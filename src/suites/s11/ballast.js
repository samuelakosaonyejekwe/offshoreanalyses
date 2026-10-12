// Ballasting and stability control: the inclination caused by a load on deck, the tank transfers (or fill and
// discharge) that bring the unit back to its target heel, trim and draft, the pumping time, the free-surface loss
// of the new tank state, and an exact re-solve of the floating equilibrium with the planned tank contents.
import { meshHydrostatics, solveEquilibrium, draftForVolume, freeSurfaceCorrection } from '../../engine/hydro/hydrostatics.js';
import { ballastPlan } from '../../engine/hydro/stability.js';
import { hullInputs, hullFrom, hullEnclosedVolume, G, DEG, need, run, round, rnd, fmt, adv, check, checkMin, sev, worst, linspace, fuelCo2Text } from './_common.js';

const TANK_COLS = [
  { key: 'name', label: 'Tank', type: 'text' }, { key: 'x', label: 'x (fwd of midship)', unit: 'm', type: 'number' }, { key: 'y', label: 'y (to port)', unit: 'm', type: 'number' },
  { key: 'zb', label: 'Tank bottom above keel', unit: 'm', type: 'number' }, { key: 'h', label: 'Tank height', unit: 'm', type: 'number' },
  { key: 'cap', label: 'Capacity', unit: 'm³', type: 'number' }, { key: 'vol', label: 'Content', unit: 'm³', type: 'number' },
  { key: 'l', label: 'Free-surface length', unit: 'm', type: 'number' }, { key: 'b', label: 'Free-surface breadth', unit: 'm', type: 'number' }, { key: 'rho', label: 'Liquid density', unit: 'kg/m³', type: 'number' },
];
const tk = (name, x, y) => ({ name, x, y, zb: 0.5, h: 8, cap: 800, vol: 400, l: 12.5, b: 8, rho: 1025 });

export default {
  id: 's11-ballast',
  suite: 11,
  title: 'Ballasting and stability control',
  summary: 'The tank-by-tank ballast transfer that corrects the heel and trim caused by a load change, how long the pumps need, and the stability of the unit with the new tank contents.',
  covers: [15],
  fidelity: 'numerical',
  standards: ['IMO MODU Code 2009 §4.10 (ballast systems of column-stabilised units)', 'DNV-OS-D101 (2021) Ch.2 Sec.4 (ballast and bilge systems)', 'IMO 2008 IS Code Part B ch. 3 (free-surface effects)', 'DNV-OS-C301 (2021) Ch.2 Sec.1', 'ISO 19901-5 (2021) (weight control)'],
  models: ['Small-angle hydrostatic influence coefficients (draft, heel and trim per unit ballast mass)', 'Bounded least-squares ballast plan with tank capacity and fill limits', 'Pump-rate limited transfer time', 'Free-surface correction of slack tanks', 'Nonlinear floating equilibrium of the closed hull with the planned tank contents (ballast-induced mass and centre-of-gravity change)'],
  limits: 'The plan is found on linear small-angle hydrostatics and then re-checked by an exact equilibrium of the hull surface; tanks are prismatic with a centroid that rises linearly with fill, all transfers run in parallel at the total pump capacity, and pipe losses, valve sequencing, pump curves and air-pipe limits are not modelled. Not applicable to a tension-leg platform, whose attitude is set by the tendons.',
  searchTerms: 'ballast system semi-submersible heel trim correction pump capacity',
  inputs: [
    ...hullInputs({ type: 'semi', only: ['semi', 'ship', 'barge', 'spar'] }),
    { key: 'loadT', label: 'Load placed on the unit', unit: 't', type: 'number', default: 400, min: 0, max: 50000, group: 'Load change (initial condition)', role: 'ic', help: 'Mass added to the upright unit (crane lift landed, module installed, bulk loaded). It heels and trims the unit before ballasting.' },
    { key: 'lx', label: 'Load x (forward of midship)', unit: 'm', type: 'number', default: 30, min: -250, max: 250, group: 'Load change (initial condition)', role: 'ic', help: 'Longitudinal position of the load.' },
    { key: 'ly', label: 'Load y (to port)', unit: 'm', type: 'number', default: 30, min: -80, max: 80, group: 'Load change (initial condition)', role: 'ic', help: 'Transverse position of the load.' },
    { key: 'lz', label: 'Load height above keel', unit: 'm', type: 'number', default: 50, min: 0, max: 250, group: 'Load change (initial condition)', role: 'ic', help: 'Height of the centre of gravity of the load.' },
    { key: 'tanks', label: 'Ballast tanks', type: 'table', group: 'Ballast system', role: 'ic', columns: TANK_COLS,
      default: [tk('P1 fwd', 45, 35), tk('P2', 15, 35), tk('P3', -15, 35), tk('P4 aft', -45, 35), tk('S1 fwd', 45, -35), tk('S2', 15, -35), tk('S3', -15, -35), tk('S4 aft', -45, -35)],
      help: 'Centroid plan position, tank bottom and height, capacity and current content, and the free-surface dimensions used while the tank is slack. The contents are part of the initial displacement.' },
    { key: 'mode', label: 'Ballast operation', type: 'select', options: [{ value: 'transfer', label: 'Internal transfer (total ballast constant)' }, { value: 'fill-discharge', label: 'Fill from and discharge to the sea' }], default: 'transfer', group: 'Ballast system', role: 'operation', help: 'Transfer keeps the displacement; fill and discharge can also restore the draft.' },
    { key: 'keepDraft', label: 'Restore the original draft (fill / discharge only)', type: 'bool', default: true, group: 'Ballast system', role: 'operation', showIf: { key: 'mode', in: ['fill-discharge'] }, help: 'Discharge ballast equal to the added load so that the unit returns to its draft.' },
    { key: 'fillMin', label: 'Minimum tank fill', unit: '-', type: 'number', default: 0.02, min: 0, max: 0.5, group: 'Ballast system', role: 'bc', help: 'Lowest permitted fill fraction (pump suction, stripping limit).' },
    { key: 'fillMax', label: 'Maximum tank fill', unit: '-', type: 'number', default: 0.98, min: 0.5, max: 1, group: 'Ballast system', role: 'bc', help: 'Highest permitted fill fraction (overflow and air-pipe margin).' },
    { key: 'pumpRate', label: 'Capacity of one ballast pump', unit: 'm³/h', type: 'number', default: 500, min: 10, max: 10000, group: 'Ballast system', role: 'bc', design: true, help: 'Rated flow of one pump at the system head.' },
    { key: 'nPumps', label: 'Pumps running', unit: '-', type: 'int', default: 2, min: 1, max: 12, group: 'Ballast system', role: 'bc', help: 'Number of pumps available for the operation (consider one out of service).' },
    { key: 'pumpHead', label: 'Pump head', unit: 'm', type: 'number', default: 30, min: 2, max: 120, group: 'Ballast system', role: 'bc', help: 'Total head for the energy estimate.' },
    { key: 'pumpEff', label: 'Pump and drive efficiency', unit: '-', type: 'number', default: 0.65, min: 0.2, max: 0.95, group: 'Ballast system', role: 'bc', help: 'Overall efficiency from electrical power to hydraulic power.' },
    { key: 'heelTarget', label: 'Target heel', unit: 'deg', type: 'number', default: 0, min: -10, max: 10, group: 'Targets and acceptance criteria', role: 'criteria', help: 'Positive = port side up.' },
    { key: 'trimTarget', label: 'Target trim', unit: 'deg', type: 'number', default: 0, min: -10, max: 10, group: 'Targets and acceptance criteria', role: 'criteria', help: 'Positive = bow down.' },
    { key: 'tolDeg', label: 'Allowable residual inclination', unit: 'deg', type: 'number', default: 0.25, min: 0.01, max: 5, group: 'Targets and acceptance criteria', role: 'criteria', help: 'Acceptable difference from the targets after ballasting.' },
    { key: 'timeMax', label: 'Allowable ballasting time', unit: 'h', type: 'number', default: 3, min: 0.05, max: 48, group: 'Targets and acceptance criteria', role: 'criteria', help: 'IMO MODU Code asks that a column-stabilised unit can go from operating to severe-storm draft within 3 h; lift operations often need far less.' },
    { key: 'gmMin', label: 'Required metacentric height after ballasting', unit: 'm', type: 'number', default: 1.0, min: 0, max: 10, group: 'Targets and acceptance criteria', role: 'criteria', help: 'Minimum GM including the free-surface correction of the tanks left slack.' },
  ],
  compute(inp) {
    const hull = hullFrom(inp); const rho = hull.rho; const mesh = hull.mesh;
    need(hull.type !== 'tlp', 'Ballast trimming of a tension-leg platform is governed by tendon tension, not by free-floating hydrostatics: use another hull type.');
    need(inp.fillMax > inp.fillMin, 'The maximum tank fill must exceed the minimum fill.');
    const rows = (inp.tanks || []).filter((t) => t.cap > 0);
    need(rows.length >= 2 || (rows.length === 1 && inp.mode === 'fill-discharge'), 'An internal transfer needs at least two ballast tanks with a capacity; add tanks or switch to fill and discharge.');
    const tanks = rows.map((t) => {
      const r = t.rho > 0 ? t.rho : 1025; const lo = t.cap * inp.fillMin; const hi = t.cap * inp.fillMax;
      need(t.vol >= 0 && t.vol <= t.cap, `Tank "${t.name}": the content (${fmt(t.vol)} m³) must lie between 0 and the capacity (${fmt(t.cap)} m³).`);
      need(t.h > 0, `Tank "${t.name}": enter a positive tank height.`);
      // the planner works on the usable range between the fill limits; contents outside it can only move inwards
      const v = Math.min(Math.max(t.vol, lo), hi);
      return { ...t, rho: r, lo, hi, plan: { name: t.name, x: t.x, y: t.y, capacity: hi - lo, volume: v - lo, rho: r } };
    });
    const outside = tanks.filter((t) => Math.abs(t.x) > hull.Lc / 2 + 1e-6 || Math.abs(t.y) > hull.Bc / 2 + 1e-6);
    // --- initial condition: upright unit plus the load -----------------------------------------------------
    const m1 = hull.mass + inp.loadT * 1e3;
    const full = hullEnclosedVolume(hull) * rho;
    need(m1 < full, `With the load the unit weighs ${fmt(m1 / 1e3)} t, more than the buoyancy of the fully immersed hull (${fmt(full / 1e3)} t).`);
    const cog1 = [0, 1, 2].map((i) => (hull.mass * hull.cog[i] + inp.loadT * 1e3 * [inp.lx, inp.ly, inp.lz][i]) / m1);
    const solve = (mass, cog, guess) => {
      const eq = run(() => solveEquilibrium({ mesh, mass, cog, rho, guess }), 'Floating equilibrium');
      need(eq.converged && Math.abs(eq.heel) < 60 * DEG && Math.abs(eq.trim) < 60 * DEG, 'No stable floating equilibrium was found: the load capsizes or submerges the unit. Reduce the load or its offset, or lower the centre of gravity.');
      return eq;
    };
    const eq1 = solve(m1, cog1, { draft: hull.T });
    const T1 = run(() => draftForVolume({ mesh, volume: m1 / rho }), 'Draft with the load');
    const hs1 = meshHydrostatics({ mesh, draft: T1, rho, cog: cog1, mass: m1 });
    const slack = (list, vols) => list.filter((t, i) => vols[i] > 1e-6 * t.cap && vols[i] < t.cap * (1 - 1e-6) && t.l > 0 && t.b > 0).map((t) => ({ l: t.l, b: t.b, rho: t.rho }));
    const fsc1 = freeSurfaceCorrection({ tanks: slack(tanks, tanks.map((t) => t.vol)), volume: hs1.volume, rho });
    const GMt1 = hs1.GMt - fsc1.GGt; const GMl1 = hs1.GMl - fsc1.GGl;
    need(GMt1 > 0 && GMl1 > 0, `The unit is unstable with the load on board (GM = ${fmt(Math.min(GMt1, GMl1))} m after free-surface correction): lower the load or press tanks full before lifting.`);
    // --- plan ------------------------------------------------------------------------------------------------
    const target = { heelDeg: inp.heelTarget, trimDeg: inp.trimTarget };
    if (inp.mode === 'fill-discharge' && inp.keepDraft) target.draft = hull.T;
    const plan = run(() => ballastPlan({ tanks: tanks.map((t) => t.plan), hydro: { displacement: m1, waterplaneArea: hs1.waterplaneArea, lcf: hs1.cofBody[0], GMt: GMt1, GMl: GMl1, rho },
      current: { draft: T1, heelDeg: eq1.heel / DEG, trimDeg: eq1.trim / DEG }, target, mode: inp.mode, pumpRate: inp.pumpRate / 3600, nPumps: inp.nPumps }), 'Ballast plan');
    const dv = plan.tanks.map((p) => p.deltaVolume);
    // --- exact re-check with the planned tank contents -------------------------------------------------------
    // contents of a prismatic tank: centroid at zb + (v/cap)·h/2, so the vertical moment changes with the fill level too
    const state = (f) => {
      let m = m1; const mom = cog1.map((c) => c * m1); const vols = [];
      tanks.forEach((t, i) => {
        const v0 = t.vol; const v = v0 + f * dv[i] + (f > 0 ? f * (t.plan.volume + t.lo - t.vol) : 0); vols.push(v);
        const zc = (q) => t.zb + (q / t.cap) * t.h / 2;
        m += t.rho * (v - v0); mom[0] += t.rho * (v - v0) * t.x; mom[1] += t.rho * (v - v0) * t.y; mom[2] += t.rho * (v * zc(v) - v0 * zc(v0));
      });
      return { m, cog: mom.map((q) => q / m), vols };
    };
    const fr = linspace(0, 1, 9); let guess = { draft: eq1.draft, heel: eq1.heel, trim: eq1.trim };
    const path = fr.map((f) => { const s = state(f); const e = solve(s.m, s.cog, guess); guess = { draft: e.draft, heel: e.heel, trim: e.trim }; return { f, ...e, m: s.m }; });
    const end = path[path.length - 1]; const s2 = state(1);
    const T2 = run(() => draftForVolume({ mesh, volume: s2.m / rho }), 'Final draft');
    const hs2 = meshHydrostatics({ mesh, draft: T2, rho, cog: s2.cog, mass: s2.m });
    const fsc2 = freeSurfaceCorrection({ tanks: slack(tanks, s2.vols), volume: hs2.volume, rho });
    const GM2 = hs2.GMt - fsc2.GGt;
    const tPump = plan.pumpTime / 3600;
    const errHeel = Math.abs(end.heel / DEG - inp.heelTarget); const errTrim = Math.abs(end.trim / DEG - inp.trimTarget);
    const energy = rho * G * plan.volumeMoved * inp.pumpHead / inp.pumpEff / 3.6e6; // kWh
    const fuelT = energy * 0.21 / 1e3; // 0.21 kg of marine gas oil per kWh of diesel-electric generation (screening)
    const checks = [
      check('Residual heel after ballasting (exact equilibrium)', errHeel, inp.tolDeg, 'deg', 'IMO MODU Code 2009 §4.10.1; operating-manual inclination limits', 'Guards against leaving the unit listed: a residual heel reduces air gap on the low side, loads cranes and process equipment sideways and eats into the stability range.'),
      check('Residual trim after ballasting (exact equilibrium)', errTrim, inp.tolDeg, 'deg', 'IMO MODU Code 2009 §4.10.1; operating-manual inclination limits', 'Guards against a residual fore-and-aft inclination, which drilling, process trains and risers tolerate even less than heel.'),
      check('Ballasting time', tPump, inp.timeMax, 'h', 'IMO MODU Code 2009 §4.10.1 (ballast pump capacity; 3 h from operating to storm draft); DNV-OS-D101 (2021) Ch.2 Sec.4', 'Guards against a ballast system too slow for the operation: the unit stays inclined, or exposed at the wrong draft, while the pumps catch up.'),
      checkMin('Metacentric height after ballasting, with free surfaces', GM2, inp.gmMin, 'm', 'IMO 2008 IS Code Part B ch. 3; IMO MODU Code 2009 §3.2', 'Guards initial stability in the new tank state: every tank left slack by the transfer raises the effective centre of gravity.'),
    ];
    const gov = worst(checks);
    const tot0 = tanks.reduce((a, t) => a + t.rho * t.vol, 0) / 1e3; const tot2 = tanks.reduce((a, t, i) => a + t.rho * s2.vols[i], 0) / 1e3;
    const notes = [
      `Load of ${fmt(inp.loadT)} t at (${fmt(inp.lx)}, ${fmt(inp.ly)}, ${fmt(inp.lz)}) m inclines the unit to ${fmt(eq1.heel / DEG)}° heel and ${fmt(eq1.trim / DEG)}° trim (positive = port up, bow down) at ${fmt(T1, 4)} m draft.`,
      `${inp.mode === 'transfer' ? 'Internal transfer' : 'Fill and discharge'} of ${fmt(plan.volumeMoved)} m³ in ${fmt(tPump)} h with ${inp.nPumps} × ${fmt(inp.pumpRate)} m³/h. The linear plan predicts ${fmt(plan.achieved.heelDeg)}° / ${fmt(plan.achieved.trimDeg)}°; the exact equilibrium of the hull gives ${fmt(end.heel / DEG)}° / ${fmt(end.trim / DEG)}°.`,
      `${plan.feasible ? 'All targets are reachable within the tank fill limits.' : 'The targets cannot be reached within the tank capacities and fill limits: the plan is the closest achievable state.'} Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
    ];
    if (outside.length) notes.push(`Tank${outside.length > 1 ? 's' : ''} ${outside.map((t) => `"${t.name}"`).join(', ')} lie outside the plan extent of the hull (${fmt(hull.Lc)} × ${fmt(hull.Bc)} m): check the tank table against the hull type.`);
    return {
      kpis: [
        { key: 'heel0', label: 'Heel before ballasting', value: round(eq1.heel / DEG, 5), unit: 'deg', explain: 'Heel caused by the load, from the exact equilibrium of the hull surface. Positive = port side up.' },
        { key: 'trim0', label: 'Trim before ballasting', value: round(eq1.trim / DEG, 5), unit: 'deg', explain: 'Trim caused by the load. Positive = bow down.' },
        { key: 'draft0', label: 'Mean draft with the load', value: round(T1, 5), unit: 'm', explain: 'Upright-equivalent draft at the displacement including the load.' },
        { key: 'GM0', label: 'GM with the load (free-surface corrected)', value: round(GMt1, 4), unit: 'm', explain: 'Transverse metacentric height the ballast plan is computed with.' },
        { key: 'volumeMoved', label: 'Ballast volume moved', value: round(plan.volumeMoved, 5), unit: 'm³', explain: 'Volume pumped: the larger of the amounts removed and added for a transfer, their sum for fill and discharge.' },
        { key: 'pumpTime', label: 'Ballasting time', value: round(tPump, 5), unit: 'h', limit: inp.timeMax, sense: 'max', explain: 'Volume moved divided by the total running pump capacity. Real operations add valve line-up and stripping time.' },
        { key: 'heelFinal', label: 'Heel after ballasting', value: round(end.heel / DEG, 5), unit: 'deg', explain: 'Heel from the exact equilibrium with the planned tank contents; the difference from the target is the error of the linear plan plus any infeasibility.' },
        { key: 'trimFinal', label: 'Trim after ballasting', value: round(end.trim / DEG, 5), unit: 'deg', explain: 'Trim from the exact equilibrium with the planned tank contents.' },
        { key: 'draftFinal', label: 'Mean draft after ballasting', value: round(T2, 5), unit: 'm', explain: 'Upright-equivalent draft after the operation; unchanged by an internal transfer.' },
        { key: 'GMfinal', label: 'GM after ballasting (free-surface corrected)', value: round(GM2, 4), unit: 'm', limit: inp.gmMin, sense: 'min', explain: 'Transverse metacentric height in the final tank state including the virtual rise of G from slack tanks.' },
        { key: 'fsc', label: 'Free-surface correction after ballasting', value: round(fsc2.GGt, 5), unit: 'm', explain: 'Virtual rise of the centre of gravity caused by the tanks left partly filled.' },
        { key: 'ballastTotal', label: 'Total ballast after the operation', value: round(tot2, 5), unit: 't', explain: `Mass of liquid in the listed tanks after the operation (${fmt(tot0)} t before).` },
        { key: 'pumpEnergy', label: 'Pumping energy', value: round(energy, 4), unit: 'kWh', explain: 'Electrical energy for the transfer at the stated head and efficiency.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Tank contents before and after', categories: tanks.map((t) => t.name), y: { label: 'Content', unit: 'm³' },
          series: [{ name: 'Before', y: rnd(tanks.map((t) => t.vol)) }, { name: 'After', y: rnd(s2.vols) }, { name: 'Capacity', y: rnd(tanks.map((t) => t.cap)) }],
          explain: 'Liquid volume in each tank before and after the operation against its capacity. Tanks driven to a fill limit show where the system runs out of authority.' },
        { type: 'line', title: 'Attitude during the operation', x: { label: 'Time', unit: 'h' }, y: { label: 'Heel / trim', unit: 'deg' },
          series: [{ name: 'Heel', x: rnd(fr.map((f) => f * tPump)), y: rnd(path.map((p) => p.heel / DEG)) }, { name: 'Trim', x: rnd(fr.map((f) => f * tPump)), y: rnd(path.map((p) => p.trim / DEG)) }],
          explain: 'Heel and trim from exact equilibria as the transfers progress in parallel. Both should move monotonically to the targets; a curve that overshoots means the tanks are being used against each other.' },
        { type: 'scatter', title: 'Tank plan and volume change', x: { label: 'x (forward)', unit: 'm' }, y: { label: 'y (port)', unit: 'm' },
          series: [{ name: 'Filled', x: tanks.filter((t, i) => dv[i] > 1e-6).map((t) => t.x), y: tanks.filter((t, i) => dv[i] > 1e-6).map((t) => t.y) }, { name: 'Emptied', x: tanks.filter((t, i) => dv[i] < -1e-6).map((t) => t.x), y: tanks.filter((t, i) => dv[i] < -1e-6).map((t) => t.y) }, { name: 'Unchanged', x: tanks.filter((t, i) => Math.abs(dv[i]) <= 1e-6).map((t) => t.x), y: tanks.filter((t, i) => Math.abs(dv[i]) <= 1e-6).map((t) => t.y) }, { name: 'Load', x: [inp.lx], y: [inp.ly] }],
          explain: 'Plan view of the tanks. Ballast should leave the tanks near the load and go to the tanks diagonally opposite; anything else indicates a tank at its limit.' },
      ],
      tables: [
        { title: 'Ballast plan', columns: [{ key: 'name', label: 'Tank' }, { key: 'v0', label: 'Before', unit: 'm³' }, { key: 'dv', label: 'Change', unit: 'm³' }, { key: 'v1', label: 'After', unit: 'm³' }, { key: 'fill', label: 'Fill after', unit: '%' }, { key: 'dm', label: 'Mass change', unit: 't' }, { key: 't', label: 'Time at one pump', unit: 'h' }],
          rows: tanks.map((t, i) => ({ name: t.name, v0: round(t.vol, 5), dv: round(s2.vols[i] - t.vol, 5), v1: round(s2.vols[i], 5), fill: round(100 * s2.vols[i] / t.cap, 4), dm: round(t.rho * (s2.vols[i] - t.vol) / 1e3, 5), t: round(Math.abs(s2.vols[i] - t.vol) / inp.pumpRate, 4) })),
          explain: 'One row per tank: the pumping instruction (volume change), the resulting fill level and the time one pump would need for that tank alone. Use it as the ballast sequence sheet.' },
      ],
      notes,
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: residual ${fmt(errHeel)}° heel and ${fmt(errTrim)}° trim after moving ${fmt(plan.volumeMoved)} m³ in ${fmt(tPump)} h; GM ${fmt(GM2)} m.`,
          !plan.feasible ? 'Free capacity where it is needed: lower the content of the tanks under the load or raise it in the opposite tanks before the operation, widen the fill limits, or add tanks further from the centre of flotation.'
            : tPump > inp.timeMax ? `Run ${Math.ceil(inp.nPumps * tPump / inp.timeMax)} pumps or raise the pump capacity to ${fmt(inp.pumpRate * tPump / inp.timeMax)} m³/h each to finish within ${fmt(inp.timeMax)} h.`
              : GM2 < inp.gmMin ? `Press up or empty the slack tanks: free surfaces cost ${fmt(fsc2.GGt)} m of GM; ${fmt(inp.gmMin - GM2)} m must be recovered.`
                : errHeel > inp.tolDeg || errTrim > inp.tolDeg ? 'The linear plan leaves a residual on the real hull: repeat the plan from the achieved attitude (a second, smaller transfer) or trim manually on the inclinometers.'
                  : `The plan meets the targets with ${fmt(inp.timeMax - tPump)} h of time margin. Pre-ballast against the load before it lands so that the unit never reaches ${fmt(Math.max(Math.abs(eq1.heel), Math.abs(eq1.trim)) / DEG)}°, and confirm tank soundings against the plan.`,
          'IMO MODU Code 2009 §4.10; DNV-OS-D101 (2021) Ch.2 Sec.4; IMO 2008 IS Code Part B ch. 3.'),
        adv('economic', 'info', `The operation takes ${fmt(tPump)} h of pump running and ${fmt(energy)} kWh; while it lasts the lift or loading operation is on hold.`,
          `Size the pumps on the critical-path operation: each halving of the ballasting time needs twice the pump capacity but returns ${fmt(tPump / 2)} h of crane or vessel time per lift.`,
          'Operational planning per DNV-ST-N001 (2021) (weather-restricted operations: planned duration plus contingency); whole-life costing per ISO 15686-5.',
          `${fuelCo2Text(fuelT)} (0.21 kg fuel per kWh of diesel-electric generation, screening value). Internal transfer avoids taking on and discharging sea water, and with it ballast-water treatment energy and invasive-species risk.`),
      ],
    };
  },
  verify: [
    { // Box barge 100 × 20 m, draft 5 m, KG 4 m: Δ = 10 250 t. A 100 t load at y = 5 m, z = 4 m (same height as G) gives
      // Δ′ = 10 350 t, T′ = 5.0488 m, KB = 2.5244, BM = 20²/(12·5.0488) = 6.6023, GM = 9.1267 − 4 = 5.1267 m,
      // tan φ = 100·5/(10 350·5.1267) = 0.009423 → φ = 0.5399° (port down, i.e. negative in the engine convention).
      // Transfer between tanks at y = ±8 m: ρ·v·16 = 100 t·5 m → v = 500/(1.025·16) = 30.49 m³; at 100 m³/h: 0.3049 h.
      name: 'Weight shift on a box barge: heel w·d/(Δ·GM) and the transfer volume that cancels it',
      source: 'Biran & López-Pulido, Ship Hydrostatics and Stability (2014) ch. 7 — hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, loadT: 100, lx: 0, ly: 5, lz: 4, mode: 'transfer', pumpRate: 100, nPumps: 1, fillMin: 0, fillMax: 1,
        tanks: [{ name: 'P', x: 0, y: 8, zb: 0, h: 4, cap: 1000, vol: 500, l: 0, b: 0, rho: 1025 }, { name: 'S', x: 0, y: -8, zb: 0, h: 4, cap: 1000, vol: 500, l: 0, b: 0, rho: 1025 }] },
      expect: [{ kpi: 'heel0', value: -0.5399, tol: 5e-3 }, { kpi: 'draft0', value: 5.0488, tol: 1e-4 }, { kpi: 'volumeMoved', value: 30.49, tol: 5e-3 }, { kpi: 'pumpTime', value: 0.3049, tol: 5e-3 }, { kpi: 'heelFinal', value: 0, abs: 0.01 }] },
    { name: 'No load: nothing to pump', source: 'Limiting case — an upright unit at its targets needs no ballast movement', inputs: { loadT: 0 },
      expect: [{ kpi: 'volumeMoved', value: 0, abs: 1e-3 }, { kpi: 'heel0', value: 0, abs: 1e-6 }, { kpi: 'heelFinal', value: 0, abs: 1e-4 }] },
  ],
};
