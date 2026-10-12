// Progressive flooding of one compartment through one opening: Torricelli inflow driven by the head difference,
// quasi-static re-equilibration of the hull at every step, time to equalise or to reach the limiting inclination,
// submergence of a second (downflooding) opening and the dewatering capacity that would be needed.
import { bodyToEarth } from '../../engine/hydro/hydrostatics.js';
import { progressiveFlooding } from '../../engine/hydro/stability.js';
import { hullInputs, hullFrom, hullPlot, G, DEG, need, run, round, rnd, fmt, adv, check, checkMin, sev, worst, thin, fuelCo2Text } from './_common.js';

export default {
  id: 's11-progressive-flooding',
  suite: 11,
  title: 'Progressive flooding in time',
  summary: 'How fast a breached compartment fills, how the draft, heel and trim develop while it does, and how much time there is before the unit reaches its limiting inclination or a second opening goes under.',
  covers: [12, 17],
  fidelity: 'numerical',
  standards: ['IMO MODU Code 2009 §3.4 (damage stability), §3.6 (watertight integrity)', 'SOLAS 2020 Ch. II-1 Reg. 7-2, 8 (intermediate stages of flooding)', 'DNV-OS-C301 (2021) Ch.2 Sec.1 [5]', 'IMO MSC.362(92) (time-domain flooding guidance)', 'NORSOK N-003 (2017) §9 (accidental actions)'],
  models: ['Torricelli orifice flow Q = Cd·A·√(2g·Δh)', 'Added-weight flooding with the floodwater centroid from the clipped compartment', 'Nonlinear floating equilibrium at each time step (quasi-static)', 'Explicit Euler time integration', 'Water-volume conservation check'],
  limits: 'Progressive flooding is one compartment and one opening to the sea, starting empty: no flow between compartments, no air compression or venting, no sloshing or transient roll (each step is a static equilibrium), explicit Euler in time. The run stops when the levels equalise, the compartment is full or no floating equilibrium is found, which is reported as loss of the unit although it can also mean capsize beyond the solver range. Final-stage survivability of several compartments is the damage stability module.',
  searchTerms: 'progressive flooding time domain damaged ship compartment orifice',
  inputs: [
    ...hullInputs({ type: 'semi', only: ['semi', 'ship', 'barge', 'spar'] }),
    { key: 'cxMin', label: 'Compartment aft end x', unit: 'm', type: 'number', default: 27.5, min: -300, max: 300, group: 'Flooded compartment', role: 'geometry', help: 'Aft boundary of the box-shaped compartment (x forward of midship).' },
    { key: 'cxMax', label: 'Compartment forward end x', unit: 'm', type: 'number', default: 42.5, min: -300, max: 300, group: 'Flooded compartment', role: 'geometry', help: 'Forward boundary.' },
    { key: 'cyMin', label: 'Compartment starboard side y', unit: 'm', type: 'number', default: 35, min: -100, max: 100, group: 'Flooded compartment', role: 'geometry', help: 'Starboard boundary (y positive to port).' },
    { key: 'cyMax', label: 'Compartment port side y', unit: 'm', type: 'number', default: 42.5, min: -100, max: 100, group: 'Flooded compartment', role: 'geometry', help: 'Port boundary.' },
    { key: 'czMin', label: 'Compartment bottom above keel', unit: 'm', type: 'number', default: 16, min: 0, max: 250, group: 'Flooded compartment', role: 'geometry', help: 'Lower watertight boundary (tank top or flat).' },
    { key: 'czMax', label: 'Compartment top above keel', unit: 'm', type: 'number', default: 26, min: 0.5, max: 250, group: 'Flooded compartment', role: 'geometry', help: 'Upper watertight boundary (deck or flat).' },
    { key: 'perm', label: 'Permeability', unit: '-', type: 'number', default: 0.95, min: 0.05, max: 1, group: 'Flooded compartment', role: 'ic', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Share of the compartment volume that can fill with water: 0.95 voids and tanks, 0.85 machinery, 0.6 stores. The compartment starts empty.' },
    { key: 'ox', label: 'Opening x', unit: 'm', type: 'number', default: 35, min: -300, max: 300, group: 'Opening to the sea', role: 'bc', help: 'Position of the centre of the breach or open valve.' },
    { key: 'oy', label: 'Opening y', unit: 'm', type: 'number', default: 42.5, min: -100, max: 100, group: 'Opening to the sea', role: 'bc', help: 'Transverse position of the opening.' },
    { key: 'oz', label: 'Opening height above keel', unit: 'm', type: 'number', default: 19, min: 0, max: 250, group: 'Opening to the sea', role: 'bc', help: 'Height of the centre of the opening; it must be below the waterline for water to enter.' },
    { key: 'oArea', label: 'Opening area', unit: 'm²', type: 'number', default: 0.3, min: 0.001, max: 50, group: 'Opening to the sea', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.5 }, design: true, help: 'Area of the breach: 0.03 m² for a 200 mm pipe, 0.3 m² for a supply-vessel bow puncture, several m² for a major collision.' },
    { key: 'oCd', label: 'Discharge coefficient', unit: '-', type: 'number', default: 0.6, min: 0.3, max: 1, group: 'Opening to the sea', role: 'bc', help: '0.6 for a sharp-edged hole; ragged damage openings are similar.' },
    { key: 'dfX', label: 'Second opening x', unit: 'm', type: 'number', default: 35, min: -300, max: 300, group: 'Downflooding point', role: 'bc', help: 'An unprotected opening (vent, door, chain locker) that would extend the flooding if it reached the water.' },
    { key: 'dfY', label: 'Second opening y', unit: 'm', type: 'number', default: 42.5, min: -100, max: 100, group: 'Downflooding point', role: 'bc', help: 'Transverse position.' },
    { key: 'dfZ', label: 'Second opening height above keel', unit: 'm', type: 'number', default: 38, min: 0, max: 300, group: 'Downflooding point', role: 'bc', help: 'Height of its lower edge.' },
    { key: 'pumpCap', label: 'Dewatering pump capacity', unit: 'm³/h', type: 'number', default: 500, min: 0, max: 20000, group: 'Response', role: 'operation', help: 'Bilge or ballast pump capacity available to the compartment, for comparison with the inflow (it is not subtracted from the inflow).' },
    { key: 'heelMax', label: 'Limiting inclination', unit: 'deg', type: 'number', default: 17, min: 2, max: 45, group: 'Acceptance criteria', role: 'criteria', help: 'IMO MODU Code: 17° for a column-stabilised unit after damage; 15° is common for surface units.' },
    { key: 'tResp', label: 'Required response time', unit: 'min', type: 'number', default: 15, min: 0, max: 720, group: 'Acceptance criteria', role: 'criteria', help: 'Time needed to detect, close valves, counter-ballast or muster before the limiting inclination or the second opening is reached.' },
    { key: 'tMaxMin', label: 'Simulated duration', unit: 'min', type: 'number', default: 60, min: 1, max: 720, group: 'Numerical controls', role: 'control', help: 'The run stops earlier when the levels equalise.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 10, min: 1, max: 120, group: 'Numerical controls', role: 'control', help: 'Explicit Euler step; halve it until the time to equalise stops changing.' },
  ],
  meshTime: { key: 'dt', label: 'Time step', levels: [20, 10, 5], qoi: 'tEqualise', isSize: true },
  compute(inp) {
    const hull = hullFrom(inp, { nTheta: 24 }); const rho = hull.rho;
    need(hull.type !== 'tlp', 'Flooding of a tension-leg platform changes tendon tension rather than attitude and is not covered by this free-floating model.');
    need(inp.cxMax > inp.cxMin && inp.cyMax > inp.cyMin && inp.czMax > inp.czMin, 'The compartment has no volume: each maximum must exceed its minimum.');
    need(inp.oz < hull.T, `The opening at ${fmt(inp.oz)} m above keel is above the waterline (draft ${fmt(hull.T)} m): no water enters. Lower the opening or deepen the draft.`);
    const tol = 1e-6;
    need(inp.ox >= inp.cxMin - tol && inp.ox <= inp.cxMax + tol && inp.oy >= inp.cyMin - tol && inp.oy <= inp.cyMax + tol && inp.oz >= inp.czMin - tol && inp.oz <= inp.czMax + tol,
      'The opening must lie on or inside the compartment it floods: move the opening onto a compartment boundary.');
    const box = { xMin: inp.cxMin, xMax: inp.cxMax, yMin: inp.cyMin, yMax: inp.cyMax, zMin: inp.czMin, zMax: inp.czMax };
    const vComp = (box.xMax - box.xMin) * (box.yMax - box.yMin) * (box.zMax - box.zMin) * inp.perm;
    need(inp.tMaxMin * 60 / inp.dt <= 4000, `The run needs ${fmt(inp.tMaxMin * 60 / inp.dt)} steps: raise the time step or shorten the duration (limit 4000 steps).`);
    const r = run(() => progressiveFlooding({ mesh: hull.mesh, mass: hull.mass, cog: hull.cog, compartment: { box, permeability: inp.perm }, opening: { point: [inp.ox, inp.oy, inp.oz], area: inp.oArea, Cd: inp.oCd }, dt: inp.dt, tMax: inp.tMaxMin * 60, rho }), 'Progressive flooding');
    const n = r.t.length; const last = n - 1;
    const incl = r.heelDeg.map((h, i) => Math.atan(Math.hypot(Math.tan(h * DEG), Math.tan(r.trimDeg[i] * DEG))) / DEG);
    const dfH = r.t.map((_, i) => bodyToEarth([inp.dfX, inp.dfY, inp.dfZ], r.draft[i], r.heelDeg[i] * DEG, r.trimDeg[i] * DEG)[2]);
    const firstAt = (arr, pred) => { const i = arr.findIndex(pred); return i < 0 ? null : r.t[i]; };
    const tHeel = firstAt(incl, (v) => v >= inp.heelMax); const tDf = firstAt(dfH, (v) => v <= 0);
    const tCrit = r.sunk ? r.t[last] : [tHeel, tDf].filter((v) => v !== null).reduce((a, v) => Math.min(a, v), Infinity);
    const horizon = inp.tMaxMin * 60;
    const inclMax = Math.max(...incl); const dfMin = Math.min(...dfH); const qMax = Math.max(...r.flow);
    // conservation: the stored floodwater must equal the time integral of the inflow (explicit Euler, same quadrature)
    let vIn = 0; for (let i = 0; i < last; i++) vIn += r.flow[i] * inp.dt;
    const balance = r.floodMass[last] / rho - Math.min(vIn, vComp);
    const tEq = r.timeToEqualise;
    const h0 = hull.T - inp.oz;
    const aPump = inp.pumpCap > 0 ? inp.pumpCap / 3600 / (inp.oCd * Math.sqrt(2 * G * h0)) : 0;
    const checks = [
      check('Maximum inclination during flooding', r.sunk ? Math.max(inclMax, inp.heelMax * 1.5) : inclMax, inp.heelMax, 'deg', 'IMO MODU Code 2009 §3.4.3 (inclination after damage); SOLAS Ch. II-1 Reg. 7-2 (intermediate stages)', 'Guards against an inclination at which evacuation, ballast systems and equipment stop working, at any stage of the flooding and not only at the end.'),
      checkMin('Time before the limiting inclination or downflooding', Number.isFinite(tCrit) ? Math.max(tCrit / 60, 1e-6) : horizon / 60, inp.tResp, 'min', 'IMO MODU Code 2009 §3.6 (watertight integrity, remote closing); emergency-response time (project criterion)', 'Guards the time the crew needs to react. If the limit is reached sooner than valves can be closed or counter-ballast applied, the damage case is not controllable.'),
      checkMin('Second opening above water', Math.max(dfMin, 0), 0, 'm', 'IMO MODU Code 2009 §3.4.3 (openings below the final waterline to be watertight)', 'Guards against progressive flooding of further spaces: once an unprotected opening is immersed the flooding no longer stops at one compartment.'),
    ];
    if (dfMin <= 0) checks[2].uc = 99;
    if (!Number.isFinite(tCrit)) checks[1].uc = 0;
    const gov = worst(checks);
    const [tt, fm, dr, he, tr, hd, fl, ic, df] = thin([r.t, r.floodMass, r.draft, r.heelDeg, r.trimDeg, r.head, r.flow, incl, dfH], 400);
    const tm = rnd(tt.map((t) => t / 60));
    const p3 = hullPlot(hull.mesh, hull.T, 'Hull with the flooded compartment', 'The hull surface coloured by still-water pressure head, with the flooded compartment drawn as a wire box. Check that the compartment sits where you intended and that the opening is below the waterline.');
    const n0 = p3.nodes.length;
    for (const z of [box.zMin, box.zMax]) for (const [x, y] of [[box.xMin, box.yMin], [box.xMax, box.yMin], [box.xMax, box.yMax], [box.xMin, box.yMax]]) { p3.nodes.push([x, y, z]); p3.nodeValues.push(0); }
    for (let i = 0; i < 4; i++) p3.lines.push([n0 + i, n0 + (i + 1) % 4], [n0 + 4 + i, n0 + 4 + (i + 1) % 4], [n0 + i, n0 + 4 + i]);
    const fuelT = r.floodMass[last] / rho * rho * G * 20 / 0.65 / 3.6e6 * 0.21 / 1e3; // dewatering against 20 m head
    return {
      kpis: [
        { key: 'floodMass', label: 'Floodwater on board at the end', value: round(r.floodMass[last] / 1e3, 5), unit: 't', explain: 'Mass of sea water in the compartment when the run stops (equalised, full, or end of the simulated time).' },
        { key: 'fillFrac', label: 'Compartment filling', value: round(r.floodMass[last] / (rho * vComp), 5), unit: '-', explain: 'Floodwater volume over the floodable volume of the compartment; 1.0 means the compartment is full.' },
        { key: 'tEqualise', label: 'Time to equalise', value: round((tEq ?? horizon) / 60, 5), unit: 'min', explain: tEq === null ? 'The flooding had not stopped at the end of the simulated time (or the unit was lost): the value shown is the simulated duration.' : 'Time until the inside level reaches the sea level or the compartment is full and the inflow stops.' },
        { key: 'equalised', label: 'Flooding stopped (1 = yes)', value: r.equalised ? 1 : 0, unit: '-', explain: 'Whether the inflow stopped within the simulated time. Zero with a lost unit means no floating equilibrium was found.' },
        { key: 'lost', label: 'Unit lost (1 = yes)', value: r.sunk ? 1 : 0, unit: '-', explain: 'No floating equilibrium exists with the floodwater on board: the unit sinks or capsizes in this damage case.' },
        { key: 'draftEnd', label: 'Draft at the end', value: round(r.draft[last], 5), unit: 'm', explain: 'Depth of the keel origin below the waterline in the last computed state.' },
        { key: 'heelEnd', label: 'Heel at the end', value: round(r.heelDeg[last], 5), unit: 'deg', explain: 'Heel in the last computed state. Positive = port side up, so port-side damage gives a negative value.' },
        { key: 'trimEnd', label: 'Trim at the end', value: round(r.trimDeg[last], 5), unit: 'deg', explain: 'Trim in the last computed state. Positive = bow down.' },
        { key: 'inclMax', label: 'Largest inclination during flooding', value: round(inclMax, 5), unit: 'deg', limit: inp.heelMax, sense: 'max', explain: 'Largest slope of the deck plane (heel and trim combined) at any time during the flooding.' },
        { key: 'tLimit', label: 'Time to the limiting inclination or downflooding', value: round((Number.isFinite(tCrit) ? tCrit : horizon) / 60, 5), unit: 'min', explain: Number.isFinite(tCrit) ? 'Time at which the limiting inclination is reached, the second opening goes under, or the unit is lost: the time available to respond.' : 'Neither limit is reached within the simulated time; the value shown is the simulated duration.' },
        { key: 'q0', label: 'Initial inflow', value: round(r.flow[0] * 3600, 5), unit: 'm³/h', explain: 'Flow through the opening at the first instant, Cd·A·√(2g·h) with h the depth of the opening below the still waterline.' },
        { key: 'qMax', label: 'Peak inflow', value: round(qMax * 3600, 5), unit: 'm³/h', explain: 'Largest flow through the opening during the run. It can exceed the initial value when the unit heels towards the damage and the opening sinks deeper. Compare with the dewatering capacity.' },
        { key: 'dfMin', label: 'Lowest height of the second opening above water', value: round(dfMin, 5), unit: 'm', limit: 0, sense: 'min', explain: 'Minimum clearance of the downflooding point above the sea surface during the flooding; negative means it is immersed.' },
        { key: 'volBalance', label: 'Water-volume conservation error', value: round(balance, 6), unit: 'm³', explain: 'Stored floodwater minus the time integral of the inflow. It should be zero to rounding; anything else means the numerical run cannot be trusted.' },
        { key: 'pumpArea', label: 'Opening area the pumps can hold', value: round(aPump, 5), unit: 'm²', explain: 'Largest opening at this depth whose initial inflow equals the dewatering capacity. Real breaches are usually far larger: pumps buy time, they rarely stop flooding.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Floodwater and inflow', x: { label: 'Time', unit: 'min' }, y: { label: 'Mass (t) / inflow (m³/min)', unit: '-' },
          series: [{ name: 'Floodwater mass (t)', x: tm, y: rnd(fm.map((v) => v / 1e3)) }, { name: 'Inflow (m³/min)', x: tm, y: rnd(fl.map((v) => v * 60)) }],
          explain: 'Water on board against time with the inflow rate. The inflow falls as the inside level rises towards the sea level; a curve that is still climbing at the end means the simulated time was too short.' },
        { type: 'line', title: 'Attitude during flooding', x: { label: 'Time', unit: 'min' }, y: { label: 'Angle', unit: 'deg' },
          series: [{ name: 'Heel', x: tm, y: rnd(he) }, { name: 'Trim', x: tm, y: rnd(tr) }, { name: 'Inclination', x: tm, y: rnd(ic) }, { name: 'Limit', x: [tm[0], tm[tm.length - 1]], y: [inp.heelMax, inp.heelMax], dash: true }],
          explain: 'Heel, trim and their combined inclination from a static equilibrium at each step. Read the time at which the inclination crosses the limit line: that is the response time available.' },
        { type: 'line', title: 'Draft, driving head and downflooding clearance', x: { label: 'Time', unit: 'min' }, y: { label: 'Height', unit: 'm' },
          series: [{ name: 'Draft', x: tm, y: rnd(dr) }, { name: 'Driving head at the opening', x: tm, y: rnd(hd) }, { name: 'Second opening above water', x: tm, y: rnd(df) }],
          explain: 'Sinkage of the unit, the head difference that drives the inflow (zero when equalised) and the clearance of the second opening. If the clearance reaches zero, flooding spreads beyond this compartment.' },
        p3,
      ],
      tables: [
        { title: 'Flooding history', columns: [{ key: 't', label: 'Time', unit: 'min' }, { key: 'm', label: 'Floodwater', unit: 't' }, { key: 'q', label: 'Inflow', unit: 'm³/h' }, { key: 'T', label: 'Draft', unit: 'm' }, { key: 'h', label: 'Heel', unit: 'deg' }, { key: 'tr', label: 'Trim', unit: 'deg' }, { key: 'df', label: 'Second opening', unit: 'm' }],
          rows: (() => { const [a, b, c, d2, e, f, g2] = thin([r.t, r.floodMass, r.flow, r.draft, r.heelDeg, r.trimDeg, dfH], 24); return a.map((t, i) => ({ t: round(t / 60, 4), m: round(b[i] / 1e3, 5), q: round(c[i] * 3600, 5), T: round(d2[i], 5), h: round(e[i], 4), tr: round(f[i], 4), df: round(g2[i], 4) })); })(),
          explain: 'Sampled states of the flooding: use a row as the intermediate stage of flooding for a stability check, or to time counter-ballasting.' },
      ],
      notes: [
        `Compartment of ${fmt(vComp)} m³ floodable volume breached by a ${fmt(inp.oArea)} m² opening ${fmt(h0)} m below the waterline: initial inflow ${fmt(r.flow[0] * 3600)} m³/h.`,
        r.sunk ? `No floating equilibrium was found after ${fmt(r.t[last] / 60)} min with ${fmt(r.floodMass[last] / 1e3)} t on board: the unit is taken as lost in this damage case.` : r.equalised ? `Flooding stops after ${fmt(tEq / 60)} min with ${fmt(r.floodMass[last] / 1e3)} t on board at ${fmt(r.heelDeg[last])}° heel and ${fmt(r.trimDeg[last])}° trim (positive = port up, bow down).` : `Flooding is still in progress after ${fmt(horizon / 60)} min (${fmt(r.floodMass[last] / 1e3)} t on board): extend the simulated duration.`,
        `Volume balance error ${fmt(balance)} m³ over ${n - 1} steps of ${fmt(inp.dt)} s. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', r.sunk ? 'fail' : sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: inclination up to ${fmt(inclMax)}° against ${fmt(inp.heelMax)}°, ${Number.isFinite(tCrit) ? `limit reached after ${fmt(tCrit / 60)} min` : 'no limit reached'}, peak inflow ${fmt(qMax * 3600)} m³/h against ${fmt(inp.pumpCap)} m³/h of pumps.`,
          r.sunk || inclMax > inp.heelMax ? `Subdivide the compartment: halving its length halves the floodwater moment; a floodable volume of about ${fmt(vComp * inp.heelMax / Math.max(inclMax, inp.heelMax))} m³ keeps the inclination near ${fmt(inp.heelMax)}°. Alternatively pre-plan counter-ballast into the diagonally opposite tank.`
            : Number.isFinite(tCrit) && tCrit / 60 < inp.tResp ? `The limit is reached in ${fmt(tCrit / 60)} min, sooner than the ${fmt(inp.tResp)} min response time: raise or make watertight the second opening, or fit remotely closed valves on this compartment.`
              : `The damage case is survivable with margin. Confirm the final condition with the damage stability module (residual righting lever and range), and keep the second opening at least ${fmt(Math.max(dfMin, 0))} m clear as computed.`,
          'IMO MODU Code 2009 §3.4 and §3.6; SOLAS 2020 Ch. II-1 Reg. 7-2 (intermediate flooding stages); DNV-OS-C301 (2021).'),
        adv('economic', 'info', `Dewatering ${fmt(r.floodMass[last] / 1e3)} t after the event takes ${fmt(inp.pumpCap > 0 ? r.floodMass[last] / rho / inp.pumpCap : 0)} h at ${fmt(inp.pumpCap)} m³/h, during which the unit is off station-keeping duty or production.`,
          'Weigh an extra watertight flat or bulkhead in the exposed columns against the downtime and salvage exposure of this damage case; subdivision is cheapest when decided before the compartment layout is frozen.',
          'Risk-based design: accidental limit state with a 10⁻⁴ annual event (NORSOK N-003 (2017) §9; DNV-OS-A101 (2021)); whole-life costing per ISO 15686-5.',
          `Dewatering against a 20 m head burns about ${fuelCo2Text(fuelT)} The larger footprint is the repair: a breached column plate replaced offshore needs a cofferdam and a support vessel for weeks.`),
      ],
    };
  },
  verify: [
    { // Box barge 100 × 20 m at 5 m draft; midship compartment 20 × 20 m over the full depth, permeability 1, opening 1 m² at the
      // keel, Cd = 0.6. Added weight: draft T = 5 + (400/2000)·w for an inside depth w, head = T − w = 5 − 0.8·w, so
      // 400·dw/dt = 0.6·1·√(2g·(5 − 0.8w)) → t = 400/(0.6·√(2g))·(2/0.8)·√5 = (400/2.6573)·5.5902 = 841.5 s = 14.02 min.
      // Final state: w = T = 6.25 m (same as the lost-buoyancy result T·L/(L − l)); floodwater 400·6.25·1.025 = 2562.5 t.
      name: 'Flooding of a midship compartment of a box barge through a keel opening',
      source: 'Torricelli orifice with added-weight sinkage, integrated by hand; final draft per Biran & López-Pulido (2014) ch. 11',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, cxMin: -10, cxMax: 10, cyMin: -10, cyMax: 10, czMin: 0, czMax: 12, perm: 1, ox: 0, oy: 0, oz: 0, oArea: 1, oCd: 0.6, dfX: 0, dfY: 0, dfZ: 12, dt: 5, tMaxMin: 30 },
      expect: [{ kpi: 'tEqualise', value: 14.02, tol: 0.02 }, { kpi: 'floodMass', value: 2562.5, tol: 2e-3 }, { kpi: 'draftEnd', value: 6.25, tol: 2e-3 }, { kpi: 'heelEnd', value: 0, abs: 1e-3 }, { kpi: 'volBalance', value: 0, abs: 1e-3 }, { kpi: 'dfMin', value: 5.75, tol: 2e-3 }] },
    { // Initial inflow of the default case: head = 22 − 19 = 3 m → Q = 0.6·0.3·√(2·9.80665·3) = 1.3808 m³/s = 4971 m³/h.
      name: 'Initial inflow equals the Torricelli value at the still-water head',
      source: 'Torricelli (orifice equation) — hand calculation',
      inputs: {},
      expect: [{ kpi: 'q0', value: 4971, tol: 2e-3 }] },
  ],
};
