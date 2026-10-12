// Damage stability by the lost-buoyancy method for a set of damage cases: final floating position, residual
// righting-lever curve, opening immersion and a probability-weighted survival index.
import { bodyToEarth, bodyIntegrals, solveEquilibrium } from '../../engine/hydro/hydrostatics.js';
import { gzCurve } from '../../engine/hydro/stability.js';
import { hullInputs, hullFrom, hullPlot, need, run, round, rnd, fmt, adv, check, checkMin, sev, worst, steelCo2Text } from './_common.js';
import { damagedEquilibrium, lostFrom, DEG } from './_damage.js';

const COLS = [{ key: 'name', label: 'Compartment', type: 'text' }, { key: 'dcase', label: 'Damage case', type: 'int' }, { key: 'prob', label: 'Case probability', type: 'number' },
  { key: 'xMin', label: 'x min', unit: 'm', type: 'number' }, { key: 'xMax', label: 'x max', unit: 'm', type: 'number' }, { key: 'yMin', label: 'y min', unit: 'm', type: 'number' }, { key: 'yMax', label: 'y max', unit: 'm', type: 'number' },
  { key: 'zMin', label: 'z min (above keel)', unit: 'm', type: 'number' }, { key: 'zMax', label: 'z max', unit: 'm', type: 'number' }, { key: 'perm', label: 'Permeability', type: 'number' }];

export default {
  id: 's11-damage-stability',
  suite: 11,
  title: 'Damage stability and compartment flooding',
  summary: 'Final floating position, residual righting levers and opening immersion of the unit after flooding of one or more compartments, for a list of damage cases with a probability-weighted survival index.',
  covers: [12],
  fidelity: 'numerical',
  standards: ['IMO MODU Code 2009 §3.4–3.6 (damage stability of column-stabilised and surface units)', 'SOLAS 1974 Ch. II-1 Part B-1 (probabilistic damage stability concept)', 'DNV-OS-C301 (2021) Ch.2 Sec.1', 'ABS MODU Rules Part 3 Ch.3 Sec.2'],
  models: ['Lost-buoyancy (constant-displacement) method', 'Nonlinear hydrostatic equilibrium of the damaged hull (Newton with relaxation fallback)', 'Residual righting-arm curve at constant displacement', 'Permeability-weighted compartment flooding', 'Probability-weighted survival index A = Σ pᵢ·sᵢ with sᵢ ∈ {0, 1}'],
  limits: 'Final-stage flooding only: compartments are open to the sea and fully equalised, with no intermediate stages, cross-flooding, cargo outflow or transient roll. Compartments are rectangular boxes in hull coordinates and must lie inside the hull (a column compartment is entered as a box of equal plan area). The residual curve is computed at the damaged trim held fixed. Survival factors are pass/fail, not the SOLAS s-factor formula, and wind heeling in the damaged condition is not applied.',
  searchTerms: 'damage stability lost buoyancy residual righting lever MODU flooding',
  inputs: [
    ...hullInputs({ type: 'semi' }),
    { key: 'compartments', label: 'Flooded compartments by damage case', type: 'table', columns: COLS, group: 'Damage cases', role: 'ic',
      default: [
        { name: 'Column C1, waterline compartment', dcase: 1, prob: 0.5, xMin: 28.35, xMax: 41.65, yMin: 28.35, yMax: 41.65, zMin: 17, zMax: 25, perm: 0.95 },
        { name: 'Pontoon P1 end tank', dcase: 2, prob: 0.35, xMin: 45, xMax: 55, yMin: 35, yMax: 43.5, zMin: 0, zMax: 9, perm: 0.95 },
        { name: 'Pontoon P1 end tank', dcase: 3, prob: 0.15, xMin: 45, xMax: 55, yMin: 35, yMax: 43.5, zMin: 0, zMax: 9, perm: 0.95 },
        { name: 'Pontoon P1 adjacent tank', dcase: 3, prob: 0.15, xMin: 45, xMax: 55, yMin: 26.5, yMax: 35, zMin: 0, zMax: 9, perm: 0.95 },
      ],
      help: 'One row per flooded space, grouped into damage cases by the case number; x forward of midship, y to port, z above keel. Permeability 0.95 for tanks and voids, 0.85 for machinery spaces, 0.60 for stores. The case probability is read from the first row of each case.' },
    { key: 'dfX', label: 'Opening x', unit: 'm', type: 'number', default: 30, min: -250, max: 250, group: 'Openings', role: 'bc', help: 'Position of the lowest opening through which progressive flooding could occur (vent, door, chain locker).' },
    { key: 'dfY', label: 'Opening y', unit: 'm', type: 'number', default: 30, min: -80, max: 80, group: 'Openings', role: 'bc', help: 'Transverse position (positive to port).' },
    { key: 'dfZ', label: 'Opening height above keel', unit: 'm', type: 'number', default: 44, min: 1, max: 250, group: 'Openings', role: 'bc', help: 'Height of the opening.' },
    { key: 'inclMax', label: 'Allowable inclination after damage', unit: 'deg', type: 'number', default: 17, min: 2, max: 35, group: 'Acceptance criteria', role: 'criteria', help: 'MODU Code column-stabilised units: 17° in any direction (with wind). Ships: 15° typical before cross-flooding. Confirm against the governing edition.' },
    { key: 'rangeMin', label: 'Required residual range of stability', unit: 'deg', type: 'number', default: 7, min: 0, max: 40, group: 'Acceptance criteria', role: 'criteria', help: 'Positive righting range beyond the damaged equilibrium (MODU Code: at least 7°).' },
    { key: 'openMargin', label: 'Required opening height above the damaged waterline', unit: 'm', type: 'number', default: 4, min: 0, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'MODU Code column-stabilised units: openings within 4 m above the final waterline must be watertight. Use 0 to require only that the opening stays dry.' },
    { key: 'gzMin', label: 'Required residual maximum righting lever', unit: 'm', type: 'number', default: 0.1, min: 0, max: 5, group: 'Acceptance criteria', role: 'criteria', help: 'Minimum residual GZ within the range (0.1 m is the SOLAS-type value).' },
    { key: 'indexReq', label: 'Required survival index', unit: '-', type: 'number', default: 0.9, min: 0, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'Minimum probability-weighted share of the listed damage cases that must be survived (1.0 = every case, the deterministic MODU approach).' },
    { key: 'dHeel', label: 'Heel step of the residual curve', unit: 'deg', type: 'number', default: 2.5, min: 1, max: 5, group: 'Numerical controls', role: 'control', help: 'Spacing of the residual righting-lever points.' },
  ],
  mesh: { key: 'dHeel', label: 'Heel step', levels: [5, 2.5, 1.25], qoi: 'rangeRes', isSize: true },
  compute(inp) {
    const hull = hullFrom(inp, { nTheta: 16 });
    need(hull.type !== 'tlp', 'Damage stability of a tension-leg platform is governed by tendon tension redistribution, not free-floating equilibrium: assess the free-floating (tow) condition with another hull type or use the tendon checks.');
    need(hull.hs.GMt > 0 && hull.hs.GMl > 0, `The intact unit is unstable (GM = ${fmt(Math.min(hull.hs.GMt, hull.hs.GMl))} m): fix the intact condition first.`);
    const rows = (inp.compartments || []).filter((r) => Number.isFinite(r.xMin));
    need(rows.length > 0, 'List at least one flooded compartment.');
    const ids = [...new Set(rows.map((r) => Math.round(r.dcase || 1)))].sort((a, b) => a - b);
    need(ids.length <= 12, 'At most 12 damage cases can be evaluated in one run.');
    const open = [inp.dfX, inp.dfY, inp.dfZ];
    const cases = ids.map((id) => {
      const mine = rows.filter((r) => Math.round(r.dcase || 1) === id);
      const comps = mine.map((r) => {
        need(r.xMax > r.xMin && r.yMax > r.yMin && r.zMax > r.zMin, `Compartment "${r.name}": each max coordinate must exceed its min.`);
        need(r.zMin >= -1e-6 && r.zMax <= hull.depth + 1e-6, `Compartment "${r.name}" extends outside the hull depth (0…${fmt(hull.depth)} m).`);
        return { name: r.name, box: { xMin: r.xMin, xMax: r.xMax, yMin: r.yMin, yMax: r.yMax, zMin: r.zMin, zMax: r.zMax }, permeability: Math.min(Math.max(r.perm ?? 0.95, 0), 1) };
      });
      const lost = lostFrom(comps);
      // engine Newton first; it stalls or finds an unstable root when a flooded space crosses the waterline, so the result is
      // accepted only if the relaxation solver agrees
      const eq = damagedEquilibrium({ mesh: hull.mesh, lost, mass: hull.mass, cog: hull.cog, rho: hull.rho, gmT: hull.hs.GMt, gmL: hull.hs.GMl, draft0: hull.T });
      let viaEngine = false;
      try { const e = solveEquilibrium({ mesh: hull.mesh, mass: hull.mass, cog: hull.cog, rho: hull.rho, lost, maxIter: 25 }); viaEngine = e.converged && Math.abs(e.heel - eq.heel) < 2e-3 && Math.abs(e.trim - eq.trim) < 2e-3; } catch { viaEngine = false; }
      const prob = Math.max(mine[0].prob ?? 1, 0);
      if (eq.capsized || !eq.converged) return { id, names: mine.map((r) => r.name), prob, capsized: true, uc: 99, incl: 90, heel: eq.heel, trim: eq.trim, draft: eq.draft, range: 0, gzMax: 0, openH: -99, lostV: 0, viaEngine, curve: null, checks: [] };
      const side = eq.heel < -1e-9 ? -1 : 1;
      const heq = Math.abs(eq.heel) / DEG;
      // the curve starts one step below the damaged equilibrium: the part from upright to there is not needed for the criteria
      const h0 = Math.max(0, Math.floor(heq / inp.dHeel - 1) * inp.dHeel);
      const heels = []; for (let h = h0; h <= Math.min(h0 + 40, 80) + 1e-9; h += inp.dHeel) heels.push(round(h, 6));
      const curve = run(() => gzCurve({ mesh: hull.mesh, mass: hull.mass, cog: hull.cog, rho: hull.rho, heelsDeg: heels, side, freeTrim: false, trim: eq.trim, lost, openings: [{ name: 'opening', point: open }] }), `Residual righting levers, case ${id}`);
      let vanish = heels[heels.length - 1];
      for (let i = 1; i < heels.length; i++) if (heels[i] > heq + 1e-6 && curve.gz[i] <= 0 && curve.gz[i - 1] > 0) { vanish = heels[i - 1] + (heels[i] - heels[i - 1]) * curve.gz[i - 1] / (curve.gz[i - 1] - curve.gz[i]); break; }
      const end = Math.min(vanish, curve.downfloodDeg ?? Infinity);
      const range = Math.max(0, end - heq);
      const gzRes = Math.max(0, ...curve.gz.filter((_, i) => curve.heelDeg[i] >= heq && curve.heelDeg[i] <= end));
      const incl = Math.acos(Math.cos(eq.heel) * Math.cos(eq.trim)) / DEG;
      const openH = bodyToEarth(open, eq.draft, eq.heel, eq.trim)[2];
      let lostV = 0; for (const c of lost) lostV += c.permeability * bodyIntegrals(c.mesh, null, eq.draft, eq.heel, eq.trim).V;
      const checks = [
        check(`Case ${id}: inclination after damage`, incl, inp.inclMax, 'deg', 'IMO MODU Code 2009 §3.4.3.1 (column-stabilised: 17°); §3.4.1 for surface units', 'Guards against a final list so large that evacuation, ballast-system operation and the remaining stability are compromised.'),
        checkMin(`Case ${id}: residual range of stability`, range, Math.max(inp.rangeMin, 1e-6), 'deg', 'IMO MODU Code 2009 §3.4.3.3 (range of at least 7°)', 'Guards the damaged unit against capsizing in wind and waves: a short residual range means a small extra heeling moment overturns it.'),
        checkMin(`Case ${id}: residual maximum righting lever`, gzRes, Math.max(inp.gzMin, 1e-6), 'm', 'SOLAS Ch. II-1 Reg. 7-2 concept (residual GZ ≥ 0.1 m); project criterion for offshore units', 'Guards against a damaged unit with almost no righting moment left, which cannot resist wind or shifting loads.'),
        check(`Case ${id}: opening clearance above the damaged waterline`, inp.openMargin, Math.max(openH, 1e-6), 'm', 'IMO MODU Code 2009 §3.4.3.2 (openings below the final waterline, and within 4 m above it, watertight)', 'Guards against progressive flooding through vents and doors, which turns a survivable damage into a loss.'),
      ];
      return { id, names: mine.map((r) => r.name), prob, capsized: false, uc: worst(checks).uc, incl, heel: eq.heel, trim: eq.trim, draft: eq.draft, range, gzMax: gzRes, openH, lostV, viaEngine, curve, checks, heq };
    });
    const pTot = cases.reduce((a, c) => a + c.prob, 0) || 1;
    const index = cases.reduce((a, c) => a + (c.uc <= 1 ? c.prob : 0), 0) / pTot;
    const g = cases.reduce((a, c) => (c.uc > a.uc ? c : a), cases[0]);
    const checks = cases.flatMap((c) => (c.capsized ? [check(`Case ${c.id}: floating equilibrium exists`, 99, 1, '-', 'IMO MODU Code 2009 §3.4', 'No stable floating position was found within 60°: the unit capsizes or sinks in this damage case.')] : c.checks));
    checks.push(checkMin('Survival index over the damage cases', index, Math.max(inp.indexReq, 1e-6), '-', 'SOLAS Ch. II-1 Part B-1 (attained index A ≥ required index R) — simplified pass/fail survival factors', 'Probability-weighted share of the damage cases survived; below the requirement the subdivision is insufficient for the assumed damage statistics.'));
    const gc = g.curve;
    return {
      kpis: [
        { key: 'draftDamaged', label: 'Draft after damage (governing case)', value: round(g.draft, 5), unit: 'm', explain: 'Draft at the keel origin in the final flooded condition: the unit sinks until the remaining intact volume again displaces its weight.' },
        { key: 'heelDamaged', label: 'Heel after damage (magnitude)', value: round(Math.abs(g.heel) / DEG, 4), unit: 'deg', explain: 'Final list about the longitudinal axis towards the damaged side.' },
        { key: 'trimDamaged', label: 'Trim after damage (magnitude)', value: round(Math.abs(g.trim) / DEG, 4), unit: 'deg', explain: 'Final inclination about the transverse axis.' },
        { key: 'inclination', label: 'Total inclination after damage', value: round(g.incl, 4), unit: 'deg', limit: inp.inclMax, sense: 'max', explain: 'Angle of the deck from horizontal combining heel and trim; the MODU Code limit for column-stabilised units is 17°.' },
        { key: 'lostVolume', label: 'Lost buoyancy volume', value: round(g.lostV, 5), unit: 'm³', explain: 'Volume of the flooded spaces below the final waterline times permeability: the buoyancy the rest of the hull has to replace.' },
        { key: 'rangeRes', label: 'Residual range of stability', value: round(g.range, 4), unit: 'deg', limit: inp.rangeMin, sense: 'min', explain: 'Heel from the damaged equilibrium to the angle where the righting lever vanishes or an opening immerses.' },
        { key: 'gzRes', label: 'Residual maximum righting lever', value: round(g.gzMax, 4), unit: 'm', limit: inp.gzMin, sense: 'min', explain: 'Largest righting lever available beyond the damaged equilibrium.' },
        { key: 'openHeight', label: 'Opening height above the damaged waterline', value: round(g.openH, 4), unit: 'm', limit: inp.openMargin, sense: 'min', explain: 'Vertical clearance of the listed opening in the final condition. Negative means it is under water and progressive flooding follows.' },
        { key: 'index', label: 'Survival index', value: round(index, 4), unit: '-', limit: inp.indexReq, sense: 'min', explain: 'Sum of the probabilities of the damage cases that meet every criterion, divided by the total probability listed.' },
        { key: 'nFail', label: 'Damage cases not survived', value: cases.filter((c) => c.uc > 1).length, unit: '-', explain: 'Number of listed cases failing at least one criterion.' },
        { key: 'ucMax', label: 'Governing utilisation', value: round(Math.min(g.uc, 99), 4), unit: '-', limit: 1, sense: 'max', explain: 'Largest criterion utilisation over all damage cases.' },
      ],
      checks,
      plots: [
        ...(gc ? [{ type: 'line', title: `Residual righting-lever curve, case ${g.id}`, x: { label: 'Heel towards the damaged side', unit: 'deg' }, y: { label: 'GZ', unit: 'm' },
          series: [{ name: 'Residual GZ', x: gc.heelDeg, y: rnd(gc.gz) }, { name: 'Damaged equilibrium', x: [round(g.heq), round(g.heq)], y: [Math.min(0, ...gc.gz.map((v) => round(v))), round(Math.max(g.gzMax, 0.1))], dash: true },
            ...(gc.downfloodDeg != null ? [{ name: 'Opening immerses', x: [round(gc.downfloodDeg), round(gc.downfloodDeg)], y: [0, round(Math.max(g.gzMax, 0.1))], dash: true }] : [])],
          explain: 'Righting lever of the flooded unit heeling towards the damage. It is negative below the damaged equilibrium (the curve starts just before it), then positive; the width and height of the positive lobe are the residual range and lever the criteria ask for.' }] : []),
        { type: 'bar', title: 'Damage cases: inclination and utilisation', categories: cases.map((c) => `Case ${c.id}`), y: { label: 'Inclination (deg) / utilisation', unit: '-' },
          series: [{ name: 'Inclination (deg)', y: cases.map((c) => round(c.incl, 4)) }, { name: 'Governing utilisation', y: cases.map((c) => round(Math.min(c.uc, 5), 4)) }],
          explain: 'Final inclination and worst criterion utilisation for each damage case. Cases with utilisation above 1.0 are not survived and count as zero in the index.' },
        hullPlot(hull.mesh, g.draft, 'Hull at the damaged draft', 'Intact hull surface coloured by pressure head at the damaged draft (upright view). Compare the waterline with the compartment heights entered to see which spaces straddle it.'),
      ],
      tables: [
        { title: 'Damage case results', columns: [{ key: 'id', label: 'Case' }, { key: 'sp', label: 'Flooded spaces' }, { key: 'p', label: 'Probability' }, { key: 'd', label: 'Draft', unit: 'm' }, { key: 'h', label: 'Heel', unit: 'deg' }, { key: 't', label: 'Trim', unit: 'deg' }, { key: 'v', label: 'Lost volume', unit: 'm³' }, { key: 'r', label: 'Residual range', unit: 'deg' }, { key: 'gz', label: 'Residual GZ', unit: 'm' }, { key: 'o', label: 'Opening height', unit: 'm' }, { key: 's', label: 'Survives' }],
          rows: cases.map((c) => ({ id: c.id, sp: c.names.join(' + '), p: round(c.prob / pTot, 4), d: round(c.draft, 4), h: round(c.heel / DEG, 4), t: round(c.trim / DEG, 4), v: round(c.lostV, 5), r: round(c.range, 4), gz: round(c.gzMax, 4), o: round(c.openH, 4), s: c.capsized ? 'capsizes' : c.uc <= 1 ? 'yes' : 'no' })),
          explain: 'One row per damage case: final floating position (heel positive = port side up, trim positive = bow down), residual stability and whether every criterion is met.' },
      ],
      notes: [
        `${cases.length} damage case(s); governing case ${g.id} (${g.names.join(' + ')}): ${g.capsized ? 'no stable equilibrium' : `draft ${fmt(g.draft)} m, heel ${fmt(g.heel / DEG)}°, trim ${fmt(g.trim / DEG)}°, inclination ${fmt(g.incl)}°`}.`,
        `Survival index ${fmt(index)} against ${fmt(inp.indexReq)} required.`,
        cases.some((c) => !c.viaEngine && !c.capsized) ? 'In at least one case the flooded space crosses the waterline, where the restoring moment is piecewise: the equilibrium was found by the relaxation solver (the Newton solver alone stalls or returns an unstable position there).' : 'Newton and relaxation solvers agree on every equilibrium.',
      ],
      advice: [
        adv('engineering', sev(g.uc), g.capsized ? `Case ${g.id} has no stable floating position.` : `Case ${g.id} governs: inclination ${fmt(g.incl)}° (limit ${fmt(inp.inclMax)}°), residual range ${fmt(g.range)}°, opening ${fmt(g.openH)} m above the damaged waterline.`,
          g.uc > 1 ? `Reduce the lost buoyancy or its lever: add a watertight flat or bulkhead to cut the largest flooded space (${fmt(g.lostV)} m³ lost now), lower KG to regain residual range, or raise and close the opening at z = ${fmt(inp.dfZ)} m.`
            : 'All listed cases are survived. Extend the list to every compartment within the damage extent and penetration of the governing code, add the wind heeling lever in the damaged condition, and check intermediate flooding stages with the progressive-flooding module.',
          'IMO MODU Code 2009 §3.4 (damage stability), §3.5 (extent of damage), §3.6 (watertight integrity).'),
        adv('economic', 'info', `Subdivision is cheap at design and almost impossible to add later: one extra watertight flat in a ${fmt(hull.type === 'semi' ? inp.Dcol : hull.Bc)} m column or tank is of the order of ${fmt(hull.type === 'semi' ? Math.PI * inp.Dcol ** 2 / 4 * 0.012 * 7.85 * 1.5 : hull.Bc * 20 * 0.012 * 7.85 * 1.5)} t of steel.`,
          'Trade extra subdivision against allowable KG: each damage case that fails forces a lower operating KG, which is paid for as lost deck payload over the whole life of the unit.',
          'Whole-life costing per ISO 15686-5; weight control per ISO 19901-5 (2021).',
          `${steelCo2Text(hull.type === 'semi' ? Math.PI * inp.Dcol ** 2 / 4 * 0.012 * 7.85 * 1.5 : hull.Bc * 20 * 0.012 * 7.85 * 1.5)} per added flat — small against the environmental consequence of losing the unit and its inventory.`),
      ],
    };
  },
  verify: [
    { // Box 100 × 20 m at 5 m draft; a full-breadth, full-depth midship compartment 20 m long is flooded (permeability 1):
      // intact waterplane length falls to 80 m → T′ = T·L/(L − l) = 5·100/80 = 6.25 m, no heel or trim; lost volume = 20·20·6.25 = 2500 m³.
      name: 'Symmetric midship flooding of a box barge: T′ = T·L/(L − l)',
      source: 'Biran & López-Pulido, Ship Hydrostatics and Stability (2014) ch. 11 — lost-buoyancy method, worked by hand',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, dfX: 0, dfY: 10, dfZ: 12, openMargin: 0,
        compartments: [{ name: 'midship hold', dcase: 1, prob: 1, xMin: -10, xMax: 10, yMin: -10, yMax: 10, zMin: 0, zMax: 12, perm: 1 }] },
      expect: [{ kpi: 'draftDamaged', value: 6.25, tol: 1e-4 }, { kpi: 'heelDamaged', value: 0, abs: 1e-3 }, { kpi: 'lostVolume', value: 2500, tol: 1e-3 }, { kpi: 'openHeight', value: 5.75, tol: 1e-3 }] },
    { // Same barge, a port-side tank 20 × 10 m over the full depth at 50 % permeability: the effective waterplane falls from 2000
      // to 1900 m², so the parallel sinkage gives T′ ≈ 10 000/1900 = 5.263 m at the flotation centre (the small heel adds < 1 %).
      name: 'Asymmetric flooding with permeability: parallel sinkage ∇/(Aw − μ·a)',
      source: 'Lost-buoyancy method, parallel sinkage (Biran & López-Pulido 2014 ch. 11) — hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, dfX: 0, dfY: 0, dfZ: 12, openMargin: 0,
        compartments: [{ name: 'port tank', dcase: 1, prob: 1, xMin: -10, xMax: 10, yMin: 0, yMax: 10, zMin: 0, zMax: 12, perm: 0.5 }] },
      expect: [{ kpi: 'draftDamaged', value: 5.263, tol: 0.01 }, { kpi: 'index', value: 1, tol: 1e-9 }] },
  ],
};
