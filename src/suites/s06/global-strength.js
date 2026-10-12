// Global strength of a pile-supported dolphin: rigid cap on nonlinear pile-head stiffnesses, pile forces and
// ISO 19902 tubular member checks.
import { DEG, need, fmt, round, cap, adv, check, ucCheck, arr, pileInputs, soilInputs, layoutInputs, capInputs, capLoadInputs, capLoad, pileProps, profileFrom, analyseDolphin, pileChecks, pileSteelT, steelCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS, CO2_CONCRETE } from './_common.js';

export default {
  id: 's06-global-strength', suite: 6, title: 'Global structural strength',
  summary: 'Cap displacement, pile axial forces, shears, bending moments and tubular member utilisation of a dolphin under berthing or mooring load with its self-weight.',
  covers: [3], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §13 (tubular members)', 'BS 6349-2:2019 §6–7', 'API RP 2GEO (2011, add. 2014) §8.5', 'EN 1993-5:2007 (steel piling)'],
  models: ['Static equilibrium of a rigid pile cap (matrix stiffness method)', 'Beam on nonlinear Winkler foundation (p–y) for the pile-head stiffness', 'Secant stiffness iteration', 'Beam-column interaction of tubular members', 'Raked-pile push–pull resolution'],
  limits: 'Rigid cap; every pile takes the secant head stiffness of the most loaded pile (updated twice), so shadowing between rows enters only through the p-multiplier. Bending below the cap is taken from the governing single-pile solution scaled by pile shear. The cap itself, pile-to-cap connection and local punching are not designed here. ' + NOT_HIFI,
  publishes: ['dolphin.pileAxialMax', 'dolphin.pileTensionMax', 'dolphin.pileShearMax', 'dolphin.pileMomentMax', 'dolphin.Kh', 'dolphin.capDisp', 'dolphin.stressPerLoad', 'dolphin.pileD', 'dolphin.pileT', 'dolphin.Lfree', 'dolphin.Lemb', 'dolphin.capWeight', 'quantities.steel_t', 'quantities.concrete_m3'],
  searchTerms: 'breasting dolphin pile group design berthing load raked piles',
  inputs: [
    ...layoutInputs(),
    ...pileInputs(),
    ...capInputs(),
    ...soilInputs(),
    ...capLoadInputs(),
    { key: 'fsAxial', label: 'Factor of safety on axial pile capacity', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Design factors', role: 'criteria', help: `API RP 2A-WSD: 1.5 extreme, 2.0 operating (on unfactored loads). ${CONFIRM}` },
    { key: 'nElem', label: 'Elements over the embedded length', unit: '-', type: 'int', default: 40, min: 10, max: 320, group: 'Numerical', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [20, 40, 80, 160], qoi: 'uCap' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const profile = profileFrom(inp, P.D);
    const load = capLoad(inp);
    need(Math.hypot(load.Fx, load.Fy) > 0 || Math.abs(load.Fz) > 0, 'No load is applied to the dolphin.');
    const res = analyseDolphin(inp, P, profile, load);
    const rows = pileChecks(res, P);
    const gov = rows.reduce((a, r) => (r.uc > a.uc ? r : a), rows[0]);
    const Nmax = Math.max(...rows.map((r) => r.axial)), Nmin = Math.min(...rows.map((r) => r.axial));
    const Vmax = Math.max(...rows.map((r) => r.shear)), Mmax = Math.max(...rows.map((r) => r.moment));
    const H = Math.hypot(load.Fx, load.Fy);
    const checks = [
      ucCheck(`Tubular member check, pile ${gov.id} (${gov.gov || 'combined'})`, gov.uc, gov.ref, 'Combined axial force and bending of the most loaded pile against the ISO 19902 strength and beam-column stability expressions. Above 1 the pile yields or buckles as a beam-column and the dolphin loses stiffness.'),
      check('Pile axial compression against soil capacity', Math.max(Nmax, 0) / 1e3, res.axCap / inp.fsAxial / 1e3, 'kN', 'API RP 2GEO (2011) §8.1; API RP 2A-WSD §6.3.4', 'Largest pile compression against shaft friction plus end bearing divided by the factor of safety. If exceeded the pile settles and sheds load to its neighbours.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const steelT = pileSteelT(inp, res.n), conc = inp.capB * inp.capL * inp.capT;
    // wall thickness that brings the member check to 1.0 (bending-governed: uc ∝ 1/t to first order)
    const tReq = inp.corr + (inp.t - inp.corr) * gov.uc;
    const dSteel = pileSteelT({ ...inp, t: tReq }, res.n) - steelT;
    const u = res.disp.u, dir = u > 0 && res.grp ? [res.disp.ux / u, res.disp.uy / u] : [0, 1];
    // 3-D view: each pile sampled along its length, deformed by the governing pile shape scaled by its shear share.
    const sol = res.sol, nz = sol.z.length, step = Math.max(1, Math.floor(nz / 14));
    const nodes = [], lines = [], vals = [], def = [];
    const yMax = Math.max(...Array.from(sol.y, Math.abs), 1e-9), scale = 0.1 * (inp.Lfree + inp.Lemb) / yMax;
    res.piles.forEach((p, k) => {
      const b = p.batter * DEG, a = p.azimuth * DEG, share = res.rows[k].shear / Math.max(Vmax, 1);
      let prev = -1;
      for (let i = 0; i < nz; i += step) {
        const s = sol.z[i] + inp.Lfree; // distance down from the head
        const x = p.x + Math.sin(b) * Math.cos(a) * s, y = p.y + Math.sin(b) * Math.sin(a) * s, z = -Math.cos(b) * s;
        nodes.push([round(x), round(y), round(z)]); vals.push(round(rows[k].uc, 3));
        def.push([round(x + dir[0] * sol.y[i] * share * scale), round(y + dir[1] * sol.y[i] * share * scale), round(z)]);
        if (prev >= 0) lines.push([prev, nodes.length - 1]);
        prev = nodes.length - 1;
      }
    });
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `${w.name}: utilisation ${fmt(w.uc)}. Cap displacement ${fmt(u * 1e3)} mm under the factored load of ${fmt(H / 1e3)} kN; pile forces ${fmt(Nmin / 1e3)} to ${fmt(Nmax / 1e3)} kN axial, ${fmt(Mmax / 1e3)} kN·m bending.`,
        ok ? `Reserve of ${fmt((1 - w.uc) * 100)} % on the governing check. Confirm with the cyclic p–y curves and the lower-bound soil strength, and check the pile-to-cap connection for ${fmt(gov.headMoment / 1e3)} kN·m.`
          : (w.name.startsWith('Tubular') ? `Increase the wall from ${fmt(inp.t)} mm to about ${fmt(tReq)} mm (bending-governed estimate), raise the steel grade, or add rake so that the load is carried in push–pull.` : w.name.startsWith('Pile axial') ? `Lengthen the piles or widen the pile spacing: the push–pull lever arm is the cheapest way to cut axial force.` : 'Increase pile diameter or embedment, or reduce the load per pile by adding piles.'),
        'ISO 19902:2020 §13; API RP 2GEO (2011) §8; BS 6349-2:2019',
        ok ? `The piles hold ${steelCo2Text(steelT)}` : `The thicker wall adds ${steelCo2Text(Math.max(dSteel, 0))}`),
      adv('economic', ok ? 'info' : 'warn', `Quantities: ${fmt(steelT)} t of pile steel in ${res.n} pile(s) and ${fmt(conc)} m³ of cap concrete (≈ ${fmt(conc * CO2_CONCRETE)} t CO₂e at the screening factor ${CO2_CONCRETE} t/m³).`,
        ok && w.uc < 0.7 ? `Utilisation is ${fmt(w.uc)}: a wall of about ${fmt(Math.max(tReq, inp.corr + 12))} mm would still pass, saving roughly ${fmt(Math.max(-dSteel, 0))} t of steel; check drivability and fatigue before reducing.` : 'Pile tonnage and the number of pile-driving set-ups drive the cost: compare fewer, larger piles against more, lighter ones at equal utilisation.',
        ECON_BASIS),
    ];
    return {
      kpis: [
        { key: 'ucMax', label: 'Governing member utilisation', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'Largest ISO 19902 utilisation over all piles (strength and beam-column stability, with the partial resistance factors). Above 1.0 the pile section is inadequate.' },
        { key: 'uCap', label: 'Cap horizontal displacement', unit: 'mm', value: u * 1e3, explain: 'Resultant horizontal movement of the pile cap under the factored load, including soil nonlinearity. Compare with the fender deflection: a stiff dolphin should move much less than the fender.' },
        { key: 'rotCap', label: 'Cap rotation (tilt)', unit: 'deg', value: res.disp.rot / DEG, explain: 'Rotation of the cap about a horizontal axis (pile-head rotation for a single pile). It changes the fender face angle and the walkway seats.' },
        { key: 'twist', label: 'Cap twist in plan', unit: 'deg', value: Math.abs(res.disp.rz) / DEG, explain: 'Rotation about the vertical axis from eccentric load, resisted by the lateral stiffness of the outer piles.' },
        { key: 'Kgroup', label: 'Dolphin lateral stiffness (secant)', unit: 'MN/m', value: res.Kgroup / 1e6, explain: 'Applied horizontal load over cap displacement at this load level. It falls as the soil yields; use it to couple the dolphin to berthing and mooring models.' },
        { key: 'Nmax', label: 'Largest pile compression', unit: 'kN', value: Nmax / 1e3, explain: 'Greatest axial compression at a pile head (weight plus push–pull from overturning). Governs axial bearing.' },
        { key: 'Nmin', label: 'Smallest pile axial force (tension −)', unit: 'kN', value: Nmin / 1e3, explain: 'Least axial force at a pile head; a negative value is tension, to be checked against pull-out capacity and the head connection.' },
        { key: 'Vmax', label: 'Largest pile shear', unit: 'kN', value: Vmax / 1e3, explain: 'Greatest lateral force at a pile head: the load each pile passes to the soil.' },
        { key: 'Mmax', label: 'Largest pile bending moment', unit: 'kN·m', value: Mmax / 1e3, explain: 'Greatest bending moment in any pile: at the cap for fixed heads, below the mudline for pinned heads.' },
        { key: 'sigmaMax', label: 'Largest combined stress', unit: 'MPa', value: gov.sigma / 1e6, limit: P.fy / 1e6, sense: 'max', explain: 'Axial plus bending stress in the governing pile under factored load, against yield. An elastic indicator; the code check above also covers buckling.' },
        { key: 'mobil', label: 'Soil resistance mobilised', unit: '-', value: cap(res.mobil, 1.5), explain: 'Peak ratio of lateral soil reaction to its ultimate value along the governing pile. Above about 0.7 near the surface the response is strongly nonlinear.' },
        { key: 'Le', label: 'Equivalent pile length to fixity', unit: 'm', value: res.Le, explain: 'Length of the cantilever (fixed-head: fixed–guided column) with the same lateral head stiffness as the pile in soil. Used as the buckling length.' },
        { key: 'steel', label: 'Pile steel', unit: 't', value: steelT, explain: 'Mass of all piles at the as-built wall thickness over free plus embedded length.' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'Dolphin piles: utilisation and deflected shape', nodes, lines, nodeValues: vals, valueLabel: 'Member utilisation', valueUnit: '-', deformed: def,
          explain: `Piles from head to toe coloured by member utilisation, with the deflected shape exaggerated ${fmt(scale)} times. Look for the pile carrying the most push–pull and for the depth at which the deflection dies out.` },
        { type: 'line', title: 'Bending moment and deflection of the governing pile', x: { label: 'Moment (kN·m) / deflection (mm)', unit: '' }, y: { label: 'Depth below mudline', unit: 'm' },
          series: [{ name: 'Moment (kN·m)', x: arr(Array.from(sol.moment, (m) => m / 1e3)), y: arr(Array.from(sol.z, (z) => -z)) }, { name: 'Deflection (mm)', x: arr(Array.from(sol.y, (y) => y * 1e3)), y: arr(Array.from(sol.z, (z) => -z)) }],
          explain: 'Moment and lateral deflection along the most loaded pile (negative depth is below the mudline). The moment peak fixes the wall thickness; the first zero of deflection shows how much embedment is working.' },
        { type: 'bar', title: 'Pile-head forces', categories: rows.map((r) => `Pile ${r.id}`), series: [{ name: 'Axial (kN, compression +)', y: rows.map((r) => round(r.axial / 1e3)) }, { name: 'Shear (kN)', y: rows.map((r) => round(r.shear / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Axial force and shear at each pile head. Opposite signs of axial force across the group show the push–pull couple resisting overturning.' },
      ],
      tables: [{ title: 'Pile forces and utilisation', columns: [{ key: 'id', label: 'Pile' }, { key: 'N', label: 'Axial (compression +)', unit: 'kN' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'Mh', label: 'Head moment', unit: 'kN·m' }, { key: 'M', label: 'Design moment', unit: 'kN·m' }, { key: 'T', label: 'Torsion', unit: 'kN·m' }, { key: 's', label: 'Combined stress', unit: 'MPa' }, { key: 'uc', label: 'Utilisation' }, { key: 'gov', label: 'Governing check' }],
        rows: rows.map((r) => ({ id: r.id, N: round(r.axial / 1e3), V: round(r.shear / 1e3), Mh: round(r.headMoment / 1e3), M: round(r.moment / 1e3), T: round(r.torsion / 1e3), s: round(r.sigma / 1e6), uc: round(r.uc), gov: r.gov })),
        explain: 'One row per pile under the factored load: head forces from the rigid-cap analysis, the design moment along the pile and the ISO 19902 utilisation with the expression that governs.' }],
      notes: [`${w.name} governs.`, `Pile heads ${res.single ? 'single pile: free head with the cap moment applied' : res.fixed ? 'moment-fixed into the cap' : 'pinned'}; equivalent length to fixity ${fmt(res.Le)} m with K = ${res.Kfac}.`,
        ...(gov.warnings?.length ? gov.warnings.slice(0, 2) : []), ...(res.grp ? [`Equilibrium error of the group solution: ${res.grp.equilibriumError.toExponential(1)}.`] : [])],
      advice,
      publish: { 'dolphin.pileAxialMax': Math.max(Nmax, 0), 'dolphin.pileTensionMax': Math.max(-Nmin, 0), 'dolphin.pileShearMax': Vmax, 'dolphin.pileMomentMax': Mmax, 'dolphin.Kh': res.Kgroup, 'dolphin.capDisp': u,
        'dolphin.stressPerLoad': H > 0 ? (Mmax / P.Z) / H : 0, 'dolphin.pileD': inp.D, 'dolphin.pileT': inp.t / 1e3, 'dolphin.Lfree': inp.Lfree, 'dolphin.Lemb': inp.Lemb, 'dolphin.capWeight': load.W, 'quantities.steel_t': steelT, 'quantities.concrete_m3': conc },
      quantities: { steel_t: round(steelT), concrete_m3: round(conc), vessel_days: { piling: round(res.n * 0.75 + 2, 3) } },
    };
  },
  verify: [
    // Pinned heads, four vertical piles at y = ±2 m: the overturning moment F·z = 1000·2 = 2000 kN·m is carried by push–pull,
    // ΔN = M·y/Σy² = 2000·2/16 = 250 kN on top of W/4 = 2000/4 = 500 kN → 750 and 250 kN; shear 1000/4 = 250 kN per pile.
    { name: 'Pinned four-pile group: N = W/n ± M·y/Σy²', source: 'Rigid-cap statics (Poulos & Davis 1980, ch. 9): moment equilibrium of a symmetric pinned-head group',
      inputs: { headFix: 'pinned', Fy: 1000, Fx: 0, zLoad: 2, ecc: 0, Mz: 0, Fup: 0, gammaF: 1, gammaG: 1, capB: 4, capL: 4, capT: 5, gammaConc: 25, Wequip: 0 },
      expect: [{ kpi: 'Nmax', value: 750, tol: 1e-4 }, { kpi: 'Nmin', value: 250, tol: 1e-4 }, { kpi: 'Vmax', value: 250, tol: 1e-4 }] },
  ],
};
