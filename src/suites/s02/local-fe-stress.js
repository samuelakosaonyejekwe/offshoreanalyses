// Local finite-element stress analysis of a deck or web plate with a circular opening (penetration, manhole,
// cope hole) under the membrane stress taken from the global frame analysis: plane-stress Q4 mesh with a
// refinement ladder, stress concentration, yield check and discretisation error indicator.
import { plateWithHoleMesh, solvePlane, zzError } from '../../engine/structural/fem2d.js';
import { steelOf, need, run, fmt, round, adv, check, meshPlot, steelCo2Text, CONFIRM, linspace } from './_common.js';

export default {
  id: 's02-local-fe-stress', suite: 2, title: 'Local finite-element stress at a plate opening',
  summary: 'Peak stress, stress concentration factor and yield utilisation around a circular opening in a deck or girder-web plate under the global membrane stresses, with mesh convergence.',
  covers: [2], fidelity: 'numerical',
  standards: ['DNV-RP-C208 (2022) §4–§5 (FE stress analysis and mesh requirements)', 'DNV-RP-C203 (2021) §4.3, App. C (hot-spot stress and SCF at openings)', 'ISO 19902:2020 §12 (structural modelling)', 'EN 1993-1-5:2006 Annex C (FE methods)'],
  models: ['Cauchy elasticity', 'Generalised Hooke’s law (plane stress)', 'Q4 incompatible-mode finite elements (QM6)', 'von Mises yielding', 'Tresca yield theory', 'Zienkiewicz–Zhu error indicator', 'Weld stress concentration models'],
  limits: 'Linear elastic plane-stress quarter model of a doubly symmetric plate with a central circular hole under uniform edge stresses; no plate bending, no edge reinforcement, no weld geometry and no plasticity. Peak stress above yield means local yielding, which needs a nonlinear analysis outside this module.',
  searchTerms: 'stress concentration circular hole finite width plate finite element mesh convergence deck opening',
  inputs: [
    { key: 'W', label: 'Plate width (across the load, y)', unit: 'mm', type: 'number', default: 2000, min: 50, max: 2e4, group: 'Plate and opening', role: 'geometry', help: 'Full width between stiffeners or flanges.' },
    { key: 'H', label: 'Plate length (along the load, x)', unit: 'mm', type: 'number', default: 4000, min: 50, max: 4e4, group: 'Plate and opening', role: 'geometry' },
    { key: 't', label: 'Plate thickness', unit: 'mm', type: 'number', default: 12, min: 2, max: 150, group: 'Plate and opening', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'd', label: 'Opening diameter', unit: 'mm', type: 'number', default: 600, min: 5, max: 1e4, group: 'Plate and opening', role: 'geometry', design: true },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Plate and opening', role: 'material' },
    { key: 'sx', label: 'Nominal membrane stress along x', unit: 'MPa', type: 'number', default: 90, min: -600, max: 600, group: 'Boundary stresses', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Gross-section stress on the far edges, from the global frame analysis (flange or deck-plate stress). Tension positive.' },
    { key: 'sy', label: 'Nominal membrane stress along y', unit: 'MPa', type: 'number', default: 20, min: -600, max: 600, group: 'Boundary stresses', role: 'bc', help: 'Transverse gross stress (biaxial deck plate); 0 for a uniaxial strip.' },
    { key: 'dSx', label: 'Nominal stress range along x (fatigue)', unit: 'MPa', type: 'number', default: 25, min: 0, max: 600, group: 'Boundary stresses', role: 'load', help: 'Cyclic part of the x stress; the hot-spot range at the opening edge is reported for use with an S–N curve.' },
    { key: 'gM', label: 'Material factor γM', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: `NORSOK N-004 / DNV-OS-C101 use 1.15; EN 1993 uses 1.0. ${CONFIRM}` },
    { key: 'nr', label: 'Elements from the opening to the plate edge', type: 'int', default: 16, min: 4, max: 64, group: 'Solver', role: 'control', help: 'Radial divisions, graded toward the opening; the same number is used along each 45° sector of the opening edge.' },
  ],
  mesh: { key: 'nr', label: 'Radial elements', levels: [8, 16, 32], qoi: 'sigMax' },
  compute(inp, ctx) {
    // SI: mm → m, MPa → Pa.
    const st = steelOf(ctx, inp.mat);
    const a = inp.d / 2e3, w = inp.W / 2e3, h = inp.H / 2e3, t = inp.t / 1e3, sx = inp.sx * 1e6, sy = inp.sy * 1e6;
    need(a < 0.9 * Math.min(w, h), 'The opening must be smaller than 90 % of the plate width and length.');
    need(Math.abs(sx) + Math.abs(sy) > 0, 'Enter a non-zero membrane stress on at least one pair of edges.');
    const mesh = run(() => plateWithHoleMesh({ a, w: h, h: w, nr: inp.nr, ntheta: inp.nr, grading: Math.max(2, Math.min(40, (Math.min(w, h) - a) / a * 2)) }), 'Mesh generation');
    const res = run(() => solvePlane({ mesh, E: st.E, nu: st.nu, thickness: t, mode: 'stress',
      fix: [{ nodes: mesh.nodeSets.symX, uy: 0 }, { nodes: mesh.nodeSets.symY, ux: 0 }],
      edgeLoads: [{ edges: mesh.edges.right, tx: sx }, { edges: mesh.edges.top, ty: sy }] }), 'Plane-stress solution');
    const S = res.stress, hole = mesh.nodeSets.hole;
    let sigMax = 0, vmMax = 0, tresca = 0, at = 0;
    const th = [], hoop = [];
    for (const n of hole) {
      const x = mesh.nodes[2 * n], y = mesh.nodes[2 * n + 1], c = x / a, s = y / a;
      // Hoop stress on the free edge: σθθ = σxx sin² + σyy cos² − 2 σxy sin cos.
      const hp = S.sxx[n] * s * s + S.syy[n] * c * c - 2 * S.sxy[n] * s * c;
      th.push(Math.atan2(y, x) * 180 / Math.PI); hoop.push(hp / 1e6);
      if (Math.abs(hp) > Math.abs(sigMax)) { sigMax = hp; at = th[th.length - 1]; }
    }
    for (let n = 0; n < mesh.nn; n++) { vmMax = Math.max(vmMax, S.vonMises[n]); tresca = Math.max(tresca, Math.max(Math.abs(S.s1[n] - S.s2[n]), Math.abs(S.s1[n]), Math.abs(S.s2[n]))); }
    const sRef = Math.max(Math.abs(sx), Math.abs(sy)), Kt = Math.abs(sigMax) / sRef;
    const kirsch = Math.max(Math.abs(3 * sx - sy), Math.abs(3 * sy - sx));
    const zz = run(() => zzError({ mesh, result: res, E: st.E, nu: st.nu, thickness: t }), 'Error indicator');
    const fyd = st.fy / inp.gM, uc = vmMax / fyd;
    const netUc = Math.max(Math.abs(sx) / (1 - a / w), Math.abs(sy) / (1 - a / h)) / fyd;
    const hot = Kt * inp.dSx;
    const order = th.map((v, i) => i).sort((p, q) => th[p] - th[q]);
    // Stress along the ligament (y axis at x = 0 when σx governs): nodes on the symY edge.
    const lig = Array.from(mesh.nodeSets.symY).sort((p, q) => mesh.nodes[2 * p + 1] - mesh.nodes[2 * q + 1]);
    const tReq = inp.t * uc, ok = uc <= 1 && netUc <= 1;
    const plateT = 4 * w * h * t * st.rho / 1e3;
    return {
      kpis: [
        { key: 'sigMax', label: 'Peak edge stress at the opening', unit: 'MPa', value: Math.abs(sigMax) / 1e6, explain: 'Largest tangential (hoop) stress on the free edge of the opening. This is the hot-spot stress for yield and for fatigue of a ground, unwelded edge.' },
        { key: 'Kt', label: 'Stress concentration factor', unit: '-', value: Kt, explain: 'Peak edge stress over the larger gross nominal stress. An infinite plate in uniaxial tension gives 3.0; a finite width raises it, biaxial tension lowers it.' },
        { key: 'vmMax', label: 'Peak von Mises stress', unit: 'MPa', value: vmMax / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Largest equivalent stress in the model. Above the design yield strength the edge yields locally; a small yielded zone redistributes, a large one needs a nonlinear check.' },
        { key: 'tresca', label: 'Peak Tresca stress', unit: 'MPa', value: tresca / 1e6, explain: 'Largest difference of principal stresses (maximum-shear criterion); slightly more conservative than von Mises.' },
        { key: 'hotSpot', label: 'Hot-spot stress range', unit: 'MPa', value: hot, explain: 'Nominal stress range times the computed concentration factor: enter it in the fatigue module with the S–N curve of the edge (C for a ground cut edge, lower for welded reinforcement).' },
        { key: 'zzEta', label: 'Discretisation error indicator η', unit: '%', value: zz.eta * 100, limit: 5, sense: 'max', explain: 'Zienkiewicz–Zhu relative error in the energy norm. Below about 5 % the mesh is adequate for stress; refine if it is higher.' },
        { key: 'kirsch', label: 'Infinite-plate (Kirsch) peak stress', unit: 'MPa', value: kirsch / 1e6, explain: 'Closed-form peak stress for an infinite plate, 3σ₁ − σ₂: the finite-element value should be a little higher because of the finite width.' },
      ],
      checks: [
        check('Local yield at the opening edge', vmMax / 1e6, fyd / 1e6, 'MPa', `von Mises, fy/γM with γM = ${inp.gM} (DNV-RP-C208 §5; NORSOK N-004 §6.1)`, 'Guards against yielding at the edge of the opening. Exceeding it means a plastic zone forms at the edge; under cyclic load this leads to low-cycle fatigue cracking.'),
        check('Net-section yield', netUc * fyd / 1e6, fyd / 1e6, 'MPa', 'Net-section stress σ·W/(W − d) ≤ fy/γM', 'Guards against yielding of the whole ligament beside the opening, which would be a real loss of member capacity, not a local peak.'),
      ],
      plots: [
        meshPlot({ title: 'von Mises stress (quarter model)', mesh, values: Array.from(S.vonMises, (v) => v / 1e6), valueLabel: 'von Mises stress', valueUnit: 'MPa', explain: 'Quarter of the plate (symmetry about both axes, opening at the origin). The peak sits on the opening edge, 90° from the larger applied stress, and decays within about one diameter.' }),
        { type: 'line', title: 'Hoop stress around the opening edge', x: { label: 'Angle from the x axis', unit: 'deg' }, y: { label: 'Edge stress', unit: 'MPa' },
          series: [{ name: 'Finite element', x: order.map((i) => round(th[i])), y: order.map((i) => round(hoop[i])) }, { name: 'Kirsch (infinite plate)', x: linspace(0, 90, 31).map((v) => round(v)), y: linspace(0, 90, 31).map((v) => round(((sx + sy) - 2 * (sx - sy) * Math.cos(2 * v * Math.PI / 180)) / 1e6)), dash: true }],
          explain: 'Tangential stress on the free edge. The maximum is at 90° for tension along x; a compressive value at 0° is expected for uniaxial load. The gap to the dashed closed form is the finite-width effect.' },
        { type: 'line', title: 'Stress across the ligament', x: { label: 'Distance from the opening centre (y)', unit: 'mm' }, y: { label: 'σxx', unit: 'MPa' },
          series: [{ name: 'σxx on x = 0', x: lig.map((n) => round(mesh.nodes[2 * n + 1] * 1e3)), y: lig.map((n) => round(S.sxx[n] / 1e6)) }, { name: 'Nominal', x: [a * 1e3, w * 1e3], y: [inp.sx, inp.sx], dash: true }],
          explain: 'Stress normal to the ligament from the opening edge to the plate edge: the peak decays quickly; the area under the curve equals the applied force (equilibrium).' },
      ],
      tables: [{ title: 'Mesh and equilibrium', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [{ q: 'Nodes', v: mesh.nn, u: '-' }, { q: 'Elements', v: mesh.ne, u: '-' }, { q: 'Applied force x (quarter)', v: round(res.appliedSum.fx / 1e3), u: 'kN' }, { q: 'Reaction x (quarter)', v: round(res.reactionSum.fx / 1e3), u: 'kN' }, { q: 'Strain energy (quarter)', v: round(res.strainEnergy), u: 'J' }, { q: 'Angle of the peak edge stress', v: round(at), u: 'deg' }],
        explain: 'Model size and the force balance of the quarter model: applied and reacting forces must cancel. Use the convergence study to confirm the peak stress is mesh-independent.' }],
      notes: [`Peak edge stress ${fmt(Math.abs(sigMax) / 1e6)} MPa at ${fmt(at)}° from the x axis: SCF ${fmt(Kt)} against ${fmt(kirsch / sRef)} for an infinite plate.`, `Mesh ${mesh.ne} Q4 elements, error indicator ${fmt(zz.eta * 100)} %.`],
      advice: [
        ok ? adv('engineering', uc > 0.9 ? 'warn' : 'ok', `Peak von Mises stress ${fmt(vmMax / 1e6)} MPa is ${fmt(uc * 100)} % of the design yield strength ${fmt(fyd / 1e6)} MPa; SCF ${fmt(Kt)}.`,
          `Keep the plate at ${inp.t} mm. Check fatigue with the hot-spot range ${fmt(hot)} MPa, and grind the cut edge if the opening is flame-cut.`, 'DNV-RP-C208 (2022) §5; DNV-RP-C203 (2021) §4.3',
          `The plate could go down to about ${fmt(tReq, 2)} mm for static strength (${steelCo2Text(plateT * (1 - uc))} saved per panel) if buckling and fatigue allow.`)
          : adv('engineering', 'fail', `Peak von Mises stress ${fmt(vmMax / 1e6)} MPa exceeds the design yield strength ${fmt(fyd / 1e6)} MPa (utilisation ${fmt(Math.max(uc, netUc))}).`,
            `Thicken the plate locally to at least ${fmt(Math.ceil(tReq), 3)} mm (insert plate of about two opening diameters), add a ring or coaming around the opening, or reduce the diameter; a nonlinear check can accept a small yielded zone.`, 'DNV-RP-C208 (2022) §5; NORSOK N-004 (2022) §6',
            `An insert plate of 2d × 2d at ${fmt(Math.ceil(tReq), 3)} mm adds about ${steelCo2Text((4 * a) ** 2 * Math.max(tReq - inp.t, 1) / 1e3 * st.rho / 1e3)} Local reinforcement uses far less steel than thickening the whole panel.`),
        adv('economic', 'info', `An insert or ring reinforcement costs fabrication hours (fit-up, full-penetration welds, NDT) rather than tonnes: the panel weighs ${fmt(plateT)} t.`,
          'Group penetrations so that one reinforced insert serves several, and fix their positions before plate cutting — late openings in a finished deck cost several times the original.', 'Design-for-fabrication practice; AACE RP 18R-97 (estimate maturity vs late change)'),
      ],
    };
  },
  verify: [
    { // Kirsch: infinite plate, uniaxial tension → σmax = 3σ. With d/W = 0.02 the finite-width correction is < 0.1 %.
      name: 'Small hole in a wide plate, uniaxial tension: σmax = 3σ (Kirsch)', source: 'Kirsch (1898); Timoshenko & Goodier, Theory of Elasticity, §35',
      inputs: { W: 5000, H: 5000, d: 100, sx: 100, sy: 0, nr: 32 }, expect: [{ kpi: 'sigMax', value: 300, tol: 0.02 }, { kpi: 'Kt', value: 3, tol: 0.02 }] },
    { // Equal biaxial tension: σmax = 3σ − σ = 2σ.
      name: 'Small hole, equal biaxial tension: σmax = 2σ', source: 'Timoshenko & Goodier, Theory of Elasticity, §35 (superposition of the Kirsch solution)',
      inputs: { W: 5000, H: 5000, d: 100, sx: 100, sy: 100, nr: 32 }, expect: [{ kpi: 'sigMax', value: 200, tol: 0.02 }] },
    { // Finite width, d/W = 0.3: Heywood's net-section factor Ktn = 2 + (1 − d/W)³ = 2.343, so σmax = Ktn σ/(1 − d/W) = 3.347σ.
      name: 'Finite-width strip, d/W = 0.3: σmax = [2 + (1 − d/W)³] σ/(1 − d/W) (Heywood)', source: 'Heywood (1952); Pilkey, Peterson’s Stress Concentration Factors, Chart 4.1',
      inputs: { W: 1000, H: 6000, d: 300, sx: 100, sy: 0, nr: 32 }, expect: [{ kpi: 'sigMax', value: 334.7, tol: 0.03 }] },
  ],
};
