// Static strength and stress analysis: one factored static load set (gravity, topsides, design wave and tabulated
// operational loads) solved by the direct stiffness method, with section forces and stresses in every member.
import { fmt, round, adv, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, solve, memberRows, byGroup, pileForces, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows, groupChecks, tableLoads, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { memberDiagram, sectionStress } from '../../engine/structural/frame3d.js';

export default {
  id: 's01-static-strength', suite: 1, title: 'Static strength and stress analysis',
  summary: 'Member axial forces, shears, moments and torsion with axial, bending, shear, von Mises and Tresca stresses for one factored static load set, checked against yield and the tubular member code.',
  covers: [2], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §12, §13.2–13.3', 'API RP 2A-WSD 22nd ed. §6.2–6.3', 'EN 1993-1-1:2005 §6.2.1 (von Mises yield criterion)'],
  models: ['Direct stiffness method', 'Euler–Bernoulli / Timoshenko beam theory', 'Saint-Venant torsion', 'Generalised Hooke\'s law (linear elasticity)', 'von Mises and Tresca yield criteria', 'Principle of virtual work (Clapeyron energy check)', 'Morison equation', 'ISO 19902 tubular member checks'],
  limits: `Linear elastic, small displacements, rigid joints and nominal (beam) stresses: no joint flexibility, hot-spot stress or plasticity — see the joint, global–local and nonlinear modules. ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'jacket static analysis member stress von Mises tubular',
  inputs: [
    ...stdInputs(),
    { key: 'loads', label: 'Additional static loads', type: 'table', group: 'Topsides and gravity loads', role: 'bc',
      columns: [{ key: 'name', label: 'Load', type: 'text' }, { key: 'z', label: 'Elevation', unit: 'm', type: 'number' }, { key: 'Fx', label: 'Fx', unit: 'kN', type: 'number' }, { key: 'Fy', label: 'Fy', unit: 'kN', type: 'number' }, { key: 'Fz', label: 'Fz (up +)', unit: 'kN', type: 'number' }],
      default: [{ name: 'Crane slewing reaction', z: 15, Fx: 300, Fy: 200, Fz: -1200 }, { name: 'Conductor and riser guide loads', z: -20, Fx: 250, Fy: 0, Fz: -400 }, { name: 'Boat landing and fender weight', z: 2, Fx: 0, Fy: 0, Fz: -600 }],
      help: 'Operational and equipment loads at their elevation (still water = 0), shared by the leg nodes nearest to that level. They take the gravity-load factor.' },
    { key: 'gM', label: 'Material factor on yield for the stress check', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Design criteria', role: 'criteria', help: 'von Mises stress is checked against fy/γM (1.15 in DNV-OS-C101 / NORSOK N-004; 1/0.6 for a basic allowable-stress check).' },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI conversion: kN, MN and mm inputs are converted inside the shared model; stresses reported in MPa.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Static' });
    st.lc.nodal = [...st.lc.nodal, ...tableLoads(J, inp.loads, inp.gG)];
    const res = solve(J, [st.lc]).Static, rows = memberRows(J, res, inp, { name: 'Static', hydro: true }), groups = byGroup(rows), dp = displacements(J, res), pf = pileForces(J, res);
    const fyd = J.steel.fy / inp.gM, gov = rows.reduce((a, r) => (r.vm > a.vm ? r : a)), govUc = rows.reduce((a, r) => (r.uc > a.uc ? r : a));
    const tresca = (r) => Math.sqrt((Math.abs(r.sigA) + r.sigB) ** 2 + 4 * r.tau ** 2);
    // Clapeyron: for a linear system the strain energy equals half the work of the applied loads. Reactions do no
    // work on rigid supports, so ½ Σ R·u over the supports is the energy stored in the foundation springs.
    const u = res.displacements, eSpring = res.reactions.reduce((s, r) => s + 0.5 * r.R.reduce((q, v, k) => q + v * -u[6 * r.node + k], 0), 0);
    // section-force and stress diagrams along the member with the highest von Mises stress
    const gm = J.orig[gov.i], dia = { s: [], N: [], M: [], V: [], vm: [] };
    let s0 = 0;
    for (const e of J.memberMap[gov.i]) {
      const d = memberDiagram(res, e, 9);
      d.s.forEach((s, k) => { const M = Math.hypot(d.My[k], d.Mz[k]), V = Math.hypot(d.Vy[k], d.Vz[k]); dia.s.push(round(s0 + s)); dia.N.push(round(d.N[k] / 1e3)); dia.M.push(round(M / 1e3)); dia.V.push(round(V / 1e3)); dia.vm.push(round(sectionStress({ section: gm.sec, N: d.N[k], Vy: V, T: d.T[k], My: d.My[k], Mz: d.Mz[k] }).vonMises / 1e6)); });
      s0 += d.s[d.s.length - 1];
    }
    const steel = steelMass(J, inp), gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a)), gSteel = groupSteel(J, gg.group), ok = gov.vm <= fyd && govUc.uc <= 1;
    const top = J.topNodes.concat(J.mudNodes).slice(0, 12);
    const checks = [
      check('von Mises stress (most stressed member)', gov.vm / 1e6, fyd / 1e6, 'MPa', `fy/γM with γM = ${inp.gM}; EN 1993-1-1 Eq. 6.1`, 'Equivalent stress from axial, bending, shear and torsion at the worst section. Above the limit the extreme fibre yields under the factored load.'),
      check('Tresca stress (same member)', tresca(gov) / 1e6, fyd / 1e6, 'MPa', 'Maximum shear stress criterion, σ1 − σ3 ≤ fy/γM', 'The Tresca criterion is up to 15 % more conservative than von Mises when shear is present; it is the yield criterion behind several pressure-vessel codes.'),
      ...groupChecks(groups, inp),
    ];
    return {
      kpis: [
        { key: 'maxVm', label: 'Highest von Mises stress', unit: 'MPa', value: gov.vm / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Largest nominal equivalent stress in any member under the factored load set. It must stay below yield divided by the material factor.' },
        { key: 'maxAxial', label: 'Highest axial stress (magnitude)', unit: 'MPa', value: Math.max(...rows.map((r) => Math.abs(r.sigA))) / 1e6, explain: 'Largest N/A in any member. Braces and legs of a jacket work mainly axially, so this is normally most of the stress.' },
        { key: 'maxBending', label: 'Highest bending stress', unit: 'MPa', value: Math.max(...rows.map((r) => r.sigB)) / 1e6, explain: 'Largest M/W in any member. High bending relative to axial stress points to frame action, local wave load on a long span or a missing brace.' },
        { key: 'maxShear', label: 'Highest shear stress', unit: 'MPa', value: Math.max(...rows.map((r) => r.tau)) / 1e6, explain: 'Peak shear stress from transverse shear and torsion. It is normally small in tubular frames; a high value flags a short, stiff member.' },
        { key: 'maxUc', label: 'Highest member utilisation', unit: '-', value: govUc.uc, limit: 1, sense: 'max', explain: 'Largest code utilisation including column buckling and beam-column interaction, which the stress check alone does not cover.' },
        { key: 'maxDisp', label: 'Largest nodal displacement', unit: 'mm', value: dp.max * 1e3, explain: 'Largest translation of any node. Useful as a sanity check of stiffness and as the convergence quantity for the mesh.' },
        { key: 'deckDisp', label: 'Lateral displacement at the jacket top', unit: 'mm', value: dp.top * 1e3, explain: 'Horizontal movement of the leg tops under the load set, including foundation flexibility.' },
        { key: 'maxRot', label: 'Largest nodal rotation', unit: 'mrad', value: dp.rotMax * 1e3, explain: 'Largest rotation of any node; it indicates member end rotations that joints and attachments have to follow.' },
        { key: 'sumV', label: 'Total vertical reaction', unit: 'MN', value: pf.sumV / 1e6, explain: 'Sum of the vertical support reactions; equals the factored gravity load less buoyancy plus the vertical wave load.' },
        { key: 'sumH', label: 'Total horizontal reaction', unit: 'MN', value: pf.sumH / 1e6, explain: 'Resultant horizontal support reaction: the factored base shear of this load set.' },
        { key: 'eqResidual', label: 'Equilibrium residual (relative)', unit: '-', value: res.equilibrium.relative < 1e-9 ? 0 : Number(res.equilibrium.relative.toPrecision(3)), explain: 'Out-of-balance force of the solution relative to the applied load. Values above about 1e-6 indicate an ill-conditioned or nearly unstable model. A residual below 1e-9 of the applied load is arithmetic round-off and is reported as 0.' },
        { key: 'energy', label: 'Strain energy', unit: 'kJ', value: res.strainEnergy / 1e3, explain: 'Elastic energy stored in the members. With the energy of the foundation springs it equals half the work of the applied loads (Clapeyron), which is an independent check of the solution.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: 'von Mises stress and deformed shape', values: nodeValuesFromMembers(J, rows.map((r) => r.vm / 1e6)), label: 'von Mises stress', unit: 'MPa', u: res.displacements, explain: 'Highest nominal von Mises stress of the members meeting at each node, on the deflected shape. Look for isolated hot members, which carry more than their neighbours.' }),
        { type: 'line', title: `Section forces along member ${gov.i} (${gm.label})`, x: { label: 'Distance from end A', unit: 'm' }, y: { label: 'Force / moment', unit: 'kN, kN·m' }, series: [{ name: 'Axial N (kN)', x: dia.s, y: dia.N }, { name: 'Shear V (kN)', x: dia.s, y: dia.V }, { name: 'Moment M (kN·m)', x: dia.s, y: dia.M }],
          explain: 'Axial force, resultant shear and resultant bending moment along the most stressed member. The moment peak shows whether the ends (frame action) or the span (local wave load) govern.' },
        { type: 'line', title: `von Mises stress along member ${gov.i}`, x: { label: 'Distance from end A', unit: 'm' }, y: { label: 'Stress', unit: 'MPa' }, series: [{ name: 'von Mises', x: dia.s, y: dia.vm }, { name: 'fy/γM', x: [dia.s[0], dia.s[dia.s.length - 1]], y: [round(fyd / 1e6), round(fyd / 1e6)], dash: true }],
          explain: 'Equivalent stress along the member against the design yield stress. A crossing of the dashed line marks where the section yields first.' },
        { type: 'bar', title: 'Stress components by member group', categories: groups.map((g) => g.label), series: [{ name: 'Axial', y: groups.map((g) => round(Math.max(...rows.filter((r) => r.group === g.group).map((r) => Math.abs(r.sigA))) / 1e6)) }, { name: 'Bending', y: groups.map((g) => round(Math.max(...rows.filter((r) => r.group === g.group).map((r) => r.sigB)) / 1e6)) }, { name: 'von Mises', y: groups.map((g) => round(Math.max(...rows.filter((r) => r.group === g.group).map((r) => r.vm)) / 1e6)) }], y: { label: 'Stress', unit: 'MPa' },
          explain: 'Largest axial, bending and equivalent stress in each member group. Groups far below the others are candidates for lighter sections.' },
      ],
      tables: [
        { title: 'Most utilised members', columns: memberCols, rows: memberTableRows(rows), explain: 'The fifteen members with the highest code utilisation, with their extreme section forces: the list for joint, local and fatigue checks.' },
        { title: 'Stress components of the most stressed members', columns: [{ key: 'i', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'sa', label: 'Axial', unit: 'MPa' }, { key: 'sb', label: 'Bending', unit: 'MPa' }, { key: 'tau', label: 'Shear', unit: 'MPa' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'tr', label: 'Tresca', unit: 'MPa' }, { key: 'r', label: 'vM / (fy/γM)' }],
          rows: rows.slice().sort((a, b) => b.vm - a.vm).slice(0, 12).map((r) => ({ i: r.i, g: r.label, sa: round(r.sigA / 1e6), sb: round(r.sigB / 1e6), tau: round(r.tau / 1e6), vm: round(r.vm / 1e6), tr: round(tresca(r) / 1e6), r: round(r.vm / fyd) })),
          explain: 'Axial (tension positive), bending and shear stress with the two equivalent stresses for the twelve most stressed members. The last column is the stress utilisation.' },
        { title: 'Displacements and rotations at leg tops and pile heads', columns: [{ key: 'n', label: 'Node' }, { key: 'z', label: 'Elevation', unit: 'm' }, { key: 'ux', label: 'ux', unit: 'mm' }, { key: 'uy', label: 'uy', unit: 'mm' }, { key: 'uz', label: 'uz', unit: 'mm' }, { key: 'rx', label: 'rx', unit: 'mrad' }, { key: 'ry', label: 'ry', unit: 'mrad' }, { key: 'rz', label: 'rz', unit: 'mrad' }],
          rows: top.map((n) => ({ n, z: round(J.model.nodes[n][2]), ux: round(u[6 * n] * 1e3), uy: round(u[6 * n + 1] * 1e3), uz: round(u[6 * n + 2] * 1e3), rx: round(u[6 * n + 3] * 1e3), ry: round(u[6 * n + 4] * 1e3), rz: round(u[6 * n + 5] * 1e3) })),
          explain: 'Global displacements and rotations of the topsides interface nodes and of the pile heads. Pile-head values are zero for rigid supports and show the foundation movement for springs.' },
        { title: 'Support reactions', columns: [{ key: 'node', label: 'Leg node' }, { key: 'axial', label: 'Vertical (compression +)', unit: 'MN' }, { key: 'shear', label: 'Shear', unit: 'MN' }, { key: 'moment', label: 'Moment', unit: 'MN·m' }], rows: pf.piles.map((p) => ({ node: p.node, axial: round(p.axial / 1e6), shear: round(p.shear / 1e6), moment: round(p.moment / 1e6) })), explain: 'Factored reactions at each leg. Negative vertical values are uplift.' },
      ],
      notes: [`Highest von Mises stress ${fmt(gov.vm / 1e6)} MPa in ${gm.label.toLowerCase()} member ${gov.i} (axial ${fmt(gov.sigA / 1e6)} MPa, bending ${fmt(gov.sigB / 1e6)} MPa) against ${fmt(fyd / 1e6)} MPa.`,
        `Highest code utilisation ${fmt(govUc.uc)} in member ${govUc.i} (${govUc.ref}). Wave by ${st.wf.label}, heading ${fmt(st.beta)}°, base shear ${fmt(st.baseShear / 1e6)} MN unfactored.`,
        `Energy check: member strain energy ${fmt(res.strainEnergy / 1e3)} kJ, foundation springs ${fmt(eSpring / 1e3)} kJ.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Stresses are within yield: von Mises ${fmt(gov.vm / 1e6)} MPa ≤ ${fmt(fyd / 1e6)} MPa and member utilisation ${fmt(govUc.uc)}.`, `Stability governs before stress in slender braces: the code utilisation (${fmt(govUc.uc)}) is ${govUc.uc > gov.vm / fyd ? 'higher' : 'lower'} than the stress ratio (${fmt(gov.vm / fyd)}). Take the member end forces into the joint check and confirm the load table with the topsides weight report.`, BASIS_MEMBER, `The ${gg.label.toLowerCase()} run at ${fmt(gg.uc)}; bringing them to 0.85 would change steel by about ${steelCo2(gSteel * (1 - gg.uc / 0.85))}`)
          : adv('engineering', 'fail', `Strength is exceeded: von Mises ${fmt(gov.vm / 1e6)} MPa against ${fmt(fyd / 1e6)} MPa, member utilisation ${fmt(govUc.uc)} (${govUc.ref}).`, `Increase the wall thickness of the ${gg.label.toLowerCase()} by about ${fmt((Math.max(gg.uc, gov.vm / fyd) / 0.9 - 1) * 100)} % or move the tabulated loads to a framed level; thicker wall adds strength without adding wave load.`, BASIS_MEMBER, `The extra wall is roughly ${steelCo2(gSteel * Math.max(Math.max(gg.uc, gov.vm / fyd) / 0.9 - 1, 0.05))}`),
        adv('economic', 'info', `${fmt(steel)} t of jacket steel carries this load set at a peak stress ratio of ${fmt(gov.vm / fyd)}.`, 'Group members into a few standard sizes sized to 0.8–0.9 utilisation: each distinct can size adds procurement and welding-procedure cost, while unused capacity is paid for in tonnes of fabricated steel and in lift weight.', BASIS_COST),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Cantilever tube, L = 95 m, D = 1.5 m, t = 45 mm, E = 210 GPa, tip load P = 10 kN:
    // I = π(1.5⁴ − 1.41⁴)/64 = 0.054485 m⁴, A = 0.20570 m²; bending δ = PL³/(3EI) = 1e4·95³/(3·2.1e11·0.054485) = 249.78 mm,
    // shear δ = PL/(G·A/2) = 1e4·95/(8.077e10·0.10285) = 0.11 mm, total 249.89 mm; σb = PL·(D/2)/I = 13.077 MPa.
    { name: 'Cantilever tube with tip load: δ = PL³/3EI, σ = PL/W', source: 'Euler–Bernoulli beam theory (Roark\'s Formulas for Stress and Strain, Table 8.1 case 1a)',
      inputs: { geom: testFrame('column'), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 0, FxTop: 10, found: 'fixed', gG: 1, gE: 1, loads: [], nSub: 2 },
      expect: [{ kpi: 'deckDisp', value: 249.89, tol: 0.002 }, { kpi: 'maxBending', value: 13.077, tol: 0.002 }, { kpi: 'sumH', value: 0.01, tol: 1e-6 }] },
  ],
};
