// Tubular joint strength: every brace-to-chord connection of the model is classified by load path (K / X / Y) and
// checked against the API RP 2A simple-joint capacity with brace axial load, in-plane and out-of-plane bending.
import { need, run, fmt, round, adv, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, solve, jointData, jointLoads, plot3d, steelMass, foundNote, testFrame } from './_model.js';
import { apiJointCheck, classifyJoint } from '../../engine/structural/joints.js';

// lowest chord load factor of the classes and of bending that enter the joint capacity
const qf = (c) => Math.min(...Object.values(c.Qf).filter(Number.isFinite), 1);

export default {
  id: 's01-joints', suite: 1, title: 'Tubular joint strength and punching shear',
  summary: 'Capacity and utilisation of every brace-to-chord joint under brace axial load and bending, with K/X/Y classification, chord load factor, equivalent punching shear stress and the 50 % brace-strength rule.',
  covers: [7], fidelity: 'analytical',
  standards: ['API RP 2A-WSD 22nd ed. §7.2–7.3 (21st ed. Supp. 2 §4.2–4.3), Eq. 4.3-1 to 4.3-5', 'ISO 19902:2020 §14.2–14.3', 'NORSOK N-004 (2022) §6.4'],
  models: ['API RP 2A simple tubular joint strength (Qu, Qf)', 'Load-path joint classification (K, X, Y)', 'Punching shear on the chord wall', 'Brace load interaction P/Pa + (M/Ma)²ipb + |M/Ma|opb', 'Direct stiffness method with rigid joints'],
  limits: `Simple (unstiffened, non-overlapping, ungrouted) joints with member end forces from a rigid-jointed frame; equal-diameter crossings of X braces, overlapping, ring-stiffened and grouted joints are not checked, and joint flexibility is ignored. Capacities are allowable (working-stress) values with the safety factor given, so loads are entered unfactored. ${CODE_CAVEAT} ${MORISON_CAVEAT}`,
  searchTerms: 'tubular joint punching shear strength API RP 2A K joint',
  inputs: [
    ...stdInputs({ crit: { gG: 1.0, gE: 1.0 } }),
    { key: 'canT', label: 'Leg joint-can wall thickness', unit: 'mm', type: 'number', default: 55, min: 0, max: 200, group: 'Joints', role: 'geometry', design: true, help: 'Thickened chord wall at the leg joints; 0 or less than the leg wall uses the leg wall. Can steel is not added to the frame stiffness.' },
    { key: 'gap', label: 'Gap between braces at K joints', unit: 'mm', type: 'number', default: 75, min: -500, max: 1000, group: 'Joints', role: 'geometry', help: 'Toe-to-toe gap along the chord; 50 mm minimum is good practice, negative values are overlaps (outside the simple-joint check).' },
    { key: 'FS', label: 'Joint safety factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 2.5, group: 'Joints', role: 'criteria', help: `API RP 2A-WSD: 1.6 for operating conditions; with the one-third increase for extreme storm 1.6/1.33 = 1.2. ${CONFIRM}` },
    { key: 'classify', label: 'Joint classification', type: 'select', default: 'auto', group: 'Joints', role: 'criteria', options: [{ value: 'auto', label: 'By load path (K / X / Y fractions)' }, { value: 'Y', label: 'All as T/Y (conservative for balanced K)' }, { value: 'X', label: 'All as X (cross) joints' }] },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxUc' },
  compute(inp, ctx) {
    // SI conversion: mm → m here and inside the shared model; forces in N.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Storm' }), res = solve(J, [st.lc]).Storm;
    const jd = jointData(J, { canT: inp.canT });
    need(jd.length, 'The model has no brace-to-chord joint (a smaller tubular ending on a larger one).');
    const jl = jointLoads(J, res, jd), fy = J.steel.fy, fu = J.steel.fu;
    // classification per planar joint
    const cls = jd.map(() => ({ Y: 1 }));
    if (inp.classify === 'auto') { const planes = {}; jd.forEach((j, k) => (planes[j.plane] || (planes[j.plane] = [])).push(k)); for (const ks of Object.values(planes)) { const c = classifyJoint(ks.map((k) => ({ P: jl[k].P, thetaDeg: Math.min(jd[k].thetaDeg, 90), side: jd[k].side }))); ks.forEach((k, i) => { cls[k] = c[i]; }); } }
    else if (inp.classify === 'X') cls.forEach((c, k) => { cls[k] = { X: 1 }; });
    const rows = jd.map((j, k) => {
      const L = jl[k], c = run(() => apiJointCheck({ D: j.D, T: j.T, d: j.d, t: j.t, thetaDeg: Math.min(j.thetaDeg, 90), fyChord: fy, fuChord: fu, fyBrace: fy, classification: cls[k], gap: inp.gap / 1e3, P: L.P, Mipb: Math.abs(L.Mipb), Mopb: Math.abs(L.Mopb), chord: L.chord, FS: inp.FS }), `Joint check at node ${j.node}`);
      const sin = Math.sin(j.thetaDeg * Math.PI / 180), vp = Math.abs(L.P) * sin / (Math.PI * j.d * j.T), vpa = c.Pa * sin / (Math.PI * j.d * j.T);
      // A chord stressed beyond its own capacity drives the chord load factor to zero or below: the joint then has no
      // capacity left, whatever the sign the formula returns.
      const dead = !(c.Pa > 0) || !(c.MaIpb > 0) || !(c.MaOpb > 0) || !(qf(c) > 0);
      if (dead) { c.Pa = Math.max(c.Pa, 0); c.ucAxial = 99; c.ucIpb = Math.max(c.ucIpb, 0); c.ucOpb = Math.max(c.ucOpb, 0); if (!c.warnings.some((w) => /Qf/.test(w))) c.warnings = [...c.warnings, 'The chord is overstressed at this joint: the chord load factor Qf is not positive.']; }
      return { k, j, L, c, uc: !dead && Number.isFinite(c.uc) ? Math.min(Math.max(c.uc, 0), 99) : 99, vp, vpa: dead ? 0 : vpa, cl: cls[k], beta: j.d / j.D, gamma: j.D / (2 * j.T), tau: j.t / j.T, minUc: !dead && c.minimumCapacity.provided > 0 ? c.minimumCapacity.required / c.minimumCapacity.provided : 99, warn: c.warnings };
    });
    const gov = rows.reduce((a, r) => (r.uc > a.uc ? r : a)), mc = rows.reduce((a, r) => (r.minUc > a.minUc ? r : a)), nFail = rows.filter((r) => r.uc > 1).length, nWarn = rows.filter((r) => r.warn.length).length;
    const nodeUc = new Array(J.model.nodes.length).fill(0); rows.forEach((r) => { nodeUc[r.j.node] = Math.max(nodeUc[r.j.node], r.uc); });
    const top = rows.slice().sort((a, b) => b.uc - a.uc), clsTxt = (c) => ['K', 'X', 'Y'].filter((q) => c[q] > 0.005).map((q) => `${q} ${fmt(c[q] * 100, 2)}%`).join(' '), steel = steelMass(J, inp);
    // joint-can steel: a can of length 2.5 D at every leg node with a brace
    const canNodes = new Set(rows.filter((r) => r.j.c.group === 'leg').map((r) => r.j.node)), legD = Math.max(...rows.map((r) => r.j.D)), legT = Math.max(...J.orig.filter((m) => m.group === 'leg').map((m) => m.sec.t), 0);
    const canPerMm = canNodes.size * Math.PI * legD * 2.5 * legD * 1e-3 * J.steel.rho / 1e3, tReq = gov.j.T * Math.sqrt(Math.max(gov.uc, 0.01) / 0.9), ok = gov.uc <= 1;
    const checks = [
      check(`Joint strength interaction — node ${gov.j.node}, brace ${gov.j.brace}`, gov.uc, 1, 'interaction', `API RP 2A-WSD Eq. 4.3-5 (${clsTxt(gov.cl)}), FS = ${inp.FS}`, 'Brace axial load and bending against the allowable joint capacity of the chord wall. If exceeded the chord wall yields and punches or pulls out locally before the members reach their strength.'),
      check('Equivalent punching shear on the chord wall (same joint)', gov.vp / 1e6, gov.vpa / 1e6, 'MPa', 'vp = P·sinθ/(π·d·T) against the same joint capacity; API RP 2A-WSD (punching shear format of earlier editions)', 'Average shear stress through the chord wall around the brace footprint from the axial load alone. It is the classical way of looking at the same failure: the brace punching through, or tearing out of, the chord.'),
      check(`Joint develops 50 % of the brace strength — node ${mc.j.node}, brace ${mc.j.brace}`, mc.c.minimumCapacity.required / 1e3, mc.c.minimumCapacity.provided / 1e3, 'kN', 'API RP 2A-WSD §4.2.3 (minimum capacity)', 'Joints should not be the weak link: each should carry at least half the effective strength of the brace it connects, whatever the computed load.'),
    ];
    return {
      kpis: [
        { key: 'maxUc', label: 'Highest joint utilisation', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'Largest interaction ratio of any brace-to-chord joint under the applied load case. Above 1.0 the chord can is too thin or the brace too large for it.' },
        { key: 'ucAxial', label: 'Axial part of the governing utilisation', unit: '-', value: Math.min(gov.c.ucAxial, 99), explain: 'P/Pa of the governing joint. In braced jackets most of the utilisation should be axial; a large bending share indicates frame action or wave load on the brace.' },
        { key: 'PaGov', label: 'Allowable axial capacity of the governing joint', unit: 'kN', value: gov.c.Pa / 1e3, explain: 'Allowable brace axial load of the governing joint for its classification, including the chord load factor and the safety factor.' },
        { key: 'punching', label: 'Equivalent punching shear stress (governing joint)', unit: 'MPa', value: gov.vp / 1e6, explain: 'Brace axial force normal to the chord divided by the footprint perimeter times chord wall thickness.' },
        { key: 'QfGov', label: 'Chord load factor Qf (governing joint)', unit: '-', value: qf(gov.c), explain: 'Reduction of joint capacity for stress already in the chord. Values well below 1 mean the leg itself is highly stressed at the joint.' },
        { key: 'nJoints', label: 'Joints checked', unit: '-', value: rows.length, explain: 'Number of brace ends framing into a larger chord in the model.' },
        { key: 'nFail', label: 'Joints above 1.0', unit: '-', value: nFail, limit: 0, sense: 'max', explain: 'Count of joints whose interaction ratio exceeds unity.' },
        { key: 'minCapUc', label: 'Worst 50 % brace-strength ratio', unit: '-', value: Math.min(mc.minUc, 99), limit: 1, sense: 'max', explain: 'Half the brace strength over the joint capacity for the weakest joint; above 1.0 the joint fails the minimum-capacity rule.' },
        { key: 'tCanReq', label: 'Chord wall needed at the governing joint (0.9 utilisation)', unit: 'mm', value: tReq * 1e3, explain: 'Joint capacity grows with T², so the wall for a 0.9 utilisation follows from the square root of the utilisation. Compare with the can thickness given.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: 'Joint utilisation', values: nodeUc, label: 'Joint utilisation', unit: '-', u: res.displacements, explain: 'Highest joint interaction ratio at each node. Hot nodes in one bay point to a load-path problem; a general high level points to a thin chord.' }),
        { type: 'bar', title: 'Most utilised joints', categories: top.slice(0, 12).map((r) => `n${r.j.node}/b${r.j.brace}`), series: [{ name: 'Axial', y: top.slice(0, 12).map((r) => round(Math.min(r.c.ucAxial, 99))) }, { name: 'In-plane bending', y: top.slice(0, 12).map((r) => round(Math.min(r.c.ucIpb, 99))) }, { name: 'Out-of-plane bending', y: top.slice(0, 12).map((r) => round(Math.min(r.c.ucOpb, 99))) }], y: { label: 'Utilisation component', unit: '-' },
          explain: 'Axial, in-plane and out-of-plane terms of the interaction equation for the twelve most utilised joints (node / brace number). The terms add as axial + in-plane² + out-of-plane.' },
        { type: 'scatter', title: 'Utilisation against diameter ratio β', x: { label: 'β = d/D', unit: '-' }, y: { label: 'Joint utilisation', unit: '-' }, series: [{ name: 'Joints', x: rows.map((r) => round(r.beta)), y: rows.map((r) => round(r.uc)) }, { name: 'Limit', x: [0.2, 1], y: [1, 1], dash: true }],
          explain: 'Each point is one joint. Small-β joints are weak in punching; points above the line fail. The formulas are valid for 0.2 ≤ β ≤ 1.0 and 10 ≤ γ ≤ 50.' },
      ],
      tables: [
        { title: 'Joint checks (most utilised)', columns: [{ key: 'n', label: 'Node' }, { key: 'z', label: 'Elevation', unit: 'm' }, { key: 'b', label: 'Brace' }, { key: 'c', label: 'Chord' }, { key: 'beta', label: 'β' }, { key: 'gamma', label: 'γ' }, { key: 'th', label: 'θ', unit: 'deg' }, { key: 'cl', label: 'Class' }, { key: 'P', label: 'P (tension +)', unit: 'kN' }, { key: 'Pa', label: 'Pa', unit: 'kN' }, { key: 'Mi', label: 'Mipb', unit: 'kN·m' }, { key: 'Mo', label: 'Mopb', unit: 'kN·m' }, { key: 'Qf', label: 'Qf' }, { key: 'vp', label: 'Punching shear', unit: 'MPa' }, { key: 'uc', label: 'Utilisation' }],
          rows: top.slice(0, 25).map((r) => ({ n: r.j.node, z: round(r.j.z, 3), b: r.j.brace, c: r.j.chord, beta: round(r.beta, 3), gamma: round(r.gamma, 3), th: round(r.j.thetaDeg, 3), cl: clsTxt(r.cl), P: round(r.L.P / 1e3), Pa: round(r.c.Pa / 1e3), Mi: round(Math.abs(r.L.Mipb) / 1e3), Mo: round(Math.abs(r.L.Mopb) / 1e3), Qf: round(qf(r.c), 3), vp: round(r.vp / 1e6, 3), uc: round(r.uc) })),
          explain: 'Geometry, classification, brace loads, allowable axial capacity and utilisation of the twenty-five most utilised joints. Chord numbers refer to the member tables of the in-place analysis.' },
        { title: 'Validity warnings', columns: [{ key: 'n', label: 'Node' }, { key: 'b', label: 'Brace' }, { key: 'w', label: 'Warning' }], rows: rows.filter((r) => r.warn.length).slice(0, 20).map((r) => ({ n: r.j.node, b: r.j.brace, w: r.warn.join('; ') })), explain: 'Joints outside the validity range of the simple-joint formulas (β, γ, θ, yield strength or gap). Their utilisation is indicative only and needs a detailed assessment.' },
      ],
      notes: [`Governing joint: node ${gov.j.node} (z = ${fmt(gov.j.z)} m), brace ${gov.j.brace} on chord ${gov.j.chord}, β = ${fmt(gov.beta)}, γ = ${fmt(gov.gamma)}, ${clsTxt(gov.cl)}; utilisation ${fmt(gov.uc)}.`,
        `${rows.length} joints checked, ${nFail} above 1.0, ${nWarn} outside the validity range. Chord wall at leg joints ${fmt(Math.max(inp.canT, legT * 1e3))} mm.`,
        `Unfactored loads with FS = ${fmt(inp.FS)}; wave by ${st.wf.label}, heading ${fmt(st.beta)}°.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `All ${rows.length} joints pass: highest utilisation ${fmt(gov.uc)} at node ${gov.j.node}; the 50 % brace-strength rule gives ${fmt(mc.minUc)}.`, gov.uc < 0.6 ? `The can could be thinned to about ${fmt(Math.max(tReq * 1e3, legT * 1e3))} mm at 0.9 utilisation, but check fatigue first: hot-spot stress, not static strength, usually sets can thickness.` : 'Keep the can thickness; confirm classification for the other headings and carry the same joints into the fatigue module.', 'API RP 2A-WSD 22nd ed. §7.3; ISO 19902:2020 §14.3', `Each millimetre of can wall on the ${canNodes.size} leg joints is about ${steelCo2(canPerMm)}`)
          : adv('engineering', 'fail', `${nFail} joint(s) fail: utilisation ${fmt(gov.uc)} at node ${gov.j.node}, brace ${gov.j.brace} (${clsTxt(gov.cl)}).`, `Increase the chord can to about ${fmt(tReq * 1e3)} mm at that joint (capacity ∝ T²), raise β with a larger brace stub, or grout / ring-stiffen an existing joint.`, 'API RP 2A-WSD 22nd ed. §7.3; ISO 19902:2020 §14.3', `Thickening the cans by ${fmt(Math.max(tReq * 1e3 - gov.j.T * 1e3, 1))} mm on ${canNodes.size} leg joints is about ${steelCo2(canPerMm * Math.max(tReq * 1e3 - gov.j.T * 1e3, 1))}`),
        adv('economic', 'info', `Joint cans are ${fmt(canPerMm * Math.max(inp.canT - legT * 1e3, 0))} t of thick, higher-quality steel over the leg wall on ${canNodes.size} nodes.`, 'Thick cans need through-thickness (Z-quality) plate, preheat and more NDT: standardise on one or two can thicknesses and keep β above about 0.4 so cans stay moderate. A simple joint is far cheaper than a ring-stiffened or grouted one, offshore or in the yard.', `${BASIS_COST}; fabrication practice per ISO 19902:2020 §19–§20`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Y joint, θ = 45°, D = 1.0 m, T = 25 mm (γ = 20), d = 0.5 m (β = 0.5), brace in compression, chord nearly unloaded (Qf ≈ 1):
    // Qu = 2.8 + (20 + 0.8γ)β^1.6 = 2.8 + 36·0.32988 = 14.676; Pa = Qu·Fy·T²/(FS·sinθ) = 14.676·355e6·0.025²/(1.2·0.70711) = 3837.5 kN.
    // In-plane: Qu = (5 + 0.7γ)β^1.2 = 19·0.43528 = 8.270, Ma = 8.270·355e6·0.025²·0.5/(1.2·0.70711) = 1081.3 kN·m.
    // Load 10 kN at the brace tip (20 m out): P = 7.071 kN, M = 200 kN·m → 7.071/3837.5 + (200/1081.3)² = 0.0361.
    { name: 'Y joint capacity Pa = Qu·Qf·Fy·T²/(FS·sinθ) and interaction', source: 'API RP 2A-WSD 21st ed. Supplement 2 / 22nd ed., Table 4.3-1 and Eq. 4.3-1, 4.3-5 (hand calculation)',
      inputs: { geom: testFrame('tjoint', { thetaDeg: 45 }), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 0.01, found: 'fixed', canT: 0, FS: 1.2, nSub: 1 },
      expect: [{ kpi: 'PaGov', value: 3837.5, tol: 0.01 }, { kpi: 'maxUc', value: 0.0361, tol: 0.03 }] },
  ],
};
