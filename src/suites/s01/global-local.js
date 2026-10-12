// Global–local analysis: member end forces from the global frame are taken to the critical tubular joints, where
// Efthymiou stress concentration factors give the hot-spot stresses around the brace–chord intersection, the
// local yielding (shakedown) margin and a Neuber estimate of the plastic strain concentration.
import { need, run, fmt, round, adv, check, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, solve, memberRows, jointData, jointLoads, plot3d, steelMass, foundNote, testFrame } from './_model.js';
import { tyJointHotSpots } from '../../engine/structural/scf.js';

export default {
  id: 's01-global-local', suite: 1, title: 'Global and local analysis — joint hot-spot stresses',
  summary: 'Transfers brace end forces of the global frame to the most loaded joints and evaluates hot-spot stresses at eight points around each intersection, the local yield margin and the boundary forces for a detailed sub-model.',
  covers: [4], fidelity: 'reduced-order',
  standards: ['DNV-RP-C203 (2021) §3.3 and Appendix B (Efthymiou SCFs)', 'ISO 19902:2020 §12.3, §16.10 (local analysis, hot-spot stress)', 'DNV-RP-C208 (2019) §4 (local strain criteria)', 'API RP 2A-WSD 22nd ed. §8.3 (SCFs)'],
  models: ['Direct stiffness method (global beam model)', 'Global–local transfer of member end forces (sub-model boundary loads)', 'Efthymiou stress concentration factors for T/Y joints', 'Hot-spot stress superposition at crown and saddle points', 'Neuber rule for local plastic strain', 'Elastic shakedown criterion'],
  limits: `Parametric SCFs replace the shell or solid sub-model: every brace is treated as a T/Y joint (conservative for balanced K action), joints are unstiffened and non-overlapping, and weld geometry, residual stress and misalignment are not modelled. Beam–shell–solid coupling and nonlinear local FE are outside the module — the boundary-force table is the hand-over to such a sub-model. ${CODE_CAVEAT} ${MORISON_CAVEAT}`,
  searchTerms: 'tubular joint hot spot stress Efthymiou stress concentration global local',
  inputs: [
    ...stdInputs(),
    { key: 'canT', label: 'Leg joint-can wall thickness', unit: 'mm', type: 'number', default: 55, min: 0, max: 200, group: 'Local model', role: 'geometry', design: true, help: 'Chord wall at leg joints; 0 uses the leg wall.' },
    { key: 'nLocal', label: 'Joints analysed locally', type: 'int', default: 10, min: 1, max: 60, group: 'Local model', role: 'control', help: 'The joints with the highest nominal brace stress are taken to the local analysis.' },
    { key: 'chordC', label: 'Chord end fixity parameter C', unit: '-', type: 'number', default: 0.7, min: 0.5, max: 1, group: 'Local model', role: 'bc', help: 'Boundary condition of the local model: 0.5 pinned chord ends, 1.0 fully fixed, 0.7 typical for jacket legs (DNV-RP-C203 Table B-1).' },
    { key: 'hsLimit', label: 'Hot-spot stress limit as a multiple of yield', unit: '-', type: 'number', default: 2, min: 1, max: 3, group: 'Design criteria', role: 'criteria', help: `Elastic shakedown: a linear-elastic hot-spot stress up to 2·fy shakes down to elastic cycling after first yield. ${CONFIRM}` },
    { key: 'epsLimit', label: 'Local plastic strain limit', unit: '%', type: 'number', default: 4, min: 0.5, max: 15, group: 'Design criteria', role: 'criteria', help: `Order of 4–5 % for welded structural steel in DNV-RP-C208 type local checks. ${CONFIRM}` },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxHotSpot' },
  compute(inp, ctx) {
    // SI conversion: mm → m; stresses reported in MPa.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Storm' }), res = solve(J, [st.lc]).Storm, rowsM = memberRows(J, res, inp, { name: 'Storm' });
    const jd = jointData(J, { canT: inp.canT });
    need(jd.length, 'The model has no brace-to-chord joint for a local analysis.');
    const jl = jointLoads(J, res, jd), fy = J.steel.fy, E = J.steel.E;
    const nom = jd.map((j, k) => ({ k, s: Math.abs(jl[k].sx) + Math.hypot(jl[k].sip, jl[k].sop) })).sort((a, b) => b.s - a.s).slice(0, inp.nLocal);
    const rows = nom.map(({ k, s }) => {
      const j = jd[k], L = jl[k], h = run(() => tyJointHotSpots({ D: j.D, T: j.T, d: j.d, t: j.t, L: j.Lc, thetaDeg: Math.min(j.thetaDeg, 90), sigmaX: L.sx, sigmaMy: L.sip, sigmaMz: L.sop, chordEnds: 'general', C: inp.chordC }), `Hot-spot stresses at node ${j.node}`);
      const hs = Math.max(h.chord.maxAbs, h.brace.maxAbs), side = h.chord.maxAbs >= h.brace.maxAbs ? 'chord' : 'brace', g = side === 'chord' ? h.chord : h.brace;
      // Neuber: σ·ε = σ_hs²/E with σ = fy on the plateau → total strain σ_hs²/(E·fy)
      const eps = hs > fy ? hs * hs / (E * fy) - fy / E : 0;
      return { j, L, h, hs, side, pt: g.governingPoint, nomS: s, scfEff: s > 0 ? hs / s : 0, eps };
    });
    const gov = rows.reduce((a, r) => (r.hs > a.hs ? r : a)), nodeHs = new Array(J.model.nodes.length).fill(0); rows.forEach((r) => { nodeHs[r.j.node] = Math.max(nodeHs[r.j.node], r.hs / 1e6); });
    const maxScf = Math.max(gov.h.scf.axial.chordSaddle, gov.h.scf.axial.chordCrown, gov.h.scf.axial.braceSaddle, gov.h.scf.axial.braceCrown);
    const ucMember = Math.max(...rowsM.map((r) => r.uc)), steel = steelMass(J, inp), lim = inp.hsLimit * fy, ok = gov.hs <= lim && gov.eps * 100 <= inp.epsLimit;
    const pts = [1, 2, 3, 4, 5, 6, 7, 8], warn = rows.filter((r) => r.h.warnings.length);
    // chord wall for the limit: chord-side SCFs scale roughly with τ·γ ∝ 1/T^1.5 at constant brace
    const tReq = gov.j.T * (gov.hs / lim) ** (1 / 1.5);
    return {
      kpis: [
        { key: 'maxHotSpot', label: 'Highest hot-spot stress', unit: 'MPa', value: gov.hs / 1e6, limit: lim / 1e6, sense: 'max', explain: 'Largest linear-elastic stress at the weld toe of any analysed joint under the factored storm load. It may exceed yield locally; up to twice yield the detail shakes down and cycles elastically afterwards.' },
        { key: 'nomStress', label: 'Nominal brace stress at that joint', unit: 'MPa', value: gov.nomS / 1e6, explain: 'Axial plus bending stress in the brace from the global beam model: what the frame analysis sees before the local concentration.' },
        { key: 'scfEff', label: 'Effective concentration (hot-spot / nominal)', unit: '-', value: gov.scfEff, explain: 'Ratio of local to global stress at the governing joint, combining the axial and bending SCFs with the actual load mix.' },
        { key: 'scfAxialMax', label: 'Largest axial SCF of the governing joint', unit: '-', value: maxScf, explain: 'Highest of the four Efthymiou axial factors (chord and brace side, crown and saddle). Values above about 6 flag a thin chord or a small diameter ratio.' },
        { key: 'scfChordSaddle', label: 'Chord saddle SCF, axial (governing joint)', unit: '-', value: gov.h.scf.axial.chordSaddle, explain: 'Efthymiou factor at the chord saddle for brace axial load: the usual peak location for T and Y joints.' },
        { key: 'plasticStrain', label: 'Local plastic strain (Neuber)', unit: '%', value: gov.eps * 100, limit: inp.epsLimit, sense: 'max', explain: 'Estimate of the plastic strain at the hot spot when the elastic stress exceeds yield, for an elastic–perfectly plastic material. Zero means the hot spot stays elastic.' },
        { key: 'ucMember', label: 'Highest global member utilisation', unit: '-', value: ucMember, limit: 1, sense: 'max', explain: 'Code utilisation from the global analysis of the same load case, for comparison with the local result.' },
        { key: 'nYield', label: 'Analysed joints with hot-spot stress above yield', unit: '-', value: rows.filter((r) => r.hs > fy).length, explain: 'Joints where first loading causes local yielding at the weld toe. This is acceptable for strength within the shakedown limit but these are the fatigue-critical locations.' },
      ],
      checks: [
        check(`Hot-spot stress — node ${gov.j.node}, brace ${gov.j.brace} (${gov.side} side, point ${gov.pt})`, gov.hs / 1e6, lim / 1e6, 'MPa', `Elastic shakedown, σhs ≤ ${inp.hsLimit}·fy; hot-spot stress per DNV-RP-C203 §3.3.2`, 'Linear-elastic peak stress at the weld toe against the shakedown limit. Beyond it the detail yields in both directions every storm cycle and fails by low-cycle fatigue or incremental collapse.'),
        check('Local plastic strain at the hot spot', gov.eps * 100, inp.epsLimit, '%', 'Neuber rule; strain limit per DNV-RP-C208 §4 (project criterion)', 'Plastic strain concentrated at the weld toe when the elastic hot-spot stress exceeds yield. Above the limit a crack can initiate in a single extreme event.'),
        check('Global member strength (same load case)', ucMember, 1, 'interaction', `${inp.code} tubular member check`, 'Member utilisation from the global beam model. Local analysis only makes sense for a frame that passes globally.'),
      ],
      plots: [
        { type: 'line', title: `Hot-spot stress around the intersection — node ${gov.j.node}, brace ${gov.j.brace}`, x: { label: 'Point (1, 5 crown; 3, 7 saddle)', unit: '-' }, y: { label: 'Hot-spot stress', unit: 'MPa' }, series: [{ name: 'Chord side', x: pts, y: Array.from(gov.h.chord.sigma, (v) => round(v / 1e6)) }, { name: 'Brace side', x: pts, y: Array.from(gov.h.brace.sigma, (v) => round(v / 1e6)) }, { name: 'Yield', x: [1, 8], y: [round(fy / 1e6), round(fy / 1e6)], dash: true }],
          explain: 'Stress at the eight standard points around the weld for the governing joint. Peaks at 3 and 7 mean saddle (axial or out-of-plane) control, peaks at 1 and 5 mean in-plane bending control; the dashed line is first yield.' },
        { type: 'bar', title: 'Nominal and hot-spot stress of the analysed joints', categories: rows.map((r) => `n${r.j.node}/b${r.j.brace}`), series: [{ name: 'Nominal', y: rows.map((r) => round(r.nomS / 1e6)) }, { name: 'Hot spot', y: rows.map((r) => round(r.hs / 1e6)) }], y: { label: 'Stress', unit: 'MPa' },
          explain: 'Global (nominal) brace stress and the local hot-spot stress for each analysed joint. The ratio of the two bars is the effective stress concentration.' },
        plot3d(J, { title: 'Localised high-stress regions', values: nodeHs, label: 'Hot-spot stress', unit: 'MPa', u: res.displacements, explain: 'Highest hot-spot stress at the nodes taken to the local analysis (zero elsewhere), on the deflected shape of the global model.' }),
      ],
      tables: [
        { title: 'Local joint results', columns: [{ key: 'n', label: 'Node' }, { key: 'b', label: 'Brace' }, { key: 'beta', label: 'β' }, { key: 'gamma', label: 'γ' }, { key: 'tau', label: 'τ' }, { key: 'th', label: 'θ', unit: 'deg' }, { key: 'sa', label: 'SCF axial chord saddle' }, { key: 'si', label: 'SCF in-plane chord crown' }, { key: 'so', label: 'SCF out-of-plane chord saddle' }, { key: 'nom', label: 'Nominal', unit: 'MPa' }, { key: 'hs', label: 'Hot spot', unit: 'MPa' }, { key: 'loc', label: 'Location' }, { key: 'eps', label: 'Plastic strain', unit: '%' }],
          rows: rows.map((r) => ({ n: r.j.node, b: r.j.brace, beta: round(r.h.parameters.beta, 3), gamma: round(r.h.parameters.gamma, 3), tau: round(r.h.parameters.tau, 3), th: round(r.j.thetaDeg, 3), sa: round(r.h.scf.axial.chordSaddle, 3), si: round(r.h.scf.ipb.chordCrown, 3), so: round(r.h.scf.opb.chordSaddle, 3), nom: round(r.nomS / 1e6), hs: round(r.hs / 1e6), loc: `${r.side}, point ${r.pt}`, eps: round(r.eps * 100, 3) })),
          explain: 'Geometry parameters, the three chord-side SCFs and the resulting stresses for each analysed joint, with the location of the peak around the weld.' },
        { title: 'Boundary forces for a detailed sub-model', columns: [{ key: 'n', label: 'Node' }, { key: 'b', label: 'Brace' }, { key: 'P', label: 'Brace axial (tension +)', unit: 'kN' }, { key: 'Mi', label: 'In-plane moment', unit: 'kN·m' }, { key: 'Mo', label: 'Out-of-plane moment', unit: 'kN·m' }, { key: 'Pc', label: 'Chord axial', unit: 'kN' }, { key: 'Mc', label: 'Chord moment', unit: 'kN·m' }],
          rows: rows.map((r) => ({ n: r.j.node, b: r.j.brace, P: round(r.L.P / 1e3), Mi: round(r.L.Mipb / 1e3), Mo: round(r.L.Mopb / 1e3), Pc: round(r.L.chord.P / 1e3), Mc: round(Math.hypot(r.L.chord.Mipb, r.L.chord.Mopb) / 1e3) })),
          explain: 'Member end forces at each analysed joint from the global model: the loads to apply at the cut boundaries of a shell or solid sub-model of that joint.' },
      ],
      notes: [`Governing joint: node ${gov.j.node}, brace ${gov.j.brace}: nominal ${fmt(gov.nomS / 1e6)} MPa → hot spot ${fmt(gov.hs / 1e6)} MPa on the ${gov.side} side at point ${gov.pt} (effective concentration ${fmt(gov.scfEff)}).`,
        `${rows.filter((r) => r.hs > fy).length} of ${rows.length} analysed joints exceed yield at the hot spot; limit ${fmt(lim / 1e6)} MPa. ${warn.length ? `${warn.length} joint(s) are outside the Efthymiou validity range: ${warn[0].h.warnings[0]}.` : 'All analysed joints are inside the Efthymiou validity range.'}`,
        `Wave by ${st.wf.label}, heading ${fmt(st.beta)}°.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Local stresses are acceptable: hot spot ${fmt(gov.hs / 1e6)} MPa ≤ ${fmt(lim / 1e6)} MPa, plastic strain ${fmt(gov.eps * 100)} %.`, gov.hs > fy ? 'The weld toe yields on first loading and then shakes down. Take the listed joints into the fatigue module, and build a shell sub-model with the boundary forces only for joints outside the parametric range or with stiffeners.' : 'The hot spots stay elastic. A sub-model is not needed for strength; use the same SCFs in the fatigue module.', 'DNV-RP-C203 (2021) §3.3; ISO 19902:2020 §16.10', `Lowering SCFs by geometry (larger β, thicker can) is more material-efficient than thicker braces: a can 5 mm thicker on one joint is about ${steelCo2(Math.PI * gov.j.D * 2.5 * gov.j.D * 0.005 * J.steel.rho / 1e3)}`)
          : adv('engineering', 'fail', `Local stress is too high at node ${gov.j.node}: hot spot ${fmt(gov.hs / 1e6)} MPa against ${fmt(lim / 1e6)} MPa, plastic strain ${fmt(gov.eps * 100)} %.`, `Thicken the chord can to about ${fmt(tReq * 1e3)} mm, increase the brace diameter towards β ≈ 0.6–0.8, or add ring stiffeners and confirm with a shell sub-model loaded by the boundary forces in the table.`, 'DNV-RP-C203 (2021) Appendix B; DNV-RP-C208 (2019) §4', `The thicker can is about ${steelCo2(Math.PI * gov.j.D * 2.5 * gov.j.D * Math.max(tReq - gov.j.T, 0.002) * J.steel.rho / 1e3)}`),
        adv('economic', 'info', `Only ${rows.length} of ${jd.length} joints need local attention at this load level; the effective concentration is ${fmt(gov.scfEff)}.`, 'Reserve shell-element sub-models for the few joints outside the parametric SCF range: each costs days of analyst time, while parametric SCFs screen the whole jacket in seconds. Offshore, a hot-spot grind or a clamp costs far more than a thicker can at fabrication.', BASIS_COST),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // T joint, θ = 90°, β = 0.5, γ = 20, τ = 0.5, α = 2·10/1.0 = 20 (no short-chord correction; the chord-fixity term carries sin²2θ = 0):
    // chord saddle SCF, axial = γ·τ^1.1·(1.11 − 3(β − 0.52)²)·sin^1.6θ = 20·0.46652·1.1088 = 10.345.
    { name: 'Efthymiou T-joint chord saddle SCF for axial load', source: 'Efthymiou (1988), OTJ\'88; DNV-RP-C203 Table B-1 Eq. (1)',
      inputs: { geom: testFrame('tjoint', { thetaDeg: 90 }), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 0.01, found: 'fixed', canT: 0, nSub: 1 },
      expect: [{ kpi: 'scfChordSaddle', value: 10.345, tol: 0.002 }] },
  ],
};
