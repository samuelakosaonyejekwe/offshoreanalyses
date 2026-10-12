// Buckling and stability of the dolphin: elastic critical load of the pile group as a sway frame on
// equivalent-fixity piles, Euler and ISO 19902 column strength, local buckling and beam-column interaction.
import { analyzeFrame } from '../../engine/structural/frame3d.js';
import { bucklingAnalysis } from '../../engine/structural/framebuckling.js';
import { isoTubularStrengths, isoTubularCheck } from '../../engine/structural/tubular.js';
import { need, run, fmt, round, cap, adv, check, ucCheck, linspace, pileInputs, soilInputs, layoutInputs, capInputs, capLoadInputs, capLoad, frameInputs, fixityLength, pileProps, profileFrom, pilesFrom, frameModel, framePlot, govUc, pileSteelT, steelCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

export default {
  id: 's06-buckling-stability', suite: 6, title: 'Structural buckling and stability',
  summary: 'Elastic critical load factor of the pile group, effective lengths, column and beam-column utilisation of the piles and local buckling of the tube wall.',
  covers: [9], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §13.2.3 (column and local buckling), §13.3.3 (beam-column)', 'EN 1993-1-1:2005 §5.2.1 (αcr and second-order effects)', 'API RP 2A-WSD 22nd ed. §6.10 (pile wall)', 'DNV-RP-C202 (2019) (shell buckling, reference)'],
  models: ['Euler buckling', 'Linear eigenvalue buckling of a 3-D frame (geometric stiffness)', 'Matrix stiffness frame analysis with Timoshenko elements', 'Beam-column theory with moment amplification', 'ISO 19902 column curve and local buckling strength', 'Equivalent depth of fixity from p–y pile-head stiffness'],
  limits: 'Piles are replaced by columns fixed at an equivalent depth derived from the secant pile-head stiffness; the soil gives no other lateral support, which is conservative for global buckling. Linear eigenvalue buckling with first-order member forces; imperfections and plasticity enter through the code column curve. Hydrostatic collapse and ring-stiffened shell buckling are not covered. ' + NOT_HIFI,
  publishes: ['dolphin.lambdaCr', 'dolphin.bucklingLength'],
  searchTerms: 'pile buckling unsupported length effective length sway frame tubular column',
  inputs: [
    ...layoutInputs(),
    ...pileInputs(),
    ...capInputs(),
    ...soilInputs(),
    ...capLoadInputs(),
    ...frameInputs(),
    { key: 'lamReq', label: 'Required elastic critical load factor', unit: '-', type: 'number', default: 3, min: 1, max: 20, group: 'Design factors', role: 'criteria', help: `Margin of the factored loads against elastic frame instability. EN 1993-1-1 §5.2.1 asks for second-order analysis below 10; this module amplifies the moments by 1/(1 − 1/αcr). ${CONFIRM}` },
  ],
  mesh: { key: 'nSub', label: 'Elements per pile', levels: [2, 4, 8], qoi: 'lambda1' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const piles = pilesFrom(inp), n = piles.length;
    const fixed = inp.headFix === 'fixed' && n > 1;
    const profile = inp.LeUser > 0 ? null : profileFrom(inp, P.D);
    const Le = fixityLength(inp, P, profile, fixed);
    const K = fixed ? 1 : 2;   // sway column: fixed–guided (K = 1) or fixed–free (K = 2)
    const load = capLoad(inp);
    need(load.Fz < 0, 'The net vertical load is upwards: no pile is in compression and there is nothing to buckle. Check the weight and uplift.');
    const fm = frameModel({ piles, P, Le, nSub: inp.nSub, fixedHead: fixed });
    fm.model.loadCases = [{ name: 'ref', nodal: [{ node: fm.centre, F: [load.Fx, load.Fy, load.Fz, load.Mx, load.My, load.Mz] }] }];
    const opt = { shearDeformation: !!inp.shear };
    const an = run(() => analyzeFrame(fm.model, opt), 'Frame analysis').results.ref;
    const bk = run(() => bucklingAnalysis(fm.model, { loadCase: 'ref', nModes: 1, ...opt }), 'Buckling analysis');
    const lam = bk.factors.find((f) => f > 0);
    need(Number.isFinite(lam) && lam > 0, 'No positive buckling factor was found: the reference loading puts no pile in compression.');
    const amp = lam > 1 ? 1 / (1 - 1 / lam) : 99;
    const S = isoTubularStrengths({ D: P.D, t: P.t, fy: P.fy, E: P.E, L: Le, K });
    const PE = Math.PI ** 2 * P.EI / (K * Le) ** 2;
    // per pile: axial force at the head and the largest first-order moment along the pile
    const rows = piles.map((p, k) => {
      let N = 0, M = 0, V = 0;
      for (let e = 0; e < inp.nSub; e++) {
        const mf = an.memberForces[k * inp.nSub + e];
        for (const end of [mf.endA, mf.endB]) { M = Math.max(M, Math.hypot(end.My, end.Mz)); V = Math.max(V, Math.hypot(end.Vy, end.Vz)); if (Math.abs(end.N) > Math.abs(N)) N = end.N; }
      }
      const Nc = Math.max(-N, 0);
      const chk = isoTubularCheck({ D: P.D, t: P.t, fy: P.fy, E: P.E, L: Le, Ky: K, Kz: K, Cmy: 0.85, Cmz: 0.85, N, My: M, V });
      const g = govUc(chk);
      return { id: k + 1, N, Nc, M, V, ucCol: Math.min(1.18 * Nc / S.Pc, 99), ucBC: g.uc, gov: g.name };
    });
    const gov = rows.reduce((a, r) => (r.ucBC > a.ucBC ? r : a), rows[0]);
    const NcMax = Math.max(...rows.map((r) => r.Nc));
    const Dt = P.D / P.t, par = P.fy * P.D / (P.E * P.t);
    const checks = [
      { name: 'Elastic critical load factor of the frame', demand: round(inp.lamReq, 4), capacity: round(cap(lam, 1e4), 4), unit: 'factor (required / achieved)', uc: Math.min(inp.lamReq / lam, 99), ref: 'EN 1993-1-1:2005 §5.2.1; linear eigenvalue buckling', explain: 'Multiplier on the factored loads at which the pile group buckles sideways as an elastic frame. Close to 1 the dolphin is unstable; below about 10 second-order moments are significant.' },
      ucCheck(`Column buckling, pile ${rows.reduce((a, r) => (r.ucCol > a.ucCol ? r : a), rows[0]).id}`, Math.max(...rows.map((r) => r.ucCol)), 'ISO 19902:2020 §13.2.3 (γR,c = 1.18)', 'Axial compression against the column buckling strength over the unsupported length to fixity. If exceeded the pile bows sideways and sheds its load.'),
      ucCheck(`Beam-column interaction, pile ${gov.id} (${gov.gov || 'combined'})`, gov.ucBC, 'ISO 19902:2020 §13.3.3', 'Compression with bending amplified by the P–δ effect, against the member stability and local strength expressions.'),
      check('Wall slenderness D/t', Dt, 120, '-', 'ISO 19902:2020 §13.2 (validity limit D/t ≤ 120); API RP 2A-WSD §6.10.6', 'Diameter-to-thickness ratio. Thin walls buckle locally before the section yields and are easily damaged in driving.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    // column curve for the plot
    const lams = linspace(0.05, 3, 40), fcRel = lams.map((l) => round(l <= 1.34 ? 1 - 0.278 * l * l : 0.9 / (l * l)));
    const tMinApi = (6.35 + P.D * 1e3 / 100);
    const steelT = pileSteelT(inp, n);
    const mode = bk.modes[bk.factors.indexOf(lam)] ?? bk.modes[0];
    let mx = 0; for (let i = 0; i < fm.model.nodes.length; i++) mx = Math.max(mx, Math.hypot(mode[6 * i], mode[6 * i + 1], mode[6 * i + 2]));
    const scale = 0.15 * Le / (mx || 1);
    const be = bk.effectiveLengths.filter((e) => e.N < 0 && e.member < n * inp.nSub).sort((a, b) => a.N - b.N)[0];
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Elastic critical factor ${fmt(lam)} (moment amplification ${fmt(amp)}); largest pile compression ${fmt(NcMax / 1e3)} kN against a column strength of ${fmt(S.Pc / 1e3)} kN over ${fmt(Le)} m to fixity (K = ${K}); D/t = ${fmt(Dt)}.`,
        ok ? `Stability is adequate. Re-check with scour at its design depth and with the cyclic soil stiffness, both of which lengthen the unsupported pile; ${par > 0.0517 ? 'the wall is too slender for full plastic bending (reduced bending strength applies).' : 'the wall is compact enough for full plastic bending.'}`
          : (w.name.startsWith('Wall') ? `Increase the wall to at least ${fmt(P.D * 1e3 / 120 + inp.corr)} mm (D/t ≤ 120); API RP 2A minimum for driving is ${fmt(tMinApi)} mm.` : `Shorten the unsupported length (scour protection, bracing below the cap or raked piles), or increase the diameter: Euler load grows with D³·t. Required elastic critical load per pile is about ${fmt(inp.lamReq * NcMax / 1e3)} kN against ${fmt(PE / 1e3)} kN.`),
        'ISO 19902:2020 §13.2.3 and §13.3.3; EN 1993-1-1:2005 §5.2',
        `Stability is governed by stiffness, not strength: a higher steel grade does not help, while diameter does. The piles hold ${steelCo2Text(steelT)}`),
      adv('economic', ok ? 'info' : 'warn', `Unsupported length ${fmt(Le)} m, slenderness λ = ${fmt(S.lambda)}; column strength is ${fmt(S.fc / P.fy * 100)} % of yield.`,
        S.lambda > 1 ? 'The piles work at a low fraction of yield in compression: steel is being bought for stiffness. Compare raked piles or a submerged bracing frame (more fabrication, less tonnage) with thicker plain piles.' : 'Column buckling costs little strength here; wall thickness can be optimised on bending and fatigue rather than on stability.',
        ECON_BASIS),
    ];
    return {
      kpis: [
        { key: 'lambda1', label: 'Elastic critical load factor', unit: '-', value: lam, limit: inp.lamReq, sense: 'min', explain: 'Factor by which the applied factored loads must be multiplied to cause elastic sway buckling of the pile group. It measures global stability and the size of second-order effects.' },
        { key: 'amp', label: 'Second-order moment amplification', unit: '-', value: Math.min(amp, 99), explain: '1/(1 − 1/αcr): how much the first-order sway moments grow because the vertical load acts on the displaced cap. Above about 1.1 second-order effects must be designed for.' },
        { key: 'Le', label: 'Length from pile head to fixity', unit: 'm', value: Le, explain: 'Unsupported column length: free length in water plus the depth of apparent fixity in the soil at the working load.' },
        { key: 'PE', label: 'Euler load of one pile', unit: 'kN', value: PE / 1e3, explain: 'π²EI/(K·L)² for the sway mode of a single pile. An upper bound that ignores imperfections and yielding.' },
        { key: 'Pc', label: 'Column buckling strength of one pile', unit: 'kN', value: S.Pc / 1e3, explain: 'Representative axial strength from the ISO 19902 column curve, including local buckling of the wall. Divide by 1.18 for the design resistance.' },
        { key: 'NcMax', label: 'Largest pile compression', unit: 'kN', value: NcMax / 1e3, explain: 'Greatest axial compression in any pile under the factored loads (weight plus push–pull from the horizontal load).' },
        { key: 'ucCol', label: 'Column buckling utilisation', unit: '-', value: Math.max(...rows.map((r) => r.ucCol)), limit: 1, sense: 'max', explain: 'Factored compression over design column strength for the most loaded pile.' },
        { key: 'ucBC', label: 'Beam-column utilisation', unit: '-', value: gov.ucBC, limit: 1, sense: 'max', explain: 'ISO 19902 interaction of compression and amplified bending for the most loaded pile. The check that normally governs a sway dolphin.' },
        { key: 'slender', label: 'Column slenderness λ', unit: '-', value: S.lambda, explain: 'Reduced slenderness (K·L/πr)·√(fyc/E). Below about 0.4 buckling is negligible; above 1.34 the pile buckles elastically.' },
        { key: 'Dt', label: 'Diameter / thickness', unit: '-', value: Dt, limit: 120, sense: 'max', explain: 'Wall slenderness of the corroded tube. Above about 60 local buckling reduces the axial strength; 120 is the limit of the code expressions.' },
        { key: 'localRed', label: 'Local buckling strength / yield', unit: '-', value: S.fyc / P.fy, explain: 'Ratio of the representative local buckling strength of the wall to yield. 1.0 means the wall is stocky enough to reach yield in compression.' },
        { key: 'Keff', label: 'Effective length factor from the eigenvalue', unit: '-', value: be ? Math.max(be.Ley, be.Lez) / Le : K, explain: 'Effective length of the most compressed pile back-calculated from the frame buckling mode, over its length to fixity. It reflects cap rotation and load sharing; compare with the idealised 1.0 (fixed–guided) or 2.0 (cantilever).' },
      ],
      checks,
      plots: [
        framePlot(fm.model, { disp: mode, scale, title: 'First buckling mode of the pile group', explain: 'Equivalent-fixity frame with the first buckling mode overlaid (arbitrary amplitude). A sideways sway of the cap with the piles bending in double curvature is the expected mode; twisting indicates an unbalanced layout.' }),
        { type: 'line', title: 'Column curve and pile operating point', x: { label: 'Slenderness λ', unit: '-' }, y: { label: 'Stress / local buckling strength', unit: '-' },
          series: [{ name: 'ISO 19902 column strength fc/fyc', x: lams.map((l) => round(l)), y: fcRel }, { name: 'Euler', x: lams.filter((l) => l > 0.75).map((l) => round(l)), y: lams.filter((l) => l > 0.75).map((l) => round(1 / (l * l))), dash: true }, { name: 'Most loaded pile', x: [round(S.lambda)], y: [round(NcMax / P.A / S.fyc)] }],
          explain: 'Design column strength against slenderness with the elastic Euler curve. The marker is the axial stress of the most compressed pile: its distance below the curve is the reserve before bending is added.' },
        { type: 'bar', title: 'Pile compression against column strength', categories: rows.map((r) => `Pile ${r.id}`), series: [{ name: 'Compression (kN)', y: rows.map((r) => round(r.Nc / 1e3)) }, { name: 'Design column strength (kN)', y: rows.map(() => round(S.Pc / 1.18 / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Axial compression in each pile next to the design buckling resistance. Piles in tension show zero.' },
      ],
      tables: [{ title: 'Pile stability checks', columns: [{ key: 'id', label: 'Pile' }, { key: 'N', label: 'Axial (tension +)', unit: 'kN' }, { key: 'M', label: 'First-order moment', unit: 'kN·m' }, { key: 'ucCol', label: 'Column UC' }, { key: 'ucBC', label: 'Beam-column UC' }, { key: 'gov', label: 'Governing expression' }],
        rows: rows.map((r) => ({ id: r.id, N: round(r.N / 1e3), M: round(r.M / 1e3), ucCol: round(r.ucCol), ucBC: round(r.ucBC), gov: r.gov })),
        explain: 'One row per pile: first-order forces from the frame analysis and the ISO 19902 column and beam-column utilisations.' }],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, `Length to fixity ${inp.LeUser > 0 ? 'as specified' : 'from the secant pile-head stiffness at the reference shear'}: ${fmt(Le)} m.`, ...S.warnings.slice(0, 2), bk.converged ? `Eigenvalue iteration converged in ${bk.iterations} iterations.` : 'Eigenvalue iteration did not fully converge: refine the mesh.'],
      advice,
      publish: { 'dolphin.lambdaCr': lam, 'dolphin.bucklingLength': K * Le },
      quantities: { steel_t: round(steelT) },
    };
  },
  verify: [
    // Single cantilever pile, L = 20 m, EI = 210e9·π/64·(1 − 0.96⁴) = 1.5530e9 N·m²: Pcr = π²EI/(2L)² = 9579.6 kN;
    // applied compression 1000 kN (cap 4×4×2.5 m at 25 kN/m³) → factor 9.580.
    { name: 'Euler buckling of a cantilever: Pcr = π²EI/(4L²)', source: 'Euler (1744); Timoshenko & Gere, Theory of Elastic Stability §2.1',
      inputs: { piles: [{ x: 0, y: 0, batter: 0, azimuth: 0 }], D: 1, t: 20, corr: 0, steel: 's355', LeUser: 20, shear: false, nSub: 8, capB: 4, capL: 4, capT: 2.5, gammaConc: 25, Wequip: 0, Fx: 0, Fy: 0, Fup: 0, Mz: 0, gammaF: 1, gammaG: 1 },
      expect: [{ kpi: 'lambda1', value: 9.5796, tol: 0.002 }, { kpi: 'PE', value: 9579.6, tol: 1e-3 }] },
    // Four vertical piles 20 m apart, heads fixed in a rigid cap: sway buckling as fixed–guided columns, Pcr = π²EI/L² each
    // = 38 318 kN; total vertical load 4000 kN (1000 kN per pile) → factor 38.3 (pile axial flexibility lets the cap rock slightly: −0.4 %).
    { name: 'Sway buckling of a rigid-cap group: Pcr = π²EI/L² per pile', source: 'Timoshenko & Gere, Theory of Elastic Stability (fixed–guided column)',
      inputs: { piles: [{ x: -10, y: -10, batter: 0, azimuth: 0 }, { x: 10, y: -10, batter: 0, azimuth: 0 }, { x: 10, y: 10, batter: 0, azimuth: 0 }, { x: -10, y: 10, batter: 0, azimuth: 0 }], headFix: 'fixed', D: 1, t: 20, corr: 0, steel: 's355', LeUser: 20, shear: false, nSub: 8, capB: 8, capL: 8, capT: 2.5, gammaConc: 25, Wequip: 0, Fx: 0, Fy: 0, Fup: 0, Mz: 0, gammaF: 1, gammaG: 1 },
      expect: [{ kpi: 'lambda1', value: 38.318, tol: 0.01 }] },
  ],
};
