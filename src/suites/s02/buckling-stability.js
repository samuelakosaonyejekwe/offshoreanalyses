// Stability of the topsides: elastic frame buckling and second-order (P–Δ) amplification of the deck frame,
// effective lengths and buckling resistance of the columns, plate and stiffened-panel buckling of the deck
// plating acting as a girder flange, and lateral-torsional buckling of the girders.
import { bucklingAnalysis, pDeltaAnalysis } from '../../engine/structural/framebuckling.js';
import { refineModel } from '../../engine/structural/frame3d.js';
import { ec3Chi, elasticCriticalMoment } from '../../engine/structural/steelcheck.js';
import { plateUltimateStrength, stiffenerColumnBuckling, stiffenerTripping, plateBucklingCoefficient } from '../../engine/structural/plates.js';
import { deckInputs, solverInput, buildDeck, twistRestrained, loadCase, equipmentOf, framePlot, need, run, fmt, round, cap, adv, check, fsCheck, steelCo2Text, CONFIRM, DECK_LIMITS } from './_common.js';

/**
 * Members for the frame stability eigenproblem. On top of the twist restraint, a welded deck plate is an in-plane
 * diaphragm: a girder or beam under it cannot bow sideways between grid points. Without that restraint the lowest
 * eigenmode of a braced deck is weak-axis bowing of the compressed chord girder of the perimeter truss, which the
 * plating prevents. The plate strip tributary to the member (width b, thickness t) is added to the weak-axis second
 * moment of area as t b³/12. Grated decks (plate thickness 0) get no such restraint and keep the mode.
 */
function stabilityMembers(deck, tPlate) {
  return twistRestrained(deck).map((m, i) => {
    const o = deck.orig[i], b = deck.trib[i];
    return tPlate > 0 && o.horizontal && o.kind === 'I' && b > 0 ? { ...m, section: { ...m.section, Iz: o.sec.Iz + (tPlate * b ** 3) / 12 } } : m;
  });
}

export default {
  id: 's02-buckling-stability', suite: 2, title: 'Buckling and stability',
  summary: 'Elastic critical load factor and sway amplification of the deck frame, column effective lengths and buckling resistance, deck plate and stiffened-panel buckling, and lateral-torsional buckling of girders.',
  covers: [5], fidelity: 'numerical',
  standards: ['EN 1993-1-1:2005 §5.2 (αcr and second-order effects), §6.3 (member buckling)', 'EN 1993-1-5:2006 §4, §5, §10 (plate buckling)', 'DNV-RP-C201 (2010) §6–§7 (stiffened panels)', 'ISO 19902:2020 §13.5 (effective lengths)', 'NORSOK N-004 (2022) §6.3'],
  models: ['Geometric stiffness formulation', 'Linear (eigenvalue) frame stability analysis', 'Second-order analysis (P–Δ and P–δ)', 'Euler buckling', 'Lateral–torsional buckling', 'Plate buckling (effective width)', 'Perry–Robertson column curves', 'Stiffener tripping'],
  limits: `${DECK_LIMITS} Frame buckling is linear eigenvalue analysis on the first-order axial forces; imperfections enter only through the code buckling curves. The deck plate is taken as an in-plane diaphragm that stops plated girders and beams twisting or bowing sideways, so αcr is the lowest mode of the legs, columns, braces and frame sway. Nonlinear post-buckling and imperfection-sensitivity analyses are outside this module.`,
  searchTerms: 'frame buckling critical load factor P-delta effective length stiffened panel buckling offshore deck',
  inputs: [
    ...deckInputs(),
    { key: 'swFac', label: 'Self-weight factor', unit: '-', type: 'number', default: 1, min: 0, max: 1.5, group: 'Permanent and equipment loads', role: 'load', help: 'Multiplies the structural self-weight inside the factored combination; 0 leaves only equipment and live load.' },
    { key: 'Hx', label: 'Lateral load along x (factored)', unit: 'kN', type: 'number', default: 1200, min: 0, max: 5e4, group: 'Environmental and interface loads', role: 'bc', help: 'Used for the second-order sway amplification.' },
    { key: 'gG', label: 'Partial factor on permanent actions', unit: '-', type: 'number', default: 1.3, min: 0.8, max: 1.6, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gQ', label: 'Partial factor on variable actions', unit: '-', type: 'number', default: 1.5, min: 0, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'alphaReq', label: 'Required critical load factor αcr', unit: '-', type: 'number', default: 10, min: 1.5, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `EN 1993-1-1 §5.2.1: first-order analysis is adequate for αcr ≥ 10 (elastic); below it use the second-order results. ${CONFIRM}` },
    { key: 'gM1', label: 'γM1 (buckling resistance)', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: `1.15 NORSOK / DNV; 1.0 EN 1993-1-1 recommended value. ${CONFIRM}` },
    { key: 'colCurve', label: 'Column buckling curve', type: 'select', default: 'a', group: 'Acceptance criteria', role: 'criteria', options: ['a0', 'a', 'b', 'c', 'd'].map((v) => ({ value: v, label: `Curve ${v}` })), help: 'EN 1993-1-1 Table 6.2: hot-finished hollow sections a; cold-formed or welded tubes c.' },
    { key: 'sStiff', label: 'Deck stiffener spacing', unit: 'mm', type: 'number', default: 700, min: 150, max: 3000, group: 'Deck plate panel', role: 'geometry', design: true },
    { key: 'Lpanel', label: 'Stiffener span between beams', unit: 'm', type: 'number', default: 3.3, min: 0.3, max: 15, group: 'Deck plate panel', role: 'geometry' },
    { key: 'hw', label: 'Stiffener web height', unit: 'mm', type: 'number', default: 240, min: 40, max: 1200, group: 'Deck plate panel', role: 'geometry' },
    { key: 'tw', label: 'Stiffener web thickness', unit: 'mm', type: 'number', default: 10, min: 4, max: 60, group: 'Deck plate panel', role: 'geometry' },
    { key: 'bf', label: 'Stiffener flange width (0 = flat bar)', unit: 'mm', type: 'number', default: 100, min: 0, max: 600, group: 'Deck plate panel', role: 'geometry' },
    { key: 'tf', label: 'Stiffener flange thickness', unit: 'mm', type: 'number', default: 12, min: 0, max: 60, group: 'Deck plate panel', role: 'geometry' },
    { key: 'sigX', label: 'In-plane compressive stress in the deck plate', unit: 'MPa', type: 'number', default: 110, min: 0, max: 600, group: 'Deck plate panel', role: 'bc', help: 'Factored membrane compression along the stiffeners when the plate acts as the flange of a girder (from the in-place or local stress analysis).' },
    { key: 'tauXY', label: 'In-plane shear stress in the deck plate', unit: 'MPa', type: 'number', default: 25, min: 0, max: 400, group: 'Deck plate panel', role: 'bc' },
    { key: 'latFrac', label: 'Laterally unrestrained fraction of a girder bay', unit: '-', type: 'number', default: 0.25, min: 0.02, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'Spacing of lateral restraints to the compression flange as a fraction of the bay length.' },
    solverInput(4),
  ],
  // One element per member cannot bow between its joints, so the ladder starts at two; a member held at both ends gains
  // little from two to four elements and most from four to eight, so the asymptotic range is the last three levels.
  mesh: { key: 'nSub', label: 'Elements per member', levels: [2, 4, 8, 16], qoi: 'alphaCr' },
  compute(inp, ctx) {
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), st = deck.steel;
    need(inp.tPlate > 0, 'The plate buckling checks need a deck plate thickness above zero.');
    const zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0), extra = [];
    deck.levelNodes.forEach((ns, k) => { for (const n of ns) extra.push({ node: n, F: [(inp.Hx * 1e3 * zs[k]) / zsum / ns.length, 0, 0, 0, 0, 0] }); });
    const grav = loadCase(deck, inp, { name: 'gravity', gDead: inp.gG * inp.swFac, gEquip: inp.gG, gLive: inp.gQ, eqs });
    const comb = loadCase(deck, inp, { name: 'design', gDead: inp.gG * inp.swFac, gEquip: inp.gG, gLive: inp.gQ, extra, eqs });
    const ref = run(() => refineModel({ ...deck.model, members: stabilityMembers(deck, inp.tPlate / 1e3), loadCases: [grav, comb] }, inp.nSub), 'Mesh refinement');
    const bk = run(() => bucklingAnalysis(ref.model, { loadCase: 'gravity', nModes: 3 }), 'Frame buckling analysis (check that the frame is stable and carries compression)');
    const pos = bk.factors.filter((f) => f > 0 && Number.isFinite(f));
    need(pos.length, 'No positive buckling load factor was found: the frame carries no compression under the gravity combination.');
    const alphaCr = Math.min(pos[0], 9999);
    // Second-order sway amplification under the combination with lateral load.
    let amp = 1, pdOk = true;
    try {
      const pd = pDeltaAnalysis(ref.model, { loadCase: 'design' }); pdOk = pd.converged;
      // Sway amplification: largest horizontal displacement, second order over first order (the engine's own ratio
      // is on the largest displacement of any kind, which gravity sag dominates in a braced deck).
      const lat = (u) => { let m = 0; for (let i = 0; i < u.length / 6; i++) m = Math.max(m, Math.hypot(u[6 * i], u[6 * i + 1])); return m; };
      const l1 = lat(pd.firstOrder.displacements), l2 = lat(pd.result.displacements);
      amp = l1 > 1e-9 ? l2 / l1 : pd.amplification;
    } catch { pdOk = false; amp = alphaCr > 1 ? 1 / (1 - 1 / alphaCr) : 99; }
    const ampEst = alphaCr > 1 ? 1 / (1 - 1 / alphaCr) : 99;
    // Columns: effective length from the critical force of each compressed vertical member, K capped at 5.
    const cols = [];
    deck.orig.forEach((o, i) => {
      if (o.horizontal || o.kind !== 'tube') return;
      let N = 0; for (const k of ref.memberMap[i]) N = Math.min(N, bk.axial[k]);
      if (N >= -1) return;
      const Pcr = alphaCr * -N, Le = Math.PI * Math.sqrt(st.E * o.sec.Iy / Pcr), K = cap(Le / o.L, 5);
      const Npl = o.sec.A * st.fy, lam = Math.sqrt(Npl / (Math.PI ** 2 * st.E * o.sec.Iy / (K * o.L) ** 2)), chi = ec3Chi(lam, inp.colCurve), NbRd = chi * Npl / inp.gM1;
      cols.push({ i, group: o.group, N: -N, L: o.L, K, lam, chi, NbRd, uc: -N / NbRd });
    });
    need(cols.length, 'No compressed column, leg or brace was found under the gravity combination.');
    const gc = cols.reduce((a, c) => (c.uc > a.uc ? c : a));
    // Deck plate between stiffeners and the stiffened panel as a column (plate acting as the girder flange).
    const b = inp.sStiff / 1e3, a = inp.Lpanel, t = inp.tPlate / 1e3, sx = inp.sigX * 1e6, tau = inp.tauXY * 1e6;
    need(a >= b, 'The stiffener span must not be shorter than the stiffener spacing.');
    const pl = run(() => plateUltimateStrength({ a, b, t, E: st.E, nu: st.nu, fy: st.fy, sigmaX: sx, tau }), 'Plate buckling');
    const kc = plateBucklingCoefficient({ aspect: a / b }).k, sigE = Math.PI ** 2 * st.E / (12 * (1 - st.nu ** 2)) * (t / b) ** 2;
    const sc = run(() => stiffenerColumnBuckling({ L: a, s: b, tp: t, hw: inp.hw / 1e3, tw: inp.tw / 1e3, bf: inp.bf / 1e3, tf: inp.tf / 1e3, E: st.E, nu: st.nu, fy: st.fy }), 'Stiffened panel buckling');
    const trip = run(() => stiffenerTripping({ L: a, hw: inp.hw / 1e3, tw: inp.tw / 1e3, bf: inp.bf / 1e3, tf: inp.tf / 1e3, E: st.E, nu: st.nu, fy: st.fy }), 'Stiffener tripping');
    // Girder lateral-torsional buckling between restraints.
    const g = deck.sections.girder, Llt = inp.latFrac * deck.dx;
    const Mcr = elasticCriticalMoment({ E: st.E, G: st.G, Iz: g.Iz, It: g.J, Iw: g.Iw, L: Llt }), lamLT = Math.sqrt(g.Wply * st.fy / Mcr);
    const checks = [
      fsCheck('Elastic frame stability (αcr)', alphaCr, inp.alphaReq, 'EN 1993-1-1 §5.2.1 Eq. (5.1)', 'Ratio of the elastic critical load of the frame to the design load. Below the required value, sway (P–Δ) effects must be included in the member forces; below about 3 the frame is close to elastic instability.'),
      check(`Column flexural buckling — ${gc.group} ${gc.i + 1}`, gc.N / 1e3, gc.NbRd / 1e3, 'kN', `EN 1993-1-1 §6.3.1, curve ${inp.colCurve}, K = ${fmt(gc.K)} from the frame buckling mode`, 'Axial compression against the buckling resistance with the effective length of the governing sway or braced mode. Exceeding it means the column bows and sheds load.'),
      check('Deck plate buckling between stiffeners', pl.utilisation * inp.gM1, 1, 'interaction', 'EN 1993-1-5 §4.4, §5, §10 Eq. (10.5), resistance / γM1', 'Compression and shear in the plate panel against its ultimate (post-buckling) strength. Above 1.0 the plate buckles between stiffeners and sheds its share of the girder flange force.'),
      check('Stiffened panel column buckling', sx / 1e6, sc.sigmaU / inp.gM1 / 1e6, 'MPa', 'DNV-RP-C201 §7; EN 1993-1-5 §4.5.3 (Perry–Robertson, effective plate flange)', 'Buckling of the stiffener with its effective plate between the deck beams. Failure removes the plate from the girder flange over a whole bay.'),
      check('Stiffener tripping', sx / 1e6, trip.sigmaTinelastic / inp.gM1 / 1e6, 'MPa', 'DNV-RP-C201 §7.5 (torsional buckling of stiffeners)', 'Sideways twisting of the stiffener about its weld line under compression. It is sudden and leaves the plate unsupported; add tripping brackets if it governs.'),
      check('Girder lateral-torsional slenderness', lamLT, 0.4, 'λ̄LT', 'EN 1993-1-1 §6.3.2.2(4): no LTB reduction for λ̄LT ≤ 0.4', 'Non-dimensional slenderness of the girder between lateral restraints. Above 0.4 the bending resistance must be reduced for lateral-torsional buckling (done in the in-place check).'),
    ];
    const hard = checks.filter((c, i) => i !== 0 && i !== 5), gov = hard.reduce((p, c) => (c.uc > p.uc ? c : p)), ok = hard.every((c) => c.uc <= 1) && alphaCr >= 3;
    const mode = bk.modes[bk.factors.indexOf(pos[0])] || bk.modes[0], nN = ref.model.nodes.length;
    const mval = new Array(nN); for (let i = 0; i < nN; i++) mval[i] = Math.hypot(mode[6 * i], mode[6 * i + 1], mode[6 * i + 2]);
    const colT = deck.orig.filter((o) => !o.horizontal).reduce((s, o) => s + o.sec.A * o.L * st.rho, 0) / 1e3;
    const lamP = linspaceLam();
    return {
      kpis: [
        { key: 'alphaCr', label: 'Elastic critical load factor αcr', unit: '-', value: alphaCr, limit: inp.alphaReq, sense: 'min', explain: 'Multiplier on the factored gravity loads at which the frame buckles elastically. Above 10 sway effects are negligible; between 3 and 10 amplify the sway moments; below 3 a full second-order analysis is essential.' },
        { key: 'amp', label: 'Second-order sway amplification', unit: '-', value: amp, explain: `Largest horizontal displacement from the P–Δ analysis over the first-order value (the estimate 1/(1 − 1/αcr) is ${fmt(ampEst)}). Multiply sway moments by this factor.` },
        { key: 'Kcol', label: 'Effective length factor of the governing column', unit: '-', value: gc.K, explain: 'Back-calculated from the frame buckling mode: Le = π√(EI/Pcr). About 0.7–1.0 for braced frames, above 1 for sway frames. Capped at 5 for lightly loaded members.' },
        { key: 'ucCol', label: 'Column buckling utilisation', unit: '-', value: round(gc.uc, 5), limit: 1, sense: 'max', explain: 'Design compression over the flexural buckling resistance of the governing leg, column or brace.' },
        { key: 'sigCrPlate', label: 'Elastic plate buckling stress', unit: 'MPa', value: kc * sigE / 1e6, explain: 'Critical compressive stress of the simply supported plate panel between stiffeners, k π²E/(12(1−ν²)) (t/b)².' },
        { key: 'ucPlate', label: 'Plate buckling utilisation', unit: '-', value: round(pl.utilisation * inp.gM1, 5), limit: 1, sense: 'max', explain: 'Combined compression and shear against the ultimate plate strength of EN 1993-1-5 with γM1.' },
        { key: 'sigUpanel', label: 'Stiffened panel ultimate stress', unit: 'MPa', value: sc.sigmaU / 1e6, explain: 'Average compressive stress at which the stiffener with its effective plate buckles as a column between the deck beams.' },
        { key: 'lamLT', label: 'Girder LTB slenderness λ̄LT', unit: '-', value: lamLT, explain: 'Slenderness for lateral-torsional buckling between restraints; below 0.4 the girder reaches its full plastic moment.' },
      ],
      checks,
      plots: [
        framePlot({ title: 'First buckling mode', ref, u: mode, values: mval, valueLabel: 'Mode amplitude (normalised)', valueUnit: '-', explain: 'Shape of the lowest elastic buckling mode. A sway of whole storeys points to missing bracing; a single bowing member points to a slender column or brace.' }),
        { type: 'line', title: 'Column buckling curve and operating points', x: { label: 'Slenderness λ̄', unit: '-' }, y: { label: 'χ = Nb / (A fy)', unit: '-' },
          series: [{ name: `Curve ${inp.colCurve}`, x: lamP, y: lamP.map((l) => round(ec3Chi(l, inp.colCurve))) }, { name: 'Euler', x: lamP.filter((l) => l >= 1), y: lamP.filter((l) => l >= 1).map((l) => round(1 / (l * l))), dash: true }, { name: 'Compressed members (N γM1 / A fy)', x: cols.map((c) => round(c.lam)), y: cols.map((c) => round(c.N * inp.gM1 / (deck.orig[c.i].sec.A * st.fy))), mode: 'markers' }],
          explain: 'Each compressed leg, column and brace plotted at its slenderness and axial stress ratio. Points must lie below the buckling curve; the gap between the curve and the Euler hyperbola is the allowance for imperfections.' },
        { type: 'bar', title: 'Stability utilisation', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => round(c.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'All stability checks side by side; the frame, member and plate modes interact, so keep each one comfortably below 1.0.' },
      ],
      tables: [{ title: 'Compressed members', columns: [{ key: 'm', label: 'Member' }, { key: 'group', label: 'Group' }, { key: 'N', label: 'Compression', unit: 'kN' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'K', label: 'K' }, { key: 'lam', label: 'λ̄' }, { key: 'chi', label: 'χ' }, { key: 'Nb', label: 'Nb,Rd', unit: 'kN' }, { key: 'uc', label: 'UC' }],
        rows: cols.sort((p, q) => q.uc - p.uc).slice(0, 40).map((c) => ({ m: c.i + 1, group: c.group, N: round(c.N / 1e3), L: round(c.L), K: round(c.K), lam: round(c.lam), chi: round(c.chi), Nb: round(c.NbRd / 1e3), uc: round(c.uc) })),
        explain: 'Compressed legs, columns and braces sorted by utilisation, with the effective length factor from the frame buckling mode. Lightly loaded members show large K values because the frame, not the member, buckles first.' }],
      notes: [`αcr = ${fmt(alphaCr)}: ${alphaCr >= 10 ? 'first-order analysis is adequate' : alphaCr >= 3 ? 'sway effects must be included (amplified moments or second-order analysis)' : 'the frame is close to elastic instability'}.`,
        `Second-order amplification ${fmt(amp)}${pdOk ? '' : ' (P–Δ iteration did not converge: estimate from αcr)'}; ${gov.name.toLowerCase()} governs the member and plate checks at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', alphaCr < inp.alphaReq || gov.uc > 0.9 ? 'warn' : 'ok', `Frame αcr = ${fmt(alphaCr)} (required ${inp.alphaReq}); ${gov.name.toLowerCase()} at ${fmt(gov.uc)}; plate panel elastic buckling stress ${fmt(kc * sigE / 1e6)} MPa.`,
          alphaCr < inp.alphaReq ? `Amplify sway moments by ${fmt(amp)} in the member checks, or raise αcr above ${inp.alphaReq} with perimeter bracing.` : `First-order member forces are adequate. Use K = ${fmt(gc.K)} for the governing column in the in-place member checks.`,
          'EN 1993-1-1:2005 §5.2.1–5.2.2; EN 1993-1-5:2006 §4', `Legs, columns and braces hold ${steelCo2Text(colT)} Bracing raises αcr with far less steel than heavier columns.`)
          : adv('engineering', 'fail', alphaCr < 3 ? `The frame is close to elastic instability: αcr = ${fmt(alphaCr)}.` : `${gov.name} fails at utilisation ${fmt(gov.uc)}.`,
            alphaCr < 3 ? 'Add vertical bracing or fix the leg bases; then re-check.' : gov.name.startsWith('Column') ? `Increase the column radius of gyration (larger diameter) or brace it at mid-height; area needed is about ${fmt((gov.uc - 1) * 100, 2)} % more at the same slenderness.` : gov.name.startsWith('Deck plate') ? `Increase the plate to about ${fmt(inp.tPlate * Math.sqrt(gov.uc), 2)} mm or reduce the stiffener spacing to about ${fmt(inp.sStiff / Math.sqrt(gov.uc), 3)} mm.` : 'Deepen the stiffener web, add a flange or fit tripping brackets at mid-span of the stiffeners.',
            'EN 1993-1-1:2005 §5.2, §6.3; EN 1993-1-5:2006 §4; DNV-RP-C201 (2010) §7', `Closer stiffeners or tripping brackets restore panel stability with a few kg/m²; 1 mm of extra plate over ${fmt(deck.planArea * zs.length)} m² is ${steelCo2Text(deck.planArea * zs.length * 7.85 / 1e3)}`),
        adv('economic', 'info', `Stability is controlled by ${fmt(colT)} t of legs, columns and braces and by the stiffener layout of ${fmt(deck.planArea * zs.length)} m² of plated deck.`,
          'Decide braced versus moment frame early: bracing saves column and girder tonnage but blocks equipment access and escape routes, which is a layout cost rather than a steel cost.', 'Concept-selection practice; AACE RP 18R-97 estimate classes'),
      ],
    };
  },
  verify: [
    { // Four slender columns (CHS 300 × 10, h = 20 m, fixed bases) under a stiff grid, 100 t at the plan centre → P = 245.17 kN each.
      // Stiff girders restrain the column heads against rotation, so the sway mode has K = 1:
      // I = π(0.3⁴ − 0.28⁴)/64 = 9.5896e-5 m⁴; Pcr = π²·210e9·9.5896e-5/20² = 496.89 kN; αcr = 496.89/245.17 = 2.027.
      name: 'Sway buckling of fixed-base columns under a stiff deck: Pcr = π²EI/h² (K = 1)', source: 'Timoshenko & Gere, Theory of Elastic Stability, §2.4 (sway column with both ends fixed against rotation)',
      inputs: { swFac: 0, gG: 1, gQ: 0, live: 0, Hx: 0, tPlate: 10, secSteel: 0, nx: 1, ny: 1, nLevels: 1, zLow: 20, cD: 300, ct: 10, bracing: 'none', columns: 'corners', nSub: 4, equip: [{ name: 'Block', W: 100, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
      expect: [{ kpi: 'alphaCr', value: 2.027, tol: 0.02 }, { kpi: 'Kcol', value: 1.0, tol: 0.02 }] },
    { // σcr = k π²E/(12(1 − ν²)) (t/b)² with k = 4 for a long simply supported plate: π²·210000/(12·0.91) = 189 800 MPa; (10/700)² → 38.73 MPa; × 4 = 154.9 MPa.
      name: 'Long simply supported plate in compression: σcr = 4 π²E/(12(1−ν²)) (t/b)²', source: 'Timoshenko & Gere, Theory of Elastic Stability, §9.2',
      inputs: { tPlate: 10, sStiff: 700, Lpanel: 3.5 }, expect: [{ kpi: 'sigCrPlate', value: 154.9, tol: 0.005 }] },
  ],
};

function linspaceLam() { const out = []; for (let i = 0; i <= 60; i++) out.push(round(0.05 * i, 4)); return out; }
