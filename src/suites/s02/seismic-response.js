// Seismic response of the topsides by response-spectrum analysis: modes of the deck frame with its equipment
// mass, ISO 19901-2 or EN 1998-1 spectrum, CQC combination, base shear, floor accelerations for equipment
// holding-down, and member checks under gravity plus earthquake.
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { responseSpectrumAnalysis, isoSeismicSpectrum, ec8Spectrum } from '../../engine/structural/framedynamics.js';
import { deckInputs, solverInput, buildDeck, loadCase, solveDeck, memberForces, memberChecks, groupSummary, dynamicModel, framePlot, equipmentOf, weights, need, run, fmt, round, adv, check, fsCheck, steelCo2Text, CONFIRM, DECK_LIMITS, G } from './_common.js';

export default {
  id: 's02-seismic-response', suite: 2, title: 'Seismic response of the topsides',
  summary: 'Response-spectrum analysis of the deck frame: modal periods and mass participation, seismic base shear and overturning, storey displacements, floor accelerations for equipment supports and member utilisation under gravity plus earthquake.',
  covers: [11], fidelity: 'numerical',
  standards: ['ISO 19901-2:2017 §7–§8 (seismic design procedures, ELE / ALE)', 'ISO 19901-3:2014 §7.4 (seismic actions on topsides)', 'API RP 2EQ (2014) §6', 'EN 1998-1:2004 §3.2.2, §4.3.3.3 (response spectrum, modal combination)', 'ISO 19902:2020 §11'],
  models: ['Eigenvalue analysis', 'Modal superposition', 'Response-spectrum method with CQC combination', 'Mass–stiffness–damping formulation', 'Lagrange’s equations of motion (consistent mass)', 'Matrix stiffness method'],
  limits: `${DECK_LIMITS} The spectrum is applied at the jacket interface as if it were the ground: the filtering and amplification of the motion by the jacket and foundation is represented only by the entered amplification factor, so an integrated jacket–topsides analysis or a deck-level floor spectrum should confirm the result. Linear modal response; ductility is represented only by the scale factor on the spectrum. Nonlinear time-history analysis is outside this module.`,
  publishes: ['topsides.seismicShear', 'topsides.f1'],
  searchTerms: 'seismic response spectrum analysis offshore topsides ISO 19901-2 floor acceleration equipment',
  inputs: [
    ...deckInputs(),
    { key: 'code', label: 'Spectrum', type: 'select', default: 'iso', group: 'Seismic action', role: 'environment', options: [{ value: 'iso', label: 'ISO 19901-2 / API RP 2EQ site spectrum' }, { value: 'ec8', label: 'EN 1998-1 elastic spectrum (type 1)' }] },
    { key: 'Sa02', label: 'Mapped spectral acceleration at 0.2 s', unit: 'g', type: 'number', default: 0.12, min: 0.01, max: 4, group: 'Seismic action', role: 'environment', showIf: { key: 'code', in: ['iso'] }, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'ISO 19901-2 Annex B maps (1000-year rock values) or the site-specific hazard study.' },
    { key: 'Sa10', label: 'Mapped spectral acceleration at 1.0 s', unit: 'g', type: 'number', default: 0.04, min: 0.005, max: 3, group: 'Seismic action', role: 'environment', showIf: { key: 'code', in: ['iso'] } },
    { key: 'Ca', label: 'Site coefficient Ca', unit: '-', type: 'number', default: 1.2, min: 0.5, max: 3, group: 'Seismic action', role: 'soil', showIf: { key: 'code', in: ['iso'] }, help: `ISO 19901-2 Table 6 by site class and Sa(0.2). ${CONFIRM}` },
    { key: 'Cv', label: 'Site coefficient Cv', unit: '-', type: 'number', default: 1.6, min: 0.5, max: 4, group: 'Seismic action', role: 'soil', showIf: { key: 'code', in: ['iso'] }, help: `ISO 19901-2 Table 7. ${CONFIRM}` },
    { key: 'ag', label: 'Design ground acceleration (type A ground)', unit: 'm/s²', type: 'number', default: 1.5, min: 0.05, max: 10, group: 'Seismic action', role: 'environment', bind: 'site.pga475', showIf: { key: 'code', in: ['ec8'] } },
    { key: 'ground', label: 'Ground type', type: 'select', default: 'C', group: 'Seismic action', role: 'soil', showIf: { key: 'code', in: ['ec8'] }, options: ['A', 'B', 'C', 'D', 'E'].map((v) => ({ value: v, label: `Type ${v}` })) },
    { key: 'scale', label: 'Scale factor on the spectrum', unit: '-', type: 'number', default: 1, min: 0.05, max: 5, group: 'Seismic action', role: 'bc', help: 'Combine here the return-period scaling (ALE to ELE: 1/Cr), the importance factor and the amplification of the motion between the seabed and the deck by the jacket (often 1.5–3 near the jacket period).' },
    { key: 'dirs', label: 'Excitation directions', type: 'select', default: 'xy', group: 'Seismic action', role: 'bc', options: [{ value: 'xy', label: 'x and y (SRSS)' }, { value: 'x', label: 'x only' }, { value: 'y', label: 'y only' }, { value: 'xyz', label: 'x, y and vertical' }] },
    { key: 'zeta', label: 'Damping ratio', unit: '%', type: 'number', default: 5, min: 0.5, max: 15, group: 'Dynamic properties', role: 'ic', help: 'ISO 19901-2 spectra are for 5 %; welded steel topsides 2–5 %.' },
    { key: 'psi', label: 'Fraction of live load present (mass and weight)', unit: '-', type: 'number', default: 0.75, min: 0, max: 1, group: 'Dynamic properties', role: 'ic', help: `ISO 19902 uses 75 % of the maximum variable actions in the seismic combination. ${CONFIRM}` },
    { key: 'nModes', label: 'Modes retained', type: 'int', default: 16, min: 1, max: 60, group: 'Dynamic properties', role: 'control' },
    { key: 'massReq', label: 'Required participating mass', unit: '%', type: 'number', default: 90, min: 50, max: 99, group: 'Acceptance criteria', role: 'criteria', help: 'EN 1998-1 §4.3.3.3.1: the modes retained should mobilise at least 90 % of the mass in each horizontal direction.' },
    { key: 'driftRatio', label: 'Storey drift limit (height / value)', unit: '-', type: 'number', default: 100, min: 50, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: `h/100 for bare steel framing, h/200 where brittle cladding, firewalls or rigid piping cross the storey. ${CONFIRM}` },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'baseShear' },
  compute(inp, ctx) {
    // SI: g → m/s², % → fraction.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), zeta = inp.zeta / 100;
    const sp = run(() => (inp.code === 'iso' ? isoSeismicSpectrum({ Sa02: inp.Sa02 * G, Sa10: inp.Sa10 * G, Ca: inp.Ca, Cv: inp.Cv, zeta, scale: inp.scale }) : (() => { const s = ec8Spectrum({ ag: inp.ag, groundType: inp.ground, type: 1, zeta }); return { Sa: (T) => inp.scale * s.Sa(T) }; })()), 'Response spectrum');
    const Sa = sp.Sa, directions = inp.dirs.split('');
    const ref = dynamicModel(deck, inp, { liveFrac: inp.psi, nSub: inp.nSub, eqs });
    const modal = run(() => modalAnalysis(ref.model, { nModes: Math.round(inp.nModes), tol: 1e-6 }), 'Modal analysis (check supports and bracing)');
    const rsa = run(() => responseSpectrumAnalysis(modal, { spectrum: Sa, directions, zeta, combination: 'CQC', directional: 'SRSS' }), 'Response-spectrum analysis');
    const nN = ref.model.nodes.length, u = rsa.displacements, T = Array.from(modal.periods);
    need(T.every((v) => Number.isFinite(v) && v > 0), 'The modal analysis returned no valid period: the frame is a mechanism or has no mass.');
    const hDirs = directions.filter((d) => d !== 'z'), Vb = Math.hypot(...hDirs.map((d) => rsa.baseShear[d]));
    const w = weights(deck, inp, { liveFrac: inp.psi, eqs }), massPart = Math.min(...hDirs.map((d) => rsa.massParticipation[d]));
    // Floor acceleration per node and horizontal direction: SRSS of modal accelerations φ Γ Sa plus the rigid
    // (missing-mass) part at the zero-period acceleration.
    const dIdx = { x: 0, y: 1, z: 2 }, acc = new Array(nN).fill(0);
    for (const d of hDirs) {
      const miss = Math.max(1 - rsa.massParticipation[d], 0) * Sa(0);
      for (let i = 0; i < nN; i++) { let s = miss * miss; for (let n = 0; n < T.length; n++) s += (modal.modes[n][6 * i + dIdx[d]] * modal.participation[n][d] * Sa(T[n])) ** 2; acc[i] = Math.hypot(acc[i], Math.sqrt(s)); }
    }
    const levels = deck.levelsZ.map((z, k) => { const ns = deck.levelNodes[k]; return { level: k + 1, z, a: Math.max(...ns.map((n) => acc[n])), d: Math.max(...ns.map((n) => Math.hypot(u[6 * n], u[6 * n + 1]))) }; });
    const storeys = levels.map((l, k) => { const h = k === 0 ? l.z : l.z - levels[k - 1].z, dr = Math.abs(l.d - (k === 0 ? 0 : levels[k - 1].d)); return { ...l, h, drift: dr }; });
    const wd = storeys.reduce((a, c) => (c.drift / c.h > a.drift / a.h ? c : a));
    // Gravity + earthquake member forces: static gravity (G + ψQ) with the spectral forces added in the unfavourable sense.
    const { ref: sref, sol } = solveDeck(deck, [loadCase(deck, inp, { gLive: inp.psi, eqs })], inp.nSub);
    const fg = memberForces(deck, sref, sol.results.design);
    const comb = fg.map((f) => {
      let N = 0, Vy = 0, Vz = 0, Tq = 0, My = 0, Mz = 0;
      for (const k of ref.memberMap[f.i]) { const q = rsa.memberForces[k]; for (const o of [0, 6]) { N = Math.max(N, q[o]); Vy = Math.max(Vy, q[o + 1]); Vz = Math.max(Vz, q[o + 2]); Tq = Math.max(Tq, q[o + 3]); My = Math.max(My, q[o + 4]); Mz = Math.max(Mz, q[o + 5]); } }
      return { ...f, N: f.N - N, Vy: f.Vy + Vy, Vz: f.Vz + Vz, T: f.T + Tq, My: f.My + My, Mz: f.Mz + Mz, Ns: N };
    });
    const rows = memberChecks(deck, comb), rowsG = memberChecks(deck, fg), gov = rows.reduce((a, c) => (c.uc > a.uc ? c : a)), groups = groupSummary(rows), ucG = Math.max(...rowsG.map((r) => r.uc));
    const eqRows = eqs.map((e) => { const a = levels[e.k].a, H = e.mass * a, M = H * e.zc, arm = Math.max(Math.min(e.a, e.b), 0.5), up = M / arm - e.W / 2; return { name: e.name, level: e.k + 1, a, H, M, up }; });
    const worstEq = eqRows.reduce((a, c) => (c.a * 1 > a.a ? c : a), eqRows[0] || { a: 0, name: '-', H: 0, up: 0 });
    const checks = [
      check('Member strength under gravity plus earthquake', gov.uc, 1, 'interaction', gov.ref, `Highest member utilisation (${gov.group.toLowerCase()} ${gov.i + 1}) with the CQC seismic forces added to gravity in the unfavourable sense. Above 1.0 the member yields or buckles in the design earthquake; for the abnormal-level event, ductile yielding may be acceptable if a nonlinear check shows no collapse.`),
      fsCheck('Participating mass of the retained modes', massPart, inp.massReq / 100, 'EN 1998-1 §4.3.3.3.1(3); ISO 19901-2 §8.3', 'Share of the mass mobilised by the retained modes in the weakest horizontal direction. Below the requirement the base shear and member forces are underestimated: retain more modes.'),
      check(`Storey drift — level ${wd.level}`, wd.drift * 1e3, wd.h / inp.driftRatio * 1e3, 'mm', `Storey height / ${inp.driftRatio} (project criterion; EN 1998-1 §4.4.3.2 gives h/200 for brittle non-structural elements)`, 'Relative seismic displacement between consecutive decks, which racks piping, firewalls and cladding and adds P–Δ moments.'),
    ];
    const strength = checks[0].uc, ok = checks.every((c) => c.uc <= 1), govC = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    const Tp = []; for (let i = 0; i <= 80; i++) Tp.push(round(0.02 * 1.07 ** i * (i === 0 ? 0 : 1), 4));
    const colT = deck.orig.filter((o) => !o.horizontal).reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'T1', label: 'Fundamental period', unit: 's', value: T[0], explain: 'Longest natural period of the topsides on a rigid interface. On a real jacket the combined system is longer; compare with the jacket period to judge interaction.' },
        { key: 'Sa1', label: 'Spectral acceleration at the fundamental period', unit: 'g', value: Sa(T[0]) / G, explain: 'Ordinate of the scaled design spectrum at the first period: the acceleration the first mode responds to.' },
        { key: 'baseShear', label: 'Seismic base shear', unit: 'kN', value: Vb / 1e3, explain: 'Horizontal force at the jacket interface from the combined modes and directions (CQC, SRSS).' },
        { key: 'shearCoef', label: 'Base shear over weight', unit: '-', value: Vb / w.total, explain: 'Seismic coefficient: base shear divided by the weight present in the earthquake. Typical values are 0.05–0.4 depending on seismicity and period.' },
        { key: 'massPart', label: 'Participating mass (weakest direction)', unit: '%', value: massPart * 100, limit: inp.massReq, sense: 'min', explain: 'Cumulative effective modal mass of the retained modes over the total mass.' },
        { key: 'topDisp', label: 'Seismic displacement of the top deck', unit: 'mm', value: levels[levels.length - 1].d * 1e3, explain: 'Peak horizontal displacement of the top level relative to the interface.' },
        { key: 'topAcc', label: 'Peak acceleration of the top deck', unit: 'g', value: levels[levels.length - 1].a / G, explain: 'Peak horizontal acceleration at the top level: the input for equipment holding-down and piping supports there.' },
        { key: 'maxUc', label: 'Highest member utilisation', unit: '-', value: round(gov.uc, 5), limit: 1, sense: 'max', explain: `Gravity plus earthquake (gravity alone gives ${fmt(ucG)}).` },
        { key: 'eqShear', label: 'Largest equipment holding-down shear', unit: 'kN', value: Math.max(0, ...eqRows.map((e) => e.H)) / 1e3, explain: 'Mass times floor acceleration of the heaviest-loaded equipment item: the horizontal force its anchor bolts and stops must carry.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Design spectrum and modal periods', x: { label: 'Period', unit: 's' }, y: { label: 'Spectral acceleration', unit: 'g' }, series: [{ name: 'Design spectrum', x: Tp, y: Tp.map((t) => round(Sa(t) / G)) }, { name: 'Modes', x: T.map((t) => round(t)), y: T.map((t) => round(Sa(t) / G)), mode: 'markers' }],
          explain: 'Scaled design spectrum with the natural periods of the topsides marked. Modes on the plateau attract the largest accelerations; long-period sway modes sit on the falling branch.' },
        framePlot({ title: 'Peak seismic displacement', ref, u, values: Array.from({ length: nN }, (_, i) => Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2]) * 1e3), valueLabel: 'Peak displacement', valueUnit: 'mm', explain: 'Envelope of peak displacements from the modal combination (all values positive, so the shape is an envelope, not a deflected shape at one instant).' }),
        { type: 'line', title: 'Floor acceleration and displacement over height', x: { label: 'Peak acceleration (g) / displacement (dm)', unit: '-' }, y: { label: 'Elevation above the interface', unit: 'm' }, series: [{ name: 'Acceleration (g)', x: [round(Sa(0) / G), ...levels.map((l) => round(l.a / G))], y: [0, ...levels.map((l) => l.z)] }, { name: 'Displacement (dm)', x: [0, ...levels.map((l) => round(l.d * 10))], y: [0, ...levels.map((l) => l.z)] }],
          explain: 'Peak horizontal acceleration and displacement of each deck. Acceleration growing with height is the whip effect equipment on the upper decks must be anchored for.' },
        { type: 'bar', title: 'Utilisation by member group', categories: groups.map((g) => g.group), series: [{ name: 'Gravity + earthquake', y: groups.map((g) => round(g.uc)) }, { name: 'Gravity', y: groups.map((g) => round(Math.max(...rowsG.filter((r) => r.group === g.group).map((r) => r.uc)))) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Governing utilisation of each member group with and without the earthquake: braces and legs normally take the seismic increment.' },
      ],
      tables: [
        { title: 'Modes', columns: [{ key: 'mode', label: 'Mode' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'Sa', label: 'Sa', unit: 'g' }, { key: 'mx', label: 'Mass x', unit: '%' }, { key: 'my', label: 'Mass y', unit: '%' }, { key: 'mz', label: 'Mass z', unit: '%' }],
          rows: T.map((t, n) => ({ mode: n + 1, T: round(t), f: round(1 / t), Sa: round(Sa(t) / G), mx: round(modal.massFraction[n].x * 100, 3), my: round(modal.massFraction[n].y * 100, 3), mz: round(modal.massFraction[n].z * 100, 3) })), explain: 'Period, spectral ordinate and effective mass fraction of each retained mode. The modes with large mass fractions carry the base shear.' },
        { title: 'Equipment seismic loads', columns: [{ key: 'name', label: 'Item' }, { key: 'level', label: 'Level' }, { key: 'a', label: 'Floor acceleration', unit: 'g' }, { key: 'H', label: 'Horizontal force', unit: 'kN' }, { key: 'M', label: 'Overturning moment', unit: 'kN·m' }, { key: 'up', label: 'Net uplift per side', unit: 'kN' }],
          rows: eqRows.map((e) => ({ name: e.name, level: e.level, a: round(e.a / G), H: round(e.H / 1e3), M: round(e.M / 1e3), up: round(Math.max(e.up, 0) / 1e3) })), explain: 'Holding-down design forces for each equipment item treated as rigid: mass × floor acceleration at its deck, the moment about its base and the net uplift on one side after deducting half the weight (zero means no uplift).' },
      ],
      notes: [`T1 = ${fmt(T[0])} s; base shear ${fmt(Vb / 1e3)} kN = ${fmt(Vb / w.total * 100)} % of the ${fmt(w.total / G / 1e3)} t present.`, `${gov.group} ${gov.i + 1} governs at ${fmt(gov.uc)} (${gov.ref}); retained modes mobilise ${fmt(massPart * 100)} % of the mass.`,
        'The spectrum acts at the jacket interface: confirm the scale factor against the jacket analysis or use a deck-level floor spectrum.'],
      advice: [
        ok ? adv('engineering', govC.uc > 0.9 ? 'warn' : 'ok', `The frame carries the design earthquake at utilisation ${fmt(strength)}; base shear coefficient ${fmt(Vb / w.total)}, top-deck acceleration ${fmt(levels[levels.length - 1].a / G)} g.`,
          `Design equipment holding-down for the tabulated floor accelerations (largest ${fmt(worstEq.a / G)} g on ${worstEq.name}) and provide ductile detailing of brace connections; confirm the result in the integrated jacket–topsides model.`, 'ISO 19901-2:2017 §8; ISO 19901-3:2014 §7.4; EN 1998-1:2004 §4.3.3.3',
          `Seismic demand scales with mass: each tonne removed from the upper decks lowers the base shear by about ${fmt(Vb / (w.total / G / 1e3) / 1e3)} kN, and the legs, columns and braces resisting it weigh ${fmt(colT)} t (${steelCo2Text(colT)}).`)
          : adv('engineering', 'fail', `${govC.name} is not met (utilisation ${fmt(govC.uc)}): base shear ${fmt(Vb / 1e3)} kN, member utilisation ${fmt(strength)}, participating mass ${fmt(massPart * 100)} %.`,
            massPart < inp.massReq / 100 && strength <= 1 ? `Retain more modes (currently ${T.length}) until at least ${inp.massReq} % of the mass participates.` : govC.name.startsWith('Storey') ? 'Add bracing between the levels concerned or stiffen the columns; then re-check drift.' : `Strengthen the ${gov.group.toLowerCase()} group by about ${fmt((strength - 1) * 100, 2)} %, add braced bays in both directions, or move heavy equipment to the lower deck to cut the overturning.`,
            'ISO 19901-2:2017 §8; ISO 19902:2020 §11; EN 1998-1:2004 §4.3.3.3', `Strengthening the lateral system by ${fmt(Math.max(strength - 1, 0.05) * 100, 2)} % is about ${steelCo2Text(Math.max(strength - 1, 0.05) * colT)} Lowering mass and centre of gravity reduces the demand without adding steel.`),
        adv('economic', 'info', `Seismic design adds little steel when it is considered from the start (lateral system ${fmt(colT)} t here), but equipment anchorage and piping flexibility for ${fmt(levels[levels.length - 1].a / G)} g affect every package on the upper deck.`,
          'Issue the floor accelerations to equipment and piping vendors with the purchase specification; retrofitting holding-down bolts and pipe guides offshore costs far more than specifying them.', 'ISO 19901-2:2017 two-level (ELE / ALE) design philosophy; whole-life costing per ISO 15686-5'),
      ],
      publish: { 'topsides.seismicShear': Vb, 'topsides.f1': 1 / T[0] },
    };
  },
  verify: [{
    // Stiff deck on four fixed-base legs CHS 300 × 10, h = 20 m (see the vibration module): k = 48EI/h³ = 120 829 N/m,
    // generalised mass M + 13/35 m_legs = 148 492 + 2 125 = 150 617 kg → T = 2π√(M*/k) = 7.015 s.
    // ISO 19901-2 long-period branch (T > 4 s), Cv = 1, Sa(1.0) = 0.2 g: Sa = 4·0.2·9.80665/7.015² = 0.1594 m/s².
    // Effective modal mass Γ²/M* with Γ = M + ½ m_legs = 151 353 kg: 152 093 kg → V = 152 093·0.1594 = 24.25 kN.
    name: 'Single-mode base shear of a stiff deck on four legs: V = M_eff · Sa(T1)', source: 'Chopra, Dynamics of Structures, §13.7; ISO 19901-2 §7.1 spectrum shape (hand calculation)',
    inputs: { code: 'iso', Sa02: 0.5, Sa10: 0.2, Ca: 1, Cv: 1, scale: 1, dirs: 'x', psi: 0, nModes: 3, live: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, zLow: 20, cD: 300, ct: 10, bracing: 'none', columns: 'corners', nSub: 4, equip: [{ name: 'Block', W: 100, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'T1', value: 7.015, tol: 0.01 }, { kpi: 'baseShear', value: 24.25, tol: 0.03 }],
  }],
};
