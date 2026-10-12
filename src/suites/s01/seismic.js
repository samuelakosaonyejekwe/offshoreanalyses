// Seismic analysis: response-spectrum analysis (ISO 19901-2 site spectrum, CQC) combined with gravity for the member
// checks, and a modal time-history analysis under a recorded or spectrum-compatible ground motion.
import { G, need, run, fmt, round, adv, check, linspace, steelCo2, CODE_CAVEAT, CONFIRM, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { geomInputs, foundationInputs, criteriaInputs, meshInputs, buildJacket, loadCase, solve, memberForces, checkMember, byGroup, pileForces, nodeValuesFromMembers, plot3d, groupChecks, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { isoSeismicSpectrum, responseSpectrumAnalysis, modalTimeHistory, spectrumCompatibleAccelerogram, responseSpectrum } from '../../engine/structural/framedynamics.js';

const DIRS = { x: ['x'], xy: ['x', 'y'], xyz: ['x', 'y', 'z'] };

export default {
  id: 's01-seismic', suite: 1, title: 'Seismic response — spectrum and time history',
  summary: 'Earthquake base shear, deck displacement, member utilisation and pile-head forces from a response-spectrum analysis, with displacement and base-shear histories from a modal time-history analysis.',
  covers: [10], fidelity: 'numerical',
  standards: ['ISO 19901-2:2017 §7–§9 (seismic design procedures, ELE / ALE)', 'ISO 19902:2020 §11', 'API RP 2EQ (2014)', 'EN 1998-1:2004 §4.3.3 (modal combination)'],
  models: ['Response spectrum analysis with CQC / SRSS modal combination', 'ISO 19901-2 site response spectrum', 'Modal superposition time history (exact piecewise-linear recurrence)', 'Spectrum-compatible synthetic accelerogram (SIMQKE procedure)', 'Eigenvalue analysis with added mass', 'ISO 19902 tubular member checks'],
  limits: `Linear elastic response on linearised pile-head springs with uniform base excitation: no soil nonlinearity, kinematic interaction, liquefaction or ductility check — the abnormal-level earthquake needs a nonlinear pushover or time-history analysis. Mapped accelerations and site coefficients are inputs to be taken from ISO 19901-2 or a site-specific hazard study. ${CODE_CAVEAT}`,
  searchTerms: 'offshore jacket seismic response spectrum analysis ISO 19901-2',
  inputs: [
    ...geomInputs(),
    { key: 'Wtop', label: 'Topsides weight', unit: 'MN', type: 'number', default: 40, min: 0, max: 2000, group: 'Mass and gravity', role: 'bc', bind: 'topsides.weight', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Acts as gravity load and, divided by g, as mass on the leg tops.' },
    { key: 'swCont', label: 'Self-weight contingency', unit: '-', type: 'number', default: 1.1, min: 0, max: 1.5, group: 'Mass and gravity', role: 'bc' },
    { key: 'addCa', label: 'Added mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 2, group: 'Mass and gravity', role: 'criteria' },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 300, group: 'Mass and gravity', role: 'environment' },
    { key: 'mgDepth', label: 'Marine growth down to depth', unit: 'm', type: 'number', default: 50, min: 0, max: 400, group: 'Mass and gravity', role: 'environment' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Mass and gravity', role: 'environment', bind: 'site.rhoW' },
    { key: 'Sa02', label: 'Mapped spectral acceleration at 0.2 s (rock)', unit: 'g', type: 'number', default: 0.5, min: 0.01, max: 3, group: 'Seismic action', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'ISO 19901-2 Annex B maps (1000-year return period) or a site-specific hazard study.' },
    { key: 'Sa10', label: 'Mapped spectral acceleration at 1.0 s (rock)', unit: 'g', type: 'number', default: 0.2, min: 0.005, max: 2, group: 'Seismic action', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'siteCa', label: 'Site coefficient Ca (short period)', unit: '-', type: 'number', default: 1.2, min: 0.8, max: 2.5, group: 'Seismic action', role: 'soil', help: `ISO 19901-2 Table 6 by site class and Sa(0.2); soil amplification. ${CONFIRM}` },
    { key: 'siteCv', label: 'Site coefficient Cv (long period)', unit: '-', type: 'number', default: 1.7, min: 0.8, max: 3.5, group: 'Seismic action', role: 'soil', help: `ISO 19901-2 Table 7. ${CONFIRM}` },
    { key: 'scale', label: 'Scale factor on the spectrum', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 3, group: 'Seismic action', role: 'criteria', help: 'ELE spectrum = ALE (1000-year) spectrum divided by the seismic reserve capacity factor Cr and adjusted by the scale factor NALE; 0.5 is a typical ELE/ALE ratio. 1.0 analyses the mapped spectrum.' },
    { key: 'zeta', label: 'Damping ratio', unit: '-', type: 'number', default: 0.05, min: 0.01, max: 0.2, group: 'Seismic action', role: 'criteria', help: '5 % of critical is the ISO 19901-2 reference.' },
    { key: 'dirs', label: 'Excitation directions', type: 'select', default: 'xyz', group: 'Seismic action', role: 'bc', options: [{ value: 'x', label: 'x only' }, { value: 'xy', label: 'x and y' }, { value: 'xyz', label: 'x, y and vertical (half spectrum)' }] },
    { key: 'comb', label: 'Modal combination', type: 'select', default: 'CQC', group: 'Seismic action', role: 'control', options: [{ value: 'CQC', label: 'CQC' }, { value: 'SRSS', label: 'SRSS' }] },
    { key: 'dirComb', label: 'Directional combination', type: 'select', default: 'SRSS', group: 'Seismic action', role: 'control', options: [{ value: 'SRSS', label: 'SRSS' }, { value: '100-30-30', label: '100 % + 30 % + 30 %' }] },
    { key: 'accel', label: 'Ground acceleration record in x (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Acceleration', yUnit: 'm/s²', group: 'Seismic action', role: 'bc', help: 'Leave empty to generate a record compatible with the spectrum above (fixed seed).' },
    { key: 'dur', label: 'Duration of the synthetic record', unit: 's', type: 'number', default: 20, min: 8, max: 60, group: 'Seismic action', role: 'control' },
    ...foundationInputs(), ...criteriaInputs({ gG: 1.1, gE: 0.9 }),
    { key: 'nModes', label: 'Modes', type: 'int', default: 12, min: 1, max: 40, group: 'Solver', role: 'control', help: 'Enough modes for about 90 % mass participation.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.02, min: 0.002, max: 0.1, group: 'Solver', role: 'control' },
    ...meshInputs({ nSub: 1, phase: false }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'baseShear' },
  meshTime: { key: 'dt', label: 'Time step', levels: [0.04, 0.02, 0.01], qoi: 'thDeckDisp', isSize: true },
  compute(inp, ctx) {
    // SI conversion: spectral accelerations g → m/s², MN → N and kg.
    const J = buildJacket(inp, ctx), mTop = inp.Wtop * 1e6 / G, nT = J.topNodes.length, dirs = DIRS[inp.dirs];
    const gi = { ...inp, FxTop: 0, MyTop: 0, windA: 0, Uw: 0, surge: 0 };
    const stat = solve(J, [loadCase('Gravity', J, gi, { fG: inp.gG, fE: 0 })]).Gravity;
    if (mTop > 0) J.model.nodeMasses = J.topNodes.map((k) => ({ node: k, m: mTop / nT }));
    const opt = { nModes: inp.nModes, waterline: 0, rhoW: inp.rhoW, Ca: inp.addCa };
    if (inp.mg > 0) opt.marineGrowth = { thickness: inp.mg / 1e3, rho: 1325, zTop: 0, zBottom: -inp.mgDepth };
    const mo = run(() => modalAnalysis(J.model, opt), 'Modal analysis');
    const sp = isoSeismicSpectrum({ Sa02: inp.Sa02 * G, Sa10: inp.Sa10 * G, Ca: inp.siteCa, Cv: inp.siteCv, zeta: inp.zeta, scale: inp.scale });
    const rsa = run(() => responseSpectrumAnalysis(mo, { spectrum: { x: sp.Sa, y: sp.Sa, z: (T) => 0.5 * sp.Sa(T) }, directions: dirs, zeta: inp.zeta, combination: inp.comb, directional: inp.dirComb }), 'Response-spectrum analysis');
    const u = rsa.displacements, deck = Math.max(...J.topNodes.map((k) => Math.hypot(u[6 * k], u[6 * k + 1]))), V = Math.hypot(rsa.baseShear.x, dirs.includes('y') ? rsa.baseShear.y : 0);
    // gravity ± factored seismic member forces (the spectrum result carries no sign: both senses of axial force are checked)
    const FG = memberForces(J, stat), rows = J.orig.map((m, i) => {
      let N = 0, My = 0, Mz = 0, Vs = 0, T = 0;
      for (const e of J.memberMap[i]) { const f = rsa.memberForces[e]; N = Math.max(N, f[0], f[6]); My = Math.max(My, f[4], f[10]); Mz = Math.max(Mz, f[5], f[11]); Vs = Math.max(Vs, Math.hypot(Math.max(f[1], f[7]), Math.max(f[2], f[8]))); T = Math.max(T, f[3], f[9]); }
      const g = FG[i], e = inp.gE, cases = [-1, 1].map((s) => ({ N: g.N + s * e * N, My: g.My + e * My, Mz: g.Mz + e * Mz, V: g.V + e * Vs, T: g.T + e * T }));
      const cs = cases.map((f) => checkMember(J, m, f, inp)), k = cs[0].uc >= cs[1].uc ? 0 : 1;
      return { i, group: m.group, label: m.label, L: m.L, N: cases[k].N, NE: N, M: Math.hypot(cases[k].My, cases[k].Mz), ME: Math.hypot(My, Mz), V: cases[k].V, T: cases[k].T, vm: cs[k].vm, uc: cs[k].uc, ref: cs[k].ref, case: 'G ± E' };
    });
    const groups = byGroup(rows), gov = rows.reduce((a, r) => (r.uc > a.uc ? r : a)), pfG = pileForces(J, stat);
    const legE = J.orig.filter((m) => J.mudNodes.includes(m.a) || J.mudNodes.includes(m.b)).map((m) => rows[m.i].NE), pileE = legE.length ? Math.max(...legE) : 0;
    // time history in x: recorded motion or a spectrum-compatible synthetic record on a fixed 0.01 s grid
    let rec;
    if (inp.accel?.x?.length > 3) {
      const tx = inp.accel.x.map(Number), ay = inp.accel.y.map(Number);
      need(tx.every((v, i) => Number.isFinite(v) && Number.isFinite(ay[i]) && (i === 0 || v > tx[i - 1])), 'The acceleration record needs increasing times and finite accelerations.');
      rec = { t: tx, acc: ay, source: 'recorded' };
    } else {
      const s = run(() => spectrumCompatibleAccelerogram({ spectrum: sp.Sa, duration: inp.dur, dt: 0.01, zeta: inp.zeta, fMin: 0.2, fMax: 25, nFreq: 60, iterations: 6, seed: 7 }), 'Synthetic accelerogram');
      rec = { t: Array.from(s.t), acc: Array.from(s.acc), source: 'synthetic', err: s.meanError };
    }
    const tEnd = rec.t[rec.t.length - 1], nSteps = Math.max(10, Math.round((tEnd + 2 * mo.periods[0]) / inp.dt));
    let ip = 0; const at = (t) => { if (t <= rec.t[0] || t >= tEnd) return 0; while (ip > 0 && rec.t[ip] > t) ip--; while (rec.t[ip + 1] < t) ip++; const w = (t - rec.t[ip]) / (rec.t[ip + 1] - rec.t[ip]); return rec.acc[ip] + w * (rec.acc[ip + 1] - rec.acc[ip]); };
    const th = run(() => modalTimeHistory(mo, { dt: inp.dt, nSteps, zeta: inp.zeta, baseAcc: { x: at }, outputs: [{ node: J.topNodes[0], dof: 'ux' }] }), 'Time-history analysis');
    const Vt = new Float64Array(nSteps + 1);
    mo.participation.forEach((p, k) => { const c = p.x * mo.omegas[k] ** 2, q = th.q[k]; for (let i = 0; i <= nSteps; i++) Vt[i] += c * q[i]; });
    const dk = th.outputs[0], thDeck = dk.reduce((s, v) => Math.max(s, Math.abs(v)), 0), thV = Vt.reduce((s, v) => Math.max(s, Math.abs(v)), 0), pga = rec.acc.reduce((s, v) => Math.max(s, Math.abs(v)), 0);
    const step = Math.max(1, Math.ceil((nSteps + 1) / 1200)), ti = [], di = [], vi = [];
    for (let i = 0; i <= nSteps; i += step) { ti.push(round(th.t[i])); di.push(round(dk[i] * 1e3)); vi.push(round(Vt[i] / 1e6)); }
    const aStep = Math.max(1, Math.ceil(rec.t.length / 1200)), Tp = linspace(0.05, 5, 60), rsRec = responseSpectrum({ acc: Float64Array.from({ length: Math.round(tEnd / 0.01) + 1 }, (_, i) => at(i * 0.01)), dt: 0.01, periods: Tp, zeta: inp.zeta });
    const steel = steelMass(J, inp), gg = groups.reduce((a, g) => (g.uc > a.uc ? g : a)), gSteel = groupSteel(J, gg.group), ok = gov.uc <= 1, mp = rsa.massParticipation.x;
    const pAx = pfG.axialMax + inp.gE * pileE, pTen = Math.max(0, inp.gE * pileE - Math.min(...pfG.piles.map((p) => p.axial)));
    return {
      kpis: [
        { key: 'baseShear', label: 'Seismic base shear (spectrum, unfactored)', unit: 'MN', value: V / 1e6, explain: 'Combined peak horizontal base shear from the response-spectrum analysis. Compare with the storm base shear: where the earthquake is smaller, seismic loading only governs local, mass-driven details.' },
        { key: 'deckDisp', label: 'Deck displacement (spectrum)', unit: 'mm', value: deck * 1e3, explain: 'Peak horizontal displacement of the leg tops relative to the ground. It sets the demand on risers, conductors and bridge bearings.' },
        { key: 'maxUc', label: 'Highest member utilisation, gravity ± earthquake', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'Code utilisation with factored gravity and the factored spectrum forces in the worse sense. Above 1.0 the member does not meet the extreme-level earthquake strength requirement.' },
        { key: 'T1', label: 'First natural period', unit: 's', value: mo.periods[0], explain: 'Fundamental period used to read the spectrum; long-period jackets sit on the falling 1/T branch.' },
        { key: 'SaT1', label: 'Spectral acceleration at T1', unit: 'g', value: sp.Sa(mo.periods[0]) / G, explain: 'Design spectral acceleration of the first mode, including site coefficients, damping and the scale factor.' },
        { key: 'massPart', label: 'Mass participation of the modes used (x)', unit: '%', value: mp * 100, limit: 90, sense: 'min', explain: 'Share of the mass captured by the modes in the analysis; below 90 % add modes, because missing mass means missing base shear.' },
        { key: 'thDeckDisp', label: 'Peak deck displacement (time history, x)', unit: 'mm', value: thDeck * 1e3, explain: 'Largest deck displacement in the time-history analysis for the x component. It should be close to the x-direction spectrum result when the record matches the spectrum.' },
        { key: 'thBaseShear', label: 'Peak base shear (time history, x)', unit: 'MN', value: thV / 1e6, explain: 'Largest instantaneous base shear in the time history; unlike the spectrum value it carries sign and timing.' },
        { key: 'pga', label: 'Peak ground acceleration of the record', unit: 'g', value: pga / G, explain: 'Largest acceleration of the ground motion used for the time history.' },
        { key: 'pileAxial', label: 'Largest pile-head compression, gravity + earthquake', unit: 'MN', value: pAx / 1e6, explain: 'Gravity reaction plus the factored seismic axial force of the leg at the mudline: the seismic demand on pile axial capacity.' },
      ],
      checks: [...groupChecks(groups, inp, 'gravity ± earthquake'),
        check('Mass participation of the modes', 0.9, Math.max(mp, 1e-6), '-', 'ISO 19901-2:2017 §9.2; EN 1998-1 §4.3.3.3.1 (90 % of the mass)', 'A response-spectrum analysis with too few modes misses part of the inertia force. Add modes until 90 % of the mass participates.')],
      plots: [
        { type: 'line', title: 'Design spectrum and spectrum of the ground motion', x: { label: 'Period', unit: 's' }, y: { label: 'Spectral acceleration', unit: 'g' }, series: [{ name: 'ISO 19901-2 site spectrum', x: Tp.map((v) => round(v)), y: Tp.map((T) => round(sp.Sa(T) / G)) }, { name: `Record (${rec.source})`, x: Tp.map((v) => round(v)), y: Array.from(rsRec.Sa, (v) => round(v / G)), dash: true }, { name: 'Modes', x: mo.periods.slice(0, 8).map((v) => round(v)), y: mo.periods.slice(0, 8).map((T) => round(sp.Sa(T) / G)) }],
          explain: 'Design spectrum with the natural periods marked on it, and the response spectrum of the ground motion used for the time history. The two curves should agree near the important periods.' },
        { type: 'line', title: 'Deck displacement history (x)', x: { label: 'Time', unit: 's' }, y: { label: 'Displacement', unit: 'mm' }, series: [{ name: 'Leg top', x: ti, y: di }], explain: 'Displacement of the leg top relative to the ground during the earthquake. Look at the number of large cycles and whether the response keeps growing, which indicates resonance with the record.' },
        { type: 'line', title: 'Base shear history (x)', x: { label: 'Time', unit: 's' }, y: { label: 'Base shear', unit: 'MN' }, series: [{ name: 'Base shear', x: ti, y: vi }], explain: 'Total horizontal force transferred to the foundation in x through the event: the foundation response history. Peak values compare with the spectrum result.' },
        { type: 'line', title: 'Ground acceleration', x: { label: 'Time', unit: 's' }, y: { label: 'Acceleration', unit: 'm/s²' }, series: [{ name: rec.source, x: rec.t.filter((v, i) => i % aStep === 0).map((v) => round(v)), y: rec.acc.filter((v, i) => i % aStep === 0).map((v) => round(v)) }], explain: 'Input motion at the mudline. A synthetic record has a build-up, a strong phase and a decay; a recorded motion should be scaled to the design spectrum before use.' },
        plot3d(J, { title: 'Utilisation under gravity ± earthquake and spectrum displacements', values: nodeValuesFromMembers(J, rows.map((r) => r.uc)), label: 'Utilisation', unit: '-', u, explain: 'Member utilisation at the nodes with the combined peak displacement shape of the spectrum analysis (all values positive, so the shape is an envelope, not a deflected state).' }),
      ],
      tables: [
        { title: 'Modes and spectral ordinates', columns: [{ key: 'k', label: 'Mode' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'Sa', label: 'Sa', unit: 'g' }, { key: 'mx', label: 'Mass x', unit: '%' }, { key: 'my', label: 'Mass y', unit: '%' }, { key: 'mz', label: 'Mass z', unit: '%' }], rows: mo.periods.map((T, k) => ({ k: k + 1, T: round(T), Sa: round(sp.Sa(T) / G), mx: round(mo.massFraction[k].x * 100, 3), my: round(mo.massFraction[k].y * 100, 3), mz: round(mo.massFraction[k].z * 100, 3) })), explain: 'Period, design spectral acceleration and effective mass of each mode: the modes with high mass on the plateau of the spectrum dominate.' },
        { title: 'Most utilised members, gravity ± earthquake', columns: [{ key: 'i', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'N', label: 'Axial combined', unit: 'kN' }, { key: 'NE', label: 'Axial seismic', unit: 'kN' }, { key: 'M', label: 'Moment combined', unit: 'kN·m' }, { key: 'ME', label: 'Moment seismic', unit: 'kN·m' }, { key: 'uc', label: 'Utilisation' }, { key: 'ref', label: 'Governing check' }],
          rows: rows.slice().sort((a, b) => b.uc - a.uc).slice(0, 15).map((r) => ({ i: r.i, g: r.label, N: round(r.N / 1e3), NE: round(r.NE / 1e3), M: round(r.M / 1e3), ME: round(r.ME / 1e3), uc: round(r.uc), ref: r.ref })), explain: 'Member forces from gravity combined with the factored spectrum forces, and the seismic share alone. Members where the seismic share dominates are the ones to detail for ductility.' },
        { title: 'Base shear by direction', columns: [{ key: 'd', label: 'Excitation' }, { key: 'vx', label: 'Base shear x', unit: 'MN' }, { key: 'vy', label: 'Base shear y', unit: 'MN' }, { key: 'vz', label: 'Vertical', unit: 'MN' }], rows: dirs.map((d) => ({ d, vx: round(rsa.perDirection[d].baseShear.x / 1e6), vy: round(rsa.perDirection[d].baseShear.y / 1e6), vz: round(rsa.perDirection[d].baseShear.z / 1e6) })), explain: 'Peak reactions for each excitation component before the directional combination.' },
      ],
      notes: [`T1 = ${fmt(mo.periods[0])} s, Sa(T1) = ${fmt(sp.Sa(mo.periods[0]) / G)} g; base shear ${fmt(V / 1e6)} MN (${fmt(V / (mo.totalMass.x * G) * 100)} % of the vibrating weight); ${fmt(mp * 100)} % mass participation with ${mo.periods.length} modes.`,
        `Time history (${rec.source} record, PGA ${fmt(pga / G)} g${rec.err != null ? `, mean spectral mismatch ${fmt(rec.err * 100)} %` : ''}): deck ${fmt(thDeck * 1e3)} mm, base shear ${fmt(thV / 1e6)} MN against ${fmt(rsa.perDirection.x.baseShear.x / 1e6)} MN from the x spectrum.`,
        `Governing member ${gov.i} (${gov.label.toLowerCase()}), ${gov.ref}. Pile-head tension ${pTen > 0 ? `${fmt(pTen / 1e6)} MN` : 'does not occur'}.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Extreme-level earthquake strength is met: utilisation ${fmt(gov.uc)} with base shear ${fmt(V / 1e6)} MN and deck displacement ${fmt(deck * 1e3)} mm.`, mp < 0.9 ? `Raise the number of modes: only ${fmt(mp * 100)} % of the mass participates.` : 'Confirm the mapped accelerations and site class with a site-specific hazard study, then demonstrate the abnormal-level earthquake by a pushover with a seismic load pattern (ductility), and check deck appurtenances for the floor accelerations.', 'ISO 19901-2:2017 §8 (ELE) and §9 (ALE); ISO 19902:2020 §11', `Seismic force is mass times acceleration: lighter topsides reduce demand directly, whereas heavier bracing adds stiffness, shortens the period and can raise the spectral acceleration. The jacket holds ${steelCo2(steel)}`)
          : adv('engineering', 'fail', `Earthquake strength is exceeded: utilisation ${fmt(gov.uc)} in member ${gov.i} (${gov.label.toLowerCase()}, ${gov.ref}).`, `Strengthen the ${gg.label.toLowerCase()} by about ${fmt((gg.uc / 0.9 - 1) * 100)} % in wall thickness, or reduce topsides mass; check whether a softer, more ductile bracing arrangement moves the period down the spectrum.`, BASIS_MEMBER, `The strengthening is roughly ${steelCo2(gSteel * Math.max(gg.uc / 0.9 - 1, 0.05))}`),
        adv('economic', 'info', `Seismic base shear is ${fmt(V / (mo.totalMass.x * G) * 100)} % of the vibrating weight for a spectrum scaled by ${fmt(inp.scale)}.`, 'A site-specific probabilistic hazard study usually pays for itself in moderate-seismicity regions: mapped values are conservative, and the spectrum scales steel in the lower bays and pile penetration directly.', `${BASIS_COST}; ISO 19901-2:2017 §6.4 (detailed seismic action procedure)`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Single-mode oscillator: 20 m cantilever tube carrying 5 MN (509 858 kg) at the tip.
    // Flexibility L³/(3EI) + L/(G·A/2) = 2.3306e-7 + 2.41e-9 m/N → k = 4.2468e6 N/m; Rayleigh generalised mass
    // M + 0.2357·mL = 509 858 + 7 612 = 517 470 kg → T1 = 2π√(517 470/4.2468e6) = 2.1933 s.
    // Spectrum (Ca = Cv = 1, scale 1, Sa(1.0) = 0.2 g): Sa = 0.2·9.80665/2.1933 = 0.8942 m/s²; effective modal mass
    // (M + 0.375·mL)²/(M + 0.2357·mL) = 521 968²/517 470 = 526 505 kg → V = 0.4708 MN. The round column has two equal
    // sway modes, so two modes are requested; CQC adds their fully correlated x contributions.
    { name: 'Single-mode cantilever with tip mass: V = M_eff·Sa(T1)', source: 'Chopra, Dynamics of Structures (4th ed.) §6.7 and §8.3 (generalised SDOF) with the ISO 19901-2 spectrum',
      inputs: { geom: testFrame('column', { H: 20 }), Wtop: 5, swCont: 0, addCa: 0, mg: 0, Sa02: 0.5, Sa10: 0.2, siteCa: 1, siteCv: 1, scale: 1, dirs: 'x', found: 'fixed', nModes: 2, nSub: 4 },
      expect: [{ kpi: 'T1', value: 2.1933, tol: 0.005 }, { kpi: 'baseShear', value: 0.4708, tol: 0.01 }] },
  ],
};
