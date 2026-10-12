// Seismic response of a piled wharf bent and liquefaction of its foundation soils: modal analysis, response
// spectrum analysis with a time-history cross-check, pile checks in the seismic combination, liquefaction
// triggering from SPT data with post-liquefaction settlement and loss of pile support, and the sliding-block
// displacement of the under-wharf slope.
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { ec8Spectrum, responseSpectrumAnalysis, spectrumCompatibleAccelerogram, modalTimeHistory } from '../../engine/structural/framedynamics.js';
import { isoTubularCheck } from '../../engine/structural/tubular.js';
import { sptProfile, liquefactionPotentialIndex, postLiquefactionSettlement } from '../../engine/geotech/liquefaction.js';
import { newmarkDisplacement } from '../../engine/geotech/slope.js';
import { G, need, run, r, f, cap, check, sfCheck, ucCheck, eng, eco, linspace, thin, bentInputs, pileInputs, bentModel, solveBent, subdivideBent, framePlot, deformScale, steelCo2, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

const SPT_DEFAULT = [
  { z: 2, N160: 24, FC: 8, gamma: 19, sus: 1 }, { z: 4, N160: 25, FC: 8, gamma: 19, sus: 1 }, { z: 6, N160: 27, FC: 10, gamma: 19.5, sus: 1 },
  { z: 9, N160: 8, FC: 60, gamma: 17.5, sus: 0 }, { z: 12, N160: 29, FC: 10, gamma: 19.5, sus: 1 }, { z: 16, N160: 32, FC: 5, gamma: 20, sus: 1 }, { z: 20, N160: 40, FC: 5, gamma: 20, sus: 1 },
];

/** Modal and response-spectrum solution of the bent with its seismic mass; extraFree lengthens the piles (lost soil support). */
function seismicBent(inp, ctx, spec, extraFree) {
  const B = bentModel({ ...inp, crane: [], Ffender: 0, Fmoor: 0, Fenv: 0 }, ctx, { scour: extraFree });
  // seismic mass carried by the bent: slab, surfacing, the ψE share of cargo and the crane share, lumped at the pile heads
  const mDeck = ((inp.slabT * 25e3 + inp.gSdl * 1e3 + inp.psiE * inp.qLive * 1e3) * B.sy * Math.max(B.deckLen, 0)) / G + inp.mCrane * 1e3;
  const heads = B.meta.headNodes, sub = subdivideBent(B.model, B.meta, inp.nSub);
  const model = { ...sub.model, nodeMasses: heads.map((n) => ({ node: n, m: mDeck / heads.length })) };
  const nFree = (model.nodes.length - B.meta.toeNodes.length) * 2;
  const modal = run(() => modalAnalysis(model, { nModes: Math.max(1, Math.min(Math.round(inp.nModes), nFree)), waterline: 0, Ca: 1 }), 'Modal analysis');
  const rsa = run(() => responseSpectrumAnalysis(modal, { spectrum: spec.design.Sa, directions: ['x'], zeta: inp.zeta / 100, combination: 'CQC' }), 'Response-spectrum analysis');
  const rsaE = inp.q > 1 ? responseSpectrumAnalysis(modal, { spectrum: spec.elastic.Sa, directions: ['x'], zeta: inp.zeta / 100, combination: 'CQC' }) : rsa;
  return { B, model, modal, rsa, rsaE, mDeck, sub, u: Math.abs(rsaE.displacements[6 * B.headSea]) };
}

export default {
  id: 's05-seismic-liquefaction', suite: 5, title: 'Seismic response and liquefaction',
  summary: 'Natural periods, spectral base shear, deck displacement and pile utilisation of the wharf bent in the design earthquake, liquefaction triggering and settlement of the foundation soils, the response with liquefied layers, and the residual slope displacement.',
  covers: [8], fidelity: 'numerical',
  standards: ['EN 1998-1:2004 §3.2.2, §4.3.3.3', 'EN 1998-5:2004 §4.1.4, Annex B (liquefaction)', 'ASCE 61-14 Seismic Design of Piers and Wharves', 'PIANC WG34 Seismic Design Guidelines for Port Structures (2001)', 'Youd et al. (2001); Idriss & Boulanger (2008) SPT-based triggering', 'ISO 19902:2020 §13 (tubular members)'],
  models: ['Modal analysis (subspace iteration, consistent mass with added water mass)', 'Response spectrum method (CQC)', 'Time-history analysis by modal superposition', 'Spectrum-compatible synthetic accelerogram', 'Liquefaction triggering correlations (simplified procedure, SPT)', 'Post-liquefaction volumetric and shear strain', 'Newmark sliding block'],
  limits: `Linear response-spectrum analysis of one bent in its own plane on equivalent-fixity piles; the behaviour factor stands in for ductility. Liquefaction is the simplified stress-based procedure with SPT correlations, and its effect on the piles is taken as loss of lateral support over the liquefied thickness: kinematic loading from lateral spreading, pore-pressure build-up in time and nonlinear soil–pile response need effective-stress site response and pushover analysis. ${NOT_HIFI}`,
  publishes: ['wharf.T1', 'wharf.seismicBaseShear', 'wharf.seismicDisplacement', 'wharf.liquefiedThickness'],
  searchTerms: 'pile supported wharf seismic design liquefaction lateral spreading ASCE 61',
  inputs: [
    ...bentInputs(), ...pileInputs(),
    { key: 'gSdl', label: 'Superimposed dead load', unit: 'kPa', type: 'number', default: 2.5, min: 0, max: 30, group: 'Seismic mass', role: 'ic' },
    { key: 'qLive', label: 'Uniform cargo load on the deck', unit: 'kPa', type: 'number', default: 30, min: 0, max: 200, group: 'Seismic mass', role: 'ic' },
    { key: 'psiE', label: 'Share of the cargo load present in the earthquake', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Seismic mass', role: 'criteria', help: `EN 1998-1 §4.2.4 with ψ2 = 0.8 for storage and φ = 0.5–1.0; ASCE 61 uses 10 % of the uniform live load. ${CONFIRM}` },
    { key: 'mCrane', label: 'Crane mass carried by this bent', unit: 't', type: 'number', default: 60, min: 0, max: 2000, group: 'Seismic mass', role: 'ic', help: 'Share of the crane mass that moves with the deck (a crane on flexible legs is partly isolated; ASCE 61 §3.5).' },
    { key: 'ag', label: 'Design ground acceleration on rock', unit: 'm/s²', type: 'number', default: 1.5, min: 0.1, max: 8, group: 'Earthquake', role: 'bc', bind: 'site.pga475', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Reference peak ground acceleration times the importance factor (475-year event for ordinary berths).' },
    { key: 'ground', label: 'Ground type', type: 'select', default: 'C', group: 'Earthquake', role: 'bc', options: ['A', 'B', 'C', 'D', 'E'].map((v) => ({ value: v, label: `Type ${v}` })), help: 'EN 1998-1 Table 3.1: A rock, B very dense, C dense or stiff, D loose to medium or soft to firm, E shallow soft layer on rock.' },
    { key: 'specType', label: 'Spectrum type', type: 'select', default: '1', group: 'Earthquake', role: 'bc', options: [{ value: '1', label: 'Type 1 (surface-wave magnitude above 5.5)' }, { value: '2', label: 'Type 2 (magnitude up to 5.5)' }] },
    { key: 'Mw', label: 'Moment magnitude', unit: '-', type: 'number', default: 7.0, min: 5, max: 9, group: 'Earthquake', role: 'bc', help: 'Used for the magnitude scaling factor of the liquefaction check.' },
    { key: 'q', label: 'Behaviour factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Earthquake', role: 'criteria', help: `1.0 elastic; 1.5 limited ductility; higher only with ductile pile-head detailing verified by pushover (ASCE 61 uses displacement-based design instead). ${CONFIRM}` },
    { key: 'zeta', label: 'Damping ratio', unit: '%', type: 'number', default: 5, min: 1, max: 20, group: 'Earthquake', role: 'criteria' },
    { key: 'spt', label: 'SPT profile below the seabed', type: 'table', group: 'Liquefaction', role: 'soil',
      columns: [{ key: 'z', label: 'Depth of the layer base', unit: 'm', type: 'number' }, { key: 'N160', label: '(N1)60', unit: '-', type: 'number' }, { key: 'FC', label: 'Fines content', unit: '%', type: 'number' }, { key: 'gamma', label: 'Total unit weight', unit: 'kN/m³', type: 'number' }, { key: 'sus', label: 'Susceptible (1 / 0)', unit: '-', type: 'number' }],
      default: SPT_DEFAULT, help: 'One row per layer from the borehole, depths increasing; each row represents the soil between the previous depth and its own. Mark clays and plastic silts as not susceptible (0).' },
    { key: 'fsLiq', label: 'Required factor of safety against liquefaction', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Liquefaction', role: 'criteria', help: `EN 1998-5 §4.1.4(11) asks for 1.25 on the stress ratio; 1.1–1.3 in port practice. ${CONFIRM}` },
    { key: 'ky', label: 'Yield acceleration of the under-wharf slope', unit: 'g', type: 'number', default: 0.15, min: 0.01, max: 1, group: 'Liquefaction', role: 'soil', help: 'Horizontal acceleration at which the slope or dyke starts to slide (pseudo-static factor of safety of 1), with residual strengths where layers liquefy.' },
    { key: 'uLimit', label: 'Allowable deck displacement', unit: 'mm', type: 'number', default: 250, min: 10, max: 2000, group: 'Acceptance', role: 'criteria', help: 'Set by crane rail gauge, joints and pile-head ductility for the performance level chosen.' },
    { key: 'dLimit', label: 'Allowable residual ground displacement', unit: 'mm', type: 'number', default: 300, min: 10, max: 3000, group: 'Acceptance', role: 'criteria', help: 'PIANC WG34 damage criteria: about 100 mm serviceable, 300 mm repairable for pile-supported wharves.' },
    { key: 'nModes', label: 'Modes', unit: '-', type: 'int', default: 6, min: 1, max: 20, group: 'Numerics', role: 'control' },
    { key: 'nSub', label: 'Elements per pile', unit: '-', type: 'int', default: 4, min: 1, max: 32, group: 'Numerics', role: 'control', help: 'Subdivision of each pile for the mass distribution in the dynamic analysis.' },
    { key: 'record', label: 'Recorded ground acceleration history (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Ground acceleration', yUnit: 'm/s²', group: 'Earthquake', role: 'bc', help: 'A site-specific or code-matched accelerogram at uniform time step. When given it replaces the synthetic spectrum-compatible record in the time-history check and in the sliding-block analysis.' },
    { key: 'duration', label: 'Duration of the synthetic record', unit: 's', type: 'number', default: 20, min: 5, max: 60, group: 'Numerics', role: 'control' },
    { key: 'dtTH', label: 'Time step of the time-history analysis', unit: 's', type: 'number', default: 0.01, min: 0.01, max: 0.08, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nSub', label: 'Elements per pile', levels: [1, 2, 4, 8], qoi: 'T1' },
  meshTime: { key: 'dtTH', label: 'Time step', levels: [0.04, 0.02, 0.01], qoi: 'uTH', isSize: true },
  compute(inp, ctx) {
    const zeta = inp.zeta / 100;
    const elastic = run(() => ec8Spectrum({ ag: inp.ag, groundType: inp.ground, type: Number(inp.specType), zeta }), 'Response spectrum');
    const design = inp.q > 1 ? ec8Spectrum({ ag: inp.ag, groundType: inp.ground, type: Number(inp.specType), zeta, q: inp.q }) : elastic;
    const spec = { elastic, design }, amax = inp.ag * elastic.S;
    // liquefaction triggering (simplified procedure) with the surface acceleration ag·S
    const rows = (inp.spt || []).map((q, i) => ({ z: Number(q.z), N160: Number(q.N160), FC: Number(q.FC) || 0, gamma: Number(q.gamma) * 1e3, sus: Number(q.sus) !== 0, i }));
    need(rows.length >= 1, 'Give at least one SPT layer.');
    rows.forEach((q, i) => { need(q.z > (i ? rows[i - 1].z : 0), `SPT row ${i + 1}: depths must increase.`); need(q.N160 > 0 && q.gamma > 10e3, `SPT row ${i + 1}: (N1)60 must be positive and the total unit weight above 10 kN/m³.`); });
    const liq = run(() => sptProfile({ points: rows.map((q) => ({ z: q.z, N160: q.sus ? q.N160 : 60, FC: q.FC, gamma: q.gamma })), waterTableDepth: 0, amax, Mw: inp.Mw }), 'Liquefaction assessment');
    const lay = rows.map((q, i) => ({ ...q, thick: q.z - (i ? rows[i - 1].z : 0), FS: q.sus ? Math.min(liq[i].FS, 5) : 5, CSR: liq[i].CSR, CRR: liq[i].CRR75, N160cs: liq[i].N160cs }));
    const sus = lay.filter((q) => q.sus);
    const crit = sus.length ? sus.reduce((m, q) => (q.FS < m.FS ? q : m), sus[0]) : null;
    const liqd = sus.filter((q) => q.FS < 1), tLiq = liqd.reduce((s, q) => s + q.thick, 0);
    const lpi = sus.length ? liquefactionPotentialIndex(sus.map((q) => q.z - q.thick / 2), sus.map((q) => q.FS)) : 0;
    const post = postLiquefactionSettlement(sus.map((q) => ({ thickness: q.thick, N160cs: q.N160cs, FS: q.FS })));
    // structure: as built, and with lateral support lost over the liquefied thickness
    const A = seismicBent(inp, ctx, spec, 0);
    const Lq = tLiq > 0 ? seismicBent(inp, ctx, spec, tLiq) : null;
    const pileRows = (S0) => {
      const st = solveBent(S0.B, inp).combo({ G: 1, Q: inp.psiE });
      return st.piles.map((p, k) => {
        let Ns = 0, Vs = 0, Ms = 0; // envelope over the elements of this pile
        for (const e of S0.sub.pileSubs[k]) { const m = S0.rsa.memberForces[e]; Ns = Math.max(Ns, m[0], m[6]); Vs = Math.max(Vs, Math.hypot(m[1], m[2]), Math.hypot(m[7], m[8])); Ms = Math.max(Ms, Math.hypot(m[4], m[5]), Math.hypot(m[10], m[11])); }
        const N = p.N - Ns, M = p.M + Ms, V = p.V + Vs; // gravity (tension +) with the seismic axial force taken as added compression
        const c = isoTubularCheck({ D: S0.B.P.D, t: S0.B.P.t, fy: S0.B.P.fy, E: S0.B.P.E, L: p.L, Ky: inp.headFix === 'pinned' ? 2.1 : 1.2, Kz: inp.headFix === 'pinned' ? 2.1 : 1.2, N, My: M, V });
        return { id: p.id, group: p.group, L: p.L, N, Nt: p.N + Ns, M, V, uc: cap(c.governing?.uc ?? Math.max(...c.checks.map((q) => q.uc))), gov: c.governing?.name ?? '' };
      });
    };
    const pA = pileRows(A), pL = Lq ? pileRows(Lq) : null;
    const govA = pA.reduce((m, p) => (p.uc > m.uc ? p : m), pA[0]), govL = pL ? pL.reduce((m, p) => (p.uc > m.uc ? p : m), pL[0]) : null;
    const T1 = A.modal.periods[0], Vb = A.rsa.baseShear.x, W = A.modal.totalMass.x * G;
    const hEff = Math.max(...pA.map((p) => p.L)), theta = W * A.u / Math.max(Vb * inp.q * hEff, 1);
    // time history: record matched to the elastic spectrum, generated on a fine step and integrated at the chosen step
    const dt0 = 0.01, stepN = Math.max(1, Math.round(inp.dtTH / dt0)), dt = stepN * dt0;
    const user = inp.record?.y?.length > 20 && inp.record?.x?.length === inp.record.y.length;
    let rec, acc;
    if (user) { // recorded motion: resampled to the analysis step by linear interpolation
      const tx = inp.record.x.map(Number), ay = inp.record.y.map(Number);
      need(tx.every((v, i) => Number.isFinite(v) && Number.isFinite(ay[i]) && (i === 0 || v > tx[i - 1])), 'The recorded acceleration history needs strictly increasing times and finite accelerations.');
      const T = tx[tx.length - 1] - tx[0]; need(T > 2 && T < 400, 'The recorded acceleration history must be between 2 s and 400 s long.');
      const resample = (step) => { const n = Math.floor(T / step) + 1, o = new Float64Array(n); let j = 0; for (let i = 0; i < n; i++) { const t = tx[0] + i * step; while (j < tx.length - 2 && tx[j + 1] < t) j++; o[i] = ay[j] + (ay[j + 1] - ay[j]) * (t - tx[j]) / (tx[j + 1] - tx[j]); } return o; };
      const fine = resample(dt0); rec = { acc: fine, pga: Math.max(...Array.from(fine, Math.abs)) }; acc = Array.from(resample(dt));
    } else {
      rec = run(() => spectrumCompatibleAccelerogram({ spectrum: elastic.Sa, duration: inp.duration, dt: dt0, zeta, nFreq: 50, iterations: 5, seed: 5 }), 'Synthetic accelerogram');
      acc = Array.from(rec.acc).filter((_, i) => i % stepN === 0);
    }
    const th = run(() => modalTimeHistory(A.modal, { dt, nSteps: acc.length - 1, zeta, baseAcc: { x: acc }, outputs: [{ node: A.B.headSea, dof: 0 }] }), 'Time-history analysis');
    let uTH = 0; for (const v of th.outputs[0]) uTH = Math.max(uTH, Math.abs(v));
    const slide = inp.ky * G < rec.pga ? newmarkDisplacement({ acc: rec.acc, dt: dt0, ky: inp.ky * G, bothPolarities: true }).displacement : 0;
    const resid = Math.max(slide, tLiq > 0 ? post.LDI : 0);
    const uGov = Lq ? Math.max(A.u, Lq.u) : A.u, ucGov = govL ? Math.max(govA.uc, govL.uc) : govA.uc;
    const checks = [
      ucCheck(`Pile ${govA.id} under gravity + earthquake (${govA.gov || 'combined'})`, govA.uc, 'ISO 19902:2020 §13.3; EN 1998-1:2004 §4.4.2', 'Axial force and bending in the seismic combination against member strength and stability. Above 1 the pile yields: acceptable only where it has been detailed and verified as a ductile hinge.'),
      ...(govL ? [ucCheck(`Pile ${govL.id} with liquefied layers (${govL.gov || 'combined'})`, govL.uc, 'EN 1998-5:2004 §5.4.2(6); ASCE 61-14 §4.7', 'The same check with lateral soil support removed over the liquefied thickness. The longer unsupported pile is more flexible and buckles at a lower load.')] : []),
      check('Deck displacement', uGov * 1e3, inp.uLimit, 'mm', 'EN 1998-1:2004 §4.3.4 (ds = q·de); PIANC WG34 (2001) Table 4.3', 'Peak sway of the deck in the design earthquake. Exceedance derails cranes, opens joints and demands more rotation at the pile heads than they can supply.'),
      sfCheck(`Liquefaction triggering${crit ? ` at ${f(crit.z - crit.thick / 2)} m` : ''}`, crit ? crit.FS : 5, inp.fsLiq, 'EN 1998-5:2004 §4.1.4, Annex B; Youd et al. (2001)', 'Cyclic resistance over cyclic stress in the most critical susceptible layer. Below the required factor the layer loses strength and stiffness during shaking: settlement, lateral spreading and loss of pile support follow.'),
      check('Residual ground displacement', resid * 1e3, inp.dLimit, 'mm', 'PIANC WG34 (2001) §4; Newmark (1965)', 'Permanent seaward movement of the slope from sliding-block analysis and, where layers liquefy, the lateral displacement index. It loads the piles kinematically and sets the repair needed.'),
      check('P–Δ sensitivity coefficient θ', theta, 0.3, '-', 'EN 1998-1:2004 §4.4.2.2', 'Overturning of the deck weight on the displaced bent relative to the seismic shear. Below 0.1 it is negligible; above 0.3 the bent is too flexible to be stable in the earthquake.'),
    ];
    const w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), ok = w.uc <= 1;
    const Ts = linspace(0.02, Math.max(4, 1.5 * T1), 60), mode = A.modal.modes[0], sc = deformScale(A.model, A.rsaE.displacements);
    const [tt, uu, aa] = thin([Array.from(th.t), Array.from(th.outputs[0], (v) => v * 1e3), acc], 500);
    const tNeed = inp.corr + (inp.t - inp.corr) * ucGov, steelBent = A.B.pileLen * A.B.P.massPerM / 1e3;
    const plateau = T1 >= elastic.TB && T1 <= elastic.TC;
    return {
      kpis: [
        { key: 'T1', label: 'Fundamental period', value: r(T1, 6), unit: 's', explain: 'Period of the sway mode of the bent with its seismic mass and the water moving with the piles. It places the wharf on the response spectrum.' },
        { key: 'Sa', label: 'Design spectral acceleration at T1', value: r(design.Sa(T1)), unit: 'm/s²', explain: 'Ordinate of the design spectrum (reduced by the behaviour factor) at the fundamental period.' },
        { key: 'Vb', label: 'Seismic base shear', value: r(Vb / 1e3), unit: 'kN', explain: `Design horizontal force on the bent from the modal combination: ${f(Vb / W * 100)} % of the seismic weight.` },
        { key: 'mass', label: 'Seismic mass of the bent', value: r(A.modal.totalMass.x / 1e3), unit: 't', explain: 'Deck, cargo share, crane share, piles and the water that moves with them. Force scales directly with it.' },
        { key: 'u', label: 'Deck displacement', value: r(A.u * 1e3), unit: 'mm', limit: inp.uLimit, sense: 'max', explain: 'Peak sway from the elastic spectrum (the equal-displacement estimate of the inelastic response), without liquefaction.' },
        { key: 'uTH', label: 'Deck displacement from the time history', value: r(uTH * 1e3), unit: 'mm', explain: `Peak sway from modal time-history analysis under ${user ? 'the recorded motion entered' : 'a synthetic record matched to the elastic spectrum'}. ${user ? 'It differs from the spectral value as far as the record departs from the design spectrum.' : 'It should agree with the spectral value to within the scatter of a single record (about ±15 %).'}` },
        { key: 'uLiq', label: 'Deck displacement with liquefied layers', value: r((Lq ? Lq.u : A.u) * 1e3), unit: 'mm', explain: 'Sway with lateral pile support removed over the liquefied thickness. Equal to the value above when nothing liquefies.' },
        { key: 'pileUc', label: 'Governing pile utilisation', value: r(ucGov), unit: '-', limit: 1, sense: 'max', explain: 'Largest ISO 19902 utilisation of any pile in the seismic combination, with and without liquefaction.' },
        { key: 'FSliq', label: 'Lowest factor of safety against liquefaction', value: r(crit ? crit.FS : 5), unit: '-', limit: inp.fsLiq, sense: 'min', explain: 'Cyclic resistance ratio over cyclic stress ratio in the weakest susceptible layer, with magnitude and overburden corrections. Reported as 5 when no layer is susceptible.' },
        { key: 'CSR', label: 'Cyclic stress ratio in the critical layer', value: r(crit ? crit.CSR : 0), unit: '-', explain: 'Earthquake shear stress over vertical effective stress, 0.65·(amax/g)·(σv/σ′v)·rd, at the base of the critical layer.' },
        { key: 'tLiq', label: 'Thickness of liquefied soil', value: r(tLiq), unit: 'm', explain: 'Total thickness of susceptible layers with a factor of safety below 1. The piles are taken as unsupported over this length.' },
        { key: 'LPI', label: 'Liquefaction potential index', value: r(lpi), unit: '-', explain: 'Depth-weighted severity over the top 20 m (Iwasaki): below 5 low, 5–15 high, above 15 severe surface effects.' },
        { key: 'settle', label: 'Post-liquefaction settlement', value: r(post.settlement * 1e3), unit: 'mm', explain: 'Free-field settlement from reconsolidation of the liquefied layers. It drags down on the piles and drops the backland relative to the piled deck.' },
        { key: 'resid', label: 'Residual ground displacement', value: r(resid * 1e3), unit: 'mm', limit: inp.dLimit, sense: 'max', explain: 'Larger of the sliding-block displacement of the slope and the lateral displacement index of the liquefied layers.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Response spectrum', x: { label: 'Period', unit: 's' }, y: { label: 'Spectral acceleration', unit: 'm/s²' },
          series: [{ name: 'Elastic', x: Ts.map((v) => r(v)), y: Ts.map((v) => r(elastic.Sa(v))) }, ...(inp.q > 1 ? [{ name: `Design (q = ${f(inp.q)})`, x: Ts.map((v) => r(v)), y: Ts.map((v) => r(design.Sa(v))), dash: true }] : []), { name: 'Fundamental period', x: [r(T1), r(T1)], y: [0, r(elastic.Sa(T1))], dash: true }, ...(Lq ? [{ name: 'With liquefied layers', x: [r(Lq.modal.periods[0]), r(Lq.modal.periods[0])], y: [0, r(elastic.Sa(Lq.modal.periods[0]))], dash: true }] : [])],
          explain: 'Spectrum for the ground type with the periods of the bent. On the plateau stiffening does not reduce force; on the falling branch flexibility lowers force but raises displacement.' },
        framePlot(A.model, { disp: mode, scale: (() => { let mx = 0; for (let i = 0; i < A.model.nodes.length; i++) mx = Math.max(mx, Math.hypot(mode[6 * i], mode[6 * i + 2])); return mx > 0 ? r(2 / mx, 3) : 1; })(), title: `First mode shape (T = ${f(T1)} s)`, values: A.model.nodes.map((_, i) => Math.abs(A.rsaE.displacements[6 * i]) * 1e3), valueLabel: 'Peak displacement', valueUnit: 'mm',
          explain: 'Sway mode of the bent, coloured by peak seismic displacement. The short landward piles restrain the deck and attract most of the shear.' }),
        { type: 'line', title: 'Deck displacement history', x: { label: 'Time', unit: 's' }, y: { label: 'Displacement (mm) / ground acceleration × 10 (m/s²)', unit: '-' },
          series: [{ name: 'Deck displacement (mm)', x: tt.map((v) => r(v)), y: uu.map((v) => r(v)) }, { name: 'Ground acceleration × 10', x: tt.map((v) => r(v)), y: aa.map((v) => r(v * 10)), dash: true }],
          explain: 'Elastic response of the deck to the ground motion used. The peak should be close to the spectral displacement; the number of large cycles indicates low-cycle demand on the pile heads.' },
        { type: 'line', title: 'Liquefaction factor of safety with depth', x: { label: 'Factor of safety', unit: '-' }, y: { label: 'Depth below the seabed', unit: 'm' },
          series: [{ name: 'FS', x: lay.map((q) => r(Math.min(q.FS, 3))), y: lay.map((q) => r(-(q.z - q.thick / 2))) }, { name: 'Required', x: [inp.fsLiq, inp.fsLiq], y: [0, r(-lay[lay.length - 1].z)], dash: true }],
          explain: 'Factor of safety of each layer (capped at 3; non-susceptible layers plot at the cap). Layers left of the dashed line are expected to liquefy in the design earthquake.' },
      ],
      tables: [
        { title: 'Liquefaction assessment by layer', columns: [{ key: 'z', label: 'Layer base', unit: 'm' }, { key: 'thick', label: 'Thickness', unit: 'm' }, { key: 'N', label: '(N1)60cs' }, { key: 'CSR', label: 'CSR' }, { key: 'CRR', label: 'CRR (M 7.5)' }, { key: 'FS', label: 'FS' }, { key: 'state', label: 'Result' }],
          rows: lay.map((q) => ({ z: r(q.z), thick: r(q.thick), N: q.sus ? r(q.N160cs, 3) : '–', CSR: r(q.CSR, 3), CRR: !q.sus ? '–' : Number.isFinite(q.CRR) ? r(q.CRR, 3) : 'too dense to liquefy', FS: q.sus ? r(q.FS, 3) : '–', state: !q.sus ? 'not susceptible' : q.FS < 1 ? 'liquefies' : q.FS < inp.fsLiq ? 'marginal' : 'stable' })),
          explain: 'One row per SPT layer: demand (CSR), resistance (CRR) and the factor of safety with magnitude scaling. Use it to target ground improvement at the layers that fail.' },
        { title: 'Pile forces in the seismic combination', columns: [{ key: 'id', label: 'Pile' }, { key: 'N', label: 'Compression', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }, { key: 'ucL', label: 'Utilisation, liquefied' }],
          rows: pA.map((p, k) => ({ id: p.group === 'raker' ? `${p.id} raker` : p.id, N: r(-p.N / 1e3), M: r(p.M / 1e3), V: r(p.V / 1e3), uc: r(p.uc, 3), ucL: pL ? r(pL[k].uc, 3) : '–' })),
          explain: 'Gravity plus design earthquake forces in each pile, with the utilisation as built and with liquefied layers. The short landward piles normally govern.' },
      ],
      notes: [`${plateau ? 'The bent sits on the spectral plateau' : T1 > elastic.TC ? 'The bent is on the descending branch of the spectrum' : 'The bent is stiffer than the plateau'}; surface acceleration ${f(amax / G)} g (soil factor ${f(elastic.S)}); modal mass participation ${f(A.rsa.massParticipation.x * 100)} %.`,
        `${tLiq > 0 ? `${f(tLiq)} m of soil liquefies: period lengthens from ${f(T1)} s to ${f(Lq.modal.periods[0])} s` : 'No susceptible layer liquefies'}; sliding-block displacement ${f(slide * 1e3)} mm for a yield acceleration of ${f(inp.ky)} g. Time-history over spectral displacement ${f(uTH / Math.max(A.u, 1e-9))}.`],
      advice: [
        ok ? eng('ok', `Period ${f(T1)} s, base shear ${f(Vb / 1e3)} kN, deck displacement ${f(uGov * 1e3)} mm, pile utilisation ${f(ucGov)}; lowest liquefaction factor ${f(crit ? crit.FS : 5)}, residual displacement ${f(resid * 1e3)} mm.`,
          'Confirm the pile-head connection can develop the moments in the table, add the kinematic pile demand from the residual ground displacement, and verify the behaviour factor by pushover before relying on it.', 'EN 1998-1:2004; EN 1998-5:2004; ASCE 61-14; PIANC WG34 (2001)',
          `A pile wall of about ${f(Math.max(tNeed, 10))} mm would reach unity in the seismic combination: ${steelCo2(steelBent * (Math.max(tNeed, 10) - inp.t) / inp.t)} per bent against the present wall.`)
          : eng('fail', `${w.name}: utilisation ${f(w.uc)}. ${tLiq > 0 ? `${f(tLiq)} m of soil liquefies (settlement ${f(post.settlement * 1e3)} mm). ` : ''}Deck displacement ${f(uGov * 1e3)} mm; residual ground displacement ${f(resid * 1e3)} mm.`,
            w.name.startsWith('Liquefaction') || w.name.includes('liquefied') || w.name.startsWith('Residual') ? `Improve the ground under and behind the wharf (vibro-compaction, stone columns or deep mixing to raise (N1)60cs above about ${f(crit ? crit.N160cs * inp.fsLiq / Math.max(crit.FS, 0.05) : 25)} in the critical layer), or design the piles for the unsupported length and kinematic loads.` : w.name.startsWith('Deck') || w.name.startsWith('P–Δ') ? 'Stiffen the bent (larger piles, shorter free length, rakers with ductile fuses) or articulate cranes and joints for the movement.' : `Increase the pile wall to about ${f(Math.ceil(tNeed))} mm, or detail the pile heads as ductile hinges and verify a higher behaviour factor by pushover.`,
            'EN 1998-5:2004 §4.1.4, §5.4; ASCE 61-14 §4; PIANC WG34 (2001)',
            `Thicker piles to ${f(Math.ceil(tNeed))} mm: ${steelCo2(steelBent * (Math.ceil(tNeed) - inp.t) / inp.t)} per bent. Ground improvement treats the cause and also removes the settlement of ${f(post.settlement * 1e3)} mm that the piles cannot.`),
        eco(ok ? 'info' : 'warn', `The design earthquake brings ${f(Vb / 1e3)} kN of base shear and ${f(resid * 1e3)} mm of residual ground movement per bent; seismic mass is ${f(A.modal.totalMass.x / 1e3)} t, of which the cargo share is ${f(inp.psiE * inp.qLive * 1e3 * A.B.sy * A.B.deckLen / G / 1e3)} t.`,
          'Agree the performance levels with the terminal operator (serviceable after the frequent event, repairable after the design event): they decide how much damage and downtime is carried as risk and how much is bought as ground improvement and steel.', `Performance-based design per PIANC WG34 and ASCE 61-14; ${ECON_BASIS}`),
      ],
      publish: { 'wharf.T1': T1, 'wharf.seismicBaseShear': Vb, 'wharf.seismicDisplacement': uGov, 'wharf.liquefiedThickness': tLiq },
    };
  },
  verify: [
    // One pile as a cantilever with a lumped deck mass: tube 1.0 × 25 mm, EI = 210e9·π/64·(1 − 0.95⁴) = 1.9121e9 N·m²;
    // L = 0.5 + 5 + 5 = 10.5 m; flexibility L³/(3EI) + L/(G·A/2) = 2.018e-7 + 3.4e-9 m/N → k = 4.873 MN/m;
    // m = 500 t → T = 2π√(m/k) = 2.013 s (pile and entrained water add under 1 %). Ground A, type 1 (TC = 0.4, TD = 2 s):
    // Sa = 2.5·ag·TC·TD/T² = 2.5·3·0.8/4.05 = 1.481 m/s² → V = m·Sa = 740 kN; u = Sa/ω² = 0.152 m.
    { name: 'Single-pile bent as a one-degree-of-freedom oscillator: T = 2π√(m/k), V = m·Sa', source: 'Chopra, Dynamics of Structures §6; EN 1998-1 Eq. 3.5 (hand calculation in the module source)',
      inputs: { nPiles: 1, overhang: 0, D: 1, t: 25, corr: 0, deckLevel: 0.5, depth: 5, fixDepth: 5, mCrane: 500, qLive: 0, gSdl: 0, ag: 3, ground: 'A', specType: '1', q: 1, spt: [{ z: 5, N160: 40, FC: 5, gamma: 20, sus: 0 }] },
      expect: [{ kpi: 'T1', value: 2.013, tol: 0.015 }, { kpi: 'Vb', value: 740, tol: 0.04 }, { kpi: 'u', value: 152, tol: 0.04 }] },
    // Cyclic stress ratio at 5 m in a submerged deposit of γ = 19 kN/m³, amax = 0.25 g on rock-like ground (S = 1):
    // σv = 95 kPa, σ′v = 95 − 9.81·5 = 45.95 kPa, rd = 1 − 0.00765·5 = 0.9618 → CSR = 0.65·0.25·(95/45.95)·0.9618 = 0.3231.
    { name: 'Cyclic stress ratio CSR = 0.65 (amax/g)(σv/σ′v) rd', source: 'Seed & Idriss (1971); Youd et al. (2001) Eq. 1–2 (hand calculation)',
      inputs: { ag: 2.4517, ground: 'A', spt: [{ z: 5, N160: 10, FC: 5, gamma: 19, sus: 1 }] }, expect: [{ kpi: 'CSR', value: 0.3231, tol: 0.005 }] },
  ],
};
