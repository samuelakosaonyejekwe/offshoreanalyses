// Mudmat–soil interaction: contact pressure with lift-off, combined V–H–M–T envelope, embedment and cyclic degradation.
import { contactPressure, mudmatEnvelope, torsionalCapacityRect } from '../../engine/geotech/mudmat.js';
import { cyclicDegradationIndex, mccUndrainedStrengthRatio } from '../../engine/geotech/soilmodels.js';
import { DEG, need, run, fmt, round, adv, check, steelCo2Text, matInputs, loadInputs, NOT_FE } from './_common.js';

export default {
  id: 's04-mudmat-soil-interaction', suite: 4, title: 'Mudmat–soil interaction',
  summary: 'Contact pressure under the mat with lift-off, effective contact area, combined V–H–M–T failure envelope, local embedment and the effect of cyclic strength degradation.',
  covers: [5], fidelity: 'reduced-order',
  standards: ['ISO 19901-4:2016 §7', 'API RP 2GEO (2011) §7 and Annex A', 'DNV-RP-C212 (2019) §5–6'],
  models: ['Rigid mat on tensionless Winkler bed (nonlinear contact with separation)', 'Combined V–H–M–T loading envelope (generic interaction)', 'Green (1954) exact V–H interaction', 'Interface adhesion', 'Seabed penetration under bearing', 'Cyclic soil degradation (Idriss, Dobry & Singh)', 'Modified Cam-Clay undrained strength ratio', 'Critical-state soil mechanics'],
  limits: `The mat is rigid and the soil a bed of independent tensionless springs: contact is a plane with lift-off, not a continuum finite-element solution. The V–H–M–T envelope is a generic interaction with exact strip limits (Green; effective width) and user exponents — it is not the published Gourvenec (2007) or Feng et al. (2014) fit and must be calibrated before final design. ${NOT_FE}`,
  publishes: ['subsea.contactPmax', 'subsea.contactFraction', 'subsea.envelopeF'],
  searchTerms: 'mudmat failure envelope VHM undrained clay contact pressure lift-off',
  inputs: [
    ...matInputs(),
    { key: 'skirted', label: 'Skirted base (can sustain transient tension)', type: 'bool', default: false, group: 'Mudmat geometry', role: 'bc', help: 'Skirted mats in clay can carry short-term suction under the heel; unskirted mats lift off.' },
    ...loadInputs({ torsion: true, T: 200 }),
    { key: 'ML', label: 'Moment giving eccentricity along L', unit: 'kN·m', type: 'number', default: 300, min: 0, max: 1e6, group: 'Design loads at mudline', role: 'bc' },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 5, min: 0.2, max: 500, group: 'Soil', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.5, min: 0, max: 50, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 6, min: 2, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'alphaInt', label: 'Interface adhesion factor', unit: '-', type: 'number', default: 0.6, min: 0.05, max: 1, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.2 }, help: 'Interface over intact strength; about 1/sensitivity for steel plate, 1.0 for a skirted base.' },
    { key: 'ks', label: 'Vertical subgrade modulus (interface stiffness)', unit: 'kPa/m', type: 'number', default: 400, min: 10, max: 1e5, group: 'Soil', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 },
      help: 'Pressure per unit settlement of the bed. For a rigid rectangle on an elastic half-space ks ≈ E / (0.9·B·(1 − ν²)); calibrate against plate-load tests.' },
    { key: 'nCyc', label: 'Equivalent number of storm cycles', unit: '-', type: 'number', default: 1000, min: 1, max: 1e7, group: 'Cyclic loading', role: 'load' },
    { key: 'tDeg', label: 'Degradation parameter t', unit: '-', type: 'number', default: 0.04, min: 0, max: 0.4, group: 'Cyclic loading', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'δ = N^−t (Idriss, Dobry & Singh 1978). Depends on cyclic strain, plasticity and OCR: about 0.01–0.05 at small strain in plastic clay, higher for low-plasticity soil. From cyclic DSS tests.' },
    { key: 'mcM', label: 'Critical-state friction ratio M', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 1.8, group: 'Critical-state cross-check', role: 'soil' },
    { key: 'mcLambda', label: 'Compression index λ', unit: '-', type: 'number', default: 0.25, min: 0.02, max: 1, group: 'Critical-state cross-check', role: 'soil' },
    { key: 'mcKappa', label: 'Swelling index κ', unit: '-', type: 'number', default: 0.05, min: 0.002, max: 0.3, group: 'Critical-state cross-check', role: 'soil' },
    { key: 'ocr', label: 'Over-consolidation ratio', unit: '-', type: 'number', default: 1.5, min: 1, max: 50, group: 'Critical-state cross-check', role: 'ic' },
    { key: 'expA', label: 'Envelope exponent a (horizontal)', unit: '-', type: 'number', default: 2, min: 0.5, max: 10, group: 'Failure envelope', role: 'criteria', help: 'Exponent on H in (h/h*)^a + (m/m*)^b + (t/t*)^c = 1. Neutral default 2 — not a published fit; calibrate against finite-element or published envelopes.' },
    { key: 'expB', label: 'Envelope exponent b (moment)', unit: '-', type: 'number', default: 2, min: 0.5, max: 10, group: 'Failure envelope', role: 'criteria' },
    { key: 'expC', label: 'Envelope exponent c (torsion)', unit: '-', type: 'number', default: 2, min: 0.5, max: 10, group: 'Failure envelope', role: 'criteria' },
    { key: 'gammaM', label: 'Material factor on strength', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Failure envelope', role: 'criteria' },
    { key: 'nCell', label: 'Contact integration cells per side', type: 'int', default: 80, min: 8, max: 640, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nCell', label: 'Contact cells per side', levels: [20, 40, 80, 160], qoi: 'pMax' },
  compute(inp) {
    // SI: kN → N, kN·m → N·m, kPa → Pa, kPa/m → Pa/m.
    const { B, L } = inp, A = B * L;
    const V = inp.V * 1e3, H = inp.H * 1e3, MB = inp.M * 1e3, ML = inp.ML * 1e3, T = inp.T * 1e3;
    const ex = MB / V, ey = ML / V;
    need(ex < B / 2 && ey < L / 2, `The load resultant falls outside the mat (eccentricities ${fmt(ex)} m and ${fmt(ey)} m): the structure overturns. Enlarge the mat or reduce the moments.`);
    need(inp.mcKappa < inp.mcLambda, 'The swelling index κ must be smaller than the compression index λ.');
    const cp = run(() => contactPressure({ B, L, V, ex, ey, n: inp.nCell }), 'Contact pressure');
    // Strength, degraded for cyclic loading, at B/4 below the base and divided by the material factor.
    const delta = cyclicDegradationIndex(inp.nCyc, inp.tDeg);
    const suStat = (inp.su0 + inp.suGrad * Math.min(B, L) / 4) * 1e3;
    const su = suStat * delta / inp.gammaM;
    const bS = Math.min(B, L), lS = Math.max(B, L);
    const Vult = (2 + Math.PI) * su * (1 + 0.2 * bS / lS) * A;
    const Hult = (inp.skirted ? 1 : inp.alphaInt) * su * A;
    // Uniaxial moment capacity: Bransby & Randolph (1998) 0.69·A·B·su with full tension; (2+π)·A·B·su/8 without.
    const Mult = (inp.skirted ? 0.69 : (2 + Math.PI) / 8) * A * B * su;
    const Tult = torsionalCapacityRect({ B: bS, L: lS, su: (inp.skirted ? 1 : inp.alphaInt) * su });
    const Mres = Math.hypot(MB, ML * B / L); // moment about the other axis scaled to the same normalisation
    const env = (o = {}) => mudmatEnvelope({ V, H, M: Mres, T, Vult, Hult, Mult, Tult, tension: inp.skirted, a: inp.expA, b: inp.expB, c: inp.expC, ...o });
    const e = env();
    const fEnv = Number.isFinite(e.f) ? e.f : 99;
    // Local embedment: depth at which bearing on the contact area balances the peak pressure (su increasing with depth).
    const Nc = (2 + Math.PI) * (1 + 0.2 * bS / lS);
    const suNeed = cp.pMax / Nc;
    const suMud = inp.su0 * 1e3 * delta, kk = inp.suGrad * 1e3 * delta;
    const zEmb = suNeed <= suMud ? 0 : kk > 0 ? (suNeed - suMud) / kk : 99;
    const ks = inp.ks * 1e3;
    const sMax = cp.pMax / ks, sMean = V / A / ks;
    const tiltX = Math.atan(cp.gx / ks) / DEG, tiltY = Math.atan(cp.gy / ks) / DEG;
    const friction = Math.min(H, Hult);
    // Pressure field on a coarse grid for the heat map (tensionless plane).
    const ng = 24, xs = [], ys = [], z = [];
    for (let i = 0; i < ng; i++) xs.push(round(-B / 2 + (i + 0.5) * B / ng));
    for (let j = 0; j < ng; j++) { ys.push(round(-L / 2 + (j + 0.5) * L / ng)); z.push(xs.map((x) => round(Math.max(cp.p0 + cp.gx * x + cp.gy * (-L / 2 + (j + 0.5) * L / ng), 0) / 1e3))); }
    // V–H and V–M sections of the envelope with the design point.
    const vv = [], hh = [], mm = [];
    for (let i = 1; i < 40; i++) { const v = i / 40; const q = mudmatEnvelope({ V: v * Vult, Vult, Hult, Mult, tension: inp.skirted }); vv.push(round(v * Vult / 1e3)); hh.push(round(q.hStar * Hult / 1e3)); mm.push(round(q.mStar * Mult / 1e3)); }
    const cyc = [], dl = [];
    for (let i = 0; i <= 20; i++) { const n = 10 ** (i * 0.3); cyc.push(round(n)); dl.push(cyclicDegradationIndex(n, inp.tDeg)); }
    // Critical-state cross-check of the strength ratio at B/4: su/σ′v from Modified Cam-Clay (isotropic approximation).
    const sigV = inp.gammaSub * 1e3 * Math.max(Math.min(B, L) / 4, 0.1);
    const mccRatio = run(() => mccUndrainedStrengthRatio({ M: inp.mcM, lambda: inp.mcLambda, kappa: inp.mcKappa, ocr: inp.ocr }), 'Critical-state cross-check');
    const profRatio = suStat / sigV;
    const checks = [
      { name: 'Combined V–H–M–T envelope', demand: round(fEnv), capacity: 1, unit: 'envelope value', uc: Math.min(fEnv, 99), ref: 'Generic interaction with Green (1954) and effective-width limits; cf. ISO 19901-4:2016 §7.3 and API RP 2GEO Annex A', explain: 'Guards against combined bearing, sliding, overturning and twisting failure. Outside the envelope the soil yields and the mat moves along the governing mechanism.' },
      check('Vertical mobilisation V/Vult', V / 1e3, Vult / 1e3, 'kN', 'Prandtl / Skempton shape factor, degraded and factored strength', 'Vertical load against the uniaxial bearing capacity. Above about 0.5 the moment and horizontal capacities fall quickly.'),
      check('Peak contact pressure against local bearing', cp.pMax / 1e3, Nc * su / 1e3, 'kPa', 'Local bearing at the most loaded corner', 'Guards against local soil failure under the toe, which shows as progressive embedment and tilt of the mat.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const advice = [];
    advice.push(adv('engineering', ok ? 'ok' : 'fail',
      `Envelope value ${fmt(fEnv)} (≤ 1 inside), load factor on H, M and T at constant V to reach failure ${fmt(Math.min(e.loadFactor, 99))}; ${fmt(cp.contactFraction * 100)} % of the base is in contact and cyclic loading leaves ${fmt(delta * 100)} % of the static strength.`,
      ok ? (cp.fullContact ? 'Full contact is maintained. Confirm the envelope exponents against a finite-element or published envelope for this aspect ratio and strength profile before final design.' : `Part of the base lifts off: ${inp.skirted ? 'verify the transient suction assumed under the heel' : 'consider skirts or re-centre the load'} and confirm the envelope exponents against finite-element analysis.`)
        : `Reduce the eccentricity (${fmt(ex)} m along B, ${fmt(ey)} m along L), add skirts to mobilise tension and full interface strength, or enlarge the mat until the envelope value drops below 1.`,
      'ISO 19901-4:2016 §7.3; API RP 2GEO (2011) Annex A; Gourvenec (2007) and Feng et al. (2014) for published envelopes',
      'Skirts raise horizontal and moment capacity with a few tonnes of plate, where a larger mat would add steel in proportion to plan area and enlarge the disturbed seabed footprint.'));
    advice.push(adv('economic', 'info', `The result rests on the cyclic degradation (strength factor ${fmt(delta)}) and the interface factor ${fmt(inp.alphaInt)}; both are test-derived and move the envelope value roughly in inverse proportion.`,
      'A set of cyclic direct-simple-shear and interface shear tests on site samples is a small cost against mat steel: a 10 % gain in design strength reduces the required plan area by about 10 %, here ' + fmt(0.1 * A * 0.16) + ' t of steel at 0.16 t/m².',
      'Value-of-information reasoning; ISO 19901-8:2014 marine soil investigations',
      steelCo2Text(0.1 * A * 0.16)));
    return {
      kpis: [
        { key: 'pMax', label: 'Peak contact pressure', unit: 'kPa', value: cp.pMax / 1e3, explain: 'Largest soil pressure under the rigid mat, at the corner nearest the load resultant. It drives local bearing and plate design.' },
        { key: 'pMean', label: 'Mean contact pressure', unit: 'kPa', value: V / (A * cp.contactFraction) / 1e3, explain: 'Vertical load over the area actually in contact.' },
        { key: 'contact', label: 'Effective contact area', unit: 'm²', value: A * cp.contactFraction, explain: 'Plan area in compression. Less than the gross area means the heel has lifted off the seabed.' },
        { key: 'contactFrac', label: 'Contact fraction', unit: '-', value: cp.contactFraction, limit: 0.75, sense: 'min', explain: 'Share of the base in contact. Below about 0.75 the foundation is moment-dominated and sensitive to small load changes.' },
        { key: 'fEnv', label: 'Failure-envelope value', unit: '-', value: fEnv, limit: 1, sense: 'max', explain: 'Value of the V–H–M–T interaction: below 1 the load point is inside the failure envelope. Capped at 99 when V exceeds Vult.' },
        { key: 'loadFactor', label: 'Load factor to the envelope', unit: '-', value: Math.min(e.loadFactor, 99), explain: 'Factor on H, M and T at constant V that brings the load point onto the envelope: the true reserve under combined loading.' },
        { key: 'delta', label: 'Cyclic degradation index', unit: '-', value: delta, explain: 'Ratio of cyclic to static strength, N^−t. Values well below 1 mean storm loading softens the soil appreciably.' },
        { key: 'friction', label: 'Mobilised base shear', unit: 'kN', value: friction / 1e3, explain: 'Horizontal force carried by interface shear, limited by adhesion on the gross area.' },
        { key: 'zEmb', label: 'Local embedment at the toe', unit: 'm', value: zEmb, explain: 'Depth to which the most loaded edge must penetrate before the soil strength balances the peak pressure. Zero when the mudline strength suffices; 99 when it never does.' },
        { key: 'sMax', label: 'Winkler settlement at the toe', unit: 'mm', value: sMax * 1e3, explain: 'Peak pressure over the subgrade modulus: elastic compression of the bed under the most loaded corner.' },
        { key: 'tilt', label: 'Mat rotation (resultant)', unit: 'deg', value: Math.hypot(tiltX, tiltY), explain: 'Tilt of the rigid mat on the Winkler bed from the pressure gradient. Compare with the installation and connector tolerances.' },
        { key: 'mccRatio', label: 'Cam-Clay su/σ′v cross-check', unit: '-', value: mccRatio, explain: 'Undrained strength ratio predicted by Modified Cam-Clay for the given M, λ, κ and OCR. Compare with the ratio implied by the strength profile to spot inconsistent soil data.' },
      ],
      checks,
      plots: [
        { type: 'heatmap', title: 'Contact pressure under the mat', xs, ys, z, zLabel: 'Pressure', zUnit: 'kPa', x: { label: 'x along B', unit: 'm' }, y: { label: 'y along L', unit: 'm' },
          explain: 'Soil pressure on the base. Zero-pressure cells have lifted off; the hot corner is where local bearing and plate bending are checked.' },
        { type: 'line', title: 'V–H section of the failure envelope', x: { label: 'Vertical load', unit: 'kN' }, y: { label: 'Horizontal capacity', unit: 'kN' },
          series: [{ name: 'Envelope (M = T = 0)', x: vv, y: hh }, { name: 'Design point', x: [V / 1e3], y: [H / 1e3] }],
          explain: 'Horizontal capacity available at each vertical load (Green 1954). The design point must sit below the curve; note the plateau at low V and the drop near Vult.' },
        { type: 'line', title: 'V–M section of the failure envelope', x: { label: 'Vertical load', unit: 'kN' }, y: { label: 'Moment capacity', unit: 'kN·m' },
          series: [{ name: 'Envelope (H = T = 0)', x: vv, y: mm }, { name: 'Design point', x: [V / 1e3], y: [Mres / 1e3] }],
          explain: 'Moment capacity against vertical load. Without tension it peaks at half the vertical capacity and vanishes at zero load — a light structure has little moment resistance.' },
        { type: 'line', title: 'Cyclic strength degradation', x: { label: 'Number of cycles', unit: '-' }, y: { label: 'Degradation index δ', unit: '-' },
          series: [{ name: 'δ = N^−t', x: cyc, y: dl }, { name: 'Design cycles', x: [inp.nCyc], y: [delta] }],
          explain: 'Loss of strength with the number of load cycles for the chosen degradation parameter. The marked point is the storm considered.' },
      ],
      tables: [{ title: 'Uniaxial capacities and state', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Design strength (degraded, factored)', su / 1e3, 'kPa'], ['Vult', Vult / 1e3, 'kN'], ['Hult', Hult / 1e3, 'kN'], ['Mult', Mult / 1e3, 'kN·m'], ['Tult', Tult / 1e3, 'kN·m'], ['v = V/Vult', e.v, '-'], ['h* at this v', e.hStar, '-'], ['m* at this v', e.mStar, '-'],
          ['Corner pressures (+B,+L)', cp.corners.pp / 1e3, 'kPa'], ['(+B,−L)', cp.corners.pm / 1e3, 'kPa'], ['(−B,+L)', cp.corners.mp / 1e3, 'kPa'], ['(−B,−L)', cp.corners.mm / 1e3, 'kPa'], ['Mean Winkler settlement', sMean * 1e3, 'mm'], ['su/σ′v implied by the profile at B/4', profRatio, '-']]
          .map(([k, v, u]) => ({ k, v: round(v), u })),
        explain: 'The capacities that normalise the envelope, the reduced limits at the present vertical mobilisation, and the corner values of the contact-pressure plane (negative means lift-off at that corner).' }],
      notes: [cp.fullContact ? 'Full contact: pressure follows V/A·(1 ± 6e/B ± 6e/L).' : `Lift-off: ${fmt((1 - cp.contactFraction) * 100)} % of the base carries no pressure.`, `Envelope exponents a, b, c = ${inp.expA}, ${inp.expB}, ${inp.expC} (generic, to be calibrated).`,
        Math.abs(mccRatio - profRatio) / mccRatio > 0.5 ? `The strength profile implies su/σ′v = ${fmt(profRatio)} at B/4 against ${fmt(mccRatio)} from Cam-Clay: review the soil parameters.` : 'Strength profile and critical-state parameters are broadly consistent.'],
      advice,
      publish: { 'subsea.contactPmax': cp.pMax, 'subsea.contactFraction': cp.contactFraction, 'subsea.envelopeF': fEnv },
      quantities: {},
    };
  },
  verify: [
    { name: 'Full contact: pmax = V/A·(1 + 6e/B)', source: 'Navier linear pressure distribution: 1000/140·(1 + 6·1/10) = 11.429 kPa',
      inputs: { B: 10, L: 14, V: 1000, M: 1000, ML: 0, H: 0, T: 0 }, expect: [{ kpi: 'pMax', value: 11.4286, tol: 1e-4 }, { kpi: 'contactFrac', value: 1, tol: 1e-6 }] },
    { name: 'Uniaxial lift-off: pmax = 2V / (3L(B/2 − e)), contact width 3(B/2 − e)', source: 'Triangular pressure block outside the kern: e = 3 m → pmax = 2000/(3·14·2) = 23.81 kPa, contact 6/10',
      inputs: { B: 10, L: 14, V: 1000, M: 3000, ML: 0, H: 0, T: 0, nCell: 320 }, expect: [{ kpi: 'pMax', value: 23.8095, tol: 0.01 }, { kpi: 'contactFrac', value: 0.6, tol: 0.01 }] },
  ],
};
