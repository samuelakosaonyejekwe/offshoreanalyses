// Structural reliability and life extension of the foundation: first-order reliability of the ultimate limit state
// at the mudline with corrosion over time, lognormal fatigue reliability, combined annual failure probability
// against the target, and the service life it supports.
import { form } from '../../engine/reliability/form.js';
import { normCdf, normInv } from '../../engine/reliability/sampling.js';
import { fatigueReliabilityLognormal, medianLifeFromCharacteristic } from '../../engine/fatigue/damage.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, tube, steelInputs, steelOf, LIMIT_CODE, CONFIRM } from './_common.js';

// Mean over characteristic (98 % annual fractile) for the load distribution families offered.
const meanFromChar = (kind, char, cov) => (kind === 'gumbel' ? char / (1 + 2.5923 * cov) : char / (1 + 2.0537 * cov));

export default {
  id: 's03-reliability-life-extension', suite: 3, title: 'Structural reliability and life extension',
  summary: 'Annual reliability index of the monopile in the extreme limit state as corrosion progresses, fatigue reliability with time, the combined annual failure probability against the code target, and the years of operation the structure supports beyond its design life.',
  covers: [16], fidelity: 'analytical',
  standards: ['IEC 61400-1:2019 Annex K and §7.6.1 (target reliability, annual pf 5·10⁻⁴)', 'DNV-ST-0126 (2021) §2.3 (safety classes, target annual pf 10⁻⁴)', 'DNV-ST-0262 (2016) and DNV-SE-0263 (lifetime extension of wind turbines)', 'ISO 2394:2015 (general principles on reliability)', 'JCSS Probabilistic Model Code (2001)', 'DNV-RP-C210 (2021) probabilistic inspection planning', 'Wirsching (1984) fatigue reliability'],
  models: ['First-order reliability method (Hasofer–Lind / Rackwitz–Fiessler) with Nataf transformation', 'Limit state: model uncertainty × yield or buckling strength × section modulus − load uncertainty × annual-maximum moment', 'Uniform corrosion wastage after coating life', 'Lognormal fatigue reliability format (Wirsching) on the Miner sum', 'Series system of ultimate and fatigue failure modes', 'Importance factors (α²) of the basic variables'],
  limits: `One governing section and two failure modes treated as independent; system effects, inspection updating, soil-capacity reliability and correlated wind–wave extremes are not modelled. Distribution types and coefficients of variation are editable defaults in the spirit of the JCSS model code and must be justified for the project. ${LIMIT_CODE}`,
  searchTerms: 'structural reliability FORM offshore wind monopile lifetime extension annual failure probability',
  inputs: [
    { key: 'D', label: 'Outer diameter at the governing section', unit: 'm', type: 'number', default: 6, min: 0.5, max: 14, group: 'Section', role: 'geometry', design: true },
    { key: 't', label: 'As-built wall thickness', unit: 'mm', type: 'number', default: 70, min: 10, max: 160, group: 'Section', role: 'geometry', design: true },
    ...steelInputs(),
    { key: 'resFactor', label: 'Characteristic buckling strength / yield', unit: '-', type: 'number', default: 0.85, min: 0.3, max: 1, group: 'Section', role: 'material', help: 'Ratio fks/fy from the buckling module for slender shells (1.0 when yield governs).' },
    { key: 'Mchar', label: 'Characteristic extreme moment (annual 98 % fractile)', unit: 'MN·m', type: 'number', default: 281, min: 1, max: 6000, group: 'Loads', role: 'bc', bind: 'owt.mudlineMoment', help: 'Fifty-year mudline moment from the extreme-loading module (unfactored).' },
    { key: 'loadDist', label: 'Distribution of the annual maximum moment', type: 'select', default: 'gumbel', group: 'Probabilistic model', role: 'load', options: [{ value: 'gumbel', label: 'Gumbel (extreme value type I)' }, { value: 'normal', label: 'Normal' }] },
    { key: 'covS', label: 'CoV of the annual maximum moment', unit: '-', type: 'number', default: 0.2, min: 0.02, max: 0.6, group: 'Probabilistic model', role: 'load', help: 'About 0.15–0.25 for storm-governed loads, lower for operation-governed loads with thrust limitation.' },
    { key: 'covXs', label: 'CoV of the load-model uncertainty', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.4, group: 'Probabilistic model', role: 'load', help: 'Aeroelastic and hydrodynamic model uncertainty (lognormal, median 1).' },
    { key: 'resDist', label: 'Distribution of the material strength', type: 'select', default: 'lognormal', group: 'Probabilistic model', role: 'material', options: [{ value: 'lognormal', label: 'Lognormal' }, { value: 'normal', label: 'Normal' }] },
    { key: 'fyBias', label: 'Mean yield / nominal yield', unit: '-', type: 'number', default: 1.12, min: 0.9, max: 1.4, group: 'Probabilistic model', role: 'material', help: 'Structural steel plate typically 1.08–1.15 (JCSS).' },
    { key: 'covFy', label: 'CoV of yield strength', unit: '-', type: 'number', default: 0.07, min: 0.01, max: 0.2, group: 'Probabilistic model', role: 'material' },
    { key: 'covXr', label: 'CoV of the resistance-model uncertainty', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.3, group: 'Probabilistic model', role: 'material', help: 'About 0.05 for yield of a compact section, 0.10–0.15 for shell buckling.' },
    { key: 'corr', label: 'Corrosion rate after coating breakdown', unit: 'mm/yr', type: 'number', default: 0.1, min: 0, max: 1, group: 'Degradation', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Each exposed face: about 0.1 mm/yr submerged with imperfect cathodic protection, 0.3 mm/yr in the splash zone (DNV-RP-0416).' },
    { key: 'tCoat', label: 'Coating / protection life', unit: 'years', type: 'number', default: 15, min: 0, max: 60, group: 'Degradation', role: 'ic' },
    { key: 'faces', label: 'Corroding faces', type: 'int', default: 2, min: 0, max: 2, group: 'Degradation', role: 'ic', help: 'Two when the inside is flooded and unprotected, one when the inside is sealed or protected.' },
    { key: 'age', label: 'Years in service', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Degradation', role: 'ic', help: 'Present age of the structure: the starting point of the remaining-life evaluation.' },
    { key: 'fatLife', label: 'Fatigue life with the design fatigue factor', unit: 'years', type: 'number', default: 80, min: 1, max: 5000, group: 'Fatigue', role: 'bc', bind: 'life.fatigueYears', help: 'From the fatigue module: years to a factored Miner sum of 1.' },
    { key: 'dff', label: 'Design fatigue factor used there', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Fatigue', role: 'criteria' },
    { key: 'mSN', label: 'S–N slope for the reliability format', unit: '-', type: 'number', default: 3, min: 3, max: 5, group: 'Fatigue', role: 'material' },
    { key: 'covB', label: 'CoV of the fatigue stress estimate', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 0.5, group: 'Fatigue', role: 'load', help: 'Uncertainty of computed stress ranges (loads, SCF, analysis); 0.15–0.25 for a design-stage estimate, about 0.10 when calibrated by strain monitoring.' },
    { key: 'designLife', label: 'Design life', unit: 'years', type: 'number', default: 25, min: 5, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'pfTarget', label: 'Target annual failure probability', unit: '-', type: 'number', default: 5e-4, min: 1e-6, max: 1e-2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-1 Annex K: 5·10⁻⁴ (β = 3.3) for unmanned turbines; DNV-ST-0126 normal safety class 10⁻⁴. ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const D = inp.D, t0 = inp.t / 1e3;
    const fyMean = steel.fy * inp.fyBias;
    const Mchar = inp.Mchar * 1e6, muS = meanFromChar(inp.loadDist, Mchar, inp.covS);
    const tAt = (yr) => Math.max(t0 - (inp.faces * inp.corr * Math.max(yr - inp.tCoat, 0)) / 1e3, 0.002);
    const det = (v) => ({ dist: 'deterministic', value: v });
    const vars = {
      fy: { dist: inp.resDist, mean: fyMean, cov: inp.covFy },
      S: { dist: inp.loadDist, mean: muS, cov: inp.covS },
      Xr: inp.covXr > 1e-6 ? { dist: 'lognormal', mean: 1, cov: inp.covXr } : det(1),
      Xs: inp.covXs > 1e-6 ? { dist: 'lognormal', mean: 1, cov: inp.covXs } : det(1),
    };
    const uls = (yr) => {
      const W = tube(D, tAt(yr)).W * inp.resFactor;
      const r = run(() => form((x) => x.Xr * x.fy * W - x.Xs * x.S, vars), 'Reliability analysis');
      need(Number.isFinite(r.beta), 'The reliability analysis did not return a finite index: check the distributions and the load level.');
      return r;
    };
    const horizon = Math.max(inp.designLife + 25, 50);
    const years = Array.from({ length: Math.round(horizon) }, (_, i) => i + 1);
    // ULS evaluated on a coarse grid in time (it changes only through corrosion) and interpolated.
    const grid = [0, inp.tCoat, ...[0.25, 0.5, 0.75, 1].map((q) => inp.tCoat + q * (horizon - inp.tCoat))].filter((v, i, a) => v <= horizon && a.indexOf(v) === i).sort((a, b) => a - b);
    const gb = grid.map((yr) => uls(yr).beta);
    const betaU = (yr) => { if (yr <= grid[0]) return gb[0]; for (let i = 1; i < grid.length; i++) if (yr <= grid[i]) return gb[i - 1] + ((gb[i] - gb[i - 1]) * (yr - grid[i - 1])) / (grid[i] - grid[i - 1]); return gb[gb.length - 1]; };
    const now = uls(inp.age), eol = uls(inp.designLife);
    // Fatigue: characteristic life = DFF × factored life; median from the scatter of the S–N curve.
    const charLife = inp.fatLife * inp.dff;
    const med = medianLifeFromCharacteristic(charLife);
    const fr = fatigueReliabilityLognormal({ medianLife: med.medianLife, m: inp.mSN, times: years, aCov: med.aCov, bCov: inp.covB });
    const pfU = years.map((y) => normCdf(-betaU(y)));
    const pfF = Array.from(fr.annualPf, (v) => (Number.isFinite(v) ? v : 0));
    const pfT = years.map((_, i) => pfU[i] + pfF[i]);
    const idx = (yr) => Math.min(Math.max(Math.round(yr) - 1, 0), years.length - 1);
    const iEnd = idx(inp.designLife), iNow = idx(Math.max(inp.age, 1));
    // First year in which the combined annual failure probability exceeds the target.
    let yTarget = horizon;
    for (let i = 0; i < years.length; i++) if (pfT[i] > inp.pfTarget) { yTarget = years[i]; break; }
    const ext = yTarget - inp.designLife, remain = Math.max(yTarget - inp.age, 0);
    const betaT = -normInv(inp.pfTarget);
    const betaEnd = -normInv(Math.min(Math.max(pfT[iEnd], 1e-15), 0.5));
    const imp = now.importance;
    const impRows = Object.entries(imp).filter(([, v]) => v > 1e-6).sort((a, b) => b[1] - a[1]);
    const names = { fy: 'Yield strength', S: 'Annual maximum moment', Xr: 'Resistance model', Xs: 'Load model' };
    const checks = [
      check('Annual failure probability at the end of the design life', pfT[iEnd], inp.pfTarget, '1/yr', 'IEC 61400-1:2019 Annex K; DNV-ST-0126 (2021) §2.3', 'Guards the code safety level: the combined annual probability of ultimate or fatigue failure in the last design year must not exceed the target for the safety class.'),
      check('Ultimate limit state today: target over achieved reliability index', betaT, Math.max(now.beta, 1e-9), 'β', 'ISO 2394:2015 §8; IEC 61400-1:2019 Annex K', 'Annual reliability of the corroded section against the extreme storm load. A lower index than the target means the section is too thin for the load or its uncertainty.'),
    ];
    checks[1].uc = Math.min(betaT / Math.max(now.beta, 1e-9), 99);
    const gov = worst(checks), ok = gov.uc <= 1;
    // Wall thickness that meets the target at the end of life in the ULS (β changes about linearly with ln t).
    const secNow = tube(D, t0);
    let tReq = t0;
    if (!ok) { for (let k = 0; k < 40; k++) { tReq *= 1.03; const W = tube(D, Math.max(tReq - (inp.faces * inp.corr * Math.max(inp.designLife - inp.tCoat, 0)) / 1e3, 0.002)).W * inp.resFactor; if (form((x) => x.Xr * x.fy * W - x.Xs * x.S, vars).beta >= betaT + 0.05) break; } }
    const dMass = (tube(D, tReq).A - secNow.A) * steel.rho * 10 / 1e3;
    return {
      kpis: [
        { key: 'betaUls', label: 'Annual reliability index, extreme limit state (today)', unit: '-', value: now.beta, limit: betaT, sense: 'min', explain: `Hasofer–Lind index of the governing section against the annual maximum moment at the present wall thickness. The target ${fmt(betaT)} corresponds to the annual failure probability ${inp.pfTarget}.` },
        { key: 'pfUls', label: 'Annual failure probability, extreme limit state (today)', unit: '1/yr', value: now.pf, explain: 'Probability that the extreme moment exceeds the section resistance in one year.' },
        { key: 'betaUlsEnd', label: 'Annual reliability index at the end of the design life', unit: '-', value: eol.beta, explain: 'Same index with the wall reduced by corrosion to the end of the design life.' },
        { key: 'pfFatEnd', label: 'Annual fatigue failure probability in the last design year', unit: '1/yr', value: pfF[iEnd], explain: 'Probability of fatigue failure occurring in the final year of the design life, given survival until then (lognormal format on the Miner sum).' },
        { key: 'betaFatEnd', label: 'Cumulative fatigue reliability index at the design life', unit: '-', value: fr.beta[iEnd], explain: 'Reliability index for no fatigue failure up to the end of the design life. Values of 2–3 are usual for details designed with a fatigue factor of 2–3.' },
        { key: 'pfTotalEnd', label: 'Combined annual failure probability at the design life', unit: '1/yr', value: pfT[iEnd], limit: inp.pfTarget, sense: 'max', explain: 'Ultimate plus fatigue annual failure probability in the last design year, the quantity compared with the code target.' },
        { key: 'betaTotalEnd', label: 'Combined annual reliability index at the design life', unit: '-', value: betaEnd, explain: 'Reliability index equivalent to the combined annual failure probability.' },
        { key: 'yearsToTarget', label: 'Service life supported at the target reliability', unit: 'years', value: yTarget, limit: inp.designLife, sense: 'min', explain: `Age at which the combined annual failure probability first exceeds the target (capped at the ${fmt(horizon)}-year horizon analysed).` },
        { key: 'extension', label: 'Possible life extension', unit: 'years', value: Math.max(ext, 0), explain: 'Years beyond the design life for which the target reliability is still met without intervention, before inspection results or monitoring data are used to update the estimate.' },
        { key: 'remaining', label: 'Remaining useful life from today', unit: 'years', value: remain, explain: 'Supported service life minus the present age.' },
        { key: 'tEnd', label: 'Wall thickness at the end of the design life', unit: 'mm', value: tAt(inp.designLife) * 1e3, explain: 'As-built wall less the corrosion wastage after the coating life on the stated number of faces.' },
        { key: 'alphaLoad', label: 'Importance of the load variables', unit: '%', value: ((imp.S || 0) + (imp.Xs || 0)) * 100, explain: 'Share of the uncertainty at the design point that comes from the load and the load model. A high share means better load data (monitoring) improves reliability more than better material data.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Annual failure probability with age', x: { label: 'Age', unit: 'years' }, y: { label: 'log10 annual failure probability', unit: '-' },
          series: [{ name: 'Extreme limit state', x: years, y: pfU.map((v) => Math.log10(Math.max(v, 1e-12))) }, { name: 'Fatigue', x: years, y: pfF.map((v) => Math.log10(Math.max(v, 1e-12))) }, { name: 'Combined', x: years, y: pfT.map((v) => Math.log10(Math.max(v, 1e-12))) },
            { name: 'Target', x: [1, horizon], y: [Math.log10(inp.pfTarget), Math.log10(inp.pfTarget)], dash: true }, { name: 'Design life', x: [inp.designLife, inp.designLife], y: [-12, Math.log10(Math.max(...pfT, inp.pfTarget)) + 0.5], dash: true }],
          explain: 'How the yearly risk of failure grows with age (logarithmic scale, floor at 10⁻¹²): the extreme limit state rises once corrosion starts, fatigue rises steadily. The crossing of the combined curve with the target is the supported service life.' },
        { type: 'bar', title: 'Importance factors at the design point', categories: impRows.map(([k]) => names[k] || k), series: [{ name: 'α², %', y: impRows.map(([, v]) => v * 100) }], x: { label: 'Basic variable' }, y: { label: 'Importance', unit: '%' },
          explain: 'Contribution of each uncertain quantity to the failure probability of the extreme limit state. Reduce the tallest bar — usually the annual maximum load — to gain reliability most efficiently.' },
        { type: 'line', title: 'Reliability index with age', x: { label: 'Age', unit: 'years' }, y: { label: 'Reliability index β', unit: '-' },
          series: [{ name: 'Extreme limit state, annual', x: years, y: years.map(betaU) }, { name: 'Fatigue, annual', x: years, y: Array.from(fr.annualBeta, (v) => (Number.isFinite(v) ? Math.min(v, 8) : 8)) }, { name: 'Target', x: [1, horizon], y: [betaT, betaT], dash: true }],
          explain: 'Annual reliability indices of the two failure modes against the target. Fatigue reliability falls throughout life; the extreme limit state is flat until the coating life ends.' },
      ],
      tables: [
        { title: 'Design point of the extreme limit state (today)', columns: [{ key: 'v', label: 'Variable' }, { key: 'm', label: 'Mean' }, { key: 'd', label: 'Design-point value' }, { key: 'a', label: 'α' }, { key: 'i', label: 'Importance α²', unit: '%' }],
          rows: Object.keys(vars).map((k) => ({ v: names[k], m: k === 'fy' ? `${fmt(fyMean / 1e6)} MPa` : k === 'S' ? `${fmt(muS / 1e6)} MN·m` : '1.0', d: k === 'fy' ? `${fmt(now.designPoint[k] / 1e6)} MPa` : k === 'S' ? `${fmt(now.designPoint[k] / 1e6)} MN·m` : fmt(now.designPoint[k]), a: round(now.alpha[k] ?? 0), i: round((imp[k] ?? 0) * 100) })),
          explain: 'Most probable combination of values at failure. Positive α marks a load-type variable, negative a resistance-type variable; the design-point values show which partial factors the case would need.' },
        { title: 'Reliability with age', columns: [{ key: 'y', label: 'Age', unit: 'years' }, { key: 't', label: 'Wall', unit: 'mm' }, { key: 'bu', label: 'β extreme' }, { key: 'pf', label: 'Annual pf fatigue' }, { key: 'pt', label: 'Annual pf combined' }, { key: 'ok', label: 'Meets target' }],
          rows: years.filter((y) => y % 5 === 0).map((y) => ({ y, t: round(tAt(y) * 1e3), bu: round(betaU(y)), pf: round(pfF[y - 1]), pt: round(pfT[y - 1]), ok: pfT[y - 1] <= inp.pfTarget ? 'yes' : 'no' })),
          explain: 'Five-yearly snapshot of wall thickness, reliability and failure probabilities: the basis of a lifetime-extension assessment and of the inspection interval.' },
      ],
      notes: [`Characteristic moment ${fmt(inp.Mchar)} MN·m → mean annual maximum ${fmt(muS / 1e6)} MN·m (${inp.loadDist}, CoV ${inp.covS}).`, `Characteristic fatigue life ${fmt(charLife)} years, median ${fmt(med.medianLife)} years.`,
        `${gov.name} governs (utilisation ${fmt(gov.uc)}). ${pfF[iEnd] > pfU[iEnd] ? 'Fatigue' : 'The extreme limit state'} dominates the risk at the end of the design life.`, `The target is met until year ${fmt(yTarget)}${yTarget >= horizon ? ' (the end of the analysed horizon)' : ''}.`],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `The target annual failure probability ${inp.pfTarget} is met through the design life (combined ${fmt(pfT[iEnd])} in year ${inp.designLife}) and until about year ${fmt(yTarget)}: ${fmt(Math.max(ext, 0))} years of extension are supported analytically.`,
          `Support a life-extension case with evidence: measured natural frequency and strain-based damage rates (to cut the stress CoV from ${inp.covB} towards 0.10), wall-thickness surveys in the splash zone, and inspection of the governing weld before year ${fmt(Math.max(Math.min(yTarget - 5, inp.designLife), 1))}.`,
          'DNV-ST-0262 (2016) §3–5; DNV-SE-0263; IEC 61400-1:2019 Annex K; DNV-RP-C210 (2021)',
          `Each year of extended operation avoids building replacement capacity: the ${fmt(secNow.A * steel.rho / 1e3)} t/m of installed steel keeps working with no new embodied carbon (replacement would be 1.9 t CO₂e per tonne at the screening factor).`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}. The target is exceeded from year ${fmt(yTarget)}.`,
            pfF[iEnd] > pfU[iEnd] ? 'Fatigue governs: reduce the stress uncertainty by monitoring, inspect and if necessary grind or repair the governing weld, or operate with load-reducing control modes in the last years.' : `The extreme limit state governs: a wall of about ${fmt(tReq * 1e3)} mm (now ${inp.t} mm) meets the target at the end of life; maintain the corrosion protection to stop the loss of ${fmt((t0 - tAt(inp.designLife)) * 1e3)} mm.`,
            'IEC 61400-1:2019 Annex K; DNV-ST-0126 (2021) §2.3; DNV-RP-0416 (corrosion protection); DNV-RP-C210 (2021)',
            `Thickening a 10 m section to ${fmt(tReq * 1e3)} mm adds ${fmt(dMass)} t — ${steelCo2Text(dMass)} Renewing anodes or coating preserves the existing wall for a small fraction of that footprint.`),
        adv('economic', ext > 0 ? 'info' : 'warn', ext > 0 ? `About ${fmt(ext)} years of operation beyond the ${inp.designLife}-year design life are supported at the target reliability, from a structure that is already paid for.` : `The reliability target is not held to the end of the ${inp.designLife}-year design life.`,
          ext > 0 ? 'Value the extension as added energy production against inspection, monitoring and re-certification cost; the load variables carry ' + fmt(((imp.S || 0) + (imp.Xs || 0)) * 100) + ' % of the uncertainty, so monitoring data are the most valuable evidence.' : 'Compare strengthening, corrosion-protection renewal and a shorter or derated operating period on net present value before committing to repair.',
          'Lifetime-extension business case (DNV-ST-0262; LCOE with extended operating life); value-of-information for structural health monitoring'),
      ],
    };
  },
  verify: [
    { name: 'Linear limit state with normal variables: β = (μR − μS)/√(σR² + σS²)', source: 'Cornell reliability index, exact for normal variables (Madsen, Krenk & Lind, Methods of Structural Safety, §4.1): W = 1.9110 m³ (D = 6 m, t = 70 mm), μR = 355 MPa × W = 678.4 MN·m, σR = 67.84 MN·m; μS = 300 MN·m (characteristic 423.2 MN·m at +2.0537σ), σS = 60 MN·m → β = 378.4/90.57 = 4.178',
      inputs: { D: 6, t: 70, steel: 's355', resFactor: 1, fyBias: 1, covFy: 0.1, resDist: 'normal', loadDist: 'normal', Mchar: 423.222, covS: 0.2, covXr: 0, covXs: 0, corr: 0 }, expect: [{ kpi: 'betaUls', value: 4.178, tol: 0.002 }] },
  ],
};
