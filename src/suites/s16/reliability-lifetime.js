// Structural integrity, reliability and lifetime assessment: first-order reliability of the ultimate limit state of
// the support structure, time-dependent fatigue reliability (lognormal format), corrosion life, component and
// system reliability of the turbine sub-assemblies and the resulting technical lifetime.
import { form } from '../../engine/reliability/form.js';
import { fatigueReliabilityLognormal } from '../../engine/fatigue/damage.js';
import { systemReliability, lifeReliability, mttf } from '../../engine/ram/rbd.js';
import { weibullFunctions } from '../../engine/ram/lifedata.js';
import { parisLawBs7910, parisLifeClosedForm } from '../../engine/fatigue/fracture.js';
import { normalQuantile } from '../../engine/aero/aep.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, LIMIT_CODE, CONFIRM, steelCo2Text, energyCo2Text } from './_common.js';

// Indicative Weibull lives to major replacement of main sub-assemblies (characteristic life η in years, shape β):
// order-of-magnitude values in the range of published offshore reliability surveys — replace with fleet data.
const COMPONENTS = [['Blades (set)', 70, 3], ['Pitch system', 28, 2], ['Main bearing', 45, 2.5], ['Gearbox', 30, 2.5], ['Generator', 45, 2], ['Power converter', 22, 1.5], ['Transformer', 55, 2.5], ['Yaw system', 60, 2], ['Array cable and termination', 80, 1.5]];
const lnPars = (mean, cov) => ({ dist: 'lognormal', mean, cov });

export default {
  id: 's16-reliability-lifetime', suite: 16, title: 'Structural reliability, component reliability and lifetime',
  summary: 'Annual failure probability of the support structure in the extreme load case, fatigue reliability over the service life, corrosion life, reliability of the turbine sub-assemblies and the technical lifetime of the system.',
  covers: [18], fidelity: 'analytical',
  standards: ['IEC 61400-1:2019 Annex K (target reliability, annual Pf 5·10⁻⁴)', 'DNV-ST-0126 (2021) §2 (safety philosophy, target annual Pf 10⁻⁴) and §4.11', 'ISO 2394:2015 General principles on reliability for structures', 'DNV-RP-C210 (2021) Probabilistic methods for planning of inspection for fatigue cracks', 'DNV-ST-0262 (2016) Lifetime extension of wind turbines', 'IEC 61078:2016 Reliability block diagrams', 'DNV-RP-0416 (2016) Corrosion protection for wind turbines'],
  models: ['First-order reliability method (HL-RF with Nataf transform)', 'Probabilistic reliability of the ultimate limit state with lognormal resistance and Gumbel or lognormal annual-maximum load', 'Fatigue damage reliability in the lognormal (Wirsching) format', 'Fracture mechanics: Paris-law crack growth (BS 7910)', 'Weibull life distributions and hazard', 'Series reliability block diagram of sub-assemblies', 'Uniform corrosion wastage against the corrosion allowance'],
  limits: `A single ultimate limit state (bending of the pile at the mudline) and a single fatigue detail represent the structure; soil, buckling and connection limit states, inspection updating and failure-assessment diagrams are not modelled; crack growth uses a constant geometry factor. Component lives are indicative Weibull parameters, independent and without repair (the O&M module treats repairable availability). ${LIMIT_CODE}`,
  publishes: ['life.systemYears', 'life.designYears', 'reliability.betaUls', 'reliability.pfFatigue'],
  searchTerms: 'offshore wind support structure reliability FORM fatigue reliability lifetime extension',
  inputs: [
    { key: 'D', label: 'Pile outer diameter at the critical section', unit: 'm', type: 'number', default: 6.5, min: 1, max: 14, group: 'Structure', role: 'geometry' },
    { key: 't', label: 'Wall thickness at the critical section', unit: 'mm', type: 'number', default: 80, min: 15, max: 160, group: 'Structure', role: 'geometry', design: true },
    { key: 'steel', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Structure', role: 'material' },
    { key: 'fyBias', label: 'Mean yield strength / nominal', unit: '-', type: 'number', default: 1.12, min: 1, max: 1.4, group: 'Resistance statistics', role: 'material', help: 'Structural steel plates: mean about 1.1–1.2 × the specified minimum (JCSS Probabilistic Model Code).' },
    { key: 'fyCov', label: 'Coefficient of variation of yield strength', unit: '-', type: 'number', default: 0.07, min: 0.01, max: 0.2, group: 'Resistance statistics', role: 'material' },
    { key: 'xrCov', label: 'Resistance model uncertainty (cov)', unit: '-', type: 'number', default: 0.05, min: 0, max: 0.3, group: 'Resistance statistics', role: 'material', help: 'Geometry and capacity-model uncertainty, lognormal with unit mean. 0 = deterministic.' },
    { key: 'Mk', label: 'Characteristic extreme mudline moment (50-year)', unit: 'MN·m', type: 'number', default: 325, min: 1, max: 5000, group: 'Load statistics', role: 'bc', bind: 'owt.mudlineMoment', help: 'From the extreme-loads module (unfactored governing case).' },
    { key: 'loadDist', label: 'Distribution of the annual maximum load', type: 'select', default: 'gumbel', group: 'Load statistics', role: 'bc', options: [{ value: 'gumbel', label: 'Gumbel (annual maxima)' }, { value: 'lognormal', label: 'Lognormal' }] },
    { key: 'loadBias', label: 'Mean annual maximum / characteristic value', unit: '-', type: 'number', default: 0.65, min: 0.2, max: 1.2, group: 'Load statistics', role: 'bc', help: 'For a Gumbel annual maximum with cov 0.2 the 98 % fractile (50-year value) is about 1.5 × the mean → 0.65.' },
    { key: 'loadCov', label: 'Coefficient of variation of the annual maximum load', unit: '-', type: 'number', default: 0.2, min: 0.02, max: 0.6, group: 'Load statistics', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Wind-dominated extremes about 0.15–0.25, wave-dominated 0.25–0.40.' },
    { key: 'xsCov', label: 'Load model uncertainty (cov)', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.4, group: 'Load statistics', role: 'bc', help: 'Aerodynamic, hydrodynamic and dynamic-response model uncertainty, lognormal with unit mean. 0 = deterministic.' },
    { key: 'fatLife', label: 'Calculated design fatigue life (with design fatigue factor)', unit: 'yr', type: 'number', default: 30, min: 1, max: 1e5, group: 'Fatigue', role: 'ic', bind: 'life.fatigueYears', help: 'From the fatigue module: life at a Miner sum of 1/DFF on the characteristic S–N curve.' },
    { key: 'dff', label: 'Design fatigue factor contained in that life', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Fatigue', role: 'criteria' },
    { key: 'mSn', label: 'S–N slope m', unit: '-', type: 'number', default: 3, min: 3, max: 5, group: 'Fatigue', role: 'material', help: 'The lognormal format assumes a single slope; 3 for the upper branch of welded-joint curves.' },
    { key: 'snCov', label: 'Scatter of the S–N intercept (cov)', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Fatigue', role: 'material', help: 'DNV-RP-C203 curves: standard deviation of log N = 0.20 → cov ≈ 0.49.' },
    { key: 'deltaCov', label: 'Scatter of the Miner sum at failure (cov)', unit: '-', type: 'number', default: 0.3, min: 0.05, max: 0.8, group: 'Fatigue', role: 'material' },
    { key: 'bCov', label: 'Fatigue stress model uncertainty (cov)', unit: '-', type: 'number', default: 0.15, min: 0.05, max: 0.5, group: 'Fatigue', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Load, stress-concentration and dynamic-response uncertainty on the stress ranges (0.15–0.30).' },
    { key: 'a0', label: 'Initial weld flaw depth', unit: 'mm', type: 'number', default: 0.5, min: 0.05, max: 10, group: 'Crack growth', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Initial condition of the fracture-mechanics check: 0.1–1 mm for welds accepted by non-destructive testing (DNV-RP-C210), or the size of a detected defect.' },
    { key: 'dsEq', label: 'Equivalent stress range at the detail', unit: 'MPa', type: 'number', default: 12, min: 1, max: 200, group: 'Crack growth', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Constant-amplitude range with the same damage as the load spectrum (m = 3), from the damage-equivalent moment of the fatigue module divided by the section modulus, with the stress concentration.' },
    { key: 'nCyc', label: 'Stress cycles per year', unit: '10⁶/yr', type: 'number', default: 6, min: 0.1, max: 60, group: 'Crack growth', role: 'bc', help: 'About the first natural frequency × seconds per year × operating fraction (0.25 Hz → 7·10⁶).' },
    { key: 'fmEnv', label: 'Crack-growth law (BS 7910 simplified, steel)', type: 'select', default: 'air', group: 'Crack growth', role: 'environment', options: [{ value: 'air', label: 'In air, or sea water with effective cathodic protection at low ΔK' }, { value: 'marine', label: 'Marine environment (unprotected, upper bound)' }], help: 'The marine law is 4.4 times faster than the air law at the same ΔK.' },
    { key: 'yGeom', label: 'Crack geometry factor Y (including weld-toe magnification)', unit: '-', type: 'number', default: 1.2, min: 0.5, max: 3, group: 'Crack growth', role: 'material', help: 'K = Y·Δσ·√(πa); 1.12 for a shallow edge crack, higher near the weld toe.' },
    { key: 'age', label: 'Present age of the asset', unit: 'yr', type: 'number', default: 0, min: 0, max: 60, group: 'Service life', role: 'ic', help: 'Initial condition: 0 for a new design, the years in service for a life-extension assessment.' },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 5, max: 60, group: 'Service life', role: 'criteria', bind: 'life.designYears' },
    { key: 'corrAllow', label: 'Corrosion allowance', unit: 'mm', type: 'number', default: 6, min: 0, max: 20, group: 'Corrosion', role: 'material' },
    { key: 'corrRate', label: 'Corrosion rate after coating breakdown', unit: 'mm/yr', type: 'number', default: 0.3, min: 0.01, max: 1, group: 'Corrosion', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'DNV-RP-0416: 0.30 mm/yr in the splash zone (temperate), 0.10 mm/yr submerged with cathodic protection failure.' },
    { key: 'coatLife', label: 'Coating life before corrosion starts', unit: 'yr', type: 'number', default: 15, min: 0, max: 40, group: 'Corrosion', role: 'material' },
    { key: 'components', label: 'Sub-assemblies: Weibull life to major replacement', type: 'table', group: 'Component reliability', role: 'material', columns: [{ key: 'name', label: 'Sub-assembly', type: 'text' }, { key: 'eta', label: 'Characteristic life η', unit: 'yr', type: 'number' }, { key: 'beta', label: 'Shape β', unit: '-', type: 'number' }],
      default: COMPONENTS.map(([name, eta, beta]) => ({ name, eta, beta })), help: 'Indicative values: replace with fleet failure data (maintenance history). β > 1 is wear-out, β = 1 random failures.' },
    { key: 'pfTarget', label: 'Target annual failure probability', unit: '-', type: 'number', default: 5e-4, min: 1e-6, max: 1e-2, group: 'Acceptance', role: 'criteria', help: `IEC 61400-1:2019 Annex K: 5·10⁻⁴ per year (β = 3.3) for unmanned wind turbines; DNV-ST-0126 normal safety class: 10⁻⁴ (β = 3.7). ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0, `Steel '${inp.steel}' is not in the materials library.`);
    const D = inp.D, t = inp.t / 1e3;
    need(t < D / 4, 'Wall thickness must be below a quarter of the diameter.');
    const W = Math.PI / 32 * (D ** 4 - (D - 2 * t) ** 4) / D;
    const Rmean = st.fy * inp.fyBias * W, Smean = inp.Mk * 1e6 * inp.loadBias;
    // ---- ULS: g = Xr·fy·W − Xs·M ----
    const vars = { fy: lnPars(st.fy * inp.fyBias, inp.fyCov), M: inp.loadDist === 'gumbel' ? { dist: 'gumbel', mean: Smean, cov: inp.loadCov } : lnPars(Smean, inp.loadCov),
      ...(inp.xrCov > 0 ? { Xr: lnPars(1, inp.xrCov) } : {}), ...(inp.xsCov > 0 ? { Xs: lnPars(1, inp.xsCov) } : {}) };
    const g = (x) => ((x.Xr ?? 1) * x.fy * W - (x.Xs ?? 1) * x.M) / Rmean;
    const fr = run(() => form(g, vars), 'First-order reliability analysis');
    need(Number.isFinite(fr.beta), 'The reliability analysis did not converge: check the load and resistance statistics.');
    const betaT = -normalQuantile(inp.pfTarget);
    // deterministic check with partial factors for comparison
    const ucDet = inp.Mk * 1e6 * 1.35 / (st.fy * W / 1.1);
    // ---- fatigue reliability ----
    const charLife = inp.fatLife * inp.dff;
    // characteristic S–N curve is the mean minus two standard deviations of log N: median life = characteristic × 10^(2·s), s from the intercept cov
    const sLog = Math.sqrt(Math.log(1 + inp.snCov ** 2)) / Math.LN10, medianLife = charLife * 10 ** (2 * sLog);
    const tEnd = Math.max(inp.designLife * 1.6, inp.age + 10), times = linspace(1, Math.ceil(tEnd), Math.ceil(tEnd));
    const fat = run(() => fatigueReliabilityLognormal({ medianLife, m: inp.mSn, times, deltaCov: inp.deltaCov, aCov: inp.snCov, bCov: inp.bCov }), 'Fatigue reliability');
    const iEnd = Math.min(times.length - 1, Math.max(0, Math.round(inp.designLife) - 1));
    const pfFatEnd = fat.pf[iEnd], betaFatEnd = fat.beta[iEnd], annPfEnd = fat.annualPf[iEnd];
    // year in which the annual fatigue failure probability first exceeds the target
    let yearFat = times[times.length - 1] + 1;
    for (let i = 0; i < times.length; i++) if (fat.annualPf[i] > inp.pfTarget) { yearFat = times[i]; break; }
    const fatReached = yearFat <= times[times.length - 1];
    // ---- fracture mechanics: Paris growth of a weld flaw to half the wall thickness (BS 7910 simplified law, threshold not credited) ----
    const law = parisLawBs7910({ environment: inp.fmEnv, curve: 'simplified', variableAmplitude: true });
    need(inp.a0 / 1e3 < 0.45 * t, 'The initial flaw is deeper than 45 % of the wall: this is a fitness-for-service case, not a design check.');
    const nCrack = parisLifeClosedForm({ a0: inp.a0 / 1e3, af: 0.5 * t, stressRange: inp.dsEq * 1e6, Y: inp.yGeom, A: law.stages[0].A, m: law.stages[0].m });
    const crackLife = nCrack / (inp.nCyc * 1e6);
    // ---- corrosion ----
    const corrLife = inp.coatLife + (inp.corrRate > 0 ? inp.corrAllow / inp.corrRate : 1e4);
    // ---- components and system ----
    const comps = (inp.components || []).filter((c) => c.eta > 0 && c.beta > 0);
    need(comps.length >= 1, 'The component table needs at least one row with a positive characteristic life and shape.');
    const block = { type: 'series', children: comps.map((c) => ({ type: 'component', id: c.name, life: { dist: 'weibull', shape: c.beta, scale: c.eta } })) };
    const Rsys = (tt) => systemReliability(block, tt);
    const rEnd = Rsys(inp.designLife), mt = run(() => mttf(block, { scale: 10 }), 'System MTTF');
    const expRepl = comps.map((c) => (inp.designLife / c.eta) ** c.beta);   // expected major replacements with minimal repair (cumulative hazard)
    const totRepl = expRepl.reduce((a, b) => a + b, 0);
    const structLife = Math.min(fatReached ? yearFat : 1e4, corrLife, inp.fatLife);   // crack-growth life is reported separately: it governs inspection, not the uninspected design life
    const remaining = Math.max(0, structLife - inp.age);
    const checks = [
      check('Annual reliability of the pile in the extreme load case', betaT, fr.beta, '-', `ISO 2394:2015; DNV-ST-0126 §2.3 / IEC 61400-1 Annex K (target annual Pf = ${inp.pfTarget})`, 'Demand is the target reliability index, capacity the index from the FORM analysis of yield at the critical section. Below the target the annual probability of structural failure exceeds what the code calibration assumes.'),
      check('Annual fatigue failure probability at the end of the design life', annPfEnd, inp.pfTarget, '1/yr', 'DNV-RP-C210 (2021) §4 (target reliability for fatigue); Wirsching (1984)', 'Probability that the fatigue detail fails in the last year of the design life against the target. Exceeding it calls for inspection, a higher design fatigue factor or weld improvement.'),
      check('Crack-growth life of an initial weld flaw against the design life', inp.designLife, Math.min(crackLife, 1e5), 'yr', 'BS 7910:2019 clause 8 (fatigue crack growth); DNV-RP-C210 (2021) §5 (fracture-mechanics model)', 'Design life against the time for the assumed initial flaw to grow to half the wall thickness. A shorter crack-growth life means the detail depends on in-service inspection to remain safe.'),
      check('Corrosion life against the design life', inp.designLife, corrLife, 'yr', 'DNV-RP-0416 (2016) §5 (corrosion allowance = rate × (design life − coating life))', 'Design life against the years until the corrosion allowance is consumed after coating breakdown. Beyond that the wall is thinner than the strength and fatigue design assumed.'),
      check('Structural lifetime against the design life', inp.designLife, structLife, 'yr', 'DNV-ST-0262 (2016) lifetime extension; DNV-ST-0126 §4.11', 'Design life against the earliest of the fatigue design life, the year the fatigue reliability falls below target and the corrosion life: the technical life of the primary structure without intervention.'),
    ];
    const gov = worst(checks);
    const ts = linspace(0, Math.ceil(tEnd), 41);
    const imp = Object.entries(fr.importance || {}).map(([k, v]) => ({ k, v }));
    const label = { fy: 'Yield strength', M: 'Annual maximum load', Xr: 'Resistance model', Xs: 'Load model' };
    const steelPerM = Math.PI * (D - t) * t * 7.85;
    // thickness for the target reliability: resistance ∝ t, β changes by about ln(ratio)/σ_ln with σ_ln ≈ combined log standard deviation
    const sigLn = Math.sqrt(Math.log(1 + inp.fyCov ** 2) + Math.log(1 + inp.loadCov ** 2) + Math.log(1 + inp.xrCov ** 2) + Math.log(1 + inp.xsCov ** 2));
    const tReq = inp.t * Math.exp((betaT - fr.beta) * sigLn);
    return {
      kpis: [
        { key: 'betaUls', label: 'Annual reliability index, ultimate limit state', value: round(fr.beta, 4), unit: '-', limit: round(betaT, 4), sense: 'min', explain: 'First-order reliability index of yielding at the critical pile section in any one year. 3.7 corresponds to an annual failure probability of 10⁻⁴.' },
        { key: 'pfUls', label: 'Annual failure probability, ultimate limit state', value: Number(fr.pf.toExponential(3)), unit: '1/yr', explain: 'Probability that the annual maximum load exceeds the resistance, Φ(−β). Compare with the target for the consequence class.' },
        { key: 'ucDet', label: 'Deterministic utilisation (γf = 1.35, γM = 1.10)', value: round(ucDet, 4), unit: '-', explain: 'Partial-factor check of the same section for comparison: a utilisation of 1.0 should correspond roughly to the target reliability if the statistics match the code calibration.' },
        { key: 'betaFatEnd', label: 'Cumulative fatigue reliability index at the design life', value: round(betaFatEnd, 4), unit: '-', explain: `Reliability index for fatigue failure at any time up to year ${inp.designLife} (cumulative failure probability ${Number(pfFatEnd.toExponential(2))}).` },
        { key: 'annPfFat', label: 'Annual fatigue failure probability in the last design year', value: Number(annPfEnd.toExponential(3)), unit: '1/yr', limit: inp.pfTarget, sense: 'max', explain: 'Increment of the fatigue failure probability in the final year of the design life: the quantity compared with the annual target.' },
        { key: 'medianLife', label: 'Median fatigue life', value: round(Math.min(medianLife, 1e6), 5), unit: 'yr', explain: `Life on the median S–N curve with a Miner sum of 1: characteristic life (${fmt(charLife)} yr) times the mean-to-characteristic factor of the S–N scatter.` },
        { key: 'yearFat', label: 'Year the fatigue reliability falls below target', value: round(fatReached ? yearFat : times[times.length - 1], 4), unit: 'yr', explain: fatReached ? 'First year in which the annual fatigue failure probability exceeds the target: the latest time for a first inspection or the end of the uninspected life.' : `Not reached within ${fmt(times[times.length - 1])} years (the value shown is the end of the evaluated period).` },
        { key: 'crackLife', label: 'Crack-growth life of the weld flaw to half thickness', value: round(Math.min(crackLife, 1e5), 5), unit: 'yr', explain: `Years for a ${inp.a0} mm flaw to grow to half the wall under the equivalent stress range by the Paris law (BS 7910 simplified law for steel, selected environment). It supports the S–N result and sets inspection intervals: a detectable crack must not reach critical size between inspections.` },
        { key: 'corrLife', label: 'Corrosion life', value: round(Math.min(corrLife, 1e4), 4), unit: 'yr', explain: 'Coating life plus the time to consume the corrosion allowance at the stated rate.' },
        { key: 'systemYears', label: 'Technical lifetime of the primary structure', value: round(Math.min(structLife, 1e4), 4), unit: 'yr', limit: inp.designLife, sense: 'min', explain: 'Earliest of fatigue design life, fatigue-reliability limit and corrosion life. Sub-assemblies of the turbine are replaceable and do not end the life of the asset.' },
        { key: 'remaining', label: 'Remaining structural life from the present age', value: round(Math.min(remaining, 1e4), 4), unit: 'yr', explain: 'Technical lifetime minus the years already in service: the basis of a lifetime-extension case.' },
        { key: 'rSys', label: 'Probability of no major replacement within the design life', value: round(rEnd, 5), unit: '-', explain: 'Series reliability of the sub-assemblies: the chance that a turbine reaches the end of the design life without any major component exchange. It is normally small — major replacements are part of the plan.' },
        { key: 'mttf', label: 'Mean time to the first major replacement', value: round(mt, 4), unit: 'yr', explain: 'Expected time until the first sub-assembly needs a major replacement on one turbine.' },
        { key: 'replacements', label: 'Expected major replacements per turbine in the design life', value: round(totRepl, 4), unit: '-', explain: 'Sum of the cumulative Weibull hazards of the sub-assemblies: the number of major exchanges to budget for, each needing a jack-up or heavy-lift vessel.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Fatigue reliability over the service life', x: { label: 'Year in service', unit: 'yr' }, y: { label: 'Reliability index', unit: '-' }, series: [{ name: 'Annual β', x: times, y: Array.from(fat.annualBeta, (v) => round(Math.min(v, 8), 4)) }, { name: 'Cumulative β', x: times, y: Array.from(fat.beta, (v) => round(Math.min(v, 8), 4)) }, { name: 'Target annual β', x: [times[0], times[times.length - 1]], y: [betaT, betaT], dash: true }, { name: 'Design life', x: [inp.designLife, inp.designLife], y: [0, 8], dash: true }],
          explain: 'Annual and cumulative fatigue reliability index against time (capped at 8). Where the annual curve crosses the target line, the detail needs inspection or the uninspected life ends.' },
        { type: 'line', title: 'Reliability of the sub-assemblies and of the turbine', x: { label: 'Year in service', unit: 'yr' }, y: { label: 'Probability of no major replacement', unit: '-' }, series: [{ name: 'All sub-assemblies (series)', x: ts.map((v) => round(v, 3)), y: ts.map((v) => round(Rsys(v), 5)) }, ...comps.slice(0, 8).map((c) => ({ name: c.name, x: ts.map((v) => round(v, 3)), y: ts.map((v) => round(lifeReliability({ dist: 'weibull', shape: c.beta, scale: c.eta }, v), 5)), dash: true }))],
          explain: 'Survival curves to the first major replacement for each sub-assembly and for the turbine as a series system. The steepest curves identify the components that drive heavy-maintenance campaigns.' },
        { type: 'bar', title: 'Importance of the random variables in the ultimate limit state', categories: imp.map((q) => label[q.k] || q.k), series: [{ name: 'α² (share of variance)', y: imp.map((q) => round(q.v, 4)) }], y: { label: 'Importance', unit: '-' },
          explain: 'Squared FORM sensitivity factors: the share of each variable in the uncertainty of the limit state. The load usually dominates — better site and load data raise reliability more than stronger steel.' },
        { type: 'bar', title: 'Expected major replacements per turbine in the design life', categories: comps.map((c) => c.name), series: [{ name: 'Expected number', y: expRepl.map((v) => round(v, 4)) }], y: { label: 'Replacements', unit: '-' },
          explain: 'Cumulative hazard of each sub-assembly over the design life. Values near or above 1 mean the component is replaced at least once per turbine.' },
      ],
      tables: [
        { title: 'Sub-assemblies', columns: [{ key: 'name', label: 'Sub-assembly' }, { key: 'eta', label: 'η', unit: 'yr' }, { key: 'beta', label: 'β', unit: '-' }, { key: 'R', label: 'Survival to the design life', unit: '-' }, { key: 'h', label: 'Hazard in the last design year', unit: '1/yr' }, { key: 'n', label: 'Expected replacements', unit: '-' }],
          rows: comps.map((c, i) => { const w = weibullFunctions(inp.designLife, { shape: c.beta, scale: c.eta }); return { name: c.name, eta: c.eta, beta: c.beta, R: round(w.reliability, 5), h: round(w.hazard, 5), n: round(expRepl[i], 4) }; }),
          explain: 'Weibull parameters, survival probability, end-of-life hazard and expected number of major replacements of each sub-assembly. Feed the hazards into the O&M and spare-parts planning.' },
        { title: 'Design point of the ultimate limit state', columns: [{ key: 'v', label: 'Variable' }, { key: 'mean', label: 'Mean' }, { key: 'dp', label: 'Design point' }, { key: 'alpha', label: 'α' }, { key: 'u', label: 'Unit' }],
          rows: Object.keys(vars).map((k) => ({ v: label[k] || k, mean: k === 'fy' ? round(st.fy * inp.fyBias / 1e6, 5) : k === 'M' ? round(Smean / 1e6, 5) : 1, dp: k === 'fy' ? round(fr.designPoint[k] / 1e6, 5) : k === 'M' ? round(fr.designPoint[k] / 1e6, 5) : round(fr.designPoint[k], 5), alpha: round(fr.alpha[k], 4), u: k === 'fy' ? 'MPa' : k === 'M' ? 'MN·m' : '-' })),
          explain: 'Most likely combination of the variables at failure. Negative α marks resistance-type variables, positive α load-type; the design-point values correspond to implied partial factors.' },
      ],
      notes: [
        `ULS: mean resistance ${fmt(Rmean / 1e6, 4)} MN·m against a mean annual maximum of ${fmt(Smean / 1e6, 4)} MN·m (${inp.loadDist}, cov ${inp.loadCov}); β = ${fmt(fr.beta)} (Pf = ${fr.pf.toExponential(2)} per year), target ${fmt(betaT)}; FORM ${fr.converged ? 'converged' : 'did not fully converge'}.`,
        `Fatigue: design life ${fmt(inp.fatLife)} yr with DFF ${inp.dff} → characteristic ${fmt(charLife)} yr, median ${fmt(medianLife)} yr; annual Pf in year ${inp.designLife}: ${annPfEnd.toExponential(2)}.`,
        `Lifetime: fatigue ${fmt(inp.fatLife)} yr, reliability limit ${fatReached ? `${fmt(yearFat)} yr` : 'beyond the evaluated period'}, corrosion ${fmt(Math.min(corrLife, 1e4))} yr → technical lifetime ${fmt(Math.min(structLife, 1e4))} yr; remaining ${fmt(Math.min(remaining, 1e4))} yr at age ${inp.age}.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Annual reliability index ${fmt(fr.beta)} against a target of ${fmt(betaT)}; annual fatigue failure probability in year ${inp.designLife} is ${annPfEnd.toExponential(2)} against ${inp.pfTarget}; technical lifetime ${fmt(Math.min(structLife, 1e4))} years.`,
          gov.uc > 1 ? (gov === checks[0] ? `Raise the resistance: a wall of about ${fmt(tReq, 3)} mm at the critical section reaches the target, or reduce load uncertainty with site-specific extremes (the load carries ${fmt((fr.importance?.M ?? 0) * 100)} % of the variance).` : gov === checks[1] ? `Plan a first inspection of the fatigue detail before year ${fmt(yearFat)}, or raise the design fatigue factor (life ${fmt(inp.fatLife * 1.5)} yr would be needed).` : gov === checks[2] ? `The assumed flaw grows through half the wall in ${fmt(crackLife)} years: reduce the stress range (thicker wall, ground weld), tighten the acceptance flaw size, or inspect at intervals of about ${fmt(crackLife / 3)} years.` : gov === checks[3] ? `Increase the corrosion allowance to ${fmt(inp.corrRate * Math.max(inp.designLife - inp.coatLife, 0))} mm or extend the coating / cathodic-protection life.` : 'Lengthen the governing life (fatigue, corrosion) or plan inspection and repair before it is reached.') : `The structure meets the reliability targets with a technical life of ${fmt(Math.min(structLife, 1e4))} years: a lifetime-extension case of about ${fmt(Math.max(0, Math.min(structLife, 1e4) - inp.designLife))} years can be examined with inspection data and measured loads.`,
          'ISO 2394:2015 and DNV-ST-0126 (2021) §2 (target reliability); DNV-RP-C210 (2021) for probabilistic fatigue inspection planning; DNV-ST-0262 (2016) for lifetime extension.',
          `Every year of extended life spreads the embodied carbon of the foundation (${fmt(steelPerM)} t of steel per metre of pile; ${steelCo2Text(steelPerM * 60).replace(/^.*?≈/, '≈')} for a 60 m pile) over more energy: life extension is the most carbon-efficient investment in an existing farm.`),
        adv('economic', 'info', `Expect about ${fmt(totRepl)} major component replacements per turbine in ${inp.designLife} years (mean time to the first: ${fmt(mt)} years); the probability of none is ${fmt(rEnd * 100)} %.`,
          `Budget heavy-maintenance campaigns for the components with the highest expected counts (${comps.map((c, i) => [c.name, expRepl[i]]).sort((a, b) => b[1] - a[1]).slice(0, 2).map((q) => `${q[0]} ${fmt(q[1])}`).join(', ')}) and align them with the remaining structural life of ${fmt(Math.min(remaining, 1e4))} years; inspection of the fatigue detail is cheaper than the design-fatigue-factor step from 2 to 3 when access exists.`,
          'Risk-based inspection per DNV-RP-G101 / DNV-RP-C210 and whole-life costing per ISO 15686-5: expected failure cost = probability × consequence, traded against inspection, repair and life-extension revenue.',
          energyCo2Text(Math.max(0, Math.min(structLife, 60) - inp.designLife) * 20000)),
      ],
      publish: { 'life.systemYears': Math.min(structLife, 1e4), 'life.designYears': inp.designLife, 'reliability.betaUls': fr.beta, 'reliability.pfFatigue': pfFatEnd },
    };
  },
  verify: [
    { name: 'Lognormal resistance and load: β = ln(R̃/S̃) / √(ln(1 + V_R²) + ln(1 + V_S²))', source: 'Closed-form reliability index of two lognormal variables (Melchers & Beck, Structural Reliability Analysis and Prediction, §1.4.3)',
      // D = 6 m, t = 60 mm: W = 1.64624 m³; mean R = 355 MPa × 1.0 × W = 584.42 MN·m (V = 0.10); mean S = 292.21 MN·m × 1.0 (V = 0.20) → mean ratio 2
      // medians 2/√1.01 and 1/√1.04 → ratio 2.02948, ln = 0.70778; σ = √(0.0099503 + 0.0392207) = 0.221746 → β = 3.1919
      inputs: { D: 6, t: 60, steel: 's355', fyBias: 1, fyCov: 0.1, xrCov: 0, Mk: 292.208, loadDist: 'lognormal', loadBias: 1, loadCov: 0.2, xsCov: 0 }, expect: [{ kpi: 'betaUls', value: 3.1919, tol: 0.003 }] },
    { name: 'Fatigue reliability in the lognormal format: β(t) = ln(T̃/t) / σ_lnT', source: 'Wirsching (1984) J. Struct. Eng. ASCE 110(10): hand calculation',
      // characteristic life 100 yr × 10^(2·s) with s = √ln(1.25)/ln10 = 0.20515 → median 257.2 yr; σ² = ln(1.09) + ln(1.25) + 9·ln(1.0625) = 0.85494 → σ = 0.92463
      // β(25 yr) = ln(257.2/25)/0.92463 = 2.521
      inputs: { fatLife: 100, dff: 1, mSn: 3, snCov: 0.5, deltaCov: 0.3, bCov: 0.25, designLife: 25 }, expect: [{ kpi: 'betaFatEnd', value: 2.521, tol: 0.005 }] },
    { name: 'Paris-law life with constant geometry factor: N = [a_f^(−½) − a_0^(−½)] / [−½·A·(Y·Δσ·√π)³]', source: 'BS 7910:2019 clause 8 simplified law for steel in air (A = 5.21·10⁻¹³ in N–mm units, m = 3); closed-form integration (hand calculation)',
      // A = 5.21e-13 mm/cycle per (N/mm^1.5)³ = 1.6475e-32 m/cycle per (Pa·√m)³ (1 MPa·√m = 31.623 N/mm^1.5)
      // a0 = 0.5 mm, af = 40 mm (half of 80 mm), Y = 1.2, Δσ = 12 MPa: N = (5.0 − 44.721)/(−0.5 × 1.6475e-32 × (2.5523e7)³) = 2.900e8 cycles → 48.33 yr at 6·10⁶ cycles a year
      inputs: { t: 80, a0: 0.5, dsEq: 12, nCyc: 6, yGeom: 1.2, fmEnv: 'air' }, expect: [{ kpi: 'crackLife', value: 48.33, tol: 0.003 }] },
  ],
};
