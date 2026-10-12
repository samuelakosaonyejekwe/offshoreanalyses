// Probabilistic investment appraisal: Latin-hypercube Monte Carlo of the after-tax project cash flow with stated
// distributions and correlation, scenario analysis, tornado, variance decomposition, Sobol indices, value at risk.
import { summarise, valueAtRisk, probabilityBelow, expectedMonetaryValue, convergence, varianceDecomposition, sobol } from '../../engine/econ/risk.js';
import { dist } from '../../engine/math/stats.js';
import { projectInputs, projectParams, dcf, fastFlow, npvDelayed, lcoeOf } from './_project.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, PCT_VALUE, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, tableRows, uniqueNames, threePoint, sampleLHS, percentiles, thinSample, irrFast, irrReport, rankCorr, range } from './_common.js';

const VARS = { price: 'Price', production: 'Production (yield)', capex: 'Capital cost', opex: 'Operating cost' };
const UNC = [
  { variable: 'price', shape: 'pert', low: 78, likely: 100, high: 122 }, { variable: 'production', shape: 'normal', low: 90, likely: 100, high: 110 },
  { variable: 'capex', shape: 'lognormal', low: 95, likely: 100, high: 118 }, { variable: 'opex', shape: 'triangular', low: 88, likely: 100, high: 130 },
];
const SCEN = [
  { name: 'Base case', price: 100, production: 100, capex: 100, opex: 100, delay: 0, life: 0, probability: 40 },
  { name: 'Low price', price: 80, production: 100, capex: 100, opex: 100, delay: 0, life: 0, probability: 15 },
  { name: 'High price', price: 115, production: 100, capex: 100, opex: 100, delay: 0, life: 0, probability: 10 },
  { name: 'Installation overrun and delay', price: 100, production: 100, capex: 118, opex: 100, delay: 1, life: 0, probability: 15 },
  { name: 'Increased maintenance cost', price: 100, production: 97, capex: 100, opex: 130, delay: 0, life: 0, probability: 10 },
  { name: 'Shortened design life', price: 100, production: 100, capex: 100, opex: 105, delay: 0, life: 20, probability: 10 },
];

export default {
  id: 's17-uncertainty-risk', suite: 17, title: 'Uncertainty, sensitivity and investment risk',
  summary: 'Propagates price, yield, cost, availability and schedule uncertainty through the after-tax cash flow to give the distributions of NPV, IRR and LCOE, P10/P50/P90 values, probability of loss, value at risk, the ranking of drivers and a probability-weighted scenario comparison.',
  covers: [7], fidelity: 'numerical',
  standards: ['AACE International RP 57R-09 (integrated cost and schedule risk analysis by Monte Carlo)', 'AACE International RP 41R-08 (range estimating)', 'ISO 31000:2018 and IEC 31010:2019 (risk assessment techniques: Monte Carlo, scenario analysis, sensitivity)', 'Saltelli et al. (2008), Global Sensitivity Analysis: The Primer', 'Rockafellar & Uryasev (2000) (conditional value at risk)'],
  models: ['Probability distribution modelling', 'Normal, lognormal, triangular, PERT and uniform distributions', 'Weibull distribution (start-up delay)', 'Beta distribution (availability loss)', 'Monte Carlo simulation with Latin hypercube sampling', 'Rank correlation by Gaussian copula', 'Scenario analysis', 'One-at-a-time sensitivity and tornado diagram', 'Variance decomposition (standardised regression coefficients)', 'Sobol first-order and total sensitivity indices', 'Value at risk and conditional value at risk', 'Expected monetary value', 'Certainty equivalent with exponential utility', 'Probabilistic investment appraisal'],
  limits: `Uncertainty enters as whole-life multipliers on price, production, capital and operating cost, an availability loss and a start-up delay; year-to-year price volatility is in the forecasting module and construction logic in the schedule-risk module. Cost and delay share one rank correlation; other inputs are independent. Results describe the ranges entered, which are the analyst's judgement. ${CUR_NOTE}`,
  publishes: ['econ.npvMean', 'econ.npvP10', 'econ.npvP50', 'econ.npvP90', 'econ.probLoss', 'econ.valueAtRisk'],
  searchTerms: 'Monte Carlo investment appraisal NPV distribution Sobol sensitivity value at risk',
  mesh: { key: 'nSamples', label: 'Monte Carlo samples', levels: [500, 1000, 2000, 4000], qoi: 'npvMean' },
  inputs: [
    ...projectInputs(),
    { key: 'uncertainty', label: 'Uncertain drivers (as % of the base value)', type: 'table', group: 'Project uncertainty', role: 'ic',
      help: 'Opening uncertainty state. Variable: price, production, capex or opex. Shape: pert, triangular, uniform, normal or lognormal. Low/high are the distribution ends, or the cumulative P10/P90 for normal and lognormal.',
      columns: [{ key: 'variable', label: 'Variable', type: 'text' }, { key: 'shape', label: 'Distribution', type: 'text' }, { key: 'low', label: 'Low', unit: '%', type: 'number' }, { key: 'likely', label: 'Most likely', unit: '%', type: 'number' }, { key: 'high', label: 'High', unit: '%', type: 'number' }], default: UNC },
    { key: 'availLossMax', label: 'Largest plausible availability loss', unit: '%', type: 'number', default: 8, min: 0, max: 60, group: 'Project uncertainty', role: 'ic', help: 'Production lost to failures and access delays beyond the base case, upper bound of a Beta distribution. 0 switches it off.' },
    { key: 'betaA', label: 'Beta shape a (availability loss)', unit: '-', type: 'number', default: 1.5, min: 0.2, max: 20, group: 'Project uncertainty', role: 'ic' },
    { key: 'betaB', label: 'Beta shape b (availability loss)', unit: '-', type: 'number', default: 5, min: 0.2, max: 20, group: 'Project uncertainty', role: 'ic', help: 'b > a skews the loss towards small values with a tail of bad years.' },
    { key: 'delayScale', label: 'Start-up delay, Weibull scale', unit: 'yr', type: 'number', default: 0.45, min: 0, max: 5, group: 'Project uncertainty', role: 'ic', help: 'Characteristic construction delay (63rd percentile). 0 switches delay off.' },
    { key: 'delayShape', label: 'Start-up delay, Weibull shape', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 6, group: 'Project uncertainty', role: 'ic' },
    { key: 'rhoCostDelay', label: 'Rank correlation of capital cost, operating cost and delay', unit: '-', type: 'number', default: 0.4, min: 0, max: 0.95, group: 'Project uncertainty', role: 'ic',
      help: 'Overruns and delays tend to occur together. Price and yield are sampled independently of them.' },
    { key: 'scenarios', label: 'Economic scenarios', type: 'table', group: 'Economic scenarios', role: 'bc', help: 'Discrete alternative futures with probabilities summing to 100 %. Values in % of base; delay in years; life in years (0 = base life).',
      columns: [{ key: 'name', label: 'Scenario', type: 'text' }, { key: 'price', label: 'Price', unit: '%', type: 'number' }, { key: 'production', label: 'Production', unit: '%', type: 'number' }, { key: 'capex', label: 'CAPEX', unit: '%', type: 'number' }, { key: 'opex', label: 'OPEX', unit: '%', type: 'number' }, { key: 'delay', label: 'Delay', unit: 'yr', type: 'number' }, { key: 'life', label: 'Life', unit: 'yr', type: 'number' }, { key: 'probability', label: 'Probability', unit: '%', type: 'number' }], default: SCEN },
    { key: 'confidence', label: 'Confidence level for value at risk', unit: '%', type: 'number', default: 95, min: 50, max: 99.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'maxProbLoss', label: 'Tolerated probability of negative NPV', unit: '%', type: 'number', default: 30, min: 0.1, max: 100, group: 'Acceptance criteria', role: 'criteria', help: 'Investment risk tolerance.' },
    { key: 'riskTolerance', label: 'Corporate risk tolerance', unit: CUR, type: 'number', default: 1.5e9, min: 1e3, max: 1e13, group: 'Acceptance criteria', role: 'bc', help: `Exponential-utility risk tolerance R: the certainty equivalent is −R ln E[exp(−NPV/R)]. Often taken near one sixth of equity. ${PLACEHOLDER}` },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '-', type: 'int', default: 2000, min: 200, max: 20000, group: 'Analysis control', role: 'control' },
    { key: 'nSobol', label: 'Sobol base sample', unit: '-', type: 'int', default: 256, min: 0, max: 4096, group: 'Analysis control', role: 'control', help: 'Model runs = base × (number of variables + 2). 0 skips the Sobol indices.' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control', help: 'Fixed seed: identical inputs reproduce identical results.' },
  ],
  compute(inp) {
    const p = projectParams(inp);
    const unc = tableRows(inp.uncertainty, 'Uncertain drivers', { low: { min: 0 }, likely: { min: 0 }, high: { min: 0 } });
    uniqueNames(unc, 'Uncertain drivers', 'variable');
    const sc = tableRows(inp.scenarios, 'Economic scenarios', { price: { min: 0 }, production: { min: 0 }, capex: { min: 0 }, opex: { min: 0 }, delay: { min: 0, max: 10 }, life: { min: 0, max: 60 }, probability: { min: 0, max: 100 } });
    uniqueNames(sc, 'Economic scenarios');
    const pSum = sc.reduce((s, x) => s + x.probability, 0);
    need(Math.abs(pSum - 100) < 0.01, `Scenario probabilities must sum to 100 % (they sum to ${fmt(pSum, 5)} %). Scenarios are exclusive alternatives, not simultaneous extremes.`);
    const specs = {};
    for (const u of unc) {
      const v = String(u.variable).trim().toLowerCase(); const shape = String(u.shape).trim().toLowerCase();
      need(VARS[v], `Uncertain drivers: variable '${u.variable}' is not one of ${Object.keys(VARS).join(', ')}.`);
      specs[v] = threePoint(shape, u.low / 100, u.likely / 100, u.high / 100, `Uncertain driver '${v}'`);
    }
    const n = Math.round(inp.nSamples);
    const hasDelay = inp.delayScale > 0; const hasAvail = inp.availLossMax > 0;
    // group A (correlated): capex, opex, delay; group B (independent): price, production, availability loss
    const A = {}; const B = {};
    for (const k of ['capex', 'opex']) if (specs[k]) A[k] = specs[k];
    if (hasDelay) A.delay = { dist: 'weibull', scale: inp.delayScale, shape: inp.delayShape };
    for (const k of ['price', 'production']) if (specs[k]) B[k] = specs[k];
    if (hasAvail) B.availLoss = { dist: 'beta', a: inp.betaA, b: inp.betaB, min: 0, max: inp.availLossMax / 100 };
    need(Object.keys(A).length + Object.keys(B).length > 0, 'No uncertain driver is defined: add at least one row or switch on delay or availability loss.');
    const sA = Object.keys(A).length ? run(() => sampleLHS(A, n, inp.seed, Object.keys(A).length > 1 ? inp.rhoCostDelay : 0), 'Sampling') : {};
    const sB = Object.keys(B).length ? run(() => sampleLHS(B, n, inp.seed + 104729, 0), 'Sampling') : {};
    const S = { ...sA, ...sB };
    const names = Object.keys(S);

    const evalNpv = (x) => {
      const f = { price: x.price ?? 1, prod: (x.production ?? 1) * (1 - (x.availLoss ?? 0)), capex: x.capex ?? 1, opex: x.opex ?? 1 };
      return npvDelayed(p, f, x.delay ?? 0);
    };
    const npvS = new Float64Array(n); const irrS = new Float64Array(n); const lcoeS = new Float64Array(n);
    let irrValid = 0;
    const buf = new Float64Array(p.nc + p.nOp + 16); const r0 = {}; const r1 = {};
    const f = { price: 1, prod: 1, capex: 1, opex: 1 };
    for (let k = 0; k < n; k++) {
      f.price = S.price ? S.price[k] : 1; f.prod = (S.production ? S.production[k] : 1) * (1 - (S.availLoss ? S.availLoss[k] : 0));
      f.capex = S.capex ? S.capex[k] : 1; f.opex = S.opex ? S.opex[k] : 1;
      const d = S.delay ? Math.min(S.delay[k], 12) : 0; const d0 = Math.floor(d); const w = d - d0;
      fastFlow(p, f, d0, buf, r0);
      let ir = irrFast(buf.subarray(0, r0.rows)); let npv = r0.npv; let lc = r0.lcoe;
      if (w > 1e-9) { // fractional delay: interpolate between the neighbouring whole-year cases
        fastFlow(p, f, d0 + 1, buf, r1);
        npv += w * (r1.npv - npv); lc += w * (r1.lcoe - lc); ir += w * (irrFast(buf.subarray(0, r1.rows)) - ir);
      }
      npvS[k] = npv; irrS[k] = ir; lcoeS[k] = lc; if (Number.isFinite(ir)) irrValid++;
    }
    const base = dcf(p); const baseIrr = irrReport(base.fcff); const baseLcoe = lcoeOf(p, base).lcoe;
    const st = summarise(npvS, { kind: 'value' });
    const P = percentiles(npvS);
    const irrOk = Float64Array.from(irrS).filter(Number.isFinite); const PI = irrOk.length > 10 ? percentiles(irrOk) : null;
    const PL = percentiles(lcoeS);
    const alpha = inp.confidence / 100;
    const vr = valueAtRisk(npvS, alpha); const vrMean = valueAtRisk(npvS, alpha, { relativeTo: 'mean' });
    const pl = probabilityBelow(npvS, 0);
    const conv = convergence(npvS, 30);
    // certainty equivalent, exponential utility (log-sum-exp for stability)
    const R = inp.riskTolerance; let mx = -Infinity; for (let k = 0; k < n; k++) mx = Math.max(mx, -npvS[k] / R);
    let acc = 0; for (let k = 0; k < n; k++) acc += Math.exp(-npvS[k] / R - mx);
    const ce = -R * (mx + Math.log(acc / n));
    // drivers: one-at-a-time between each input's cumulative P10 and P90, variance shares, rank correlations, Sobol
    const all = { ...A, ...B };
    const q = (nm, u) => dist(all[nm]).icdf(u);
    const neutral = (nm) => (nm === 'delay' || nm === 'availLoss' ? 0 : 1);
    const label = (nm) => VARS[nm] || (nm === 'delay' ? 'Start-up delay' : 'Availability loss');
    const tor = names.map((nm) => { const xb = Object.fromEntries(names.map((m) => [m, neutral(m)])); const lo = evalNpv({ ...xb, [nm]: q(nm, 0.1) }); const hi = evalNpv({ ...xb, [nm]: q(nm, 0.9) }); return { name: label(nm), low: Math.min(lo, hi), high: Math.max(lo, hi) }; })
      .sort((a, b) => (b.high - b.low) - (a.high - a.low));
    const vd = names.length > 1 ? run(() => varianceDecomposition(S, npvS), 'Variance decomposition') : { varianceShare: { [names[0]]: 1 }, r2: 1 };
    const rc = Object.fromEntries(names.map((nm) => [nm, rankCorr(S[nm], npvS)]));
    let sob = null;
    if (inp.nSobol >= 32 && names.length > 1) sob = run(() => sobol({ model: evalNpv, vars: all, n: Math.round(inp.nSobol), seed: inp.seed }), 'Sobol indices');
    const rank = names.map((nm) => ({ nm, share: vd.varianceShare[nm] ?? 0 })).sort((a, b) => b.share - a.share);
    // scenarios
    const scen = sc.map((s) => {
      const pp = s.life > 0 ? projectParams({ ...inp, opLife: Math.round(s.life) }) : p;
      const f = { price: s.price / 100, prod: s.production / 100, capex: s.capex / 100, opex: s.opex / 100 };
      const cf = dcf({ ...pp, delay: Math.round(s.delay) }, f);
      const ir = irrFast(cf.fcff);
      return { ...s, npv: npvDelayed(pp, f, s.delay), irr: Number.isFinite(ir) ? ir : null, lcoe: lcoeOf(pp, cf).lcoe };
    });
    const emv = expectedMonetaryValue(scen.map((s) => ({ p: s.probability / 100, value: s.npv })));

    const terms = p.real ? 'real' : 'nominal';
    const kpis = [
      { key: 'npvBase', label: `Deterministic NPV (${terms}, after tax)`, value: mil(base.npv), unit: MU, explain: 'NPV with every driver at its base value. It is usually above the mean of the distribution because downside ranges are wider than upside ones.' },
      { key: 'npvMean', label: 'Expected NPV', value: mil(st.mean), unit: MU, limit: 0, sense: 'min', explain: 'Mean of the simulated NPV distribution: the risk-neutral value of the project given the uncertainty entered.' },
      { key: 'npvStd', label: 'Standard deviation of NPV', value: mil(st.std), unit: MU, explain: 'Spread of the NPV distribution; an uncertainty indicator to set against the expected value.' },
      { key: 'npvP10', label: 'NPV P10 (cumulative, downside)', value: mil(P(10)), unit: MU, explain: 'NPV not exceeded with 10 % probability: the downside case. In the exceedance convention used for resources this same number is called P90.' },
      { key: 'npvP50', label: 'NPV P50 (median)', value: mil(P(50)), unit: MU, explain: 'NPV with an even chance of being exceeded.' },
      { key: 'npvP90', label: 'NPV P90 (cumulative, upside)', value: mil(P(90)), unit: MU, explain: 'NPV not exceeded with 90 % probability: the upside case (exceedance P10).' },
      { key: 'probLoss', label: 'Probability of negative NPV', value: round(100 * pl.probability), unit: '%', limit: inp.maxProbLoss, sense: 'max', explain: `Share of outcomes that destroy value at the hurdle rate (standard error ${fmt(100 * pl.standardError)} points). Compared with the tolerated probability of loss.` },
      { key: 'var', label: `Value at risk at ${inp.confidence} %`, value: mil(vr.var), unit: MU, explain: `Loss (negative NPV) not exceeded with ${inp.confidence} % confidence. Zero or negative means even this bad outcome still has positive NPV.` },
      { key: 'cvar', label: `Conditional value at risk at ${inp.confidence} %`, value: mil(vr.cvar), unit: MU, explain: `Average loss in the worst ${fmt(100 - inp.confidence)} % of outcomes: how bad things are when they are bad. More informative than VaR for skewed distributions.` },
      { key: 'varMean', label: 'Shortfall from the mean at that confidence', value: mil(vrMean.var), unit: MU, explain: 'How far below the expected NPV the outcome can fall at the chosen confidence: the earnings-at-risk form of VaR.' },
      { key: 'ce', label: 'Risk-adjusted value (certainty equivalent)', value: mil(ce), unit: MU, explain: 'The certain amount a decision maker with the stated risk tolerance would exchange for this uncertain NPV. The gap to the expected NPV is the risk premium.' },
      { key: 'emv', label: 'Expected NPV over the scenarios', value: mil(emv.emv), unit: MU, explain: 'Probability-weighted NPV of the discrete scenarios (expected monetary value); an independent cross-check of the simulation mean.' },
      { key: 'lcoeP50', label: 'LCOE P50', value: round(PL(50)), unit: `${CUR}/MWh`, explain: `Median real pre-tax levelised cost; cumulative P10 ${fmt(PL(10))} and P90 ${fmt(PL(90))} (P90 is the high-cost case). Deterministic value ${fmt(baseLcoe)}.` },
      { key: 'mcError', label: 'Monte Carlo standard error of the mean', value: mil(st.standardError), unit: MU, explain: 'Sampling error of the expected NPV at this sample size. Differences between options smaller than about twice this are not resolved.' },
    ];
    if (PI) kpis.splice(7, 0, { key: 'irrP50', label: 'IRR P50', value: round(100 * PI(50)), unit: '%/yr', explain: `Median economic IRR in ${terms} terms; cumulative P10 ${pct(PI(10))} and P90 ${pct(PI(90))}. Defined for ${pct(irrValid / n)} of samples; where no sign change exists IRR is absent and NPV governs.` });
    const checks = [
      check('Probability of loss ≤ tolerance', 100 * pl.probability, inp.maxProbLoss, '%', 'Investment risk tolerance (corporate criterion)', 'Guards against sanctioning a project with a higher chance of destroying value than the investor accepts; above 1.0 the risk must be reduced or transferred.'),
      check('Conditional value at risk ≤ corporate risk tolerance', Math.max(vr.cvar, 0), inp.riskTolerance, CUR, 'Rockafellar & Uryasev (2000); corporate risk appetite', 'Guards against a tail loss larger than the organisation can absorb. Above 1.0 the average loss in the worst cases exceeds the stated risk tolerance and the exposure should be shared, insured or reduced.'),
    ];
    const cdfN = 50;
    const plots = [
      { type: 'hist', title: 'Distribution of NPV', values: thinSample(npvS, 1500).map(mil), bins: 40, x: { label: 'NPV', unit: MU }, y: { label: 'Frequency' }, explain: 'Simulated NPV outcomes. The area left of zero is the probability of loss; a long left tail signals delay and overrun risk.' },
      { type: 'line', title: 'Cumulative probability of NPV', x: { label: 'NPV', unit: MU }, y: { label: 'Probability not exceeded', unit: '%' },
        series: [{ name: 'NPV', x: range(cdfN + 1, (k) => mil(P(100 * k / cdfN))), y: range(cdfN + 1, (k) => round(100 * k / cdfN)) }, { name: 'Break-even', x: [0, 0], y: [0, 100], dash: true }, { name: 'Deterministic', x: [mil(base.npv), mil(base.npv)], y: [0, 100], dash: true }],
        explain: 'Read the probability of loss where the curve crosses zero, and P10/P50/P90 at 10, 50 and 90 % (cumulative convention). The deterministic case is rarely the median.' },
      { type: 'tornado', title: 'NPV with each driver at its P10 and P90', base: mil(evalNpv(Object.fromEntries(names.map((m) => [m, neutral(m)])))), items: tor.map((t) => ({ name: t.name, low: mil(t.low), high: mil(t.high) })), x: { label: 'NPV', unit: MU },
        explain: 'One driver at a time moved between its cumulative P10 and P90 with the others at base. Bars are ordered by swing: the top ones deserve the engineering and commercial effort.' },
      { type: 'bar', title: 'Share of NPV variance by driver', categories: rank.map((r) => label(r.nm)), x: { label: 'Driver' }, y: { label: 'Share of variance', unit: '%' },
        series: [{ name: 'Regression (SRC²)', y: rank.map((r) => round(100 * r.share, 4)) }, ...(sob ? [{ name: 'Sobol first order', y: rank.map((r) => round(100 * Math.max(sob.first[r.nm], 0), 4)) }, { name: 'Sobol total effect', y: rank.map((r) => round(100 * Math.max(sob.total[r.nm], 0), 4)) }] : [])],
        explain: 'How much of the NPV variance each driver explains. Where the Sobol total effect exceeds the first-order index, the driver acts through interactions with others (for example delay with price).' },
      { type: 'bar', title: 'NPV by scenario', categories: [...scen.map((s) => s.name), 'Expected'], x: { label: 'Scenario' }, y: { label: 'NPV', unit: MU }, series: [{ name: 'NPV', y: [...scen.map((s) => mil(s.npv)), mil(emv.emv)] }],
        explain: 'NPV of each discrete scenario and their probability-weighted mean. Scenarios show coherent stories that a distribution hides.' },
      { type: 'line', title: 'Monte Carlo convergence', x: { label: 'Samples', unit: '-' }, y: { label: 'Running mean NPV', unit: MU }, series: [{ name: 'Running mean', x: conv.n, y: conv.runningMean.map(mil) }], explain: 'Running mean of NPV against sample count; the tail should be flat to within the standard error quoted.' },
    ];
    const tables = [
      { title: 'Economic percentiles (cumulative convention)', columns: [{ key: 'm', label: 'Measure' }, { key: 'p10', label: 'P10' }, { key: 'p50', label: 'P50' }, { key: 'p90', label: 'P90' }, { key: 'mean', label: 'Mean' }, { key: 'unit', label: 'Unit' }],
        rows: [{ m: 'NPV', p10: mil(P(10)), p50: mil(P(50)), p90: mil(P(90)), mean: mil(st.mean), unit: MU }, ...(PI ? [{ m: 'IRR', p10: round(100 * PI(10), 4), p50: round(100 * PI(50), 4), p90: round(100 * PI(90), 4), mean: round(100 * irrOk.reduce((a, c) => a + c, 0) / irrOk.length, 4), unit: '%/yr' }] : []),
          { m: 'LCOE', p10: round(PL(10), 4), p50: round(PL(50), 4), p90: round(PL(90), 4), mean: round(Array.from(lcoeS).reduce((a, c) => a + c, 0) / n, 4), unit: `${CUR}/MWh` }],
        explain: 'P10, P50 and P90 as cumulative (non-exceedance) percentiles for every measure: P10 is the low value, P90 the high value. For NPV and IRR low is bad; for LCOE high is bad.' },
      { title: 'Driver ranking', columns: [{ key: 'name', label: 'Driver' }, { key: 'dist', label: 'Distribution' }, { key: 'p10', label: 'Input P10' }, { key: 'p90', label: 'Input P90' }, { key: 'rc', label: 'Rank correlation with NPV' }, { key: 'src', label: 'Variance share', unit: '%' }, { key: 's1', label: 'Sobol first', unit: '%' }, { key: 'st', label: 'Sobol total', unit: '%' }],
        rows: rank.map((r) => ({ name: label(r.nm), dist: all[r.nm].dist, p10: round(q(r.nm, 0.1), 4), p90: round(q(r.nm, 0.9), 4), rc: round(rc[r.nm], 3), src: round(100 * r.share, 3), s1: sob ? round(100 * sob.first[r.nm], 3) : '–', st: sob ? round(100 * sob.total[r.nm], 3) : '–' })),
        explain: 'Each uncertain input with its distribution, its P10–P90 range and three measures of its influence on NPV. Inputs are multipliers on the base value, except delay (years) and availability loss (fraction).' },
      { title: 'Scenario analysis', columns: [{ key: 'name', label: 'Scenario' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'npv', label: 'NPV', unit: MU }, { key: 'irr', label: 'IRR', unit: '%/yr' }, { key: 'lcoe', label: 'LCOE', unit: `${CUR}/MWh` }, { key: 'd', label: 'NPV change from base', unit: MU }],
        rows: scen.map((s) => ({ name: s.name, p: s.probability, npv: mil(s.npv), irr: s.irr != null ? round(100 * s.irr, 4) : '–', lcoe: round(s.lcoe, 4), d: mil(s.npv - base.npv) })), explain: 'Each named scenario with its probability and results; the last column is the value at stake if that future occurs.' },
    ];
    const top = tor[0];
    const riskPrem = st.mean - ce;
    const advice = [
      adv('economic', pl.probability * 100 <= inp.maxProbLoss && st.mean > 0 ? 'ok' : 'warn', `Expected NPV ${fM(st.mean)} against ${fM(base.npv)} deterministic; cumulative P10 / P50 / P90 = ${fM(P(10))} / ${fM(P(50))} / ${fM(P(90))}; probability of loss ${pct(pl.probability)}; CVaR at ${inp.confidence} % ${fM(vr.cvar)}; certainty equivalent ${fM(ce)} (risk premium ${fM(riskPrem)}).`,
        pl.probability * 100 > inp.maxProbLoss ? `The chance of loss exceeds the ${inp.maxProbLoss} % tolerance. Reduce the top driver's range (contracts, hedges, fixed-price EPC) or re-scope until the P10 NPV of ${fM(P(10))} is acceptable.` : `Risk is inside tolerance. Report the P10 beside the P50 and hold the gap between deterministic and expected NPV (${fM(base.npv - st.mean)}) as the expected cost of risk.`,
        'Probabilistic investment appraisal, AACE RP 57R-09 and IEC 31010:2019 §B.5 (Monte Carlo); CVaR per Rockafellar & Uryasev (2000); percentiles stated as cumulative non-exceedance.'),
      adv('engineering', 'info', `"${top.name}" swings NPV by ${fM(top.high - top.low)} between its P10 and P90${rank.length > 1 ? ` and "${label(rank[0].nm)}" explains ${pct(rank[0].share)} of the variance (linear model R² ${fmt(vd.r2)})` : ''}.`,
        `Reduce that range at its technical source: yield uncertainty with longer site measurements and wake validation, capital cost with firm quantities and vessel-day estimates, availability with the reliability analysis, delay with weather-window and float analysis. A one-year delay alone changes NPV by ${fM(npvDelayed(p, null, 1) - base.npv)}.`,
        'Global sensitivity analysis, Saltelli et al. (2008); ISO 31000:2018 §6.4 (risk analysis feeding treatment).',
        `Carbon-price and embodied-carbon exposure belong in the ranges: a 60 000 t steel scope embodies about ${fmt(60000 * CO2_STEEL / 1e3)} kt CO₂e (screening factor ${CO2_STEEL} t CO₂e/t), so a 50 ${CUR}/t swing in a carbon border price moves capital cost by about ${fM(60000 * CO2_STEEL * 50)}. Include it as a cost driver where it applies.`),
    ];
    if (Math.abs(emv.emv - st.mean) > 3 * st.standardError + 0.15 * st.std) advice.push(adv('economic', 'info', `The scenario expectation (${fM(emv.emv)}) and the simulation mean (${fM(st.mean)}) differ by ${fM(Math.abs(emv.emv - st.mean))}.`, 'Reconcile the scenario probabilities with the distributions: both describe the same uncertainty and should give similar expectations.', 'Scenario analysis versus simulation, IEC 31010:2019 §B.2 and §B.5.'));
    return { kpis, checks, plots, tables, advice,
      notes: [PCT_VALUE, `${terms === 'real' ? 'Real' : 'Nominal'} terms at ${pct(p.rate)}; after tax. ${n} Latin-hypercube samples, seed ${inp.seed}.`, `Correlation: capital cost, operating cost and delay share a rank correlation of ${inp.rhoCostDelay}; price, yield and availability loss are independent. Extremes are therefore not assumed to coincide.`, baseIrr.text],
      publish: { 'econ.npvMean': st.mean, 'econ.npvP10': P(10), 'econ.npvP50': P(50), 'econ.npvP90': P(90), 'econ.probLoss': pl.probability, 'econ.valueAtRisk': vr.var } };
  },
  verify: [
    // only price uncertain, no tax or inflation: NPV = PV(revenue) f − PV(opex) − capex with f ~ N(1, σ), σ = 0.1/1.2816 = 0.078028
    // PV(revenue) = 200 × (P/A, 8 %, 20) = 1963.629 M; NPV(1) = −1000 + 150 × 9.818147 = 472.722 M; std = 1963.629 × 0.078028 = 153.218 M;
    // P10 = 472.722 − 1.2816 × 153.218 = 276.36 M; P(NPV < 0) = Φ(−3.0853) = 0.10 %
    { name: 'Linear propagation of a normal price factor', source: 'Closed form: NPV is linear in the price factor, so mean = NPV(1), σ = PV(revenue) σ_f, P10 = mean − 1.2816 σ',
      inputs: { basis: 'real', inflation: 0, sunk: 0, aep: 2000, degradation: 0, constructionYears: 1, opLife: 20, capex: 1e9, opexFixed: 5e7, opexVar: 0, decom: 0, price: 100, taxRate: 0, royalty: 0, discountRate: 8,
        uncertainty: [{ variable: 'price', shape: 'normal', low: 90, likely: 100, high: 110 }], availLossMax: 0, delayScale: 0, nSamples: 4000, nSobol: 0 },
      expect: [{ kpi: 'npvBase', value: 472.722, tol: 1e-5 }, { kpi: 'npvMean', value: 472.722, tol: 0.002 }, { kpi: 'npvStd', value: 153.218, tol: 0.01 }, { kpi: 'npvP10', value: 276.36, tol: 0.01 }, { kpi: 'probLoss', value: 0.10, abs: 0.1 }] },
    // scenario expectation: 50 % × NPV(price 100 %) + 50 % × NPV(price 90 %) = 472.722 − 0.5 × 0.1 × 1963.629 = 374.541 M
    { name: 'Expected monetary value of two price scenarios', source: 'Hand calculation: 0.5 × 472.722 + 0.5 × (472.722 − 196.363) = 374.541 M',
      inputs: { basis: 'real', inflation: 0, sunk: 0, aep: 2000, degradation: 0, constructionYears: 1, opLife: 20, capex: 1e9, opexFixed: 5e7, opexVar: 0, decom: 0, price: 100, taxRate: 0, royalty: 0, discountRate: 8, nSamples: 500, nSobol: 0,
        scenarios: [{ name: 'Base', price: 100, production: 100, capex: 100, opex: 100, delay: 0, life: 0, probability: 50 }, { name: 'Low', price: 90, production: 100, capex: 100, opex: 100, delay: 0, life: 0, probability: 50 }] },
      expect: [{ kpi: 'emv', value: 374.541, tol: 1e-5 }] },
  ],
};
