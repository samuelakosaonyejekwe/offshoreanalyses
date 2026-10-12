// Probabilistic cost estimate: range estimating of the cost lines (Latin hypercube with rank correlation), a
// discrete risk register, contingency at a chosen confidence level and Bayesian updating with received quotations.
import { AACE_CLASSES } from '../../engine/econ/cost.js';
import { riskRegister, summarise, convergence } from '../../engine/econ/risk.js';
import { bayesNormalUpdate, priorFromEstimateClass } from '../../engine/econ/forecast.js';
import { percentileSorted, dist } from '../../engine/math/stats.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, PCT_COST, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, tableRows, uniqueNames, threePoint, SHAPES, thinSample, range, sampleLHS } from './_common.js';

const LINES = [
  { name: 'Turbines and equipment', low: 5.4e8, likely: 5.8e8, high: 6.9e8 }, { name: 'Foundations (steel, fabrication)', low: 2.0e8, likely: 2.4e8, high: 3.3e8 },
  { name: 'Cables (supply)', low: 1.0e8, likely: 1.15e8, high: 1.5e8 }, { name: 'Electrical substations', low: 1.4e8, likely: 1.55e8, high: 2.0e8 },
  { name: 'Marine installation', low: 1.5e8, likely: 2.0e8, high: 3.4e8 }, { name: 'Engineering, management, owner', low: 1.9e8, likely: 2.2e8, high: 2.8e8 },
];
const RISKS = [
  { name: 'Unexpected soil conditions at foundations', probability: 15, low: 1.0e7, likely: 2.5e7, high: 7.0e7 },
  { name: 'Cable damage during lay or burial', probability: 20, low: 5.0e6, likely: 1.5e7, high: 4.5e7 },
  { name: 'Late delivery of installation vessel', probability: 25, low: 8.0e6, likely: 2.0e7, high: 6.0e7 },
  { name: 'Grid connection delay', probability: 10, low: 1.0e7, likely: 3.0e7, high: 9.0e7 },
  { name: 'Favourable steel price at award (opportunity)', probability: 30, low: -3.0e7, likely: -1.5e7, high: -5.0e6 },
];

export default {
  id: 's17-cost-risk', suite: 17, title: 'Probabilistic cost estimate and contingency',
  summary: 'Combines three-point cost ranges, their correlation and a risk register into the distribution of total CAPEX, the contingency for a chosen confidence level and a Bayesian update from received quotations.',
  covers: [2, 7], fidelity: 'numerical',
  standards: ['AACE International RP 41R-08 (risk analysis and contingency determination using range estimating)', 'AACE International RP 57R-09 (integrated cost and schedule risk analysis)', 'AACE International RP 18R-97 (estimate classes)', 'ISO 31000:2018 (risk management)'],
  models: ['Probabilistic cost estimating', 'Triangular, PERT, uniform, normal and lognormal distributions', 'Latin hypercube sampling with rank correlation (one-factor Gaussian copula)', 'Risk register (probability × impact) simulation', 'Expected monetary value', 'Bayesian updating of cost forecasts (conjugate normal)', 'Monte Carlo convergence diagnostics'],
  limits: `Cost lines share one rank correlation; risks are independent of each other and of the line ranges. Schedule-driven cost is represented only through the register (see the schedule-risk module). Results describe the uncertainty the user has entered, not the market. ${CUR_NOTE}`,
  publishes: ['econ.capexP50', 'econ.capexP80', 'econ.capexP90', 'econ.contingency'],
  searchTerms: 'range estimating contingency Monte Carlo cost risk analysis P80',
  mesh: { key: 'nSamples', label: 'Monte Carlo samples', levels: [1000, 2000, 4000, 8000], qoi: 'mean' },
  inputs: [
    { key: 'aaceClass', label: 'Estimate class (AACE 18R-97)', type: 'select', default: '4', group: 'Opening uncertainty state', role: 'ic', options: [5, 4, 3, 2, 1].map((c) => ({ value: String(c), label: `Class ${c} — ${AACE_CLASSES[c].usage}` })),
      help: 'Used for the prior of the Bayesian update and as a plausibility check on the entered ranges.' },
    { key: 'lines', label: 'Cost lines (three-point ranges)', type: 'table', group: 'Cost estimate ranges', role: 'cost',
      help: `Low and high are the ends of the chosen distribution (or the cumulative P10 and P90 for normal and lognormal). ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Cost line', type: 'text' }, { key: 'low', label: 'Low', unit: CUR, type: 'number' }, { key: 'likely', label: 'Most likely', unit: CUR, type: 'number' }, { key: 'high', label: 'High', unit: CUR, type: 'number' }], default: LINES },
    { key: 'shape', label: 'Distribution of the cost lines', type: 'select', default: 'pert', group: 'Cost estimate ranges', role: 'control', options: SHAPES },
    { key: 'rho', label: 'Rank correlation between cost lines', unit: '-', type: 'number', default: 0.3, min: 0, max: 0.95, group: 'Cost estimate ranges', role: 'control',
      help: 'Common market drivers make lines move together. Zero understates the spread of the total; one would make every extreme coincide.' },
    { key: 'risks', label: 'Risk register', type: 'table', group: 'Failure and risk costs', role: 'cost',
      help: `Discrete events: probability of occurrence and three-point impact if they occur (negative for opportunities). ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Risk', type: 'text' }, { key: 'probability', label: 'Probability', unit: '%', type: 'number' }, { key: 'low', label: 'Impact low', unit: CUR, type: 'number' }, { key: 'likely', label: 'Impact likely', unit: CUR, type: 'number' }, { key: 'high', label: 'Impact high', unit: CUR, type: 'number' }], default: RISKS },
    { key: 'quotes', label: 'Received quotations for the whole scope', type: 'table', group: 'Bayesian update', role: 'cost',
      help: `Independent contractor or supplier totals on the same scope and price base. Leave empty to skip the update. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'amount', label: 'Quoted total', unit: CUR, type: 'number' }],
      default: [{ name: 'Bidder 1', amount: 1.69e9 }, { name: 'Bidder 2', amount: 1.78e9 }, { name: 'Bidder 3', amount: 1.61e9 }] },
    { key: 'quoteCov', label: 'Scatter of a single quotation', unit: '%', type: 'number', default: 8, min: 0.5, max: 60, group: 'Bayesian update', role: 'control', help: 'Standard deviation of one quotation about the true cost, as a percentage.' },
    { key: 'confidence', label: 'Confidence level for the budget', unit: '%', type: 'number', default: 80, min: 50, max: 99, group: 'Acceptance criteria', role: 'criteria', help: 'Cumulative percentile of total cost the budget should cover (P80 is common for sanction).' },
    { key: 'budgetCap', label: 'Available capital budget', unit: CUR, type: 'number', default: 1.9e9, min: 1, max: 1e12, group: 'Financial constraints', role: 'bc' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '-', type: 'int', default: 4000, min: 200, max: 50000, group: 'Analysis control', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Analysis control', role: 'control', help: 'Fixed seed: the same inputs always give the same result.' },
  ],
  compute(inp) {
    const lines = tableRows(inp.lines, 'Cost lines', { low: {}, likely: {}, high: {} });
    uniqueNames(lines, 'Cost lines');
    const risks = tableRows(inp.risks, 'Risk register', { probability: { min: 0, max: 100 }, low: {}, likely: {}, high: {} }, { min: 0 });
    uniqueNames(risks, 'Risk register');
    const quotes = tableRows(inp.quotes, 'Quotations', { amount: { min: 1e-9 } }, { min: 0 });
    const n = Math.round(inp.nSamples);
    const items = lines.map((l) => ({ name: l.name, dist: threePoint(inp.shape, l.low, l.likely, l.high, `Cost line '${l.name}'`) }));
    const base = lines.reduce((s, l) => s + l.likely, 0);
    need(base > 0, 'The most-likely values of the cost lines must sum to a positive base estimate.');
    // cost lines: Latin hypercube with one common rank correlation (one-factor Gaussian copula)
    const cols = run(() => sampleLHS(Object.fromEntries(items.map((i) => [i.name, i.dist])), n, inp.seed, lines.length > 1 ? inp.rho : 0), 'Range estimate');
    const lineTotal = new Float64Array(n);
    for (const i of items) { const c = cols[i.name]; for (let k = 0; k < n; k++) lineTotal[k] += c[k]; }
    // the register uses an independent stream so that events are not tied to the line ranges
    const rr = risks.length ? run(() => riskRegister({ risks: risks.map((r) => ({ name: r.name, probability: r.probability / 100, impact: r.low === r.high ? r.likely : threePoint('triangular', r.low, r.likely, r.high, `Risk '${r.name}'`) })), n, seed: inp.seed + 7919 }), 'Risk register') : null;
    const total = new Float64Array(n);
    for (let k = 0; k < n; k++) total[k] = lineTotal[k] + (rr ? rr.samples[k] : 0);
    const st = summarise(total, { kind: 'cost' });
    const sorted = Float64Array.from(total).sort();
    const P = (q) => percentileSorted(sorted, q);
    const pConf = P(inp.confidence);
    const cont = pConf - base;
    let over = 0; for (let k = 0; k < n; k++) if (total[k] > inp.budgetCap) over++;
    const pOver = over / n;
    const conv = convergence(total, 30);
    // contribution of each line and of the register to the spread of the total
    // analytic P10/P90 of each line for the tornado
    const drivers = items.map((it, j) => { const d = dist(it.dist); return { name: it.name, low: d.icdf(0.1), high: d.icdf(0.9), likely: lines[j].likely }; });
    const cls = Number(inp.aaceClass);
    let post = null;
    if (quotes.length) {
      const prior = priorFromEstimateClass({ estimate: base, cls });
      post = run(() => bayesNormalUpdate({ priorMean: prior.priorMean, priorStd: prior.priorStd, observations: quotes.map((q) => q.amount), obsStd: inp.quoteCov / 100 * base }), 'Bayesian update');
      post.prior = prior;
    }

    const kpis = [
      { key: 'base', label: 'Base estimate (sum of most-likely values)', value: mil(base), unit: MU, explain: 'Deterministic estimate before contingency. Because cost ranges are skewed to the high side, this figure usually sits below the P50 of the distribution.' },
      { key: 'mean', label: 'Expected total cost', value: mil(st.mean), unit: MU, explain: 'Mean of the simulated total including the risk register. The difference from the base estimate is the expected cost of the uncertainty that has been described.' },
      { key: 'p10', label: 'P10 total cost (cumulative)', value: mil(P(10)), unit: MU, explain: 'Cost not exceeded with 10 % probability: the optimistic end of the range.' },
      { key: 'p50', label: 'P50 total cost (cumulative)', value: mil(P(50)), unit: MU, explain: 'Cost not exceeded with 50 % probability: an even chance of over- or under-running. A budget at P50 overruns half the time.' },
      { key: 'p80', label: 'P80 total cost (cumulative)', value: mil(P(80)), unit: MU, explain: 'Cost not exceeded with 80 % probability: a common sanction budget level.' },
      { key: 'p90', label: 'P90 total cost (cumulative)', value: mil(P(90)), unit: MU, explain: 'Cost not exceeded with 90 % probability: the conservative, high-cost case used for funding headroom and lender stress cases.' },
      { key: 'budget', label: `Budget at P${inp.confidence}`, value: mil(pConf), unit: MU, limit: mil(inp.budgetCap), sense: 'max', explain: `Total cost at the chosen cumulative confidence level; the probability of exceeding it is ${100 - inp.confidence} %.` },
      { key: 'contingency', label: `Contingency to reach P${inp.confidence}`, value: mil(cont), unit: MU, explain: 'Amount to add to the base estimate to reach the chosen confidence. It is expected to be spent in part; it is not a reserve for scope change.' },
      { key: 'contPct', label: 'Contingency as share of base', value: round(100 * cont / base), unit: '%', explain: 'Risk-derived contingency percentage for comparison with the flat allowance used in the deterministic estimate.' },
      { key: 'riskEV', label: 'Expected value of the risk register', value: mil(rr ? rr.expectedValue : 0), unit: MU, explain: 'Sum of probability × mean impact over the register (expected monetary value). Opportunities reduce it.' },
      { key: 'pOver', label: 'Probability of exceeding the budget', value: round(100 * pOver), unit: '%', limit: 100 - inp.confidence, sense: 'max', explain: 'Share of simulated outcomes above the available capital. It should not exceed the overrun probability implied by the chosen confidence level.' },
      { key: 'mcError', label: 'Monte Carlo standard error of the mean', value: round(100 * conv.relativeError, 3), unit: '%', explain: 'Sampling error of the expected cost with the chosen number of samples; it should be well below the differences being decided on.' },
    ];
    if (post) kpis.push(
      { key: 'postMean', label: 'Posterior expected cost after quotations', value: mil(post.mean), unit: MU, explain: `Estimate updated with ${quotes.length} quotation(s): a precision-weighted average of the class-${cls} prior and the quotations. The prior carries ${pct(post.priorWeight)} of the weight.` },
      { key: 'postP90', label: 'Posterior P90 (cumulative)', value: mil(post.p90), unit: MU, explain: 'Cost not exceeded with 90 % probability after the update. The narrowing against the prior range is the value of the market information received.' });

    const checks = [
      check(`P${inp.confidence} cost within available budget`, pConf, inp.budgetCap, CUR, 'AACE RP 41R-08 (contingency at a stated confidence level)', 'Guards against a budget that will be overrun more often than the decision maker accepts; above 1.0 the funding does not cover the chosen confidence.'),
      check('Base estimate within available budget', base, inp.budgetCap, CUR, 'Financial constraint (affordability limit)', 'Guards against a project that is unaffordable even before any allowance for uncertainty.'),
    ];
    const cdfN = 60;
    const plots = [
      { type: 'hist', title: 'Distribution of total cost', values: thinSample(total, 1500).map(mil), bins: 40, x: { label: 'Total cost', unit: MU }, y: { label: 'Frequency' },
        explain: 'Simulated outcomes of total CAPEX. The long right tail is what contingency is for; compare where the budget falls in the distribution.' },
      { type: 'line', title: 'Cumulative probability (S-curve)', x: { label: 'Total cost', unit: MU }, y: { label: 'Probability not exceeded', unit: '%' },
        series: [{ name: 'Total cost', x: range(cdfN + 1, (k) => mil(P(100 * k / cdfN))), y: range(cdfN + 1, (k) => round(100 * k / cdfN)) }, { name: 'Base estimate', x: [mil(base), mil(base)], y: [0, 100], dash: true },
          { name: `P${inp.confidence}`, x: [mil(pConf), mil(pConf)], y: [0, 100], dash: true }, { name: 'Budget', x: [mil(inp.budgetCap), mil(inp.budgetCap)], y: [0, 100], dash: true }],
        explain: 'Read across from a confidence level to the cost that delivers it, or up from the budget to the probability of staying inside it. The height at which the base estimate cuts the curve is its true confidence.' },
      { type: 'tornado', title: 'Cost lines: P10 to P90 range', base: mil(base), items: drivers.map((d) => ({ name: d.name, low: mil(base - d.likely + d.low), high: mil(base - d.likely + d.high) })).sort((a, b) => (b.high - b.low) - (a.high - a.low)),
        x: { label: 'Total cost with one line at its P10 / P90', unit: MU }, explain: 'Each line moved from its cumulative P10 to its P90 with the others at their most-likely value. The widest bars are the lines whose definition or quotation would reduce uncertainty most.' },
      { type: 'line', title: 'Monte Carlo convergence', x: { label: 'Samples', unit: '-' }, y: { label: 'Running mean', unit: MU }, series: [{ name: 'Running mean', x: conv.n, y: conv.runningMean.map(mil) }],
        explain: 'Running mean of the total against sample count. A flat tail shows the sample is large enough; a drifting tail calls for more samples.' },
    ];
    const tables = [
      { title: 'Percentiles of total cost', columns: [{ key: 'p', label: 'Cumulative percentile' }, { key: 'cost', label: 'Total cost', unit: MU }, { key: 'cont', label: 'Contingency on base', unit: MU }, { key: 'pc', label: 'Contingency', unit: '%' }],
        rows: [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95].map((q) => ({ p: `P${q}`, cost: mil(P(q)), cont: mil(P(q) - base), pc: round(100 * (P(q) - base) / base, 3) })),
        explain: 'Total cost and the contingency it implies at each cumulative confidence level (Pxx = not exceeded with xx % probability).' },
      { title: 'Risk register ranking', columns: [{ key: 'name', label: 'Risk' }, { key: 'ev', label: 'Expected value', unit: MU }, { key: 'share', label: 'Share of register', unit: '%' }],
        rows: (rr ? rr.byRisk : []).map((b) => ({ name: b.name, ev: mil(b.expected), share: round(100 * b.share, 3) })),
        explain: 'Risks ranked by probability × mean impact. Mitigation spend is justified up to the expected value it removes.' },
    ];
    if (post) tables.push({ title: 'Bayesian update with quotations', columns: [{ key: 'item', label: 'Item' }, { key: 'mean', label: 'Mean', unit: MU }, { key: 'sd', label: 'Standard deviation', unit: MU }],
      rows: [{ item: `Prior from Class ${cls} band`, mean: mil(post.prior.priorMean), sd: mil(post.prior.priorStd) }, { item: `Quotations (n = ${quotes.length})`, mean: mil(quotes.reduce((s, q) => s + q.amount, 0) / quotes.length), sd: mil(inp.quoteCov / 100 * base) }, { item: 'Posterior', mean: mil(post.mean), sd: mil(post.std) }],
      explain: 'How the estimate moves and tightens when quotations arrive: the posterior lies between prior and quotations, closer to whichever is more precise.' });

    const widest = plots[2].items[0];
    const baseConf = 100 * (sorted.findIndex((v) => v >= base) < 0 ? 1 : sorted.findIndex((v) => v >= base) / n);
    const advice = [
      adv('economic', pConf <= inp.budgetCap ? 'ok' : 'fail', `Base ${fM(base)}; P50 ${fM(P(50))}, P80 ${fM(P(80))}, P90 ${fM(P(90))} (cumulative). The base estimate has only a ${fmt(baseConf)} % chance of not being exceeded. Contingency to P${inp.confidence} is ${fM(cont)} (${pct(cont / base)} of base).`,
        pConf <= inp.budgetCap ? `Set the budget at ${fM(pConf)} and hold ${fM(P(90) - pConf)} as management reserve to P90; release contingency as lines are contracted.` : `The P${inp.confidence} cost exceeds the ${fM(inp.budgetCap)} available by ${fM(pConf - inp.budgetCap)} (overrun probability ${pct(pOver)}): add funding, cut scope or mitigate the top risks before sanction.`,
        'AACE RP 41R-08 (range estimating, contingency at a confidence level) and RP 57R-09; percentiles reported as cumulative non-exceedance.'),
      adv('engineering', 'info', `The widest line is "${widest.name}" (${fM((widest.high - widest.low) * 1e6)} between its P10 and P90)${rr ? `; the largest register item is "${rr.byRisk[0].name}" with an expected value of ${fM(rr.byRisk[0].expected)}` : ''}.`,
        'Narrow that range with engineering definition: firm quantities from the structural and installation analyses, a site investigation for ground risk, or a weather-window study for marine durations. Up to the expected value removed is worth spending on it.',
        'Value of information and risk treatment, ISO 31000:2018 §6.5; AACE RP 18R-97 (accuracy improves with definition).',
        `Cost uncertainty in steel-heavy lines is largely quantity uncertainty: a 5 % tonnage growth on a 60 000 t scope is 3000 t, about ${fmt(3000 * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e/t, so freezing quantities early limits both the cost tail and the embodied carbon.`),
    ];
    if (post) advice.push(adv('economic', 'info', `Quotations average ${fM(quotes.reduce((s, q) => s + q.amount, 0) / quotes.length)}; the updated expected cost is ${fM(post.mean)} with a standard deviation of ${fM(post.std)} (prior ${fM(post.prior.priorStd)}).`,
      Math.abs(post.mean - st.mean) > 2 * post.std ? 'The market evidence differs from the simulated estimate by more than two posterior standard deviations: re-base the line ranges on the quotations.' : 'The quotations are consistent with the estimate; narrow the line ranges accordingly and re-run to release contingency.',
      'Bayesian updating of cost forecasts; AACE RP 18R-97 class bands read as an 80 % interval.'));
    return { kpis, checks, plots, tables, advice,
      notes: [PCT_COST, `${n} Latin-hypercube samples, seed ${inp.seed}; cost lines share a rank correlation of ${inp.rho}, risks are independent events.`, 'Ranges and probabilities are the estimator’s judgement entered above; no default represents market data.'],
      publish: { 'econ.capexP50': P(50), 'econ.capexP80': P(80), 'econ.capexP90': P(90), 'econ.contingency': cont } };
  },
  verify: [
    // Triangular(80, 100, 150): mean = (80 + 100 + 150)/3 = 110; F(mode) = 20²/(70·20) = 0.2857 < 0.5 → P50 = 150 − sqrt(0.5·70·50) = 108.167;
    // P90 = 150 − sqrt(0.1·70·50) = 131.292
    { name: 'Single triangular line (80, 100, 150)', source: 'Closed-form triangular distribution: mean 110, P50 = 150 − √1750 = 108.167, P90 = 150 − √350 = 131.292',
      inputs: { lines: [{ name: 'A', low: 80e6, likely: 100e6, high: 150e6 }], shape: 'triangular', risks: [], quotes: [], rho: 0, nSamples: 8000 },
      expect: [{ kpi: 'mean', value: 110, tol: 0.002 }, { kpi: 'p50', value: 108.167, tol: 0.003 }, { kpi: 'p90', value: 131.292, tol: 0.003 }, { kpi: 'base', value: 100, tol: 1e-9 }] },
    // register only adds its expected value to the mean: 0.25 × 40 = 10 on a fixed base of 100 → mean 110 (uniform(99.999, 100.001) line)
    { name: 'Risk register expected value: 25 % × 40', source: 'Expected monetary value Σ p × impact = 0.25 × 40 = 10',
      inputs: { lines: [{ name: 'A', low: 99.999e6, likely: 100e6, high: 100.001e6 }], shape: 'uniform', risks: [{ name: 'R', probability: 25, low: 40e6, likely: 40e6, high: 40e6 }], quotes: [], nSamples: 8000 },
      expect: [{ kpi: 'riskEV', value: 10, tol: 1e-9 }, { kpi: 'mean', value: 110, tol: 0.01 }] },
  ],
};
