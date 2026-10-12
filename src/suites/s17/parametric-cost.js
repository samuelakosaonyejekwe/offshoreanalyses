// Early-phase (top-down) cost estimating: analogous estimate with cost–capacity scaling and index escalation,
// regression cost-estimating relationship from comparable projects, factorial estimate and learning curves.
import { capacityScaling, scalingExponent, indexFactor, factorialEstimate, serialProductionSavings, crawfordCurve, wrightCurve, estimateRange } from '../../engine/econ/cost.js';
import { regressionCostModel, referenceClassUplift } from '../../engine/econ/forecast.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, CO2_STEEL, need, run, round, fmt, mil, fM, pct, adv, check, cost, tableRows, uniqueNames, seriesXY, linspace, range } from './_common.js';

const COMPS = [
  { name: 'Comparable A', capacity: 210, cost: 0.86e9 }, { name: 'Comparable B', capacity: 320, cost: 1.18e9 }, { name: 'Comparable C', capacity: 400, cost: 1.52e9 },
  { name: 'Comparable D', capacity: 480, cost: 1.63e9 }, { name: 'Comparable E', capacity: 590, cost: 2.05e9 }, { name: 'Comparable F', capacity: 720, cost: 2.28e9 },
  { name: 'Comparable G', capacity: 860, cost: 2.87e9 }, { name: 'Comparable H', capacity: 1100, cost: 3.34e9 },
];
const UPLIFT = [{ pOverrun: 50, uplift: 15 }, { pOverrun: 40, uplift: 22 }, { pOverrun: 30, uplift: 30 }, { pOverrun: 20, uplift: 42 }, { pOverrun: 10, uplift: 65 }];

export default {
  id: 's17-parametric-cost', suite: 17, title: 'Parametric, analogous and factorial cost estimating',
  summary: 'Early-phase CAPEX from a reference project scaled by capacity and cost index, a regression on comparable projects, a factorial estimate on equipment cost and the learning-curve cost of repeated units.',
  covers: [2], fidelity: 'screening',
  standards: ['AACE International RP 18R-97 (Class 5 and Class 4 estimating methods)', 'AACE International RP 17R-97', 'AACE International RP 58R-10 (escalation by indices)', 'Williams (1947) six-tenths rule; Lang (1948) factors', 'Wright (1936) and Crawford (1944) learning curves', 'Flyvbjerg (2006) reference-class forecasting'],
  models: ['Analogous cost estimating', 'Parametric cost-estimating relationship (log–log regression)', 'Equipment cost–capacity scaling with exponent', 'Construction cost indices', 'Factorial (Lang) estimation', 'Learning and experience curves (Wright cumulative average, Crawford unit)', 'Reference-class uplift for optimism bias'],
  limits: `Screening methods for concept selection (AACE Class 5–4): valid only inside the capacity range and vintage of the reference data, and only for like-for-like scope, water depth and contract strategy. No quantities are taken off. ${CUR_NOTE}`,
  publishes: ['econ.capexParametric'],
  searchTerms: 'cost capacity scaling exponent learning curve offshore wind cost estimating relationship',
  inputs: [
    { key: 'capacity', label: 'Capacity of the new project', unit: 'MW', type: 'number', default: 500, min: 0.1, max: 50000, group: 'Project information', role: 'operation', bind: 'energy.ratedPower', design: true, help: 'Any consistent capacity measure (MW, t/day, m³/h) may be used if the reference data use the same one.' },
    { key: 'estYear', label: 'Price year of the new estimate', unit: 'yr', type: 'number', default: 2026, min: 1950, max: 2100, group: 'Opening financial state', role: 'ic' },
    { key: 'budgetCap', label: 'Available capital budget', unit: CUR, type: 'number', default: 2.4e9, min: 1, max: 1e12, group: 'Financial constraints', role: 'bc' },

    cost({ key: 'refCost', label: 'Reference project outturn cost', unit: CUR, default: 0.9e9, min: 1, max: 1e12, group: 'Analogous estimate', help: 'Actual final cost of a completed comparable project, not its sanction budget.' }),
    { key: 'refCapacity', label: 'Reference project capacity', unit: 'MW', type: 'number', default: 320, min: 0.1, max: 50000, group: 'Analogous estimate', role: 'operation' },
    { key: 'refYear', label: 'Price year of the reference cost', unit: 'yr', type: 'number', default: 2019, min: 1950, max: 2100, group: 'Analogous estimate', role: 'ic' },
    cost({ key: 'exponent', label: 'Cost–capacity exponent', unit: '-', default: 0.85, min: 0.2, max: 1.2, group: 'Analogous estimate', uncertain: { dist: 'uniform', min: 0.7, max: 0.95 },
      help: '0.6 is the classic equipment value; multi-unit offshore developments scale closer to 0.8–0.95 because units are repeated rather than enlarged.' }),
    cost({ key: 'locationFactor', label: 'Location factor', unit: '-', default: 1.0, min: 0.3, max: 3, group: 'Analogous estimate', help: 'Labour, logistics and regulatory differences between the two sites.' }),
    cost({ key: 'complexityFactor', label: 'Complexity factor', unit: '-', default: 1.05, min: 0.3, max: 3, group: 'Analogous estimate', help: 'Water depth, distance to shore, soil and metocean severity relative to the reference.' }),
    { key: 'index', label: 'Cost index history', type: 'series', group: 'Cost indices', role: 'cost', xLabel: 'Year', yLabel: 'Index', xUnit: 'yr', yUnit: '-',
      help: `Composite construction or equipment cost index covering both price years. ${PLACEHOLDER}`,
      default: { x: [2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026], y: [100, 102.1, 105.4, 107.6, 108.2, 114.9, 128.6, 133.2, 135.8, 138.9, 142.4] } },

    { key: 'comparables', label: 'Comparable completed projects', type: 'table', group: 'Regression (cost-estimating relationship)', role: 'cost',
      help: `Outturn costs normalised to the estimate price year; at least four, segmented to one project type and period. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Project', type: 'text' }, { key: 'capacity', label: 'Capacity', unit: 'MW', type: 'number' }, { key: 'cost', label: 'Outturn cost', unit: CUR, type: 'number' }], default: COMPS },

    cost({ key: 'equipmentCost', label: 'Delivered major equipment cost', unit: CUR, default: 6.6e8, min: 0, max: 1e12, group: 'Factorial estimate' }),
    cost({ key: 'langFactor', label: 'Installed-cost factor on equipment', unit: '-', default: 2.7, min: 1, max: 10, group: 'Factorial estimate', help: 'Total installed cost / delivered equipment cost. Published Lang factors (3.1–4.7) are for onshore process plant; use project history offshore.' }),

    cost({ key: 'firstUnit', label: 'First-unit cost of the repeated item', unit: CUR, default: 7.2e6, min: 0, max: 1e11, group: 'Learning curve', help: 'For example the first foundation fabricated or installed.' }),
    { key: 'units', label: 'Number of identical units', unit: '-', type: 'int', default: 34, min: 1, max: 2000, group: 'Learning curve', role: 'operation', design: true },
    cost({ key: 'learningRate', label: 'Learning rate', unit: '%', default: 8, min: 0, max: 40, group: 'Learning curve', help: 'Cost reduction per doubling of cumulative units.' }),
    { key: 'learningModel', label: 'Learning model', type: 'select', default: 'crawford', group: 'Learning curve', role: 'control', options: [{ value: 'crawford', label: 'Crawford (unit cost)' }, { value: 'wright', label: 'Wright (cumulative average)' }] },

    { key: 'pOverrun', label: 'Acceptable probability of overrun', unit: '%', type: 'number', default: 30, min: 10, max: 50, group: 'Reference-class uplift', role: 'criteria', help: 'Risk appetite of the decision maker: 50 % is a P50 budget, 20 % a P80 budget.' },
    { key: 'upliftTable', label: 'Reference-class uplift table', type: 'table', group: 'Reference-class uplift', role: 'cost',
      help: `Required uplift on the estimate for each acceptable overrun probability, from the outturn distribution of a comparable class of projects. ${PLACEHOLDER}`,
      columns: [{ key: 'pOverrun', label: 'Acceptable overrun probability', unit: '%', type: 'number' }, { key: 'uplift', label: 'Uplift', unit: '%', type: 'number' }], default: UPLIFT },
  ],
  compute(inp) {
    const idx = seriesXY(inp.index, 'Cost index history', 2);
    need(idx.y.every((v) => v > 0), 'Cost index values must be positive.');
    for (let k = 1; k < idx.x.length; k++) need(idx.x[k] > idx.x[k - 1], 'Cost index years must be ascending.');
    const lo = idx.x[0]; const hi = idx.x[idx.x.length - 1];
    need(inp.refYear >= lo && inp.refYear <= hi && inp.estYear >= lo && inp.estYear <= hi, `Both price years (${inp.refYear} and ${inp.estYear}) must lie inside the cost index history (${lo}–${hi}). Extend the index or use the forecasting module to project it.`);
    const comps = tableRows(inp.comparables, 'Comparable projects', { capacity: { min: 1e-9 }, cost: { min: 1e-9 } }, { min: 4 });
    uniqueNames(comps, 'Comparable projects');
    need(new Set(comps.map((c) => c.capacity)).size >= 3, 'Comparable projects: at least three different capacities are needed to fit a scaling law.');
    const up = tableRows(inp.upliftTable, 'Uplift table', { pOverrun: { min: 0.1, max: 99 }, uplift: { min: -50, max: 500 } }, { min: 2 });

    // (1) analogous estimate: scale by capacity, escalate by index, adjust for location and complexity
    const ix = run(() => indexFactor({ years: idx.x, values: idx.y }, inp.refYear, inp.estYear), 'Cost index');
    const analog = run(() => capacityScaling({ C1: inp.refCost, S1: inp.refCapacity, S2: inp.capacity, n: inp.exponent, index1: ix.from, index2: ix.to, locationFactor: inp.locationFactor, complexityFactor: inp.complexityFactor }), 'Capacity scaling');
    // (2) regression CER: ln C = b0 + b1 ln S; 80 % prediction interval read as cumulative P10–P90
    const reg = run(() => regressionCostModel({ X: comps.map((c) => [c.capacity]), y: comps.map((c) => c.cost), logLog: true }), 'Regression');
    const pred = reg.predict([inp.capacity], 0.8);
    const cmin = Math.min(...comps.map((c) => c.capacity)); const cmax = Math.max(...comps.map((c) => c.capacity));
    const extrap = inp.capacity < cmin || inp.capacity > cmax;
    const ends = comps.slice().sort((a, b) => a.capacity - b.capacity);
    const nEnds = scalingExponent({ C1: ends[0].cost, S1: ends[0].capacity, C2: ends[ends.length - 1].cost, S2: ends[ends.length - 1].capacity });
    // (3) factorial estimate on equipment cost (factor supplied by the user)
    const fac = inp.equipmentCost > 0 ? run(() => factorialEstimate({ method: 'lang', equipment: [{ name: 'Major equipment', cost: inp.equipmentCost }], plantType: 'project', factors: { project: inp.langFactor } }), 'Factorial estimate').total : 0;
    // (4) learning on repeated units
    const lr = inp.learningRate / 100;
    const learn = run(() => serialProductionSavings({ firstUnitCost: inp.firstUnit, n: inp.units, learningRate: lr, model: inp.learningModel }), 'Learning curve');
    const curve = inp.learningModel === 'wright' ? wrightCurve : crawfordCurve;
    const step = Math.max(1, Math.ceil(inp.units / 200));
    const un = range(Math.ceil(inp.units / step), (k) => Math.min(inp.units, 1 + k * step));
    const unitCost = un.map((u) => curve({ T1: inp.firstUnit, n: u, learningRate: lr }).unit);

    const methods = [['Analogous (scaled reference)', analog], ['Regression on comparables', pred.mean]];
    if (fac > 0) methods.push(['Factorial on equipment', fac]);
    const vals = methods.map((m) => m[1]);
    const central = vals.reduce((a, b) => a + b, 0) / vals.length;
    const spread = (Math.max(...vals) - Math.min(...vals)) / central;
    const tbl = { acceptableOverrunProbability: up.map((r) => r.pOverrun / 100), uplift: up.map((r) => r.uplift / 100) };
    const rc = run(() => referenceClassUplift({ estimate: central, acceptableOverrunProbability: inp.pOverrun / 100, table: tbl }), 'Reference-class uplift (the acceptable overrun probability must lie inside the table)');
    const rng = estimateRange(central, 5);

    const kpis = [
      { key: 'analog', label: 'Analogous estimate', value: mil(analog), unit: MU, explain: `Reference cost scaled by (capacity ratio)^${inp.exponent}, escalated by the index ratio ${fmt(ix.factor, 4)} and adjusted for location and complexity. Only as good as the likeness of the reference project.` },
      { key: 'regression', label: 'Regression estimate (expected)', value: mil(pred.mean), unit: MU, explain: 'Cost predicted by the log–log relationship fitted to the comparable projects at the new capacity. It averages over project-specific luck that a single reference carries.' },
      { key: 'regLow', label: 'Regression P10 (cumulative)', value: mil(pred.piLow), unit: MU, explain: 'Lower end of the 80 % prediction interval for one new project: the cost not exceeded with 10 % probability, given the scatter of the comparables.' },
      { key: 'regHigh', label: 'Regression P90 (cumulative)', value: mil(pred.piHigh), unit: MU, explain: 'Upper end of the 80 % prediction interval: the cost not exceeded with 90 % probability. This is a statistical band around the fit, not an allowance for scope differences.' },
      { key: 'fitExponent', label: 'Fitted cost–capacity exponent', value: round(reg.beta[1]), unit: '-', explain: 'Slope of ln(cost) on ln(capacity) in the comparables. Below 1 there are economies of scale; compare it with the exponent assumed for the analogous estimate.' },
      { key: 'r2', label: 'Regression R²', value: round(reg.r2), unit: '-', explain: 'Share of the variance of ln(cost) explained by capacity. A low value means capacity alone is a poor cost driver for this data set and further drivers (depth, distance) are needed.' },
      { key: 'indexFactor', label: 'Escalation factor between price years', value: round(ix.factor, 6), unit: '-', explain: `Index at ${inp.estYear} over index at ${inp.refYear}. Mixing price years without this correction is one of the commonest errors in analogous estimates.` },
      { key: 'factorial', label: 'Factorial estimate', value: mil(fac), unit: MU, explain: 'Delivered equipment cost multiplied by the installed-cost factor. Quick, but the factor must come from projects of the same kind.' },
      { key: 'central', label: 'Central early estimate (mean of methods)', value: mil(central), unit: MU, limit: mil(inp.budgetCap), sense: 'max', explain: 'Average of the independent methods; treat it as a Class 5 figure with the wide accuracy band of that class.' },
      { key: 'spread', label: 'Spread between methods', value: round(spread * 100), unit: '%', explain: 'Range of the methods relative to their mean. A large spread says the reference data disagree and no single figure should be quoted without its range.' },
      { key: 'uplifted', label: `Budget at ${inp.pOverrun} % overrun probability`, value: mil(rc.upliftedEstimate), unit: MU, explain: 'Central estimate plus the reference-class uplift for the chosen risk appetite: the outside-view correction for optimism bias.' },
      { key: 'learnTotal', label: `Cost of ${inp.units} repeated units`, value: mil(learn.total), unit: MU, explain: 'Total cost of the series with learning; compare with the no-learning figure in the saving KPI.' },
      { key: 'learnSaving', label: 'Saving from learning', value: mil(learn.saving), unit: MU, explain: `Difference from ${inp.units} × first-unit cost. Realised only if the units really are identical and built without long breaks.` },
      { key: 'learnAvg', label: 'Average unit cost', value: mil(learn.average), unit: MU, explain: 'Series total divided by the number of units: the rate to carry into a bottom-up estimate instead of the first-unit cost.' },
    ];
    const checks = [
      check('Central estimate within budget', central, inp.budgetCap, CUR, 'Financial constraint (affordability limit)', 'Guards against carrying a concept forward that is not affordable even on its central early estimate.'),
      check('Risk-adjusted budget within available capital', rc.upliftedEstimate, inp.budgetCap, CUR, 'Flyvbjerg (2006) reference-class forecasting', 'Guards against optimism bias: if the uplifted figure exceeds the capital available, the chance of overrun is higher than the decision maker accepts.'),
      check('Agreement between methods (spread ≤ 30 %)', spread * 100, 30, '%', 'AACE RP 18R-97 (Class 5 low-side/high-side accuracy)', 'Guards against quoting one method when the others disagree; a spread beyond the class accuracy means the reference data are not alike.'),
    ];
    const sx = linspace(Math.min(cmin, inp.capacity) * 0.8, Math.max(cmax, inp.capacity) * 1.15, 60);
    const plots = [
      { type: 'line', title: 'Cost against capacity', x: { label: 'Capacity', unit: 'MW' }, y: { label: 'Cost', unit: MU },
        series: [{ name: 'Regression fit', x: sx.map((v) => round(v)), y: sx.map((s) => mil(reg.predict([s], 0.8).mean)) }, { name: 'Regression P10', x: sx.map((v) => round(v)), y: sx.map((s) => mil(reg.predict([s], 0.8).piLow)), dash: true },
          { name: 'Regression P90', x: sx.map((v) => round(v)), y: sx.map((s) => mil(reg.predict([s], 0.8).piHigh)), dash: true },
          { name: `Scaled reference (n = ${inp.exponent})`, x: sx.map((v) => round(v)), y: sx.map((s) => mil(analog * (s / inp.capacity) ** inp.exponent)) },
          { name: 'Comparables', x: comps.map((c) => c.capacity), y: comps.map((c) => mil(c.cost)) }, { name: 'New project', x: [inp.capacity], y: [mil(central)] }],
        explain: 'Comparable projects, the fitted relationship with its 80 % prediction band (cumulative P10 and P90) and the scaled reference. Look at whether the new capacity sits inside the data and whether the two laws diverge there.' },
      { type: 'line', title: 'Learning curve of the repeated unit', x: { label: 'Unit number', unit: '-' }, y: { label: 'Unit cost', unit: MU },
        series: [{ name: 'Unit cost', x: un, y: unitCost.map(mil) }, { name: 'First-unit cost', x: [1, inp.units], y: [mil(inp.firstUnit), mil(inp.firstUnit)], dash: true }],
        explain: 'Cost of each successive unit. The steep early drop is where serial fabrication pays; the curve flattens, so extra units late in the series add little further saving.' },
      { type: 'line', title: 'Cost index', x: { label: 'Year', unit: '' }, y: { label: 'Index', unit: '-' }, series: [{ name: 'Index', x: idx.x, y: idx.y }, { name: 'Price years', x: [inp.refYear, inp.estYear], y: [ix.from, ix.to] }],
        explain: 'The escalation index with the two price years marked; the ratio of the marked values is the factor applied to the reference cost.' },
      { type: 'bar', title: 'Estimates by method', categories: [...methods.map((m) => m[0]), 'Central', 'Risk-adjusted', 'Budget'], x: { label: '' }, y: { label: 'CAPEX', unit: MU }, series: [{ name: 'CAPEX', y: [...vals, central, rc.upliftedEstimate, inp.budgetCap].map(mil) }],
        explain: 'Each independent method beside the central value, the reference-class budget and the capital available. Agreement between the first bars is the measure of confidence.' },
    ];
    const tables = [
      { title: 'Estimates by method', columns: [{ key: 'm', label: 'Method' }, { key: 'v', label: 'Estimate', unit: MU }, { key: 'd', label: 'Difference from central', unit: '%' }, { key: 'basis', label: 'Basis' }],
        rows: [{ m: methods[0][0], v: mil(analog), d: round(100 * (analog / central - 1), 3), basis: `${fM(inp.refCost)} × (${inp.capacity}/${inp.refCapacity})^${inp.exponent} × ${fmt(ix.factor, 4)} × ${inp.locationFactor} × ${inp.complexityFactor}` },
          { m: methods[1][0], v: mil(pred.mean), d: round(100 * (pred.mean / central - 1), 3), basis: `ln C = ${fmt(reg.beta[0], 4)} + ${fmt(reg.beta[1], 4)} ln S, R² = ${fmt(reg.r2)}, n = ${comps.length}${extrap ? ' (extrapolated)' : ''}` },
          ...(fac > 0 ? [{ m: 'Factorial on equipment', v: mil(fac), d: round(100 * (fac / central - 1), 3), basis: `${fM(inp.equipmentCost)} × ${inp.langFactor}` }] : []),
          { m: 'AACE Class 5 range of the central value', v: mil(central), d: 0, basis: `${fM(rng.wide.low)} to ${fM(rng.wide.high)} (wide band)` }],
        explain: 'Each method with its arithmetic, so the estimate can be reproduced by hand and the weakest assumption identified.' },
      { title: 'Comparable projects and residuals', columns: [{ key: 'name', label: 'Project' }, { key: 'capacity', label: 'Capacity', unit: 'MW' }, { key: 'cost', label: 'Outturn cost', unit: MU }, { key: 'fit', label: 'Fitted', unit: MU }, { key: 'res', label: 'Residual', unit: '%' }, { key: 'perMW', label: 'Cost per MW', unit: `${MU}/MW` }],
        rows: comps.map((c) => { const f = reg.predict([c.capacity]).mean; return { name: c.name, capacity: c.capacity, cost: mil(c.cost), fit: mil(f), res: round(100 * (c.cost / f - 1), 3), perMW: round(c.cost / c.capacity / 1e6, 4) }; }),
        explain: 'The data behind the regression with each project’s departure from the fit. A large residual flags a project whose scope or period differs and should be segmented out.' },
    ];
    const advice = [
      adv('engineering', extrap ? 'warn' : 'info', `The fitted exponent is ${fmt(reg.beta[1])} (two-point check between the smallest and largest comparable: ${fmt(nEnds)}) against ${inp.exponent} assumed; the new capacity of ${inp.capacity} MW is ${extrap ? 'outside' : 'inside'} the ${fmt(cmin)}–${fmt(cmax)} MW range of the data.`,
        extrap ? 'Add comparables that bracket the new capacity or rely on the bottom-up estimate: a power law extrapolated beyond its data is not evidence.' : `Capacity explains ${pct(reg.r2)} of the cost variance; add the next technical driver (water depth, distance to shore, unit rating) before using the relationship beyond concept selection.`,
        'Parametric estimating, AACE RP 18R-97 (Class 5/4 methods); Williams (1947) cost–capacity relation.',
        `Serial production of ${inp.units} identical units saves ${fM(learn.saving)} (${pct(learn.savingFraction)}) of fabrication effort; standardising the design also cuts rework and scrap — each 1 % of steel not wasted on a 1000 t unit is about ${fmt(10 * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t CO₂e/t).`),
      adv('economic', spread > 0.3 || rc.upliftedEstimate > inp.budgetCap ? 'warn' : 'ok', `Methods give ${methods.map((m) => `${fM(m[1])} (${m[0].toLowerCase()})`).join(', ')}: central ${fM(central)}, spread ${pct(spread)}. With a ${inp.pOverrun} % accepted chance of overrun the budget should be ${fM(rc.upliftedEstimate)} (uplift ${pct(rc.uplift)}).`,
        rc.upliftedEstimate > inp.budgetCap ? `The risk-adjusted figure exceeds the ${fM(inp.budgetCap)} available by ${fM(rc.upliftedEstimate - inp.budgetCap)}: reduce capacity, phase the development or accept a higher overrun probability explicitly.` : 'Carry the central value as a Class 5 estimate with its range into concept selection and replace it with the bottom-up estimate as definition matures.',
        'AACE RP 18R-97 Table 1 (Class 5: −20 % to −50 % low, +30 % to +100 % high); reference-class forecasting, Flyvbjerg (2006).'),
    ];
    return { kpis, checks, plots, tables, advice,
      notes: [`All estimates are in ${inp.estYear} money; the reference cost was moved from ${inp.refYear} with the index ratio ${fmt(ix.factor, 4)}.`,
        'Regression P10 and P90 are cumulative percentiles (cost not exceeded with 10 % and 90 % probability) of the prediction for a single new project.',
        'Indices, factors, exponents and comparable costs are indicative placeholders to be replaced by the estimator’s own data.'],
      publish: { 'econ.capexParametric': central } };
  },
  verify: [
    // C2 = 100 × 2^0.6 = 151.572 (index ratio 1, factors 1)
    { name: 'Six-tenths rule: doubling capacity multiplies cost by 2^0.6', source: 'Williams (1947): C2 = C1 (S2/S1)^0.6 = 100 × 1.515717 = 151.5717',
      inputs: { refCost: 100e6, refCapacity: 100, capacity: 200, exponent: 0.6, refYear: 2020, estYear: 2020, locationFactor: 1, complexityFactor: 1 }, expect: [{ kpi: 'analog', value: 151.5717, tol: 1e-5 }, { kpi: 'indexFactor', value: 1, tol: 1e-9 }] },
    // index 107.6 (2019) → 142.4 (2026): factor 1.323420
    { name: 'Index escalation between price years', source: 'Ratio of tabulated index values 142.4 / 107.6 = 1.323420', inputs: { refYear: 2019, estYear: 2026 }, expect: [{ kpi: 'indexFactor', value: 1.32342, tol: 1e-5 }] },
    // Crawford, 10 % learning: b = log2(0.9) = −0.152003; unit costs 1, 0.9, 3^b = 0.846206, 0.81 → total 3.556206
    { name: 'Crawford unit curve, four units at 10 % learning', source: 'Hand calculation: T1 (1 + 0.9 + 3^log2(0.9) + 0.81) = 3.556206 T1', inputs: { firstUnit: 1e6, units: 4, learningRate: 10, learningModel: 'crawford' }, expect: [{ kpi: 'learnTotal', value: 3.556206, tol: 1e-5 }] },
    // exact power law C = 2e6 S^0.8 in the comparables → fitted exponent 0.8, R² = 1, prediction at 500 MW = 2e6 × 500^0.8 = 288.540 M
    { name: 'Regression recovers an exact power law', source: 'Synthetic data C = 2e6 × S^0.8: slope 0.8 and C(500) = 288.540 M',
      inputs: { capacity: 500, comparables: [100, 200, 400, 800, 1600].map((s, k) => ({ name: `P${k}`, capacity: s, cost: 2e6 * s ** 0.8 })) }, expect: [{ kpi: 'fitExponent', value: 0.8, tol: 1e-6 }, { kpi: 'regression', value: 288.540, tol: 1e-4 }] },
  ],
};
