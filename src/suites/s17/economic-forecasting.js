// Economic forecasting and data-driven models: time-series forecasts of a price or index with hold-out
// back-testing (Holt–Winters, ARIMA with drift, seasonal ARIMA, local-level Kalman filter, vector autoregression),
// stochastic price processes with a percentile cone (geometric Brownian motion, Schwartz mean reversion, jump
// diffusion), cost-estimating relationships (multivariate regression, gradient boosting, random forest, Gaussian
// process, neural network) compared by cross-validation, partial pooling of sparse benchmark data, a Poisson
// regression of maintenance demand and carbon-price scenario paths.
import { arima } from '../../engine/math/timeseries.js';
import { makeRng } from '../../engine/math/rng.js';
import { gbmSimulate, gbmEstimate, schwartzSimulate, schwartzEstimate, jumpDiffusionSimulate, scenarioCone, seasonalDecompose, errorMetrics, holtWintersForecast, arimaForecast, backtest, kalmanLocalLevel, var1Fit,
  regressionCostModel, glmFit, gradientBoosting, randomForest, mlpRegression, hierarchicalPooling, carbonPricePaths } from '../../engine/econ/forecast.js';
import { carbonCost } from '../../engine/econ/carbon.js';
import { gaussianProcess } from '../../engine/reliability/surrogate.js';
import { CUR, MU, CUR_NOTE, PLACEHOLDER, PCT_CUM, need, run, round, fmt, mil, fM, pct, adv, check, checkMin, cost, tableRows, uniqueNames, seriesXY, percentiles, range, sum } from './_common.js';

// Illustrative quarterly price history and a correlated driver index: a rising level, a seasonal pattern and
// mean-reverting noise from a fixed seed. Placeholders for the user's own records.
const HIST = (() => {
  const g = makeRng(2026); const price = []; const driver = []; let a = 0; let b = 0;
  for (let t = 0; t < 48; t++) {
    a = 0.72 * a + 0.085 * g.normal(); b = 0.6 * b + 0.05 * g.normal();
    price.push(Number(((62 + 0.45 * t) * (1 + [0.07, -0.03, -0.06, 0.02][t % 4]) * Math.exp(a)).toFixed(2)));
    driver.push(Number((100 * Math.exp(0.004 * t + 0.6 * a + b)).toFixed(2)));
  }
  return { x: range(48), price, driver };
})();
// Illustrative completed-project records: cost rises with capacity, water depth and distance to shore.
const PROJECTS = (() => {
  const g = makeRng(404); const rows = [];
  for (let k = 0; k < 36; k++) {
    const mw = Math.round(150 + 950 * g.next()); const depth = Math.round(15 + 45 * g.next()); const dist = Math.round(10 + 110 * g.next());
    const c = 2.6e6 * mw ** 0.93 * (1 + 0.011 * (depth - 30)) * (1 + 0.0024 * dist) * Math.exp(0.07 * g.normal());
    rows.push({ name: `Project ${String(k + 1).padStart(2, '0')}`, mw, depth, dist, cost: Math.round(c / 1e5) * 1e5 });
  }
  return rows;
})();
const GROUPS = [
  { name: 'Fixed-bottom wind', mean: 12, variance: 180, n: 14 }, { name: 'Floating wind', mean: 31, variance: 420, n: 4 },
  { name: 'Subsea tie-back', mean: 18, variance: 260, n: 9 }, { name: 'Export cable and grid', mean: 24, variance: 350, n: 6 },
];
const INTERVENTIONS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10].map((age, k) => ({ age, count: [21, 19, 24, 26, 25, 31, 33, 38, 41, 47][k] }));
const TS_MODELS = { holt: 'Holt–Winters exponential smoothing', arima: 'ARIMA', sarima: 'Seasonal ARIMA', kalman: 'Local-level Kalman filter' };
const COST_MODELS = { regression: 'Multivariate regression', boosting: 'Gradient boosting', forest: 'Random forest', gp: 'Gaussian process', mlp: 'Neural network' };

/** Lag-0 and lag-1 autocovariance of the first differences, for the moment estimate of the local-level variances. */
function diffMoments(y) {
  const d = []; for (let t = 1; t < y.length; t++) d.push(y[t] - y[t - 1]);
  const m = sum(d) / d.length; let g0 = 0; let g1 = 0;
  for (let t = 0; t < d.length; t++) { g0 += (d[t] - m) ** 2; if (t > 0) g1 += (d[t] - m) * (d[t - 1] - m); }
  return { g0: g0 / d.length, g1: g1 / d.length };
}

export default {
  id: 's17-economic-forecasting', suite: 17, title: 'Price, cost and maintenance-demand forecasting',
  summary: 'Forecasts a price or cost index from its history with back-tested time-series models and a stochastic percentile cone, converts it into revenue at risk, and fits data-driven cost and maintenance-demand models with honest out-of-sample error.',
  covers: [11], fidelity: 'numerical',
  standards: ['Hyndman & Athanasopoulos, Forecasting: Principles and Practice, 3rd ed., ch. 5, 8 and 9 (accuracy evaluation, exponential smoothing, ARIMA)', 'Box, Jenkins & Reinsel, Time Series Analysis, 4th ed., ch. 5, 7 and 9', 'Harvey (1989), Forecasting, Structural Time Series Models and the Kalman Filter', 'Schwartz (1997), J. Finance 52(3) (mean-reverting commodity prices); Merton (1976) jump diffusion', 'ISPA/SCEA Parametric Estimating Handbook, 4th ed. (cost-estimating relationships and their validation)', 'AACE International RP 58R-10 (escalation estimating using indices and forecasts)', 'Gelman et al., Bayesian Data Analysis, 3rd ed., ch. 5 (hierarchical models)'],
  models: ['Autoregressive integrated moving average (ARIMA) with drift', 'Seasonal ARIMA', 'Holt–Winters exponential smoothing', 'State-space local-level model with Kalman filtering', 'Vector autoregression', 'Classical seasonal decomposition', 'Geometric Brownian motion', 'Mean-reverting (Schwartz one-factor) price process', 'Jump diffusion', 'Hold-out back-testing against a naive benchmark', 'Multivariate regression cost-estimating relationship with prediction interval', 'Generalised linear regression (Poisson)', 'Gradient boosting', 'Random forest', 'Gaussian process regression', 'Neural-network regression', 'Bayesian hierarchical (partial-pooling) model', 'Carbon-price scenario paths'],
  limits: `Forecasts extrapolate the statistical structure of the history supplied: structural breaks, policy changes and market redesign are outside any of these models and must enter as scenarios. The tree, Gaussian-process and neural-network models are fitted to a small table and are reported with their cross-validated error; they are not to be used outside the range of the data. The Gaussian-process length scales are fitted once to all records and held in the cross-validation folds, which flatters its error slightly. The seasonal ARIMA limits are those of the differenced series and understate the true uncertainty. ${CUR_NOTE}`,
  publishes: ['econ.priceForecast', 'econ.priceP10', 'econ.priceP90', 'econ.priceVolatility', 'econ.revenueForecast', 'econ.capexDataDriven'],
  searchTerms: 'electricity price forecasting ARIMA mean reversion cost estimating relationship machine learning offshore',
  mesh: { key: 'nPaths', label: 'Simulated price paths', levels: [500, 1000, 2000, 4000], qoi: 'meanEnd' },
  inputs: [
    { key: 'firstYear', label: 'Year of the first observation', unit: 'yr', type: 'int', default: 2014, min: 1950, max: 2100, group: 'Opening economic state', role: 'ic' },
    { key: 'basis', label: 'Terms of the history', type: 'select', default: 'real', group: 'Opening economic state', role: 'ic', options: [{ value: 'real', label: 'Real (constant money)' }, { value: 'nominal', label: 'Nominal (money of the day)' }],
      help: 'Forecasts are in the same terms as the history. Deflate a nominal series before using the result with a real discount rate.' },
    { key: 'perYear', label: 'Observations per year', type: 'select', default: '4', group: 'Opening economic state', role: 'ic', options: [{ value: '1', label: 'Annual' }, { value: '4', label: 'Quarterly' }, { value: '12', label: 'Monthly' }], help: 'Sets the seasonal period and the time step of the stochastic processes.' },
    { key: 'history', label: 'Price or index history', type: 'series', group: 'Price history', role: 'cost', xLabel: 'Period', yLabel: 'Price', xUnit: '-', yUnit: `${CUR}/MWh`, default: { x: HIST.x, y: HIST.price },
      help: `Equally spaced observations, oldest first, at least 16. ${PLACEHOLDER}` },
    { key: 'driver', label: 'Related driver series (optional)', type: 'series', group: 'Price history', role: 'cost', xLabel: 'Period', yLabel: 'Driver index', xUnit: '-', yUnit: '-', default: { x: HIST.x, y: HIST.driver },
      help: `A second series observed at the same dates (fuel price, steel or vessel-rate index) for the vector autoregression. Leave empty to skip. ${PLACEHOLDER}` },

    { key: 'horizon', label: 'Forecast horizon', unit: 'periods', type: 'int', default: 12, min: 1, max: 120, group: 'Time-series models', role: 'control' },
    { key: 'holdout', label: 'Hold-out length for back-testing', unit: 'periods', type: 'int', default: 8, min: 2, max: 60, group: 'Time-series models', role: 'control', help: 'The last observations withheld from fitting and used to score each model.' },
    { key: 'model', label: 'Forecast model', type: 'select', default: 'auto', group: 'Time-series models', role: 'control',
      options: [{ value: 'auto', label: 'Best on the hold-out test' }, ...Object.entries(TS_MODELS).map(([value, label]) => ({ value, label }))] },
    { key: 'arP', label: 'ARIMA order p', unit: '-', type: 'int', default: 1, min: 0, max: 4, group: 'Time-series models', role: 'control' },
    { key: 'arD', label: 'ARIMA order d', unit: '-', type: 'int', default: 1, min: 0, max: 2, group: 'Time-series models', role: 'control' },
    { key: 'arQ', label: 'ARIMA order q', unit: '-', type: 'int', default: 0, min: 0, max: 3, group: 'Time-series models', role: 'control' },
    { key: 'drift', label: 'ARIMA with drift', type: 'bool', default: true, group: 'Time-series models', role: 'control', help: 'Estimates the mean of the differenced series, so that a trending history is extrapolated rather than held flat.' },

    { key: 'process', label: 'Stochastic price process', type: 'select', default: 'auto', group: 'Stochastic price cone', role: 'control',
      options: [{ value: 'auto', label: 'Mean reverting if the history supports it, else random walk' }, { value: 'gbm', label: 'Geometric Brownian motion' }, { value: 'schwartz', label: 'Mean reverting (Schwartz)' }, { value: 'jump', label: 'Jump diffusion' }] },
    { key: 'jumpRate', label: 'Jump frequency', unit: '1/yr', type: 'number', default: 0.5, min: 0, max: 20, group: 'Stochastic price cone', role: 'ic', showIf: { key: 'process', in: ['jump'] } },
    { key: 'jumpMean', label: 'Mean jump size (log)', unit: '%', type: 'number', default: -8, min: -90, max: 200, group: 'Stochastic price cone', role: 'ic', showIf: { key: 'process', in: ['jump'] } },
    { key: 'jumpStd', label: 'Jump size standard deviation (log)', unit: '%', type: 'number', default: 12, min: 0, max: 200, group: 'Stochastic price cone', role: 'ic', showIf: { key: 'process', in: ['jump'] } },
    { key: 'nPaths', label: 'Simulated price paths', unit: '-', type: 'int', default: 1000, min: 200, max: 20000, group: 'Stochastic price cone', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 17, min: 1, max: 1e6, group: 'Stochastic price cone', role: 'control', help: 'Fixed seed: identical inputs reproduce identical results.' },

    { key: 'volume', label: 'Annual volume sold', unit: 'GWh/yr', type: 'number', default: 2150, min: 0.001, max: 2e5, group: 'Revenue boundary', role: 'operation', bind: 'energy.aep', help: 'Net annual production from the energy-yield suite when available.' },
    { key: 'contractShare', label: 'Share of volume under a fixed-price contract', unit: '%', type: 'number', default: 70, min: 0, max: 100, group: 'Revenue boundary', role: 'bc', design: true, help: 'Revenue boundary: the contracted share earns the contract price, the remainder the forecast market price.' },
    cost({ key: 'contractPrice', label: 'Contract price', unit: `${CUR}/MWh`, default: 102, min: 0, max: 1e5, group: 'Revenue boundary' }),
    cost({ key: 'breakEven', label: 'Break-even price of the project', unit: `${CUR}/MWh`, default: 82, min: 0, max: 1e5, group: 'Revenue boundary', bind: 'econ.breakEvenTariff', help: 'Benchmark the downside price is checked against; from the levelised-cost module when available.' }),

    { key: 'projects', label: 'Completed projects', type: 'table', group: 'Data-driven cost model', role: 'cost', help: `Final cost of comparable projects in the price year of the estimate, with the three drivers. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Project', type: 'text' }, { key: 'mw', label: 'Capacity', unit: 'MW', type: 'number' }, { key: 'depth', label: 'Water depth', unit: 'm', type: 'number' }, { key: 'dist', label: 'Distance to shore', unit: 'km', type: 'number' }, { key: 'cost', label: 'Final cost', unit: CUR, type: 'number' }], default: PROJECTS },
    { key: 'newMW', label: 'New project: capacity', unit: 'MW', type: 'number', default: 500, min: 0.1, max: 20000, group: 'Data-driven cost model', role: 'operation', bind: 'energy.ratedPower' },
    { key: 'newDepth', label: 'New project: water depth', unit: 'm', type: 'number', default: 38, min: 0.1, max: 3000, group: 'Data-driven cost model', role: 'environment', bind: 'site.depth' },
    { key: 'newDist', label: 'New project: distance to shore', unit: 'km', type: 'number', default: 65, min: 0.1, max: 1000, group: 'Data-driven cost model', role: 'geometry', bind: 'cable.length' },
    { key: 'costForm', label: 'Form of the regression', type: 'select', default: 'loglog', group: 'Data-driven cost model', role: 'control', options: [{ value: 'loglog', label: 'Power law (log–log)' }, { value: 'linear', label: 'Linear' }] },
    { key: 'costModel', label: 'Cost model reported', type: 'select', default: 'auto', group: 'Data-driven cost model', role: 'control', options: [{ value: 'auto', label: 'Lowest cross-validated error' }, ...Object.entries(COST_MODELS).map(([value, label]) => ({ value, label }))] },
    { key: 'folds', label: 'Cross-validation folds', unit: '-', type: 'int', default: 3, min: 2, max: 10, group: 'Data-driven cost model', role: 'control' },

    { key: 'groups', label: 'Observed cost growth by project class', type: 'table', group: 'Pooling of sparse benchmarks', role: 'cost', help: `Mean and variance of estimate-to-final cost growth in percent, and the number of projects observed, by class. ${PLACEHOLDER}`,
      columns: [{ key: 'name', label: 'Class', type: 'text' }, { key: 'mean', label: 'Mean growth', unit: '%', type: 'number' }, { key: 'variance', label: 'Variance', unit: '%²', type: 'number' }, { key: 'n', label: 'Projects', type: 'number' }], default: GROUPS },

    { key: 'interventions', label: 'Maintenance interventions by asset age', type: 'table', group: 'Maintenance-demand model', role: 'ic', help: 'Recorded offshore interventions per year against the age of the asset: the opening maintenance state from which demand is projected.',
      columns: [{ key: 'age', label: 'Asset age', unit: 'yr', type: 'number' }, { key: 'count', label: 'Interventions', unit: '1/yr', type: 'number' }], default: INTERVENTIONS },
    { key: 'glmYears', label: 'Years of maintenance demand to project', unit: 'yr', type: 'int', default: 5, min: 1, max: 30, group: 'Maintenance-demand model', role: 'control' },
    cost({ key: 'interventionCost', label: 'Cost per intervention', unit: CUR, default: 1.4e5, min: 0, max: 1e9, group: 'Maintenance-demand model', help: 'Vessel, crew and parts for one offshore visit.' }),

    { key: 'emissions', label: 'Emissions exposed to a carbon price', unit: 't CO₂e/yr', type: 'number', default: 11500, min: 0, max: 1e9, group: 'Carbon-price scenarios', role: 'operation', help: 'Operational emissions that carry, or may come to carry, a carbon price.' },
    cost({ key: 'carbonStart', label: 'Carbon price today', unit: `${CUR}/t CO₂e`, default: 85, min: 0, max: 5000, group: 'Carbon-price scenarios' }),
    cost({ key: 'carbonLow', label: 'Carbon price growth, low scenario', unit: '%/yr', default: 0, min: -20, max: 50, group: 'Carbon-price scenarios' }),
    cost({ key: 'carbonMid', label: 'Carbon price growth, central scenario', unit: '%/yr', default: 4, min: -20, max: 50, group: 'Carbon-price scenarios' }),
    cost({ key: 'carbonHigh', label: 'Carbon price growth, high scenario', unit: '%/yr', default: 9, min: -20, max: 50, group: 'Carbon-price scenarios' }),
    { key: 'carbonYears', label: 'Years of carbon cost considered', unit: 'yr', type: 'int', default: 25, min: 1, max: 60, group: 'Carbon-price scenarios', role: 'bc', bind: 'life.designYears' },
    { key: 'discountRate', label: 'Discount rate', unit: '%/yr', type: 'number', default: 6, min: -5, max: 50, group: 'Carbon-price scenarios', role: 'bc', help: `In the terms of the history. ${PLACEHOLDER}` },
  ],
  compute(inp) {
    const y = seriesXY(inp.history, 'Price history', 16).y;
    need(y.every((v) => v > 0), 'Price history: every value must be positive (the stochastic processes work on the logarithm of the price).');
    const n = y.length; const m = Number(inp.perYear); const dt = 1 / m; const h = Math.round(inp.horizon);
    const seasonal = m > 1 && n >= 2 * m + 4;
    const hold = Math.round(inp.holdout);
    need(hold <= n - Math.max(12, seasonal ? 2 * m + 4 : 12), `Hold-out of ${hold} periods leaves too little history to fit: with ${n} observations use at most ${n - Math.max(12, seasonal ? 2 * m + 4 : 12)}.`);

    /* ---- seasonal pattern and seasonally adjusted series ---- */
    const sIdx = seasonal ? seasonalDecompose(y, m).seasonalIndices : null;
    const season = (t) => (sIdx ? sIdx[((t % m) + m) % m] : 0);
    let adj = y.map((v, t) => v - season(t));
    const adjusted = adj.every((v) => v > 0);
    if (!adjusted) adj = y.slice(); // additive adjustment would turn a low price negative: keep the raw series

    /* ---- (A) time-series models, each as (training series, steps ahead) → forecast ---- */
    const fitArima = (x) => arima(x, { p: inp.arP, d: inp.arD, q: inp.arQ, includeMean: inp.arD === 0 ? true : !!inp.drift });
    const models = {
      holt: (x, k) => holtWintersForecast({ x, horizon: k, ...(m > 1 && x.length >= 2 * m + 2 ? { season: m, seasonal: 'additive' } : {}) }).forecast,
      // the engine's ARIMA wrapper drops the drift term whenever d ≥ 1, so the time-series library is called directly
      arima: (x, k) => Array.from(fitArima(x).forecast(k).mean),
      ...(seasonal ? { sarima: (x, k) => arimaForecast({ x, p: inp.arP, d: 0, q: inp.arQ, horizon: k, season: m, D: 1 }).mean } : {}),
      kalman: (x, k) => {
        const si = m > 1 && x.length >= 2 * m ? seasonalDecompose(x, m).seasonalIndices : null;
        const xa = si ? x.map((v, t) => v - si[t % m]) : x;
        // local level: Var(Δy) = q + 2r and Cov(Δy_t, Δy_t−1) = −r (Harvey 1989, §2.4), solved for q and r
        const { g0, g1 } = diffMoments(xa); const floor = 1e-9 * (g0 + 1e-12) + 1e-12;
        const r = Math.max(-g1, floor); const qv = Math.max(g0 - 2 * r, floor);
        const f = kalmanLocalLevel({ y: xa, q: qv, r, horizon: k }).forecast;
        return f.map((v, j) => v + (si ? si[(x.length + j) % m] : 0));
      },
    };
    const tests = Object.keys(models).map((key) => {
      try {
        const bt = backtest({ x: y, holdout: hold, forecast: models[key] });
        need(Number.isFinite(bt.metrics.rmse), 'non-finite forecast');
        return { key, ok: true, rmse: bt.metrics.rmse, mape: bt.metrics.mape, bias: bt.metrics.bias, skill: bt.skill, naive: bt.naive.rmse, forecast: bt.forecast };
      } catch (e) { return { key, ok: false, error: e.message }; }
    });
    const okTests = tests.filter((t) => t.ok).sort((a, b) => a.rmse - b.rmse);
    need(okTests.length > 0, `No time-series model could be fitted to the history: ${tests.map((t) => `${TS_MODELS[t.key]}: ${t.error}`).join('; ')}`);
    let chosen = okTests[0];
    if (inp.model !== 'auto') {
      const t = tests.find((x) => x.key === inp.model);
      need(t, `${TS_MODELS[inp.model]} needs seasonal data: at least ${2 * m + 4} observations at more than one observation per year.`);
      need(t.ok, `${TS_MODELS[inp.model]} could not be fitted to this history (${t.error}). Choose another model or "best on the hold-out test".`);
      chosen = t;
    }
    const fc = run(() => Array.from(models[chosen.key](y, h)), `${TS_MODELS[chosen.key]} forecast`);
    need(fc.every(Number.isFinite), `${TS_MODELS[chosen.key]} returned a non-finite forecast; choose another model.`);
    let lim = null; let arFit = null;
    try { arFit = fitArima(y); if (chosen.key === 'arima') { const f = arFit.forecast(h); lim = { lo: Array.from(f.lo95), hi: Array.from(f.hi95) }; } } catch (e) { arFit = null; }
    const naiveRmse = chosen.naive;

    /* ---- (B) stochastic process and percentile cone ---- */
    const ge = run(() => gbmEstimate(adj, dt), 'Volatility estimate');
    let se = null; try { se = schwartzEstimate(adj, dt); } catch (e) { se = null; }
    const mr = !!(se && se.meanReverting && Number.isFinite(se.kappa) && se.kappa > 0);
    let proc = inp.process === 'auto' ? (mr ? 'schwartz' : 'gbm') : inp.process;
    need(proc !== 'schwartz' || mr, 'The history does not show mean reversion (the fitted autoregressive slope is outside 0–1), so the Schwartz process cannot be estimated. Use the random walk or "auto".');
    const S0 = adj[n - 1]; const nPaths = Math.round(inp.nPaths);
    const sim = run(() => (proc === 'schwartz' ? schwartzSimulate({ S0, kappa: se.kappa, mu: se.mu, sigma: se.sigma, dt, steps: h, paths: nPaths, seed: inp.seed })
      : proc === 'jump' ? jumpDiffusionSimulate({ S0, mu: ge.mu, sigma: ge.sigma, lambda: inp.jumpRate, jumpMean: inp.jumpMean / 100, jumpStd: inp.jumpStd / 100, T: h * dt, steps: h, paths: nPaths, seed: inp.seed })
        : gbmSimulate({ S0, mu: ge.mu, sigma: ge.sigma, T: h * dt, steps: h, paths: nPaths, seed: inp.seed })), 'Price simulation');
    const cone = scenarioCone(sim.paths, [10, 50, 90]);
    const sAdd = (k) => (adjusted ? season(n - 1 + k) : 0); // step k of a path is period n − 1 + k
    const band = (key) => cone.bands[key].map((v, k) => v + sAdd(k));
    const b10 = band('p10'); const b50 = band('p50'); const b90 = band('p90'); const bMean = cone.mean.map((v, k) => v + sAdd(k));

    /* ---- revenue: contracted share at the contract price, remainder at the forecast price ---- */
    const share = inp.contractShare / 100; const vol = inp.volume * 1e3; // GWh → MWh
    const ny = Math.min(m, h);
    const blended = (price) => share * inp.contractPrice + (1 - share) * price;
    const revModel = vol * blended(sum(fc.slice(0, ny)) / ny);
    const yrAvg = Float64Array.from(sim.paths, (p) => { let s = 0; for (let k = 1; k <= ny; k++) s += p[k] + sAdd(k); return s / ny; });
    const PR = percentiles(yrAvg);
    const rev = (pp) => vol * blended(PR(pp));
    const p10Avg = sum(b10.slice(1)) / h;

    /* ---- (C) vector autoregression with the driver series ---- */
    let varFit = null; let varFc = null;
    const drv = (inp.driver?.y || []).map(Number);
    if (drv.length > 0) {
      need(drv.length === n && drv.every(Number.isFinite), `Driver series: it must have the same ${n} observations as the price history (it has ${drv.length}), or be left empty.`);
      varFit = run(() => var1Fit(adj.map((v, t) => [v, drv[t]])), 'Vector autoregression');
      varFc = varFit.forecast(h).map((row, k) => row[0] + sAdd(k + 1));
    }

    /* ---- (D) data-driven cost model ---- */
    const pr = tableRows(inp.projects, 'Completed projects', { mw: { min: 1e-9 }, depth: { min: 0 }, dist: { min: 0 }, cost: { min: 1e-9 } }, { min: 8 });
    uniqueNames(pr, 'Completed projects');
    const X = pr.map((r) => [r.mw, r.depth, r.dist]); const yc = pr.map((r) => r.cost);
    const logLog = inp.costForm === 'loglog';
    need(!logLog || X.every((r) => r.every((v) => v > 0)), 'Completed projects: the power-law regression needs positive capacity, depth and distance in every row; use the linear form otherwise.');
    const xNew = [inp.newMW, inp.newDepth, inp.newDist];
    const reg = run(() => regressionCostModel({ X, y: yc, logLog }), 'Cost regression');
    const meanCost = sum(yc) / yc.length;
    // Gaussian process: the observation noise is the scatter left by the regression (in log units for the power law,
    // so multiplied by the mean cost), otherwise the process would interpolate the noise of the records. Length scales
    // are fitted once on all records and held for the cross-validation folds.
    const gpNoise = Math.max((logLog ? reg.sigma * meanCost : reg.sigma) ** 2, 1e-10 * meanCost ** 2);
    let gpAll = null;
    const gpFit = (A, b) => {
      if (A === X) { gpAll = gpAll || gaussianProcess(X, yc, { kernel: 'matern52', noise: gpNoise }); return gpAll; }
      gpAll = gpAll || gaussianProcess(X, yc, { kernel: 'matern52', noise: gpNoise });
      return gaussianProcess(A, b, { kernel: 'matern52', noise: gpNoise, lengthScale: gpAll.params.lengthScale, optimise: false });
    };
    const fitters = {
      regression: (A, b) => { const f = regressionCostModel({ X: A, y: b, logLog }); return (x) => f.predict(x, 0.8).mean; },
      boosting: (A, b) => gradientBoosting({ X: A, y: b, trees: 100, learningRate: 0.1, maxDepth: 2 }).predict,
      forest: (A, b) => randomForest({ X: A, y: b, trees: 30, maxDepth: 4, features: 2, seed: inp.seed }).predict,
      gp: (A, b) => { const f = gpFit(A, b); return (x) => f.predict(x).mean; },
      mlp: (A, b) => mlpRegression({ X: A, y: b, hidden: 4, epochs: 300, seed: inp.seed }).predict,
    };
    const K = Math.min(Math.round(inp.folds), Math.floor(pr.length / 2));
    const costRes = Object.keys(fitters).map((key) => {
      try {
        const cv = new Array(pr.length).fill(NaN);
        for (let f = 0; f < K; f++) {
          const tr = []; const te = []; for (let i = 0; i < pr.length; i++) (i % K === f ? te : tr).push(i);
          const pred = fitters[key](tr.map((i) => X[i]), tr.map((i) => yc[i]));
          for (const i of te) cv[i] = pred(X[i]);
        }
        const em = errorMetrics(cv, yc);
        const full = fitters[key](X, yc)(xNew);
        need(Number.isFinite(em.rmse) && Number.isFinite(full), 'non-finite prediction');
        return { key, ok: true, rmse: em.rmse, mape: em.mape, bias: em.bias, pred: full, cv };
      } catch (e) { return { key, ok: false, error: e.message }; }
    });
    const regP = reg.predict(xNew, 0.8);
    const okCost = costRes.filter((c) => c.ok).sort((a, b) => a.rmse - b.rmse);
    need(okCost.length > 0, 'No cost model could be fitted to the project table.');
    let bestCost = okCost[0];
    if (inp.costModel !== 'auto') { const c = costRes.find((x) => x.key === inp.costModel); need(c.ok, `${COST_MODELS[inp.costModel]} could not be fitted to the project table (${c.error}).`); bestCost = c; }
    const lo = [0, 1, 2].map((j) => Math.min(...X.map((r) => r[j]))); const hi = [0, 1, 2].map((j) => Math.max(...X.map((r) => r[j])));
    const outside = ['capacity', 'water depth', 'distance'].filter((_, j) => xNew[j] < lo[j] || xNew[j] > hi[j]);
    // marginal cost of each driver from the regression, by central difference at the new project
    const marg = [0, 1, 2].map((j) => { const d = 0.01 * xNew[j]; const a = xNew.slice(); const b = xNew.slice(); a[j] -= d; b[j] += d; return (reg.predict(b, 0.8).mean - reg.predict(a, 0.8).mean) / (2 * d); });

    /* ---- (E) partial pooling of sparse benchmark groups ---- */
    const grp = tableRows(inp.groups, 'Cost growth by class', { mean: {}, variance: { min: 0 }, n: { min: 1 } }, { min: 2 });
    uniqueNames(grp, 'Cost growth by class');
    const pool = run(() => hierarchicalPooling(grp.map((g) => ({ name: g.name, mean: g.mean, variance: g.variance, n: g.n }))), 'Partial pooling');

    /* ---- (F) Poisson regression of maintenance demand on age ---- */
    const iv = tableRows(inp.interventions, 'Maintenance interventions', { age: { min: 0 }, count: { min: 0 } }, { min: 3 });
    const glm = run(() => glmFit({ X: iv.map((r) => [r.age]), y: iv.map((r) => r.count), family: 'poisson' }), 'Maintenance-demand model');
    need(glm.converged && glm.beta.every(Number.isFinite), 'The maintenance-demand model did not converge: check the intervention counts.');
    const lastAge = Math.max(...iv.map((r) => r.age));
    const nextAges = range(Math.round(inp.glmYears), (k) => lastAge + 1 + k);
    const nextCounts = nextAges.map((a) => glm.predict([a]));
    const ivNext = sum(nextCounts);

    /* ---- (G) carbon-price scenarios ---- */
    const cy = Math.round(inp.carbonYears); const dr = inp.discountRate / 100;
    const paths = carbonPricePaths({ years: cy, scenarios: { low: { start: inp.carbonStart, growth: inp.carbonLow / 100 }, central: { start: inp.carbonStart, growth: inp.carbonMid / 100 }, high: { start: inp.carbonStart, growth: inp.carbonHigh / 100 } } });
    const em = [0, ...new Array(cy).fill(inp.emissions)];
    const cPV = Object.fromEntries(['low', 'central', 'high'].map((k) => [k, carbonCost({ emissions_t: em, price: [paths[k][0], ...paths[k]], rate: dr }).pv]));

    const unitP = `${CUR}/MWh`;
    const perLabel = m === 1 ? 'year' : m === 4 ? 'quarter' : 'month';
    const kpis = [
      { key: 'fcNext', label: `Forecast for the next ${perLabel}`, value: round(fc[0]), unit: unitP, explain: `Central forecast of ${TS_MODELS[chosen.key]} one period ahead, in ${inp.basis} terms. The last observation was ${fmt(y[n - 1])}.` },
      { key: 'fcEnd', label: `Forecast at the horizon (${h} periods)`, value: round(fc[h - 1]), unit: unitP, explain: `Central forecast ${fmt(h * dt)} years ahead from the chosen model. Its reliability falls with distance: compare with the percentile cone, which widens with the horizon.` },
      { key: 'rmse', label: 'Hold-out forecast error (RMSE)', value: round(chosen.rmse), unit: unitP, explain: `Root-mean-square error of the chosen model on the last ${hold} observations, which it did not see. The honest measure of accuracy; in-sample fit is always better.` },
      { key: 'mape', label: 'Hold-out mean absolute percentage error', value: round(100 * chosen.mape), unit: '%', explain: 'Average size of the hold-out error relative to the observation. Under about 10 % is useful for budgeting; revenue cases need the distribution, not this average.' },
      { key: 'skill', label: 'Skill against the naive forecast', value: round(chosen.skill), unit: '-', limit: 0, sense: 'min', explain: 'One minus the ratio of the model error to the error of "no change from the last value". Zero or negative means the model adds nothing over the naive forecast and should not be relied on.' },
      { key: 'volAnnual', label: 'Annual price volatility', value: round(100 * ge.sigma), unit: '%/√yr', explain: `Standard deviation of the logarithmic price change, annualised${adjusted && sIdx ? ' after removing the seasonal pattern' : ''}. The input the option and risk analyses need; standard error ${pct(ge.seSigma)}.` },
      { key: 'driftAnnual', label: 'Annual price drift', value: round(100 * ge.mu), unit: '%/yr', explain: `Mean growth rate of the price over the history in ${inp.basis} terms. With only ${fmt(n * dt)} years of data it is poorly determined: a trend seen in the past is weak evidence of a trend to come.` },
    ];
    if (mr) kpis.push(
      { key: 'halfLife', label: 'Half-life of price shocks', value: round(se.halfLife), unit: 'yr', explain: 'Time for a deviation from the long-run level to halve under the mean-reverting model. Short half-lives mean today’s price says little about prices a few years out.' },
      { key: 'longRun', label: 'Long-run price level', value: round(Math.exp(se.alpha)), unit: unitP, explain: 'Median level the mean-reverting process pulls towards. If it sits below the break-even price, a currently high price is not a basis for investment.' });
    kpis.push(
      { key: 'meanEnd', label: 'Simulated mean price at the horizon', value: round(bMean[h]), unit: unitP, explain: `Mean of ${nPaths} simulated paths of the ${proc === 'schwartz' ? 'mean-reverting' : proc === 'jump' ? 'jump-diffusion' : 'random-walk'} process at the horizon. Its sampling error falls with the square root of the number of paths.` },
      { key: 'p10End', label: 'Price P10 at the horizon (cumulative, downside)', value: round(b10[h]), unit: unitP, explain: 'Price not exceeded with 10 % probability at the horizon: the low case. In the exceedance convention of resource estimates the same number is called P90.' },
      { key: 'p50End', label: 'Price P50 at the horizon (median)', value: round(b50[h]), unit: unitP, explain: 'Price with an even chance of being exceeded at the horizon.' },
      { key: 'p90End', label: 'Price P90 at the horizon (cumulative, upside)', value: round(b90[h]), unit: unitP, explain: 'Price not exceeded with 90 % probability at the horizon: the high case (exceedance P10).' },
      { key: 'revNext', label: 'Expected revenue, next year', value: mil(revModel), unit: `${MU}/yr`, explain: `Volume times the blend of ${inp.contractShare} % at the contract price and the rest at the model forecast averaged over the next ${ny} period(s).` },
      { key: 'revP10', label: 'Revenue P10, next year (cumulative, downside)', value: mil(rev(10)), unit: `${MU}/yr`, explain: 'Annual revenue not exceeded with 10 % probability given the simulated price: the figure a lender sizes debt on. The contracted share is what holds it up.' },
      { key: 'revP90', label: 'Revenue P90, next year (cumulative, upside)', value: mil(rev(90)), unit: `${MU}/yr`, explain: 'Annual revenue not exceeded with 90 % probability: the high case. The spread between P10 and P90 is the merchant exposure.' },
    );
    if (varFc) kpis.push({ key: 'varEnd', label: 'Forecast at the horizon from the vector autoregression', value: round(varFc[h - 1]), unit: unitP, explain: `Price forecast when the driver series is modelled jointly with it. The driver’s lagged effect on the price is ${fmt(varFit.A[0][1])} per unit; a value near zero means the driver adds no information.` });
    kpis.push(
      { key: 'costPred', label: `Predicted cost of the new project (${COST_MODELS[bestCost.key]})`, value: mil(bestCost.pred), unit: MU, explain: `Estimate from the ${inp.costModel === 'auto' ? 'model with the lowest cross-validated error' : 'selected model'} for ${fmt(inp.newMW)} MW in ${fmt(inp.newDepth)} m of water, ${fmt(inp.newDist)} km from shore. An early (Class 5) figure: it knows only what the comparables share with this project.` },
      { key: 'costReg', label: 'Regression estimate', value: mil(regP.mean), unit: MU, explain: `Multivariate ${logLog ? 'power-law' : 'linear'} regression on capacity, depth and distance (R² ${fmt(reg.r2)}). The transparent model: every coefficient can be read and challenged.` },
      { key: 'costRegP10', label: 'Regression P10 (cumulative)', value: mil(regP.piLow), unit: MU, explain: 'Lower end of the 80 % prediction interval for one new project: the cost not exceeded with 10 % probability given the scatter of the comparables.' },
      { key: 'costRegP90', label: 'Regression P90 (cumulative)', value: mil(regP.piHigh), unit: MU, explain: 'Upper end of the 80 % prediction interval: the cost not exceeded with 90 % probability. Scope differences from the comparables are not inside this band.' },
      { key: 'cvError', label: 'Cross-validated error of the reported cost model', value: round(100 * bestCost.rmse / meanCost), unit: '% of mean cost', explain: `Out-of-sample RMSE over ${K} folds relative to the mean project cost. This, not the fit to the training data, is the accuracy to expect on a new project inside the data range.` },
      { key: 'grandMean', label: 'Pooled mean cost growth across classes', value: round(pool.grandMean), unit: '%', explain: 'Precision-weighted mean of the class means in the hierarchical model: the best estimate for a class with no data of its own.' },
      { key: 'tau', label: 'Between-class standard deviation', value: round(Math.sqrt(pool.tau2)), unit: '%', explain: 'How much the true class means differ from each other. Zero means the classes are indistinguishable and the pooled mean applies to all; large means each class keeps its own estimate.' },
      { key: 'interventionsNext', label: `Expected interventions over the next ${Math.round(inp.glmYears)} years`, value: round(ivNext), unit: '-', explain: `Poisson regression of interventions on age, projected from age ${fmt(lastAge + 1)}. Demand changes by ${pct(Math.exp(glm.beta[1]) - 1)} per year of age.` },
      { key: 'interventionCost', label: 'Cost of that maintenance demand', value: mil(ivNext * inp.interventionCost), unit: MU, explain: 'Projected interventions times the cost per intervention, undiscounted: the maintenance budget implied by the ageing trend in the records.' },
      { key: 'carbonPV', label: 'Present value of carbon cost, central scenario', value: mil(cPV.central), unit: MU, explain: `Emissions of ${fmt(inp.emissions)} t CO₂e/yr priced along the central path over ${cy} years and discounted at ${inp.discountRate} %. Low and high scenarios give ${fM(cPV.low)} and ${fM(cPV.high)}.` },
    );

    const rowsPerDriver = pr.length / 3;
    const checks = [
      check('Forecast error ≤ naive benchmark', chosen.rmse, Math.max(naiveRmse, 1e-12), unitP, 'Hyndman & Athanasopoulos, Forecasting: Principles and Practice, §5.2 (benchmark methods)',
        'Guards against relying on a model that is no better than assuming no change. Above 1.0 the model has no demonstrated skill on unseen data and its forecast should not carry a decision.'),
      checkMin('Downside price (cumulative P10, horizon average) ≥ break-even price', blended(p10Avg), Math.max(inp.breakEven, 1e-9), unitP, 'Lender downside case on realised price (project-finance practice)',
        'Guards against a revenue stream that fails to cover cost in the low-price case. Above 1.0 the blended realised price in the P10 case is below break-even: raise the contracted share or lower the cost.'),
      checkMin('Project records per cost driver ≥ 10', rowsPerDriver, 10, 'rows per driver', 'ISPA/SCEA Parametric Estimating Handbook, 4th ed., ch. 3 (sample size for cost-estimating relationships)',
        'Guards against fitting flexible models to too little data, which reproduces noise rather than cost behaviour. Above 1.0 only the regression, with its prediction interval, should be used.'),
    ];

    const tHist = range(n); const tFc = range(h, (k) => n + k); const tCone = range(h + 1, (k) => n - 1 + k);
    const tail = Math.max(0, n - Math.max(4 * m, 24));
    const plots = [
      { type: 'line', title: 'History and forecast', x: { label: 'Period', unit: '-' }, y: { label: 'Price', unit: unitP },
        series: [{ name: 'History', x: tHist, y: y.map((v) => round(v)) }, { name: `Forecast (${TS_MODELS[chosen.key]})`, x: tFc, y: fc.map((v) => round(v)) }, { name: 'Hold-out test', x: range(hold, (k) => n - hold + k), y: chosen.forecast.map((v) => round(v)), dash: true },
          ...(lim ? [{ name: '95 % lower', x: tFc, y: lim.lo.map((v) => round(v)), dash: true }, { name: '95 % upper', x: tFc, y: lim.hi.map((v) => round(v)), dash: true }] : []),
          ...(varFc ? [{ name: 'Vector autoregression', x: tFc, y: varFc.map((v) => round(v)), dash: true }] : [])],
        explain: 'The observed series, the forecast and, over the last observations, what the model would have predicted without seeing them. How far the dashed hold-out line strays from the history is the error to expect.' },
      { type: 'line', title: 'Stochastic price cone', x: { label: 'Period', unit: '-' }, y: { label: 'Price', unit: unitP },
        series: [{ name: 'History', x: tHist.slice(tail), y: y.slice(tail).map((v) => round(v)) }, { name: 'P90 (cumulative)', x: tCone, y: b90.map((v) => round(v)), dash: true }, { name: 'P50', x: tCone, y: b50.map((v) => round(v)) }, { name: 'P10 (cumulative)', x: tCone, y: b10.map((v) => round(v)), dash: true },
          { name: 'Break-even', x: [tHist[tail], n - 1 + h], y: [inp.breakEven, inp.breakEven], dash: true }],
        explain: 'Simulated percentiles of the price (cumulative: P10 is the low line). A random walk widens without limit; a mean-reverting price settles into a band. Look for where the P10 line crosses the break-even price.' },
      { type: 'bar', title: 'Hold-out error by forecast model', categories: [...okTests.map((t) => TS_MODELS[t.key]), 'Naive (no change)'], x: { label: 'Model' }, y: { label: 'RMSE', unit: unitP }, series: [{ name: 'RMSE', y: [...okTests.map((t) => round(t.rmse)), round(naiveRmse)] }],
        explain: 'Out-of-sample error of every model that could be fitted, lowest first, beside the naive benchmark. A model taller than the last bar has no skill on this series.' },
      { type: 'scatter', title: 'Cost model: predicted against actual (cross-validated)', x: { label: 'Actual cost', unit: MU }, y: { label: 'Predicted cost', unit: MU },
        series: [{ name: COST_MODELS[bestCost.key], x: yc.map(mil), y: bestCost.cv.map(mil) }, { name: 'Perfect prediction', x: [mil(Math.min(...yc)), mil(Math.max(...yc))], y: [mil(Math.min(...yc)), mil(Math.max(...yc))], dash: true }],
        explain: 'Each completed project predicted by a model that was not shown it. Points on the line are predicted exactly; the vertical scatter is the error to expect, and a tilt away from the line is bias at small or large projects.' },
      { type: 'bar', title: 'Cross-validated error by cost model', categories: okCost.map((c) => COST_MODELS[c.key]), x: { label: 'Model' }, y: { label: 'RMSE', unit: '% of mean cost' }, series: [{ name: 'RMSE', y: okCost.map((c) => round(100 * c.rmse / meanCost)) }],
        explain: 'Out-of-sample error of each cost model. Where the flexible models do no better than the regression, the regression is preferred because it can be explained and audited.' },
      { type: 'line', title: 'Maintenance demand against asset age', x: { label: 'Asset age', unit: 'yr' }, y: { label: 'Interventions', unit: '1/yr' },
        series: [{ name: 'Recorded', x: iv.map((r) => r.age), y: iv.map((r) => r.count) }, { name: 'Poisson regression', x: [...iv.map((r) => r.age), ...nextAges], y: [...iv.map((r) => round(glm.predict([r.age]))), ...nextCounts.map((v) => round(v))], dash: true }],
        explain: 'Recorded interventions and the fitted ageing trend extended over the projection years. A rising curve is the signature of wear-out; a flat one means age does not explain the demand.' },
      { type: 'line', title: 'Carbon-price scenarios', x: { label: 'Year', unit: '-' }, y: { label: 'Carbon price', unit: `${CUR}/t CO₂e` },
        series: ['low', 'central', 'high'].map((k) => ({ name: k[0].toUpperCase() + k.slice(1), x: range(cy, (j) => j + 1), y: paths[k].map((v) => round(v)) })),
        explain: 'The three assumed carbon-price paths. They are scenarios, not forecasts: the spread between them is the policy risk carried by emitting activities.' },
    ];
    if (sIdx) plots.push({ type: 'bar', title: 'Seasonal pattern', categories: range(m, (k) => (m === 4 ? `Q${k + 1}` : `M${k + 1}`)), x: { label: 'Position in the year' }, y: { label: 'Seasonal component', unit: unitP }, series: [{ name: 'Seasonal index', y: sIdx.map((v) => round(v)) }],
      explain: 'Average deviation of each position in the year from the trend (additive). Revenue timing and maintenance outages should be planned against this pattern.' });

    const tables = [
      { title: 'Forecast model comparison (hold-out test)', columns: [{ key: 'model', label: 'Model' }, { key: 'rmse', label: 'RMSE', unit: unitP }, { key: 'mape', label: 'MAPE', unit: '%' }, { key: 'bias', label: 'Bias', unit: unitP }, { key: 'skill', label: 'Skill vs naive' }, { key: 'note', label: 'Note' }],
        rows: tests.map((t) => (t.ok ? { model: TS_MODELS[t.key], rmse: round(t.rmse, 5), mape: round(100 * t.mape, 4), bias: round(t.bias, 4), skill: round(t.skill, 4), note: t.key === chosen.key ? 'used' : '' } : { model: TS_MODELS[t.key], rmse: '–', mape: '–', bias: '–', skill: '–', note: `not fitted: ${t.error}` })),
        explain: 'Every candidate scored on the same withheld observations. Bias is the mean of forecast minus actual: a large bias with a small RMSE points to a missed trend.' },
      { title: 'Forecast and percentile cone', columns: [{ key: 't', label: 'Period' }, { key: 'year', label: 'Year' }, { key: 'fc', label: 'Model forecast', unit: unitP }, { key: 'p10', label: 'P10', unit: unitP }, { key: 'p50', label: 'P50', unit: unitP }, { key: 'p90', label: 'P90', unit: unitP }, { key: 'mean', label: 'Simulated mean', unit: unitP }],
        rows: fc.map((v, k) => ({ t: n + k, year: round(inp.firstYear + (n + k) * dt, 6), fc: round(v, 5), p10: round(b10[k + 1], 5), p50: round(b50[k + 1], 5), p90: round(b90[k + 1], 5), mean: round(bMean[k + 1], 5) })),
        explain: 'The central forecast of the time-series model beside the percentiles of the simulated process, period by period. Percentiles are cumulative: P10 is the low value and P90 the high value.' },
      { title: 'Estimated process parameters', columns: [{ key: 'item', label: 'Parameter' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
        rows: [{ item: 'Volatility (log returns)', value: round(100 * ge.sigma, 5), unit: '%/√yr' }, { item: 'Drift', value: round(100 * ge.mu, 5), unit: '%/yr' },
          { item: 'Mean reversion found', value: mr ? 'yes' : 'no', unit: '' }, ...(mr ? [{ item: 'Reversion speed κ', value: round(se.kappa, 5), unit: '1/yr' }, { item: 'Half-life', value: round(se.halfLife, 5), unit: 'yr' }, { item: 'Long-run median', value: round(Math.exp(se.alpha), 5), unit: unitP }, { item: 'Volatility of the mean-reverting process', value: round(100 * se.sigma, 5), unit: '%/√yr' }] : []),
          ...(arFit ? [{ item: `ARIMA(${inp.arP},${inp.arD},${inp.arQ}) AR coefficients`, value: arFit.phi.map((v) => fmt(v, 4)).join(', ') || '–', unit: '' }, { item: 'ARIMA MA coefficients', value: arFit.theta.map((v) => fmt(v, 4)).join(', ') || '–', unit: '' }, { item: 'ARIMA drift per period', value: round(arFit.mean, 5), unit: unitP }, { item: 'ARIMA stationary', value: arFit.stationary ? 'yes' : 'no', unit: '' }] : []),
          ...(varFit ? [{ item: 'VAR: price on lagged price', value: round(varFit.A[0][0], 5), unit: '' }, { item: 'VAR: price on lagged driver', value: round(varFit.A[0][1], 5), unit: '' }, { item: 'VAR: driver on lagged price', value: round(varFit.A[1][0], 5), unit: '' }, { item: 'VAR: driver on lagged driver', value: round(varFit.A[1][1], 5), unit: '' }] : []),
          { item: 'Process simulated', value: proc === 'schwartz' ? 'Schwartz mean reversion' : proc === 'jump' ? 'Jump diffusion' : 'Geometric Brownian motion', unit: '' }],
        explain: 'The fitted parameters behind the forecast and the cone, for audit and for use in the option and risk modules.' },
      { title: 'Cost models compared', columns: [{ key: 'model', label: 'Model' }, { key: 'rmse', label: 'CV RMSE', unit: '% of mean' }, { key: 'mape', label: 'CV MAPE', unit: '%' }, { key: 'pred', label: 'Prediction for the new project', unit: MU }, { key: 'note', label: 'Note' }],
        rows: costRes.map((c) => (c.ok ? { model: COST_MODELS[c.key], rmse: round(100 * c.rmse / meanCost, 4), mape: round(100 * c.mape, 4), pred: mil(c.pred), note: c.key === bestCost.key ? 'reported' : '' } : { model: COST_MODELS[c.key], rmse: '–', mape: '–', pred: '–', note: `not fitted: ${c.error}` })),
        explain: 'Each model with its out-of-sample error and its estimate for the new project. Disagreement between models of similar error is itself a measure of the uncertainty of an early estimate.' },
      { title: 'Cost growth by class: raw and pooled', columns: [{ key: 'name', label: 'Class' }, { key: 'n', label: 'Projects' }, { key: 'raw', label: 'Observed mean', unit: '%' }, { key: 'pooled', label: 'Pooled estimate', unit: '%' }, { key: 'shrink', label: 'Weight on the common mean', unit: '%' }, { key: 'std', label: 'Standard error', unit: '%' }],
        rows: pool.groups.map((g, k) => ({ name: g.name, n: grp[k].n, raw: round(g.raw, 5), pooled: round(g.pooled, 5), shrink: round(100 * g.shrinkage, 4), std: round(g.std, 4) })),
        explain: 'Each class mean pulled towards the common mean in proportion to how little data it rests on. Use the pooled column, not the raw one, as the benchmark for a class with few projects.' },
    ];

    const span = rev(90) - rev(10);
    const advice = [
      adv('economic', chosen.skill > 0 && blended(p10Avg) >= inp.breakEven ? 'ok' : 'warn',
        `${TS_MODELS[chosen.key]} forecasts ${fmt(fc[h - 1])} ${unitP} at ${fmt(h * dt)} years with a hold-out error of ${pct(chosen.mape)} (skill ${fmt(chosen.skill)} against the naive forecast). Simulated price at the horizon, cumulative P10 / P50 / P90: ${fmt(b10[h])} / ${fmt(b50[h])} / ${fmt(b90[h])}. Next-year revenue ${fM(revModel)} expected, ${fM(rev(10))} at P10 and ${fM(rev(90))} at P90.`,
        chosen.skill <= 0 ? 'The model does not beat "no change": use the last value with the percentile cone, and treat any trend as a scenario rather than a forecast.'
          : blended(p10Avg) < inp.breakEven ? `The blended P10 price of ${fmt(blended(p10Avg))} is below the break-even of ${fmt(inp.breakEven)}: the contracted share would need to be ${inp.contractPrice > p10Avg ? pct(Math.min(1, Math.max(0, (inp.breakEven - p10Avg) / (inp.contractPrice - p10Avg)))) : 'above 100 %'} at ${fmt(inp.contractPrice)} ${unitP} to cover it.`
            : `Downside price covers break-even with ${inp.contractShare} % contracted. The ${fM(span)} P10–P90 revenue spread is the merchant exposure: each 10 points of contracted share removes about ${fM(0.1 * vol * (PR(90) - PR(10)))} of it.`,
        `Out-of-sample evaluation, Hyndman & Athanasopoulos (2021) ch. 5; lender downside price cases in project finance (Yescombe 2014, ch. 9); escalation and price forecasts per AACE RP 58R-10. ${PCT_CUM}`),
      adv('engineering', outside.length || rowsPerDriver < 10 ? 'warn' : 'info',
        `The cost of the new project is driven by ${fmt(inp.newMW)} MW, ${fmt(inp.newDepth)} m depth and ${fmt(inp.newDist)} km distance: at this point one more megawatt adds ${fM(marg[0])}, one more metre of water ${fM(marg[1])} and one more kilometre from shore ${fM(marg[2])}. Estimate ${fM(bestCost.pred)} (${COST_MODELS[bestCost.key]}, cross-validated error ${pct(bestCost.rmse / meanCost)}); regression 80 % interval ${fM(regP.piLow)} to ${fM(regP.piHigh)}.${outside.length ? ` The new project lies outside the data range in ${outside.join(' and ')}.` : ''}`,
        `${outside.length ? 'Do not use the tree, Gaussian-process or neural-network estimate outside the data range: they cannot extrapolate. Use the regression and widen its interval. ' : ''}Firm up the driver the estimate is most sensitive to with its engineering source (bathymetry and foundation concept for depth, the cable route for distance), then replace this early figure with the quantity-based estimate. Maintenance demand is rising ${pct(Math.exp(glm.beta[1]) - 1)} per year of age: ${fmt(ivNext)} interventions over the next ${Math.round(inp.glmYears)} years.`,
        'Cost-estimating relationships and their validation, ISPA/SCEA Parametric Estimating Handbook, 4th ed., ch. 3; AACE RP 18R-97 (Class 5 and 4 estimates by parametric models).',
        `Operational emissions of ${fmt(inp.emissions)} t CO₂e/yr carry a carbon cost with a present value of ${fM(cPV.low)} to ${fM(cPV.high)} across the scenarios (${fM(cPV.central)} central) over ${cy} years. Each tonne per year abated is worth about ${fmt(inp.emissions > 0 ? cPV.central / inp.emissions : 0)} ${CUR} in present value on the central path, which is the budget available for lower-carbon vessels and logistics.`),
    ];
    if (rowsPerDriver < 10) advice.push(adv('economic', 'warn', `Only ${pr.length} project records support three cost drivers (${fmt(rowsPerDriver)} per driver).`, 'Rely on the regression and its prediction interval; treat the machine-learning estimates as a cross-check until at least 30 comparable records are available.', 'ISPA/SCEA Parametric Estimating Handbook, 4th ed., ch. 3 (data sufficiency).'));

    const notes = [PCT_CUM, `${inp.basis === 'real' ? 'Real' : 'Nominal'} terms, as the history. ${n} observations at ${m} per year from ${inp.firstYear}; forecast ${h} periods; hold-out ${hold} periods.`,
      `Forecast model: ${TS_MODELS[chosen.key]}${inp.model === 'auto' ? ' (lowest hold-out error)' : ' (selected)'}. Price process: ${proc === 'schwartz' ? 'mean reverting' : proc === 'jump' ? 'jump diffusion' : 'geometric Brownian motion'}, ${nPaths} paths, seed ${inp.seed}.`,
      sIdx ? (adjusted ? 'Seasonality was removed before estimating volatility and added back to the cone.' : 'The seasonal pattern is too large relative to the price level to remove additively; volatility includes it.') : 'No seasonal component (annual data or too short a history).'];
    return { kpis, checks, plots, tables, advice, notes,
      publish: { 'econ.priceForecast': fc[h - 1], 'econ.priceP10': b10[h], 'econ.priceP90': b90[h], 'econ.priceVolatility': ge.sigma, 'econ.revenueForecast': revModel, 'econ.capexDataDriven': bestCost.pred } };
  },
  verify: [
    // Holt's linear method reproduces a straight line exactly: y = 50 + 2t for t = 0…23, so 6 steps ahead is 50 + 2 × 29 = 108
    { name: 'Linear trend is extrapolated exactly', source: 'Limiting case: Holt’s linear exponential smoothing has zero error on a noise-free straight line, forecast = 50 + 2 (23 + 6) = 108',
      inputs: { perYear: '1', history: { x: range(24), y: range(24, (t) => 50 + 2 * t) }, driver: { x: [], y: [] }, model: 'holt', horizon: 6, holdout: 4 },
      expect: [{ kpi: 'fcEnd', value: 108, tol: 1e-4 }, { kpi: 'fcNext', value: 98, tol: 1e-4 }] },
    // 24 annual log returns alternating +0.1 and −0.1: mean 0, sample standard deviation 0.1 √(24/23) = 0.1021508
    { name: 'Volatility from log returns', source: 'Hull, Options, Futures and Other Derivatives, ch. 15: σ = sd(ln S_t+1/S_t)/√Δt = 0.1 √(24/23) = 10.21508 % for returns alternating ±0.1 at Δt = 1',
      inputs: { perYear: '1', history: { x: range(25), y: range(25, (t) => 100 * Math.exp(t % 2 ? 0.1 : 0)) }, driver: { x: [], y: [] }, holdout: 4, horizon: 4 },
      expect: [{ kpi: 'volAnnual', value: 10.21508, tol: 1e-5 }] },
    // exact linear cost surface: C = 100 M + 2 M per MW + 3 M per m + 1 M per km → 100 + 1000 + 120 + 60 = 1280 M at (500, 40, 60)
    { name: 'Regression recovers an exact linear cost relationship', source: 'Hand calculation: C = 100 + 2 × 500 + 3 × 40 + 1 × 60 = 1280 M',
      inputs: { costForm: 'linear', costModel: 'regression', newMW: 500, newDepth: 40, newDist: 60,
        projects: [[200, 20, 30], [300, 35, 80], [450, 25, 20], [600, 50, 110], [750, 30, 60], [900, 45, 40], [350, 55, 95], [520, 18, 70], [680, 42, 15], [820, 28, 100], [260, 48, 50], [980, 38, 85]].map(([mw, depth, dist], k) => ({ name: `P${k}`, mw, depth, dist, cost: 1e8 + 2e6 * mw + 3e6 * depth + 1e6 * dist })) },
      expect: [{ kpi: 'costReg', value: 1280, tol: 1e-6 }, { kpi: 'costPred', value: 1280, tol: 1e-6 }] },
    // three classes with standard error² = 4/4 = 1: fixed-effect mean 12, Q = 4 + 4 + 0 = 8, τ² = (Q − 2)/(3 − 1) = 3
    { name: 'Partial pooling of three equal-precision classes', source: 'DerSimonian & Laird (1986) moment estimator by hand: means 10, 14, 12 with s²/n = 1 give mean 12, τ² = (8 − 2)/2 = 3, τ = 1.7320508',
      inputs: { groups: [{ name: 'a', mean: 10, variance: 4, n: 4 }, { name: 'b', mean: 14, variance: 4, n: 4 }, { name: 'c', mean: 12, variance: 4, n: 4 }] },
      expect: [{ kpi: 'grandMean', value: 12, tol: 1e-9 }, { kpi: 'tau', value: 1.7320508, tol: 1e-6 }] },
    // counts exactly exp(0.5 + 0.3 age) for ages 0…5: the next two years give e^2.3 + e^2.6 = 9.97418 + 13.46374 = 23.43792
    { name: 'Poisson regression recovers an exact log-linear trend', source: 'Maximum likelihood on noise-free data returns the generating coefficients: e^(0.5 + 0.3 × 6) + e^(0.5 + 0.3 × 7) = 23.43792',
      inputs: { glmYears: 2, interventions: range(6, (a) => ({ age: a, count: Math.exp(0.5 + 0.3 * a) })) }, expect: [{ kpi: 'interventionsNext', value: 23.43792, tol: 1e-5 }] },
  ],
};
