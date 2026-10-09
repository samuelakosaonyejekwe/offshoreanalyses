// Price and cost forecasting: stochastic processes (GBM, Ornstein–Uhlenbeck / Schwartz one-factor, jump
// diffusion) with seeded simulation and parameter estimation, scenario cones, seasonal decomposition, wrappers
// over the shared time-series library with hold-out back-testing, regression cost models with intervals,
// generalised linear models, vector autoregression, Kalman local-level filter, Gaussian-process and tree-ensemble
// regressors, a small neural network, Bayesian updating, reference-class uplift, escalation, exchange-rate and
// carbon-price scenarios.
// Data-driven models are only as good as the data given to them: every fit returns its error measures, and the
// back-test function should be used before a forecast is relied upon.

import { makeRng } from '../math/rng.js';
import { betaI } from '../math/special.js';
import { solveDense, inverse } from '../math/linalg.js';
import { expSmoothing, arima, kalmanFilter } from '../math/timeseries.js';
import { gaussianProcess } from '../reliability/surrogate.js';
import { requireNumber, requireSeries } from './tvm.js';
import { normInv, quantile } from './risk.js';
import { costAtCumulativeVolume, AACE_CLASSES } from './cost.js';

// ---------------------------------------------------------------------------------------------------------------
// Stochastic processes
// ---------------------------------------------------------------------------------------------------------------

/**
 * Geometric Brownian motion paths by the exact discretisation S_(t+Δ) = S_t exp[(μ − σ²/2)Δ + σ√Δ z]. Seeded.
 * Ref: Hull, Options, Futures and Other Derivatives, ch. 14; Dixit & Pindyck (1994), Investment under
 * Uncertainty, ch. 3.
 * @param {{S0: number, mu: number, sigma: number, T: number, steps: number, paths?: number, seed?: number}} o start
 *   value (currency), drift (1/yr), volatility (1/√yr), horizon (yr), number of steps and paths
 * @returns {{t: number[], paths: Float64Array[]}} each path has steps + 1 values
 */
export function gbmSimulate({ S0, mu, sigma, T, steps, paths = 1000, seed = 1 }) {
  requireNumber('S0', S0, { gt: 0 }); requireNumber('mu', mu); requireNumber('sigma', sigma, { min: 0 });
  requireNumber('T', T, { gt: 0 }); requireNumber('steps', steps, { min: 1, integer: true }); requireNumber('paths', paths, { min: 1, integer: true });
  const dt = T / steps; const drift = (mu - 0.5 * sigma * sigma) * dt; const vol = sigma * Math.sqrt(dt);
  const rng = makeRng(seed); const out = [];
  for (let p = 0; p < paths; p++) {
    const a = new Float64Array(steps + 1); a[0] = S0;
    for (let k = 1; k <= steps; k++) a[k] = a[k - 1] * Math.exp(drift + vol * rng.normal());
    out.push(a);
  }
  return { t: Array.from({ length: steps + 1 }, (_, k) => k * dt), paths: out };
}

/**
 * GBM parameters from a price history by the moments of log returns: σ = sd(r)/√Δ, μ = mean(r)/Δ + σ²/2.
 * Ref: Hull, Options, Futures and Other Derivatives, ch. 15 (estimating volatility from historical data).
 * @param {number[]} prices positive prices at equal spacing
 * @param {number} dt spacing in years (e.g. 1/12 for monthly)
 * @returns {{mu: number, sigma: number, seSigma: number, n: number}} per year; seSigma ≈ σ/√(2n)
 */
export function gbmEstimate(prices, dt) {
  requireSeries('prices', prices, { minLength: 3 }); requireNumber('dt', dt, { gt: 0 });
  const n = prices.length - 1; let m = 0; const r = new Float64Array(n);
  for (let k = 0; k < n; k++) { if (!(prices[k] > 0 && prices[k + 1] > 0)) throw new RangeError('prices must be positive'); r[k] = Math.log(prices[k + 1] / prices[k]); m += r[k]; }
  m /= n;
  let v = 0; for (let k = 0; k < n; k++) v += (r[k] - m) ** 2;
  const sigma = Math.sqrt(v / (n - 1) / dt);
  return { mu: m / dt + 0.5 * sigma * sigma, sigma, seSigma: sigma / Math.sqrt(2 * n), n };
}

/**
 * Ornstein–Uhlenbeck (mean-reverting) paths dx = κ(θ − x)dt + σ dW by the exact transition
 * x_(t+Δ) = θ + (x_t − θ)e^(−κΔ) + σ √[(1 − e^(−2κΔ))/(2κ)] z. Seeded.
 * Ref: Uhlenbeck, G.E. & Ornstein, L.S. (1930), Physical Review 36; Dixit & Pindyck (1994), ch. 3.
 * @param {{x0: number, kappa: number, theta: number, sigma: number, dt: number, steps: number, paths?: number,
 *   seed?: number}} o start, reversion speed (1/yr), long-run mean, volatility (unit/√yr), step (yr)
 * @returns {{t: number[], paths: Float64Array[]}}
 */
export function ouSimulate({ x0, kappa, theta, sigma, dt, steps, paths = 1, seed = 1 }) {
  requireNumber('x0', x0); requireNumber('kappa', kappa, { gt: 0 }); requireNumber('theta', theta);
  requireNumber('sigma', sigma, { min: 0 }); requireNumber('dt', dt, { gt: 0 }); requireNumber('steps', steps, { min: 1, integer: true });
  const e = Math.exp(-kappa * dt); const sd = sigma * Math.sqrt((1 - e * e) / (2 * kappa));
  const rng = makeRng(seed); const out = [];
  for (let p = 0; p < paths; p++) {
    const a = new Float64Array(steps + 1); a[0] = x0;
    for (let k = 1; k <= steps; k++) a[k] = theta + (a[k - 1] - theta) * e + sd * rng.normal();
    out.push(a);
  }
  return { t: Array.from({ length: steps + 1 }, (_, k) => k * dt), paths: out };
}

/**
 * OU parameters from an equally spaced history by AR(1) regression x_(t+1) = a + b x_t + ε with the exact
 * discretisation: κ = −ln b / Δ, θ = a/(1 − b), σ = sd(ε) √[2κ/(1 − b²)].
 * Ref: Dixit & Pindyck (1994), Investment under Uncertainty, ch. 3 (estimation of a mean-reverting process).
 * @param {number[]} x observations at equal spacing
 * @param {number} dt spacing (yr)
 * @returns {{kappa: number, theta: number, sigma: number, halfLife: number, b: number, meanReverting: boolean,
 *   residualStd: number}} meanReverting is false (κ, θ, σ NaN) when the AR(1) slope is outside (0, 1)
 */
export function ouEstimate(x, dt) {
  requireSeries('x', x, { minLength: 4 }); requireNumber('dt', dt, { gt: 0 });
  const n = x.length - 1; let sx = 0; let sy = 0; let sxx = 0; let sxy = 0;
  for (let k = 0; k < n; k++) { sx += x[k]; sy += x[k + 1]; sxx += x[k] * x[k]; sxy += x[k] * x[k + 1]; }
  const b = (n * sxy - sx * sy) / (n * sxx - sx * sx); const a = (sy - b * sx) / n;
  let sse = 0; for (let k = 0; k < n; k++) sse += (x[k + 1] - a - b * x[k]) ** 2;
  const sd = Math.sqrt(sse / (n - 2));
  if (!(b > 0 && b < 1)) return { kappa: NaN, theta: NaN, sigma: NaN, halfLife: NaN, b, meanReverting: false, residualStd: sd };
  const kappa = -Math.log(b) / dt;
  return { kappa, theta: a / (1 - b), sigma: sd * Math.sqrt(2 * kappa / (1 - b * b)), halfLife: Math.LN2 / kappa, b,
    meanReverting: true, residualStd: sd };
}

/**
 * Schwartz one-factor commodity price model dS = κ(μ − ln S)S dt + σS dW: ln S is an OU process with long-run
 * mean α = μ − σ²/(2κ); simulated exactly in the log. Seeded.
 * Ref: Schwartz, E.S. (1997), The stochastic behavior of commodity prices, J. Finance 52(3), Model 1.
 * @param {{S0: number, kappa: number, mu: number, sigma: number, dt: number, steps: number, paths?: number,
 *   seed?: number}} o spot (currency), reversion speed (1/yr), long-run log-price level μ, volatility (1/√yr)
 * @returns {{t: number[], paths: Float64Array[], longRunMedian: number}} prices; longRunMedian = e^α
 */
export function schwartzSimulate({ S0, kappa, mu, sigma, dt, steps, paths = 1, seed = 1 }) {
  requireNumber('S0', S0, { gt: 0 }); requireNumber('kappa', kappa, { gt: 0 });
  const alpha = mu - sigma * sigma / (2 * kappa);
  const r = ouSimulate({ x0: Math.log(S0), kappa, theta: alpha, sigma, dt, steps, paths, seed });
  for (const p of r.paths) for (let k = 0; k < p.length; k++) p[k] = Math.exp(p[k]);
  return { t: r.t, paths: r.paths, longRunMedian: Math.exp(alpha) };
}

/**
 * Schwartz one-factor parameters from a price history (OU estimation on ln S, then μ = α + σ²/(2κ)).
 * Ref: Schwartz (1997), J. Finance 52(3).
 * @param {number[]} prices positive prices at equal spacing
 * @param {number} dt spacing (yr)
 * @returns {{kappa: number, mu: number, alpha: number, sigma: number, halfLife: number, meanReverting: boolean}}
 */
export function schwartzEstimate(prices, dt) {
  requireSeries('prices', prices, { minLength: 4 });
  const e = ouEstimate(Array.from(prices, (p) => { if (!(p > 0)) throw new RangeError('prices must be positive'); return Math.log(p); }), dt);
  return { kappa: e.kappa, alpha: e.theta, mu: e.theta + e.sigma * e.sigma / (2 * e.kappa), sigma: e.sigma, halfLife: e.halfLife,
    meanReverting: e.meanReverting };
}

/**
 * Merton jump-diffusion paths: GBM plus compound-Poisson lognormal jumps, with the drift compensated so that
 * E[S_t] = S0 e^(μt). Seeded.
 * Ref: Merton, R.C. (1976), Option pricing when underlying stock returns are discontinuous, J. Financial
 * Economics 3(1–2).
 * @param {{S0: number, mu: number, sigma: number, lambda: number, jumpMean: number, jumpStd: number, T: number,
 *   steps: number, paths?: number, seed?: number}} o lambda = jumps per year; jump size ln J ~ N(jumpMean, jumpStd²)
 * @returns {{t: number[], paths: Float64Array[]}}
 */
export function jumpDiffusionSimulate({ S0, mu, sigma, lambda, jumpMean, jumpStd, T, steps, paths = 1000, seed = 1 }) {
  requireNumber('S0', S0, { gt: 0 }); requireNumber('lambda', lambda, { min: 0 }); requireNumber('jumpStd', jumpStd, { min: 0 });
  requireNumber('T', T, { gt: 0 }); requireNumber('steps', steps, { min: 1, integer: true });
  const dt = T / steps; const k = Math.exp(jumpMean + 0.5 * jumpStd * jumpStd) - 1;
  const drift = (mu - lambda * k - 0.5 * sigma * sigma) * dt; const vol = sigma * Math.sqrt(dt);
  const rng = makeRng(seed); const out = []; const L = Math.exp(-lambda * dt);
  for (let p = 0; p < paths; p++) {
    const a = new Float64Array(steps + 1); a[0] = S0;
    for (let s = 1; s <= steps; s++) {
      let nJ = 0; let prod = rng.next();
      while (prod > L && nJ < 1000) { nJ++; prod *= rng.next(); }       // Knuth's Poisson sampler
      let jump = 0; for (let j = 0; j < nJ; j++) jump += jumpMean + jumpStd * rng.normal();
      a[s] = a[s - 1] * Math.exp(drift + vol * rng.normal() + jump);
    }
    out.push(a);
  }
  return { t: Array.from({ length: steps + 1 }, (_, s) => s * dt), paths: out };
}

/**
 * Scenario cone from simulated paths: cumulative percentiles (default P10/P50/P90) at every time step. The
 * convention is cumulative: P90 is the value not exceeded with 90 % probability.
 * Ref: fan-chart presentation of simulated paths.
 * @param {ArrayLike<number>[]} paths simulated paths of equal length
 * @param {number[]} [percentiles] cumulative percentiles in percent
 * @returns {{bands: Record<string, number[]>, mean: number[], convention: string}} bands keyed 'p10', 'p50', …
 */
export function scenarioCone(paths, percentiles = [10, 50, 90]) {
  if (!Array.isArray(paths) || paths.length < 2) throw new RangeError('paths must contain at least 2 paths');
  const m = paths[0].length; const bands = {}; const mean = new Array(m).fill(0);
  for (const p of percentiles) bands[`p${p}`] = new Array(m);
  const col = new Float64Array(paths.length);
  for (let k = 0; k < m; k++) {
    for (let i = 0; i < paths.length; i++) { col[i] = paths[i][k]; mean[k] += col[i] / paths.length; }
    const s = Float64Array.from(col).sort();
    for (const p of percentiles) bands[`p${p}`][k] = quantile(s, p / 100);
  }
  return { bands, mean, convention: 'Cumulative: pXX is the value not exceeded with XX % probability.' };
}

/**
 * Analytic lognormal cone around a forward curve: P_p(t) = F(t) exp(−σ²t/2 + z_p σ √t), so the forward is the
 * mean of the distribution at each date.
 * Ref: Black, F. (1976), The pricing of commodity contracts, J. Financial Economics 3(1–2) (lognormal forward).
 * @param {{forward: number[], times: number[], sigma: number, percentiles?: number[]}} o forward prices (currency),
 *   their times (yr), volatility (1/√yr), cumulative percentiles in percent
 * @returns {{bands: Record<string, number[]>, convention: string}}
 */
export function forwardCurveCone({ forward, times, sigma, percentiles = [10, 50, 90] }) {
  requireSeries('forward', forward); requireSeries('times', times); requireNumber('sigma', sigma, { min: 0 });
  if (forward.length !== times.length) throw new RangeError('forward and times must have the same length');
  const bands = {};
  for (const p of percentiles) {
    const z = normInv(requireNumber('percentile', p, { gt: 0, lt: 100 }) / 100);
    bands[`p${p}`] = Array.from(forward, (f, k) => f * Math.exp(-0.5 * sigma * sigma * times[k] + z * sigma * Math.sqrt(Math.max(0, times[k]))));
  }
  return { bands, convention: 'Cumulative: pXX is the value not exceeded with XX % probability.' };
}

// ---------------------------------------------------------------------------------------------------------------
// Time-series wrappers and back-testing
// ---------------------------------------------------------------------------------------------------------------

/**
 * Classical additive seasonal decomposition: centred moving-average trend, seasonal component as the mean
 * detrended value per position (normalised to sum to zero) and remainder.
 * Ref: Hyndman & Athanasopoulos, Forecasting: Principles and Practice, ch. 3 (classical decomposition).
 * @param {number[]} x series
 * @param {number} period season length in samples (>= 2; at least two full seasons of data)
 * @returns {{trend: (number | null)[], seasonal: number[], remainder: (number | null)[], seasonalIndices: number[]}}
 */
export function seasonalDecompose(x, period) {
  requireSeries('x', x); requireNumber('period', period, { min: 2, integer: true });
  const n = x.length;
  if (n < 2 * period) throw new RangeError('x must cover at least two full seasons');
  const trend = new Array(n).fill(null); const half = Math.floor(period / 2);
  for (let t = half; t < n - half; t++) {
    let s = 0;
    if (period % 2) for (let k = -half; k <= half; k++) s += x[t + k];
    else { for (let k = -half + 1; k < half; k++) s += x[t + k]; s += 0.5 * (x[t - half] + x[t + half]); }
    trend[t] = s / period;
  }
  const sum = new Array(period).fill(0); const cnt = new Array(period).fill(0);
  for (let t = 0; t < n; t++) if (trend[t] !== null) { sum[t % period] += x[t] - trend[t]; cnt[t % period]++; }
  const idx = sum.map((s, k) => (cnt[k] ? s / cnt[k] : 0));
  const mean = idx.reduce((a, v) => a + v, 0) / period;
  const seasonalIndices = idx.map((v) => v - mean);
  const seasonal = Array.from({ length: n }, (_, t) => seasonalIndices[t % period]);
  return { trend, seasonal, remainder: trend.map((tr, t) => (tr === null ? null : x[t] - tr - seasonal[t])), seasonalIndices };
}

/**
 * Forecast error metrics on paired predictions and observations: MAPE (mean |error|/|obs|), RMSE, MAE and bias
 * (mean of prediction − observation).
 * Ref: Hyndman & Koehler (2006), Another look at measures of forecast accuracy, Int. J. Forecasting 22(4).
 * @param {ArrayLike<number>} pred predictions
 * @param {ArrayLike<number>} obs observations (same length; MAPE skips zero observations)
 * @returns {{mape: number, rmse: number, mae: number, bias: number, n: number}} mape as a fraction
 */
export function errorMetrics(pred, obs) {
  if (pred.length !== obs.length || pred.length === 0) throw new RangeError('pred and obs must have equal non-zero length');
  let ape = 0; let na = 0; let se = 0; let ae = 0; let b = 0; const n = pred.length;
  for (let k = 0; k < n; k++) {
    const e = pred[k] - obs[k]; se += e * e; ae += Math.abs(e); b += e;
    if (obs[k] !== 0) { ape += Math.abs(e / obs[k]); na++; }
  }
  return { mape: na ? ape / na : NaN, rmse: Math.sqrt(se / n), mae: ae / n, bias: b / n, n };
}

/**
 * Holt / Holt–Winters forecast (wrapper over the shared exponential-smoothing routine).
 * Ref: Holt (1957); Winters (1960), Management Science 6(3).
 * @param {{x: number[], horizon: number, season?: number, trend?: boolean, seasonal?: 'additive' | 'multiplicative',
 *   alpha?: number, beta?: number, gamma?: number}} o omitted smoothing constants are estimated
 * @returns {{forecast: number[], fitted: number[], rmse: number, params: {alpha: number, beta: number | null,
 *   gamma: number | null}}}
 */
export function holtWintersForecast({ x, horizon, season, trend = true, seasonal, alpha, beta, gamma }) {
  requireSeries('x', x, { minLength: 3 }); requireNumber('horizon', horizon, { min: 1, integer: true });
  const opts = { trend };
  if (season) opts.season = season;
  if (seasonal) opts.seasonal = seasonal;
  if (alpha !== undefined) opts.alpha = alpha;
  if (beta !== undefined) opts.beta = beta;
  if (gamma !== undefined) opts.gamma = gamma;
  const m = expSmoothing(x, opts);
  return { forecast: Array.from(m.forecast(horizon)), fitted: Array.from(m.fitted), rmse: m.rmse,
    params: { alpha: m.alpha, beta: m.beta, gamma: m.gamma } };
}

/**
 * ARIMA(p, d, q) forecast with 95 % limits, optionally seasonal: with `season` and D = 1 the series is first
 * seasonally differenced (x_t − x_(t−s)), the ARIMA is fitted to the result and the forecast is integrated back
 * (limits are then reported for the mean path only as mean ± the differenced-series limits).
 * Ref: Box, Jenkins & Reinsel, Time Series Analysis, 4th ed., ch. 5, 7 and 9 (seasonal models).
 * @param {{x: number[], p?: number, d?: number, q?: number, horizon: number, season?: number, D?: number}} o
 * @returns {{mean: number[], lo95: number[], hi95: number[], aic: number, sigma2: number, stationary: boolean}}
 */
export function arimaForecast({ x, p = 1, d = 0, q = 0, horizon, season, D = 0 }) {
  requireSeries('x', x, { minLength: 6 }); requireNumber('horizon', horizon, { min: 1, integer: true });
  if (!season || !D) {
    const m = arima(x, { p, d, q }); const f = m.forecast(horizon);
    return { mean: Array.from(f.mean), lo95: Array.from(f.lo95), hi95: Array.from(f.hi95), aic: m.aic, sigma2: m.sigma2, stationary: m.stationary };
  }
  requireNumber('season', season, { min: 2, integer: true });
  if (D !== 1) throw new RangeError('only D = 0 or 1 is supported');
  if (x.length < 2 * season + 4) throw new RangeError('x is too short for seasonal differencing');
  const w = []; for (let t = season; t < x.length; t++) w.push(x[t] - x[t - season]);
  const m = arima(w, { p, d, q }); const f = m.forecast(horizon);
  const ext = Array.from(x); const mean = []; const lo95 = []; const hi95 = [];
  for (let h = 0; h < horizon; h++) {
    const v = f.mean[h] + ext[ext.length - season];
    ext.push(v); mean.push(v); lo95.push(v - (f.mean[h] - f.lo95[h])); hi95.push(v + (f.hi95[h] - f.mean[h]));
  }
  return { mean, lo95, hi95, aic: m.aic, sigma2: m.sigma2, stationary: m.stationary };
}

/**
 * Hold-out back-test of any forecasting method: the last `holdout` points are withheld, the method is fitted to
 * the rest and its forecast is scored against the withheld data; a naive last-value forecast is scored alongside
 * as the benchmark to beat.
 * Ref: Hyndman & Athanasopoulos, Forecasting: Principles and Practice, ch. 5 (evaluating forecast accuracy).
 * @param {{x: number[], holdout: number, forecast: (train: number[], h: number) => ArrayLike<number>}} o
 * @returns {{metrics: object, naive: object, skill: number, forecast: number[], actual: number[]}} skill =
 *   1 − RMSE / RMSE_naive (positive = better than naive)
 */
export function backtest({ x, holdout, forecast }) {
  requireSeries('x', x, { minLength: 4 }); requireNumber('holdout', holdout, { min: 1, max: x.length - 2, integer: true });
  const train = Array.from(x).slice(0, x.length - holdout); const actual = Array.from(x).slice(x.length - holdout);
  const f = Array.from(forecast(train, holdout));
  if (f.length !== holdout) throw new RangeError('forecast must return `holdout` values');
  const metrics = errorMetrics(f, actual);
  const naive = errorMetrics(new Array(holdout).fill(train[train.length - 1]), actual);
  return { metrics, naive, skill: naive.rmse > 0 ? 1 - metrics.rmse / naive.rmse : NaN, forecast: f, actual };
}

/**
 * Local-level (random walk plus noise) state-space model filtered with the shared Kalman filter; the forecast is
 * the last filtered level with variance growing by q per step.
 * Ref: Harvey, A.C. (1989), Forecasting, Structural Time Series Models and the Kalman Filter, ch. 2; Durbin &
 * Koopman, Time Series Analysis by State Space Methods, ch. 2.
 * @param {{y: (number | null)[], q: number, r: number, horizon?: number}} o observations (null = missing), level
 *   noise variance q and observation noise variance r (unit²)
 * @returns {{level: number[], variance: number[], forecast: number[], forecastStd: number[], logLik: number}}
 */
export function kalmanLocalLevel({ y, q, r, horizon = 0 }) {
  if (!Array.isArray(y) || y.length < 2) throw new RangeError('y must have at least 2 observations');
  requireNumber('q', q, { min: 0 }); requireNumber('r', r, { gt: 0 });
  const first = y.find((v) => v !== null && Number.isFinite(v));
  const res = kalmanFilter({ F: [[1]], H: [[1]], Q: [[q]], R: [[r]], x0: [first], P0: [[1e6 * (r + q)]] }, y);
  const level = res.x.map((v) => v[0]); const variance = res.P.map((P) => P[0][0]);
  const last = level[level.length - 1]; const pv = variance[variance.length - 1];
  return { level, variance, forecast: new Array(horizon).fill(last),
    forecastStd: Array.from({ length: horizon }, (_, h) => Math.sqrt(pv + (h + 1) * q + r)), logLik: res.logLik };
}

/**
 * Vector autoregression of order 1, y_t = c + A y_(t−1) + ε, by equation-wise least squares, with iterated
 * multi-step forecasts.
 * Ref: Lütkepohl, H. (2005), New Introduction to Multiple Time Series Analysis, ch. 3.
 * @param {number[][]} series rows = time, columns = variables (T × k, T > k + 2)
 * @returns {{c: number[], A: number[][], residualStd: number[], stable: boolean | null,
 *   forecast: (h: number) => number[][]}} forecast(h) returns h rows
 */
export function var1Fit(series) {
  if (!Array.isArray(series) || series.length < 4) throw new RangeError('series needs at least 4 rows');
  const k = series[0].length; const T = series.length - 1;
  if (T <= k + 1) throw new RangeError('series is too short for the number of variables');
  const G = Array.from({ length: k + 1 }, () => new Array(k + 1).fill(0));
  const rhs = Array.from({ length: k }, () => new Array(k + 1).fill(0));
  for (let t = 0; t < T; t++) {
    const z = [1, ...series[t]];
    for (let a = 0; a <= k; a++) { for (let b = 0; b <= k; b++) G[a][b] += z[a] * z[b]; for (let e = 0; e < k; e++) rhs[e][a] += z[a] * series[t + 1][e]; }
  }
  const c = []; const A = [];
  for (let e = 0; e < k; e++) { const beta = solveDense(G.map((r) => r.slice()), rhs[e]); c.push(beta[0]); A.push(Array.from(beta).slice(1)); }
  const residualStd = new Array(k).fill(0);
  for (let t = 0; t < T; t++) for (let e = 0; e < k; e++) {
    let p = c[e]; for (let j = 0; j < k; j++) p += A[e][j] * series[t][j];
    residualStd[e] += (series[t + 1][e] - p) ** 2;
  }
  for (let e = 0; e < k; e++) residualStd[e] = Math.sqrt(residualStd[e] / (T - k - 1));
  // stability by the induced infinity norm (sufficient, not necessary)
  const normInf = Math.max(...A.map((r) => r.reduce((s, v) => s + Math.abs(v), 0)));
  const last = series[series.length - 1];
  return { c, A, residualStd, stable: normInf < 1 ? true : null,
    forecast: (h) => { const out = []; let y = Array.from(last); for (let s = 0; s < h; s++) { y = c.map((ci, e) => ci + A[e].reduce((a, v, j) => a + v * y[j], 0)); out.push(y); } return out; } };
}

// ---------------------------------------------------------------------------------------------------------------
// Regression cost models
// ---------------------------------------------------------------------------------------------------------------

function tCdf(t, dof) {
  const x = dof / (dof + t * t);
  const tail = 0.5 * betaI(dof / 2, 0.5, x);
  return t >= 0 ? 1 - tail : tail;
}

/**
 * Quantile of Student's t distribution by bisection on the exact cdf (regularised incomplete beta function).
 * Ref: Abramowitz & Stegun 26.7.1.
 * @param {number} p probability in (0, 1)
 * @param {number} dof degrees of freedom (> 0)
 * @returns {number}
 */
export function tQuantile(p, dof) {
  requireNumber('p', p, { gt: 0, lt: 1 }); requireNumber('dof', dof, { gt: 0 });
  let lo = -1e3; let hi = 1e3;
  for (let it = 0; it < 200; it++) { const mid = 0.5 * (lo + hi); if (tCdf(mid, dof) < p) lo = mid; else hi = mid; }
  return 0.5 * (lo + hi);
}

/**
 * Regression-based cost model (cost-estimating relationship): ordinary least squares of cost on one or more
 * drivers, with coefficient standard errors and, for a new point, the confidence interval of the mean and the
 * prediction interval of a single outcome. With logLog the model is fitted as ln C = β0 + Σ βj ln xj (power law)
 * and results are transformed back (intervals become multiplicative; the back-transformed mean is the median).
 * Ref: Draper & Smith, Applied Regression Analysis, 3rd ed., ch. 3 and 5; ISPA/SCEA Parametric Estimating
 * Handbook, 4th ed., ch. 3.
 * @param {{X: number[][], y: number[], logLog?: boolean}} o rows of driver values and observed costs (n > p + 1)
 * @returns {{beta: number[], se: number[], r2: number, adjR2: number, sigma: number, dof: number,
 *   predict: (x: number[], level?: number) => {mean: number, ciLow: number, ciHigh: number, piLow: number,
 *   piHigh: number}}} beta[0] is the intercept; level = two-sided confidence (default 0.9)
 */
export function regressionCostModel({ X, y, logLog = false }) {
  if (!Array.isArray(X) || X.length !== y.length) throw new RangeError('X must have one row per y');
  const n = y.length; const p = X[0].length + 1;
  if (n <= p) throw new RangeError('need more observations than coefficients');
  const tx = (v) => { if (logLog && !(v > 0)) throw new RangeError('logLog needs positive values'); return logLog ? Math.log(v) : v; };
  const A = X.map((r) => [1, ...Array.from(r, tx)]); const yy = Array.from(y, tx);
  const G = Array.from({ length: p }, (_, a) => Array.from({ length: p }, (_, b) => A.reduce((s, r) => s + r[a] * r[b], 0)));
  const Gi = inverse(G);
  const Aty = Array.from({ length: p }, (_, a) => A.reduce((s, r, i) => s + r[a] * yy[i], 0));
  const beta = Gi.map((row) => row.reduce((s, v, j) => s + v * Aty[j], 0));
  const mean = yy.reduce((s, v) => s + v, 0) / n; let sse = 0; let sst = 0;
  for (let i = 0; i < n; i++) { const f = A[i].reduce((s, v, j) => s + v * beta[j], 0); sse += (yy[i] - f) ** 2; sst += (yy[i] - mean) ** 2; }
  const dof = n - p; const s2 = sse / dof;
  const predict = (x, level = 0.9) => {
    if (!x || x.length !== p - 1) throw new RangeError(`x must have ${p - 1} value(s)`);
    const a = [1, ...Array.from(x, tx)];
    const m = a.reduce((s, v, j) => s + v * beta[j], 0);
    let h = 0; for (let i = 0; i < p; i++) for (let j = 0; j < p; j++) h += a[i] * Gi[i][j] * a[j];
    const t = tQuantile(0.5 + level / 2, dof);
    const ci = t * Math.sqrt(s2 * h); const pi = t * Math.sqrt(s2 * (1 + h));
    const back = (v) => (logLog ? Math.exp(v) : v);
    return { mean: back(m), ciLow: back(m - ci), ciHigh: back(m + ci), piLow: back(m - pi), piHigh: back(m + pi) };
  };
  return { beta, se: Gi.map((r, j) => Math.sqrt(s2 * r[j])), r2: sst > 0 ? 1 - sse / sst : NaN,
    adjR2: sst > 0 ? 1 - (sse / dof) / (sst / (n - 1)) : NaN, sigma: Math.sqrt(s2), dof, predict };
}

/**
 * Generalised linear model by iteratively reweighted least squares: 'gaussian' (identity link), 'poisson' (log
 * link, e.g. failure counts) or 'binomial' (logit link, e.g. probability of overrun).
 * Ref: McCullagh & Nelder (1989), Generalized Linear Models, 2nd ed., §2.5 (IRLS).
 * @param {{X: number[][], y: number[], family?: 'gaussian' | 'poisson' | 'binomial', maxIter?: number}} o rows of
 *   regressors (intercept added) and responses (counts for poisson, 0/1 or proportions for binomial)
 * @returns {{beta: number[], deviance: number, iterations: number, converged: boolean, predict: (x: number[]) => number}}
 *   predict returns the mean response
 */
export function glmFit({ X, y, family = 'gaussian', maxIter = 50 }) {
  if (!Array.isArray(X) || X.length !== y.length || X.length < 2) throw new RangeError('X must have one row per y');
  const n = y.length; const p = X[0].length + 1; const A = X.map((r) => [1, ...r]);
  const link = { gaussian: { inv: (e) => e, dmu: () => 1, v: () => 1 },
    poisson: { inv: (e) => Math.exp(Math.min(e, 700)), dmu: (m) => m, v: (m) => m },
    binomial: { inv: (e) => 1 / (1 + Math.exp(-e)), dmu: (m) => m * (1 - m), v: (m) => m * (1 - m) } }[family];
  if (!link) throw new RangeError("family must be 'gaussian', 'poisson' or 'binomial'");
  let beta = new Array(p).fill(0);
  if (family === 'poisson') beta[0] = Math.log(Math.max(1e-9, y.reduce((s, v) => s + v, 0) / n));
  let converged = false; let it = 0;
  for (; it < maxIter && !converged; it++) {
    const G = Array.from({ length: p }, () => new Array(p).fill(0)); const b = new Array(p).fill(0);
    for (let i = 0; i < n; i++) {
      const eta = A[i].reduce((s, v, j) => s + v * beta[j], 0); const mu = link.inv(eta);
      const d = Math.max(1e-12, link.dmu(mu)); const w = d * d / Math.max(1e-12, link.v(mu));
      const z = eta + (y[i] - mu) / d;
      for (let a = 0; a < p; a++) { b[a] += w * A[i][a] * z; for (let c = 0; c < p; c++) G[a][c] += w * A[i][a] * A[i][c]; }
    }
    const next = Array.from(solveDense(G, b));
    let diff = 0; for (let j = 0; j < p; j++) diff = Math.max(diff, Math.abs(next[j] - beta[j]));
    beta = next; converged = diff < 1e-10 * (1 + Math.max(...beta.map(Math.abs)));
  }
  const predict = (x) => link.inv([1, ...x].reduce((s, v, j) => s + v * beta[j], 0));
  let dev = 0;
  for (let i = 0; i < n; i++) {
    const mu = predict(X[i]);
    if (family === 'gaussian') dev += (y[i] - mu) ** 2;
    else if (family === 'poisson') dev += 2 * ((y[i] > 0 ? y[i] * Math.log(y[i] / mu) : 0) - (y[i] - mu));
    else dev += -2 * (y[i] * Math.log(Math.max(mu, 1e-300)) + (1 - y[i]) * Math.log(Math.max(1 - mu, 1e-300)));
  }
  return { beta, deviance: dev, iterations: it, converged, predict };
}

/**
 * Gaussian-process regression of cost (or any response) on drivers, returning mean and standard deviation of the
 * prediction (wrapper over the shared surrogate).
 * Ref: Rasmussen & Williams (2006), Gaussian Processes for Machine Learning, ch. 2.
 * @param {{X: number[][], y: number[], kernel?: 'rbf' | 'matern52', noise?: number}} o noise = observation variance
 * @returns {{predict: (x: number[]) => {mean: number, std: number}, logLik: number, params: object}}
 */
export function gpRegression({ X, y, kernel = 'matern52', noise }) {
  const gp = gaussianProcess(X, y, noise !== undefined ? { kernel, noise } : { kernel });
  return { predict: (x) => gp.predict(x), logLik: gp.logLik, params: gp.params };
}

// --- regression trees (shared by boosting and the forest) ---
function buildTree(X, y, idx, depth, maxDepth, minLeaf, pickFeatures) {
  let mean = 0; for (const i of idx) mean += y[i]; mean /= idx.length;
  if (depth >= maxDepth || idx.length < 2 * minLeaf) return { value: mean };
  let best = null; let sst = 0; for (const i of idx) sst += (y[i] - mean) ** 2;
  for (const f of pickFeatures(X[0].length)) {
    const order = idx.slice().sort((a, b) => X[a][f] - X[b][f]);
    let sl = 0; let sl2 = 0; let st = 0; let st2 = 0;
    for (const i of order) { st += y[i]; st2 += y[i] * y[i]; }
    for (let k = 0; k < order.length - 1; k++) {
      const v = y[order[k]]; sl += v; sl2 += v * v;
      if (k + 1 < minLeaf || order.length - k - 1 < minLeaf || X[order[k]][f] === X[order[k + 1]][f]) continue;
      const nl = k + 1; const nr = order.length - nl;
      const sse = (sl2 - sl * sl / nl) + (st2 - sl2 - (st - sl) ** 2 / nr);
      if (!best || sse < best.sse) best = { sse, f, thr: 0.5 * (X[order[k]][f] + X[order[k + 1]][f]), left: order.slice(0, nl), right: order.slice(nl) };
    }
  }
  if (!best || best.sse >= sst - 1e-12 * Math.max(1, sst)) return { value: mean };
  return { f: best.f, thr: best.thr, left: buildTree(X, y, best.left, depth + 1, maxDepth, minLeaf, pickFeatures),
    right: buildTree(X, y, best.right, depth + 1, maxDepth, minLeaf, pickFeatures) };
}
const treePredict = (node, x) => { let nd = node; while (nd.value === undefined) nd = x[nd.f] <= nd.thr ? nd.left : nd.right; return nd.value; };
const checkXY = (X, y) => {
  if (!Array.isArray(X) || X.length !== y.length || X.length < 4) throw new RangeError('X must have one row per y (at least 4 rows)');
};

/**
 * Gradient-boosted regression trees with squared-error loss: each shallow tree is fitted to the current residuals
 * and added with a shrinkage factor.
 * Ref: Friedman, J.H. (2001), Greedy function approximation: a gradient boosting machine, Annals of Statistics 29(5).
 * @param {{X: number[][], y: number[], trees?: number, learningRate?: number, maxDepth?: number, minLeaf?: number}} o
 * @returns {{predict: (x: number[]) => number, trainRmse: number, history: number[]}} history = training RMSE
 *   after each tree
 */
export function gradientBoosting({ X, y, trees = 100, learningRate = 0.1, maxDepth = 2, minLeaf = 2 }) {
  checkXY(X, y); requireNumber('learningRate', learningRate, { gt: 0, max: 1 });
  const n = y.length; const f0 = y.reduce((s, v) => s + v, 0) / n;
  const pred = new Float64Array(n).fill(f0); const res = new Float64Array(n); const forest = []; const history = [];
  const all = Array.from({ length: n }, (_, i) => i); const feats = (d) => Array.from({ length: d }, (_, j) => j);
  for (let m = 0; m < trees; m++) {
    for (let i = 0; i < n; i++) res[i] = y[i] - pred[i];
    const tree = buildTree(X, res, all, 0, maxDepth, minLeaf, feats);
    forest.push(tree);
    let se = 0; for (let i = 0; i < n; i++) { pred[i] += learningRate * treePredict(tree, X[i]); se += (y[i] - pred[i]) ** 2; }
    history.push(Math.sqrt(se / n));
  }
  return { predict: (x) => forest.reduce((s, t) => s + learningRate * treePredict(t, x), f0), trainRmse: history[history.length - 1], history };
}

/**
 * Random-forest regression: bootstrap-aggregated regression trees with a random feature subset at each split;
 * the out-of-bag RMSE is an honest error estimate. Seeded.
 * Ref: Breiman, L. (2001), Random forests, Machine Learning 45(1).
 * @param {{X: number[][], y: number[], trees?: number, maxDepth?: number, minLeaf?: number, features?: number,
 *   seed?: number}} o features = number tried per split (default ⌈d/3⌉)
 * @returns {{predict: (x: number[]) => number, oobRmse: number}}
 */
export function randomForest({ X, y, trees = 100, maxDepth = 6, minLeaf = 2, features, seed = 1 }) {
  checkXY(X, y);
  const n = y.length; const d = X[0].length; const mtry = features ?? Math.max(1, Math.ceil(d / 3));
  const rng = makeRng(seed); const forest = [];
  const oobSum = new Float64Array(n); const oobCnt = new Float64Array(n);
  const pick = (dd) => { const a = Array.from({ length: dd }, (_, j) => j); for (let i = dd - 1; i > 0; i--) { const r = Math.floor(rng.next() * (i + 1)); [a[i], a[r]] = [a[r], a[i]]; } return a.slice(0, mtry); };
  for (let m = 0; m < trees; m++) {
    const inBag = new Uint8Array(n); const idx = [];
    for (let i = 0; i < n; i++) { const r = Math.floor(rng.next() * n); idx.push(r); inBag[r] = 1; }
    const tree = buildTree(X, y, idx, 0, maxDepth, minLeaf, pick);
    forest.push(tree);
    for (let i = 0; i < n; i++) if (!inBag[i]) { oobSum[i] += treePredict(tree, X[i]); oobCnt[i]++; }
  }
  let se = 0; let c = 0;
  for (let i = 0; i < n; i++) if (oobCnt[i] > 0) { se += (y[i] - oobSum[i] / oobCnt[i]) ** 2; c++; }
  return { predict: (x) => forest.reduce((s, t) => s + treePredict(t, x), 0) / forest.length, oobRmse: c ? Math.sqrt(se / c) : NaN };
}

/**
 * Small feed-forward neural network for regression: one tanh hidden layer, linear output, inputs and output
 * standardised, trained by full-batch gradient descent with momentum. Seeded. Intended for smooth low-dimensional
 * relationships with enough trustworthy data; always compare with the regression model on a hold-out set.
 * Ref: Rumelhart, Hinton & Williams (1986), Learning representations by back-propagating errors, Nature 323.
 * @param {{X: number[][], y: number[], hidden?: number, epochs?: number, learningRate?: number, momentum?: number,
 *   seed?: number}} o
 * @returns {{predict: (x: number[]) => number, trainRmse: number, history: number[]}} history = training RMSE (in
 *   the units of y) sampled every 10 epochs
 */
export function mlpRegression({ X, y, hidden = 6, epochs = 2000, learningRate = 0.05, momentum = 0.9, seed = 1 }) {
  checkXY(X, y);
  const n = y.length; const d = X[0].length; const rng = makeRng(seed);
  const mx = new Array(d).fill(0); const sx = new Array(d).fill(0);
  for (let j = 0; j < d; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j] / n; for (let i = 0; i < n; i++) sx[j] += (X[i][j] - mx[j]) ** 2; sx[j] = Math.sqrt(sx[j] / n) || 1; }
  const my = y.reduce((s, v) => s + v, 0) / n; const sy = Math.sqrt(y.reduce((s, v) => s + (v - my) ** 2, 0) / n) || 1;
  const Z = X.map((r) => r.map((v, j) => (v - mx[j]) / sx[j])); const t = y.map((v) => (v - my) / sy);
  const W1 = Array.from({ length: hidden }, () => Array.from({ length: d + 1 }, () => (rng.next() - 0.5) / Math.sqrt(d)));
  const W2 = Array.from({ length: hidden + 1 }, () => (rng.next() - 0.5) / Math.sqrt(hidden));
  const V1 = W1.map((r) => r.map(() => 0)); const V2 = W2.map(() => 0);
  const hOut = new Float64Array(hidden); const history = [];
  const fwd = (z) => { let o = W2[hidden]; for (let h = 0; h < hidden; h++) { let a = W1[h][d]; for (let j = 0; j < d; j++) a += W1[h][j] * z[j]; hOut[h] = Math.tanh(a); o += W2[h] * hOut[h]; } return o; };
  for (let e = 0; e < epochs; e++) {
    const g1 = W1.map((r) => r.map(() => 0)); const g2 = W2.map(() => 0); let se = 0;
    for (let i = 0; i < n; i++) {
      const err = fwd(Z[i]) - t[i]; se += err * err;
      for (let h = 0; h < hidden; h++) {
        g2[h] += err * hOut[h];
        const dh = err * W2[h] * (1 - hOut[h] * hOut[h]);
        for (let j = 0; j < d; j++) g1[h][j] += dh * Z[i][j];
        g1[h][d] += dh;
      }
      g2[hidden] += err;
    }
    for (let h = 0; h <= hidden; h++) { V2[h] = momentum * V2[h] - learningRate * g2[h] / n; W2[h] += V2[h]; }
    for (let h = 0; h < hidden; h++) for (let j = 0; j <= d; j++) { V1[h][j] = momentum * V1[h][j] - learningRate * g1[h][j] / n; W1[h][j] += V1[h][j]; }
    if (e % 10 === 0 || e === epochs - 1) history.push(Math.sqrt(se / n) * sy);
  }
  const predict = (x) => fwd(x.map((v, j) => (v - mx[j]) / sx[j])) * sy + my;
  let se = 0; for (let i = 0; i < n; i++) se += (predict(X[i]) - y[i]) ** 2;
  return { predict, trainRmse: Math.sqrt(se / n), history };
}

// ---------------------------------------------------------------------------------------------------------------
// Bayesian updating, reference classes, escalation and scenarios
// ---------------------------------------------------------------------------------------------------------------

/**
 * Bayesian updating of a cost estimate with the conjugate normal model: prior N(m0, s0²), n independent quotes
 * with known sampling standard deviation s; posterior precision = 1/s0² + n/s², posterior mean is the
 * precision-weighted average.
 * Ref: Gelman et al., Bayesian Data Analysis, 3rd ed., §2.5 (normal model with known variance).
 * @param {{priorMean: number, priorStd: number, observations: number[], obsStd: number}} o currency
 * @returns {{mean: number, std: number, priorWeight: number, predictiveStd: number, p10: number, p50: number,
 *   p90: number, convention: string}} posterior of the true cost; predictiveStd applies to the next quote
 */
export function bayesNormalUpdate({ priorMean, priorStd, observations, obsStd }) {
  requireNumber('priorMean', priorMean); requireNumber('priorStd', priorStd, { gt: 0 }); requireNumber('obsStd', obsStd, { gt: 0 });
  const n = observations.length; const ybar = n ? requireSeries('observations', observations).reduce((s, v) => s + v, 0) / n : 0;
  const pp = 1 / (priorStd * priorStd); const pl = n / (obsStd * obsStd);
  const mean = (pp * priorMean + pl * ybar) / (pp + pl); const std = Math.sqrt(1 / (pp + pl));
  return { mean, std, priorWeight: pp / (pp + pl), predictiveStd: Math.sqrt(std * std + obsStd * obsStd),
    p10: mean + normInv(0.1) * std, p50: mean, p90: mean + normInv(0.9) * std,
    convention: 'Cumulative: P90 is the cost not exceeded with 90 % probability.' };
}

/**
 * Normal prior for a cost estimate from its AACE class: the class's wide accuracy band is read as an 80 %
 * interval (P10–P90), giving std = (high − low)/(2 × 1.2816) and a mean shifted to the middle of the band
 * (estimates are skewed to overrun).
 * Ref: AACE International RP 18R-97 (accuracy ranges stated at about 80 % confidence); the normal approximation
 * of the band is a modelling choice.
 * @param {{estimate: number, cls: 1 | 2 | 3 | 4 | 5, band?: 'wide' | 'narrow'}} o currency
 * @returns {{priorMean: number, priorStd: number, low: number, high: number}}
 */
export function priorFromEstimateClass({ estimate, cls, band = 'wide' }) {
  requireNumber('estimate', estimate, { gt: 0 });
  const c = AACE_CLASSES[cls];
  if (!c) throw new RangeError('cls must be 1, 2, 3, 4 or 5');
  const k = band === 'narrow' ? 0 : 1;
  const low = estimate * (1 + c.low[k]); const high = estimate * (1 + c.high[k]);
  return { priorMean: 0.5 * (low + high), priorStd: (high - low) / (2 * normInv(0.9)), low, high };
}

/**
 * Hierarchical (partial-pooling) estimate of group means with the normal–normal model fitted by empirical Bayes:
 * between-group variance τ² by the method of moments, each group mean shrunk towards the grand mean with weight
 * τ²/(τ² + s_j²/n_j). Useful for pooling sparse cost or failure data across similar projects or assets.
 * Ref: Gelman et al., Bayesian Data Analysis, 3rd ed., ch. 5 (hierarchical normal model); DerSimonian & Laird
 * (1986), Controlled Clinical Trials 7(3) for the moment estimator.
 * @param {{name: string, mean: number, variance: number, n: number}[]} groups sample mean, sample variance and size
 * @returns {{grandMean: number, tau2: number, groups: {name: string, raw: number, pooled: number, shrinkage: number,
 *   std: number}[]}} shrinkage = weight on the grand mean (0 = no pooling, 1 = complete pooling)
 */
export function hierarchicalPooling(groups) {
  if (!Array.isArray(groups) || groups.length < 2) throw new RangeError('at least two groups are required');
  const se2 = groups.map((g) => requireNumber('variance', g.variance, { min: 0 }) / requireNumber('n', g.n, { min: 1 }) || 1e-300);
  const w = se2.map((v) => 1 / v); const sw = w.reduce((s, v) => s + v, 0);
  const fixed = groups.reduce((s, g, j) => s + w[j] * g.mean, 0) / sw;
  const Q = groups.reduce((s, g, j) => s + w[j] * (g.mean - fixed) ** 2, 0);
  const c = sw - w.reduce((s, v) => s + v * v, 0) / sw;
  const tau2 = Math.max(0, (Q - (groups.length - 1)) / c);
  const wr = se2.map((v) => 1 / (v + tau2)); const swr = wr.reduce((s, v) => s + v, 0);
  const grandMean = groups.reduce((s, g, j) => s + wr[j] * g.mean, 0) / swr;
  return { grandMean, tau2, groups: groups.map((g, j) => {
    const B = tau2 > 0 ? se2[j] / (se2[j] + tau2) : 1;
    return { name: g.name, raw: g.mean, pooled: B * grandMean + (1 - B) * g.mean, shrinkage: B, std: Math.sqrt((1 - B) * se2[j] + B * B / swr) };
  }) };
}

/**
 * Editable reference-class uplift table: required uplift on the base estimate (fraction) against the acceptable
 * chance of cost overrun. The numbers are indicative placeholders with a plausible shape; they must be replaced
 * by the distribution of overruns of a genuinely comparable reference class.
 * Ref: method of Flyvbjerg, B. (2006), From Nobel Prize to project management: getting risks right, Project
 * Management Journal 37(3); values are placeholders, not taken from any published reference class.
 */
export const defaultOptimismBias = {
  source: 'indicative placeholder — replace with the overrun distribution of a comparable reference class',
  basisYear: null,
  acceptableOverrunProbability: [0.5, 0.4, 0.3, 0.2, 0.1],
  uplift: [0.15, 0.22, 0.30, 0.42, 0.65],
};

/**
 * Reference-class forecasting uplift: the uplift for an acceptable probability of overrun, interpolated in the
 * table, and the uplifted estimate.
 * Ref: Flyvbjerg (2006), Project Management Journal 37(3); Kahneman & Tversky (1979) (outside view).
 * @param {{estimate: number, acceptableOverrunProbability: number, table?: object}} o currency, probability within
 *   the table's range
 * @returns {{uplift: number, upliftedEstimate: number, usesPlaceholderTable: boolean}}
 */
export function referenceClassUplift({ estimate, acceptableOverrunProbability, table = defaultOptimismBias }) {
  requireNumber('estimate', estimate, { min: 0 });
  const ps = table.acceptableOverrunProbability; const us = table.uplift;
  const p = requireNumber('acceptableOverrunProbability', acceptableOverrunProbability, { min: Math.min(...ps), max: Math.max(...ps) });
  const pts = ps.map((v, k) => [v, us[k]]).sort((a, b) => a[0] - b[0]);
  let u = pts[pts.length - 1][1];
  for (let k = 1; k < pts.length; k++) if (p <= pts[k][0]) { u = pts[k - 1][1] + (pts[k][1] - pts[k - 1][1]) * (p - pts[k - 1][0]) / (pts[k][0] - pts[k - 1][0]); break; }
  return { uplift: u, upliftedEstimate: estimate * (1 + u), usesPlaceholderTable: table === defaultOptimismBias };
}

/**
 * Learning-curve cost projection along a path of cumulative volume.
 * Ref: Wright (1936); Boston Consulting Group (1968).
 * @param {{c0: number, q0: number, volumes: number[], learningRate: number}} o reference cost at cumulative volume
 *   q0 and the future cumulative volumes
 * @returns {number[]} cost at each cumulative volume (currency per unit)
 */
export function learningCurveProjection({ c0, q0, volumes, learningRate }) {
  return Array.from(requireSeries('volumes', volumes), (q) => costAtCumulativeVolume({ c0, q0, q, learningRate }));
}

/**
 * Escalation-index forecast from a history: compound annual growth ('cagr', between first and last point) or
 * least-squares linear trend, extrapolated over the horizon.
 * Ref: AACE International RP 58R-10 (escalation estimating using indices and forecasts).
 * @param {{years: number[], values: number[], horizon: number, method?: 'cagr' | 'linear'}} o index history at
 *   annual spacing and number of future years
 * @returns {{years: number[], values: number[], rate: number}} forecast years and index values; rate = implied
 *   annual growth (cagr) or slope/last value (linear)
 */
export function escalationForecast({ years, values, horizon, method = 'cagr' }) {
  requireSeries('years', years, { minLength: 2 }); requireSeries('values', values, { minLength: 2 });
  requireNumber('horizon', horizon, { min: 1, integer: true });
  const n = years.length; const y1 = years[n - 1]; const v1 = values[n - 1];
  const outY = Array.from({ length: horizon }, (_, k) => y1 + k + 1);
  if (method === 'cagr') {
    if (!(values[0] > 0 && v1 > 0)) throw new RangeError('index values must be positive');
    const g = Math.pow(v1 / values[0], 1 / (y1 - years[0])) - 1;
    return { years: outY, values: outY.map((y) => v1 * Math.pow(1 + g, y - y1)), rate: g };
  }
  if (method !== 'linear') throw new RangeError("method must be 'cagr' or 'linear'");
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0;
  for (let k = 0; k < n; k++) { sx += years[k]; sy += values[k]; sxx += years[k] * years[k]; sxy += years[k] * values[k]; }
  const b = (n * sxy - sx * sy) / (n * sxx - sx * sx); const a = (sy - b * sx) / n;
  return { years: outY, values: outY.map((y) => a + b * y), rate: b / v1 };
}

/**
 * Exchange-rate scenarios: the forward path from covered interest parity F_t = S [(1 + r_d)/(1 + r_f)]^t and a
 * lognormal uncertainty cone around it.
 * Ref: covered interest-rate parity, Hull, Options, Futures and Other Derivatives, ch. 5.
 * @param {{spot: number, domesticRate: number, foreignRate: number, sigma: number, years: number,
 *   percentiles?: number[]}} o spot = domestic currency per unit of foreign currency; annual rates (fractions);
 *   volatility (1/√yr)
 * @returns {{t: number[], forward: number[], bands: Record<string, number[]>, convention: string}}
 */
export function fxScenarios({ spot, domesticRate, foreignRate, sigma, years, percentiles = [10, 50, 90] }) {
  requireNumber('spot', spot, { gt: 0 }); requireNumber('domesticRate', domesticRate, { gt: -1 });
  requireNumber('foreignRate', foreignRate, { gt: -1 }); requireNumber('years', years, { min: 1, integer: true });
  const t = Array.from({ length: years + 1 }, (_, k) => k);
  const forward = t.map((k) => spot * Math.pow((1 + domesticRate) / (1 + foreignRate), k));
  return { t, forward, ...forwardCurveCone({ forward, times: t, sigma, percentiles }) };
}

/**
 * Carbon-price scenario paths: each scenario grows from its start price at its own annual rate, bounded by an
 * optional floor and cap. Start prices and growth rates are scenario assumptions supplied by the caller (or by
 * live data).
 * Ref: scenario construction for internal carbon pricing; no price data are embedded.
 * @param {{years: number, scenarios: Record<string, {start: number, growth: number}>, floor?: number, cap?: number}} o
 *   prices in currency per t CO2e, growth as a fraction per year
 * @returns {Record<string, number[]>} price per year (index 0 = first year) for each scenario
 */
export function carbonPricePaths({ years, scenarios, floor = 0, cap = Infinity }) {
  requireNumber('years', years, { min: 1, integer: true });
  const out = {};
  for (const name of Object.keys(scenarios || {})) {
    const s = scenarios[name];
    requireNumber(`${name}.start`, s.start, { min: 0 }); requireNumber(`${name}.growth`, s.growth, { gt: -1 });
    out[name] = Array.from({ length: years }, (_, k) => Math.min(cap, Math.max(floor, s.start * Math.pow(1 + s.growth, k))));
  }
  if (Object.keys(out).length === 0) throw new RangeError('scenarios must name at least one path');
  return out;
}
