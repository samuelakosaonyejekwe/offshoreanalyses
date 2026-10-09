/**
 * Time-series analysis: autocorrelation, differencing, smoothing, Holt–Winters exponential smoothing,
 * ARIMA(p, d, q) by conditional sum of squares, and the linear Kalman filter.
 *
 * Series are equally spaced samples (plain or typed arrays); results carry the units of the series [u].
 */
import { levenbergMarquardt, nelderMead } from './optim.js';
import { inverse } from './linalg.js';

function needSeries(x, name, nMin = 1) {
  if (!x || typeof x.length !== 'number' || x.length < nMin) throw new RangeError(`${name}: x needs at least ${nMin} value(s)`);
}

function meanOf(x) {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i];
  return s / x.length;
}

/**
 * Sample autocorrelation function r_k = c_k / c_0 with c_k = (1/n) Σ (x_t - x̄)(x_{t+k} - x̄).
 * Reference: Box, Jenkins & Reinsel, "Time Series Analysis", 4th ed., §2.1.5. Validity: 0 <= maxLag < n,
 * non-constant series.
 * @param {ArrayLike<number>} x series [u]
 * @param {number} maxLag largest lag [samples]
 * @returns {number[]} r_0..r_maxLag (r_0 = 1) [-]
 */
export function acf(x, maxLag) {
  needSeries(x, 'acf', 2);
  const n = x.length;
  if (!Number.isInteger(maxLag) || maxLag < 0 || maxLag >= n) throw new RangeError('acf: maxLag must be an integer within 0..n-1');
  const m = meanOf(x);
  let c0 = 0;
  for (let i = 0; i < n; i++) c0 += (x[i] - m) * (x[i] - m);
  const r = new Array(maxLag + 1);
  r[0] = 1;
  for (let k = 1; k <= maxLag; k++) {
    let s = 0;
    for (let i = 0; i + k < n; i++) s += (x[i] - m) * (x[i + k] - m);
    r[k] = c0 > 0 ? s / c0 : 0;
  }
  return r;
}

/** Durbin–Levinson recursion on autocorrelations r[0..K]; returns { pacf, phi (AR coefficients of order K) }. */
function durbinLevinson(r, K) {
  const pac = new Array(K + 1);
  pac[0] = 1;
  let phi = [];
  let v = 1;
  for (let k = 1; k <= K; k++) {
    let num = r[k];
    for (let j = 1; j < k; j++) num -= phi[j - 1] * r[k - j];
    const a = v > 1e-300 ? num / v : 0;
    const next = new Array(k);
    for (let j = 1; j < k; j++) next[j - 1] = phi[j - 1] - a * phi[k - j - 1];
    next[k - 1] = a;
    phi = next;
    v *= 1 - a * a;
    pac[k] = a;
  }
  return { pacf: pac, phi };
}

/**
 * Sample partial autocorrelation function by the Durbin–Levinson recursion.
 * Reference: Box, Jenkins & Reinsel, 4th ed., §3.2.6 and App. A3.2. Validity: 0 <= maxLag < n.
 * @param {ArrayLike<number>} x series [u]
 * @param {number} maxLag largest lag [samples]
 * @returns {number[]} φ_00..φ_KK (first entry 1 by convention) [-]
 */
export function pacf(x, maxLag) {
  return durbinLevinson(acf(x, maxLag), maxLag).pacf;
}

/**
 * Yule–Walker estimate of AR(p) coefficients (x_t - μ = Σ φ_i (x_{t-i} - μ) + e_t).
 * Reference: Box, Jenkins & Reinsel, 4th ed., §3.2.2.
 * @param {ArrayLike<number>} x series [u]
 * @param {number} p order
 * @returns {number[]} φ_1..φ_p [-]
 */
export function yuleWalker(x, p) {
  if (p === 0) return [];
  return durbinLevinson(acf(x, p), p).phi;
}

/**
 * d-th order backward difference ∇^d x (length n - d).
 * @param {ArrayLike<number>} x series [u]
 * @param {number} [d=1] order of differencing (integer >= 0)
 * @returns {number[]} differenced series [u]
 */
export function difference(x, d = 1) {
  needSeries(x, 'difference');
  if (!Number.isInteger(d) || d < 0 || d >= x.length) throw new RangeError('difference: d must be an integer within 0..n-1');
  let w = Array.from(x);
  for (let k = 0; k < d; k++) {
    const nw = new Array(w.length - 1);
    for (let i = 1; i < w.length; i++) nw[i - 1] = w[i] - w[i - 1];
    w = nw;
  }
  return w;
}

/**
 * Centred moving average of window w (same length as the input; the window is truncated symmetrically near the
 * ends so no phase shift is introduced). For even w the classical 2×w centred average (half weights at the two
 * extreme points) is used.
 * Reference: Brockwell & Davis, "Introduction to Time Series and Forecasting", §1.5.
 * @param {ArrayLike<number>} x series [u]
 * @param {number} w window length [samples], integer >= 1
 * @returns {number[]} smoothed series [u]
 */
export function movingAverage(x, w) {
  needSeries(x, 'movingAverage');
  if (!Number.isInteger(w) || w < 1) throw new RangeError('movingAverage: w must be a positive integer');
  const n = x.length;
  const out = new Array(n);
  const cs = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) cs[i + 1] = cs[i] + x[i];
  const even = w % 2 === 0;
  const hFull = even ? w / 2 : (w - 1) / 2;
  for (let i = 0; i < n; i++) {
    const h = Math.min(hFull, i, n - 1 - i);
    if (even && h === hFull && h > 0) {
      const inner = cs[i + h] - cs[i - h + 1];
      out[i] = (inner + 0.5 * (x[i - h] + x[i + h])) / (2 * h);
    } else {
      out[i] = (cs[i + h + 1] - cs[i - h]) / (2 * h + 1);
    }
  }
  return out;
}

/**
 * Exponential smoothing: simple (level), Holt (level + linear trend) or Holt–Winters (level + trend + season).
 *
 * The model is selected by the options: `season` > 1 activates the seasonal component (with trend unless
 * `trend: false`); `beta` (or `trend: true`) activates the trend. Any active smoothing parameter that is omitted
 * is estimated by minimising the one-step-ahead sum of squared errors (bounded Nelder–Mead on [0, 1]).
 * Recursions (additive): l_t = α(x_t - s_{t-m}) + (1-α)(l_{t-1} + b_{t-1}); b_t = β(l_t - l_{t-1}) + (1-β)b_{t-1};
 * s_t = γ(x_t - l_t) + (1-γ)s_{t-m}; multiplicative seasonality replaces differences by ratios.
 * Reference: Holt (1957); Winters (1960), Manage. Sci. 6, 324–342; Hyndman & Athanasopoulos, "Forecasting:
 * Principles and Practice", ch. 8. Validity: n >= 2 (>= 2 full seasons for the seasonal model); positive data
 * for multiplicative seasonality.
 * @param {ArrayLike<number>} x series [u]
 * @param {{ alpha?: number, beta?: number, gamma?: number, season?: number, trend?: boolean,
 *   seasonal?: 'additive' | 'multiplicative' }} [opts] smoothing constants in [0, 1]; season length [samples]
 * @returns {{ fitted: number[], forecast: (h: number) => number[], level: number, trend: number,
 *   seasonals: number[], alpha: number, beta: number | null, gamma: number | null, sse: number, rmse: number }}
 *   fitted: one-step-ahead predictions aligned with x [u]; forecast(h): next h values [u]; final states;
 *   parameters used; sum of squared one-step errors [u²] and its root mean [u]
 */
export function expSmoothing(x, opts = {}) {
  needSeries(x, 'expSmoothing', 2);
  const n = x.length;
  const m = opts.season && opts.season > 1 ? Math.round(opts.season) : 0;
  const useTrend = opts.trend !== undefined ? !!opts.trend : opts.beta !== undefined || m > 0;
  const mult = opts.seasonal === 'multiplicative';
  if (m && n < 2 * m) throw new RangeError('expSmoothing: x must span at least two full seasons');
  for (const k of ['alpha', 'beta', 'gamma']) {
    if (opts[k] !== undefined && opts[k] !== null && !(opts[k] >= 0 && opts[k] <= 1)) throw new RangeError(`expSmoothing: ${k} must be within 0..1`);
  }
  const run = (alpha, beta, gamma, keep) => {
    let level;
    let trend = 0;
    let seas = null;
    let start;
    if (m) {
      let m1 = 0;
      let m2 = 0;
      for (let i = 0; i < m; i++) {
        m1 += x[i];
        m2 += x[m + i];
      }
      m1 /= m;
      m2 /= m;
      level = m1;
      trend = useTrend ? (m2 - m1) / m : 0;
      seas = new Array(m);
      for (let i = 0; i < m; i++) {
        const line = m1 + trend * (i - (m - 1) / 2); // trend line through the first season
        seas[i] = mult ? x[i] / line : x[i] - line;
      }
      // start the recursion at the beginning with the level moved back half a season
      level = m1 - (trend * (m - 1)) / 2 - trend;
      start = 0;
    } else if (useTrend) {
      level = x[0];
      trend = x[1] - x[0];
      start = 1;
    } else {
      level = x[0];
      start = 1;
    }
    const fitted = keep ? new Array(n) : null;
    if (keep) for (let i = 0; i < start; i++) fitted[i] = x[i];
    let sse = 0;
    for (let t = start; t < n; t++) {
      const s = m ? seas[t % m] : mult ? 1 : 0;
      const base = level + trend;
      const pred = m ? (mult ? base * s : base + s) : base;
      const err = x[t] - pred;
      sse += err * err;
      if (keep) fitted[t] = pred;
      const prev = level;
      if (m) {
        level = alpha * (mult ? x[t] / s : x[t] - s) + (1 - alpha) * base;
        seas[t % m] = gamma * (mult ? x[t] / level : x[t] - level) + (1 - gamma) * s;
      } else {
        level = alpha * x[t] + (1 - alpha) * base;
      }
      if (useTrend) trend = beta * (level - prev) + (1 - beta) * trend;
    }
    return { sse, fitted, level, trend, seas, count: n - start };
  };
  const given = { alpha: opts.alpha ?? null, beta: useTrend ? (opts.beta ?? null) : 0, gamma: m ? (opts.gamma ?? null) : 0 };
  const free = ['alpha', 'beta', 'gamma'].filter((k) => given[k] === null);
  const par = { ...given };
  if (free.length) {
    const start0 = { alpha: 0.3, beta: 0.1, gamma: 0.1 };
    const obj = (v) => {
      free.forEach((k, i) => {
        par[k] = v[i];
      });
      const s = run(par.alpha, par.beta, par.gamma, false).sse;
      return Number.isFinite(s) ? s : 1e300;
    };
    let best = null;
    for (const scale of [1, 2.2]) {
      const r = nelderMead(
        obj,
        free.map((k) => Math.min(start0[k] * scale, 0.9)),
        { bounds: free.map(() => [0, 1]), tol: 1e-10, maxIter: 600, step: 0.15 },
      );
      if (!best || r.f < best.f) best = r;
    }
    free.forEach((k, i) => {
      par[k] = best.x[i];
    });
  }
  const res = run(par.alpha, par.beta, par.gamma, true);
  const seasonals = m ? res.seas.slice() : [];
  const forecast = (h) => {
    if (!Number.isInteger(h) || h < 1) throw new RangeError('expSmoothing.forecast: h must be a positive integer');
    const out = new Array(h);
    for (let k = 1; k <= h; k++) {
      const base = res.level + (useTrend ? k * res.trend : 0);
      if (m) {
        const s = seasonals[(n + k - 1) % m];
        out[k - 1] = mult ? base * s : base + s;
      } else out[k - 1] = base;
    }
    return out;
  };
  return {
    fitted: res.fitted,
    forecast,
    level: res.level,
    trend: useTrend ? res.trend : 0,
    seasonals,
    alpha: par.alpha,
    beta: useTrend ? par.beta : null,
    gamma: m ? par.gamma : null,
    sse: res.sse,
    rmse: Math.sqrt(res.sse / Math.max(res.count, 1)),
  };
}

/**
 * ARIMA(p, d, q) estimated by conditional sum of squares (CSS).
 *
 * Model for w_t = ∇^d x_t: (w_t - μ) = Σ φ_i (w_{t-i} - μ) + e_t + Σ θ_j e_{t-j}. Residuals are computed
 * recursively from t = p with pre-sample innovations set to zero and the sum of squares is minimised by
 * Levenberg–Marquardt, started from the Yule–Walker AR estimate. σ² = SSE / (number of residuals);
 * AIC = N ln σ² + N(1 + ln 2π) + 2(k + 1), k = number of estimated coefficients.
 * Forecast variances use the ψ-weights of the integrated model, Var[e_h] = σ² Σ_{j<h} ψ_j², with 95 % limits
 * mean ± 1.959964·sd (Gaussian innovations).
 * Reference: Box, Jenkins & Reinsel, "Time Series Analysis", 4th ed., §7.1.2–7.1.3 (conditional sum of squares)
 * and §5.1–5.2 (forecasts, ψ-weights). Validity: n - d > p + q + 1; the estimate is not constrained to the
 * stationary/invertible region — check `stationary`.
 * @param {ArrayLike<number>} x series [u]
 * @param {{ p?: number, d?: number, q?: number, includeMean?: boolean }} [order] orders (default 1, 0, 0);
 *   includeMean: estimate μ (default true when d = 0, false otherwise)
 * @returns {{ phi: number[], theta: number[], sigma2: number, aic: number,
 *   forecast: (h: number) => { mean: number[], lo95: number[], hi95: number[], se: number[] }, mean: number,
 *   residuals: number[], se: number[], converged: boolean, stationary: boolean, p: number, d: number, q: number }}
 *   AR and MA coefficients [-], innovation variance [u²], AIC [-], forecast function (arrays of length h, [u]),
 *   estimated mean μ of the differenced series [u], CSS residuals [u], standard errors of [phi, theta, μ]
 */
export function arima(x, order = {}) {
  needSeries(x, 'arima', 3);
  const p = order.p ?? 1;
  const d = order.d ?? 0;
  const q = order.q ?? 0;
  for (const [k, v] of [['p', p], ['d', d], ['q', q]]) {
    if (!Number.isInteger(v) || v < 0) throw new RangeError(`arima: ${k} must be a non-negative integer`);
  }
  const w = difference(x, d);
  const n = w.length;
  const includeMean = order.includeMean ?? d === 0;
  const k = p + q + (includeMean ? 1 : 0);
  if (n - p <= k) throw new RangeError('arima: series too short for the requested orders');
  const N = n - p;
  const resid = (par, out) => {
    const mu = includeMean ? par[p + q] : 0;
    const e = out || new Float64Array(N);
    for (let t = p; t < n; t++) {
      let v = w[t] - mu;
      for (let i = 0; i < p; i++) v -= par[i] * (w[t - 1 - i] - mu);
      for (let j = 0; j < q; j++) {
        const idx = t - 1 - j - p;
        if (idx >= 0) v -= par[p + j] * e[idx];
      }
      // guard against explosive non-invertible trial values
      e[t - p] = Math.abs(v) < 1e150 ? v : 1e150;
    }
    return e;
  };
  const wMean = meanOf(w);
  const start = new Array(k).fill(0);
  if (p > 0) {
    const yw = yuleWalker(w, p);
    for (let i = 0; i < p; i++) start[i] = Number.isFinite(yw[i]) ? yw[i] : 0;
  }
  if (includeMean) start[p + q] = wMean;
  let par = start;
  let converged = true;
  let se = new Array(k).fill(0);
  if (k > 0) {
    const fit = levenbergMarquardt((v) => resid(v), start, { tol: 1e-12, maxIter: 300 });
    par = fit.x;
    converged = fit.converged;
    se = fit.se;
  }
  const e = resid(par);
  let sse = 0;
  for (let i = 0; i < N; i++) sse += e[i] * e[i];
  const sigma2 = sse / N;
  const phi = par.slice(0, p);
  const theta = par.slice(p, p + q);
  const mu = includeMean ? par[p + q] : 0;
  const aic = N * Math.log(sigma2) + N * (1 + Math.log(2 * Math.PI)) + 2 * (k + 1);

  // AR polynomial of the integrated model: (1 - Σ φ_i B^i)(1 - B)^d = 1 - Σ φ*_i B^i
  let poly = [1, ...phi.map((v) => -v)];
  for (let i = 0; i < d; i++) {
    const np = new Array(poly.length + 1).fill(0);
    for (let j = 0; j < poly.length; j++) {
      np[j] += poly[j];
      np[j + 1] -= poly[j];
    }
    poly = np;
  }
  const phiStar = poly.slice(1).map((v) => -v);

  // last value of each differencing level, for re-integration
  const lastLevels = [];
  {
    let lev = Array.from(x);
    for (let i = 0; i < d; i++) {
      lastLevels.push(lev[lev.length - 1]);
      const nl = new Array(lev.length - 1);
      for (let j = 1; j < lev.length; j++) nl[j - 1] = lev[j] - lev[j - 1];
      lev = nl;
    }
  }

  const forecast = (h) => {
    if (!Number.isInteger(h) || h < 1) throw new RangeError('arima.forecast: h must be a positive integer');
    const wf = new Array(h);
    for (let s = 0; s < h; s++) {
      let v = mu;
      for (let i = 0; i < p; i++) {
        const idx = n + s - 1 - i;
        const wv = idx >= n ? wf[idx - n] : idx >= 0 ? w[idx] : mu;
        v += phi[i] * (wv - mu);
      }
      for (let j = 0; j < q; j++) {
        const idx = n + s - 1 - j; // innovation index in w-time
        if (idx < n && idx - p >= 0) v += theta[j] * e[idx - p];
      }
      wf[s] = v;
    }
    const last = lastLevels.slice();
    const meanF = new Array(h);
    for (let s = 0; s < h; s++) {
      let v = wf[s];
      for (let lvl = d - 1; lvl >= 0; lvl--) {
        v = last[lvl] + v;
        last[lvl] = v;
      }
      meanF[s] = v;
    }
    const psi = new Array(h).fill(0);
    psi[0] = 1;
    for (let j = 1; j < h; j++) {
      let v = j <= q ? theta[j - 1] : 0;
      for (let i = 1; i <= Math.min(j, phiStar.length); i++) v += phiStar[i - 1] * psi[j - i];
      psi[j] = v;
    }
    const seF = new Array(h);
    const lo95 = new Array(h);
    const hi95 = new Array(h);
    let acc = 0;
    for (let s = 0; s < h; s++) {
      acc += psi[s] * psi[s];
      seF[s] = Math.sqrt(sigma2 * acc);
      lo95[s] = meanF[s] - 1.959963984540054 * seF[s];
      hi95[s] = meanF[s] + 1.959963984540054 * seF[s];
    }
    return { mean: meanF, lo95, hi95, se: seF };
  };

  // stationarity check through the reflection coefficients of φ (Schur–Cohn / step-down recursion)
  let stationary = true;
  {
    let a = phi.slice();
    for (let m = p; m >= 1 && stationary; m--) {
      const km = a[m - 1];
      if (!(Math.abs(km) < 1)) {
        stationary = false;
        break;
      }
      const prev = new Array(m - 1);
      for (let i = 0; i < m - 1; i++) prev[i] = (a[i] + km * a[m - 2 - i]) / (1 - km * km);
      a = prev;
    }
  }

  return { phi, theta, sigma2, aic, forecast, mean: mu, residuals: Array.from(e), se, converged, stationary, p, d, q };
}

// ---------------------------------------------------------------------------------------------------------------
// Kalman filter
// ---------------------------------------------------------------------------------------------------------------

function mmul(A, B) {
  const n = A.length;
  const m = B[0].length;
  const k = B.length;
  const C = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = new Array(m).fill(0);
    const Ai = A[i];
    for (let l = 0; l < k; l++) {
      const a = Ai[l];
      if (a === 0) continue;
      const Bl = B[l];
      for (let j = 0; j < m; j++) row[j] += a * Bl[j];
    }
    C[i] = row;
  }
  return C;
}

function mT(A) {
  const n = A.length;
  const m = A[0].length;
  const T = new Array(m);
  for (let j = 0; j < m; j++) {
    T[j] = new Array(n);
    for (let i = 0; i < n; i++) T[j][i] = A[i][j];
  }
  return T;
}

function asMatrix(M, name) {
  if (typeof M === 'number') return [[M]];
  if (!Array.isArray(M) && !ArrayBuffer.isView(M)) throw new RangeError(`kalmanFilter: ${name} must be a matrix`);
  if (typeof M[0] === 'number') return [Array.from(M)];
  return M.map((r) => Array.from(r));
}

/**
 * Linear discrete-time Kalman filter.
 *
 * State model x_k = F x_{k-1} + w, w ~ N(0, Q); observation z_k = H x_k + v, v ~ N(0, R). Each step performs the
 * time update (x⁻ = F x, P⁻ = F P Fᵀ + Q) followed by the measurement update with gain K = P⁻Hᵀ(H P⁻ Hᵀ + R)⁻¹
 * and the Joseph-form covariance update P = (I - K H) P⁻ (I - K H)ᵀ + K R Kᵀ. x0, P0 are the state estimate and
 * covariance before the first observation. An observation that is null or contains NaN is treated as missing
 * (time update only).
 * Reference: Kalman (1960), J. Basic Eng. 82, 35–45; Simon, "Optimal State Estimation" (2006), §5.1.
 * Validity: small dense systems; R (and S) positive definite.
 * @param {{ F: number[][], H: number[][], Q: number[][], R: number[][], x0: number[], P0: number[][] }} model
 *   state transition (n×n), observation matrix (m×n), process noise covariance (n×n), measurement noise
 *   covariance (m×m), initial state (n) and covariance (n×n)
 * @param {(number | number[] | null)[]} observations measurement vectors (scalars allowed when m = 1)
 * @returns {{ x: number[][], P: number[][][], innovations: number[][], xPred: number[][], S: number[][][],
 *   logLik: number }} filtered states per step, filtered covariances, innovations z - H x⁻, predicted states,
 *   innovation covariances, Gaussian log-likelihood of the observations
 */
export function kalmanFilter(model, observations) {
  if (!model) throw new RangeError('kalmanFilter: model is required');
  const F = asMatrix(model.F, 'F');
  const H = asMatrix(model.H, 'H');
  const n = F.length;
  const m = H.length;
  const Q = asMatrix(model.Q, 'Q');
  const R = asMatrix(model.R, 'R');
  if (F[0].length !== n) throw new RangeError('kalmanFilter: F must be square');
  if (H[0].length !== n) throw new RangeError('kalmanFilter: H must have as many columns as states');
  if (Q.length !== n || Q[0].length !== n) throw new RangeError('kalmanFilter: Q must be n×n');
  if (R.length !== m || R[0].length !== m) throw new RangeError('kalmanFilter: R must be m×m');
  let x = typeof model.x0 === 'number' ? [model.x0] : Array.from(model.x0 || []);
  if (x.length !== n) throw new RangeError('kalmanFilter: x0 must have one entry per state');
  let P = asMatrix(model.P0, 'P0');
  if (P.length !== n || P[0].length !== n) throw new RangeError('kalmanFilter: P0 must be n×n');
  if (!observations || typeof observations.length !== 'number') throw new RangeError('kalmanFilter: observations must be an array');
  const Ft = mT(F);
  const Ht = mT(H);
  const xs = [];
  const Ps = [];
  const inns = [];
  const xPred = [];
  const Ss = [];
  let logLik = 0;
  for (let kStep = 0; kStep < observations.length; kStep++) {
    // time update
    const xp = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) xp[i] += F[i][j] * x[j];
    const Pp = mmul(mmul(F, P), Ft);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Pp[i][j] += Q[i][j];
    xPred.push(xp);
    let z = observations[kStep];
    if (typeof z === 'number') z = [z];
    let missing = z === null || z === undefined;
    if (!missing) {
      if (z.length !== m) throw new RangeError(`kalmanFilter: observation ${kStep} must have ${m} entries`);
      for (let i = 0; i < m; i++) if (!Number.isFinite(z[i])) missing = true;
    }
    if (missing) {
      x = xp;
      P = Pp;
      xs.push(x.slice());
      Ps.push(P.map((r) => r.slice()));
      inns.push(new Array(m).fill(NaN));
      Ss.push(null);
      continue;
    }
    const y = new Array(m);
    for (let i = 0; i < m; i++) {
      let s = z[i];
      for (let j = 0; j < n; j++) s -= H[i][j] * xp[j];
      y[i] = s;
    }
    const PHt = mmul(Pp, Ht);
    const S = mmul(H, PHt);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) S[i][j] += R[i][j];
    let Sinv;
    let detS;
    if (m === 1) {
      if (!(S[0][0] > 0)) throw new RangeError('kalmanFilter: innovation variance is not positive');
      Sinv = [[1 / S[0][0]]];
      detS = S[0][0];
    } else {
      Sinv = inverse(S);
      detS = detSmall(S);
    }
    const K = mmul(PHt, Sinv);
    x = xp.slice();
    for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) x[i] += K[i][j] * y[j];
    const KH = mmul(K, H);
    const A = new Array(n);
    for (let i = 0; i < n; i++) {
      A[i] = new Array(n);
      for (let j = 0; j < n; j++) A[i][j] = (i === j ? 1 : 0) - KH[i][j];
    }
    const Pn = mmul(mmul(A, Pp), mT(A));
    const KRKt = mmul(mmul(K, R), mT(K));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Pn[i][j] += KRKt[i][j];
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const s = 0.5 * (Pn[i][j] + Pn[j][i]);
        Pn[i][j] = s;
        Pn[j][i] = s;
      }
    }
    P = Pn;
    let quad = 0;
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) quad += y[i] * Sinv[i][j] * y[j];
    logLik += -0.5 * (m * Math.log(2 * Math.PI) + Math.log(detS) + quad);
    xs.push(x.slice());
    Ps.push(P.map((r) => r.slice()));
    inns.push(y);
    Ss.push(S);
  }
  return { x: xs, P: Ps, innovations: inns, xPred, S: Ss, logLik };
}

/** Determinant of a small dense matrix by Gaussian elimination with partial pivoting. */
function detSmall(M) {
  const n = M.length;
  const A = M.map((r) => r.slice());
  let det = 1;
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
    if (A[piv][c] === 0) return 0;
    if (piv !== c) {
      const t = A[piv];
      A[piv] = A[c];
      A[c] = t;
      det = -det;
    }
    det *= A[c][c];
    for (let r = c + 1; r < n; r++) {
      const f = A[r][c] / A[c][c];
      for (let j = c; j < n; j++) A[r][j] -= f * A[c][j];
    }
  }
  return det;
}
