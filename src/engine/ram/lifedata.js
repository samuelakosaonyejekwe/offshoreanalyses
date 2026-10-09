// Life-data analysis: distribution fitting with right-censored data, probability plotting, hazard models,
// renewal theory, repairable-system trend analysis and spares provisioning.
// Times may be in any consistent unit.

import { erfc, erfinv, gamma } from '../math/special.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function isSeries(x) {
  return Array.isArray(x) || ArrayBuffer.isView(x);
}
const Phi = (z) => 0.5 * erfc(-z / Math.SQRT2);
const PhiInv = (p) => Math.SQRT2 * erfinv(2 * p - 1);
const phi = (z) => Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);

function prepare(times, censored) {
  need(isSeries(times) && times.length > 0, 'times must be a non-empty array');
  if (censored !== undefined) need(isSeries(censored) && censored.length === times.length, 'censored must match times in length');
  let r = 0;
  for (let i = 0; i < times.length; i++) {
    need(times[i] > 0 && Number.isFinite(times[i]), `times[${i}] must be positive and finite`);
    if (!censored || !censored[i]) r++;
  }
  need(r > 0, 'at least one failure (uncensored time) is required');
  return r;
}

/**
 * Maximum-likelihood fit of a two-parameter Weibull distribution to complete or right-censored life data.
 * The shape β solves  Σ tᵢ^β·ln tᵢ / Σ tᵢ^β − 1/β − (1/r)·Σ_failures ln tᵢ = 0  (sums over all units, r failures);
 * the scale is η = (Σ tᵢ^β / r)^(1/β).
 * Reference: Cohen (1965) Technometrics 7, 579–588; Meeker & Escobar, Statistical Methods for Reliability Data
 * (1998); IEC 61649.
 * @param {object} o
 * @param {ArrayLike<number>} o.times failure or suspension times [time]
 * @param {ArrayLike<boolean>} [o.censored] true where the unit was suspended (right-censored)
 * @returns {{shape:number, scale:number, logLikelihood:number, failures:number, n:number, mean:number,
 *   seShape:number, seScale:number}} standard errors from the observed information matrix
 */
export function fitWeibullMle({ times, censored }) {
  const r = prepare(times, censored);
  const n = times.length;
  let tMax = 0;
  let sumLnFail = 0;
  for (let i = 0; i < n; i++) {
    if (times[i] > tMax) tMax = times[i];
    if (!censored || !censored[i]) sumLnFail += Math.log(times[i]);
  }
  const x = Float64Array.from(times, (t) => t / tMax);
  const lnx = x.map(Math.log);
  const meanLn = sumLnFail / r - Math.log(tMax);
  const g = (b) => {
    let s0 = 0;
    let s1 = 0;
    for (let i = 0; i < n; i++) {
      const p = x[i] ** b;
      s0 += p;
      s1 += p * lnx[i];
    }
    return s1 / s0 - 1 / b - meanLn;
  };
  let lo = 1e-3;
  let hi = 1;
  while (g(hi) < 0 && hi < 1e4) hi *= 2;
  need(g(hi) >= 0, 'times: Weibull shape did not converge (degenerate data)');
  for (let i = 0; i < 200 && hi - lo > 1e-13 * hi; i++) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) < 0) lo = mid; else hi = mid;
  }
  const shape = 0.5 * (lo + hi);
  let s0 = 0;
  for (let i = 0; i < n; i++) s0 += x[i] ** shape;
  const scale = tMax * (s0 / r) ** (1 / shape);
  const ll = (b, e) => {
    let l = 0;
    for (let i = 0; i < n; i++) {
      const z = (times[i] / e) ** b;
      l -= z;
      if (!censored || !censored[i]) l += Math.log(b / e) + (b - 1) * Math.log(times[i] / e);
    }
    return l;
  };
  // Observed information by central differences in (ln β, ln η), mapped back by the delta method.
  const h = 1e-4;
  const f = (u, v) => ll(shape * Math.exp(u), scale * Math.exp(v));
  const f0 = f(0, 0);
  const huu = (f(h, 0) - 2 * f0 + f(-h, 0)) / (h * h);
  const hvv = (f(0, h) - 2 * f0 + f(0, -h)) / (h * h);
  const huv = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h);
  const det = huu * hvv - huv * huv;
  const seShape = det > 0 ? shape * Math.sqrt(-hvv / det) : NaN;
  const seScale = det > 0 ? scale * Math.sqrt(-huu / det) : NaN;
  return { shape, scale, logLikelihood: f0, failures: r, n, mean: scale * gamma(1 + 1 / shape), seShape, seScale };
}

/**
 * Maximum-likelihood fit of an exponential distribution with right censoring: λ = r / Σ tᵢ (total time on test).
 * Reference: Epstein & Sobel (1953); Meeker & Escobar.
 * @param {object} o
 * @param {ArrayLike<number>} o.times [time]
 * @param {ArrayLike<boolean>} [o.censored]
 * @returns {{rate:number, mtbf:number, failures:number, totalTime:number, logLikelihood:number}}
 */
export function fitExponentialMle({ times, censored }) {
  const r = prepare(times, censored);
  let T = 0;
  for (let i = 0; i < times.length; i++) T += times[i];
  const rate = r / T;
  return { rate, mtbf: 1 / rate, failures: r, totalTime: T, logLikelihood: r * Math.log(rate) - rate * T };
}

/**
 * Maximum-likelihood fit of a lognormal distribution with right censoring, by the EM algorithm on ln t
 * (censored observations are replaced by their conditional first and second moments each iteration).
 * Reference: Dempster, Laird & Rubin (1977) JRSS B 39, 1–38; Meeker & Escobar.
 * @param {object} o
 * @param {ArrayLike<number>} o.times [time]
 * @param {ArrayLike<boolean>} [o.censored]
 * @returns {{mu:number, sigma:number, median:number, mean:number, failures:number, iterations:number}} mu and sigma
 *   of ln t
 */
export function fitLognormalMle({ times, censored }) {
  const r = prepare(times, censored);
  const n = times.length;
  const y = Float64Array.from(times, Math.log);
  let mu = 0;
  for (let i = 0; i < n; i++) mu += y[i];
  mu /= n;
  let s2 = 0;
  for (let i = 0; i < n; i++) s2 += (y[i] - mu) ** 2;
  let sigma = Math.sqrt(s2 / n) || 1;
  need(r >= 2 || sigma > 0, 'times: at least two distinct failures are required');
  let it = 0;
  for (; it < 20000; it++) {
    let m1 = 0;
    let m2 = 0;
    for (let i = 0; i < n; i++) {
      if (censored && censored[i]) {
        const z = (y[i] - mu) / sigma;
        const surv = 1 - Phi(z);
        const hz = surv > 1e-300 ? phi(z) / surv : z;
        m1 += mu + sigma * hz;
        m2 += mu * mu + sigma * sigma + sigma * (y[i] + mu) * hz;
      } else {
        m1 += y[i];
        m2 += y[i] * y[i];
      }
    }
    const muNew = m1 / n;
    const sigNew = Math.sqrt(Math.max(m2 / n - muNew * muNew, 1e-300));
    const d = Math.abs(muNew - mu) + Math.abs(sigNew - sigma);
    mu = muNew;
    sigma = sigNew;
    if (d < 1e-12 * (1 + Math.abs(mu))) break;
  }
  return { mu, sigma, median: Math.exp(mu), mean: Math.exp(mu + 0.5 * sigma * sigma), failures: r, iterations: it };
}

/**
 * Median ranks with Johnson's adjustment for suspensions and Benard's approximation F = (j − 0.3)/(n + 0.4).
 * Reference: Johnson (1964) The Statistical Treatment of Fatigue Experiments; Benard & Bos-Levenbach (1953);
 * IEC 61649; Abernethy, The New Weibull Handbook.
 * @param {object} o
 * @param {ArrayLike<number>} o.times [time]
 * @param {ArrayLike<boolean>} [o.censored]
 * @returns {{t:number[], rank:number[], F:number[]}} failure times in ascending order with adjusted rank and
 *   median-rank estimate of the unreliability
 */
export function medianRanks({ times, censored }) {
  prepare(times, censored);
  const n = times.length;
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => times[a] - times[b] || ((censored && censored[a]) ? 1 : 0) - ((censored && censored[b]) ? 1 : 0));
  const t = [];
  const rank = [];
  const F = [];
  let prev = 0;
  for (let k = 0; k < n; k++) {
    const i = order[k];
    if (censored && censored[i]) continue;
    const remaining = n - k;
    const adj = prev + (n + 1 - prev) / (1 + remaining);
    prev = adj;
    t.push(times[i]);
    rank.push(adj);
    F.push((adj - 0.3) / (n + 0.4));
  }
  return { t, rank, F };
}

/**
 * Median-rank regression (probability-plot fit) with the plot coordinates.
 *  Weibull:     x = ln t, y = ln(−ln(1 − F));  β = slope, η = exp(−intercept/β)
 *  lognormal:   x = ln t, y = Φ⁻¹(F);          σ = 1/slope, μ = −intercept/slope
 *  exponential: x = t,    y = −ln(1 − F);      λ = slope of the line through the origin
 * Regression of y on x (least squares).
 * Reference: IEC 61649; Abernethy, The New Weibull Handbook.
 * @param {object} o
 * @param {ArrayLike<number>} o.times [time]
 * @param {ArrayLike<boolean>} [o.censored]
 * @param {string} [o.dist='weibull'] 'weibull' | 'lognormal' | 'exponential'
 * @returns {{dist:string, shape?:number, scale?:number, mu?:number, sigma?:number, rate?:number, r2:number,
 *   plot: {x:number[], y:number[], t:number[], F:number[]}}}
 */
export function medianRankRegression({ times, censored, dist = 'weibull' }) {
  const mr = medianRanks({ times, censored });
  const m = mr.t.length;
  need(m >= 2, 'at least two failures are required for regression');
  const x = [];
  const y = [];
  for (let i = 0; i < m; i++) {
    if (dist === 'weibull') { x.push(Math.log(mr.t[i])); y.push(Math.log(-Math.log(1 - mr.F[i]))); }
    else if (dist === 'lognormal') { x.push(Math.log(mr.t[i])); y.push(PhiInv(mr.F[i])); }
    else if (dist === 'exponential') { x.push(mr.t[i]); y.push(-Math.log(1 - mr.F[i])); }
    else throw new RangeError(`dist: unknown distribution "${dist}"`);
  }
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0; let syy = 0;
  for (let i = 0; i < m; i++) { sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; syy += y[i] * y[i]; }
  const plot = { x, y, t: mr.t, F: mr.F };
  if (dist === 'exponential') {
    const rate = sxy / sxx;
    let sse = 0;
    for (let i = 0; i < m; i++) sse += (y[i] - rate * x[i]) ** 2;
    const sst = syy - (sy * sy) / m;
    return { dist, rate, r2: sst > 0 ? 1 - sse / sst : 1, plot };
  }
  const den = m * sxx - sx * sx;
  need(den > 0, 'times: failure times must not all be equal');
  const slope = (m * sxy - sx * sy) / den;
  const intercept = (sy - slope * sx) / m;
  const r2 = (m * sxy - sx * sy) ** 2 / (den * (m * syy - sy * sy));
  if (dist === 'weibull') return { dist, shape: slope, scale: Math.exp(-intercept / slope), r2, plot };
  return { dist, sigma: 1 / slope, mu: -intercept / slope, r2, plot };
}

/**
 * Weibull reliability functions at time t: R = exp(−(t/η)^β), F = 1 − R, pdf, hazard h = (β/η)(t/η)^(β−1).
 * @param {number} t time [time]
 * @param {{shape:number, scale:number}} par
 * @returns {{reliability:number, cdf:number, pdf:number, hazard:number}}
 */
export function weibullFunctions(t, { shape, scale }) {
  need(t >= 0, 't must be non-negative');
  need(shape > 0 && scale > 0, 'shape and scale must be positive');
  const z = (t / scale) ** shape;
  const R = Math.exp(-z);
  const hazard = t === 0 ? (shape < 1 ? Infinity : shape === 1 ? 1 / scale : 0) : (shape / scale) * (t / scale) ** (shape - 1);
  return { reliability: R, cdf: 1 - R, pdf: Number.isFinite(hazard) ? hazard * R : Infinity, hazard };
}

/**
 * B-life: the time by which a fraction p of the population has failed, B_p = η·(−ln(1 − p))^(1/β)
 * (B10: p = 0.10).
 * Reference: IEC 61649.
 * @param {{shape:number, scale:number}} par Weibull parameters
 * @param {number} [p=0.10] fraction failed [-]
 * @returns {number} B-life [time]
 */
export function bLife({ shape, scale }, p = 0.10) {
  need(shape > 0 && scale > 0, 'shape and scale must be positive');
  need(p > 0 && p < 1, 'p must be in (0, 1)');
  return scale * (-Math.log(1 - p)) ** (1 / shape);
}

/**
 * Bathtub hazard as a competing-risk composite of Weibull modes (infant mortality β < 1, random β = 1,
 * wear-out β > 1): h(t) = Σ (βᵢ/ηᵢ)(t/ηᵢ)^(βᵢ−1), R(t) = exp(−Σ (t/ηᵢ)^βᵢ).
 * Reference: Meeker & Escobar, Statistical Methods for Reliability Data (competing risks).
 * @param {number} t time [time], > 0
 * @param {{shape:number, scale:number}[]} modes
 * @returns {{hazard:number, reliability:number, contributions:number[]}} hazard [1/time]
 */
export function bathtubHazard(t, modes) {
  need(t > 0, 't must be positive');
  need(Array.isArray(modes) && modes.length > 0, 'modes must be a non-empty array');
  let h = 0;
  let H = 0;
  const contributions = modes.map((m) => {
    need(m.shape > 0 && m.scale > 0, 'mode shape and scale must be positive');
    const c = (m.shape / m.scale) * (t / m.scale) ** (m.shape - 1);
    h += c;
    H += (t / m.scale) ** m.shape;
    return c;
  });
  return { hazard: h, reliability: Math.exp(-H), contributions };
}

/**
 * Renewal function M(t) = E[number of renewals in (0, t]] by discretising the renewal equation
 * M(t) = F(t) + ∫₀ᵗ M(t − x) dF(x) (Riemann–Stieltjes scheme, second-order accurate).
 * Reference: Xie (1989) "On the solution of renewal-type integral equations", Commun. Statist. – Simul. 18(1),
 * 281–293; Barlow & Proschan, Mathematical Theory of Reliability.
 * @param {object} o
 * @param {function(number):number} o.cdf life distribution F(t)
 * @param {number} o.tMax largest time [time]
 * @param {number} [o.steps=1000] number of intervals [-]
 * @returns {{t: Float64Array, M: Float64Array}} renewal function on the grid
 */
export function renewalFunction({ cdf, tMax, steps = 1000 }) {
  need(typeof cdf === 'function', 'cdf must be a function');
  need(tMax > 0, 'tMax must be positive');
  need(Number.isInteger(steps) && steps > 0, 'steps must be a positive integer');
  const h = tMax / steps;
  const t = new Float64Array(steps + 1);
  const F = new Float64Array(steps + 1);
  const M = new Float64Array(steps + 1);
  for (let i = 0; i <= steps; i++) { t[i] = i * h; F[i] = cdf(t[i]); }
  const dF1 = F[1] - F[0];
  for (let i = 1; i <= steps; i++) {
    let s = F[i] + 0.5 * M[i - 1] * dF1;
    for (let j = 2; j <= i; j++) s += 0.5 * (M[i - j + 1] + M[i - j]) * (F[j] - F[j - 1]);
    M[i] = s / (1 - 0.5 * dF1);
  }
  return { t, M };
}

/**
 * Renewal function of a Weibull life distribution at time t (convenience wrapper).
 * @param {object} o
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.t [time]
 * @param {number} [o.steps=1000]
 * @returns {number} expected number of renewals [-]
 */
export function weibullRenewals({ shape, scale, t, steps = 1000 }) {
  need(shape > 0 && scale > 0, 'shape and scale must be positive');
  if (t === 0) return 0;
  const r = renewalFunction({ cdf: (x) => 1 - Math.exp(-((x / scale) ** shape)), tMax: t, steps });
  return r.M[steps];
}

/**
 * Crow–AMSAA (power-law non-homogeneous Poisson process) fit for a repairable system:
 * E[N(t)] = λ·t^β, intensity u(t) = λ·β·t^(β−1). Maximum likelihood:
 *   time-terminated at T:   β = n / Σ ln(T/tᵢ),            λ = n / T^β
 *   failure-terminated:     β = n / Σ_{i<n} ln(t_n/tᵢ),    λ = n / t_n^β
 * β < 1: reliability growth; β > 1: deterioration.
 * Reference: Crow (1974) "Reliability analysis for complex repairable systems"; MIL-HDBK-189C; IEC 61164.
 * @param {object} o
 * @param {ArrayLike<number>} o.times cumulative failure times of one system, ascending [time]
 * @param {number} [o.endTime] end of observation T [time]; omitted = failure-terminated at the last failure
 * @returns {{beta:number, lambda:number, intensity:number, instantaneousMtbf:number, cumulativeMtbf:number,
 *   n:number, endTime:number, expected:(function(number):number)}} intensity and MTBF at the end time
 */
export function crowAmsaa({ times, endTime }) {
  need(isSeries(times) && times.length >= 2, 'times must contain at least two failures');
  for (let i = 0; i < times.length; i++) {
    need(times[i] > 0, `times[${i}] must be positive`);
    if (i > 0) need(times[i] >= times[i - 1], 'times must be ascending cumulative times');
  }
  const n = times.length;
  const T = endTime === undefined ? times[n - 1] : endTime;
  need(T >= times[n - 1], 'endTime must not precede the last failure');
  let s = 0;
  for (let i = 0; i < n; i++) s += Math.log(T / times[i]);
  need(s > 0, 'times: failures must not all coincide with the end time');
  const beta = n / s;
  const lambda = n / T ** beta;
  const intensity = lambda * beta * T ** (beta - 1);
  return { beta, lambda, intensity, instantaneousMtbf: 1 / intensity, cumulativeMtbf: T / n, n, endTime: T, expected: (t) => lambda * t ** beta };
}

/**
 * Laplace trend test for a repairable system: U = (Σtᵢ/n − T/2)/(T·√(1/(12n))) is asymptotically standard
 * normal under a homogeneous Poisson process. U > 0: failures concentrate late (deterioration); U < 0: improvement.
 * For failure-terminated data the last failure is the end time and is excluded from the sum.
 * Reference: MIL-HDBK-189C; Ascher & Feingold, Repairable Systems Reliability (1984).
 * @param {object} o
 * @param {ArrayLike<number>} o.times cumulative failure times [time]
 * @param {number} [o.endTime] end of observation [time]; omitted = failure-terminated
 * @returns {{U:number, pValue:number, trend:string}} two-sided p-value; trend 'deteriorating' | 'improving' |
 *   'none' at the 5 % level
 */
export function laplaceTest({ times, endTime }) {
  need(isSeries(times) && times.length >= 2, 'times must contain at least two failures');
  let n = times.length;
  let T = endTime;
  if (T === undefined) {
    T = times[n - 1];
    n -= 1;
  }
  need(T > 0, 'endTime must be positive');
  let s = 0;
  for (let i = 0; i < n; i++) s += times[i];
  const U = (s / n - T / 2) / (T * Math.sqrt(1 / (12 * n)));
  const pValue = erfc(Math.abs(U) / Math.SQRT2);
  return { U, pValue, trend: pValue < 0.05 ? (U > 0 ? 'deteriorating' : 'improving') : 'none' };
}

/**
 * MTBF, MTTR and availability from operating records or from summary values.
 *   inherent availability A = MTBF/(MTBF + MTTR); operational availability includes the logistic delay.
 * Reference: IEC 60050-192; MIL-HDBK-338B.
 * @param {object} o
 * @param {ArrayLike<number>} [o.upTimes] observed times between failures [time]
 * @param {ArrayLike<number>} [o.downTimes] observed repair (down) times [time]
 * @param {number} [o.mtbf] [time] (instead of upTimes)
 * @param {number} [o.mttr] [time] (instead of downTimes)
 * @param {number} [o.logisticDelay=0] mean waiting time for spares/access per failure [time]
 * @returns {{mtbf:number, mttr:number, failureRate:number, repairRate:number, availability:number,
 *   operationalAvailability:number}}
 */
export function availabilityMetrics({ upTimes, downTimes, mtbf, mttr, logisticDelay = 0 }) {
  const mean = (a, name) => {
    need(isSeries(a) && a.length > 0, `${name} must be a non-empty array`);
    let s = 0;
    for (let i = 0; i < a.length; i++) { need(a[i] >= 0, `${name} must be non-negative`); s += a[i]; }
    return s / a.length;
  };
  const up = mtbf === undefined ? mean(upTimes, 'upTimes') : mtbf;
  const down = mttr === undefined ? mean(downTimes, 'downTimes') : mttr;
  need(up > 0, 'mtbf must be positive');
  need(down >= 0, 'mttr must be non-negative');
  need(logisticDelay >= 0, 'logisticDelay must be non-negative');
  return {
    mtbf: up, mttr: down, failureRate: 1 / up, repairRate: down > 0 ? 1 / down : Infinity,
    availability: up / (up + down), operationalAvailability: up / (up + down + logisticDelay),
  };
}

/**
 * Spare-parts quantity for a service level: the smallest s with P[Poisson(λ·n·T) <= s] >= level, where the
 * demand during the replenishment (or mission) time T comes from n units failing at rate λ each.
 * Reference: Sherbrooke, Optimal Inventory Modeling of Systems (2nd ed.); MIL-HDBK-338B.
 * @param {object} o
 * @param {number} o.failureRate failure rate per unit [1/time]
 * @param {number} [o.units=1] number of units in service [-]
 * @param {number} o.period replenishment lead time or support period [time]
 * @param {number} [o.serviceLevel=0.95] required probability of no stock-out [-]
 * @returns {{spares:number, expectedDemand:number, achievedLevel:number}}
 */
export function sparesForServiceLevel({ failureRate, units = 1, period, serviceLevel = 0.95 }) {
  need(failureRate >= 0, 'failureRate must be non-negative');
  need(units > 0, 'units must be positive');
  need(period >= 0, 'period must be non-negative');
  need(serviceLevel > 0 && serviceLevel < 1, 'serviceLevel must be in (0, 1)');
  const m = failureRate * units * period;
  need(m < 700, 'expected demand is too large for the direct Poisson sum');
  let term = Math.exp(-m);
  let cum = term;
  let s = 0;
  while (cum < serviceLevel && s < 100000) {
    s++;
    term *= m / s;
    cum += term;
  }
  return { spares: s, expectedDemand: m, achievedLevel: Math.min(1, cum) };
}
