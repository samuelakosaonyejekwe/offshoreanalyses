/**
 * Descriptive statistics, probability distributions, regression and goodness-of-fit.
 *
 * All functions accept plain arrays or typed arrays. Quantities are unit-agnostic: results carry the units of the
 * data (moments) or are dimensionless (probabilities, correlation, skewness, kurtosis).
 */
import { erfc, lgamma, gammaP, gammaQ, betaI } from './special.js';
import { leastSquares, inverse } from './linalg.js';

const SQRT2 = Math.SQRT2;
const SQRT2PI = Math.sqrt(2 * Math.PI);
/** Euler–Mascheroni constant. */
export const EULER_GAMMA = 0.5772156649015329;

function need(arr, name, nMin = 1) {
  if (!arr || typeof arr.length !== 'number' || arr.length < nMin) {
    throw new RangeError(`${name}: needs an array with at least ${nMin} value(s)`);
  }
}

// ---------------------------------------------------------------------------------------------------------------
// Descriptive statistics
// ---------------------------------------------------------------------------------------------------------------

/**
 * Arithmetic mean. Reference: definition. Validity: n >= 1.
 * @param {ArrayLike<number>} arr data [u]
 * @returns {number} mean [u]
 */
export function mean(arr) {
  need(arr, 'mean');
  let s = 0;
  for (let i = 0; i < arr.length; i++) s += arr[i];
  return s / arr.length;
}

/**
 * Unbiased sample variance (divisor n-1; n for `population = true`). Validity: n >= 2 (n = 1 returns 0).
 * @param {ArrayLike<number>} arr data [u]
 * @param {boolean} [population=false] use divisor n
 * @returns {number} variance [u²]
 */
export function variance(arr, population = false) {
  need(arr, 'variance');
  const n = arr.length;
  if (n === 1) return 0;
  const m = mean(arr);
  let s = 0;
  for (let i = 0; i < n; i++) {
    const d = arr[i] - m;
    s += d * d;
  }
  return s / (population ? n : n - 1);
}

/**
 * Sample standard deviation (divisor n-1 unless `population`).
 * @param {ArrayLike<number>} arr data [u]
 * @param {boolean} [population=false] use divisor n
 * @returns {number} standard deviation [u]
 */
export function std(arr, population = false) {
  return Math.sqrt(variance(arr, population));
}

/**
 * Minimum (loop based, safe for very long arrays).
 * @param {ArrayLike<number>} arr data [u]
 * @returns {number} smallest value [u]
 */
export function min(arr) {
  need(arr, 'min');
  let m = Infinity;
  for (let i = 0; i < arr.length; i++) if (arr[i] < m) m = arr[i];
  return m;
}

/**
 * Maximum (loop based, safe for very long arrays).
 * @param {ArrayLike<number>} arr data [u]
 * @returns {number} largest value [u]
 */
export function max(arr) {
  need(arr, 'max');
  let m = -Infinity;
  for (let i = 0; i < arr.length; i++) if (arr[i] > m) m = arr[i];
  return m;
}

function sortedCopy(arr) {
  return Float64Array.from(arr).sort();
}

/**
 * Percentile of an already ascending-sorted array, linear interpolation between order statistics
 * (Hyndman & Fan 1996 type 7, the spreadsheet / NumPy default).
 * @param {ArrayLike<number>} sorted ascending data [u]
 * @param {number} p percentile in 0..100 [%]
 * @returns {number} percentile [u]
 */
export function percentileSorted(sorted, p) {
  need(sorted, 'percentile');
  if (!(p >= 0 && p <= 100)) throw new RangeError('percentile: p must be within 0..100');
  const h = ((sorted.length - 1) * p) / 100;
  const i = Math.floor(h);
  const j = Math.min(i + 1, sorted.length - 1);
  return sorted[i] + (h - i) * (sorted[j] - sorted[i]);
}

/**
 * Percentile with linear interpolation (Hyndman & Fan 1996, type 7).
 * @param {ArrayLike<number>} arr data (any order) [u]
 * @param {number} p percentile in 0..100 [%]
 * @returns {number} percentile [u]
 */
export function percentile(arr, p) {
  need(arr, 'percentile');
  return percentileSorted(sortedCopy(arr), p);
}

/**
 * Median (50th percentile).
 * @param {ArrayLike<number>} arr data [u]
 * @returns {number} median [u]
 */
export function median(arr) {
  return percentile(arr, 50);
}

/**
 * Sample skewness g1 = m3 / m2^1.5 (biased moment estimator, as in most statistical packages by default).
 * Validity: n >= 2 and non-zero variance.
 * @param {ArrayLike<number>} arr data [u]
 * @returns {number} skewness [-]
 */
export function skewness(arr) {
  need(arr, 'skewness', 2);
  const n = arr.length;
  const m = mean(arr);
  let m2 = 0;
  let m3 = 0;
  for (let i = 0; i < n; i++) {
    const d = arr[i] - m;
    m2 += d * d;
    m3 += d * d * d;
  }
  m2 /= n;
  m3 /= n;
  return m2 > 0 ? m3 / Math.pow(m2, 1.5) : 0;
}

/**
 * Sample excess kurtosis g2 = m4 / m2² - 3 (0 for a normal distribution). Pass `excess = false` for the Pearson
 * kurtosis m4 / m2² (3 for a normal distribution). Validity: n >= 2 and non-zero variance.
 * @param {ArrayLike<number>} arr data [u]
 * @param {boolean} [excess=true] subtract 3
 * @returns {number} kurtosis [-]
 */
export function kurtosis(arr, excess = true) {
  need(arr, 'kurtosis', 2);
  const n = arr.length;
  const m = mean(arr);
  let m2 = 0;
  let m4 = 0;
  for (let i = 0; i < n; i++) {
    const d2 = (arr[i] - m) * (arr[i] - m);
    m2 += d2;
    m4 += d2 * d2;
  }
  m2 /= n;
  m4 /= n;
  if (!(m2 > 0)) return 0;
  return m4 / (m2 * m2) - (excess ? 3 : 0);
}

/**
 * Histogram with equal-width bins over [min, max] (or explicit ascending edges). The last bin is closed on the
 * right. Values outside explicit edges are ignored.
 * @param {ArrayLike<number>} arr data [u]
 * @param {number|number[]} [bins=10] number of bins, or array of bin edges [u]
 * @returns {{ edges: number[], counts: number[], centres: number[], density: number[] }}
 *   edges (nBins+1) [u], counts per bin [-], bin centres [u], probability density (integrates to 1) [1/u]
 */
export function histogram(arr, bins = 10) {
  need(arr, 'histogram');
  let edges;
  if (typeof bins === 'number') {
    if (!(bins >= 1) || !Number.isInteger(bins)) throw new RangeError('histogram: bins must be a positive integer');
    let lo = min(arr);
    let hi = max(arr);
    if (hi === lo) {
      lo -= 0.5;
      hi += 0.5;
    }
    edges = new Array(bins + 1);
    for (let i = 0; i <= bins; i++) edges[i] = lo + ((hi - lo) * i) / bins;
  } else {
    edges = Array.from(bins);
    if (edges.length < 2) throw new RangeError('histogram: bins needs at least two edges');
    for (let i = 1; i < edges.length; i++) {
      if (!(edges[i] > edges[i - 1])) throw new RangeError('histogram: bins edges must be strictly ascending');
    }
  }
  const nb = edges.length - 1;
  const counts = new Array(nb).fill(0);
  const lo = edges[0];
  const hi = edges[nb];
  const uniform = typeof bins === 'number';
  let total = 0;
  for (let i = 0; i < arr.length; i++) {
    const v = arr[i];
    if (!(v >= lo && v <= hi)) continue;
    let k;
    if (uniform) {
      k = Math.min(nb - 1, Math.floor(((v - lo) / (hi - lo)) * nb));
    } else {
      let a = 0;
      let b = nb;
      while (b - a > 1) {
        const m = (a + b) >> 1;
        if (v >= edges[m]) a = m;
        else b = m;
      }
      k = a;
    }
    counts[k]++;
    total++;
  }
  const centres = new Array(nb);
  const density = new Array(nb);
  for (let k = 0; k < nb; k++) {
    centres[k] = 0.5 * (edges[k] + edges[k + 1]);
    density[k] = total > 0 ? counts[k] / (total * (edges[k + 1] - edges[k])) : 0;
  }
  return { edges, counts, centres, density };
}

/**
 * Unbiased sample covariance (divisor n-1). Validity: equal lengths, n >= 2.
 * @param {ArrayLike<number>} a first variable [ua]
 * @param {ArrayLike<number>} b second variable [ub]
 * @returns {number} covariance [ua·ub]
 */
export function covariance(a, b) {
  need(a, 'covariance', 2);
  if (!b || b.length !== a.length) throw new RangeError('covariance: a and b must have the same length');
  const ma = mean(a);
  const mb = mean(b);
  let s = 0;
  for (let i = 0; i < a.length; i++) s += (a[i] - ma) * (b[i] - mb);
  return s / (a.length - 1);
}

/**
 * Pearson product-moment correlation coefficient. Validity: equal lengths, n >= 2; returns 0 when either variable
 * is constant.
 * @param {ArrayLike<number>} a first variable [ua]
 * @param {ArrayLike<number>} b second variable [ub]
 * @returns {number} correlation in [-1, 1] [-]
 */
export function correlation(a, b) {
  need(a, 'correlation', 2);
  if (!b || b.length !== a.length) throw new RangeError('correlation: a and b must have the same length');
  const ma = mean(a);
  const mb = mean(b);
  let sab = 0;
  let saa = 0;
  let sbb = 0;
  for (let i = 0; i < a.length; i++) {
    const da = a[i] - ma;
    const db = b[i] - mb;
    sab += da * db;
    saa += da * da;
    sbb += db * db;
  }
  const den = Math.sqrt(saa * sbb);
  return den > 0 ? sab / den : 0;
}

/**
 * Fractional ranks (1-based, ties receive the average rank).
 * @param {ArrayLike<number>} arr data [u]
 * @returns {Float64Array} ranks 1..n [-]
 */
export function ranks(arr) {
  need(arr, 'ranks');
  const n = arr.length;
  const idx = new Uint32Array(n);
  for (let i = 0; i < n; i++) idx[i] = i;
  idx.sort((i, j) => arr[i] - arr[j]);
  const r = new Float64Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    while (j + 1 < n && arr[idx[j + 1]] === arr[idx[i]]) j++;
    const avg = 0.5 * (i + j) + 1;
    for (let k = i; k <= j; k++) r[idx[k]] = avg;
    i = j + 1;
  }
  return r;
}

/**
 * Spearman rank correlation coefficient (Pearson correlation of average ranks).
 * @param {ArrayLike<number>} a first variable [ua]
 * @param {ArrayLike<number>} b second variable [ub]
 * @returns {number} rank correlation in [-1, 1] [-]
 */
export function rankCorrelation(a, b) {
  need(a, 'rankCorrelation', 2);
  if (!b || b.length !== a.length) throw new RangeError('rankCorrelation: a and b must have the same length');
  return correlation(ranks(a), ranks(b));
}

// ---------------------------------------------------------------------------------------------------------------
// Standard normal helpers
// ---------------------------------------------------------------------------------------------------------------

/**
 * Standard normal probability density φ(z).
 * @param {number} z standard normal variate [-]
 * @returns {number} density [-]
 */
export function normPdf(z) {
  return Math.exp(-0.5 * z * z) / SQRT2PI;
}

/**
 * Standard normal cumulative distribution Φ(z) = erfc(-z/√2)/2 (accurate in both tails).
 * @param {number} z standard normal variate [-]
 * @returns {number} probability [-]
 */
export function normCdf(z) {
  return 0.5 * erfc(-z / SQRT2);
}

const NA = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239];
const NB = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
const NC = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
const ND = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];

/**
 * Inverse standard normal distribution Φ⁻¹(p).
 * Algorithm: P. J. Acklam's rational approximation (relative error 1.15e-9) followed by Halley refinement steps
 * on Φ evaluated with erfc, giving close to machine precision including the far tails.
 * Validity: 0 <= p <= 1 (returns ∓Infinity at the ends).
 * @param {number} p probability [-]
 * @returns {number} standard normal quantile [-]
 */
export function normInv(p) {
  if (!(p >= 0 && p <= 1)) throw new RangeError('normInv: p must be within 0..1');
  if (p === 0) return -Infinity;
  if (p === 1) return Infinity;
  const pLow = 0.02425;
  let x;
  if (p < pLow) {
    const q = Math.sqrt(-2 * Math.log(p));
    x = (((((NC[0] * q + NC[1]) * q + NC[2]) * q + NC[3]) * q + NC[4]) * q + NC[5]) / ((((ND[0] * q + ND[1]) * q + ND[2]) * q + ND[3]) * q + 1);
  } else if (p <= 1 - pLow) {
    const q = p - 0.5;
    const r = q * q;
    x = ((((((NA[0] * r + NA[1]) * r + NA[2]) * r + NA[3]) * r + NA[4]) * r + NA[5]) * q) / (((((NB[0] * r + NB[1]) * r + NB[2]) * r + NB[3]) * r + NB[4]) * r + 1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((NC[0] * q + NC[1]) * q + NC[2]) * q + NC[3]) * q + NC[4]) * q + NC[5]) / ((((ND[0] * q + ND[1]) * q + ND[2]) * q + ND[3]) * q + 1);
  }
  for (let it = 0; it < 2; it++) {
    // work on the smaller tail to keep relative precision
    const e = x < 0 ? 0.5 * erfc(-x / SQRT2) - p : -(0.5 * erfc(x / SQRT2) - (1 - p));
    const u = e * SQRT2PI * Math.exp(0.5 * x * x);
    if (!Number.isFinite(u)) break;
    x -= u / (1 + 0.5 * x * u);
  }
  return x;
}

// ---------------------------------------------------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------------------------------------------------

function gammaFn(x) {
  return Math.exp(lgamma(x));
}

function num(spec, key, fallback) {
  const v = spec[key];
  if (v === undefined || v === null) {
    if (fallback !== undefined) return fallback;
    throw new RangeError(`dist(${spec.dist}): ${key} is required`);
  }
  if (typeof v !== 'number' || Number.isNaN(v)) throw new RangeError(`dist(${spec.dist}): ${key} must be a number`);
  return v;
}

function has(spec, key) {
  return spec[key] !== undefined && spec[key] !== null;
}

/** Standard deviation requested through { std } or { mean, cov }. Returns undefined when neither is present. */
function stdOf(spec) {
  if (has(spec, 'std')) return num(spec, 'std');
  if (has(spec, 'cov')) return Math.abs(num(spec, 'mean')) * num(spec, 'cov');
  return undefined;
}

function checkP(p, name) {
  if (!(p >= 0 && p <= 1)) throw new RangeError(`dist(${name}).icdf: p must be within 0..1`);
}

/** Safeguarded Newton inversion of a continuous cdf on a finite bracket [lo, hi]. */
function invertCdf(cdf, pdf, p, lo, hi, x0) {
  let a = lo;
  let b = hi;
  let x = x0 > a && x0 < b ? x0 : 0.5 * (a + b);
  for (let it = 0; it < 200; it++) {
    const f = cdf(x) - p;
    if (f === 0) return x;
    if (f > 0) b = x;
    else a = x;
    const d = pdf(x);
    let xn = d > 0 && Number.isFinite(d) ? x - f / d : NaN;
    if (!(xn > a && xn < b)) xn = 0.5 * (a + b);
    if (Math.abs(xn - x) <= 4e-16 * Math.max(Math.abs(x), 1e-300) || b - a <= 4e-16 * Math.max(Math.abs(a), Math.abs(b), 1e-300)) return xn;
    x = xn;
  }
  return x;
}

/** Gamma(k, 1) deviate, Marsaglia & Tsang (2000) "A simple method for generating gamma variables". */
function sampleGammaUnit(k, rng) {
  if (k < 1) {
    const u = 1 - rng.next();
    return sampleGammaUnit(k + 1, rng) * Math.pow(u, 1 / k);
  }
  const d = k - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x;
    let v;
    do {
      x = rng.normal();
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;
    const u = 1 - rng.next();
    if (u < 1 - 0.0331 * x * x * x * x) return d * v;
    if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
  }
}

function openUniform(rng) {
  let u = rng.next();
  while (u <= 0) u = rng.next();
  return u;
}

function finish(name, spec, o) {
  const d = {
    name,
    spec,
    pdf: o.pdf,
    cdf: o.cdf,
    icdf: o.icdf,
    mean: o.mean,
    std: o.std,
    variance: o.std * o.std,
    support: o.support,
    params: o.params || {},
    sample: o.sample || ((rng) => o.icdf(openUniform(rng))),
    /** x for a standard normal variate z (isoprobabilistic transform x = F⁻¹(Φ(z))). */
    fromZ:
      o.fromZ ||
      ((z) => {
        // use the survival side for the upper tail to keep precision
        const p = normCdf(z);
        return o.icdf(Math.min(Math.max(p, 1e-300), 1 - 1.1e-16));
      }),
    /** Standard normal variate for x (z = Φ⁻¹(F(x))). */
    toZ:
      o.toZ ||
      ((x) => {
        const p = o.cdf(x);
        return normInv(Math.min(Math.max(p, 1e-300), 1 - 1.1e-16));
      }),
  };
  return d;
}

function makeNormal(spec) {
  const m = num(spec, 'mean');
  const s = stdOf(spec);
  if (!(s > 0)) throw new RangeError('dist(normal): std must be > 0');
  return finish('normal', spec, {
    pdf: (x) => normPdf((x - m) / s) / s,
    cdf: (x) => normCdf((x - m) / s),
    icdf: (p) => {
      checkP(p, 'normal');
      return m + s * normInv(p);
    },
    mean: m,
    std: s,
    support: [-Infinity, Infinity],
    params: { mean: m, std: s },
    sample: (rng) => m + s * rng.normal(),
    fromZ: (z) => m + s * z,
    toZ: (x) => (x - m) / s,
  });
}

function makeLognormal(spec) {
  let mu;
  let sig;
  if (has(spec, 'mu') && has(spec, 'sigma') && !has(spec, 'mean')) {
    mu = num(spec, 'mu');
    sig = num(spec, 'sigma');
  } else {
    const m = num(spec, 'mean');
    const s = stdOf(spec);
    if (!(m > 0)) throw new RangeError('dist(lognormal): mean must be > 0');
    if (!(s > 0)) throw new RangeError('dist(lognormal): std must be > 0');
    const v2 = (s / m) * (s / m);
    sig = Math.sqrt(Math.log1p(v2));
    mu = Math.log(m) - 0.5 * sig * sig;
  }
  if (!(sig > 0)) throw new RangeError('dist(lognormal): sigma must be > 0');
  const m = Math.exp(mu + 0.5 * sig * sig);
  const s = m * Math.sqrt(Math.expm1(sig * sig));
  return finish('lognormal', spec, {
    pdf: (x) => (x > 0 ? normPdf((Math.log(x) - mu) / sig) / (x * sig) : 0),
    cdf: (x) => (x > 0 ? normCdf((Math.log(x) - mu) / sig) : 0),
    icdf: (p) => {
      checkP(p, 'lognormal');
      return p === 0 ? 0 : Math.exp(mu + sig * normInv(p));
    },
    mean: m,
    std: s,
    support: [0, Infinity],
    params: { mu, sigma: sig },
    sample: (rng) => Math.exp(mu + sig * rng.normal()),
    fromZ: (z) => Math.exp(mu + sig * z),
    toZ: (x) => (Math.log(x) - mu) / sig,
  });
}

function makeUniform(spec) {
  let a;
  let b;
  if (has(spec, 'min') && has(spec, 'max')) {
    a = num(spec, 'min');
    b = num(spec, 'max');
  } else {
    const m = num(spec, 'mean');
    const s = stdOf(spec);
    if (!(s > 0)) throw new RangeError('dist(uniform): min and max (or mean and std/cov) are required');
    a = m - Math.sqrt(3) * s;
    b = m + Math.sqrt(3) * s;
  }
  if (!(b > a)) throw new RangeError('dist(uniform): max must be > min');
  const w = b - a;
  return finish('uniform', spec, {
    pdf: (x) => (x >= a && x <= b ? 1 / w : 0),
    cdf: (x) => (x <= a ? 0 : x >= b ? 1 : (x - a) / w),
    icdf: (p) => {
      checkP(p, 'uniform');
      return a + w * p;
    },
    mean: 0.5 * (a + b),
    std: w / Math.sqrt(12),
    support: [a, b],
    params: { min: a, max: b },
  });
}

function makeTriangular(spec) {
  let a;
  let c;
  let b;
  if (has(spec, 'min') && has(spec, 'max')) {
    a = num(spec, 'min');
    b = num(spec, 'max');
    c = num(spec, 'mode', 0.5 * (a + b));
  } else {
    const m = num(spec, 'mean');
    const s = stdOf(spec);
    if (!(s > 0)) throw new RangeError('dist(triangular): min, mode and max (or mean and std/cov) are required');
    a = m - Math.sqrt(6) * s;
    b = m + Math.sqrt(6) * s;
    c = m;
  }
  if (!(b > a)) throw new RangeError('dist(triangular): max must be > min');
  if (!(c >= a && c <= b)) throw new RangeError('dist(triangular): mode must lie within [min, max]');
  const w = b - a;
  const fc = (c - a) / w;
  return finish('triangular', spec, {
    pdf: (x) => {
      if (x < a || x > b) return 0;
      if (x < c) return (2 * (x - a)) / (w * (c - a));
      if (x > c) return (2 * (b - x)) / (w * (b - c));
      return 2 / w;
    },
    cdf: (x) => {
      if (x <= a) return 0;
      if (x >= b) return 1;
      if (x <= c) return ((x - a) * (x - a)) / (w * (c - a));
      return 1 - ((b - x) * (b - x)) / (w * (b - c));
    },
    icdf: (p) => {
      checkP(p, 'triangular');
      if (p < fc) return a + Math.sqrt(p * w * (c - a));
      return b - Math.sqrt((1 - p) * w * (b - c));
    },
    mean: (a + b + c) / 3,
    std: Math.sqrt((a * a + b * b + c * c - a * b - a * c - b * c) / 18),
    support: [a, b],
    params: { min: a, mode: c, max: b },
  });
}

function betaCore(name, spec, al, be, lo, hi) {
  if (!(al > 0)) throw new RangeError(`dist(${name}): shape a must be > 0`);
  if (!(be > 0)) throw new RangeError(`dist(${name}): shape b must be > 0`);
  if (!(hi > lo)) throw new RangeError(`dist(${name}): max must be > min`);
  const w = hi - lo;
  const lnB = lgamma(al) + lgamma(be) - lgamma(al + be);
  const pdf01 = (t) => {
    if (t < 0 || t > 1) return 0;
    if (t === 0) return al < 1 ? Infinity : al === 1 ? be : 0;
    if (t === 1) return be < 1 ? Infinity : be === 1 ? al : 0;
    return Math.exp((al - 1) * Math.log(t) + (be - 1) * Math.log1p(-t) - lnB);
  };
  const cdf01 = (t) => (t <= 0 ? 0 : t >= 1 ? 1 : betaI(al, be, t));
  const m01 = al / (al + be);
  const s01 = Math.sqrt((al * be) / ((al + be) * (al + be) * (al + be + 1)));
  return finish(name, spec, {
    pdf: (x) => pdf01((x - lo) / w) / w,
    cdf: (x) => cdf01((x - lo) / w),
    icdf: (p) => {
      checkP(p, name);
      if (p === 0) return lo;
      if (p === 1) return hi;
      return lo + w * invertCdf(cdf01, pdf01, p, 0, 1, m01);
    },
    mean: lo + w * m01,
    std: w * s01,
    support: [lo, hi],
    params: { a: al, b: be, min: lo, max: hi },
    sample: (rng) => {
      const g1 = sampleGammaUnit(al, rng);
      const g2 = sampleGammaUnit(be, rng);
      return lo + (w * g1) / (g1 + g2);
    },
  });
}

function makeBeta(spec) {
  const lo = num(spec, 'min', 0);
  const hi = num(spec, 'max', 1);
  let al;
  let be;
  if (has(spec, 'a') && has(spec, 'b')) {
    al = num(spec, 'a');
    be = num(spec, 'b');
  } else {
    // method of moments on [min, max]
    const m = (num(spec, 'mean') - lo) / (hi - lo);
    const s = stdOf(spec) / (hi - lo);
    const t = (m * (1 - m)) / (s * s) - 1;
    if (!(m > 0 && m < 1 && t > 0)) throw new RangeError('dist(beta): mean/std are not attainable on [min, max]');
    al = m * t;
    be = (1 - m) * t;
  }
  return betaCore('beta', spec, al, be, lo, hi);
}

function makePert(spec) {
  const a = num(spec, 'min');
  const b = num(spec, 'max');
  const c = num(spec, 'mode');
  const lam = num(spec, 'lambda', 4);
  if (!(b > a)) throw new RangeError('dist(pert): max must be > min');
  if (!(c >= a && c <= b)) throw new RangeError('dist(pert): mode must lie within [min, max]');
  if (!(lam > 0)) throw new RangeError('dist(pert): lambda must be > 0');
  const al = 1 + (lam * (c - a)) / (b - a);
  const be = 1 + (lam * (b - c)) / (b - a);
  const d = betaCore('pert', spec, al, be, a, b);
  d.params = { min: a, mode: c, max: b, lambda: lam, a: al, b: be };
  return d;
}

function weibullShapeFromCov(cov) {
  // cov² = Γ(1+2/k)/Γ(1+1/k)² - 1, monotonically decreasing in k
  const f = (k) => Math.exp(lgamma(1 + 2 / k) - 2 * lgamma(1 + 1 / k)) - 1 - cov * cov;
  let a = 0.05;
  let b = 500;
  if (f(a) < 0 || f(b) > 0) throw new RangeError('dist(weibull): cov is outside the attainable range');
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (a + b);
    if (f(m) > 0) a = m;
    else b = m;
    if (b - a < 1e-14 * b) break;
  }
  return 0.5 * (a + b);
}

function makeWeibull(spec) {
  const loc = num(spec, 'loc', 0);
  let k;
  let lam;
  if (has(spec, 'scale') && has(spec, 'shape')) {
    k = num(spec, 'shape');
    lam = num(spec, 'scale');
  } else {
    const m = num(spec, 'mean') - loc;
    const s = stdOf(spec);
    if (!(m > 0 && s > 0)) throw new RangeError('dist(weibull): scale and shape (or mean > loc and std/cov) are required');
    k = weibullShapeFromCov(s / m);
    lam = m / gammaFn(1 + 1 / k);
  }
  if (!(k > 0)) throw new RangeError('dist(weibull): shape must be > 0');
  if (!(lam > 0)) throw new RangeError('dist(weibull): scale must be > 0');
  const g1 = gammaFn(1 + 1 / k);
  const g2 = gammaFn(1 + 2 / k);
  return finish('weibull', spec, {
    pdf: (x) => {
      const z = (x - loc) / lam;
      if (z < 0) return 0;
      if (z === 0) return k < 1 ? Infinity : k === 1 ? 1 / lam : 0;
      return (k / lam) * Math.pow(z, k - 1) * Math.exp(-Math.pow(z, k));
    },
    cdf: (x) => {
      const z = (x - loc) / lam;
      return z <= 0 ? 0 : -Math.expm1(-Math.pow(z, k));
    },
    icdf: (p) => {
      checkP(p, 'weibull');
      return loc + lam * Math.pow(-Math.log1p(-p), 1 / k);
    },
    mean: loc + lam * g1,
    std: lam * Math.sqrt(Math.max(g2 - g1 * g1, 0)),
    support: [loc, Infinity],
    params: { scale: lam, shape: k, loc },
    fromZ: (z) => {
      // -ln(1 - Φ(z)) evaluated on the survival side
      const q = normCdf(-z);
      return loc + lam * Math.pow(-Math.log(Math.max(q, 1e-300)), 1 / k);
    },
  });
}

function makeGumbel(spec) {
  let loc;
  let sc;
  if (has(spec, 'loc') && has(spec, 'scale')) {
    loc = num(spec, 'loc');
    sc = num(spec, 'scale');
  } else {
    const m = num(spec, 'mean');
    const s = stdOf(spec);
    if (!(s > 0)) throw new RangeError('dist(gumbel): loc and scale (or mean and std/cov) are required');
    sc = (s * Math.sqrt(6)) / Math.PI;
    loc = m - EULER_GAMMA * sc;
  }
  if (!(sc > 0)) throw new RangeError('dist(gumbel): scale must be > 0');
  return finish('gumbel', spec, {
    pdf: (x) => {
      const z = (x - loc) / sc;
      return Math.exp(-z - Math.exp(-z)) / sc;
    },
    cdf: (x) => Math.exp(-Math.exp(-(x - loc) / sc)),
    icdf: (p) => {
      checkP(p, 'gumbel');
      return loc - sc * Math.log(-Math.log(p));
    },
    mean: loc + EULER_GAMMA * sc,
    std: (Math.PI * sc) / Math.sqrt(6),
    support: [-Infinity, Infinity],
    params: { loc, scale: sc },
    fromZ: (z) => {
      // -ln Φ(z): for large z use -ln(1-q) = -log1p(-q) with q the survival probability
      const q = normCdf(-z);
      const nl = q < 0.5 ? -Math.log1p(-q) : -Math.log(Math.max(normCdf(z), 1e-300));
      return loc - sc * Math.log(Math.max(nl, 1e-300));
    },
  });
}

function makeGev(spec) {
  const loc = num(spec, 'loc');
  const sc = num(spec, 'scale');
  const xi = num(spec, 'shape', 0);
  if (!(sc > 0)) throw new RangeError('dist(gev): scale must be > 0');
  if (Math.abs(xi) < 1e-12) {
    const g = makeGumbel({ dist: 'gumbel', loc, scale: sc });
    g.name = 'gev';
    g.spec = spec;
    g.params = { loc, scale: sc, shape: 0 };
    return g;
  }
  let m = Infinity;
  let s = Infinity;
  if (xi < 1) {
    const g1 = gammaFn(1 - xi);
    m = loc + (sc * (g1 - 1)) / xi;
    if (xi < 0.5) s = (sc * Math.sqrt(Math.max(gammaFn(1 - 2 * xi) - g1 * g1, 0))) / Math.abs(xi);
  }
  const bound = loc - sc / xi;
  return finish('gev', spec, {
    pdf: (x) => {
      const t = 1 + (xi * (x - loc)) / sc;
      if (t <= 0) return 0;
      const lt = Math.log(t);
      return Math.exp(-(1 / xi + 1) * lt - Math.exp(-lt / xi)) / sc;
    },
    cdf: (x) => {
      const t = 1 + (xi * (x - loc)) / sc;
      if (t <= 0) return xi > 0 ? 0 : 1;
      return Math.exp(-Math.pow(t, -1 / xi));
    },
    icdf: (p) => {
      checkP(p, 'gev');
      if (p === 0) return xi > 0 ? bound : -Infinity;
      if (p === 1) return xi > 0 ? Infinity : bound;
      return loc + (sc * Math.expm1(-xi * Math.log(-Math.log(p)))) / xi;
    },
    mean: m,
    std: s,
    support: xi > 0 ? [bound, Infinity] : [-Infinity, bound],
    params: { loc, scale: sc, shape: xi },
    fromZ: (z) => {
      // -ln Φ(z) from the survival probability for the upper tail
      const q = normCdf(-z);
      const nl = q < 0.5 ? -Math.log1p(-q) : -Math.log(Math.max(normCdf(z), 1e-300));
      return loc + (sc * Math.expm1(-xi * Math.log(Math.max(nl, 1e-300)))) / xi;
    },
  });
}

function makeGpd(spec) {
  const loc = num(spec, 'loc', 0);
  const sc = num(spec, 'scale');
  const xi = num(spec, 'shape', 0);
  if (!(sc > 0)) throw new RangeError('dist(gpd): scale must be > 0');
  const zero = Math.abs(xi) < 1e-12;
  const upper = xi < 0 && !zero ? loc - sc / xi : Infinity;
  return finish('gpd', spec, {
    pdf: (x) => {
      const z = (x - loc) / sc;
      if (z < 0) return 0;
      if (zero) return Math.exp(-z) / sc;
      const t = 1 + xi * z;
      if (t <= 0) return 0;
      return Math.pow(t, -1 / xi - 1) / sc;
    },
    cdf: (x) => {
      const z = (x - loc) / sc;
      if (z <= 0) return 0;
      if (zero) return -Math.expm1(-z);
      const t = 1 + xi * z;
      if (t <= 0) return 1;
      return -Math.expm1(-Math.log(t) / xi);
    },
    icdf: (p) => {
      checkP(p, 'gpd');
      if (p === 1) return upper;
      const l = -Math.log1p(-p);
      return zero ? loc + sc * l : loc + (sc * Math.expm1(xi * l)) / xi;
    },
    mean: xi < 1 ? loc + sc / (1 - xi) : Infinity,
    std: xi < 0.5 ? sc / ((1 - xi) * Math.sqrt(1 - 2 * xi)) : Infinity,
    support: [loc, upper],
    params: { loc, scale: sc, shape: xi },
    fromZ: (z) => {
      const l = -Math.log(Math.max(normCdf(-z), 1e-300)); // -ln(1 - Φ(z))
      return zero ? loc + sc * l : loc + (sc * Math.expm1(xi * l)) / xi;
    },
  });
}

function makeExponential(spec) {
  const rate = has(spec, 'rate') ? num(spec, 'rate') : 1 / num(spec, 'mean');
  if (!(rate > 0) || !Number.isFinite(rate)) throw new RangeError('dist(exponential): rate must be > 0');
  return finish('exponential', spec, {
    pdf: (x) => (x < 0 ? 0 : rate * Math.exp(-rate * x)),
    cdf: (x) => (x <= 0 ? 0 : -Math.expm1(-rate * x)),
    icdf: (p) => {
      checkP(p, 'exponential');
      return -Math.log1p(-p) / rate;
    },
    mean: 1 / rate,
    std: 1 / rate,
    support: [0, Infinity],
    params: { rate },
    fromZ: (z) => -Math.log(Math.max(normCdf(-z), 1e-300)) / rate,
  });
}

function makeGamma(spec) {
  let k;
  let th;
  if (has(spec, 'shape') && has(spec, 'scale')) {
    k = num(spec, 'shape');
    th = num(spec, 'scale');
  } else {
    const m = num(spec, 'mean');
    const s = stdOf(spec);
    if (!(m > 0 && s > 0)) throw new RangeError('dist(gamma): shape and scale (or mean > 0 and std/cov) are required');
    k = (m * m) / (s * s);
    th = (s * s) / m;
  }
  if (!(k > 0)) throw new RangeError('dist(gamma): shape must be > 0');
  if (!(th > 0)) throw new RangeError('dist(gamma): scale must be > 0');
  const lg = lgamma(k);
  const pdf1 = (t) => {
    if (t < 0) return 0;
    if (t === 0) return k < 1 ? Infinity : k === 1 ? 1 : 0;
    return Math.exp((k - 1) * Math.log(t) - t - lg);
  };
  const cdf1 = (t) => (t <= 0 ? 0 : gammaP(k, t));
  return finish('gamma', spec, {
    pdf: (x) => pdf1(x / th) / th,
    cdf: (x) => cdf1(x / th),
    icdf: (p) => {
      checkP(p, 'gamma');
      if (p === 0) return 0;
      if (p === 1) return Infinity;
      // Wilson–Hilferty starting value, then safeguarded Newton
      const z = normInv(p);
      const c = 1 / (9 * k);
      let x0 = k * Math.pow(Math.max(1 - c + z * Math.sqrt(c), 1e-3), 3);
      let hi = Math.max(2 * x0, k + 10 * Math.sqrt(k) + 20);
      while (cdf1(hi) < p && hi < 1e300) hi *= 2;
      if (!(x0 > 0 && x0 < hi)) x0 = 0.5 * hi;
      return th * invertCdf(cdf1, pdf1, p, 0, hi, x0);
    },
    mean: k * th,
    std: Math.sqrt(k) * th,
    support: [0, Infinity],
    params: { shape: k, scale: th },
    sample: (rng) => th * sampleGammaUnit(k, rng),
    toZ: (x) => {
      const t = x / th;
      const p = cdf1(t);
      return p < 0.5 ? normInv(Math.max(p, 1e-300)) : -normInv(Math.max(gammaQ(k, t), 1e-300));
    },
  });
}

function makeRayleigh(spec) {
  const sg = has(spec, 'sigma') ? num(spec, 'sigma') : num(spec, 'mean') / Math.sqrt(Math.PI / 2);
  if (!(sg > 0)) throw new RangeError('dist(rayleigh): sigma must be > 0');
  return finish('rayleigh', spec, {
    pdf: (x) => (x < 0 ? 0 : (x / (sg * sg)) * Math.exp((-x * x) / (2 * sg * sg))),
    cdf: (x) => (x <= 0 ? 0 : -Math.expm1((-x * x) / (2 * sg * sg))),
    icdf: (p) => {
      checkP(p, 'rayleigh');
      return sg * Math.sqrt(-2 * Math.log1p(-p));
    },
    mean: sg * Math.sqrt(Math.PI / 2),
    std: sg * Math.sqrt(2 - Math.PI / 2),
    support: [0, Infinity],
    params: { sigma: sg },
    fromZ: (z) => sg * Math.sqrt(-2 * Math.log(Math.max(normCdf(-z), 1e-300))),
  });
}

function makeDeterministic(spec) {
  const v = has(spec, 'value') ? num(spec, 'value') : num(spec, 'mean');
  return finish('deterministic', spec, {
    pdf: (x) => (x === v ? Infinity : 0),
    cdf: (x) => (x >= v ? 1 : 0),
    icdf: (p) => {
      checkP(p, 'deterministic');
      return v;
    },
    mean: v,
    std: 0,
    support: [v, v],
    params: { value: v },
    sample: () => v,
    fromZ: () => v,
    toZ: () => 0,
  });
}

const MAKERS = {
  normal: makeNormal,
  gaussian: makeNormal,
  lognormal: makeLognormal,
  uniform: makeUniform,
  triangular: makeTriangular,
  pert: makePert,
  beta: makeBeta,
  weibull: makeWeibull,
  gumbel: makeGumbel,
  gev: makeGev,
  gpd: makeGpd,
  exponential: makeExponential,
  gamma: makeGamma,
  rayleigh: makeRayleigh,
  deterministic: makeDeterministic,
  constant: makeDeterministic,
};

/** Names of the distributions understood by {@link dist}. */
export const DISTRIBUTIONS = ['normal', 'lognormal', 'uniform', 'triangular', 'pert', 'beta', 'weibull', 'gumbel', 'gev', 'gpd', 'exponential', 'gamma', 'rayleigh', 'deterministic'];

/**
 * Build a probability distribution object from a declarative specification.
 *
 * Supported `spec.dist` and parameters (x carries the units of the variable, [u]):
 *  - `normal` { mean, std }
 *  - `lognormal` { mean, std } — mean and standard deviation of X itself (not of ln X); alternatively { mu, sigma }
 *    of ln X
 *  - `uniform` { min, max }
 *  - `triangular` { min, mode, max }
 *  - `pert` { min, mode, max, lambda = 4 } — beta-PERT: Beta(1 + λ(mode-min)/(max-min), 1 + λ(max-mode)/(max-min))
 *    scaled to [min, max] (Malcolm et al. 1959; Vose 2008 "Risk Analysis" §III)
 *  - `beta` { a, b, min = 0, max = 1 }
 *  - `weibull` { scale, shape, loc = 0 } — F = 1 - exp(-((x-loc)/scale)^shape)
 *  - `gumbel` { loc, scale } — maxima, F = exp(-exp(-(x-loc)/scale))
 *  - `gev` { loc, scale, shape ξ } — F = exp(-(1 + ξ(x-loc)/scale)^(-1/ξ)); ξ > 0 heavy (Fréchet) tail, ξ < 0 bounded
 *    (Weibull-type) upper tail, ξ = 0 Gumbel (Coles 2001 convention; Hosking's k = -ξ)
 *  - `gpd` { loc, scale, shape ξ } — F = 1 - (1 + ξ(x-loc)/scale)^(-1/ξ), same sign convention as `gev`
 *  - `exponential` { rate } [1/u]
 *  - `gamma` { shape, scale }
 *  - `rayleigh` { sigma }
 *  - `deterministic` { value }
 *
 * Convenience: `{ mean, cov }` may replace `{ mean, std }` (std = |mean|·cov), and the two-parameter families
 * `uniform`, `triangular` (symmetric), `beta`, `weibull`, `gumbel`, `gamma` accept { mean, std | cov } and are then
 * moment-matched; `exponential`, `rayleigh` and `deterministic` accept { mean }.
 *
 * icdf is exact: closed form everywhere except `beta`/`pert`/`gamma` (safeguarded Newton on the regularised
 * incomplete beta/gamma functions) and the normal family (Acklam approximation + Halley refinement).
 *
 * References: Johnson, Kotz & Balakrishnan, "Continuous Univariate Distributions" (1994/95); Coles, "An
 * Introduction to Statistical Modeling of Extreme Values" (2001).
 * Validity: parameters within their natural ranges, otherwise RangeError naming the parameter. `mean`/`std` are
 * Infinity where the moment does not exist (GEV/GPD with ξ >= 1 or 1/2).
 *
 * @param {object} spec distribution specification, see above
 * @returns {{ pdf: (x: number) => number, cdf: (x: number) => number, icdf: (p: number) => number, mean: number,
 *   std: number, sample: (rng: { next: () => number, normal: () => number }) => number, variance: number,
 *   support: [number, number], params: object, name: string, spec: object, fromZ: (z: number) => number,
 *   toZ: (x: number) => number }}
 *   pdf [1/u], cdf [-], icdf [u], mean [u], std [u], sample(rng) [u]; extras: variance [u²], support bounds [u],
 *   resolved native parameters, isoprobabilistic transforms fromZ (x = F⁻¹(Φ(z))) and toZ (z = Φ⁻¹(F(x))).
 */
export function dist(spec) {
  if (typeof spec === 'number') return makeDeterministic({ dist: 'deterministic', value: spec });
  if (!spec || typeof spec !== 'object') throw new RangeError('dist: spec must be an object');
  if (typeof spec.pdf === 'function' && typeof spec.icdf === 'function') return spec; // already a distribution
  const name = String(spec.dist === undefined ? '' : spec.dist).toLowerCase();
  const maker = MAKERS[name];
  if (!maker) throw new RangeError(`dist: unknown distribution "${spec.dist}"`);
  return maker(spec);
}

// ---------------------------------------------------------------------------------------------------------------
// Regression
// ---------------------------------------------------------------------------------------------------------------

/**
 * Ordinary least-squares straight line y = intercept + slope·x with standard errors.
 * Reference: Draper & Smith, "Applied Regression Analysis", 3rd ed., ch. 1. Validity: n >= 2, x not constant
 * (standard errors need n >= 3).
 * @param {ArrayLike<number>} x abscissae [ux]
 * @param {ArrayLike<number>} y ordinates [uy]
 * @returns {{ slope: number, intercept: number, r2: number, seSlope: number, seIntercept: number, n: number,
 *   sigma: number }} slope [uy/ux], intercept [uy], coefficient of determination [-], standard errors of slope
 *   [uy/ux] and intercept [uy], sample size, residual standard deviation [uy]
 */
export function linearRegression(x, y) {
  need(x, 'linearRegression: x', 2);
  if (!y || y.length !== x.length) throw new RangeError('linearRegression: x and y must have the same length');
  const n = x.length;
  const mx = mean(x);
  const my = mean(y);
  let sxx = 0;
  let sxy = 0;
  let syy = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx;
    const dy = y[i] - my;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  if (!(sxx > 0)) throw new RangeError('linearRegression: x must not be constant');
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  const sse = Math.max(syy - slope * sxy, 0);
  const r2 = syy > 0 ? 1 - sse / syy : 1;
  const s2 = n > 2 ? sse / (n - 2) : 0;
  return {
    slope,
    intercept,
    r2,
    seSlope: Math.sqrt(s2 / sxx),
    seIntercept: Math.sqrt(s2 * (1 / n + (mx * mx) / sxx)),
    n,
    sigma: Math.sqrt(s2),
  };
}

/**
 * Least-squares polynomial fit (Householder QR on the Vandermonde matrix of the centred and scaled abscissae,
 * coefficients converted back to powers of x).
 * Validity: number of points > degree; degree <= ~12 for well-spread data.
 * @param {ArrayLike<number>} x abscissae [ux]
 * @param {ArrayLike<number>} y ordinates [uy]
 * @param {number} degree polynomial degree [-]
 * @returns {number[]} coefficients, highest degree first [uy/ux^k]
 */
export function polyfit(x, y, degree) {
  need(x, 'polyfit: x');
  if (!y || y.length !== x.length) throw new RangeError('polyfit: x and y must have the same length');
  if (!Number.isInteger(degree) || degree < 0) throw new RangeError('polyfit: degree must be a non-negative integer');
  const n = x.length;
  if (n < degree + 1) throw new RangeError('polyfit: needs at least degree + 1 points');
  const xm = mean(x);
  let xs = 0;
  for (let i = 0; i < n; i++) xs = Math.max(xs, Math.abs(x[i] - xm));
  if (!(xs > 0)) xs = 1;
  const A = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = new Array(degree + 1);
    const t = (x[i] - xm) / xs;
    let p = 1;
    for (let k = 0; k <= degree; k++) {
      row[k] = p;
      p *= t;
    }
    A[i] = row;
  }
  const sol = leastSquares(A, Array.from(y)).x; // ascending powers of t
  // expand Σ a_k ((x - xm)/xs)^k into powers of x
  const asc = new Array(degree + 1).fill(0);
  let basis = [1]; // ((x - xm)/xs)^k, ascending powers of x
  for (let k = 0; k <= degree; k++) {
    for (let j = 0; j < basis.length; j++) asc[j] += sol[k] * basis[j];
    const nb = new Array(basis.length + 1).fill(0);
    for (let j = 0; j < basis.length; j++) {
      nb[j] -= (basis[j] * xm) / xs;
      nb[j + 1] += basis[j] / xs;
    }
    basis = nb;
  }
  return asc.reverse();
}

/**
 * Evaluate a polynomial by Horner's rule.
 * @param {ArrayLike<number>} c coefficients, highest degree first
 * @param {number|ArrayLike<number>} x evaluation point(s) [ux]
 * @returns {number|number[]} value(s) [uy]
 */
export function polyval(c, x) {
  const one = (v) => {
    let s = 0;
    for (let i = 0; i < c.length; i++) s = s * v + c[i];
    return s;
  };
  if (typeof x === 'number') return one(x);
  return Array.from(x, one);
}

/**
 * Multiple linear regression y = β0 + Σ βj·xj by Householder QR, with coefficient standard errors
 * se = sqrt(diag((AᵀA)⁻¹)·SSE/(n - p)).
 * Reference: Draper & Smith, "Applied Regression Analysis", 3rd ed., ch. 5. Validity: n > number of regressors + 1
 * and linearly independent columns.
 * @param {ArrayLike<ArrayLike<number>>} X regressor rows (n × k), intercept column added automatically
 * @param {ArrayLike<number>} y response (n) [uy]
 * @returns {{ beta: number[], r2: number, residuals: number[], se: number[], sigma: number, adjR2: number }}
 *   coefficients (beta[0] = intercept), coefficient of determination [-], residuals y - ŷ [uy], standard errors
 *   of beta, residual standard deviation [uy], adjusted R² [-]
 */
export function multipleRegression(X, y) {
  need(X, 'multipleRegression: X');
  if (!y || y.length !== X.length) throw new RangeError('multipleRegression: X and y must have the same number of rows');
  const n = X.length;
  const k = X[0].length;
  const p = k + 1;
  if (n < p) throw new RangeError('multipleRegression: needs at least as many rows as coefficients');
  const A = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = new Array(p);
    row[0] = 1;
    for (let j = 0; j < k; j++) row[j + 1] = X[i][j];
    A[i] = row;
  }
  const beta = Array.from(leastSquares(A, Array.from(y)).x);
  const residuals = new Array(n);
  const my = mean(y);
  let sse = 0;
  let sst = 0;
  for (let i = 0; i < n; i++) {
    let f = 0;
    for (let j = 0; j < p; j++) f += A[i][j] * beta[j];
    residuals[i] = y[i] - f;
    sse += residuals[i] * residuals[i];
    sst += (y[i] - my) * (y[i] - my);
  }
  const AtA = new Array(p);
  for (let a = 0; a < p; a++) {
    AtA[a] = new Array(p).fill(0);
    for (let b = 0; b <= a; b++) {
      let s = 0;
      for (let i = 0; i < n; i++) s += A[i][a] * A[i][b];
      AtA[a][b] = s;
    }
  }
  for (let a = 0; a < p; a++) for (let b = a + 1; b < p; b++) AtA[a][b] = AtA[b][a];
  const s2 = n > p ? sse / (n - p) : 0;
  const inv = inverse(AtA);
  const se = new Array(p);
  for (let j = 0; j < p; j++) se[j] = Math.sqrt(Math.max(inv[j][j] * s2, 0));
  const r2 = sst > 0 ? 1 - sse / sst : 1;
  return { beta, r2, residuals, se, sigma: Math.sqrt(s2), adjR2: n > p ? 1 - ((1 - r2) * (n - 1)) / (n - p) : r2 };
}

/**
 * Agreement metrics between predictions and observations.
 * rmse = sqrt(mean((pred-obs)²)); mae = mean|pred-obs|; bias = mean(pred - obs); nrmse = rmse / (max(obs) -
 * min(obs)) (NaN when the observations are constant); r2 = 1 - SSE/SST (can be negative); maxAbs = max|pred-obs|.
 * @param {ArrayLike<number>} pred predictions [u]
 * @param {ArrayLike<number>} obs observations [u]
 * @returns {{ rmse: number, mae: number, bias: number, nrmse: number, r2: number, maxAbs: number, n: number }}
 *   rmse, mae, bias, maxAbs in [u]; nrmse, r2 dimensionless; n sample size
 */
export function metrics(pred, obs) {
  need(pred, 'metrics: pred');
  if (!obs || obs.length !== pred.length) throw new RangeError('metrics: pred and obs must have the same length');
  const n = pred.length;
  const mo = mean(obs);
  let sse = 0;
  let sae = 0;
  let sb = 0;
  let sst = 0;
  let mx = 0;
  for (let i = 0; i < n; i++) {
    const e = pred[i] - obs[i];
    sse += e * e;
    sae += Math.abs(e);
    sb += e;
    sst += (obs[i] - mo) * (obs[i] - mo);
    if (Math.abs(e) > mx) mx = Math.abs(e);
  }
  const rmse = Math.sqrt(sse / n);
  const range = max(obs) - min(obs);
  return {
    rmse,
    mae: sae / n,
    bias: sb / n,
    nrmse: range > 0 ? rmse / range : NaN,
    r2: sst > 0 ? 1 - sse / sst : sse === 0 ? 1 : NaN,
    maxAbs: mx,
    n,
  };
}

/**
 * One-sample Kolmogorov–Smirnov test against a fully specified continuous distribution.
 * D = sup|Fn - F|; p-value from the Kolmogorov limiting distribution with Stephens' finite-sample correction
 * λ = (√n + 0.12 + 0.11/√n)·D, Q(λ) = 2 Σ (-1)^(j-1) exp(-2 j² λ²).
 * Reference: Stephens (1970), J. R. Stat. Soc. B 32, 115–122; Press et al., "Numerical Recipes" §14.3.
 * Validity: n >= ~5; p is not valid when the distribution parameters were estimated from the same sample.
 * @param {ArrayLike<number>} sample observations [u]
 * @param {(x: number) => number} cdf hypothesised cumulative distribution function
 * @returns {{ D: number, p: number, n: number }} KS statistic [-], p-value [-], sample size
 */
export function ksTest(sample, cdf) {
  need(sample, 'ksTest: sample');
  if (typeof cdf !== 'function') throw new RangeError('ksTest: cdf must be a function');
  const s = sortedCopy(sample);
  const n = s.length;
  let D = 0;
  for (let i = 0; i < n; i++) {
    const f = cdf(s[i]);
    const d = Math.max(f - i / n, (i + 1) / n - f);
    if (d > D) D = d;
  }
  const sn = Math.sqrt(n);
  const lam = (sn + 0.12 + 0.11 / sn) * D;
  let p;
  if (lam < 0.2) p = 1;
  else {
    let sum = 0;
    let sign = 1;
    for (let j = 1; j <= 100; j++) {
      const term = sign * Math.exp(-2 * j * j * lam * lam);
      sum += term;
      if (Math.abs(term) < 1e-16) break;
      sign = -sign;
    }
    p = Math.min(Math.max(2 * sum, 0), 1);
  }
  return { D, p, n };
}
