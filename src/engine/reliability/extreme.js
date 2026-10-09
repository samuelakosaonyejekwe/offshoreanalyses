/**
 * Extreme-value statistics for metocean and load data: distribution fitting (Gumbel, Weibull, GEV, GPD,
 * lognormal, normal), annual-maxima and peaks-over-threshold extraction, return values with bootstrap
 * confidence bands, IFORM environmental contours and scatter diagrams.
 *
 * Every fit result carries a `dist` field and the parameter names of math/stats.js `dist()`, so it can be passed
 * straight to `returnValue` or `dist`. Shape convention for GEV and GPD: ξ > 0 heavy (Fréchet / Pareto) tail,
 * ξ < 0 bounded upper tail — the same convention as `dist({ dist: 'gev' | 'gpd' })`; Hosking's k equals -ξ.
 */
import { makeRng } from '../math/rng.js';
import { dist, normInv, percentileSorted, EULER_GAMMA } from '../math/stats.js';
import { lgamma } from '../math/special.js';
import { brent } from '../math/roots.js';
import { goldenSection, nelderMead, levenbergMarquardt } from '../math/optim.js';

function clean(x, name, nMin) {
  if (!x || typeof x.length !== 'number') throw new RangeError(`${name}: x must be an array`);
  const a = Float64Array.from(x);
  for (let i = 0; i < a.length; i++) if (!Number.isFinite(a[i])) throw new RangeError(`${name}: x[${i}] is not a finite number`);
  if (a.length < nMin) throw new RangeError(`${name}: needs at least ${nMin} values`);
  return a;
}

function gam(x) {
  return Math.exp(lgamma(x));
}

function logLikOf(spec, x) {
  const d = dist(spec);
  let s = 0;
  for (let i = 0; i < x.length; i++) s += Math.log(d.pdf(x[i]));
  return s;
}

/**
 * Sample L-moments from the unbiased probability-weighted moments b_r.
 * λ1 = b0, λ2 = 2b1 - b0, λ3 = 6b2 - 6b1 + b0, λ4 = 20b3 - 30b2 + 12b1 - b0.
 * Reference: Hosking (1990), J. R. Stat. Soc. B 52, 105–124; Hosking & Wallis (1997), "Regional Frequency
 * Analysis", §2.7. Validity: n >= 4.
 * @param {ArrayLike<number>} x sample [u]
 * @returns {{ l1: number, l2: number, l3: number, l4: number, t3: number, t4: number }} L-moments [u] and the
 *   L-skewness / L-kurtosis ratios [-]
 */
export function lMoments(x) {
  const s = clean(x, 'lMoments', 4).sort();
  const n = s.length;
  let b0 = 0;
  let b1 = 0;
  let b2 = 0;
  let b3 = 0;
  for (let i = 0; i < n; i++) {
    const v = s[i];
    b0 += v;
    b1 += (v * i) / (n - 1);
    b2 += (v * i * (i - 1)) / ((n - 1) * (n - 2));
    b3 += (v * i * (i - 1) * (i - 2)) / ((n - 1) * (n - 2) * (n - 3));
  }
  b0 /= n;
  b1 /= n;
  b2 /= n;
  b3 /= n;
  const l1 = b0;
  const l2 = 2 * b1 - b0;
  const l3 = 6 * b2 - 6 * b1 + b0;
  const l4 = 20 * b3 - 30 * b2 + 12 * b1 - b0;
  return { l1, l2, l3, l4, t3: l2 !== 0 ? l3 / l2 : 0, t4: l2 !== 0 ? l4 / l2 : 0 };
}

function meanStd(x) {
  const n = x.length;
  let m = 0;
  for (let i = 0; i < n; i++) m += x[i];
  m /= n;
  let v = 0;
  for (let i = 0; i < n; i++) v += (x[i] - m) * (x[i] - m);
  return { mean: m, std: n > 1 ? Math.sqrt(v / (n - 1)) : 0 };
}

/**
 * Fit a Gumbel (type I, maxima) distribution F = exp(-exp(-(x - loc)/scale)).
 *  - 'mle': maximum likelihood; the scale solves s = x̄ - Σ x e^(-x/s) / Σ e^(-x/s) (Brent), loc = -s ln mean(e^(-x/s)).
 *    Standard errors from the asymptotic covariance: Var(loc) = 1.10867 s²/n, Var(scale) = 0.60793 s²/n.
 *  - 'mom': method of moments, scale = √6·s_x/π, loc = x̄ - γ·scale.
 *  - 'lmom': L-moments, scale = λ2/ln 2, loc = λ1 - γ·scale.
 * Reference: Gumbel (1958), "Statistics of Extremes"; Johnson, Kotz & Balakrishnan (1995), ch. 22 (MLE and
 * its asymptotic covariance); Hosking (1990). Validity: n >= 3 (4 for 'lmom'), non-constant sample.
 * @param {ArrayLike<number>} x sample of maxima [u]
 * @param {{ method?: 'mle' | 'mom' | 'lmom' }} [opts] estimation method (default 'mle')
 * @returns {{ dist: 'gumbel', loc: number, scale: number, se?: { loc: number, scale: number }, logLik: number,
 *   method: string, n: number }} parameters [u], standard errors (MLE only) [u], log-likelihood [-]
 */
export function fitGumbel(x, opts = {}) {
  const method = opts.method ?? 'mle';
  const a = clean(x, 'fitGumbel', 3);
  const n = a.length;
  const { mean: m, std: s } = meanStd(a);
  if (!(s > 0)) throw new RangeError('fitGumbel: x must not be constant');
  let loc;
  let scale;
  let se;
  if (method === 'mom') {
    scale = (Math.sqrt(6) * s) / Math.PI;
    loc = m - EULER_GAMMA * scale;
  } else if (method === 'lmom') {
    const lm = lMoments(a);
    scale = lm.l2 / Math.LN2;
    loc = lm.l1 - EULER_GAMMA * scale;
  } else if (method === 'mle') {
    let xmin = Infinity;
    for (let i = 0; i < n; i++) if (a[i] < xmin) xmin = a[i];
    const f = (sc) => {
      let num = 0;
      let den = 0;
      for (let i = 0; i < n; i++) {
        const e = Math.exp(-(a[i] - xmin) / sc);
        num += a[i] * e;
        den += e;
      }
      return sc - m + num / den;
    };
    const s0 = (Math.sqrt(6) * s) / Math.PI;
    let lo = s0 / 4;
    let hi = s0 * 4;
    for (let k = 0; k < 60 && f(lo) > 0; k++) lo /= 4;
    for (let k = 0; k < 60 && f(hi) < 0; k++) hi *= 4;
    scale = brent(f, lo, hi, { tol: 1e-13 * s0, maxIter: 200 });
    let den = 0;
    for (let i = 0; i < n; i++) den += Math.exp(-(a[i] - xmin) / scale);
    loc = xmin - scale * Math.log(den / n);
    se = { loc: scale * Math.sqrt(1.10867 / n), scale: scale * Math.sqrt(0.60793 / n) };
  } else throw new RangeError("fitGumbel: method must be 'mle', 'mom' or 'lmom'");
  const out = { dist: 'gumbel', loc, scale, logLik: logLikOf({ dist: 'gumbel', loc, scale }, a), method, n };
  if (se) out.se = se;
  return out;
}

function weibull2Mle(y, kStart) {
  // y > 0; returns { scale, shape } by solving h(k) = Σ y^k ln y / Σ y^k - 1/k - mean(ln y) = 0.
  // h is increasing with h' = Var_w(ln y) + 1/k² > 0 (weights y^k), so a safeguarded Newton iteration converges
  // from any positive start; the usual start is k0 = (π/√6)/sd(ln y).
  const n = y.length;
  let ymax = 0;
  for (let i = 0; i < n; i++) if (y[i] > ymax) ymax = y[i];
  const ln = new Float64Array(n);
  let mln = 0;
  for (let i = 0; i < n; i++) {
    ln[i] = Math.log(y[i] / ymax);
    mln += ln[i];
  }
  mln /= n;
  let k = kStart;
  if (!(k > 0)) {
    let v = 0;
    for (let i = 0; i < n; i++) v += (ln[i] - mln) * (ln[i] - mln);
    v = Math.sqrt(v / Math.max(n - 1, 1));
    k = v > 0 ? 1.2825498301618641 / v : 1;
  }
  let lo = 0;
  let hi = Infinity;
  let den = 0;
  for (let it = 0; it < 100; it++) {
    let s0 = 0;
    let s1 = 0;
    let s2 = 0;
    for (let i = 0; i < n; i++) {
      const l = ln[i];
      const e = Math.exp(k * l);
      s0 += e;
      s1 += e * l;
      s2 += e * l * l;
    }
    den = s0;
    const m1 = s1 / s0;
    const h = m1 - 1 / k - mln;
    if (h > 0) hi = k;
    else lo = k;
    const dh = s2 / s0 - m1 * m1 + 1 / (k * k);
    let kn = k - h / dh;
    if (!(kn > lo && kn < hi)) kn = Number.isFinite(hi) ? 0.5 * (lo + hi) : 2 * k;
    if (Math.abs(kn - k) <= 1e-13 * k) {
      k = kn;
      break;
    }
    k = kn;
    if (!(k < 1e6)) throw new RangeError('fitWeibull: shape estimate diverged (sample is nearly constant)');
  }
  // Σ (y/ymax)^k at the final k
  den = 0;
  let sl = 0;
  for (let i = 0; i < n; i++) {
    den += Math.exp(k * ln[i]);
    sl += ln[i];
  }
  const scale = ymax * Math.pow(den / n, 1 / k);
  // log-likelihood: n ln k - n k ln(scale) + (k - 1) Σ ln y - Σ (y/scale)^k, with Σ (y/scale)^k = n at the optimum
  const logLik = n * Math.log(k) - n * k * Math.log(scale) + (k - 1) * (sl + n * Math.log(ymax)) - n;
  return { scale, shape: k, logLik };
}

/**
 * Fit a Weibull distribution F = 1 - exp(-((x - loc)/scale)^shape) by maximum likelihood.
 * Two-parameter (loc = 0): the shape solves the profile likelihood equation (safeguarded Newton),
 * scale = (mean x^k)^(1/k).
 * Three-parameter: the location is found by maximising the profile log-likelihood over loc < min(x) (log-spaced
 * scan followed by golden-section refinement), with the two-parameter MLE applied to x - loc.
 * Reference: Johnson, Kotz & Balakrishnan (1994), ch. 21; Cohen (1965), Technometrics 7, 579–588;
 * Smith (1985), Biometrika 72, 67–90 (regularity requires shape > 2; the likelihood is unbounded for shape < 1,
 * in which case the scan keeps loc at least 1e-6 of the data range below the minimum).
 * Validity: x > 0 for the two-parameter fit; n >= 5 (>= 10 recommended for three parameters).
 * @param {ArrayLike<number>} x sample [u]
 * @param {{ threeParam?: boolean }} [opts] estimate the location parameter (default false)
 * @returns {{ dist: 'weibull', scale: number, shape: number, loc: number, logLik: number, n: number }}
 *   scale [u], shape [-], loc [u], log-likelihood [-]
 */
export function fitWeibull(x, opts = {}) {
  const a = clean(x, 'fitWeibull', 5);
  const n = a.length;
  let xmin = Infinity;
  let xmax = -Infinity;
  for (let i = 0; i < n; i++) {
    if (a[i] < xmin) xmin = a[i];
    if (a[i] > xmax) xmax = a[i];
  }
  if (!(xmax > xmin)) throw new RangeError('fitWeibull: x must not be constant');
  if (!opts.threeParam) {
    if (!(xmin > 0)) throw new RangeError('fitWeibull: x must be > 0 for the two-parameter fit');
    const r = weibull2Mle(a);
    return { dist: 'weibull', scale: r.scale, shape: r.shape, loc: 0, logLik: r.logLik, n };
  }
  const range = xmax - xmin;
  const y = new Float64Array(n);
  let kWarm;
  const prof = (t) => {
    const loc = xmin - Math.exp(t);
    for (let i = 0; i < n; i++) y[i] = a[i] - loc;
    try {
      const r = weibull2Mle(y, kWarm);
      if (!Number.isFinite(r.logLik)) return null;
      kWarm = r.shape;
      return r;
    } catch (e) {
      return null;
    }
  };
  const tLo = Math.log(1e-6 * range);
  const tHi = Math.log(20 * range);
  const nGrid = 33;
  let best = -Infinity;
  let ib = -1;
  const ts = new Array(nGrid);
  for (let i = 0; i < nGrid; i++) {
    ts[i] = tLo + ((tHi - tLo) * i) / (nGrid - 1);
    const r = prof(ts[i]);
    const v = r ? r.logLik : -Infinity;
    if (v > best) {
      best = v;
      ib = i;
    }
  }
  if (ib < 0) throw new Error('fitWeibull: three-parameter likelihood could not be evaluated');
  const a0 = ts[Math.max(ib - 1, 0)];
  const b0 = ts[Math.min(ib + 1, nGrid - 1)];
  const gs = goldenSection(
    (t) => {
      const r = prof(t);
      return r ? -r.logLik : Infinity;
    },
    a0,
    b0,
    { tol: 1e-6 },
  );
  const tBest = -gs.f >= best ? gs.x : ts[ib];
  const r = prof(tBest);
  return { dist: 'weibull', scale: r.scale, shape: r.shape, loc: xmin - Math.exp(tBest), logLik: r.logLik, n };
}

/**
 * Fit a generalised extreme-value distribution by L-moments.
 * In Hosking's notation (k = -ξ): τ3 = 2(1 - 3^-k)/(1 - 2^-k) - 3 is solved for k (starting from Hosking's
 * approximation k ≈ 7.8590c + 2.9554c², c = 2/(3 + τ3) - ln 2/ln 3, refined by Brent/secant to full precision),
 * then scale = λ2 k / ((1 - 2^-k) Γ(1 + k)) and loc = λ1 - scale (1 - Γ(1 + k))/k.
 * Reference: Hosking, Wallis & Wood (1985), Technometrics 27, 251–261; Hosking & Wallis (1997), App. A.7.
 * Validity: n >= 4 (>= 20 recommended); -1 < k, i.e. ξ < 1.
 * @param {ArrayLike<number>} x sample of block maxima [u]
 * @returns {{ dist: 'gev', loc: number, scale: number, shape: number, logLik: number, n: number }}
 *   loc [u], scale [u], shape ξ [-] (ξ > 0 heavy tail), log-likelihood [-]
 */
export function fitGev(x) {
  const a = clean(x, 'fitGev', 4);
  const lm = lMoments(a);
  if (!(lm.l2 > 0)) throw new RangeError('fitGev: x must not be constant');
  const t3 = lm.t3;
  const c = 2 / (3 + t3) - Math.LN2 / Math.log(3);
  let k = 7.859 * c + 2.9554 * c * c;
  const ratio = (kk) => (Math.abs(kk) < 1e-9 ? Math.log(3) / Math.LN2 : Math.expm1(-kk * Math.log(3)) / Math.expm1(-kk * Math.LN2));
  const f = (kk) => 2 * ratio(kk) - 3 - t3;
  // secant refinement of Hosking's approximation
  let k0 = k;
  let f0 = f(k0);
  let k1 = k + (Math.abs(k) > 1e-3 ? 1e-3 * k : 1e-3);
  let f1 = f(k1);
  for (let it = 0; it < 50 && Math.abs(f1) > 1e-15 && f1 !== f0; it++) {
    const k2 = k1 - (f1 * (k1 - k0)) / (f1 - f0);
    k0 = k1;
    f0 = f1;
    k1 = k2;
    f1 = f(k1);
    if (!Number.isFinite(k1)) break;
  }
  if (Number.isFinite(k1) && Math.abs(f1) < 1e-9) k = k1;
  if (!(k > -1)) throw new RangeError('fitGev: sample L-skewness is outside the range of the GEV family (ξ >= 1)');
  let scale;
  let loc;
  if (Math.abs(k) < 1e-8) {
    scale = lm.l2 / Math.LN2;
    loc = lm.l1 - EULER_GAMMA * scale;
    k = 0;
  } else {
    const g = gam(1 + k);
    scale = (lm.l2 * k) / (-Math.expm1(-k * Math.LN2) * g);
    loc = lm.l1 - (scale * (1 - g)) / k;
  }
  const spec = { dist: 'gev', loc, scale, shape: k === 0 ? 0 : -k };
  return { ...spec, logLik: logLikOf(spec, a), n: a.length };
}

/**
 * Fit a generalised Pareto distribution F = 1 - (1 + ξ y/scale)^(-1/ξ) to threshold excesses y >= 0.
 *  - 'lmom' (default): ξ = 2 - λ1/λ2, scale = λ1 (1 - ξ) (Hosking & Wallis 1987).
 *  - 'mle': maximum likelihood by Nelder–Mead on (ln scale, ξ) started from the L-moment estimate.
 * Reference: Hosking & Wallis (1987), Technometrics 29, 339–349; Coles (2001), §4.3. Validity: n >= 4 (>= 30
 * recommended); L-moment estimator requires ξ < 1, MLE is regular for ξ > -1/2.
 * @param {ArrayLike<number>} excess excesses over the threshold (value - threshold) [u]
 * @param {{ method?: 'lmom' | 'mle' }} [opts] estimation method
 * @returns {{ dist: 'gpd', loc: 0, scale: number, shape: number, logLik: number, method: string, n: number }}
 *   scale [u], shape ξ [-], log-likelihood [-]. Set `loc` to the threshold before calling `returnValue`.
 */
export function fitGpd(excess, opts = {}) {
  const method = opts.method ?? 'lmom';
  const a = clean(excess, 'fitGpd', 4);
  const n = a.length;
  for (let i = 0; i < n; i++) if (a[i] < 0) throw new RangeError('fitGpd: excesses must be >= 0 (subtract the threshold first)');
  const lm = lMoments(a);
  if (!(lm.l2 > 0)) throw new RangeError('fitGpd: excesses must not be constant');
  let shape = 2 - lm.l1 / lm.l2;
  let scale = lm.l1 * (1 - shape);
  if (method === 'mle') {
    let ymax = 0;
    for (let i = 0; i < n; i++) if (a[i] > ymax) ymax = a[i];
    const nll = (p) => {
      const sc = Math.exp(p[0]);
      const xi = p[1];
      if (xi < 0 && ymax >= -sc / xi) return 1e300;
      let s = n * p[0];
      if (Math.abs(xi) < 1e-9) {
        for (let i = 0; i < n; i++) s += a[i] / sc;
      } else {
        for (let i = 0; i < n; i++) s += (1 / xi + 1) * Math.log1p((xi * a[i]) / sc);
      }
      return s;
    };
    let start = [Math.log(scale), shape];
    if (!(nll(start) < 1e299)) start = [Math.log(lm.l1), 0.01];
    const r = nelderMead(nll, start, { tol: 1e-12, maxIter: 2000, step: [0.2, 0.1] });
    scale = Math.exp(r.x[0]);
    shape = r.x[1];
  } else if (method !== 'lmom') throw new RangeError("fitGpd: method must be 'lmom' or 'mle'");
  const spec = { dist: 'gpd', loc: 0, scale, shape };
  return { ...spec, logLik: logLikOf(spec, a), method, n };
}

/**
 * Fit a lognormal distribution: mu and sigma are the mean and (n-1) standard deviation of ln x.
 * Reference: Johnson, Kotz & Balakrishnan (1994), ch. 14. Validity: x > 0, n >= 2.
 * @param {ArrayLike<number>} x sample [u]
 * @returns {{ dist: 'lognormal', mu: number, sigma: number, mean: number, std: number, logLik: number, n: number }}
 *   mu, sigma of ln X [-]; implied mean and standard deviation of X [u]; log-likelihood [-]
 */
export function fitLognormal(x) {
  const a = clean(x, 'fitLognormal', 2);
  const ln = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    if (!(a[i] > 0)) throw new RangeError('fitLognormal: x must be > 0');
    ln[i] = Math.log(a[i]);
  }
  const { mean: mu, std: sigma } = meanStd(ln);
  if (!(sigma > 0)) throw new RangeError('fitLognormal: x must not be constant');
  const m = Math.exp(mu + 0.5 * sigma * sigma);
  const s = m * Math.sqrt(Math.expm1(sigma * sigma));
  return { dist: 'lognormal', mu, sigma, mean: m, std: s, logLik: logLikOf({ dist: 'lognormal', mu, sigma }, a), n: a.length };
}

/**
 * Fit a normal distribution: sample mean and (n-1) standard deviation.
 * @param {ArrayLike<number>} x sample [u]
 * @returns {{ dist: 'normal', mean: number, std: number, mu: number, sigma: number, logLik: number, n: number }}
 *   mean and standard deviation [u] (mu, sigma are aliases); log-likelihood [-]
 */
export function fitNormal(x) {
  const a = clean(x, 'fitNormal', 2);
  const { mean: m, std: s } = meanStd(a);
  if (!(s > 0)) throw new RangeError('fitNormal: x must not be constant');
  return { dist: 'normal', mean: m, std: s, mu: m, sigma: s, logLik: logLikOf({ dist: 'normal', mean: m, std: s }, a), n: a.length };
}

const MS_PER_HOUR = 3600e3;
const MS_PER_YEAR = 365.2425 * 24 * MS_PER_HOUR;

function toMs(t, name, i) {
  let v;
  if (typeof t === 'number') v = t;
  else if (t instanceof Date) v = t.getTime();
  else if (typeof t === 'string') v = Date.parse(t);
  else v = NaN;
  if (!Number.isFinite(v)) throw new RangeError(`${name}: times[${i}] is not a valid time`);
  return v;
}

/**
 * Annual maxima of a time series, grouped by UTC calendar year.
 * Reference: block-maxima approach, Coles (2001), §3.1.
 * @param {ArrayLike<number | string | Date>} times sample times, ms since the Unix epoch or ISO-8601 strings
 * @param {ArrayLike<number>} values observations [u] (non-finite values are skipped)
 * @returns {{ years: number[], maxima: number[], counts: number[], timesOfMax: number[] }} calendar years
 *   (ascending), maximum per year [u], number of valid observations per year, time of each maximum [ms]
 */
export function annualMaxima(times, values) {
  if (!times || !values || times.length !== values.length) throw new RangeError('annualMaxima: times and values must have the same length');
  const map = new Map();
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (!Number.isFinite(v)) continue;
    const t = toMs(times[i], 'annualMaxima', i);
    const y = new Date(t).getUTCFullYear();
    const e = map.get(y);
    if (!e) map.set(y, { max: v, count: 1, t });
    else {
      e.count++;
      if (v > e.max) {
        e.max = v;
        e.t = t;
      }
    }
  }
  const years = Array.from(map.keys()).sort((a, b) => a - b);
  return {
    years,
    maxima: years.map((y) => map.get(y).max),
    counts: years.map((y) => map.get(y).count),
    timesOfMax: years.map((y) => map.get(y).t),
  };
}

/**
 * Peaks over threshold with runs declustering: exceedances closer than `minSeparationHours` belong to the same
 * storm cluster and only the cluster maximum is retained, giving approximately independent peaks.
 * ratePerYear = number of peaks / record length in years (365.2425-day years, first to last sample).
 * Reference: Coles (2001), §4.4 and §5.3.2 (declustering). Validity: threshold high enough for
 * the GPD limit and low enough to keep >~ 30 peaks; record length > 0.
 * @param {ArrayLike<number | string | Date>} times sample times, ms epoch or ISO strings (any order)
 * @param {ArrayLike<number>} values observations [u]
 * @param {{ threshold: number, minSeparationHours?: number }} opts threshold [u]; minimum separation between
 *   clusters [h] (default 48)
 * @returns {{ peaks: number[], times: number[], ratePerYear: number, excesses: number[], threshold: number,
 *   years: number }} cluster maxima [u], their times [ms], mean number of peaks per year [1/yr], peaks minus
 *   threshold [u], record length [yr]
 */
export function peaksOverThreshold(times, values, opts = {}) {
  if (!times || !values || times.length !== values.length) throw new RangeError('peaksOverThreshold: times and values must have the same length');
  const threshold = opts.threshold;
  if (typeof threshold !== 'number' || Number.isNaN(threshold)) throw new RangeError('peaksOverThreshold: threshold is required');
  const sepMs = (opts.minSeparationHours ?? 48) * MS_PER_HOUR;
  if (!(sepMs >= 0)) throw new RangeError('peaksOverThreshold: minSeparationHours must be >= 0');
  const n = values.length;
  if (n < 2) throw new RangeError('peaksOverThreshold: needs at least two samples');
  const t = new Float64Array(n);
  let sorted = true;
  for (let i = 0; i < n; i++) {
    t[i] = toMs(times[i], 'peaksOverThreshold', i);
    if (i > 0 && t[i] < t[i - 1]) sorted = false;
  }
  let order = null;
  if (!sorted) {
    order = new Uint32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
    order.sort((i, j) => t[i] - t[j]);
  }
  const peaks = [];
  const pt = [];
  let last = -Infinity;
  let cur = -Infinity;
  let curT = 0;
  let open = false;
  let tMin = Infinity;
  let tMax = -Infinity;
  for (let k = 0; k < n; k++) {
    const i = order ? order[k] : k;
    const ti = t[i];
    if (ti < tMin) tMin = ti;
    if (ti > tMax) tMax = ti;
    const v = values[i];
    if (!(v > threshold)) continue;
    if (open && ti - last > sepMs) {
      peaks.push(cur);
      pt.push(curT);
      open = false;
    }
    if (!open) {
      open = true;
      cur = v;
      curT = ti;
    } else if (v > cur) {
      cur = v;
      curT = ti;
    }
    last = ti;
  }
  if (open) {
    peaks.push(cur);
    pt.push(curT);
  }
  const years = (tMax - tMin) / MS_PER_YEAR;
  if (!(years > 0)) throw new RangeError('peaksOverThreshold: the record has zero duration');
  return { peaks, times: pt, ratePerYear: peaks.length / years, excesses: peaks.map((p) => p - threshold), threshold, years };
}

/**
 * Return value: the level exceeded on average once every T years, x_T = F⁻¹(1 - 1/(λT)), where λ is the mean
 * number of events per year described by the fitted distribution (1 for annual maxima, the peak rate for a
 * peaks-over-threshold fit, 2922 for 3-hourly sea states, ...).
 * Reference: Coles (2001), §3.1.3 and §4.3.3. Validity: λT > 1.
 * @param {object} fit distribution specification `{ dist, ...params }` (any fit result of this module)
 * @param {number} T return period [yr]
 * @param {{ ratePerYear?: number }} [opts] events per year λ (default 1)
 * @returns {number} return value [u]
 */
export function returnValue(fit, T, opts = {}) {
  if (!fit || typeof fit !== 'object' || !fit.dist) throw new RangeError('returnValue: fit must be { dist, ...params }');
  const rate = opts.ratePerYear ?? 1;
  if (!(T > 0)) throw new RangeError('returnValue: T must be > 0');
  if (!(rate > 0)) throw new RangeError('returnValue: ratePerYear must be > 0');
  if (!(rate * T > 1)) throw new RangeError('returnValue: ratePerYear·T must exceed 1');
  return dist(fit).icdf(1 - 1 / (rate * T));
}

/**
 * Non-parametric bootstrap confidence band of a return value: the sample is resampled with replacement, refitted
 * and the return value recomputed; the band is the percentile interval.
 * Reference: Efron & Tibshirani (1993), "An Introduction to the Bootstrap", ch. 13; Coles (2001), §2.6.
 * Validity: i.i.d. sample (block maxima or declustered peaks); resamples for which the fit fails are skipped.
 * @param {ArrayLike<number>} x sample [u]
 * @param {(x: Float64Array) => object} fitFn fitting function returning `{ dist, ...params }` (e.g. fitGumbel)
 * @param {number} T return period [yr]
 * @param {{ n?: number, seed?: number, ratePerYear?: number, level?: number }} [opts] resamples (500); seed (1);
 *   events per year (1); confidence level (0.95)
 * @returns {{ lo: number, median: number, hi: number, estimate: number, nValid: number, values: Float64Array }}
 *   lower / median / upper percentile of the bootstrap return values [u], the point estimate from the original
 *   sample [u], number of successful resamples, sorted bootstrap values [u]
 */
export function bootstrapReturnValue(x, fitFn, T, opts = {}) {
  const a = clean(x, 'bootstrapReturnValue', 3);
  if (typeof fitFn !== 'function') throw new RangeError('bootstrapReturnValue: fitFn must be a function');
  const nb = opts.n ?? 500;
  if (!Number.isInteger(nb) || nb < 10) throw new RangeError('bootstrapReturnValue: n must be an integer >= 10');
  const level = opts.level ?? 0.95;
  if (!(level > 0 && level < 1)) throw new RangeError('bootstrapReturnValue: level must be within (0, 1)');
  const rvOpts = { ratePerYear: opts.ratePerYear ?? 1 };
  const estimate = returnValue(fitFn(a), T, rvOpts);
  const rng = makeRng(opts.seed ?? 1);
  const n = a.length;
  const res = new Float64Array(n);
  const vals = [];
  for (let b = 0; b < nb; b++) {
    for (let i = 0; i < n; i++) res[i] = a[rng.int(n)];
    try {
      const v = returnValue(fitFn(res), T, rvOpts);
      if (Number.isFinite(v)) vals.push(v);
    } catch (e) {
      // degenerate resample: skip
    }
  }
  if (vals.length < 5) throw new Error('bootstrapReturnValue: the fit failed on almost every resample');
  const sorted = Float64Array.from(vals).sort();
  const tail = (1 - level) / 2;
  return {
    lo: percentileSorted(sorted, 100 * tail),
    median: percentileSorted(sorted, 50),
    hi: percentileSorted(sorted, 100 * (1 - tail)),
    estimate,
    nValid: vals.length,
    values: sorted,
  };
}

/** Weighted linear least squares y ≈ c0 + c1·f. Returns { c0, c1, sse }. */
function wls(f, y, w) {
  let sw = 0;
  let sf = 0;
  let sy = 0;
  let sff = 0;
  let sfy = 0;
  for (let i = 0; i < f.length; i++) {
    sw += w[i];
    sf += w[i] * f[i];
    sy += w[i] * y[i];
    sff += w[i] * f[i] * f[i];
    sfy += w[i] * f[i] * y[i];
  }
  const den = sw * sff - sf * sf;
  let c0;
  let c1;
  if (Math.abs(den) < 1e-14 * Math.max(sw * sff, 1e-300)) {
    c1 = 0;
    c0 = sy / sw;
  } else {
    c1 = (sw * sfy - sf * sy) / den;
    c0 = (sy - c1 * sf) / sw;
  }
  let sse = 0;
  for (let i = 0; i < f.length; i++) {
    const r = y[i] - c0 - c1 * f[i];
    sse += w[i] * r * r;
  }
  return { c0, c1, sse };
}

/**
 * Environmental contour of significant wave height and peak period by the inverse first-order reliability
 * method (IFORM).
 *
 * Joint model (DNV-RP-C205 §3.6.3, conditional modelling approach after Haver):
 *  - Hs: three-parameter Weibull marginal fitted by maximum likelihood;
 *  - Tp | Hs: lognormal with μ(h) = a0 + a1 h^a2 and σ(h) = b0 + b1 exp(b2 h) (mean and standard deviation of
 *    ln Tp), fitted by weighted non-linear least squares to the per-bin statistics of ln Tp.
 * The contour is the image of the circle of radius β = Φ⁻¹(1 - 1/N) in standard normal space, N = number of sea
 * states in T years, through the Rosenblatt transformation: h = F_Hs⁻¹(Φ(β cos θ)), tp = exp(μ(h) + σ(h) β sin θ).
 * Reference: Winterstein, Ude, Cornell, Bjerager & Haver (1993), "Environmental parameters for extreme response:
 * inverse FORM with omission factors", ICOSSAR-93; Haver (1985), Coastal Eng. 9, 357–370; DNV-RP-C205
 * §3.6.3 (joint Hs–Tp model) and §3.7 (environmental contours). Validity: >= 200 sea states spanning at least ~4 Hs bins of 10+ samples; hs > 0, tp > 0;
 * serial correlation of sea states is neglected (as usual for IFORM contours).
 *
 * @param {{ hs: ArrayLike<number>, tp: ArrayLike<number> }} data simultaneous Hs [m] and Tp [s]
 * @param {number} T return period [yr]
 * @param {{ seaStateHours?: number, nPoints?: number, hsBin?: number, minCount?: number }} [opts] sea-state
 *   duration [h] (1); number of contour points (72); Hs bin width for the conditional fit [m] (range/15);
 *   minimum samples per bin (10, sparse bins are merged)
 * @returns {{ hs: number[], tp: number[], fit: { weibull: { scale: number, shape: number, loc: number },
 *   mu: { a0: number, a1: number, a2: number }, sigma: { b0: number, b1: number, b2: number }, beta: number,
 *   nSeaStates: number, bins: { hs: number, mu: number, sigma: number, count: number }[] }, hsMax: number,
 *   tpAtHsMax: number }} contour points (θ = 0 first, counter-clockwise, not repeated at the end) Hs [m] and
 *   Tp [s]; fitted model; largest Hs on the contour [m] with its Tp [s]
 */
export function iformContour(data, T, opts = {}) {
  if (!data || !data.hs || !data.tp || data.hs.length !== data.tp.length) throw new RangeError('iformContour: hs and tp must be arrays of equal length');
  const n = data.hs.length;
  if (n < 50) throw new RangeError('iformContour: at least 50 sea states are required');
  if (!(T > 0)) throw new RangeError('iformContour: T must be > 0');
  const seaStateHours = opts.seaStateHours ?? 1;
  if (!(seaStateHours > 0)) throw new RangeError('iformContour: seaStateHours must be > 0');
  const nPoints = opts.nPoints ?? 72;
  if (!Number.isInteger(nPoints) || nPoints < 4) throw new RangeError('iformContour: nPoints must be an integer >= 4');
  const hs = Float64Array.from(data.hs);
  const lt = new Float64Array(n);
  let hmin = Infinity;
  let hmax = -Infinity;
  for (let i = 0; i < n; i++) {
    if (!(hs[i] > 0) || !Number.isFinite(hs[i])) throw new RangeError(`iformContour: hs[${i}] must be > 0`);
    if (!(data.tp[i] > 0) || !Number.isFinite(data.tp[i])) throw new RangeError(`iformContour: tp[${i}] must be > 0`);
    lt[i] = Math.log(data.tp[i]);
    if (hs[i] < hmin) hmin = hs[i];
    if (hs[i] > hmax) hmax = hs[i];
  }
  const wb = fitWeibull(hs, { threeParam: true });

  // per-bin statistics of ln Tp
  const order = new Uint32Array(n);
  for (let i = 0; i < n; i++) order[i] = i;
  order.sort((i, j) => hs[i] - hs[j]);
  const binW = opts.hsBin ?? (hmax - hmin) / 15;
  if (!(binW > 0)) throw new RangeError('iformContour: hs must not be constant');
  const minCount = opts.minCount ?? 10;
  const bins = [];
  let start = 0;
  let edge = hmin + binW;
  for (let k = 0; k < n; k++) {
    const last = k === n - 1;
    const nextBeyond = last || hs[order[k + 1]] >= edge;
    if (nextBeyond && (k - start + 1 >= minCount || last)) {
      let c = 0;
      let sh = 0;
      let sl = 0;
      for (let q = start; q <= k; q++) {
        sh += hs[order[q]];
        sl += lt[order[q]];
        c++;
      }
      const mh = sh / c;
      const ml = sl / c;
      let v = 0;
      for (let q = start; q <= k; q++) v += (lt[order[q]] - ml) * (lt[order[q]] - ml);
      bins.push({ hs: mh, mu: ml, sigma: c > 1 ? Math.sqrt(v / (c - 1)) : 0, count: c });
      start = k + 1;
    }
    if (!last) while (hs[order[k + 1]] >= edge) edge += binW;
  }
  // merge an under-populated last bin into its neighbour
  if (bins.length > 1 && bins[bins.length - 1].count < Math.max(3, minCount / 2)) {
    const b = bins.pop();
    const a = bins[bins.length - 1];
    const c = a.count + b.count;
    const mu = (a.mu * a.count + b.mu * b.count) / c;
    const varPooled = ((a.count - 1) * a.sigma * a.sigma + (b.count - 1) * b.sigma * b.sigma + a.count * (a.mu - mu) ** 2 + b.count * (b.mu - mu) ** 2) / (c - 1);
    bins[bins.length - 1] = { hs: (a.hs * a.count + b.hs * b.count) / c, mu, sigma: Math.sqrt(varPooled), count: c };
  }
  if (bins.length < 3) throw new RangeError('iformContour: not enough populated Hs bins for the conditional Tp model (need >= 3)');
  const hb = bins.map((b) => b.hs);
  const wts = bins.map((b) => b.count);
  const sq = wts.map(Math.sqrt);

  // μ(h) = a0 + a1 h^a2: profile over a2 (linear in a0, a1), then Levenberg–Marquardt polish
  const muB = bins.map((b) => b.mu);
  let bestMu = null;
  for (let a2 = 0.05; a2 <= 3.0001; a2 += 0.05) {
    const r = wls(
      hb.map((h) => Math.pow(h, a2)),
      muB,
      wts,
    );
    if (!bestMu || r.sse < bestMu.sse) bestMu = { ...r, a2 };
  }
  let muP = [bestMu.c0, bestMu.c1, bestMu.a2];
  try {
    const lmr = levenbergMarquardt((p) => hb.map((h, i) => sq[i] * (p[0] + p[1] * Math.pow(h, p[2]) - muB[i])), muP, {
      bounds: [[-Infinity, Infinity], [-Infinity, Infinity], [0.01, 5]],
      maxIter: 100,
    });
    if (lmr.cost * 2 <= bestMu.sse && lmr.x.every(Number.isFinite)) muP = lmr.x;
  } catch (e) {
    // keep the profile solution
  }
  // σ(h) = b0 + b1 exp(b2 h), b2 <= 0
  const sgB = bins.map((b) => b.sigma);
  let bestSg = null;
  for (let k = 0; k <= 60; k++) {
    const b2 = -0.01 * Math.pow(300, k / 60); // -0.01 .. -3
    const r = wls(
      hb.map((h) => Math.exp(b2 * h)),
      sgB,
      wts,
    );
    if (!bestSg || r.sse < bestSg.sse) bestSg = { ...r, b2 };
  }
  let sgP = [bestSg.c0, bestSg.c1, bestSg.b2];
  try {
    const lmr = levenbergMarquardt((p) => hb.map((h, i) => sq[i] * (p[0] + p[1] * Math.exp(p[2] * h) - sgB[i])), sgP, {
      bounds: [[-Infinity, Infinity], [-Infinity, Infinity], [-10, 0]],
      maxIter: 100,
    });
    if (lmr.cost * 2 <= bestSg.sse && lmr.x.every(Number.isFinite)) sgP = lmr.x;
  } catch (e) {
    // keep the profile solution
  }
  let sgMinObs = Infinity;
  for (const s of sgB) if (s > 0 && s < sgMinObs) sgMinObs = s;
  const sgFloor = Number.isFinite(sgMinObs) ? 0.25 * sgMinObs : 1e-3;
  const muFn = (h) => muP[0] + muP[1] * Math.pow(Math.max(h, 0), muP[2]);
  const sgFn = (h) => Math.max(sgP[0] + sgP[1] * Math.exp(sgP[2] * h), sgFloor);

  const nSeaStates = (T * 365.25 * 24) / seaStateHours;
  if (!(nSeaStates > 1)) throw new RangeError('iformContour: T must cover more than one sea state');
  const beta = normInv(1 - 1 / nSeaStates);
  const hsD = dist({ dist: 'weibull', scale: wb.scale, shape: wb.shape, loc: wb.loc });
  const outH = new Array(nPoints);
  const outT = new Array(nPoints);
  let hsMax = -Infinity;
  let tpAt = NaN;
  for (let k = 0; k < nPoints; k++) {
    const th = (2 * Math.PI * k) / nPoints;
    const h = hsD.fromZ(beta * Math.cos(th));
    const tp = Math.exp(muFn(h) + sgFn(h) * beta * Math.sin(th));
    outH[k] = h;
    outT[k] = tp;
    if (h > hsMax) {
      hsMax = h;
      tpAt = tp;
    }
  }
  return {
    hs: outH,
    tp: outT,
    fit: {
      weibull: { scale: wb.scale, shape: wb.shape, loc: wb.loc },
      mu: { a0: muP[0], a1: muP[1], a2: muP[2] },
      sigma: { b0: sgP[0], b1: sgP[1], b2: sgP[2] },
      beta,
      nSeaStates,
      bins,
    },
    hsMax,
    tpAtHsMax: tpAt,
  };
}

/**
 * Hs–Tp scatter diagram (joint occurrence table).
 * Reference: joint frequency table of sea states as used for long-term wave statistics, DNV-RP-C205 §3.6.
 * @param {ArrayLike<number>} hs significant wave heights [m]
 * @param {ArrayLike<number>} tp peak periods [s]
 * @param {{ hsBin?: number, tpBin?: number }} [opts] class widths: Hs [m] (0.5), Tp [s] (1)
 * @returns {{ hsEdges: number[], tpEdges: number[], counts: number[][], prob: number[][], n: number }}
 *   class boundaries (bins are [edge_i, edge_{i+1})), counts[iHs][iTp], probabilities (sum 1), samples counted
 */
export function scatterDiagram(hs, tp, opts = {}) {
  if (!hs || !tp || hs.length !== tp.length) throw new RangeError('scatterDiagram: hs and tp must have the same length');
  const hsBin = opts.hsBin ?? 0.5;
  const tpBin = opts.tpBin ?? 1;
  if (!(hsBin > 0)) throw new RangeError('scatterDiagram: hsBin must be > 0');
  if (!(tpBin > 0)) throw new RangeError('scatterDiagram: tpBin must be > 0');
  let hLo = Infinity;
  let hHi = -Infinity;
  let tLo = Infinity;
  let tHi = -Infinity;
  let n = 0;
  for (let i = 0; i < hs.length; i++) {
    if (!Number.isFinite(hs[i]) || !Number.isFinite(tp[i])) continue;
    n++;
    if (hs[i] < hLo) hLo = hs[i];
    if (hs[i] > hHi) hHi = hs[i];
    if (tp[i] < tLo) tLo = tp[i];
    if (tp[i] > tHi) tHi = tp[i];
  }
  if (n === 0) throw new RangeError('scatterDiagram: no finite samples');
  const h0 = Math.floor(hLo / hsBin);
  const t0 = Math.floor(tLo / tpBin);
  const nh = Math.floor(hHi / hsBin) - h0 + 1;
  const nt = Math.floor(tHi / tpBin) - t0 + 1;
  const hsEdges = Array.from({ length: nh + 1 }, (_, i) => (h0 + i) * hsBin);
  const tpEdges = Array.from({ length: nt + 1 }, (_, i) => (t0 + i) * tpBin);
  const counts = Array.from({ length: nh }, () => new Array(nt).fill(0));
  for (let i = 0; i < hs.length; i++) {
    if (!Number.isFinite(hs[i]) || !Number.isFinite(tp[i])) continue;
    const ih = Math.min(Math.max(Math.floor(hs[i] / hsBin) - h0, 0), nh - 1);
    const it = Math.min(Math.max(Math.floor(tp[i] / tpBin) - t0, 0), nt - 1);
    counts[ih][it]++;
  }
  const prob = counts.map((row) => row.map((c) => c / n));
  return { hsEdges, tpEdges, counts, prob, n };
}
