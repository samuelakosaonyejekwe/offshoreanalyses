// Uncertainty and risk: deterministic sensitivity (tornado, spider, data tables, scenarios), seeded Monte Carlo /
// Latin-hypercube simulation with rank correlations, percentile reporting with explicit conventions, VaR / CVaR,
// EMV, variance-based importance, convergence diagnostics, risk-register quantification and value of information.

import { makeRng } from '../math/rng.js';
import { dist } from '../math/stats.js';
import { solveDense } from '../math/linalg.js';
import { latinHypercube, sobolIndices } from '../reliability/sampling.js';
import { requireNumber, requireSeries } from './tvm.js';
import { buildProjectCashflow } from './cashflow.js';

// ---------------------------------------------------------------------------------------------------------------
// Small statistics kept local so that the percentile definition is fixed and documented
// ---------------------------------------------------------------------------------------------------------------

/**
 * Quantile of an ascending-sorted sample by linear interpolation between order statistics (Hyndman–Fan type 7,
 * the spreadsheet PERCENTILE.INC definition).
 * Ref: Hyndman, R.J. & Fan, Y. (1996), Sample quantiles in statistical packages, The American Statistician 50(4).
 * @param {ArrayLike<number>} sorted ascending values
 * @param {number} q probability 0–1 (cumulative, non-exceedance)
 * @returns {number} same unit as the sample
 */
export function quantile(sorted, q) {
  requireNumber('q', q, { min: 0, max: 1 });
  const n = sorted.length;
  if (n === 0) throw new RangeError('sample must not be empty');
  const h = (n - 1) * q; const lo = Math.floor(h); const hi = Math.min(n - 1, lo + 1);
  return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
}

const SQRT2PI = Math.sqrt(2 * Math.PI);

/**
 * Standard normal density φ(x).
 * Ref: Abramowitz & Stegun, Handbook of Mathematical Functions, 26.2.1.
 * @param {number} x standard normal variate (-)
 * @returns {number}
 */
export function normPdf(x) { return Math.exp(-0.5 * x * x) / SQRT2PI; }

/**
 * Standard normal cumulative distribution Φ(x) to double precision: power series
 * Φ = 1/2 + φ(x) Σ x^(2k+1) / (1·3·5···(2k+1)) for |x| < 5, Laplace continued fraction for the tails.
 * Ref: Abramowitz & Stegun 26.2.10–26.2.14; Marsaglia, G. (2004), Evaluating the normal distribution,
 * J. Statistical Software 11(4).
 * @param {number} x standard normal variate (-)
 * @returns {number} probability 0–1
 */
export function normCdf(x) {
  if (Number.isNaN(x)) return NaN;
  const ax = Math.abs(x);
  let upper;
  if (ax < 5) {
    let s = ax; let t = ax;
    for (let i = 1; i < 400; i++) { t *= ax * ax / (2 * i + 1); s += t; if (t < s * 1e-17) break; }
    upper = 0.5 - s * normPdf(ax);
  } else if (ax > 40) upper = 0;
  else {
    let f = ax;
    for (let k = 80; k >= 1; k--) f = ax + k / f;
    upper = normPdf(ax) / f;
  }
  return x > 0 ? 1 - upper : upper;
}

/**
 * Inverse standard normal cumulative distribution: rational starting approximation refined by Newton steps on
 * the exact Φ, so the result is accurate to double precision.
 * Ref: Acklam, P.J. (2003), An algorithm for computing the inverse normal cumulative distribution function;
 * refinement by Newton's method.
 * @param {number} p probability in (0, 1)
 * @returns {number} standard normal variate (-)
 */
export function normInv(p) {
  requireNumber('p', p, { gt: 0, lt: 1 });
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  let x;
  if (p < 0.02425 || p > 1 - 0.02425) {
    const q = Math.sqrt(-2 * Math.log(p < 0.5 ? p : 1 - p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    if (p > 0.5) x = -x;
  } else {
    const q = p - 0.5; const r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  for (let it = 0; it < 4; it++) {
    const e = normCdf(x) - p;
    const dx = e / normPdf(x);
    if (!Number.isFinite(dx)) break;
    x -= dx / (1 + 0.5 * x * dx);
    if (Math.abs(dx) < 1e-15 * (1 + Math.abs(x))) break;
  }
  return x;
}

function ranks(a) {
  const n = a.length;
  const idx = Array.from({ length: n }, (_, i) => i).sort((p, q) => a[p] - a[q]);
  const r = new Float64Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    while (j + 1 < n && a[idx[j + 1]] === a[idx[i]]) j++;
    const avg = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) r[idx[k]] = avg;
    i = j + 1;
  }
  return r;
}

function pearson(a, b) {
  const n = a.length; let ma = 0; let mb = 0;
  for (let k = 0; k < n; k++) { ma += a[k]; mb += b[k]; }
  ma /= n; mb /= n;
  let sab = 0; let saa = 0; let sbb = 0;
  for (let k = 0; k < n; k++) { const x = a[k] - ma; const y = b[k] - mb; sab += x * y; saa += x * x; sbb += y * y; }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : 0;
}

/**
 * Spearman rank correlation of two samples (average ranks for ties).
 * Ref: Spearman, C. (1904), Am. J. Psychology 15(1).
 * @param {ArrayLike<number>} a sample
 * @param {ArrayLike<number>} b sample of the same length
 * @returns {number} −1 … 1
 */
export function spearman(a, b) {
  if (a.length !== b.length || a.length < 2) throw new RangeError('samples must have equal length >= 2');
  return pearson(ranks(a), ranks(b));
}

/**
 * Summary statistics of a sample with percentiles reported under both conventions, stated explicitly:
 * `cumulative.pXX` is the value NOT exceeded with XX % probability (cost convention); `exceedance.PXX` is the
 * value exceeded with XX % probability (resource/energy convention, so exceedance.P90 = cumulative.p10).
 * Ref: SPE-PRMS (2018) for the exceedance convention; AACE RP 41R-08 for the cumulative cost convention.
 * @param {ArrayLike<number>} samples values (any unit)
 * @param {{kind?: 'resource' | 'cost' | 'value'}} [o] selects which convention `headline` uses (resource and
 *   value → exceedance, cost → cumulative)
 * @returns {{n: number, mean: number, std: number, standardError: number, min: number, max: number,
 *   cumulative: Record<string, number>, exceedance: Record<string, number>, headline: {P10: number, P50: number,
 *   P90: number, convention: string}, convention: string}}
 */
export function summarise(samples, { kind = 'value' } = {}) {
  const n = samples.length;
  if (n < 2) throw new RangeError('samples must contain at least 2 values');
  const sorted = Float64Array.from(samples).sort();
  let mean = 0; for (let k = 0; k < n; k++) mean += sorted[k]; mean /= n;
  let v = 0; for (let k = 0; k < n; k++) v += (sorted[k] - mean) ** 2;
  const std = Math.sqrt(v / (n - 1));
  const cumulative = {}; const exceedance = {};
  for (const p of [1, 5, 10, 20, 25, 50, 75, 80, 90, 95, 99]) {
    cumulative[`p${p}`] = quantile(sorted, p / 100);
    exceedance[`P${p}`] = quantile(sorted, 1 - p / 100);
  }
  const exc = kind !== 'cost';
  const headline = exc
    ? { P10: exceedance.P10, P50: exceedance.P50, P90: exceedance.P90,
      convention: 'Exceedance: P90 is the value exceeded with 90 % probability (low case), P10 the high case.' }
    : { P10: cumulative.p10, P50: cumulative.p50, P90: cumulative.p90,
      convention: 'Cumulative: P90 is the cost not exceeded with 90 % probability (high case), P10 the low case.' };
  return { n, mean, std, standardError: std / Math.sqrt(n), min: sorted[0], max: sorted[n - 1], cumulative, exceedance,
    headline, convention: 'cumulative.pXX = non-exceedance percentile; exceedance.PXX = value exceeded with XX % probability.' };
}

// ---------------------------------------------------------------------------------------------------------------
// Deterministic sensitivity
// ---------------------------------------------------------------------------------------------------------------

/**
 * One-at-a-time sensitivity: each input is moved to its low and high value with all others at base; the result is
 * tornado data sorted by swing.
 * Ref: Clemen & Reilly, Making Hard Decisions, ch. 5 (tornado diagrams).
 * @param {{model: (x: object) => number, base: Record<string, number>, ranges: Record<string, {low: number,
 *   high: number}>}} o model of named inputs; ranges give absolute low/high input values
 * @returns {{base: number, items: {name: string, low: number, high: number, lowInput: number, highInput: number,
 *   swing: number}[]}} low/high are the model outputs at the low/high input
 */
export function oneAtATime({ model, base, ranges }) {
  if (typeof model !== 'function') throw new RangeError('model must be a function');
  const b = model({ ...base });
  const items = [];
  for (const name of Object.keys(ranges)) {
    const r = ranges[name];
    requireNumber(`ranges.${name}.low`, r.low); requireNumber(`ranges.${name}.high`, r.high);
    const low = model({ ...base, [name]: r.low }); const high = model({ ...base, [name]: r.high });
    items.push({ name, low, high, lowInput: r.low, highInput: r.high, swing: Math.abs(high - low) });
  }
  items.sort((p, q) => q.swing - p.swing);
  return { base: b, items };
}

/**
 * Spider-plot data: model output against a common set of multipliers applied to each selected input in turn.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 18 (sensitivity graphs).
 * @param {{model: (x: object) => number, base: Record<string, number>, names: string[], factors?: number[]}} o
 *   factors = multipliers on the base value (default 0.7 … 1.3)
 * @returns {{base: number, series: {name: string, x: number[], y: number[]}[]}} x = factor, y = output
 */
export function spider({ model, base, names, factors = [0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3] }) {
  requireSeries('factors', factors);
  const series = names.map((name) => {
    requireNumber(`base.${name}`, base[name]);
    return { name, x: Array.from(factors), y: Array.from(factors, (f) => model({ ...base, [name]: base[name] * f })) };
  });
  return { base: model({ ...base }), series };
}

/**
 * Two-way data table: model output on the grid of two inputs, z[iy][ix].
 * Ref: spreadsheet two-variable data table; Blank & Tarquin ch. 18.
 * @param {{model: (x: object) => number, base: object, x: {name: string, values: number[]}, y: {name: string,
 *   values: number[]}}} o
 * @returns {{xs: number[], ys: number[], z: number[][]}}
 */
export function dataTable2({ model, base, x, y }) {
  requireSeries('x.values', x.values); requireSeries('y.values', y.values);
  const z = Array.from(y.values, (yv) => Array.from(x.values, (xv) => model({ ...base, [x.name]: xv, [y.name]: yv })));
  return { xs: Array.from(x.values), ys: Array.from(y.values), z };
}

/**
 * Scenario manager: evaluates the model for the base case and for named scenarios given as input overrides.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 18 (scenario analysis).
 * @param {{model: (x: object) => number | object, base: object, scenarios: Record<string, object>}} o e.g.
 *   scenarios = { low: { price: 40 }, high: { price: 80 } }
 * @returns {{name: string, inputs: object, result: number | object}[]} base first
 */
export function scenarios({ model, base, scenarios: sc }) {
  const out = [{ name: 'base', inputs: { ...base }, result: model({ ...base }) }];
  for (const name of Object.keys(sc || {})) {
    const inputs = { ...base, ...sc[name] };
    out.push({ name, inputs, result: model(inputs) });
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Sampling with rank correlation
// ---------------------------------------------------------------------------------------------------------------

function choleskyLower(C, what) {
  const k = C.length; const L = Array.from({ length: k }, () => new Array(k).fill(0));
  for (let i = 0; i < k; i++) {
    for (let j = 0; j <= i; j++) {
      let s = C[i][j];
      for (let m = 0; m < j; m++) s -= L[i][m] * L[j][m];
      if (i === j) {
        if (!(s > 1e-12)) throw new RangeError(`${what} is not positive definite`);
        L[i][i] = Math.sqrt(s);
      } else L[i][j] = s / L[j][j];
    }
  }
  return L;
}

/**
 * Impose a target Spearman rank-correlation matrix on independent sample columns by the Iman–Conover
 * restricted-pairing method. Only the pairing is changed: every marginal sample is preserved exactly.
 * Ref: Iman, R.L. & Conover, W.J. (1982), A distribution-free approach to inducing rank correlation among input
 * variables, Communications in Statistics B11(3).
 * @param {Record<string, Float64Array>} samples columns of equal length n
 * @param {{a: string, b: string, rho: number}[]} rankCorr target rank correlations between named columns
 * @param {number} [seed] seed for the score permutation (default 1)
 * @returns {Record<string, Float64Array>} re-paired columns
 */
export function imposeRankCorrelation(samples, rankCorr, seed = 1) {
  const names = Object.keys(samples); const k = names.length;
  if (!rankCorr || rankCorr.length === 0 || k < 2) return samples;
  const n = samples[names[0]].length;
  const C = Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => (i === j ? 1 : 0)));
  for (const c of rankCorr) {
    const i = names.indexOf(c.a); const j = names.indexOf(c.b);
    if (i < 0 || j < 0 || i === j) throw new RangeError(`rankCorr refers to unknown or identical variables '${c.a}', '${c.b}'`);
    requireNumber('rankCorr.rho', c.rho, { gt: -1, lt: 1 });
    C[i][j] = c.rho; C[j][i] = c.rho;
  }
  const P = choleskyLower(C, 'rank-correlation matrix');
  const rng = makeRng(seed);
  const scores = new Float64Array(n);
  for (let i = 0; i < n; i++) scores[i] = normInv((i + 1) / (n + 1));
  const M = [];
  for (let j = 0; j < k; j++) {
    const col = Float64Array.from(scores);
    for (let i = n - 1; i > 0; i--) { const r = Math.floor(rng.next() * (i + 1)); const t = col[i]; col[i] = col[r]; col[r] = t; }
    M.push(col);
  }
  const T = Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => (i === j ? 1 : pearson(M[i], M[j]))));
  const Q = choleskyLower(T, 'score correlation matrix');
  const Mstar = Array.from({ length: k }, () => new Float64Array(n));
  const y = new Float64Array(k);
  for (let r = 0; r < n; r++) {
    for (let i = 0; i < k; i++) { let s = M[i][r]; for (let m = 0; m < i; m++) s -= Q[i][m] * y[m]; y[i] = s / Q[i][i]; }
    for (let i = 0; i < k; i++) { let s = 0; for (let m = 0; m <= i; m++) s += P[i][m] * y[m]; Mstar[i][r] = s; }
  }
  const out = {};
  for (let j = 0; j < k; j++) {
    const sorted = Float64Array.from(samples[names[j]]).sort();
    const rk = ranks(Mstar[j]);
    const col = new Float64Array(n);
    for (let r = 0; r < n; r++) col[r] = sorted[Math.min(n - 1, Math.round(rk[r]) - 1)];
    out[names[j]] = col;
  }
  return out;
}

/**
 * Draw input samples for named distributions by Latin hypercube (default) or plain Monte Carlo, seeded, with
 * optional Spearman rank correlations imposed by Iman–Conover.
 * Ref: McKay, Beckman & Conover (1979), Technometrics 21(2) (LHS); Iman & Conover (1982).
 * @param {{vars: Record<string, object>, n?: number, seed?: number, method?: 'lhs' | 'mc',
 *   rankCorr?: {a: string, b: string, rho: number}[]}} o vars = { name: distribution spec of math/stats dist() }
 * @returns {{samples: Record<string, Float64Array>, names: string[], means: number[], n: number, seed: number,
 *   method: string}}
 */
export function sampleInputs({ vars, n = 2000, seed = 1, method = 'lhs', rankCorr }) {
  requireNumber('n', n, { min: 2, integer: true }); requireNumber('seed', seed);
  const names = Object.keys(vars || {});
  if (names.length === 0) throw new RangeError('vars must name at least one distribution');
  const dists = names.map((name) => dist(vars[name]));
  let samples;
  if (method === 'lhs') {
    const lh = latinHypercube(vars, n, { seed });
    samples = {};
    for (const name of names) samples[name] = Float64Array.from(lh[name]);
  } else if (method === 'mc') {
    const rng = makeRng(seed);
    samples = {};
    for (const name of names) samples[name] = new Float64Array(n);
    for (let k = 0; k < n; k++) for (let j = 0; j < names.length; j++) samples[names[j]][k] = dists[j].sample(rng);
  } else throw new RangeError("method must be 'lhs' or 'mc'");
  if (rankCorr && rankCorr.length) samples = imposeRankCorrelation(samples, rankCorr, seed + 7919);
  return { samples, names, means: dists.map((d) => d.mean), n, seed, method };
}

/**
 * Monte Carlo / LHS propagation of input distributions through a model returning one number or an object of
 * named numbers. Deterministic for a given seed.
 * Ref: Vose, D., Risk Analysis: A Quantitative Guide, 3rd ed., ch. 4 and 7.
 * @param {{model: (x: object) => number | Record<string, number>, vars: Record<string, object>, fixed?: object,
 *   n?: number, seed?: number, method?: 'lhs' | 'mc', rankCorr?: object[]}} o fixed = constant inputs merged
 *   into every evaluation
 * @returns {{outputs: Record<string, Float64Array>, inputs: Record<string, Float64Array>,
 *   stats: Record<string, object>, n: number, seed: number, method: string}} stats as from summarise(); a scalar
 *   model's output is named 'value'
 */
export function simulate({ model, vars, fixed, n = 2000, seed = 1, method = 'lhs', rankCorr }) {
  if (typeof model !== 'function') throw new RangeError('model must be a function');
  const s = sampleInputs({ vars, n, seed, method, rankCorr });
  const outputs = {};
  for (let k = 0; k < n; k++) {
    const x = { ...(fixed || {}) };
    for (const name of s.names) x[name] = s.samples[name][k];
    let r = model(x);
    if (typeof r === 'number') r = { value: r };
    for (const key of Object.keys(r)) {
      if (typeof r[key] !== 'number') continue;
      (outputs[key] ??= new Float64Array(n).fill(NaN))[k] = r[key];
    }
  }
  const stats = {};
  for (const key of Object.keys(outputs)) {
    const finite = outputs[key].filter((v) => Number.isFinite(v));
    stats[key] = finite.length >= 2 ? { ...summarise(finite), nValid: finite.length } : { nValid: finite.length };
  }
  return { outputs, inputs: s.samples, stats, n, seed, method };
}

// ---------------------------------------------------------------------------------------------------------------
// Risk measures
// ---------------------------------------------------------------------------------------------------------------

/**
 * Probability that a sampled outcome falls below a threshold (default: P[NPV < 0]) with its binomial standard error.
 * Ref: Vose, Risk Analysis, ch. 4.
 * @param {ArrayLike<number>} samples outcomes (e.g. NPV, currency)
 * @param {number} [threshold] default 0
 * @returns {{probability: number, standardError: number}}
 */
export function probabilityBelow(samples, threshold = 0) {
  const n = samples.length;
  if (n === 0) throw new RangeError('samples must not be empty');
  let c = 0; for (let k = 0; k < n; k++) if (samples[k] < threshold) c++;
  const p = c / n;
  return { probability: p, standardError: Math.sqrt(p * (1 - p) / n) };
}

/**
 * Value at risk and conditional value at risk (expected shortfall) of sampled outcomes at confidence α.
 * Losses are L = reference − outcome (reference 0 by default, or the sample mean with relativeTo 'mean');
 * VaR_α is the α-quantile of L and CVaR_α = E[L | L >= VaR_α]. Positive numbers are losses.
 * Ref: Rockafellar, R.T. & Uryasev, S. (2000), Optimization of conditional value-at-risk, J. Risk 2(3);
 * Jorion, Value at Risk, 3rd ed., ch. 5.
 * @param {ArrayLike<number>} samples outcomes where larger is better (e.g. NPV, currency)
 * @param {number} [alpha] confidence level in (0, 1), default 0.95
 * @param {{relativeTo?: 'zero' | 'mean'}} [o]
 * @returns {{alpha: number, var: number, cvar: number, reference: number}} currency
 */
export function valueAtRisk(samples, alpha = 0.95, { relativeTo = 'zero' } = {}) {
  requireNumber('alpha', alpha, { gt: 0, lt: 1 });
  const n = samples.length;
  if (n < 2) throw new RangeError('samples must contain at least 2 values');
  let ref = 0;
  if (relativeTo === 'mean') { for (let k = 0; k < n; k++) ref += samples[k]; ref /= n; }
  const L = new Float64Array(n);
  for (let k = 0; k < n; k++) L[k] = ref - samples[k];
  L.sort();
  const v = quantile(L, alpha);
  let s = 0; let c = 0;
  for (let k = n - 1; k >= 0 && L[k] >= v; k--) { s += L[k]; c++; }
  return { alpha, var: v, cvar: c ? s / c : v, reference: ref };
}

/**
 * Closed-form VaR and CVaR for a normally distributed outcome X ~ N(mean, std): VaR_α = −mean + std z_α,
 * CVaR_α = −mean + std φ(z_α) / (1 − α).
 * Ref: Rockafellar & Uryasev (2000); McNeil, Frey & Embrechts, Quantitative Risk Management, ch. 2.
 * @param {{mean: number, std: number, alpha?: number}} o outcome mean and standard deviation (currency), confidence
 * @returns {{var: number, cvar: number, z: number}}
 */
export function normalVaR({ mean, std, alpha = 0.95 }) {
  requireNumber('mean', mean); requireNumber('std', std, { min: 0 }); requireNumber('alpha', alpha, { gt: 0, lt: 1 });
  const z = normInv(alpha);
  return { var: -mean + std * z, cvar: -mean + std * normPdf(z) / (1 - alpha), z };
}

/**
 * Expected monetary value Σ p_i v_i of discrete outcomes; probabilities must sum to 1.
 * Ref: Clemen & Reilly, Making Hard Decisions, ch. 4.
 * @param {{p: number, value: number}[]} outcomes probabilities (0–1) and values (currency)
 * @returns {{emv: number, std: number}}
 */
export function expectedMonetaryValue(outcomes) {
  if (!Array.isArray(outcomes) || outcomes.length === 0) throw new RangeError('outcomes must be a non-empty array');
  let sp = 0; let m = 0; let m2 = 0;
  for (const o of outcomes) {
    requireNumber('p', o.p, { min: 0, max: 1 }); requireNumber('value', o.value);
    sp += o.p; m += o.p * o.value; m2 += o.p * o.value * o.value;
  }
  if (Math.abs(sp - 1) > 1e-9) throw new RangeError('outcome probabilities must sum to 1');
  return { emv: m, std: Math.sqrt(Math.max(0, m2 - m * m)) };
}

/**
 * Convergence diagnostics of a Monte Carlo sample: running mean and standard error of the mean at increasing
 * sample sizes, in the order the samples were generated.
 * Ref: Vose, Risk Analysis, ch. 7 (how many iterations).
 * @param {ArrayLike<number>} samples values in generation order
 * @param {number} [points] number of checkpoints (default 40)
 * @returns {{n: number[], runningMean: number[], standardError: number[], relativeError: number}} relativeError =
 *   final standard error / |final mean|
 */
export function convergence(samples, points = 40) {
  const N = samples.length;
  if (N < 2) throw new RangeError('samples must contain at least 2 values');
  const n = []; const runningMean = []; const standardError = [];
  let s = 0; let s2 = 0; let next = 1;
  const step = Math.max(1, Math.floor(N / points));
  for (let k = 0; k < N; k++) {
    s += samples[k]; s2 += samples[k] * samples[k];
    const c = k + 1;
    if (c >= 2 && (c === N || c >= next)) {
      const m = s / c; const v = Math.max(0, (s2 - c * m * m) / (c - 1));
      n.push(c); runningMean.push(m); standardError.push(Math.sqrt(v / c));
      next = c + step;
    }
  }
  const last = n.length - 1;
  return { n, runningMean, standardError,
    relativeError: runningMean[last] !== 0 ? standardError[last] / Math.abs(runningMean[last]) : NaN };
}

// ---------------------------------------------------------------------------------------------------------------
// Importance measures
// ---------------------------------------------------------------------------------------------------------------

/**
 * "What matters" ranking: Spearman rank correlation of every sampled input with an output, sorted by absolute
 * value, with each input's share of the summed squared correlations.
 * Ref: Saltelli et al., Global Sensitivity Analysis: The Primer, ch. 1 (rank-based measures).
 * @param {Record<string, ArrayLike<number>>} inputs sampled inputs
 * @param {ArrayLike<number>} output sampled output (same length; non-finite entries are skipped)
 * @returns {{name: string, rankCorrelation: number, share: number}[]}
 */
export function whatMatters(inputs, output) {
  const keep = [];
  for (let k = 0; k < output.length; k++) if (Number.isFinite(output[k])) keep.push(k);
  if (keep.length < 3) throw new RangeError('output needs at least 3 finite samples');
  const y = keep.map((k) => output[k]);
  const res = Object.keys(inputs).map((name) => ({ name, rankCorrelation: spearman(keep.map((k) => inputs[name][k]), y), share: 0 }));
  let s = 0; for (const r of res) s += r.rankCorrelation ** 2;
  for (const r of res) r.share = s > 0 ? r.rankCorrelation ** 2 / s : 0;
  return res.sort((p, q) => Math.abs(q.rankCorrelation) - Math.abs(p.rankCorrelation));
}

/**
 * Variance decomposition by standardised regression coefficients: output regressed linearly on the sampled
 * inputs; SRC_i² is the fraction of output variance attributed to input i (exact for a linear model with
 * independent inputs; R² shows how much the linear model explains).
 * Ref: Saltelli et al., Global Sensitivity Analysis: The Primer, ch. 1.
 * @param {Record<string, ArrayLike<number>>} inputs sampled inputs
 * @param {ArrayLike<number>} output sampled output
 * @returns {{src: Record<string, number>, varianceShare: Record<string, number>, r2: number}}
 */
export function varianceDecomposition(inputs, output) {
  const names = Object.keys(inputs); const k = names.length; const n = output.length;
  if (n <= k + 1) throw new RangeError('not enough samples for the regression');
  const z = (a) => {
    let m = 0; for (let i = 0; i < n; i++) m += a[i]; m /= n;
    let v = 0; for (let i = 0; i < n; i++) v += (a[i] - m) ** 2;
    const s = Math.sqrt(v / (n - 1)) || 1;
    return Float64Array.from(a, (x) => (x - m) / s);
  };
  const X = names.map((nm) => z(inputs[nm])); const y = z(output);
  const A = Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => { let s = 0; for (let r = 0; r < n; r++) s += X[i][r] * X[j][r]; return s; }));
  const b = Array.from({ length: k }, (_, i) => { let s = 0; for (let r = 0; r < n; r++) s += X[i][r] * y[r]; return s; });
  const beta = solveDense(A, b);
  let sse = 0; let sst = 0;
  for (let r = 0; r < n; r++) { let p = 0; for (let i = 0; i < k; i++) p += beta[i] * X[i][r]; sse += (y[r] - p) ** 2; sst += y[r] ** 2; }
  const src = {}; const varianceShare = {};
  names.forEach((nm, i) => { src[nm] = beta[i]; varianceShare[nm] = beta[i] ** 2; });
  return { src, varianceShare, r2: sst > 0 ? 1 - sse / sst : NaN };
}

/**
 * Sobol first-order and total-effect indices of a scalar model (Saltelli sampling, shared library routine).
 * Ref: Sobol', I.M. (2001), Math. Comput. Simul. 55; Saltelli et al. (2010), Comput. Phys. Commun. 181(2).
 * @param {{model: (x: object) => number, vars: Record<string, object>, n?: number, seed?: number}} o n = base
 *   sample size (model runs = n (k + 2))
 * @returns {{first: Record<string, number>, total: Record<string, number>, ranking: string[]}}
 */
export function sobol({ model, vars, n = 1024, seed = 1 }) {
  const r = sobolIndices(model, vars, { n, seed });
  const ranking = Object.keys(r.total).sort((a, b) => r.total[b] - r.total[a]);
  return { first: r.first, total: r.total, ranking };
}

// ---------------------------------------------------------------------------------------------------------------
// Stochastic cash-flow simulation
// ---------------------------------------------------------------------------------------------------------------

const FACTOR_VARS = { priceFactor: 'price', productionFactor: 'production', capexFactor: 'capex', opexFactor: 'opex' };

/**
 * Monte Carlo / LHS simulation of the project cash-flow model. Recognised variable names: priceFactor,
 * productionFactor, capexFactor, opexFactor (multipliers), discountRate, inflation, delayYears (start-up delay;
 * fractional delays interpolate between the neighbouring whole-year models). Any other variable is passed to
 * `apply(inputs, x)` which must return modified inputs.
 * Outputs: npv, irr (NaN when not unique), capex (total), lcoe (pre-tax levelised cost = PV of capex, opex and
 * decommissioning / PV of production at the discount rate), delayYears.
 * Percentile conventions are stated in the result; costs are reported under both.
 * Ref: Vose, Risk Analysis, ch. 4; Savage, S.L. (2009), The Flaw of Averages (why distributions, not point values).
 * @param {{inputs: object, vars: Record<string, object>, apply?: (inputs: object, x: object) => object, n?: number,
 *   seed?: number, method?: 'lhs' | 'mc', rankCorr?: object[], alpha?: number}} o inputs as for
 *   buildProjectCashflow; alpha = VaR confidence (default 0.95)
 * @returns {{outputs: Record<string, Float64Array>, inputs: Record<string, Float64Array>, npv: object, irr: object,
 *   capex: object, lcoe: object, probNpvNegative: number, var: object, ranking: object[], convergence: object,
 *   conventions: Record<string, string>, n: number, seed: number}}
 */
export function simulateCashflow({ inputs, vars, apply, n = 1000, seed = 1, method = 'lhs', rankCorr, alpha = 0.95 }) {
  const evalOnce = (x, delay) => {
    let inp = { ...inputs, factors: { ...(inputs.factors || {}) }, timeline: { ...inputs.timeline } };
    const custom = {};
    for (const k of Object.keys(x)) {
      if (FACTOR_VARS[k]) inp.factors[FACTOR_VARS[k]] = (inp.factors[FACTOR_VARS[k]] ?? 1) * x[k];
      else if (k === 'discountRate') inp.discountRate = x[k];
      else if (k === 'inflation') inp.inflation = x[k];
      else if (k !== 'delayYears') custom[k] = x[k];
    }
    if (delay !== undefined) inp.timeline.startDelayYears = (inputs.timeline.startDelayYears ?? 0) + delay;
    if (apply) inp = apply(inp, x);
    else if (Object.keys(custom).length) throw new RangeError(`no apply() given for variables: ${Object.keys(custom).join(', ')}`);
    const m = buildProjectCashflow(inp);
    const dec = inp.decommissioning?.fund ? m.table.decomFund : m.table.decomSpend;
    let pvDec = 0; for (let t = 0; t < m.n; t++) pvDec += dec[t] * m.table.discountFactor[t];
    return { npv: m.metrics.npv, irr: m.metrics.irr ?? NaN, capex: m.totals.capex,
      lcoe: m.totals.pvProduction > 0 ? (m.totals.pvCapex + m.totals.pvOpex + pvDec) / m.totals.pvProduction : NaN };
  };
  const model = (x) => {
    if (x.delayYears === undefined) return { ...evalOnce(x), delayYears: 0 };
    const d = Math.max(0, x.delayYears); const lo = Math.floor(d); const w = d - lo;
    const a = evalOnce(x, lo);
    if (w < 1e-9) return { ...a, delayYears: d };
    const b = evalOnce(x, lo + 1);
    const mix = (p, q) => p + w * (q - p);
    return { npv: mix(a.npv, b.npv), irr: mix(a.irr, b.irr), capex: mix(a.capex, b.capex), lcoe: mix(a.lcoe, b.lcoe), delayYears: d };
  };
  const sim = simulate({ model, vars, n, seed, method, rankCorr });
  const finite = (a) => a.filter((v) => Number.isFinite(v));
  const sm = (key, kind) => { const f = finite(sim.outputs[key]); return f.length >= 2 ? summarise(f, { kind }) : null; };
  return { outputs: sim.outputs, inputs: sim.inputs, npv: sm('npv', 'value'), irr: sm('irr', 'value'),
    capex: sm('capex', 'cost'), lcoe: sm('lcoe', 'cost'),
    probNpvNegative: probabilityBelow(sim.outputs.npv, 0).probability, var: valueAtRisk(sim.outputs.npv, alpha),
    ranking: whatMatters(sim.inputs, sim.outputs.npv), convergence: convergence(sim.outputs.npv),
    conventions: { npv: 'exceedance (P90 = value exceeded with 90 % probability)',
      cost: 'both reported: cumulative.pXX (not exceeded with XX %) and exceedance.PXX (exceeded with XX %)' },
    n, seed };
}

// ---------------------------------------------------------------------------------------------------------------
// Risk register and value of information
// ---------------------------------------------------------------------------------------------------------------

/**
 * Risk-register quantification: each risk occurs with its probability and, if it occurs, has an impact drawn from
 * its distribution (threats positive cost, opportunities negative). Returns the distribution of total impact and
 * the risk-adjusted contingency at P50 / P80 / P90 (cumulative convention).
 * Ref: AACE International RP 44R-08 (risk analysis and contingency determination using expected value);
 * ISO 31010:2019 Annex B (Monte Carlo simulation).
 * @param {{risks: {name: string, probability: number, impact: number | object}[], n?: number, seed?: number}} o
 *   impact = fixed currency amount or distribution spec (currency)
 * @returns {{expectedValue: number, mean: number, std: number, p50: number, p80: number, p90: number,
 *   byRisk: {name: string, expected: number, share: number}[], samples: Float64Array, convention: string}}
 */
export function riskRegister({ risks, n = 5000, seed = 1 }) {
  if (!Array.isArray(risks) || risks.length === 0) throw new RangeError('risks must be a non-empty array');
  requireNumber('n', n, { min: 2, integer: true });
  const rng = makeRng(seed);
  const ds = risks.map((r) => {
    requireNumber(`probability of '${r.name}'`, r.probability, { min: 0, max: 1 });
    return typeof r.impact === 'number' ? { mean: r.impact, sample: () => r.impact } : dist(r.impact);
  });
  const total = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    let s = 0;
    for (let j = 0; j < risks.length; j++) {
      const u = rng.next();
      const x = ds[j].sample(rng);
      if (u < risks[j].probability) s += x;
    }
    total[k] = s;
  }
  const st = summarise(total, { kind: 'cost' });
  let ev = 0;
  const byRisk = risks.map((r, j) => { const e = r.probability * ds[j].mean; ev += e; return { name: r.name, expected: e, share: 0 }; });
  for (const b of byRisk) b.share = ev !== 0 ? b.expected / ev : 0;
  byRisk.sort((a, b) => Math.abs(b.expected) - Math.abs(a.expected));
  return { expectedValue: ev, mean: st.mean, std: st.std, p50: st.cumulative.p50, p80: st.cumulative.p80,
    p90: st.cumulative.p90, byRisk, samples: total,
    convention: 'Cumulative: Pxx is the total risk cost not exceeded with xx % probability.' };
}

/**
 * Expected value of perfect information for a discrete decision: EVPI = Σ_s p_s max_a v(a, s) − max_a Σ_s p_s v(a, s).
 * Ref: Raiffa & Schlaifer (1961), Applied Statistical Decision Theory; Clemen & Reilly, Making Hard Decisions, ch. 12.
 * @param {{p: number, payoffs: number[]}[]} states state probabilities (sum 1) and payoff of each action (currency)
 * @returns {{evpi: number, evWithoutInformation: number, evWithPerfectInformation: number, bestAction: number}}
 */
export function evpi(states) {
  if (!Array.isArray(states) || states.length === 0) throw new RangeError('states must be a non-empty array');
  const na = states[0].payoffs.length;
  let sp = 0; let withPi = 0; const ev = new Array(na).fill(0);
  for (const s of states) {
    requireNumber('p', s.p, { min: 0, max: 1 }); requireSeries('payoffs', s.payoffs);
    if (s.payoffs.length !== na) throw new RangeError('every state needs one payoff per action');
    sp += s.p; withPi += s.p * Math.max(...s.payoffs);
    for (let a = 0; a < na; a++) ev[a] += s.p * s.payoffs[a];
  }
  if (Math.abs(sp - 1) > 1e-9) throw new RangeError('state probabilities must sum to 1');
  const best = ev.indexOf(Math.max(...ev));
  return { evpi: withPi - ev[best], evWithoutInformation: ev[best], evWithPerfectInformation: withPi, bestAction: best };
}

/**
 * Expected value of imperfect (sample) information by Bayesian pre-posterior analysis: for each possible signal
 * the posterior over states is formed, the best action chosen, and the result averaged over the signal
 * probabilities. EVII = that expectation − best prior expected value (0 <= EVII <= EVPI).
 * Ref: Raiffa & Schlaifer (1961); Clemen & Reilly, Making Hard Decisions, ch. 12.
 * @param {{states: {p: number, payoffs: number[]}[], likelihood: number[][]}} o likelihood[s][k] = P(signal k |
 *   state s), rows sum to 1
 * @returns {{evii: number, evpi: number, evWithInformation: number, evWithoutInformation: number,
 *   signals: {probability: number, posterior: number[], bestAction: number, value: number}[]}}
 */
export function evii({ states, likelihood }) {
  const base = evpi(states);
  const ns = states.length; const na = states[0].payoffs.length;
  if (!Array.isArray(likelihood) || likelihood.length !== ns) throw new RangeError('likelihood needs one row per state');
  const nk = likelihood[0].length;
  for (const row of likelihood) {
    let s = 0; for (const v of row) s += requireNumber('likelihood', v, { min: 0, max: 1 });
    if (row.length !== nk || Math.abs(s - 1) > 1e-9) throw new RangeError('likelihood rows must have equal length and sum to 1');
  }
  const signals = []; let withInfo = 0;
  for (let k = 0; k < nk; k++) {
    let pk = 0; for (let s = 0; s < ns; s++) pk += states[s].p * likelihood[s][k];
    const post = states.map((st, s) => (pk > 0 ? st.p * likelihood[s][k] / pk : 0));
    let bestA = 0; let bestV = -Infinity;
    for (let a = 0; a < na; a++) { let v = 0; for (let s = 0; s < ns; s++) v += post[s] * states[s].payoffs[a]; if (v > bestV) { bestV = v; bestA = a; } }
    if (pk > 0) withInfo += pk * bestV;
    signals.push({ probability: pk, posterior: post, bestAction: bestA, value: pk > 0 ? bestV : 0 });
  }
  return { evii: withInfo - base.evWithoutInformation, evpi: base.evpi, evWithInformation: withInfo,
    evWithoutInformation: base.evWithoutInformation, signals };
}
