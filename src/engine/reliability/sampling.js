/**
 * Sampling-based uncertainty propagation: crude Monte Carlo, Latin hypercube sampling with Iman–Conover rank
 * correlation, variance-based (Sobol) sensitivity indices by the Saltelli scheme, and one-at-a-time tornado
 * sensitivities.
 *
 * Random variables are declared as `{ name: distSpec }` with the specifications understood by `dist()` of
 * math/stats.js; models receive a plain object `{ name: value }`.
 *
 * Correlation input (`corr`) is accepted in any of these forms:
 *   - an object keyed by a pair of names joined with ',', '|', ':' or '~':  { 'Hs,Tp': 0.6 }
 *   - a nested object: { Hs: { Tp: 0.6 } }
 *   - an array of triplets: [['Hs', 'Tp', 0.6]] or [{ a: 'Hs', b: 'Tp', rho: 0.6 }]
 *   - a full symmetric matrix (number[][]) ordered like Object.keys(vars)
 */
import { makeRng } from '../math/rng.js';
import { dist, normCdf, normInv, percentileSorted } from '../math/stats.js';
import { cholesky } from '../math/linalg.js';
import { gaussLegendre } from '../math/quad.js';
import { brent } from '../math/roots.js';

/**
 * Resolve the variable declarations into names and distribution objects.
 * @param {Object<string, object | number>} vars { name: distSpec } (a bare number is a deterministic value)
 * @returns {{ names: string[], dists: object[] }}
 */
export function resolveVars(vars) {
  if (!vars || typeof vars !== 'object') throw new RangeError('vars: must be an object of { name: distSpec }');
  const names = Object.keys(vars);
  if (names.length === 0) throw new RangeError('vars: at least one variable is required');
  const dists = names.map((k) => {
    try {
      return dist(vars[k]);
    } catch (e) {
      throw new RangeError(`vars.${k}: ${e.message}`);
    }
  });
  return { names, dists };
}

/**
 * Build the full correlation matrix from any of the accepted `corr` forms.
 * @param {string[]} names variable names (matrix order)
 * @param {object | any[] | undefined} corr correlation specification, see the module header
 * @returns {number[][] | null} symmetric matrix with unit diagonal, or null when no off-diagonal term is given
 */
export function corrMatrix(names, corr) {
  if (corr === undefined || corr === null) return null;
  const n = names.length;
  const R = Array.from({ length: n }, (_, i) => {
    const row = new Array(n).fill(0);
    row[i] = 1;
    return row;
  });
  const index = new Map(names.map((k, i) => [k, i]));
  let any = false;
  const set = (a, b, rho) => {
    if (!index.has(a)) throw new RangeError(`corr: unknown variable "${a}"`);
    if (!index.has(b)) throw new RangeError(`corr: unknown variable "${b}"`);
    if (typeof rho !== 'number' || !(rho >= -1 && rho <= 1)) throw new RangeError(`corr: correlation of ${a} and ${b} must be within -1..1`);
    const i = index.get(a);
    const j = index.get(b);
    if (i === j) return;
    R[i][j] = rho;
    R[j][i] = rho;
    if (rho !== 0) any = true;
  };
  if (Array.isArray(corr)) {
    if (corr.length === n && corr.every((r) => (Array.isArray(r) || ArrayBuffer.isView(r)) && r.length === n && typeof r[0] === 'number')) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < i; j++) {
          if (Math.abs(corr[i][j] - corr[j][i]) > 1e-9) throw new RangeError('corr: matrix must be symmetric');
          set(names[i], names[j], corr[i][j]);
        }
      }
    } else {
      for (const item of corr) {
        if (Array.isArray(item)) set(item[0], item[1], item[2]);
        else if (item && typeof item === 'object') set(item.a ?? item.i ?? item.x, item.b ?? item.j ?? item.y, item.rho ?? item.value ?? item.r);
        else throw new RangeError('corr: array entries must be [a, b, rho] or { a, b, rho }');
      }
    }
  } else if (typeof corr === 'object') {
    for (const key of Object.keys(corr)) {
      const v = corr[key];
      if (v && typeof v === 'object') {
        for (const k2 of Object.keys(v)) set(key, k2, v[k2]);
      } else {
        let parts = null;
        for (const sep of [',', '|', ':', '~']) {
          const pp = key.split(sep);
          if (pp.length === 2 && index.has(pp[0].trim()) && index.has(pp[1].trim())) {
            parts = pp;
            break;
          }
        }
        if (!parts) throw new RangeError(`corr: cannot interpret key "${key}" as a pair of variable names`);
        set(parts[0].trim(), parts[1].trim(), v);
      }
    }
  } else throw new RangeError('corr: must be an object, an array of triplets or a matrix');
  return any ? R : null;
}

let glCache = null;

/**
 * Nataf model: correlation ρ0 of the underlying standard normal pair that reproduces the linear correlation ρ of
 * two variables with given marginals. Closed forms for normal/lognormal pairs, otherwise the defining integral
 * ρ = E[(X_i-μ_i)(X_j-μ_j)]/(σ_i σ_j) is evaluated by 48×48-point Gauss–Legendre quadrature on [-8, 8]² and
 * solved by Brent's method.
 * Reference: Liu & Der Kiureghian (1986), Probab. Eng. Mech. 1(2), 105–112. Validity: marginals with finite
 * variance; |ρ| must be attainable for the pair.
 * @param {object} di distribution object of the first variable
 * @param {object} dj distribution object of the second variable
 * @param {number} rho target linear correlation of the physical variables [-]
 * @returns {number} normal-space correlation ρ0 [-]
 */
export function natafRho(di, dj, rho) {
  if (rho === 0) return 0;
  if (!(di.std > 0) || !(dj.std > 0) || !Number.isFinite(di.std) || !Number.isFinite(dj.std)) return rho;
  const ni = di.name;
  const nj = dj.name;
  if (ni === 'normal' && nj === 'normal') return rho;
  const cv = (d) => d.std / d.mean;
  if (ni === 'lognormal' && nj === 'lognormal') {
    const a = cv(di);
    const b = cv(dj);
    const arg = 1 + rho * a * b;
    if (!(arg > 0)) throw new RangeError('corr: correlation not attainable for this lognormal pair');
    const r0 = Math.log(arg) / Math.sqrt(Math.log1p(a * a) * Math.log1p(b * b));
    if (!(Math.abs(r0) < 1)) throw new RangeError('corr: correlation not attainable for this lognormal pair');
    return r0;
  }
  if ((ni === 'normal' && nj === 'lognormal') || (ni === 'lognormal' && nj === 'normal')) {
    const c = cv(ni === 'lognormal' ? di : dj);
    const r0 = (rho * c) / Math.sqrt(Math.log1p(c * c));
    if (!(Math.abs(r0) < 1)) throw new RangeError('corr: correlation not attainable for this normal–lognormal pair');
    return r0;
  }
  if (!glCache) {
    const gl = gaussLegendre(48);
    const z = new Float64Array(48);
    const w = new Float64Array(48);
    for (let k = 0; k < 48; k++) {
      z[k] = 8 * gl.x[k];
      w[k] = (8 * gl.w[k] * Math.exp(-0.5 * z[k] * z[k])) / Math.sqrt(2 * Math.PI);
    }
    glCache = { z, w };
  }
  const { z, w } = glCache;
  const hi = new Float64Array(48);
  for (let k = 0; k < 48; k++) hi[k] = (di.fromZ(z[k]) - di.mean) / di.std;
  const f = (r0) => {
    const s = Math.sqrt(1 - r0 * r0);
    let acc = 0;
    for (let k = 0; k < 48; k++) {
      if (w[k] < 1e-18) continue;
      let inner = 0;
      for (let l = 0; l < 48; l++) {
        if (w[l] < 1e-18) continue;
        inner += w[l] * (dj.fromZ(r0 * z[k] + s * z[l]) - dj.mean);
      }
      acc += w[k] * hi[k] * inner;
    }
    return acc / dj.std - rho;
  };
  const lo = f(-0.9999);
  const up = f(0.9999);
  if (!(lo < 0 && up > 0)) throw new RangeError(`corr: correlation ${rho} is not attainable for a ${ni}–${nj} pair`);
  return brent(f, -0.9999, 0.9999, { tol: 1e-10, maxIter: 100 });
}

/**
 * Isoprobabilistic transformation between independent standard normal space u and physical space x.
 * Independent variables: x_i = F_i⁻¹(Φ(u_i)). Correlated variables (Nataf / Gaussian copula): z = L u with
 * L Lᵀ = R0, x_i = F_i⁻¹(Φ(z_i)), where R0 is the Nataf-adjusted normal-space correlation matrix (or the given
 * matrix itself when `nataf` is false, i.e. a Gaussian copula parameterised directly).
 * Reference: Liu & Der Kiureghian (1986); Der Kiureghian & Liu (1986), J. Eng. Mech. 112(1), 85–104.
 * @param {string[]} names variable names
 * @param {object[]} dists distribution objects
 * @param {number[][] | null} R correlation matrix of the physical variables (null: independent)
 * @param {{ nataf?: boolean }} [opts] apply the Nataf correlation correction (default true)
 * @returns {{ n: number, L: number[][] | null, R0: number[][] | null, toX: (u: ArrayLike<number>, out?: object) => object,
 *   toU: (x: object) => number[] }}
 */
export function makeTransform(names, dists, R, opts = {}) {
  const n = names.length;
  let L = null;
  let R0 = null;
  if (R) {
    R0 = R.map((r) => r.slice());
    if (opts.nataf !== false) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < i; j++) {
          if (R[i][j] !== 0) {
            const r0 = natafRho(dists[i], dists[j], R[i][j]);
            R0[i][j] = r0;
            R0[j][i] = r0;
          }
        }
      }
    }
    try {
      L = cholesky(R0);
    } catch (e) {
      throw new RangeError('corr: correlation matrix is not positive definite');
    }
  }
  const z = new Float64Array(n);
  const toX = (u, out = {}) => {
    if (L) {
      for (let i = 0; i < n; i++) {
        let s = 0;
        const Li = L[i];
        for (let j = 0; j <= i; j++) s += Li[j] * u[j];
        z[i] = s;
      }
      for (let i = 0; i < n; i++) out[names[i]] = dists[i].fromZ(z[i]);
    } else {
      for (let i = 0; i < n; i++) out[names[i]] = dists[i].fromZ(u[i]);
    }
    return out;
  };
  const toU = (x) => {
    const u = new Array(n);
    for (let i = 0; i < n; i++) u[i] = dists[i].std > 0 ? dists[i].toZ(x[names[i]]) : 0;
    if (L) {
      for (let i = 0; i < n; i++) {
        let s = u[i];
        for (let j = 0; j < i; j++) s -= L[i][j] * u[j];
        u[i] = s / L[i][i];
      }
    }
    return u;
  };
  return { n, L, R0, toX, toU };
}

function checkN(n, name) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError(`${name}: n must be a positive integer`);
}

/**
 * Crude Monte Carlo simulation of a limit-state (or any scalar) function.
 * Independent inputs are drawn from their marginals; with `corr` the Nataf / Gaussian-copula transformation is
 * used. pf = P[g < 0] = N_f / n with the 95 % Wilson score interval.
 * Reference: Rubinstein & Kroese, "Simulation and the Monte Carlo Method"; Wilson (1927), JASA 22, 209–212.
 * Validity: pf resolvable only when n·pf >~ 10.
 * The object passed to `g` is reused between calls — copy it if it must be kept.
 * @param {(x: Object<string, number>) => number} g limit state / response function (failure when g < 0)
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ n?: number, seed?: number, corr?: object, keepInputs?: boolean, nataf?: boolean }} [opts]
 *   n samples (10000); seed (1); corr correlations; keepInputs: return the sampled inputs (default n <= 200000)
 * @returns {{ samples: Float64Array, inputs?: Object<string, Float64Array>, mean: number, std: number, pf: number,
 *   ciPf: [number, number], covPf: number, nFail: number, n: number, beta: number }}
 *   g values; inputs per variable; mean and standard deviation of g; failure probability; its 95 % interval;
 *   coefficient of variation of the pf estimator; failures counted; generalised reliability index -Φ⁻¹(pf)
 */
export function monteCarlo(g, vars, opts = {}) {
  if (typeof g !== 'function') throw new RangeError('monteCarlo: g must be a function');
  const { names, dists } = resolveVars(vars);
  const n = opts.n ?? 10000;
  checkN(n, 'monteCarlo');
  const rng = makeRng(opts.seed ?? 1);
  const d = names.length;
  const R = corrMatrix(names, opts.corr);
  const tr = R ? makeTransform(names, dists, R, { nataf: opts.nataf }) : null;
  const keep = opts.keepInputs ?? n <= 200000;
  const inputs = keep ? Object.fromEntries(names.map((k) => [k, new Float64Array(n)])) : undefined;
  const cols = keep ? names.map((k) => inputs[k]) : null;
  const samples = new Float64Array(n);
  const x = {};
  const u = new Float64Array(d);
  let s1 = 0;
  let nFail = 0;
  for (let k = 0; k < n; k++) {
    if (tr) {
      for (let i = 0; i < d; i++) u[i] = rng.normal();
      tr.toX(u, x);
    } else {
      for (let i = 0; i < d; i++) x[names[i]] = dists[i].sample(rng);
    }
    if (keep) for (let i = 0; i < d; i++) cols[i][k] = x[names[i]];
    const v = g(x);
    samples[k] = v;
    s1 += v;
    if (v < 0) nFail++;
  }
  const mean = s1 / n;
  let s2 = 0;
  for (let k = 0; k < n; k++) s2 += (samples[k] - mean) * (samples[k] - mean);
  const std = n > 1 ? Math.sqrt(s2 / (n - 1)) : 0;
  const pf = nFail / n;
  const zc = 1.959963984540054;
  const den = 1 + (zc * zc) / n;
  const centre = (pf + (zc * zc) / (2 * n)) / den;
  const halfW = (zc * Math.sqrt((pf * (1 - pf)) / n + (zc * zc) / (4 * n * n))) / den;
  const out = {
    samples,
    mean,
    std,
    pf,
    ciPf: [Math.max(0, centre - halfW), Math.min(1, centre + halfW)],
    covPf: nFail > 0 ? Math.sqrt((1 - pf) / (n * pf)) : Infinity,
    nFail,
    n,
    beta: pf > 0 && pf < 1 ? -normInv(pf) : pf === 0 ? Infinity : -Infinity,
  };
  if (keep) out.inputs = inputs;
  return out;
}

/** Solve L·X = B for lower-triangular L where B is given by rows (in place on a copy). */
function lowerInverse(L) {
  const n = L.length;
  const inv = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let c = 0; c < n; c++) {
    for (let i = c; i < n; i++) {
      let s = i === c ? 1 : 0;
      for (let j = c; j < i; j++) s -= L[i][j] * inv[j][c];
      inv[i][c] = s / L[i][i];
    }
  }
  return inv;
}

function argsort(a) {
  const idx = new Uint32Array(a.length);
  for (let i = 0; i < idx.length; i++) idx[i] = i;
  idx.sort((i, j) => a[i] - a[j]);
  return idx;
}

/**
 * Latin hypercube sample: each marginal is divided into n equiprobable strata, one value is drawn uniformly
 * inside each stratum and mapped through the inverse cdf; the strata are randomly paired between variables.
 * With `corr` the pairing is rearranged by the Iman–Conover method so that the sample rank (Spearman)
 * correlation approximates the target while every marginal stays exactly stratified.
 * Reference: McKay, Beckman & Conover (1979), Technometrics 21, 239–245; Iman & Conover (1982), Commun. Stat.
 * Simul. Comput. 11, 311–334. Validity: n > number of variables for the correlated variant.
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {number} n sample size
 * @param {{ seed?: number, corr?: object }} [opts] seed (1); target rank correlations
 * @returns {Object<string, Float64Array>} one column of n values per variable
 */
export function latinHypercube(vars, n, opts = {}) {
  const { names, dists } = resolveVars(vars);
  checkN(n, 'latinHypercube');
  const rng = makeRng(opts.seed ?? 1);
  const d = names.length;
  const R = corrMatrix(names, opts.corr);
  const cols = new Array(d);
  const perm = new Uint32Array(n);
  for (let i = 0; i < d; i++) {
    for (let k = 0; k < n; k++) perm[k] = k;
    rng.shuffle(perm);
    const col = new Float64Array(n);
    const di = dists[i];
    for (let k = 0; k < n; k++) {
      let p = (perm[k] + rng.next()) / n;
      if (p <= 0) p = 0.5 / n;
      col[k] = di.icdf(p);
    }
    cols[i] = col;
  }
  if (R) {
    if (n <= d) throw new RangeError('latinHypercube: n must exceed the number of variables when corr is given');
    // van der Waerden scores, independently permuted per column
    const M = new Array(d);
    const base = new Float64Array(n);
    for (let k = 0; k < n; k++) base[k] = normInv((k + 1) / (n + 1));
    for (let i = 0; i < d; i++) {
      const c = Float64Array.from(base);
      rng.shuffle(c);
      M[i] = c;
    }
    // sample correlation T of the scores (columns have zero mean and equal variance)
    let ss = 0;
    for (let k = 0; k < n; k++) ss += base[k] * base[k];
    const T = Array.from({ length: d }, () => new Array(d).fill(0));
    for (let i = 0; i < d; i++) {
      for (let j = 0; j <= i; j++) {
        let s = 0;
        for (let k = 0; k < n; k++) s += M[i][k] * M[j][k];
        T[i][j] = s / ss;
        T[j][i] = T[i][j];
      }
    }
    let P;
    let Q;
    try {
      P = cholesky(R);
      Q = cholesky(T);
    } catch (e) {
      throw new RangeError('latinHypercube: corr is not positive definite');
    }
    const Qi = lowerInverse(Q);
    // S = P Q⁻¹ ; M* = M Sᵀ  (row vector form: m*_i = Σ_j S_ij m_j)
    const S = Array.from({ length: d }, () => new Array(d).fill(0));
    for (let i = 0; i < d; i++) for (let j = 0; j < d; j++) for (let l = 0; l < d; l++) S[i][j] += P[i][l] * Qi[l][j];
    const target = new Float64Array(n);
    for (let i = 0; i < d; i++) {
      target.fill(0);
      for (let j = 0; j < d; j++) {
        const sij = S[i][j];
        if (sij === 0) continue;
        const Mj = M[j];
        for (let k = 0; k < n; k++) target[k] += sij * Mj[k];
      }
      const order = argsort(target);
      const sorted = Float64Array.from(cols[i]).sort();
      for (let r = 0; r < n; r++) cols[i][order[r]] = sorted[r];
    }
  }
  const out = {};
  for (let i = 0; i < d; i++) out[names[i]] = cols[i];
  return out;
}

function summarise(arr) {
  const n = arr.length;
  const s = Float64Array.from(arr).sort();
  let m = 0;
  for (let i = 0; i < n; i++) m += s[i];
  m /= n;
  let v = 0;
  for (let i = 0; i < n; i++) v += (s[i] - m) * (s[i] - m);
  return {
    mean: m,
    std: n > 1 ? Math.sqrt(v / (n - 1)) : 0,
    p10: percentileSorted(s, 10),
    p50: percentileSorted(s, 50),
    p90: percentileSorted(s, 90),
    min: s[0],
    max: s[n - 1],
    p05: percentileSorted(s, 5),
    p95: percentileSorted(s, 95),
  };
}

/**
 * Propagate input uncertainty through a model with several named outputs.
 * Reference: McKay et al. (1979) for the default Latin hypercube design. Validity: the model must return the
 * same output names for every sample (a bare number is stored under the name `value`).
 * @param {(x: Object<string, number>) => Object<string, number> | number} model model evaluated per sample
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ n?: number, seed?: number, method?: 'lhs' | 'mc', corr?: object }} [opts] n (1000); seed (1);
 *   sampling method (Latin hypercube by default); optional correlations (rank correlation for 'lhs', Nataf for 'mc')
 * @returns {{ outputs: Object<string, Float64Array>, inputs: Object<string, Float64Array>,
 *   stats: Object<string, { mean: number, std: number, p10: number, p50: number, p90: number, min: number,
 *   max: number, p05: number, p95: number }> }}
 */
export function propagate(model, vars, opts = {}) {
  if (typeof model !== 'function') throw new RangeError('propagate: model must be a function');
  const n = opts.n ?? 1000;
  checkN(n, 'propagate');
  const method = opts.method ?? 'lhs';
  const { names, dists } = resolveVars(vars);
  let inputs;
  if (method === 'lhs') inputs = latinHypercube(vars, n, { seed: opts.seed ?? 1, corr: opts.corr });
  else if (method === 'mc') {
    const rng = makeRng(opts.seed ?? 1);
    const d = names.length;
    const R = corrMatrix(names, opts.corr);
    const tr = R ? makeTransform(names, dists, R) : null;
    inputs = Object.fromEntries(names.map((k) => [k, new Float64Array(n)]));
    const u = new Float64Array(d);
    const x = {};
    for (let k = 0; k < n; k++) {
      if (tr) {
        for (let i = 0; i < d; i++) u[i] = rng.normal();
        tr.toX(u, x);
        for (let i = 0; i < d; i++) inputs[names[i]][k] = x[names[i]];
      } else {
        for (let i = 0; i < d; i++) inputs[names[i]][k] = dists[i].sample(rng);
      }
    }
  } else throw new RangeError("propagate: method must be 'lhs' or 'mc'");
  const outputs = {};
  let outNames = null;
  for (let k = 0; k < n; k++) {
    const x = {};
    for (let i = 0; i < names.length; i++) x[names[i]] = inputs[names[i]][k];
    let r = model(x);
    if (typeof r === 'number') r = { value: r };
    if (!r || typeof r !== 'object') throw new RangeError('propagate: model must return an object of named outputs');
    if (!outNames) {
      outNames = Object.keys(r).filter((o) => typeof r[o] === 'number');
      for (const o of outNames) outputs[o] = new Float64Array(n);
    }
    for (const o of outNames) outputs[o][k] = r[o];
  }
  const stats = {};
  for (const o of outNames) stats[o] = summarise(outputs[o]);
  return { outputs, inputs, stats };
}

function firstPrimes(count) {
  const p = [];
  for (let c = 2; p.length < count; c++) {
    let ok = true;
    for (let i = 0; i < p.length && p[i] * p[i] <= c; i++) {
      if (c % p[i] === 0) {
        ok = false;
        break;
      }
    }
    if (ok) p.push(c);
  }
  return p;
}

/**
 * Randomised Halton low-discrepancy sequence (independent random digit permutations per digit position and
 * dimension, which keeps every one-dimensional projection uniformly distributed and removes the correlation
 * between high-dimension pairs of the plain sequence).
 * Reference: Halton (1960), Numer. Math. 2, 84–90; Owen (2017), "A randomized Halton algorithm in R",
 * arXiv:1706.02808.
 * @param {number} n number of points
 * @param {number} dim number of dimensions
 * @param {{ next: () => number, shuffle: (a: any) => any }} rng seeded generator
 * @returns {Float64Array[]} dim columns of n values in (0, 1)
 */
export function scrambledHalton(n, dim, rng) {
  const primes = firstPrimes(dim);
  const cols = new Array(dim);
  for (let j = 0; j < dim; j++) {
    const b = primes[j];
    const K = Math.ceil(Math.log(4294967296 * 1024) / Math.log(b));
    const perms = new Array(K);
    for (let k = 0; k < K; k++) {
      const p = new Uint16Array(b);
      for (let i = 0; i < b; i++) p[i] = i;
      rng.shuffle(p);
      perms[k] = p;
    }
    const col = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let idx = i;
      let f = 1 / b;
      let v = 0;
      for (let k = 0; k < K; k++) {
        const digit = idx % b;
        idx = (idx - digit) / b;
        v += perms[k][digit] * f;
        f /= b;
      }
      col[i] = Math.min(Math.max(v, 1e-12), 1 - 1e-12);
    }
    cols[j] = col;
  }
  return cols;
}

/**
 * Variance-based global sensitivity indices (Sobol) for independent inputs by the Saltelli sampling scheme:
 * two base matrices A and B (n rows) and the d hybrid matrices AB_i (A with column i taken from B), i.e.
 * n·(d + 2) model evaluations. Estimators: first order S_i = mean(f_B·(f_ABi - f_A)) / V (Saltelli et al. 2010,
 * Table 2 (b)); total effect ST_i = mean((f_A - f_ABi)²) / (2V) (Jansen 1999). The base sample is a randomised
 * Halton sequence in 2d dimensions.
 * Reference: Sobol' (1993); Saltelli et al. (2010), Comput. Phys. Commun. 181, 259–270; Jansen (1999), Comput.
 * Phys. Commun. 117, 35–43. Validity: independent inputs, finite output variance; sampling error ~ O(1/√n) or
 * better (indices of non-influential inputs scatter around 0 and may be slightly negative).
 * @param {(x: Object<string, number>) => number} model scalar model (the argument object is reused)
 * @param {Object<string, object>} vars { name: distSpec }
 * @param {{ n?: number, seed?: number }} [opts] base sample size n (1024); seed (1)
 * @returns {{ first: Object<string, number>, total: Object<string, number>, variance: number, mean: number,
 *   evaluations: number }} first-order and total-effect indices [-], output variance and mean, model calls
 */
export function sobolIndices(model, vars, opts = {}) {
  if (typeof model !== 'function') throw new RangeError('sobolIndices: model must be a function');
  const { names, dists } = resolveVars(vars);
  const n = opts.n ?? 1024;
  checkN(n, 'sobolIndices');
  if (n < 8) throw new RangeError('sobolIndices: n must be at least 8');
  const d = names.length;
  const rng = makeRng(opts.seed ?? 1);
  const U = scrambledHalton(n, 2 * d, rng);
  const A = new Array(d);
  const B = new Array(d);
  for (let i = 0; i < d; i++) {
    A[i] = new Float64Array(n);
    B[i] = new Float64Array(n);
    const di = dists[i];
    for (let k = 0; k < n; k++) {
      A[i][k] = di.icdf(U[i][k]);
      B[i][k] = di.icdf(U[d + i][k]);
    }
  }
  const x = {};
  const evalWith = (k, swap) => {
    for (let i = 0; i < d; i++) x[names[i]] = i === swap ? B[i][k] : A[i][k];
    return model(x);
  };
  const fA = new Float64Array(n);
  const fB = new Float64Array(n);
  let s = 0;
  for (let k = 0; k < n; k++) {
    fA[k] = evalWith(k, -1);
    for (let i = 0; i < d; i++) x[names[i]] = B[i][k];
    fB[k] = model(x);
    s += fA[k] + fB[k];
  }
  const f0 = s / (2 * n);
  let V = 0;
  for (let k = 0; k < n; k++) V += (fA[k] - f0) * (fA[k] - f0) + (fB[k] - f0) * (fB[k] - f0);
  V /= 2 * n - 1;
  const first = {};
  const total = {};
  for (let i = 0; i < d; i++) {
    let sf = 0;
    let st = 0;
    if (dists[i].std > 0) {
      for (let k = 0; k < n; k++) {
        const fab = evalWith(k, i);
        sf += fB[k] * (fab - fA[k]);
        st += (fA[k] - fab) * (fA[k] - fab);
      }
    }
    first[names[i]] = V > 0 ? sf / n / V : 0;
    total[names[i]] = V > 0 ? st / (2 * n) / V : 0;
  }
  return { first, total, variance: V, mean: f0, evaluations: n * (d + 2) };
}

/**
 * One-at-a-time (tornado) sensitivity: every variable is moved to its low and high value in turn while the
 * others stay at their base value.
 * Reference: standard deterministic sensitivity practice, e.g. Clemen, "Making Hard Decisions", ch. 5.
 * @param {(x: Object<string, number>) => number} model scalar model
 * @param {Object<string, { low: number, base: number, high: number }>} vars ranges per variable
 * @returns {{ name: string, low: number, high: number, swing: number, base: number, xLow: number, xHigh: number }[]}
 *   model output at the low and high input value, swing = |high - low|, base-case output and the input values
 *   used; sorted by decreasing swing
 */
export function tornado(model, vars) {
  if (typeof model !== 'function') throw new RangeError('tornado: model must be a function');
  if (!vars || typeof vars !== 'object') throw new RangeError('tornado: vars must be an object');
  const names = Object.keys(vars);
  const baseX = {};
  for (const k of names) {
    const v = vars[k];
    if (!v || typeof v.low !== 'number' || typeof v.high !== 'number') throw new RangeError(`tornado: vars.${k} needs numeric low and high`);
    baseX[k] = typeof v.base === 'number' ? v.base : 0.5 * (v.low + v.high);
  }
  const one = (x) => {
    const r = model(x);
    return typeof r === 'number' ? r : NaN;
  };
  const base = one({ ...baseX });
  const out = names.map((k) => {
    const low = one({ ...baseX, [k]: vars[k].low });
    const high = one({ ...baseX, [k]: vars[k].high });
    return { name: k, low, high, swing: Math.abs(high - low), base, xLow: vars[k].low, xHigh: vars[k].high };
  });
  out.sort((a, b) => b.swing - a.swing);
  return out;
}

export { normCdf, normInv };
