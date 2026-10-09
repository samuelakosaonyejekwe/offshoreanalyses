/**
 * Surrogate modelling and Bayesian tools: Gaussian-process regression (kriging), random-walk Metropolis MCMC and
 * Bayesian optimisation with the expected-improvement acquisition function.
 */
import { makeRng } from '../math/rng.js';
import { normCdf, normPdf } from '../math/stats.js';
import { cholesky } from '../math/linalg.js';
import { nelderMead } from '../math/optim.js';

const SQRT5 = Math.sqrt(5);

function toRows(X, name) {
  if (!X || typeof X.length !== 'number' || X.length === 0) throw new RangeError(`${name}: X must be a non-empty array of rows`);
  return Array.from(X, (r, i) => {
    const row = typeof r === 'number' ? [r] : Array.from(r);
    for (const v of row) if (!Number.isFinite(v)) throw new RangeError(`${name}: X[${i}] contains a non-finite value`);
    return row;
  });
}

/**
 * Gaussian-process regression with a constant (sample-mean) trend and an anisotropic stationary kernel.
 *
 * Inputs are scaled to the unit cube and the response is standardised. Kernels (r² = Σ ((x_d - x'_d)/ℓ_d)²):
 * 'rbf' k = exp(-r²/2); 'matern52' k = (1 + √5 r + 5r²/3) exp(-√5 r). The covariance is σ²(K + g I) with nugget
 * g = noise / var(y); σ² is profiled analytically (σ̂² = yᵀ(K+gI)⁻¹y/n) and the length scales maximise the
 * concentrated log marginal likelihood (bounded Nelder–Mead on ln ℓ, two starts; one length scale per input for
 * d <= 8, a common one otherwise). All solves use the Cholesky factorisation.
 * Prediction: mean = m + k*ᵀ(K+gI)⁻¹(y - m), variance of the latent function σ²(1 - k*ᵀ(K+gI)⁻¹k*).
 * Reference: Rasmussen & Williams (2006), "Gaussian Processes for Machine Learning", Alg. 2.1 and §5.4.1;
 * Jones, Schonlau & Welch (1998), J. Global Optim. 13, 455–492 (concentrated likelihood).
 * Validity: n >= 2 distinct points, n <= ~500 (O(n³) per likelihood evaluation), smooth responses.
 *
 * @param {ArrayLike<ArrayLike<number> | number>} X training inputs, n rows of d values (numbers allowed for d = 1)
 * @param {ArrayLike<number>} y training responses (n) [uy]
 * @param {{ kernel?: 'rbf' | 'matern52', lengthScale?: number | number[], noise?: number, optimise?: boolean }} [opts]
 *   kernel ('rbf'); lengthScale in the units of X (start value when optimising, fixed value otherwise; default
 *   0.3 × the range of each input); noise: observation noise variance [uy²] (1e-8); optimise length scales (true)
 * @returns {{ predict: (x: number | number[] | number[][]) => { mean: number, std: number } | { mean: number, std: number }[],
 *   logLik: number, params: { kernel: string, lengthScale: number[], signalVariance: number, noise: number,
 *   mean: number }, n: number, dim: number }} predictor (a single point gives one result, an array of points an
 *   array of results; std is the posterior standard deviation of the latent function [uy]); log marginal
 *   likelihood of y [-]; fitted hyper-parameters in physical units
 */
export function gaussianProcess(X, y, opts = {}) {
  const rows = toRows(X, 'gaussianProcess');
  const n = rows.length;
  const d = rows[0].length;
  if (!y || y.length !== n) throw new RangeError('gaussianProcess: y must have one value per row of X');
  if (n < 2) throw new RangeError('gaussianProcess: at least two training points are required');
  for (const r of rows) if (r.length !== d) throw new RangeError('gaussianProcess: all rows of X must have the same length');
  const kernel = opts.kernel ?? 'rbf';
  if (kernel !== 'rbf' && kernel !== 'matern52') throw new RangeError("gaussianProcess: kernel must be 'rbf' or 'matern52'");
  const noise = opts.noise ?? 1e-8;
  if (!(noise >= 0)) throw new RangeError('gaussianProcess: noise must be >= 0');
  const optimise = opts.optimise ?? true;

  // scaling
  const xlo = new Float64Array(d).fill(Infinity);
  const xr = new Float64Array(d);
  for (let j = 0; j < d; j++) {
    let hi = -Infinity;
    for (let i = 0; i < n; i++) {
      if (rows[i][j] < xlo[j]) xlo[j] = rows[i][j];
      if (rows[i][j] > hi) hi = rows[i][j];
    }
    xr[j] = hi > xlo[j] ? hi - xlo[j] : 1;
  }
  const Z = new Float64Array(n * d);
  for (let i = 0; i < n; i++) for (let j = 0; j < d; j++) Z[i * d + j] = (rows[i][j] - xlo[j]) / xr[j];
  let ym = 0;
  for (let i = 0; i < n; i++) {
    if (!Number.isFinite(y[i])) throw new RangeError(`gaussianProcess: y[${i}] is not finite`);
    ym += y[i];
  }
  ym /= n;
  let ys = 0;
  for (let i = 0; i < n; i++) ys += (y[i] - ym) * (y[i] - ym);
  ys = Math.sqrt(ys / (n - 1));
  if (!(ys > 0)) ys = 1;
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) t[i] = (y[i] - ym) / ys;
  const nugget = Math.max(noise / (ys * ys), 1e-10);

  const kfun = (r2) => {
    if (kernel === 'rbf') return Math.exp(-0.5 * r2);
    const r = Math.sqrt(r2);
    return (1 + SQRT5 * r + (5 / 3) * r2) * Math.exp(-SQRT5 * r);
  };
  const ard = d <= 8;
  const expand = (p) => {
    const inv = new Float64Array(d);
    for (let j = 0; j < d; j++) inv[j] = Math.exp(-(ard ? p[j] : p[0]));
    return inv;
  };

  const factor = (invL) => {
    let jitter = nugget;
    for (let attempt = 0; attempt < 6; attempt++) {
      const K = new Array(n);
      for (let i = 0; i < n; i++) K[i] = new Array(n);
      for (let i = 0; i < n; i++) {
        K[i][i] = 1 + jitter;
        for (let k = 0; k < i; k++) {
          let r2 = 0;
          for (let j = 0; j < d; j++) {
            const q = (Z[i * d + j] - Z[k * d + j]) * invL[j];
            r2 += q * q;
          }
          const v = kfun(r2);
          K[i][k] = v;
          K[k][i] = v;
        }
      }
      try {
        const L = cholesky(K);
        // α = K⁻¹ t by forward and back substitution
        const w = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          let s = t[i];
          const Li = L[i];
          for (let k = 0; k < i; k++) s -= Li[k] * w[k];
          w[i] = s / Li[i];
        }
        let quad = 0;
        let logDet = 0;
        for (let i = 0; i < n; i++) {
          quad += w[i] * w[i];
          logDet += Math.log(L[i][i]);
        }
        const alpha = new Float64Array(n);
        for (let i = n - 1; i >= 0; i--) {
          let s = w[i];
          for (let k = i + 1; k < n; k++) s -= L[k][i] * alpha[k];
          alpha[i] = s / L[i][i];
        }
        const s2 = quad / n;
        return { L, alpha, s2, nll: 0.5 * n * Math.log(s2) + logDet, jitter, logDet };
      } catch (e) {
        jitter *= 100;
      }
    }
    return null;
  };

  const nPar = ard ? d : 1;
  let p0 = new Array(nPar).fill(Math.log(0.3));
  if (opts.lengthScale !== undefined) {
    const ls = typeof opts.lengthScale === 'number' ? new Array(d).fill(opts.lengthScale) : Array.from(opts.lengthScale);
    if (ls.length !== d) throw new RangeError('gaussianProcess: lengthScale must be a number or one value per input');
    for (const v of ls) if (!(v > 0)) throw new RangeError('gaussianProcess: lengthScale must be > 0');
    if (ard) p0 = ls.map((v, j) => Math.log(v / xr[j]));
    else p0 = [ls.reduce((s, v, j) => s + Math.log(v / xr[j]), 0) / d];
  }
  let pBest = p0;
  if (optimise) {
    const obj = (p) => {
      const f = factor(expand(p));
      return f && Number.isFinite(f.nll) ? f.nll : 1e300;
    };
    const bounds = new Array(nPar).fill(0).map(() => [Math.log(0.01), Math.log(30)]);
    const clampP = (p) => p.map((v) => Math.min(Math.max(v, bounds[0][0]), bounds[0][1]));
    const starts = [clampP(p0)];
    starts.push(clampP(p0.map((v) => v + Math.log(4))));
    let best = null;
    for (const st of starts) {
      const r = nelderMead(obj, st, { bounds, tol: 1e-7, maxIter: 60 + 60 * nPar, step: 0.7 });
      if (!best || r.f < best.f) best = r;
    }
    pBest = best.x;
  }
  const invL = expand(pBest);
  const fit = factor(invL);
  if (!fit) throw new Error('gaussianProcess: covariance matrix is not positive definite (duplicate points? increase noise)');
  const { L, alpha, s2 } = fit;
  const logLik = -0.5 * n * Math.log(2 * Math.PI * s2) - fit.logDet - 0.5 * n - n * Math.log(ys);

  const kv = new Float64Array(n);
  const one = (x) => {
    const xv = typeof x === 'number' ? [x] : x;
    if (!xv || xv.length !== d) throw new RangeError(`gaussianProcess.predict: x must have ${d} value(s)`);
    let m = 0;
    for (let i = 0; i < n; i++) {
      let r2 = 0;
      for (let j = 0; j < d; j++) {
        const q = ((xv[j] - xlo[j]) / xr[j] - Z[i * d + j]) * invL[j];
        r2 += q * q;
      }
      kv[i] = kfun(r2);
      m += kv[i] * alpha[i];
    }
    // v = L⁻¹ k*
    let vv = 0;
    for (let i = 0; i < n; i++) {
      let s = kv[i];
      const Li = L[i];
      for (let k = 0; k < i; k++) s -= Li[k] * kv[k];
      kv[i] = s / Li[i];
      vv += kv[i] * kv[i];
    }
    const variance = s2 * Math.max(1 - vv, 0);
    return { mean: ym + ys * m, std: ys * Math.sqrt(variance) };
  };
  const predict = (x) => {
    if (Array.isArray(x) && x.length > 0 && typeof x[0] !== 'number') return x.map(one);
    if (d === 1 && Array.isArray(x) && x.length !== 1) return x.map(one);
    return one(x);
  };
  return {
    predict,
    logLik,
    params: {
      kernel,
      lengthScale: Array.from(invL, (v, j) => xr[j] / v),
      signalVariance: s2 * ys * ys,
      noise: fit.jitter * s2 * ys * ys,
      mean: ym,
    },
    n,
    dim: d,
  };
}

/**
 * Random-walk Metropolis sampler with independent Gaussian proposals per coordinate.
 * Reference: Metropolis et al. (1953), J. Chem. Phys. 21, 1087–1092; Gelman et al., "Bayesian Data Analysis",
 * 3rd ed., §11.2 (target acceptance 0.2–0.5). Validity: logPost finite at x0; return -Infinity outside the
 * support.
 * @param {(x: number[]) => number} logPost log of the (unnormalised) posterior density
 * @param {ArrayLike<number>} x0 starting point
 * @param {{ n?: number, stepSizes?: number | number[], seed?: number, burn?: number, thin?: number }} [opts]
 *   n retained samples (5000); proposal standard deviations (default 0.1·max(|x0|, 1)); seed (1); burn-in
 *   iterations discarded (n/5); thinning interval (1)
 * @returns {{ chain: number[][], acceptance: number, logPost: number[], mean: number[], std: number[] }}
 *   retained samples, acceptance rate over all iterations [-], log-posterior of each retained sample, posterior
 *   mean and standard deviation per coordinate
 */
export function metropolis(logPost, x0, opts = {}) {
  if (typeof logPost !== 'function') throw new RangeError('metropolis: logPost must be a function');
  const d = x0.length;
  if (!(d >= 1)) throw new RangeError('metropolis: x0 must have at least one element');
  const n = opts.n ?? 5000;
  if (!Number.isInteger(n) || n < 1) throw new RangeError('metropolis: n must be a positive integer');
  const burn = opts.burn ?? Math.floor(n / 5);
  const thin = opts.thin ?? 1;
  if (!Number.isInteger(burn) || burn < 0) throw new RangeError('metropolis: burn must be a non-negative integer');
  if (!Number.isInteger(thin) || thin < 1) throw new RangeError('metropolis: thin must be a positive integer');
  let steps;
  if (opts.stepSizes === undefined) steps = Array.from(x0, (v) => 0.1 * Math.max(Math.abs(v), 1));
  else if (typeof opts.stepSizes === 'number') steps = new Array(d).fill(opts.stepSizes);
  else steps = Array.from(opts.stepSizes);
  if (steps.length !== d) throw new RangeError('metropolis: stepSizes must have one entry per coordinate');
  for (const s of steps) if (!(s > 0)) throw new RangeError('metropolis: stepSizes must be > 0');
  const rng = makeRng(opts.seed ?? 1);
  let x = Array.from(x0);
  let lp = logPost(x);
  if (!Number.isFinite(lp)) throw new RangeError('metropolis: logPost(x0) must be finite');
  const chain = new Array(n);
  const lps = new Array(n);
  const cand = new Array(d);
  let accepted = 0;
  const total = burn + n * thin;
  let kept = 0;
  for (let it = 0; it < total; it++) {
    for (let j = 0; j < d; j++) cand[j] = x[j] + steps[j] * rng.normal();
    const lc = logPost(cand);
    if (lc >= lp || Math.log(1 - rng.next()) < lc - lp) {
      x = cand.slice();
      lp = lc;
      accepted++;
    }
    if (it >= burn && (it - burn) % thin === thin - 1) {
      chain[kept] = x.slice();
      lps[kept] = lp;
      kept++;
    }
  }
  const mean = new Array(d).fill(0);
  const std = new Array(d).fill(0);
  for (let k = 0; k < n; k++) for (let j = 0; j < d; j++) mean[j] += chain[k][j] / n;
  for (let k = 0; k < n; k++) for (let j = 0; j < d; j++) std[j] += (chain[k][j] - mean[j]) ** 2;
  for (let j = 0; j < d; j++) std[j] = n > 1 ? Math.sqrt(std[j] / (n - 1)) : 0;
  return { chain, acceptance: accepted / total, logPost: lps, mean, std };
}

/**
 * Expected improvement (for minimisation) of a Gaussian prediction over the incumbent best value:
 * EI = (f_min - μ) Φ(z) + s φ(z), z = (f_min - μ)/s.
 * Reference: Jones, Schonlau & Welch (1998), eq. 15.
 * @param {number} mean predicted mean [uy]
 * @param {number} std predicted standard deviation [uy]
 * @param {number} fBest incumbent minimum [uy]
 * @returns {number} expected improvement [uy]
 */
export function expectedImprovement(mean, std, fBest) {
  if (!(std > 0)) return Math.max(fBest - mean, 0);
  const z = (fBest - mean) / std;
  return (fBest - mean) * normCdf(z) + std * normPdf(z);
}

/**
 * Bayesian (efficient global) optimisation: a Gaussian-process surrogate (Matérn 5/2) is fitted to all
 * evaluations and the next point maximises the expected improvement, searched by seeded random candidates plus
 * local candidates around the incumbent and polished by bounded Nelder–Mead.
 * Reference: Jones, Schonlau & Welch (1998), J. Global Optim. 13, 455–492; Mockus (1975).
 * Validity: expensive, smooth, low-dimensional objectives (d <= ~10); finite bounds; minimisation.
 * @param {(x: number[]) => number} f objective to minimise
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ nInit?: number, nIter?: number, seed?: number, kernel?: 'rbf' | 'matern52' }} [opts] size of the
 *   initial Latin hypercube design (max(5, 2d + 1)); number of adaptive evaluations (25); seed (1)
 * @returns {{ x: number[], f: number, history: number[], X: number[][], y: number[] }} best point and value,
 *   best-so-far objective after each evaluation, all evaluated points and their values
 */
export function bayesOptimise(f, bounds, opts = {}) {
  if (typeof f !== 'function') throw new RangeError('bayesOptimise: f must be a function');
  if (!Array.isArray(bounds) || bounds.length === 0) throw new RangeError('bayesOptimise: bounds are required');
  const d = bounds.length;
  for (let j = 0; j < d; j++) {
    const b = bounds[j];
    if (!b || !Number.isFinite(b[0]) || !Number.isFinite(b[1]) || !(b[1] > b[0])) throw new RangeError(`bayesOptimise: bounds[${j}] must be finite with lo < hi`);
  }
  const nInit = opts.nInit ?? Math.max(5, 2 * d + 1);
  const nIter = opts.nIter ?? 25;
  if (!Number.isInteger(nInit) || nInit < 2) throw new RangeError('bayesOptimise: nInit must be an integer >= 2');
  if (!Number.isInteger(nIter) || nIter < 0) throw new RangeError('bayesOptimise: nIter must be a non-negative integer');
  const kernel = opts.kernel ?? 'matern52';
  const rng = makeRng(opts.seed ?? 1);
  const X = [];
  const y = [];
  const history = [];
  let fBest = Infinity;
  let xBest = null;
  const evaluate = (x) => {
    const v = f(x.slice());
    if (!Number.isFinite(v)) throw new Error('bayesOptimise: objective returned a non-finite value');
    X.push(x.slice());
    y.push(v);
    if (v < fBest) {
      fBest = v;
      xBest = x.slice();
    }
    history.push(fBest);
  };
  // Latin hypercube initial design
  const perms = [];
  for (let j = 0; j < d; j++) perms.push(rng.shuffle(Array.from({ length: nInit }, (_, i) => i)));
  for (let i = 0; i < nInit; i++) {
    const x = new Array(d);
    for (let j = 0; j < d; j++) x[j] = bounds[j][0] + ((perms[j][i] + rng.next()) / nInit) * (bounds[j][1] - bounds[j][0]);
    evaluate(x);
  }
  const span = bounds.map((b) => b[1] - b[0]);
  const tooClose = (x) => {
    for (const p of X) {
      let m = 0;
      for (let j = 0; j < d; j++) m = Math.max(m, Math.abs(p[j] - x[j]) / span[j]);
      if (m < 1e-7) return true;
    }
    return false;
  };
  for (let it = 0; it < nIter; it++) {
    let ySpan = 0;
    for (const v of y) ySpan = Math.max(ySpan, Math.abs(v - fBest));
    let next = null;
    if (ySpan > 0) {
      const gp = gaussianProcess(X, y, { kernel, noise: 1e-10 * ySpan * ySpan });
      const negEi = (x) => {
        const p = gp.predict(x);
        return -expectedImprovement(p.mean, p.std, fBest);
      };
      const nCand = Math.min(3000, 400 * d);
      const cands = [];
      for (let c = 0; c < nCand; c++) {
        const x = new Array(d);
        if (c % 4 === 3) {
          // local candidate around the incumbent
          const sc = c % 8 === 3 ? 0.02 : 0.1;
          for (let j = 0; j < d; j++) x[j] = Math.min(Math.max(xBest[j] + sc * span[j] * rng.normal(), bounds[j][0]), bounds[j][1]);
        } else {
          for (let j = 0; j < d; j++) x[j] = bounds[j][0] + span[j] * rng.next();
        }
        cands.push({ x, v: negEi(x) });
      }
      cands.sort((a, b) => a.v - b.v);
      let best = cands[0];
      for (let c = 0; c < Math.min(3, cands.length); c++) {
        const r = nelderMead(negEi, cands[c].x, { bounds, tol: 1e-9, maxIter: 40 * d + 20, step: span.map((s) => 0.02 * s) });
        if (r.f < best.v) best = { x: r.x, v: r.f };
      }
      if (!tooClose(best.x)) next = best.x;
    }
    if (!next) {
      // surrogate gave no new information: explore at random
      next = new Array(d);
      for (let j = 0; j < d; j++) next[j] = bounds[j][0] + span[j] * rng.next();
    }
    evaluate(next);
  }
  return { x: xBest, f: fBest, history, X, y };
}
