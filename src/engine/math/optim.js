/**
 * Optimisation: one-dimensional search, local derivative-free and quasi-Newton methods, non-linear least squares,
 * population-based global methods, multi-objective NSGA-II, linear and mixed-integer linear programming.
 *
 * All methods minimise. Bounds are arrays of [lo, hi] pairs, one per variable (±Infinity allowed for the local
 * methods). Stochastic methods draw exclusively from the seeded generator of rng.js and are therefore
 * reproducible for a given seed.
 */
import { makeRng } from './rng.js';
import { solveDense, inverse } from './linalg.js';
import { brent } from './roots.js';

function checkFn(f, name) {
  if (typeof f !== 'function') throw new RangeError(`${name}: objective must be a function`);
}

function normBounds(bounds, n, name, finite = false) {
  if (bounds === undefined || bounds === null) {
    if (finite) throw new RangeError(`${name}: bounds are required`);
    return null;
  }
  if (!Array.isArray(bounds) || bounds.length !== n) throw new RangeError(`${name}: bounds must have one [lo, hi] pair per variable`);
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const b = bounds[i];
    lo[i] = b && b[0] !== undefined && b[0] !== null ? b[0] : -Infinity;
    hi[i] = b && b[1] !== undefined && b[1] !== null ? b[1] : Infinity;
    if (!(hi[i] >= lo[i])) throw new RangeError(`${name}: bounds[${i}] must satisfy lo <= hi`);
    if (finite && !(Number.isFinite(lo[i]) && Number.isFinite(hi[i]))) throw new RangeError(`${name}: bounds[${i}] must be finite`);
  }
  return { lo, hi };
}

function clampInto(x, b) {
  if (!b) return x;
  for (let i = 0; i < x.length; i++) {
    if (x[i] < b.lo[i]) x[i] = b.lo[i];
    else if (x[i] > b.hi[i]) x[i] = b.hi[i];
  }
  return x;
}

// ---------------------------------------------------------------------------------------------------------------
// One-dimensional
// ---------------------------------------------------------------------------------------------------------------

/**
 * Golden-section search for the minimum of a unimodal function on [a, b].
 * Reference: Kiefer (1953); Press et al., "Numerical Recipes" §10.2. Validity: f unimodal on [a, b]; otherwise a
 * local minimum is returned.
 * @param {(x: number) => number} f objective
 * @param {number} a lower end of the interval
 * @param {number} b upper end of the interval
 * @param {{ tol?: number, maxIter?: number }} [opts] absolute interval tolerance (default 1e-8·(1+|a|+|b|))
 * @returns {{ x: number, f: number, iterations: number }} minimiser, minimum, iterations used
 */
export function goldenSection(f, a, b, opts = {}) {
  checkFn(f, 'goldenSection');
  if (!(Number.isFinite(a) && Number.isFinite(b))) throw new RangeError('goldenSection: a and b must be finite');
  if (a > b) [a, b] = [b, a];
  const tol = opts.tol ?? 1e-8 * (1 + Math.abs(a) + Math.abs(b));
  const maxIter = opts.maxIter ?? 200;
  const invPhi = (Math.sqrt(5) - 1) / 2;
  let c = b - invPhi * (b - a);
  let d = a + invPhi * (b - a);
  let fc = f(c);
  let fd = f(d);
  let it = 0;
  while (b - a > tol && it < maxIter) {
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - invPhi * (b - a);
      fc = f(c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + invPhi * (b - a);
      fd = f(d);
    }
    it++;
  }
  return fc < fd ? { x: c, f: fc, iterations: it } : { x: d, f: fd, iterations: it };
}

/**
 * Find x in [a, b] with f(x) = target (Brent's method on f - target).
 * Reference: Brent (1973), "Algorithms for Minimization without Derivatives", ch. 4. Validity: f continuous and
 * f - target changes sign over [a, b]; otherwise a RangeError is thrown.
 * @param {(x: number) => number} f function
 * @param {number} target required function value
 * @param {number} a lower end of the search interval
 * @param {number} b upper end of the search interval
 * @returns {number} x such that f(x) = target
 */
export function goalSeek(f, target, a, b) {
  checkFn(f, 'goalSeek');
  const h = (x) => f(x) - target;
  const fa = h(a);
  const fb = h(b);
  if (fa === 0) return a;
  if (fb === 0) return b;
  if (!(fa * fb < 0)) throw new RangeError('goalSeek: f(a) - target and f(b) - target must have opposite signs');
  return brent(h, a, b, { tol: 1e-13, maxIter: 300 });
}

// ---------------------------------------------------------------------------------------------------------------
// Nelder–Mead
// ---------------------------------------------------------------------------------------------------------------

/**
 * Nelder–Mead downhill simplex with optional box bounds (trial points are projected onto the box).
 * Coefficients: reflection 1, expansion 2, contraction 0.5, shrink 0.5.
 * Reference: Nelder & Mead (1965), Comput. J. 7, 308–313; Lagarias et al. (1998), SIAM J. Optim. 9, 112–147.
 * Validity: smooth or mildly noisy objectives, n <= ~20.
 * @param {(x: number[]) => number} f objective
 * @param {ArrayLike<number>} x0 starting point
 * @param {{ bounds?: [number, number][], tol?: number, maxIter?: number, step?: number | number[] }} [opts]
 *   tol: relative spread of function values over the simplex combined with a √tol simplex-size test (1e-8);
 *   step: initial simplex edge (scalar or per variable; default 5 % of |x0|, or of the bound range)
 * @returns {{ x: number[], f: number, iterations: number, converged: boolean, evaluations: number }}
 */
export function nelderMead(f, x0, opts = {}) {
  checkFn(f, 'nelderMead');
  const n = x0.length;
  if (!(n >= 1)) throw new RangeError('nelderMead: x0 must have at least one element');
  const b = normBounds(opts.bounds, n, 'nelderMead');
  const tol = opts.tol ?? 1e-8;
  const maxIter = opts.maxIter ?? 2000;
  let evals = 0;
  const F = (x) => {
    evals++;
    const v = f(x);
    return Number.isNaN(v) ? Infinity : v;
  };
  const pts = new Array(n + 1);
  const fv = new Float64Array(n + 1);
  pts[0] = clampInto(Array.from(x0), b);
  for (let i = 0; i < n; i++) {
    const p = pts[0].slice();
    let h;
    if (opts.step !== undefined) h = typeof opts.step === 'number' ? opts.step : opts.step[i];
    else if (b && Number.isFinite(b.lo[i]) && Number.isFinite(b.hi[i])) h = 0.05 * (b.hi[i] - b.lo[i]);
    else h = p[i] !== 0 ? 0.05 * Math.abs(p[i]) : 0.00025 * 100;
    if (h === 0) h = 1e-3;
    if (b && p[i] + h > b.hi[i]) h = -Math.abs(h);
    p[i] += h;
    clampInto(p, b);
    if (p[i] === pts[0][i]) p[i] = pts[0][i] - 0.5 * h; // degenerate against a bound
    clampInto(p, b);
    pts[i + 1] = p;
  }
  for (let i = 0; i <= n; i++) fv[i] = F(pts[i]);
  const order = new Array(n + 1);
  const cen = new Array(n);
  const trial = (coef, worst) => {
    const p = new Array(n);
    for (let j = 0; j < n; j++) p[j] = cen[j] + coef * (worst[j] - cen[j]);
    return clampInto(p, b);
  };
  let it = 0;
  let converged = false;
  for (; it < maxIter; it++) {
    for (let i = 0; i <= n; i++) order[i] = i;
    order.sort((i, j) => fv[i] - fv[j]);
    const ib = order[0];
    const iw = order[n];
    const isw = order[n - 1];
    const fb = fv[ib];
    let size = 0;
    let scale = 0;
    for (let j = 0; j < n; j++) {
      scale = Math.max(scale, Math.abs(pts[ib][j]));
      for (let i = 0; i <= n; i++) size = Math.max(size, Math.abs(pts[i][j] - pts[ib][j]));
    }
    const spread = fv[iw] - fb;
    if ((spread <= tol * (Math.abs(fb) + tol) && size <= Math.sqrt(tol) * (1 + scale)) || size <= 1e-15 * (1 + scale)) {
      converged = true;
      break;
    }
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let i = 0; i <= n; i++) if (i !== iw) s += pts[i][j];
      cen[j] = s / n;
    }
    const worst = pts[iw];
    const xr = trial(-1, worst);
    const fr = F(xr);
    if (fr < fb) {
      const xe = trial(-2, worst);
      const fe = F(xe);
      if (fe < fr) {
        pts[iw] = xe;
        fv[iw] = fe;
      } else {
        pts[iw] = xr;
        fv[iw] = fr;
      }
    } else if (fr < fv[isw]) {
      pts[iw] = xr;
      fv[iw] = fr;
    } else {
      const outside = fr < fv[iw];
      const xc = trial(outside ? -0.5 : 0.5, worst);
      const fc = F(xc);
      if (fc < (outside ? fr : fv[iw])) {
        pts[iw] = xc;
        fv[iw] = fc;
      } else {
        for (let i = 0; i <= n; i++) {
          if (i === ib) continue;
          for (let j = 0; j < n; j++) pts[i][j] = pts[ib][j] + 0.5 * (pts[i][j] - pts[ib][j]);
          fv[i] = F(pts[i]);
        }
      }
    }
  }
  let ib = 0;
  for (let i = 1; i <= n; i++) if (fv[i] < fv[ib]) ib = i;
  return { x: pts[ib].slice(), f: fv[ib], iterations: it, converged, evaluations: evals };
}

// ---------------------------------------------------------------------------------------------------------------
// Levenberg–Marquardt
// ---------------------------------------------------------------------------------------------------------------

/**
 * Levenberg–Marquardt non-linear least squares, minimising cost = ½ Σ r_i(x)², with Marquardt diagonal scaling
 * and projection onto optional box bounds.
 * The parameter covariance is (JᵀJ)⁻¹·σ² with σ² = Σ r²/(m - n) evaluated at the solution (zero when m <= n).
 * Reference: Marquardt (1963), J. SIAM 11, 431–441; Moré (1978), "The Levenberg–Marquardt algorithm:
 * implementation and theory". Validity: m residuals >= 1, smooth residuals; forward-difference Jacobian when
 * `jac` is omitted.
 * @param {(x: number[]) => ArrayLike<number>} residuals residual vector function (length m)
 * @param {ArrayLike<number>} x0 starting parameters (length n)
 * @param {{ bounds?: [number, number][], tol?: number, maxIter?: number, jac?: (x: number[]) => number[][] }} [opts]
 *   tol: relative step / relative cost-reduction / gradient tolerance (1e-10); jac returns the m×n Jacobian
 * @returns {{ x: number[], cost: number, iterations: number, converged: boolean, covariance: number[][],
 *   residuals: number[], sigma2: number, se: number[] }} solution, ½Σr², iterations, convergence flag, covariance
 *   matrix (n×n), residual vector at the solution, residual variance, parameter standard errors
 */
export function levenbergMarquardt(residuals, x0, opts = {}) {
  checkFn(residuals, 'levenbergMarquardt');
  const n = x0.length;
  if (!(n >= 1)) throw new RangeError('levenbergMarquardt: x0 must have at least one element');
  const b = normBounds(opts.bounds, n, 'levenbergMarquardt');
  const tol = opts.tol ?? 1e-10;
  const maxIter = opts.maxIter ?? 200;
  let x = clampInto(Array.from(x0), b);
  let r = Array.from(residuals(x));
  const m = r.length;
  if (!(m >= 1)) throw new RangeError('levenbergMarquardt: residuals must return a non-empty array');
  const half = (v) => {
    let s = 0;
    for (let i = 0; i < v.length; i++) s += v[i] * v[i];
    return Number.isNaN(s) ? Infinity : 0.5 * s;
  };
  let cost = half(r);
  if (!Number.isFinite(cost)) throw new RangeError('levenbergMarquardt: residuals are not finite at x0');

  const jacobian = (xx, rr) => {
    if (opts.jac) return opts.jac(xx);
    const J = new Array(m);
    for (let i = 0; i < m; i++) J[i] = new Float64Array(n);
    const xp = xx.slice();
    for (let j = 0; j < n; j++) {
      let h = 1e-7 * Math.max(Math.abs(xx[j]), 1e-3);
      if (b && xx[j] + h > b.hi[j]) h = -h;
      xp[j] = xx[j] + h;
      const rp = residuals(xp);
      xp[j] = xx[j];
      for (let i = 0; i < m; i++) J[i][j] = (rp[i] - rr[i]) / h;
    }
    return J;
  };
  const normal = (J, rr) => {
    const A = new Array(n);
    const g = new Array(n).fill(0);
    for (let a = 0; a < n; a++) A[a] = new Array(n).fill(0);
    for (let i = 0; i < m; i++) {
      const Ji = J[i];
      const ri = rr[i];
      for (let a = 0; a < n; a++) {
        const ja = Ji[a];
        if (ja === 0) continue;
        g[a] += ja * ri;
        const Aa = A[a];
        for (let c = 0; c <= a; c++) Aa[c] += ja * Ji[c];
      }
    }
    for (let a = 0; a < n; a++) for (let c = a + 1; c < n; c++) A[a][c] = A[c][a];
    return { A, g };
  };

  let lambda = 1e-3;
  let it = 0;
  let converged = false;
  let J = jacobian(x, r);
  while (it < maxIter) {
    it++;
    const { A, g } = normal(J, r);
    // projected gradient test
    let gmax = 0;
    for (let j = 0; j < n; j++) {
      const atLo = b && x[j] <= b.lo[j] && g[j] > 0;
      const atHi = b && x[j] >= b.hi[j] && g[j] < 0;
      if (!atLo && !atHi) gmax = Math.max(gmax, Math.abs(g[j]));
    }
    if (gmax <= 1e-14 * (1 + cost) || cost <= 1e-300) {
      converged = true;
      break;
    }
    let accepted = false;
    let done = false;
    for (let inner = 0; inner < 40; inner++) {
      const M = new Array(n);
      for (let a = 0; a < n; a++) {
        M[a] = A[a].slice();
        M[a][a] += lambda * (A[a][a] > 0 ? A[a][a] : 1) + 1e-300;
      }
      let delta;
      try {
        delta = solveDense(
          M,
          g.map((v) => -v),
        );
      } catch (e) {
        lambda *= 10;
        continue;
      }
      const xn = new Array(n);
      let bad = false;
      for (let j = 0; j < n; j++) {
        xn[j] = x[j] + delta[j];
        if (!Number.isFinite(xn[j])) bad = true;
      }
      if (bad) {
        lambda *= 10;
        continue;
      }
      clampInto(xn, b);
      const rn = Array.from(residuals(xn));
      const cn = half(rn);
      let step = 0;
      let xs = 0;
      for (let j = 0; j < n; j++) {
        step = Math.max(step, Math.abs(xn[j] - x[j]));
        xs = Math.max(xs, Math.abs(x[j]));
      }
      if (cn < cost) {
        const red = cost - cn;
        const smallStep = step <= tol * (xs + tol);
        const smallRed = red <= tol * cost;
        x = xn;
        r = rn;
        cost = cn;
        lambda = Math.max(lambda / 10, 1e-15);
        accepted = true;
        if (smallStep || smallRed) done = true;
        break;
      }
      if (step <= tol * (xs + tol) * 1e-2 || lambda > 1e15) {
        done = true; // no further progress possible
        break;
      }
      lambda *= 10;
    }
    if (accepted) J = jacobian(x, r);
    if (done || !accepted) {
      converged = true;
      break;
    }
  }
  const { A } = normal(J, r);
  const sigma2 = m > n ? (2 * cost) / (m - n) : 0;
  let covariance;
  try {
    const inv = inverse(A);
    covariance = inv.map((row) => Array.from(row, (v) => v * sigma2));
  } catch (e) {
    covariance = Array.from({ length: n }, () => new Array(n).fill(Infinity));
  }
  const se = covariance.map((row, i) => Math.sqrt(Math.max(row[i], 0)));
  return { x, cost, iterations: it, converged, covariance, residuals: r, sigma2, se };
}

// ---------------------------------------------------------------------------------------------------------------
// BFGS
// ---------------------------------------------------------------------------------------------------------------

/**
 * Quasi-Newton BFGS minimisation with Armijo backtracking; box bounds are handled by a projected active-set
 * scheme (variables at a bound whose gradient points outward are frozen for the step, trial points are
 * projected onto the box).
 * Reference: Nocedal & Wright, "Numerical Optimization", 2nd ed., §6.1 (BFGS) and §16.7 / Bertsekas (1982)
 * (gradient projection). Validity: smooth objectives; central-difference gradient when `grad` is omitted.
 * @param {(x: number[]) => number} f objective
 * @param {ArrayLike<number>} x0 starting point
 * @param {{ grad?: (x: number[]) => ArrayLike<number>, bounds?: [number, number][], tol?: number, maxIter?: number }} [opts]
 *   tol: infinity norm of the projected gradient (default 1e-8)
 * @returns {{ x: number[], f: number, iterations: number, converged: boolean, gradNorm: number, evaluations: number }}
 */
export function bfgs(f, x0, opts = {}) {
  checkFn(f, 'bfgs');
  const n = x0.length;
  if (!(n >= 1)) throw new RangeError('bfgs: x0 must have at least one element');
  const b = normBounds(opts.bounds, n, 'bfgs');
  const tol = opts.tol ?? 1e-8;
  const maxIter = opts.maxIter ?? 500;
  let evals = 0;
  const F = (x) => {
    evals++;
    const v = f(x);
    return Number.isNaN(v) ? Infinity : v;
  };
  const numGrad = (x) => {
    const g = new Array(n);
    const xp = x.slice();
    for (let j = 0; j < n; j++) {
      const h = 6e-6 * Math.max(Math.abs(x[j]), 1e-2);
      let hp = h;
      let hm = h;
      if (b) {
        hp = Math.min(h, b.hi[j] - x[j]);
        hm = Math.min(h, x[j] - b.lo[j]);
        if (hp + hm === 0) {
          g[j] = 0;
          continue;
        }
      }
      xp[j] = x[j] + hp;
      const fp = hp > 0 ? F(xp) : NaN;
      xp[j] = x[j] - hm;
      const fm = hm > 0 ? F(xp) : NaN;
      xp[j] = x[j];
      if (hp > 0 && hm > 0) g[j] = (fp - fm) / (hp + hm);
      else if (hp > 0) g[j] = (fp - F(x)) / hp;
      else g[j] = (F(x) - fm) / hm;
    }
    return g;
  };
  const G = opts.grad ? (x) => Array.from(opts.grad(x)) : numGrad;
  let x = clampInto(Array.from(x0), b);
  let fx = F(x);
  if (!Number.isFinite(fx)) throw new RangeError('bfgs: objective is not finite at x0');
  let g = G(x);
  const H = new Array(n);
  const resetH = () => {
    for (let i = 0; i < n; i++) {
      H[i] = new Float64Array(n);
      H[i][i] = 1;
    }
  };
  resetH();
  const active = new Uint8Array(n);
  const d = new Array(n);
  let it = 0;
  let converged = false;
  let pg = Infinity;
  let fresh = true;
  for (; it < maxIter; it++) {
    pg = 0;
    for (let j = 0; j < n; j++) {
      active[j] = b && ((x[j] <= b.lo[j] && g[j] > 0) || (x[j] >= b.hi[j] && g[j] < 0)) ? 1 : 0;
      if (!active[j]) pg = Math.max(pg, Math.abs(g[j]));
    }
    if (pg <= tol) {
      converged = true;
      break;
    }
    let slope = 0;
    for (let i = 0; i < n; i++) {
      if (active[i]) {
        d[i] = 0;
        continue;
      }
      let s = 0;
      for (let j = 0; j < n; j++) if (!active[j]) s -= H[i][j] * g[j];
      d[i] = s;
      slope += s * g[i];
    }
    if (!(slope < 0)) {
      resetH();
      fresh = true;
      slope = 0;
      for (let i = 0; i < n; i++) {
        d[i] = active[i] ? 0 : -g[i];
        slope += d[i] * g[i];
      }
    }
    // initial step: unit for a built-up Hessian, scaled for steepest descent
    let t = 1;
    if (fresh) {
      let dn = 0;
      let xn = 0;
      for (let i = 0; i < n; i++) {
        dn = Math.max(dn, Math.abs(d[i]));
        xn = Math.max(xn, Math.abs(x[i]));
      }
      t = Math.min(1, Math.max(xn, 1) / dn);
    }
    let xn;
    let fn = Infinity;
    let ok = false;
    for (let ls = 0; ls < 60; ls++) {
      xn = new Array(n);
      for (let i = 0; i < n; i++) xn[i] = x[i] + t * d[i];
      clampInto(xn, b);
      fn = F(xn);
      let dec = 0;
      for (let i = 0; i < n; i++) dec += g[i] * (xn[i] - x[i]);
      if (fn <= fx + 1e-4 * dec && dec < 0) {
        ok = true;
        break;
      }
      t *= ls < 3 ? 0.5 : 0.25;
    }
    if (!ok) {
      if (!fresh) {
        resetH();
        fresh = true;
        continue; // retry along steepest descent
      }
      converged = pg <= Math.max(1e3 * tol, 1e-5) * (1 + Math.abs(fx));
      break;
    }
    const gn = G(xn);
    const s = new Array(n);
    const y = new Array(n);
    let sy = 0;
    let ss = 0;
    let yy = 0;
    for (let i = 0; i < n; i++) {
      s[i] = xn[i] - x[i];
      y[i] = gn[i] - g[i];
      sy += s[i] * y[i];
      ss += s[i] * s[i];
      yy += y[i] * y[i];
    }
    if (sy > 1e-10 * Math.sqrt(ss * yy) && sy > 0) {
      if (fresh) {
        const sc = sy / yy; // initial Hessian scaling (Nocedal & Wright eq. 6.20)
        for (let i = 0; i < n; i++) H[i][i] = sc;
        fresh = false;
      }
      const Hy = new Array(n);
      let yHy = 0;
      for (let i = 0; i < n; i++) {
        let a = 0;
        for (let j = 0; j < n; j++) a += H[i][j] * y[j];
        Hy[i] = a;
        yHy += a * y[i];
      }
      const rho = 1 / sy;
      const c = (1 + yHy * rho) * rho;
      for (let i = 0; i < n; i++) {
        const Hi = H[i];
        for (let j = 0; j < n; j++) Hi[j] += c * s[i] * s[j] - rho * (Hy[i] * s[j] + s[i] * Hy[j]);
      }
    }
    const df = fx - fn;
    x = xn;
    fx = fn;
    g = gn;
    if (df <= 1e-16 * (1 + Math.abs(fx)) && Math.sqrt(ss) <= 1e-14 * (1 + Math.sqrt(x.reduce((a, v) => a + v * v, 0)))) {
      converged = true;
      break;
    }
  }
  return { x, f: fx, iterations: it, converged, gradNorm: pg, evaluations: evals };
}

// ---------------------------------------------------------------------------------------------------------------
// Population-based global methods
// ---------------------------------------------------------------------------------------------------------------

function randomPoint(b, n, rng) {
  const x = new Array(n);
  for (let j = 0; j < n; j++) x[j] = b.lo[j] + (b.hi[j] - b.lo[j]) * rng.next();
  return x;
}

function evalSafe(f, x) {
  const v = f(x);
  return Number.isNaN(v) ? Infinity : v;
}

/**
 * Differential evolution, DE/rand/1/bin with reflection at the bounds.
 * Reference: Storn & Price (1997), J. Global Optim. 11, 341–359. Validity: finite bounds; continuous variables.
 * @param {(x: number[]) => number} f objective
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ pop?: number, F?: number, CR?: number, maxGen?: number, seed?: number, tol?: number }} [opts]
 *   pop population size (15·n, at least 5); F differential weight (0.7); CR crossover rate (0.9); maxGen (200);
 *   seed (1); tol: stop when max(f) - min(f) over the population <= tol·(1 + |min f|) (1e-12)
 * @returns {{ x: number[], f: number, generations: number, history: number[], evaluations: number }}
 *   history[k] = best objective after generation k
 */
export function differentialEvolution(f, bounds, opts = {}) {
  checkFn(f, 'differentialEvolution');
  const n = bounds.length;
  const b = normBounds(bounds, n, 'differentialEvolution', true);
  const np = Math.max(5, Math.round(opts.pop ?? 15 * n));
  const Fw = opts.F ?? 0.7;
  const CR = opts.CR ?? 0.9;
  const maxGen = opts.maxGen ?? 200;
  const tol = opts.tol ?? 1e-12;
  const rng = makeRng(opts.seed ?? 1);
  const P = new Array(np);
  const fv = new Float64Array(np);
  let evals = 0;
  for (let i = 0; i < np; i++) {
    // stratified initial population
    P[i] = new Array(n);
  }
  for (let j = 0; j < n; j++) {
    const perm = rng.shuffle(Array.from({ length: np }, (_, i) => i));
    for (let i = 0; i < np; i++) P[i][j] = b.lo[j] + ((perm[i] + rng.next()) / np) * (b.hi[j] - b.lo[j]);
  }
  for (let i = 0; i < np; i++) {
    fv[i] = evalSafe(f, P[i]);
    evals++;
  }
  const history = [];
  const trial = new Array(n);
  let gen = 0;
  for (; gen < maxGen; gen++) {
    for (let i = 0; i < np; i++) {
      let r1;
      let r2;
      let r3;
      do r1 = rng.int(np);
      while (r1 === i);
      do r2 = rng.int(np);
      while (r2 === i || r2 === r1);
      do r3 = rng.int(np);
      while (r3 === i || r3 === r1 || r3 === r2);
      const jr = rng.int(n);
      for (let j = 0; j < n; j++) {
        if (j === jr || rng.next() < CR) {
          let v = P[r1][j] + Fw * (P[r2][j] - P[r3][j]);
          if (v < b.lo[j]) v = Math.min(b.hi[j], 2 * b.lo[j] - v);
          else if (v > b.hi[j]) v = Math.max(b.lo[j], 2 * b.hi[j] - v);
          trial[j] = v;
        } else trial[j] = P[i][j];
      }
      const ft = evalSafe(f, trial);
      evals++;
      if (ft <= fv[i]) {
        fv[i] = ft;
        for (let j = 0; j < n; j++) P[i][j] = trial[j];
      }
    }
    let fmin = Infinity;
    let fmax = -Infinity;
    for (let i = 0; i < np; i++) {
      if (fv[i] < fmin) fmin = fv[i];
      if (fv[i] > fmax) fmax = fv[i];
    }
    history.push(fmin);
    if (fmax - fmin <= tol * (1 + Math.abs(fmin))) {
      gen++;
      break;
    }
  }
  let ib = 0;
  for (let i = 1; i < np; i++) if (fv[i] < fv[ib]) ib = i;
  return { x: P[ib].slice(), f: fv[ib], generations: gen, history, evaluations: evals };
}

/**
 * Particle swarm optimisation (global-best topology, constriction-type coefficients, velocity clamping,
 * absorbing walls).
 * Reference: Kennedy & Eberhart (1995), Proc. IEEE ICNN; Clerc & Kennedy (2002), IEEE Trans. Evol. Comput. 6,
 * 58–73 (w = 0.7298, c1 = c2 = 1.49618). Validity: finite bounds.
 * @param {(x: number[]) => number} f objective
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ particles?: number, maxIter?: number, seed?: number, w?: number, c1?: number, c2?: number }} [opts]
 *   particles (max(30, 10·n)); maxIter (200); seed (1); inertia w; cognitive c1; social c2
 * @returns {{ x: number[], f: number, history: number[], evaluations: number }} history[k] = best after iteration k
 */
export function particleSwarm(f, bounds, opts = {}) {
  checkFn(f, 'particleSwarm');
  const n = bounds.length;
  const b = normBounds(bounds, n, 'particleSwarm', true);
  const np = Math.max(2, Math.round(opts.particles ?? Math.max(30, 10 * n)));
  const maxIter = opts.maxIter ?? 200;
  const w = opts.w ?? 0.7298;
  const c1 = opts.c1 ?? 1.49618;
  const c2 = opts.c2 ?? 1.49618;
  const rng = makeRng(opts.seed ?? 1);
  const X = new Array(np);
  const V = new Array(np);
  const PB = new Array(np);
  const fpb = new Float64Array(np);
  let gb = null;
  let fgb = Infinity;
  let evals = 0;
  for (let i = 0; i < np; i++) {
    X[i] = randomPoint(b, n, rng);
    V[i] = new Array(n);
    for (let j = 0; j < n; j++) V[i][j] = (rng.next() - 0.5) * 0.2 * (b.hi[j] - b.lo[j]);
    PB[i] = X[i].slice();
    fpb[i] = evalSafe(f, X[i]);
    evals++;
    if (fpb[i] < fgb) {
      fgb = fpb[i];
      gb = X[i].slice();
    }
  }
  const history = [];
  for (let it = 0; it < maxIter; it++) {
    for (let i = 0; i < np; i++) {
      const x = X[i];
      const v = V[i];
      for (let j = 0; j < n; j++) {
        const range = b.hi[j] - b.lo[j];
        let vj = w * v[j] + c1 * rng.next() * (PB[i][j] - x[j]) + c2 * rng.next() * (gb[j] - x[j]);
        const vmax = 0.5 * range;
        if (vj > vmax) vj = vmax;
        else if (vj < -vmax) vj = -vmax;
        let xj = x[j] + vj;
        if (xj < b.lo[j]) {
          xj = b.lo[j];
          vj = 0;
        } else if (xj > b.hi[j]) {
          xj = b.hi[j];
          vj = 0;
        }
        x[j] = xj;
        v[j] = vj;
      }
      const fx = evalSafe(f, x);
      evals++;
      if (fx < fpb[i]) {
        fpb[i] = fx;
        for (let j = 0; j < n; j++) PB[i][j] = x[j];
      }
    }
    for (let i = 0; i < np; i++) {
      if (fpb[i] < fgb) {
        fgb = fpb[i];
        gb = PB[i].slice();
      }
    }
    history.push(fgb);
  }
  return { x: gb, f: fgb, history, evaluations: evals };
}

/** Simulated binary crossover (Deb & Agrawal 1995) of one gene pair within [lo, hi]. */
function sbxGene(p1, p2, lo, hi, eta, rng) {
  if (Math.abs(p1 - p2) < 1e-14) return [p1, p2];
  const y1 = Math.min(p1, p2);
  const y2 = Math.max(p1, p2);
  const u = rng.next();
  const child = (beta) => {
    const alpha = 2 - Math.pow(beta, -(eta + 1));
    return u <= 1 / alpha ? Math.pow(u * alpha, 1 / (eta + 1)) : Math.pow(1 / (2 - u * alpha), 1 / (eta + 1));
  };
  const range = y2 - y1;
  const bq1 = child(1 + (2 * (y1 - lo)) / range);
  const bq2 = child(1 + (2 * (hi - y2)) / range);
  let c1 = 0.5 * (y1 + y2 - bq1 * range);
  let c2 = 0.5 * (y1 + y2 + bq2 * range);
  c1 = Math.min(Math.max(c1, lo), hi);
  c2 = Math.min(Math.max(c2, lo), hi);
  return rng.next() < 0.5 ? [c2, c1] : [c1, c2];
}

/** Polynomial mutation (Deb & Goyal 1996) of one gene within [lo, hi]. */
function polyMutate(x, lo, hi, eta, rng) {
  const range = hi - lo;
  if (!(range > 0)) return x;
  const d1 = (x - lo) / range;
  const d2 = (hi - x) / range;
  const u = rng.next();
  const mp = 1 / (eta + 1);
  let dq;
  if (u < 0.5) dq = Math.pow(2 * u + (1 - 2 * u) * Math.pow(1 - d1, eta + 1), mp) - 1;
  else dq = 1 - Math.pow(2 * (1 - u) + 2 * (u - 0.5) * Math.pow(1 - d2, eta + 1), mp);
  return Math.min(Math.max(x + dq * range, lo), hi);
}

/**
 * Real-coded genetic algorithm: binary tournament selection, simulated binary crossover (SBX, η = 10),
 * polynomial mutation (η = 20, index annealed to 60 for late-stage refinement) and elitism of the two best
 * individuals. Variables flagged in `integer` are kept integral.
 * Reference: Deb, "Multi-Objective Optimization using Evolutionary Algorithms" (2001), §4.2; Goldberg (1989).
 * Validity: finite bounds.
 * @param {(x: number[]) => number} f objective
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ pop?: number, generations?: number, seed?: number, mutation?: number, crossover?: number,
 *   integer?: boolean[] }} [opts] pop (max(40, 10·n), made even); generations (200); seed (1); mutation
 *   probability per gene (max(1/n, 0.1)); crossover probability per pair (0.9)
 * @returns {{ x: number[], f: number, history: number[], evaluations: number }} history[k] = best after generation k
 */
export function geneticAlgorithm(f, bounds, opts = {}) {
  checkFn(f, 'geneticAlgorithm');
  const n = bounds.length;
  const b = normBounds(bounds, n, 'geneticAlgorithm', true);
  let np = Math.max(4, Math.round(opts.pop ?? Math.max(40, 10 * n)));
  if (np % 2) np++;
  const gens = opts.generations ?? 200;
  const pm = opts.mutation ?? Math.max(1 / n, 0.1);
  const pc = opts.crossover ?? 0.9;
  const integer = opts.integer || [];
  const rng = makeRng(opts.seed ?? 1);
  const fix = (x) => {
    for (let j = 0; j < n; j++) {
      if (integer[j]) x[j] = Math.min(Math.max(Math.round(x[j]), Math.ceil(b.lo[j])), Math.floor(b.hi[j]));
    }
    return x;
  };
  let P = new Array(np);
  let fv = new Float64Array(np);
  let evals = 0;
  for (let i = 0; i < np; i++) {
    P[i] = fix(randomPoint(b, n, rng));
    fv[i] = evalSafe(f, P[i]);
    evals++;
  }
  const history = [];
  const pick = () => {
    const i = rng.int(np);
    const j = rng.int(np);
    return fv[i] <= fv[j] ? i : j;
  };
  for (let g = 0; g < gens; g++) {
    const order = Array.from({ length: np }, (_, i) => i).sort((i, j) => fv[i] - fv[j]);
    const Q = [P[order[0]].slice(), P[order[1]].slice()];
    const fq = [fv[order[0]], fv[order[1]]];
    const etaM = 20 + (40 * g) / Math.max(gens - 1, 1);
    while (Q.length < np) {
      const a = P[pick()];
      const c = P[pick()];
      const c1 = a.slice();
      const c2 = c.slice();
      if (rng.next() < pc) {
        for (let j = 0; j < n; j++) {
          if (rng.next() < 0.5 || n === 1) {
            const pair = sbxGene(a[j], c[j], b.lo[j], b.hi[j], 10, rng);
            c1[j] = pair[0];
            c2[j] = pair[1];
          }
        }
      }
      for (let j = 0; j < n; j++) {
        if (rng.next() < pm) c1[j] = integer[j] ? c1[j] + (rng.next() < 0.5 ? -1 : 1) * (1 + rng.int(Math.max(1, Math.round(0.1 * (b.hi[j] - b.lo[j]))))) : polyMutate(c1[j], b.lo[j], b.hi[j], etaM, rng);
        if (rng.next() < pm) c2[j] = integer[j] ? c2[j] + (rng.next() < 0.5 ? -1 : 1) * (1 + rng.int(Math.max(1, Math.round(0.1 * (b.hi[j] - b.lo[j]))))) : polyMutate(c2[j], b.lo[j], b.hi[j], etaM, rng);
      }
      clampInto(c1, b);
      clampInto(c2, b);
      fix(c1);
      fix(c2);
      Q.push(c1);
      fq.push(evalSafe(f, c1));
      evals++;
      if (Q.length < np) {
        Q.push(c2);
        fq.push(evalSafe(f, c2));
        evals++;
      }
    }
    P = Q;
    fv = Float64Array.from(fq);
    let fb = Infinity;
    for (let i = 0; i < np; i++) if (fv[i] < fb) fb = fv[i];
    history.push(fb);
  }
  let ib = 0;
  for (let i = 1; i < np; i++) if (fv[i] < fv[ib]) ib = i;
  return { x: P[ib].slice(), f: fv[ib], history, evaluations: evals };
}

/**
 * Simulated annealing with Metropolis acceptance, geometric cooling and adaptive Gaussian step lengths
 * (step scale tuned every 50 moves towards ~40 % acceptance, as in Corana et al.).
 * Reference: Kirkpatrick, Gelatt & Vecchi (1983), Science 220, 671–680; Corana et al. (1987), ACM TOMS 13,
 * 262–280. Validity: finite bounds.
 * @param {(x: number[]) => number} f objective
 * @param {ArrayLike<number>} x0 starting point (projected onto the bounds)
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ T0?: number, cooling?: number, maxIter?: number, seed?: number }} [opts]
 *   T0 initial temperature in objective units (default: mean |Δf| of 20 random trial moves); cooling: temperature
 *   factor applied every iteration (default chosen so that T ends at 1e-10·T0); maxIter (20000); seed (1)
 * @returns {{ x: number[], f: number, history: number[], evaluations: number }}
 *   history: best objective sampled at up to 200 equally spaced iterations
 */
export function simulatedAnnealing(f, x0, bounds, opts = {}) {
  checkFn(f, 'simulatedAnnealing');
  const n = x0.length;
  const b = normBounds(bounds, n, 'simulatedAnnealing', true);
  const maxIter = opts.maxIter ?? 20000;
  const rng = makeRng(opts.seed ?? 1);
  let x = clampInto(Array.from(x0), b);
  let fx = evalSafe(f, x);
  let evals = 1;
  let best = x.slice();
  let fbest = fx;
  let scale = 0.1;
  const cand = new Array(n);
  const propose = () => {
    for (let j = 0; j < n; j++) {
      const range = b.hi[j] - b.lo[j];
      let v = x[j] + scale * range * rng.normal();
      // reflect into the box
      for (let k = 0; k < 8 && (v < b.lo[j] || v > b.hi[j]); k++) v = v < b.lo[j] ? 2 * b.lo[j] - v : 2 * b.hi[j] - v;
      cand[j] = Math.min(Math.max(v, b.lo[j]), b.hi[j]);
    }
  };
  let T0 = opts.T0;
  if (T0 === undefined) {
    let s = 0;
    let c = 0;
    for (let k = 0; k < 20; k++) {
      propose();
      const fc = evalSafe(f, cand);
      evals++;
      if (Number.isFinite(fc) && Number.isFinite(fx)) {
        s += Math.abs(fc - fx);
        c++;
      }
      if (fc < fbest) {
        fbest = fc;
        best = cand.slice();
      }
    }
    T0 = c > 0 && s > 0 ? s / c : 1;
  }
  if (!(T0 > 0)) throw new RangeError('simulatedAnnealing: T0 must be > 0');
  const cooling = opts.cooling ?? Math.pow(1e-10, 1 / Math.max(maxIter, 1));
  if (!(cooling > 0 && cooling < 1)) throw new RangeError('simulatedAnnealing: cooling must be within (0, 1)');
  let T = T0;
  const history = [];
  const every = Math.max(1, Math.floor(maxIter / 200));
  let acc = 0;
  for (let it = 0; it < maxIter; it++) {
    propose();
    const fc = evalSafe(f, cand);
    evals++;
    const df = fc - fx;
    if (df <= 0 || rng.next() < Math.exp(-df / T)) {
      for (let j = 0; j < n; j++) x[j] = cand[j];
      fx = fc;
      acc++;
      if (fx < fbest) {
        fbest = fx;
        best = x.slice();
      }
    }
    if ((it + 1) % 50 === 0) {
      const rate = acc / 50;
      if (rate > 0.5) scale = Math.min(scale * 1.5, 0.5);
      else if (rate < 0.3) scale = Math.max(scale / 1.5, 1e-12);
      acc = 0;
      // restart from the best point when the walk has drifted far above it
      if (fx > fbest + 10 * T) {
        x = best.slice();
        fx = fbest;
      }
    }
    T *= cooling;
    if ((it + 1) % every === 0) history.push(fbest);
  }
  return { x: best, f: fbest, history, evaluations: evals };
}

// ---------------------------------------------------------------------------------------------------------------
// NSGA-II
// ---------------------------------------------------------------------------------------------------------------

/**
 * Fast non-dominated sort (Deb et al. 2002, §III-A).
 * @param {ArrayLike<number>[]} F objective vectors (all minimised)
 * @returns {{ rank: Int32Array, fronts: number[][] }} rank per individual (0 = non-dominated) and index lists
 */
export function nonDominatedSort(F) {
  const N = F.length;
  const m = N ? F[0].length : 0;
  const dominatedBy = new Int32Array(N);
  const dominates = new Array(N);
  const rank = new Int32Array(N);
  const fronts = [[]];
  for (let p = 0; p < N; p++) dominates[p] = [];
  for (let p = 0; p < N; p++) {
    const fp = F[p];
    for (let q = p + 1; q < N; q++) {
      const fq = F[q];
      let pb = false;
      let qb = false;
      for (let k = 0; k < m; k++) {
        if (fp[k] < fq[k]) pb = true;
        else if (fp[k] > fq[k]) qb = true;
      }
      if (pb && !qb) {
        dominates[p].push(q);
        dominatedBy[q]++;
      } else if (qb && !pb) {
        dominates[q].push(p);
        dominatedBy[p]++;
      }
    }
  }
  for (let p = 0; p < N; p++) if (dominatedBy[p] === 0) fronts[0].push(p);
  let i = 0;
  while (fronts[i].length) {
    const nextFront = [];
    for (const p of fronts[i]) {
      for (const q of dominates[p]) {
        if (--dominatedBy[q] === 0) {
          rank[q] = i + 1;
          nextFront.push(q);
        }
      }
    }
    fronts.push(nextFront);
    i++;
  }
  fronts.pop();
  return { rank, fronts };
}

/**
 * Crowding distance of the members of one front (Deb et al. 2002, §III-B).
 * @param {ArrayLike<number>[]} F objective vectors of the whole population
 * @param {number[]} front indices of the front members
 * @param {Float64Array} out crowding distance written at the member indices
 */
export function crowdingDistance(F, front, out) {
  const l = front.length;
  for (const i of front) out[i] = 0;
  if (l === 0) return out;
  const m = F[front[0]].length;
  const idx = front.slice();
  for (let k = 0; k < m; k++) {
    idx.sort((a, c) => F[a][k] - F[c][k]);
    const fmin = F[idx[0]][k];
    const fmax = F[idx[l - 1]][k];
    out[idx[0]] = Infinity;
    out[idx[l - 1]] = Infinity;
    if (!(fmax > fmin)) continue;
    for (let i = 1; i < l - 1; i++) out[idx[i]] += (F[idx[i + 1]][k] - F[idx[i - 1]][k]) / (fmax - fmin);
  }
  return out;
}

/**
 * NSGA-II multi-objective optimisation: fast non-dominated sorting, crowding-distance diversity preservation,
 * crowded binary tournament, SBX crossover (pc = 0.9, η = 15) and polynomial mutation (pm = 1/n, η = 20).
 * Reference: Deb, Pratap, Agarwal & Meyarivan (2002), IEEE Trans. Evol. Comput. 6, 182–197. Validity: finite
 * bounds, 2–4 objectives (crowding distance loses resolving power beyond that).
 * @param {(x: number[]) => ArrayLike<number>} fs objective vector function (all objectives minimised)
 * @param {[number, number][]} bounds finite [lo, hi] per variable
 * @param {{ pop?: number, generations?: number, seed?: number }} [opts] pop (100, made a multiple of 2);
 *   generations (200); seed (1)
 * @returns {{ front: { x: number[], f: number[] }[], evaluations: number }} non-dominated set of the final
 *   population, sorted by the first objective
 */
export function nsga2(fs, bounds, opts = {}) {
  checkFn(fs, 'nsga2');
  const n = bounds.length;
  const b = normBounds(bounds, n, 'nsga2', true);
  let np = Math.max(4, Math.round(opts.pop ?? 100));
  if (np % 2) np++;
  const gens = opts.generations ?? 200;
  const rng = makeRng(opts.seed ?? 1);
  const pm = 1 / n;
  let evals = 0;
  const evalF = (x) => {
    evals++;
    return Array.from(fs(x), (v) => (Number.isNaN(v) ? Infinity : v));
  };
  let X = new Array(np);
  let F = new Array(np);
  for (let i = 0; i < np; i++) {
    X[i] = randomPoint(b, n, rng);
    F[i] = evalF(X[i]);
  }
  let { rank, fronts } = nonDominatedSort(F);
  let crowd = new Float64Array(np);
  for (const fr of fronts) crowdingDistance(F, fr, crowd);
  const tournament = () => {
    const i = rng.int(np);
    const j = rng.int(np);
    if (rank[i] !== rank[j]) return rank[i] < rank[j] ? i : j;
    if (crowd[i] !== crowd[j]) return crowd[i] > crowd[j] ? i : j;
    return rng.next() < 0.5 ? i : j;
  };
  for (let g = 0; g < gens; g++) {
    const XA = X.slice();
    const FA = F.slice();
    while (XA.length < 2 * np) {
      const a = X[tournament()];
      const c = X[tournament()];
      const c1 = a.slice();
      const c2 = c.slice();
      if (rng.next() < 0.9) {
        for (let j = 0; j < n; j++) {
          if (rng.next() < 0.5) {
            const pair = sbxGene(a[j], c[j], b.lo[j], b.hi[j], 15, rng);
            c1[j] = pair[0];
            c2[j] = pair[1];
          }
        }
      }
      for (let j = 0; j < n; j++) {
        if (rng.next() < pm) c1[j] = polyMutate(c1[j], b.lo[j], b.hi[j], 20, rng);
        if (rng.next() < pm) c2[j] = polyMutate(c2[j], b.lo[j], b.hi[j], 20, rng);
      }
      XA.push(c1, c2);
      FA.push(evalF(c1), evalF(c2));
    }
    const sorted = nonDominatedSort(FA);
    const cd = new Float64Array(2 * np);
    const keep = [];
    for (const fr of sorted.fronts) {
      crowdingDistance(FA, fr, cd);
      if (keep.length + fr.length <= np) {
        for (const i of fr) keep.push(i);
      } else {
        const rest = fr.slice().sort((i, j) => cd[j] - cd[i]);
        for (let k = 0; keep.length < np; k++) keep.push(rest[k]);
      }
      if (keep.length >= np) break;
    }
    X = keep.map((i) => XA[i]);
    F = keep.map((i) => FA[i]);
    const r2 = new Int32Array(np);
    const c2 = new Float64Array(np);
    for (let k = 0; k < np; k++) {
      r2[k] = sorted.rank[keep[k]];
      c2[k] = cd[keep[k]];
    }
    rank = r2;
    crowd = c2;
  }
  ({ rank, fronts } = nonDominatedSort(F));
  const front = fronts[0].map((i) => ({ x: X[i].slice(), f: F[i].slice() }));
  front.sort((p, q) => p.f[0] - q.f[0]);
  return { front, evaluations: evals };
}

// ---------------------------------------------------------------------------------------------------------------
// Linear programming
// ---------------------------------------------------------------------------------------------------------------

const LP_EPS = 1e-9;

/** Pivot the dense tableau (rows T[0..m-1], objective row `z`) on (pr, pc). */
function pivot(T, z, basis, pr, pc, ncol) {
  const row = T[pr];
  const pv = row[pc];
  for (let j = 0; j <= ncol; j++) row[j] /= pv;
  row[pc] = 1;
  for (let i = 0; i < T.length; i++) {
    if (i === pr) continue;
    const Ti = T[i];
    const fct = Ti[pc];
    if (fct === 0) continue;
    for (let j = 0; j <= ncol; j++) Ti[j] -= fct * row[j];
    Ti[pc] = 0;
  }
  const fz = z[pc];
  if (fz !== 0) {
    for (let j = 0; j <= ncol; j++) z[j] -= fz * row[j];
    z[pc] = 0;
  }
  basis[pr] = pc;
}

/**
 * Primal simplex iterations with Bland's anti-cycling rule (smallest-index entering variable among those with a
 * negative reduced cost; smallest-index leaving variable among the minimum-ratio ties).
 * @returns {'optimal' | 'unbounded' | 'iteration-limit'}
 */
function simplexIterate(T, z, basis, ncol, nEnter, maxIter) {
  const m = T.length;
  for (let it = 0; it < maxIter; it++) {
    let pc = -1;
    for (let j = 0; j < nEnter; j++) {
      if (z[j] < -LP_EPS) {
        pc = j;
        break;
      }
    }
    if (pc < 0) return 'optimal';
    let pr = -1;
    let best = Infinity;
    for (let i = 0; i < m; i++) {
      const a = T[i][pc];
      if (a > LP_EPS) {
        const ratio = T[i][ncol] / a;
        if (ratio < best - 1e-12 || (Math.abs(ratio - best) <= 1e-12 && pr >= 0 && basis[i] < basis[pr])) {
          best = ratio;
          pr = i;
        }
      }
    }
    if (pr < 0) return 'unbounded';
    pivot(T, z, basis, pr, pc, ncol);
  }
  return 'iteration-limit';
}

/**
 * Linear programme: minimise cᵀx subject to Aub·x <= bub, Aeq·x = beq and lo <= x <= hi.
 * Method: dense two-phase primal simplex with Bland's rule. Variables are shifted to their finite lower bound
 * (mirrored when only the upper bound is finite, split when free); finite upper bounds become rows.
 * Default bounds are [0, +Infinity) for every variable (the usual LP convention).
 * Reference: Dantzig (1963), "Linear Programming and Extensions"; Bland (1977), Math. Oper. Res. 2, 103–107;
 * Chvátal (1983), "Linear Programming", ch. 8. Validity: dense problems up to a few hundred rows and columns.
 * @param {{ c: number[], Aub?: number[][], bub?: number[], Aeq?: number[][], beq?: number[],
 *   bounds?: [number, number][] }} problem objective coefficients; inequality and equality systems; bounds
 * @returns {{ x: number[] | null, f: number, status: 'optimal' | 'infeasible' | 'unbounded' }}
 *   x is null and f is +Infinity (infeasible) or -Infinity (unbounded) when no optimum exists
 */
export function linprog(problem) {
  if (!problem || !problem.c) throw new RangeError('linprog: c is required');
  const c = problem.c;
  const n = c.length;
  const Aub = problem.Aub || [];
  const bub = problem.bub || [];
  const Aeq = problem.Aeq || [];
  const beq = problem.beq || [];
  if (Aub.length !== bub.length) throw new RangeError('linprog: Aub and bub must have the same number of rows');
  if (Aeq.length !== beq.length) throw new RangeError('linprog: Aeq and beq must have the same number of rows');
  for (const row of Aub) if (row.length !== n) throw new RangeError('linprog: every Aub row must have c.length entries');
  for (const row of Aeq) if (row.length !== n) throw new RangeError('linprog: every Aeq row must have c.length entries');
  const bnd = problem.bounds;
  if (bnd && bnd.length !== n) throw new RangeError('linprog: bounds must have one [lo, hi] pair per variable');

  // x_j = off_j + sgn_j * y_col(j) [ - y_col2(j) ]
  const off = new Float64Array(n);
  const sgn = new Float64Array(n);
  const col = new Int32Array(n);
  const col2 = new Int32Array(n).fill(-1);
  const ubRows = []; // [column, upper value]
  let ny = 0;
  for (let j = 0; j < n; j++) {
    let lo = 0;
    let hi = Infinity;
    if (bnd && bnd[j]) {
      lo = bnd[j][0] === null || bnd[j][0] === undefined ? -Infinity : bnd[j][0];
      hi = bnd[j][1] === null || bnd[j][1] === undefined ? Infinity : bnd[j][1];
    }
    if (Number.isNaN(lo) || Number.isNaN(hi)) throw new RangeError(`linprog: bounds[${j}] must be numbers`);
    if (lo > hi + 1e-12) return { x: null, f: Infinity, status: 'infeasible' };
    if (Number.isFinite(lo)) {
      off[j] = lo;
      sgn[j] = 1;
      col[j] = ny++;
      if (Number.isFinite(hi)) ubRows.push([col[j], Math.max(hi - lo, 0)]);
    } else if (Number.isFinite(hi)) {
      off[j] = hi;
      sgn[j] = -1;
      col[j] = ny++;
    } else {
      off[j] = 0;
      sgn[j] = 1;
      col[j] = ny++;
      col2[j] = ny++;
    }
  }
  const mUb = Aub.length + ubRows.length;
  const m = mUb + Aeq.length;
  // rows in y-space
  const rows = new Array(m);
  const rhs = new Float64Array(m);
  const isLe = new Uint8Array(m);
  const fill = (src, bval, i, le) => {
    const r = new Float64Array(ny);
    let shift = 0;
    for (let j = 0; j < n; j++) {
      const a = src[j];
      if (a === 0) continue;
      shift += a * off[j];
      r[col[j]] += a * sgn[j];
      if (col2[j] >= 0) r[col2[j]] -= a;
    }
    rows[i] = r;
    rhs[i] = bval - shift;
    isLe[i] = le;
  };
  let ri = 0;
  for (let i = 0; i < Aub.length; i++) fill(Aub[i], bub[i], ri++, 1);
  for (const [cj, u] of ubRows) {
    const r = new Float64Array(ny);
    r[cj] = 1;
    rows[ri] = r;
    rhs[ri] = u;
    isLe[ri] = 1;
    ri++;
  }
  for (let i = 0; i < Aeq.length; i++) fill(Aeq[i], beq[i], ri++, 0);

  // columns: ny structural | mUb slacks | artificials
  const needArt = new Uint8Array(m);
  let nArt = 0;
  for (let i = 0; i < m; i++) {
    if (!isLe[i] || rhs[i] < 0) {
      needArt[i] = 1;
      nArt++;
    }
  }
  const nReal = ny + mUb;
  const ncol = nReal + nArt;
  const T = new Array(m);
  const basis = new Int32Array(m);
  let ai = 0;
  let si = 0;
  for (let i = 0; i < m; i++) {
    const t = new Float64Array(ncol + 1);
    const neg = rhs[i] < 0 ? -1 : 1;
    // row equilibration keeps the fixed pivot tolerance meaningful
    let scale = 0;
    for (let j = 0; j < ny; j++) scale = Math.max(scale, Math.abs(rows[i][j]));
    if (!(scale > 0)) scale = 1;
    const mult = neg / scale;
    for (let j = 0; j < ny; j++) t[j] = rows[i][j] * mult;
    t[ncol] = rhs[i] * mult;
    if (isLe[i]) {
      t[ny + si] = mult;
      if (!needArt[i]) basis[i] = ny + si;
      si++;
    }
    if (needArt[i]) {
      t[nReal + ai] = 1;
      basis[i] = nReal + ai;
      ai++;
    }
    T[i] = t;
  }
  const maxIter = 200 * (m + ncol) + 2000;
  if (nArt > 0) {
    // phase 1: minimise the sum of artificial variables
    const z = new Float64Array(ncol + 1);
    for (let i = 0; i < m; i++) {
      if (!needArt[i]) continue;
      for (let j = 0; j < nReal; j++) z[j] -= T[i][j];
      z[ncol] -= T[i][ncol];
    }
    const st = simplexIterate(T, z, basis, ncol, nReal, maxIter);
    let scaleB = 1;
    for (let i = 0; i < m; i++) scaleB = Math.max(scaleB, Math.abs(T[i][ncol]));
    if (st !== 'optimal' || -z[ncol] > 1e-7 * scaleB) return { x: null, f: Infinity, status: 'infeasible' };
    // drive remaining (degenerate) artificials out of the basis
    for (let i = 0; i < m; i++) {
      if (basis[i] < nReal) continue;
      let pc = -1;
      for (let j = 0; j < nReal; j++) {
        if (Math.abs(T[i][j]) > 1e-8) {
          pc = j;
          break;
        }
      }
      if (pc >= 0) pivot(T, z, basis, i, pc, ncol);
      else T[i][ncol] = 0; // redundant constraint row
    }
  }
  // phase 2
  const cy = new Float64Array(ncol);
  let cost0 = 0;
  for (let j = 0; j < n; j++) {
    cost0 += c[j] * off[j];
    cy[col[j]] += c[j] * sgn[j];
    if (col2[j] >= 0) cy[col2[j]] -= c[j];
  }
  const z = new Float64Array(ncol + 1);
  for (let j = 0; j < nReal; j++) z[j] = cy[j];
  for (let i = 0; i < m; i++) {
    const cb = basis[i] < nReal ? cy[basis[i]] : 0;
    if (cb === 0) continue;
    for (let j = 0; j <= ncol; j++) z[j] -= cb * T[i][j];
  }
  for (let i = 0; i < m; i++) if (basis[i] < nReal) z[basis[i]] = 0;
  const st = simplexIterate(T, z, basis, ncol, nReal, maxIter);
  if (st === 'unbounded') return { x: null, f: -Infinity, status: 'unbounded' };
  if (st !== 'optimal') throw new Error('linprog: simplex iteration limit reached');
  const y = new Float64Array(ncol);
  for (let i = 0; i < m; i++) y[basis[i]] = T[i][ncol];
  const x = new Array(n);
  let fval = 0;
  for (let j = 0; j < n; j++) {
    let v = off[j] + sgn[j] * y[col[j]];
    if (col2[j] >= 0) v -= y[col2[j]];
    x[j] = v;
    fval += c[j] * v;
  }
  return { x, f: fval, status: 'optimal' };
}

/**
 * Mixed-integer linear programme by depth-first branch and bound on the LP relaxation (branching on the most
 * fractional integer variable, the child nearer to the relaxed value explored first, bounding by the incumbent).
 * Reference: Land & Doig (1960), Econometrica 28, 497–520; Wolsey (1998), "Integer Programming", ch. 7.
 * Validity: small to medium problems (tens of integer variables); integer variables should have finite bounds
 * for guaranteed termination. Default bounds [0, +Infinity).
 * @param {{ c: number[], Aub?: number[][], bub?: number[], Aeq?: number[][], beq?: number[],
 *   bounds?: [number, number][], integer?: boolean[], maxNodes?: number }} problem as {@link linprog} plus the
 *   integrality flags; maxNodes caps the search (default 100000)
 * @returns {{ x: number[] | null, f: number, status: 'optimal' | 'infeasible' | 'unbounded' | 'node-limit',
 *   nodes: number }} best integer-feasible solution, its objective, status and the number of LP nodes solved
 */
export function milp(problem) {
  if (!problem || !problem.c) throw new RangeError('milp: c is required');
  const n = problem.c.length;
  const integer = problem.integer || new Array(n).fill(true);
  if (integer.length !== n) throw new RangeError('milp: integer must have one flag per variable');
  const maxNodes = problem.maxNodes ?? 100000;
  const base = [];
  for (let j = 0; j < n; j++) {
    const bj = problem.bounds && problem.bounds[j] ? problem.bounds[j] : [0, Infinity];
    let lo = bj[0] === null || bj[0] === undefined ? -Infinity : bj[0];
    let hi = bj[1] === null || bj[1] === undefined ? Infinity : bj[1];
    if (integer[j]) {
      lo = Math.ceil(lo - 1e-9);
      hi = Math.floor(hi + 1e-9);
    }
    base.push([lo, hi]);
  }
  const TOL = 1e-6;
  let best = null;
  let fbest = Infinity;
  let nodes = 0;
  let rootStatus = null;
  let limited = false;
  const stack = [base];
  while (stack.length) {
    if (nodes >= maxNodes) {
      limited = true;
      break;
    }
    const bounds = stack.pop();
    nodes++;
    const res = linprog({ c: problem.c, Aub: problem.Aub, bub: problem.bub, Aeq: problem.Aeq, beq: problem.beq, bounds });
    if (rootStatus === null) rootStatus = res.status;
    if (res.status === 'unbounded') {
      if (nodes === 1) return { x: null, f: -Infinity, status: 'unbounded', nodes };
      continue;
    }
    if (res.status !== 'optimal') continue;
    if (res.f >= fbest - 1e-9 * (1 + Math.abs(fbest))) continue;
    let jb = -1;
    let worst = TOL;
    for (let j = 0; j < n; j++) {
      if (!integer[j]) continue;
      const fr = Math.abs(res.x[j] - Math.round(res.x[j]));
      if (fr > worst) {
        worst = fr;
        jb = j;
      }
    }
    if (jb < 0) {
      best = res.x.map((v, j) => (integer[j] ? Math.round(v) : v));
      fbest = res.f;
      continue;
    }
    const v = res.x[jb];
    const down = bounds.map((p) => p.slice());
    const up = bounds.map((p) => p.slice());
    down[jb][1] = Math.floor(v);
    up[jb][0] = Math.ceil(v);
    // push the farther child first so the nearer one is explored next
    if (v - Math.floor(v) < 0.5) stack.push(up, down);
    else stack.push(down, up);
  }
  if (best === null) return { x: null, f: Infinity, status: limited ? 'node-limit' : 'infeasible', nodes };
  let f = 0;
  for (let j = 0; j < n; j++) f += problem.c[j] * best[j];
  return { x: best, f, status: limited ? 'node-limit' : 'optimal', nodes };
}
