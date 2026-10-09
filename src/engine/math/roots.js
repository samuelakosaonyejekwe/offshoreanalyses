/**
 * Root finding for scalar equations, nonlinear systems and polynomials.
 *
 * All routines are unit-agnostic: tolerances on x are in the units of x, residual tolerances in the units of f.
 *
 * References: Brent, "Algorithms for Minimization without Derivatives" (1973), ch. 4 (zeroin);
 * Dennis & Schnabel, "Numerical Methods for Unconstrained Optimization and Nonlinear Equations" (1983), ch. 6
 * (damped Newton with backtracking); Aberth, Math. Comp. 27 (1973) 339–344 (simultaneous polynomial roots).
 */

import { solveDense, pinv, matVec } from './linalg.js';

const EPS = 2.220446049250313e-16;

function checkFn(f, name) {
  if (typeof f !== 'function') throw new RangeError(`${name} must be a function`);
}

function evalFinite(f, x, name) {
  const v = f(x);
  if (Number.isNaN(v)) throw new RangeError(`${name}: f(${x}) is NaN`);
  return v;
}

/**
 * Root of f on a bracketing interval by Brent's method (inverse quadratic interpolation / secant / bisection).
 * Reference: Brent (1973) ch. 4. Validity: f continuous on [a, b] with f(a)·f(b) <= 0; throws RangeError if the
 * interval does not bracket a sign change.
 * @param {(x: number) => number} f function
 * @param {number} a interval end [x units]
 * @param {number} b other interval end [x units]
 * @param {{ tol?: number, maxIter?: number }} [opts] tol: absolute tolerance on the root [x units] (default 1e-12);
 *   maxIter: iteration cap (default 200)
 * @returns {number} root [x units]
 */
export function brent(f, a, b, { tol = 1e-12, maxIter = 200 } = {}) {
  checkFn(f, 'brent: f');
  if (!Number.isFinite(a) || !Number.isFinite(b)) throw new RangeError('brent: a and b must be finite');
  let fa = evalFinite(f, a, 'brent');
  let fb = evalFinite(f, b, 'brent');
  if (fa === 0) return a;
  if (fb === 0) return b;
  if (fa > 0 === fb > 0) {
    throw new RangeError(`brent: no sign change on [a, b] (f(${a}) = ${fa}, f(${b}) = ${fb})`);
  }
  let c = a;
  let fc = fa;
  let d = b - a;
  let e = d;
  for (let iter = 0; iter < maxIter; iter++) {
    if (fb > 0 === fc > 0) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }
    if (Math.abs(fc) < Math.abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const tol1 = 2 * EPS * Math.abs(b) + 0.5 * tol;
    const xm = 0.5 * (c - b);
    if (Math.abs(xm) <= tol1 || fb === 0) return b;
    if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
      const s = fb / fa;
      let p;
      let q;
      if (a === c) {
        p = 2 * xm * s;
        q = 1 - s;
      } else {
        const qq = fa / fc;
        const r = fb / fc;
        p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1));
        q = (qq - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) q = -q;
      p = Math.abs(p);
      const min1 = 3 * xm * q - Math.abs(tol1 * q);
      const min2 = Math.abs(e * q);
      if (2 * p < Math.min(min1, min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (Math.abs(d) > tol1) b += d;
    else b += xm >= 0 ? tol1 : -tol1;
    fb = evalFinite(f, b, 'brent');
  }
  throw new Error(`brent: no convergence in ${maxIter} iterations`);
}

/**
 * Root of f on a bracketing interval by bisection.
 * Validity: f(a)·f(b) <= 0 (throws RangeError otherwise); f need not be continuous (returns the sign-change point).
 * @param {(x: number) => number} f function
 * @param {number} a interval end [x units]
 * @param {number} b other interval end [x units]
 * @param {{ tol?: number, maxIter?: number }} [opts] tol: absolute tolerance on the root (default 1e-12);
 *   maxIter: iteration cap (default 200)
 * @returns {number} root [x units]
 */
export function bisect(f, a, b, { tol = 1e-12, maxIter = 200 } = {}) {
  checkFn(f, 'bisect: f');
  if (!Number.isFinite(a) || !Number.isFinite(b)) throw new RangeError('bisect: a and b must be finite');
  let fa = evalFinite(f, a, 'bisect');
  const fb = evalFinite(f, b, 'bisect');
  if (fa === 0) return a;
  if (fb === 0) return b;
  if (fa > 0 === fb > 0) {
    throw new RangeError(`bisect: no sign change on [a, b] (f(${a}) = ${fa}, f(${b}) = ${fb})`);
  }
  let lo = a;
  let hi = b;
  for (let iter = 0; iter < maxIter; iter++) {
    const mid = lo + 0.5 * (hi - lo);
    if (Math.abs(hi - lo) <= tol + 2 * EPS * Math.abs(mid) || mid === lo || mid === hi) return mid;
    const fm = evalFinite(f, mid, 'bisect');
    if (fm === 0) return mid;
    if (fm > 0 === fa > 0) {
      lo = mid;
      fa = fm;
    } else {
      hi = mid;
    }
  }
  return lo + 0.5 * (hi - lo);
}

/**
 * Newton–Raphson iteration for a scalar equation, with a central-difference derivative when df is omitted.
 * Validity: f smooth near the root and x0 inside the basin of attraction; throws Error on a zero derivative or
 * when maxIter is exhausted (use brent with a bracket for guaranteed convergence).
 * @param {(x: number) => number} f function
 * @param {number} x0 starting point [x units]
 * @param {{ df?: (x: number) => number, tol?: number, maxIter?: number }} [opts] df: derivative; tol: relative
 *   step tolerance |dx| <= tol·(1 + |x|) (default 1e-12); maxIter (default 50)
 * @returns {number} root [x units]
 */
export function newton(f, x0, { df, tol = 1e-12, maxIter = 50 } = {}) {
  checkFn(f, 'newton: f');
  if (!Number.isFinite(x0)) throw new RangeError('newton: x0 must be finite');
  let x = x0;
  for (let iter = 0; iter < maxIter; iter++) {
    const fx = f(x);
    if (fx === 0) return x;
    if (!Number.isFinite(fx)) throw new Error(`newton: f(${x}) is not finite`);
    let d;
    if (df) {
      d = df(x);
    } else {
      const h = 6.0554544523933395e-6 * Math.max(Math.abs(x), 1); // eps^(1/3)
      d = (f(x + h) - f(x - h)) / (2 * h);
    }
    if (d === 0 || !Number.isFinite(d)) throw new Error(`newton: derivative is zero or not finite at x = ${x}`);
    const dx = fx / d;
    x -= dx;
    if (Math.abs(dx) <= tol * (1 + Math.abs(x))) return x;
  }
  throw new Error(`newton: no convergence in ${maxIter} iterations`);
}

/**
 * Expand an interval geometrically around x0 until it brackets a sign change of f.
 * Reference: Press et al., "Numerical Recipes" §9.1 (zbrac), with optional hard limits.
 * Validity: f defined on the searched range; throws Error if no sign change is found within maxIter expansions.
 * @param {(x: number) => number} f function
 * @param {number} x0 starting point [x units]
 * @param {{ factor?: number, maxIter?: number, min?: number, max?: number, step?: number }} [opts] factor: growth
 *   factor per expansion (default 1.6); maxIter (default 60); min/max: hard limits of the search range; step:
 *   initial half-width (default 5 % of max(|x0|, 1))
 * @returns {[number, number]} [a, b] with a < b and f(a)·f(b) <= 0
 */
export function findBracket(f, x0, { factor = 1.6, maxIter = 60, min = -Infinity, max = Infinity, step } = {}) {
  checkFn(f, 'findBracket: f');
  if (!Number.isFinite(x0)) throw new RangeError('findBracket: x0 must be finite');
  if (!(factor > 1)) throw new RangeError('findBracket: factor must be > 1');
  if (!(min < max)) throw new RangeError('findBracket: min must be < max');
  const x = Math.min(Math.max(x0, min), max);
  const h = step !== undefined ? Math.abs(step) : 0.05 * Math.max(Math.abs(x), 1);
  if (!(h > 0)) throw new RangeError('findBracket: step must be non-zero');
  let a = Math.max(x - h, min);
  let b = Math.min(x + h, max);
  let fa = evalFinite(f, a, 'findBracket');
  let fb = evalFinite(f, b, 'findBracket');
  for (let iter = 0; iter < maxIter; iter++) {
    if (fa === 0 || fb === 0 || fa > 0 !== fb > 0) return [a, b];
    const canA = a > min;
    const canB = b < max;
    if (!canA && !canB) break;
    const w = b - a;
    if (canA && (!canB || Math.abs(fa) < Math.abs(fb))) {
      a = Math.max(a - factor * w, min);
      fa = evalFinite(f, a, 'findBracket');
    } else {
      b = Math.min(b + factor * w, max);
      fb = evalFinite(f, b, 'findBracket');
    }
  }
  if (fa === 0 || fb === 0 || fa > 0 !== fb > 0) return [a, b];
  throw new Error(`findBracket: no sign change found in [${a}, ${b}]`);
}

function norm2(v) {
  let s = 0;
  for (let i = 0; i < v.length; i++) s += v[i] * v[i];
  return Math.sqrt(s);
}

/**
 * Solve the nonlinear system F(x) = 0 by Newton's method with a backtracking (Armijo) line search on ½‖F‖².
 * The Jacobian is analytic when `jac` is given, otherwise forward finite differences. A singular Jacobian falls
 * back to a pseudo-inverse (minimum-norm) step. Reference: Dennis & Schnabel (1983) §6.3.
 * Validity: F continuously differentiable, square system. Does not throw on non-convergence: check `converged`.
 * @param {(x: number[]) => ArrayLike<number>} F residual function, n -> n
 * @param {ArrayLike<number>} x0 starting point, length n
 * @param {{ jac?: (x: number[]) => number[][], tol?: number, maxIter?: number, lineSearch?: boolean,
 *   fdStep?: number }} [opts] tol: converged when ‖F‖₂ <= tol·max(1, ‖F(x0)‖₂) or when a full Newton step is
 *   smaller than tol·(1 + ‖x‖₂) (default 1e-10); maxIter (default 50); lineSearch (default true); fdStep: relative
 *   finite-difference step (default sqrt(eps))
 * @returns {{ x: number[], iterations: number, residual: number, converged: boolean }} x: solution; residual:
 *   final ‖F(x)‖₂ [units of F]
 */
export function newtonSystem(F, x0, { jac, tol = 1e-10, maxIter = 50, lineSearch = true, fdStep } = {}) {
  checkFn(F, 'newtonSystem: F');
  if (!x0 || !(x0.length >= 1)) throw new RangeError('newtonSystem: x0 must be a non-empty vector');
  const n = x0.length;
  let x = Array.from(x0);
  let Fx = Array.from(F(x.slice()));
  if (Fx.length !== n) throw new RangeError('newtonSystem: F must return a vector of the same length as x0');
  let res = norm2(Fx);
  if (!Number.isFinite(res)) throw new RangeError('newtonSystem: F(x0) is not finite');
  const target = tol * Math.max(1, res);
  const rel = fdStep !== undefined ? fdStep : Math.sqrt(EPS);
  let iterations = 0;
  let converged = res <= target;
  while (!converged && iterations < maxIter) {
    iterations++;
    let J;
    if (jac) {
      J = jac(x.slice());
    } else {
      J = new Array(n);
      for (let i = 0; i < n; i++) J[i] = new Array(n);
      const xp = x.slice();
      for (let j = 0; j < n; j++) {
        const h = rel * Math.max(Math.abs(x[j]), 1);
        xp[j] = x[j] + h;
        const Fp = F(xp.slice());
        xp[j] = x[j];
        for (let i = 0; i < n; i++) J[i][j] = (Fp[i] - Fx[i]) / h;
      }
    }
    let dx;
    try {
      dx = solveDense(J, Fx);
    } catch (err) {
      dx = matVec(pinv(J), Fx);
    }
    for (let i = 0; i < n; i++) dx[i] = -dx[i];
    let lam = 1;
    let xn;
    let Fn;
    let resN;
    for (let ls = 0; ; ls++) {
      xn = new Array(n);
      for (let i = 0; i < n; i++) xn[i] = x[i] + lam * dx[i];
      Fn = Array.from(F(xn.slice()));
      resN = norm2(Fn);
      // Armijo on ½‖F‖² along the Newton direction: slope = −‖F‖²
      if (!lineSearch) break;
      if (Number.isFinite(resN) && resN * resN <= res * res * (1 - 2e-4 * lam)) break;
      if (ls >= 40) break;
      lam *= 0.5;
    }
    if (!Number.isFinite(resN)) break;
    const stepNorm = lam * norm2(dx);
    x = xn;
    Fx = Fn;
    res = resN;
    if (res <= target) converged = true;
    else if (lam === 1 && stepNorm <= tol * (1 + norm2(x))) converged = true;
    else if (stepNorm === 0) break;
  }
  return { x, iterations, residual: res, converged };
}

/**
 * Fixed-point (Picard) iteration x ← x + relax·(g(x) − x) for a scalar or a vector unknown.
 * Validity: g a contraction near the solution (|g'| < 1), or made so by under-relaxation 0 < relax < 1.
 * Does not throw on non-convergence: check `converged`.
 * @param {(x: any) => any} g iteration map (number -> number, or number[] -> number[])
 * @param {number | ArrayLike<number>} x0 starting value
 * @param {{ tol?: number, maxIter?: number, relax?: number }} [opts] tol: converged when ‖g(x) − x‖₂ <=
 *   tol·(1 + ‖x‖₂) (default 1e-10); maxIter (default 200); relax: relaxation factor (default 1)
 * @returns {{ x: number | number[], iterations: number, converged: boolean }} x has the type of x0
 */
export function fixedPoint(g, x0, { tol = 1e-10, maxIter = 200, relax = 1 } = {}) {
  checkFn(g, 'fixedPoint: g');
  if (!(relax > 0)) throw new RangeError('fixedPoint: relax must be > 0');
  const scalar = typeof x0 === 'number';
  let x = scalar ? [x0] : Array.from(x0);
  const n = x.length;
  for (let iter = 1; iter <= maxIter; iter++) {
    const gx = scalar ? [g(x[0])] : g(x.slice());
    let d2 = 0;
    let x2 = 0;
    const xn = new Array(n);
    for (let i = 0; i < n; i++) {
      const d = gx[i] - x[i];
      if (!Number.isFinite(d)) throw new Error('fixedPoint: iteration produced a non-finite value');
      d2 += d * d;
      xn[i] = x[i] + relax * d;
      x2 += xn[i] * xn[i];
    }
    x = xn;
    if (Math.sqrt(d2) <= tol * (1 + Math.sqrt(x2))) {
      return { x: scalar ? x[0] : x, iterations: iter, converged: true };
    }
  }
  return { x: scalar ? x[0] : x, iterations: maxIter, converged: false };
}

/**
 * All complex roots of a real-coefficient polynomial by the Aberth–Ehrlich simultaneous iteration.
 * coeffs = [c0, c1, ..., cn] represents c0·xⁿ + c1·xⁿ⁻¹ + ... + cn (highest degree first). Leading zeros are
 * dropped; exact zero roots are deflated first. Roots are returned sorted by real part, then imaginary part;
 * imaginary parts below ~6e-8·|z| are snapped to zero (resolution limit of a double real root).
 * Validity: degree >= 1 after trimming; multiple roots of multiplicity m are accurate to about eps^(1/m).
 * @param {ArrayLike<number>} coeffs polynomial coefficients, highest degree first
 * @returns {{ re: number[], im: number[] }} real and imaginary parts of the roots
 */
export function polyRoots(coeffs) {
  if (!coeffs || !(coeffs.length >= 1)) throw new RangeError('polyRoots: coeffs must be a non-empty array');
  const c = Array.from(coeffs);
  for (const v of c) if (!Number.isFinite(v)) throw new RangeError('polyRoots: coeffs must be finite');
  while (c.length > 0 && c[0] === 0) c.shift();
  if (c.length < 2) throw new RangeError('polyRoots: polynomial must have degree >= 1');
  let nZero = 0;
  while (c.length > 1 && c[c.length - 1] === 0) {
    c.pop();
    nZero++;
  }
  const n = c.length - 1;
  const zr = new Float64Array(n);
  const zi = new Float64Array(n);
  if (n === 1) {
    zr[0] = -c[1] / c[0];
  } else if (n >= 2) {
    // normalise to monic and scale the variable by the geometric-mean root modulus to avoid overflow
    const lead = c[0];
    const rho = Math.pow(Math.abs(c[n] / lead), 1 / n);
    const a = new Float64Array(n + 1);
    let sc = 1;
    for (let k = 0; k <= n; k++) {
      a[k] = c[k] / lead / sc; // p(rho·w)/rhoⁿ has coefficients c_k / rho^k
      sc *= rho;
    }
    // starting points: circle of radius from the Cauchy bound of the scaled polynomial, slightly rotated
    let bound = 0;
    for (let k = 1; k <= n; k++) bound = Math.max(bound, Math.pow(Math.abs(a[k]), 1 / k));
    const r0 = Math.min(Math.max(bound, 0.5), 2);
    for (let k = 0; k < n; k++) {
      const th = (2 * Math.PI * k) / n + 0.4;
      zr[k] = r0 * Math.cos(th);
      zi[k] = r0 * Math.sin(th);
    }
    const maxIter = 500;
    for (let iter = 0; iter < maxIter; iter++) {
      let moved = 0;
      for (let k = 0; k < n; k++) {
        const x = zr[k];
        const y = zi[k];
        // Horner for p and p'
        let pr = a[0];
        let pi = 0;
        let dr = 0;
        let di = 0;
        for (let j = 1; j <= n; j++) {
          const tr = dr * x - di * y + pr;
          const ti = dr * y + di * x + pi;
          dr = tr;
          di = ti;
          const qr = pr * x - pi * y + a[j];
          const qi = pr * y + pi * x;
          pr = qr;
          pi = qi;
        }
        if (pr === 0 && pi === 0) continue;
        // u = p / p'
        let ur;
        let ui;
        const dm = dr * dr + di * di;
        if (dm === 0) {
          ur = 1e-3 * (1 + Math.abs(x));
          ui = 1e-3 * (1 + Math.abs(y));
        } else {
          ur = (pr * dr + pi * di) / dm;
          ui = (pi * dr - pr * di) / dm;
        }
        // s = Σ_{j≠k} 1 / (z_k − z_j)
        let sr = 0;
        let si = 0;
        for (let j = 0; j < n; j++) {
          if (j === k) continue;
          const ex = x - zr[j];
          const ey = y - zi[j];
          const em = ex * ex + ey * ey;
          if (em === 0) continue;
          sr += ex / em;
          si -= ey / em;
        }
        // w = u / (1 − u·s)
        const gr = 1 - (ur * sr - ui * si);
        const gi = -(ur * si + ui * sr);
        const gm = gr * gr + gi * gi;
        let wr = ur;
        let wi = ui;
        if (gm > 0) {
          wr = (ur * gr + ui * gi) / gm;
          wi = (ui * gr - ur * gi) / gm;
        }
        zr[k] = x - wr;
        zi[k] = y - wi;
        const rel = Math.hypot(wr, wi) / Math.max(Math.hypot(zr[k], zi[k]), 1e-300);
        if (rel > moved) moved = rel;
      }
      if (moved <= 4 * EPS) break;
    }
    for (let k = 0; k < n; k++) {
      zr[k] *= rho;
      zi[k] *= rho;
    }
  }
  const roots = [];
  for (let k = 0; k < n; k++) {
    let im = zi[k];
    if (Math.abs(im) <= 6e-8 * Math.hypot(zr[k], zi[k])) im = 0;
    roots.push([zr[k], im]);
  }
  for (let k = 0; k < nZero; k++) roots.push([0, 0]);
  roots.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
  return { re: roots.map((r) => r[0]), im: roots.map((r) => r[1]) };
}
