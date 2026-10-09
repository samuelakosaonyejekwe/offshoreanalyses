/**
 * Numerical quadrature: sampled-data rules, Gauss–Legendre nodes and adaptive Gauss–Kronrod integration.
 *
 * Unit-agnostic: the result carries the units of y times the units of x.
 *
 * References: Davis & Rabinowitz, "Methods of Numerical Integration", 2nd ed. (trapezoid/Simpson §2.1–2.2,
 * Gauss–Legendre §2.7); Piessens et al., "QUADPACK" (1983) (G7–K15 pair, QAG global adaptive strategy and error
 * estimate).
 */

const EPS = 2.220446049250313e-16;

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

function checkSamples(y, x, name) {
  if (!isVec(y)) throw new RangeError(`${name}: y must be an array`);
  if (isVec(x)) {
    if (x.length !== y.length) throw new RangeError(`${name}: x and y must have the same length`);
    return true;
  }
  if (x !== undefined && (typeof x !== 'number' || !Number.isFinite(x))) {
    throw new RangeError(`${name}: second argument must be an abscissa array or a finite spacing dx`);
  }
  return false;
}

/**
 * Composite trapezoidal rule for sampled data (second order).
 * Validity: any ordered abscissae (non-uniform allowed); fewer than two samples integrate to 0.
 * @param {ArrayLike<number>} y samples [y units]
 * @param {ArrayLike<number> | number} [x] abscissae (same length as y) or uniform spacing dx (default 1) [x units]
 * @returns {number} ∫ y dx [y units · x units]
 */
export function trapz(y, x = 1) {
  const nonUniform = checkSamples(y, x, 'trapz');
  const n = y.length;
  if (n < 2) return 0;
  let s = 0;
  if (nonUniform) {
    for (let i = 1; i < n; i++) s += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    return s;
  }
  for (let i = 1; i < n - 1; i++) s += y[i];
  return x * (s + 0.5 * (y[0] + y[n - 1]));
}

/**
 * Cumulative trapezoidal integral; out[0] = 0 and out[k] = ∫ from x[0] to x[k].
 * @param {ArrayLike<number>} y samples [y units]
 * @param {ArrayLike<number> | number} [x] abscissae (same length as y) or uniform spacing dx (default 1) [x units]
 * @returns {Float64Array} running integral, same length as y [y units · x units]
 */
export function cumtrapz(y, x = 1) {
  const nonUniform = checkSamples(y, x, 'cumtrapz');
  const n = y.length;
  const out = new Float64Array(n);
  for (let i = 1; i < n; i++) {
    const h = nonUniform ? x[i] - x[i - 1] : x;
    out[i] = out[i - 1] + 0.5 * h * (y[i] + y[i - 1]);
  }
  return out;
}

/**
 * Composite Simpson rule for sampled data (fourth order).
 * Uniform spacing: Simpson 1/3 panels; when the number of intervals is odd the last three intervals use the
 * Simpson 3/8 rule so the order is preserved. If an abscissa array is passed instead of dx, the non-uniform
 * (three-point parabola) Simpson rule is used, with the standard end correction for an odd number of intervals.
 * Validity: at least two samples (two samples fall back to the trapezoid; fewer integrate to 0).
 * @param {ArrayLike<number>} y samples [y units]
 * @param {number | ArrayLike<number>} [dx] uniform spacing (default 1) or abscissae [x units]
 * @returns {number} ∫ y dx [y units · x units]
 */
export function simpson(y, dx = 1) {
  const nonUniform = checkSamples(y, dx, 'simpson');
  const n = y.length;
  if (n < 2) return 0;
  if (n === 2) return 0.5 * (nonUniform ? dx[1] - dx[0] : dx) * (y[0] + y[1]);
  const m = n - 1; // intervals
  if (!nonUniform) {
    let s = 0;
    const mEven = m % 2 === 0 ? m : m - 3; // intervals handled by the 1/3 rule
    if (mEven > 0) {
      let odd = 0;
      let even = 0;
      for (let i = 1; i < mEven; i += 2) odd += y[i];
      for (let i = 2; i < mEven; i += 2) even += y[i];
      s += (dx / 3) * (y[0] + 4 * odd + 2 * even + y[mEven]);
    }
    if (mEven !== m) {
      s += (3 * dx / 8) * (y[m - 3] + 3 * y[m - 2] + 3 * y[m - 1] + y[m]);
    }
    return s;
  }
  const x = dx;
  let s = 0;
  const last = m % 2 === 0 ? m : m - 1;
  for (let i = 0; i < last; i += 2) {
    const h0 = x[i + 1] - x[i];
    const h1 = x[i + 2] - x[i + 1];
    if (h0 === 0 || h1 === 0) {
      s += 0.5 * h0 * (y[i] + y[i + 1]) + 0.5 * h1 * (y[i + 1] + y[i + 2]);
      continue;
    }
    const hs = h0 + h1;
    s += (hs / 6) * ((2 - h1 / h0) * y[i] + ((hs * hs) / (h0 * h1)) * y[i + 1] + (2 - h0 / h1) * y[i + 2]);
  }
  if (last !== m) {
    const h1 = x[m] - x[m - 1];
    const h0 = x[m - 1] - x[m - 2];
    if (h0 === 0 || h1 === 0) {
      s += 0.5 * h1 * (y[m] + y[m - 1]);
    } else {
      const alpha = (2 * h1 * h1 + 3 * h1 * h0) / (6 * (h0 + h1));
      const beta = (h1 * h1 + 3 * h1 * h0) / (6 * h0);
      const eta = (h1 * h1 * h1) / (6 * h0 * (h0 + h1));
      s += alpha * y[m] + beta * y[m - 1] - eta * y[m - 2];
    }
  }
  return s;
}

const glCache = new Map();

/**
 * Gauss–Legendre nodes and weights on [−1, 1] (exact for polynomials of degree <= 2n − 1).
 * Nodes are the roots of the Legendre polynomial Pₙ found by Newton iteration on the three-term recurrence;
 * weights w = 2 / ((1 − x²) Pₙ'(x)²). Reference: Davis & Rabinowitz §2.7.
 * Validity: integer 1 <= n <= 64 as specified; larger n (up to 1024) is accepted and remains accurate.
 * @param {number} n number of points [-]
 * @returns {{ x: number[], w: number[] }} nodes (ascending) and weights (Σw = 2)
 */
export function gaussLegendre(n) {
  if (!Number.isInteger(n) || n < 1 || n > 1024) throw new RangeError('gaussLegendre: n must be an integer in 1..1024');
  const hit = glCache.get(n);
  if (hit) return { x: hit.x.slice(), w: hit.w.slice() };
  const x = new Array(n);
  const w = new Array(n);
  const half = (n + 1) >> 1;
  for (let i = 0; i < half; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let pp = 1;
    for (let it = 0; it < 100; it++) {
      let p1 = 1;
      let p2 = 0;
      for (let j = 1; j <= n; j++) {
        const p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      pp = (n * (z * p1 - p2)) / (z * z - 1);
      const dz = p1 / pp;
      z -= dz;
      if (Math.abs(dz) <= 4 * EPS) break;
    }
    // re-evaluate the derivative at the converged node for the weight
    let p1 = 1;
    let p2 = 0;
    for (let j = 1; j <= n; j++) {
      const p3 = p2;
      p2 = p1;
      p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
    }
    pp = (n * (z * p1 - p2)) / (z * z - 1);
    const wi = 2 / ((1 - z * z) * pp * pp);
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = wi;
    w[n - 1 - i] = wi;
  }
  if (n % 2 === 1) x[half - 1] = 0;
  glCache.set(n, { x: x.slice(), w: w.slice() });
  return { x, w };
}

// Gauss(7)–Kronrod(15) pair, QUADPACK qk15. XGK[1], XGK[3], XGK[5], XGK[7] are the 7-point Gauss nodes.
const XGK = [
  0.991455371120812639206854697526329, 0.949107912342758524526189684047851, 0.864864423359769072789712788640926,
  0.741531185599394439863864773280788, 0.58608723546769113029414483825873, 0.405845151377397166906606412076961,
  0.207784955007898467600689403773245, 0,
];
const WGK = [
  0.02293532201052922496373200805897, 0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
  0.140653259715525918745189590510238, 0.16900472663926790282658342659855, 0.190350578064785409913256402421014,
  0.204432940075298892414161999234649, 0.209482141084727828012999174891714,
];
const WG = [
  0.129484966168869693270611432679082, 0.27970539148927666790146777142378, 0.381830050505118944950369775488975,
  0.417959183673469387755102040816327,
];

const k15out = { result: 0, err: 0, resabs: 0 };

/**
 * One G7–K15 evaluation on [a, b]; writes result, error estimate and ∫|f| into k15out.
 * fv1/fv2 are caller-owned scratch buffers so that nested integrals (f calling integrate) stay re-entrant.
 */
function k15(g, a, b, fv1, fv2) {
  const c = 0.5 * (a + b);
  const h = 0.5 * (b - a);
  const fc = g(c);
  let resg = WG[3] * fc;
  let resk = WGK[7] * fc;
  let resabs = Math.abs(resk);
  for (let j = 0; j < 7; j++) {
    const dx = h * XGK[j];
    const f1 = g(c - dx);
    const f2 = g(c + dx);
    fv1[j] = f1;
    fv2[j] = f2;
    resk += WGK[j] * (f1 + f2);
    resabs += WGK[j] * (Math.abs(f1) + Math.abs(f2));
    if (j % 2 === 1) resg += WG[(j - 1) >> 1] * (f1 + f2);
  }
  const mean = 0.5 * resk;
  let resasc = WGK[7] * Math.abs(fc - mean);
  for (let j = 0; j < 7; j++) resasc += WGK[j] * (Math.abs(fv1[j] - mean) + Math.abs(fv2[j] - mean));
  const ah = Math.abs(h);
  resabs *= ah;
  resasc *= ah;
  let err = Math.abs((resk - resg) * h);
  if (resasc !== 0 && err !== 0) err = resasc * Math.min(1, Math.pow((200 * err) / resasc, 1.5));
  if (resabs > Number.MIN_VALUE / (50 * EPS)) err = Math.max(50 * EPS * resabs, err);
  k15out.result = resk * h;
  k15out.err = err;
  k15out.resabs = resabs;
}

/**
 * Definite integral of f over [a, b] by globally adaptive Gauss–Kronrod (G7–K15) quadrature (QUADPACK QAG
 * strategy: always bisect the sub-interval with the largest error estimate).
 * Infinite limits are mapped to a finite interval (x = a + t/(1−t) for a semi-infinite range, x = t/(1−t²) for the
 * whole line). Integrable end-point singularities are handled by the bisection (end points are never evaluated).
 * Validity: f finite at interior points; returns the best estimate if the tolerance cannot be reached within
 * maxDepth bisections of any sub-interval (no error is thrown for that case).
 * @param {(x: number) => number} f integrand
 * @param {number} a lower limit (may be −Infinity) [x units]
 * @param {number} b upper limit (may be +Infinity) [x units]
 * @param {{ tol?: number, maxDepth?: number, absTol?: number }} [opts] tol: relative tolerance on the integral
 *   (default 1e-9); maxDepth: maximum bisection depth of a sub-interval (default 40); absTol: additional absolute
 *   tolerance (default 0)
 * @returns {number} ∫ₐᵇ f(x) dx [f units · x units]
 */
export function integrate(f, a, b, { tol = 1e-9, maxDepth = 40, absTol = 0 } = {}) {
  if (typeof f !== 'function') throw new RangeError('integrate: f must be a function');
  if (Number.isNaN(a) || Number.isNaN(b)) throw new RangeError('integrate: a and b must not be NaN');
  if (!(tol > 0)) throw new RangeError('integrate: tol must be > 0');
  if (a === b) return 0;
  if (a > b) return -integrate(f, b, a, { tol, maxDepth, absTol });
  let g;
  let lo;
  let hi;
  const check = (x, v) => {
    if (Number.isNaN(v)) throw new RangeError(`integrate: f(${x}) is NaN`);
    return v;
  };
  if (Number.isFinite(a) && Number.isFinite(b)) {
    g = (x) => check(x, f(x));
    lo = a;
    hi = b;
  } else if (Number.isFinite(a)) {
    g = (t) => {
      const u = 1 - t;
      const v = f(a + t / u) / (u * u);
      return Number.isFinite(v) ? v : 0;
    };
    lo = 0;
    hi = 1;
  } else if (Number.isFinite(b)) {
    g = (t) => {
      const u = 1 - t;
      const v = f(b - t / u) / (u * u);
      return Number.isFinite(v) ? v : 0;
    };
    lo = 0;
    hi = 1;
  } else {
    g = (t) => {
      const u = 1 - t * t;
      const v = (f(t / u) * (1 + t * t)) / (u * u);
      return Number.isFinite(v) ? v : 0;
    };
    lo = -1;
    hi = 1;
  }
  const fv1 = new Float64Array(7);
  const fv2 = new Float64Array(7);
  const maxIntervals = 20000;
  let cap = 64;
  let ia = new Float64Array(cap);
  let ib = new Float64Array(cap);
  let ir = new Float64Array(cap);
  let ie = new Float64Array(cap);
  let id = new Int32Array(cap);
  k15(g, lo, hi, fv1, fv2);
  ia[0] = lo;
  ib[0] = hi;
  ir[0] = k15out.result;
  ie[0] = k15out.err;
  id[0] = 0;
  let count = 1;
  let total = ir[0];
  let errSum = ie[0];
  let resabs = k15out.resabs;
  for (;;) {
    const target = Math.max(tol * Math.abs(total), absTol, 50 * EPS * resabs);
    if (errSum <= target || count >= maxIntervals) break;
    // interval with the largest error that may still be split
    let k = -1;
    let emax = 0;
    for (let i = 0; i < count; i++) {
      if (ie[i] > emax && id[i] < maxDepth) {
        emax = ie[i];
        k = i;
      }
    }
    if (k < 0) break;
    const a1 = ia[k];
    const b2 = ib[k];
    const mid = 0.5 * (a1 + b2);
    if (!(mid > a1 && mid < b2)) {
      id[k] = maxDepth;
      continue;
    }
    if (count === cap) {
      cap *= 2;
      const grow = (arr, Ctor) => {
        const nu = new Ctor(cap);
        nu.set(arr);
        return nu;
      };
      ia = grow(ia, Float64Array);
      ib = grow(ib, Float64Array);
      ir = grow(ir, Float64Array);
      ie = grow(ie, Float64Array);
      id = grow(id, Int32Array);
    }
    const depth = id[k] + 1;
    k15(g, a1, mid, fv1, fv2);
    const r1 = k15out.result;
    const e1 = k15out.err;
    k15(g, mid, b2, fv1, fv2);
    const r2 = k15out.result;
    const e2 = k15out.err;
    ib[k] = mid;
    ir[k] = r1;
    ie[k] = e1;
    id[k] = depth;
    ia[count] = mid;
    ib[count] = b2;
    ir[count] = r2;
    ie[count] = e2;
    id[count] = depth;
    count++;
    // re-sum to avoid drift from repeated updates
    total = 0;
    errSum = 0;
    for (let i = 0; i < count; i++) {
      total += ir[i];
      errSum += ie[i];
    }
  }
  return total;
}
