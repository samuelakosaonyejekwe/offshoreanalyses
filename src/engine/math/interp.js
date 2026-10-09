/**
 * Interpolation, grids and resampling for tabulated data.
 *
 * Unit-agnostic: results carry the units of ys; abscissae must be in consistent units and strictly ascending.
 *
 * References: de Boor, "A Practical Guide to Splines" (2001), ch. IV (cubic spline end conditions);
 * Fritsch & Carlson, SIAM J. Numer. Anal. 17 (1980) 238–246 and Fritsch & Butland, SIAM J. Sci. Stat. Comput. 5
 * (1984) 300–304 (monotone piecewise cubic Hermite, PCHIP).
 */

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

function checkTable(xs, ys, name, minPoints = 1) {
  if (!isVec(xs) || !isVec(ys)) throw new RangeError(`${name}: xs and ys must be arrays`);
  if (xs.length !== ys.length) throw new RangeError(`${name}: xs and ys must have the same length`);
  if (xs.length < minPoints) throw new RangeError(`${name}: xs must have at least ${minPoints} point(s)`);
  for (let i = 1; i < xs.length; i++) {
    if (!(xs[i] > xs[i - 1])) throw new RangeError(`${name}: xs must be strictly ascending (violated at index ${i})`);
  }
}

/** Index i of the interval [xs[i], xs[i+1]] containing x, clamped to 0..n−2 (binary search). */
function locate(xs, x) {
  const n = xs.length;
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid;
    else hi = mid;
  }
  return lo;
}

function checkExtrap(mode, name, allowed) {
  if (!allowed.includes(mode)) {
    throw new RangeError(`${name}: extrapolate must be one of ${allowed.map((s) => `'${s}'`).join(', ')}`);
  }
}

function linearEval(xs, ys, x, extrapolate) {
  const n = xs.length;
  if (Number.isNaN(x)) return NaN;
  if (n === 1) {
    if (extrapolate === 'error' && x !== xs[0]) throw new RangeError(`interp1: x = ${x} is outside the table`);
    return ys[0];
  }
  if (x < xs[0] || x > xs[n - 1]) {
    if (extrapolate === 'error') {
      throw new RangeError(`interp1: x = ${x} is outside the table range [${xs[0]}, ${xs[n - 1]}]`);
    }
    if (extrapolate === 'clamp') return x < xs[0] ? ys[0] : ys[n - 1];
  }
  const i = locate(xs, x);
  const t = (x - xs[i]) / (xs[i + 1] - xs[i]);
  return ys[i] + t * (ys[i + 1] - ys[i]);
}

/**
 * Piecewise-linear interpolation in a table.
 * Validity: xs strictly ascending (checked); outside the table the behaviour follows `extrapolate`.
 * @param {ArrayLike<number>} xs abscissae, ascending [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {number} x query point [x units]
 * @param {{ extrapolate?: 'clamp' | 'linear' | 'error' }} [opts] 'clamp' (default) holds the end value, 'linear'
 *   extends the end segment, 'error' throws RangeError
 * @returns {number} interpolated value [y units]
 */
export function interp1(xs, ys, x, { extrapolate = 'clamp' } = {}) {
  checkTable(xs, ys, 'interp1');
  checkExtrap(extrapolate, 'interp1', ['clamp', 'linear', 'error']);
  return linearEval(xs, ys, x, extrapolate);
}

/**
 * Build a reusable piecewise-linear interpolant (table validated once).
 * @param {ArrayLike<number>} xs abscissae, strictly ascending [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {{ extrapolate?: 'clamp' | 'linear' | 'error' }} [opts] as for interp1 (default 'clamp')
 * @returns {(x: number) => number} interpolant returning [y units]
 */
export function makeLinear(xs, ys, { extrapolate = 'clamp' } = {}) {
  checkTable(xs, ys, 'makeLinear');
  checkExtrap(extrapolate, 'makeLinear', ['clamp', 'linear', 'error']);
  const X = Float64Array.from(xs);
  const Y = Float64Array.from(ys);
  return (x) => linearEval(X, Y, x, extrapolate);
}

/** Evaluate a piecewise cubic Hermite (values Y, slopes D) at x; order 0, 1 or 2 selects value/derivatives. */
function hermiteEval(X, Y, D, x, order, extrapolate, name) {
  const n = X.length;
  if (Number.isNaN(x)) return NaN;
  if (n === 1) return order === 0 ? Y[0] : 0;
  if (x < X[0] || x > X[n - 1]) {
    if (extrapolate === 'error') {
      throw new RangeError(`${name}: x = ${x} is outside the table range [${X[0]}, ${X[n - 1]}]`);
    }
    const e = x < X[0] ? 0 : n - 1;
    if (extrapolate === 'clamp') return order === 0 ? Y[e] : 0;
    if (extrapolate === 'linear') {
      if (order === 0) return Y[e] + D[e] * (x - X[e]);
      return order === 1 ? D[e] : 0;
    }
    // 'extend': fall through and evaluate the end cubic outside its interval
  }
  const i = locate(X, x);
  const h = X[i + 1] - X[i];
  const t = x - X[i];
  const delta = (Y[i + 1] - Y[i]) / h;
  const c2 = (3 * delta - 2 * D[i] - D[i + 1]) / h;
  const c3 = (D[i] + D[i + 1] - 2 * delta) / (h * h);
  if (order === 0) return Y[i] + t * (D[i] + t * (c2 + t * c3));
  if (order === 1) return D[i] + t * (2 * c2 + 3 * t * c3);
  return 2 * c2 + 6 * t * c3;
}

/**
 * Shape-preserving piecewise cubic Hermite interpolant (PCHIP): monotone data give a monotone curve, no
 * overshoot. Interior slopes are the weighted harmonic means of Fritsch & Butland (1984); end slopes use the
 * three-point non-centred formula with the monotonicity limiter of Fritsch & Carlson (1980). C¹ continuous.
 * Validity: xs strictly ascending, at least one point (two points give a straight line).
 * @param {ArrayLike<number>} xs abscissae [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {{ extrapolate?: 'clamp' | 'linear' | 'extend' | 'error' }} [opts] behaviour outside the table: 'clamp'
 *   (default) holds the end value, 'linear' continues with the end slope, 'extend' evaluates the end cubic,
 *   'error' throws RangeError
 * @returns {((x: number) => number) & { d1: (x: number) => number }} interpolant [y units]; fn.d1 is its first
 *   derivative [y units / x units]
 */
export function makePchip(xs, ys, { extrapolate = 'clamp' } = {}) {
  checkTable(xs, ys, 'makePchip');
  checkExtrap(extrapolate, 'makePchip', ['clamp', 'linear', 'extend', 'error']);
  const n = xs.length;
  const X = Float64Array.from(xs);
  const Y = Float64Array.from(ys);
  const D = new Float64Array(n);
  if (n === 2) {
    D[0] = (Y[1] - Y[0]) / (X[1] - X[0]);
    D[1] = D[0];
  } else if (n > 2) {
    const h = new Float64Array(n - 1);
    const del = new Float64Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
      h[i] = X[i + 1] - X[i];
      del[i] = (Y[i + 1] - Y[i]) / h[i];
    }
    for (let i = 1; i < n - 1; i++) {
      if (del[i - 1] * del[i] > 0) {
        const w1 = 2 * h[i] + h[i - 1];
        const w2 = h[i] + 2 * h[i - 1];
        D[i] = (w1 + w2) / (w1 / del[i - 1] + w2 / del[i]);
      } else {
        D[i] = 0;
      }
    }
    const endSlope = (h0, h1, d0, d1) => {
      let d = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
      if (d * d0 <= 0) d = 0;
      else if (d0 * d1 <= 0 && Math.abs(d) > 3 * Math.abs(d0)) d = 3 * d0;
      return d;
    };
    D[0] = endSlope(h[0], h[1], del[0], del[1]);
    D[n - 1] = endSlope(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
  }
  const fn = (x) => hermiteEval(X, Y, D, x, 0, extrapolate, 'makePchip');
  fn.d1 = (x) => hermiteEval(X, Y, D, x, 1, extrapolate, 'makePchip');
  return fn;
}

/**
 * C² cubic spline interpolant.
 * End conditions (de Boor ch. IV): 'natural' (zero second derivative at both ends), 'not-a-knot' (third
 * derivative continuous across the first and last interior knots; reproduces cubics exactly), or 'clamped'
 * (end slopes d0 and dn prescribed). The tridiagonal system for the knot second derivatives is solved by the
 * Thomas algorithm (diagonally dominant).
 * Validity: xs strictly ascending, at least two points (two points give a straight line; 'not-a-knot' with
 * three points gives the interpolating parabola).
 * @param {ArrayLike<number>} xs abscissae [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {{ bc?: 'natural' | 'not-a-knot' | 'clamped', d0?: number, dn?: number,
 *   extrapolate?: 'linear' | 'clamp' | 'extend' | 'error' }} [opts] bc: end condition (default 'natural');
 *   d0, dn: end slopes for 'clamped' [y units / x units] (default 0); extrapolate: outside the table 'linear'
 *   (default) continues with the end slope and zero curvature, 'clamp' holds the end value, 'extend' evaluates
 *   the end cubic, 'error' throws RangeError
 * @returns {((x: number) => number) & { d1: (x: number) => number, d2: (x: number) => number }} spline value
 *   [y units]; fn.d1 first derivative [y/x]; fn.d2 second derivative [y/x²]
 */
export function makeSpline(xs, ys, { bc = 'natural', d0 = 0, dn = 0, extrapolate = 'linear' } = {}) {
  checkTable(xs, ys, 'makeSpline', 2);
  checkExtrap(extrapolate, 'makeSpline', ['linear', 'clamp', 'extend', 'error']);
  if (bc !== 'natural' && bc !== 'not-a-knot' && bc !== 'clamped') {
    throw new RangeError("makeSpline: bc must be 'natural', 'not-a-knot' or 'clamped'");
  }
  const n = xs.length;
  const X = Float64Array.from(xs);
  const Y = Float64Array.from(ys);
  const M = new Float64Array(n); // second derivatives at the knots
  const h = new Float64Array(n - 1);
  const del = new Float64Array(n - 1);
  for (let i = 0; i < n - 1; i++) {
    h[i] = X[i + 1] - X[i];
    del[i] = (Y[i + 1] - Y[i]) / h[i];
  }
  if (n === 2) {
    if (bc === 'clamped') {
      // single cubic Hermite segment: 2h M0 + h M1 = 6(δ − d0); h M0 + 2h M1 = 6(dn − δ)
      const r0 = (6 * (del[0] - d0)) / h[0];
      const r1 = (6 * (dn - del[0])) / h[0];
      M[0] = (2 * r0 - r1) / 3;
      M[1] = (2 * r1 - r0) / 3;
    }
  } else if (bc === 'not-a-knot' && n === 3) {
    const m = (2 * (del[1] - del[0])) / (h[0] + h[1]);
    M[0] = m;
    M[1] = m;
    M[2] = m;
  } else {
    const lo = new Float64Array(n);
    const di = new Float64Array(n);
    const up = new Float64Array(n);
    const rhs = new Float64Array(n);
    for (let i = 1; i < n - 1; i++) {
      lo[i] = h[i - 1];
      di[i] = 2 * (h[i - 1] + h[i]);
      up[i] = h[i];
      rhs[i] = 6 * (del[i] - del[i - 1]);
    }
    let first = 0;
    let last = n - 1;
    if (bc === 'natural') {
      di[0] = 1;
      di[n - 1] = 1;
    } else if (bc === 'clamped') {
      di[0] = 2 * h[0];
      up[0] = h[0];
      rhs[0] = 6 * (del[0] - d0);
      lo[n - 1] = h[n - 2];
      di[n - 1] = 2 * h[n - 2];
      rhs[n - 1] = 6 * (dn - del[n - 2]);
    } else {
      // not-a-knot: eliminate M0 = ((h0+h1) M1 − h0 M2)/h1 and M_{n−1} likewise; the reduced system stays
      // tridiagonal and diagonally dominant.
      first = 1;
      last = n - 2;
      const h0 = h[0];
      const h1 = h[1];
      di[1] = ((h0 + h1) * (h0 + 2 * h1)) / h1;
      up[1] = (h1 * h1 - h0 * h0) / h1;
      const ha = h[n - 3];
      const hb = h[n - 2];
      di[n - 2] = ((ha + hb) * (hb + 2 * ha)) / ha;
      lo[n - 2] = (ha * ha - hb * hb) / ha;
    }
    // Thomas algorithm on rows first..last
    for (let i = first + 1; i <= last; i++) {
      const w = lo[i] / di[i - 1];
      di[i] -= w * up[i - 1];
      rhs[i] -= w * rhs[i - 1];
    }
    M[last] = rhs[last] / di[last];
    for (let i = last - 1; i >= first; i--) M[i] = (rhs[i] - up[i] * M[i + 1]) / di[i];
    if (bc === 'not-a-knot') {
      M[0] = ((h[0] + h[1]) * M[1] - h[0] * M[2]) / h[1];
      M[n - 1] = ((h[n - 3] + h[n - 2]) * M[n - 2] - h[n - 2] * M[n - 3]) / h[n - 3];
    }
  }
  const slopeAt = (e) => {
    // end slopes of the spline itself
    if (e === 0) return del[0] - (h[0] * (2 * M[0] + M[1])) / 6;
    return del[n - 2] + (h[n - 2] * (M[n - 2] + 2 * M[n - 1])) / 6;
  };
  const s0 = slopeAt(0);
  const sn = slopeAt(1);
  const evalOrder = (x, order) => {
    if (Number.isNaN(x)) return NaN;
    if (x < X[0] || x > X[n - 1]) {
      if (extrapolate === 'error') {
        throw new RangeError(`makeSpline: x = ${x} is outside the table range [${X[0]}, ${X[n - 1]}]`);
      }
      const left = x < X[0];
      if (extrapolate === 'clamp') return order === 0 ? (left ? Y[0] : Y[n - 1]) : 0;
      if (extrapolate === 'linear') {
        const s = left ? s0 : sn;
        if (order === 0) return left ? Y[0] + s * (x - X[0]) : Y[n - 1] + s * (x - X[n - 1]);
        return order === 1 ? s : 0;
      }
    }
    const i = locate(X, x);
    const hi = h[i];
    const t = x - X[i];
    const b = del[i] - (hi * (2 * M[i] + M[i + 1])) / 6;
    const c = 0.5 * M[i];
    const d = (M[i + 1] - M[i]) / (6 * hi);
    if (order === 0) return Y[i] + t * (b + t * (c + t * d));
    if (order === 1) return b + t * (2 * c + 3 * t * d);
    return 2 * c + 6 * t * d;
  };
  const fn = (x) => evalOrder(x, 0);
  fn.d1 = (x) => evalOrder(x, 1);
  fn.d2 = (x) => evalOrder(x, 2);
  return fn;
}

/**
 * Bilinear interpolation on a rectilinear grid, clamped to the grid edges.
 * Validity: xs and ys strictly ascending; z has ys.length rows of xs.length values (z[iy][ix]). A grid with a
 * single point along an axis is treated as constant along that axis.
 * @param {ArrayLike<number>} xs grid abscissae along x [x units]
 * @param {ArrayLike<number>} ys grid abscissae along y [y units]
 * @param {ArrayLike<number>[]} z grid values, z[iy][ix] [z units]
 * @param {number} x query x [x units]
 * @param {number} y query y [y units]
 * @returns {number} interpolated value [z units]
 */
export function interp2(xs, ys, z, x, y) {
  if (!isVec(xs) || !isVec(ys) || xs.length < 1 || ys.length < 1) throw new RangeError('interp2: xs and ys must be non-empty arrays');
  if (!Array.isArray(z) || z.length !== ys.length) throw new RangeError('interp2: z must have ys.length rows');
  const nx = xs.length;
  const ny = ys.length;
  if (!isVec(z[0]) || z[0].length !== nx) throw new RangeError('interp2: each row of z must have xs.length values');
  if (Number.isNaN(x) || Number.isNaN(y)) return NaN;
  let ix = 0;
  let tx = 0;
  if (nx > 1) {
    const xc = Math.min(Math.max(x, xs[0]), xs[nx - 1]);
    ix = locate(xs, xc);
    const hx = xs[ix + 1] - xs[ix];
    if (!(hx > 0)) throw new RangeError('interp2: xs must be strictly ascending');
    tx = (xc - xs[ix]) / hx;
  }
  let iy = 0;
  let ty = 0;
  if (ny > 1) {
    const yc = Math.min(Math.max(y, ys[0]), ys[ny - 1]);
    iy = locate(ys, yc);
    const hy = ys[iy + 1] - ys[iy];
    if (!(hy > 0)) throw new RangeError('interp2: ys must be strictly ascending');
    ty = (yc - ys[iy]) / hy;
  }
  const ix1 = nx > 1 ? ix + 1 : ix;
  const iy1 = ny > 1 ? iy + 1 : iy;
  const z00 = z[iy][ix];
  const z01 = z[iy][ix1];
  const z10 = z[iy1][ix];
  const z11 = z[iy1][ix1];
  return (1 - ty) * ((1 - tx) * z00 + tx * z01) + ty * ((1 - tx) * z10 + tx * z11);
}

/**
 * n equally spaced points from a to b inclusive (end point exact).
 * @param {number} a first value
 * @param {number} b last value
 * @param {number} n number of points (>= 1; n = 1 returns [a])
 * @returns {Float64Array} grid
 */
export function linspace(a, b, n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('linspace: n must be an integer >= 1');
  const out = new Float64Array(n);
  if (n === 1) {
    out[0] = a;
    return out;
  }
  const step = (b - a) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = a + i * step;
  out[n - 1] = b;
  return out;
}

/**
 * n logarithmically spaced points from 10^a to 10^b inclusive.
 * @param {number} a exponent of the first value
 * @param {number} b exponent of the last value
 * @param {number} n number of points (>= 1)
 * @returns {Float64Array} grid
 */
export function logspace(a, b, n) {
  const e = linspace(a, b, n);
  for (let i = 0; i < n; i++) e[i] = Math.pow(10, e[i]);
  return e;
}

/**
 * Resample a tabulated curve onto n equally spaced abscissae spanning [xs[0], xs[last]] by linear interpolation.
 * Validity: xs strictly ascending, at least two points; n >= 2.
 * @param {ArrayLike<number>} xs abscissae [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {number} n number of output points
 * @returns {{ x: Float64Array, y: Float64Array }} uniform abscissae and interpolated ordinates
 */
export function resample(xs, ys, n) {
  checkTable(xs, ys, 'resample', 2);
  if (!Number.isInteger(n) || n < 2) throw new RangeError('resample: n must be an integer >= 2');
  const x = linspace(xs[0], xs[xs.length - 1], n);
  const y = new Float64Array(n);
  let j = 0;
  const last = xs.length - 2;
  for (let i = 0; i < n; i++) {
    const xi = x[i];
    while (j < last && xs[j + 1] < xi) j++;
    const t = (xi - xs[j]) / (xs[j + 1] - xs[j]);
    y[i] = ys[j] + t * (ys[j + 1] - ys[j]);
  }
  y[0] = ys[0];
  y[n - 1] = ys[ys.length - 1];
  return { x, y };
}

/**
 * Reduce a long series to at most maxPoints samples for plotting while keeping its envelope: the first and last
 * samples are kept and the interior is split into equal-count buckets, each contributing its minimum and maximum
 * sample in their original order. Series already short enough are returned as copies.
 * Validity: xs and ys of equal length (xs need not be uniform); maxPoints >= 4.
 * @param {ArrayLike<number>} xs abscissae [x units]
 * @param {ArrayLike<number>} ys ordinates [y units]
 * @param {number} maxPoints maximum number of output samples
 * @returns {{ x: Float64Array, y: Float64Array }} decimated series (length <= maxPoints)
 */
export function decimate(xs, ys, maxPoints) {
  if (!isVec(xs) || !isVec(ys) || xs.length !== ys.length) throw new RangeError('decimate: xs and ys must be arrays of equal length');
  if (!Number.isInteger(maxPoints) || maxPoints < 4) throw new RangeError('decimate: maxPoints must be an integer >= 4');
  const n = xs.length;
  if (n <= maxPoints) return { x: Float64Array.from(xs), y: Float64Array.from(ys) };
  const nb = (maxPoints - 2) >> 1;
  const idx = new Int32Array(2 * nb + 2);
  let cnt = 0;
  idx[cnt++] = 0;
  const inner = n - 2;
  for (let b = 0; b < nb; b++) {
    const s = 1 + Math.floor((b * inner) / nb);
    const e = 1 + Math.floor(((b + 1) * inner) / nb);
    if (e <= s) continue;
    let imin = s;
    let imax = s;
    for (let i = s + 1; i < e; i++) {
      if (ys[i] < ys[imin]) imin = i;
      if (ys[i] > ys[imax]) imax = i;
    }
    if (imin === imax) {
      idx[cnt++] = imin;
    } else if (imin < imax) {
      idx[cnt++] = imin;
      idx[cnt++] = imax;
    } else {
      idx[cnt++] = imax;
      idx[cnt++] = imin;
    }
  }
  idx[cnt++] = n - 1;
  const x = new Float64Array(cnt);
  const y = new Float64Array(cnt);
  for (let k = 0; k < cnt; k++) {
    x[k] = xs[idx[k]];
    y[k] = ys[idx[k]];
  }
  return { x, y };
}
