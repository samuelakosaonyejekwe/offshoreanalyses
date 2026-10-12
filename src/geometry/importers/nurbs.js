// Non-uniform rational B-spline evaluation for the CAD text readers (DXF SPLINE, IGES 126 / 128, STEP
// B_SPLINE_CURVE_WITH_KNOTS). Algorithms A2.1 (FindSpan), A2.2 (BasisFuns), A4.1 (CurvePoint) and A4.3
// (SurfacePoint) of Piegl & Tiller, "The NURBS Book", 2nd ed. (1997).
//
// Conventions: degree p, n + 1 control points, knot vector U of length n + p + 2 (non-decreasing).
// Control points are flat x, y, z; weights (optional) are one per control point.

/** Knot span index i such that U[i] ≤ u < U[i+1] (the last non-empty span for u = U[n+1]). */
export function findSpan(n, p, u, U) {
  if (u >= U[n + 1]) { let i = n; while (i > p && U[i] === U[i + 1]) i--; return i; }
  if (u <= U[p]) { let i = p; while (i < n && U[i + 1] === U[i]) i++; return i; }
  let lo = p, hi = n + 1, mid = (lo + hi) >> 1;
  while (u < U[mid] || u >= U[mid + 1]) { if (u < U[mid]) hi = mid; else lo = mid; mid = (lo + hi) >> 1; }
  return mid;
}

/** The p + 1 non-zero basis functions at u, written into N. */
export function basisFuns(span, u, p, U, N, left, right) {
  N[0] = 1;
  for (let j = 1; j <= p; j++) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;
    let saved = 0;
    for (let r = 0; r < j; r++) {
      const d = right[r + 1] + left[j - r];
      const t = d !== 0 ? N[r] / d : 0;
      N[r] = saved + right[r + 1] * t;
      saved = left[j - r] * t;
    }
    N[j] = saved;
  }
}

/** Validate / repair the inputs of a curve; returns null when the data cannot define a curve. */
function curveSetup(p, U, P, W) {
  const n = P.length / 3 - 1;
  if (!(p >= 1) || n < p || U.length < n + p + 2) return null;
  return { n, lo: U[p], hi: U[n + 1], rational: !!(W && W.length === n + 1) };
}

/**
 * Sample a (rational) B-spline curve at `count + 1` parameters, evenly spaced within each knot span so that
 * knots (where the curve may have corners) are always hit.
 * @param {number} p degree @param {ArrayLike<number>} U knots @param {ArrayLike<number>} P control points (flat xyz)
 * @param {ArrayLike<number>|null} W weights @param {number} [perSpan=8] segments per non-empty knot span
 * @param {number} [u0] @param {number} [u1] parameter range (defaults to the full curve)
 * @returns {Float64Array|null} flat xyz, or null when the definition is invalid
 */
export function sampleCurve(p, U, P, W, perSpan = 8, u0, u1) {
  const s = curveSetup(p, U, P, W);
  if (!s) return null;
  const lo = u0 !== undefined && u0 === u0 ? Math.max(s.lo, Math.min(u0, s.hi)) : s.lo;
  const hi = u1 !== undefined && u1 === u1 ? Math.max(s.lo, Math.min(u1, s.hi)) : s.hi;
  const params = [lo];
  for (let i = p; i <= s.n; i++) {
    const a = Math.max(U[i], lo), b = Math.min(U[i + 1], hi);
    if (!(b > a)) continue;
    for (let k = 1; k <= perSpan; k++) params.push(a + ((b - a) * k) / perSpan);
  }
  if (params.length < 2) return null;
  const out = new Float64Array(params.length * 3);
  const N = new Float64Array(p + 1), L = new Float64Array(p + 1), R = new Float64Array(p + 1);
  params.forEach((u, k) => {
    const span = findSpan(s.n, p, u, U);
    basisFuns(span, u, p, U, N, L, R);
    let x = 0, y = 0, z = 0, w = 0;
    for (let j = 0; j <= p; j++) {
      const i = span - p + j, c = s.rational ? N[j] * W[i] : N[j];
      x += c * P[3 * i]; y += c * P[3 * i + 1]; z += c * P[3 * i + 2]; w += c;
    }
    if (!s.rational || w === 0) w = 1;
    out[3 * k] = x / w; out[3 * k + 1] = y / w; out[3 * k + 2] = z / w;
  });
  return out;
}

/**
 * Tessellate a (rational) B-spline surface into triangles.
 * Control net: (nu + 1) × (nv + 1) points with index iu + (nu + 1)·iv (u varies fastest, as in IGES 128).
 * @returns {{ coords: Float64Array, indices: Uint32Array, nu: number, nv: number }|null}
 */
export function tessellateSurface(pu, pv, U, V, P, W, nu, nv, perSpan = 6, range) {
  if (!(pu >= 1 && pv >= 1) || nu < pu || nv < pv || U.length < nu + pu + 2 || V.length < nv + pv + 2 || P.length < 3 * (nu + 1) * (nv + 1)) return null;
  const grid = (p, K, n, lo, hi) => {
    const a0 = lo !== undefined && lo === lo ? Math.max(K[p], lo) : K[p], b0 = hi !== undefined && hi === hi ? Math.min(K[n + 1], hi) : K[n + 1];
    const out = [a0];
    for (let i = p; i <= n; i++) { const a = Math.max(K[i], a0), b = Math.min(K[i + 1], b0); if (!(b > a)) continue; for (let k = 1; k <= perSpan; k++) out.push(a + ((b - a) * k) / perSpan); }
    return out;
  };
  const us = grid(pu, U, nu, range && range[0], range && range[1]), vs = grid(pv, V, nv, range && range[2], range && range[3]);
  if (us.length < 2 || vs.length < 2) return null;
  const rational = !!(W && W.length === (nu + 1) * (nv + 1));
  const Nu = new Float64Array(pu + 1), Nv = new Float64Array(pv + 1), L = new Float64Array(Math.max(pu, pv) + 1), R = new Float64Array(Math.max(pu, pv) + 1);
  const coords = new Float64Array(us.length * vs.length * 3);
  const spansU = us.map((u) => findSpan(nu, pu, u, U));
  const basisU = us.map((u, k) => { basisFuns(spansU[k], u, pu, U, Nu, L, R); return Float64Array.from(Nu); });
  vs.forEach((v, jv) => {
    const sv = findSpan(nv, pv, v, V);
    basisFuns(sv, v, pv, V, Nv, L, R);
    us.forEach((u, ju) => {
      const su = spansU[ju], bu = basisU[ju];
      let x = 0, y = 0, z = 0, w = 0;
      for (let b = 0; b <= pv; b++) for (let a = 0; a <= pu; a++) {
        const i = su - pu + a + (nu + 1) * (sv - pv + b);
        const c = bu[a] * Nv[b] * (rational ? W[i] : 1);
        x += c * P[3 * i]; y += c * P[3 * i + 1]; z += c * P[3 * i + 2]; w += c;
      }
      if (!rational || w === 0) w = 1;
      const o = 3 * (ju + us.length * jv);
      coords[o] = x / w; coords[o + 1] = y / w; coords[o + 2] = z / w;
    });
  });
  const m = us.length, q = vs.length;
  const indices = new Uint32Array((m - 1) * (q - 1) * 6);
  let k = 0;
  for (let j = 0; j + 1 < q; j++) for (let i = 0; i + 1 < m; i++) {
    const a = i + m * j, b = a + 1, c = a + m + 1, d = a + m;
    indices[k++] = a; indices[k++] = b; indices[k++] = c; indices[k++] = a; indices[k++] = c; indices[k++] = d;
  }
  return { coords, indices, nu: m, nv: q };
}

/** Points of a circular arc: centre c, radius r, from angle a0 sweeping `sweep` radians, in the plane (ex, ey). */
export function arcPoints(c, ex, ey, r, a0, sweep, stepDeg = 5) {
  const n = Math.max(2, Math.ceil(Math.abs(sweep) / ((stepDeg * Math.PI) / 180)));
  const out = new Float64Array((n + 1) * 3);
  for (let k = 0; k <= n; k++) {
    const t = a0 + (sweep * k) / n, cs = r * Math.cos(t), sn = r * Math.sin(t);
    out[3 * k] = c[0] + cs * ex[0] + sn * ey[0]; out[3 * k + 1] = c[1] + cs * ex[1] + sn * ey[1]; out[3 * k + 2] = c[2] + cs * ex[2] + sn * ey[2];
  }
  return out;
}
