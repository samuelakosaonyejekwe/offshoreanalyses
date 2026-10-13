// Tessellation of boundary-representation faces, shared by the B-rep text readers (acis.js, sab.js through
// acis.js, parasolid.js).
//
// A reader supplies, per face, a surface (one of the analytic surfaces below, with its 3-D ↔ parameter maps) and
// the boundary loops as closed 3-D polylines built from the polylines of the shared edges. The face is
// triangulated in the parameter plane of its surface between those loops: the loop polygons are ear-clipped and
// interior edges that span more than one angular step of the surface are bisected until none is left. Boundary
// segments are never subdivided, so two faces that share an edge use exactly the same points along it and the
// body comes out watertight when the file's topology is.
//
// Parameters are chosen so that ∂u × ∂v points away from the axis / centre (plane: along its normal);
// `inward` records that the surface normal of the file is the opposite.

import { V3, Welder } from './common.js';
import { findSpan, basisFuns } from './nurbs.js';

const TAU = 2 * Math.PI;
const nearest = (u, ref) => u + TAU * Math.round((ref - u) / TAU);

function anyPerp(n) {
  const s = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  return V3.norm(V3.sub(s, V3.scale(n, V3.dot(s, n))));
}
function perpTo(x, n) {
  if (!x || !(V3.len(x) > 0)) return anyPerp(n);
  const p = V3.sub(x, V3.scale(n, V3.dot(x, n)));
  return V3.len(p) > 1e-9 * V3.len(x) ? V3.norm(p) : anyPerp(n);
}

export { anyPerp, perpTo };

/** Angular step [rad] for a curve of radius r so that the sag stays below `defl`, at most `angStep`. */
export function stepper(defl, angStep) {
  const lowest = (0.25 * Math.PI) / 180;
  return (r) => {
    let s = angStep;
    if (r > defl) s = Math.min(s, 2 * Math.acos(1 - defl / r));
    return Math.max(s, lowest);
  };
}

/** Plane through O with unit normal N; `xdir` (optional) fixes the direction of u. */
export function planeSurface(O, N, xdir) {
  const e1 = perpTo(xdir || null, N), e2 = V3.cross(N, e1);
  return {
    type: 'plane', inward: false, perU: false, perV: false,
    toParam: (p) => { const d = V3.sub(p, O); return [V3.dot(d, e1), V3.dot(d, e2)]; },
    fromParam: (u, v) => [O[0] + u * e1[0] + v * e2[0], O[1] + u * e1[1] + v * e2[1], O[2] + u * e1[2] + v * e2[2]],
    scales: () => ({ su: 1, sv: 1, stepU: Infinity, stepV: Infinity }),
  };
}

/**
 * Cylinder or cone about the axis (C, A): the section at axial distance v from C is an ellipse with semi-axes
 * R0 + v·tan along X and ratio·(R0 + v·tan) along A × X.
 * @param {boolean} inward surface normal of the file points towards the axis
 */
export function coneSurface(C, A, X0, R0, ratio, tan, inward, stepFor) {
  const X = perpTo(X0, A), Y = V3.cross(A, X);
  const eps = 1e-9 * R0;
  return {
    type: Math.abs(tan) < 1e-12 ? 'cylinder' : 'cone', inward, perU: true, perV: false, straightV: true,
    apex: Math.abs(tan) > 1e-12 ? -R0 / tan : null,
    toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y) / ratio; return [Math.hypot(x, y) > eps ? Math.atan2(y, x) : NaN, V3.dot(d, A)]; },
    fromParam: (u, v) => { const r = R0 + v * tan, cx = r * Math.cos(u), cy = r * ratio * Math.sin(u); return [C[0] + cx * X[0] + cy * Y[0] + v * A[0], C[1] + cx * X[1] + cy * Y[1] + v * A[1], C[2] + cx * X[2] + cy * Y[2] + v * A[2]]; },
    radiusAt: (v) => Math.abs(R0 + v * tan) * Math.max(1, ratio),
    // generators are straight: only the angle a triangle spans matters
    scales: (rmax) => { const st = stepFor(rmax); return { su: 1 / st, sv: 1 / (8 * rmax * st), stepU: st, stepV: Infinity, squeeze: 0.6 }; },
  };
}

/** Sphere of radius R about C; Z is the pole axis, X0 the direction of longitude 0. */
export function sphereSurface(C, R, Z, X0, inward, stepFor) {
  const X = perpTo(X0, Z), Y = V3.cross(Z, X);
  return {
    type: 'sphere', inward, perU: true, perV: false, poles: true, C, R,
    toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y), z = V3.dot(d, Z); const h = Math.hypot(x, y); return [h > 1e-9 * R ? Math.atan2(y, x) : NaN, Math.atan2(z, h)]; },
    fromParam: (u, v) => { const c = R * Math.cos(v), cx = c * Math.cos(u), cy = c * Math.sin(u), cz = R * Math.sin(v); return [C[0] + cx * X[0] + cy * Y[0] + cz * Z[0], C[1] + cx * X[1] + cy * Y[1] + cz * Z[1], C[2] + cx * X[2] + cy * Y[2] + cz * Z[2]]; },
    scales: () => { const st = stepFor(R); return { su: 1 / st, sv: 1 / st, stepU: st, stepV: st }; },
  };
}

/**
 * Torus about the axis (C, N) with major radius R and tube radius r. A ring torus (R > r) is periodic both ways.
 * Of a spindle torus (0 < R < r, the tube crosses the axis) the outer sheet is read: it is bounded by the two
 * points where it meets the axis, which are treated like the poles of a sphere; the inner ("lemon") sheet is not
 * handled and faces on it are refused.
 */
export function torusSurface(C, N, R, r, X0, inward, stepFor) {
  if (!(r > 0) || !(R > 0)) return null;
  const X = perpTo(X0, N), Y = V3.cross(N, X);
  const fromParam = (u, v) => { const c = R + r * Math.cos(v), cx = c * Math.cos(u), cy = c * Math.sin(u), cz = r * Math.sin(v); return [C[0] + cx * X[0] + cy * Y[0] + cz * N[0], C[1] + cx * X[1] + cy * Y[1] + cz * N[1], C[2] + cx * X[2] + cy * Y[2] + cz * N[2]]; };
  // around the axis the facets may not sag more than they do around the tube, or the chords of a slender
  // torus would cut through it
  const scales = () => { const b = stepFor(r), a = Math.min(stepFor(R + r), 2 * Math.acos(1 - (r * (1 - Math.cos(b / 2))) / (R + r))); return { su: 1 / a, sv: 1 / b, stepU: a, stepV: b }; };
  if (R > r) {
    return {
      type: 'torus', inward, perU: true, perV: true, C, R, r, fromParam, scales,
      toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y); return [Math.atan2(y, x), Math.atan2(V3.dot(d, N), Math.hypot(x, y) - R)]; },
    };
  }
  if (!(R < r)) return null;
  const vPole = Math.acos(-R / r);
  return {
    type: 'torus', inward, perU: true, perV: false, poleRows: [-vPole, vPole], poleMid: 0, C, R, r, fromParam, scales,
    toParam: (p) => {
      const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y), z = V3.dot(d, N), h = Math.hypot(x, y);
      // a point that fits the inner sheet better than the outer one cannot be placed
      if (Math.abs(Math.hypot(h + R, z) - r) < Math.abs(Math.hypot(h - R, z) - r) - 1e-9 * r) return [NaN, NaN];
      return [h > 1e-9 * r ? Math.atan2(y, x) : NaN, Math.atan2(z, h - R)];
    },
  };
}

/**
 * B-spline surface (rational when `W` is given) as a face surface. Both parameters are rescaled to 0 … 2π so
 * that a surface closed in a direction can be treated like the periodic analytic surfaces. A surface closed in v
 * only, or collapsed to a point along an edge u = const, is transposed (which reverses ∂u × ∂v: `inward`).
 * @param {{ pu, pv, U, V, nu, nv, P, W }} net degrees, full knot vectors, numbers of control points and the control
 *   net with index iu + nu·iv
 * @param {number} defl largest sag of a facet [model units] @param {number} angStep largest turn per facet [rad]
 * @returns {object|null} null when the net is inconsistent or collapsed in a way that is not handled
 */
export function splineSurface(net, defl, angStep, swapped = false) {
  const { pu, pv, U, V, nu, nv, P, W } = net;
  if (!(pu >= 1 && pv >= 1 && nu > pu && nv > pv) || U.length < nu + pu + 1 || V.length < nv + pv + 1 || P.length < 3 * nu * nv) return null;
  const ua = U[pu], ub = U[nu], va = V[pv], vb = V[nv];
  if (!(ub > ua) || !(vb > va)) return null;
  const Nu = new Float64Array(pu + 1), Nv = new Float64Array(pv + 1), L = new Float64Array(Math.max(pu, pv) + 1), R = new Float64Array(Math.max(pu, pv) + 1);
  const at = (u, v) => {
    const su = findSpan(nu - 1, pu, u, U);
    basisFuns(su, u, pu, U, Nu, L, R);
    const sv = findSpan(nv - 1, pv, v, V);
    basisFuns(sv, v, pv, V, Nv, L, R);
    let x = 0, y = 0, z = 0, w = 0;
    for (let b = 0; b <= pv; b++) for (let a = 0; a <= pu; a++) {
      const i = su - pu + a + nu * (sv - pv + b);
      const c = Nu[a] * Nv[b] * (W ? W[i] : 1);
      x += c * P[3 * i]; y += c * P[3 * i + 1]; z += c * P[3 * i + 2]; w += c;
    }
    if (!W || w === 0) w = 1;
    return [x / w, y / w, z / w];
  };
  const spans = (K, p, n) => { let c = 0; for (let i = p; i < n; i++) if (K[i + 1] > K[i]) c++; return c; };
  const gu = Math.max(8, Math.min(48, 4 * spans(U, pu, nu))), gv = Math.max(8, Math.min(48, 4 * spans(V, pv, nv)));
  const G = [];
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let j = 0; j <= gv; j++) for (let i = 0; i <= gu; i++) {
    const p = at(ua + ((ub - ua) * i) / gu, va + ((vb - va) * j) / gv);
    if (!(p[0] === p[0] && p[1] === p[1] && p[2] === p[2])) return null;
    G.push(p);
    for (let k = 0; k < 3; k++) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
  }
  const size = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  if (!(size > 0)) return null;
  const g = (i, j) => G[i + (gu + 1) * j];
  const d3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  const same = 1e-6 * size;
  let closedU = true, closedV = true, rowLo = true, rowHi = true, colLo = true, colHi = true;
  for (let j = 0; j <= gv; j++) { if (d3(g(0, j), g(gu, j)) > same) closedU = false; if (d3(g(0, j), g(0, 0)) > same) colLo = false; if (d3(g(gu, j), g(gu, 0)) > same) colHi = false; }
  for (let i = 0; i <= gu; i++) { if (d3(g(i, 0), g(i, gv)) > same) closedV = false; if (d3(g(i, 0), g(0, 0)) > same) rowLo = false; if (d3(g(i, gv), g(0, gv)) > same) rowHi = false; }
  if ((colLo || colHi) || (closedV && !closedU)) {
    if (swapped || ((colLo || colHi) && (rowLo || rowHi))) return null;
    const Q = new Float64Array(P.length), Wq = W ? new Float64Array(W.length) : null;
    for (let iu = 0; iu < nu; iu++) for (let iv = 0; iv < nv; iv++) { const from = iu + nu * iv, to = iv + nv * iu; Q[3 * to] = P[3 * from]; Q[3 * to + 1] = P[3 * from + 1]; Q[3 * to + 2] = P[3 * from + 2]; if (Wq) Wq[to] = W[from]; }
    return splineSurface({ pu: pv, pv: pu, U: V, V: U, nu: nv, nv: nu, P: Q, W: Wq }, defl, angStep, true);
  }
  // facets per direction from the sag and the turn of the sampled iso-lines
  const need = (along) => {
    const n = along ? gu : gv, m = along ? gv : gu;
    let sag = 0, turn = 0;
    for (let j = 0; j <= m; j++) for (let i = 1; i < n; i++) {
      const a = along ? g(i - 1, j) : g(j, i - 1), b = along ? g(i, j) : g(j, i), c = along ? g(i + 1, j) : g(j, i + 1);
      sag = Math.max(sag, d3(b, [(a[0] + c[0]) / 2, (a[1] + c[1]) / 2, (a[2] + c[2]) / 2]));
      const e1 = V3.sub(b, a), e2 = V3.sub(c, b), l = V3.len(e1) * V3.len(e2);
      if (l > 1e-12 * size * size) turn = Math.max(turn, Math.acos(Math.max(-1, Math.min(1, V3.dot(e1, e2) / l))));
    }
    return Math.max(1, Math.min(512, Math.ceil(n * Math.max(Math.sqrt(sag / (4 * defl)), turn / angStep) - 1e-9)));
  };
  const stepU = TAU / need(true), stepV = TAU / need(false);
  const uOf = (x) => { let f = x / TAU; if (closedU) f -= Math.floor(f); return ua + (ub - ua) * Math.max(0, Math.min(1, f)); };
  const vOf = (y) => { let f = y / TAU; if (closedV) f -= Math.floor(f); return va + (vb - va) * Math.max(0, Math.min(1, f)); };
  const fromParam = (x, y) => at(uOf(x), vOf(y));
  const hx = 1e-5 * TAU;
  const toParam = (p) => {
    let bi = 0, bd = Infinity;
    for (let k = 0; k < G.length; k++) { const q = G[k]; const d = (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2; if (d < bd) { bd = d; bi = k; } }
    let x = (TAU * (bi % (gu + 1))) / gu, y = (TAU * Math.floor(bi / (gu + 1))) / gv;
    const lim = (t, closed) => (closed ? t : Math.max(0, Math.min(TAU, t)));
    let S = fromParam(x, y), Su, Sv;
    const derivs = () => {
      const x0 = Math.max(0, x - hx), x1 = Math.min(TAU, x + hx), y0 = Math.max(0, y - hx), y1 = Math.min(TAU, y + hx);
      Su = V3.scale(V3.sub(fromParam(x1, y), fromParam(x0, y)), 1 / (x1 - x0));
      Sv = V3.scale(V3.sub(fromParam(x, y1), fromParam(x, y0)), 1 / (y1 - y0));
    };
    for (let it = 0; it < 24; it++) {
      derivs();
      const r = V3.sub(p, S), a = V3.dot(Su, Su), b = V3.dot(Su, Sv), c = V3.dot(Sv, Sv), e = V3.dot(Su, r), f = V3.dot(Sv, r);
      const det = a * c - b * b;
      let dx, dy;
      if (det > 1e-12 * a * c && det > 0) { dx = (e * c - f * b) / det; dy = (a * f - b * e) / det; } else if (c > a) { dx = 0; dy = c > 0 ? f / c : 0; } else { dy = 0; dx = a > 0 ? e / a : 0; }
      // damped: a step is taken only as far as it brings the point closer
      let k = 1, nx, ny, Sn, dn;
      const d0 = V3.len(r);
      for (;;) { nx = lim(x + k * dx, closedU); ny = lim(y + k * dy, closedV); Sn = fromParam(nx, ny); dn = d3(Sn, p); if (dn <= d0 || k < 1e-3) break; k /= 2; }
      if (dn > d0) break;
      const moved = Math.abs(nx - x) + Math.abs(ny - y);
      x = nx; y = ny; S = Sn;
      if (moved < 1e-10) break;
    }
    derivs();
    const dist = d3(S, p);
    if (closedU) x -= TAU * Math.floor(x / TAU);
    if (closedV) y -= TAU * Math.floor(y / TAU);
    // on a collapsed row the u parameter means nothing
    return [V3.len(Su) * TAU < 1e-5 * size ? NaN : x, y, dist];
  };
  return {
    type: 'spline', inward: swapped, perU: closedU, perV: closedV, size,
    // rescaled knots, where the surface may have creases
    knotsU: [...new Set(U.slice(pu, nu + 1))].map((k) => (TAU * (k - ua)) / (ub - ua)), knotsV: [...new Set(V.slice(pv, nv + 1))].map((k) => (TAU * (k - va)) / (vb - va)),
    poleRows: rowLo || rowHi ? [rowLo ? 0 : null, rowHi ? TAU : null] : null, poleMid: Math.PI,
    toParam, fromParam,
    scales: () => ({ su: 1 / stepU, sv: 1 / stepV, stepU, stepV }),
  };
}

// ---------------------------------------------------------------------------------------------------------
// Face domain in parameter space
// ---------------------------------------------------------------------------------------------------------

const areaOf = (L) => { let a = 0; for (let i = 0, j = L.length - 1; i < L.length; j = i++) a += L[j].x * L[i].y - L[i].x * L[j].y; return a / 2; };
const meanOf = (L, k) => { let s = 0; for (const q of L) s += q[k]; return s / L.length; };
const shifted = (L, dx, dy = 0) => (dx || dy ? L.map((q) => ({ x: q.x + dx, y: q.y + dy, p: q.p })) : L);

/**
 * Boundary loops of a face in the parameter plane of a periodic surface, unwrapped. A vertex at a pole or apex
 * (where u is undefined) becomes the stretch of the degenerate row between the u of its neighbours.
 * @returns {{ pts: Array<{x, y, p}>, wu: number, wv: number }|null}
 */
function unwrapLoop(surf, pts3, stepU, flip = false) {
  const raw = pts3.map((p) => { const [u, v] = surf.toParam(p); return { u, v, p, sing: u !== u }; });
  if (raw.some((q) => q.v !== q.v)) return OFF_SHEET;
  const m = raw.length;
  if (raw.every((q) => q.sing)) return null;
  const out = [];
  const firstSing = raw.findIndex((q, i) => q.sing && !raw[(i + 1) % m].sing);
  const fixV = (v, ref) => (surf.perV ? nearest(v, ref) : v);
  const fixU = (u, ref) => (surf.perU ? nearest(u, ref) : u);
  if (firstSing < 0) {
    let lu = raw[0].u, lv = raw[0].v;
    out.push({ x: lu, y: lv, p: raw[0].p });
    for (let i = 1; i < m; i++) { lu = fixU(raw[i].u, lu); lv = fixV(raw[i].v, lv); out.push({ x: lu, y: lv, p: raw[i].p }); }
    return { pts: out, wu: surf.perU ? Math.round((nearest(raw[0].u, lu) - raw[0].u) / TAU) : 0, wv: surf.perV ? Math.round((nearest(raw[0].v, lv) - raw[0].v) / TAU) : 0 };
  }
  // start just after a degenerate vertex so that the loop closes along the degenerate row without winding
  const seq = [];
  for (let k = 1; k <= m; k++) seq.push(raw[(firstSing + k) % m]);
  let lu = seq[0].u, lv = seq[0].v;
  out.push({ x: lu, y: lv, p: seq[0].p });
  for (let i = 1; i < m; i++) {
    const q = seq[i];
    if (!q.sing) { lu = fixU(q.u, lu); lv = fixV(q.v, lv); out.push({ x: lu, y: lv, p: q.p }); continue; }
    out.push({ x: lu, y: q.v, p: q.p });
    let nx = null;
    for (let k = i + 1; k < m && !nx; k++) if (!seq[k].sing) nx = seq[k];
    let un = nx ? fixU(nx.u, lu) : out[0].x;
    // a boundary that runs straight through a pole continues on the opposite meridian; the face lies to the
    // left of its boundary, which decides the way round the pole
    if (nx && (surf.poles || surf.poleRows) && Math.abs(Math.abs(un - lu) - Math.PI) < 1e-6) un = lu + (q.v > (surf.poles ? 0 : surf.poleMid ?? Math.PI) ? -1 : 1) * (flip ? -1 : 1) * Math.PI;
    // the degenerate row is cut into steps (all at the same 3-D point) so that no triangle fans across it
    if (un !== lu) { const k = Math.max(1, Math.ceil(Math.abs(un - lu) / (stepU || TAU) - 1e-9)); for (let j = 1; j <= k; j++) out.push({ x: lu + ((un - lu) * j) / k, y: q.v, p: q.p }); }
    lu = un;
  }
  return { pts: out, wu: 0, wv: 0 };
}

/**
 * Polygons (even–odd rule) bounding a face on a periodic surface.
 * @param {boolean} flip face normal opposes ∂u × ∂v, so loops run clockwise in the parameter plane
 * @returns {Array<Array<{x, y, p}>>|string} loops, or the reason the domain could not be built
 */
function periodicDomain(surf, loops, flip, sc) {
  const orient = flip ? -1 : 1;
  const mk = (x, y) => ({ x, y, p: surf.fromParam(x, y) });
  const copy = (q, dx, dy) => ({ x: q.x + dx, y: q.y + dy, p: q.p });
  const divs = (a, b) => Math.max(1, Math.ceil(Math.max(Math.abs(b.x - a.x) / sc.stepU, Number.isFinite(sc.stepV) ? Math.abs(b.y - a.y) / sc.stepV : 0) - 1e-9));
  // interior points of the straight parameter-space line a → b (shared by both sides of a seam)
  const between = (a, b) => { const k = divs(a, b), o = []; for (let i = 1; i < k; i++) o.push(mk(a.x + ((b.x - a.x) * i) / k, a.y + ((b.y - a.y) * i) / k)); return o; };
  const wind = loops.filter((l) => l.wu || l.wv), flat = loops.filter((l) => !l.wu && !l.wv);
  if (wind.some((l) => l.wu && l.wv)) return 'a boundary loop winds around both directions of the surface';

  // full parameter rectangle of a closed surface, cut where no hole lies
  const rectangle = (holes) => {
    const place = (list, k, lo) => list.map((h) => { const d = -TAU * Math.floor((meanOf(h, k) - lo) / TAU); return k === 'x' ? shifted(h, d, 0) : shifted(h, 0, d); });
    const clear = (hs, k, lo) => hs.every((h) => h.every((q) => q[k] > lo + 1e-6 && q[k] < lo + TAU - 1e-6));
    let u0 = null, hs = holes;
    for (let c = 0; c < 16 && u0 === null; c++) { const lo = -Math.PI + 0.0137 + (c * Math.PI) / 8; const t = place(holes, 'x', lo); if (clear(t, 'x', lo)) { u0 = lo; hs = t; } }
    if (u0 === null) return 'its holes leave no meridian free for the seam';
    if (surf.poles) {
      const s = mk(u0, -Math.PI / 2), nth = mk(u0, Math.PI / 2);
      const side = between(s, nth);
      const row = between(s, { x: u0 + TAU, y: s.y });
      return [[s, ...row.map((q) => ({ x: q.x, y: q.y, p: s.p })), copy(s, TAU, 0), ...side.map((q) => copy(q, TAU, 0)), copy(nth, TAU, 0), ...row.slice().reverse().map((q) => ({ x: q.x, y: nth.y, p: nth.p })), nth, ...side.slice().reverse()], ...hs];
    }
    let v0 = null;
    for (let c = 0; c < 16 && v0 === null; c++) { const lo = -Math.PI + 0.0137 + (c * Math.PI) / 8; const t = place(hs, 'y', lo); if (clear(t, 'y', lo)) { v0 = lo; hs = t; } }
    if (v0 === null) return 'its holes leave no parallel free for the seam';
    const o = mk(u0, v0), bu = between(o, { x: u0 + TAU, y: v0 }), bv = between(o, { x: u0, y: v0 + TAU });
    return [[o, ...bu, copy(o, TAU, 0), ...bv.map((q) => copy(q, TAU, 0)), copy(o, TAU, TAU), ...bu.slice().reverse().map((q) => copy(q, 0, TAU)), copy(o, 0, TAU), ...bv.slice().reverse()], ...hs];
  };

  if (!wind.length) {
    if (!flat.length) return surf.poles || surf.perV ? rectangle([]) : 'the face has no boundary on an unbounded surface';
    const withArea = flat.map((l) => ({ l, a: areaOf(l.pts) }));
    if (withArea.some((q) => q.a * orient > 0) || !(surf.poles || surf.perV)) {
      withArea.sort((p, q) => Math.abs(q.a) - Math.abs(p.a));
      const ref = withArea[0].l.pts, mu = meanOf(ref, 'x'), mv = meanOf(ref, 'y');
      return withArea.map((q) => shifted(q.l.pts, TAU * Math.round((mu - meanOf(q.l.pts, 'x')) / TAU), surf.perV ? TAU * Math.round((mv - meanOf(q.l.pts, 'y')) / TAU) : 0));
    }
    return rectangle(flat.map((l) => l.pts));
  }

  // work with "along" = the winding direction; swap coordinates for loops that wind in v
  const alongV = !!wind[0].wv;
  if (wind.some((l) => !!l.wv !== alongV)) return 'its boundary loops wind around different directions of the surface';
  const sw = (L) => (alongV ? L.map((q) => ({ x: q.y, y: q.x, p: q.p })) : L);
  const back = (polys) => (alongV ? polys.map(sw) : polys);
  const mkA = (x, y) => (alongV ? { x, y, p: surf.fromParam(y, x) } : mk(x, y));
  const betweenA = (a, b) => (alongV ? between({ x: a.y, y: a.x }, { x: b.y, y: b.x }).map((q) => ({ x: q.y, y: q.x, p: q.p })) : between(a, b));
  let A = sw(wind[0].pts);
  const wA = alongV ? wind[0].wv : wind[0].wu;
  const crossPeriodic = alongV ? surf.perU : surf.perV;
  // side of loop A on which the face lies, as a sign on the cross coordinate
  const side = Math.sign(wA) * (alongV ? -1 : 1) * orient;
  // The seam starts at a vertex of loop A and runs into the face. Starting where A reaches furthest into the
  // face keeps the seam clear of A itself (a stepped boundary would otherwise be crossed); of those vertices
  // the first that leaves every inner loop on one side of the seam is taken.
  let lo = 0, holes = null;
  const flatA = flat.map((l) => sw(l.pts));
  const order = A.map((_, k) => k).sort((p, q) => side * (A[q].y - A[p].y) || p - q);
  for (const k of order) {
    if (holes) break;
    const start = A[k].x;
    lo = Math.min(start, start + TAU * wA);
    const hs = flatA.map((h) => shifted(h, -TAU * Math.floor((meanOf(h, 'x') - lo) / TAU), 0));
    if (hs.every((h) => h.every((q) => q.x > lo + 1e-6 && q.x < lo + TAU - 1e-6))) { holes = hs; if (k) A = A.slice(k).concat(A.slice(0, k).map((q) => copy(q, TAU * wA, 0))); }
  }
  if (!holes) return 'its inner loops leave no place for the seam of the surface';
  const a0 = A[0], aEnd = copy(a0, TAU * wA, 0);

  if (wind.length === 1) {
    if (alongV) return 'a single boundary loop winds around the tube of the surface';
    let vp;
    if (surf.poles) vp = side * Math.PI / 2;
    else if (surf.poleRows) { vp = side > 0 ? surf.poleRows[1] : surf.poleRows[0]; if (vp === null) return 'a single closed boundary leaves the face unbounded along the surface'; } else if (surf.apex !== null && surf.apex !== undefined) vp = surf.apex;
    else return 'a single closed boundary leaves the face unbounded along the axis';
    const pole = mkA(a0.x, vp);
    const seam = surf.straightV ? [] : betweenA(a0, pole);
    const poleEnd = copy(pole, TAU * wA, 0);
    const row = betweenA(poleEnd, pole).map((q) => ({ x: q.x, y: q.y, p: pole.p }));
    return back([[...A, aEnd, ...seam.map((q) => copy(q, TAU * wA, 0)), poleEnd, ...row, pole, ...seam.slice().reverse()], ...holes]);
  }
  if (wind.length !== 2) return `${wind.length} boundary loops wind around the surface`;
  let B = sw(wind[1].pts);
  const wB = alongV ? wind[1].wv : wind[1].wu;
  if (Math.sign(wB) === Math.sign(wA)) B = B.slice().reverse();
  // start B at its vertex nearest (in angle) to the start of A, so that the seam is as short as possible
  let best = 0, bd = Infinity;
  for (let k = 0; k < B.length; k++) { const d = Math.abs(nearest(B[k].x, a0.x) - a0.x); if (d < bd) { bd = d; best = k; } }
  const base = B[best].x;
  // B runs against A: rebuild it continuously from its new start
  const Br = [];
  for (let k = 0; k < B.length; k++) { const q = B[(best + k) % B.length]; let x = q.x; if (k) { x = nearest(q.x, Br[k - 1].x); } Br.push({ x, y: q.y, p: q.p }); }
  let dx = nearest(base, a0.x) - base, dy = 0;
  if (crossPeriodic) { const ca = meanOf(A, 'y'), cb = meanOf(Br, 'y'); let d = side * (cb - ca); d -= TAU * Math.floor(d / TAU); dy = ca + side * d - cb; }
  const Bs = shifted(Br, dx, dy);
  const b0 = Bs[0];
  const seam = betweenA(a0, b0);
  return back([[...A, aEnd, ...seam.map((q) => copy(q, TAU * wA, 0)), ...Bs.map((q) => copy(q, TAU * wA, 0)), b0, ...seam.slice().reverse()], ...holes]);
}

/**
 * Triangulate a polygon with holes by ear clipping, keeping every boundary vertex (collinear ones included) and
 * every boundary segment — which is what lets neighbouring faces share their boundary points. Holes are joined
 * to the outer ring by a bridge between two mutually visible vertices. O(n²).
 * @param {Array<{X: number, Y: number}>} V vertices @param {number[][]} rings vertex indices per ring, outer first
 * @returns {number[]|null} counter-clockwise triangles, or null when the rings touch or cross (a hole cannot be
 *   joined, or no ear is left in a part with area)
 */
function triangulatePolygon(V, rings) {
  const cross = (a, b, c) => (V[b].X - V[a].X) * (V[c].Y - V[a].Y) - (V[b].Y - V[a].Y) * (V[c].X - V[a].X);
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const r of rings) for (const i of r) { if (V[i].X < x0) x0 = V[i].X; if (V[i].X > x1) x1 = V[i].X; if (V[i].Y < y0) y0 = V[i].Y; if (V[i].Y > y1) y1 = V[i].Y; }
  const eps = 1e-13 * Math.max(x1 - x0, y1 - y0) ** 2;
  const areaOfRing = (r) => { let a = 0; for (let i = 0, j = r.length - 1; i < r.length; j = i++) a += V[r[j]].X * V[r[i]].Y - V[r[i]].X * V[r[j]].Y; return a; };
  const list = (r, ccw) => {
    const idx = (areaOfRing(r) > 0) === ccw ? r : r.slice().reverse();
    const nodes = idx.map((i) => ({ i, prev: null, next: null }));
    nodes.forEach((n, k) => { n.next = nodes[(k + 1) % nodes.length]; n.prev = nodes[(k + nodes.length - 1) % nodes.length]; });
    return nodes[0];
  };
  const each = (start, fn) => { let n = start; do { fn(n); n = n.next; } while (n !== start); };
  // p–q is blocked by the edge a–b when they cross or touch (end points they share apart)
  const near = 1e-9 * Math.max(x1 - x0, y1 - y0);
  const onSeg = (p, q, a) => {
    const ex = V[q].X - V[p].X, ey = V[q].Y - V[p].Y, l = Math.hypot(ex, ey);
    if (!(l > 0)) return false;
    const t = ((V[a].X - V[p].X) * ex + (V[a].Y - V[p].Y) * ey) / l;
    return t >= -near && t <= l + near && Math.abs(cross(p, q, a)) <= near * l;
  };
  const blocked = (p, q, a, b) => {
    const shareA = a === p || a === q, shareB = b === p || b === q;
    if (shareA && shareB) return false;
    if ((!shareA && onSeg(p, q, a)) || (!shareB && onSeg(p, q, b))) return true;
    if (shareA || shareB) return false;
    if (onSeg(a, b, p) || onSeg(a, b, q)) return true;
    const d1 = cross(p, q, a), d2 = cross(p, q, b), d3 = cross(a, b, p), d4 = cross(a, b, q);
    return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
  };
  // b lies in the interior angle of the (counter-clockwise) polygon at node a
  const inCone = (a, b) => (cross(a.prev.i, a.i, a.next.i) > eps
    ? cross(a.i, b, a.prev.i) > eps && cross(b, a.i, a.next.i) > eps
    : !(cross(a.i, b, a.next.i) >= -eps && cross(b, a.i, a.prev.i) >= -eps));
  let outer = list(rings[0], true);
  const holes = rings.slice(1).map((r) => list(r, false));
  const rightmost = (h) => { let m = h; each(h, (n) => { if (V[n.i].X > V[m.i].X) m = n; }); return m; };
  const pending = holes.map(rightmost).sort((a, b) => V[b.i].X - V[a.i].X);
  for (let hi = 0; hi < pending.length; hi++) {
    const h = pending[hi];
    const cand = [];
    each(outer, (n) => cand.push(n));
    const d2 = (n) => (V[n.i].X - V[h.i].X) ** 2 + (V[n.i].Y - V[h.i].Y) ** 2;
    cand.sort((a, b) => d2(a) - d2(b));
    let found = null;
    for (const o of cand) {
      if (!inCone(o, h.i)) continue;
      let ok = true;
      const test = (n) => { if (ok && blocked(h.i, o.i, n.i, n.next.i)) ok = false; };
      each(outer, test);
      for (let k = hi; k < pending.length && ok; k++) each(pending[k], test);
      if (ok) { found = o; break; }
    }
    if (!found) return null;
    // … o → h → (hole) → h′ → o′ → o.next …
    const o2 = { i: found.i, prev: null, next: found.next }, h2 = { i: h.i, prev: h.prev, next: o2 };
    found.next.prev = o2; o2.prev = h2;
    h.prev.next = h2;
    found.next = h; h.prev = found;
  }
  const T = [];
  let count = 0;
  each(outer, () => count++);
  const inTri = (a, b, c, p, slack) => cross(a, b, p) >= slack && cross(b, c, p) >= slack && cross(c, a, p) >= slack;
  const samePt = (a, b) => V[a].X === V[b].X && V[a].Y === V[b].Y;
  const isEar = (n, slack) => {
    const a = n.prev.i, b = n.i, c = n.next.i;
    if (!(cross(a, b, c) > eps)) return false;
    for (let p = n.next.next; p !== n.prev; p = p.next) {
      const q = p.i;
      if (q === a || q === b || q === c || samePt(q, a) || samePt(q, b) || samePt(q, c)) continue;
      if (inTri(a, b, c, q, slack)) return false;
    }
    return true;
  };
  let n = outer, idle = 0, slack = -eps;
  while (count > 3) {
    if (isEar(n, slack)) {
      T.push(n.prev.i, n.i, n.next.i);
      n.prev.next = n.next; n.next.prev = n.prev;
      n = n.next; count--; idle = 0; slack = -eps;
      continue;
    }
    n = n.next;
    if (++idle > count) {
      // no ear under the strict test (points on the border of a candidate count as inside): accept points on
      // the border, then give up on what is left (collinear remainders)
      if (slack < 0) { slack = eps; idle = 0; continue; }
      // what is left is either a sliver without area (collinear points) or a part that is not a simple polygon
      let left = 0;
      each(n, (q) => { left += V[q.i].X * V[q.next.i].Y - V[q.next.i].X * V[q.i].Y; });
      if (Math.abs(left) > 1e-7 * Math.abs(areaOfRing(rings[0]))) return null;
      break;
    }
  }
  if (count === 3 && cross(n.prev.i, n.i, n.next.i) > eps) T.push(n.prev.i, n.i, n.next.i);
  // the triangles must cover the area between the rings exactly; they do not when a ring crosses itself or another
  let want = 0, got = 0;
  rings.forEach((r, k) => { want += (k ? -1 : 1) * Math.abs(areaOfRing(r)); });
  for (let t = 0; t < T.length; t += 3) got += cross(T[t], T[t + 1], T[t + 2]);
  if (Math.abs(got - want) > 1e-6 * Math.abs(want)) return null;
  return T;
}

/**
 * Lawson edge flips towards the constrained Delaunay triangulation of the given vertices: an interior edge is
 * flipped when the opposite vertex of one of its triangles lies inside the circumcircle of the other and the two
 * triangles form a strictly convex quadrilateral. Boundary edges (`fixed`) are never flipped.
 * @param {number[]} T triangles (counter-clockwise), changed in place @param {Array<{X, Y}>} V vertices
 */
function flipToDelaunay(T, V, fixed, key) {
  const area2 = (a, b, c) => (V[b].X - V[a].X) * (V[c].Y - V[a].Y) - (V[b].Y - V[a].Y) * (V[c].X - V[a].X);
  let scale = 0;
  for (let t = 0; t < T.length; t += 3) scale = Math.max(scale, Math.abs(area2(T[t], T[t + 1], T[t + 2])));
  const eps = 1e-12 * scale;
  // d lies inside the circumcircle of the counter-clockwise triangle (a, b, c): the angles that see the edge
  // a–b from c and from d add up to more than π (written with sines and cosines, so it is free of scale)
  const inCircle = (a, b, c, d) => {
    const ax = V[a].X - V[c].X, ay = V[a].Y - V[c].Y, bx = V[b].X - V[c].X, by = V[b].Y - V[c].Y;
    const px = V[b].X - V[d].X, py = V[b].Y - V[d].Y, qx = V[a].X - V[d].X, qy = V[a].Y - V[d].Y;
    const s1 = ax * by - ay * bx, c1 = ax * bx + ay * by, s2 = px * qy - py * qx, c2 = px * qx + py * qy;
    return s1 * c2 + c1 * s2 < -1e-12 * Math.sqrt((ax * ax + ay * ay) * (bx * bx + by * by) * (px * px + py * py) * (qx * qx + qy * qy));
  };
  // neighbour across the edge that starts at corner k of triangle t
  const N = new Int32Array(T.length).fill(-1);
  const open = new Map();
  for (let e = 0; e < T.length; e++) {
    const t = e - (e % 3), kk = key(T[e], T[t + (e + 1) % 3]);
    const o = open.get(kk);
    if (o === undefined) open.set(kk, e); else { N[e] = (o - (o % 3)) / 3; N[o] = t / 3; open.delete(kk); }
  }
  const stack = [];
  for (let e = 0; e < T.length; e++) if (N[e] > (e - (e % 3)) / 3) stack.push(e);
  const relink = (n, from, to, a, b) => { if (n < 0) return; for (let q = 0; q < 3; q++) if (N[3 * n + q] === from && T[3 * n + q] === a && T[3 * n + (q + 1) % 3] === b) { N[3 * n + q] = to; return; } };
  // an ear-clipped strip needs of the order of n² flips; the cap only stops a run that numerical noise keeps alive
  for (let guard = 4e7; stack.length && guard > 0; guard--) {
    const e = stack.pop(), k = e % 3, t = (e - k) / 3, u = N[e];
    if (u < 0) continue;
    const a = T[3 * t + k], b = T[3 * t + (k + 1) % 3], c = T[3 * t + (k + 2) % 3];
    if (fixed.has(key(a, b))) continue;
    let k2 = -1;
    for (let q = 0; q < 3; q++) if (T[3 * u + q] === b && T[3 * u + (q + 1) % 3] === a) k2 = q;
    if (k2 < 0) continue;
    const d = T[3 * u + (k2 + 2) % 3];
    if (!inCircle(a, b, c, d) || !(area2(a, d, c) > eps && area2(d, b, c) > eps)) continue;
    const nBC = N[3 * t + (k + 1) % 3], nCA = N[3 * t + (k + 2) % 3], nAD = N[3 * u + (k2 + 1) % 3], nDB = N[3 * u + (k2 + 2) % 3];
    T[3 * t] = a; T[3 * t + 1] = d; T[3 * t + 2] = c;
    N[3 * t] = nAD; N[3 * t + 1] = u; N[3 * t + 2] = nCA;
    T[3 * u] = d; T[3 * u + 1] = b; T[3 * u + 2] = c;
    N[3 * u] = nDB; N[3 * u + 1] = nBC; N[3 * u + 2] = t;
    relink(nAD, u, t, d, a);
    relink(nBC, t, u, c, b);
    stack.push(3 * t, 3 * t + 2, 3 * u, 3 * u + 1);
  }
}

const OFF_SHEET = 'part of its boundary lies on a sheet of the surface that is not handled';

/** Reason returned when the facets of the boundary loops of a face touch or cross: a finer tessellation separates them. */
export const CROSSING = 'its boundary loops touch or cross at this facet size';

/**
 * Triangulate a face domain and lift it to 3-D.
 * @param {Array<Array<{x, y, p}>>} polys loops in the parameter plane with the 3-D point of every vertex
 * @param {{ su: number, sv: number, stepU: number, stepV: number }} sc scaling used for the triangulation and
 *   the largest parameter span one triangle edge may have (Infinity: no limit in that direction)
 * @returns {{ pts: number[][], indices: number[] }|string|null} CROSSING when the loops cannot be told apart
 */
function liftDomain(surf, polys, sc) {
  const rings = [];
  for (const L of polys) {
    const kept = [];
    for (const q of L) {
      const X = q.x * sc.su, Y = q.y * sc.sv;
      if (kept.length && kept[kept.length - 1].X === X && kept[kept.length - 1].Y === Y) continue;
      kept.push({ X, Y, x: q.x, y: q.y, p: q.p });
    }
    if (kept.length > 1 && kept[0].X === kept[kept.length - 1].X && kept[0].Y === kept[kept.length - 1].Y) kept.pop();
    if (kept.length < 3) continue;
    let a = 0;
    for (let i = 0, j = kept.length - 1; i < kept.length; j = i++) a += kept[j].X * kept[i].Y - kept[i].X * kept[j].Y;
    rings.push({ kept, area: Math.abs(a) });
  }
  if (!rings.length) return null;
  // the ring of largest area is the outer boundary, the others are holes in it
  let outer = 0;
  rings.forEach((r, i) => { if (r.area > rings[outer].area) outer = i; });
  if (!(rings[outer].area > 0)) return null;
  const order = [rings[outer], ...rings.filter((_, i) => i !== outer)];
  const V = [], idx = [];
  const boundary = new Set();
  const key = (a, b) => (a < b ? a * 4194304 + b : b * 4194304 + a);
  for (const r of order) {
    const base = V.length;
    r.kept.forEach((q, i) => { V.push(q); boundary.add(key(base + i, base + (i + 1) % r.kept.length)); });
    idx.push(r.kept.map((_, i) => base + i));
  }
  let T = triangulatePolygon(V, idx);
  if (!T) return CROSSING;
  if (!T.length) return null;
  {
    // Every boundary segment must be the edge of a triangle, or the face would not meet its neighbours. Ear
    // clipping may bridge over boundary vertices that lie on a straight stretch of a ring; a triangle whose
    // edge runs along such a stretch is cut at those vertices.
    const used = () => { const u = new Set(); for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) u.add(key(T[t + k], T[t + (k + 1) % 3])); return u; };
    // a segment of a collapsed row (pole, apex) has no length in space and need not be an edge
    const real = (k) => V[Math.floor(k / 4194304)].p !== V[k % 4194304].p;
    let have = used();
    let missing = false;
    for (const k of boundary) if (!have.has(k) && real(k)) { missing = true; break; }
    if (missing) {
      const ringOf = new Int32Array(V.length), posOf = new Int32Array(V.length);
      idx.forEach((r, ri) => r.forEach((v, k) => { ringOf[v] = ri; posOf[v] = k; }));
      const N = [];
      for (let t = 0; t < T.length; t += 3) {
        let tris = [[T[t], T[t + 1], T[t + 2]]];
        for (let e = 0; e < 3; e++) {
          const p = T[t + e], q = T[t + (e + 1) % 3];
          if (boundary.has(key(p, q)) || ringOf[p] !== ringOf[q]) continue;
          const ring = idx[ringOf[p]], n = ring.length;
          // the ring vertices between p and q, walking either way, when all of them lie on the segment p–q
          let chain = null;
          for (const dir of [1, -1]) {
            const steps = (((posOf[q] - posOf[p]) * dir) % n + n) % n;
            if (steps < 2 || steps > n - 2) continue;
            const mids = [];
            let ok = true;
            const ex = V[q].X - V[p].X, ey = V[q].Y - V[p].Y, l2 = ex * ex + ey * ey;
            for (let m = 1; m < steps && ok; m++) {
              const v = ring[(((posOf[p] + dir * m) % n) + n) % n];
              const f = ((V[v].X - V[p].X) * ex + (V[v].Y - V[p].Y) * ey) / l2, off = Math.abs((V[v].X - V[p].X) * ey - (V[v].Y - V[p].Y) * ex);
              if (!(f > 0 && f < 1) || off > 1e-9 * l2 || (mids.length && f <= mids[mids.length - 1].f) || V[v].p === V[p].p || V[v].p === V[q].p) ok = false; else mids.push({ v, f });
            }
            if (ok && mids.length) { chain = [p, ...mids.map((x) => x.v), q]; break; }
          }
          if (!chain) continue;
          const k = tris.findIndex((x) => (x[0] === p && x[1] === q) || (x[1] === p && x[2] === q) || (x[2] === p && x[0] === q));
          const x = tris[k], o = x[0] === p ? x[2] : x[1] === p ? x[0] : x[1];
          tris.splice(k, 1, ...chain.slice(0, -1).map((v, m) => [v, chain[m + 1], o]));
        }
        for (const x of tris) N.push(x[0], x[1], x[2]);
      }
      T = N;
      have = used();
      for (const k of boundary) if (!have.has(k) && real(k)) return CROSSING;
    }
  }
  flipToDelaunay(T, V, boundary, key);
  const limU = Number.isFinite(sc.stepU) ? sc.stepU * (1 + 1e-9) : Infinity, limV = Number.isFinite(sc.stepV) ? sc.stepV * (1 + 1e-9) : Infinity;
  if (limU !== Infinity || limV !== Infinity) {
    const long = (a, b) => Math.abs(V[a].x - V[b].x) > limU || Math.abs(V[a].y - V[b].y) > limV;
    for (let pass = 0; pass < 14 && T.length < 600000; pass++) {
      const mid = new Map();
      const m = (a, b) => {
        const k = key(a, b);
        if (boundary.has(k) || !long(a, b)) return -1;
        let v = mid.get(k);
        if (v === undefined) { const x = (V[a].x + V[b].x) / 2, y = (V[a].y + V[b].y) / 2; v = V.length; V.push({ X: x * sc.su, Y: y * sc.sv, x, y, p: surf.fromParam(x, y) }); mid.set(k, v); }
        return v;
      };
      const N = [];
      for (let t = 0; t < T.length; t += 3) {
        const v = [T[t], T[t + 1], T[t + 2]];
        const e = [m(v[0], v[1]), m(v[1], v[2]), m(v[2], v[0])];
        const n = (e[0] >= 0) + (e[1] >= 0) + (e[2] >= 0);
        if (n === 0) N.push(v[0], v[1], v[2]);
        else if (n === 3) N.push(v[0], e[0], e[2], e[0], v[1], e[1], e[2], e[1], v[2], e[0], e[1], e[2]);
        else if (n === 1) { const k = e[0] >= 0 ? 0 : e[1] >= 0 ? 1 : 2; const a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3]; N.push(a, e[k], c, e[k], b, c); } else {
          // two edges cut: k is the corner between them
          const k = e[0] < 0 ? 2 : e[1] < 0 ? 0 : 1;
          const a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3], mab = e[k], mca = e[(k + 2) % 3];
          N.push(a, mab, mca, mab, b, c, mab, c, mca);
        }
      }
      if (!mid.size) break;
      T = N;
      flipToDelaunay(T, V, boundary, key);
    }
  }
  return { pts: V.map((q) => q.p), indices: T };
}

/**
 * Triangulate one face.
 * @param {object} surf surface from planeSurface / coneSurface / sphereSurface / torusSurface
 * @param {number[][][]} loops boundary loops as closed 3-D polylines (no repeated end point); points are shared
 *   by reference with the neighbouring faces
 * @param {boolean} flip the face normal opposes ∂u × ∂v of the surface parametrisation
 * @param {number} [tol] largest distance a boundary point may lie off the surface; beyond it the face is refused
 * @returns {{ pts: number[][], indices: number[] }|string} triangles (counter-clockwise about ∂u × ∂v), or the
 *   reason the face could not be triangulated
 */
export function tessellateFace(surf, loops, flip, tol = Infinity) {
  let polys, sc;
  if (tol !== Infinity && surf.type !== 'spline') {
    // the boundary of a face lies on its surface; when it does not, the surface was not understood
    let worst = 0;
    for (const L of loops) for (const p of L) { const [u, v] = surf.toParam(p); if (u === u && v === v) { const q = surf.fromParam(u, v); worst = Math.max(worst, Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2])); } }
    if (worst > tol) return `its boundary lies ${worst.toPrecision(2)} model units off the surface`;
  }
  if (surf.type === 'plane') {
    sc = surf.scales();
    polys = loops.map((L) => L.map((p) => { const [x, y] = surf.toParam(p); return { x, y, p }; }));
    if (!polys.length) return 'a planar face has no boundary loops';
  } else if (surf.type === 'spline') {
    sc = surf.scales();
    const mk = (x, y) => ({ x, y, p: surf.fromParam(x, y) });
    if (!loops.length) {
      // no boundary: the whole patch
      // the outline runs through the knot lines, cut to the facet size in between
      const cuts = (knots, step) => { const o = [0]; for (let k = 1; k < knots.length; k++) { const a = knots[k - 1], b = knots[k], n = Math.max(1, Math.ceil((b - a) / step - 1e-9)); for (let i = 1; i <= n; i++) o.push(i === n ? b : a + ((b - a) * i) / n); } o[o.length - 1] = TAU; return o; };
      const us = cuts(surf.knotsU, sc.stepU), vs = cuts(surf.knotsV, sc.stepV), ring = [];
      for (let i = 0; i + 1 < us.length; i++) ring.push(mk(us[i], 0));
      for (let j = 0; j + 1 < vs.length; j++) ring.push(mk(TAU, vs[j]));
      for (let i = us.length - 1; i > 0; i--) ring.push(mk(us[i], TAU));
      for (let j = vs.length - 1; j > 0; j--) ring.push(mk(0, vs[j]));
      polys = [ring];
    } else {
      let worst = 0;
      const tp = surf.toParam;
      const probe = { perU: surf.perU, perV: surf.perV, poleRows: surf.poleRows, poleMid: surf.poleMid, toParam: (p) => { const r = tp(p); if (r[2] > worst) worst = r[2]; return r; } };
      const un = loops.map((L) => unwrapLoop(probe, L, sc.stepU, flip)).filter(Boolean);
      if (un.includes(OFF_SHEET)) return OFF_SHEET;
      if (worst > tol) return `its boundary lies ${worst.toPrecision(2)} model units off the B-spline surface stored for it`;
      if (surf.perU || surf.perV) {
        polys = periodicDomain(surf, un, flip, sc);
        if (typeof polys === 'string') return polys;
      } else polys = un.map((l) => l.pts);
    }
  } else {
    let rmax = 0;
    if (surf.radiusAt) { for (const L of loops) for (const p of L) rmax = Math.max(rmax, surf.radiusAt(surf.toParam(p)[1])); if (!(rmax > 0)) return 'a conical face has no boundary'; }
    sc = surf.scales(rmax);
    const un = loops.map((L) => unwrapLoop(surf, L, sc.stepU, flip)).filter(Boolean);
    if (un.includes(OFF_SHEET)) return OFF_SHEET;
    polys = periodicDomain(surf, un, flip, sc);
    if (typeof polys === 'string') return polys;
  }
  if (sc.squeeze) {
    let v0 = Infinity, v1 = -Infinity;
    for (const L of polys) for (const q of L) { if (q.y < v0) v0 = q.y; if (q.y > v1) v1 = q.y; }
    if (v1 > v0) sc.sv = sc.squeeze / (v1 - v0);
  }
  return liftDomain(surf, polys, sc) || 'its boundary is degenerate';
}

// ---------------------------------------------------------------------------------------------------------
// Mesh of a body: welding, closing, orientation
// ---------------------------------------------------------------------------------------------------------

/**
 * Close T-junctions: wherever a vertex lies on an edge that only one triangle uses, that triangle is split there.
 * They arise where the triangulation of one face subdivided a shared boundary segment and its neighbour did not.
 * @param {Float64Array} P welded coordinates @param {number[]} T triangles (flat) @param {number[]} F face id per triangle
 */
function closeTJunctions(P, T, F, tol) {
  for (let pass = 0; pass < 12; pass++) {
    const nV = P.length / 3;
    const use = new Map();
    for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) {
      const a = T[t + k], b = T[t + (k + 1) % 3];
      const key = a < b ? a * nV + b : b * nV + a;
      const e = use.get(key);
      if (e) e.n++; else use.set(key, { n: 1, t, k });
    }
    const open = [];
    const bv = new Set();
    for (const e of use.values()) if (e.n === 1) { open.push(e); bv.add(T[e.t + e.k]); bv.add(T[e.t + (e.k + 1) % 3]); }
    if (!open.length) return;
    const cand = [...bv];
    const split = new Map(); // triangle offset → { k, mids }
    for (const e of open) {
      if (split.has(e.t)) continue;
      const a = T[e.t + e.k], b = T[e.t + (e.k + 1) % 3];
      const ax = P[3 * a], ay = P[3 * a + 1], az = P[3 * a + 2];
      const ex = P[3 * b] - ax, ey = P[3 * b + 1] - ay, ez = P[3 * b + 2] - az;
      const l2 = ex * ex + ey * ey + ez * ez;
      if (!(l2 > 0)) continue;
      const mids = [];
      for (const c of cand) {
        if (c === a || c === b) continue;
        const dx = P[3 * c] - ax, dy = P[3 * c + 1] - ay, dz = P[3 * c + 2] - az;
        const f = (dx * ex + dy * ey + dz * ez) / l2;
        if (!(f > 0 && f < 1)) continue;
        if (Math.hypot(dx - f * ex, dy - f * ey, dz - f * ez) <= tol && f * Math.sqrt(l2) > tol && (1 - f) * Math.sqrt(l2) > tol) mids.push([f, c]);
      }
      if (mids.length) { mids.sort((p, q) => p[0] - q[0]); split.set(e.t, { k: e.k, mids: mids.map((q) => q[1]) }); }
    }
    if (!split.size) return;
    const nT = [], nF = [];
    for (let t = 0; t < T.length; t += 3) {
      const s = split.get(t);
      if (!s) { nT.push(T[t], T[t + 1], T[t + 2]); nF.push(F[t / 3]); continue; }
      const a = T[t + s.k], b = T[t + (s.k + 1) % 3], o = T[t + (s.k + 2) % 3];
      const chain = [a, ...s.mids, b];
      for (let i = 0; i + 1 < chain.length; i++) { nT.push(chain[i], chain[i + 1], o); nF.push(F[t / 3]); }
    }
    T.length = 0; F.length = 0;
    for (const v of nT) T.push(v);
    for (const v of nF) F.push(v);
  }
}

/** Edge use of a triangle mesh: `closed` when every edge is used by exactly two triangles, `oriented` when in addition the two run it in opposite directions. */
function edgeUse(T, nV) {
  const use = new Map();
  let oriented = true;
  for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) {
    const a = T[t + k], b = T[t + (k + 1) % 3];
    const key = a < b ? a * nV + b : b * nV + a, dir = a < b ? 1 : 4;
    const u = use.get(key) || 0;
    if (u & (dir === 1 ? 3 : 12)) oriented = false;
    use.set(key, u + dir);
  }
  for (const n of use.values()) if (n !== 5) return { closed: (n === 2 || n === 8) && [...use.values()].every((m) => m === 5 || m === 2 || m === 8), oriented: false };
  return { closed: use.size > 0, oriented: oriented && use.size > 0 };
}

function volumeOf(P, T) {
  let s = 0;
  for (let t = 0; t < T.length; t += 3) {
    const a = 3 * T[t], b = 3 * T[t + 1], c = 3 * T[t + 2];
    s += P[a] * (P[b + 1] * P[c + 2] - P[b + 2] * P[c + 1]) - P[a + 1] * (P[b] * P[c + 2] - P[b + 2] * P[c]) + P[a + 2] * (P[b] * P[c + 1] - P[b + 1] * P[c]);
  }
  return s / 6;
}

/**
 * Collects the triangles of the faces of one body and finishes them into one mesh.
 */
export class BodyMesh {
  /** @param {number} diag extent of the body (sets the welding tolerance) */
  constructor(diag) {
    this.diag = diag;
    this.weld = new Welder(diag * 1e-9, 4096);
    this.T = []; this.F = []; this.names = [];
  }
  /**
   * Add the triangles of a face; `flip` reverses their winding. Returns the number of triangles kept.
   * @param {number[][]} pts @param {ArrayLike<number>} idx @param {boolean} flip @param {string} name group name
   */
  addFace(pts, idx, flip, name) {
    const ids = pts.map((p) => this.weld.add(p[0], p[1], p[2]));
    let n = 0;
    for (let i = 0; i < idx.length; i += 3) {
      const a = ids[idx[i]], b = ids[idx[i + (flip ? 2 : 1)]], c = ids[idx[i + (flip ? 1 : 2)]];
      if (a === b || b === c || a === c) continue;
      this.T.push(a, b, c); this.F.push(this.names.length); n++;
    }
    if (n) this.names.push(name);
    return n;
  }
  get empty() { return !this.T.length; }
  /**
   * @param {{ apply?: Function, mirror?: boolean }|null} tf placement of the body
   * @param {number} scale factor applied to the coordinates after the placement
   * @returns {{ coords: Float64Array, indices: Uint32Array, groups: object[], closed: boolean, oriented: boolean, reoriented: boolean, volume: number|null }}
   *   `closed`: every edge is shared by two triangles; `oriented`: and they agree on the outward side; `volume`
   *   of the returned coordinates, null for an open mesh
   */
  finish(tf, scale = 1) {
    const { T, F } = this;
    const P = this.weld.finish();
    closeTJunctions(P, T, F, this.diag * 1e-7);
    const { closed, oriented } = edgeUse(T, P.length / 3);
    let reoriented = false, volume = null;
    if (closed) {
      volume = volumeOf(P, T);
      if (volume < 0) { for (let t = 0; t < T.length; t += 3) { const x = T[t + 1]; T[t + 1] = T[t + 2]; T[t + 2] = x; } reoriented = true; volume = -volume; }
    }
    const mirror = !!(tf && tf.mirror);
    const coords = new Float64Array(P.length);
    for (let i = 0; i < P.length; i += 3) { const p = tf && tf.apply ? tf.apply([P[i], P[i + 1], P[i + 2]]) : [P[i], P[i + 1], P[i + 2]]; coords[i] = p[0] * scale; coords[i + 1] = p[1] * scale; coords[i + 2] = p[2] * scale; }
    // triangles grouped by face, in face order
    const order = Array.from({ length: F.length }, (_, i) => i).sort((a, b) => F[a] - F[b] || a - b);
    const indices = new Uint32Array(T.length);
    const groups = [];
    order.forEach((t, k) => {
      indices[3 * k] = T[3 * t]; indices[3 * k + 1] = T[3 * t + (mirror ? 2 : 1)]; indices[3 * k + 2] = T[3 * t + (mirror ? 1 : 2)];
      const g = groups[groups.length - 1];
      if (g && g.face === F[t]) g.count += 3; else groups.push({ name: this.names[F[t]], start: 3 * k, count: 3, face: F[t] });
    });
    for (const g of groups) delete g.face;
    // the volume of the body as placed (a placement may scale it)
    if (closed) volume = Math.abs(volumeOf(coords, indices));
    return { coords, indices, groups, closed, oriented, reoriented, volume };
  }
}
