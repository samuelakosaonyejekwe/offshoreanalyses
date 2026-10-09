/**
 * Geometric measurement of imported models: areas, volumes, centroids, inertia, bounding boxes, sections,
 * hydrostatic cuts, wall thickness, route and bathymetry statistics.
 *
 * All lengths in metres; surfaces are triangle meshes (positions x, y, z per vertex; indices three per face),
 * oriented outward where a volume is involved. z is up. Grid: z[iy·nx + ix] at (x0 + ix·dx, y0 + iy·dy).
 *
 * References: Mirtich, "Fast and accurate computation of polyhedral mass properties", J. Graphics Tools 1(2)
 * 1996 (divergence-theorem integrals); Tonon, "Explicit exact formulas for the 3-D tetrahedron inertia tensor",
 * J. Math. Stat. 1(1) 2005; Möller & Trumbore, "Fast, minimum storage ray/triangle intersection", 1997.
 */
import { eigSym } from '../engine/math/eigen.js';
import { makeRng } from '../engine/math/rng.js';
import { toSurface, bounds } from './topology.js';

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const unit = (a, name = 'direction') => {
  const l = Math.hypot(a[0], a[1], a[2]);
  if (!(l > 0)) throw new RangeError(`${name}: must be a non-zero vector`);
  return [a[0] / l, a[1] / l, a[2] / l];
};

/** Right-handed in-plane basis (u, v) with u × v = n; u is the global x axis projected into the plane when possible. */
export function planeBasis(normal) {
  const n = unit(normal, 'normal');
  const seed = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  const d = dot(seed, n);
  const u = unit([seed[0] - d * n[0], seed[1] - d * n[1], seed[2] - d * n[2]]);
  return { u, v: cross(n, u), n };
}

/** Total area of a triangle mesh [m²]. */
export function surfaceArea(positions, indices) {
  const P = positions;
  let s = 0;
  for (let f = 0; f < indices.length; f += 3) {
    const a = 3 * indices[f], b = 3 * indices[f + 1], c = 3 * indices[f + 2];
    const ux = P[b] - P[a], uy = P[b + 1] - P[a + 1], uz = P[b + 2] - P[a + 2];
    const vx = P[c] - P[a], vy = P[c + 1] - P[a + 1], vz = P[c + 2] - P[a + 2];
    s += Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
  }
  return s / 2;
}

/**
 * Enclosed volume by the divergence theorem [m³]: positive for a closed, outward-oriented surface.
 * Meaningless for open surfaces (the result then depends on the origin).
 */
export function signedVolume(positions, indices) {
  const P = positions;
  // integrate about the first vertex to limit cancellation for geometry far from the origin
  const ox = P[0] || 0, oy = P[1] || 0, oz = P[2] || 0;
  let s = 0;
  for (let f = 0; f < indices.length; f += 3) {
    const a = 3 * indices[f], b = 3 * indices[f + 1], c = 3 * indices[f + 2];
    const ax = P[a] - ox, ay = P[a + 1] - oy, az = P[a + 2] - oz;
    const bx = P[b] - ox, by = P[b + 1] - oy, bz = P[b + 2] - oz;
    const cx = P[c] - ox, cy = P[c + 1] - oy, cz = P[c + 2] - oz;
    s += ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
  }
  return s / 6;
}

/**
 * Mass properties of the solid enclosed by a closed outward-oriented surface (uniform density).
 * @param {{ density?: number }} [opts] density [kg/m³], default 1 (then mass = volume, inertia in m⁵)
 * @returns {{ volume: number, mass: number, centroid: number[], inertia: number[][], inertiaOrigin: number[][] }}
 *   volume [m³], mass [kg], centroid [m], inertia tensor about the centroid and about the global origin [kg·m²]
 *   (rows/columns x, y, z; products of inertia carry the negative sign: Ixy = −∫xy dm).
 */
export function inertiaSolid(positions, indices, { density = 1 } = {}) {
  const P = positions;
  const o = [P[0] || 0, P[1] || 0, P[2] || 0];
  let vol = 0;
  const m1 = [0, 0, 0];
  const c2 = [0, 0, 0, 0, 0, 0]; // ∫xx, yy, zz, xy, yz, zx about o
  const pr = [[0, 0], [1, 1], [2, 2], [0, 1], [1, 2], [2, 0]];
  for (let f = 0; f < indices.length; f += 3) {
    const a = [0, 0, 0], b = [0, 0, 0], c = [0, 0, 0];
    for (let k = 0; k < 3; k++) { a[k] = P[3 * indices[f] + k] - o[k]; b[k] = P[3 * indices[f + 1] + k] - o[k]; c[k] = P[3 * indices[f + 2] + k] - o[k]; }
    const d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    vol += d / 6;
    for (let k = 0; k < 3; k++) m1[k] += d / 24 * (a[k] + b[k] + c[k]);
    for (let q = 0; q < 6; q++) {
      const i = pr[q][0], j = pr[q][1];
      const si = a[i] + b[i] + c[i], sj = a[j] + b[j] + c[j];
      c2[q] += d / 120 * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + si * sj);
    }
  }
  if (!(Math.abs(vol) > 0)) throw new RangeError('positions/indices: the surface encloses no volume');
  const cg = [m1[0] / vol, m1[1] / vol, m1[2] / vol];
  // central second moments, then tensor
  const cc = c2.map((v, q) => v - vol * cg[pr[q][0]] * cg[pr[q][1]]);
  const tensor = (s, V, g) => {
    // s: second moments about a point; shift not needed here, g unused
    return [
      [s[1] + s[2], -s[3], -s[5]],
      [-s[3], s[0] + s[2], -s[4]],
      [-s[5], -s[4], s[0] + s[1]],
    ].map((r) => r.map((x) => x * density));
  };
  const c = [cg[0] + o[0], cg[1] + o[1], cg[2] + o[2]];
  const so = cc.map((v, q) => v + vol * c[pr[q][0]] * c[pr[q][1]]);
  return { volume: vol, mass: density * vol, centroid: c, inertia: tensor(cc), inertiaOrigin: tensor(so) };
}

/**
 * Mass properties of a thin shell of uniform thickness lying on the surface (plating, hull shell).
 * @param {{ arealDensity?: number, thickness?: number, density?: number }} [opts] mass per area [kg/m²]
 *   (default 1), or thickness [m] × density [kg/m³]
 * @returns {{ area: number, mass: number, centroid: number[], inertia: number[][], inertiaOrigin: number[][] }}
 *   area [m²], mass [kg], area centroid [m], inertia about centroid / origin [kg·m²]. Thin-shell assumption:
 *   thickness ≪ radius of curvature; the t³ term is neglected.
 */
export function inertiaShell(positions, indices, { arealDensity, thickness, density } = {}) {
  const rho = arealDensity ?? (thickness !== undefined && density !== undefined ? thickness * density : 1);
  const P = positions;
  let A = 0;
  const m1 = [0, 0, 0];
  const s = [0, 0, 0, 0, 0, 0];
  const pr = [[0, 0], [1, 1], [2, 2], [0, 1], [1, 2], [2, 0]];
  for (let f = 0; f < indices.length; f += 3) {
    const a = [0, 0, 0], b = [0, 0, 0], c = [0, 0, 0];
    for (let k = 0; k < 3; k++) { a[k] = P[3 * indices[f] + k]; b[k] = P[3 * indices[f + 1] + k]; c[k] = P[3 * indices[f + 2] + k]; }
    const n = cross(sub(b, a), sub(c, a));
    const ar = Math.hypot(n[0], n[1], n[2]) / 2;
    A += ar;
    for (let k = 0; k < 3; k++) m1[k] += ar * (a[k] + b[k] + c[k]) / 3;
    for (let q = 0; q < 6; q++) {
      const i = pr[q][0], j = pr[q][1];
      s[q] += ar / 12 * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + (a[i] + b[i] + c[i]) * (a[j] + b[j] + c[j]));
    }
  }
  if (!(A > 0)) throw new RangeError('positions/indices: the surface has no area');
  const cg = m1.map((v) => v / A);
  const sc = s.map((v, q) => v - A * cg[pr[q][0]] * cg[pr[q][1]]);
  const tensor = (m) => [[m[1] + m[2], -m[3], -m[5]], [-m[3], m[0] + m[2], -m[4]], [-m[5], -m[4], m[0] + m[1]]].map((r) => r.map((x) => x * rho));
  return { area: A, mass: rho * A, centroid: cg, inertia: tensor(sc), inertiaOrigin: tensor(s) };
}

/**
 * Principal moments and axes of a symmetric 3×3 inertia tensor.
 * @returns {{ moments: number[], axes: number[][], radiiOfGyration?: number[] }} moments ascending; axes[k] is
 *   the unit vector of moment k (right-handed set); radii of gyration [m] when `mass` is given.
 */
export function principalAxes(inertia, mass) {
  const e = eigSym(inertia.map((r) => r.slice()));
  const axes = [0, 1, 2].map((k) => unit(e.vectors[k].slice(0, 3)));
  const t = cross(axes[0], axes[1]);
  if (dot(t, axes[2]) < 0) axes[2] = axes[2].map((x) => -x);
  const out = { moments: e.values.slice(0, 3), axes };
  if (mass > 0) out.radiiOfGyration = out.moments.map((m) => Math.sqrt(Math.max(0, m) / mass));
  return out;
}

/** Axis-aligned bounding box: { min, max, size, center, diagonal } [m]. */
export function aabb(positions) {
  const b = bounds(positions);
  return { min: b.min, max: b.max, size: [b.max[0] - b.min[0], b.max[1] - b.min[1], b.max[2] - b.min[2]], center: [(b.min[0] + b.max[0]) / 2, (b.min[1] + b.max[1]) / 2, (b.min[2] + b.max[2]) / 2], diagonal: b.diagonal };
}

/**
 * Oriented bounding box aligned with the principal axes of the vertex cloud (PCA). It is tight for
 * box-like and prismatic bodies but is not the minimum-volume box in general.
 * @returns {{ center: number[], axes: number[][], halfSizes: number[], size: number[], volume: number }}
 *   axes sorted by decreasing extent.
 */
export function obb(positions) {
  const n = (positions.length / 3) | 0;
  if (n === 0) throw new RangeError('positions: empty');
  const m = [0, 0, 0];
  for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) m[k] += positions[3 * i + k];
  for (let k = 0; k < 3; k++) m[k] /= n;
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < n; i++) {
    const d = [positions[3 * i] - m[0], positions[3 * i + 1] - m[1], positions[3 * i + 2] - m[2]];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += d[a] * d[b] / n;
  }
  const e = eigSym(C);
  let axes = [2, 1, 0].map((k) => unit(e.vectors[k].slice(0, 3)));
  if (dot(cross(axes[0], axes[1]), axes[2]) < 0) axes[2] = axes[2].map((x) => -x);
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) {
    const p = [positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]];
    for (let k = 0; k < 3; k++) { const t = dot(p, axes[k]); if (t < lo[k]) lo[k] = t; if (t > hi[k]) hi[k] = t; }
  }
  let half = [0, 1, 2].map((k) => (hi[k] - lo[k]) / 2);
  const mid = [0, 1, 2].map((k) => (hi[k] + lo[k]) / 2);
  const center = [0, 1, 2].map((j) => mid[0] * axes[0][j] + mid[1] * axes[1][j] + mid[2] * axes[2][j]);
  const order = [0, 1, 2].sort((a, b) => half[b] - half[a]);
  axes = order.map((k) => axes[k]); half = order.map((k) => half[k]);
  if (dot(cross(axes[0], axes[1]), axes[2]) < 0) axes[2] = axes[2].map((x) => -x);
  return { center, axes, halfSizes: half, size: half.map((h) => 2 * h), volume: 8 * half[0] * half[1] * half[2] };
}

/**
 * Projected (silhouette) area along a direction — wind or current exposed area — by rasterising the triangles
 * onto the plane normal to the direction and counting covered cells.
 * Accuracy: each cell is tested at its centre, so the error is bounded by about (silhouette perimeter × cell
 * size) / 2 and is typically a fifth of that; at the default 1024 cells along the longer side a compact body
 * is within ~0.1–0.3 %. Double the resolution to halve the error (first-order method).
 * @param {number[]} direction viewing direction (need not be unit)
 * @param {{ resolution?: number }} [opts] cells along the longer side of the projected bounding box (default 1024, max 4096)
 * @returns {{ area: number, cellSize: number, width: number, height: number, centroid: number[], errorBound: number }}
 *   area [m²]; projected extents [m]; centroid of the silhouette in 3-D on the plane through the origin [m];
 *   errorBound = cell size × silhouette pixel perimeter / 2 [m²].
 */
export function projectedArea(positions, indices, direction, { resolution = 1024 } = {}) {
  const { u, v } = planeBasis(direction);
  const res = Math.max(8, Math.min(4096, Math.round(resolution)));
  const n = (positions.length / 3) | 0;
  const X = new Float64Array(n), Y = new Float64Array(n);
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (let i = 0; i < n; i++) {
    const px = positions[3 * i], py = positions[3 * i + 1], pz = positions[3 * i + 2];
    const x = px * u[0] + py * u[1] + pz * u[2], y = px * v[0] + py * v[1] + pz * v[2];
    X[i] = x; Y[i] = y;
    if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
  }
  const w = x1 - x0, h = y1 - y0;
  if (!(w > 0) || !(h > 0)) return { area: 0, cellSize: 0, width: Math.max(w, 0) || 0, height: Math.max(h, 0) || 0, centroid: [0, 0, 0], errorBound: 0 };
  const d = Math.max(w, h) / res;
  const nx = Math.max(1, Math.ceil(w / d)), ny = Math.max(1, Math.ceil(h / d));
  const bm = new Uint8Array(nx * ny);
  for (let f = 0; f < indices.length; f += 3) {
    const a = indices[f], b = indices[f + 1], c = indices[f + 2];
    const ys = [Y[a], Y[b], Y[c]], xs = [X[a], X[b], X[c]];
    const ymin = Math.min(ys[0], ys[1], ys[2]), ymax = Math.max(ys[0], ys[1], ys[2]);
    const j0 = Math.max(0, Math.ceil((ymin - y0) / d - 0.5)), j1 = Math.min(ny - 1, Math.floor((ymax - y0) / d - 0.5));
    for (let j = j0; j <= j1; j++) {
      const yc = y0 + (j + 0.5) * d;
      let xa = Infinity, xb = -Infinity;
      for (let k = 0; k < 3; k++) {
        const k2 = (k + 1) % 3;
        const ya = ys[k], yb = ys[k2];
        if ((ya <= yc && yb >= yc) || (yb <= yc && ya >= yc)) {
          const x = ya === yb ? xs[k] : xs[k] + (yc - ya) / (yb - ya) * (xs[k2] - xs[k]);
          if (x < xa) xa = x; if (x > xb) xb = x;
          if (ya === yb) { if (xs[k2] < xa) xa = xs[k2]; if (xs[k2] > xb) xb = xs[k2]; }
        }
      }
      if (xa > xb) continue;
      const i0 = Math.max(0, Math.ceil((xa - x0) / d - 0.5)), i1 = Math.min(nx - 1, Math.floor((xb - x0) / d - 0.5));
      for (let i = i0; i <= i1; i++) bm[j * nx + i] = 1;
    }
  }
  let count = 0, sx = 0, sy = 0, per = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    if (!bm[j * nx + i]) continue;
    count++; sx += i; sy += j;
    if (i === 0 || !bm[j * nx + i - 1]) per++;
    if (i === nx - 1 || !bm[j * nx + i + 1]) per++;
    if (j === 0 || !bm[(j - 1) * nx + i]) per++;
    if (j === ny - 1 || !bm[(j + 1) * nx + i]) per++;
  }
  const cx = count ? x0 + (sx / count + 0.5) * d : 0, cy = count ? y0 + (sy / count + 0.5) * d : 0;
  return {
    area: count * d * d, cellSize: d, width: w, height: h,
    centroid: [cx * u[0] + cy * v[0], cx * u[1] + cy * v[1], cx * u[2] + cy * v[2]],
    errorBound: per * d * d / 2,
  };
}

/** Area, centroid and second moments of a closed 2-D polygon (signed: counter-clockwise positive). */
function polygonMoments(xy) {
  const n = xy.length / 2;
  let A = 0, Sx = 0, Sy = 0, Ixx = 0, Iyy = 0, Ixy = 0;
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    const x0 = xy[2 * i], y0 = xy[2 * i + 1], x1 = xy[2 * j], y1 = xy[2 * j + 1];
    const c = x0 * y1 - x1 * y0;
    A += c; Sx += (x0 + x1) * c; Sy += (y0 + y1) * c;
    Ixx += (y0 * y0 + y0 * y1 + y1 * y1) * c;
    Iyy += (x0 * x0 + x0 * x1 + x1 * x1) * c;
    Ixy += (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * c;
  }
  return { A: A / 2, Sx: Sx / 6, Sy: Sy / 6, Ixx: Ixx / 12, Iyy: Iyy / 12, Ixy: Ixy / 24 };
}

function pointInLoop(x, y, L) {
  let inside = false;
  const n = L.length / 2;
  for (let i = 0, j = n - 1; i < n; j = i++) {
    const xi = L[2 * i], yi = L[2 * i + 1], xj = L[2 * j], yj = L[2 * j + 1];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/**
 * Area properties of a planar region bounded by closed loops. Holes are recognised by nesting (even–odd rule),
 * so loop winding does not matter.
 * @param {ArrayLike<number>[]} loops packed x, y per loop [m] (first point not repeated)
 * @returns {{ area: number, centroid: number[], Ixx: number, Iyy: number, Ixy: number, I1: number, I2: number,
 *   thetaDeg: number, perimeter: number, rx: number, ry: number }} area [m²]; centroid [m]; second moments about
 *   the centroidal axes [m⁴] (Ixx = ∫y² dA, Iyy = ∫x² dA, Ixy = ∫xy dA); principal moments I1 ≥ I2 and the angle
 *   of axis 1 from x [deg]; perimeter of all loops [m]; radii of gyration [m].
 */
export function sectionProperties(loops) {
  if (!Array.isArray(loops) || !loops.length) throw new RangeError('loops: at least one loop is required');
  const clean = loops.map((L) => {
    let n = (L.length / 2) | 0;
    if (n > 1 && L[0] === L[2 * n - 2] && L[1] === L[2 * n - 1]) n--;
    return Float64Array.from(L).subarray(0, 2 * n);
  }).filter((L) => L.length >= 6);
  let A = 0, Sx = 0, Sy = 0, Ixx = 0, Iyy = 0, Ixy = 0, per = 0;
  clean.forEach((L, li) => {
    const m = polygonMoments(L);
    let depth = 0;
    clean.forEach((M, mi) => { if (mi !== li && pointInLoop(L[0], L[1], M)) depth++; });
    const s = (depth % 2 === 0 ? 1 : -1) * Math.sign(m.A || 1);
    A += s * m.A; Sx += s * m.Sx; Sy += s * m.Sy; Ixx += s * m.Ixx; Iyy += s * m.Iyy; Ixy += s * m.Ixy;
    const n = L.length / 2;
    for (let i = 0; i < n; i++) { const j = (i + 1) % n; per += Math.hypot(L[2 * j] - L[2 * i], L[2 * j + 1] - L[2 * i + 1]); }
  });
  if (!(A > 0)) return { area: 0, centroid: [0, 0], Ixx: 0, Iyy: 0, Ixy: 0, I1: 0, I2: 0, thetaDeg: 0, perimeter: per, rx: 0, ry: 0 };
  const cx = Sx / A, cy = Sy / A;
  const ixx = Ixx - A * cy * cy, iyy = Iyy - A * cx * cx, ixy = Ixy - A * cx * cy;
  const avg = (ixx + iyy) / 2, r = Math.hypot((ixx - iyy) / 2, ixy);
  // axis of the larger principal second moment, measured from x
  const theta = 0.5 * Math.atan2(-2 * ixy, ixx - iyy);
  return { area: A, centroid: [cx, cy], Ixx: ixx, Iyy: iyy, Ixy: ixy, I1: avg + r, I2: avg - r, thetaDeg: theta * 180 / Math.PI, perimeter: per, rx: Math.sqrt(ixx / A), ry: Math.sqrt(iyy / A) };
}

/**
 * Cut a surface with a plane and return the section loops with their area properties (waterplane, frame
 * section, pile cross-section). A vertex exactly on the plane is treated as lying on the positive side, so
 * closed surfaces always give closed loops. For outward-oriented surfaces outer loops run counter-clockwise
 * seen from the tip of the normal.
 * @param {{ point: number[], normal: number[] }} plane a point on the plane [m] and its normal
 * @returns {{ loops: Float64Array[], loops2d: Float64Array[], closed: boolean[], basis: { u, v, n }, origin: number[],
 *   area: number, centroid: number[], centroid2d: number[], Iuu: number, Ivv: number, Iuv: number, I1: number, I2: number,
 *   perimeter: number }} loops as packed 3-D points and as (u, v) coordinates from `origin` (the plane point);
 *   area [m²], centroid in 3-D and in plane coordinates, second moments about centroidal axes parallel to u and v
 *   [m⁴] (Iuu = ∫v² dA resists bending about the u axis). Properties are taken from closed loops only.
 */
export function planeCut(positions, indices, plane) {
  if (!plane || !plane.point || !plane.normal) throw new RangeError('plane: { point, normal } is required');
  const { u, v, n } = planeBasis(plane.normal);
  const o = plane.point;
  const P = positions;
  const nV = (P.length / 3) | 0;
  const d = new Float64Array(nV);
  for (let i = 0; i < nV; i++) d[i] = (P[3 * i] - o[0]) * n[0] + (P[3 * i + 1] - o[1]) * n[1] + (P[3 * i + 2] - o[2]) * n[2];
  const pts = new Map(); // edge key → [u, v]
  const nextOf = new Map();
  const ekey = (a, b) => (a < b ? a * nV + b : b * nV + a);
  const point = (a, b) => {
    const k = ekey(a, b);
    if (!pts.has(k)) {
      const t = d[a] / (d[a] - d[b]);
      const x = P[3 * a] + t * (P[3 * b] - P[3 * a]) - o[0], y = P[3 * a + 1] + t * (P[3 * b + 1] - P[3 * a + 1]) - o[1], z = P[3 * a + 2] + t * (P[3 * b + 2] - P[3 * a + 2]) - o[2];
      pts.set(k, [x * u[0] + y * u[1] + z * u[2], x * v[0] + y * v[1] + z * v[2]]);
    }
    return k;
  };
  for (let f = 0; f < indices.length; f += 3) {
    const t = [indices[f], indices[f + 1], indices[f + 2]];
    const s = [d[t[0]] >= 0, d[t[1]] >= 0, d[t[2]] >= 0];
    const up = s[0] + s[1] + s[2];
    if (up === 0 || up === 3) continue;
    // rotate so that corner 0 is the lone corner on its side
    let k = 0;
    while (!(s[k] !== s[(k + 1) % 3] && s[k] !== s[(k + 2) % 3])) k++;
    const a = t[k], b = t[(k + 1) % 3], c = t[(k + 2) % 3];
    const e1 = point(a, b), e2 = point(c, a);
    // walking with the material on the left seen from +n: the lone corner above → from edge a–b to edge c–a
    if (s[k]) { pushNext(nextOf, e1, e2); } else { pushNext(nextOf, e2, e1); }
  }
  const loops2d = [], closed = [];
  for (const [start] of nextOf) {
    for (;;) {
      const l0 = nextOf.get(start);
      if (!l0 || !l0.length) break;
      const keys = [start];
      let cur = l0.pop(), ok = false;
      for (let guard = 0; guard <= pts.size; guard++) {
        if (cur === start) { ok = true; break; }
        keys.push(cur);
        const nx = nextOf.get(cur);
        if (!nx || !nx.length) break;
        cur = nx.pop();
      }
      const L = new Float64Array(2 * keys.length);
      keys.forEach((k, i) => { const p = pts.get(k); L[2 * i] = p[0]; L[2 * i + 1] = p[1]; });
      loops2d.push(L); closed.push(ok);
    }
  }
  let A = 0, Sx = 0, Sy = 0, Ixx = 0, Iyy = 0, Ixy = 0, per = 0;
  loops2d.forEach((L, i) => {
    const m = L.length / 2;
    for (let q = 0; q + (closed[i] ? 0 : 1) < m; q++) { const j = (q + 1) % m; per += Math.hypot(L[2 * j] - L[2 * q], L[2 * j + 1] - L[2 * q + 1]); }
    if (!closed[i]) return;
    const pm = polygonMoments(L);
    A += pm.A; Sx += pm.Sx; Sy += pm.Sy; Ixx += pm.Ixx; Iyy += pm.Iyy; Ixy += pm.Ixy;
  });
  // an inward-oriented surface gives negative loops; report magnitudes
  const sgn = A < 0 ? -1 : 1;
  A *= sgn; Sx *= sgn; Sy *= sgn; Ixx *= sgn; Iyy *= sgn; Ixy *= sgn;
  const cu = A > 0 ? Sx / A : 0, cv = A > 0 ? Sy / A : 0;
  const Iuu = Ixx - A * cv * cv, Ivv = Iyy - A * cu * cu, Iuv = Ixy - A * cu * cv;
  const avg = (Iuu + Ivv) / 2, r = Math.hypot((Iuu - Ivv) / 2, Iuv);
  const to3 = (x, y) => [o[0] + x * u[0] + y * v[0], o[1] + x * u[1] + y * v[1], o[2] + x * u[2] + y * v[2]];
  const loops = loops2d.map((L) => { const out = new Float64Array(L.length / 2 * 3); for (let i = 0; i < L.length / 2; i++) out.set(to3(L[2 * i], L[2 * i + 1]), 3 * i); return out; });
  return { loops, loops2d, closed, basis: { u, v, n }, origin: [o[0], o[1], o[2]], area: A, centroid: to3(cu, cv), centroid2d: [cu, cv], Iuu, Ivv, Iuv, I1: avg + r, I2: avg - r, perimeter: per };
}

function pushNext(map, a, b) { const l = map.get(a); if (l) l.push(b); else map.set(a, [b]); }

/**
 * Volume and centroid of the part of a closed outward-oriented surface on the negative side of a plane
 * (displaced volume below a waterplane, tank contents below a liquid level). The faces are clipped to the
 * half-space and integrated about a point on the plane, so the cap contributes nothing and need not be built.
 * @param {{ point: number[], normal: number[] }} plane normal points to the excluded side (up for a waterplane)
 * @returns {{ volume: number, centroid: number[], wettedArea: number }} volume [m³], its centroid (centre of
 *   buoyancy) [m], area of the clipped surface (wetted surface) [m²]
 */
export function volumeBelowPlane(positions, indices, plane) {
  if (!plane || !plane.point || !plane.normal) throw new RangeError('plane: { point, normal } is required');
  const n = unit(plane.normal, 'plane.normal');
  const o = plane.point;
  const P = positions;
  let vol = 0, area = 0;
  const m = [0, 0, 0];
  const poly = [];
  const addTri = (a, b, c) => {
    const dd = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    vol += dd / 6;
    for (let k = 0; k < 3; k++) m[k] += dd / 24 * (a[k] + b[k] + c[k]);
    const cr = cross(sub(b, a), sub(c, a));
    area += Math.hypot(cr[0], cr[1], cr[2]) / 2;
  };
  for (let f = 0; f < indices.length; f += 3) {
    const p = [], dist = [];
    for (let k = 0; k < 3; k++) {
      const i = 3 * indices[f + k];
      const q = [P[i] - o[0], P[i + 1] - o[1], P[i + 2] - o[2]];
      p.push(q); dist.push(dot(q, n));
    }
    if (dist[0] >= 0 && dist[1] >= 0 && dist[2] >= 0) continue;
    if (dist[0] <= 0 && dist[1] <= 0 && dist[2] <= 0) { addTri(p[0], p[1], p[2]); continue; }
    poly.length = 0;
    for (let k = 0; k < 3; k++) {
      const k2 = (k + 1) % 3;
      if (dist[k] <= 0) poly.push(p[k]);
      if ((dist[k] < 0 && dist[k2] > 0) || (dist[k] > 0 && dist[k2] < 0)) {
        const t = dist[k] / (dist[k] - dist[k2]);
        poly.push([p[k][0] + t * (p[k2][0] - p[k][0]), p[k][1] + t * (p[k2][1] - p[k][1]), p[k][2] + t * (p[k2][2] - p[k][2])]);
      }
    }
    for (let k = 1; k + 1 < poly.length; k++) addTri(poly[0], poly[k], poly[k + 1]);
  }
  const c = Math.abs(vol) > 0 ? [m[0] / vol + o[0], m[1] / vol + o[1], m[2] / vol + o[2]] : [o[0], o[1], o[2]];
  return { volume: vol, centroid: c, wettedArea: area };
}

/**
 * Wall-thickness estimate of a closed, outward-oriented thin-walled body by sampled inward rays: from the
 * centroid of sampled faces a ray is cast against the outward normal and the distance to the first face hit is
 * recorded. For a solid body this is the local chord, not a wall thickness. Cost is samples × faces.
 * @param {{ samples?: number, seed?: number }} [opts] number of sampled faces (default 256, area-weighted,
 *   deterministic for a seed)
 * @returns {{ min: number, max: number, mean: number, median: number, p05: number, p95: number, samples: number,
 *   misses: number, values: Float32Array, faces: Uint32Array }} thickness statistics [m]; `misses` counts rays that
 *   hit nothing (open surface or reversed normals)
 */
export function wallThickness(positions, indices, { samples = 256, seed = 1 } = {}) {
  const P = positions;
  const F = (indices.length / 3) | 0;
  if (F === 0) throw new RangeError('indices: empty surface');
  const rng = makeRng(seed);
  const cum = new Float64Array(F);
  let tot = 0;
  for (let f = 0; f < F; f++) {
    const a = 3 * indices[3 * f], b = 3 * indices[3 * f + 1], c = 3 * indices[3 * f + 2];
    const cr = cross([P[b] - P[a], P[b + 1] - P[a + 1], P[b + 2] - P[a + 2]], [P[c] - P[a], P[c + 1] - P[a + 1], P[c + 2] - P[a + 2]]);
    tot += Math.hypot(cr[0], cr[1], cr[2]) / 2;
    cum[f] = tot;
  }
  const ns = Math.min(Math.max(1, samples | 0), 100000);
  const vals = [], faces = [];
  let misses = 0;
  for (let s = 0; s < ns; s++) {
    const r = rng.next() * tot;
    let lo = 0, hi = F - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (cum[mid] < r) lo = mid + 1; else hi = mid; }
    const f = lo;
    const a = 3 * indices[3 * f], b = 3 * indices[3 * f + 1], c = 3 * indices[3 * f + 2];
    const e1 = [P[b] - P[a], P[b + 1] - P[a + 1], P[b + 2] - P[a + 2]], e2 = [P[c] - P[a], P[c + 1] - P[a + 1], P[c + 2] - P[a + 2]];
    const nr = cross(e1, e2);
    const nl = Math.hypot(nr[0], nr[1], nr[2]);
    if (!(nl > 0)) continue;
    const dir = [-nr[0] / nl, -nr[1] / nl, -nr[2] / nl];
    const org = [(P[a] + P[b] + P[c]) / 3, (P[a + 1] + P[b + 1] + P[c + 1]) / 3, (P[a + 2] + P[b + 2] + P[c + 2]) / 3];
    let best = Infinity;
    for (let g = 0; g < F; g++) {
      if (g === f) continue;
      const t = rayTriangle(org, dir, P, indices[3 * g], indices[3 * g + 1], indices[3 * g + 2]);
      if (t > 1e-9 && t < best) best = t;
    }
    if (Number.isFinite(best)) { vals.push(best); faces.push(f); } else misses++;
  }
  const sorted = vals.slice().sort((x, y) => x - y);
  const q = (p) => (sorted.length ? sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(p * (sorted.length - 1))))] : NaN);
  return {
    min: sorted.length ? sorted[0] : NaN, max: sorted.length ? sorted[sorted.length - 1] : NaN,
    mean: sorted.length ? sorted.reduce((x, y) => x + y, 0) / sorted.length : NaN,
    median: q(0.5), p05: q(0.05), p95: q(0.95), samples: sorted.length, misses,
    values: Float32Array.from(vals), faces: Uint32Array.from(faces),
  };
}

/**
 * Ray–triangle intersection (Möller–Trumbore). Returns the ray parameter t (distance for a unit direction) or
 * Infinity when the ray misses. Both faces of the triangle are hit.
 */
export function rayTriangle(org, dir, P, ia, ib, ic) {
  const a = 3 * ia, b = 3 * ib, c = 3 * ic;
  const e1x = P[b] - P[a], e1y = P[b + 1] - P[a + 1], e1z = P[b + 2] - P[a + 2];
  const e2x = P[c] - P[a], e2y = P[c + 1] - P[a + 1], e2z = P[c + 2] - P[a + 2];
  const px = dir[1] * e2z - dir[2] * e2y, py = dir[2] * e2x - dir[0] * e2z, pz = dir[0] * e2y - dir[1] * e2x;
  const det = e1x * px + e1y * py + e1z * pz;
  if (Math.abs(det) < 1e-300) return Infinity;
  const inv = 1 / det;
  const tx = org[0] - P[a], ty = org[1] - P[a + 1], tz = org[2] - P[a + 2];
  const u = (tx * px + ty * py + tz * pz) * inv;
  if (u < 0 || u > 1) return Infinity;
  const qx = ty * e1z - tz * e1y, qy = tz * e1x - tx * e1z, qz = tx * e1y - ty * e1x;
  const v = (dir[0] * qx + dir[1] * qy + dir[2] * qz) * inv;
  if (v < 0 || u + v > 1) return Infinity;
  const t = (e2x * qx + e2y * qy + e2z * qz) * inv;
  return t >= 0 ? t : Infinity;
}

/**
 * Route geometry of a 3-D polyline against chainage.
 * @param {ArrayLike<number>} points packed x, y, z [m]
 * @returns {{ length: number, horizontalLength: number, chainage: Float64Array, curvature: Float64Array,
 *   inclinationDeg: Float64Array, headingDeg: Float64Array, minRadius: number, maxInclinationDeg: number,
 *   zMin: number, zMax: number }} chainage along the line [m] per point; curvature 1/R [1/m] per point from the
 *   circle through each point and its neighbours (0 at the ends); inclination above horizontal [deg] and heading
 *   clockwise from +y (north) [deg] of the segment leaving each point (last point repeats the previous segment);
 *   minimum radius of curvature [m] (Infinity for a straight line).
 */
export function polylineMeasure(points) {
  const n = (points.length / 3) | 0;
  if (n < 2) throw new RangeError('points: at least two points are required');
  const ch = new Float64Array(n), kap = new Float64Array(n), inc = new Float64Array(n), hd = new Float64Array(n);
  let hl = 0, zmin = Infinity, zmax = -Infinity, maxInc = 0;
  for (let i = 0; i < n; i++) { const z = points[3 * i + 2]; if (z < zmin) zmin = z; if (z > zmax) zmax = z; }
  for (let i = 0; i + 1 < n; i++) {
    const dx = points[3 * i + 3] - points[3 * i], dy = points[3 * i + 4] - points[3 * i + 1], dz = points[3 * i + 5] - points[3 * i + 2];
    const h = Math.hypot(dx, dy), s = Math.hypot(h, dz);
    ch[i + 1] = ch[i] + s; hl += h;
    inc[i] = s > 0 ? Math.atan2(dz, h) * 180 / Math.PI : 0;
    hd[i] = h > 0 ? (Math.atan2(dx, dy) * 180 / Math.PI + 360) % 360 : (i > 0 ? hd[i - 1] : 0);
    if (Math.abs(inc[i]) > maxInc) maxInc = Math.abs(inc[i]);
  }
  inc[n - 1] = inc[n - 2]; hd[n - 1] = hd[n - 2];
  let kmax = 0;
  for (let i = 1; i + 1 < n; i++) {
    const a = [points[3 * i] - points[3 * i - 3], points[3 * i + 1] - points[3 * i - 2], points[3 * i + 2] - points[3 * i - 1]];
    const b = [points[3 * i + 3] - points[3 * i], points[3 * i + 4] - points[3 * i + 1], points[3 * i + 5] - points[3 * i + 2]];
    const c = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    const la = Math.hypot(a[0], a[1], a[2]), lb = Math.hypot(b[0], b[1], b[2]), lc = Math.hypot(c[0], c[1], c[2]);
    const cr = cross(a, b);
    const den = la * lb * lc;
    // Menger curvature: 4 × triangle area / product of the sides
    kap[i] = den > 0 ? 2 * Math.hypot(cr[0], cr[1], cr[2]) / den : 0;
    if (kap[i] > kmax) kmax = kap[i];
  }
  return { length: ch[n - 1], horizontalLength: hl, chainage: ch, curvature: kap, inclinationDeg: inc, headingDeg: hd, minRadius: kmax > 0 ? 1 / kmax : Infinity, maxInclinationDeg: maxInc, zMin: zmin, zMax: zmax };
}

/**
 * Bilinear sample of a grid at (x, y) [m]. Returns NaN outside the grid or where any of the four surrounding
 * nodes has no data.
 */
export function gridSample(grid, x, y) {
  const fx = (x - grid.x0) / grid.dx, fy = (y - grid.y0) / grid.dy;
  const eps = 1e-9;
  if (!(fx >= -eps) || !(fy >= -eps) || fx > grid.nx - 1 + eps || fy > grid.ny - 1 + eps) return NaN;
  const i = Math.min(Math.max(0, Math.floor(fx)), Math.max(0, grid.nx - 2)), j = Math.min(Math.max(0, Math.floor(fy)), Math.max(0, grid.ny - 2));
  const tx = grid.nx > 1 ? Math.min(1, Math.max(0, fx - i)) : 0, ty = grid.ny > 1 ? Math.min(1, Math.max(0, fy - j)) : 0;
  const i1 = Math.min(i + 1, grid.nx - 1), j1 = Math.min(j + 1, grid.ny - 1);
  const z00 = grid.z[j * grid.nx + i], z10 = grid.z[j * grid.nx + i1], z01 = grid.z[j1 * grid.nx + i], z11 = grid.z[j1 * grid.nx + i1];
  return (1 - tx) * (1 - ty) * z00 + tx * (1 - ty) * z10 + (1 - tx) * ty * z01 + tx * ty * z11;
}

/**
 * Bathymetry / terrain statistics of a grid.
 * @returns {{ nodes: number, valid: number, noData: number, zMin: number, zMax: number, zMean: number, zStd: number,
 *   depthMin: number, depthMax: number, slopeMeanDeg: number, slopeMaxDeg: number, slopeP95Deg: number,
 *   roughness: number, extentX: number, extentY: number, area: number }} elevations [m, z up]; depths = −z
 *   (positive down) so depthMax = −zMin; slope from central differences [deg]; roughness = standard deviation of
 *   the difference between each node and the mean of its eight neighbours [m]; plan extents [m] and area [m²].
 */
export function gridStats(grid) {
  const { nx, ny, z } = grid;
  if (!(nx > 0) || !(ny > 0) || !z || z.length !== nx * ny) throw new RangeError('grid: nx·ny must match z.length');
  let valid = 0, zmin = Infinity, zmax = -Infinity, sum = 0, sum2 = 0;
  for (let i = 0; i < z.length; i++) { const v = z[i]; if (Number.isFinite(v)) { valid++; sum += v; sum2 += v * v; if (v < zmin) zmin = v; if (v > zmax) zmax = v; } }
  const mean = valid ? sum / valid : NaN;
  const slopes = [];
  let rs = 0, rs2 = 0, rn = 0;
  const dx = Math.abs(grid.dx), dy = Math.abs(grid.dy);
  for (let j = 1; j + 1 < ny; j++) for (let i = 1; i + 1 < nx; i++) {
    const c = j * nx + i;
    const zc = z[c], zw = z[c - 1], ze = z[c + 1], zs = z[c - nx], zn = z[c + nx];
    if (!Number.isFinite(zc + zw + ze + zs + zn)) continue;
    slopes.push(Math.atan(Math.hypot((ze - zw) / (2 * dx), (zn - zs) / (2 * dy))) * 180 / Math.PI);
    const corners = z[c - nx - 1] + z[c - nx + 1] + z[c + nx - 1] + z[c + nx + 1];
    if (!Number.isFinite(corners)) continue;
    const r = zc - (zw + ze + zs + zn + corners) / 8;
    rs += r; rs2 += r * r; rn++;
  }
  slopes.sort((a, b) => a - b);
  const sm = slopes.length ? slopes.reduce((a, b) => a + b, 0) / slopes.length : NaN;
  return {
    nodes: nx * ny, valid, noData: nx * ny - valid, zMin: valid ? zmin : NaN, zMax: valid ? zmax : NaN, zMean: mean,
    zStd: valid ? Math.sqrt(Math.max(0, sum2 / valid - mean * mean)) : NaN,
    depthMin: valid ? -zmax : NaN, depthMax: valid ? -zmin : NaN,
    slopeMeanDeg: sm, slopeMaxDeg: slopes.length ? slopes[slopes.length - 1] : NaN, slopeP95Deg: slopes.length ? slopes[Math.round(0.95 * (slopes.length - 1))] : NaN,
    roughness: rn ? Math.sqrt(Math.max(0, rs2 / rn - (rs / rn) ** 2)) : NaN,
    extentX: (nx - 1) * dx, extentY: (ny - 1) * dy, area: (nx - 1) * dx * (ny - 1) * dy,
  };
}

/**
 * Seabed profile under a route: bilinear samples of the grid along a polyline.
 * @param {object} grid grid model
 * @param {ArrayLike<number>} points packed route points; x, y used for position, z (if `stride` is 3) reported as lineZ
 * @param {{ step?: number, stride?: 2|3 }} [opts] sample spacing along the plan route [m] (default: the smaller
 *   grid spacing); stride of `points` (default 3)
 * @returns {{ chainage: Float64Array, x: Float64Array, y: Float64Array, z: Float64Array, lineZ: Float64Array,
 *   gap: Float64Array, slopeDeg: Float64Array, outside: number }} plan chainage [m], sample positions, seabed
 *   elevation [m] (NaN outside the grid), route elevation and gap = lineZ − z [m] (free-span height where
 *   positive), along-route seabed slope [deg], number of samples outside the grid.
 */
export function gridProfile(grid, points, { step, stride = 3 } = {}) {
  const n = (points.length / stride) | 0;
  if (n < 2) throw new RangeError('points: at least two points are required');
  const ds = step ?? Math.min(Math.abs(grid.dx), Math.abs(grid.dy));
  if (!(ds > 0)) throw new RangeError('step: must be > 0');
  const ch = [], xs = [], ys = [], lz = [];
  let s0 = 0;
  for (let i = 0; i + 1 < n; i++) {
    const ax = points[stride * i], ay = points[stride * i + 1], bx = points[stride * i + stride], by = points[stride * i + stride + 1];
    const az = stride > 2 ? points[stride * i + 2] : NaN, bz = stride > 2 ? points[stride * i + stride + 2] : NaN;
    const L = Math.hypot(bx - ax, by - ay);
    const m = Math.max(1, Math.ceil(L / ds - 1e-9));
    for (let k = i === 0 ? 0 : 1; k <= m; k++) {
      const t = k / m;
      ch.push(s0 + t * L); xs.push(ax + t * (bx - ax)); ys.push(ay + t * (by - ay)); lz.push(az + t * (bz - az));
    }
    s0 += L;
  }
  const N = ch.length;
  const z = new Float64Array(N), gap = new Float64Array(N), slope = new Float64Array(N);
  let outside = 0;
  for (let i = 0; i < N; i++) { z[i] = gridSample(grid, xs[i], ys[i]); if (Number.isNaN(z[i])) outside++; gap[i] = lz[i] - z[i]; }
  for (let i = 0; i < N; i++) {
    const a = Math.max(0, i - 1), b = Math.min(N - 1, i + 1);
    slope[i] = ch[b] > ch[a] ? Math.atan((z[b] - z[a]) / (ch[b] - ch[a])) * 180 / Math.PI : 0;
  }
  return { chainage: Float64Array.from(ch), x: Float64Array.from(xs), y: Float64Array.from(ys), z, lineZ: Float64Array.from(lz), gap, slopeDeg: slope, outside };
}

/**
 * Headline measurements of any model kind, for the geometry workbench summary panel.
 * @returns {object} kind-specific fields: surface → { area, volume (NaN if the result is not meaningful because
 *   the model is not passed through a watertightness check here — use topology.inspect), centroid, bbox };
 *   volume → { volume, boundaryArea, bbox }; frame → { totalLength, steelVolume, bbox }; polyline → polylineMeasure;
 *   grid → gridStats; points → { count, bbox }; section → sectionProperties.
 */
export function measure(model) {
  if (!model) throw new RangeError('model: required');
  switch (model.kind) {
    case 'surface': {
      const area = surfaceArea(model.positions, model.indices);
      const volume = signedVolume(model.positions, model.indices);
      let centroid = null;
      try { centroid = (Math.abs(volume) > 1e-12 * area ** 1.5 ? inertiaSolid(model.positions, model.indices) : inertiaShell(model.positions, model.indices)).centroid; } catch { centroid = null; }
      return { kind: 'surface', area, volume, centroid, bbox: aabb(model.positions), faces: model.indices.length / 3, vertices: model.positions.length / 3 };
    }
    case 'volume': {
      const s = toSurface(model);
      return { kind: 'volume', volume: signedVolume(s.positions, s.indices), boundaryArea: surfaceArea(s.positions, s.indices), bbox: aabb(model.positions), nodes: model.positions.length / 3, cells: model.cells.reduce((t, b) => t + b.conn.length / ({ tet4: 4, tet10: 10, hex8: 8, hex20: 20, wedge6: 6, pyr5: 5 }[b.type] || 1), 0) };
    }
    case 'frame': {
      let L = 0, steel = 0;
      for (const m of model.members) {
        const l = Math.hypot(model.nodes[3 * m.a] - model.nodes[3 * m.b], model.nodes[3 * m.a + 1] - model.nodes[3 * m.b + 1], model.nodes[3 * m.a + 2] - model.nodes[3 * m.b + 2]);
        L += l;
        const s = m.section;
        if (s) steel += l * (s.type === 'tube' ? Math.PI * (s.D - s.t) * s.t : s.A || 0);
      }
      return { kind: 'frame', totalLength: L, steelVolume: steel, bbox: aabb(model.nodes), nodes: model.nodes.length / 3, members: model.members.length };
    }
    case 'polyline': return { kind: 'polyline', ...polylineMeasure(model.points), bbox: aabb(model.points) };
    case 'grid': return { kind: 'grid', ...gridStats(model) };
    case 'points': return { kind: 'points', count: model.positions.length / 3, bbox: aabb(model.positions) };
    case 'section': return { kind: 'section', ...sectionProperties(model.loops) };
    default: throw new RangeError(`model.kind: '${model.kind}' is not recognised`);
  }
}
