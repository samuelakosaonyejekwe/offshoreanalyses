/**
 * Mesh generators: closed, outward-oriented triangle surfaces of common offshore bodies with named face groups,
 * a quality-refined Delaunay mesher for planar sections, and conversions between point clouds, grids and surfaces.
 *
 * Conventions: lengths in metres, z up, waterline at z = 0 for floating bodies. Surfaces follow CONTRACT §4:
 * `{ kind: 'surface', name, positions: Float32Array, indices: Uint32Array, groups: [{ name, start, count }] }`
 * with group ranges in index units. Curved bodies are inscribed polyhedra: their area and volume converge to the
 * analytic values from below at second order in the facet size.
 *
 * References: Ruppert, "A Delaunay refinement algorithm for quality 2-dimensional mesh generation",
 * J. Algorithms 18 (1995); Shewchuk, "Delaunay refinement algorithms for triangular mesh generation",
 * Comput. Geom. 22 (2002); Bowyer (1981) and Watson (1981), Comput. J. 24; Shepard, "A two-dimensional
 * interpolation function for irregularly-spaced data", Proc. ACM Nat. Conf. 1968; Wigley hull: Journée,
 * "Experiments and calculations on 4 Wigley hull forms in head waves", TU Delft report 909 (1992).
 */

const pos = (v, name) => { if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be > 0`); return v; };
const count = (v, name, min = 1) => { const n = Math.round(v); if (!(n >= min)) throw new RangeError(`${name}: must be at least ${min}`); return n; };
const lin = (a, b, n) => { const out = new Float64Array(n + 1); for (let i = 0; i <= n; i++) out[i] = i === n ? b : a + (b - a) * i / n; return out; };
const cells = (length, size, explicit, min = 1) => (explicit !== undefined ? count(explicit, 'divisions', min) : Math.max(min, Math.round(Math.abs(length) / size)));

/** Accumulates vertices and grouped triangles and emits a surface model with contiguous groups. */
class Builder {
  constructor() { this.p = []; this.f = []; this.g = []; this.names = []; this.byName = new Map(); this.keys = new Map(); }

  vertex(x, y, z) { this.p.push(x, y, z); return this.p.length / 3 - 1; }

  /** Vertex shared through a caller-chosen key, so that adjoining patches use the same index. */
  shared(key, x, y, z) {
    let i = this.keys.get(key);
    if (i === undefined) { i = this.vertex(x, y, z); this.keys.set(key, i); }
    return i;
  }

  group(name) {
    let i = this.byName.get(name);
    if (i === undefined) { i = this.names.length; this.names.push(name); this.byName.set(name, i); }
    return i;
  }

  tri(a, b, c, g) { if (a !== b && b !== c && c !== a) { this.f.push(a, b, c); this.g.push(g); } }

  quad(a, b, c, d, g) { this.tri(a, b, c, g); this.tri(a, c, d, g); }

  build(name, extra) {
    const F = this.g.length;
    const order = Array.from({ length: F }, (_, i) => i).sort((i, j) => this.g[i] - this.g[j] || i - j);
    const indices = new Uint32Array(3 * F);
    const groups = [];
    let cur = -1;
    order.forEach((f, n) => {
      indices[3 * n] = this.f[3 * f]; indices[3 * n + 1] = this.f[3 * f + 1]; indices[3 * n + 2] = this.f[3 * f + 2];
      if (this.g[f] !== cur) { cur = this.g[f]; groups.push({ name: this.names[cur], start: 3 * n, count: 0 }); }
      groups[groups.length - 1].count += 3;
    });
    return { kind: 'surface', name, positions: Float32Array.from(this.p), indices, groups, ...extra };
  }
}

/** Lattice box faces between coordinate arrays; `map` warps the lattice, `skipTop(i, j)` leaves top cells open. */
function latticeBox(b, xs, ys, zs, names, { prefix = '', map = null, skipTop = null } = {}) {
  const nx = xs.length - 1, ny = ys.length - 1, nz = zs.length - 1;
  const v = (i, j, k) => {
    const key = `${prefix}${i},${j},${k}`;
    const have = b.keys.get(key);
    if (have !== undefined) return have;
    const q = map ? map(xs[i], ys[j], zs[k]) : [xs[i], ys[j], zs[k]];
    return b.shared(key, q[0], q[1], q[2]);
  };
  const g = Object.fromEntries(Object.entries(names).map(([k, n]) => [k, b.group(n)]));
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    b.quad(v(i, j, 0), v(i, j + 1, 0), v(i + 1, j + 1, 0), v(i + 1, j, 0), g.zm);
    if (!skipTop || !skipTop(i, j)) b.quad(v(i, j, nz), v(i + 1, j, nz), v(i + 1, j + 1, nz), v(i, j + 1, nz), g.zp);
  }
  for (let i = 0; i < nx; i++) for (let k = 0; k < nz; k++) {
    b.quad(v(i, 0, k), v(i + 1, 0, k), v(i + 1, 0, k + 1), v(i, 0, k + 1), g.ym);
    b.quad(v(i, ny, k), v(i, ny, k + 1), v(i + 1, ny, k + 1), v(i + 1, ny, k), g.yp);
  }
  for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) {
    b.quad(v(0, j, k), v(0, j, k + 1), v(0, j + 1, k + 1), v(0, j + 1, k), g.xm);
    b.quad(v(nx, j, k), v(nx, j + 1, k), v(nx, j + 1, k + 1), v(nx, j, k + 1), g.xp);
  }
  return v;
}

/**
 * Surface of revolution about the z axis. The profile is a list of { r, z, group } points walked counter-clockwise
 * in the (r, z) half-plane (bottom of the axis → outside → top of the axis); `group` names the strip that starts
 * at that point. Points with r = 0 become poles. A profile that does not touch the axis must be `closed`.
 */
function revolve(b, profile, n, { closed = false, centre = [0, 0, 0], prefix = 'r' } = {}) {
  const rings = profile.map((p, pi) => {
    if (!(p.r > 0)) return [b.shared(`${prefix}${pi}`, centre[0], centre[1], centre[2] + p.z)];
    const ring = new Array(n);
    for (let i = 0; i < n; i++) {
      const a = 2 * Math.PI * i / n;
      ring[i] = b.shared(`${prefix}${pi},${i}`, centre[0] + p.r * Math.cos(a), centre[1] + p.r * Math.sin(a), centre[2] + p.z);
    }
    return ring;
  });
  const m = profile.length;
  for (let s = 0; s < (closed ? m : m - 1); s++) {
    const r0 = rings[s], r1 = rings[(s + 1) % m];
    const g = b.group(profile[s].group || 'surface');
    for (let i = 0; i < n; i++) {
      const j = (i + 1) % n;
      const a = r0[i % r0.length], bb = r0[j % r0.length], c = r1[j % r1.length], d = r1[i % r1.length];
      b.tri(a, bb, c, g);
      b.tri(a, c, d, g);
    }
  }
  return rings;
}

/** Insert intermediate points so that no profile strip is longer than `size`. */
function refineProfile(profile, size, closed = false) {
  const out = [];
  const m = profile.length;
  for (let s = 0; s < (closed ? m : m - 1); s++) {
    const p = profile[s], q = profile[(s + 1) % m];
    const k = Math.max(1, Math.round(Math.hypot(q.r - p.r, q.z - p.z) / size));
    for (let i = 0; i < k; i++) out.push({ r: p.r + (q.r - p.r) * i / k, z: p.z + (q.z - p.z) * i / k, group: p.group });
  }
  if (!closed) out.push(profile[m - 1]);
  return out;
}

const segmentsFor = (radius, size, explicit) => (explicit !== undefined ? count(explicit, 'segments', 3) : Math.max(8, Math.round(2 * Math.PI * radius / size)));

/** Move a body built along +z from the origin so that its axis runs from point a towards point b. */
function placeAlong(model, a, b) {
  const d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
  const L = Math.hypot(d[0], d[1], d[2]);
  if (!(L > 0)) throw new RangeError('a/b: the two points must differ');
  const w = [d[0] / L, d[1] / L, d[2] / L];
  const seed = Math.abs(w[2]) < 0.999 ? [0, 0, 1] : [1, 0, 0];
  let u = [seed[1] * w[2] - seed[2] * w[1], seed[2] * w[0] - seed[0] * w[2], seed[0] * w[1] - seed[1] * w[0]];
  const ul = Math.hypot(u[0], u[1], u[2]);
  u = [u[0] / ul, u[1] / ul, u[2] / ul];
  const v = [w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0]];
  const P = model.positions;
  for (let i = 0; i < P.length; i += 3) {
    const x = P[i], y = P[i + 1], z = P[i + 2];
    P[i] = a[0] + x * u[0] + y * v[0] + z * w[0]; P[i + 1] = a[1] + x * u[1] + y * v[1] + z * w[1]; P[i + 2] = a[2] + x * u[2] + y * v[2] + z * w[2];
  }
  return model;
}

/**
 * Rectangular box.
 * @param {object} o { size: [lx, ly, lz] [m], center = [0, 0, 0], divisions?: [nx, ny, nz], elementSize? [m], name }
 * @returns {object} surface with groups x-, x+, y-, y+, z-, z+ (12 triangles when undivided)
 */
export function box({ size, center = [0, 0, 0], divisions, elementSize, name = 'box' } = {}) {
  if (!size || size.length !== 3) throw new RangeError('size: [lx, ly, lz] is required');
  size.forEach((s, i) => pos(s, `size[${i}]`));
  const h = elementSize ? pos(elementSize, 'elementSize') : Infinity;
  const n = [0, 1, 2].map((k) => cells(size[k], h, divisions ? divisions[k] : undefined));
  const b = new Builder();
  latticeBox(b, lin(center[0] - size[0] / 2, center[0] + size[0] / 2, n[0]), lin(center[1] - size[1] / 2, center[1] + size[1] / 2, n[1]), lin(center[2] - size[2] / 2, center[2] + size[2] / 2, n[2]),
    { xm: 'x-', xp: 'x+', ym: 'y-', yp: 'y+', zm: 'z-', zp: 'z+' });
  return b.build(name);
}

/**
 * Closed frustum of a cone (or a cone when one radius is zero) with its axis along z.
 * @param {object} o { radiusBottom, radiusTop [m] (one may be 0), height [m], base = [0, 0, 0] (centre of the bottom),
 *   segments? (around), elementSize? [m] (sets segments, height and cap divisions), heightDivisions?, capRings?, name }
 * @returns {object} surface with groups bottom, side, top
 */
export function cone({ radiusBottom, radiusTop = 0, height, base = [0, 0, 0], segments, elementSize, heightDivisions, capRings, name = 'cone' } = {}) {
  if (!(radiusBottom >= 0) || !(radiusTop >= 0) || !(radiusBottom + radiusTop > 0)) throw new RangeError('radiusBottom/radiusTop: must be >= 0 and not both zero');
  pos(height, 'height');
  const R = Math.max(radiusBottom, radiusTop);
  const h = elementSize ? pos(elementSize, 'elementSize') : Infinity;
  const n = segmentsFor(R, h, segments);
  const slant = Math.hypot(height, radiusBottom - radiusTop);
  const nh = cells(slant, h, heightDivisions), nb = cells(radiusBottom, h, capRings), nt = cells(radiusTop, h, capRings);
  const prof = [];
  for (let i = 0; i < nb && radiusBottom > 0; i++) prof.push({ r: radiusBottom * i / nb, z: 0, group: 'bottom' });
  for (let i = 0; i < nh; i++) prof.push({ r: radiusBottom + (radiusTop - radiusBottom) * i / nh, z: height * i / nh, group: 'side' });
  for (let i = 0; i < nt && radiusTop > 0; i++) prof.push({ r: radiusTop * (1 - i / nt), z: height, group: 'top' });
  prof.push({ r: 0, z: height, group: 'top' });
  const b = new Builder();
  revolve(b, prof, n, { centre: base });
  return b.build(name);
}

/**
 * Closed circular cylinder with its axis along z.
 * @param {object} o { radius [m], height [m], base = [0, 0, 0], segments?, elementSize?, heightDivisions?, capRings?, name }
 * @returns {object} surface with groups bottom, side, top
 */
export function cylinder({ radius, height, base, segments, elementSize, heightDivisions, capRings, name = 'cylinder' } = {}) {
  pos(radius, 'radius');
  return cone({ radiusBottom: radius, radiusTop: radius, height, base, segments, elementSize, heightDivisions, capRings, name });
}

/**
 * Sphere by recursive subdivision of an icosahedron (near-uniform triangles).
 * @param {object} o { radius [m], center = [0, 0, 0], subdivisions = 3 (0…7; 20·4ˢ triangles), elementSize? [m], name }
 * @returns {object} surface with one group 'surface'
 */
export function sphere({ radius, center = [0, 0, 0], subdivisions = 3, elementSize, name = 'sphere' } = {}) {
  pos(radius, 'radius');
  let s = subdivisions;
  // icosahedron edge = 1.0515 R, halved by each subdivision
  if (elementSize) s = Math.max(0, Math.min(7, Math.round(Math.log2(1.0515 * radius / pos(elementSize, 'elementSize')))));
  if (!(s >= 0 && s <= 7) || !Number.isInteger(s)) throw new RangeError('subdivisions: integer 0…7');
  const t = (1 + Math.sqrt(5)) / 2;
  let V = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]];
  let F = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
  const norm = (p) => { const l = Math.hypot(p[0], p[1], p[2]); return [p[0] / l, p[1] / l, p[2] / l]; };
  V = V.map(norm);
  for (let it = 0; it < s; it++) {
    const mid = new Map();
    const mp = (a, b) => {
      const key = a < b ? a * 1e9 + b : b * 1e9 + a;
      let i = mid.get(key);
      if (i === undefined) { i = V.length; V.push(norm([V[a][0] + V[b][0], V[a][1] + V[b][1], V[a][2] + V[b][2]])); mid.set(key, i); }
      return i;
    };
    const nf = [];
    for (const [a, b, c] of F) { const ab = mp(a, b), bc = mp(b, c), ca = mp(c, a); nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]); }
    F = nf;
  }
  const b = new Builder();
  const g = b.group('surface');
  for (const p of V) b.vertex(center[0] + radius * p[0], center[1] + radius * p[1], center[2] + radius * p[2]);
  for (const f of F) b.tri(f[0], f[1], f[2], g);
  return b.build(name);
}

/**
 * Closed tube (member) between two points; a different end radius gives a tapered member.
 * @param {object} o { a: [x, y, z], b: [x, y, z] [m], radius [m], radiusB? [m] (at b; default = radius), segments?,
 *   elementSize?, heightDivisions?, capRings?, name }
 * @returns {object} surface with groups bottom (at a), side, top (at b)
 */
export function tube({ a, b, radius, radiusB, segments, elementSize, heightDivisions, capRings, name = 'tube' } = {}) {
  if (!a || !b) throw new RangeError('a/b: end points are required');
  const L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const m = cone({ radiusBottom: pos(radius, 'radius'), radiusTop: radiusB ?? radius, height: pos(L, 'a/b distance'), segments, elementSize, heightDivisions, capRings, name });
  return placeAlong(m, a, b);
}

/**
 * Box-shaped barge with optional raked ends, centred on the origin in plan, keel at z = −draft.
 * Volume = L·B·D − B·(rakeBow + rakeStern)·rakeHeight / 2.
 * @param {object} o { length, beam, depth [m], draft = depth / 2 [m], rakeBow = 0, rakeStern = 0 [m] (horizontal
 *   length of the raked part of the bottom), rakeHeight = depth / 2 [m], elementSize? [m], divisions?: [nx, ny, nz], name }
 * @returns {object} surface with groups stern, bow, port, starboard, bottom, deck (bow at +x, port at +y)
 */
export function barge({ length, beam, depth, draft, rakeBow = 0, rakeStern = 0, rakeHeight, elementSize, divisions, name = 'barge' } = {}) {
  pos(length, 'length'); pos(beam, 'beam'); pos(depth, 'depth');
  const T = draft ?? depth / 2;
  const hr = rakeHeight ?? depth / 2;
  if (!(rakeBow >= 0) || !(rakeStern >= 0) || rakeBow + rakeStern >= length) throw new RangeError('rakeBow/rakeStern: must be >= 0 and together shorter than the length');
  const raked = rakeBow > 0 || rakeStern > 0;
  if (raked && !(hr > 0 && hr <= depth)) throw new RangeError('rakeHeight: must be within (0, depth]');
  const h = elementSize ? pos(elementSize, 'elementSize') : Math.min(length, beam, depth);
  const nx = cells(length, h, divisions ? divisions[0] : undefined), ny = cells(beam, h, divisions ? divisions[1] : undefined);
  let zs;
  if (raked && hr < depth) {
    const n1 = cells(hr, h, divisions ? Math.max(1, Math.round(divisions[2] * hr / depth)) : undefined);
    const n2 = cells(depth - hr, h, divisions ? Math.max(1, divisions[2] - n1) : undefined);
    zs = Float64Array.from([...lin(-T, -T + hr, n1), ...Array.from(lin(-T + hr, depth - T, n2)).slice(1)]);
  } else zs = lin(-T, depth - T, cells(depth, h, divisions ? divisions[2] : undefined));
  const map = raked ? (x, y, z) => {
    const f = Math.max(0, 1 - (z + T) / hr);
    const xs = -length / 2 + rakeStern * f, xb = length / 2 - rakeBow * f;
    return [xs + (x + length / 2) / length * (xb - xs), y, z];
  } : null;
  const b = new Builder();
  latticeBox(b, lin(-length / 2, length / 2, nx), lin(-beam / 2, beam / 2, ny), zs, { xm: 'stern', xp: 'bow', ym: 'starboard', yp: 'port', zm: 'bottom', zp: 'deck' }, { map });
  return b.build(name);
}

/**
 * Wigley parabolic hull, y = ±(B/2)·(1 − (2x/L)²)·(1 − (z/T)²) for −T ≤ z ≤ 0, closed by the waterplane at z = 0.
 * Analytic values: displaced volume 4/9·L·B·T, waterplane area 2/3·L·B, centre of buoyancy at z = −3T/8.
 * @param {object} o { length, beam, draft [m], nx = 40 (along the length), nz = 12 (over the draft), closed = true
 *   (false leaves the waterplane open, as for a panel model of the wetted hull), name }
 * @returns {object} surface with groups hull-starboard (y < 0), hull-port (y > 0) and deck
 */
export function wigleyHull({ length, beam, draft, nx = 40, nz = 12, closed = true, name = 'wigley' } = {}) {
  pos(length, 'length'); pos(beam, 'beam'); pos(draft, 'draft');
  count(nx, 'nx', 2); count(nz, 'nz', 1);
  const b = new Builder();
  const half = (i, j) => {
    const x = -length / 2 + length * i / nx, z = -draft + draft * j / nz;
    const xi = 2 * x / length, zeta = z / draft;
    return [x, i === 0 || i === nx || j === 0 ? 0 : beam / 2 * (1 - xi * xi) * (1 - zeta * zeta), z];
  };
  const v = (side, i, j) => {
    const p = half(i, j);
    return p[1] === 0 ? b.shared(`c${i},${j}`, p[0], 0, p[2]) : b.shared(`${side}${i},${j}`, p[0], side === 'p' ? p[1] : -p[1], p[2]);
  };
  const gs = b.group('hull-starboard'), gp = b.group('hull-port'), gd = b.group('deck');
  for (let i = 0; i < nx; i++) for (let j = 0; j < nz; j++) {
    // the diagonal always runs to the centreline corner nearer the end, so that no triangle lies wholly on the
    // centreplane (it would be shared by both sides and collapse the stem)
    if (2 * i < nx) {
      b.quad(v('p', i, j), v('p', i, j + 1), v('p', i + 1, j + 1), v('p', i + 1, j), gp);
      b.quad(v('s', i, j), v('s', i + 1, j), v('s', i + 1, j + 1), v('s', i, j + 1), gs);
    } else {
      b.quad(v('p', i, j + 1), v('p', i + 1, j + 1), v('p', i + 1, j), v('p', i, j), gp);
      b.quad(v('s', i + 1, j), v('s', i + 1, j + 1), v('s', i, j + 1), v('s', i, j), gs);
    }
  }
  if (closed) for (let i = 0; i < nx; i++) b.quad(v('s', i, nz), v('s', i + 1, nz), v('p', i + 1, nz), v('p', i, nz), gd);
  return b.build(name);
}

/**
 * Spar hull: vertical cylinder with an optional reduced upper column joined by a conical taper (Hywind-type).
 * @param {object} o { diameter [m], draft [m], freeboard = 10 [m], upperDiameter? [m], taperTop?, taperBottom? [m]
 *   (elevations of the taper ends, z up, required with upperDiameter), segments?, elementSize? [m], name }
 * @returns {object} surface with groups keel, lower, taper, upper, top (keel at z = −draft)
 */
export function spar({ diameter, draft, freeboard = 10, upperDiameter, taperTop, taperBottom, segments, elementSize, name = 'spar' } = {}) {
  pos(diameter, 'diameter'); pos(draft, 'draft');
  if (!(freeboard >= 0)) throw new RangeError('freeboard: must be >= 0');
  const R = diameter / 2;
  const h = elementSize ? pos(elementSize, 'elementSize') : diameter / 4;
  let prof;
  if (upperDiameter !== undefined && upperDiameter !== diameter) {
    pos(upperDiameter, 'upperDiameter');
    if (!(taperBottom > -draft && taperTop > taperBottom && taperTop < freeboard)) throw new RangeError('taperBottom/taperTop: need −draft < taperBottom < taperTop < freeboard');
    prof = [{ r: 0, z: -draft, group: 'keel' }, { r: R, z: -draft, group: 'lower' }, { r: R, z: taperBottom, group: 'taper' }, { r: upperDiameter / 2, z: taperTop, group: 'upper' }, { r: upperDiameter / 2, z: freeboard, group: 'top' }, { r: 0, z: freeboard, group: 'top' }];
  } else {
    if (!(freeboard > -draft)) throw new RangeError('freeboard: must be above the keel');
    prof = [{ r: 0, z: -draft, group: 'keel' }, { r: R, z: -draft, group: 'lower' }, { r: R, z: freeboard, group: 'top' }, { r: 0, z: freeboard, group: 'top' }];
  }
  const b = new Builder();
  revolve(b, refineProfile(prof, h), segmentsFor(R, h, segments));
  return b.build(name);
}

/**
 * Monopile: straight or tapered pile, solid envelope (default) or hollow with a wall thickness.
 * The hollow pile is a closed surface of genus 1 (V − E + F = 0): outer wall, inner wall and two end annuli.
 * @param {object} o { diameter [m] (at the toe), zTop, zToe [m] (z up; e.g. interface level and pile toe),
 *   topDiameter? [m], taperTop?, taperBottom? [m] (elevations of the conical section, required with topDiameter),
 *   wallThickness? [m], segments?, elementSize? [m], name }
 * @returns {object} surface with groups toe, shaft, taper, upper, top and, when hollow, inner
 */
export function monopile({ diameter, zTop, zToe, topDiameter, taperTop, taperBottom, wallThickness, segments, elementSize, name = 'monopile' } = {}) {
  pos(diameter, 'diameter');
  if (!(zTop > zToe)) throw new RangeError('zTop/zToe: zTop must be above zToe');
  const R = diameter / 2;
  const h = elementSize ? pos(elementSize, 'elementSize') : diameter / 4;
  const tapered = topDiameter !== undefined && topDiameter !== diameter;
  if (tapered) {
    pos(topDiameter, 'topDiameter');
    if (!(taperBottom > zToe && taperTop > taperBottom && taperTop < zTop)) throw new RangeError('taperBottom/taperTop: need zToe < taperBottom < taperTop < zTop');
  }
  const Rt = tapered ? topDiameter / 2 : R;
  const outer = tapered
    ? [{ r: R, z: zToe, group: 'shaft' }, { r: R, z: taperBottom, group: 'taper' }, { r: Rt, z: taperTop, group: 'upper' }, { r: Rt, z: zTop, group: 'top' }]
    : [{ r: R, z: zToe, group: 'shaft' }, { r: R, z: zTop, group: 'top' }];
  const b = new Builder();
  const n = segmentsFor(R, h, segments);
  if (wallThickness !== undefined) {
    const t = pos(wallThickness, 'wallThickness');
    if (!(t < Math.min(R, Rt))) throw new RangeError('wallThickness: must be smaller than the radius');
    const inner = outer.slice().reverse().map((p) => ({ r: p.r - t, z: p.z, group: 'inner' }));
    inner[inner.length - 1].group = 'toe';
    revolve(b, refineProfile([...outer, ...inner], h, true), n, { closed: true });
  } else {
    revolve(b, refineProfile([{ r: 0, z: zToe, group: 'toe' }, ...outer, { r: 0, z: zTop, group: 'top' }], h), n);
  }
  return b.build(name);
}

/**
 * Twin-pontoon semi-submersible: two rectangular pontoons along x with circular columns standing on them. The
 * column feet are imprinted in the pontoon decks, so the hull is one watertight surface without internal faces.
 * Volume = 2·Lp·Bp·Hp + N·A·Hc with A the area of the column polygon (→ π D²/4).
 * @param {object} o { pontoonLength, pontoonWidth, pontoonHeight [m], pontoonSpacing [m] (centre to centre),
 *   columnDiameter [m] (< pontoonWidth), columnHeight [m] (above the pontoon deck), columnsPerPontoon = 2,
 *   columnSpacing? [m] (centre to centre along x; default spreads them over 70 % of the length), draft [m]
 *   (keel at z = −draft), segments = 16 (around a column, multiple of 4), elementSize? [m], name }
 * @returns {object} surface with groups pontoons, columns, column-tops. The two pontoon-and-column hulls are
 *   separate closed components (bracing and deck are not modelled), so V − E + F = 4.
 */
export function semiSubmersible({ pontoonLength, pontoonWidth, pontoonHeight, pontoonSpacing, columnDiameter, columnHeight, columnsPerPontoon = 2, columnSpacing, draft, segments = 16, elementSize, name = 'semi-submersible' } = {}) {
  const Lp = pos(pontoonLength, 'pontoonLength'), Bp = pos(pontoonWidth, 'pontoonWidth'), Hp = pos(pontoonHeight, 'pontoonHeight');
  const S = pos(pontoonSpacing, 'pontoonSpacing'), D = pos(columnDiameter, 'columnDiameter'), Hc = pos(columnHeight, 'columnHeight');
  const T = draft ?? Hp + Hc / 2;
  const nc = count(columnsPerPontoon, 'columnsPerPontoon');
  if (!(D < Bp)) throw new RangeError('columnDiameter: must be smaller than the pontoon width');
  if (!(S > Bp)) throw new RangeError('pontoonSpacing: pontoons overlap');
  if (segments % 4 !== 0 || segments < 8) throw new RangeError('segments: a multiple of 4, at least 8');
  const m = segments / 4;
  const sx = columnSpacing ?? (nc > 1 ? 0.7 * Lp / (nc - 1) : 0);
  const xc = Array.from({ length: nc }, (_, i) => (i - (nc - 1) / 2) * sx);
  if (xc[0] - Bp / 2 < -Lp / 2 - 1e-9 || (nc > 1 && sx < Bp - 1e-9)) throw new RangeError('columnSpacing: column patches overlap or overhang the pontoon ends');
  const h = elementSize ? pos(elementSize, 'elementSize') : Bp / m;
  // x lattice: free stretches between column patches, each patch m cells of the pontoon width
  const xs = [-Lp / 2];
  const patchStart = [];
  const stretch = (to) => { const from = xs[xs.length - 1]; if (to - from < 1e-9 * Lp) { xs[xs.length - 1] = to; return; } const k = Math.max(1, Math.round((to - from) / h)); for (let i = 1; i <= k; i++) xs.push(i === k ? to : from + (to - from) * i / k); };
  for (const x of xc) {
    stretch(x - Bp / 2);
    patchStart.push(xs.length - 1);
    for (let i = 1; i <= m; i++) xs.push(i === m ? x + Bp / 2 : x - Bp / 2 + Bp * i / m);
  }
  stretch(Lp / 2);
  const zs = lin(-T, -T + Hp, Math.max(1, Math.round(Hp / h)));
  const nzc = Math.max(1, Math.round(Hc / h)), ncap = Math.max(1, Math.round(D / 2 / h));
  const b = new Builder();
  const gCol = b.group('columns'), gPon = b.group('pontoons');
  [-S / 2, S / 2].forEach((yc, pi) => {
    const ys = lin(yc - Bp / 2, yc + Bp / 2, m);
    const inPatch = (i) => patchStart.some((s) => i >= s && i < s + m);
    const v = latticeBox(b, Float64Array.from(xs), ys, zs, { xm: 'pontoons', xp: 'pontoons', ym: 'pontoons', yp: 'pontoons', zm: 'pontoons', zp: 'pontoons' }, { prefix: `p${pi}:`, skipTop: (i) => inPatch(i) });
    const kTop = zs.length - 1;
    xc.forEach((x0, ci) => {
      const s = patchStart[ci];
      // square patch boundary, counter-clockwise from the corner at (+x, −y): matches column angles −45° + 90°·k/m
      const sq = [];
      for (let k = 0; k < m; k++) sq.push(v(s + m, k, kTop));
      for (let k = 0; k < m; k++) sq.push(v(s + m - k, m, kTop));
      for (let k = 0; k < m; k++) sq.push(v(s, m - k, kTop));
      for (let k = 0; k < m; k++) sq.push(v(s + k, 0, kTop));
      const ring = (z, r, tag) => Array.from({ length: 4 * m }, (_, k) => {
        const a = -Math.PI / 4 + Math.PI / 2 * k / m;
        return b.shared(`c${pi},${ci},${tag},${k}`, x0 + r * Math.cos(a), yc + r * Math.sin(a), z);
      });
      let prev = ring(-T + Hp, D / 2, 'z0');
      for (let k = 0; k < 4 * m; k++) { const k1 = (k + 1) % (4 * m); b.quad(sq[k], sq[k1], prev[k1], prev[k], gPon); }
      for (let q = 1; q <= nzc; q++) {
        const next = ring(-T + Hp + Hc * q / nzc, D / 2, `z${q}`);
        for (let k = 0; k < 4 * m; k++) { const k1 = (k + 1) % (4 * m); b.quad(prev[k], prev[k1], next[k1], next[k], gCol); }
        prev = next;
      }
      const gTop = b.group('column-tops');
      for (let q = 1; q <= ncap; q++) {
        const zt = -T + Hp + Hc;
        const next = q === ncap ? [b.shared(`c${pi},${ci},pole`, x0, yc, zt)] : ring(zt, D / 2 * (1 - q / ncap), `t${q}`);
        for (let k = 0; k < 4 * m; k++) {
          const k1 = (k + 1) % (4 * m);
          b.tri(prev[k], prev[k1], next[k1 % next.length], gTop);
          b.tri(prev[k], next[k1 % next.length], next[k % next.length], gTop);
        }
        prev = next;
      }
    });
  });
  return b.build(name);
}

/* --------------------------------------------------------------------------- planar Delaunay triangulation */

const orient = (ax, ay, bx, by, cx, cy) => (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);

/** Incremental Delaunay triangulation with constrained (blocking) sub-segments — Bowyer–Watson kernel. */
class Triangulation {
  constructor(minX, minY, maxX, maxY) {
    const s = Math.max(maxX - minX, maxY - minY, 1e-12), cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
    this.x = [cx - 30 * s, cx + 30 * s, cx];
    this.y = [cy - 20 * s, cy - 20 * s, cy + 40 * s];
    this.tv = [0, 1, 2]; this.tn = [-1, -1, -1]; this.alive = [1]; this.free = [];
    this.vt = [0, 0, 0];
    this.seg = new Map();
    this.eps = 1e-14 * s * s;
    this.last = 0;
    this.mark = [0]; this.stamp = 0;
  }

  segKey(a, b) { return a < b ? a * 67108864 + b : b * 67108864 + a; }

  inCircle(t, px, py) {
    const a = this.tv[3 * t], b = this.tv[3 * t + 1], c = this.tv[3 * t + 2];
    const ax = this.x[a] - px, ay = this.y[a] - py, bx = this.x[b] - px, by = this.y[b] - py, cx = this.x[c] - px, cy = this.y[c] - py;
    return (ax * ax + ay * ay) * (bx * cy - by * cx) - (bx * bx + by * by) * (ax * cy - ay * cx) + (cx * cx + cy * cy) * (ax * by - ay * bx) > 0;
  }

  locate(px, py) {
    let t = this.alive[this.last] ? this.last : this.alive.findIndex((a) => a);
    for (let guard = 0; guard < this.alive.length + 16; guard++) {
      let moved = false;
      for (let k = 0; k < 3; k++) {
        const a = this.tv[3 * t + (k + 1) % 3], b = this.tv[3 * t + (k + 2) % 3];
        if (orient(this.x[a], this.y[a], this.x[b], this.y[b], px, py) < 0 && this.tn[3 * t + k] >= 0) { t = this.tn[3 * t + k]; moved = true; break; }
      }
      if (!moved) return t;
    }
    // walk cycled (degenerate input): fall back to exhaustive search
    for (let q = 0; q < this.alive.length; q++) {
      if (!this.alive[q]) continue;
      let inside = true;
      for (let k = 0; k < 3 && inside; k++) {
        const a = this.tv[3 * q + (k + 1) % 3], b = this.tv[3 * q + (k + 2) % 3];
        if (orient(this.x[a], this.y[a], this.x[b], this.y[b], px, py) < 0) inside = false;
      }
      if (inside) return q;
    }
    return t;
  }

  /**
   * Insert a point. Returns { v, tris } or { blocked: [[a, b], …] } when the cavity is bounded by a sub-segment
   * that the point lies on the wrong side of or (with `guard`) encroaches, or null for a degenerate position.
   */
  insert(px, py, seed, guard = false) {
    const stamp = ++this.stamp;
    const cav = [seed];
    this.mark[seed] = stamp;
    const grow = (from) => {
      for (let q = from; q < cav.length; q++) {
        const t = cav[q];
        for (let k = 0; k < 3; k++) {
          const n = this.tn[3 * t + k];
          if (n < 0 || this.mark[n] === stamp) continue;
          if (this.seg.size && this.seg.has(this.segKey(this.tv[3 * t + (k + 1) % 3], this.tv[3 * t + (k + 2) % 3]))) continue;
          if (this.inCircle(n, px, py)) { this.mark[n] = stamp; cav.push(n); }
        }
      }
    };
    grow(0);
    let edges;
    for (let pass = 0; ; pass++) {
      edges = [];
      let grew = false;
      const before = cav.length;
      const blocked = [];
      for (const t of cav) {
        for (let k = 0; k < 3; k++) {
          const n = this.tn[3 * t + k];
          if (n >= 0 && this.mark[n] === stamp) continue;
          const a = this.tv[3 * t + (k + 1) % 3], b = this.tv[3 * t + (k + 2) % 3];
          const o = orient(this.x[a], this.y[a], this.x[b], this.y[b], px, py);
          const isSeg = this.seg.size > 0 && this.seg.has(this.segKey(a, b));
          if (o <= this.eps) {
            // the point is not strictly inside this edge: widen the cavity unless a sub-segment or the hull stops it
            if (isSeg) { blocked.push([a, b]); continue; }
            if (n < 0 || pass > 64) return null;
            this.mark[n] = stamp; cav.push(n); grew = true;
          } else {
            if (guard && isSeg && (this.x[a] - px) * (this.x[b] - px) + (this.y[a] - py) * (this.y[b] - py) < 0) blocked.push([a, b]);
            edges.push(a, b, n);
          }
        }
      }
      if (blocked.length) return { blocked };
      if (!grew) break;
      grow(before);
    }
    const v = this.x.length;
    this.x.push(px); this.y.push(py); this.vt.push(-1);
    for (const t of cav) { this.alive[t] = 0; this.free.push(t); }
    const startAt = new Map(), endAt = new Map();
    const tris = [];
    for (let e = 0; e < edges.length; e += 3) {
      const a = edges[e], b = edges[e + 1], n = edges[e + 2];
      let t = this.free.pop();
      if (t === undefined) { t = this.alive.length; this.alive.push(1); this.tv.push(0, 0, 0); this.tn.push(-1, -1, -1); this.mark.push(0); }
      this.alive[t] = 1;
      this.tv[3 * t] = a; this.tv[3 * t + 1] = b; this.tv[3 * t + 2] = v;
      this.tn[3 * t + 2] = n;
      if (n >= 0) {
        for (let k = 0; k < 3; k++) {
          const na = this.tv[3 * n + (k + 1) % 3], nb = this.tv[3 * n + (k + 2) % 3];
          if (na === b && nb === a) { this.tn[3 * n + k] = t; break; }
        }
      }
      this.vt[a] = t; this.vt[b] = t;
      startAt.set(a, t); endAt.set(b, t);
      tris.push(t);
    }
    this.vt[v] = tris[0];
    for (const t of tris) {
      this.tn[3 * t] = startAt.get(this.tv[3 * t + 1]) ?? -1;
      this.tn[3 * t + 1] = endAt.get(this.tv[3 * t]) ?? -1;
    }
    this.last = tris[0];
    return { v, tris };
  }

  /** Triangle having the directed edge a → b, or −1 when a and b are not joined. */
  findEdge(a, b) {
    const t0 = this.vt[a];
    let t = t0;
    for (let guard = 0; guard < 10000; guard++) {
      let k = 0;
      while (k < 3 && this.tv[3 * t + k] !== a) k++;
      if (k === 3) return -1;
      if (this.tv[3 * t + (k + 1) % 3] === b) return t;
      t = this.tn[3 * t + (k + 1) % 3];
      if (t < 0 || t === t0) break;
    }
    // the fan was cut by the hull: sweep the other way
    t = t0;
    for (let guard = 0; guard < 10000; guard++) {
      let k = 0;
      while (k < 3 && this.tv[3 * t + k] !== a) k++;
      if (k === 3) return -1;
      if (this.tv[3 * t + (k + 1) % 3] === b) return t;
      t = this.tn[3 * t + (k + 2) % 3];
      if (t < 0 || t === t0) break;
    }
    return -1;
  }
}

/**
 * Delaunay triangulation of a planar point set (unconstrained).
 * @param {ArrayLike<number>} xy packed x, y [m]; coincident points are triangulated once
 * @returns {Uint32Array} triangles, counter-clockwise, indexing the input points. Slivers along the convex hull
 *   whose circumcircle is more than about 30 times the size of the point set may be left out.
 */
export function delaunay(xy) {
  const n = (xy.length / 2) | 0;
  if (n < 3) return new Uint32Array(0);
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (let i = 0; i < n; i++) { const x = xy[2 * i], y = xy[2 * i + 1]; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
  const T = new Triangulation(x0, y0, x1, y1);
  const id = [0, 0, 0];
  const seen = new Map();
  // insert along a coarse grid sweep so that the walking point location stays short
  const order = Array.from({ length: n }, (_, i) => i);
  const g = Math.max(1, Math.round(Math.sqrt(n) / 2));
  const cell = (i) => { const cx = Math.min(g - 1, Math.floor((xy[2 * i] - x0) / ((x1 - x0) || 1) * g)); const cy = Math.min(g - 1, Math.floor((xy[2 * i + 1] - y0) / ((y1 - y0) || 1) * g)); return cx * g + (cx % 2 ? g - 1 - cy : cy); };
  const ck = order.map(cell);
  order.sort((a, b) => ck[a] - ck[b]);
  for (const i of order) {
    const key = `${xy[2 * i]},${xy[2 * i + 1]}`;
    if (seen.has(key)) continue;
    const r = T.insert(xy[2 * i], xy[2 * i + 1], T.locate(xy[2 * i], xy[2 * i + 1]));
    if (r && r.v !== undefined) { seen.set(key, r.v); id[r.v] = i; }
  }
  const out = [];
  for (let t = 0; t < T.alive.length; t++) {
    if (!T.alive[t]) continue;
    const a = T.tv[3 * t], b = T.tv[3 * t + 1], c = T.tv[3 * t + 2];
    if (a > 2 && b > 2 && c > 2) out.push(id[a], id[b], id[c]);
  }
  return Uint32Array.from(out);
}

function inLoops(loops, x, y) {
  let inside = false;
  for (const L of loops) {
    const n = L.length / 2;
    for (let i = 0, j = n - 1; i < n; j = i++) {
      const xi = L[2 * i], yi = L[2 * i + 1], xj = L[2 * j], yj = L[2 * j + 1];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
    }
  }
  return inside;
}

/**
 * Quality triangulation of a planar section (Ruppert's Delaunay refinement). Loops bound the material by the
 * even–odd rule, so holes are simply further loops; loop winding does not matter. The boundary is respected
 * exactly; boundary segments are subdivided as the refinement requires. The result is a Delaunay triangulation
 * of its own vertices: no vertex lies inside the circumcircle of any triangle.
 *
 * Validity: loops must be simple and must not touch or cross each other. Angles of the input below about 30°
 * cannot be improved; triangles pinned in such corners are left as they are and reported in `stats`.
 *
 * @param {object} section CONTRACT §4 section model ({ loops: [Float64Array(2n)] })
 * @param {{ size?: number, minAngle?: number, maxVertices?: number }} [opts] target element edge length [m]
 *   (default: quality refinement only); smallest allowed angle [deg], default 25, at most 33; safety cap on the
 *   number of vertices, default 200 000
 * @returns {object} surface model in the plane z = 0 (counter-clockwise triangles, normal +z) with extra fields
 *   `xy` (Float64Array, full-precision coordinates), `boundary` (Uint32Array, vertex pairs of boundary edges) and
 *   `stats: { vertices, triangles, minAngleDeg, maxEdge, area, truncated }`
 */
export function triangulateSection(section, { size, minAngle = 25, maxVertices = 200000 } = {}) {
  if (!section || !Array.isArray(section.loops) || !section.loops.length) throw new RangeError('section.loops: at least one loop is required');
  if (size !== undefined) pos(size, 'size');
  if (!(minAngle >= 0 && minAngle <= 33)) throw new RangeError('minAngle: must be within 0…33 degrees');
  const loops = section.loops.map((L) => {
    let n = (L.length / 2) | 0;
    if (n > 1 && L[0] === L[2 * n - 2] && L[1] === L[2 * n - 1]) n--;
    const pts = [];
    for (let i = 0; i < n; i++) {
      if (!Number.isFinite(L[2 * i] + L[2 * i + 1])) throw new RangeError('section.loops: coordinates must be finite');
      if (pts.length && pts[pts.length - 2] === L[2 * i] && pts[pts.length - 1] === L[2 * i + 1]) continue;
      pts.push(L[2 * i], L[2 * i + 1]);
    }
    return pts;
  }).filter((L) => L.length >= 6);
  if (!loops.length) throw new RangeError('section.loops: every loop has fewer than three distinct points');
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const L of loops) for (let i = 0; i < L.length; i += 2) { if (L[i] < x0) x0 = L[i]; if (L[i] > x1) x1 = L[i]; if (L[i + 1] < y0) y0 = L[i + 1]; if (L[i + 1] > y1) y1 = L[i + 1]; }
  // work about the centre of the section to keep the in-circle determinants well conditioned
  const ox = (x0 + x1) / 2, oy = (y0 + y1) / 2;
  const local = loops.map((L) => L.map((v, i) => v - (i % 2 ? oy : ox)));
  const span = Math.max(x1 - x0, y1 - y0);
  const T = new Triangulation(x0 - ox, y0 - oy, x1 - ox, y1 - oy);
  const vtype = [0, 0, 0]; // 0 input corner, 1 on a segment, 2 free
  const vseg = [-1, -1, -1];
  const byCoord = new Map();
  const addVertex = (x, y, type, segId) => {
    const key = `${x},${y}`;
    const have = byCoord.get(key);
    if (have !== undefined) return have;
    const r = T.insert(x, y, T.locate(x, y));
    if (!r || r.v === undefined) throw new RangeError('section.loops: loops touch or repeat a point');
    byCoord.set(key, r.v); vtype[r.v] = type; vseg[r.v] = segId;
    return r.v;
  };
  // input segments, pre-divided to the target size; each records the corners it runs between
  const segEnds = [];
  const pending = [];
  for (const L of local) {
    const n = L.length / 2;
    const corner = [];
    for (let i = 0; i < n; i++) corner.push(addVertex(L[2 * i], L[2 * i + 1], 0, -1));
    for (let i = 0; i < n; i++) {
      const j = (i + 1) % n;
      const id = segEnds.length;
      segEnds.push([corner[i], corner[j]]);
      const len = Math.hypot(L[2 * j] - L[2 * i], L[2 * j + 1] - L[2 * i + 1]);
      const k = size ? Math.max(1, Math.ceil(len / size - 1e-9)) : 1;
      let prev = corner[i];
      for (let q = 1; q <= k; q++) {
        const v = q === k ? corner[j] : addVertex(L[2 * i] + (L[2 * j] - L[2 * i]) * q / k, L[2 * i + 1] + (L[2 * j + 1] - L[2 * i + 1]) * q / k, 1, id);
        pending.push([prev, v, id]);
        prev = v;
      }
    }
  }
  for (const [a, b, id] of pending) T.seg.set(T.segKey(a, b), id);
  const segQueue = pending.map(([a, b]) => [a, b, false]);
  const triQueue = [];
  let truncated = false;

  const encroaches = (a, b, v) => v > 2 && (T.x[a] - T.x[v]) * (T.x[b] - T.x[v]) + (T.y[a] - T.y[v]) * (T.y[b] - T.y[v]) < 0;
  const third = (t, a, b) => { for (let k = 0; k < 3; k++) { const v = T.tv[3 * t + k]; if (v !== a && v !== b) return v; } return -1; };

  const splitSegment = (a, b, force) => {
    const key = T.segKey(a, b);
    const id = T.seg.get(key);
    if (id === undefined) return;
    const t = T.findEdge(a, b);
    if (!force && t >= 0) {
      const u = T.findEdge(b, a);
      if (!encroaches(a, b, third(t, a, b)) && (u < 0 || !encroaches(a, b, third(u, a, b)))) return;
    }
    if (T.x.length >= maxVertices) { truncated = true; return; }
    const ax = T.x[a], ay = T.y[a], bx = T.x[b], by = T.y[b];
    const len = Math.hypot(bx - ax, by - ay);
    let f = 0.5;
    // concentric shells about an input corner: neighbouring segments are cut at equal distances from it, which
    // stops two segments meeting at a small angle from splitting each other for ever (Ruppert 1995, §6)
    const ca = vtype[a] === 0, cb = vtype[b] === 0;
    if (ca !== cb) {
      const unit = 2 ** Math.round(Math.log2(len / 2 / span)) * span;
      f = ca ? unit / len : 1 - unit / len;
      if (!(f > 0.2 && f < 0.8)) f = 0.5;
    }
    const mx = ax + f * (bx - ax), my = ay + f * (by - ay);
    T.seg.delete(key);
    const seed = t >= 0 ? t : T.locate(mx, my);
    const r = T.insert(mx, my, seed);
    if (!r || r.v === undefined) { T.seg.set(key, id); return; }
    vtype[r.v] = 1; vseg[r.v] = id;
    T.seg.set(T.segKey(a, r.v), id); T.seg.set(T.segKey(r.v, b), id);
    segQueue.push([a, r.v, false], [r.v, b, false]);
    afterInsert(r);
  };
  const afterInsert = (r) => {
    for (const t of r.tris) {
      const a = T.tv[3 * t], b = T.tv[3 * t + 1];
      if (T.seg.has(T.segKey(a, b)) && encroaches(a, b, r.v)) segQueue.push([a, b, true]);
      triQueue.push(t);
    }
  };
  const drainSegments = () => {
    for (let guard = 0; segQueue.length && guard < 50 * maxVertices; guard++) {
      const [a, b, force] = segQueue.pop();
      splitSegment(a, b, force);
    }
    segQueue.length = 0;
  };
  drainSegments();

  const sinMin = Math.sin(minAngle * Math.PI / 180);
  const rMax = size ? 0.75 * size : Infinity;
  const domain = new Map();
  const isInside = (t) => {
    const a = T.tv[3 * t], b = T.tv[3 * t + 1], c = T.tv[3 * t + 2];
    if (a < 3 || b < 3 || c < 3) return false;
    return inLoops(local, (T.x[a] + T.x[b] + T.x[c]) / 3, (T.y[a] + T.y[b] + T.y[c]) / 3);
  };
  // a skinny triangle whose short edge joins two segments at equal distance from their common corner is the
  // unavoidable consequence of a small input angle and is not refined
  const pinned = (p, q) => {
    if (vtype[p] !== 1 || vtype[q] !== 1 || vseg[p] === vseg[q]) return false;
    const sp = segEnds[vseg[p]], sq = segEnds[vseg[q]];
    const c = sp[0] === sq[0] || sp[0] === sq[1] ? sp[0] : sp[1] === sq[0] || sp[1] === sq[1] ? sp[1] : -1;
    if (c < 0) return false;
    const dp = Math.hypot(T.x[p] - T.x[c], T.y[p] - T.y[c]), dq = Math.hypot(T.x[q] - T.x[c], T.y[q] - T.y[c]);
    return Math.abs(dp - dq) <= 1e-9 * (dp + dq);
  };
  for (let t = 0; t < T.alive.length; t++) if (T.alive[t]) triQueue.push(t);
  for (let guard = 0; triQueue.length && guard < 60 * maxVertices; guard++) {
    const t = triQueue.pop();
    if (!T.alive[t]) continue;
    const a = T.tv[3 * t], b = T.tv[3 * t + 1], c = T.tv[3 * t + 2];
    const sig = a * 1e12 + b * 1e6 + c;
    let ins = domain.get(t);
    if (!ins || ins[0] !== sig) { ins = [sig, isInside(t)]; domain.set(t, ins); }
    if (!ins[1]) continue;
    const ax = T.x[a], ay = T.y[a], bx = T.x[b] - ax, by = T.y[b] - ay, cx = T.x[c] - ax, cy = T.y[c] - ay;
    const d = 2 * (bx * cy - by * cx);
    if (!(d > 0)) continue;
    const b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const ux = (cy * b2 - by * c2) / d, uy = (bx * c2 - cx * b2) / d;
    const R = Math.hypot(ux, uy);
    const lab = Math.sqrt(b2), lac = Math.sqrt(c2), lbc = Math.hypot(cx - bx, cy - by);
    const lmin = Math.min(lab, lac, lbc);
    const tooBig = R > rMax;
    let skinny = lmin / (2 * R) < sinMin;
    if (skinny && !tooBig) {
      const [p, q] = lmin === lab ? [a, b] : lmin === lac ? [a, c] : [b, c];
      if (pinned(p, q)) skinny = false;
    }
    if (!tooBig && !skinny) continue;
    if (T.x.length >= maxVertices) { truncated = true; break; }
    const r = T.insert(ax + ux, ay + uy, t, true);
    if (!r) continue;
    if (r.blocked) {
      for (const [p, q] of r.blocked) segQueue.push([p, q, true]);
      drainSegments();
      if (T.alive[t] && T.tv[3 * t] === a && T.tv[3 * t + 1] === b && T.tv[3 * t + 2] === c) triQueue.push(t);
      continue;
    }
    vtype[r.v] = 2; vseg[r.v] = -1;
    afterInsert(r);
    drainSegments();
  }
  if (triQueue.length) truncated = true;

  const map = new Int32Array(T.x.length).fill(-1);
  const idx = [], xy = [];
  let minSin = 1, maxEdge = 0, area = 0;
  for (let t = 0; t < T.alive.length; t++) {
    if (!T.alive[t] || !isInside(t)) continue;
    const tri = [T.tv[3 * t], T.tv[3 * t + 1], T.tv[3 * t + 2]];
    for (const v of tri) { if (map[v] < 0) { map[v] = xy.length / 2; xy.push(T.x[v] + ox, T.y[v] + oy); } idx.push(map[v]); }
    const [a, b, c] = tri;
    const l = [Math.hypot(T.x[b] - T.x[a], T.y[b] - T.y[a]), Math.hypot(T.x[c] - T.x[b], T.y[c] - T.y[b]), Math.hypot(T.x[a] - T.x[c], T.y[a] - T.y[c])];
    const A = orient(T.x[a], T.y[a], T.x[b], T.y[b], T.x[c], T.y[c]) / 2;
    area += A;
    maxEdge = Math.max(maxEdge, l[0], l[1], l[2]);
    // smallest angle is opposite the shortest edge: sin θ = 2A / (product of the two longer edges)
    l.sort((p, q) => p - q);
    minSin = Math.min(minSin, 2 * A / (l[1] * l[2]));
  }
  const boundary = [];
  for (const key of T.seg.keys()) {
    const a = Math.floor(key / 67108864), b = key - a * 67108864;
    if (map[a] >= 0 && map[b] >= 0) boundary.push(map[a], map[b]);
  }
  const nV = xy.length / 2;
  const positions = new Float32Array(3 * nV);
  for (let i = 0; i < nV; i++) { positions[3 * i] = xy[2 * i]; positions[3 * i + 1] = xy[2 * i + 1]; }
  return {
    kind: 'surface', name: section.name || 'section', positions, indices: Uint32Array.from(idx),
    xy: Float64Array.from(xy), boundary: Uint32Array.from(boundary),
    stats: { vertices: nV, triangles: idx.length / 3, minAngleDeg: Math.asin(Math.min(1, Math.max(0, minSin))) * 180 / Math.PI, maxEdge, area, truncated },
  };
}

/* ------------------------------------------------------------------------------------ points, grids, surfaces */

/**
 * Grid a scattered point cloud by inverse-distance weighting (Shepard): z = Σ wᵢ zᵢ / Σ wᵢ with wᵢ = dᵢ^−power
 * over the nearest `neighbours` points within `radius`. A node that coincides with a point takes its value.
 * @param {object} points CONTRACT §4 points model (or { positions })
 * @param {{ dx: number, dy?: number, power?: number, neighbours?: number, radius?: number,
 *   bounds?: { x0, y0, x1, y1 }, name?: string }} opts node spacing [m] (dy defaults to dx); power, default 2;
 *   neighbours, default 12; search radius [m], default 5 × the larger spacing (nodes with no point inside it
 *   get NaN); plan extent, default the bounding box of the points
 * @returns {object} grid model (z up, row-major z[iy·nx + ix])
 */
export function pointsToGrid(points, { dx, dy, power = 2, neighbours = 12, radius, bounds, name } = {}) {
  const P = points && (points.positions || points);
  if (!P || P.length < 3) throw new RangeError('points: at least one point is required');
  pos(dx, 'dx');
  const sy = dy ?? dx;
  pos(sy, 'dy');
  pos(power, 'power');
  const kN = count(neighbours, 'neighbours');
  const n = (P.length / 3) | 0;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (let i = 0; i < n; i++) { const x = P[3 * i], y = P[3 * i + 1]; if (!Number.isFinite(x + y + P[3 * i + 2])) continue; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
  if (!(x0 <= x1)) throw new RangeError('points: no finite points');
  const gx0 = bounds ? bounds.x0 : x0, gy0 = bounds ? bounds.y0 : y0, gx1 = bounds ? bounds.x1 : x1, gy1 = bounds ? bounds.y1 : y1;
  const nx = Math.max(1, Math.floor((gx1 - gx0) / dx + 1e-9) + 1), ny = Math.max(1, Math.floor((gy1 - gy0) / sy + 1e-9) + 1);
  if (nx * ny > 5e7) throw new RangeError('dx/dy: the grid would have more than 50 million nodes');
  const rad = radius ?? 5 * Math.max(dx, sy);
  pos(rad, 'radius');
  // bucket the points on a hash grid sized so that a bucket holds a few points
  const bs = Math.max(Math.sqrt((x1 - x0 + dx) * (y1 - y0 + sy) / n) * 1.5, Math.min(dx, sy) * 0.5);
  const bnx = Math.max(1, Math.ceil((x1 - x0) / bs) + 1), bny = Math.max(1, Math.ceil((y1 - y0) / bs) + 1);
  const start = new Uint32Array(bnx * bny + 1);
  const bucket = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    const x = P[3 * i], y = P[3 * i + 1];
    if (!Number.isFinite(x + y + P[3 * i + 2])) continue;
    bucket[i] = Math.floor((y - y0) / bs) * bnx + Math.floor((x - x0) / bs);
    start[bucket[i] + 1]++;
  }
  for (let q = 0; q < bnx * bny; q++) start[q + 1] += start[q];
  const fill = start.slice(0, bnx * bny);
  const order = new Uint32Array(start[bnx * bny]);
  for (let i = 0; i < n; i++) if (bucket[i] >= 0) order[fill[bucket[i]]++] = i;
  const z = new Float32Array(nx * ny);
  const bestD = new Float64Array(kN), bestZ = new Float64Array(kN);
  const maxRing = Math.ceil(rad / bs) + 1;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const gx = gx0 + i * dx, gy = gy0 + j * sy;
    const cx = Math.floor((gx - x0) / bs), cy = Math.floor((gy - y0) / bs);
    let found = 0, worst = rad * rad;
    for (let ring = 0; ring <= maxRing; ring++) {
      // once k points are held, a ring further out than the worst of them cannot improve the set
      if (found === kN && (ring - 1) * bs > Math.sqrt(worst)) break;
      for (let by = cy - ring; by <= cy + ring; by++) {
        if (by < 0 || by >= bny) continue;
        const edge = by === cy - ring || by === cy + ring;
        for (let bx = cx - ring; bx <= cx + ring; bx += edge ? 1 : 2 * ring || 1) {
          if (bx < 0 || bx >= bnx) continue;
          for (let q = start[by * bnx + bx]; q < start[by * bnx + bx + 1]; q++) {
            const p = order[q];
            const d2 = (P[3 * p] - gx) ** 2 + (P[3 * p + 1] - gy) ** 2;
            if (d2 > rad * rad || (found === kN && d2 >= worst)) continue;
            let s = found < kN ? found++ : kN - 1;
            while (s > 0 && bestD[s - 1] > d2) { bestD[s] = bestD[s - 1]; bestZ[s] = bestZ[s - 1]; s--; }
            bestD[s] = d2; bestZ[s] = P[3 * p + 2];
            if (found === kN) worst = bestD[kN - 1];
          }
        }
      }
    }
    if (!found) { z[j * nx + i] = NaN; continue; }
    if (bestD[0] === 0) { z[j * nx + i] = bestZ[0]; continue; }
    let sw = 0, swz = 0;
    for (let s = 0; s < found; s++) { const w = bestD[s] ** (-power / 2); sw += w; swz += w * bestZ[s]; }
    z[j * nx + i] = swz / sw;
  }
  return { kind: 'grid', name: name || (points.name ? `${points.name} (gridded)` : 'grid'), nx, ny, x0: gx0, y0: gy0, dx, dy: sy, z, ...(points.crs ? { crs: points.crs } : {}) };
}

/**
 * Triangulated surface of a grid (two triangles per cell, normals up for positive spacings). Cells touching a
 * no-data node are left out and unused nodes are dropped.
 * @param {object} grid CONTRACT §4 grid model
 * @param {{ stride?: number, name?: string }} [opts] keep every `stride`-th node (default 1)
 * @returns {object} surface model with `gridNode` (Uint32Array): the grid node index of every vertex
 */
export function gridToSurface(grid, { stride = 1, name } = {}) {
  const { nx, ny } = grid;
  if (!(nx > 0) || !(ny > 0) || !grid.z || grid.z.length !== nx * ny) throw new RangeError('grid: nx·ny must match z.length');
  const s = count(stride, 'stride');
  const mx = Math.floor((nx - 1) / s) + 1, my = Math.floor((ny - 1) / s) + 1;
  const id = new Int32Array(mx * my).fill(-1);
  const P = [], src = [], idx = [];
  const v = (i, j) => {
    const q = j * mx + i;
    if (id[q] < 0) { id[q] = src.length; const g = j * s * nx + i * s; src.push(g); P.push(grid.x0 + i * s * grid.dx, grid.y0 + j * s * grid.dy, grid.z[g]); }
    return id[q];
  };
  const ok = (i, j) => Number.isFinite(grid.z[j * s * nx + i * s]);
  const ccw = grid.dx * grid.dy > 0;
  for (let j = 0; j + 1 < my; j++) for (let i = 0; i + 1 < mx; i++) {
    if (!ok(i, j) || !ok(i + 1, j) || !ok(i + 1, j + 1) || !ok(i, j + 1)) continue;
    const a = v(i, j), b = v(i + 1, j), c = v(i + 1, j + 1), d = v(i, j + 1);
    if (ccw) idx.push(a, b, c, a, c, d); else idx.push(a, c, b, a, d, c);
  }
  return { kind: 'surface', name: name || grid.name || 'grid surface', positions: Float32Array.from(P), indices: Uint32Array.from(idx), gridNode: Uint32Array.from(src) };
}
