/**
 * Mesh refinement, coarsening, smoothing and resampling.
 *
 *  - uniform refinement: triangle 1→4, quadrilateral 1→4, tetrahedron 1→8 (Bey), hexahedron 1→8, wedge 1→8,
 *    pyramid 1→6 pyramids + 4 tetrahedra; conforming across mixed cell types because shared edges and faces get
 *    the same new nodes
 *  - `refineTo`: longest-edge bisection (Rivara) of a surface down to a target edge length, or repeated uniform
 *    refinement of a volume mesh; member subdivision for frames; resampling for polylines and grids
 *  - `coarsenTo`: quadric-error edge-collapse decimation of a surface (Garland & Heckbert)
 *  - `smooth`: Laplacian / Taubin smoothing with fixed boundaries
 *
 * Every function returns new arrays and leaves its input untouched. Surface groups are carried through.
 *
 * References: Bey, "Tetrahedral grid refinement", Computing 55 (1995); Rivara, "Algorithms for refining
 * triangular grids suitable for adaptive and multigrid techniques", IJNME 20 (1984); Garland & Heckbert,
 * "Surface simplification using quadric error metrics", SIGGRAPH 97; Taubin, "A signal processing approach to
 * fair surface design", SIGGRAPH 95.
 */
import { CELL_INFO, CELL_EDGES, buildEdges, faceGroupIds, groupsFromIds } from './topology.js';
import { gridSample } from './measure.js';

const LIMIT = 2e7;

/** Shared mid-nodes: one new node per distinct set of parent nodes, placed at their mean (or projected). */
class NodePool {
  constructor(positions, project) {
    this.src = positions;
    this.n0 = (positions.length / 3) | 0;
    this.extra = [];
    this.edge = new Map();
    this.multi = new Map();
    this.project = project;
  }

  coord(i, k) { return i < this.n0 ? this.src[3 * i + k] : this.extra[3 * (i - this.n0) + k]; }

  add(x, y, z) {
    if (this.project) { const q = this.project(x, y, z); x = q[0]; y = q[1]; z = q[2]; }
    this.extra.push(x, y, z);
    return this.n0 + this.extra.length / 3 - 1;
  }

  mid(a, b) {
    const key = a < b ? a * 4294967296 + b : b * 4294967296 + a;
    let i = this.edge.get(key);
    if (i === undefined) {
      i = this.add((this.coord(a, 0) + this.coord(b, 0)) / 2, (this.coord(a, 1) + this.coord(b, 1)) / 2, (this.coord(a, 2) + this.coord(b, 2)) / 2);
      this.edge.set(key, i);
    }
    return i;
  }

  centre(ids) {
    const key = ids.slice().sort((p, q) => p - q).join(',');
    let i = this.multi.get(key);
    if (i === undefined) {
      let x = 0, y = 0, z = 0;
      for (const v of ids) { x += this.coord(v, 0); y += this.coord(v, 1); z += this.coord(v, 2); }
      i = this.add(x / ids.length, y / ids.length, z / ids.length);
      this.multi.set(key, i);
    }
    return i;
  }

  positions() {
    const out = new this.src.constructor(this.src.length + this.extra.length);
    out.set(this.src); out.set(this.extra, this.src.length);
    return out;
  }
}

/**
 * Split every triangle into four (new nodes at the edge mid-points).
 * @param {ArrayLike<number>} positions packed x, y, z [m]
 * @param {ArrayLike<number>} indices three vertex ids per triangle
 * @param {(x: number, y: number, z: number) => number[]} [project] moves each new node onto the true geometry
 * @returns {{ positions, indices: Uint32Array }} child triangles 4f…4f+3 belong to parent f
 */
export function refineTriangles(positions, indices, project) {
  const pool = new NodePool(positions, project);
  const out = new Uint32Array(indices.length * 4);
  for (let f = 0, o = 0; f < indices.length; f += 3) {
    const a = indices[f], b = indices[f + 1], c = indices[f + 2];
    const ab = pool.mid(a, b), bc = pool.mid(b, c), ca = pool.mid(c, a);
    out[o++] = a; out[o++] = ab; out[o++] = ca;
    out[o++] = ab; out[o++] = b; out[o++] = bc;
    out[o++] = ca; out[o++] = bc; out[o++] = c;
    out[o++] = ab; out[o++] = bc; out[o++] = ca;
  }
  return { positions: pool.positions(), indices: out };
}

/**
 * Split every quadrilateral into four (new nodes at the edge mid-points and at the face centre).
 * @param {ArrayLike<number>} quads four vertex ids per quadrilateral, in order around the face
 * @returns {{ positions, quads: Uint32Array }}
 */
export function refineQuads(positions, quads, project) {
  const pool = new NodePool(positions, project);
  const out = new Uint32Array(quads.length * 4);
  for (let f = 0, o = 0; f < quads.length; f += 4) {
    const v = [quads[f], quads[f + 1], quads[f + 2], quads[f + 3]];
    const m = [pool.mid(v[0], v[1]), pool.mid(v[1], v[2]), pool.mid(v[2], v[3]), pool.mid(v[3], v[0])];
    const c = pool.centre(v);
    for (let k = 0; k < 4; k++) { out[o++] = v[k]; out[o++] = m[k]; out[o++] = c; out[o++] = m[(k + 3) % 4]; }
  }
  return { positions: pool.positions(), quads: out };
}

// hexahedron corner k sits at lattice position HEX_AT[k] of the unit cube
const HEX_AT = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

function refineCell(pool, type, c, push) {
  const m = (i, j) => pool.mid(c[i], c[j]);
  if (type === 'tet4') {
    const m01 = m(0, 1), m02 = m(0, 2), m03 = m(0, 3), m12 = m(1, 2), m13 = m(1, 3), m23 = m(2, 3);
    push('tet4', [c[0], m01, m02, m03]); push('tet4', [m01, c[1], m12, m13]); push('tet4', [m02, m12, c[2], m23]); push('tet4', [m03, m13, m23, c[3]]);
    // inner octahedron cut along the m02–m13 diagonal (Bey's regular refinement)
    push('tet4', [m01, m02, m03, m13]); push('tet4', [m01, m02, m13, m12]); push('tet4', [m02, m03, m13, m23]); push('tet4', [m02, m12, m23, m13]);
  } else if (type === 'hex8') {
    const node = (i, j, k) => {
      const ids = [];
      for (let q = 0; q < 8; q++) {
        const a = HEX_AT[q];
        if ((i === 1 || a[0] === i / 2) && (j === 1 || a[1] === j / 2) && (k === 1 || a[2] === k / 2)) ids.push(c[q]);
      }
      return ids.length === 1 ? ids[0] : ids.length === 2 ? pool.mid(ids[0], ids[1]) : pool.centre(ids);
    };
    const L = [];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) L[(i * 3 + j) * 3 + k] = node(i, j, k);
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) push('hex8', HEX_AT.map((a) => L[((i + a[0]) * 3 + j + a[1]) * 3 + k + a[2]]));
  } else if (type === 'wedge6') {
    const b = [c[0], c[1], c[2]], t = [c[3], c[4], c[5]];
    const mid3 = [m(0, 3), m(1, 4), m(2, 5)];
    const eb = [m(0, 1), m(1, 2), m(2, 0)], et = [m(3, 4), m(4, 5), m(5, 3)];
    const em = [pool.centre([c[0], c[1], c[3], c[4]]), pool.centre([c[1], c[2], c[4], c[5]]), pool.centre([c[2], c[0], c[5], c[3]])];
    const layer = (lo, elo, hi, ehi) => {
      push('wedge6', [lo[0], elo[0], elo[2], hi[0], ehi[0], ehi[2]]);
      push('wedge6', [elo[0], lo[1], elo[1], ehi[0], hi[1], ehi[1]]);
      push('wedge6', [elo[2], elo[1], lo[2], ehi[2], ehi[1], hi[2]]);
      push('wedge6', [elo[0], elo[1], elo[2], ehi[0], ehi[1], ehi[2]]);
    };
    layer(b, eb, mid3, em);
    layer(mid3, em, t, et);
  } else if (type === 'pyr5') {
    const e = [m(0, 1), m(1, 2), m(2, 3), m(3, 0)], a = [m(0, 4), m(1, 4), m(2, 4), m(3, 4)];
    const cb = pool.centre([c[0], c[1], c[2], c[3]]);
    push('pyr5', [c[0], e[0], cb, e[3], a[0]]); push('pyr5', [e[0], c[1], e[1], cb, a[1]]);
    push('pyr5', [cb, e[1], c[2], e[2], a[2]]); push('pyr5', [e[3], cb, e[2], c[3], a[3]]);
    push('pyr5', [a[0], a[1], a[2], a[3], c[4]]);
    push('pyr5', [a[0], a[3], a[2], a[1], cb]);
    for (let k = 0; k < 4; k++) push('tet4', [e[k], a[k], a[(k + 1) % 4], cb]);
  } else throw new RangeError(`cells.type: '${type}' cannot be refined`);
}

/**
 * Uniform refinement of a volume mesh: every cell is split into eight (a pyramid into six pyramids and four
 * tetrahedra). Quadratic cells (tet10, hex20) are refined on their corner nodes and come out linear.
 * @param {object} model CONTRACT §4 volume model
 * @param {{ project?: Function }} [opts]
 * @returns {object} volume model; every block carries `parent` (Uint32Array: running index of the parent cell).
 *   Existing nodes keep their numbers, so node sets stay valid; `sets.cells` are mapped to the child cells.
 */
export function refineVolume(model, { project } = {}) {
  const pool = new NodePool(model.positions, project);
  const blocks = new Map();
  let parent = 0, total = 0;
  const push = (type, conn) => {
    let b = blocks.get(type);
    if (!b) { b = { conn: [], parent: [] }; blocks.set(type, b); }
    for (const v of conn) b.conn.push(v);
    b.parent.push(parent);
    total++;
  };
  for (const block of model.cells || []) {
    const info = CELL_INFO[block.type];
    if (!info) throw new RangeError(`cells.type: unsupported cell type '${block.type}'`);
    const lin = block.type === 'tet10' ? 'tet4' : block.type === 'hex20' ? 'hex8' : block.type;
    const nc = (block.conn.length / info.nodes) | 0;
    if (total + 8 * nc > LIMIT) throw new RangeError('model: refinement would exceed 20 million cells');
    const c = new Array(info.corners);
    for (let i = 0; i < nc; i++, parent++) {
      for (let k = 0; k < info.corners; k++) c[k] = block.conn[i * info.nodes + k];
      refineCell(pool, lin, c, push);
    }
  }
  const cells = [];
  const childOf = [];
  let base = 0;
  for (const [type, b] of blocks) {
    cells.push({ type, conn: Uint32Array.from(b.conn), parent: Uint32Array.from(b.parent) });
    b.parent.forEach((p, i) => { (childOf[p] || (childOf[p] = [])).push(base + i); });
    base += b.parent.length;
  }
  const out = { ...model, positions: pool.positions(), cells };
  delete out.nodeIds;
  if (model.sets && model.sets.cells && !ArrayBuffer.isView(model.sets.cells)) {
    const mapped = {};
    for (const [name, ids] of Object.entries(model.sets.cells)) mapped[name] = Uint32Array.from(Array.from(ids).flatMap((p) => childOf[p] || []));
    out.sets = { ...model.sets, cells: mapped };
  }
  return out;
}

/**
 * Divide every frame member into equal parts.
 * @param {object} frame CONTRACT §4 frame model
 * @param {{ parts?: number, maxLength?: number }} opts number of parts per member, or the longest allowed
 *   element [m] (each member is then cut into the smallest number of equal parts not exceeding it)
 * @returns {object} frame with `memberParent` (Uint32Array: index of the original member of every element);
 *   original nodes keep their numbers and section, material and group are inherited
 */
export function subdivideMembers(frame, { parts, maxLength } = {}) {
  if (!(parts >= 1) && !(maxLength > 0)) throw new RangeError('parts/maxLength: give parts >= 1 or maxLength > 0');
  const N = frame.nodes;
  const extra = [];
  const n0 = (N.length / 3) | 0;
  const members = [], parent = [];
  frame.members.forEach((m, mi) => {
    const dx = N[3 * m.b] - N[3 * m.a], dy = N[3 * m.b + 1] - N[3 * m.a + 1], dz = N[3 * m.b + 2] - N[3 * m.a + 2];
    const k = parts >= 1 ? Math.round(parts) : Math.max(1, Math.ceil(Math.hypot(dx, dy, dz) / maxLength - 1e-9));
    let prev = m.a;
    for (let q = 1; q <= k; q++) {
      let next = m.b;
      if (q < k) { next = n0 + extra.length / 3; extra.push(N[3 * m.a] + dx * q / k, N[3 * m.a + 1] + dy * q / k, N[3 * m.a + 2] + dz * q / k); }
      members.push({ ...m, a: prev, b: next });
      parent.push(mi);
      prev = next;
    }
  });
  const nodes = new Float64Array(N.length + extra.length);
  nodes.set(N); nodes.set(extra, N.length);
  const out = { ...frame, nodes, members, memberParent: Uint32Array.from(parent) };
  delete out.memberOrigin;
  return out;
}

/**
 * Resample a polyline at equal spacing along its length (the end points are kept; corners are cut where they
 * fall between samples). Numeric attributes with one value per point are interpolated linearly.
 * @param {object} polyline CONTRACT §4 polyline model
 * @param {{ spacing?: number, count?: number }} opts target spacing [m] or number of points (>= 2)
 * @returns {object} polyline with `attrs.kp` set to the chainage [m] when the input had none
 */
export function resamplePolyline(polyline, { spacing, count } = {}) {
  const P = polyline.points, n = (P.length / 3) | 0;
  if (n < 2) throw new RangeError('polyline.points: at least two points are required');
  const s = new Float64Array(n);
  for (let i = 1; i < n; i++) s[i] = s[i - 1] + Math.hypot(P[3 * i] - P[3 * i - 3], P[3 * i + 1] - P[3 * i - 2], P[3 * i + 2] - P[3 * i - 1]);
  const L = s[n - 1];
  let m;
  if (count >= 2) m = Math.round(count);
  else if (spacing > 0) m = Math.max(1, Math.round(L / spacing)) + 1;
  else throw new RangeError('spacing/count: give spacing > 0 or count >= 2');
  const pts = new Float64Array(3 * m);
  const seg = new Uint32Array(m), frac = new Float64Array(m);
  let j = 0;
  for (let i = 0; i < m; i++) {
    const t = i === m - 1 ? L : L * i / (m - 1);
    while (j < n - 2 && s[j + 1] < t) j++;
    const d = s[j + 1] - s[j];
    const f = d > 0 ? Math.min(1, Math.max(0, (t - s[j]) / d)) : 0;
    seg[i] = j; frac[i] = f;
    for (let k = 0; k < 3; k++) pts[3 * i + k] = P[3 * j + k] + f * (P[3 * j + 3 + k] - P[3 * j + k]);
  }
  const attrs = {};
  for (const [key, a] of Object.entries(polyline.attrs || {})) {
    if (a && a.length === n && typeof a[0] === 'number') {
      const o = new Float64Array(m);
      for (let i = 0; i < m; i++) o[i] = a[seg[i]] + frac[i] * (a[seg[i] + 1] - a[seg[i]]);
      attrs[key] = o;
    }
  }
  if (!attrs.kp) { attrs.kp = new Float64Array(m); for (let i = 0; i < m; i++) attrs.kp[i] = i === m - 1 ? L : L * i / (m - 1); }
  return { ...polyline, points: pts, attrs };
}

/**
 * Resample a grid to a new spacing over the same extent: bilinear interpolation when refining, the mean of the
 * valid source nodes inside each new cell when coarsening (so that sub-grid relief is averaged, not aliased).
 * @param {object} grid CONTRACT §4 grid model
 * @param {{ dx?: number, dy?: number, factor?: number }} opts new spacing [m] (dy defaults to dx with the sign
 *   of the original), or `factor` > 1 to refine / < 1 to coarsen both directions
 * @returns {object} grid model
 */
export function resampleGrid(grid, { dx, dy, factor } = {}) {
  let ndx, ndy;
  if (factor > 0) { ndx = grid.dx / factor; ndy = grid.dy / factor; }
  else if (dx > 0) { ndx = Math.sign(grid.dx) * dx; ndy = Math.sign(grid.dy) * (dy > 0 ? dy : dx); }
  else throw new RangeError('dx/factor: give dx > 0 or factor > 0');
  const ex = (grid.nx - 1) * grid.dx, ey = (grid.ny - 1) * grid.dy;
  const nx = Math.max(1, Math.floor(ex / ndx + 1e-9) + 1), ny = Math.max(1, Math.floor(ey / ndy + 1e-9) + 1);
  if (nx * ny > 5e7) throw new RangeError('dx/factor: the grid would have more than 50 million nodes');
  const z = new Float32Array(nx * ny);
  const rx = Math.abs(ndx / grid.dx), ry = Math.abs(ndy / grid.dy);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const x = grid.x0 + i * ndx, y = grid.y0 + j * ndy;
    if (rx <= 1 && ry <= 1) { z[j * nx + i] = gridSample(grid, x, y); continue; }
    const fi = i * rx, fj = j * ry;
    const i0 = Math.max(0, Math.ceil(fi - rx / 2 - 1e-9)), i1 = Math.min(grid.nx - 1, Math.floor(fi + rx / 2 + 1e-9));
    const j0 = Math.max(0, Math.ceil(fj - ry / 2 - 1e-9)), j1 = Math.min(grid.ny - 1, Math.floor(fj + ry / 2 + 1e-9));
    let sum = 0, cnt = 0;
    for (let b = j0; b <= j1; b++) for (let a = i0; a <= i1; a++) { const v = grid.z[b * grid.nx + a]; if (Number.isFinite(v)) { sum += v; cnt++; } }
    z[j * nx + i] = cnt ? sum / cnt : gridSample(grid, x, y);
  }
  return { ...grid, nx, ny, dx: ndx, dy: ndy, z };
}

/**
 * One level of uniform refinement of any model kind: surface 1→4, volume 1→8, frame members halved, polyline
 * and grid spacing halved.
 * @param {object} model CONTRACT §4 model
 * @param {{ project?: Function }} [opts] `project(x, y, z) -> [x, y, z]` snaps new nodes to the true geometry
 * @returns {object} refined model of the same kind
 */
export function refine(model, opts = {}) {
  if (!model) throw new RangeError('model: required');
  switch (model.kind) {
    case 'surface': {
      if (model.indices.length * 4 > 3 * LIMIT) throw new RangeError('model: refinement would exceed 20 million triangles');
      const r = refineTriangles(model.positions, model.indices, opts.project);
      const out = { ...model, positions: r.positions, indices: r.indices };
      delete out.faceCell;
      if (model.groups) out.groups = model.groups.map((g) => ({ ...g, start: g.start * 4, count: g.count * 4 }));
      return out;
    }
    case 'volume': return refineVolume(model, opts);
    case 'frame': return subdivideMembers(model, { parts: 2 });
    case 'polyline': return resamplePolyline(model, { count: 2 * ((model.points.length / 3) | 0) - 1 });
    case 'grid': return resampleGrid(model, { factor: 2 });
    default: throw new RangeError(`model.kind: a ${model.kind} model cannot be refined`);
  }
}

function edgeStats(positions, pairs) {
  let max = 0, sum = 0;
  const n = pairs.length / 2;
  for (let e = 0; e < n; e++) {
    const a = 3 * pairs[2 * e], b = 3 * pairs[2 * e + 1];
    const l = Math.hypot(positions[a] - positions[b], positions[a + 1] - positions[b + 1], positions[a + 2] - positions[b + 2]);
    sum += l; if (l > max) max = l;
  }
  return { max, mean: n ? sum / n : 0 };
}

/** Longest-edge bisection of a triangle surface until no edge is longer than `h`. */
function bisectSurface(model, h, project) {
  let P = Array.from(model.positions);
  let idx = Array.from(model.indices);
  const g0 = model.groups ? faceGroupIds(model) : null;
  let gid = g0 ? Array.from(g0.ids) : null;
  const len2 = (a, b) => (P[3 * a] - P[3 * b]) ** 2 + (P[3 * a + 1] - P[3 * b + 1]) ** 2 + (P[3 * a + 2] - P[3 * b + 2]) ** 2;
  const h2 = h * h * (1 + 1e-9);
  for (let pass = 0; pass < 60; pass++) {
    const F = idx.length / 3;
    const nV = P.length / 3;
    const key = (a, b) => (a < b ? a * nV + b : b * nV + a);
    const marked = new Map();
    for (let f = 0; f < F; f++) for (let k = 0; k < 3; k++) {
      const a = idx[3 * f + k], b = idx[3 * f + (k + 1) % 3];
      if (len2(a, b) > h2) marked.set(key(a, b), -1);
    }
    if (!marked.size) break;
    // closure: a triangle that is cut must also have its longest edge cut, which bounds the smallest angle
    for (let changed = true; changed;) {
      changed = false;
      for (let f = 0; f < F; f++) {
        let any = false, lk = 0, ll = -1;
        for (let k = 0; k < 3; k++) {
          const a = idx[3 * f + k], b = idx[3 * f + (k + 1) % 3];
          if (marked.has(key(a, b))) any = true;
          const l = len2(a, b);
          if (l > ll) { ll = l; lk = k; }
        }
        if (!any) continue;
        const kk = key(idx[3 * f + lk], idx[3 * f + (lk + 1) % 3]);
        if (!marked.has(kk)) { marked.set(kk, -1); changed = true; }
      }
    }
    if (F + 3 * marked.size > LIMIT) throw new RangeError('targetSize: refinement would exceed 20 million triangles');
    const mid = (a, b) => {
      const k = key(a, b);
      let m = marked.get(k);
      if (m === undefined) return -1;
      if (m < 0) {
        let x = (P[3 * a] + P[3 * b]) / 2, y = (P[3 * a + 1] + P[3 * b + 1]) / 2, z = (P[3 * a + 2] + P[3 * b + 2]) / 2;
        if (project) { const q = project(x, y, z); x = q[0]; y = q[1]; z = q[2]; }
        m = P.length / 3; P.push(x, y, z); marked.set(k, m);
      }
      return m;
    };
    const next = [], ng = gid ? [] : null;
    const emit = (f, a, b, c) => { next.push(a, b, c); if (ng) ng.push(gid[f]); };
    for (let f = 0; f < F; f++) {
      const v = [idx[3 * f], idx[3 * f + 1], idx[3 * f + 2]];
      const m = [mid(v[0], v[1]), mid(v[1], v[2]), mid(v[2], v[0])];
      const cnt = (m[0] >= 0) + (m[1] >= 0) + (m[2] >= 0);
      if (cnt === 0) emit(f, v[0], v[1], v[2]);
      else if (cnt === 3) { emit(f, v[0], m[0], m[2]); emit(f, m[0], v[1], m[1]); emit(f, m[2], m[1], v[2]); emit(f, m[0], m[1], m[2]); }
      else if (cnt === 1) {
        const k = m[0] >= 0 ? 0 : m[1] >= 0 ? 1 : 2;
        const a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3];
        emit(f, a, m[k], c); emit(f, m[k], b, c);
      } else {
        // rotate so that the uncut edge is c–a
        const k = m[0] < 0 ? 1 : m[1] < 0 ? 2 : 0;
        const a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3], mab = m[k], mbc = m[(k + 1) % 3];
        emit(f, mab, b, mbc);
        if (len2(a, mbc) <= len2(mab, c)) { emit(f, a, mab, mbc); emit(f, a, mbc, c); } else { emit(f, a, mab, c); emit(f, mab, mbc, c); }
      }
    }
    idx = next; gid = ng;
  }
  const out = { ...model, positions: new model.positions.constructor(P), indices: Uint32Array.from(idx) };
  delete out.faceCell;
  if (gid) { const r = groupsFromIds(out.indices, Int32Array.from(gid), g0.names); out.indices = r.indices; out.groups = r.groups; }
  return out;
}

/**
 * Refine a model until its elements are no larger than a target size.
 *  surface  — longest-edge bisection until every edge ≤ targetSize (conforming, keeps watertightness and groups)
 *  volume   — uniform refinement until the longest cell edge ≤ targetSize
 *  frame    — members cut into equal parts no longer than targetSize
 *  polyline — resampled at targetSize; grid — resampled to targetSize when that is finer than the present spacing
 * @param {object} model CONTRACT §4 model
 * @param {number} targetSize largest allowed element edge [m]
 * @param {{ project?: Function }} [opts] snaps new surface / volume nodes to the true geometry
 * @returns {object} refined model (the input itself when it already meets the target)
 */
export function refineTo(model, targetSize, opts = {}) {
  if (!model) throw new RangeError('model: required');
  if (!(targetSize > 0)) throw new RangeError('targetSize: must be > 0');
  switch (model.kind) {
    case 'surface': return bisectSurface(model, targetSize, opts.project);
    case 'volume': {
      let m = model;
      for (let guard = 0; guard < 12; guard++) {
        const pairs = [];
        for (const b of m.cells) {
          const info = CELL_INFO[b.type];
          const E = CELL_EDGES[b.type];
          for (let c = 0; c < b.conn.length; c += info.nodes) for (const [i, j] of E) pairs.push(b.conn[c + i], b.conn[c + j]);
        }
        if (edgeStats(m.positions, pairs).max <= targetSize * (1 + 1e-9)) return m;
        m = refineVolume(m, opts);
      }
      return m;
    }
    case 'frame': return subdivideMembers(model, { maxLength: targetSize });
    case 'polyline': return resamplePolyline(model, { spacing: targetSize });
    case 'grid': return Math.max(Math.abs(model.dx), Math.abs(model.dy)) <= targetSize ? model : resampleGrid(model, { dx: targetSize });
    default: throw new RangeError(`model.kind: a ${model.kind} model cannot be refined`);
  }
}

/* ------------------------------------------------------------------------------------------------ decimation */

class Heap {
  constructor() { this.a = []; }

  push(item) {
    const a = this.a;
    let i = a.length;
    a.push(item);
    while (i > 0) { const p = (i - 1) >> 1; if (a[p].cost <= item.cost) break; a[i] = a[p]; i = p; }
    a[i] = item;
  }

  pop() {
    const a = this.a;
    const top = a[0], last = a.pop();
    if (a.length) {
      let i = 0;
      for (;;) {
        let c = 2 * i + 1;
        if (c >= a.length) break;
        if (c + 1 < a.length && a[c + 1].cost < a[c].cost) c++;
        if (a[c].cost >= last.cost) break;
        a[i] = a[c]; i = c;
      }
      a[i] = last;
    }
    return top;
  }

  get size() { return this.a.length; }
}

/**
 * Quadric-error-metric decimation of a triangle surface (Garland & Heckbert 1997): edges are collapsed in order
 * of increasing error against the planes of the original faces. Collapses that would fold a face over, pinch
 * the surface (link condition) or move a boundary are refused or penalised, so a watertight manifold input
 * stays watertight and manifold, and open boundaries and group borders keep their outline.
 *
 * Vertices are expected to be shared between faces (weld STL-style input first).
 *
 * @param {object} model CONTRACT §4 surface model
 * @param {number} targetCount number of triangles wanted (the result may stop above it when no further collapse
 *   is admissible)
 * @param {{ borderWeight?: number, maxNormalTurnDeg?: number }} [opts] penalty factor for moving boundary and
 *   group-border edges (default 100); largest rotation of a face normal allowed in one collapse (default 80°)
 * @returns {object} surface model with groups preserved and `decimation: { from, to, maxError }` (largest
 *   accepted quadric error, m² × m² area weighting)
 */
export function decimate(model, targetCount, { borderWeight = 100, maxNormalTurnDeg = 80 } = {}) {
  if (!model || model.kind !== 'surface') throw new RangeError('model: a surface model is required');
  if (!(targetCount >= 0)) throw new RangeError('targetCount: must be >= 0');
  const nV = (model.positions.length / 3) | 0;
  const F0 = (model.indices.length / 3) | 0;
  if (F0 <= targetCount) return model;
  const P = Float64Array.from(model.positions);
  const T = Int32Array.from(model.indices);
  const alive = new Uint8Array(F0).fill(1);
  const gid = faceGroupIds(model).ids;
  const vf = Array.from({ length: nV }, () => []);
  for (let f = 0; f < F0; f++) for (let k = 0; k < 3; k++) vf[T[3 * f + k]].push(f);
  const Q = new Float64Array(10 * nV);
  const addPlane = (v, a, b, c, d, w) => {
    const q = 10 * v;
    Q[q] += w * a * a; Q[q + 1] += w * a * b; Q[q + 2] += w * a * c; Q[q + 3] += w * a * d;
    Q[q + 4] += w * b * b; Q[q + 5] += w * b * c; Q[q + 6] += w * b * d; Q[q + 7] += w * c * c; Q[q + 8] += w * c * d; Q[q + 9] += w * d * d;
  };
  const normal = (f, out) => {
    const a = 3 * T[3 * f], b = 3 * T[3 * f + 1], c = 3 * T[3 * f + 2];
    const ux = P[b] - P[a], uy = P[b + 1] - P[a + 1], uz = P[b + 2] - P[a + 2], vx = P[c] - P[a], vy = P[c + 1] - P[a + 1], vz = P[c + 2] - P[a + 2];
    out[0] = uy * vz - uz * vy; out[1] = uz * vx - ux * vz; out[2] = ux * vy - uy * vx;
    return out;
  };
  const nrm = [0, 0, 0];
  for (let f = 0; f < F0; f++) {
    normal(f, nrm);
    const l = Math.hypot(nrm[0], nrm[1], nrm[2]);
    if (!(l > 0)) continue;
    const a = nrm[0] / l, b = nrm[1] / l, c = nrm[2] / l, v0 = 3 * T[3 * f];
    const d = -(a * P[v0] + b * P[v0 + 1] + c * P[v0 + 2]);
    for (let k = 0; k < 3; k++) addPlane(T[3 * f + k], a, b, c, d, l / 2);
  }
  // boundary edges and edges between groups: a plane through the edge, square to the face, holds the outline
  const ed = buildEdges(model.indices, nV);
  for (let e = 0; e < ed.edgeCount; e++) {
    const s = ed.edgeFaceStart[e], n = ed.edgeFaceStart[e + 1] - s;
    const f0 = (ed.edgeHalfEdges[s] / 3) | 0;
    const border = n === 1 || (n === 2 && gid[f0] !== gid[(ed.edgeHalfEdges[s + 1] / 3) | 0]) || n > 2;
    if (!border) continue;
    const a = ed.edgeVerts[2 * e], b = ed.edgeVerts[2 * e + 1];
    for (let q = s; q < s + n; q++) {
      normal((ed.edgeHalfEdges[q] / 3) | 0, nrm);
      const ex = P[3 * b] - P[3 * a], ey = P[3 * b + 1] - P[3 * a + 1], ez = P[3 * b + 2] - P[3 * a + 2];
      let px = ey * nrm[2] - ez * nrm[1], py = ez * nrm[0] - ex * nrm[2], pz = ex * nrm[1] - ey * nrm[0];
      const l = Math.hypot(px, py, pz);
      if (!(l > 0)) continue;
      px /= l; py /= l; pz /= l;
      const d = -(px * P[3 * a] + py * P[3 * a + 1] + pz * P[3 * a + 2]);
      const w = borderWeight * (ex * ex + ey * ey + ez * ez);
      addPlane(a, px, py, pz, d, w); addPlane(b, px, py, pz, d, w);
    }
  }
  const stamp = new Uint32Array(nV);
  const heap = new Heap();
  let areaSum = 0;
  for (let f = 0; f < F0; f++) { normal(f, nrm); areaSum += Math.hypot(nrm[0], nrm[1], nrm[2]) / 2; }
  const tie = 1e-7 * areaSum / F0;
  const q = new Float64Array(10);
  const evalQ = (x, y, z) => q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y + 2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
  const consider = (a, b) => {
    for (let k = 0; k < 10; k++) q[k] = Q[10 * a + k] + Q[10 * b + k];
    const det = q[0] * (q[4] * q[7] - q[5] * q[5]) - q[1] * (q[1] * q[7] - q[5] * q[2]) + q[2] * (q[1] * q[5] - q[4] * q[2]);
    const scale = Math.abs(q[0]) + Math.abs(q[4]) + Math.abs(q[7]);
    let x, y, z, cost = Infinity;
    if (Math.abs(det) > 1e-9 * scale * scale * scale) {
      x = -(q[3] * (q[4] * q[7] - q[5] * q[5]) - q[1] * (q[6] * q[7] - q[5] * q[8]) + q[2] * (q[6] * q[5] - q[4] * q[8])) / det;
      y = -(q[0] * (q[6] * q[7] - q[8] * q[5]) - q[3] * (q[1] * q[7] - q[5] * q[2]) + q[2] * (q[1] * q[8] - q[6] * q[2])) / det;
      z = -(q[0] * (q[4] * q[8] - q[5] * q[6]) - q[1] * (q[1] * q[8] - q[6] * q[2]) + q[3] * (q[1] * q[5] - q[4] * q[2])) / det;
      // an optimum far from the edge is the mark of a nearly singular quadric: fall back to the edge itself
      const mx = (P[3 * a] + P[3 * b]) / 2, my = (P[3 * a + 1] + P[3 * b + 1]) / 2, mz = (P[3 * a + 2] + P[3 * b + 2]) / 2;
      const el = Math.hypot(P[3 * a] - P[3 * b], P[3 * a + 1] - P[3 * b + 1], P[3 * a + 2] - P[3 * b + 2]);
      if (Math.hypot(x - mx, y - my, z - mz) <= 2 * el) cost = evalQ(x, y, z);
    }
    if (!Number.isFinite(cost)) {
      for (const t of [0.5, 0, 1]) {
        const cx = P[3 * a] + t * (P[3 * b] - P[3 * a]), cy = P[3 * a + 1] + t * (P[3 * b + 1] - P[3 * a + 1]), cz = P[3 * a + 2] + t * (P[3 * b + 2] - P[3 * a + 2]);
        const c = evalQ(cx, cy, cz);
        if (c < cost) { cost = c; x = cx; y = cy; z = cz; }
      }
    }
    // a small length term orders the collapses on flat regions (where every error is zero) shortest edge first
    const e2 = (P[3 * a] - P[3 * b]) ** 2 + (P[3 * a + 1] - P[3 * b + 1]) ** 2 + (P[3 * a + 2] - P[3 * b + 2]) ** 2;
    heap.push({ cost: Math.max(cost, 0) + tie * e2, err: Math.max(cost, 0), a, b, sa: stamp[a], sb: stamp[b], x, y, z });
  };
  for (let e = 0; e < ed.edgeCount; e++) consider(ed.edgeVerts[2 * e], ed.edgeVerts[2 * e + 1]);
  const cosMin = Math.cos(maxNormalTurnDeg * Math.PI / 180);
  const n0 = [0, 0, 0], n1 = [0, 0, 0];
  const ring = (v, set) => { for (const f of vf[v]) for (let k = 0; k < 3; k++) { const u = T[3 * f + k]; if (u !== v) set.add(u); } return set; };
  let faces = F0, maxError = 0;
  while (faces > targetCount && heap.size) {
    const c = heap.pop();
    const { a, b } = c;
    if (stamp[a] !== c.sa || stamp[b] !== c.sb) continue;
    // link condition: the two ends may share only the vertices opposite the edge
    const shared = vf[a].filter((f) => T[3 * f] === b || T[3 * f + 1] === b || T[3 * f + 2] === b);
    if (!shared.length || shared.length > 2) continue;
    const ra = ring(a, new Set()), rb = ring(b, new Set());
    let common = 0;
    for (const u of ra) if (rb.has(u)) common++;
    if (common !== shared.length) continue;
    // a closed fan must keep at least a tetrahedron's worth of faces
    if (shared.length === 2 && ra.size + rb.size - common - 2 < 3) continue;
    let ok = true;
    for (const v of [a, b]) {
      for (const f of vf[v]) {
        if (shared.includes(f)) continue;
        normal(f, n0);
        const keep = [P[3 * v], P[3 * v + 1], P[3 * v + 2]];
        P[3 * v] = c.x; P[3 * v + 1] = c.y; P[3 * v + 2] = c.z;
        normal(f, n1);
        P[3 * v] = keep[0]; P[3 * v + 1] = keep[1]; P[3 * v + 2] = keep[2];
        const l0 = Math.hypot(n0[0], n0[1], n0[2]), l1 = Math.hypot(n1[0], n1[1], n1[2]);
        if (!(l1 > 1e-12 * l0) || (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]) < cosMin * l0 * l1) { ok = false; break; }
      }
      if (!ok) break;
    }
    if (!ok) continue;
    for (const f of shared) {
      alive[f] = 0; faces--;
      for (let k = 0; k < 3; k++) { const u = T[3 * f + k]; if (u !== a && u !== b) vf[u] = vf[u].filter((g) => g !== f); }
    }
    const fa = vf[a].filter((f) => alive[f]);
    for (const f of vf[b]) {
      if (!alive[f]) continue;
      for (let k = 0; k < 3; k++) if (T[3 * f + k] === b) T[3 * f + k] = a;
      fa.push(f);
    }
    vf[a] = fa; vf[b] = [];
    for (let k = 0; k < 10; k++) Q[10 * a + k] += Q[10 * b + k];
    P[3 * a] = c.x; P[3 * a + 1] = c.y; P[3 * a + 2] = c.z;
    stamp[a]++; stamp[b]++;
    if (c.err > maxError) maxError = c.err;
    for (const u of ring(a, new Set())) consider(a, u);
  }
  const map = new Int32Array(nV).fill(-1);
  const pos = [], idx = [], ids = [];
  for (let f = 0; f < F0; f++) {
    if (!alive[f]) continue;
    for (let k = 0; k < 3; k++) {
      const v = T[3 * f + k];
      if (map[v] < 0) { map[v] = pos.length / 3; pos.push(P[3 * v], P[3 * v + 1], P[3 * v + 2]); }
      idx.push(map[v]);
    }
    ids.push(gid[f]);
  }
  const out = { ...model, positions: new model.positions.constructor(pos), indices: Uint32Array.from(idx), decimation: { from: F0, to: idx.length / 3, maxError } };
  delete out.faceCell;
  if (model.groups) { const r = groupsFromIds(out.indices, Int32Array.from(ids), model.groups.map((g) => g.name)); out.indices = r.indices; out.groups = r.groups; }
  return out;
}

/**
 * Coarsen a model to about a target element count.
 *  surface  — quadric-error decimation to `targetCount` triangles (see `decimate`)
 *  points   — every k-th point so that about `targetCount` remain
 *  polyline — resampled to `targetCount` points; grid — resampled so that about `targetCount` nodes remain
 * Volume and frame meshes cannot be coarsened here (re-mesh the geometry instead).
 * @param {object} model CONTRACT §4 model
 * @param {number} targetCount wanted number of triangles / points / nodes
 * @param {object} [opts] passed to `decimate` for surfaces
 * @returns {object} coarsened model (the input itself when it is already at or below the target)
 */
export function coarsenTo(model, targetCount, opts = {}) {
  if (!model) throw new RangeError('model: required');
  if (!(targetCount >= 1)) throw new RangeError('targetCount: must be >= 1');
  switch (model.kind) {
    case 'surface': return decimate(model, targetCount, opts);
    case 'points': {
      const n = (model.positions.length / 3) | 0;
      if (n <= targetCount) return model;
      const m = Math.round(targetCount);
      const keep = new Uint32Array(m);
      for (let i = 0; i < m; i++) keep[i] = Math.min(n - 1, Math.floor(i * n / m));
      const pos = new model.positions.constructor(3 * m);
      for (let i = 0; i < m; i++) { pos[3 * i] = model.positions[3 * keep[i]]; pos[3 * i + 1] = model.positions[3 * keep[i] + 1]; pos[3 * i + 2] = model.positions[3 * keep[i] + 2]; }
      const out = { ...model, positions: pos };
      if (model.attrs) { out.attrs = {}; for (const [k, a] of Object.entries(model.attrs)) out.attrs[k] = a && a.length === n ? a.constructor.from(keep, (i) => a[i]) : a; }
      return out;
    }
    case 'polyline': return (model.points.length / 3) <= targetCount ? model : resamplePolyline(model, { count: Math.max(2, Math.round(targetCount)) });
    case 'grid': {
      const n = model.nx * model.ny;
      if (n <= targetCount) return model;
      return resampleGrid(model, { factor: Math.sqrt(targetCount / n) });
    }
    default: throw new RangeError(`model.kind: a ${model.kind} model cannot be coarsened; re-mesh it at a larger element size`);
  }
}

/**
 * Laplacian smoothing of a triangle surface: every free vertex moves towards the mean of its edge neighbours.
 * With `taubin` (default) each step λ is followed by a negative step μ, giving the transfer function
 * (1 − λk)(1 − μk) per pass. The default μ = −λ makes it 1 − λ²k²: high-frequency noise is damped while smooth
 * shape (small k) is left almost untouched and nothing is amplified, so the body does not shrink as it does
 * under plain Laplacian smoothing. Boundary vertices and vertices on group borders stay where they are.
 * @param {object} model CONTRACT §4 surface model (shared vertices)
 * @param {{ iterations?: number, lambda?: number, taubin?: boolean, mu?: number, fixBoundary?: boolean,
 *   fixGroupBorders?: boolean, fixed?: ArrayLike<number> }} [opts] iterations (default 10); step λ in (0, 1]
 *   (default 0.5); μ ≤ −λ (default −λ; Taubin's pass-band variant is μ = −λ / (1 − k_PB·λ) with k_PB ≈ 0.1,
 *   which slightly amplifies wavelengths longer than the pass band); further vertex ids to hold in `fixed`
 * @returns {object} surface model with new positions (connectivity unchanged)
 */
export function smooth(model, { iterations = 10, lambda = 0.5, taubin = true, mu, fixBoundary = true, fixGroupBorders = true, fixed } = {}) {
  if (!model || model.kind !== 'surface') throw new RangeError('model: a surface model is required');
  if (!(lambda > 0 && lambda <= 1)) throw new RangeError('lambda: must be within (0, 1]');
  if (mu === undefined) mu = -lambda;
  if (taubin && !(mu <= -lambda)) throw new RangeError('mu: must not be above −lambda');
  const nV = (model.positions.length / 3) | 0;
  const ed = buildEdges(model.indices, nV);
  const hold = new Uint8Array(nV);
  if (fixed) for (const v of fixed) hold[v] = 1;
  const gid = faceGroupIds(model).ids;
  for (let e = 0; e < ed.edgeCount; e++) {
    const s = ed.edgeFaceStart[e], n = ed.edgeFaceStart[e + 1] - s;
    const border = (fixBoundary && n !== 2) || (fixGroupBorders && n === 2 && gid[(ed.edgeHalfEdges[s] / 3) | 0] !== gid[(ed.edgeHalfEdges[s + 1] / 3) | 0]);
    if (border) { hold[ed.edgeVerts[2 * e]] = 1; hold[ed.edgeVerts[2 * e + 1]] = 1; }
  }
  const deg = new Uint32Array(nV);
  for (let e = 0; e < ed.edgeCount; e++) { deg[ed.edgeVerts[2 * e]]++; deg[ed.edgeVerts[2 * e + 1]]++; }
  let P = Float64Array.from(model.positions);
  const acc = new Float64Array(3 * nV);
  const step = (f) => {
    acc.fill(0);
    for (let e = 0; e < ed.edgeCount; e++) {
      const a = 3 * ed.edgeVerts[2 * e], b = 3 * ed.edgeVerts[2 * e + 1];
      acc[a] += P[b]; acc[a + 1] += P[b + 1]; acc[a + 2] += P[b + 2];
      acc[b] += P[a]; acc[b + 1] += P[a + 1]; acc[b + 2] += P[a + 2];
    }
    const N = new Float64Array(3 * nV);
    for (let v = 0; v < nV; v++) {
      for (let k = 0; k < 3; k++) N[3 * v + k] = hold[v] || !deg[v] ? P[3 * v + k] : P[3 * v + k] + f * (acc[3 * v + k] / deg[v] - P[3 * v + k]);
    }
    P = N;
  };
  for (let it = 0; it < iterations; it++) { step(lambda); if (taubin) step(mu); }
  return { ...model, positions: new model.positions.constructor(P) };
}
