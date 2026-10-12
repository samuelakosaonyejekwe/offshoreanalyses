/**
 * Mesh topology: vertex welding, edge–face adjacency, boundary loops, manifoldness, connected components,
 * orientation, Euler characteristic and a one-call inspection for every geometry model kind (CONTRACT §4).
 *
 * Conventions
 *  - Triangle surfaces: `positions` (x, y, z per vertex) [m], `indices` (three vertex ids per face).
 *    Half-edge h = 3·f + k runs from corner k to corner (k + 1) mod 3 of face f.
 *  - Outward orientation: counter-clockwise seen from outside, so the divergence-theorem volume is positive.
 *  - Volume cells use the Gmsh / Abaqus / CGNS corner ordering (see CELL_FACES); mid-side nodes of tet10 and
 *    hex20 are ignored for topology (corner nodes only).
 *  - Grid: z[iy·nx + ix] at x = x0 + ix·dx, y = y0 + iy·dy; NaN marks no data.
 *
 * References: Edelsbrunner, "Geometry and Topology for Mesh Generation" (2001) ch. 3 (Euler–Poincaré formula
 * χ = V − E + F = 2C − 2g − b); Botsch et al., "Polygon Mesh Processing" (2010) ch. 1–2 (manifold conditions).
 */

/** Number of corner nodes and total nodes per cell type. */
export const CELL_INFO = {
  tet4: { corners: 4, nodes: 4 }, tet10: { corners: 4, nodes: 10 },
  hex8: { corners: 8, nodes: 8 }, hex20: { corners: 8, nodes: 20 },
  wedge6: { corners: 6, nodes: 6 }, pyr5: { corners: 5, nodes: 5 },
};

/** Outward-oriented faces by local corner index. */
export const CELL_FACES = {
  tet4: [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
  hex8: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
  wedge6: [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
  pyr5: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
};
CELL_FACES.tet10 = CELL_FACES.tet4;
CELL_FACES.hex20 = CELL_FACES.hex8;

/** Edges by local corner index. */
export const CELL_EDGES = {
  tet4: [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]],
  hex8: [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]],
  wedge6: [[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3], [0, 3], [1, 4], [2, 5]],
  pyr5: [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 4], [2, 4], [3, 4]],
};
CELL_EDGES.tet10 = CELL_EDGES.tet4;
CELL_EDGES.hex20 = CELL_EDGES.hex8;

const NONE = 0xffffffff;

function requireSurface(model, name = 'model') {
  if (!model || !model.positions || !model.indices) throw new RangeError(`${name}: a triangle surface with positions and indices is required`);
  if (model.indices.length % 3 !== 0) throw new RangeError(`${name}.indices: length must be a multiple of 3`);
}

/** Axis-aligned bounds of a packed coordinate array. Returns { min: [x,y,z], max: [x,y,z], diagonal } [m]. */
export function bounds(positions, dim = 3) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += dim) {
    for (let k = 0; k < dim; k++) {
      const v = positions[i + k];
      if (v < min[k]) min[k] = v;
      if (v > max[k]) max[k] = v;
    }
  }
  if (dim < 3) { min[2] = 0; max[2] = 0; }
  if (!(min[0] <= max[0])) return { min: [0, 0, 0], max: [0, 0, 0], diagonal: 0 };
  return { min, max, diagonal: Math.hypot(max[0] - min[0], max[1] - min[1], max[2] - min[2]) };
}

/**
 * Merge vertices closer than `tol` using a spatial hash with cell size `tol` (27-cell neighbourhood search).
 * Greedy: each vertex joins the first earlier vertex within the tolerance, so chains longer than `tol` are not
 * collapsed transitively. tol = 0 merges bit-identical coordinates only.
 * @param {ArrayLike<number>} positions packed x, y, z [m]
 * @param {number} tol merge distance [m], >= 0
 * @returns {{ positions: Float32Array|Float64Array, map: Uint32Array, count: number, merged: number, rep: Uint32Array }}
 *   new packed positions (same array type), old→new vertex map, new vertex count, number of vertices removed,
 *   and for every new vertex the first old vertex it represents.
 */
export function weldVertices(positions, tol = 0) {
  if (!(tol >= 0)) throw new RangeError('tol: must be >= 0');
  const n = (positions.length / 3) | 0;
  const map = new Uint32Array(n);
  const rep = new Uint32Array(n);
  const head = new Map();
  const next = new Int32Array(n).fill(-1);
  const b = bounds(positions);
  const cell = tol > 0 ? tol : 1;
  const inv = 1 / cell;
  const ox = b.min[0], oy = b.min[1], oz = b.min[2];
  const tol2 = tol * tol;
  const key = (ix, iy, iz) => (Math.imul(ix, 73856093) ^ Math.imul(iy, 19349663) ^ Math.imul(iz, 83492791)) | 0;
  let count = 0;
  for (let i = 0; i < n; i++) {
    const x = positions[3 * i], y = positions[3 * i + 1], z = positions[3 * i + 2];
    let found = -1;
    let k0;
    if (tol > 0) {
      const ix = Math.floor((x - ox) * inv), iy = Math.floor((y - oy) * inv), iz = Math.floor((z - oz) * inv);
      k0 = key(ix, iy, iz);
      search:
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        let j = head.get(key(ix + dx, iy + dy, iz + dz));
        while (j !== undefined && j >= 0) {
          const r = rep[j];
          const ex = positions[3 * r] - x, ey = positions[3 * r + 1] - y, ez = positions[3 * r + 2] - z;
          if (ex * ex + ey * ey + ez * ez <= tol2) { found = j; break search; }
          j = next[j];
        }
      }
    } else {
      // exact match: hash the raw coordinate values (−0 and +0 treated alike)
      k0 = key((x * 1e6) | 0, (y * 1e6) | 0, (z * 1e6) | 0) ^ Math.imul(Math.fround(x + y * 3 + z * 7) * 1e3 | 0, 2654435761);
      let j = head.get(k0);
      while (j !== undefined && j >= 0) {
        const r = rep[j];
        if (positions[3 * r] === x && positions[3 * r + 1] === y && positions[3 * r + 2] === z) { found = j; break; }
        j = next[j];
      }
    }
    if (found >= 0) { map[i] = found; continue; }
    const h = head.get(k0);
    next[count] = h === undefined ? -1 : h;
    head.set(k0, count);
    rep[count] = i;
    map[i] = count++;
  }
  const out = new positions.constructor(3 * count);
  for (let j = 0; j < count; j++) {
    const r = rep[j];
    out[3 * j] = positions[3 * r]; out[3 * j + 1] = positions[3 * r + 1]; out[3 * j + 2] = positions[3 * r + 2];
  }
  return { positions: out, map, count, merged: n - count, rep: rep.slice(0, count) };
}

/**
 * Weld a surface model: merges vertices within `tol` and remaps the faces (faces are not removed).
 * @param {object} model surface model
 * @param {number} [tol] merge distance [m]; default 1e-6 × bounding-box diagonal
 * @returns {{ model: object, merged: number, map: Uint32Array }}
 */
export function weld(model, tol) {
  requireSurface(model);
  const t = tol ?? 1e-6 * bounds(model.positions).diagonal;
  const w = weldVertices(model.positions, t);
  const indices = new Uint32Array(model.indices.length);
  for (let i = 0; i < indices.length; i++) indices[i] = w.map[model.indices[i]];
  return { model: { ...model, positions: w.positions, indices }, merged: w.merged, map: w.map };
}

/**
 * Edge table of a triangle mesh.
 * @param {ArrayLike<number>} indices three vertex ids per face
 * @param {number} nVerts vertex count (ids must be < nVerts)
 * @returns {{ edgeCount: number, edgeVerts: Uint32Array, faceEdges: Uint32Array, edgeFaceStart: Uint32Array,
 *   edgeHalfEdges: Uint32Array, edgeFaceCount: (e: number) => number }}
 *   edgeVerts[2e], [2e+1] = end vertices (lower id first); faceEdges[3f+k] = edge id of half-edge k (0xFFFFFFFF
 *   for a collapsed edge with both ends on the same vertex); the half-edges of edge e are
 *   edgeHalfEdges[edgeFaceStart[e] .. edgeFaceStart[e+1]) (face = h / 3 | 0, corner = h % 3).
 */
export function buildEdges(indices, nVerts) {
  const F = (indices.length / 3) | 0;
  const map = new Map();
  const faceEdges = new Uint32Array(3 * F);
  const ev = [];
  const cnt = [];
  let E = 0;
  for (let f = 0; f < F; f++) {
    for (let k = 0; k < 3; k++) {
      const a = indices[3 * f + k], b = indices[3 * f + (k + 1) % 3];
      if (a === b) { faceEdges[3 * f + k] = NONE; continue; }
      const key = a < b ? a * nVerts + b : b * nVerts + a;
      let id = map.get(key);
      if (id === undefined) { id = E++; map.set(key, id); ev.push(a < b ? a : b, a < b ? b : a); cnt.push(0); }
      cnt[id]++;
      faceEdges[3 * f + k] = id;
    }
  }
  const edgeFaceStart = new Uint32Array(E + 1);
  for (let e = 0; e < E; e++) edgeFaceStart[e + 1] = edgeFaceStart[e] + cnt[e];
  const fill = edgeFaceStart.slice(0, E);
  const edgeHalfEdges = new Uint32Array(edgeFaceStart[E]);
  for (let h = 0; h < 3 * F; h++) {
    const e = faceEdges[h];
    if (e !== NONE) edgeHalfEdges[fill[e]++] = h;
  }
  return {
    edgeCount: E, edgeVerts: Uint32Array.from(ev), faceEdges, edgeFaceStart, edgeHalfEdges,
    edgeFaceCount: (e) => edgeFaceStart[e + 1] - edgeFaceStart[e],
  };
}

/**
 * Boundary loops: chains of edges used by exactly one face, in the direction of that face's half-edge.
 * @returns {{ loops: number[][], closed: boolean[], edgeCount: number }} vertex ids per loop; `closed[i]` is
 *   false where a chain could not be closed (non-manifold boundary).
 */
export function boundaryLoops(indices, nVerts, edges = buildEdges(indices, nVerts)) {
  const out = new Map();
  let nb = 0;
  for (let e = 0; e < edges.edgeCount; e++) {
    if (edges.edgeFaceStart[e + 1] - edges.edgeFaceStart[e] !== 1) continue;
    const h = edges.edgeHalfEdges[edges.edgeFaceStart[e]];
    const f = (h / 3) | 0, k = h % 3;
    const a = indices[3 * f + k], b = indices[3 * f + (k + 1) % 3];
    const list = out.get(a);
    if (list) list.push(b); else out.set(a, [b]);
    nb++;
  }
  const loops = [], closed = [];
  for (const [start] of out) {
    for (;;) {
      const first = out.get(start);
      if (!first || !first.length) break;
      const loop = [start];
      let v = first.pop();
      let ok = false;
      for (let guard = 0; guard <= nb; guard++) {
        if (v === start) { ok = true; break; }
        loop.push(v);
        const nx = out.get(v);
        if (!nx || !nx.length) break;
        v = nx.pop();
      }
      loops.push(loop); closed.push(ok);
    }
  }
  return { loops, closed, edgeCount: nb };
}

/** Edges shared by more than two faces. Returns [[i, j], ...] vertex pairs. */
export function nonManifoldEdges(indices, nVerts, edges = buildEdges(indices, nVerts)) {
  const res = [];
  for (let e = 0; e < edges.edgeCount; e++) {
    if (edges.edgeFaceStart[e + 1] - edges.edgeFaceStart[e] > 2) res.push([edges.edgeVerts[2 * e], edges.edgeVerts[2 * e + 1]]);
  }
  return res;
}

function ufFind(p, i) { while (p[i] !== i) { p[i] = p[p[i]]; i = p[i]; } return i; }

/**
 * Corner fans: for every face corner the id of the edge-connected fan of faces it belongs to around its vertex.
 * A manifold vertex has one fan. Internal helper shared with the repair tools.
 */
export function cornerFans(indices, nVerts, edges = buildEdges(indices, nVerts)) {
  const n = indices.length;
  const p = new Uint32Array(n);
  for (let i = 0; i < n; i++) p[i] = i;
  for (let e = 0; e < edges.edgeCount; e++) {
    const s = edges.edgeFaceStart[e], t = edges.edgeFaceStart[e + 1];
    if (t - s < 2) continue;
    const u = edges.edgeVerts[2 * e], v = edges.edgeVerts[2 * e + 1];
    let cu = -1, cv = -1;
    for (let q = s; q < t; q++) {
      const f = (edges.edgeHalfEdges[q] / 3) | 0;
      for (let k = 0; k < 3; k++) {
        const c = 3 * f + k;
        if (indices[c] === u) { if (cu < 0) cu = c; else p[ufFind(p, c)] = ufFind(p, cu); }
        else if (indices[c] === v) { if (cv < 0) cv = c; else p[ufFind(p, c)] = ufFind(p, cv); }
      }
    }
  }
  for (let i = 0; i < n; i++) p[i] = ufFind(p, i);
  return p;
}

/** Vertices whose incident faces form more than one edge-connected fan (bow-tie vertices). Returns vertex ids. */
export function nonManifoldVertices(indices, nVerts, edges = buildEdges(indices, nVerts)) {
  const fans = cornerFans(indices, nVerts, edges);
  const first = new Int32Array(nVerts).fill(-1);
  const flag = new Uint8Array(nVerts);
  for (let c = 0; c < indices.length; c++) {
    const v = indices[c];
    if (first[v] < 0) first[v] = fans[c];
    else if (first[v] !== fans[c]) flag[v] = 1;
  }
  const res = [];
  for (let v = 0; v < nVerts; v++) if (flag[v]) res.push(v);
  return res;
}

/**
 * Connected components of the faces.
 * @param {'edge'|'vertex'} [by] faces are connected through shared edges (default) or shared vertices
 * @returns {{ count: number, faceComponent: Int32Array, sizes: number[] }}
 */
export function connectedComponents(indices, nVerts, by = 'edge', edges) {
  const F = (indices.length / 3) | 0;
  const p = new Uint32Array(F);
  for (let i = 0; i < F; i++) p[i] = i;
  if (by === 'vertex') {
    const first = new Int32Array(nVerts).fill(-1);
    for (let c = 0; c < indices.length; c++) {
      const v = indices[c], f = (c / 3) | 0;
      if (first[v] < 0) first[v] = f; else p[ufFind(p, f)] = ufFind(p, first[v]);
    }
  } else {
    const ed = edges || buildEdges(indices, nVerts);
    for (let e = 0; e < ed.edgeCount; e++) {
      const s = ed.edgeFaceStart[e], t = ed.edgeFaceStart[e + 1];
      const f0 = (ed.edgeHalfEdges[s] / 3) | 0;
      for (let q = s + 1; q < t; q++) p[ufFind(p, (ed.edgeHalfEdges[q] / 3) | 0)] = ufFind(p, f0);
    }
  }
  const faceComponent = new Int32Array(F);
  const ids = new Map();
  const sizes = [];
  for (let f = 0; f < F; f++) {
    const r = ufFind(p, f);
    let id = ids.get(r);
    if (id === undefined) { id = sizes.length; ids.set(r, id); sizes.push(0); }
    faceComponent[f] = id; sizes[id]++;
  }
  return { count: sizes.length, faceComponent, sizes };
}

/** Faces using the same three vertices as an earlier face (any orientation). Returns face ids. */
export function duplicateFaces(indices, nVerts) {
  const F = (indices.length / 3) | 0;
  const seen = new Set();
  const res = [];
  for (let f = 0; f < F; f++) {
    let a = indices[3 * f], b = indices[3 * f + 1], c = indices[3 * f + 2], t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    const key = `${a},${b},${c}`;
    if (seen.has(key)) res.push(f); else seen.add(key);
  }
  return res;
}

/** Twice the triangle area. */
function area2(P, a, b, c) {
  const ux = P[3 * b] - P[3 * a], uy = P[3 * b + 1] - P[3 * a + 1], uz = P[3 * b + 2] - P[3 * a + 2];
  const vx = P[3 * c] - P[3 * a], vy = P[3 * c + 1] - P[3 * a + 1], vz = P[3 * c + 2] - P[3 * a + 2];
  return Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
}

/**
 * Faces with a repeated vertex or with a height below `relTol` × longest edge (collinear corners).
 * @param {number} [relTol] slenderness limit, default 1e-7
 * @returns {number[]} face ids
 */
export function degenerateFaces(positions, indices, relTol = 1e-7) {
  const res = [];
  const F = (indices.length / 3) | 0;
  for (let f = 0; f < F; f++) {
    const a = indices[3 * f], b = indices[3 * f + 1], c = indices[3 * f + 2];
    if (a === b || b === c || c === a) { res.push(f); continue; }
    const l2 = Math.max(dist2(positions, a, b), dist2(positions, b, c), dist2(positions, c, a));
    const A2 = area2(positions, a, b, c);
    if (!(A2 > relTol * l2)) res.push(f);
  }
  return res;
}

function dist2(P, a, b) {
  const x = P[3 * a] - P[3 * b], y = P[3 * a + 1] - P[3 * b + 1], z = P[3 * a + 2] - P[3 * b + 2];
  return x * x + y * y + z * z;
}

/**
 * Make the face winding consistent by flood fill across two-face edges, then turn every closed component
 * outward (positive divergence-theorem volume). Open components keep the orientation of their seed face.
 * @returns {{ indices: Uint32Array, flipped: number, orientable: boolean, components: number,
 *   componentVolume: number[], componentClosed: boolean[] }} new index array, number of faces whose winding
 *   changed, whether a consistent orientation exists (false for Möbius-like or non-manifold conflicts).
 */
export function orientConsistently(positions, indices, nVerts = (positions.length / 3) | 0) {
  const F = (indices.length / 3) | 0;
  const out = Uint32Array.from(indices);
  const edges = buildEdges(indices, nVerts);
  const comp = new Int32Array(F).fill(-1);
  const flip = new Uint8Array(F);
  let orientable = true;
  const stack = [];
  const vol = [], closed = [];
  let nc = 0;
  const dir = (h) => indices[h] < indices[3 * ((h / 3) | 0) + (h % 3 + 1) % 3] ? 1 : 0;
  for (let seed = 0; seed < F; seed++) {
    if (comp[seed] >= 0) continue;
    comp[seed] = nc; stack.push(seed);
    while (stack.length) {
      const f = stack.pop();
      for (let k = 0; k < 3; k++) {
        const e = edges.faceEdges[3 * f + k];
        if (e === NONE) continue;
        const s = edges.edgeFaceStart[e];
        if (edges.edgeFaceStart[e + 1] - s !== 2) continue;
        const h0 = edges.edgeHalfEdges[s], h1 = edges.edgeHalfEdges[s + 1];
        const hOther = ((h0 / 3) | 0) === f && h0 === 3 * f + k ? h1 : h0;
        const g = (hOther / 3) | 0;
        if (g === f) continue;
        // consistent neighbours traverse the shared edge in opposite directions
        const same = dir(3 * f + k) === dir(hOther);
        const want = same ? 1 - flip[f] : flip[f];
        if (comp[g] < 0) { comp[g] = nc; flip[g] = want; stack.push(g); }
        else if (flip[g] !== want) orientable = false;
      }
    }
    nc++;
  }
  for (let c = 0; c < nc; c++) { vol.push(0); closed.push(true); }
  for (let e = 0; e < edges.edgeCount; e++) {
    if (edges.edgeFaceStart[e + 1] - edges.edgeFaceStart[e] !== 2) {
      for (let q = edges.edgeFaceStart[e]; q < edges.edgeFaceStart[e + 1]; q++) closed[comp[(edges.edgeHalfEdges[q] / 3) | 0]] = false;
    }
  }
  const P = positions;
  for (let f = 0; f < F; f++) {
    const a = 3 * indices[3 * f], b = 3 * indices[3 * f + 1], c = 3 * indices[3 * f + 2];
    const v = (P[a] * (P[b + 1] * P[c + 2] - P[b + 2] * P[c + 1]) - P[a + 1] * (P[b] * P[c + 2] - P[b + 2] * P[c]) + P[a + 2] * (P[b] * P[c + 1] - P[b + 1] * P[c])) / 6;
    vol[comp[f]] += flip[f] ? -v : v;
  }
  let flipped = 0;
  for (let f = 0; f < F; f++) {
    let fl = flip[f];
    if (closed[comp[f]] && vol[comp[f]] < 0) fl = 1 - fl;
    if (fl) { const t = out[3 * f + 1]; out[3 * f + 1] = out[3 * f + 2]; out[3 * f + 2] = t; flipped++; }
  }
  return { indices: out, flipped, orientable, components: nc, componentVolume: vol.map((v, c) => (closed[c] ? Math.abs(v) : v)), componentClosed: closed };
}

/**
 * Euler characteristic χ = V − E + F and genus. For an orientable surface with C components and b boundary
 * loops, χ = 2C − 2g − b (Euler–Poincaré), so g = (2C − b − χ) / 2 (total over all components).
 * @returns {{ V: number, E: number, F: number, chi: number, components: number, boundaryLoops: number, genus: number }}
 *   V counts referenced vertices only; genus is NaN when the formula does not give a non-negative integer
 *   (non-manifold input).
 */
export function eulerCharacteristic(indices, nVerts, edges = buildEdges(indices, nVerts)) {
  const used = new Uint8Array(nVerts);
  let V = 0;
  for (let i = 0; i < indices.length; i++) if (!used[indices[i]]) { used[indices[i]] = 1; V++; }
  const F = (indices.length / 3) | 0;
  const E = edges.edgeCount;
  const chi = V - E + F;
  const components = connectedComponents(indices, nVerts, 'vertex').count;
  const b = boundaryLoops(indices, nVerts, edges).loops.length;
  const g = (2 * components - b - chi) / 2;
  return { V, E, F, chi, components, boundaryLoops: b, genus: Number.isInteger(g) && g >= 0 ? g : NaN };
}

/**
 * Boundary faces of a volume mesh: cell faces that belong to exactly one cell, outward oriented.
 * @param {object} model volume model
 * @returns {{ tris: Uint32Array, quads: Uint32Array, triCell: Uint32Array, quadCell: Uint32Array,
 *   interiorFaces: number, overSharedFaces: number }} node ids of boundary triangles / quadrilaterals, the global
 *   cell index (running through model.cells in order) owning each, the number of faces shared by two cells and
 *   the number shared by more than two (non-conforming input).
 */
export function boundaryFaces(model) {
  let total = 0;
  for (const block of model.cells || []) {
    const info = CELL_INFO[block.type];
    if (!info) throw new RangeError(`cells.type: unsupported cell type '${block.type}'`);
    total += ((block.conn.length / info.nodes) | 0) * CELL_FACES[block.type].length;
  }
  // open-addressing hash of faces keyed by their sorted node ids (typed arrays only: meshes of millions of cells)
  let size = 16;
  while (size < 2 * total) size *= 2;
  const mask = size - 1;
  const table = new Int32Array(size).fill(-1);
  const sorted = new Uint32Array(4 * total), ids = new Uint32Array(4 * total), owner = new Uint32Array(total), count = new Uint8Array(total), nodes = new Uint8Array(total);
  let nRec = 0, cellBase = 0;
  const q = [0, 0, 0, 0];
  for (const block of model.cells || []) {
    const ft = CELL_FACES[block.type];
    const nn = CELL_INFO[block.type].nodes;
    const nc = (block.conn.length / nn) | 0;
    for (let c = 0; c < nc; c++) {
      for (let f = 0; f < ft.length; f++) {
        const lf = ft[f], m = lf.length;
        for (let k = 0; k < m; k++) q[k] = block.conn[c * nn + lf[k]];
        let s0 = q[0], s1 = q[1], s2 = q[2], s3 = m === 4 ? q[3] : NONE, t;
        if (s0 > s1) { t = s0; s0 = s1; s1 = t; }
        if (s2 > s3) { t = s2; s2 = s3; s3 = t; }
        if (s0 > s2) { t = s0; s0 = s2; s2 = t; }
        if (s1 > s3) { t = s1; s1 = s3; s3 = t; }
        if (s1 > s2) { t = s1; s1 = s2; s2 = t; }
        let h = (Math.imul(s0, 73856093) ^ Math.imul(s1, 19349663) ^ Math.imul(s2, 83492791) ^ Math.imul(s3, 2654435761)) & mask;
        for (;;) {
          const r = table[h];
          if (r < 0) {
            table[h] = nRec;
            sorted[4 * nRec] = s0; sorted[4 * nRec + 1] = s1; sorted[4 * nRec + 2] = s2; sorted[4 * nRec + 3] = s3;
            for (let k = 0; k < m; k++) ids[4 * nRec + k] = q[k];
            owner[nRec] = cellBase + c; count[nRec] = 1; nodes[nRec] = m;
            nRec++;
            break;
          }
          if (sorted[4 * r] === s0 && sorted[4 * r + 1] === s1 && sorted[4 * r + 2] === s2 && sorted[4 * r + 3] === s3) { if (count[r] < 255) count[r]++; break; }
          h = (h + 1) & mask;
        }
      }
    }
    cellBase += nc;
  }
  let nt = 0, nq = 0, interior = 0, over = 0;
  for (let r = 0; r < nRec; r++) { if (count[r] === 1) { if (nodes[r] === 3) nt++; else nq++; } else if (count[r] === 2) interior++; else over++; }
  const tris = new Uint32Array(3 * nt), quads = new Uint32Array(4 * nq), triCell = new Uint32Array(nt), quadCell = new Uint32Array(nq);
  nt = 0; nq = 0;
  for (let r = 0; r < nRec; r++) {
    if (count[r] !== 1) continue;
    if (nodes[r] === 3) { tris[3 * nt] = ids[4 * r]; tris[3 * nt + 1] = ids[4 * r + 1]; tris[3 * nt + 2] = ids[4 * r + 2]; triCell[nt++] = owner[r]; }
    else { quads[4 * nq] = ids[4 * r]; quads[4 * nq + 1] = ids[4 * r + 1]; quads[4 * nq + 2] = ids[4 * r + 2]; quads[4 * nq + 3] = ids[4 * r + 3]; quadCell[nq++] = owner[r]; }
  }
  return { tris, quads, triCell, quadCell, interiorFaces: interior, overSharedFaces: over };
}

/**
 * Signed volume of one cell from its corner nodes (sum of tetrahedra fanned from the first corner over the
 * outward faces). Positive for correctly ordered cells.
 */
export function cellVolume(positions, conn, offset, type) {
  const ft = CELL_FACES[type];
  const P = positions;
  const o = 3 * conn[offset];
  let v = 0;
  for (let q = 0; q < ft.length; q++) {
    const lf = ft[q];
    for (let t = 1; t + 1 < lf.length; t++) {
      const a = 3 * conn[offset + lf[0]], b = 3 * conn[offset + lf[t]], c = 3 * conn[offset + lf[t + 1]];
      const ax = P[a] - P[o], ay = P[a + 1] - P[o + 1], az = P[a + 2] - P[o + 2];
      const bx = P[b] - P[o], by = P[b + 1] - P[o + 1], bz = P[b + 2] - P[o + 2];
      const cx = P[c] - P[o], cy = P[c + 1] - P[o + 1], cz = P[c + 2] - P[o + 2];
      v += ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    }
  }
  return v / 6;
}

/**
 * Triangle surface of any model that has one: a surface is returned as is, a volume gives its boundary
 * (quadrilaterals split on the 0–2 diagonal), a grid gives its height-field surface (cells touching NaN skipped).
 * @returns {object|null} surface model, with `faceCell` (Uint32Array) for volumes; null for other kinds.
 */
export function toSurface(model) {
  if (!model) return null;
  if (model.kind === 'surface') return model;
  if (model.kind === 'volume') {
    const bf = boundaryFaces(model);
    const nt = bf.tris.length / 3, nq = bf.quads.length / 4;
    const indices = new Uint32Array(3 * (nt + 2 * nq));
    const faceCell = new Uint32Array(nt + 2 * nq);
    indices.set(bf.tris);
    faceCell.set(bf.triCell);
    for (let q = 0; q < nq; q++) {
      const a = bf.quads[4 * q], b = bf.quads[4 * q + 1], c = bf.quads[4 * q + 2], d = bf.quads[4 * q + 3];
      indices.set([a, b, c, a, c, d], 3 * (nt + 2 * q));
      faceCell[nt + 2 * q] = faceCell[nt + 2 * q + 1] = bf.quadCell[q];
    }
    return { kind: 'surface', name: model.name, positions: model.positions, indices, faceCell };
  }
  if (model.kind === 'grid') {
    const { nx, ny } = model;
    const positions = new Float32Array(3 * nx * ny);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const p = 3 * (j * nx + i);
      positions[p] = model.x0 + i * model.dx; positions[p + 1] = model.y0 + j * model.dy; positions[p + 2] = model.z[j * nx + i];
    }
    const idx = [];
    // counter-clockwise seen from +z when dx·dy > 0
    const ccw = model.dx * model.dy > 0;
    for (let j = 0; j + 1 < ny; j++) for (let i = 0; i + 1 < nx; i++) {
      const a = j * nx + i, b = a + 1, c = a + nx + 1, d = a + nx;
      if (Number.isNaN(model.z[a]) || Number.isNaN(model.z[b]) || Number.isNaN(model.z[c]) || Number.isNaN(model.z[d])) continue;
      if (ccw) idx.push(a, b, c, a, c, d); else idx.push(a, c, b, a, d, c);
    }
    return { kind: 'surface', name: model.name, positions, indices: Uint32Array.from(idx) };
  }
  return null;
}

function segTri(P, p, q, a, b, c, eps) {
  // strict interior intersection of segment p–q with triangle a, b, c (Möller–Trumbore)
  const e1x = P[b] - P[a], e1y = P[b + 1] - P[a + 1], e1z = P[b + 2] - P[a + 2];
  const e2x = P[c] - P[a], e2y = P[c + 1] - P[a + 1], e2z = P[c + 2] - P[a + 2];
  const dx = P[q] - P[p], dy = P[q + 1] - P[p + 1], dz = P[q + 2] - P[p + 2];
  const px = dy * e2z - dz * e2y, py = dz * e2x - dx * e2z, pz = dx * e2y - dy * e2x;
  const det = e1x * px + e1y * py + e1z * pz;
  const scale = Math.sqrt((e1x * e1x + e1y * e1y + e1z * e1z) * (e2x * e2x + e2y * e2y + e2z * e2z) * (dx * dx + dy * dy + dz * dz));
  if (Math.abs(det) <= 1e-12 * scale) return false;
  const inv = 1 / det;
  const tx = P[p] - P[a], ty = P[p + 1] - P[a + 1], tz = P[p + 2] - P[a + 2];
  const u = (tx * px + ty * py + tz * pz) * inv;
  if (u < eps || u > 1 - eps) return false;
  const qx = ty * e1z - tz * e1y, qy = tz * e1x - tx * e1z, qz = tx * e1y - ty * e1x;
  const v = (dx * qx + dy * qy + dz * qz) * inv;
  if (v < eps || u + v > 1 - eps) return false;
  const t = (e2x * qx + e2y * qy + e2z * qz) * inv;
  return t > eps && t < 1 - eps;
}

/**
 * Self-intersection search: pairs of faces that do not share a vertex and whose interiors cross
 * (an edge of one passes strictly through the other). Uniform-grid broad phase.
 * Coplanar overlaps and intersections between faces that share a vertex or an edge are not detected.
 * @param {{ maxPairs?: number }} [opts] stop after this many pairs (default 1000)
 * @returns {{ pairs: number[][], truncated: boolean }} face-id pairs
 */
export function selfIntersections(positions, indices, { maxPairs = 1000 } = {}) {
  const F = (indices.length / 3) | 0;
  const pairs = [];
  if (F < 2) return { pairs, truncated: false };
  const bb = bounds(positions);
  const fb = new Float64Array(6 * F);
  let ext = 0;
  for (let f = 0; f < F; f++) {
    let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let k = 0; k < 3; k++) {
      const v = 3 * indices[3 * f + k];
      x0 = Math.min(x0, positions[v]); x1 = Math.max(x1, positions[v]);
      y0 = Math.min(y0, positions[v + 1]); y1 = Math.max(y1, positions[v + 1]);
      z0 = Math.min(z0, positions[v + 2]); z1 = Math.max(z1, positions[v + 2]);
    }
    fb.set([x0, y0, z0, x1, y1, z1], 6 * f);
    ext += Math.max(x1 - x0, y1 - y0, z1 - z0);
  }
  // cubic cells about twice the mean face size, whatever the proportions of the body (a long slender hull must
  // not put all its faces into a few slab-shaped cells); at most 1024 cells along any axis
  const h = Math.max(2 * ext / F, bb.diagonal / 1024, 1e-300);
  const nx = Math.max(1, Math.ceil((bb.max[0] - bb.min[0]) / h)), ny = Math.max(1, Math.ceil((bb.max[1] - bb.min[1]) / h)), nz = Math.max(1, Math.ceil((bb.max[2] - bb.min[2]) / h));
  const cells = new Map();
  for (let f = 0; f < F; f++) {
    const i0 = Math.min(nx - 1, Math.floor((fb[6 * f] - bb.min[0]) / h)), i1 = Math.min(nx - 1, Math.floor((fb[6 * f + 3] - bb.min[0]) / h));
    const j0 = Math.min(ny - 1, Math.floor((fb[6 * f + 1] - bb.min[1]) / h)), j1 = Math.min(ny - 1, Math.floor((fb[6 * f + 4] - bb.min[1]) / h));
    const k0 = Math.min(nz - 1, Math.floor((fb[6 * f + 2] - bb.min[2]) / h)), k1 = Math.min(nz - 1, Math.floor((fb[6 * f + 5] - bb.min[2]) / h));
    for (let i = i0; i <= i1; i++) for (let j = j0; j <= j1; j++) for (let k = k0; k <= k1; k++) {
      const key = (i * ny + j) * nz + k;
      const l = cells.get(key);
      if (l) l.push(f); else cells.set(key, [f]);
    }
  }
  const seen = new Set();
  const eps = 1e-9;
  for (const list of cells.values()) {
    for (let i = 0; i < list.length; i++) {
      const f = list[i];
      for (let j = i + 1; j < list.length; j++) {
        const g = list[j];
        if (fb[6 * f] > fb[6 * g + 3] || fb[6 * g] > fb[6 * f + 3] || fb[6 * f + 1] > fb[6 * g + 4] || fb[6 * g + 1] > fb[6 * f + 4] || fb[6 * f + 2] > fb[6 * g + 5] || fb[6 * g + 2] > fb[6 * f + 5]) continue;
        const fa = indices[3 * f], fbv = indices[3 * f + 1], fc = indices[3 * f + 2];
        const ga = indices[3 * g], gb = indices[3 * g + 1], gc = indices[3 * g + 2];
        if (fa === ga || fa === gb || fa === gc || fbv === ga || fbv === gb || fbv === gc || fc === ga || fc === gb || fc === gc) continue;
        const key = f < g ? f * F + g : g * F + f;
        if (seen.has(key)) continue;
        seen.add(key);
        const A = [3 * fa, 3 * fbv, 3 * fc], B = [3 * ga, 3 * gb, 3 * gc];
        let hit = false;
        for (let k = 0; k < 3 && !hit; k++) {
          hit = segTri(positions, A[k], A[(k + 1) % 3], B[0], B[1], B[2], eps) || segTri(positions, B[k], B[(k + 1) % 3], A[0], A[1], A[2], eps);
        }
        if (hit) {
          pairs.push([Math.min(f, g), Math.max(f, g)]);
          if (pairs.length >= maxPairs) return { pairs, truncated: true };
        }
      }
    }
  }
  return { pairs, truncated: false };
}

const CAP = 5000;
const issue = (list, severity, code, count, message, extra) => { if (count > 0) list.push({ code, count, severity, message, ...extra }); };
const plural = (n, s, p = `${s}s`) => `${n} ${n === 1 ? s : p}`;

function inspectSurface(model, opts) {
  const issues = [];
  const P = model.positions;
  const nV0 = (P.length / 3) | 0;
  const F = (model.indices.length / 3) | 0;
  let bad = 0;
  for (let i = 0; i < P.length; i++) if (!Number.isFinite(P[i])) bad++;
  issue(issues, 'error', 'non-finite-coordinates', bad, `${plural(bad, 'coordinate')} are NaN or infinite.`);
  let oob = 0;
  for (let i = 0; i < model.indices.length; i++) if (model.indices[i] >= nV0) oob++;
  issue(issues, 'error', 'index-out-of-range', oob, `${plural(oob, 'face corner')} refer to vertices that do not exist.`);
  const bb = bounds(P);
  if (oob || bad || F === 0) {
    if (F === 0) issue(issues, 'error', 'empty', 1, 'The surface has no faces.');
    return { ok: false, watertight: false, issues, stats: { vertices: nV0, faces: F, bbox: bb } };
  }
  const tol = opts.tolerance ?? 1e-6 * bb.diagonal;
  const w = opts.weld === false ? null : weldVertices(P, tol);
  let idx = model.indices, pos = P, rep = null;
  if (w && w.merged > 0) {
    idx = new Uint32Array(model.indices.length);
    for (let i = 0; i < idx.length; i++) idx[i] = w.map[model.indices[i]];
    pos = w.positions; rep = w.rep;
    issue(issues, 'info', 'unwelded-vertices', w.merged, `${plural(w.merged, 'vertex', 'vertices')} coincide with another vertex within ${tol.toExponential(1)} m (typical of STL); topology was evaluated after welding them.`);
  }
  const nV = (pos.length / 3) | 0;
  const orig = (v) => (rep ? rep[v] : v);
  const deg = degenerateFaces(pos, idx);
  issue(issues, 'error', 'degenerate-faces', deg.length, `${plural(deg.length, 'face')} have zero area (repeated or collinear corners).`, { elements: deg.slice(0, CAP) });
  const dup = duplicateFaces(idx, nV);
  issue(issues, 'error', 'duplicate-faces', dup.length, `${plural(dup.length, 'face')} repeat another face.`, { elements: dup.slice(0, CAP) });
  const edges = buildEdges(idx, nV);
  const loops = boundaryLoops(idx, nV, edges);
  const bEdges = [];
  for (const l of loops.loops) for (let i = 0; i < l.length && bEdges.length < CAP; i++) bEdges.push([orig(l[i]), orig(l[(i + 1) % l.length])]);
  issue(issues, 'warn', 'open-boundary', loops.edgeCount, `${plural(loops.edgeCount, 'boundary edge')} in ${plural(loops.loops.length, 'loop')}: the surface is open (holes or free edges).`, { edges: bEdges, loops: loops.loops.length });
  const nme = nonManifoldEdges(idx, nV, edges);
  issue(issues, 'error', 'non-manifold-edges', nme.length, `${plural(nme.length, 'edge')} are shared by more than two faces.`, { edges: nme.slice(0, CAP).map(([a, b]) => [orig(a), orig(b)]) });
  const nmv = nonManifoldVertices(idx, nV, edges);
  issue(issues, 'error', 'non-manifold-vertices', nmv.length, `${plural(nmv.length, 'vertex', 'vertices')} join otherwise separate parts of the surface (bow-tie).`, { vertices: nmv.slice(0, CAP).map(orig) });
  const used = new Uint8Array(nV);
  for (let i = 0; i < idx.length; i++) used[idx[i]] = 1;
  let unref = 0;
  for (let v = 0; v < nV; v++) if (!used[v]) unref++;
  issue(issues, 'info', 'unreferenced-vertices', unref, `${plural(unref, 'vertex', 'vertices')} are not used by any face.`);
  const ori = orientConsistently(pos, idx, nV);
  if (!ori.orientable) issue(issues, 'error', 'not-orientable', 1, 'Face windings cannot be made consistent (one-sided or non-manifold surface).');
  else issue(issues, 'warn', 'inconsistent-orientation', ori.flipped, `${plural(ori.flipped, 'face')} are wound against the outward-normal convention (reversed normals).`);
  const comps = connectedComponents(idx, nV, 'edge', edges);
  issue(issues, 'info', 'multiple-components', comps.count > 1 ? comps.count : 0, `The surface has ${comps.count} separate parts.`);
  const eu = eulerCharacteristic(idx, nV, edges);
  let nSelf = 0;
  const doSelf = opts.selfIntersections ?? F <= 60000;
  if (doSelf) {
    const si = selfIntersections(pos, idx);
    nSelf = si.pairs.length;
    issue(issues, 'error', 'self-intersections', nSelf, `${si.truncated ? 'At least ' : ''}${plural(nSelf, 'pair')} of faces cross each other.`, { elements: [...new Set(si.pairs.flat())].slice(0, CAP) });
  }
  const closed = loops.edgeCount === 0;
  const manifold = nme.length === 0 && nmv.length === 0;
  const reasons = [];
  if (!closed) reasons.push(`open: ${plural(loops.edgeCount, 'boundary edge')} in ${plural(loops.loops.length, 'loop')}`);
  if (nme.length) reasons.push(`${plural(nme.length, 'non-manifold edge')}`);
  if (nmv.length) reasons.push(`${plural(nmv.length, 'non-manifold vertex', 'non-manifold vertices')}`);
  if (!ori.orientable) reasons.push('not orientable');
  if (deg.length) reasons.push(`${plural(deg.length, 'degenerate face')}`);
  if (dup.length) reasons.push(`${plural(dup.length, 'duplicate face')}`);
  if (nSelf) reasons.push(`${plural(nSelf, 'self-intersecting face pair')}`);
  const watertight = reasons.length === 0;
  let area = 0;
  for (let f = 0; f < F; f++) area += area2(pos, idx[3 * f], idx[3 * f + 1], idx[3 * f + 2]) / 2;
  const volume = ori.componentVolume.reduce((s, v, c) => s + (ori.componentClosed[c] ? v : 0), 0);
  return {
    ok: !issues.some((i) => i.severity === 'error'),
    watertight, reasons, issues,
    stats: {
      vertices: nV0, weldedVertices: nV, faces: F, edges: eu.E, eulerCharacteristic: eu.chi, genus: eu.genus,
      components: comps.count, boundaryLoops: loops.loops.length, boundaryEdges: loops.edgeCount, closed, manifold,
      orientable: ori.orientable, outward: ori.orientable && ori.flipped === 0, area, volume: closed ? volume : NaN,
      bbox: bb, selfIntersectionsChecked: doSelf,
    },
  };
}

function inspectVolume(model) {
  const issues = [];
  const P = model.positions;
  const nV = (P.length / 3) | 0;
  let nCells = 0, inverted = 0, zero = 0, oob = 0, volume = 0;
  const byType = {};
  const used = new Uint8Array(nV);
  const invertedIds = [];
  for (const block of model.cells || []) {
    const info = CELL_INFO[block.type];
    if (!info) { issue(issues, 'error', 'unsupported-cell-type', 1, `Cell type '${block.type}' is not supported.`); continue; }
    const nc = (block.conn.length / info.nodes) | 0;
    if (block.conn.length % info.nodes) issue(issues, 'error', 'bad-connectivity-length', 1, `Connectivity of the ${block.type} block is not a multiple of ${info.nodes}.`);
    byType[block.type] = (byType[block.type] || 0) + nc;
    for (let c = 0; c < nc; c++) {
      let ok = true;
      for (let k = 0; k < info.nodes; k++) { const v = block.conn[c * info.nodes + k]; if (v >= nV) { ok = false; oob++; } else used[v] = 1; }
      if (!ok) continue;
      const v = cellVolume(P, block.conn, c * info.nodes, block.type);
      volume += v;
      if (v < 0) { inverted++; if (invertedIds.length < CAP) invertedIds.push(nCells + c); } else if (v === 0 || !Number.isFinite(v)) zero++;
    }
    nCells += nc;
  }
  issue(issues, 'error', 'index-out-of-range', oob, `${plural(oob, 'cell node')} refer to nodes that do not exist.`);
  issue(issues, 'error', 'inverted-cells', inverted, `${plural(inverted, 'cell')} have negative volume (node order reversed).`, { elements: invertedIds });
  issue(issues, 'error', 'zero-volume-cells', zero, `${plural(zero, 'cell')} have zero volume.`);
  let unref = 0;
  for (let v = 0; v < nV; v++) if (!used[v]) unref++;
  issue(issues, 'info', 'unreferenced-nodes', unref, `${plural(unref, 'node')} are not used by any cell.`);
  if (nCells === 0) issue(issues, 'error', 'empty', 1, 'The volume mesh has no cells.');
  let boundary = null;
  if (!oob && nCells) {
    const bf = boundaryFaces(model);
    issue(issues, 'error', 'over-shared-faces', bf.overSharedFaces, `${plural(bf.overSharedFaces, 'face')} are shared by more than two cells.`);
    const surf = toSurface(model);
    const edges = buildEdges(surf.indices, nV);
    // a conforming mesh has a closed two-manifold boundary. A hanging node leaves the unmatched coarse and fine
    // faces back to back, so their edges carry one, three or four boundary faces instead of two
    let odd = 0;
    for (let e = 0; e < edges.edgeCount; e++) if (edges.edgeFaceStart[e + 1] - edges.edgeFaceStart[e] !== 2) odd++;
    issue(issues, 'error', 'non-conforming-boundary', odd, `${plural(odd, 'edge')} of the mesh boundary do not join exactly two boundary faces: cells do not match face to face (hanging nodes, mixed face types, or cells touching only along an edge).`);
    boundary = { triangles: bf.tris.length / 3, quads: bf.quads.length / 4, interiorFaces: bf.interiorFaces, components: connectedComponents(surf.indices, nV, 'edge', edges).count };
  }
  return {
    ok: !issues.some((i) => i.severity === 'error'), watertight: !!boundary && !issues.some((i) => i.code === 'non-conforming-boundary'), issues,
    stats: { nodes: nV, cells: nCells, cellsByType: byType, volume, boundary, bbox: bounds(P), sets: Object.keys(model.sets || {}) },
  };
}

function inspectFrame(model, opts) {
  const issues = [];
  const N = model.nodes, nN = (N.length / 3) | 0, members = model.members || [];
  const bb = bounds(N);
  const tol = opts.tolerance ?? 1e-6 * (bb.diagonal || 1);
  let zeroLen = 0, oob = 0, dupM = 0;
  const deg = new Uint32Array(nN);
  const seen = new Set();
  const p = new Uint32Array(nN);
  for (let i = 0; i < nN; i++) p[i] = i;
  const zeroIds = [], dupIds = [];
  let total = 0, lmin = Infinity, lmax = 0;
  members.forEach((m, i) => {
    if (!(m.a < nN) || !(m.b < nN) || m.a < 0 || m.b < 0) { oob++; return; }
    const L = Math.sqrt(dist2(N, m.a, m.b));
    if (m.a === m.b || L <= tol) { zeroLen++; zeroIds.push(i); return; }
    total += L; lmin = Math.min(lmin, L); lmax = Math.max(lmax, L);
    const key = m.a < m.b ? `${m.a},${m.b}` : `${m.b},${m.a}`;
    if (seen.has(key)) { dupM++; dupIds.push(i); } else seen.add(key);
    deg[m.a]++; deg[m.b]++;
    p[ufFind(p, m.a)] = ufFind(p, m.b);
  });
  issue(issues, 'error', 'index-out-of-range', oob, `${plural(oob, 'member')} refer to nodes that do not exist.`);
  issue(issues, 'error', 'zero-length-members', zeroLen, `${plural(zeroLen, 'member')} have zero length.`, { elements: zeroIds });
  issue(issues, 'error', 'duplicate-members', dupM, `${plural(dupM, 'member')} connect the same two nodes as another member.`, { elements: dupIds });
  const w = weldVertices(N, tol);
  issue(issues, 'warn', 'coincident-nodes', w.merged, `${plural(w.merged, 'node')} coincide with another node within ${tol.toExponential(1)} m: members meeting there are not connected.`);
  let unused = 0, free = 0;
  const roots = new Set();
  for (let i = 0; i < nN; i++) {
    if (deg[i] === 0) unused++; else { roots.add(ufFind(p, i)); if (deg[i] === 1) free++; }
  }
  issue(issues, 'info', 'unreferenced-nodes', unused, `${plural(unused, 'node')} are not used by any member.`);
  issue(issues, 'info', 'free-ends', free, `${plural(free, 'node')} have a single member (supports, tips or dangling members — check each is intended).`);
  issue(issues, 'warn', 'disconnected-parts', roots.size > 1 ? roots.size : 0, `The frame is in ${roots.size} unconnected parts.`);
  const noSec = members.filter((m) => !m.section).length;
  issue(issues, 'info', 'members-without-section', noSec, `${plural(noSec, 'member')} have no section assigned.`);
  return {
    ok: !issues.some((i) => i.severity === 'error'), watertight: false, issues,
    stats: { nodes: nN, members: members.length, components: roots.size, freeEnds: free, totalLength: total, minLength: Number.isFinite(lmin) ? lmin : 0, maxLength: lmax, bbox: bb },
  };
}

function segmentsCross(ax, ay, bx, by, cx, cy, dx, dy) {
  const o = (px, py, qx, qy, rx, ry) => (qx - px) * (ry - py) - (qy - py) * (rx - px);
  const d1 = o(cx, cy, dx, dy, ax, ay), d2 = o(cx, cy, dx, dy, bx, by), d3 = o(ax, ay, bx, by, cx, cy), d4 = o(ax, ay, bx, by, dx, dy);
  return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}

function inspectOther(model) {
  const issues = [];
  const stats = {};
  if (model.kind === 'polyline') {
    const P = model.points, n = (P.length / 3) | 0;
    let bad = 0, dup = 0, len = 0;
    for (let i = 0; i < P.length; i++) if (!Number.isFinite(P[i])) bad++;
    for (let i = 1; i < n; i++) { const d = Math.sqrt(dist2(P, i, i - 1)); len += d; if (d === 0) dup++; }
    issue(issues, 'error', 'non-finite-coordinates', bad, `${plural(bad, 'coordinate')} are NaN or infinite.`);
    issue(issues, 'warn', 'repeated-points', dup, `${plural(dup, 'point')} repeat the previous point.`);
    if (n < 2) issue(issues, 'error', 'too-few-points', 1, 'A polyline needs at least two points.');
    Object.assign(stats, { points: n, length: len, bbox: bounds(P) });
  } else if (model.kind === 'points') {
    const P = model.positions, n = (P.length / 3) | 0;
    let bad = 0;
    for (let i = 0; i < P.length; i++) if (!Number.isFinite(P[i])) bad++;
    issue(issues, 'error', 'non-finite-coordinates', bad, `${plural(bad, 'coordinate')} are NaN or infinite.`);
    if (n === 0) issue(issues, 'error', 'empty', 1, 'The point set is empty.');
    const dupes = bad || n === 0 ? 0 : weldVertices(P, 0).merged;
    issue(issues, 'info', 'duplicate-points', dupes, `${plural(dupes, 'point')} repeat another point exactly.`);
    Object.assign(stats, { points: n, bbox: bounds(P) });
  } else if (model.kind === 'grid') {
    const n = model.nx * model.ny;
    if (!model.z || model.z.length !== n) issue(issues, 'error', 'bad-grid-size', 1, `The grid has ${model.z ? model.z.length : 0} values but nx·ny = ${n}.`);
    if (!(model.dx !== 0) || !(model.dy !== 0)) issue(issues, 'error', 'zero-spacing', 1, 'Grid spacing dx or dy is zero.');
    let nodata = 0, zmin = Infinity, zmax = -Infinity;
    for (let i = 0; i < (model.z ? model.z.length : 0); i++) { const v = model.z[i]; if (Number.isFinite(v)) { if (v < zmin) zmin = v; if (v > zmax) zmax = v; } else nodata++; }
    issue(issues, 'info', 'no-data-cells', nodata, `${plural(nodata, 'grid node')} have no data.`);
    if (nodata === n) issue(issues, 'error', 'empty', 1, 'The grid has no valid values.');
    Object.assign(stats, { nx: model.nx, ny: model.ny, nodes: n, noData: nodata, zMin: zmin, zMax: zmax, crs: model.crs || null });
  } else if (model.kind === 'section') {
    let pts = 0, open = 0, tiny = 0, cross = 0, area = 0;
    const loops = model.loops || [];
    for (const L of loops) {
      const n = (L.length / 2) | 0;
      pts += n;
      if (n < 3) { tiny++; continue; }
      let a = 0;
      for (let i = 0; i < n; i++) { const j = (i + 1) % n; a += L[2 * i] * L[2 * j + 1] - L[2 * j] * L[2 * i + 1]; }
      area += Math.abs(a) / 2 * (loops.indexOf(L) === 0 ? 1 : 0);
      if (L[0] === L[2 * n - 2] && L[1] === L[2 * n - 1]) open++; // repeated closing point
      if (n <= 2000) {
        for (let i = 0; i < n; i++) for (let j = i + 2; j < n; j++) {
          if (i === 0 && j === n - 1) continue;
          const i2 = (i + 1) % n, j2 = (j + 1) % n;
          if (segmentsCross(L[2 * i], L[2 * i + 1], L[2 * i2], L[2 * i2 + 1], L[2 * j], L[2 * j + 1], L[2 * j2], L[2 * j2 + 1])) cross++;
        }
      }
    }
    if (!loops.length) issue(issues, 'error', 'empty', 1, 'The section has no loops.');
    issue(issues, 'error', 'degenerate-loops', tiny, `${plural(tiny, 'loop')} have fewer than three points.`);
    issue(issues, 'info', 'repeated-closing-point', open, `${plural(open, 'loop')} repeat their first point at the end.`);
    issue(issues, 'error', 'self-intersecting-loops', cross, `${plural(cross, 'pair')} of loop segments cross.`);
    Object.assign(stats, { loops: loops.length, points: pts, outerArea: area });
  } else {
    issue(issues, 'error', 'unknown-kind', 1, `Model kind '${model.kind}' is not recognised.`);
  }
  return { ok: !issues.some((i) => i.severity === 'error'), watertight: false, issues, stats };
}

/**
 * Inspect any geometry model: validity, connectivity and, for surfaces, the watertightness verdict.
 * @param {object} model CONTRACT §4 model
 * @param {{ tolerance?: number, weld?: boolean, selfIntersections?: boolean }} [opts] weld tolerance [m]
 *   (default 1e-6 × bounding-box diagonal); weld = false inspects the connectivity exactly as stored;
 *   self-intersection search defaults to on for surfaces up to 60 000 faces.
 * @returns {{ ok: boolean, watertight: boolean, reasons?: string[], issues: Array<{ code: string, count: number,
 *   severity: 'error'|'warn'|'info', message: string, edges?: number[][], elements?: number[], vertices?: number[] }>,
 *   stats: object }} `ok` is true when no issue has severity 'error'. `edges`, `elements` and `vertices` (first
 *   5000) index the model as given and can be passed to the viewer's highlight.
 */
export function inspect(model, opts = {}) {
  if (!model || typeof model !== 'object') throw new RangeError('model: a geometry model is required');
  let res;
  if (model.kind === 'surface') { requireSurface(model); res = inspectSurface(model, opts); }
  else if (model.kind === 'volume') res = inspectVolume(model);
  else if (model.kind === 'frame') res = inspectFrame(model, opts);
  else res = inspectOther(model);
  res.kind = model.kind;
  res.name = model.name;
  return res;
}

/**
 * Watertightness verdict for a surface (or the boundary of a volume mesh).
 * @returns {{ watertight: boolean, reasons: string[] }}
 */
export function watertight(model, opts = {}) {
  const s = toSurface(model);
  if (!s) return { watertight: false, reasons: [`a ${model && model.kind} model has no surface`] };
  const r = inspectSurface(s, opts);
  return { watertight: r.watertight, reasons: r.reasons || r.issues.map((i) => i.message) };
}

/** Per-face group id (−1 where a face is in no group) and the group names. */
export function faceGroupIds(model) {
  const F = (model.indices.length / 3) | 0;
  const ids = new Int32Array(F).fill(-1);
  const names = (model.groups || []).map((g) => g.name);
  (model.groups || []).forEach((g, gi) => {
    // group ranges are in index units (start, count over the index array), as in WebGL draw ranges
    const f0 = (g.start / 3) | 0, f1 = Math.min(F, ((g.start + g.count) / 3) | 0);
    for (let f = f0; f < f1; f++) ids[f] = gi;
  });
  return { ids, names };
}

/**
 * Rebuild contiguous groups from per-face group ids: faces are stably reordered by group.
 * @returns {{ indices: Uint32Array, groups: Array<{ name: string, start: number, count: number }>, order: Uint32Array }}
 *   `order[newFace] = oldFace`; group start/count are in index units.
 */
export function groupsFromIds(indices, ids, names) {
  const F = (indices.length / 3) | 0;
  const order = new Uint32Array(F);
  for (let f = 0; f < F; f++) order[f] = f;
  let sorted = true;
  for (let f = 1; f < F; f++) if (ids[f] < ids[f - 1]) { sorted = false; break; }
  const key = (g) => (g < 0 ? names.length : g);
  if (!sorted) {
    const arr = Array.from(order);
    arr.sort((a, b) => key(ids[a]) - key(ids[b]) || a - b);
    order.set(arr);
  } else if (F && ids[0] < 0 && ids[F - 1] >= 0) {
    const arr = Array.from(order);
    arr.sort((a, b) => key(ids[a]) - key(ids[b]) || a - b);
    order.set(arr);
  }
  const out = new Uint32Array(indices.length);
  const groups = [];
  let cur = -2;
  for (let n = 0; n < F; n++) {
    const f = order[n];
    out[3 * n] = indices[3 * f]; out[3 * n + 1] = indices[3 * f + 1]; out[3 * n + 2] = indices[3 * f + 2];
    const g = ids[f];
    if (g >= 0) {
      if (g !== cur) groups.push({ name: names[g], start: 3 * n, count: 0 });
      groups[groups.length - 1].count += 3;
    }
    cur = g;
  }
  return { indices: out, groups, order };
}
