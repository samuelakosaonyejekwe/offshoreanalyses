// Shared fixture builders and assertions for the geometry importer tests.
import assert from 'node:assert/strict';

export const enc = (s) => new TextEncoder().encode(s);

export function close(a, b, tol = 1e-9, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

export function closeArr(a, b, tol = 1e-9, msg = '') {
  assert.equal(a.length, b.length, `${msg} length`);
  for (let i = 0; i < a.length; i++) close(a[i], b[i], tol, `${msg}[${i}]`);
}

/** Unit cube [0,1]³ as 8 vertices and 12 outward-oriented triangles. */
export const CUBE_V = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
export const CUBE_Q = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
export const CUBE_T = CUBE_Q.flatMap((q) => [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]);

/** Enclosed volume of an indexed triangle surface by the divergence theorem (positive when outward). */
export function surfaceVolume(m) {
  const p = m.positions, t = m.indices;
  let v = 0;
  for (let i = 0; i < t.length; i += 3) {
    const a = 3 * t[i], b = 3 * t[i + 1], c = 3 * t[i + 2];
    v += p[a] * (p[b + 1] * p[c + 2] - p[b + 2] * p[c + 1]) - p[a + 1] * (p[b] * p[c + 2] - p[b + 2] * p[c]) + p[a + 2] * (p[b] * p[c + 1] - p[b + 1] * p[c]);
  }
  return v / 6;
}

/** Total area of an indexed triangle surface. */
export function surfaceArea(m) {
  const p = m.positions, t = m.indices;
  let s = 0;
  for (let i = 0; i < t.length; i += 3) {
    const a = 3 * t[i], b = 3 * t[i + 1], c = 3 * t[i + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2];
    const vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    s += Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx) / 2;
  }
  return s;
}

export function binaryStl(tris, verts, header = '') {
  const out = new Uint8Array(84 + 50 * tris.length);
  out.set(enc(header).subarray(0, 80));
  const dv = new DataView(out.buffer);
  dv.setUint32(80, tris.length, true);
  tris.forEach((t, i) => t.forEach((k, j) => verts[k].forEach((c, a) => dv.setFloat32(84 + 50 * i + 12 + 12 * j + 4 * a, c, true))));
  return out;
}

export function concat(...parts) {
  const n = parts.reduce((s, p) => s + p.length, 0);
  const out = new Uint8Array(n);
  let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}

/** Every importer must return these. */
export function assertContract(r, kinds) {
  assert.ok(Array.isArray(r.models) && r.models.length > 0, 'models');
  assert.ok(Array.isArray(r.warnings), 'warnings');
  assert.ok(r.meta && r.meta.counts && r.meta.bbox, 'meta.counts and meta.bbox');
  assert.ok(r.meta.elementTypes && typeof r.meta.elementTypes === 'object', 'meta.elementTypes');
  assert.ok(r.meta.skipped && typeof r.meta.skipped === 'object', 'meta.skipped');
  assert.ok('units' in r.meta && 'crs' in r.meta, 'meta.units and meta.crs');
  for (const v of [...r.meta.bbox.min, ...r.meta.bbox.max]) assert.ok(Number.isFinite(v), 'finite bbox');
  if (kinds) for (const m of r.models) assert.ok(kinds.includes(m.kind), `unexpected model kind ${m.kind}`);
  for (const m of r.models) {
    assert.equal(typeof m.name, 'string');
    if (m.kind === 'surface') { assert.ok(m.positions instanceof Float32Array); assert.ok(m.indices instanceof Uint32Array); assert.equal(m.indices.length % 3, 0); const n = m.positions.length / 3; for (const i of m.indices) assert.ok(i < n, 'index in range'); }
    if (m.kind === 'points') assert.ok(m.positions instanceof Float32Array);
    if (m.kind === 'polyline') assert.ok(m.points instanceof Float64Array);
    if (m.kind === 'frame') { assert.ok(m.nodes instanceof Float64Array); const n = m.nodes.length / 3; for (const mb of m.members) assert.ok(mb.a < n && mb.b < n && mb.a !== mb.b, 'member ends'); }
    if (m.kind === 'grid') { assert.ok(m.z instanceof Float32Array); assert.equal(m.z.length, m.nx * m.ny); assert.ok(m.dx > 0 && m.dy > 0); }
    if (m.kind === 'section') for (const l of m.loops) assert.ok(l instanceof Float64Array && l.length >= 6);
    if (m.kind === 'volume') for (const c of m.cells) assert.ok(c.conn instanceof Uint32Array);
  }
}

// ---- finite-element reference shapes (canonical corner order: see ELEM in importers/common.js) ----------

export const REF = {
  tet: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
  hex: [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],
  wedge: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]],
  pyr: [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]],
  tri: [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
  quad: [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
};

/** Canonical mid-side edges (pairs of corner indices) in canonical node order. */
export const CANON_EDGES = {
  tri6: [[0, 1], [1, 2], [2, 0]],
  quad8: [[0, 1], [1, 2], [2, 3], [3, 0]],
  tet10: [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]],
  hex20: [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]],
};

/** Corner coordinates shifted by `d`, followed by the midpoints of `edges` in the order given. */
export function withMids(corners, edges, d = [0, 0, 0]) {
  const c = corners.map((p) => [p[0] + d[0], p[1] + d[1], p[2] + d[2]]);
  return c.concat(edges.map(([a, b]) => [(c[a][0] + c[b][0]) / 2, (c[a][1] + c[b][1]) / 2, (c[a][2] + c[b][2]) / 2]));
}

const APEX = { tet4: [0, 1, 2, 3], tet10: [0, 1, 2, 3], hex8: [0, 1, 3, 4], hex20: [0, 1, 3, 4], wedge6: [0, 1, 2, 3], pyr5: [0, 1, 3, 4] };
const NODES = { tet4: 4, tet10: 10, hex8: 8, hex20: 20, wedge6: 6, pyr5: 5, tri3: 3, tri6: 6, quad4: 4, quad8: 8, quad9: 9 };

/**
 * Check every cell/shell block of a model: corner orientation gives a positive volume (right-handed canonical
 * ordering) and every mid-side node sits at the midpoint of the canonical corner pair.
 */
export function checkElements(model) {
  const pos = model.positions;
  const o = model.origin || [0, 0, 0];
  const P = (i) => [pos[3 * i] + o[0], pos[3 * i + 1] + o[1], pos[3 * i + 2] + o[2]];
  const blocks = model.kind === 'volume' ? model.cells : model.elements || [];
  for (const b of blocks) {
    const n = NODES[b.type];
    assert.ok(n, `known block type ${b.type}`);
    assert.equal(b.conn.length % n, 0);
    for (let e = 0; e < b.conn.length / n; e++) {
      const c = Array.from(b.conn.subarray(e * n, (e + 1) * n), P);
      if (APEX[b.type]) {
        const [a, bb, cc, d] = APEX[b.type].map((k) => c[k]);
        const u = [bb[0] - a[0], bb[1] - a[1], bb[2] - a[2]], v = [cc[0] - a[0], cc[1] - a[1], cc[2] - a[2]], w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        const vol = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
        assert.ok(vol > 0, `${b.type} element ${e} has positive orientation (got ${vol})`);
      }
      const edges = CANON_EDGES[b.type];
      if (edges) edges.forEach(([i, j], k) => { for (let a = 0; a < 3; a++) close(c[b.type === 'tri6' ? 3 + k : b.type === 'quad8' ? 4 + k : b.type === 'tet10' ? 4 + k : 8 + k][a], (c[i][a] + c[j][a]) / 2, 1e-6, `${b.type} mid-side node ${k}`); });
    }
  }
}
