import test from 'node:test';
import assert from 'node:assert/strict';
import {
  bounds, weldVertices, weld, buildEdges, boundaryLoops, nonManifoldEdges, nonManifoldVertices, connectedComponents,
  duplicateFaces, degenerateFaces, orientConsistently, eulerCharacteristic, boundaryFaces, cellVolume, toSurface,
  selfIntersections, inspect, watertight, faceGroupIds, groupsFromIds,
} from '../../src/mesh/topology.js';
import { box, sphere, cylinder, monopile } from '../../src/mesh/generate.js';
import { close, hexBlock } from './helpers.js';

const soup = (m) => {
  // STL-style: every triangle has its own three vertices
  const P = new Float32Array(m.indices.length * 3);
  for (let i = 0; i < m.indices.length; i++) P.set(m.positions.subarray(3 * m.indices[i], 3 * m.indices[i] + 3), 3 * i);
  return { kind: 'surface', name: 'soup', positions: P, indices: Uint32Array.from({ length: m.indices.length }, (_, i) => i) };
};
const codes = (r) => r.issues.map((i) => i.code);

test('vertex welding: exact, tolerance-based, and the old → new map', () => {
  const P = Float64Array.from([0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1e-9, 0, 5, 5, 5]);
  const exact = weldVertices(P, 0);
  assert.equal(exact.count, 4); assert.equal(exact.merged, 1);
  assert.deepEqual(Array.from(exact.map), [0, 1, 0, 2, 3]);
  const tol = weldVertices(P, 1e-6);
  assert.equal(tol.count, 3);
  assert.deepEqual(Array.from(tol.map), [0, 1, 0, 1, 2]);
  assert.deepEqual(Array.from(tol.rep), [0, 1, 4]);
  assert.ok(tol.positions instanceof Float64Array);
  assert.throws(() => weldVertices(P, -1), RangeError);
  // an STL-style triangle soup of a box welds back to its 8 corners
  const s = soup(box({ size: [1, 2, 3] }));
  assert.equal(s.positions.length / 3, 36);
  const w = weld(s);
  assert.equal(w.model.positions.length / 3, 8); assert.equal(w.merged, 28);
  assert.deepEqual(bounds(w.model.positions).max, [0.5, 1, 1.5]);
});

test('edge–face adjacency of a cube: 18 edges, each shared by two faces', () => {
  const m = box({ size: [1, 1, 1] });
  const e = buildEdges(m.indices, 8);
  assert.equal(e.edgeCount, 18);
  for (let i = 0; i < e.edgeCount; i++) {
    assert.equal(e.edgeFaceCount(i), 2);
    assert.ok(e.edgeVerts[2 * i] < e.edgeVerts[2 * i + 1]);
    // the two half-edges run in opposite directions on an oriented surface
    const [h0, h1] = [e.edgeHalfEdges[e.edgeFaceStart[i]], e.edgeHalfEdges[e.edgeFaceStart[i] + 1]];
    assert.equal(m.indices[h0], m.indices[3 * ((h1 / 3) | 0) + (h1 % 3 + 1) % 3]);
  }
  assert.equal(boundaryLoops(m.indices, 8).loops.length, 0);
});

test('Euler characteristic: sphere 2, torus-like hollow pile 0, disc 1, two components 4', () => {
  const s = sphere({ radius: 1, subdivisions: 2 });
  const eu = eulerCharacteristic(s.indices, s.positions.length / 3);
  assert.deepEqual([eu.V, eu.E, eu.F, eu.chi, eu.genus, eu.boundaryLoops], [162, 480, 320, 2, 0, 0]);
  const t = monopile({ diameter: 2, zTop: 5, zToe: 0, wallThickness: 0.2, segments: 12, elementSize: 1 });
  const et = eulerCharacteristic(t.indices, t.positions.length / 3);
  assert.equal(et.chi, 0); assert.equal(et.genus, 1);
  const b = box({ size: [1, 1, 1], divisions: [2, 2, 2] });
  const top = b.groups.find((g) => g.name === 'z+');
  const disc = b.indices.slice(top.start, top.start + top.count);
  const ed = eulerCharacteristic(disc, b.positions.length / 3);
  assert.equal(ed.chi, 1); assert.equal(ed.boundaryLoops, 1); assert.equal(ed.genus, 0);
  const two = Uint32Array.from([...s.indices, ...Array.from(s.indices, (v) => v + 162)]);
  assert.equal(eulerCharacteristic(two, 324).chi, 4);
});

test('boundary loops: a removed face leaves one closed three-edge loop running against the face winding', () => {
  const m = box({ size: [1, 1, 1] });
  const idx = m.indices.slice(3);
  const bl = boundaryLoops(idx, 8);
  assert.equal(bl.loops.length, 1); assert.equal(bl.edgeCount, 3); assert.deepEqual(bl.closed, [true]);
  assert.deepEqual(bl.loops[0].slice().sort(), Array.from(m.indices.slice(0, 3)).sort());
  // open cylinder side: two loops of n edges
  const c = cylinder({ radius: 1, height: 2, segments: 12 });
  const side = c.groups.find((g) => g.name === 'side');
  const b2 = boundaryLoops(c.indices.slice(side.start, side.start + side.count), c.positions.length / 3);
  assert.equal(b2.loops.length, 2);
  assert.deepEqual(b2.loops.map((l) => l.length), [12, 12]);
});

test('non-manifold edges and bow-tie vertices are found', () => {
  // three triangles on one edge (a fin on a book spine)
  const fin = Uint32Array.from([0, 1, 2, 0, 1, 3, 0, 1, 4]);
  assert.deepEqual(nonManifoldEdges(fin, 5), [[0, 1]]);
  // two triangles touching at a single vertex
  const bow = Uint32Array.from([0, 1, 2, 2, 3, 4]);
  assert.deepEqual(nonManifoldVertices(bow, 5), [2]);
  assert.deepEqual(nonManifoldEdges(bow, 5), []);
  assert.deepEqual(nonManifoldVertices(box({ size: [1, 1, 1] }).indices, 8), []);
});

test('connected components by edge and by vertex', () => {
  const bow = Uint32Array.from([0, 1, 2, 2, 3, 4, 5, 6, 7]);
  assert.equal(connectedComponents(bow, 8, 'edge').count, 3);
  const cv = connectedComponents(bow, 8, 'vertex');
  assert.equal(cv.count, 2);
  assert.deepEqual(cv.sizes, [2, 1]);
  assert.deepEqual(Array.from(cv.faceComponent), [0, 0, 1]);
});

test('degenerate and duplicate faces', () => {
  const P = Float64Array.from([0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0]);
  const idx = Uint32Array.from([0, 1, 3, 0, 1, 2, 1, 1, 3, 3, 1, 0, 0, 1, 3]);
  assert.deepEqual(degenerateFaces(P, idx), [1, 2]); // collinear, repeated vertex
  assert.deepEqual(duplicateFaces(idx, 4), [3, 4]); // same vertices in any order
});

test('orientation: scrambled windings become consistent and outward; volume sign follows', () => {
  const m = sphere({ radius: 1, subdivisions: 2 });
  const idx = Uint32Array.from(m.indices);
  for (const f of [0, 5, 17, 100, 250]) { const t = idx[3 * f + 1]; idx[3 * f + 1] = idx[3 * f + 2]; idx[3 * f + 2] = t; }
  const o = orientConsistently(m.positions, idx);
  assert.equal(o.flipped, 5); assert.ok(o.orientable); assert.equal(o.components, 1);
  assert.deepEqual(Array.from(o.indices), Array.from(m.indices));
  // a wholly inside-out sphere is turned outward
  const inv = Uint32Array.from(m.indices);
  for (let f = 0; f < inv.length; f += 3) { const t = inv[f + 1]; inv[f + 1] = inv[f + 2]; inv[f + 2] = t; }
  const o2 = orientConsistently(m.positions, inv);
  assert.equal(o2.flipped, 320);
  assert.ok(o2.componentVolume[0] > 4);
  // Möbius strip: no consistent orientation exists
  const n = 8;
  const mob = [];
  for (let i = 0; i < n; i++) { const a = 2 * i, b = 2 * i + 1, c = i === n - 1 ? 1 : 2 * i + 2, d = i === n - 1 ? 0 : 2 * i + 3; mob.push(a, b, c, b, d, c); }
  const P = new Float64Array(6 * n);
  for (let i = 0; i < n; i++) { const t = 2 * Math.PI * i / n; for (const s of [0, 1]) { const w = (s - 0.5) * 0.4; P.set([(1 + w * Math.cos(t / 2)) * Math.cos(t), (1 + w * Math.cos(t / 2)) * Math.sin(t), w * Math.sin(t / 2)], 3 * (2 * i + s)); } }
  assert.equal(orientConsistently(P, Uint32Array.from(mob)).orientable, false);
});

test('volume meshes: boundary faces, cell volumes and inspection of a hexahedral block', () => {
  const m = hexBlock(3, 2);
  const bf = boundaryFaces(m);
  assert.equal(bf.quads.length / 4, 54); assert.equal(bf.tris.length, 0);
  assert.equal(bf.interiorFaces, 3 * 9 * 2); assert.equal(bf.overSharedFaces, 0);
  close(cellVolume(m.positions, m.cells[0].conn, 0, 'hex8'), 8 / 27, 1e-12);
  const s = toSurface(m);
  assert.equal(s.indices.length / 3, 108); assert.equal(s.faceCell.length, 108);
  const r = inspect(m);
  assert.ok(r.ok && r.watertight);
  close(r.stats.volume, 8, 1e-12);
  assert.deepEqual(r.stats.cellsByType, { hex8: 27 });
  assert.equal(r.stats.boundary.quads, 54);
  assert.ok(watertight(m).watertight);
  // reversed node order → inverted cell; unknown cell type → error
  const tet = { kind: 'volume', positions: Float64Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 9, 9, 9]), cells: [{ type: 'tet4', conn: Uint32Array.from([0, 2, 1, 3]) }] };
  const rt = inspect(tet);
  assert.ok(!rt.ok); assert.ok(codes(rt).includes('inverted-cells')); assert.ok(codes(rt).includes('unreferenced-nodes'));
  close(cellVolume(tet.positions, Uint32Array.from([0, 1, 2, 3]), 0, 'tet4'), 1 / 6, 1e-15);
  assert.ok(codes(inspect({ kind: 'volume', positions: tet.positions, cells: [{ type: 'hex27', conn: new Uint32Array(27) }] })).includes('unsupported-cell-type'));
  // a hanging node: one coarse cell beside two half-height cells leaves free edges on the boundary
  const P = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 2, 0, 0, 2, 1, 0, 2, 0, 0.5, 2, 1, 0.5, 1, 0, 0.5, 1, 1, 0.5, 2, 0, 1, 2, 1, 1];
  const hang = { kind: 'volume', positions: Float64Array.from(P), cells: [{ type: 'hex8', conn: Uint32Array.from([0, 1, 2, 3, 4, 5, 6, 7, 1, 8, 9, 2, 12, 10, 11, 13, 12, 10, 11, 13, 5, 14, 15, 6]) }] };
  const rh = inspect(hang);
  assert.ok(codes(rh).includes('non-conforming-boundary')); assert.equal(rh.watertight, false);
});

test('self-intersection: two interpenetrating boxes are reported, a clean sphere is not', () => {
  // offset so that no edge of one box passes exactly through an edge of the other (touching is not reported)
  const a = box({ size: [2, 2, 2] }), b = box({ size: [2, 2, 2], center: [1, 0.7, 0.4] });
  const P = Float32Array.from([...a.positions, ...b.positions]);
  const idx = Uint32Array.from([...a.indices, ...Array.from(b.indices, (v) => v + 8)]);
  const si = selfIntersections(P, idx);
  assert.ok(si.pairs.length >= 6);
  const r = inspect({ kind: 'surface', positions: P, indices: idx });
  assert.ok(codes(r).includes('self-intersections')); assert.equal(r.watertight, false);
  const s = sphere({ radius: 1, subdivisions: 3 });
  assert.equal(selfIntersections(s.positions, s.indices).pairs.length, 0);
});

test('inspect surface: a triangle soup is welded for the verdict; issues carry codes, counts and messages', () => {
  const r = inspect(soup(sphere({ radius: 2, subdivisions: 2 })));
  assert.ok(r.ok && r.watertight);
  assert.deepEqual(codes(r), ['unwelded-vertices']);
  assert.equal(r.issues[0].count, 960 - 162);
  assert.equal(typeof r.issues[0].message, 'string');
  assert.equal(r.stats.eulerCharacteristic, 2); assert.equal(r.stats.weldedVertices, 162);
  assert.equal(r.kind, 'surface');
  // unwelded inspection sees 320 separate triangles
  const raw = inspect(soup(sphere({ radius: 2, subdivisions: 2 })), { weld: false, selfIntersections: false });
  assert.equal(raw.watertight, false); assert.equal(raw.stats.components, 320);
  assert.ok(!inspect({ kind: 'surface', positions: new Float32Array(9), indices: Uint32Array.from([0, 1, 7]) }).ok);
  assert.ok(codes(inspect({ kind: 'surface', positions: new Float32Array(0), indices: new Uint32Array(0) })).includes('empty'));
  assert.throws(() => inspect(null), RangeError);
  assert.throws(() => inspect({ kind: 'surface', positions: new Float32Array(9) }), RangeError);
});

test('inspect frame: dangling, zero-length and duplicate members, coincident nodes, disconnected parts', () => {
  const nodes = Float64Array.from([0, 0, 0, 0, 0, 10, 10, 0, 10, 10, 0, 0, 10, 0, 1e-9, 50, 0, 0, 50, 0, 5, 99, 99, 99]);
  const members = [{ a: 0, b: 1, section: { type: 'tube', D: 1, t: 0.02 } }, { a: 1, b: 2 }, { a: 2, b: 3 }, { a: 2, b: 1 }, { a: 3, b: 3 }, { a: 3, b: 4 }, { a: 5, b: 6 }, { a: 0, b: 40 }];
  const r = inspect({ kind: 'frame', name: 'f', nodes, members });
  const by = Object.fromEntries(r.issues.map((i) => [i.code, i]));
  assert.equal(by['duplicate-members'].count, 1); assert.deepEqual(by['duplicate-members'].elements, [3]);
  assert.equal(by['zero-length-members'].count, 2); assert.deepEqual(by['zero-length-members'].elements, [4, 5]);
  assert.equal(by['index-out-of-range'].count, 1);
  assert.equal(by['coincident-nodes'].count, 1);
  assert.equal(by['disconnected-parts'].count, 2);
  assert.equal(by['unreferenced-nodes'].count, 2);
  assert.equal(by['free-ends'].count, 4);
  assert.equal(by['members-without-section'].count, 7);
  assert.ok(!r.ok); assert.equal(r.watertight, false);
  close(r.stats.totalLength, 10 + 10 + 10 + 10 + 5, 1e-9);
  const good = inspect({ kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 10, 0, 10]), members: [{ a: 0, b: 1, section: { A: 1 } }, { a: 1, b: 2, section: { A: 1 } }] });
  assert.ok(good.ok); assert.equal(good.stats.components, 1); assert.equal(good.stats.freeEnds, 2);
});

test('inspect polyline, points, grid and section', () => {
  const pl = inspect({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 3, 4, 0, 3, 4, 0, 3, 4, 12]) });
  assert.ok(pl.ok); assert.equal(pl.issues[0].code, 'repeated-points'); close(pl.stats.length, 17, 1e-12);
  assert.ok(!inspect({ kind: 'polyline', points: Float64Array.from([0, 0, 0]) }).ok);
  const pts = inspect({ kind: 'points', positions: Float32Array.from([0, 0, 0, 1, 1, 1, 0, 0, 0]) });
  assert.equal(pts.issues[0].code, 'duplicate-points'); assert.equal(pts.stats.points, 3);
  assert.ok(!inspect({ kind: 'points', positions: Float32Array.from([0, NaN, 0]) }).ok);
  const g = inspect({ kind: 'grid', nx: 2, ny: 2, x0: 0, y0: 0, dx: 1, dy: 1, z: Float32Array.from([-5, NaN, -7, -6]) });
  assert.ok(g.ok); assert.equal(g.stats.noData, 1); assert.equal(g.stats.zMin, -7);
  assert.ok(!inspect({ kind: 'grid', nx: 2, ny: 2, x0: 0, y0: 0, dx: 1, dy: 1, z: new Float32Array(3) }).ok);
  const sec = inspect({ kind: 'section', loops: [Float64Array.from([0, 0, 2, 0, 2, 1, 0, 1])] });
  assert.ok(sec.ok); close(sec.stats.outerArea, 2, 1e-12);
  const bow = inspect({ kind: 'section', loops: [Float64Array.from([0, 0, 2, 2, 2, 0, 0, 2])] });
  assert.ok(codes(bow).includes('self-intersecting-loops'));
  assert.ok(codes(inspect({ kind: 'hologram' })).includes('unknown-kind'));
});

test('face groups survive a round trip through per-face ids', () => {
  const m = box({ size: [1, 1, 1], divisions: [2, 1, 1] });
  const g = faceGroupIds(m);
  assert.equal(g.names.length, 6); assert.ok(g.ids.every((i) => i >= 0));
  // reverse the face order: groups are rebuilt contiguous
  const F = m.indices.length / 3;
  const idx = new Uint32Array(m.indices.length), ids = new Int32Array(F);
  for (let f = 0; f < F; f++) { idx.set(m.indices.subarray(3 * (F - 1 - f), 3 * (F - f)), 3 * f); ids[f] = g.ids[F - 1 - f]; }
  const r = groupsFromIds(idx, ids, g.names);
  assert.deepEqual(r.groups.map((x) => [x.name, x.count]), m.groups.map((x) => [x.name, x.count]));
  assert.deepEqual(Array.from(r.indices).sort(), Array.from(m.indices).sort());
});
