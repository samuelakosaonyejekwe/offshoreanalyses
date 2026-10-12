import test from 'node:test';
import assert from 'node:assert/strict';
import {
  repair, removeDegenerateFaces, removeDuplicateFaces, removeUnreferencedVertices, fillHoles, fixOrientation,
  splitNonManifoldVertices, dropSmallComponents, healFrame,
} from '../../src/mesh/repair.js';
import { inspect } from '../../src/mesh/topology.js';
import { box, sphere, cylinder } from '../../src/mesh/generate.js';
import { close, rel, hexBlock } from './helpers.js';

const codes = (r) => r.issues.map((i) => i.code);
const acts = (r) => Object.fromEntries(r.actions.map((a) => [a.code, a.count]));

/** A sphere damaged the way exported CAD and scan meshes usually are. */
function broken() {
  const s = sphere({ radius: 2, subdivisions: 3 });
  const F = s.indices.length / 3;
  const P = [], I = [];
  // triangle soup (no shared vertices), with two faces missing, some reversed
  const skip = new Set([10, 400]);
  const flip = new Set([3, 50, 51, 52, 700]);
  for (let f = 0; f < F; f++) {
    if (skip.has(f)) continue;
    const v = [s.indices[3 * f], s.indices[3 * f + 1], s.indices[3 * f + 2]];
    if (flip.has(f)) v.reverse();
    for (const q of v) { I.push(P.length / 3); P.push(s.positions[3 * q], s.positions[3 * q + 1], s.positions[3 * q + 2]); }
  }
  const tri = (a, b, c) => { for (const q of [a, b, c]) { I.push(P.length / 3); P.push(...q); } };
  const p0 = Array.from(s.positions.subarray(3 * s.indices[0], 3 * s.indices[0] + 3)), p1 = Array.from(s.positions.subarray(3 * s.indices[1], 3 * s.indices[1] + 3));
  // an exact copy of face 0, a collapsed face on one of its edges, a vertex nothing uses, and a speck of debris far away
  tri(...[0, 1, 2].map((k) => Array.from(s.positions.subarray(3 * s.indices[k], 3 * s.indices[k] + 3))));
  tri(p0, p1, p1);
  P.push(9, 9, 9);
  tri([20, 0, 0], [20.01, 0, 0], [20, 0.01, 0]);
  return { model: { kind: 'surface', name: 'damaged', positions: Float32Array.from(P), indices: Uint32Array.from(I), groups: [{ name: 'hull', start: 0, count: 3 * (F - 2) }] }, exact: s };
}

test('a deliberately broken mesh is diagnosed', () => {
  const { model } = broken();
  const r = inspect(model);
  assert.equal(r.ok, false); assert.equal(r.watertight, false);
  const by = Object.fromEntries(r.issues.map((i) => [i.code, i.count]));
  assert.ok(by['unwelded-vertices'] > 3000);
  assert.equal(by['duplicate-faces'], 1);
  assert.equal(by['degenerate-faces'], 1);
  assert.equal(by['unreferenced-vertices'], 1);
  assert.equal(by['open-boundary'], 6 + 3); // two triangular holes and the free debris triangle
  assert.equal(by['non-manifold-edges'], 3); // the edges of the doubled face
  assert.equal(by['multiple-components'], 2);
  assert.ok(by['inconsistent-orientation'] >= 5);
  assert.ok(r.reasons.some((t) => t.startsWith('open')));
  assert.equal(r.issues.find((i) => i.code === 'open-boundary').edges.length, 9);
});

test('…and repaired: watertight, outward, analytic volume restored, every action reported', () => {
  const { model, exact } = broken();
  const before = Float32Array.from(model.positions);
  const fixed = repair(model);
  assert.deepEqual(Array.from(model.positions), Array.from(before), 'input untouched');
  const a = acts(fixed);
  assert.ok(a['welded-vertices'] > 3000);
  assert.equal(a['removed-degenerate-faces'], 1);
  assert.equal(a['removed-duplicate-faces'], 1);
  assert.equal(a['dropped-small-components'], 1);
  assert.equal(a['filled-holes'], 2);
  assert.equal(a['flipped-faces'], 5);
  assert.equal(a['removed-unreferenced-vertices'], 4); // the stray vertex and the three debris corners
  assert.ok(fixed.actions.every((x) => typeof x.message === 'string' && x.count > 0));
  const r = inspect(fixed.model);
  assert.ok(r.ok, codes(r).join(','));
  assert.ok(r.watertight); assert.ok(r.stats.outward);
  assert.deepEqual(codes(r), []);
  assert.equal(r.stats.eulerCharacteristic, 2);
  assert.equal(fixed.model.indices.length / 3, exact.indices.length / 3);
  assert.equal(fixed.model.positions.length, exact.positions.length);
  rel(r.stats.volume, inspect(exact).stats.volume, 1e-6);
  // groups survive: the original hull plus the patches
  assert.deepEqual(fixed.model.groups.map((g) => g.name), ['hull', 'filled-holes']);
  assert.equal(fixed.model.groups[1].count, 6);
  // repairing a clean mesh does nothing
  assert.deepEqual(repair(exact).actions, []);
});

test('individual surface repairs', () => {
  const b = box({ size: [1, 1, 1] });
  const dup = { ...b, indices: Uint32Array.from([...b.indices, b.indices[1], b.indices[2], b.indices[0]]) };
  assert.equal(removeDuplicateFaces(dup).removed, 1);
  const deg = { ...b, indices: Uint32Array.from([...b.indices, 0, 0, 1]) };
  const rd = removeDegenerateFaces(deg);
  assert.equal(rd.removed, 1); assert.equal(rd.model.indices.length, 36);
  assert.deepEqual(rd.model.groups.map((g) => g.count), b.groups.map((g) => g.count));
  const extra = { ...b, positions: Float32Array.from([...b.positions, 7, 7, 7]) };
  const ru = removeUnreferencedVertices(extra);
  assert.equal(ru.removed, 1); assert.equal(ru.map[8], 0xffffffff);
  const flipped = { ...b, indices: Uint32Array.from(b.indices) };
  [flipped.indices[1], flipped.indices[2]] = [flipped.indices[2], flipped.indices[1]];
  const fo = fixOrientation(flipped);
  assert.equal(fo.flipped, 1); assert.deepEqual(Array.from(fo.model.indices), Array.from(b.indices));
  // two cubes sharing one corner vertex: the bow-tie is separated
  const c2 = box({ size: [1, 1, 1], center: [1, 1, 1] });
  const corner = Array.from({ length: 8 }, (_, v) => Array.from(c2.positions.subarray(3 * v, 3 * v + 3)));
  const shared = corner.findIndex((p) => p.every((x) => x === 0.5));
  const P = Float32Array.from([...b.positions, ...c2.positions]);
  let target = -1;
  for (let v = 0; v < 8; v++) if (b.positions[3 * v] === 0.5 && b.positions[3 * v + 1] === 0.5 && b.positions[3 * v + 2] === 0.5) target = v;
  const idx = Uint32Array.from([...b.indices, ...Array.from(c2.indices, (v) => (v === shared ? target : v + 8))]);
  const tie = { kind: 'surface', positions: P, indices: idx };
  assert.ok(codes(inspect(tie, { weld: false })).includes('non-manifold-vertices'));
  const sp = splitNonManifoldVertices(tie);
  assert.equal(sp.split, 1);
  assert.ok(!codes(inspect(sp.model, { weld: false })).includes('non-manifold-vertices'));
});

test('hole filling: small holes are closed to the right volume, large ones are left and reported', () => {
  const c = cylinder({ radius: 1, height: 2, segments: 24, capRings: 1 });
  const top = c.groups.find((g) => g.name === 'top');
  const open = { ...c, indices: Uint32Array.from([...c.indices.subarray(0, top.start), ...c.indices.subarray(top.start + top.count)]), groups: c.groups.filter((g) => g !== top) };
  assert.equal(inspect(open).watertight, false);
  const f = fillHoles(open);
  assert.equal(f.filled, 1); assert.equal(f.addedFaces, 24);
  const r = inspect(fixOrientation(f.model).model);
  assert.ok(r.watertight);
  rel(r.stats.volume, inspect(c).stats.volume, 1e-6);
  assert.equal(f.model.groups[f.model.groups.length - 1].name, 'filled-holes');
  const left = repair(open, { maxHoleEdges: 8 });
  assert.equal(acts(left)['holes-left-open'], 1);
  assert.equal(inspect(left.model).watertight, false);
  // a quadrilateral hole takes two triangles
  const b = box({ size: [1, 1, 1] });
  const zp = b.groups.find((g) => g.name === 'z+');
  const holed = { kind: 'surface', positions: b.positions, indices: Uint32Array.from([...b.indices.subarray(0, zp.start), ...b.indices.subarray(zp.start + zp.count)]) };
  const fq = fillHoles(holed);
  assert.equal(fq.addedFaces, 2);
  assert.ok(inspect(fq.model).watertight);
  close(inspect(fq.model).stats.volume, 1, 1e-6);
});

test('small components are dropped by face count or by area share', () => {
  const s = sphere({ radius: 1, subdivisions: 2 });
  const d = sphere({ radius: 0.001, subdivisions: 0, center: [5, 0, 0] });
  const n = s.positions.length / 3;
  const both = { kind: 'surface', positions: Float32Array.from([...s.positions, ...d.positions]), indices: Uint32Array.from([...s.indices, ...Array.from(d.indices, (v) => v + n)]) };
  const r = dropSmallComponents(both);
  assert.equal(r.removedComponents, 1); assert.equal(r.removedFaces, 20);
  assert.equal(dropSmallComponents(both, { minAreaFraction: 1e-9 }).removedComponents, 0);
  assert.equal(dropSmallComponents(s).removedComponents, 0);
});

test('frame healing: near-coincident nodes merged, members split at intermediate nodes, duplicates removed', () => {
  // a portal whose beam passes over an unconnected mid-column, with a doubled node and a repeated member
  const nodes = Float64Array.from([
    0, 0, 0, 0, 0, 10, 20, 0, 10, 20, 0, 0, // 0–3 portal
    10, 0, 0, 10, 0, 10, // 4, 5 mid column; node 5 lies on the beam 1–2
    0, 0, 10.0004, // 6 duplicates node 1 within tolerance
    50, 50, 50, // 7 unused
  ]);
  const sec = { type: 'tube', D: 0.5, t: 0.02 };
  const members = [{ a: 0, b: 1, section: sec, group: 'legs' }, { a: 6, b: 2, section: sec, group: 'beam' }, { a: 2, b: 3, section: sec, group: 'legs' }, { a: 4, b: 5, section: sec, group: 'legs' }, { a: 2, b: 1, section: sec, group: 'beam' }, { a: 3, b: 3 }];
  const frame = { kind: 'frame', name: 'portal', nodes, members, sets: { nodes: { supports: Uint32Array.from([0, 3, 4]) }, members: { beam: Uint32Array.from([1]) } } };
  const before = inspect(frame, { tolerance: 1e-3 });
  assert.ok(codes(before).includes('coincident-nodes') && codes(before).includes('disconnected-parts') && codes(before).includes('zero-length-members'));
  const h = repair(frame);
  const a = acts(h);
  assert.equal(a['merged-nodes'], 1);
  assert.equal(a['split-members'], 2); // both copies of the beam are cut at the column head
  assert.equal(a['removed-zero-length'], 1);
  assert.equal(a['removed-duplicates'], 2);
  assert.equal(a['removed-unused-nodes'], 1);
  assert.equal(h.model.nodes.length / 3, 6);
  assert.equal(h.model.members.length, 5);
  assert.deepEqual(Array.from(h.model.memberOrigin), [0, 1, 1, 2, 3]);
  assert.ok(h.model.members.every((m) => m.section === sec));
  assert.equal(h.flags.isolatedParts, 1);
  assert.equal(h.flags.freeEndNodes.length, 3);
  assert.deepEqual(h.flags.overlappingMembers, [[1, 4], [1, 4]]);
  assert.deepEqual(Array.from(h.model.sets.members.beam), [1, 2]);
  assert.equal(h.model.sets.nodes.supports.length, 3);
  const after = inspect(h.model);
  assert.ok(after.ok); assert.equal(after.stats.components, 1);
  close(after.stats.totalLength, 10 + 20 + 10 + 10, 1e-3);
  assert.equal(frame.members.length, 6, 'input untouched');
  // splitting can be switched off
  assert.equal(healFrame(frame, { split: false }).model.members.length, 4);
});

test('volume repair: inverted cells reordered, zero-volume cells and unused nodes removed, sets remapped', () => {
  const m = hexBlock(2);
  const conn = Uint32Array.from(m.cells[0].conn);
  // cell 3 turned inside out (bottom and top faces swapped)
  for (let k = 0; k < 4; k++) { const t = conn[24 + k]; conn[24 + k] = conn[28 + k]; conn[28 + k] = t; }
  const P = Float64Array.from([...m.positions, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);
  const bad = { kind: 'volume', positions: P, cells: [{ type: 'hex8', conn }, { type: 'tet4', conn: Uint32Array.from([27, 28, 29, 30]) }], sets: { nodes: { corner: Uint32Array.from([0, 26, 30]) }, cells: { last: Uint32Array.from([7, 8]) } } };
  const r0 = inspect(bad);
  assert.ok(codes(r0).includes('inverted-cells') && codes(r0).includes('zero-volume-cells'));
  const fixed = repair(bad);
  assert.deepEqual(acts(fixed), { 'reordered-inverted-cells': 1, 'removed-zero-volume-cells': 1, 'removed-unreferenced-nodes': 4 });
  const r1 = inspect(fixed.model);
  assert.ok(r1.ok && r1.watertight);
  close(r1.stats.volume, 1, 1e-12);
  assert.equal(fixed.model.cells.length, 1);
  assert.deepEqual(Array.from(fixed.model.sets.nodes.corner), [0, 26]);
  assert.deepEqual(Array.from(fixed.model.sets.cells.last), [7]);
  const tet = repair({ kind: 'volume', positions: Float64Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]), cells: [{ type: 'tet4', conn: Uint32Array.from([0, 2, 1, 3]) }] });
  assert.ok(inspect(tet.model).ok);
});

test('polyline, points, grid and section clean-up', () => {
  const pl = repair({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 0, 0, 0, 1, 0, 0, NaN, 0, 0, 2, 0, 0]), attrs: { kp: Float64Array.from([0, 0, 1, 1.5, 2]) } });
  assert.equal(pl.model.points.length, 9); assert.deepEqual(Array.from(pl.model.attrs.kp), [0, 1, 2]);
  const pts = repair({ kind: 'points', positions: Float32Array.from([0, 0, 0, 1, 1, 1, 0, 0, 0, Infinity, 0, 0]) });
  assert.equal(pts.model.positions.length, 6); assert.equal(pts.actions[0].count, 2);
  const g = repair({ kind: 'grid', nx: 3, ny: 3, x0: 0, y0: 0, dx: 1, dy: 1, z: Float32Array.from([1, 1, 1, 1, NaN, 1, 1, 1, 1]) });
  assert.equal(g.model.z[4], 1); assert.equal(g.actions[0].code, 'filled-no-data');
  const sec = repair({ kind: 'section', loops: [Float64Array.from([0, 0, 0, 1, 1, 1, 1, 0, 0, 0]), Float64Array.from([0.2, 0.2, 0.8, 0.2, 0.8, 0.8, 0.2, 0.8]), Float64Array.from([0, 0, 1, 1])] });
  assert.equal(sec.model.loops.length, 2);
  assert.deepEqual(Array.from(sec.model.loops[0]), [1, 0, 1, 1, 0, 1, 0, 0]); // closing point dropped, made counter-clockwise
  assert.deepEqual(Array.from(sec.model.loops[1]), [0.2, 0.8, 0.8, 0.8, 0.8, 0.2, 0.2, 0.2]); // hole made clockwise
  assert.throws(() => repair({ kind: 'hologram' }), RangeError);
  assert.throws(() => repair(null), RangeError);
});
