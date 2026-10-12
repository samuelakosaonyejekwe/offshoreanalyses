import test from 'node:test';
import assert from 'node:assert/strict';
import { refineTriangles, refineQuads, refineVolume, subdivideMembers, resamplePolyline, resampleGrid, refine, refineTo, decimate, coarsenTo, smooth } from '../../src/mesh/refine.js';
import { inspect, cellVolume, CELL_INFO } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { assess } from '../../src/mesh/quality.js';
import { box, sphere, cylinder, wigleyHull } from '../../src/mesh/generate.js';
import { orderFromErrors } from '../../src/mesh/convergence.js';
import { makeRng } from '../../src/engine/math/rng.js';
import { close, rel, hexBlock } from './helpers.js';

const volumes = (m) => { const out = []; for (const b of m.cells) { const nn = CELL_INFO[b.type].nodes; for (let o = 0; o < b.conn.length; o += nn) out.push(cellVolume(m.positions, b.conn, o, b.type)); } return out; };
const sum = (a) => a.reduce((s, v) => s + v, 0);
const single = (type, pts, conn) => ({ kind: 'volume', positions: Float64Array.from(pts.flat()), cells: [{ type, conn: Uint32Array.from(conn) }] });
const S3 = Math.sqrt(3);

test('triangle 1 → 4: counts V + E, 4F; area unchanged; closed surface stays closed', () => {
  const b = box({ size: [1, 2, 3] });
  const r = refineTriangles(b.positions, b.indices);
  assert.equal(r.indices.length / 3, 48);
  assert.equal(r.positions.length / 3, 8 + 18);
  rel(surfaceArea(r.positions, r.indices), 22, 1e-6);
  const m = refine(refine(b));
  assert.equal(m.indices.length / 3, 12 * 16);
  const i = inspect(m);
  assert.ok(i.watertight && i.stats.outward); assert.equal(i.stats.eulerCharacteristic, 2);
  rel(i.stats.volume, 6, 1e-6);
  // groups scale with their faces and stay contiguous
  assert.deepEqual(m.groups.map((g) => [g.name, g.count / 3]), b.groups.map((g) => [g.name, 16 * g.count / 3]));
  assert.equal(m.groups.reduce((s, g) => s + g.count, 0), m.indices.length);
  assert.equal(b.indices.length, 36, 'input untouched');
});

test('quadrilateral 1 → 4: 4 children, shared mid-side nodes, area preserved', () => {
  const P = Float64Array.from([0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0]);
  const r = refineQuads(P, Uint32Array.from([0, 1, 4, 3, 1, 2, 5, 4]));
  assert.equal(r.quads.length / 4, 8);
  assert.equal(r.positions.length / 3, 6 + 7 + 2); // 7 edges, 2 centres
  let area = 0;
  for (let q = 0; q < r.quads.length; q += 4) {
    const [a, b, c, d] = r.quads.subarray(q, q + 4);
    const x = (i) => r.positions[3 * i], y = (i) => r.positions[3 * i + 1];
    const A = 0.5 * ((x(c) - x(a)) * (y(d) - y(b)) - (x(d) - x(b)) * (y(c) - y(a)));
    close(A, 0.25, 1e-12); area += A;
  }
  close(area, 2, 1e-12);
});

test('tetrahedron 1 → 8 (Bey): eight positive children of volume V/8 each, conforming after two levels', () => {
  const m = single('tet4', [[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]], [0, 1, 2, 3]);
  const r = refineVolume(m);
  assert.equal(r.cells.length, 1); assert.equal(r.cells[0].type, 'tet4');
  const v = volumes(r);
  assert.equal(v.length, 8); assert.equal(r.positions.length / 3, 10);
  v.forEach((x) => close(x, 8 / 3 / 8, 1e-12));
  assert.deepEqual(Array.from(r.cells[0].parent), new Array(8).fill(0));
  const r2 = refineVolume(r);
  const i2 = inspect(r2);
  assert.equal(i2.stats.cells, 64); assert.ok(i2.ok && i2.watertight);
  close(i2.stats.volume, 8 / 3, 1e-12);
  assert.equal(i2.stats.boundary.triangles, 4 * 16);
  // Bey's refinement keeps the children in a bounded number of shape classes
  assert.ok(assess(r2).summary.scaledJacobian.min > 0.4);
});

test('hexahedron 1 → 8: 27 nodes from one cell; a block refines to a conforming block', () => {
  const one = hexBlock(1);
  const r = refineVolume(one);
  assert.equal(r.positions.length / 3, 27);
  const v = volumes(r);
  assert.equal(v.length, 8);
  v.forEach((x) => close(x, 1 / 8, 1e-12));
  const blk = refine(hexBlock(2));
  const i = inspect(blk);
  assert.equal(i.stats.cells, 64); assert.equal(i.stats.nodes, 125);
  assert.ok(i.ok && i.watertight);
  assert.equal(i.stats.boundary.quads, 6 * 16); assert.equal(i.stats.boundary.interiorFaces, 3 * 16 * 3);
  close(i.stats.volume, 1, 1e-12);
  assert.equal(assess(blk).grade, 'A');
  // quadratic cells are refined on their corners
  const h20 = { kind: 'volume', positions: one.positions, cells: [{ type: 'hex20', conn: Uint32Array.from([...one.cells[0].conn, ...new Array(12).fill(0)]) }] };
  assert.equal(refineVolume(h20).cells[0].type, 'hex8');
});

test('wedge 1 → 8 and pyramid 1 → 6 pyramids + 4 tetrahedra: positive children, volume conserved', () => {
  const w = refineVolume(single('wedge6', [[0, 0, 0], [1, 0, 0], [0.5, S3 / 2, 0], [0, 0, 1], [1, 0, 1], [0.5, S3 / 2, 1]], [0, 1, 2, 3, 4, 5]));
  const vw = volumes(w);
  assert.equal(vw.length, 8);
  vw.forEach((x) => close(x, S3 / 4 / 8, 1e-12));
  assert.ok(inspect(refineVolume(w)).watertight);
  const py = single('pyr5', [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]], [0, 1, 2, 3, 4]);
  const p = refineVolume(py);
  assert.deepEqual(p.cells.map((b) => [b.type, b.conn.length / CELL_INFO[b.type].nodes]), [['pyr5', 6], ['tet4', 4]]);
  const vp = volumes(p);
  assert.ok(vp.every((x) => x > 0));
  close(sum(vp), 1 / 3, 1e-12);
  vp.slice(0, 6).forEach((x) => close(x, 1 / 24, 1e-12));
  vp.slice(6).forEach((x) => close(x, 1 / 48, 1e-12));
  const p2 = inspect(refineVolume(p));
  assert.ok(p2.ok && p2.watertight); close(p2.stats.volume, 1 / 3, 1e-12);
  assert.throws(() => refineVolume({ kind: 'volume', positions: py.positions, cells: [{ type: 'hex27', conn: new Uint32Array(27) }] }), RangeError);
});

test('mixed mesh stays conforming: a pyramid capping a hexahedron, cell sets follow their children', () => {
  const P = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0.5, 0.5, 2]];
  const m = { kind: 'volume', positions: Float64Array.from(P.flat()), cells: [{ type: 'hex8', conn: Uint32Array.from([0, 1, 2, 3, 4, 5, 6, 7]) }, { type: 'pyr5', conn: Uint32Array.from([4, 5, 6, 7, 8]) }], sets: { nodes: { apex: Uint32Array.from([8]) }, cells: { cap: Uint32Array.from([1]) } } };
  const r = refineVolume(m);
  const i = inspect(r);
  assert.ok(i.ok && i.watertight, i.issues.map((x) => x.code).join(','));
  close(i.stats.volume, 1 + 1 / 3, 1e-12);
  assert.deepEqual(i.stats.cellsByType, { hex8: 8, pyr5: 6, tet4: 4 });
  assert.equal(r.sets.cells.cap.length, 10);
  assert.deepEqual(Array.from(r.sets.nodes.apex), [8]);
  close(r.positions[3 * 8 + 2], 2, 0);
});

test('refinement with projection: sphere area and volume errors fall by four per level', () => {
  const R = 1.5;
  const onSphere = (x, y, z) => { const l = Math.hypot(x, y, z); return [R * x / l, R * y / l, R * z / l]; };
  let m = { ...sphere({ radius: R, subdivisions: 1 }), positions: Float64Array.from(sphere({ radius: R, subdivisions: 1 }).positions) };
  const hs = [], eV = [], eA = [];
  for (let lvl = 0; lvl < 4; lvl++) {
    hs.push(2 ** -lvl);
    eV.push(4 / 3 * Math.PI * R ** 3 - signedVolume(m.positions, m.indices));
    eA.push(4 * Math.PI * R * R - surfaceArea(m.positions, m.indices));
    m = refine(m, { project: onSphere });
  }
  const oV = orderFromErrors(hs, eV), oA = orderFromErrors(hs, eA);
  close(oV.pairwise[0], 2, 0.02); close(oA.pairwise[0], 2, 0.02);
  // without projection the faceted volume does not change at all
  const flat = refine(sphere({ radius: R, subdivisions: 1 }));
  rel(signedVolume(flat.positions, flat.indices), signedVolume(sphere({ radius: R, subdivisions: 1 }).positions, sphere({ radius: R, subdivisions: 1 }).indices), 1e-6);
});

test('refineTo on a surface: every edge within the target, still watertight, same volume, groups kept', () => {
  const b = box({ size: [4, 2, 1] });
  const m = refineTo(b, 0.5);
  const q = assess(m, { panel: false });
  assert.ok(q.summary.size.max <= 0.5 + 1e-6);
  const i = inspect(m);
  assert.ok(i.watertight && i.stats.outward); assert.equal(i.stats.eulerCharacteristic, 2);
  rel(i.stats.volume, 8, 1e-6);
  assert.deepEqual(m.groups.map((g) => g.name), b.groups.map((g) => g.name));
  assert.equal(m.groups.reduce((s, g) => s + g.count, 0), m.indices.length);
  // longest-edge bisection never makes the smallest angle worse than half the original one
  assert.ok(q.summary.minAngle.min >= 0.5 * assess(b, { panel: false }).summary.minAngle.min - 1e-6);
  assert.equal(refineTo(m, 0.5).indices.length, m.indices.length, 'already fine enough');
  assert.throws(() => refineTo(b, 0), RangeError);
});

test('refineTo on volumes, frames, polylines and grids', () => {
  const v = refineTo(hexBlock(1, 2), 0.6);
  assert.equal(inspect(v).stats.cells, 64); // 2 m → 1 m → 0.5 m
  const frame = { kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 4, 3, 10]), members: [{ a: 0, b: 1, section: { type: 'tube', D: 1, t: 0.03 }, group: 'leg' }, { a: 1, b: 2, group: 'brace' }] };
  const f = refineTo(frame, 2.5);
  assert.equal(f.members.length, 4 + 2);
  assert.deepEqual(Array.from(f.memberParent), [0, 0, 0, 0, 1, 1]);
  assert.equal(f.nodes.length / 3, 3 + 3 + 1);
  assert.deepEqual(Array.from(f.nodes.subarray(0, 9)), Array.from(frame.nodes));
  assert.ok(f.members.slice(0, 4).every((m) => m.group === 'leg' && m.section.D === 1));
  const fi = inspect(f);
  assert.equal(fi.stats.components, 1); close(fi.stats.totalLength, 15, 1e-12); close(fi.stats.maxLength, 2.5, 1e-12);
  assert.equal(subdivideMembers(frame, { parts: 3 }).members.length, 6);
  assert.throws(() => subdivideMembers(frame, {}), RangeError);
  const pl = refineTo({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 10, 0, 0, 10, 10, 0]), attrs: { kp: Float64Array.from([100, 110, 120]), temp: Float64Array.from([4, 6, 10]) } }, 2);
  assert.equal(pl.points.length / 3, 11);
  close(pl.points[3 * 5], 10, 1e-12); close(pl.points[3 * 5 + 1], 0, 1e-12); // the corner falls on a sample
  close(pl.attrs.kp[10], 120, 1e-12); close(pl.attrs.temp[7], 6 + 0.4 * 4, 1e-12);
  const g = { kind: 'grid', nx: 5, ny: 3, x0: 0, y0: 0, dx: 10, dy: 10, z: Float32Array.from({ length: 15 }, (_, i) => -(i % 5) * 2 - Math.floor(i / 5)) };
  const fine = refineTo(g, 5);
  assert.equal(fine.nx, 9); assert.equal(fine.ny, 5);
  close(fine.z[1], -1, 1e-6); close(fine.z[9 + 1], -1.5, 1e-6); // bilinear on a plane is exact
  assert.equal(refineTo(g, 20), g);
  assert.throws(() => refine({ kind: 'points', positions: new Float32Array(3) }), RangeError);
});

test('grid coarsening averages the cells it replaces; polyline resampling keeps the end points', () => {
  const nx = 9, ny = 9;
  const z = new Float32Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) z[j * nx + i] = (i + j) % 2 ? 1 : -1; // checkerboard
  const g = { kind: 'grid', nx, ny, x0: 0, y0: 0, dx: 1, dy: 1, z };
  const c = resampleGrid(g, { factor: 0.5 });
  assert.equal(c.nx, 5); assert.equal(c.dx, 2);
  // interior coarse node: mean of the 3 × 3 block around it (5 × (−1) + 4 × (+1)) / 9
  close(c.z[1 * 5 + 1], -1 / 9, 1e-6);
  const cg = coarsenTo(g, 25);
  assert.equal(cg.nx * cg.ny, 25);
  const pl = resamplePolyline({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 3, 4, 0, 3, 4, 5]) }, { count: 5 });
  assert.deepEqual(Array.from(pl.points.subarray(0, 3)), [0, 0, 0]); assert.deepEqual(Array.from(pl.points.subarray(12)), [3, 4, 5]);
  assert.deepEqual(Array.from(pl.attrs.kp), [0, 2.5, 5, 7.5, 10]);
  assert.equal(coarsenTo({ kind: 'polyline', points: new Float64Array(300).map((_, i) => i) }, 10).points.length, 30);
  const pts = coarsenTo({ kind: 'points', positions: new Float32Array(3000).map((_, i) => i), attrs: { intensity: new Uint8Array(1000).map((_, i) => i % 256) } }, 100);
  assert.equal(pts.positions.length, 300); assert.equal(pts.attrs.intensity.length, 100); assert.equal(pts.positions[3], 30);
  assert.throws(() => coarsenTo(hexBlock(2), 4), RangeError);
});

test('quadric decimation keeps a sphere closed and its volume within 1 % at a fifth of the faces', () => {
  const s = sphere({ radius: 2, subdivisions: 4 });
  const v0 = signedVolume(s.positions, s.indices);
  const d = coarsenTo(s, 1000);
  assert.equal(d.indices.length / 3, 1000);
  assert.deepEqual(d.decimation.from, 5120);
  const i = inspect(d);
  assert.ok(i.watertight && i.stats.outward, (i.reasons || []).join(';'));
  assert.equal(i.stats.eulerCharacteristic, 2);
  rel(i.stats.volume, v0, 0.01);
  rel(i.stats.area, surfaceArea(s.positions, s.indices), 0.01);
  assert.ok(assess(d, { panel: false }).summary.minAngle.min > 15);
  // harder: a twentieth
  const d2 = decimate(s, 256);
  rel(inspect(d2).stats.volume, v0, 0.04);
  assert.ok(inspect(d2).watertight);
  assert.equal(decimate(s, 1e6), s);
  assert.equal(s.indices.length / 3, 5120, 'input untouched');
});

test('decimation of flat-faced and grouped bodies: exact volume, sharp edges and groups preserved', () => {
  const b = refineTo(box({ size: [4, 2, 1] }), 0.25);
  const d = decimate(b, 60);
  const i = inspect(d);
  assert.ok(i.watertight);
  rel(i.stats.volume, 8, 1e-5); rel(i.stats.area, 28, 1e-5);
  assert.ok(d.indices.length / 3 <= 60);
  assert.deepEqual(d.groups.map((g) => g.name).sort(), b.groups.map((g) => g.name).sort());
  // every face still lies in the plane of its group
  const top = d.groups.find((g) => g.name === 'z+');
  for (let k = top.start; k < top.start + top.count; k++) close(d.positions[3 * d.indices[k] + 2], 0.5, 1e-5);
  // a cylinder keeps its rim: cap area stays that of the polygon
  const c = cylinder({ radius: 1, height: 3, segments: 48, heightDivisions: 12, capRings: 4 });
  const dc = decimate(c, 400);
  assert.ok(inspect(dc).watertight);
  rel(signedVolume(dc.positions, dc.indices), signedVolume(c.positions, c.indices), 0.01);
  // an open surface keeps its boundary outline
  const hull = wigleyHull({ length: 100, beam: 10, draft: 6, nx: 40, nz: 10, closed: false });
  const dh = decimate(hull, 400);
  const ih = inspect(dh);
  assert.equal(ih.stats.boundaryLoops, 1);
  let zmax = -Infinity;
  for (let k = 2; k < dh.positions.length; k += 3) zmax = Math.max(zmax, dh.positions[k]);
  close(zmax, 0, 1e-3); // held by a penalty, not a constraint: sub-millimetre drift on a 100 m hull
  assert.throws(() => decimate({ kind: 'frame' }, 10), RangeError);
});

test('Taubin smoothing removes noise from a sphere without shrinking it; plain Laplacian shrinks', () => {
  const R = 1;
  const s = sphere({ radius: R, subdivisions: 3 });
  const rng = makeRng(11);
  const noisy = { ...s, positions: Float64Array.from(s.positions, (v) => v + 0.02 * (rng.next() - 0.5)) };
  const rough = (m) => { let e = 0; for (let i = 0; i < m.positions.length; i += 3) e += (Math.hypot(m.positions[i], m.positions[i + 1], m.positions[i + 2]) - R) ** 2; return Math.sqrt(e / (m.positions.length / 3)); };
  const t = smooth(noisy, { iterations: 20 });
  assert.ok(rough(t) < 0.5 * rough(noisy), `${rough(t)} vs ${rough(noisy)}`);
  rel(signedVolume(t.positions, t.indices), signedVolume(s.positions, s.indices), 0.01);
  const lap = smooth(noisy, { iterations: 10, taubin: false });
  assert.ok(signedVolume(lap.positions, lap.indices) < 0.97 * signedVolume(s.positions, s.indices));
  assert.equal(t.indices, noisy.indices);
  // boundaries and group borders are held: a box keeps its corners and edges exactly
  const b = refineTo(box({ size: [2, 2, 2] }), 0.5);
  const sb = smooth(b, { iterations: 5 });
  rel(signedVolume(sb.positions, sb.indices), 8, 1e-6);
  const free = smooth(b, { iterations: 5, fixGroupBorders: false, taubin: false });
  assert.ok(signedVolume(free.positions, free.indices) < 7.9);
  const fixed = smooth(noisy, { iterations: 3, fixed: [0, 1, 2] });
  assert.deepEqual(Array.from(fixed.positions.subarray(0, 9)), Array.from(noisy.positions.subarray(0, 9)));
  assert.throws(() => smooth(s, { lambda: 2 }), RangeError);
  assert.throws(() => smooth(s, { mu: -0.1 }), RangeError);
});
