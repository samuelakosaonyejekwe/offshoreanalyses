import test from 'node:test';
import assert from 'node:assert/strict';
import { box, cylinder, cone, sphere, tube, barge, wigleyHull, semiSubmersible, spar, monopile, triangulateSection, delaunay, pointsToGrid, gridToSurface } from '../../src/mesh/generate.js';
import { inspect, faceGroupIds } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume, inertiaSolid, volumeBelowPlane, planeCut, gridSample } from '../../src/mesh/measure.js';
import { orderFromErrors } from '../../src/mesh/convergence.js';
import { makeRng } from '../../src/engine/math/rng.js';
import { close, rel, closeV, circle } from './helpers.js';

const checkClosed = (m, chi = 2) => {
  const r = inspect(m);
  assert.equal(r.stats.eulerCharacteristic, chi, `${m.name}: V − E + F`);
  assert.ok(r.watertight, `${m.name}: watertight (${(r.reasons || []).join('; ')})`);
  assert.ok(r.stats.outward, `${m.name}: outward normals`);
  assert.ok(r.ok, `${m.name}: no errors`);
  const g = faceGroupIds(m);
  assert.ok(g.ids.every((i) => i >= 0), `${m.name}: every face belongs to a named group`);
  assert.equal(m.groups.reduce((s, x) => s + x.count, 0), m.indices.length, `${m.name}: groups cover the index array`);
  assert.ok(m.positions instanceof Float32Array && m.indices instanceof Uint32Array);
  return r;
};

test('box: V − E + F = 2, exact area, volume and inertia (I = m (b² + c²) / 12)', () => {
  const m = box({ size: [2, 3, 4], center: [1, -2, 0.5], elementSize: 0.5 });
  const r = checkClosed(m);
  assert.deepEqual(m.groups.map((g) => g.name).sort(), ['x+', 'x-', 'y+', 'y-', 'z+', 'z-']);
  rel(r.stats.area, 2 * (6 + 8 + 12), 1e-6);
  const ip = inertiaSolid(m.positions, m.indices);
  rel(ip.volume, 24, 1e-6);
  closeV(ip.centroid, [1, -2, 0.5], 1e-5);
  rel(ip.inertia[0][0], 24 * (9 + 16) / 12, 1e-5);
  rel(ip.inertia[1][1], 24 * (4 + 16) / 12, 1e-5);
  rel(ip.inertia[2][2], 24 * (4 + 9) / 12, 1e-5);
  close(ip.inertia[0][1], 0, 1e-4); close(ip.inertia[1][2], 0, 1e-4);
  assert.equal(box({ size: [1, 1, 1] }).indices.length / 3, 12);
  assert.throws(() => box({ size: [1, 0, 1] }), RangeError);
});

test('cylinder: inscribed prism values exactly, analytic values at second order', () => {
  const R = 1.5, H = 4;
  const hs = [], eV = [], eA = [];
  for (const n of [16, 32, 64]) {
    const m = cylinder({ radius: R, height: H, segments: n, heightDivisions: n / 8, capRings: n / 8 });
    const r = checkClosed(m);
    assert.deepEqual(m.groups.map((g) => g.name), ['bottom', 'side', 'top']);
    // exact for the n-gon prism
    rel(r.stats.volume, n / 2 * R * R * Math.sin(2 * Math.PI / n) * H, 1e-5);
    rel(r.stats.area, n * R * R * Math.sin(2 * Math.PI / n) + 2 * n * R * Math.sin(Math.PI / n) * H, 1e-5);
    hs.push(1 / n); eV.push(Math.PI * R * R * H - r.stats.volume); eA.push(2 * Math.PI * R * (R + H) - r.stats.area);
  }
  close(orderFromErrors(hs, eV).p, 2, 0.05, 'volume order');
  close(orderFromErrors(hs, eA).p, 2, 0.05, 'area order');
  const ip = inertiaSolid(cylinder({ radius: R, height: H, segments: 256 }).positions, cylinder({ radius: R, height: H, segments: 256 }).indices);
  const M = Math.PI * R * R * H;
  rel(ip.inertia[2][2], M * R * R / 2, 2e-3, 'Izz = m R² / 2');
  rel(ip.inertia[0][0], M * (3 * R * R + H * H) / 12, 2e-3, 'Ixx = m (3R² + h²) / 12');
  closeV(ip.centroid, [0, 0, H / 2], 1e-5);
});

test('cone and frustum: volume π h (R² + R r + r²) / 3, centroid at h / 4 for a full cone', () => {
  const c = cone({ radiusBottom: 2, radiusTop: 0, height: 6, segments: 128, heightDivisions: 8 });
  checkClosed(c);
  const ip = inertiaSolid(c.positions, c.indices);
  rel(ip.volume, Math.PI * 4 * 6 / 3, 1e-3);
  close(ip.centroid[2], 6 / 4, 1e-3);
  const f = cone({ radiusBottom: 2, radiusTop: 1, height: 3, segments: 128 });
  checkClosed(f);
  rel(signedVolume(f.positions, f.indices), Math.PI * 3 * (4 + 2 + 1) / 3, 1e-3);
  checkClosed(cone({ radiusBottom: 0, radiusTop: 1, height: 1, segments: 12 }));
  assert.throws(() => cone({ radiusBottom: 0, radiusTop: 0, height: 1 }), RangeError);
});

test('sphere: 20·4ˢ triangles, area 4πR², volume 4πR³/3 and inertia 2mR²/5 converge at second order', () => {
  const R = 2;
  const hs = [], eV = [], eA = [], eI = [];
  for (const s of [1, 2, 3, 4]) {
    const m = sphere({ radius: R, subdivisions: s });
    assert.equal(m.indices.length / 3, 20 * 4 ** s);
    const r = checkClosed(m);
    const ip = inertiaSolid(m.positions, m.indices);
    hs.push(2 ** -s); eV.push(4 / 3 * Math.PI * R ** 3 - r.stats.volume); eA.push(4 * Math.PI * R * R - r.stats.area);
    eI.push(0.4 * (4 / 3 * Math.PI * R ** 3) * R * R - ip.inertia[0][0]);
    closeV(ip.centroid, [0, 0, 0], 1e-5);
  }
  assert.ok(eV.every((e) => e > 0) && eA.every((e) => e > 0), 'inscribed polyhedron is smaller than the sphere');
  close(orderFromErrors(hs, eV).pairwise[0], 2, 0.05, 'volume order');
  close(orderFromErrors(hs, eA).pairwise[0], 2, 0.05, 'area order');
  close(orderFromErrors(hs, eI).pairwise[0], 2, 0.06, 'inertia order');
  assert.ok(eV[3] / (4 / 3 * Math.PI * R ** 3) < 3e-3);
  assert.equal(sphere({ radius: 1, elementSize: 0.13 }).indices.length / 3, 20 * 4 ** 3);
});

test('tube between two points: axis, length and volume', () => {
  const a = [1, 2, 3], b = [4, 6, 15];
  const m = tube({ a, b, radius: 0.5, segments: 96 });
  checkClosed(m);
  const ip = inertiaSolid(m.positions, m.indices);
  rel(ip.volume, Math.PI * 0.25 * 13, 1e-3);
  closeV(ip.centroid, [2.5, 4, 9], 1e-4);
  const t = tube({ a: [0, 0, 0], b: [0, 0, -10], radius: 1, radiusB: 0.5, segments: 64 });
  checkClosed(t);
  rel(signedVolume(t.positions, t.indices), Math.PI * 10 * (1 + 0.5 + 0.25) / 3, 2e-3);
  assert.throws(() => tube({ a: [0, 0, 0], b: [0, 0, 0], radius: 1 }), RangeError);
});

test('barge: box and raked volumes are exact; displacement below the waterline', () => {
  const m = barge({ length: 100, beam: 30, depth: 8, draft: 5, elementSize: 4 });
  checkClosed(m);
  assert.deepEqual(m.groups.map((g) => g.name).sort(), ['bottom', 'bow', 'deck', 'port', 'starboard', 'stern']);
  rel(signedVolume(m.positions, m.indices), 24000, 1e-6);
  const wet = volumeBelowPlane(m.positions, m.indices, { point: [0, 0, 0], normal: [0, 0, 1] });
  rel(wet.volume, 100 * 30 * 5, 1e-6);
  close(wet.centroid[2], -2.5, 1e-5);
  rel(wet.wettedArea, 100 * 30 + 2 * 5 * 130, 1e-6);
  const r = barge({ length: 100, beam: 30, depth: 8, draft: 5, rakeBow: 10, rakeStern: 6, rakeHeight: 4, elementSize: 4 });
  checkClosed(r);
  rel(signedVolume(r.positions, r.indices), 24000 - 30 * 16 * 4 / 2, 1e-6);
  assert.throws(() => barge({ length: 10, beam: 5, depth: 2, rakeBow: 6, rakeStern: 6 }), RangeError);
});

test('Wigley hull: volume 4/9 LBT, waterplane 2/3 LB, KB at −3T/8, second-order convergence', () => {
  const L = 100, B = 10, T = 6.25;
  const hs = [], eV = [];
  for (const n of [20, 40, 80]) {
    const m = wigleyHull({ length: L, beam: B, draft: T, nx: n, nz: n / 4 });
    const r = checkClosed(m);
    hs.push(1 / n); eV.push(4 / 9 * L * B * T - r.stats.volume);
  }
  close(orderFromErrors(hs, eV).p, 2, 0.1);
  const m = wigleyHull({ length: L, beam: B, draft: T, nx: 160, nz: 40 });
  const ip = inertiaSolid(m.positions, m.indices);
  rel(ip.volume, 4 / 9 * L * B * T, 1e-3);
  close(ip.centroid[2], -3 * T / 8, 5e-3);
  close(ip.centroid[0], 0, 1e-4); close(ip.centroid[1], 0, 1e-4);
  const wp = planeCut(m.positions, m.indices, { point: [0, 0, -1e-4], normal: [0, 0, 1] });
  rel(wp.area, 2 / 3 * L * B, 1e-3);
  const open = inspect(wigleyHull({ length: L, beam: B, draft: T, nx: 20, nz: 6, closed: false }));
  assert.equal(open.watertight, false);
  assert.equal(open.stats.boundaryLoops, 1);
});

test('spar and monopile: closed bodies of revolution; hollow pile has genus 1', () => {
  // an n-sided body of revolution has n·sin(2π/n)/(2π) of the volume of the round one, section by section
  const n = 48, k = n * Math.sin(2 * Math.PI / n) / (2 * Math.PI);
  const s = spar({ diameter: 9.4, draft: 120, freeboard: 10, upperDiameter: 6.5, taperTop: -4, taperBottom: -12, segments: n, elementSize: 4 });
  checkClosed(s);
  assert.deepEqual(s.groups.map((g) => g.name), ['keel', 'lower', 'taper', 'upper', 'top']);
  const exact = Math.PI * (4.7 ** 2 * 108 + 8 * (4.7 ** 2 + 4.7 * 3.25 + 3.25 ** 2) / 3 + 3.25 ** 2 * 14);
  rel(signedVolume(s.positions, s.indices), k * exact, 1e-5);
  const wet = volumeBelowPlane(s.positions, s.indices, { point: [0, 0, 0], normal: [0, 0, 1] });
  rel(wet.volume, k * (exact - Math.PI * 3.25 ** 2 * 10), 1e-5);
  const p = monopile({ diameter: 8, zTop: 20, zToe: -60, segments: n, elementSize: 2 });
  checkClosed(p);
  rel(signedVolume(p.positions, p.indices), k * Math.PI * 16 * 80, 1e-5);
  const h = monopile({ diameter: 8, zTop: 20, zToe: -60, wallThickness: 0.08, segments: n, elementSize: 2 });
  const r = checkClosed(h, 0);
  assert.equal(r.stats.genus, 1);
  rel(r.stats.volume, k * Math.PI * (16 - 3.92 ** 2) * 80, 1e-4, 'steel volume');
  assert.ok(h.groups.some((g) => g.name === 'inner'));
  const t = monopile({ diameter: 8, topDiameter: 6, taperBottom: -5, taperTop: 5, zTop: 20, zToe: -60, segments: n, elementSize: 2 });
  checkClosed(t);
  rel(signedVolume(t.positions, t.indices), k * Math.PI * (16 * 55 + 10 * (16 + 12 + 9) / 3 + 9 * 15), 1e-5);
  assert.throws(() => monopile({ diameter: 8, zTop: -60, zToe: 20 }), RangeError);
});

test('semi-submersible: two watertight hulls without internal faces, volume of pontoons plus columns', () => {
  const m = semiSubmersible({ pontoonLength: 100, pontoonWidth: 16, pontoonHeight: 8, pontoonSpacing: 60, columnDiameter: 12, columnHeight: 30, draft: 20, segments: 32 });
  const r = checkClosed(m, 4);
  assert.equal(r.stats.components, 2);
  assert.deepEqual(m.groups.map((g) => g.name).sort(), ['column-tops', 'columns', 'pontoons']);
  const colArea = 32 / 2 * 36 * Math.sin(2 * Math.PI / 32);
  rel(r.stats.volume, 2 * 100 * 16 * 8 + 4 * 30 * colArea, 1e-5);
  const wet = volumeBelowPlane(m.positions, m.indices, { point: [0, 0, 0], normal: [0, 0, 1] });
  rel(wet.volume, 2 * 100 * 16 * 8 + 4 * 12 * colArea, 1e-5);
  const wp = planeCut(m.positions, m.indices, { point: [0, 0, 0], normal: [0, 0, 1] });
  assert.equal(wp.loops.length, 4);
  rel(wp.area, 4 * colArea, 1e-5);
  assert.throws(() => semiSubmersible({ pontoonLength: 100, pontoonWidth: 10, pontoonHeight: 8, pontoonSpacing: 60, columnDiameter: 12, columnHeight: 30 }), RangeError);
});

const emptyCircles = (xy, tris, tol = 1e-9) => {
  const n = xy.length / 2;
  let worst = 0;
  for (let t = 0; t < tris.length; t += 3) {
    const [a, b, c] = [tris[t], tris[t + 1], tris[t + 2]];
    const ax = xy[2 * a], ay = xy[2 * a + 1], bx = xy[2 * b] - ax, by = xy[2 * b + 1] - ay, cx = xy[2 * c] - ax, cy = xy[2 * c + 1] - ay;
    const d = 2 * (bx * cy - by * cx);
    assert.ok(d > 0, 'counter-clockwise triangle');
    const ux = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d, uy = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d;
    const R = Math.hypot(ux, uy);
    for (let v = 0; v < n; v++) {
      if (v === a || v === b || v === c) continue;
      const dist = Math.hypot(xy[2 * v] - ax - ux, xy[2 * v + 1] - ay - uy);
      if (R - dist > worst * R) worst = (R - dist) / R;
    }
  }
  assert.ok(worst <= tol, `a vertex lies ${worst} radii inside a circumcircle`);
};

test('Delaunay triangulation of random points: empty circumcircles and 2n − 2 − h triangles', () => {
  const rng = makeRng(7);
  const n = 400;
  const xy = new Float64Array(2 * n);
  for (let i = 0; i < 2 * n; i++) xy[i] = rng.next() * 10;
  // corners of the square make the convex hull known: h = 4
  xy.set([-1, -1, 11, -1, 11, 11, -1, 11]);
  const tris = delaunay(xy);
  assert.equal(tris.length / 3, 2 * n - 2 - 4);
  emptyCircles(xy, tris);
  assert.equal(delaunay([0, 0, 1, 0]).length, 0);
  // four cocircular points still give two triangles; a repeated point is used once
  assert.equal(delaunay([0, 0, 1, 0, 1, 1, 0, 1, 1, 1]).length / 3, 2);
});

test('section meshing: Delaunay empty-circle property, minimum angle, size, area and boundary conformity', () => {
  const sq = triangulateSection({ kind: 'section', name: 'plate', loops: [Float64Array.from([0, 0, 2, 0, 2, 1, 0, 1])] }, { size: 0.1, minAngle: 28 });
  emptyCircles(sq.xy, sq.indices);
  assert.ok(sq.stats.minAngleDeg >= 28, `min angle ${sq.stats.minAngleDeg}`);
  assert.ok(sq.stats.maxEdge <= 0.15);
  close(sq.stats.area, 2, 1e-12);
  assert.equal(sq.stats.truncated, false);
  const r = inspect(sq, { selfIntersections: false });
  assert.equal(r.stats.eulerCharacteristic, 1); // a disc
  assert.equal(r.stats.boundaryLoops, 1);
  assert.equal(r.stats.boundaryEdges, sq.boundary.length / 2);
  // every boundary vertex lies on the outline
  for (const v of sq.boundary) { const x = sq.xy[2 * v], y = sq.xy[2 * v + 1]; assert.ok(Math.min(Math.abs(x), Math.abs(x - 2), Math.abs(y), Math.abs(y - 1)) < 1e-12); }
  assert.ok(sq.positions.every((v, i) => i % 3 !== 2 || v === 0));
});

test('section meshing: tube with a hole (genus via Euler characteristic 0) and area of the annulus polygon', () => {
  const n = 48;
  const m = triangulateSection({ kind: 'section', loops: [circle(1, n), circle(0.9, n)] }, { size: 0.05 });
  emptyCircles(m.xy, m.indices);
  const r = inspect(m, { selfIntersections: false });
  assert.equal(r.stats.eulerCharacteristic, 0);
  assert.equal(r.stats.boundaryLoops, 2);
  close(m.stats.area, n / 2 * Math.sin(2 * Math.PI / n) * (1 - 0.81), 1e-12);
  assert.ok(m.stats.minAngleDeg >= 25);
  // hole winding does not matter
  const rev = circle(0.9, n).reverse();
  const hole = Float64Array.from({ length: 2 * n }, (_, i) => rev[i % 2 ? i - 1 : i + 1]);
  close(triangulateSection({ loops: [circle(1, n), hole] }, { size: 0.05 }).stats.area, m.stats.area, 1e-12);
});

test('section meshing: without a size only quality is enforced; small input angles end the refinement cleanly', () => {
  const coarse = triangulateSection({ loops: [Float64Array.from([0, 0, 1, 0, 1, 1, 0, 1])] });
  assert.equal(coarse.stats.triangles, 2);
  const L = triangulateSection({ loops: [Float64Array.from([0, 0, 2, 0, 2, 0.2, 0.2, 0.2, 0.2, 3, 0, 3, 0, 0])] }, { minAngle: 30 });
  assert.ok(L.stats.minAngleDeg >= 30 && !L.stats.truncated);
  close(L.stats.area, 0.4 + 0.2 * 2.8, 1e-12);
  emptyCircles(L.xy, L.indices);
  // a 5.7° wedge cannot be improved at its tip, but the mesher must stop and respect the area
  const wedge = triangulateSection({ loops: [Float64Array.from([0, 0, 10, 0, 10, 1])] }, { size: 0.5 });
  assert.equal(wedge.stats.truncated, false);
  close(wedge.stats.area, 5, 1e-10);
  close(wedge.stats.minAngleDeg, Math.atan(0.1) * 180 / Math.PI, 0.5);
  assert.throws(() => triangulateSection({ loops: [] }), RangeError);
  assert.throws(() => triangulateSection({ loops: [Float64Array.from([0, 0, 1, 0, 1, 1])] }, { minAngle: 40 }), RangeError);
});

test('points to grid by inverse-distance weighting: exact at data points, bounded by the data, NaN beyond the radius', () => {
  const rng = makeRng(3);
  const n = 4000;
  const P = new Float32Array(3 * n);
  const plane = (x, y) => -20 - 0.02 * x + 0.01 * y;
  for (let i = 0; i < n; i++) { const x = rng.next() * 100, y = rng.next() * 50; P.set([x, y, plane(x, y)], 3 * i); }
  // pin the four corners so the grid extent is known
  P.set([0, 0, plane(0, 0), 100, 0, plane(100, 0), 100, 50, plane(100, 50), 0, 50, plane(0, 50)]);
  const g = pointsToGrid({ kind: 'points', name: 'survey', positions: P }, { dx: 5 });
  assert.equal(g.kind, 'grid'); assert.equal(g.nx, 21); assert.equal(g.ny, 11);
  close(g.z[0], plane(0, 0), 1e-5, 'node on a data point');
  let worst = 0;
  for (let j = 0; j < g.ny; j++) for (let i = 0; i < g.nx; i++) {
    const z = g.z[j * g.nx + i];
    assert.ok(z >= -22.01 && z <= -19.49, 'IDW never leaves the range of the data');
    worst = Math.max(worst, Math.abs(z - plane(i * 5, j * 5)));
  }
  // IDW is not linear-exact; with ~40 points per cell the error stays a fraction of the fall across one cell (0.1 m)
  assert.ok(worst < 0.08, `worst error on a plane ${worst}`);
  const sparse = pointsToGrid({ positions: Float32Array.from([0, 0, 1, 100, 0, 3]) }, { dx: 10, radius: 15, bounds: { x0: 0, y0: 0, x1: 100, y1: 0 } });
  assert.equal(sparse.nx, 11); assert.equal(sparse.ny, 1);
  assert.equal(sparse.z[0], 1); assert.equal(sparse.z[10], 3);
  assert.ok(Number.isNaN(sparse.z[5]));
  // two points at equal distance: the mean, whatever the power
  close(pointsToGrid({ positions: Float32Array.from([0, 0, 1, 100, 0, 3]) }, { dx: 50, radius: 60, power: 3 }).z[1], 2, 1e-6);
  assert.throws(() => pointsToGrid({ positions: P }, { dx: 0 }), RangeError);
});

test('grid to surface: two triangles per valid cell, upward normals, grid node map', () => {
  const grid = { kind: 'grid', name: 'g', nx: 4, ny: 3, x0: 10, y0: 20, dx: 2, dy: 2, z: Float32Array.from([0, 0, 0, 0, 0, 1, NaN, 0, 0, 0, 0, 0]) };
  const s = gridToSurface(grid);
  assert.equal(s.indices.length / 3, 2 * (6 - 4));
  assert.equal(s.positions.length / 3, 6); // the two cells left share an edge; unused nodes are dropped
  const [a, b, c] = s.indices;
  const u = [s.positions[3 * b] - s.positions[3 * a], s.positions[3 * b + 1] - s.positions[3 * a + 1]], v = [s.positions[3 * c] - s.positions[3 * a], s.positions[3 * c + 1] - s.positions[3 * a + 1]];
  assert.ok(u[0] * v[1] - u[1] * v[0] > 0);
  for (let i = 0; i < s.gridNode.length; i++) close(s.positions[3 * i + 2], grid.z[s.gridNode[i]], 0);
  grid.z[6] = 2;
  const full = gridToSurface(grid);
  assert.equal(full.indices.length / 3, 12);
  close(surfaceArea(gridToSurface({ ...grid, z: new Float32Array(12) }).positions, gridToSurface({ ...grid, z: new Float32Array(12) }).indices), 24, 1e-6);
  assert.equal(gridToSurface({ ...grid, nx: 5, ny: 5, z: new Float32Array(25) }, { stride: 2 }).indices.length / 3, 8);
  close(gridSample(grid, 12, 22), 1, 1e-12);
});
