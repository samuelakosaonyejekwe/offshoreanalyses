import test from 'node:test';
import assert from 'node:assert/strict';
import {
  planeBasis, surfaceArea, signedVolume, inertiaSolid, inertiaShell, principalAxes, aabb, obb, projectedArea,
  sectionProperties, planeCut, volumeBelowPlane, wallThickness, rayTriangle, polylineMeasure, gridSample, gridStats,
  gridProfile, measure,
} from '../../src/mesh/measure.js';
import { box, sphere, cylinder, monopile, tube } from '../../src/mesh/generate.js';
import { close, rel, closeV, circle, hexBlock } from './helpers.js';

test('area and volume of a box far from the origin (cancellation guarded)', () => {
  const m = box({ size: [3, 4, 5], center: [1000, -2000, 500] });
  rel(surfaceArea(m.positions, m.indices), 94, 1e-6);
  rel(signedVolume(m.positions, m.indices), 60, 1e-5);
  const flipped = Uint32Array.from(m.indices);
  for (let f = 0; f < flipped.length; f += 3) { const t = flipped[f + 1]; flipped[f + 1] = flipped[f + 2]; flipped[f + 2] = t; }
  rel(signedVolume(m.positions, flipped), -60, 1e-5);
});

test('solid inertia: parallel-axis theorem, density scaling and tensor sign convention', () => {
  const m = box({ size: [2, 4, 6], center: [3, 0, 0] });
  const ip = inertiaSolid(m.positions, m.indices, { density: 1025 });
  const M = 1025 * 48;
  rel(ip.mass, M, 1e-6);
  rel(ip.inertia[0][0], M * (16 + 36) / 12, 1e-5);
  rel(ip.inertiaOrigin[1][1], M * (4 + 36) / 12 + M * 9, 1e-5);
  rel(ip.inertiaOrigin[2][2], M * (4 + 16) / 12 + M * 9, 1e-5);
  // off-centre box rotated products: Ixy = −∫xy dm
  const o = box({ size: [2, 2, 2], center: [1, 1, 0] });
  rel(inertiaSolid(o.positions, o.indices).inertiaOrigin[0][1], -8 * 1 * 1, 1e-5);
  assert.throws(() => inertiaSolid(new Float32Array(9), Uint32Array.from([0, 1, 2])), RangeError);
});

test('thin spherical shell: area 4πR², I = 2/3 m R²', () => {
  const R = 3;
  const m = sphere({ radius: R, subdivisions: 5 });
  const sh = inertiaShell(m.positions, m.indices, { thickness: 0.02, density: 7850 });
  rel(sh.area, 4 * Math.PI * R * R, 1e-3);
  rel(sh.mass, 0.02 * 7850 * 4 * Math.PI * R * R, 1e-3);
  rel(sh.inertia[2][2], 2 / 3 * sh.mass * R * R, 2e-3);
  closeV(sh.centroid, [0, 0, 0], 1e-4);
});

test('principal axes of a tilted member recover its axis and radii of gyration', () => {
  const a = [0, 0, 0], b = [6, 8, 0];
  const m = tube({ a, b, radius: 0.5, segments: 64 });
  const ip = inertiaSolid(m.positions, m.indices);
  const pa = principalAxes(ip.inertia, ip.mass);
  assert.ok(pa.moments[0] < pa.moments[1]);
  // smallest moment about the member axis (0.6, 0.8, 0)
  close(Math.abs(pa.axes[0][0] * 0.6 + pa.axes[0][1] * 0.8), 1, 1e-3);
  rel(pa.radiiOfGyration[0], 0.5 / Math.SQRT2, 2e-3);
  rel(pa.radiiOfGyration[2], Math.sqrt(0.25 / 4 + 100 / 12), 2e-3);
});

test('bounding boxes: axis-aligned and principal-axis oriented', () => {
  const m = tube({ a: [0, 0, 0], b: [6, 8, 0], radius: 0.5, segments: 64 });
  const bb = aabb(m.positions);
  close(bb.size[2], 1, 1e-5);
  closeV(bb.center, [3, 4, 0], 0.5);
  const ob = obb(m.positions);
  close(ob.size[0], 10, 1e-3);
  close(ob.size[1], 1, 0.02); close(ob.size[2], 1, 0.02);
  assert.ok(ob.volume < (bb.size[0] * bb.size[1] * bb.size[2]) / 3);
  assert.throws(() => obb(new Float32Array(0)), RangeError);
});

test('projected (silhouette) area of a sphere and of a box seen along a diagonal', () => {
  const s = sphere({ radius: 2, subdivisions: 4 });
  const pa = projectedArea(s.positions, s.indices, [1, 2, 3], { resolution: 512 });
  rel(pa.area, Math.PI * 4, 0.01);
  assert.ok(pa.errorBound > Math.abs(pa.area - Math.PI * 4) / 4);
  const b = box({ size: [1, 1, 1] });
  // unit cube along (1,1,1): regular hexagon of area √3
  rel(projectedArea(b.positions, b.indices, [1, 1, 1]).area, Math.sqrt(3), 0.01);
  rel(projectedArea(b.positions, b.indices, [0, 0, 5]).area, 1, 0.005);
});

test('section properties: rectangle, hollow box by nesting, rotated rectangle principal axes', () => {
  const r = sectionProperties([Float64Array.from([0, 0, 4, 0, 4, 2, 0, 2])]);
  close(r.area, 8, 1e-12); closeV(r.centroid, [2, 1], 1e-12);
  close(r.Ixx, 4 * 8 / 12, 1e-12); close(r.Iyy, 2 * 64 / 12, 1e-12); close(r.Ixy, 0, 1e-12);
  close(r.perimeter, 12, 1e-12); close(r.rx, Math.sqrt(1 / 3), 1e-12);
  close(r.I1, 128 / 12, 1e-12); close(Math.abs(r.thetaDeg) % 180, 90, 1e-9);
  // hollow: outer 4 × 2, hole 2 × 1 (given clockwise or counter-clockwise alike), closing point repeated
  const hollow = sectionProperties([Float64Array.from([0, 0, 4, 0, 4, 2, 0, 2, 0, 0]), Float64Array.from([1, 0.5, 1, 1.5, 3, 1.5, 3, 0.5])]);
  close(hollow.area, 6, 1e-12);
  close(hollow.Ixx, (4 * 8 - 2 * 1) / 12, 1e-12);
  close(hollow.Iyy, (2 * 64 - 1 * 8) / 12, 1e-12);
  // circular tube against π (R⁴ − r⁴) / 4 with polygon correction
  const n = 720;
  const t = sectionProperties([circle(1, n), circle(0.8, n)]);
  rel(t.Ixx, Math.PI * (1 - 0.8 ** 4) / 4, 1e-4);
  rel(t.area, Math.PI * (1 - 0.64), 1e-4);
  assert.throws(() => sectionProperties([]), RangeError);
});

test('plane cut of a cylinder: circular section properties and loop orientation', () => {
  const R = 1.5;
  const m = cylinder({ radius: R, height: 10, segments: 360, base: [2, 3, -5] });
  const cut = planeCut(m.positions, m.indices, { point: [0, 0, 0.123], normal: [0, 0, 2] });
  assert.equal(cut.loops.length, 1); assert.deepEqual(cut.closed, [true]);
  rel(cut.area, Math.PI * R * R, 1e-4);
  closeV(cut.centroid, [2, 3, 0.123], 1e-4);
  rel(cut.Iuu, Math.PI * R ** 4 / 4, 3e-4); rel(cut.Ivv, Math.PI * R ** 4 / 4, 3e-4);
  rel(cut.perimeter, 2 * Math.PI * R, 1e-4);
  // inclined cut through the axis: ellipse of area π R² / cos θ
  const inc = planeCut(m.positions, m.indices, { point: [2, 3, 0], normal: [Math.sin(0.5), 0, Math.cos(0.5)] });
  rel(inc.area, Math.PI * R * R / Math.cos(0.5), 1e-3);
  // a plane that misses gives no loops
  assert.equal(planeCut(m.positions, m.indices, { point: [0, 0, 50], normal: [0, 0, 1] }).loops.length, 0);
  assert.throws(() => planeCut(m.positions, m.indices, { point: [0, 0, 0] }), RangeError);
});

test('volume below a plane: hemisphere 2πR³/3 with centroid 3R/8 below the cut; heeled box', () => {
  const R = 2;
  const s = sphere({ radius: R, subdivisions: 5 });
  const h = volumeBelowPlane(s.positions, s.indices, { point: [0, 0, 0], normal: [0, 0, 1] });
  rel(h.volume, 2 / 3 * Math.PI * R ** 3, 2e-3);
  close(h.centroid[2], -3 * R / 8, 2e-3);
  rel(h.wettedArea, 2 * Math.PI * R * R, 2e-3);
  // spherical cap of height hc: V = π hc² (3R − hc) / 3
  const cap = volumeBelowPlane(s.positions, s.indices, { point: [0, 0, -1.2], normal: [0, 0, 1] });
  rel(cap.volume, Math.PI * 0.8 ** 2 * (3 * R - 0.8) / 3, 5e-3);
  // unit-width box heeled so that the waterplane passes through its centre: half the volume whatever the angle
  const b = box({ size: [4, 2, 2] });
  for (const ang of [0, 0.2, 0.7]) rel(volumeBelowPlane(b.positions, b.indices, { point: [0, 0, 0], normal: [0, Math.sin(ang), Math.cos(ang)] }).volume, 8, 1e-6);
  const heel = volumeBelowPlane(b.positions, b.indices, { point: [0, 0, 0], normal: [0, Math.sin(0.2), Math.cos(0.2)] });
  // wedge transfer moves the centre of buoyancy towards the immersed side: y_B = −(b²/12T)·tan θ about the box
  close(heel.centroid[1], -(4 / 12) * Math.tan(0.2), 1e-5);
  assert.equal(volumeBelowPlane(b.positions, b.indices, { point: [0, 0, -5], normal: [0, 0, 1] }).volume, 0);
});

test('wall thickness of a hollow pile by inward rays', () => {
  const m = monopile({ diameter: 6, zTop: 10, zToe: 0, wallThickness: 0.25, segments: 96, elementSize: 1 });
  const w = wallThickness(m.positions, m.indices, { samples: 200, seed: 4 });
  assert.equal(w.misses, 0);
  close(w.median, 0.25, 0.01);
  assert.ok(w.min > 0.2);
  close(rayTriangle([0.2, 0.2, 5], [0, 0, -1], Float64Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0]), 0, 1, 2), 5, 1e-12);
  assert.equal(rayTriangle([2, 2, 5], [0, 0, -1], Float64Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0]), 0, 1, 2), Infinity);
});

test('polyline: length, inclination and heading against chainage, radius of a circular arc', () => {
  const p = polylineMeasure(Float64Array.from([0, 0, 0, 30, 40, 0, 30, 40, -50, 30, 140, -50]));
  close(p.length, 200, 1e-12); close(p.horizontalLength, 150, 1e-12);
  closeV(Array.from(p.chainage), [0, 50, 100, 200], 1e-12);
  closeV(Array.from(p.inclinationDeg), [0, -90, 0, 0], 1e-12);
  close(p.headingDeg[0], Math.atan2(30, 40) * 180 / Math.PI, 1e-12); // from north, clockwise
  close(p.headingDeg[2], 0, 1e-12);
  close(p.maxInclinationDeg, 90, 1e-12); close(p.zMin, -50, 0);
  const R = 500, n = 90;
  const arc = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) arc.set([R * Math.cos(i * 0.01), R * Math.sin(i * 0.01), -100], 3 * i);
  const a = polylineMeasure(arc);
  rel(a.minRadius, R, 1e-6);
  rel(a.length, R * 0.89, 1e-4);
  rel(a.curvature[40], 1 / R, 1e-6);
  assert.equal(polylineMeasure(Float64Array.from([0, 0, 0, 1, 0, 0, 2, 0, 0])).minRadius, Infinity);
  // constant 5° slope
  const s = polylineMeasure(Float64Array.from([0, 0, 0, 100, 0, -100 * Math.tan(5 * Math.PI / 180)]));
  close(s.inclinationDeg[0], -5, 1e-9);
  assert.throws(() => polylineMeasure(Float64Array.from([0, 0, 0])), RangeError);
});

const slope = (nx = 21, ny = 11) => {
  const z = new Float32Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) z[j * nx + i] = -30 - 0.1 * i * 5;
  return { kind: 'grid', name: 'slope', nx, ny, x0: 1000, y0: 2000, dx: 5, dy: 5, z };
};

test('grid statistics of a plane sloping 1 : 10', () => {
  const g = slope();
  const s = gridStats(g);
  assert.equal(s.nodes, 231); assert.equal(s.valid, 231); assert.equal(s.noData, 0);
  close(s.zMax, -30, 1e-5); close(s.zMin, -40, 1e-5); close(s.depthMax, 40, 1e-5); close(s.zMean, -35, 1e-4);
  close(s.slopeMeanDeg, Math.atan(0.1) * 180 / Math.PI, 1e-4); close(s.slopeMaxDeg, s.slopeMeanDeg, 1e-4);
  close(s.roughness, 0, 1e-5);
  close(s.extentX, 100, 0); close(s.extentY, 50, 0); close(s.area, 5000, 0);
  g.z[5] = NaN;
  assert.equal(gridStats(g).noData, 1);
  close(gridSample(slope(), 1012.5, 2007), -31.25, 1e-5);
  assert.ok(Number.isNaN(gridSample(slope(), 999, 2007)));
  assert.throws(() => gridStats({ nx: 2, ny: 2, z: new Float32Array(3) }), RangeError);
});

test('grid profile under a route: chainage, seabed level, gap and slope along the line', () => {
  const g = slope();
  const route = Float64Array.from([1000, 2025, -29, 1060, 2025, -36, 1060, 2045, -36, 1200, 2045, -36]);
  const p = gridProfile(g, route);
  close(p.chainage[p.chainage.length - 1], 60 + 20 + 140, 1e-9);
  assert.equal(p.chainage.length, 12 + 4 + 28 + 1);
  close(p.z[0], -30, 1e-5); close(p.gap[0], 1, 1e-5);
  // at x = 1060 the seabed is at −36: the route touches down
  close(p.z[12], -36, 1e-5); close(p.gap[12], 0, 1e-5);
  close(p.slopeDeg[5], -Math.atan(0.1) * 180 / Math.PI, 1e-4);
  close(p.slopeDeg[14], 0, 1e-4); // running along the contour
  assert.equal(p.outside, 20); // beyond x = 1100
  assert.ok(Number.isNaN(p.z[p.z.length - 1]));
  const coarse = gridProfile(g, Float64Array.from([1000, 2000, 1100, 2050]), { step: 50, stride: 2 });
  assert.equal(coarse.chainage.length, 4);
  assert.ok(Number.isNaN(coarse.gap[0]));
  assert.throws(() => gridProfile(g, Float64Array.from([0, 0, 0])), RangeError);
});

test('measure(): headline numbers for every model kind', () => {
  const b = measure(box({ size: [1, 2, 3] }));
  rel(b.area, 22, 1e-6); rel(b.volume, 6, 1e-6); closeV(b.centroid, [0, 0, 0], 1e-6); assert.equal(b.faces, 12);
  const v = measure(hexBlock(2, 3));
  close(v.volume, 27, 1e-12); close(v.boundaryArea, 54, 1e-12); assert.equal(v.cells, 8);
  const f = measure({ kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 3, 4, 10]), members: [{ a: 0, b: 1, section: { type: 'tube', D: 1, t: 0.05 } }, { a: 1, b: 2, section: { A: 0.2 } }] });
  close(f.totalLength, 15, 1e-12); close(f.steelVolume, 10 * Math.PI * 0.95 * 0.05 + 5 * 0.2, 1e-12);
  close(measure({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 3, 4, 0]) }).length, 5, 1e-12);
  assert.equal(measure(slope()).kind, 'grid');
  assert.equal(measure({ kind: 'points', positions: new Float32Array(30) }).count, 10);
  close(measure({ kind: 'section', loops: [circle(1, 360)] }).area, Math.PI, 1e-3);
  assert.throws(() => measure({ kind: 'hologram' }), RangeError);
  closeV(planeBasis([0, 0, 3]).u, [1, 0, 0], 1e-12); closeV(planeBasis([0, 0, 3]).v, [0, 1, 0], 1e-12);
});
