import test from 'node:test';
import assert from 'node:assert/strict';
import { triQuality, quadQuality, cellQuality, panelFrequencyLimit, panelCheck, assess, QUALITY_LIMITS } from '../../src/mesh/quality.js';
import { box, sphere, barge, wigleyHull } from '../../src/mesh/generate.js';
import { close, rel, hexBlock } from './helpers.js';

const P = (pts) => Float64Array.from(pts.flat());
const S3 = Math.sqrt(3);

test('equilateral triangle is ideal on every metric', () => {
  const q = triQuality(P([[0, 0, 0], [2, 0, 0], [1, S3, 0]]), 0, 1, 2);
  close(q.aspectRatio, 1, 1e-12); close(q.minAngle, 60, 1e-9); close(q.maxAngle, 60, 1e-9);
  close(q.skewness, 0, 1e-12); close(q.scaledJacobian, 1, 1e-12); close(q.size, 2, 1e-12); close(q.area, S3, 1e-12);
});

test('right isosceles and sliver triangles: hand-calculated values', () => {
  const r = triQuality(P([[0, 0, 0], [1, 0, 0], [0, 1, 0]]), 0, 1, 2);
  close(r.minAngle, 45, 1e-9); close(r.maxAngle, 90, 1e-9);
  close(r.skewness, 0.25, 1e-12); // max((90 − 60)/120, (60 − 45)/60)
  close(r.aspectRatio, Math.SQRT2 * (2 + Math.SQRT2) / (4 * S3 * 0.5), 1e-12);
  close(r.scaledJacobian, 2 / S3 * Math.SQRT1_2, 1e-12); // sin 45° · 2/√3
  const s = triQuality(P([[0, 0, 0], [10, 0, 0], [5, 0.1, 0]]), 0, 1, 2);
  assert.ok(s.minAngle < 1.2 && s.maxAngle > 177 && s.skewness > 0.97 && s.aspectRatio > 50);
  const d = triQuality(P([[0, 0, 0], [1, 0, 0], [2, 0, 0]]), 0, 1, 2);
  assert.equal(d.aspectRatio, Infinity); assert.equal(d.scaledJacobian, 0); assert.equal(d.skewness, 1);
});

test('square is ideal; rectangle, parallelogram and warped quadrilaterals', () => {
  const sq = quadQuality(P([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]), 0, 1, 2, 3);
  close(sq.aspectRatio, 1, 1e-12); close(sq.minAngle, 90, 1e-9); close(sq.maxAngle, 90, 1e-9);
  close(sq.skewness, 0, 1e-12); close(sq.scaledJacobian, 1, 1e-12); close(sq.warpage, 0, 1e-9); close(sq.area, 1, 1e-12);
  const rect = quadQuality(P([[0, 0, 0], [4, 0, 0], [4, 1, 0], [0, 1, 0]]), 0, 1, 2, 3);
  close(rect.aspectRatio, 4, 1e-12); close(rect.scaledJacobian, 1, 1e-12);
  const par = quadQuality(P([[0, 0, 0], [1, 0, 0], [1.5, S3 / 2, 0], [0.5, S3 / 2, 0]]), 0, 1, 2, 3);
  close(par.minAngle, 60, 1e-9); close(par.maxAngle, 120, 1e-9); close(par.skewness, 1 / 3, 1e-9); close(par.scaledJacobian, S3 / 2, 1e-12);
  const warped = quadQuality(P([[0, 0, 0], [1, 0, 0], [1, 1, 0.3], [0, 1, 0]]), 0, 1, 2, 3);
  assert.ok(warped.warpage > 10);
  // re-entrant (arrow-head) quadrilateral has a negative corner Jacobian
  assert.ok(quadQuality(P([[0, 0, 0], [1, 0, 0], [0.2, 0.2, 0], [0, 1, 0]]), 0, 1, 2, 3).scaledJacobian < 0);
});

test('regular tetrahedron and cube are ideal; prism and pyramid with equilateral faces too', () => {
  const tet = P([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]);
  const q = cellQuality(tet, [0, 1, 2, 3], 0, 'tet4');
  assert.ok(q.volume > 0);
  close(q.aspectRatio, 1, 1e-12); close(q.minAngle, Math.acos(1 / 3) * 180 / Math.PI, 1e-9); close(q.maxAngle, q.minAngle, 1e-9);
  close(q.skewness, 0, 1e-12); close(q.scaledJacobian, 1, 1e-12); close(q.volume, 8 / 3, 1e-12);
  const hex = P([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]);
  const h = cellQuality(hex, [0, 1, 2, 3, 4, 5, 6, 7], 0, 'hex8');
  close(h.aspectRatio, 1, 1e-12); close(h.minAngle, 90, 1e-9); close(h.maxAngle, 90, 1e-9); close(h.skewness, 0, 1e-12); close(h.scaledJacobian, 1, 1e-12);
  const w = cellQuality(P([[0, 0, 0], [1, 0, 0], [0.5, S3 / 2, 0], [0, 0, 1], [1, 0, 1], [0.5, S3 / 2, 1]]), [0, 1, 2, 3, 4, 5], 0, 'wedge6');
  close(w.scaledJacobian, 1, 1e-12); close(w.skewness, 0, 1e-9); close(w.aspectRatio, 1, 1e-12);
  const py = cellQuality(P([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, Math.SQRT1_2]]), [0, 1, 2, 3, 4], 0, 'pyr5');
  close(py.scaledJacobian, 1, 1e-12); close(py.skewness, 0, 1e-9);
  assert.throws(() => cellQuality(hex, [0, 1, 2, 3, 4, 5, 6, 7], 0, 'hex27'), RangeError);
});

test('distorted cells: corner tetrahedron, sheared and inverted hexahedra', () => {
  // corner of a cube: three mutually perpendicular unit edges
  const c = cellQuality(P([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]), [0, 1, 2, 3], 0, 'tet4');
  close(c.volume, 1 / 6, 1e-12);
  // R = √3/2, r = 1/(3 + √3): radius ratio R / 3r
  close(c.aspectRatio, (S3 / 2) * (3 + S3) / 3, 1e-12);
  close(c.skewness, 1 - (1 / 6) / (8 * S3 / 27 * (S3 / 2) ** 3), 1e-12); // = 0.5
  close(c.minAngle, Math.acos(1 / S3) * 180 / Math.PI, 1e-9); close(c.maxAngle, 90, 1e-9);
  close(c.scaledJacobian, Math.SQRT2 / 2, 1e-12); // √2 · 6V / (√2·√2·1) at the worst corner
  // a flat tetrahedron
  assert.ok(cellQuality(P([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.3, 0.3, 0.001]]), [0, 1, 2, 3], 0, 'tet4').skewness > 0.99);
  // hexahedron sheared by 45° in x
  const sh = cellQuality(P([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 1], [2, 0, 1], [2, 1, 1], [1, 1, 1]]), [0, 1, 2, 3, 4, 5, 6, 7], 0, 'hex8');
  close(sh.scaledJacobian, Math.SQRT1_2, 1e-12); close(sh.minAngle, 45, 1e-9); close(sh.maxAngle, 135, 1e-9); close(sh.skewness, 0.5, 1e-9);
  close(sh.aspectRatio, Math.SQRT2, 1e-12); close(sh.volume, 1, 1e-12);
  const inv = cellQuality(P([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]), [4, 5, 6, 7, 0, 1, 2, 3], 0, 'hex8');
  assert.ok(inv.volume < 0 && inv.scaledJacobian <= -1 + 1e-12);
});

test('assess: a regular hexahedral block grades A with a single-bin histogram', () => {
  const r = assess(hexBlock(4));
  assert.equal(r.grade, 'A');
  assert.equal(r.summary.elements, 64); assert.deepEqual(r.summary.byType, { hex8: 64 });
  assert.equal(r.perElement.scaledJacobian.length, 64);
  assert.ok(r.perElement.rating.every((v) => v === 0));
  assert.equal(r.histogram.metric, 'scaledJacobian');
  assert.equal(r.histogram.counts[r.histogram.counts.length - 1], 64);
  assert.equal(r.histogram.edges.length, r.histogram.counts.length + 1);
  close(r.summary.size.mean, 0.25, 1e-6); close(r.summary.size.ratio, 1, 1e-6);
  assert.equal(r.summary.poor + r.summary.marginal + r.summary.inverted, 0);
  assert.match(r.diagnosis.join(' '), /All elements are within the recommended limits/);
});

test('assess: surface grading, per-element ratings and plain-language diagnosis', () => {
  const s = assess(sphere({ radius: 1, subdivisions: 3 }));
  assert.equal(s.grade, 'A');
  assert.ok(s.summary.minAngle.min > 50);
  assert.ok(s.panel && s.panel.normalsOutward);
  // one long thin face on a stretched box
  const b = box({ size: [100, 1, 1] });
  const r = assess(b, { panel: false });
  assert.equal(r.grade, 'D'); assert.equal(r.panel, undefined);
  assert.ok(r.summary.poor >= 8);
  assert.ok(r.perElement.rating.some((v) => v === 2));
  assert.match(r.diagnosis.join(' '), /corner angle below the limit/);
  assert.match(r.diagnosis.join(' '), /too stretched/);
  assert.ok(r.summary.worst.value < 0.05);
  // limits can be overridden per family
  const loose = assess(b, { panel: false, limits: { tri: { aspectRatio: [1e3, 1e4], minAngle: [0.1, 0.01], maxAngle: [179.9, 179.99], skewness: [0.999, 0.9999], scaledJacobian: [0.001, 0.0001] } } });
  assert.equal(loose.grade, 'A');
  assert.deepEqual(QUALITY_LIMITS.tri.minAngle, [20, 10]);
  // an inverted cell forces grade D
  const blk = hexBlock(2);
  const conn = Uint32Array.from(blk.cells[0].conn);
  for (let k = 0; k < 4; k++) { const t = conn[k]; conn[k] = conn[k + 4]; conn[k + 4] = t; }
  const bad = assess({ ...blk, cells: [{ type: 'hex8', conn }] });
  assert.equal(bad.grade, 'D'); assert.equal(bad.summary.inverted, 1); assert.equal(bad.summary.worst.element, 0);
  assert.match(bad.diagnosis.join(' '), /inverted/);
  // quadrilateral panels supplied alongside triangles
  const q = assess({ kind: 'surface', positions: Float32Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]), indices: new Uint32Array(0), quads: Uint32Array.from([0, 1, 2, 3]) }, { panel: false });
  assert.deepEqual(q.summary.byType, { quad4: 1 }); assert.equal(q.grade, 'A');
});

test('assess: frames, polylines, grids and kinds without elements', () => {
  const f = assess({ kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 0, 0, 10.05]), members: [{ a: 0, b: 1 }, { a: 1, b: 2 }] });
  assert.equal(f.grade, 'C'); close(f.summary.size.ratio, 200, 1e-2);
  assert.match(f.diagnosis.join(' '), /conditioning/);
  assert.equal(assess({ kind: 'polyline', points: Float64Array.from([0, 0, 0, 1, 0, 0, 2, 0, 0]) }).grade, 'A');
  assert.equal(assess({ kind: 'grid', nx: 10, ny: 10, dx: 5, dy: 25, x0: 0, y0: 0, z: new Float32Array(100) }).grade, 'B');
  assert.equal(assess({ kind: 'points', positions: new Float32Array(9) }).grade, '-');
  assert.match(assess({ kind: 'section', loops: [] }).diagnosis[0], /triangulate/);
  assert.throws(() => assess({ kind: 'hologram' }), RangeError);
});

test('panels per wavelength: λ = 7h and the deep- and finite-depth dispersion relation', () => {
  const d = panelFrequencyLimit(2);
  close(d.minWavelength, 14, 1e-12);
  // deep water: T = sqrt(2πλ/g)
  close(d.minPeriod, Math.sqrt(2 * Math.PI * 14 / 9.80665), 1e-12);
  close(d.maxOmega * d.minPeriod, 2 * Math.PI, 1e-12);
  close(d.maxFrequencyHz, 1 / d.minPeriod, 1e-12);
  // shallow water (kd → 0): ω = k sqrt(g d)
  const s = panelFrequencyLimit(100, { depth: 2, panelsPerWavelength: 10 });
  rel(s.maxOmega, 2 * Math.PI / 1000 * Math.sqrt(9.80665 * 2), 1e-4);
  assert.ok(panelFrequencyLimit(2, { depth: 10 }).maxOmega < d.maxOmega);
  assert.throws(() => panelFrequencyLimit(0), RangeError);
});

test('panel mesh check: normals, closure, wetted panels and adequacy for a wave period', () => {
  const m = barge({ length: 100, beam: 30, depth: 8, draft: 5, elementSize: 4 });
  const p = panelCheck(m, { minPeriod: 4 });
  assert.ok(p.closed && p.normalsOutward && p.normalsConsistent);
  assert.equal(p.reversedPanels, 0);
  assert.equal(p.panels, m.indices.length / 3);
  assert.ok(p.panelsAboveWaterline > 0 && p.wettedPanels + p.panelsAboveWaterline === p.panels);
  rel(p.panelSizeP95, Math.hypot(4, 4), 0.05); // diagonal of a 4 m square cell
  close(p.minWavelength, 7 * p.panelSizeP95, 1e-12);
  // 4 s wave in deep water: λ = g T² / 2π = 24.98 m → 3.57 m panels
  rel(p.requiredPanelSize, 9.80665 * 16 / (2 * Math.PI) / 7, 1e-9);
  assert.equal(p.adequate, false);
  assert.match(p.notes.join(' '), /too large for 4 s waves/);
  assert.equal(panelCheck(m, { minPeriod: 8 }).adequate, true);
  // inside-out hull
  const inv = { ...m, indices: Uint32Array.from(m.indices) };
  for (let f = 0; f < inv.indices.length; f += 3) { const t = inv.indices[f + 1]; inv.indices[f + 1] = inv.indices[f + 2]; inv.indices[f + 2] = t; }
  const pi = panelCheck(inv);
  assert.equal(pi.normalsOutward, false); assert.equal(pi.normalsConsistent, true); assert.equal(pi.reversedPanels, pi.panels);
  // wetted hull cut at the waterline: open, which is expected
  const open = panelCheck(wigleyHull({ length: 100, beam: 10, draft: 6, closed: false }));
  assert.equal(open.closed, false); assert.equal(open.panelsAboveWaterline, 0);
  assert.match(open.notes.join(' '), /open/);
  assert.throws(() => panelCheck({ kind: 'frame' }), RangeError);
});
