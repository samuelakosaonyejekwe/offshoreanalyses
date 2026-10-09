import test from 'node:test';
import assert from 'node:assert/strict';
import { tzClay, tzSand, qz, hyperbolicCurve, scaleCurve, sampleCurve, API_TZ_CLAY, API_QZ } from '../../src/engine/geotech/tz.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('API clay t–z: tabulated points, residual plateau, odd symmetry', () => {
  const D = 1.5, tmax = 80e3;
  const c = tzClay({ tmax, D, residualRatio: 0.8 });
  const pts = [[0.0016, 0.30], [0.0031, 0.50], [0.0057, 0.75], [0.0080, 0.90], [0.0100, 1.00], [0.0200, 0.80], [0.05, 0.80]];
  for (const [r, f] of pts) close(c.p(r * D), f * tmax);
  close(c.p(0.0008 * D), 0.15 * tmax); // linear from the origin
  close(c.p(0.015 * D), 0.9 * tmax); // midway down the softening branch
  close(c.p(-0.0031 * D), -0.5 * tmax);
  close(c.k(0.001 * D), (0.3 * tmax) / (0.0016 * D));
  assert.ok(c.k(0.015 * D) < 0);
  assert.equal(c.k(0.03 * D), 0);
  assert.deepEqual(API_TZ_CLAY.zOverD, [0, 0.0016, 0.0031, 0.0057, 0.0080, 0.0100, 0.0200]);
  assert.throws(() => tzClay({ tmax, D, residualRatio: 1.2 }), RangeError);
});

test('API sand t–z: linear to tmax at 2.5 mm', () => {
  const c = tzSand({ tmax: 60e3 });
  close(c.p(0.00125), 30e3); close(c.p(0.0025), 60e3); close(c.p(0.02), 60e3); close(c.p(-0.02), -60e3);
  close(tzSand({ tmax: 60e3, zPeak: 0.00254 }).p(0.00127), 30e3);
});

test('API Q–z: tabulated points, full mobilisation at 0.1 D, no uplift resistance', () => {
  const D = 2, Qp = 10e6;
  const c = qz({ Qp, D });
  const pts = [[0.002, 0.25], [0.013, 0.50], [0.042, 0.75], [0.073, 0.90], [0.100, 1.00], [0.3, 1.0]];
  for (const [r, f] of pts) close(c.p(r * D), f * Qp);
  assert.equal(c.p(-0.01), 0); assert.equal(c.k(-0.01), 0);
  assert.deepEqual(API_QZ.qOverQp, [0, 0.25, 0.50, 0.75, 0.90, 1.00]);
});

test('hyperbolic curve: initial stiffness, asymptote, Rf cap; scaling and sampling', () => {
  const h = hyperbolicCurve({ pmax: 100e3, ki: 50e6 });
  close(h.k(0), 50e6);
  close(h.p(0.002), 0.002 / (1 / 50e6 + 0.002 / 100e3)); // = 50 kPa: half at z = pmax/ki
  close(h.p(0.002), 50e3);
  assert.ok(h.p(10) < 100e3 && h.p(10) > 99.9e3);
  const hr = hyperbolicCurve({ pmax: 100e3, ki: 50e6, Rf: 0.9 });
  close(hr.p(1), 100e3); assert.equal(hr.k(1), 0);
  assert.equal(hyperbolicCurve({ pmax: 1e6, ki: 1e8, odd: false }).p(-0.01), 0);
  const s = scaleCurve(h, 0.7);
  close(s.p(0.002), 35e3); close(s.pu, 70e3);
  const smp = sampleCurve(h, { zMax: 0.01, n: 5 });
  assert.deepEqual(smp.z.map((v) => +v.toFixed(4)), [0, 0.0025, 0.005, 0.0075, 0.01]);
});
