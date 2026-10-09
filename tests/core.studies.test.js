import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gci3, defaults, fromSI, convergence, sensitivity, monteCarlo, advise, calibrate, validation } from '../src/core/studies.js';

const ctx = { site: {}, bus: {}, material: () => null, progress() {} };
// y = a·x² + c/n²: second-order "discretisation error" with a known exact value
const mod = {
  id: 's00-t', inputs: [
    { key: 'a', type: 'number', default: 2, min: 0.1, max: 10, role: 'material', group: 'g', label: 'a', design: true },
    { key: 'x', type: 'number', default: 3, min: 0, max: 10, role: 'bc', group: 'g', label: 'x', unit: 'mm', bind: 'k.x' },
    { key: 'n', type: 'int', default: 8, min: 1, max: 1e4, role: 'control', group: 'g', label: 'n' }],
  mesh: { key: 'n', levels: [4, 8, 16, 32], qoi: 'y' },
  compute: (v) => ({ kpis: [{ key: 'y', label: 'y', value: v.a * v.x * v.x + 5 / (v.n * v.n), limit: 30, sense: 'max' }] }),
};

test('GCI: Celik et al. (2008) Table 1 worked example', () => {
  // φ1 = 6.063, φ2 = 5.972, φ3 = 5.863, r21 = 1.5, r32 = 1.333 → p = 1.53, φext = 6.1685, GCI_fine = 2.2 %
  const r = gci3(1, 1.5, 1.5 * 1.333, 6.063, 5.972, 5.863);
  assert.ok(Math.abs(r.p - 1.53) < 0.02, `p = ${r.p}`);
  assert.ok(Math.abs(r.fExt - 6.1685) < 0.002, `fExt = ${r.fExt}`);
  assert.ok(Math.abs(r.gciFine - 0.022) < 0.001, `GCI = ${r.gciFine}`);
});
test('GCI recovers order and exact value from synthetic data', () => {
  const f = (h) => 10 + 3 * h ** 2;
  const r = gci3(0.25, 0.5, 1, f(0.25), f(0.5), f(1));
  assert.ok(Math.abs(r.p - 2) < 1e-9); assert.ok(Math.abs(r.fExt - 10) < 1e-9); assert.equal(r.behaviour, 'monotone');
});
test('defaults: bus value is converted to the field unit and marked linked', () => {
  const d = defaults(mod, { 'k.x': 0.004 });
  assert.equal(d.values.x, 4); assert.equal(d.linked.x, 'k.x');
  assert.equal(defaults(mod, { 'k.x': 0.004 }, { x: 7 }).values.x, 7);
  assert.equal(fromSI(2.5e6, 'MPa'), 2.5);
});
test('convergence study: second order, exact extrapolation, coarsest admissible level', async () => {
  const r = await convergence(mod, defaults(mod).values, ctx, { tol: 0.001 });
  assert.ok(Math.abs(r.observedOrder - 2) < 1e-6); assert.ok(Math.abs(r.fExt - 18) < 1e-9);
  assert.equal(r.selected, 32); // error 5/n²/18: 16 → 0.108 %, 32 → 0.027 %
  assert.equal((await convergence(mod, defaults(mod).values, ctx, { tol: 0.002 })).selected, 16);
});
test('sensitivity ranks x above a (quadratic vs linear)', async () => {
  const r = await sensitivity(mod, defaults(mod).values, ctx, { pct: 10 });
  const swing = (k) => { const row = r.rows.find((q) => q.key === k); return Math.abs(row.out.high.y - row.out.low.y); };
  assert.ok(swing('x') > swing('a'));
});
test('Monte Carlo mean of a linear response equals the analytic mean; seed reproducible', async () => {
  const lin = { ...mod, compute: (v) => ({ kpis: [{ key: 'y', label: 'y', value: 4 * v.a }] }) };
  const o = { n: 400, seed: 3, vars: { a: { dist: 'normal', cov: 0.1 } } };
  const r1 = await monteCarlo(lin, defaults(lin).values, ctx, o), r2 = await monteCarlo(lin, defaults(lin).values, ctx, o);
  assert.ok(Math.abs(r1.stats.y.mean - 8) < 0.01); assert.ok(Math.abs(r1.stats.y.std - 0.8) < 0.03);
  assert.equal(r1.stats.y.p90, r2.stats.y.p90);
});
test('advisor finds the limiting value of a design variable', async () => {
  const v = { ...defaults(mod).values, x: 4.2 }; // y = 35.4 > 30 → fail; a must fall to (30 − 5/64)/17.64
  const r = await advise(mod, v, ctx);
  assert.equal(r.status, 'fail');
  assert.ok(Math.abs(r.options[0].to - (30 - 5 / 64) / 17.64) < 2e-3);
});
test('calibration recovers a parameter; validation reports metrics', async () => {
  const data = [1, 2, 3, 4].map((x) => ({ set: { x }, kpi: 'y', value: 3.3 * x * x + 5 / 64 }));
  const r = await calibrate(mod, defaults(mod).values, ctx, { params: [{ key: 'a', min: 1, max: 6 }], data });
  assert.ok(Math.abs(r.params[0].fitted - 3.3) < 1e-3);
  const val = await validation(mod, { ...defaults(mod).values, a: 3.3 }, ctx, { data });
  assert.ok(val.pass); assert.ok(val.metrics.rmse < 1e-9);
  await assert.rejects(calibrate(mod, defaults(mod).values, ctx, { params: [{ key: 'a', min: 1, max: 6 }, { key: 'x', min: 0, max: 9 }], data: data.slice(0, 1) }), /identifiable/);
});
