import { test } from 'node:test';
import assert from 'node:assert/strict';
import { monteCarlo, latinHypercube, propagate, sobolIndices, tornado, natafRho, corrMatrix, scrambledHalton } from '../../src/engine/reliability/sampling.js';
import { dist, normCdf, correlation, rankCorrelation, mean, std } from '../../src/engine/math/stats.js';
import { makeRng } from '../../src/engine/math/rng.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);

const RS = { R: { dist: 'normal', mean: 10, std: 1.5 }, S: { dist: 'normal', mean: 6, std: 1.2 } };

test('monteCarlo: linear normal R - S, pf = Phi(-(muR - muS)/sqrt(sR² + sS²)) within the reported interval', () => {
  const beta = 4 / Math.sqrt(1.5 ** 2 + 1.2 ** 2); // 2.0825
  const pfExact = normCdf(-beta); // 0.01865
  const r = monteCarlo((x) => x.R - x.S, RS, { n: 200000, seed: 4 });
  assert.equal(r.samples.length, 200000);
  assert.equal(r.inputs.R.length, 200000);
  close(r.mean, 4, 0.02);
  close(r.std, Math.sqrt(1.5 ** 2 + 1.2 ** 2), 0.01);
  assert.ok(r.ciPf[0] <= pfExact && pfExact <= r.ciPf[1], `pf ${r.pf}, ci ${r.ciPf}, exact ${pfExact}`);
  close(r.pf, pfExact, 4 * Math.sqrt((pfExact * (1 - pfExact)) / 200000));
  close(r.beta, beta, 0.03);
  close(r.covPf, Math.sqrt((1 - r.pf) / (200000 * r.pf)), 1e-12);
  // deterministic for a seed, different for another
  const again = monteCarlo((x) => x.R - x.S, RS, { n: 2000, seed: 4 });
  const again2 = monteCarlo((x) => x.R - x.S, RS, { n: 2000, seed: 4 });
  assert.equal(again.pf, again2.pf);
  assert.deepEqual(Array.from(again.samples.slice(0, 5)), Array.from(again2.samples.slice(0, 5)));
  assert.throws(() => monteCarlo((x) => x.R, RS, { n: 0 }), RangeError);
  assert.throws(() => monteCarlo((x) => x.R, { R: { dist: 'normal', mean: 1, std: -1 } }), /vars\.R/);
});

test('monteCarlo: correlated normals — var(R - S) = sR² + sS² - 2 rho sR sS', () => {
  const rho = 0.6;
  const sd = Math.sqrt(1.5 ** 2 + 1.2 ** 2 - 2 * rho * 1.5 * 1.2);
  const r = monteCarlo((x) => x.R - x.S, RS, { n: 100000, seed: 2, corr: { 'R,S': rho } });
  close(r.std, sd, 0.012);
  close(correlation(r.inputs.R, r.inputs.S), rho, 0.01);
  const pfExact = normCdf(-4 / sd);
  close(r.pf, pfExact, 5 * Math.sqrt(pfExact / 100000));
});

test('corrMatrix: accepted input forms give the same matrix', () => {
  const names = ['a', 'b', 'c'];
  const ref = [[1, 0.3, 0], [0.3, 1, -0.2], [0, -0.2, 1]];
  assert.deepEqual(corrMatrix(names, { 'a,b': 0.3, 'b|c': -0.2 }), ref);
  assert.deepEqual(corrMatrix(names, { a: { b: 0.3 }, c: { b: -0.2 } }), ref);
  assert.deepEqual(corrMatrix(names, [['a', 'b', 0.3], { a: 'b', b: 'c', rho: -0.2 }]), ref);
  assert.deepEqual(corrMatrix(names, ref), ref);
  assert.equal(corrMatrix(names, undefined), null);
  assert.throws(() => corrMatrix(names, { 'a,z': 0.1 }), RangeError);
  assert.throws(() => corrMatrix(names, { 'a,b': 1.5 }), RangeError);
});

test('natafRho: closed forms — uniform pair rho0 = 2 sin(pi rho / 6); lognormal pair; normal–lognormal', () => {
  const u1 = dist({ dist: 'uniform', min: 0, max: 1 });
  const u2 = dist({ dist: 'uniform', min: -3, max: 5 });
  for (const rho of [0.3, 0.7, -0.5]) close(natafRho(u1, u2, rho), 2 * Math.sin((Math.PI * rho) / 6), 2e-5, `uniform rho=${rho}`);
  const l1 = dist({ dist: 'lognormal', mean: 10, cov: 0.3 });
  const l2 = dist({ dist: 'lognormal', mean: 4, cov: 0.5 });
  close(natafRho(l1, l2, 0.4), Math.log(1 + 0.4 * 0.3 * 0.5) / Math.sqrt(Math.log(1.09) * Math.log(1.25)), 1e-12);
  const nn = dist({ dist: 'normal', mean: 0, std: 1 });
  close(natafRho(nn, l2, 0.4), (0.4 * 0.5) / Math.sqrt(Math.log(1.25)), 1e-12);
  assert.equal(natafRho(nn, nn, 0.4), 0.4);
  // numerical route, exponential pair: Liu & Der Kiureghian give rho0/rho = 1.229 - 0.367 rho + 0.153 rho² (max error 1.5 %)
  const e = dist({ dist: 'exponential', rate: 1 });
  for (const rho of [0.3, 0.6]) {
    const f = 1.229 - 0.367 * rho + 0.153 * rho * rho;
    assert.ok(Math.abs(natafRho(e, e, rho) / (rho * f) - 1) < 0.02, `exponential rho=${rho}`);
  }
});

test('monteCarlo: Nataf transformation reproduces the target linear correlation of non-normal marginals', () => {
  const vars = { a: { dist: 'exponential', rate: 0.5 }, b: { dist: 'uniform', min: 0, max: 10 }, c: { dist: 'gumbel', loc: 3, scale: 1 } };
  const r = monteCarlo((x) => x.a + x.b + x.c, vars, { n: 100000, seed: 12, corr: { 'a,b': 0.5, 'b,c': -0.3 } });
  close(correlation(r.inputs.a, r.inputs.b), 0.5, 0.015);
  close(correlation(r.inputs.b, r.inputs.c), -0.3, 0.015);
  close(correlation(r.inputs.a, r.inputs.c), 0, 0.015);
  // marginals are untouched
  close(mean(r.inputs.a), 2, 0.03);
  close(std(r.inputs.a), 2, 0.04);
  close(mean(r.inputs.c), 3 + 0.5772156649, 0.02);
});

test('latinHypercube: exactly one sample per equiprobable stratum; mean error far below crude Monte Carlo', () => {
  const vars = { x: { dist: 'normal', mean: 5, std: 2 }, y: { dist: 'weibull', scale: 3, shape: 1.8 } };
  const n = 500;
  const s = latinHypercube(vars, n, { seed: 3 });
  for (const k of ['x', 'y']) {
    const d = dist(vars[k]);
    const seen = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      const stratum = Math.min(n - 1, Math.floor(d.cdf(s[k][i]) * n));
      seen[stratum]++;
    }
    assert.ok(seen.every((c) => c === 1), `${k}: every stratum holds exactly one sample`);
    // stratification: |mean error| << sigma/sqrt(n) = 0.089 for x
    close(mean(s[k]), d.mean, (0.1 * d.std) / Math.sqrt(n), k);
  }
  // independent pairing: small spurious correlation
  assert.ok(Math.abs(correlation(s.x, s.y)) < 0.15);
  assert.deepEqual(Array.from(latinHypercube(vars, 20, { seed: 3 }).x), Array.from(latinHypercube(vars, 20, { seed: 3 }).x));
});

test('latinHypercube: Iman–Conover induces the target rank correlation and keeps the marginals', () => {
  const vars = { hs: { dist: 'weibull', scale: 2, shape: 1.5 }, tp: { dist: 'lognormal', mean: 9, std: 2 }, u: { dist: 'uniform', min: 0, max: 1 } };
  const n = 2000;
  const s = latinHypercube(vars, n, { seed: 9, corr: { 'hs,tp': 0.7, 'tp,u': -0.4 } });
  close(rankCorrelation(s.hs, s.tp), 0.7, 0.03);
  close(rankCorrelation(s.tp, s.u), -0.4, 0.03);
  close(rankCorrelation(s.hs, s.u), 0, 0.05);
  const plain = latinHypercube(vars, n, { seed: 9 });
  for (const k of ['hs', 'tp', 'u']) {
    const d = dist(vars[k]);
    const seen = new Uint8Array(n);
    for (let i = 0; i < n; i++) seen[Math.min(n - 1, Math.floor(d.cdf(s[k][i]) * n))]++;
    assert.ok(seen.every((c) => c === 1), `${k} remains stratified`);
    close(mean(s[k]), mean(plain[k]), 1e-9 * Math.abs(mean(plain[k])) + 1e-12);
  }
});

test('propagate: linear model of normals — output mean/std/percentiles match the closed form', () => {
  const vars = { a: { dist: 'normal', mean: 2, std: 0.5 }, b: { dist: 'normal', mean: -1, std: 1 } };
  const model = (x) => ({ y: 3 + 2 * x.a - x.b, twice: 2 * x.a });
  for (const method of ['lhs', 'mc']) {
    const n = method === 'lhs' ? 4000 : 40000;
    const r = propagate(model, vars, { n, seed: 5, method });
    assert.equal(r.outputs.y.length, n);
    assert.equal(r.inputs.a.length, n);
    const sd = Math.sqrt(4 * 0.25 + 1);
    close(r.stats.y.mean, 8, 0.03, method);
    close(r.stats.y.std, sd, 0.03, method);
    close(r.stats.y.p50, 8, 0.05, method);
    close(r.stats.y.p10, 8 - 1.2815515655 * sd, 0.07, method);
    close(r.stats.y.p90, 8 + 1.2815515655 * sd, 0.07, method);
    assert.ok(r.stats.y.min < r.stats.y.p10 && r.stats.y.max > r.stats.y.p90);
    close(r.stats.twice.mean, 4, 0.02, method);
  }
  // a model returning a bare number is stored as `value`
  assert.ok('value' in propagate((x) => x.a, vars, { n: 10 }).stats);
  assert.throws(() => propagate(model, vars, { method: 'xyz' }), RangeError);
});

test('scrambledHalton: uniform one-dimensional projections and low integration error', () => {
  const pts = scrambledHalton(4096, 6, makeRng(1));
  for (const col of pts) {
    close(mean(col), 0.5, 2e-3);
    // integral of x² over [0,1] = 1/3
    let s = 0;
    for (const v of col) s += v * v;
    close(s / 4096, 1 / 3, 2e-3);
  }
  // integral of the product of the 6 coordinates = 2^-6, error far below the Monte Carlo 1/sqrt(n)
  let prod = 0;
  for (let i = 0; i < 4096; i++) prod += pts[0][i] * pts[1][i] * pts[2][i] * pts[3][i] * pts[4][i] * pts[5][i];
  close(prod / 4096, 1 / 64, 6e-4);
  assert.ok(Math.abs(correlation(pts[4], pts[5])) < 0.03);
});

test('sobolIndices: Ishigami function (a = 7, b = 0.1), analytical S1 = 0.3139, S2 = 0.4424, S3 = 0, ST1 = 0.5576, ST3 = 0.2437', () => {
  const a = 7;
  const b = 0.1;
  const ishigami = (x) => Math.sin(x.x1) + a * Math.sin(x.x2) ** 2 + b * x.x3 ** 4 * Math.sin(x.x1);
  const u = { dist: 'uniform', min: -Math.PI, max: Math.PI };
  const r = sobolIndices(ishigami, { x1: u, x2: u, x3: u }, { n: 8192, seed: 1 });
  // analytical variance: a²/8 + b pi^4/5 + b² pi^8/18 + 1/2 = 13.8446
  const V = (a * a) / 8 + (b * Math.PI ** 4) / 5 + (b * b * Math.PI ** 8) / 18 + 0.5;
  close(r.variance, V, 0.15);
  close(r.mean, a / 2, 0.03);
  close(r.first.x1, 0.3139, 0.02, 'S1');
  close(r.first.x2, 0.4424, 0.02, 'S2');
  close(r.first.x3, 0, 0.02, 'S3');
  close(r.total.x1, 0.5576, 0.02, 'ST1');
  close(r.total.x2, 0.4424, 0.02, 'ST2');
  close(r.total.x3, 0.2437, 0.02, 'ST3');
  assert.equal(r.evaluations, 8192 * 5);
  // reproducible
  const again = sobolIndices(ishigami, { x1: u, x2: u, x3: u }, { n: 256, seed: 1 });
  assert.equal(again.first.x1, sobolIndices(ishigami, { x1: u, x2: u, x3: u }, { n: 256, seed: 1 }).first.x1);
});

test('sobolIndices: additive linear model — S_i = ST_i = a_i² var_i / V, deterministic input has zero index', () => {
  const vars = {
    p: { dist: 'normal', mean: 1, std: 2 },
    q: { dist: 'uniform', min: 0, max: 6 },
    r: { dist: 'lognormal', mean: 3, std: 1 },
    k: { dist: 'deterministic', value: 4 },
  };
  const model = (x) => 3 * x.p - 2 * x.q + x.r + x.k;
  const vp = 9 * 4;
  const vq = 4 * 3;
  const vr = 1;
  const V = vp + vq + vr;
  const s = sobolIndices(model, vars, { n: 4096, seed: 3 });
  close(s.first.p, vp / V, 0.015);
  close(s.first.q, vq / V, 0.015);
  close(s.first.r, vr / V, 0.01);
  close(s.total.p, vp / V, 0.015);
  close(s.total.q, vq / V, 0.015);
  assert.equal(s.first.k, 0);
  assert.equal(s.total.k, 0);
  close(s.first.p + s.first.q + s.first.r, 1, 0.02);
});

test('tornado: one-at-a-time swings sorted by magnitude', () => {
  const model = (x) => x.a * 10 + x.b - 0.1 * x.c;
  const t = tornado(model, { a: { low: 0.9, base: 1, high: 1.2 }, b: { low: 0, base: 5, high: 20 }, c: { low: 0, base: 10, high: 20 } });
  assert.deepEqual(t.map((e) => e.name), ['b', 'a', 'c']);
  // base output = 10 + 5 - 1 = 14 ; b: low 9, high 29 ; a: low 13, high 16 ; c: low 15, high 13
  close(t[0].low, 9, 1e-12);
  close(t[0].high, 29, 1e-12);
  close(t[0].swing, 20, 1e-12);
  close(t[1].swing, 3, 1e-12);
  close(t[2].low, 15, 1e-12);
  close(t[2].high, 13, 1e-12);
  close(t[2].swing, 2, 1e-12);
  close(t[0].base, 14, 1e-12);
  assert.throws(() => tornado(model, { a: { base: 1 } }), RangeError);
});
