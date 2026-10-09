import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gaussianProcess, metropolis, bayesOptimise, expectedImprovement } from '../../src/engine/reliability/surrogate.js';
import { makeRng } from '../../src/engine/math/rng.js';
import { normCdf, normPdf, correlation } from '../../src/engine/math/stats.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);

test('gaussianProcess: interpolates its training points with vanishing variance (both kernels)', () => {
  const X = [0, 0.7, 1.5, 2.1, 3.0, 3.8, 4.4, 5.2, 6.0].map((v) => [v]);
  const y = X.map(([v]) => Math.sin(v) + 0.1 * v);
  for (const kernel of ['rbf', 'matern52']) {
    const gp = gaussianProcess(X, y, { kernel });
    X.forEach((x, i) => {
      const p = gp.predict(x);
      close(p.mean, y[i], 2e-4, `${kernel} mean at training point ${i}`);
      assert.ok(p.std < 2e-3, `${kernel} std at training point ${i}: ${p.std}`);
    });
    // between the points the true function is recovered and the predicted uncertainty is honest
    for (const x of [0.35, 1.1, 2.6, 4.9, 5.6]) {
      const p = gp.predict([x]);
      const truth = Math.sin(x) + 0.1 * x;
      close(p.mean, truth, kernel === 'rbf' ? 0.02 : 0.06, `${kernel} at ${x}`);
      assert.ok(Math.abs(p.mean - truth) < 4 * p.std + 1e-3, `${kernel} error within 4 std at ${x}`);
    }
    // far from the data the std grows towards the prior and the mean reverts to the constant trend
    const far = gp.predict([30]);
    assert.ok(far.std > 50 * gp.predict([2.6]).std);
    close(far.mean, gp.params.mean, 1e-3);
    assert.ok(Number.isFinite(gp.logLik));
    assert.equal(gp.params.lengthScale.length, 1);
  }
});

test('gaussianProcess: two-dimensional smooth function, vectorised prediction and scalar convenience', () => {
  const rng = makeRng(10);
  const f = (a, b) => Math.exp(-((a - 1) ** 2) / 2) * Math.cos(b) + 0.2 * a;
  const X = Array.from({ length: 60 }, () => [rng.uniform(-1, 3), rng.uniform(0, 4)]);
  const y = X.map(([a, b]) => f(a, b));
  const gp = gaussianProcess(X, y);
  assert.equal(gp.dim, 2);
  assert.equal(gp.params.lengthScale.length, 2);
  const test = Array.from({ length: 200 }, () => [rng.uniform(-0.8, 2.8), rng.uniform(0.2, 3.8)]);
  const pred = gp.predict(test);
  assert.equal(pred.length, 200);
  let mx = 0;
  for (let i = 0; i < 200; i++) mx = Math.max(mx, Math.abs(pred[i].mean - f(test[i][0], test[i][1])));
  assert.ok(mx < 0.03, `max error ${mx}`);
  close(gp.predict(X[7]).mean, y[7], 1e-4);
  // 1-D convenience: scalar inputs
  const g1 = gaussianProcess([0, 1, 2, 3], [0, 1, 4, 9]);
  close(g1.predict(2).mean, 4, 1e-3);
  assert.equal(g1.predict([0.5, 1.5]).length, 2);
});

test('gaussianProcess: marginal-likelihood fit beats a poor fixed length scale; noise smooths noisy data', () => {
  const X = Array.from({ length: 25 }, (_, i) => [i * 0.25]);
  const y = X.map(([v]) => Math.sin(2 * v));
  const opt = gaussianProcess(X, y);
  const fixed = gaussianProcess(X, y, { lengthScale: 0.05, optimise: false });
  close(fixed.params.lengthScale[0], 0.05, 1e-12);
  assert.ok(opt.logLik > fixed.logLik + 5, `logLik ${opt.logLik} vs ${fixed.logLik}`);
  // the too-short length scale reverts to the mean between points; the optimised one does not
  const truth = Math.sin(2 * 1.125);
  assert.ok(Math.abs(opt.predict([1.125]).mean - truth) < 1e-3);
  assert.ok(Math.abs(fixed.predict([1.125]).mean - truth) > 0.05);
  // noisy observations: with the known noise variance the GP regresses rather than interpolates
  const rng = makeRng(4);
  const Xn = Array.from({ length: 80 }, (_, i) => [i * 0.075]);
  const yn = Xn.map(([v]) => Math.sin(v) + 0.2 * rng.normal());
  const gn = gaussianProcess(Xn, yn, { noise: 0.04 });
  let rmsFit = 0;
  let rmsTruth = 0;
  Xn.forEach(([v], i) => {
    const m = gn.predict([v]).mean;
    rmsFit += (m - yn[i]) ** 2;
    rmsTruth += (m - Math.sin(v)) ** 2;
  });
  rmsFit = Math.sqrt(rmsFit / 80);
  rmsTruth = Math.sqrt(rmsTruth / 80);
  assert.ok(rmsFit > 0.12, `residual rms ${rmsFit} ~ noise level 0.2`);
  assert.ok(rmsTruth < 0.1, `error to the true function ${rmsTruth}`);
  assert.throws(() => gaussianProcess([[1]], [1]), RangeError);
  assert.throws(() => gaussianProcess([[1], [2]], [1]), RangeError);
  assert.throws(() => gaussianProcess([[1], [2]], [1, 2], { kernel: 'nope' }), RangeError);
});

test('metropolis: samples a correlated bivariate normal — means, standard deviations and correlation', () => {
  // target: x ~ N(1, 2²), y | x ~ N(-1 + 0.5 (x - 1), 0.75) -> sd_y = sqrt(0.25·4 + 0.75) = 1.3229, rho = 1/1.3229
  const logPost = ([x, y]) => -0.5 * ((x - 1) / 2) ** 2 - (0.5 * (y + 1 - 0.5 * (x - 1)) ** 2) / 0.75;
  const r = metropolis(logPost, [0, 0], { n: 60000, stepSizes: [1.6, 1.0], seed: 12, burn: 2000 });
  assert.equal(r.chain.length, 60000);
  assert.ok(r.acceptance > 0.2 && r.acceptance < 0.6, `acceptance ${r.acceptance}`);
  close(r.mean[0], 1, 0.1);
  close(r.mean[1], -1, 0.07);
  close(r.std[0], 2, 0.1);
  close(r.std[1], Math.sqrt(1.75), 0.07);
  close(correlation(r.chain.map((c) => c[0]), r.chain.map((c) => c[1])), 1 / Math.sqrt(1.75), 0.03);
  const again = metropolis(logPost, [0, 0], { n: 200, stepSizes: [1.6, 1.0], seed: 12 });
  assert.deepEqual(again.chain[199], metropolis(logPost, [0, 0], { n: 200, stepSizes: [1.6, 1.0], seed: 12 }).chain[199]);
});

test('metropolis: respects support constraints (exponential posterior, mean 1/rate) and Bayesian update of a normal mean', () => {
  const r = metropolis(([x]) => (x < 0 ? -Infinity : -2 * x), [1], { n: 40000, stepSizes: 1.0, seed: 3 });
  assert.ok(r.chain.every((c) => c[0] >= 0));
  close(r.mean[0], 0.5, 0.03);
  // conjugate check: prior N(0, 10²), 5 observations with sigma 1 -> posterior N(n ybar/(n + 0.01), 1/(n + 0.01))
  const obs = [2.1, 1.7, 2.4, 1.9, 2.2];
  const lp = ([m]) => -0.5 * (m / 10) ** 2 - 0.5 * obs.reduce((s, v) => s + (v - m) ** 2, 0);
  const p = metropolis(lp, [0], { n: 40000, stepSizes: 0.8, seed: 8 });
  const prec = 5 + 0.01;
  close(p.mean[0], obs.reduce((s, v) => s + v, 0) / prec, 0.02);
  close(p.std[0], 1 / Math.sqrt(prec), 0.02);
  assert.throws(() => metropolis(() => -Infinity, [0]), RangeError);
});

test('expectedImprovement: closed form', () => {
  // mean = fBest: EI = s·phi(0)
  close(expectedImprovement(2, 0.5, 2), 0.5 * normPdf(0), 1e-15);
  close(expectedImprovement(1, 0.5, 2), 1 * normCdf(2) + 0.5 * normPdf(2), 1e-15);
  assert.equal(expectedImprovement(3, 0, 2), 0);
  assert.equal(expectedImprovement(1, 0, 2), 1);
});

test('bayesOptimise: Forrester et al. (2008) function (6x-2)² sin(12x-4), minimum f = -6.0207 at x = 0.7572', () => {
  const f = ([x]) => (6 * x - 2) ** 2 * Math.sin(12 * x - 4);
  const r = bayesOptimise(f, [[0, 1]], { nInit: 5, nIter: 15, seed: 1 });
  close(r.x[0], 0.7572, 0.005);
  close(r.f, -6.0207, 0.01);
  assert.equal(r.history.length, 20);
  for (let i = 1; i < r.history.length; i++) assert.ok(r.history[i] <= r.history[i - 1]);
  assert.equal(r.X.length, 20);
  assert.deepEqual(bayesOptimise(f, [[0, 1]], { nInit: 5, nIter: 15, seed: 1 }).x, r.x);
});

test('bayesOptimise: Branin function, global minimum 0.397887, within 45 evaluations', () => {
  const branin = ([x, y]) => {
    const b = 5.1 / (4 * Math.PI ** 2);
    const c = 5 / Math.PI;
    return (y - b * x * x + c * x - 6) ** 2 + 10 * (1 - 1 / (8 * Math.PI)) * Math.cos(x) + 10;
  };
  const bounds = [[-5, 10], [0, 15]];
  const r = bayesOptimise(branin, bounds, { nInit: 10, nIter: 35, seed: 2 });
  assert.ok(r.f < 0.397887 + 0.05, `best ${r.f}`);
  assert.ok(r.x[0] >= -5 && r.x[0] <= 10 && r.x[1] >= 0 && r.x[1] <= 15);
  for (const p of r.X) assert.ok(p[0] >= -5 && p[0] <= 10 && p[1] >= 0 && p[1] <= 15);
  assert.throws(() => bayesOptimise(branin, [[0, Infinity], [0, 1]]), RangeError);
});
