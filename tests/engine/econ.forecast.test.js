import test from 'node:test';
import assert from 'node:assert/strict';
import * as F from '../../src/engine/econ/forecast.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const sd = (a) => { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); };

test('GBM: simulated terminal moments match E[S_T] = S0 e^(μT) and sd[ln S_T] = σ√T; estimation recovers σ', () => {
  const g = F.gbmSimulate({ S0: 100, mu: 0.08, sigma: 0.25, T: 2, steps: 24, paths: 20000, seed: 7 });
  assert.equal(g.t.length, 25); assert.equal(g.paths.length, 20000); assert.equal(g.paths[0][0], 100);
  const ST = g.paths.map((p) => p[24]);
  rel(mean(ST), 100 * Math.exp(0.16), 0.01);
  const lr = ST.map((v) => Math.log(v / 100));
  rel(sd(lr), 0.25 * Math.SQRT2, 0.02);
  close(mean(lr), (0.08 - 0.03125) * 2, 0.01);
  // deterministic series: prices growing 1 % per month exactly → σ = 0, μ = 12 ln 1.01
  const flat = F.gbmEstimate(Array.from({ length: 25 }, (_, i) => 50 * 1.01 ** i), 1 / 12);
  close(flat.sigma, 0, 1e-12); close(flat.mu, 12 * Math.log(1.01), 1e-12); assert.equal(flat.n, 24);
  // alternating log returns ±0.1 (12 of each): sample sd = 0.1 √(24/23) per step → σ = 0.35386 per √yr at monthly spacing
  const alt = [100]; for (let i = 0; i < 24; i++) alt.push(alt[i] * Math.exp(i % 2 ? -0.1 : 0.1));
  const e = F.gbmEstimate(alt, 1 / 12);
  close(e.sigma, 0.1 * Math.sqrt(24 / 23) * Math.sqrt(12), 1e-12);
  close(e.mu, 0.5 * e.sigma ** 2, 1e-12); close(e.seSigma, e.sigma / Math.sqrt(48), 1e-12);
  // estimate from one long simulated path
  const one = F.gbmSimulate({ S0: 100, mu: 0.05, sigma: 0.3, T: 40, steps: 4800, paths: 1, seed: 11 }).paths[0];
  rel(F.gbmEstimate(Array.from(one), 1 / 120).sigma, 0.3, 0.05);
  assert.deepEqual(F.gbmSimulate({ S0: 1, mu: 0, sigma: 0.2, T: 1, steps: 4, paths: 2, seed: 5 }).paths, F.gbmSimulate({ S0: 1, mu: 0, sigma: 0.2, T: 1, steps: 4, paths: 2, seed: 5 }).paths);
});

test('Ornstein–Uhlenbeck: stationary mean θ and variance σ²/(2κ); AR(1) estimation on an exact deterministic decay', () => {
  const o = F.ouSimulate({ x0: 10, kappa: 2, theta: 4, sigma: 1.5, dt: 0.05, steps: 400, paths: 3000, seed: 3 });
  const end = o.paths.map((p) => p[400]);
  close(mean(end), 4, 0.05); rel(sd(end), 1.5 / 2, 0.04);
  // E[x_t] = θ + (x0 − θ) e^(−κt): at t = 0.5, 4 + 6 e^(−1) = 6.2073
  close(mean(o.paths.map((p) => p[10])), 4 + 6 * Math.exp(-1), 0.05);
  // noise-free decay x_t = θ + (x0 − θ) b^t with b = e^(−κΔ): the regression is exact
  const x = Array.from({ length: 30 }, (_, i) => 5 + 3 * Math.exp(-0.8 * 0.25 * i));
  const e = F.ouEstimate(x, 0.25);
  close(e.kappa, 0.8, 1e-8); close(e.theta, 5, 1e-8); close(e.b, Math.exp(-0.2), 1e-10);
  close(e.halfLife, Math.LN2 / 0.8, 1e-8); assert.equal(e.meanReverting, true); close(e.sigma, 0, 1e-6);
  // a trending series is not mean reverting
  const t = F.ouEstimate(Array.from({ length: 30 }, (_, i) => 1.1 ** i), 1);
  assert.equal(t.meanReverting, false); assert.ok(Number.isNaN(t.kappa));
  // estimate from a long simulated path
  const p = F.ouSimulate({ x0: 4, kappa: 2, theta: 4, sigma: 1.5, dt: 0.02, steps: 40000, paths: 1, seed: 9 }).paths[0];
  const est = F.ouEstimate(Array.from(p), 0.02);
  rel(est.kappa, 2, 0.15); close(est.theta, 4, 0.1); rel(est.sigma, 1.5, 0.03);
});

test('Schwartz one-factor model: long-run median exp(μ − σ²/2κ); estimation is OU on the log price', () => {
  const s = F.schwartzSimulate({ S0: 60, kappa: 1.5, mu: Math.log(80), sigma: 0.3, dt: 1 / 12, steps: 240, paths: 2000, seed: 2 });
  const alpha = Math.log(80) - 0.09 / 3;
  close(s.longRunMedian, Math.exp(alpha), 1e-9);
  const lnEnd = s.paths.map((p) => Math.log(p[240]));
  close(mean(lnEnd), alpha, 0.02); rel(sd(lnEnd), 0.3 / Math.sqrt(3), 0.05);
  const prices = Array.from({ length: 40 }, (_, i) => Math.exp(4 + 0.5 * Math.exp(-0.6 * 0.5 * i)));
  const e = F.schwartzEstimate(prices, 0.5);
  close(e.kappa, 0.6, 1e-7); close(e.alpha, 4, 1e-7); close(e.mu, 4 + e.sigma ** 2 / 1.2, 1e-9); assert.equal(e.meanReverting, true);
});

test('Merton jump diffusion: compensated drift keeps E[S_T] = S0 e^(μT); without jumps it is GBM', () => {
  const j = F.jumpDiffusionSimulate({ S0: 100, mu: 0.05, sigma: 0.15, lambda: 2, jumpMean: -0.1, jumpStd: 0.2, T: 1, steps: 50, paths: 40000, seed: 4 });
  rel(mean(j.paths.map((p) => p[50])), 100 * Math.exp(0.05), 0.015);
  // total log variance = σ²T + λT (m² + s²) = 0.0225 + 2 × 0.05 = 0.1225
  const lr = j.paths.map((p) => Math.log(p[50] / 100));
  rel(sd(lr) ** 2, 0.1225, 0.06);
  const g = F.jumpDiffusionSimulate({ S0: 100, mu: 0.05, sigma: 0.15, lambda: 0, jumpMean: 0, jumpStd: 0, T: 1, steps: 10, paths: 20000, seed: 4 });
  rel(sd(g.paths.map((p) => Math.log(p[10] / 100))), 0.15, 0.02);
});

test('scenario cones: percentiles of paths; lognormal forward cone has the forward as its mean and P50 = F e^(−σ²t/2)', () => {
  // five constant paths 1…5: P50 = 3, mean = 3; cumulative convention P10 < P90
  const cone = F.scenarioCone([1, 2, 3, 4, 5].map((v) => [v, 2 * v]));
  close(cone.bands.p50[0], 3, 1e-12); close(cone.bands.p50[1], 6, 1e-12); close(cone.mean[1], 6, 1e-12);
  assert.ok(cone.bands.p10[0] < cone.bands.p90[0]);
  const fc = F.forwardCurveCone({ forward: [50, 55], times: [1, 4], sigma: 0.3 });
  close(fc.bands.p50[0], 50 * Math.exp(-0.045), 1e-9); close(fc.bands.p50[1], 55 * Math.exp(-0.18), 1e-9);
  // z90 = 1.281552: P90 = F exp(−σ²t/2 + 1.281552 σ √t)
  close(fc.bands.p90[1], 55 * Math.exp(-0.18 + 1.2815516 * 0.6), 1e-5);
  close(fc.bands.p10[1] * fc.bands.p90[1], fc.bands.p50[1] ** 2, 1e-8);
  // FX: covered interest parity F_t = S ((1 + rd)/(1 + rf))^t
  const fx = F.fxScenarios({ spot: 1.25, domesticRate: 0.05, foreignRate: 0.02, sigma: 0.1, years: 3 });
  const k = fx.t.indexOf(3);
  close(fx.forward[k], 1.25 * (1.05 / 1.02) ** 3, 1e-12);
  close(fx.bands.p50[k], fx.forward[k] * Math.exp(-0.015), 1e-9);
});

test('seasonal decomposition recovers a linear trend plus a zero-sum seasonal pattern exactly', () => {
  const pat = [3, -1, -2, 0];
  const x = Array.from({ length: 24 }, (_, t) => 10 + 0.5 * t + pat[t % 4]);
  const d = F.seasonalDecompose(x, 4);
  for (let k = 0; k < 4; k++) close(d.seasonalIndices[k], pat[k], 1e-10);
  close(d.trend[10], 10 + 0.5 * 10, 1e-10); close(d.remainder[10], 0, 1e-10);
  assert.equal(d.trend[0], null);
  close(d.seasonalIndices.reduce((s, v) => s + v, 0), 0, 1e-12);
  assert.throws(() => F.seasonalDecompose([1, 2, 3, 4, 5], 4), RangeError);
});

test('error metrics by hand: predictions [110, 90, 100] against [100, 100, 100]', () => {
  const m = F.errorMetrics([110, 90, 100], [100, 100, 100]);
  close(m.mape, 0.2 / 3, 1e-12); close(m.rmse, Math.sqrt(200 / 3), 1e-12); close(m.mae, 20 / 3, 1e-12); close(m.bias, 0, 1e-12);
  assert.equal(m.n, 3);
  close(F.errorMetrics([1, 2], [0, 4]).mape, 0.5, 1e-12); // the zero observation is skipped
});

test('Holt and Holt–Winters continue an exact linear trend and an exact additive season', () => {
  const lin = Array.from({ length: 20 }, (_, t) => 5 + 2 * t);
  const h = F.holtWintersForecast({ x: lin, horizon: 3 });
  for (let k = 0; k < 3; k++) close(h.forecast[k], 5 + 2 * (20 + k), 1e-3);
  assert.ok(h.rmse < 1e-3);
  const pat = [4, -4, 2, -2];
  const sea = Array.from({ length: 32 }, (_, t) => 50 + t + pat[t % 4]);
  const w = F.holtWintersForecast({ x: sea, horizon: 4, season: 4, seasonal: 'additive' });
  for (let k = 0; k < 4; k++) close(w.forecast[k], 50 + 32 + k + pat[k], 0.3);
});

test('ARIMA forecast: with d = 1 the trend continues only when the drift is estimated (includeMean)', () => {
  // x_t = 10 + 2t + small bounded wiggle: differences have mean ≈ 2
  const x = Array.from({ length: 40 }, (_, t) => 10 + 2 * t + Math.sin(1.7 * t));
  const flat = F.arimaForecast({ x, p: 1, d: 1, q: 0, horizon: 10 });
  const drift = F.arimaForecast({ x, p: 1, d: 1, q: 0, horizon: 10, includeMean: true });
  close(flat.drift, 0, 1e-12);
  close(drift.drift, 2, 0.05);
  // ten steps ahead the true level is 10 + 2 × 49 = 108; the drift model follows it, the default flattens out
  close(drift.mean[9], 108, 1.5);
  assert.ok(flat.mean[9] < 95, `default d = 1 forecast flattens: ${flat.mean[9]}`);
  close(drift.mean[9] - drift.mean[8], 2, 0.1);
  assert.ok(drift.lo95[9] < drift.mean[9] && drift.hi95[9] > drift.mean[9]);
  // an exact random walk with drift (constant increments 3): ARIMA(0,1,0) with mean reproduces the line
  const line = Array.from({ length: 12 }, (_, t) => 7 + 3 * t);
  const l = F.arimaForecast({ x: line, p: 0, d: 1, q: 0, horizon: 4, includeMean: true });
  for (let k = 0; k < 4; k++) close(l.mean[k], 7 + 3 * (12 + k), 1e-8);
  close(l.drift, 3, 1e-10);
  // d = 0 keeps its default (mean estimated): a stationary AR(1) forecast reverts to the series mean
  const st = Array.from({ length: 60 }, (_, t) => 20 + 3 * Math.sin(0.9 * t));
  const s0 = F.arimaForecast({ x: st, p: 1, d: 0, q: 0, horizon: 40 });
  close(s0.drift, mean(st), 0.5);
  const off = F.arimaForecast({ x: st, p: 1, d: 0, q: 0, horizon: 3, includeMean: false });
  assert.equal(off.drift, 0);
  // seasonal differencing with drift: season 4, level rising 1 per step → seasonal difference 4 exactly
  const pat = [3, -1, -2, 0];
  const sx = Array.from({ length: 28 }, (_, t) => 10 + t + pat[t % 4]);
  const sf = F.arimaForecast({ x: sx, p: 0, d: 0, q: 0, horizon: 4, season: 4, D: 1 });
  for (let k = 0; k < 4; k++) close(sf.mean[k], 10 + 28 + k + pat[k], 1e-8);
  assert.throws(() => F.arimaForecast({ x, horizon: 2, includeMean: 'yes' }), RangeError);
  assert.throws(() => F.arimaForecast({ x, horizon: 2, season: 4, D: 2 }), RangeError);
});

test('back-test: a perfect forecaster has zero error and skill 1 against the naive last-value forecast', () => {
  const x = Array.from({ length: 12 }, (_, t) => 3 * t);
  const b = F.backtest({ x, holdout: 3, forecast: (train, h) => Array.from({ length: h }, (_, k) => 3 * (train.length + k)) });
  close(b.metrics.rmse, 0, 1e-12); close(b.skill, 1, 1e-12);
  // naive forecast 24 against 27, 30, 33: RMSE = √((9 + 36 + 81)/3) = 6.4807
  close(b.naive.rmse, Math.sqrt(42), 1e-12);
  assert.deepEqual(b.actual, [27, 30, 33]);
  assert.throws(() => F.backtest({ x, holdout: 3, forecast: () => [1] }), RangeError);
});

test('local-level Kalman filter: steady-state gain for q = r is the golden ratio conjugate; forecast variance grows by q', () => {
  // P_pred = P + q, K = P_pred/(P_pred + r), P = (1 − K) P_pred; with q = r = 1 the fixed point is K = (√5 − 1)/2
  const y = Array.from({ length: 60 }, (_, t) => (t % 2 ? 1 : -1));
  const k = F.kalmanLocalLevel({ y, q: 1, r: 1, horizon: 3 });
  const K = (Math.sqrt(5) - 1) / 2;
  close(k.variance[59], K * 1, 1e-9); // P = K r
  close(k.forecast[0], k.level[59], 1e-12); close(k.forecast[2], k.level[59], 1e-12);
  close(k.forecastStd[1] ** 2 - k.forecastStd[0] ** 2, 1, 1e-9);
  // constant observations: the level converges to the constant; a missing value is bridged
  const c = F.kalmanLocalLevel({ y: [5, 5, null, 5, 5, 5, 5, 5], q: 0.01, r: 1 });
  close(c.level[7], 5, 1e-6);
  assert.ok(Number.isFinite(c.logLik));
});

test('VAR(1): coefficients recovered exactly from a noise-free linear recursion; forecasts iterate the recursion', () => {
  const A = [[0.5, 0.2], [-0.1, 0.7]]; const c = [1, 2];
  const rows = [[3, -1]];
  for (let t = 1; t < 30; t++) { const p = rows[t - 1]; rows.push([c[0] + A[0][0] * p[0] + A[0][1] * p[1], c[1] + A[1][0] * p[0] + A[1][1] * p[1]]); }
  const v = F.var1Fit(rows.slice(0, 12));
  for (let i = 0; i < 2; i++) { close(v.c[i], c[i], 1e-7); for (let j = 0; j < 2; j++) close(v.A[i][j], A[i][j], 1e-7); }
  const f = v.forecast(2);
  close(f[0][0], rows[12][0], 1e-6); close(f[1][1], rows[13][1], 1e-6);
  assert.equal(v.stable, true);
});

test('Student t quantiles: t(0.975; 10) = 2.2281, t(0.95; 5) = 2.0150, large dof → normal 1.9600', () => {
  close(F.tQuantile(0.975, 10), 2.2281, 1e-4); close(F.tQuantile(0.95, 5), 2.0150, 1e-4);
  close(F.tQuantile(0.975, 1e6), 1.95996, 1e-4); close(F.tQuantile(0.5, 7), 0, 1e-6);
  close(F.tQuantile(0.025, 10), -2.2281, 1e-4);
  // dof = 1 is the Cauchy distribution: quantile tan(π(p − ½)) → t(0.75; 1) = 1
  close(F.tQuantile(0.75, 1), 1, 1e-6);
});

test('regression cost model: exact line and exact power law; prediction interval wider than confidence interval', () => {
  // textbook data x = 1…5, y = 2.2, 4.1, 6.2, 7.9, 10.1: slope = 19.6/10 = 1.96, intercept = 6.1 − 1.96 × 3 = 0.22
  const m = F.regressionCostModel({ X: [[1], [2], [3], [4], [5]], y: [2.2, 4.1, 6.2, 7.9, 10.1] });
  close(m.beta[1], 1.96, 1e-10); close(m.beta[0], 0.22, 1e-10); assert.equal(m.dof, 3);
  // SSE = Σy² − b0 Σy − b1 Σxy = 224.91 − 0.22 × 30.5 − 1.96 × 111.1 = 0.444 → σ = √(0.444/3) = 0.12166... check via r²
  const sst = [2.2, 4.1, 6.2, 7.9, 10.1].reduce((s, v) => s + (v - 6.1) ** 2, 0);
  const sse = [1, 2, 3, 4, 5].reduce((s, x, i) => s + ([2.2, 4.1, 6.2, 7.9, 10.1][i] - 0.22 - 1.96 * x) ** 2, 0);
  close(m.r2, 1 - sse / sst, 1e-10); close(m.sigma, Math.sqrt(sse / 3), 1e-10);
  close(m.se[1], Math.sqrt(sse / 3 / 10), 1e-10);
  const p = m.predict([3]);
  close(p.mean, 6.1, 1e-10);
  // at the mean of x: CI half-width = t(0.95; 3) σ/√5, PI half-width = t σ √(1 + 1/5)
  const tq = F.tQuantile(0.95, 3);
  close(p.ciHigh - p.mean, tq * m.sigma / Math.sqrt(5), 1e-8); close(p.piHigh - p.mean, tq * m.sigma * Math.sqrt(1.2), 1e-8);
  // power law C = 3 x^0.7 fitted in logs is exact
  const xs = [1, 2, 4, 8, 16, 32];
  const pl = F.regressionCostModel({ X: xs.map((x) => [x]), y: xs.map((x) => 3 * x ** 0.7), logLog: true });
  close(pl.beta[1], 0.7, 1e-10); close(Math.exp(pl.beta[0]), 3, 1e-9); close(pl.predict([10]).mean, 3 * 10 ** 0.7, 1e-8);
  assert.throws(() => F.regressionCostModel({ X: [[1], [2]], y: [1, 2] }), RangeError);
});

test('GLM by IRLS: gaussian equals least squares; poisson and binomial recover exact log-linear and logistic means', () => {
  const g = F.glmFit({ X: [[1], [2], [3], [4], [5]], y: [2.2, 4.1, 6.2, 7.9, 10.1] });
  close(g.beta[0], 0.22, 1e-8); close(g.beta[1], 1.96, 1e-8);
  const xs = [0, 1, 2, 3, 4, 5];
  const po = F.glmFit({ X: xs.map((x) => [x]), y: xs.map((x) => Math.exp(0.3 + 0.4 * x)), family: 'poisson' });
  close(po.beta[0], 0.3, 1e-6); close(po.beta[1], 0.4, 1e-6); assert.ok(po.converged && po.deviance < 1e-8);
  close(po.predict([6]), Math.exp(2.7), 1e-4);
  const bi = F.glmFit({ X: xs.map((x) => [x]), y: xs.map((x) => 1 / (1 + Math.exp(-(-2 + 0.8 * x)))), family: 'binomial' });
  close(bi.beta[0], -2, 1e-5); close(bi.beta[1], 0.8, 1e-5); close(bi.predict([2.5]), 0.5, 1e-5);
});

test('machine-learning regressors reproduce simple relationships (GP interpolation, boosting, forest, MLP)', () => {
  const X = Array.from({ length: 25 }, (_, i) => [i / 4]);
  const y = X.map(([x]) => Math.sin(x) + 0.5 * x);
  const gp = F.gpRegression({ X, y });
  const p = gp.predict([2.1]);
  close(p.mean, Math.sin(2.1) + 1.05, 0.02); assert.ok(p.std >= 0 && p.std < 0.1);
  assert.ok(gp.predict([40]).std > 5 * p.std, 'uncertainty grows away from the data');
  const step = X.map(([x]) => (x < 3 ? 1 : 5));
  const gb = F.gradientBoosting({ X, y: step, trees: 200, learningRate: 0.2 });
  close(gb.predict([1]), 1, 0.05); close(gb.predict([5]), 5, 0.05);
  assert.ok(gb.history[gb.history.length - 1] < gb.history[0] && gb.trainRmse < 0.05);
  const rf = F.randomForest({ X, y: step, trees: 60, seed: 2 });
  close(rf.predict([1]), 1, 0.3); close(rf.predict([5.5]), 5, 0.3); assert.ok(rf.oobRmse >= 0 && rf.oobRmse < 1.5);
  const mlp = F.mlpRegression({ X, y: X.map(([x]) => 2 + 3 * x), hidden: 4, epochs: 3000, seed: 1 });
  close(mlp.predict([3]), 11, 0.5); assert.ok(mlp.trainRmse < 0.3);
});

test('Bayesian cost updating: conjugate normal posterior, estimate-class prior and hierarchical shrinkage', () => {
  // prior N(100, 20²), three quotes 110, 120, 115 with sd 10: precision 1/400 + 3/100 = 0.0325
  //   mean = (100/400 + 345/100)/0.0325 = 113.846, sd = 1/√0.0325 = 5.547
  const b = F.bayesNormalUpdate({ priorMean: 100, priorStd: 20, observations: [110, 120, 115], obsStd: 10 });
  close(b.mean, 3.7 / 0.0325, 1e-10); close(b.std, 1 / Math.sqrt(0.0325), 1e-10);
  close(b.priorWeight, 0.0025 / 0.0325, 1e-12); close(b.predictiveStd, Math.sqrt(1 / 0.0325 + 100), 1e-10);
  close(b.p50, b.mean, 1e-10); close(b.p90 - b.p10, 2 * 1.2815516 * b.std, 1e-6);
  const pr = F.priorFromEstimateClass({ estimate: 100, cls: 3 });
  rel(pr.priorStd, (pr.high - pr.low) / (2 * 1.2816), 1e-4); close(pr.priorMean, 0.5 * (pr.low + pr.high), 1e-9);
  assert.ok(pr.low < 100 && pr.high > 100);
  assert.ok(F.priorFromEstimateClass({ estimate: 100, cls: 5 }).priorStd > pr.priorStd, 'a class 5 estimate is less certain');
  // identical group means: no between-group variance, complete pooling
  const same = F.hierarchicalPooling([{ name: 'a', mean: 10, variance: 4, n: 4 }, { name: 'b', mean: 10, variance: 4, n: 4 }]);
  close(same.tau2, 0, 1e-12); close(same.groups[0].pooled, 10, 1e-12); close(same.groups[0].shrinkage, 1, 1e-12);
  // widely separated precise groups: almost no pooling
  const far = F.hierarchicalPooling([{ name: 'a', mean: 0, variance: 1, n: 100 }, { name: 'b', mean: 100, variance: 1, n: 100 }]);
  assert.ok(far.groups[0].shrinkage < 0.01 && Math.abs(far.groups[0].pooled) < 1); close(far.grandMean, 50, 1e-9);
});

test('reference-class uplift, learning curve, escalation and carbon price paths', () => {
  const tb = { probabilities: [0.5, 0.2], uplifts: [0.1, 0.3] };
  const table = Object.keys(F.defaultOptimismBias).length ? F.defaultOptimismBias : tb;
  const u = F.referenceClassUplift({ estimate: 200, acceptableOverrunProbability: 0.5 });
  close(u.upliftedEstimate, 200 * (1 + u.uplift), 1e-9); assert.equal(u.usesPlaceholderTable, true);
  assert.ok(F.referenceClassUplift({ estimate: 200, acceptableOverrunProbability: 0.2 }).uplift > u.uplift, 'less tolerance for overrun needs more uplift');
  assert.ok(table);
  // Wright's law, 15 % learning rate (85 % progress ratio): cost falls to 85 % per doubling
  const lc = F.learningCurveProjection({ c0: 100, q0: 10, volumes: [10, 20, 40, 15], learningRate: 0.15 });
  close(lc[0], 100, 1e-9); close(lc[1], 85, 1e-9); close(lc[2], 72.25, 1e-9); close(lc[3], 100 * 1.5 ** Math.log2(0.85), 1e-9);
  // index 100 → 121 over two years: CAGR 10 %; linear trend through 100, 110, 120: slope 10
  const e = F.escalationForecast({ years: [2020, 2021, 2022], values: [100, 110, 121], horizon: 2 });
  close(e.rate, 0.1, 1e-12); assert.deepEqual(e.years, [2023, 2024]); close(e.values[1], 121 * 1.21, 1e-9);
  const l = F.escalationForecast({ years: [2020, 2021, 2022], values: [100, 110, 120], horizon: 2, method: 'linear' });
  close(l.values[0], 130, 1e-9); close(l.values[1], 140, 1e-9); close(l.rate, 10 / 120, 1e-12);
  const cp = F.carbonPricePaths({ years: 4, scenarios: { low: { start: 50, growth: 0 }, high: { start: 80, growth: 0.25 } }, cap: 130 });
  assert.deepEqual(cp.low, [50, 50, 50, 50]); close(cp.high[1], 100, 1e-12); close(cp.high[2], 125, 1e-12); close(cp.high[3], 130, 1e-12);
});
