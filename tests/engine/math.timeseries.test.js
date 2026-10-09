import { test } from 'node:test';
import assert from 'node:assert/strict';
import { acf, pacf, difference, movingAverage, expSmoothing, arima, kalmanFilter, yuleWalker } from '../../src/engine/math/timeseries.js';
import { makeRng } from '../../src/engine/math/rng.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);

function simulateArma({ phi = [], theta = [], mu = 0, sigma = 1, n, seed }) {
  const rng = makeRng(seed);
  const burn = 500;
  const x = new Array(n + burn).fill(0);
  const e = new Array(n + burn).fill(0);
  for (let t = 0; t < n + burn; t++) {
    e[t] = sigma * rng.normal();
    let v = e[t];
    for (let i = 0; i < phi.length; i++) if (t - 1 - i >= 0) v += phi[i] * x[t - 1 - i];
    for (let j = 0; j < theta.length; j++) if (t - 1 - j >= 0) v += theta[j] * e[t - 1 - j];
    x[t] = v;
  }
  return x.slice(burn).map((v) => v + mu);
}

test('acf: hand calculation for [1, 2, 3, 4, 5] (r1 = 0.4, r2 = -0.1)', () => {
  // deviations -2..2, c0 = 10; lag 1: 2 + 0 + 0 + 2 = 4; lag 2: 0 - 1 + 0 = -1
  const r = acf([1, 2, 3, 4, 5], 2);
  close(r[0], 1, 1e-15);
  close(r[1], 0.4, 1e-15);
  close(r[2], -0.1, 1e-15);
  assert.throws(() => acf([1, 2, 3], 3), RangeError);
});

test('acf / pacf of a simulated AR(1) with phi = 0.7: r_k = phi^k, PACF cuts off after lag 1', () => {
  const x = simulateArma({ phi: [0.7], n: 20000, seed: 11 });
  const r = acf(x, 5);
  for (let k = 1; k <= 5; k++) close(r[k], 0.7 ** k, 0.03, `r${k}`);
  const p = pacf(x, 5);
  close(p[1], 0.7, 0.02);
  for (let k = 2; k <= 5; k++) assert.ok(Math.abs(p[k]) < 0.03, `pacf ${k} = ${p[k]}`); // ~ 4/sqrt(n)
});

test('difference and movingAverage', () => {
  assert.deepEqual(difference([1, 4, 9, 16, 25], 1), [3, 5, 7, 9]);
  assert.deepEqual(difference([1, 4, 9, 16, 25], 2), [2, 2, 2]);
  assert.deepEqual(difference([1, 2], 0), [1, 2]);
  // a centred moving average reproduces a straight line exactly, including at the truncated ends
  const line = Array.from({ length: 9 }, (_, i) => 2 + 3 * i);
  movingAverage(line, 5).forEach((v, i) => close(v, line[i], 1e-12));
  movingAverage(line, 4).forEach((v, i) => close(v, line[i], 1e-12));
  // 2x4 centred average removes a period-4 seasonal pattern
  const seas = Array.from({ length: 16 }, (_, i) => 10 + [3, -1, -4, 2][i % 4]);
  const ma = movingAverage(seas, 4);
  for (let i = 2; i < 14; i++) close(ma[i], 10, 1e-12);
  assert.deepEqual(movingAverage([1, 2, 6], 3), [1, 3, 6]);
});

test('expSmoothing: simple smoothing recursion by hand (alpha = 0.5)', () => {
  // l0 = 3; l1 = 0.5*5 + 0.5*3 = 4; l2 = 0.5*9 + 0.5*4 = 6.5; fitted = one-step-ahead predictions
  const r = expSmoothing([3, 5, 9], { alpha: 0.5 });
  assert.deepEqual(r.fitted, [3, 3, 4]);
  close(r.level, 6.5, 1e-14);
  assert.deepEqual(r.forecast(2), [6.5, 6.5]);
});

test('expSmoothing: Holt reproduces and extrapolates an exact linear trend', () => {
  const x = Array.from({ length: 20 }, (_, i) => 5 + 2 * i);
  const r = expSmoothing(x, { alpha: 0.4, beta: 0.3 });
  r.fitted.forEach((v, i) => close(v, x[i], 1e-10));
  const f = r.forecast(3);
  close(f[0], 45, 1e-10);
  close(f[2], 49, 1e-10);
  close(r.trend, 2, 1e-12);
});

test('expSmoothing: Holt–Winters additive recovers level, trend and season of a noiseless seasonal series', () => {
  const s = [4, -2, -5, 3];
  const gen = (i) => 20 + 0.5 * i + s[i % 4];
  const x = Array.from({ length: 40 }, (_, i) => gen(i));
  const r = expSmoothing(x, { alpha: 0.3, beta: 0.1, gamma: 0.2, season: 4 });
  const f = r.forecast(8);
  for (let k = 0; k < 8; k++) close(f[k], gen(40 + k), 1e-8, `h=${k + 1}`);
  close(r.trend, 0.5, 1e-9);
  assert.ok(r.sse < 1e-16);
});

test('expSmoothing: Holt–Winters with estimated parameters forecasts a noisy seasonal series; multiplicative variant', () => {
  const rng = makeRng(8);
  const s = [6, -1, -7, 2, 3, -3];
  const gen = (i) => 50 + 0.3 * i + s[i % 6];
  const x = Array.from({ length: 120 }, (_, i) => gen(i) + 0.5 * rng.normal());
  const r = expSmoothing(x, { season: 6 });
  for (const k of ['alpha', 'beta', 'gamma']) assert.ok(r[k] >= 0 && r[k] <= 1);
  const f = r.forecast(12);
  for (let k = 0; k < 12; k++) close(f[k], gen(120 + k), 1.5, `h=${k + 1}`);
  assert.ok(r.rmse < 0.8, `rmse ${r.rmse}`);
  // multiplicative seasonality: x = (100 + 2 i) * factor
  const fac = [1.2, 0.8, 0.9, 1.1];
  const gm = (i) => (100 + 2 * i) * fac[i % 4];
  const xm = Array.from({ length: 48 }, (_, i) => gm(i));
  const rm = expSmoothing(xm, { season: 4, seasonal: 'multiplicative' });
  const fm = rm.forecast(4);
  for (let k = 0; k < 4; k++) assert.ok(Math.abs(fm[k] / gm(48 + k) - 1) < 0.02, `mult h=${k + 1}: ${fm[k]} vs ${gm(48 + k)}`);
});

test('arima: AR(2) coefficient recovery (phi = 0.6, -0.3; mean 5; sigma 2) from 8000 seeded samples', () => {
  const x = simulateArma({ phi: [0.6, -0.3], mu: 5, sigma: 2, n: 8000, seed: 21 });
  const fit = arima(x, { p: 2, d: 0, q: 0 });
  assert.ok(fit.converged);
  assert.ok(fit.stationary);
  // asymptotic se of phi ≈ sqrt((1 - phi2²)/n) = 0.0107
  close(fit.phi[0], 0.6, 0.035);
  close(fit.phi[1], -0.3, 0.035);
  close(fit.sigma2, 4, 0.2);
  close(fit.mean, 5, 0.15);
  close(fit.se[0], Math.sqrt((1 - 0.09) / 8000), 0.003);
  // Yule–Walker agrees
  const yw = yuleWalker(x, 2);
  close(yw[0], fit.phi[0], 0.01);
  close(yw[1], fit.phi[1], 0.01);
  // long-horizon forecast tends to the mean with the process variance sigma²·Σpsi² = gamma0
  const f = fit.forecast(60);
  close(f.mean[59], fit.mean, 1e-6);
  // gamma0 of AR(2): sigma² (1 - phi2) / ((1 + phi2)((1 - phi2)² - phi1²))
  const g0 = (fit.sigma2 * (1 - fit.phi[1])) / ((1 + fit.phi[1]) * ((1 - fit.phi[1]) ** 2 - fit.phi[0] ** 2));
  close(f.se[59], Math.sqrt(g0), 1e-6);
  close(f.se[0], Math.sqrt(fit.sigma2), 1e-12);
  close(f.hi95[0] - f.lo95[0], 2 * 1.959964 * Math.sqrt(fit.sigma2), 1e-5);
  // one-step forecast equals the AR recursion on the last two observations
  const n = x.length;
  close(f.mean[0], fit.mean + fit.phi[0] * (x[n - 1] - fit.mean) + fit.phi[1] * (x[n - 2] - fit.mean), 1e-10);
});

test('arima: MA(1) and ARMA(1,1) parameter recovery by conditional sum of squares', () => {
  const x = simulateArma({ theta: [0.5], n: 8000, seed: 5 });
  const ma = arima(x, { p: 0, d: 0, q: 1 });
  close(ma.theta[0], 0.5, 0.04);
  close(ma.sigma2, 1, 0.06);
  // MA(1) forecasts: beyond one step the forecast is the mean and the variance is sigma²(1 + theta²)
  const f = ma.forecast(3);
  close(f.mean[2], ma.mean, 1e-12);
  close(f.se[1], Math.sqrt(ma.sigma2 * (1 + ma.theta[0] ** 2)), 1e-12);
  const y = simulateArma({ phi: [0.7], theta: [0.4], n: 10000, seed: 6 });
  const arma = arima(y, { p: 1, d: 0, q: 1 });
  close(arma.phi[0], 0.7, 0.04);
  close(arma.theta[0], 0.4, 0.05);
  // the correct model has a lower AIC than white noise around a mean
  assert.ok(arma.aic < arima(y, { p: 0, d: 0, q: 0 }).aic);
});

test('arima: ARIMA(0,1,0) random walk — forecast is the last value and se grows as sigma·sqrt(h)', () => {
  const rng = makeRng(31);
  const x = [10];
  for (let i = 1; i < 3000; i++) x.push(x[i - 1] + 0.5 * rng.normal());
  const fit = arima(x, { p: 0, d: 1, q: 0 });
  close(fit.sigma2, 0.25, 0.02);
  const f = fit.forecast(9);
  for (let h = 1; h <= 9; h++) {
    close(f.mean[h - 1], x[x.length - 1], 1e-12);
    close(f.se[h - 1], Math.sqrt(fit.sigma2 * h), 1e-12);
  }
});

test('arima: ARIMA(1,1,0) on an integrated AR(1) recovers phi and integrates the forecast', () => {
  const w = simulateArma({ phi: [0.5], n: 6000, seed: 14 });
  const x = [0];
  for (let i = 0; i < w.length; i++) x.push(x[i] + w[i]);
  const fit = arima(x, { p: 1, d: 1, q: 0 });
  close(fit.phi[0], 0.5, 0.04);
  const f = fit.forecast(2);
  const n = x.length;
  const d1 = x[n - 1] - x[n - 2];
  close(f.mean[0], x[n - 1] + fit.phi[0] * d1, 1e-10);
  close(f.mean[1], x[n - 1] + fit.phi[0] * d1 + fit.phi[0] ** 2 * d1, 1e-10);
  // psi weights of ARIMA(1,1,0): psi1 = 1 + phi
  close(f.se[1], Math.sqrt(fit.sigma2 * (1 + (1 + fit.phi[0]) ** 2)), 1e-10);
  assert.throws(() => arima([1, 2, 3], { p: 2, d: 0, q: 2 }), RangeError);
});

test('kalmanFilter: scalar constant — the estimate is the Bayesian running mean and P = 1/(1/P0 + k/R)', () => {
  const z = [1.2, 0.8, 1.1, 0.9, 1.3];
  const R = 0.5;
  const P0 = 4;
  const r = kalmanFilter({ F: [[1]], H: [[1]], Q: [[0]], R: [[R]], x0: [0], P0: [[P0]] }, z);
  let sum = 0;
  for (let k = 0; k < z.length; k++) {
    sum += z[k];
    const prec = 1 / P0 + (k + 1) / R;
    close(r.P[k][0][0], 1 / prec, 1e-12);
    close(r.x[k][0], sum / R / prec, 1e-12);
  }
  close(r.innovations[0][0], 1.2, 1e-15);
});

test('kalmanFilter: constant-velocity tracking from noisy positions', () => {
  const dt = 1;
  const rng = makeRng(17);
  const sigmaZ = 2;
  const n = 300;
  const truthPos = (k) => 5 + 1.5 * (k + 1) * dt;
  const obs = Array.from({ length: n }, (_, k) => [truthPos(k) + sigmaZ * rng.normal()]);
  const q = 1e-6;
  const model = {
    F: [[1, dt], [0, 1]],
    H: [[1, 0]],
    Q: [[(q * dt ** 3) / 3, (q * dt ** 2) / 2], [(q * dt ** 2) / 2, q * dt]],
    R: [[sigmaZ ** 2]],
    x0: [0, 0],
    P0: [[100, 0], [0, 100]],
  };
  const r = kalmanFilter(model, obs);
  assert.equal(r.x.length, n);
  const last = r.x[n - 1];
  close(last[1], 1.5, 0.02); // velocity
  close(last[0], truthPos(n - 1), 4 * Math.sqrt(r.P[n - 1][0][0]));
  // the filter is far more accurate than the raw measurements
  let seF = 0;
  let seZ = 0;
  for (let k = 100; k < n; k++) {
    seF += (r.x[k][0] - truthPos(k)) ** 2;
    seZ += (obs[k][0] - truthPos(k)) ** 2;
  }
  assert.ok(seF < 0.1 * seZ, `filter mse ${seF / 200} vs measurement mse ${seZ / 200}`);
  // covariance stays symmetric and decreases
  close(r.P[n - 1][0][1], r.P[n - 1][1][0], 1e-12);
  assert.ok(r.P[n - 1][0][0] < r.P[10][0][0]);
  // innovations are white with variance ≈ S: normalised innovation squared averages to 1
  let nis = 0;
  for (let k = 50; k < n; k++) nis += r.innovations[k][0] ** 2 / r.S[k][0][0];
  close(nis / (n - 50), 1, 0.25);
  assert.ok(Number.isFinite(r.logLik));
  // without process noise and with a diffuse prior the filter equals batch least squares of a straight line:
  // variance of the fitted value at the last of n points = sigma² [1/n + (n - tbar)²/Sxx] = sigma² 2(2n-1)/(n(n+1))
  const exact = kalmanFilter({ ...model, Q: [[0, 0], [0, 0]], P0: [[1e8, 0], [0, 1e8]] }, obs);
  close(exact.P[n - 1][0][0], (sigmaZ ** 2 * 2 * (2 * n - 1)) / (n * (n + 1)), 1e-6);
});

test('kalmanFilter: missing observations only propagate the state; two-sensor update', () => {
  const model = { F: [[1, 1], [0, 1]], H: [[1, 0]], Q: [[0, 0], [0, 0]], R: [[1]], x0: [0, 1], P0: [[1, 0], [0, 1]] };
  const r = kalmanFilter(model, [null, [NaN], 3.1]);
  assert.deepEqual(r.x[0], [1, 1]);
  assert.deepEqual(r.x[1], [2, 1]);
  assert.ok(Number.isNaN(r.innovations[0][0]));
  // P after two predictions: [[1 + 4, 2], [2, 1]] -> wait: P1 = F P0 Fᵀ = [[2,1],[1,1]], P2 = [[5,2],[2,1]]
  assert.deepEqual(r.P[1], [[5, 2], [2, 1]]);
  // two independent sensors of the same scalar state with R1 = 1, R2 = 4 and a diffuse prior:
  // posterior ≈ inverse-variance weighted mean (4*z1 + z2)/5
  const two = kalmanFilter({ F: [[1]], H: [[1], [1]], Q: [[0]], R: [[1, 0], [0, 4]], x0: [0], P0: [[1e9]] }, [[10, 20]]);
  close(two.x[0][0], (4 * 10 + 20) / 5, 1e-6);
  close(two.P[0][0][0], 0.8, 1e-6);
  assert.throws(() => kalmanFilter({ F: [[1]], H: [[1]], Q: [[0]], R: [[1]], x0: [0, 0], P0: [[1]] }, [1]), RangeError);
});
