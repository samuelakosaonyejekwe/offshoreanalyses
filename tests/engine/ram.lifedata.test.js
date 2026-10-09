import test from 'node:test';
import assert from 'node:assert/strict';
import {
  fitWeibullMle, fitExponentialMle, fitLognormalMle, medianRanks, medianRankRegression, weibullFunctions, bLife,
  bathtubHazard, renewalFunction, weibullRenewals, crowAmsaa, laplaceTest, availabilityMetrics, sparesForServiceLevel,
} from '../../src/engine/ram/lifedata.js';
import { makeRng } from '../../src/engine/math/rng.js';
import { gamma } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

function weibullSample(n, shape, scale, seed, censorAt) {
  const rng = makeRng(seed);
  const times = [];
  const censored = [];
  for (let i = 0; i < n; i++) {
    const t = scale * (-Math.log(1 - rng.next())) ** (1 / shape);
    if (censorAt !== undefined && t > censorAt) { times.push(censorAt); censored.push(true); } else { times.push(t); censored.push(false); }
  }
  return { times, censored };
}

test('Weibull MLE recovers the parameters from right-censored data (β = 2, η = 1000, 24 % suspended at 1200)', () => {
  const data = weibullSample(3000, 2, 1000, 5, 1200);
  const nCens = data.censored.filter(Boolean).length;
  close(nCens / 3000, Math.exp(-1.44), 0.1);
  const fit = fitWeibullMle(data);
  close(fit.shape, 2, 0.04);
  close(fit.scale, 1000, 0.03);
  assert.equal(fit.failures, 3000 - nCens);
  // Estimates lie within about three standard errors of the truth.
  assert.ok(Math.abs(fit.shape - 2) < 3.5 * fit.seShape, `shape ${fit.shape} ± ${fit.seShape}`);
  assert.ok(Math.abs(fit.scale - 1000) < 3.5 * fit.seScale, `scale ${fit.scale} ± ${fit.seScale}`);
  // Ignoring the censoring flag (treating suspensions as failures) biases the scale low.
  const naive = fitWeibullMle({ times: data.times });
  assert.ok(naive.scale < 0.93 * 1000);
  // The MLE satisfies the likelihood equation: η^β = Σ t^β / r.
  let s = 0;
  for (const t of data.times) s += t ** fit.shape;
  close(fit.scale ** fit.shape, s / fit.failures, 1e-9);
});

test('Weibull MLE on complete data: likelihood equations hold; heavy censoring still recovers the shape', () => {
  const data = weibullSample(500, 0.8, 200, 9);
  const fit = fitWeibullMle({ times: data.times });
  close(fit.shape, 0.8, 0.08);
  close(fit.scale, 200, 0.12);
  // Score equation for β: Σ t^β ln t / Σ t^β − 1/β − mean(ln t) = 0
  let s0 = 0; let s1 = 0; let sl = 0;
  for (const t of data.times) { s0 += t ** fit.shape; s1 += t ** fit.shape * Math.log(t); sl += Math.log(t); }
  assert.ok(Math.abs(s1 / s0 - 1 / fit.shape - sl / 500) < 1e-9);
  const heavy = weibullSample(4000, 3, 100, 21, 70); // about 70 % suspended
  const f2 = fitWeibullMle(heavy);
  close(f2.shape, 3, 0.08);
  close(f2.scale, 100, 0.05);
});

test('exponential MLE with censoring: λ = failures / total time on test', () => {
  const r = fitExponentialMle({ times: [100, 200, 300, 400, 500], censored: [false, false, true, false, true] });
  close(r.rate, 3 / 1500, 1e-14);
  close(r.mtbf, 500, 1e-12);
});

test('lognormal MLE with right censoring (EM) recovers μ and σ', () => {
  const rng = makeRng(17);
  const times = [];
  const censored = [];
  for (let i = 0; i < 4000; i++) {
    const t = Math.exp(5 + 0.6 * rng.normal());
    if (t > 200) { times.push(200); censored.push(true); } else { times.push(t); censored.push(false); }
  }
  const fit = fitLognormalMle({ times, censored });
  close(fit.mu, 5, 0.01);
  close(fit.sigma, 0.6, 0.05);
  // Complete data: closed-form MLE (mean and population standard deviation of ln t).
  const c = fitLognormalMle({ times: [10, 20, 40] });
  close(c.mu, Math.log(20), 1e-10);
  close(c.sigma, Math.sqrt((2 * Math.log(2) ** 2) / 3), 1e-8);
});

test('median ranks with Johnson adjustment for suspensions: hand calculation', () => {
  // n = 4: 10 (F), 20 (S), 30 (F), 40 (F). Adjusted ranks: 1; 1 + (5 − 1)/(1 + 2) = 2.3333; 2.3333 + (5 − 2.3333)/(1 + 1) = 3.6667
  // Benard: (1 − 0.3)/4.4 = 0.15909; (2.3333 − 0.3)/4.4 = 0.46212; (3.6667 − 0.3)/4.4 = 0.76515
  const r = medianRanks({ times: [30, 10, 40, 20], censored: [false, false, false, true] });
  assert.deepEqual(r.t, [10, 30, 40]);
  close(r.rank[1], 7 / 3, 1e-12);
  close(r.rank[2], 11 / 3, 1e-12);
  close(r.F[0], 0.7 / 4.4, 1e-12);
  close(r.F[1], (7 / 3 - 0.3) / 4.4, 1e-12);
  close(r.F[2], (11 / 3 - 0.3) / 4.4, 1e-12);
});

test('median-rank regression: exact recovery from quantile-perfect data and probability-plot coordinates', () => {
  // Times placed exactly at the Benard median-rank quantiles of Weibull(β = 1.5, η = 500) plot on a straight line.
  const n = 12;
  const times = [];
  for (let i = 1; i <= n; i++) times.push(500 * (-Math.log(1 - (i - 0.3) / (n + 0.4))) ** (1 / 1.5));
  const w = medianRankRegression({ times });
  close(w.shape, 1.5, 1e-9);
  close(w.scale, 500, 1e-9);
  close(w.r2, 1, 1e-9);
  close(w.plot.x[0], Math.log(times[0]), 1e-14);
  close(w.plot.y[0], Math.log(-Math.log(1 - 0.7 / 12.4)), 1e-12);
  // Random censored sample: regression and MLE agree within sampling error.
  const data = weibullSample(400, 2, 1000, 31, 1300);
  const mrr = medianRankRegression(data);
  const mle = fitWeibullMle(data);
  close(mrr.shape, mle.shape, 0.12);
  close(mrr.scale, mle.scale, 0.05);
  const e = medianRankRegression({ times: weibullSample(300, 1, 50, 2).times, dist: 'exponential' });
  close(e.rate, 1 / 50, 0.15);
});

test('Weibull functions, B10 life and bathtub composite hazard', () => {
  const f = weibullFunctions(500, { shape: 2, scale: 1000 });
  close(f.reliability, Math.exp(-0.25), 1e-14);
  close(f.hazard, (2 / 1000) * 0.5, 1e-14);
  close(f.pdf, f.hazard * f.reliability, 1e-14);
  // B10 = η·(−ln 0.9)^(1/β) = 1000·0.105361^0.5 = 324.59
  close(bLife({ shape: 2, scale: 1000 }), 324.59, 1e-4);
  close(weibullFunctions(bLife({ shape: 2, scale: 1000 }), { shape: 2, scale: 1000 }).cdf, 0.1, 1e-12);
  // Bathtub: decreasing early, flat mid-life, increasing late.
  const modes = [{ shape: 0.5, scale: 5000 }, { shape: 1, scale: 2000 }, { shape: 4, scale: 3000 }];
  const h = [10, 500, 1200, 4000].map((t) => bathtubHazard(t, modes).hazard);
  assert.ok(h[0] > h[1] && h[3] > h[2] && h[3] > h[1]);
  close(bathtubHazard(1000, modes).reliability, Math.exp(-(0.2 ** 0.5 + 0.5 + (1 / 3) ** 4)), 1e-12);
});

test('renewal function: M(t) = λt for exponential lives; Weibull asymptote t/μ + (CoV² − 1)/2', () => {
  const r = renewalFunction({ cdf: (t) => 1 - Math.exp(-t / 100), tMax: 500, steps: 500 });
  close(r.M[500], 5, 1e-4);
  close(r.M[100], 1, 1e-4);
  const shape = 3;
  const scale = 100;
  const mean = scale * gamma(1 + 1 / shape);
  const cv2 = gamma(1 + 2 / shape) / gamma(1 + 1 / shape) ** 2 - 1;
  close(weibullRenewals({ shape, scale, t: 1000, steps: 2000 }), 1000 / mean + (cv2 - 1) / 2, 2e-3);
  // Short times: M(t) ≈ F(t).
  close(weibullRenewals({ shape, scale, t: 20 }), 1 - Math.exp(-(0.2 ** 3)), 1e-3);
});

test('Crow–AMSAA fit recovers the power-law NHPP parameters; Laplace test detects the trend', () => {
  // Simulate N(t) with E[N] = λ·t^β by transforming a unit-rate Poisson process: t_i = (Γ_i/λ)^(1/β).
  const sim = (beta, lambda, T, seed) => {
    const rng = makeRng(seed);
    const times = [];
    let g = 0;
    for (;;) {
      g -= Math.log(1 - rng.next());
      const t = (g / lambda) ** (1 / beta);
      if (t > T) break;
      times.push(t);
    }
    return times;
  };
  const times = sim(1.6, 0.002, 5000, 4);
  const fit = crowAmsaa({ times, endTime: 5000 });
  close(fit.beta, 1.6, 0.08);
  close(fit.expected(5000), times.length, 1e-10);
  close(fit.lambda * 5000 ** fit.beta, times.length, 1e-10);
  close(fit.intensity, (fit.beta * times.length) / 5000, 1e-10);
  const lap = laplaceTest({ times, endTime: 5000 });
  assert.ok(lap.U > 3);
  assert.equal(lap.trend, 'deteriorating');
  // Reliability growth (β < 1).
  const growth = sim(0.6, 2, 5000, 8);
  assert.ok(crowAmsaa({ times: growth, endTime: 5000 }).beta < 0.75);
  assert.equal(laplaceTest({ times: growth, endTime: 5000 }).trend, 'improving');
  // Hand calculation, failure-terminated: times 100, 300, 600 → β = 3/(ln 6 + ln 2) = 1.2073, λ = 3/600^β
  const h = crowAmsaa({ times: [100, 300, 600] });
  close(h.beta, 3 / (Math.log(6) + Math.log(2)), 1e-12);
  close(h.lambda, 3 / 600 ** h.beta, 1e-12);
  // Laplace, hand: times 20, 50, 90 observed to T = 100: U = (160/3 − 50)/(100·√(1/36)) = 0.2
  close(laplaceTest({ times: [20, 50, 90], endTime: 100 }).U, 0.2, 1e-12);
});

test('MTBF, MTTR, availability and Poisson spares', () => {
  const a = availabilityMetrics({ upTimes: [900, 1100, 1000], downTimes: [8, 12, 10], logisticDelay: 40 });
  close(a.mtbf, 1000, 1e-12);
  close(a.mttr, 10, 1e-12);
  close(a.availability, 1000 / 1010, 1e-12);
  close(a.operationalAvailability, 1000 / 1050, 1e-12);
  // Demand mean 2: Poisson CDF 0.9473 at 4, 0.9834 at 5 → 5 spares for 95 %.
  const s = sparesForServiceLevel({ failureRate: 1e-4, units: 10, period: 2000, serviceLevel: 0.95 });
  close(s.expectedDemand, 2, 1e-12);
  assert.equal(s.spares, 5);
  close(s.achievedLevel, Math.exp(-2) * (1 + 2 + 2 + 4 / 3 + 2 / 3 + 4 / 15), 1e-12);
  assert.equal(sparesForServiceLevel({ failureRate: 1e-4, units: 10, period: 2000, serviceLevel: 0.94 }).spares, 4);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => fitWeibullMle({ times: [1, 2], censored: [true, true] }), /failure/);
  assert.throws(() => fitWeibullMle({ times: [1, -2] }), /times\[1\]/);
  assert.throws(() => crowAmsaa({ times: [5, 3] }), /ascending/);
  assert.throws(() => bLife({ shape: 2, scale: 1 }, 1.5), /p must/);
});
