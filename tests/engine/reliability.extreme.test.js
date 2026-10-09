import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  lMoments,
  fitGumbel,
  fitWeibull,
  fitGev,
  fitGpd,
  fitLognormal,
  fitNormal,
  annualMaxima,
  peaksOverThreshold,
  returnValue,
  bootstrapReturnValue,
  iformContour,
  scatterDiagram,
} from '../../src/engine/reliability/extreme.js';
import { dist, normInv } from '../../src/engine/math/stats.js';
import { makeRng } from '../../src/engine/math/rng.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);

function draw(spec, n, seed) {
  const d = dist(spec);
  const rng = makeRng(seed);
  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) x[i] = d.sample(rng);
  return x;
}

test('lMoments: hand calculation for [1, 2, 3, 4] and population values of an exponential sample', () => {
  // b0 = 2.5 ; b1 = (2·1/3 + 3·2/3 + 4·1)/4 = 5/3 ; l2 = 2 b1 - b0 = 5/6 ; symmetric -> l3 = 0
  const lm = lMoments([4, 2, 1, 3]);
  close(lm.l1, 2.5, 1e-14);
  close(lm.l2, 5 / 6, 1e-14);
  close(lm.l3, 0, 1e-14);
  // exponential(mean 2): l1 = 2, l2 = 1, t3 = 1/3, t4 = 1/6 (Hosking 1990, Table 1)
  const e = lMoments(draw({ dist: 'exponential', rate: 0.5 }, 100000, 3));
  close(e.l1, 2, 0.02);
  close(e.l2, 1, 0.012);
  close(e.t3, 1 / 3, 0.006);
  close(e.t4, 1 / 6, 0.006);
});

test('fitGumbel: parameter recovery (loc 12, scale 2.5) by MLE, moments and L-moments from 20000 seeded samples', () => {
  const x = draw({ dist: 'gumbel', loc: 12, scale: 2.5 }, 20000, 101);
  const mle = fitGumbel(x);
  const mom = fitGumbel(x, { method: 'mom' });
  const lmo = fitGumbel(x, { method: 'lmom' });
  for (const f of [mle, mom, lmo]) {
    assert.equal(f.dist, 'gumbel');
    close(f.loc, 12, 0.08, f.method);
    close(f.scale, 2.5, 0.06, f.method);
  }
  // MLE: asymptotic standard errors, parameters within 4 se, and the highest likelihood of the three
  close(mle.se.loc, 2.5 * Math.sqrt(1.10867 / 20000), 5e-4);
  close(mle.se.scale, 2.5 * Math.sqrt(0.60793 / 20000), 5e-4);
  assert.ok(Math.abs(mle.loc - 12) < 4 * mle.se.loc && Math.abs(mle.scale - 2.5) < 4 * mle.se.scale);
  assert.ok(mle.logLik >= mom.logLik && mle.logLik >= lmo.logLik);
  // the MLE satisfies its likelihood equations: mean(exp(-z)) = 1 and mean(z) - mean(z exp(-z)) = 1
  let s1 = 0;
  let s2 = 0;
  for (const v of x) {
    const z = (v - mle.loc) / mle.scale;
    s1 += Math.exp(-z);
    s2 += z - z * Math.exp(-z);
  }
  close(s1 / x.length, 1, 1e-9);
  close(s2 / x.length, 1, 1e-8);
  assert.throws(() => fitGumbel([1, 1, 1, 1]), RangeError);
  assert.throws(() => fitGumbel(x, { method: 'abc' }), RangeError);
});

test('fitWeibull: two-parameter MLE recovery (scale 3, shape 1.8) and score equation', () => {
  const x = draw({ dist: 'weibull', scale: 3, shape: 1.8 }, 20000, 7);
  const f = fitWeibull(x);
  close(f.scale, 3, 0.04);
  close(f.shape, 1.8, 0.03);
  assert.equal(f.loc, 0);
  // logLik reported equals the direct sum of log densities
  const d = dist(f);
  let ll = 0;
  for (const v of x) ll += Math.log(d.pdf(v));
  close(f.logLik, ll, 1e-6 * Math.abs(ll));
  assert.throws(() => fitWeibull([-1, 2, 3, 4, 5]), RangeError);
});

test('fitWeibull: three-parameter MLE recovery (scale 2.2, shape 1.5, loc 0.6), typical of an Hs marginal', () => {
  const truth = { dist: 'weibull', scale: 2.2, shape: 1.5, loc: 0.6 };
  const x = draw(truth, 20000, 19);
  const f = fitWeibull(x, { threeParam: true });
  close(f.loc, 0.6, 0.03);
  close(f.scale, 2.2, 0.06);
  close(f.shape, 1.5, 0.04);
  // better likelihood than the two-parameter fit and than the true parameters
  assert.ok(f.logLik > fitWeibull(x).logLik);
  const dt = dist(truth);
  let llTrue = 0;
  for (const v of x) llTrue += Math.log(dt.pdf(v));
  assert.ok(f.logLik >= llTrue - 1e-6);
  // a larger shape (regular case, shape > 2)
  const y = draw({ dist: 'weibull', scale: 5, shape: 3, loc: 10 }, 20000, 23);
  const g = fitWeibull(y, { threeParam: true });
  close(g.loc, 10, 0.2);
  close(g.scale, 5, 0.2);
  close(g.shape, 3, 0.15);
});

test('fitGev: L-moment recovery for heavy (xi = 0.15), bounded (xi = -0.2) and Gumbel (xi = 0) tails', () => {
  for (const xi of [0.15, -0.2, 0]) {
    const x = draw({ dist: 'gev', loc: 20, scale: 3, shape: xi }, 40000, 55);
    const f = fitGev(x);
    assert.equal(f.dist, 'gev');
    close(f.shape, xi, 0.02, `xi=${xi} shape`);
    close(f.loc, 20, 0.08, `xi=${xi} loc`);
    close(f.scale, 3, 0.08, `xi=${xi} scale`);
    // the fitted distribution reproduces the sample L-moments exactly (definition of the estimator)
    const d = dist(f);
    const lm = lMoments(x);
    // l1 = mean
    close(d.mean, lm.l1, 1e-8 * Math.abs(lm.l1) + 1e-9, `xi=${xi} mean`);
    // l2 = ∫ x (2F - 1) dF, by the mid-point rule on the quantile function
    let l2 = 0;
    const n = 400000;
    for (let i = 0; i < n; i++) {
      const p = (i + 0.5) / n;
      l2 += d.icdf(p) * (2 * p - 1);
    }
    close(l2 / n, lm.l2, 2e-3 * lm.l2, `xi=${xi} l2`);
  }
});

test('fitGpd: recovery of scale 1.5 and shape (0.1 and -0.2) by L-moments and maximum likelihood', () => {
  for (const xi of [0.1, -0.2]) {
    const x = draw({ dist: 'gpd', loc: 0, scale: 1.5, shape: xi }, 30000, 77);
    const lm = fitGpd(x);
    const ml = fitGpd(x, { method: 'mle' });
    for (const f of [lm, ml]) {
      close(f.scale, 1.5, 0.05, `${f.method} xi=${xi}`);
      close(f.shape, xi, 0.025, `${f.method} xi=${xi}`);
    }
    assert.ok(ml.logLik >= lm.logLik - 1e-6);
  }
  assert.throws(() => fitGpd([-1, 1, 2, 3]), RangeError);
});

test('fitLognormal / fitNormal: parameter recovery and consistency with dist()', () => {
  const x = draw({ dist: 'lognormal', mean: 9, std: 2.5 }, 50000, 5);
  const f = fitLognormal(x);
  const sig = Math.sqrt(Math.log(1 + (2.5 / 9) ** 2));
  close(f.sigma, sig, 0.004);
  close(f.mu, Math.log(9) - 0.5 * sig * sig, 0.005);
  close(f.mean, 9, 0.05);
  close(dist(f).icdf(0.5), Math.exp(f.mu), 1e-9);
  close(dist({ dist: 'lognormal', mu: f.mu, sigma: f.sigma }).mean, f.mean, 1e-9);
  const g = fitNormal(draw({ dist: 'normal', mean: -3, std: 0.7 }, 50000, 6));
  close(g.mean, -3, 0.01);
  close(g.std, 0.7, 0.01);
  assert.equal(g.mu, g.mean);
});

test('returnValue: analytic quantiles — Gumbel x_T = loc - scale ln(-ln(1 - 1/T)), GPD with event rate, Weibull sea states', () => {
  const gum = { dist: 'gumbel', loc: 10, scale: 2 };
  for (const T of [10, 100, 10000]) close(returnValue(gum, T), 10 - 2 * Math.log(-Math.log(1 - 1 / T)), 1e-10);
  // POT: x_T = u + (sigma/xi) ((lambda T)^xi - 1)
  const gpd = { dist: 'gpd', loc: 6, scale: 1.2, shape: -0.1 };
  close(returnValue(gpd, 100, { ratePerYear: 4.5 }), 6 + (1.2 / -0.1) * (Math.pow(450, -0.1) - 1), 1e-10);
  // exponential excess limit
  close(returnValue({ dist: 'gpd', loc: 6, scale: 1.2, shape: 0 }, 50, { ratePerYear: 3 }), 6 + 1.2 * Math.log(150), 1e-10);
  // 3-hourly sea states (2922 per year), Weibull: x_T = loc + scale (ln(N))^(1/shape)
  const N = 2922 * 100;
  close(returnValue({ dist: 'weibull', scale: 2, shape: 1.4, loc: 0.5 }, 100, { ratePerYear: 2922 }), 0.5 + 2 * Math.pow(Math.log(N), 1 / 1.4), 1e-9);
  // GEV 100-year value from a fit to synthetic annual maxima equals the true quantile within sampling error
  const truth = { dist: 'gev', loc: 8, scale: 1.2, shape: -0.1 };
  const fit = fitGev(draw(truth, 5000, 9));
  close(returnValue(fit, 100), dist(truth).icdf(0.99), 0.15);
  assert.throws(() => returnValue(gum, 1), RangeError);
  assert.throws(() => returnValue({ loc: 1 }, 10), RangeError);
});

test('annualMaxima: grouping by UTC year from ISO strings and epoch milliseconds', () => {
  const times = ['2019-03-01T00:00:00Z', '2019-11-20T06:00:00Z', '2020-01-05T00:00:00Z', '2020-12-31T23:00:00Z', '2021-06-15T12:00:00Z', '2021-06-16T12:00:00Z'];
  const values = [3.1, 5.4, 2.2, 6.8, NaN, 4.0];
  const r = annualMaxima(times, values);
  assert.deepEqual(r.years, [2019, 2020, 2021]);
  assert.deepEqual(r.maxima, [5.4, 6.8, 4.0]);
  assert.deepEqual(r.counts, [2, 2, 1]);
  assert.equal(r.timesOfMax[1], Date.parse('2020-12-31T23:00:00Z'));
  const ms = times.map((t) => Date.parse(t));
  assert.deepEqual(annualMaxima(ms, values).maxima, r.maxima);
  assert.throws(() => annualMaxima(['not a date'], [1]), RangeError);
});

test('peaksOverThreshold: declustering keeps one peak per storm (hand-built hourly series)', () => {
  const h = 3600e3;
  const t0 = Date.UTC(2020, 0, 1);
  const n = 24 * 366; // one leap year of hourly data
  const times = Array.from({ length: n }, (_, i) => t0 + i * h);
  const values = new Array(n).fill(1);
  // storm 1: hours 100..110 peaking at 7.5 ; a second burst 20 h later (same storm for a 48 h window) peaking at 6
  for (let i = 100; i <= 110; i++) values[i] = 5 + 2.5 * (1 - Math.abs(i - 105) / 5);
  values[130] = 6;
  // storm 2: hour 400, peak 9 ; storm 3: hours 2000-2002, peak 5.5
  values[400] = 9;
  values[2000] = 5.2;
  values[2001] = 5.5;
  values[2002] = 5.1;
  const r = peaksOverThreshold(times, values, { threshold: 5, minSeparationHours: 48 });
  assert.deepEqual(r.peaks, [7.5, 9, 5.5]);
  assert.deepEqual(r.times, [t0 + 105 * h, t0 + 400 * h, t0 + 2001 * h]);
  close(r.excesses[0], 2.5, 1e-12);
  // record length (n-1) hours -> rate = 3 / years
  close(r.ratePerYear, 3 / (((n - 1) * h) / (365.2425 * 24 * h)), 1e-9);
  // a shorter separation window splits the first storm in two
  assert.deepEqual(peaksOverThreshold(times, values, { threshold: 5, minSeparationHours: 12 }).peaks, [7.5, 6, 9, 5.5]);
  // unsorted input gives the same result
  const idx = makeRng(1).shuffle(Array.from({ length: n }, (_, i) => i));
  const r2 = peaksOverThreshold(idx.map((i) => times[i]), idx.map((i) => values[i]), { threshold: 5 });
  assert.deepEqual(r2.peaks, r.peaks);
  assert.throws(() => peaksOverThreshold(times, values, {}), RangeError);
});

test('peaksOverThreshold + fitGpd + returnValue: 100-year value of a synthetic storm record matches the generating model', () => {
  // storms arrive as a Poisson process (12 per year) with GPD peak excesses over u = 4 (scale 1.1, xi = -0.08)
  const rng = makeRng(404);
  const years = 600; // long synthetic record so that the sampling error of the 100-yr value is ~0.2
  const yearMs = 365.2425 * 24 * 3600e3;
  const gp = dist({ dist: 'gpd', loc: 4, scale: 1.1, shape: -0.08 });
  const times = [0];
  const values = [0];
  let t = 0;
  for (;;) {
    t += rng.exponential(12) * yearMs + 5 * 24 * 3600e3; // at least 5 days apart
    if (t > years * yearMs) break;
    times.push(t);
    values.push(gp.icdf(rng.next()));
  }
  times.push(years * yearMs);
  values.push(0);
  const pot = peaksOverThreshold(times, values, { threshold: 4, minSeparationHours: 48 });
  assert.equal(pot.peaks.length, times.length - 2);
  const fit = fitGpd(pot.excesses, { method: 'mle' });
  const x100 = returnValue({ ...fit, loc: 4 }, 100, { ratePerYear: pot.ratePerYear });
  const truth = gp.icdf(1 - 1 / (pot.ratePerYear * 100));
  close(x100, truth, 0.6, `100-yr ${x100} vs ${truth}`);
});

test('bootstrapReturnValue: percentile band brackets the estimate and has the asymptotic width', () => {
  const truth = { dist: 'gumbel', loc: 10, scale: 2 };
  const x = draw(truth, 200, 33);
  const b = bootstrapReturnValue(x, (s) => fitGumbel(s, { method: 'mom' }), 100, { n: 400, seed: 2 });
  assert.ok(b.lo < b.median && b.median < b.hi);
  assert.ok(b.lo < b.estimate && b.estimate < b.hi);
  const exact = dist(truth).icdf(0.99); // 19.2
  assert.ok(b.lo < exact && exact < b.hi, `band [${b.lo}, ${b.hi}] vs exact ${exact}`);
  // delta-method standard error of the MLE quantile: scale sqrt((1.10867 + 0.51404 y + 0.60793 y²)/n), y = 4.6
  // = 2 sqrt(16.34/200) = 0.57 (MLE); the moment estimator is less efficient, so the band is somewhat wider
  const width = b.hi - b.lo;
  assert.ok(width > 2 * 1.96 * 0.5 && width < 2 * 1.96 * 1.3, `width ${width}`);
  assert.equal(b.nValid, 400);
  const again = bootstrapReturnValue(x, (s) => fitGumbel(s, { method: 'mom' }), 100, { n: 400, seed: 2 });
  assert.equal(again.lo, b.lo);
});

test('iformContour: recovers the generating Hs–Tp model and the analytic contour extremes', () => {
  // generating model (Haver-type): Hs ~ Weibull(scale 2.0, shape 1.45, loc 0.5);
  // ln Tp | Hs ~ N(mu(h), sigma(h)²), mu = 1.3 + 0.55 h^0.4, sigma = 0.06 + 0.25 exp(-0.35 h)
  const n = 60000;
  const rng = makeRng(2718);
  const hsD = dist({ dist: 'weibull', scale: 2.0, shape: 1.45, loc: 0.5 });
  const mu = (h) => 1.3 + 0.55 * Math.pow(h, 0.4);
  const sg = (h) => 0.06 + 0.25 * Math.exp(-0.35 * h);
  const hs = new Float64Array(n);
  const tp = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    hs[i] = hsD.icdf(rng.next());
    tp[i] = Math.exp(mu(hs[i]) + sg(hs[i]) * rng.normal());
  }
  const T = 50;
  const r = iformContour({ hs, tp }, T, { seaStateHours: 3, nPoints: 72 });
  assert.equal(r.hs.length, 72);
  assert.equal(r.tp.length, 72);
  // marginal
  close(r.fit.weibull.loc, 0.5, 0.03);
  close(r.fit.weibull.scale, 2.0, 0.04);
  close(r.fit.weibull.shape, 1.45, 0.03);
  // conditional model: compare the fitted functions with the truth over the populated range
  const muFit = (h) => r.fit.mu.a0 + r.fit.mu.a1 * Math.pow(h, r.fit.mu.a2);
  const sgFit = (h) => r.fit.sigma.b0 + r.fit.sigma.b1 * Math.exp(r.fit.sigma.b2 * h);
  for (const h of [1, 2, 3, 5, 7]) {
    close(muFit(h), mu(h), 0.02, `mu(${h})`);
    close(sgFit(h), sg(h), 0.015, `sigma(${h})`);
  }
  // reliability index: N = 50 yr of 3-hour sea states
  const N = (T * 365.25 * 24) / 3;
  close(r.fit.beta, normInv(1 - 1 / N), 1e-12);
  // the largest Hs on the contour is the marginal N-sea-state quantile (theta = 0)
  const hsExact = hsD.icdf(1 - 1 / N);
  close(r.hsMax, hsExact, 0.03 * hsExact, 'hsMax');
  assert.equal(r.hs[0], r.hsMax);
  close(r.tpAtHsMax, Math.exp(mu(hsExact)), 0.05 * Math.exp(mu(hsExact)), 'tp at hsMax');
  // every contour point lies on the beta-circle of the true model (within fitting error)
  for (let k = 0; k < 72; k += 6) {
    const u1 = hsD.toZ(r.hs[k]);
    const u2 = (Math.log(r.tp[k]) - mu(r.hs[k])) / sg(r.hs[k]);
    close(Math.hypot(u1, u2), r.fit.beta, 0.35, `point ${k}`);
  }
  // a longer return period gives a contour that encloses the shorter one in Hs
  const r2 = iformContour({ hs, tp }, 1, { seaStateHours: 3 });
  assert.ok(r2.hsMax < r.hsMax);
  assert.throws(() => iformContour({ hs: [1, 2], tp: [1, 2] }, 50), RangeError);
});

test('scatterDiagram: hand-checked table', () => {
  const hs = [0.2, 0.4, 0.6, 1.2, 1.3, 2.9];
  const tp = [4.5, 4.9, 5.1, 6.0, 6.9, 9.99];
  const r = scatterDiagram(hs, tp, { hsBin: 1, tpBin: 2 });
  assert.deepEqual(r.hsEdges, [0, 1, 2, 3]);
  assert.deepEqual(r.tpEdges, [4, 6, 8, 10]);
  assert.deepEqual(r.counts, [[3, 0, 0], [0, 2, 0], [0, 0, 1]]);
  close(r.prob[0][0], 0.5, 1e-15);
  let s = 0;
  for (const row of r.prob) for (const p of row) s += p;
  close(s, 1, 1e-14);
  const d = scatterDiagram(hs, tp);
  assert.equal(d.hsEdges[1] - d.hsEdges[0], 0.5);
  assert.equal(d.tpEdges[1] - d.tpEdges[0], 1);
  assert.equal(d.n, 6);
});
