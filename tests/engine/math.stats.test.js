import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  mean,
  variance,
  std,
  min,
  max,
  median,
  percentile,
  skewness,
  kurtosis,
  histogram,
  covariance,
  correlation,
  rankCorrelation,
  dist,
  normCdf,
  normInv,
  linearRegression,
  polyfit,
  polyval,
  multipleRegression,
  metrics,
  ksTest,
} from '../../src/engine/math/stats.js';
import { makeRng } from '../../src/engine/math/rng.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} expected ${b}, got ${a} (rel ${tol})`);

test('descriptive statistics: hand calculation on [2, 4, 4, 4, 5, 5, 7, 9]', () => {
  const x = [2, 4, 4, 4, 5, 5, 7, 9];
  // sum = 40, mean = 5; squared deviations 9+1+1+1+0+0+4+16 = 32
  assert.equal(mean(x), 5);
  close(variance(x), 32 / 7, 1e-14);
  close(variance(x, true), 4, 1e-14);
  close(std(x, true), 2, 1e-14);
  assert.equal(min(x), 2);
  assert.equal(max(x), 9);
  assert.equal(median(x), 4.5);
  // type-7 percentile: h = 7 * 0.25 = 1.75 -> 4 ; p = 90: h = 6.3 -> 7 + 0.3 * 2 = 7.6
  close(percentile(x, 25), 4, 1e-14);
  close(percentile(x, 90), 7.6, 1e-12);
  assert.equal(percentile(x, 0), 2);
  assert.equal(percentile(x, 100), 9);
  // third central moment = (-27 -1 -1 -1 + 0 + 0 + 8 + 64)/8 = 5.25 ; m2^1.5 = 8
  close(skewness(x), 5.25 / 8, 1e-14);
  // fourth central moment = (81 + 1 + 1 + 1 + 0 + 0 + 16 + 256)/8 = 44.5 ; /16 = 2.78125
  close(kurtosis(x, false), 2.78125, 1e-14);
  close(kurtosis(x), -0.21875, 1e-14);
  assert.throws(() => mean([]), RangeError);
  assert.throws(() => percentile(x, 101), RangeError);
});

test('histogram: counts, edges and density', () => {
  const h = histogram([0, 0.1, 0.4, 0.5, 0.9, 1], 2);
  assert.deepEqual(h.edges, [0, 0.5, 1]);
  assert.deepEqual(h.counts, [3, 3]);
  close(h.density[0] * 0.5 + h.density[1] * 0.5, 1, 1e-14);
  const h2 = histogram([1, 2, 2, 3, 10], [0, 2, 4]);
  assert.deepEqual(h2.counts, [1, 3]);
});

test('covariance, correlation and rank correlation', () => {
  const a = [1, 2, 3, 4, 5];
  const b = [2, 4, 6, 8, 10];
  close(covariance(a, b), 5, 1e-14);
  close(correlation(a, b), 1, 1e-14);
  close(correlation(a, [5, 4, 3, 2, 1]), -1, 1e-14);
  // monotone but non-linear relation: Spearman = 1, Pearson < 1
  const c = a.map((v) => Math.exp(v));
  close(rankCorrelation(a, c), 1, 1e-14);
  assert.ok(correlation(a, c) < 0.95);
  // textbook Spearman example with d^2 sum: ranks (1,2,3,4,5) vs (2,1,4,3,5): rho = 1 - 6*4/(5*24) = 0.8
  close(rankCorrelation([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]), 0.8, 1e-14);
});

test('standard normal: tabulated values and inverse round trip into the far tails', () => {
  close(normCdf(0), 0.5, 1e-15);
  close(normCdf(1.959963984540054), 0.975, 1e-12);
  close(normCdf(-3), 1.3498980316300946e-3, 1e-14);
  rel(normCdf(-6), 9.865876450376946e-10, 1e-9);
  close(normInv(0.975), 1.959963984540054, 1e-10);
  close(normInv(0.5), 0, 1e-15);
  rel(normInv(1e-10), -6.361340902404056, 1e-9);
  for (const p of [1e-12, 1e-6, 0.01, 0.2, 0.5, 0.7, 0.99, 1 - 1e-9]) rel(normCdf(normInv(p)), p, 1e-9, `p=${p}`);
  assert.equal(normInv(0), -Infinity);
  assert.throws(() => normInv(1.1), RangeError);
});

const SPECS = [
  { dist: 'normal', mean: 3, std: 2 },
  { dist: 'lognormal', mean: 10, std: 3 },
  { dist: 'uniform', min: -1, max: 4 },
  { dist: 'triangular', min: 1, mode: 2, max: 6 },
  { dist: 'pert', min: 10, mode: 14, max: 30 },
  { dist: 'beta', a: 2.5, b: 4, min: 1, max: 3 },
  { dist: 'weibull', scale: 2.2, shape: 1.6, loc: 0.4 },
  { dist: 'gumbel', loc: 5, scale: 1.5 },
  { dist: 'gev', loc: 5, scale: 1.5, shape: 0.15 },
  { dist: 'gev', loc: 5, scale: 1.5, shape: -0.2 },
  { dist: 'gpd', loc: 1, scale: 2, shape: 0.1 },
  { dist: 'gpd', loc: 1, scale: 2, shape: -0.25 },
  { dist: 'exponential', rate: 0.4 },
  { dist: 'gamma', shape: 3.5, scale: 1.2 },
  { dist: 'gamma', shape: 0.7, scale: 2 },
  { dist: 'rayleigh', sigma: 1.7 },
];
const label = (s) => `${s.dist}${s.shape !== undefined && (s.dist === 'gev' || s.dist === 'gpd') ? `(xi=${s.shape})` : ''}`;

/** composite Simpson on [a, b] with n (even) intervals — independent of the engine quadrature */
function simpson(f, a, b, n) {
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * f(a + i * h);
  return (s * h) / 3;
}

test('dist: cdf(icdf(p)) = p for every distribution', () => {
  for (const spec of SPECS) {
    const d = dist(spec);
    for (const p of [1e-6, 0.001, 0.05, 0.3, 0.5, 0.8, 0.99, 0.999999]) {
      const x = d.icdf(p);
      assert.ok(Number.isFinite(x), `${label(spec)} icdf(${p})`);
      close(d.cdf(x), p, 1e-10 + 1e-9 * p, `${label(spec)} p=${p}`);
    }
    assert.throws(() => d.icdf(-0.1), RangeError);
    assert.throws(() => d.icdf(1.5), RangeError);
  }
});

test('dist: pdf integrates to 1 and to the cdf, pdf = d(cdf)/dx', () => {
  for (const spec of SPECS) {
    const d = dist(spec);
    const lo = d.icdf(1e-7);
    const hi = d.icdf(1 - 1e-7);
    // integrate between quantiles that avoid endpoint singularities (gamma with shape < 1)
    const a = d.icdf(0.02);
    const b = d.icdf(0.98);
    close(simpson(d.pdf, a, b, 4000), 0.96, 2e-6, `${label(spec)} integral of pdf between the 2 % and 98 % quantiles`);
    if (!(spec.dist === 'gamma' && spec.shape < 1)) {
      close(simpson(d.pdf, lo, hi, 40000), 1, 5e-5, `${label(spec)} total probability`);
    }
    for (const p of [0.1, 0.5, 0.9]) {
      const x = d.icdf(p);
      const h = 1e-5 * (b - a);
      close((d.cdf(x + h) - d.cdf(x - h)) / (2 * h), d.pdf(x), 1e-6 * Math.max(1, d.pdf(x)), `${label(spec)} pdf vs dcdf at p=${p}`);
    }
  }
});

test('dist: analytic mean/std agree with numerical integration of x·pdf', () => {
  for (const spec of SPECS) {
    const d = dist(spec);
    // E[X] = ∫ icdf(p) dp, E[X²] = ∫ icdf(p)² dp (substitution p = F(x)), mid-point rule on a graded grid
    let m1 = 0;
    let m2 = 0;
    // (fewer nodes for the families whose icdf is a Newton iteration on the incomplete beta/gamma function)
    const n = ['beta', 'pert', 'gamma'].includes(spec.dist) ? 20000 : 200000;
    for (let i = 0; i < n; i++) {
      const x = d.icdf((i + 0.5) / n);
      m1 += x;
      m2 += x * x;
    }
    m1 /= n;
    const sd = Math.sqrt(m2 / n - m1 * m1);
    const heavy = spec.dist === 'gev' || spec.dist === 'gpd' || spec.dist === 'lognormal';
    rel(m1, d.mean, heavy ? 2e-3 : 2e-4, `${label(spec)} mean`);
    rel(sd, d.std, heavy ? 3e-2 : 2e-3, `${label(spec)} std`);
  }
});

test('dist: sample moments from 2e5 seeded draws match the analytic mean and std', () => {
  for (const spec of SPECS) {
    const d = dist(spec);
    const rng = makeRng(20240601);
    const n = 200000;
    let s = 0;
    let s2 = 0;
    const [lo, hi] = d.support;
    for (let i = 0; i < n; i++) {
      const x = d.sample(rng);
      assert.ok(x >= lo && x <= hi);
      s += x;
      s2 += x * x;
    }
    const m = s / n;
    const sd = Math.sqrt(s2 / n - m * m);
    // mean within 5 standard errors; std within 2 % (4 % for the heavy-tailed GEV/GPD with xi > 0)
    close(m, d.mean, (5 * d.std) / Math.sqrt(n), `${label(spec)} sample mean`);
    const heavy = (spec.dist === 'gev' || spec.dist === 'gpd') && spec.shape > 0;
    rel(sd, d.std, heavy ? 0.04 : 0.02, `${label(spec)} sample std`);
  }
});

test('dist: samples pass a Kolmogorov–Smirnov test against their own cdf (gamma and beta use rejection samplers)', () => {
  for (const spec of [SPECS[5], SPECS[13], SPECS[14], SPECS[4], SPECS[0]]) {
    const d = dist(spec);
    const rng = makeRng(77);
    const x = Array.from({ length: 20000 }, () => d.sample(rng));
    const ks = ksTest(x, d.cdf);
    assert.ok(ks.p > 0.001, `${label(spec)} KS p = ${ks.p}, D = ${ks.D}`);
  }
  // and the test rejects a wrong distribution
  const rng = makeRng(1);
  const x = Array.from({ length: 2000 }, () => rng.normal());
  assert.ok(ksTest(x, dist({ dist: 'normal', mean: 0.2, std: 1 }).cdf).p < 1e-6);
  assert.ok(ksTest(x, normCdf).p > 0.01);
});

test('dist: closed-form parameter relations', () => {
  // lognormal parameterised by the mean and std of X: sigma² = ln(1 + V²), mu = ln(m) - sigma²/2
  const ln = dist({ dist: 'lognormal', mean: 10, std: 3 });
  close(ln.params.sigma, Math.sqrt(Math.log(1.09)), 1e-14);
  close(ln.params.mu, Math.log(10) - 0.5 * Math.log(1.09), 1e-14);
  close(ln.mean, 10, 1e-12);
  close(ln.std, 3, 1e-12);
  close(ln.icdf(0.5), 10 / Math.sqrt(1.09), 1e-12); // median = m / sqrt(1 + V²)
  // { mean, cov } convenience
  const n2 = dist({ dist: 'normal', mean: 200, cov: 0.1 });
  close(n2.std, 20, 1e-12);
  close(dist({ dist: 'lognormal', mean: 50, cov: 0.2 }).std, 10, 1e-10);
  // beta-PERT: mean = (min + 4 mode + max)/6, alpha = 1 + 4 (mode-min)/(max-min)
  const pert = dist({ dist: 'pert', min: 10, mode: 14, max: 30 });
  close(pert.mean, (10 + 4 * 14 + 30) / 6, 1e-12);
  close(pert.params.a, 1.8, 1e-12);
  close(pert.params.b, 4.2, 1e-12);
  // PERT pdf peaks at the mode
  assert.ok(pert.pdf(14) > pert.pdf(13.5) && pert.pdf(14) > pert.pdf(14.5));
  // triangular: mean = (a + b + c)/3, cdf(mode) = (c-a)/(b-a)
  const tri = dist({ dist: 'triangular', min: 1, mode: 2, max: 6 });
  close(tri.mean, 3, 1e-14);
  close(tri.cdf(2), 0.2, 1e-14);
  // Gumbel: mean = loc + 0.5772 scale, std = pi scale / sqrt(6)
  const gu = dist({ dist: 'gumbel', loc: 5, scale: 1.5 });
  close(gu.mean, 5 + 0.5772156649015329 * 1.5, 1e-13);
  close(gu.std, (Math.PI * 1.5) / Math.sqrt(6), 1e-13);
  // Weibull with shape 1 is exponential; with shape 2 Rayleigh (scale = sigma sqrt 2)
  const w1 = dist({ dist: 'weibull', scale: 2.5, shape: 1 });
  const ex = dist({ dist: 'exponential', rate: 0.4 });
  close(w1.cdf(3), ex.cdf(3), 1e-14);
  const w2 = dist({ dist: 'weibull', scale: 1.7 * Math.SQRT2, shape: 2 });
  const ra = dist({ dist: 'rayleigh', sigma: 1.7 });
  close(w2.icdf(0.9), ra.icdf(0.9), 1e-12);
  close(w2.mean, ra.mean, 1e-12);
  // gamma with integer shape 2: cdf = 1 - e^-t (1 + t)
  const ga = dist({ dist: 'gamma', shape: 2, scale: 3 });
  close(ga.cdf(4.5), 1 - Math.exp(-1.5) * 2.5, 1e-12);
  // beta(1,1) is uniform
  close(dist({ dist: 'beta', a: 1, b: 1 }).cdf(0.3), 0.3, 1e-12);
  // deterministic
  const de = dist({ dist: 'deterministic', value: 4.2 });
  assert.equal(de.icdf(0.3), 4.2);
  assert.equal(de.std, 0);
  assert.equal(de.sample(makeRng(1)), 4.2);
});

test('dist: GEV sign convention (xi > 0 heavy upper tail, xi < 0 bounded above, xi -> 0 Gumbel) and GPD limits', () => {
  const heavy = dist({ dist: 'gev', loc: 0, scale: 1, shape: 0.2 });
  const light = dist({ dist: 'gev', loc: 0, scale: 1, shape: -0.2 });
  const gum = dist({ dist: 'gumbel', loc: 0, scale: 1 });
  assert.ok(heavy.icdf(0.9999) > gum.icdf(0.9999) && gum.icdf(0.9999) > light.icdf(0.9999));
  assert.equal(light.support[1], 5); // loc - scale/xi
  assert.equal(heavy.support[0], -5);
  assert.equal(light.cdf(6), 1);
  assert.equal(heavy.cdf(-6), 0);
  const near = dist({ dist: 'gev', loc: 0, scale: 1, shape: 1e-7 });
  close(near.icdf(0.99), gum.icdf(0.99), 1e-5);
  close(dist({ dist: 'gev', loc: 0, scale: 1, shape: 0 }).cdf(1.3), gum.cdf(1.3), 1e-15);
  // GEV mean = loc + scale (Γ(1-xi) - 1)/xi : xi = 0.5 -> Γ(0.5) = sqrt(pi)
  close(dist({ dist: 'gev', loc: 1, scale: 2, shape: 0.5 }).mean, 1 + 2 * (Math.sqrt(Math.PI) - 1) / 0.5, 1e-10);
  assert.equal(dist({ dist: 'gev', loc: 1, scale: 2, shape: 0.5 }).std, Infinity);
  // GPD with xi = 0 is exponential; xi = -1 is uniform on [loc, loc + scale]
  close(dist({ dist: 'gpd', loc: 0, scale: 2.5, shape: 0 }).cdf(3), 1 - Math.exp(-1.2), 1e-14);
  close(dist({ dist: 'gpd', loc: 1, scale: 2, shape: -1 }).cdf(2), 0.5, 1e-14);
  close(dist({ dist: 'gpd', loc: 1, scale: 2, shape: 0.25 }).mean, 1 + 2 / 0.75, 1e-14);
});

test('dist: moment-matched construction from { mean, std } reproduces the moments', () => {
  for (const name of ['uniform', 'triangular', 'gumbel', 'gamma', 'weibull', 'beta']) {
    const spec = name === 'beta' ? { dist: name, mean: 0.3, std: 0.1 } : { dist: name, mean: 12, std: 3 };
    const d = dist(spec);
    rel(d.mean, spec.mean, 1e-9, name);
    rel(d.std, spec.std, 1e-9, name);
  }
  rel(dist({ dist: 'weibull', mean: 12, cov: 0.25 }).std, 3, 1e-9);
  close(dist({ dist: 'exponential', mean: 4 }).params.rate, 0.25, 1e-15);
});

test('dist: invalid parameters throw RangeError naming the argument', () => {
  assert.throws(() => dist({ dist: 'normal', mean: 0, std: -1 }), /std/);
  assert.throws(() => dist({ dist: 'lognormal', mean: -1, std: 1 }), /mean/);
  assert.throws(() => dist({ dist: 'uniform', min: 2, max: 1 }), /max/);
  assert.throws(() => dist({ dist: 'triangular', min: 0, mode: 5, max: 1 }), /mode/);
  assert.throws(() => dist({ dist: 'weibull', scale: 1, shape: 0 }), /shape/);
  assert.throws(() => dist({ dist: 'gev', loc: 0, scale: 0, shape: 0 }), /scale/);
  assert.throws(() => dist({ dist: 'nope' }), RangeError);
  assert.throws(() => dist({ dist: 'normal', std: 1 }), /mean/);
});

test('dist: isoprobabilistic transforms fromZ/toZ are consistent and stay finite in the tails', () => {
  for (const spec of SPECS) {
    const d = dist(spec);
    for (const z of [-3, -1, 0, 0.5, 2, 4]) {
      const x = d.fromZ(z);
      close(d.toZ(x), z, 1e-7, `${label(spec)} z=${z}`);
    }
    assert.ok(Number.isFinite(d.fromZ(7.5)), label(spec));
    assert.ok(Number.isFinite(d.fromZ(-7.5)), label(spec));
  }
});

test('linearRegression: exact line and textbook standard errors', () => {
  const e = linearRegression([0, 1, 2, 3], [1, 3, 5, 7]);
  close(e.slope, 2, 1e-14);
  close(e.intercept, 1, 1e-14);
  close(e.r2, 1, 1e-14);
  // x = 1..5, y = 2, 4, 5, 4, 5: Sxx = 10, Sxy = 6, slope 0.6, intercept 2.2, SSE = 2.4, s² = 0.8
  const r = linearRegression([1, 2, 3, 4, 5], [2, 4, 5, 4, 5]);
  close(r.slope, 0.6, 1e-14);
  close(r.intercept, 2.2, 1e-14);
  close(r.r2, 0.6, 1e-14);
  close(r.seSlope, Math.sqrt(0.8 / 10), 1e-14);
  close(r.seIntercept, Math.sqrt(0.8 * (1 / 5 + 9 / 10)), 1e-14);
  assert.throws(() => linearRegression([1, 1, 1], [1, 2, 3]), RangeError);
});

test('polyfit / polyval: recovers a cubic exactly, highest degree first', () => {
  const c = [2, -3, 0.5, 7];
  const x = [-2, -1, 0, 0.5, 1, 2, 3, 4.5];
  const y = x.map((v) => polyval(c, v));
  const fit = polyfit(x, y, 3);
  for (let i = 0; i < 4; i++) close(fit[i], c[i], 1e-10);
  close(polyval([1, 0, -1], 3), 8, 1e-14);
  assert.deepEqual(polyval([1, 1], [1, 2]), [2, 3]);
  // badly offset abscissae (years) remain well conditioned
  const yrs = [2001, 2002, 2003, 2004, 2005, 2006];
  const f2 = polyfit(yrs, yrs.map((v) => 3 * (v - 2000) ** 2 + 1), 2);
  close(polyval(f2, 2003.5), 3 * 3.5 ** 2 + 1, 1e-6);
});

test('multipleRegression: exact plane and noisy recovery with standard errors', () => {
  const X = [[1, 2], [2, 1], [3, 5], [4, 3], [5, 8], [6, 1]];
  const y = X.map(([a, b]) => 1.5 + 2 * a - 0.5 * b);
  const r = multipleRegression(X, y);
  close(r.beta[0], 1.5, 1e-10);
  close(r.beta[1], 2, 1e-10);
  close(r.beta[2], -0.5, 1e-10);
  close(r.r2, 1, 1e-12);
  const rng = makeRng(3);
  const Xn = [];
  const yn = [];
  for (let i = 0; i < 2000; i++) {
    const a = rng.next() * 10;
    const b = rng.next() * 4;
    Xn.push([a, b]);
    yn.push(4 + 0.7 * a - 1.2 * b + 0.3 * rng.normal());
  }
  const rn = multipleRegression(Xn, yn);
  const truth = [4, 0.7, -1.2];
  for (let j = 0; j < 3; j++) {
    assert.ok(rn.se[j] > 0);
    assert.ok(Math.abs(rn.beta[j] - truth[j]) < 4 * rn.se[j]);
  }
  // se of a slope ≈ sigma / (sd_x sqrt(n)) = 0.3 / (10/sqrt(12) * sqrt(2000))
  close(rn.se[1], 0.3 / ((10 / Math.sqrt(12)) * Math.sqrt(2000)), 3e-4);
  assert.equal(rn.residuals.length, 2000);
});

test('metrics: hand calculation', () => {
  const m = metrics([1, 2, 3, 5], [1, 2, 4, 4]);
  // errors 0, 0, -1, 1
  close(m.rmse, Math.sqrt(0.5), 1e-14);
  close(m.mae, 0.5, 1e-14);
  close(m.bias, 0, 1e-14);
  close(m.maxAbs, 1, 1e-14);
  close(m.nrmse, Math.sqrt(0.5) / 3, 1e-14);
  // SST of obs (mean 2.75) = 3.0625 + 0.5625 + 1.5625 + 1.5625 = 6.75
  close(m.r2, 1 - 2 / 6.75, 1e-14);
  assert.equal(m.n, 4);
  close(metrics([2, 3], [1, 2]).bias, 1, 1e-14);
});
