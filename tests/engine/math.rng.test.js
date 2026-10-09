import { test } from 'node:test';
import assert from 'node:assert/strict';
import { makeRng } from '../../src/engine/math/rng.js';

test('rng: same seed gives the same stream, different seeds differ', () => {
  const a = makeRng(123);
  const b = makeRng(123);
  const c = makeRng(124);
  let same = 0;
  for (let i = 0; i < 1000; i++) {
    const va = a.next();
    assert.equal(va, b.next());
    if (va === c.next()) same++;
    assert.ok(va >= 0 && va < 1);
  }
  assert.ok(same < 3);
  // fractional and large seeds are valid and distinct
  assert.notEqual(makeRng(0.5).next(), makeRng(0.25).next());
  assert.notEqual(makeRng(2 ** 40 + 1).next(), makeRng(1).next());
  assert.throws(() => makeRng(NaN), RangeError);
});

test('rng: uniform moments (mean 1/2, variance 1/12) and chi-square uniformity over 100 bins, n = 2e5', () => {
  const r = makeRng(7);
  const n = 200000;
  const bins = new Array(100).fill(0);
  let s = 0;
  let s2 = 0;
  for (let i = 0; i < n; i++) {
    const u = r.next();
    s += u;
    s2 += u * u;
    bins[Math.floor(u * 100)]++;
  }
  const m = s / n;
  // standard error of the mean = sqrt(1/12/n) = 6.5e-4 -> 5 sigma
  assert.ok(Math.abs(m - 0.5) < 3.3e-3);
  assert.ok(Math.abs(s2 / n - m * m - 1 / 12) < 1e-3);
  let chi2 = 0;
  for (const c of bins) chi2 += (c - n / 100) ** 2 / (n / 100);
  // chi-square with 99 dof: mean 99, sd 14.1; accept within ~4.5 sd
  assert.ok(chi2 > 40 && chi2 < 165, `chi2 = ${chi2}`);
});

test('rng: serial correlation of successive uniforms is negligible (|r1| < 5/sqrt(n))', () => {
  const r = makeRng(99);
  const n = 200000;
  let prev = r.next();
  let sxy = 0;
  for (let i = 0; i < n; i++) {
    const u = r.next();
    sxy += (prev - 0.5) * (u - 0.5);
    prev = u;
  }
  assert.ok(Math.abs((sxy / n) * 12) < 5 / Math.sqrt(n));
});

test('rng: normal deviates have mean 0, variance 1, skewness 0, kurtosis 3 and P(|z| > 1.96) = 0.05', () => {
  const r = makeRng(2024);
  const n = 200000;
  let s1 = 0;
  let s2 = 0;
  let s3 = 0;
  let s4 = 0;
  let tail = 0;
  for (let i = 0; i < n; i++) {
    const z = r.normal();
    s1 += z;
    s2 += z * z;
    s3 += z * z * z;
    s4 += z * z * z * z;
    if (Math.abs(z) > 1.959964) tail++;
  }
  assert.ok(Math.abs(s1 / n) < 0.012);
  assert.ok(Math.abs(s2 / n - 1) < 0.016);
  assert.ok(Math.abs(s3 / n) < 0.03);
  assert.ok(Math.abs(s4 / n - 3) < 0.11);
  assert.ok(Math.abs(tail / n - 0.05) < 0.0025);
});

test('rng: int(n) covers 0..n-1 uniformly and shuffle is an unbiased permutation', () => {
  const r = makeRng(5);
  const counts = new Array(6).fill(0);
  for (let i = 0; i < 60000; i++) counts[r.int(6)]++;
  for (const c of counts) assert.ok(Math.abs(c - 10000) < 500);
  assert.throws(() => r.int(0), RangeError);
  // each of the 6 permutations of [0,1,2] appears with probability 1/6
  const seen = new Map();
  for (let i = 0; i < 60000; i++) {
    const key = r.shuffle([0, 1, 2]).join('');
    seen.set(key, (seen.get(key) || 0) + 1);
  }
  assert.equal(seen.size, 6);
  for (const c of seen.values()) assert.ok(Math.abs(c - 10000) < 500);
  const arr = [1, 2, 3, 4, 5, 6, 7];
  assert.deepEqual(r.shuffle(arr.slice()).sort(), arr);
});
