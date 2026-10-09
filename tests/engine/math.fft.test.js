import test from 'node:test';
import assert from 'node:assert/strict';
import { fft, ifft, psdWelch, spectralMoments, synthesizeFromSpectrum, windowFunction } from '../../src/engine/math/fft.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

/** Linear congruential generator (Numerical Recipes constants), uniform in [0, 1). */
function lcg(seed) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** Direct O(N²) evaluation of the DFT definition. */
function dft(re, im) {
  const n = re.length;
  const R = new Float64Array(n);
  const I = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      const a = (-2 * Math.PI * ((k * j) % n)) / n;
      R[k] += re[j] * Math.cos(a) - im[j] * Math.sin(a);
      I[k] += re[j] * Math.sin(a) + im[j] * Math.cos(a);
    }
  }
  return { re: R, im: I };
}

function trapzArea(f, s) {
  let a = 0;
  for (let k = 1; k < f.length; k++) a += 0.5 * (s[k] + s[k - 1]) * (f[k] - f[k - 1]);
  return a;
}

function variance(x) {
  let m = 0;
  for (let i = 0; i < x.length; i++) m += x[i];
  m /= x.length;
  let v = 0;
  for (let i = 0; i < x.length; i++) v += (x[i] - m) * (x[i] - m);
  return v / x.length;
}

test('fft: hand-computed 4-point transform of [1, 2, 3, 4] = [10, −2+2i, −2, −2−2i]', () => {
  const { re, im } = fft([1, 2, 3, 4]);
  const er = [10, -2, -2, -2];
  const ei = [0, 2, 0, -2];
  for (let k = 0; k < 4; k++) {
    close(re[k], er[k], 1e-14);
    close(im[k], ei[k], 1e-14);
  }
  assert.ok(re instanceof Float64Array && im instanceof Float64Array);
});

test('fft: matches the DFT definition for power-of-two, composite, odd and prime lengths (complex input)', () => {
  const rnd = lcg(42);
  for (const n of [1, 2, 3, 8, 12, 15, 64, 97, 100, 127]) {
    const xr = Array.from({ length: n }, () => rnd() - 0.5);
    const xi = Array.from({ length: n }, () => rnd() - 0.5);
    const a = fft(xr, xi);
    const b = dft(xr, xi);
    for (let k = 0; k < n; k++) {
      close(a.re[k], b.re[k], 1e-11, `n=${n} re[${k}]`);
      close(a.im[k], b.im[k], 1e-11, `n=${n} im[${k}]`);
    }
  }
});

test('fft: impulse -> flat spectrum, constant -> DC only; inputs are not modified', () => {
  for (const n of [16, 30]) {
    const imp = new Array(n).fill(0);
    imp[0] = 1;
    const a = fft(imp);
    for (let k = 0; k < n; k++) {
      close(a.re[k], 1, 1e-13);
      close(a.im[k], 0, 1e-13);
    }
    assert.equal(imp[0], 1);
    const c = fft(new Float64Array(n).fill(3));
    close(c.re[0], 3 * n, 1e-11);
    for (let k = 1; k < n; k++) close(Math.hypot(c.re[k], c.im[k]), 0, 1e-11);
  }
  assert.throws(() => fft([]), RangeError);
  assert.throws(() => fft([1, 2], [1]), RangeError);
});

test('ifft(fft(x)) = x for lengths 1024 (radix-2) and 1000, 1023 (Bluestein)', () => {
  const rnd = lcg(7);
  for (const n of [1024, 1000, 1023]) {
    const xr = Float64Array.from({ length: n }, () => rnd() - 0.5);
    const xi = Float64Array.from({ length: n }, () => rnd() - 0.5);
    const F = fft(xr, xi);
    const b = ifft(F.re, F.im);
    let e = 0;
    for (let k = 0; k < n; k++) e = Math.max(e, Math.abs(b.re[k] - xr[k]), Math.abs(b.im[k] - xi[k]));
    assert.ok(e < 1e-12, `n=${n}: round-trip error ${e}`);
  }
});

test('fft: Parseval — Σ|x|² = (1/N) Σ|X|² for N = 4096 and N = 3000', () => {
  const rnd = lcg(2718);
  for (const n of [4096, 3000]) {
    const x = Float64Array.from({ length: n }, () => rnd() - 0.5);
    const { re, im } = fft(x);
    let et = 0;
    let ef = 0;
    for (let k = 0; k < n; k++) {
      et += x[k] * x[k];
      ef += re[k] * re[k] + im[k] * im[k];
    }
    close(ef / n / et, 1, 1e-12, `n=${n}`);
    // real input: Hermitian symmetry X[N−k] = conj(X[k])
    close(re[n - 5], re[5], 1e-10);
    close(im[n - 5], -im[5], 1e-10);
  }
});

test('fft: single-tone detection — A cos(2π k0 n/N + φ) gives |X[k0]| = A N/2 and phase φ, other bins empty', () => {
  for (const [n, k0] of [[512, 37], [600, 41]]) {
    const A = 2.5;
    const phi = 0.7;
    const x = Float64Array.from({ length: n }, (_, j) => A * Math.cos((2 * Math.PI * k0 * j) / n + phi));
    const { re, im } = fft(x);
    let kmax = 0;
    let vmax = 0;
    for (let k = 0; k <= n / 2; k++) {
      const m = Math.hypot(re[k], im[k]);
      if (m > vmax) {
        vmax = m;
        kmax = k;
      }
      if (k !== k0) assert.ok(m < 1e-9 * n, `leakage in bin ${k}: ${m}`);
    }
    assert.equal(kmax, k0);
    close(vmax, (A * n) / 2, 1e-9 * n);
    close(Math.atan2(im[k0], re[k0]), phi, 1e-11);
  }
});

test('fft: linearity and the shift theorem x[n−m] <-> X[k] e^{−2πikm/N} (N = 90)', () => {
  const rnd = lcg(5);
  const n = 90;
  const m = 7;
  const x = Array.from({ length: n }, () => rnd());
  const y = Array.from({ length: n }, () => rnd());
  const X = fft(x);
  const Y = fft(y);
  const Z = fft(x.map((v, i) => 2 * v - 3 * y[i]));
  const S = fft(x.map((_, i) => x[(i - m + n) % n]));
  for (let k = 0; k < n; k++) {
    close(Z.re[k], 2 * X.re[k] - 3 * Y.re[k], 1e-11);
    close(Z.im[k], 2 * X.im[k] - 3 * Y.im[k], 1e-11);
    const a = (-2 * Math.PI * k * m) / n;
    close(S.re[k], X.re[k] * Math.cos(a) - X.im[k] * Math.sin(a), 1e-11);
    close(S.im[k], X.re[k] * Math.sin(a) + X.im[k] * Math.cos(a), 1e-11);
  }
});

test('windowFunction: Hann is periodic with mean 1/2 and mean square 3/8; unknown name throws', () => {
  const w = windowFunction('hann', 256);
  assert.equal(w[0], 0);
  close(w[128], 1, 1e-15);
  let m = 0;
  let m2 = 0;
  for (const v of w) {
    m += v;
    m2 += v * v;
  }
  close(m / 256, 0.5, 1e-14);
  close(m2 / 256, 0.375, 1e-14);
  assert.deepEqual(Array.from(windowFunction('rect', 3)), [1, 1, 1]);
  assert.throws(() => windowFunction('kaiser', 8), RangeError);
});

test('psdWelch: white noise — area under the PSD equals the variance, level is flat at σ²/(fs/2)', () => {
  const rnd = lcg(1234);
  const n = 1 << 16;
  const fs = 20;
  // sum of 12 uniforms minus 6: unit-variance, nearly Gaussian
  const x = Float64Array.from({ length: n }, () => {
    let s = -6;
    for (let k = 0; k < 12; k++) s += rnd();
    return 3 * s + 5; // σ = 3, mean 5 (removed by detrending)
  });
  const v = variance(x);
  close(v, 9, 0.2);
  const { f, psd, nSegments } = psdWelch(x, fs, { nfft: 1024 });
  assert.equal(f.length, 513);
  assert.equal(f[0], 0);
  close(f[512], fs / 2, 1e-12);
  close(f[1], fs / 1024, 1e-15);
  assert.equal(nSegments, 127);
  const area = trapzArea(f, psd);
  close(area / v, 1, 0.03, 'area = variance');
  // flat level: mean of the interior bins against σ²/(fs/2)
  let lvl = 0;
  for (let k = 10; k < 500; k++) lvl += psd[k];
  lvl /= 490;
  close(lvl / (v / (fs / 2)), 1, 0.03, 'white level');
});

test('psdWelch: sine tone — peak at the tone frequency, area = A²/2 for every window', () => {
  const fs = 100;
  const f0 = 12.5;
  const A = 3;
  const n = 20000;
  const x = Float64Array.from({ length: n }, (_, k) => A * Math.sin((2 * Math.PI * f0 * k) / fs + 0.3));
  for (const window of ['hann', 'hamming', 'blackman', 'rect']) {
    const { f, psd } = psdWelch(x, fs, { nfft: 2000, window });
    let kmax = 0;
    for (let k = 1; k < psd.length; k++) if (psd[k] > psd[kmax]) kmax = k;
    close(f[kmax], f0, fs / 2000 / 2, `${window}: peak frequency`);
    // the periodogram is evaluated on bins, so the area is the plain sum × Δf
    let area = 0;
    for (let k = 0; k < psd.length; k++) area += psd[k];
    area *= fs / 2000;
    close(area, (A * A) / 2, 1e-6 * A * A, `${window}: area`);
  }
});

test('psdWelch: non power-of-two nfft, odd nfft, overlap and short-record handling', () => {
  const rnd = lcg(9);
  const x = Float64Array.from({ length: 5000 }, () => rnd() - 0.5);
  const v = variance(x);
  for (const [nfft, overlap] of [[500, 0.5], [333, 0.75], [512, 0]]) {
    const { f, psd } = psdWelch(x, 4, { nfft, overlap });
    assert.equal(f.length, Math.floor(nfft / 2) + 1);
    let area = 0;
    for (let k = 0; k < psd.length; k++) area += psd[k];
    area *= 4 / nfft;
    close(area / v, 1, 0.08, `nfft=${nfft}`);
  }
  // record shorter than nfft: a single segment of the full length
  const s = psdWelch(x.subarray(0, 300), 4, { nfft: 1024 });
  assert.equal(s.nSegments, 1);
  assert.equal(s.f.length, 151);
  assert.throws(() => psdWelch(x, 0), RangeError);
  assert.throws(() => psdWelch(x, 4, { overlap: 1 }), RangeError);
  assert.throws(() => psdWelch([1], 4), RangeError);
});

test('psdWelch: detrend none keeps the mean as a DC line (mean² = ∫ of the DC lobe)', () => {
  const x = new Float64Array(4096).fill(2);
  const a = psdWelch(x, 1, { nfft: 256, detrend: 'none', window: 'rect' });
  // rectangular window: all power in bin 0, psd[0]·Δf = mean²
  close(a.psd[0] / 256, 4, 1e-10);
  const b = psdWelch(x, 1, { nfft: 256 });
  for (const p of b.psd) close(p, 0, 1e-20);
});

test('spectralMoments: box spectrum S = 2 on [0.1, 0.3] Hz — m_n = 2 (0.3^(n+1) − 0.1^(n+1))/(n+1)', () => {
  const n = 2001;
  const f = new Float64Array(n);
  const S = new Float64Array(n).fill(2);
  for (let i = 0; i < n; i++) f[i] = 0.1 + (0.2 * i) / (n - 1);
  const m = spectralMoments(f, S);
  assert.deepEqual(Object.keys(m), ['m0', 'm1', 'm2', 'm4']);
  close(m.m0, 0.4, 1e-13);
  close(m.m1, 0.3 ** 2 - 0.1 ** 2, 1e-13);
  close(m.m2, (2 * (0.3 ** 3 - 0.1 ** 3)) / 3, 1e-8);
  close(m.m4, (2 * (0.3 ** 5 - 0.1 ** 5)) / 5, 1e-9);
  // derived periods: Tz = sqrt(m0/m2), T1 = m0/m1
  close(Math.sqrt(m.m0 / m.m2), Math.sqrt(0.4 / ((2 * 0.026) / 3)), 1e-6);
  const neg = spectralMoments(f, S, [-1, 3]);
  close(neg.m_1, 2 * Math.log(3), 1e-6);
  close(neg.m3, (2 * (0.3 ** 4 - 0.1 ** 4)) / 4, 1e-9);
  assert.throws(() => spectralMoments([0.2, 0.1], [1, 1]), RangeError);
  assert.throws(() => spectralMoments([0.1, 0.2], [1]), RangeError);
});

test('synthesizeFromSpectrum: deterministic per seed, variance = m0, spectrum recovered by psdWelch', () => {
  // narrow-banded Gaussian-shaped spectrum centred on 0.1 Hz, irregular frequency grid to avoid repetition
  const nf = 300;
  const f = new Float64Array(nf);
  const S = new Float64Array(nf);
  const rnd = lcg(31);
  let fk = 0.04;
  for (let i = 0; i < nf; i++) {
    fk += 0.0004 * (0.5 + rnd());
    f[i] = fk;
    S[i] = 5 * Math.exp(-(((fk - 0.1) / 0.015) ** 2));
  }
  const m0 = spectralMoments(f, S, [0]).m0;
  const a = synthesizeFromSpectrum(f, S, { duration: 36000, dt: 0.5, seed: 11 });
  assert.equal(a.t.length, 72001);
  assert.equal(a.t[0], 0);
  close(a.t[72000], 36000, 1e-9);
  const b = synthesizeFromSpectrum(f, S, { duration: 36000, dt: 0.5, seed: 11 });
  assert.deepEqual(a.x, b.x, 'same seed -> identical series');
  const c = synthesizeFromSpectrum(f, S, { duration: 36000, dt: 0.5, seed: 12 });
  assert.notEqual(a.x[100], c.x[100]);
  close(variance(a.x) / m0, 1, 0.02, 'variance = m0');
  // significant height from the series: 4σ = 4 sqrt(m0)
  close((4 * Math.sqrt(variance(c.x))) / (4 * Math.sqrt(m0)), 1, 0.02);
  const { f: fw, psd } = psdWelch(a.x, 2, { nfft: 4096 });
  let kmax = 0;
  for (let k = 1; k < psd.length; k++) if (psd[k] > psd[kmax]) kmax = k;
  close(fw[kmax], 0.1, 0.004, 'spectral peak');
  close(trapzArea(fw, psd) / m0, 1, 0.03, 'area of the recovered spectrum');
  // spot check against the direct sum at one time instant is implied by determinism; check a single component
  const one = synthesizeFromSpectrum([0.1, 0.2, 0.3], [0, 8, 0], { duration: 100, dt: 0.01, seed: 3 });
  // amplitude sqrt(2·S·Δf) = sqrt(2·8·0.1)
  let peak = 0;
  for (const v of one.x) peak = Math.max(peak, Math.abs(v));
  close(peak, Math.sqrt(1.6), 1e-4);
});

test('synthesizeFromSpectrum: Rayleigh amplitudes keep the expected variance; argument checks', () => {
  const nf = 400;
  const f = Float64Array.from({ length: nf }, (_, i) => 0.05 + 0.00037 * i);
  const S = Float64Array.from(f, (v) => Math.exp(-(((v - 0.12) / 0.03) ** 2)));
  const m0 = spectralMoments(f, S, [0]).m0;
  let acc = 0;
  for (let seed = 1; seed <= 8; seed++) {
    const r = synthesizeFromSpectrum(f, S, { duration: 2700, dt: 1, seed, randomAmplitude: true });
    acc += variance(r.x);
  }
  close(acc / 8 / m0, 1, 0.1);
  assert.throws(() => synthesizeFromSpectrum(f, S, { duration: 0, dt: 1 }), RangeError);
  assert.throws(() => synthesizeFromSpectrum(f, S, { duration: 10, dt: -1 }), RangeError);
  assert.throws(() => synthesizeFromSpectrum([0.1, 0.1], [1, 1], { duration: 10, dt: 1 }), RangeError);
  assert.throws(() => synthesizeFromSpectrum([0.1, 0.2], [1, -1], { duration: 10, dt: 1 }), RangeError);
});
