/**
 * Discrete Fourier transform and spectral analysis of sampled signals.
 *
 * Transform convention: X[k] = Σ_{n=0}^{N−1} x[n] e^{−2πi k n / N} (no scaling); the inverse carries the 1/N.
 * Frequencies are in Hz, sampling rates in Hz, spectra are one-sided variance densities in (signal units)²/Hz.
 *
 * References: Cooley & Tukey, Math. Comp. 19 (1965) 297–301 (radix-2 FFT); Bluestein, IEEE Trans. Audio
 * Electroacoust. 18 (1970) 451–455 (chirp-z transform for arbitrary length); Welch, IEEE Trans. Audio
 * Electroacoust. 15 (1967) 70–73 (averaged modified periodograms); Harris, Proc. IEEE 66 (1978) 51–83 (windows);
 * DNV-RP-C205 (2019) §3.3.2 (random-phase superposition of harmonic wave components, spectral moments §3.5.3).
 */

import { makeRng } from './rng.js';

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

const twiddleCache = new Map();

/** cos/sin tables for a radix-2 transform of length n: w[k] = e^{−2πik/n}, k < n/2. */
function twiddles(n) {
  let tw = twiddleCache.get(n);
  if (!tw) {
    const half = n >> 1;
    const c = new Float64Array(half);
    const s = new Float64Array(half);
    for (let k = 0; k < half; k++) {
      const a = (-2 * Math.PI * k) / n;
      c[k] = Math.cos(a);
      s[k] = Math.sin(a);
    }
    tw = { c, s };
    if (twiddleCache.size > 24) twiddleCache.clear();
    twiddleCache.set(n, tw);
  }
  return tw;
}

/** In-place iterative radix-2 forward transform; n must be a power of two. */
function fftRadix2(re, im) {
  const n = re.length;
  if (n < 2) return;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i];
      re[i] = re[j];
      re[j] = t;
      t = im[i];
      im[i] = im[j];
      im[j] = t;
    }
  }
  const { c, s } = twiddles(n);
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1;
    const step = n / len;
    for (let i = 0; i < n; i += len) {
      for (let k = 0, q = 0; k < half; k++, q += step) {
        const wr = c[q];
        const wi = s[q];
        const a = i + k;
        const b = a + half;
        const xr = re[b] * wr - im[b] * wi;
        const xi = re[b] * wi + im[b] * wr;
        re[b] = re[a] - xr;
        im[b] = im[a] - xi;
        re[a] += xr;
        im[a] += xi;
      }
    }
  }
}

/** Forward transform of arbitrary length by Bluestein's chirp-z algorithm (three radix-2 transforms). */
function fftBluestein(re, im) {
  const n = re.length;
  let m = 1;
  while (m < 2 * n - 1) m <<= 1;
  const wr = new Float64Array(n);
  const wi = new Float64Array(n);
  // chirp e^{−iπk²/n}; k² is reduced modulo 2n to keep the phase accurate for large k
  for (let k = 0; k < n; k++) {
    const k2 = (k * k) % (2 * n);
    const a = (Math.PI * k2) / n;
    wr[k] = Math.cos(a);
    wi[k] = -Math.sin(a);
  }
  const ar = new Float64Array(m);
  const ai = new Float64Array(m);
  const br = new Float64Array(m);
  const bi = new Float64Array(m);
  for (let k = 0; k < n; k++) {
    ar[k] = re[k] * wr[k] - im[k] * wi[k];
    ai[k] = re[k] * wi[k] + im[k] * wr[k];
  }
  br[0] = wr[0];
  bi[0] = -wi[0];
  for (let k = 1; k < n; k++) {
    br[k] = wr[k];
    bi[k] = -wi[k];
    br[m - k] = wr[k];
    bi[m - k] = -wi[k];
  }
  fftRadix2(ar, ai);
  fftRadix2(br, bi);
  for (let k = 0; k < m; k++) {
    const pr = ar[k] * br[k] - ai[k] * bi[k];
    const pi = ar[k] * bi[k] + ai[k] * br[k];
    // conjugate so that a second forward transform performs the inverse
    ar[k] = pr;
    ai[k] = -pi;
  }
  fftRadix2(ar, ai);
  const inv = 1 / m;
  for (let k = 0; k < n; k++) {
    const cr = ar[k] * inv;
    const ci = -ai[k] * inv;
    re[k] = cr * wr[k] - ci * wi[k];
    im[k] = cr * wi[k] + ci * wr[k];
  }
}

function prepare(re, im, name) {
  if (!isVec(re) || re.length < 1) throw new RangeError(`${name}: re must be a non-empty array`);
  const n = re.length;
  if (im !== undefined && im !== null && (!isVec(im) || im.length !== n)) {
    throw new RangeError(`${name}: im must have the same length as re`);
  }
  const r = Float64Array.from(re);
  const i = im ? Float64Array.from(im) : new Float64Array(n);
  return [r, i, n];
}

function transform(r, i, n) {
  if ((n & (n - 1)) === 0) fftRadix2(r, i);
  else fftBluestein(r, i);
}

/**
 * Forward discrete Fourier transform X[k] = Σ x[n] e^{−2πikn/N} of a real or complex sequence of any length:
 * iterative radix-2 Cooley–Tukey for powers of two, Bluestein's chirp-z algorithm otherwise (O(N log N) always).
 * Inputs are not modified.
 * Validity: N >= 1.
 * @param {ArrayLike<number>} re real part, length N
 * @param {ArrayLike<number>} [im] imaginary part, length N (zero if omitted)
 * @returns {{ re: Float64Array, im: Float64Array }} transform, length N, bin k at frequency k·fs/N
 */
export function fft(re, im) {
  const [r, i, n] = prepare(re, im, 'fft');
  transform(r, i, n);
  return { re: r, im: i };
}

/**
 * Inverse discrete Fourier transform x[n] = (1/N) Σ X[k] e^{+2πikn/N}; ifft(fft(x)) = x.
 * Validity: N >= 1, any length.
 * @param {ArrayLike<number>} re real part of the spectrum, length N
 * @param {ArrayLike<number>} [im] imaginary part of the spectrum, length N (zero if omitted)
 * @returns {{ re: Float64Array, im: Float64Array }} sequence, length N
 */
export function ifft(re, im) {
  const [r, i, n] = prepare(re, im, 'ifft');
  for (let k = 0; k < n; k++) i[k] = -i[k];
  transform(r, i, n);
  const inv = 1 / n;
  for (let k = 0; k < n; k++) {
    r[k] *= inv;
    i[k] *= -inv;
  }
  return { re: r, im: i };
}

/**
 * Window coefficients (DFT-periodic form, Harris 1978): w[k], k = 0..n−1.
 * @param {'hann' | 'hamming' | 'blackman' | 'rect' | 'boxcar' | 'none'} name window name
 * @param {number} n window length
 * @returns {Float64Array} window, length n
 */
export function windowFunction(name, n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('windowFunction: n must be an integer >= 1');
  const w = new Float64Array(n);
  const a = (2 * Math.PI) / n;
  switch (name) {
    case 'hann':
    case 'hanning':
      for (let k = 0; k < n; k++) w[k] = 0.5 - 0.5 * Math.cos(a * k);
      break;
    case 'hamming':
      for (let k = 0; k < n; k++) w[k] = 0.54 - 0.46 * Math.cos(a * k);
      break;
    case 'blackman':
      for (let k = 0; k < n; k++) w[k] = 0.42 - 0.5 * Math.cos(a * k) + 0.08 * Math.cos(2 * a * k);
      break;
    case 'rect':
    case 'boxcar':
    case 'none':
      w.fill(1);
      break;
    default:
      throw new RangeError("windowFunction: name must be 'hann', 'hamming', 'blackman' or 'rect'");
  }
  if (n === 1) w[0] = 1;
  return w;
}

/**
 * One-sided power spectral density by Welch's method (averaged, windowed, overlapping periodograms).
 * Scaling is "density": ∫ psd df over 0..fs/2 equals the variance of the (detrended) signal, i.e. psd is in
 * (units of x)²/Hz. Each segment has its mean removed before windowing (set detrend: 'none' to keep it).
 * If the record is shorter than nfft a single segment of the full record length is used.
 * Validity: fs > 0, at least 2 samples, 0 <= overlap < 1; any nfft (not restricted to powers of two).
 * @param {ArrayLike<number>} x signal samples [signal units]
 * @param {number} fs sampling frequency [Hz]
 * @param {{ nfft?: number, overlap?: number, window?: string, detrend?: 'mean' | 'none' }} [opts] nfft: segment
 *   length in samples (default 1024); overlap: fraction of a segment shared with the next (default 0.5); window:
 *   'hann' (default), 'hamming', 'blackman' or 'rect'; detrend: per-segment mean removal (default 'mean')
 * @returns {{ f: Float64Array, psd: Float64Array, nSegments: number }} f: frequencies 0..fs/2 [Hz] (spacing
 *   fs/nfft); psd: one-sided density [units²/Hz]; nSegments: number of averaged segments
 */
export function psdWelch(x, fs, { nfft = 1024, overlap = 0.5, window = 'hann', detrend = 'mean' } = {}) {
  if (!isVec(x) || x.length < 2) throw new RangeError('psdWelch: x must have at least 2 samples');
  if (!(fs > 0) || !Number.isFinite(fs)) throw new RangeError('psdWelch: fs must be > 0');
  if (!Number.isInteger(nfft) || nfft < 2) throw new RangeError('psdWelch: nfft must be an integer >= 2');
  if (!(overlap >= 0 && overlap < 1)) throw new RangeError('psdWelch: overlap must be in [0, 1)');
  if (detrend !== 'mean' && detrend !== 'none') throw new RangeError("psdWelch: detrend must be 'mean' or 'none'");
  const N = x.length;
  const L = Math.min(nfft, N);
  const w = windowFunction(window, L);
  let w2 = 0;
  for (let k = 0; k < L; k++) w2 += w[k] * w[k];
  const step = Math.max(1, Math.floor(L * (1 - overlap)));
  const nSeg = Math.floor((N - L) / step) + 1;
  const nf = (L >> 1) + 1;
  const psd = new Float64Array(nf);
  const re = new Float64Array(L);
  const im = new Float64Array(L);
  for (let s = 0; s < nSeg; s++) {
    const o = s * step;
    let mean = 0;
    if (detrend === 'mean') {
      for (let k = 0; k < L; k++) mean += x[o + k];
      mean /= L;
    }
    for (let k = 0; k < L; k++) {
      re[k] = (x[o + k] - mean) * w[k];
      im[k] = 0;
    }
    transform(re, im, L);
    for (let k = 0; k < nf; k++) psd[k] += re[k] * re[k] + im[k] * im[k];
  }
  const scale = 1 / (fs * w2 * nSeg);
  const f = new Float64Array(nf);
  for (let k = 0; k < nf; k++) {
    // fold negative frequencies onto the positive ones (not DC, nor Nyquist when L is even)
    const twoSided = k === 0 || (L % 2 === 0 && k === nf - 1);
    psd[k] *= twoSided ? scale : 2 * scale;
    f[k] = (k * fs) / L;
  }
  return { f, psd, nSegments: nSeg };
}

/**
 * Spectral moments m_n = ∫ fⁿ S(f) df of a one-sided spectrum by the trapezoidal rule, with f in Hz (so that
 * the mean zero-crossing period is Tz = sqrt(m0/m2) in seconds and Hs = 4·sqrt(m0); DNV-RP-C205 §3.5.3).
 * Moments about angular frequency follow as m_n(ω) = (2π)ⁿ m_n(f).
 * Validity: f ascending and non-negative, S on the same grid. For negative orders the point f = 0 is skipped.
 * @param {ArrayLike<number>} f frequencies [Hz]
 * @param {ArrayLike<number>} S one-sided spectral density [units²/Hz]
 * @param {number[]} [orders] moment orders (default [0, 1, 2, 4])
 * @returns {Object<string, number>} moments keyed 'm0', 'm1', 'm2', 'm4', ... [units²·Hzⁿ]; a negative order −k is
 *   keyed 'm_k' (e.g. 'm_1' for the energy-period moment)
 */
export function spectralMoments(f, S, orders = [0, 1, 2, 4]) {
  if (!isVec(f) || !isVec(S) || f.length !== S.length) throw new RangeError('spectralMoments: f and S must be arrays of equal length');
  if (!Array.isArray(orders)) throw new RangeError('spectralMoments: orders must be an array of numbers');
  const n = f.length;
  for (let i = 1; i < n; i++) {
    if (!(f[i] >= f[i - 1])) throw new RangeError(`spectralMoments: f must be ascending (violated at index ${i})`);
  }
  if (n > 0 && f[0] < 0) throw new RangeError('spectralMoments: f must be non-negative');
  const out = {};
  for (const p of orders) {
    if (!Number.isFinite(p)) throw new RangeError('spectralMoments: orders must be finite numbers');
    let m = 0;
    let prev = 0;
    for (let i = 0; i < n; i++) {
      let g;
      if (p === 0) g = S[i];
      else if (f[i] === 0) g = 0;
      else g = Math.pow(f[i], p) * S[i];
      if (i > 0) m += 0.5 * (f[i] - f[i - 1]) * (g + prev);
      prev = g;
    }
    out[p < 0 ? `m_${-p}` : `m${p}`] = m;
  }
  return out;
}

/**
 * Synthesise a stationary Gaussian-like time series from a one-sided spectrum by superposition of harmonic
 * components with random phases (DNV-RP-C205 §3.3.2): x(t) = Σ a_i cos(2π f_i t + φ_i), a_i = sqrt(2 S_i Δf_i),
 * φ_i uniform in [0, 2π) from the seeded generator, Δf_i the trapezoidal band width of sample i. The variance of
 * the series over a long record equals m0 = ∫ S df. With `randomAmplitude` the amplitudes are Rayleigh
 * distributed with E[a_i²] = 2 S_i Δf_i.
 * Validity: f ascending [Hz], S >= 0. On an equally spaced frequency grid the series repeats after 1/Δf: use
 * Δf <= 1/duration or unequal spacing. Cost O(samples × components).
 * @param {ArrayLike<number>} f frequencies [Hz]
 * @param {ArrayLike<number>} S one-sided spectral density [units²/Hz]
 * @param {{ duration: number, dt: number, seed?: number, randomAmplitude?: boolean }} opts duration: record length
 *   [s]; dt: time step [s]; seed: generator seed (default 1); randomAmplitude: Rayleigh amplitudes (default false)
 * @returns {{ t: Float64Array, x: Float64Array }} t: times 0, dt, ... <= duration [s]; x: signal [units]
 */
export function synthesizeFromSpectrum(f, S, { duration, dt, seed = 1, randomAmplitude = false } = {}) {
  if (!isVec(f) || !isVec(S) || f.length !== S.length || f.length < 2) {
    throw new RangeError('synthesizeFromSpectrum: f and S must be arrays of equal length >= 2');
  }
  if (!(duration > 0) || !Number.isFinite(duration)) throw new RangeError('synthesizeFromSpectrum: duration must be > 0');
  if (!(dt > 0) || !Number.isFinite(dt)) throw new RangeError('synthesizeFromSpectrum: dt must be > 0');
  const nf = f.length;
  for (let i = 1; i < nf; i++) {
    if (!(f[i] > f[i - 1])) throw new RangeError(`synthesizeFromSpectrum: f must be strictly ascending (violated at index ${i})`);
  }
  const rng = makeRng(seed);
  const nt = Math.floor(duration / dt + 1e-9) + 1;
  const t = new Float64Array(nt);
  const x = new Float64Array(nt);
  for (let k = 0; k < nt; k++) t[k] = k * dt;
  const RESYNC = 2048;
  for (let i = 0; i < nf; i++) {
    const df = i === 0 ? 0.5 * (f[1] - f[0]) : i === nf - 1 ? 0.5 * (f[nf - 1] - f[nf - 2]) : 0.5 * (f[i + 1] - f[i - 1]);
    // draw the random numbers for every component so the phases do not depend on which S are zero
    const phi = 2 * Math.PI * rng.next();
    const ray = randomAmplitude ? -Math.log(1 - rng.next()) : 1;
    const Si = S[i];
    if (!(Si >= 0)) throw new RangeError(`synthesizeFromSpectrum: S[${i}] must be >= 0`);
    if (Si === 0) continue;
    const amp = Math.sqrt(2 * Si * df * ray);
    const w = 2 * Math.PI * f[i];
    const cd = Math.cos(w * dt);
    const sd = Math.sin(w * dt);
    // phasor rotation with periodic exact re-synchronisation (keeps round-off drift below ~1e-12)
    for (let k0 = 0; k0 < nt; k0 += RESYNC) {
      let c = amp * Math.cos(w * t[k0] + phi);
      let s = amp * Math.sin(w * t[k0] + phi);
      const kEnd = Math.min(nt, k0 + RESYNC);
      for (let k = k0; k < kEnd; k++) {
        x[k] += c;
        const cn = c * cd - s * sd;
        s = s * cd + c * sd;
        c = cn;
      }
    }
  }
  return { t, x };
}
