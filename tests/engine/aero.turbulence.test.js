import test from 'node:test';
import assert from 'node:assert/strict';
import {
  TURBULENCE_CLASSES, WIND_CLASSES, NWP_EXPONENT, turbulenceScale, ntmSigma, etmSigma, ewm, eogGust, eogSpeed, edc, ecd,
  ewsIncrement, eventSeries, kaimalParams, kaimalSpectrum, vonKarmanSpectrum, coherenceIec, synthesizeSeries,
  synthesizeWindField, rotorEffectiveWind, rotorAveragingFilter, rotorEffectiveTurbulence, rotorEffectiveFromPoint,
  offshoreTurbulence, peakFactor, gustFactor,
} from '../../src/engine/aero/turbulence.js';
import { psdWelch } from '../../src/engine/math/fft.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${tol})`);
const variance = (x) => { let m = 0; for (const v of x) m += v; m /= x.length; let s = 0; for (const v of x) s += (v - m) ** 2; return s / x.length; };
const corr = (a, b) => {
  let ma = 0; let mb = 0; for (let i = 0; i < a.length; i++) { ma += a[i]; mb += b[i]; } ma /= a.length; mb /= a.length;
  let sab = 0; let saa = 0; let sbb = 0;
  for (let i = 0; i < a.length; i++) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) ** 2; sbb += (b[i] - mb) ** 2; }
  return sab / Math.sqrt(saa * sbb);
};
// integral of f(x) on a log grid (trapezoid in ln x)
function logIntegral(f, lo, hi, n = 20000) {
  let s = 0; let xp = lo; let fp = f(lo) * lo;
  for (let i = 1; i <= n; i++) { const x = lo * Math.pow(hi / lo, i / n); const fx = f(x) * x; s += 0.5 * (fx + fp) * Math.log(x / xp); xp = x; fp = fx; }
  return s;
}

test('IEC 61400-1 class tables: Iref A+ 0.18, A 0.16, B 0.14, C 0.12; Vref 50 / 42.5 / 37.5; shear exponents', () => {
  assert.deepEqual(TURBULENCE_CLASSES, { 'A+': 0.18, A: 0.16, B: 0.14, C: 0.12 });
  assert.deepEqual(WIND_CLASSES, { I: 50, II: 42.5, III: 37.5 });
  assert.equal(NWP_EXPONENT.onshore, 0.2); assert.equal(NWP_EXPONENT.offshore, 0.14);
});

test('turbulence scale: 0.7 z below 60 m, 42 m above', () => {
  close(turbulenceScale(30), 21, 1e-12); assert.equal(turbulenceScale(60), 42); assert.equal(turbulenceScale(150), 42);
});

test('NTM by hand: Iref 0.16, Vhub 15 m/s -> sigma1 = 0.16 (0.75*15 + 5.6) = 2.696 m/s', () => {
  close(ntmSigma({ Vhub: 15, Iref: 0.16 }), 2.696, 1e-12);
  close(ntmSigma({ Vhub: 15, Iref: 0.16 }) / 15, 0.17973, 1e-5);
  assert.throws(() => ntmSigma({ Vhub: 10, Iref: 1.2 }), RangeError);
});

test('ETM by hand: c = 2, Iref 0.16, Vave 10, Vhub 15 -> 2*0.16*(0.072*8*3.5 + 10) = 3.84512 m/s', () => {
  close(etmSigma({ Vhub: 15, Iref: 0.16, Vave: 10 }), 3.84512, 1e-10);
  assert.ok(etmSigma({ Vhub: 15, Iref: 0.16, Vave: 10 }) > ntmSigma({ Vhub: 15, Iref: 0.16 }));
});

test('EWM: Ve50 = 1.4 Vref (z/zhub)^0.11, Ve1 = 0.8 Ve50, turbulent V50 = Vref with sigma = 0.11 Vhub', () => {
  const h = ewm({ Vref: 50, z: 90, hubHeight: 90 });
  close(h.Ve50, 70, 1e-12); close(h.Ve1, 56, 1e-12); close(h.V50, 50, 1e-12); close(h.V1, 40, 1e-12); close(h.sigma50, 5.5, 1e-12);
  close(ewm({ Vref: 50, z: 45, hubHeight: 90 }).Ve50, 70 * Math.pow(0.5, 0.11), 1e-12);
  close(ewm({ Vref: 50, z: 45, hubHeight: 90 }).Ve50, 64.861, 2e-3);
});

test('EOG: gust magnitude by hand and Mexican-hat extrema (T = 10.5 s)', () => {
  const s1 = ntmSigma({ Vhub: 25, Iref: 0.16 });                 // 3.896
  // min{1.35 (56 - 25) = 41.85 ; 3.3 * 3.896 / (1 + 0.1*126/42) = 9.8898}
  const vg = eogGust({ Vhub: 25, sigma1: s1, D: 126, lambda1: 42, Ve1: 56 });
  close(vg, 3.3 * 3.896 / 1.3, 1e-9); close(vg, 9.8898, 1e-4);
  close(eogGust({ Vhub: 55, sigma1: s1, D: 126, lambda1: 42, Ve1: 56 }), 1.35, 1e-12);
  const o = { Vhub: 25, Vgust: vg };
  assert.equal(eogSpeed(-1, o), 25); assert.equal(eogSpeed(11, o), 25);
  close(eogSpeed(0, o), 25, 1e-12); close(eogSpeed(10.5, o), 25, 1e-9);
  close(eogSpeed(5.25, o), 25 + 0.74 * vg, 1e-12);                // peak: sin(3pi/2) = -1, 1 - cos(pi) = 2
  // the peak is the maximum, the shape is symmetric about T/2 and has two equal dips before and after
  let mx = -Infinity; let mn = Infinity; let tMn = 0;
  for (let i = 0; i <= 10500; i++) {
    const t = i / 1000; const v = eogSpeed(t, o);
    close(v, eogSpeed(10.5 - t, o), 1e-9);
    if (v > mx) mx = v;
    if (v < mn) { mn = v; tMn = t; }
  }
  close(mx, 25 + 0.74 * vg, 1e-9);
  // dips: stationary points of sin(3 pi x)(1 - cos 2 pi x), 3 cos(3 pi x)(1 - cos 2 pi x) + 2 sin(3 pi x) sin(2 pi x) = 0,
  // at x = t/T = 0.2332 (and 1 - x): depth 0.37 * 0.724 = 0.268 Vgust
  close((25 - mn) / vg, 0.37 * Math.sin(3 * Math.PI * 0.2332) * (1 - Math.cos(2 * Math.PI * 0.2332)), 1e-4);
  close(tMn / 10.5, 0.2332, 2e-3);
  close((25 - mn) / vg, 0.268, 1e-3);
});

test('EDC and ECD by hand', () => {
  // theta_e = 4 atan(2.096 / (10 * 1.3)) = 4 * 9.159 deg = 36.64 deg
  const d = edc({ Vhub: 10, sigma1: 2.096, D: 126, lambda1: 42 });
  close(d.thetaEDeg, 4 * Math.atan(2.096 / 13) * 180 / Math.PI, 1e-12); close(d.thetaEDeg, 36.64, 0.01);
  close(d.directionDeg(3), 0.5 * d.thetaEDeg, 1e-12); close(d.directionDeg(6), d.thetaEDeg, 1e-12); assert.equal(d.directionDeg(-1), 0);
  const g = ecd({ Vhub: 12 });
  assert.equal(g.Vcg, 15); assert.equal(g.T, 10); close(g.thetaCgDeg, 60, 1e-12);
  close(g.speed(5), 12 + 7.5, 1e-12); close(g.speed(10), 27, 1e-12); close(g.speed(50), 27, 1e-12);
  close(g.directionDeg(5), 30, 1e-12);
  assert.equal(ecd({ Vhub: 3 }).thetaCgDeg, 180);
});

test('EWS increment by hand at the blade tip and mid-event', () => {
  // (D/2)/D * (2.5 + 0.2*6.4*sigma1*(D/Lambda1)^0.25) * (1 - cos(pi)) = 2.5 + 1.28 * 2 * 3^0.25
  close(ewsIncrement({ offset: 63, t: 6, sigma1: 2, D: 126, lambda1: 42 }), 2.5 + 1.28 * 2 * Math.pow(3, 0.25), 1e-12);
  assert.equal(ewsIncrement({ offset: 63, t: 13, sigma1: 2, D: 126, lambda1: 42 }), 0);
  close(ewsIncrement({ offset: -63, t: 6, sigma1: 2, D: 126, lambda1: 42 }), -(2.5 + 1.28 * 2 * Math.pow(3, 0.25)), 1e-12);
  const s = eventSeries({ fn: (t) => eogSpeed(t, { Vhub: 10, Vgust: 5 }), dt: 0.25, duration: 30, tStart: 10 });
  assert.equal(s.v[0], 10); close(s.v[61], 10 + 0.74 * 5, 1e-9);
});

test('Kaimal parameters and spectrum: scales 8.1/2.7/0.66 Lambda1, sigma ratios 1/0.8/0.5, integral equals sigma^2', () => {
  const k = kaimalParams({ hubHeight: 90, sigma1: 2 });
  close(k.u.L, 340.2, 1e-9); close(k.v.L, 113.4, 1e-9); close(k.w.L, 27.72, 1e-9); close(k.Lc, 340.2, 1e-9);
  assert.equal(k.v.sigma, 1.6); assert.equal(k.w.sigma, 1.0);
  for (const c of [k.u, k.v, k.w]) {
    const I = logIntegral((f) => kaimalSpectrum(f, { Vhub: 10, sigma: c.sigma, L: c.L }), 1e-7, 1e4);
    rel(I, c.sigma ** 2, 2e-3);
  }
  // f S / sigma^2 = 4 n / (1 + 6 n)^(5/3) at n = f L / V = 1
  close(1 * kaimalSpectrum(10 / 340.2, { Vhub: 10, sigma: 2, L: 340.2 }) * (10 / 340.2) / 4, 4 / Math.pow(7, 5 / 3), 1e-12);
});

test('von Karman isotropic spectra integrate to sigma^2 and satisfy the isotropy relation at low frequency', () => {
  for (const comp of ['u', 'v']) {
    const I = logIntegral((f) => vonKarmanSpectrum(f, { V: 12, sigma: 1.5, L: 150, component: comp }), 1e-7, 1e4);
    rel(I, 2.25, 3e-3, comp);
  }
  // S_v(0) = S_u(0)/2 for isotropic turbulence; -5/3 inertial-range ratio S_v/S_u -> 4/3
  rel(vonKarmanSpectrum(0, { V: 12, sigma: 1.5, L: 150, component: 'v' }), 0.5 * vonKarmanSpectrum(0, { V: 12, sigma: 1.5, L: 150 }), 1e-12);
  rel(vonKarmanSpectrum(50, { V: 12, sigma: 1.5, L: 150, component: 'w' }) / vonKarmanSpectrum(50, { V: 12, sigma: 1.5, L: 150 }), 4 / 3, 1e-3);
});

test('IEC coherence by hand and its decay with separation and frequency', () => {
  // f = 0.1 Hz, r = 20 m, V = 10 m/s, Lc = 340.2 m: exp(-12 sqrt(0.2^2 + (0.12*20/340.2)^2)) = 0.09058
  close(coherenceIec(0.1, 20, { Vhub: 10, Lc: 340.2 }), Math.exp(-12 * Math.sqrt(0.04 + (2.4 / 340.2) ** 2)), 1e-15);
  close(coherenceIec(0.1, 20, { Vhub: 10, Lc: 340.2 }), 0.09058, 1e-4);
  assert.equal(coherenceIec(0.5, 0, { Vhub: 10, Lc: 340.2 }), 1);
  assert.ok(coherenceIec(0.05, 60, { Vhub: 10, Lc: 340.2 }) < coherenceIec(0.05, 20, { Vhub: 10, Lc: 340.2 }));
  assert.ok(coherenceIec(0.2, 20, { Vhub: 10, Lc: 340.2 }) < coherenceIec(0.05, 20, { Vhub: 10, Lc: 340.2 }));
});

test('synthesizeSeries: variance equals the sum of S df exactly; Welch spectrum reproduces the Kaimal target', () => {
  const p = { Vhub: 10, sigma: 1.8, L: 340.2 };
  const s = synthesizeSeries({ spectrum: (f) => kaimalSpectrum(f, p), duration: 4096, dt: 0.25, seed: 7 });
  assert.equal(s.x.length, 16384);
  rel(variance(s.x), s.variance, 1e-9);
  let m = 0; for (const v of s.x) m += v;
  close(m / s.x.length, 0, 1e-9);
  // band-limited analytic variance: integral of the Kaimal spectrum from df/2 to Nyquist (closed form)
  const cum = (f) => p.sigma ** 2 * (1 - Math.pow(1 + 6 * f * p.L / p.Vhub, -2 / 3));
  rel(s.variance, cum(2) - cum(0.5 / 4096), 0.02);
  const w = psdWelch(s.x, 4, { nfft: 1024 });
  // band-averaged comparison over a decade and a half
  for (const [f0, f1] of [[0.02, 0.05], [0.05, 0.15], [0.15, 0.5], [0.5, 1.5]]) {
    let a = 0; let b = 0; let n = 0;
    for (let i = 0; i < w.f.length; i++) if (w.f[i] >= f0 && w.f[i] < f1) { a += w.psd[i]; b += kaimalSpectrum(w.f[i], p); n++; }
    assert.ok(n > 5);
    rel(a / n, b / n, 0.15, `band ${f0}-${f1} Hz`);
  }
  // seeded: reproducible, and different seeds differ
  const s2 = synthesizeSeries({ spectrum: (f) => kaimalSpectrum(f, p), duration: 4096, dt: 0.25, seed: 7 });
  assert.equal(s2.x[1234], s.x[1234]);
  const s3 = synthesizeSeries({ spectrum: (f) => kaimalSpectrum(f, p), duration: 4096, dt: 0.25, seed: 8 });
  assert.notEqual(s3.x[1234], s.x[1234]);
});

test('Veers synthesis on a 5x5 grid: target variance, Kaimal spectrum, and correlation decaying with separation as predicted by the coherence model', () => {
  const Vhub = 10; const sigma1 = ntmSigma({ Vhub, Iref: 0.14 });
  const fld = synthesizeWindField({ Vhub, sigma1, hubHeight: 90, width: 120, height: 120, ny: 5, nz: 5, duration: 2048, dt: 0.25, seed: 11 });
  assert.equal(fld.u.length, 25); assert.equal(fld.points.length, 25);
  // first point has unit Cholesky amplitude: its variance is exactly the (scaled) band variance = sigma1^2
  rel(variance(fld.u[0]), sigma1 ** 2, 1e-9);
  let mv = 0; for (const u of fld.u) mv += variance(u) / 25;
  rel(mv, sigma1 ** 2, 0.08, 'mean point variance');
  assert.ok(fld.scale > 1 && fld.scale < 1.12, `scale ${fld.scale}`);
  // spectrum at the centre point
  const w = psdWelch(fld.u[12], 4, { nfft: 512 });
  let a = 0; let b = 0;
  for (let i = 0; i < w.f.length; i++) if (w.f[i] >= 0.03 && w.f[i] < 1) { a += w.psd[i]; b += fld.scale ** 2 * kaimalSpectrum(w.f[i], { Vhub, sigma: sigma1, L: 340.2 }); }
  rel(a, b, 0.2);
  // predicted correlation rho(r) = sum S Coh df / sum S df
  const df = 1 / 2048;
  const rho = (r) => { let n = 0; let d = 0; for (let m = 1; m < 4096; m++) { const S = kaimalSpectrum(m * df, { Vhub, sigma: sigma1, L: 340.2 }); n += S * coherenceIec(m * df, r, { Vhub, Lc: 340.2 }); d += S; } return n / d; };
  // average the sample correlation over all pairs at each lateral separation
  const byR = new Map();
  for (let i = 0; i < 25; i++) for (let j = i + 1; j < 25; j++) {
    const r = Math.round(Math.hypot(fld.points[i].y - fld.points[j].y, fld.points[i].z - fld.points[j].z));
    if (![30, 60, 90, 120].includes(r)) continue;
    const e = byR.get(r) || [0, 0]; e[0] += corr(fld.u[i], fld.u[j]); e[1]++; byR.set(r, e);
  }
  const c = [30, 60, 90, 120].map((r) => byR.get(r)[0] / byR.get(r)[1]);
  assert.ok(c[0] > c[1] && c[1] > c[2] && c[2] > c[3], `correlation decays: ${c}`);
  [30, 60, 90, 120].forEach((r, k) => close(c[k], rho(r), 0.08, `rho(${r} m)`));
  // deterministic for a seed
  const f2 = synthesizeWindField({ Vhub, sigma1, hubHeight: 90, width: 120, height: 120, ny: 5, nz: 5, duration: 2048, dt: 0.25, seed: 11 });
  assert.equal(f2.u[7][100], fld.u[7][100]);
  assert.throws(() => synthesizeWindField({ Vhub, sigma1, hubHeight: 90, width: 1, height: 1, ny: 12, nz: 12, duration: 10, dt: 1, seed: 1 }), RangeError);
});

test('rotor-effective wind: disc averaging reduces the variance; grid average and averaging filter agree', () => {
  const Vhub = 12; const sigma1 = 1.8; const R = 63;
  const fld = synthesizeWindField({ Vhub, sigma1, hubHeight: 90, width: 126, height: 126, ny: 7, nz: 7, duration: 2048, dt: 0.5, seed: 5 });
  const re = rotorEffectiveWind(fld, { R });
  assert.ok(re.nPoints >= 29 && re.nPoints < 49);
  let m = 0; for (const v of re.u) m += v; m /= re.u.length;
  close(m, Vhub, 1e-9);
  const sGrid = Math.sqrt(variance(re.u));
  assert.ok(sGrid < 0.8 * sigma1);
  const rt = rotorEffectiveTurbulence({ Vhub, sigma1, hubHeight: 90, R, duration: 2048, dt: 0.5, seed: 5 });
  rel(Math.sqrt(variance(rt.u)), rt.sigmaEff, 1e-9);
  rel(rt.sigmaEff, sGrid, 0.15);
  // with shear the cubic (energy) mean exceeds the linear mean
  const re3 = rotorEffectiveWind(fld, { R, shearExp: 0.2, power: 3 });
  const re1 = rotorEffectiveWind(fld, { R, shearExp: 0.2 });
  let m3 = 0; let m1 = 0; for (let i = 0; i < re3.u.length; i++) { m3 += re3.u[i]; m1 += re1.u[i]; }
  assert.ok(m3 > m1);
});

test('rotorAveragingFilter: unity for a point, monotonic low-pass for a disc, high-frequency roll-off', () => {
  const f = [0, 0.01, 0.03, 0.1, 0.3, 1];
  const h0 = rotorAveragingFilter(f, { R: 0, Vhub: 10, Lc: 340.2 });
  for (const v of h0) close(v, 1, 1e-12);
  const h = rotorAveragingFilter(f, { R: 63, Vhub: 10, Lc: 340.2 });
  for (let i = 1; i < h.length; i++) assert.ok(h[i] < h[i - 1]);
  assert.ok(h[0] < 1 && h[0] > 0.6);
  assert.ok(h[5] < 0.03);
  // filtering a point series keeps the mean and reduces the variance by the filter-weighted spectrum
  const p = synthesizeSeries({ spectrum: (fr) => kaimalSpectrum(fr, { Vhub: 10, sigma: 1.5, L: 340.2 }), duration: 1024, dt: 0.5, seed: 3 });
  const u = Float64Array.from(p.x, (v) => v + 10);
  const r = rotorEffectiveFromPoint({ u, dt: 0.5, R: 63, Vhub: 10, Lc: 340.2 });
  let m = 0; for (const v of r) m += v;
  close(m / r.length, 10, 1e-9);
  const df = 1 / 1024; const fr = []; for (let k = 1; k < 1024; k++) fr.push(k * df);
  const H2 = rotorAveragingFilter(fr, { R: 63, Vhub: 10, Lc: 340.2 });
  let expect = 0; fr.forEach((x, k) => { expect += kaimalSpectrum(x, { Vhub: 10, sigma: 1.5, L: 340.2 }) * H2[k] * df; });
  rel(variance(r), expect, 1e-6);
});

test('offshore turbulence (IEC 61400-3 Charnock roughness): self-consistent roughness and sigma1', () => {
  const o = offshoreTurbulence({ Vhub: 15, hubHeight: 100, I15: 0.06 });
  // z0 = Ac u*^2 / g with u* = kappa V / ln(z/z0)
  rel(o.z0, 0.011 * o.uStar ** 2 / 9.80665, 1e-9);
  rel(o.uStar, 0.4 * 15 / Math.log(100 / o.z0), 1e-12);
  close(o.sigma1, 15 / Math.log(100 / o.z0) + 1.28 * 1.44 * 0.06, 1e-12);
  assert.ok(o.z0 > 1e-5 && o.z0 < 5e-3);
  // sea roughness grows with wind speed, so the roughness-driven part of the intensity rises with wind
  assert.ok(offshoreTurbulence({ Vhub: 25, hubHeight: 100, I15: 0 }).TI > offshoreTurbulence({ Vhub: 10, hubHeight: 100, I15: 0 }).TI);
});

test('peak factor (Davenport) and gust factor by hand', () => {
  // nu T = 100: sqrt(2 ln 100) = 3.0349, + 0.5772/3.0349 = 3.2251
  close(peakFactor(0.5, 200), 3.2251, 1e-4);
  close(gustFactor(0.12, 3.5), 1.42, 1e-12);
  assert.throws(() => peakFactor(0.001, 10), RangeError);
});
