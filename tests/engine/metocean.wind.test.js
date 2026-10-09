import test from 'node:test';
import assert from 'node:assert/strict';
import * as Wd from '../../src/engine/metocean/wind.js';
import { makeRng } from '../../src/engine/math/rng.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('power-law and logarithmic profiles: reference value, round trip and hand values', () => {
  assert.equal(Wd.powerLawProfile(10, 25), 25);
  assert.ok(rel(Wd.powerLawProfile(80, 25, 10, 0.12), 25 * 8 ** 0.12) < 1e-15);
  // round trip 10 m -> 90 m -> 10 m
  const up = Wd.powerLawProfile(90, 20, 10, 0.14);
  assert.ok(rel(Wd.powerLawProfile(10, up, 90, 0.14), 20) < 1e-14);
  assert.equal(Wd.logProfile(10, 25), 25);
  assert.ok(rel(Wd.logProfile(100, 25, 10, 0.001), (25 * Math.log(1e5)) / Math.log(1e4)) < 1e-15); // 31.25
  const lu = Wd.logProfile(120, 18, 10, 0.0005);
  assert.ok(rel(Wd.logProfile(10, lu, 120, 0.0005), 18) < 1e-14);
  // friction velocity reproduces the profile: U(z) = (u*/κ) ln(z/z0)
  const us = Wd.frictionVelocity(25, 10, 0.001);
  assert.ok(rel((us / 0.4) * Math.log(100 / 0.001), Wd.logProfile(100, 25, 10, 0.001)) < 1e-14);
  // equivalent shear exponent makes the two profiles agree at both heights
  const al = Wd.equivalentShearExponent(10, 100, 0.001);
  assert.ok(rel(Wd.powerLawProfile(100, 25, 10, al), Wd.logProfile(100, 25, 10, 0.001)) < 1e-14);
  assert.throws(() => Wd.logProfile(0.0001, 25, 10, 0.001), RangeError);
});

test('Frøya profile (ISO 19901-1 Eqs. A.4–A.7): hand calculation for U0 = 30 m/s', () => {
  const U0 = 30;
  assert.equal(Wd.froyaWind({ U0 }).u, 30); // z = 10 m, t = 1 h is the reference condition
  const C = 0.0573 * Math.sqrt(1 + 0.15 * 30); // 0.13438
  assert.ok(Math.abs(C - 0.13438) < 1e-5);
  const a = Wd.froyaWind({ U0, z: 50 });
  assert.ok(rel(a.Uz, 30 * (1 + C * Math.log(5))) < 1e-14); // 36.49 m/s
  assert.ok(Math.abs(a.Uz - 36.49) < 0.01);
  const b = Wd.froyaWind({ U0, z: 10, t: 3 });
  const Iu = 0.06 * (1 + 0.043 * 30); // 0.1374
  assert.ok(rel(b.Iu, Iu) < 1e-14);
  assert.ok(rel(b.u, 30 * (1 - 0.41 * Iu * Math.log(3 / 3600))) < 1e-14); // 41.98 m/s
  assert.ok(Math.abs(b.u - 41.98) < 0.01);
  assert.ok(rel(Wd.gustFactor({ U0, z: 10, t: 3 }), b.u / 30) < 1e-14);
  const c = Wd.froyaWind({ U0, z: 50, t: 60 });
  assert.ok(rel(c.Iu, Iu * 5 ** -0.22) < 1e-14);
  assert.ok(rel(c.u, a.Uz * (1 - 0.41 * c.Iu * Math.log(60 / 3600))) < 1e-14);
  assert.throws(() => Wd.froyaWind({ U0, t: 7200 }), RangeError);
});

test('Frøya time-averaging conversion round-trips (u(z,t) -> U0 -> u(z,t)) and between conditions', () => {
  for (const U0 of [8, 20, 35, 50]) {
    for (const [z, t] of [[10, 3600], [10, 3], [30, 600], [120, 60]]) {
      const u = Wd.froyaWind({ U0, z, t }).u;
      assert.ok(rel(Wd.froyaReferenceSpeed({ u, z, t }), U0) < 1e-12);
    }
  }
  const fwd = Wd.convertWindSpeed({ u: 40, zFrom: 10, tFrom: 600, zTo: 80, tTo: 3 });
  const back = Wd.convertWindSpeed({ u: fwd.u, zFrom: 80, tFrom: 3, zTo: 10, tTo: 600 });
  assert.ok(rel(back.u, 40) < 1e-11);
  assert.ok(fwd.u > 40 && fwd.U0 < 40);
});

test('Frøya spectrum (ISO 19901-1 Eqs. A.8–A.9): hand calculation and height/speed scaling', () => {
  // U0 = 10 m/s, z = 10 m, f = 0.01 Hz: f̃ = 1.72, S = 320 / (1 + 1.72^0.468)^(5/(3·0.468))
  const hand = 320 / Math.pow(1 + Math.pow(1.72, 0.468), 5 / (3 * 0.468));
  assert.ok(rel(Wd.froyaSpectrum(0.01, { U0: 10, z: 10 }), hand) < 1e-14);
  assert.ok(Math.abs(hand - 16.77) < 0.02);
  // low-frequency plateau: S -> 320 (U0/10)² (z/10)^0.45
  assert.ok(rel(Wd.froyaSpectrum(1e-12, { U0: 20, z: 40 }), 320 * 4 * 4 ** 0.45) < 1e-3);
  // high-frequency tail falls as f^(-5/3)
  // (approached slowly because n = 0.468: local slope = −(5/3) f̃ⁿ/(1 + f̃ⁿ))
  const r = Wd.froyaSpectrum(2e6, { U0: 20, z: 40 }) / Wd.froyaSpectrum(1e6, { U0: 20, z: 40 });
  assert.ok(Math.abs(Math.log2(r) + 5 / 3) < 0.002);
  const ft = 172 * 0.1 * 4 ** (2 / 3) * 2 ** -0.75; const h = 1e-4;
  const slope = Math.log(Wd.froyaSpectrum(0.1 * (1 + h), { U0: 20, z: 40 }) / Wd.froyaSpectrum(0.1 * (1 - h), { U0: 20, z: 40 })) / Math.log((1 + h) / (1 - h));
  assert.ok(Math.abs(slope + ((5 / 3) * ft ** 0.468) / (1 + ft ** 0.468)) < 1e-6);
});

test('IEC 61400-1 normal and extreme turbulence models: hand values', () => {
  // class A at 15 m/s: σ1 = 0.16 (0.75·15 + 5.6) = 2.696 m/s
  const n = Wd.iecNtm(15, Wd.IEC_IREF.A);
  assert.ok(rel(n.sigma1, 2.696) < 1e-13 && rel(n.TI, 2.696 / 15) < 1e-13);
  // ETM, class IB (Vave = 10 m/s), Vhub = 20: σ1 = 2·0.14·(0.072·(5 + 3)(10 − 4) + 10) = 3.76768
  const e = Wd.iecEtm({ Vhub: 20, Iref: Wd.IEC_IREF.B, Vave: 0.2 * Wd.IEC_VREF.I });
  assert.ok(rel(e.sigma1, 2 * 0.14 * (0.072 * 8 * 6 + 10)) < 1e-14);
  assert.ok(e.sigma1 > Wd.iecNtm(20, 0.14).sigma1);
  assert.equal(Wd.iecTurbulenceScale(30), 21);
  assert.equal(Wd.iecTurbulenceScale(100), 42);
});

test('Kaimal spectrum (IEC 61400-1 Annex C) integrates to σ_k² for each component', () => {
  const o = { Vhub: 12, sigma1: 2.0, z: 90 };
  for (const [comp, fac] of [['u', 1], ['v', 0.8], ['w', 0.5]]) {
    // integrate in ln f from 1e-7 to 1e5 Hz (trapezoid) and add the analytic f^(-5/3) tail
    let s = 0; const n = 20000; const a = Math.log(1e-7); const b = Math.log(1e5);
    for (let i = 0; i <= n; i++) {
      const f = Math.exp(a + ((b - a) * i) / n);
      s += (i === 0 || i === n ? 0.5 : 1) * f * Wd.kaimalSpectrum(f, { ...o, component: comp }) * ((b - a) / n);
    }
    s += 1.5 * 1e5 * Wd.kaimalSpectrum(1e5, { ...o, component: comp });
    assert.ok(rel(s, (fac * 2.0) ** 2) < 1e-4, comp);
  }
  // S(0) = 4 σ² L/V with L = 8.1·42 m
  assert.ok(rel(Wd.kaimalSpectrum(0, o), (4 * 4 * 8.1 * 42) / 12) < 1e-14);
});

test('IEC 61400-1 extreme wind speed model and extreme operating gust', () => {
  const w = Wd.iecExtremeWind({ Vref: 50, zhub: 100 });
  assert.ok(rel(w.Ve50, 70) < 1e-14 && rel(w.Ve1, 56) < 1e-14 && rel(w.V1, 40) < 1e-14);
  assert.ok(rel(w.sigma1, 5.5) < 1e-14);
  assert.ok(rel(Wd.iecExtremeWind({ Vref: 50, zhub: 100, z: 50 }).Ve50, 70 * 0.5 ** 0.11) < 1e-14);
  // EOG at Vhub = 25 m/s, class IA, D = 120 m, hub 100 m
  const sigma1 = Wd.iecNtm(25, 0.16).sigma1; // 3.896
  const gust = Wd.iecEog({ Vhub: 25, Ve1: 56, sigma1, D: 120, zhub: 100 });
  const v2 = (3.3 * sigma1) / (1 + (0.1 * 120) / 42);
  assert.ok(rel(gust.Vgust, Math.min(1.35 * (56 - 25), v2)) < 1e-14); // 10.0 m/s, turbulence branch governs
  assert.ok(rel(gust.Vgust, v2) < 1e-14);
  // gust starts and ends at Vhub, dips before and after the peak; the peak is at mid-duration:
  // sin(3π/2)(1 − cos π) = −2  ->  V = Vhub + 0.74 Vgust
  const at = (t) => Wd.iecEog({ Vhub: 25, Ve1: 56, sigma1, D: 120, zhub: 100, t }).V;
  assert.ok(Math.abs(at(0) - 25) < 1e-12 && Math.abs(at(10.5) - 25) < 1e-12 && at(20) === 25);
  assert.ok(rel(at(5.25), 25 + 0.74 * gust.Vgust) < 1e-13);
  assert.ok(rel(gust.Vmax, 25 + 0.74 * gust.Vgust) < 1e-9);
  assert.ok(gust.Vmin < 25 && at(2) < 25);
});

test('wind force F = ½ ρ Cs A U² with API RP 2A shape coefficients; air density from the gas law', () => {
  assert.ok(rel(Wd.airDensity(), 101325 / (287.05 * 288.15)) < 1e-15);
  assert.ok(Math.abs(Wd.airDensity() - 1.225) < 5e-4);
  assert.ok(Wd.airDensity({ temperatureC: 30 }) < Wd.airDensity({ temperatureC: 0 }));
  assert.ok(Wd.airDensity({ temperatureC: 25, relativeHumidity: 1 }) < Wd.airDensity({ temperatureC: 25 }));
  // cylinder: ½·1.225·0.5·20 m²·(40 m/s)² = 9800 N
  const f = Wd.windForce({ U: 40, A: 20, Cs: Wd.API_SHAPE_COEFFICIENTS.cylindricalSections, rho: 1.225 });
  assert.ok(rel(f.F, 9800) < 1e-14 && rel(f.q, 980) < 1e-14);
  assert.deepEqual({ ...Wd.API_SHAPE_COEFFICIENTS }, { beams: 1.5, buildingSides: 1.5, cylindricalSections: 0.5, overallPlatformProjectedArea: 1.0 });
  // height coefficient is equivalent to using the profile speed at the centroid
  const Ch = Wd.heightCoefficient(60, 10, 0.12);
  const direct = Wd.windForce({ U: Wd.powerLawProfile(60, 40, 10, 0.12), A: 20, Cs: 1.5 }).F;
  assert.ok(rel(Wd.windForce({ U: 40, A: 20, Cs: 1.5, Ch }).F, direct) < 1e-14);
  assert.ok(rel(Wd.windForce({ U: 40, A: 20, alpha: Math.PI / 6 }).F, 0.5 * 980 * 20) < 1e-13);
});

test('Weibull wind climate: moments, power density for k = 2, mean/std round trip', () => {
  const A = 10; const k = 2;
  const w = Wd.weibullWind({ A, k, U: 8, rho: 1.225 });
  assert.ok(rel(w.mean, (A * Math.sqrt(Math.PI)) / 2) < 1e-12); // Rayleigh mean
  assert.ok(rel(w.powerDensity, 0.5 * 1.225 * 1000 * 0.75 * Math.sqrt(Math.PI)) < 1e-12); // Γ(2.5) = 3√π/4
  assert.ok(rel(w.cdf, 1 - Math.exp(-0.64)) < 1e-14);
  assert.ok(rel(w.pdf, 0.16 * Math.exp(-0.64)) < 1e-14);
  assert.ok(rel(w.median, A * Math.sqrt(Math.LN2)) < 1e-14 && rel(w.mode, A * Math.SQRT1_2) < 1e-14);
  for (const [a, kk] of [[8.5, 1.6], [11, 2.3], [6, 3.5]]) {
    const m = Wd.weibullWind({ A: a, k: kk });
    const f = Wd.weibullFromMeanStd(m.mean, m.std);
    assert.ok(rel(f.A, a) < 1e-9 && rel(f.k, kk) < 1e-9);
  }
});

test('Weibull maximum-likelihood fit recovers the parameters of a synthetic sample', () => {
  const rng = makeRng(42);
  const A = 9.5; const k = 2.2; const n = 40000;
  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) x[i] = A * Math.pow(-Math.log(1 - rng.next()), 1 / k);
  const f = Wd.weibullFitMle(x);
  // asymptotic standard errors: se(k) ≈ 0.78 k/√n = 0.0086, se(A) ≈ 1.05 A/(k √n) = 0.023
  assert.ok(Math.abs(f.k - k) < 0.035 && Math.abs(f.A - A) < 0.09, `A=${f.A} k=${f.k}`);
  assert.equal(f.n, n);
  // likelihood equation is satisfied: perturbing either parameter lowers the log-likelihood
  const ll = (a, kk) => { let s = 0; for (let i = 0; i < n; i++) s += Math.log(kk / a) + (kk - 1) * Math.log(x[i] / a) - Math.pow(x[i] / a, kk); return s; };
  assert.ok(rel(f.logLik, ll(f.A, f.k)) < 1e-10);
  assert.ok(ll(f.A * 1.01, f.k) < f.logLik && ll(f.A * 0.99, f.k) < f.logLik && ll(f.A, f.k * 1.01) < f.logLik && ll(f.A, f.k * 0.99) < f.logLik);
  assert.equal(Wd.weibullFitMle([0, 5, 7, -1, 9]).nExcluded, 2);
  assert.throws(() => Wd.weibullFitMle([3]), RangeError);
});
