import test from 'node:test';
import assert from 'node:assert/strict';
import {
  spectralMoment, spectralParameters, narrowBandDamageRate, wirschingLightFactor, dirlikParameters, dirlikDamageRate,
  tovoBenasciuttiDamageRate, alpha075DamageRate, spectralDamage, bimodalDamageRate, stressSpectrum, longTermSpectralDamage,
} from '../../src/engine/fatigue/spectral.js';
import { rainflowFourPoint } from '../../src/engine/fatigue/rainflow.js';
import { minerSumCycles } from '../../src/engine/fatigue/damage.js';
import { makeRng } from '../../src/engine/math/rng.js';
import { gamma } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// Flat (band-limited white) PSD of level S0 between f1 and f2.
function band(f1, f2, S0, n = 2001) {
  const f = new Float64Array(n);
  const S = new Float64Array(n);
  for (let i = 0; i < n; i++) { f[i] = f1 + ((f2 - f1) * i) / (n - 1); S[i] = S0; }
  return { f, S };
}

// Gaussian process from a PSD: random-phase harmonics with one randomly placed frequency per bin, which avoids
// the periodicity of an equally spaced set.
function synthesize(f1, f2, S0, { nComp, duration, dt, seed }) {
  const rng = makeRng(seed);
  const df = (f2 - f1) / nComp;
  const amp = Math.sqrt(2 * S0 * df);
  const w = new Float64Array(nComp);
  const ph = new Float64Array(nComp);
  for (let k = 0; k < nComp; k++) { w[k] = 2 * Math.PI * (f1 + (k + rng.next()) * df); ph[k] = 2 * Math.PI * rng.next(); }
  const n = Math.round(duration / dt);
  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const t = i * dt;
    let s = 0;
    for (let k = 0; k < nComp; k++) s += Math.cos(w[k] * t + ph[k]);
    x[i] = amp * s;
  }
  return x;
}

test('spectral moments and bandwidth parameters of a flat band: closed form m_n = S0·(f2^(n+1) − f1^(n+1))/(n+1)', () => {
  const { f, S } = band(0.5, 1.5, 4);
  const mn = (n) => (4 * (1.5 ** (n + 1) - 0.5 ** (n + 1))) / (n + 1);
  for (const n of [0, 1, 2, 4]) close(spectralMoment(f, S, n), mn(n), 1e-6, `m${n}`);
  const p = spectralParameters(f, S);
  close(p.nu0, Math.sqrt(mn(2) / mn(0)), 1e-6);
  close(p.nuP, Math.sqrt(mn(4) / mn(2)), 1e-6);
  close(p.alpha1, mn(1) / Math.sqrt(mn(0) * mn(2)), 1e-6);
  close(p.alpha2, mn(2) / Math.sqrt(mn(0) * mn(4)), 1e-6);
  close(p.epsilon, Math.sqrt(1 - p.alpha2 ** 2), 1e-12);
});

test('narrow-band closed form, m = 3: d = ν/a·(2√2σ)³·Γ(2.5), and bilinear form against direct integration of the Rayleigh range density', () => {
  const sn1 = { m1: 3, a1: 1e12, m2: null, a2: null, sKnee: 0 };
  close(narrowBandDamageRate({ sigma: 10, rate: 0.2, sn: sn1 }), (0.2 * (2 * Math.SQRT2 * 10) ** 3 * gamma(2.5)) / 1e12, 1e-12);
  // Bilinear curve continuous at S1 = 30: a2 = a1·S1^(m2 − m1).
  const sn2 = { m1: 3, a1: 1e12, m2: 5, a2: 1e12 * 30 ** 2, sKnee: 30 };
  const sigma = 8;
  // Rayleigh density of ranges: p(S) = S/(4σ²)·exp(−S²/(8σ²)).
  let sum = 0;
  const h = 0.002;
  for (let s = h / 2; s < 400; s += h) {
    const N = s < 30 ? sn2.a2 * s ** -5 : sn2.a1 * s ** -3;
    sum += ((s / (4 * sigma * sigma)) * Math.exp((-s * s) / (8 * sigma * sigma)) * h) / N;
  }
  close(narrowBandDamageRate({ sigma, rate: 1, sn: sn2 }), sum, 1e-6);
});

test('Wirsching–Light factor: hand value for m = 3, ε = 0.6', () => {
  // a = 0.926 − 0.099 = 0.827, b = 4.761 − 2.323 = 2.438, λ = 0.827 + 0.173·0.4^2.438
  close(wirschingLightFactor(0.6, 3), 0.827 + 0.173 * 0.4 ** 2.438, 1e-12);
  assert.equal(wirschingLightFactor(0, 3), 1);
});

test('Dirlik: parameters written out from the moments, density integrates to one, closed-form damage equals integration', () => {
  const { f, S } = band(0.1, 2.0, 50);
  const p = spectralParameters(f, S);
  const d = dirlikParameters(p);
  const g = p.m2 / Math.sqrt(p.m0 * p.m4);
  const xm = (p.m1 / p.m0) * Math.sqrt(p.m2 / p.m4);
  const D1 = (2 * (xm - g * g)) / (1 + g * g);
  const R = (g - xm - D1 * D1) / (1 - g - D1 + D1 * D1);
  const D2 = (1 - g - D1 + D1 * D1) / (1 - R);
  const D3 = 1 - D1 - D2;
  const Q = (1.25 * (g - D3 - D2 * R)) / D1;
  close(d.D1, D1, 1e-12); close(d.D2, D2, 1e-12); close(d.D3, D3, 1e-12); close(d.Q, Q, 1e-12); close(d.R, R, 1e-12);
  let area = 0;
  let dmg = 0;
  const h = 0.005;
  const sn = { m1: 4, a1: 1e14, m2: null, a2: null, sKnee: 0 };
  for (let s = h / 2; s < 40 * p.sigma; s += h) { const pdf = d.pdf(s); area += pdf * h; dmg += (pdf * h * s ** 4) / 1e14; }
  close(area, 1, 1e-6);
  close(dirlikDamageRate({ moments: p, sn }), p.nuP * dmg, 1e-5);
  // Bilinear path (numerical) collapses to the closed form when both branches coincide.
  const same = { m1: 4, a1: 1e14, m2: 4, a2: 1e14, sKnee: 3 * p.sigma };
  close(dirlikDamageRate({ moments: p, sn: same }), dirlikDamageRate({ moments: p, sn }), 2e-4);
});

test('narrow spectrum: Dirlik, Tovo–Benasciutti, α0.75 and Wirsching–Light all approach the narrow-band damage', () => {
  const { f, S } = band(0.98, 1.02, 2500);
  const sn = { m1: 3, a1: 1e12, m2: null, a2: null, sKnee: 0 };
  const nb = spectralDamage({ f, S, sn, method: 'narrow-band' }).damageRate;
  close(spectralDamage({ f, S, sn, method: 'dirlik' }).damageRate, nb, 0.03, 'Dirlik');
  close(tovoBenasciuttiDamageRate({ f, S, sn }).damageRate, nb, 0.01, 'TB');
  close(alpha075DamageRate({ f, S, sn }), nb, 0.01, 'alpha075');
  close(spectralDamage({ f, S, sn, method: 'wirsching-light' }).damageRate, nb, 0.01, 'WL');
});

test('narrow-band damage against rainflow counting of a synthesized narrow-band Gaussian process (m = 3)', () => {
  const f1 = 0.95;
  const f2 = 1.05;
  const S0 = 1000; // variance 100 → σ = 10
  const duration = 12000;
  const x = synthesize(f1, f2, S0, { nComp: 60, duration, dt: 0.04, seed: 7 });
  const sn = { m1: 3, a1: 1e12, m2: null, a2: null, sKnee: 0 };
  const timeDomain = minerSumCycles(rainflowFourPoint(x, { residue: 'half' }).cycles, sn);
  const { f, S } = band(f1, f2, S0);
  const nb = spectralDamage({ f, S, sn, duration, method: 'narrow-band' });
  close(timeDomain, nb.damage, 0.06, 'rainflow vs narrow-band');
  // Cycle count: ν0·T.
  close(rainflowFourPoint(x, { residue: 'half' }).cycles.reduce((a, c) => a + c.count, 0), nb.cycles, 0.02);
});

test('broad-band spectrum: Dirlik within 10 % of time-domain rainflow damage (Dirlik 1985 reports a few per cent); narrow-band is conservative', () => {
  const f1 = 0.1;
  const f2 = 2.0;
  const S0 = 50;
  const duration = 6000;
  const sn = { m1: 3, a1: 1e12, m2: null, a2: null, sKnee: 0 };
  const { f, S } = band(f1, f2, S0);
  const par = spectralParameters(f, S);
  assert.ok(par.alpha2 < 0.8, 'spectrum must be broad-banded');
  let timeDomain = 0;
  for (const seed of [11, 12]) {
    const x = synthesize(f1, f2, S0, { nComp: 160, duration, dt: 0.025, seed });
    timeDomain += minerSumCycles(rainflowFourPoint(x, { residue: 'half' }).cycles, sn) / 2;
  }
  const dk = spectralDamage({ f, S, sn, duration, method: 'dirlik' }).damage;
  const nb = spectralDamage({ f, S, sn, duration, method: 'narrow-band' }).damage;
  const tb = spectralDamage({ f, S, sn, duration, method: 'tovo-benasciutti' }).damage;
  const wl = spectralDamage({ f, S, sn, duration, method: 'wirsching-light' }).damage;
  close(dk, timeDomain, 0.10, 'Dirlik vs rainflow');
  close(tb, timeDomain, 0.10, 'Tovo-Benasciutti vs rainflow');
  assert.ok(nb > timeDomain, 'narrow-band must be conservative');
  assert.ok(wl < nb && wl > 0.6 * nb);
});

test('Jiao–Moan / DNV-OS-E301 dual narrow-band correction: hand value and limits', () => {
  const sn = { m1: 3, a1: 1e12, m2: null, a2: null, sKnee: 0 };
  // λL = λW = 0.5, νL = 0.01, νW = 0.1, δW = 0.1:
  // νY = √0.00505 = 0.071063, νe = √(0.25e-4 + 0.25e-4) = 0.0070711,
  // bracket = 0 + √(π/4)·3·Γ(2)/Γ(2.5) = 2, ρ = (νe/νY)·2 + (νW/νY)·0.5^1.5 = 0.19901 + 0.49752 = 0.69653
  const r = bimodalDamageRate({ sigmaLF: 10, nuLF: 0.01, sigmaWF: 10, nuWF: 0.1, sn });
  close(r.nuY, 0.071063, 1e-4);
  close(r.rho, 0.69653, 1e-4);
  close(r.sigmaY, Math.sqrt(200), 1e-12);
  close(r.damageRate, r.rho * r.combinedSpectrumRate, 1e-12);
  close(r.combinedSpectrumRate, (r.nuY * (2 * Math.SQRT2 * r.sigmaY) ** 3 * gamma(2.5)) / 1e12, 1e-12);
  // Bounds: simple summation <= dual narrow-band <= combined spectrum for this case.
  const sum = bimodalDamageRate({ sigmaLF: 10, nuLF: 0.01, sigmaWF: 10, nuWF: 0.1, sn, method: 'simple-summation' }).damageRate;
  assert.ok(sum < r.damageRate && r.damageRate < r.combinedSpectrumRate);
  // No wave-frequency part: everything reduces to the low-frequency narrow-band damage.
  const lf = bimodalDamageRate({ sigmaLF: 10, nuLF: 0.01, sigmaWF: 0, nuWF: 0.1, sn });
  close(lf.damageRate, lf.rateLF, 1e-12);
});

test('long-term damage over a scatter diagram: probability weighting and |H|^m scaling', () => {
  const { f, S } = band(0.05, 0.3, 2);
  const sn = { m1: 3, a1: 1e30, m2: null, a2: null, sKnee: 0 };
  const H = Float64Array.from(f, () => 5e6);
  const year = 31557600;
  const one = longTermSpectralDamage({ states: [{ probability: 1, f, waveSpectrum: S, transfer: H }], sn, duration: year });
  const split = longTermSpectralDamage({ states: [{ probability: 0.3, f, waveSpectrum: S, transfer: H }, { probability: 0.7, f, waveSpectrum: S, transfer: H }], sn, duration: year });
  close(split.damage, one.damage, 1e-12);
  close(split.contributions[0].damage / split.damage, 0.3, 1e-12);
  const H2 = Float64Array.from(f, () => 1e7);
  const twice = longTermSpectralDamage({ states: [{ probability: 1, f, waveSpectrum: S, transfer: H2 }], sn, duration: year });
  close(twice.damage / one.damage, 8, 1e-9);
  close(one.life, year / one.damage, 1e-12);
  close(stressSpectrum([2, 3], [10, 100])[1], 30000, 1e-14);
  // Direct stress PSD gives the same answer as wave spectrum × transfer function.
  const direct = longTermSpectralDamage({ states: [{ probability: 1, f, S: stressSpectrum(S, H) }], sn, duration: year });
  close(direct.damage, one.damage, 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => spectralMoment([0, 1], [1], 0), RangeError);
  assert.throws(() => spectralDamage({ f: [0, 1], S: [1, 1], sn: { m1: 3, a1: 1 }, method: 'nope' }), /method/);
  assert.throws(() => wirschingLightFactor(1.2, 3), /epsilon/);
});
