import test from 'node:test';
import assert from 'node:assert/strict';
import * as S from '../../src/engine/metocean/spectra.js';

const g = 9.80665;
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const Hs = 4; const Tp = 10; const wp = (2 * Math.PI) / Tp;
const range = { wMin: 0.2 * wp, wMax: 12 * wp, n: 8000, tail: true };

test('Pierson–Moskowitz: m0 = Hs²/16, peak at ωp, closed-form Tz/Tp = 0.7104 and Tm01/Tp = 0.7718', () => {
  const f = (w) => S.piersonMoskowitz(w, Hs, Tp);
  const m = S.spectrumMoments(f, range);
  assert.ok(rel(m.m0, (Hs * Hs) / 16) < 1e-4);
  const p = S.spectralParameters(m);
  assert.ok(rel(p.Hm0, Hs) < 1e-4);
  // closed form for A ω⁻⁵ exp(−B ω⁻⁴): m0/m2 = 1/sqrt(π B), B = 1.25 ωp⁴  ->  Tz/Tp = (1.25 π)^(-1/4)
  assert.ok(rel(p.Tz / Tp, Math.pow(1.25 * Math.PI, -0.25)) < 1e-4);
  assert.ok(Math.abs(p.Tz / Tp - 0.7104) < 5e-4);
  assert.ok(Math.abs(p.Tm01 / Tp - 0.7718) < 5e-4);
  assert.ok(Math.abs(p.Te / Tp - 0.8572) < 5e-4);
  // peak value S(ωp) = (5/16) Hs² e^(-5/4) / ωp and location
  assert.ok(rel(f(wp), ((5 / 16) * Hs * Hs * Math.exp(-1.25)) / wp) < 1e-14);
  assert.ok(f(wp) > f(0.99 * wp) && f(wp) > f(1.01 * wp));
  assert.equal(f(0), 0);
});

test('JONSWAP (DNV-RP-C205 §3.5.5): integrates to Hs²/16 within 0.5 % for γ = 1…5 and peaks at Tp', () => {
  for (const gam of [1, 2, 3.3, 5]) {
    const f = (w) => S.jonswap(w, { Hs, Tp, gamma: gam });
    const m = S.spectrumMoments(f, range);
    assert.ok(rel(m.m0, (Hs * Hs) / 16) < 5e-3, `gamma ${gam}`);
    let best = 0; let wBest = 0;
    for (let w = 0.5 * wp; w < 2 * wp; w += wp / 2000) { const v = f(w); if (v > best) { best = v; wBest = w; } }
    assert.ok(rel((2 * Math.PI) / wBest, Tp) < 1e-3, `peak for gamma ${gam}`);
    // peak ordinate = Aγ γ × PM peak
    assert.ok(rel(f(wp), (1 - 0.287 * Math.log(gam)) * gam * S.piersonMoskowitz(wp, Hs, Tp)) < 1e-14);
  }
  assert.ok(rel(S.jonswap(0.7, { Hs, Tp, gamma: 1 }), S.piersonMoskowitz(0.7, Hs, Tp)) < 1e-15);
});

test('JONSWAP: Tz/Tp and Tm01/Tp from numerical moments match the DNV-RP-C205 §3.5.5.4 fits', () => {
  for (const gam of [1, 2, 3.3, 5, 7]) {
    const p = S.spectralParameters(S.spectrumMoments((w) => S.jonswap(w, { Hs, Tp, gamma: gam }), range));
    const r = S.jonswapPeriodRatios(gam);
    assert.ok(Math.abs(p.Tz / Tp - r.TzOverTp) < 3e-3, `Tz gamma ${gam}`);
    assert.ok(Math.abs(p.Tm01 / Tp - r.Tm01OverTp) < 3e-3, `Tm01 gamma ${gam}`);
  }
  // hand evaluation of the fit at γ = 3.3
  assert.ok(Math.abs(S.jonswapPeriodRatios(3.3).TzOverTp - (0.6673 + 0.05037 * 3.3 - 0.00623 * 10.89 + 0.0003341 * 35.937)) < 1e-15);
  assert.throws(() => S.jonswapPeriodRatios(0.5), RangeError);
});

test('JONSWAP γ selection from Tp/√Hs (DNV-RP-C205 §3.5.5.5)', () => {
  assert.equal(S.jonswapGamma(4, 7), 5); // 3.5
  assert.equal(S.jonswapGamma(4, 10), 1); // 5.0
  assert.ok(rel(S.jonswapGamma(4, 8), Math.exp(5.75 - 1.15 * 4)) < 1e-15); // 4.0 -> 3.158
  // continuous at both ends
  assert.ok(Math.abs(Math.exp(5.75 - 1.15 * 3.6) - 5) < 0.02 && Math.abs(Math.exp(5.75 - 1.15 * 5) - 1) < 1e-12);
  const auto = S.jonswap(0.8, { Hs: 4, Tp: 8 });
  assert.ok(rel(auto, S.jonswap(0.8, { Hs: 4, Tp: 8, gamma: Math.exp(1.15) })) < 1e-14);
});

test('TMA depth factor: deep limit 1, shallow limit ω² d / 2g, monotonic; reduces Hs in finite depth', () => {
  assert.ok(Math.abs(S.tmaFactor(1.0, 2000) - 1) < 1e-12);
  const w = 0.05; const d = 2;
  assert.ok(rel(S.tmaFactor(w, d), (w * w * d) / (2 * g)) < 2e-3);
  let prev = 0;
  for (let ww = 0.05; ww < 3; ww += 0.05) { const v = S.tmaFactor(ww, 20); assert.ok(v >= prev && v <= 1); prev = v; }
  // Kitaigorodskii factor from its definition k⁻³ ∂k/∂ω, by finite differences
  const kOf = (om, dd) => { let k = (om * om) / g; for (let i = 0; i < 200; i++) k = (om * om) / (g * Math.tanh(k * dd)); return k; };
  const om = 0.6; const dd = 15; const h = 1e-5;
  const fin = (kOf(om + h, dd) - kOf(om - h, dd)) / (2 * h) / kOf(om, dd) ** 3;
  const deep = ((2 * om) / g) / ((om * om) / g) ** 3;
  assert.ok(rel(S.tmaFactor(om, dd), fin / deep) < 1e-6);
  const mDeep = S.spectrumMoments((x) => S.jonswap(x, { Hs, Tp, gamma: 3.3 }), range);
  const mTma = S.spectrumMoments((x) => S.tma(x, { Hs, Tp, gamma: 3.3, d: 15 }), range);
  assert.ok(mTma.m0 < 0.8 * mDeep.m0 && mTma.m0 > 0.2 * mDeep.m0);
});

test('Bretschneider / ITTC two-parameter spectrum returns the stated Hs and period (Tp, Tz or T1)', () => {
  const r = { wMin: 0.05, wMax: 30, n: 20000, tail: true };
  const pz = S.spectralParameters(S.spectrumMoments((w) => S.bretschneider(w, { Hs, T: 8, periodType: 'Tz' }), r));
  assert.ok(rel(pz.Hm0, Hs) < 1e-6 && rel(pz.Tz, 8) < 1e-6);
  const p1 = S.spectralParameters(S.spectrumMoments((w) => S.bretschneider(w, { Hs, T: 8, periodType: 'T1' }), r));
  assert.ok(rel(p1.Hm0, Hs) < 1e-6 && rel(p1.Tm01, 8) < 1e-6);
  assert.ok(rel(S.bretschneider(0.9, { Hs, T: 8 }), S.piersonMoskowitz(0.9, Hs, 8)) < 1e-14);
  // ITTC coefficients A = 4π³Hs²/Tz⁴ = 173 Hs²/T1⁴-type check: B = 16π³/Tz⁴ and, for T1, B ≈ 691/T1⁴
  const w = 0.7;
  const B = (16 * Math.PI ** 3) / 8 ** 4;
  assert.ok(rel(S.bretschneider(w, { Hs, T: 8, periodType: 'Tz' }), ((B * Hs * Hs) / 4 / w ** 5) * Math.exp(-B / w ** 4)) < 1e-13);
  const B1 = -Math.log(S.bretschneider(1, { Hs, T: 8, periodType: 'T1' }) / S.bretschneider(1e6, { Hs, T: 8, periodType: 'T1' }) / 1e30);
  assert.ok(Math.abs(B1 * 8 ** 4 - 691) < 1.5);
});

test('Ochi–Hubble six-parameter spectrum: each part integrates to Hs_j²/16 and peaks at ωp_j; λ = 1 is PM', () => {
  const o = { Hs1: 3, Tp1: 15, lambda1: 3, Hs2: 2, Tp2: 7, lambda2: 1.5 };
  const m = S.spectrumMoments((w) => S.ochiHubble(w, o), { wMin: 0.05, wMax: 30, n: 20000 });
  assert.ok(rel(16 * m.m0, 9 + 4) < 1e-6);
  const one = (w) => S.ochiHubble(w, { ...o, Hs2: 0 });
  const w1 = (2 * Math.PI) / 15;
  assert.ok(one(w1) > one(0.995 * w1) && one(w1) > one(1.005 * w1));
  assert.ok(rel(S.ochiHubble(0.8, { Hs1: 4, Tp1: 10, lambda1: 1, Hs2: 0, Tp2: 5, lambda2: 1 }), S.piersonMoskowitz(0.8, 4, 10)) < 1e-13);
});

test('directional spreading functions integrate to one; cos² peak value 2/π', () => {
  const integ = (f) => { let s = 0; const n = 20000; for (let i = 0; i < n; i++) s += f(-Math.PI + (2 * Math.PI * (i + 0.5)) / n) * ((2 * Math.PI) / n); return s; };
  for (const n of [0, 2, 4, 7.5, 20]) assert.ok(Math.abs(integ((t) => S.spreadingCosN(t, 0.4, n)) - 1) < 1e-6);
  for (const s of [1, 5, 10, 40]) assert.ok(Math.abs(integ((t) => S.spreadingCos2s(t, -2.9, s)) - 1) < 1e-6);
  assert.ok(rel(S.spreadingCosN(0.4, 0.4, 2), 2 / Math.PI) < 1e-13);
  assert.equal(S.spreadingCosN(0.4 + 2, 0.4, 2), 0);
  // cos-2s with s = 1: D = (1/π) cos²(θ/2)
  assert.ok(rel(S.spreadingCos2s(1, 0, 1), Math.cos(0.5) ** 2 / Math.PI) < 1e-13);
});

test('tabulated moments agree with function moments; bandwidth parameters of a narrow spectrum vanish', () => {
  const f = (w) => S.jonswap(w, { Hs, Tp, gamma: 3.3 });
  const n = 6000; const om = new Float64Array(n); const sp = new Float64Array(n);
  for (let i = 0; i < n; i++) { om[i] = 0.2 * wp + (i * 11.8 * wp) / (n - 1); sp[i] = f(om[i]); }
  const a = S.tabulatedMoments(om, sp);
  const b = S.spectrumMoments(f, { wMin: 0.2 * wp, wMax: 12 * wp, n: 8000 });
  for (const key of ['mMinus1', 'm0', 'm1', 'm2', 'm4']) assert.ok(rel(a[key], b[key]) < 1e-4, key);
  const w0 = 0.8; const sg = 0.004;
  const nb = S.spectralParameters(S.spectrumMoments((w) => Math.exp(-0.5 * ((w - w0) / sg) ** 2), { wMin: 0.7, wMax: 0.9, n: 4000 }));
  assert.ok(nb.epsilon < 0.02 && nb.nu < 0.01);
  assert.ok(rel(nb.Tz, (2 * Math.PI) / w0) < 1e-4 && rel(nb.Te, (2 * Math.PI) / w0) < 1e-4);
  assert.throws(() => S.tabulatedMoments([1, 1], [0, 0]), RangeError);
});

test('Rayleigh statistics: exceedance, most probable maximum and expected maximum vs the asymptotic formula', () => {
  assert.ok(rel(S.rayleighExceedance(Hs, Hs), Math.exp(-2)) < 1e-15); // 13.5 % of waves exceed Hs
  assert.ok(rel(S.mostProbableMaxHeight(1, 1000), Math.sqrt(0.5 * Math.log(1000))) < 1e-15); // 1.858 Hs
  // the mode of the maximum satisfies N·P(H > h) = 1 to leading order
  assert.ok(rel(1000 * S.rayleighExceedance(S.mostProbableMaxHeight(2, 1000), 2), 1) < 1e-12);
  for (const [N, tol] of [[100, 1e-2], [1000, 5e-3], [10000, 3e-3]]) {
    const ex = S.expectedMaxHeight(1, N, 'exact');
    const as = S.expectedMaxHeight(1, N);
    assert.ok(rel(as, ex) < tol, `N=${N}: ${as} vs ${ex}`);
    assert.ok(ex > S.mostProbableMaxHeight(1, N));
  }
  // exact value for N = 1 is the Rayleigh mean, Hs sqrt(π/8)
  assert.ok(rel(S.expectedMaxHeight(1, 1, 'exact'), Math.sqrt(Math.PI / 8)) < 1e-9);
  // exact value for N = 2: E = 2 E1 − E[min] = Hs sqrt(π/8)(2 − 1/√2)
  assert.ok(rel(S.expectedMaxHeight(1, 2, 'exact'), Math.sqrt(Math.PI / 8) * (2 - Math.SQRT1_2)) < 1e-9);
  const m = S.maxWaveHeight({ Hs: 10, Tz: 10.8, duration: 10800 });
  assert.ok(rel(m.N, 1000) < 1e-12 && rel(m.mostProbableOverHs, 1.8585) < 1e-4);
});

test('Forristall (2000) crest distribution: Rayleigh limit, hand calculation of 2D and 3D parameters', () => {
  // vanishing steepness in deep water: α = 0.3536 = 1/√8, β = 2  ->  Rayleigh crests P = exp(−8 x²/Hs²)
  const lin = S.forristallCrest({ Hs: 1e-6, T1: 10, d: Infinity });
  assert.ok(Math.abs(lin.alpha - 0.3536) < 1e-8 && Math.abs(lin.beta - 2) < 1e-6);
  // Hs = 10 m, T1 = 10 s, d = 50 m
  const Hs1 = 10; const T1 = 10; const d = 50;
  let k = (2 * Math.PI / T1) ** 2 / g;
  for (let i = 0; i < 200; i++) k = (2 * Math.PI / T1) ** 2 / (g * Math.tanh(k * d));
  const S1 = (2 * Math.PI * Hs1) / (g * T1 * T1); // 0.06407
  const Ur = Hs1 / (k * k * d ** 3);
  const p2 = S.forristallCrest({ Hs: Hs1, T1, d, type: '2D', N: 1000 });
  const p3 = S.forristallCrest({ Hs: Hs1, T1, d, type: '3D', N: 1000 });
  assert.ok(Math.abs(S1 - 0.06407) < 1e-5);
  assert.ok(rel(p2.S1, S1) < 1e-13 && rel(p2.Ur, Ur) < 1e-9);
  assert.ok(rel(p2.alpha, 0.3536 + 0.2892 * S1 + 0.1060 * Ur) < 1e-12);
  assert.ok(rel(p2.beta, 2 - 2.1597 * S1 + 0.0968 * Ur * Ur) < 1e-12);
  assert.ok(rel(p3.alpha, 0.3536 + 0.2568 * S1 + 0.0800 * Ur) < 1e-12);
  assert.ok(rel(p3.beta, 2 - 1.7912 * S1 - 0.5302 * Ur + 0.284 * Ur * Ur) < 1e-12);
  // second-order crests exceed the linear (Rayleigh) crest at the same probability; 2D above 3D in this regime
  const linMax = 0.3536 * Hs1 * Math.sqrt(Math.log(1000));
  assert.ok(p2.mostProbableMax > p3.mostProbableMax && p3.mostProbableMax > linMax);
  assert.ok(rel(S.forristallExceedance(p3.mostProbableMax, p3), 1 / 1000) < 1e-10);
  assert.ok(p3.expectedMax > p3.mostProbableMax);
});

test('near-bed velocity (DNV-RP-F109 §3.4.3): narrow-band limit equals the linear wave of height Hs', () => {
  const d = 30; const T = 9; const w0 = (2 * Math.PI) / T; const sg = 0.002;
  const HsN = 3;
  const m0 = (HsN * HsN) / 16;
  const spec = (w) => (m0 / (sg * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * ((w - w0) / sg) ** 2);
  const r = S.nearBedVelocity(spec, { d, wMin: w0 - 0.03, wMax: w0 + 0.03, n: 2000 });
  let k = (w0 * w0) / g;
  for (let i = 0; i < 300; i++) k = (w0 * w0) / (g * Math.tanh(k * d));
  assert.ok(rel(r.Us, (Math.PI * HsN) / (T * Math.sinh(k * d))) < 1e-3);
  assert.ok(rel(r.Tu, T) < 1e-3);
});

test('near-bed velocity for JONSWAP: Tu exceeds Tz of the surface, Us falls with depth, deep water gives zero', () => {
  const f = (w) => S.jonswap(w, { Hs: 5, Tp: 11, gamma: 3.3 });
  const lim = { wMin: 0.1, wMax: 6 };
  const a = S.nearBedVelocity(f, { d: 30, ...lim });
  const b = S.nearBedVelocity(f, { d: 60, ...lim });
  const tz = 11 * S.jonswapPeriodRatios(3.3).TzOverTp;
  assert.ok(a.Tu > tz && a.Tu < 11 * 1.1);
  assert.ok(b.Us < a.Us && a.Us > 0.3 && a.Us < 2.5);
  assert.ok(S.nearBedVelocity(f, { d: 3000, ...lim }).Us < 1e-12);
  // independent quadrature (midpoint rule) of ∫ (ω/sinh kd)² S dω
  let M0 = 0; const n = 40000; const dw = (lim.wMax - lim.wMin) / n;
  for (let i = 0; i < n; i++) {
    const w = lim.wMin + (i + 0.5) * dw;
    let k = (w * w) / g / Math.sqrt(Math.tanh((w * w * 30) / g));
    for (let j = 0; j < 60; j++) { const th = Math.tanh(k * 30); k -= (k * th - (w * w) / g) / (th + k * 30 * (1 - th * th)); }
    M0 += (w / Math.sinh(k * 30)) ** 2 * f(w) * dw;
  }
  assert.ok(rel(a.Us, 2 * Math.sqrt(M0)) < 1e-5);
});

test('random sea-state synthesis is deterministic for a seed and reproduces the spectral variance', () => {
  const f = (w) => S.jonswap(w, { Hs, Tp, gamma: 3.3 });
  const o = { wMin: 0.3 * wp, wMax: 5 * wp, duration: 3600, dt: 0.5, nFreq: 300 };
  const a = S.synthesizeSeaState(f, { ...o, seed: 7 });
  const b = S.synthesizeSeaState(f, { ...o, seed: 7 });
  const c = S.synthesizeSeaState(f, { ...o, seed: 8 });
  assert.equal(a.eta.length, 7201);
  assert.deepEqual(Array.from(a.eta.slice(0, 50)), Array.from(b.eta.slice(0, 50)));
  assert.notEqual(a.eta[10], c.eta[10]);
  let v = 0; let mean = 0;
  for (let i = 0; i < a.eta.length; i++) { mean += a.eta[i]; v += a.eta[i] ** 2; }
  mean /= a.eta.length; v = v / a.eta.length - mean * mean;
  // variance equals Σ a_i²/2 (band-limited m0) up to the finite-record sampling error
  let target = 0;
  for (let i = 0; i < a.components.amp.length; i++) target += 0.5 * a.components.amp[i] ** 2;
  assert.ok(rel(target, (Hs * Hs) / 16) < 0.02);
  assert.ok(rel(v, target) < 0.06);
  // the recurrence used for speed reproduces the direct sum
  const j = 5000; let direct = 0;
  for (let i = 0; i < a.components.amp.length; i++) direct += a.components.amp[i] * Math.cos(a.components.omega[i] * a.t[j] + a.components.phase[i]);
  assert.ok(Math.abs(direct - a.eta[j]) < 1e-9);
});

test('regular-wave equivalents: energy-equivalent height Hs/√2 and design-wave period range', () => {
  const r = S.regularWaveEquivalents({ Hs: 8, Tp: 13, Tz: 10, duration: 10800 });
  assert.ok(rel(r.energy.H, 8 / Math.SQRT2) < 1e-15);
  assert.ok(rel(r.design.N, 1080) < 1e-12);
  assert.ok(rel(r.design.H, 8 * Math.sqrt(0.5 * Math.log(1080))) < 1e-14);
  assert.ok(rel(r.design.TLow, Math.sqrt(6.5 * r.design.H)) < 1e-14 && rel(r.design.THigh, Math.sqrt(11 * r.design.H)) < 1e-14);
});
