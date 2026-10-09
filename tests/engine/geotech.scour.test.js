import test from 'node:test';
import assert from 'node:assert/strict';
import * as s from '../../src/engine/geotech/scour.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('Sumer et al. (1992) pile scour in waves: zero at KC = 6, S/D = 0.446 at KC = 20, 1.3 for KC -> infinity', () => {
  assert.equal(s.pileScour({ KC: 6 }).SoverD, 0);
  assert.equal(s.pileScour({ KC: 3 }).SoverD, 0);
  // 1.3 (1 - exp(-0.03·14)) = 1.3 · 0.34295
  close(s.pileScour({ KC: 20 }).SoverD, 0.4458, 2e-4);
  close(s.pileScour({ KC: 1e6 }).SoverD, 1.3, 1e-9);
  close(s.pileScour().SoverD, 1.3, 1e-14);
  close(s.keuleganCarpenter(1.5, 10, 1.2), 12.5, 1e-12);
});

test('steady current: 1.3 D mean, 2.0 D with one standard deviation (sigma = 0.7)', () => {
  close(s.pileScour({ Ucw: 1 }).SoverD, 1.3, 1e-14);
  close(s.pileScour({ Ucw: 1, nSigma: 1 }).SoverD, 2.0, 1e-14);
});

test('Sumer & Fredsoe (2001) combined waves and current: KC = 10, Ucw = 0.5 -> S/D = 0.995 (hand)', () => {
  const r = s.pileScour({ KC: 10, Ucw: 0.5 });
  // A = 0.03 + 0.75·0.5^2.6 = 0.1537; B = 6 exp(-2.35) = 0.5722
  close(r.A, 0.1537, 2e-4);
  close(r.B, 0.5722, 2e-4);
  close(r.SoverD, 0.9948, 3e-4);
  // wave-only limit and monotonic increase with current
  close(s.pileScour({ KC: 10, Ucw: 0 }).SoverD, 1.3 * (1 - Math.exp(-0.12)), 1e-12);
  assert.ok(s.pileScour({ KC: 10, Ucw: 0.7 }).SoverD > r.SoverD);
});

test('Rudolph & Bos coefficients reduce to the wave-only form at Ucw = 0', () => {
  close(s.pileScourRudolphBos({ KC: 9, Ucw: 0 }).SoverD, s.pileScour({ KC: 9 }).SoverD, 1e-12);
});

test('Breusers et al. (1977): S = 1.5 D tanh(h/D) -> 1.446 D at h/D = 2, 1.5 D in deep water', () => {
  close(s.scourBreusers({ D: 2, h: 4 }), 1.5 * 2 * Math.tanh(2), 1e-14);
  close(s.scourBreusers({ D: 2, h: 4 }) / 2, 1.446, 1e-3);
  close(s.scourBreusers({ D: 1, h: 100 }), 1.5, 1e-9);
});

test('HEC-18 CSU equation by hand: a = 2 m, y1 = 5 m, V1 = 2 m/s, circular pier, K3 = 1.1 -> 3.54 m', () => {
  const r = s.scourHec18({ a: 2, y1: 5, V1: 2 });
  close(r.Fr, 0.28563, 1e-4);
  // 2.0·1.0·1.0·1.1·2·(2.5^0.35 = 1.3781)·(0.28563^0.43 = 0.58343)
  close(r.ys, 3.5377, 5e-4);
  // skewed rectangular pier L/a = 4 at 30 deg: K2 = (cos30 + 4 sin30)^0.65
  const k = s.scourHec18({ a: 2, y1: 5, V1: 2, nose: 'square', skew: Math.PI / 6, L: 8 });
  close(k.K2, (Math.cos(Math.PI / 6) + 2) ** 0.65, 1e-12);
  // cap at 2.4 a for round-nosed aligned piers at low Froude number
  assert.equal(s.scourHec18({ a: 0.5, y1: 30, V1: 3 }).ys, 1.2);
});

test('Sheppard-Melville (HEC-18): clear-water peak at V = Vc equals live-bed start; f-functions by hand', () => {
  const o = { a: 2, y1: 6, Vc: 0.5, D50: 0.0005 };
  const below = s.scourSheppardMelville({ ...o, V1: 0.5 * (1 - 1e-9) });
  const above = s.scourSheppardMelville({ ...o, V1: 0.5 });
  close(below.ys, above.ys, 1e-6);
  close(below.f1, Math.tanh(3 ** 0.4), 1e-12);
  const ad = 4000;
  close(below.f3, ad / (0.4 * ad ** 1.2 + 10.6 * ad ** -0.13), 1e-12);
  close(above.ys, 2.5 * below.f1 * below.f3 * 2, 1e-9);
  // beyond the live-bed peak velocity the depth is 2.2 f1 a
  close(s.scourSheppardMelville({ ...o, V1: 10 }).ys, 2.2 * below.f1 * 2, 1e-12);
  assert.equal(s.scourSheppardMelville({ ...o, V1: 0.1 }).ys, 0);
});

test('time scale of scour (Sumer et al. 1992): T* by hand for current and waves; exponential development', () => {
  const c = s.scourTimeScale({ mode: 'current', D: 1, theta: 0.1, d: 0.0002, delta: 2 });
  close(c.Tstar, (2 / 2000) * 0.1 ** -2.2, 1e-12);
  close(c.Tstar, 0.15849, 1e-4);
  close(c.T, c.Tstar / Math.sqrt(9.80665 * 1.65 * 0.0002 ** 3), 1e-12);
  const w = s.scourTimeScale({ mode: 'waves', D: 1, theta: 0.2, d: 0.0002, KC: 10 });
  close(w.Tstar, 1e-6 * 50 ** 3, 1e-12);
  close(s.scourDevelopment(2, 3600, 3600), 2 * (1 - 1 / Math.E), 1e-14);
});

test('pipeline scour: S/D = 0.1 sqrt(KC) = 0.4 at KC = 16; current 0.625 exp(-0.6 e/D); onset criterion by hand', () => {
  close(s.pipelineScour({ mode: 'waves', KC: 16 }), 0.4, 1e-14);
  close(s.pipelineScour({ mode: 'current' }), 0.625, 1e-14);
  close(s.pipelineScour({ mode: 'current', gapRatio: 1 }), 0.625 * Math.exp(-0.6), 1e-14);
  const o = s.pipelineScourOnset({ D: 0.5, e: 0.05, n: 0.4, s: 2.65, U: 1.5 });
  // Ucr^2 = g D (1 - n)(s - 1)·0.025 exp(9 sqrt(0.1)) = 9.80665·0.5·0.6·1.65·0.025·17.219
  close(o.Ucr ** 2, 9.80665 * 0.5 * 0.6 * 1.65 * 0.025 * Math.exp(9 * Math.sqrt(0.1)), 1e-12);
  close(o.Ucr, 1.446, 2e-3);
  assert.equal(o.onset, true);
});

test('Soulsby (1997) critical Shields parameter for 0.2 mm sand: D* = 4.12, theta_cr = 0.055', () => {
  const r = s.criticalShields({ d: 0.0002 });
  close(r.Dstar, 4.12, 2e-3);
  close(r.thetaCr, 0.3 / (1 + 1.2 * r.Dstar) + 0.055 * (1 - Math.exp(-0.02 * r.Dstar)), 1e-12);
  close(r.thetaCr, 0.0548, 3e-3);
  // coarse limit 0.055
  close(s.criticalShields({ d: 0.2 }).thetaCr, 0.055, 2e-3);
});

test('bed shear stress (Soulsby): current drag coefficient and wave friction factor by hand; combined mean and maximum', () => {
  const d = 0.001, h = 20;
  const z0 = 2.5 * d / 30;
  const r = s.bedShearStress({ d, h, U: 1, Uw: 0.8, T: 10 });
  const Cd = (0.4 / (Math.log(h / z0) - 1)) ** 2;
  close(r.tauC, 1025 * Cd, 1e-12);
  const fw = 1.39 * (0.8 * 10 / (2 * Math.PI) / z0) ** -0.52;
  close(r.tauW, 0.5 * 1025 * fw * 0.64, 1e-12);
  const tm = r.tauC * (1 + 1.2 * (r.tauW / (r.tauC + r.tauW)) ** 3.2);
  close(r.tauMean, tm, 1e-12);
  close(r.tauMax, tm + r.tauW, 1e-12);
  // waves at right angles
  close(s.bedShearStress({ d, h, U: 1, Uw: 0.8, T: 10, phi: Math.PI / 2 }).tauMax, Math.hypot(tm, r.tauW), 1e-12);
});

test('rock sizing by the Shields approach: the returned stone is exactly at threshold and grows with amplification', () => {
  const o = { h: 25, U: 1.2, Uw: 1.0, T: 9, amplification: 4 };
  const r = s.rockSizeShields(o);
  close(r.tauMax, r.tauCr, 1e-6);
  assert.ok(r.d50 > 0.01 && r.d50 < 1);
  assert.ok(s.rockSizeShields({ ...o, amplification: 2 }).d50 < r.d50);
});

test('Terzaghi filter rules', () => {
  const ok = s.filterCriteria({ D15f: 0.004, D50f: 0.01, D15b: 0.0005, D50b: 0.001, D85b: 0.002 });
  assert.equal(ok.ok, true);
  close(ok.retention.ratio, 2, 1e-12);
  close(ok.permeability.ratio, 8, 1e-12);
  const bad = s.filterCriteria({ D15f: 0.02, D50f: 0.05, D15b: 0.0005, D50b: 0.001, D85b: 0.002 });
  assert.equal(bad.retention.ok, false);
  assert.equal(bad.uniformity.ok, false);
});

test('armour take-off and scour effect on overburden (hand)', () => {
  const a = s.armourTakeoff({ D: 6, extentFactor: 2.5, thickness: 1.5, porosity: 0.4, rhoS: 2650 });
  close(a.area, Math.PI * (225 - 9), 1e-12);
  close(a.mass, Math.PI * 216 * 1.5 * 0.6 * 2650, 1e-12);
  const e = s.scourEffect({ embedment: 30, generalScour: 1, localScour: 3, gammaSub: 9000, reductionDepth: 12 });
  assert.equal(e.embedment, 26);
  assert.equal(e.sigmaV(3.9), 0);
  close(e.sigmaV(16), 9000 * 15, 1e-12);            // fully recovered at 4 + 12 m
  close(e.sigmaV(10), 9000 * 15 * 6 / 12, 1e-12);   // linear in between
  close(e.sigmaV(25), 9000 * 24, 1e-12);
  close(s.globalScour({ Dref: 2, kG: 0.4, localScour: 2.6 }).total, 3.4, 1e-12);
  close(s.largeStructureScour({ Lref: 30, kS: 0.05, KC: 0.5, DoverL: 0.2 }).S, 1.5, 1e-12);
});
