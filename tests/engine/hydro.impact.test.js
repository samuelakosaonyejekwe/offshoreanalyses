import test from 'node:test';
import assert from 'node:assert/strict';
import * as im from '../../src/engine/hydro/impact.js';
import { besselJ0, besselJ1 } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);
const DEG = Math.PI / 180;

test('Wedge water entry: von Kármán F = pi rho V³ t / tan² beta; Wagner / von Kármán = pi²/4', () => {
  const r = im.wedgeWaterEntry({ deadriseDeg: 20, V: 4, t: 0.05, rho: 1000 });
  const tb = Math.tan(20 * DEG);
  close(r.vonKarman.force, (Math.PI * 1000 * 64 * 0.05) / (tb * tb), 1e-12);
  close(r.wagner.force / r.vonKarman.force, (Math.PI * Math.PI) / 4, 1e-12);
  close(r.ratio, 2.4674, 1e-4);
  close(r.wagner.c / r.vonKarman.c, Math.PI / 2, 1e-12);
  // force = d/dt (added mass × V): finite-difference check of the Wagner added mass
  const dt = 1e-6;
  const a1 = im.wedgeWaterEntry({ deadriseDeg: 20, V: 4, t: 0.05 + dt, rho: 1000 }).wagner.addedMass;
  const a0 = im.wedgeWaterEntry({ deadriseDeg: 20, V: 4, t: 0.05 - dt, rho: 1000 }).wagner.addedMass;
  close(((a1 - a0) / (2 * dt)) * 4, r.wagner.force, 1e-6);
  // Wagner peak pressure coefficient pi²/(4 tan² beta): 79.4 at 10° deadrise
  close(im.wedgeWaterEntry({ deadriseDeg: 10, V: 1, t: 0.01 }).wagner.Cp, 79.36, 1e-3);
  assert.ok(im.wedgeWaterEntry({ deadriseDeg: 20, V: 4, t: 0.05, halfBeam: 0.5 }).chineWetted);
});

test('DNV-RP-C205 slamming pressure coefficients and Campbell–Weynberg cylinder slam coefficient', () => {
  assert.equal(im.slammingPressureCoefficient({ shape: 'cylinder' }), 5.15);
  close(im.slammingPressureCoefficient({ shape: 'flat' }), 2 * Math.PI, 1e-15);
  close(im.slammingPressureCoefficient({ shape: 'wedge', deadriseDeg: 30 }), 2.5 / Math.tan(30 * DEG) ** 1.1, 1e-15);
  close(im.slammingPressureCoefficient({ shape: 'wedge', deadriseDeg: 45 }), 2.5, 1e-12);
  assert.throws(() => im.slammingPressureCoefficient({ shape: 'wedge', deadriseDeg: 10 }), RangeError);
  close(im.cylinderSlamCoefficient({ s: 0, D: 1 }), 5.15, 1e-15);
  close(im.cylinderSlamCoefficient({ s: 1, D: 1 }), 5.15 * (1 / 20 + 0.107), 1e-12);
  close(im.cylinderSlamCoefficient({ s: 0.1, D: 2 }), 5.15 * (2 / 3.9 + 0.00535), 1e-12);
});

test('Ochi slamming probability: exp(-(d²/2m0r + v²/2m0v)); threshold 0.093 sqrt(gL) = 12 ft/s for a 520 ft ship', () => {
  const L = 520 * 0.3048;
  const r = im.slammingProbability({ m0r: 4, m0v: 3, d: 5, L });
  close(r.vThreshold, 12 * 0.3048, 5e-3);
  close(r.probability, Math.exp(-(25 / 8 + r.vThreshold ** 2 / 6)), 1e-12);
  close(r.ratePerHour, (3600 / (2 * Math.PI)) * Math.sqrt(3 / 4) * r.probability, 1e-12);
  close(im.slammingProbability({ m0r: 4, m0v: 3, d: 5, vThreshold: 0 }).probability, Math.exp(-25 / 8), 1e-12);
  const w = im.deckWetnessProbability({ m0r: 4, m2r: 3, freeboard: 6 });
  close(w.probability, Math.exp(-36 / 8), 1e-12);
});

test('Dam-break (Ritter) deck flow: h = 4/9 h0 and u = 2/3 sqrt(g h0) at the edge; front at 2 sqrt(g h0) t; volume', () => {
  const g = 9.80665; const h0 = 3; const c0 = Math.sqrt(g * h0); const t = 1.5;
  const e = im.damBreakDeckFlow({ h0, x: 0, t });
  close(e.h, (4 / 9) * h0, 1e-12); close(e.u, (2 / 3) * c0, 1e-12);
  close(e.frontPosition, 2 * c0 * t, 1e-12);
  assert.ok(!im.damBreakDeckFlow({ h0, x: 2.1 * c0 * t, t }).arrived);
  // water on the deck side: ∫ h dx from 0 to the front = (8/27) h0 c0 t; flux through the edge = h u = (8/27) h0 c0
  let vol = 0; const n = 20000; const xf = 2 * c0 * t;
  for (let i = 0; i < n; i++) vol += im.damBreakDeckFlow({ h0, x: ((i + 0.5) * xf) / n, t }).h * (xf / n);
  close(vol, (8 / 27) * h0 * c0 * t, 1e-6);
  close(e.h * e.u, (8 / 27) * h0 * c0, 1e-12);
  close(e.momentumFlux, 1025 * e.h * e.u * e.u, 1e-12); close(e.stagnationPressure, 0.5 * 1025 * e.u * e.u, 1e-12);
});

test('Sloshing: rectangular w² = g (pi/L) tanh(pi h/L); cylindrical uses the first zero of J1\' (1.8412)', () => {
  const r = im.sloshingPeriods({ shape: 'rectangular', L: 10, h: 5 });
  // hand: g·0.314159·tanh(1.5708) = 9.80665·0.314159·0.917152 = 2.82563 -> w = 1.68096 rad/s -> T = 3.7379 s
  close(r.periods[0], 3.7379, 1e-4);
  close(r.omegas[1] ** 2, 9.80665 * (2 * Math.PI / 10) * Math.tanh(Math.PI), 1e-12);
  const c = im.sloshingPeriods({ shape: 'cylindrical', R: 5, h: 50 });
  close(c.omegas[0] ** 2, (9.80665 * 1.841183781) / 5, 1e-9);
  // the tabulated constants are zeros of J1'(x) = J0(x) - J1(x)/x
  for (const x of [1.841183781, 5.331442774, 8.536316366]) near(besselJ0(x) - besselJ1(x) / x, 0, 1e-7);
  // shallow-fill limit: T = 2L / sqrt(g h)
  close(im.sloshingPeriods({ shape: 'rectangular', L: 40, h: 0.5 }).periods[0], 80 / Math.sqrt(9.80665 * 0.5), 2e-3);
  const chk = im.sloshingResonanceCheck({ sloshingPeriod: 8, motionPeriods: [6, 9, 12] });
  assert.deepEqual(chk.resonant, [false, true, false]);
});

test('Wave-in-deck silhouette force and VIM screening (hand calculations)', () => {
  const r = im.waveInDeckForce({ crestElevation: 16, deckBottom: 14, width: 30, u: 6, U: 1, Cd: 2.5, waveKinematicsFactor: 0.88, currentBlockageFactor: 0.8 });
  close(r.area, 60, 1e-12); close(r.force, 0.5 * 1025 * 2.5 * (0.88 * 6 + 0.8) ** 2 * 60, 1e-12); close(r.elevation, 15, 1e-12);
  assert.equal(im.waveInDeckForce({ crestElevation: 12, deckBottom: 14, width: 30, u: 6, Cd: 2.5 }).force, 0);
  const v = im.vimScreening({ U: 1.2, Tn: 150, D: 30, AoverD: 0.5 });
  close(v.Ur, 6, 1e-12); assert.ok(v.lockIn); close(v.amplitude, 15, 1e-12); close(v.UrResonance, 5, 1e-12);
  close(v.sheddingPeriod, 30 / (0.2 * 1.2), 1e-12);
  assert.ok(!im.vimScreening({ U: 0.4, Tn: 150, D: 30 }).lockIn);
});
