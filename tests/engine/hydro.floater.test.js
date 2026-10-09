import test from 'node:test';
import assert from 'node:assert/strict';
import * as fl from '../../src/engine/hydro/floater.js';
import { rectangleAddedMassCoefficient } from '../../src/engine/hydro/morison.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);
const rho = 1025; const g = 9.80665;

function semi(extra = {}) {
  const members = [];
  for (const y of [30, -30]) members.push({ type: 'rect', p1: [-45, y, -20], p2: [45, y, -20], b: 14, h: 8, Cd: 2 });
  for (const x of [30, -30]) for (const y of [30, -30]) members.push({ type: 'cylinder', p1: [x, y, -16], p2: [x, y, 15], D: 12, Cd: 1 });
  return fl.buildFloater({ members, cog: [0, 0, -8], gyradii: [32, 32, 38], rho, g, nStrips: 24, ...extra });
}

test('Spar (vertical cylinder) heave natural period T = 2 pi sqrt((m + a)/(rho g Aw)), a = half-disc 4/3 rho r³', () => {
  const D = 10; const draft = 60;
  const model = fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 10], D, end1: 'half-disc' }], cog: [0, 0, -40], gyradii: [20, 20, 4], rho, g });
  const Aw = (Math.PI * D * D) / 4; const m = rho * Aw * draft; const a = (4 / 3) * rho * (D / 2) ** 3;
  close(model.volume, Aw * draft, 1e-12); close(model.mass, m, 1e-12);
  close(model.A[2][2], a, 1e-12); close(model.C[2][2], rho * g * Aw, 1e-12);
  const np = fl.naturalPeriods(model);
  close(np.uncoupled[2], 2 * Math.PI * Math.sqrt((m + a) / (rho * g * Aw)), 1e-12);
  // heave is uncoupled for an axisymmetric spar: the coupled eigen-solution contains the same period
  assert.ok(np.coupled.some((T) => Math.abs(T - np.uncoupled[2]) < 1e-6 * T));
  // pitch stiffness = rho g V GM with GM = KB + BM - KG (hand)
  const GM = -30 + (D * D) / (16 * draft) + 40;
  close(model.C[4][4], rho * g * Aw * draft * GM, 1e-9);
  // surge added mass = rho pi r² × draft (Ca = 1)
  close(model.A[0][0], rho * Aw * draft, 1e-12);
});

test('Spar heave excitation: bottom pressure rho g Aw exp(-kT) corrected by the base added mass (hand formula)', () => {
  const D = 10; const draft = 60;
  const model = fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 10], D, end1: 'half-disc' }], cog: [0, 0, -40], gyradii: [20, 20, 4], rho, g, nStrips: 200 });
  const w = 0.4; const k = w * w / g; const Aw = (Math.PI * D * D) / 4;
  const r = fl.floaterRAO({ model, omega: [w], headingDeg: 0, viscous: false });
  const expected = rho * g * Aw * Math.exp(-k * draft) - (4 / 3) * rho * 125 * w * w * Math.exp(-k * draft);
  close(r.force.re[2][0], expected, 2e-5);
  near(r.force.im[2][0], 0, 1e-6 * expected);
  // surge force on a deep slender column: (1 + Ca) rho A ∫ a_x dz = 2 rho A g (1 - exp(-kT)), 90° ahead of the crest
  close(r.force.im[0][0], 2 * rho * Aw * g * (1 - Math.exp(-k * draft)), 2e-5);
});

test('Semi-submersible: heave cancellation period exists and matches the pontoon/column force balance', () => {
  const model = semi();
  const omega = Array.from({ length: 400 }, (_, i) => 0.15 + (i * 0.35) / 399);
  const r = fl.floaterRAO({ model, omega, headingDeg: 90, viscous: false });
  // beam sea, symmetric hull: heave force is real; locate its long-period sign change
  let wc = null;
  for (let i = 1; i < omega.length; i++) {
    if (r.force.re[2][i - 1] > 0 && r.force.re[2][i] <= 0) { wc = omega[i - 1] + ((omega[i] - omega[i - 1]) * r.force.re[2][i - 1]) / (r.force.re[2][i - 1] - r.force.re[2][i]); break; }
  }
  assert.ok(wc !== null, 'no cancellation found');
  // hand balance (Faltinsen 1990, ch. 3): 4 rho g A_c exp(-k z_c) = omega² · 2 (rho V_p + a_33) exp(-k z_p),
  // column base at 16 m, pontoon centre at 20 m, omega² = g k  ->  4 A_c exp(4k) = 2 (V_p + a_33/rho) k
  const Ac = Math.PI * 36; const Vp = 90 * 14 * 8;
  const a33 = rectangleAddedMassCoefficient(14 / 8) * Math.PI * 49 * 90;
  let lo = 1e-4; let hi = 0.02;
  for (let i = 0; i < 80; i++) { const k = 0.5 * (lo + hi); if (2 * (Vp + a33) * k - 4 * Ac * Math.exp(4 * k) < 0) lo = k; else hi = k; }
  const wHand = Math.sqrt(g * lo);
  close(wc, wHand, 3e-3, 'cancellation frequency');
  const Tc = (2 * Math.PI) / wc;
  assert.ok(Tc > 18 && Tc < 30, `cancellation period ${Tc}`);
  // the heave RAO has a minimum close to zero there
  const rc = fl.floaterRAO({ model, omega: [wc], headingDeg: 90, viscous: false });
  assert.ok(rc.amp[2][0] < 0.01);
});

test('Semi-submersible: long-wave limits (heave -> 1, surge -> 1, pitch -> wave slope k) and natural periods', () => {
  const model = semi();
  const np = fl.naturalPeriods(model);
  assert.ok(np.uncoupled[2] > 20 && np.uncoupled[2] < 30, `heave period ${np.uncoupled[2]}`);
  assert.equal(np.uncoupled[0], null);
  const w = (2 * Math.PI) / 300;
  const r = fl.floaterRAO({ model, omega: [w], headingDeg: 0, viscous: false, radiationDampingRatio: 0 });
  close(r.amp[2][0], 1, 0.02, 'heave'); close(r.amp[0][0], 1, 0.03, 'surge');
  close(r.amp[4][0], (w * w) / g, 0.05, 'pitch');
});

test('Stochastic drag linearisation converges and reduces the resonant heave response', () => {
  const model = semi();
  const Tn = fl.naturalPeriods(model).uncoupled[2];
  const omega = Array.from({ length: 160 }, (_, i) => 0.1 + (i * 1.4) / 159);
  omega.push((2 * Math.PI) / Tn); omega.sort((a, b) => a - b);
  const i0 = omega.indexOf((2 * Math.PI) / Tn);
  const sea = { Hs: 8, Tp: 16, gamma: 2 };
  const a = fl.floaterRAO({ model, omega, headingDeg: 90, viscous: false, radiationDampingRatio: 0.005 });
  const b = fl.floaterRAO({ model, omega, headingDeg: 90, seaState: sea, radiationDampingRatio: 0.005 });
  assert.ok(b.converged && b.iterations >= 1);
  assert.ok(b.amp[2][i0] < 0.6 * a.amp[2][i0], `${b.amp[2][i0]} vs ${a.amp[2][i0]}`);
  assert.ok(b.viscousDamping[2][2] > 0 && b.viscousDamping[1][1] > 0);
  // drag damping grows with sea state (sigma of the relative velocity)
  const c = fl.floaterRAO({ model, omega, headingDeg: 90, seaState: { Hs: 2, Tp: 16, gamma: 2 }, radiationDampingRatio: 0.005 });
  assert.ok(c.viscousDamping[2][2] < b.viscousDamping[2][2]);
});

test('Point motion: rigid-body kinematics s = xi + theta × r; acceleration = omega² × displacement', () => {
  const z = () => [0]; const rao = { omega: [0.5], re: [z(), z(), z(), [0.05], [0.1], [0.02]], im: [z(), z(), z(), z(), z(), z()] };
  const p = fl.pointMotion({ rao, point: [10, 4, 3] });
  near(p.re[2][0], 0.05 * 4 - 0.1 * 10, 1e-12);
  near(p.re[0][0], 0.1 * 3 - 0.02 * 4, 1e-12);
  near(p.re[1][0], 0.02 * 10 - 0.05 * 3, 1e-12);
  close(p.accAmp[2][0], 0.25 * 0.8, 1e-12);
});

test('Response statistics: unit RAO returns the wave statistics (significant amplitude Hs/2, MPM = sqrt(2 m0 ln N))', () => {
  const omega = Array.from({ length: 3000 }, (_, i) => 0.05 + (i * 4.95) / 2999);
  const S = fl.waveSpectrum(omega, { Hs: 6, Tp: 12, gamma: 1 });
  const st = fl.responseStatistics({ omega, amp: omega.map(() => 1), S, duration: 10800 });
  close(st.significant, 3, 2e-3);
  close(st.mpm, Math.sqrt(2 * st.m0 * Math.log(10800 / st.Tz)), 1e-12);
  // Pierson–Moskowitz closed-form moments: m2 = m0 wp² sqrt(5 pi)/2, so Tz/Tp = 0.7104; tail truncated at 5 rad/s
  close(st.Tz / 12, 0.7104, 0.02);
});

test('Air gap of a fixed deck point: relative motion = wave, sigma = Hs/4', () => {
  const stiff = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? 1e16 : 0)));
  const model = semi({ Kext: stiff });
  const omega = Array.from({ length: 300 }, (_, i) => 0.2 + (i * 1.3) / 299);
  const rao = fl.floaterRAO({ model, omega, headingDeg: 30, viscous: false });
  const S = fl.waveSpectrum(omega, { Hs: 10, Tp: 14 });
  const ag = fl.airGap({ rao, point: [20, 10, 15], staticAirGap: 15, S });
  const sig = fl.responseStatistics({ omega, amp: omega.map(() => 1), S }).sigma;
  close(ag.sigma, sig, 1e-6);
  close(ag.minAirGap, 15 - ag.mpmRelative, 1e-12);
  const amp = fl.airGap({ rao, point: [20, 10, 15], staticAirGap: 15, S, waveFactor: 1.2 });
  close(amp.sigma, 1.2 * sig, 1e-6);
});

test('Box barge, Jensen–Mansour–Olsen (2004): long-wave limits and an independently evaluated point', () => {
  const L = 100; const B = 20; const T = 5;
  const lw = fl.boxBargeRAO({ L, B, T, omega: [0.05], headingDeg: 180, GMt: 4 });
  close(lw.heave[0], 1, 1e-3); close(lw.pitch[0], 0.05 * 0.05 / g, 2e-3);
  close(lw.heaveNaturalPeriod, 2 * Math.PI * Math.sqrt(10 / g), 1e-12);
  // hand evaluation at omega = 0.6 rad/s, head sea, zero speed (alpha = 1)
  const w = 0.6; const k = w * w / g; const ke = k;
  const A = 2 * Math.sin(0.5 * k * B) * Math.exp(-k * T);
  const f = Math.sqrt((1 - k * T) ** 2 + (A * A / (k * B)) ** 2);
  const kappa = Math.exp(-ke * T);
  const F = kappa * f * (2 / (ke * L)) * Math.sin(ke * L / 2);
  const Gp = kappa * f * (24 / ((ke * L) ** 2 * L)) * (Math.sin(ke * L / 2) - (ke * L / 2) * Math.cos(ke * L / 2));
  const eta = 1 / Math.sqrt((1 - 2 * k * T) ** 2 + (A * A / (k * B)) ** 2);
  const r = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 180 });
  close(r.heave[0], eta * F, 1e-12); close(r.pitch[0], eta * Gp, 1e-12);
  assert.equal(r.roll, null);
  // beam sea: no pitch; heave forcing without the length (Smith) reduction
  const beam = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90 });
  near(beam.pitch[0], 0, 1e-12); close(beam.heave[0], eta * f, 1e-9);
});

test('Box barge roll (1-DOF Froude–Krylov): quasi-static roll equals the wave slope in long beam waves; resonance at Tn', () => {
  const L = 100; const B = 20; const T = 5; const GMt = 3;
  const w = 0.05;
  const r = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90, GMt, kxx: 8 });
  const wn0 = Math.sqrt(g * GMt) / 8;
  // wave slope k times the 1-DOF dynamic amplification (0.5 %); remaining terms are of order kT, kB
  close(r.roll[0], (w * w / g) / (1 - (w / wn0) ** 2), 3e-3);
  close(r.rollNaturalPeriod, 2 * Math.PI * 8 / Math.sqrt(g * GMt), 1e-12);
  near(fl.boxBargeRAO({ L, B, T, omega: [0.5], headingDeg: 180, GMt }).roll[0], 0, 1e-12);
  // dynamic amplification at resonance = 1/(2 zeta) of the quasi-static response to the same moment
  const wn = (2 * Math.PI) / r.rollNaturalPeriod;
  const z5 = fl.boxBargeRAO({ L, B, T, omega: [wn], headingDeg: 90, GMt, kxx: 8, rollDampingRatio: 0.05 }).roll[0];
  const z10 = fl.boxBargeRAO({ L, B, T, omega: [wn], headingDeg: 90, GMt, kxx: 8, rollDampingRatio: 0.10 }).roll[0];
  close(z5 / z10, 2, 1e-9);
});

test('Operability and limiting Hs: linear scaling with Hs', () => {
  const omega = Array.from({ length: 400 }, (_, i) => 0.1 + (i * 2.9) / 399);
  const amp = omega.map(() => 0.5);
  const lim = fl.limitingHs({ omega, amp, limit: 1.0, Tp: [8, 12] });
  // significant amplitude = 0.5 × Hs/2 -> Hs_lim = 4 m (JONSWAP m0 normalisation within 1 %)
  close(lim.HsLimit[0], 4, 0.015); close(lim.HsLimit[1], 4, 0.015);
  const op = fl.operability({
    seaStates: [{ Hs: 1, Tp: 8, prob: 0.5 }, { Hs: 3, Tp: 8, prob: 0.3 }, { Hs: 5, Tp: 10, prob: 0.2 }],
    evaluate: (Hs) => ({ heave: 0.25 * Hs, roll: 0.5 * Hs }), limits: { heave: 1.0, roll: 1.8 },
  });
  close(op.operability, 0.8, 1e-12); assert.equal(op.perState[2].governing, 'roll');
  close(op.governingCounts.roll, 0.2, 1e-12);
});
