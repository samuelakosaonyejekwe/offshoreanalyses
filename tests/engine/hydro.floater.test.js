import test from 'node:test';
import assert from 'node:assert/strict';
import * as fl from '../../src/engine/hydro/floater.js';
import { rectangleAddedMassCoefficient, rectangleRollInertiaCoefficient, finiteCylinderAddedMassCoefficient } from '../../src/engine/hydro/morison.js';

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

test('Wave kinematics are integrated exactly along each strip: forces do not depend on the number of strips', () => {
  const D = 8; const draft = 40; const Aw = (Math.PI * D * D) / 4; const ma = (2 / 3) * Math.PI * rho * 64;
  const spar = (nStrips, depth) => fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 4], D, end1: 'hemisphere' }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho, g, nStrips, ...(depth ? { depth } : {}) });
  close(spar(12).A[2][2], ma, 1e-12);
  // deep water, kT = 4.5 (T = 6 s): the heave force (rho g A - omega² m_a) exp(-kT) is 1 % of rho g A; centre sampling
  // of 12 strips misses it by (k dl)²/24 of rho g A (1 - exp(-kT)), a third of the force
  for (const T of [25, 10, 6, 5]) {
    const w = (2 * Math.PI) / T; const k = w * w / g;
    const expected = (rho * g * Aw - w * w * ma) * Math.exp(-k * draft);
    for (const n of [1, 12, 50]) {
      const r = fl.floaterRAO({ model: spar(n), omega: [w], viscous: false });
      close(r.force.re[2][0], expected, 1e-9, `heave force, T = ${T} s, ${n} strips`);
      near(r.force.im[2][0], 0, 1e-9 * rho * g * Aw);
      // surge: (1 + Ca) rho A ∫ a_x dz = 2 rho g A (1 - exp(-kT)), 90° ahead of the crest
      close(r.force.im[0][0], 2 * rho * g * Aw * (1 - Math.exp(-k * draft)), 1e-12, `surge force, T = ${T} s, ${n} strips`);
    }
  }
  // finite depth 60 m: bottom pressure rho g cosh k(d - T)/cosh kd, vertical acceleration -omega² sinh k(d - T)/sinh kd;
  // horizontal ∫ cosh k(z + d)/sinh kd dz = (1 - sinh k(d - T)/sinh kd)/k
  const d = 60; const w = 0.5; let k = w * w / g;
  for (let i = 0; i < 50; i++) k = (w * w) / (g * Math.tanh(k * d));
  for (const n of [1, 12]) {
    const r = fl.floaterRAO({ model: spar(n, d), omega: [w], viscous: false });
    close(r.force.re[2][0], rho * g * Aw * Math.cosh(k * (d - draft)) / Math.cosh(k * d) - w * w * ma * Math.sinh(k * (d - draft)) / Math.sinh(k * d), 1e-9, `finite depth, ${n} strips`);
    close(r.force.im[0][0], 2 * rho * Aw * w * w * (1 - Math.sinh(k * (d - draft)) / Math.sinh(k * d)) / k, 1e-9);
  }
  // an inclined submerged brace in an oblique wave: one strip and sixty give the same force (moments use the strip-centre lever)
  const brace = (nStrips) => fl.buildFloater({ members: [{ type: 'cylinder', p1: [-10, -4, -25], p2: [12, 9, -6], D: 2 }], cog: [0, 0, -15], gyradii: [5, 5, 5], rho, g, nStrips });
  const a = fl.floaterRAO({ model: brace(1), omega: [0.9], headingDeg: 35, viscous: false, radiation: 'none' });
  const b = fl.floaterRAO({ model: brace(60), omega: [0.9], headingDeg: 35, viscous: false, radiation: 'none' });
  for (let i = 0; i < 3; i++) { close(a.force.re[i][0], b.force.re[i][0], 1e-10); close(a.force.im[i][0], b.force.im[i][0], 1e-10); }
});

test('Free end of a column: the transverse added mass loses that of a finite cylinder (DNV-RP-C205 Table D-2) at the end', () => {
  const D = 8; const draft = 40; const Aw = (Math.PI * D * D) / 4;
  const mk = (free1) => fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 4], D, free1 }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho, g });
  const strip = rho * Aw * draft;
  // the column and its image in the still-water plane: length 80 m, diameter 8 m -> C_A = 1 - 0.36/10
  close(finiteCylinderAddedMassCoefficient(10), 0.964, 1e-12);
  const loss = (1 - 0.964) * strip;
  const m0 = mk(false); const m1 = mk(true);
  close(m0.A[0][0], strip, 1e-12); close(m1.A[0][0], strip - loss, 1e-12); close(m1.A[1][1], strip - loss, 1e-12);
  close(m1.A[2][2], 0, 1e-12, 'no axial change');
  // the loss sits at the keel: lever 40 m in the coupling and the pitch inertia
  close(m1.A[4][4], m0.A[4][4] - loss * draft * draft, 1e-12); close(m1.A[0][4], m0.A[0][4] + loss * draft, 1e-12);
  assert.ok(m0.A[4][4] / m1.A[4][4] > 1.1);
  // and in the excitation: long-wave surge force omega² (rho V + A11), with the keel term at exp(-kT)
  const w = 0.3; const k = w * w / g;
  const f = (m) => fl.floaterRAO({ model: m, omega: [w], viscous: false }).force.im[0][0];
  close(f(m0) - f(m1), w * w * loss * Math.exp(-k * draft), 1e-9);
  // a dry end is not corrected
  close(fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 4], D, free2: true }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho, g }).A[0][0], strip, 1e-12);
});

test('Radiation damping from the Haskind relation: closed forms for a vertical column, symmetry, resonance amplitude', () => {
  const D = 8; const draft = 40; const Aw = (Math.PI * D * D) / 4;
  const mk = (depth) => fl.buildFloater({ members: [{ type: 'cylinder', p1: [0, 0, -draft], p2: [0, 0, 4], D, end1: 'hemisphere' }], cog: [0, 0, -25], gyradii: [12, 12, 3], rho, g, ...(depth ? { depth } : {}) });
  const model = mk();
  for (const w of [0.3, 0.48, 0.8]) {
    const r = fl.floaterRAO({ model, omega: [w], headingDeg: 0, viscous: false });
    const X = (i) => ({ re: r.force.re[i][0], im: r.force.im[i][0] });
    const Bm = fl.radiationDamping({ model, omega: [w] })[0];
    const c = w ** 3 / (rho * g ** 3);
    // heave force independent of heading: ∫|X3|² = 2 pi |X3|²; surge force ∝ cos(beta): ∫ = pi |X1|²
    close(Bm[2][2], c * (X(2).re ** 2 + X(2).im ** 2) / 2, 1e-10, 'B33');
    close(Bm[0][0], c * (X(0).re ** 2 + X(0).im ** 2) / 4, 1e-10, 'B11');
    close(Bm[1][1], Bm[0][0], 1e-10, 'B22'); close(Bm[3][3], Bm[4][4], 1e-10, 'B44');
    close(Bm[0][4], c * (X(0).re * X(4).re + X(0).im * X(4).im) / 4, 1e-10, 'B15');
    close(Bm[1][3], -Bm[0][4], 1e-10, 'B24 = -B15 for a body of revolution');
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) near(Bm[i][j], Bm[j][i], 1e-9 * Bm[4][4]);
    near(Bm[0][2], 0, 1e-9 * Math.sqrt(Bm[0][0] * Bm[2][2])); near(Bm[5][5], 0, 1e-12 * Bm[4][4]);
    assert.deepEqual(r.radiationDamping[0], Bm);
    // the heading count does not matter once the phase variation is resolved
    close(fl.radiationDamping({ model, omega: [w], nHeadings: 64 })[0][0][0], Bm[0][0], 1e-10);
  }
  // long waves: X3 -> rho g A_wp, so B33 -> rho A_wp² omega³ / (2 g), the damping of a small heaving body
  const wl = 0.02;
  close(fl.radiationDamping({ model, omega: [wl] })[0][2][2], rho * Aw * Aw * wl ** 3 / (2 * g), 1e-2);
  // finite depth: B33 = k |X3|² / (4 rho g c_g), c_g = (omega/2k)(1 + 2kd/sinh 2kd)
  const d = 60; const w = 0.5; let k = w * w / g;
  for (let i = 0; i < 50; i++) k = (w * w) / (g * Math.tanh(k * d));
  const shallow = mk(d);
  const rs = fl.floaterRAO({ model: shallow, omega: [w], viscous: false });
  const cg = (w / (2 * k)) * (1 + 2 * k * d / Math.sinh(2 * k * d));
  close(rs.radiationDamping[0][2][2], k * (rs.force.re[2][0] ** 2 + rs.force.im[2][0] ** 2) / (4 * rho * g * cg), 1e-10);
  // heave resonance with radiation damping alone: |z| = |X3| / (omega B33) = 2 rho g³ / (omega⁴ |X3|)
  const wn = (2 * Math.PI) / fl.naturalPeriods(model).uncoupled[2];
  const res = fl.floaterRAO({ model, omega: [wn], viscous: false });
  close(res.amp[2][0], 2 * rho * g ** 3 / (wn ** 4 * Math.abs(res.force.re[2][0])), 1e-6);
  // 'none' leaves the resonance undamped; an added fraction of critical gives |z| = |X3| / (2 zeta C33) without radiation
  const zeta = 0.03;
  const rz = fl.floaterRAO({ model, omega: [wn], viscous: false, radiation: 'none', extraDampingRatio: zeta });
  assert.equal(rz.radiationDamping, null);
  close(rz.amp[2][0], Math.abs(rz.force.re[2][0]) / (2 * zeta * model.C[2][2]), 1e-6);
  // the former name of the added ratio is still accepted, and both together are refused
  const old = fl.floaterRAO({ model, omega: [wn], viscous: false, radiation: 'none', radiationDampingRatio: zeta });
  close(old.amp[2][0], rz.amp[2][0], 1e-12);
  assert.throws(() => fl.floaterRAO({ model, omega: [wn], extraDampingRatio: 0.01, radiationDampingRatio: 0.01 }), /not both/);
  assert.throws(() => fl.floaterRAO({ model, omega: [wn], radiation: 'ratio' }), /haskind/);
});

test('Radiation damping of a four-column unit: symmetric, positive semi-definite, converged in the number of headings', () => {
  const model = semi();
  const omega = [0.2, 0.35, 0.6, 0.9];
  const Bm = fl.radiationDamping({ model, omega });
  Bm.forEach((M6, iw) => {
    for (let i = 0; i < 6; i++) {
      assert.ok(M6[i][i] >= 0, `B${i + 1}${i + 1} at ${omega[iw]} rad/s`);
      for (let j = 0; j < 6; j++) { near(M6[i][j], M6[j][i], 1e-9 * Math.sqrt(M6[i][i] * M6[j][j]) + 1e-6); assert.ok(M6[i][j] ** 2 <= M6[i][i] * M6[j][j] * (1 + 1e-9) + 1e-6); }
    }
  });
  // doubling the heading count changes nothing: the default resolves the phase differences between the columns
  const fine = fl.radiationDamping({ model, omega, nHeadings: 720 });
  for (let iw = 0; iw < omega.length; iw++) for (const i of [0, 2, 4, 5]) close(Bm[iw][i][i], fine[iw][i][i], 1e-8, `B${i + 1}${i + 1} at ${omega[iw]} rad/s`);
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
  // beam sea: no pitch; the bottom pressure keeps its depth decay exp(-kT) and is averaged over the beam instead of the length
  const beam = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90 });
  near(beam.pitch[0], 0, 1e-12); close(beam.heave[0], eta * f * Math.exp(-k * T) * Math.sin(k * B / 2) / (k * B / 2), 1e-9);
  // the limits of the stated validity follow the heave natural period of the same equation of motion
  close(beam.validPeriodMin.heave, 2 * beam.heaveNaturalPeriod, 1e-12); assert.equal(beam.validPeriodMin.roll, null);
});

test('Box barge heave force is the bottom-pressure integral at any heading: exp(-kT) sinc(kx L/2) sinc(ky B/2)', () => {
  const L = 100; const B = 20; const T = 5;
  // long beam waves: heave -> exp(-kT)(1 - kT)/(1 - 2kT) = 1 + O((kT)²), not 1 + kT (kT = 0.08: 1.011, not 1.08)
  const w = Math.sqrt(0.08 * g / T); const k = w * w / g;
  const beam = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90 }).heave[0];
  const A = 2 * Math.sin(0.5 * k * B) * Math.exp(-k * T); const d = A * A / (k * B);
  close(beam, Math.exp(-k * T) * (Math.sin(k * B / 2) / (k * B / 2)) * Math.hypot(1 - k * T, d) / Math.hypot(1 - 2 * k * T, d), 1e-12);
  assert.ok(Math.abs(beam - 1) < 0.02, `beam-sea heave at kT = 0.08: ${beam}`);
  // a free body follows a very long wave from any direction
  for (const h of [0, 45, 90, 135, 180]) close(fl.boxBargeRAO({ L, B, T, omega: [0.03], headingDeg: h }).heave[0], 1, 1e-3, `heading ${h}`);
  // oblique sea: the ratio to head sea at the same frequency is the ratio of the two pressure integrals
  const w2 = 0.5; const k2 = w2 * w2 / g; const bt = 3 * Math.PI / 4;
  const sinc = (x) => Math.sin(x) / x;
  const obl = fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 135 }); const head = fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 180 });
  close(obl.heave[0] / head.heave[0], Math.abs(sinc(0.5 * k2 * L * Math.cos(bt)) * sinc(0.5 * k2 * B * Math.sin(bt)) / sinc(0.5 * k2 * L)), 1e-9);
  // heading symmetry of a box
  close(fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 45 }).heave[0], obl.heave[0], 1e-12);
  close(fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 45 }).pitch[0], obl.pitch[0], 1e-12);
});

test('Box barge roll: quasi-static roll equals the wave slope in long beam waves; added inertia, natural period, damping', () => {
  const L = 100; const B = 20; const T = 5; const GMt = 3;
  const disp = rho * L * B * T; const C44 = disp * g * GMt;
  // strip added inertia of the section 20 × 5 m: half the double body 20 × 10 m, DNV-RP-C205 Table D-1, ratio 2 -> 0.15
  const A44 = 0.5 * 0.15 * rho * Math.PI * 10 ** 4 * L;
  const w = 0.05;
  const r = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90, GMt, kxx: 8 });
  close(r.rollAddedInertia, A44, 1e-12);
  const wn0 = Math.sqrt(g * GMt) / 8;
  // wave slope k times the 1-DOF dynamic amplification (0.5 %), less the added-inertia moment of the relative roll
  // (0.2 %) and the first-order depth decay of the Froude–Krylov moment about G (0.19 %). The latter by hand, from the
  // pressure on the bottom (e^{-kT}) and on the two sides (lever z) and the horizontal force at the lever z_G:
  // M_FK = rho g L k [B³/12 (1 - kT) - B T²/2 (1 - 2kT/3) - z_G B T (1 - kT/2)], whose leading term is C44 k.
  const k = w * w / g; const zG = T / 2 + B * B / (12 * T) - GMt - T;
  const decay = k * T * (B ** 3 / 12 - B * T * T / 3 - zG * B * T / 2) / (B ** 3 / 12 - B * T * T / 2 - zG * B * T);
  close(B ** 3 / 12 - B * T * T / 2 - zG * B * T, B * T * GMt, 1e-12);
  // damping: 5 % of critical plus the (here negligible) radiation damping; what remains is of order (kT)²
  const amp = 1 / Math.hypot(1 - (w / wn0) ** 2, 2 * 0.05 * w / wn0 + w * r.rollRadiationDamping[0] / C44);
  close(r.roll[0], k * (1 - decay - w * w * A44 / C44) * amp, 1e-5);
  close(r.rollNaturalPeriod, 2 * Math.PI * 8 / Math.sqrt(g * GMt), 1e-12);
  close(r.validPeriodMin.roll, 1.3 * r.rollNaturalPeriod, 1e-12);
  near(fl.boxBargeRAO({ L, B, T, omega: [0.5], headingDeg: 180, GMt }).roll[0], 0, 1e-12);
  // gyradius of the mass alone: the strip added inertia is added, I44 = Δ kxxDry² + A44
  const dry = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90, GMt, kxxDry: 7 });
  close(dry.rollNaturalPeriod, 2 * Math.PI * Math.sqrt((disp * 49 + A44) / C44), 1e-12);
  assert.throws(() => fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90, GMt, kxx: 8, kxxDry: 7 }), /not both/);
  // a deep narrow section uses the draught as the longer half-side: 4 × 10 m -> double body 4 × 20 m, ratio 5
  close(fl.boxBargeRAO({ L, B: 4, T: 10, omega: [w], headingDeg: 90, GMt: 0.5 }).rollAddedInertia, 0.5 * rectangleRollInertiaCoefficient(5) * rho * Math.PI * 10 ** 4 * L, 1e-12);
  // without radiation damping the amplification at resonance is 1/(2 zeta) of the quasi-static response to the same moment
  const wn = (2 * Math.PI) / r.rollNaturalPeriod;
  const at = (z, rad) => fl.boxBargeRAO({ L, B, T, omega: [wn], headingDeg: 90, GMt, kxx: 8, rollDampingRatio: z, rollRadiation: rad });
  close(at(0.05, false).roll[0] / at(0.10, false).roll[0], 2, 1e-9);
  assert.equal(at(0.05, false).rollRadiationDamping[0], 0);
  // with it the damping is B_rad + 2 zeta sqrt(I44 C44): the two responses are in the inverse ratio of the totals
  const Brad = at(0.05, true).rollRadiationDamping[0]; const Bc = 2 * Math.sqrt(disp * 64 * C44);
  assert.ok(Brad > 0);
  close(at(0.05, true).roll[0] / at(0.10, true).roll[0], (Brad + 0.10 * Bc) / (Brad + 0.05 * Bc), 1e-9);
  close(at(0.05, true).roll[0] / at(0.05, false).roll[0], 0.05 * Bc / (Brad + 0.05 * Bc), 1e-9);
});

test('Box barge roll radiation damping satisfies the Haskind relation with the roll moment of the model', () => {
  const L = 100; const B = 20; const T = 5; const GMt = 3;
  const disp = rho * L * B * T; const C44 = disp * g * GMt;
  const A44 = 0.5 * 0.15 * rho * Math.PI * 10 ** 4 * L;
  // long waves: M(beta) -> C44 k sin(beta) (1 - omega² A44/C44), so ∫|M|² d beta = pi (…)² and
  // B_rad -> omega³ (C44 k)² (1 - omega² A44/C44)² / (4 rho g³)
  const w = 0.04; const k = w * w / g;
  const r = fl.boxBargeRAO({ L, B, T, omega: [w], headingDeg: 90, GMt, kxx: 8 });
  close(r.rollRadiationDamping[0], w ** 3 * (C44 * k * (1 - w * w * A44 / C44)) ** 2 / (4 * rho * g ** 3), 5e-3);
  // any frequency: the moment at a heading is recovered from the undamped off-resonance response, M = phi (C44 - omega² I44),
  // and its square integrated over all headings gives the damping, omega³/(4 pi rho g³) ∫|M|² d beta
  const w2 = 0.7; const I44 = disp * 64; const n = 720; let sum = 0;
  for (let j = 0; j < n; j++) {
    const phi = fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 360 * (j + 0.5) / n, GMt, kxx: 8, rollDampingRatio: 0, rollRadiation: false }).roll[0];
    sum += (phi * (C44 - w2 * w2 * I44)) ** 2;
  }
  const expected = w2 ** 3 / (4 * Math.PI * rho * g ** 3) * sum * 2 * Math.PI / n;
  close(fl.boxBargeRAO({ L, B, T, omega: [w2], headingDeg: 90, GMt, kxx: 8 }).rollRadiationDamping[0], expected, 1e-9);
  // resonance with radiation damping alone: |phi| = |M| / (omega_n B_rad)
  const wn = Math.sqrt(C44 / I44);
  const und = fl.boxBargeRAO({ L, B, T, omega: [wn * 0.5], headingDeg: 90, GMt, kxx: 8, rollDampingRatio: 0, rollRadiation: false });
  assert.ok(und.roll[0] > 0);
  const res = fl.boxBargeRAO({ L, B, T, omega: [wn], headingDeg: 90, GMt, kxx: 8, rollDampingRatio: 0 });
  const off = fl.boxBargeRAO({ L, B, T, omega: [wn], headingDeg: 90, GMt, kxx: 9, rollDampingRatio: 0, rollRadiation: false });
  const M = off.roll[0] * Math.abs(C44 - wn * wn * disp * 81); // same moment, read from a detuned undamped oscillator
  close(res.roll[0], M / (wn * res.rollRadiationDamping[0]), 1e-6);
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
