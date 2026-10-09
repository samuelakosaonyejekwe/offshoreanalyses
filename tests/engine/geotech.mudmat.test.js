import test from 'node:test';
import assert from 'node:assert/strict';
import * as m from '../../src/engine/geotech/mudmat.js';

const close = (a, e, rel = 1e-9, msg = '') => assert.ok(Math.abs(a - e) <= rel * Math.max(Math.abs(e), 1e-300), `${msg} got ${a}, expected ${e}`);

test('Green (1954) V-H interaction: 2 + pi at H = 0, 1 + pi/2 at full sliding; inverse round trip', () => {
  close(m.greenVH(0), 2 + Math.PI, 1e-14);
  close(m.greenVH(1), 1 + Math.PI / 2, 1e-14);
  for (const h of [0.2, 0.5, 0.9]) close(m.greenHofV(m.greenVH(h) / (2 + Math.PI)), h, 1e-9);
  assert.equal(m.greenHofV(0.3), 1);
});

test('plastic torsion of a square base: closed form vs numerical integration of r dA', () => {
  const B = 2, L = 2, su = 1;
  // mean distance from the centre of a unit quadrant = [sqrt(2) + ln(1 + sqrt(2))]/3 = 0.76520
  close(m.torsionalCapacityRect({ B, L, su }), 4 * (Math.SQRT2 + Math.log(1 + Math.SQRT2)) / 3, 1e-12);
  const B2 = 3, L2 = 7;
  let s = 0;
  const n = 600;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) s += Math.hypot(-B2 / 2 + (i + 0.5) * B2 / n, -L2 / 2 + (j + 0.5) * L2 / n);
  close(m.torsionalCapacityRect({ B: B2, L: L2, su: 5 }), 5 * s * (B2 / n) * (L2 / n), 1e-5);
});

test('torsion-equivalent horizontal load: H for T = 0, 4T/L for H = 0', () => {
  close(m.torsionEquivalentH({ H: 100, T: 0, Leff: 5 }), 100, 1e-14);
  close(m.torsionEquivalentH({ H: 0, T: 50, Leff: 5 }), 40, 1e-14);
  close(m.torsionEquivalentH({ H: 30, T: 50, Leff: 5 }), 20 + Math.hypot(30, 20), 1e-14);
});

test('classical capacity, central vertical load on clay: square mat 6.17 su A (Brinch Hansen)', () => {
  const r = m.mudmatClassicalCapacity({ B: 6, L: 6, V: 1e5, su: 5e3, gamma: 5e3 });
  close(r.Vult, (2 + Math.PI) * 1.2 * 5e3 * 36, 1e-12);
  close(r.ucBearing, 1e5 / r.Vult, 1e-12);
  assert.equal(r.ucSliding, 0);
});

test('classical capacity: moments shrink the effective area by hand (e = M/V)', () => {
  const r = m.mudmatClassicalCapacity({ B: 6, L: 8, V: 4e5, MB: 2e5, ML: 4e5, su: 8e3, gamma: 5e3 });
  close(r.eB, 0.5, 1e-12);
  close(r.eL, 1.0, 1e-12);
  close(r.Aeff, 5 * 6, 1e-12);
  // Hansen additive form on the 5 x 6 m effective area: 5.14 su (1 + 0.2·5/6)
  close(r.Vult, (2 + Math.PI) * 8e3 * (1 + 0.2 * 5 / 6) * 30, 1e-12);
  // sliding with torsion uses the equivalent load on the effective area
  const t = m.mudmatClassicalCapacity({ B: 6, L: 8, V: 4e5, HB: 3e4, T: 6e4, su: 8e3, gamma: 5e3 });
  close(t.Hequiv, 2 * 6e4 / 8 + Math.hypot(3e4, 2 * 6e4 / 8), 1e-12);
  close(t.ucSliding, t.Hequiv / (8e3 * 48), 1e-12);
  assert.throws(() => m.mudmatClassicalCapacity({ B: 6, L: 8, V: 1e5, MB: 4e5, su: 8e3, gamma: 5e3 }), /eccentricity/);
});

test('generic envelope: uniaxial points lie on the surface; v = 0.5 gives the peak no-tension moment', () => {
  const base = { Vult: 1000, Hult: 200, Mult: 500 };
  close(m.mudmatEnvelope({ ...base, V: 300, H: 200 }).f, 1, 1e-12);            // pure sliding below v = 0.5
  close(m.mudmatEnvelope({ ...base, V: 500, M: 500 }).f, 1, 1e-12);            // m* = 4 v (1 - v) = 1
  close(m.mudmatEnvelope({ ...base, V: 250, M: 375 }).f, 1, 1e-12);            // m* = 0.75
  const r = m.mudmatEnvelope({ ...base, V: 500, H: 100, M: 250 });
  close(r.f, 0.25 + 0.25, 1e-12);
  close(r.loadFactor, Math.SQRT2, 1e-9);
  assert.equal(m.mudmatEnvelope({ ...base, V: 1200 }).inside, false);
  // Green strip result reproduced through the envelope: V = 0.8 Vult
  close(m.mudmatEnvelope({ ...base, V: 800, H: 200 * m.greenHofV(0.8) }).f, 1, 1e-9);
});

test('skirt penetration in clay by hand; self-weight penetration closed form for uniform strength without tip', () => {
  const r = m.skirtPenetrationClay({ z: 0.5, wallLength: 40, t: 0.01, su0: 4e3, k: 2e3, alpha: 0.5, Nc: 7.5, gamma: 6e3 });
  close(r.friction, 0.5 * (4e3 + 0.5e3) * 0.5 * 40 * 2, 1e-12);
  close(r.tip, (7.5 * 5e3 + 3e3) * 0.4, 1e-12);
  // uniform su, t = 0: W = alpha su 2 P z  ->  z = W / (alpha su 2 P)
  const s = m.selfWeightPenetration({ W: 30e3, skirtHeight: 1, wallLength: 40, t: 0, su0: 5e3, alpha: 0.5 });
  close(s.z, 30e3 / (0.5 * 5e3 * 80), 1e-8);
  assert.equal(s.full, false);
  close(s.ballastRequired, 0.5 * 5e3 * 80 - 30e3, 1e-12);
  assert.equal(m.selfWeightPenetration({ W: 1e9, skirtHeight: 1, wallLength: 40, t: 0.01, su0: 5e3 }).full, true);
});

test('CPT skirt resistance: linear qc profile integrates exactly', () => {
  const r = m.skirtPenetrationCpt({ z: 2, qc: (z) => 1e6 * z, wallLength: 10, t: 0.02, kp: 0.4, kf: 0.03 });
  close(r.tip, 0.4 * 0.2 * 2e6, 1e-12);
  close(r.friction, 2 * 10 * 0.03 * 0.5 * 1e6 * 4, 1e-12);
});

test('contact pressure: full contact V/A (1 ± 6e/B); uniaxial lift-off 2V / (3 L (B/2 - e))', () => {
  const full = m.contactPressure({ B: 4, L: 6, V: 2400, ex: 0.5 });
  close(full.pMax, 100 * (1 + 0.75), 1e-12);
  close(full.pMin, 100 * 0.25, 1e-12);
  assert.equal(full.fullContact, true);
  const lift = m.contactPressure({ B: 4, L: 6, V: 2400, ex: 1.2, n: 400 });
  close(lift.pMax, 2 * 2400 / (3 * 6 * (2 - 1.2)), 2e-3);
  close(lift.contactFraction, 3 * 0.8 / 4, 5e-3);
  assert.equal(lift.fullContact, false);
});

test('contact pressure, bi-axial lift-off: equilibrium of the integrated pressure field', () => {
  const B = 5, L = 8, V = 1e6, ex = 1.3, ey = 2.0;
  const r = m.contactPressure({ B, L, V, ex, ey, n: 300 });
  let F = 0, Mx = 0, My = 0;
  const n = 900;
  for (let i = 0; i < n; i++) {
    const x = -B / 2 + (i + 0.5) * B / n;
    for (let j = 0; j < n; j++) {
      const y = -L / 2 + (j + 0.5) * L / n;
      const p = Math.max(r.p0 + r.gx * x + r.gy * y, 0) * (B / n) * (L / n);
      F += p; Mx += p * x; My += p * y;
    }
  }
  close(F, V, 2e-3);
  close(Mx / F, ex, 3e-3);
  close(My / F, ey, 3e-3);
  assert.ok(r.contactFraction < 0.6);
});

test('landing on a linear soil spring: energy balance and velocity-limit round trip', () => {
  const p = { mass: 50e3, addedMass: 30e3, W: 400e3, k: 2e7 };
  const r = m.landingImpact({ ...p, v: 0.5 });
  // 0.5 M v^2 + W d = 0.5 k d^2
  close(0.5 * 80e3 * 0.25 + 400e3 * r.delta, 0.5 * 2e7 * r.delta ** 2, 1e-12);
  // zero velocity: sudden release gives twice the static load
  close(m.landingImpact({ ...p, v: 0 }).daf, 2, 1e-12);
  const v = m.landingVelocityLimit({ ...p, Qallow: r.Fmax });
  close(v, 0.5, 1e-10);
  assert.equal(m.landingVelocityLimit({ ...p, Qallow: 700e3 }), 0);
});

test('breakout force and on-bottom stability by hand', () => {
  const b = m.breakoutForce({ W: 100e3, A: 20, su: 5e3, breakoutFactor: 0.5, Nc: 6, skirtFriction: 10e3 });
  close(b.F, 100e3 + 10e3 + 0.5 * 6 * 5e3 * 20, 1e-12);
  const s = m.onBottomStability({ W: 200e3, FH: 40e3, FL: 50e3, arm: 2, B: 6, mu: 0.5 });
  close(s.fsSliding, 0.5 * 150e3 / 40e3, 1e-12);
  close(s.fsOverturning, 150e3 * 3 / 80e3, 1e-12);
  close(s.fsUplift, 4, 1e-12);
});

test('required mat size solver returns a mat that just satisfies the governing check', () => {
  const loads = { V: 6e5, HB: 4e4, MB: 2e5, su: 6e3, gamma: 5e3, aspect: 1.5 };
  const r = m.requiredMatSize(loads);
  close(Math.max(r.ucBearing, r.ucSliding), 1, 1e-6);
  const smaller = m.mudmatClassicalCapacity({ ...loads, B: r.B * 0.98, L: r.L * 0.98 });
  assert.ok(Math.max(smaller.ucBearing, smaller.ucSliding) > 1);
});
