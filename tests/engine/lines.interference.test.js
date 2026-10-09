import test from 'node:test';
import assert from 'node:assert/strict';
import {
  currentDeflection, huseWake, blevinsWakeDrag, clearanceUnderCurrent, impactVelocityEstimate,
} from '../../src/engine/lines/interference.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const rho = 1025;

test('constant tension, uniform current: parabola with y_mid = f L^2 / (8 T)', () => {
  const L = 300, T = 2e6, d = 0.5, Cd = 1.0, U = 1.0;
  const f = 0.5 * rho * Cd * d * U * U; // 256.25 N/m
  const r = currentDeflection({ L, w: 0, d, Cd, current: U, topTension: T, n: 100 });
  assert.ok(rel(r.yMax, f * L * L / (8 * T)) < 1e-12); // 1.4414 m
  assert.ok(rel(r.yMax, 1.44140625) < 1e-9);
  assert.ok(rel(r.totalDrag, f * L) < 1e-12);
  assert.ok(Math.abs(r.y[0]) < 1e-12 && Math.abs(r.y[100]) < 1e-12);
  assert.ok(rel(Math.tan(r.bottomAngle), f * L / (2 * T)) < 1e-12);
});

test('top-tensioned line with weight: deflection matches the closed-form logarithmic solution', () => {
  // Te = Tb + w s, uniform f: y(s) = -f s/w + (c + f Tb/w) ln(1 + w s/Tb)/w, c = f L/ln(1 + w L/Tb) - f Tb/w
  const L = 1000, w = 800, Ttop = 1.5e6, d = 0.4, Cd = 1.2, U = 0.8;
  const f = 0.5 * rho * Cd * d * U * U, Tb = Ttop - w * L;
  const c = f * L / Math.log(1 + w * L / Tb) - f * Tb / w;
  const exact = (s) => -f * s / w + (c + f * Tb / w) * Math.log(1 + w * s / Tb) / w;
  const r = currentDeflection({ L, w, d, Cd, current: U, topTension: Ttop, n: 400 });
  for (const i of [50, 150, 200, 320]) assert.ok(rel(r.y[i], exact(r.s[i])) < 2e-5);
  // the maximum sits below mid-height because the tension is lower there
  assert.ok(r.sAtMax < L / 2);
  // refinement: trapezoidal rule converges at second order
  const e = (n) => { const q = currentDeflection({ L, w, d, Cd, current: U, topTension: Ttop, n }); return Math.abs(q.y[n / 2] - exact(L / 2)); };
  const p = Math.log2(e(50) / e(100));
  assert.ok(Math.abs(p - 2) < 0.1);
});

test('free-hanging line: uniform line hangs straight at tan(a) = f/w; with an end weight the closed form holds', () => {
  const L = 200, w = 300, d = 0.2, Cd = 1.2, U = 0.6;
  const f = 0.5 * rho * Cd * d * U * U;
  const r = currentDeflection({ L, w, d, Cd, current: U, mode: 'free-hanging' });
  assert.ok(rel(r.y[0], f * L / w) < 1e-10);
  assert.ok(rel(Math.tan(-r.topAngle), f / w) < 1e-10 && rel(Math.tan(-r.bottomAngle), f / w) < 1e-10);
  // end weight W: y(0) = f [L/w - (W/w^2) ln(1 + w L/W)]
  const W = 2e4;
  const r2 = currentDeflection({ L, w, d, Cd, current: U, mode: 'free-hanging', endWeight: W, n: 800 });
  assert.ok(rel(r2.y[0], f * (L / w - (W / (w * w)) * Math.log(1 + w * L / W))) < 1e-5);
  assert.ok(r2.y[0] < r.y[0]);
  // sheared current profile: zero current -> zero deflection below the sheared layer is still carried by the line
  const r3 = currentDeflection({ L, w, d, Cd, current: { s: [0, 100, 200], speed: [0, 0, 1] }, mode: 'free-hanging', endWeight: W });
  assert.ok(Math.abs(r3.y[0] - r3.y[100]) < 1e-12); // no load below 100 m: that part hangs vertically
  assert.throws(() => currentDeflection({ L, w, d, Cd, current: U, topTension: 100 }), /topTension/);
});

test('Huse wake: hand calculation and half-width property', () => {
  // Cd = 1, D = 1, x = 6: xs = 10, b = 0.25 sqrt(10) = 0.7906, u/U0 = sqrt(0.1) = 0.3162
  const a = huseWake({ x: 6, D: 1, Cd: 1, U0: 2 });
  assert.ok(rel(a.xs, 10) < 1e-14 && rel(a.halfWidth, 0.25 * Math.sqrt(10)) < 1e-14);
  assert.ok(rel(a.deficit, 2 * Math.sqrt(0.1)) < 1e-14);
  assert.ok(rel(a.dragFactor, (1 - Math.sqrt(0.1)) ** 2) < 1e-14); // 0.4675
  // at y = b the deficit is half the centreline value (exp(-0.693) = 0.5)
  assert.ok(rel(huseWake({ x: 6, y: a.halfWidth, D: 1, Cd: 1, U0: 2 }).deficit, a.deficit / 2) < 3e-4);
  // deficit decays downstream, shielding vanishes far away
  assert.ok(huseWake({ x: 60, D: 1, Cd: 1 }).deficit < huseWake({ x: 6, D: 1, Cd: 1 }).deficit);
  assert.ok(huseWake({ x: 1e6, D: 1, Cd: 1 }).dragFactor > 0.99);
  assert.throws(() => huseWake({ x: -1, D: 1, Cd: 1 }), /x/);
});

test('Blevins wake drag: Cd = Cd0 [1 - sqrt(Cd D/x) exp(-4.5 y^2/(Cd D x))]^2 (hand calculation)', () => {
  // x = 10 D, Cd = 1: r = sqrt(0.1) -> Cd = (1 - 0.31623)^2 = 0.46754
  assert.ok(rel(blevinsWakeDrag({ x: 10, D: 1, Cd: 1 }).Cd, 0.467544) < 1e-5);
  // off the centreline by y = 1: r = 0.31623 exp(-0.45) = 0.20164 -> 0.63738
  assert.ok(rel(blevinsWakeDrag({ x: 10, y: 1, D: 1, Cd: 1 }).Cd, (1 - Math.sqrt(0.1) * Math.exp(-0.45)) ** 2) < 1e-12);
  assert.ok(rel(blevinsWakeDrag({ x: 10, D: 1, Cd: 1, Cd0: 1.2 }).Cd, 1.2 * 0.467544) < 1e-5);
});

test('two risers in line: wake shielding reduces the clearance; small spacing is flagged as collision', () => {
  const riser = { L: 300, w: 800, d: 0.4, Cd: 1.2, topTension: 2e6 };
  const none = clearanceUnderCurrent({ upstream: riser, downstream: riser, current: 0.8, spacing: 6, wakeModel: 'none' });
  // identical lines without shielding deflect equally: clearance = spacing - D everywhere
  for (const i of [0, 77, 200]) assert.ok(Math.abs(none.clearance[i] - 5.6) < 1e-9);
  const huse = clearanceUnderCurrent({ upstream: riser, downstream: riser, current: 0.8, spacing: 6 });
  assert.ok(huse.converged);
  assert.ok(huse.minClearance < 5.6 - 0.3 && huse.minClearance > 4.5 && !huse.collision);
  assert.ok(huse.downstream.yMax < huse.upstream.yMax);
  // consistency: converged shielding factor equals the Huse factor at the converged spacing
  const i = 100;
  const x = 6 + huse.downstream.y[i] - huse.upstream.y[i];
  assert.ok(Math.abs(huse.dragFactor[i] - huseWake({ x, D: 0.4, Cd: 1.2 }).dragFactor) < 1e-4);
  const blev = clearanceUnderCurrent({ upstream: riser, downstream: riser, current: 0.8, spacing: 6, wakeModel: 'blevins' });
  assert.ok(blev.minClearance < 5.6);
  // with a lighter-tensioned upstream riser and 0.9 m spacing the upstream riser reaches the downstream one
  const close = clearanceUnderCurrent({ upstream: { ...riser, topTension: 1.0e6 }, downstream: riser, current: 1.0, spacing: 0.9 });
  assert.equal(close.collision, true);
  assert.ok(close.minClearance <= 0);
});

test('impact velocity estimate: v = sqrt((2 f g - k g^2)/m) (hand calculation)', () => {
  // f = 300 N/m, g = 2 m, k = 100 N/m2, m = 400 kg/m -> E = 600 - 200 = 400 J/m, v = sqrt(2) m/s
  const r = impactVelocityEstimate({ force: 300, gap: 2, mass: 400, stiffness: 100 });
  assert.ok(rel(r.velocity, Math.SQRT2) < 1e-14 && r.energy === 400 && r.reaches);
  // restoring force stops the line before contact: f g < k g^2 / 2
  assert.equal(impactVelocityEstimate({ force: 300, gap: 7, mass: 400, stiffness: 100 }).reaches, false);
  assert.throws(() => impactVelocityEstimate({ force: 1, gap: 1, mass: 0 }), /mass/);
});
