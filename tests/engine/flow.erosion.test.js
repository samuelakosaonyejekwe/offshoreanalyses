import test from 'node:test';
import assert from 'node:assert/strict';
import * as E from '../../src/engine/flow/erosion.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('API RP 14E erosional velocity: C = 100, ρ = 100 kg/m³ (6.243 lb/ft³) → 40.02 ft/s = 12.20 m/s', () => {
  // native units: Ve = 100/√6.2428 = 40.023 ft/s × 0.3048 = 12.199 m/s
  rel(E.erosionalVelocity({ rho: 100 }), 100 / Math.sqrt(100 / 16.018463) * 0.3048, 1e-6);
  rel(E.erosionalVelocity({ rho: 100 }), 12.199, 1e-4);
  rel(E.erosionalVelocity({ rho: 400, C: 150 }), 1.5 * 12.199 / 2, 1e-4);
  assert.throws(() => E.erosionalVelocity({ rho: 0 }), RangeError);
});

test('DNV-RP-O501 impact-angle function and erosion ratio', () => {
  // normal impact: 0.6 × 1^0.6 × (1 − e^(−10π)) = 0.6; grazing: 0
  close(E.impactAngleFunction(Math.PI / 2), 0.6, 1e-9); close(E.impactAngleFunction(0), 0, 1e-15);
  // 30°: 0.6 × [0.5 + 7.2 × 0.25]^0.6 × (1 − e^(−10.472)) = 0.6 × 2.3^0.6 = 0.9891
  rel(E.impactAngleFunction(Math.PI / 6), 0.6 * 2.3 ** 0.6 * (1 - Math.exp(-20 * Math.PI / 6)), 1e-12);
  rel(E.impactAngleFunction(Math.PI / 6), 0.9891, 1e-3);
  let max = 0; let at = 0;
  for (let d = 1; d <= 90; d++) { const f = E.impactAngleFunction(d * Math.PI / 180); if (f > max) { max = f; at = d; } }
  assert.ok(max > 0.95 && max <= 1.01 && at > 20 && at < 45, `maximum ${max} at ${at}°`);
  // E = K U^n F(α): 30 m/s at 30° → 2e-9 × 30^2.6 × 0.9891 = 1.370e-5 kg/kg
  rel(E.impactErosionRatio({ Up: 30, alpha: Math.PI / 6 }), 2e-9 * 30 ** 2.6 * 0.9891, 1e-3);
  rel(E.impactErosionRatio({ Up: 60, alpha: Math.PI / 6 }) / E.impactErosionRatio({ Up: 30, alpha: Math.PI / 6 }), 2 ** 2.6, 1e-12);
});

test('straight-pipe erosion: 2.5e-5 U^2.6 D^−2 ṁ', () => {
  // 20 m/s, D = 0.1 m, 0.01 kg/s: 2.5e-5 × 2412.7 × 100 × 0.01 = 0.0603 mm/yr
  rel(E.erosionStraightPipe({ Up: 20, D: 0.1, mSand: 0.01 }), 2.5e-5 * 20 ** 2.6 / 0.01 * 0.01, 1e-12);
  rel(E.erosionStraightPipe({ Up: 20, D: 0.1, mSand: 0.01 }), 0.0603, 2e-3);
  close(E.erosionStraightPipe({ Up: 20, D: 0.1, mSand: 0 }), 0, 1e-15);
});

test('bend erosion (DNV-RP-O501) by hand: gas with sand in a 1.5 D bend', () => {
  // Up = 20 m/s, D = 0.1 m, 250 µm sand at 0.001 kg/s, gas ρ = 50 kg/m³, μ = 1.5e-5 Pa·s, Rc = 1.5
  const o = { Up: 20, D: 0.1, mSand: 0.001, dp: 250e-6, rhoM: 50, muM: 1.5e-5 };
  const r = E.erosionBend(o);
  const alpha = Math.atan(1 / Math.sqrt(3)); // 30°
  close(r.alpha, alpha, 1e-12);
  const A = 50 ** 2 * Math.tan(alpha) * 20 * 0.1 / (2650 * 1.5e-5);
  rel(r.A, A, 1e-12);
  let gc = 50 / (2650 * (1.88 * Math.log(A) - 6.04)); if (!(gc > 0 && gc < 0.1)) gc = 0.1;
  rel(r.gammaC, gc, 1e-12);
  const gamma = 250e-6 / 0.1; const G = gamma < gc ? gamma / gc : 1;
  rel(r.G, G, 1e-12);
  const F = E.impactAngleFunction(alpha); rel(r.F, F, 1e-12);
  const At = Math.PI * 0.01 / (4 * Math.sin(alpha));
  const rate = 2e-9 * F * 20 ** 2.6 / (7800 * At) * G * 2.5 * 1 * 0.001;
  rel(r.rate / rate, 3.15e10, 2e-3); // m/s → mm/year
  // proportional to sand rate and to the geometry factor; steeper than quadratic in velocity
  rel(E.erosionBend({ ...o, mSand: 0.002 }).rate, 2 * r.rate, 1e-12); rel(E.erosionBend({ ...o, GF: 3 }).rate, 3 * r.rate, 1e-12);
  assert.ok(E.erosionBend({ ...o, Up: 40 }).rate > 2 ** 2.6 * r.rate * 0.5);
  // a larger bend radius gives a shallower impact angle
  assert.ok(E.erosionBend({ ...o, Rc: 5 }).alpha < r.alpha);
});

test('Salama correlation and the limiting velocity for an allowable rate', () => {
  // W = 10 kg/day, Vm = 20 m/s, d = 250 µm, D = 100 mm, ρm = 100 kg/m³: 10 × 400 × 250/(5.5 × 10⁴ × 100) = 0.1818 mm/yr
  const s = E.erosionSalama({ mSand: 10 / 86400, Vm: 20, dp: 250e-6, D: 0.1, rhoM: 100 });
  rel(s, 10 * 400 * 250 / (5.5 * 1e4 * 100), 1e-9); rel(s, 0.1818, 1e-3);
  // fixed sand rate, Salama: ER ∝ V² → V = 20 √(0.1/0.1818) = 14.83 m/s
  const v = E.maxVelocityForErosion({ model: 'salama', allowable: 0.1, mSand: 10 / 86400, dp: 250e-6, D: 0.1, rhoM: 100 });
  rel(v.v, 20 * Math.sqrt(0.1 / s), 1e-4); rel(v.rate, 0.1, 1e-4);
  // fixed concentration: sand rate ∝ V → ER ∝ V³
  const c = 10 / 86400 / (20 * Math.PI * 0.01 / 4);
  const vc = E.maxVelocityForErosion({ model: 'salama', allowable: 0.1, sandConcentration: c, dp: 250e-6, D: 0.1, rhoM: 100 });
  rel(vc.v, 20 * (0.1 / s) ** (1 / 3), 1e-4);
  // straight pipe at fixed rate: V = (allowable D²/(2.5e-5 ṁ))^(1/2.6); a synergy factor lowers the limit
  const st = E.maxVelocityForErosion({ model: 'straight', allowable: 0.05, mSand: 0.01, D: 0.1 });
  rel(st.v, (0.05 * 0.01 / (2.5e-5 * 0.01)) ** (1 / 2.6), 1e-4);
  assert.ok(E.maxVelocityForErosion({ model: 'straight', allowable: 0.05, mSand: 0.01, D: 0.1, interactionFactor: 2 }).v < st.v);
  // never reached within the search range
  close(E.maxVelocityForErosion({ model: 'straight', allowable: 1e6, mSand: 1e-6, D: 0.5, vMaxSearch: 50 }).v, 50, 1e-9);
  const b = E.maxVelocityForErosion({ allowable: 0.1, mSand: 0.001, D: 0.1, dp: 250e-6, rhoM: 50, muM: 1.5e-5 });
  rel(E.erosionBend({ Up: b.v, D: 0.1, mSand: 0.001, dp: 250e-6, rhoM: 50, muM: 1.5e-5 }).rate, 0.1, 1e-3);
});

test('sand transport: Oroskar–Turian critical velocity and Durand deposit velocity', () => {
  // 300 µm sand, D = 0.2 m, C = 1 %, water: √(g d (s − 1)) = √(9.80665 × 3e-4 × 1.65) = 0.06967 m/s
  const o = { d: 300e-6, D: 0.2, C: 0.01, rhoL: 1000, muL: 1e-3 };
  const w = Math.sqrt(9.80665 * 300e-6 * 1.65);
  const NRe = 0.2 * 1000 * w / 1e-3;
  const vc = w * 1.85 * 0.01 ** 0.1536 * 0.99 ** 0.3564 * (300e-6 / 0.2) ** -0.378 * NRe ** 0.09 * 0.96 ** 0.30;
  rel(E.criticalVelocityOroskarTurian(o), vc, 1e-12);
  assert.ok(vc > 0.5 && vc < 3, `critical velocity ${vc} m/s`);
  // larger pipes and denser slurries need more velocity
  assert.ok(E.criticalVelocityOroskarTurian({ ...o, D: 0.4 }) > vc && E.criticalVelocityOroskarTurian({ ...o, C: 0.05 }) > vc);
  // Durand: F_L = 1.3, D = 0.2 m, s = 2.65: 1.3 √(2 × 9.80665 × 0.2 × 1.65) = 3.307 m/s
  rel(E.depositVelocityDurand({ D: 0.2, rhoL: 1000, FL: 1.3 }), 1.3 * Math.sqrt(2 * 9.80665 * 0.2 * 1.65), 1e-12);
  rel(E.depositVelocityDurand({ D: 0.2, rhoL: 1000, FL: 1.3 }), 3.307, 1e-3);
});
