import test from 'node:test';
import assert from 'node:assert/strict';
import * as I from '../../src/engine/pipeline/integrity.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300), `${a} vs ${b}`);

// worked example used throughout: 24 in x 12.7 mm, X52 (SMYS 359 MPa), defect 5 mm deep x 200 mm long
const ex = { D: 0.6096, t: 0.0127, d: 0.005, L: 0.2, smys: 359e6 };

test('ASME B31G original: worked example by hand', () => {
  // z = 0.2^2/(0.6096*0.0127) = 0.04/0.00774192 = 5.16668 ; M = sqrt(1 + 0.8*5.16668) = 2.26569
  // d/t = 0.393701 ; (2/3)d/t = 0.262467
  // SF = 1.1*359*(1 - 0.262467)/(1 - 0.262467/2.26569) = 394.9*0.737533/0.884156 = 329.41 MPa
  // PF = 2*329.41*0.0127/0.6096 = 13.726 MPa ; PS = PF/1.39 = 9.874 MPa
  const r = I.b31gOriginal(ex);
  close(r.z, 5.16668, 1e-5);
  close(r.M, 2.26569, 1e-5);
  close(r.failureStress, 329.41e6, 1e-4);
  close(r.failurePressure, 13.726e6, 1e-4);
  close(r.safePressure, 13.726e6 / 1.39, 1e-4);
  // long defect (z > 20): SF = 1.1 SMYS (1 - d/t)
  const long = I.b31gOriginal({ ...ex, L: 0.5 });
  assert.ok(long.z > 20);
  close(long.failureStress, 394.9e6 * (1 - 0.393701), 1e-5);
  // no defect: flow stress
  close(I.b31gOriginal({ ...ex, d: 0 }).failureStress, 394.9e6, 1e-12);
  assert.equal(I.b31gOriginal({ ...ex, d: 0.011 }).withinDepthLimits, false);
  assert.throws(() => I.b31gOriginal({ ...ex, d: 0.02 }), RangeError);
});

test('Modified B31G (0.85 dL): worked example by hand', () => {
  // M = sqrt(1 + 0.6275*5.16668 - 0.003375*5.16668^2) = sqrt(1 + 3.24209 - 0.09009) = sqrt(4.15200) = 2.03764
  // 0.85 d/t = 0.334646 ; SF = (359 + 69)*(1 - 0.334646)/(1 - 0.334646/2.03764) = 428*0.665354/0.835768 = 340.73 MPa
  // PF = 2*340.73*0.0127/0.6096 = 14.197 MPa
  const r = I.b31gModified(ex);
  close(r.M, 2.03764, 1e-5);
  close(r.flowStress, 428e6, 1e-12);
  close(r.failureStress, 340.73e6, 1e-4);
  close(r.failurePressure, 14.197e6, 1e-4);
  // z > 50 branch: M = 0.032 z + 3.3 ; the two forms meet at z = 50 (4.8926 vs 4.9000)
  const L50 = Math.sqrt(50 * ex.D * ex.t);
  close(I.b31gModified({ ...ex, L: L50 * 1.0000001 }).M, 4.9, 1e-6);
  close(I.b31gModified({ ...ex, L: L50 }).M, Math.sqrt(1 + 31.375 - 8.4375), 1e-9);
  close(I.b31gModified({ ...ex, L: 1.2 }).M, 0.032 * 1.44 / (ex.D * ex.t) + 3.3, 1e-12);
});

test('effective area (RSTRENG): uniform profile equals the closed form; a local deep pit is found by the sub-length search', () => {
  // flat-bottomed defect: A/A0 = d/t for every sub-length, the full length governs (largest M)
  const x = Array.from({ length: 21 }, (_, i) => i * 0.01);
  const flat = I.effectiveArea({ D: ex.D, t: ex.t, x, depth: x.map(() => 0.005), smys: ex.smys });
  const M = Math.sqrt(1 + 0.6275 * 5.16668 - 0.003375 * 5.16668 ** 2);
  close(flat.failureStress, 428e6 * (1 - 0.393701) / (1 - 0.393701 / M), 1e-5);
  close(flat.effectiveLength, 0.2, 1e-12);
  // profile: long shallow corrosion (1 mm over 1 m) with a short 6 mm pit in the middle
  const xs = Array.from({ length: 101 }, (_, i) => i * 0.01);
  const dep = xs.map((v) => (Math.abs(v - 0.5) <= 0.0301 ? 0.006 : 0.001));
  const r = I.effectiveArea({ D: ex.D, t: ex.t, x: xs, depth: dep, smys: ex.smys });
  // less conservative than treating the whole length at the maximum depth with 0.85 dL
  const mod = I.b31gModified({ ...ex, d: 0.006, L: 1.0 });
  assert.ok(r.failureStress > mod.failureStress);
  // not above any single sub-length evaluated by hand: the pit alone (0.06 m at 6 mm)
  const zp = 0.06 ** 2 / (ex.D * ex.t); const Mp = Math.sqrt(1 + 0.6275 * zp - 0.003375 * zp * zp);
  const pit = 428e6 * (1 - 6 / 12.7) / (1 - (6 / 12.7) / Mp);
  assert.ok(r.failureStress <= pit * (1 + 1e-12));
  assert.ok(r.effectiveLength >= 0.06 && r.effectiveLength <= 1.0);
  close(r.maxDepth, 0.006, 1e-12);
});

test('DNV-RP-F101 single defect: capacity equation by hand, Part B and Part A formats', () => {
  // Q = sqrt(1 + 0.31*L^2/(D t)) = sqrt(1 + 0.31*5.16668) = sqrt(2.60167) = 1.61297
  // Pf = 2*0.0127*455e6/(0.6096 - 0.0127)*(1 - 0.393701)/(1 - 0.393701/1.61297) = 19.3617e6*0.606299/0.755915 = 15.530 MPa
  const r = I.dnvF101Single({ ...ex, fu: 455e6 });
  close(r.Q, 1.61297, 1e-5);
  close(r.failurePressure, 15.530e6, 1e-4);
  close(r.safePressure, 0.9 * 0.72 * r.failurePressure, 1e-12);
  assert.equal(r.allowablePartA, null);
  // Part A with caller factors: (d/t)* = 0.393701 + 1.0*0.08 = 0.473701 ; gamma_d*(d/t)* = 0.61581
  const a = I.dnvF101Single({ ...ex, fu: 455e6, gammaM: 0.74, gammaD: 1.3, epsD: 1.0, stdDepth: 0.08 });
  const rs = 1.3 * (0.393701 + 0.08);
  close(a.allowablePartA, 0.74 * 19.3617e6 * (1 - rs) / (1 - rs / 1.61297), 1e-4);
  assert.equal(I.dnvF101Single({ ...ex, d: 0.011, fu: 455e6, gammaM: 0.74, gammaD: 1.3, epsD: 1, stdDepth: 0.08 }).allowablePartA, 0);
  // defect-free pipe: burst = 2 t fu/(D - t)
  close(I.dnvF101Single({ ...ex, d: 0, fu: 455e6 }).failurePressure, 2 * 0.0127 * 455e6 / 0.5969, 1e-12);
  // very long defect: Q -> large, capacity -> remaining ligament
  close(I.dnvF101Single({ ...ex, L: 500, fu: 455e6 }).failurePressure, 19.3617e6 * 0.606299, 2e-3);
});

test('remaining life: depth-limited and pressure-limited cases, re-inspection and de-rating', () => {
  // low MAOP: the 80 % depth limit governs: (0.8*12.7 - 5)/0.4 mm/yr = 12.9 years
  const a = I.remainingLife({ assess: I.b31gModified, defect: ex, maop: 1e6, depthRate: 0.0004 });
  assert.equal(a.governing, 'depth');
  close(a.remainingLife, 12.9, 1e-6);
  close(a.reinspectionInterval, 6.45, 1e-6);
  // MAOP close to today's safe pressure: pressure governs, and the safe pressure at that time equals the MAOP
  const now = I.b31gModified(ex).safePressure;
  const b = I.remainingLife({ assess: I.b31gModified, defect: ex, maop: 0.9 * now, depthRate: 0.0004, lengthRate: 0.005, interval: 3 });
  assert.equal(b.governing, 'pressure');
  assert.ok(b.remainingLife > 0 && b.remainingLife < 12.9);
  const atEnd = I.b31gModified({ ...ex, d: ex.d + 0.0004 * b.remainingLife, L: ex.L + 0.005 * b.remainingLife }).safePressure;
  close(atEnd, 0.9 * now, 1e-6);
  close(b.deratedPressure, I.b31gModified({ ...ex, d: 0.0062, L: 0.215 }).safePressure, 1e-12);
  assert.ok(b.deratedPressure < now);
  // already unacceptable
  const c = I.remainingLife({ assess: I.b31gModified, defect: ex, maop: 2 * now, depthRate: 0.0004 });
  assert.equal(c.remainingLife, 0); assert.equal(c.acceptableNow, false);
  // caller growth law
  const d = I.remainingLife({ assess: I.b31gModified, defect: ex, maop: 1e6, depthRate: (y) => 0.005 + 0.001 * Math.sqrt(y) });
  close(d.remainingLife, ((0.8 * 0.0127 - 0.005) / 0.001) ** 2, 1e-5);
});

test('dent strain per ASME B31.8 Appendix R: hand calculation', () => {
  // t = 10 mm, R0 = 0.3 m, R1 = -0.2 m, R2 = -1.0 m, d = 20 mm, L = 0.4 m
  // eps1 = 0.005*(1/0.3 + 1/0.2) = 0.0416667 ; eps2 = 0.005 ; eps3 = 0.5*(0.05)^2 = 0.00125
  // eps_i = sqrt(eps1^2 - eps1*0.00625 + 0.00625^2) = sqrt(1.73611e-3 - 2.60417e-4 + 3.90625e-5) = 0.038920
  // eps_o = sqrt(eps1^2 + eps1*(-0.00375) + 0.00375^2) = sqrt(1.73611e-3 - 1.5625e-4 + 1.40625e-5) = 0.039924
  const r = I.dentStrain({ t: 0.01, R0: 0.3, R1: -0.2, R2: -1.0, d: 0.02, L: 0.4 });
  close(r.eps1, 0.0416667, 1e-5); close(r.eps2, 0.005, 1e-12); close(r.eps3, 0.00125, 1e-12);
  close(r.epsInside, 0.038920, 1e-4);
  close(r.epsOutside, 0.039924, 1e-4);
  close(r.uc, r.epsOutside / 0.06, 1e-12);
  assert.equal(I.dentDepthCheck({ d: 0.03, D: 0.6096 }).acceptable, true);       // 4.9 % < 6 %
  assert.equal(I.dentDepthCheck({ d: 0.03, D: 0.6096, onWeld: true }).acceptable, false);
  assert.equal(I.dentDepthCheck({ d: 0.005, D: 0.6096, withStressRaiser: true }).acceptable, false);
});

test('wall-loss extreme value: Gumbel parameters recovered from exact quantiles; return level formula', () => {
  // maxima placed at the plotting-position quantiles of Gumbel(loc = 2 mm, scale = 0.4 mm)
  const n = 400;
  const x = Array.from({ length: n }, (_, i) => 0.002 - 0.0004 * Math.log(-Math.log((i + 0.5) / n)));
  const r = I.wallLossExtreme({ blockMaxima: x, returnPeriod: 100, t: 0.0127 });
  close(r.scale, 0.0004, 0.03);
  close(r.loc, 0.002, 0.01);
  close(r.predictedMax, r.loc - r.scale * Math.log(-Math.log(0.99)), 1e-12);
  close(r.predictedMax, 0.002 + 0.0004 * 4.60015, 0.02);
  close(r.remainingWall, 0.0127 - r.predictedMax, 1e-12);
  assert.ok(r.exceedanceOfWall >= 0 && r.exceedanceOfWall < 1e-6);
  assert.throws(() => I.wallLossExtreme({ blockMaxima: [1, 2, 3], returnPeriod: 10 }), RangeError);
});

test('fitness-for-service summary and crack hook', () => {
  const f = I.fitnessForService({ ...ex, maop: 9.9e6, smts: 455e6 });
  assert.equal(f.levels.length, 3);
  assert.equal(f.levels[0].acceptable, false);        // B31G original: 9.874 MPa < 9.9 MPa
  assert.equal(f.levels[1].acceptable, true);         // Modified B31G: 10.21 MPa
  assert.equal(f.acceptable, true);
  assert.match(f.governingLevel, /Modified/);
  const bad = I.fitnessForService({ ...ex, d: 0.0105, maop: 5e6 });
  assert.equal(bad.acceptable, false);
  assert.match(bad.action, /80 %/);
  const c = I.crackAssessment({ flaw: { a: 0.003 }, assess: (fl) => ({ acceptable: fl.a < 0.005, uc: fl.a / 0.005 }) });
  assert.equal(c.acceptable, true); close(c.uc, 0.6, 1e-12);
  assert.throws(() => I.crackAssessment({ flaw: {}, assess: null }), RangeError);
});
