import test from 'node:test';
import assert from 'node:assert/strict';
import {
  armourStressRange, meanStressCorrection, snCycles, genericArmourSN, minerDamage, armourFatigueCoefficients,
  armourFatigueDamage, pressureCycleStressRange, wireSlipPerCycle, archardWear, sheathCreepIntoGaps, sheathStrainCheck, endFittingAnchorage,
  multiaxialCriterion,
} from '../../src/engine/xsection/flexfatigue.js';
import { buildFlexibleStack, exampleFlexibleDesign, axisymmetricResponse } from '../../src/engine/xsection/flexible.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);
const DEG = Math.PI / 180;
const wire = { E: 2.0e11, A: 5e-5, R: 0.1, alpha: 30 * DEG, thickness: 0.004, width: 0.0125 };

test('armour stress range: stick below slip, capped at twice the friction stress after slip (hand values)', () => {
  const f = 500; // sigma_f,max = pi R f/(2 sin a A) = 3.1416 MPa (see helix test)
  const small = armourStressRange({ ...wire, frictionPerLength: f, kappaRange: 1e-4 });
  assert.equal(small.slipped, false);
  close(small.frictionRange, 2.0e11 * 0.75 * 0.1 * 1e-4); // 1.5 MPa
  const big = armourStressRange({ ...wire, frictionPerLength: f, kappaRange: 0.05 });
  assert.equal(big.slipped, true);
  close(big.frictionRange, 2 * Math.PI * 1e6, 1e-9);
  close(big.normalBendingRange, 2.0e11 * 0.002 * 0.75 ** 2 * 0.05); // E t/2 cos^4 dk = 11.25 MPa
  close(big.transverseBendingRange, 2.0e11 * 0.00625 * 1.25 * Math.cos(30 * DEG) * 0.05);
  close(big.rangeExtrados, big.frictionRange + big.normalBendingRange);
  close(big.governing, Math.max(big.rangeExtrados, big.transverseBendingRange));
});

test('Goodman and Gerber mean-stress corrections', () => {
  close(meanStressCorrection({ stressRange: 100e6, meanStress: 350e6, uts: 1400e6 }).equivalentRange, 100e6 / 0.75);
  close(meanStressCorrection({ stressRange: 100e6, meanStress: 350e6, uts: 1400e6, method: 'gerber' }).equivalentRange, 100e6 / (1 - 1 / 16));
  close(meanStressCorrection({ stressRange: 100e6, meanStress: -200e6, uts: 1400e6 }).equivalentRange, 100e6);
  assert.throws(() => meanStressCorrection({ stressRange: 1, meanStress: 2e9, uts: 1e9 }), /meanStress/);
});

test('S-N evaluation: single slope, bilinear continuity at the knee, cut-off', () => {
  // N = 1e12 at 100 MPa with m = 3  ->  a = 1e12 * (1e8)^3 = 1e36
  const c = { a: 1e36, m: 3 };
  close(snCycles({ stressRange: 100e6, curve: c }).cycles, 1e12, 1e-9);
  close(snCycles({ stressRange: 200e6, curve: c }).cycles, 1.25e11, 1e-9);
  const bi = { a: 1e36, m: 3, m2: 5, nKnee: 1e13 };
  const sKnee = (1e36 / 1e13) ** (1 / 3);
  close(snCycles({ stressRange: sKnee, curve: bi }).cycles, 1e13, 1e-6);
  close(snCycles({ stressRange: sKnee / 2, curve: bi }).cycles, 1e13 * 32, 1e-6);
  assert.equal(snCycles({ stressRange: 10e6, curve: { ...c, limit: 20e6 } }).cycles, Infinity);
});

test('generic armour curve passes through the textbook anchor points and is labelled generic', () => {
  const uts = 1200e6;
  const c = genericArmourSN({ uts });
  assert.equal(c.generic, true); assert.match(c.label, /GENERIC/);
  close(snCycles({ stressRange: 2 * 0.9 * uts, curve: c }).cycles, 1e3, 1e-6);
  close(snCycles({ stressRange: 2 * 0.5 * uts * 1.0001, curve: c }).cycles, 1e6, 2e-3);
  assert.equal(snCycles({ stressRange: 0.99 * uts, curve: c }).cycles, Infinity);
  // corrosive environments need an explicit knock-down and have no endurance limit
  assert.throws(() => genericArmourSN({ uts, environment: 'sour' }), /environmentFactor/);
  const sour = genericArmourSN({ uts, environment: 'sour', environmentFactor: 2 });
  assert.equal(sour.limit, undefined);
  close(snCycles({ stressRange: 0.9 * uts, curve: sour }).cycles, 1e3, 1e-6);
});

test('Miner damage by hand, design fatigue factor 10', () => {
  const curve = { a: 1e36, m: 3 };
  const r = minerDamage({ bins: [{ stressRange: 100e6, cycles: 1e9 }, { stressRange: 200e6, cycles: 2.5e8 }], curve });
  close(r.damage, 1e9 / 1e12 + 2.5e8 / 1.25e11, 1e-9); // 0.001 + 0.002
  close(r.life, 1 / 0.003, 1e-9);
  close(r.lifeWithDff, 1 / 0.03, 1e-9);
  // SCF 1.2 with m = 3 multiplies damage by 1.728; corrosion factor 2 doubles it again
  const s = minerDamage({ bins: [{ stressRange: 100e6, cycles: 1e9 }], curve, scf: 1.2, corrosionFactor: 2 });
  close(s.damage, 0.001 * 1.728 * 2, 1e-9);
});

test('armour fatigue over load bins reproduces a hand-built bin', () => {
  const curve = { a: 1e36, m: 3 };
  const r = armourFatigueDamage({
    wire, bins: [{ tension: 4e5, curvatureRange: 0.05, cycles: 1e6 }],
    meanStressPerTension: 500, frictionPerTension: 1.25e-3, curve, uts: 1400e6,
  });
  // T = 400 kN -> mean 200 MPa, f = 500 N/m -> range = 2 pi MPa + 11.25 MPa (previous test)
  const range = (2 * Math.PI + 11.25) * 1e6;
  close(r.stress[0].meanStress, 200e6); close(r.stress[0].frictionPerLength, 500);
  const transverse = 2.0e11 * 0.00625 * 1.25 * Math.cos(30 * DEG) * 0.05;
  const gov = Math.max(range, transverse);
  close(r.stress[0].range, gov, 1e-9);
  const eq = gov / (1 - 200 / 1400);
  close(r.damage, 1e6 / (1e36 * eq ** -3), 1e-9);
});

test('fatigue coefficients from the flexible-pipe model are consistent with a direct solution', () => {
  const st = buildFlexibleStack(exampleFlexibleDesign('riser-6in'));
  const c = armourFatigueCoefficients(st, 6, { tensionLow: 2e5, tensionHigh: 6e5, mu: 0.1, loads: { pInt: 1e7 } });
  const r = axisymmetricResponse(st, { tension: 4.5e5, pInt: 1e7 });
  close(c.meanStress0 + c.meanStressPerTension * 4.5e5, r.layers[6].wireStress, 1e-8);
  const L = st.layers[6];
  const f = 0.1 * (r.layers[6].contactIn + r.layers[6].contactOut) * 2 * Math.PI * L.R * Math.cos(L.alpha) / L.n;
  close(c.friction0 + c.frictionPerTension * 4.5e5, f, 1e-8);
  assert.ok(c.meanStressPerTension > 0 && c.frictionPerTension > 0);
  const d = armourFatigueDamage({ ...c, wire: c.wire, bins: [{ tension: 4e5, curvatureRange: 0.02, cycles: 5e6 }], curve: genericArmourSN({ uts: L.uts, environment: 'sweet', environmentFactor: 1.5 }), uts: L.uts });
  assert.ok(d.damage > 0 && Number.isFinite(d.lifeWithDff));
});

test('pressure-armour stress range for a pressure cycle is linear in the pressure step', () => {
  const st = buildFlexibleStack(exampleFlexibleDesign('riser-6in'));
  const a = pressureCycleStressRange(st, 2, { pLow: 5e6, pHigh: 2e7, loads: { tension: 3e5 } });
  const b = pressureCycleStressRange(st, 2, { pLow: 5e6, pHigh: 1.25e7, loads: { tension: 3e5 } });
  close(a.stressRange, 2 * b.stressRange, 1e-8);
  close(a.meanStress, 0.5 * (a.stressLow + a.stressHigh));
  // thin-ring bound: the armour cannot carry more hoop load than the whole pressure step p r_i
  const L = st.layers[2];
  const hoopMetal = L.nA * Math.sin(L.alpha) ** 2 / (2 * Math.PI * L.R * Math.cos(L.alpha));
  assert.ok(a.stressRange < 1.5e7 * st.layers[1].rIn / hoopMetal && a.stressRange > 0.5 * 1.5e7 * st.layers[1].rIn / hoopMetal);
});

test('wire slip equals the integral of the stick strain along the wire; Archard wear by hand', () => {
  const R = 0.1, alpha = 30 * DEG, dk = 0.05;
  // integrate strain cos^2(a) R k cos(psi) over a quarter pitch, ds = R dpsi/sin(a)
  let slip = 0; const m = 100000;
  for (let i = 0; i < m; i++) slip += 0.75 * R * dk * Math.cos((i + 0.5) / m * Math.PI / 2) * R / 0.5 * (Math.PI / 2 / m);
  const s = wireSlipPerCycle({ R, alpha, kappaRange: dk });
  close(s.slipRange, slip, 1e-8);
  close(s.slipRange, 0.01 * 0.75 * 0.05 / 0.5);
  const w = archardWear({ wearCoefficient: 1e-15, contactPressure: 5e6, slidingPerCycle: s.slidingPerCycle, cycles: 1e8, thickness: 0.004 });
  close(w.wearDepth, 1e-15 * 5e6 * 2 * 7.5e-4 * 1e8); // 0.75 mm
  close(w.stressFactor, 0.004 / (0.004 - 7.5e-4));
  assert.equal(archardWear({ wearCoefficient: 1e-13, contactPressure: 5e6, slidingPerCycle: 1e-3, cycles: 1e10, thickness: 0.004 }).wornThrough, true);
});

test('sheath creep volume balance and strain checks', () => {
  const c = sheathCreepIntoGaps({ gapWidth: 0.002, intrusionDepth: 0.003, pitch: 0.012, sheathThickness: 0.008 });
  close(c.thicknessLoss, 0.002 * 0.003 / 0.012); // 0.5 mm
  close(c.reduction, 0.0625); close(c.utilisation, 0.0625 / 0.3);
  const s = sheathStrainCheck({ rOut: 0.1, kappaMax: 0.5, strainLimit: 0.077, kappaRange: 0.1, strainRangeLimit: 0.02 });
  close(s.strain, 0.05); close(s.utilisation, 0.05 / 0.077); close(s.rangeUtilisation, 0.5);
});

test('end-fitting anchorage: bond capacity and wedge statics', () => {
  const a = endFittingAnchorage({ bondStress: 10e6, perimeter: 0.032, embeddedLength: 0.25, wireArea: 48e-6, uts: 1400e6 });
  close(a.bondCapacity, 10e6 * 0.032 * 0.25); // 80 kN
  close(a.efficiency, 80e3 / (48e-6 * 1400e6));
  // wedge: clamping force per unit pull 1/tan(beta + phi); check against force triangle (N from axial equilibrium)
  const beta = 10 * DEG, mu = 0.2;
  const N = 1 / (Math.sin(beta) + mu * Math.cos(beta));
  const radial = N * Math.cos(beta) - mu * N * Math.sin(beta);
  const w = endFittingAnchorage({ bondStress: 10e6, perimeter: 0.032, embeddedLength: 0.25, wireArea: 48e-6, uts: 1400e6, coneHalfAngle: beta, muCone: mu, muWire: 0.1 });
  close(w.wedgeAmplification, radial, 1e-12);
  close(w.capacity, w.bondCapacity / (1 - 0.1 * radial), 1e-12);
  assert.equal(w.selfReleasing, beta > Math.atan(mu));
});

test('multiaxial criteria reduce to known uniaxial and torsion values', () => {
  const sa = 200e6;
  const uni = [sa, 0, 0, 0, 0, 0], tor = [0, 0, 0, 150e6, 0, 0];
  // Sines, uniaxial: sqrt(J2) = sa/sqrt(3); mean hydrostatic = sm/3
  const s = multiaxialCriterion({ criterion: 'sines', amplitude: uni, mean: [90e6, 0, 0, 0, 0, 0], k: 0.5, limit: 200e6 });
  close(s.shearTerm, sa / Math.sqrt(3)); close(s.normalTerm, 30e6); close(s.utilisation, (sa / Math.sqrt(3) + 15e6) / 200e6);
  // Crossland in torsion: hydrostatic stress zero -> equivalent = tau
  close(multiaxialCriterion({ criterion: 'crossland', amplitude: tor, k: 0.4, limit: 200e6 }).equivalent, 150e6);
  // Dang Van, uniaxial fully reversed: tau_a = sa/2, sigma_H,max = sa/3. With a = 3(t/f - 1/2), b = t the criterion
  // is exactly met at sa = f (bending limit)
  const t = 180e6, f = 300e6, aDV = 3 * (t / f - 0.5);
  close(multiaxialCriterion({ criterion: 'dangvan', amplitude: [f, 0, 0, 0, 0, 0], k: aDV, limit: t }).utilisation, 1, 1e-12);
  close(multiaxialCriterion({ criterion: 'dangvan', amplitude: [0, 0, 0, t, 0, 0], k: aDV, limit: t }).utilisation, 1, 1e-12);
  // Findley in pure torsion: max over planes of tau + k sigma_n = tau sqrt(1 + k^2)
  const k = 0.3;
  close(multiaxialCriterion({ criterion: 'findley', amplitude: tor, k, limit: 200e6, nGrid: 180 }).equivalent, 150e6 * Math.sqrt(1 + k * k), 1e-4);
  // Findley uniaxial: max of (sa/2)(sin 2x + k (1 + cos 2x)) = (sa/2)(k + sqrt(1 + k^2))
  close(multiaxialCriterion({ criterion: 'findley', amplitude: uni, k, limit: 200e6, nGrid: 180 }).equivalent, sa / 2 * (k + Math.sqrt(1 + k * k)), 1e-4);
  assert.throws(() => multiaxialCriterion({ criterion: 'x', amplitude: uni, k: 1, limit: 1 }), /criterion/);
});
