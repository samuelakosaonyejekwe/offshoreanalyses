import test from 'node:test';
import assert from 'node:assert/strict';
import {
  interlockedRingStiffness, solidRingStiffness, ringElasticCollapse, combinedElasticCollapse, ovalisedRingCollapse,
  collapseAssessment, birdcagingCheck, wireFoundationBuckling, lateralBucklingScreen, overbendCheck,
} from '../../src/engine/xsection/flexcollapse.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);

test('ring elastic collapse p_cr = 3 EI/R^3 = (E/4)(t/R)^3 for a solid ring (Levy; Timoshenko & Gere Ch. 7)', () => {
  const E = 2e11, t = 0.006, R = 0.08;
  const EI = solidRingStiffness({ E, thickness: t }).EI;
  close(EI, E * t ** 3 / 12);
  close(ringElasticCollapse({ EI, R }).pCr, E / 4 * (t / R) ** 3, 1e-12);
  // long tube: plane-strain factor gives the classical 2E/(1-nu^2) (t/D)^3
  const EIps = solidRingStiffness({ E, thickness: t, nu: 0.3 }).EI;
  close(ringElasticCollapse({ EI: EIps, R }).pCr, 2 * E / (1 - 0.09) * (t / (2 * R)) ** 3, 1e-12);
});

test('I-equivalent thickness of an interlocked profile', () => {
  const r = interlockedRingStiffness({ E: 1.93e11, Imin: 6e-11, pitch: 0.016, nTendons: 1, K: 1 });
  close(r.EI, 1.93e11 * 6e-11 / 0.016);
  close(1.93e11 * r.tEq ** 3 / 12, r.EI, 1e-12);
  close(interlockedRingStiffness({ E: 1.93e11, Imin: 6e-11, pitch: 0.016, nTendons: 2, K: 0.8 }).EI, 1.6 * r.EI, 1e-12);
});

test('combined resistance is the sum of the ring stiffness terms', () => {
  const c = combinedElasticCollapse({ rings: [{ EI: 700, R: 0.08 }, { EI: 9000, R: 0.095 }] });
  close(c.pCr, 3 * 700 / 0.08 ** 3 + 3 * 9000 / 0.095 ** 3, 1e-12);
});

test('Timoshenko ovalised ring: limits for ovality -> 0 and first-yield equation satisfied', () => {
  const R = 0.08, t = 0.006, E = 2e11;
  const pCr = E / 4 * (t / R) ** 3; // 21.1 MPa
  // yield governs
  const low = ovalisedRingCollapse({ pCr, yieldStress: 200e6, R, thickness: t, ovality: 0 });
  close(low.pCollapse, Math.min(200e6 * t / R, pCr), 1e-9);
  close(low.pCollapse, 15e6, 1e-9);
  // elastic buckling governs
  const high = ovalisedRingCollapse({ pCr, yieldStress: 600e6, R, thickness: t, ovality: 0 });
  close(high.pCollapse, pCr, 1e-9);
  // finite ovality: root satisfies p R/t + 6 p R w0 /(t^2 (1 - p/pcr)) = sigma_y
  const sy = 350e6, ov = 0.004;
  const r = ovalisedRingCollapse({ pCr, yieldStress: sy, R, thickness: t, ovality: ov });
  const p = r.pCollapse, w0 = ov * R;
  close(p * R / t + 6 * p * R * w0 / (t * t * (1 - p / pCr)), sy, 1e-9);
  assert.ok(p < Math.min(sy * t / R, pCr));
  // monotonic in ovality
  assert.ok(ovalisedRingCollapse({ pCr, yieldStress: sy, R, thickness: t, ovality: 0.01 }).pCollapse < p);
});

test('wet collapse (carcass alone) is below dry collapse (carcass + pressure armour); depth and utilisation', () => {
  const carcass = { EI: 800, R: 0.08, thickness: 0.004, yieldStress: 400e6 };
  const pa = { EI: 6000, R: 0.095 };
  const wet = collapseAssessment({ carcass, pressureArmour: pa, annulus: 'flooded', ovality: 0.002, pDesign: 2e6 });
  const dry = collapseAssessment({ carcass, pressureArmour: pa, annulus: 'dry', ovality: 0.002 });
  close(wet.pElastic, 3 * 800 / 0.08 ** 3, 1e-12);
  assert.ok(dry.pCollapse > wet.pCollapse);
  close(wet.collapseDepth, wet.pCollapse / (1025 * 9.80665), 1e-12);
  close(wet.utilisation, 2e6 / (0.67 * wet.pCollapse), 1e-12);
  const bent = collapseAssessment({ carcass, annulus: 'flooded', curvatureFactor: 0.9 });
  close(bent.pCollapse, 0.9 * wet.pCollapse, 1e-12);
  assert.throws(() => collapseAssessment({ carcass, curvatureFactor: 1.2 }), /curvatureFactor/);
});

test('bird-caging: outward pressure C tan^2(a)/(2 pi R^2) against tape hoop strength, by hand', () => {
  const alpha = 35 * Math.PI / 180, R = 0.1;
  // single layer, tape 1.5 mm at 600 MPa
  const r = birdcagingCheck({ compression: 1e5, layers: [{ R, alpha }], tape: { strength: 600e6, thickness: 0.0015, R: 0.103 } });
  close(r.demandPressure, 1e5 * Math.tan(alpha) ** 2 / (2 * Math.PI * 0.01), 1e-12);
  close(r.capacityPressure, 600e6 * 0.0015 / 0.103, 1e-12);
  close(r.criticalCompression, r.capacityPressure * 2 * Math.PI * 0.01 / Math.tan(alpha) ** 2, 1e-12);
  close(r.utilisation, r.demandPressure / r.capacityPressure, 1e-12);
  // dry annulus: external pressure supports; flooded: it does not
  const dry = birdcagingCheck({ compression: 1e5, layers: [{ R, alpha }], pExt: 5e6 });
  close(dry.capacityPressure, 5e6);
  const wet = birdcagingCheck({ compression: 1e5, layers: [{ R, alpha }], pExt: 5e6, annulusFlooded: true });
  assert.equal(wet.capacityPressure, 0);
  assert.equal(wet.utilisation, Infinity);
  // two equal layers at the same radius share the load: same demand as one layer
  const two = birdcagingCheck({ compression: 1e5, layers: [{ R, alpha }, { R, alpha }], pExt: 5e6 });
  close(two.demandPressure, r.demandPressure, 1e-12);
});

test('wire on elastic foundation: 2 sqrt(k EI) is the minimum of P(l) = EI (pi/l)^2 + k (l/pi)^2', () => {
  const E = 2.07e11, In = 0.012 * 0.004 ** 3 / 12, k = 4e6;
  const r = wireFoundationBuckling({ E, In, k });
  let min = Infinity, lMin = 0;
  for (let l = 0.01; l < 1; l += 1e-5) {
    const P = E * In * (Math.PI / l) ** 2 + k * (l / Math.PI) ** 2;
    if (P < min) { min = P; lMin = l; }
  }
  close(r.Pcr, min, 1e-6);
  close(r.halfWavelength, lMin, 1e-3);
  // foundation modulus from a restraining ring
  const f = wireFoundationBuckling({ E, In, restraintE: 1e9, restraintThickness: 0.008, restraintR: 0.11, R: 0.1, alpha: 0.6, n: 40 });
  close(f.k, 1e9 * 0.008 / 0.11 ** 2 * 2 * Math.PI * 0.1 * Math.cos(0.6) / 40, 1e-12);
  close(f.layerCompression, 40 * f.Pcr * Math.cos(0.6), 1e-12);
});

test('lateral buckling screen is the Euler load over the stated buckling length', () => {
  const E = 2.07e11, It = 0.004 * 0.012 ** 3 / 12, R = 0.1, alpha = 0.6, n = 40;
  const r = lateralBucklingScreen({ E, It, R, alpha, n, compression: 1e5 });
  const Lw = Math.hypot(2 * Math.PI * R, 2 * Math.PI * R / Math.tan(alpha)); // wire length per pitch
  close(r.bucklingLength, Lw, 1e-12);
  close(r.Pcr, Math.PI ** 2 * E * It / Lw ** 2, 1e-12);
  close(r.utilisation, 1e5 / (n * r.Pcr * Math.cos(alpha)), 1e-12);
  close(lateralBucklingScreen({ E, It, R, alpha, n, lengthFactor: 0.5 }).Pcr, 4 * r.Pcr, 1e-12);
});

test('over-bending check', () => {
  const r = overbendCheck({ curvature: 0.5, lockingRadius: 1.5 });
  close(r.utilisation, 0.75); assert.equal(r.locked, false); close(r.marginRadius, 0.5);
  assert.equal(overbendCheck({ curvature: 1, lockingRadius: 1.5 }).locked, true);
  assert.throws(() => overbendCheck({ curvature: -1, lockingRadius: 1 }), /curvature/);
});
