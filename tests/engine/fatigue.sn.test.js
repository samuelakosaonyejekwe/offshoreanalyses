import test from 'node:test';
import assert from 'node:assert/strict';
import {
  SN_CURVE_NAMES, snCurve, snCurveFromDetailCategory, thicknessFactor, meanStressFactor, snConstants, cyclesToFailure,
  stressRangeAtCycles, weibullClosedFormDamage, weibullScaleFromMax, allowableStressRangeWeibull, dffUtilisation,
  mooringTnApi, mooringCyclesApi, mooringSnDnv, mooringCyclesDnv, mooringNominalArea, mooringSnConstants,
  chainBreakingLoad, cyclesFromConstants,
} from '../../src/engine/fatigue/sn.js';

const MPa = 1e6;
const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('DNV-RP-C203 curve D in air: hand-calculated lives on both branches', () => {
  const D = snCurve('D');
  // Upper branch: N = 10^12.164 / 100^3 = 1.4588e6
  close(cyclesToFailure({ curve: D, stressRange: 100 * MPa }), 10 ** 12.164 / 1e6, 1e-12);
  close(cyclesToFailure({ curve: D, stressRange: 100 * MPa }), 1.4588e6, 1e-4);
  // Lower branch (below the 52.63 MPa knee): N = 10^15.606 / 40^5 = 3.9418e7
  close(cyclesToFailure({ curve: D, stressRange: 40 * MPa }), 10 ** 15.606 / 40 ** 5, 1e-12);
  close(cyclesToFailure({ curve: D, stressRange: 40 * MPa }), 3.9418e7, 1e-4);
});

test('slope change: the two branches meet at the knee and the knee matches the tabulated fatigue limit', () => {
  for (const name of SN_CURVE_NAMES) {
    const c = snCurve(name);
    const sn = snConstants(c);
    const above = cyclesFromConstants(sn, sn.sKnee * (1 + 1e-9));
    const below = cyclesFromConstants(sn, sn.sKnee * (1 - 1e-9));
    close(above, below, 1e-6, name);
    // Tabulated fatigue limit at 1e7 cycles agrees with the intersection of the tabulated branches.
    close(c.kneeStress, c.fatigueLimitTabulated, 5e-3, `${name} fatigue limit`);
    close(c.kneeCycles, 1e7, 0.015, `${name} knee cycles`);
    // Seawater with cathodic protection: slope change at 1e6 cycles (1.8e6 for T).
    close(snCurve(name, 'seawater-cp').kneeCycles, name === 'T' ? 1.8e6 : 1e6, 0.015, `${name} CP knee`);
  }
});

test('environment: seawater with cathodic protection and free corrosion, curve D at 100 MPa', () => {
  close(cyclesToFailure({ curve: snCurve('D', 'seawater-cp'), stressRange: 100 * MPa }), 10 ** 11.764 / 1e6, 1e-12);
  close(cyclesToFailure({ curve: snCurve('D', 'free-corrosion'), stressRange: 100 * MPa }), 10 ** 11.687 / 1e6, 1e-12);
  // Free corrosion is single-slope: still m = 3 at a low stress range.
  close(cyclesToFailure({ curve: snCurve('D', 'free-corrosion'), stressRange: 20 * MPa }), 10 ** 11.687 / 8000, 1e-12);
  // CP and air coincide on the lower branch.
  close(cyclesToFailure({ curve: snCurve('D', 'seawater-cp'), stressRange: 30 * MPa }), cyclesToFailure({ curve: snCurve('D'), stressRange: 30 * MPa }), 1e-12);
});

test('thickness correction: D curve, t = 50 mm, k = 0.2 → life divided by 2^0.6', () => {
  const D = snCurve('D');
  close(thicknessFactor(0.05, 0.025, 0.2), 2 ** 0.2, 1e-14);
  assert.equal(thicknessFactor(0.02, 0.025, 0.2), 1);
  const n = cyclesToFailure({ curve: D, stressRange: 100 * MPa, thickness: 0.05 });
  close(n, 10 ** 12.164 / 1e6 / 2 ** 0.6, 1e-12);
  close(n, 9.6245e5, 1e-4);
  // Tubular T curve refers to 16 mm.
  assert.equal(snCurve('T').tRef, 0.016);
});

test('stress concentration factor multiplies the stress range; inverse function round-trips', () => {
  const E = snCurve('E');
  close(cyclesToFailure({ curve: E, stressRange: 50 * MPa, scf: 2 }), cyclesToFailure({ curve: E, stressRange: 100 * MPa }), 1e-12);
  for (const N of [1e4, 1e6, 5e7, 1e9]) {
    const s = stressRangeAtCycles({ curve: E, cycles: N, thickness: 0.04, scf: 1.3 });
    close(cyclesToFailure({ curve: E, stressRange: s, thickness: 0.04, scf: 1.3 }), N, 1e-10);
  }
});

test('detail-category relationship: curve D is the 90 MPa category at 2e6 cycles', () => {
  const D = snCurve('D');
  close(stressRangeAtCycles({ curve: D, cycles: 2e6 }), 90 * MPa, 2e-3);
  const fat90 = snCurveFromDetailCategory({ category: 90 * MPa });
  close(fat90.logA1, 12.164, 1e-3);
  close(fat90.logA2, 15.606, 1e-3);
  // F3 is the 56 MPa category, W3 the 36 MPa category.
  close(stressRangeAtCycles({ curve: snCurve('F3'), cycles: 2e6 }), 56 * MPa, 3e-3);
  close(stressRangeAtCycles({ curve: snCurve('W3'), cycles: 2e6 }), 36 * MPa, 3e-3);
});

test('mean-stress reduction factor (DNV-RP-C203 §2.5)', () => {
  assert.equal(meanStressFactor({ stressMax: 100, stressMin: 10 }), 1);
  close(meanStressFactor({ stressMax: -10, stressMin: -100 }), 0.6, 1e-14);
  close(meanStressFactor({ stressMax: 50, stressMin: -50 }), 0.8, 1e-14);
  close(meanStressFactor({ stressMax: 50, stressMin: -50, detail: 'welded' }), 0.9, 1e-14);
  close(meanStressFactor({ stressMax: -10, stressMin: -100, detail: 'welded' }), 0.8, 1e-14);
});

test('Weibull closed form, single slope: D = n0/a·q^m·Γ(1 + m/h) written out for h = 1, m = 3 (Γ(4) = 6)', () => {
  const sn = { m1: 3, a1: 1e12 * MPa ** 3, m2: null, a2: null, sKnee: 0 };
  const d = weibullClosedFormDamage({ sn, scale: 10 * MPa, shape: 1, cycles: 1e8 });
  close(d, (1e8 * 1000 * 6) / 1e12, 1e-10);
  close(weibullScaleFromMax(200 * MPa, 1, 1e8), (200 * MPa) / Math.log(1e8), 1e-14);
});

test('allowable stress range for a Weibull long-term distribution reproduces the target damage', () => {
  const curve = snCurve('F');
  const r = allowableStressRangeWeibull({ curve, shape: 0.9, cycles: 1e8, dff: 3, thickness: 0.04 });
  const sn = snConstants(curve, { thickness: 0.04 });
  close(weibullClosedFormDamage({ sn, scale: r.scale, shape: 0.9, cycles: 1e8 }), 1 / 3, 1e-9);
  close(r.maxStressRange, r.scale * Math.log(1e8) ** (1 / 0.9), 1e-12);
  // A higher DFF must lower the allowable range; for a single-slope curve it scales with DFF^(-1/m).
  const free = snCurve('F', 'free-corrosion');
  const a1 = allowableStressRangeWeibull({ curve: free, shape: 1, cycles: 1e8, dff: 1 }).maxStressRange;
  const a10 = allowableStressRangeWeibull({ curve: free, shape: 1, cycles: 1e8, dff: 10 }).maxStressRange;
  close(a10 / a1, 10 ** (-1 / 3), 1e-8);
  close(dffUtilisation(0.2, 3), 0.6, 1e-12);
});

test('mooring T–N curves: API RP 2SK tension-ratio form', () => {
  assert.deepEqual(mooringTnApi('chain-stud'), { M: 3.0, K: 1000 });
  // Studless chain, R = 0.1: N = 316 / 0.1^3 = 3.16e5
  close(mooringCyclesApi({ type: 'chain-studless', tensionRange: 1e6, referenceBreakingStrength: 1e7 }), 3.16e5, 1e-12);
  // Spiral strand, Lm = 0.3: K = 10^(3.25 − 3.43·0.3) = 10^2.221, N = K / 0.1^5.05
  close(mooringCyclesApi({ type: 'rope-spiral-strand', tensionRange: 1e6, referenceBreakingStrength: 1e7, meanLoadRatio: 0.3 }), 10 ** 2.221 * 10 ** 5.05, 1e-10);
  close(mooringTnApi('rope-six-strand', 0.2).K, 10 ** (3.20 - 0.558), 1e-12);
});

test('mooring S–N curves: DNV-OS-E301 nominal-stress form with corrosion allowance', () => {
  assert.deepEqual(mooringSnDnv('chain-studless'), { aD: 6.0e10, m: 3.0 });
  assert.deepEqual(mooringSnDnv('chain-stud'), { aD: 1.2e11, m: 3.0 });
  assert.deepEqual(mooringSnDnv('rope-stranded'), { aD: 3.4e14, m: 4.0 });
  assert.deepEqual(mooringSnDnv('rope-spiral'), { aD: 1.7e17, m: 4.8 });
  // 100 mm studless chain, ΔT = 1 MN: A = 2·π·0.1²/4 = 0.015708 m², s = 63.662 MPa, N = 6e10/63.662³ = 2.3255e5
  const r = mooringCyclesDnv({ type: 'chain-studless', tensionRange: 1e6, diameter: 0.1 });
  close(r.stressRange, 63.662e6, 1e-5);
  close(r.cycles, 6e10 / (1 / (2 * Math.PI * 0.0025)) ** 3, 1e-12);
  close(r.cycles, 2.3255e5, 1e-4);
  // 10 mm allowance on diameter: area ratio 0.81, life ratio 0.81³.
  const c = mooringCyclesDnv({ type: 'chain-studless', tensionRange: 1e6, diameter: 0.1, corrosionAllowance: 0.01 });
  close(c.cycles / r.cycles, 0.81 ** 3, 1e-12);
  close(mooringNominalArea({ type: 'rope-spiral', diameter: 0.1 }).area, Math.PI * 0.0025, 1e-14);
  // Neutral constants reproduce the same life from a stress in Pa.
  close(cyclesFromConstants(mooringSnConstants('chain-studless'), r.stressRange), r.cycles, 1e-10);
});

test('chain break load: R3, 76 mm → 0.0223·76²·(44 − 0.08·76) = 4884 kN', () => {
  close(chainBreakingLoad({ diameter: 0.076, grade: 'R3' }), 0.0223 * 76 * 76 * (44 - 6.08) * 1000, 1e-12);
  close(chainBreakingLoad({ diameter: 0.076, grade: 'R3' }), 4884e3, 1e-3);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => snCurve('Z'), RangeError);
  assert.throws(() => snCurve('D', 'vacuum'), RangeError);
  assert.throws(() => cyclesToFailure({ curve: snCurve('D'), stressRange: -1 }), /stressRange/);
  assert.throws(() => mooringCyclesDnv({ type: 'chain-studless', tensionRange: 1, diameter: 0.1, corrosionAllowance: 0.2 }), /corrosionAllowance/);
});
