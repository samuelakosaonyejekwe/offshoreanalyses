import test from 'node:test';
import assert from 'node:assert/strict';
import {
  SN_CURVE_NAMES, SN_JOINT_CURVE_NAMES, SN_EDITIONS, snCurveLabel, snCurve, snCurveFromDetailCategory, thicknessFactor, meanStressFactor, snConstants, cyclesToFailure,
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

// ---------------------------------------------------------------------------------------------------------------
// Curve provenance: every curve carries the constants of the standard and edition it is named after.
// Sources (quotes in src/library/factor-evidence.json): DNV-RP-C203 October 2011 Tables 2-1 to 2-3; DNVGL-RP-0005
// June 2014 Table 2-1; DNVGL-RP-C203 April 2016 and September 2019 Tables 2-1 to 2-4 and §2.4.3; API RP 2A-WSD
// 21st ed. Supplement 3 Table 5.5.1-1 and §5.5.2; ISO 19902:2007 Table 16.11-1 and §16.11.4.
// ---------------------------------------------------------------------------------------------------------------

test('curve lists, labels and the default edition', () => {
  assert.deepEqual([...SN_CURVE_NAMES], ['B1', 'B2', 'C', 'C1', 'C2', 'D', 'E', 'F', 'F1', 'F3', 'G', 'W1', 'W2', 'W3', 'T']);
  assert.deepEqual([...SN_JOINT_CURVE_NAMES], ['API-WJ', 'API-CJ', 'ISO-TJ', 'ISO-CJ']);
  assert.deepEqual([...SN_EDITIONS], ['2011', '2014', '2016']);
  assert.equal(snCurve('T').edition, '2016');
  assert.equal(snCurve('T').standard, 'DNV-RP-C203');
  assert.equal(snCurve('ISO-TJ').edition, null);
  assert.equal(snCurveLabel('T'), 'DNV-RP-C203 T');
  assert.equal(snCurveLabel('API-WJ'), 'API RP 2A-WSD WJ');
  assert.equal(snCurveLabel('ISO-CJ'), 'ISO 19902 CJ');
  assert.throws(() => snCurve('T', 'air', { edition: '2024' }), /edition/);
  assert.throws(() => snCurve('X'), /unknown S-N curve/);
});

test('DNV T curve, 2016 / 2019 editions: 12.48 / 16.13, 67.09 MPa, 16 mm, k = 0.25, in all three environments', () => {
  const air = snCurve('T'); const cp = snCurve('T', 'seawater-cp'); const fc = snCurve('T', 'free-corrosion');
  assert.deepEqual([air.m1, air.logA1, air.m2, air.logA2, air.k, air.kScfAbove10, air.tRef], [3, 12.48, 5, 16.13, 0.25, null, 0.016]);
  assert.equal(air.fatigueLimitTabulated, 67.09 * MPa);
  assert.deepEqual([cp.m1, cp.logA1, cp.m2, cp.logA2, cp.k, cp.tRef], [3, 12.18, 5, 16.13, 0.25, 0.016]);
  assert.deepEqual([fc.m1, fc.logA1, fc.m2, fc.k, fc.tRef], [3, 12.03, null, 0.25, 0.016]);
  // air, 100 MPa: N = 10^12.48/100³ = 3.01995e6 at t <= 16 mm; at 32 mm the range is 100·2^0.25, N = 3.01995e6/2^0.75 = 1.79567e6
  close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.016 }), 3.01995e6, 1e-5);
  close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.032 }), 1.79567e6, 1e-5);
  // an SCF above 10 does not change the exponent in these editions: 10 MPa × 12 = 120 MPa at 32 mm
  close(cyclesToFailure({ curve: air, stressRange: 10 * MPa, scf: 12, thickness: 0.032 }), 10 ** 12.48 / (120 * 2 ** 0.25) ** 3, 1e-12);
  // cathodic protection: 100 MPa is above the slope change (94.4 MPa at 1.8e6): N = 10^12.18/100³ = 1.51356e6;
  // 50 MPa is below it: N = 10^16.13/50⁵ = 4.31668e7
  close(cyclesToFailure({ curve: cp, stressRange: 100 * MPa }), 1.51356e6, 1e-5);
  close(cyclesToFailure({ curve: cp, stressRange: 50 * MPa }), 4.31668e7, 1e-5);
  close(cp.kneeCycles, 1.8e6, 2e-3);
  // free corrosion: N = 10^12.03/100³ = 1.07152e6, single slope
  close(cyclesToFailure({ curve: fc, stressRange: 100 * MPa }), 1.07152e6, 1e-5);
  close(cyclesToFailure({ curve: fc, stressRange: 10 * MPa }), 1.07152e9, 1e-5);
});

test('DNV T curve, 2011 and 2014 editions: equal to D with 32 mm and k = 0.25 (0.30 for SCF > 10)', () => {
  for (const edition of ['2011', '2014']) {
    const air = snCurve('T', 'air', { edition }); const D = snCurve('D', 'air', { edition });
    assert.deepEqual([air.m1, air.logA1, air.m2, air.logA2, air.k, air.kScfAbove10, air.tRef], [3, 12.164, 5, 15.606, 0.25, 0.30, 0.032]);
    assert.deepEqual([air.logA1, air.logA2, air.fatigueLimitTabulated], [D.logA1, D.logA2, D.fatigueLimitTabulated]);
    assert.equal(air.fatigueLimitTabulated, 52.63 * MPa);
    assert.equal(snCurve('T', 'seawater-cp', { edition }).logA1, 11.764);
    assert.equal(snCurve('T', 'free-corrosion', { edition }).logA1, 11.687);
    close(snCurve('T', 'seawater-cp', { edition }).kneeCycles, 1e6, 3e-3);
    // 100 MPa: N = 10^12.164/100³ = 1.45881e6 up to 32 mm; at 64 mm N = 1.45881e6/2^0.75 = 8.67416e5
    close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.032 }), 1.45881e6, 1e-5);
    close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.064 }), 8.67416e5, 1e-5);
    // SCF = 12 > 10: exponent 0.30 → range 10·12·2^0.30 = 147.74 MPa, N = 10^12.164/147.74³ = 4.52407e5
    close(cyclesToFailure({ curve: air, stressRange: 10 * MPa, scf: 12, thickness: 0.064 }), 4.52407e5, 1e-5);
    // SCF = 10 keeps 0.25; an explicit k wins over both
    close(cyclesToFailure({ curve: air, stressRange: 12 * MPa, scf: 10, thickness: 0.064 }), 10 ** 12.164 / (120 * 2 ** 0.25) ** 3, 1e-12);
    close(cyclesToFailure({ curve: air, stressRange: 10 * MPa, scf: 12, thickness: 0.064, k: 0.2 }), 10 ** 12.164 / (120 * 2 ** 0.2) ** 3, 1e-12);
    close(cyclesToFailure({ curve: snCurve('T', 'seawater-cp', { edition }), stressRange: 100 * MPa }), 5.80764e5, 1e-5);
    close(cyclesToFailure({ curve: snCurve('T', 'free-corrosion', { edition }), stressRange: 100 * MPa }), 4.86407e5, 1e-5);
  }
  // the 2016 curve is higher at every thickness: 2.07 times the life at 16 mm, 10^(0.316 − 0.75 log 2) = 1.23 times at 32 mm and above
  const r = (t) => cyclesToFailure({ curve: snCurve('T'), stressRange: 100 * MPa, thickness: t }) / cyclesToFailure({ curve: snCurve('T', 'air', { edition: '2011' }), stressRange: 100 * MPa, thickness: t });
  close(r(0.016), 10 ** 0.316, 1e-12); close(r(0.032), 10 ** 0.316 / 2 ** 0.75, 1e-12); close(r(0.064), 10 ** 0.316 / 2 ** 0.75, 1e-12);
});

test('thickness exponents of C, C1, C2 by edition and environment; life at 100 MPa and t = 50 mm', () => {
  const k = (name, env, edition) => snCurve(name, env, { edition }).k;
  for (const env of ['air', 'seawater-cp']) {
    assert.deepEqual(['C', 'C1', 'C2'].map((n) => k(n, env, '2016')), [0.05, 0.10, 0.15]);
    assert.deepEqual(['C', 'C1', 'C2'].map((n) => k(n, env, '2014')), [0.05, 0.10, 0.15]);
    assert.deepEqual(['C', 'C1', 'C2'].map((n) => k(n, env, '2011')), [0.15, 0.15, 0.15]);
  }
  // free corrosion: 0.15 in every edition
  for (const edition of SN_EDITIONS) assert.deepEqual(['C', 'C1', 'C2'].map((n) => k(n, 'free-corrosion', edition)), [0.15, 0.15, 0.15]);
  // the other curves do not depend on the edition
  for (const name of ['B1', 'B2', 'D', 'E', 'F', 'F1', 'F3', 'G', 'W1', 'W2', 'W3']) {
    for (const env of ['air', 'seawater-cp', 'free-corrosion']) assert.deepEqual(snCurve(name, env, { edition: '2011' }), { ...snCurve(name, env), edition: '2011' });
  }
  assert.deepEqual(['B1', 'B2', 'D', 'E', 'F', 'F1', 'F3', 'G', 'W1', 'W2', 'W3'].map((n) => snCurve(n).k), [0, 0, 0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]);
  const life = (name, env, edition) => cyclesToFailure({ curve: snCurve(name, env, { edition }), stressRange: 100 * MPa, thickness: 0.05 });
  // C: 10^12.592/100³ = 3.90841e6; 2016: /2^(3·0.05) = 3.52246e6; 2011: /2^(3·0.15) = 2.86112e6
  close(life('C', 'air', '2016'), 3.52246e6, 1e-5); close(life('C', 'air', '2011'), 2.86112e6, 1e-5);
  // C1: 10^12.449/100³ = 2.81190e6; 2016: /2^0.30 = 2.28397e6; 2011: /2^0.45 = 2.05843e6
  close(life('C1', 'air', '2016'), 2.28397e6, 1e-5); close(life('C1', 'air', '2011'), 2.05843e6, 1e-5);
  // C2: 10^12.301/100³/2^0.45 = 1.46398e6 in both
  close(life('C2', 'air', '2016'), 1.46398e6, 1e-5); close(life('C2', 'air', '2011'), 1.46398e6, 1e-5);
  // free corrosion: C 10^12.115/100³/2^0.45 = 9.53974e5, C1 10^11.972/100³/2^0.45 = 6.86336e5
  close(life('C', 'free-corrosion', '2016'), 9.53974e5, 1e-5); close(life('C1', 'free-corrosion', '2016'), 6.86336e5, 1e-5);
});

test('every DNV curve: tabulated constants in air, with cathodic protection and in free corrosion; life at 100 MPa', () => {
  // [log a1 air, log a2, fatigue limit MPa, log a1 CP, log a free corrosion], m1 = 4 for B1 and B2, 3 otherwise
  const table = {
    B1: [15.117, 17.146, 106.97, 14.917, 12.436], B2: [14.885, 16.856, 93.59, 14.685, 12.262], C: [12.592, 16.320, 73.10, 12.192, 12.115],
    C1: [12.449, 16.081, 65.50, 12.049, 11.972], C2: [12.301, 15.835, 58.48, 11.901, 11.824], D: [12.164, 15.606, 52.63, 11.764, 11.687],
    E: [12.010, 15.350, 46.78, 11.610, 11.533], F: [11.855, 15.091, 41.52, 11.455, 11.378], F1: [11.699, 14.832, 36.84, 11.299, 11.222],
    F3: [11.546, 14.576, 32.75, 11.146, 11.068], G: [11.398, 14.330, 29.24, 10.998, 10.921], W1: [11.261, 14.101, 26.32, 10.861, 10.784],
    W2: [11.107, 13.845, 23.39, 10.707, 10.630], W3: [10.970, 13.617, 21.05, 10.570, 10.493],
  };
  for (const [name, [a1, a2, lim, cp, fc]] of Object.entries(table)) {
    const m1 = name.startsWith('B') ? 4 : 3;
    const air = snCurve(name);
    assert.deepEqual([air.m1, air.logA1, air.m2, air.logA2, air.tRef], [m1, a1, 5, a2, 0.025], name);
    assert.equal(air.fatigueLimitTabulated, lim * MPa);
    assert.equal(snCurve(name, 'seawater-cp').logA1, cp);
    assert.equal(snCurve(name, 'seawater-cp').logA2, a2);
    assert.deepEqual([snCurve(name, 'free-corrosion').m1, snCurve(name, 'free-corrosion').logA1, snCurve(name, 'free-corrosion').m2], [3, fc, null]);
    // 100 MPa at the reference thickness is on the first branch of every curve except B1 (limit 106.97 MPa)
    const first = 100 > air.kneeStress / MPa;
    close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.025 }), first ? 10 ** a1 / 100 ** m1 : 10 ** a2 / 100 ** 5, 1e-12, name);
    close(cyclesToFailure({ curve: snCurve(name, 'free-corrosion'), stressRange: 100 * MPa }), 10 ** fc / 1e6, 1e-12, name);
  }
  // spot values by hand: F in air 10^11.855/100³ = 7.1614e5; B1 in air (second branch) 10^17.146/100⁵ = 1.3996e7
  close(cyclesToFailure({ curve: snCurve('F'), stressRange: 100 * MPa }), 7.1614e5, 1e-4);
  close(cyclesToFailure({ curve: snCurve('B1'), stressRange: 100 * MPa }), 1.3996e7, 1e-4);
});

test('API RP 2A-WSD WJ and ISO 19902 TJ: 12.48 / 16.13, 16 mm, 0.25; environment rules of each standard', () => {
  for (const name of ['API-WJ', 'ISO-TJ']) {
    const air = snCurve(name);
    assert.deepEqual([air.m1, air.logA1, air.m2, air.logA2, air.k, air.tRef, air.fatigueLimitTabulated], [3, 12.48, 5, 16.13, 0.25, 0.016, null]);
    close(air.kneeCycles, 1e7, 0.015);
    // 100 MPa: 3.01995e6 at 16 mm, 1.79567e6 at 32 mm; 40 MPa (second branch): 10^16.13/40⁵ = 1.31735e8
    close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.016 }), 3.01995e6, 1e-5);
    close(cyclesToFailure({ curve: air, stressRange: 100 * MPa, thickness: 0.032 }), 1.79567e6, 1e-5);
    close(cyclesToFailure({ curve: air, stressRange: 40 * MPa }), 1.31735e8, 1e-5);
    // free corrosion: a factor of 3 on N on the m = 3 line, no change of slope → 10^12.48/3/100³ = 1.00665e6
    const fc = snCurve(name, 'free-corrosion');
    assert.equal(fc.m2, null);
    close(cyclesToFailure({ curve: fc, stressRange: 100 * MPa }), 1.00665e6, 1e-5);
    close(cyclesToFailure({ curve: fc, stressRange: 10 * MPa }), 1.00665e9, 1e-5);
    // the same record whatever DNV edition is asked for
    assert.deepEqual(snCurve(name, 'air', { edition: '2011' }), air);
  }
  // cathodic protection: ISO tabulates 12.18; API halves the life of the m = 3 branch, 10^12.48/2/100³ = 1.50998e6
  assert.equal(snCurve('ISO-TJ', 'seawater-cp').logA1, 12.18);
  close(cyclesToFailure({ curve: snCurve('ISO-TJ', 'seawater-cp'), stressRange: 100 * MPa }), 1.51356e6, 1e-5);
  close(cyclesToFailure({ curve: snCurve('API-WJ', 'seawater-cp'), stressRange: 100 * MPa }), 1.50998e6, 1e-5);
  for (const name of ['API-WJ', 'ISO-TJ']) {
    const cp = snCurve(name, 'seawater-cp');
    assert.equal(cp.logA2, 16.13);
    // slope change moves to 94 MPa at 1.8e6 cycles (API Figure 5.5-1; ISO Table 16.11-1)
    close(cp.kneeStress, 94 * MPa, 6e-3); close(cp.kneeCycles, 1.8e6, 8e-3);
    close(cyclesToFailure({ curve: cp, stressRange: 40 * MPa }), 1.31735e8, 1e-5);
  }
  // the DNV T curve of the 2016 edition is the ISO TJ curve in air and with cathodic protection, not in free corrosion
  for (const env of ['air', 'seawater-cp']) {
    const a = snCurve('T', env); const b = snCurve('ISO-TJ', env);
    assert.deepEqual([a.m1, a.logA1, a.m2, a.logA2, a.k, a.tRef], [b.m1, b.logA1, b.m2, b.logA2, b.k, b.tRef]);
  }
  close(snCurve('T', 'free-corrosion').logA1 - snCurve('ISO-TJ', 'free-corrosion').logA1, 12.03 - (12.48 - Math.log10(3)), 1e-9);
});

test('cast-joint curves: API CJ 15.17 (m = 4) / 17.21 (m = 5), ISO CJ 15.17 single slope; 38 mm, 0.15; air only', () => {
  const api = snCurve('API-CJ'); const iso = snCurve('ISO-CJ');
  assert.deepEqual([api.m1, api.logA1, api.m2, api.logA2, api.k, api.tRef], [4, 15.17, 5, 17.21, 0.15, 0.038]);
  assert.deepEqual([iso.m1, iso.logA1, iso.m2, iso.logA2, iso.k, iso.tRef], [4, 15.17, null, null, 0.15, 0.038]);
  close(api.kneeCycles, 1e7, 0.03);
  // 150 MPa: N = 10^15.17/150⁴ = 2.92170e6 up to 38 mm; at 76 mm N = 2.92170e6/2^(4·0.15) = 1.92760e6
  for (const c of [api, iso]) {
    close(cyclesToFailure({ curve: c, stressRange: 150 * MPa, thickness: 0.038 }), 2.92170e6, 1e-5);
    close(cyclesToFailure({ curve: c, stressRange: 150 * MPa, thickness: 0.076 }), 1.92760e6, 1e-5);
  }
  // 80 MPa: API second branch 10^17.21/80⁵ = 4.94937e7; ISO stays on m = 4, 10^15.17/80⁴ = 3.61110e7
  close(cyclesToFailure({ curve: api, stressRange: 80 * MPa }), 4.94937e7, 1e-5);
  close(cyclesToFailure({ curve: iso, stressRange: 80 * MPa }), 3.61110e7, 1e-5);
  for (const name of ['API-CJ', 'ISO-CJ']) for (const env of ['seawater-cp', 'free-corrosion']) assert.throws(() => snCurve(name, env), /air only/);
});
