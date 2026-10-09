import test from 'node:test';
import assert from 'node:assert/strict';
import {
  minerSum, minerSumCycles, minerBlocks, minerWeibull, weibullHistogram, fatigueLife, combineDamage,
  mansonHalfordDamage, extrapolateHistogram, fatigueReliabilityLognormal, medianLifeFromCharacteristic,
} from '../../src/engine/fatigue/damage.js';
import { snCurve, snConstants } from '../../src/engine/fatigue/sn.js';
import { gamma } from '../../src/engine/math/special.js';

const MPa = 1e6;
const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('Palmgren–Miner sum, hand calculation on curve D (air)', () => {
  const sn = snConstants(snCurve('D'));
  // 100 MPa: N = 1.4588e6; 60 MPa: N = 10^12.164/60³ = 6.7538e6; 40 MPa (below knee): N = 3.9418e7.
  const r = minerSum({ ranges: [100 * MPa, 60 * MPa, 40 * MPa], counts: [1e5, 1e6, 1e7], sn });
  const expected = 1e5 / (10 ** 12.164 / 1e6) + 1e6 / (10 ** 12.164 / 216000) + 1e7 / (10 ** 15.606 / 40 ** 5);
  close(r.damage, expected, 1e-12);
  close(r.damage, 0.06855 + 0.14806 + 0.25369, 1e-3);
  close(minerSumCycles([{ range: 100 * MPa, count: 1e5 }, { range: 60 * MPa, count: 1e6 }, { range: 40 * MPa, count: 1e7 }], sn), expected, 1e-12);
});

test('block programme: repetitions to failure', () => {
  const sn = { m1: 3, a1: 1e12 * MPa ** 3, m2: null, a2: null, sKnee: 0 };
  // Per repetition: 1000/1e6 + 8000/8e6... N(100) = 1e6, N(50) = 8e6 → D = 0.001 + 0.001 = 0.002
  const r = minerBlocks({ blocks: [{ stressRange: 100 * MPa, cycles: 1000 }, { stressRange: 50 * MPa, cycles: 8000 }], sn, repetitions: 100 });
  close(r.damagePerRepetition, 0.002, 1e-12);
  close(r.damage, 0.2, 1e-12);
  close(r.repetitionsToFailure, 500, 1e-12);
});

test('Weibull closed-form damage, single slope: equals n0/a·q^m·Γ(1 + m/h) and the numerical histogram sum', () => {
  const sn = snConstants(snCurve('F', 'free-corrosion'));
  const q = 12 * MPa;
  const h = 0.85;
  const n0 = 7e7;
  const cf = minerWeibull({ sn, shape: h, cycles: n0, scale: q }).damage;
  close(cf, (n0 * 12 ** 3 * gamma(1 + 3 / h)) / 10 ** 11.378, 1e-10);
  const hist = weibullHistogram({ scale: q, shape: h, cycles: n0, bins: 40000 });
  close(hist.counts.reduce((a, b) => a + b, 0), n0, 1e-9);
  close(minerSum({ ...hist, sn }).damage, cf, 1e-4);
});

test('Weibull closed-form damage, bilinear curve (incomplete gamma functions) against the numerical histogram sum', () => {
  for (const [name, env, q, h] of [['D', 'air', 9, 1.0], ['F3', 'seawater-cp', 14, 0.8], ['C', 'air', 20, 1.2], ['T', 'air', 6, 0.9]]) {
    const sn = snConstants(snCurve(name, env), { thickness: 0.04 });
    const cf = minerWeibull({ sn, shape: h, cycles: 1e8, scale: q * MPa }).damage;
    const hist = weibullHistogram({ scale: q * MPa, shape: h, cycles: 1e8, bins: 60000 });
    close(minerSum({ ...hist, sn }).damage, cf, 2e-4, `${name} ${env}`);
  }
  // Maximum-range form: q = Δσ0/(ln n0)^(1/h).
  const sn = snConstants(snCurve('D'));
  const a = minerWeibull({ sn, shape: 1, cycles: 1e8, maxStressRange: 200 * MPa });
  close(a.scale, (200 * MPa) / Math.log(1e8), 1e-12);
});

test('life, DFF utilisation, accumulated-damage start state and remaining life', () => {
  // 0.02 per year → characteristic life 50 yr; DFF 3 → design life 16.67 yr.
  const r = fatigueLife({ damageRate: 0.02, designLife: 25, dff: 3 });
  close(r.characteristicLife, 50, 1e-12);
  close(r.designLife, 50 / 3, 1e-12);
  close(r.utilisation, 1.5, 1e-12);
  assert.equal(r.ok, false);
  // 0.15 already consumed: remaining = (1/3 − 0.15)/0.02 = 9.167 yr.
  const s = fatigueLife({ damageRate: 0.02, dff: 3, accumulatedDamage: 0.15, designLife: 5 });
  close(s.remainingLife, (1 / 3 - 0.15) / 0.02, 1e-12);
  close(s.utilisation, 3 * (0.15 + 0.1), 1e-12);
  assert.equal(s.ok, true);
});

test('combination of damage from several sources is a plain sum, independent of order', () => {
  const src = [{ name: 'wave', damage: 0.05, dff: 3 }, { name: 'VIV', damage: 0.01, dff: 10 }, { name: 'installation', damage: 0.02 }, { name: 'transport', damage: 0.005, dff: 2 }];
  const a = combineDamage(src);
  const b = combineDamage(src.slice().reverse());
  close(a.damage, 0.085, 1e-12);
  close(a.factoredDamage, 0.15 + 0.1 + 0.02 + 0.01, 1e-12);
  close(a.factoredDamage, b.factoredDamage, 1e-14);
  assert.equal(a.ok, true);
  close(a.shares.reduce((s, x) => s + x.share, 0), 1, 1e-12);
});

test('Manson–Halford damage curve, two-level test: n2/N2 = 1 − (n1/N1)^((N1/N2)^0.4)', () => {
  // High → low: N1 = 1e3, N2 = 1e5, n1/N1 = 0.5 → remaining fraction 1 − 0.5^(0.01^0.4) = 1 − 0.5^0.158489 = 0.10403
  const remaining = 1 - 0.5 ** (0.01 ** 0.4);
  close(remaining, 0.10403, 1e-4);
  const just = mansonHalfordDamage([{ cycles: 500, cyclesToFailure: 1e3 }, { cycles: remaining * 1e5 * 1.001, cyclesToFailure: 1e5 }]);
  const notYet = mansonHalfordDamage([{ cycles: 500, cyclesToFailure: 1e3 }, { cycles: remaining * 1e5 * 0.999, cyclesToFailure: 1e5 }]);
  assert.equal(just.failed, true);
  assert.equal(notYet.failed, false);
  assert.ok(just.minerSum < 1, 'high-low sequence fails below a Miner sum of 1');
  // Low → high survives beyond a Miner sum of 1.
  const lowHigh = mansonHalfordDamage([{ cycles: 5e4, cyclesToFailure: 1e5 }, { cycles: 600, cyclesToFailure: 1e3 }]);
  assert.equal(lowHigh.failed, false);
  assert.ok(lowHigh.minerSum > 1);
});

test('load-history extrapolation: counts scale with duration and the Weibull fit recovers the parent distribution', () => {
  const sn = { m1: 3, a1: 1e12 * MPa ** 3, m2: null, a2: null, sKnee: 0 };
  const hist = weibullHistogram({ scale: 10 * MPa, shape: 1.1, cycles: 1e6, bins: 4000 });
  const r = extrapolateHistogram({ ranges: hist.ranges, counts: hist.counts, measuredDuration: 30, targetDuration: 3000, sn });
  close(r.cycles, 1e8, 1e-9);
  close(r.shape, 1.1, 0.01);
  close(r.scale, 10 * MPa, 0.01);
  close(r.maxRangeExpected, 10 * MPa * Math.log(1e8) ** (1 / 1.1), 0.02);
  close(r.damageScaled, 100 * minerSum({ ...hist, sn }).damage, 1e-10);
  close(r.damageWeibull, r.damageScaled, 0.02);
});

test('lognormal (Wirsching) fatigue reliability: hand calculation', () => {
  // σ² = ln(1.09) + ln(1.25) + 9·ln(1.0625) = 0.086178 + 0.223144 + 0.545622 = 0.854944, σ = 0.924632
  const r = fatigueReliabilityLognormal({ medianLife: 200, m: 3, times: [20, 200], deltaCov: 0.3, aCov: 0.5, bCov: 0.25 });
  close(r.sigmaLnT, 0.924632, 1e-5);
  close(r.beta[0], Math.log(10) / 0.924632, 1e-5);
  close(r.beta[0], 2.4903, 1e-4);
  close(r.pf[0], 0.00638, 0.01);
  assert.ok(Math.abs(r.beta[1]) < 1e-12);
  close(r.pf[1], 0.5, 1e-12);
  // Annual probability is the increment of the cumulative probability; reliability index decreases with time.
  const t = Array.from({ length: 30 }, (_, i) => i + 1);
  const s = fatigueReliabilityLognormal({ medianLife: 200, m: 3, times: t });
  close(s.annualPf.reduce((a, b) => a + b, 0), s.pf[29], 1e-10);
  for (let i = 1; i < 30; i++) assert.ok(s.beta[i] < s.beta[i - 1]);
  // Stress bias B̃ = 1.2 shortens the median time to failure by 1.2³.
  const b = fatigueReliabilityLognormal({ medianLife: 200, m: 3, times: [20], bMedian: 1.2 });
  close(b.medianTimeToFailure, 200 / 1.728, 1e-12);
});

test('median life from the characteristic (mean − 2 s) S–N curve: factor 10^(2·0.2) = 2.512', () => {
  const r = medianLifeFromCharacteristic(100);
  close(r.medianLife, 251.19, 1e-4);
  close(r.aCov, Math.sqrt(Math.exp((0.2 * Math.LN10) ** 2) - 1), 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => minerSum({ ranges: [1], counts: [1, 2], sn: {} }), RangeError);
  assert.throws(() => fatigueLife({ damageRate: -1 }), /damageRate/);
  assert.throws(() => fatigueReliabilityLognormal({ medianLife: 0, m: 3, times: [1] }), /medianLife/);
});
