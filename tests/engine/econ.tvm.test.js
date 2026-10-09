import test from 'node:test';
import assert from 'node:assert/strict';
import * as tvm from '../../src/engine/econ/tvm.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('standard factors match the compound-interest tables (Blank & Tarquin, 10 % / 8 % / 6 % tables)', () => {
  close(tvm.factorAP(0.10, 5), 0.26380, 2e-5);
  close(tvm.factorPA(0.08, 10), 6.7101, 1e-5);
  close(tvm.factorFP(0.06, 10), 1.7908, 3e-5);
  close(tvm.factorPG(0.10, 5), 6.8618, 1e-5);
  close(tvm.factorAF(0.10, 5), 0.16380, 2e-5);
  close(tvm.factorFA(0.10, 5), 6.1051, 1e-5);
  close(tvm.factorAG(0.10, 5), 1.8101, 2e-5);
  close(tvm.factorPF(0.10, 5), 0.6209, 2e-5);
});

test('factor identities: A/P = A/F + i, P/G = (P/A)(A/G), F/G = (P/G)(F/P)', () => {
  for (const [i, n] of [[0.03, 7], [0.12, 25], [0.2, 3]]) {
    close(tvm.factorAP(i, n), tvm.factorAF(i, n) + i, 1e-12);
    close(tvm.factorPG(i, n), tvm.factorPA(i, n) * tvm.factorAG(i, n), 1e-10);
    close(tvm.factorFG(i, n), tvm.factorPG(i, n) * tvm.factorFP(i, n), 1e-12);
    close(tvm.factorPA(i, n) * tvm.factorAP(i, n), 1, 1e-12);
  }
});

test('zero-interest limits: P/A = n, A/P = 1/n, P/G = n(n-1)/2, A/G = (n-1)/2', () => {
  close(tvm.factorPA(0, 8), 8); close(tvm.factorAP(0, 8), 1 / 8); close(tvm.factorFA(0, 8), 8);
  close(tvm.factorPG(0, 8), 28); close(tvm.factorAG(0, 8), 3.5);
  close(tvm.factorPG(1e-10, 8), 28, 1e-6);
});

test('gradient factor equals the brute-force present value of 0, G, 2G, ...', () => {
  const i = 0.07; const n = 12; let pv = 0;
  for (let k = 1; k <= n; k++) pv += (k - 1) * Math.pow(1 + i, -k);
  close(tvm.factorPG(i, n), pv, 1e-11);
});

test('geometric gradient: hand value for i = 10 %, g = 5 %, n = 3 and the g = i limit', () => {
  // 1/1.1 + 1.05/1.21 + 1.1025/1.331
  close(tvm.factorPGeometric(0.10, 0.05, 3), 1 / 1.1 + 1.05 / 1.21 + 1.1025 / 1.331, 1e-12);
  close(tvm.factorPGeometric(0.08, 0.08, 6), 6 / 1.08, 1e-12);
  close(tvm.pvGrowingAnnuity({ A1: 100, i: 0.10, g: 0.05, n: 3 }), 260.5184, 1e-6);
});

test('compounding: 12 % nominal monthly is 12.6825 % effective, continuous 12.7497 %', () => {
  close(tvm.effectiveRate(0.12, 12), 0.126825, 1e-6);
  close(tvm.effectiveRate(0.12, Infinity), Math.exp(0.12) - 1, 1e-14);
  close(tvm.nominalRate(tvm.effectiveRate(0.12, 12), 12), 0.12, 1e-12);
  close(tvm.nominalRate(0.10, Infinity), Math.log(1.1), 1e-14);
  close(tvm.compoundAmount({ P: 1000, r: 0.08, t: 3, m: 4 }), 1000 * Math.pow(1.02, 12), 1e-12);
  close(tvm.fvContinuous(1000, 0.05, 10), 1648.7212707, 1e-9);
  close(tvm.pvContinuous(tvm.fvContinuous(250, 0.07, 4), 0.07, 4), 250, 1e-12);
  close(Math.pow(1 + tvm.periodicRate(0.10, 12), 12), 1.10, 1e-12);
});

test('single sums, annuities (ordinary and due), perpetuities', () => {
  close(tvm.futureValue(1000, 0.06, 10), 1790.8477, 1e-7);
  close(tvm.presentValue(1790.8477, 0.06, 10), 1000, 1e-7);
  close(tvm.pvAnnuity({ A: 100, i: 0.08, n: 10 }), 671.00814, 1e-7);
  close(tvm.pvAnnuity({ A: 100, i: 0.08, n: 10, due: true }), 671.00814 * 1.08, 1e-7);
  close(tvm.fvAnnuity({ A: 100, i: 0.10, n: 5 }), 610.51, 1e-9);
  close(tvm.annuityPayment({ P: 1000, i: 0.10, n: 5 }), 263.79748, 1e-7);
  close(tvm.pvPerpetuity(50, 0.05), 1000);
  close(tvm.pvGrowingPerpetuity(50, 0.08, 0.03), 1000);
  assert.throws(() => tvm.pvGrowingPerpetuity(50, 0.03, 0.03), RangeError);
  // capitalised cost: 1000 now, 50 a year for ever, renewal of 200 every 10 years at 5 %
  close(tvm.capitalisedCost({ firstCost: 1000, annualCost: 50, renewalCost: 200, life: 10, i: 0.05 }),
    1000 + 50 / 0.05 + 200 * (0.05 / (Math.pow(1.05, 10) - 1)) / 0.05, 1e-12);
});

test('equivalent annual worth and cost by hand', () => {
  // NPV = -1000 + 400 (P/A, 10 %, 3) = -5.2592; AW = NPV (A/P, 10 %, 3)
  const pa = (1 - Math.pow(1.1, -3)) / 0.1;
  close(tvm.equivalentAnnualWorth([-1000, 400, 400, 400], 0.10), (-1000 + 400 * pa) / pa, 1e-12);
  // EAC = (10000 - 2000)(A/P, 10 %, 5) + 2000 x 0.10 + 500 = 2810.38
  const e = tvm.equivalentAnnualCost({ capex: 10000, annualCost: 500, salvage: 2000, i: 0.10, n: 5 });
  close(e.capitalRecovery, 2310.38, 1e-5); close(e.eac, 2810.38, 1e-5);
});

test('Fisher relation is exact and reversible; inflate/deflate are inverse', () => {
  close(tvm.realToNominal(0.05, 0.03), 0.0815, 1e-14);
  close(tvm.nominalToReal(0.0815, 0.03), 0.05, 1e-13);
  close(tvm.deflate(tvm.inflate(123, 0.025, 7), 0.025, 7), 123, 1e-13);
  close(tvm.inflate(100, 0.02, 10), 121.899442, 1e-8);
});

test('discount factors: end-year, mid-year and time-varying rates', () => {
  close(tvm.discountFactor(0.10, 3), 1 / 1.331, 1e-14);
  close(tvm.discountFactor(0.10, 3, { midYear: true }), Math.pow(1.1, -2.5), 1e-14);
  close(tvm.discountFactor(0.10, 0, { midYear: true }), 1);
  const df = tvm.discountFactors([0, 0.05, 0.10, 0.20], 4);
  close(df[3], 1 / (1.05 * 1.10 * 1.20), 1e-14);
  const mid = tvm.discountFactors(0.10, 4, { midYear: true });
  close(mid[2], Math.pow(1.1, -1.5), 1e-14);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => tvm.factorPA(-1, 5), /i must be/);
  assert.throws(() => tvm.factorAP(0.1, 0), /n must be/);
  assert.throws(() => tvm.requireNumber('x', NaN), /x must be a finite number/);
  assert.throws(() => tvm.requireSeries('s', [1, 'a']), /s\[1\]/);
});
