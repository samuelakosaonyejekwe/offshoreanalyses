import test from 'node:test';
import assert from 'node:assert/strict';
import * as dep from '../../src/engine/econ/depreciation.js';

const close = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `got ${a}, expected ${b}`);
const sum = (a) => a.reduce((s, v) => s + v, 0);

test('straight line: (basis - salvage) / life each year, book value ends at salvage', () => {
  const r = dep.straightLine({ basis: 10000, salvage: 1000, life: 5 });
  assert.deepEqual(r.charges, [1800, 1800, 1800, 1800, 1800]);
  close(r.bookValues[4], 1000);
});

test('double declining balance with switch to straight line (hand schedule 4000, 2400, 1440, 1080, 1080)', () => {
  const r = dep.decliningBalance({ basis: 10000, life: 5, factor: 2 });
  const expected = [4000, 2400, 1440, 1080, 1080];
  r.charges.forEach((c, k) => close(c, expected[k], 1e-12));
  assert.equal(r.switchYear, 4);
  close(r.total, 10000);
});

test('declining balance without switch never goes below salvage', () => {
  const r = dep.decliningBalance({ basis: 1000, salvage: 100, life: 5, factor: 2, switchToSL: false });
  // 400, 240, 144, 86.4, then limited to reach salvage 100: 129.6 - 100 = 29.6
  [400, 240, 144, 86.4, 29.6].forEach((v, k) => close(r.charges[k], v, 1e-10));
  close(r.bookValues[4], 100);
});

test('sum-of-years-digits: 5/15, 4/15, ... of the depreciable amount', () => {
  const r = dep.sumOfYearsDigits({ basis: 16000, salvage: 1000, life: 5 });
  assert.deepEqual(r.charges.map((v) => Math.round(v)), [5000, 4000, 3000, 2000, 1000]);
});

test('units of production follows output and stops at the depreciable amount', () => {
  const r = dep.unitsOfProduction({ basis: 1000, salvage: 100, units: [10, 30, 60], totalUnits: 100 });
  [90, 270, 540].forEach((v, k) => close(r.charges[k], v));
  const capped = dep.unitsOfProduction({ basis: 1000, units: [80, 80], totalUnits: 100 });
  close(sum(capped.charges), 1000);
});

test('MACRS 5-year half-year percentages are the IRS Pub. 946 Table A-1 values', () => {
  assert.deepEqual(Array.from(dep.MACRS_HALF_YEAR[5]), [20.00, 32.00, 19.20, 11.52, 11.52, 5.76]);
  assert.deepEqual(Array.from(dep.MACRS_HALF_YEAR[3]), [33.33, 44.45, 14.81, 7.41]);
  assert.deepEqual(Array.from(dep.MACRS_HALF_YEAR[7]), [14.29, 24.49, 17.49, 12.49, 8.93, 8.92, 8.93, 4.46]);
  const r = dep.macrs({ basis: 100000, recoveryPeriod: 5 });
  [20000, 32000, 19200, 11520, 11520, 5760].forEach((v, k) => close(r.charges[k], v, 1e-12));
});

test('every MACRS table sums to 100 % and has recovery period + 1 entries', () => {
  for (const n of [3, 5, 7, 10, 15, 20]) {
    assert.equal(dep.MACRS_HALF_YEAR[n].length, n + 1);
    close(sum(dep.MACRS_HALF_YEAR[n]), 100, 1e-9);
  }
  assert.throws(() => dep.macrs({ basis: 1, recoveryPeriod: 4 }), RangeError);
});

test('MACRS 5-year table is reproduced by 200 % DB with half-year convention and SL switch', () => {
  // year 1: half of 40 %; then 40 % of the balance until straight line over the remaining 2.5 years is larger
  let bv = 100; const out = [];
  out.push(bv * 0.4 / 2); bv -= out[0];
  for (let y = 2; y <= 5; y++) { const remaining = 5.5 - (y - 1); const c = Math.max(bv * 0.4, bv / remaining); out.push(c); bv -= c; }
  out.push(bv);
  out.forEach((v, k) => close(v, dep.MACRS_HALF_YEAR[5][k], 1e-3));
});

test('writing-down allowance pool: 18 % reducing balance with a first-year allowance', () => {
  const r = dep.writingDownAllowance({ additions: [1000, 0, 0], rate: 0.18 });
  [180, 147.6, 121.032].forEach((v, k) => close(r.charges[k], v, 1e-12));
  close(r.residualPool, 1000 * Math.pow(0.82, 3));
  const f = dep.writingDownAllowance({ additions: [1000, 500], rate: 0.20, firstYearAllowance: 0.5 });
  // year 1: FYA 500, pool 500, WDA 100 -> 600; year 2: FYA 250, pool 400 + 250 = 650, WDA 130 -> 380
  close(f.charges[0], 600); close(f.charges[1], 380);
});

test('bonus depreciation front-loads the first year and still totals the basis', () => {
  const r = dep.bonusDepreciation({ basis: 1000, bonusFraction: 0.4, method: { method: 'sl', life: 3 } });
  [600, 200, 200].forEach((v, k) => close(r.charges[k], v));
  close(dep.depreciationSchedule({ method: 'macrs', recoveryPeriod: 5, basis: 1000, bonusFraction: 0.5 })[0], 500 + 100);
});

test('dispatcher and tax-shield present value by hand', () => {
  assert.deepEqual(dep.depreciationSchedule({ method: 'schedule', basis: 200, fractions: [0.5, 0.3, 0.2] }), [100, 60, 40]);
  assert.deepEqual(dep.depreciationSchedule({ method: 'none', basis: 200 }), []);
  assert.throws(() => dep.depreciationSchedule({ method: 'magic', basis: 1 }), RangeError);
  // 200 a year for 5 years, 30 % tax, 10 %: 60 x (P/A, 10 %, 5) = 60 x 3.790787
  const ts = dep.taxShieldPV({ charges: [200, 200, 200, 200, 200], taxRate: 0.3, rate: 0.10 });
  close(ts.pv, 60 * 3.7907868, 1e-7); close(ts.undiscounted, 300);
});
