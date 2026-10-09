import test from 'node:test';
import assert from 'node:assert/strict';
import * as fin from '../../src/engine/econ/finance.js';
import { npv, irr } from '../../src/engine/econ/cashflow.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const sum = (a) => a.reduce((s, v) => s + v, 0);

test('WACC by hand: 60/40, kE 12 %, kD 6 %, tax 25 % -> 9.0 % post-tax, 12.0 % pre-tax, 9.6 % vanilla', () => {
  const w = fin.wacc({ equity: 60, debt: 40, costOfEquity: 0.12, costOfDebt: 0.06, taxRate: 0.25 });
  close(w.postTax, 0.6 * 0.12 + 0.4 * 0.06 * 0.75, 1e-15); close(w.postTax, 0.09, 1e-12);
  close(w.preTax, 0.12, 1e-12); close(w.vanilla, 0.096, 1e-12);
});

test('CAPM with country risk premium and Hamada lever/unlever round trip', () => {
  close(fin.capm({ riskFree: 0.03, beta: 1.2, marketPremium: 0.05, countryRiskPremium: 0.02 }), 0.11, 1e-14);
  close(fin.capm({ riskFree: 0.03, beta: 1.2, marketPremium: 0.05, countryRiskPremium: 0.02, lambda: 0.5 }), 0.10, 1e-14);
  const bu = fin.unleverBeta({ betaEquity: 1.2, debtToEquity: 0.5, taxRate: 0.25 });
  close(bu, 1.2 / 1.375, 1e-14);
  close(fin.releverBeta({ betaAsset: bu, debtToEquity: 0.5, taxRate: 0.25 }), 1.2, 1e-14);
});

test('annuity loan: level payment 129.50 on 1000 at 5 % over 10 years, principal sums to the loan', () => {
  const s = fin.debtSchedule({ type: 'annuity', principal: 1000, rate: 0.05, tenor: 10 });
  close(s.payment[0], 1000 * 0.05 / (1 - 1.05 ** -10), 1e-12); close(s.payment[0], 129.50457, 1e-7);
  close(sum(s.principal), 1000, 1e-12); assert.equal(s.closing[9], 0);
  for (let k = 0; k < 10; k++) close(s.payment[k], s.payment[0], 1e-10);
  close(s.interest[0], 50); close(s.totalInterest, 10 * s.payment[0] - 1000, 1e-10);
  close(npv([0, ...s.payment], 0.05), 1000, 1e-12);
});

test('equal-principal, bullet and grace-period schedules', () => {
  const e = fin.debtSchedule({ type: 'equalPrincipal', principal: 1000, rate: 0.06, tenor: 5 });
  assert.deepEqual(e.principal.map(Math.round), [200, 200, 200, 200, 200]);
  close(e.interest[2], 600 * 0.06);
  const b = fin.debtSchedule({ type: 'bullet', principal: 1000, rate: 0.06, tenor: 4 });
  assert.deepEqual(b.principal, [0, 0, 0, 1000]); close(b.payment[3], 1060);
  const g = fin.debtSchedule({ type: 'annuity', principal: 1000, rate: 0.05, tenor: 6, grace: 2 });
  assert.equal(g.principal[0], 0); assert.equal(g.principal[1], 0); close(g.payment[1], 50);
  close(g.payment[2], 1000 * 0.05 / (1 - 1.05 ** -4), 1e-12); close(sum(g.principal), 1000, 1e-12);
});

test('sculpted debt gives a constant DSCR and repays exactly; sizing and principal modes are inverse', () => {
  const cfads = [120, 150, 180, 160, 140, 130, 110, 100];
  const s = fin.sculptedDebt({ cfads, rate: 0.06, dscr: 1.35 });
  const cov = fin.dscr(cfads, s.payment);
  for (const v of cov.dscr) close(v, 1.35, 1e-12);
  close(cov.min, 1.35, 1e-12); close(cov.average, 1.35, 1e-12);
  close(s.principalAmount, npv([0, ...cfads], 0.06) / 1.35, 1e-12);
  close(sum(s.principal), s.principalAmount, 1e-10); assert.equal(s.closing[7], 0);
  const back = fin.sculptedDebt({ cfads, rate: 0.06, principal: s.principalAmount });
  close(back.dscr, 1.35, 1e-12);
  // for sculpted debt the loan-life cover ratio equals the DSCR in every period
  const r = fin.coverageRatios({ cfads, opening: s.opening, rate: 0.06 });
  for (const v of r.llcr) close(v, 1.35, 1e-10);
});

test('LLCR and PLCR by hand for a two-year loan inside a three-year project', () => {
  const r = fin.coverageRatios({ cfads: [100, 100, 100], opening: [150, 80], rate: 0.10 });
  close(r.llcr[0], (100 / 1.1 + 100 / 1.21) / 150, 1e-13);
  close(r.plcr[0], (100 / 1.1 + 100 / 1.21 + 100 / 1.331) / 150, 1e-13);
  close(r.llcr[1], (100 / 1.1) / 80, 1e-13);
});

test('debt capacity: annuity profile sized so that the minimum DSCR equals the target; gearing cap binds', () => {
  const cfads = [150, 140, 120, 160, 170];
  const c = fin.debtCapacity({ cfads, rate: 0.07, tenor: 5, targetDscr: 1.3, profile: 'annuity' });
  const s = fin.debtSchedule({ type: 'annuity', principal: c.debt, rate: 0.07, tenor: 5 });
  close(fin.dscr(cfads, s.payment).min, 1.3, 1e-12);
  close(c.debt, 120 / 1.3 * (1 - 1.07 ** -5) / 0.07, 1e-12);
  const capped = fin.debtCapacity({ cfads, rate: 0.07, tenor: 5, targetDscr: 1.3, profile: 'sculpted', maxGearing: 0.5, totalFunding: 600 });
  assert.equal(capped.binding, 'gearing'); close(capped.debt, 300);
});

test('DSRA: six months of the next debt service, released at maturity', () => {
  const d = fin.dsra({ debtService: [100, 120, 80], monthsCover: 6 });
  close(d.initial, 50); assert.deepEqual(d.target, [60, 40, 0]); assert.deepEqual(d.movement, [10, -20, -40]);
  close(d.initial + sum(d.movement), 0, 1e-12);
});

test('construction financing: interest during construction capitalised (hand: 5.0 on the first draw)', () => {
  const c = fin.constructionFinancing({ costs: [100, 100], gearing: 0.5, rate: 0.10 });
  assert.deepEqual(c.debtDraw, [50, 52.5]); assert.deepEqual(c.idc, [0, 5]);
  close(c.facility, 102.5); close(c.totalEquity, 102.5); close(c.totalUses, 205);
  const f = fin.constructionFinancing({ costs: [100, 100], gearing: 0.5, rate: 0.10, arrangementFeeRate: 0.02, commitmentFeeRate: 0.01 });
  // fees are consistent with the facility they are charged on (fixed point)
  close(f.fees[0], 0.02 * f.facility, 1e-10);
  close(f.fees[1], 0.01 * (f.facility - f.balance[0]), 1e-10);
  close(f.facility, sum(f.debtDraw), 1e-10);
});

const fcff = [-600, -400, 180, 190, 200, 200, 200, 190, 180, 170, 160, 150];

test('project-finance model: gearing case — funding balances, debt repaid, equity IRR above project IRR', () => {
  const r = fin.financeProject({ fcff, firstOperatingRow: 2, debt: { gearing: 0.6, rate: 0.05, tenor: 8, type: 'annuity' }, taxRate: 0.25 });
  close(r.debt + r.equity, 1000 + r.idc + r.fees, 1e-9);
  close(sum(r.principal), r.debt, 1e-9);
  assert.ok(r.projectIrr > 0.05 && r.equityIrr > r.projectIrr);
  // equity flow identity in an operating year
  close(r.equityCashflows[4], fcff[4] + 0.25 * r.interest[4] - r.debtService[4], 1e-10);
  close(r.dscr[4], fcff[4] / r.debtService[4], 1e-12);
  close(r.equityIrr, irr(r.equityCashflows).irr, 1e-12);
});

test('project-finance model: sized to a target DSCR with sculpted repayment, DSCR constant and LLCR equal to it', () => {
  const r = fin.financeProject({ fcff, firstOperatingRow: 2, debt: { targetDscr: 1.4, maxGearing: 0.95, rate: 0.05, tenor: 9, type: 'sculpted' }, taxRate: 0.2 });
  assert.equal(r.binding, 'dscr');
  close(r.minDscr, 1.4, 1e-7); close(r.averageDscr, 1.4, 1e-7); close(r.llcr, 1.4, 1e-7);
  assert.ok(r.plcr > r.llcr);
  const capped = fin.financeProject({ fcff, firstOperatingRow: 2, debt: { targetDscr: 1.05, maxGearing: 0.5, rate: 0.05, tenor: 9, type: 'sculpted' } });
  assert.equal(capped.binding, 'gearing'); close(capped.gearing, 0.5); assert.ok(capped.minDscr > 1.05);
});

test('dividend lock-up traps cash while DSCR is below the trigger and releases it afterwards', () => {
  const weak = [-1000, 125, 125, 300, 300, 300];
  const base = { fcff: weak, firstOperatingRow: 1, debt: { gearing: 0.5, rate: 0.05, tenor: 5, type: 'annuity' } };
  const free = fin.financeProject(base);
  const lock = fin.financeProject({ ...base, debt: { ...base.debt, lockupDscr: 1.2 } });
  assert.ok(free.dscr[1] < 1.2 && free.dscr[3] > 1.2);
  assert.equal(lock.lockedUp[1], true); assert.equal(lock.equityCashflows[1], 0);
  close(sum(lock.equityCashflows), sum(free.equityCashflows), 1e-10);       // cash is deferred, not lost
  assert.ok(lock.equityIrr <= free.equityIrr);
});

test('DSRA in the model is funded by equity at financial close and returned', () => {
  const a = fin.financeProject({ fcff, firstOperatingRow: 2, debt: { gearing: 0.6, rate: 0.05, tenor: 8 } });
  const b = fin.financeProject({ fcff, firstOperatingRow: 2, debt: { gearing: 0.6, rate: 0.05, tenor: 8, dsraMonths: 6 } });
  close(b.dsraBalance[1], 0.5 * b.debtService[2], 1e-12);
  close(sum(b.equityCashflows), sum(a.equityCashflows), 1e-9);
  assert.ok(b.equityIrr < a.equityIrr);
});

test('equity cash flows, adjusted present value and refinancing by hand', () => {
  const e = fin.equityCashflows({ fcff: [-100, 60, 60], drawdown: [50, 0, 0], interest: [0, 5, 2.5], principal: [0, 25, 25], taxRate: 0.2, equityRate: 0.12, projectRate: 0.08 });
  assert.deepEqual(e.equityCashflows, [-50, 60 - 4 - 25, 60 - 2 - 25]);
  const a = fin.adjustedPresentValue({ fcff: [-100, 60, 60], unleveredRate: 0.10, interest: [0, 5, 2.5], taxRate: 0.2, debtRate: 0.05, issueCosts: 1 });
  close(a.apv, (-100 + 60 / 1.1 + 60 / 1.21) + (1 / 1.05 + 0.5 / 1.05 ** 2) - 1, 1e-13);
  const old = fin.debtSchedule({ principal: 1000, rate: 0.08, tenor: 5 });
  const rf = fin.refinance({ balance: 1000, remainingService: old.payment, newRate: 0.05, newTenor: 5, feeRate: 0.01, discountRate: 0.05 });
  close(rf.pvNewService, 1000, 1e-10); close(rf.gain, npv([0, ...old.payment], 0.05) - 1000 - 10, 1e-10);
  assert.ok(rf.gain > 0);
});

test('lease versus buy: net advantage to leasing by hand (no tax, 8 %, four payments in advance)', () => {
  const r = fin.leaseVsBuy({ price: 1000, leasePayments: [260, 260, 260, 260], debtRate: 0.08 });
  const pv = 260 * (1 + (1 - 1.08 ** -3) / 0.08);
  close(r.nal, 1000 - pv, 1e-12); assert.equal(r.decision, 'lease');
  const t = fin.leaseVsBuy({ price: 1000, leasePayments: [320, 320, 320, 320], depreciation: [250, 250, 250, 250], taxRate: 0.3, debtRate: 0.08, inAdvance: false });
  const rr = 0.08 * 0.7; const pa = (1 - (1 + rr) ** -4) / rr;
  close(t.nal, 1000 - 320 * 0.7 * pa - 0.3 * 250 * pa, 1e-12); assert.equal(t.decision, 'buy');
});

test('inflation-linked debt reduces to the real annuity at zero inflation and scales with the index', () => {
  const real = fin.debtSchedule({ principal: 1000, rate: 0.03, tenor: 5 });
  const zero = fin.inflationLinkedDebt({ principal: 1000, realRate: 0.03, tenor: 5, inflation: 0 });
  zero.payment.forEach((v, k) => close(v, real.payment[k], 1e-12));
  const idx = fin.inflationLinkedDebt({ principal: 1000, realRate: 0.03, tenor: 5, inflation: 0.02 });
  close(idx.payment[2], real.payment[2] * 1.02 ** 3, 1e-12); assert.equal(idx.closing[4], 0);
});
