import test from 'node:test';
import assert from 'node:assert/strict';
import * as cf from '../../src/engine/econ/cashflow.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const flows = [-1000, 300, 400, 500, 600];

test('NPV of a hand-computed cash flow at 10 %', () => {
  // 300/1.1 + 400/1.21 + 500/1.331 + 600/1.4641 - 1000 = 388.771
  close(cf.npv(flows, 0.10), -1000 + 300 / 1.1 + 400 / 1.21 + 500 / 1.331 + 600 / 1.4641, 1e-13);
  close(cf.npv(flows, 0.10), 388.7712, 1e-6);
  close(cf.npv(flows, 0), 800);
});

test('NPV with mid-year convention and with time-varying rates', () => {
  close(cf.npv([-100, 110], 0.10, { midYear: true }), -100 + 110 / Math.sqrt(1.1), 1e-13);
  close(cf.npv([-100, 60, 60], [0, 0.05, 0.10]), -100 + 60 / 1.05 + 60 / (1.05 * 1.10), 1e-13);
});

test('IRR: unique root, closed form for a two-period annuity', () => {
  close(cf.irr([-1000, 1100]).irr, 0.10, 1e-10);
  // 60x^2 + 60x - 100 = 0, x = 1/(1+r)
  const x = (-60 + Math.sqrt(3600 + 24000)) / 120;
  const r = cf.irr([-100, 60, 60]);
  assert.equal(r.status, 'unique');
  close(r.irr, 1 / x - 1, 1e-10);
  assert.equal(r.norstromUnique, true);
  close(cf.npv([-100, 60, 60], r.irr), 0, 1e-9);
});

test('classic multiple-IRR case (-1000, 2500, -1540) gives 10 % and 40 % with status multiple', () => {
  const r = cf.irr([-1000, 2500, -1540]);
  assert.equal(r.status, 'multiple');
  assert.equal(r.irr, null);
  assert.equal(r.roots.length, 2);
  close(r.roots[0], 0.10, 1e-9); close(r.roots[1], 0.40, 1e-9);
  assert.equal(r.descartesSignChanges, 2);
  assert.equal(r.norstromUnique, false);
});

test('IRR status none: no sign change, and sign changes without a real root', () => {
  assert.equal(cf.irr([-100, -50, -20]).status, 'none');
  // 1 - 2x + 2x^2 has no real root although the flows change sign twice
  const r = cf.irr([1, -2, 2]);
  assert.equal(r.status, 'none'); assert.equal(r.descartesSignChanges, 2);
});

test('IRR finds negative rates and a tangent (double) root', () => {
  close(cf.irr([-100, 50, 30]).irr, cf.irr([-100, 50, 30]).roots[0]);
  assert.ok(cf.irr([-100, 50, 30]).irr < 0);
  close(cf.npv([-100, 50, 30], cf.irr([-100, 50, 30]).irr), 0, 1e-8);
  const d = cf.irr([-100, 200, -100]);            // -(1 - x)^2: double root at r = 0
  assert.equal(d.roots.length, 1); close(d.roots[0], 0, 1e-6);
});

test('IRR of a 30-year project matches the annuity relation', () => {
  const a = [-1000, ...new Array(30).fill(120)];
  const r = cf.irr(a).irr;
  close(120 * (1 - Math.pow(1 + r, -30)) / r, 1000, 1e-9);
});

test('MIRR by hand: FV of inflows at 12 %, PV of outflows at 10 %', () => {
  const fv = 300 * 1.12 ** 3 + 400 * 1.12 ** 2 + 500 * 1.12 + 600;      // 2083.2384
  close(fv, 2083.2384, 1e-9);
  close(cf.mirr(flows, 0.10, 0.12), Math.pow(fv / 1000, 0.25) - 1, 1e-13);
  close(cf.mirr(flows, 0.10, 0.12), 0.201392, 1e-5);
});

test('payback with interpolation: simple 2.6 years, discounted 3.0513 years', () => {
  const p = cf.payback(flows, 0.10);
  close(p.simple, 2 + 300 / 500, 1e-12);
  const cum3 = -1000 + 300 / 1.1 + 400 / 1.21 + 500 / 1.331;
  close(p.discounted, 3 + -cum3 / (600 / 1.4641), 1e-12);
  assert.equal(cf.payback([-100, 10, 10]).simple, null);
  assert.equal(cf.payback([-100, 150, -80]).staysPositive, false);
});

test('profitability index, benefit-cost ratios and ROI by hand', () => {
  close(cf.profitabilityIndex(flows, 0.10).pi, 1.3887712, 1e-6);
  const pa = (1 - 1.1 ** -5) / 0.1;
  const b = cf.benefitCostRatio({ benefits: [0, 500, 500, 500, 500, 500], investment: [1000], om: [0, 100, 100, 100, 100, 100], rate: 0.10 });
  close(b.conventional, 500 * pa / (1000 + 100 * pa), 1e-12);
  close(b.modified, (500 * pa - 100 * pa) / 1000, 1e-12);
  close(cf.roi(flows).roi, 0.8); close(cf.roi(flows).annualised, 0.2);
});

test('metric set is consistent with the individual functions', () => {
  const m = cf.cashflowMetrics(flows, 0.10);
  close(m.npv, cf.npv(flows, 0.10)); close(m.irr, cf.irr(flows).irr); assert.equal(m.irrStatus, 'unique');
  close(m.eaw, m.npv * 0.1 / (1 - 1.1 ** -4), 1e-12);
});

test('incremental analysis: challenger accepted only when the increment earns the MARR', () => {
  const pa = (1 - 1.1 ** -5) / 0.1;
  const A = { name: 'A', cashflows: [-1000, 300, 300, 300, 300, 300] };
  const B = { name: 'B', cashflows: [-1500, 430, 430, 430, 430, 430] };
  const r = cf.incrementalAnalysis({ alternatives: [B, A], marr: 0.10 });
  assert.equal(r.selected, 'A');
  close(r.steps[1].incrementalNpv, -500 + 130 * pa, 1e-12);          // -7.20: reject B
  assert.ok(r.steps[1].incrementalIrr < 0.10);
  const B2 = { name: 'B', cashflows: [-1500, 440, 440, 440, 440, 440] };
  assert.equal(cf.incrementalAnalysis({ alternatives: [A, B2], marr: 0.10 }).selected, 'B');
  assert.equal(cf.incrementalAnalysis({ alternatives: [{ name: 'X', cashflows: [-100, 50, 50] }], marr: 0.10 }).selected, null);
});

test('unequal lives: equivalent annual worth and least-common-multiple NPV agree', () => {
  const r = cf.compareUnequalLives({ alternatives: [{ name: 'A', cashflows: [-100, 70, 70] }, { name: 'B', cashflows: [-150, 70, 70, 70] }], rate: 0.10 });
  assert.equal(r.studyPeriod, 6);
  const eawA = 70 - 100 * 0.1 / (1 - 1.1 ** -2); const eawB = 70 - 150 * 0.1 / (1 - 1.1 ** -3);
  close(r.results[0].eaw, eawA, 1e-12); close(r.results[1].eaw, eawB, 1e-12);
  const pa6 = (1 - 1.1 ** -6) / 0.1;
  close(r.results[0].npvCommonMultiple, eawA * pa6, 1e-11); close(r.results[1].npvCommonMultiple, eawB * pa6, 1e-11);
  assert.equal(r.bestByEaw, 'A');
});

test('currency conversion uses cross rates and enforces base-date consistency', () => {
  const rates = { EUR: 1, USD: 1.10, GBP: 0.85 };
  close(cf.convertCurrency({ amount: 110, from: 'USD', to: 'GBP', rates }).amount, 85, 1e-12);
  assert.throws(() => cf.convertCurrency({ amount: 1, from: 'USD', to: 'GBP', rates, amountDate: '2020-01-01', rateDate: '2023-01-01' }), /differs/);
  assert.ok(cf.convertCurrency({ amount: 1, from: 'USD', to: 'GBP', rates, amountDate: '2020-01-01', rateDate: '2023-01-01', allowDateMismatch: true }).warning);
  assert.throws(() => cf.convertCurrency({ amount: 1, from: 'USD', to: 'JPY', rates }), RangeError);
});

const simple = { timeline: { constructionYears: 1, operatingYears: 3 }, discountRate: 0.10,
  capex: [{ name: 'plant', amount: 1000 }], opex: { fixed: 100 }, production: 100, price: 5 };

test('project model, no tax or inflation: FCFF = (-1000, 400, 400, 400), flows in the intended periods', () => {
  const m = cf.buildProjectCashflow(simple);
  assert.deepEqual(m.fcff, [-1000, 400, 400, 400]);
  assert.deepEqual(m.table.phase, ['construction', 'operation', 'operation', 'operation']);
  assert.deepEqual(m.table.capex, [1000, 0, 0, 0]);
  assert.deepEqual(m.table.revenue, [0, 500, 500, 500]);
  close(m.metrics.npv, -1000 + 400 * (1 - 1.1 ** -3) / 0.1, 1e-12);
  close(m.table.cumulativeDiscounted[3], m.metrics.npv, 1e-12);
  assert.equal(m.rows.length, 4); assert.equal(m.rows[2].ebitda, 400);
});

test('project model tax: straight-line depreciation with loss carry-forward (hand: tax 0, 0, 50)', () => {
  const m = cf.buildProjectCashflow({ ...simple, tax: { rate: 0.25, depreciation: { method: 'sl', life: 2 } } });
  assert.deepEqual(m.table.depreciation, [0, 500, 500, 0]);
  assert.deepEqual(m.table.taxableIncome, [0, -100, -100, 400]);
  assert.deepEqual(m.table.tax, [0, 0, 0, 50]);              // 400 less 200 of losses brought forward, at 25 %
  assert.deepEqual(m.fcff, [-1000, 400, 400, 350]);
  const noLcf = cf.buildProjectCashflow({ ...simple, tax: { rate: 0.25, lossCarryForward: false, depreciation: { method: 'sl', life: 2 } } });
  assert.deepEqual(noLcf.table.tax, [0, 0, 0, 100]);
});

test('project model: investment tax credit and allowance reduce tax as specified', () => {
  const m = cf.buildProjectCashflow({ ...simple, tax: { rate: 0.25, investmentTaxCredit: 0.1 } });
  // no depreciation: tax 100 a year; credit 100 offsets the first year
  assert.deepEqual(m.table.tax, [0, 0, 100, 100]);
  const a = cf.buildProjectCashflow({ ...simple, tax: { rate: 0.25, investmentAllowance: 0.2 } });
  assert.deepEqual(a.table.tax, [0, 50, 100, 100]);           // (400 - 200) x 25 %
});

test('project model: real and nominal modes give the same NPV with Fisher-consistent rates', () => {
  const inp = { timeline: { constructionYears: 2, operatingYears: 12 }, inflation: 0.025, discountRate: 0.09,
    capex: [{ name: 'a', amount: 1000, phasing: [0.4, 0.6], escalation: 0.01 }], opex: { fixed: 30, variablePerUnit: 0.5, escalation: 0.005 },
    production: 110, price: { type: 'indexed', price: 4 }, royaltyRate: 0.05,
    tax: { rate: 0.28, depreciation: { method: 'macrs', recoveryPeriod: 7 } }, workingCapital: { fractionOfRevenue: 0.08 },
    salvage: 50, decommissioning: { cost: 120, fund: { rate: 0.03 } } };
  const nominal = cf.buildProjectCashflow(inp);
  const real = cf.buildProjectCashflow({ ...inp, mode: 'real', discountRate: 1.09 / 1.025 - 1 });
  close(real.metrics.npv, nominal.metrics.npv, 1e-10);
  close(real.table.fcff[5] * 1.025 ** 5, nominal.table.fcff[5], 1e-12);
  // a fixed nominal price loses real value at the inflation rate
  const fixed = cf.buildProjectCashflow({ ...inp, mode: 'real', discountRate: 0.06, price: { type: 'fixed', price: 4 } });
  close(fixed.table.price[4] / fixed.table.price[3], 1 / 1.025, 1e-12);
});

test('project model reconciliation: FCFF identity, working capital nets to zero, escrow matures at the decommissioning cost', () => {
  const m = cf.buildProjectCashflow({ timeline: { constructionYears: 2, operatingYears: 10 }, inflation: 0.02, discountRate: 0.08,
    capex: [{ amount: 800, phasing: [0.5, 0.5] }, { amount: 100, year: 6, category: 'replacement' }], opex: { fixed: 40 }, production: 100, price: 4,
    tax: { rate: 0.25, depreciation: { method: 'sl', life: 5 } }, workingCapital: { fractionOfRevenue: 0.1 },
    decommissioning: { cost: 150, fund: { rate: 0.04 } }, salvage: 10 });
  const t = m.table;
  for (let k = 0; k < m.n; k++) close(t.fcff[k], t.ebitda[k] - t.tax[k] - t.capex[k] - t.deltaWC[k] - t.decomFund[k] + t.salvage[k], 1e-10, `row ${k}`);
  close(t.deltaWC.reduce((s, v) => s + v, 0), 0, 1e-10);
  const dy = m.timeline.decomYear; assert.equal(dy, 12);
  close(t.decomSpend[dy], 150 * 1.02 ** 12, 1e-12);
  let fv = 0; for (let k = 0; k < m.n; k++) fv += t.decomFund[k] * 1.04 ** (dy - k);
  close(fv, t.decomSpend[dy], 1e-10);
  assert.equal(t.decomFundBalance[dy], 0);
  assert.equal(t.fcff[dy], 0);                      // paid from the escrow, not from project cash
  close(t.capex[6], 100 * 1.02 ** 6, 1e-12);
  close(m.totals.capexByCategory.replacement, 100 * 1.02 ** 6, 1e-12);
});

test('project model without an escrow: decommissioning is a cash outflow in its year', () => {
  const m = cf.buildProjectCashflow({ ...simple, decommissioning: { cost: 200 } });
  assert.deepEqual(m.fcff, [-1000, 400, 400, 400, -200]);
  assert.equal(m.table.phase[4], 'decommissioning');
});

test('price contracts: two-way CfD pays the strike whatever the market; PPA floor and cap; merchant series', () => {
  const base = { ...simple, timeline: { constructionYears: 1, operatingYears: 3 } };
  const cfd = cf.buildProjectCashflow({ ...base, price: { type: 'cfd', strike: 6, strikeIndexation: 0, reference: [2, 6, 9] } });
  assert.deepEqual(cfd.table.revenue, [0, 600, 600, 600]);
  assert.deepEqual(cfd.table.settlement, [0, 400, 0, -300]);
  const oneWay = cf.buildProjectCashflow({ ...base, price: { type: 'cfd', strike: 6, strikeIndexation: 0, reference: [2, 6, 9], twoWay: false } });
  assert.deepEqual(oneWay.table.price, [0, 6, 6, 9]);
  const basis = cf.buildProjectCashflow({ ...base, price: { type: 'cfd', strike: 6, strikeIndexation: 0, reference: [5, 5, 5], capture: [4, 4, 4] } });
  assert.deepEqual(basis.table.price, [0, 5, 5, 5]);            // capture 4 + (6 - 5)
  const ppa = cf.buildProjectCashflow({ ...base, price: { type: 'ppa', floor: 4, cap: 7, market: [2, 6, 9], indexation: 0 } });
  assert.deepEqual(ppa.table.price, [0, 4, 6, 7]);
  const term = cf.buildProjectCashflow({ ...base, price: { type: 'ppa', price: 5, indexation: 0, market: [2, 6, 9], termYears: 2 } });
  assert.deepEqual(term.table.price, [0, 5, 5, 9]);
  const mer = cf.buildProjectCashflow({ ...base, price: { type: 'merchant', series: [3, 4, 5] } });
  assert.deepEqual(mer.table.revenue, [0, 300, 400, 500]);
});

test('sunk cost is excluded from cash flow but keeps its depreciation tax shield', () => {
  const inp = { ...simple, capex: [{ name: 'plant', amount: 1000 }, { name: 'studies already paid', amount: 300, sunk: true }],
    tax: { rate: 0.25, depreciation: { method: 'sl', life: 3 } } };
  const m = cf.buildProjectCashflow(inp);
  assert.equal(m.table.capex[0], 1000);
  close(m.table.depreciation[1], 1300 / 3, 1e-12);
  const incl = cf.buildProjectCashflow({ ...inp, excludeSunk: false });
  assert.equal(incl.table.capex[0], 1300);
  close(m.metrics.npv - incl.metrics.npv, 300, 1e-10);
});

test('royalty, start-up delay and phasing validation', () => {
  const r = cf.buildProjectCashflow({ ...simple, royaltyRate: 0.1 });
  assert.deepEqual(r.table.ebitda, [0, 350, 350, 350]);
  const d = cf.buildProjectCashflow({ ...simple, timeline: { constructionYears: 1, operatingYears: 3, startDelayYears: 1 }, delayCostPerYear: 20 });
  assert.deepEqual(d.fcff, [-1000, -20, 400, 400, 400]);
  assert.throws(() => cf.buildProjectCashflow({ ...simple, capex: [{ amount: 100, phasing: [0.5, 0.4] }] }), /sum to 1/);
  assert.throws(() => cf.buildProjectCashflow({ ...simple, discountRate: undefined }), /discountRate/);
});

test('break-even price, volume, tariff and capex by goal seek return NPV = 0', () => {
  const inp = { ...simple, inflation: 0.02, tax: { rate: 0.2, depreciation: { method: 'sl', life: 3 } } };
  const p = cf.breakEven({ inputs: inp, variable: 'price' });
  close(cf.buildProjectCashflow({ ...inp, price: { type: 'indexed', price: p.value } }).metrics.npv, 0, 1e-7);
  const tar = cf.breakEven({ inputs: inp, variable: 'tariff' });
  assert.ok(tar.value > p.value);                         // a flat nominal tariff must start higher than an indexed one
  const v = cf.breakEven({ inputs: inp, variable: 'volume' });
  close(cf.buildProjectCashflow({ ...inp, factors: { production: v.factor } }).metrics.npv, 0, 1e-7);
  const c = cf.breakEven({ inputs: inp, variable: 'capex' });
  close(cf.buildProjectCashflow({ ...inp, factors: { capex: c.factor } }).metrics.npv, 0, 1e-7);
  // untaxed, uninflated case has a closed form: p x 100 - 100 = 1000 (A/P, 10 %, 3)
  close(cf.breakEven({ inputs: simple, variable: 'price' }).value, (1000 * 0.1 / (1 - 1.1 ** -3) + 100) / 100, 1e-9);
});
