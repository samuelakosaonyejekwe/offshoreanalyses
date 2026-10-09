// Cost of capital, debt schedules, project-finance coverage ratios and debt sizing, equity cash flows, adjusted
// present value and lease-versus-buy. Rates are fractions per period; debt-schedule arrays are indexed by period
// 1..n at array index 0..n-1 unless stated.

import { requireNumber, requireSeries, factorAP, factorPA } from './tvm.js';
import { npv, irr } from './cashflow.js';

// ---------------------------------------------------------------------------------------------------------------
// Cost of capital
// ---------------------------------------------------------------------------------------------------------------

/**
 * Weighted average cost of capital. vanilla = wE kE + wD kD; postTax = wE kE + wD kD (1 − T);
 * preTax = postTax / (1 − T) (the pre-tax rate that gives the same value on pre-tax flows for a perpetuity).
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 19.
 * @param {{equity: number, debt: number, costOfEquity: number, costOfDebt: number, taxRate?: number}} o equity and
 *   debt are market values (currency) or weights; costs are fractions per year; taxRate fraction 0–1 (< 1)
 * @returns {{weightEquity: number, weightDebt: number, vanilla: number, postTax: number, preTax: number}}
 */
export function wacc({ equity, debt, costOfEquity, costOfDebt, taxRate = 0 }) {
  requireNumber('equity', equity, { min: 0 }); requireNumber('debt', debt, { min: 0 });
  requireNumber('costOfEquity', costOfEquity); requireNumber('costOfDebt', costOfDebt);
  requireNumber('taxRate', taxRate, { min: 0, lt: 1 });
  const v = equity + debt;
  if (!(v > 0)) throw new RangeError('equity + debt must be > 0');
  const we = equity / v; const wd = debt / v;
  const postTax = we * costOfEquity + wd * costOfDebt * (1 - taxRate);
  return { weightEquity: we, weightDebt: wd, vanilla: we * costOfEquity + wd * costOfDebt, postTax,
    preTax: postTax / (1 - taxRate) };
}

/**
 * CAPM cost of equity kE = rf + β (market premium) + λ (country risk premium) + size premium.
 * Ref: Sharpe (1964), J. Finance 19(3); country-risk extension after Damodaran (2003), Measuring company exposure
 * to country risk.
 * @param {{riskFree: number, beta: number, marketPremium: number, countryRiskPremium?: number, lambda?: number,
 *   sizePremium?: number}} o fractions per year; beta and lambda dimensionless (lambda = exposure to country risk)
 * @returns {number} cost of equity (fraction per year)
 */
export function capm({ riskFree, beta, marketPremium, countryRiskPremium = 0, lambda = 1, sizePremium = 0 }) {
  requireNumber('riskFree', riskFree); requireNumber('beta', beta); requireNumber('marketPremium', marketPremium);
  requireNumber('countryRiskPremium', countryRiskPremium); requireNumber('lambda', lambda);
  return riskFree + beta * marketPremium + lambda * countryRiskPremium + requireNumber('sizePremium', sizePremium);
}

/**
 * Unlever an equity beta with the Hamada equation βU = βL / [1 + (1 − T) D/E] (debt beta zero).
 * Ref: Hamada, R.S. (1972), J. Finance 27(2).
 * @param {{betaEquity: number, debtToEquity: number, taxRate?: number}} o D/E ratio (>= 0), tax rate fraction
 * @returns {number} asset (unlevered) beta
 */
export function unleverBeta({ betaEquity, debtToEquity, taxRate = 0 }) {
  requireNumber('betaEquity', betaEquity); requireNumber('debtToEquity', debtToEquity, { min: 0 });
  requireNumber('taxRate', taxRate, { min: 0, max: 1 });
  return betaEquity / (1 + (1 - taxRate) * debtToEquity);
}

/**
 * Relever an asset beta with the Hamada equation βL = βU [1 + (1 − T) D/E].
 * Ref: Hamada, R.S. (1972), J. Finance 27(2).
 * @param {{betaAsset: number, debtToEquity: number, taxRate?: number}} o D/E ratio (>= 0), tax rate fraction
 * @returns {number} equity (levered) beta
 */
export function releverBeta({ betaAsset, debtToEquity, taxRate = 0 }) {
  requireNumber('betaAsset', betaAsset); requireNumber('debtToEquity', debtToEquity, { min: 0 });
  requireNumber('taxRate', taxRate, { min: 0, max: 1 });
  return betaAsset * (1 + (1 - taxRate) * debtToEquity);
}

// ---------------------------------------------------------------------------------------------------------------
// Debt schedules
// ---------------------------------------------------------------------------------------------------------------

/**
 * Debt repayment schedule. 'annuity' = level payment, 'equalPrincipal' = constant amortisation, 'bullet' =
 * interest only with repayment at maturity. `grace` periods at the start are interest-only and count inside the
 * tenor.
 * Ref: Yescombe, E.R., Principles of Project Finance, 2nd ed.
 * @param {{type?: 'annuity' | 'equalPrincipal' | 'bullet', principal: number, rate: number, tenor: number,
 *   grace?: number}} o principal (currency), rate per period (fraction, >= 0), tenor (periods, integer),
 *   grace (periods, integer < tenor)
 * @returns {{opening: number[], interest: number[], principal: number[], payment: number[], closing: number[],
 *   totalInterest: number, totalPrincipal: number}} currency per period
 */
export function debtSchedule({ type = 'annuity', principal, rate, tenor, grace = 0 }) {
  requireNumber('principal', principal, { min: 0 }); requireNumber('rate', rate, { min: 0 });
  requireNumber('tenor', tenor, { min: 1, integer: true });
  requireNumber('grace', grace, { min: 0, max: tenor - 1, integer: true });
  const n = tenor - grace;
  const level = type === 'annuity' ? principal * factorAP(rate, n) : 0;
  const out = { opening: [], interest: [], principal: [], payment: [], closing: [], totalInterest: 0, totalPrincipal: 0 };
  let bal = principal;
  for (let k = 1; k <= tenor; k++) {
    const int = bal * rate;
    let pr = 0;
    if (k > grace) {
      if (type === 'annuity') pr = k === tenor ? bal : level - int;
      else if (type === 'equalPrincipal') pr = k === tenor ? bal : principal / n;
      else if (type === 'bullet') pr = k === tenor ? bal : 0;
      else throw new RangeError("type must be 'annuity', 'equalPrincipal' or 'bullet'");
    }
    out.opening.push(bal); out.interest.push(int); out.principal.push(pr); out.payment.push(int + pr);
    bal -= pr; out.closing.push(Math.abs(bal) < 1e-9 * Math.max(1, principal) ? 0 : bal);
    out.totalInterest += int; out.totalPrincipal += pr;
  }
  return out;
}

/**
 * Debt sculpted to cash flow available for debt service: debt service_k = CFADS_k / DSCR for a constant DSCR.
 * Give `dscr` to size the debt (principal = PV of the debt service at the debt rate) or `principal` to find the
 * constant DSCR that repays it exactly (DSCR = PV(CFADS) / principal).
 * Ref: Yescombe, Principles of Project Finance (sculpted repayment); Gatti, Project Finance in Theory and
 * Practice.
 * @param {{cfads: number[], rate: number, dscr?: number, principal?: number}} o CFADS per repayment period
 *   (currency, periods 1..n), debt rate per period (fraction), target DSCR (> 0) or principal (currency)
 * @returns {{principalAmount: number, dscr: number, opening: number[], interest: number[], principal: number[],
 *   payment: number[], closing: number[], negativeAmortisation: boolean}}
 */
export function sculptedDebt({ cfads, rate, dscr, principal }) {
  requireSeries('cfads', cfads); requireNumber('rate', rate, { min: 0 });
  let pvCfads = 0;
  for (let k = 0; k < cfads.length; k++) pvCfads += cfads[k] * Math.pow(1 + rate, -(k + 1));
  let cover = dscr; let debt = principal;
  if (dscr !== undefined) { requireNumber('dscr', dscr, { gt: 0 }); debt = pvCfads / dscr; }
  else if (principal !== undefined) { requireNumber('principal', principal, { gt: 0 }); cover = pvCfads / principal; }
  else throw new RangeError('either dscr or principal is required');
  const out = { principalAmount: debt, dscr: cover, opening: [], interest: [], principal: [], payment: [], closing: [],
    negativeAmortisation: false };
  let bal = debt;
  for (let k = 0; k < cfads.length; k++) {
    const pay = cfads[k] / cover; const int = bal * rate; const pr = pay - int;
    if (pr < -1e-9 * Math.max(1, debt)) out.negativeAmortisation = true;
    out.opening.push(bal); out.interest.push(int); out.principal.push(pr); out.payment.push(pay);
    bal -= pr; out.closing.push(Math.abs(bal) < 1e-9 * Math.max(1, debt) ? 0 : bal);
  }
  return out;
}

/**
 * Debt-service coverage ratio per period, DSCR_k = CFADS_k / (interest_k + principal_k); periods with no debt
 * service are reported as null.
 * Ref: Yescombe, Principles of Project Finance.
 * @param {number[]} cfads cash flow available for debt service per period (currency)
 * @param {number[]} debtService interest + principal per period (currency)
 * @returns {{dscr: (number | null)[], min: number, average: number, minPeriod: number}} minPeriod is 1-based
 */
export function dscr(cfads, debtService) {
  requireSeries('cfads', cfads); requireSeries('debtService', debtService);
  const out = []; let min = Infinity; let sum = 0; let cnt = 0; let minPeriod = 0;
  for (let k = 0; k < debtService.length; k++) {
    if (debtService[k] > 0) {
      const v = (k < cfads.length ? cfads[k] : 0) / debtService[k];
      out.push(v); sum += v; cnt++;
      if (v < min) { min = v; minPeriod = k + 1; }
    } else out.push(null);
  }
  return { dscr: out, min: cnt ? min : NaN, average: cnt ? sum / cnt : NaN, minPeriod };
}

/**
 * Loan-life and project-life coverage ratios at the start of each period:
 * LLCR_k = [PV of CFADS over the remaining loan life + DSRA balance] / debt outstanding, PLCR_k likewise over the
 * remaining project life, both discounted at the debt rate.
 * Ref: Yescombe, Principles of Project Finance; Gatti, Project Finance in Theory and Practice.
 * @param {{cfads: number[], opening: number[], rate: number, dsra?: number[]}} o CFADS for periods 1..N (project
 *   life), opening debt balance for periods 1..n (loan life, n <= N), debt rate (fraction), DSRA balance at the
 *   start of each period (currency, optional)
 * @returns {{llcr: (number | null)[], plcr: (number | null)[]}} one value per loan period
 */
export function coverageRatios({ cfads, opening, rate, dsra }) {
  requireSeries('cfads', cfads); requireSeries('opening', opening); requireNumber('rate', rate, { gt: -1 });
  const n = opening.length; const N = cfads.length;
  const llcr = []; const plcr = [];
  for (let k = 0; k < n; k++) {
    if (!(opening[k] > 0)) { llcr.push(null); plcr.push(null); continue; }
    let l = 0; let p = 0;
    for (let j = k; j < N; j++) {
      const v = cfads[j] * Math.pow(1 + rate, -(j - k + 1));
      if (j < n) l += v;
      p += v;
    }
    const res = dsra ? dsra[k] ?? 0 : 0;
    llcr.push((l + res) / opening[k]); plcr.push((p + res) / opening[k]);
  }
  return { llcr, plcr };
}

/**
 * Debt capacity for a target minimum DSCR and tenor. 'sculpted': PV of CFADS / DSCR. Other profiles: the largest
 * principal whose schedule keeps CFADS_k / debt service_k >= target in every period (schedules scale linearly
 * with principal). Optionally capped by a maximum gearing on total funding.
 * Ref: Yescombe, Principles of Project Finance (debt sizing).
 * @param {{cfads: number[], rate: number, tenor: number, targetDscr: number, profile?: string, grace?: number,
 *   maxGearing?: number, totalFunding?: number}} o CFADS per period (currency), rate (fraction), tenor (periods),
 *   profile 'sculpted' | 'annuity' | 'equalPrincipal' | 'bullet', gearing cap (fraction) on totalFunding (currency)
 * @returns {{debt: number, byDscr: number, byGearing: number, binding: 'dscr' | 'gearing', gearing: number | null}}
 */
export function debtCapacity({ cfads, rate, tenor, targetDscr, profile = 'sculpted', grace = 0, maxGearing, totalFunding }) {
  requireSeries('cfads', cfads); requireNumber('targetDscr', targetDscr, { gt: 0 });
  requireNumber('tenor', tenor, { min: 1, max: cfads.length, integer: true });
  const cf = Array.from(cfads).slice(0, tenor);
  let byDscr;
  if (profile === 'sculpted') byDscr = sculptedDebt({ cfads: cf, rate, dscr: targetDscr }).principalAmount;
  else {
    const unit = debtSchedule({ type: profile, principal: 1, rate, tenor, grace });
    byDscr = Infinity;
    for (let k = 0; k < tenor; k++) if (unit.payment[k] > 0) byDscr = Math.min(byDscr, cf[k] / (targetDscr * unit.payment[k]));
    byDscr = Math.max(0, byDscr);
  }
  const byGearing = maxGearing !== undefined && totalFunding !== undefined ? maxGearing * totalFunding : Infinity;
  const debt = Math.min(byDscr, byGearing);
  return { debt, byDscr, byGearing, binding: byDscr <= byGearing ? 'dscr' : 'gearing',
    gearing: totalFunding ? debt / totalFunding : null };
}

/**
 * Debt-service reserve account: target balance at the end of each period = `monthsCover`/12 × the next period's
 * debt service, and the cash movement needed to hold it (positive = funding, negative = release).
 * Ref: Yescombe, Principles of Project Finance (reserve accounts).
 * @param {{debtService: number[], monthsCover?: number, periodsPerYear?: number}} o debt service per period
 *   (currency), cover in months (default 6), periods per year (default 1)
 * @returns {{initial: number, target: number[], movement: number[]}} initial = balance required before period 1
 */
export function dsra({ debtService, monthsCover = 6, periodsPerYear = 1 }) {
  requireSeries('debtService', debtService); requireNumber('monthsCover', monthsCover, { min: 0 });
  const f = monthsCover / 12 * periodsPerYear;
  const n = debtService.length;
  const target = []; const movement = [];
  const initial = f * debtService[0];
  let prev = initial;
  for (let k = 0; k < n; k++) {
    const t = k + 1 < n ? f * debtService[k + 1] : 0;
    target.push(t); movement.push(t - prev); prev = t;
  }
  return { initial, target, movement };
}

/**
 * Financing during construction: debt drawdowns with interest during construction (IDC) capitalised, commitment
 * fee on the undrawn facility and an arrangement (up-front) fee on the facility, funded pro rata at the gearing.
 * Interest accrues on the balance at the start of each period (draws are taken at the end of the period row,
 * consistent with end-of-period cash-flow timing).
 * Ref: Yescombe, Principles of Project Finance.
 * @param {{costs: number[], gearing: number, rate: number, commitmentFeeRate?: number, arrangementFeeRate?: number}} o
 *   construction costs per period (currency, positive), gearing = debt share of total uses (0–1), rates per
 *   period (fractions); the arrangement fee is charged once in the first period
 * @returns {{debtDraw: number[], equityDraw: number[], idc: number[], fees: number[], balance: number[],
 *   facility: number, totalIdc: number, totalFees: number, totalUses: number, totalEquity: number}}
 */
export function constructionFinancing({ costs, gearing, rate, commitmentFeeRate = 0, arrangementFeeRate = 0 }) {
  requireSeries('costs', costs); requireNumber('gearing', gearing, { min: 0, lt: 1 });
  requireNumber('rate', rate, { min: 0 });
  const n = costs.length;
  let facility = 0; let res = null;
  for (let it = 0; it < 60; it++) {
    const debtDraw = []; const equityDraw = []; const idc = []; const fees = []; const balance = [];
    let bal = 0;
    for (let t = 0; t < n; t++) {
      const i = rate * bal;
      const fee = (t === 0 ? arrangementFeeRate * facility : 0) + commitmentFeeRate * Math.max(0, facility - bal) * (t > 0 ? 1 : 0);
      const uses = costs[t] + i + fee;
      const d = gearing * uses;
      debtDraw.push(d); equityDraw.push(uses - d); idc.push(i); fees.push(fee);
      bal += d; balance.push(bal);
    }
    res = { debtDraw, equityDraw, idc, fees, balance };
    if (Math.abs(bal - facility) <= 1e-12 * Math.max(1, bal)) { facility = bal; break; }
    facility = bal;
  }
  const sum = (a) => a.reduce((s, v) => s + v, 0);
  return { ...res, facility, totalIdc: sum(res.idc), totalFees: sum(res.fees),
    totalUses: sum(res.debtDraw) + sum(res.equityDraw), totalEquity: sum(res.equityDraw) };
}

/**
 * Integrated project-finance model on top of an unlevered after-tax cash flow (FCFF): construction funding with
 * IDC and fees, term debt (given gearing, or sized to a target DSCR up to a gearing cap), coverage ratios, DSRA,
 * distribution lock-up and equity cash flows.
 * Assumptions stated plainly: CFADS is the FCFF of the operating rows (tax computed without interest); the
 * interest tax shield taxRate × interest is added to the equity flow in the period it arises, i.e. taxable
 * capacity is assumed. Rows before `firstOperatingRow` (at least row 0) are funding rows; debt service starts in
 * the next row.
 * Ref: Yescombe, Principles of Project Finance; Gatti, Project Finance in Theory and Practice.
 * @param {object} o
 * @param {number[]} o.fcff unlevered free cash flow per row (currency), row 0 = time zero
 * @param {number} [o.firstOperatingRow] default: first row with positive FCFF
 * @param {{gearing?: number, targetDscr?: number, maxGearing?: number, rate: number, tenor: number,
 *   type?: 'annuity' | 'equalPrincipal' | 'bullet' | 'sculpted', grace?: number, arrangementFeeRate?: number,
 *   commitmentFeeRate?: number, dsraMonths?: number, lockupDscr?: number}} o.debt give `gearing` (0–1) or
 *   `targetDscr`; rate per year (fraction), tenor years from the first operating row
 * @param {number} [o.taxRate] fraction used for the interest tax shield (default 0)
 * @param {number} [o.equityDiscountRate] cost of equity for equity NPV (fraction, default debt.rate)
 * @returns {{gearing: number, debt: number, equity: number, idc: number, fees: number, binding: string,
 *   cfads: number[], interest: number[], principal: number[], debtService: number[], balance: number[],
 *   dscr: (number | null)[], minDscr: number, averageDscr: number, llcr: number | null, plcr: number | null,
 *   dsraBalance: number[], lockedUp: boolean[], equityCashflows: number[], equityIrr: number | null,
 *   equityIrrStatus: string, equityNpv: number, projectIrr: number | null}} all series are per model row
 */
export function financeProject({ fcff, firstOperatingRow, debt, taxRate = 0, equityDiscountRate }) {
  requireSeries('fcff', fcff, { minLength: 2 });
  if (!debt) throw new RangeError('debt terms are required');
  requireNumber('debt.rate', debt.rate, { min: 0 }); requireNumber('debt.tenor', debt.tenor, { min: 1, integer: true });
  requireNumber('taxRate', taxRate, { min: 0, max: 1 });
  const N = fcff.length;
  let s = firstOperatingRow ?? fcff.findIndex((v) => v > 0);
  if (s < 0) throw new RangeError('fcff has no positive (operating) row');
  s = Math.max(1, s);
  if (s >= N) throw new RangeError('firstOperatingRow is beyond the cash-flow horizon');
  const tenor = Math.min(debt.tenor, N - s);
  const type = debt.type ?? 'annuity';
  const costs = []; for (let t = 0; t < s; t++) costs.push(Math.max(0, -fcff[t]));
  const cfads = Array.from(fcff, (v, t) => (t >= s ? v : 0));
  const opCfads = cfads.slice(s);
  const months = debt.dsraMonths ?? 0;

  const evaluate = (g) => {
    const cons = constructionFinancing({ costs, gearing: g, rate: debt.rate,
      commitmentFeeRate: debt.commitmentFeeRate ?? 0, arrangementFeeRate: debt.arrangementFeeRate ?? 0 });
    const B = cons.facility;
    const sch = B > 0
      ? (type === 'sculpted'
        ? sculptedDebt({ cfads: opCfads.slice(0, tenor), rate: debt.rate, principal: B })
        : debtSchedule({ type, principal: B, rate: debt.rate, tenor, grace: debt.grace ?? 0 }))
      : { opening: new Array(tenor).fill(0), interest: new Array(tenor).fill(0), principal: new Array(tenor).fill(0),
        payment: new Array(tenor).fill(0), closing: new Array(tenor).fill(0) };
    return { cons, B, sch };
  };
  const capacityGap = (g) => {
    const { B, cons } = evaluate(g);
    const cap = debtCapacity({ cfads: opCfads, rate: debt.rate, tenor, targetDscr: debt.targetDscr,
      profile: type, grace: debt.grace ?? 0 }).debt;
    return { gap: B - cap, cons };
  };

  let g; let binding = 'gearing';
  if (debt.gearing !== undefined) g = requireNumber('debt.gearing', debt.gearing, { min: 0, lt: 1 });
  else if (debt.targetDscr !== undefined) {
    const gMax = requireNumber('debt.maxGearing', debt.maxGearing ?? 0.9, { min: 0, lt: 1 });
    if (capacityGap(gMax).gap <= 0) g = gMax;
    else {
      let lo = 0; let hi = gMax;
      for (let it = 0; it < 100; it++) { const mid = 0.5 * (lo + hi); if (capacityGap(mid).gap > 0) hi = mid; else lo = mid; }
      g = 0.5 * (lo + hi); binding = 'dscr';
    }
  } else throw new RangeError('debt.gearing or debt.targetDscr is required');

  const { cons, B, sch } = evaluate(g);
  const interest = new Array(N).fill(0); const principal = new Array(N).fill(0); const debtService = new Array(N).fill(0);
  const balance = new Array(N).fill(0); const dsraBalance = new Array(N).fill(0); const lockedUp = new Array(N).fill(false);
  for (let t = 0; t < s; t++) balance[t] = cons.balance[t];
  for (let k = 0; k < tenor; k++) {
    interest[s + k] = sch.interest[k]; principal[s + k] = sch.principal[k]; debtService[s + k] = sch.payment[k];
    balance[s + k] = sch.closing[k];
  }
  const cov = dscr(opCfads, sch.payment);
  const res = B > 0 ? dsra({ debtService: sch.payment, monthsCover: months }) : { initial: 0, target: [], movement: [] };
  const ratios = B > 0 ? coverageRatios({ cfads: opCfads, opening: sch.opening, rate: debt.rate }) : { llcr: [null], plcr: [null] };
  const equityCashflows = new Array(N).fill(0);
  for (let t = 0; t < s; t++) equityCashflows[t] = -cons.equityDraw[t] + (fcff[t] > 0 ? fcff[t] : 0);
  equityCashflows[s - 1] -= res.initial;
  dsraBalance[s - 1] = res.initial;
  let trapped = 0;
  const dscrRow = new Array(N).fill(null);
  for (let t = s; t < N; t++) {
    const k = t - s;
    const move = k < tenor ? res.movement[k] ?? 0 : 0;
    dsraBalance[t] = k < tenor ? res.target[k] ?? 0 : 0;
    let cash = cfads[t] + taxRate * interest[t] - debtService[t] - move;
    if (k < tenor) dscrRow[t] = cov.dscr[k];
    if (debt.lockupDscr !== undefined && k < tenor && dscrRow[t] !== null && dscrRow[t] < debt.lockupDscr && cash > 0) {
      trapped += cash; cash = 0; lockedUp[t] = true;
    } else if (trapped > 0) { cash += trapped; trapped = 0; }
    equityCashflows[t] = cash;
  }
  if (trapped > 0) equityCashflows[N - 1] += trapped;
  const eIrr = irr(equityCashflows);
  const pIrr = irr(Array.from(fcff));
  return { gearing: g, debt: B, equity: cons.totalEquity, idc: cons.totalIdc, fees: cons.totalFees, binding,
    cfads, interest, principal, debtService, balance, dscr: dscrRow, minDscr: cov.min, averageDscr: cov.average,
    llcr: ratios.llcr[0], plcr: ratios.plcr[0], dsraBalance, lockedUp, equityCashflows, equityIrr: eIrr.irr,
    equityIrrStatus: eIrr.status, equityNpv: npv(equityCashflows, equityDiscountRate ?? debt.rate),
    projectIrr: pIrr.irr };
}

/**
 * Refinancing hook: replace the remaining debt at a chosen period by a new loan (new rate, tenor and profile,
 * with an up-front fee) and value the change in debt service.
 * Ref: Yescombe, Principles of Project Finance (refinancing).
 * @param {{balance: number, remainingService: number[], newRate: number, newTenor: number, newType?: string,
 *   feeRate?: number, discountRate: number}} o outstanding balance (currency), remaining debt service of the old
 *   loan per period, new terms, fee as a fraction of the balance, discount rate for the gain (fraction)
 * @returns {{schedule: object, fee: number, pvOldService: number, pvNewService: number, gain: number}} gain > 0
 *   means the refinancing adds value at the discount rate
 */
export function refinance({ balance, remainingService, newRate, newTenor, newType = 'annuity', feeRate = 0, discountRate }) {
  requireNumber('balance', balance, { gt: 0 }); requireSeries('remainingService', remainingService);
  requireNumber('discountRate', discountRate, { gt: -1 });
  const schedule = debtSchedule({ type: newType, principal: balance, rate: newRate, tenor: newTenor });
  const fee = feeRate * balance;
  const pvOld = npv([0, ...remainingService], discountRate);
  const pvNew = npv([0, ...schedule.payment], discountRate);
  return { schedule, fee, pvOldService: pvOld, pvNewService: pvNew, gain: pvOld - pvNew - fee };
}

/**
 * Equity cash flows from an unlevered cash flow and a debt schedule:
 * equity_t = FCFF_t + drawdown_t − interest_t (1 − T) − principal_t; returns equity and project IRR/NPV.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 19 (flow-to-equity method).
 * @param {{fcff: number[], drawdown: number[], interest: number[], principal: number[], taxRate?: number,
 *   equityRate: number, projectRate: number}} o per-row series (currency) and discount rates (fractions)
 * @returns {{equityCashflows: number[], equityIrr: number | null, equityNpv: number, projectIrr: number | null,
 *   projectNpv: number, equityIrrStatus: string}}
 */
export function equityCashflows({ fcff, drawdown, interest, principal, taxRate = 0, equityRate, projectRate }) {
  requireSeries('fcff', fcff, { minLength: 2 });
  const at = (a, t) => (a && t < a.length ? a[t] : 0);
  const eq = Array.from(fcff, (v, t) => v + at(drawdown, t) - at(interest, t) * (1 - taxRate) - at(principal, t));
  const e = irr(eq); const p = irr(Array.from(fcff));
  return { equityCashflows: eq, equityIrr: e.irr, equityIrrStatus: e.status, equityNpv: npv(eq, equityRate),
    projectIrr: p.irr, projectNpv: npv(fcff, projectRate) };
}

/**
 * Adjusted present value: base-case NPV of the unlevered flows at the unlevered cost of equity, plus the PV of
 * interest tax shields at the cost of debt, less issue costs, plus other financing side effects.
 * Ref: Myers, S.C. (1974), J. Finance 29(1); Brealey, Myers & Allen ch. 19.
 * @param {{fcff: number[], unleveredRate: number, interest: number[], taxRate: number, debtRate: number,
 *   issueCosts?: number, otherSideEffects?: number}} o per-row series (currency), fractions, currency
 * @returns {{baseNpv: number, pvTaxShield: number, apv: number}} currency
 */
export function adjustedPresentValue({ fcff, unleveredRate, interest, taxRate, debtRate, issueCosts = 0, otherSideEffects = 0 }) {
  requireNumber('taxRate', taxRate, { min: 0, max: 1 });
  const baseNpv = npv(fcff, unleveredRate);
  const pvTaxShield = npv(Array.from(requireSeries('interest', interest), (v) => v * taxRate), debtRate);
  return { baseNpv, pvTaxShield, apv: baseNpv + pvTaxShield - issueCosts + otherSideEffects };
}

/**
 * Lease versus buy: net advantage to leasing NAL = purchase price − PV[after-tax lease payments + lost
 * depreciation tax shields] − PV[after-tax salvage forgone] + PV[after-tax owner-only costs avoided], all
 * discounted at the after-tax cost of debt. Lease when NAL > 0.
 * Ref: Myers, Dill & Bautista (1976), J. Finance 31(3); Brealey, Myers & Allen ch. 25.
 * @param {{price: number, leasePayments: number[], inAdvance?: boolean, depreciation?: number[], taxRate?: number,
 *   debtRate: number, salvage?: number, ownerCosts?: number[]}} o price (currency), lease payments per period,
 *   payments at period start when inAdvance (default true), depreciation charges of ownership per period
 *   (periods 1..n), pre-tax cost of debt (fraction), after-tax salvage value at the end (currency), owner-only
 *   operating costs per period (pre-tax)
 * @returns {{nal: number, pvLease: number, pvDepreciationShield: number, pvSalvage: number, pvOwnerCosts: number,
 *   discountRate: number, decision: 'lease' | 'buy'}}
 */
export function leaseVsBuy({ price, leasePayments, inAdvance = true, depreciation = [], taxRate = 0, debtRate,
  salvage = 0, ownerCosts = [] }) {
  requireNumber('price', price, { min: 0 }); requireSeries('leasePayments', leasePayments);
  requireNumber('taxRate', taxRate, { min: 0, max: 1 }); requireNumber('debtRate', debtRate, { gt: -1 });
  const r = debtRate * (1 - taxRate);
  const n = leasePayments.length;
  let pvLease = 0; let pvDep = 0; let pvOwn = 0;
  for (let k = 0; k < n; k++) pvLease += leasePayments[k] * (1 - taxRate) * Math.pow(1 + r, -(inAdvance ? k : k + 1));
  for (let k = 0; k < depreciation.length; k++) pvDep += taxRate * depreciation[k] * Math.pow(1 + r, -(k + 1));
  for (let k = 0; k < ownerCosts.length; k++) pvOwn += ownerCosts[k] * (1 - taxRate) * Math.pow(1 + r, -(k + 1));
  const pvSalvage = salvage * Math.pow(1 + r, -n);
  const nal = price - pvLease - pvDep - pvSalvage + pvOwn;
  return { nal, pvLease, pvDepreciationShield: pvDep, pvSalvage, pvOwnerCosts: pvOwn, discountRate: r,
    decision: nal > 0 ? 'lease' : 'buy' };
}

/**
 * Inflation-linked (capital-indexed) annuity debt: a real annuity at the real rate whose every cash flow and
 * balance is scaled by the cumulative inflation index.
 * Ref: Yescombe, Principles of Project Finance (inflation-indexed financing).
 * @param {{principal: number, realRate: number, tenor: number, inflation: number | number[]}} o principal
 *   (currency), real rate per period (fraction), tenor (periods), inflation per period (fraction or series)
 * @returns {{index: number[], interest: number[], principal: number[], payment: number[], closing: number[],
 *   indexation: number[]}} nominal currency per period; indexation = uplift of the opening balance in the period
 */
export function inflationLinkedDebt({ principal, realRate, tenor, inflation }) {
  const real = debtSchedule({ type: 'annuity', principal, rate: realRate, tenor });
  const out = { index: [], interest: [], principal: [], payment: [], closing: [], indexation: [] };
  let idx = 1;
  for (let k = 0; k < tenor; k++) {
    const f = typeof inflation === 'number' ? inflation : inflation[Math.min(k, inflation.length - 1)];
    requireNumber('inflation', f, { gt: -1 });
    out.indexation.push(real.opening[k] * idx * f);
    idx *= 1 + f;
    out.index.push(idx); out.interest.push(real.interest[k] * idx); out.principal.push(real.principal[k] * idx);
    out.payment.push(real.payment[k] * idx); out.closing.push(real.closing[k] * idx);
  }
  return out;
}

/**
 * Level annual payment of an annuity loan, and the PV of an annuity per unit payment (convenience re-exports of
 * the capital-recovery and present-worth factors for finance callers).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} principal currency
 * @param {number} rate fraction per period
 * @param {number} tenor periods
 * @returns {{payment: number, pvPerUnit: number}}
 */
export function annuityLoanPayment(principal, rate, tenor) {
  return { payment: requireNumber('principal', principal) * factorAP(rate, tenor), pvPerUnit: factorPA(rate, tenor) };
}
