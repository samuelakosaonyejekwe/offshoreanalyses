// Discounted cash-flow measures (NPV, IRR with multiple-root diagnostics, MIRR, PI, B/C, ROI, payback), break-even
// goal seek, incremental analysis of mutually exclusive alternatives and the project cash-flow model builder.
// Conventions: cashflows[k] is the net flow of period k (inflow positive), cashflows[0] at time zero.

import { brent } from '../math/roots.js';
import { requireNumber, requireSeries, discountFactors, factorAP } from './tvm.js';
import { depreciationSchedule } from './depreciation.js';

// ---------------------------------------------------------------------------------------------------------------
// Present-value measures
// ---------------------------------------------------------------------------------------------------------------

/**
 * Net present value NPV = Σ CF_k × DF_k. End-of-period discounting by default; with midYear the flows of
 * periods k >= 1 are discounted at k - 0.5. `rate` may be a vector (rate[k] applies between k-1 and k).
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 2 and 5.
 * @param {number[]} cashflows currency per period
 * @param {number | number[]} rate discount rate per period (fraction, > -1) or one per period
 * @param {{midYear?: boolean}} [o]
 * @returns {number} currency at time zero
 */
export function npv(cashflows, rate, o = {}) {
  requireSeries('cashflows', cashflows);
  const df = discountFactors(rate, cashflows.length, o);
  let s = 0;
  for (let k = 0; k < cashflows.length; k++) s += cashflows[k] * df[k];
  return s;
}

/**
 * Discounted and cumulative discounted series of a cash flow.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 5.
 * @param {number[]} cashflows currency per period
 * @param {number | number[]} rate discount rate (fraction) or per-period vector
 * @param {{midYear?: boolean}} [o]
 * @returns {{factors: number[], discounted: number[], cumulative: number[], cumulativeDiscounted: number[]}}
 */
export function discountSeries(cashflows, rate, o = {}) {
  requireSeries('cashflows', cashflows);
  const df = discountFactors(rate, cashflows.length, o);
  const discounted = []; const cumulative = []; const cumulativeDiscounted = [];
  let c = 0; let cd = 0;
  for (let k = 0; k < cashflows.length; k++) {
    c += cashflows[k]; cd += cashflows[k] * df[k];
    discounted.push(cashflows[k] * df[k]); cumulative.push(c); cumulativeDiscounted.push(cd);
  }
  return { factors: Array.from(df), discounted, cumulative, cumulativeDiscounted };
}

/**
 * NPV profile: NPV evaluated on a list of discount rates (for plotting and visual IRR checks).
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 5.
 * @param {number[]} cashflows currency per period
 * @param {number[]} rates discount rates (fractions, > -1)
 * @returns {{rates: number[], npv: number[]}}
 */
export function npvProfile(cashflows, rates) {
  requireSeries('rates', rates);
  return { rates: Array.from(rates), npv: Array.from(rates, (r) => npv(cashflows, r)) };
}

function signChanges(a) {
  let n = 0; let prev = 0;
  for (const v of a) {
    const s = Math.sign(v);
    if (s === 0) continue;
    if (prev !== 0 && s !== prev) n++;
    prev = s;
  }
  return n;
}

// Sign-preserving evaluation of Σ CF_k x^k for x in (0, ∞), mapped from u in (0, 1) with x = u / (1 - u).
// For x > 1 the polynomial is divided by x^n so that long series do not overflow.
function pvPoly(cf, u) {
  const n = cf.length - 1;
  const x = u / (1 - u);
  let s = 0;
  if (x <= 1) { for (let k = n; k >= 0; k--) s = s * x + cf[k]; return s; }
  const y = 1 / x;
  for (let k = 0; k <= n; k++) s = s * y + cf[k];
  return s;
}

/**
 * Internal rate(s) of return: every real root r in (-1, ∞) of Σ CF_k (1 + r)^-k = 0, found by sign-change
 * bracketing on a dense grid plus tangency (double-root) detection, each refined with Brent's method.
 * Diagnostics: Descartes' rule (number of sign changes in the flows bounds the number of roots r > -1) and
 * Norstrom's criterion (cumulative flow series starting negative with exactly one sign change and a non-zero
 * final value guarantees a unique positive IRR).
 * Roots closer than about 1e-3 apart in rate may be reported as one.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 7 (multiple rates of return, Descartes and Norstrom tests);
 * Norstrom, C.J. (1972), J. Financial and Quantitative Analysis 7(3).
 * @param {number[]} cashflows currency per period (at least two entries)
 * @param {{minRate?: number, maxRate?: number, grid?: number}} [o] search window in rate (defaults -0.9999, 9999)
 *   and number of grid intervals (default 4000)
 * @returns {{status: 'unique' | 'multiple' | 'none', irr: number | null, roots: number[], descartesSignChanges: number,
 *   norstromSignChanges: number, norstromUnique: boolean, message: string}} irr is set only when status is 'unique';
 *   rates are fractions per period
 */
export function irr(cashflows, o = {}) {
  requireSeries('cashflows', cashflows, { minLength: 2 });
  const cf = Array.from(cashflows);
  const descartes = signChanges(cf);
  const cum = []; let c = 0;
  for (const v of cf) { c += v; cum.push(c); }
  const norstromSignChanges = signChanges(cum);
  const firstNonZero = cum.find((v) => v !== 0) ?? 0;
  const norstromUnique = firstNonZero < 0 && norstromSignChanges === 1 && cum[cum.length - 1] !== 0;
  const roots = [];
  if (descartes > 0) {
    const minRate = o.minRate ?? -0.9999; const maxRate = o.maxRate ?? 9999;
    requireNumber('minRate', minRate, { gt: -1 }); requireNumber('maxRate', maxRate, { gt: minRate });
    const N = o.grid ?? 4000;
    const uLo = 1 / (2 + maxRate); const uHi = 1 / (2 + minRate);   // u = x/(1+x), x = 1/(1+r)
    const f = (u) => pvPoly(cf, u);
    let scale = 0;
    for (const v of cf) scale = Math.max(scale, Math.abs(v));
    const toRate = (u) => (1 - u) / u - 1;
    const push = (u) => {
      const r = toRate(u);
      if (!roots.some((q) => Math.abs(q - r) <= 1e-9 * (1 + Math.abs(r)))) roots.push(r);
    };
    let u0 = uLo; let f0 = f(u0); let fm1 = NaN;
    for (let k = 1; k <= N; k++) {
      const u1 = uLo + (uHi - uLo) * k / N;
      const f1 = f(u1);
      if (f0 === 0) push(u0);
      else if (f0 * f1 < 0) push(brent(f, u0, u1, { tol: 1e-15, maxIter: 300 }));
      else if (k >= 2 && Number.isFinite(fm1) && f0 * fm1 > 0 && f0 * f1 > 0
        && Math.abs(f0) < Math.abs(fm1) && Math.abs(f0) < Math.abs(f1)) {
        // local extremum of |f| without a sign change: test for a tangent (double) root by golden search
        let a = u0 - (uHi - uLo) / N; let b = u1;
        const g = (u) => Math.abs(f(u));
        const gr = 0.6180339887498949;
        let x1 = b - gr * (b - a); let x2 = a + gr * (b - a); let g1 = g(x1); let g2 = g(x2);
        for (let it = 0; it < 80; it++) {
          if (g1 < g2) { b = x2; x2 = x1; g2 = g1; x1 = b - gr * (b - a); g1 = g(x1); }
          else { a = x1; x1 = x2; g1 = g2; x2 = a + gr * (b - a); g2 = g(x2); }
        }
        if (Math.min(g1, g2) <= 1e-10 * scale) push(0.5 * (a + b));
      }
      fm1 = f0; u0 = u1; f0 = f1;
    }
    if (f0 === 0) push(u0);
    roots.sort((p, q) => p - q);
  }
  const status = roots.length === 0 ? 'none' : roots.length === 1 ? 'unique' : 'multiple';
  const message = status === 'unique' ? 'Single internal rate of return.'
    : status === 'multiple'
      ? `Non-conventional cash flow: ${roots.length} internal rates of return; use NPV or MIRR for the decision.`
      : 'No internal rate of return in (-100 %, ∞): the flows never change sign or NPV does not cross zero.';
  return { status, irr: status === 'unique' ? roots[0] : null, roots, descartesSignChanges: descartes,
    norstromSignChanges, norstromUnique, message };
}

/**
 * Modified internal rate of return: negative flows are discounted to time 0 at the finance rate, positive flows
 * compounded to period n at the reinvestment rate; MIRR = (FV_pos / -PV_neg)^(1/n) - 1.
 * Ref: Lin, S.A.Y. (1976), The modified internal rate of return and investment criterion, Eng. Economist 21(4);
 * Brealey, Myers & Allen ch. 5.
 * @param {number[]} cashflows currency per period
 * @param {number} financeRate rate applied to outflows (fraction, > -1)
 * @param {number} reinvestRate rate applied to inflows (fraction, > -1)
 * @returns {number} MIRR (fraction per period); NaN when there are no outflows or no inflows
 */
export function mirr(cashflows, financeRate, reinvestRate) {
  requireSeries('cashflows', cashflows, { minLength: 2 });
  requireNumber('financeRate', financeRate, { gt: -1 }); requireNumber('reinvestRate', reinvestRate, { gt: -1 });
  const n = cashflows.length - 1;
  let pvNeg = 0; let fvPos = 0;
  for (let k = 0; k <= n; k++) {
    if (cashflows[k] < 0) pvNeg += cashflows[k] * Math.pow(1 + financeRate, -k);
    else fvPos += cashflows[k] * Math.pow(1 + reinvestRate, n - k);
  }
  if (pvNeg === 0 || fvPos === 0) return NaN;
  return Math.pow(fvPos / -pvNeg, 1 / n) - 1;
}

// Splits a net flow into investment (leading negative flows before the first positive flow) and the rest.
function splitInvestment(cashflows, investment) {
  const n = cashflows.length;
  const inv = new Array(n).fill(0);
  if (investment) {
    requireSeries('investment', investment);
    for (let k = 0; k < Math.min(n, investment.length); k++) inv[k] = investment[k];
  } else {
    for (let k = 0; k < n && cashflows[k] <= 0; k++) inv[k] = -cashflows[k];
  }
  return inv;
}

/**
 * Profitability index PI = PV(net flows excluding investment) / PV(investment) = 1 + NPV / PV(investment).
 * When `investment` is not given, the leading run of non-positive flows is taken as the investment.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 5 (capital rationing).
 * @param {number[]} cashflows net currency per period
 * @param {number} rate discount rate (fraction)
 * @param {{investment?: number[], midYear?: boolean}} [o] investment = positive capital outlays per period
 * @returns {{pi: number, pvInvestment: number, pvReturns: number, npv: number}}
 */
export function profitabilityIndex(cashflows, rate, o = {}) {
  const inv = splitInvestment(cashflows, o.investment);
  const pvInvestment = npv(inv, rate, o);
  const total = npv(cashflows, rate, o);
  const pvReturns = total + pvInvestment;
  return { pi: pvInvestment > 0 ? pvReturns / pvInvestment : NaN, pvInvestment, pvReturns, npv: total };
}

/**
 * Benefit–cost ratios. Conventional: PV(B) / [PV(I) - PV(S) + PV(O&M)]. Modified: [PV(B) - PV(O&M)] / [PV(I) - PV(S)].
 * Ref: Blank & Tarquin, Engineering Economy, ch. 9.
 * @param {{benefits: number[], investment: number[], om?: number[], salvage?: number[], rate: number}} o all series
 *   are positive magnitudes per period (currency), rate is the discount rate (fraction)
 * @returns {{conventional: number, modified: number, pvBenefits: number, pvInvestment: number, pvOm: number,
 *   pvSalvage: number, netBenefit: number}}
 */
export function benefitCostRatio({ benefits, investment, om = [0], salvage = [0], rate }) {
  const pvBenefits = npv(benefits, rate); const pvInvestment = npv(investment, rate);
  const pvOm = npv(om, rate); const pvSalvage = npv(salvage, rate);
  const netInv = pvInvestment - pvSalvage;
  return { conventional: pvBenefits / (netInv + pvOm), modified: (pvBenefits - pvOm) / netInv, pvBenefits,
    pvInvestment, pvOm, pvSalvage, netBenefit: pvBenefits - pvOm - netInv };
}

/**
 * Return on investment (undiscounted): ROI = (total net return - investment) / investment, and its simple
 * annual average over the number of periods after time zero.
 * Ref: standard accounting definition; Blank & Tarquin, Engineering Economy, ch. 7.
 * @param {number[]} cashflows net currency per period
 * @param {{investment?: number[]}} [o] capital outlays per period (positive); default leading non-positive flows
 * @returns {{roi: number, annualised: number, investment: number, netReturn: number}} fractions and currency
 */
export function roi(cashflows, o = {}) {
  requireSeries('cashflows', cashflows, { minLength: 2 });
  const inv = splitInvestment(cashflows, o.investment);
  let investment = 0; let total = 0;
  for (let k = 0; k < cashflows.length; k++) { investment += inv[k]; total += cashflows[k]; }
  const r = investment > 0 ? total / investment : NaN;
  return { roi: r, annualised: r / (cashflows.length - 1), investment, netReturn: total + investment };
}

/**
 * Simple and discounted payback periods: the first time the (discounted) cumulative flow reaches zero, with
 * linear interpolation inside the period in which it turns non-negative. null when never recovered.
 * A later return to negative (e.g. decommissioning) is flagged by `staysPositive`.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 13.
 * @param {number[]} cashflows currency per period
 * @param {number} [rate] discount rate for the discounted payback (fraction, default 0)
 * @param {{midYear?: boolean}} [o]
 * @returns {{simple: number | null, discounted: number | null, staysPositive: boolean}} periods from time zero
 */
export function payback(cashflows, rate = 0, o = {}) {
  requireSeries('cashflows', cashflows);
  const df = discountFactors(rate, cashflows.length, o);
  const first = (useDf) => {
    let c = 0;
    for (let k = 0; k < cashflows.length; k++) {
      const v = cashflows[k] * (useDf ? df[k] : 1);
      const prev = c; c += v;
      if (prev < 0 && c >= 0) return k - 1 + (-prev) / v;
      if (k === 0 && c >= 0 && v > 0) return 0;
    }
    return null;
  };
  let total = 0; let staysPositive = true; let crossed = false;
  for (let k = 0; k < cashflows.length; k++) {
    total += cashflows[k];
    if (total >= 0) crossed = true; else if (crossed) staysPositive = false;
  }
  return { simple: first(false), discounted: first(true), staysPositive: crossed && staysPositive };
}

/**
 * The standard metric set for one cash-flow series.
 * Ref: definitions as on npv, irr, mirr, profitabilityIndex, payback, roi.
 * @param {number[]} cashflows currency per period
 * @param {number | number[]} rate discount rate (fraction)
 * @param {{financeRate?: number, reinvestRate?: number, midYear?: boolean, investment?: number[]}} [o] MIRR rates
 *   default to the (first) discount rate
 * @returns {{npv: number, irr: number | null, irrStatus: string, irrRoots: number[], irrDiagnostics: object,
 *   mirr: number, pi: number, roi: number, paybackSimple: number | null, paybackDiscounted: number | null,
 *   eaw: number}} eaw = equivalent annual worth (currency per period)
 */
export function cashflowMetrics(cashflows, rate, o = {}) {
  const r0 = typeof rate === 'number' ? rate : rate[Math.min(1, rate.length - 1)];
  const v = npv(cashflows, rate, o);
  const ir = irr(cashflows);
  const pb = payback(cashflows, rate, o);
  const n = cashflows.length - 1;
  return { npv: v, irr: ir.irr, irrStatus: ir.status, irrRoots: ir.roots, irrDiagnostics: ir,
    mirr: mirr(cashflows, o.financeRate ?? r0, o.reinvestRate ?? r0),
    pi: profitabilityIndex(cashflows, rate, o).pi, roi: roi(cashflows, o).roi,
    paybackSimple: pb.simple, paybackDiscounted: pb.discounted,
    eaw: n > 0 ? v * factorAP(r0, n) : v };
}

// ---------------------------------------------------------------------------------------------------------------
// Goal seek and break-even
// ---------------------------------------------------------------------------------------------------------------

/**
 * Solve f(x) = target for a scalar x by expanding the bracket [lo, hi] geometrically (within [min, max]) until
 * the sign changes, then Brent's method. Throws when no sign change is found.
 * Ref: Brent, R.P. (1973), Algorithms for Minimization without Derivatives, ch. 4.
 * @param {(x: number) => number} f function of the unknown
 * @param {{target?: number, lo: number, hi: number, min?: number, max?: number, tol?: number}} o bracket guess and limits
 * @returns {number} x such that f(x) = target
 */
export function solveFor(f, { target = 0, lo, hi, min = -Infinity, max = Infinity, tol = 1e-12 }) {
  requireNumber('lo', lo); requireNumber('hi', hi, { gt: lo });
  const g = (x) => f(x) - target;
  let a = lo; let b = hi; let ga = g(a); let gb = g(b);
  for (let it = 0; it < 80 && ga * gb > 0; it++) {
    const w = b - a;
    if (Math.abs(ga) < Math.abs(gb)) { a = Math.max(min, a - w); ga = g(a); }
    else { b = Math.min(max, b + w); gb = g(b); }
    if (a === min && b === max) break;
  }
  if (ga === 0) return a;
  if (gb === 0) return b;
  if (!(ga * gb < 0)) throw new RangeError('solveFor: target not reachable inside the search limits');
  return brent(g, a, b, { tol, maxIter: 300 });
}

/**
 * Break-even of the project cash-flow model: the value of one driver that makes a metric equal a target
 * (default NPV = 0).
 * variable: 'price' — constant real (inflation-indexed) unit price replacing the price contract;
 * 'tariff' — flat nominal unit price; 'volume' — multiplier on the production schedule;
 * 'capex' — multiplier on all capex; 'opex' — multiplier on all opex.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 13 (break-even analysis).
 * @param {{inputs: object, variable: 'price' | 'tariff' | 'volume' | 'capex' | 'opex', metric?: string,
 *   target?: number, lo?: number, hi?: number}} o inputs as for buildProjectCashflow; metric is a key of
 *   model.metrics (default 'npv')
 * @returns {{variable: string, value: number, factor: number | null, metric: string, target: number}} value is the
 *   break-even price (currency per unit) for 'price'/'tariff', otherwise the multiplier; for 'volume' and 'capex'
 *   multiply the base schedule/estimate by `factor`
 */
export function breakEven({ inputs, variable, metric = 'npv', target = 0, lo, hi }) {
  const run = (x) => {
    const inp = { ...inputs, factors: { ...(inputs.factors || {}) } };
    if (variable === 'price') inp.price = { type: 'indexed', price: x };
    else if (variable === 'tariff') inp.price = { type: 'fixed', price: x };
    else if (variable === 'volume') inp.factors.production = x;
    else if (variable === 'capex') inp.factors.capex = x;
    else if (variable === 'opex') inp.factors.opex = x;
    else throw new RangeError("variable must be 'price', 'tariff', 'volume', 'capex' or 'opex'");
    if (variable === 'price' || variable === 'tariff') inp.factors.price = 1;
    const m = buildProjectCashflow(inp).metrics[metric];
    if (typeof m !== 'number' || !Number.isFinite(m)) throw new RangeError(`metric '${metric}' is not available`);
    return m;
  };
  const isPrice = variable === 'price' || variable === 'tariff';
  const x = solveFor(run, { target, lo: lo ?? (isPrice ? 0 : 0.01), hi: hi ?? (isPrice ? 100 : 2), min: 0, max: 1e12 });
  return { variable, value: x, factor: isPrice ? null : x, metric, target };
}

// ---------------------------------------------------------------------------------------------------------------
// Mutually exclusive alternatives
// ---------------------------------------------------------------------------------------------------------------

/**
 * Incremental (challenger–defender) analysis of mutually exclusive alternatives with equal lives. Alternatives
 * are ordered by first cost; each challenger is accepted when the incremental flow (challenger − defender) has
 * NPV >= 0 at the MARR (equivalent to incremental IRR >= MARR for a conventional increment; the incremental
 * IRR status is reported so that non-conventional increments are visible).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 8 (incremental ROR analysis).
 * @param {{alternatives: {name: string, cashflows: number[]}[], marr: number, includeDoNothing?: boolean}} o
 *   marr = minimum attractive rate of return (fraction); do-nothing is the first defender unless disabled
 * @returns {{selected: string | null, steps: {defender: string, challenger: string, incrementalNpv: number,
 *   incrementalIrr: number | null, irrStatus: string, accepted: boolean}[]}}
 */
export function incrementalAnalysis({ alternatives, marr, includeDoNothing = true }) {
  if (!Array.isArray(alternatives) || alternatives.length === 0) throw new RangeError('alternatives must be a non-empty array');
  requireNumber('marr', marr, { gt: -1 });
  const n = alternatives[0].cashflows.length;
  for (const a of alternatives) {
    requireSeries(`cashflows of ${a.name}`, a.cashflows, { minLength: 2 });
    if (a.cashflows.length !== n) throw new RangeError('alternatives must have equal lives; use compareUnequalLives');
  }
  const sorted = alternatives.slice().sort((p, q) => -p.cashflows[0] - -q.cashflows[0]);
  let defender = includeDoNothing ? { name: 'Do nothing', cashflows: new Array(n).fill(0) } : sorted.shift();
  const steps = [];
  for (const ch of sorted) {
    const delta = ch.cashflows.map((v, k) => v - defender.cashflows[k]);
    const dn = npv(delta, marr);
    const di = irr(delta);
    const accepted = dn >= 0;
    steps.push({ defender: defender.name, challenger: ch.name, incrementalNpv: dn, incrementalIrr: di.irr,
      irrStatus: di.status, accepted });
    if (accepted) defender = ch;
  }
  return { selected: defender.name === 'Do nothing' ? null : defender.name, steps };
}

function gcd(a, b) { while (b) { [a, b] = [b, a % b]; } return a; }

/**
 * Comparison of alternatives with unequal lives under the repeatability assumption: equivalent annual worth of
 * one life cycle, and NPV over the least common multiple of lives with like-for-like renewal.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 5 (LCM) and ch. 6 (annual worth).
 * @param {{alternatives: {name: string, cashflows: number[]}[], rate: number}} o cashflows over one life cycle
 *   (life = cashflows.length - 1 periods), rate = MARR (fraction)
 * @returns {{studyPeriod: number, results: {name: string, life: number, npvOneCycle: number, eaw: number,
 *   npvCommonMultiple: number}[], bestByEaw: string}}
 */
export function compareUnequalLives({ alternatives, rate }) {
  if (!Array.isArray(alternatives) || alternatives.length === 0) throw new RangeError('alternatives must be a non-empty array');
  requireNumber('rate', rate, { gt: -1 });
  let lcm = 1;
  for (const a of alternatives) {
    requireSeries(`cashflows of ${a.name}`, a.cashflows, { minLength: 2 });
    const life = a.cashflows.length - 1;
    lcm = lcm * life / gcd(lcm, life);
  }
  const results = alternatives.map((a) => {
    const life = a.cashflows.length - 1;
    const one = npv(a.cashflows, rate);
    let rep = 0;
    for (let c = 0; c < lcm / life; c++) rep += one * Math.pow(1 + rate, -c * life);
    return { name: a.name, life, npvOneCycle: one, eaw: one * factorAP(rate, life), npvCommonMultiple: rep };
  });
  let best = results[0];
  for (const r of results) if (r.eaw > best.eaw) best = r;
  return { studyPeriod: lcm, results, bestByEaw: best.name };
}

// ---------------------------------------------------------------------------------------------------------------
// Currency conversion
// ---------------------------------------------------------------------------------------------------------------

/**
 * Convert an amount between currencies using a rate table quoted as units of each currency per one unit of a
 * common reference currency, with a base-date consistency check: the amount's price date and the rate date must
 * agree within `toleranceDays`, otherwise a RangeError is thrown (or a warning returned with allowDateMismatch).
 * Ref: consistent-base-date rule for cost data, AACE RP 58R-10 (escalation and currency in estimates).
 * @param {{amount: number, from: string, to: string, rates: Record<string, number>, amountDate?: string,
 *   rateDate?: string, toleranceDays?: number, allowDateMismatch?: boolean}} o rates[code] = units of `code` per
 *   reference unit (> 0); dates are ISO strings
 * @returns {{amount: number, rate: number, warning: string | null}} converted amount (currency `to`) and the
 *   cross rate (units of `to` per unit of `from`)
 */
export function convertCurrency({ amount, from, to, rates, amountDate, rateDate, toleranceDays = 31,
  allowDateMismatch = false }) {
  requireNumber('amount', amount);
  if (!rates || !(rates[from] > 0)) throw new RangeError(`rates has no positive entry for '${from}'`);
  if (!(rates[to] > 0)) throw new RangeError(`rates has no positive entry for '${to}'`);
  let warning = null;
  if (amountDate && rateDate) {
    const d = Math.abs(Date.parse(amountDate) - Date.parse(rateDate)) / 86400000;
    if (!Number.isFinite(d)) throw new RangeError('amountDate and rateDate must be ISO dates');
    if (d > toleranceDays) {
      const msg = `exchange-rate date ${rateDate} differs from the amount's price date ${amountDate} by ${Math.round(d)} days`;
      if (!allowDateMismatch) throw new RangeError(msg);
      warning = msg;
    }
  }
  const rate = rates[to] / rates[from];
  return { amount: amount * rate, rate, warning };
}

// ---------------------------------------------------------------------------------------------------------------
// Project cash-flow model
// ---------------------------------------------------------------------------------------------------------------

const seriesAt = (s, k, name) => {
  if (typeof s === 'number') return s;
  if (!s || !s.length) throw new RangeError(`${name} must be a number or a non-empty series`);
  return k < s.length ? s[k] : s[s.length - 1];
};

/**
 * Annual project cash-flow model (free cash flow to firm, unlevered, after tax).
 *
 * Timeline: rows t = 0 … ; construction occupies rows 0 … nc-1, an optional start-up delay follows, operations
 * run for `operatingYears` rows, decommissioning is spent in row `decomYear` (default: the row after the last
 * operating year). Row 0 is time zero for discounting.
 *
 * Money inputs (capex amounts, opex, merchant prices, salvage, decommissioning cost) are in constant money of the
 * price base date (row `priceBaseOffset`, default row 0) and are inflated internally with `inflation` plus any
 * item-specific real `escalation`. The model is always evaluated in nominal money so that depreciation, loss
 * carry-forward and fixed-price contracts are treated correctly, then — in 'real' mode — every money column is
 * deflated to base-date money. `discountRate` must be stated in the same terms as `mode`; with
 * r_nominal = (1 + r_real)(1 + inflation) − 1 both modes give the same NPV.
 *
 * Price contracts (`price`): a number or { type: 'indexed', price, indexation? } — constant real price escalated
 * at `indexation` (default inflation); { type: 'fixed', price } — flat nominal price; { type: 'merchant', series,
 * nominal? } — market series per operating year; { type: 'cfd', strike, strikeIndexation?, reference, capture?,
 * twoWay?, termYears? } — contract for difference: realised = capture + (strike − reference) for a two-way CfD,
 * capture + max(0, strike − reference) for one-way, market after the term; { type: 'ppa', price? | floor?, cap?,
 * market, discount?, contractedFraction?, indexation?, termYears? } — fixed-price PPA, or a floor/cap on the
 * discounted market price, for the contracted share of output.
 *
 * Tax: taxable income = EBITDA − depreciation − investment allowance − decommissioning deduction + salvage; losses
 * are carried forward (and optionally carried back against tax already paid); an investment tax credit is offset
 * against tax payable and carried forward. Sunk capex items are excluded from cash flow (when `excludeSunk`) but
 * still depreciate, because their tax shield is a future cash flow.
 *
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 6 (cash-flow rules: incremental after-tax
 * flows, sunk costs, working capital, consistent inflation treatment); Blank & Tarquin ch. 14 and 17.
 *
 * @param {object} inp
 * @param {{constructionYears: number, operatingYears: number, startDelayYears?: number, decomYear?: number,
 *   baseYear?: number}} inp.timeline integers (years); baseYear is the calendar year of row 0
 * @param {'nominal' | 'real'} [inp.mode] reporting basis (default 'nominal')
 * @param {number} [inp.inflation] general inflation (fraction per year, default 0)
 * @param {number} [inp.priceBaseOffset] row of the price base date relative to row 0 (default 0; may be negative)
 * @param {number | number[]} inp.discountRate fraction per year, in the terms of `mode`
 * @param {boolean} [inp.midYear] mid-year discounting (default false)
 * @param {{name?: string, category?: string, amount: number, phasing?: number[], startYear?: number, year?: number,
 *   escalation?: number, depreciable?: boolean, sunk?: boolean}[]} inp.capex amounts in base-date money (sunk
 *   items: historical cost, not inflated); phasing fractions sum to 1 over consecutive rows from startYear
 * @param {boolean} [inp.excludeSunk] default true
 * @param {{fixed?: number | number[], variablePerUnit?: number | number[], escalation?: number}} [inp.opex]
 *   fixed per operating year, variable per unit of production (base-date money), real escalation (fraction/yr)
 * @param {number | number[]} inp.production units per operating year
 * @param {number | object} inp.price see above (currency per unit)
 * @param {number} [inp.royaltyRate] fraction of revenue
 * @param {{rate?: number, lossCarryForward?: boolean, lossCarryBack?: boolean, depreciation?: object,
 *   investmentAllowance?: number, investmentTaxCredit?: number}} [inp.tax] depreciation = method spec for
 *   depreciationSchedule without `basis`
 * @param {{fractionOfRevenue?: number}} [inp.workingCapital] recovered in the last operating year
 * @param {number} [inp.salvage] residual value received in the last operating year (base-date money, taxable)
 * @param {{cost: number, escalation?: number, fund?: {rate?: number, startOperatingYear?: number},
 *   deductible?: 'spend' | 'contribution' | 'none'}} [inp.decommissioning] cost in base-date money; with `fund`
 *   level nominal contributions accrue in an escrow earning `rate` so that it equals the cost at decomYear
 * @param {number} [inp.delayCostPerYear] standby cost per start-up delay year (base-date money)
 * @param {{price?: number, production?: number, capex?: number, opex?: number}} [inp.factors] multipliers used by
 *   break-even, sensitivity and Monte Carlo studies (default 1)
 * @returns {{mode: string, n: number, rows: object[], table: Record<string, number[]>, fcff: number[],
 *   metrics: object, totals: object, timeline: object}} table columns (currency per year unless stated): t,
 *   calendarYear, phase, production (units), price (currency/unit), revenue, royalty, opex, ebitda, depreciation,
 *   taxableIncome, tax, capex, deltaWC, decomFund, decomFundBalance, decomSpend, salvage, fcff, cumulative,
 *   discountFactor (-), discounted, cumulativeDiscounted
 */
export function buildProjectCashflow(inp) {
  if (!inp || !inp.timeline) throw new RangeError('timeline is required');
  const tl = inp.timeline;
  const nc = requireNumber('timeline.constructionYears', tl.constructionYears, { min: 0, integer: true });
  const no = requireNumber('timeline.operatingYears', tl.operatingYears, { min: 1, integer: true });
  const nd = requireNumber('timeline.startDelayYears', tl.startDelayYears ?? 0, { min: 0, integer: true });
  const opStart = nc + nd;
  const opEnd = opStart + no - 1;
  const decom = inp.decommissioning && inp.decommissioning.cost ? inp.decommissioning : null;
  const decomYear = decom ? requireNumber('timeline.decomYear', tl.decomYear ?? opEnd + 1, { min: 0, integer: true }) : null;
  const n = Math.max(opEnd, decomYear ?? 0) + 1;
  const mode = inp.mode ?? 'nominal';
  if (mode !== 'nominal' && mode !== 'real') throw new RangeError("mode must be 'nominal' or 'real'");
  const infl = requireNumber('inflation', inp.inflation ?? 0, { gt: -1 });
  const base = requireNumber('priceBaseOffset', inp.priceBaseOffset ?? 0);
  const fac = { price: 1, production: 1, capex: 1, opex: 1, ...(inp.factors || {}) };
  const idx = new Float64Array(n);
  for (let t = 0; t < n; t++) idx[t] = Math.pow(1 + infl, t - base);
  const esc = (rate, t) => Math.pow(1 + (rate ?? 0), t - base);
  const col = () => new Array(n).fill(0);

  // --- capital expenditure and depreciable cohorts (nominal) ---
  const capex = col(); const depBasis = col(); const capexByCategory = {};
  if (!Array.isArray(inp.capex)) throw new RangeError('capex must be an array of items');
  for (const item of inp.capex) {
    requireNumber(`capex amount (${item.name ?? 'item'})`, item.amount, { min: 0 });
    let start = item.year ?? item.startYear ?? 0;
    let ph = item.phasing;
    if (item.year !== undefined || !ph) {
      ph = item.year !== undefined || nc <= 1 ? [1] : new Array(nc).fill(1 / nc);
    }
    requireSeries('capex phasing', ph);
    let sum = 0; for (const f of ph) sum += f;
    if (Math.abs(sum - 1) > 1e-6) throw new RangeError(`capex phasing of '${item.name ?? 'item'}' must sum to 1`);
    requireNumber('capex startYear', start, { min: 0, integer: true });
    if (start + ph.length > n) throw new RangeError(`capex item '${item.name ?? 'item'}' extends beyond the timeline`);
    for (let j = 0; j < ph.length; j++) {
      const t = start + j;
      const nominal = item.sunk ? item.amount * ph[j]
        : item.amount * fac.capex * ph[j] * idx[t] * esc(item.escalation, t);
      if (item.depreciable !== false) depBasis[t] += nominal;
      if (item.sunk && inp.excludeSunk !== false) continue;
      capex[t] += nominal;
      const cat = item.category ?? 'uncategorised';
      capexByCategory[cat] = (capexByCategory[cat] ?? 0) + nominal;
    }
  }

  // --- production, price, revenue, opex (nominal) ---
  const production = col(); const price = col(); const revenue = col(); const royalty = col(); const opex = col();
  const settlement = col(); const phase = new Array(n).fill('post');
  const P = inp.price;
  const pSpec = typeof P === 'number' ? { type: 'indexed', price: P } : P;
  if (!pSpec || !pSpec.type) throw new RangeError('price must be a number or a contract object with a type');
  const ox = inp.opex || {};
  for (let t = 0; t < n; t++) {
    phase[t] = t < nc ? 'construction' : t < opStart ? 'delay' : t <= opEnd ? 'operation' : t === decomYear ? 'decommissioning' : 'post';
    if (t >= nc && t < opStart) opex[t] = (inp.delayCostPerYear ?? 0) * idx[t];
    if (t < opStart || t > opEnd) continue;
    const k = t - opStart;
    production[t] = seriesAt(inp.production, k, 'production') * fac.production;
    const mkt = (s) => seriesAt(s, k, 'price series') * (pSpec.nominal ? 1 : idx[t]);
    const index = (v, rate) => v * Math.pow(1 + (rate ?? infl), t - base);
    const inTerm = pSpec.termYears === undefined || k < pSpec.termYears;
    let p;
    switch (pSpec.type) {
      case 'fixed': p = requireNumber('price.price', pSpec.price); break;
      case 'indexed': p = index(requireNumber('price.price', pSpec.price), pSpec.indexation); break;
      case 'merchant': p = mkt(pSpec.series); break;
      case 'cfd': {
        const ref = mkt(pSpec.reference);
        const cap = pSpec.capture !== undefined ? mkt(pSpec.capture) : ref;
        if (inTerm) {
          const strike = index(requireNumber('price.strike', pSpec.strike), pSpec.strikeIndexation);
          const diff = pSpec.twoWay === false ? Math.max(0, strike - ref) : strike - ref;
          p = cap + diff;
          settlement[t] = diff * production[t] * fac.price;
        } else p = cap;
        break;
      }
      case 'ppa': {
        const m = mkt(pSpec.market ?? pSpec.price ?? 0);
        if (inTerm) {
          let c;
          if (pSpec.price !== undefined) c = index(pSpec.price, pSpec.indexation);
          else {
            c = m * (1 - (pSpec.discount ?? 0));
            if (pSpec.floor !== undefined) c = Math.max(c, index(pSpec.floor, pSpec.indexation));
            if (pSpec.cap !== undefined) c = Math.min(c, index(pSpec.cap, pSpec.indexation));
          }
          const share = pSpec.contractedFraction ?? 1;
          p = share * c + (1 - share) * m;
        } else p = m;
        break;
      }
      default: throw new RangeError(`price.type '${pSpec.type}' is not supported`);
    }
    price[t] = p * fac.price;
    revenue[t] = price[t] * production[t];
    royalty[t] = (inp.royaltyRate ?? 0) * revenue[t];
    const fixed = ox.fixed !== undefined ? seriesAt(ox.fixed, k, 'opex.fixed') : 0;
    const varU = ox.variablePerUnit !== undefined ? seriesAt(ox.variablePerUnit, k, 'opex.variablePerUnit') : 0;
    opex[t] = (fixed + varU * production[t]) * fac.opex * idx[t] * esc(ox.escalation, t);
  }

  // --- depreciation by cohort; starts when the asset enters service ---
  const tax = inp.tax || {};
  const taxRate = requireNumber('tax.rate', tax.rate ?? 0, { min: 0, max: 1 });
  const depreciation = col(); const allowance = col(); const credit = col();
  let undepreciated = 0;
  for (let y = 0; y < n; y++) {
    if (depBasis[y] <= 0) continue;
    const s = Math.max(y, opStart);
    if (tax.depreciation && tax.depreciation.method !== 'none') {
      const ch = depreciationSchedule({ ...tax.depreciation, basis: depBasis[y] });
      let used = 0;
      for (let j = 0; j < ch.length && s + j < n; j++) { depreciation[s + j] += ch[j]; used += ch[j]; }
      undepreciated += depBasis[y] - used;
    }
    if (s < n) {
      allowance[s] += (tax.investmentAllowance ?? 0) * depBasis[y];
      credit[s] += (tax.investmentTaxCredit ?? 0) * depBasis[y];
    }
  }

  // --- working capital, salvage, decommissioning ---
  const deltaWC = col(); const salvage = col(); const decomSpend = col(); const decomFund = col();
  const decomFundBalance = col(); const decomDeduction = col();
  const wcFrac = inp.workingCapital?.fractionOfRevenue ?? 0;
  let wcPrev = 0;
  for (let t = opStart; t <= opEnd; t++) {
    const wc = t === opEnd ? 0 : wcFrac * revenue[t];
    deltaWC[t] = wc - wcPrev; wcPrev = wc;
  }
  salvage[opEnd] = (inp.salvage ?? 0) * idx[opEnd];
  let decomNominal = 0;
  if (decom) {
    requireNumber('decommissioning.cost', decom.cost, { min: 0 });
    decomNominal = decom.cost * idx[decomYear] * esc(decom.escalation, decomYear);
    const how = decom.deductible ?? (decom.fund ? 'contribution' : 'spend');
    if (decom.fund) {
      const fr = requireNumber('decommissioning.fund.rate', decom.fund.rate ?? 0, { gt: -1 });
      const s0 = opStart + (decom.fund.startOperatingYear ?? 0);
      const s1 = Math.min(opEnd, decomYear);
      if (s0 > s1) throw new RangeError('decommissioning.fund.startOperatingYear leaves no contribution years');
      let ann = 0;
      for (let t = s0; t <= s1; t++) ann += Math.pow(1 + fr, decomYear - t);
      const A = decomNominal / ann;
      let bal = 0;
      for (let t = 0; t < n; t++) {
        bal *= t > 0 ? 1 + fr : 1;
        if (t >= s0 && t <= s1) { decomFund[t] = A; bal += A; }
        if (t === decomYear) bal -= decomNominal;
        decomFundBalance[t] = Math.abs(bal) < 1e-9 * Math.max(1, decomNominal) ? 0 : bal;
        if (how === 'contribution') decomDeduction[t] = decomFund[t];
      }
      decomSpend[decomYear] = decomNominal;      // paid from the escrow: no project cash flow in that year
    } else decomSpend[decomYear] = decomNominal;
    if (how === 'spend') decomDeduction[decomYear] = decomNominal;
  }

  // --- tax ---
  const ebitda = col(); const taxableIncome = col(); const taxPaid = col(); const fcff = col();
  let lossPool = 0; let creditPool = 0; let cumTax = 0;
  const lcf = tax.lossCarryForward !== false;
  for (let t = 0; t < n; t++) {
    ebitda[t] = revenue[t] - royalty[t] - opex[t];
    const ti = ebitda[t] - depreciation[t] - allowance[t] - decomDeduction[t] + salvage[t];
    taxableIncome[t] = ti;
    let due = 0;
    if (ti < 0) {
      let loss = -ti;
      if (tax.lossCarryBack && cumTax > 0) {
        const refund = Math.min(cumTax, taxRate * loss);
        due = -refund;
        loss -= taxRate > 0 ? refund / taxRate : 0;
      }
      if (lcf) lossPool += loss;
    } else {
      const use = lcf ? Math.min(lossPool, ti) : 0;
      lossPool -= use;
      due = taxRate * (ti - use);
    }
    creditPool += credit[t];
    if (due > 0 && creditPool > 0) { const c = Math.min(creditPool, due); due -= c; creditPool -= c; }
    taxPaid[t] = due; cumTax += due;
    const decomCash = decom && decom.fund ? 0 : decomSpend[t];
    fcff[t] = ebitda[t] - due - capex[t] - deltaWC[t] - decomFund[t] - decomCash + salvage[t];
  }

  // --- real mode: deflate every money column to base-date money ---
  const money = { price, revenue, royalty, opex, ebitda, depreciation, taxableIncome, tax: taxPaid, capex, deltaWC,
    decomFund, decomFundBalance, decomSpend, salvage, settlement, fcff };
  if (mode === 'real') {
    for (const key of Object.keys(money)) { const a = money[key]; for (let t = 0; t < n; t++) a[t] /= idx[t]; }
    for (const k of Object.keys(capexByCategory)) capexByCategory[k] = NaN;   // mixed-date sum is meaningless in real terms
  }

  const rate = inp.discountRate;
  if (rate === undefined) throw new RangeError('discountRate is required');
  const ds = discountSeries(fcff, rate, { midYear: !!inp.midYear });
  const table = { t: Array.from({ length: n }, (_, t) => t),
    calendarYear: Array.from({ length: n }, (_, t) => (tl.baseYear ?? 0) + t), phase, production, ...money,
    cumulative: ds.cumulative, discountFactor: ds.factors, discounted: ds.discounted,
    cumulativeDiscounted: ds.cumulativeDiscounted };
  const keys = Object.keys(table);
  const rows = [];
  for (let t = 0; t < n; t++) { const r = {}; for (const k of keys) r[k] = table[k][t]; rows.push(r); }
  const sum = (a) => a.reduce((s, v) => s + v, 0);
  const metrics = cashflowMetrics(fcff, rate, { midYear: !!inp.midYear, investment: capex,
    financeRate: inp.financeRate, reinvestRate: inp.reinvestRate });
  const totals = { capex: sum(capex), opex: sum(opex), revenue: sum(revenue), royalty: sum(royalty), tax: sum(taxPaid),
    production: sum(production), decommissioning: mode === 'real' ? sum(decomSpend) : decomNominal,
    pvCapex: npv(capex, rate, { midYear: !!inp.midYear }), pvOpex: npv(opex, rate, { midYear: !!inp.midYear }),
    pvRevenue: npv(revenue, rate, { midYear: !!inp.midYear }),
    pvProduction: npv(production, rate, { midYear: !!inp.midYear }),
    undepreciatedBasis: undepreciated, unusedLosses: lossPool, unusedCredits: creditPool,
    capexByCategory: mode === 'real' ? null : capexByCategory };
  return { mode, n, rows, table, fcff, metrics, totals,
    timeline: { constructionYears: nc, startDelayYears: nd, operatingYears: no, firstOperatingRow: opStart,
      lastOperatingRow: opEnd, decomYear } };
}
