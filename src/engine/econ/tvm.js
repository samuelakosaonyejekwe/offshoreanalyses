// Time value of money: compounding, the standard engineering-economy factors, annuities, gradients,
// equivalent annual worth and real/nominal (Fisher) conversions.
// Money is a plain number in one project currency; rates are fractions per period (0.08 = 8 %); n is in periods
// (years unless the caller uses another consistent period).

/**
 * Validate a finite number against optional bounds; throws RangeError naming the argument. Returns the value.
 * Ref: house input-validation rule (docs/CONTRACT.md §1).
 * @param {string} name argument name used in the message
 * @param {number} v value (any unit)
 * @param {{min?: number, max?: number, gt?: number, lt?: number, integer?: boolean}} [o] inclusive min/max, exclusive gt/lt
 * @returns {number} v
 */
export function requireNumber(name, v, o = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new RangeError(`${name} must be a finite number`);
  if (o.min !== undefined && v < o.min) throw new RangeError(`${name} must be >= ${o.min}`);
  if (o.max !== undefined && v > o.max) throw new RangeError(`${name} must be <= ${o.max}`);
  if (o.gt !== undefined && !(v > o.gt)) throw new RangeError(`${name} must be > ${o.gt}`);
  if (o.lt !== undefined && !(v < o.lt)) throw new RangeError(`${name} must be < ${o.lt}`);
  if (o.integer && !Number.isInteger(v)) throw new RangeError(`${name} must be an integer`);
  return v;
}

/**
 * Validate an array of finite numbers; throws RangeError naming the argument. Returns the array.
 * Ref: house input-validation rule (docs/CONTRACT.md §1).
 * @param {string} name argument name
 * @param {ArrayLike<number>} a values (any unit)
 * @param {{minLength?: number}} [o]
 * @returns {ArrayLike<number>} a
 */
export function requireSeries(name, a, o = {}) {
  if (!a || typeof a.length !== 'number') throw new RangeError(`${name} must be an array of numbers`);
  if (a.length < (o.minLength ?? 1)) throw new RangeError(`${name} must have at least ${o.minLength ?? 1} element(s)`);
  for (let k = 0; k < a.length; k++) {
    if (typeof a[k] !== 'number' || !Number.isFinite(a[k])) throw new RangeError(`${name}[${k}] must be a finite number`);
  }
  return a;
}

const rate = (i) => requireNumber('i', i, { gt: -1 });
const periods = (n) => requireNumber('n', n, { min: 0 });

// ---------------------------------------------------------------------------------------------------------------
// Standard factors (ANSI Z94 functional notation)
// ---------------------------------------------------------------------------------------------------------------

/**
 * Single-payment compound-amount factor (F/P, i, n) = (1 + i)^n. Dimensionless.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2; Park, Contemporary Engineering Economics, ch. 3.
 * @param {number} i interest rate per period (fraction, > -1)
 * @param {number} n number of periods (>= 0, may be fractional)
 * @returns {number}
 */
export function factorFP(i, n) { return Math.pow(1 + rate(i), periods(n)); }

/**
 * Single-payment present-worth factor (P/F, i, n) = (1 + i)^-n. Dimensionless.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number}
 */
export function factorPF(i, n) { return Math.pow(1 + rate(i), -periods(n)); }

/**
 * Uniform-series present-worth factor (P/A, i, n) = [(1+i)^n - 1] / [i (1+i)^n]; equals n when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number} dimensionless
 */
export function factorPA(i, n) {
  rate(i); periods(n);
  if (Math.abs(i) < 1e-12) return n;
  return -Math.expm1(-n * Math.log1p(i)) / i;
}

/**
 * Capital-recovery factor (A/P, i, n) = i (1+i)^n / [(1+i)^n - 1]; equals 1/n when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (> 0)
 * @returns {number} dimensionless (per period)
 */
export function factorAP(i, n) { requireNumber('n', n, { gt: 0 }); return 1 / factorPA(i, n); }

/**
 * Uniform-series compound-amount factor (F/A, i, n) = [(1+i)^n - 1] / i; equals n when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number} dimensionless
 */
export function factorFA(i, n) {
  rate(i); periods(n);
  if (Math.abs(i) < 1e-12) return n;
  return Math.expm1(n * Math.log1p(i)) / i;
}

/**
 * Sinking-fund factor (A/F, i, n) = i / [(1+i)^n - 1]; equals 1/n when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (> 0)
 * @returns {number} dimensionless (per period)
 */
export function factorAF(i, n) { requireNumber('n', n, { gt: 0 }); return 1 / factorFA(i, n); }

/**
 * Arithmetic-gradient present-worth factor (P/G, i, n) = [(1+i)^n - i n - 1] / [i² (1+i)^n].
 * The gradient series is 0, G, 2G, ..., (n-1)G at periods 1..n. Equals n(n-1)/2 when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number} dimensionless
 */
export function factorPG(i, n) {
  rate(i); periods(n);
  if (Math.abs(i) < 1e-9) return n * (n - 1) / 2;
  const f = Math.pow(1 + i, n);
  return (f - i * n - 1) / (i * i * f);
}

/**
 * Arithmetic-gradient uniform-series factor (A/G, i, n) = 1/i - n / [(1+i)^n - 1]; equals (n-1)/2 when i = 0.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (> 0)
 * @returns {number} dimensionless
 */
export function factorAG(i, n) {
  rate(i); requireNumber('n', n, { gt: 0 });
  if (Math.abs(i) < 1e-9) return (n - 1) / 2;
  return 1 / i - n / Math.expm1(n * Math.log1p(i));
}

/**
 * Arithmetic-gradient future-worth factor (F/G, i, n) = (P/G, i, n) (F/P, i, n).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number} dimensionless
 */
export function factorFG(i, n) { return factorPG(i, n) * factorFP(i, n); }

/**
 * Geometric-gradient present-worth factor (P/A1, g, i, n): PV of A1, A1(1+g), ..., A1(1+g)^(n-1) at periods 1..n
 * per unit A1. = [1 - ((1+g)/(1+i))^n] / (i - g), or n/(1+i) when g = i.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} g growth rate per period (fraction, > -1)
 * @param {number} n periods (>= 0)
 * @returns {number} dimensionless
 */
export function factorPGeometric(i, g, n) {
  rate(i); requireNumber('g', g, { gt: -1 }); periods(n);
  if (Math.abs(i - g) < 1e-12) return n / (1 + i);
  return (1 - Math.pow((1 + g) / (1 + i), n)) / (i - g);
}

/**
 * The six standard engineering-economy factors plus the gradient factors for one (i, n) pair. Dimensionless.
 * Ref: Blank & Tarquin, Engineering Economy, compound-interest factor tables.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} n periods (> 0)
 * @returns {{FP: number, PF: number, AF: number, FA: number, AP: number, PA: number, PG: number, AG: number}}
 */
export function factorTable(i, n) {
  return { FP: factorFP(i, n), PF: factorPF(i, n), AF: factorAF(i, n), FA: factorFA(i, n), AP: factorAP(i, n),
    PA: factorPA(i, n), PG: factorPG(i, n), AG: factorAG(i, n) };
}

// ---------------------------------------------------------------------------------------------------------------
// Compounding and rate conversions
// ---------------------------------------------------------------------------------------------------------------

/**
 * Future value of a single sum with discrete compounding m times per year: F = P (1 + r/m)^(m t). Currency.
 * Ref: Park, Contemporary Engineering Economics, ch. 4.
 * @param {{P: number, r: number, t: number, m?: number}} o P present sum (currency), r nominal annual rate
 *   (fraction), t years (>= 0), m compounding periods per year (default 1)
 * @returns {number} future value (currency)
 */
export function compoundAmount({ P, r, t, m = 1 }) {
  requireNumber('P', P); requireNumber('t', t, { min: 0 }); requireNumber('m', m, { gt: 0 });
  requireNumber('r', r, { gt: -m });
  return P * Math.pow(1 + r / m, m * t);
}

/**
 * Future value with continuous compounding F = P e^(r t). Currency.
 * Ref: Park, Contemporary Engineering Economics, ch. 4.
 * @param {number} P present sum (currency)
 * @param {number} r nominal annual rate, continuously compounded (fraction)
 * @param {number} t years (>= 0)
 * @returns {number} future value (currency)
 */
export function fvContinuous(P, r, t) {
  requireNumber('P', P); requireNumber('r', r); requireNumber('t', t, { min: 0 });
  return P * Math.exp(r * t);
}

/**
 * Present value with continuous discounting P = F e^(-r t). Currency.
 * Ref: Park, Contemporary Engineering Economics, ch. 4.
 * @param {number} F future sum (currency)
 * @param {number} r continuously compounded annual rate (fraction)
 * @param {number} t years (>= 0)
 * @returns {number} present value (currency)
 */
export function pvContinuous(F, r, t) {
  requireNumber('F', F); requireNumber('r', r); requireNumber('t', t, { min: 0 });
  return F * Math.exp(-r * t);
}

/**
 * Effective annual rate from a nominal annual rate compounded m times per year: (1 + r/m)^m - 1;
 * m = Infinity gives continuous compounding e^r - 1.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 4.
 * @param {number} r nominal annual rate (fraction)
 * @param {number} [m] compounding periods per year (> 0 or Infinity, default 1)
 * @returns {number} effective annual rate (fraction)
 */
export function effectiveRate(r, m = 1) {
  requireNumber('r', r);
  if (m === Infinity) return Math.expm1(r);
  requireNumber('m', m, { gt: 0 });
  return Math.pow(1 + r / m, m) - 1;
}

/**
 * Nominal annual rate (compounded m times per year) equivalent to an effective annual rate:
 * m [(1 + ia)^(1/m) - 1]; m = Infinity gives ln(1 + ia).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 4.
 * @param {number} ia effective annual rate (fraction, > -1)
 * @param {number} [m] compounding periods per year (> 0 or Infinity, default 1)
 * @returns {number} nominal annual rate (fraction)
 */
export function nominalRate(ia, m = 1) {
  requireNumber('ia', ia, { gt: -1 });
  if (m === Infinity) return Math.log1p(ia);
  requireNumber('m', m, { gt: 0 });
  return m * (Math.pow(1 + ia, 1 / m) - 1);
}

/**
 * Effective rate per sub-period equivalent to an effective annual rate: (1 + ia)^(1/m) - 1.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 4.
 * @param {number} ia effective annual rate (fraction, > -1)
 * @param {number} m sub-periods per year (> 0)
 * @returns {number} effective rate per sub-period (fraction)
 */
export function periodicRate(ia, m) {
  requireNumber('ia', ia, { gt: -1 }); requireNumber('m', m, { gt: 0 });
  return Math.pow(1 + ia, 1 / m) - 1;
}

// ---------------------------------------------------------------------------------------------------------------
// Single sums, annuities, perpetuities
// ---------------------------------------------------------------------------------------------------------------

/**
 * Future value of a single sum F = P (1 + i)^n. Currency.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} P present sum (currency)
 * @param {number} i rate per period (fraction)
 * @param {number} n periods
 * @returns {number} currency
 */
export function futureValue(P, i, n) { return requireNumber('P', P) * factorFP(i, n); }

/**
 * Present value of a single sum P = F (1 + i)^-n. Currency.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {number} F future sum (currency)
 * @param {number} i rate per period (fraction)
 * @param {number} n periods
 * @returns {number} currency
 */
export function presentValue(F, i, n) { return requireNumber('F', F) * factorPF(i, n); }

/**
 * Present value of a level annuity of n payments A. Ordinary (end of period) by default; an annuity due
 * (start of period) is the ordinary value times (1 + i).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2; Brealey, Myers & Allen, Principles of Corporate Finance, ch. 2.
 * @param {{A: number, i: number, n: number, due?: boolean}} o A payment per period (currency), i rate per period
 *   (fraction), n number of payments, due true for payments at period start
 * @returns {number} currency
 */
export function pvAnnuity({ A, i, n, due = false }) {
  return requireNumber('A', A) * factorPA(i, n) * (due ? 1 + i : 1);
}

/**
 * Future value (at the end of period n) of a level annuity of n payments A, ordinary or due.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {{A: number, i: number, n: number, due?: boolean}} o as for pvAnnuity
 * @returns {number} currency
 */
export function fvAnnuity({ A, i, n, due = false }) {
  return requireNumber('A', A) * factorFA(i, n) * (due ? 1 + i : 1);
}

/**
 * Level payment that repays (ordinary annuity) a present sum P over n periods: A = P (A/P, i, n); for an annuity
 * due the payment is divided by (1 + i).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 2.
 * @param {{P: number, i: number, n: number, due?: boolean}} o P present sum (currency), i rate (fraction), n periods
 * @returns {number} payment per period (currency)
 */
export function annuityPayment({ P, i, n, due = false }) {
  return requireNumber('P', P) * factorAP(i, n) / (due ? 1 + i : 1);
}

/**
 * Present value of a growing annuity: first payment A1 at period 1, growing at g per period for n payments.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 2; Blank & Tarquin ch. 2.
 * @param {{A1: number, i: number, g: number, n: number, due?: boolean}} o A1 first payment (currency), i discount
 *   rate, g growth rate (fractions per period), n payments, due true when A1 is paid at time 0
 * @returns {number} currency
 */
export function pvGrowingAnnuity({ A1, i, g, n, due = false }) {
  return requireNumber('A1', A1) * factorPGeometric(i, g, n) * (due ? 1 + i : 1);
}

/**
 * Present value of a level perpetuity A / i (first payment one period from now).
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 2.
 * @param {number} A payment per period (currency)
 * @param {number} i rate per period (fraction, > 0)
 * @returns {number} currency
 */
export function pvPerpetuity(A, i) { requireNumber('i', i, { gt: 0 }); return requireNumber('A', A) / i; }

/**
 * Present value of a growing perpetuity A1 / (i - g) (Gordon growth); requires i > g.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 2.
 * @param {number} A1 payment at the end of period 1 (currency)
 * @param {number} i discount rate per period (fraction)
 * @param {number} g growth rate per period (fraction, < i)
 * @returns {number} currency
 */
export function pvGrowingPerpetuity(A1, i, g) {
  requireNumber('A1', A1); requireNumber('i', i); requireNumber('g', g);
  if (!(i > g)) throw new RangeError('i must be greater than g for a growing perpetuity');
  return A1 / (i - g);
}

/**
 * Capitalised cost: present value of an asset renewed for ever — first cost plus PV of a perpetual annual cost and
 * of renewals every `life` periods.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 5 (capitalised cost).
 * @param {{firstCost: number, annualCost?: number, renewalCost?: number, life?: number, i: number}} o currency,
 *   currency per period, currency per renewal, periods between renewals, rate per period (> 0)
 * @returns {number} currency
 */
export function capitalisedCost({ firstCost, annualCost = 0, renewalCost = 0, life = 1, i }) {
  requireNumber('firstCost', firstCost); requireNumber('i', i, { gt: 0 });
  const renew = renewalCost ? renewalCost * factorAF(i, requireNumber('life', life, { gt: 0 })) / i : 0;
  return firstCost + annualCost / i + renew;
}

// ---------------------------------------------------------------------------------------------------------------
// Equivalent annual worth / cost
// ---------------------------------------------------------------------------------------------------------------

/**
 * Equivalent (uniform) annual worth of a cash-flow series: AW = NPV (A/P, i, n), n = last period index.
 * cashflows[0] occurs at time 0, cashflows[k] at the end of period k.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 6.
 * @param {number[]} cashflows currency per period (inflows positive)
 * @param {number} i rate per period (fraction, > -1)
 * @returns {number} currency per period
 */
export function equivalentAnnualWorth(cashflows, i) {
  requireSeries('cashflows', cashflows, { minLength: 2 }); rate(i);
  let pv = 0;
  for (let k = 0; k < cashflows.length; k++) pv += cashflows[k] * Math.pow(1 + i, -k);
  return pv * factorAP(i, cashflows.length - 1);
}

/**
 * Equivalent annual cost of owning and operating an asset: capital recovery (P - S)(A/P) + S i, plus annual cost.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 6 (capital recovery).
 * @param {{capex: number, annualCost?: number, salvage?: number, i: number, n: number}} o capex first cost
 *   (currency, positive), annualCost end-of-period cost (currency/period), salvage value at n (currency), i rate
 *   per period (fraction), n life (periods, > 0)
 * @returns {{capitalRecovery: number, annualCost: number, eac: number}} currency per period
 */
export function equivalentAnnualCost({ capex, annualCost = 0, salvage = 0, i, n }) {
  requireNumber('capex', capex); requireNumber('annualCost', annualCost); requireNumber('salvage', salvage);
  const capitalRecovery = (capex - salvage) * factorAP(i, n) + salvage * i;
  return { capitalRecovery, annualCost, eac: capitalRecovery + annualCost };
}

// ---------------------------------------------------------------------------------------------------------------
// Inflation and discounting conventions
// ---------------------------------------------------------------------------------------------------------------

/**
 * Nominal (market) rate from a real rate and inflation, exact Fisher relation (1 + in) = (1 + ir)(1 + f).
 * Ref: Fisher (1930), The Theory of Interest; Blank & Tarquin, Engineering Economy, ch. 14.
 * @param {number} real real rate (fraction, > -1)
 * @param {number} inflation inflation rate (fraction, > -1)
 * @returns {number} nominal rate (fraction)
 */
export function realToNominal(real, inflation) {
  requireNumber('real', real, { gt: -1 }); requireNumber('inflation', inflation, { gt: -1 });
  return (1 + real) * (1 + inflation) - 1;
}

/**
 * Real rate from a nominal rate and inflation, exact Fisher relation ir = (1 + in)/(1 + f) - 1.
 * Ref: Fisher (1930), The Theory of Interest; Blank & Tarquin, Engineering Economy, ch. 14.
 * @param {number} nominal nominal rate (fraction, > -1)
 * @param {number} inflation inflation rate (fraction, > -1)
 * @returns {number} real rate (fraction)
 */
export function nominalToReal(nominal, inflation) {
  requireNumber('nominal', nominal, { gt: -1 }); requireNumber('inflation', inflation, { gt: -1 });
  return (1 + nominal) / (1 + inflation) - 1;
}

/**
 * Convert an amount in base-date (constant) money to then-current (nominal) money t periods later.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 14.
 * @param {number} amount constant-money amount (currency of base date)
 * @param {number} inflation inflation rate per period (fraction, > -1)
 * @param {number} t periods after the base date (may be negative)
 * @returns {number} nominal amount (currency of date t)
 */
export function inflate(amount, inflation, t) {
  requireNumber('amount', amount); requireNumber('inflation', inflation, { gt: -1 }); requireNumber('t', t);
  return amount * Math.pow(1 + inflation, t);
}

/**
 * Convert a nominal amount at period t to base-date (constant) money.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 14.
 * @param {number} amount nominal amount (currency of date t)
 * @param {number} inflation inflation rate per period (fraction, > -1)
 * @param {number} t periods after the base date
 * @returns {number} constant-money amount
 */
export function deflate(amount, inflation, t) { return inflate(amount, inflation, -t); }

/**
 * Discount factor for a flow in period k. End-of-period convention: (1+i)^-k. Mid-period convention: flows of
 * period k >= 1 are taken at k - 0.5 (period 0 stays at time 0).
 * Ref: SPE PRMS guidance on mid-year discounting; Brealey, Myers & Allen ch. 2.
 * @param {number} i rate per period (fraction, > -1)
 * @param {number} k period index (>= 0)
 * @param {{midYear?: boolean}} [o]
 * @returns {number} dimensionless
 */
export function discountFactor(i, k, { midYear = false } = {}) {
  rate(i); requireNumber('k', k, { min: 0 });
  const t = midYear && k > 0 ? k - 0.5 : k;
  return Math.pow(1 + i, -t);
}

/**
 * Discount factors for periods 0..n-1 under a constant rate or a time-varying rate vector. With a vector,
 * rates[k] applies between period k-1 and k (rates[0] is unused); mid-period flows take half of period k's rate
 * compounding: DF_k(mid) = DF_(k-1) (1 + r_k)^-0.5.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 3 (term structure).
 * @param {number | number[]} rates rate per period (fraction) or one rate per period index
 * @param {number} n number of periods (>= 1)
 * @param {{midYear?: boolean}} [o]
 * @returns {Float64Array} dimensionless factors, length n
 */
export function discountFactors(rates, n, { midYear = false } = {}) {
  requireNumber('n', n, { min: 1, integer: true });
  const df = new Float64Array(n);
  const vec = typeof rates !== 'number';
  if (vec) {
    requireSeries('rates', rates);
    if (rates.length < n) throw new RangeError('rates must have one entry per period');
  } else rate(rates);
  let end = 1;
  df[0] = 1;
  for (let k = 1; k < n; k++) {
    const r = vec ? rates[k] : rates;
    if (!(r > -1)) throw new RangeError('rates must be > -1');
    df[k] = midYear ? end / Math.sqrt(1 + r) : end / (1 + r);
    end /= 1 + r;
  }
  return df;
}
