// Tax and book depreciation schedules and the present value of the depreciation tax shield.
// All schedules return annual charges in the currency of the basis; year index 0 is the first year of service.

import { requireNumber, requireSeries } from './tvm.js';

/**
 * MACRS general depreciation system percentage tables, half-year convention (percent of unadjusted basis per
 * recovery year). 3/5/7/10-year classes use 200 % declining balance, 15/20-year classes 150 %, switching to
 * straight line. Values are the published table values, not recomputed.
 * Ref: IRS Publication 946, Appendix A, Table A-1.
 */
export const MACRS_HALF_YEAR = Object.freeze({
  3: Object.freeze([33.33, 44.45, 14.81, 7.41]),
  5: Object.freeze([20.00, 32.00, 19.20, 11.52, 11.52, 5.76]),
  7: Object.freeze([14.29, 24.49, 17.49, 12.49, 8.93, 8.92, 8.93, 4.46]),
  10: Object.freeze([10.00, 18.00, 14.40, 11.52, 9.22, 7.37, 6.55, 6.55, 6.56, 6.55, 3.28]),
  15: Object.freeze([5.00, 9.50, 8.55, 7.70, 6.93, 6.23, 5.90, 5.90, 5.91, 5.90, 5.91, 5.90, 5.91, 5.90, 5.91, 2.95]),
  20: Object.freeze([3.750, 7.219, 6.677, 6.177, 5.713, 5.285, 4.888, 4.522, 4.462, 4.461, 4.462, 4.461, 4.462, 4.461,
    4.462, 4.461, 4.462, 4.461, 4.462, 4.461, 2.231]),
});

function finish(charges, basis) {
  const bookValues = new Array(charges.length);
  let bv = basis;
  for (let k = 0; k < charges.length; k++) { bv -= charges[k]; bookValues[k] = bv; }
  return { charges, bookValues, total: basis - bv };
}

/**
 * Straight-line depreciation: (basis - salvage) / life each year.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 16.
 * @param {{basis: number, salvage?: number, life: number}} o basis and salvage in currency, life in years (integer >= 1)
 * @returns {{charges: number[], bookValues: number[], total: number}} currency per year, end-of-year book values
 */
export function straightLine({ basis, salvage = 0, life }) {
  requireNumber('basis', basis, { min: 0 }); requireNumber('salvage', salvage, { min: 0 });
  requireNumber('life', life, { min: 1, integer: true });
  const d = (basis - salvage) / life;
  return finish(new Array(life).fill(d), basis);
}

/**
 * Declining-balance depreciation at rate factor/life on the book value (factor 2 = double declining balance),
 * never depreciating below salvage. With switchToSL the charge switches to straight line over the remaining life
 * as soon as that is larger, so the asset is fully written down to salvage by the end of its life.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 16 (DB, DDB and switching).
 * @param {{basis: number, salvage?: number, life: number, factor?: number, switchToSL?: boolean}} o currency,
 *   currency, years (integer), DB multiple (0 < factor <= life), switch flag (default true)
 * @returns {{charges: number[], bookValues: number[], total: number, switchYear: number | null}} switchYear is
 *   the 1-based year in which straight line takes over (null if never)
 */
export function decliningBalance({ basis, salvage = 0, life, factor = 2, switchToSL = true }) {
  requireNumber('basis', basis, { min: 0 }); requireNumber('salvage', salvage, { min: 0 });
  requireNumber('life', life, { min: 1, integer: true }); requireNumber('factor', factor, { gt: 0, max: life });
  const d = factor / life;
  const charges = [];
  let bv = basis;
  let switchYear = null;
  for (let k = 0; k < life; k++) {
    let c = Math.max(0, Math.min(bv * d, bv - salvage));
    if (switchToSL) {
      const sl = Math.max(0, (bv - salvage) / (life - k));
      if (sl > c + 1e-12 * Math.max(1, basis) || switchYear !== null) {
        if (switchYear === null) switchYear = k + 1;
        c = sl;
      }
    }
    charges.push(c);
    bv -= c;
  }
  return { ...finish(charges, basis), switchYear };
}

/**
 * Sum-of-years-digits depreciation: year k (1-based) charge = (life - k + 1) / [life (life + 1) / 2] (basis - salvage).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 16.
 * @param {{basis: number, salvage?: number, life: number}} o currency, currency, years (integer >= 1)
 * @returns {{charges: number[], bookValues: number[], total: number}}
 */
export function sumOfYearsDigits({ basis, salvage = 0, life }) {
  requireNumber('basis', basis, { min: 0 }); requireNumber('salvage', salvage, { min: 0 });
  requireNumber('life', life, { min: 1, integer: true });
  const syd = life * (life + 1) / 2;
  const charges = [];
  for (let k = 1; k <= life; k++) charges.push((life - k + 1) / syd * (basis - salvage));
  return finish(charges, basis);
}

/**
 * Units-of-production depreciation: charge in a year = (basis - salvage) × units produced / total units.
 * totalUnits defaults to the sum of the supplied production (e.g. reserves, operating hours).
 * Ref: Blank & Tarquin, Engineering Economy, ch. 16 (unit-of-production); IAS 16 §62.
 * @param {{basis: number, salvage?: number, units: number[], totalUnits?: number}} o currency, currency,
 *   production per year (any unit), total recoverable units (same unit)
 * @returns {{charges: number[], bookValues: number[], total: number}}
 */
export function unitsOfProduction({ basis, salvage = 0, units, totalUnits }) {
  requireNumber('basis', basis, { min: 0 }); requireNumber('salvage', salvage, { min: 0 });
  requireSeries('units', units);
  let sum = 0;
  for (const u of units) { if (u < 0) throw new RangeError('units must be >= 0'); sum += u; }
  const tot = totalUnits ?? sum;
  requireNumber('totalUnits', tot, { gt: 0 });
  let left = basis - salvage;
  const charges = [];
  for (const u of units) { const c = Math.min(left, (basis - salvage) * u / tot); charges.push(c); left -= c; }
  return finish(charges, basis);
}

/**
 * MACRS (GDS, half-year convention) depreciation from the IRS percentage table; salvage is ignored by MACRS.
 * Ref: IRS Publication 946, Appendix A, Table A-1.
 * @param {{basis: number, recoveryPeriod: 3 | 5 | 7 | 10 | 15 | 20}} o currency, property class in years
 * @returns {{charges: number[], bookValues: number[], total: number, percentages: number[]}} recoveryPeriod + 1 years
 */
export function macrs({ basis, recoveryPeriod }) {
  requireNumber('basis', basis, { min: 0 });
  const pct = MACRS_HALF_YEAR[recoveryPeriod];
  if (!pct) throw new RangeError('recoveryPeriod must be one of 3, 5, 7, 10, 15, 20');
  return { ...finish(pct.map((p) => basis * p / 100), basis), percentages: pct.slice() };
}

/**
 * Writing-down allowance on a reducing-balance pool (UK capital-allowance style): each year the allowance is
 * rate × pool after additions; an optional first-year allowance is taken on each addition in its year and only
 * the remainder enters the pool. The pool is never fully exhausted — the residual is reported.
 * Ref: UK Capital Allowances Act 2001, Part 2 (main/special-rate pools); the rate is an argument because it is
 * set by the Finance Act in force.
 * @param {{additions: number[], rate: number, firstYearAllowance?: number, years?: number}} o additions = qualifying
 *   expenditure per year (currency), rate = WDA fraction per year (0–1), firstYearAllowance = fraction of each
 *   addition deducted immediately (0–1), years = schedule length (default additions.length)
 * @returns {{charges: number[], pool: number[], total: number, residualPool: number}} allowances per year and
 *   closing pool per year (currency)
 */
export function writingDownAllowance({ additions, rate, firstYearAllowance = 0, years }) {
  requireSeries('additions', additions); requireNumber('rate', rate, { min: 0, max: 1 });
  requireNumber('firstYearAllowance', firstYearAllowance, { min: 0, max: 1 });
  const n = years ?? additions.length;
  requireNumber('years', n, { min: 1, integer: true });
  const charges = []; const pool = [];
  let p = 0; let total = 0;
  for (let k = 0; k < n; k++) {
    const add = k < additions.length ? additions[k] : 0;
    const fya = add * firstYearAllowance;
    p += add - fya;
    const wda = p * rate;
    p -= wda;
    charges.push(fya + wda); pool.push(p); total += fya + wda;
  }
  return { charges, pool, total, residualPool: p };
}

/**
 * Bonus (first-year) depreciation: bonusFraction of the basis is expensed in year 1 and the remainder follows
 * the supplied method.
 * Ref: US IRC §168(k) (bonus fraction is statutory and time-dependent, hence an argument).
 * @param {{basis: number, bonusFraction: number, method?: object}} o currency, fraction 0–1, and a method spec as
 *   accepted by depreciationSchedule (without basis) applied to the remaining basis
 * @returns {{charges: number[], bookValues: number[], total: number, bonus: number}}
 */
export function bonusDepreciation({ basis, bonusFraction, method = { method: 'sl', life: 1 } }) {
  requireNumber('basis', basis, { min: 0 }); requireNumber('bonusFraction', bonusFraction, { min: 0, max: 1 });
  const bonus = basis * bonusFraction;
  const rest = depreciationSchedule({ ...method, basis: basis - bonus });
  const charges = rest.slice();
  if (charges.length === 0) charges.push(0);
  charges[0] += bonus;
  return { ...finish(charges, basis), bonus };
}

/**
 * Dispatcher returning the annual charges for a method specification.
 * method: 'sl' | 'db' | 'ddb' | 'syd' | 'uop' | 'macrs' | 'wda' | 'schedule' | 'none'.
 * 'schedule' uses `fractions` (of basis per year); 'wda' needs `rate` and `years`; 'uop' needs `units`.
 * Ref: methods as cited on the individual functions.
 * @param {object} spec { method, basis (currency), life?, salvage?, factor?, switchToSL?, recoveryPeriod?, rate?,
 *   years?, firstYearAllowance?, units?, totalUnits?, fractions?, bonusFraction? }
 * @returns {number[]} charges per year (currency), year 0 = first year of service
 */
export function depreciationSchedule(spec) {
  if (!spec || typeof spec !== 'object') throw new RangeError('depreciation spec must be an object');
  const { method = 'sl', basis } = spec;
  requireNumber('basis', basis, { min: 0 });
  if (spec.bonusFraction) {
    const { bonusFraction, ...rest } = spec;
    delete rest.basis;
    return bonusDepreciation({ basis, bonusFraction, method: rest }).charges;
  }
  switch (method) {
    case 'none': return [];
    case 'sl': return straightLine(spec).charges;
    case 'db': return decliningBalance({ factor: 1.5, ...spec }).charges;
    case 'ddb': return decliningBalance({ ...spec, factor: 2 }).charges;
    case 'syd': return sumOfYearsDigits(spec).charges;
    case 'uop': return unitsOfProduction(spec).charges;
    case 'macrs': return macrs({ basis, recoveryPeriod: spec.recoveryPeriod ?? spec.life }).charges;
    case 'wda': return writingDownAllowance({ additions: [basis], rate: spec.rate,
      firstYearAllowance: spec.firstYearAllowance ?? 0, years: spec.years ?? spec.life ?? 25 }).charges;
    case 'schedule': {
      requireSeries('fractions', spec.fractions);
      return Array.from(spec.fractions, (f) => basis * f);
    }
    default: throw new RangeError(`method '${method}' is not a known depreciation method`);
  }
}

/**
 * Present value of the depreciation tax shield: Σ taxRate × charge_k / (1 + rate)^(startPeriod + k).
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 6.
 * @param {{charges: number[], taxRate: number, rate: number, startPeriod?: number}} o charges per year (currency),
 *   tax rate (fraction 0–1), discount rate (fraction), period at which the first charge is taken (default 1)
 * @returns {{pv: number, undiscounted: number}} currency
 */
export function taxShieldPV({ charges, taxRate, rate, startPeriod = 1 }) {
  requireSeries('charges', charges); requireNumber('taxRate', taxRate, { min: 0, max: 1 });
  requireNumber('rate', rate, { gt: -1 }); requireNumber('startPeriod', startPeriod, { min: 0 });
  let pv = 0; let sum = 0;
  for (let k = 0; k < charges.length; k++) {
    pv += taxRate * charges[k] * Math.pow(1 + rate, -(startPeriod + k));
    sum += taxRate * charges[k];
  }
  return { pv, undiscounted: sum };
}
