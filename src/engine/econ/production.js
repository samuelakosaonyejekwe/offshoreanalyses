// Production and revenue volume models: Arps decline curves, plateau/decline field profiles, multi-well build-up
// under a facility constraint, associated water and gas streams, gas sales shrinkage, renewable generation
// profiles and tariff/throughput services.
// Rates are per time unit (the same unit as 1/Di, normally years); volumes are rate × time.

import { levenbergMarquardt } from '../math/optim.js';
import { requireNumber, requireSeries } from './tvm.js';
import { normInv } from './risk.js';

const arpsArgs = (qi, Di, b) => {
  requireNumber('qi', qi, { gt: 0 }); requireNumber('Di', Di, { gt: 0 }); requireNumber('b', b, { min: 0, max: 2 });
};
const EXP_TOL = 1e-9;

/**
 * Arps decline rate q(t): exponential (b = 0) q = qi e^(−Di t); hyperbolic (0 < b < 1, extended to b <= 2)
 * q = qi (1 + b Di t)^(−1/b); harmonic (b = 1) q = qi / (1 + Di t).
 * Ref: Arps, J.J. (1945), Analysis of decline curves, Trans. AIME 160, 228–247.
 * @param {{qi: number, Di: number, b?: number, t: number}} o initial rate (volume/time), nominal initial decline
 *   (1/time), decline exponent (-), time (same unit as 1/Di, >= 0)
 * @returns {number} rate (volume/time)
 */
export function arpsRate({ qi, Di, b = 0, t }) {
  arpsArgs(qi, Di, b); requireNumber('t', t, { min: 0 });
  return b < EXP_TOL ? qi * Math.exp(-Di * t) : qi * Math.pow(1 + b * Di * t, -1 / b);
}

/**
 * Arps cumulative production to time t. Exponential: (qi − q)/Di; hyperbolic: qi^b [qi^(1−b) − q^(1−b)] /
 * [(1 − b) Di]; harmonic: (qi/Di) ln(qi/q).
 * Ref: Arps (1945), Trans. AIME 160.
 * @param {{qi: number, Di: number, b?: number, t: number}} o as for arpsRate
 * @returns {number} cumulative volume
 */
export function arpsCumulative({ qi, Di, b = 0, t }) {
  const q = arpsRate({ qi, Di, b, t });
  if (b < EXP_TOL) return (qi - q) / Di;
  if (Math.abs(b - 1) < EXP_TOL) return qi / Di * Math.log(qi / q);
  return Math.pow(qi, b) / ((1 - b) * Di) * (Math.pow(qi, 1 - b) - Math.pow(q, 1 - b));
}

/**
 * Time for the rate to fall to a limit (economic limit): exponential ln(qi/q)/Di; otherwise [(qi/q)^b − 1]/(b Di).
 * Ref: Arps (1945), Trans. AIME 160.
 * @param {{qi: number, Di: number, b?: number, qLimit: number}} o qLimit in (0, qi]
 * @returns {number} time (same unit as 1/Di)
 */
export function arpsTimeToRate({ qi, Di, b = 0, qLimit }) {
  arpsArgs(qi, Di, b); requireNumber('qLimit', qLimit, { gt: 0, max: qi });
  return b < EXP_TOL ? Math.log(qi / qLimit) / Di : (Math.pow(qi / qLimit, b) - 1) / (b * Di);
}

/**
 * Estimated ultimate recovery of the decline segment down to an economic-limit rate, with the time to reach it.
 * For b >= 1 the cumulative is unbounded as the limit tends to zero, so a positive limit is required.
 * Ref: Arps (1945), Trans. AIME 160.
 * @param {{qi: number, Di: number, b?: number, qLimit: number}} o
 * @returns {{eur: number, timeToLimit: number}} volume and time
 */
export function arpsEur({ qi, Di, b = 0, qLimit }) {
  const t = arpsTimeToRate({ qi, Di, b, qLimit });
  return { eur: arpsCumulative({ qi, Di, b, t }), timeToLimit: t };
}

/**
 * Conversion between the nominal (instantaneous) decline Di and the effective decline over one time unit
 * De = (qi − q(1))/qi. For b = 0: De = 1 − e^(−Di); generally De = 1 − (1 + b Di)^(−1/b).
 * Ref: Arps (1945), Trans. AIME 160 (nominal and effective decline).
 * @param {{Di?: number, De?: number, b?: number}} o give one of Di (1/time) or De (fraction per time unit, 0–1)
 * @returns {{Di: number, De: number}}
 */
export function declineConversion({ Di, De, b = 0 }) {
  requireNumber('b', b, { min: 0, max: 2 });
  if (Di !== undefined) {
    requireNumber('Di', Di, { gt: 0 });
    return { Di, De: b < EXP_TOL ? 1 - Math.exp(-Di) : 1 - Math.pow(1 + b * Di, -1 / b) };
  }
  requireNumber('De', De, { gt: 0, lt: 1 });
  return { Di: b < EXP_TOL ? -Math.log(1 - De) : (Math.pow(1 - De, -b) - 1) / b, De };
}

/**
 * Fit an Arps decline to rate–time data by nonlinear least squares (Levenberg–Marquardt) on the rates.
 * model 'exponential' fixes b = 0, 'harmonic' fixes b = 1, 'hyperbolic' fits b in [0, bMax].
 * Ref: Arps (1945); Marquardt, D.W. (1963), J. SIAM 11(2).
 * @param {{t: number[], q: number[], model?: 'exponential' | 'hyperbolic' | 'harmonic', bMax?: number}} o times and
 *   positive rates (at least 3 points)
 * @returns {{qi: number, Di: number, b: number, rmse: number, r2: number, converged: boolean}}
 */
export function arpsFit({ t, q, model = 'hyperbolic', bMax = 1 }) {
  requireSeries('t', t, { minLength: 3 }); requireSeries('q', q, { minLength: 3 });
  if (t.length !== q.length) throw new RangeError('t and q must have the same length');
  for (const v of q) if (!(v > 0)) throw new RangeError('q must be positive');
  // starting values from a log-linear (exponential) fit
  const n = t.length; let st = 0; let sy = 0; let stt = 0; let sty = 0;
  for (let k = 0; k < n; k++) { const y = Math.log(q[k]); st += t[k]; sy += y; stt += t[k] * t[k]; sty += t[k] * y; }
  const slope = (n * sty - st * sy) / (n * stt - st * st);
  const D0 = Math.max(1e-6, -slope); const q0 = Math.exp((sy - slope * st) / n);
  const fixedB = model === 'exponential' ? 0 : model === 'harmonic' ? 1 : null;
  if (fixedB === null && model !== 'hyperbolic') throw new RangeError("model must be 'exponential', 'hyperbolic' or 'harmonic'");
  const pred = (p, tk) => { const b = fixedB ?? p[2]; return b < EXP_TOL ? p[0] * Math.exp(-p[1] * tk) : p[0] * Math.pow(1 + b * p[1] * tk, -1 / b); };
  const res = (p) => Array.from(t, (tk, k) => pred(p, tk) - q[k]);
  const run = (x0) => levenbergMarquardt(res, x0, { bounds: fixedB === null
    ? [[1e-12, Infinity], [1e-9, Infinity], [0, bMax]] : [[1e-12, Infinity], [1e-9, Infinity]], tol: 1e-14, maxIter: 400 });
  let best = null;
  const starts = fixedB === null ? [[q0, D0, 0.5 * bMax], [q0, 1.5 * D0, 0.9 * bMax], [q0, D0, 0.1 * bMax]] : [[q0, D0], [q0, 2 * D0]];
  for (const s of starts) { const r = run(s); if (!best || r.cost < best.cost) best = r; }
  const p = best.x; let sse = 0; let sst = 0; const mean = q.reduce((a, v) => a + v, 0) / n;
  for (let k = 0; k < n; k++) { sse += (pred(p, t[k]) - q[k]) ** 2; sst += (q[k] - mean) ** 2; }
  return { qi: p[0], Di: p[1], b: fixedB ?? p[2], rmse: Math.sqrt(sse / n), r2: sst > 0 ? 1 - sse / sst : NaN, converged: !!best.converged };
}

/**
 * Field production profile: linear ramp-up to plateau, plateau, then Arps decline to an economic limit; returns
 * volumes per year by exact integration of each segment within each year.
 * Ref: plateau-and-decline profile construction, Jahn, Cook & Graham, Hydrocarbon Exploration and Production,
 * 2nd ed.; Arps (1945) for the decline.
 * @param {{plateauRate: number, rampYears?: number, plateauYears: number, Di: number, b?: number, years: number,
 *   qLimit?: number, uptime?: number}} o plateau rate (volume/year), durations (years), nominal decline (1/year),
 *   horizon (integer years), economic-limit rate (volume/year), uptime fraction applied to all volumes
 * @returns {{annual: number[], cumulative: number[], endRate: number[], total: number, lastEconomicYear: number}}
 *   annual volumes; endRate is the instantaneous rate at each year end (after uptime)
 */
export function plateauDeclineProfile({ plateauRate, rampYears = 0, plateauYears, Di, b = 0, years, qLimit = 0, uptime = 1 }) {
  requireNumber('plateauRate', plateauRate, { gt: 0 }); requireNumber('rampYears', rampYears, { min: 0 });
  requireNumber('plateauYears', plateauYears, { min: 0 }); requireNumber('years', years, { min: 1, integer: true });
  requireNumber('uptime', uptime, { min: 0, max: 1 }); requireNumber('qLimit', qLimit, { min: 0, lt: plateauRate });
  arpsArgs(plateauRate, Di, b);
  const tDecl = rampYears + plateauYears;
  const tEnd = qLimit > 0 ? tDecl + arpsTimeToRate({ qi: plateauRate, Di, b, qLimit }) : Infinity;
  const cum = (t) => {                                          // cumulative volume to time t
    const tt = Math.min(t, tEnd);
    let v = 0;
    if (rampYears > 0) { const tr = Math.min(tt, rampYears); v += 0.5 * plateauRate * tr * tr / rampYears; }
    if (tt > rampYears) v += plateauRate * (Math.min(tt, tDecl) - rampYears);
    if (tt > tDecl) v += arpsCumulative({ qi: plateauRate, Di, b, t: tt - tDecl });
    return v;
  };
  const rateAt = (t) => (t >= tEnd ? 0 : t < rampYears ? plateauRate * t / rampYears : t <= tDecl ? plateauRate
    : arpsRate({ qi: plateauRate, Di, b, t: t - tDecl }));
  const annual = []; const cumulative = []; const endRate = []; let last = 0;
  for (let y = 1; y <= years; y++) {
    const v = (cum(y) - cum(y - 1)) * uptime;
    annual.push(v); cumulative.push(cum(y) * uptime); endRate.push(rateAt(y) * uptime);
    if (v > 0) last = y;
  }
  return { annual, cumulative, endRate, total: cumulative[years - 1], lastEconomicYear: last };
}

/**
 * Multi-well phased development against a facility capacity. Each well follows an Arps decline expressed on its
 * own cumulative production (rate–cumulative form), so production curtailed by the facility limit is deferred,
 * not lost: when total potential exceeds capacity all wells are choked back pro rata.
 * Rate–cumulative forms: b = 0: q = qi − Di Np; b = 1: q = qi exp(−Di Np/qi); otherwise
 * q = [qi^(1−b) − (1 − b) Di Np / qi^b]^(1/(1−b)).
 * Ref: Arps (1945) rate–cumulative relations; constrained-profile construction as in Jahn, Cook & Graham.
 * @param {{wellsOnline: number[], well: {qi: number, Di: number, b?: number}, facilityCapacity?: number,
 *   uptime?: number, stepsPerYear?: number}} o wellsOnline[y] = wells started in year y (1 per entry), well type
 *   curve (volume/year, 1/year), capacity (volume/year), uptime fraction, integration steps per year (default 24)
 * @returns {{potential: number[], constrained: number[], deferred: number[], wellCount: number[],
 *   cumulative: number[], capacityUtilisation: number[]}} annual volumes; deferred = potential − constrained
 */
export function multiWellProfile({ wellsOnline, well, facilityCapacity = Infinity, uptime = 1, stepsPerYear = 24 }) {
  requireSeries('wellsOnline', wellsOnline); requireNumber('uptime', uptime, { min: 0, max: 1 });
  requireNumber('stepsPerYear', stepsPerYear, { min: 1, integer: true });
  const { qi, Di, b = 0 } = well || {};
  arpsArgs(qi, Di, b);
  if (!(facilityCapacity > 0)) throw new RangeError('facilityCapacity must be > 0');
  const rateOfCum = (Np) => {
    if (b < EXP_TOL) return Math.max(0, qi - Di * Np);
    if (Math.abs(b - 1) < EXP_TOL) return qi * Math.exp(-Di * Np / qi);
    const base = Math.pow(qi, 1 - b) - (1 - b) * Di * Np / Math.pow(qi, b);
    return b < 1 ? (base > 0 ? Math.pow(base, 1 / (1 - b)) : 0) : Math.pow(base, 1 / (1 - b));
  };
  const years = wellsOnline.length; const dt = 1 / stepsPerYear;
  const groups = [];                 // wells started in the same year share one cumulative
  const potential = []; const constrained = []; const deferred = []; const wellCount = []; const cumulative = [];
  const capacityUtilisation = [];
  // unconstrained reference: each group declines on its own time
  const free = [];
  let total = 0; let count = 0;
  for (let y = 0; y < years; y++) {
    const add = requireNumber('wellsOnline', wellsOnline[y], { min: 0 });
    if (add > 0) { groups.push({ n: add, Np: 0 }); free.push({ n: add, Np: 0 }); count += add; }
    let pot = 0; let con = 0;
    for (let s = 0; s < stepsPerYear; s++) {
      let qFree = 0;
      for (const g of free) { const q = rateOfCum(g.Np) * uptime; qFree += g.n * q; g.Np += q * dt; }
      let qPot = 0;
      for (const g of groups) { g.q = rateOfCum(g.Np) * uptime; qPot += g.n * g.q; }
      const f = qPot > facilityCapacity ? facilityCapacity / qPot : 1;
      for (const g of groups) g.Np += g.q * f * dt;
      pot += qFree * dt; con += qPot * f * dt;
    }
    total += con;
    potential.push(pot); constrained.push(con); deferred.push(pot - con); wellCount.push(count); cumulative.push(total);
    capacityUtilisation.push(Number.isFinite(facilityCapacity) ? con / facilityCapacity : NaN);
  }
  return { potential, constrained, deferred, wellCount, cumulative, capacityUtilisation };
}

/**
 * Associated streams from an oil profile with water-cut and gas–oil-ratio trends given as tables per year:
 * water = oil × wc / (1 − wc); gas = oil × GOR; liquids = oil + water.
 * Ref: definitions of water cut and GOR, SPE Petroleum Engineering Handbook, vol. I.
 * @param {{oil: number[], waterCut?: number | number[], gor?: number | number[]}} o oil volume per year, water cut
 *   as a fraction of total liquid (0 <= wc < 1), GOR in gas volume per oil volume (consistent units)
 * @returns {{oil: number[], water: number[], gas: number[], liquids: number[]}} per year
 */
export function associatedStreams({ oil, waterCut = 0, gor = 0 }) {
  requireSeries('oil', oil);
  const at = (s, k) => (typeof s === 'number' ? s : s[Math.min(k, s.length - 1)]);
  const water = []; const gas = []; const liquids = [];
  for (let k = 0; k < oil.length; k++) {
    const wc = requireNumber('waterCut', at(waterCut, k), { min: 0, lt: 1 });
    const g = requireNumber('gor', at(gor, k), { min: 0 });
    const w = oil[k] * wc / (1 - wc);
    water.push(w); gas.push(oil[k] * g); liquids.push(oil[k] + w);
  }
  return { oil: Array.from(oil), water, gas, liquids };
}

/**
 * Gas sales volume after fuel, flare and processing shrinkage: sales = gross × (1 − fuel − flare) × (1 − shrinkage).
 * Ref: gas balance accounting (wellhead to sales), SPE Petroleum Engineering Handbook, vol. III.
 * @param {{gross: number[], fuelFraction?: number, flareFraction?: number, shrinkage?: number}} o gross gas per
 *   year (volume), fractions of gross (0–1), processing shrinkage (0–1)
 * @returns {{sales: number[], fuel: number[], flare: number[], shrink: number[]}} per year
 */
export function gasSales({ gross, fuelFraction = 0, flareFraction = 0, shrinkage = 0 }) {
  requireSeries('gross', gross); requireNumber('fuelFraction', fuelFraction, { min: 0, max: 1 });
  requireNumber('flareFraction', flareFraction, { min: 0, max: 1 }); requireNumber('shrinkage', shrinkage, { min: 0, max: 1 });
  if (fuelFraction + flareFraction > 1) throw new RangeError('fuelFraction + flareFraction must be <= 1');
  const sales = []; const fuel = []; const flare = []; const shrink = [];
  for (const g of gross) {
    const net = g * (1 - fuelFraction - flareFraction);
    fuel.push(g * fuelFraction); flare.push(g * flareFraction); shrink.push(net * shrinkage); sales.push(net * (1 - shrinkage));
  }
  return { sales, fuel, flare, shrink };
}

/**
 * Renewable generation profile by year: energy = gross AEP × availability × (1 − electrical losses) ×
 * (1 − curtailment) × (1 − degradation)^(year − 1), with P50 and exceedance cases from a normal distribution of
 * annual yield (Pxx = P50 (1 − z_xx σ), z_90 = 1.2816; P90 is the energy exceeded with 90 % probability) and an
 * optional monthly shape.
 * Ref: energy-yield assessment practice (normally distributed long-term net yield); exceedance convention as
 * used in energy-yield reports.
 * @param {{capacityMW?: number, capacityFactor?: number, grossAep_MWh?: number, availability?: number,
 *   lossesFrac?: number, curtailment?: number | number[], degradation?: number, years: number, sigma?: number,
 *   exceedance?: number[], monthlyShape?: number[]}} o either grossAep_MWh or capacityMW × capacityFactor × 8760;
 *   fractions 0–1; degradation per year; sigma = relative standard deviation of long-term net yield (default 0);
 *   exceedance levels in percent (default [50, 75, 90, 99]); monthlyShape = 12 weights
 * @returns {{p50: number[], cases: Record<string, number[]>, netCapacityFactor: number | null,
 *   monthlyP50Year1: number[] | null, convention: string}} MWh per year
 */
export function renewableGeneration({ capacityMW, capacityFactor, grossAep_MWh, availability = 1, lossesFrac = 0,
  curtailment = 0, degradation = 0, years, sigma = 0, exceedance = [50, 75, 90, 99], monthlyShape }) {
  requireNumber('years', years, { min: 1, integer: true }); requireNumber('availability', availability, { min: 0, max: 1 });
  requireNumber('lossesFrac', lossesFrac, { min: 0, max: 1 }); requireNumber('degradation', degradation, { min: 0, lt: 1 });
  requireNumber('sigma', sigma, { min: 0, lt: 1 });
  const gross = grossAep_MWh ?? requireNumber('capacityMW', capacityMW, { gt: 0 }) * requireNumber('capacityFactor', capacityFactor, { gt: 0, max: 1 }) * 8760;
  requireNumber('grossAep_MWh', gross, { min: 0 });
  const p50 = [];
  for (let y = 0; y < years; y++) {
    const c = requireNumber('curtailment', typeof curtailment === 'number' ? curtailment : curtailment[Math.min(y, curtailment.length - 1)], { min: 0, max: 1 });
    p50.push(gross * availability * (1 - lossesFrac) * (1 - c) * Math.pow(1 - degradation, y));
  }
  const cases = {};
  for (const p of exceedance) {
    requireNumber('exceedance', p, { gt: 0, lt: 100 });
    const f = Math.max(0, 1 - normInv(p / 100) * sigma);
    cases[`P${p}`] = p50.map((v) => v * f);
  }
  let monthly = null;
  if (monthlyShape) {
    if (monthlyShape.length !== 12) throw new RangeError('monthlyShape must have 12 values');
    const s = monthlyShape.reduce((a, v) => a + requireNumber('monthlyShape', v, { min: 0 }), 0);
    if (!(s > 0)) throw new RangeError('monthlyShape must not be all zero');
    monthly = monthlyShape.map((w) => p50[0] * w / s);
  }
  return { p50, cases, netCapacityFactor: capacityMW ? p50[0] / (capacityMW * 8760) : null, monthlyP50Year1: monthly,
    convention: 'Exceedance: Pxx is the annual energy exceeded with xx % probability (P90 < P50).' };
}

/**
 * Tariff / throughput service revenue (port throughput, pipeline transport): revenue = max(volume, take-or-pay
 * volume) × tariff, with tariff escalation and a capacity cap on volume.
 * Ref: ship-or-pay / take-or-pay tariff mechanics as used in transportation agreements.
 * @param {{volumes: number[], tariff: number, escalation?: number, capacity?: number, takeOrPayFraction?: number,
 *   contractedVolume?: number}} o volumes per year (units), tariff (currency/unit, year-1 value), escalation
 *   (fraction/year), capacity (units/year), take-or-pay as a fraction of contractedVolume
 * @returns {{billedVolume: number[], revenue: number[], tariff: number[], shortfallPayments: number[]}}
 */
export function throughputRevenue({ volumes, tariff, escalation = 0, capacity = Infinity, takeOrPayFraction = 0, contractedVolume = 0 }) {
  requireSeries('volumes', volumes); requireNumber('tariff', tariff, { min: 0 }); requireNumber('escalation', escalation, { gt: -1 });
  requireNumber('takeOrPayFraction', takeOrPayFraction, { min: 0, max: 1 });
  const billedVolume = []; const revenue = []; const tar = []; const shortfallPayments = [];
  const floor = takeOrPayFraction * contractedVolume;
  volumes.forEach((v, k) => {
    const vol = Math.min(requireNumber('volumes', v, { min: 0 }), capacity);
    const t = tariff * Math.pow(1 + escalation, k);
    const billed = Math.max(vol, floor);
    billedVolume.push(billed); tar.push(t); revenue.push(billed * t); shortfallPayments.push((billed - vol) * t);
  });
  return { billedVolume, revenue, tariff: tar, shortfallPayments };
}

/**
 * Ramp-up factors for the first years of operation (learning to nameplate): linear from `start` to 1 over
 * `rampYears`, by year-average.
 * Ref: commissioning ramp-up allowance; generic linear model.
 * @param {{years: number, rampYears: number, start?: number}} o counts in years; start fraction 0–1
 * @returns {number[]} multiplier per year (0–1)
 */
export function rampUpFactors({ years, rampYears, start = 0 }) {
  requireNumber('years', years, { min: 1, integer: true }); requireNumber('rampYears', rampYears, { min: 0 });
  requireNumber('start', start, { min: 0, max: 1 });
  const f = (t) => (rampYears <= 0 || t >= rampYears ? 1 : start + (1 - start) * t / rampYears);
  const out = [];
  for (let y = 0; y < years; y++) {
    // exact average of the piecewise-linear factor over [y, y + 1]
    const a = Math.min(y, rampYears); const b = Math.min(y + 1, rampYears);
    const rampPart = rampYears > 0 ? (b - a) * 0.5 * (f(a) + (b >= rampYears ? 1 : f(b))) : 0;
    out.push(rampPart + (y + 1 - Math.max(y, Math.min(y + 1, rampYears))));
  }
  return out;
}
