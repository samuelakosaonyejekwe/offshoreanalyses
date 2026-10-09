// Maintenance optimisation: replacement policies, imperfect repair, condition-based maintenance, inspection
// intervals, downtime economics and an offshore-wind operations-and-maintenance simulator.
// Times in any consistent unit unless a name says hours or years; costs in any one currency.

import { gamma, gammaP } from '../math/special.js';
import { makeRng } from '../math/rng.js';
import { renewalFunction } from './lifedata.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function positive(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

// ∫₀ᵀ exp(−(t/η)^β) dt = (η/β)·γ(1/β, (T/η)^β), γ the lower incomplete gamma function.
function weibullMeanUpTo(T, shape, scale) {
  return (scale / shape) * gamma(1 / shape) * gammaP(1 / shape, (T / scale) ** shape);
}

function minimise1d(f, lo, hi, n = 240) {
  // Log-spaced scan followed by golden-section refinement around the best point.
  let best = lo;
  let fBest = Infinity;
  let kBest = 0;
  const r = Math.log(hi / lo) / n;
  const xs = [];
  for (let k = 0; k <= n; k++) {
    const x = lo * Math.exp(r * k);
    xs.push(x);
    const v = f(x);
    if (v < fBest) { fBest = v; best = x; kBest = k; }
  }
  let a = xs[Math.max(0, kBest - 1)];
  let b = xs[Math.min(n, kBest + 1)];
  const g = (Math.sqrt(5) - 1) / 2;
  let c = b - g * (b - a);
  let d = a + g * (b - a);
  let fc = f(c);
  let fd = f(d);
  for (let i = 0; i < 80 && b - a > 1e-10 * b; i++) {
    if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); }
    else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
  }
  const x = 0.5 * (a + b);
  const fx = f(x);
  return fx <= fBest ? { x, f: fx, atUpperEnd: kBest === n } : { x: best, f: fBest, atUpperEnd: kBest === n };
}

/**
 * Age-replacement cost rate for a Weibull life: replace preventively at age T (cost c_p) or on failure (c_f):
 *   C(T) = [c_p·R(T) + c_f·F(T)] / ∫₀ᵀ R(t) dt.
 * Reference: Barlow & Hunter (1960) Operations Research 8, 90–100; Barlow & Proschan, Mathematical Theory of
 * Reliability (1965).
 * @param {object} o
 * @param {number} o.interval replacement age T [time]
 * @param {number} o.shape Weibull shape β [-]
 * @param {number} o.scale Weibull scale η [time]
 * @param {number} o.costPreventive c_p [currency]
 * @param {number} o.costFailure c_f [currency]
 * @returns {{costRate:number, preventiveShare:number, correctiveShare:number, meanCycle:number,
 *   reliabilityAtInterval:number}} costRate [currency/time]; the two shares add up to the cost rate
 */
export function ageReplacementCostRate({ interval, shape, scale, costPreventive, costFailure }) {
  positive(interval, 'interval');
  positive(shape, 'shape');
  positive(scale, 'scale');
  need(costPreventive >= 0 && costFailure >= 0, 'costs must be non-negative');
  const R = Math.exp(-((interval / scale) ** shape));
  const cyc = weibullMeanUpTo(interval, shape, scale);
  return {
    costRate: (costPreventive * R + costFailure * (1 - R)) / cyc, preventiveShare: (costPreventive * R) / cyc,
    correctiveShare: (costFailure * (1 - R)) / cyc, meanCycle: cyc, reliabilityAtInterval: R,
  };
}

/**
 * Optimal age-replacement interval for a Weibull life (minimum long-run cost rate). A finite optimum exists
 * only for an increasing hazard (β > 1) and c_f > c_p; otherwise run-to-failure is optimal and the interval is
 * returned as Infinity with cost rate c_f/MTTF.
 * Reference: Barlow & Proschan (1965).
 * @param {object} o
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.costPreventive [currency]
 * @param {number} o.costFailure [currency]
 * @returns {{interval:number, costRate:number, runToFailureCostRate:number, saving:number}} saving = relative
 *   reduction of the cost rate against run-to-failure [-]
 */
export function optimalAgeReplacement({ shape, scale, costPreventive, costFailure }) {
  positive(shape, 'shape');
  positive(scale, 'scale');
  need(costPreventive > 0 && costFailure > 0, 'costs must be positive');
  const rtf = costFailure / (scale * gamma(1 + 1 / shape));
  if (shape <= 1 || costFailure <= costPreventive) return { interval: Infinity, costRate: rtf, runToFailureCostRate: rtf, saving: 0 };
  const f = (T) => ageReplacementCostRate({ interval: T, shape, scale, costPreventive, costFailure }).costRate;
  const r = minimise1d(f, scale * 1e-3, scale * 20);
  if (r.atUpperEnd || r.f >= rtf) return { interval: Infinity, costRate: rtf, runToFailureCostRate: rtf, saving: 0 };
  return { interval: r.x, costRate: r.f, runToFailureCostRate: rtf, saving: 1 - r.f / rtf };
}

/**
 * Block-replacement cost rate: all units are replaced every T regardless of age, failures in between are
 * renewed: C(T) = [c_p + c_f·M(T)]/T with M the renewal function.
 * Reference: Barlow & Proschan (1965).
 * @param {object} o
 * @param {number} o.interval T [time]
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.costPreventive [currency]
 * @param {number} o.costFailure [currency]
 * @param {number} [o.steps=600] discretisation of the renewal equation [-]
 * @returns {{costRate:number, expectedFailures:number}}
 */
export function blockReplacementCostRate({ interval, shape, scale, costPreventive, costFailure, steps = 600 }) {
  positive(interval, 'interval');
  positive(shape, 'shape');
  positive(scale, 'scale');
  const M = renewalFunction({ cdf: (x) => 1 - Math.exp(-((x / scale) ** shape)), tMax: interval, steps }).M[steps];
  return { costRate: (costPreventive + costFailure * M) / interval, expectedFailures: M };
}

/**
 * Optimal block-replacement interval for a Weibull life by evaluating C(T) on the renewal-function grid.
 * @param {object} o
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.costPreventive [currency]
 * @param {number} o.costFailure [currency]
 * @param {number} [o.maxInterval] search limit [time], default 3·scale
 * @param {number} [o.steps=1500]
 * @returns {{interval:number, costRate:number, expectedFailures:number, curve: {t: Float64Array, costRate: Float64Array}}}
 *   interval = Infinity when the cost rate is still falling at the search limit
 */
export function optimalBlockReplacement({ shape, scale, costPreventive, costFailure, maxInterval, steps = 1500 }) {
  positive(shape, 'shape');
  positive(scale, 'scale');
  need(costPreventive > 0 && costFailure > 0, 'costs must be positive');
  const tMax = maxInterval === undefined ? 3 * scale : maxInterval;
  const { t, M } = renewalFunction({ cdf: (x) => 1 - Math.exp(-((x / scale) ** shape)), tMax, steps });
  const costRate = new Float64Array(steps + 1);
  costRate[0] = Infinity;
  let k = 1;
  for (let i = 1; i <= steps; i++) {
    costRate[i] = (costPreventive + costFailure * M[i]) / t[i];
    if (costRate[i] < costRate[k]) k = i;
  }
  const curve = { t, costRate };
  if (k === steps) return { interval: Infinity, costRate: costRate[steps], expectedFailures: M[steps], curve };
  return { interval: t[k], costRate: costRate[k], expectedFailures: M[k], curve };
}

/**
 * Preventive versus corrective cost-rate curves for age replacement over a range of intervals.
 * @param {object} o
 * @param {ArrayLike<number>} o.intervals [time]
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.costPreventive [currency]
 * @param {number} o.costFailure [currency]
 * @returns {{intervals:number[], preventive:number[], corrective:number[], total:number[]}} [currency/time]
 */
export function costRateCurves({ intervals, shape, scale, costPreventive, costFailure }) {
  const out = { intervals: [], preventive: [], corrective: [], total: [] };
  for (let i = 0; i < intervals.length; i++) {
    const r = ageReplacementCostRate({ interval: intervals[i], shape, scale, costPreventive, costFailure });
    out.intervals.push(intervals[i]);
    out.preventive.push(r.preventiveShare);
    out.corrective.push(r.correctiveShare);
    out.total.push(r.costRate);
  }
  return out;
}

/**
 * Imperfect repair by the Kijima type I virtual-age model with a Weibull baseline: after the n-th repair the
 * virtual age is v_n = v_(n−1) + q·x_n (q = 0: as good as new, q = 1: as bad as old). Monte Carlo estimate of the
 * expected number of failures in a horizon.
 * Reference: Kijima (1989) J. Appl. Probab. 26, 89–102.
 * @param {object} o
 * @param {number} o.shape [-]
 * @param {number} o.scale [time]
 * @param {number} o.q repair effectiveness parameter in [0, 1]
 * @param {number} o.horizon [time]
 * @param {number} [o.n=5000] number of histories [-]
 * @param {number} [o.seed=1]
 * @returns {{expectedFailures:number, stdError:number}}
 */
export function kijimaExpectedFailures({ shape, scale, q, horizon, n = 5000, seed = 1 }) {
  positive(shape, 'shape');
  positive(scale, 'scale');
  need(q >= 0 && q <= 1, 'q must be in [0, 1]');
  positive(horizon, 'horizon');
  const rng = makeRng(seed);
  let s = 0;
  let s2 = 0;
  for (let k = 0; k < n; k++) {
    let t = 0;
    let v = 0;
    let c = 0;
    for (;;) {
      // Conditional Weibull given survival to virtual age v.
      const x = scale * ((v / scale) ** shape - Math.log(1 - rng.next())) ** (1 / shape) - v;
      t += x;
      if (t > horizon || c > 1e6) break;
      c++;
      v += q * x;
    }
    s += c;
    s2 += c * c;
  }
  const mean = s / n;
  return { expectedFailures: mean, stdError: Math.sqrt(Math.max(0, s2 / n - mean * mean) / n) };
}

// Marsaglia–Tsang gamma variate with unit scale.
function gammaVariate(rng, k) {
  if (k < 1) return gammaVariate(rng, k + 1) * (1 - rng.next()) ** (1 / k);
  const d = k - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  for (;;) {
    let x;
    let v;
    do {
      x = rng.normal();
      v = 1 + c * x;
    } while (v <= 0);
    v = v * v * v;
    const u = 1 - rng.next();
    if (Math.log(u) < 0.5 * x * x + d - d * v + d * Math.log(v)) return d * v;
  }
}

/**
 * Condition-based maintenance on a stationary gamma degradation process, by Monte Carlo renewal-reward.
 * Degradation X(t) has independent increments ~ Gamma(shapeRate·Δt, scale). The item is inspected every
 * `inspectionInterval`; it is replaced preventively when X >= preventiveThreshold at an inspection and
 * correctively when X >= failureLevel (found at the next inspection unless failures are self-announcing, in which
 * case the crossing time is interpolated within the interval). Cost rate = E[cycle cost]/E[cycle length].
 * Reference: van Noortwijk (2009) Reliability Engineering & System Safety 94, 2–21; Grall, Bérenguer &
 * Dieulle (2002) RESS 76, 167–180.
 * @param {object} o
 * @param {number} o.shapeRate gamma shape per unit time [1/time]
 * @param {number} o.scale gamma scale [degradation unit]; mean rate = shapeRate·scale
 * @param {number} o.failureLevel L [degradation unit]
 * @param {number} o.preventiveThreshold M <= L [degradation unit]
 * @param {number} o.inspectionInterval τ [time]
 * @param {number} o.costInspection [currency]
 * @param {number} o.costPreventive [currency]
 * @param {number} o.costFailure [currency]
 * @param {number} [o.downtimeCostRate=0] cost per unit time spent failed and undetected [currency/time]
 * @param {boolean} [o.selfAnnouncing=false]
 * @param {number} [o.n=4000] number of renewal cycles [-]
 * @param {number} [o.seed=1]
 * @returns {{costRate:number, meanCycleLength:number, failureFraction:number, inspectionsPerCycle:number}}
 */
export function cbmGammaProcess({ shapeRate, scale, failureLevel, preventiveThreshold, inspectionInterval, costInspection, costPreventive,
  costFailure, downtimeCostRate = 0, selfAnnouncing = false, n = 4000, seed = 1 }) {
  positive(shapeRate, 'shapeRate');
  positive(scale, 'scale');
  positive(failureLevel, 'failureLevel');
  need(preventiveThreshold > 0 && preventiveThreshold <= failureLevel, 'preventiveThreshold must be in (0, failureLevel]');
  positive(inspectionInterval, 'inspectionInterval');
  need(costInspection >= 0 && costPreventive >= 0 && costFailure >= 0, 'costs must be non-negative');
  const rng = makeRng(seed);
  const k = shapeRate * inspectionInterval;
  let cost = 0;
  let len = 0;
  let fails = 0;
  let insp = 0;
  for (let c = 0; c < n; c++) {
    let x = 0;
    let t = 0;
    for (;;) {
      const dx = scale * gammaVariate(rng, k);
      const prev = x;
      x += dx;
      t += inspectionInterval;
      if (x >= failureLevel) {
        fails++;
        const frac = dx > 0 ? (failureLevel - prev) / dx : 1;
        if (selfAnnouncing) {
          t -= inspectionInterval * (1 - frac);
        } else {
          cost += costInspection + downtimeCostRate * inspectionInterval * (1 - frac);
          insp++;
        }
        cost += costFailure;
        break;
      }
      cost += costInspection;
      insp++;
      if (x >= preventiveThreshold) { cost += costPreventive; break; }
    }
    len += t;
  }
  return { costRate: cost / len, meanCycleLength: len / n, failureFraction: fails / n, inspectionsPerCycle: insp / n };
}

/**
 * Optimise the preventive threshold (and optionally the inspection interval) of {@link cbmGammaProcess} on a grid
 * with common random numbers.
 * @param {object} o arguments of cbmGammaProcess without preventiveThreshold, plus:
 * @param {number[]} [o.thresholds] candidate thresholds; default 20 values from 5 % to 100 % of failureLevel
 * @param {number[]} [o.intervals] candidate inspection intervals; default [o.inspectionInterval]
 * @returns {{preventiveThreshold:number, inspectionInterval:number, costRate:number, grid: {threshold:number,
 *   interval:number, costRate:number, failureFraction:number}[]}}
 */
export function optimiseCbmThreshold(o) {
  const thresholds = o.thresholds || Array.from({ length: 20 }, (_, i) => (o.failureLevel * (i + 1)) / 20);
  const intervals = o.intervals || [o.inspectionInterval];
  let best = null;
  const grid = [];
  for (const interval of intervals) {
    for (const threshold of thresholds) {
      const r = cbmGammaProcess({ ...o, preventiveThreshold: threshold, inspectionInterval: interval });
      const row = { threshold, interval, costRate: r.costRate, failureFraction: r.failureFraction };
      grid.push(row);
      if (!best || row.costRate < best.costRate) best = row;
    }
  }
  return { preventiveThreshold: best.threshold, inspectionInterval: best.interval, costRate: best.costRate, grid };
}

function integrate(f, a, b, n = 400) {
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += (i % 2 === 1 ? 4 : 2) * f(a + i * h);
  return (s * h) / 3;
}

/**
 * Delay-time inspection model (perfect inspections every T). Defects arise as a Poisson process of rate λ and
 * each becomes a failure after a delay time h with distribution F. The fraction of defects that end as
 * failures is b(T) = (1/T)·∫₀ᵀ F(h) dh; then
 *   cost rate     C(T) = [λT·(b·c_f + (1 − b)·c_r) + c_i]/(T + d_i)
 *   downtime rate D(T) = [λT·b·d_f + d_i]/(T + d_i).
 * Reference: Christer & Waller (1984) J. Oper. Res. Soc. 35, 401–406; Christer (1999) JORS 50, 1120–1137.
 * @param {object} o
 * @param {number} o.interval inspection interval T [time]
 * @param {number} o.defectRate λ [1/time]
 * @param {function(number):number} o.delayCdf delay-time distribution F(h)
 * @param {number} [o.costInspection=0] c_i [currency]
 * @param {number} [o.costFailure=0] c_f, breakdown repair [currency]
 * @param {number} [o.costDefectRepair=0] c_r, repair of a defect found at inspection [currency]
 * @param {number} [o.downtimeInspection=0] d_i [time]
 * @param {number} [o.downtimeFailure=0] d_f [time]
 * @returns {{failureFraction:number, costRate:number, downtimeRate:number, failuresPerInterval:number}}
 */
export function delayTimeModel({ interval, defectRate, delayCdf, costInspection = 0, costFailure = 0, costDefectRepair = 0, downtimeInspection = 0, downtimeFailure = 0 }) {
  positive(interval, 'interval');
  need(defectRate >= 0, 'defectRate must be non-negative');
  need(typeof delayCdf === 'function', 'delayCdf must be a function');
  const b = integrate(delayCdf, 0, interval) / interval;
  const nDef = defectRate * interval;
  const cyc = interval + downtimeInspection;
  return {
    failureFraction: b, costRate: (nDef * (b * costFailure + (1 - b) * costDefectRepair) + costInspection) / cyc,
    downtimeRate: (nDef * b * downtimeFailure + downtimeInspection) / cyc, failuresPerInterval: nDef * b,
  };
}

/**
 * Inspection interval that minimises the delay-time cost rate (or downtime rate).
 * @param {object} o arguments of {@link delayTimeModel} without `interval`, plus:
 * @param {number} o.minInterval [time]
 * @param {number} o.maxInterval [time]
 * @param {string} [o.objective='cost'] 'cost' | 'downtime'
 * @returns {{interval:number, costRate:number, downtimeRate:number, failureFraction:number, atBound:boolean}}
 */
export function optimalInspectionInterval(o) {
  const { minInterval, maxInterval, objective = 'cost' } = o;
  positive(minInterval, 'minInterval');
  need(maxInterval > minInterval, 'maxInterval must exceed minInterval');
  need(objective === 'cost' || objective === 'downtime', 'objective must be "cost" or "downtime"');
  const f = (T) => {
    const r = delayTimeModel({ ...o, interval: T });
    return objective === 'cost' ? r.costRate : r.downtimeRate;
  };
  const r = minimise1d(f, minInterval, maxInterval);
  const res = delayTimeModel({ ...o, interval: r.x });
  const atBound = r.x <= minInterval * 1.001 || r.x >= maxInterval * 0.999;
  return { interval: r.x, costRate: res.costRate, downtimeRate: res.downtimeRate, failureFraction: res.failureFraction, atBound };
}

/**
 * Expected annual downtime and downtime cost of a repairable item or fleet.
 *   downtime per failure = MTTR + logistic delay;  failures per year = λ·(operating hours, ≈ 8760·A)
 * Reference: IEC 60300-3-3 (life-cycle costing); elementary renewal argument.
 * @param {object} o
 * @param {number} o.failureRate failures per operating hour [1/h]
 * @param {number} o.mttr active repair time [h]
 * @param {number} [o.logisticDelay=0] waiting time for access, spares and crew [h]
 * @param {number} [o.units=1] number of identical items [-]
 * @param {number} [o.downtimeCostPerHour=0] lost production and penalties [currency/h]
 * @param {number} [o.repairCost=0] direct cost per repair [currency]
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{availability:number, failuresPerYear:number, downtimeHoursPerYear:number, downtimeCostPerYear:number,
 *   repairCostPerYear:number, totalCostPerYear:number}} totals for all units
 */
export function expectedDowntime({ failureRate, mttr, logisticDelay = 0, units = 1, downtimeCostPerHour = 0, repairCost = 0, hoursPerYear = 8760 }) {
  need(failureRate >= 0, 'failureRate must be non-negative');
  need(mttr >= 0 && logisticDelay >= 0, 'mttr and logisticDelay must be non-negative');
  need(units > 0, 'units must be positive');
  const down = mttr + logisticDelay;
  const A = 1 / (1 + failureRate * down);
  const failures = failureRate * hoursPerYear * A * units;
  const hours = failures * down;
  return {
    availability: A, failuresPerYear: failures, downtimeHoursPerYear: hours, downtimeCostPerYear: hours * downtimeCostPerHour,
    repairCostPerYear: failures * repairCost, totalCostPerYear: hours * downtimeCostPerHour + failures * repairCost,
  };
}

/**
 * Default normalised power curve: cubic between cut-in and rated wind speed, 1 up to cut-out, 0 elsewhere.
 * @param {number} wind hub-height wind speed [m/s]
 * @param {object} [o]
 * @param {number} [o.cutIn=3] [m/s]
 * @param {number} [o.rated=12] [m/s]
 * @param {number} [o.cutOut=25] [m/s]
 * @returns {number} power / rated power [-]
 */
export function genericPowerFraction(wind, { cutIn = 3, rated = 12, cutOut = 25 } = {}) {
  if (!(wind >= cutIn) || wind > cutOut) return 0;
  if (wind >= rated) return 1;
  return (wind ** 3 - cutIn ** 3) / (rated ** 3 - cutIn ** 3);
}

/**
 * Offshore-wind operations-and-maintenance simulator (time-stepping Monte Carlo).
 * Each turbine fails as a Poisson process in operating time with the sum of its component failure rates; a
 * failure stops the turbine (unless the repair class says otherwise), waits for the lead time (spares and, for
 * on-demand vessels, mobilisation), then needs `repairHours` of work that only progresses in time steps where the
 * vessel's limits on significant wave height and wind speed are met and a vessel of that class is free.
 * The weather series is repeated cyclically over the simulated period.
 * Outputs follow the usual definitions: time-based availability = turbine up-time / total time; energy-based
 * availability = energy produced / energy that would have been produced with no downtime.
 * Reference: model structure after Dinwoodie, Endrerud, Hofmann, Martin & Sperstad (2015) "Reference cases for
 * verification of operation and maintenance simulation models for offshore wind farms", Wind Engineering 39(1),
 * 1–14; Carroll, McDonald & McMillan (2016) Wind Energy 19, 1107–1119 for typical failure rates and repair times
 * (supplied by the caller).
 * Limits: no travel time or shift pattern, no opportunistic grouping, repairs are pre-emptible in bad weather,
 * weather persistence comes only from the supplied series.
 * @param {object} o
 * @param {number} o.turbines number of turbines [-]
 * @param {number} o.ratedPower rated power per turbine [W]
 * @param {number} o.years simulated period [yr]
 * @param {{hs: ArrayLike<number>, wind: ArrayLike<number>, dtHours:number}} o.weather significant wave height [m],
 *   hub-height wind speed [m/s], time step [h]
 * @param {{name:string, failureRate:number, repairHours:number, vessel:string, partsCost?:number,
 *   leadTimeHours?:number, technicians?:number, stopsTurbine?:boolean}[]} o.components failureRate per turbine
 *   per operating year [1/yr]
 * @param {Object<string, {hsLimit?:number, windLimit?:number, count?:number, dayRate?:number, charter?:string,
 *   mobilisationHours?:number, mobilisationCost?:number}>} o.vessels limits [m], [m/s]; count (default
 *   unlimited); dayRate [currency/day]; charter 'permanent' (paid every day) or 'on-demand' (paid for working days
 *   plus mobilisationCost per job, with mobilisationHours added to the lead time)
 * @param {function(number):number} [o.powerFraction] normalised power curve; default {@link genericPowerFraction}
 * @param {number} [o.energyPricePerMWh=0] [currency/MWh]
 * @param {number} [o.labourRatePerHour=0] cost per technician-hour [currency/h]
 * @param {number} [o.simulations=10] number of Monte Carlo histories [-]
 * @param {number} [o.seed=1]
 * @returns {{timeAvailability:number, energyAvailability:number, energyProducedMWhPerYear:number,
 *   energyLostMWhPerYear:number, capacityFactor:number, failuresPerYear:number, meanDowntimeHours:number,
 *   costPerYear: {vessels:number, parts:number, labour:number, lostRevenue:number, total:number},
 *   vesselUtilisation: Object<string, number>, byComponent: {name:string, failuresPerYear:number,
 *   meanDowntimeHours:number}[], accessibility: Object<string, number>, availabilityStd:number}}
 *   availabilities and utilisation are fractions [-] (utilisation is in vessel-equivalents when no count is given);
 *   meanDowntimeHours is from failure to completed repair; costs [currency/yr], means over the simulations
 */
export function simulateWindOM({ turbines, ratedPower, years, weather, components, vessels, powerFraction, energyPricePerMWh = 0,
  labourRatePerHour = 0, simulations = 10, seed = 1 }) {
  need(Number.isInteger(turbines) && turbines > 0, 'turbines must be a positive integer');
  positive(ratedPower, 'ratedPower');
  positive(years, 'years');
  need(weather && weather.dtHours > 0 && weather.hs && weather.wind && weather.hs.length === weather.wind.length && weather.hs.length > 0,
    'weather must provide hs, wind (equal length) and dtHours');
  need(Array.isArray(components) && components.length > 0, 'components must be a non-empty array');
  need(vessels && typeof vessels === 'object', 'vessels must be an object');
  need(Number.isInteger(simulations) && simulations > 0, 'simulations must be a positive integer');
  const dt = weather.dtHours;
  const nW = weather.hs.length;
  const nSteps = Math.round((years * 8760) / dt);
  const pf = powerFraction || ((w) => genericPowerFraction(w));
  const vNames = Object.keys(vessels);
  let rateTotal = 0;
  components.forEach((c, i) => {
    need(c.failureRate >= 0, `components[${i}].failureRate must be non-negative`);
    need(c.repairHours > 0, `components[${i}].repairHours must be positive`);
    need(vessels[c.vessel] !== undefined, `components[${i}].vessel "${c.vessel}" is not defined in vessels`);
    rateTotal += c.failureRate;
  });
  const cumRate = [];
  components.reduce((s, c) => { cumRate.push(s + c.failureRate); return s + c.failureRate; }, 0);
  const lamHour = rateTotal / 8760;
  // Weather-dependent quantities are the same for every history.
  const power = new Float64Array(nW);
  const access = {};
  const accessibility = {};
  for (const v of vNames) access[v] = new Uint8Array(nW);
  for (let i = 0; i < nW; i++) {
    power[i] = Math.min(1, Math.max(0, pf(weather.wind[i])));
    for (const v of vNames) {
      const lim = vessels[v];
      access[v][i] = (lim.hsLimit === undefined || weather.hs[i] <= lim.hsLimit) && (lim.windLimit === undefined || weather.wind[i] <= lim.windLimit) ? 1 : 0;
    }
  }
  for (const v of vNames) accessibility[v] = access[v].reduce((s, x) => s + x, 0) / nW;

  const rng = makeRng(seed);
  const sum = { eProd: 0, ePot: 0, fails: 0, repairs: 0, downH: 0, parts: 0, labour: 0, vessel: 0 };
  const busy = Object.fromEntries(vNames.map((v) => [v, 0]));
  const compFails = new Float64Array(components.length);
  const compDown = new Float64Array(components.length);
  const compRepairs = new Float64Array(components.length);
  let avSum = 0;
  let avSum2 = 0;
  const drawLife = () => (lamHour > 0 ? -Math.log(1 - rng.next()) / lamHour : Infinity);
  for (let s = 0; s < simulations; s++) {
    const life = new Float64Array(turbines);
    const down = new Uint8Array(turbines);
    for (let k = 0; k < turbines; k++) life[k] = drawLife();
    const queues = Object.fromEntries(vNames.map((v) => [v, []]));
    let upSteps = 0;
    for (let step = 0; step < nSteps; step++) {
      const w = step % nW;
      const tNow = step * dt;
      const p = power[w];
      let upCount = 0;
      for (let k = 0; k < turbines; k++) {
        if (down[k]) continue;
        upCount++;
        life[k] -= dt;
        if (life[k] <= 0) {
          const u = rng.next() * rateTotal;
          let ci = 0;
          while (ci < components.length - 1 && u >= cumRate[ci]) ci++;
          const c = components[ci];
          const vs = vessels[c.vessel];
          const onDemand = vs.charter === 'on-demand';
          const lead = Math.max(c.leadTimeHours || 0, onDemand ? vs.mobilisationHours || 0 : 0);
          const stops = c.stopsTurbine !== false;
          queues[c.vessel].push({ turbine: k, comp: ci, readyAt: tNow + dt + lead, remaining: c.repairHours, start: tNow + dt, stops });
          if (stops) down[k] = 1; else life[k] = drawLife();
          sum.fails++;
          compFails[ci]++;
          sum.parts += c.partsCost || 0;
          sum.labour += (c.technicians || 0) * c.repairHours * labourRatePerHour;
          if (onDemand) sum.vessel += vs.mobilisationCost || 0;
        }
      }
      upSteps += upCount;
      sum.ePot += turbines * p * dt;
      sum.eProd += upCount * p * dt;
      for (const v of vNames) {
        const q = queues[v];
        if (q.length === 0 || !access[v][w]) continue;
        const cap = vessels[v].count === undefined ? Infinity : vessels[v].count;
        let served = 0;
        for (let j = 0; j < q.length && served < cap; j++) {
          const job = q[j];
          if (job.readyAt > tNow) continue;
          served++;
          job.remaining -= dt;
          busy[v] += dt;
          if (job.remaining <= 1e-9) {
            const dur = tNow + dt - job.start;
            sum.downH += dur;
            sum.repairs++;
            compDown[job.comp] += dur;
            compRepairs[job.comp]++;
            if (job.stops) { down[job.turbine] = 0; life[job.turbine] = drawLife(); }
            q.splice(j, 1);
            j--;
          }
        }
      }
    }
    const av = upSteps / (nSteps * turbines);
    avSum += av;
    avSum2 += av * av;
  }
  const totYears = simulations * ((nSteps * dt) / 8760);
  const vesselUtilisation = {};
  for (const v of vNames) {
    const vs = vessels[v];
    const count = vs.count === undefined ? null : vs.count;
    const days = (nSteps * dt) / 24;
    if (vs.charter === 'on-demand') sum.vessel += (vs.dayRate || 0) * (busy[v] / 24);
    else sum.vessel += (vs.dayRate || 0) * (count === null ? 1 : count) * days * simulations;
    vesselUtilisation[v] = busy[v] / (simulations * nSteps * dt * (count === null ? 1 : count));
  }
  const eProd = (sum.eProd * ratedPower) / 1e6 / totYears;
  const ePot = (sum.ePot * ratedPower) / 1e6 / totYears;
  const lostRevenue = (ePot - eProd) * energyPricePerMWh;
  const cost = { vessels: sum.vessel / totYears, parts: sum.parts / totYears, labour: sum.labour / totYears, lostRevenue, total: 0 };
  cost.total = cost.vessels + cost.parts + cost.labour + cost.lostRevenue;
  const meanAv = avSum / simulations;
  return {
    timeAvailability: meanAv, energyAvailability: sum.ePot > 0 ? sum.eProd / sum.ePot : 1,
    energyProducedMWhPerYear: eProd, energyLostMWhPerYear: ePot - eProd,
    capacityFactor: eProd / ((turbines * ratedPower * 8760) / 1e6), failuresPerYear: sum.fails / totYears,
    meanDowntimeHours: sum.repairs > 0 ? sum.downH / sum.repairs : 0,
    costPerYear: cost, vesselUtilisation,
    byComponent: components.map((c, i) => ({ name: c.name, failuresPerYear: compFails[i] / totYears, meanDowntimeHours: compRepairs[i] > 0 ? compDown[i] / compRepairs[i] : 0 })),
    accessibility, availabilityStd: Math.sqrt(Math.max(0, avSum2 / simulations - meanAv * meanAv)),
  };
}
