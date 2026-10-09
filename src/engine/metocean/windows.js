// Weather-window, operability and persistence analysis on metocean time series; Markov-chain sea-state model;
// Monte Carlo duration of a sequence of weather-restricted marine operations.
//
// Conventions: a series is an object of equally long arrays keyed by variable name (e.g. { hs, tp, wind }) sampled
// at a constant interval dtHours (h). Durations are in hours. Optional time stamps are epoch milliseconds (UTC).

import { makeRng } from '../math/rng.js';
import { percentile } from '../math/stats.js';

function reqPos(name, v) {
  if (typeof v !== 'number' || !(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a positive finite number`);
}

function seriesLength(series, limits) {
  let n = -1;
  for (const lim of limits) {
    const a = series[lim.key];
    if (!a || typeof a.length !== 'number') throw new RangeError(`series has no variable '${lim.key}'`);
    if (n < 0) n = a.length; else if (a.length !== n) throw new RangeError('series variables must have the same length');
  }
  if (n < 1) throw new RangeError('limits must name at least one non-empty series variable');
  return n;
}

/**
 * Operation reference period and classification of a marine operation:
 *   TR = TPOP + TC   (planned operation period plus contingency time).
 * An operation is treated as weather-restricted when TR and TPOP do not exceed the stated thresholds, otherwise it
 * is weather-unrestricted and must be designed for seasonal extremes.
 * Reference: DNV-ST-N001 (2021) Sec. 2, operational criteria (TR ≤ 96 h and TPOP ≤ 72 h for weather-restricted operations; the
 * contingency time covers uncertainty in TPOP and unforeseen events and must be justified case by case, so it is
 * an explicit argument here).
 * @param {object} o
 * @param {number} o.plannedHours planned operation period TPOP, h
 * @param {number} o.contingencyHours contingency time TC, h
 * @param {number} [o.maxReferenceHours=96] upper limit of TR for a weather-restricted operation, h
 * @param {number} [o.maxPlannedHours=72] upper limit of TPOP for a weather-restricted operation, h
 * @returns {{TR:number,TPOP:number,TC:number,weatherRestricted:boolean}} hours
 */
export function operationReferencePeriod({ plannedHours, contingencyHours, maxReferenceHours = 96, maxPlannedHours = 72 }) {
  reqPos('plannedHours', plannedHours);
  if (!(contingencyHours >= 0)) throw new RangeError('contingencyHours must be non-negative');
  const TR = plannedHours + contingencyHours;
  return { TR, TPOP: plannedHours, TC: contingencyHours, weatherRestricted: TR <= maxReferenceHours && plannedHours <= maxPlannedHours };
}

/**
 * Forecast (operational) limit from the design limit, OPWF = α · OPLIM.
 * Reference: DNV-ST-N001 (2021) Sec. 2, alpha-factor tables: α depends on the planned operation period,
 * the design wave height and the forecast/monitoring level (roughly 0.5–1.0 for waves, 0.8 typical for wind).
 * The tabulated values are not reproduced here; pass the value read from the applicable table.
 * @param {number} designLimit design limit OPLIM (e.g. Hs in m, wind in m/s)
 * @param {number} alpha alpha factor, 0 < α ≤ 1
 * @returns {number} forecast limit OPWF, same unit as designLimit
 */
export function forecastLimit(designLimit, alpha) {
  if (!(alpha > 0 && alpha <= 1)) throw new RangeError('alpha must lie in (0, 1]');
  if (!Number.isFinite(designLimit)) throw new RangeError('designLimit must be a finite number');
  return alpha * designLimit;
}

/**
 * Workability mask of a multi-variable series: a step is workable when every limit is satisfied,
 *   min ≤ value ≤ alpha · max   (either bound optional; alpha defaults to 1 and scales the upper bound only).
 * Missing values (NaN) are not workable.
 * Reference: DNV-ST-N001 (2021) Sec. 2 (operational criteria).
 * @param {Object<string, ArrayLike<number>>} series variables sampled at a constant interval
 * @param {{key:string,max?:number,min?:number,alpha?:number}[]} limits limits in the units of each variable
 * @returns {Uint8Array} 1 where all limits are met, 0 otherwise
 */
export function workabilityMask(series, limits) {
  if (!Array.isArray(limits) || limits.length === 0) throw new RangeError('limits must be a non-empty array');
  const n = seriesLength(series, limits);
  const mask = new Uint8Array(n).fill(1);
  for (const lim of limits) {
    const a = series[lim.key];
    const hi = lim.max === undefined || lim.max === null ? Infinity : forecastLimit(lim.max, lim.alpha ?? 1);
    const lo = lim.min === undefined || lim.min === null ? -Infinity : lim.min;
    if (!(hi >= lo)) throw new RangeError(`limit on '${lim.key}' has min above the (alpha-reduced) max`);
    for (let i = 0; i < n; i++) {
      const v = a[i];
      if (!(v >= lo && v <= hi)) mask[i] = 0;
    }
  }
  return mask;
}

/**
 * Uninterrupted weather windows in a workability mask. A run of m consecutive workable steps is taken to last
 * m · dtHours (each sample represents its interval).
 * Reference: persistence (duration) statistics; DNV-ST-N001 (2021) Sec. 2 (weather windows).
 * @param {ArrayLike<number>} mask workability mask (truthy = workable)
 * @param {number} dtHours sampling interval, h
 * @param {number} [minHours=0] only windows at least this long are returned, h
 * @returns {{startIndex:number,endIndex:number,startHours:number,durationHours:number}[]} endIndex is inclusive;
 *   startHours measured from the first sample
 */
export function findWindows(mask, dtHours, minHours = 0) {
  reqPos('dtHours', dtHours);
  const out = [];
  const n = mask.length;
  let i = 0;
  while (i < n) {
    if (mask[i]) {
      let j = i;
      while (j + 1 < n && mask[j + 1]) j++;
      const dur = (j - i + 1) * dtHours;
      if (dur >= minHours - 1e-9) out.push({ startIndex: i, endIndex: j, startHours: i * dtHours, durationHours: dur });
      i = j + 1;
    } else i++;
  }
  return out;
}

// Remaining workable duration from each step (steps), and waiting steps to the next start with >= need steps.
function startability(mask, need, circular) {
  const n = mask.length;
  const remaining = new Int32Array(n);
  const wait = new Int32Array(n).fill(-1);
  const passes = circular ? 2 : 1;
  let run = 0;
  for (let p = 0; p < passes; p++) {
    for (let i = n - 1; i >= 0; i--) {
      run = mask[i] ? Math.min(run + 1, n) : 0;
      remaining[i] = run;
    }
  }
  let w = -1;
  for (let p = 0; p < passes; p++) {
    for (let i = n - 1; i >= 0; i--) {
      if (remaining[i] >= need) w = 0;
      else if (w >= 0) w++;
      wait[i] = w;
    }
  }
  return { remaining, wait };
}

const SEASONS = ['DJF', 'MAM', 'JJA', 'SON'];
const SEASON_OF_MONTH = [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0];

/**
 * Operability analysis of a time series against operational limits and a required uninterrupted duration.
 *   workable fraction   = share of time with all limits satisfied;
 *   startable fraction  = share of time at which an operation needing `requiredHours` of uninterrupted workable
 *                         weather could start (this is the operability relevant for planning);
 *   waiting time        = time from each step to the next startable step (zero if startable), averaged over all
 *                         steps from which a window is reached before the end of the record.
 * Monthly and seasonal breakdowns are given when time stamps are supplied.
 * Reference: DNV-ST-N001 (2021) Sec. 2 (TR, alpha factor, weather windows); DNV-RP-N103 (2017) Sec. 9
 * (weather criteria and availability analysis).
 * Validity: constant sampling interval; gaps must be represented by NaN. The record should span whole years for
 * unbiased annual figures.
 * @param {object} o
 * @param {Object<string, ArrayLike<number>>} o.series variables
 * @param {{key:string,max?:number,min?:number,alpha?:number}[]} o.limits operational limits (alpha reduces max)
 * @param {number} o.dtHours sampling interval, h
 * @param {number} [o.requiredHours=dtHours] required uninterrupted duration (operation reference period TR), h
 * @param {ArrayLike<number>} [o.timeMs] epoch milliseconds (UTC) of each sample, for monthly statistics
 * @returns {{n:number,requiredSteps:number,workableFraction:number,startableFraction:number,
 *   meanWaitingHours:number,maxWaitingHours:number,p90WaitingHours:number,
 *   windows:{startIndex:number,endIndex:number,startHours:number,durationHours:number}[],
 *   monthly:{month:number,n:number,workableFraction:number,startableFraction:number,meanWaitingHours:number}[],
 *   seasonal:{season:string,n:number,workableFraction:number,startableFraction:number,meanWaitingHours:number}[],
 *   mask:Uint8Array}}
 *   fractions 0–1; waiting times h (NaN when no window exists); windows are those lasting at least requiredHours;
 *   month is 1–12; monthly and seasonal are empty arrays without timeMs.
 */
export function operability({ series, limits, dtHours, requiredHours, timeMs }) {
  reqPos('dtHours', dtHours);
  const req = requiredHours ?? dtHours;
  reqPos('requiredHours', req);
  const mask = workabilityMask(series, limits);
  const n = mask.length;
  if (timeMs && timeMs.length !== n) throw new RangeError('timeMs must have the same length as the series');
  const need = Math.max(1, Math.ceil(req / dtHours - 1e-9));
  const { remaining, wait } = startability(mask, need, false);
  let work = 0; let start = 0; let wSum = 0; let wN = 0; let wMax = 0;
  const waits = [];
  const mN = new Float64Array(12); const mW = new Float64Array(12); const mS = new Float64Array(12);
  const mWait = new Float64Array(12); const mWaitN = new Float64Array(12);
  for (let i = 0; i < n; i++) {
    const st = remaining[i] >= need ? 1 : 0;
    work += mask[i]; start += st;
    if (wait[i] >= 0) { wSum += wait[i]; wN++; waits.push(wait[i] * dtHours); if (wait[i] > wMax) wMax = wait[i]; }
    if (timeMs) {
      const mo = new Date(timeMs[i]).getUTCMonth();
      if (Number.isNaN(mo)) throw new RangeError('timeMs contains an invalid time stamp');
      mN[mo]++; mW[mo] += mask[i]; mS[mo] += st;
      if (wait[i] >= 0) { mWait[mo] += wait[i]; mWaitN[mo]++; }
    }
  }
  const monthly = []; const seasonal = [];
  if (timeMs) {
    const sN = [0, 0, 0, 0]; const sW = [0, 0, 0, 0]; const sS = [0, 0, 0, 0]; const sWait = [0, 0, 0, 0]; const sWaitN = [0, 0, 0, 0];
    for (let m = 0; m < 12; m++) {
      monthly.push({ month: m + 1, n: mN[m], workableFraction: mN[m] ? mW[m] / mN[m] : NaN,
        startableFraction: mN[m] ? mS[m] / mN[m] : NaN, meanWaitingHours: mWaitN[m] ? (mWait[m] / mWaitN[m]) * dtHours : NaN });
      const s = SEASON_OF_MONTH[m];
      sN[s] += mN[m]; sW[s] += mW[m]; sS[s] += mS[m]; sWait[s] += mWait[m]; sWaitN[s] += mWaitN[m];
    }
    for (let s = 0; s < 4; s++) {
      seasonal.push({ season: SEASONS[s], n: sN[s], workableFraction: sN[s] ? sW[s] / sN[s] : NaN,
        startableFraction: sN[s] ? sS[s] / sN[s] : NaN, meanWaitingHours: sWaitN[s] ? (sWait[s] / sWaitN[s]) * dtHours : NaN });
    }
  }
  return {
    n, requiredSteps: need, workableFraction: work / n, startableFraction: start / n,
    meanWaitingHours: wN ? (wSum / wN) * dtHours : NaN, maxWaitingHours: wN ? wMax * dtHours : NaN,
    p90WaitingHours: wN ? percentile(waits, 90) : NaN,
    windows: findWindows(mask, dtHours, need * dtHours), monthly, seasonal, mask,
  };
}

/**
 * Persistence (duration) statistics of calms or storms: for each threshold, the events are the uninterrupted runs
 * with the variable below (calms) or above (storms) the threshold. For each duration D the table gives
 *   eventExceedance[D]  = fraction of events lasting at least D          (non-exceedance = 1 − this value);
 *   timeFraction[D]     = fraction of the whole record occupied by events lasting at least D.
 * Reference: Graham (1982) "The parameterisation and prediction of wave height and
 * wind speed persistence statistics for oil industry operational planning purposes", Coastal Eng. 6.
 * Validity: constant sampling interval; events cut by the ends of the record are counted with their observed length.
 * @param {object} o
 * @param {ArrayLike<number>} o.values variable (e.g. Hs in m)
 * @param {number} o.dtHours sampling interval, h
 * @param {number[]} o.thresholds thresholds in the unit of the variable
 * @param {number[]} o.durations durations, h
 * @param {'below'|'above'} [o.sense='below'] calms (value ≤ threshold) or storms (value > threshold)
 * @returns {{threshold:number,occurrence:number,nEvents:number,meanDurationHours:number,maxDurationHours:number,
 *   eventExceedance:number[],timeFraction:number[]}[]} one row per threshold; occurrence = fraction of time in
 *   the state; arrays follow `durations`
 */
export function persistenceTable({ values, dtHours, thresholds, durations, sense = 'below' }) {
  reqPos('dtHours', dtHours);
  if (sense !== 'below' && sense !== 'above') throw new RangeError("sense must be 'below' or 'above'");
  const n = values.length;
  if (n < 1) throw new RangeError('values must not be empty');
  const rows = [];
  const mask = new Uint8Array(n);
  for (const th of thresholds) {
    let occ = 0;
    for (let i = 0; i < n; i++) {
      const v = values[i];
      mask[i] = sense === 'below' ? (v <= th ? 1 : 0) : (v > th ? 1 : 0);
      occ += mask[i];
    }
    const ev = findWindows(mask, dtHours, 0);
    let sum = 0; let mx = 0;
    for (const e of ev) { sum += e.durationHours; if (e.durationHours > mx) mx = e.durationHours; }
    const ex = []; const tf = [];
    for (const D of durations) {
      let c = 0; let tt = 0;
      for (const e of ev) if (e.durationHours >= D - 1e-9) { c++; tt += e.durationHours; }
      ex.push(ev.length ? c / ev.length : 0);
      tf.push(tt / (n * dtHours));
    }
    rows.push({ threshold: th, occurrence: occ / n, nEvents: ev.length, meanDurationHours: ev.length ? sum / ev.length : 0,
      maxDurationHours: mx, eventExceedance: ex, timeFraction: tf });
  }
  return rows;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Markov-chain sea-state model                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

function stateOf(v, edges) {
  // edges ascending, length nStates + 1; values outside are assigned to the end bins
  const ns = edges.length - 1;
  if (!(v >= edges[0])) return v < edges[0] ? 0 : -1;
  if (v >= edges[ns]) return ns - 1;
  let lo = 0; let hi = ns;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (edges[m] <= v) lo = m; else hi = m; }
  return lo;
}

/**
 * Stationary distribution π of a row-stochastic transition matrix (π P = π, Σ π = 1) by power iteration on the
 * lazy chain ½(I + P), which converges for periodic chains as well.
 * Reference: Norris (1997) Markov Chains §1.7–1.8.
 * Validity: irreducible chain (a unique stationary distribution); otherwise the result depends on the uniform
 * starting vector.
 * @param {number[][]} P transition matrix, P[i][j] = probability of moving from state i to state j
 * @param {object} [o]
 * @param {number} [o.tol=1e-13] convergence tolerance on the largest change
 * @param {number} [o.maxIter=100000]
 * @returns {number[]} stationary probabilities
 */
export function markovStationary(P, { tol = 1e-13, maxIter = 100000 } = {}) {
  const n = P.length;
  if (n < 1) throw new RangeError('P must not be empty');
  let pi = new Array(n).fill(1 / n);
  let nx = new Array(n).fill(0);
  for (let it = 0; it < maxIter; it++) {
    nx.fill(0);
    for (let i = 0; i < n; i++) {
      const pii = pi[i];
      if (pii === 0) continue;
      const row = P[i];
      for (let j = 0; j < n; j++) nx[j] += pii * row[j];
    }
    let diff = 0; let s = 0;
    for (let j = 0; j < n; j++) { nx[j] = 0.5 * (nx[j] + pi[j]); s += nx[j]; }
    for (let j = 0; j < n; j++) { nx[j] /= s; diff = Math.max(diff, Math.abs(nx[j] - pi[j])); }
    const t = pi; pi = nx; nx = t;
    if (diff < tol) break;
  }
  return pi;
}

/**
 * Fit a first-order Markov chain to a binned time series: maximum-likelihood transition probabilities
 * P[i][j] = N_ij / Σ_j N_ij from the counts of consecutive-step transitions. Pairs containing NaN are skipped.
 * States never left in the record are made absorbing-free by giving them the empirical occupancy distribution.
 * Reference: Anastasiou & Tsekos (1996) "Persistence statistics of marine environmental parameters from Markov
 * theory", Applied Ocean Research 18; Anderson & Goodman (1957) for the estimator.
 * Validity: the sampling interval of the fitted chain is that of the data; first-order chains under-predict long
 * persistence unless the interval is comparable to the weather decorrelation time (3–6 h is customary).
 * @param {ArrayLike<number>} values series (e.g. Hs in m)
 * @param {number[]} edges ascending bin edges, length nStates + 1, in the unit of the variable; values beyond the
 *   outer edges fall in the end bins
 * @returns {{edges:number[],nStates:number,counts:number[][],P:number[][],occupancy:number[],stationary:number[]}}
 *   counts of transitions, transition matrix, empirical state occupancy and stationary distribution
 */
export function markovFit(values, edges) {
  if (!Array.isArray(edges) || edges.length < 2) throw new RangeError('edges must have at least two entries');
  for (let i = 1; i < edges.length; i++) if (!(edges[i] > edges[i - 1])) throw new RangeError('edges must be strictly ascending');
  const ns = edges.length - 1;
  const counts = Array.from({ length: ns }, () => new Array(ns).fill(0));
  const occ = new Array(ns).fill(0);
  let prev = -1; let tot = 0;
  for (let i = 0; i < values.length; i++) {
    const s = stateOf(values[i], edges);
    if (s >= 0) { occ[s]++; tot++; if (prev >= 0) counts[prev][s]++; }
    prev = s;
  }
  if (tot < 2) throw new RangeError('values must contain at least two valid samples');
  const occupancy = occ.map((c) => c / tot);
  const P = counts.map((row) => {
    let s = 0;
    for (let j = 0; j < ns; j++) s += row[j];
    return s > 0 ? row.map((c) => c / s) : occupancy.slice();
  });
  return { edges: edges.slice(), nStates: ns, counts, P, occupancy, stationary: markovStationary(P) };
}

/**
 * Simulate a synthetic series from a fitted Markov chain (deterministic for a seed).
 * Values are placed at the bin mid-point or drawn uniformly within the bin.
 * Reference: as `markovFit`.
 * @param {{edges:number[],P:number[][],stationary:number[]}} model result of `markovFit`
 * @param {object} o
 * @param {number} o.nSteps number of steps (e.g. years × 8766 / dtHours)
 * @param {number} [o.seed=1] seed
 * @param {number} [o.startState] initial state index (default: drawn from the stationary distribution)
 * @param {'mid'|'uniform'} [o.sample='uniform'] value within the bin
 * @returns {{states:Int32Array,values:Float64Array}} state index and value (unit of the fitted variable) per step
 */
export function markovSimulate(model, { nSteps, seed = 1, startState, sample = 'uniform' }) {
  if (!Number.isInteger(nSteps) || nSteps < 1) throw new RangeError('nSteps must be a positive integer');
  if (sample !== 'mid' && sample !== 'uniform') throw new RangeError("sample must be 'mid' or 'uniform'");
  const { edges, P } = model;
  const ns = P.length;
  const rng = makeRng(seed);
  // cumulative rows for inverse-transform sampling
  const cum = P.map((row) => { let s = 0; return row.map((p) => (s += p)); });
  const draw = (c) => {
    const u = rng.next() * c[ns - 1];
    let lo = 0; let hi = ns - 1;
    while (lo < hi) { const m = (lo + hi) >> 1; if (c[m] > u) hi = m; else lo = m + 1; }
    return lo;
  };
  let s;
  if (startState === undefined || startState === null) {
    let a = 0;
    s = draw(model.stationary.map((p) => (a += p)));
  } else {
    if (!Number.isInteger(startState) || startState < 0 || startState >= ns) throw new RangeError('startState is out of range');
    s = startState;
  }
  const states = new Int32Array(nSteps); const values = new Float64Array(nSteps);
  for (let i = 0; i < nSteps; i++) {
    if (i > 0) s = draw(cum[s]);
    states[i] = s;
    const f = sample === 'mid' ? 0.5 : rng.next();
    values[i] = edges[s] + f * (edges[s + 1] - edges[s]);
  }
  return { states, values };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Monte Carlo campaign duration                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Monte Carlo duration of a campaign of sequential weather-restricted operations, by replaying the operations
 * against the (historical or synthetic) series from randomly drawn start times.
 * Each operation waits for the first step from which its limits (with alpha factor) hold without interruption for
 * its reference period TR = netHours + contingencyHours, then takes netHours to execute; the next operation may
 * start immediately afterwards. The series is treated as periodic (it should span whole years).
 *   weather downtime = completion time − Σ net durations.
 * Reference: DNV-ST-N001 (2021) Sec. 2 (TR, alpha factor); DNV-RP-N103 (2017) Sec. 9 (time-domain replay of
 * operations against a weather record).
 * Validity: operations cannot be suspended once started; technical downtime, mobilisation and learning effects
 * are not modelled; the spread reflects only the variability contained in the series supplied.
 * @param {object} o
 * @param {Object<string, ArrayLike<number>>} o.series variables at constant interval
 * @param {number} o.dtHours sampling interval, h
 * @param {{name?:string,limits:{key:string,max?:number,min?:number,alpha?:number}[],netHours:number,
 *   contingencyHours?:number}[]} o.operations operations in execution order
 * @param {number} [o.nTrials=1000] number of random start times
 * @param {number} [o.seed=1] seed
 * @param {ArrayLike<number>} [o.timeMs] epoch ms of the samples; needed for startMonths
 * @param {number[]} [o.startMonths] restrict campaign starts to these months (1–12)
 * @returns {{nTrials:number,netHours:number,completionHours:{p10:number,p50:number,p90:number,mean:number,min:number,max:number},
 *   downtimeHours:{p10:number,p50:number,p90:number,mean:number},
 *   operations:{name:string,requiredHours:number,startableFraction:number,meanWaitHours:number,p90WaitHours:number}[],
 *   samples:Float64Array}} hours; P10/P50/P90 are non-exceedance percentiles of the completion time (P90 is the
 *   conservative value); samples holds the completion time of every trial, h
 */
export function campaignMonteCarlo({ series, dtHours, operations, nTrials = 1000, seed = 1, timeMs, startMonths }) {
  reqPos('dtHours', dtHours);
  if (!Array.isArray(operations) || operations.length === 0) throw new RangeError('operations must be a non-empty array');
  if (!Number.isInteger(nTrials) || nTrials < 1) throw new RangeError('nTrials must be a positive integer');
  const prepared = operations.map((op, idx) => {
    reqPos(`operations[${idx}].netHours`, op.netHours);
    const cont = op.contingencyHours ?? 0;
    if (!(cont >= 0)) throw new RangeError(`operations[${idx}].contingencyHours must be non-negative`);
    const mask = workabilityMask(series, op.limits);
    const need = Math.max(1, Math.ceil((op.netHours + cont) / dtHours - 1e-9));
    const { remaining, wait } = startability(mask, need, true);
    if (wait[0] < 0) throw new Error(`operation '${op.name ?? idx + 1}' never finds a weather window of ${op.netHours + cont} h in the series`);
    let st = 0;
    for (let i = 0; i < mask.length; i++) if (remaining[i] >= need) st++;
    return { name: op.name ?? `operation ${idx + 1}`, wait, netSteps: Math.max(1, Math.ceil(op.netHours / dtHours - 1e-9)),
      netHours: op.netHours, requiredHours: op.netHours + cont, startableFraction: st / mask.length, waits: new Float64Array(nTrials) };
  });
  const n = prepared[0].wait.length;
  let candidates = null;
  if (startMonths && startMonths.length) {
    if (!timeMs || timeMs.length !== n) throw new RangeError('timeMs (same length as the series) is required with startMonths');
    const ok = new Set(startMonths);
    candidates = [];
    for (let i = 0; i < n; i++) if (ok.has(new Date(timeMs[i]).getUTCMonth() + 1)) candidates.push(i);
    if (candidates.length === 0) throw new RangeError('no samples fall in startMonths');
  }
  const rng = makeRng(seed);
  const samples = new Float64Array(nTrials);
  const down = new Float64Array(nTrials);
  let net = 0;
  for (const p of prepared) net += p.netHours;
  for (let tr = 0; tr < nTrials; tr++) {
    const s0 = candidates ? candidates[rng.int(candidates.length)] : rng.int(n);
    let steps = 0; // elapsed steps since the campaign start
    let hours = 0;
    for (const p of prepared) {
      const w = p.wait[(s0 + steps) % n];
      p.waits[tr] = w * dtHours;
      steps += w + p.netSteps;
      hours += w * dtHours + p.netHours;
    }
    samples[tr] = hours;
    down[tr] = hours - net;
  }
  const arr = Array.from(samples); const darr = Array.from(down);
  let mean = 0; let dmean = 0; let mn = Infinity; let mx = -Infinity;
  for (let i = 0; i < nTrials; i++) { mean += arr[i]; dmean += darr[i]; if (arr[i] < mn) mn = arr[i]; if (arr[i] > mx) mx = arr[i]; }
  return {
    nTrials, netHours: net,
    completionHours: { p10: percentile(arr, 10), p50: percentile(arr, 50), p90: percentile(arr, 90), mean: mean / nTrials, min: mn, max: mx },
    downtimeHours: { p10: percentile(darr, 10), p50: percentile(darr, 50), p90: percentile(darr, 90), mean: dmean / nTrials },
    operations: prepared.map((p) => {
      const w = Array.from(p.waits);
      let m = 0;
      for (let i = 0; i < w.length; i++) m += w[i];
      return { name: p.name, requiredHours: p.requiredHours, startableFraction: p.startableFraction, meanWaitHours: m / w.length, p90WaitHours: percentile(w, 90) };
    }),
    samples,
  };
}
