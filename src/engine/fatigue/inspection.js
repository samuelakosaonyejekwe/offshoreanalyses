// Inspection reliability: probability of detection, Bayesian updating of fatigue crack growth after inspections,
// risk-based inspection planning and value of information.
//
// SI throughout (crack sizes in m). The probabilistic crack-growth model is the one-dimensional Paris law with a
// constant geometry factor, a(N) in closed form, which keeps the Monte Carlo fast enough for interactive use.

import { makeRng } from '../math/rng.js';
import { parisCrackSizeClosedForm } from './fracture.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

/**
 * Probability of detection of a crack of size a: POD(a) = 1 − 1/(1 + (a/X0)^b).
 * Reference: DNV-RP-C210 (probabilistic methods for planning of inspection for fatigue cracks), POD curve form.
 * @param {number} a crack size (depth for EC/MPI/ACFM/UT, length for visual) [m]
 * @param {object} o
 * @param {number} o.x0 distribution parameter X0 (median detectable size) [m]
 * @param {number} o.b distribution parameter b [-]
 * @returns {number} probability of detection [-]
 */
export function pod(a, { x0, b }) {
  need(a >= 0, 'a must be non-negative');
  need(x0 > 0, 'x0 must be positive');
  need(b > 0, 'b must be positive');
  if (a === 0) return 0;
  if (!Number.isFinite(a)) return 1;
  return 1 - 1 / (1 + (a / x0) ** b);
}

/**
 * POD parameters quoted from DNV-RP-C210 for common methods. These presets are convenience values: check them
 * against the edition governing the project before use, and prefer passing x0 and b explicitly.
 *  EC / MPI / ACFM (crack depth): good conditions above water X0 = 0.40 mm, b = 1.43; normal conditions above
 *  water 0.45 mm, 0.90; below water 1.16 mm, 0.90. Ultrasonic testing (depth) 0.41 mm, 0.642.
 *  Close visual (crack length): easy access 15.78 mm, 1.079; moderate 37.15 mm, 1.095; difficult 83.03 mm, 1.079.
 * @param {string} method 'ec-good' | 'ec-normal' | 'ec-underwater' | 'ut' | 'visual-easy' | 'visual-moderate' |
 *   'visual-difficult'
 * @returns {{x0:number, b:number, measure:string}} x0 [m], b [-], measure = 'depth' or 'length'
 */
export function podPreset(method) {
  const t = {
    'ec-good': [0.40e-3, 1.43, 'depth'],
    'ec-normal': [0.45e-3, 0.90, 'depth'],
    'ec-underwater': [1.16e-3, 0.90, 'depth'],
    ut: [0.41e-3, 0.642, 'depth'],
    'visual-easy': [15.78e-3, 1.079, 'length'],
    'visual-moderate': [37.15e-3, 1.095, 'length'],
    'visual-difficult': [83.03e-3, 1.079, 'length'],
  }[method];
  need(t !== undefined, `method: unknown POD preset "${method}"`);
  return { x0: t[0], b: t[1], measure: t[2] };
}

function sampler(spec, name) {
  if (typeof spec === 'number') return () => spec;
  need(spec && typeof spec === 'object', `${name} must be a number or a distribution`);
  const { dist = 'deterministic' } = spec;
  if (dist === 'deterministic') return () => spec.value;
  const mean = spec.mean;
  const std = spec.std !== undefined ? spec.std : spec.cov !== undefined ? spec.cov * Math.abs(mean) : undefined;
  if (dist === 'exponential') {
    need(mean > 0, `${name}.mean must be positive`);
    return (rng) => -mean * Math.log(1 - rng.next());
  }
  need(std !== undefined && std >= 0, `${name} needs std or cov`);
  if (dist === 'normal') return (rng) => mean + std * rng.normal();
  if (dist === 'lognormal') {
    need(mean > 0, `${name}.mean must be positive`);
    const s2 = Math.log(1 + (std / mean) ** 2);
    const mu = Math.log(mean) - 0.5 * s2;
    const s = Math.sqrt(s2);
    return (rng) => Math.exp(mu + s * rng.normal());
  }
  throw new RangeError(`${name}.dist: unsupported distribution "${dist}"`);
}

/**
 * Monte Carlo population of fatigue cracks for inspection updating. Each sample carries an initial size a0, a
 * Paris coefficient C, a stress-modelling factor and a critical size; crack size at any time follows the closed
 * form for a constant geometry factor, with the limit state g = a_crit − a(t).
 * Random variables are given as a number or { dist: 'lognormal' | 'normal' | 'exponential' | 'deterministic',
 * mean, std | cov, value }.
 * Reference: Madsen, Skjong & Kirkemo (1987) "Probabilistic fatigue analysis of offshore structures —
 * reliability updating through inspection results"; DNV-RP-C210 (fracture-mechanics based planning).
 * Typical inputs (DNV-RP-C210): a0 exponential with mean 0.11 mm; ln C normal with the BS 7910 mean and
 * mean + 2 SD constants; stress factor lognormal with CoV 0.1–0.3.
 * @param {object} o
 * @param {number|object} o.a0 initial crack depth [m]
 * @param {number|object} o.C Paris coefficient in SI units [m/cycle/(Pa·√m)^m]
 * @param {number} o.m Paris exponent [-]
 * @param {number} o.stressRange equivalent constant-amplitude stress range (E[S^m])^(1/m) [Pa]
 * @param {number|object} [o.stressFactor=1] multiplicative uncertainty on stress [-]
 * @param {number} o.Y geometry factor [-]
 * @param {number} o.cyclesPerYear [cycles/yr]
 * @param {number|object} o.aCrit critical crack depth (e.g. wall thickness) [m]
 * @param {number} [o.n=20000] number of samples [-]
 * @param {number} [o.seed=1]
 * @returns {{n:number, m:number, a0: Float64Array, k: Float64Array, aCrit: Float64Array, weights: Float64Array,
 *   timeToFailure: Float64Array}} k = per-year growth constant C·(Y·S·√π)^m·ν (SI); weights start at 1;
 *   timeToFailure in years
 */
export function crackPopulation({ a0, C, m, stressRange, stressFactor = 1, Y, cyclesPerYear, aCrit, n = 20000, seed = 1 }) {
  need(m > 0, 'm must be positive');
  need(stressRange > 0, 'stressRange must be positive');
  need(Y > 0, 'Y must be positive');
  need(cyclesPerYear > 0, 'cyclesPerYear must be positive');
  need(Number.isInteger(n) && n > 0, 'n must be a positive integer');
  const rng = makeRng(seed);
  const sA = sampler(a0, 'a0');
  const sC = sampler(C, 'C');
  const sS = sampler(stressFactor, 'stressFactor');
  const sCr = sampler(aCrit, 'aCrit');
  const pop = { n, m, a0: new Float64Array(n), k: new Float64Array(n), aCrit: new Float64Array(n), weights: new Float64Array(n).fill(1), timeToFailure: new Float64Array(n) };
  const e = 1 - m / 2;
  for (let i = 0; i < n; i++) {
    const a = Math.max(sA(rng), 1e-9);
    const c = Math.max(sC(rng), 0);
    const s = Math.max(sS(rng), 0) * stressRange;
    const ac = Math.max(sCr(rng), 1e-9);
    const k = c * (Y * s * Math.sqrt(Math.PI)) ** m * cyclesPerYear;
    pop.a0[i] = a;
    pop.k[i] = k;
    pop.aCrit[i] = ac;
    if (a >= ac) pop.timeToFailure[i] = 0;
    else if (k === 0) pop.timeToFailure[i] = Infinity;
    else pop.timeToFailure[i] = Math.abs(e) < 1e-12 ? Math.log(ac / a) / k : (ac ** e - a ** e) / (e * k);
  }
  return pop;
}

function sizeAt(pop, i, t) {
  if (t >= pop.timeToFailure[i]) return Infinity;
  return parisCrackSizeClosedForm({ a0: pop.a0[i], cycles: t, stressRange: 1, Y: 1 / Math.sqrt(Math.PI), A: pop.k[i], m: pop.m });
}

/**
 * Crack sizes of the population at time t (Infinity for members that have already failed).
 * @param {object} pop population from {@link crackPopulation}
 * @param {number} t time [yr]
 * @returns {Float64Array} crack size [m]
 */
export function crackSizesAt(pop, t) {
  need(t >= 0, 't must be non-negative');
  const out = new Float64Array(pop.n);
  for (let i = 0; i < pop.n; i++) out[i] = sizeAt(pop, i, t);
  return out;
}

/**
 * Bayesian update of the crack population for an inspection outcome, by likelihood weighting:
 *  no detection:               w ← w·[1 − POD(a(t))], and members that would already have failed get w = 0;
 *  detection with measurement: w ← w·POD(a(t))·φ((a_m − a(t))/σ_m);
 *  detection without sizing:   w ← w·POD(a(t)).
 * The component is known to have survived to the inspection, so failed members are excluded in every case.
 * Reference: Madsen et al. (1987); DNV-RP-C210 (updating by inspection events).
 * @param {object} pop population (not modified)
 * @param {object} o
 * @param {number} o.time inspection time [yr]
 * @param {{x0:number, b:number}} o.pod POD parameters
 * @param {boolean} o.detected inspection outcome
 * @param {number} [o.measuredSize] measured crack size [m]
 * @param {number} [o.sizingStd] standard deviation of the sizing error [m] (required with measuredSize)
 * @returns {object} new population sharing the samples with updated weights
 */
export function updateAfterInspection(pop, { time, pod: podPar, detected, measuredSize, sizingStd }) {
  need(time >= 0, 'time must be non-negative');
  need(typeof detected === 'boolean', 'detected must be a boolean');
  const w = new Float64Array(pop.n);
  if (measuredSize !== undefined) need(sizingStd > 0, 'sizingStd must be positive when measuredSize is given');
  for (let i = 0; i < pop.n; i++) {
    if (pop.weights[i] === 0 || time >= pop.timeToFailure[i]) continue;
    const a = sizeAt(pop, i, time);
    const p = pod(a, podPar);
    let like;
    if (!detected) like = 1 - p;
    else if (measuredSize !== undefined) {
      const z = (measuredSize - a) / sizingStd;
      like = p * Math.exp(-0.5 * z * z);
    } else like = p;
    w[i] = pop.weights[i] * like;
  }
  return { ...pop, weights: w };
}

/**
 * Probability that an inspection at a given time finds a crack, given survival to that time and the information
 * already carried by the weights.
 * @param {object} pop population
 * @param {number} time [yr]
 * @param {{x0:number, b:number}} podPar
 * @returns {number} probability of detection at the inspection given survival [-]
 */
export function detectionProbability(pop, time, podPar) {
  let num = 0;
  let den = 0;
  for (let i = 0; i < pop.n; i++) {
    if (pop.weights[i] === 0 || time >= pop.timeToFailure[i]) continue;
    den += pop.weights[i];
    num += pop.weights[i] * pod(sizeAt(pop, i, time), podPar);
  }
  return den > 0 ? num / den : 0;
}

/**
 * Failure probability of the (weighted) population over time: cumulative P_f(t) = P[T_f <= t] and annual
 * probability ΔP_f in the year ending at t, conditional on the information carried by the weights.
 * @param {object} pop population
 * @param {ArrayLike<number>} times [yr]
 * @param {object} [o]
 * @param {number} [o.survivedTo=0] condition on survival to this time [yr]
 * @returns {{pf: Float64Array, annualPf: Float64Array, effectiveSamples:number}}
 */
export function populationFailureProbability(pop, times, { survivedTo = 0 } = {}) {
  let sw = 0;
  let sw2 = 0;
  for (let i = 0; i < pop.n; i++) {
    if (pop.timeToFailure[i] <= survivedTo) continue;
    sw += pop.weights[i];
    sw2 += pop.weights[i] ** 2;
  }
  const P = (t) => {
    if (!(sw > 0) || t <= survivedTo) return 0;
    let s = 0;
    for (let i = 0; i < pop.n; i++) {
      const tf = pop.timeToFailure[i];
      if (tf > survivedTo && tf <= t) s += pop.weights[i];
    }
    return s / sw;
  };
  const pf = new Float64Array(times.length);
  const annualPf = new Float64Array(times.length);
  for (let j = 0; j < times.length; j++) {
    pf[j] = P(times[j]);
    annualPf[j] = pf[j] - P(times[j] - 1);
  }
  return { pf, annualPf, effectiveSamples: sw2 > 0 ? (sw * sw) / sw2 : 0 };
}

/**
 * Risk-based inspection plan: inspections are placed just before the annual failure probability would exceed
 * the target, assuming each inspection finds nothing (the usual planning assumption), and the population is
 * updated after each.
 * Reference: DNV-RP-C210 (inspection planning for a target annual failure probability); Straub &
 * Faber (2005) Structural Safety 27, 335–355.
 * Monte Carlo resolution limits the smallest target that can be resolved to roughly 10/n.
 * @param {object} o
 * @param {object} o.population prior population from {@link crackPopulation}
 * @param {{x0:number, b:number}} o.pod POD parameters of the planned method
 * @param {number} o.targetAnnualPf acceptable annual failure probability [-]
 * @param {number} o.serviceLife planning horizon [yr]
 * @param {number} [o.maxInspections=50]
 * @returns {{inspectionTimes:number[], table: {year:number, annualPf:number, cumulativePf:number,
 *   inspection:boolean}[], maxAnnualPf:number, meetsTarget:boolean}} table has one row per year; probabilities are
 *   conditional on survival and on no detection at the inspections so far
 */
export function planInspections({ population, pod: podPar, targetAnnualPf, serviceLife, maxInspections = 50 }) {
  need(targetAnnualPf > 0 && targetAnnualPf < 1, 'targetAnnualPf must be in (0, 1)');
  need(serviceLife >= 1, 'serviceLife must be at least 1 year');
  let pop = population;
  let last = 0;
  const inspectionTimes = [];
  const table = [];
  let maxAnnual = 0;
  let meets = true;
  const years = Math.floor(serviceLife);
  for (let y = 1; y <= years; y++) {
    let r = populationFailureProbability(pop, [y], { survivedTo: last });
    let annual = annualGiven(pop, y, last);
    let inspected = false;
    if (annual > targetAnnualPf && y - 1 > last && inspectionTimes.length < maxInspections) {
      pop = updateAfterInspection(pop, { time: y - 1, pod: podPar, detected: false });
      last = y - 1;
      inspectionTimes.push(y - 1);
      inspected = true;
      r = populationFailureProbability(pop, [y], { survivedTo: last });
      annual = annualGiven(pop, y, last);
    }
    if (annual > targetAnnualPf) meets = false;
    if (annual > maxAnnual) maxAnnual = annual;
    table.push({ year: y, annualPf: annual, cumulativePf: r.pf[0], inspection: inspected });
  }
  return { inspectionTimes, table, maxAnnualPf: maxAnnual, meetsTarget: meets };
}

// Probability of failure in (y − 1, y] given survival to y − 1 and the inspection information up to `last`.
function annualGiven(pop, y, last) {
  let alive = 0;
  let fail = 0;
  const start = Math.max(y - 1, last);
  for (let i = 0; i < pop.n; i++) {
    const tf = pop.timeToFailure[i];
    if (tf <= start) continue;
    alive += pop.weights[i];
    if (tf <= y) fail += pop.weights[i];
  }
  return alive > 0 ? fail / alive : 0;
}

/**
 * Value of information of an inspection in a two-state decision problem (damaged / intact) with a repair action.
 * Without inspection the cheaper of "do nothing" (expected failure cost) and "repair anyway" is chosen; with
 * inspection, the cheaper action is chosen for each outcome (indication / no indication). VoI = expected cost without − expected cost with (inspection cost
 * excluded); the inspection is worthwhile when VoI exceeds its cost.
 * Reference: Raiffa & Schlaifer (1961) pre-posterior analysis; Straub (2014) Structural Safety 49, 75–85.
 * @param {object} o
 * @param {number} o.pDamage prior probability that the detail is damaged [-]
 * @param {number} o.pFailGivenDamage probability of failure in the period if damaged and not repaired [-]
 * @param {number} [o.pFailGivenIntact=0] [-]
 * @param {number} o.pod probability of an indication given damage [-]
 * @param {number} [o.falseAlarm=0] probability of an indication given no damage [-]
 * @param {number} o.costFailure [currency]
 * @param {number} o.costRepair [currency]
 * @param {number} [o.costInspection=0] [currency]
 * @param {number} [o.pFailAfterRepair=0] residual failure probability after repair [-]
 * @returns {{costWithout:number, costWith:number, valueOfInformation:number, netValue:number, worthwhile:boolean,
 *   actionWithout:string, posteriorDamageGivenIndication:number, posteriorDamageGivenClear:number}}
 */
export function valueOfInformation({ pDamage, pFailGivenDamage, pFailGivenIntact = 0, pod: pd, falseAlarm = 0, costFailure, costRepair, costInspection = 0, pFailAfterRepair = 0 }) {
  for (const [v, name] of [[pDamage, 'pDamage'], [pFailGivenDamage, 'pFailGivenDamage'], [pFailGivenIntact, 'pFailGivenIntact'], [pd, 'pod'], [falseAlarm, 'falseAlarm'], [pFailAfterRepair, 'pFailAfterRepair']]) {
    need(v >= 0 && v <= 1, `${name} must be in [0, 1]`);
  }
  need(costFailure >= 0 && costRepair >= 0 && costInspection >= 0, 'costs must be non-negative');
  const pfPrior = pDamage * pFailGivenDamage + (1 - pDamage) * pFailGivenIntact;
  const nothing = pfPrior * costFailure;
  const repair = costRepair + pFailAfterRepair * costFailure;
  const costWithout = Math.min(nothing, repair);
  const pInd = pDamage * pd + (1 - pDamage) * falseAlarm;
  const postInd = pInd > 0 ? (pDamage * pd) / pInd : 0;
  const postClear = pInd < 1 ? (pDamage * (1 - pd)) / (1 - pInd) : 0;
  const act = (pDam) => Math.min((pDam * pFailGivenDamage + (1 - pDam) * pFailGivenIntact) * costFailure, repair);
  const costWith = pInd * act(postInd) + (1 - pInd) * act(postClear);
  const voi = costWithout - costWith;
  return {
    costWithout, costWith, valueOfInformation: voi, netValue: voi - costInspection, worthwhile: voi > costInspection,
    actionWithout: nothing <= repair ? 'do-nothing' : 'repair', posteriorDamageGivenIndication: postInd, posteriorDamageGivenClear: postClear,
  };
}
