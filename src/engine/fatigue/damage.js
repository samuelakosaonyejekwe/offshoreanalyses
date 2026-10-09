// Cumulative fatigue damage: Palmgren–Miner summation, closed-form long-term damage, life and reliability.
// S–N constants are passed in the neutral form { m1, a1, m2, a2, sKnee } of fatigue/sn.js (stress in Pa).

import { erfc, erfinv, gamma, lgamma } from '../math/special.js';
import { cyclesFromConstants, weibullClosedFormDamage, weibullScaleFromMax } from './sn.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function isSeries(x) {
  return Array.isArray(x) || ArrayBuffer.isView(x);
}
const phi = (x) => 0.5 * erfc(-x / Math.SQRT2);

/**
 * Palmgren–Miner sum D = Σ nᵢ/Nᵢ for a stress-range histogram.
 * Reference: Palmgren (1924), Miner (1945); DNV-RP-C203 Eq. (2.2.1).
 * @param {object} o
 * @param {ArrayLike<number>} o.ranges stress ranges [Pa]
 * @param {ArrayLike<number>} o.counts number of cycles in each class [cycles]
 * @param {object} o.sn S–N constants
 * @returns {{damage:number, contributions: Float64Array}} damage [-] and per-class contributions [-]
 */
export function minerSum({ ranges, counts, sn }) {
  need(isSeries(ranges) && isSeries(counts) && ranges.length === counts.length, 'ranges and counts must be arrays of equal length');
  const contributions = new Float64Array(ranges.length);
  let damage = 0;
  for (let i = 0; i < ranges.length; i++) {
    need(counts[i] >= 0, `counts[${i}] must be non-negative`);
    const d = counts[i] === 0 || ranges[i] === 0 ? 0 : counts[i] / cyclesFromConstants(sn, ranges[i]);
    contributions[i] = d;
    damage += d;
  }
  return { damage, contributions };
}

/**
 * Miner sum of counted rainflow cycles ({range, count}).
 * @param {object[]} cycles output of fatigue/rainflow.js [range in Pa]
 * @param {object} sn S–N constants
 * @returns {number} damage [-]
 */
export function minerSumCycles(cycles, sn) {
  let d = 0;
  for (const c of cycles) if (c.range > 0) d += c.count / cyclesFromConstants(sn, c.range);
  return d;
}

/**
 * Miner sum of a block loading programme: each block has a constant stress range, a cycle count per repetition
 * and optionally its own S–N constants (e.g. different hot spots or environments).
 * @param {object} o
 * @param {{stressRange:number, cycles:number, sn?:object, label?:string}[]} o.blocks stressRange [Pa], cycles
 *   per repetition [cycles]
 * @param {object} [o.sn] default S–N constants
 * @param {number} [o.repetitions=1] number of times the programme is applied [-]
 * @returns {{damage:number, damagePerRepetition:number, repetitionsToFailure:number, blocks: {label:(string|undefined),
 *   damage:number}[]}}
 */
export function minerBlocks({ blocks, sn, repetitions = 1 }) {
  need(Array.isArray(blocks), 'blocks must be an array');
  need(repetitions >= 0, 'repetitions must be non-negative');
  let per = 0;
  const out = [];
  blocks.forEach((b, i) => {
    const c = b.sn || sn;
    need(c, `blocks[${i}] has no S-N constants`);
    need(b.cycles >= 0 && b.stressRange >= 0, `blocks[${i}] must have non-negative stressRange and cycles`);
    const d = b.stressRange === 0 ? 0 : b.cycles / cyclesFromConstants(c, b.stressRange);
    per += d;
    out.push({ label: b.label, damage: d * repetitions });
  });
  return { damage: per * repetitions, damagePerRepetition: per, repetitionsToFailure: per > 0 ? 1 / per : Infinity, blocks: out };
}

/**
 * Miner damage of a two-parameter Weibull long-term stress-range distribution (closed form, single or bilinear
 * S–N). Either the scale q or the maximum range Δσ0 in `cycles` cycles is given.
 * Reference: DNV-RP-C203 §5.1, Eq. (5.1.1)–(5.1.2).
 * @param {object} o
 * @param {object} o.sn S–N constants
 * @param {number} o.shape Weibull shape h [-]
 * @param {number} o.cycles number of cycles n0 [cycles]
 * @param {number} [o.scale] Weibull scale q [Pa]
 * @param {number} [o.maxStressRange] Δσ0, range exceeded once in n0 cycles [Pa]
 * @returns {{damage:number, scale:number}} damage [-], scale [Pa]
 */
export function minerWeibull({ sn, shape, cycles, scale, maxStressRange }) {
  need(scale !== undefined || maxStressRange !== undefined, 'scale or maxStressRange is required');
  const q = scale !== undefined ? scale : weibullScaleFromMax(maxStressRange, shape, cycles);
  return { damage: weibullClosedFormDamage({ sn, scale: q, shape, cycles }), scale: q };
}

/**
 * Discretise a Weibull long-term distribution into a histogram (equal-width classes, class mid-point ranges).
 * Used to cross-check closed forms and to feed block-based tools.
 * @param {object} o
 * @param {number} o.scale q [Pa]
 * @param {number} o.shape h [-]
 * @param {number} o.cycles n0 [cycles]
 * @param {number} [o.bins=200] number of classes [-]
 * @param {number} [o.maxRange] upper limit [Pa]; default = range exceeded with probability 1e-12
 * @returns {{ranges: Float64Array, counts: Float64Array}}
 */
export function weibullHistogram({ scale, shape, cycles, bins = 200, maxRange }) {
  need(scale > 0 && shape > 0, 'scale and shape must be positive');
  need(cycles >= 0, 'cycles must be non-negative');
  need(Number.isInteger(bins) && bins > 0, 'bins must be a positive integer');
  const top = maxRange === undefined ? scale * (-Math.log(1e-12)) ** (1 / shape) : maxRange;
  const w = top / bins;
  const ranges = new Float64Array(bins);
  const counts = new Float64Array(bins);
  const sf = (s) => Math.exp(-((s / scale) ** shape));
  for (let i = 0; i < bins; i++) {
    ranges[i] = (i + 0.5) * w;
    counts[i] = cycles * (sf(i * w) - sf((i + 1) * w));
  }
  return { ranges, counts };
}

/**
 * Fatigue life, utilisation and remaining life from a damage rate, with an accumulated-damage start state.
 * life = (η − D0)/(DFF·ḋ) measured from now, where η is the allowable Miner sum.
 * Reference: DNV-RP-C203; DNV-ST-0126 (design fatigue factors).
 * @param {object} o
 * @param {number} o.damageRate characteristic damage per year (without DFF) [1/yr]
 * @param {number} [o.designLife] intended service life from now [yr]
 * @param {number} [o.dff=1] design fatigue factor [-]
 * @param {number} [o.accumulatedDamage=0] characteristic damage already consumed (without DFF) [-]
 * @param {number} [o.allowableDamage=1] Miner sum at failure [-]
 * @returns {{characteristicLife:number, designLife:number, remainingLife:number, damageAtDesignLife:(number|null),
 *   utilisation:(number|null), ok:(boolean|null)}} lives in years: characteristicLife = η/ḋ from new,
 *   designLife = η/(DFF·ḋ) from new, remainingLife from now including DFF and the consumed damage;
 *   utilisation = DFF·(D0 + ḋ·designLife)/η
 */
export function fatigueLife({ damageRate, designLife, dff = 1, accumulatedDamage = 0, allowableDamage = 1 }) {
  need(damageRate >= 0, 'damageRate must be non-negative');
  need(dff > 0, 'dff must be positive');
  need(accumulatedDamage >= 0, 'accumulatedDamage must be non-negative');
  need(allowableDamage > 0, 'allowableDamage must be positive');
  const inf = damageRate === 0;
  const characteristicLife = inf ? Infinity : allowableDamage / damageRate;
  const remaining = inf ? Infinity : Math.max(0, (allowableDamage / dff - accumulatedDamage) / damageRate);
  let damageAtDesignLife = null;
  let utilisation = null;
  let ok = null;
  if (designLife !== undefined) {
    need(designLife >= 0, 'designLife must be non-negative');
    damageAtDesignLife = accumulatedDamage + damageRate * designLife;
    utilisation = (dff * damageAtDesignLife) / allowableDamage;
    ok = utilisation <= 1;
  }
  return { characteristicLife, designLife: characteristicLife / dff, remainingLife: remaining, damageAtDesignLife, utilisation, ok };
}

/**
 * Combine damage from several independent sources or phases (wave, VIV, installation, transport, …) by linear,
 * sequence-insensitive summation; each source may carry its own design fatigue factor.
 * Reference: DNV-RP-C203 §2.2; DNV-ST-F101 / DNV-RP-F204 (phase-wise fatigue budgets).
 * @param {{name:string, damage:number, dff?:number}[]} sources damage of each source over its exposure [-]
 * @param {object} [o]
 * @param {number} [o.allowableDamage=1] [-]
 * @returns {{damage:number, factoredDamage:number, utilisation:number, ok:boolean, shares: {name:string,
 *   damage:number, factoredDamage:number, share:number}[]}}
 */
export function combineDamage(sources, { allowableDamage = 1 } = {}) {
  need(Array.isArray(sources), 'sources must be an array');
  need(allowableDamage > 0, 'allowableDamage must be positive');
  let d = 0;
  let fd = 0;
  for (const s of sources) {
    need(s.damage >= 0, `sources["${s.name}"].damage must be non-negative`);
    const dff = s.dff === undefined ? 1 : s.dff;
    need(dff > 0, `sources["${s.name}"].dff must be positive`);
    d += s.damage;
    fd += s.damage * dff;
  }
  const shares = sources.map((s) => {
    const f = s.damage * (s.dff === undefined ? 1 : s.dff);
    return { name: s.name, damage: s.damage, factoredDamage: f, share: fd > 0 ? f / fd : 0 };
  });
  return { damage: d, factoredDamage: fd, utilisation: fd / allowableDamage, ok: fd <= allowableDamage, shares };
}

/**
 * Nonlinear sequence-dependent damage by the Manson–Halford damage-curve approach. Moving from a block with
 * life N₁ to one with life N₂ the consumed life fraction r transforms as r → r^((N₁/N₂)^0.4); failure when the
 * fraction reaches 1. High-to-low sequences fail earlier than Miner, low-to-high later.
 * Reference: Manson & Halford (1981) Int. J. Fracture 17(2), 169–192.
 * @param {{cycles:number, cyclesToFailure:number}[]} blocks blocks in the order applied
 * @returns {{lifeFraction:number, failed:boolean, failedAtBlock:(number|null), minerSum:number}}
 */
export function mansonHalfordDamage(blocks) {
  need(Array.isArray(blocks) && blocks.length > 0, 'blocks must be a non-empty array');
  let r = 0;
  let miner = 0;
  let prevN = null;
  for (let i = 0; i < blocks.length; i++) {
    const { cycles, cyclesToFailure: N } = blocks[i];
    need(cycles >= 0, `blocks[${i}].cycles must be non-negative`);
    need(N > 0, `blocks[${i}].cyclesToFailure must be positive`);
    if (prevN !== null && r > 0) r = r ** ((prevN / N) ** 0.4);
    r += cycles / N;
    miner += cycles / N;
    prevN = N;
    if (r >= 1) return { lifeFraction: r, failed: true, failedAtBlock: i, minerSum: miner };
  }
  return { lifeFraction: r, failed: false, failedAtBlock: null, minerSum: miner };
}

/**
 * Extrapolate a measured range histogram to a longer exposure. Counts scale with duration; the tail is described
 * by a two-parameter Weibull fitted by the method of moments to the measured ranges, giving the expected largest
 * range in the target exposure and the closed-form damage of the fitted distribution.
 * Reference: DNV-RP-C203 §5.1 (Weibull long-term distribution); Johannesson (2006) for extrapolation practice.
 * Validity: measurement period representative of the long-term loading (stationary mix of conditions).
 * @param {object} o
 * @param {ArrayLike<number>} o.ranges class ranges [Pa]
 * @param {ArrayLike<number>} o.counts measured counts [cycles]
 * @param {number} o.measuredDuration [s]
 * @param {number} o.targetDuration [s]
 * @param {object} [o.sn] S–N constants; when given, damage is returned
 * @returns {{counts: Float64Array, cycles:number, shape:number, scale:number, maxRangeMeasured:number,
 *   maxRangeExpected:number, damageScaled:(number|null), damageWeibull:(number|null)}} ranges in Pa
 */
export function extrapolateHistogram({ ranges, counts, measuredDuration, targetDuration, sn }) {
  need(isSeries(ranges) && isSeries(counts) && ranges.length === counts.length && ranges.length > 0, 'ranges and counts must be non-empty arrays of equal length');
  need(measuredDuration > 0 && targetDuration > 0, 'durations must be positive');
  const factor = targetDuration / measuredDuration;
  const scaled = new Float64Array(counts.length);
  let n = 0;
  let s1 = 0;
  let s2 = 0;
  let maxMeasured = 0;
  for (let i = 0; i < counts.length; i++) {
    scaled[i] = counts[i] * factor;
    n += counts[i];
    s1 += counts[i] * ranges[i];
    s2 += counts[i] * ranges[i] * ranges[i];
    if (counts[i] > 0 && ranges[i] > maxMeasured) maxMeasured = ranges[i];
  }
  need(n > 0 && s1 > 0, 'counts must contain cycles');
  const mean = s1 / n;
  const cv2 = Math.max(s2 / n / (mean * mean) - 1, 1e-8);
  // Method of moments: Γ(1+2/h)/Γ(1+1/h)² = 1 + CoV², solved by bisection (monotonic decreasing in h).
  let lo = 0.05;
  let hi = 200;
  for (let i = 0; i < 200; i++) {
    const h = Math.sqrt(lo * hi);
    const r = Math.exp(lgamma(1 + 2 / h) - 2 * lgamma(1 + 1 / h)) - 1;
    if (r > cv2) lo = h; else hi = h;
  }
  const shape = Math.sqrt(lo * hi);
  const scale = mean / gamma(1 + 1 / shape);
  const nTarget = n * factor;
  return {
    counts: scaled, cycles: nTarget, shape, scale, maxRangeMeasured: maxMeasured,
    maxRangeExpected: nTarget > 1 ? scale * Math.log(nTarget) ** (1 / shape) : maxMeasured,
    damageScaled: sn ? minerSum({ ranges, counts: scaled, sn }).damage : null,
    damageWeibull: sn ? weibullClosedFormDamage({ sn, scale, shape, cycles: nTarget }) : null,
  };
}

/**
 * Fatigue reliability by the lognormal format (Wirsching): time to failure T = Δ·T_c/B^m, where T_c is the
 * life computed with the median S–N curve (Miner sum 1), Δ the Miner sum at failure, B the stress-modelling
 * error and the scatter of the S–N intercept a is carried as its own lognormal factor. All variables lognormal:
 *   σ_lnT² = ln[(1 + CΔ²)(1 + C_a²)(1 + C_B²)^(m²)],   β(t) = ln(T̃/t)/σ_lnT,   P_f(t) = Φ(−β(t)).
 * The annual failure probability is the increment of P_f in the year ending at t; the hazard is that increment
 * divided by the survival probability at the start of the year.
 * Reference: Wirsching (1984) J. Struct. Eng. ASCE 110(10), 2340–2356; Wirsching & Chen (1988).
 * Typical values (Wirsching 1984): Δ̃ = 1.0, CΔ = 0.30; B̃ = 0.7–1.0, C_B = 0.15–0.30; C_a = 0.4–0.7 (welded joints).
 * Validity: single-slope S–N curve.
 * @param {object} o
 * @param {number} o.medianLife life from the median S–N curve with unit Miner sum and unbiased stress [yr]
 * @param {number} o.m S–N inverse slope [-]
 * @param {ArrayLike<number>} o.times evaluation times [yr], positive
 * @param {number} [o.deltaMedian=1] median Miner sum at failure [-]
 * @param {number} [o.deltaCov=0.3] [-]
 * @param {number} [o.aCov=0.5] coefficient of variation of the S–N intercept [-]
 * @param {number} [o.bMedian=1] median stress bias (actual/computed) [-]
 * @param {number} [o.bCov=0.25] [-]
 * @returns {{medianTimeToFailure:number, sigmaLnT:number, beta: Float64Array, pf: Float64Array,
 *   annualPf: Float64Array, hazard: Float64Array, annualBeta: Float64Array}} times in yr, probabilities [-]
 */
export function fatigueReliabilityLognormal({ medianLife, m, times, deltaMedian = 1, deltaCov = 0.3, aCov = 0.5, bMedian = 1, bCov = 0.25 }) {
  need(medianLife > 0, 'medianLife must be positive');
  need(m > 0, 'm must be positive');
  need(isSeries(times), 'times must be an array');
  need(deltaMedian > 0 && bMedian > 0, 'deltaMedian and bMedian must be positive');
  need(deltaCov >= 0 && aCov >= 0 && bCov >= 0, 'coefficients of variation must be non-negative');
  const sig = Math.sqrt(Math.log(1 + deltaCov ** 2) + Math.log(1 + aCov ** 2) + m * m * Math.log(1 + bCov ** 2));
  need(sig > 0, 'at least one coefficient of variation must be positive');
  const tMed = (deltaMedian * medianLife) / bMedian ** m;
  const n = times.length;
  const beta = new Float64Array(n);
  const pf = new Float64Array(n);
  const annualPf = new Float64Array(n);
  const hazard = new Float64Array(n);
  const annualBeta = new Float64Array(n);
  const P = (t) => (t <= 0 ? 0 : phi(-Math.log(tMed / t) / sig));
  for (let i = 0; i < n; i++) {
    const t = times[i];
    need(t > 0, `times[${i}] must be positive`);
    beta[i] = Math.log(tMed / t) / sig;
    pf[i] = P(t);
    const prev = P(t - 1);
    annualPf[i] = Math.max(0, pf[i] - prev);
    hazard[i] = prev < 1 ? annualPf[i] / (1 - prev) : 1;
    annualBeta[i] = annualPf[i] > 0 && annualPf[i] < 1 ? -Math.SQRT2 * erfinv(2 * annualPf[i] - 1) : annualPf[i] <= 0 ? Infinity : -Infinity;
  }
  return { medianTimeToFailure: tMed, sigmaLnT: sig, beta, pf, annualPf, hazard, annualBeta };
}

/**
 * Median-to-design life ratio implied by a characteristic S–N curve placed k standard deviations below the mean
 * in log N: T_median = T_characteristic·10^(k·s).
 * Reference: DNV-RP-C203 §2.4.1 (design curve = mean − 2·s_logN, s_logN = 0.20 for the tabulated curves).
 * @param {number} characteristicLife life from the design S–N curve [yr]
 * @param {object} [o]
 * @param {number} [o.sLogN=0.20] standard deviation of log10 N [-]
 * @param {number} [o.k=2] number of standard deviations [-]
 * @returns {{medianLife:number, aCov:number}} median life [yr] and the equivalent CoV of the intercept a [-]
 */
export function medianLifeFromCharacteristic(characteristicLife, { sLogN = 0.20, k = 2 } = {}) {
  need(characteristicLife > 0, 'characteristicLife must be positive');
  need(sLogN >= 0, 'sLogN must be non-negative');
  const sLn = sLogN * Math.LN10;
  return { medianLife: characteristicLife * 10 ** (k * sLogN), aCov: Math.sqrt(Math.exp(sLn * sLn) - 1) };
}
