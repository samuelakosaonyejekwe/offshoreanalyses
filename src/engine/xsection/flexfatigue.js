/**
 * Fatigue, wear and end-fitting checks for the armour layers and polymer sheaths of unbonded flexible pipes.
 *
 * Stress model per bending cycle (Saevik 2011; Feret & Bournazel 1987, see helix.js): the wire axial stress
 * follows the pipe in stick until the friction capacity is exhausted (full slip), after which only the local wire
 * bending stresses keep growing with curvature. Stresses are in Pa, S-N curves are N = a S^(-m) with S in Pa.
 */

import { criticalCurvature, wireCurvatureChange } from './helix.js';
import { axisymmetricResponse } from './flexible.js';

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function fin(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Stress ranges in a tensile-armour wire for one bending cycle of curvature range dKappa (peak to peak).
 *  - friction (stick) stress range at the extrados/intrados: min(E cos^2(a) R dKappa, 2 sigma_f,max), with
 *    sigma_f,max = pi R f /(2 sin(a) A) (Masing rule: the range response is twice the amplitude response);
 *  - local bending ranges: E (t/2) |dKn| at psi = 0 and E (b/2) |dKt| at psi = pi/2 (helix.js
 *    `wireCurvatureChange`).
 * Two hot spots are reported: extrados (friction + normal-curvature bending) and neutral axis (transverse
 * bending only, friction stress zero there).
 * Reference: Saevik (2011) Marine Structures 24, 97-111.
 * Validity: constant contact pressure through the cycle; R dKappa << 1.
 * @param {object} o
 * @param {number} o.E Pa
 * @param {number} o.A wire area, m^2
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.thickness radial wire thickness, m
 * @param {number} o.width wire width, m
 * @param {number} o.frictionPerLength friction per unit wire length, N/m
 * @param {number} o.kappaRange curvature range, 1/m
 * @param {'loxodromic'|'geodesic'} [o.curve='loxodromic']
 * @returns {{frictionRange:number, stickRange:number, slipped:boolean, normalBendingRange:number,
 *   transverseBendingRange:number, rangeExtrados:number, rangeNeutralAxis:number, governing:number}} Pa
 */
export function armourStressRange({ E, A, R, alpha, thickness, width, frictionPerLength, kappaRange, curve = 'loxodromic' }) {
  pos(thickness, 'thickness'); pos(width, 'width'); nonneg(kappaRange, 'kappaRange');
  const cc = criticalCurvature({ E, A, R, alpha, frictionPerLength });
  const stickRange = E * Math.cos(alpha) ** 2 * R * kappaRange;
  const frictionRange = Math.min(stickRange, 2 * cc.sigmaFrictionMax);
  const k0 = wireCurvatureChange({ alpha, kappa: kappaRange, psi: 0, curve });
  const k90 = wireCurvatureChange({ alpha, kappa: kappaRange, psi: Math.PI / 2, curve });
  const normalBendingRange = E * thickness / 2 * Math.abs(k0.dKn);
  const transverseBendingRange = E * width / 2 * Math.abs(k90.dKt);
  const rangeExtrados = frictionRange + normalBendingRange;
  return {
    frictionRange, stickRange, slipped: stickRange > 2 * cc.sigmaFrictionMax,
    normalBendingRange, transverseBendingRange, rangeExtrados, rangeNeutralAxis: transverseBendingRange,
    governing: Math.max(rangeExtrados, transverseBendingRange),
  };
}

/**
 * Mean-stress correction of a stress range to the equivalent fully reversed range.
 *   Goodman: S_eq = S / (1 - sigma_m/sigma_u);   Gerber: S_eq = S / (1 - (sigma_m/sigma_u)^2).
 * Compressive mean stress is not credited (treated as zero).
 * Reference: Goodman (1899), Gerber (1874); e.g. Dowling, Mechanical Behavior of Materials, §9.7.
 * Validity: 0 <= sigma_m < sigma_u.
 * @param {object} o
 * @param {number} o.stressRange Pa
 * @param {number} o.meanStress Pa
 * @param {number} o.uts ultimate tensile strength, Pa
 * @param {'goodman'|'gerber'|'none'} [o.method='goodman']
 * @returns {{equivalentRange:number, factor:number}} Pa, -
 */
export function meanStressCorrection({ stressRange, meanStress, uts, method = 'goodman' }) {
  nonneg(stressRange, 'stressRange'); fin(meanStress, 'meanStress');
  if (method === 'none') return { equivalentRange: stressRange, factor: 1 };
  pos(uts, 'uts');
  const r = Math.max(meanStress, 0) / uts;
  if (r >= 1) throw new RangeError('meanStress must be below uts');
  let factor;
  if (method === 'goodman') factor = 1 / (1 - r);
  else if (method === 'gerber') factor = 1 / (1 - r * r);
  else throw new RangeError("method must be 'goodman', 'gerber' or 'none'");
  return { equivalentRange: stressRange * factor, factor };
}

/**
 * Cycles to failure from an S-N curve N = a S^(-m) (S = stress range in Pa), optionally bilinear (second slope m2
 * beyond nKnee cycles, continuous at the knee) and with a cut-off range `limit` below which life is infinite.
 * Reference: Basquin relation; format of DNV-RP-C203 §2.4.
 * @param {object} o
 * @param {number} o.stressRange Pa
 * @param {{a:number, m:number, m2?:number, nKnee?:number, limit?:number}} o.curve a in Pa^m cycles
 * @returns {{cycles:number}} cycles to failure (Infinity below the cut-off or at zero range)
 */
export function snCycles({ stressRange, curve }) {
  nonneg(stressRange, 'stressRange');
  if (!curve) throw new RangeError('curve is required');
  pos(curve.a, 'curve.a'); pos(curve.m, 'curve.m');
  if (stressRange === 0 || (curve.limit !== undefined && stressRange <= curve.limit)) return { cycles: Infinity };
  // work in logarithms: a is of order 1e100 for S in Pa
  const la = Math.log10(curve.a), ls = Math.log10(stressRange);
  let logN = la - curve.m * ls;
  if (curve.m2 !== undefined && curve.nKnee !== undefined) {
    const lk = Math.log10(pos(curve.nKnee, 'curve.nKnee'));
    if (logN > lk) {
      const lsKnee = (la - lk) / curve.m;
      logN = lk + pos(curve.m2, 'curve.m2') * (lsKnee - ls);
    }
  }
  return { cycles: 10 ** logN };
}

/**
 * GENERIC, UNQUALIFIED armour-wire S-N curve for demonstration only. It is the textbook estimate for polished
 * steel specimens in air (Shigley's Mechanical Engineering Design, §6-7/6-8): fully reversed strength
 * f x S_ut at 10^3 cycles (f = 0.9) and endurance limit 0.5 S_ut (capped at 700 MPa) at 10^6 cycles, joined by a
 * Basquin line; expressed here in stress range (twice the amplitude). `marin` multiplies the endurance limit for
 * surface, size and reliability effects (Shigley Eq. 6-18).
 * Corrosive annulus environments: there is no generic curve that can be justified without test data. For
 * 'sweet' (CO2) and 'sour' (H2S) annulus conditions the caller must supply environmentFactor >= 1, a knock-down
 * applied to the stress range axis (i.e. strengths divided by it), and the endurance limit is removed (no fatigue
 * limit in corrosion fatigue). Project fatigue design must use S-N data from tests on the actual wire in the
 * actual annulus environment (API Spec 17J §5.3.4 / API RP 17B).
 * @param {object} o
 * @param {number} o.uts wire ultimate tensile strength, Pa
 * @param {'air'|'sweet'|'sour'} [o.environment='air']
 * @param {number} [o.environmentFactor] stress knock-down >= 1, required for 'sweet' and 'sour'
 * @param {number} [o.f=0.9] fraction of S_ut at 10^3 cycles
 * @param {number} [o.marin=1] combined Marin modification factor on the endurance limit
 * @returns {{a:number, m:number, limit?:number, generic:true, label:string}} curve for {@link snCycles}
 */
export function genericArmourSN({ uts, environment = 'air', environmentFactor, f = 0.9, marin = 1 }) {
  pos(uts, 'uts'); pos(f, 'f'); pos(marin, 'marin');
  let kd = 1;
  if (environment === 'sweet' || environment === 'sour') {
    if (!(environmentFactor >= 1)) throw new RangeError(`environmentFactor (>= 1) is required for a ${environment} annulus: no generic value exists`);
    kd = environmentFactor;
  } else if (environment !== 'air') throw new RangeError("environment must be 'air', 'sweet' or 'sour'");
  const s3 = 2 * f * uts / kd; // range at 1e3 cycles
  const s6 = 2 * marin * Math.min(0.5 * uts, 700e6) / kd; // range at 1e6 cycles
  const m = 3 / Math.log10(s3 / s6);
  const a = 10 ** (6 + m * Math.log10(s6));
  const curve = { a, m, generic: true, label: `GENERIC unqualified ${environment} curve (textbook estimate, demonstration only)` };
  if (environment === 'air') curve.limit = s6;
  return curve;
}

/**
 * Palmgren-Miner damage for a table of stress-range bins.
 * Each bin: { stressRange [Pa], cycles, meanStress [Pa] (optional) }. The range is multiplied by the stress
 * concentration factor, mean-stress corrected, then N from the S-N curve; cycles to failure are divided by
 * corrosionFactor (>= 1, life reduction for corrosion fatigue when the curve was obtained in air).
 * Service life = exposure / (damage x DFF); API Spec 17J requires a design fatigue factor of 10 on the armour
 * service life.
 * Reference: Palmgren (1924), Miner (1945); API Spec 17J §5.3.4.
 * @param {object} o
 * @param {{stressRange:number, cycles:number, meanStress?:number}[]} o.bins cycles counted over `exposureYears`
 * @param {object} o.curve S-N curve for {@link snCycles}
 * @param {number} [o.uts] Pa, needed for a mean-stress correction
 * @param {'goodman'|'gerber'|'none'} [o.meanStressMethod='none']
 * @param {number} [o.scf=1] stress concentration factor on the range (e.g. at the end-fitting entrance)
 * @param {number} [o.corrosionFactor=1] divisor on cycles to failure
 * @param {number} [o.dff=10] design fatigue factor
 * @param {number} [o.exposureYears=1] duration represented by the bin cycle counts, years
 * @returns {{damage:number, damagePerYear:number, life:number, lifeWithDff:number, bins:{equivalentRange:number,
 *   cyclesToFailure:number, damage:number}[]}} -, 1/yr, years, years
 */
export function minerDamage({ bins, curve, uts, meanStressMethod = 'none', scf = 1, corrosionFactor = 1, dff = 10, exposureYears = 1 }) {
  if (!Array.isArray(bins)) throw new RangeError('bins must be an array');
  pos(scf, 'scf'); pos(dff, 'dff'); pos(exposureYears, 'exposureYears');
  if (!(corrosionFactor >= 1)) throw new RangeError('corrosionFactor must be >= 1');
  let damage = 0;
  const rows = bins.map((b, i) => {
    nonneg(b.cycles, `bins[${i}].cycles`);
    const eq = meanStressCorrection({ stressRange: scf * nonneg(b.stressRange, `bins[${i}].stressRange`), meanStress: b.meanStress ?? 0, uts, method: meanStressMethod }).equivalentRange;
    const N = snCycles({ stressRange: eq, curve }).cycles / corrosionFactor;
    const d = b.cycles / N;
    damage += d;
    return { equivalentRange: eq, cyclesToFailure: N, damage: d };
  });
  const perYear = damage / exposureYears;
  return { damage, damagePerYear: perYear, life: perYear > 0 ? 1 / perYear : Infinity, lifeWithDff: perYear > 0 ? 1 / (perYear * dff) : Infinity, bins: rows };
}

/**
 * Linearised fatigue coefficients of one tensile-armour layer from two axisymmetric solutions of the flexible
 * pipe (flexible.js): wire mean stress and friction per unit wire length as affine functions of effective
 * tension, sigma = s0 + s1 T, f = f0 + f1 T, with f = mu (q_in + q_out) 2 pi R cos(alpha)/n.
 * Validity: contact state unchanged between the two tensions.
 * @param {object} stack layer stack from buildFlexibleStack
 * @param {number} layerIndex index of the tensile-armour layer
 * @param {object} o
 * @param {number} o.tensionLow effective tension, N
 * @param {number} o.tensionHigh effective tension, N
 * @param {number} [o.mu=0.1] friction coefficient
 * @param {object} [o.loads] other loads (pInt, pExt, annulus ...)
 * @returns {{meanStress0:number, meanStressPerTension:number, friction0:number, frictionPerTension:number,
 *   wire:{E:number, A:number, R:number, alpha:number, thickness:number, width:number|null}}} Pa, Pa/N, N/m,
 *   (N/m)/N; wire data for {@link armourFatigueDamage} (per-wire friction uses the layer wire count, or one
 *   equivalent wire of the whole layer area when the count is unknown - the stresses are identical)
 */
export function armourFatigueCoefficients(stack, layerIndex, { tensionLow, tensionHigh, mu = 0.1, loads = {} }) {
  const L = stack.layers[layerIndex];
  if (!L || L.kind !== 'helical') throw new RangeError('layerIndex must refer to a helical layer');
  fin(tensionLow, 'tensionLow'); fin(tensionHigh, 'tensionHigh'); nonneg(mu, 'mu');
  if (tensionHigh === tensionLow) throw new RangeError('tensionHigh must differ from tensionLow');
  const n = L.n ?? 1;
  const at = (T) => {
    const r = axisymmetricResponse(stack, { ...loads, tension: T }).layers[layerIndex];
    const q = Math.max(r.contactIn, 0) + Math.max(r.contactOut, 0);
    return { s: r.wireStress, f: mu * q * 2 * Math.PI * L.R * Math.cos(L.alpha) / n };
  };
  const a = at(tensionLow), b = at(tensionHigh), dT = tensionHigh - tensionLow;
  const s1 = (b.s - a.s) / dT, f1 = (b.f - a.f) / dT;
  return {
    meanStress0: a.s - s1 * tensionLow, meanStressPerTension: s1, friction0: a.f - f1 * tensionLow, frictionPerTension: f1,
    wire: { E: L.E, A: L.nA / n, R: L.R, alpha: L.alpha, thickness: L.wireThickness ?? L.t, width: L.wireWidth ?? null },
  };
}

/**
 * Pressure-armour (or any helical layer) stress range for an internal-pressure cycle, from two axisymmetric
 * solutions of the flexible pipe: range = |sigma(pHigh) - sigma(pLow)|, mean = their average. Feed the result to
 * {@link minerDamage} with the pressure-armour S-N curve. Bending-induced ovalisation stresses of the pressure
 * armour profile are not modelled.
 * Reference: API RP 17B fatigue methodology (pressure cycling of the pressure armour); flexible.js model.
 * @param {object} stack layer stack from buildFlexibleStack
 * @param {number} layerIndex index of the helical layer
 * @param {object} o
 * @param {number} o.pLow bore pressure at the bottom of the cycle, Pa
 * @param {number} o.pHigh bore pressure at the top of the cycle, Pa
 * @param {object} [o.loads] other loads held constant (tension, pExt, annulus ...)
 * @returns {{stressRange:number, meanStress:number, stressLow:number, stressHigh:number}} Pa
 */
export function pressureCycleStressRange(stack, layerIndex, { pLow, pHigh, loads = {} }) {
  const L = stack.layers[layerIndex];
  if (!L || L.kind !== 'helical') throw new RangeError('layerIndex must refer to a helical layer');
  nonneg(pLow, 'pLow'); nonneg(pHigh, 'pHigh');
  const lo = axisymmetricResponse(stack, { ...loads, pInt: pLow }).layers[layerIndex].wireStress;
  const hi = axisymmetricResponse(stack, { ...loads, pInt: pHigh }).layers[layerIndex].wireStress;
  return { stressRange: Math.abs(hi - lo), meanStress: 0.5 * (hi + lo), stressLow: lo, stressHigh: hi };
}

/**
 * Tensile-armour fatigue damage over a table of (tension, curvature range, cycles) load bins.
 * For each bin: mean stress = meanStress0 + meanStressPerTension T; friction f = friction0 + frictionPerTension T
 * (>= 0); stress range from {@link armourStressRange} (governing hot spot); mean-stress correction; Miner sum via
 * {@link minerDamage}.
 * Reference: API RP 17B §8 fatigue methodology (global load bins -> local wire stress -> Miner), Saevik (2011).
 * @param {object} o
 * @param {{E:number, A:number, R:number, alpha:number, thickness:number, width:number}} o.wire SI
 * @param {{tension:number, curvatureRange:number, cycles:number}[]} o.bins N, 1/m, cycles in exposureYears
 * @param {number} [o.meanStress0=0] Pa
 * @param {number} o.meanStressPerTension Pa/N
 * @param {number} [o.friction0=0] N/m
 * @param {number} o.frictionPerTension (N/m)/N
 * @param {object} o.curve S-N curve
 * @param {number} [o.uts] Pa
 * @param {'goodman'|'gerber'|'none'} [o.meanStressMethod='goodman']
 * @param {number} [o.scf=1]
 * @param {number} [o.corrosionFactor=1]
 * @param {number} [o.dff=10]
 * @param {number} [o.exposureYears=1]
 * @param {'loxodromic'|'geodesic'} [o.wireCurve='loxodromic']
 * @returns {object} minerDamage result plus stress rows [{ meanStress, frictionPerLength, range, slipped }]
 */
export function armourFatigueDamage({ wire, bins, meanStress0 = 0, meanStressPerTension, friction0 = 0, frictionPerTension, curve, uts, meanStressMethod = 'goodman', scf = 1, corrosionFactor = 1, dff = 10, exposureYears = 1, wireCurve = 'loxodromic' }) {
  if (!wire) throw new RangeError('wire is required');
  if (!Array.isArray(bins)) throw new RangeError('bins must be an array');
  fin(meanStressPerTension, 'meanStressPerTension'); fin(frictionPerTension, 'frictionPerTension');
  const stress = [];
  const sBins = bins.map((b, i) => {
    fin(b.tension, `bins[${i}].tension`);
    const f = Math.max(friction0 + frictionPerTension * b.tension, 0);
    const r = armourStressRange({ ...wire, frictionPerLength: f, kappaRange: b.curvatureRange, curve: wireCurve });
    const meanStress = meanStress0 + meanStressPerTension * b.tension;
    stress.push({ meanStress, frictionPerLength: f, range: r.governing, slipped: r.slipped });
    return { stressRange: r.governing, cycles: b.cycles, meanStress };
  });
  return { ...minerDamage({ bins: sBins, curve, uts, meanStressMethod, scf, corrosionFactor, dff, exposureYears }), stress };
}

/**
 * Relative slip between an armour wire and the adjacent layer per bending cycle.
 * Derivation: a fully slipped wire carries no bending-induced axial strain, so its displacement relative to the
 * pipe is the integral of the stick strain cos^2(a) R kappa cos(psi) along the wire (ds = R dpsi / sin(a)):
 * delta(psi) = R^2 cos^2(a) kappa sin(psi) / sin(a), maximum at the neutral axis. For a curvature range dKappa the
 * slip range is R^2 cos^2(a) dKappa / sin(a) and the sliding distance per cycle twice that (out and back). The
 * stick phase before slip onset is ignored, which overestimates sliding (conservative for wear).
 * Reference: Feret & Bournazel (1987).
 * @param {object} o
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.kappaRange curvature range, 1/m
 * @returns {{slipRange:number, slidingPerCycle:number}} m, m
 */
export function wireSlipPerCycle({ R, alpha, kappaRange }) {
  pos(R, 'R'); nonneg(kappaRange, 'kappaRange');
  if (!(alpha > 0 && alpha < Math.PI / 2)) throw new RangeError('alpha must lie strictly between 0 and pi/2 rad');
  const slipRange = R * R * Math.cos(alpha) ** 2 * kappaRange / Math.sin(alpha);
  return { slipRange, slidingPerCycle: 2 * slipRange };
}

/**
 * Archard wear of an armour wire / anti-wear tape interface and the resulting stress increase.
 * Archard: worn volume V = k F s, i.e. wear depth h = k p s with k the dimensional wear coefficient (m^2/N =
 * 1/Pa), p the contact pressure on the wire footprint and s the accumulated sliding distance. The wire stress
 * rises by the area-loss factor t/(t - h).
 * Reference: Archard (1953) J. Appl. Phys. 24, 981-988. k must come from tests on the actual material pair
 * (reported dimensional values for lubricated/dry steel pairs span roughly 1e-18 to 1e-13 m^2/N).
 * @param {object} o
 * @param {number} o.wearCoefficient k, m^2/N
 * @param {number} o.contactPressure p on the wire footprint (smeared pressure / fill fraction), Pa
 * @param {number} o.slidingPerCycle m
 * @param {number} o.cycles number of cycles
 * @param {number} o.thickness original wire thickness, m
 * @returns {{wearDepth:number, remainingThickness:number, stressFactor:number, wornThrough:boolean}} m, m, -
 */
export function archardWear({ wearCoefficient, contactPressure, slidingPerCycle, cycles, thickness }) {
  nonneg(wearCoefficient, 'wearCoefficient'); nonneg(contactPressure, 'contactPressure'); nonneg(slidingPerCycle, 'slidingPerCycle');
  nonneg(cycles, 'cycles'); pos(thickness, 'thickness');
  const wearDepth = wearCoefficient * contactPressure * slidingPerCycle * cycles;
  const rem = thickness - wearDepth;
  return { wearDepth, remainingThickness: Math.max(rem, 0), stressFactor: rem > 0 ? thickness / rem : Infinity, wornThrough: rem <= 0 };
}

/**
 * Loss of pressure-sheath thickness by creep of the polymer into the gaps of the adjacent armour layer (volume
 * balance): material filling a gap of width g to a depth d over each pitch Lp of the layer is drawn from the
 * sheath wall, dt = phi g d / Lp, with phi the fraction of the gap section actually filled (shape factor, 1 for a
 * rectangular plug). The intrusion depth depends on temperature, pressure and time and is an input (API Spec 17J
 * limits the reduction of pressure-sheath wall thickness due to creep into gaps to 30 % - confirm against the
 * governing edition).
 * @param {object} o
 * @param {number} o.gapWidth g, m
 * @param {number} o.intrusionDepth d, m
 * @param {number} o.pitch armour pitch Lp, m
 * @param {number} o.sheathThickness original thickness, m
 * @param {number} [o.fillShape=1] phi, -
 * @param {number} [o.allowableReduction=0.3] allowable fractional thickness reduction
 * @returns {{thicknessLoss:number, reduction:number, utilisation:number}} m, -, -
 */
export function sheathCreepIntoGaps({ gapWidth, intrusionDepth, pitch, sheathThickness, fillShape = 1, allowableReduction = 0.3 }) {
  nonneg(gapWidth, 'gapWidth'); nonneg(intrusionDepth, 'intrusionDepth'); pos(pitch, 'pitch'); pos(sheathThickness, 'sheathThickness');
  pos(fillShape, 'fillShape'); pos(allowableReduction, 'allowableReduction');
  if (gapWidth > pitch) throw new RangeError('gapWidth cannot exceed pitch');
  const thicknessLoss = fillShape * gapWidth * intrusionDepth / pitch;
  const reduction = thicknessLoss / sheathThickness;
  return { thicknessLoss, reduction, utilisation: reduction / allowableReduction };
}

/**
 * Polymer sheath bending strain check: maximum strain r_out kappa_max against the allowable strain, and strain
 * range r_out dKappa for a fatigue assessment against a caller-supplied allowable range.
 * Reference: API Spec 17J §5.3.1 (allowable sheath strains, see flexible.js SHEATH_STRAIN_LIMITS).
 * @param {object} o
 * @param {number} o.rOut sheath outer radius, m
 * @param {number} o.kappaMax maximum curvature, 1/m
 * @param {number} o.strainLimit allowable strain, -
 * @param {number} [o.kappaRange] cyclic curvature range, 1/m
 * @param {number} [o.strainRangeLimit] allowable cyclic strain range, -
 * @returns {{strain:number, utilisation:number, strainRange:number|null, rangeUtilisation:number|null}}
 */
export function sheathStrainCheck({ rOut, kappaMax, strainLimit, kappaRange, strainRangeLimit }) {
  pos(rOut, 'rOut'); nonneg(kappaMax, 'kappaMax'); pos(strainLimit, 'strainLimit');
  const strain = rOut * kappaMax;
  const strainRange = kappaRange === undefined ? null : rOut * nonneg(kappaRange, 'kappaRange');
  return {
    strain, utilisation: strain / strainLimit, strainRange,
    rangeUtilisation: strainRange !== null && strainRangeLimit !== undefined ? strainRange / pos(strainRangeLimit, 'strainRangeLimit') : null,
  };
}

/**
 * Armour-wire anchorage capacity in an epoxy-filled end fitting.
 *  - Bond: pull-out force per wire P_bond = tau_b x perimeter x embedded length (uniform bond stress).
 *  - Cone wedge: the epoxy plug of half-angle beta is drawn into the conical end-fitting body by the armour pull
 *    P; with friction angle phi = atan(mu_cone) on the cone surface, wedge statics give a radial (clamping) force
 *    P / tan(beta + phi) on the plug, which adds a friction capacity mu_wire x clamping force shared by the wires.
 *    The wedge is self-releasing on unloading when beta > phi.
 *  - Anchorage efficiency = capacity / (A sigma_u); the wire, not the anchorage, should govern (efficiency >= 1).
 * The stress concentration at the end-fitting entrance is applied as `scf` in {@link minerDamage}. Crimped /
 * swaged seal rings are not modelled: their check is contact-pressure retention of the seal against the bore
 * pressure with the supplier's qualified crimp parameters.
 * Reference: wedge statics (e.g. Meriam & Kraige, Statics, Ch. 6); bond anchorage as for grouted tendons. Bond
 * stress and friction coefficients must come from anchorage tests.
 * @param {object} o
 * @param {number} o.bondStress tau_b, Pa
 * @param {number} o.perimeter wire perimeter, m
 * @param {number} o.embeddedLength m
 * @param {number} o.wireArea m^2
 * @param {number} o.uts wire UTS, Pa
 * @param {number} [o.nWires=1] wires sharing the wedge clamping force
 * @param {number} [o.coneHalfAngle] beta, rad (omit to ignore the wedge action)
 * @param {number} [o.muCone=0] friction coefficient epoxy/body
 * @param {number} [o.muWire=0] friction coefficient wire/epoxy
 * @returns {{bondCapacity:number, wedgeAmplification:number, wedgeCapacityFraction:number, capacity:number,
 *   wireStrength:number, efficiency:number, selfReleasing:boolean|null}} N per wire; clamping force per unit pull;
 *   friction capacity from wedging as a fraction of the pull; N; N; -
 */
export function endFittingAnchorage({ bondStress, perimeter, embeddedLength, wireArea, uts, nWires = 1, coneHalfAngle, muCone = 0, muWire = 0 }) {
  nonneg(bondStress, 'bondStress'); pos(perimeter, 'perimeter'); pos(embeddedLength, 'embeddedLength'); pos(wireArea, 'wireArea');
  pos(uts, 'uts'); pos(nWires, 'nWires'); nonneg(muCone, 'muCone'); nonneg(muWire, 'muWire');
  const bondCapacity = bondStress * perimeter * embeddedLength;
  let amp = 0, selfReleasing = null;
  if (coneHalfAngle !== undefined) {
    if (!(coneHalfAngle > 0 && coneHalfAngle < Math.PI / 2)) throw new RangeError('coneHalfAngle must be in (0, pi/2)');
    const phi = Math.atan(muCone);
    if (coneHalfAngle + phi >= Math.PI / 2) throw new RangeError('coneHalfAngle + friction angle must be below pi/2');
    amp = 1 / Math.tan(coneHalfAngle + phi);
    selfReleasing = coneHalfAngle > phi;
  }
  // friction from wedging grows with the pull itself: P_cap = P_bond + muWire amp P_cap  ->  P_bond/(1 - muWire amp)
  const frac = muWire * amp;
  const capacity = frac >= 1 ? Infinity : bondCapacity / (1 - frac);
  const wireStrength = wireArea * uts;
  return { bondCapacity, wedgeAmplification: amp, wedgeCapacityFraction: frac, capacity, wireStrength, efficiency: capacity / wireStrength, selfReleasing };
}

function devAndHydro(s, name) {
  if (!Array.isArray(s) || s.length !== 6 || !s.every(Number.isFinite)) throw new RangeError(`${name} must be [xx, yy, zz, xy, yz, zx] in Pa`);
  const h = (s[0] + s[1] + s[2]) / 3;
  const j2 = ((s[0] - s[1]) ** 2 + (s[1] - s[2]) ** 2 + (s[2] - s[0]) ** 2) / 6 + s[3] ** 2 + s[4] ** 2 + s[5] ** 2;
  return { h, j2 };
}

/** Principal stresses of a symmetric tensor [xx, yy, zz, xy, yz, zx], descending (trigonometric solution). */
function principal(s) {
  const [a, b, c, d, e, f] = s;
  const m = (a + b + c) / 3;
  const p2 = (a - m) ** 2 + (b - m) ** 2 + (c - m) ** 2 + 2 * (d * d + e * e + f * f);
  if (p2 === 0) return [m, m, m];
  const p = Math.sqrt(p2 / 6);
  const B = [(a - m) / p, (b - m) / p, (c - m) / p, d / p, e / p, f / p];
  const detB = B[0] * (B[1] * B[2] - B[4] * B[4]) - B[3] * (B[3] * B[2] - B[4] * B[5]) + B[5] * (B[3] * B[4] - B[1] * B[5]);
  const phi = Math.acos(Math.min(1, Math.max(-1, detB / 2))) / 3;
  const e1 = m + 2 * p * Math.cos(phi), e3 = m + 2 * p * Math.cos(phi + 2 * Math.PI / 3);
  return [e1, 3 * m - e1 - e3, e3];
}

/**
 * Multiaxial high-cycle fatigue criteria for proportional loading, simple forms. Stress tensors are given as
 * [xx, yy, zz, xy, yz, zx] in Pa: `amplitude` (alternating part) and `mean`.
 *  - Sines (1955):      sqrt(J2,a) + k sigma_H,mean <= lambda
 *  - Crossland (1956):  sqrt(J2,a) + k sigma_H,max  <= lambda
 *  - Dang Van (1973), proportional-loading form: tau_a + a sigma_H,max <= b with tau_a the Tresca shear amplitude
 *  - Findley (1959) critical plane: max over planes of (tau_a + k sigma_n,max) <= f, searched on a grid of plane
 *    normals (sigma_n,max = normal stress from mean + |normal stress amplitude|).
 * Material constants are arguments. Usual calibrations from the fully reversed torsion limit t and bending limit
 * f_1: Dang Van a = 3 (t/f_1 - 1/2), b = t; Crossland k = 3 t/f_1 - sqrt(3), lambda = t; Sines lambda = t with k
 * from a repeated-bending limit.
 * Reference: Socie & Marquis, Multiaxial Fatigue (SAE, 2000), Ch. 5; Papadopoulos et al. (1997) Int. J. Fatigue
 * 19, 219-235.
 * Validity: proportional (in-phase) loading, infinite-life assessment.
 * @param {object} o
 * @param {'sines'|'crossland'|'dangvan'|'findley'} o.criterion
 * @param {number[]} o.amplitude stress amplitude tensor, Pa
 * @param {number[]} [o.mean] mean stress tensor, Pa (default zero)
 * @param {number} o.k hydrostatic / normal stress sensitivity (a for Dang Van), -
 * @param {number} o.limit right-hand side (lambda, b or f), Pa
 * @param {number} [o.nGrid=90] angular divisions for the Findley plane search
 * @returns {{equivalent:number, utilisation:number, shearTerm:number, normalTerm:number}} Pa, -, Pa, Pa
 */
export function multiaxialCriterion({ criterion, amplitude, mean = [0, 0, 0, 0, 0, 0], k, limit, nGrid = 90 }) {
  fin(k, 'k'); pos(limit, 'limit');
  const A = devAndHydro(amplitude, 'amplitude'), M = devAndHydro(mean, 'mean');
  let shearTerm, normalTerm;
  if (criterion === 'sines') { shearTerm = Math.sqrt(A.j2); normalTerm = M.h; }
  else if (criterion === 'crossland') { shearTerm = Math.sqrt(A.j2); normalTerm = M.h + Math.abs(A.h); }
  else if (criterion === 'dangvan') {
    const pr = principal(amplitude);
    shearTerm = (pr[0] - pr[2]) / 2; normalTerm = M.h + Math.abs(A.h);
  } else if (criterion === 'findley') {
    if (!(Number.isInteger(nGrid) && nGrid >= 8)) throw new RangeError('nGrid must be an integer >= 8');
    let best = -Infinity; shearTerm = 0; normalTerm = 0;
    const trac = (s, n) => [s[0] * n[0] + s[3] * n[1] + s[5] * n[2], s[3] * n[0] + s[1] * n[1] + s[4] * n[2], s[5] * n[0] + s[4] * n[1] + s[2] * n[2]];
    for (let i = 0; i <= nGrid; i++) {
      const th = Math.PI / 2 * i / nGrid; // polar angle, upper hemisphere
      const nPhi = i === 0 ? 1 : 4 * nGrid;
      for (let j = 0; j < nPhi; j++) {
        const ph = 2 * Math.PI * j / nPhi;
        const n = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
        const ta = trac(amplitude, n), tm = trac(mean, n);
        const sna = ta[0] * n[0] + ta[1] * n[1] + ta[2] * n[2];
        const snm = tm[0] * n[0] + tm[1] * n[1] + tm[2] * n[2];
        const tau = Math.sqrt(Math.max(ta[0] ** 2 + ta[1] ** 2 + ta[2] ** 2 - sna * sna, 0));
        const v = tau + k * (snm + Math.abs(sna));
        if (v > best) { best = v; shearTerm = tau; normalTerm = snm + Math.abs(sna); }
      }
    }
  } else throw new RangeError("criterion must be 'sines', 'crossland', 'dangvan' or 'findley'");
  const equivalent = shearTerm + k * normalTerm;
  return { equivalent, utilisation: equivalent / limit, shearTerm, normalTerm };
}
