// S–N and T–N fatigue resistance curves and their evaluation.
//
// Stress ranges are passed in Pa (SI) everywhere in the public API. The published curve constants are defined for
// stress ranges in MPa; the conversion happens inside this file (unit trap: log a is tied to MPa).

import { gamma, gammaP, gammaQ } from '../math/special.js';

const MPA = 1e6;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function positive(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

// DNV-RP-C203 Table 2-1 (air): [m1, log a1 (N <= 1e7), log a2 (m2 = 5, N > 1e7), fatigue limit at 1e7 [MPa],
// thickness exponent k, structural stress concentration embedded in the detail relative to curve D].
const AIR = {
  B1: [4.0, 15.117, 17.146, 106.97, 0.0, null],
  B2: [4.0, 14.885, 16.856, 93.59, 0.0, null],
  C: [3.0, 12.592, 16.320, 73.10, 0.05, null],
  C1: [3.0, 12.449, 16.081, 65.50, 0.10, null],
  C2: [3.0, 12.301, 15.835, 58.48, 0.15, null],
  D: [3.0, 12.164, 15.606, 52.63, 0.20, 1.00],
  E: [3.0, 12.010, 15.350, 46.78, 0.20, 1.13],
  F: [3.0, 11.855, 15.091, 41.52, 0.25, 1.27],
  F1: [3.0, 11.699, 14.832, 36.84, 0.25, 1.43],
  F3: [3.0, 11.546, 14.576, 32.75, 0.25, 1.61],
  G: [3.0, 11.398, 14.330, 29.24, 0.25, 1.80],
  W1: [3.0, 11.261, 14.101, 26.32, 0.25, 2.00],
  W2: [3.0, 11.107, 13.845, 23.39, 0.25, 2.25],
  W3: [3.0, 10.970, 13.617, 21.05, 0.25, 2.50],
  // Tubular joints, 2016 and later editions: reference thickness 16 mm.
  T: [3.0, 12.48, 16.13, 67.09, 0.25, 1.00],
};

// DNV-RP-C203 Table 2-2 (seawater with cathodic protection): log a1 for the upper branch (N <= 1e6, 1.8e6 for T);
// the lower branch (m2 = 5) is the same as in air.
const CP_LOGA1 = {
  B1: 14.917, B2: 14.685, C: 12.192, C1: 12.049, C2: 11.901, D: 11.764, E: 11.610, F: 11.455, F1: 11.299,
  F3: 11.146, G: 10.998, W1: 10.861, W2: 10.707, W3: 10.570, T: 12.18,
};

// DNV-RP-C203 Table 2-4 (free corrosion): single slope m = 3 for every curve, no fatigue limit.
const FREE_LOGA = {
  B1: 12.436, B2: 12.262, C: 12.115, C1: 11.972, C2: 11.824, D: 11.687, E: 11.533, F: 11.378, F1: 11.222,
  F3: 11.068, G: 10.921, W1: 10.784, W2: 10.630, W3: 10.493, T: 12.03,
};

/** Names of the DNV-RP-C203 curves available from {@link snCurve}. */
export const SN_CURVE_NAMES = Object.freeze(Object.keys(AIR));

/** Environments accepted by {@link snCurve}. */
export const SN_ENVIRONMENTS = Object.freeze(['air', 'seawater-cp', 'free-corrosion']);

/**
 * DNV-RP-C203 S–N curve record.
 * Reference: DNV-RP-C203 (2016–2021 edition family) Tables 2-1 (air), 2-2 (seawater with cathodic protection)
 * and 2-4 (free corrosion). log N = log a − m·log(Δσ·(t/tRef)^k), Δσ in MPa.
 * Validity: welded and base-material steel details with yield strength up to 960 MPa in air (550 MPa in seawater).
 * @param {string} name curve designation, one of SN_CURVE_NAMES
 * @param {string} [environment='air'] 'air' | 'seawater-cp' | 'free-corrosion'
 * @returns {{name:string, environment:string, m1:number, logA1:number, m2:(number|null), logA2:(number|null),
 *   kneeStress:number, kneeCycles:number, fatigueLimitTabulated:(number|null), k:number, tRef:number,
 *   scfRelativeToD:(number|null)}} kneeStress and fatigueLimitTabulated in Pa, tRef in m, kneeCycles in cycles;
 *   logA values refer to stress in MPa. kneeStress is the exact intersection of the two branches (the tabulated
 *   fatigue limit agrees with it to the rounding of the table).
 */
export function snCurve(name, environment = 'air') {
  const row = AIR[name];
  need(row !== undefined, `name: unknown S-N curve "${name}"`);
  need(SN_ENVIRONMENTS.includes(environment), `environment: unknown environment "${environment}"`);
  const [m1Air, logA1Air, logA2, limit, k, scf] = row;
  const tRef = name === 'T' ? 0.016 : 0.025;
  if (environment === 'free-corrosion') {
    return {
      name, environment, m1: 3.0, logA1: FREE_LOGA[name], m2: null, logA2: null, kneeStress: 0,
      kneeCycles: Infinity, fatigueLimitTabulated: null, k, tRef, scfRelativeToD: scf,
    };
  }
  const logA1 = environment === 'air' ? logA1Air : CP_LOGA1[name];
  const m2 = 5.0;
  const logKnee = (logA2 - logA1) / (m2 - m1Air);
  return {
    name, environment, m1: m1Air, logA1, m2, logA2, kneeStress: 10 ** logKnee * MPA,
    kneeCycles: 10 ** (logA1 - m1Air * logKnee), fatigueLimitTabulated: environment === 'air' ? limit * MPA : null,
    k, tRef, scfRelativeToD: scf,
  };
}

/**
 * Build a bilinear S–N curve from a detail category (stress range at 2·10⁶ cycles), the relationship used by
 * DNV-RP-C203 (the curve letter corresponds to the stress range at 2e6 cycles), EN 1993-1-9 and the IIW FAT classes.
 * Derived, not tabulated: a1 = 2e6·Δσc^m1; the second branch starts at kneeCycles with slope m2 and is continuous.
 * @param {object} o
 * @param {number} o.category stress range at 2e6 cycles [Pa]
 * @param {number} [o.m1=3] upper-branch inverse slope [-]
 * @param {number|null} [o.m2=5] lower-branch inverse slope [-]; null for a single-slope curve
 * @param {number} [o.kneeCycles=1e7] cycles at the slope change [cycles]
 * @param {number} [o.k=0] thickness exponent [-]
 * @param {number} [o.tRef=0.025] reference thickness [m]
 * @param {string} [o.name='custom']
 * @returns {object} curve record with the same fields as {@link snCurve}
 */
export function snCurveFromDetailCategory({ category, m1 = 3, m2 = 5, kneeCycles = 1e7, k = 0, tRef = 0.025, name = 'custom' }) {
  positive(category, 'category');
  positive(m1, 'm1');
  positive(tRef, 'tRef');
  const logA1 = Math.log10(2e6) + m1 * Math.log10(category / MPA);
  if (m2 === null || m2 === undefined) {
    return { name, environment: 'custom', m1, logA1, m2: null, logA2: null, kneeStress: 0, kneeCycles: Infinity,
      fatigueLimitTabulated: null, k, tRef, scfRelativeToD: null };
  }
  need(m2 > m1, 'm2 must exceed m1');
  positive(kneeCycles, 'kneeCycles');
  const logKnee = (logA1 - Math.log10(kneeCycles)) / m1;
  const logA2 = Math.log10(kneeCycles) + m2 * logKnee;
  return { name, environment: 'custom', m1, logA1, m2, logA2, kneeStress: 10 ** logKnee * MPA, kneeCycles,
    fatigueLimitTabulated: null, k, tRef, scfRelativeToD: null };
}

/**
 * Thickness correction factor on stress, (t/tRef)^k for t > tRef, 1 otherwise.
 * Reference: DNV-RP-C203 §2.4.3.
 * @param {number} thickness plate (or brace wall) thickness [m]
 * @param {number} tRef reference thickness [m]
 * @param {number} k thickness exponent [-]
 * @returns {number} factor on stress range [-], >= 1
 */
export function thicknessFactor(thickness, tRef, k) {
  positive(thickness, 'thickness');
  positive(tRef, 'tRef');
  need(k >= 0, 'k must be non-negative');
  return thickness > tRef ? (thickness / tRef) ** k : 1;
}

/**
 * Mean-stress reduction factor on the stress range for details without significant residual stress.
 * fm = (σt + c·|σc|)/(σt + |σc|) with c = 0.6 for base material and c = 0.8 for welded connections after
 * post-weld heat treatment or with documented low residual stress. fm = 1 when the cycle is fully tensile.
 * Reference: DNV-RP-C203 §2.5. Not applicable to as-welded details (use 1.0).
 * @param {object} o
 * @param {number} o.stressMax maximum stress of the cycle, tension positive [Pa]
 * @param {number} o.stressMin minimum stress of the cycle [Pa]
 * @param {string} [o.detail='base'] 'base' (c = 0.6) or 'welded' (c = 0.8)
 * @returns {number} reduction factor on stress range [-], between c and 1
 */
export function meanStressFactor({ stressMax, stressMin, detail = 'base' }) {
  need(Number.isFinite(stressMax) && Number.isFinite(stressMin), 'stressMax and stressMin must be finite');
  need(stressMax >= stressMin, 'stressMax must not be less than stressMin');
  need(detail === 'base' || detail === 'welded', 'detail must be "base" or "welded"');
  const c = detail === 'base' ? 0.6 : 0.8;
  const tens = Math.max(stressMax, 0) - Math.max(stressMin, 0);
  const comp = Math.max(-stressMin, 0) - Math.max(-stressMax, 0);
  if (tens + comp === 0) return 1;
  return (tens + c * comp) / (tens + comp);
}

/**
 * Stress ratio R = σmin/σmax of a cycle.
 * @param {number} stressMax [Pa]
 * @param {number} stressMin [Pa]
 * @returns {number} R [-]
 */
export function stressRatio(stressMax, stressMin) {
  need(stressMax !== 0, 'stressMax must be non-zero');
  return stressMin / stressMax;
}

function stressModifier(curve, { thickness, scf = 1, meanStressFactor: fm = 1, k } = {}) {
  positive(scf, 'scf');
  need(fm > 0 && fm <= 1, 'meanStressFactor must be in (0, 1]');
  const kk = k === undefined ? curve.k : k;
  const tf = thickness === undefined ? 1 : thicknessFactor(thickness, curve.tRef, kk);
  return scf * fm * tf;
}

/**
 * S–N constants expressed for the nominal stress range in Pa, with stress concentration, mean-stress factor and
 * thickness correction folded in: N = a1·S^(−m1) for S >= sKnee, N = a2·S^(−m2) below.
 * This is the neutral form consumed by the damage and spectral engines.
 * @param {object} curve record from {@link snCurve} or {@link snCurveFromDetailCategory}
 * @param {object} [o]
 * @param {number} [o.thickness] thickness for the size correction [m]; omitted = no correction
 * @param {number} [o.scf=1] stress concentration factor applied to the nominal range [-]
 * @param {number} [o.meanStressFactor=1] reduction factor on the range, see {@link meanStressFactor} [-]
 * @param {number} [o.k] override of the thickness exponent [-]
 * @returns {{m1:number, a1:number, m2:(number|null), a2:(number|null), sKnee:number}} a in Pa^m·cycles, sKnee in Pa
 *   (0 for single-slope curves)
 */
export function snConstants(curve, o = {}) {
  const f = stressModifier(curve, o);
  const a1 = 10 ** curve.logA1 * (MPA / f) ** curve.m1;
  if (curve.m2 === null) return { m1: curve.m1, a1, m2: null, a2: null, sKnee: 0 };
  return { m1: curve.m1, a1, m2: curve.m2, a2: 10 ** curve.logA2 * (MPA / f) ** curve.m2, sKnee: curve.kneeStress / f };
}

/**
 * Cycles to failure from constants in the neutral form of {@link snConstants}.
 * @param {{m1:number,a1:number,m2:(number|null),a2:(number|null),sKnee:number}} sn
 * @param {number} stressRange [Pa] (same unit as the constants)
 * @returns {number} cycles to failure; Infinity for a zero range
 */
export function cyclesFromConstants(sn, stressRange) {
  need(stressRange >= 0, 'stressRange must be non-negative');
  if (stressRange === 0) return Infinity;
  if (sn.m2 !== null && sn.m2 !== undefined && stressRange < sn.sKnee) return sn.a2 * stressRange ** -sn.m2;
  return sn.a1 * stressRange ** -sn.m1;
}

/**
 * Characteristic (mean minus two standard deviations) cycles to failure N(Δσ).
 * Reference: DNV-RP-C203 Eq. (2.4.1) and (2.4.3).
 * @param {object} o
 * @param {object} o.curve record from {@link snCurve}
 * @param {number} o.stressRange nominal (or hot-spot) stress range [Pa]
 * @param {number} [o.thickness] thickness [m]; omitted = no size correction
 * @param {number} [o.scf=1] stress concentration factor [-]
 * @param {number} [o.meanStressFactor=1] [-]
 * @param {number} [o.k] override of the thickness exponent [-]
 * @returns {number} cycles to failure [cycles]
 */
export function cyclesToFailure({ curve, stressRange, thickness, scf, meanStressFactor: fm, k }) {
  return cyclesFromConstants(snConstants(curve, { thickness, scf, meanStressFactor: fm, k }), stressRange);
}

/**
 * Stress range that gives a stated number of cycles to failure (inverse of {@link cyclesToFailure}).
 * @param {object} o
 * @param {object} o.curve
 * @param {number} o.cycles target cycles [cycles]
 * @param {number} [o.thickness] [m]
 * @param {number} [o.scf=1] [-]
 * @param {number} [o.meanStressFactor=1] [-]
 * @param {number} [o.k] [-]
 * @returns {number} nominal stress range [Pa]
 */
export function stressRangeAtCycles({ curve, cycles, thickness, scf, meanStressFactor: fm, k }) {
  positive(cycles, 'cycles');
  const sn = snConstants(curve, { thickness, scf, meanStressFactor: fm, k });
  const s1 = (sn.a1 / cycles) ** (1 / sn.m1);
  if (sn.m2 === null || s1 >= sn.sKnee) return s1;
  return (sn.a2 / cycles) ** (1 / sn.m2);
}

/**
 * Closed-form Miner damage for a two-parameter Weibull long-term stress-range distribution,
 * Q(Δσ) = exp(−(Δσ/q)^h), against a single- or two-slope S–N curve:
 *   D = n0·[ q^m1/a1·Γ(1 + m1/h; (S1/q)^h) + q^m2/a2·γ(1 + m2/h; (S1/q)^h) ]
 * with Γ(·;·) the upper and γ(·;·) the lower incomplete gamma function and S1 the stress at the slope change.
 * Reference: DNV-RP-C203 §5.1 (Eq. 5.1.1 form).
 * @param {object} o
 * @param {{m1:number,a1:number,m2:(number|null),a2:(number|null),sKnee:number}} o.sn constants, see {@link snConstants}
 * @param {number} o.scale Weibull scale parameter q [Pa]
 * @param {number} o.shape Weibull shape parameter h [-], typically 0.5–1.5
 * @param {number} o.cycles total number of cycles n0 [cycles]
 * @returns {number} Miner damage [-]
 */
export function weibullClosedFormDamage({ sn, scale, shape, cycles }) {
  positive(scale, 'scale');
  positive(shape, 'shape');
  need(cycles >= 0, 'cycles must be non-negative');
  const g1 = 1 + sn.m1 / shape;
  if (sn.m2 === null || sn.m2 === undefined || !(sn.sKnee > 0)) {
    return (cycles * scale ** sn.m1 * gamma(g1)) / sn.a1;
  }
  const x = (sn.sKnee / scale) ** shape;
  const g2 = 1 + sn.m2 / shape;
  const upper = gamma(g1) * gammaQ(g1, x);
  const lower = gamma(g2) * gammaP(g2, x);
  return cycles * ((scale ** sn.m1 * upper) / sn.a1 + (scale ** sn.m2 * lower) / sn.a2);
}

/**
 * Weibull scale parameter from the stress range exceeded once in n0 cycles: q = Δσ0/(ln n0)^(1/h).
 * Reference: DNV-RP-C203 Eq. (5.1.2).
 * @param {number} maxStressRange Δσ0 [Pa]
 * @param {number} shape h [-]
 * @param {number} cycles n0 [cycles], > 1
 * @returns {number} scale q [Pa]
 */
export function weibullScaleFromMax(maxStressRange, shape, cycles) {
  positive(maxStressRange, 'maxStressRange');
  positive(shape, 'shape');
  need(cycles > 1, 'cycles must exceed 1');
  return maxStressRange / Math.log(cycles) ** (1 / shape);
}

/**
 * Allowable maximum stress range Δσ0 (exceeded once in the reference number of cycles) for a target fatigue life
 * under a Weibull long-term distribution, so that damage × DFF = allowable Miner sum.
 * Solved from the closed form of {@link weibullClosedFormDamage}; reproduces the tables of DNV-RP-C203 §5.
 * @param {object} o
 * @param {object} o.curve S–N curve record
 * @param {number} o.shape Weibull shape parameter h [-]
 * @param {number} o.cycles number of cycles in the design life n0 [cycles]
 * @param {number} [o.dff=1] design fatigue factor [-]
 * @param {number} [o.allowableDamage=1] Miner sum at failure [-]
 * @param {number} [o.thickness] [m]
 * @param {number} [o.scf=1] [-]
 * @param {number} [o.meanStressFactor=1] [-]
 * @returns {{maxStressRange:number, scale:number, damage:number}} Δσ0 [Pa], Weibull scale q [Pa], damage at the
 *   solution without the DFF [-]
 */
export function allowableStressRangeWeibull({ curve, shape, cycles, dff = 1, allowableDamage = 1, thickness, scf, meanStressFactor: fm }) {
  positive(shape, 'shape');
  need(cycles > 1, 'cycles must exceed 1');
  positive(dff, 'dff');
  positive(allowableDamage, 'allowableDamage');
  const sn = snConstants(curve, { thickness, scf, meanStressFactor: fm });
  const target = allowableDamage / dff;
  const f = (lnq) => weibullClosedFormDamage({ sn, scale: Math.exp(lnq), shape, cycles }) - target;
  let lo = Math.log(1e-3);
  let hi = Math.log(1e12);
  need(f(lo) < 0 && f(hi) > 0, 'cycles: no allowable stress range exists in the search interval');
  for (let i = 0; i < 200 && hi - lo > 1e-13; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) hi = mid; else lo = mid;
  }
  const scale = Math.exp(0.5 * (lo + hi));
  return { maxStressRange: scale * Math.log(cycles) ** (1 / shape), scale, damage: target };
}

/**
 * Fatigue utilisation including the design fatigue factor: η = D·DFF (acceptable when <= 1).
 * Reference: DNV-RP-C203 (the DFF is chosen by the caller from the governing design standard).
 * @param {number} damage Miner sum over the design life [-]
 * @param {number} dff design fatigue factor [-], typically 1–10
 * @returns {number} utilisation [-]
 */
export function dffUtilisation(damage, dff) {
  need(damage >= 0, 'damage must be non-negative');
  positive(dff, 'dff');
  return damage * dff;
}

// ---------------------------------------------------------------------------------------------------------------
// Mooring line T–N / S–N curves
// ---------------------------------------------------------------------------------------------------------------

/**
 * API RP 2SK T–N parameters, N·R^M = K with R = tension range / reference breaking strength.
 * Reference: API RP 2SK (3rd ed.) T–N curves. Lm is the ratio of mean load to reference breaking
 * strength (wire rope only). For chain the reference breaking strength is that of ORQ common chain link of the
 * same diameter, whatever the actual grade.
 * @param {string} type 'chain-stud' | 'chain-studless' | 'connecting-link' | 'rope-six-strand' | 'rope-spiral-strand'
 * @param {number} [meanLoadRatio=0.3] Lm [-], used by the two wire-rope curves
 * @returns {{M:number, K:number}} slope M [-] and intercept K [cycles]
 */
export function mooringTnApi(type, meanLoadRatio = 0.3) {
  need(meanLoadRatio >= 0 && meanLoadRatio < 1, 'meanLoadRatio must be in [0, 1)');
  switch (type) {
    case 'chain-stud': return { M: 3.0, K: 1000 };
    case 'chain-studless': return { M: 3.0, K: 316 };
    case 'connecting-link': return { M: 3.0, K: 178 };
    case 'rope-six-strand': return { M: 4.09, K: 10 ** (3.20 - 2.79 * meanLoadRatio) };
    case 'rope-spiral-strand': return { M: 5.05, K: 10 ** (3.25 - 3.43 * meanLoadRatio) };
    default: throw new RangeError(`type: unknown API mooring component "${type}"`);
  }
}

/**
 * Cycles to failure of a mooring component by the API RP 2SK tension-range convention.
 * @param {object} o
 * @param {string} o.type see {@link mooringTnApi}
 * @param {number} o.tensionRange double-amplitude tension range [N]
 * @param {number} o.referenceBreakingStrength reference breaking strength [N] (ORQ value for chain; after any
 *   corrosion/wear reduction chosen by the caller)
 * @param {number} [o.meanLoadRatio=0.3] mean load / reference breaking strength [-]
 * @returns {number} cycles to failure [cycles]
 */
export function mooringCyclesApi({ type, tensionRange, referenceBreakingStrength, meanLoadRatio = 0.3 }) {
  need(tensionRange >= 0, 'tensionRange must be non-negative');
  positive(referenceBreakingStrength, 'referenceBreakingStrength');
  if (tensionRange === 0) return Infinity;
  const { M, K } = mooringTnApi(type, meanLoadRatio);
  return K * (tensionRange / referenceBreakingStrength) ** -M;
}

/**
 * DNV-OS-E301 S–N parameters for mooring components, nc(s) = aD·s^(−m) with s the nominal stress range in MPa.
 * Reference: DNV-OS-E301 Ch.2 Sec.2 [6.2], Table "S-N fatigue curve parameters".
 * @param {string} type 'chain-stud' | 'chain-studless' | 'rope-stranded' | 'rope-spiral'
 * @returns {{aD:number, m:number}} intercept aD [cycles·MPa^m] and slope m [-]
 */
export function mooringSnDnv(type) {
  switch (type) {
    case 'chain-stud': return { aD: 1.2e11, m: 3.0 };
    case 'chain-studless': return { aD: 6.0e10, m: 3.0 };
    case 'rope-stranded': return { aD: 3.4e14, m: 4.0 };
    case 'rope-spiral': return { aD: 1.7e17, m: 4.8 };
    default: throw new RangeError(`type: unknown DNV mooring component "${type}"`);
  }
}

/**
 * Nominal cross-section used with the DNV-OS-E301 mooring S–N curves: 2·πd²/4 for chain (both legs of the link),
 * πd²/4 for wire rope, with the diameter reduced by a corrosion/wear allowance.
 * Reference: DNV-OS-E301 Ch.2 Sec.2 [6.2] (the standard takes half of the corrosion allowance for fatigue; the
 * caller decides how much to pass).
 * @param {object} o
 * @param {string} o.type see {@link mooringSnDnv}
 * @param {number} o.diameter nominal chain bar or rope diameter [m]
 * @param {number} [o.corrosionAllowance=0] total reduction of diameter [m]
 * @returns {{area:number, diameter:number}} nominal area [m²] and reduced diameter [m]
 */
export function mooringNominalArea({ type, diameter, corrosionAllowance = 0 }) {
  positive(diameter, 'diameter');
  need(corrosionAllowance >= 0 && corrosionAllowance < diameter, 'corrosionAllowance must be in [0, diameter)');
  mooringSnDnv(type);
  const d = diameter - corrosionAllowance;
  const single = (Math.PI * d * d) / 4;
  return { area: type.startsWith('chain') ? 2 * single : single, diameter: d };
}

/**
 * Cycles to failure of a mooring component by the DNV-OS-E301 nominal-stress convention.
 * @param {object} o
 * @param {string} o.type see {@link mooringSnDnv}
 * @param {number} o.tensionRange double-amplitude tension range [N]
 * @param {number} o.diameter nominal diameter [m]
 * @param {number} [o.corrosionAllowance=0] reduction of diameter [m]
 * @returns {{cycles:number, stressRange:number}} cycles to failure [cycles] and nominal stress range [Pa]
 */
export function mooringCyclesDnv({ type, tensionRange, diameter, corrosionAllowance = 0 }) {
  need(tensionRange >= 0, 'tensionRange must be non-negative');
  const { aD, m } = mooringSnDnv(type);
  const { area } = mooringNominalArea({ type, diameter, corrosionAllowance });
  const s = tensionRange / area;
  return { cycles: s === 0 ? Infinity : aD * (s / MPA) ** -m, stressRange: s };
}

/**
 * Neutral S–N constants (stress in Pa) for a DNV-OS-E301 mooring curve, for use with the damage/spectral engines.
 * @param {string} type see {@link mooringSnDnv}
 * @returns {{m1:number, a1:number, m2:null, a2:null, sKnee:number}}
 */
export function mooringSnConstants(type) {
  const { aD, m } = mooringSnDnv(type);
  return { m1: m, a1: aD * MPA ** m, m2: null, a2: null, sKnee: 0 };
}

/**
 * Catalogue break load of offshore mooring chain, BL = c·d²·(44 − 0.08·d) with d in mm and BL in kN.
 * Reference: IACS UR W22 / DNV-OS-E302 break-load formula; grade coefficients c: ORQ 0.0211, R3 0.0223,
 * R3S 0.0249, R4 0.0274, R4S 0.0304, R5 0.0320. Pass `coefficient` to override.
 * Validity: stud-link and studless common link, d ≈ 50–200 mm.
 * @param {object} o
 * @param {number} o.diameter nominal bar diameter [m]
 * @param {string} [o.grade='R3'] 'ORQ' | 'R3' | 'R3S' | 'R4' | 'R4S' | 'R5'
 * @param {number} [o.coefficient] override of c [kN/mm² scale as in the formula]
 * @returns {number} minimum breaking load [N]
 */
export function chainBreakingLoad({ diameter, grade = 'R3', coefficient }) {
  positive(diameter, 'diameter');
  const table = { ORQ: 0.0211, R3: 0.0223, R3S: 0.0249, R4: 0.0274, R4S: 0.0304, R5: 0.0320 };
  const c = coefficient === undefined ? table[grade] : coefficient;
  need(c !== undefined && c > 0, `grade: unknown chain grade "${grade}"`);
  const d = diameter * 1000;
  need(d < 550, 'diameter is outside the range of the break-load formula');
  return c * d * d * (44 - 0.08 * d) * 1000;
}
