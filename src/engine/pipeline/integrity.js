/**
 * Pipeline integrity assessment: corrosion metal-loss defects (ASME B31G, Modified B31G, effective area /
 * RSTRENG, DNV-RP-F101), remaining life and re-inspection, dents (ASME B31.8 Appendix R strains and depth
 * limits), wall-loss extreme-value extrapolation, a hook for crack assessment and a fitness-for-service
 * summary.
 *
 * Notes on scope:
 *  - Gouges and other mechanical metal loss are sharper than corrosion and may hide cracking or cold work;
 *    the corrosion equations apply only after dressing to a smooth contour and crack inspection.
 *  - Pressure reversal: a defect that survived a hydrotest can fail below the test pressure on
 *    re-pressurisation after ductile tearing during the test; remaining-strength results do not cover it.
 * All quantities SI.
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

function defectArgs({ D, t, d, L, smys }) {
  pos('D', D); pos('t', t); nonNeg('d', d); nonNeg('L', L); pos('smys', smys);
  if (d >= t) throw new RangeError('d must be smaller than t');
}

/**
 * ASME B31G (original) remaining strength of a corroded pipe.
 * z = L²/(D t), M = √(1 + 0.8 z), flow stress 1.1 SMYS;
 * z ≤ 20: S_F = S_flow (1 − ⅔ d/t)/(1 − ⅔ (d/t)/M) (parabolic area);  z > 20: S_F = S_flow (1 − d/t);
 * P_F = 2 S_F t/D;  P_S = P_F/SF.
 *
 * Reference: ASME B31G-2012 §2.2(a) (Level 1, original B31G).
 * Validity: 10 % ≤ d/t ≤ 80 % (shallower pits are acceptable without assessment; deeper must be repaired);
 * internal pressure loading only; blunt defects in ductile line pipe. The validity flag is returned.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.d maximum defect depth [m]
 * @param {number} o.L axial defect length [m]
 * @param {number} o.smys [Pa]
 * @param {number} [o.safetyFactor=1.39] SF on failure pressure [-] (not below 1.25; 1/0.72 = 1.39 is customary)
 * @returns {{z:number, M:number, flowStress:number, failureStress:number, failurePressure:number, safePressure:number,
 *   depthRatio:number, withinDepthLimits:boolean}} stresses and pressures [Pa]
 */
export function b31gOriginal({ D, t, d, L, smys, safetyFactor = 1.39 }) {
  defectArgs({ D, t, d, L, smys }); pos('safetyFactor', safetyFactor);
  const z = (L * L) / (D * t);
  const M = Math.sqrt(1 + 0.8 * z);
  const flowStress = 1.1 * smys;
  const r = d / t;
  const failureStress = z <= 20 ? (flowStress * (1 - (2 / 3) * r)) / (1 - ((2 / 3) * r) / M) : flowStress * (1 - r);
  const failurePressure = (2 * failureStress * t) / D;
  return { z, M, flowStress, failureStress, failurePressure, safePressure: failurePressure / safetyFactor, depthRatio: r,
    withinDepthLimits: r <= 0.8 };
}

/** Three-term Folias factor of Modified B31G. */
function foliasModified(z) {
  return z <= 50 ? Math.sqrt(1 + 0.6275 * z - 0.003375 * z * z) : 0.032 * z + 3.3;
}

/**
 * Modified B31G (0.85 dL) remaining strength.
 * Flow stress SMYS + 69 MPa; M = √(1 + 0.6275 z − 0.003375 z²) for z ≤ 50, 0.032 z + 3.3 for z > 50;
 * S_F = S_flow (1 − 0.85 d/t)/(1 − 0.85 (d/t)/M); P_F = 2 S_F t/D.
 *
 * Reference: ASME B31G-2012 §2.2(b); Kiefner & Vieth (1989), AGA PR 3-805.
 * Validity: d/t ≤ 80 %; SMYS up to 483 MPa (X70) for the SMYS + 69 MPa flow stress.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.d [m]
 * @param {number} o.L [m]
 * @param {number} o.smys [Pa]
 * @param {number} [o.safetyFactor=1.39] [-]
 * @returns {{z:number, M:number, flowStress:number, failureStress:number, failurePressure:number, safePressure:number,
 *   depthRatio:number, withinDepthLimits:boolean}}
 */
export function b31gModified({ D, t, d, L, smys, safetyFactor = 1.39 }) {
  defectArgs({ D, t, d, L, smys }); pos('safetyFactor', safetyFactor);
  const z = (L * L) / (D * t);
  const M = foliasModified(z);
  const flowStress = smys + 69e6;
  const r = d / t;
  const failureStress = (flowStress * (1 - 0.85 * r)) / (1 - (0.85 * r) / M);
  const failurePressure = (2 * failureStress * t) / D;
  return { z, M, flowStress, failureStress, failurePressure, safePressure: failurePressure / safetyFactor, depthRatio: r,
    withinDepthLimits: r <= 0.8 };
}

/**
 * Effective-area method (RSTRENG) for a measured river-bottom profile: every contiguous sub-length of the
 * profile is assessed with its own metal-loss area and Folias factor, and the lowest failure stress governs.
 * S_F = S_flow (1 − A/A_0)/(1 − (A/A_0)/M(L_ij)), A = trapezoidal area of the depth profile over L_ij,
 * A_0 = L_ij t, flow stress SMYS + 69 MPa, M as Modified B31G.
 *
 * Reference: ASME B31G-2012 §2.3 (Level 2, effective area); Kiefner & Vieth (1989).
 * Validity: as Modified B31G; profile points ordered along the pipe axis.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {ArrayLike<number>} o.x axial positions of the depth readings [m], ascending
 * @param {ArrayLike<number>} o.depth defect depth at each position [m]
 * @param {number} o.smys [Pa]
 * @param {number} [o.safetyFactor=1.39] [-]
 * @returns {{failureStress:number, failurePressure:number, safePressure:number, effectiveLength:number,
 *   effectiveArea:number, iStart:number, iEnd:number, maxDepth:number}} governing sub-length [m], area [m²]
 */
export function effectiveArea({ D, t, x, depth, smys, safetyFactor = 1.39 }) {
  pos('D', D); pos('t', t); pos('smys', smys); pos('safetyFactor', safetyFactor);
  const n = x.length;
  if (n < 2 || depth.length !== n) throw new RangeError('x and depth must have the same length (at least 2)');
  let maxDepth = 0;
  for (let i = 0; i < n; i++) {
    nonNeg('depth', depth[i]);
    if (depth[i] >= t) throw new RangeError('depth must be smaller than t everywhere');
    if (i > 0 && !(x[i] > x[i - 1])) throw new RangeError('x must be strictly ascending');
    if (depth[i] > maxDepth) maxDepth = depth[i];
  }
  const flow = smys + 69e6;
  let best = { failureStress: flow, effectiveLength: 0, effectiveArea: 0, iStart: 0, iEnd: 0 };
  for (let i = 0; i < n - 1; i++) {
    let A = 0;
    for (let j = i + 1; j < n; j++) {
      A += 0.5 * (depth[j] + depth[j - 1]) * (x[j] - x[j - 1]);
      const Lij = x[j] - x[i];
      const ratio = A / (Lij * t);
      const M = foliasModified((Lij * Lij) / (D * t));
      const sf = (flow * (1 - ratio)) / (1 - ratio / M);
      if (sf < best.failureStress) best = { failureStress: sf, effectiveLength: Lij, effectiveArea: A, iStart: i, iEnd: j };
    }
  }
  const failurePressure = (2 * best.failureStress * t) / D;
  return { ...best, failurePressure, safePressure: failurePressure / safetyFactor, maxDepth };
}

/**
 * DNV-RP-F101 single longitudinal corrosion defect under internal pressure.
 * Capacity: P_f = 2 t σ_u/(D − t) · (1 − d/t)/(1 − d/(t Q)),  Q = √(1 + 0.31 (L/√(D t))²).
 * Part B (allowable stress design): P_sw = F · P_f, F = F1·F2 with modelling factor F1 = 0.9 and the
 *   operational usage (design) factor F2.
 * Part A (calibrated partial factors), when gammaM, gammaD and epsD are supplied:
 *   p_corr = γ_m 2 t f_u/(D − t) · (1 − γ_d (d/t)*)/(1 − γ_d (d/t)* /Q),  (d/t)* = (d/t)_meas + ε_d StD[d/t];
 *   if γ_d (d/t)* ≥ 1 the allowable pressure is zero. The Part A factors depend on safety class, inspection
 *   method and sizing accuracy (RP Tables 3-2 to 3-8); they are not tabulated here and must be supplied.
 *
 * Reference: DNV-RP-F101 (2019) §3.7.3.1 (Part A, Eq. 3.2) and §8.2 (Part B).
 * Validity: d/t ≤ 0.85; ductile line pipe up to X80; no significant axial compression or bending
 * (those need the combined-loading equations of the RP, not included).
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.d [m]
 * @param {number} o.L [m]
 * @param {number} o.fu tensile strength to be used in design (SMTS, de-rated for temperature) [Pa]
 * @param {number} [o.F2=0.72] usage factor [-] (Part B)
 * @param {number} [o.F1=0.9] modelling factor [-] (Part B)
 * @param {number} [o.gammaM] partial factor for the prediction model and safety class [-] (Part A)
 * @param {number} [o.gammaD] partial factor for corrosion depth [-] (Part A)
 * @param {number} [o.epsD] fractile factor for corrosion depth [-] (Part A)
 * @param {number} [o.stdDepth=0] standard deviation of the measured d/t [-] (Part A)
 * @returns {{Q:number, failurePressure:number, safePressure:number, allowablePartA:(number|null), depthRatio:number,
 *   withinDepthLimits:boolean}} pressures [Pa] (Part B safe working pressure in `safePressure`)
 */
export function dnvF101Single({ D, t, d, L, fu, F2 = 0.72, F1 = 0.9, gammaM, gammaD, epsD, stdDepth = 0 }) {
  pos('D', D); pos('t', t); nonNeg('d', d); nonNeg('L', L); pos('fu', fu); pos('F1', F1); pos('F2', F2);
  if (d >= t) throw new RangeError('d must be smaller than t');
  const Q = Math.sqrt(1 + 0.31 * (L / Math.sqrt(D * t)) ** 2);
  const r = d / t;
  const failurePressure = ((2 * t * fu) / (D - t)) * ((1 - r) / (1 - r / Q));
  let allowablePartA = null;
  if (gammaM !== undefined || gammaD !== undefined || epsD !== undefined) {
    pos('gammaM', gammaM); pos('gammaD', gammaD); nonNeg('epsD', epsD); nonNeg('stdDepth', stdDepth);
    const rs = gammaD * (r + epsD * stdDepth);
    allowablePartA = rs >= 1 ? 0 : gammaM * ((2 * t * fu) / (D - t)) * ((1 - rs) / (1 - rs / Q));
  }
  return { Q, failurePressure, safePressure: F1 * F2 * failurePressure, allowablePartA, depthRatio: r, withinDepthLimits: r <= 0.85 };
}

/**
 * Remaining life of a corrosion defect: time until the safe pressure of the growing defect falls to the
 * operating pressure (MAOP), or the depth reaches the limiting fraction of the wall, whichever is first;
 * with the recommended re-inspection interval and the de-rated pressure at the end of a chosen interval.
 *
 * Growth: linear rates, or caller functions d(time), L(time) (e.g. a power law or a rate fitted between two
 * inspections).
 *
 * Reference: assessment loop of ASME B31G-2012 §1.6 / DNV-RP-F101 §3.11 with any of the capacity models of
 * this file; re-inspection at a fraction of the remaining life (0.5 is customary, cf. ASME B31.8S §7).
 *
 * @param {object} o
 * @param {function(object):{safePressure:number}} o.assess capacity model, e.g. b31gModified or (a) => dnvF101Single({...a, fu})
 * @param {object} o.defect arguments for `assess` at time zero (D, t, d, L, smys, ...)
 * @param {number} o.maop operating pressure to sustain [Pa]
 * @param {number|function(number):number} o.depthRate depth growth rate [m/year], or function d(tYears) [m]
 * @param {number|function(number):number} [o.lengthRate=0] length growth rate [m/year], or function L(tYears) [m]
 * @param {number} [o.depthLimit=0.8] limiting d/t [-]
 * @param {number} [o.horizon=100] longest time considered [years]
 * @param {number} [o.reinspectFraction=0.5] fraction of remaining life for re-inspection [-]
 * @param {number} [o.interval] an inspection interval [years] for which the de-rated pressure is wanted
 * @returns {{remainingLife:number, governing:'pressure'|'depth'|'none', reinspectionInterval:number,
 *   safePressureNow:number, acceptableNow:boolean, deratedPressure:(number|null)}}
 *   times [years] (remainingLife = horizon when nothing governs); pressures [Pa]
 */
export function remainingLife({ assess, defect, maop, depthRate, lengthRate = 0, depthLimit = 0.8, horizon = 100, reinspectFraction = 0.5, interval }) {
  if (typeof assess !== 'function') throw new RangeError('assess must be a function');
  pos('maop', maop); pos('horizon', horizon); pos('depthLimit', depthLimit);
  const dOf = typeof depthRate === 'function' ? depthRate : (tt) => defect.d + nonNeg('depthRate', depthRate) * tt;
  const LOf = typeof lengthRate === 'function' ? lengthRate : (tt) => defect.L + nonNeg('lengthRate', lengthRate) * tt;
  const state = (tt) => {
    const d = dOf(tt);
    const ratio = d / defect.t;
    if (ratio >= depthLimit) return { fail: 'depth', p: 0 };
    const p = assess({ ...defect, d, L: LOf(tt) }).safePressure;
    return { fail: p < maop ? 'pressure' : null, p };
  };
  const now = state(0);
  let life = horizon;
  let governing = 'none';
  if (now.fail) { life = 0; governing = now.fail; } else {
    const N = 400;
    let lo = 0;
    let hi = null;
    for (let i = 1; i <= N; i++) {
      const tt = (horizon * i) / N;
      const s = state(tt);
      if (s.fail) { hi = tt; governing = s.fail; break; }
      lo = tt;
    }
    if (hi !== null) {
      for (let i = 0; i < 50; i++) {
        const mid = 0.5 * (lo + hi);
        const s = state(mid);
        if (s.fail) { hi = mid; governing = s.fail; } else lo = mid;
      }
      life = lo;
    }
  }
  return { remainingLife: life, governing, reinspectionInterval: life * reinspectFraction, safePressureNow: now.p, acceptableNow: !now.fail,
    deratedPressure: interval === undefined ? null : state(pos('interval', interval)).p };
}

/**
 * Strain in a plain dent from its measured curvatures.
 * ε1 = (t/2)(1/R0 − 1/R1)   circumferential bending (R1 negative for a re-entrant dent)
 * ε2 = −t/(2 R2)            longitudinal bending (R2 negative for a re-entrant profile)
 * ε3 = ½ (d/L)²             longitudinal extensional
 * ε_i = [ε1² − ε1(ε2 + ε3) + (ε2 + ε3)²]^½ (inside surface), ε_o = [ε1² + ε1(−ε2 + ε3) + (−ε2 + ε3)²]^½ (outside)
 * Acceptance: max(ε_i, ε_o) ≤ 6 %.
 *
 * Reference: ASME B31.8-2018 Nonmandatory Appendix R, R-2 and §851.4.1.
 *
 * @param {object} o
 * @param {number} o.t wall thickness [m]
 * @param {number} o.R0 initial pipe surface radius (D/2) [m]
 * @param {number} o.R1 radius of curvature in the transverse plane through the dent [m], negative if re-entrant
 * @param {number} o.R2 radius of curvature in the longitudinal plane through the dent [m], negative if re-entrant
 * @param {number} o.d dent depth [m]
 * @param {number} o.L dent length [m]
 * @param {number} [o.limit=0.06] strain limit [-]
 * @returns {{eps1:number, eps2:number, eps3:number, epsInside:number, epsOutside:number, epsMax:number, limit:number, uc:number}}
 */
export function dentStrain({ t, R0, R1, R2, d, L, limit = 0.06 }) {
  pos('t', t); pos('R0', R0); fin('R1', R1); fin('R2', R2); nonNeg('d', d); pos('L', L);
  if (R1 === 0 || R2 === 0) throw new RangeError('R1 and R2 must be non-zero');
  const eps1 = (t / 2) * (1 / R0 - 1 / R1);
  const eps2 = -t / (2 * R2);
  const eps3 = 0.5 * (d / L) ** 2;
  const a = eps2 + eps3;
  const b = -eps2 + eps3;
  const epsInside = Math.sqrt(eps1 * eps1 - eps1 * a + a * a);
  const epsOutside = Math.sqrt(eps1 * eps1 + eps1 * b + b * b);
  const epsMax = Math.max(epsInside, epsOutside);
  return { eps1, eps2, eps3, epsInside, epsOutside, epsMax, limit, uc: epsMax / limit };
}

/**
 * Dent depth screening: plain dents up to 6 % of the outer diameter; dents affecting a ductile girth or seam
 * weld up to 2 % of the outer diameter (and 4 % strain); dents with gouges, cracks or corrosion are
 * injurious regardless of depth and need a specific assessment (e.g. the EPRG dent–gouge model, not included).
 *
 * Reference: ASME B31.8-2018 §851.4.1; ASME B31.4 §451.6.2.4.
 *
 * @param {object} o
 * @param {number} o.d dent depth [m]
 * @param {number} o.D outer diameter [m]
 * @param {boolean} [o.onWeld=false] dent affects a weld
 * @param {boolean} [o.withStressRaiser=false] gouge, groove, crack or metal loss in the dent
 * @returns {{depthRatio:number, limit:number, acceptable:boolean, note:string}}
 */
export function dentDepthCheck({ d, D, onWeld = false, withStressRaiser = false }) {
  nonNeg('d', d); pos('D', D);
  const limit = onWeld ? 0.02 : 0.06;
  const depthRatio = d / D;
  if (withStressRaiser) return { depthRatio, limit: 0, acceptable: false, note: 'dent with a stress raiser is injurious; repair or assess specifically' };
  return { depthRatio, limit, acceptable: depthRatio <= limit, note: onWeld ? 'dent on weld: 2 % OD' : 'plain dent: 6 % OD' };
}

/**
 * Extreme-value extrapolation of wall loss from a partially inspected pipe. The deepest pit of each
 * inspected block (equal areas) is fitted with a Gumbel distribution (method of moments:
 * scale = s √6/π, loc = mean − 0.5772 scale) and the deepest pit expected over the whole pipe is the
 * characteristic largest value for a return period T = total area / block area:
 * x_T = loc − scale ln(−ln(1 − 1/T)).
 *
 * Reference: Gumbel (1958), Statistics of Extremes; ASTM G46 §8 / Kowaka (1994) for the pitting application.
 * Validity: blocks of equal size from a homogeneous corrosion population; at least 5 block maxima.
 *
 * @param {object} o
 * @param {ArrayLike<number>} o.blockMaxima deepest wall loss in each inspected block [m]
 * @param {number} o.returnPeriod T = total area (or length) / area (or length) of one block [-], > 1
 * @param {number} [o.t] wall thickness [m] (for the remaining-wall output)
 * @returns {{loc:number, scale:number, mean:number, std:number, predictedMax:number, remainingWall:(number|null),
 *   exceedanceOfWall:(number|null)}} depths [m]; exceedanceOfWall probability that the pipe-wide maximum exceeds t
 */
export function wallLossExtreme({ blockMaxima, returnPeriod, t }) {
  const n = blockMaxima.length;
  if (!(n >= 5)) throw new RangeError('blockMaxima must hold at least 5 values');
  if (!(returnPeriod > 1)) throw new RangeError('returnPeriod must exceed 1');
  let mean = 0;
  for (let i = 0; i < n; i++) mean += nonNeg('blockMaxima', blockMaxima[i]);
  mean /= n;
  let v = 0;
  for (let i = 0; i < n; i++) v += (blockMaxima[i] - mean) ** 2;
  const std = Math.sqrt(v / (n - 1));
  if (!(std > 0)) throw new RangeError('blockMaxima must not all be equal');
  const scale = (std * Math.sqrt(6)) / Math.PI;
  const loc = mean - 0.5772156649 * scale;
  const predictedMax = loc - scale * Math.log(-Math.log(1 - 1 / returnPeriod));
  let remainingWall = null;
  let exceed = null;
  if (t !== undefined) {
    pos('t', t);
    remainingWall = t - predictedMax;
    // maximum over T blocks: F_T(x) = exp(−T e^(−(x − loc)/scale))
    exceed = 1 - Math.exp(-returnPeriod * Math.exp(-(t - loc) / scale));
  }
  return { loc, scale, mean, std, predictedMax, remainingWall, exceedanceOfWall: exceed };
}

/**
 * Hook for crack-like flaws: passes the flaw to a caller-supplied fracture assessment (e.g. a failure
 * assessment diagram per BS 7910 / API 579 from the fracture engine) and normalises the result.
 *
 * Reference: none (dispatcher); corrosion and dent rules of this file do not apply to cracks.
 *
 * @param {object} o
 * @param {object} o.flaw flaw description forwarded unchanged
 * @param {function(object):{acceptable:boolean, uc?:number}} o.assess caller's fracture assessment
 * @returns {{acceptable:boolean, uc:(number|null), detail:object}}
 */
export function crackAssessment({ flaw, assess }) {
  if (typeof assess !== 'function') throw new RangeError('assess must be a function (fracture engine)');
  const r = assess(flaw);
  if (!r || typeof r.acceptable !== 'boolean') throw new RangeError('assess must return { acceptable: boolean, ... }');
  return { acceptable: r.acceptable, uc: Number.isFinite(r.uc) ? r.uc : null, detail: r };
}

/**
 * Fitness-for-service summary of a corrosion defect across the assessment levels available here:
 * Level 1 original B31G and Modified B31G, Level 2 effective area (when a profile is given), and
 * DNV-RP-F101 Part B (when the tensile strength is given). Reports each safe pressure against the MAOP and
 * the least conservative acceptable level.
 *
 * Reference: ASME B31G-2012 §1.7 (level hierarchy); DNV-RP-F101 Part B.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.d [m]
 * @param {number} o.L [m]
 * @param {number} o.smys [Pa]
 * @param {number} o.maop [Pa]
 * @param {number} [o.smts] [Pa] (enables the DNV-RP-F101 Part B line)
 * @param {{x:ArrayLike<number>, depth:ArrayLike<number>}} [o.profile] river-bottom profile (enables Level 2)
 * @param {number} [o.safetyFactor=1.39] [-]
 * @param {number} [o.F2=0.72] [-]
 * @returns {{levels:Array<{name:string, safePressure:number, uc:number, acceptable:boolean}>, acceptable:boolean,
 *   governingLevel:(string|null), depthRatio:number, action:string}}
 */
export function fitnessForService({ D, t, d, L, smys, maop, smts, profile, safetyFactor = 1.39, F2 = 0.72 }) {
  pos('maop', maop);
  const levels = [];
  const add = (name, p) => levels.push({ name, safePressure: p, uc: maop / p, acceptable: p >= maop });
  add('Level 1 — B31G original', b31gOriginal({ D, t, d, L, smys, safetyFactor }).safePressure);
  add('Level 1 — Modified B31G', b31gModified({ D, t, d, L, smys, safetyFactor }).safePressure);
  if (profile) add('Level 2 — effective area', effectiveArea({ D, t, x: profile.x, depth: profile.depth, smys, safetyFactor }).safePressure);
  if (smts !== undefined) add('DNV-RP-F101 Part B', dnvF101Single({ D, t, d, L, fu: smts, F2 }).safePressure);
  const depthRatio = d / t;
  const ok = levels.filter((l) => l.acceptable);
  const acceptable = ok.length > 0 && depthRatio <= 0.8;
  let action;
  if (depthRatio > 0.8) action = 'depth exceeds 80 % of the wall: repair or replace';
  else if (acceptable) action = 'acceptable at the current MAOP; set the re-inspection interval from the corrosion rate';
  else action = 'not acceptable at the current MAOP: de-rate to the highest safe pressure, repair, or carry out a Level 3 assessment';
  return { levels, acceptable, governingLevel: ok.length ? ok[0].name : null, depthRatio, action };
}
