/**
 * Steady-state current rating of power cables (100 % load factor).
 *
 * References
 *  - IEC 60287-1-1:2006 (+A1:2014) Electric cables — Calculation of the current rating — Part 1-1: current rating
 *    equations and calculation of losses.
 *  - IEC 60287-2-1:2006 Part 2-1: calculation of thermal resistance.
 *  - Neher & McGrath (1957), AIEE Trans. 76, "The calculation of the temperature rise and load capability of cable
 *    systems" (backfill geometric factor).
 *  - CIGRE TB 880 (2022) notes that the IEC 60287 armour-loss formula overestimates the loss of three-core
 *    wire-armoured submarine cables; a user reduction factor is exposed for that reason.
 *
 * Scope: single-core cables (single, trefoil, flat; solid, single-point or cross bonding; no magnetic armour in
 * AC), three-core cables with separately sheathed cores (SL type) and wire armour, and DC cables.
 * Not covered: three-core cables with a common sheath, pipe-type cables, magnetic armour on single-core AC cables,
 * steel-tape armour.
 *
 * SI units throughout: m, Ω/m, W/m, K·m/W, V, A; temperatures in °C where the name ends in C.
 */

import {
  dcResistance20, conductorDiameter, sheathResistance20, armourResistance20, armourWireCount,
  skinProximityCoefficients, metalProperties, insulationProperties,
} from './conductors.js';

const TWO_PI = 2 * Math.PI;

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Conductor AC resistance                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Skin-effect factor ys.
 * xs² = 8πf/R′·10⁻⁷·ks; ys = xs⁴/(192 + 0.8·xs⁴) for xs ≤ 2.8;
 * ys = −0.136 − 0.0177·xs + 0.0563·xs² for 2.8 < xs ≤ 3.8; ys = 0.354·xs − 0.733 for xs > 3.8.
 * Reference: IEC 60287-1-1 §2.1.2 (the two upper ranges were added by Amendment 1:2014).
 * @param {object} o
 * @param {number} o.Rdc DC resistance of the conductor at operating temperature R′, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} [o.ks=1] skin-effect coefficient (IEC 60287-1-1 Table 2)
 * @returns {{xs:number, ys:number}} dimensionless
 */
export function skinEffect({ Rdc, f, ks = 1 }) {
  pos(Rdc, 'Rdc'); nonNeg(f, 'f'); pos(ks, 'ks');
  const xs2 = ((8 * Math.PI * f) / Rdc) * 1e-7 * ks;
  const xs = Math.sqrt(xs2);
  const xs4 = xs2 * xs2;
  let ys;
  if (xs <= 2.8) ys = xs4 / (192 + 0.8 * xs4);
  else if (xs <= 3.8) ys = -0.136 - 0.0177 * xs + 0.0563 * xs2;
  else ys = 0.354 * xs - 0.733;
  return { xs, ys };
}

/**
 * Proximity-effect factor yp for circular conductors.
 * xp² = 8πf/R′·10⁻⁷·kp, F = xp⁴/(192 + 0.8·xp⁴);
 * three-core cables / three single-core cables: yp = F·(dc/s)²·[0.312·(dc/s)² + 1.18/(F + 0.27)];
 * two-core cables / two single-core cables: yp = 2.9·F·(dc/s)².
 * Reference: IEC 60287-1-1 §2.1.4.1 and §2.1.3. Valid for xp ≤ 2.8.
 * @param {object} o
 * @param {number} o.Rdc DC resistance at operating temperature, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} [o.kp=1] proximity-effect coefficient (IEC 60287-1-1 Table 2)
 * @param {number} o.dc conductor diameter, m
 * @param {number} o.s distance between conductor axes, m
 * @param {2|3} [o.conductors=3] number of conductors in the group
 * @returns {{xp:number, yp:number}} dimensionless
 */
export function proximityEffect({ Rdc, f, kp = 1, dc, s, conductors = 3 }) {
  pos(Rdc, 'Rdc'); nonNeg(f, 'f'); pos(kp, 'kp'); pos(dc, 'dc'); pos(s, 's');
  if (s < dc) throw new RangeError('s: axial spacing cannot be smaller than the conductor diameter dc');
  if (conductors !== 2 && conductors !== 3) throw new RangeError('conductors: must be 2 or 3');
  const xp2 = ((8 * Math.PI * f) / Rdc) * 1e-7 * kp;
  const xp4 = xp2 * xp2;
  const F = xp4 / (192 + 0.8 * xp4);
  const r2 = (dc / s) * (dc / s);
  const yp = conductors === 2 ? 2.9 * F * r2 : F * r2 * (0.312 * r2 + 1.18 / (F + 0.27));
  return { xp: Math.sqrt(xp2), yp };
}

/**
 * AC resistance of a conductor at operating temperature, R = R′·(1 + ys + yp).
 * Reference: IEC 60287-1-1 §2.1.
 * @param {object} o
 * @param {number} o.R20 DC resistance at 20 °C, Ω/m
 * @param {number} o.thetaC conductor temperature, °C
 * @param {number} o.alpha20 temperature coefficient at 20 °C, 1/K
 * @param {number} o.f frequency, Hz (0 for DC)
 * @param {number} [o.ks=1] skin-effect coefficient
 * @param {number} [o.kp=1] proximity-effect coefficient
 * @param {number} [o.dc] conductor diameter, m (needed for proximity effect)
 * @param {number} [o.s] axial spacing of conductors, m (omit for an isolated conductor: yp = 0)
 * @param {2|3} [o.conductors=3] conductors in the group
 * @returns {{Rdc:number, ys:number, yp:number, Rac:number}} Ω/m and dimensionless factors
 */
export function acResistance({ R20, thetaC, alpha20, f, ks = 1, kp = 1, dc, s, conductors = 3 }) {
  pos(R20, 'R20');
  const Rdc = R20 * (1 + alpha20 * (thetaC - 20));
  if (!(Rdc > 0)) throw new RangeError('thetaC: gives a non-positive resistance');
  if (!(f > 0)) return { Rdc, ys: 0, yp: 0, Rac: Rdc };
  const { ys } = skinEffect({ Rdc, f, ks });
  const yp = s ? proximityEffect({ Rdc, f, kp, dc, s, conductors }).yp : 0;
  return { Rdc, ys, yp, Rac: Rdc * (1 + ys + yp) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Dielectric loss                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Capacitance per core of a circular-conductor cable, C = ε/(18·ln(Di/dc))·10⁻⁹ F/m.
 * Reference: IEC 60287-1-1 §2.2. (The constant 18 is the rounded value of 1/(2πε0)·10⁻⁹.)
 * @param {object} o
 * @param {number} o.epsR relative permittivity of the insulation
 * @param {number} o.Di external diameter of the insulation (excluding screen), m
 * @param {number} o.dc diameter of the conductor including its screen, m
 * @returns {number} capacitance, F/m
 */
export function capacitance({ epsR, Di, dc }) {
  pos(epsR, 'epsR'); pos(Di, 'Di'); pos(dc, 'dc');
  if (Di <= dc) throw new RangeError('Di: must exceed dc');
  return (epsR / (18 * Math.log(Di / dc))) * 1e-9;
}

/**
 * Dielectric loss per unit length per phase, Wd = ω·C·U0²·tanδ.
 * Reference: IEC 60287-1-1 §2.2.
 * @param {object} o
 * @param {number} o.f frequency, Hz
 * @param {number} o.C capacitance per core, F/m
 * @param {number} o.U0 voltage to earth, V
 * @param {number} o.tanDelta loss factor of the insulation
 * @returns {number} dielectric loss, W/m
 */
export function dielectricLoss({ f, C, U0, tanDelta }) {
  nonNeg(f, 'f'); nonNeg(C, 'C'); nonNeg(U0, 'U0'); nonNeg(tanDelta, 'tanDelta');
  return TWO_PI * f * C * U0 * U0 * tanDelta;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Sheath and armour loss factors                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Eddy-current loss factor λ1″ of the sheath of a single-core cable (sheaths single-point or cross bonded).
 * λ1″ = Rs/R·[gs·λ0·(1 + Δ1 + Δ2) + (β1·ts)⁴/12] with ts in metres (the standard writes (β1·ts)⁴/(12·10¹²)
 * with ts in mm), gs = 1 + (ts/Ds)^1.74·(β1·Ds − 1.6) with Ds in metres, β1 = sqrt(4πω/(10⁷·ρs)),
 * m = ω/Rs·10⁻⁷, and λ0, Δ1, Δ2 for trefoil, flat-centre, flat-outer-leading and flat-outer-lagging cables.
 * Reference: IEC 60287-1-1 §2.3.6.1. The Δ exponents are reproduced from the standard; confirm against the
 * current edition before using the result for a final design.
 * @param {object} o
 * @param {number} o.R conductor AC resistance at operating temperature, Ω/m
 * @param {number} o.Rs sheath resistance at operating temperature, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} o.d mean sheath diameter, m
 * @param {number} o.s axial spacing of the cables, m
 * @param {number} o.thickness sheath thickness ts, m
 * @param {number} o.outerDiameter external sheath diameter Ds, m
 * @param {number} o.rhoS sheath resistivity at operating temperature, Ω·m
 * @param {'trefoil'|'flatCentre'|'flatOuterLeading'|'flatOuterLagging'} o.position cable position
 * @returns {{lambda1Eddy:number, lambda0:number, delta1:number, delta2:number, gs:number, m:number}}
 */
export function sheathEddyLossFactor({ R, Rs, f, d, s, thickness, outerDiameter, rhoS, position }) {
  pos(R, 'R'); pos(Rs, 'Rs'); pos(f, 'f'); pos(d, 'd'); pos(s, 's'); pos(thickness, 'thickness');
  pos(outerDiameter, 'outerDiameter'); pos(rhoS, 'rhoS');
  const w = TWO_PI * f;
  const m = (w / Rs) * 1e-7;
  const q = d / (2 * s);
  const mm = (m * m) / (1 + m * m);
  let lambda0; let delta1; let delta2;
  switch (position) {
    case 'trefoil':
      lambda0 = 3 * mm * q * q;
      delta1 = (1.14 * Math.pow(m, 2.45) + 0.33) * Math.pow(q, 0.92 * m + 1.66);
      delta2 = 0;
      break;
    case 'flatCentre':
      lambda0 = 6 * mm * q * q;
      delta1 = 0.86 * Math.pow(m, 3.08) * Math.pow(q, 1.4 * m + 0.7);
      delta2 = 0;
      break;
    case 'flatOuterLeading':
      lambda0 = 1.5 * mm * q * q;
      delta1 = 4.7 * Math.pow(m, 0.7) * Math.pow(q, 0.16 * m + 2);
      delta2 = 21 * Math.pow(m, 3.3) * Math.pow(q, 1.47 * m + 5.06);
      break;
    case 'flatOuterLagging':
      lambda0 = 1.5 * mm * q * q;
      delta1 = (-0.74 * (m + 2) * Math.sqrt(m)) / (2 + (m - 0.3) * (m - 0.3)) * Math.pow(q, m + 1);
      delta2 = 0.92 * Math.pow(m, 3.7) * Math.pow(q, m + 2);
      break;
    default:
      throw new RangeError('position: must be trefoil, flatCentre, flatOuterLeading or flatOuterLagging');
  }
  // The standard allows Δ1 = Δ2 = 0 when m ≤ 0.1 (lead sheaths); the full expressions are evaluated regardless.
  const beta1 = Math.sqrt((4 * Math.PI * w) / (1e7 * rhoS));
  const gs = 1 + Math.pow(thickness / outerDiameter, 1.74) * (beta1 * outerDiameter - 1.6);
  const lambda1Eddy = (Rs / R) * (gs * lambda0 * (1 + delta1 + delta2) + Math.pow(beta1 * thickness, 4) / 12);
  return { lambda1Eddy, lambda0, delta1, delta2, gs, m };
}

/**
 * Sheath loss factor λ1 = λ1′ + λ1″ of single-core cables.
 *
 * Solid (both-ends) bonding, λ1″ = 0 (IEC neglects eddy loss here except for large segmental conductors):
 *  - trefoil, or two cables:  λ1′ = Rs/R · 1/(1 + (Rs/X)²),  X = 2ω·10⁻⁷·ln(2s/d)                (§2.3.1)
 *  - flat, regularly transposed: same with X1 = 2ω·10⁻⁷·ln(2·∛2·s/d)                            (§2.3.2)
 *  - flat, not transposed: the three-cable expressions with P = X + Xm, Q = X − Xm/3,
 *    Xm = 2ω·10⁻⁷·ln 2 for the leading outer, centre and lagging outer cable                     (§2.3.3)
 * Single-point or cross bonding: λ1′ = 0 and λ1″ from sheathEddyLossFactor when sheath dimensions are given
 * (§2.3.6.1); for flat formation the largest of the three positions is returned as `lambda1`.
 *
 * @param {object} o
 * @param {number} o.R conductor AC resistance at operating temperature, Ω/m
 * @param {number} o.Rs sheath (or sheath ∥ non-magnetic armour) resistance at operating temperature, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} o.d mean diameter of the sheath, m
 * @param {number} o.s axial spacing between adjacent cables, m
 * @param {'trefoil'|'flat'|'pair'} [o.formation='trefoil']
 * @param {'solid'|'singlePoint'|'cross'} [o.bonding='solid']
 * @param {boolean} [o.transposed=false] regular transposition of a flat formation
 * @param {{thickness:number, outerDiameter:number, rhoS:number}} [o.sheath] sheath data for the eddy-current
 *        term (m, m, Ω·m); when omitted with single-point/cross bonding λ1″ is taken as 0
 * @returns {{lambda1:number, lambda1Circulating:number, lambda1Eddy:number, X:number,
 *            perCable?:{leading:number, centre:number, lagging:number}}}
 *          loss factors (ratio of sheath loss to conductor loss) and sheath reactance X in Ω/m
 */
export function sheathLossSingleCore({ R, Rs, f, d, s, formation = 'trefoil', bonding = 'solid', transposed = false, sheath }) {
  pos(R, 'R'); pos(Rs, 'Rs'); pos(f, 'f'); pos(d, 'd'); pos(s, 's');
  if (s < d) throw new RangeError('s: axial spacing cannot be smaller than the mean sheath diameter d');
  if (!['trefoil', 'flat', 'pair'].includes(formation)) throw new RangeError('formation: must be trefoil, flat or pair');
  if (!['solid', 'singlePoint', 'cross'].includes(bonding)) throw new RangeError('bonding: must be solid, singlePoint or cross');
  const w = TWO_PI * f;
  const X = 2 * w * 1e-7 * Math.log((2 * s) / d);
  if (bonding !== 'solid') {
    if (!sheath) return { lambda1: 0, lambda1Circulating: 0, lambda1Eddy: 0, X };
    const base = { R, Rs, f, d, s, thickness: sheath.thickness, outerDiameter: sheath.outerDiameter, rhoS: sheath.rhoS };
    if (formation === 'flat') {
      const leading = sheathEddyLossFactor({ ...base, position: 'flatOuterLeading' }).lambda1Eddy;
      const centre = sheathEddyLossFactor({ ...base, position: 'flatCentre' }).lambda1Eddy;
      const lagging = sheathEddyLossFactor({ ...base, position: 'flatOuterLagging' }).lambda1Eddy;
      const e = Math.max(leading, centre, lagging);
      return { lambda1: e, lambda1Circulating: 0, lambda1Eddy: e, X, perCable: { leading, centre, lagging } };
    }
    const e = sheathEddyLossFactor({ ...base, position: 'trefoil' }).lambda1Eddy;
    return { lambda1: e, lambda1Circulating: 0, lambda1Eddy: e, X };
  }
  if (formation !== 'flat') {
    const l = (Rs / R) / (1 + (Rs / X) * (Rs / X));
    return { lambda1: l, lambda1Circulating: l, lambda1Eddy: 0, X };
  }
  if (transposed) {
    const X1 = 2 * w * 1e-7 * Math.log((2 * Math.cbrt(2) * s) / d);
    const l = (Rs / R) / (1 + (Rs / X1) * (Rs / X1));
    return { lambda1: l, lambda1Circulating: l, lambda1Eddy: 0, X: X1 };
  }
  const Xm = 2 * w * 1e-7 * Math.LN2;
  const P = X + Xm;
  const Q = X - Xm / 3;
  const Rs2 = Rs * Rs;
  const a = (0.75 * P * P) / (Rs2 + P * P) + (0.25 * Q * Q) / (Rs2 + Q * Q);
  const b = (2 * Rs * P * Q * Xm) / (Math.sqrt(3) * (Rs2 + Q * Q) * (Rs2 + P * P));
  const leading = (Rs / R) * (a + b);
  const lagging = (Rs / R) * (a - b);
  const centre = (Rs / R) * (Q * Q) / (Rs2 + Q * Q);
  const l = Math.max(leading, centre, lagging);
  return { lambda1: l, lambda1Circulating: l, lambda1Eddy: 0, X, perCable: { leading, centre, lagging } };
}

/**
 * Sheath loss factor of a three-core cable with separate lead sheaths (SL type), armoured:
 * λ1′ = Rs/R · 1.5/(1 + (Rs/X)²), X = 2ω·10⁻⁷·ln(2s/d).
 * Reference: IEC 60287-1-1 §2.3.10. (Unarmoured SL cables use the trefoil single-core expression, factor 1.)
 * @param {object} o
 * @param {number} o.R conductor AC resistance at operating temperature, Ω/m
 * @param {number} o.Rs resistance of one sheath at operating temperature, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} o.d mean sheath diameter, m
 * @param {number} o.s distance between conductor axes, m
 * @param {boolean} [o.armoured=true]
 * @returns {{lambda1:number, X:number}} loss factor and sheath reactance (Ω/m)
 */
export function sheathLossThreeCoreSL({ R, Rs, f, d, s, armoured = true }) {
  pos(R, 'R'); pos(Rs, 'Rs'); pos(f, 'f'); pos(d, 'd'); pos(s, 's');
  const X = 2 * TWO_PI * f * 1e-7 * Math.log((2 * s) / d);
  const lambda1 = ((Rs / R) * (armoured ? 1.5 : 1)) / (1 + (Rs / X) * (Rs / X));
  return { lambda1, X };
}

/**
 * Armour loss factor λ2 of a three-core cable with steel-wire armour (round conductors):
 * λ2 = 1.23·RA/R·(2c/dA)²·1/((2.77·RA·10⁶/ω)² + 1),
 * multiplied by (1 − R/Rs·λ1′) for SL-type cables (screening by the sheath currents), where λ1′ is the
 * circulating-current factor of §2.3.1 (without the 1.5 multiplier).
 * Reference: IEC 60287-1-1 §2.4.2.3.1 and §2.4.2.5.
 * IEC 60287 is known to overestimate this loss for three-core submarine cables (CIGRE TB 880); `reductionFactor`
 * (≤ 1, project decision supported by measurement or FEA) scales the result.
 * @param {object} o
 * @param {number} o.R conductor AC resistance at operating temperature, Ω/m
 * @param {number} o.RA armour resistance at operating temperature, Ω/m
 * @param {number} o.f frequency, Hz
 * @param {number} o.c distance between the axis of a conductor and the cable centre, m
 * @param {number} o.dA mean diameter of the armour, m
 * @param {number} [o.Rs] sheath resistance for the SL-type screening factor, Ω/m
 * @param {number} [o.lambda1Trefoil] λ1′ = Rs/R/(1 + (Rs/X)²) for the SL-type screening factor
 * @param {number} [o.reductionFactor=1] user factor in (0, 1]
 * @returns {{lambda2:number, lambda2Unscreened:number, screeningFactor:number}}
 */
export function armourLossThreeCoreSteelWire({ R, RA, f, c, dA, Rs, lambda1Trefoil, reductionFactor = 1 }) {
  pos(R, 'R'); pos(RA, 'RA'); pos(f, 'f'); pos(c, 'c'); pos(dA, 'dA');
  if (!(reductionFactor > 0 && reductionFactor <= 1)) throw new RangeError('reductionFactor: must be in (0, 1]');
  const w = TWO_PI * f;
  const k = (2.77 * RA * 1e6) / w;
  const lambda2Unscreened = 1.23 * (RA / R) * Math.pow((2 * c) / dA, 2) / (k * k + 1);
  let screeningFactor = 1;
  if (Rs !== undefined && lambda1Trefoil !== undefined) screeningFactor = 1 - (R / Rs) * lambda1Trefoil;
  return { lambda2: lambda2Unscreened * screeningFactor * reductionFactor, lambda2Unscreened, screeningFactor };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Thermal resistances                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Thermal resistance of a concentric cylindrical layer, T = ρT/(2π)·ln(1 + 2t/d).
 * Used for T1 of single-core cables and screened cores (§4.1.2.1), T2 (§4.1.3.1) and T3 (§4.1.4.1) of
 * IEC 60287-2-1.
 * @param {object} o
 * @param {number} o.rhoT thermal resistivity of the layer, K·m/W
 * @param {number} o.innerDiameter diameter under the layer, m
 * @param {number} o.thickness layer thickness, m
 * @returns {number} thermal resistance per unit length, K·m/W
 */
export function layerThermalResistance({ rhoT, innerDiameter, thickness }) {
  nonNeg(rhoT, 'rhoT'); pos(innerDiameter, 'innerDiameter'); nonNeg(thickness, 'thickness');
  return (rhoT / TWO_PI) * Math.log(1 + (2 * thickness) / innerDiameter);
}

/**
 * Geometric factor Ḡ for the fillers and bedding of SL and SA type cables (curve fit of IEC 60287-2-1 Fig. 6,
 * §4.1.3.2): Ḡ = 2π·(0.00022619 + 2.11429·X − 20.4762·X²) for 0 < X ≤ 0.03 and
 * Ḡ = 2π·(0.0142108 + 1.17533·X − 4.49737·X² + 10.6352·X³) for 0.03 < X ≤ 0.15,
 * where X is the thickness of material between sheaths and armour as a fraction of the sheath outer diameter.
 * Then T2 = ρT/(6π)·Ḡ.
 * @param {number} X thickness ratio (0 … 0.15)
 * @returns {number} geometric factor Ḡ
 */
export function geometricFactorSL(X) {
  if (!(X >= 0 && X <= 0.15)) throw new RangeError('X: thickness ratio must be within 0–0.15 (supply T2 directly otherwise)');
  if (X === 0) return 0;
  if (X <= 0.03) return TWO_PI * (0.00022619 + 2.11429 * X - 20.4762 * X * X);
  return TWO_PI * (0.0142108 + 1.17533 * X - 4.49737 * X * X + 10.6352 * X * X * X);
}

/**
 * External thermal resistance of a single isolated buried cable (Kennelly),
 * T4 = ρT/(2π)·ln(u + sqrt(u² − 1)), u = 2L/De.
 * Reference: IEC 60287-2-1 §4.2.2.
 * @param {object} o
 * @param {number} o.rhoT soil thermal resistivity, K·m/W
 * @param {number} o.depth depth of the cable axis below the surface L, m
 * @param {number} o.De external diameter of the cable, m
 * @returns {number} T4, K·m/W
 */
export function t4Buried({ rhoT, depth, De }) {
  pos(rhoT, 'rhoT'); pos(depth, 'depth'); pos(De, 'De');
  const u = (2 * depth) / De;
  if (u < 1) throw new RangeError('depth: cable axis must be at least De/2 below the surface');
  return (rhoT / TWO_PI) * Math.log(u + Math.sqrt(u * u - 1));
}

/**
 * External thermal resistance of one cable of a group of equally loaded, identical buried cables by the method
 * of images: T4 = ρT/(2π)·{ln(u + sqrt(u² − 1)) + Σ ln(d′pk/dpk)}.
 * Reference: IEC 60287-2-1 §4.2.3.2 (superposition; valid when cables are not touching).
 * @param {object} o
 * @param {number} o.rhoT soil thermal resistivity, K·m/W
 * @param {{x:number, depth:number}[]} o.positions cable axes: horizontal coordinate and depth, m
 * @param {number} o.index index of the cable being rated
 * @param {number} o.De external cable diameter, m
 * @returns {number} T4 including mutual heating, K·m/W
 */
export function t4Group({ rhoT, positions, index, De }) {
  if (!Array.isArray(positions) || positions.length < 1) throw new RangeError('positions: must be a non-empty array');
  if (!(index >= 0 && index < positions.length)) throw new RangeError('index: out of range');
  const p = positions[index];
  let sum = t4Buried({ rhoT, depth: p.depth, De }) * (TWO_PI / rhoT);
  for (let k = 0; k < positions.length; k++) {
    if (k === index) continue;
    const q = positions[k];
    const dx = p.x - q.x;
    const d = Math.hypot(dx, p.depth - q.depth);
    if (!(d > 0)) throw new RangeError('positions: two cables coincide');
    const dImg = Math.hypot(dx, p.depth + q.depth);
    sum += Math.log(dImg / d);
  }
  return (rhoT / TWO_PI) * sum;
}

/**
 * External thermal resistance of three equally loaded single-core cables in touching trefoil, buried:
 * metallic-sheathed cables T4 = 1.5·ρT/π·(ln(2u) − 0.630); part-metallic covered or non-metallic sheathed
 * cables T4 = ρT/(2π)·(ln(2u) + 2·ln(u)). L is measured to the centre of the trefoil group and De is the
 * diameter of one cable. Reference: IEC 60287-2-1 §4.2.4.3.2 / §4.2.4.3.3.
 * @param {object} o
 * @param {number} o.rhoT soil thermal resistivity, K·m/W
 * @param {number} o.depth depth of the trefoil centre, m
 * @param {number} o.De diameter of one cable, m
 * @param {boolean} [o.metallicSheath=true]
 * @returns {number} T4, K·m/W
 */
export function t4TrefoilTouching({ rhoT, depth, De, metallicSheath = true }) {
  pos(rhoT, 'rhoT'); pos(depth, 'depth'); pos(De, 'De');
  const u = (2 * depth) / De;
  if (u < 1) throw new RangeError('depth: too shallow for the cable diameter');
  return metallicSheath
    ? ((1.5 * rhoT) / Math.PI) * (Math.log(2 * u) - 0.630)
    : (rhoT / TWO_PI) * (Math.log(2 * u) + 2 * Math.log(u));
}

/**
 * External thermal resistance of three equally loaded metallic-sheathed single-core cables laid flat and
 * touching: T4 = ρT·(0.475·ln(2u) − 0.346), valid for u ≥ 5.
 * Reference: IEC 60287-2-1 §4.2.4.2 (three cables, flat touching formation).
 * @param {object} o
 * @param {number} o.rhoT soil thermal resistivity, K·m/W
 * @param {number} o.depth depth of the cable axes, m
 * @param {number} o.De cable diameter, m
 * @returns {number} T4, K·m/W
 */
export function t4FlatTouching({ rhoT, depth, De }) {
  pos(rhoT, 'rhoT'); pos(depth, 'depth'); pos(De, 'De');
  const u = (2 * depth) / De;
  if (u < 5) throw new RangeError('depth: formula valid for u = 2L/De ≥ 5');
  return rhoT * (0.475 * Math.log(2 * u) - 0.346);
}

/**
 * Temperature rise at a buried point caused by other buried line heat sources (image method):
 * Δθ = Σ ρT/(2π)·Wk·ln(d′pk/dpk).
 * Reference: IEC 60287-2-1 §4.2.3.1 (unequally loaded cables).
 * @param {object} o
 * @param {number} o.rhoT soil thermal resistivity, K·m/W
 * @param {{x:number, depth:number}} o.target point of interest, m
 * @param {{x:number, depth:number, W:number}[]} o.sources line sources: position (m) and heat output (W/m)
 * @returns {number} temperature rise, K
 */
export function mutualHeatingRise({ rhoT, target, sources }) {
  pos(rhoT, 'rhoT');
  let sum = 0;
  for (const q of sources) {
    const dx = target.x - q.x;
    const d = Math.hypot(dx, target.depth - q.depth);
    if (!(d > 0)) throw new RangeError('sources: a source coincides with the target');
    sum += q.W * Math.log(Math.hypot(dx, target.depth + q.depth) / d);
  }
  return (rhoT / TWO_PI) * sum;
}

/**
 * Geometric factor of a rectangular backfill or duct bank (Neher–McGrath) and the resulting correction to T4.
 * Equivalent radius rb: ln(rb) = 0.5·(x/y)·(4/π − x/y)·ln(1 + y²/x²) + ln(x/2), x the shorter and y the longer
 * side (valid for y/x < 3); Gb = ln(u + sqrt(u² − 1)), u = LG/rb.
 * T4 is first computed as if all the soil had the backfill resistivity ρc, then the correction
 * ΔT4 = N/(2π)·(ρe − ρc)·Gb is added (N loaded cables in the backfill).
 * Reference: IEC 60287-2-1 §4.2.7.4; Neher & McGrath (1957).
 * @param {object} o
 * @param {number} o.width backfill width, m
 * @param {number} o.height backfill height, m
 * @param {number} o.depthCentre depth of the backfill centre LG, m
 * @param {number} [o.rhoSoil] native soil resistivity ρe, K·m/W
 * @param {number} [o.rhoBackfill] backfill resistivity ρc, K·m/W
 * @param {number} [o.nCables=1] number of loaded cables in the backfill N
 * @returns {{rb:number, Gb:number, deltaT4:number}} equivalent radius (m), geometric factor, correction (K·m/W;
 *          0 when resistivities are not given)
 */
export function backfillCorrection({ width, height, depthCentre, rhoSoil, rhoBackfill, nCables = 1 }) {
  pos(width, 'width'); pos(height, 'height'); pos(depthCentre, 'depthCentre');
  const x = Math.min(width, height);
  const y = Math.max(width, height);
  if (y / x >= 3) throw new RangeError('width/height: equivalent-radius formula valid for aspect ratio below 3');
  const r = x / y;
  const rb = Math.exp(0.5 * r * (4 / Math.PI - r) * Math.log(1 + (y * y) / (x * x)) + Math.log(x / 2));
  const u = depthCentre / rb;
  if (u < 1) throw new RangeError('depthCentre: backfill must lie below the surface');
  const Gb = Math.log(u + Math.sqrt(u * u - 1));
  const deltaT4 = rhoSoil !== undefined && rhoBackfill !== undefined
    ? (nCables / TWO_PI) * (rhoSoil - rhoBackfill) * Gb : 0;
  return { rb, Gb, deltaT4 };
}

/**
 * Thermal resistance of the air space between a cable and its duct or pipe,
 * T4′ = U/(1 + 0.1·(V + Y·θm)·De) with De in mm (converted internally from metres).
 * Reference: IEC 60287-2-1 §4.2.7.2 and Table 4. Constants known with certainty are provided through `kind`
 * (metallic conduit U = 5.2, V = 1.4, Y = 0.011; plastic duct U = 1.87, V = 0.312, Y = 0.0037); other duct types
 * must be supplied as U, V, Y. Valid for 25 mm ≤ De ≤ 100 mm (single cable) per the standard.
 * @param {object} o
 * @param {number} o.De cable (or equivalent group) diameter, m
 * @param {number} o.thetaMeanC mean temperature of the medium in the duct, °C
 * @param {'metallic'|'plastic'} [o.kind='plastic']
 * @param {number} [o.U] @param {number} [o.V] @param {number} [o.Y] overrides of the constants
 * @returns {number} T4′, K·m/W
 */
export function ductAirResistance({ De, thetaMeanC, kind = 'plastic', U, V, Y }) {
  pos(De, 'De');
  const table = { metallic: [5.2, 1.4, 0.011], plastic: [1.87, 0.312, 0.0037] };
  if (!table[kind]) throw new RangeError('kind: must be metallic or plastic (or give U, V, Y)');
  const u = U ?? table[kind][0];
  const v = V ?? table[kind][1];
  const y = Y ?? table[kind][2];
  return u / (1 + 0.1 * (v + y * thetaMeanC) * De * 1e3);
}

/**
 * Heat-dissipation coefficient h for cables in free air, h = Z/De^g + E (De in metres), W/(m²·K^(5/4)).
 * Reference: IEC 60287-2-1 §4.2.1.1 and Table 2 (black surfaces, De ≤ 0.15 m). Constants provided for a single
 * cable (Z 0.21, E 3.94, g 0.60), three cables in trefoil (0.96, 1.25, 0.20) and three cables touching,
 * horizontal (0.62, 1.95, 0.25), all clear of a vertical wall; other arrangements through Z, E, g.
 * @param {object} o
 * @param {number} o.De cable diameter, m
 * @param {'single'|'trefoil'|'flatTouching'} [o.arrangement='single']
 * @param {number} [o.Z] @param {number} [o.E] @param {number} [o.g] overrides of the constants
 * @returns {number} h, W/(m²·K^(5/4))
 */
export function freeAirDissipationCoefficient({ De, arrangement = 'single', Z, E, g }) {
  pos(De, 'De');
  const table = { single: [0.21, 3.94, 0.60], trefoil: [0.96, 1.25, 0.20], flatTouching: [0.62, 1.95, 0.25] };
  if (!table[arrangement]) throw new RangeError('arrangement: must be single, trefoil or flatTouching (or give Z, E, g)');
  const z = Z ?? table[arrangement][0];
  const e = E ?? table[arrangement][1];
  const gg = g ?? table[arrangement][2];
  return z / Math.pow(De, gg) + e;
}

/**
 * Surface temperature rise of a cable in free air, with solar radiation, from the surface energy balance
 * π·De·h·Δθs^(5/4) = W + σ·De·H, which is the fixed point of the iteration in IEC 60287-2-1 §4.2.1.
 * @param {object} o
 * @param {number} o.W total heat leaving the cable, W/m
 * @param {number} o.De cable diameter, m
 * @param {number} o.h dissipation coefficient, W/(m²·K^(5/4))
 * @param {number} [o.absorptivity=0] solar absorption coefficient σ of the surface (0.4 PE, 0.6 PVC, 0.8 bitumen/jute)
 * @param {number} [o.solarIntensity=0] solar radiation intensity H, W/m² (1000 W/m² is the usual design value)
 * @returns {number} Δθs surface rise above ambient, K
 */
export function freeAirSurfaceRise({ W, De, h, absorptivity = 0, solarIntensity = 0 }) {
  nonNeg(W, 'W'); pos(De, 'De'); pos(h, 'h'); nonNeg(absorptivity, 'absorptivity'); nonNeg(solarIntensity, 'solarIntensity');
  return Math.pow((W + absorptivity * De * solarIntensity) / (Math.PI * De * h), 0.8);
}

const SIGMA_SB = 5.670374419e-8; // Stefan–Boltzmann constant, W/(m²·K⁴)

/**
 * Cable inside a vertical J-tube or I-tube air section: steady energy balance between cable surface, tube and
 * ambient. Radiation across the annulus uses the exact grey-body exchange between long concentric cylinders,
 * q = σ·π·De·(Ts⁴ − Tt⁴)/(1/εc + (De/Dti)·(1/εt − 1)); convection in the annulus and outside the tube, and
 * external radiation, use coefficients supplied by the caller (no correlation is built in). The tube wall is
 * conductive (ρ of steel ≈ 0.02 K·m/W, negligible) and absorbs solar radiation on its projected area.
 * Solved by nested bisection. First-principles model — not an IEC 60287 clause.
 * @param {object} o
 * @param {number} o.W heat leaving the cable, W/m
 * @param {number} o.De cable diameter, m
 * @param {number} o.tubeInnerDiameter m
 * @param {number} o.tubeOuterDiameter m
 * @param {number} o.ambientC ambient air temperature, °C
 * @param {number} [o.hGap=3] convective coefficient cable surface → tube across the air gap, referred to the
 *        cable surface, W/(m²·K) (typical 1–5 for a sealed tube)
 * @param {number} [o.hOut=10] convective coefficient tube → ambient air, W/(m²·K) (5 still air … 30 windy)
 * @param {number} [o.emissivityCable=0.9] @param {number} [o.emissivityTubeInner=0.8]
 * @param {number} [o.emissivityTubeOuter=0.8]
 * @param {number} [o.absorptivity=0.6] solar absorptivity of the tube outer surface
 * @param {number} [o.solarIntensity=0] W/m²
 * @returns {{surfaceC:number, tubeC:number, rise:number, T4eff:number}} cable-surface and tube temperatures (°C),
 *          cable-surface rise above ambient (K) and effective external resistance rise/W (K·m/W, Infinity if W = 0)
 */
export function jTubeSurfaceTemperature({ W, De, tubeInnerDiameter, tubeOuterDiameter, ambientC, hGap = 3, hOut = 10,
  emissivityCable = 0.9, emissivityTubeInner = 0.8, emissivityTubeOuter = 0.8, absorptivity = 0.6, solarIntensity = 0 }) {
  nonNeg(W, 'W'); pos(De, 'De'); pos(tubeInnerDiameter, 'tubeInnerDiameter'); pos(tubeOuterDiameter, 'tubeOuterDiameter');
  if (tubeInnerDiameter <= De) throw new RangeError('tubeInnerDiameter: must exceed the cable diameter De');
  if (tubeOuterDiameter < tubeInnerDiameter) throw new RangeError('tubeOuterDiameter: must be ≥ tubeInnerDiameter');
  nonNeg(hGap, 'hGap'); pos(hOut, 'hOut'); nonNeg(solarIntensity, 'solarIntensity');
  const Ta = ambientC + 273.15;
  const qIn = W + absorptivity * tubeOuterDiameter * solarIntensity;
  // tube: heat qIn leaves by convection + radiation to ambient
  const outFlux = (Tt) => Math.PI * tubeOuterDiameter * (hOut * (Tt - Ta) + emissivityTubeOuter * SIGMA_SB * (Tt ** 4 - Ta ** 4));
  let lo = Ta; let hi = Ta + 1;
  while (outFlux(hi) < qIn) hi = Ta + (hi - Ta) * 2;
  for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (outFlux(mid) < qIn) lo = mid; else hi = mid; }
  const Tt = 0.5 * (lo + hi);
  const fRad = 1 / (1 / emissivityCable + (De / tubeInnerDiameter) * (1 / emissivityTubeInner - 1));
  const gapFlux = (Ts) => Math.PI * De * (hGap * (Ts - Tt) + fRad * SIGMA_SB * (Ts ** 4 - Tt ** 4));
  lo = Tt; hi = Tt + 1;
  while (gapFlux(hi) < W) hi = Tt + (hi - Tt) * 2;
  for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (gapFlux(mid) < W) lo = mid; else hi = mid; }
  const Ts = 0.5 * (lo + hi);
  return { surfaceC: Ts - 273.15, tubeC: Tt - 273.15, rise: Ts - Ta, T4eff: W > 0 ? (Ts - Ta) / W : Infinity };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Rating equation                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Permissible current from the IEC 60287 rating equation
 *   I = sqrt[(Δθ − Wd·[0.5·T1 + n·(T2 + T3 + v·T4)] + (v − 1)·Δθx)
 *            / (R·T1 + n·R·(1 + λ1)·T2 + n·R·(1 + λ1 + λ2)·(T3 + v·T4))].
 * With v = 1 this is the buried-cable equation without drying out (§1.4.1.1); with v = ρdry/ρwet and Δθx the
 * critical soil temperature rise it is the partial-drying-out equation (§1.4.2). For DC cables set Wd, λ1, λ2
 * to zero and R to the DC resistance (§1.4.1.2).
 * Reference: IEC 60287-1-1 §1.4.
 * @param {object} o
 * @param {number} o.dTheta permissible conductor temperature rise above ambient, K
 * @param {number} o.R conductor AC (or DC) resistance at maximum temperature, Ω/m
 * @param {number} o.T1 @param {number} [o.T2=0] @param {number} [o.T3=0] @param {number} o.T4 thermal resistances, K·m/W
 * @param {number} [o.n=1] number of load-carrying conductors in the cable
 * @param {number} [o.Wd=0] dielectric loss per phase, W/m
 * @param {number} [o.lambda1=0] @param {number} [o.lambda2=0] sheath and armour loss factors
 * @param {number} [o.v=1] ratio of dry to moist soil thermal resistivity
 * @param {number} [o.dThetaX=0] critical temperature rise of the soil (θx − θa), K
 * @returns {number} current, A (0 when the dielectric loss alone uses up the permissible rise)
 */
export function ratingEquation({ dTheta, R, T1, T2 = 0, T3 = 0, T4, n = 1, Wd = 0, lambda1 = 0, lambda2 = 0, v = 1, dThetaX = 0 }) {
  pos(R, 'R'); nonNeg(T1, 'T1'); nonNeg(T2, 'T2'); nonNeg(T3, 'T3'); nonNeg(T4, 'T4'); pos(n, 'n');
  if (!Number.isFinite(dTheta)) throw new RangeError('dTheta: must be finite');
  if (!(v >= 1)) throw new RangeError('v: must be ≥ 1');
  const num = dTheta - Wd * (0.5 * T1 + n * (T2 + T3 + v * T4)) + (v - 1) * dThetaX;
  const den = R * T1 + n * R * (1 + lambda1) * T2 + n * R * (1 + lambda1 + lambda2) * (T3 + v * T4);
  if (!(den > 0)) throw new RangeError('T1..T4: total thermal resistance must be positive');
  return num > 0 ? Math.sqrt(num / den) : 0;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Cable description                                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Build the radial geometry, electrical resistances at 20 °C and internal thermal resistances of a cable from a
 * layer description.
 *
 * Layer order (inside → out): conductor, conductor screen, insulation, insulation screen (incl. bedding tapes),
 * metallic sheath/screen, core jacket, [three cores laid up], bedding (and fillers), armour, serving.
 * Three-core cables are of the separately screened/sheathed type: laid-up diameter = (1 + 2/√3)·core diameter
 * (three touching circles), conductor-axis radius c = Dcore/√3 and axial spacing s = Dcore.
 *
 * Thermal resistances (IEC 60287-2-1): T1 = ρT/(2π)·ln(1 + 2·t1/dc) with t1 between conductor and sheath
 * (§4.1.2.1); T2 single-core concentric layer (§4.1.3.1) or, three-core SL type, ρT/(6π)·Ḡ (§4.1.3.2) plus one
 * third of the per-core jacket resistance; T3 concentric layer over the armour (§4.1.4.1). Without armour the
 * bedding and serving are both part of T3 and T2 = 0. Any of T1, T2, T3 may be overridden.
 *
 * @param {object} c cable description
 * @param {1|3} [c.cores=1]
 * @param {object} c.conductor { material: 'Cu'|'Al', area (m²), diameter? (m), R20? (Ω/m), construction?
 *        (key of SKIN_PROXIMITY), ks?, kp?, fillFactor?, alpha20? }
 * @param {number} [c.conductorScreenThickness=0] m
 * @param {object} c.insulation { material? (key of INSULATION), thickness (m), rhoT?, epsR?, tanDelta?,
 *        heatCapacity? }
 * @param {number} [c.insulationScreenThickness=0] m (everything between insulation and sheath)
 * @param {object} [c.sheath] { material, thickness (m), R20? (Ω/m), rho20?, alpha20? }
 * @param {object} [c.coreJacket] { thickness (m), rhoT }
 * @param {object} [c.bedding] { thickness (m), rhoT }
 * @param {object} [c.armour] { material, wireDiameter (m), nWires?, coverage?, layFactor?, R20? (Ω/m), alpha20?,
 *        magnetic? (default true for steel) }
 * @param {object} [c.serving] { thickness (m), rhoT }
 * @param {number} [c.T1] @param {number} [c.T2] @param {number} [c.T3] overrides, K·m/W
 * @returns {object} { cores, dc, dScreen, Di, dUnderSheath, sheathMeanDiameter, sheathOuterDiameter, coreDiameter,
 *          laidUpDiameter, c, s, armourMeanDiameter, armourOuterDiameter, De (all m), R20 (Ω/m), alpha20, ks, kp,
 *          Rs20, alphaS, rhoS20, RA20, alphaA (Ω/m, 1/K, Ω·m), armourMagnetic, T1, T2, T3 (K·m/W), epsR, tanDelta,
 *          C (F/m), layers: [{ name, innerDiameter, outerDiameter, rhoT, heatCapacity, area }] }
 */
export function cableGeometry(c) {
  if (!c || !c.conductor || !c.insulation) throw new RangeError('cable: conductor and insulation are required');
  const cores = c.cores ?? 1;
  if (cores !== 1 && cores !== 3) throw new RangeError('cable.cores: must be 1 or 3');
  const cond = c.conductor;
  pos(cond.area, 'cable.conductor.area');
  const cm = metalProperties(cond.material);
  const dc = cond.diameter ?? conductorDiameter({ area: cond.area, fillFactor: cond.fillFactor ?? 0.9 });
  const R20 = cond.R20 ?? dcResistance20({ material: cond.material, area: cond.area }).R20;
  const sp = skinProximityCoefficients(cond.construction ?? 'roundStranded');
  const ks = cond.ks ?? sp.ks;
  const kp = cond.kp ?? sp.kp;

  const tcs = c.conductorScreenThickness ?? 0;
  const tis = c.insulationScreenThickness ?? 0;
  const ins = c.insulation;
  pos(ins.thickness, 'cable.insulation.thickness');
  const im = ins.material ? insulationProperties(ins.material) : {};
  const rhoTi = ins.rhoT ?? im.rhoT;
  const epsR = ins.epsR ?? im.epsR;
  const tanDelta = ins.tanDelta ?? im.tanDelta;
  if (!(rhoTi > 0)) throw new RangeError('cable.insulation.rhoT: required when no material is given');

  const dScreen = dc + 2 * tcs;
  const Di = dScreen + 2 * ins.thickness;
  const dUnderSheath = Di + 2 * tis;
  const layers = [
    { name: 'conductor', innerDiameter: 0, outerDiameter: dc, rhoT: 0, heatCapacity: cm.heatCapacity, area: cond.area },
    { name: 'insulation', innerDiameter: dc, outerDiameter: dUnderSheath, rhoT: rhoTi,
      heatCapacity: ins.heatCapacity ?? im.heatCapacity ?? 2.4e6, area: (Math.PI / 4) * (dUnderSheath ** 2 - dc ** 2) },
  ];

  let d = dUnderSheath;
  let sheathMeanDiameter = null; let sheathOuterDiameter = d; let Rs20 = null; let alphaS = 0; let rhoS20 = null;
  let sheathThickness = 0;
  if (c.sheath) {
    const sh = c.sheath;
    pos(sh.thickness, 'cable.sheath.thickness');
    const smat = metalProperties(sh.material);
    sheathThickness = sh.thickness;
    sheathMeanDiameter = d + sh.thickness;
    sheathOuterDiameter = d + 2 * sh.thickness;
    rhoS20 = sh.rho20 ?? smat.rho20;
    Rs20 = sh.R20 ?? sheathResistance20({ material: sh.material, meanDiameter: sheathMeanDiameter, thickness: sh.thickness, rho20: rhoS20 });
    alphaS = sh.alpha20 ?? smat.alpha20;
    layers.push({ name: 'sheath', innerDiameter: d, outerDiameter: sheathOuterDiameter, rhoT: 0,
      heatCapacity: smat.heatCapacity, area: Math.PI * sheathMeanDiameter * sh.thickness });
    d = sheathOuterDiameter;
  }
  let Tjacket = 0;
  if (c.coreJacket && c.coreJacket.thickness > 0) {
    Tjacket = layerThermalResistance({ rhoT: c.coreJacket.rhoT, innerDiameter: d, thickness: c.coreJacket.thickness });
    const od = d + 2 * c.coreJacket.thickness;
    layers.push({ name: 'coreJacket', innerDiameter: d, outerDiameter: od, rhoT: c.coreJacket.rhoT,
      heatCapacity: c.coreJacket.heatCapacity ?? 2.4e6, area: (Math.PI / 4) * (od * od - d * d) });
    d = od;
  }
  const coreDiameter = d;
  const laidUpDiameter = cores === 3 ? coreDiameter * (1 + 2 / Math.sqrt(3)) : coreDiameter;
  const cRadius = cores === 3 ? coreDiameter / Math.sqrt(3) : 0;
  d = laidUpDiameter;

  const bed = c.bedding && c.bedding.thickness > 0 ? c.bedding : null;
  const serv = c.serving && c.serving.thickness > 0 ? c.serving : null;
  let T2 = 0; let T3 = 0;
  let armourMeanDiameter = null; let armourOuterDiameter = null; let RA20 = null; let alphaA = 0; let armourMagnetic = false;
  const outer = [];
  if (c.armour) {
    const ar = c.armour;
    pos(ar.wireDiameter, 'cable.armour.wireDiameter');
    const tb = bed ? bed.thickness : 0;
    if (cores === 3) {
      // Fillers between the cores are treated with the bedding through the geometric factor Ḡ.
      if (c.T2 === undefined) {
        T2 = ((bed ? bed.rhoT : 0) / (6 * Math.PI)) * geometricFactorSL(tb / coreDiameter) + Tjacket / 3;
      }
    } else {
      T2 = (bed ? layerThermalResistance({ rhoT: bed.rhoT, innerDiameter: d, thickness: tb }) : 0) + Tjacket;
    }
    if (bed) {
      // For three cores the bedding/filler region is everything inside the armour that is not a core.
      const od = d + 2 * tb;
      const area = (Math.PI / 4) * (od * od - cores * coreDiameter * coreDiameter);
      outer.push({ name: 'bedding', innerDiameter: d, outerDiameter: od, rhoT: bed.rhoT,
        heatCapacity: bed.heatCapacity ?? 2.0e6, area });
      d = od;
    }
    const amat = metalProperties(ar.material);
    armourMeanDiameter = d + ar.wireDiameter;
    armourOuterDiameter = d + 2 * ar.wireDiameter;
    const nWires = ar.nWires ?? armourWireCount({ meanDiameter: armourMeanDiameter, wireDiameter: ar.wireDiameter, coverage: ar.coverage ?? 0.95 });
    RA20 = ar.R20 ?? armourResistance20({ material: ar.material, wireDiameter: ar.wireDiameter, nWires, layFactor: ar.layFactor ?? 1, rho20: ar.rho20 });
    alphaA = ar.alpha20 ?? amat.alpha20;
    armourMagnetic = ar.magnetic ?? (ar.material === 'steel');
    outer.push({ name: 'armour', innerDiameter: d, outerDiameter: armourOuterDiameter, rhoT: 0,
      heatCapacity: amat.heatCapacity, area: nWires * (Math.PI / 4) * ar.wireDiameter ** 2 });
    d = armourOuterDiameter;
    if (serv) {
      T3 = layerThermalResistance({ rhoT: serv.rhoT, innerDiameter: d, thickness: serv.thickness });
      const od = d + 2 * serv.thickness;
      outer.push({ name: 'serving', innerDiameter: d, outerDiameter: od, rhoT: serv.rhoT,
        heatCapacity: serv.heatCapacity ?? 2.0e6, area: (Math.PI / 4) * (od * od - d * d) });
      d = od;
    }
  } else {
    // No armour: jacket, bedding and serving are a single covering in series (T3); T2 = 0.
    T3 = cores === 1 ? Tjacket : Tjacket / 3;
    for (const [name, lay] of [['bedding', bed], ['serving', serv]]) {
      if (!lay) continue;
      T3 += layerThermalResistance({ rhoT: lay.rhoT, innerDiameter: d, thickness: lay.thickness });
      const od = d + 2 * lay.thickness;
      outer.push({ name, innerDiameter: d, outerDiameter: od, rhoT: lay.rhoT, heatCapacity: lay.heatCapacity ?? 2.0e6,
        area: (Math.PI / 4) * (od * od - d * d) });
      d = od;
    }
  }
  const De = d;
  const T1 = c.T1 ?? layerThermalResistance({ rhoT: rhoTi, innerDiameter: dc, thickness: tcs + ins.thickness + tis });
  if (c.T2 !== undefined) T2 = c.T2;
  if (c.T3 !== undefined) T3 = c.T3;
  const C = epsR ? capacitance({ epsR, Di, dc: dScreen }) : 0;
  return {
    cores, dc, dScreen, Di, dUnderSheath, sheathMeanDiameter, sheathOuterDiameter, sheathThickness, coreDiameter,
    laidUpDiameter, c: cRadius, s: coreDiameter, armourMeanDiameter, armourOuterDiameter, De,
    R20, alpha20: cond.alpha20 ?? cm.alpha20, ks, kp, Rs20, alphaS, rhoS20, RA20, alphaA, armourMagnetic,
    T1, T2, T3, epsR: epsR ?? 0, tanDelta: tanDelta ?? 0, C,
    layers: [...layers, ...outer],
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Installation: external thermal environment                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

function groupPositions(inst, De) {
  const L = inst.depth;
  const formation = inst.formation ?? 'single';
  const s = inst.spacing ?? De;
  if (s < De * (1 - 1e-9)) throw new RangeError('installation.spacing: cannot be smaller than the cable diameter');
  if (inst.positions) return { positions: inst.positions, touching: false, s };
  switch (formation) {
    case 'single': return { positions: [{ x: 0, depth: L }], touching: false, s };
    case 'pair': return { positions: [{ x: -s / 2, depth: L }, { x: s / 2, depth: L }], touching: s <= De * (1 + 1e-9), s };
    case 'flat': return { positions: [{ x: -s, depth: L }, { x: 0, depth: L }, { x: s, depth: L }], touching: s <= De * (1 + 1e-9), s };
    case 'trefoil': {
      // L is the depth of the centroid of the group.
      const h = s / Math.sqrt(3);
      return { positions: [{ x: 0, depth: L - h }, { x: -s / 2, depth: L + h / 2 }, { x: s / 2, depth: L + h / 2 }],
        touching: s <= De * (1 + 1e-9), s };
    }
    default: throw new RangeError('installation.formation: must be single, pair, flat or trefoil');
  }
}

/**
 * External thermal model of an installation. Returns a function giving the cable-surface temperature rise above
 * ambient for a total heat output W (W/m per cable) together with the equivalent T4 where it is constant.
 *
 * installation.type:
 *  - 'buried'  : depth, rhoSoil, formation ('single' | 'pair' | 'flat' | 'trefoil'), spacing (axial, default
 *                touching), positions (explicit, overrides formation), metallicSheath, dryZone { rhoDry,
 *                criticalRiseK }, backfill { width, height, depthCentre, rhoBackfill }, externalSources
 *                [{ x, depth, W }] (other heat sources, x relative to the rated cable/group centre).
 *                T4 by Kennelly / image method / IEC touching-trefoil and touching-flat formulas.
 *  - 'seabed'  : cable lying on the seabed in water, T4 = 1/(π·De·h) with h (W/(m²·K)) supplied.
 *  - 'air'     : IEC 60287-2-1 §4.2.1 free air; arrangement, Z/E/g, absorptivity, solarIntensity.
 *  - 'jtube'   : see jTubeSurfaceTemperature (tubeInnerDiameter, tubeOuterDiameter, hGap, hOut, emissivities,
 *                absorptivity, solarIntensity).
 *  - 'duct'    : buried duct (e.g. HDD landfall): duct { innerDiameter, outerDiameter, rhoDuct, fill: 'air' |
 *                'solid', kind / U, V, Y for air, rhoFill for bentonite or grout }, plus the buried inputs.
 *  - 'custom'  : T4 given directly.
 * @param {object} inst installation description (ambientC required)
 * @param {{De:number, sheathMeanDiameter:?number}} geom cable geometry
 * @returns {{kind:string, T4:?number, v:number, dThetaX:number, extRise:number, rise:(W:number, meanC?:number)=>number,
 *            nonlinear:boolean, detail:object}}
 */
export function externalThermalModel(inst, geom) {
  if (!inst || !Number.isFinite(inst.ambientC)) throw new RangeError('installation.ambientC: required');
  const De = geom.De;
  const type = inst.type ?? 'buried';
  const linear = (T4, extra = {}) => ({ kind: type, T4, v: 1, dThetaX: 0, extRise: 0, nonlinear: false, detail: {}, rise: (W) => W * T4, ...extra });

  if (type === 'custom') { nonNeg(inst.T4, 'installation.T4'); return linear(inst.T4); }
  if (type === 'seabed') { pos(inst.h, 'installation.h'); return linear(1 / (Math.PI * De * inst.h)); }
  if (type === 'air') {
    const h = freeAirDissipationCoefficient({ De, arrangement: inst.arrangement ?? 'single', Z: inst.Z, E: inst.E, g: inst.g });
    const absorptivity = inst.absorptivity ?? 0; const solarIntensity = inst.solarIntensity ?? 0;
    return { kind: type, T4: null, v: 1, dThetaX: 0, extRise: 0, nonlinear: true, detail: { h },
      rise: (W) => freeAirSurfaceRise({ W, De, h, absorptivity, solarIntensity }) };
  }
  if (type === 'jtube') {
    return { kind: type, T4: null, v: 1, dThetaX: 0, extRise: 0, nonlinear: true, detail: {},
      rise: (W) => jTubeSurfaceTemperature({ ...inst, W, De }).rise };
  }
  if (type !== 'buried' && type !== 'duct') throw new RangeError('installation.type: must be buried, duct, seabed, air, jtube or custom');

  pos(inst.depth, 'installation.depth'); pos(inst.rhoSoil, 'installation.rhoSoil');
  // With a thermal backfill T4 is evaluated with the backfill resistivity and corrected (Neher–McGrath).
  const rho = inst.backfill ? inst.backfill.rhoBackfill : inst.rhoSoil;
  pos(rho, 'installation.backfill.rhoBackfill');
  const Dext = type === 'duct' ? inst.duct.outerDiameter : De;
  if (type === 'duct') {
    if (!inst.duct) throw new RangeError('installation.duct: required for type "duct"');
    pos(inst.duct.innerDiameter, 'installation.duct.innerDiameter'); pos(inst.duct.outerDiameter, 'installation.duct.outerDiameter');
    if (inst.duct.innerDiameter <= De) throw new RangeError('installation.duct.innerDiameter: must exceed the cable diameter');
  }
  const { positions, touching } = groupPositions(inst, Dext);
  const formation = inst.positions ? 'custom' : (inst.formation ?? 'single');
  let T4soil; let hottest = 0;
  if (formation === 'trefoil' && touching && type === 'buried') {
    T4soil = t4TrefoilTouching({ rhoT: rho, depth: inst.depth, De, metallicSheath: inst.metallicSheath ?? (geom.sheathMeanDiameter != null) });
  } else if (formation === 'flat' && touching && type === 'buried' && (2 * inst.depth) / De >= 5
    && (inst.metallicSheath ?? (geom.sheathMeanDiameter != null))) {
    T4soil = t4FlatTouching({ rhoT: rho, depth: inst.depth, De });
    hottest = 1;
  } else {
    T4soil = -Infinity;
    for (let i = 0; i < positions.length; i++) {
      const t = t4Group({ rhoT: rho, positions, index: i, De: Dext });
      if (t > T4soil) { T4soil = t; hottest = i; }
    }
  }
  let deltaBackfill = 0;
  if (inst.backfill) {
    deltaBackfill = backfillCorrection({ ...inst.backfill, rhoSoil: inst.rhoSoil, nCables: inst.backfill.nCables ?? positions.length }).deltaT4;
    T4soil += deltaBackfill;
  }
  let extRise = 0;
  if (inst.externalSources && inst.externalSources.length) {
    extRise = mutualHeatingRise({ rhoT: inst.rhoSoil, target: positions[hottest], sources: inst.externalSources });
  }
  let v = 1; let dThetaX = 0;
  if (inst.dryZone) {
    pos(inst.dryZone.rhoDry, 'installation.dryZone.rhoDry'); pos(inst.dryZone.criticalRiseK, 'installation.dryZone.criticalRiseK');
    v = inst.dryZone.rhoDry / rho;
    if (v < 1) throw new RangeError('installation.dryZone.rhoDry: must be ≥ the moist soil resistivity');
    dThetaX = inst.dryZone.criticalRiseK;
  }
  const detail = { positions, hottest, T4soil, deltaBackfill };
  if (type === 'buried') {
    return { kind: type, T4: T4soil, v, dThetaX, extRise, nonlinear: false, detail, rise: (W) => W * T4soil };
  }
  // Duct: T4 = T4′ (cable → duct) + T4″ (duct wall) + T4‴ (outside the duct).
  const duct = inst.duct;
  const T4wall = (duct.rhoDuct ?? 3.5) / TWO_PI * Math.log(duct.outerDiameter / duct.innerDiameter);
  const solidFill = duct.fill === 'solid';
  if (solidFill) pos(duct.rhoFill, 'installation.duct.rhoFill');
  const t4In = (meanC) => (solidFill
    ? (duct.rhoFill / TWO_PI) * Math.log(duct.innerDiameter / De)
    : ductAirResistance({ De, thetaMeanC: meanC, kind: duct.kind ?? 'plastic', U: duct.U, V: duct.V, Y: duct.Y }));
  detail.T4wall = T4wall;
  detail.t4In = t4In;
  // Dry-zone factor applies to the soil part only.
  return {
    kind: type, T4: null, v, dThetaX, extRise, nonlinear: !solidFill, detail,
    ductInner: (meanC) => t4In(meanC) + T4wall,
    rise: (W, meanC = inst.ambientC + 30) => W * (t4In(meanC) + T4wall + T4soil),
    T4soil,
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Complete rating                                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

function lossFactors({ g, sys, inst, R, thetaSheathC, thetaArmourC, s }) {
  // Returns λ1, λ2 for the present temperatures.
  if (sys.dc || !(sys.f > 0)) return { lambda1: 0, lambda2: 0, lambda1Circulating: 0, lambda1Eddy: 0, Rs: null, RA: null };
  const f = sys.f;
  let Rs = g.Rs20 != null ? g.Rs20 * (1 + g.alphaS * (thetaSheathC - 20)) : null;
  const RA = g.RA20 != null ? g.RA20 * (1 + g.alphaA * (thetaArmourC - 20)) : null;
  let lambda1 = 0; let lambda2 = 0; let lambda1Circulating = 0; let lambda1Eddy = 0;
  if (g.cores === 3) {
    if (Rs != null) {
      const sl = sheathLossThreeCoreSL({ R, Rs, f, d: g.sheathMeanDiameter, s: g.s, armoured: RA != null });
      lambda1 = lambda1Circulating = sl.lambda1;
      if (RA != null && g.armourMagnetic) {
        const l1t = (Rs / R) / (1 + (Rs / sl.X) * (Rs / sl.X));
        lambda2 = armourLossThreeCoreSteelWire({ R, RA, f, c: g.c, dA: g.armourMeanDiameter, Rs, lambda1Trefoil: l1t,
          reductionFactor: sys.armourLossReduction ?? 1 }).lambda2;
      }
    } else if (RA != null && g.armourMagnetic) {
      lambda2 = armourLossThreeCoreSteelWire({ R, RA, f, c: g.c, dA: g.armourMeanDiameter, reductionFactor: sys.armourLossReduction ?? 1 }).lambda2;
    }
    if (sys.lambda2 !== undefined) lambda2 = sys.lambda2;
    return { lambda1, lambda2, lambda1Circulating, lambda1Eddy, Rs, RA };
  }
  // Single-core cables.
  if (RA != null && g.armourMagnetic) {
    throw new RangeError('cable.armour: magnetic armour on single-core AC cables is not covered (IEC 60287-1-1 §2.4.2.1); give system.lambda2 and a non-magnetic armour');
  }
  const formation = inst.formation ?? 'single';
  if (Rs != null && formation !== 'single' && !inst.positions) {
    const bonding = inst.bonding ?? 'solid';
    let Re = Rs; let d = g.sheathMeanDiameter; let share = 1;
    if (RA != null && bonding === 'solid') {
      // Non-magnetic armour bonded with the sheath: parallel resistance and rms diameter (IEC 60287-1-1 §2.4.1).
      Re = (Rs * RA) / (Rs + RA);
      d = Math.sqrt((g.sheathMeanDiameter ** 2 + g.armourMeanDiameter ** 2) / 2);
      share = RA / (Rs + RA);
    }
    const sl = sheathLossSingleCore({ R, Rs: Re, f, d, s, formation, bonding, transposed: inst.transposed ?? false,
      sheath: g.sheathThickness > 0 ? { thickness: g.sheathThickness, outerDiameter: g.sheathOuterDiameter,
        rhoS: g.rhoS20 * (1 + g.alphaS * (thetaSheathC - 20)) } : undefined });
    lambda1 = sl.lambda1 * share;
    lambda2 = sl.lambda1 * (1 - share);
    lambda1Circulating = sl.lambda1Circulating * share;
    lambda1Eddy = sl.lambda1Eddy * share;
  }
  if (sys.lambda1 !== undefined) lambda1 = sys.lambda1;
  if (sys.lambda2 !== undefined) lambda2 = sys.lambda2;
  Rs = Rs ?? null;
  return { lambda1, lambda2, lambda1Circulating, lambda1Eddy, Rs, RA };
}

function electricalState({ g, sys, inst, thetaC, thetaSheathC, thetaArmourC }) {
  const dc = !!sys.dc || !(sys.f > 0);
  // Axial spacing for proximity effect: cores of a three-core cable, or neighbouring single-core cables.
  let s;
  if (g.cores === 3) s = g.s;
  else if ((inst.formation ?? 'single') !== 'single' && !inst.positions) s = inst.spacing ?? (inst.type === 'duct' ? inst.duct.outerDiameter : g.De);
  const twoCond = (inst.formation === 'pair');
  const r = acResistance({ R20: g.R20, thetaC, alpha20: g.alpha20, f: dc ? 0 : sys.f, ks: g.ks, kp: g.kp, dc: g.dc, s,
    conductors: twoCond ? 2 : 3 });
  const lf = lossFactors({ g, sys, inst, R: r.Rac, thetaSheathC, thetaArmourC, s: s ?? g.De });
  let Wd = 0;
  if (!dc && sys.includeDielectric !== false) {
    const U0 = sys.U0 ?? (sys.U !== undefined ? sys.U / Math.sqrt(3) : 0);
    Wd = dielectricLoss({ f: sys.f, C: g.C, U0, tanDelta: sys.tanDelta ?? g.tanDelta });
  }
  return { ...r, ...lf, Wd };
}

function validateSystem(sys) {
  if (!sys) throw new RangeError('system: required');
  if (!sys.dc && !(sys.f > 0)) throw new RangeError('system.f: frequency required for AC (or set system.dc = true)');
  if (!Number.isFinite(sys.thetaMaxC)) throw new RangeError('system.thetaMaxC: maximum conductor temperature required');
}

/**
 * Conductor, sheath, armour and surface temperatures and the loss breakdown of a cable carrying a given current.
 * Resistances of conductor, sheath and armour follow their own temperatures (fixed-point iteration);
 * θs = θc − (I²R + 0.5·Wd)·T1, θA = θs − [I²R·(1 + λ1) + Wd]·n·T2, θe = θA − [I²R·(1 + λ1 + λ2) + Wd]·n·T3,
 * and the surface rise follows the installation model (with the two-zone drying-out correction when active).
 * Reference: IEC 60287-1-1 §1.4, §2; IEC 60287-2-1.
 * @param {object} o
 * @param {object} o.cable cable description (see cableGeometry) or a geometry object returned by it
 * @param {object} o.installation installation (see externalThermalModel)
 * @param {object} o.system { f (Hz), U (line-to-line, V) or U0 (V), dc?, thetaMaxC, tanDelta?, lambda1?, lambda2?,
 *        armourLossReduction?, includeDielectric? }
 * @param {number} o.I current, A
 * @returns {object} { I, conductorC, sheathC, armourC, surfaceC, ambientC, R, Rdc, ys, yp, lambda1, lambda2,
 *          Wd (W/m per core), losses: { conductor, dielectric, sheath, armour, total } (W/m per cable),
 *          T1, T2, T3, T4eff (K·m/W), dried (boolean), iterations }
 */
export function temperaturesAtCurrent({ cable, installation, system, I }) {
  validateSystem(system); nonNeg(I, 'I');
  const g = cable.layers ? cable : cableGeometry(cable);
  const ext = externalThermalModel(installation, g);
  const n = g.cores;
  let thetaC = installation.ambientC + 30; let thetaS = thetaC - 5; let thetaA = thetaS - 2; let thetaE = thetaA - 2;
  let st; let Wtot = 0; let dried = false; let surfRise = 0; let it = 0;
  for (; it < 200; it++) {
    st = electricalState({ g, sys: system, inst: installation, thetaC, thetaSheathC: thetaS, thetaArmourC: thetaA });
    const Wc = I * I * st.Rac;
    Wtot = n * (Wc * (1 + st.lambda1 + st.lambda2) + st.Wd);
    if (ext.kind === 'duct') {
      const soil = soilRise(Wtot, ext);
      dried = soil.dried;
      const ductOuterRise = soil.rise + ext.extRise;
      // Mean temperature of the duct filling = average of cable surface (previous iterate) and duct inner wall.
      const wallC = installation.ambientC + ductOuterRise + Wtot * ext.detail.T4wall;
      surfRise = ductOuterRise + Wtot * ext.ductInner(0.5 * (thetaE + wallC));
    } else if (ext.nonlinear) {
      surfRise = ext.rise(Wtot);
    } else {
      const soil = soilRise(Wtot, ext);
      dried = soil.dried;
      surfRise = soil.rise + ext.extRise;
    }
    const newE = installation.ambientC + surfRise;
    const newA = newE + (Wc * (1 + st.lambda1 + st.lambda2) + st.Wd) * n * g.T3;
    const newS = newA + (Wc * (1 + st.lambda1) + st.Wd) * n * g.T2;
    const newC = newS + (Wc + 0.5 * st.Wd) * g.T1;
    if (!Number.isFinite(newC) || newC > 5000) throw new RangeError('I: no thermal equilibrium (thermal runaway) at this current');
    const err = Math.abs(newC - thetaC);
    thetaC = newC; thetaS = newS; thetaA = newA; thetaE = newE;
    if (err < 1e-9) break;
  }
  const Wc = I * I * st.Rac;
  return {
    I, conductorC: thetaC, sheathC: thetaS, armourC: thetaA, surfaceC: thetaE, ambientC: installation.ambientC,
    R: st.Rac, Rdc: st.Rdc, ys: st.ys, yp: st.yp, lambda1: st.lambda1, lambda2: st.lambda2,
    lambda1Circulating: st.lambda1Circulating, lambda1Eddy: st.lambda1Eddy, Wd: st.Wd,
    losses: { conductor: n * Wc, dielectric: n * st.Wd, sheath: n * Wc * st.lambda1, armour: n * Wc * st.lambda2, total: Wtot },
    T1: g.T1, T2: g.T2, T3: g.T3, T4eff: Wtot > 0 ? (surfRise - (ext.extRise || 0)) / Wtot : (ext.T4 ?? NaN),
    externalRise: ext.extRise || 0, dried, iterations: it + 1,
  };
}

function soilRise(W, ext) {
  // Two-zone drying-out: when the wet-soil surface rise exceeds the critical rise the dry zone forms and the
  // rise becomes v·W·T4 − (v − 1)·Δθx (IEC 60287-1-1 §1.4.2 rearranged).
  const T4 = ext.kind === 'duct' ? ext.T4soil : ext.T4;
  const wet = W * T4;
  if (ext.v > 1 && wet > ext.dThetaX) return { rise: ext.v * wet - (ext.v - 1) * ext.dThetaX, dried: true };
  return { rise: wet, dried: false };
}

/**
 * Continuous current rating of a cable for an installation.
 * The current is found so that the conductor reaches thetaMaxC, with conductor, sheath and armour resistances at
 * their own converged temperatures. For installations with constant thermal resistances this is the fixed point
 * of the IEC 60287-1-1 §1.4 rating equation (see ratingEquation); free-air, J-tube and air-filled-duct
 * installations are non-linear and are solved on the same energy balance by bisection.
 * DC cables (system.dc): no skin/proximity effect, no dielectric, sheath or armour loss; an optional limit on the
 * temperature drop across the insulation (system.dTInsulationMaxK, used for mass-impregnated and DC-XLPE cables to
 * limit stress inversion) gives a second rating I = sqrt(ΔT/(R·T1)) and the lower one governs.
 * Flat formations with solid bonding: the sheath-loss factor of the worst cable is combined with the external
 * resistance of the hottest (centre) cable, which is conservative.
 * Reference: IEC 60287-1-1 §1.4; IEC 60287-2-1.
 * @param {object} o
 * @param {object} o.cable cable description (see cableGeometry) or geometry object
 * @param {object} o.installation installation (see externalThermalModel)
 * @param {object} o.system { f, U or U0, dc?, thetaMaxC, dTInsulationMaxK?, tanDelta?, lambda1?, lambda2?,
 *        armourLossReduction?, includeDielectric? }
 * @returns {object} { I (A), governedBy: 'conductorTemperature' | 'insulationTemperatureDrop', IThermal,
 *          IInsulationDrop, ...fields of temperaturesAtCurrent at the rating, energyBalanceResidualK }
 */
export function rateCable({ cable, installation, system }) {
  validateSystem(system);
  const g = cable.layers ? cable : cableGeometry(cable);
  const target = system.thetaMaxC;
  if (!installation || !Number.isFinite(installation.ambientC)) throw new RangeError('installation.ambientC: required');
  if (!(target > installation.ambientC)) throw new RangeError('system.thetaMaxC: must exceed the ambient temperature');
  const tempAt = (I) => temperaturesAtCurrent({ cable: g, installation, system, I }).conductorC;
  let IThermal = 0;
  if (tempAt(0) < target) {
    let lo = 0; let hi = 100;
    // Bracket: temperature rises monotonically with current; runaway counts as "too hot".
    const hot = (I) => { try { return tempAt(I) >= target; } catch (e) { if (e instanceof RangeError && /thermal/.test(e.message)) return true; throw e; } };
    while (!hot(hi)) { lo = hi; hi *= 2; if (hi > 1e6) throw new RangeError('installation: rating exceeds 1 MA, check inputs'); }
    for (let i = 0; i < 60 && hi - lo > 1e-9 * hi; i++) {
      const mid = 0.5 * (lo + hi);
      if (hot(mid)) hi = mid; else lo = mid;
    }
    IThermal = 0.5 * (lo + hi);
  }
  let I = IThermal; let governedBy = 'conductorTemperature'; let IInsulationDrop = Infinity;
  if (system.dc && system.dTInsulationMaxK !== undefined) {
    pos(system.dTInsulationMaxK, 'system.dTInsulationMaxK');
    // Drop across insulation = I²·R(θc)·T1 with θc from the thermal solution at that current.
    const drop = (x) => { const t = temperaturesAtCurrent({ cable: g, installation, system, I: x }); return x * x * t.R * g.T1; };
    if (drop(IThermal) > system.dTInsulationMaxK) {
      let lo = 0; let hi = IThermal;
      for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (drop(mid) > system.dTInsulationMaxK) hi = mid; else lo = mid; }
      IInsulationDrop = 0.5 * (lo + hi);
      I = IInsulationDrop; governedBy = 'insulationTemperatureDrop';
    }
  }
  const t = temperaturesAtCurrent({ cable: g, installation, system, I });
  // Independent check of the heat balance through the thermal network.
  const n = g.cores;
  const Wc = I * I * t.R;
  const drops = (Wc + 0.5 * t.Wd) * g.T1 + (Wc * (1 + t.lambda1) + t.Wd) * n * g.T2
    + (Wc * (1 + t.lambda1 + t.lambda2) + t.Wd) * n * g.T3 + (t.surfaceC - t.ambientC);
  return { ...t, I, governedBy, IThermal, IInsulationDrop, De: g.De,
    energyBalanceResidualK: t.conductorC - t.ambientC - drops };
}

/**
 * Rating as a function of one installation parameter (burial depth, soil resistivity, ambient temperature, ...).
 * @param {object} o
 * @param {object} o.cable @param {object} o.installation @param {object} o.system as for rateCable
 * @param {string} o.parameter key of the installation object to vary (e.g. 'depth' in m, 'rhoSoil' in K·m/W)
 * @param {number[]} o.values values of that parameter
 * @returns {{parameter:string, values:number[], I:number[], surfaceC:number[]}} ratings in A, surface temperatures in °C
 */
export function sweepRating({ cable, installation, system, parameter, values }) {
  if (typeof parameter !== 'string') throw new RangeError('parameter: must be a key of the installation object');
  if (!Array.isArray(values) && !ArrayBuffer.isView(values)) throw new RangeError('values: must be an array');
  const g = cable.layers ? cable : cableGeometry(cable);
  const I = []; const surfaceC = [];
  for (const v of values) {
    const r = rateCable({ cable: g, installation: { ...installation, [parameter]: v }, system });
    I.push(r.I); surfaceC.push(r.surfaceC);
  }
  return { parameter, values: Array.from(values), I, surfaceC };
}

/**
 * Rating of each section of a cable route and identification of the limiting section (thermal hotspot).
 * Sections are treated as long compared with the longitudinal heat-flow length (no axial conduction credit).
 * @param {object} o
 * @param {object} o.cable @param {object} o.system as for rateCable
 * @param {{name:string, length:number, installation:object}[]} o.sections route sections (length in m)
 * @param {number} [o.I] if given, conductor temperature of each section at this current is also returned, A
 * @returns {{sections:{name:string, length:number, kpStart:number, I:number, conductorAtLoadC:?number}[],
 *            limitingIndex:number, limitingName:string, I:number}} route rating in A
 */
export function routeHotspot({ cable, system, sections, I }) {
  if (!Array.isArray(sections) || sections.length === 0) throw new RangeError('sections: must be a non-empty array');
  const g = cable.layers ? cable : cableGeometry(cable);
  let kp = 0; let limitingIndex = 0; let Imin = Infinity;
  const out = sections.map((sct, i) => {
    const r = rateCable({ cable: g, installation: sct.installation, system });
    const row = { name: sct.name ?? `Section ${i + 1}`, length: sct.length ?? 0, kpStart: kp, I: r.I,
      conductorAtLoadC: I !== undefined ? safeTemp(g, sct.installation, system, I) : null };
    kp += sct.length ?? 0;
    if (r.I < Imin) { Imin = r.I; limitingIndex = i; }
    return row;
  });
  return { sections: out, limitingIndex, limitingName: out[limitingIndex].name, I: Imin };
}

function safeTemp(g, installation, system, I) {
  try { return temperaturesAtCurrent({ cable: g, installation, system, I }).conductorC; } catch (e) {
    if (e instanceof RangeError) return Infinity; throw e;
  }
}
