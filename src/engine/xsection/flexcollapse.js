/**
 * Hydrostatic collapse, radial (bird-caging) and lateral instability checks for unbonded flexible pipes.
 *
 * All ring formulas are per unit pipe length: EI in N m (N m^2 per metre), pressures in Pa. Radii are mean ring
 * radii. Screening-level models: they do not replace the manufacturer's validated collapse and buckling
 * methodology required by API Spec 17J.
 */

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function angle(v, name = 'alpha') {
  if (!(typeof v === 'number' && v > 0 && v < Math.PI / 2)) throw new RangeError(`${name} must lie strictly between 0 and pi/2 rad`);
  return v;
}

/**
 * Equivalent ring bending stiffness per unit length of an interlocked (carcass or pressure-armour) layer and the
 * thickness of the solid ring with the same stiffness (I-equivalence):
 *     EI_eq = K n E I_min / L_p,     t_eq = (12 EI_eq / E)^(1/3)
 * with n the number of profiles (tendons) in the layer, I_min the smallest second moment of area of the profile
 * cross-section, L_p the pitch and K a factor accounting for lay angle and overlap/fill of the profile.
 * Reference: API RP 17B §8.4 (collapse; equivalent ring stiffness of the carcass). K has no universal value:
 * it is profile specific (argument, default 1).
 * @param {object} o
 * @param {number} o.E Young's modulus, Pa
 * @param {number} o.Imin smallest second moment of the profile section, m^4
 * @param {number} o.pitch layer pitch L_p, m
 * @param {number} [o.nTendons=1] number of profiles wound in the layer
 * @param {number} [o.K=1] lay-angle / overlap factor, -
 * @returns {{EI:number, tEq:number}} N m (per unit length), m
 */
export function interlockedRingStiffness({ E, Imin, pitch, nTendons = 1, K = 1 }) {
  pos(E, 'E'); pos(Imin, 'Imin'); pos(pitch, 'pitch'); pos(nTendons, 'nTendons'); pos(K, 'K');
  const EI = K * nTendons * E * Imin / pitch;
  return { EI, tEq: Math.cbrt(12 * EI / E) };
}

/**
 * Ring bending stiffness per unit length of a solid (or fill-reduced) layer: EI = F E t^3 / 12, optionally with
 * the plane-strain factor 1/(1 - nu^2) for a long tube.
 * Reference: Timoshenko & Gere, Theory of Elastic Stability, §7.2-7.5.
 * @param {object} o
 * @param {number} o.E Pa
 * @param {number} o.thickness m
 * @param {number} [o.fill=1] metal fraction of the layer along the pipe axis, -
 * @param {number} [o.nu] if given, the plane-strain factor is applied
 * @returns {{EI:number}} N m
 */
export function solidRingStiffness({ E, thickness, fill = 1, nu }) {
  pos(E, 'E'); pos(thickness, 'thickness');
  if (!(fill > 0 && fill <= 1)) throw new RangeError('fill must be in (0, 1]');
  const ps = nu === undefined ? 1 : 1 / (1 - nu * nu);
  return { EI: fill * E * thickness ** 3 / 12 * ps };
}

/**
 * Elastic critical pressure of a circular ring under uniform external (fluid) pressure: p_cr = 3 EI / R^3.
 * Reference: Levy (1884); Timoshenko & Gere, Theory of Elastic Stability, Ch. 7 (buckling of a circular ring).
 * Validity: thin ring, pressure remaining normal to the deformed ring, mode n = 2.
 * @param {object} o
 * @param {number} o.EI ring bending stiffness per unit length, N m
 * @param {number} o.R mean ring radius, m
 * @returns {{pCr:number}} Pa
 */
export function ringElasticCollapse({ EI, R }) {
  pos(EI, 'EI'); pos(R, 'R');
  return { pCr: 3 * EI / R ** 3 };
}

/**
 * Combined elastic collapse resistance of concentric rings that ovalise together (carcass plus pressure armour
 * with a dry annulus): the rings share the n = 2 mode, so their critical pressures add,
 *     p_cr = sum_i 3 EI_i / R_i^3.
 * Reference: API RP 17B §8.4 (sum of the layer ring stiffnesses); Timoshenko & Gere §7.
 * Validity: layers remain in contact and do not transfer shear (unbonded); pressure referred to each ring's own
 * radius (radius differences between the layers neglected in the load transfer).
 * @param {object} o
 * @param {{EI:number, R:number}[]} o.rings ring stiffness (N m) and radius (m) of each supporting layer
 * @returns {{pCr:number, contributions:number[]}} Pa
 */
export function combinedElasticCollapse({ rings }) {
  if (!Array.isArray(rings) || rings.length === 0) throw new RangeError('rings must be a non-empty array');
  const contributions = rings.map((r, i) => 3 * pos(r.EI, `rings[${i}].EI`) / pos(r.R, `rings[${i}].R`) ** 3);
  return { pCr: contributions.reduce((a, c) => a + c, 0), contributions };
}

/**
 * Collapse pressure of an initially oval ring at first yield (Timoshenko).
 * Derivation: a ring with initial radial deviation w0 cos(2 theta) under pressure p develops the bending moment
 * M = p R w0 / (1 - p/p_cr) (amplified ovality). First yield at the extreme fibre:
 *     p R / a + p R w0 / (Z (1 - p/p_cr)) = sigma_y
 * with a the membrane area per unit length and Z the section modulus per unit length (solid ring: a = t,
 * Z = t^2/6). With p_y = sigma_y a / R and g = a w0 / Z (solid ring: 6 w0/t) this is the quadratic
 *     p^2 - [p_y + (1 + g) p_cr] p + p_y p_cr = 0,
 * whose smaller root is the collapse pressure. As w0 -> 0 the roots are p_y and p_cr, i.e. min(yield, elastic).
 * Ovality is defined as (Dmax - Dmin)/(Dmax + Dmin) (API Spec 17J), so w0 = ovality x R.
 * Reference: Timoshenko & Gere, Theory of Elastic Stability, Ch. 7 (tube with initial ellipticity under external
 * pressure).
 * Validity: thin ring, elastic up to first yield (conservative for strain-hardening carcass strip).
 * @param {object} o
 * @param {number} o.pCr elastic critical pressure, Pa
 * @param {number} o.yieldStress Pa
 * @param {number} o.R mean radius, m
 * @param {number} o.thickness membrane thickness a (area per unit length), m
 * @param {number} [o.sectionModulus] Z per unit length, m^2 (default thickness^2/6)
 * @param {number} [o.ovality=0.002] (Dmax - Dmin)/(Dmax + Dmin), - (API RP 17B suggests 0.2 % as a minimum when no
 *   data exist - confirm for the project)
 * @returns {{pCollapse:number, pYield:number, pCr:number, w0:number}} Pa, Pa, Pa, m
 */
export function ovalisedRingCollapse({ pCr, yieldStress, R, thickness, sectionModulus, ovality = 0.002 }) {
  pos(pCr, 'pCr'); pos(yieldStress, 'yieldStress'); pos(R, 'R'); pos(thickness, 'thickness'); nonneg(ovality, 'ovality');
  const Z = sectionModulus === undefined ? thickness * thickness / 6 : pos(sectionModulus, 'sectionModulus');
  const w0 = ovality * R;
  const pY = yieldStress * thickness / R;
  const g = thickness * w0 / Z;
  const b = pY + (1 + g) * pCr;
  // numerically stable smaller root of p^2 - b p + pY pCr = 0
  const pCollapse = 2 * pY * pCr / (b + Math.sqrt(Math.max(b * b - 4 * pY * pCr, 0)));
  return { pCollapse, pYield: pY, pCr, w0 };
}

/**
 * Collapse assessment of a rough-bore flexible pipe.
 *  - 'flooded' (wet collapse, breached outer sheath): the pressure acts on the pressure sheath and is resisted by
 *    the carcass alone.
 *  - 'dry': the pressure on the outer sheath is transmitted through the armour layers; carcass and pressure
 *    armour ovalise together and their elastic resistances add ({@link combinedElasticCollapse}). The yield term
 *    conservatively uses the carcass section only.
 * The ovalised-ring formula ({@link ovalisedRingCollapse}) is applied to the elastic pressure, then multiplied by
 * the curvature knock-down factor (bent pipe: argument, 1 = straight). Collapse depth = p/(rho g); utilisation =
 * design pressure/(factor x collapse pressure) with the API 17J carcass factor as argument.
 * Reference: API RP 17B §8.4; API Spec 17J Table 6.
 * @param {object} o
 * @param {{EI:number, R:number, thickness:number, sectionModulus?:number, yieldStress:number}} o.carcass ring
 *   stiffness N m, radius m, membrane thickness m, section modulus m^2, yield Pa
 * @param {{EI:number, R:number}} [o.pressureArmour] supporting layer for the dry case
 * @param {'dry'|'flooded'} [o.annulus='flooded']
 * @param {number} [o.ovality=0.002]
 * @param {number} [o.curvatureFactor=1] knock-down for bending, (0, 1]
 * @param {number} [o.pDesign] external design pressure, Pa (for utilisation)
 * @param {number} [o.factor=0.67] permissible utilisation (see flexible.js carcassBucklingFactor)
 * @param {number} [o.rho=1025] seawater density, kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @returns {{pElastic:number, pCollapse:number, collapseDepth:number, allowableDepth:number,
 *   utilisation:number|null, mode:string}} Pa, Pa, m, m, -
 */
export function collapseAssessment({ carcass, pressureArmour, annulus = 'flooded', ovality = 0.002, curvatureFactor = 1, pDesign, factor = 0.67, rho = 1025, g = 9.80665 }) {
  if (!carcass) throw new RangeError('carcass is required');
  if (annulus !== 'dry' && annulus !== 'flooded') throw new RangeError("annulus must be 'dry' or 'flooded'");
  if (!(curvatureFactor > 0 && curvatureFactor <= 1)) throw new RangeError('curvatureFactor must be in (0, 1]');
  pos(factor, 'factor'); pos(rho, 'rho'); pos(g, 'g');
  const rings = [{ EI: carcass.EI, R: carcass.R }];
  if (annulus === 'dry' && pressureArmour) rings.push(pressureArmour);
  const pElastic = combinedElasticCollapse({ rings }).pCr;
  const ov = ovalisedRingCollapse({ pCr: pElastic, yieldStress: carcass.yieldStress, R: carcass.R, thickness: carcass.thickness, sectionModulus: carcass.sectionModulus, ovality });
  const pCollapse = ov.pCollapse * curvatureFactor;
  return {
    pElastic, pCollapse,
    collapseDepth: pCollapse / (rho * g),
    allowableDepth: factor * pCollapse / (rho * g),
    utilisation: pDesign === undefined ? null : nonneg(pDesign, 'pDesign') / (factor * pCollapse),
    mode: annulus === 'flooded' ? 'wet collapse (carcass alone)' : 'dry collapse (carcass + pressure armour)',
  };
}

/**
 * Radial instability (bird-caging) screening of the tensile armour under axial compression.
 * Derivation: a helical layer carrying a compressive axial force C (N) has wire force Fw = C/(n cos(alpha)); the
 * wire curvature sin^2(alpha)/R turns it into an outward smeared pressure (helix.js contact-pressure formula with
 * the sign reversed)
 *     q = n Fw sin^2(alpha)/(2 pi R^2 cos(alpha)) = C tan^2(alpha) / (2 pi R^2).
 * Inner layers bear on the outer one, so the pressure demand at the outermost armour radius Ro is
 * sum_i q_i R_i / Ro. It is resisted by the hoop strength of the high-strength tape, sigma_tape t_tape / R_tape,
 * by the outer sheath, sigma_sheath t_sheath / R_sheath, and - only while the annulus is dry - by the net external
 * pressure on the outer sheath. With a flooded annulus the external pressure gives no support.
 * The compression is shared between the layers in proportion to their axial stiffness weights.
 * Reference: API RP 17B §8 (bird-caging failure mode); statics as derived above.
 * Validity: strength (not stability) criterion for the restraining layers; ignores friction and wire bending
 * stiffness, so it is conservative for short compressed lengths.
 * @param {object} o
 * @param {number} o.compression total compressive force in the tensile armour (true wall compression), N (>= 0)
 * @param {{R:number, alpha:number, weight?:number}[]} o.layers armour layers inner to outer: radius m, lay angle
 *   rad, relative axial stiffness weight (default n E A cos^3 equal -> 1)
 * @param {{strength:number, thickness:number, R:number}} [o.tape] hoop strength Pa, total thickness m, radius m
 * @param {{strength:number, thickness:number, R:number}} [o.sheath] outer sheath yield Pa, thickness m, radius m
 * @param {number} [o.pExt=0] external pressure, Pa
 * @param {number} [o.pAnnulus=0] annulus pressure, Pa
 * @param {boolean} [o.annulusFlooded=false]
 * @returns {{demandPressure:number, capacityPressure:number, utilisation:number, criticalCompression:number,
 *   tapeCapacity:number, sheathCapacity:number, pressureSupport:number}} Pa, Pa, -, N, Pa, Pa, Pa
 */
export function birdcagingCheck({ compression, layers, tape, sheath, pExt = 0, pAnnulus = 0, annulusFlooded = false }) {
  nonneg(compression, 'compression'); nonneg(pExt, 'pExt'); nonneg(pAnnulus, 'pAnnulus');
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  const wSum = layers.reduce((a, L) => a + (L.weight ?? 1), 0);
  const Ro = layers[layers.length - 1].R;
  let perUnit = 0;
  layers.forEach((L, i) => {
    pos(L.R, `layers[${i}].R`); angle(L.alpha, `layers[${i}].alpha`);
    perUnit += (L.weight ?? 1) / wSum * Math.tan(L.alpha) ** 2 / (2 * Math.PI * L.R * L.R) * L.R / Ro;
  });
  const tapeCapacity = tape ? pos(tape.strength, 'tape.strength') * pos(tape.thickness, 'tape.thickness') / pos(tape.R, 'tape.R') : 0;
  const sheathCapacity = sheath ? pos(sheath.strength, 'sheath.strength') * pos(sheath.thickness, 'sheath.thickness') / pos(sheath.R, 'sheath.R') : 0;
  const pressureSupport = annulusFlooded ? 0 : Math.max(pExt - pAnnulus, 0);
  const capacityPressure = tapeCapacity + sheathCapacity + pressureSupport;
  const demandPressure = compression * perUnit;
  return {
    demandPressure, capacityPressure,
    utilisation: capacityPressure > 0 ? demandPressure / capacityPressure : (demandPressure > 0 ? Infinity : 0),
    criticalCompression: capacityPressure / perUnit,
    tapeCapacity, sheathCapacity, pressureSupport,
  };
}

/**
 * Buckling load of a tensile-armour wire treated as a beam on an elastic (Winkler) radial foundation.
 * Derivation: E I w'''' + P w'' + k w = 0 with w = sin(pi x / l) gives P(l) = E I (pi/l)^2 + k (l/pi)^2, minimum
 *     P_cr = 2 sqrt(k E I)   at half-wavelength  l = pi (E I / k)^(1/4)
 * (Hetenyi, Beams on Elastic Foundation, Ch. VII). The foundation is the hoop stiffness of the restraining layer
 * (tape and/or sheath): a ring of modulus Er and thickness tr at radius Rr gives a pressure Er tr / Rr^2 per unit
 * radial displacement, acting on the tributary width of one wire 2 pi R cos(alpha)/n:
 *     k = Er tr / Rr^2 x 2 pi R cos(alpha) / n     (N/m per m of wire per m displacement).
 * The initial wire curvature sin^2(alpha)/R does not change this bifurcation load to first order; it produces the
 * pre-buckling outward load checked in {@link birdcagingCheck}.
 * Validity: wire long compared with the half-wavelength; elastic restraint; friction ignored. Screening only.
 * @param {object} o
 * @param {number} o.E wire modulus, Pa
 * @param {number} o.In wire second moment for bending out of the layer (b t^3/12), m^4
 * @param {number} [o.k] foundation modulus, N/m^2; if omitted built from the restraint fields below
 * @param {number} [o.restraintE] Pa
 * @param {number} [o.restraintThickness] m
 * @param {number} [o.restraintR] m
 * @param {number} [o.R] armour helix radius, m
 * @param {number} [o.alpha] lay angle, rad
 * @param {number} [o.n] wires in the layer
 * @returns {{Pcr:number, halfWavelength:number, k:number, layerCompression:number|null}} N per wire, m, N/m^2,
 *   N axial compression of the whole layer n Pcr cos(alpha) (null when n/alpha not given)
 */
export function wireFoundationBuckling({ E, In, k, restraintE, restraintThickness, restraintR, R, alpha, n }) {
  pos(E, 'E'); pos(In, 'In');
  let kf = k;
  if (kf === undefined) {
    pos(restraintE, 'restraintE'); pos(restraintThickness, 'restraintThickness'); pos(restraintR, 'restraintR');
    pos(R, 'R'); angle(alpha); pos(n, 'n');
    kf = restraintE * restraintThickness / restraintR ** 2 * 2 * Math.PI * R * Math.cos(alpha) / n;
  } else pos(kf, 'k');
  const Pcr = 2 * Math.sqrt(kf * E * In);
  return {
    Pcr, halfWavelength: Math.PI * (E * In / kf) ** 0.25, k: kf,
    layerCompression: n !== undefined && alpha !== undefined ? n * Pcr * Math.cos(alpha) : null,
  };
}

/**
 * Lateral (in-surface) buckling screening of tensile-armour wires under compression with cyclic bending.
 * Screening formula (derived, not a published capacity model): the wire is idealised as a pin-ended Euler column
 * bending about its strong axis within the cylindrical surface, radially restrained by the adjacent layers and
 * free of lateral friction (repeated bending cycles in a flooded annulus progressively remove lateral friction
 * restraint):
 *     P_cr = pi^2 E I_t / L_b^2,     L_b = lengthFactor x pitch / cos(alpha)   (wire length per pitch)
 *     C_cr = n P_cr cos(alpha)   (layer compression)
 * lengthFactor is an assumption of this screen (default 1: one buckle per pitch). Published lateral-buckling
 * tests and models (e.g. Ostergaard, Lyckegaard & Andreasen 2012, Marine Structures 27) show capacity depends on
 * curvature amplitude, number of cycles, friction and tape restraint, none of which is represented here; use the
 * result only to rank designs or to flag cases needing detailed analysis.
 * @param {object} o
 * @param {number} o.E Pa
 * @param {number} o.It wire second moment for in-surface bending (t b^3/12), m^4
 * @param {number} o.R helix radius, m
 * @param {number} o.alpha lay angle, rad
 * @param {number} o.n wires in the layer
 * @param {number} [o.lengthFactor=1] buckling length / wire length per pitch (suggested study range 0.5-2)
 * @param {number} [o.compression] layer compression for the utilisation, N
 * @returns {{Pcr:number, bucklingLength:number, layerCompression:number, utilisation:number|null}} N, m, N, -
 */
export function lateralBucklingScreen({ E, It, R, alpha, n, lengthFactor = 1, compression }) {
  pos(E, 'E'); pos(It, 'It'); pos(R, 'R'); angle(alpha); pos(n, 'n'); pos(lengthFactor, 'lengthFactor');
  const Lb = lengthFactor * 2 * Math.PI * R / Math.sin(alpha);
  const Pcr = Math.PI ** 2 * E * It / (Lb * Lb);
  const layerCompression = n * Pcr * Math.cos(alpha);
  return { Pcr, bucklingLength: Lb, layerCompression, utilisation: compression === undefined ? null : nonneg(compression, 'compression') / layerCompression };
}

/**
 * Over-bending check against unlocking of the interlocked carcass / pressure armour. Beyond the locking radius
 * the profile gaps are closed (intrados) or the hooks fully engaged (extrados); further bending loads the
 * interlock directly and can unlock the profile. Utilisation = locking radius / bend radius = kappa x rho_lock.
 * A related mechanism, upheaval of armour wires through gaps in a damaged tape or sheath, is not modelled: it
 * requires the local restraint geometry.
 * Reference: API RP 17B §8 (over-bending: unlocking of interlocked layers); locking radius from flexible.js
 * `minimumBendRadius`.
 * @param {object} o
 * @param {number} o.curvature applied curvature, 1/m
 * @param {number} o.lockingRadius m
 * @returns {{utilisation:number, locked:boolean, marginRadius:number}} -, -, m (bend radius minus locking
 *   radius; Infinity for a straight pipe)
 */
export function overbendCheck({ curvature, lockingRadius }) {
  nonneg(curvature, 'curvature'); pos(lockingRadius, 'lockingRadius');
  const utilisation = curvature * lockingRadius;
  return { utilisation, locked: utilisation >= 1, marginRadius: curvature > 0 ? 1 / curvature - lockingRadius : Infinity };
}
