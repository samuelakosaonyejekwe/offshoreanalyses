/**
 * Submarine power and telecom cable mechanics: cross-section builder, handling limits, laying tension,
 * on-bottom stability, thermo-mechanical thrust, and telecom laying kinematics.
 *
 * Helical-layer stiffness terms are the standard forms (see helix.js for the derivation):
 * axial n E A cos^3(a), torsion n E A R^2 sin^2(a) cos(a), coupling lay n E A R sin(a) cos^2(a).
 */

const DEG = Math.PI / 180;

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
 * GENERIC handbook material constants for cable construction (room temperature, round numbers from standard
 * references such as Kaye & Laby and polymer supplier ranges). Indicative only: replace with project data.
 * E in Pa, rho in kg/m^3, alphaT in 1/K.
 */
export const GENERIC_CABLE_MATERIALS = Object.freeze({
  copper: Object.freeze({ E: 1.15e11, rho: 8900, alphaT: 1.7e-5, generic: true }),
  aluminium: Object.freeze({ E: 6.9e10, rho: 2700, alphaT: 2.3e-5, generic: true }),
  lead: Object.freeze({ E: 1.6e10, rho: 11340, alphaT: 2.9e-5, generic: true }),
  steel: Object.freeze({ E: 2.0e11, rho: 7850, alphaT: 1.2e-5, generic: true }),
  xlpe: Object.freeze({ E: 2.0e8, rho: 920, alphaT: 2.0e-4, generic: true }),
  epr: Object.freeze({ E: 2.0e7, rho: 1300, alphaT: 2.0e-4, generic: true }),
  massImpregnated: Object.freeze({ E: 5.0e8, rho: 1100, alphaT: 1.0e-4, generic: true }),
  pe: Object.freeze({ E: 8.0e8, rho: 950, alphaT: 2.0e-4, generic: true }),
  polypropyleneYarn: Object.freeze({ E: 5.0e7, rho: 600, alphaT: 1.0e-4, generic: true }),
});

/**
 * GENERIC repeatered/repeaterless telecom cable types (order-of-magnitude properties of the usual 17 mm
 * lightweight family, rounded; NOT manufacturer data - replace with the datasheet of the actual cable):
 * LW lightweight, LWP lightweight protected, SA single armoured, DA double armoured.
 * outerDiameter m; massAir kg/m; weightWater N/m; ntts, nots, npts N (nominal transient / operating / permanent
 * tensile strengths).
 */
export const GENERIC_TELECOM_CABLES = Object.freeze({
  LW: Object.freeze({ label: 'GENERIC lightweight (indicative)', generic: true, outerDiameter: 0.017, massAir: 0.55, weightWater: 3.0, ntts: 70e3, nots: 50e3, npts: 30e3 }),
  LWP: Object.freeze({ label: 'GENERIC lightweight protected (indicative)', generic: true, outerDiameter: 0.023, massAir: 0.85, weightWater: 4.2, ntts: 70e3, nots: 50e3, npts: 30e3 }),
  SA: Object.freeze({ label: 'GENERIC single armoured (indicative)', generic: true, outerDiameter: 0.028, massAir: 1.7, weightWater: 11, ntts: 220e3, nots: 150e3, npts: 100e3 }),
  DA: Object.freeze({ label: 'GENERIC double armoured (indicative)', generic: true, outerDiameter: 0.038, massAir: 3.6, weightWater: 25, ntts: 400e3, nots: 280e3, npts: 180e3 }),
});

function buildConcentric(layers, r0, tag) {
  // returns per-layer records and the outer radius
  let r = r0;
  const out = [];
  layers.forEach((L, i) => {
    const name = L.name ?? L.type ?? `${tag}[${i}]`;
    const E = nonneg(L.E ?? 0, `${tag}[${i}].E`), rho = nonneg(L.rho ?? 0, `${tag}[${i}].rho`);
    if (L.type === 'conductor') {
      if (r !== 0) throw new RangeError(`${tag}[${i}]: a conductor must be the innermost layer`);
      const A = pos(L.area, `${tag}[${i}].area`);
      const ff = L.fillFactor ?? 1;
      if (!(ff > 0 && ff <= 1)) throw new RangeError(`${tag}[${i}].fillFactor must be in (0, 1]`);
      const ro = Math.sqrt(A / (Math.PI * ff));
      const Isolid = Math.PI * ro ** 4 / 4;
      let EImin = 0;
      if (L.strandDiameter) { const ds = pos(L.strandDiameter, `${tag}[${i}].strandDiameter`); EImin = E * (A / (Math.PI * ds * ds / 4)) * Math.PI * ds ** 4 / 64; }
      out.push({ name, type: 'conductor', ri: 0, ro, area: A, mass: rho * A, EA: E * A, GJ: 0, EImax: E * Isolid * ff, EImin, K12: 0, E, material: L.material });
      r = ro;
    } else if (L.type === 'armour') {
      if (!(L.alphaDeg > 0 && L.alphaDeg < 90)) throw new RangeError(`${tag}[${i}].alphaDeg must be in (0, 90)`);
      const a = L.alphaDeg * DEG, c = Math.cos(a), s = Math.sin(a), lay = L.lay ?? 1;
      if (lay !== 1 && lay !== -1) throw new RangeError(`${tag}[${i}].lay must be +1 or -1`);
      let t, Aw, width, Iw;
      if (L.wireDiameter !== undefined) { t = pos(L.wireDiameter, `${tag}[${i}].wireDiameter`); Aw = Math.PI * t * t / 4; width = t; Iw = Math.PI * t ** 4 / 64; }
      else { t = pos(L.wireThickness, `${tag}[${i}].wireThickness`); width = pos(L.wireWidth, `${tag}[${i}].wireWidth`); Aw = t * width; Iw = width * t ** 3 / 12; }
      const R = r + t / 2;
      const nMax = 2 * Math.PI * R * c / width;
      const n = L.nWires ?? Math.floor((L.fill ?? 0.95) * nMax + 1e-9);
      if (!(n >= 1 && n <= nMax * (1 + 1e-9))) throw new RangeError(`${tag}[${i}]: wire count does not fit on the circumference`);
      const k = n * E * Aw;
      out.push({
        name, type: 'armour', ri: r, ro: r + t, R, alpha: a, lay, nWires: n, wireArea: Aw, area: n * Aw / c, mass: rho * n * Aw / c,
        EA: k * c ** 3, GJ: k * R * R * s * s * c, K12: lay * k * R * s * c * c, EImax: k * R * R * c ** 3 / 2, EImin: n * E * Iw * c, E, uts: L.uts, material: L.material,
      });
      r += t;
    } else {
      const t = pos(L.thickness, `${tag}[${i}].thickness`);
      const ri = r, ro = r + t, A = Math.PI * (ro * ro - ri * ri), I = Math.PI * (ro ** 4 - ri ** 4) / 4;
      const G = E / (2 * (1 + (L.nu ?? 0.4)));
      out.push({ name, type: L.type ?? 'solid', ri, ro, area: A, mass: rho * A, EA: E * A, GJ: G * 2 * I, K12: 0, EImax: E * I, EImin: E * I, E, material: L.material });
      r = ro;
    }
  });
  return { layers: out, ro: r };
}

/**
 * Build a submarine cable cross-section and its mechanical properties.
 *
 * Definition: { cores: 1 | 2 | 3, core: [layers of one core, conductor outwards], coreLayLength [m] (3-core /
 * bipole lay-up), fillerDensity [kg/m^3], fibreUnit { mass [kg/m] }, outer: [layers over the assembled cores] }.
 * Layer types:
 *   { type: 'conductor', area [m^2 metal], fillFactor (metal/envelope, e.g. ~0.9 compacted stranded), E, rho,
 *     strandDiameter [m] (optional, for the lower-bound bending stiffness) }
 *   { type: 'armour', wireDiameter [m] | (wireWidth, wireThickness) [m], nWires | fill (default 0.95),
 *     alphaDeg, lay (+1/-1), E, rho, uts }
 *   any other type (screen, insulation, sheath, bedding, serving ...): { thickness [m], E, rho, nu }
 * Geometry of the lay-up (exact for touching circles): three cores of diameter d sit on a pitch radius d/sqrt(3)
 * inside an envelope d (1 + 2/sqrt(3)); two cores on radius d/2 inside 2d. Core lay angle tan(a) = 2 pi Rc/L.
 *
 * Returned bending stiffness is a bracket: EImin (all strands and wires slip: each bends about its own axis) and
 * EImax (solid conductor, armour in full stick). Real cables lie between and are temperature and curvature
 * dependent (bitumen, polymers): use measured values for dynamic analysis.
 * Reference: CIGRE TB 623 (2015) mechanical properties in intent; helix formulas per helix.js.
 * @param {object} def cable definition (see above)
 * @param {object} [o]
 * @param {number} [o.rhoSeawater=1025] kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @returns {object} { outerDiameter [m], coreDiameter [m], mass [kg/m], weightAir, submergedWeight [N/m],
 *   specificGravity [-], EA (torque-free), EAfixed [N], GJ (tension-free), GJfixed [N m^2], coupling K12 [N m],
 *   torquePerTension [m], torqueBalanced, EImin, EImax [N m^2], conductorArea [m^2 total], layers: [...] }
 */
export function buildCable(def, { rhoSeawater = 1025, g = 9.80665 } = {}) {
  if (!def || !Array.isArray(def.core) || def.core.length === 0) throw new RangeError('def.core must be a non-empty array of layers');
  pos(rhoSeawater, 'rhoSeawater'); pos(g, 'g');
  const nCores = def.cores ?? 1;
  if (![1, 2, 3].includes(nCores)) throw new RangeError('cores must be 1, 2 or 3');
  const core = buildConcentric(def.core, 0, 'core');
  const dCore = 2 * core.ro;
  let mass = 0, K11 = 0, K12 = 0, K22 = 0, EImin = 0, EImax = 0, condA = 0, r0;
  const sum = (key) => core.layers.reduce((a, L) => a + L[key], 0);
  const records = [];
  if (nCores === 1) {
    mass += sum('mass'); K11 += sum('EA'); K12 += sum('K12'); K22 += sum('GJ'); EImin += sum('EImin'); EImax += sum('EImax');
    r0 = core.ro;
    for (const L of core.layers) records.push({ ...L, group: 'core' });
  } else {
    const Rc = nCores === 3 ? dCore / Math.sqrt(3) : dCore / 2;
    const env = nCores === 3 ? dCore * (1 + 2 / Math.sqrt(3)) : 2 * dCore;
    const layLength = pos(def.coreLayLength, 'coreLayLength');
    const a = Math.atan(2 * Math.PI * Rc / layLength), c = Math.cos(a), s = Math.sin(a), lay = def.coreLay ?? 1;
    const EAc = sum('EA'), EIcMin = sum('EImin'), EIcMax = sum('EImax'), GJc = sum('GJ');
    mass += nCores * sum('mass') / c;
    K11 += nCores * EAc * c ** 3; K12 += lay * nCores * EAc * Rc * s * c * c; K22 += nCores * (EAc * Rc * Rc * s * s * c + GJc * c);
    EImin += nCores * EIcMin * c; EImax += nCores * (EIcMax * c + EAc * Rc * Rc * c ** 3 / 2);
    const fillArea = Math.PI * env * env / 4 - nCores * Math.PI * dCore * dCore / 4;
    mass += nonneg(def.fillerDensity ?? 0, 'fillerDensity') * fillArea + nonneg(def.fibreUnit?.mass ?? 0, 'fibreUnit.mass');
    r0 = env / 2;
    for (const L of core.layers) records.push({ ...L, group: 'core', count: nCores, layAngle: a, pitchRadius: Rc });
  }
  for (const L of core.layers) if (L.type === 'conductor') condA += nCores * L.area;
  const outer = buildConcentric(def.outer ?? [], r0, 'outer');
  // buildConcentric starts conductors at r = 0 only; outer layers start from the assembled radius
  for (const L of outer.layers) {
    mass += L.mass; K11 += L.EA; K12 += L.K12; K22 += L.GJ; EImin += L.EImin; EImax += L.EImax;
    records.push({ ...L, group: 'outer' });
  }
  const OD = 2 * outer.ro;
  const disp = rhoSeawater * Math.PI * OD * OD / 4;
  return {
    outerDiameter: OD, coreDiameter: dCore, mass, weightAir: mass * g, submergedWeight: (mass - disp) * g, specificGravity: mass / disp,
    EA: K22 > 0 ? K11 - K12 * K12 / K22 : K11, EAfixed: K11, GJ: K22 - K12 * K12 / K11, GJfixed: K22, coupling: K12,
    torquePerTension: K12 / K11, torqueBalanced: Math.abs(K12) < 0.02 * Math.sqrt(K11 * K22),
    EImin, EImax, conductorArea: condA, cores: nCores, layers: records,
  };
}

/**
 * Maximum allowable pulling tension.
 *  - conductor-based (pulling eye / stocking on the conductors): T = sigma_c x total conductor area, with
 *    sigma_c = 50 N/mm^2 for copper and 30 N/mm^2 for aluminium (long-standing cable-pulling practice quoted in
 *    cable standards and manufacturer guidance; overridable);
 *  - armour-based (pulling by the armour): T = utilisation x sum over armour layers of n A sigma_u cos(alpha);
 *    the utilisation is a required argument because CIGRE TB 623 leaves the permissible armour stress to the
 *    tensile-bending type test of the actual cable.
 * @param {object} cable result of {@link buildCable}
 * @param {object} [o]
 * @param {'copper'|'aluminium'} [o.conductor='copper']
 * @param {number} [o.conductorStress] allowable conductor stress, Pa (default 50e6 Cu / 30e6 Al)
 * @param {number} [o.armourUtilisation] allowable fraction of armour UTS (omit to skip the armour limit)
 * @returns {{conductorLimit:number, armourLimit:number|null, armourStrength:number|null}} N
 */
export function cablePullingTension(cable, { conductor = 'copper', conductorStress, armourUtilisation } = {}) {
  if (conductor !== 'copper' && conductor !== 'aluminium') throw new RangeError("conductor must be 'copper' or 'aluminium'");
  const sc = conductorStress === undefined ? (conductor === 'copper' ? 50e6 : 30e6) : pos(conductorStress, 'conductorStress');
  let strength = null;
  for (const L of cable.layers) if (L.type === 'armour' && L.uts) strength = (strength ?? 0) + L.nWires * L.wireArea * L.uts * Math.cos(L.alpha);
  let armourLimit = null;
  if (armourUtilisation !== undefined) {
    if (!(armourUtilisation > 0 && armourUtilisation <= 1)) throw new RangeError('armourUtilisation must be in (0, 1]');
    if (strength === null) throw new RangeError('cable has no armour layer with uts');
    armourLimit = armourUtilisation * strength;
  }
  return { conductorLimit: sc * cable.conductorArea, armourLimit, armourStrength: strength };
}

/**
 * Handling limits: minimum bend radius as multiples of the overall diameter, sidewall pressure on a bend, and
 * tensioner grip.
 *  - MBR = factor x OD: installation (under tension) typically 15-20 x OD, static/installed 10-12 x OD
 *    (manufacturer-specific; defaults are the conservative ends of those ranges).
 *  - Sidewall pressure on a chute / sheave / bend: SWP = T / R (N per metre of cable), against an allowable from
 *    the cable's crush test (argument).
 *  - Tensioner: squeeze q per track per metre must not exceed the allowable; grip length L = SF T /(mu n q).
 * Reference: CIGRE TB 623 (2015) handling parameters in intent; capstan/sidewall statics.
 * @param {object} o
 * @param {number} o.outerDiameter m
 * @param {number} [o.installFactor=20] MBR / OD during installation
 * @param {number} [o.staticFactor=12] MBR / OD installed
 * @param {number} [o.tension] N
 * @param {number} [o.bendRadius] actual radius, m
 * @param {number} [o.sidewallLimit] allowable sidewall pressure, N/m
 * @param {number} [o.squeeze] tensioner squeeze per track, N/m
 * @param {number} [o.squeezeLimit] allowable squeeze, N/m
 * @param {number} [o.nTracks=2]
 * @param {number} [o.mu] pad friction coefficient
 * @param {number} [o.safetyFactor=1.5] on grip
 * @returns {{mbrInstall:number, mbrStatic:number, bendUtilisation:number|null, sidewallPressure:number|null,
 *   sidewallUtilisation:number|null, minRadiusForSidewall:number|null, squeezeUtilisation:number|null,
 *   gripLength:number|null}} m, m, -, N/m, -, m, -, m
 */
export function cableHandlingLimits({ outerDiameter, installFactor = 20, staticFactor = 12, tension, bendRadius, sidewallLimit, squeeze, squeezeLimit, nTracks = 2, mu, safetyFactor = 1.5 }) {
  pos(outerDiameter, 'outerDiameter'); pos(installFactor, 'installFactor'); pos(staticFactor, 'staticFactor');
  const out = { mbrInstall: installFactor * outerDiameter, mbrStatic: staticFactor * outerDiameter, bendUtilisation: null, sidewallPressure: null, sidewallUtilisation: null, minRadiusForSidewall: null, squeezeUtilisation: null, gripLength: null };
  if (bendRadius !== undefined) out.bendUtilisation = out.mbrInstall / pos(bendRadius, 'bendRadius');
  if (tension !== undefined) {
    nonneg(tension, 'tension');
    if (bendRadius !== undefined) out.sidewallPressure = tension / bendRadius;
    if (sidewallLimit !== undefined) {
      out.minRadiusForSidewall = tension / pos(sidewallLimit, 'sidewallLimit');
      if (out.sidewallPressure !== null) out.sidewallUtilisation = out.sidewallPressure / sidewallLimit;
    }
  }
  if (squeeze !== undefined) {
    pos(squeeze, 'squeeze');
    if (squeezeLimit !== undefined) out.squeezeUtilisation = squeeze / pos(squeezeLimit, 'squeezeLimit');
    if (tension !== undefined && mu !== undefined) out.gripLength = pos(safetyFactor, 'safetyFactor') * tension / (pos(mu, 'mu') * pos(nTracks, 'nTracks') * squeeze);
  }
  return out;
}

/**
 * Impact / crush classification against caller-supplied energy bands (e.g. from the cable's impact type test or
 * a dropped-object study). Returns the first band whose maxEnergy is not exceeded.
 * @param {object} o
 * @param {number} o.energy impact energy, J
 * @param {{maxEnergy:number, label:string}[]} o.bands ascending energy bands, J
 * @returns {{label:string, band:number, exceedsAll:boolean}}
 */
export function cableImpactClass({ energy, bands }) {
  nonneg(energy, 'energy');
  if (!Array.isArray(bands) || bands.length === 0) throw new RangeError('bands must be a non-empty array');
  for (let i = 0; i < bands.length; i++) {
    if (i > 0 && !(bands[i].maxEnergy > bands[i - 1].maxEnergy)) throw new RangeError('bands must be ascending in maxEnergy');
    if (energy <= bands[i].maxEnergy) return { label: bands[i].label, band: i, exceedsAll: false };
  }
  return { label: 'exceeds all bands', band: bands.length, exceedsAll: true };
}

/**
 * Coiling ability. Laying a cable into a fixed coil (tank or basket without a rotating turntable) imposes one
 * full turn of twist per coil turn: twist rate = 2 pi /(2 pi Rcoil) = 1/Rcoil rad/m (exact geometry). A cable
 * with armour wound in one direction accepts this twist in the lay-opening direction (the armour bird-cages
 * slightly and recovers); a torque-balanced cable with cross-wound armour layers resists twist in both directions
 * and must be stored on a turntable or reel. The armour strain induced by the twist at zero axial strain is
 * lay R sin(a) cos(a) / Rcoil (helix kinematics).
 * Reference: CIGRE TB 623 (2015) coiling test rationale; helix.js Eq. (K1).
 * @param {object} cable result of {@link buildCable}
 * @param {object} o
 * @param {number} o.coilRadius m
 * @returns {{coilable:boolean, twistRate:number, armour:{name:string, strain:number}[], note:string}} rad/m
 */
export function cableCoiling(cable, { coilRadius }) {
  pos(coilRadius, 'coilRadius');
  const arm = cable.layers.filter((L) => L.type === 'armour');
  const lays = new Set(arm.map((L) => L.lay));
  const coilable = arm.length > 0 && lays.size === 1;
  const twistRate = 1 / coilRadius;
  return {
    coilable, twistRate,
    armour: arm.map((L) => ({ name: L.name, strain: -Math.abs(L.R * Math.sin(L.alpha) * Math.cos(L.alpha) * twistRate) })),
    note: arm.length === 0 ? 'Unarmoured cable: coiling governed by the core twist tolerance.'
      : coilable ? `Unidirectional armour: coil in the lay-opening direction (${arm[0].lay > 0 ? 'against a right-hand' : 'against a left-hand'} lay); armour strain listed is the opening (compressive) strain.`
        : 'Cross-wound (torque-balanced) armour: not coilable, use a turntable or reel.',
  };
}

/**
 * Laying tension.
 * Static catenary (exact for an inextensible cable without bending stiffness or current): top tension
 * T = H + w d with H the bottom (horizontal) tension, w the submerged weight per unit length and d the water
 * depth; departure angle from the horizontal cos(theta) = H/T; suspended length s = sqrt(d (d + 2 H/w));
 * layback x = (H/w) acosh(1 + w d/H); radius at touchdown R = H/w.
 * Test/design tension in the form used for the tensile-bending test of submarine cables for laying depths up to
 * 500 m (CIGRE Electra 171 / TB 623, as recalled): T = dynamicFactor x w d + H with dynamicFactor = 1.3 covering
 * vessel-motion dynamics. For deeper water the recommendation adds an explicit dynamic term from vessel heave;
 * here any dynamic amplification beyond the factor must be supplied by a dynamic analysis.
 * @param {object} o
 * @param {number} o.submergedWeight w, N/m
 * @param {number} o.depth d, m
 * @param {number} o.bottomTension H, N (> 0)
 * @param {number} [o.dynamicFactor=1.3]
 * @returns {{staticTopTension:number, designTension:number, departureAngle:number, suspendedLength:number,
 *   layback:number, touchdownRadius:number}} N, N, rad from horizontal, m, m, m
 */
export function cableLayTension({ submergedWeight, depth, bottomTension, dynamicFactor = 1.3 }) {
  const w = pos(submergedWeight, 'submergedWeight'), d = pos(depth, 'depth'), H = pos(bottomTension, 'bottomTension');
  pos(dynamicFactor, 'dynamicFactor');
  const T = H + w * d;
  return {
    staticTopTension: T, designTension: dynamicFactor * w * d + H, departureAngle: Math.acos(H / T),
    suspendedLength: Math.sqrt(d * (d + 2 * H / w)), layback: H / w * Math.acosh(1 + w * d / H), touchdownRadius: H / w,
  };
}

/**
 * On-bottom stability of an unburied cable in steady near-bed flow and its specific gravity.
 * Lateral balance per unit length (Coulomb friction): mu (w - F_L) >= SF x F_D, with
 * F_D = 0.5 rho C_D D U^2 and F_L = 0.5 rho C_L D U^2. Required submerged weight w_req = F_L + SF F_D/mu.
 * Specific gravity = mass /(rho pi D^2/4).
 * Reference: Morison drag/lift on a cylinder on the seabed; absolute lateral static stability in the sense of
 * DNV-RP-F109 (whose calibrated peak-load coefficients should be used for design). C_D and C_L are required
 * arguments: they depend on Reynolds number, roughness, embedment and wave/current ratio (order 0.7-1.2 and
 * 0.8-1.0 for a cylinder resting on the bed in steady flow).
 * @param {object} o
 * @param {number} o.submergedWeight N/m
 * @param {number} o.outerDiameter m
 * @param {number} o.velocity near-bed velocity normal to the cable, m/s
 * @param {number} o.mu lateral friction coefficient
 * @param {number} o.Cd drag coefficient
 * @param {number} o.Cl lift coefficient
 * @param {number} [o.rho=1025] kg/m^3
 * @param {number} [o.safetyFactor=1.1]
 * @param {number} [o.mass] kg/m (for the specific gravity)
 * @returns {{drag:number, lift:number, resistance:number, utilisation:number, requiredWeight:number,
 *   criticalVelocity:number, specificGravity:number|null, stable:boolean}} N/m, N/m, N/m, -, N/m, m/s
 */
export function cableOnBottomStability({ submergedWeight, outerDiameter, velocity, mu, Cd, Cl, rho = 1025, safetyFactor = 1.1, mass }) {
  fin(submergedWeight, 'submergedWeight'); pos(outerDiameter, 'outerDiameter'); nonneg(velocity, 'velocity'); pos(mu, 'mu');
  pos(Cd, 'Cd'); nonneg(Cl, 'Cl'); pos(rho, 'rho'); pos(safetyFactor, 'safetyFactor');
  const q = 0.5 * rho * outerDiameter * velocity * velocity;
  const drag = Cd * q, lift = Cl * q;
  const resistance = mu * Math.max(submergedWeight - lift, 0);
  return {
    drag, lift, resistance,
    utilisation: resistance > 0 ? safetyFactor * drag / resistance : (drag > 0 ? Infinity : 0),
    requiredWeight: lift + safetyFactor * drag / mu,
    criticalVelocity: submergedWeight > 0 ? Math.sqrt(submergedWeight / (0.5 * rho * outerDiameter * (Cl + safetyFactor * Cd / mu))) : 0,
    specificGravity: mass === undefined ? null : pos(mass, 'mass') / (rho * Math.PI * outerDiameter ** 2 / 4),
    stable: safetyFactor * drag <= resistance,
  };
}

/**
 * Thermo-mechanical behaviour of the conductor.
 *  - Rigid (fully restrained) installation: thrust F = E_eff A alpha dT, with E_eff the effective modulus of the
 *    stranded conductor (lower than the metal modulus; argument).
 *  - Flexible installation (snaking): if the cable can move sideways the thermal excess length is absorbed in a
 *    sinusoidal snake y = a sin(2 pi x/lambda), whose arc-length excess is (pi a/lambda)^2 to second order. With a
 *    residual thrust Fres the strain to absorb is alpha dT - Fres/(E_eff A), so
 *    a = (lambda/pi) sqrt(alpha dT - Fres/(E_eff A)) when starting straight, or the increase from an initial
 *    amplitude a0: a = sqrt(a0^2 + (lambda/pi)^2 strain). Peak curvature of the snake = a (2 pi/lambda)^2.
 * Reference: thermal expansion of a restrained bar; sine-wave arc length (cable cleat/clamp design uses
 * these forces).
 * @param {object} o
 * @param {number} o.area conductor area, m^2
 * @param {number} o.Eeff effective conductor modulus, Pa
 * @param {number} o.alphaT expansion coefficient, 1/K
 * @param {number} o.deltaT conductor temperature rise, K
 * @param {number} [o.snakeWavelength] lambda, m
 * @param {number} [o.initialAmplitude=0] a0, m
 * @param {number} [o.residualThrust=0] N
 * @returns {{thrust:number, freeStrain:number, snakeAmplitude:number|null, snakeCurvature:number|null}} N, -,
 *   m, 1/m
 */
export function conductorThermalThrust({ area, Eeff, alphaT, deltaT, snakeWavelength, initialAmplitude = 0, residualThrust = 0 }) {
  pos(area, 'area'); pos(Eeff, 'Eeff'); pos(alphaT, 'alphaT'); fin(deltaT, 'deltaT'); nonneg(initialAmplitude, 'initialAmplitude'); nonneg(residualThrust, 'residualThrust');
  const freeStrain = alphaT * deltaT;
  let amp = null, curv = null;
  if (snakeWavelength !== undefined) {
    pos(snakeWavelength, 'snakeWavelength');
    const strain = Math.max(freeStrain - residualThrust / (Eeff * area), 0);
    amp = Math.sqrt(initialAmplitude ** 2 + (snakeWavelength / Math.PI) ** 2 * strain);
    curv = amp * (2 * Math.PI / snakeWavelength) ** 2;
  }
  return { thrust: Eeff * area * freeStrain, freeStrain, snakeAmplitude: amp, snakeCurvature: curv };
}

/**
 * Telecom cable tension utilisation against its rated strengths:
 * NTTS (nominal transient tensile strength: short-duration peak, e.g. recovery dynamics),
 * NOTS (nominal operating tensile strength: sustained during an operation, typically up to 48 h),
 * NPTS (nominal permanent tensile strength: residual tension for the system life).
 * Reference: ITU-T G.972 (definitions of terms for optical fibre submarine cable systems); values from the cable
 * datasheet.
 * @param {object} o
 * @param {number} [o.transientTension] N
 * @param {number} [o.operatingTension] N
 * @param {number} [o.permanentTension] N
 * @param {number} o.ntts N
 * @param {number} [o.nots] N
 * @param {number} [o.npts] N
 * @returns {{transient:number|null, operating:number|null, permanent:number|null, governing:number}}
 *   utilisations, -
 */
export function telecomTensionUtilisation({ transientTension, operatingTension, permanentTension, ntts, nots, npts }) {
  pos(ntts, 'ntts');
  const u = (T, cap, n) => (T === undefined || cap === undefined ? null : nonneg(T, n) / pos(cap, n));
  const r = { transient: u(transientTension, ntts, 'transientTension'), operating: u(operatingTension, nots, 'operatingTension'), permanent: u(permanentTension, npts, 'permanentTension') };
  return { ...r, governing: Math.max(r.transient ?? 0, r.operating ?? 0, r.permanent ?? 0) };
}

/**
 * Steady-state sinking of a cable during laying (straight-line solution of Zajac).
 * Derivation: at constant ship speed V the suspended cable is a straight line at angle alpha to the horizontal
 * when the weight component normal to the cable balances the normal drag,
 *     w cos(alpha) = 0.5 rho C_D D (V sin(alpha))^2.
 * For small angles alpha ~ H/V with the hydrodynamic constant H = sqrt(2 w /(rho C_D D)) (rad m/s; usually quoted
 * in degree-knots). The exact angle is cos(alpha) = (sqrt(1 + 4 k^2) - 1)/(2 k), k = rho C_D D V^2/(2 w).
 * The sinking speed normal to the surface is V tan(alpha) (touchdown lags the ship by depth/tan(alpha)).
 * Reference: Zajac (1957), Bell System Technical Journal 36, 1129-1207. C_D is the normal drag coefficient of
 * the cable (argument).
 * @param {object} o
 * @param {number} o.weightWater w, N/m
 * @param {number} o.outerDiameter D, m
 * @param {number} o.shipSpeed V, m/s
 * @param {number} o.Cd normal drag coefficient
 * @param {number} [o.rho=1025] kg/m^3
 * @param {number} [o.depth] m (for the touchdown lag)
 * @returns {{hydrodynamicConstant:number, hydrodynamicConstantDegKnots:number, sinkAngle:number,
 *   sinkAngleSmall:number, sinkAngleDeg:number, layback:number|null, sinkTime:number|null}} rad m/s; deg knots;
 *   rad (exact); rad (H/V); deg; m; s
 */
export function telecomSinkAngle({ weightWater, outerDiameter, shipSpeed, Cd, rho = 1025, depth }) {
  const w = pos(weightWater, 'weightWater'), D = pos(outerDiameter, 'outerDiameter'), V = pos(shipSpeed, 'shipSpeed');
  pos(Cd, 'Cd'); pos(rho, 'rho');
  const H = Math.sqrt(2 * w / (rho * Cd * D));
  const k = rho * Cd * D * V * V / (2 * w);
  const alpha = Math.acos((Math.sqrt(1 + 4 * k * k) - 1) / (2 * k));
  return {
    hydrodynamicConstant: H, hydrodynamicConstantDegKnots: H / DEG / 0.514444, sinkAngle: alpha, sinkAngleSmall: H / V, sinkAngleDeg: alpha / DEG,
    layback: depth === undefined ? null : pos(depth, 'depth') / Math.tan(alpha),
    sinkTime: depth === undefined ? null : depth / (V * Math.tan(alpha)),
  };
}

/**
 * Slack control during laying (exact straight-line geometry). The ship advances horizontally at Vs; the cable
 * leaves at sink angle alpha; the seabed descends ahead at slope angle beta (negative for a rising seabed).
 * Per unit ship advance the touchdown point travels sin(a)/sin(a + b) along the seabed and the suspended length
 * grows by sin(b)/sin(a + b), so laying with bottom slack e (cable length per seabed length = 1 + e) requires
 *     Vc/Vs = [(1 + e) sin(alpha) + sin(beta)] / sin(alpha + beta),
 * which reduces to 1 + e on a flat seabed. Conversely e = [(Vc/Vs) sin(a + b) - sin(b)]/sin(a) - 1.
 * Reference: Zajac (1957) §III (laying over a sloping bottom); derived from the triangle ship - touchdown -
 * seabed.
 * Validity: steady straight-line configuration, alpha + beta in (0, pi).
 * @param {object} o
 * @param {number} o.sinkAngle alpha, rad
 * @param {number} [o.seabedSlope=0] beta, rad (positive: seabed descending in the lay direction)
 * @param {number} o.shipSpeed Vs, m/s
 * @param {number} [o.slack] required bottom slack e (fraction) -> pay-out speed
 * @param {number} [o.payoutSpeed] Vc, m/s -> resulting slack (give exactly one of slack / payoutSpeed)
 * @returns {{payoutSpeed:number, slack:number, slackPercent:number, speedRatio:number, zeroSlackRatio:number}}
 *   m/s, -, %, Vc/Vs, Vc/Vs for e = 0
 */
export function telecomSlack({ sinkAngle, seabedSlope = 0, shipSpeed, slack, payoutSpeed }) {
  const a = sinkAngle, b = fin(seabedSlope, 'seabedSlope');
  if (!(a > 0 && a < Math.PI / 2)) throw new RangeError('sinkAngle must be in (0, pi/2)');
  if (!(a + b > 0 && a + b < Math.PI)) throw new RangeError('sinkAngle + seabedSlope must be in (0, pi)');
  pos(shipSpeed, 'shipSpeed');
  if ((slack === undefined) === (payoutSpeed === undefined)) throw new RangeError('give exactly one of slack or payoutSpeed');
  const sab = Math.sin(a + b), zero = (Math.sin(a) + Math.sin(b)) / sab;
  let e, ratio;
  if (slack !== undefined) { e = fin(slack, 'slack'); ratio = ((1 + e) * Math.sin(a) + Math.sin(b)) / sab; }
  else { ratio = pos(payoutSpeed, 'payoutSpeed') / shipSpeed; e = (ratio * sab - Math.sin(b)) / Math.sin(a) - 1; }
  return { payoutSpeed: ratio * shipSpeed, slack: e, slackPercent: 100 * e, speedRatio: ratio, zeroSlackRatio: zero };
}
