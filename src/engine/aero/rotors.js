/**
 * Rotor and turbine definitions: data structure, the NREL 5-MW reference turbine planform, geometric scaling to
 * other ratings, an optimum (Betz / Glauert-Schmitz) rotor designer and tubular tower / monopile property helpers.
 *
 * Rotor data structure (plain object, SI units, degrees only where the field name ends in `Deg`):
 *   name, label                          provenance text
 *   blades                               number of blades [-]
 *   hubRadius, tipRadius                 [m] measured along the blade (pitch axis)
 *   preconeDeg, tiltDeg                  [deg] magnitudes; upwind rotor, both move the blade tips away from the tower
 *   overhang                             [m] yaw axis to hub centre along the shaft
 *   hubHeight                            [m]
 *   stations[]                           { r [m], chord [m], twistDeg [deg], thickness [-], airfoil (id), polar }
 *   polars                               { id: polar } (see polar.js)
 *   ratedPower                           [W] electrical
 *   ratedWind, cutIn, cutOut             [m/s]
 *   minRotorSpeed, ratedRotorSpeed, maxRotorSpeed   [rad/s]
 *   gearboxRatio                         [-] generator speed / rotor speed
 *   drivetrainEfficiency                 [-] electrical power / shaft power
 *   rotorInertia, generatorInertia       [kg m2] rotor about the shaft; generator about the high-speed shaft
 *   drivetrainInertia                    [kg m2] total referred to the low-speed shaft
 *   rotorMass, nacelleMass, rnaMass      [kg]
 *   tower                                { height (base to top) [m], baseElevation [m], baseDiameter, topDiameter,
 *                                          baseThickness, topThickness [m], density [kg/m3], E [Pa], mass [kg] }
 */

import { genericPolar, cylinderPolar } from './polar.js';

const D2R = Math.PI / 180;
const RPM = Math.PI / 30;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/**
 * Validate a rotor definition, attach a polar to every station and add derived quantities.
 * @param {object} def rotor definition (fields listed in the module header); `stations` ascending in r with
 *   hubRadius < r < tipRadius; each station needs `polar` or an `airfoil` id present in `def.polars`
 * @returns {object} the rotor: a shallow copy of `def` with `stations[i].polar` resolved, `rotorRadius` =
 *   tipRadius cos(precone) [m], `sweptArea` [m2] and `solidity` (blade planform area / swept area) [-]
 */
export function createRotor(def) {
  need(def && Number.isInteger(def.blades) && def.blades >= 1, 'blades', 'integer >= 1 required');
  need(def.tipRadius > 0, 'tipRadius', 'must be > 0');
  need(def.hubRadius >= 0 && def.hubRadius < def.tipRadius, 'hubRadius', 'must be in [0, tipRadius)');
  need(Array.isArray(def.stations) && def.stations.length >= 2, 'stations', 'at least two stations required');
  const polars = def.polars || {};
  let prev = def.hubRadius;
  let planform = 0;
  const stations = def.stations.map((s, i) => {
    need(s.r > prev - 1e-12 && s.r < def.tipRadius + 1e-12, `stations[${i}].r`, 'must ascend between hub and tip radius');
    need(s.chord > 0, `stations[${i}].chord`, 'must be > 0');
    const polar = s.polar || polars[s.airfoil];
    need(polar && typeof polar.clAt === 'function', `stations[${i}]`, `no polar for airfoil "${s.airfoil}"`);
    prev = s.r;
    return { r: s.r, chord: s.chord, twistDeg: s.twistDeg ?? 0, thickness: s.thickness ?? polar.thickness, airfoil: s.airfoil ?? polar.name, polar };
  });
  for (let i = 0; i <= stations.length; i++) {
    const r0 = i === 0 ? def.hubRadius : stations[i - 1].r;
    const r1 = i === stations.length ? def.tipRadius : stations[i].r;
    const c0 = i === 0 ? stations[0].chord : stations[i - 1].chord;
    const c1 = i === stations.length ? stations[i - 1].chord : stations[i].chord;
    planform += 0.5 * (c0 + c1) * (r1 - r0);
  }
  const cone = (def.preconeDeg ?? 0) * D2R;
  const rotorRadius = def.tipRadius * Math.cos(cone);
  const sweptArea = Math.PI * rotorRadius * rotorRadius;
  return {
    ...def, preconeDeg: def.preconeDeg ?? 0, tiltDeg: def.tiltDeg ?? 0, stations, polars,
    rotorRadius, sweptArea, solidity: def.blades * planform / sweptArea,
  };
}

/**
 * Section properties of a tapered circular steel tube (tower, monopile, transition piece) with linearly varying
 * outer diameter and wall thickness. A = pi (D - t) t,  I = pi/64 (D^4 - (D - 2t)^4).
 * @param {object} o
 * @param {number} o.zBase elevation of the lower end [m]
 * @param {number} o.zTop elevation of the upper end [m]
 * @param {number} o.dBase outer diameter at the lower end [m]
 * @param {number} o.dTop outer diameter at the upper end [m]
 * @param {number} o.tBase wall thickness at the lower end [m]
 * @param {number} o.tTop wall thickness at the upper end [m]
 * @param {number} [o.n=21] number of stations (>= 2)
 * @param {number} [o.density=7850] effective density incl. secondary steel [kg/m3]
 * @param {number} [o.E=210e9] Young's modulus [Pa]
 * @returns {{z:number[],D:number[],t:number[],area:number[],I:number[],massPerLength:number[],EI:number[],
 *   mass:number,cogZ:number}} stations: elevation [m], diameter [m], thickness [m], steel area [m2], second moment of
 *   area [m4], mass per length [kg/m], bending stiffness [N m2]; total mass [kg] (Simpson-exact for the linear taper)
 *   and elevation of the centre of mass [m]
 */
export function tubularDistribution(o) {
  need(o && o.zTop > o.zBase, 'zTop', 'must exceed zBase');
  need(o.dBase > 0 && o.dTop > 0, 'dBase/dTop', 'must be > 0');
  need(o.tBase > 0 && o.tTop > 0 && 2 * o.tBase < o.dBase && 2 * o.tTop < o.dTop, 'tBase/tTop', 'must be > 0 and < D/2');
  const n = o.n ?? 21;
  need(Number.isInteger(n) && n >= 2, 'n', 'integer >= 2 required');
  const rho = o.density ?? 7850;
  const E = o.E ?? 210e9;
  const out = { z: [], D: [], t: [], area: [], I: [], massPerLength: [], EI: [], mass: 0, cogZ: 0 };
  const sec = (s) => {
    const D = o.dBase + (o.dTop - o.dBase) * s;
    const t = o.tBase + (o.tTop - o.tBase) * s;
    return { D, t, A: Math.PI * (D - t) * t, I: Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4) };
  };
  for (let i = 0; i < n; i++) {
    const s = i / (n - 1);
    const q = sec(s);
    out.z.push(o.zBase + (o.zTop - o.zBase) * s);
    out.D.push(q.D); out.t.push(q.t); out.area.push(q.A); out.I.push(q.I);
    out.massPerLength.push(rho * q.A); out.EI.push(E * q.I);
  }
  // A(s) is quadratic and s A(s) cubic in s: Simpson's rule is exact for both
  const L = o.zTop - o.zBase;
  const a0 = sec(0).A; const a1 = sec(0.5).A; const a2 = sec(1).A;
  out.mass = rho * L * (a0 + 4 * a1 + a2) / 6;
  out.cogZ = o.zBase + L * (4 * 0.5 * a1 + a2) / (a0 + 4 * a1 + a2);
  return out;
}

const NREL5_R = [2.8667, 5.6, 8.3333, 11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333];
const NREL5_CHORD = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419];
const NREL5_TWIST = [13.308, 13.308, 13.308, 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106];
const NREL5_AF = ['Cylinder1', 'Cylinder1', 'Cylinder2', 'DU40', 'DU35', 'DU35', 'DU30', 'DU25', 'DU25', 'DU21', 'DU21',
  'NACA64', 'NACA64', 'NACA64', 'NACA64', 'NACA64', 'NACA64'];
const NREL5_THICK = { Cylinder1: 1, Cylinder2: 1, DU40: 0.40, DU35: 0.35, DU30: 0.30, DU25: 0.25, DU21: 0.21, NACA64: 0.18 };

/**
 * Polar set for the NREL 5-MW section ids built from GENERIC parametric polars keyed by relative thickness
 * (DU40 40 %, DU35 35 %, DU30 30 %, DU25 25 %, DU21 21 %, NACA64-618 18 %) and constant-drag cylinders
 * (Cylinder1 Cd = 0.50, Cylinder2 Cd = 0.35 as in the reference definition). These are NOT the published DU/NACA
 * tables — pass measured polars through `nrel5mw({ polars })` for design work.
 * @param {object} [o]
 * @param {number} [o.aspectRatio=17] blade aspect ratio for the Viterna extrapolation [-]
 * @param {object} [o.overrides] per-id overrides of genericPolar arguments, e.g. { DU21: { clMax: 1.3 } }
 * @returns {Object<string, object>} polars by section id
 */
export function nrel5mwGenericPolars(o = {}) {
  const ar = o.aspectRatio ?? 17;
  const ov = o.overrides || {};
  const out = {
    Cylinder1: cylinderPolar({ cd: 0.5, name: 'Cylinder1' }),
    Cylinder2: cylinderPolar({ cd: 0.35, name: 'Cylinder2' }),
  };
  for (const id of ['DU40', 'DU35', 'DU30', 'DU25', 'DU21', 'NACA64']) {
    out[id] = genericPolar({ thickness: NREL5_THICK[id], aspectRatio: ar, name: `generic-${Math.round(NREL5_THICK[id] * 100)} (${id} position)`, ...(ov[id] || {}) });
  }
  return out;
}

/**
 * NREL offshore 5-MW baseline wind turbine: gross properties and blade aerodynamic planform at the 17 analysis
 * nodes. Reference: Jonkman, J., Butterfield, S., Musial, W. and Scott, G. (2009), "Definition of a 5-MW Reference
 * Wind Turbine for Offshore System Development", NREL/TP-500-38060 (Tables 1-1, 3-1, 6-1, 7-1).
 * The airfoil polars default to generic parametric polars by thickness (see nrel5mwGenericPolars): performance
 * computed with them is indicative only.
 * @param {object} [o]
 * @param {Object<string, object>} [o.polars] polars keyed by Cylinder1, Cylinder2, DU40, DU35, DU30, DU25, DU21, NACA64
 * @returns {object} rotor (see createRotor)
 */
export function nrel5mw(o = {}) {
  const usingGeneric = !o.polars;
  const polars = o.polars || nrel5mwGenericPolars();
  // The report quotes base/top wall thicknesses of 27/19 mm; its distributed tower properties (5590.87 kg/m and
  // EI = 614.34e9 N m2 at the base with rho = 8500 kg/m3, E = 210 GPa) and the integrated mass of 347,460 kg
  // correspond to those thicknesses increased by 30 % (35.1/24.7 mm). The effective values are used for mass and
  // stiffness; the nominal ones are kept for reference.
  const tower = { height: 87.6, baseElevation: 0, baseDiameter: 6.0, topDiameter: 3.87, baseThickness: 0.027 * 1.3, topThickness: 0.019 * 1.3,
    nominalBaseThickness: 0.027, nominalTopThickness: 0.019, density: 8500, E: 210e9, mass: 347460 };
  return createRotor({
    name: 'NREL 5-MW reference',
    label: usingGeneric
      ? 'NREL 5-MW planform and gross properties (Jonkman et al. 2009) with GENERIC parametric polars by thickness'
      : 'NREL 5-MW planform and gross properties (Jonkman et al. 2009) with user-supplied polars',
    blades: 3, hubRadius: 1.5, tipRadius: 63, preconeDeg: 2.5, tiltDeg: 5, overhang: 5, hubHeight: 90,
    stations: NREL5_R.map((r, i) => ({ r, chord: NREL5_CHORD[i], twistDeg: NREL5_TWIST[i], thickness: NREL5_THICK[NREL5_AF[i]], airfoil: NREL5_AF[i] })),
    polars,
    ratedPower: 5e6, ratedWind: 11.4, cutIn: 3, cutOut: 25,
    minRotorSpeed: 6.9 * RPM, ratedRotorSpeed: 12.1 * RPM, maxRotorSpeed: 12.1 * RPM,
    gearboxRatio: 97, drivetrainEfficiency: 0.944,
    generatorInertia: 534.116, drivetrainInertia: 43784700, rotorInertia: 43784700 - 534.116 * 97 * 97,
    rotorMass: 110000, nacelleMass: 240000, rnaMass: 350000, bladeMass: 17740, bladeLength: 61.5,
    tower,
  });
}

/**
 * Headline public parameters of other open reference turbines, recorded for comparison with the scaled rotors
 * produced by `scaleRotor` (which are NOT these machines).
 *  - IEA 15-MW: Gaertner et al. (2020), "Definition of the IEA 15-Megawatt Offshore Reference Wind Turbine",
 *    NREL/TP-5000-75698.
 *  - DTU 10-MW: Bak et al. (2013), "Description of the DTU 10 MW Reference Wind Turbine", DTU Wind Energy
 *    Report-I-0092.
 * Units: ratedPower [W], rotorDiameter, hubHeight [m], ratedWind, cutIn, cutOut [m/s], ratedRotorSpeedRpm [rpm],
 * rnaMass [kg].
 */
export const PUBLIC_REFERENCE_TURBINES = {
  'NREL-5MW': { ratedPower: 5e6, rotorDiameter: 126, hubHeight: 90, ratedWind: 11.4, cutIn: 3, cutOut: 25, ratedRotorSpeedRpm: 12.1, rnaMass: 350000 },
  'DTU-10MW': { ratedPower: 10e6, rotorDiameter: 178.3, hubHeight: 119, ratedWind: 11.4, cutIn: 4, cutOut: 25, ratedRotorSpeedRpm: 9.6 },
  'IEA-15MW': { ratedPower: 15e6, rotorDiameter: 240, hubHeight: 150, ratedWind: 10.59, cutIn: 3, cutOut: 25, ratedRotorSpeedRpm: 7.56, rnaMass: 1017000 },
};

/**
 * Geometric scaling of a reference rotor to another rating ("generic scaled" rotor — NOT the IEA/DTU reference
 * machine of that rating). Classical similarity at constant tip speed and constant specific power:
 * lengths x s with s = sqrt(P/Pref); rotor speeds x 1/s; masses x s^massExponent; mass moments of inertia
 * x s^(massExponent + 2); wall thicknesses x s. Rated, cut-in and cut-out wind speeds, polars, twist, efficiency
 * and gearbox ratio are kept. The cube law (massExponent = 3) is the geometric-similarity upper bound; built blades
 * follow exponents of roughly 2.0 - 2.5.
 * @param {object} ref reference rotor (createRotor result)
 * @param {object} o
 * @param {number} o.ratedPower target rated power [W]
 * @param {number} [o.massExponent=3] mass scaling exponent [-]
 * @param {number} [o.hubHeight] hub height [m]; default keeps the reference blade-tip clearance to the base datum
 * @param {string} [o.name]
 * @returns {object} scaled rotor with `scaleFactor` s and `scaledFrom`
 */
export function scaleRotor(ref, o) {
  need(o && o.ratedPower > 0, 'ratedPower', 'must be > 0');
  need(ref && ref.ratedPower > 0, 'ref.ratedPower', 'reference rotor needs a rated power');
  const s = Math.sqrt(o.ratedPower / ref.ratedPower);
  const me = o.massExponent ?? 3;
  const sm = Math.pow(s, me);
  const si = Math.pow(s, me + 2);
  const sc = (v, f) => (v === undefined ? undefined : v * f);
  const hubHeight = o.hubHeight ?? (ref.hubHeight !== undefined ? ref.hubHeight - ref.tipRadius + ref.tipRadius * s : undefined);
  const tower = ref.tower ? {
    ...ref.tower,
    height: hubHeight !== undefined && ref.hubHeight !== undefined ? ref.tower.height + (hubHeight - ref.hubHeight) : ref.tower.height * s,
    baseDiameter: ref.tower.baseDiameter * s, topDiameter: ref.tower.topDiameter * s,
    baseThickness: ref.tower.baseThickness * s, topThickness: ref.tower.topThickness * s, mass: undefined,
  } : undefined;
  if (tower) {
    tower.mass = tubularDistribution({ zBase: 0, zTop: tower.height, dBase: tower.baseDiameter, dTop: tower.topDiameter,
      tBase: tower.baseThickness, tTop: tower.topThickness, density: tower.density, n: 3 }).mass;
  }
  const mw = (o.ratedPower / 1e6).toFixed(1).replace(/\.0$/, '');
  return createRotor({
    ...ref,
    name: o.name || `Generic ${mw}-MW (scaled from ${ref.name})`,
    label: `generic scaled rotor: ${ref.name} scaled geometrically by ${s.toFixed(4)} at constant tip speed — not a published reference turbine. Base: ${ref.label}`,
    scaledFrom: ref.name, scaleFactor: s,
    hubRadius: ref.hubRadius * s, tipRadius: ref.tipRadius * s, overhang: sc(ref.overhang, s), hubHeight,
    stations: ref.stations.map((st) => ({ ...st, r: st.r * s, chord: st.chord * s })),
    ratedPower: o.ratedPower,
    minRotorSpeed: sc(ref.minRotorSpeed, 1 / s), ratedRotorSpeed: sc(ref.ratedRotorSpeed, 1 / s), maxRotorSpeed: sc(ref.maxRotorSpeed, 1 / s),
    rotorInertia: sc(ref.rotorInertia, si), generatorInertia: sc(ref.generatorInertia, si), drivetrainInertia: sc(ref.drivetrainInertia, si),
    rotorMass: sc(ref.rotorMass, sm), nacelleMass: sc(ref.nacelleMass, sm), rnaMass: sc(ref.rnaMass, sm),
    bladeMass: sc(ref.bladeMass, sm), bladeLength: sc(ref.bladeLength, s),
    tower,
  });
}

/**
 * Local optimum blade geometry for an ideal rotor at local speed ratio lambda_r = Omega r / V.
 *  - with wake rotation (Glauert optimum, Schmitz form): phi = (2/3) atan(1/lambda_r),
 *    c = 8 pi r (1 - cos phi) / (B Cl);
 *  - without wake rotation (Betz optimum): phi = atan(2/(3 lambda_r)), c = 8 pi r sin(phi) / (3 B Cl lambda_r).
 * Reference: Manwell, McGowan and Rogers, Wind Energy Explained, 2nd ed., sections 3.7 and 3.9; Gasch and Twele,
 * Wind Power Plants, ch. 5 (Schmitz).
 * @param {object} o
 * @param {number} o.r local radius [m]
 * @param {number} o.lambdaR local speed ratio [-] (> 0)
 * @param {number} o.blades number of blades [-]
 * @param {number} o.clDesign design lift coefficient [-]
 * @param {boolean} [o.wakeRotation=true]
 * @returns {{phi:number, chord:number}} inflow angle [rad], chord [m]
 */
export function optimalSection(o) {
  need(o && o.lambdaR > 0, 'lambdaR', 'must be > 0');
  need(o.clDesign > 0, 'clDesign', 'must be > 0');
  need(o.blades >= 1, 'blades', 'must be >= 1');
  if (o.wakeRotation ?? true) {
    const phi = (2 / 3) * Math.atan(1 / o.lambdaR);
    return { phi, chord: 8 * Math.PI * o.r * (1 - Math.cos(phi)) / (o.blades * o.clDesign) };
  }
  const phi = Math.atan(2 / (3 * o.lambdaR));
  return { phi, chord: 8 * Math.PI * o.r * Math.sin(phi) / (3 * o.blades * o.clDesign * o.lambdaR) };
}

/**
 * Design an optimum rotor (chord and twist distributions) for a design tip-speed ratio and design lift coefficient;
 * twist = phi - alphaDesign at zero pitch. Suitable as a starting planform for tidal and custom rotors (no tip-loss
 * or structural adjustment of the planform is made; `maxChord` optionally clips the root).
 * @param {object} o
 * @param {number} o.tipRadius [m]
 * @param {number} [o.hubRadius=0.05*tipRadius] [m]
 * @param {number} o.blades [-]
 * @param {number} o.tsr design tip-speed ratio [-]
 * @param {number} o.clDesign design lift coefficient [-]
 * @param {number} o.alphaDesignDeg angle of attack giving clDesign on the chosen polar [deg]
 * @param {object} o.polar polar used at every station (createPolar result, +-180 deg)
 * @param {number} [o.nStations=20] number of equal-width annuli (stations at annulus mid-radius)
 * @param {boolean} [o.wakeRotation=true] Glauert/Schmitz (true) or Betz (false) optimum
 * @param {number} [o.maxChord] chord limit [m]
 * @param {string} [o.name]
 * @returns {object} rotor (createRotor) with `designTsr`
 */
export function designOptimalRotor(o) {
  need(o && o.tipRadius > 0, 'tipRadius', 'must be > 0');
  need(o.tsr > 0, 'tsr', 'must be > 0');
  need(o.polar && typeof o.polar.clAt === 'function', 'polar', 'polar object required');
  const hub = o.hubRadius ?? 0.05 * o.tipRadius;
  const n = o.nStations ?? 20;
  need(Number.isInteger(n) && n >= 2, 'nStations', 'integer >= 2 required');
  const stations = [];
  for (let i = 0; i < n; i++) {
    const r = hub + (o.tipRadius - hub) * (i + 0.5) / n;
    const sec = optimalSection({ r, lambdaR: o.tsr * r / o.tipRadius, blades: o.blades, clDesign: o.clDesign, wakeRotation: o.wakeRotation });
    stations.push({
      r, chord: o.maxChord ? Math.min(o.maxChord, sec.chord) : sec.chord,
      twistDeg: sec.phi / D2R - o.alphaDesignDeg, thickness: o.polar.thickness, airfoil: o.polar.name, polar: o.polar,
    });
  }
  return createRotor({
    name: o.name || `Optimum rotor (TSR ${o.tsr}, ${o.blades} blades)`,
    label: `ideal ${(o.wakeRotation ?? true) ? 'Glauert/Schmitz' : 'Betz'} optimum planform for TSR ${o.tsr}, Cl ${o.clDesign}; polar: ${o.polar.label}`,
    blades: o.blades, hubRadius: hub, tipRadius: o.tipRadius, preconeDeg: 0, tiltDeg: 0, stations, designTsr: o.tsr,
  });
}

/**
 * Tower mass and stiffness distribution of a rotor's tower definition (wrapper of tubularDistribution).
 * @param {object} rotor rotor with a `tower` record
 * @param {number} [n=21] number of stations
 * @returns {object} see tubularDistribution; elevations start at tower.baseElevation
 */
export function towerDistribution(rotor, n = 21) {
  const t = rotor && rotor.tower;
  need(t && t.height > 0, 'rotor.tower', 'tower definition required');
  const z0 = t.baseElevation ?? 0;
  return tubularDistribution({ zBase: z0, zTop: z0 + t.height, dBase: t.baseDiameter, dTop: t.topDiameter,
    tBase: t.baseThickness, tTop: t.topThickness, n, density: t.density ?? 7850, E: t.E ?? 210e9 });
}
