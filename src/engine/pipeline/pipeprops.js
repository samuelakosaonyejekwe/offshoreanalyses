/**
 * Pipe section properties, layered coatings, weights, buoyancy, line-pipe material grades and
 * standard pipe sizes.
 *
 * All quantities SI (m, kg, N, Pa); temperatures in °C where the argument name ends in `C`.
 */

const G = 9.80665;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Geometric properties of a circular hollow section (exact annulus formulas, no thin-wall approximation).
 *
 * Reference: elementary section geometry.
 * Validity: 0 < t < D/2.
 *
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @returns {{D:number, t:number, Di:number, As:number, Ai:number, Ae:number, I:number, J:number, Z:number, rGyr:number, DtRatio:number}}
 *   Di internal diameter [m]; As steel area [m²]; Ai bore area [m²]; Ae area enclosed by the steel outer
 *   surface [m²]; I second moment of area [m⁴]; J polar moment [m⁴]; Z elastic section modulus [m³];
 *   rGyr radius of gyration [m]; DtRatio D/t [-]
 */
export function pipeSection({ D, t }) {
  pos('D', D); pos('t', t);
  if (t >= D / 2) throw new RangeError('t must be smaller than D/2');
  const Di = D - 2 * t;
  const As = (Math.PI / 4) * (D * D - Di * Di);
  const Ai = (Math.PI / 4) * Di * Di;
  const Ae = (Math.PI / 4) * D * D;
  const I = (Math.PI / 64) * (D ** 4 - Di ** 4);
  return { D, t, Di, As, Ai, Ae, I, J: 2 * I, Z: (2 * I) / D, rGyr: Math.sqrt(I / As), DtRatio: D / t };
}

/**
 * Weights per unit length of a coated pipe: steel, concentric coating layers (anti-corrosion, insulation,
 * concrete weight coat with water absorption), contents and marine growth; buoyancy, submerged weight for
 * the empty, water-filled and operating conditions, and specific gravity.
 *
 * Water absorption of a layer is the mass of absorbed water per unit dry mass of that layer (concrete weight
 * coating is typically specified with 2–5 % by mass); it adds to the weight in the submerged conditions only.
 * Buoyancy is the weight of water displaced by the outermost surface (including marine growth if given).
 *
 * Reference: Archimedes; layer areas from annulus geometry.
 * Validity: any concentric build-up; field joints are not modelled (apply an equivalent density if needed).
 *
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t steel wall thickness [m]
 * @param {number} [o.rhoSteel=7850] steel density [kg/m³]
 * @param {Array<{name?:string, thickness:number, density:number, absorption?:number}>} [o.coatings=[]]
 *   layers from the steel outward: thickness [m], dry density [kg/m³], absorption [kg water / kg dry layer]
 * @param {number} [o.rhoContents=0] density of the operating contents [kg/m³]
 * @param {number} [o.rhoWater=1025] seawater density [kg/m³]
 * @param {number} [o.rhoFlood=o.rhoWater] density of the flooding / hydrotest water [kg/m³]
 * @param {{thickness:number, density:number}} [o.marineGrowth] marine growth layer: thickness [m], density [kg/m³]
 * @param {number} [o.g=9.80665] gravitational acceleration [m/s²]
 * @returns {{Dtotal:number, Dhydro:number, massSteel:number, massCoatings:number, massLayers:number[],
 *   massAbsorbed:number, massContents:number, massFlood:number, massMarineGrowth:number, massDisplaced:number,
 *   weightDry:number, buoyancy:number, submergedEmpty:number, submergedFlooded:number, submergedOperating:number,
 *   sgEmpty:number, sgFlooded:number, sgOperating:number, massOperating:number}}
 *   Dtotal outer diameter over coatings [m]; Dhydro hydrodynamic diameter incl. marine growth [m];
 *   mass* in kg/m; weightDry weight in air, empty, without absorbed water or marine growth [N/m];
 *   buoyancy [N/m]; submerged* [N/m] (negative = buoyant); sg* = weight in the condition / buoyancy [-];
 *   massOperating structural + contents mass, without hydrodynamic added mass [kg/m]
 */
export function pipeWeights({ D, t, rhoSteel = 7850, coatings = [], rhoContents = 0, rhoWater = 1025, rhoFlood,
  marineGrowth, g = G }) {
  const sec = pipeSection({ D, t });
  pos('rhoSteel', rhoSteel); pos('rhoWater', rhoWater); nonNeg('rhoContents', rhoContents); pos('g', g);
  const rhoFl = rhoFlood === undefined ? rhoWater : pos('rhoFlood', rhoFlood);
  const massSteel = rhoSteel * sec.As;
  let Do = D;
  let massCoatings = 0;
  let massAbsorbed = 0;
  const massLayers = [];
  coatings.forEach((c, i) => {
    nonNeg(`coatings[${i}].thickness`, c.thickness); pos(`coatings[${i}].density`, c.density);
    const absorption = c.absorption === undefined ? 0 : nonNeg(`coatings[${i}].absorption`, c.absorption);
    const Dn = Do + 2 * c.thickness;
    const m = c.density * (Math.PI / 4) * (Dn * Dn - Do * Do);
    massLayers.push(m);
    massCoatings += m;
    massAbsorbed += m * absorption;
    Do = Dn;
  });
  const Dtotal = Do;
  let massMarineGrowth = 0;
  let Dhydro = Dtotal;
  if (marineGrowth && marineGrowth.thickness > 0) {
    pos('marineGrowth.density', marineGrowth.density);
    Dhydro = Dtotal + 2 * marineGrowth.thickness;
    massMarineGrowth = marineGrowth.density * (Math.PI / 4) * (Dhydro * Dhydro - Dtotal * Dtotal);
  }
  const massContents = rhoContents * sec.Ai;
  const massFlood = rhoFl * sec.Ai;
  const massDisplaced = rhoWater * (Math.PI / 4) * Dhydro * Dhydro;
  const buoyancy = massDisplaced * g;
  const wetBase = massSteel + massCoatings + massAbsorbed + massMarineGrowth;
  const weightDry = (massSteel + massCoatings) * g;
  const wEmpty = wetBase * g;
  const wFlood = (wetBase + massFlood) * g;
  const wOper = (wetBase + massContents) * g;
  return {
    Dtotal, Dhydro, massSteel, massCoatings, massLayers, massAbsorbed, massContents, massFlood, massMarineGrowth,
    massDisplaced, weightDry, buoyancy,
    submergedEmpty: wEmpty - buoyancy, submergedFlooded: wFlood - buoyancy, submergedOperating: wOper - buoyancy,
    sgEmpty: wEmpty / buoyancy, sgFlooded: wFlood / buoyancy, sgOperating: wOper / buoyancy,
    massOperating: wetBase + massContents,
  };
}

/** Specified minimum yield / tensile strength [Pa] by grade. API 5L (46th ed.) PSL 2 Table 7 and DNV-ST-F101 Table 7-5. */
const GRADES = {
  // API 5L / ISO 3183 PSL 2
  B: [245e6, 415e6], X42: [290e6, 415e6], X46: [320e6, 435e6], X52: [360e6, 460e6], X56: [390e6, 490e6],
  X60: [415e6, 520e6], X65: [450e6, 535e6], X70: [485e6, 570e6], X80: [555e6, 625e6],
  // DNV-ST-F101 C-Mn line pipe
  'DNV 245': [245e6, 370e6], 'DNV 290': [290e6, 415e6], 'DNV 320': [320e6, 435e6], 'DNV 360': [360e6, 460e6],
  'DNV 390': [390e6, 490e6], 'DNV 415': [415e6, 520e6], 'DNV 450': [450e6, 535e6], 'DNV 485': [485e6, 570e6],
  'DNV 555': [555e6, 625e6],
};

/**
 * Line-pipe grade record: SMYS and SMTS with common elastic constants.
 *
 * Reference: API 5L PSL 2 Table 7 (B, X42–X80); DNV-ST-F101 Table 7-5 (DNV 245–555).
 * The elastic constants are the customary carbon-steel design values and are returned for convenience.
 *
 * @param {string} grade e.g. 'X65' or 'DNV 450' (case-insensitive; 'L450' and '450' are accepted as DNV 450)
 * @returns {{grade:string, smys:number, smts:number, E:number, nu:number, alpha:number, rho:number}}
 *   smys, smts [Pa]; E [Pa]; nu [-]; alpha thermal expansion [1/K]; rho [kg/m³]
 */
export function steelGrade(grade) {
  if (typeof grade !== 'string') throw new RangeError('grade must be a string such as "X65" or "DNV 450"');
  let key = grade.trim().toUpperCase().replace(/^API\s*5L\s*/, '').replace(/\s+/g, ' ');
  const m = key.match(/^(?:DNV\s*|L)?(\d{3})$/);
  if (m) key = `DNV ${m[1]}`;
  const rec = GRADES[key];
  if (!rec) throw new RangeError(`grade "${grade}" is not in the library (${Object.keys(GRADES).join(', ')})`);
  return { grade: key, smys: rec[0], smts: rec[1], E: 207e9, nu: 0.3, alpha: 11.7e-6, rho: 7850 };
}

/** @returns {string[]} names of the grades in the library */
export function listSteelGrades() {
  return Object.keys(GRADES);
}

/**
 * Temperature de-rating of yield (and tensile) strength for C-Mn and 13Cr line pipe.
 *
 * Default curve: zero up to 50 °C, rising linearly to 30 MPa at 100 °C and to 70 MPa at 200 °C
 * (DNV-ST-F101 (2021) Fig. 5-2, C-Mn/13Cr line, read from the figure). Duplex grades de-rate from 20 °C and
 * are not built in: pass the project curve through `curve`.
 *
 * Validity: up to the last point of the curve (200 °C for the default); above it throws.
 *
 * @param {number} tempC design temperature [°C]
 * @param {object} [o]
 * @param {Array<[number, number]>} [o.curve] override: ascending [temperature °C, de-rating Pa] points,
 *   linearly interpolated, zero below the first point
 * @returns {number} de-rating value f_temp [Pa] to subtract from SMYS / SMTS
 */
export function temperatureDerating(tempC, { curve = [[50, 0], [100, 30e6], [200, 70e6]] } = {}) {
  if (!Number.isFinite(tempC)) throw new RangeError('tempC must be finite');
  if (!Array.isArray(curve) || curve.length < 2) throw new RangeError('curve must hold at least two [tempC, Pa] points');
  if (tempC <= curve[0][0]) return curve[0][1];
  const last = curve[curve.length - 1];
  if (tempC > last[0]) throw new RangeError(`tempC exceeds the de-rating curve (${last[0]} °C)`);
  for (let i = 1; i < curve.length; i++) {
    if (tempC <= curve[i][0]) {
      const [x0, y0] = curve[i - 1];
      const [x1, y1] = curve[i];
      return y0 + ((y1 - y0) * (tempC - x0)) / (x1 - x0);
    }
  }
  return last[1];
}

/**
 * Fabrication factor α_fab for collapse and propagation capacity.
 *
 * Reference: DNV-ST-F101 (2021) Table 5-5: seamless 1.00; UO, TRB, ERW and HFW 0.93; UOE 0.85.
 *
 * @param {'seamless'|'UO'|'TRB'|'ERW'|'HFW'|'UOE'} process pipe manufacturing process
 * @returns {number} α_fab [-]
 */
export function fabricationFactor(process) {
  const p = String(process).toUpperCase();
  if (p === 'SEAMLESS' || p === 'SMLS') return 1.0;
  if (p === 'UO' || p === 'TRB' || p === 'ERW' || p === 'HFW') return 0.93;
  if (p === 'UOE') return 0.85;
  throw new RangeError('process must be one of seamless, UO, TRB, ERW, HFW, UOE');
}

/**
 * Material strength factor α_U.
 *
 * Reference: DNV-ST-F101 (2021) Table 5-3: 0.96 normally, 1.00 with supplementary requirement U.
 *
 * @param {boolean} [supplementaryU=false] true when supplementary requirement U is specified
 * @returns {number} α_U [-]
 */
export function materialStrengthFactor(supplementaryU = false) {
  return supplementaryU ? 1.0 : 0.96;
}

/**
 * Characteristic material strengths f_y = (SMYS − f_y,temp)·α_U and f_u = (SMTS − f_u,temp)·α_U.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.4 and 5.5.
 *
 * @param {object} o
 * @param {number} o.smys specified minimum yield strength [Pa]
 * @param {number} o.smts specified minimum tensile strength [Pa]
 * @param {number} [o.tempC=20] design temperature [°C]
 * @param {number} [o.alphaU=0.96] material strength factor [-]
 * @param {number} [o.fyTemp] yield de-rating override [Pa] (default from temperatureDerating)
 * @param {number} [o.fuTemp] tensile de-rating override [Pa] (default equal to the yield de-rating)
 * @param {Array<[number, number]>} [o.curve] de-rating curve override, see temperatureDerating
 * @returns {{fy:number, fu:number, fyTemp:number, fuTemp:number, alphaU:number}} strengths [Pa]
 */
export function characteristicStrength({ smys, smts, tempC = 20, alphaU = 0.96, fyTemp, fuTemp, curve }) {
  pos('smys', smys); pos('smts', smts); pos('alphaU', alphaU);
  const dy = fyTemp === undefined ? temperatureDerating(tempC, curve ? { curve } : undefined) : nonNeg('fyTemp', fyTemp);
  const du = fuTemp === undefined ? dy : nonNeg('fuTemp', fuTemp);
  return { fy: (smys - dy) * alphaU, fu: (smts - du) * alphaU, fyTemp: dy, fuTemp: du, alphaU };
}

// ASME B36.10M: outside diameter [in] by NPS (exact up to NPS 12; equal to the NPS from 14 upward).
const NPS_OD_IN = [
  ['1/8', 0.405], ['1/4', 0.540], ['3/8', 0.675], ['1/2', 0.840], ['3/4', 1.050], ['1', 1.315], ['1 1/4', 1.660],
  ['1 1/2', 1.900], ['2', 2.375], ['2 1/2', 2.875], ['3', 3.500], ['3 1/2', 4.000], ['4', 4.500], ['5', 5.563],
  ['6', 6.625], ['8', 8.625], ['10', 10.750], ['12', 12.750],
];
for (const n of [14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48]) NPS_OD_IN.push([String(n), n]);

// ASME B36.10M wall thickness [in] by NPS: { STD, XS, 40, 80, 160 } (blank where the schedule does not exist).
const SCHEDULES_IN = {
  '1/2': { STD: 0.109, XS: 0.147, 40: 0.109, 80: 0.147, 160: 0.188 },
  '3/4': { STD: 0.113, XS: 0.154, 40: 0.113, 80: 0.154, 160: 0.219 },
  '1': { STD: 0.133, XS: 0.179, 40: 0.133, 80: 0.179, 160: 0.250 },
  '1 1/4': { STD: 0.140, XS: 0.191, 40: 0.140, 80: 0.191, 160: 0.250 },
  '1 1/2': { STD: 0.145, XS: 0.200, 40: 0.145, 80: 0.200, 160: 0.281 },
  '2': { STD: 0.154, XS: 0.218, 40: 0.154, 80: 0.218, 160: 0.344 },
  '2 1/2': { STD: 0.203, XS: 0.276, 40: 0.203, 80: 0.276, 160: 0.375 },
  '3': { STD: 0.216, XS: 0.300, 40: 0.216, 80: 0.300, 160: 0.438 },
  '3 1/2': { STD: 0.226, XS: 0.318, 40: 0.226, 80: 0.318 },
  '4': { STD: 0.237, XS: 0.337, 40: 0.237, 80: 0.337, 160: 0.531 },
  '5': { STD: 0.258, XS: 0.375, 40: 0.258, 80: 0.375, 160: 0.625 },
  '6': { STD: 0.280, XS: 0.432, 40: 0.280, 80: 0.432, 160: 0.719 },
  '8': { STD: 0.322, XS: 0.500, 40: 0.322, 80: 0.500, 160: 0.906 },
  '10': { STD: 0.365, XS: 0.500, 40: 0.365, 80: 0.594, 160: 1.125 },
  '12': { STD: 0.375, XS: 0.500, 40: 0.406, 80: 0.688, 160: 1.312 },
  '14': { STD: 0.375, XS: 0.500, 40: 0.438, 80: 0.750, 160: 1.406 },
  '16': { STD: 0.375, XS: 0.500, 40: 0.500, 80: 0.844, 160: 1.594 },
  '18': { STD: 0.375, XS: 0.500, 40: 0.562, 80: 0.938, 160: 1.781 },
  '20': { STD: 0.375, XS: 0.500, 40: 0.594, 80: 1.031, 160: 1.969 },
  '24': { STD: 0.375, XS: 0.500, 40: 0.688, 80: 1.219, 160: 2.344 },
};
for (const n of [22, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48]) SCHEDULES_IN[String(n)] = { STD: 0.375, XS: 0.500 };

const IN = 0.0254;

/**
 * Outside diameter for a nominal pipe size.
 *
 * Reference: ASME B36.10M (the inch outside diameters are exact; converted with 25.4 mm/in).
 * Validity: NPS 1/8 to 48 (even sizes above 24).
 *
 * @param {string|number} nps nominal pipe size, e.g. '12', 12, '1 1/2', 0.5
 * @returns {number} outside diameter [m]
 */
export function npsOuterDiameter(nps) {
  const key = normaliseNps(nps);
  const row = NPS_OD_IN.find((r) => r[0] === key);
  if (!row) throw new RangeError(`nps "${nps}" is not a standard size`);
  return row[1] * IN;
}

function normaliseNps(nps) {
  if (typeof nps === 'number') {
    const fr = { 0.125: '1/8', 0.25: '1/4', 0.375: '3/8', 0.5: '1/2', 0.75: '3/4', 1.25: '1 1/4', 1.5: '1 1/2', 2.5: '2 1/2', 3.5: '3 1/2' };
    return fr[nps] || String(nps);
  }
  return String(nps).trim().replace(/"/g, '').replace(/-/g, ' ').replace(/\s+/g, ' ');
}

/**
 * Standard pipe table: outside diameter by NPS and wall thickness by schedule.
 *
 * Reference: ASME B36.10M. Wall thicknesses are held for schedules STD, XS, 40, 80 and 160 from NPS 1/2 to 24
 * and STD / XS up to NPS 48; other schedules are not included (line pipe to API 5L may be ordered in any
 * thickness, so treat the table as a starting point).
 *
 * @param {object} [o]
 * @param {Array<string|number>} [o.sizes] subset of NPS designations (default: all)
 * @param {Array<string|number>} [o.schedules=['STD','XS','40','80','160']] schedules to list
 * @returns {Array<{nps:string, D:number, walls:Array<{schedule:string, t:number, Di:number, massSteel:number}>}>}
 *   D [m]; t [m]; Di [m]; massSteel bare steel mass at 7850 kg/m³ [kg/m]
 */
export function standardPipeTable({ sizes, schedules = ['STD', 'XS', '40', '80', '160'] } = {}) {
  const keys = sizes ? sizes.map(normaliseNps) : NPS_OD_IN.map((r) => r[0]);
  return keys.map((k) => {
    const D = npsOuterDiameter(k);
    const sch = SCHEDULES_IN[k] || {};
    const walls = [];
    for (const s of schedules) {
      const tin = sch[String(s).toUpperCase()];
      if (tin === undefined) continue;
      const t = tin * IN;
      walls.push({ schedule: String(s).toUpperCase(), t, Di: D - 2 * t, massSteel: 7850 * Math.PI * (D - t) * t });
    }
    return { nps: k, D, walls };
  });
}

/**
 * Wall thickness for an NPS / schedule pair.
 *
 * Reference: ASME B36.10M (subset described under standardPipeTable).
 *
 * @param {string|number} nps nominal pipe size
 * @param {string|number} schedule 'STD', 'XS', 40, 80 or 160
 * @returns {number} wall thickness [m]
 */
export function scheduleWallThickness(nps, schedule) {
  const k = normaliseNps(nps);
  const sch = SCHEDULES_IN[k];
  const tin = sch && sch[String(schedule).toUpperCase()];
  if (tin === undefined) throw new RangeError(`schedule ${schedule} is not tabulated for nps ${nps}`);
  return tin * IN;
}
