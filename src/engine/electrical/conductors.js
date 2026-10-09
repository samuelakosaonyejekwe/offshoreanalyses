/**
 * Conductor and cable-material data for power-cable calculations.
 *
 * Sources:
 *  - IEC 60228:2004 Table 2 (class 2 stranded conductors, maximum DC resistance at 20 °C, plain annealed copper
 *    and aluminium).
 *  - IEC 60287-1-1:2006 Table 1 (electrical resistivities and temperature coefficients), Table 2 (skin- and
 *    proximity-effect coefficients ks, kp), Table 3 (relative permittivity and loss factor).
 *  - IEC 60287-2-1:2006 Table 1 (thermal resistivities of materials).
 *  - IEC 60853-2:1989 Table E1 (volumetric heat capacities).
 * Everything is SI: Ω/m, m, m², Ω·m, K·m/W, J/(m³·K).
 */

/** Standard nominal cross-sections of IEC 60228, mm². */
const SIZES_MM2 = [16, 25, 35, 50, 70, 95, 120, 150, 185, 240, 300, 400, 500, 630, 800, 1000, 1200, 1400, 1600,
  1800, 2000, 2500];

/** IEC 60228 Table 2, maximum DC resistance at 20 °C, Ω/km, plain copper wires, class 2. */
const R20_CU_OHM_PER_KM = [1.15, 0.727, 0.524, 0.387, 0.268, 0.193, 0.153, 0.124, 0.0991, 0.0754, 0.0601, 0.0470,
  0.0366, 0.0283, 0.0221, 0.0176, 0.0151, 0.0129, 0.0113, 0.0101, 0.0090, 0.0072];

/** IEC 60228 Table 2, maximum DC resistance at 20 °C, Ω/km, aluminium, class 2. */
const R20_AL_OHM_PER_KM = [1.91, 1.20, 0.868, 0.641, 0.443, 0.320, 0.253, 0.206, 0.164, 0.125, 0.100, 0.0778,
  0.0605, 0.0469, 0.0367, 0.0291, 0.0247, 0.0212, 0.0186, 0.0165, 0.0149, 0.0127];

/**
 * Electrical properties of metals used in cables at 20 °C (IEC 60287-1-1 Table 1) with the adiabatic
 * short-circuit constants of IEC 60949 Table I and volumetric heat capacities of IEC 60853-2 Table E1.
 *  rho20  electrical resistivity, Ω·m
 *  alpha20 temperature coefficient of resistance at 20 °C, 1/K
 *  K      IEC 60949 constant, A·s^0.5/mm²
 *  beta   reciprocal of the temperature coefficient at 0 °C, K
 *  heatCapacity volumetric specific heat, J/(m³·K)
 */
export const METALS = Object.freeze({
  Cu: Object.freeze({ name: 'Copper', rho20: 1.7241e-8, alpha20: 3.93e-3, K: 226, beta: 234.5, heatCapacity: 3.45e6 }),
  Al: Object.freeze({ name: 'Aluminium', rho20: 2.8264e-8, alpha20: 4.03e-3, K: 148, beta: 228, heatCapacity: 2.5e6 }),
  lead: Object.freeze({ name: 'Lead or lead alloy', rho20: 21.4e-8, alpha20: 4.0e-3, K: 41, beta: 230, heatCapacity: 1.45e6 }),
  steel: Object.freeze({ name: 'Steel', rho20: 13.8e-8, alpha20: 4.5e-3, K: 78, beta: 202, heatCapacity: 3.8e6 }),
  bronze: Object.freeze({ name: 'Bronze', rho20: 3.5e-8, alpha20: 3.0e-3, K: 180, beta: 313, heatCapacity: 3.4e6 }),
  stainless: Object.freeze({ name: 'Stainless steel', rho20: 70e-8, alpha20: 0, K: NaN, beta: NaN, heatCapacity: 3.8e6 }),
});

/**
 * Insulation and covering materials.
 *  rhoT  thermal resistivity, K·m/W (IEC 60287-2-1 Table 1)
 *  rhoTHighVoltage  value for the higher voltage class where the standard distinguishes one, K·m/W
 *  epsR  relative permittivity (IEC 60287-1-1 Table 3)
 *  tanDelta  loss factor for cables above 18/30 (36) kV; tanDeltaMv for cables up to 18/30 (36) kV
 *  heatCapacity  volumetric specific heat, J/(m³·K) (IEC 60853-2 Table E1)
 *  maxConductorC / maxShortCircuitC  usual maximum conductor temperatures in service and after a short circuit, °C
 *     (IEC 60502-2 / IEC 60840 values for XLPE, EPR and PVC; 70/130 °C is for low-density PE; 55 °C for paper is
 *     the usual mass-impregnated HVDC limit and no short-circuit value is given; all are editable project data)
 */
export const INSULATION = Object.freeze({
  xlpe: Object.freeze({ name: 'XLPE (unfilled)', rhoT: 3.5, epsR: 2.5, tanDelta: 0.001, tanDeltaMv: 0.004,
    heatCapacity: 2.4e6, maxConductorC: 90, maxShortCircuitC: 250 }),
  pe: Object.freeze({ name: 'Polyethylene', rhoT: 3.5, epsR: 2.3, tanDelta: 0.001, tanDeltaMv: 0.001,
    heatCapacity: 2.4e6, maxConductorC: 70, maxShortCircuitC: 130 }),
  epr: Object.freeze({ name: 'EPR', rhoT: 5.0, rhoTLowVoltage: 3.5, epsR: 3.0, tanDelta: 0.005, tanDeltaMv: 0.020,
    heatCapacity: 2.0e6, maxConductorC: 90, maxShortCircuitC: 250 }),
  pvc: Object.freeze({ name: 'PVC', rhoT: 6.0, rhoTLowVoltage: 5.0, epsR: 8.0, tanDelta: 0.1, tanDeltaMv: 0.1,
    heatCapacity: 1.7e6, maxConductorC: 70, maxShortCircuitC: 160 }),
  paper: Object.freeze({ name: 'Mass-impregnated paper (solid type)', rhoT: 6.0, epsR: 4.0, tanDelta: 0.01,
    tanDeltaMv: 0.01, heatCapacity: 2.0e6, maxConductorC: 55, maxShortCircuitC: NaN }),
});

/**
 * Thermal resistivities of protective coverings and duct materials, K·m/W (IEC 60287-2-1 Table 1).
 */
export const COVERING_RHO_T = Object.freeze({
  juteAndFibrous: 6.0, rubberSandwich: 6.0, polychloroprene: 5.5, pvcUpTo35kV: 5.0, pvcAbove35kV: 6.0, pe: 3.5,
  ductConcrete: 1.0, ductFibre: 4.8, ductEarthenware: 1.2, ductPvc: 6.0, ductPe: 3.5,
});

/**
 * Skin- and proximity-effect coefficients, IEC 60287-1-1:2006 Table 2.
 * `impregnated` distinguishes dried-and-impregnated (paper) conductors from extruded-insulation conductors.
 * The `milliken*2014` entries are the values introduced by Amendment 1:2014 for segmental copper conductors;
 * they are provided for convenience and should be confirmed against the cable maker's data.
 */
export const SKIN_PROXIMITY = Object.freeze({
  roundStranded: Object.freeze({ ks: 1, kp: 1 }),
  roundStrandedImpregnated: Object.freeze({ ks: 1, kp: 0.8 }),
  roundCompact: Object.freeze({ ks: 1, kp: 1 }),
  roundCompactImpregnated: Object.freeze({ ks: 1, kp: 0.8 }),
  sectorShaped: Object.freeze({ ks: 1, kp: 1 }),
  sectorShapedImpregnated: Object.freeze({ ks: 1, kp: 0.8 }),
  milliken: Object.freeze({ ks: 0.435, kp: 0.37 }),
  millikenInsulatedWires2014: Object.freeze({ ks: 0.35, kp: 0.20 }),
  millikenBareUnidirectional2014: Object.freeze({ ks: 0.62, kp: 0.37 }),
  millikenBareBidirectional2014: Object.freeze({ ks: 0.80, kp: 0.37 }),
});

function metal(material) {
  const m = METALS[material];
  if (!m) throw new RangeError(`material: unknown metal "${material}" (use one of ${Object.keys(METALS).join(', ')})`);
  return m;
}

/**
 * List of the standard IEC 60228 class-2 nominal cross-sections covered by the resistance table.
 * @returns {number[]} cross-sections, m²
 */
export function standardSizes() {
  return SIZES_MM2.map((s) => s * 1e-6);
}

/**
 * Maximum DC resistance at 20 °C of an IEC 60228 class-2 stranded conductor.
 * Tabulated sizes (16–2500 mm²) return the standard value exactly. Other sizes (including 3000 mm², which is
 * outside the IEC 60228 table) are computed as R = ρ20·kLay/A with a documented lay/stranding allowance.
 * Reference: IEC 60228:2004 Table 2.
 * @param {object} o
 * @param {'Cu'|'Al'} o.material conductor metal
 * @param {number} o.area nominal cross-section, m² (16e-6 … 3000e-6)
 * @param {number} [o.layFactor=1.02] resistance allowance for stranding lay and tolerance when the size is not
 *        tabulated (typical 1.02–1.04; the tabulated values are 1.02–1.13 times ρ20/A_nominal, the higher ratios at
 *        small sizes where the nominal area exceeds the metallic area)
 * @returns {{R20:number, tabulated:boolean}} R20 in Ω/m; tabulated=true when taken from IEC 60228
 */
export function dcResistance20({ material, area, layFactor = 1.02 }) {
  if (material !== 'Cu' && material !== 'Al') throw new RangeError('material: IEC 60228 table covers "Cu" and "Al" only');
  if (!(area > 0)) throw new RangeError('area: must be positive (m²)');
  if (!(layFactor >= 1)) throw new RangeError('layFactor: must be ≥ 1');
  const mm2 = area * 1e6;
  const table = material === 'Cu' ? R20_CU_OHM_PER_KM : R20_AL_OHM_PER_KM;
  for (let i = 0; i < SIZES_MM2.length; i++) {
    if (Math.abs(SIZES_MM2[i] - mm2) < 1e-6 * SIZES_MM2[i]) return { R20: table[i] * 1e-3, tabulated: true };
  }
  if (mm2 < 16 || mm2 > 3000) throw new RangeError('area: outside 16–3000 mm²');
  return { R20: (metal(material).rho20 * layFactor) / area, tabulated: false };
}

/**
 * DC resistance at temperature θ: R = R20·[1 + α20·(θ − 20)].
 * Reference: IEC 60287-1-1 §2.1.1.
 * @param {object} o
 * @param {number} o.R20 resistance at 20 °C, Ω/m (or Ω)
 * @param {number} o.thetaC temperature, °C
 * @param {number} [o.alpha20] temperature coefficient at 20 °C, 1/K (overrides material)
 * @param {string} [o.material] metal key of METALS (used when alpha20 is omitted)
 * @returns {number} resistance at θ, same unit as R20
 */
export function resistanceAtTemperature({ R20, thetaC, alpha20, material }) {
  if (!(R20 >= 0)) throw new RangeError('R20: must be ≥ 0');
  if (!Number.isFinite(thetaC)) throw new RangeError('thetaC: must be finite');
  const a = alpha20 ?? metal(material).alpha20;
  return R20 * (1 + a * (thetaC - 20));
}

/**
 * Estimated diameter of a stranded conductor from its cross-section: d = sqrt(4A/(π·fill)).
 * Geometric identity; the fill factor is the ratio of metal area to the circumscribed circle area.
 * @param {object} o
 * @param {number} o.area metallic cross-section, m²
 * @param {number} [o.fillFactor=0.9] typical 0.76 (uncompacted concentric strand) … 0.92 (compacted), 0.95–0.97
 *        (water-blocked compacted / profiled wire), 1.0 solid
 * @returns {number} conductor diameter, m
 */
export function conductorDiameter({ area, fillFactor = 0.9 }) {
  if (!(area > 0)) throw new RangeError('area: must be positive (m²)');
  if (!(fillFactor > 0 && fillFactor <= 1)) throw new RangeError('fillFactor: must be in (0, 1]');
  return Math.sqrt((4 * area) / (Math.PI * fillFactor));
}

/**
 * Resistance per unit length at 20 °C of a tubular metallic sheath: Rs = ρ/(π·d·t) (thin-wall, mean diameter).
 * Reference: IEC 60287-1-1 §2.3 (resistance of sheath from its cross-section).
 * @param {object} o
 * @param {string} o.material metal key of METALS
 * @param {number} o.meanDiameter mean sheath diameter, m
 * @param {number} o.thickness sheath thickness, m
 * @param {number} [o.rho20] resistivity override, Ω·m
 * @returns {number} resistance at 20 °C, Ω/m
 */
export function sheathResistance20({ material, meanDiameter, thickness, rho20 }) {
  if (!(meanDiameter > 0)) throw new RangeError('meanDiameter: must be positive (m)');
  if (!(thickness > 0)) throw new RangeError('thickness: must be positive (m)');
  const rho = rho20 ?? metal(material).rho20;
  return rho / (Math.PI * meanDiameter * thickness);
}

/**
 * Resistance per unit cable length at 20 °C of a wire armour layer: RA = ρ·kLay/(n·π·df²/4).
 * The lay factor accounts for the longer path of helically applied wires,
 * kLay = sqrt(1 + (π·dA/p)²) for lay length p; it is exposed as an argument.
 * @param {object} o
 * @param {string} o.material metal key of METALS
 * @param {number} o.wireDiameter wire diameter df, m
 * @param {number} o.nWires number of wires
 * @param {number} [o.layFactor=1.0] ratio of wire length to cable length (≥ 1; typical 1.02–1.10)
 * @param {number} [o.rho20] resistivity override, Ω·m
 * @returns {number} resistance at 20 °C, Ω/m
 */
export function armourResistance20({ material, wireDiameter, nWires, layFactor = 1.0, rho20 }) {
  if (!(wireDiameter > 0)) throw new RangeError('wireDiameter: must be positive (m)');
  if (!(nWires >= 1)) throw new RangeError('nWires: must be ≥ 1');
  if (!(layFactor >= 1)) throw new RangeError('layFactor: must be ≥ 1');
  const rho = rho20 ?? metal(material).rho20;
  return (rho * layFactor) / (nWires * Math.PI * wireDiameter * wireDiameter / 4);
}

/**
 * Number of round armour wires that fit on a pitch circle with a given coverage.
 * Geometric: n = floor(coverage·π·dA/df).
 * @param {object} o
 * @param {number} o.meanDiameter armour pitch-circle diameter dA, m
 * @param {number} o.wireDiameter wire diameter, m
 * @param {number} [o.coverage=0.95] fraction of the circumference covered (0–1)
 * @returns {number} number of wires (integer)
 */
export function armourWireCount({ meanDiameter, wireDiameter, coverage = 0.95 }) {
  if (!(meanDiameter > 0)) throw new RangeError('meanDiameter: must be positive (m)');
  if (!(wireDiameter > 0)) throw new RangeError('wireDiameter: must be positive (m)');
  if (!(coverage > 0 && coverage <= 1)) throw new RangeError('coverage: must be in (0, 1]');
  return Math.max(1, Math.floor((coverage * Math.PI * meanDiameter) / wireDiameter));
}

/**
 * Skin- and proximity-effect coefficients for a conductor construction.
 * Reference: IEC 60287-1-1:2006 Table 2 (see SKIN_PROXIMITY).
 * @param {string} type key of SKIN_PROXIMITY
 * @returns {{ks:number, kp:number}} dimensionless coefficients
 */
export function skinProximityCoefficients(type) {
  const c = SKIN_PROXIMITY[type];
  if (!c) throw new RangeError(`type: unknown conductor construction "${type}"`);
  return { ks: c.ks, kp: c.kp };
}

/**
 * Property record of a metal (resistivity, temperature coefficient, IEC 60949 constants, heat capacity).
 * @param {string} material key of METALS
 * @returns {{name:string, rho20:number, alpha20:number, K:number, beta:number, heatCapacity:number}}
 *          Ω·m, 1/K, A·s^0.5/mm², K, J/(m³·K)
 */
export function metalProperties(material) {
  return metal(material);
}

/**
 * Property record of an insulation material.
 * @param {string} material key of INSULATION
 * @returns {object} see INSULATION
 */
export function insulationProperties(material) {
  const m = INSULATION[material];
  if (!m) throw new RangeError(`material: unknown insulation "${material}" (use one of ${Object.keys(INSULATION).join(', ')})`);
  return m;
}
