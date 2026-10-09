/**
 * Annulus environment of unbonded flexible pipes: gas permeation through the pressure sheath, pressure build-up
 * and venting, composition and sour/sweet classification, flooding, annulus volume testing and rapid gas
 * decompression screening.
 *
 * Units: permeability coefficients Pe in mol m /(m^2 s Pa) (amount of gas x thickness / area / time / pressure
 * difference); use {@link permeabilityFromCm3} to convert the customary cm^3(STP) cm /(cm^2 s bar). Permeability
 * data are strongly material-, grade- and plasticiser-dependent: no built-in values are provided.
 */

/** Molar gas constant, J/(mol K) (CODATA 2018, exact). */
export const R_GAS = 8.314462618;
/** ISO 15156-2 / NACE MR0175 threshold for sour service: H2S partial pressure 0.3 kPa (0.05 psi), Pa. */
export const H2S_SOUR_LIMIT = 300;

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
 * Convert a permeability coefficient from cm^3(STP) cm /(cm^2 s bar) to mol m /(m^2 s Pa).
 * 1 cm^3(STP) = 1/22414 mol (ideal gas at 273.15 K, 101.325 kPa); 1 cm /(cm^2 bar) = 1e-2 m /(1e-4 m^2 1e5 Pa).
 * @param {number} value coefficient in cm^3(STP) cm /(cm^2 s bar)
 * @returns {{Pe:number}} mol m /(m^2 s Pa)
 */
export function permeabilityFromCm3(value) {
  nonneg(value, 'value');
  return { Pe: value * 1e-3 / 22414 };
}

/**
 * Temperature dependence of a permeability coefficient: Pe(T) = Pe0 exp(-Ea /(R T)).
 * Reference: Arrhenius form for permeation in polymers (e.g. Crank & Park, Diffusion in Polymers, Ch. 2).
 * Pe0 and Ea are material/gas-pair data from permeation tests (arguments; activation energies for CH4, CO2,
 * H2S, H2O in PE/PA/PVDF are typically of order 20-60 kJ/mol).
 * @param {object} o
 * @param {number} o.Pe0 pre-exponential coefficient, mol m /(m^2 s Pa)
 * @param {number} o.Ea activation energy, J/mol
 * @param {number} o.T absolute temperature, K
 * @returns {{Pe:number}} mol m /(m^2 s Pa)
 */
export function arrheniusPermeability({ Pe0, Ea, T }) {
  nonneg(Pe0, 'Pe0'); nonneg(Ea, 'Ea'); pos(T, 'T');
  return { Pe: Pe0 * Math.exp(-Ea / (R_GAS * T)) };
}

/**
 * Steady radial permeation of gas components through a cylindrical sheath:
 *     n_dot_i = Pe_i (f_bore,i - f_annulus,i) 2 pi L / ln(ro/ri)       (mol/s)
 * with f the fugacity (partial pressure for an ideal gas) of component i on each side. Each component may give Pe
 * directly or (Pe0, Ea) for {@link arrheniusPermeability} at temperature T.
 * Reference: Fick's law / solution-diffusion model integrated over a cylinder wall (Crank, The Mathematics of
 * Diffusion, §5.2); API RP 17B §6 (annulus permeation calculation) in intent.
 * Validity: steady state, constant coefficient through the wall (evaluate Pe at the mean sheath temperature),
 * no shielding by the armour layers (conservative: shielding by the pressure armour reduces the effective area).
 * @param {object} o
 * @param {{name:string, Pe?:number, Pe0?:number, Ea?:number, fugacityBore:number, fugacityAnnulus?:number}[]}
 *   o.components fugacities in Pa
 * @param {number} o.ri sheath inner radius, m
 * @param {number} o.ro sheath outer radius, m
 * @param {number} o.length pipe length, m
 * @param {number} [o.T] sheath temperature, K (needed with Pe0/Ea)
 * @param {number} [o.areaFactor=1] effective-area factor for shielding by the pressure armour, (0, 1]
 * @returns {{rates:{name:string, Pe:number, molarRate:number}[], totalMolarRate:number,
 *   moleFractions:Record<string, number>}} mol/s; mole fractions of the permeate stream
 */
export function annulusPermeation({ components, ri, ro, length, T, areaFactor = 1 }) {
  if (!Array.isArray(components) || components.length === 0) throw new RangeError('components must be a non-empty array');
  pos(ri, 'ri'); pos(ro, 'ro'); pos(length, 'length');
  if (!(ro > ri)) throw new RangeError('ro must exceed ri');
  if (!(areaFactor > 0 && areaFactor <= 1)) throw new RangeError('areaFactor must be in (0, 1]');
  const geom = 2 * Math.PI * length / Math.log(ro / ri) * areaFactor;
  let total = 0;
  const rates = components.map((c, i) => {
    const Pe = c.Pe !== undefined ? nonneg(c.Pe, `components[${i}].Pe`) : arrheniusPermeability({ Pe0: c.Pe0, Ea: c.Ea, T }).Pe;
    const df = nonneg(c.fugacityBore, `components[${i}].fugacityBore`) - nonneg(c.fugacityAnnulus ?? 0, `components[${i}].fugacityAnnulus`);
    const molarRate = Pe * df * geom;
    total += molarRate;
    return { name: c.name ?? `component ${i + 1}`, Pe, molarRate };
  });
  const moleFractions = {};
  for (const r of rates) moleFractions[r.name] = total !== 0 ? r.molarRate / total : 0;
  return { rates, totalMolarRate: total, moleFractions };
}

/**
 * Free (gas) volume of the annulus: sum over the layers between the pressure sheath and the outer sheath of
 * (1 - fill) x layer annulus area x length, plus any explicit gap volume.
 * @param {object} o
 * @param {{rIn:number, rOut:number, fill:number}[]} o.layers annulus layers: radii m, solid fraction 0..1
 * @param {number} o.length m
 * @param {number} [o.extraAreaVoid=0] additional void area (gaps between layers), m^2
 * @returns {{voidArea:number, volume:number, annulusArea:number, voidFraction:number}} m^2, m^3, m^2, -
 */
export function annulusFreeVolume({ layers, length, extraAreaVoid = 0 }) {
  if (!Array.isArray(layers)) throw new RangeError('layers must be an array');
  pos(length, 'length'); nonneg(extraAreaVoid, 'extraAreaVoid');
  let area = 0, voidA = extraAreaVoid;
  layers.forEach((L, i) => {
    pos(L.rIn, `layers[${i}].rIn`); pos(L.rOut, `layers[${i}].rOut`);
    if (!(L.rOut > L.rIn)) throw new RangeError(`layers[${i}].rOut must exceed rIn`);
    if (!(L.fill >= 0 && L.fill <= 1)) throw new RangeError(`layers[${i}].fill must be in [0, 1]`);
    const a = Math.PI * (L.rOut ** 2 - L.rIn ** 2);
    area += a; voidA += (1 - L.fill) * a;
  });
  return { voidArea: voidA, volume: voidA * length, annulusArea: area, voidFraction: area > 0 ? voidA / (area + extraAreaVoid) : 1 };
}

/**
 * Annulus pressure build-up between vent-valve openings (ideal gas, constant temperature and inflow):
 *     dp/dt = n_dot R T / V,   t_vent = (p_vent - p_0) V / (n_dot R T),
 * and the time-averaged vent flow at standard conditions Q_std = n_dot R T_std / p_std.
 * Reference: ideal-gas law; API Spec 17J requires the annulus vent system to open before the annulus pressure
 * can lift the outer sheath (vent opening pressures are typically a few bar above ambient - argument).
 * Validity: permeation rate unaffected by the (small) annulus back-pressure; no condensation.
 * @param {object} o
 * @param {number} o.molarRate total permeation rate into the annulus, mol/s
 * @param {number} o.volume annulus free volume, m^3
 * @param {number} o.T annulus gas temperature, K
 * @param {number} [o.p0=101325] pressure after venting (vent re-seat), Pa
 * @param {number} o.pVent vent-valve opening pressure, Pa (absolute)
 * @param {number} [o.Tstd=288.15] K
 * @param {number} [o.pStd=101325] Pa
 * @returns {{pressureRate:number, timeToVent:number, ventInterval:number, ventRateStd:number,
 *   molesPerVent:number}} Pa/s, s, s, m^3(std)/s, mol
 */
export function annulusPressureBuildUp({ molarRate, volume, T, p0 = 101325, pVent, Tstd = 288.15, pStd = 101325 }) {
  nonneg(molarRate, 'molarRate'); pos(volume, 'volume'); pos(T, 'T'); pos(p0, 'p0'); pos(pVent, 'pVent'); pos(Tstd, 'Tstd'); pos(pStd, 'pStd');
  if (!(pVent > p0)) throw new RangeError('pVent must exceed p0');
  const rate = molarRate * R_GAS * T / volume;
  const t = rate > 0 ? (pVent - p0) / rate : Infinity;
  return { pressureRate: rate, timeToVent: t, ventInterval: t, ventRateStd: molarRate * R_GAS * Tstd / pStd, molesPerVent: (pVent - p0) * volume / (R_GAS * T) };
}

/**
 * Saturation pressure of water by the Antoine equation log10(p[mmHg]) = 8.07131 - 1730.63/(233.426 + T[degC]),
 * valid 1-100 degC (constants as tabulated in the DDBST / Dortmund Data Bank and widely reproduced).
 * @param {number} TC temperature, degC
 * @returns {{pSat:number}} Pa
 */
export function waterSaturationPressure(TC) {
  fin(TC, 'TC');
  if (TC < 1 || TC > 100) throw new RangeError('TC must be between 1 and 100 degC');
  return { pSat: 10 ** (8.07131 - 1730.63 / (233.426 + TC)) * 133.322368 };
}

/**
 * Steady-state annulus composition and classification. With continuous permeation and venting the annulus gas
 * composition tends to the mole fractions of the permeate stream; partial pressure = mole fraction x annulus
 * pressure. Classification: 'sour' when p_H2S >= 0.3 kPa (ISO 15156-2 / NACE MR0175 threshold,
 * {@link H2S_SOUR_LIMIT}), else 'sweet' when CO2 is present, else 'benign'. Water condenses when its partial
 * pressure exceeds the saturation pressure at the annulus temperature.
 * Consequences: a condensed-water or seawater-flooded annulus with CO2/H2S is a corrosion-fatigue environment
 * for the armour wires - S-N data for that environment are required (flexfatigue.js), and sour conditions call
 * for wires qualified to ISO 15156 (SSC/HIC).
 * @param {object} o
 * @param {Record<string, number>} o.moleFractions by component name; the names 'H2S', 'CO2', 'H2O' are recognised
 * @param {number} o.pAnnulus annulus absolute pressure, Pa
 * @param {number} [o.TC] annulus temperature, degC (for the condensation check)
 * @param {number} [o.sourLimit=300] Pa
 * @returns {{partialPressures:Record<string, number>, classification:'sour'|'sweet'|'benign', h2sRatio:number,
 *   waterCondenses:boolean|null}} Pa; p_H2S / sour limit
 */
export function annulusClassification({ moleFractions, pAnnulus, TC, sourLimit = H2S_SOUR_LIMIT }) {
  if (!moleFractions || typeof moleFractions !== 'object') throw new RangeError('moleFractions must be an object');
  pos(pAnnulus, 'pAnnulus'); pos(sourLimit, 'sourLimit');
  const partialPressures = {};
  for (const [k, v] of Object.entries(moleFractions)) partialPressures[k] = nonneg(v, `moleFractions.${k}`) * pAnnulus;
  const pH2S = partialPressures.H2S ?? 0, pCO2 = partialPressures.CO2 ?? 0;
  const classification = pH2S >= sourLimit ? 'sour' : pCO2 > 0 ? 'sweet' : 'benign';
  const waterCondenses = TC === undefined || partialPressures.H2O === undefined ? null : partialPressures.H2O > waterSaturationPressure(TC).pSat;
  return { partialPressures, classification, h2sRatio: pH2S / sourLimit, waterCondenses };
}

/**
 * pH estimate of condensed water in the annulus (hook). For pure water in equilibrium with CO2 the carbonic-acid
 * balance [H+]^2 = K1 KH pCO2 gives pH = (pK1 + pKH - log10(pCO2/atm))/2, with pK1 = 6.35 and
 * pKH = 1.47 (KH in mol/(L atm)) at 25 degC (Stumm & Morgan, Aquatic Chemistry, carbonate equilibria). In a confined
 * annulus the large steel-surface-to-water-volume ratio supersaturates the water with dissolved iron and raises
 * the pH markedly (values around 5.5-6.5 are reported); that effect is not predictable from first principles here
 * and enters as `confinementShift` (argument) or a measured `pHOverride`.
 * @param {object} o
 * @param {number} o.pCO2 CO2 partial pressure, Pa
 * @param {number} [o.pK1=6.35] first dissociation constant of carbonic acid
 * @param {number} [o.pKH=1.47] Henry constant exponent
 * @param {number} [o.confinementShift=0] pH increase from iron supersaturation in confined water
 * @param {number} [o.pHOverride] measured/assumed pH returned as is
 * @returns {{pH:number, pHPureWater:number}}
 */
export function annulusPH({ pCO2, pK1 = 6.35, pKH = 1.47, confinementShift = 0, pHOverride }) {
  pos(pCO2, 'pCO2'); nonneg(confinementShift, 'confinementShift');
  const pure = Math.min(0.5 * (pK1 + pKH - Math.log10(pCO2 / 101325)), 7);
  return { pH: pHOverride !== undefined ? fin(pHOverride, 'pHOverride') : Math.min(pure + confinementShift, 14), pHPureWater: pure };
}

/**
 * Outer-sheath damage / seawater flooding scenario: water volume and mass entering the annulus free volume, the
 * annulus pressure (hydrostatic at the breach), and the flags that change the design checks.
 * Consequences reported: wet collapse (carcass alone resists, flexcollapse.js), loss of external-pressure
 * support against bird-caging, corrosion-fatigue S-N curves for the armour, added weight.
 * @param {object} o
 * @param {number} o.freeVolume annulus free volume over the flooded length, m^3
 * @param {number} o.floodedLength m
 * @param {number} o.breachDepth water depth at the damage, m
 * @param {number} [o.floodedFraction=1] fraction of the free volume filled
 * @param {number} [o.rho=1025] kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @param {number} [o.pAtm=101325] Pa
 * @returns {{waterVolume:number, addedMassPerLength:number, addedWeightPerLength:number, annulusPressure:number,
 *   consequences:string[]}} m^3, kg/m, N/m (in air; submerged the water is neutrally buoyant), Pa abs
 */
export function annulusFlooding({ freeVolume, floodedLength, breachDepth, floodedFraction = 1, rho = 1025, g = 9.80665, pAtm = 101325 }) {
  nonneg(freeVolume, 'freeVolume'); pos(floodedLength, 'floodedLength'); nonneg(breachDepth, 'breachDepth'); pos(rho, 'rho'); pos(g, 'g');
  if (!(floodedFraction >= 0 && floodedFraction <= 1)) throw new RangeError('floodedFraction must be in [0, 1]');
  const waterVolume = freeVolume * floodedFraction;
  const m = rho * waterVolume / floodedLength;
  return {
    waterVolume, addedMassPerLength: m, addedWeightPerLength: m * g, annulusPressure: pAtm + rho * g * breachDepth,
    consequences: [
      'External pressure acts on the pressure sheath: check wet collapse (carcass alone).',
      'No external-pressure support of the tensile armour: re-check radial and lateral buckling under compression.',
      'Seawater with permeated CO2/H2S: use corrosion-fatigue S-N data for the armour wires.',
      'Vent system and annulus monitoring no longer indicate permeation: plan annulus testing/inspection.',
    ],
  };
}

/**
 * Annulus free-volume test evaluation.
 *  - 'pressure': a known amount of gas (volume V_inj at standard conditions) raises the annulus pressure from p1
 *    to p2 at temperature T: V_free = n R T/(p2 - p1), n = pStd V_inj/(R Tstd).
 *  - 'vacuum': the evacuated annulus at p1 is connected to a reference vessel V_ref at p_ref and equalises at p2
 *    (Boyle's law, isothermal): V_free = V_ref (p_ref - p2)/(p2 - p1).
 * Flooded fraction = 1 - V_free/V_design; a leak is estimated from a pressure drift dp/dt at constant temperature
 * as n_dot = (dp/dt) V_free/(R T).
 * Reference: ideal-gas law (standard annulus vacuum / nitrogen volume test practice for flexible risers).
 * @param {object} o
 * @param {'pressure'|'vacuum'} o.method
 * @param {number} o.p1 initial annulus pressure, Pa abs
 * @param {number} o.p2 final annulus pressure, Pa abs
 * @param {number} [o.T=288.15] gas temperature, K
 * @param {number} [o.injectedStdVolume] m^3 at (Tstd, pStd), for 'pressure'
 * @param {number} [o.refVolume] m^3, for 'vacuum'
 * @param {number} [o.refPressure] Pa abs, for 'vacuum'
 * @param {number} [o.designVolume] m^3 dry design free volume
 * @param {number} [o.pressureDrift] dp/dt during a hold period, Pa/s
 * @param {number} [o.Tstd=288.15]
 * @param {number} [o.pStd=101325]
 * @returns {{freeVolume:number, floodedFraction:number|null, leakMolarRate:number|null,
 *   leakStdRate:number|null}} m^3, -, mol/s, m^3(std)/s
 */
export function annulusVolumeTest({ method, p1, p2, T = 288.15, injectedStdVolume, refVolume, refPressure, designVolume, pressureDrift, Tstd = 288.15, pStd = 101325 }) {
  nonneg(p1, 'p1'); pos(p2, 'p2'); pos(T, 'T');
  if (!(p2 > p1)) throw new RangeError('p2 must exceed p1');
  let V;
  if (method === 'pressure') {
    const n = pStd * pos(injectedStdVolume, 'injectedStdVolume') / (R_GAS * Tstd);
    V = n * R_GAS * T / (p2 - p1);
  } else if (method === 'vacuum') {
    pos(refVolume, 'refVolume'); pos(refPressure, 'refPressure');
    if (!(refPressure > p2)) throw new RangeError('refPressure must exceed p2');
    V = refVolume * (refPressure - p2) / (p2 - p1);
  } else throw new RangeError("method must be 'pressure' or 'vacuum'");
  const leak = pressureDrift === undefined ? null : fin(pressureDrift, 'pressureDrift') * V / (R_GAS * T);
  return {
    freeVolume: V,
    floodedFraction: designVolume === undefined ? null : 1 - V / pos(designVolume, 'designVolume'),
    leakMolarRate: leak, leakStdRate: leak === null ? null : leak * R_GAS * Tstd / pStd,
  };
}

/**
 * Rapid gas decompression (blistering) risk screening of a polymer sheath: compares the bore depressurisation
 * rate and pressure drop with limits for which the material has been qualified (blistering tests per API Spec
 * 17J are run at a stated pressure, temperature, gas composition, number of cycles and decompression
 * rate - those qualified values are the arguments).
 * @param {object} o
 * @param {number} o.pressureDrop depressurisation step, Pa
 * @param {number} o.duration time taken, s
 * @param {number} o.qualifiedRate decompression rate covered by the qualification test, Pa/s
 * @param {number} [o.qualifiedPressure] maximum pressure covered by the qualification, Pa
 * @param {number} [o.initialPressure] pressure before blow-down, Pa
 * @returns {{rate:number, rateUtilisation:number, pressureUtilisation:number|null, risk:'within qualification'|
 *   'outside qualification'}} Pa/s
 */
export function rgdScreening({ pressureDrop, duration, qualifiedRate, qualifiedPressure, initialPressure }) {
  nonneg(pressureDrop, 'pressureDrop'); pos(duration, 'duration'); pos(qualifiedRate, 'qualifiedRate');
  const rate = pressureDrop / duration;
  const pu = qualifiedPressure !== undefined && initialPressure !== undefined ? nonneg(initialPressure, 'initialPressure') / pos(qualifiedPressure, 'qualifiedPressure') : null;
  const ru = rate / qualifiedRate;
  return { rate, rateUtilisation: ru, pressureUtilisation: pu, risk: ru <= 1 && (pu === null || pu <= 1) ? 'within qualification' : 'outside qualification' };
}
