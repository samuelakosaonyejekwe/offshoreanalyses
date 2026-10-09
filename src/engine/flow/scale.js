// Mineral-scale screening: aqueous speciation helpers, activity coefficients, temperature-dependent solubility
// products, saturation indices for carbonate and sulphate scales, seawater/formation-water mixing, classical water
// stability indices and a first-order deposition-rate screen.
// Concentrations: molality (mol/kg water) inside the thermodynamic functions; water analyses in mg/L are converted
// by `toMolal`. T in K, pressures in Pa.

import { brineDensity } from './water.js';

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/** Ion data: molar mass (g/mol), charge, and ion-size parameter å (Å, Kielland 1937) for extended Debye–Hückel. */
export const IONS = Object.freeze({
  Na: Object.freeze({ M: 22.98977, z: 1, a: 4.0 }),
  K: Object.freeze({ M: 39.0983, z: 1, a: 3.0 }),
  Ca: Object.freeze({ M: 40.078, z: 2, a: 6.0 }),
  Mg: Object.freeze({ M: 24.305, z: 2, a: 8.0 }),
  Ba: Object.freeze({ M: 137.327, z: 2, a: 5.0 }),
  Sr: Object.freeze({ M: 87.62, z: 2, a: 5.0 }),
  Fe: Object.freeze({ M: 55.845, z: 2, a: 6.0 }),
  Cl: Object.freeze({ M: 35.453, z: -1, a: 3.0 }),
  SO4: Object.freeze({ M: 96.0626, z: -2, a: 4.0 }),
  HCO3: Object.freeze({ M: 61.0168, z: -1, a: 4.0 }),
  CO3: Object.freeze({ M: 60.0089, z: -2, a: 4.5 }),
});

/**
 * Composition of standard seawater (salinity 35) in mg per kg of solution, after the reference composition of
 * Millero, Feistel, Wright & McDougall, Deep-Sea Res. I 55 (2008) 50–72 (rounded). Multiply by the density
 * (≈ 1.025 kg/L) for mg/L. Provided as an editable default for mixing studies.
 */
export const SEAWATER_MG_PER_KG = Object.freeze({ Na: 10781, K: 399, Mg: 1284, Ca: 412, Sr: 7.9, Cl: 19353,
  SO4: 2712, HCO3: 105 });

/**
 * Convert a water analysis in mg/L to molalities (mol/kg water): m_i = c_i/(1000 M_i) / [(ρ − TDS)/1000],
 * with ρ the solution density in g/L (from the McCain brine correlation at 288.7 K when not given) and TDS the sum
 * of the listed ions in g/L.
 * Reference: definition of molality. @param {object} o
 * @param {Object<string, number>} o.mgL ion concentrations, mg/L, keyed as `IONS`
 * @param {number} [o.density] solution density, kg/m³
 * @returns {{ molal: Object<string, number>, tds: number, density: number, chargeBalance: number }} molalities;
 *   TDS in mg/L; density kg/m³; charge imbalance (Σ cations − Σ anions)/(Σ cations + Σ anions) in equivalents
 */
export function toMolal({ mgL, density }) {
  if (!mgL || typeof mgL !== 'object') throw new RangeError('mgL must be an object of ion concentrations');
  let tds = 0;
  for (const [ion, c] of Object.entries(mgL)) {
    if (!IONS[ion]) throw new RangeError(`ion ${ion} is not in the ion table`);
    tds += need(`mgL.${ion}`, c, 0, 5e5);
  }
  let rho = density;
  if (rho == null) {
    // salinity in g/kg needs the density itself: two fixed-point passes are ample
    rho = 1000;
    for (let pass = 0; pass < 3; pass++) rho = brineDensity({ T: 288.7, S: Math.min(tds / rho, 300), P: 101325 }).rhoSc;
  }
  const kgWaterPerL = (rho - tds / 1000) / 1000;
  if (!(kgWaterPerL > 0)) throw new RangeError('TDS exceeds the solution density');
  const molal = {};
  let cat = 0;
  let an = 0;
  for (const [ion, c] of Object.entries(mgL)) {
    molal[ion] = c / (1000 * IONS[ion].M) / kgWaterPerL;
    const eq = molal[ion] * Math.abs(IONS[ion].z);
    if (IONS[ion].z > 0) cat += eq; else an += eq;
  }
  return { molal, tds, density: rho, chargeBalance: cat + an > 0 ? (cat - an) / (cat + an) : 0 };
}

/**
 * Ionic strength I = ½ Σ m_i z_i² (molal).
 * Reference: Lewis & Randall (1921). @param {Object<string, number>} molal molalities keyed as `IONS`
 * @returns {number} mol/kg
 */
export function ionicStrength(molal) {
  let I = 0;
  for (const [ion, m] of Object.entries(molal)) {
    if (!IONS[ion]) throw new RangeError(`ion ${ion} is not in the ion table`);
    I += 0.5 * need(`molal.${ion}`, m, 0, 100) * IONS[ion].z ** 2;
  }
  return I;
}

/**
 * Static dielectric constant of water at atmospheric pressure, ε = 87.740 − 0.40008 t + 9.398e-4 t² − 1.410e-6 t³
 * (t °C). Reference: Malmberg & Maryott, J. Res. NBS 56 (1956) 1–8. Validity: 0–100 °C.
 * @param {number} T K @returns {number} relative permittivity, –
 */
export function waterDielectric(T) {
  need('T', T, 273.15, 373.16);
  const t = T - 273.15;
  return 87.740 - 0.40008 * t + 9.398e-4 * t * t - 1.410e-6 * t ** 3;
}

/**
 * Debye–Hückel parameters from first principles:
 * A = 1.82483e6 √ρ/(εT)^1.5 (kg^½·mol^−½), B = 50.2916 √ρ/(εT)^½ (Å⁻¹·kg^½·mol^−½), ρ in g/cm³.
 * At 25 °C: A = 0.5108 (log10 basis), B = 0.3286.
 * Reference: Debye & Hückel (1923); constants as in Truesdell & Jones, J. Res. USGS 2 (1974) 233 (WATEQ).
 * Validity: 0–100 °C (dielectric fit), pure-water density at 1 atm.
 * @param {number} T K @returns {{ A: number, B: number }}
 */
export function debyeHuckelParameters(T) {
  const eps = waterDielectric(T);
  const rho = brineDensity({ T, S: 0, P: 101325 }).rho / 1000;
  return { A: 1.82483e6 * Math.sqrt(rho) / (eps * T) ** 1.5, B: 50.2916 * Math.sqrt(rho) / Math.sqrt(eps * T) };
}

/**
 * Single-ion activity coefficient.
 *  - 'davies': log γ = −A z² [√I/(1 + √I) − 0.3 I]  (Davies 1962; usable to I ≈ 0.5 mol/kg).
 *  - 'extendedDH': log γ = −A z² √I/(1 + B å √I)  (extended Debye–Hückel; I ≲ 0.1 mol/kg).
 * Seawater (I ≈ 0.7) and formation brines (I up to several mol/kg) are beyond the formal range of both; a Pitzer
 * ion-interaction model is needed for quantitative work there. Results at high ionic strength are screening values
 * only (the Davies form is used with the ionic strength as given, no cap).
 * Reference: Davies, "Ion Association" (1962); Debye & Hückel (1923); Kielland, J. Am. Chem. Soc. 59 (1937) 1675.
 * @param {object} o
 * @param {number} o.z ion charge @param {number} o.I ionic strength, mol/kg @param {number} [o.T=298.15] K
 * @param {'davies'|'extendedDH'} [o.model='davies'] @param {number} [o.a] ion-size parameter, Å (extendedDH)
 * @returns {number} activity coefficient, –
 */
export function activityCoefficient({ z, I, T = 298.15, model = 'davies', a }) {
  need('z', z, -6, 6);
  need('I', I, 0, 20);
  const { A, B } = debyeHuckelParameters(T);
  const sq = Math.sqrt(I);
  if (model === 'davies') return 10 ** (-A * z * z * (sq / (1 + sq) - 0.3 * I));
  if (model === 'extendedDH') return 10 ** (-A * z * z * sq / (1 + B * need('a', a, 1, 15) * sq));
  throw new RangeError(`model = ${model} is not one of davies, extendedDH`);
}

const LOGK = {
  // log10 K = A + B T + C/T + D log10 T + E/T²  (T in K)
  calcite: { c: [-171.9065, -0.077993, 2839.319, 71.595, 0], range: [273.15, 363.15], M: 100.0869,
    ref: 'Plummer & Busenberg, Geochim. Cosmochim. Acta 46 (1982) 1011' },
  barite: { c: [136.035, 0, -7680.41, -48.595, 0], range: [273.15, 573.15], M: 233.39,
    ref: 'Langmuir & Melchior, Geochim. Cosmochim. Acta 49 (1985) 2423 (after Blount 1977)' },
  celestite: { c: [-14805.9622, -2.4660924, 756968.533, 5436.3588, -40553604], range: [273.15, 363.15], M: 183.68,
    ref: 'Reardon & Armstrong, Geochim. Cosmochim. Acta 51 (1987) 63' },
  gypsum: { c: [68.2401, 0, -3221.51, -25.0627, 0], range: [273.15, 383.15], M: 172.17,
    ref: 'Langmuir & Melchior (1985)' },
  anhydrite: { c: [197.52, 0, -8669.8, -69.835, 0], range: [273.15, 473.15], M: 136.14,
    ref: 'Langmuir & Melchior (1985)' },
  // carbonate system (Plummer & Busenberg 1982): Henry constant in mol/(kg·atm)
  KH: { c: [108.3865, 0.01985076, -6919.53, -40.45154, 669365], range: [273.15, 523.15] },
  K1: { c: [-356.3094, -0.06091964, 21834.37, 126.8339, -1684915], range: [273.15, 523.15] },
  K2: { c: [-107.8871, -0.03252849, 5151.79, 38.92561, -563713.9], range: [273.15, 523.15] },
};

function logKfit(key, T) {
  const k = LOGK[key];
  need('T', T, k.range[0], k.range[1] + 1e-9);
  const [A, B, C, D, E] = k.c;
  return A + B * T + C / T + D * Math.log10(T) + E / (T * T);
}

/**
 * Thermodynamic solubility product of a scale mineral, log10 K_sp(T) = A + B T + C/T + D log10 T + E/T²:
 *   calcite   CaCO3:        −171.9065 − 0.077993 T + 2839.319/T + 71.595 log T            (25 °C: −8.48)
 *   barite    BaSO4:        136.035 − 7680.41/T − 48.595 log T                           (25 °C: −9.97)
 *   celestite SrSO4:        −14805.9622 − 2.4660924 T + 756968.533/T − 40553604/T² + 5436.3588 log T   (−6.63)
 *   gypsum    CaSO4·2H2O:   68.2401 − 3221.51/T − 25.0627 log T                          (25 °C: −4.58)
 *   anhydrite CaSO4:        197.52 − 8669.8/T − 69.835 log T                             (25 °C: −4.36)
 * References: Plummer & Busenberg (1982); Langmuir & Melchior (1985); Reardon & Armstrong (1987) — the analytic
 * expressions as compiled in the USGS PHREEQC database (phreeqc.dat).
 * Validity: saturation pressure of water; calcite and celestite 0–90 °C, barite to 300 °C, gypsum to 110 °C,
 * anhydrite to 200 °C. An optional pressure correction ln(K_P/K) = −ΔV_r (P − P_0)/(R T) is applied when the
 * molar volume change of dissolution `dVr` (m³/mol, negative for these salts) is supplied.
 * @param {object} o
 * @param {'calcite'|'barite'|'celestite'|'gypsum'|'anhydrite'} o.mineral
 * @param {number} o.T K @param {number} [o.P=101325] Pa @param {number} [o.dVr] m³/mol
 * @returns {{ logK: number, K: number, molarMass: number }} molar mass of the solid in kg/mol
 */
export function solubilityProduct({ mineral, T, P = 101325, dVr }) {
  const k = LOGK[mineral];
  if (!k || !k.M) throw new RangeError(`mineral = ${mineral} is not one of calcite, barite, celestite, gypsum, anhydrite`);
  let logK = logKfit(mineral, T);
  if (dVr != null) logK += -need('dVr', dVr, -1e-3, 1e-3) * (need('P', P, 0, 2e8) - 101325) / (8.314462618 * T) / Math.LN10;
  return { logK, K: 10 ** logK, molarMass: k.M / 1000 };
}

/**
 * Carbonate-system equilibrium constants (molal scale, CO2 fugacity in atm):
 *   log K_H = 108.3865 + 0.01985076 T − 6919.53/T − 40.45154 log T + 669365/T²          (CO2(g) ⇌ CO2(aq); −1.468 at 25 °C)
 *   log K_1 = −356.3094 − 0.06091964 T + 21834.37/T + 126.8339 log T − 1684915/T²        (−6.352)
 *   log K_2 = −107.8871 − 0.03252849 T + 5151.79/T + 38.92561 log T − 563713.9/T²        (−10.329)
 *   log K_w = −4470.99/T + 6.0875 − 0.01706 T                                             (−13.995; Harned & Owen)
 * Reference: Plummer & Busenberg, Geochim. Cosmochim. Acta 46 (1982) 1011–1040; Harned & Owen, "The Physical
 * Chemistry of Electrolytic Solutions" (1958). Validity: 0–250 °C (K_w fit 0–60 °C).
 * @param {number} T K
 * @returns {{ logKH: number, logK1: number, logK2: number, logKw: number }}
 */
export function carbonateConstants(T) {
  return { logKH: logKfit('KH', T), logK1: logKfit('K1', T), logK2: logKfit('K2', T),
    logKw: -4470.99 / T + 6.0875 - 0.01706 * T };
}

/**
 * pH of water in equilibrium with a CO2 partial pressure for a given bicarbonate alkalinity, from K_1, K_H and the
 * charge balance [HCO3⁻] = Alk + [H⁺]: a_H² + γ Alk a_H − K_1 K_H p_CO2 = 0. Limits: a_H = K_1 K_H p_CO2/(γ Alk)
 * for buffered water and a_H = √(K_1 K_H p_CO2) for condensed water without alkalinity.
 * Reference: carbonate equilibria with the constants of `carbonateConstants`; CO2 treated as an ideal gas
 * (multiply p_CO2 by a fugacity coefficient for high pressure).
 * @param {object} o
 * @param {number} o.T K @param {number} o.pCO2 CO2 partial pressure, Pa
 * @param {number} [o.mHCO3=0] bicarbonate alkalinity, mol/kg @param {number} [o.I=0] ionic strength, mol/kg
 * @returns {number} pH (activity scale)
 */
export function phFromCO2({ T, pCO2, mHCO3 = 0, I = 0 }) {
  need('pCO2', pCO2, 1e-3, 1e9);
  need('mHCO3', mHCO3, 0, 10);
  const k = carbonateConstants(T);
  const aCO2 = 10 ** k.logKH * (pCO2 / 101325);
  const K1 = 10 ** k.logK1;
  const g1 = activityCoefficient({ z: 1, I, T });
  // charge balance [HCO3⁻] = Alk + [H⁺] (carbonate and hydroxide neglected) in K1 = a_H γ [HCO3⁻]/a_CO2 gives
  // a_H² + γ Alk a_H − K1 a_CO2 = 0
  const b = mHCO3 * g1;
  const aH = 0.5 * (-b + Math.sqrt(b * b + 4 * K1 * aCO2));
  return -Math.log10(aH);
}

/**
 * Saturation index and ratio of calcium carbonate (calcite): SI = log10(a_Ca a_CO3/K_sp), with
 * a_CO3 = K_2 a_HCO3/a_H. The pH is either given or computed from the CO2 partial pressure (`phFromCO2`).
 * Reference: Plummer & Busenberg (1982) constants; Davies activity coefficients. Validity: screening accuracy at
 * oilfield ionic strengths (see `activityCoefficient`); ion pairing (CaHCO3⁺, CaSO4°) is neglected.
 * @param {object} o
 * @param {Object<string, number>} o.molal water composition, mol/kg (must contain Ca and HCO3)
 * @param {number} o.T K @param {number} [o.pH] @param {number} [o.pCO2] Pa (one of pH and pCO2)
 * @returns {{ SI: number, SR: number, pH: number, I: number, logKsp: number }}
 */
export function calciteSaturation({ molal, T, pH, pCO2 }) {
  if ((pH == null) === (pCO2 == null)) throw new RangeError('give exactly one of pH and pCO2');
  if (!(molal?.Ca > 0) || !(molal?.HCO3 > 0)) throw new RangeError('molal.Ca and molal.HCO3 must be positive');
  const I = ionicStrength(molal);
  const ph = pH ?? phFromCO2({ T, pCO2, mHCO3: molal.HCO3, I });
  need('pH', ph, 0, 14);
  const k = carbonateConstants(T);
  const g1 = activityCoefficient({ z: 1, I, T });
  const g2 = activityCoefficient({ z: 2, I, T });
  const logIAP = Math.log10(g2 * molal.Ca) + k.logK2 + Math.log10(g1 * molal.HCO3) + ph;
  const logKsp = solubilityProduct({ mineral: 'calcite', T }).logK;
  return { SI: logIAP - logKsp, SR: 10 ** (logIAP - logKsp), pH: ph, I, logKsp };
}

const SULPHATES = { barite: 'Ba', celestite: 'Sr', gypsum: 'Ca', anhydrite: 'Ca' };

/**
 * Saturation index, saturation ratio and precipitation potential of a sulphate scale (barite, celestite, gypsum,
 * anhydrite): SI = log10(γ_M m_M γ_SO4 m_SO4/K_sp). The precipitation potential is the amount x that must
 * precipitate to return to saturation at constant activity coefficients:
 *   (m_M − x)(m_SO4 − x) γ_M γ_SO4 = K_sp  →  x = ½[(m_M + m_SO4) − √((m_M − m_SO4)² + 4 K_sp/(γ_M γ_SO4))].
 * Each mineral is treated independently (no competition for sulphate, water activity = 1 for gypsum).
 * Reference: `solubilityProduct`; Davies activity coefficients (screening at oilfield ionic strength).
 * @param {object} o
 * @param {'barite'|'celestite'|'gypsum'|'anhydrite'} o.mineral
 * @param {Object<string, number>} o.molal water composition, mol/kg
 * @param {number} o.T K @param {number} [o.P=101325] Pa @param {number} [o.dVr] m³/mol
 * @returns {{ SI: number, SR: number, precipitate: number, mass: number, I: number }} precipitate in mol/kg water;
 *   mass in kg of solid per kg of water (≈ per litre; multiply by the water density for kg/m³)
 */
export function sulphateSaturation({ mineral, molal, T, P = 101325, dVr }) {
  const cation = SULPHATES[mineral];
  if (!cation) throw new RangeError(`mineral = ${mineral} is not one of barite, celestite, gypsum, anhydrite`);
  const mM = molal?.[cation] ?? 0;
  const mS = molal?.SO4 ?? 0;
  const I = ionicStrength(molal);
  const ksp = solubilityProduct({ mineral, T, P, dVr });
  const g2 = activityCoefficient({ z: 2, I, T });
  if (!(mM > 0) || !(mS > 0)) return { SI: -Infinity, SR: 0, precipitate: 0, mass: 0, I };
  const SR = g2 * mM * g2 * mS / ksp.K;
  const x = Math.max(0, 0.5 * ((mM + mS) - Math.sqrt((mM - mS) ** 2 + 4 * ksp.K / (g2 * g2))));
  return { SI: Math.log10(SR), SR, precipitate: x, mass: x * ksp.molarMass, I };
}

/**
 * Mix two waters on a mass-of-water basis (conservative ions): m_mix = f m_a + (1 − f) m_b.
 * Reference: mass balance. @param {Object<string, number>} a molalities of water A
 * @param {Object<string, number>} b molalities of water B @param {number} f mass fraction of water A, –
 * @returns {Object<string, number>} mixed molalities
 */
export function mixWaters(a, b, f) {
  need('f', f, 0, 1);
  const out = {};
  for (const ion of new Set([...Object.keys(a), ...Object.keys(b)])) out[ion] = f * (a[ion] ?? 0) + (1 - f) * (b[ion] ?? 0);
  return out;
}

/**
 * Sulphate-scaling tendency when injected seawater mixes with formation water, over the whole range of mixing
 * ratios: saturation index and precipitation potential of barite, celestite and gypsum/anhydrite at each ratio,
 * and the worst ratio (largest precipitated mass) per mineral.
 * Reference: mixing-line method (e.g. Oddo & Tomson practice; Jordan & Mackay, SPE scale papers) with the
 * equilibrium data of `solubilityProduct`.
 * Validity: equilibrium (no kinetics), independent minerals, Davies activity model — a screening calculation; use
 * a Pitzer-based scale package for design decisions.
 * @param {object} o
 * @param {Object<string, number>} o.seawater molalities @param {Object<string, number>} o.formationWater molalities
 * @param {number} o.T K @param {number} [o.P=101325] Pa @param {number} [o.n=51] number of mixing ratios
 * @param {string[]} [o.minerals=['barite','celestite','anhydrite']]
 * @param {number} [o.rhoWater=1000] kg of water per m³ of mixed brine, for the mass potential
 * @returns {{ fraction: number[], results: Object<string, { SI: number[], mass: number[], worstFraction: number,
 *   worstMass: number, maxSI: number }> }} fraction = seawater mass fraction; mass in kg solid per m³ of brine
 */
export function seawaterMixingScale({ seawater, formationWater, T, P = 101325, n = 51,
  minerals = ['barite', 'celestite', 'anhydrite'], rhoWater = 1000 }) {
  need('n', n, 3, 2001);
  need('rhoWater', rhoWater, 500, 1300);
  const fraction = Array.from({ length: n }, (_, i) => i / (n - 1));
  const results = {};
  for (const mineral of minerals) {
    const SI = [];
    const mass = [];
    let worst = 0;
    for (let i = 0; i < n; i++) {
      const r = sulphateSaturation({ mineral, molal: mixWaters(seawater, formationWater, fraction[i]), T, P });
      SI.push(r.SI);
      mass.push(r.mass * rhoWater);
      if (mass[i] > mass[worst]) worst = i;
    }
    results[mineral] = { SI, mass, worstFraction: fraction[worst], worstMass: mass[worst],
      maxSI: Math.max(...SI.filter(Number.isFinite), -Infinity) };
  }
  return { fraction, results };
}

/**
 * Langelier saturation index and Ryznar stability index:
 *   pH_s = (9.3 + A + B) − (C + D),  A = (log10 TDS − 1)/10,  B = −13.12 log10(T) + 34.55,
 *   C = log10(Ca as mg/L CaCO3) − 0.4,  D = log10(alkalinity as mg/L CaCO3);
 *   LSI = pH − pH_s (positive: scale-forming),  RSI = 2 pH_s − pH (< 6 scale-forming, > 7 corrosive).
 * Reference: Langelier, J. AWWA 28 (1936) 1500; Ryznar, J. AWWA 36 (1944) 472; the algebraic form of the Langelier
 * tables is the usual one of water-treatment practice. Validity: fresh to brackish water, TDS < ~10 000 mg/L,
 * 0–90 °C; use Stiff–Davis for brines.
 * @param {object} o
 * @param {number} o.pH @param {number} o.T K @param {number} o.tds total dissolved solids, mg/L
 * @param {number} o.calcium calcium hardness, mg/L as CaCO3 @param {number} o.alkalinity mg/L as CaCO3
 * @returns {{ pHs: number, LSI: number, RSI: number }}
 */
export function langelierRyznar({ pH, T, tds, calcium, alkalinity }) {
  need('pH', pH, 0, 14);
  need('T', T, 273.15, 373.15);
  need('tds', tds, 1, 1e5);
  need('calcium', calcium, 1e-3, 1e6);
  need('alkalinity', alkalinity, 1e-3, 1e6);
  const A = (Math.log10(tds) - 1) / 10;
  const B = -13.12 * Math.log10(T) + 34.55;
  const C = Math.log10(calcium) - 0.4;
  const D = Math.log10(alkalinity);
  const pHs = 9.3 + A + B - (C + D);
  return { pHs, LSI: pH - pHs, RSI: 2 * pHs - pH };
}

/**
 * Stiff–Davis stability index for oilfield brines: SDI = pH − pCa − pAlk − K, with pCa = −log10[Ca²⁺ mol/L],
 * pAlk = −log10[total alkalinity eq/L] and K an empirical constant read from the Stiff–Davis chart as a function
 * of ionic strength and temperature (argument — read it from the published chart).
 * Reference: Stiff & Davis, Trans. AIME 195 (1952) 213–216.
 * @param {object} o
 * @param {number} o.pH @param {number} o.Ca calcium, mol/L @param {number} o.alkalinity eq/L @param {number} o.K
 * @returns {{ SDI: number, pHs: number }} positive SDI: CaCO3 scale-forming
 */
export function stiffDavis({ pH, Ca, alkalinity, K }) {
  need('pH', pH, 0, 14);
  need('Ca', Ca, 1e-9, 10);
  need('alkalinity', alkalinity, 1e-9, 10);
  need('K', K, 0, 6);
  const pHs = K - Math.log10(Ca) - Math.log10(alkalinity);
  return { SDI: pH - pHs, pHs };
}

/**
 * Scale deposition-rate screen from supersaturation with surface-reaction kinetics of order n:
 * R = k (SR^(1/ν) − 1)^n for SR > 1 (kg/(m²·s)), ν = number of ions per formula unit (2 for these salts);
 * thickness growth = R/ρ_scale. The rate constant k is an argument (strongly temperature-, surface- and
 * inhibitor-dependent). With an inhibitor, the rate is set to zero when the dose is at or above the minimum
 * inhibitor concentration (MIC) supplied by the user.
 * Reference: surface-controlled crystal growth law (Nancollas & Reddy, J. Colloid Interface Sci. 37 (1971) 824:
 * second order for CaCO3; n = 1 gives the first-order screen).
 * @param {object} o
 * @param {number} o.SR saturation ratio IAP/K_sp @param {number} o.k rate constant, kg/(m²·s)
 * @param {number} [o.order=1] @param {number} [o.nu=2] @param {number} [o.rhoScale=2710] kg/m³ (calcite 2710, barite 4480)
 * @param {number} [o.inhibitorDose] mg/L @param {number} [o.mic] minimum inhibitor concentration, mg/L
 * @returns {{ rate: number, thicknessRate: number, inhibited: boolean }} kg/(m²·s), m/s
 */
export function scaleDepositionRate({ SR, k, order = 1, nu = 2, rhoScale = 2710, inhibitorDose, mic }) {
  need('SR', SR, 0, 1e12);
  need('k', k, 0, 1e3);
  need('order', order, 0.5, 5);
  need('nu', nu, 1, 6);
  need('rhoScale', rhoScale, 500, 1e4);
  const inhibited = inhibitorDose != null && mic != null && inhibitorDose >= mic;
  const rate = SR > 1 && !inhibited ? k * (SR ** (1 / nu) - 1) ** order : 0;
  return { rate, thicknessRate: rate / rhoScale, inhibited };
}
