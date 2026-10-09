// Internal corrosion screening: CO2 corrosion (de Waard family), inhibitor availability, corrosion allowance,
// sour-service and top-of-line screening, oxygen corrosion by mass-transfer limit, and electrochemical kinetics
// (Butler–Volmer, Tafel, mixed potential, Faraday conversion).
// Units: T in K, partial pressures and total pressure in Pa (the published correlations use bar internally),
// corrosion rates in mm/year, current density A/m², potentials V.

import { brent } from '../math/roots.js';

const R = 8.314462618;
const F = 96485.33212;
const SEC_PER_YEAR = 365.25 * 86400;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * de Waard & Milliams CO2 corrosion rate of carbon steel (equation of the 1975 nomogram):
 * log10 V = 5.8 − 1710/T + 0.67 log10 p_CO2  (V mm/year, T K, p_CO2 bar).
 * Reference: de Waard & Milliams, Corrosion 31(5) (1975) 177–181. Validity: wet CO2, no protective film,
 * T up to about 60–70 °C (above the scaling temperature the rate is over-predicted), p_CO2 up to ~10 bar;
 * a worst-case "nomogram" rate before correction factors.
 * @param {object} o
 * @param {number} o.T temperature, K @param {number} o.pCO2 CO2 partial pressure (or fugacity), Pa
 * @returns {number} corrosion rate, mm/year
 */
export function deWaardMilliams({ T, pCO2 }) {
  need('T', T, 273.15, 473.15);
  need('pCO2', pCO2, 1, 1e8);
  return 10 ** (5.8 - 1710 / T + 0.67 * Math.log10(pCO2 / 1e5));
}

/**
 * Fugacity coefficient of CO2 in the gas phase used with the de Waard correlations:
 * log10 φ = (0.0031 − 1.4/T) P  (P total pressure in bar, T in K).
 * Reference: de Waard, Lotz & Milliams, Corrosion 47(12) (1991) 976–985, as revised in de Waard & Lotz,
 * NACE Corrosion/93 paper 69. Validity: up to about 250 bar.
 * @param {object} o
 * @param {number} o.T K @param {number} o.P total pressure, Pa
 * @returns {number} fugacity coefficient, –
 */
export function co2FugacityCoefficient({ T, P }) {
  need('T', T, 273.15, 473.15);
  need('P', P, 0, 5e7);
  return 10 ** ((0.0031 - 1.4 / T) * (P / 1e5));
}

/**
 * pH of pure water saturated with CO2 (condensed water, no dissolved iron or bicarbonate).
 *  - '1991': pH = 3.71 + 0.00417 t − 0.5 log10 p_CO2      (de Waard, Lotz & Milliams 1991)
 *  - '1995': pH = 3.82 + 0.00384 t − 0.5 log10 f_CO2      (de Waard, Lotz & Dugstad 1995)
 * (t in °C, pressure in bar.)
 * Reference: as cited. Validity: about 5–90 °C, p_CO2 up to ~10 bar.
 * @param {object} o
 * @param {number} o.T K @param {number} o.pCO2 CO2 partial pressure or fugacity, Pa
 * @param {'1991'|'1995'} [o.version='1995']
 * @returns {number} pH
 */
export function phCO2Saturated({ T, pCO2, version = '1995' }) {
  need('T', T, 273.15, 473.15);
  need('pCO2', pCO2, 1, 1e8);
  const t = T - 273.15;
  const lp = Math.log10(pCO2 / 1e5);
  if (version === '1991') return 3.71 + 0.00417 * t - 0.5 * lp;
  if (version === '1995') return 3.82 + 0.00384 * t - 0.5 * lp;
  throw new RangeError(`version = ${version} is not one of 1991, 1995`);
}

/**
 * de Waard, Lotz & Milliams (1991) CO2 corrosion rate with correction factors:
 *   V = V_nomogram(f_CO2) · F_scale · F_pH · F_glycol,
 *   V_nomogram: log V = 5.8 − 1710/T + 0.67 log f_CO2,  f_CO2 = φ p_CO2 (`co2FugacityCoefficient`);
 *   protective scale (high temperature): log F_scale = 2400/T − 0.6 log f_CO2 − 6.7, capped at 1
 *     (i.e. effective above the scaling temperature T_scale = 2400/(6.7 + 0.6 log f_CO2));
 *   pH: with pH_sat = min(1.36 + 1307/T − 0.17 log f_CO2, 5.4 − 0.66 log f_CO2) (saturation with FeCO3 / Fe3O4),
 *     log F_pH = 0.32 (pH_sat − pH_act) for pH_sat > pH_act, log F_pH = −0.13 (pH_act − pH_sat)^1.6 otherwise;
 *     applied only when the actual pH is given;
 *   glycol: log F_glycol = 1.6 (log W − 2), W = water content of the water–glycol mixture in wt %.
 * Reference: de Waard, Lotz & Milliams, Corrosion 47(12) (1991) 976–985. Validity: carbon steel, sweet service,
 * T ≤ ~140 °C, f_CO2 ≤ ~10 bar; factors for crude-oil wetting and condensation are left to the user
 * (`extraFactor`).
 * @param {object} o
 * @param {number} o.T K @param {number} o.pCO2 CO2 partial pressure, Pa
 * @param {number} [o.P] total pressure, Pa (for the fugacity coefficient; φ = 1 when omitted)
 * @param {number} [o.pH] actual pH of the water phase
 * @param {number} [o.waterInGlycol] water content of the aqueous glycol phase, wt % (100 = no glycol)
 * @param {number} [o.extraFactor=1] additional multiplicative factor
 * @returns {{ rate: number, nomogram: number, fCO2: number, Fscale: number, FpH: number, Fglycol: number,
 *   Tscale: number, pHsat: number }} rates mm/year, fugacity Pa, scaling temperature K
 */
export function deWaardLotzMilliams({ T, pCO2, P, pH, waterInGlycol, extraFactor = 1 }) {
  need('T', T, 273.15, 473.15);
  need('pCO2', pCO2, 1, 1e8);
  need('extraFactor', extraFactor, 0, 10);
  const phi = P != null ? co2FugacityCoefficient({ T, P }) : 1;
  const fbar = phi * pCO2 / 1e5;
  const lf = Math.log10(fbar);
  const nomogram = 10 ** (5.8 - 1710 / T + 0.67 * lf);
  const Fscale = Math.min(1, 10 ** (2400 / T - 0.6 * lf - 6.7));
  const pHsat = Math.min(1.36 + 1307 / T - 0.17 * lf, 5.4 - 0.66 * lf);
  let FpH = 1;
  if (pH != null) {
    need('pH', pH, 0, 14);
    FpH = pHsat > pH ? 10 ** (0.32 * (pHsat - pH)) : 10 ** (-0.13 * (pH - pHsat) ** 1.6);
  }
  let Fglycol = 1;
  if (waterInGlycol != null) {
    need('waterInGlycol', waterInGlycol, 0.1, 100);
    Fglycol = 10 ** (1.6 * (Math.log10(waterInGlycol) - 2));
  }
  return { rate: nomogram * Fscale * FpH * Fglycol * extraFactor, nomogram, fCO2: fbar * 1e5, Fscale, FpH, Fglycol,
    Tscale: 2400 / (6.7 + 0.6 * lf), pHsat };
}

/**
 * de Waard, Lotz & Dugstad (1995) resistance model with flow dependence:
 *   1/V = 1/V_r + 1/V_m,
 *   reaction:      log V_r = 4.93 − 1119/T + 0.58 log f_CO2 − 0.34 (pH_act − pH_CO2),
 *   mass transfer: V_m = 2.45 U^0.8/d^0.2 · f_CO2,
 * with f_CO2 in bar, liquid velocity U in m/s, hydraulic diameter d in m, pH_CO2 = 3.82 + 0.00384 t − 0.5 log f_CO2
 * (pure CO2-saturated water). The form quoted here is the base equation (best-fit constants for normalised
 * St 52 steel); the steel-composition (Cr, C) factors of the paper are not included. Protective-scale and glycol
 * factors of the 1991 model may be applied on top through `extraFactor`.
 * Reference: de Waard, Lotz & Dugstad, NACE Corrosion/95 paper 128.
 * Validity: 20–90 °C, U 1.5–13 m/s in the test data, f_CO2 0.3–20 bar, no protective film.
 * @param {object} o
 * @param {number} o.T K @param {number} o.pCO2 CO2 partial pressure, Pa
 * @param {number} o.U liquid velocity, m/s @param {number} o.d hydraulic diameter, m
 * @param {number} [o.P] total pressure, Pa (fugacity coefficient; φ = 1 when omitted)
 * @param {number} [o.pH] actual pH (default: pH of CO2-saturated water, i.e. no pH correction)
 * @param {number} [o.extraFactor=1]
 * @returns {{ rate: number, Vr: number, Vm: number, pHCO2: number, fCO2: number }} mm/year; fugacity Pa
 */
export function deWaardLotzDugstad({ T, pCO2, U, d, P, pH, extraFactor = 1 }) {
  need('T', T, 273.15, 473.15);
  need('pCO2', pCO2, 1, 1e8);
  need('U', U, 1e-6, 100);
  need('d', d, 1e-4, 10);
  need('extraFactor', extraFactor, 0, 10);
  const phi = P != null ? co2FugacityCoefficient({ T, P }) : 1;
  const fbar = phi * pCO2 / 1e5;
  const pHCO2 = phCO2Saturated({ T, pCO2: fbar * 1e5, version: '1995' });
  const pHact = pH ?? pHCO2;
  need('pH', pHact, 0, 14);
  const Vr = 10 ** (4.93 - 1119 / T + 0.58 * Math.log10(fbar) - 0.34 * (pHact - pHCO2));
  const Vm = 2.45 * U ** 0.8 / d ** 0.2 * fbar;
  return { rate: extraFactor / (1 / Vr + 1 / Vm), Vr, Vm, pHCO2, fCO2: fbar * 1e5 };
}

/**
 * Inhibited corrosion rate by the inhibitor-availability model:
 * CR = A · CR_inhibited + (1 − A) · CR_uninhibited, A = fraction of time the inhibitor is available at the
 * required dose.
 * Reference: NORSOK M-001 (materials selection), inhibitor-availability concept. Typical design values:
 * A = 0.90–0.95 (0.99 with redundancy and monitoring), CR_inhibited ≈ 0.1 mm/year — both are arguments.
 * @param {object} o
 * @param {number} o.uninhibited mm/year @param {number} [o.inhibited=0.1] mm/year @param {number} [o.availability=0.9] –
 * @returns {number} time-averaged corrosion rate, mm/year
 */
export function inhibitedRate({ uninhibited, inhibited = 0.1, availability = 0.9 }) {
  need('uninhibited', uninhibited, 0, 1e3);
  need('inhibited', inhibited, 0, 1e3);
  need('availability', availability, 0, 1);
  return availability * inhibited + (1 - availability) * uninhibited;
}

/**
 * Corrosion allowance and remaining life: CA = CR · life; life = (allowance)/CR.
 * Reference: definition; many company specifications limit the allowance for carbon steel to 6–10 mm, above which
 * a corrosion-resistant alloy or cladding is selected (limit taken as an argument).
 * @param {object} o
 * @param {number} o.rate corrosion rate, mm/year @param {number} o.life design life, years
 * @param {number} [o.maxAllowance=8] largest practicable allowance, mm
 * @param {number} [o.available] allowance actually provided, mm (to compute the achieved life)
 * @returns {{ required: number, feasible: boolean, lifeAchieved: number|null }} mm, –, years
 */
export function corrosionAllowance({ rate, life, maxAllowance = 8, available }) {
  need('rate', rate, 0, 1e3);
  need('life', life, 0, 200);
  need('maxAllowance', maxAllowance, 0, 100);
  const required = rate * life;
  return { required, feasible: required <= maxAllowance,
    lifeAchieved: available != null ? (rate > 0 ? need('available', available, 0, 100) / rate : Infinity) : null };
}

/**
 * Sour-service and corrosion-regime screening.
 * Sour service for sulphide stress cracking: p_H2S ≥ 0.3 kPa (0.05 psi) per ISO 15156-2 / NACE MR0175 (clause 7).
 * Corrosion regime from the CO2/H2S partial-pressure ratio (rule of thumb after Pots et al., NACE Corrosion/2002
 * paper 02235): ratio > 500 CO2-dominated (sweet), 20–500 mixed, < 20 H2S-dominated (sour, iron-sulphide films).
 * The CO2 correlations in this file are not valid in the sour regime.
 * @param {object} o
 * @param {number} o.pH2S Pa @param {number} o.pCO2 Pa
 * @returns {{ sourService: boolean, ratio: number, regime: 'sweet'|'mixed'|'sour'|'none' }}
 */
export function sourServiceScreening({ pH2S, pCO2 }) {
  need('pH2S', pH2S, 0, 1e8);
  need('pCO2', pCO2, 0, 1e8);
  const ratio = pH2S > 0 ? pCO2 / pH2S : Infinity;
  const regime = pH2S === 0 && pCO2 === 0 ? 'none' : ratio > 500 ? 'sweet' : ratio >= 20 ? 'mixed' : 'sour';
  return { sourService: pH2S >= 300, ratio, regime };
}

/**
 * Top-of-line corrosion screening for wet-gas lines: TLC is considered likely when the water condensation rate on
 * the upper pipe wall exceeds a critical value. The condensation rate is system-specific (heat loss through the
 * wall and water content of the gas) and is an input; the threshold is an argument with the commonly quoted range 0.15–0.25 g/(m²·s) (Gunaltun, Supriyatman & Achmad, NACE Corrosion/99 paper 36).
 * @param {object} o
 * @param {number} o.condensationRate kg/(m²·s) @param {number} [o.threshold=2.5e-4] kg/(m²·s)
 * @returns {{ atRisk: boolean, ratio: number }}
 */
export function topOfLineScreening({ condensationRate, threshold = 2.5e-4 }) {
  need('condensationRate', condensationRate, 0, 1);
  need('threshold', threshold, 1e-9, 1);
  return { atRisk: condensationRate >= threshold, ratio: condensationRate / threshold };
}

/**
 * Mass-transfer coefficient in turbulent pipe flow by Berger & Hau (1977): Sh = 0.0165 Re^0.86 Sc^0.33,
 * k_m = Sh D_diff/d.
 * Reference: Berger & Hau, Int. J. Heat Mass Transfer 20 (1977) 1185–1194. Validity: 8e3 < Re < 2e5,
 * 1000 < Sc < 6000, smooth straight pipe.
 * @param {object} o
 * @param {number} o.Re – @param {number} o.Sc Schmidt number μ/(ρ D), – @param {number} o.Ddiff m²/s
 * @param {number} o.d pipe diameter, m
 * @returns {{ km: number, Sh: number }} m/s
 */
export function massTransferBergerHau({ Re, Sc, Ddiff, d }) {
  need('Re', Re, 1, 1e8);
  need('Sc', Sc, 1e-2, 1e6);
  need('Ddiff', Ddiff, 1e-12, 1e-6);
  need('d', d, 1e-4, 10);
  const Sh = 0.0165 * Re ** 0.86 * Sc ** 0.33;
  return { km: Sh * Ddiff / d, Sh };
}

/**
 * Oxygen corrosion rate of carbon steel under mass-transfer (diffusion) control of oxygen reduction:
 * limiting current i_lim = 4 F k_m c_O2 (O2 + 2H2O + 4e⁻ → 4OH⁻), corrosion rate from Faraday's law for Fe → Fe²⁺.
 * Equivalent: 2 mol Fe per mol O2 → CR = 2 k_m c_O2 M_Fe/ρ_Fe.
 * Reference: diffusion-limited cathodic reaction (e.g. Fontana, "Corrosion Engineering" ch. 9); the oxygen
 * concentration and mass-transfer coefficient are inputs (k_m e.g. from `massTransferBergerHau`; Oldfield &
 * Todd and Berger–Hau type models differ only in how k_m is obtained).
 * @param {object} o
 * @param {number} o.cO2 dissolved oxygen, kg/m³ (1 ppm by mass in water ≈ 1e-3 kg/m³)
 * @param {number} o.km mass-transfer coefficient, m/s
 * @returns {{ rate: number, iLim: number }} mm/year and A/m²
 */
export function oxygenCorrosionRate({ cO2, km }) {
  need('cO2', cO2, 0, 1);
  need('km', km, 0, 1);
  const iLim = 4 * F * km * cO2 / 0.0319988;
  return { rate: currentToRate({ i: iLim }), iLim };
}

/**
 * Faraday conversion of a corrosion current density to a penetration rate: CR = i M/(n F ρ).
 * For iron (M = 55.845 g/mol, n = 2, ρ = 7870 kg/m³): 1 A/m² = 1.16 mm/year.
 * Reference: Faraday's laws; ASTM G102.
 * @param {object} o
 * @param {number} o.i current density, A/m²
 * @param {number} [o.M=0.055845] molar mass, kg/mol @param {number} [o.n=2] electrons per atom
 * @param {number} [o.rho=7870] density, kg/m³
 * @returns {number} mm/year
 */
export function currentToRate({ i, M = 0.055845, n = 2, rho = 7870 }) {
  need('i', i, 0, 1e7);
  need('M', M, 1e-3, 1);
  need('n', n, 1, 8);
  need('rho', rho, 100, 3e4);
  return i * M / (n * F * rho) * SEC_PER_YEAR * 1e3;
}

/**
 * Butler–Volmer current density of a single electrode reaction:
 * i = i_0 {exp[α_a n F η/(R T)] − exp[−α_c n F η/(R T)]}, η = E − E_eq.
 * Reference: Bard & Faulkner, "Electrochemical Methods" (2001) ch. 3.
 * @param {object} o
 * @param {number} o.i0 exchange current density, A/m² @param {number} o.eta overpotential, V
 * @param {number} [o.alphaA=0.5] @param {number} [o.alphaC=0.5] transfer coefficients
 * @param {number} [o.n=1] electrons in the rate-determining step @param {number} [o.T=298.15] K
 * @returns {number} net (anodic-positive) current density, A/m²
 */
export function butlerVolmer({ i0, eta, alphaA = 0.5, alphaC = 0.5, n = 1, T = 298.15 }) {
  need('i0', i0, 0, 1e9);
  need('eta', eta, -5, 5);
  need('alphaA', alphaA, 0, 4);
  need('alphaC', alphaC, 0, 4);
  need('T', T, 200, 800);
  const k = n * F / (R * T);
  return i0 * (Math.exp(alphaA * k * eta) - Math.exp(-alphaC * k * eta));
}

/**
 * Tafel slope b = 2.303 R T/(α n F) (V per decade) and the Tafel overpotential η = b log10(i/i_0).
 * At 25 °C and α n = 0.5: b = 0.1183 V/decade.
 * Reference: Tafel (1905); high-overpotential limit of the Butler–Volmer equation.
 * @param {object} o
 * @param {number} [o.alpha=0.5] @param {number} [o.n=1] @param {number} [o.T=298.15] K
 * @param {number} [o.i] current density, A/m² @param {number} [o.i0] exchange current density, A/m²
 * @returns {{ b: number, eta: number|null }} V/decade, V
 */
export function tafel({ alpha = 0.5, n = 1, T = 298.15, i, i0 }) {
  need('alpha', alpha, 1e-3, 4);
  need('T', T, 200, 800);
  const b = Math.LN10 * R * T / (alpha * n * F);
  return { b, eta: i != null && i0 != null ? b * Math.log10(need('i', i, 1e-30, 1e9) / need('i0', i0, 1e-30, 1e9)) : null };
}

function anodicCurrent(a, E) {
  return a.i0 * 10 ** ((E - a.E0) / a.b);
}

function cathodicCurrent(c, E) {
  // activation branch in series with the diffusion limit: 1/i = 1/i_act + 1/i_lim
  const act = c.i0 * 10 ** (-(E - c.E0) / c.b);
  return c.iLim != null && Number.isFinite(c.iLim) ? 1 / (1 / act + 1 / c.iLim) : act;
}

function checkBranch(name, r, lim) {
  if (!r) throw new RangeError(`${name} reaction parameters are required`);
  need(`${name}.E0`, r.E0, -5, 5);
  need(`${name}.i0`, r.i0, 1e-30, 1e9);
  need(`${name}.b`, r.b, 1e-4, 5);
  if (lim && r.iLim != null) need(`${name}.iLim`, r.iLim, 1e-30, Infinity);
}

/**
 * Mixed-potential (corrosion) solution for one anodic and one cathodic reaction:
 *   anodic (metal dissolution, Tafel):   i_a = i_0a · 10^((E − E_a)/b_a)
 *   cathodic (Tafel with diffusion limit): 1/i_c = 1/[i_0c · 10^(−(E − E_c)/b_c)] + 1/i_lim
 * The corrosion potential satisfies i_a(E_corr) = i_c(E_corr); solved by Brent between E_a − 1 V and E_c + 1 V.
 * Without a diffusion limit the closed form is E_corr = [E_a/b_a + E_c/b_c + log10(i_0c/i_0a)]/(1/b_a + 1/b_c);
 * with i_lim far below the activation current the result tends to i_corr = i_lim.
 * Reference: Wagner & Traud, Z. Elektrochem. 44 (1938) 391; Stern & Geary, J. Electrochem. Soc. 104 (1957) 56.
 * Validity: both reactions far from their own equilibrium potentials (reverse reactions neglected), uniform
 * corrosion, no ohmic drop.
 * @param {object} o
 * @param {{ E0: number, i0: number, b: number }} o.anodic equilibrium potential V, exchange current density A/m²,
 *   Tafel slope V/decade
 * @param {{ E0: number, i0: number, b: number, iLim?: number }} o.cathodic as above plus limiting current A/m²
 * @returns {{ Ecorr: number, icorr: number, rate: number }} V, A/m², mm/year for iron (see `currentToRate`)
 */
export function mixedPotential({ anodic, cathodic }) {
  checkBranch('anodic', anodic, false);
  checkBranch('cathodic', cathodic, true);
  if (!(cathodic.E0 > anodic.E0)) throw new RangeError('cathodic.E0 must exceed anodic.E0 for spontaneous corrosion');
  const f = (E) => Math.log(anodicCurrent(anodic, E)) - Math.log(cathodicCurrent(cathodic, E));
  const Ecorr = brent(f, anodic.E0 - 1, cathodic.E0 + 1, { tol: 1e-14 });
  const icorr = anodicCurrent(anodic, Ecorr);
  return { Ecorr, icorr, rate: currentToRate({ i: icorr }) };
}

/**
 * Polarisation curve (Evans diagram data) for the two-reaction system of `mixedPotential`.
 * Reference: as `mixedPotential`.
 * @param {object} o
 * @param {{ E0: number, i0: number, b: number }} o.anodic @param {{ E0: number, i0: number, b: number, iLim?: number }} o.cathodic
 * @param {number} [o.Emin] V (default E_corr − 0.3) @param {number} [o.Emax] V (default E_corr + 0.3)
 * @param {number} [o.n=121] points
 * @returns {{ E: number[], iAnodic: number[], iCathodic: number[], iNet: number[], Ecorr: number, icorr: number }}
 *   potentials V; partial and net (anodic-positive) current densities A/m²
 */
export function polarisationCurve({ anodic, cathodic, Emin, Emax, n = 121 }) {
  const mp = mixedPotential({ anodic, cathodic });
  const lo = Emin ?? mp.Ecorr - 0.3;
  const hi = Emax ?? mp.Ecorr + 0.3;
  need('n', n, 2, 1e5);
  if (!(hi > lo)) throw new RangeError('Emax must exceed Emin');
  const E = [];
  const ia = [];
  const ic = [];
  const net = [];
  for (let k = 0; k < n; k++) {
    const e = lo + (hi - lo) * k / (n - 1);
    E.push(e);
    ia.push(anodicCurrent(anodic, e));
    ic.push(cathodicCurrent(cathodic, e));
    net.push(ia[k] - ic[k]);
  }
  return { E, iAnodic: ia, iCathodic: ic, iNet: net, Ecorr: mp.Ecorr, icorr: mp.icorr };
}
