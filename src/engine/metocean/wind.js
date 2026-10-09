// Wind models: mean profiles, the ISO 19901-1 / NORSOK N-003 (Frøya) profile, gust and spectrum, IEC 61400-1
// turbulence and extreme-wind models, wind force on bodies, air density and Weibull wind-climate helpers.
//
// Conventions: z is height above mean sea level (m), wind speeds in m/s, frequency f in Hz, wind spectra are
// one-sided densities in (m/s)²/Hz.

import { gamma as gammaFn } from '../math/special.js';

/** Dry-air density at 15 °C and 101325 Pa (ISO standard atmosphere), kg/m³. */
export const RHO_AIR = 1.225;
/** von Kármán constant, -. */
export const KARMAN = 0.4;

function reqPos(name, v) {
  if (typeof v !== 'number' || !(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a positive finite number`);
}
function reqNonNeg(name, v) {
  if (typeof v !== 'number' || !(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a non-negative finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Mean profiles                                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Power-law wind profile, U(z) = Uref (z/zRef)^α.
 * Reference: DNV-RP-C205 (2019) §2.3.2.11 (α ≈ 0.12 open sea, 0.14–0.16 open country); IEC 61400-1 (2019) §6.3.2.2
 * normal wind profile with α = 0.2 (IEC 61400-3-1 uses 0.14 offshore).
 * Validity: neutral atmosphere, heights within the surface layer (up to roughly 100–200 m).
 * @param {number} z height, m
 * @param {number} Uref wind speed at zRef, m/s
 * @param {number} [zRef=10] reference height, m
 * @param {number} [alpha=0.12] shear exponent, -
 * @returns {number} wind speed at z, m/s
 */
export function powerLawProfile(z, Uref, zRef = 10, alpha = 0.12) {
  reqPos('z', z); reqNonNeg('Uref', Uref); reqPos('zRef', zRef);
  return Uref * Math.pow(z / zRef, alpha);
}

/**
 * Logarithmic wind profile, U(z) = Uref ln(z/z0) / ln(zRef/z0) = (u* / κ) ln(z/z0), κ = 0.4.
 * Reference: DNV-RP-C205 (2019) §2.3.2.4–2.3.2.6 (z0 ≈ 0.0001–0.01 m for open sea).
 * Validity: neutral stratification, z > z0, surface layer.
 * @param {number} z height, m
 * @param {number} Uref wind speed at zRef, m/s
 * @param {number} [zRef=10] reference height, m
 * @param {number} [z0=0.0002] roughness length, m
 * @returns {number} wind speed at z, m/s
 */
export function logProfile(z, Uref, zRef = 10, z0 = 0.0002) {
  reqPos('z', z); reqNonNeg('Uref', Uref); reqPos('zRef', zRef); reqPos('z0', z0);
  if (!(z > z0) || !(zRef > z0)) throw new RangeError('z and zRef must exceed z0');
  return (Uref * Math.log(z / z0)) / Math.log(zRef / z0);
}

/**
 * Friction velocity of the logarithmic profile, u* = κ U(z) / ln(z/z0).
 * Reference: DNV-RP-C205 (2019) §2.3.2.4.
 * @param {number} U wind speed at height z, m/s
 * @param {number} z height, m
 * @param {number} z0 roughness length, m
 * @returns {number} friction velocity, m/s
 */
export function frictionVelocity(U, z, z0) {
  reqNonNeg('U', U); reqPos('z', z); reqPos('z0', z0);
  if (!(z > z0)) throw new RangeError('z must exceed z0');
  return (KARMAN * U) / Math.log(z / z0);
}

/**
 * Shear exponent of the power law equivalent to a logarithmic profile between two heights,
 * α = ln[ln(z2/z0)/ln(z1/z0)] / ln(z2/z1).
 * Reference: DNV-RP-C205 (2019) §2.3.2.12.
 * @param {number} z1 lower height, m
 * @param {number} z2 upper height, m
 * @param {number} z0 roughness length, m
 * @returns {number} α, -
 */
export function equivalentShearExponent(z1, z2, z0) {
  reqPos('z1', z1); reqPos('z2', z2); reqPos('z0', z0);
  if (!(z1 > z0) || !(z2 > z1)) throw new RangeError('heights must satisfy z0 < z1 < z2');
  return Math.log(Math.log(z2 / z0) / Math.log(z1 / z0)) / Math.log(z2 / z1);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* ISO 19901-1 / NORSOK N-003 (Frøya) wind                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

const T0 = 3600;

/**
 * Frøya (NPD) wind profile with gust-averaging conversion for strong maritime winds:
 *   U(z) = U0 [1 + C ln(z/10)],        C = 0.0573 √(1 + 0.15 U0)
 *   Iu(z) = 0.06 [1 + 0.043 U0] (z/10)^(−0.22)
 *   u(z, t) = U(z) [1 − 0.41 Iu(z) ln(t/t0)],   t ≤ t0 = 3600 s
 * U0 is the 1-hour mean wind speed at 10 m; u(z, t) is the speed averaged over t with the same annual probability
 * of exceedance.
 * Reference: ISO 19901-1 (2015) §A.7.3, Eqs. A.4–A.7; NORSOK N-003 (2017) §6.3.2; Andersen & Løvseth (1992/2006);
 * DNV-RP-C205 (2019) §2.3.2.13.
 * Validity: near-neutral strong winds over open sea (U0 above about 10 m/s), 3 s ≤ t ≤ 3600 s, 10 m ≤ z ≲ 100 m.
 * @param {object} o
 * @param {number} o.U0 1-hour mean wind speed at 10 m, m/s
 * @param {number} [o.z=10] height above sea level, m
 * @param {number} [o.t=3600] averaging period, s (≤ 3600)
 * @returns {{u:number,Uz:number,Iu:number,C:number,gustFactor:number,sigma:number}} u speed averaged over t at z,
 *   m/s; Uz 1-hour mean at z, m/s; Iu turbulence intensity (-); C profile coefficient (-); gustFactor = u/Uz (-);
 *   sigma = Iu·Uz standard deviation of the along-wind speed, m/s.
 */
export function froyaWind({ U0, z = 10, t = T0 }) {
  reqNonNeg('U0', U0); reqPos('z', z); reqPos('t', t);
  if (t > T0) throw new RangeError('t must not exceed 3600 s');
  const C = 0.0573 * Math.sqrt(1 + 0.15 * U0);
  const Uz = U0 * (1 + C * Math.log(z / 10));
  const Iu = 0.06 * (1 + 0.043 * U0) * Math.pow(z / 10, -0.22);
  const gf = 1 - 0.41 * Iu * Math.log(t / T0);
  return { u: Uz * gf, Uz, Iu, C, gustFactor: gf, sigma: Iu * Uz };
}

/**
 * Inverse of `froyaWind`: the 1-hour mean wind speed at 10 m that gives a stated speed u at height z averaged
 * over t (bisection; u(z, t) increases monotonically with U0).
 * Reference and validity: see `froyaWind`.
 * @param {object} o
 * @param {number} o.u wind speed at height z averaged over t, m/s
 * @param {number} [o.z=10] height, m
 * @param {number} [o.t=3600] averaging period, s
 * @returns {number} U0, m/s
 */
export function froyaReferenceSpeed({ u, z = 10, t = T0 }) {
  reqNonNeg('u', u);
  if (u === 0) return 0;
  let lo = 0; let hi = Math.max(1, u);
  while (froyaWind({ U0: hi, z, t }).u < u) { hi *= 2; if (hi > 1e4) throw new RangeError('u cannot be matched by the Frøya profile'); }
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (lo + hi);
    if (froyaWind({ U0: m, z, t }).u < u) lo = m; else hi = m;
    if (hi - lo <= 1e-14 * hi) break;
  }
  return 0.5 * (lo + hi);
}

/**
 * Convert a wind speed between heights and averaging periods through the Frøya model.
 * Reference and validity: see `froyaWind`.
 * @param {object} o
 * @param {number} o.u known wind speed, m/s
 * @param {number} o.zFrom height of the known speed, m
 * @param {number} o.tFrom averaging period of the known speed, s
 * @param {number} o.zTo target height, m
 * @param {number} o.tTo target averaging period, s
 * @returns {{u:number,U0:number}} u converted speed, m/s; U0 implied 1-hour mean at 10 m, m/s
 */
export function convertWindSpeed({ u, zFrom, tFrom, zTo, tTo }) {
  const U0 = froyaReferenceSpeed({ u, z: zFrom, t: tFrom });
  return { u: froyaWind({ U0, z: zTo, t: tTo }).u, U0 };
}

/**
 * Gust factor relative to the 1-hour mean at the same height, G = u(z, t)/U(z) = 1 − 0.41 Iu(z) ln(t/3600).
 * Reference and validity: see `froyaWind` (ISO 19901-1 Eq. A.6).
 * @param {object} o
 * @param {number} o.U0 1-hour mean wind speed at 10 m, m/s
 * @param {number} [o.z=10] height, m
 * @param {number} [o.t=3] gust averaging period, s
 * @returns {number} gust factor, -
 */
export function gustFactor({ U0, z = 10, t = 3 }) { return froyaWind({ U0, z, t }).gustFactor; }

/**
 * Frøya (NPD) spectrum of the longitudinal wind-speed fluctuations at a point:
 *   S(f) = 320 (U0/10)² (z/10)^0.45 / (1 + f̃ⁿ)^(5/(3n)),   f̃ = 172 f (z/10)^(2/3) (U0/10)^(−0.75),   n = 0.468.
 * Reference: ISO 19901-1 (2015) §A.7.4, Eqs. A.8–A.9; NORSOK N-003 (2017) §6.3.2; DNV-RP-C205 (2019) §2.3.4.12.
 * Validity: strong winds over open sea, 0.00167 Hz ≤ f ≤ 0.5 Hz; particularly suited to low-frequency excitation
 * of compliant and moored structures.
 * @param {number} f frequency, Hz
 * @param {object} o
 * @param {number} o.U0 1-hour mean wind speed at 10 m, m/s
 * @param {number} [o.z=10] height, m
 * @returns {number} S(f), (m/s)²/Hz
 */
export function froyaSpectrum(f, { U0, z = 10 }) {
  reqPos('U0', U0); reqPos('z', z);
  if (!(f > 0)) return 0;
  const n = 0.468;
  const ft = 172 * f * Math.pow(z / 10, 2 / 3) * Math.pow(U0 / 10, -0.75);
  return (320 * (U0 / 10) ** 2 * Math.pow(z / 10, 0.45)) / Math.pow(1 + Math.pow(ft, n), 5 / (3 * n));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* IEC 61400-1                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/** Reference turbulence intensities Iref at 15 m/s per IEC 61400-1 (2019) Table 1 (A+ added in edition 4). */
export const IEC_IREF = Object.freeze({ 'A+': 0.18, A: 0.16, B: 0.14, C: 0.12 });
/** Reference wind speeds Vref (10-min mean, 50-year, hub height) per IEC 61400-1 (2019) Table 1, m/s. */
export const IEC_VREF = Object.freeze({ I: 50, II: 42.5, III: 37.5 });

/**
 * Longitudinal turbulence scale parameter Λ1 = 0.7 z for z ≤ 60 m, 42 m for z ≥ 60 m.
 * Reference: IEC 61400-1 (2019) §6.3.1, Eq. 5.
 * @param {number} z hub height, m
 * @returns {number} Λ1, m
 */
export function iecTurbulenceScale(z) { reqPos('z', z); return z <= 60 ? 0.7 * z : 42; }

/**
 * Normal turbulence model: representative (90 % quantile) standard deviation of the longitudinal wind speed,
 * σ1 = Iref (0.75 Vhub + b), b = 5.6 m/s.
 * Reference: IEC 61400-1 (2019) §6.3.2.3, Eq. 10.
 * Validity: standard wind-turbine classes; site-specific offshore turbulence is usually lower.
 * @param {number} Vhub 10-min mean wind speed at hub height, m/s
 * @param {number} Iref reference turbulence intensity (see `IEC_IREF`), -
 * @returns {{sigma1:number,TI:number}} sigma1 m/s; TI = σ1/Vhub (-)
 */
export function iecNtm(Vhub, Iref) {
  reqPos('Vhub', Vhub); reqPos('Iref', Iref);
  const sigma1 = Iref * (0.75 * Vhub + 5.6);
  return { sigma1, TI: sigma1 / Vhub };
}

/**
 * Extreme turbulence model, σ1 = c Iref [0.072 (Vave/c + 3)(Vhub/c − 4) + 10], c = 2 m/s.
 * Reference: IEC 61400-1 (2019) §6.3.3.4, Eq. 19 (Vave = 0.2 Vref for the standard classes).
 * @param {object} o
 * @param {number} o.Vhub 10-min mean wind speed at hub height, m/s
 * @param {number} o.Iref reference turbulence intensity, -
 * @param {number} o.Vave annual mean wind speed at hub height, m/s
 * @returns {{sigma1:number,TI:number}} sigma1 m/s; TI (-)
 */
export function iecEtm({ Vhub, Iref, Vave }) {
  reqPos('Vhub', Vhub); reqPos('Iref', Iref); reqPos('Vave', Vave);
  const c = 2;
  const sigma1 = c * Iref * (0.072 * (Vave / c + 3) * (Vhub / c - 4) + 10);
  return { sigma1, TI: sigma1 / Vhub };
}

/**
 * Kaimal spectrum of a turbulence velocity component,
 *   f S_k(f) / σ_k² = 4 f L_k/Vhub / (1 + 6 f L_k/Vhub)^(5/3),
 * with L_k = 8.1 Λ1, 2.7 Λ1, 0.66 Λ1 and σ_k = σ1, 0.8 σ1, 0.5 σ1 for the longitudinal (u), lateral (v) and
 * vertical (w) components. The spectrum integrates exactly to σ_k² over 0 ≤ f < ∞.
 * Reference: IEC 61400-1 (2019) Annex C, Eq. C.14 and Table C.1; Kaimal et al. (1972).
 * @param {number} f frequency, Hz
 * @param {object} o
 * @param {number} o.Vhub 10-min mean wind speed at hub height, m/s
 * @param {number} o.sigma1 longitudinal standard deviation, m/s
 * @param {number} o.z hub height (for Λ1), m
 * @param {'u'|'v'|'w'} [o.component='u']
 * @returns {number} S_k(f), (m/s)²/Hz
 */
export function kaimalSpectrum(f, { Vhub, sigma1, z, component = 'u' }) {
  reqPos('Vhub', Vhub); reqNonNeg('sigma1', sigma1);
  const lam = iecTurbulenceScale(z);
  let L; let sig;
  if (component === 'u') { L = 8.1 * lam; sig = sigma1; }
  else if (component === 'v') { L = 2.7 * lam; sig = 0.8 * sigma1; }
  else if (component === 'w') { L = 0.66 * lam; sig = 0.5 * sigma1; }
  else throw new RangeError("component must be 'u', 'v' or 'w'");
  if (!(f >= 0)) return 0;
  const r = L / Vhub;
  return (sig * sig * 4 * r) / Math.pow(1 + 6 * f * r, 5 / 3);
}

/**
 * Extreme wind speed model (EWM): 50-year and 1-year extreme wind speeds as functions of height.
 *   steady:    Ve50(z) = 1.4 Vref (z/zhub)^0.11,  Ve1(z) = 0.8 Ve50(z)          (3-s gust)
 *   turbulent: V50(z) = Vref (z/zhub)^0.11,       V1(z) = 0.8 V50(z),  σ1 = 0.11 Vhub   (10-min mean)
 * Reference: IEC 61400-1 (2019) §6.3.3.2, Eqs. 13–17.
 * @param {object} o
 * @param {number} o.Vref reference wind speed (10-min mean, 50-year, hub height), m/s
 * @param {number} o.zhub hub height, m
 * @param {number} [o.z=zhub] height, m
 * @returns {{Ve50:number,Ve1:number,V50:number,V1:number,sigma1:number}} speeds m/s; sigma1 for the turbulent
 *   model with Vhub = V50(zhub) = Vref, m/s
 */
export function iecExtremeWind({ Vref, zhub, z = zhub }) {
  reqPos('Vref', Vref); reqPos('zhub', zhub); reqPos('z', z);
  const sh = Math.pow(z / zhub, 0.11);
  return { Ve50: 1.4 * Vref * sh, Ve1: 0.8 * 1.4 * Vref * sh, V50: Vref * sh, V1: 0.8 * Vref * sh, sigma1: 0.11 * Vref };
}

/**
 * Extreme operating gust (EOG): gust magnitude and the wind-speed time history at hub height.
 *   Vgust = min{ 1.35 (Ve1 − Vhub);  3.3 σ1 / (1 + 0.1 D/Λ1) }
 *   V(t) = V − 0.37 Vgust sin(3π t/T) [1 − cos(2π t/T)]  for 0 ≤ t ≤ T,  V otherwise;   T = 10.5 s
 * Reference: IEC 61400-1 (2019) §6.3.3.3, Eqs. 17–18.
 * @param {object} o
 * @param {number} o.Vhub 10-min mean wind speed at hub height, m/s
 * @param {number} o.Ve1 1-year extreme (3-s gust) wind speed at hub height, m/s
 * @param {number} o.sigma1 NTM standard deviation at Vhub, m/s
 * @param {number} o.D rotor diameter, m
 * @param {number} o.zhub hub height, m
 * @param {number} [o.t] time from the start of the gust, s (when given, `V` is returned for that instant)
 * @param {number} [o.T=10.5] gust duration, s
 * @returns {{Vgust:number,T:number,V:number,Vmax:number,Vmin:number}} Vgust m/s; T s; V wind speed at t (Vhub
 *   when t is not given), m/s; Vmax, Vmin extremes of the gust history, m/s
 */
export function iecEog({ Vhub, Ve1, sigma1, D, zhub, t, T = 10.5 }) {
  reqPos('Vhub', Vhub); reqPos('Ve1', Ve1); reqNonNeg('sigma1', sigma1); reqPos('D', D); reqPos('T', T);
  const lam = iecTurbulenceScale(zhub);
  const Vgust = Math.max(0, Math.min(1.35 * (Ve1 - Vhub), (3.3 * sigma1) / (1 + (0.1 * D) / lam)));
  const shape = (tt) => (tt < 0 || tt > T ? 0 : -0.37 * Math.sin((3 * Math.PI * tt) / T) * (1 - Math.cos((2 * Math.PI * tt) / T)));
  let mx = 0; let mn = 0;
  for (let i = 0; i <= 2000; i++) { const s = shape((T * i) / 2000); if (s > mx) mx = s; if (s < mn) mn = s; }
  const V = t === undefined || t === null ? Vhub : Vhub + Vgust * shape(t);
  return { Vgust, T, V, Vmax: Vhub + Vgust * mx, Vmin: Vhub + Vgust * mn };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Wind force and air properties                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Shape coefficients Cs for perpendicular wind approach.
 * Reference: API RP 2A-WSD (2014) §5.3.2.4 (22nd ed.; §2.3.2e in the 21st ed.).
 */
export const API_SHAPE_COEFFICIENTS = Object.freeze({
  beams: 1.5,
  buildingSides: 1.5,
  cylindricalSections: 0.5,
  overallPlatformProjectedArea: 1.0,
});

/**
 * Density of air from the ideal-gas law, ρ = (p − 0.378 e)/(Rd T), with Rd = 287.05 J/(kg·K) and vapour pressure
 * e = RH · es(T), es from the Tetens formula es = 610.78 exp[17.27 T_C/(T_C + 237.3)] Pa.
 * Reference: ideal-gas law (ISO 2533 standard atmosphere constants); Tetens (1930) / Monteith & Unsworth (2008).
 * Validity: −40 °C to +50 °C, pressures near atmospheric.
 * @param {object} [o]
 * @param {number} [o.temperatureC=15] air temperature, °C
 * @param {number} [o.pressure=101325] absolute pressure, Pa
 * @param {number} [o.relativeHumidity=0] relative humidity, 0–1
 * @returns {number} air density, kg/m³
 */
export function airDensity({ temperatureC = 15, pressure = 101325, relativeHumidity = 0 } = {}) {
  reqPos('pressure', pressure);
  if (!(temperatureC > -273.15)) throw new RangeError('temperatureC must exceed absolute zero');
  if (!(relativeHumidity >= 0 && relativeHumidity <= 1)) throw new RangeError('relativeHumidity must lie in [0, 1]');
  const e = relativeHumidity * 610.78 * Math.exp((17.27 * temperatureC) / (temperatureC + 237.3));
  return (pressure - 0.378 * e) / (287.05 * (temperatureC + 273.15));
}

/**
 * Wind force on a body, F = ½ ρ Cs Ch A U² (optionally reduced by the angle of incidence for a flat face,
 * F = ½ ρ Cs Ch A U² sin α, α = angle between the wind and the surface).
 * Reference: API RP 2A-WSD (2014) §5.3.2.2, Eq. 5.6 (F = (ρ/2) u² Cs A); DNV-RP-C205 (2019) §5.3.1
 * (FW = C q S sin α). The height coefficient Ch is 1 when U is taken at the elevation of the area; use
 * `heightCoefficient` to refer a 10-m speed to that elevation.
 * Validity: steady (quasi-static) wind load; no shielding, solidification or dynamic amplification.
 * @param {object} o
 * @param {number} o.U wind speed (appropriate gust duration), m/s
 * @param {number} o.A projected area normal to the wind, m²
 * @param {number} [o.Cs=1] shape coefficient (see `API_SHAPE_COEFFICIENTS`), -
 * @param {number} [o.Ch=1] height coefficient, -
 * @param {number} [o.rho=1.225] air density, kg/m³
 * @param {number} [o.alpha=π/2] angle between wind direction and the exposed surface, rad
 * @returns {{F:number,q:number}} F force N; q = ½ ρ U² basic pressure Pa
 */
export function windForce({ U, A, Cs = 1, Ch = 1, rho = RHO_AIR, alpha = Math.PI / 2 }) {
  reqNonNeg('U', U); reqNonNeg('A', A); reqNonNeg('Cs', Cs); reqNonNeg('Ch', Ch); reqPos('rho', rho);
  const q = 0.5 * rho * U * U;
  return { F: q * Cs * Ch * A * Math.sin(alpha), q };
}

/**
 * Height coefficient referring a wind pressure at the reference height to height z through a power-law profile,
 * Ch = [U(z)/U(zRef)]² = (z/zRef)^(2α).
 * Reference: derived from the profile (API RP 2A-WSD §5.3.2.1 applies the profile directly; tabulated height
 * coefficients of mobile-unit codes are step approximations of the same quantity and are not reproduced here).
 * @param {number} z height of the area centroid, m
 * @param {number} [zRef=10] reference height, m
 * @param {number} [alpha=0.12] shear exponent, -
 * @returns {number} Ch, -
 */
export function heightCoefficient(z, zRef = 10, alpha = 0.12) {
  reqPos('z', z); reqPos('zRef', zRef);
  return Math.pow(z / zRef, 2 * alpha);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Weibull wind climate                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Two-parameter Weibull wind-speed distribution: density, exceedance and moments.
 *   f(U) = (k/A)(U/A)^(k−1) exp[−(U/A)^k],  F(U) = 1 − exp[−(U/A)^k],  mean = A Γ(1 + 1/k),
 *   mean cube = A³ Γ(1 + 3/k),  mean power density = ½ ρ A³ Γ(1 + 3/k).
 * Reference: IEC 61400-1 (2019) §6.3.2.1 (Rayleigh is k = 2); DNV-RP-C205 (2019) §2.3.1.
 * @param {object} o
 * @param {number} o.A scale parameter, m/s
 * @param {number} o.k shape parameter, -
 * @param {number} [o.U] wind speed at which pdf/cdf are evaluated, m/s
 * @param {number} [o.rho=1.225] air density, kg/m³
 * @returns {{pdf:number,cdf:number,mean:number,std:number,meanCube:number,powerDensity:number,mode:number,
 *   median:number}} pdf s/m; cdf (-); mean, std, mode, median m/s; meanCube m³/s³; powerDensity W/m²
 */
export function weibullWind({ A, k, U, rho = RHO_AIR }) {
  reqPos('A', A); reqPos('k', k);
  const g1 = gammaFn(1 + 1 / k); const g2 = gammaFn(1 + 2 / k); const g3 = gammaFn(1 + 3 / k);
  let pdf = NaN; let cdf = NaN;
  if (U !== undefined && U !== null) {
    if (U < 0) { pdf = 0; cdf = 0; } else {
      const x = Math.pow(U / A, k);
      pdf = U === 0 ? (k === 1 ? 1 / A : k < 1 ? Infinity : 0) : (k / A) * Math.pow(U / A, k - 1) * Math.exp(-x);
      cdf = 1 - Math.exp(-x);
    }
  }
  return {
    pdf, cdf, mean: A * g1, std: A * Math.sqrt(Math.max(0, g2 - g1 * g1)), meanCube: A ** 3 * g3,
    powerDensity: 0.5 * rho * A ** 3 * g3, mode: k > 1 ? A * Math.pow((k - 1) / k, 1 / k) : 0,
    median: A * Math.pow(Math.LN2, 1 / k),
  };
}

/**
 * Weibull parameters from the mean and standard deviation (method of moments): k solves
 * Γ(1 + 2/k)/Γ(1 + 1/k)² − 1 = (σ/μ)², then A = μ/Γ(1 + 1/k).
 * Reference: method of moments; Justus et al. (1978) J. Appl. Meteorol. 17.
 * Validity: coefficient of variation between about 0.03 and 4 (k from 0.3 to 50).
 * @param {number} mean mean wind speed, m/s
 * @param {number} std standard deviation, m/s
 * @returns {{A:number,k:number}} A m/s; k (-)
 */
export function weibullFromMeanStd(mean, std) {
  reqPos('mean', mean); reqPos('std', std);
  const target = (std / mean) ** 2;
  const f = (k) => gammaFn(1 + 2 / k) / gammaFn(1 + 1 / k) ** 2 - 1 - target;
  let lo = 0.3; let hi = 50;
  if (f(lo) < 0 || f(hi) > 0) throw new RangeError('std/mean is outside the range representable by a Weibull shape of 0.3–50');
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (lo + hi);
    if (f(m) > 0) lo = m; else hi = m;
    if (hi - lo < 1e-13 * hi) break;
  }
  const k = 0.5 * (lo + hi);
  return { A: mean / gammaFn(1 + 1 / k), k };
}

/**
 * Weibull parameters by maximum likelihood: k solves Σ x^k ln x / Σ x^k − 1/k − mean(ln x) = 0 and
 * A = (mean(x^k))^(1/k). Zero and negative values (calms) are excluded and reported.
 * Reference: Cohen (1965) Technometrics 7(4); Seguro & Lambert (2000) J. Wind Eng. Ind. Aerodyn. 85.
 * Validity: at least two distinct positive values.
 * @param {ArrayLike<number>} speeds wind-speed sample, m/s
 * @returns {{A:number,k:number,n:number,nExcluded:number,logLik:number}} A m/s; k (-); n values used; nExcluded
 *   non-positive values dropped; logLik log-likelihood at the optimum
 */
export function weibullFitMle(speeds) {
  let n = 0; let sumLn = 0; let mx = 0;
  for (let i = 0; i < speeds.length; i++) {
    const x = speeds[i];
    if (x > 0 && Number.isFinite(x)) { n++; sumLn += Math.log(x); if (x > mx) mx = x; }
  }
  if (n < 2) throw new RangeError('speeds must contain at least two positive values');
  const meanLn = sumLn / n;
  // scale by the maximum so that x^k cannot overflow
  const h = (k) => {
    let a = 0; let b = 0;
    for (let i = 0; i < speeds.length; i++) {
      const x = speeds[i];
      if (x > 0 && Number.isFinite(x)) { const p = Math.pow(x / mx, k); a += p * Math.log(x); b += p; }
    }
    return a / b - 1 / k - meanLn;
  };
  let lo = 0.02; let hi = 2;
  while (h(hi) < 0) { hi *= 2; if (hi > 1e4) throw new RangeError('speeds have no spread; Weibull shape is unbounded'); }
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (lo + hi);
    if (h(m) < 0) lo = m; else hi = m;
    if (hi - lo < 1e-13 * hi) break;
  }
  const k = 0.5 * (lo + hi);
  let s = 0;
  for (let i = 0; i < speeds.length; i++) { const x = speeds[i]; if (x > 0 && Number.isFinite(x)) s += Math.pow(x / mx, k); }
  const A = mx * Math.pow(s / n, 1 / k);
  const logLik = n * Math.log(k) - n * k * Math.log(A) + (k - 1) * sumLn - n;
  return { A, k, n, nExcluded: speeds.length - n, logLik };
}
