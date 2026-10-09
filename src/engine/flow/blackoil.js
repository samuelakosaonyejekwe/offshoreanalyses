// Black-oil PVT correlations. Public API is SI: P absolute in Pa, T in K, gas–oil ratios in Sm³/Sm³, rates in
// Sm³/s at stock-tank conditions, densities kg/m³, viscosities Pa·s, interfacial tension N/m.
// The published correlations are evaluated internally in their original oilfield units (psia, °F, scf/STB, cP)
// with the published coefficients; conversions use exact factors.
// Standard conditions for all volumetric ratios: 60 °F (288.7056 K) and 101 325 Pa.

import { brineDensity, brineViscosity, seawaterViscosity } from './water.js';

/** Universal gas constant, J/(mol·K). */
export const R_GAS = 8.314462618;
/** Standard temperature (60 °F), K. */
export const T_SC = 288.70555555555555;
/** Standard pressure, Pa. */
export const P_SC = 101325;
/** Molar mass of dry air, kg/mol. */
export const M_AIR = 0.0289647;

const PSI = 6894.757293168;
const SCF_PER_STB = 5.614583333333333; // (Sm³/Sm³) → scf/STB, 1 bbl = 5.614583 ft³
const RHO_WATER_SC = 999.016; // kg/m³ at 60 °F (62.366 lb/ft³), basis of oil specific gravity

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}
const toF = (T) => T * 1.8 - 459.67;
const toR = (T) => T * 1.8;

/**
 * Oil specific gravity (water = 1 at 60 °F) from API gravity, γo = 141.5/(131.5 + API).
 * Reference: API gravity definition. Validity: API > −131.5 (practically 5–70).
 * @param {number} api API gravity, °API
 * @returns {number} specific gravity, –
 */
export function apiToSg(api) {
  need('api', api, 0, 120);
  return 141.5 / (131.5 + api);
}

/**
 * API gravity from specific gravity, API = 141.5/γo − 131.5. Reference: API gravity definition.
 * @param {number} sg oil specific gravity, – (0.55–1.08)
 * @returns {number} °API
 */
export function sgToApi(sg) {
  need('sg', sg, 0.5, 1.1);
  return 141.5 / sg - 131.5;
}

/**
 * Gas density at standard conditions from gas gravity (ideal gas, Z = 1): ρ = γg M_air P_sc/(R T_sc).
 * Reference: ideal-gas law. @param {number} gammaG gas specific gravity (air = 1), –
 * @returns {number} kg/Sm³
 */
export function gasDensitySc(gammaG) {
  need('gammaG', gammaG, 0.05, 3);
  return gammaG * M_AIR * P_SC / (R_GAS * T_SC);
}

/**
 * Pseudo-critical properties of a natural gas from its gravity with sour-gas correction.
 * Sutton (1985): Ppc = 756.8 − 131.0 γ − 3.6 γ² psia, Tpc = 169.2 + 349.5 γ − 74.0 γ² °R (applied to the total gas
 * gravity). Wichert & Aziz (1972): ε = 120 (A^0.9 − A^1.6) + 15 (B^0.5 − B⁴) °R with A = yH2S + yCO2, B = yH2S;
 * T'pc = Tpc − ε, P'pc = Ppc T'pc/(Tpc + B (1 − B) ε).
 * References: Sutton, SPE 14265 (1985); Wichert & Aziz, Hydrocarbon Processing 51(5) (1972).
 * Validity: 0.57 < γ < 1.68; yCO2 ≤ 0.544, yH2S ≤ 0.738.
 * @param {object} o
 * @param {number} o.gammaG gas specific gravity (air = 1), –
 * @param {number} [o.yCO2=0] CO2 mole fraction, –
 * @param {number} [o.yH2S=0] H2S mole fraction, –
 * @returns {{ Tpc: number, Ppc: number, TpcUncorrected: number, PpcUncorrected: number, epsilon: number }}
 *   corrected pseudo-critical temperature K and pressure Pa, Sutton values before correction, ε in K
 */
export function gasPseudoCriticals({ gammaG, yCO2 = 0, yH2S = 0 }) {
  need('gammaG', gammaG, 0.5, 2.0);
  need('yCO2', yCO2, 0, 0.8);
  need('yH2S', yH2S, 0, 0.8);
  const Ppc = 756.8 - 131.0 * gammaG - 3.6 * gammaG * gammaG;
  const Tpc = 169.2 + 349.5 * gammaG - 74.0 * gammaG * gammaG;
  const A = yCO2 + yH2S;
  const B = yH2S;
  const eps = A > 0 ? 120 * (A ** 0.9 - A ** 1.6) + 15 * (Math.sqrt(B) - B ** 4) : 0;
  const Tc = Tpc - eps;
  const Pc = Ppc * Tc / (Tpc + B * (1 - B) * eps);
  return { Tpc: Tc / 1.8, Ppc: Pc * PSI, TpcUncorrected: Tpc / 1.8, PpcUncorrected: Ppc * PSI, epsilon: eps / 1.8 };
}

const DAK = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210];

function dakZofRho(rr, Tr) {
  const A = DAK;
  const r2 = rr * rr;
  const c1 = A[0] + A[1] / Tr + A[2] / Tr ** 3 + A[3] / Tr ** 4 + A[4] / Tr ** 5;
  const c2 = A[5] + A[6] / Tr + A[7] / (Tr * Tr);
  const c3 = A[8] * (A[6] / Tr + A[7] / (Tr * Tr));
  return 1 + c1 * rr + c2 * r2 - c3 * rr ** 5 + A[9] * (1 + A[10] * r2) * (r2 / Tr ** 3) * Math.exp(-A[10] * r2);
}

function zDak(Tr, Pr) {
  // Unknown: reduced density ρr = 0.27 Ppr/(Z Tpr). g(ρr) = ρr Z(ρr) − 0.27 Ppr/Tpr = 0, smallest positive root.
  const target = 0.27 * Pr / Tr;
  const g = (r) => r * dakZofRho(r, Tr) - target;
  let r = target; // ideal-gas start
  let ok = false;
  for (let it = 0; it < 60; it++) {
    const h = 1e-6 * Math.max(r, 1e-3);
    const gr = g(r);
    const dg = (g(r + h) - g(r - h)) / (2 * h);
    if (!(dg > 0)) break;
    let rn = r - gr / dg;
    if (rn <= 0) rn = 0.5 * r;
    if (Math.abs(rn - r) < 1e-13 * Math.max(1, r)) { r = rn; ok = true; break; }
    r = rn;
  }
  if (!ok || !(r > 0) || r > 4) {
    // bracket the lowest root on a grid, then bisect
    let a = 0;
    let b = 0;
    const n = 400;
    for (let i = 1; i <= n; i++) {
      b = 4 * i / n;
      if (g(b) > 0) break;
      a = b;
    }
    for (let it = 0; it < 80; it++) {
      const m = 0.5 * (a + b);
      if (g(m) > 0) b = m; else a = m;
    }
    r = 0.5 * (a + b);
  }
  return target / r;
}

function zHallYarborough(Tr, Pr) {
  const t = 1 / Tr;
  const a = 0.06125 * Pr * t * Math.exp(-1.2 * (1 - t) ** 2);
  const b = 14.76 * t - 9.76 * t * t + 4.58 * t ** 3;
  const c = 90.7 * t - 242.2 * t * t + 42.4 * t ** 3;
  const d = 2.18 + 2.82 * t;
  const F = (y) => -a + (y + y * y + y ** 3 - y ** 4) / (1 - y) ** 3 - b * y * y + c * y ** d;
  const dF = (y) => (1 + 4 * y + 4 * y * y - 4 * y ** 3 + y ** 4) / (1 - y) ** 4 - 2 * b * y + c * d * y ** (d - 1);
  let y = Math.min(0.5, Math.max(1e-6, a)); // ideal-gas start: Z = 1 → y = a
  let ok = false;
  for (let it = 0; it < 100; it++) {
    let yn = y - F(y) / dF(y);
    if (!(yn > 0)) yn = 0.5 * y;
    if (yn >= 1) yn = 0.5 * (y + 1);
    if (Math.abs(yn - y) < 1e-14) { y = yn; ok = true; break; }
    y = yn;
  }
  if (!ok) {
    let lo = 1e-12;
    let hi = 0.999;
    for (let it = 0; it < 100; it++) {
      const m = 0.5 * (lo + hi);
      if (F(m) > 0) hi = m; else lo = m;
    }
    y = 0.5 * (lo + hi);
  }
  return a / y;
}

/**
 * Gas compressibility factor from reduced conditions.
 * `method = 'dak'`: Dranchuk & Abou-Kassem (1975) eleven-constant fit of the Standing–Katz chart
 * (A1…A11 = 0.3265, −1.0700, −0.5339, 0.01569, −0.05165, 0.5475, −0.7361, 0.1844, 0.1056, 0.6134, 0.7210).
 * `method = 'hallYarborough'`: Hall & Yarborough (1973) hard-sphere form.
 * References: Dranchuk & Abou-Kassem, J. Can. Pet. Technol. 14(3) (1975); Hall & Yarborough, Oil & Gas J. (18 June
 * 1973). Validity: 1.0 ≤ Tpr ≤ 3.0, 0.2 ≤ Ppr ≤ 30 (DAK); 1.2 ≤ Tpr ≤ 3, Ppr ≤ 24 (H–Y). Evaluation is permitted
 * down to Tpr = 0.92 for marching solvers, with degraded accuracy.
 * @param {object} o
 * @param {number} o.Tpr pseudo-reduced temperature, –
 * @param {number} o.Ppr pseudo-reduced pressure, –
 * @param {'dak'|'hallYarborough'} [o.method='dak']
 * @returns {number} Z, –
 */
export function zFactorReduced({ Tpr, Ppr, method = 'dak' }) {
  need('Tpr', Tpr, 0.92, 5);
  need('Ppr', Ppr, 0, 40);
  if (Ppr === 0) return 1;
  if (method === 'dak') return zDak(Tpr, Ppr);
  if (method === 'hallYarborough') return zHallYarborough(Tpr, Ppr);
  throw new RangeError(`method = ${method} is not one of dak, hallYarborough`);
}

/**
 * Gas compressibility factor from gas gravity (Sutton pseudo-criticals with Wichert–Aziz correction, then DAK or
 * Hall–Yarborough). See `gasPseudoCriticals` and `zFactorReduced` for references and validity.
 * Temperatures below Tpr = 0.92 are evaluated at Tpr = 0.92 (the correlation has no meaning there).
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.gammaG gas gravity (air = 1), –
 * @param {number} [o.yCO2=0] mole fraction CO2
 * @param {number} [o.yH2S=0] mole fraction H2S
 * @param {'dak'|'hallYarborough'} [o.method='dak']
 * @returns {number} Z, –
 */
export function zFactor({ P, T, gammaG, yCO2 = 0, yH2S = 0, method = 'dak' }) {
  need('P', P, 0, 3e8);
  need('T', T, 100, 1000);
  const pc = gasPseudoCriticals({ gammaG, yCO2, yH2S });
  return zFactorReduced({ Tpr: Math.max(0.92, Math.min(5, T / pc.Tpc)), Ppr: Math.min(40, P / pc.Ppc), method });
}

/**
 * Real-gas density ρ = P M/(Z R T) with M = γg M_air.
 * Reference: real-gas law. @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.gammaG gas gravity, –
 * @param {number} [o.Z] compressibility factor; computed by `zFactor` when omitted
 * @param {number} [o.yCO2=0] @param {number} [o.yH2S=0]
 * @returns {number} density, kg/m³
 */
export function gasDensity({ P, T, gammaG, Z, yCO2 = 0, yH2S = 0 }) {
  const z = Z ?? zFactor({ P, T, gammaG, yCO2, yH2S });
  return P * gammaG * M_AIR / (z * R_GAS * T);
}

/**
 * Gas formation volume factor Bg = (P_sc/T_sc) Z T/P (Z_sc = 1).
 * Reference: real-gas law. @param {object} o
 * @param {number} o.P Pa @param {number} o.T K @param {number} o.Z compressibility factor, –
 * @returns {number} m³/Sm³
 */
export function gasFvf({ P, T, Z }) {
  need('P', P, 1, 3e8);
  need('T', T, 100, 1000);
  need('Z', Z, 0.05, 5);
  return P_SC * Z * T / (T_SC * P);
}

/**
 * Gas viscosity by Lee, Gonzalez & Eakin (1966):
 * μ = 1e-4 K exp(X ρ^Y) cP, K = (9.4 + 0.02 M) T^1.5/(209 + 19 M + T), X = 3.5 + 986/T + 0.01 M,
 * Y = 2.4 − 0.2 X, with T in °R, ρ in g/cm³, M in g/mol.
 * Reference: Lee, Gonzalez & Eakin, "The viscosity of natural gases", JPT 18(8) (1966) 997–1000.
 * Validity: 311–444 K, 0.7–55 MPa, hydrocarbon gases with ≤ 3.2 mol % CO2; standard deviation 2.7 %.
 * @param {object} o
 * @param {number} o.T temperature, K
 * @param {number} o.rho gas density, kg/m³
 * @param {number} o.M molar mass, kg/mol
 * @returns {number} viscosity, Pa·s
 */
export function gasViscosityLee({ T, rho, M }) {
  need('T', T, 100, 1000);
  need('rho', rho, 0, 800);
  need('M', M, 0.002, 0.2);
  const Mg = M * 1000;
  const TR = toR(T);
  const K = (9.4 + 0.02 * Mg) * TR ** 1.5 / (209 + 19 * Mg + TR);
  const X = 3.5 + 986 / TR + 0.01 * Mg;
  const Y = 2.4 - 0.2 * X;
  return 1e-4 * K * Math.exp(X * (rho / 1000) ** Y) * 1e-3;
}

function vbGasGravity(gammaG, api, Tsep, Psep) {
  // Vazquez–Beggs normalisation of gas gravity to a 100 psig separator.
  if (Tsep == null || Psep == null) return gammaG;
  return gammaG * (1 + 5.912e-5 * api * toF(Tsep) * Math.log10(Psep / PSI / 114.7));
}

function glasoRsPbStar(api, gammaG, TF) {
  return (rs) => (rs / gammaG) ** 0.816 * TF ** 0.172 / api ** 0.989;
}

/**
 * Solution gas–oil ratio at pressure P (saturated oil).
 * `standing`: Rs = γg [(p/18.2 + 1.4)·10^(0.0125 API − 0.00091 T)]^1.2048 (p psia, T °F, Rs scf/STB); the exponent
 *   is evaluated as 1/0.83 so that it inverts the bubble-point equation exactly.
 * `vazquezBeggs`: Rs = C1 γgs p^C2 exp(C3 API/(T + 460)); API ≤ 30: C1 = 0.0362, C2 = 1.0937, C3 = 25.7240;
 *   API > 30: C1 = 0.0178, C2 = 1.1870, C3 = 23.931; γgs = γg [1 + 5.912e-5 API Tsep log10(psep/114.7)].
 * `glaso`: inverse of log pb = 1.7669 + 1.7447 log pb* − 0.30218 (log pb*)²,
 *   pb* = (Rs/γg)^0.816 T^0.172/API^0.989.
 * References: Standing (1947) API Drilling & Production Practice; Vazquez & Beggs, JPT 32(6) (1980);
 * Glasø, JPT 32(5) (1980). Validity (approx.): Standing 130–7000 psia, 100–258 °F, 16.5–63.8 °API,
 * γg 0.59–0.95; Vazquez–Beggs 50–5250 psia, 70–295 °F, 16–58 °API; Glasø 165–7142 psia, 80–280 °F, 22–48 °API.
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.api oil gravity, °API
 * @param {number} o.gammaG gas gravity (air = 1), –
 * @param {'standing'|'vazquezBeggs'|'glaso'} [o.method='standing']
 * @param {number} [o.Tsep] separator temperature, K (Vazquez–Beggs gravity correction; omit for none)
 * @param {number} [o.Psep] separator absolute pressure, Pa
 * @returns {number} Rs, Sm³/Sm³
 */
export function solutionGor({ P, T, api, gammaG, method = 'standing', Tsep, Psep }) {
  need('P', P, 0, 2e8);
  need('T', T, 200, 500);
  need('api', api, 5, 80);
  need('gammaG', gammaG, 0.5, 2);
  const p = P / PSI;
  const TF = toF(T);
  let rs;
  if (method === 'standing') {
    // exponent 1/0.83 (published rounded as 1.2048) keeps Rs the exact inverse of the bubble-point equation
    rs = gammaG * ((p / 18.2 + 1.4) * 10 ** (0.0125 * api - 0.00091 * TF)) ** (1 / 0.83);
  } else if (method === 'vazquezBeggs') {
    const [C1, C2, C3] = api <= 30 ? [0.0362, 1.0937, 25.7240] : [0.0178, 1.1870, 23.931];
    rs = C1 * vbGasGravity(gammaG, api, Tsep, Psep) * p ** C2 * Math.exp(C3 * api / (TF + 460));
  } else if (method === 'glaso') {
    const lp = Math.log10(Math.max(p, 1e-9));
    const disc = 1.7447 * 1.7447 - 4 * 0.30218 * (lp - 1.7669);
    if (disc < 0) throw new RangeError(`P = ${P} Pa is above the range of the Glasø bubble-point correlation`);
    const x = (1.7447 - Math.sqrt(disc)) / (2 * 0.30218);
    const pbStar = 10 ** x;
    rs = gammaG * (pbStar * api ** 0.989 / Math.max(TF, 1) ** 0.172) ** (1 / 0.816);
  } else {
    throw new RangeError(`method = ${method} is not one of standing, vazquezBeggs, glaso`);
  }
  return Math.max(0, rs) / SCF_PER_STB;
}

/**
 * Bubble-point pressure for a given solution GOR (inverse of `solutionGor`, same correlations and references).
 * Standing: pb = 18.2 [(Rs/γg)^0.83 · 10^(0.00091 T − 0.0125 API) − 1.4].
 * @param {object} o
 * @param {number} o.Rs solution GOR at the bubble point, Sm³/Sm³
 * @param {number} o.T temperature, K
 * @param {number} o.api °API
 * @param {number} o.gammaG gas gravity, –
 * @param {'standing'|'vazquezBeggs'|'glaso'} [o.method='standing']
 * @param {number} [o.Tsep] K @param {number} [o.Psep] Pa
 * @returns {number} bubble-point pressure, Pa (absolute)
 */
export function bubblePoint({ Rs, T, api, gammaG, method = 'standing', Tsep, Psep }) {
  need('Rs', Rs, 0, 1e4);
  need('T', T, 200, 500);
  need('api', api, 5, 80);
  need('gammaG', gammaG, 0.5, 2);
  const rs = Rs * SCF_PER_STB;
  const TF = toF(T);
  let pb;
  if (method === 'standing') {
    pb = 18.2 * ((rs / gammaG) ** 0.83 * 10 ** (0.00091 * TF - 0.0125 * api) - 1.4);
  } else if (method === 'vazquezBeggs') {
    const [C1, C2, C3] = api <= 30 ? [0.0362, 1.0937, 25.7240] : [0.0178, 1.1870, 23.931];
    pb = (rs / (C1 * vbGasGravity(gammaG, api, Tsep, Psep) * Math.exp(C3 * api / (TF + 460)))) ** (1 / C2);
  } else if (method === 'glaso') {
    if (rs <= 0) return P_SC;
    const x = Math.log10(glasoRsPbStar(api, gammaG, Math.max(TF, 1))(rs));
    pb = 10 ** (1.7669 + 1.7447 * x - 0.30218 * x * x);
  } else {
    throw new RangeError(`method = ${method} is not one of standing, vazquezBeggs, glaso`);
  }
  return Math.max(pb * PSI, P_SC);
}

/**
 * Oil formation volume factor at saturation for a given dissolved GOR.
 * `standing`: Bo = 0.9759 + 0.000120 [Rs (γg/γo)^0.5 + 1.25 T]^1.2.
 * `vazquezBeggs`: Bo = 1 + C1 Rs + (T − 60)(API/γgs)(C2 + C3 Rs); API ≤ 30: 4.677e-4, 1.751e-5, −1.811e-8;
 *   API > 30: 4.670e-4, 1.100e-5, 1.337e-9.
 * `glaso`: log(Bo − 1) = −6.58511 + 2.91329 log B* − 0.27683 (log B*)², B* = Rs (γg/γo)^0.526 + 0.968 T.
 * (Rs scf/STB, T °F.) References and validity as `solutionGor`.
 * @param {object} o
 * @param {number} o.Rs dissolved GOR, Sm³/Sm³
 * @param {number} o.T temperature, K
 * @param {number} o.api °API
 * @param {number} o.gammaG gas gravity, –
 * @param {'standing'|'vazquezBeggs'|'glaso'} [o.method='standing']
 * @param {number} [o.Tsep] K @param {number} [o.Psep] Pa
 * @returns {number} Bo, m³/Sm³
 */
export function oilFvfSaturated({ Rs, T, api, gammaG, method = 'standing', Tsep, Psep }) {
  need('Rs', Rs, 0, 1e4);
  need('T', T, 200, 500);
  need('api', api, 5, 80);
  need('gammaG', gammaG, 0.5, 2);
  const rs = Rs * SCF_PER_STB;
  const TF = toF(T);
  const go = apiToSg(api);
  if (method === 'standing') {
    return 0.9759 + 0.000120 * Math.max(0, rs * Math.sqrt(gammaG / go) + 1.25 * TF) ** 1.2;
  }
  if (method === 'vazquezBeggs') {
    const [C1, C2, C3] = api <= 30 ? [4.677e-4, 1.751e-5, -1.811e-8] : [4.670e-4, 1.100e-5, 1.337e-9];
    return 1 + C1 * rs + (TF - 60) * (api / vbGasGravity(gammaG, api, Tsep, Psep)) * (C2 + C3 * rs);
  }
  if (method === 'glaso') {
    const bs = rs * (gammaG / go) ** 0.526 + 0.968 * TF;
    if (bs <= 0) return 1;
    const x = Math.log10(bs);
    return 1 + 10 ** (-6.58511 + 2.91329 * x - 0.27683 * x * x);
  }
  throw new RangeError(`method = ${method} is not one of standing, vazquezBeggs, glaso`);
}

/**
 * Isothermal compressibility of undersaturated oil (Vazquez & Beggs 1980):
 * co = (−1433 + 5 Rsb + 17.2 T − 1180 γgs + 12.61 API)/(1e5 p)  (1/psi; Rsb scf/STB, T °F, p psia).
 * Reference: Vazquez & Beggs, JPT 32(6) (1980). Validity: P ≥ Pb, ranges as `solutionGor`.
 * The result is floored at 1e-11 1/Pa where the numerator turns negative (heavy, cold, gas-free oil).
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.Rsb solution GOR at the bubble point, Sm³/Sm³
 * @param {number} o.api °API
 * @param {number} o.gammaG gas gravity (separator-normalised if available), –
 * @returns {number} compressibility, 1/Pa
 */
export function oilCompressibility({ P, T, Rsb, api, gammaG }) {
  need('P', P, 1, 3e8);
  const A = vbCoNumerator({ T, Rsb, api, gammaG });
  return Math.max(A / (1e5 * (P / PSI)) / PSI, 1e-11);
}

function vbCoNumerator({ T, Rsb, api, gammaG }) {
  return -1433 + 5 * Rsb * SCF_PER_STB + 17.2 * toF(T) - 1180 * gammaG + 12.61 * api;
}

/**
 * Dead-oil viscosity by Beggs & Robinson (1975): μod = 10^X − 1 cP, X = Y T^−1.163, Y = 10^Z,
 * Z = 3.0324 − 0.02023 API (T °F).
 * Reference: Beggs & Robinson, JPT 27(9) (1975) 1140–1141. Validity: 294–419 K (70–295 °F), 16–58 °API.
 * The correlation rises very steeply below 294 K; treat seabed-temperature values as indicative only.
 * @param {object} o
 * @param {number} o.T temperature, K (> 255.4 K so that T °F > 0)
 * @param {number} o.api °API
 * @returns {number} viscosity, Pa·s
 */
export function deadOilViscosity({ T, api }) {
  need('T', T, 256, 500);
  need('api', api, 5, 80);
  const X = 10 ** (3.0324 - 0.02023 * api) * toF(T) ** -1.163;
  return (10 ** X - 1) * 1e-3;
}

/**
 * Saturated (live) oil viscosity by Beggs & Robinson (1975): μo = A μod^B,
 * A = 10.715 (Rs + 100)^−0.515, B = 5.44 (Rs + 150)^−0.338 (Rs scf/STB, μ cP).
 * Reference: Beggs & Robinson (1975). Validity: Rs 20–2070 scf/STB, 0–5250 psig, 70–295 °F.
 * @param {object} o
 * @param {number} o.muDead dead-oil viscosity, Pa·s
 * @param {number} o.Rs dissolved GOR, Sm³/Sm³
 * @returns {number} viscosity, Pa·s
 */
export function liveOilViscosity({ muDead, Rs }) {
  need('muDead', muDead, 1e-6, 1e4);
  need('Rs', Rs, 0, 1e4);
  const rs = Rs * SCF_PER_STB;
  const A = 10.715 * (rs + 100) ** -0.515;
  const B = 5.44 * (rs + 150) ** -0.338;
  return A * (muDead * 1e3) ** B * 1e-3;
}

/**
 * Undersaturated oil viscosity by Vazquez & Beggs (1980): μo = μob (p/pb)^m,
 * m = 2.6 p^1.187 exp(−11.513 − 8.98e-5 p) (p psia).
 * Reference: Vazquez & Beggs, JPT 32(6) (1980). Validity: 141–9515 psia, 0.117–148 cP.
 * @param {object} o
 * @param {number} o.muOb viscosity at the bubble point, Pa·s
 * @param {number} o.P absolute pressure, Pa (≥ Pb)
 * @param {number} o.Pb bubble-point pressure, Pa
 * @returns {number} viscosity, Pa·s
 */
export function undersaturatedOilViscosity({ muOb, P, Pb }) {
  need('muOb', muOb, 1e-6, 1e4);
  need('Pb', Pb, 1, 3e8);
  need('P', P, Pb, 3e8);
  const p = P / PSI;
  const m = 2.6 * p ** 1.187 * Math.exp(-11.513 - 8.98e-5 * p);
  return muOb * (P / Pb) ** m;
}

/**
 * Gas–oil interfacial tension by the Baker & Swerdloff (1956) charts in the form fitted by Beggs:
 * dead oil σ68 = 39 − 0.2571 API, σ100 = 37.5 − 0.2571 API dyn/cm (linear interpolation between 68 and 100 °F,
 * no extrapolation outside), multiplied by the dissolved-gas factor C = 1 − 0.024 p^0.45 (p psia).
 * Reference: Baker & Swerdloff, Oil & Gas J. (2 Jan 1956); Beggs in Bradley (ed.), Petroleum Engineering Handbook
 * (SPE 1987) ch. 22. Validity: the factor reaches zero near 27 MPa; the result is floored at 1 mN/m.
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.api °API
 * @returns {number} interfacial tension, N/m
 */
export function gasOilIft({ P, T, api }) {
  need('P', P, 0, 3e8);
  need('T', T, 200, 500);
  need('api', api, 5, 80);
  const TF = Math.min(100, Math.max(68, toF(T)));
  const s68 = 39 - 0.2571 * api;
  const s100 = 37.5 - 0.2571 * api;
  const sd = s68 + (TF - 68) * (s100 - s68) / 32;
  const C = 1 - 0.024 * (P / PSI) ** 0.45;
  return Math.max(sd * C, 1) * 1e-3;
}

/**
 * Gas–water interfacial tension: σ(74 °F) = 75 − 1.108 p^0.349, σ(280 °F) = 53 − 0.1048 p^0.637 dyn/cm (p psia),
 * linear interpolation in temperature, clamped to the two curves outside 74–280 °F.
 * Reference: Hough, Rzasa & Wood (1951) data as fitted by Beggs (Petroleum Engineering Handbook, SPE 1987, ch. 22).
 * Validity: up to about 55 MPa; floored at 1 mN/m.
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @returns {number} interfacial tension, N/m
 */
export function gasWaterIft({ P, T }) {
  need('P', P, 0, 3e8);
  need('T', T, 200, 500);
  const p = P / PSI;
  const TF = Math.min(280, Math.max(74, toF(T)));
  const s74 = 75 - 1.108 * p ** 0.349;
  const s280 = 53 - 0.1048 * p ** 0.637;
  return Math.max(s74 + (TF - 74) * (s280 - s74) / 206, 1) * 1e-3;
}

/**
 * Oil PVT state at (P, T) for a stock-tank oil with total produced GOR `Rsb`: bubble point, dissolved GOR,
 * formation volume factor (exponential compressibility correction above the bubble point,
 * Bo = Bob (pb/p)^A with A the Vazquez–Beggs numerator ×1e-5), density by material balance
 * ρo = (ρ_sto + Rs ρ_g,sc)/Bo, and viscosity (Beggs–Robinson; Vazquez–Beggs above the bubble point).
 * References: as the individual correlations.
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} o.api °API
 * @param {number} o.gammaG gas gravity, –
 * @param {number} o.Rsb total (produced) solution GOR, Sm³/Sm³
 * @param {'standing'|'vazquezBeggs'|'glaso'} [o.method='standing']
 * @param {number} [o.Tsep] K @param {number} [o.Psep] Pa
 * @param {number} [o.Pb] measured bubble point at T, Pa (overrides the correlation; Rs is then scaled so that
 *   Rs(Pb) = Rsb)
 * @returns {{ Pb: number, Rs: number, Bo: number, rho: number, mu: number, co: number, saturated: boolean }}
 *   Pa, Sm³/Sm³, m³/Sm³, kg/m³, Pa·s, 1/Pa
 */
export function oilProperties({ P, T, api, gammaG, Rsb, method = 'standing', Tsep, Psep, Pb }) {
  need('P', P, 1, 3e8);
  need('Rsb', Rsb, 0, 1e4);
  const common = { T, api, gammaG, method, Tsep, Psep };
  const pbCorr = bubblePoint({ Rs: Rsb, ...common });
  const pb = Pb ?? pbCorr;
  const saturated = P < pb;
  let Rs = Rsb;
  if (saturated) {
    // with a measured bubble point the correlation is shifted along the pressure axis
    const pEval = Pb != null ? P * pbCorr / pb : P;
    Rs = Math.min(Rsb, solutionGor({ P: pEval, ...common }));
  }
  const rhoSto = apiToSg(api) * RHO_WATER_SC;
  const gGas = method === 'vazquezBeggs' ? vbGasGravity(gammaG, api, Tsep, Psep) : gammaG;
  const coNum = vbCoNumerator({ T, Rsb, api, gammaG: gGas });
  const co = Math.max(coNum / (1e5 * (P / PSI)) / PSI, 1e-11);
  let Bo = oilFvfSaturated({ Rs, ...common });
  const muDead = deadOilViscosity({ T, api });
  let mu = liveOilViscosity({ muDead, Rs });
  if (!saturated && P > pb) {
    Bo *= (pb / P) ** Math.max(coNum * 1e-5, 0);
    mu = undersaturatedOilViscosity({ muOb: mu, P, Pb: pb });
  }
  const rho = (rhoSto + Rs * gasDensitySc(gammaG)) / Bo;
  return { Pb: pb, Rs, Bo, rho, mu, co, saturated };
}

/**
 * Gas PVT state at (P, T) from gas gravity: Z (DAK or Hall–Yarborough with Sutton + Wichert–Aziz), density,
 * formation volume factor and Lee–Gonzalez–Eakin viscosity. References: as the individual correlations.
 * @param {object} o
 * @param {number} o.P Pa @param {number} o.T K @param {number} o.gammaG –
 * @param {number} [o.yCO2=0] @param {number} [o.yH2S=0]
 * @param {'dak'|'hallYarborough'} [o.zMethod='dak']
 * @returns {{ Z: number, rho: number, Bg: number, mu: number, M: number, Tpr: number, Ppr: number }}
 *   –, kg/m³, m³/Sm³, Pa·s, kg/mol, –, –
 */
export function gasProperties({ P, T, gammaG, yCO2 = 0, yH2S = 0, zMethod = 'dak' }) {
  need('P', P, 1, 3e8);
  const pc = gasPseudoCriticals({ gammaG, yCO2, yH2S });
  const Z = zFactor({ P, T, gammaG, yCO2, yH2S, method: zMethod });
  const M = gammaG * M_AIR;
  const rho = P * M / (Z * R_GAS * T);
  return { Z, rho, Bg: gasFvf({ P, T, Z }), mu: gasViscosityLee({ T, rho: Math.min(rho, 800), M }), M,
    Tpr: T / pc.Tpc, Ppr: P / pc.Ppc };
}

/**
 * Water / brine PVT hook used by `fluidAt`: McCain density and formation volume factor; viscosity by Sharqawy
 * et al. below 311 K and McCain (with pressure factor) above. See water.js for references and validity.
 * @param {object} o
 * @param {number} o.P Pa @param {number} o.T K
 * @param {number} [o.S=35] dissolved solids, g/kg
 * @returns {{ rho: number, rhoSc: number, Bw: number, mu: number }} kg/m³, kg/m³, m³/Sm³, Pa·s
 */
export function waterProperties({ P, T, S = 35 }) {
  const d = brineDensity({ T, S, P });
  const mu = T < 311 ? seawaterViscosity({ T: Math.max(T, 273.15), S: Math.min(S, 150) })
    : brineViscosity({ T, S: Math.min(S, 260), P });
  return { rho: d.rho, rhoSc: d.rhoSc, Bw: d.Bw, mu };
}

/**
 * In-situ phase rates and properties of a black-oil stream at (P, T) — the single entry point for multiphase
 * hydraulics.
 *
 * Stream definition (stock-tank basis): oil rate `qo`, produced gas–oil ratio `GOR` (total gas at standard
 * conditions per stock-tank oil), water cut `WC` = qw/(qo + qw). Alternatively give `qgSc` and/or `qwSc`
 * explicitly (they override GOR and WC; needed for gas–water systems with qo = 0).
 * Gas dissolved in water, water vapour in gas and the change of free-gas gravity with pressure are neglected.
 * Liquid-mixture viscosity and interfacial tension are in-situ volume-fraction weighted (no emulsion model).
 * References: correlations named in `oilProperties`, `gasProperties`, `waterProperties`, `gasOilIft`,
 * `gasWaterIft`. Validity: the intersection of the correlation ranges; departures are listed in `warnings`.
 *
 * @param {object} o
 * @param {number} o.P absolute pressure, Pa
 * @param {number} o.T temperature, K
 * @param {number} [o.qo=0] stock-tank oil rate, Sm³/s
 * @param {number} [o.GOR=0] produced gas–oil ratio, Sm³/Sm³
 * @param {number} [o.WC=0] water cut at stock-tank conditions, – (0 ≤ WC < 1)
 * @param {number} [o.qgSc] total gas rate at standard conditions, Sm³/s (overrides GOR·qo)
 * @param {number} [o.qwSc] water rate at standard conditions, Sm³/s (overrides WC)
 * @param {number} o.gammaG gas gravity (air = 1), –
 * @param {number} [o.api=35] oil gravity, °API
 * @param {number} [o.salinity=35] water dissolved solids, g/kg
 * @param {number} [o.yCO2=0] @param {number} [o.yH2S=0] gas mole fractions
 * @param {'standing'|'vazquezBeggs'|'glaso'} [o.method='standing'] Rs/Pb/Bo correlation
 * @param {'dak'|'hallYarborough'} [o.zMethod='dak']
 * @param {number} [o.Pb] measured bubble point at T, Pa
 * @param {number} [o.Tsep] K @param {number} [o.Psep] Pa (Vazquez–Beggs gas-gravity normalisation)
 * @param {number} [o.D] pipe inner diameter, m — when given, superficial velocities are returned
 * @returns {{
 *   P: number, T: number, Pb: number, Rs: number, Bo: number, Bg: number, Bw: number, Z: number,
 *   qOil: number, qGas: number, qWater: number, qLiquid: number, qTotal: number,
 *   mOil: number, mGas: number, mWater: number, mTotal: number,
 *   rhoO: number, rhoG: number, rhoW: number, rhoL: number, rhoNoSlip: number,
 *   muO: number, muG: number, muW: number, muL: number, muNoSlip: number,
 *   sigmaOG: number, sigmaWG: number, sigmaL: number,
 *   lambdaL: number, waterFraction: number, qgFreeSc: number,
 *   vsl?: number, vsg?: number, vm?: number, warnings: string[] }}
 *   Pressures Pa; Rs Sm³/Sm³; Bo, Bg, Bw m³/Sm³; in-situ volumetric rates q* m³/s; mass rates m* kg/s;
 *   densities kg/m³; viscosities Pa·s; interfacial tensions N/m; lambdaL = no-slip liquid holdup (–);
 *   waterFraction = in-situ water volume fraction of the liquid (–); qgFreeSc = free gas at standard
 *   conditions, Sm³/s; vsl, vsg, vm = superficial liquid, gas and mixture velocities, m/s.
 */
export function fluidAt({ P, T, qo = 0, GOR = 0, WC = 0, qgSc, qwSc, gammaG, api = 35, salinity = 35,
  yCO2 = 0, yH2S = 0, method = 'standing', zMethod = 'dak', Pb, Tsep, Psep, D }) {
  need('P', P, 1e3, 3e8);
  need('T', T, 256, 500);
  need('qo', qo, 0, 1e3);
  need('GOR', GOR, 0, 1e6);
  need('WC', WC, 0, 1);
  const warnings = [];
  const qg = qgSc ?? qo * GOR;
  let qw = qwSc;
  if (qw == null) {
    if (WC >= 1) throw new RangeError('WC = 1 leaves the water rate undefined; give qwSc explicitly');
    qw = qo * WC / (1 - WC);
  }
  need('qgSc', qg, 0, 1e6);
  need('qwSc', qw, 0, 1e3);

  const gas = gasProperties({ P, T, gammaG, yCO2, yH2S, zMethod });
  if (gas.Tpr < 1.05) warnings.push(`Tpr = ${gas.Tpr.toFixed(3)} is below the Z-factor correlation range (1.05)`);
  const wat = waterProperties({ P, T, S: salinity });
  const rhoGsc = gasDensitySc(gammaG);

  let oil = null;
  let qgFree = qg;
  if (qo > 0) {
    oil = oilProperties({ P, T, api, gammaG, Rsb: qg / qo, method, Tsep, Psep, Pb });
    qgFree = Math.max(0, qg - qo * oil.Rs);
    if (T < 294) warnings.push('T is below the Beggs–Robinson viscosity range (294 K); oil viscosity is indicative');
    if (api < 16 || api > 58) warnings.push('API gravity is outside the black-oil correlation range (16–58)');
  } else {
    // dead-oil properties are still reported so that callers always get finite numbers
    oil = oilProperties({ P, T, api, gammaG, Rsb: 0, method: 'standing' });
  }

  const qOil = qo * oil.Bo;
  const qGas = qgFree * gas.Bg;
  const qWater = qw * wat.Bw;
  const qLiquid = qOil + qWater;
  const qTotal = qLiquid + qGas;
  const fw = qLiquid > 0 ? qWater / qLiquid : 0;
  const sigmaOG = gasOilIft({ P, T, api });
  const sigmaWG = gasWaterIft({ P, T });
  const rhoL = (1 - fw) * oil.rho + fw * wat.rho;
  const muL = (1 - fw) * oil.mu + fw * wat.mu;
  const sigmaL = (1 - fw) * sigmaOG + fw * sigmaWG;
  const lambdaL = qTotal > 0 ? qLiquid / qTotal : 0;
  const out = {
    P, T, Pb: oil.Pb, Rs: qo > 0 ? oil.Rs : 0, Bo: oil.Bo, Bg: gas.Bg, Bw: wat.Bw, Z: gas.Z,
    qOil, qGas, qWater, qLiquid, qTotal,
    mOil: qOil * oil.rho, mGas: qgFree * rhoGsc, mWater: qw * wat.rhoSc,
    mTotal: qOil * oil.rho + qgFree * rhoGsc + qw * wat.rhoSc,
    rhoO: oil.rho, rhoG: gas.rho, rhoW: wat.rho, rhoL,
    rhoNoSlip: lambdaL * rhoL + (1 - lambdaL) * gas.rho,
    muO: oil.mu, muG: gas.mu, muW: wat.mu, muL,
    muNoSlip: lambdaL * muL + (1 - lambdaL) * gas.mu,
    sigmaOG, sigmaWG, sigmaL, lambdaL, waterFraction: fw, qgFreeSc: qgFree, warnings,
  };
  if (D != null) {
    need('D', D, 1e-3, 10);
    const A = Math.PI * D * D / 4;
    out.vsl = qLiquid / A;
    out.vsg = qGas / A;
    out.vm = qTotal / A;
  }
  return out;
}
