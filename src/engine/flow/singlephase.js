// Single-phase pipeline hydraulics: friction factors, liquid and gas line pressure drop, classical gas
// transmission equations, coupled pressure–temperature marching, line sizing, machinery power, and water hammer
// (Joukowsky and method of characteristics).
// SI units: P absolute in Pa, T in K, lengths m, mass rate kg/s, volumetric rate m³/s, head m.

import { brent } from '../math/roots.js';
import { segmentTemperature } from './thermalpipe.js';
import { zFactor, gasViscosityLee, M_AIR, R_GAS, T_SC, P_SC } from './blackoil.js';
import { erosionalVelocity } from './erosion.js';

const G0 = 9.80665;
const PSI = 6894.757293168;
const INCH = 0.0254;
const MILE = 1609.344;
const FT3 = 0.028316846592;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * Reynolds number Re = ρ v D/μ. Reference: definition.
 * @param {object} o
 * @param {number} o.rho kg/m³ @param {number} o.v m/s @param {number} o.D m @param {number} o.mu Pa·s
 * @returns {number} –
 */
export function reynolds({ rho, v, D, mu }) {
  need('rho', rho, 1e-6, 1e5);
  need('D', D, 1e-6, 100);
  need('mu', mu, 1e-9, 1e6);
  return rho * Math.abs(need('v', v, -1e4, 1e4)) * D / mu;
}

function colebrook(Re, rr) {
  // Newton iteration on x = 1/√f: x + 2 log10(rr/3.7 + 2.51 x/Re) = 0, started from Haaland
  let x = -1.8 * Math.log10((rr / 3.7) ** 1.11 + 6.9 / Re);
  for (let it = 0; it < 50; it++) {
    const a = rr / 3.7 + 2.51 * x / Re;
    const F = x + 2 * Math.log10(a);
    const dF = 1 + 2 * 2.51 / (Re * a * Math.LN10);
    const dx = F / dF;
    x -= dx;
    if (Math.abs(dx) < 1e-14 * Math.abs(x)) break;
  }
  return 1 / (x * x);
}

/**
 * Darcy–Weisbach friction factor.
 *  - 'colebrook' (default): laminar f = 64/Re for Re < 2300; otherwise the Colebrook–White equation
 *    1/√f = −2 log10[(ε/D)/3.7 + 2.51/(Re √f)] solved to machine precision by Newton iteration.
 *  - 'haaland': 1/√f = −1.8 log10[((ε/D)/3.7)^1.11 + 6.9/Re] (explicit, within about 2 % of Colebrook).
 *  - 'churchill': all-regime f = 8[(8/Re)^12 + (A + B)^−1.5]^(1/12), A = {2.457 ln[1/((7/Re)^0.9 + 0.27 ε/D)]}^16,
 *    B = (37 530/Re)^16 — continuous through transition, preferred inside marching solvers.
 *  - 'laminar': 64/Re.
 * References: Colebrook, J. ICE 11 (1939) 133; Haaland, J. Fluids Eng. 105 (1983) 89; Churchill, Chem. Eng. 84(24)
 * (1977) 91; Moody, Trans. ASME 66 (1944) 671.
 * Validity: Newtonian fluid, circular pipe, ε/D ≤ 0.05. The 2300–4000 transition band is not physically resolved.
 * @param {object} o
 * @param {number} o.Re Reynolds number, – (> 0)
 * @param {number} [o.relRough=0] relative roughness ε/D, –
 * @param {'colebrook'|'haaland'|'churchill'|'laminar'} [o.method='colebrook']
 * @returns {number} Darcy friction factor, – (Fanning = Darcy/4)
 */
export function frictionFactor({ Re, relRough = 0, method = 'colebrook' }) {
  need('Re', Re, 1e-9, 1e12);
  need('relRough', relRough, 0, 0.1);
  switch (method) {
    case 'laminar': return 64 / Re;
    case 'colebrook': return Re < 2300 ? 64 / Re : colebrook(Re, relRough);
    case 'haaland': return Re < 2300 ? 64 / Re : (-1.8 * Math.log10((relRough / 3.7) ** 1.11 + 6.9 / Re)) ** -2;
    case 'churchill': {
      const A = (2.457 * Math.log(1 / ((7 / Re) ** 0.9 + 0.27 * relRough))) ** 16;
      const B = (37530 / Re) ** 16;
      return 8 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12);
    }
    default: throw new RangeError(`method = ${method} is not one of colebrook, haaland, churchill, laminar`);
  }
}

/**
 * Minor (fitting) pressure loss ΔP = ΣK ρ v²/2 and the equivalent pipe length L_eq = ΣK D/f.
 * Reference: Crane Technical Paper 410; Idelchik, "Handbook of Hydraulic Resistance".
 * @param {object} o
 * @param {number} o.K sum of resistance coefficients, – @param {number} o.rho kg/m³ @param {number} o.v m/s
 * @param {number} [o.D] m @param {number} [o.f] Darcy friction factor (both needed for the equivalent length)
 * @returns {{ dP: number, head: number, Leq: number|null }} Pa, m of fluid, m
 */
export function minorLoss({ K, rho, v, D, f }) {
  need('K', K, 0, 1e6);
  need('rho', rho, 1e-6, 1e5);
  need('v', v, -1e4, 1e4);
  return { dP: K * rho * v * v / 2, head: K * v * v / (2 * G0), Leq: D != null && f != null ? K * D / f : null };
}

function normaliseProfile(profile, L, dz) {
  if (profile) {
    if (!Array.isArray(profile) || profile.length < 2) throw new RangeError('profile needs at least two points');
    for (let i = 1; i < profile.length; i++) {
      if (!(profile[i].x > profile[i - 1].x)) throw new RangeError(`profile x must increase strictly (point ${i})`);
    }
    return profile.map((p) => ({ x: p.x, z: p.z ?? 0 }));
  }
  need('L', L, 1e-9, 1e8);
  return [{ x: 0, z: 0 }, { x: L, z: dz ?? 0 }];
}

/**
 * Pressure drop of an incompressible (liquid) line along an elevation profile:
 * ΔP = f (L/D) ρ v²/2 + ΣK ρ v²/2 + ρ g Δz (Darcy–Weisbach + minor losses + static head).
 * Reference: Darcy–Weisbach; Bernoulli with losses. Validity: constant density and viscosity, full pipe.
 * @param {object} o
 * @param {number} o.Q volumetric flow rate, m³/s
 * @param {number} o.D inner diameter, m
 * @param {number} o.rho kg/m³ @param {number} o.mu Pa·s
 * @param {number} [o.rough=4.5e-5] absolute roughness, m
 * @param {{ x: number, z?: number }[]} [o.profile] route points (distance m, elevation m), or give L and dz
 * @param {number} [o.L] length, m @param {number} [o.dz=0] outlet minus inlet elevation, m
 * @param {number} [o.K=0] sum of minor-loss coefficients
 * @param {number} [o.Pin] inlet pressure, Pa (to report the pressure profile)
 * @param {string} [o.method='colebrook'] friction-factor method
 * @returns {{ dP: number, dPfriction: number, dPminor: number, dPelevation: number, v: number, Re: number,
 *   f: number, gradient: number, x: number[], P: number[]|null }}
 *   pressure drops Pa (positive = loss in the flow direction), friction gradient Pa/m
 */
export function liquidLinePressureDrop({ Q, D, rho, mu, rough = 4.5e-5, profile, L, dz = 0, K = 0, Pin, method = 'colebrook' }) {
  need('Q', Q, -1e4, 1e4);
  need('D', D, 1e-4, 50);
  need('rough', rough, 0, 0.05 * D);
  const pts = normaliseProfile(profile, L, dz);
  const v = Q / (Math.PI * D * D / 4);
  const Re = reynolds({ rho, v, D, mu });
  const f = Re > 0 ? frictionFactor({ Re, relRough: rough / D, method }) : 0;
  const grad = f * rho * v * Math.abs(v) / (2 * D);
  const len = pts[pts.length - 1].x - pts[0].x;
  const dPf = grad * len;
  const dPm = K * rho * v * Math.abs(v) / 2;
  const dPe = rho * G0 * (pts[pts.length - 1].z - pts[0].z);
  let P = null;
  if (Pin != null) {
    // minor losses are distributed uniformly along the route for the reported profile
    P = pts.map((p) => Pin - (grad + dPm / len) * (p.x - pts[0].x) - rho * G0 * (p.z - pts[0].z));
  }
  return { dP: dPf + dPm + dPe, dPfriction: dPf, dPminor: dPm, dPelevation: dPe, v, Re, f, gradient: grad,
    x: pts.map((p) => p.x), P };
}

function elevationTerms(M, Z, T, L, dz) {
  // s = 2 g M Δz/(Z R T); equivalent length L_e = L (e^s − 1)/s
  const s = 2 * G0 * M * dz / (Z * R_GAS * T);
  const Le = Math.abs(s) < 1e-9 ? L * (1 + s / 2) : L * Math.expm1(s) / s;
  return { s, Le, es: Math.exp(s) };
}

/**
 * Steady isothermal gas flow through a pipe (general flow equation, derived from the momentum balance
 * dP/dx = −f ρ v²/(2D) − ρ g dz/dx with ρ = P M/(Z R T) at average Z and T, kinetic-energy term neglected):
 *   ṁ = (π/4) E √[ (P1² − e^s P2²) M D⁵ / (f Z R T L_e) ],  s = 2 g M Δz/(Z R T),  L_e = L (e^s − 1)/s.
 * The Darcy friction factor is taken from `f`, or computed from Reynolds number and roughness (iterating on the
 * flow rate), or set by a classical transmission equation through `equation`:
 *   'aga' — AGA fully turbulent: 1/√f_Fanning = 4 log10(3.7 D/ε) (rough-pipe law, no drag factor);
 *   'weymouth' — f = 0.032/D_in^(1/3) (D in inches), the friction law underlying the Weymouth equation.
 * Either the outlet pressure `P2` or the mass rate `mDot` must be given; the other is returned.
 * Reference: e.g. Menon, "Gas Pipeline Hydraulics" (2005) ch. 2; AGA NB-13 (1965); Weymouth, Trans. ASME 34 (1912).
 * Validity: steady, isothermal, single-phase gas; velocity well below sonic (kinetic term ignored).
 * @param {object} o
 * @param {number} o.P1 inlet pressure, Pa
 * @param {number} [o.P2] outlet pressure, Pa @param {number} [o.mDot] mass rate, kg/s
 * @param {number} o.D inner diameter, m @param {number} o.L length, m
 * @param {number} o.T average flowing temperature, K
 * @param {number} o.Z average compressibility factor, –
 * @param {number} [o.M] molar mass, kg/mol (or give gammaG) @param {number} [o.gammaG] gas gravity, –
 * @param {number} [o.dz=0] outlet minus inlet elevation, m
 * @param {number} [o.E=1] pipeline efficiency factor, –
 * @param {number} [o.f] Darcy friction factor
 * @param {'general'|'aga'|'weymouth'} [o.equation='general']
 * @param {number} [o.rough=1.5e-5] absolute roughness, m @param {number} [o.mu] gas viscosity, Pa·s ('general')
 * @returns {{ mDot: number, P2: number, f: number, qStd: number, v1: number, Re: number|null }}
 *   kg/s, Pa, –, standard volumetric rate Sm³/s (60 °F, 101.325 kPa), inlet velocity m/s
 */
export function gasFlowGeneral({ P1, P2, mDot, D, L, T, Z, M, gammaG, dz = 0, E = 1, f, equation = 'general',
  rough = 1.5e-5, mu }) {
  need('P1', P1, 1, 1e9);
  need('D', D, 1e-4, 20);
  need('L', L, 1e-6, 1e8);
  need('T', T, 50, 2000);
  need('Z', Z, 0.05, 5);
  need('E', E, 1e-3, 1.5);
  const Mm = M ?? (need('gammaG', gammaG, 0.05, 5) * M_AIR);
  need('M', Mm, 1e-3, 1);
  if ((P2 == null) === (mDot == null)) throw new RangeError('give exactly one of P2 and mDot');
  const { Le, es } = elevationTerms(Mm, Z, T, L, dz);
  const A = Math.PI * D * D / 4;
  const fFixed = f ?? (equation === 'aga' ? 4 / (4 * Math.log10(3.7 * D / need('rough', rough, 1e-9, 0.05 * D))) ** 2
    : equation === 'weymouth' ? 0.032 / Math.cbrt(D / INCH) : null);
  if (fFixed == null && equation !== 'general') throw new RangeError(`equation = ${equation} is not one of general, aga, weymouth`);
  const fOf = (m) => {
    if (fFixed != null) return fFixed;
    need('mu', mu, 1e-7, 1e-2);
    const Re = Math.max(4 * Math.abs(m) / (Math.PI * D * mu), 1e-6);
    return frictionFactor({ Re, relRough: rough / D });
  };
  const flow = (p2, ff) => {
    const d = P1 * P1 - es * p2 * p2;
    return Math.sign(d) * A * E * Math.sqrt(Math.abs(d) * Mm * D / (ff * Z * R_GAS * T * Le));
  };
  let m = mDot;
  let p2 = P2;
  let ff;
  if (mDot == null) {
    need('P2', P2, 0, 1e9);
    ff = fFixed ?? 0.015;
    m = flow(P2, ff);
    for (let it = 0; it < 100 && fFixed == null; it++) {
      const fn = fOf(m);
      const mn = flow(P2, fn);
      const done = Math.abs(mn - m) <= 1e-13 * Math.abs(mn);
      ff = fn;
      m = mn;
      if (done) break;
    }
  } else {
    need('mDot', mDot, 0, 1e7);
    ff = fOf(mDot);
    const p2sq = (P1 * P1 - (mDot / (A * E)) ** 2 * ff * Z * R_GAS * T * Le / (Mm * D)) / es;
    if (!(p2sq > 0)) throw new RangeError('the specified flow rate cannot be delivered: outlet pressure would be ≤ 0');
    p2 = Math.sqrt(p2sq);
  }
  const rho1 = P1 * Mm / (Z * R_GAS * T);
  return { mDot: m, P2: p2, f: ff, qStd: m / (Mm * P_SC / (R_GAS * T_SC)), v1: m / (rho1 * A),
    Re: mu ? 4 * Math.abs(m) / (Math.PI * D * mu) : null };
}

// Published oilfield-unit forms: Q (scf/d) = C E (Tb/Pb)^a [(P1² − e^s P2²)/(γ^b T Le Z)]^m D^n,
// P psia, T °R, L miles, D inches.
const TRANSMISSION = {
  weymouth: { C: 433.5, a: 1, b: 1, m: 0.5, n: 2.667 },
  panhandleA: { C: 435.87, a: 1.0788, b: 0.8539, m: 0.5394, n: 2.6182 },
  panhandleB: { C: 737, a: 1.02, b: 0.961, m: 0.51, n: 2.53 },
};

/**
 * Classical gas-transmission equations in their published oilfield-unit form
 *   Q = C E (T_b/P_b)^a [(P1² − e^s P2²)/(γ^b T L_e Z)]^m D^n   (Q scf/d, P psia, T °R, L miles, D inches):
 *   Weymouth:    C = 433.5,  a = 1,      b = 1,      m = 0.5,    n = 2.667
 *   Panhandle A: C = 435.87, a = 1.0788, b = 0.8539, m = 0.5394, n = 2.6182
 *   Panhandle B: C = 737,    a = 1.02,   b = 0.961,  m = 0.51,   n = 2.53
 * with the elevation correction s = 2 g M Δz/(Z R T) (= 0.0375 γ ΔH/(T Z) in ft, °R) and L_e = L (e^s − 1)/s.
 * Give either `P2` (returns the flow) or `qStd` (returns the outlet pressure).
 * References: Weymouth (1912); Panhandle Eastern Pipe Line Co. (1940s) and the modified Panhandle (1956) as given
 * in GPSA Engineering Data Book §17 and Menon (2005) ch. 2.
 * Validity: fully turbulent, large-diameter, high-pressure transmission lines; Weymouth is conservative for
 * D > ~0.3 m; Panhandle A for Re 5e6–1.1e7, Panhandle B for Re 4e6–4e7. Efficiency E typically 0.88–0.95.
 * @param {object} o
 * @param {'weymouth'|'panhandleA'|'panhandleB'} o.equation
 * @param {number} o.P1 Pa @param {number} [o.P2] Pa @param {number} [o.qStd] Sm³/s at (Tb, Pb)
 * @param {number} o.D m @param {number} o.L m @param {number} o.T K @param {number} o.Z – @param {number} o.gammaG –
 * @param {number} [o.dz=0] m @param {number} [o.E=1] –
 * @param {number} [o.Tb=288.7056] base temperature, K @param {number} [o.Pb=101325] base pressure, Pa
 * @returns {{ qStd: number, P2: number, mDot: number }} Sm³/s at base conditions, Pa, kg/s
 */
export function gasTransmission({ equation, P1, P2, qStd, D, L, T, Z, gammaG, dz = 0, E = 1, Tb = T_SC, Pb = P_SC }) {
  const c = TRANSMISSION[equation];
  if (!c) throw new RangeError(`equation = ${equation} is not one of weymouth, panhandleA, panhandleB`);
  need('P1', P1, 1, 1e9);
  need('D', D, 1e-3, 20);
  need('L', L, 1e-6, 1e8);
  need('T', T, 50, 2000);
  need('Z', Z, 0.05, 5);
  need('gammaG', gammaG, 0.05, 5);
  need('E', E, 1e-3, 1.5);
  if ((P2 == null) === (qStd == null)) throw new RangeError('give exactly one of P2 and qStd');
  const { Le, es } = elevationTerms(gammaG * M_AIR, Z, T, L, dz);
  const coef = c.C * E * ((Tb * 1.8) / (Pb / PSI)) ** c.a * (D / INCH) ** c.n;
  const denom = gammaG ** c.b * (T * 1.8) * (Le / MILE) * Z;
  const toSI = FT3 / 86400; // scf/d → m³/s at base conditions
  let q = qStd;
  let p2 = P2;
  if (qStd == null) {
    need('P2', P2, 0, 1e9);
    const d = (P1 * P1 - es * P2 * P2) / (PSI * PSI);
    q = Math.sign(d) * coef * (Math.abs(d) / denom) ** c.m * toSI;
  } else {
    need('qStd', qStd, 0, 1e9);
    const d = (qStd / toSI / coef) ** (1 / c.m) * denom * PSI * PSI;
    const p2sq = (P1 * P1 - d) / es;
    if (!(p2sq > 0)) throw new RangeError('the specified flow rate cannot be delivered: outlet pressure would be ≤ 0');
    p2 = Math.sqrt(p2sq);
  }
  return { qStd: q, P2: p2, mDot: q * gammaG * M_AIR * Pb / (R_GAS * Tb) };
}

/**
 * Joule–Thomson coefficient of a natural gas from the Z-factor correlation,
 * μ_JT = (∂T/∂P)_h = R T² (∂Z/∂T)_P / (P M cp), with the derivative by central difference of `zFactor`.
 * Reference: thermodynamic identity μ_JT = [T(∂v/∂T)_P − v]/cp with v = Z R T/(P M).
 * Validity: as the Z-factor correlation (Tpr > 1.05).
 * @param {object} o
 * @param {number} o.P Pa @param {number} o.T K @param {number} o.gammaG – @param {number} o.cp J/(kg·K)
 * @param {number} [o.yCO2=0] @param {number} [o.yH2S=0]
 * @returns {number} K/Pa
 */
export function gasJouleThomson({ P, T, gammaG, cp, yCO2 = 0, yH2S = 0 }) {
  need('cp', cp, 100, 1e5);
  const h = 0.5;
  const dZ = (zFactor({ P, T: T + h, gammaG, yCO2, yH2S }) - zFactor({ P, T: T - h, gammaG, yCO2, yH2S })) / (2 * h);
  return R_GAS * T * T * dZ / (P * gammaG * M_AIR * cp);
}

/**
 * Non-isothermal steady marching solver for a single-phase line along an elevation profile, coupling the momentum
 * balance and the energy balance segment by segment.
 * Each profile interval is split into `nSub` steps; on every step the properties are evaluated at the step-average
 * pressure and temperature (two fixed-point passes), the pressure change is
 *   ΔP = −[f ρ v²/(2D) + ρ g sinθ] Δx / (1 − ρ v²/P)        (gas; the denominator is the kinetic-energy term)
 *   ΔP = −[f ρ v²/(2D) + ρ g sinθ] Δx                        (liquid)
 * and the temperature follows the analytic exponential solution of `segmentTemperature` with the local U-value,
 * the Joule–Thomson coefficient and the elevation term.
 * Gas: Z by `zFactor`, viscosity by Lee–Gonzalez–Eakin, μ_JT by `gasJouleThomson` (or as given). Liquid: constant
 * ρ, μ (or μ(T) function), μ_JT = −1/(ρ cp).
 * Reference: conservation of momentum and energy for 1-D pipe flow (e.g. Bird, Stewart & Lightfoot ch. 15);
 * friction by Churchill (1977) by default for continuity across regimes.
 * Validity: single phase, subsonic (throws when ρv²/P reaches 0.5), U referred to the inner diameter.
 * @param {object} o
 * @param {{ type: 'gas', gammaG: number, cp: number, yCO2?: number, yH2S?: number, muJT?: number } |
 *         { type: 'liquid', rho: number, mu: number|function(number):number, cp: number }} o.fluid
 * @param {number} o.mDot mass flow rate, kg/s
 * @param {number} o.Pin inlet pressure, Pa @param {number} o.Tin inlet temperature, K
 * @param {number} o.D inner diameter, m @param {number} [o.rough=4.5e-5] m
 * @param {{ x: number, z?: number }[]} o.profile route points
 * @param {number|number[]} o.U overall heat-transfer coefficient on πD, W/(m²·K); scalar or per profile interval
 * @param {number|number[]} o.Tamb ambient temperature, K; scalar or per profile interval
 * @param {number} [o.nSub=10] steps per profile interval
 * @param {string} [o.method='churchill'] friction-factor method
 * @returns {{ x: number[], z: number[], P: number[], T: number[], v: number[], rho: number[], Pout: number,
 *   Tout: number, dP: number }} values at the profile points; Pa, K, m/s, kg/m³
 */
export function lineProfile({ fluid, mDot, Pin, Tin, D, rough = 4.5e-5, profile, U, Tamb, nSub = 10, method = 'churchill' }) {
  need('mDot', mDot, 1e-9, 1e6);
  need('Pin', Pin, 1, 1e9);
  need('Tin', Tin, 50, 2000);
  need('D', D, 1e-4, 20);
  need('nSub', nSub, 1, 1e5);
  const pts = normaliseProfile(profile);
  if (!fluid || (fluid.type !== 'gas' && fluid.type !== 'liquid')) throw new RangeError("fluid.type must be 'gas' or 'liquid'");
  const gas = fluid.type === 'gas';
  const A = Math.PI * D * D / 4;
  const G = mDot / A;
  const props = (P, T) => {
    if (gas) {
      const Z = zFactor({ P, T, gammaG: fluid.gammaG, yCO2: fluid.yCO2 ?? 0, yH2S: fluid.yH2S ?? 0 });
      const M = fluid.gammaG * M_AIR;
      const rho = P * M / (Z * R_GAS * T);
      const muJT = fluid.muJT ?? gasJouleThomson({ P, T, gammaG: fluid.gammaG, cp: fluid.cp, yCO2: fluid.yCO2 ?? 0, yH2S: fluid.yH2S ?? 0 });
      return { rho, mu: gasViscosityLee({ T, rho: Math.min(rho, 800), M }), muJT };
    }
    const mu = typeof fluid.mu === 'function' ? fluid.mu(T) : fluid.mu;
    return { rho: fluid.rho, mu, muJT: -1 / (fluid.rho * fluid.cp) };
  };
  const out = { x: [pts[0].x], z: [pts[0].z], P: [Pin], T: [Tin], v: [], rho: [] };
  let P = Pin;
  let T = Tin;
  let p0 = props(P, T);
  out.v.push(G / p0.rho);
  out.rho.push(p0.rho);
  for (let i = 1; i < pts.length; i++) {
    const Lseg = pts[i].x - pts[i - 1].x;
    const len = Math.hypot(Lseg, pts[i].z - pts[i - 1].z); // length along the pipe
    const ds = len / nSub;
    const sinT = (pts[i].z - pts[i - 1].z) / len;
    const Ui = Array.isArray(U) ? U[i - 1] : U;
    const Ta = Array.isArray(Tamb) ? Tamb[i - 1] : Tamb;
    for (let k = 0; k < nSub; k++) {
      let Pn = P;
      let Tn = T;
      for (let pass = 0; pass < 3; pass++) {
        const pr = props(0.5 * (P + Pn), 0.5 * (T + Tn));
        const v = G / pr.rho;
        const Re = pr.rho * v * D / pr.mu;
        const f = frictionFactor({ Re, relRough: rough / D, method });
        let grad = -(f * pr.rho * v * v / (2 * D) + pr.rho * G0 * sinT);
        if (gas) {
          const kin = pr.rho * v * v / (0.5 * (P + Pn));
          if (kin >= 0.5) throw new Error('gas velocity approaches sonic conditions; the line is choked for this flow rate');
          grad /= (1 - kin);
        }
        Pn = P + grad * ds;
        if (!(Pn > 0)) throw new Error('pressure falls to zero along the line: the flow rate cannot be delivered');
        Tn = segmentTemperature({ Tin: T, Tamb: Ta, L: ds, mDot, cp: fluid.cp, U: Ui, D, muJT: pr.muJT, dPdx: grad,
          dzdx: sinT }).Tout;
      }
      P = Pn;
      T = Tn;
    }
    p0 = props(P, T);
    out.x.push(pts[i].x);
    out.z.push(pts[i].z);
    out.P.push(P);
    out.T.push(T);
    out.v.push(G / p0.rho);
    out.rho.push(p0.rho);
  }
  return { ...out, Pout: P, Tout: T, dP: Pin - P };
}

/**
 * Noise-related velocity limit for gas piping in the form V = a (1/ρ)^n, capped at `cap`.
 * The defaults a = 175, n = 0.43, cap = 60 m/s follow the form used in NORSOK P-001 (process design) for gas lines;
 * they are exposed as arguments because project specifications differ — confirm against the governing standard.
 * @param {object} o
 * @param {number} o.rho gas density, kg/m³
 * @param {number} [o.a=175] @param {number} [o.n=0.43] @param {number} [o.cap=60] m/s
 * @returns {number} velocity limit, m/s
 */
export function noiseVelocityLimit({ rho, a = 175, n = 0.43, cap = 60 }) {
  need('rho', rho, 1e-3, 2000);
  need('a', a, 1, 1e4);
  need('n', n, 0, 2);
  need('cap', cap, 0.1, 1e3);
  return Math.min(cap, a * (1 / rho) ** n);
}

/**
 * Velocity screening of a line against the API RP 14E erosional velocity and the gas noise limit (hooks into
 * erosion.js and `noiseVelocityLimit`).
 * Reference: API RP 14E §2.5.a; see `noiseVelocityLimit`.
 * @param {object} o
 * @param {number} o.v actual (mixture) velocity, m/s @param {number} o.rho (mixture) density, kg/m³
 * @param {number} [o.C=100] API RP 14E C-factor in its native units, (lb/ft³)^½·ft/s
 * @param {boolean} [o.gas=false] also apply the noise limit
 * @returns {{ erosional: number, erosionalRatio: number, noise: number|null, noiseRatio: number|null, ok: boolean }}
 *   limits m/s and utilisation ratios v/limit
 */
export function velocityChecks({ v, rho, C = 100, gas = false }) {
  need('v', v, 0, 1e4);
  const ve = erosionalVelocity({ rho, C });
  const vn = gas ? noiseVelocityLimit({ rho }) : null;
  return { erosional: ve, erosionalRatio: v / ve, noise: vn, noiseRatio: vn ? v / vn : null,
    ok: v <= ve && (vn == null || v <= vn) };
}

/**
 * Minimum inner diameter of a liquid or (locally incompressible) gas line for a maximum velocity and/or a maximum
 * frictional pressure gradient, and the smallest candidate size that satisfies both.
 * Reference: Darcy–Weisbach with Colebrook–White; velocity criterion v = Q/A.
 * @param {object} o
 * @param {number} o.Q volumetric flow rate at flowing conditions, m³/s
 * @param {number} o.rho kg/m³ @param {number} o.mu Pa·s
 * @param {number} [o.vMax] maximum velocity, m/s
 * @param {number} [o.gradMax] maximum friction gradient, Pa/m
 * @param {number} [o.rough=4.5e-5] m
 * @param {number[]} [o.candidates] available inner diameters, m
 * @returns {{ Dvelocity: number|null, Dgradient: number|null, Drequired: number, governing: string,
 *   selected: number|null, v: number, gradient: number }} diameters m; velocity and gradient at the selected (or
 *   required) diameter
 */
export function sizeLine({ Q, rho, mu, vMax, gradMax, rough = 4.5e-5, candidates }) {
  need('Q', Q, 1e-12, 1e4);
  need('rho', rho, 1e-6, 1e5);
  need('mu', mu, 1e-9, 1e6);
  if (vMax == null && gradMax == null) throw new RangeError('give vMax and/or gradMax');
  const gradOf = (D) => {
    const v = Q / (Math.PI * D * D / 4);
    const Re = rho * v * D / mu;
    return frictionFactor({ Re, relRough: Math.min(rough / D, 0.1) }) * rho * v * v / (2 * D);
  };
  const Dv = vMax != null ? Math.sqrt(4 * Q / (Math.PI * need('vMax', vMax, 1e-6, 1e4))) : null;
  let Dg = null;
  if (gradMax != null) {
    need('gradMax', gradMax, 1e-9, 1e9);
    // the friction gradient falls steeply with D; search on ln D between 1 mm and 20 m
    const lo = 1e-3;
    const hi = 20;
    if (gradOf(hi) > gradMax) throw new RangeError('gradMax cannot be met with a diameter up to 20 m');
    Dg = gradOf(lo) <= gradMax ? lo
      : Math.exp(brent((lnD) => Math.log(gradOf(Math.exp(lnD)) / gradMax), Math.log(lo), Math.log(hi), { tol: 1e-12 }));
  }
  const Dreq = Math.max(Dv ?? 0, Dg ?? 0);
  let selected = null;
  if (candidates) {
    const sorted = candidates.slice().sort((a, b) => a - b);
    selected = sorted.find((d) => d >= Dreq * (1 - 1e-12)
      && (gradMax == null || gradOf(d) <= gradMax * (1 + 1e-9))) ?? null;
  }
  const Duse = selected ?? Dreq;
  return { Dvelocity: Dv, Dgradient: Dg, Drequired: Dreq, governing: (Dv ?? 0) >= (Dg ?? 0) ? 'velocity' : 'gradient',
    selected, v: Q / (Math.PI * Duse * Duse / 4), gradient: gradOf(Duse) };
}

/**
 * Pump hydraulic and shaft power: P_h = ρ g H Q = ΔP Q, P_shaft = P_h/η.
 * Reference: energy equation for an incompressible fluid.
 * @param {object} o
 * @param {number} o.Q m³/s @param {number} [o.dP] differential pressure, Pa
 * @param {number} [o.head] differential head, m (with rho, alternative to dP) @param {number} [o.rho] kg/m³
 * @param {number} [o.efficiency=0.75] –
 * @returns {{ hydraulic: number, shaft: number, dP: number }} W, W, Pa
 */
export function pumpPower({ Q, dP, head, rho, efficiency = 0.75 }) {
  need('Q', Q, 0, 1e4);
  need('efficiency', efficiency, 1e-3, 1);
  const d = dP ?? need('rho', rho, 1e-3, 1e5) * G0 * need('head', head, -1e6, 1e6);
  return { hydraulic: d * Q, shaft: d * Q / efficiency, dP: d };
}

/**
 * Compressor polytropic head, discharge temperature and power:
 *   (n − 1)/n = (k − 1)/(k η_p),  H_p = Z_avg R T1/M · n/(n − 1) · [(P2/P1)^((n−1)/n) − 1]  (J/kg),
 *   T2 = T1 (P2/P1)^((n−1)/n),  gas power = ṁ H_p/η_p,  shaft power = gas power/η_mech.
 * Reference: GPSA Engineering Data Book §13 (polytropic calculation); ASME PTC 10 (ideal-gas exponent form).
 * Validity: single section without intercooling, average Z, constant isentropic exponent k.
 * @param {object} o
 * @param {number} o.mDot kg/s @param {number} o.P1 Pa @param {number} o.P2 Pa @param {number} o.T1 K
 * @param {number} o.M molar mass, kg/mol @param {number} o.k isentropic exponent cp/cv, –
 * @param {number} [o.Z=1] average compressibility factor
 * @param {number} [o.polytropicEfficiency=0.78] @param {number} [o.mechanicalEfficiency=0.98]
 * @returns {{ head: number, T2: number, n: number, gasPower: number, shaftPower: number }} J/kg, K, –, W, W
 */
export function compressorPolytropic({ mDot, P1, P2, T1, M, k, Z = 1, polytropicEfficiency = 0.78, mechanicalEfficiency = 0.98 }) {
  need('mDot', mDot, 0, 1e6);
  need('P1', P1, 1, 1e9);
  need('P2', P2, P1, 1e9);
  need('T1', T1, 50, 2000);
  need('M', M, 1e-3, 1);
  need('k', k, 1.0001, 2);
  need('Z', Z, 0.05, 5);
  need('polytropicEfficiency', polytropicEfficiency, 0.05, 1);
  need('mechanicalEfficiency', mechanicalEfficiency, 0.05, 1);
  const e = (k - 1) / (k * polytropicEfficiency);
  const r = P2 / P1;
  const head = Z * R_GAS * T1 / M * (r ** e - 1) / e;
  const gasPower = mDot * head / polytropicEfficiency;
  return { head, T2: T1 * r ** e, n: 1 / (1 - e), gasPower, shaftPower: gasPower / mechanicalEfficiency };
}

/**
 * Pressure-wave speed in a liquid-filled elastic pipe (thin-walled):
 *   a = √(K/ρ) / √(1 + (K/E)(D/e) c1),
 * c1 = 1 − ν/2 (anchored at the upstream end only), 1 − ν² (anchored against axial movement throughout),
 * 1 (expansion joints throughout).
 * Reference: Wylie & Streeter, "Fluid Transients in Systems" (1993) ch. 1; Korteweg (1878).
 * Validity: D/e > ~25 (thin wall), single-phase liquid without free gas.
 * @param {object} o
 * @param {number} o.K fluid bulk modulus, Pa @param {number} o.rho fluid density, kg/m³
 * @param {number} o.D inner diameter, m @param {number} o.e wall thickness, m
 * @param {number} o.E Young's modulus of the wall, Pa @param {number} [o.nu=0.3] Poisson's ratio
 * @param {'upstream'|'throughout'|'expansionJoints'} [o.restraint='throughout']
 * @returns {number} wave speed, m/s
 */
export function waveSpeed({ K, rho, D, e, E, nu = 0.3, restraint = 'throughout' }) {
  need('K', K, 1e3, 1e12);
  need('rho', rho, 1e-3, 3e4);
  need('D', D, 1e-4, 50);
  need('e', e, 1e-6, 10);
  need('E', E, 1e5, 1e13);
  need('nu', nu, 0, 0.5);
  const c1 = restraint === 'upstream' ? 1 - nu / 2 : restraint === 'throughout' ? 1 - nu * nu
    : restraint === 'expansionJoints' ? 1 : NaN;
  if (Number.isNaN(c1)) throw new RangeError(`restraint = ${restraint} is not one of upstream, throughout, expansionJoints`);
  return Math.sqrt(K / rho) / Math.sqrt(1 + (K / E) * (D / e) * c1);
}

/**
 * Joukowsky surge for a rapid velocity change: ΔP = ρ a Δv, Δh = a Δv/g; applies in full when the closure time is
 * shorter than the pipe period 2L/a. For slower linear closure a first estimate is ΔP·(2L/a)/t_c.
 * Reference: Joukowsky (1898); Wylie & Streeter (1993) ch. 1.
 * @param {object} o
 * @param {number} o.rho kg/m³ @param {number} o.a wave speed, m/s @param {number} o.dv velocity change, m/s
 * @param {number} [o.L] pipe length, m @param {number} [o.tClose] closure time, s
 * @returns {{ dP: number, dHead: number, period: number|null, rapid: boolean|null, dPslow: number|null }}
 *   Pa, m, s (2L/a), whether the closure is rapid, reduced estimate for slow closure Pa
 */
export function joukowsky({ rho, a, dv, L, tClose }) {
  need('rho', rho, 1e-3, 3e4);
  need('a', a, 1, 1e4);
  need('dv', dv, -1e3, 1e3);
  const dP = rho * a * Math.abs(dv);
  const period = L != null ? 2 * L / a : null;
  const rapid = period != null && tClose != null ? tClose <= period : null;
  return { dP, dHead: a * Math.abs(dv) / G0, period, rapid,
    dPslow: rapid === false ? dP * period / tClose : rapid === true ? dP : null };
}

/**
 * Water-hammer transient in a single uniform pipe by the method of characteristics (fixed grid, Δt = Δx/a).
 * Compatibility equations (Wylie & Streeter):
 *   C⁺: H_P = C_P − B Q_P,  C_P = H_A + B Q_A − R Q_A|Q_A|
 *   C⁻: H_P = C_M + B Q_P,  C_M = H_B − B Q_B + R Q_B|Q_B|,  B = a/(gA),  R = f Δx/(2 g D A²).
 * Boundaries: upstream reservoir of head `H0` (constant or function of time); downstream either a valve
 * discharging to a fixed back head with opening law τ(t) (Q = τ Q0 √(ΔH/ΔH0)), or a dead end (Q = 0).
 * Default valve law: τ = (1 − t/t_c)^m for t < t_c, 0 afterwards.
 * Initial condition: steady flow Q0 with the friction grade line (valve), or rest at `Hinit` (dead end).
 * Friction: quasi-steady Darcy–Weisbach with constant f (unsteady friction and column separation are not
 * modelled; heads below vapour head are flagged in `cavitation`).
 * Reference: Wylie & Streeter, "Fluid Transients in Systems" (1993) ch. 3; Chaudhry, "Applied Hydraulic
 * Transients" (2014) ch. 3.
 * @param {object} o
 * @param {number} o.L pipe length, m @param {number} o.D inner diameter, m @param {number} o.a wave speed, m/s
 * @param {number} [o.f=0] Darcy friction factor
 * @param {number|function(number):number} o.H0 upstream reservoir head, m (or function of time, s)
 * @param {number} [o.Q0=0] initial flow rate, m³/s
 * @param {'valve'|'deadEnd'} [o.downstream='valve']
 * @param {number} [o.Hback=0] head downstream of the valve, m
 * @param {number} [o.tClose=0] valve closure time, s (0 = instantaneous)
 * @param {number} [o.m=1] exponent of the closure law
 * @param {function(number):number} [o.tau] custom relative opening τ(t) ∈ [0, 1] (overrides tClose, m)
 * @param {number} [o.Hinit] initial head for a dead-end line at rest, m (default H0 at t = 0)
 * @param {number} [o.n=50] number of reaches @param {number} o.tEnd simulated time, s
 * @param {number} [o.Hvapour] vapour head, m (to flag column separation)
 * @param {number} [o.maxOutput=2000] maximum stored time levels
 * @returns {{ t: number[], Hend: number[], Qend: number[], Hmid: number[], Qin: number[], x: number[],
 *   Hmax: number[], Hmin: number[], HmaxAll: number, HminAll: number, dt: number, cavitation: boolean }}
 *   Hend, Qend = head and flow at the downstream end; Hmid at mid-length; Qin at the reservoir; envelopes
 *   Hmax(x), Hmin(x) in m
 */
export function waterHammerMOC({ L, D, a, f = 0, H0, Q0 = 0, downstream = 'valve', Hback = 0, tClose = 0, m = 1, tau,
  Hinit, n = 50, tEnd, Hvapour, maxOutput = 2000 }) {
  need('L', L, 1e-3, 1e7);
  need('D', D, 1e-4, 50);
  need('a', a, 1, 1e4);
  need('f', f, 0, 1);
  need('Q0', Q0, -1e4, 1e4);
  need('n', n, 2, 1e5);
  need('tEnd', tEnd, 0, 1e7);
  need('tClose', tClose, 0, 1e7);
  if (downstream !== 'valve' && downstream !== 'deadEnd') throw new RangeError(`downstream = ${downstream} is not one of valve, deadEnd`);
  const N = Math.round(n);
  const A = Math.PI * D * D / 4;
  const dx = L / N;
  const dt = dx / a;
  const B = a / (G0 * A);
  const Rf = f * dx / (2 * G0 * D * A * A);
  const Hres = typeof H0 === 'function' ? H0 : () => H0;
  const H0v = need('H0', Hres(0), -1e6, 1e6);
  let H = new Float64Array(N + 1);
  let Q = new Float64Array(N + 1);
  const valve = downstream === 'valve';
  if (valve) {
    for (let i = 0; i <= N; i++) { H[i] = H0v - i * Rf * Q0 * Math.abs(Q0); Q[i] = Q0; }
  } else {
    H.fill(Hinit ?? H0v);
  }
  const dH0 = H[N] - Hback;
  if (valve && Q0 !== 0 && !(dH0 > 0)) throw new RangeError('initial head across the valve must be positive: reduce Q0 or friction');
  const tauOf = tau ?? ((t) => (tClose <= 0 ? (t > 0 ? 0 : 1) : t >= tClose ? 0 : (1 - t / tClose) ** m));
  const Hmax = Float64Array.from(H);
  const Hmin = Float64Array.from(H);
  const steps = Math.ceil(tEnd / dt - 1e-9);
  const stride = Math.max(1, Math.ceil(steps / maxOutput));
  const mid = Math.round(N / 2);
  const out = { t: [0], Hend: [H[N]], Qend: [Q[N]], Hmid: [H[mid]], Qin: [Q[0]] };
  let Hn = new Float64Array(N + 1);
  let Qn = new Float64Array(N + 1);
  for (let s = 1; s <= steps; s++) {
    const t = s * dt;
    for (let i = 1; i < N; i++) {
      const Cp = H[i - 1] + B * Q[i - 1] - Rf * Q[i - 1] * Math.abs(Q[i - 1]);
      const Cm = H[i + 1] - B * Q[i + 1] + Rf * Q[i + 1] * Math.abs(Q[i + 1]);
      Hn[i] = 0.5 * (Cp + Cm);
      Qn[i] = (Cp - Cm) / (2 * B);
    }
    // upstream reservoir
    const Cm0 = H[1] - B * Q[1] + Rf * Q[1] * Math.abs(Q[1]);
    Hn[0] = Hres(t);
    Qn[0] = (Hn[0] - Cm0) / B;
    // downstream boundary
    const CpN = H[N - 1] + B * Q[N - 1] - Rf * Q[N - 1] * Math.abs(Q[N - 1]);
    if (!valve) {
      Qn[N] = 0;
      Hn[N] = CpN;
    } else {
      const tv = Math.min(1, Math.max(0, tauOf(t)));
      if (tv === 0 || Q0 === 0) {
        Qn[N] = 0;
        Hn[N] = CpN;
      } else {
        const Cv = (tv * Q0) ** 2 / (2 * dH0);
        const head = CpN - Hback;
        // forward flow when the available head is positive, reverse flow otherwise
        Qn[N] = head >= 0 ? -B * Cv + Math.sqrt((B * Cv) ** 2 + 2 * Cv * head)
          : B * Cv - Math.sqrt((B * Cv) ** 2 - 2 * Cv * head);
        Hn[N] = CpN - B * Qn[N];
      }
    }
    for (let i = 0; i <= N; i++) {
      if (Hn[i] > Hmax[i]) Hmax[i] = Hn[i];
      if (Hn[i] < Hmin[i]) Hmin[i] = Hn[i];
    }
    [H, Hn] = [Hn, H];
    [Q, Qn] = [Qn, Q];
    if (s % stride === 0 || s === steps) {
      out.t.push(t);
      out.Hend.push(H[N]);
      out.Qend.push(Q[N]);
      out.Hmid.push(H[mid]);
      out.Qin.push(Q[0]);
    }
  }
  let HmaxAll = -Infinity;
  let HminAll = Infinity;
  for (let i = 0; i <= N; i++) { HmaxAll = Math.max(HmaxAll, Hmax[i]); HminAll = Math.min(HminAll, Hmin[i]); }
  return { ...out, x: Array.from({ length: N + 1 }, (_, i) => i * dx), Hmax: Array.from(Hmax), Hmin: Array.from(Hmin),
    HmaxAll, HminAll, dt, cavitation: Hvapour != null && HminAll < Hvapour };
}
