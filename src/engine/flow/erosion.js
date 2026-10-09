// Erosion and solids-transport screening: API RP 14E erosional velocity, DNV-RP-O501 particle-impact erosion for
// straight pipe and bends, Salama (2000), allowable-velocity inversion and minimum transport velocity for sand.
// SI units unless a correlation's native units are stated; erosion rates are returned in mm/year.

import { brent } from '../math/roots.js';

const G0 = 9.80665;
const SEC_PER_YEAR = 3600 * 24 * 365;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * API RP 14E erosional velocity V_e = C/√ρ_m with V_e in ft/s and ρ_m in lb/ft³; converted exactly to SI:
 * V_e [m/s] = C · 0.3048 · √(16.018463/ρ_m[kg/m³]) = 1.2199 C/√ρ_m.
 * Guidance in API RP 14E §2.5.a: C = 100 for continuous and 125 for intermittent service in solids-free fluids;
 * 150–200 (continuous) and up to 250 (intermittent) where corrosion is controlled or corrosion-resistant alloys
 * are used; reduce C when solids are present.
 * Reference: API RP 14E (5th ed. 1991, reaffirmed) §2.5.a. Validity: empirical screening only — it contains no
 * sand rate, particle size or geometry.
 * @param {object} o
 * @param {number} o.rho gas/liquid mixture density at flowing conditions, kg/m³
 * @param {number} [o.C=100] empirical constant in native units, (lb/ft³)^½·ft/s
 * @returns {number} erosional velocity, m/s
 */
export function erosionalVelocity({ rho, C = 100 }) {
  need('rho', rho, 1e-3, 3e4);
  need('C', C, 1, 1000);
  return C * 0.3048 * Math.sqrt(16.01846337396014 / rho);
}

/**
 * Ductile-material impact-angle function of DNV-RP-O501 (2015 edition):
 * F(α) = 0.6 [sin α + 7.2 (sin α − sin²α)]^0.6 [1 − exp(−20 α)], α in radians, 0 ≤ α ≤ π/2 (maximum ≈ 1 near 35°,
 * 0.6 at normal impact).
 * Reference: DNV-RP-O501 "Managing sand production and erosion" (2015), impact-angle function for ductile materials.
 * @param {number} alpha impact angle, rad
 * @returns {number} F, –
 */
export function impactAngleFunction(alpha) {
  need('alpha', alpha, 0, Math.PI / 2 + 1e-12);
  const s = Math.sin(alpha);
  return 0.6 * (s + 7.2 * (s - s * s)) ** 0.6 * (1 - Math.exp(-20 * alpha));
}

/**
 * Generic particle-impact erosion relation of DNV-RP-O501: mass of wall material removed per mass of sand,
 * E = K U_p^n F(α), with the steel-grade constants K = 2.0e-9 (m/s)^−n and n = 2.6 as defaults.
 * Reference: DNV-RP-O501 (2015), material-constants table (steel grades). Other materials: supply K and n from that table.
 * @param {object} o
 * @param {number} o.Up particle impact velocity, m/s
 * @param {number} o.alpha impact angle, rad
 * @param {number} [o.K=2.0e-9] material constant, (m/s)^−n
 * @param {number} [o.n=2.6] velocity exponent
 * @returns {number} erosion ratio, kg wall/kg sand
 */
export function impactErosionRatio({ Up, alpha, K = 2.0e-9, n = 2.6 }) {
  need('Up', Up, 0, 1e3);
  need('K', K, 1e-15, 1e-3);
  need('n', n, 1, 5);
  return K * Up ** n * impactAngleFunction(alpha);
}

/**
 * Sand erosion rate in smooth, straight steel pipe (DNV-RP-O501): E_L = 2.5e-5 U_p^2.6 D^−2 ṁ_p  [mm/year]
 * with U_p in m/s (taken equal to the mixture velocity), D in m, ṁ_p in kg/s.
 * Reference: DNV-RP-O501 (2015), smooth and straight pipes. Validity: steel, straight runs away from
 * disturbances; usually negligible compared with bends.
 * @param {object} o
 * @param {number} o.Up mixture (particle) velocity, m/s @param {number} o.D inner diameter, m
 * @param {number} o.mSand sand mass rate, kg/s
 * @returns {number} erosion rate, mm/year
 */
export function erosionStraightPipe({ Up, D, mSand }) {
  need('Up', Up, 0, 1e3);
  need('D', D, 1e-3, 10);
  need('mSand', mSand, 0, 1e3);
  return 2.5e-5 * Up ** 2.6 * mSand / (D * D);
}

/**
 * Sand erosion rate in a pipe bend (DNV-RP-O501 bend model):
 *   characteristic impact angle α = arctan(1/√(2 R_c)), R_c = bend radius in pipe diameters;
 *   A = ρ_m² tan α U_p D/(ρ_p μ_m);  critical particle size ratio γ_c = ρ_m/(ρ_p [1.88 ln A − 6.04]),
 *   limited to 0.1 when the expression is ≤ 0 or ≥ 0.1;  γ = d_p/D;
 *   particle-size correction G = γ/γ_c for γ < γ_c, otherwise 1;
 *   impacted area A_t = π D²/(4 sin α);
 *   E = K F(α) U_p^n/(ρ_t A_t) · G · C_1 · GF · ṁ_p   [m/s], C_1 = 2.5, GF = geometry factor;
 * converted to mm/year with 3.15e10 (1000 mm/m × 365 d).
 * Reference: DNV-RP-O501 (2015), pipe-bend model. Validity: steel-type ductile materials, bend radius 0.5–50 D,
 * particles 20–5000 µm, U_p up to 200 m/s; stated to be conservative (model uncertainty about a factor 3).
 * Geometry factor GF: 1 for a bend preceded by ≥ 10 D of straight pipe; DNV recommends 2 (bends closer together,
 * same plane) or 3 (out of plane) — supply as argument.
 * @param {object} o
 * @param {number} o.Up particle velocity (mixture velocity), m/s
 * @param {number} o.D inner diameter, m
 * @param {number} o.mSand sand mass rate, kg/s
 * @param {number} o.dp particle diameter, m
 * @param {number} o.rhoM fluid mixture density, kg/m³ @param {number} o.muM fluid mixture viscosity, Pa·s
 * @param {number} [o.rhoP=2650] particle density, kg/m³
 * @param {number} [o.rhoT=7800] target (wall) material density, kg/m³
 * @param {number} [o.Rc=1.5] bend radius of curvature in pipe diameters
 * @param {number} [o.GF=1] geometry factor @param {number} [o.C1=2.5] model factor
 * @param {number} [o.K=2.0e-9] @param {number} [o.n=2.6]
 * @returns {{ rate: number, alpha: number, A: number, gammaC: number, G: number, F: number }} erosion rate
 *   mm/year; impact angle rad; dimensionless group A; critical size ratio; size correction; F(α)
 */
export function erosionBend({ Up, D, mSand, dp, rhoM, muM, rhoP = 2650, rhoT = 7800, Rc = 1.5, GF = 1, C1 = 2.5,
  K = 2.0e-9, n = 2.6 }) {
  need('Up', Up, 0, 1e3);
  need('D', D, 1e-3, 10);
  need('mSand', mSand, 0, 1e3);
  need('dp', dp, 1e-7, 0.05);
  need('rhoM', rhoM, 1e-3, 3e3);
  need('muM', muM, 1e-7, 10);
  need('rhoP', rhoP, 100, 2e4);
  need('rhoT', rhoT, 100, 3e4);
  need('Rc', Rc, 0.5, 50);
  need('GF', GF, 0.1, 10);
  const alpha = Math.atan(1 / Math.sqrt(2 * Rc));
  if (Up === 0) return { rate: 0, alpha, A: 0, gammaC: 0.1, G: 1, F: impactAngleFunction(alpha) };
  const A = rhoM * rhoM * Math.tan(alpha) * Up * D / (rhoP * muM);
  let gammaC = rhoM / (rhoP * (1.88 * Math.log(A) - 6.04));
  if (!(gammaC > 0) || gammaC >= 0.1) gammaC = 0.1;
  const gamma = dp / D;
  const G = gamma < gammaC ? gamma / gammaC : 1;
  const At = Math.PI * D * D / (4 * Math.sin(alpha));
  const F = impactAngleFunction(alpha);
  const rate = K * F * Up ** n / (rhoT * At) * G * C1 * GF * mSand * 1000 * SEC_PER_YEAR;
  return { rate, alpha, A, gammaC, G, F };
}

/**
 * Salama (2000) sand-erosion correlation for multiphase flow:
 * ER = W V_m² d/(S_m D² ρ_m)  [mm/year], with W in kg/day, V_m in m/s, d in µm, D in mm, ρ_m in kg/m³.
 * The geometry constant S_m is an argument; Salama gives 5.5 for pipe bends (and a larger value for blinded
 * tees) — confirm against the paper before relying on a default.
 * Reference: Salama, M.M., "An alternative to API 14E erosional velocity limits for sand-laden fluids", J. Energy
 * Resources Technology 122 (2000) 71–77. Validity: calibrated on gas–liquid tests in elbows; sand ≤ ~400 µm.
 * @param {object} o
 * @param {number} o.mSand sand mass rate, kg/s @param {number} o.Vm mixture velocity, m/s
 * @param {number} o.dp particle diameter, m @param {number} o.D inner diameter, m
 * @param {number} o.rhoM mixture density, kg/m³ @param {number} [o.Sm=5.5] geometry constant
 * @returns {number} erosion rate, mm/year
 */
export function erosionSalama({ mSand, Vm, dp, D, rhoM, Sm = 5.5 }) {
  need('mSand', mSand, 0, 1e3);
  need('Vm', Vm, 0, 1e3);
  need('dp', dp, 1e-7, 0.05);
  need('D', D, 1e-3, 10);
  need('rhoM', rhoM, 1e-3, 3e3);
  need('Sm', Sm, 0.1, 1e4);
  return (mSand * 86400) * Vm * Vm * (dp * 1e6) / (Sm * (D * 1e3) ** 2 * rhoM);
}

/**
 * Maximum mixture velocity for an allowable erosion rate, by inverting a rate model at fixed sand concentration
 * (sand mass rate proportional to velocity) or at fixed sand mass rate.
 * Reference: inversion (Brent) of `erosionBend`, `erosionStraightPipe` or `erosionSalama`.
 * @param {object} o
 * @param {'bend'|'straight'|'salama'} [o.model='bend']
 * @param {number} o.allowable allowable erosion rate, mm/year
 * @param {number} [o.mSand] fixed sand rate, kg/s — or —
 * @param {number} [o.sandConcentration] sand mass per mixture volume, kg/m³ (then ṁ_p = c·V·πD²/4)
 * @param {number} o.D m; remaining model arguments (dp, rhoM, muM, Rc, GF, Sm …) are passed through
 * @param {number} [o.vMaxSearch=300] upper bound of the search, m/s
 * @param {number} [o.interactionFactor=1] erosion–corrosion synergy multiplier applied to the erosion rate (≥ 1;
 *   project-specific, no default correlation is implied)
 * @returns {{ v: number, rate: number }} limiting velocity m/s (vMaxSearch if never reached) and rate there mm/yr
 */
export function maxVelocityForErosion({ model = 'bend', allowable, mSand, sandConcentration, vMaxSearch = 300,
  interactionFactor = 1, ...args }) {
  need('allowable', allowable, 1e-9, 1e6);
  need('interactionFactor', interactionFactor, 1, 100);
  if ((mSand == null) === (sandConcentration == null)) throw new RangeError('give exactly one of mSand and sandConcentration');
  const D = need('D', args.D, 1e-3, 10);
  const rateAt = (v) => {
    const ms = mSand ?? sandConcentration * v * Math.PI * D * D / 4;
    let r;
    if (model === 'bend') r = erosionBend({ ...args, Up: v, mSand: ms }).rate;
    else if (model === 'straight') r = erosionStraightPipe({ Up: v, D, mSand: ms });
    else if (model === 'salama') r = erosionSalama({ ...args, Vm: v, mSand: ms });
    else throw new RangeError(`model = ${model} is not one of bend, straight, salama`);
    return r * interactionFactor;
  };
  if (rateAt(vMaxSearch) <= allowable) return { v: vMaxSearch, rate: rateAt(vMaxSearch) };
  const v = brent((x) => rateAt(x) - allowable, 1e-6, vMaxSearch, { tol: 1e-9 });
  return { v, rate: rateAt(v) };
}

/**
 * Critical (minimum transport) velocity to keep sand suspended and avoid a stationary bed in horizontal pipe,
 * Oroskar & Turian (1980):
 *   V_c/√(g d (s − 1)) = 1.85 C^0.1536 (1 − C)^0.3564 (d/D)^−0.378 Ñ_Re^0.09 x^0.30,
 *   Ñ_Re = D ρ_l √(g d (s − 1))/μ_l,  s = ρ_s/ρ_l,
 * C = solids volume fraction, x = fraction of turbulent eddies with velocity above the settling velocity
 * (≈ 0.96 for most slurries, taken as an argument).
 * Reference: Oroskar & Turian, AIChE J. 26(4) (1980) 550–558. Validity: C 0.01–0.5 (extrapolated with caution to
 * the dilute sand loads of production systems), Newtonian carrier liquid, horizontal flow.
 * @param {object} o
 * @param {number} o.d particle diameter, m @param {number} o.D pipe inner diameter, m
 * @param {number} o.C solids volume fraction, –
 * @param {number} o.rhoL liquid density, kg/m³ @param {number} o.muL liquid viscosity, Pa·s
 * @param {number} [o.rhoS=2650] solids density, kg/m³ @param {number} [o.x=0.96]
 * @param {number} [o.g=9.80665]
 * @returns {number} critical velocity, m/s
 */
export function criticalVelocityOroskarTurian({ d, D, C, rhoL, muL, rhoS = 2650, x = 0.96, g = G0 }) {
  need('d', d, 1e-7, 0.1);
  need('D', D, 1e-3, 10);
  need('C', C, 1e-9, 0.6);
  need('rhoL', rhoL, 1, 3e3);
  need('muL', muL, 1e-6, 10);
  need('rhoS', rhoS, rhoL * (1 + 1e-9), 3e4);
  need('x', x, 1e-3, 1);
  const vs = Math.sqrt(g * d * (rhoS / rhoL - 1));
  const NRe = D * rhoL * vs / muL;
  return 1.85 * vs * C ** 0.1536 * (1 - C) ** 0.3564 * (d / D) ** -0.378 * NRe ** 0.09 * x ** 0.30;
}

/**
 * Durand limit-deposit velocity V_D = F_L √(2 g D (s − 1)).
 * F_L is read from the Durand–Condolios chart as a function of particle size and concentration and is taken as an
 * argument (typically 0.8–1.5; about 1.3–1.4 for uniform sand coarser than ~0.5 mm at a few per cent by volume).
 * Reference: Durand, R., Proc. Minnesota Int. Hydraulics Convention (1953) 89–103.
 * @param {object} o
 * @param {number} o.D m @param {number} o.rhoL kg/m³ @param {number} [o.rhoS=2650] kg/m³
 * @param {number} o.FL Durand factor, – @param {number} [o.g=9.80665]
 * @returns {number} deposit velocity, m/s
 */
export function depositVelocityDurand({ D, rhoL, rhoS = 2650, FL, g = G0 }) {
  need('D', D, 1e-3, 10);
  need('rhoL', rhoL, 1, 3e3);
  need('rhoS', rhoS, rhoL * (1 + 1e-9), 3e4);
  need('FL', FL, 0.1, 3);
  return FL * Math.sqrt(2 * g * D * (rhoS / rhoL - 1));
}
