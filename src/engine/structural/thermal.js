/**
 * Heat transfer and fire engineering for structural members.
 *
 *  - One-dimensional transient conduction through multilayer plane walls and cylindrical shells by an implicit
 *    finite-volume scheme (vertex-centred control volumes, θ-method in time, Picard iteration for radiation and
 *    temperature-dependent properties), with temperature, flux, convection and radiation boundary conditions.
 *  - Steady multilayer resistance network (plane and cylindrical) and lumped-capacitance cooling.
 *  - Nominal fire curves of EN 1991-1-2 §3.2 (standard ISO 834, external, hydrocarbon).
 *  - Carbon steel thermal properties and reduction factors of EN 1993-1-2 §3 (Table 3.1), temperature rise of
 *    unprotected and insulated steel members by the step equations of §4.2.5, critical temperature (§4.2.4),
 *    passive fire protection thickness sizing.
 *  - Thermal strain and stress of restrained members.
 *
 * Units: SI (m, s, W, J, kg). Temperatures are in °C wherever the argument or field name ends in `C`; temperature
 * differences are in K. Heat fluxes are positive into the body.
 */
import { solveTridiagonal } from '../math/linalg.js';
import { brent } from '../math/roots.js';
import { erf } from '../math/special.js';

/** Stefan–Boltzmann constant [W/(m²·K⁴)] (EN 1991-1-2 §3.1(6) uses 5.67e-8). */
export const STEFAN_BOLTZMANN = 5.67e-8;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
const val = (v, x) => (typeof v === 'function' ? v(x) : v);

/* ------------------------------------------------------------------ fire curves */

/**
 * Standard temperature–time curve (ISO 834 / EN 1991-1-2 Eq. 3.4): θg = 20 + 345 log10(8 t + 1), t in minutes.
 * @param {number} t time from ignition [s] (>= 0)
 * @returns {number} gas temperature [°C]
 */
export function isoFireC(t) {
  nonNeg('t', t);
  return 20 + 345 * Math.log10(8 * t / 60 + 1);
}

/**
 * Hydrocarbon fire curve (EN 1991-1-2 Eq. 3.6): θg = 20 + 1080 (1 − 0.325 e^(−0.167 t) − 0.675 e^(−2.5 t)), t in minutes.
 * @param {number} t time [s]
 * @returns {number} gas temperature [°C]
 */
export function hydrocarbonFireC(t) {
  nonNeg('t', t);
  const m = t / 60;
  return 20 + 1080 * (1 - 0.325 * Math.exp(-0.167 * m) - 0.675 * Math.exp(-2.5 * m));
}

/**
 * External fire curve (EN 1991-1-2 Eq. 3.5): θg = 20 + 660 (1 − 0.687 e^(−0.32 t) − 0.313 e^(−3.8 t)), t in minutes.
 * @param {number} t time [s]
 * @returns {number} gas temperature [°C]
 */
export function externalFireC(t) {
  nonNeg('t', t);
  const m = t / 60;
  return 20 + 660 * (1 - 0.687 * Math.exp(-0.32 * m) - 0.313 * Math.exp(-3.8 * m));
}

/**
 * Gas temperature of a named nominal fire curve, with the coefficient of convective heat transfer that
 * EN 1991-1-2 §3.2 attaches to it (25 W/m²K standard and external, 50 W/m²K hydrocarbon).
 * @param {'iso834'|'hydrocarbon'|'external'} curve
 * @param {number} t time [s]
 * @returns {{gasC:number, alphaC:number}} gas temperature [°C], convection coefficient [W/(m²·K)]
 */
export function fireCurve(curve, t) {
  if (curve === 'iso834') return { gasC: isoFireC(t), alphaC: 25 };
  if (curve === 'hydrocarbon') return { gasC: hydrocarbonFireC(t), alphaC: 50 };
  if (curve === 'external') return { gasC: externalFireC(t), alphaC: 25 };
  throw new RangeError(`curve "${curve}" is not one of iso834, hydrocarbon, external`);
}

function resolveFire(fire, alphaC) {
  if (typeof fire === 'function') {
    if (alphaC === undefined) throw new RangeError('alphaC is required when fire is a function');
    return { gas: fire, alphaC: pos('alphaC', alphaC) };
  }
  const a = fireCurve(fire, 0).alphaC;
  return { gas: (t) => fireCurve(fire, t).gasC, alphaC: alphaC === undefined ? a : pos('alphaC', alphaC) };
}

/* ------------------------------------------------------------------ carbon steel at elevated temperature */

/**
 * Specific heat of carbon steel (EN 1993-1-2 §3.4.1.2):
 * 20–600 °C: 425 + 0.773 θ − 1.69e-3 θ² + 2.22e-6 θ³; 600–735 °C: 666 + 13002/(738 − θ);
 * 735–900 °C: 545 + 17820/(θ − 731); 900–1200 °C: 650.
 * @param {number} thetaC steel temperature [°C], 20–1200 (clamped outside)
 * @returns {number} ca [J/(kg·K)]
 */
export function steelSpecificHeat(thetaC) {
  const th = Math.min(1200, Math.max(20, fin('thetaC', thetaC)));
  if (th < 600) return 425 + 0.773 * th - 1.69e-3 * th * th + 2.22e-6 * th ** 3;
  if (th < 735) return 666 + 13002 / (738 - th);
  if (th < 900) return 545 + 17820 / (th - 731);
  return 650;
}

/**
 * Thermal conductivity of carbon steel (EN 1993-1-2 §3.4.1.3): 54 − 3.33e-2 θ below 800 °C, 27.3 above.
 * @param {number} thetaC steel temperature [°C], 20–1200 (clamped outside)
 * @returns {number} λa [W/(m·K)]
 */
export function steelConductivity(thetaC) {
  const th = Math.min(1200, Math.max(20, fin('thetaC', thetaC)));
  return th < 800 ? 54 - 3.33e-2 * th : 27.3;
}

/**
 * Relative thermal elongation of carbon steel from 20 °C (EN 1993-1-2 §3.4.1.1):
 * 20–750 °C: 1.2e-5 θ + 0.4e-8 θ² − 2.416e-4; 750–860 °C: 1.1e-2; 860–1200 °C: 2e-5 θ − 6.2e-3.
 * @param {number} thetaC steel temperature [°C], 20–1200
 * @returns {number} Δl/l [-]
 */
export function steelThermalStrain(thetaC) {
  fin('thetaC', thetaC);
  if (thetaC < 20 || thetaC > 1200) throw new RangeError('thetaC must lie in [20, 1200]');
  if (thetaC < 750) return 1.2e-5 * thetaC + 0.4e-8 * thetaC * thetaC - 2.416e-4;
  if (thetaC <= 860) return 1.1e-2;
  return 2e-5 * thetaC - 6.2e-3;
}

/* EN 1993-1-2 Table 3.1: reduction factors for carbon steel relative to the values at 20 °C. */
const T31 = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200];
const KY = [1, 1, 1, 1, 1, 0.78, 0.47, 0.23, 0.11, 0.06, 0.04, 0.02, 0];
const KP = [1, 1, 0.807, 0.613, 0.42, 0.36, 0.18, 0.075, 0.05, 0.0375, 0.025, 0.0125, 0];
const KE = [1, 1, 0.9, 0.8, 0.7, 0.6, 0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225, 0];

function tab(ys, th) {
  if (th <= T31[0]) return ys[0];
  if (th >= 1200) return 0;
  let i = 1;
  while (T31[i] < th) i++;
  const f = (th - T31[i - 1]) / (T31[i] - T31[i - 1]);
  return ys[i - 1] + f * (ys[i] - ys[i - 1]);
}

/**
 * Reduction factors for the stress–strain relationship of carbon steel at elevated temperature, EN 1993-1-2
 * Table 3.1 with linear interpolation: ky,θ (effective yield strength), kp,θ (proportional limit), kE,θ (slope of
 * the linear elastic range).
 * @param {number} thetaC steel temperature [°C] (values below 20 °C return 1; at and above 1200 °C return 0)
 * @returns {{ky:number, kp:number, kE:number}} reduction factors [-]
 */
export function steelReductionFactors(thetaC) {
  fin('thetaC', thetaC);
  return { ky: tab(KY, thetaC), kp: tab(KP, thetaC), kE: tab(KE, thetaC) };
}

/**
 * Critical temperature of a carbon steel member for a degree of utilisation μ0 at time t = 0.
 * method 'formula': EN 1993-1-2 Eq. (4.22) θa,cr = 39.19 ln[1/(0.9674 μ0^3.833) − 1] + 482 (μ0 ≥ 0.013; members
 * without instability phenomena and with uniform temperature). method 'table': the temperature at which
 * ky,θ of Table 3.1 equals μ0 (exact inverse of the tabulated strength reduction).
 * @param {number} mu0 degree of utilisation Efi,d/Rfi,d,0 [-], 0.013 ≤ μ0 ≤ 1
 * @param {'formula'|'table'} [method='formula']
 * @returns {number} θa,cr [°C]
 */
export function criticalTemperatureC(mu0, method = 'formula') {
  if (!(mu0 >= 0.013 && mu0 <= 1)) throw new RangeError('mu0 must lie in [0.013, 1]');
  if (method === 'formula') return 39.19 * Math.log(1 / (0.9674 * mu0 ** 3.833) - 1) + 482;
  if (method !== 'table') throw new RangeError('method must be "formula" or "table"');
  if (mu0 >= 1) return 400;
  for (let i = 4; i < T31.length - 1; i++) {
    if (mu0 <= KY[i] && mu0 >= KY[i + 1]) return T31[i] + (KY[i] - mu0) / (KY[i] - KY[i + 1]) * (T31[i + 1] - T31[i]);
  }
  return 1200;
}

/**
 * Section factor Am/V of a circular hollow section heated on its outer surface only: π D / [π (D − t) t].
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @returns {number} Am/V [1/m]
 */
export function tubeSectionFactor({ D, t }) {
  pos('D', D); pos('t', t);
  if (t > D / 2) throw new RangeError('t must not exceed D/2');
  return D / ((D - t) * t);
}

/**
 * Temperature of an unprotected steel member in a fire by the step equation of EN 1993-1-2 §4.2.5.1:
 * Δθa,t = ksh (Am/V)/(ca ρa) ḣnet Δt, with the net heat flux of EN 1991-1-2 §3.1
 * ḣnet = αc (θg − θa) + Φ εm εf σ [(θg + 273)⁴ − (θa + 273)⁴] and ca(θa) of §3.4.1.2. Δt must not exceed 5 s.
 * @param {object} o
 * @param {number} o.sectionFactor Am/V [1/m] (not less than 10 1/m in the standard)
 * @param {'iso834'|'hydrocarbon'|'external'|function(number):number} o.fire nominal curve name or gas temperature [°C] as a function of time [s]
 * @param {number} o.duration exposure time [s]
 * @param {number} [o.dt=5] time step [s] (≤ 5)
 * @param {number} [o.ksh=1] correction factor for the shadow effect [-]
 * @param {number} [o.alphaC] convection coefficient [W/(m²·K)] (default: that of the named curve; required for a function)
 * @param {number} [o.emissivity=0.7] surface emissivity of the member εm [-] (0.7 carbon steel, EN 1993-1-2 §2.2(2))
 * @param {number} [o.emissivityFire=1] emissivity of the fire εf [-]
 * @param {number} [o.phi=1] configuration factor Φ [-]
 * @param {number} [o.rho=7850] steel density [kg/m³]
 * @param {number} [o.theta0C=20] initial steel temperature [°C]
 * @returns {{t:Float64Array, gasC:Float64Array, steelC:Float64Array, finalC:number, maxC:number}} time [s], temperatures [°C]
 */
export function steelTemperatureUnprotected({
  sectionFactor, fire, duration, dt = 5, ksh = 1, alphaC, emissivity = 0.7, emissivityFire = 1, phi = 1, rho = 7850, theta0C = 20,
}) {
  pos('sectionFactor', sectionFactor); pos('duration', duration); pos('dt', dt); pos('ksh', ksh); pos('rho', rho);
  if (dt > 5) throw new RangeError('dt must not exceed 5 s (EN 1993-1-2 §4.2.5.1(4))');
  const f = resolveFire(fire, alphaC);
  const n = Math.ceil(duration / dt - 1e-9);
  const h = duration / n;
  const t = new Float64Array(n + 1); const gasC = new Float64Array(n + 1); const steelC = new Float64Array(n + 1);
  steelC[0] = theta0C; gasC[0] = f.gas(0);
  const rad = phi * emissivity * emissivityFire * STEFAN_BOLTZMANN;
  let maxC = theta0C;
  for (let i = 1; i <= n; i++) {
    t[i] = i * h;
    const g = f.gas(t[i]);
    const a = steelC[i - 1];
    const hnet = f.alphaC * (g - a) + rad * ((g + 273) ** 4 - (a + 273) ** 4);
    steelC[i] = a + ksh * sectionFactor / (steelSpecificHeat(a) * rho) * hnet * h;
    gasC[i] = g;
    if (steelC[i] > maxC) maxC = steelC[i];
  }
  return { t, gasC, steelC, finalC: steelC[n], maxC };
}

/**
 * Temperature of a steel member insulated by fire protection material, EN 1993-1-2 §4.2.5.2 Eq. (4.27):
 * Δθa,t = λp (Ap/V) (θg − θa) / [dp ca ρa (1 + φ/3)] Δt − (e^(φ/10) − 1) Δθg (not negative while the gas temperature rises),
 * φ = cp ρp dp (Ap/V)/(ca ρa). Δt must not exceed 30 s. Protection properties are taken as constant.
 * @param {object} o
 * @param {number} o.sectionFactor Ap/V of the protected member [1/m]
 * @param {number} o.dp protection thickness [m]
 * @param {number} o.lambdaP thermal conductivity of the protection [W/(m·K)]
 * @param {number} [o.rhoP=0] density of the protection [kg/m³] (0 neglects its heat capacity, φ = 0)
 * @param {number} [o.cP=0] specific heat of the protection [J/(kg·K)]
 * @param {'iso834'|'hydrocarbon'|'external'|function(number):number} o.fire curve name or gas temperature [°C] as a function of time [s]
 * @param {number} o.duration exposure time [s]
 * @param {number} [o.dt=30] time step [s] (≤ 30)
 * @param {number} [o.rho=7850] steel density [kg/m³]
 * @param {number} [o.theta0C=20] initial steel temperature [°C]
 * @returns {{t:Float64Array, gasC:Float64Array, steelC:Float64Array, finalC:number, maxC:number}} time [s], temperatures [°C]
 */
export function steelTemperatureProtected({ sectionFactor, dp, lambdaP, rhoP = 0, cP = 0, fire, duration, dt = 30, rho = 7850, theta0C = 20 }) {
  pos('sectionFactor', sectionFactor); pos('dp', dp); pos('lambdaP', lambdaP); nonNeg('rhoP', rhoP); nonNeg('cP', cP);
  pos('duration', duration); pos('dt', dt); pos('rho', rho);
  if (dt > 30) throw new RangeError('dt must not exceed 30 s (EN 1993-1-2 §4.2.5.2(3))');
  const gas = typeof fire === 'function' ? fire : resolveFire(fire).gas;
  const n = Math.ceil(duration / dt - 1e-9);
  const h = duration / n;
  const t = new Float64Array(n + 1); const gasC = new Float64Array(n + 1); const steelC = new Float64Array(n + 1);
  steelC[0] = theta0C; gasC[0] = gas(0);
  let maxC = theta0C;
  for (let i = 1; i <= n; i++) {
    t[i] = i * h;
    const g = gas(t[i]);
    const a = steelC[i - 1];
    const ca = steelSpecificHeat(a);
    const fi = cP * rhoP * dp * sectionFactor / (ca * rho);
    const dg = g - gasC[i - 1];
    let d = lambdaP * sectionFactor * (g - a) / (dp * ca * rho * (1 + fi / 3)) * h - (Math.exp(fi / 10) - 1) * dg;
    if (d < 0 && dg > 0) d = 0;
    steelC[i] = a + d; gasC[i] = g;
    if (steelC[i] > maxC) maxC = steelC[i];
  }
  return { t, gasC, steelC, finalC: steelC[n], maxC };
}

/**
 * First time at which a temperature history reaches a limit (linear interpolation between steps).
 * @param {{t:ArrayLike<number>, steelC:ArrayLike<number>}} history result of steelTemperatureUnprotected / steelTemperatureProtected
 * @param {number} limitC temperature limit [°C]
 * @returns {number} time [s], or Infinity when the limit is not reached within the history
 */
export function timeToTemperature(history, limitC) {
  fin('limitC', limitC);
  const { t, steelC } = history;
  if (steelC[0] >= limitC) return t[0];
  for (let i = 1; i < t.length; i++) {
    if (steelC[i] >= limitC) return t[i - 1] + (limitC - steelC[i - 1]) / (steelC[i] - steelC[i - 1]) * (t[i] - t[i - 1]);
  }
  return Infinity;
}

/**
 * Passive fire protection thickness that keeps the steel temperature at or below a limit for the required
 * duration, by root finding on the insulated-member step equation of EN 1993-1-2 §4.2.5.2.
 * @param {object} o
 * @param {number} o.sectionFactor Ap/V [1/m]
 * @param {number} o.lambdaP protection conductivity [W/(m·K)]
 * @param {number} [o.rhoP=0] protection density [kg/m³]
 * @param {number} [o.cP=0] protection specific heat [J/(kg·K)]
 * @param {'iso834'|'hydrocarbon'|'external'|function(number):number} o.fire fire curve
 * @param {number} o.duration required fire resistance time [s]
 * @param {number} o.limitC limiting (critical) steel temperature [°C]
 * @param {number} [o.dpMin=1e-4] smallest thickness considered [m]
 * @param {number} [o.dpMax=0.5] largest thickness considered [m]
 * @returns {{dp:number, steelC:number, required:boolean}} thickness [m] and the steel temperature reached with it [°C]; when dpMin
 *   already satisfies the limit, dp = dpMin and `required` is false (check the unprotected member separately)
 */
export function pfpThickness({ sectionFactor, lambdaP, rhoP = 0, cP = 0, fire, duration, limitC, dpMin = 1e-4, dpMax = 0.5 }) {
  pos('dpMin', dpMin); pos('dpMax', dpMax); fin('limitC', limitC);
  if (!(dpMax > dpMin)) throw new RangeError('dpMax must exceed dpMin');
  const run = (dp) => steelTemperatureProtected({ sectionFactor, dp, lambdaP, rhoP, cP, fire, duration, dt: Math.min(30, duration / 20) }).maxC;
  const lo = run(dpMin);
  if (lo <= limitC) return { dp: dpMin, steelC: lo, required: false };
  const hi = run(dpMax);
  if (hi > limitC) throw new RangeError('limitC cannot be met with a protection thickness up to dpMax');
  const dp = brent((d) => run(d) - limitC, dpMin, dpMax, { tol: 1e-9 });
  return { dp, steelC: run(dp), required: true };
}

/**
 * Strength and stiffness of carbon steel at temperature: fy,θ = ky,θ fy, Ea,θ = kE,θ E (EN 1993-1-2 §3.2.1).
 * @param {object} o
 * @param {number} o.fy yield strength at 20 °C [Pa]
 * @param {number} [o.E=2.1e11] Young's modulus at 20 °C [Pa]
 * @param {number} o.thetaC steel temperature [°C]
 * @returns {{fy:number, fp:number, E:number, ky:number, kp:number, kE:number}} reduced properties [Pa] and factors [-]
 */
export function steelAtTemperature({ fy, E = 2.1e11, thetaC }) {
  pos('fy', fy); pos('E', E);
  const k = steelReductionFactors(thetaC);
  return { fy: k.ky * fy, fp: k.kp * fy, E: k.kE * E, ...k };
}

/* ------------------------------------------------------------------ thermal stress */

/**
 * Stress in a uniformly heated prismatic member whose expansion is restrained. With an axial restraint of
 * stiffness ks in series with the member (stiffness EA/L) the restrained fraction is r = ks/(ks + EA/L) and
 * σ = −r E εth (compression for heating), εth = α ΔT or a supplied free thermal strain. The elastic value is
 * capped at the yield strength when `fy` is given.
 * @param {object} o
 * @param {number} o.E Young's modulus at the member temperature [Pa]
 * @param {number} [o.alpha=1.2e-5] coefficient of thermal expansion [1/K] (ignored when `thermalStrain` is given)
 * @param {number} [o.dT] temperature rise [K]
 * @param {number} [o.thermalStrain] free thermal strain [-] (e.g. steelThermalStrain(θ) − steelThermalStrain(θ0))
 * @param {number} [o.restraint=1] degree of restraint r [-] in [0, 1]; or give A, L and ks
 * @param {number} [o.A] cross-section area [m²]
 * @param {number} [o.L] member length [m]
 * @param {number} [o.ks] axial stiffness of the restraint [N/m]
 * @param {number} [o.fy] yield strength at the member temperature [Pa]
 * @returns {{stress:number, elasticStress:number, restraint:number, freeStrain:number, force:(number|undefined), yielded:boolean}}
 *   stress [Pa] (negative = compression), free strain [-], axial force [N] when A is given
 */
export function restrainedThermalStress({ E, alpha = 1.2e-5, dT, thermalStrain, restraint = 1, A, L, ks, fy }) {
  pos('E', E);
  let eth;
  if (thermalStrain !== undefined) eth = fin('thermalStrain', thermalStrain);
  else eth = alpha * fin('dT', dT);
  let r = restraint;
  if (ks !== undefined) {
    nonNeg('ks', ks); pos('A', A); pos('L', L);
    r = ks === Infinity ? 1 : ks / (ks + E * A / L);
  }
  if (!(r >= 0 && r <= 1)) throw new RangeError('restraint must lie in [0, 1]');
  const elastic = -r * E * eth;
  let stress = elastic; let yielded = false;
  if (fy !== undefined && Math.abs(elastic) > fy) { stress = Math.sign(elastic) * fy; yielded = true; }
  return { stress, elasticStress: elastic, restraint: r, freeStrain: eth, force: A === undefined ? undefined : stress * A, yielded };
}

/* ------------------------------------------------------------------ steady conduction and lumped cooling */

/**
 * Steady one-dimensional heat flow through layers in series with surface film coefficients (thermal resistance
 * network; Incropera et al., Fundamentals of Heat and Mass Transfer, ch. 3).
 * Plane wall, per unit area: R = 1/h1 + Σ t/k + 1/h2. Cylinder, per unit length: R' = 1/(2π r1 h1) + Σ ln(ro/ri)/(2π k) + 1/(2π rn h2).
 * @param {object} o
 * @param {'plane'|'cylinder'} [o.geometry='plane']
 * @param {Array<{thickness:number, k:number}>} o.layers from the inner (left) side outward: thickness [m], conductivity [W/(m·K)]
 * @param {number} [o.rInner] inner radius [m] (cylinder)
 * @param {number} [o.hInner=Infinity] film coefficient on the inner / left surface [W/(m²·K)]
 * @param {number} [o.hOuter=Infinity] film coefficient on the outer / right surface [W/(m²·K)]
 * @param {number} o.TinnerC fluid (or surface when h = Infinity) temperature on the inner / left side [°C]
 * @param {number} o.TouterC fluid temperature on the outer / right side [°C]
 * @returns {{q:number, resistance:number, U:number, interfaceC:number[]}} heat flow inner → outer [W/m² plane, W/m cylinder], total
 *   resistance [m²K/W or mK/W], overall coefficient U [W/(m²·K)] (cylinder: referred to the outer surface), temperatures at the
 *   inner surface, each layer interface and the outer surface [°C]
 */
export function steadyMultilayer({ geometry = 'plane', layers, rInner, hInner = Infinity, hOuter = Infinity, TinnerC, TouterC }) {
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  fin('TinnerC', TinnerC); fin('TouterC', TouterC);
  if (!(hInner > 0) || !(hOuter > 0)) throw new RangeError('hInner and hOuter must be positive');
  const cyl = geometry === 'cylinder';
  if (!cyl && geometry !== 'plane') throw new RangeError('geometry must be "plane" or "cylinder"');
  if (cyl) pos('rInner', rInner);
  const res = [];
  let r = cyl ? rInner : 0;
  res.push(cyl ? 1 / (2 * Math.PI * r * hInner) : 1 / hInner);
  layers.forEach((l, i) => {
    pos(`layers[${i}].thickness`, l.thickness); pos(`layers[${i}].k`, l.k);
    res.push(cyl ? Math.log((r + l.thickness) / r) / (2 * Math.PI * l.k) : l.thickness / l.k);
    r += l.thickness;
  });
  res.push(cyl ? 1 / (2 * Math.PI * r * hOuter) : 1 / hOuter);
  const R = res.reduce((s, v) => s + v, 0);
  const q = (TinnerC - TouterC) / R;
  const interfaceC = [];
  let T = TinnerC - q * res[0];
  interfaceC.push(T);
  for (let i = 1; i <= layers.length; i++) { T -= q * res[i]; interfaceC.push(T); }
  return { q, resistance: R, U: cyl ? 1 / (R * 2 * Math.PI * r) : 1 / R, interfaceC };
}

/**
 * Lumped-capacitance (Newtonian) cooling or heating: T(t) = T∞ + (T0 − T∞) e^(−t/τ), τ = m cp/(h A)
 * (Incropera ch. 5.1). Valid for Biot number Bi = h Lc/k < 0.1, reported when k and Lc are given.
 * @param {object} o
 * @param {number} o.T0C initial temperature [°C]
 * @param {number} o.TinfC ambient temperature [°C]
 * @param {number} o.h heat transfer coefficient [W/(m²·K)]
 * @param {number} o.area exposed surface area [m²]
 * @param {number} o.mass mass [kg]
 * @param {number} o.cp specific heat [J/(kg·K)]
 * @param {number} [o.t=0] time [s]
 * @param {number} [o.k] conductivity of the body [W/(m·K)] (for the Biot number)
 * @param {number} [o.Lc] characteristic length V/A [m] (for the Biot number)
 * @returns {{TC:number, tau:number, biot:(number|undefined), valid:(boolean|undefined), timeTo:function(number):number}} temperature [°C],
 *   time constant [s], Biot number [-], and timeTo(TC) → time [s] to reach a temperature between T0 and T∞
 */
export function lumpedCooling({ T0C, TinfC, h, area, mass, cp, t = 0, k, Lc }) {
  fin('T0C', T0C); fin('TinfC', TinfC); pos('h', h); pos('area', area); pos('mass', mass); pos('cp', cp); nonNeg('t', t);
  const tau = mass * cp / (h * area);
  const biot = k !== undefined && Lc !== undefined ? h * pos('Lc', Lc) / pos('k', k) : undefined;
  const timeTo = (TC) => {
    const ratio = (TC - TinfC) / (T0C - TinfC);
    if (!(ratio > 0 && ratio <= 1)) throw new RangeError('TC must lie between TinfC (exclusive) and T0C');
    return -tau * Math.log(ratio);
  };
  return { TC: TinfC + (T0C - TinfC) * Math.exp(-t / tau), tau, biot, valid: biot === undefined ? undefined : biot < 0.1, timeTo };
}

/* ------------------------------------------------------------------ 1-D transient conduction */

function checkBc(name, bc) {
  if (!bc) return { type: 'insulated' };
  const types = ['insulated', 'temperature', 'flux', 'convection'];
  if (!types.includes(bc.type)) throw new RangeError(`${name}.type must be one of ${types.join(', ')}`);
  if (bc.type === 'temperature' && bc.TC === undefined) throw new RangeError(`${name}.TC is required`);
  if (bc.type === 'flux' && bc.q === undefined) throw new RangeError(`${name}.q is required`);
  if (bc.type === 'convection') {
    if (bc.TinfC === undefined) throw new RangeError(`${name}.TinfC is required`);
    if (!(bc.h >= 0) && typeof bc.h !== 'function') throw new RangeError(`${name}.h must be non-negative`);
    if (bc.emissivity !== undefined && !(bc.emissivity >= 0 && bc.emissivity <= 1)) throw new RangeError(`${name}.emissivity must lie in [0, 1]`);
  }
  return bc;
}

/**
 * One-dimensional transient heat conduction through a multilayer plane wall (per unit area) or a multilayer
 * cylindrical shell / solid cylinder (per unit length): ρ c ∂T/∂t = (1/r^m) ∂/∂x (r^m k ∂T/∂x) + q̇, m = 0 or 1.
 * Vertex-centred finite volumes; the conductance between neighbouring nodes is k/Δx (plane) or 2πk/ln(r₂/r₁)
 * (cylinder, exact for steady radial flow); θ-method in time (θ = 1 backward Euler, unconditionally stable and
 * monotone; θ = 0.5 Crank–Nicolson, second order). Radiation and temperature-dependent properties are treated by
 * Picard iteration within each step. Reference: Patankar, Numerical Heat Transfer and Fluid Flow (1980) ch. 4.
 * Boundary conditions (left = inner / x = 0, right = outer):
 *   { type: 'insulated' } | { type: 'temperature', TC } | { type: 'flux', q } (W/m² into the body) |
 *   { type: 'convection', h, TinfC, emissivity?, TradC? } with q = h (T∞ − Ts) + ε σ (Trad⁴ − Ts⁴) (absolute temperatures;
 *   Trad defaults to T∞). TC, q, h, TinfC, TradC may be numbers or functions of time [s].
 * @param {object} o
 * @param {'plane'|'cylinder'} [o.geometry='plane']
 * @param {Array<{thickness:number, k:(number|function), rho:(number|function), cp:(number|function), n?:number, qdot?:number}>} o.layers
 *   from the left / inner side: thickness [m], conductivity [W/(m·K)], density [kg/m³], specific heat [J/(kg·K)] (numbers or
 *   functions of temperature [°C]), number of cells n (default: about 80 cells shared in proportion to thickness, at least 4 each),
 *   volumetric heat source [W/m³]
 * @param {number} [o.rInner=0] inner radius [m] (cylinder; 0 = solid cylinder, left boundary is then the symmetry axis)
 * @param {number|function(number):number} o.T0C initial temperature [°C], uniform or a function of the coordinate [m]
 * @param {object} [o.left] boundary condition at the left / inner surface (default insulated)
 * @param {object} [o.right] boundary condition at the right / outer surface (default insulated)
 * @param {number} o.tEnd duration [s]
 * @param {number} o.dt time step [s]
 * @param {number} [o.theta=1] time weighting [-] in [0.5, 1]
 * @param {number} [o.nOut=50] number of stored profiles (besides the initial one)
 * @param {number} [o.tol=1e-6] Picard convergence tolerance on temperature [K]
 * @param {number} [o.maxIter=50] Picard iteration limit per step
 * @returns {{x:Float64Array, t:Float64Array, TC:Float64Array[], finalC:Float64Array, qLeft:Float64Array, qRight:Float64Array,
 *   energy:Float64Array, interfaces:number[]}} node coordinates (x or r) [m], stored times [s], temperature profiles [°C], heat flux into the
 *   body at the left / right surface at the stored times [W/m²], stored heat relative to 0 °C per unit area [J/m²] (plane) or per unit
 *   length [J/m] (cylinder) evaluated with the properties at the current temperature, node indices of the layer interfaces
 */
export function conduction1d({
  geometry = 'plane', layers, rInner = 0, T0C, left, right, tEnd, dt, theta = 1, nOut = 50, tol = 1e-6, maxIter = 50,
}) {
  const cyl = geometry === 'cylinder';
  if (!cyl && geometry !== 'plane') throw new RangeError('geometry must be "plane" or "cylinder"');
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  pos('tEnd', tEnd); pos('dt', dt); nonNeg('rInner', rInner);
  if (!(theta >= 0.5 && theta <= 1)) throw new RangeError('theta must lie in [0.5, 1]');
  if (!Number.isInteger(nOut) || nOut < 1) throw new RangeError('nOut must be a positive integer');
  const bcL = checkBc('left', left); const bcR = checkBc('right', right);
  const axis = cyl && rInner === 0;
  if (axis && bcL.type !== 'insulated') throw new RangeError('left must be insulated for a solid cylinder (symmetry axis)');
  let total = 0;
  layers.forEach((l, i) => {
    pos(`layers[${i}].thickness`, l.thickness); total += l.thickness;
    for (const p of ['k', 'rho', 'cp']) if (typeof l[p] !== 'function') pos(`layers[${i}].${p}`, l[p]);
  });
  // grid
  const segLayer = []; const xs = [cyl ? rInner : 0]; const interfaces = [0];
  layers.forEach((l, j) => {
    const n = l.n !== undefined ? l.n : Math.max(4, Math.round(80 * l.thickness / total));
    if (!Number.isInteger(n) || n < 1) throw new RangeError(`layers[${j}].n must be a positive integer`);
    const x0 = xs[xs.length - 1];
    for (let i = 1; i <= n; i++) { xs.push(x0 + l.thickness * i / n); segLayer.push(j); }
    interfaces.push(xs.length - 1);
  });
  const N = xs.length; const nSeg = N - 1;
  const x = Float64Array.from(xs);
  // geometric factors: conductance = k·gk, capacity of the left / right half of a segment = ρc·vl, ρc·vr
  const gk = new Float64Array(nSeg); const vl = new Float64Array(nSeg); const vr = new Float64Array(nSeg);
  for (let s = 0; s < nSeg; s++) {
    const a = x[s]; const b = x[s + 1]; const m = 0.5 * (a + b);
    if (cyl) {
      gk[s] = a > 0 ? 2 * Math.PI / Math.log(b / a) : 2 * Math.PI * m / (b - a);
      vl[s] = Math.PI * (m * m - a * a); vr[s] = Math.PI * (b * b - m * m);
    } else { gk[s] = 1 / (b - a); vl[s] = m - a; vr[s] = b - m; }
  }
  const aL = cyl ? 2 * Math.PI * x[0] : 1; const aR = cyl ? 2 * Math.PI * x[N - 1] : 1;
  const nonlinear = layers.some((l) => typeof l.k === 'function' || typeof l.rho === 'function' || typeof l.cp === 'function')
    || [bcL, bcR].some((b) => b.type === 'convection' && b.emissivity > 0);
  let T = new Float64Array(N);
  for (let i = 0; i < N; i++) T[i] = fin('T0C', val(T0C, x[i]));
  const G = new Float64Array(nSeg); const C = new Float64Array(N); const Q = new Float64Array(N);
  // assemble conductances, capacities and sources at temperature field Tk
  const assemble = (Tk) => {
    C.fill(0); Q.fill(0);
    for (let s = 0; s < nSeg; s++) {
      const l = layers[segLayer[s]];
      const tm = 0.5 * (Tk[s] + Tk[s + 1]);
      G[s] = val(l.k, tm) * gk[s];
      const rc = val(l.rho, tm) * val(l.cp, tm);
      C[s] += rc * vl[s]; C[s + 1] += rc * vr[s];
      if (l.qdot) { Q[s] += l.qdot * vl[s]; Q[s + 1] += l.qdot * vr[s]; }
    }
  };
  // boundary as q·A = hEff·A·(Tref − Ts) + qFix·A ; returns [hA, hA·Tref + qA] at time t for surface temperature Ts
  const bcTerms = (bc, area, t, Ts) => {
    if (bc.type === 'flux') return [0, val(bc.q, t) * area];
    if (bc.type === 'convection') {
      const Tinf = val(bc.TinfC, t); const h = val(bc.h, t);
      let hA = h * area; let rhs = h * area * Tinf;
      if (bc.emissivity > 0) {
        const Tr = (bc.TradC === undefined ? Tinf : val(bc.TradC, t)) + 273.15; const Tsk = Ts + 273.15;
        const hr = bc.emissivity * STEFAN_BOLTZMANN * (Tr * Tr + Tsk * Tsk) * (Tr + Tsk);
        hA += hr * area; rhs += hr * area * (Tr - 273.15);
      }
      return [hA, rhs];
    }
    return [0, 0];
  };
  const surfaceFlux = (bc, area, t, Tprof, node, nb, seg) => {
    // prescribed temperature: conduction flux into the first cell (lags the true wall flux by the heat stored in the half cell)
    if (bc.type === 'temperature') return area > 0 ? G[seg] * (Tprof[node] - Tprof[nb]) / area : 0;
    const [hA, rhs] = bcTerms(bc, area, t, Tprof[node]);
    return area > 0 ? (rhs - hA * Tprof[node]) / area : 0;
  };
  const energyOf = (Tk) => { let e = 0; for (let i = 0; i < N; i++) e += C[i] * Tk[i]; return e; };
  const nSteps = Math.ceil(tEnd / dt - 1e-9);
  const h = tEnd / nSteps;
  const every = Math.max(1, Math.floor(nSteps / nOut));
  const tOut = [0]; const prof = [Float64Array.from(T)]; const qL = []; const qR = []; const en = [];
  assemble(T);
  if (bcL.type === 'temperature') T[0] = val(bcL.TC, 0);
  if (bcR.type === 'temperature') T[N - 1] = val(bcR.TC, 0);
  prof[0] = Float64Array.from(T);
  qL.push(surfaceFlux(bcL, aL, 0, T, 0, 1, 0)); qR.push(surfaceFlux(bcR, aR, 0, T, N - 1, N - 2, nSeg - 1)); en.push(energyOf(T));
  const lower = new Float64Array(N); const diag = new Float64Array(N); const upper = new Float64Array(N); const rhs = new Float64Array(N);
  const expl = new Float64Array(N);
  for (let step = 1; step <= nSteps; step++) {
    const t0 = (step - 1) * h; const t1 = step * h;
    // explicit part (1 − θ)·L(T⁰) with the old-time coefficients
    if (theta < 1) {
      assemble(T);
      const [hAl, rl] = bcTerms(bcL, aL, t0, T[0]); const [hAr, rr] = bcTerms(bcR, aR, t0, T[N - 1]);
      for (let i = 0; i < N; i++) {
        let f = Q[i];
        if (i > 0) f += G[i - 1] * (T[i - 1] - T[i]);
        if (i < N - 1) f += G[i] * (T[i + 1] - T[i]);
        expl[i] = f;
      }
      expl[0] += rl - hAl * T[0]; expl[N - 1] += rr - hAr * T[N - 1];
    }
    let Tk = Float64Array.from(T);
    let Tn = Tk;
    for (let it = 0; it < (nonlinear ? maxIter : 1); it++) {
      assemble(Tk);
      const [hAl, rl] = bcTerms(bcL, aL, t1, Tk[0]); const [hAr, rr] = bcTerms(bcR, aR, t1, Tk[N - 1]);
      for (let i = 0; i < N; i++) {
        const gl = i > 0 ? G[i - 1] : 0; const gr = i < N - 1 ? G[i] : 0;
        lower[i] = -theta * gl; upper[i] = -theta * gr;
        diag[i] = C[i] / h + theta * (gl + gr);
        rhs[i] = C[i] / h * T[i] + theta * Q[i] + (1 - theta) * expl[i];
      }
      diag[0] += theta * hAl; rhs[0] += theta * rl;
      diag[N - 1] += theta * hAr; rhs[N - 1] += theta * rr;
      if (bcL.type === 'temperature') { diag[0] = 1; upper[0] = 0; rhs[0] = val(bcL.TC, t1); }
      if (bcR.type === 'temperature') { diag[N - 1] = 1; lower[N - 1] = 0; rhs[N - 1] = val(bcR.TC, t1); }
      Tn = solveTridiagonal(lower, diag, upper, rhs);
      if (!nonlinear) break;
      let err = 0;
      for (let i = 0; i < N; i++) err = Math.max(err, Math.abs(Tn[i] - Tk[i]));
      Tk = Float64Array.from(Tn);
      if (err < tol) break;
      if (it === maxIter - 1) throw new Error(`conduction1d: Picard iteration did not converge at t = ${t1} s (reduce dt)`);
    }
    T = Float64Array.from(Tn);
    if (step % every === 0 || step === nSteps) {
      assemble(T);
      tOut.push(t1); prof.push(Float64Array.from(T));
      qL.push(surfaceFlux(bcL, aL, t1, T, 0, 1, 0)); qR.push(surfaceFlux(bcR, aR, t1, T, N - 1, N - 2, nSeg - 1)); en.push(energyOf(T));
    }
  }
  return {
    x, t: Float64Array.from(tOut), TC: prof, finalC: T, qLeft: Float64Array.from(qL), qRight: Float64Array.from(qR),
    energy: Float64Array.from(en), interfaces,
  };
}

/**
 * Temperature in a semi-infinite solid after a step change of its surface temperature:
 * T(x, t) = Ts + (T0 − Ts) erf(x / (2 √(α t))), surface flux q = k (Ts − T0)/√(π α t) (Carslaw & Jaeger, Conduction of
 * Heat in Solids, §2.4). Valid while the thermal penetration depth ≈ 4√(αt) is less than the body thickness.
 * @param {object} o
 * @param {number} o.x depth below the surface [m]
 * @param {number} o.t time [s] (> 0)
 * @param {number} o.alpha thermal diffusivity k/(ρ c) [m²/s]
 * @param {number} o.T0C initial temperature [°C]
 * @param {number} o.TsC surface temperature for t > 0 [°C]
 * @param {number} [o.k] conductivity [W/(m·K)] for the surface flux
 * @returns {{TC:number, flux:(number|undefined)}} temperature [°C], surface heat flux into the solid [W/m²]
 */
export function semiInfiniteStep({ x, t, alpha, T0C, TsC, k }) {
  nonNeg('x', x); pos('t', t); pos('alpha', alpha);
  fin('T0C', T0C); fin('TsC', TsC);
  return {
    TC: TsC + (T0C - TsC) * erf(x / (2 * Math.sqrt(alpha * t))),
    flux: k === undefined ? undefined : k * (TsC - T0C) / Math.sqrt(Math.PI * alpha * t),
  };
}
