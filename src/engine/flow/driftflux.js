/**
 * Homogeneous and drift-flux closures for gas–liquid pipe flow: void fraction / holdup, mixture
 * properties, slip-ratio correlations and two-phase friction multipliers.
 *
 * Conventions used throughout the flow engines:
 *   vsl, vsg  superficial velocities [m/s]        vm = vsl + vsg mixture velocity [m/s]
 *   theta     pipe inclination from horizontal [rad], positive for upward flow
 *   holdup    in-situ liquid volume fraction [-]  alpha = 1 - holdup, in-situ gas (void) fraction [-]
 *   lambdaL   no-slip (input) liquid fraction vsl / vm [-]
 */

/** Standard gravity [m/s²]. */
export const G = 9.80665;

/**
 * Validate a numeric argument.
 * @param {string} name argument name used in the error message
 * @param {number} v value
 * @param {{min?: number, max?: number, exclusiveMin?: boolean}} [lim] inclusive bounds (min exclusive on request)
 * @returns {number} v
 * @throws {RangeError} when v is not finite or outside the bounds
 */
export function checkArg(name, v, { min = -Infinity, max = Infinity, exclusiveMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (exclusiveMin && v <= min)) {
    throw new RangeError(`${name} = ${v} is outside the valid range ${exclusiveMin ? '(' : '['}${min}, ${max}]`);
  }
  return v;
}

/**
 * Darcy (Moody) friction factor: 64/Re laminar (Re < 2000), Colebrook–White turbulent (Re > 4000),
 * linear blend between the two in the transition band.
 * Reference: Colebrook (1939), J. Inst. Civ. Eng. 11, 133–156.
 * @param {number} Re Reynolds number [-], >= 0
 * @param {number} [relRough=0] relative roughness ε/D [-], 0 … 0.05
 * @returns {number} Darcy friction factor [-] (0 when Re = 0)
 */
export function darcyFrictionFactor(Re, relRough = 0) {
  checkArg('Re', Re, { min: 0 });
  checkArg('relRough', relRough, { min: 0, max: 0.1 });
  return darcyFrictionFactorFast(Re, relRough);
}

/**
 * darcyFrictionFactor without argument validation, for inner loops of the transient solver.
 * @param {number} Re Reynolds number [-], >= 0
 * @param {number} relRough relative roughness ε/D [-]
 * @returns {number} Darcy friction factor [-]
 */
export function darcyFrictionFactorFast(Re, relRough) {
  if (Re === 0) return 0;
  if (Re <= 2000) return 64 / Re;
  const R = Re >= 4000 ? Re : 4000;
  // Haaland start value, then Newton iterations on x = 1/√f in x + 2 log10(a + b x) = 0.
  const a = relRough / 3.7;
  const b = 2.51 / R;
  let x = -1.8 * Math.log10((relRough / 3.7) ** 1.11 + 6.9 / R);
  for (let i = 0; i < 4; i++) {
    const q = a + b * x;
    x -= (x + 2 * Math.log10(q)) / (1 + 2 * b / (q * Math.LN10));
  }
  const turb = 1 / (x * x);
  if (Re >= 4000) return turb;
  const w = (Re - 2000) / 2000;
  return (1 - w) * (64 / 2000) + w * turb;
}

/**
 * Slip mixture density ρ = H ρL + (1 − H) ρG.
 * @param {number} holdup liquid volume fraction [-], 0 … 1
 * @param {number} rhoL liquid density [kg/m³]
 * @param {number} rhoG gas density [kg/m³]
 * @returns {number} mixture density [kg/m³]
 */
export function mixtureDensity(holdup, rhoL, rhoG) {
  checkArg('holdup', holdup, { min: 0, max: 1 });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  return holdup * rhoL + (1 - holdup) * rhoG;
}

/**
 * Two-phase mixture viscosity by a named mixing rule.
 *  'linear'          μ = H μL + (1 − H) μG, volume weighted (Dukler et al. 1964 when H is the no-slip fraction)
 *  'mcadams'         1/μ = x/μG + (1 − x)/μL (McAdams et al. 1942)
 *  'cicchitti'       μ = x μG + (1 − x) μL (Cicchitti et al. 1960)
 *  'beattie-whalley' μ = μL (1 − β)(1 + 2.5 β) + μG β, β = homogeneous void fraction (Beattie & Whalley 1982)
 *  'arrhenius'       μ = μL^H · μG^(1 − H) (used by Hagedorn & Brown 1965)
 * @param {object} o
 * @param {'linear'|'mcadams'|'cicchitti'|'beattie-whalley'|'arrhenius'} [o.rule='linear']
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.holdup] liquid volume fraction [-] ('linear', 'arrhenius'; for 'beattie-whalley' pass the no-slip fraction)
 * @param {number} [o.quality] gas mass fraction x [-] ('mcadams', 'cicchitti')
 * @returns {number} mixture viscosity [Pa·s]
 */
export function mixtureViscosity({ rule = 'linear', muL, muG, holdup, quality }) {
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('muG', muG, { min: 0, exclusiveMin: true });
  switch (rule) {
    case 'linear':
      checkArg('holdup', holdup, { min: 0, max: 1 });
      return holdup * muL + (1 - holdup) * muG;
    case 'arrhenius':
      checkArg('holdup', holdup, { min: 0, max: 1 });
      return muL ** holdup * muG ** (1 - holdup);
    case 'beattie-whalley': {
      checkArg('holdup', holdup, { min: 0, max: 1 });
      const beta = 1 - holdup;
      return muL * (1 - beta) * (1 + 2.5 * beta) + muG * beta;
    }
    case 'mcadams':
      checkArg('quality', quality, { min: 0, max: 1 });
      return 1 / (quality / muG + (1 - quality) / muL);
    case 'cicchitti':
      checkArg('quality', quality, { min: 0, max: 1 });
      return quality * muG + (1 - quality) * muL;
    default:
      throw new RangeError(`rule = ${rule} is not a known viscosity mixing rule`);
  }
}

/**
 * Homogeneous (no-slip) model: both phases travel at the mixture velocity.
 * Reference: Wallis (1969), One-dimensional Two-phase Flow, ch. 2.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vsg superficial gas velocity [m/s], >= 0 (vsl + vsg > 0)
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @returns {{vm: number, lambdaL: number, alpha: number, holdup: number, rhoNs: number, muNs: number,
 *   massFlux: number, quality: number}} vm [m/s]; lambdaL = holdup, alpha [-]; rhoNs [kg/m³]; muNs [Pa·s]
 *   (volume-weighted); massFlux [kg/(m²·s)]; quality = gas mass fraction [-]
 */
export function homogeneous({ vsl, vsg, rhoL, rhoG, muL, muG }) {
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vsg', vsg, { min: 0 });
  const vm = vsl + vsg;
  checkArg('vsl + vsg', vm, { min: 0, exclusiveMin: true });
  const lambdaL = vsl / vm;
  const rhoNs = mixtureDensity(lambdaL, rhoL, rhoG);
  const muNs = mixtureViscosity({ rule: 'linear', muL, muG, holdup: lambdaL });
  const massFlux = rhoL * vsl + rhoG * vsg;
  return { vm, lambdaL, holdup: lambdaL, alpha: 1 - lambdaL, rhoNs, muNs, massFlux, quality: rhoG * vsg / massFlux };
}

/**
 * Terminal rise velocity of a bubble swarm in an infinite medium, v∞ = 1.53 [g σ (ρL − ρG)/ρL²]^(1/4).
 * Reference: Harmathy (1960), AIChE J. 6(2), 281–288.
 * @param {object} o
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³], < rhoL
 * @param {number} o.sigma surface tension [N/m]
 * @returns {number} rise velocity [m/s]
 */
export function harmathyVelocity({ rhoL, rhoG, sigma }) {
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, max: rhoL });
  checkArg('sigma', sigma, { min: 0, exclusiveMin: true });
  return 1.53 * (G * sigma * (rhoL - rhoG) / (rhoL * rhoL)) ** 0.25;
}

/**
 * Taylor-bubble rise velocity in stagnant liquid in a vertical pipe, v = 0.35 √(g D (ρL − ρG)/ρL).
 * Reference: Dumitrescu (1943); Davies & Taylor (1950), Proc. R. Soc. A 200, 375–390; Nicklin et al. (1962).
 * Valid for inertia-dominated bubbles (Eötvös number > ~100, low liquid viscosity).
 * @param {object} o
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.rhoL=1] liquid density [kg/m³]
 * @param {number} [o.rhoG=0] gas density [kg/m³]
 * @returns {number} rise velocity [m/s]
 */
export function taylorBubbleVelocity({ D, rhoL = 1, rhoG = 0 }) {
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, max: rhoL });
  return 0.35 * Math.sqrt(G * D * (rhoL - rhoG) / rhoL);
}

/**
 * Bendiksen (1984) distribution coefficient and drift velocity of long bubbles in inclined pipes:
 *   Fr = vm/√(gD) < 3.5 :  C0 = 1.05 + 0.15 sin²θ,  vd = (0.35 sinθ + 0.54 cosθ) √(gD)
 *   Fr ≥ 3.5           :  C0 = 1.20,                vd = 0.35 sinθ √(gD)
 * Reference: Bendiksen (1984), Int. J. Multiphase Flow 10(4), 467–483.
 * Validity: 0 ≤ θ ≤ 90° (air–water, D = 19–50 mm in the experiments). For θ < 0 the same expression is
 * evaluated and the result is flagged `extrapolated` — at small downward angles and low velocity the real
 * bubble can turn and travel upstream, which this closure does not represent.
 * @param {object} o
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad], −π/2 … π/2
 * @param {number} [o.froudeSwitch=3.5] Froude number of the switch [-]
 * @returns {{C0: number, vd: number, Fr: number, extrapolated: boolean}} C0 [-], vd [m/s], Fr [-]
 */
export function bendiksen({ vm, D, theta = 0, froudeSwitch = 3.5 }) {
  checkArg('vm', vm);
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  const sgd = Math.sqrt(G * D);
  const Fr = Math.abs(vm) / sgd;
  const s = Math.sin(theta);
  const c = Math.cos(theta);
  if (Fr < froudeSwitch) return { C0: 1.05 + 0.15 * s * s, vd: (0.35 * s + 0.54 * c) * sgd, Fr, extrapolated: theta < 0 };
  return { C0: 1.2, vd: 0.35 * s * sgd, Fr, extrapolated: theta < 0 };
}

/**
 * Hasan–Kabir drift-flux parameters for upward flow in vertical and deviated wells.
 *   bubbly  C0 = 1.2,  vd = Harmathy v∞
 *   slug    C0 = 1.2,  vd = 0.35 √(gDΔρ/ρL) · √(sinθ) (1 + cosθ)^1.2
 *   churn   C0 = 1.15, vd as slug
 *   annular C0 = 1.0,  vd = 0
 * Reference: Hasan & Kabir (1988), SPE Prod. Eng. 3(2), 263–272 and 3(4), 474–482; Hasan, Kabir & Sayarpour
 * (2010), J. Pet. Sci. Eng. 72, 42–49 (parameter table). Validity: 0 < θ ≤ 90° (upward).
 * @param {object} o
 * @param {'bubbly'|'slug'|'churn'|'annular'} o.regime
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=π/2] inclination from horizontal [rad], 0 < θ ≤ π/2
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.sigma surface tension [N/m]
 * @returns {{C0: number, vd: number}} C0 [-], vd [m/s]
 */
export function hasanKabir({ regime, D, theta = Math.PI / 2, rhoL, rhoG, sigma }) {
  checkArg('theta', theta, { min: 0, max: Math.PI / 2, exclusiveMin: true });
  const incl = Math.sqrt(Math.sin(theta)) * (1 + Math.cos(theta)) ** 1.2;
  switch (regime) {
    case 'bubbly': return { C0: 1.2, vd: harmathyVelocity({ rhoL, rhoG, sigma }) };
    case 'slug': return { C0: 1.2, vd: taylorBubbleVelocity({ D, rhoL, rhoG }) * incl };
    case 'churn': return { C0: 1.15, vd: taylorBubbleVelocity({ D, rhoL, rhoG }) * incl };
    case 'annular': return { C0: 1.0, vd: 0 };
    default: throw new RangeError(`regime = ${regime} is not a Hasan–Kabir regime`);
  }
}

/**
 * Regime-dependent drift-flux closure (distribution coefficient C0 and drift velocity vd).
 *   'bubbly'      C0 = 1.2, vd = 1.53 [gσΔρ/ρL²]^¼ sinθ          (Zuber & Findlay 1965; Harmathy 1960)
 *   'slug'        C0 = 1.2, vd = 0.35 √(gDΔρ/ρL) sinθ             (Nicklin et al. 1962)
 *   'bendiksen'   inclination-dependent C0 and vd with the Fr = 3.5 switch (Bendiksen 1984)
 *   'annular'     C0 = 1.0, vd = 0 (gas core and film travel together to first order; override with o.C0 / o.vd)
 *   'homogeneous' C0 = 1.0, vd = 0
 *   'hasan-kabir-bubbly' | '-slug' | '-churn' | '-annular'  (see hasanKabir)
 * The sinθ factor on the vertical rise velocities projects buoyancy on the pipe axis; it is exact at θ = 90°
 * and an engineering approximation elsewhere — use 'bendiksen' for inclined slug flow.
 * @param {object} o
 * @param {string} o.regime closure name (above)
 * @param {number} [o.vm] mixture velocity [m/s] ('bendiksen' only)
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=π/2] inclination [rad]
 * @param {number} [o.rhoL] liquid density [kg/m³]
 * @param {number} [o.rhoG] gas density [kg/m³]
 * @param {number} [o.sigma] surface tension [N/m] (bubbly closures)
 * @param {number} [o.C0] override of the distribution coefficient [-] (typ. 1.0–1.3)
 * @param {number} [o.vd] override of the drift velocity [m/s]
 * @returns {{C0: number, vd: number}} C0 [-], vd [m/s]
 */
export function driftClosure({ regime, vm, D, theta = Math.PI / 2, rhoL, rhoG, sigma, C0, vd }) {
  let r;
  switch (regime) {
    case 'bubbly': r = { C0: 1.2, vd: harmathyVelocity({ rhoL, rhoG, sigma }) * Math.sin(theta) }; break;
    case 'slug': r = { C0: 1.2, vd: taylorBubbleVelocity({ D, rhoL, rhoG }) * Math.sin(theta) }; break;
    case 'bendiksen': { const b = bendiksen({ vm, D, theta }); r = { C0: b.C0, vd: b.vd }; break; }
    case 'annular':
    case 'homogeneous': r = { C0: 1, vd: 0 }; break;
    default:
      if (typeof regime === 'string' && regime.startsWith('hasan-kabir-')) {
        r = hasanKabir({ regime: regime.slice(12), D, theta, rhoL, rhoG, sigma });
      } else throw new RangeError(`regime = ${regime} is not a known drift-flux closure`);
  }
  if (C0 !== undefined) r.C0 = checkArg('C0', C0, { min: 0.5, max: 2 });
  if (vd !== undefined) r.vd = checkArg('vd', vd);
  return r;
}

/**
 * Zuber–Findlay drift-flux relation vG = C0 vm + vd, hence α = vsg / (C0 vm + vd).
 * Reference: Zuber & Findlay (1965), J. Heat Transfer 87(4), 453–468.
 * The void fraction is clipped to [0, 1]; `clipped` reports when the closure was outside its physical range
 * (e.g. counter-current conditions).
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s]
 * @param {number} o.vsg superficial gas velocity [m/s], >= 0
 * @param {number} o.C0 distribution coefficient [-]
 * @param {number} o.vd drift velocity [m/s]
 * @returns {{alpha: number, holdup: number, vG: number, vL: number, vm: number, slipRatio: number,
 *   clipped: boolean}} alpha, holdup [-]; in-situ phase velocities vG, vL [m/s]; slipRatio = vG/vL [-]
 */
export function zuberFindlay({ vsl, vsg, C0, vd }) {
  checkArg('vsl', vsl);
  checkArg('vsg', vsg, { min: 0 });
  checkArg('C0', C0, { min: 0, exclusiveMin: true });
  checkArg('vd', vd);
  const vm = vsl + vsg;
  if (vsg === 0) return { alpha: 0, holdup: 1, vG: C0 * vm + vd, vL: vsl, vm, slipRatio: vsl !== 0 ? (C0 * vm + vd) / vsl : Infinity, clipped: false };
  const vGraw = C0 * vm + vd;
  let alpha = vGraw > 0 ? vsg / vGraw : 1;
  let clipped = false;
  if (alpha > 1 || alpha < 0) { alpha = Math.min(1, Math.max(0, alpha)); clipped = true; }
  const vG = vsg / alpha;
  const vL = alpha < 1 ? vsl / (1 - alpha) : 0;
  return { alpha, holdup: 1 - alpha, vG, vL, vm, slipRatio: vL !== 0 ? vG / vL : Infinity, clipped };
}

/**
 * Holdup from a named drift-flux closure in one call (closure + Zuber–Findlay).
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s]
 * @param {number} o.vsg superficial gas velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=π/2] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} [o.sigma] surface tension [N/m]
 * @param {string} [o.regime='bendiksen'] closure name (see driftClosure)
 * @param {number} [o.C0] override [-]
 * @param {number} [o.vd] override [m/s]
 * @returns {{alpha: number, holdup: number, vG: number, vL: number, vm: number, slipRatio: number,
 *   clipped: boolean, C0: number, vd: number, rhoM: number}} as zuberFindlay plus C0 [-], vd [m/s] and the
 *   slip mixture density rhoM [kg/m³]
 */
export function driftFluxHoldup({ vsl, vsg, D, theta = Math.PI / 2, rhoL, rhoG, sigma, regime = 'bendiksen', C0, vd }) {
  const c = driftClosure({ regime, vm: vsl + vsg, D, theta, rhoL, rhoG, sigma, C0, vd });
  const z = zuberFindlay({ vsl, vsg, C0: c.C0, vd: c.vd });
  return { ...z, C0: c.C0, vd: c.vd, rhoM: mixtureDensity(z.holdup, rhoL, rhoG) };
}

/**
 * Void fraction from quality and slip ratio S = vG/vL:  α = 1 / (1 + S (1 − x)/x · ρG/ρL).
 * Exact kinematic identity (mass balance of both phases).
 * @param {object} o
 * @param {number} o.quality gas mass fraction x [-], 0 … 1
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} [o.S=1] slip ratio [-], > 0
 * @returns {number} void fraction α [-]
 */
export function voidFromSlip({ quality, rhoL, rhoG, S = 1 }) {
  checkArg('quality', quality, { min: 0, max: 1 });
  checkArg('S', S, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  if (quality === 0) return 0;
  return 1 / (1 + S * (1 - quality) / quality * rhoG / rhoL);
}

/**
 * Chisholm slip ratio S = √(1 − x (1 − ρL/ρG)).
 * Reference: Chisholm (1973), Int. J. Heat Mass Transfer 16, 347–358.
 * @param {object} o
 * @param {number} o.quality gas mass fraction x [-]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @returns {number} slip ratio vG/vL [-]
 */
export function slipRatioChisholm({ quality, rhoL, rhoG }) {
  checkArg('quality', quality, { min: 0, max: 1 });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  return Math.sqrt(1 - quality * (1 - rhoL / rhoG));
}

/**
 * Zivi slip ratio S = (ρL/ρG)^(1/3), from minimum kinetic-energy flux of an ideal annular flow.
 * Reference: Zivi (1964), J. Heat Transfer 86, 247–252.
 * @param {object} o
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @returns {number} slip ratio vG/vL [-]
 */
export function slipRatioZivi({ rhoL, rhoG }) {
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  return Math.cbrt(rhoL / rhoG);
}

/** Frictional gradient of one phase flowing alone in the pipe [Pa/m] and its flow type. */
function aloneGradient(v, rho, mu, D, relRough) {
  const Re = rho * v * D / mu;
  return { Re, turbulent: Re > 2000, grad: darcyFrictionFactor(Re, relRough) * rho * v * v / (2 * D) };
}

/**
 * Lockhart–Martinelli two-phase frictional pressure gradient with Chisholm's C-coefficient fit:
 *   X² = (dP/dx)L / (dP/dx)G,  φL² = 1 + C/X + 1/X²,  φG² = 1 + C X + X²,
 *   C = 20 (turbulent–turbulent), 12 (viscous liquid – turbulent gas), 10 (turbulent liquid – viscous gas),
 *   5 (viscous–viscous); a phase is "turbulent" when its superficial Reynolds number exceeds 2000.
 * Holdup: Wallis fit of the Lockhart–Martinelli curve α = (1 + X^0.8)^(−0.378), and Butterworth's form
 *   (1 − α)/α = 0.28 [(1 − x)/x]^0.64 (ρG/ρL)^0.36 (μL/μG)^0.07.
 * References: Lockhart & Martinelli (1949), Chem. Eng. Prog. 45, 39–48; Chisholm (1967), Int. J. Heat Mass
 * Transfer 10, 1767–1778; Wallis (1969) p. 51; Butterworth (1975), Int. J. Multiphase Flow 1, 845–850.
 * Validity: horizontal adiabatic flow, low pressure, small diameter (≤ ~50 mm) in the original data.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.roughness=0] absolute wall roughness [m]
 * @returns {{X: number, C: number, phiL2: number, phiG2: number, gradL: number, gradG: number,
 *   gradFriction: number, alphaWallis: number, holdupWallis: number, alphaButterworth: number,
 *   holdupButterworth: number}} X, C, multipliers [-]; gradients [Pa/m] (pressure loss per metre); void
 *   fractions and holdups [-]
 */
export function lockhartMartinelli({ vsl, vsg, D, rhoL, rhoG, muL, muG, roughness = 0 }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('roughness', roughness, { min: 0 });
  const l = aloneGradient(vsl, rhoL, muL, D, roughness / D);
  const g = aloneGradient(vsg, rhoG, muG, D, roughness / D);
  const X = Math.sqrt(l.grad / g.grad);
  const C = l.turbulent ? (g.turbulent ? 20 : 10) : (g.turbulent ? 12 : 5);
  const phiL2 = 1 + C / X + 1 / (X * X);
  const phiG2 = 1 + C * X + X * X;
  const alphaWallis = (1 + X ** 0.8) ** -0.378;
  const mL = rhoL * vsl;
  const mG = rhoG * vsg;
  const alphaButterworth = 1 / (1 + 0.28 * (mL / mG) ** 0.64 * (rhoG / rhoL) ** 0.36 * (muL / muG) ** 0.07);
  return {
    X, C, phiL2, phiG2, gradL: l.grad, gradG: g.grad, gradFriction: phiL2 * l.grad,
    alphaWallis, holdupWallis: 1 - alphaWallis, alphaButterworth, holdupButterworth: 1 - alphaButterworth,
  };
}

/**
 * Friedel two-phase friction multiplier (liquid-only basis):
 *   φLO² = E + 3.24 F H / (Fr^0.045 We^0.035)
 *   E = (1 − x)² + x² (ρL fGO)/(ρG fLO),  F = x^0.78 (1 − x)^0.224,
 *   H = (ρL/ρG)^0.91 (μG/μL)^0.19 (1 − μG/μL)^0.7,
 *   Fr = G²/(g D ρh²),  We = G² D/(σ ρh),  1/ρh = x/ρG + (1 − x)/ρL.
 * Reference: Friedel (1979), European Two-phase Flow Group Meeting, Ispra, paper E2 (form for horizontal and
 * vertical upward flow). Validity: μL/μG < 1000; recommended for D > ~4 mm.
 * @param {object} o
 * @param {number} o.massFlux total mass flux G [kg/(m²·s)], > 0
 * @param {number} o.quality gas mass fraction x [-], 0 … 1
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s], < muL
 * @param {number} o.sigma surface tension [N/m]
 * @param {number} [o.roughness=0] absolute wall roughness [m]
 * @returns {{phiLo2: number, gradLo: number, gradFriction: number, E: number, F: number, H: number, Fr: number,
 *   We: number, rhoH: number}} multiplier [-]; liquid-only and two-phase frictional gradients [Pa/m];
 *   intermediate groups [-]; homogeneous density rhoH [kg/m³]
 */
export function friedel({ massFlux, quality, D, rhoL, rhoG, muL, muG, sigma, roughness = 0 }) {
  checkArg('massFlux', massFlux, { min: 0, exclusiveMin: true });
  checkArg('quality', quality, { min: 0, max: 1 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('muG', muG, { min: 0, max: muL });
  checkArg('sigma', sigma, { min: 0, exclusiveMin: true });
  const x = quality;
  const fLo = darcyFrictionFactor(massFlux * D / muL, roughness / D);
  const fGo = darcyFrictionFactor(massFlux * D / muG, roughness / D);
  const rhoH = 1 / (x / rhoG + (1 - x) / rhoL);
  const E = (1 - x) ** 2 + x * x * (rhoL * fGo) / (rhoG * fLo);
  const F = x ** 0.78 * (1 - x) ** 0.224;
  const H = (rhoL / rhoG) ** 0.91 * (muG / muL) ** 0.19 * (1 - muG / muL) ** 0.7;
  const Fr = massFlux * massFlux / (G * D * rhoH * rhoH);
  const We = massFlux * massFlux * D / (sigma * rhoH);
  const phiLo2 = E + 3.24 * F * H / (Fr ** 0.045 * We ** 0.035);
  const gradLo = fLo * massFlux * massFlux / (2 * D * rhoL);
  return { phiLo2, gradLo, gradFriction: phiLo2 * gradLo, E, F, H, Fr, We, rhoH };
}
