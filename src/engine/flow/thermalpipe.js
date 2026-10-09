// Pipeline heat transfer: film coefficients, multilayer overall heat-transfer coefficient, burial, steady
// temperature profile, insulation sizing, shutdown cooldown (lumped and 1-D radial transient conduction), warm-up,
// pipe-in-pipe annulus and direct electrical heating demand.
// SI units: temperatures in K (differences in K), lengths m, h and U in W/(m²·K), k in W/(m·K), heat flow per unit
// length W/m, time s.

import { solveTridiagonal } from '../math/linalg.js';
import { brent } from '../math/roots.js';

const SIGMA_SB = 5.670374419e-8; // Stefan–Boltzmann constant, W/(m²·K⁴)
const G0 = 9.80665;

function need(name, v, lo = -Infinity, hi = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < lo || v > hi) {
    throw new RangeError(`${name} = ${v} is outside the valid range [${lo}, ${hi}]`);
  }
  return v;
}

/**
 * Nusselt number for flow inside a circular pipe.
 *  - 'gnielinski' (default for turbulent flow): Nu = (f/8)(Re − 1000) Pr / [1 + 12.7 √(f/8)(Pr^(2/3) − 1)] with the
 *    Petukhov smooth-pipe friction factor f = (0.790 ln Re − 1.64)^−2 (or `f` if supplied);
 *    valid 3000 < Re < 5e6, 0.5 < Pr < 2000.
 *  - 'dittusBoelter': Nu = 0.023 Re^0.8 Pr^n, n = 0.3 for a cooled fluid, 0.4 for a heated fluid;
 *    valid Re > 1e4, 0.6 < Pr < 160, L/D > 10.
 *  - 'siederTate': laminar, Nu = 1.86 (Re Pr D/L)^(1/3) (μ/μw)^0.14, not less than 3.66.
 *  - 'laminar': fully developed, constant wall temperature, Nu = 3.66.
 *  - 'auto': 3.66 for Re < 2300, Gnielinski above 3000, linear blend in between.
 * References: Gnielinski, Int. Chem. Eng. 16 (1976) 359; Dittus & Boelter (1930); Sieder & Tate, Ind. Eng. Chem.
 * 28 (1936) 1429; Incropera et al., "Fundamentals of Heat and Mass Transfer", ch. 8.
 * @param {object} o
 * @param {number} o.Re Reynolds number, – @param {number} o.Pr Prandtl number, –
 * @param {'auto'|'gnielinski'|'dittusBoelter'|'siederTate'|'laminar'} [o.method='auto']
 * @param {boolean} [o.heating=false] fluid is being heated (Dittus–Boelter exponent)
 * @param {number} [o.f] Darcy friction factor for Gnielinski (rough pipes)
 * @param {number} [o.LoverD] length-to-diameter ratio (Sieder–Tate)
 * @param {number} [o.muRatio=1] bulk-to-wall viscosity ratio μ/μw (Sieder–Tate)
 * @returns {number} Nu = h D/k, –
 */
export function nusseltInternal({ Re, Pr, method = 'auto', heating = false, f, LoverD, muRatio = 1 }) {
  need('Re', Re, 0, 1e9);
  need('Pr', Pr, 1e-3, 1e6);
  const gn = (re) => {
    const ff = f ?? (0.790 * Math.log(re) - 1.64) ** -2;
    return (ff / 8) * (re - 1000) * Pr / (1 + 12.7 * Math.sqrt(ff / 8) * (Pr ** (2 / 3) - 1));
  };
  switch (method) {
    case 'laminar': return 3.66;
    case 'siederTate': {
      need('LoverD', LoverD, 1e-6, 1e9);
      need('muRatio', muRatio, 1e-3, 1e3);
      return Math.max(3.66, 1.86 * Math.cbrt(Re * Pr / LoverD) * muRatio ** 0.14);
    }
    case 'dittusBoelter': return 0.023 * Re ** 0.8 * Pr ** (heating ? 0.4 : 0.3);
    case 'gnielinski':
      if (Re <= 1000) throw new RangeError(`Re = ${Re} is below the range of the Gnielinski correlation`);
      return gn(Re);
    case 'auto':
      if (Re < 2300) return 3.66;
      if (Re >= 3000) return gn(Re);
      return 3.66 + (gn(3000) - 3.66) * (Re - 2300) / 700;
    default: throw new RangeError(`method = ${method} is not a known internal-flow correlation`);
  }
}

/**
 * Internal film coefficient h = Nu k/D from fluid properties and velocity (uses `nusseltInternal`).
 * @param {object} o
 * @param {number} o.v mean velocity, m/s @param {number} o.D inner diameter, m
 * @param {number} o.rho kg/m³ @param {number} o.mu Pa·s @param {number} o.k W/(m·K) @param {number} o.cp J/(kg·K)
 * @param {string} [o.method='auto'] @param {boolean} [o.heating=false] @param {number} [o.L] pipe length, m
 * @returns {{ h: number, Nu: number, Re: number, Pr: number }} h in W/(m²·K)
 */
export function internalFilmCoefficient({ v, D, rho, mu, k, cp, method = 'auto', heating = false, L }) {
  need('D', D, 1e-4, 20);
  need('rho', rho, 1e-3, 3e4);
  need('mu', mu, 1e-7, 1e4);
  need('k', k, 1e-4, 1e3);
  need('cp', cp, 1, 1e5);
  const Re = rho * Math.abs(need('v', v, -1e3, 1e3)) * D / mu;
  const Pr = cp * mu / k;
  const Nu = nusseltInternal({ Re, Pr, method, heating, LoverD: L != null ? L / D : undefined });
  return { h: Nu * k / D, Nu, Re, Pr };
}

/**
 * External forced-convection film coefficient for a cylinder in cross-flow (Churchill & Bernstein 1977):
 * Nu = 0.3 + 0.62 Re^½ Pr^⅓ / [1 + (0.4/Pr)^⅔]^¼ · [1 + (Re/282 000)^⅝]^⅘.
 * Reference: Churchill & Bernstein, J. Heat Transfer 99 (1977) 300–306. Validity: Re·Pr > 0.2; properties at the
 * film temperature. Natural convection is not included (use a floor value for still water).
 * @param {object} o
 * @param {number} o.v current velocity normal to the pipe, m/s @param {number} o.D outer diameter, m
 * @param {number} o.rho kg/m³ @param {number} o.mu Pa·s @param {number} o.k W/(m·K) @param {number} o.cp J/(kg·K)
 * @returns {{ h: number, Nu: number, Re: number, Pr: number }} h in W/(m²·K)
 */
export function externalFilmCoefficient({ v, D, rho, mu, k, cp }) {
  need('D', D, 1e-4, 20);
  need('rho', rho, 1e-3, 3e4);
  need('mu', mu, 1e-7, 1e2);
  need('k', k, 1e-4, 1e3);
  need('cp', cp, 1, 1e5);
  const Re = rho * Math.abs(need('v', v, -100, 100)) * D / mu;
  const Pr = cp * mu / k;
  const Nu = 0.3 + 0.62 * Math.sqrt(Re) * Math.cbrt(Pr) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25
    * (1 + (Re / 282000) ** 0.625) ** 0.8;
  return { h: Nu * k / D, Nu, Re, Pr };
}

/**
 * Equivalent external heat-transfer coefficient (referred to the pipe outer surface) of a buried or partially
 * buried pipe.
 * Fully buried (H ≥ D/2): conduction shape factor of an isothermal cylinder below an isothermal plane,
 *   R' = cosh⁻¹(2H'/D)/(2π k_soil) per unit length, so h = 2 k_soil/[D cosh⁻¹(2H'/D)],
 *   where the finite seabed film coefficient is represented by a fictitious extra cover H' = H + k_soil/h_sea
 *   (exact for h_sea → ∞).
 * Partially buried (−D/2 < H < D/2): perimeter-weighted interpolation between the exposed-pipe coefficient and the
 *   just-buried value, h = φ h_exposed + (1 − φ) h_buried(H = D/2), with exposed fraction φ = arccos(2H/D)/π. This
 *   interpolation is an engineering approximation (no closed form exists); expect ±20 %.
 * Exposed (H ≤ −D/2): h = h_exposed.
 * Reference: Carslaw & Jaeger, "Conduction of Heat in Solids" (shape factor of a buried cylinder); Incropera et al.
 * ch. 4 (conduction shape factors).
 * @param {object} o
 * @param {number} o.D outer diameter of the outermost coating, m
 * @param {number} o.H depth of the pipe centreline below the seabed, m (negative: above the seabed)
 * @param {number} o.kSoil soil thermal conductivity, W/(m·K)
 * @param {number} [o.hSea=Infinity] seabed-to-water film coefficient, W/(m²·K)
 * @param {number} [o.hExposed] film coefficient of the exposed pipe surface, W/(m²·K) (needed when H < D/2)
 * @returns {{ h: number, exposedFraction: number, resistance: number }} h in W/(m²·K) on πD, exposed perimeter
 *   fraction, thermal resistance per unit length in m·K/W
 */
export function burialFilmCoefficient({ D, H, kSoil, hSea = Infinity, hExposed }) {
  need('D', D, 1e-3, 20);
  need('H', H, -1e3, 1e4);
  need('kSoil', kSoil, 1e-2, 20);
  need('hSea', hSea, 1e-3, Infinity);
  const extra = Number.isFinite(hSea) ? kSoil / hSea : 0;
  const buried = (depth) => {
    const arg = 2 * (depth + extra) / D;
    return arg > 1 ? 2 * kSoil / (D * Math.acosh(arg)) : Infinity;
  };
  let h;
  let phi = 0;
  if (H >= D / 2) {
    h = buried(H);
  } else {
    if (hExposed == null) throw new RangeError('hExposed is required for a partially buried or exposed pipe');
    need('hExposed', hExposed, 1e-3, 1e6);
    phi = H <= -D / 2 ? 1 : Math.acos(2 * H / D) / Math.PI;
    const hb = Math.min(buried(D / 2), hExposed);
    h = phi * hExposed + (1 - phi) * hb;
  }
  return { h, exposedFraction: phi, resistance: 1 / (h * Math.PI * D) };
}

/**
 * Thickness of a concentric soil annulus that has the same steady thermal resistance as the burial shape factor,
 * for use as an extra layer in the radial transient solver: ln(r_eq/r_o) = cosh⁻¹(2H/D) → r_eq = r_o exp(cosh⁻¹(2H/D)).
 * Reference: equivalence of steady resistances (the transient response of the true half-space is only approximated).
 * @param {object} o
 * @param {number} o.D outer diameter, m @param {number} o.H centreline depth, m (H > D/2)
 * @returns {number} equivalent soil layer thickness, m
 */
export function equivalentSoilThickness({ D, H }) {
  need('D', D, 1e-3, 20);
  need('H', H, D / 2 * (1 + 1e-9), 1e4);
  return 0.5 * D * (Math.exp(Math.acosh(2 * H / D)) - 1);
}

/**
 * Overall heat-transfer coefficient of a multilayer pipe (steady radial conduction in series with the films):
 *   1/(U_ref D_ref) = 1/(h_i D_i) + Σ_j ln(D_j,out/D_j,in)/(2 k_j) + 1/(h_o D_o).
 * Reference: Incropera et al. ch. 3 (radial systems); industry practice refers U to the steel inner or outer diameter.
 * @param {object} o
 * @param {number} o.Di inner diameter of the innermost layer, m
 * @param {{ t: number, k: number }[]} o.layers layers from inside out: thickness m, conductivity W/(m·K)
 * @param {number} [o.hi=Infinity] internal film coefficient, W/(m²·K)
 * @param {number} [o.ho=Infinity] external film coefficient (on the outermost diameter), W/(m²·K)
 * @param {number} [o.Dref] reference diameter for U, m (default Di)
 * @returns {{ U: number, Dref: number, Do: number, UA: number, resistance: number,
 *   resistances: { internal: number, layers: number[], external: number } }}
 *   U in W/(m²·K) on π·Dref; UA = heat loss per unit length per kelvin, W/(m·K); resistances in m·K/W
 */
export function overallU({ Di, layers = [], hi = Infinity, ho = Infinity, Dref }) {
  need('Di', Di, 1e-4, 20);
  need('hi', hi, 1e-6, Infinity);
  need('ho', ho, 1e-6, Infinity);
  let d = Di;
  const Rl = [];
  layers.forEach((l, j) => {
    need(`layers[${j}].t`, l.t, 0, 20);
    need(`layers[${j}].k`, l.k, 1e-5, 1e4);
    const dOut = d + 2 * l.t;
    Rl.push(Math.log(dOut / d) / (2 * Math.PI * l.k));
    d = dOut;
  });
  const Rin = Number.isFinite(hi) ? 1 / (hi * Math.PI * Di) : 0;
  const Rout = Number.isFinite(ho) ? 1 / (ho * Math.PI * d) : 0;
  const Rtot = Rin + Rl.reduce((s, r) => s + r, 0) + Rout;
  if (!(Rtot > 0)) throw new RangeError('total thermal resistance is zero: give at least one layer or film coefficient');
  const ref = Dref ?? Di;
  need('Dref', ref, 1e-4, 50);
  return { U: 1 / (Rtot * Math.PI * ref), Dref: ref, Do: d, UA: 1 / Rtot, resistance: Rtot,
    resistances: { internal: Rin, layers: Rl, external: Rout } };
}

/**
 * Steady temperature change of a flowing fluid over one pipe segment with constant properties (analytic solution
 * of the energy equation):
 *   dT/dx = −(T − T_amb)/λ + μ_JT dP/dx − (g/cp) dz/dx,  λ = ṁ cp/(U π D)
 *   T(L) = T∞ + (T_in − T∞) exp(−L/λ),  T∞ = T_amb + λ (μ_JT dP/dx − (g/cp) dz/dx).
 * The elevation term is the potential-energy change; kinetic-energy change is neglected. For an incompressible
 * liquid μ_JT = −1/(ρ cp) (frictional heating).
 * Reference: Coulter & Bardon (1979) pipeline temperature equation; Alves, Alhanati & Shoham, SPE Prod. Eng. 7
 * (1992) 363 (general form).
 * @param {object} o
 * @param {number} o.Tin inlet temperature, K @param {number} o.Tamb ambient temperature, K
 * @param {number} o.L segment length, m @param {number} o.mDot mass flow rate, kg/s (> 0)
 * @param {number} o.cp fluid heat capacity, J/(kg·K)
 * @param {number} o.U overall heat-transfer coefficient, W/(m²·K), referred to diameter `D`
 * @param {number} o.D reference diameter of U, m
 * @param {number} [o.muJT=0] Joule–Thomson coefficient, K/Pa
 * @param {number} [o.dPdx=0] pressure gradient along the flow, Pa/m (negative for a pressure drop)
 * @param {number} [o.dzdx=0] elevation gradient, m/m (positive uphill)
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{ Tout: number, Tinf: number, lambda: number, heatLoss: number }} K, K, relaxation length m, heat
 *   transferred to the surroundings over the segment W
 */
export function segmentTemperature({ Tin, Tamb, L, mDot, cp, U, D, muJT = 0, dPdx = 0, dzdx = 0, g = G0 }) {
  need('Tin', Tin, 1, 5000);
  need('Tamb', Tamb, 1, 5000);
  need('L', L, 0, 1e7);
  need('mDot', mDot, 1e-12, 1e7);
  need('cp', cp, 1, 1e5);
  need('U', U, 0, 1e6);
  need('D', D, 1e-4, 50);
  const source = muJT * dPdx - (g / cp) * dzdx; // K/m
  if (U === 0) {
    const Tout = Tin + source * L;
    return { Tout, Tinf: Infinity, lambda: Infinity, heatLoss: 0 };
  }
  const lambda = mDot * cp / (U * Math.PI * D);
  const Tinf = Tamb + lambda * source;
  const e = Math.exp(-L / lambda);
  const Tout = Tinf + (Tin - Tinf) * e;
  // ∫ U π D (T − T_amb) dx
  const heatLoss = U * Math.PI * D * ((Tinf - Tamb) * L + (Tin - Tinf) * lambda * (1 - e));
  return { Tout, Tinf, lambda, heatLoss };
}

/**
 * Steady temperature profile along a pipeline route by marching `segmentTemperature` over the profile points.
 * Per-point arrays may be given for ambient temperature, U and pressure (for the Joule–Thomson term).
 * Reference: as `segmentTemperature`. Validity: constant cp and μ_JT per call (re-evaluate outside for strongly
 * varying properties, as the non-isothermal hydraulic solver does).
 * @param {object} o
 * @param {number} o.Tin inlet temperature, K
 * @param {number} o.mDot kg/s @param {number} o.cp J/(kg·K)
 * @param {number|number[]} o.U W/(m²·K) on πD; scalar or one value per segment
 * @param {number} o.D reference diameter of U, m
 * @param {{ x: number, z?: number }[]} o.profile route points: distance m, elevation m
 * @param {number|number[]} o.Tamb ambient temperature, K; scalar or one value per segment
 * @param {number} [o.muJT=0] K/Pa
 * @param {number[]} [o.P] pressure at each profile point, Pa
 * @returns {{ x: number[], T: number[], Tarrival: number, heatLoss: number }} K; total heat loss W
 */
export function temperatureProfile({ Tin, mDot, cp, U, D, profile, Tamb, muJT = 0, P }) {
  if (!Array.isArray(profile) || profile.length < 2) throw new RangeError('profile needs at least two points');
  if (P && P.length !== profile.length) throw new RangeError('P must have one value per profile point');
  const x = [profile[0].x];
  const T = [Tin];
  let heat = 0;
  for (let i = 1; i < profile.length; i++) {
    const L = profile[i].x - profile[i - 1].x;
    if (!(L > 0)) throw new RangeError(`profile x must increase strictly (point ${i})`);
    const seg = segmentTemperature({
      Tin: T[i - 1], L, mDot, cp, D, muJT,
      Tamb: Array.isArray(Tamb) ? Tamb[i - 1] : Tamb,
      U: Array.isArray(U) ? U[i - 1] : U,
      dPdx: P ? (P[i] - P[i - 1]) / L : 0,
      dzdx: ((profile[i].z ?? 0) - (profile[i - 1].z ?? 0)) / L,
    });
    x.push(profile[i].x);
    T.push(seg.Tout);
    heat += seg.heatLoss;
  }
  return { x, T, Tarrival: T[T.length - 1], heatLoss: heat };
}

/**
 * Insulation thickness of one layer needed to reach a target overall U-value or a target arrival temperature.
 * Solved with Brent's method on the layer thickness between 0 and `tMax`.
 * Reference: `overallU` and `segmentTemperature` (horizontal line, no Joule–Thomson term for the arrival target).
 * @param {object} o
 * @param {number} o.Di m
 * @param {{ t: number, k: number }[]} o.layers all layers; the thickness of layer `index` is the unknown
 * @param {number} o.index index of the insulation layer in `layers`
 * @param {number} [o.hi=Infinity] @param {number} [o.ho=Infinity] W/(m²·K)
 * @param {number} [o.Dref] reference diameter of U, m (default Di)
 * @param {number} [o.Utarget] target U, W/(m²·K) — give this or the arrival-temperature set below
 * @param {number} [o.Tarrival] target arrival temperature, K
 * @param {number} [o.Tin] K @param {number} [o.Tamb] K @param {number} [o.L] m @param {number} [o.mDot] kg/s
 * @param {number} [o.cp] J/(kg·K)
 * @param {number} [o.tMax=1] upper bound of the search, m
 * @returns {{ t: number, U: number, Do: number }} required thickness m (0 if the bare pipe already meets the
 *   target), resulting U and outer diameter
 */
export function insulationThickness({ Di, layers, index, hi = Infinity, ho = Infinity, Dref, Utarget, Tarrival,
  Tin, Tamb, L, mDot, cp, tMax = 1 }) {
  if (!Array.isArray(layers) || !(index >= 0 && index < layers.length)) throw new RangeError('index must address a layer');
  const uOf = (t) => overallU({ Di, hi, ho, Dref, layers: layers.map((l, j) => (j === index ? { ...l, t } : l)) });
  let target = Utarget;
  if (target == null) {
    need('Tarrival', Tarrival, 1, 5000);
    need('Tin', Tin, 1, 5000);
    need('Tamb', Tamb, 1, 5000);
    if (!((Tin - Tamb) * (Tarrival - Tamb) > 0) || Math.abs(Tarrival - Tamb) >= Math.abs(Tin - Tamb)) {
      throw new RangeError('Tarrival must lie strictly between Tamb and Tin');
    }
    // T_arr = T_amb + (T_in − T_amb) exp(−U π D L/(ṁ cp))
    target = -Math.log((Tarrival - Tamb) / (Tin - Tamb)) * need('mDot', mDot, 1e-9, 1e7) * need('cp', cp, 1, 1e5)
      / (Math.PI * (Dref ?? Di) * need('L', L, 1e-6, 1e7));
  }
  need('Utarget', target, 1e-9, 1e6);
  const f = (t) => uOf(t).U - target;
  if (f(0) <= 0) { const u = uOf(0); return { t: 0, U: u.U, Do: u.Do }; }
  if (f(tMax) > 0) throw new RangeError(`target cannot be met with a thickness up to tMax = ${tMax} m`);
  const t = brent(f, 0, tMax, { tol: 1e-10 });
  const u = uOf(t);
  return { t, U: u.U, Do: u.Do };
}

/**
 * Single-capacity (lumped) shutdown cooldown: all heat capacities at the contents temperature,
 *   T(t) = T_amb + (T_0 − T_amb) exp(−t/τ),  τ = Σ(m' c)/(U π D),
 * and the time to cool to a threshold, t = τ ln[(T_0 − T_amb)/(T_th − T_amb)].
 * Reference: Newton's law of cooling (Biot number ≪ 1); a conservative-fast estimate when insulation heat capacity
 * is excluded from Σ m'c.
 * @param {object} o
 * @param {number} o.T0 initial temperature, K @param {number} o.Tamb K
 * @param {number} o.U W/(m²·K) on πD @param {number} o.D reference diameter of U, m
 * @param {number} o.heatCapacity Σ(m' c): heat capacity per unit length of contents + wall, J/(m·K)
 * @param {number} [o.Tthreshold] temperature of interest (hydrate or wax appearance), K
 * @param {number[]} [o.t] times at which to evaluate, s
 * @returns {{ tau: number, timeTo: number|null, T: number[]|null }} time constant s; time to threshold s (Infinity
 *   if never reached, 0 if already below); temperatures at `t`
 */
export function cooldownLumped({ T0, Tamb, U, D, heatCapacity, Tthreshold, t }) {
  need('T0', T0, 1, 5000);
  need('Tamb', Tamb, 1, 5000);
  need('U', U, 1e-9, 1e6);
  need('D', D, 1e-4, 50);
  need('heatCapacity', heatCapacity, 1e-6, 1e12);
  const tau = heatCapacity / (U * Math.PI * D);
  let timeTo = null;
  if (Tthreshold != null) {
    const r = (Tthreshold - Tamb) / (T0 - Tamb);
    timeTo = r >= 1 ? 0 : r <= 0 ? Infinity : -tau * Math.log(r);
  }
  return { tau, timeTo, T: t ? Array.from(t, (ti) => Tamb + (T0 - Tamb) * Math.exp(-ti / tau)) : null };
}

function buildRadialNetwork({ Di, contents, layers, hi, ho, cellsPerLayer }) {
  // Finite-volume cells in each layer (equal radial thickness); node at the cell mid-radius.
  const cells = [];
  let r = Di / 2;
  layers.forEach((l, j) => {
    need(`layers[${j}].t`, l.t, 1e-9, 50);
    need(`layers[${j}].k`, l.k, 1e-5, 1e4);
    need(`layers[${j}].rho`, l.rho, 1e-3, 3e4);
    need(`layers[${j}].cp`, l.cp, 1, 1e5);
    const nc = Math.max(1, Math.round(l.n ?? cellsPerLayer));
    const dr = l.t / nc;
    for (let c = 0; c < nc; c++) {
      const rin = r + c * dr;
      const rout = rin + dr;
      cells.push({ rin, rout, rm: 0.5 * (rin + rout), k: l.k, C: l.rho * l.cp * Math.PI * (rout * rout - rin * rin), layer: j });
    }
    r += l.t;
  });
  const ro = r;
  const n = cells.length;
  // conductances per unit length between consecutive nodes: [contents–cell0, cell0–cell1, …, cell(n−1)–ambient]
  const G = new Float64Array(n + 1);
  const half = (ra, rb, k) => Math.log(rb / ra) / (2 * Math.PI * k);
  const Rin = (Number.isFinite(hi) ? 1 / (hi * Math.PI * Di) : 0) + half(cells[0].rin, cells[0].rm, cells[0].k);
  G[0] = 1 / Rin;
  for (let i = 1; i < n; i++) {
    G[i] = 1 / (half(cells[i - 1].rm, cells[i - 1].rout, cells[i - 1].k) + half(cells[i].rin, cells[i].rm, cells[i].k));
  }
  G[n] = 1 / (half(cells[n - 1].rm, cells[n - 1].rout, cells[n - 1].k) + (Number.isFinite(ho) ? 1 / (ho * Math.PI * 2 * ro) : 0));
  const Cc = contents ? contents.rho * contents.cp * Math.PI * Di * Di / 4 : 0;
  return { cells, G, n, ro, Cc };
}

function steadyWallProfile(net, Tfluid, Tamb) {
  // series network: temperature drops in proportion to resistance
  const { G, n } = net;
  let Rtot = 0;
  for (let i = 0; i <= n; i++) Rtot += 1 / G[i];
  const q = (Tfluid - Tamb) / Rtot;
  const T = new Float64Array(n);
  let Tprev = Tfluid;
  for (let i = 0; i < n; i++) { Tprev -= q / G[i]; T[i] = Tprev; }
  return { T, q, Rtot };
}

/**
 * Transient radial conduction through the pipe wall, coatings and (optionally) an equivalent soil annulus, with
 * the contents as a well-mixed lump — shutdown cooldown — or held at a fixed fluid temperature — warm-up.
 *
 * Discretisation: finite volumes in radius (`cellsPerLayer` cells of equal thickness per layer, or `n` per layer),
 * cylindrical half-cell resistances ln(r_b/r_a)/(2πk) between nodes, film resistances at both ends, θ-method in
 * time (θ = 1 backward Euler, unconditionally stable and first-order; θ = 0.5 Crank–Nicolson, second-order).
 * With `cellsPerLayer = 1` this is the classical lumped multi-capacity model (one node per layer).
 * The scheme conserves energy exactly: the change of stored heat equals the time-integrated ambient loss.
 * Initial condition: `initial = 'steady'` — steady conduction profile between `T0` (fluid) and `Tamb`;
 * `'uniform'` — everything at `T0`; `'ambient'` — wall at `Tamb`, contents at `T0`; or an array of cell temperatures.
 * Reference: Carslaw & Jaeger ch. 7; Patankar, "Numerical Heat Transfer and Fluid Flow" ch. 4 (finite-volume
 * discretisation). Validity: axisymmetric conduction, constant properties, no phase change or natural-convection
 * model in the contents (represent it through `hi`); buried pipes are approximated by an annulus
 * (`equivalentSoilThickness`), which reproduces the steady loss but only approximates the soil transient.
 *
 * @param {object} o
 * @param {number} o.Di inner diameter, m
 * @param {{ rho: number, cp: number }} [o.contents] contents density kg/m³ and heat capacity J/(kg·K) (effective
 *   values for a multiphase hold-up); required for mode 'cooldown'
 * @param {{ t: number, k: number, rho: number, cp: number, n?: number }[]} o.layers from inside out
 * @param {number} [o.hi=Infinity] contents-to-wall film coefficient, W/(m²·K)
 * @param {number} [o.ho=Infinity] outer film coefficient, W/(m²·K)
 * @param {number} o.Tamb ambient temperature, K
 * @param {number} o.T0 initial contents temperature (cooldown) or fixed fluid temperature (warm-up), K
 * @param {'cooldown'|'warmup'} [o.mode='cooldown']
 * @param {'steady'|'uniform'|'ambient'|number[]} [o.initial] default 'steady' for cooldown, 'ambient' for warm-up
 * @param {number} o.tEnd end time, s @param {number} o.dt time step, s
 * @param {number} [o.cellsPerLayer=8] @param {number} [o.theta=1]
 * @param {number} [o.maxOutput=600] maximum number of stored time levels
 * @returns {{ t: number[], Tcontents: number[], TwallInner: number[], Touter: number[], qLoss: number[],
 *   qFluid: number[], r: number[], Tfinal: number[], energyStored0: number, energyStored: number,
 *   energyLost: number, energyIn: number, Usteady: number }}
 *   temperatures K; qLoss = heat flow to ambient W/m; qFluid = heat flow from the contents/fluid into the wall W/m;
 *   r = cell-node radii m; energies J/m relative to T_amb (energyIn = heat supplied by the fluid in warm-up mode);
 *   Usteady = steady U on πDi, W/(m²·K)
 */
export function radialTransient({ Di, contents, layers, hi = Infinity, ho = Infinity, Tamb, T0, mode = 'cooldown',
  initial, tEnd, dt, cellsPerLayer = 8, theta = 1, maxOutput = 600 }) {
  need('Di', Di, 1e-4, 20);
  need('Tamb', Tamb, 1, 5000);
  need('T0', T0, 1, 5000);
  need('tEnd', tEnd, 0, 1e10);
  need('dt', dt, 1e-9, 1e9);
  need('theta', theta, 0.5, 1);
  need('hi', hi, 1e-9, Infinity);
  need('ho', ho, 1e-9, Infinity);
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must contain at least one layer');
  if (mode !== 'cooldown' && mode !== 'warmup') throw new RangeError(`mode = ${mode} is not one of cooldown, warmup`);
  const cool = mode === 'cooldown';
  if (cool) {
    if (!contents) throw new RangeError('contents (rho, cp) are required for a cooldown');
    need('contents.rho', contents.rho, 1e-3, 3e4);
    need('contents.cp', contents.cp, 1, 1e5);
  }
  const net = buildRadialNetwork({ Di, contents: cool ? contents : null, layers, hi, ho, cellsPerLayer });
  const { cells, G, n, Cc } = net;
  const steady = steadyWallProfile(net, T0, Tamb);
  const init = initial ?? (cool ? 'steady' : 'ambient');
  let Tw;
  if (Array.isArray(init) || ArrayBuffer.isView(init)) {
    if (init.length !== n) throw new RangeError(`initial must have ${n} cell temperatures`);
    Tw = Float64Array.from(init);
  } else if (init === 'steady') Tw = steady.T.slice();
  else if (init === 'uniform') Tw = new Float64Array(n).fill(T0);
  else if (init === 'ambient') Tw = new Float64Array(n).fill(Tamb);
  else throw new RangeError(`initial = ${init} is not one of steady, uniform, ambient or an array`);

  // unknown vector: [contents (cooldown only), cell 0 … cell n−1]
  const off = cool ? 1 : 0;
  const N = n + off;
  let T = new Float64Array(N);
  if (cool) T[0] = T0;
  for (let i = 0; i < n; i++) T[i + off] = Tw[i];
  const C = new Float64Array(N);
  if (cool) C[0] = Cc;
  for (let i = 0; i < n; i++) C[i + off] = cells[i].C;
  // link k joins unknown k−1 and k (k = 1…N−1); gL[0] is the fixed-temperature link on the left (warm-up only),
  // gR is the ambient link on the right
  const gLink = new Float64Array(N + 1);
  for (let i = 0; i <= n; i++) gLink[i + off] = G[i];
  // in cooldown gLink[0] is unused (adiabatic centre); in warm-up gLink[0] = G[0] couples to the fixed fluid T
  const Tleft = T0;
  const lower = new Float64Array(N);
  const diag = new Float64Array(N);
  const upper = new Float64Array(N);
  const rhs = new Float64Array(N);
  const energy = (v) => { let e = 0; for (let i = 0; i < N; i++) e += C[i] * (v[i] - Tamb); return e; };
  const qOut = (v) => gLink[N] * (v[N - 1] - Tamb);
  const qIn = (v) => (cool ? G[0] * (v[0] - v[1]) : G[0] * (Tleft - v[0]));
  const innerSurface = (v) => {
    // wall temperature at the inner surface: from the heat flow through the film
    const q = qIn(v);
    const Tf = cool ? v[0] : Tleft;
    return Number.isFinite(hi) ? Tf - q / (hi * Math.PI * Di) : Tf;
  };
  const outerSurface = (v) => (Number.isFinite(ho) ? Tamb + qOut(v) / (ho * Math.PI * 2 * net.ro) : Tamb);

  const nSteps = Math.max(1, Math.ceil(tEnd / dt - 1e-9));
  const h = tEnd > 0 ? tEnd / nSteps : dt;
  const stride = Math.max(1, Math.ceil(nSteps / maxOutput));
  const out = { t: [0], Tcontents: [cool ? T[0] : Tleft], TwallInner: [innerSurface(T)], Touter: [outerSurface(T)],
    qLoss: [qOut(T)], qFluid: [qIn(T)] };
  const E0 = energy(T);
  let lost = 0;
  let supplied = 0;
  for (let s = 1; s <= (tEnd > 0 ? nSteps : 0); s++) {
    for (let i = 0; i < N; i++) {
      const gl = i === 0 ? (cool ? 0 : gLink[0]) : gLink[i];
      const gr = gLink[i + 1];
      const Tl = i === 0 ? Tleft : T[i - 1];
      const Tr = i === N - 1 ? Tamb : T[i + 1];
      const cdt = C[i] / h;
      lower[i] = i === 0 ? 0 : -theta * gl;
      upper[i] = i === N - 1 ? 0 : -theta * gr;
      diag[i] = cdt + theta * (gl + gr);
      rhs[i] = cdt * T[i] + (1 - theta) * (gl * (Tl - T[i]) + gr * (Tr - T[i]));
      if (i === 0 && !cool) rhs[i] += theta * gl * Tleft;
      if (i === N - 1) rhs[i] += theta * gr * Tamb;
    }
    const Tn = solveTridiagonal(lower, diag, upper, rhs);
    lost += h * (theta * qOut(Tn) + (1 - theta) * qOut(T));
    if (!cool) supplied += h * (theta * qIn(Tn) + (1 - theta) * qIn(T));
    T = Float64Array.from(Tn);
    if (s % stride === 0 || s === nSteps) {
      out.t.push(s * h);
      out.Tcontents.push(cool ? T[0] : Tleft);
      out.TwallInner.push(innerSurface(T));
      out.Touter.push(outerSurface(T));
      out.qLoss.push(qOut(T));
      out.qFluid.push(qIn(T));
    }
  }
  return {
    ...out, r: cells.map((c) => c.rm), Tfinal: Array.from(T.subarray(off)),
    energyStored0: E0, energyStored: energy(T), energyLost: lost, energyIn: supplied,
    Usteady: 1 / (steady.Rtot * Math.PI * Di),
  };
}

/**
 * Lumped multi-capacity cooldown: one thermal node per layer plus the contents (see `radialTransient`, which this
 * calls with one cell per layer). Reference and validity as `radialTransient`; adequate when each layer has a small
 * Biot number, otherwise refine with `radialTransient`.
 * @param {object} o same arguments as `radialTransient` (mode fixed to 'cooldown', `cellsPerLayer` to 1)
 * @returns {object} as `radialTransient`
 */
export function cooldownMultiCapacity(o) {
  return radialTransient({ ...o, mode: 'cooldown', cellsPerLayer: 1, layers: o.layers.map((l) => ({ ...l, n: 1 })) });
}

/**
 * Time at which a monotonically evolving temperature history first crosses a threshold ("no-touch time" to the
 * hydrate equilibrium temperature or the wax appearance temperature), by linear interpolation between samples.
 * Reference: definition. @param {object} o
 * @param {number[]} o.t times, s @param {number[]} o.T temperatures, K @param {number} o.threshold K
 * @returns {number} time s; 0 if the first sample is already beyond the threshold; Infinity if never crossed
 */
export function timeToTemperature({ t, T, threshold }) {
  if (!t || !T || t.length !== T.length || t.length === 0) throw new RangeError('t and T must be equal-length arrays');
  need('threshold', threshold, 0, 5000);
  const cooling = T[T.length - 1] <= T[0];
  const beyond = (v) => (cooling ? v <= threshold : v >= threshold);
  if (beyond(T[0])) return 0;
  for (let i = 1; i < T.length; i++) {
    if (beyond(T[i])) return t[i - 1] + (t[i] - t[i - 1]) * (threshold - T[i - 1]) / (T[i] - T[i - 1]);
  }
  return Infinity;
}

/**
 * Heat transfer across the gas-filled annulus of a pipe-in-pipe: radiation between concentric grey surfaces,
 * conduction through the gas and natural convection by the Raithby–Hollands effective conductivity.
 *   q_rad = σ π D_1 (T_1⁴ − T_2⁴)/[1/ε_1 + (D_1/D_2)(1/ε_2 − 1)]
 *   q_cond = 2π k_eff (T_1 − T_2)/ln(D_2/D_1)
 *   k_eff/k = 0.386 [Pr/(0.861 + Pr)]^¼ Ra_c^¼ (not less than 1),
 *   Ra_c = [ln(D_2/D_1)]⁴ Ra_L / [L³ (D_1^(−3/5) + D_2^(−3/5))⁵], L = (D_2 − D_1)/2, Ra_L = g β ΔT L³ ρ² cp/(μ k).
 * Reference: Raithby & Hollands (1975) as given in Incropera et al. ch. 9 (valid 100 < Ra_c < 1e7); grey-body
 * concentric-cylinder exchange, Incropera et al. ch. 13.
 * @param {object} o
 * @param {number} o.D1 outer diameter of the inner pipe (or its coating), m
 * @param {number} o.D2 inner diameter of the outer pipe, m
 * @param {number} o.T1 inner surface temperature, K @param {number} o.T2 outer surface temperature, K
 * @param {number} [o.eps1=0.8] @param {number} [o.eps2=0.8] surface emissivities, –
 * @param {{ k: number, rho: number, mu: number, cp: number, beta?: number }} o.gas annulus gas properties at the
 *   mean temperature (β = 1/T_mean when omitted, ideal gas)
 * @param {number} [o.g=9.80665]
 * @returns {{ q: number, qRad: number, qCond: number, kEff: number, RaC: number, hEquivalent: number,
 *   kEquivalent: number }} heat flows W/m; hEquivalent = q/(π D_1 ΔT) W/(m²·K); kEquivalent = conductivity of a
 *   solid layer giving the same total heat flow, W/(m·K)
 */
export function pipeInPipeAnnulus({ D1, D2, T1, T2, eps1 = 0.8, eps2 = 0.8, gas, g = G0 }) {
  need('D1', D1, 1e-3, 20);
  need('D2', D2, D1 * (1 + 1e-9), 40);
  need('T1', T1, 1, 3000);
  need('T2', T2, 1, 3000);
  need('eps1', eps1, 1e-3, 1);
  need('eps2', eps2, 1e-3, 1);
  if (!gas) throw new RangeError('gas properties are required');
  const dT = T1 - T2;
  const qRad = SIGMA_SB * Math.PI * D1 * (T1 ** 4 - T2 ** 4) / (1 / eps1 + (D1 / D2) * (1 / eps2 - 1));
  const L = 0.5 * (D2 - D1);
  const beta = gas.beta ?? 2 / (T1 + T2);
  const Pr = gas.cp * gas.mu / gas.k;
  const RaL = g * beta * Math.abs(dT) * L ** 3 * gas.rho * gas.rho * gas.cp / (gas.mu * gas.k);
  const ln = Math.log(D2 / D1);
  const RaC = ln ** 4 * RaL / (L ** 3 * (D1 ** -0.6 + D2 ** -0.6) ** 5);
  const kEff = gas.k * Math.max(1, 0.386 * (Pr / (0.861 + Pr)) ** 0.25 * RaC ** 0.25);
  const qCond = 2 * Math.PI * kEff * dT / ln;
  const q = qRad + qCond;
  const hEq = dT !== 0 ? q / (Math.PI * D1 * dT) : 0;
  return { q, qRad, qCond, kEff, RaC, hEquivalent: hEq, kEquivalent: dT !== 0 ? q * ln / (2 * Math.PI * dT) : kEff };
}

/**
 * Direct electrical heating (or any distributed heating) power demand.
 * Steady: heat to hold the contents at `Ttarget`, q' = U π D (T_target − T_amb).
 * Warm-up within `warmupTime`: additionally the sensible heat Σ(m'c)(T_target − T_start)/t, with the mean loss
 * during the ramp taken at the average of start and target temperature (linear ramp).
 * Electrical input = thermal demand / efficiency (fraction of the supplied power dissipated in the pipe).
 * Reference: energy balance; system-specific efficiency per the DEH supplier (typically 0.5–0.8 for open-loop DEH).
 * @param {object} o
 * @param {number} o.U W/(m²·K) on πD @param {number} o.D m @param {number} o.L heated length, m
 * @param {number} o.Ttarget K @param {number} o.Tamb K
 * @param {number} [o.efficiency=1] –
 * @param {number} [o.heatCapacity] Σ(m'c) per unit length, J/(m·K) (warm-up)
 * @param {number} [o.Tstart] start temperature, K (default Tamb) @param {number} [o.warmupTime] s
 * @returns {{ holdPerLength: number, holdPower: number, warmupPerLength: number|null, warmupPower: number|null }}
 *   electrical power per unit length W/m and total W for holding and for warm-up
 */
export function heatingPower({ U, D, L, Ttarget, Tamb, efficiency = 1, heatCapacity, Tstart, warmupTime }) {
  need('U', U, 0, 1e6);
  need('D', D, 1e-4, 50);
  need('L', L, 0, 1e7);
  need('Ttarget', Ttarget, 1, 5000);
  need('Tamb', Tamb, 1, 5000);
  need('efficiency', efficiency, 1e-3, 1);
  const ua = U * Math.PI * D;
  const hold = ua * (Ttarget - Tamb) / efficiency;
  let warm = null;
  if (warmupTime != null) {
    need('warmupTime', warmupTime, 1e-6, 1e10);
    need('heatCapacity', heatCapacity, 0, 1e12);
    const Ts = Tstart ?? Tamb;
    warm = (heatCapacity * (Ttarget - Ts) / warmupTime + ua * (0.5 * (Ttarget + Ts) - Tamb)) / efficiency;
  }
  return { holdPerLength: hold, holdPower: hold * L, warmupPerLength: warm, warmupPower: warm == null ? null : warm * L };
}

/**
 * External thermal resistance of a single buried cable or pipe and the corresponding steady heat loss:
 * T4 = (ρ_T/2π) ln(u + √(u² − 1)), u = 2L/D_e  (identical to the cosh⁻¹ shape factor), q' = ΔT/T4.
 * Reference: IEC 60287-2-1 (external thermal resistance of a single isolated buried cable); Carslaw & Jaeger.
 * Validity: isothermal ground surface, homogeneous soil, u > 1.
 * @param {object} o
 * @param {number} o.depth burial depth to the axis, m @param {number} o.De external diameter, m
 * @param {number} o.rhoT soil thermal resistivity, K·m/W (= 1/k_soil)
 * @param {number} [o.dT] surface temperature rise above ambient soil, K
 * @returns {{ T4: number, q: number|null }} thermal resistance K·m/W and heat loss W/m (when dT is given)
 */
export function buriedThermalResistance({ depth, De, rhoT, dT }) {
  need('De', De, 1e-4, 20);
  need('depth', depth, De / 2 * (1 + 1e-9), 1e4);
  need('rhoT', rhoT, 1e-3, 100);
  const u = 2 * depth / De;
  const T4 = rhoT / (2 * Math.PI) * Math.log(u + Math.sqrt(u * u - 1));
  return { T4, q: dT == null ? null : dT / T4 };
}
