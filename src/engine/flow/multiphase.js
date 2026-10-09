/**
 * Steady-state gas–liquid pressure-gradient models and a marching solver along an elevation profile:
 * Beggs & Brill (1973) with the Payne et al. (1979) corrections, homogeneous and drift-flux gradients,
 * mechanistic stratified (two-fluid) and slug (unit-cell) models, turndown curve and hydraulic capacity.
 *
 * Gradients named grad* are pressure LOSS per metre in the flow direction [Pa/m] (dP/ds = −gradTotal).
 * theta is the pipe inclination from horizontal [rad], positive for upward flow.
 */
import { brent } from '../math/roots.js';
import { G, checkArg, darcyFrictionFactor, driftFluxHoldup } from './driftflux.js';
import { beggsBrillPattern, flowRegime, stratifiedGeometry, stratifiedForces } from './flowregime.js';
import { slugUnitCell } from './slug.js';

// Beggs & Brill horizontal holdup HL(0) = a λ^b / Fr^c.
const BB_HOLDUP = {
  segregated: [0.98, 0.4846, 0.0868],
  intermittent: [0.845, 0.5351, 0.0173],
  distributed: [1.065, 0.5824, 0.0609],
};
// Inclination correction C = (1 − λ) ln(d λ^e Nlv^f Fr^g); distributed uphill has no correction.
const BB_UPHILL = {
  segregated: [0.011, -3.768, 3.539, -1.614],
  intermittent: [2.96, 0.305, -0.4473, 0.0978],
  distributed: null,
};
const BB_DOWNHILL = [4.70, -0.3692, 0.1244, -0.5056];

/**
 * Erosional velocity of API RP 14E: ve = C/√ρ with ve in ft/s, ρ in lb/ft³; in SI ve = 1.21990 C/√ρ
 * (0.3048 · √16.01846). C = 100 for continuous and 125 for intermittent solids-free service.
 * Reference: API RP 14E (1991, R2019) §2.5.a, Eq. 2.14. Empirical screening limit, not a physical erosion model.
 * @param {number} rho gas/liquid mixture density at flowing conditions (no-slip) [kg/m³]
 * @param {number} [C=100] empirical constant in the API field-unit form [-]
 * @returns {number} erosional velocity [m/s]
 */
export function erosionalVelocity(rho, C = 100) {
  checkArg('rho', rho, { min: 0, exclusiveMin: true });
  checkArg('C', C, { min: 0, exclusiveMin: true });
  return 0.3048 * Math.sqrt(16.018463) * C / Math.sqrt(rho);
}

function singlePhase({ v, rho, mu, D, theta, roughness, P, isGas }) {
  const gradFriction = darcyFrictionFactor(rho * v * D / mu, roughness / D) * rho * v * v / (2 * D);
  const gradGravity = rho * G * Math.sin(theta);
  const Ek = isGas && P ? rho * v * v / P : 0;
  return { gradFriction, gradGravity, Ek };
}

function finish(r, { vsl, vsg, rhoL, rhoG, P, erosionalC }) {
  const vm = vsl + vsg;
  const lambdaL = vsl / vm;
  const rhoNs = lambdaL * rhoL + (1 - lambdaL) * rhoG;
  const H = r.holdup;
  const rhoM = H * rhoL + (1 - H) * rhoG;
  let Ek = r.Ek;
  if (Ek === undefined) Ek = P ? rhoM * vm * vsg / P : 0;
  if (Ek >= 0.95) throw new Error(`Acceleration term Ek = ${Ek.toFixed(2)}: flow is near choking; the steady gradient is not valid`);
  const base = r.gradFriction + r.gradGravity;
  const gradTotal = base / (1 - Ek);
  const ve = erosionalVelocity(rhoNs, erosionalC);
  return {
    method: r.method, regime: r.regime, holdup: H, lambdaL, vm, vsl, vsg, vL: H > 0 ? vsl / H : 0, vG: H < 1 ? vsg / (1 - H) : 0,
    rhoM, rhoNs, gradFriction: r.gradFriction, gradGravity: r.gradGravity, gradAcceleration: gradTotal - base, gradTotal, Ek,
    erosionalVelocity: ve, erosionalRatio: vm / ve, detail: r.detail ?? null,
  };
}

/**
 * Beggs & Brill two-phase pressure gradient for any inclination, with the Payne et al. corrections.
 *   pattern from λ = vsl/vm and Fr = vm²/(gD) (see beggsBrillPattern)
 *   HL(0) = a λ^b / Fr^c ≥ λ:  segregated 0.98/0.4846/0.0868, intermittent 0.845/0.5351/0.0173,
 *                              distributed 1.065/0.5824/0.0609
 *   C = (1 − λ) ln(d λ^e Nlv^f Fr^g) ≥ 0,  Nlv = vsl (ρL/(g σ))^¼
 *       uphill segregated 0.011/−3.768/3.539/−1.614, intermittent 2.96/0.305/−0.4473/0.0978, distributed C = 0;
 *       downhill (all patterns) 4.70/−0.3692/0.1244/−0.5056
 *   ψ = 1 + C [sin(1.8θ) − ⅓ sin³(1.8θ)],  HL(θ) = HL(0) ψ; transition pattern: A HL,seg + (1 − A) HL,int with
 *       A = (L3 − Fr)/(L3 − L2)
 *   Payne: HL(θ) × 0.924 uphill, × 0.685 downhill; normalising friction factor fn from the Moody diagram at the
 *       actual roughness (smooth pipe in the original)
 *   friction: ftp = fn e^s, s = ln y / (−0.0523 + 3.182 ln y − 0.8725 ln²y + 0.01853 ln⁴y), y = λ/HL(θ)²,
 *       s = ln(2.2 y − 1.2) for 1 < y < 1.2;  gradFriction = ftp ρn vm²/(2D)
 *   gravity: ρs g sinθ, ρs = ρL HL + ρG (1 − HL);  acceleration: Ek = ρs vm vsg/P, total = (f + g)/(1 − Ek)
 * References: Beggs & Brill (1973), J. Pet. Technol. 25(5), 607–617; Payne, Palmer, Brill & Beggs (1979),
 * J. Pet. Technol. 31(9), 1198–1208; Brill & Mukherjee (1999), Multiphase Flow in Wells, SPE Monograph 17, §4.
 * Validity: 25–38 mm air–water test data, all inclinations; widely used outside that range. With Payne's
 * downhill factor the holdup of an almost gas-free downhill flow does not tend to 1 — single-phase input
 * (vsg = 0 or vsl = 0) is handled exactly as single-phase Darcy flow.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vsg superficial gas velocity [m/s], >= 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} o.sigma surface tension [N/m]
 * @param {number} [o.P] absolute pressure [Pa] (acceleration term; omitted → no acceleration)
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {boolean} [o.payne=true] apply the Payne et al. holdup factors
 * @param {number} [o.erosionalC=100] API RP 14E constant [-]
 * @returns {object} common gradient result (see pressureGradient) with detail = {pattern, holdup0, C, psi, Fr,
 *   Nlv, L1, L2, L3, L4, fn, ftp, s, y} [-]
 */
export function beggsBrill({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, sigma, P, roughness = 0, payne = true, erosionalC = 100 }) {
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vsg', vsg, { min: 0 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('muG', muG, { min: 0, exclusiveMin: true });
  checkArg('roughness', roughness, { min: 0 });
  const vm = vsl + vsg;
  checkArg('vsl + vsg', vm, { min: 0, exclusiveMin: true });
  const common = { vsl, vsg, rhoL, rhoG, P, erosionalC };
  if (vsg === 0 || vsl === 0) {
    const gas = vsl === 0;
    const sp = singlePhase({ v: vm, rho: gas ? rhoG : rhoL, mu: gas ? muG : muL, D, theta, roughness, P, isGas: gas });
    return finish({ method: 'beggs-brill', regime: gas ? 'single-phase-gas' : 'single-phase-liquid', holdup: gas ? 0 : 1, ...sp }, common);
  }
  checkArg('sigma', sigma, { min: 0, exclusiveMin: true });
  const pat = beggsBrillPattern({ vsl, vsg, D });
  const lam = pat.lambdaL;
  const Fr = pat.Fr;
  const Nlv = vsl * (rhoL / (G * sigma)) ** 0.25;
  const trig = Math.sin(1.8 * theta);
  const incl = trig - trig ** 3 / 3;
  const forPattern = (name) => {
    const [a, b, c] = BB_HOLDUP[name];
    const holdup0 = Math.min(1, Math.max(lam, a * lam ** b / Fr ** c));
    let C = 0;
    if (theta !== 0) {
      const k = theta > 0 ? BB_UPHILL[name] : BB_DOWNHILL;
      if (k) C = Math.max(0, (1 - lam) * Math.log(k[0] * lam ** k[1] * Nlv ** k[2] * Fr ** k[3]));
    }
    const psi = 1 + C * incl;
    return { holdup0, C, psi, holdup: holdup0 * psi };
  };
  let h;
  if (pat.pattern === 'transition') {
    const seg = forPattern('segregated');
    const int = forPattern('intermittent');
    const A = (pat.L3 - Fr) / (pat.L3 - pat.L2);
    h = {
      holdup0: A * seg.holdup0 + (1 - A) * int.holdup0, C: A * seg.C + (1 - A) * int.C, psi: A * seg.psi + (1 - A) * int.psi,
      holdup: A * seg.holdup + (1 - A) * int.holdup,
    };
  } else h = forPattern(pat.pattern);
  let HL = h.holdup;
  if (payne && theta !== 0) HL *= theta > 0 ? 0.924 : 0.685;
  // Uphill the gas cannot lag the liquid, so the no-slip value is a physical lower bound.
  if (theta > 0) HL = Math.max(HL, lam);
  HL = Math.min(1, Math.max(1e-9, HL));
  const rhoN = lam * rhoL + (1 - lam) * rhoG;
  const muN = lam * muL + (1 - lam) * muG;
  const fn = darcyFrictionFactor(rhoN * vm * D / muN, roughness / D);
  const y = lam / (HL * HL);
  let s;
  if (y > 1 && y < 1.2) s = Math.log(2.2 * y - 1.2);
  else {
    const ly = Math.log(y);
    s = ly === 0 ? 0 : ly / (-0.0523 + 3.182 * ly - 0.8725 * ly * ly + 0.01853 * ly ** 4);
  }
  const ftp = fn * Math.exp(s);
  const rhoS = HL * rhoL + (1 - HL) * rhoG;
  return finish({
    method: 'beggs-brill', regime: pat.pattern, holdup: HL, gradFriction: ftp * rhoN * vm * vm / (2 * D),
    gradGravity: rhoS * G * Math.sin(theta),
    detail: { pattern: pat.pattern, holdup0: h.holdup0, C: h.C, psi: h.psi, Fr, Nlv, L1: pat.L1, L2: pat.L2, L3: pat.L3, L4: pat.L4, fn, ftp, s, y },
  }, common);
}

/**
 * Steady stratified flow from the two-fluid momentum balances (liquid level solved from the combined balance,
 * pressure gradient from the sum of both phase balances).
 * Reference: Taitel & Dukler (1976), AIChE J. 22(1), 47–55; interfacial-friction options of stratifiedForces.
 * When several levels satisfy the balance (upward inclination) the lowest is returned.
 * Validity: stratified smooth/wavy flow, |θ| ≲ 30°.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {'smooth'|'constant'|'andritsos'} [o.interfacial='smooth']
 * @param {number} [o.fi] constant interfacial Fanning factor [-]
 * @returns {{hL: number, holdup: number, vL: number, vG: number, tauWL: number, tauWG: number, tauI: number,
 *   gradFriction: number, gradGravity: number, gradTotal: number, residual: number}} hL = level/D [-];
 *   holdup [-]; in-situ velocities [m/s]; shear stresses [Pa]; gradients [Pa/m]; momentum residual [Pa/m]
 */
export function stratifiedTwoFluid({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, roughness = 0, interfacial = 'smooth', fi }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0, exclusiveMin: true });
  const fluid = { D, theta, rhoL, rhoG, muL, muG, roughness, interfacial, fi };
  const at = (h) => {
    const H = stratifiedGeometry(h).holdup;
    return stratifiedForces({ ...fluid, h, vL: vsl / H, vG: vsg / (1 - H) });
  };
  const f = (h) => at(h).residual;
  const N = 300;
  let hPrev = 1e-6;
  let hL = null;
  if (f(hPrev) >= 0) hL = hPrev;
  for (let k = 1; k < N && hL === null; k++) {
    const h = 0.5 * (1 - Math.cos(Math.PI * k / N));
    if (h <= hPrev) continue;
    if (f(h) >= 0) hL = brent(f, hPrev, h, { tol: 1e-13 });
    hPrev = h;
  }
  if (hL === null) hL = brent(f, hPrev, 1 - 1e-9, { tol: 1e-13 });
  const r = at(hL);
  return {
    hL, holdup: r.holdup, vL: vsl / r.holdup, vG: vsg / (1 - r.holdup), tauWL: r.tauWL, tauWG: r.tauWG, tauI: r.tauI,
    gradFriction: r.gradFriction, gradGravity: r.gradGravity, gradTotal: r.gradFriction + r.gradGravity, residual: r.residual,
  };
}

/**
 * Two-phase pressure gradient and holdup at a point by a named method.
 *   'beggs-brill'  Beggs & Brill (1973) + Payne et al. (1979)
 *   'homogeneous'  no slip; f(Re = ρn vm D/μn) ρn vm²/(2D) + ρn g sinθ
 *   'drift-flux'   holdup from a drift-flux closure (o.closure, default 'bendiksen'; o.C0/o.vd overrides);
 *                  friction f(Re = ρm vm D/μm) ρm vm²/(2D) on the slip mixture density and holdup-weighted
 *                  viscosity — the same closure set as the transient solver
 *   'mechanistic'  flow pattern from flowRegime(); stratified → two-fluid balance; intermittent with
 *                  |θ| ≤ 15° → slug unit cell; slug/churn at steeper angles → drift flux (Bendiksen);
 *                  bubble → drift flux (Harmathy); dispersed bubble → homogeneous; annular → Beggs & Brill
 *                  (no mechanistic annular model is implemented; `detail.fallback` says so)
 * Acceleration: Ek = ρm vm vsg / P (isothermal gas expansion at constant mass flux), total = (f + g)/(1 − Ek).
 * @param {object} o
 * @param {'beggs-brill'|'homogeneous'|'drift-flux'|'mechanistic'} [o.method='beggs-brill']
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vsg superficial gas velocity [m/s], >= 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.sigma] surface tension [N/m] (Beggs–Brill, bubbly closure, vertical map)
 * @param {number} [o.P] absolute pressure [Pa] (acceleration term)
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {string} [o.closure='bendiksen'] drift-flux closure name
 * @param {number} [o.C0] drift-flux C0 override [-]
 * @param {number} [o.vd] drift-flux vd override [m/s]
 * @param {string} [o.interfacial] stratified/film interfacial friction option
 * @param {number} [o.fi] constant interfacial Fanning factor [-]
 * @param {number} [o.slugLength] slug-body length for the unit cell [m]
 * @param {boolean} [o.payne=true] Payne corrections in Beggs–Brill
 * @param {number} [o.erosionalC=100] API RP 14E constant [-]
 * @returns {{method: string, regime: string, holdup: number, lambdaL: number, vm: number, vsl: number,
 *   vsg: number, vL: number, vG: number, rhoM: number, rhoNs: number, gradFriction: number,
 *   gradGravity: number, gradAcceleration: number, gradTotal: number, Ek: number, erosionalVelocity: number,
 *   erosionalRatio: number, detail: object|null}} holdup, lambdaL, Ek, erosionalRatio = vm/ve [-]; velocities
 *   [m/s]; densities [kg/m³]; gradients = pressure loss per metre [Pa/m]
 */
export function pressureGradient(o) {
  const { method = 'beggs-brill', vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, sigma, P, roughness = 0, closure = 'bendiksen',
    C0, vd, interfacial, fi, slugLength, payne = true, erosionalC = 100 } = o;
  if (method === 'beggs-brill') return beggsBrill({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, sigma, P, roughness, payne, erosionalC });
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vsg', vsg, { min: 0 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, exclusiveMin: true });
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('muG', muG, { min: 0, exclusiveMin: true });
  const vm = vsl + vsg;
  checkArg('vsl + vsg', vm, { min: 0, exclusiveMin: true });
  const common = { vsl, vsg, rhoL, rhoG, P, erosionalC };
  if (vsg === 0 || vsl === 0) {
    const gas = vsl === 0;
    const sp = singlePhase({ v: vm, rho: gas ? rhoG : rhoL, mu: gas ? muG : muL, D, theta, roughness, P, isGas: gas });
    return finish({ method, regime: gas ? 'single-phase-gas' : 'single-phase-liquid', holdup: gas ? 0 : 1, ...sp }, common);
  }
  const sinT = Math.sin(theta);
  const mixture = (H, regime, detail) => {
    const rho = H * rhoL + (1 - H) * rhoG;
    const mu = H * muL + (1 - H) * muG;
    return {
      method, regime, holdup: H, detail,
      gradFriction: darcyFrictionFactor(rho * vm * D / mu, roughness / D) * rho * vm * vm / (2 * D), gradGravity: rho * G * sinT,
    };
  };
  const drift = (name, regime, detail) => {
    const d = driftFluxHoldup({ vsl, vsg, D, theta, rhoL, rhoG, sigma, regime: name, C0, vd });
    return mixture(d.holdup, regime, { ...detail, C0: d.C0, vd: d.vd, clipped: d.clipped });
  };
  switch (method) {
    case 'homogeneous': return finish(mixture(vsl / vm, 'homogeneous', null), common);
    case 'drift-flux': return finish(drift(closure, `drift-flux (${closure})`, {}), common);
    case 'mechanistic': {
      const fr = flowRegime({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, sigma });
      const regime = fr.regime;
      if (regime.startsWith('stratified')) {
        const s = stratifiedTwoFluid({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, roughness, interfacial: interfacial ?? (regime === 'stratified-wavy' ? 'andritsos' : 'smooth'), fi });
        return finish({ method, regime, holdup: s.holdup, gradFriction: s.gradFriction, gradGravity: s.gradGravity, detail: { model: 'two-fluid stratified', hL: s.hL } }, common);
      }
      if (regime === 'intermittent' || regime === 'slug' || regime === 'churn') {
        if (Math.abs(theta) <= 15 * Math.PI / 180) {
          const u = slugUnitCell({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, roughness, slugLength, interfacial, fi });
          if (u.feasible) {
            return finish({
              method, regime, holdup: u.holdup, gradFriction: u.gradFriction, gradGravity: u.gradGravity,
              detail: { model: 'slug unit cell', vt: u.vt, slugHoldup: u.slugHoldup, filmHoldup: u.filmHoldup, slugFraction: u.slugFraction, frequency: u.frequency },
            }, common);
          }
        }
        return finish(drift('bendiksen', regime, { model: 'drift flux (Bendiksen)' }), common);
      }
      if (regime === 'bubble') return finish(drift('bubbly', regime, { model: 'drift flux (Harmathy)' }), common);
      if (regime === 'dispersed-bubble') return finish(mixture(vsl / vm, regime, { model: 'homogeneous' }), common);
      const b = beggsBrill({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, sigma, P, roughness, payne, erosionalC });
      return { ...b, method, regime, detail: { ...b.detail, model: 'Beggs & Brill', fallback: 'no mechanistic annular model' } };
    }
    default: throw new RangeError(`method = ${method} is not a known pressure-gradient method`);
  }
}

/**
 * Convenience wrapper: local gradient from mass rates and a fluid-property callback.
 * @param {object} o any option of pressureGradient (method, D, theta, roughness, …) plus
 * @param {(P: number, T: number) => {rhoL: number, rhoG: number, muL: number, muG: number, sigma: number,
 *   quality?: number}} o.fluidAt properties at pressure P [Pa] and temperature T [K]: densities [kg/m³],
 *   viscosities [Pa·s], surface tension [N/m], optional gas mass fraction `quality` [-] (flash result)
 * @param {number} o.P absolute pressure [Pa]
 * @param {number} o.T temperature [K]
 * @param {number} [o.massRateL] liquid mass rate [kg/s]
 * @param {number} [o.massRateG] gas mass rate [kg/s]
 * @param {number} [o.massRate] total mass rate [kg/s] (with `gasMassFraction` or fluidAt().quality)
 * @param {number} [o.gasMassFraction] gas mass fraction [-]
 * @returns {object} pressureGradient result plus `fluid` (the callback result) and massRateL, massRateG [kg/s]
 */
export function gradientAt({ fluidAt, P, T, massRateL, massRateG, massRate, gasMassFraction, ...rest }) {
  if (typeof fluidAt !== 'function') throw new RangeError('fluidAt must be a function (P, T) => properties');
  checkArg('P', P, { min: 0, exclusiveMin: true });
  const fl = fluidAt(P, T);
  let mL = massRateL;
  let mG = massRateG;
  if (massRate !== undefined) {
    const x = fl.quality ?? gasMassFraction;
    checkArg('gasMassFraction', x, { min: 0, max: 1 });
    mG = massRate * x; mL = massRate * (1 - x);
  }
  checkArg('massRateL', mL, { min: 0 });
  checkArg('massRateG', mG, { min: 0 });
  const A = Math.PI * rest.D * rest.D / 4;
  const r = pressureGradient({ ...rest, ...fl, P, vsl: mL / (fl.rhoL * A), vsg: mG / (fl.rhoG * A) });
  return { ...r, fluid: fl, massRateL: mL, massRateG: mG };
}

function buildStations(profile, maxStep) {
  if (!Array.isArray(profile) || profile.length < 2) throw new RangeError('profile must list at least two {x, z} points');
  let total = 0;
  for (let i = 0; i < profile.length; i++) {
    checkArg(`profile[${i}].x`, profile[i].x);
    checkArg(`profile[${i}].z`, profile[i].z);
    if (i > 0) {
      const len = Math.hypot(profile[i].x - profile[i - 1].x, profile[i].z - profile[i - 1].z);
      if (!(len > 0)) throw new RangeError(`profile[${i}] coincides with the previous point`);
      total += len;
    }
  }
  const step = maxStep ?? total / 100;
  checkArg('maxStep', step, { min: 0, exclusiveMin: true });
  const st = [{ s: 0, x: profile[0].x, z: profile[0].z, theta: 0, pointIndex: 0 }];
  let s = 0;
  for (let i = 1; i < profile.length; i++) {
    const dx = profile[i].x - profile[i - 1].x;
    const dz = profile[i].z - profile[i - 1].z;
    const len = Math.hypot(dx, dz);
    const theta = Math.asin(Math.max(-1, Math.min(1, dz / len)));
    const n = Math.min(5000, Math.max(1, Math.ceil(len / step - 1e-9)));
    for (let k = 1; k <= n; k++) {
      st.push({ s: s + len * k / n, x: profile[i - 1].x + dx * k / n, z: profile[i - 1].z + dz * k / n, theta, pointIndex: k === n ? i : -1 });
    }
    s += len;
  }
  st[0].theta = st[1].theta;
  if (st.length > 20001) throw new RangeError('maxStep gives more than 20000 stations');
  return { stations: st, length: total };
}

function temperatureFunction(temperature, profile, stations) {
  if (typeof temperature === 'function') return temperature;
  if (typeof temperature === 'number') { checkArg('temperature', temperature, { min: 0, exclusiveMin: true }); return () => temperature; }
  if (Array.isArray(temperature) || ArrayBuffer.isView(temperature)) {
    if (temperature.length !== profile.length) throw new RangeError('temperature array must have one value per profile point');
    const sPts = stations.filter((p) => p.pointIndex >= 0).map((p) => p.s);
    return (s) => {
      let i = 1;
      while (i < sPts.length - 1 && s > sPts[i]) i++;
      const w = Math.min(1, Math.max(0, (s - sPts[i - 1]) / (sPts[i] - sPts[i - 1])));
      return temperature[i - 1] * (1 - w) + temperature[i] * w;
    };
  }
  throw new RangeError('temperature must be a number [K], an array per profile point or a function of chainage');
}

/**
 * Steady two-phase flow along an elevation profile by marching the pressure gradient (classical 4th-order
 * Runge–Kutta on dP/ds, properties re-evaluated at the local pressure and temperature through `fluidAt`).
 *
 * Boundary specification — exactly two of {inlet pressure, outlet pressure, flow rate}:
 *   pIn + rate   march downstream from the inlet;
 *   pOut + rate  march upstream from the outlet (backward integration of the same ODE — the converged result of
 *                a shooting iteration on the inlet pressure);
 *   pIn + pOut   the total mass rate is solved (Brent) at the stated gas mass fraction. The ΔP–rate curve of an
 *                uphill line has a minimum, so two rates can satisfy the same pressures; the friction-dominated
 *                (highest) root is returned and `multipleSolutions` is set.
 * Any other combination throws an Error naming the over- or under-constraint.
 *
 * @param {object} o
 * @param {{x: number, z: number}[]} o.profile horizontal distance x [m] and elevation z [m, up] of the route
 *   points, inlet first
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {(P: number, T: number) => object} o.fluidAt property callback (see gradientAt)
 * @param {number|number[]|((s: number) => number)} o.temperature fluid temperature [K]: constant, one value per
 *   profile point (linear in chainage between points) or a function of chainage s [m]
 * @param {number} [o.pIn] inlet absolute pressure [Pa]
 * @param {number} [o.pOut] outlet absolute pressure [Pa]
 * @param {number} [o.massRateL] liquid mass rate [kg/s] (no interphase mass transfer)
 * @param {number} [o.massRateG] gas mass rate [kg/s]
 * @param {number} [o.massRate] total mass rate [kg/s], split by fluidAt().quality or `gasMassFraction`
 * @param {number} [o.gasMassFraction] gas mass fraction [-] (needed for the rate solve unless fluidAt returns quality)
 * @param {string} [o.method='beggs-brill'] gradient method (see pressureGradient) and its options (closure, C0, vd,
 *   interfacial, fi, slugLength, payne, erosionalC) passed through
 * @param {number} [o.maxStep=length/100] largest integration step [m]
 * @param {[number, number]} [o.rateBracket] search range of the total mass rate for the rate solve [kg/s]
 * @returns {{stations: object[], pIn: number, pOut: number, dP: number, dPFriction: number, dPGravity: number,
 *   dPAcceleration: number, massRateL: number, massRateG: number, massRate: number, length: number,
 *   liquidInventory: number, meanHoldup: number, maxErosionalRatio: number, solvedFor: string,
 *   multipleSolutions: boolean}} stations: per station {s, x, z [m]; theta [rad]; P [Pa]; T [K]; regime;
 *   holdup, lambdaL [-]; vsl, vsg, vm, vL, vG [m/s]; rhoM [kg/m³]; gradFriction, gradGravity, gradAcceleration,
 *   gradTotal [Pa/m]; erosionalRatio [-]}; pressures and pressure-drop components [Pa] (dP = pIn − pOut);
 *   mass rates [kg/s]; length [m]; liquidInventory [m³]; solvedFor ∈ 'pOut' | 'pIn' | 'massRate'
 */
export function marchPipeline(o) {
  const { profile, D, roughness = 0, fluidAt, temperature, pIn, pOut, massRateL, massRateG, massRate, gasMassFraction,
    maxStep, rateBracket, ...methodOptions } = o;
  checkArg('D', D, { min: 0, exclusiveMin: true });
  if (typeof fluidAt !== 'function') throw new RangeError('fluidAt must be a function (P, T) => properties');
  const hasPhaseRates = massRateL !== undefined || massRateG !== undefined;
  if (hasPhaseRates && massRate !== undefined) throw new Error('Over-constrained flow rate: give either massRate or massRateL/massRateG, not both.');
  if (hasPhaseRates && (massRateL === undefined || massRateG === undefined)) throw new Error('Under-constrained flow rate: both massRateL and massRateG are needed.');
  const hasRate = hasPhaseRates || massRate !== undefined;
  const nGiven = (pIn !== undefined) + (pOut !== undefined) + hasRate;
  if (nGiven === 3) throw new Error('Over-constrained hydraulic boundary conditions: inlet pressure, outlet pressure and flow rate were all specified. Specify exactly two; the third is the solution.');
  if (nGiven < 2) throw new Error('Under-constrained hydraulic boundary conditions: specify exactly two of inlet pressure, outlet pressure and flow rate.');
  if (pIn !== undefined) checkArg('pIn', pIn, { min: 0, exclusiveMin: true });
  if (pOut !== undefined) checkArg('pOut', pOut, { min: 0, exclusiveMin: true });
  const { stations, length } = buildStations(profile, maxStep);
  const Tat = temperatureFunction(temperature, profile, stations);
  const A = Math.PI * D * D / 4;
  const n = stations.length;

  const local = (P, s, theta, rate) => {
    if (!(P > 0)) throw new Error(`Pressure falls to zero at chainage ${s.toFixed(1)} m: the flow rate cannot be delivered with the stated boundary pressure.`);
    return gradientAt({ ...methodOptions, D, roughness, theta, fluidAt, P, T: Tat(s), ...rate });
  };
  // forward = true: from the inlet with P0 = pIn; false: from the outlet with P0 = pOut.
  const march = (P0, forward, rate) => {
    const P = new Float64Array(n);
    if (forward) {
      P[0] = P0;
      for (let i = 0; i < n - 1; i++) {
        const a = stations[i];
        const b = stations[i + 1];
        const h = b.s - a.s;
        const f = (p, s) => -local(p, s, b.theta, rate).gradTotal;
        const k1 = f(P[i], a.s);
        const k2 = f(P[i] + h * k1 / 2, a.s + h / 2);
        const k3 = f(P[i] + h * k2 / 2, a.s + h / 2);
        const k4 = f(P[i] + h * k3, b.s);
        P[i + 1] = P[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
      }
    } else {
      P[n - 1] = P0;
      for (let i = n - 1; i > 0; i--) {
        const a = stations[i - 1];
        const b = stations[i];
        const h = b.s - a.s;
        const f = (p, s) => local(p, s, b.theta, rate).gradTotal;
        const k1 = f(P[i], b.s);
        const k2 = f(P[i] + h * k1 / 2, b.s - h / 2);
        const k3 = f(P[i] + h * k2 / 2, b.s - h / 2);
        const k4 = f(P[i] + h * k3, a.s);
        P[i - 1] = P[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
      }
    }
    if (!(P[0] > 0) || !(P[n - 1] > 0) || !Number.isFinite(P[0] + P[n - 1])) throw new Error('Pressure falls to zero along the line: the flow rate cannot be delivered with the stated boundary pressure.');
    return P;
  };

  let rate = hasPhaseRates ? { massRateL, massRateG } : { massRate, gasMassFraction };
  let P;
  let solvedFor;
  let multipleSolutions = false;
  if (hasRate) {
    if (hasPhaseRates) { checkArg('massRateL', massRateL, { min: 0 }); checkArg('massRateG', massRateG, { min: 0 }); }
    else checkArg('massRate', massRate, { min: 0, exclusiveMin: true });
    if (pIn !== undefined) { P = march(pIn, true, rate); solvedFor = 'pOut'; } else { P = march(pOut, false, rate); solvedFor = 'pIn'; }
  } else {
    if (!(pIn > 0 && pOut > 0)) throw new RangeError('pIn and pOut must be positive');
    solvedFor = 'massRate';
    // Residual: inlet pressure needed for rate m (marching upstream from the outlet) minus the available one.
    const need = (m) => {
      try { return march(pOut, false, { massRate: m, gasMassFraction })[0] - pIn; } catch (e) { if (e instanceof RangeError) throw e; return Infinity; }
    };
    const fl = fluidAt(pIn, Tat(0));
    const x = fl.quality ?? gasMassFraction;
    checkArg('gasMassFraction', x, { min: 0, max: 1 });
    const rhoH = 1 / (x / fl.rhoG + (1 - x) / fl.rhoL);
    const lo = rateBracket ? rateBracket[0] : 1e-3 * rhoH * A;
    const hi = rateBracket ? rateBracket[1] : 100 * rhoH * A;
    checkArg('rateBracket[0]', lo, { min: 0, exclusiveMin: true });
    checkArg('rateBracket[1]', hi, { min: lo, exclusiveMin: true });
    const N = 48;
    const ms = Array.from({ length: N + 1 }, (_, i) => lo * (hi / lo) ** (i / N));
    const rs = ms.map(need);
    let k = -1;
    let crossings = 0;
    for (let i = 0; i < N; i++) {
      if ((rs[i] < 0) !== (rs[i + 1] < 0)) { crossings++; if (rs[i] < 0 && rs[i + 1] >= 0) k = i; }
    }
    if (k < 0) {
      if (rs.every((r) => r >= 0)) throw new Error('No flow solution: the available pressure difference is below the minimum (static head plus friction) of the line at every rate in the search range.');
      throw new Error('No flow solution within the rate search range: the available pressure difference exceeds the line resistance at the upper bracket; widen rateBracket.');
    }
    multipleSolutions = crossings > 1;
    let b = ms[k + 1];
    // Shrink a bracket whose upper end failed (non-deliverable rate) until the residual is finite.
    while (!Number.isFinite(need(b))) b = Math.sqrt(ms[k] * b);
    const m = Math.exp(brent((lm) => need(Math.exp(lm)), Math.log(ms[k]), Math.log(b), { tol: 1e-12 }));
    rate = { massRate: m, gasMassFraction };
    P = march(pOut, false, rate);
  }

  const out = new Array(n);
  let inventory = 0;
  let dPf = 0;
  let dPg = 0;
  let dPa = 0;
  let maxEr = 0;
  let prev = null;
  let mL = 0;
  let mG = 0;
  for (let i = 0; i < n; i++) {
    const s = stations[i];
    const g = local(P[i], s.s, s.theta, rate);
    out[i] = {
      s: s.s, x: s.x, z: s.z, theta: s.theta, P: P[i], T: Tat(s.s), regime: g.regime, holdup: g.holdup, lambdaL: g.lambdaL,
      vsl: g.vsl, vsg: g.vsg, vm: g.vm, vL: g.vL, vG: g.vG, rhoM: g.rhoM, gradFriction: g.gradFriction, gradGravity: g.gradGravity,
      gradAcceleration: g.gradAcceleration, gradTotal: g.gradTotal, erosionalRatio: g.erosionalRatio,
    };
    if (i === 0) { mL = g.massRateL; mG = g.massRateG; }
    if (g.erosionalRatio > maxEr) maxEr = g.erosionalRatio;
    if (i > 0) {
      // The inclination changes at route points: integrate each step with the upstream end evaluated at the
      // inclination of this step.
      const a = stations[i - 1];
      const up = a.theta === s.theta && prev ? prev : local(P[i - 1], a.s, s.theta, rate);
      const h = s.s - a.s;
      inventory += 0.5 * (up.holdup + g.holdup) * A * h;
      dPf += 0.5 * (up.gradFriction + g.gradFriction) * h;
      dPg += 0.5 * (up.gradGravity + g.gradGravity) * h;
      dPa += 0.5 * (up.gradAcceleration + g.gradAcceleration) * h;
    }
    prev = g;
  }
  return {
    stations: out, pIn: P[0], pOut: P[n - 1], dP: P[0] - P[n - 1], dPFriction: dPf, dPGravity: dPg, dPAcceleration: dPa,
    massRateL: mL, massRateG: mG, massRate: mL + mG, length, liquidInventory: inventory, meanHoldup: inventory / (A * length),
    maxErosionalRatio: maxEr, solvedFor, multipleSolutions,
  };
}

/**
 * Turndown curve: pressure drop against total mass rate at a fixed outlet pressure and gas mass fraction. In
 * an uphill line the curve has a minimum separating the gravity-dominated branch (liquid accumulates as the
 * rate falls; operation is unstable and prone to terrain slugging) from the friction-dominated branch.
 * @param {object} o marchPipeline arguments without pIn and without rates, plus
 * @param {number} o.pOut outlet pressure [Pa]
 * @param {number[]} o.massRates total mass rates to evaluate [kg/s], ascending
 * @param {number} [o.gasMassFraction] gas mass fraction [-] (unless fluidAt returns quality)
 * @returns {{massRate: number[], dP: number[], pIn: number[], dPFriction: number[], dPGravity: number[],
 *   liquidInventory: number[], minimum: {massRate: number, dP: number, index: number}, hasInteriorMinimum: boolean}}
 *   rates [kg/s]; pressures [Pa]; inventories [m³]; the sampled minimum of dP and whether it lies inside the range
 */
export function turndownCurve({ massRates, pIn, massRate, massRateL, massRateG, ...rest }) {
  if (pIn !== undefined || massRate !== undefined || massRateL !== undefined || massRateG !== undefined) {
    throw new Error('Over-constrained turndown specification: give the outlet pressure and the list massRates only.');
  }
  if (!Array.isArray(massRates) || massRates.length < 2) throw new RangeError('massRates must list at least two rates');
  const res = { massRate: [], dP: [], pIn: [], dPFriction: [], dPGravity: [], liquidInventory: [] };
  let iMin = 0;
  massRates.forEach((m, i) => {
    const r = marchPipeline({ ...rest, massRate: m });
    res.massRate.push(m); res.dP.push(r.dP); res.pIn.push(r.pIn); res.dPFriction.push(r.dPFriction); res.dPGravity.push(r.dPGravity);
    res.liquidInventory.push(r.liquidInventory);
    if (r.dP < res.dP[iMin]) iMin = i;
  });
  return { ...res, minimum: { massRate: res.massRate[iMin], dP: res.dP[iMin], index: iMin }, hasInteriorMinimum: iMin > 0 && iMin < massRates.length - 1 };
}

/**
 * Hydraulic capacity: the largest total mass rate the line delivers between a maximum inlet pressure and a
 * minimum arrival pressure (friction-dominated root of the rate solve in marchPipeline).
 * @param {object} o marchPipeline arguments without rates, plus
 * @param {number} o.pIn maximum inlet pressure [Pa]
 * @param {number} o.pOut minimum outlet (arrival) pressure [Pa]
 * @param {number} [o.gasMassFraction] gas mass fraction [-]
 * @returns {{massRate: number, massRateL: number, massRateG: number, maxErosionalRatio: number,
 *   multipleSolutions: boolean, result: object}} rates [kg/s]; erosional ratio [-]; full marchPipeline result
 */
export function hydraulicCapacity(o) {
  if (o.massRate !== undefined || o.massRateL !== undefined || o.massRateG !== undefined) {
    throw new Error('Over-constrained capacity specification: give both pressures and no flow rate.');
  }
  const r = marchPipeline(o);
  return { massRate: r.massRate, massRateL: r.massRateL, massRateG: r.massRateG, maxErosionalRatio: r.maxErosionalRatio, multipleSolutions: r.multipleSolutions, result: r };
}
