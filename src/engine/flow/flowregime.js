/**
 * Gas–liquid flow-pattern prediction: Taitel & Dukler (1976) mechanistic map for horizontal and slightly
 * inclined pipes, Taitel, Barnea & Dukler (1980) map for upward vertical flow, the Beggs & Brill (1973)
 * horizontal pattern boundaries, and a regime-map generator for plotting.
 *
 * Also the dimensional shear/gravity balance of a stratified layer, shared by the two-fluid, slug-film and
 * stability models.
 *
 * theta is the pipe inclination from horizontal [rad], positive for upward flow.
 */
import { brent } from '../math/roots.js';
import { G, checkArg, darcyFrictionFactor } from './driftflux.js';

const QUARTER_PI = Math.PI / 4;

/**
 * Dimensionless geometry of stratified flow in a circular pipe (lengths scaled by D, areas by D²,
 * velocities by the superficial velocity of the phase).
 * Reference: Taitel & Dukler (1976), AIChE J. 22(1), 47–55, Eqs. (11)–(18).
 * @param {number} h dimensionless liquid level hL/D [-], 0 < h < 1
 * @returns {{AL: number, AG: number, SL: number, SG: number, Si: number, dAL: number, uL: number, uG: number,
 *   DL: number, DG: number, holdup: number}} areas ÃL, ÃG; wetted perimeters S̃L, S̃G; interface width S̃i;
 *   dÃL/dh̃; velocity ratios ũL = vL/vsl, ũG = vG/vsg; hydraulic diameters D̃L, D̃G; holdup = ÃL/(π/4) — all [-]
 */
export function stratifiedGeometry(h) {
  checkArg('h', h, { min: 0, max: 1 });
  const c = 2 * h - 1;
  const ac = Math.acos(c);
  const s = Math.sqrt(Math.max(0, 1 - c * c));
  const AL = 0.25 * (Math.PI - ac + c * s);
  const AG = 0.25 * (ac - c * s);
  const SL = Math.PI - ac;
  return {
    AL, AG, SL, SG: ac, Si: s, dAL: s, uL: QUARTER_PI / AL, uG: QUARTER_PI / AG,
    DL: 4 * AL / SL, DG: 4 * AG / (ac + s), holdup: AL / QUARTER_PI,
  };
}

/**
 * Residual of the dimensionless combined momentum balance of equilibrium stratified flow,
 *   X² [(ũL D̃L)^(−n) ũL² S̃L/ÃL] − [(ũG D̃G)^(−m) ũG² (S̃G/ÃG + S̃i/ÃL + S̃i/ÃG)] − 4Y,
 * which is zero at the equilibrium level (fi = fG, interface velocity ≪ gas velocity).
 * Reference: Taitel & Dukler (1976), Eq. (9). Y is defined here with theta positive upward, i.e.
 * Y = −(ρL − ρG) g sinθ / |dP/dx|SG (the paper's α is positive downward).
 * @param {number} h dimensionless level hL/D [-]
 * @param {number} X2 Lockhart–Martinelli parameter squared [-]
 * @param {number} Y inclination parameter [-]
 * @param {number} [n=0.2] liquid friction exponent [-]
 * @param {number} [m=0.2] gas friction exponent [-]
 * @returns {number} residual [-]
 */
export function taitelDuklerResidual(h, X2, Y, n = 0.2, m = 0.2) {
  const g = stratifiedGeometry(h);
  return X2 * (g.uL * g.DL) ** -n * g.uL * g.uL * g.SL / g.AL
    - (g.uG * g.DG) ** -m * g.uG * g.uG * (g.SG / g.AG + g.Si / g.AL + g.Si / g.AG) - 4 * Y;
}

/**
 * Equilibrium stratified liquid level from X and Y. When several roots exist (upward inclination) the
 * lowest level is returned, which is the physically observed thin-film solution (Barnea 1987).
 * @param {object} o
 * @param {number} o.X2 Lockhart–Martinelli parameter squared [-], > 0
 * @param {number} [o.Y=0] inclination parameter [-] (see taitelDuklerResidual)
 * @param {number} [o.n=0.2] liquid friction exponent [-]
 * @param {number} [o.m=0.2] gas friction exponent [-]
 * @returns {number} hL/D [-]
 */
export function equilibriumLevel({ X2, Y = 0, n = 0.2, m = 0.2 }) {
  checkArg('X2', X2, { min: 0, exclusiveMin: true });
  checkArg('Y', Y);
  const f = (h) => taitelDuklerResidual(h, X2, Y, n, m);
  const N = 400;
  let hPrev = 1e-7;
  let fPrev = f(hPrev);
  if (!(fPrev > 0)) return hPrev;
  for (let k = 1; k < N; k++) {
    const h = 0.5 * (1 - Math.cos(Math.PI * k / N));
    const fh = f(h);
    if (fh <= 0) return brent(f, hPrev, h, { tol: 1e-14 });
    hPrev = h; fPrev = fh;
  }
  return brent(f, hPrev, 1 - 1e-9, { tol: 1e-14 });
}

/**
 * Dimensionless liquid level hL/D that gives a stated stratified holdup (inverse of the circular-segment area).
 * @param {number} holdup liquid area fraction [-], 0 < holdup < 1
 * @returns {number} hL/D [-]
 */
export function levelFromHoldup(holdup) {
  checkArg('holdup', holdup, { min: 0, max: 1 });
  if (holdup <= 0) return 0;
  if (holdup >= 1) return 1;
  return brent((h) => stratifiedGeometry(h).holdup - holdup, 0, 1, { tol: 1e-14 });
}

/**
 * Wall and interfacial shear of a stratified gas–liquid layer and the residual of the combined (pressure
 * eliminated) momentum balance of the two-fluid model,
 *   R = −τwL SL/AL + τwG SG/AG + τi Si (1/AL + 1/AG) − (ρL − ρG) g sinθ      [Pa/m]
 * (R = 0 for steady, fully developed stratified flow). Shear stresses are signed with the phase velocity:
 * τw = f ρ v|v|/2, τi = fi ρG (vG − vL)|vG − vL|/2, Fanning factors from the Colebrook/laminar law on the
 * hydraulic diameters DL = 4AL/SL (open channel) and DG = 4AG/(SG + Si) (closed duct).
 * Interfacial friction: 'smooth' fi = fG (Taitel & Dukler 1976); 'constant' fi = o.fi (0.0142 was used for
 * wavy flow by Shoham & Taitel 1984, AIChE J. 30(3), 377–385, after Cohen & Hanratty 1968); 'andritsos'
 * fi/fG = 1 + 15 √(hL/D) (vsg/vsg,t − 1) for vsg > vsg,t = 5 √(ρG0/ρG) m/s (Andritsos & Hanratty 1987,
 * AIChE J. 33(3), 444–454; ρG0 = gas density at atmospheric pressure).
 * Reference for the balance: Taitel & Dukler (1976), Eqs. (1)–(3).
 * @param {object} o
 * @param {number} o.h liquid level hL/D [-], 0 < h < 1
 * @param {number} o.vL in-situ liquid velocity [m/s] (signed)
 * @param {number} o.vG in-situ gas velocity [m/s] (signed)
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {'smooth'|'constant'|'andritsos'} [o.interfacial='smooth']
 * @param {number} [o.fi=0.0142] constant interfacial Fanning factor [-] ('constant')
 * @param {number} [o.rhoG0=1.2] gas density at atmospheric pressure [kg/m³] ('andritsos')
 * @returns {{residual: number, tauWL: number, tauWG: number, tauI: number, fL: number, fG: number, fI: number,
 *   holdup: number, gradFriction: number, gradGravity: number, AL: number, AG: number, SL: number, SG: number,
 *   Si: number}} residual [Pa/m]; shear stresses [Pa]; Fanning factors [-]; holdup [-]; gradFriction =
 *   (τwL SL + τwG SG)/A and gradGravity = ρmix g sinθ, pressure loss per metre [Pa/m]; areas [m²]; perimeters [m]
 */
export function stratifiedForces({ h, vL, vG, D, theta = 0, rhoL, rhoG, muL, muG, roughness = 0,
  interfacial = 'smooth', fi = 0.0142, rhoG0 = 1.2 }) {
  checkArg('h', h, { min: 0, max: 1, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  const g = stratifiedGeometry(h);
  const D2 = D * D;
  const AL = g.AL * D2;
  const AG = g.AG * D2;
  const SL = g.SL * D;
  const SG = g.SG * D;
  const Si = g.Si * D;
  const dL = g.DL * D;
  const dG = g.DG * D;
  const fL = darcyFrictionFactor(rhoL * Math.abs(vL) * dL / muL, Math.min(0.05, roughness / dL)) / 4;
  const fG = darcyFrictionFactor(rhoG * Math.abs(vG) * dG / muG, Math.min(0.05, roughness / dG)) / 4;
  let fI;
  if (interfacial === 'smooth') fI = fG;
  else if (interfacial === 'constant') fI = checkArg('fi', fi, { min: 0 });
  else if (interfacial === 'andritsos') {
    const vsg = Math.abs(vG) * (1 - g.holdup);
    const vsgT = 5 * Math.sqrt(rhoG0 / rhoG);
    fI = fG * (vsg > vsgT ? 1 + 15 * Math.sqrt(h) * (vsg / vsgT - 1) : 1);
  } else throw new RangeError(`interfacial = ${interfacial} is not a known interfacial-friction option`);
  const tauWL = fL * rhoL * vL * Math.abs(vL) / 2;
  const tauWG = fG * rhoG * vG * Math.abs(vG) / 2;
  const dv = vG - vL;
  const tauI = fI * rhoG * dv * Math.abs(dv) / 2;
  const sinT = Math.sin(theta);
  const residual = -tauWL * SL / AL + tauWG * SG / AG + tauI * Si * (1 / AL + 1 / AG) - (rhoL - rhoG) * G * sinT;
  const A = QUARTER_PI * D2;
  return {
    residual, tauWL, tauWG, tauI, fL, fG, fI, holdup: g.holdup,
    gradFriction: (tauWL * SL + tauWG * SG) / A, gradGravity: (g.holdup * rhoL + (1 - g.holdup) * rhoG) * G * sinT,
    AL, AG, SL, SG, Si,
  };
}

function fanningCoefficients(Re) {
  // Smooth-pipe Blasius-type law used by Taitel & Dukler: f = C Re^(−n).
  return Re > 2000 ? { C: 0.046, n: 0.2 } : { C: 16, n: 1 };
}

/**
 * Dimensionless groups of the Taitel–Dukler map.
 *   X² = |dP/dx|SL / |dP/dx|SG,  Y = −(ρL − ρG) g sinθ / |dP/dx|SG,
 *   F = √(ρG/(ρL − ρG)) vsg / √(D g cosθ),  K² = F² ReSL,  T² = |dP/dx|SL / ((ρL − ρG) g cosθ)
 * with superficial gradients from f = C Re^(−n): C = 0.046, n = 0.2 turbulent; C = 16, n = 1 laminar
 * (Fanning), the regime of each phase being decided on its superficial Reynolds number (2000).
 * Reference: Taitel & Dukler (1976), Eqs. (7)–(10), (24)–(31).
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @returns {{X2: number, X: number, Y: number, F: number, K: number, T: number, n: number, m: number,
 *   gradSL: number, gradSG: number, ReSL: number, ReSG: number}} groups [-]; superficial frictional gradients
 *   [Pa/m]; F, K, T are Infinity when cosθ <= 0
 */
export function taitelDuklerParameters({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, max: rhoL, exclusiveMin: true });
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('muG', muG, { min: 0, exclusiveMin: true });
  const ReSL = rhoL * vsl * D / muL;
  const ReSG = rhoG * vsg * D / muG;
  const cl = fanningCoefficients(ReSL);
  const cg = fanningCoefficients(ReSG);
  const gradSL = 4 * cl.C * ReSL ** -cl.n / D * rhoL * vsl * vsl / 2;
  const gradSG = 4 * cg.C * ReSG ** -cg.n / D * rhoG * vsg * vsg / 2;
  const dRho = rhoL - rhoG;
  const cosT = Math.cos(theta);
  const X2 = gradSL / gradSG;
  const Y = -dRho * G * Math.sin(theta) / gradSG;
  const ok = cosT > 1e-12;
  const F = ok ? Math.sqrt(rhoG / dRho) * vsg / Math.sqrt(D * G * cosT) : Infinity;
  const K = ok ? F * Math.sqrt(ReSL) : Infinity;
  const T = ok ? Math.sqrt(gradSL / (dRho * G * cosT)) : Infinity;
  return { X2, X: Math.sqrt(X2), Y, F, K, T, n: cl.n, m: cg.n, gradSL, gradSG, ReSL, ReSG };
}

/**
 * Taitel & Dukler (1976) flow pattern for horizontal and slightly inclined pipes.
 * Steps: (1) equilibrium stratified level from the combined momentum balance; (2) stratified is unstable when
 * F² [ũG² (dÃL/dh̃L) / ((1 − h̃L)² ÃG)] ≥ 1 (Kelvin–Helmholtz with the (1 − h̃L) finite-wave factor);
 * (3) if unstable: annular when h̃L < annularLevel, else dispersed bubble when
 * T² ≥ 8 ÃG / (S̃i ũL² (ũL D̃L)^(−n)), else intermittent; (4) if stable: wavy when
 * K ≥ 2 / (√ũL ũG √s) with Jeffreys' sheltering coefficient s = 0.01, else smooth. For downward flow the
 * film is also taken as wavy when its Froude number vL/√(g hL) exceeds 1.5 (Barnea, Shoham & Taitel 1982,
 * Chem. Eng. Sci. 37(5), 735–740).
 * Reference: Taitel & Dukler (1976), AIChE J. 22(1), 47–55. Validity: |θ| ≲ 10°, fully developed flow,
 * low-viscosity liquids; larger angles are evaluated but flagged by flowRegime().
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.annularLevel=0.5] h̃L below which unstable stratified flow becomes annular [-]
 *   (0.5 Taitel & Dukler; 0.35 Barnea et al. 1982)
 * @param {number} [o.sheltering=0.01] Jeffreys sheltering coefficient s [-]
 * @returns {{regime: 'stratified-smooth'|'stratified-wavy'|'intermittent'|'annular'|'dispersed-bubble',
 *   hL: number, holdupStratified: number, X: number, Y: number, F: number, K: number, T: number,
 *   stratifiedStabilityRatio: number, dispersedRatio: number, wavyRatio: number, residual: number}}
 *   hL = equilibrium hL/D [-]; holdupStratified [-]; ratios are (left side)/(right side) of each criterion
 *   (transition at 1) [-]; residual of the momentum balance at hL [-]
 */
export function taitelDukler({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, annularLevel = 0.5, sheltering = 0.01 }) {
  checkArg('annularLevel', annularLevel, { min: 0, max: 1, exclusiveMin: true });
  checkArg('sheltering', sheltering, { min: 0, exclusiveMin: true });
  const p = taitelDuklerParameters({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG });
  const hL = equilibriumLevel(p);
  const g = stratifiedGeometry(hL);
  const C2 = 1 - hL;
  const stratifiedStabilityRatio = p.F * p.F * g.uG * g.uG * g.dAL / (C2 * C2 * g.AG);
  const dispersedRatio = p.T * p.T * g.Si * g.uL * g.uL * (g.uL * g.DL) ** -p.n / (8 * g.AG);
  let wavyRatio = p.K * Math.sqrt(g.uL) * g.uG * Math.sqrt(sheltering) / 2;
  if (theta < 0) {
    const frFilm = vsl * g.uL / Math.sqrt(G * hL * D);
    wavyRatio = Math.max(wavyRatio, frFilm / 1.5);
  }
  let regime;
  if (stratifiedStabilityRatio >= 1) {
    if (hL < annularLevel) regime = 'annular';
    else regime = dispersedRatio >= 1 ? 'dispersed-bubble' : 'intermittent';
  } else regime = wavyRatio >= 1 ? 'stratified-wavy' : 'stratified-smooth';
  return {
    regime, hL, holdupStratified: g.holdup, X: p.X, Y: p.Y, F: p.F, K: p.K, T: p.T,
    stratifiedStabilityRatio, dispersedRatio, wavyRatio, residual: taitelDuklerResidual(hL, p.X2, p.Y, p.n, p.m),
  };
}

/**
 * Taitel, Barnea & Dukler (1980) flow pattern for upward vertical flow.
 *   annular            vsg ρG^½ / [σ g (ρL − ρG)]^¼ ≥ 3.1 (droplet-suspension, Turner type)
 *   dispersed bubble   vm ≥ 4.0 D^0.429 (σ/ρL)^0.089 / νL^0.072 · [g (ρL − ρG)/ρL]^0.446  and  α_ns ≤ 0.52
 *   bubble             vsl ≥ 3.0 vsg − 1.15 [g (ρL − ρG) σ/ρL²]^¼ (α < 0.25), only when
 *                      [ρL² g D² / ((ρL − ρG) σ)]^¼ ≥ 4.36
 *   churn              slug flow that has not developed: L/D < lE/D = 40.6 (vm/√(gD) + 0.22)
 *   slug               otherwise
 * Reference: Taitel, Barnea & Dukler (1980), AIChE J. 26(3), 345–354, Eqs. (5), (7), (13), (14), (26), (30).
 * Validity: upward vertical (use for θ ≳ 70°), low-viscosity liquid.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.sigma surface tension [N/m]
 * @param {number} [o.LoverD=Infinity] distance from the pipe entrance in diameters [-]; Infinity = developed
 *   flow (no churn region)
 * @returns {{regime: 'bubble'|'dispersed-bubble'|'slug'|'churn'|'annular', annularRatio: number,
 *   dispersedRatio: number, bubbleMargin: number, bubblePossible: boolean, entranceLengthOverD: number}}
 *   annularRatio = vsg/vsg,annular [-]; dispersedRatio = vm/vm,dispersed [-]; bubbleMargin = vsl − vsl,transition
 *   [m/s]; entranceLengthOverD = lE/D [-]
 */
export function taitelBarneaDuklerVertical({ vsl, vsg, D, rhoL, rhoG, muL, sigma, LoverD = Infinity }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, max: rhoL, exclusiveMin: true });
  checkArg('muL', muL, { min: 0, exclusiveMin: true });
  checkArg('sigma', sigma, { min: 0, exclusiveMin: true });
  if (!(LoverD > 0)) throw new RangeError(`LoverD = ${LoverD} must be positive`);
  const dRho = rhoL - rhoG;
  const vm = vsl + vsg;
  const annularRatio = vsg * Math.sqrt(rhoG) / (3.1 * (sigma * G * dRho) ** 0.25);
  const vmDispersed = 4.0 * D ** 0.429 * (sigma / rhoL) ** 0.089 / (muL / rhoL) ** 0.072 * (G * dRho / rhoL) ** 0.446;
  const dispersedRatio = vm / vmDispersed;
  const bubblePossible = (rhoL * rhoL * G * D * D / (dRho * sigma)) ** 0.25 >= 4.36;
  const bubbleMargin = vsl - (3.0 * vsg - 1.15 * (G * dRho * sigma / (rhoL * rhoL)) ** 0.25);
  const entranceLengthOverD = 40.6 * (vm / Math.sqrt(G * D) + 0.22);
  let regime;
  if (annularRatio >= 1) regime = 'annular';
  else if (dispersedRatio >= 1 && vsg / vm <= 0.52) regime = 'dispersed-bubble';
  else if (bubblePossible && bubbleMargin >= 0 && dispersedRatio < 1) regime = 'bubble';
  else regime = LoverD < entranceLengthOverD ? 'churn' : 'slug';
  return { regime, annularRatio, dispersedRatio, bubbleMargin, bubblePossible, entranceLengthOverD };
}

/**
 * Beggs & Brill (1973) horizontal flow-pattern boundaries and pattern.
 *   L1 = 316 λ^0.302, L2 = 0.0009252 λ^(−2.4684), L3 = 0.10 λ^(−1.4516), L4 = 0.5 λ^(−6.738), Fr = vm²/(gD)
 *   segregated   λ < 0.01 and Fr < L1, or λ ≥ 0.01 and Fr < L2
 *   transition   λ ≥ 0.01 and L2 ≤ Fr ≤ L3
 *   intermittent 0.01 ≤ λ < 0.4 and L3 < Fr ≤ L1, or λ ≥ 0.4 and L3 < Fr ≤ L4
 *   distributed  λ < 0.4 and Fr ≥ L1, or λ ≥ 0.4 and Fr > L4
 * Reference: Beggs & Brill (1973), J. Pet. Technol. 25(5), 607–617, with the revised boundaries of Brill &
 * Beggs, Two-Phase Flow in Pipes (1991). Validity: 25–38 mm air–water data; widely extrapolated.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vsg superficial gas velocity [m/s], >= 0
 * @param {number} o.D pipe inner diameter [m]
 * @returns {{pattern: 'segregated'|'transition'|'intermittent'|'distributed', lambdaL: number, Fr: number,
 *   L1: number, L2: number, L3: number, L4: number}} all [-]
 */
export function beggsBrillPattern({ vsl, vsg, D }) {
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vsg', vsg, { min: 0 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  const vm = vsl + vsg;
  checkArg('vsl + vsg', vm, { min: 0, exclusiveMin: true });
  const lam = vsl / vm;
  const Fr = vm * vm / (G * D);
  const L1 = 316 * lam ** 0.302;
  const L2 = 0.0009252 * lam ** -2.4684;
  const L3 = 0.10 * lam ** -1.4516;
  const L4 = 0.5 * lam ** -6.738;
  let pattern;
  if ((lam < 0.01 && Fr < L1) || (lam >= 0.01 && Fr < L2)) pattern = 'segregated';
  else if (lam >= 0.01 && Fr >= L2 && Fr <= L3) pattern = 'transition';
  else if ((lam >= 0.01 && lam < 0.4 && Fr > L3 && Fr <= L1) || (lam >= 0.4 && Fr > L3 && Fr <= L4)) pattern = 'intermittent';
  else pattern = 'distributed';
  return { pattern, lambdaL: lam, Fr, L1, L2, L3, L4 };
}

const DEG = Math.PI / 180;

/**
 * Flow pattern for any inclination by combining the two mechanistic maps:
 *   |θ| ≤ nearHorizontalDeg (10°)      Taitel & Dukler (1976) at the actual angle
 *   θ ≥ nearVerticalDeg (70°)          Taitel, Barnea & Dukler (1980)
 *   in between, upward                 both maps are evaluated; the Taitel–Dukler result is used below the
 *                                      mid-angle and the vertical result above it (`interpolated: true`)
 *   θ < −nearHorizontalDeg             Taitel & Dukler momentum balance at the actual angle (`extrapolated`)
 * This is a documented engineering blend, not Barnea's (1987) unified model.
 * @param {object} o as taitelDukler plus
 * @param {number} o.sigma surface tension [N/m]
 * @param {number} [o.LoverD] entrance length in diameters for the churn criterion [-]
 * @param {number} [o.nearHorizontalDeg=10] upper angle of the near-horizontal range [deg]
 * @param {number} [o.nearVerticalDeg=70] lower angle of the near-vertical range [deg]
 * @returns {{regime: string, model: 'taitel-dukler'|'taitel-barnea-dukler', interpolated: boolean,
 *   extrapolated: boolean, horizontal: object|null, vertical: object|null}} regime name; the detailed results
 *   of the maps that were evaluated
 */
export function flowRegime({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, sigma, LoverD, annularLevel, sheltering,
  nearHorizontalDeg = 10, nearVerticalDeg = 70 }) {
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  const deg = theta / DEG;
  const td = () => taitelDukler({ vsl, vsg, D, theta, rhoL, rhoG, muL, muG, annularLevel, sheltering });
  const tbd = () => taitelBarneaDuklerVertical({ vsl, vsg, D, rhoL, rhoG, muL, sigma, LoverD });
  if (Math.abs(deg) <= nearHorizontalDeg) {
    const h = td();
    return { regime: h.regime, model: 'taitel-dukler', interpolated: false, extrapolated: false, horizontal: h, vertical: null };
  }
  if (deg >= nearVerticalDeg) {
    const v = tbd();
    return { regime: v.regime, model: 'taitel-barnea-dukler', interpolated: false, extrapolated: false, horizontal: null, vertical: v };
  }
  if (deg < 0) {
    const h = td();
    return { regime: h.regime, model: 'taitel-dukler', interpolated: false, extrapolated: true, horizontal: h, vertical: null };
  }
  const h = td();
  const v = tbd();
  const useVertical = deg >= 0.5 * (nearHorizontalDeg + nearVerticalDeg);
  return {
    regime: useVertical ? v.regime : h.regime, model: useVertical ? 'taitel-barnea-dukler' : 'taitel-dukler',
    interpolated: true, extrapolated: false, horizontal: h, vertical: v,
  };
}

/**
 * Flow-regime map on a logarithmic (vsg, vsl) grid for plotting, with the operating point classified.
 * @param {object} o fluid and pipe arguments of the selected model (D, theta, rhoL, rhoG, muL, muG, sigma, …) plus
 * @param {'unified'|'taitel-dukler'|'vertical'|'beggs-brill'} [o.model='unified']
 * @param {[number, number]} [o.vsgRange=[0.01, 100]] gas superficial-velocity range [m/s]
 * @param {[number, number]} [o.vslRange=[0.001, 10]] liquid superficial-velocity range [m/s]
 * @param {number} [o.nx=60] grid points along vsg [-]
 * @param {number} [o.ny=60] grid points along vsl [-]
 * @param {{vsg: number, vsl: number}} [o.point] operating point [m/s]
 * @returns {{xs: number[], ys: number[], z: number[][], legend: string[], regimes: string[][],
 *   point: {vsg: number, vsl: number, regime: string}|null}} xs = vsg [m/s], ys = vsl [m/s], z[iy][ix] = index
 *   into legend [-], regimes[iy][ix] = name
 */
export function regimeMap({ model = 'unified', vsgRange = [0.01, 100], vslRange = [0.001, 10], nx = 60, ny = 60, point, ...fluid }) {
  checkArg('nx', nx, { min: 2, max: 400 });
  checkArg('ny', ny, { min: 2, max: 400 });
  checkArg('vsgRange[0]', vsgRange[0], { min: 0, exclusiveMin: true });
  checkArg('vsgRange[1]', vsgRange[1], { min: vsgRange[0], exclusiveMin: true });
  checkArg('vslRange[0]', vslRange[0], { min: 0, exclusiveMin: true });
  checkArg('vslRange[1]', vslRange[1], { min: vslRange[0], exclusiveMin: true });
  const classify = (vsg, vsl) => {
    switch (model) {
      case 'unified': return flowRegime({ ...fluid, vsl, vsg }).regime;
      case 'taitel-dukler': return taitelDukler({ ...fluid, vsl, vsg }).regime;
      case 'vertical': return taitelBarneaDuklerVertical({ ...fluid, vsl, vsg }).regime;
      case 'beggs-brill': return beggsBrillPattern({ D: fluid.D, vsl, vsg }).pattern;
      default: throw new RangeError(`model = ${model} is not a known flow-regime model`);
    }
  };
  const axis = (r, n) => Array.from({ length: n }, (_, i) => r[0] * (r[1] / r[0]) ** (i / (n - 1)));
  const xs = axis(vsgRange, nx);
  const ys = axis(vslRange, ny);
  const legend = [];
  const z = [];
  const regimes = [];
  for (let iy = 0; iy < ny; iy++) {
    const row = new Array(nx);
    const names = new Array(nx);
    for (let ix = 0; ix < nx; ix++) {
      const r = classify(xs[ix], ys[iy]);
      let k = legend.indexOf(r);
      if (k < 0) { legend.push(r); k = legend.length - 1; }
      row[ix] = k; names[ix] = r;
    }
    z.push(row); regimes.push(names);
  }
  const pt = point ? { vsg: point.vsg, vsl: point.vsl, regime: classify(point.vsg, point.vsl) } : null;
  return { xs, ys, z, legend, regimes, point: pt };
}
