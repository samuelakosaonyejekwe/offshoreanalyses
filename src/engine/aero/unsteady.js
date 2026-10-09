/**
 * Unsteady rotor aerodynamics: dynamic inflow, dynamic stall, tower shadow, shear/veer, rotational sampling and a
 * time-marching unsteady BEM for one rigid rotor.
 */

import { bemSolve, bemElement, towerPotentialFactor } from './bem.js';
import { polarLinearFit } from './polar.js';

const D2R = Math.PI / 180;
const TWO_PI = 2 * Math.PI;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/**
 * Time constants of the Oye dynamic-inflow model:
 *   tau1 = 1.1 / (1 - 1.3 a) * R / V0,   tau2 = (0.39 - 0.26 (r/R)^2) tau1,   a limited to 0.5.
 * Reference: Oye, S. (1990); Hansen, M.O.L., Aerodynamics of Wind Turbines, 2nd ed., Eqs. 9.11-9.12.
 * @param {object} o
 * @param {number} o.a axial induction factor [-]
 * @param {number} o.R rotor radius [m]
 * @param {number} o.V0 free-stream speed [m/s] (> 0)
 * @param {number} o.r local radius [m]
 * @returns {{tau1:number, tau2:number}} time constants [s]
 */
export function oyeTimeConstants(o) {
  need(o && o.V0 > 0, 'V0', 'must be > 0');
  need(o.R > 0 && o.r >= 0 && o.r <= o.R * (1 + 1e-9), 'r', 'must be in [0, R]');
  const a = Math.min(0.5, Math.max(0, o.a));
  const tau1 = 1.1 / (1 - 1.3 * a) * o.R / o.V0;
  return { tau1, tau2: (0.39 - 0.26 * (o.r / o.R) ** 2) * tau1 };
}

/**
 * One step of the Oye dynamic-inflow filter for an induced velocity:
 *   y + tau1 dy/dt = x + k tau1 dx/dt,    W + tau2 dW/dt = y,    k = 0.6
 * integrated analytically over the step (Hansen, Aerodynamics of Wind Turbines, Eqs. 9.13-9.15):
 *   H = x_i + k tau1 (x_i - x_{i-1})/dt,  y_i = H + (y_{i-1} - H) exp(-dt/tau1),  W_i = y_i + (W_{i-1} - y_i) exp(-dt/tau2).
 * @param {{xPrev:number, y:number, W:number}} state filter state, mutated: previous quasi-steady value, intermediate
 *   and filtered induced velocity [m/s]; initialise all three with the quasi-steady value
 * @param {object} o
 * @param {number} o.x quasi-steady induced velocity at the new time [m/s]
 * @param {number} o.dt time step [s]
 * @param {number} o.tau1 [s]
 * @param {number} o.tau2 [s]
 * @param {number} [o.k=0.6] lead constant [-]
 * @returns {object} the state (W is the filtered induced velocity [m/s])
 */
export function oyeDynamicInflowStep(state, o) {
  const k = o.k ?? 0.6;
  const H = o.x + k * o.tau1 * (o.x - state.xPrev) / o.dt;
  state.y = H + (state.y - H) * Math.exp(-o.dt / o.tau1);
  state.W = state.y + (state.W - state.y) * Math.exp(-o.dt / o.tau2);
  state.xPrev = o.x;
  return state;
}

/**
 * Rate of change of the mean axial induced velocity from the Pitt-Peters dynamic-inflow model reduced to its uniform
 * (axial) state: apparent mass of an impermeable disc m_a = (8/3) rho R^3,
 *   m_a dv/dt + 2 rho pi R^2 v |V - v| = T.
 * Reference: Pitt, D.M. and Peters, D.A. (1981), "Theoretical prediction of dynamic-inflow derivatives", Vertica 5.
 * The steady state reproduces momentum theory, T = 2 rho A v (V - v).
 * @param {object} o
 * @param {number} o.v mean axial induced velocity at the disc [m/s]
 * @param {number} o.thrust rotor thrust [N]
 * @param {number} o.V free-stream speed [m/s]
 * @param {number} o.R rotor radius [m]
 * @param {number} [o.rho=1.225] density [kg/m3]
 * @returns {number} dv/dt [m/s2]
 */
export function pittPetersRate(o) {
  need(o && o.R > 0, 'R', 'must be > 0');
  const rho = o.rho ?? 1.225;
  return (o.thrust - 2 * rho * Math.PI * o.R * o.R * o.v * Math.abs(o.V - o.v)) / ((8 / 3) * rho * o.R ** 3);
}

/**
 * Kirchhoff/Helmholtz trailing-edge separation description of a static polar:
 *   Cl,st = m (alpha - alpha0) ((1 + sqrt(f))/2)^2   ->   f(alpha) = (2 sqrt(Cl,st / (m (alpha - alpha0))) - 1)^2
 * (f = 1 attached, f = 0 fully separated, reached when Cl,st falls to a quarter of the linear value), and the
 * fully-separated lift Cl,fs = (Cl,st - m (alpha - alpha0) f)/(1 - f), with Cl,fs = Cl,st/2 where f = 1.
 * Reference: Hansen, M.H., Gaunaa, M. and Madsen, H.Aa. (2004), "A Beddoes-Leishman type dynamic stall model in
 * state-space and indicial formulations", Riso-R-1354(EN), section 2.
 * @param {object} polar static polar (polar.js), extrapolated to +-180 deg
 * @param {object} [o]
 * @param {number} [o.alpha0] zero-lift angle [rad]; default fitted
 * @param {number} [o.slope] attached lift slope [1/rad]; default fitted
 * @returns {{alpha0:number, slope:number, cd0:number, clMax:number, fAt:function(number):number, clFsAt:function(number):number}}
 *   zero-lift angle [rad], slope [1/rad], drag at alpha0 [-], static maximum lift [-], evaluators for alpha in rad
 */
export function createSeparationModel(polar, o = {}) {
  const fit = polarLinearFit(polar);
  const alpha0 = o.alpha0 ?? fit.alpha0;
  const m = o.slope ?? fit.slope;
  need(m > 0, 'slope', 'must be > 0');
  const n = 1441;                         // 0.25 deg grid over +-180 deg
  const f = new Float64Array(n);
  const fs = new Float64Array(n);
  let clMax = -Infinity;
  for (let i = 0; i < n; i++) {
    const al = -Math.PI + TWO_PI * i / (n - 1);
    const cl = polar.clAt(al);
    if (Math.abs(al) < 0.6 && cl > clMax) clMax = cl;
    const lin = m * (al - alpha0);
    let fi;
    if (Math.abs(al - alpha0) < 1e-6) fi = 1;
    else {
      const x = cl / lin;
      fi = x >= 1 ? 1 : (x <= 0.25 ? 0 : (2 * Math.sqrt(x) - 1) ** 2);
    }
    // beyond 60 deg from alpha0 the flow is treated as fully separated
    if (Math.abs(al - alpha0) > Math.PI / 3) fi = 0;
    f[i] = fi;
    fs[i] = fi > 1 - 1e-6 ? 0.5 * cl : (cl - lin * fi) / (1 - fi);
  }
  const at = (arr) => (alpha) => {
    let a = alpha;
    if (a > Math.PI || a < -Math.PI) a -= TWO_PI * Math.floor((a + Math.PI) / TWO_PI);
    const u = (a + Math.PI) / TWO_PI * (n - 1);
    const i = Math.min(n - 2, Math.max(0, Math.floor(u)));
    return arr[i] + (u - i) * (arr[i + 1] - arr[i]);
  };
  return { alpha0, slope: m, cd0: polar.cdAt(alpha0), clMax, fAt: at(f), clFsAt: at(fs) };
}

/**
 * One step of Oye's dynamic-stall model (trailing-edge separation lag only):
 *   Cl = fs Cl,inv + (1 - fs) Cl,fs,   dfs/dt = (f_st(alpha) - fs)/tau,   tau = A c / W,   Cl,inv = m (alpha - alpha0),
 * integrated analytically over the step. Reference: Oye, S. (1991), "Dynamic stall simulated as time lag of
 * separation"; Hansen, Aerodynamics of Wind Turbines, section 9.3 (A about 4).
 * @param {{fs:number}} state separation state [-], mutated (initialise with sep.fAt(alpha))
 * @param {object} o
 * @param {object} o.sep separation model (createSeparationModel)
 * @param {number} o.alpha angle of attack [rad]
 * @param {number} o.W relative speed [m/s]
 * @param {number} o.chord chord [m]
 * @param {number} o.dt time step [s]
 * @param {number} [o.A=4] time-constant factor [-]
 * @returns {object} state with `cl` [-]
 */
export function oyeDynamicStallStep(state, o) {
  const tau = (o.A ?? 4) * o.chord / Math.max(o.W, 1e-6);
  const fst = o.sep.fAt(o.alpha);
  state.fs = fst + (state.fs - fst) * Math.exp(-o.dt / tau);
  state.cl = state.fs * o.sep.slope * (o.alpha - o.sep.alpha0) + (1 - state.fs) * o.sep.clFsAt(o.alpha);
  return state;
}

/**
 * Initial state of the Beddoes-Leishman model in steady flow at angle alpha.
 * @param {object} sep separation model
 * @param {number} alpha angle of attack [rad]
 * @param {object} [c] constants (A1, A2 as in beddoesLeishmanStep)
 * @returns {{x1:number,x2:number,x3:number,x4:number,xv:number,cvPrev:number,tauV:number,alphaPrev:number,cl:number,cd:number}} state
 */
export function beddoesLeishmanInit(sep, alpha, c = {}) {
  const A1 = c.A1 ?? 0.165;
  const A2 = c.A2 ?? 0.335;
  return { x1: A1 * alpha, x2: A2 * alpha, x3: sep.slope * (alpha - sep.alpha0), x4: sep.fAt(alpha), xv: 0, cvPrev: 0, tauV: 0, alphaPrev: alpha, cl: 0, cd: 0 };
}

/**
 * One step of a Beddoes-Leishman type dynamic-stall model in the incompressible state-space form of Hansen, Gaunaa
 * and Madsen (2004, Riso-R-1354), integrated with exact exponentials over the step for frozen U:
 *   states x1, x2: attached-flow (Theodorsen) lag with the Jones two-term approximation of the Wagner function,
 *                  phi(s) = 1 - A1 exp(-b1 s) - A2 exp(-b2 s), A1 = 0.165, A2 = 0.335, b1 = 0.0455, b2 = 0.3,
 *                  s = 2 U t / c;  alpha_E = alpha (1 - A1 - A2) + x1 + x2;
 *   state x3:      leading-edge pressure lag of the unsteady lift, time constant Tp * c/(2U);
 *   state x4:      trailing-edge separation point lag, time constant Tf * c/(2U), driven by f_st(x3/m + alpha0);
 *   Cl = m (alpha_E - alpha0) x4 + Cl,fs(alpha_E) (1 - x4) + pi (c/2U) dalpha/dt            (last term: added mass)
 *   Cd = Cd,st(alpha_E) + (alpha - alpha_E) Cl + (Cd,st(alpha_E) - Cd0) [ (sqrt(f_st(alpha_E)) - sqrt(x4))/2 - (f_st(alpha_E) - x4)/4 ].
 * Optional fifth state (leading-edge vortex lift, simplified from Leishman and Beddoes 1989): while the lagged lift
 * x3 exceeds the critical value `cn1` the vortex-lift increment Cv = m (alpha_E - alpha0) (1 - ((1 + sqrt(x4))/2)^2)
 * feeds  xv_n = xv_{n-1} exp(-ds/Tv) + (Cv_n - Cv_{n-1}) exp(-ds/(2 Tv))  for a travel time of Tvl semi-chords,
 * after which xv only decays. Tv and Tvl are airfoil-dependent fitted constants (Leishman and Beddoes quote Tv
 * about 6 and Tvl about 7-11 for the NACA 0012); they are arguments and the vortex term is off by default.
 * Limits: incompressible, plunge and pitch not distinguished (alpha taken at the three-quarter chord), no
 * leading-edge-stall criterion, constants not tuned to a particular airfoil.
 * @param {object} state state from beddoesLeishmanInit, mutated
 * @param {object} o
 * @param {object} o.sep separation model (createSeparationModel)
 * @param {object} o.polar static polar (for Cd)
 * @param {number} o.alpha angle of attack at the new time [rad]
 * @param {number} o.W relative speed [m/s]
 * @param {number} o.chord chord [m]
 * @param {number} o.dt time step [s]
 * @param {number} [o.A1=0.165] @param {number} [o.A2=0.335] @param {number} [o.b1=0.0455] @param {number} [o.b2=0.3]
 * @param {number} [o.Tp=1.5] pressure-lag constant [semi-chords] (typical 0.8 - 2)
 * @param {number} [o.Tf=6] separation-lag constant [semi-chords] (typical 3 - 10)
 * @param {boolean} [o.vortexLift=false]
 * @param {number} [o.Tv=6] vortex decay constant [semi-chords]
 * @param {number} [o.Tvl=11] vortex travel time [semi-chords]
 * @param {number} [o.cn1] critical lift for vortex onset [-]; default the static maximum lift
 * @returns {object} state with `cl`, `cd` [-]
 */
export function beddoesLeishmanStep(state, o) {
  const A1 = o.A1 ?? 0.165; const A2 = o.A2 ?? 0.335; const b1 = o.b1 ?? 0.0455; const b2 = o.b2 ?? 0.3;
  const sep = o.sep;
  const Tu = 0.5 * o.chord / Math.max(o.W, 1e-6);
  const ds = o.dt / Tu;
  const aDot = (o.alpha - state.alphaPrev) / o.dt;
  const aMid = 0.5 * (o.alpha + state.alphaPrev);
  const e1 = Math.exp(-b1 * ds); const e2 = Math.exp(-b2 * ds);
  state.x1 = state.x1 * e1 + A1 * aMid * (1 - e1);
  state.x2 = state.x2 * e2 + A2 * aMid * (1 - e2);
  const aE = o.alpha * (1 - A1 - A2) + state.x1 + state.x2;
  const clCirc = sep.slope * (aE - sep.alpha0);
  const clAdd = Math.PI * Tu * aDot;
  const ep = Math.exp(-ds / (o.Tp ?? 1.5));
  state.x3 = state.x3 * ep + (clCirc + clAdd) * (1 - ep);
  const fLag = sep.fAt(state.x3 / sep.slope + sep.alpha0);
  const ef = Math.exp(-ds / (o.Tf ?? 6));
  state.x4 = state.x4 * ef + fLag * (1 - ef);
  let cl = clCirc * state.x4 + sep.clFsAt(aE) * (1 - state.x4) + clAdd;
  if (o.vortexLift) {
    const Tv = o.Tv ?? 6; const Tvl = o.Tvl ?? 11;
    const cv = clCirc * (1 - (0.5 * (1 + Math.sqrt(Math.max(0, state.x4)))) ** 2);
    if (Math.abs(state.x3) > (o.cn1 ?? sep.clMax)) state.tauV += ds; else state.tauV = 0;
    state.xv *= Math.exp(-ds / Tv);
    if (state.tauV > 0 && state.tauV < Tvl) state.xv += (cv - state.cvPrev) * Math.exp(-0.5 * ds / Tv);
    state.cvPrev = cv;
    cl += state.xv;
  }
  const fE = sep.fAt(aE);
  const cdSt = o.polar.cdAt(aE);
  state.cd = cdSt + (o.alpha - aE) * cl + (cdSt - sep.cd0) * (0.5 * (Math.sqrt(fE) - Math.sqrt(Math.max(0, state.x4))) - 0.25 * (fE - state.x4));
  state.cl = cl;
  state.alphaE = aE;
  state.alphaPrev = o.alpha;
  return state;
}

/**
 * Tower shadow of an upwind rotor: streamwise speed factor from the potential flow about the tower (see
 * bem.js towerPotentialFactor).
 * @param {object} o
 * @param {number} o.x streamwise distance of the blade element from the tower axis [m]
 * @param {number} o.y lateral distance from the tower axis [m]
 * @param {number} o.towerDiameter local tower diameter [m]
 * @returns {number} u/U [-]
 */
export function towerShadowPotential(o) {
  return towerPotentialFactor(o.x, o.y, 0.5 * o.towerDiameter);
}

/**
 * Wind speed and direction at height z from hub-height values: power law V = Vhub (z/zhub)^alpha (IEC 61400-1
 * normal wind profile) or logarithmic law V = Vhub ln(z/z0)/ln(zhub/z0) when `z0` is given; linear veer with height.
 * @param {object} o
 * @param {number} o.Vhub hub-height speed [m/s]
 * @param {number} o.z height [m] (> 0)
 * @param {number} o.hubHeight [m]
 * @param {number} [o.shearExp=0.2] power-law exponent [-] (IEC: 0.2 onshore, 0.14 offshore)
 * @param {number} [o.z0] roughness length [m]; selects the logarithmic law
 * @param {number} [o.veerDegPerM=0] direction change per metre of height [deg/m]
 * @returns {{speed:number, directionDeg:number}} speed [m/s], direction relative to hub height [deg]
 */
export function shearedWind(o) {
  need(o && o.z > 0 && o.hubHeight > 0, 'z/hubHeight', 'must be > 0');
  let speed;
  if (o.z0 !== undefined) {
    need(o.z0 > 0 && o.z0 < o.hubHeight, 'z0', 'must be in (0, hubHeight)');
    speed = o.Vhub * Math.log(Math.max(o.z, o.z0 * 1.0001) / o.z0) / Math.log(o.hubHeight / o.z0);
  } else {
    speed = o.Vhub * Math.pow(o.z / o.hubHeight, o.shearExp ?? 0.2);
  }
  return { speed, directionDeg: (o.veerDegPerM ?? 0) * (o.z - o.hubHeight) };
}

/**
 * Discrete Fourier harmonics of one period of a uniformly sampled periodic signal:
 * y(psi) = mean + sum_n amp_n cos(n psi - phase_n).
 * @param {number[]} y samples over one period (uniform, end point excluded)
 * @param {number} [nMax=8] highest harmonic
 * @returns {{mean:number, harmonics:Array<{n:number, amp:number, phase:number}>}} amplitudes in the units of y, phase [rad]
 */
export function harmonicAmplitudes(y, nMax = 8) {
  const N = y.length;
  need(N >= 2, 'y', 'at least two samples required');
  let mean = 0;
  for (let i = 0; i < N; i++) mean += y[i];
  mean /= N;
  const harmonics = [];
  for (let n = 1; n <= Math.min(nMax, Math.floor((N - 1) / 2)); n++) {
    let c = 0; let s = 0;
    for (let i = 0; i < N; i++) { const p = TWO_PI * n * i / N; c += y[i] * Math.cos(p); s += y[i] * Math.sin(p); }
    harmonics.push({ n, amp: 2 * Math.hypot(c, s) / N, phase: Math.atan2(s, c) });
  }
  return { mean, harmonics };
}

/**
 * Rotational sampling of a steady, spatially non-uniform inflow (shear, tower shadow, tilt, yaw): per-blade loads
 * versus azimuth from the sector BEM and their harmonics. A single blade sees multiples of the rotation frequency
 * (1P, 2P, ...); summed over B equally spaced identical blades only multiples of B survive in the rotor loads (3P,
 * 6P for three blades), which is verified exactly by the phase sum.
 * Limits: quasi-steady induction per sector (no dynamic inflow/stall), rigid blades.
 * @param {object} o bemSolve options (rotor, U, omega, pitchDeg, rho, shearExp, towerShadow, yawDeg, ...)
 * @param {number} [o.nAzimuth=72] azimuth positions (multiple of the blade number)
 * @param {number} [o.nHarmonics=9]
 * @returns {{psi:number[], bladeFlapMoment:number[], bladeThrust:number[], rotorThrust:number[], rotorTorque:number[],
 *   bladeFlapHarmonics:object, rotorThrustHarmonics:object, rotorTorqueHarmonics:object, frequency1P:number}}
 *   azimuth [rad]; blade root flapwise moment [N m] and thrust [N] versus azimuth; instantaneous rotor thrust [N]
 *   and torque [N m] versus azimuth of blade 1; harmonic sets (see harmonicAmplitudes, n = multiples of 1P);
 *   rotation frequency [Hz]
 */
export function rotationalSampling(o) {
  const B = o.rotor.blades;
  let nAz = o.nAzimuth ?? 72;
  nAz = Math.ceil(nAz / B) * B;
  const s = bemSolve({ ...o, nSectors: nAz });
  const az = s.azimuthal;
  const rotorThrust = new Array(nAz).fill(0);
  const rotorTorque = new Array(nAz).fill(0);
  for (let i = 0; i < nAz; i++) {
    for (let b = 0; b < B; b++) {
      const j = (i + b * nAz / B) % nAz;
      rotorThrust[i] += az.bladeThrust[j];
      rotorTorque[i] += az.bladeTorque[j];
    }
  }
  const nH = o.nHarmonics ?? 9;
  return {
    psi: az.psi, bladeFlapMoment: az.flapMoment, bladeThrust: az.bladeThrust, rotorThrust, rotorTorque,
    bladeFlapHarmonics: harmonicAmplitudes(az.flapMoment, nH),
    rotorThrustHarmonics: harmonicAmplitudes(rotorThrust, nH),
    rotorTorqueHarmonics: harmonicAmplitudes(rotorTorque, nH),
    frequency1P: o.omega / TWO_PI,
  };
}

function seriesAt(v, k, t, fallback) {
  if (v === undefined || v === null) return fallback;
  if (typeof v === 'number') return v;
  if (typeof v === 'function') return v(t);
  return v[Math.min(k, v.length - 1)];
}

/**
 * Time-marching unsteady BEM of one rigid rotor with prescribed inflow, rotor speed and pitch histories.
 * Per blade element and time step: local wind (rotor-effective speed with optional power-law shear and potential
 * tower shadow, or a caller-supplied field) -> quasi-steady induction from the steady element equations
 * (bem.js) -> Oye dynamic-inflow filter on the axial and tangential induced velocities -> angle of attack ->
 * optional dynamic stall (Oye or Beddoes-Leishman) -> loads. Blades are treated independently (no annulus
 * averaging), yaw and skewed wake are not included; tilt and precone enter the kinematics as in bemSolve.
 * @param {object} o
 * @param {object} o.rotor rotor
 * @param {number} o.dt time step [s] (resolve the blade passage: about 0.02 - 0.05 s for multi-MW rotors)
 * @param {number} o.nSteps number of steps
 * @param {number|number[]|function(number):number} o.wind rotor-effective hub-height speed [m/s]: constant, series
 *   at t = k dt, or function of time [s]
 * @param {function(number, number, number):number} [o.windField] local streamwise speed (t [s], y lateral [m],
 *   z above hub [m]) -> [m/s]; overrides `wind`, shear and tower shadow
 * @param {number|number[]|function(number):number} o.omega rotor speed [rad/s]
 * @param {number|number[]|function(number):number} [o.pitchDeg=0] collective pitch [deg]
 * @param {number} [o.rho=1.225] [kg/m3]
 * @param {number} [o.shearExp=0] @param {number} [o.hubHeight]
 * @param {object|boolean} [o.towerShadow=false] as in bemSolve
 * @param {number} [o.tiltDeg] @param {number} [o.preconeDeg]
 * @param {'oye'|'none'} [o.dynamicInflow='oye']
 * @param {'none'|'oye'|'bl'} [o.dynamicStall='none']
 * @param {object} [o.stall] constants forwarded to the dynamic-stall step (A, Tp, Tf, vortexLift, ...)
 * @param {number} [o.azimuth0=0] azimuth of blade 1 at t = 0 [rad]
 * @returns {{t:Float64Array, thrust:Float64Array, torque:Float64Array, power:Float64Array, azimuth:Float64Array,
 *   flapMoment:Float64Array[], edgeMoment:Float64Array[], wind:Float64Array}} time [s]; rotor thrust [N], torque
 *   [N m], aerodynamic power [W]; azimuth of blade 1 [rad]; per-blade aerodynamic root flapwise/edgewise moments
 *   [N m] (index [blade][step]); hub wind speed used [m/s]
 */
export function unsteadyBem(o) {
  const rotor = o.rotor;
  need(rotor && rotor.stations, 'rotor', 'rotor required');
  need(o.dt > 0, 'dt', 'must be > 0');
  need(Number.isInteger(o.nSteps) && o.nSteps >= 1, 'nSteps', 'integer >= 1 required');
  const rho = o.rho ?? 1.225;
  const B = rotor.blades;
  const st = rotor.stations;
  const n = st.length;
  const Rtip = rotor.tipRadius; const Rhub = rotor.hubRadius;
  const tilt = (o.tiltDeg ?? rotor.tiltDeg ?? 0) * D2R;
  const cone = (o.preconeDeg ?? rotor.preconeDeg ?? 0) * D2R;
  const sT = Math.sin(tilt); const cT = Math.cos(tilt); const sc = Math.sin(cone); const cc = Math.cos(cone);
  const shear = o.shearExp ?? 0;
  const hubHt = o.hubHeight ?? rotor.hubHeight;
  if (shear !== 0) need(hubHt > 0, 'hubHeight', 'required when shearExp != 0');
  let ts = null;
  if (o.towerShadow) {
    ts = o.towerShadow === true
      ? { topDiameter: rotor.tower?.topDiameter, baseDiameter: rotor.tower?.baseDiameter, height: rotor.tower?.height, overhang: rotor.overhang, hubAboveTowerTop: 0 }
      : { hubAboveTowerTop: 0, ...o.towerShadow };
    need(ts.topDiameter > 0 && ts.baseDiameter > 0 && ts.height > 0 && ts.overhang > 0, 'towerShadow', 'topDiameter, baseDiameter, height and overhang required');
  }
  const dynIn = (o.dynamicInflow ?? 'oye') === 'oye';
  const stallMode = o.dynamicStall ?? 'none';
  need(['none', 'oye', 'bl'].includes(stallMode), 'dynamicStall', "must be 'none', 'oye' or 'bl'");
  const seps = stallMode === 'none' ? null : st.map((s) => (s.polar.clAt(0.1) === 0 && s.polar.clAt(-0.1) === 0 ? null : createSeparationModel(s.polar)));
  const N = o.nSteps + 1;
  const res = {
    t: new Float64Array(N), thrust: new Float64Array(N), torque: new Float64Array(N), power: new Float64Array(N),
    azimuth: new Float64Array(N), wind: new Float64Array(N), flapMoment: [], edgeMoment: [],
  };
  for (let b = 0; b < B; b++) { res.flapMoment.push(new Float64Array(N)); res.edgeMoment.push(new Float64Array(N)); }
  const ax = []; const tn = []; const ds = [];
  for (let i = 0; i < B * n; i++) { ax.push({ xPrev: 0, y: 0, W: 0 }); tn.push({ xPrev: 0, y: 0, W: 0 }); ds.push(null); }
  const el = { Vx: 0, Vy: 0, r: 0, chord: 0, theta: 0, polar: null, blades: B, tipRadius: Rtip, hubRadius: Rhub, rho };
  const out = {};
  const stallArgs = { ...(o.stall || {}), sep: null, polar: null, alpha: 0, W: 0, chord: 0, dt: o.dt };
  const Np = new Float64Array(n); const Tp = new Float64Array(n);
  let psi1 = o.azimuth0 ?? 0;
  for (let k = 0; k < N; k++) {
    const t = k * o.dt;
    const U = o.windField ? o.windField(t, 0, 0) : seriesAt(o.wind, k, t, 0);
    const om = seriesAt(o.omega, k, t, 0);
    const pitch = seriesAt(o.pitchDeg, k, t, 0) * D2R;
    need(U >= 0 && om >= 0, 'wind/omega', 'must be >= 0');
    let T = 0; let Q = 0;
    for (let b = 0; b < B; b++) {
      const psi = psi1 + TWO_PI * b / B;
      const sa = Math.sin(psi); const ca = Math.cos(psi);
      for (let i = 0; i < n; i++) {
        const s = st[i];
        const r = s.r;
        const z = r * cc * ca * cT + r * sc * sT;
        let V;
        if (o.windField) V = o.windField(t, -r * cc * sa, z);
        else {
          V = U;
          if (shear !== 0) V *= Math.pow(Math.max(1e-3, 1 + z / hubHt), shear);
          if (ts && -z - ts.hubAboveTowerTop > 0) {
            const zb = -z - ts.hubAboveTowerTop;
            const rad = 0.5 * (ts.topDiameter + (ts.baseDiameter - ts.topDiameter) * Math.min(1, zb / ts.height));
            V *= towerPotentialFactor(ts.overhang * cT + r * (sc * cT - cc * ca * sT), r * cc * sa, rad);
          }
        }
        el.Vx = V * (sT * ca * sc + cT * cc);
        el.Vy = V * sT * sa + om * r * cc;
        el.r = r; el.chord = s.chord; el.theta = s.twistDeg * D2R + pitch; el.polar = s.polar;
        bemElement(el, out);
        const id = b * n + i;
        let wA = out.a * el.Vx;
        let wT = out.ap * el.Vy;
        if (dynIn && out.converged) {
          const A = ax[id]; const Tt = tn[id];
          if (k === 0) { A.xPrev = wA; A.y = wA; A.W = wA; Tt.xPrev = wT; Tt.y = wT; Tt.W = wT; } else {
            const aa = Math.min(0.5, Math.max(0, out.a));
            const tau1 = 1.1 / (1 - 1.3 * aa) * Rtip / Math.max(U, 0.1);
            const tau2 = (0.39 - 0.26 * (r / Rtip) ** 2) * tau1;
            const kk = 0.6;
            let H = wA + kk * tau1 * (wA - A.xPrev) / o.dt;
            A.y = H + (A.y - H) * Math.exp(-o.dt / tau1); A.W = A.y + (A.W - A.y) * Math.exp(-o.dt / tau2); A.xPrev = wA;
            H = wT + kk * tau1 * (wT - Tt.xPrev) / o.dt;
            Tt.y = H + (Tt.y - H) * Math.exp(-o.dt / tau1); Tt.W = Tt.y + (Tt.W - Tt.y) * Math.exp(-o.dt / tau2); Tt.xPrev = wT;
            wA = A.W; wT = Tt.W;
          }
        }
        const vx = el.Vx - wA; const vy = el.Vy + wT;
        const phi = Math.atan2(vx, vy);
        const W = Math.hypot(vx, vy);
        const alpha = phi - el.theta;
        let cl = s.polar.clAt(alpha); let cd = s.polar.cdAt(alpha);
        const sep = seps ? seps[i] : null;
        if (sep) {
          if (stallMode === 'oye') {
            if (ds[id] === null) ds[id] = { fs: sep.fAt(alpha), cl: 0 };
            stallArgs.sep = sep; stallArgs.alpha = alpha; stallArgs.W = W; stallArgs.chord = s.chord;
            cl = oyeDynamicStallStep(ds[id], stallArgs).cl;
          } else {
            if (ds[id] === null) ds[id] = beddoesLeishmanInit(sep, alpha, o.stall);
            stallArgs.sep = sep; stallArgs.polar = s.polar; stallArgs.alpha = alpha; stallArgs.W = W; stallArgs.chord = s.chord;
            const d = beddoesLeishmanStep(ds[id], stallArgs);
            cl = d.cl; cd = d.cd;
          }
        }
        const q = 0.5 * rho * W * W * s.chord;
        const sp = Math.sin(phi); const cp = Math.cos(phi);
        Np[i] = (cl * cp + cd * sp) * q;
        Tp[i] = (cl * sp - cd * cp) * q;
      }
      let T1 = 0; let Q1 = 0; let Mf = 0; let Me = 0; let rP = Rhub; let nP = 0; let tP = 0;
      for (let i = 0; i <= n; i++) {
        const r = i < n ? st[i].r : Rtip; const nn = i < n ? Np[i] : 0; const tt = i < n ? Tp[i] : 0; const dr = r - rP;
        T1 += 0.5 * (nn + nP) * dr; Q1 += 0.5 * (tt * r + tP * rP) * dr;
        Mf += 0.5 * (nn * (r - Rhub) + nP * (rP - Rhub)) * dr; Me += 0.5 * (tt * (r - Rhub) + tP * (rP - Rhub)) * dr;
        rP = r; nP = nn; tP = tt;
      }
      T += T1 * cc; Q += Q1 * cc;
      res.flapMoment[b][k] = Mf; res.edgeMoment[b][k] = Me;
    }
    res.t[k] = t; res.thrust[k] = T; res.torque[k] = Q; res.power[k] = Q * om; res.azimuth[k] = psi1; res.wind[k] = U;
    psi1 += om * o.dt;
  }
  return res;
}
