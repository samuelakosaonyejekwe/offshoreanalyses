/**
 * Variable-speed, variable-pitch turbine control: steady operating schedule, region-2 torque law, gain-scheduled PI
 * collective pitch, pitch actuator, supervisory sequences, yaw logic and a 1-DOF closed-loop simulation.
 *
 * The control structure follows the baseline controller of the NREL 5-MW definition (Jonkman et al. 2009,
 * NREL/TP-500-38060, section 7): generator-torque regions 1, 1.5, 2, 2.5, 3 scheduled on low-pass-filtered generator
 * speed and a gain-scheduled PI pitch loop on generator-speed error. Gains are computed from the caller's own
 * Cp/Ct surfaces; nothing here is tuned to a specific machine.
 *
 * Units: pitch angles in rad unless the name ends in `Deg`; "Lss" = low-speed (rotor) shaft, "Hss" = generator shaft.
 */

import { surfaceLookup } from './bem.js';

const D2R = Math.PI / 180;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/**
 * Region-2 optimal torque gain: generator torque T = k omega^2 keeps the rotor at the optimum tip-speed ratio,
 *   k = 0.5 rho pi R^5 Cp,max / lambda_opt^3   (low-speed shaft),   k_hss = k / N^3   (generator shaft, ideal gearbox).
 * Reference: Jonkman et al. (2009), section 7.2; Burton et al., Wind Energy Handbook, section 8.3.2. Mechanical
 * losses are not included: with a drivetrain efficiency below 1 applied to electrical power the gain is unchanged.
 * @param {object} o
 * @param {number} o.R rotor radius [m]
 * @param {number} o.cpMax maximum power coefficient [-]
 * @param {number} o.tsrOpt optimum tip-speed ratio [-]
 * @param {number} [o.rho=1.225] air density [kg/m3]
 * @param {number} [o.gearboxRatio=1] N [-]
 * @returns {{kLss:number, kHss:number, kHssRpm:number}} gains [N m/(rad/s)^2] on rotor and generator shaft, and the
 *   generator-shaft gain per rpm^2 [N m/rpm^2]
 */
export function region2Gain(o) {
  need(o && o.R > 0, 'R', 'must be > 0');
  need(o.cpMax > 0 && o.cpMax < 16 / 27 + 1e-9, 'cpMax', 'must be in (0, 16/27]');
  need(o.tsrOpt > 0, 'tsrOpt', 'must be > 0');
  const N = o.gearboxRatio ?? 1;
  const kLss = 0.5 * (o.rho ?? 1.225) * Math.PI * o.R ** 5 * o.cpMax / o.tsrOpt ** 3;
  const kHss = kLss / N ** 3;
  return { kLss, kHss, kHssRpm: kHss * (Math.PI / 30) ** 2 };
}

/**
 * PI collective-pitch gains by pole placement on the 1-DOF rotor model linearised at rated speed
 * (Jonkman et al. 2009, Eqs. 7-13 and 7-14; Hansen et al. 2005, Riso-R-1500):
 *   Kp = 2 I Omega0 zeta omega_n / (N (-dP/dtheta)),   Ki = I Omega0 omega_n^2 / (N (-dP/dtheta)),
 * acting on generator-speed error [rad/s] and giving pitch [rad]. Worked check: I = 43,784,700 kg m2,
 * Omega0 = 12.1 rpm, N = 97, zeta = 0.7, omega_n = 0.6 rad/s, dP/dtheta = -25.52e6 W/rad give Kp = 0.01883 s and
 * Ki = 0.008069 (the published NREL 5-MW baseline values).
 * @param {object} o
 * @param {number} o.inertia drivetrain inertia referred to the low-speed shaft [kg m2]
 * @param {number} o.omegaRated rated rotor speed [rad/s]
 * @param {number} o.dPdTheta sensitivity of aerodynamic power to collective pitch at the operating point [W/rad] (< 0)
 * @param {number} [o.gearboxRatio=1] N [-]
 * @param {number} [o.zeta=0.7] target damping ratio [-] (0.6 - 0.7 recommended)
 * @param {number} [o.omegaN=0.6] target natural frequency [rad/s] (0.6 for bottom-fixed; lower for floaters)
 * @returns {{kp:number, ki:number}} proportional gain [s] and integral gain [-]
 */
export function piPitchGains(o) {
  need(o && o.inertia > 0, 'inertia', 'must be > 0');
  need(o.omegaRated > 0, 'omegaRated', 'must be > 0');
  need(o.dPdTheta < 0, 'dPdTheta', 'must be < 0 (power falls as the blade pitches to feather)');
  const N = o.gearboxRatio ?? 1;
  const z = o.zeta ?? 0.7;
  const wn = o.omegaN ?? 0.6;
  need(z > 0 && wn > 0, 'zeta/omegaN', 'must be > 0');
  return { kp: 2 * o.inertia * o.omegaRated * z * wn / (N * -o.dPdTheta), ki: o.inertia * o.omegaRated * wn * wn / (N * -o.dPdTheta) };
}

/**
 * Gain-correction (scheduling) factor GK(theta) = 1 / (1 + theta/thetaK), where thetaK is the pitch angle at which
 * the pitch sensitivity dP/dtheta has doubled from its value at zero pitch (Jonkman et al. 2009, Eq. 7-18).
 * @param {number} theta collective pitch [rad]
 * @param {number} thetaK [rad] (> 0)
 * @returns {number} GK [-]
 */
export function gainCorrection(theta, thetaK) {
  need(thetaK > 0, 'thetaK', 'must be > 0');
  return 1 / (1 + Math.max(0, theta) / thetaK);
}

function aeroPower(surface, rho, U, omega, pitchDeg) {
  const R = surface.R;
  const c = surfaceLookup(surface, omega * R / Math.max(U, 1e-6), pitchDeg);
  const q = 0.5 * rho * Math.PI * R * R * U * U;
  return { power: c.cp * q * U, thrust: c.ct * q, torque: c.cq * q * R, cp: c.cp, ct: c.ct };
}

function pitchForPower(surface, rho, U, omega, target, fineDeg) {
  let lo = fineDeg;
  const pMax = surface.pitchDeg[surface.pitchDeg.length - 1];
  if (aeroPower(surface, rho, U, omega, lo).power <= target) return lo;
  let hi = lo;
  // march to the first pitch that drops below the target (power is not monotonic far beyond feather)
  do { lo = hi; hi = Math.min(pMax, hi + 1); } while (aeroPower(surface, rho, U, omega, hi).power > target && hi < pMax);
  if (aeroPower(surface, rho, U, omega, hi).power > target) return pMax;
  for (let i = 0; i < 50; i++) {
    const mid = 0.5 * (lo + hi);
    if (aeroPower(surface, rho, U, omega, mid).power > target) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

/**
 * Optimum of a Cp surface: grid maximum refined by a parabola through the neighbouring nodes in each direction.
 * @param {object} surface cpSurface result (bem.js)
 * @param {object} [o]
 * @param {number} [o.pitchDeg] restrict the search to this pitch [deg]
 * @returns {{cpMax:number, tsrOpt:number, pitchOptDeg:number}} [-], [-], [deg]
 */
export function surfaceOptimum(surface, o = {}) {
  const T = surface.tsr; const P = surface.pitchDeg;
  if (o.pitchDeg !== undefined) {
    let best = -Infinity; let lb = T[0];
    const n = 400;
    for (let i = 0; i <= n; i++) {
      const l = T[0] + (T[T.length - 1] - T[0]) * i / n;
      const c = surfaceLookup(surface, l, o.pitchDeg).cp;
      if (c > best) { best = c; lb = l; }
    }
    return { cpMax: best, tsrOpt: lb, pitchOptDeg: o.pitchDeg };
  }
  let bi = 0; let bj = 0; let best = -Infinity;
  for (let j = 0; j < P.length; j++) for (let i = 0; i < T.length; i++) if (surface.cp[j][i] > best) { best = surface.cp[j][i]; bi = i; bj = j; }
  const vertex = (x0, x1, x2, y0, y1, y2) => {
    const d = (x0 - x1) * (x0 - x2) * (x1 - x2);
    const A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d;
    const Bq = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d;
    if (!(A < 0)) return [x1, y1];
    const xv = -Bq / (2 * A);
    const C = y1 - A * x1 * x1 - Bq * x1;
    return [xv, A * xv * xv + Bq * xv + C];
  };
  let tsrOpt = T[bi]; let pitchOpt = P[bj]; let cp = best;
  if (bi > 0 && bi < T.length - 1) {
    const v = vertex(T[bi - 1], T[bi], T[bi + 1], surface.cp[bj][bi - 1], best, surface.cp[bj][bi + 1]);
    tsrOpt = v[0]; cp = Math.max(cp, v[1]);
  }
  if (bj > 0 && bj < P.length - 1) {
    const v = vertex(P[bj - 1], P[bj], P[bj + 1], surface.cp[bj - 1][bi], best, surface.cp[bj + 1][bi]);
    pitchOpt = v[0]; cp = Math.max(cp, v[1]);
  }
  return { cpMax: cp, tsrOpt, pitchOptDeg: pitchOpt };
}

/**
 * Steady operating schedule versus wind speed from Cp/Ct surfaces.
 * Regions: 1.5 (rotor speed held at its minimum), 2 (optimum tip-speed ratio, fine pitch), 2.5 (rotor speed held
 * at rated, power still below rated), 3 (rated power held by pitching to feather at rated speed). Optional peak
 * shaving: wherever the thrust exceeds `maxThrust` the pitch is increased until the thrust equals the limit (power
 * is given up around rated wind speed).
 * @param {object} o
 * @param {object} o.surface cpSurface result (bem.js), computed for the same rotor
 * @param {number} o.ratedPower electrical rated power [W]
 * @param {number} o.ratedRotorSpeed [rad/s]
 * @param {number} [o.minRotorSpeed=0] [rad/s]
 * @param {number} [o.efficiency=1] drivetrain efficiency [-]
 * @param {number} [o.gearboxRatio=1] [-]
 * @param {number} [o.rho=1.225] [kg/m3]
 * @param {number} [o.tsrOpt] @param {number} [o.pitchFineDeg] default: optimum of the surface
 * @param {number[]} [o.speeds] wind speeds [m/s]; default cutIn..cutOut step 0.5
 * @param {number} [o.cutIn=3] @param {number} [o.cutOut=25] [m/s]
 * @param {number} [o.maxThrust] thrust limit for peak shaving [N]
 * @returns {{U:number[],omega:number[],pitchDeg:number[],power:number[],aeroPower:number[],thrust:number[],
 *   aeroTorque:number[],genTorque:number[],cp:number[],ct:number[],tsr:number[],region:string[],ratedWind:number,
 *   cpMax:number,tsrOpt:number,pitchFineDeg:number}} rotor speed [rad/s], pitch [deg], electrical and aerodynamic
 *   power [W], thrust [N], aerodynamic torque on the rotor shaft and generator torque on the generator shaft [N m],
 *   coefficients, region, rated wind speed [m/s]
 */
export function operatingSchedule(o) {
  const sf = o.surface;
  need(sf && sf.cp && sf.R > 0, 'surface', 'cpSurface result required');
  need(o.ratedPower > 0, 'ratedPower', 'must be > 0');
  need(o.ratedRotorSpeed > 0, 'ratedRotorSpeed', 'must be > 0');
  const rho = o.rho ?? 1.225;
  const eta = o.efficiency ?? 1;
  const N = o.gearboxRatio ?? 1;
  const wMin = o.minRotorSpeed ?? 0;
  const wR = o.ratedRotorSpeed;
  const cutIn = o.cutIn ?? 3; const cutOut = o.cutOut ?? 25;
  const opt = surfaceOptimum(sf, o.pitchFineDeg !== undefined ? { pitchDeg: o.pitchFineDeg } : {});
  const tsrOpt = o.tsrOpt ?? opt.tsrOpt;
  const fine = o.pitchFineDeg ?? opt.pitchOptDeg;
  const R = sf.R;
  const Paero = o.ratedPower / eta;
  const speeds = o.speeds ?? (() => { const a = []; for (let u = cutIn; u <= cutOut + 1e-9; u += 0.5) a.push(u); return a; })();
  const wAt = (U) => Math.min(wR, Math.max(wMin, tsrOpt * U / R));
  const out = { U: [], omega: [], pitchDeg: [], power: [], aeroPower: [], thrust: [], aeroTorque: [], genTorque: [], cp: [], ct: [], tsr: [], region: [],
    ratedWind: NaN, cpMax: surfaceLookup(sf, tsrOpt, fine).cp, tsrOpt, pitchFineDeg: fine };
  {
    const g = (U) => aeroPower(sf, rho, U, wAt(U), fine).power - Paero;
    let lo = Math.max(cutIn, 0.5);
    if (g(lo) < 0) {
      let x = lo;
      while (x < cutOut && g(Math.min(cutOut, x + 0.5)) < 0) x += 0.5;
      if (x < cutOut) {
        lo = x; let hi = Math.min(cutOut, x + 0.5);
        for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (g(mid) < 0) lo = mid; else hi = mid; }
        out.ratedWind = 0.5 * (lo + hi);
      }
    }
  }
  for (const U of speeds) {
    let w = 0; let p = fine; let s = { power: 0, thrust: 0, torque: 0, cp: 0, ct: 0 }; let region = 'off';
    if (U >= cutIn - 1e-9 && U <= cutOut + 1e-9 && U > 0) {
      w = wAt(U);
      const wOpt = tsrOpt * U / R;
      region = wOpt < wMin ? '1.5' : (wOpt > wR ? '2.5' : '2');
      s = aeroPower(sf, rho, U, w, fine);
      if (s.power > Paero) { region = '3'; p = pitchForPower(sf, rho, U, w, Paero, fine); s = aeroPower(sf, rho, U, w, p); }
      if (o.maxThrust !== undefined && s.thrust > o.maxThrust) {
        let lo = p; let hi = sf.pitchDeg[sf.pitchDeg.length - 1];
        for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (aeroPower(sf, rho, U, w, mid).thrust > o.maxThrust) lo = mid; else hi = mid; }
        p = 0.5 * (lo + hi); s = aeroPower(sf, rho, U, w, p);
        region += ' (peak-shaved)';
      }
    }
    out.U.push(U); out.omega.push(w); out.pitchDeg.push(p); out.region.push(region);
    out.aeroPower.push(s.power); out.power.push(Math.max(0, eta * s.power)); out.thrust.push(s.thrust);
    out.aeroTorque.push(s.torque); out.genTorque.push(w > 0 ? Math.max(0, s.power) / (w * N) : 0);
    out.cp.push(s.cp); out.ct.push(s.ct); out.tsr.push(U > 0 ? w * R / U : 0);
  }
  return out;
}

/**
 * Sensitivity of aerodynamic power to collective pitch along the region-3 operating line (rated speed, rated
 * power), from the Cp surface by central differences, with the linear fit dP/dtheta = s0 (1 + theta/thetaK) (theta
 * measured from fine pitch) used for gain scheduling (Jonkman et al. 2009, Fig. 7-3: the published NREL 5-MW values are s0 = -25.52e6 W/rad and
 * thetaK = 6.302336 deg).
 * @param {object} o
 * @param {object} o.surface cpSurface result
 * @param {number} o.ratedPowerAero aerodynamic (shaft) rated power [W]
 * @param {number} o.ratedRotorSpeed [rad/s]
 * @param {number} [o.rho=1.225]
 * @param {number} [o.pitchFineDeg=0] fine pitch [deg]
 * @param {number[]} [o.speeds] wind speeds above rated [m/s]; default rated..25 in 1 m/s steps
 * @param {number} [o.stepDeg=1] half-width of the central difference [deg]
 * @returns {{U:number[], pitchDeg:number[], dPdTheta:number[], dPdTheta0:number, thetaK:number, thetaKDeg:number, ratedWind:number}}
 *   operating points, sensitivities [W/rad], fitted value at zero pitch [W/rad], thetaK [rad] and [deg]
 */
export function pitchSensitivity(o) {
  const sf = o.surface;
  const rho = o.rho ?? 1.225;
  const fine = o.pitchFineDeg ?? 0;
  const h = o.stepDeg ?? 1;
  const w = o.ratedRotorSpeed;
  need(o.ratedPowerAero > 0 && w > 0, 'ratedPowerAero/ratedRotorSpeed', 'must be > 0');
  let lo = 1; let hi = 40;
  need(aeroPower(sf, rho, hi, w, fine).power > o.ratedPowerAero, 'ratedPowerAero', 'not reached below 40 m/s on this surface');
  for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (aeroPower(sf, rho, mid, w, fine).power < o.ratedPowerAero) lo = mid; else hi = mid; }
  const Ur = 0.5 * (lo + hi);
  const speeds = o.speeds ?? (() => { const a = [Ur]; for (let u = Math.ceil(Ur + 0.5); u <= 25; u += 1) a.push(u); return a; })();
  const res = { U: [], pitchDeg: [], dPdTheta: [], dPdTheta0: 0, thetaK: 0, thetaKDeg: 0, ratedWind: Ur };
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0;
  for (const U of speeds) {
    const p = pitchForPower(sf, rho, U, w, o.ratedPowerAero, fine);
    const d = (aeroPower(sf, rho, U, w, p + h).power - aeroPower(sf, rho, U, w, p - h).power) / (2 * h * D2R);
    res.U.push(U); res.pitchDeg.push(p); res.dPdTheta.push(d);
    const x = (p - fine) * D2R;       // pitch measured from fine pitch, as used by the gain schedule
    sx += x; sy += d; sxx += x * x; sxy += x * d;
  }
  const n = speeds.length;
  need(n >= 2, 'speeds', 'at least two operating points required');
  const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const s0 = (sy - slope * sx) / n;
  res.dPdTheta0 = s0;
  res.thetaK = s0 / slope;
  res.thetaKDeg = res.thetaK / D2R;
  return res;
}

/**
 * Derived constants of the five-region generator-torque law (Jonkman et al. 2009, section 7.2; all on the
 * generator shaft): region 1.5 is a linear ramp from zero torque at the cut-in generator speed to the region-2
 * curve at (1 + region15Span) times that speed; region 2.5 is a linear ramp of an induction-machine-like slope with
 * synchronous speed ratedGenSpeed*region3Fraction/(1 + slip) that meets rated torque at region3Fraction of rated
 * speed; region 3 holds constant power (or constant torque).
 * @param {object} c
 * @param {number} c.k region-2 gain on the generator shaft [N m/(rad/s)^2]
 * @param {number} c.ratedGenSpeed rated generator speed [rad/s]
 * @param {number} c.ratedPowerMech rated mechanical generator power, electrical/efficiency [W]
 * @param {number} [c.minGenSpeed=0] cut-in generator speed [rad/s]
 * @param {number} [c.region15Span=0.3] [-]
 * @param {number} [c.slipPercent=10] region-2.5 slip [%]
 * @param {number} [c.region3Fraction=0.99] fraction of rated speed at which region 3 starts [-]
 * @returns {{k:number, rtGnSp:number, rtPwr:number, rtTq:number, ctInSp:number, rgn2Sp:number, sySp:number,
 *   slope15:number, slope25:number, trGnSp:number}} speeds [rad/s], power [W], torque [N m], slopes [N m s/rad]
 */
export function torqueLawConstants(c) {
  need(c && c.k >= 0, 'k', 'must be >= 0');
  need(c.ratedGenSpeed > 0 && c.ratedPowerMech > 0, 'ratedGenSpeed/ratedPowerMech', 'must be > 0');
  const rtGnSp = (c.region3Fraction ?? 0.99) * c.ratedGenSpeed;
  const rtPwr = c.ratedPowerMech;
  const rtTq = rtPwr / rtGnSp;
  const ctInSp = c.minGenSpeed ?? 0;
  const rgn2Sp = ctInSp * (1 + (c.region15Span ?? 0.3));
  const sySp = rtGnSp / (1 + 0.01 * (c.slipPercent ?? 10));
  const slope15 = rgn2Sp > ctInSp ? c.k * rgn2Sp * rgn2Sp / (rgn2Sp - ctInSp) : 0;
  const slope25 = rtTq / (rtGnSp - sySp);
  const disc = slope25 * (slope25 - 4 * c.k * sySp);
  const trGnSp = c.k === 0 ? sySp : (disc >= 0 ? (slope25 - Math.sqrt(disc)) / (2 * c.k) : rtGnSp);
  return { k: c.k, rtGnSp, rtPwr, rtTq, ctInSp, rgn2Sp, sySp, slope15, slope25, trGnSp };
}

/**
 * Generator torque demand from the five-region law.
 * @param {object} tc constants from torqueLawConstants
 * @param {number} genSpeed filtered generator speed [rad/s]
 * @param {object} [o]
 * @param {number} [o.pitch=0] collective pitch [rad]
 * @param {number} [o.region3MinPitch=0.0174533] pitch above which region 3 is forced [rad] (1 deg)
 * @param {boolean} [o.constantPower=true] region 3 holds power (true) or rated torque (false)
 * @returns {{torque:number, region:string}} generator-shaft torque [N m], region label
 */
export function generatorTorqueLaw(tc, genSpeed, o = {}) {
  const w = genSpeed;
  if (w >= tc.rtGnSp || (o.pitch ?? 0) >= (o.region3MinPitch ?? D2R)) {
    return { torque: (o.constantPower ?? true) ? tc.rtPwr / Math.max(w, 0.2 * tc.rtGnSp) : tc.rtTq, region: '3' };
  }
  if (w <= tc.ctInSp) return { torque: 0, region: '1' };
  if (w < tc.rgn2Sp) return { torque: tc.slope15 * (w - tc.ctInSp), region: '1.5' };
  if (w < tc.trGnSp) return { torque: tc.k * w * w, region: '2' };
  return { torque: Math.max(0, tc.slope25 * (w - tc.sySp)), region: '2.5' };
}

/**
 * One step of a second-order pitch actuator with position and rate limits:
 *   theta'' + 2 zeta w theta' + w^2 theta = w^2 theta_cmd, |theta'| <= rateLimit, min <= theta <= max,
 * integrated with the exact discrete solution of the unsaturated oscillator over sub-steps (semi-implicit Euler
 * with dt/8 sub-steps), then clipped.
 * @param {{theta:number, rate:number}} state pitch [rad] and pitch rate [rad/s], mutated
 * @param {object} o
 * @param {number} o.cmd commanded pitch [rad]
 * @param {number} o.dt step [s]
 * @param {number} [o.omegaAct=6.283] actuator natural frequency [rad/s] (typically 1 - 2 Hz)
 * @param {number} [o.zetaAct=0.7] actuator damping ratio [-]
 * @param {number} [o.rateLimit=0.1396] [rad/s] (8 deg/s)
 * @param {number} [o.min=0] [rad] @param {number} [o.max=1.5708] [rad]
 * @returns {object} state
 */
export function pitchActuatorStep(state, o) {
  const w = o.omegaAct ?? 2 * Math.PI;
  const z = o.zetaAct ?? 0.7;
  const rl = o.rateLimit ?? 8 * D2R;
  const lo = o.min ?? 0; const hi = o.max ?? Math.PI / 2;
  const nSub = Math.max(1, Math.ceil(o.dt * w * 8));
  const h = o.dt / nSub;
  for (let i = 0; i < nSub; i++) {
    const acc = w * w * (o.cmd - state.theta) - 2 * z * w * state.rate;
    state.rate += h * acc;
    if (state.rate > rl) state.rate = rl; else if (state.rate < -rl) state.rate = -rl;
    state.theta += h * state.rate;
    if (state.theta > hi) { state.theta = hi; if (state.rate > 0) state.rate = 0; } else if (state.theta < lo) { state.theta = lo; if (state.rate < 0) state.rate = 0; }
  }
  return state;
}

/**
 * Controller configuration for a rotor from its Cp surface: region-2 gain, torque-law constants, pitch sensitivity,
 * PI gains and gain-schedule angle, with NREL-baseline-style limits as editable defaults.
 * @param {object} o
 * @param {object} o.rotor rotor with ratedPower [W], ratedRotorSpeed, minRotorSpeed [rad/s], gearboxRatio,
 *   drivetrainEfficiency, drivetrainInertia [kg m2]
 * @param {object} o.surface cpSurface result for the rotor
 * @param {number} [o.rho=1.225]
 * @param {number} [o.zeta=0.7] @param {number} [o.omegaN=0.6] pitch-loop design targets [-], [rad/s]
 * @param {number} [o.pitchFineDeg] fine pitch [deg]; default optimum of the surface, not below 0
 * @param {number} [o.pitchRateDeg=8] pitch rate limit [deg/s] (NREL 5-MW: 8)
 * @param {number} [o.maxTorqueFactor=1.1] torque limit / rated torque [-] (NREL 5-MW: 1.1)
 * @param {number} [o.torqueRate] generator torque rate limit [N m/s]; default rated torque per 2.87 s
 *   (the NREL 5-MW ratio 43,094 N m / 15,000 N m/s)
 * @param {number} [o.filterHz=0.25] corner frequency of the generator-speed low-pass filter [Hz] (NREL 5-MW: 0.25)
 * @returns {object} controller configuration for createController, with `design` (cpMax, tsrOpt, dPdTheta0, ...)
 */
export function controllerFromRotor(o) {
  const r = o.rotor; const sf = o.surface;
  need(r && r.ratedPower > 0 && r.ratedRotorSpeed > 0, 'rotor', 'ratedPower and ratedRotorSpeed required');
  need(r.drivetrainInertia > 0, 'rotor.drivetrainInertia', 'must be > 0');
  const rho = o.rho ?? 1.225;
  const N = r.gearboxRatio ?? 1;
  const eta = r.drivetrainEfficiency ?? 1;
  let opt = surfaceOptimum(sf);
  let fine = o.pitchFineDeg ?? Math.max(0, opt.pitchOptDeg);
  if (fine !== opt.pitchOptDeg) opt = surfaceOptimum(sf, { pitchDeg: fine });
  const g = region2Gain({ R: sf.R, cpMax: opt.cpMax, tsrOpt: opt.tsrOpt, rho, gearboxRatio: N });
  const Pmech = r.ratedPower / eta;
  const sens = pitchSensitivity({ surface: sf, ratedPowerAero: Pmech, ratedRotorSpeed: r.ratedRotorSpeed, rho, pitchFineDeg: fine });
  const gains = piPitchGains({ inertia: r.drivetrainInertia, omegaRated: r.ratedRotorSpeed, dPdTheta: sens.dPdTheta0, gearboxRatio: N, zeta: o.zeta, omegaN: o.omegaN });
  const ratedGenSpeed = r.ratedRotorSpeed * N;
  const ratedTorque = Pmech / ratedGenSpeed;
  return {
    gearboxRatio: N, efficiency: eta,
    k: g.kHss, ratedGenSpeed, ratedPowerMech: Pmech, ratedTorque, minGenSpeed: (r.minRotorSpeed ?? 0) * N,
    region15Span: 0.3, slipPercent: 10, region3Fraction: 0.99, constantPower: true,
    maxTorque: (o.maxTorqueFactor ?? 1.1) * ratedTorque, torqueRate: o.torqueRate ?? ratedTorque / 2.87,
    kp: gains.kp, ki: gains.ki, thetaK: sens.thetaK,
    pitchMin: fine * D2R, pitchMax: 90 * D2R, pitchRate: (o.pitchRateDeg ?? 8) * D2R, filterHz: o.filterHz ?? 0.25,
    towerDamperGain: 0, cutIn: r.cutIn, cutOut: r.cutOut,
    design: { cpMax: opt.cpMax, tsrOpt: opt.tsrOpt, pitchFineDeg: fine, dPdTheta0: sens.dPdTheta0, thetaKDeg: sens.thetaKDeg,
      ratedWind: sens.ratedWind, zeta: o.zeta ?? 0.7, omegaN: o.omegaN ?? 0.6, kHssRpm: g.kHssRpm },
  };
}

/**
 * Discrete baseline controller with a simple supervisory state machine.
 *
 * Modes: 'parked' (pitch at feather, no torque) -> 'startup' (pitch ramps to fine at `startupPitchRate`; generator
 * connects when the speed passes the cut-in generator speed) -> 'production' (torque law + gain-scheduled PI pitch
 * with integrator clamping as anti-windup; optional tower fore-aft damper adding towerDamperGain * tower-top
 * velocity to the pitch demand) -> 'shutdown' (normal stop: pitch to feather at the normal rate, torque law kept
 * until the speed drops below cut-in generator speed) or 'estop' (emergency stop: pitch to feather at `estopPitchRate`,
 * generator torque removed at once when `gridLoss`, mechanical brake `brakeTorque` on the generator shaft below
 * `brakeGenSpeed`). 'estop' is latched; 'shutdown' returns to 'parked' at standstill. Automatic transitions:
 * wind above cutOut or below cutIn (when `windSpeed` is supplied, compared against `cfg.cutOut`/`cfg.cutIn`), and
 * over-speed above `overspeedFactor` x rated (when given) -> estop. This is a documented minimal logic, not a
 * certified safety system.
 * @param {object} cfg configuration (controllerFromRotor result or equivalent): k, ratedGenSpeed, ratedPowerMech,
 *   minGenSpeed, region15Span, slipPercent, region3Fraction, constantPower, maxTorque [N m], torqueRate [N m/s],
 *   kp [s], ki [-], thetaK [rad], pitchMin, pitchMax [rad], pitchRate [rad/s], filterHz [Hz],
 *   towerDamperGain [rad s/m], startupPitchRate [rad/s] (default 2 deg/s), estopPitchRate [rad/s] (default pitchRate),
 *   brakeTorque [N m] (default 0), brakeGenSpeed [rad/s] (default 0), overspeedFactor [-] (default none),
 *   cutIn, cutOut [m/s]
 * @param {object} [init]
 * @param {'parked'|'startup'|'production'} [init.mode='production']
 * @param {number} [init.pitch] initial pitch [rad] (seeds the integrator); default pitchMin
 * @param {number} [init.genSpeed] initial generator speed [rad/s]; default rated
 * @param {number} [init.genTorque] initial torque [N m]; default from the torque law
 * @returns {{step:function(object):object, state:object, cfg:object}} controller; `step({ dt [s], genSpeed [rad/s],
 *   pitch [rad] (measured), towerTopVel [m/s]?, windSpeed [m/s]?, command?: 'start'|'stop'|'estop', gridLoss?: boolean })`
 *   returns { genTorque [N m], pitchCmd [rad], brakeTorque [N m], mode, region, genSpeedFiltered [rad/s] }
 */
export function createController(cfg, init = {}) {
  need(cfg && cfg.kp >= 0 && cfg.ki >= 0, 'cfg.kp/ki', 'must be >= 0');
  need(cfg.pitchMax > cfg.pitchMin, 'cfg.pitchMax', 'must exceed pitchMin');
  const tc = torqueLawConstants(cfg);
  const pRate = cfg.pitchRate ?? 8 * D2R;
  const st = {
    mode: init.mode ?? 'production',
    genSpeedF: init.genSpeed ?? cfg.ratedGenSpeed,
    pitchCmd: init.pitch ?? (init.mode === 'parked' ? cfg.pitchMax : cfg.pitchMin),
    integ: init.pitch ?? cfg.pitchMin,
    torque: 0, gridLost: false, time: 0, region: '',
  };
  st.torque = init.genTorque ?? (st.mode === 'production' ? Math.min(cfg.maxTorque ?? Infinity, generatorTorqueLaw(tc, st.genSpeedF, { pitch: st.pitchCmd, constantPower: cfg.constantPower }).torque) : 0);
  const res = { genTorque: 0, pitchCmd: 0, brakeTorque: 0, mode: st.mode, region: '', genSpeedFiltered: 0 };
  function slew(target, rate, dt) {
    const d = target - st.pitchCmd;
    const m = rate * dt;
    st.pitchCmd += d > m ? m : (d < -m ? -m : d);
  }
  function torqueTo(target, dt) {
    const tq = Math.min(cfg.maxTorque ?? Infinity, Math.max(0, target));
    const m = (cfg.torqueRate ?? Infinity) * dt;
    const d = tq - st.torque;
    st.torque += d > m ? m : (d < -m ? -m : d);
  }
  function step(inp) {
    const dt = inp.dt;
    need(dt > 0, 'dt', 'must be > 0');
    st.time += dt;
    const alpha = Math.exp(-2 * Math.PI * dt * (cfg.filterHz ?? 0.25));
    st.genSpeedF = (1 - alpha) * inp.genSpeed + alpha * st.genSpeedF;
    const w = st.genSpeedF;
    if (inp.gridLoss) st.gridLost = true;
    // supervisory transitions
    if (inp.command === 'estop' || (st.gridLost && st.mode !== 'parked')) st.mode = 'estop';
    else if (cfg.overspeedFactor && inp.genSpeed > cfg.overspeedFactor * cfg.ratedGenSpeed && st.mode !== 'parked') st.mode = 'estop';
    else if (inp.command === 'stop' && (st.mode === 'production' || st.mode === 'startup')) st.mode = 'shutdown';
    else if (inp.command === 'start' && st.mode === 'parked') st.mode = 'startup';
    if (inp.windSpeed !== undefined && st.mode === 'production' && cfg.cutOut !== undefined && inp.windSpeed > cfg.cutOut) st.mode = 'shutdown';
    let brake = 0;
    let region = '';
    if (st.mode === 'production') {
      const law = generatorTorqueLaw(tc, w, { pitch: st.pitchCmd, constantPower: cfg.constantPower, region3MinPitch: cfg.pitchMin + D2R });
      region = law.region;
      torqueTo(law.torque, dt);
      const gk = gainCorrection(inp.pitch - cfg.pitchMin, cfg.thetaK);
      const err = w - cfg.ratedGenSpeed;
      st.integ += gk * cfg.ki * err * dt;
      if (st.integ > cfg.pitchMax) st.integ = cfg.pitchMax; else if (st.integ < cfg.pitchMin) st.integ = cfg.pitchMin;
      let target = gk * cfg.kp * err + st.integ;
      if (cfg.towerDamperGain && inp.towerTopVel !== undefined) target += cfg.towerDamperGain * inp.towerTopVel;
      if (target > cfg.pitchMax) target = cfg.pitchMax; else if (target < cfg.pitchMin) target = cfg.pitchMin;
      slew(target, pRate, dt);
    } else if (st.mode === 'startup') {
      slew(cfg.pitchMin, cfg.startupPitchRate ?? 2 * D2R, dt);
      st.integ = st.pitchCmd;
      torqueTo(0, dt);
      region = 'startup';
      if (w > tc.ctInSp && st.pitchCmd <= cfg.pitchMin + 1e-9) st.mode = 'production';
      else if (w > Math.max(tc.rgn2Sp, tc.ctInSp)) { st.mode = 'production'; }
    } else if (st.mode === 'shutdown') {
      slew(cfg.pitchMax, pRate, dt);
      const law = generatorTorqueLaw(tc, w, { pitch: 0, constantPower: false });
      torqueTo(w > tc.ctInSp ? Math.min(law.torque, tc.rtTq) : 0, dt);
      region = 'shutdown';
      if (w < 0.02 * cfg.ratedGenSpeed && st.pitchCmd >= cfg.pitchMax - 1e-9) st.mode = 'parked';
    } else if (st.mode === 'estop') {
      slew(cfg.pitchMax, cfg.estopPitchRate ?? pRate, dt);
      if (st.gridLost) st.torque = 0; else torqueTo(0, dt);
      if (inp.genSpeed < (cfg.brakeGenSpeed ?? 0)) brake = cfg.brakeTorque ?? 0;
      region = 'estop';
    } else {
      slew(cfg.pitchMax, pRate, dt);
      torqueTo(0, dt);
      st.integ = cfg.pitchMin;
      region = 'parked';
    }
    st.region = region;
    res.genTorque = st.torque; res.pitchCmd = st.pitchCmd; res.brakeTorque = brake; res.mode = st.mode; res.region = region; res.genSpeedFiltered = w;
    return res;
  }
  return { step, state: st, cfg, torqueConstants: tc };
}

/**
 * One step of a dead-band yaw controller with delay: yawing starts when the low-pass-filtered yaw error has exceeded
 * the dead band continuously for `delay` seconds and then proceeds at constant rate until the filtered error has
 * changed sign or fallen below `stopDeg`. All thresholds are arguments (typical: dead band 8 - 15 deg, delay
 * 30 - 120 s, rate 0.3 - 0.5 deg/s).
 * @param {{errorF:number, timer:number, yawing:number}} state filtered error [deg], timer [s], direction (-1, 0, 1); mutated
 * @param {object} o
 * @param {number} o.errorDeg instantaneous yaw misalignment, wind direction minus nacelle direction [deg]
 * @param {number} o.dt [s]
 * @param {number} [o.deadbandDeg=8] @param {number} [o.delay=30] [s] @param {number} [o.rateDegPerS=0.3]
 * @param {number} [o.stopDeg=1] @param {number} [o.filterTau=10] error low-pass time constant [s]
 * @returns {{yawRateDegPerS:number, yawing:boolean}} nacelle yaw rate command [deg/s]
 */
export function yawControllerStep(state, o) {
  const tau = o.filterTau ?? 10;
  const db = o.deadbandDeg ?? 8;
  const a = Math.exp(-o.dt / tau);
  state.errorF = a * (state.errorF ?? 0) + (1 - a) * o.errorDeg;
  state.timer = state.timer ?? 0;
  state.yawing = state.yawing ?? 0;
  if (state.yawing === 0) {
    if (Math.abs(state.errorF) > db) state.timer += o.dt; else state.timer = 0;
    if (state.timer >= (o.delay ?? 30)) { state.yawing = Math.sign(state.errorF); state.timer = 0; }
  } else if (state.errorF * state.yawing <= (o.stopDeg ?? 1)) {
    state.yawing = 0;
  }
  return { yawRateDegPerS: state.yawing * (o.rateDegPerS ?? 0.3), yawing: state.yawing !== 0 };
}

/**
 * Floating-platform negative-damping check for the pitch controller: above rated the constant-power pitch loop
 * makes thrust fall with relative wind speed; to avoid driving platform pitch unstable the closed-loop rotor-speed
 * natural frequency is placed below the platform-pitch natural frequency (Larsen, T.J. and Hanson, T.D. (2007),
 * J. Phys. Conf. Ser. 75 012073; Jonkman (2008) reduced omega_n from 0.6 to 0.2 rad/s for the OC3-Hywind spar).
 * @param {object} o
 * @param {number} o.omegaN controller natural frequency [rad/s]
 * @param {number} o.platformPitchFreqHz platform pitch natural frequency [Hz]
 * @param {number} [o.margin=1] required ratio platform frequency / controller frequency [-]
 * @returns {{ok:boolean, ratio:number, omegaNMax:number}} pass flag, frequency ratio [-], largest admissible omega_n [rad/s]
 */
export function floatingDetuneCheck(o) {
  need(o && o.omegaN > 0 && o.platformPitchFreqHz > 0, 'omegaN/platformPitchFreqHz', 'must be > 0');
  const wp = 2 * Math.PI * o.platformPitchFreqHz;
  const m = o.margin ?? 1;
  return { ok: o.omegaN * m < wp, ratio: wp / o.omegaN, omegaNMax: wp / m };
}

/**
 * Closed-loop simulation of the 1-DOF rigid drivetrain with pitch actuator and the baseline controller:
 *   I dOmega/dt = Q_aero(V, Omega, theta) - N (T_gen + T_brake),  Q_aero from the Cq surface,
 * explicit mid-point integration of the rotor equation with the controller sampled once per step.
 * Limits: rigid shaft and tower, rotor-effective wind, quasi-steady aerodynamics from the surfaces.
 * @param {object} o
 * @param {object} o.surface cpSurface result
 * @param {object} o.cfg controller configuration (controllerFromRotor)
 * @param {number} o.inertia drivetrain inertia on the low-speed shaft [kg m2]
 * @param {number|number[]|function(number):number} o.wind rotor-effective wind speed [m/s]
 * @param {number} o.dt [s] @param {number} o.duration [s]
 * @param {number} [o.rho=1.225]
 * @param {number} [o.omega0] initial rotor speed [rad/s]; default rated
 * @param {number} [o.pitch0Deg] initial pitch [deg]; default fine pitch
 * @param {'parked'|'startup'|'production'} [o.mode='production'] initial controller mode
 * @param {object} [o.actuator] pitchActuatorStep constants (omegaAct, zetaAct)
 * @param {function(number):object} [o.events] t [s] -> { command?, gridLoss? } passed to the controller
 * @returns {{t:Float64Array, wind:Float64Array, rotorSpeed:Float64Array, pitchDeg:Float64Array, genTorque:Float64Array,
 *   power:Float64Array, thrust:Float64Array, aeroTorque:Float64Array, mode:string[]}} rotor speed [rad/s], generator
 *   torque [N m], electrical power [W], thrust [N], aerodynamic torque [N m]
 */
export function simulateClosedLoop(o) {
  need(o && o.dt > 0 && o.duration > 0, 'dt/duration', 'must be > 0');
  need(o.inertia > 0, 'inertia', 'must be > 0');
  const cfg = o.cfg; const sf = o.surface;
  const rho = o.rho ?? 1.225;
  const N = cfg.gearboxRatio ?? 1;
  const eta = cfg.efficiency ?? 1;
  const n = Math.round(o.duration / o.dt);
  const mode0 = o.mode ?? 'production';
  let om = o.omega0 ?? cfg.ratedGenSpeed / N;
  const act = { theta: (o.pitch0Deg !== undefined ? o.pitch0Deg * D2R : (mode0 === 'parked' ? cfg.pitchMax : cfg.pitchMin)), rate: 0 };
  const ctrl = createController(cfg, { mode: mode0, pitch: act.theta, genSpeed: om * N });
  const windAt = (k, t) => (typeof o.wind === 'number' ? o.wind : (typeof o.wind === 'function' ? o.wind(t) : o.wind[Math.min(k, o.wind.length - 1)]));
  const out = { t: new Float64Array(n + 1), wind: new Float64Array(n + 1), rotorSpeed: new Float64Array(n + 1), pitchDeg: new Float64Array(n + 1),
    genTorque: new Float64Array(n + 1), power: new Float64Array(n + 1), thrust: new Float64Array(n + 1), aeroTorque: new Float64Array(n + 1), mode: [] };
  const actArgs = { omegaAct: 2 * Math.PI, zetaAct: 0.7, ...(o.actuator || {}), cmd: 0, dt: o.dt, rateLimit: Math.max(cfg.pitchRate, cfg.estopPitchRate ?? 0), min: Math.min(0, cfg.pitchMin), max: cfg.pitchMax };
  const inp = { dt: o.dt, genSpeed: 0, pitch: 0, command: undefined, gridLoss: false, windSpeed: undefined };
  for (let k = 0; k <= n; k++) {
    const t = k * o.dt;
    const V = windAt(k, t);
    const ev = o.events ? (o.events(t) || {}) : {};
    inp.genSpeed = om * N; inp.pitch = act.theta; inp.command = ev.command; inp.gridLoss = !!ev.gridLoss;
    const c = ctrl.step(inp);
    const a0 = aeroPower(sf, rho, V, Math.max(om, 1e-4), act.theta / D2R);
    out.t[k] = t; out.wind[k] = V; out.rotorSpeed[k] = om; out.pitchDeg[k] = act.theta / D2R; out.genTorque[k] = c.genTorque;
    out.power[k] = eta * c.genTorque * om * N; out.thrust[k] = a0.thrust; out.aeroTorque[k] = a0.torque; out.mode.push(c.mode);
    if (k === n) break;
    const load = N * (c.genTorque + (om > 0 ? c.brakeTorque : 0));
    const k1 = (a0.torque - load) / o.inertia;
    const omMid = Math.max(0, om + 0.5 * o.dt * k1);
    const k2 = (aeroPower(sf, rho, V, Math.max(omMid, 1e-4), act.theta / D2R).torque - load) / o.inertia;
    om = Math.max(0, om + o.dt * k2);
    actArgs.cmd = c.pitchCmd;
    pitchActuatorStep(act, actArgs);
  }
  return out;
}
