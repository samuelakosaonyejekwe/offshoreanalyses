// Tidal-stream resource and turbine performance: harmonic synthesis, velocity statistics and profiles, turbine
// power/thrust with an injected rotor model or an actuator-disc fallback, cavitation, channel blockage,
// array wakes, bidirectional yaw loss, annual energy, extreme loads, flux resource and a 0-D tidal-range model.
//
// Units: SI; velocities m/s, time s (unless a name ends in Hours), power W, energy Wh where the name ends in Wh.

import { farmFlow } from './wake.js';

const G = 9.80665;
const DEG = Math.PI / 180;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Angular speeds of the principal tidal constituents in degrees per mean solar hour
 * (Doodson/Schureman astronomical speeds as tabulated by NOAA CO-OPS and the IHO constituent list).
 */
export const TIDAL_SPEEDS_DEG_PER_HOUR = Object.freeze({
  M2: 28.9841042,
  S2: 30.0,
  N2: 28.4397295,
  K2: 30.0821373,
  K1: 15.0410686,
  O1: 13.9430356,
  P1: 14.9589314,
  Q1: 13.3986609,
  M4: 57.9682084,
});

/**
 * Period of a tidal constituent, T = 360 / speed (hours).
 * @param {string} name constituent name (M2, S2, N2, K2, K1, O1, P1, Q1, M4)
 * @returns {{periodHours:number, period:number, omega:number}} period (h, s) and angular frequency (rad/s)
 */
export function constituentPeriod(name) {
  const sp = TIDAL_SPEEDS_DEG_PER_HOUR[name];
  need(sp !== undefined, `constituent '${name}' is not known`);
  const periodHours = 360 / sp;
  return { periodHours, period: periodHours * 3600, omega: (sp * DEG) / 3600 };
}

/**
 * Harmonic synthesis of the along-stream tidal velocity,
 *   u(t) = u0 + Σ_i a_i cos(ω_i t − φ_i),
 * with exact constituent speeds. Nodal corrections (f, u) and equilibrium arguments are not applied: phases are
 * relative to the chosen time origin, which is adequate for statistics over a year.
 * @param {object} o
 * @param {{name:string, amplitude:number, phaseDeg?:number}[]} o.constituents amplitude (m/s), phase lag (deg)
 * @param {number[]|Float64Array} [o.t] times (s); default a uniform series from `duration` and `dt`
 * @param {number} [o.duration=31536000] series length (s) when t is not given
 * @param {number} [o.dt=600] time step (s) when t is not given
 * @param {number} [o.mean=0] residual (mean) current u0 (m/s)
 * @returns {{t:Float64Array, u:Float64Array}} time (s) and signed velocity (m/s, positive = flood)
 */
export function synthesizeTide({ constituents, t, duration = 365 * 86400, dt = 600, mean = 0 }) {
  need(Array.isArray(constituents) && constituents.length > 0, 'constituents must be a non-empty array');
  let tt;
  if (t) tt = Float64Array.from(t);
  else {
    finitePos(dt, 'dt');
    finitePos(duration, 'duration');
    const n = Math.floor(duration / dt) + 1;
    need(n <= 5e6, 'duration/dt gives too many samples');
    tt = new Float64Array(n);
    for (let i = 0; i < n; i++) tt[i] = i * dt;
  }
  const u = new Float64Array(tt.length).fill(mean);
  for (const c of constituents) {
    const { omega } = constituentPeriod(c.name);
    need(Number.isFinite(c.amplitude) && c.amplitude >= 0, `amplitude of ${c.name} must be >= 0`);
    const ph = (c.phaseDeg ?? 0) * DEG;
    for (let i = 0; i < tt.length; i++) u[i] += c.amplitude * Math.cos(omega * tt[i] - ph);
  }
  return { t: tt, u };
}

/**
 * Two-constituent (M2 + S2) estimator from mean spring and mean neap peak velocities:
 *   a_M2 = (V_spring + V_neap)/2,  a_S2 = (V_spring − V_neap)/2,
 * spring–neap period 360 / (30 − 28.9841042) h = 14.765 days, and the slowly varying peak envelope
 *   U_peak(t) = √(a_M2² + a_S2² + 2 a_M2 a_S2 cos((ω_S2 − ω_M2) t)).
 * @param {object} o
 * @param {number} o.springPeak mean spring peak velocity (m/s)
 * @param {number} o.neapPeak mean neap peak velocity (m/s)
 * @returns {{constituents:{name:string, amplitude:number, phaseDeg:number}[], springNeapPeriod:number,
 *   envelope:(t:number)=>number}} constituents for `synthesizeTide` (in phase at t = 0, i.e. springs), period (s)
 */
export function springNeap({ springPeak, neapPeak }) {
  finitePos(springPeak, 'springPeak');
  need(Number.isFinite(neapPeak) && neapPeak >= 0 && neapPeak <= springPeak, 'neapPeak must be within [0, springPeak]');
  const a1 = 0.5 * (springPeak + neapPeak);
  const a2 = 0.5 * (springPeak - neapPeak);
  const dw = constituentPeriod('S2').omega - constituentPeriod('M2').omega;
  return {
    constituents: [
      { name: 'M2', amplitude: a1, phaseDeg: 0 },
      { name: 'S2', amplitude: a2, phaseDeg: 0 },
    ],
    springNeapPeriod: (2 * Math.PI) / dw,
    envelope: (t) => Math.sqrt(a1 * a1 + a2 * a2 + 2 * a1 * a2 * Math.cos(dw * t)),
  };
}

/**
 * Speed statistics of a velocity series: exceedance curve and histogram of |u|.
 * @param {number[]|Float64Array} u velocity samples (m/s), equally spaced in time
 * @param {{binWidth?:number}} [opts] histogram bin width (m/s), default 0.1
 * @returns {{speed:Float64Array, exceedance:Float64Array, bins:{u:number, prob:number}[], mean:number,
 *   max:number, cubeRootMeanCube:number}} sorted speeds (descending) with the fraction of time each is exceeded
 */
export function velocityDistribution(u, opts = {}) {
  need(u && u.length > 0, 'u must be a non-empty series');
  const bw = opts.binWidth ?? 0.1;
  finitePos(bw, 'binWidth');
  const n = u.length;
  const speed = new Float64Array(n);
  let s = 0;
  let s3 = 0;
  for (let i = 0; i < n; i++) {
    speed[i] = Math.abs(u[i]);
    s += speed[i];
    s3 += speed[i] ** 3;
  }
  speed.sort();
  speed.reverse();
  const exceedance = new Float64Array(n);
  for (let i = 0; i < n; i++) exceedance[i] = i / n;
  const nb = Math.floor(speed[0] / bw) + 1;
  const counts = new Float64Array(nb);
  for (let i = 0; i < n; i++) counts[Math.min(nb - 1, Math.floor(speed[i] / bw))]++;
  const bins = Array.from(counts, (c, i) => ({ u: (i + 0.5) * bw, prob: c / n }));
  return { speed, exceedance, bins, mean: s / n, max: speed[0], cubeRootMeanCube: Math.cbrt(s3 / n) };
}

/**
 * Vertical profile of a tidal current.
 *  - 'power': u(z) = uRef (z/zRef)^(1/n) with n = 7 (or 10); for zRef = depth, the depth mean is n/(n+1) uRef.
 *  - 'log':   u(z) = uRef ln(z/z0) / ln(zRef/z0) (rough-bed logarithmic layer).
 * z is height above the seabed. Reference: Soulsby (1997) Dynamics of Marine Sands (current profiles).
 * @param {object} o
 * @param {'power'|'log'} [o.model='power']
 * @param {number} o.uRef velocity at the reference height (m/s)
 * @param {number} o.zRef reference height above seabed (m)
 * @param {number} [o.n=7] power-law denominator (-)
 * @param {number} [o.z0] bed roughness length (m), required for 'log'
 * @returns {(z:number)=>number} velocity (m/s) at height z (m) above the seabed
 */
export function currentProfile({ model = 'power', uRef, zRef, n = 7, z0 }) {
  need(Number.isFinite(uRef), 'uRef must be finite');
  finitePos(zRef, 'zRef');
  if (model === 'power') {
    finitePos(n, 'n');
    return (z) => (z <= 0 ? 0 : uRef * Math.pow(z / zRef, 1 / n));
  }
  need(model === 'log', "model must be 'power' or 'log'");
  finitePos(z0, 'z0');
  need(z0 < zRef, 'z0 must be below zRef');
  const L = Math.log(zRef / z0);
  return (z) => (z <= z0 ? 0 : (uRef * Math.log(z / z0)) / L);
}

/**
 * Rotor-disc averages of a sheared current: area-mean velocity and the power-weighted velocity
 * (cube root of the area-mean of u³), by horizontal strips.
 * @param {object} o
 * @param {(z:number)=>number} o.profile velocity (m/s) at height z (m) above seabed
 * @param {number} o.hubHeight hub height above seabed (m)
 * @param {number} o.D rotor diameter (m)
 * @param {number} [o.nSeg=60] strips
 * @returns {{mean:number, cubeMean:number, powerVelocity:number, hub:number}} mean (m/s), mean of u³ (m³/s³),
 *   its cube root (m/s), hub velocity (m/s)
 */
export function rotorAveragedVelocity({ profile, hubHeight, D, nSeg = 60 }) {
  finitePos(D, 'D');
  need(hubHeight >= D / 2, 'hubHeight must be at least the rotor radius');
  need(typeof profile === 'function', 'profile must be a function');
  const R = D / 2;
  const below = (t) => Math.PI / 2 + Math.asin(t) + t * Math.sqrt(Math.max(0, 1 - t * t));
  let s1 = 0;
  let s3 = 0;
  for (let i = 0; i < nSeg; i++) {
    const t0 = -1 + (2 * i) / nSeg;
    const t1 = -1 + (2 * (i + 1)) / nSeg;
    const w = (below(t1) - below(t0)) / Math.PI;
    const v = profile(hubHeight + R * 0.5 * (t0 + t1));
    s1 += v * w;
    s3 += v * v * v * w;
  }
  return { mean: s1, cubeMean: s3, powerVelocity: Math.cbrt(s3), hub: profile(hubHeight) };
}

function inductionFromCp(Cp) {
  if (Cp <= 0) return 0;
  if (Cp >= 16 / 27) return 1 / 3;
  let lo = 0;
  let hi = 1 / 3;
  for (let i = 0; i < 60; i++) {
    const a = 0.5 * (lo + hi);
    if (4 * a * (1 - a) * (1 - a) < Cp) lo = a;
    else hi = a;
  }
  return 0.5 * (lo + hi);
}

/**
 * Ideal actuator-disc rotor performance (fallback when no rotor model is injected):
 * Cp = η 4 a (1 − a)², Ct = 4 a (1 − a). a = 1/3 is the Lanchester–Betz optimum (Cp = 16/27 η, Ct = 8/9).
 * @param {{a?:number, efficiency?:number}} [o] axial induction (-), hydrodynamic-to-electrical efficiency (-)
 * @returns {(U:number)=>{Cp:number, Ct:number}} rotor callback usable in `tidalPowerCurve`
 */
export function actuatorRotor({ a = 1 / 3, efficiency = 0.75 } = {}) {
  need(a > 0 && a < 0.5, 'a must be within (0, 0.5)');
  need(efficiency > 0 && efficiency <= 1, 'efficiency must be within (0, 1]');
  const Cp = efficiency * 4 * a * (1 - a) * (1 - a);
  const Ct = 4 * a * (1 - a);
  return () => ({ Cp, Ct });
}

/**
 * Adapter turning a blade-element-momentum solver into a rotor callback at a fixed tip-speed ratio and pitch
 * (variable-speed operation below rated). The solver must accept { rotor, U, omega, pitchDeg, rho } and return
 * { Cp, Ct } (or { power, thrust }).
 * @param {object} o
 * @param {Function} o.solve BEM solver
 * @param {object} o.rotor rotor definition passed through to the solver
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.tsr design tip-speed ratio λ = Ω R / U (-)
 * @param {number} [o.pitchDeg=0] blade pitch (deg)
 * @param {number} [o.rho=1025] water density (kg/m³)
 * @returns {(U:number)=>{Cp:number, Ct:number}}
 */
export function rotorFromSolver({ solve, rotor, D, tsr, pitchDeg = 0, rho = 1025 }) {
  need(typeof solve === 'function', 'solve must be a function');
  finitePos(D, 'D');
  finitePos(tsr, 'tsr');
  const A = (Math.PI * D * D) / 4;
  return (U) => {
    const r = solve({ rotor, U, omega: (tsr * U) / (D / 2), pitchDeg, rho });
    const Cp = Number.isFinite(r.Cp) ? r.Cp : r.power / (0.5 * rho * A * U ** 3);
    const Ct = Number.isFinite(r.Ct) ? r.Ct : r.thrust / (0.5 * rho * A * U * U);
    return { Cp, Ct };
  };
}

/**
 * Tidal-turbine power and thrust curve. Below rated the rotor callback supplies Cp and Ct; power is capped at
 * the rated value; above rated the required Cp = P_rated/(½ ρ A U³) is converted to a thrust coefficient with the
 * actuator-disc relation (pitch-to-feather shedding: a solves 4a(1 − a)² = Cp/η, Ct = 4a(1 − a)), limited to the
 * below-rated Ct. Below cut-in the rotor is taken as parked with zero power and `parkedCt`.
 * @param {object} o
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.ratedPower rated electrical power (W)
 * @param {(U:number)=>{Cp:number, Ct:number}} [o.rotor] rotor performance callback; default `actuatorRotor()`
 * @param {number} [o.rho=1025] sea-water density (kg/m³)
 * @param {number} [o.cutIn=0.7] cut-in speed (m/s)
 * @param {number} [o.uMax=5] highest tabulated speed (m/s)
 * @param {number} [o.step=0.1] table step (m/s)
 * @param {number} [o.efficiency=0.75] efficiency used to recover the hydrodynamic Cp above rated (-)
 * @param {number} [o.parkedCt=0] thrust coefficient below cut-in (-)
 * @returns {{u:number[], p:number[], cp:number[], ct:number[], thrust:number[], ratedSpeed:number, D:number,
 *   rho:number}} m/s, W, -, -, N; ratedSpeed is the first tabulated speed at rated power (NaN if never reached)
 */
export function tidalPowerCurve({ D, ratedPower, rotor, rho = 1025, cutIn = 0.7, uMax = 5, step = 0.1, efficiency = 0.75, parkedCt = 0 }) {
  finitePos(D, 'D');
  finitePos(ratedPower, 'ratedPower');
  finitePos(rho, 'rho');
  finitePos(step, 'step');
  need(cutIn >= 0 && uMax > cutIn, 'uMax must exceed cutIn >= 0');
  const perf = rotor ?? actuatorRotor({ efficiency });
  const A = (Math.PI * D * D) / 4;
  const out = { u: [], p: [], cp: [], ct: [], thrust: [], ratedSpeed: NaN, D, rho };
  const n = Math.round(uMax / step);
  for (let i = 0; i <= n; i++) {
    const U = i * step;
    let P = 0;
    let ct = parkedCt;
    let cp = 0;
    if (U >= cutIn && U > 0) {
      const r = perf(U);
      need(r && Number.isFinite(r.Cp) && Number.isFinite(r.Ct), 'rotor callback must return finite Cp and Ct');
      const q = 0.5 * rho * A * U ** 3;
      P = r.Cp * q;
      ct = r.Ct;
      if (P >= ratedPower) {
        P = ratedPower;
        if (Number.isNaN(out.ratedSpeed)) out.ratedSpeed = U;
        const a = inductionFromCp(P / q / efficiency);
        ct = Math.min(r.Ct, 4 * a * (1 - a));
      }
      cp = P / q;
    }
    out.u.push(U);
    out.p.push(P);
    out.cp.push(cp);
    out.ct.push(ct);
    out.thrust.push(0.5 * rho * A * ct * U * U);
  }
  return out;
}

/**
 * Blade-tip cavitation check at top dead centre.
 *   σ = (p_atm + ρ g h_tip − p_v) / (½ ρ W²),  W² = (U (1 − a))² + (Ω R)²,  h_tip = hub submergence − R.
 * Cavitation inception is predicted where σ < −Cp,min of the blade section (minimum pressure coefficient from
 * the caller's section data). Also returns the largest tip speed for which σ = −Cp,min × margin.
 * Reference: Batten et al. (2008) Renewable Energy 33:1085–1096; Carlton, Marine Propellers and Propulsion ch. 9.
 * @param {object} o
 * @param {number} o.U inflow speed (m/s)
 * @param {number} o.tipSpeed blade tip speed Ω R (m/s)
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.hubDepth hub submergence below the instantaneous free surface (m)
 * @param {number} o.cpMin minimum section pressure coefficient (negative, -)
 * @param {number} [o.rho=1025] sea-water density (kg/m³)
 * @param {number} [o.patm=101325] atmospheric pressure (Pa)
 * @param {number} [o.pv=1700] vapour pressure (Pa); ≈ 1230 Pa at 10 °C, 1700 Pa at 15 °C, 2340 Pa at 20 °C
 * @param {number} [o.a=0] axial induction at the tip (-); 0 is conservative
 * @param {number} [o.margin=1] required ratio σ / (−Cp,min) (-)
 * @returns {{sigma:number, sigmaRequired:number, cavitates:boolean, relativeSpeed:number, tipDepth:number,
 *   maxTipSpeed:number, utilisation:number}} utilisation = required/available cavitation number
 */
export function cavitationCheck({ U, tipSpeed, D, hubDepth, cpMin, rho = 1025, patm = 101325, pv = 1700, a = 0, margin = 1 }) {
  finitePos(D, 'D');
  finitePos(rho, 'rho');
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  need(Number.isFinite(tipSpeed) && tipSpeed >= 0, 'tipSpeed must be >= 0');
  need(cpMin < 0, 'cpMin must be negative');
  const tipDepth = hubDepth - D / 2;
  need(tipDepth >= 0, 'hubDepth must be at least the rotor radius (tip submerged)');
  const ua = U * (1 - a);
  const W2 = ua * ua + tipSpeed * tipSpeed;
  need(W2 > 0, 'U and tipSpeed must not both be zero');
  const dp = patm + rho * G * tipDepth - pv;
  const sigma = dp / (0.5 * rho * W2);
  const sigmaRequired = -cpMin * margin;
  const Wmax2 = dp / (0.5 * rho * sigmaRequired);
  return {
    sigma,
    sigmaRequired,
    cavitates: sigma < sigmaRequired,
    relativeSpeed: Math.sqrt(W2),
    tipDepth,
    maxTipSpeed: Wmax2 > ua * ua ? Math.sqrt(Wmax2 - ua * ua) : 0,
    utilisation: sigmaRequired / sigma,
  };
}

/**
 * Actuator disc in a channel of finite cross-section with a deformable free surface (linear momentum theory).
 * Notation (velocities scaled by the upstream speed u, depths by the upstream depth h): turbine-plane velocity
 * α2, wake velocity α4 at the pressure-equalisation section, bypass velocity β4, blockage ratio B = A_rotor /
 * A_channel, Froude number Fr = u/√(g h).
 *   bypass Bernoulli:     h4/h = 1 − (Fr²/2)(β4² − 1)
 *   continuity:           h4/h = B α2/α4 + (1 − B α2)/β4
 *   thrust (core Bernoulli): Ct = T/(½ ρ A u²) = β4² − α4²
 *   momentum (1 → 4):     (1 − (h4/h)²)/(2 Fr²) − (B/2) Ct = B α2 (α4 − 1) + (1 − B α2)(β4 − 1)
 *   power:                Cp = α2 Ct
 * For Fr → 0 the hydrostatic term tends to (β4² − 1)/2 and the system is the rigid-lid model of Garrett &
 * Cummins (2007), whose optimum is α4 = 1/3, Cp,max = (16/27)(1 − B)^−2. The far-downstream depth change after
 * mixing follows from the overall momentum balance:
 *   ½ x³ − (3/2) x² + (1 − Fr² + ½ Ct B Fr²) x − ½ Ct B Fr² = 0,  x = Δh/h.
 * References: Garrett & Cummins (2007) J. Fluid Mech. 588:243–251; Houlsby, Draper & Oldfield (2008) OUEL
 * Report 2296/08; Whelan, Graham & Peiró (2009) J. Fluid Mech. 624:281–291.
 * Valid for subcritical flow, uniform inflow, a disc spanning a fraction B of the section; no solution is
 * returned (RangeError) when the bypass flow would become critical.
 * @param {object} o
 * @param {number} o.B blockage ratio (-), 0 <= B < 1
 * @param {number} o.alpha4 wake velocity ratio u4/u (-), 0 < α4 < 1
 * @param {number} [o.Fr=0] upstream Froude number (-)
 * @returns {{alpha2:number, alpha4:number, beta4:number, Ct:number, Cp:number, h4OverH:number, dhOverH:number,
 *   CtLocal:number, CpLocal:number}} CtLocal, CpLocal are referred to the disc velocity (÷ α2², ÷ α2³)
 */
export function channelBlockage({ B, alpha4, Fr = 0 }) {
  need(Number.isFinite(B) && B >= 0 && B < 1, 'B must be within [0, 1)');
  need(alpha4 > 0 && alpha4 < 1, 'alpha4 must be within (0, 1)');
  need(Number.isFinite(Fr) && Fr >= 0 && Fr < 1, 'Fr must be within [0, 1)');
  const F2 = Fr * Fr;
  let alpha2;
  let beta4;
  if (B < 1e-12) {
    beta4 = 1;
    alpha2 = 0.5 * (1 + alpha4);
  } else {
    const h4 = (b) => 1 - 0.5 * F2 * (b * b - 1);
    const a2 = (b) => (h4(b) - 1 / b) / (B * (1 / alpha4 - 1 / b));
    const res = (b) => {
      const al = a2(b);
      const e = b * b - 1;
      return 0.5 * e * (1 - 0.25 * F2 * e) - 0.5 * B * (b * b - alpha4 * alpha4) - B * al * (alpha4 - 1) - (1 - B * al) * (b - 1);
    };
    let lo = 1 + 1e-10;
    let hi = NaN;
    let fl = res(lo);
    for (let d = 1e-8; d < 50; d *= 1.15) {
      const b = 1 + d;
      if (h4(b) <= 0) break;
      const f = res(b);
      if (fl < 0 !== f < 0) {
        hi = b;
        break;
      }
      lo = b;
      fl = f;
    }
    need(Number.isFinite(hi), 'no subcritical solution for this combination of B, Fr and alpha4');
    for (let i = 0; i < 100; i++) {
      const mid = 0.5 * (lo + hi);
      const f = res(mid);
      if (fl < 0 === f < 0) {
        lo = mid;
        fl = f;
      } else hi = mid;
    }
    beta4 = 0.5 * (lo + hi);
    alpha2 = a2(beta4);
    need(alpha2 > alpha4 - 1e-9 && alpha2 <= 1 + 1e-9, 'no physical solution (turbine-plane velocity out of range)');
  }
  const Ct = beta4 * beta4 - alpha4 * alpha4;
  const Cp = alpha2 * Ct;
  // downstream mixed depth change: smallest positive root of the cubic (Newton from x = 0)
  let x = 0;
  if (F2 > 0 && B > 0) {
    const c0 = 0.5 * Ct * B * F2;
    for (let i = 0; i < 60; i++) {
      const f = 0.5 * x ** 3 - 1.5 * x * x + (1 - F2 + c0) * x - c0;
      const df = 1.5 * x * x - 3 * x + (1 - F2 + c0);
      const dx = f / df;
      x -= dx;
      if (Math.abs(dx) < 1e-15) break;
    }
  }
  return {
    alpha2,
    alpha4,
    beta4,
    Ct,
    Cp,
    h4OverH: 1 - 0.5 * F2 * (beta4 * beta4 - 1),
    dhOverH: x,
    CtLocal: Ct / (alpha2 * alpha2),
    CpLocal: Cp / alpha2 ** 3,
  };
}

/**
 * Maximum power coefficient of an actuator disc in a blocked channel, by maximising `channelBlockage` over the
 * wake velocity ratio α4 (golden-section search). For Fr = 0 this reproduces Garrett & Cummins (2007):
 * Cp,max = (16/27)(1 − B)^−2 at α4 = 1/3.
 * @param {object} o
 * @param {number} o.B blockage ratio (-)
 * @param {number} [o.Fr=0] Froude number (-)
 * @returns {{CpMax:number, alpha4:number, Ct:number, alpha2:number, beta4:number, garrettCummins:number,
 *   dhOverH:number}} garrettCummins is the rigid-lid closed form for comparison
 */
export function channelBlockageLimit({ B, Fr = 0 }) {
  const f = (a4) => {
    try {
      return channelBlockage({ B, alpha4: a4, Fr }).Cp;
    } catch {
      return -Infinity;
    }
  };
  const gr = (Math.sqrt(5) - 1) / 2;
  let a = 0.02;
  let b = 0.98;
  let c = b - gr * (b - a);
  let d = a + gr * (b - a);
  let fc = f(c);
  let fd = f(d);
  for (let i = 0; i < 80; i++) {
    if (fc > fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - gr * (b - a);
      fc = f(c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + gr * (b - a);
      fd = f(d);
    }
  }
  const a4 = 0.5 * (a + b);
  const r = channelBlockage({ B, alpha4: a4, Fr });
  return { CpMax: r.Cp, alpha4: a4, Ct: r.Ct, alpha2: r.alpha2, beta4: r.beta4, garrettCummins: 16 / 27 / (1 - B) ** 2, dhOverH: r.dhOverH };
}

/**
 * Blockage correction of rotor coefficients measured or computed in a confined channel to unbounded flow by the
 * equal-thrust, equal-disc-velocity argument (Barnsley & Wellicome 1990; Bahaj et al. 2007): the unbounded
 * free-stream speed U_F that gives the same thrust at the same disc velocity satisfies
 *   U_F/U = α2 + Ct/(4 α2)    (from T = 2 ρ A u_t (U_F − u_t)),
 * and Ct_free = Ct (U/U_F)², Cp_free = Cp (U/U_F)³, λ_free = λ (U/U_F).
 * The channel state (α2) is found from `channelBlockage` by matching the given Ct.
 * @param {object} o
 * @param {number} o.Ct thrust coefficient referred to the channel upstream speed (-)
 * @param {number} o.Cp power coefficient referred to the channel upstream speed (-)
 * @param {number} o.B blockage ratio (-)
 * @param {number} [o.Fr=0] Froude number (-)
 * @returns {{speedRatio:number, CtFree:number, CpFree:number, alpha2:number, alpha4:number}} speedRatio = U/U_F
 */
export function blockageCorrection({ Ct, Cp, B, Fr = 0 }) {
  finitePos(Ct, 'Ct');
  need(Number.isFinite(Cp), 'Cp must be finite');
  // Ct(α4) decreases monotonically from its maximum at α4 → 0 to 0 at α4 → 1
  let lo = 0.02;
  let hi = 0.999;
  const g = (a4) => channelBlockage({ B, alpha4: a4, Fr }).Ct - Ct;
  need(g(lo) > 0, 'Ct is too large for the momentum model at this blockage');
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) > 0) lo = mid;
    else hi = mid;
  }
  const r = channelBlockage({ B, alpha4: 0.5 * (lo + hi), Fr });
  const uf = r.alpha2 + Ct / (4 * r.alpha2);
  const s = 1 / uf;
  return { speedRatio: s, CtFree: Ct * s * s, CpFree: Cp * s ** 3, alpha2: r.alpha2, alpha4: r.alpha4 };
}

/**
 * Tidal-array flow for one flow state using the Gaussian wake model of wake.js with tidal defaults.
 * `flowDirDeg` is the direction the current flows TOWARD (oceanographic convention, clockwise from north).
 * Turbine `hubHeight` is the height above the seabed; with `reflect` the seabed acts as a mirror. The free
 * surface and channel blockage are not modelled here (see `channelBlockage`).
 * Wake growth: by default from the ambient turbulence intensity with the wind-derived Niayifar & Porté-Agel fit;
 * supply `kStar` from site/flume data where available.
 * @param {object} o
 * @param {object[]} o.turbines as for `farmFlow` ({ x, y, D, hubHeight, power, ct })
 * @param {number} o.U undisturbed current speed at hub height (m/s)
 * @param {number} o.flowDirDeg direction of flow, toward (deg)
 * @param {number} [o.TI=0.1] ambient turbulence intensity (-), typically 0.08–0.15 at tidal sites
 * @param {object} [o.wake] further `farmFlow` options
 * @returns {object} result of `farmFlow`
 */
export function tidalArrayFlow({ turbines, U, flowDirDeg, TI = 0.1, wake = {} }) {
  need(Number.isFinite(flowDirDeg), 'flowDirDeg must be finite');
  return farmFlow({ model: 'gaussian', reflect: true, ...wake, turbines, U, TI, wdDeg: flowDirDeg + 180 });
}

/**
 * Fixed-heading (non-yawing, bidirectional) turbine in a tide whose flood and ebb directions are not exactly
 * opposed. The power in each phase is scaled by cos^p of the misalignment between the flow and the rotor axis
 * (p = 3 from the axial-velocity argument; measured values for yawed rotors are 2–3, hence an argument).
 * @param {object} o
 * @param {number} o.floodDirDeg flood flow direction, toward (deg)
 * @param {number} o.ebbDirDeg ebb flow direction, toward (deg)
 * @param {number} [o.headingDeg] rotor-axis heading (deg); default the energy-optimal heading
 * @param {number} [o.floodEnergy=1] energy available on the flood with perfect alignment (any unit)
 * @param {number} [o.ebbEnergy=1] energy available on the ebb with perfect alignment (same unit)
 * @param {number} [o.exponent=3] cosine exponent p (-)
 * @returns {{headingDeg:number, floodMisalignmentDeg:number, ebbMisalignmentDeg:number, floodFactor:number,
 *   ebbFactor:number, energy:number, loss:number, asymmetryDeg:number}} loss = 1 − energy/(flood + ebb);
 *   asymmetryDeg = departure of the two directions from 180° opposition
 */
export function bidirectionalYield({ floodDirDeg, ebbDirDeg, headingDeg, floodEnergy = 1, ebbEnergy = 1, exponent = 3 }) {
  need(Number.isFinite(floodDirDeg) && Number.isFinite(ebbDirDeg), 'flow directions must be finite');
  need(floodEnergy >= 0 && ebbEnergy >= 0 && floodEnergy + ebbEnergy > 0, 'energies must be >= 0 and not both zero');
  finitePos(exponent, 'exponent');
  const wrap = (a) => ((((a + 180) % 360) + 360) % 360) - 180;
  const evalH = (h) => {
    const gf = wrap(floodDirDeg - h);
    const ge = wrap(ebbDirDeg - (h + 180));
    const ff = Math.abs(gf) < 90 ? Math.cos(gf * DEG) ** exponent : 0;
    const fe = Math.abs(ge) < 90 ? Math.cos(ge * DEG) ** exponent : 0;
    return { gf, ge, ff, fe, e: floodEnergy * ff + ebbEnergy * fe };
  };
  let h = headingDeg;
  if (h === undefined) {
    let best = -1;
    const h0 = floodDirDeg;
    const span = wrap(ebbDirDeg - 180 - floodDirDeg);
    for (let i = 0; i <= 2000; i++) {
      const cand = h0 + (span * i) / 2000;
      const e = evalH(cand).e;
      if (e > best) {
        best = e;
        h = cand;
      }
    }
  }
  const r = evalH(h);
  return {
    headingDeg: h,
    floodMisalignmentDeg: r.gf,
    ebbMisalignmentDeg: r.ge,
    floodFactor: r.ff,
    ebbFactor: r.fe,
    energy: r.e,
    loss: 1 - r.e / (floodEnergy + ebbEnergy),
    asymmetryDeg: wrap(ebbDirDeg - 180 - floodDirDeg),
  };
}

function interpTable(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return x < xs[0] ? 0 : ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid;
    else hi = mid;
  }
  return ys[lo] + ((x - xs[lo]) / (xs[hi] - xs[lo])) * (ys[hi] - ys[lo]);
}

/**
 * Annual energy of a tidal turbine from a velocity time series or a speed distribution.
 * Speeds above the last tabulated value use the last power (rated).
 * @param {object} o
 * @param {{u:number[], p:number[]}} o.curve power curve (m/s, W)
 * @param {number[]|Float64Array} [o.speeds] velocity series (m/s, sign ignored), equally spaced
 * @param {{u:number, prob:number}[]} [o.distribution] speed bins with probabilities (alternative to `speeds`)
 * @param {number} [o.availability=1] technical availability (-)
 * @param {number} [o.velocityFactor=1] multiplier on every speed (e.g. rotor-average / reference ratio) (-)
 * @param {number} [o.hoursPerYear=8760]
 * @param {number} [o.ratedPower] rated power (W); default the maximum of the curve
 * @returns {{aepWh:number, meanPower:number, capacityFactor:number, ratedFraction:number}} ratedFraction =
 *   fraction of time at rated power (-)
 */
export function tidalEnergy({ curve, speeds, distribution, availability = 1, velocityFactor = 1, hoursPerYear = 8760, ratedPower }) {
  need(curve && curve.u && curve.p && curve.u.length > 1, 'curve must be a table { u: [], p: [] }');
  need(availability >= 0 && availability <= 1, 'availability must be within [0, 1]');
  need(speeds || distribution, 'either speeds or distribution must be given');
  const rated = ratedPower ?? Math.max(...curve.p);
  let mean = 0;
  let atRated = 0;
  if (speeds) {
    need(speeds.length > 0, 'speeds must be non-empty');
    for (let i = 0; i < speeds.length; i++) {
      const p = interpTable(curve.u, curve.p, Math.abs(speeds[i]) * velocityFactor);
      mean += p;
      if (p >= rated * (1 - 1e-9)) atRated++;
    }
    mean /= speeds.length;
    atRated /= speeds.length;
  } else {
    let ps = 0;
    for (const b of distribution) {
      need(b.prob >= 0, 'distribution probabilities must be >= 0');
      const p = interpTable(curve.u, curve.p, b.u * velocityFactor);
      mean += p * b.prob;
      ps += b.prob;
      if (p >= rated * (1 - 1e-9)) atRated += b.prob;
    }
    finitePos(ps, 'sum of distribution probabilities');
    mean /= ps;
    atRated /= ps;
  }
  mean *= availability;
  return { aepWh: mean * hoursPerYear, meanPower: mean, capacityFactor: mean / rated, ratedFraction: atRated };
}

/** Linear dispersion ω² = g k tanh(k h): wavenumber by Newton iteration. */
function waveNumber(omega, h) {
  const k0 = (omega * omega) / G;
  let k = k0 / Math.sqrt(Math.tanh(k0 * h));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(k * h);
    const f = G * k * th - omega * omega;
    const df = G * th + G * k * h * (1 - th * th);
    const dk = f / df;
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

function orbitalAmplitude({ H, T, depth }, zAboveBed) {
  const omega = (2 * Math.PI) / T;
  const k = waveNumber(omega, depth);
  const kh = k * depth;
  // cosh(k z)/sinh(k h) written with decaying exponentials to stay finite in deep water
  const ratio = (Math.exp(k * (zAboveBed - depth)) + Math.exp(-k * (zAboveBed + depth))) / (1 - Math.exp(-2 * kh));
  return { u: ((Math.PI * H) / T) * ratio, k, omega };
}

/**
 * Extreme and wave-induced cyclic thrust on a tidal rotor.
 * Design velocity U = spring peak + storm-surge current + wave orbital velocity amplitude at hub level (linear
 * wave theory, u_w = (π H/T) cosh(k z)/sinh(k h), z above bed), thrust T = ½ ρ A Ct U². The wave-induced thrust
 * range about the mean current U_c is ΔT = ½ ρ A Ct [(U_c + u_w)² − (U_c − u_w)²] = 2 ρ A Ct U_c u_w (u_w ≤ U_c),
 * with one cycle per wave period — the hook for wave–current fatigue counting (quasi-steady Ct, no dynamic
 * inflow or added-mass effects).
 * @param {object} o
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.Ct thrust coefficient in the governing state (-)
 * @param {number} o.springPeak maximum spring-tide velocity at hub height (m/s)
 * @param {number} [o.surgeCurrent=0] storm-surge / wind-driven current (m/s)
 * @param {{H:number, T:number, depth:number}} [o.wave] wave height (m), period (s), water depth (m)
 * @param {number} [o.hubHeight] hub height above seabed (m), required with `wave`
 * @param {number} [o.rho=1025] density (kg/m³)
 * @param {number} [o.duration=3600] duration for the cycle count (s)
 * @returns {{designVelocity:number, waveVelocity:number, thrust:number, meanThrust:number, thrustRange:number,
 *   cycles:number}} m/s, m/s, N, N, N, number of wave cycles in `duration`
 */
export function extremeThrust({ D, Ct, springPeak, surgeCurrent = 0, wave, hubHeight, rho = 1025, duration = 3600 }) {
  finitePos(D, 'D');
  need(Ct >= 0, 'Ct must be >= 0');
  need(springPeak >= 0 && surgeCurrent >= 0, 'velocities must be >= 0');
  const A = (Math.PI * D * D) / 4;
  let uw = 0;
  let cycles = 0;
  if (wave) {
    finitePos(wave.H, 'wave.H');
    finitePos(wave.T, 'wave.T');
    finitePos(wave.depth, 'wave.depth');
    need(hubHeight > 0 && hubHeight <= wave.depth, 'hubHeight must be within (0, wave.depth]');
    uw = orbitalAmplitude(wave, hubHeight).u;
    cycles = duration / wave.T;
  }
  const Uc = springPeak + surgeCurrent;
  const q = 0.5 * rho * A * Ct;
  const lowSq = Uc >= uw ? (Uc - uw) ** 2 : -((Uc - uw) ** 2);
  return {
    designVelocity: Uc + uw,
    waveVelocity: uw,
    thrust: q * (Uc + uw) ** 2,
    meanThrust: q * Uc * Uc,
    thrustRange: q * ((Uc + uw) ** 2 - lowSq),
    cycles,
  };
}

/**
 * Base shear and overturning moment of a tidal-turbine support column: rotor thrust at hub height plus Morison
 * drag on the column, f(z) = ½ ρ Cd D_c u(z)|u(z)|, integrated from the seabed to the top of the column.
 * With a wave, the drag uses current + orbital velocity at the crest phase and the inertia-force amplitude
 * (ρ Cm π D_c²/4 · ω u_w(z)) is reported separately because it is 90° out of phase with the drag maximum.
 * Reference: Morison et al. (1950); DNV-RP-C205 §6.
 * @param {object} o
 * @param {number} o.thrust rotor thrust (N)
 * @param {number} o.hubHeight hub height above seabed (m)
 * @param {number} o.columnDiameter column diameter (m)
 * @param {(z:number)=>number} o.profile current velocity (m/s) at height z above seabed
 * @param {number} [o.columnTop=hubHeight] top of the column above seabed (m)
 * @param {number} [o.Cd=1.0] drag coefficient (-); DNV-RP-C205 gives 0.65 (smooth) to 1.05 (rough) for cylinders
 * @param {number} [o.Cm=2.0] inertia coefficient (-)
 * @param {{H:number, T:number, depth:number}} [o.wave]
 * @param {number} [o.rho=1025]
 * @param {number} [o.nSeg=100] integration segments
 * @returns {{baseShear:number, baseMoment:number, columnDrag:number, columnDragMoment:number,
 *   inertiaShear:number, inertiaMoment:number}} N and N·m about the seabed
 */
export function supportLoads({ thrust, hubHeight, columnDiameter, profile, columnTop, Cd = 1.0, Cm = 2.0, wave, rho = 1025, nSeg = 100 }) {
  need(Number.isFinite(thrust), 'thrust must be finite');
  finitePos(hubHeight, 'hubHeight');
  finitePos(columnDiameter, 'columnDiameter');
  need(typeof profile === 'function', 'profile must be a function');
  const top = columnTop ?? hubHeight;
  finitePos(top, 'columnTop');
  const dz = top / nSeg;
  let Fd = 0;
  let Md = 0;
  let Fi = 0;
  let Mi = 0;
  for (let i = 0; i < nSeg; i++) {
    const z = (i + 0.5) * dz;
    let u = profile(z);
    if (wave) {
      const w = orbitalAmplitude(wave, Math.min(z, wave.depth));
      u += w.u;
      const fi = ((rho * Cm * Math.PI * columnDiameter ** 2) / 4) * w.omega * w.u * dz;
      Fi += fi;
      Mi += fi * z;
    }
    const f = 0.5 * rho * Cd * columnDiameter * u * Math.abs(u) * dz;
    Fd += f;
    Md += f * z;
  }
  return {
    baseShear: thrust + Fd,
    baseMoment: thrust * hubHeight + Md,
    columnDrag: Fd,
    columnDragMoment: Md,
    inertiaShear: Fi,
    inertiaMoment: Mi,
  };
}

/**
 * Kinetic power density statistics of a current series, P/A = ½ ρ |u|³ (W/m²).
 * @param {number[]|Float64Array} u velocity samples (m/s)
 * @param {number} [rho=1025] density (kg/m³)
 * @returns {{mean:number, max:number, meanSpeed:number, maxSpeed:number, cubeRootMeanCube:number}} mean and
 *   maximum power density (W/m²), speeds (m/s)
 */
export function kineticPowerDensity(u, rho = 1025) {
  need(u && u.length > 0, 'u must be a non-empty series');
  finitePos(rho, 'rho');
  let s = 0;
  let s3 = 0;
  let mx = 0;
  for (let i = 0; i < u.length; i++) {
    const a = Math.abs(u[i]);
    s += a;
    s3 += a * a * a;
    if (a > mx) mx = a;
  }
  const m3 = s3 / u.length;
  return { mean: 0.5 * rho * m3, max: 0.5 * rho * mx ** 3, meanSpeed: s / u.length, maxSpeed: mx, cubeRootMeanCube: Math.cbrt(m3) };
}

/**
 * Flux-method resource estimate: mean kinetic power through a channel cross-section times the significant
 * impact factor (SIF), the fraction that can be extracted without significant environmental or economic effect.
 * Reference: Black & Veatch (2005) Phase II UK Tidal Stream Energy Resource Assessment (SIF ≈ 0.1–0.2 for most
 * sites, up to 0.5 for sea lochs); the SIF is site-specific and is an argument without default.
 * @param {object} o
 * @param {number} o.meanPowerDensity mean kinetic power density over the section (W/m²)
 * @param {number} o.area channel cross-section area (m²)
 * @param {number} o.sif significant impact factor (-)
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{fluxPower:number, extractablePower:number, extractableWh:number}} W, W, Wh/yr
 */
export function fluxResource({ meanPowerDensity, area, sif, hoursPerYear = 8760 }) {
  need(meanPowerDensity >= 0, 'meanPowerDensity must be >= 0');
  finitePos(area, 'area');
  need(sif >= 0 && sif <= 1, 'sif must be within [0, 1]');
  const fluxPower = meanPowerDensity * area;
  return { fluxPower, extractablePower: fluxPower * sif, extractableWh: fluxPower * sif * hoursPerYear };
}

/**
 * Tidal-range (barrage/lagoon) ebb-generation scheme, 0-D basin model.
 * Sea level ζ(t) = (R/2) cos(2π t/T); basin of constant plan area A with level Z:
 *   holding   — gates closed until the head H = Z − ζ reaches `startHead`;
 *   generating — Q_t = Cd_t A_t √(2 g H) (orifice law), P = η ρ g Q_t H, capped at `ratedPower` by throttling Q_t;
 *                ends when H < `minHead`;
 *   filling   — once ζ > Z the sluices pass Q_s = Cd_s A_s √(2 g (ζ − Z)) until the tide turns.
 * dZ/dt = ∓Q/A, explicit time stepping. The theoretical maximum potential energy per tide released from high
 * water to low water is E = ½ ρ g A R². Reference: Prandle (1984) Adv. Water Resources 7:21–27; Baker (1991)
 * Tidal Power. A screening model: no basin hypsometry, no pumping, no two-way generation, no turbine hill chart.
 * @param {object} o
 * @param {number} o.basinArea basin plan area A (m²)
 * @param {number} o.range tidal range R (m)
 * @param {number} o.turbineArea total turbine flow area A_t (m²)
 * @param {number} o.sluiceArea total sluice area A_s (m²)
 * @param {number} [o.period=44714.16] tidal period (s), M2
 * @param {number} [o.cdTurbine=1.0] turbine discharge coefficient (-)
 * @param {number} [o.cdSluice=1.0] sluice discharge coefficient (-)
 * @param {number} [o.startHead] head to start generating (m); default R/2
 * @param {number} [o.minHead=1] minimum generating head (m)
 * @param {number} [o.efficiency=0.85] water-to-wire efficiency (-)
 * @param {number} [o.ratedPower=Infinity] installed capacity (W)
 * @param {number} [o.rho=1025]
 * @param {number} [o.dt=60] time step (s)
 * @param {number} [o.cycles=6] simulated tides; results are taken from the last one
 * @returns {{energyPerCycleWh:number, potentialPerCycleWh:number, efficiency:number, meanPower:number,
 *   peakPower:number, annualWh:number, generatingFraction:number, t:Float64Array, sea:Float64Array,
 *   basin:Float64Array, power:Float64Array}} series cover the last tide (s, m, m, W)
 */
export function tidalRangeEbb({
  basinArea,
  range,
  turbineArea,
  sluiceArea,
  period = 44714.16,
  cdTurbine = 1.0,
  cdSluice = 1.0,
  startHead,
  minHead = 1,
  efficiency = 0.85,
  ratedPower = Infinity,
  rho = 1025,
  dt = 60,
  cycles = 6,
}) {
  finitePos(basinArea, 'basinArea');
  finitePos(range, 'range');
  finitePos(turbineArea, 'turbineArea');
  finitePos(sluiceArea, 'sluiceArea');
  finitePos(period, 'period');
  finitePos(dt, 'dt');
  need(Number.isInteger(cycles) && cycles >= 1, 'cycles must be an integer >= 1');
  need(efficiency > 0 && efficiency <= 1, 'efficiency must be within (0, 1]');
  const hs = startHead ?? range / 2;
  need(hs > minHead && minHead >= 0, 'startHead must exceed minHead >= 0');
  const nPer = Math.round(period / dt);
  const h = period / nPer;
  const w = (2 * Math.PI) / period;
  let Z = range / 2;
  let mode = 'hold';
  const t = new Float64Array(nPer);
  const sea = new Float64Array(nPer);
  const basin = new Float64Array(nPer);
  const power = new Float64Array(nPer);
  let e = 0;
  let peak = 0;
  let gen = 0;
  for (let c = 0; c < cycles; c++) {
    const last = c === cycles - 1;
    for (let i = 0; i < nPer; i++) {
      const time = (c * nPer + i) * h;
      const z = 0.5 * range * Math.cos(w * time);
      const falling = Math.sin(w * time) > 0;
      const H = Z - z;
      if (mode === 'hold' && H >= hs) mode = 'generate';
      else if (mode === 'generate' && H < minHead) mode = 'wait';
      else if (mode === 'wait' && z > Z) mode = 'fill';
      else if (mode === 'fill' && (z <= Z && falling)) mode = 'hold';
      let P = 0;
      if (mode === 'generate') {
        let Q = cdTurbine * turbineArea * Math.sqrt(2 * G * H);
        P = efficiency * rho * G * Q * H;
        if (P > ratedPower) {
          P = ratedPower;
          Q = P / (efficiency * rho * G * H);
        }
        Z -= (Q * h) / basinArea;
      } else if (mode === 'fill' && z > Z) {
        const Q = cdSluice * sluiceArea * Math.sqrt(2 * G * (z - Z));
        Z = Math.min(z, Z + (Q * h) / basinArea);
      }
      if (last) {
        t[i] = i * h;
        sea[i] = z;
        basin[i] = Z;
        power[i] = P;
        e += P * h;
        if (P > peak) peak = P;
        if (P > 0) gen++;
      }
    }
  }
  const potential = 0.5 * rho * G * basinArea * range * range;
  return {
    energyPerCycleWh: e / 3600,
    potentialPerCycleWh: potential / 3600,
    efficiency: e / potential,
    meanPower: e / period,
    peakPower: peak,
    annualWh: (e / 3600) * ((8760 * 3600) / period),
    generatingFraction: gen / nPer,
    t,
    sea,
    basin,
    power,
  };
}
