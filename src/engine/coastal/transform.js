// Nearshore wave transformation: shoaling, refraction, breaking, surf-zone energy balance, set-up,
// diffraction, reflection and wave-current interaction. SI units, radians unless the name ends in Deg.

import { brent } from '../math/roots.js';

const G = 9.80665;
const TWO_PI = 2 * Math.PI;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/** Linear dispersion ω² = g k tanh(kh) solved by Newton iteration on kh (private helper). */
function waveNumber(T, h, g) {
  const omega = TWO_PI / T;
  const k0h = (omega * omega * h) / g;
  // Starting guess that is exact in both the deep and shallow limits.
  let kh = k0h / Math.sqrt(Math.tanh(k0h));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(kh);
    const f = kh * th - k0h;
    const df = th + kh * (1 - th * th);
    const d = f / df;
    kh -= d;
    if (Math.abs(d) < 1e-14 * kh) break;
  }
  return kh / h;
}

/**
 * Linear (Airy) wave properties at a given depth.
 * Reference: Dean & Dalrymple (1991) "Water Wave Mechanics for Engineers and Scientists", ch. 3–4.
 * Validity: small-amplitude waves, any relative depth.
 * @param {object} o
 * @param {number} o.T wave period [s]
 * @param {number} o.h water depth [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{k:number,L:number,C:number,n:number,Cg:number,L0:number,C0:number,omega:number}}
 *   k [rad/m], L [m], C [m/s], n = Cg/C [-], Cg [m/s], deep-water L0 [m] and C0 [m/s], omega [rad/s]
 */
export function linearWave({ T, h, g = G }) {
  pos(T, 'T');
  pos(h, 'h');
  const k = waveNumber(T, h, g);
  const kh = k * h;
  const omega = TWO_PI / T;
  const C = omega / k;
  const n = kh > 350 ? 0.5 : 0.5 * (1 + (2 * kh) / Math.sinh(2 * kh));
  const L0 = (g * T * T) / TWO_PI;
  return { k, L: TWO_PI / k, C, n, Cg: n * C, L0, C0: L0 / T, omega };
}

/**
 * Linear shoaling coefficient Ks = H/H0 = sqrt(Cg0/Cg) for normal incidence from deep water.
 * Reference: conservation of energy flux, Dean & Dalrymple (1991) §4.8.
 * Validity: linear waves, no dissipation, mild slope.
 * @param {object} o
 * @param {number} o.T period [s]
 * @param {number} o.h depth [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {number} Ks [-]
 */
export function shoalingCoefficient({ T, h, g = G }) {
  const w = linearWave({ T, h, g });
  return Math.sqrt((0.5 * w.C0) / w.Cg);
}

/**
 * Refraction over straight, parallel depth contours by Snell's law, from deep water.
 * sin θ / C = sin θ0 / C0,  Kr = sqrt(cos θ0 / cos θ).
 * Reference: Dean & Dalrymple (1991) §4.8.2.
 * Validity: straight parallel contours, linear waves, |θ0| < 90°.
 * @param {object} o
 * @param {number} o.T period [s]
 * @param {number} o.h depth [m]
 * @param {number} o.theta0Deg deep-water angle between wave crests and contours (= ray to shore-normal) [deg]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{thetaDeg:number, theta:number, Kr:number}} local angle [deg], [rad] and refraction coefficient [-]
 */
export function refractionSnell({ T, h, theta0Deg, g = G }) {
  need(Number.isFinite(theta0Deg) && Math.abs(theta0Deg) < 90, 'theta0Deg must lie in (-90, 90)');
  const w = linearWave({ T, h, g });
  const t0 = (theta0Deg * Math.PI) / 180;
  const theta = Math.asin((Math.sin(t0) * w.C) / w.C0);
  return { theta, thetaDeg: (theta * 180) / Math.PI, Kr: Math.sqrt(Math.cos(t0) / Math.cos(theta)) };
}

/**
 * Combined linear shoaling and refraction of a deep-water wave along a profile of still-water depths
 * (straight parallel contours). Optionally caps the height with a depth-limited breaker index.
 * Reference: Dean & Dalrymple (1991) §4.8; cap is the McCowan-type criterion H <= gammaB h.
 * Validity: linear, no dissipation seaward of breaking; depths > 0.
 * @param {object} o
 * @param {number} o.H0 deep-water wave height [m]
 * @param {number} o.T period [s]
 * @param {number} [o.theta0Deg] deep-water incidence angle from shore-normal [deg]
 * @param {ArrayLike<number>} o.depths still-water depths [m]
 * @param {number} [o.gammaB] if given, H is limited to gammaB*h [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Ks:Float64Array,Kr:Float64Array,H:Float64Array,thetaDeg:Float64Array,L:Float64Array,broken:Uint8Array}}
 *   coefficients [-], heights [m], angles [deg], wavelengths [m], flag 1 where the breaker cap is active
 */
export function transformProfile({ H0, T, theta0Deg = 0, depths, gammaB, g = G }) {
  pos(H0, 'H0');
  need(depths && depths.length > 0, 'depths must be a non-empty array');
  const n = depths.length;
  const Ks = new Float64Array(n);
  const Kr = new Float64Array(n);
  const H = new Float64Array(n);
  const thetaDeg = new Float64Array(n);
  const L = new Float64Array(n);
  const broken = new Uint8Array(n);
  const t0 = (theta0Deg * Math.PI) / 180;
  need(Math.abs(theta0Deg) < 90, 'theta0Deg must lie in (-90, 90)');
  for (let i = 0; i < n; i++) {
    const w = linearWave({ T, h: depths[i], g });
    const th = Math.asin((Math.sin(t0) * w.C) / w.C0);
    Ks[i] = Math.sqrt((0.5 * w.C0) / w.Cg);
    Kr[i] = Math.sqrt(Math.cos(t0) / Math.cos(th));
    thetaDeg[i] = (th * 180) / Math.PI;
    L[i] = w.L;
    let Hi = H0 * Ks[i] * Kr[i];
    if (gammaB !== undefined && Hi > gammaB * depths[i]) {
      Hi = gammaB * depths[i];
      broken[i] = 1;
    }
    H[i] = Hi;
  }
  return { Ks, Kr, H, thetaDeg, L, broken };
}

/**
 * Miche (1944) limiting wave height in finite depth: H_max = 0.142 L tanh(2πh/L).
 * Reference: Miche (1944), Annales des Ponts et Chaussées 114.
 * Validity: regular waves on a horizontal bed; tends to 0.142 L (deep) and 0.89 h (shallow).
 * @param {object} o
 * @param {number} o.T period [s]
 * @param {number} o.h depth [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Hmax:number, L:number, steepnessMax:number}} [m], [m], [-]
 */
export function breakingMiche({ T, h, g = G }) {
  const w = linearWave({ T, h, g });
  const s = 0.142 * Math.tanh(w.k * h);
  return { Hmax: s * w.L, L: w.L, steepnessMax: s };
}

/**
 * McCowan (1894) solitary-wave breaking limit Hb = 0.78 hb.
 * Reference: McCowan (1894), Phil. Mag. 38.
 * Validity: shallow water, flat bed; the index is exposed as an argument.
 * @param {object} o
 * @param {number} o.h depth [m]
 * @param {number} [o.index] breaker index Hb/hb [-], default 0.78
 * @returns {number} Hb [m]
 */
export function breakingMcCowan({ h, index = 0.78 }) {
  pos(h, 'h');
  return index * h;
}

/**
 * Weggel (1972) breaker index on a plane slope: Hb/hb = b − a Hb/(g T²),
 * a = 43.75 (1 − e^(−19 m)), b = 1.56/(1 + e^(−19.5 m)), solved explicitly for Hb at depth hb.
 * Reference: Weggel (1972), J. Waterways Harbors Coastal Eng. Div. ASCE 98(WW4); SPM (1984) Eq. 2-92..2-94.
 * Validity: regular waves, beach slope m <= about 0.1.
 * @param {object} o
 * @param {number} o.h depth at breaking [m]
 * @param {number} o.T period [s]
 * @param {number} o.slope beach slope tan β [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Hb:number, index:number, a:number, b:number}} Hb [m], Hb/hb [-], a [-], b [-]
 */
export function breakingWeggel({ h, T, slope, g = G }) {
  pos(h, 'h');
  pos(T, 'T');
  need(Number.isFinite(slope) && slope >= 0, 'slope must be >= 0');
  const a = 43.75 * (1 - Math.exp(-19 * slope));
  const b = 1.56 / (1 + Math.exp(-19.5 * slope));
  const Hb = (b * h) / (1 + (a * h) / (g * T * T));
  return { Hb, index: Hb / h, a, b };
}

/**
 * Goda (1975) breaker index: Hb/L0 = A {1 − exp[−1.5 (π hb/L0)(1 + K tan^(4/3) β)]}, A = 0.17, K = 15.
 * (Goda 2010 revised K to 11; pass slopeFactor = 11 for the revised form.)
 * Reference: Goda (1975) Rep. Port Harbour Res. Inst. 14(3); Goda (2000) "Random Seas and Design of Maritime
 * Structures" Eq. 3.21.
 * Validity: regular-wave breaking on uniform slopes 1/100 to 1/10; A = 0.12–0.18 describes the range for
 * individual random waves.
 * @param {object} o
 * @param {number} o.h depth at breaking [m]
 * @param {number} o.T period [s]
 * @param {number} o.slope beach slope tan β [-]
 * @param {number} [o.A] coefficient [-], default 0.17
 * @param {number} [o.slopeFactor] K [-], default 15
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Hb:number, index:number, L0:number}} Hb [m], Hb/hb [-], L0 [m]
 */
export function breakingGoda({ h, T, slope, A = 0.17, slopeFactor = 15, g = G }) {
  pos(h, 'h');
  pos(T, 'T');
  need(Number.isFinite(slope) && slope >= 0, 'slope must be >= 0');
  const L0 = (g * T * T) / TWO_PI;
  const Hb = A * L0 * (1 - Math.exp(-1.5 * ((Math.PI * h) / L0) * (1 + slopeFactor * Math.pow(slope, 4 / 3))));
  return { Hb, index: Hb / h, L0 };
}

/**
 * Goda's approximate formulas for significant and maximum wave height of random waves in and near the surf
 * zone on a uniform slope.
 *   h/L0 >= 0.2: H1/3 = Ks H0';  Hmax = 1.8 Ks H0'
 *   h/L0 <  0.2: H1/3 = min{β0 H0' + β1 h, βmax H0', Ks H0'};  Hmax = min{β0* H0' + β1* h, βmax* H0', 1.8 Ks H0'}
 *   β0 = 0.028 s^−0.38 exp(20 tan^1.5 θ), β1 = 0.52 exp(4.2 tan θ), βmax = max{0.92, 0.32 s^−0.29 exp(2.4 tan θ)}
 *   β0* = 0.052 s^−0.38 exp(20 tan^1.5 θ), β1* = 0.63 exp(3.8 tan θ), βmax* = max{1.65, 0.53 s^−0.29 exp(2.4 tan θ)}
 * with s = H0'/L0.
 * Reference: Goda (1975); Goda (2000) Eq. 3.25–3.26.
 * Validity: slopes 1/100 to 1/10, H0'/L0 about 0.005–0.08. Goda uses a nonlinear shoaling coefficient; the
 * linear Ks is used here unless `Ks` is supplied.
 * @param {object} o
 * @param {number} o.H0 equivalent deep-water significant wave height H0' [m]
 * @param {number} o.T significant wave period [s]
 * @param {number} o.h still-water depth [m]
 * @param {number} o.slope seabed slope tan θ [-]
 * @param {number} [o.Ks] shoaling coefficient override [-]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{Hs:number,Hmax:number,Ks:number,beta0:number,beta1:number,betaMax:number,surfZone:boolean}}
 *   heights [m], coefficients [-], surfZone true when a breaking branch governs H1/3
 */
export function godaSurfZoneHeights({ H0, T, h, slope, Ks, g = G }) {
  pos(H0, 'H0');
  pos(T, 'T');
  pos(h, 'h');
  need(Number.isFinite(slope) && slope >= 0, 'slope must be >= 0');
  const L0 = (g * T * T) / TWO_PI;
  const ks = Ks !== undefined ? Ks : shoalingCoefficient({ T, h, g });
  const s = H0 / L0;
  const e20 = Math.exp(20 * Math.pow(slope, 1.5));
  const beta0 = 0.028 * Math.pow(s, -0.38) * e20;
  const beta1 = 0.52 * Math.exp(4.2 * slope);
  const betaMax = Math.max(0.92, 0.32 * Math.pow(s, -0.29) * Math.exp(2.4 * slope));
  const b0s = 0.052 * Math.pow(s, -0.38) * e20;
  const b1s = 0.63 * Math.exp(3.8 * slope);
  const bMs = Math.max(1.65, 0.53 * Math.pow(s, -0.29) * Math.exp(2.4 * slope));
  let Hs = ks * H0;
  let Hmax = 1.8 * ks * H0;
  let surfZone = false;
  if (h / L0 < 0.2) {
    const br = Math.min(beta0 * H0 + beta1 * h, betaMax * H0);
    if (br < Hs) {
      Hs = br;
      surfZone = true;
    }
    Hmax = Math.min(b0s * H0 + b1s * h, bMs * H0, Hmax);
  }
  return { Hs, Hmax, Ks: ks, beta0, beta1, betaMax, surfZone };
}

/**
 * Battjes & Stive (1985) calibration of the Battjes–Janssen breaker parameter: γ = 0.5 + 0.4 tanh(33 s0),
 * s0 = Hrms0/L0 (deep-water steepness based on Hrms and peak period).
 * Reference: Battjes & Stive (1985), J. Geophys. Res. 90(C5).
 * @param {number} s0 deep-water steepness [-]
 * @returns {number} γ [-]
 */
export function battjesStiveGamma(s0) {
  need(Number.isFinite(s0) && s0 >= 0, 's0 must be >= 0');
  return 0.5 + 0.4 * Math.tanh(33 * s0);
}

/**
 * Fraction of broken waves Qb from the Battjes–Janssen (1978) implicit relation (1 − Qb)/(−ln Qb) = (Hrms/Hm)².
 * Reference: Battjes & Janssen (1978), Proc. 16th ICCE, Eq. 8.
 * @param {number} ratio Hrms/Hm [-] (>= 0)
 * @returns {number} Qb in [0, 1]
 */
export function battjesJanssenQb(ratio) {
  need(Number.isFinite(ratio) && ratio >= 0, 'ratio must be >= 0');
  if (ratio >= 1) return 1;
  const b2 = ratio * ratio;
  if (b2 < 0.0025) return Math.exp(-1 / b2); // below 1e-170: the root is its own asymptote
  const lo = Math.exp(-1 / b2);
  const eps = Math.min(0.5, 1 - b2);
  const f = (q) => 1 - q + b2 * Math.log(q);
  return brent(f, lo, 1 - eps, { tol: 1e-14 });
}

/**
 * One-dimensional cross-shore energy-flux balance for random waves with Battjes–Janssen (1978) bore
 * dissipation, bottom friction and wave set-up from the cross-shore momentum balance:
 *   d(E Cg cos θ)/dx = −(Db + Df),  E = ρ g Hrms²/8
 *   Db = (α/4) ρ g fp Qb Hm²,  Hm = (0.88/k) tanh(γ k h/0.88)
 *   Df = ρ fw/(16 √π) (2π fp Hrms / sinh kh)³  (time mean of τ·u with τ = ½ ρ fw u|u| and Rayleigh heights;
 *        identical to Thornton & Guza 1983 Eq. 40 with cf = fw)
 *   dSxx/dx + ρ g (d + η) dη/dx = 0,  Sxx = E [n (1 + cos²θ) − ½]   (Longuet-Higgins & Stewart 1964)
 * Direction follows Snell's law with the local depth including set-up. Marched shoreward with a Heun
 * (predictor–corrector) step on the supplied grid.
 * Reference: Battjes & Janssen (1978) Proc. 16th ICCE; Longuet-Higgins & Stewart (1964) Deep-Sea Res. 11;
 * Thornton & Guza (1983) J. Geophys. Res. 88(C10).
 * Validity: narrow-banded sea, straight parallel contours, mild slopes, no currents, no reflection.
 * @param {object} o
 * @param {ArrayLike<number>} o.x cross-shore coordinates, increasing shoreward [m]
 * @param {ArrayLike<number>} o.depth still-water depth at each x (positive down; <= 0 is dry land) [m]
 * @param {number} o.Hrms0 root-mean-square wave height at x[0] [m]
 * @param {number} o.Tp peak period [s]
 * @param {number} [o.theta0Deg] wave angle from shore-normal at x[0] [deg]
 * @param {number} [o.gamma] breaker parameter γ [-], default 0.8 (Battjes & Janssen 1978)
 * @param {number} [o.alpha] dissipation coefficient α [-], default 1
 * @param {number} [o.fw] wave friction factor [-], default 0 (no bottom friction)
 * @param {number} [o.rho] water density [kg/m³]
 * @param {number} [o.eta0] mean water level at x[0] [m]; default is the non-breaking set-down
 *   −k Hrms²/(8 sinh 2kh)
 * @param {boolean} [o.setup] include the set-up equation (default true)
 * @param {number} [o.hMin] total depth below which the march stops [m], default 0.01
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{x:Float64Array,Hrms:Float64Array,eta:Float64Array,h:Float64Array,thetaDeg:Float64Array,
 *   Qb:Float64Array,Db:Float64Array,Df:Float64Array,energyFlux:Float64Array,Sxx:Float64Array,
 *   setupMax:number,nWet:number}}
 *   Hrms [m], eta mean water level [m], h total depth [m], theta [deg], Qb [-], Db and Df [W/m²],
 *   energyFlux = E Cg cos θ [W/m], Sxx [N/m], setupMax [m], nWet = number of nodes computed
 */
export function energyFluxBalance({
  x, depth, Hrms0, Tp, theta0Deg = 0, gamma = 0.8, alpha = 1, fw = 0, rho = 1025, eta0, setup = true,
  hMin = 0.01, g = G,
}) {
  need(x && depth && x.length === depth.length && x.length >= 2, 'x and depth must be arrays of equal length >= 2');
  pos(Hrms0, 'Hrms0');
  pos(Tp, 'Tp');
  pos(gamma, 'gamma');
  need(alpha >= 0, 'alpha must be >= 0');
  need(fw >= 0, 'fw must be >= 0');
  need(Math.abs(theta0Deg) < 90, 'theta0Deg must lie in (-90, 90)');
  need(depth[0] > hMin, 'depth[0] must be wet');
  const n = x.length;
  for (let i = 1; i < n; i++) need(x[i] > x[i - 1], 'x must be strictly increasing');
  const out = {
    x: Float64Array.from(x), Hrms: new Float64Array(n), eta: new Float64Array(n), h: new Float64Array(n),
    thetaDeg: new Float64Array(n), Qb: new Float64Array(n), Db: new Float64Array(n), Df: new Float64Array(n),
    energyFlux: new Float64Array(n), Sxx: new Float64Array(n), setupMax: 0, nWet: 0,
  };
  const fp = 1 / Tp;
  const sqrtPi = Math.sqrt(Math.PI);
  const ws = linearWave({ T: Tp, h: depth[0], g });

  // Local state for a given energy flux F and total depth h.
  const st = { H: 0, theta: 0, Qb: 0, Db: 0, Df: 0, Sxx: 0, D: 0 };
  const evalState = (F, h) => {
    const w = linearWave({ T: Tp, h, g });
    const sinT = Math.max(-1, Math.min(1, snell * w.C));
    const cosT = Math.sqrt(1 - sinT * sinT);
    const E = Math.max(F, 0) / (w.Cg * cosT);
    let H = Math.sqrt((8 * E) / (rho * g));
    const kh = w.k * h;
    const Hm = (0.88 / w.k) * Math.tanh((gamma * kh) / 0.88);
    // Rayleigh distribution truncated at Hm cannot carry Hrms > Hm.
    if (alpha > 0 && H > Hm) H = Hm;
    const Ecl = (rho * g * H * H) / 8;
    const Qb = alpha > 0 ? battjesJanssenQb(H / Hm) : 0;
    st.H = H;
    st.theta = Math.asin(sinT);
    st.Qb = Qb;
    st.Db = 0.25 * alpha * rho * g * fp * Qb * Hm * Hm;
    st.Df = fw > 0 && kh < 30 ? ((rho * fw) / (16 * sqrtPi)) * Math.pow((TWO_PI * fp * H) / Math.sinh(kh), 3) : 0;
    st.D = st.Db + st.Df;
    st.Sxx = Ecl * (w.n * (1 + cosT * cosT) - 0.5);
    st.F = Ecl * w.Cg * cosT;
    return st;
  };
  const store = (i, h, eta) => {
    out.Hrms[i] = st.H;
    out.eta[i] = eta;
    out.h[i] = h;
    out.thetaDeg[i] = (st.theta * 180) / Math.PI;
    out.Qb[i] = st.Qb;
    out.Db[i] = st.Db;
    out.Df[i] = st.Df;
    out.energyFlux[i] = st.F;
    out.Sxx[i] = st.Sxx;
  };

  let eta = eta0 !== undefined ? eta0 : -(ws.k * Hrms0 * Hrms0) / (8 * Math.sinh(2 * ws.k * depth[0]));
  if (!setup) eta = eta0 !== undefined ? eta0 : 0;
  let h = depth[0] + eta;
  need(h > hMin, 'depth[0] + eta0 must be wet');
  const w0 = linearWave({ T: Tp, h, g });
  const snell = Math.sin((theta0Deg * Math.PI) / 180) / w0.C;
  const cos0 = Math.cos((theta0Deg * Math.PI) / 180);
  let F = ((rho * g * Hrms0 * Hrms0) / 8) * w0.Cg * cos0;
  evalState(F, h);
  F = st.F;
  store(0, h, eta);
  let nWet = 1;
  let setupMax = eta;
  for (let i = 1; i < n; i++) {
    const dx = x[i] - x[i - 1];
    const Dp = st.D;
    const Sp = st.Sxx;
    const hp = h;
    let etaN = eta;
    let hN = depth[i] + etaN;
    if (hN <= hMin) break;
    let FN = Math.max(F - Dp * dx, 0);
    let ok = true;
    for (let it = 0; it < 8; it++) {
      evalState(FN, hN);
      FN = Math.max(F - 0.5 * (Dp + st.D) * dx, 0);
      if (setup) {
        etaN = eta - (st.Sxx - Sp) / (rho * g * 0.5 * (hp + hN));
        hN = depth[i] + etaN;
        if (hN <= hMin) {
          ok = false;
          break;
        }
      }
    }
    if (!ok) break;
    evalState(FN, hN);
    F = st.F;
    eta = etaN;
    h = hN;
    store(i, h, eta);
    if (eta > setupMax) setupMax = eta;
    nWet++;
  }
  out.setupMax = setupMax;
  out.nWet = nWet;
  return out;
}

/**
 * Fresnel integrals C(x) = ∫0^x cos(π t²/2) dt and S(x) = ∫0^x sin(π t²/2) dt.
 * Power series for |x| <= 2, continued fraction of the complex error function beyond.
 * Reference: Abramowitz & Stegun (1964) §7.3 (series 7.3.11, 7.3.13; relation to erf 7.3.22);
 * Press et al. (2007) "Numerical Recipes" §6.8.
 * Accuracy: about 1e-12 absolute.
 * @param {number} x argument [-]
 * @returns {{C:number, S:number}}
 */
export function fresnel(x) {
  need(Number.isFinite(x) || Math.abs(x) === Infinity, 'x must be a number');
  const ax = Math.abs(x);
  let C;
  let S;
  if (ax === Infinity) {
    C = 0.5;
    S = 0.5;
  } else if (ax <= 2) {
    const t = 0.5 * Math.PI * ax * ax;
    const t2 = t * t;
    let termC = ax; // (−1)^n t^(2n) x / (2n)!
    let termS = (ax * t); // (−1)^n t^(2n+1) x / (2n+1)!
    C = termC;
    S = termS / 3;
    for (let n = 1; n < 80; n++) {
      termC *= -t2 / ((2 * n - 1) * (2 * n));
      termS *= -t2 / ((2 * n) * (2 * n + 1));
      const dc = termC / (4 * n + 1);
      const ds = termS / (4 * n + 3);
      C += dc;
      S += ds;
      if (Math.abs(dc) < 1e-17 && Math.abs(ds) < 1e-17) break;
    }
  } else {
    // Modified Lentz evaluation of the continued fraction for the complex error function
    // erfc(z) e^(z²), z = √π (1 − i) x / 2, from which C + iS = (1+i)/2 · [1 − e^(iπx²/2) · (x − ix) · h].
    const pix2 = Math.PI * ax * ax;
    let br = 1;
    const bi = -pix2;
    let den = br * br + bi * bi;
    let dr = br / den;
    let di = -bi / den;
    let hr = dr;
    let hi = di;
    let cr = 1e300;
    let ci = 0;
    let nn = -1;
    for (let kk = 2; kk < 400; kk++) {
      nn += 2;
      const an = -nn * (nn + 1);
      br += 4;
      // d = 1/(an d + b)
      const tr = an * dr + br;
      const ti = an * di + bi;
      den = tr * tr + ti * ti;
      dr = tr / den;
      di = -ti / den;
      // c = b + an / c
      den = cr * cr + ci * ci;
      const ncr = br + (an * cr) / den;
      const nci = bi - (an * ci) / den;
      cr = ncr;
      ci = nci;
      const delr = cr * dr - ci * di;
      const deli = cr * di + ci * dr;
      const nhr = hr * delr - hi * deli;
      hi = hr * deli + hi * delr;
      hr = nhr;
      if (Math.abs(delr - 1) + Math.abs(deli) < 1e-16) break;
    }
    // h <- (ax − i ax) h
    const gr = ax * (hr + hi);
    const gi = ax * (hi - hr);
    const cs = Math.cos(0.5 * pix2);
    const sn = Math.sin(0.5 * pix2);
    const er = 1 - (cs * gr - sn * gi);
    const ei = -(cs * gi + sn * gr);
    C = 0.5 * (er - ei);
    S = 0.5 * (er + ei);
  }
  return x < 0 ? { C: -C, S: -S } : { C, S };
}

// f(σ) = (1+i)/2 ∫_{−∞}^{σ} exp(−iπt²/2) dt, written to a two-element output.
function sommerfeldF(sigma, out) {
  const fr = fresnel(sigma);
  const a = fr.C + 0.5;
  const b = -(fr.S + 0.5);
  out[0] = 0.5 * (a - b);
  out[1] = 0.5 * (a + b);
}

const tmpA = [0, 0];
const tmpB = [0, 0];

function semiInfiniteField(x, y, k, theta0) {
  const r = Math.hypot(x, y);
  let theta = Math.atan2(y, x);
  if (theta < 0) theta += TWO_PI;
  const root = 2 * Math.sqrt((k * r) / Math.PI);
  sommerfeldF(root * Math.sin(0.5 * (theta - theta0)), tmpA);
  sommerfeldF(-root * Math.sin(0.5 * (theta + theta0)), tmpB);
  const p1 = -k * r * Math.cos(theta - theta0);
  const p2 = -k * r * Math.cos(theta + theta0);
  const c1 = Math.cos(p1);
  const s1 = Math.sin(p1);
  const c2 = Math.cos(p2);
  const s2 = Math.sin(p2);
  return [tmpA[0] * c1 - tmpA[1] * s1 + tmpB[0] * c2 - tmpB[1] * s2, tmpA[0] * s1 + tmpA[1] * c1 + tmpB[0] * s2 + tmpB[1] * c2];
}

/**
 * Diffraction coefficient behind (and around) a thin, rigid, fully reflecting semi-infinite breakwater in
 * constant depth: Sommerfeld's solution as applied to water waves by Penney & Price (1952),
 *   F(r,θ) = f(σ) e^(−ikr cos(θ−θ0)) + f(σ') e^(−ikr cos(θ+θ0)),
 *   σ = 2 √(kr/π) sin[(θ−θ0)/2],  σ' = −2 √(kr/π) sin[(θ+θ0)/2],
 *   f(σ) = (1+i)/2 ∫_{−∞}^{σ} e^(−iπt²/2) dt,  Kd = |F|.
 * Geometry: the breakwater lies on y = 0, x >= 0 with its tip at the origin; θ is measured from the
 * breakwater (+x axis) anticlockwise; the incident wave travels in direction θ0 (so it arrives from the
 * y < 0 side and the geometric shadow is 0 < θ < θ0).
 * Reference: Penney & Price (1952), Phil. Trans. R. Soc. A 244; Sommerfeld (1896).
 * Validity: linear waves, uniform depth, 0 < θ0 < 180°.
 * @param {object} o
 * @param {number} o.x coordinate along the breakwater from the tip [m]
 * @param {number} o.y coordinate normal to the breakwater, positive on the lee side [m]
 * @param {number} o.L local wavelength [m]
 * @param {number} [o.theta0Deg] wave travel direction measured from the breakwater [deg], default 90
 * @returns {{Kd:number, re:number, im:number, phase:number}} Kd [-], complex amplitude relative to a unit
 *   incident wave with zero phase at the tip, phase [rad]
 */
export function diffractionSemiInfinite({ x, y, L, theta0Deg = 90 }) {
  pos(L, 'L');
  need(Number.isFinite(x) && Number.isFinite(y), 'x and y must be finite');
  need(theta0Deg > 0 && theta0Deg < 180, 'theta0Deg must lie in (0, 180)');
  const F = semiInfiniteField(x, y, TWO_PI / L, (theta0Deg * Math.PI) / 180);
  return { Kd: Math.hypot(F[0], F[1]), re: F[0], im: F[1], phase: Math.atan2(F[1], F[0]) };
}

/**
 * Diffraction through a breakwater gap by superposition of the two semi-infinite solutions,
 * F = F_right + F_left − F_incident (Penney & Price 1952).
 * Geometry: breakwater arms on y = 0 for |x| >= B/2, lee side y > 0, wave travel direction θ0 from the +x axis.
 * Reference: Penney & Price (1952) §6; SPM (1984) §2.IV.
 * Validity: gap width B larger than about one wavelength (the two tips are treated as independent);
 * evaluated on the lee side y >= 0 only.
 * @param {object} o
 * @param {number} o.x alongshore coordinate from the gap centre [m]
 * @param {number} o.y distance behind the gap (>= 0) [m]
 * @param {number} o.L wavelength [m]
 * @param {number} o.B gap width [m]
 * @param {number} [o.theta0Deg] wave travel direction from the +x axis [deg], default 90 (normal incidence)
 * @returns {{Kd:number, re:number, im:number}} Kd [-] and complex amplitude (zero incident phase at gap centre)
 */
export function diffractionGap({ x, y, L, B, theta0Deg = 90 }) {
  pos(L, 'L');
  pos(B, 'B');
  need(Number.isFinite(x) && Number.isFinite(y) && y >= 0, 'y must be >= 0 (lee side)');
  need(theta0Deg > 0 && theta0Deg < 180, 'theta0Deg must lie in (0, 180)');
  const k = TWO_PI / L;
  const t0 = (theta0Deg * Math.PI) / 180;
  const ph = k * 0.5 * B * Math.cos(t0);
  const cp = Math.cos(ph);
  const sp = Math.sin(ph);
  const a = semiInfiniteField(x - 0.5 * B, y, k, t0);
  const ar = a[0] * cp + a[1] * sp; // × e^(−i ph)
  const ai = a[1] * cp - a[0] * sp;
  const b = semiInfiniteField(-(x + 0.5 * B), y, k, Math.PI - t0);
  const br = b[0] * cp - b[1] * sp; // × e^(+i ph)
  const bi = b[1] * cp + b[0] * sp;
  const pi = -k * (x * Math.cos(t0) + y * Math.sin(t0));
  const re = ar + br - Math.cos(pi);
  const im = ai + bi - Math.sin(pi);
  return { Kd: Math.hypot(re, im), re, im };
}

/**
 * Wave reflection coefficient of a sloping structure, Seelig & Ahrens (1981): Kr = a ξ² / (b + ξ²).
 * Published coefficient pairs: smooth impermeable slopes a = 1.0, b = 5.5; conservative estimate for
 * rough permeable (rubble-mound) slopes a = 0.6, b = 6.6.
 * Reference: Seelig & Ahrens (1981), CERC Technical Paper 81-1; Seelig (1983).
 * Validity: non-overtopped slopes; ξ based on the incident wave at the toe.
 * @param {object} o
 * @param {number} o.xi surf-similarity (Iribarren) number [-]
 * @param {number} [o.a] coefficient [-], default 1.0
 * @param {number} [o.b] coefficient [-], default 5.5
 * @returns {number} Kr [-]
 */
export function reflectionSeeligAhrens({ xi, a = 1.0, b = 5.5 }) {
  need(Number.isFinite(xi) && xi >= 0, 'xi must be >= 0');
  pos(a, 'a');
  pos(b, 'b');
  return (a * xi * xi) / (b + xi * xi);
}

/**
 * Linear waves meeting a uniform collinear current: Doppler-shifted dispersion (ω − kU)² = g k tanh(kh) and
 * wave-action conservation (E/σ)(U + Cgr) = const, referenced to the same depth without current.
 * Blocking occurs when the opposing current reaches the relative group velocity (U = −Cgr).
 * Reference: Longuet-Higgins & Stewart (1961) J. Fluid Mech. 10; Bretherton & Garrett (1968);
 * Peregrine (1976) Adv. Appl. Mech. 16.
 * Validity: steady, depth-uniform current, linear waves, no breaking or reflection; the amplification is
 * singular at the blocking point.
 * @param {object} o
 * @param {number} o.T absolute (fixed-frame) period [s]
 * @param {number} o.h depth [m]
 * @param {number} o.U current speed, positive following the waves, negative opposing [m/s]
 * @param {number} [o.H0] wave height at the same depth without current [m]
 * @param {number} [o.g] gravity [m/s²]
 * @returns {{blocked:boolean,k:number,L:number,sigma:number,Tr:number,Cr:number,Cgr:number,Ca:number,
 *   Cga:number,heightRatio:number,H:number,Ublock:number}}
 *   k [rad/m], L [m], sigma relative frequency [rad/s], Tr relative period [s], relative and absolute phase
 *   and group speeds [m/s], H/H0 [-], H [m] (NaN if H0 omitted), Ublock = current that just blocks this
 *   period [m/s] (negative). When blocked the wave fields are NaN.
 */
export function waveCurrent({ T, h, U, H0, g = G }) {
  pos(T, 'T');
  pos(h, 'h');
  need(Number.isFinite(U), 'U must be finite');
  const omega = TWO_PI / T;
  const w0 = linearWave({ T, h, g });
  const sig = (k) => Math.sqrt(g * k * Math.tanh(k * h));
  const cgr = (k) => {
    const kh = k * h;
    const nn = kh > 350 ? 0.5 : 0.5 * (1 + (2 * kh) / Math.sinh(2 * kh));
    return (nn * sig(k)) / k;
  };
  // Blocking current: maximise over k the opposing speed that still admits a solution, i.e. solve
  // σ(k) − k Cgr(k) = ω (tangency of the Doppler line with the dispersion curve), then Ublock = −Cgr.
  const tang = (k) => sig(k) - k * cgr(k) - omega;
  let kHi = (4 * omega * omega) / g; // deep-water tangency point; finite depth pushes it higher
  for (let i = 0; i < 200 && tang(kHi) <= 0; i++) kHi *= 2;
  const kb = brent(tang, w0.k, kHi, { tol: 1e-13 });
  const Ublock = -cgr(kb);
  const nan = Number.NaN;
  if (U < Ublock) {
    return { blocked: true, k: nan, L: nan, sigma: nan, Tr: nan, Cr: nan, Cgr: nan, Ca: nan, Cga: nan, heightRatio: nan, H: nan, Ublock };
  }
  let k;
  if (U === 0) k = w0.k;
  else if (U > 0) k = brent((kk) => sig(kk) + kk * U - omega, 1e-12 * w0.k, w0.k, { tol: 1e-14 });
  else k = U === Ublock ? kb : brent((kk) => sig(kk) + kk * U - omega, w0.k, kb, { tol: 1e-14 });
  const sigma = sig(k);
  const Cgr = cgr(k);
  const heightRatio = Math.sqrt(((w0.Cg / omega) * sigma) / (U + Cgr));
  return {
    blocked: false, k, L: TWO_PI / k, sigma, Tr: TWO_PI / sigma, Cr: sigma / k, Cgr, Ca: omega / k, Cga: U + Cgr,
    heightRatio, H: H0 !== undefined ? H0 * heightRatio : nan, Ublock,
  };
}
