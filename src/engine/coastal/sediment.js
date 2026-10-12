// Sediment properties, threshold of motion, bed shear stress, transport formulae and rock sizing under
// currents. SI units; transport rates are volumetric (solid volume, no pores) per unit width [m²/s].

import { brent } from '../math/roots.js';

const G = 9.80665;
const KAPPA = 0.4;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function pos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Dimensionless grain size D* = d [g (s − 1)/ν²]^(1/3), s = ρs/ρ.
 * Reference: van Rijn (1984a) J. Hydraul. Eng. 110(10); Soulsby (1997) "Dynamics of Marine Sands" Eq. 75.
 * @param {object} o
 * @param {number} o.d grain diameter [m]
 * @param {number} [o.rhoS] sediment density [kg/m³], default 2650
 * @param {number} [o.rho] water density [kg/m³], default 1025
 * @param {number} [o.nu] kinematic viscosity [m²/s], default 1.0e-6 (use 1.36e-6 for sea water at 10 °C)
 * @param {number} [o.g] gravity [m/s²]
 * @returns {number} D* [-]
 */
export function dimensionlessGrainSize({ d, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  pos(d, 'd');
  pos(nu, 'nu');
  need(rhoS > rho && rho > 0, 'rhoS must exceed rho');
  return d * Math.cbrt((g * (rhoS / rho - 1)) / (nu * nu));
}

/**
 * Settling velocity of natural sand, Soulsby (1997): ws = (ν/d) [√(10.36² + 1.049 D*³) − 10.36].
 * Reference: Soulsby (1997) Eq. 102. Validity: natural (irregular) grains, all D*, low concentration.
 * @param {object} o as dimensionlessGrainSize
 * @returns {number} ws [m/s]
 */
export function settlingVelocitySoulsby({ d, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  const Ds = dimensionlessGrainSize({ d, rhoS, rho, nu, g });
  return (nu / d) * (Math.sqrt(10.36 * 10.36 + 1.049 * Ds ** 3) - 10.36);
}

/**
 * Settling velocity, van Rijn (1984b):
 *   d < 100 µm:        ws = (s − 1) g d² / (18 ν)                       (Stokes)
 *   100 µm–1000 µm:    ws = (10 ν/d) [√(1 + 0.01 (s − 1) g d³/ν²) − 1]
 *   d > 1000 µm:       ws = 1.1 √((s − 1) g d)
 * Reference: van Rijn (1984b) J. Hydraul. Eng. 110(11). Validity: natural sand, clear water.
 * @param {object} o as dimensionlessGrainSize
 * @returns {number} ws [m/s]
 */
export function settlingVelocityVanRijn({ d, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  pos(d, 'd');
  pos(nu, 'nu');
  need(rhoS > rho && rho > 0, 'rhoS must exceed rho');
  const s1 = rhoS / rho - 1;
  if (d < 1e-4) return (s1 * g * d * d) / (18 * nu);
  if (d <= 1e-3) return ((10 * nu) / d) * (Math.sqrt(1 + (0.01 * s1 * g * d ** 3) / (nu * nu)) - 1);
  return 1.1 * Math.sqrt(s1 * g * d);
}

/**
 * Threshold Shields parameter, Soulsby & Whitehouse (1997):
 *   θcr = 0.30/(1 + 1.2 D*) + 0.055 [1 − exp(−0.020 D*)].
 * Reference: Soulsby & Whitehouse (1997) Proc. Pacific Coasts and Ports '97; Soulsby (1997) Eq. 77.
 * Validity: non-cohesive grains on a flat horizontal bed, currents and/or waves.
 * @param {number} Dstar dimensionless grain size [-]
 * @returns {number} θcr [-]
 */
export function criticalShields(Dstar) {
  pos(Dstar, 'Dstar');
  return 0.3 / (1 + 1.2 * Dstar) + 0.055 * (1 - Math.exp(-0.02 * Dstar));
}

/**
 * Threshold bed shear stress τcr = θcr g (ρs − ρ) d and the matching shear velocity.
 * Reference: Soulsby (1997) §6.4.
 * @param {object} o as dimensionlessGrainSize
 * @returns {{Dstar:number, thetaCr:number, tauCr:number, uStarCr:number}} [-], [-], [Pa], [m/s]
 */
export function criticalShearStress({ d, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  const Dstar = dimensionlessGrainSize({ d, rhoS, rho, nu, g });
  const thetaCr = criticalShields(Dstar);
  const tauCr = thetaCr * g * (rhoS - rho) * d;
  return { Dstar, thetaCr, tauCr, uStarCr: Math.sqrt(tauCr / rho) };
}

/**
 * Shields parameter θ = τ / (g (ρs − ρ) d).
 * Reference: Shields (1936); Soulsby (1997) §6.4.
 * @param {object} o
 * @param {number} o.tau bed shear stress (skin friction) [Pa]
 * @param {number} o.d grain diameter [m]
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.g] [m/s²]
 * @returns {number} θ [-]
 */
export function shieldsParameter({ tau, d, rhoS = 2650, rho = 1025, g = G }) {
  need(Number.isFinite(tau) && tau >= 0, 'tau must be >= 0');
  pos(d, 'd');
  need(rhoS > rho, 'rhoS must exceed rho');
  return tau / (g * (rhoS - rho) * d);
}

/**
 * Mobility check: compares the Shields parameter with the Soulsby–Whitehouse threshold.
 * Reference: Soulsby (1997) §6.4.
 * @param {object} o
 * @param {number} o.tau skin-friction bed shear stress (maximum over the wave cycle for waves) [Pa]
 * @param {number} o.d grain diameter [m]
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.nu] [m²/s]
 * @param {number} [o.g] [m/s²]
 * @returns {{theta:number, thetaCr:number, ratio:number, mobile:boolean, Dstar:number}} ratio = θ/θcr [-]
 */
export function mobilityCheck({ tau, d, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  const c = criticalShearStress({ d, rhoS, rho, nu, g });
  const theta = shieldsParameter({ tau, d, rhoS, rho, g });
  return { theta, thetaCr: c.thetaCr, ratio: theta / c.thetaCr, mobile: theta > c.thetaCr, Dstar: c.Dstar };
}

/**
 * Grain (skin-friction) roughness: Nikuradse ks = 2.5 d50, z0 = ks/30 = d50/12.
 * Reference: Soulsby (1997) §3.3. Validity: hydraulically rough flat bed.
 * @param {number} d50 median grain diameter [m]
 * @returns {{ks:number, z0:number}} [m]
 */
export function grainRoughness(d50) {
  pos(d50, 'd50');
  return { ks: 2.5 * d50, z0: d50 / 12 };
}

/**
 * Form-drag roughness of ripples or dunes: z0f = ar Δ²/λ, and the equivalent ks = 30 z0f.
 * Reference: Soulsby (1997) §7.3 (ar = 1 typical; published range 0.3–3).
 * @param {object} o
 * @param {number} o.height bedform height Δ [m]
 * @param {number} o.length bedform wavelength λ [m]
 * @param {number} [o.ar] coefficient [-], default 1
 * @returns {{z0:number, ks:number}} [m]
 */
export function bedformRoughness({ height, length, ar = 1 }) {
  pos(height, 'height');
  pos(length, 'length');
  pos(ar, 'ar');
  const z0 = (ar * height * height) / length;
  return { z0, ks: 30 * z0 };
}

/**
 * Current-ripple dimensions on sand: λr = 1000 d50, Δr = λr/7.
 * Reference: Soulsby (1997) §7.2. Validity: d50 < 0.8 mm, flow above threshold and below wash-out.
 * @param {number} d50 [m]
 * @returns {{length:number, height:number, valid:boolean}} [m]
 */
export function currentRipples(d50) {
  pos(d50, 'd50');
  const length = 1000 * d50;
  return { length, height: length / 7, valid: d50 < 0.8e-3 };
}

/**
 * Bed shear stress under a steady current from the depth-averaged velocity.
 *   'log':       CD = [κ/(ln(h/z0) − 1)]², κ = 0.40, τ = ρ CD U²     (Soulsby 1997 §3.3)
 *   'manning':   τ = ρ g n² U² / h^(1/3)                               (Manning–Strickler)
 *   'colebrook': C = 18 log10(12 h/ks), τ = ρ g U²/C²                 (Colebrook–White, rough turbulent)
 * Reference: Soulsby (1997) §3.3; van Rijn (1984a).
 * Validity: steady uniform rough turbulent flow, h >> z0 (h/z0 > e is enforced).
 * @param {object} o
 * @param {number} o.U depth-averaged current speed [m/s]
 * @param {number} o.h depth [m]
 * @param {'log'|'manning'|'colebrook'} [o.method] default 'log'
 * @param {number} [o.z0] bed roughness length [m] ('log')
 * @param {number} [o.n] Manning coefficient [s/m^(1/3)] ('manning')
 * @param {number} [o.ks] Nikuradse roughness [m] ('colebrook')
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.g] [m/s²]
 * @returns {{tau:number, uStar:number, Cd:number, chezy:number}} [Pa], [m/s], drag coefficient [-],
 *   Chézy coefficient [m^(1/2)/s]
 */
export function bedShearCurrent({ U, h, method = 'log', z0, n, ks, rho = 1025, g = G }) {
  need(Number.isFinite(U), 'U must be finite');
  pos(h, 'h');
  let Cd;
  if (method === 'log') {
    pos(z0, 'z0');
    need(h / z0 > Math.E * 1.0001, 'h/z0 must exceed e for the depth-averaged log law');
    Cd = (KAPPA / (Math.log(h / z0) - 1)) ** 2;
  } else if (method === 'manning') {
    pos(n, 'n');
    Cd = (g * n * n) / Math.cbrt(h);
  } else if (method === 'colebrook') {
    pos(ks, 'ks');
    need((12 * h) / ks > 1.0001, '12 h/ks must exceed 1');
    const C = 18 * Math.log10((12 * h) / ks);
    Cd = g / (C * C);
  } else throw new RangeError("method must be 'log', 'manning' or 'colebrook'");
  const tau = rho * Cd * U * U;
  return { tau, uStar: Math.sqrt(tau / rho), Cd, chezy: Math.sqrt(g / Cd) };
}

/**
 * Near-bed wave orbital velocity amplitude and excursion from linear theory:
 * Uw = π H /(T sinh kh), A = Uw T/(2π).
 * Reference: Soulsby (1997) §4.4 (monochromatic waves). Validity: linear waves.
 * @param {object} o
 * @param {number} o.H wave height [m]
 * @param {number} o.T period [s]
 * @param {number} o.h depth [m]
 * @param {number} [o.g] [m/s²]
 * @returns {{Uw:number, A:number, k:number}} [m/s], [m], [rad/m]
 */
export function waveOrbital({ H, T, h, g = G }) {
  pos(H, 'H');
  pos(T, 'T');
  pos(h, 'h');
  const k0h = ((2 * Math.PI / T) ** 2 * h) / g;
  let kh = k0h / Math.sqrt(Math.tanh(k0h));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(kh);
    const dd = (kh * th - k0h) / (th + kh * (1 - th * th));
    kh -= dd;
    if (Math.abs(dd) < 1e-14 * kh) break;
  }
  const Uw = kh > 700 ? 0 : (Math.PI * H) / (T * Math.sinh(kh));
  return { Uw, A: (Uw * T) / (2 * Math.PI), k: kh / h };
}

/**
 * Rough-bed wave friction factor.
 *   'soulsby': fw = 1.39 (A/z0)^−0.52                                   (Soulsby 1997 Eq. 62a)
 *   'swart':   fw = exp[5.213 (A/ks)^−0.194 − 5.977], fw = 0.3 for A/ks <= 1.57   (Swart 1974)
 * Reference: Soulsby (1997) §4.5; Swart (1974) Delft Hydraulics Lab. Publ. 131.
 * Validity: rough turbulent oscillatory boundary layer.
 * @param {object} o
 * @param {number} o.A orbital excursion amplitude [m]
 * @param {'soulsby'|'swart'} [o.method] default 'soulsby'
 * @param {number} [o.z0] roughness length [m] ('soulsby'; taken as ks/30 if only ks is given)
 * @param {number} [o.ks] Nikuradse roughness [m] ('swart'; taken as 30 z0 if only z0 is given)
 * @returns {number} fw [-]
 */
export function waveFrictionFactor({ A, method = 'soulsby', z0, ks }) {
  pos(A, 'A');
  const zz = z0 !== undefined ? z0 : ks !== undefined ? ks / 30 : undefined;
  const kk = ks !== undefined ? ks : z0 !== undefined ? 30 * z0 : undefined;
  pos(zz, 'z0 or ks');
  if (method === 'soulsby') return 1.39 * Math.pow(A / zz, -0.52);
  if (method === 'swart') return A / kk <= 1.57 ? 0.3 : Math.exp(5.213 * Math.pow(A / kk, -0.194) - 5.977);
  throw new RangeError("method must be 'soulsby' or 'swart'");
}

/**
 * Amplitude of the bed shear stress under waves: τw = ½ ρ fw Uw².
 * Reference: Soulsby (1997) §4.5.
 * @param {object} o
 * @param {number} o.Uw orbital velocity amplitude [m/s]
 * @param {number} o.T wave period [s]
 * @param {'soulsby'|'swart'} [o.method]
 * @param {number} [o.z0] [m]
 * @param {number} [o.ks] [m]
 * @param {number} [o.rho] [kg/m³]
 * @returns {{tau:number, fw:number, A:number}} [Pa], [-], [m]
 */
export function bedShearWave({ Uw, T, method = 'soulsby', z0, ks, rho = 1025 }) {
  need(Number.isFinite(Uw) && Uw >= 0, 'Uw must be >= 0');
  pos(T, 'T');
  if (Uw === 0) return { tau: 0, fw: 0, A: 0 };
  const A = (Uw * T) / (2 * Math.PI);
  const fw = waveFrictionFactor({ A, method, z0, ks });
  return { tau: 0.5 * rho * fw * Uw * Uw, fw, A };
}

/**
 * Combined wave–current bed shear stress, Soulsby (1995) "DATA2" parameterisation:
 *   τm = τc [1 + 1.2 (τw/(τc + τw))^3.2]
 *   τmax = [(τm + τw cos φ)² + (τw sin φ)²]^(1/2)
 * with τc, τw the current-alone and wave-alone stresses and φ the angle between them.
 * Reference: Soulsby (1995) in "Advances in Coastal Morphodynamics" (MAST G8-M); Soulsby (1997) Eq. 69–70.
 * Validity: rough turbulent flow; fitted to 131 laboratory and field data.
 * @param {object} o
 * @param {number} o.tauC current-alone bed shear stress [Pa]
 * @param {number} o.tauW wave-alone bed shear-stress amplitude [Pa]
 * @param {number} [o.phiDeg] angle between current and wave directions [deg], default 0
 * @returns {{tauMean:number, tauMax:number}} [Pa]
 */
export function combinedShearSoulsby({ tauC, tauW, phiDeg = 0 }) {
  need(Number.isFinite(tauC) && tauC >= 0, 'tauC must be >= 0');
  need(Number.isFinite(tauW) && tauW >= 0, 'tauW must be >= 0');
  if (tauC + tauW === 0) return { tauMean: 0, tauMax: 0 };
  const tauMean = tauC * (1 + 1.2 * Math.pow(tauW / (tauC + tauW), 3.2));
  const phi = (phiDeg * Math.PI) / 180;
  return { tauMean, tauMax: Math.hypot(tauMean + tauW * Math.cos(phi), tauW * Math.sin(phi)) };
}

/**
 * Meyer-Peter & Müller (1948) bed load: Φ = 8 (θ − θcr)^1.5, θcr = 0.047, qb = Φ [g (s − 1) d³]^(1/2).
 * Reference: Meyer-Peter & Müller (1948) Proc. 2nd IAHR Congress; Soulsby (1997) §9.2.
 * Validity: coarse sand and gravel, θ (skin friction) below about 0.25; coefficient and threshold exposed
 * (Wong & Parker 2006 re-analysis: 3.97 and 0.0495).
 * @param {object} o
 * @param {number} o.theta skin-friction Shields parameter [-]
 * @param {number} o.d grain diameter [m]
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.thetaCr] default 0.047
 * @param {number} [o.coef] default 8
 * @param {number} [o.g] [m/s²]
 * @returns {{qb:number, phi:number}} bed-load transport [m²/s], dimensionless rate [-]
 */
export function meyerPeterMuller({ theta, d, rhoS = 2650, rho = 1025, thetaCr = 0.047, coef = 8, g = G }) {
  need(Number.isFinite(theta) && theta >= 0, 'theta must be >= 0');
  pos(d, 'd');
  need(rhoS > rho, 'rhoS must exceed rho');
  const phi = theta > thetaCr ? coef * Math.pow(theta - thetaCr, 1.5) : 0;
  return { qb: phi * Math.sqrt(g * (rhoS / rho - 1) * d ** 3), phi };
}

/**
 * Engelund & Hansen (1967) total load: f Φ = 0.1 θ^(5/2) with f = 2 u*²/U², i.e.
 *   qt = 0.05 U² θ^(5/2) [g (s − 1) d50³]^(1/2) / u*² = 0.05 U⁵ / [(s − 1)² g^(1/2) d50 C³]  (C = Chézy).
 * Reference: Engelund & Hansen (1967) "A monograph on sediment transport in alluvial streams";
 * Soulsby (1997) §10.3. Validity: rivers and tidal flows with dunes, d50 > 0.15 mm; θ and u* are
 * total-stress quantities; no threshold.
 * @param {object} o
 * @param {number} o.U depth-averaged velocity [m/s]
 * @param {number} o.tau total bed shear stress [Pa]
 * @param {number} o.d50 [m]
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.g] [m/s²]
 * @returns {{qt:number, phi:number, theta:number}} total load [m²/s], Φ [-], θ [-]
 */
export function engelundHansen({ U, tau, d50, rhoS = 2650, rho = 1025, g = G }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  need(Number.isFinite(tau) && tau >= 0, 'tau must be >= 0');
  pos(d50, 'd50');
  need(rhoS > rho, 'rhoS must exceed rho');
  if (U === 0 || tau === 0) return { qt: 0, phi: 0, theta: 0 };
  const s1 = rhoS / rho - 1;
  const theta = tau / (rho * g * s1 * d50);
  const f = (2 * tau) / (rho * U * U);
  const phi = (0.1 * Math.pow(theta, 2.5)) / f;
  return { qt: phi * Math.sqrt(g * s1 * d50 ** 3), phi, theta };
}

function vanRijnT({ U, h, d50, d90, rhoS, rho, nu, g }) {
  const s1 = rhoS / rho - 1;
  const Dstar = d50 * Math.cbrt((g * s1) / (nu * nu));
  const thetaCr = criticalShields(Dstar);
  const uCr2 = thetaCr * g * s1 * d50;
  need((12 * h) / (3 * d90) > 1.0001, 'h must be much larger than d90');
  const Cg = 18 * Math.log10((12 * h) / (3 * d90));
  const uGrain = (Math.sqrt(g) / Cg) * U;
  return { s1, Dstar, thetaCr, uCr: Math.sqrt(uCr2), T: Math.max(0, (uGrain * uGrain - uCr2) / uCr2), Cgrain: Cg };
}

/**
 * van Rijn (1984a) bed load: qb = 0.053 [(s − 1) g]^0.5 d50^1.5 T^2.1 / D*^0.3,
 *   T = (u*'² − u*cr²)/u*cr²,  u*' = (g^0.5/C') U,  C' = 18 log10(12 h/(3 d90)).
 * The threshold u*cr is taken from the Soulsby–Whitehouse curve (an algebraic fit of the Shields curve
 * that van Rijn used graphically).
 * Reference: van Rijn (1984a) J. Hydraul. Eng. 110(10). Validity: 0.2–2 mm, T < 3 (the formula
 * over-predicts for larger T), h > 1 m in the calibration.
 * @param {object} o
 * @param {number} o.U depth-averaged velocity [m/s]
 * @param {number} o.h depth [m]
 * @param {number} o.d50 [m]
 * @param {number} [o.d90] [m], default 2 d50
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.nu] [m²/s]
 * @param {number} [o.g] [m/s²]
 * @returns {{qb:number, T:number, Dstar:number, uStarGrain:number, uStarCr:number}} bed load [m²/s],
 *   transport stage [-], D* [-], grain shear velocity and threshold [m/s]
 */
export function vanRijnBedLoad({ U, h, d50, d90, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  pos(h, 'h');
  pos(d50, 'd50');
  const D90 = d90 !== undefined ? d90 : 2 * d50;
  const r = vanRijnT({ U, h, d50, d90: D90, rhoS, rho, nu, g });
  const qb = (0.053 * Math.sqrt(r.s1 * g) * Math.pow(d50, 1.5) * Math.pow(r.T, 2.1)) / Math.pow(r.Dstar, 0.3);
  return { qb, T: r.T, Dstar: r.Dstar, uStarGrain: (Math.sqrt(g) / r.Cgrain) * U, uStarCr: r.uCr };
}

/**
 * van Rijn (1984b) suspended load by the reference-concentration method:
 *   ca = 0.015 (d50/a) T^1.5 / D*^0.3   (volumetric, at height a above the bed)
 *   ds = d50 [1 + 0.011 (σs − 1)(T − 25)] for T < 25 (representative suspended size), else d50
 *   Z = ws/(β κ u*),  β = 1 + 2 (ws/u*)²  (limited to 2),  φ = 2.5 (ws/u*)^0.8 (ca/c0)^0.4,  c0 = 0.65
 *   Z' = Z + φ,  F = [(a/h)^Z' − (a/h)^1.2] / [(1 − a/h)^Z' (1.2 − Z')],  qs = F U h ca
 * with the overall shear velocity u* = (g^0.5/C) U, C = 18 log10(12 h/ks).
 * Reference: van Rijn (1984b) J. Hydraul. Eng. 110(11).
 * Validity: 0.1–0.5 mm sand, 0.01 <= a/h <= 0.1, 0.1 < ws/u* < 1 for β.
 * @param {object} o
 * @param {number} o.U depth-averaged velocity [m/s]
 * @param {number} o.h depth [m]
 * @param {number} o.d50 [m]
 * @param {number} [o.d90] [m], default 2 d50
 * @param {number} [o.ks] overall (form) roughness height [m], default 3 d90
 * @param {number} [o.a] reference level [m], default max(ks, 0.01 h)
 * @param {number} [o.sigmaS] geometric standard deviation of the bed material [-], default 2.0
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.nu] [m²/s]
 * @param {number} [o.g] [m/s²]
 * @returns {{qs:number, ca:number, a:number, F:number, Z:number, Zprime:number, ds:number, ws:number,
 *   uStar:number, T:number}} suspended load [m²/s], reference concentration [-], a [m], F [-], suspension
 *   numbers [-], ds [m], ws [m/s], u* [m/s], T [-]
 */
export function vanRijnSuspendedLoad({ U, h, d50, d90, ks, a, sigmaS = 2.0, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  pos(h, 'h');
  pos(d50, 'd50');
  need(sigmaS >= 1, 'sigmaS must be >= 1');
  const D90 = d90 !== undefined ? d90 : 2 * d50;
  const KS = ks !== undefined ? ks : 3 * D90;
  pos(KS, 'ks');
  const aa = a !== undefined ? a : Math.max(KS, 0.01 * h);
  need(aa > 0 && aa < h, 'a must lie in (0, h)');
  const r = vanRijnT({ U, h, d50, d90: D90, rhoS, rho, nu, g });
  const ds = r.T < 25 ? d50 * (1 + 0.011 * (sigmaS - 1) * (r.T - 25)) : d50;
  const ws = settlingVelocityVanRijn({ d: ds, rhoS, rho, nu, g });
  const C = 18 * Math.log10((12 * h) / KS);
  const uStar = (Math.sqrt(g) / C) * U;
  if (r.T === 0 || uStar === 0) return { qs: 0, ca: 0, a: aa, F: 0, Z: Infinity, Zprime: Infinity, ds, ws, uStar, T: r.T };
  const ca = Math.min((0.015 * (d50 / aa) * Math.pow(r.T, 1.5)) / Math.pow(r.Dstar, 0.3), 0.65);
  const beta = Math.min(2, 1 + 2 * (ws / uStar) ** 2);
  const Z = ws / (beta * KAPPA * uStar);
  const Zp = Z + 2.5 * Math.pow(ws / uStar, 0.8) * Math.pow(ca / 0.65, 0.4);
  const ah = aa / h;
  let F;
  if (Math.abs(Zp - 1.2) < 1e-9) F = (-Math.pow(ah, 1.2) * Math.log(ah)) / Math.pow(1 - ah, 1.2);
  else F = (Math.pow(ah, Zp) - Math.pow(ah, 1.2)) / (Math.pow(1 - ah, Zp) * (1.2 - Zp));
  return { qs: F * U * h * ca, ca, a: aa, F, Z, Zprime: Zp, ds, ws, uStar, T: r.T };
}

function vanRijnUcr(h, d50, d90, strict = true) {
  if (strict) need(d50 >= 0.1e-3 && d50 <= 2e-3, 'd50 must lie in [0.1 mm, 2 mm] for the van Rijn threshold velocity');
  return d50 <= 0.5e-3 ? 0.19 * Math.pow(d50, 0.1) * Math.log10((4 * h) / d90) : 8.5 * Math.pow(d50, 0.6) * Math.log10((4 * h) / d90);
}

/**
 * van Rijn (1984) simplified current-alone transport formulae:
 *   qb = 0.005 U h {(U − Ucr)/[(s − 1) g d50]^0.5}^2.4 (d50/h)^1.2
 *   qs = 0.012 U h {(U − Ucr)/[(s − 1) g d50]^0.5}^2.4 (d50/h) D*^−0.6
 *   Ucr = 0.19 d50^0.1 log10(4h/d90) for 0.1 <= d50 <= 0.5 mm;  8.5 d50^0.6 log10(4h/d90) for 0.5–2 mm
 * Reference: van Rijn (1984c) J. Hydraul. Eng. 110(12); Soulsby (1997) Eq. 133a–d.
 * Validity: h = 1–20 m, U = 0.5–5 m/s, d50 = 0.1–2 mm, fresh water at 15 °C in the calibration.
 * By default a grain size outside 0.1–2 mm is rejected. With `strict: false` the formulae are evaluated anyway —
 * the fine-sand threshold expression below 0.1 mm and the coarse-sand expression above 2 mm — and the result
 * carries `outsideValidity: true`: an extrapolation beyond the calibration data (cohesion below about 0.06 mm
 * and the gravel regime above 2 mm are not represented), to be used for screening only. 4h/d90 must exceed 1
 * in either mode.
 * @param {object} o
 * @param {number} o.U [m/s]
 * @param {number} o.h [m]
 * @param {number} o.d50 [m]
 * @param {number} [o.d90] [m], default 2 d50
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.nu] [m²/s]
 * @param {number} [o.g] [m/s²]
 * @param {boolean} [o.strict] reject d50 outside 0.1–2 mm (default true); false computes with a warning flag
 * @returns {{qb:number, qs:number, qt:number, Ucr:number, outsideValidity:boolean}} [m²/s], threshold velocity
 *   [m/s]; outsideValidity is true when d50 lies outside 0.1–2 mm (possible only with strict: false)
 */
export function vanRijnSimplified({ U, h, d50, d90, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G, strict = true }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  pos(h, 'h');
  pos(d50, 'd50');
  const D90 = d90 !== undefined ? d90 : 2 * d50;
  pos(D90, 'd90');
  const s1 = rhoS / rho - 1;
  const outsideValidity = !(d50 >= 0.1e-3 && d50 <= 2e-3);
  need(4 * h > D90, 'h must exceed d90/4 for the van Rijn threshold velocity');
  const Ucr = vanRijnUcr(h, d50, D90, strict);
  if (U <= Ucr) return { qb: 0, qs: 0, qt: 0, Ucr, outsideValidity };
  const Dstar = dimensionlessGrainSize({ d: d50, rhoS, rho, nu, g });
  const m = Math.pow((U - Ucr) / Math.sqrt(s1 * g * d50), 2.4);
  const qb = 0.005 * U * h * m * Math.pow(d50 / h, 1.2);
  const qs = 0.012 * U * h * m * (d50 / h) * Math.pow(Dstar, -0.6);
  return { qb, qs, qt: qb + qs, Ucr, outsideValidity };
}

/**
 * Soulsby–van Rijn total load for combined waves and currents:
 *   qt = As U [(U² + (0.018/CD) Urms²)^0.5 − Ucr]^2.4 (1 − 1.6 tan β)
 *   Asb = 0.005 h (d50/h)^1.2 / [(s − 1) g d50]^1.2,  Ass = 0.012 d50 D*^−0.6 / [(s − 1) g d50]^1.2
 *   CD = [0.40/(ln(h/z0) − 1)]², z0 = 0.006 m (rippled bed), Ucr as van Rijn (1984).
 * Reference: Soulsby (1997) Eq. 136a–g. Validity: rippled sand beds, h = 1–20 m, U = 0.5–5 m/s,
 * d50 = 0.1–2 mm; z0 should stay 0.006 m (it is part of the calibration).
 * @param {object} o
 * @param {number} o.U depth-averaged current speed [m/s]
 * @param {number} o.Urms root-mean-square wave orbital velocity at the bed [m/s]
 * @param {number} o.h depth [m]
 * @param {number} o.d50 [m]
 * @param {number} [o.d90] [m], default 2 d50
 * @param {number} [o.bedSlope] tan β of the bed in the streamwise direction, positive up-slope [-]
 * @param {number} [o.z0] [m], default 0.006
 * @param {number} [o.rhoS] [kg/m³]
 * @param {number} [o.rho] [kg/m³]
 * @param {number} [o.nu] [m²/s]
 * @param {number} [o.g] [m/s²]
 * @returns {{qt:number, qb:number, qs:number, Ucr:number, Cd:number}} [m²/s], [m/s], [-]
 */
export function soulsbyVanRijn({ U, Urms, h, d50, d90, bedSlope = 0, z0 = 0.006, rhoS = 2650, rho = 1025, nu = 1.0e-6, g = G }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  need(Number.isFinite(Urms) && Urms >= 0, 'Urms must be >= 0');
  pos(h, 'h');
  pos(d50, 'd50');
  pos(z0, 'z0');
  need(h / z0 > Math.E * 1.0001, 'h/z0 must exceed e');
  const D90 = d90 !== undefined ? d90 : 2 * d50;
  const s1 = rhoS / rho - 1;
  const Ucr = vanRijnUcr(h, d50, D90);
  const Cd = (0.4 / (Math.log(h / z0) - 1)) ** 2;
  const Ue = Math.sqrt(U * U + (0.018 / Cd) * Urms * Urms);
  if (Ue <= Ucr || U === 0) return { qt: 0, qb: 0, qs: 0, Ucr, Cd };
  const Dstar = dimensionlessGrainSize({ d: d50, rhoS, rho, nu, g });
  const den = Math.pow(s1 * g * d50, 1.2);
  const Asb = (0.005 * h * Math.pow(d50 / h, 1.2)) / den;
  const Ass = (0.012 * d50 * Math.pow(Dstar, -0.6)) / den;
  const common = U * Math.pow(Ue - Ucr, 2.4) * (1 - 1.6 * bedSlope);
  return { qt: (Asb + Ass) * common, qb: Asb * common, qs: Ass * common, Ucr, Cd };
}

/**
 * Rouse suspended-sediment concentration profile for a parabolic eddy diffusivity:
 *   c(z)/ca = [((h − z)/z) (a/(h − a))]^Z,  Z = ws/(κ u*)  (Rouse number), κ = 0.40.
 * Reference: Rouse (1937) Trans. ASCE 102; Soulsby (1997) §8.4. Validity: steady uniform current,
 * a <= z < h.
 * @param {object} o
 * @param {number} o.z height above the bed [m]
 * @param {number} o.h depth [m]
 * @param {number} o.a reference height [m]
 * @param {number} o.ca reference concentration at z = a [any unit; the result has the same unit]
 * @param {number} o.ws settling velocity [m/s]
 * @param {number} o.uStar shear velocity [m/s]
 * @returns {{c:number, rouse:number}} concentration [unit of ca], Rouse number [-]
 */
export function rouseProfile({ z, h, a, ca, ws, uStar }) {
  pos(h, 'h');
  need(a > 0 && a < h, 'a must lie in (0, h)');
  need(z >= a && z <= h, 'z must lie in [a, h]');
  pos(ws, 'ws');
  pos(uStar, 'uStar');
  need(Number.isFinite(ca) && ca >= 0, 'ca must be >= 0');
  const rouse = ws / (KAPPA * uStar);
  return { c: ca * Math.pow(((h - z) / z) * (a / (h - a)), rouse), rouse };
}

/**
 * Isbash (1936) stone size against current: U = C [2 g Δ D]^(1/2)  →  D = U² / (2 g Δ C²).
 * C = 0.86 for exposed stones (high turbulence), 1.20 for embedded stones (low turbulence).
 * Reference: Isbash (1936) Trans. 2nd Congress on Large Dams; SPM (1984) §7.III.7; USACE EM 1110-2-1601.
 * Validity: U is the velocity acting on the stone; no depth or slope effect.
 * @param {object} o
 * @param {number} o.U velocity at the stone [m/s]
 * @param {number} o.delta relative buoyant density Δ [-]
 * @param {number} [o.C] Isbash coefficient [-], default 0.86
 * @param {number} [o.g] [m/s²]
 * @returns {number} stone diameter D (equivalent sphere, ≈ D50) [m]
 */
export function isbashStoneSize({ U, delta, C = 0.86, g = G }) {
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  pos(delta, 'delta');
  pos(C, 'C');
  return (U * U) / (2 * g * delta * C * C);
}

/**
 * Shields-based rock size against a steady current: solves θ = τ/(ρ g Δ d) = ψcr with
 * τ = ρ g U²/C², C = 18 log10(12 h/ks), ks = ksFactor·d, i.e. d = U²/(C(d)² Δ ψcr), optionally divided by a
 * side-slope reduction ksl = [1 − sin²α/sin²φ]^(1/2).
 * Reference: Shields (1936); Rock Manual (CIRIA C683, 2007) §5.2.1.3 and §5.2.3.1 (ψcr = 0.03–0.035 for
 * no movement, 0.05–0.055 for limited movement; ks = 1–3 Dn50). Validity: uniform rough turbulent flow,
 * h/d large enough for a log profile.
 * @param {object} o
 * @param {number} o.U depth-averaged velocity [m/s]
 * @param {number} o.h depth [m]
 * @param {number} o.delta Δ [-]
 * @param {number} [o.psiCr] critical Shields parameter [-], default 0.035
 * @param {number} [o.ksFactor] ks/d [-], default 2
 * @param {number} [o.slopeDeg] side-slope angle α [deg], default 0
 * @param {number} [o.reposeDeg] angle of repose φ [deg], default 40
 * @param {number} [o.g] [m/s²]
 * @returns {{d:number, chezy:number, tau:number, ksl:number}} stone size [m], Chézy coefficient
 *   [m^(1/2)/s], bed shear stress per unit density τ/ρ [m²/s²], slope factor [-]
 */
export function shieldsRockSize({ U, h, delta, psiCr = 0.035, ksFactor = 2, slopeDeg = 0, reposeDeg = 40, g = G }) {
  pos(U, 'U');
  pos(h, 'h');
  pos(delta, 'delta');
  pos(psiCr, 'psiCr');
  pos(ksFactor, 'ksFactor');
  need(slopeDeg >= 0 && slopeDeg < reposeDeg, 'slopeDeg must lie in [0, reposeDeg)');
  const sa = Math.sin((slopeDeg * Math.PI) / 180);
  const sp = Math.sin((reposeDeg * Math.PI) / 180);
  const ksl = Math.sqrt(1 - (sa * sa) / (sp * sp));
  const chezy = (d) => 18 * Math.log10((12 * h) / (ksFactor * d));
  const f = (d) => d - (U * U) / (chezy(d) ** 2 * delta * psiCr * ksl);
  const dHi = (12 * h) / (ksFactor * 2); // 12h/ks = 2: the log law is meaningless beyond this
  need(f(dHi) > 0, 'U is too large for the depth h: the required stone approaches the water depth');
  const d = brent(f, 1e-9, dHi, { tol: 1e-12 });
  const C = chezy(d);
  return { d, chezy: C, tau: (g * U * U) / (C * C), ksl };
}
