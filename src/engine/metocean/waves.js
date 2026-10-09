// Regular-wave theories: linear (Airy), Stokes 2nd and 5th order, cnoidal, solitary; dispersion, energy,
// breaking limits, regime classification, shoaling and refraction.
//
// Conventions: x horizontal in the direction of propagation (m), z vertical, positive upwards from the still-water
// level (m), d still-water depth (m, positive; Infinity allowed where stated), t time (s), H wave height (m),
// T wave period as seen by a fixed observer (s). Phase theta = k x - omega t, crest at theta = 0.
// Accelerations returned as ax, az are local (Eulerian, d/dt at a fixed point) unless a field name says "Total".

import { brent } from '../math/roots.js';

/** Standard gravity, m/s². */
export const G = 9.80665;
/** Default sea-water density, kg/m³. */
export const RHO_SEA = 1025;

const TWO_PI = 2 * Math.PI;

function reqPos(name, v) {
  if (!(typeof v === 'number') || !(v > 0) || Number.isNaN(v)) throw new RangeError(`${name} must be a positive number`);
}
function reqFinitePos(name, v) {
  if (!(typeof v === 'number') || !(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a positive finite number`);
}
function reqFinite(name, v) {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new RangeError(`${name} must be a finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Dispersion                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

// Intrinsic angular frequency sqrt(g k tanh kd) and its k-derivative (the relative group velocity).
function sigmaOf(k, d, g) {
  const kd = k * d;
  const th = Number.isFinite(kd) ? Math.tanh(kd) : 1;
  const sig = Math.sqrt(g * k * th);
  let n;
  if (!Number.isFinite(kd) || kd > 350) n = 0.5;
  else if (kd < 1e-8) n = 1;
  else n = 0.5 * (1 + (2 * kd) / Math.sinh(2 * kd));
  return { sig, cg: k > 0 ? (n * sig) / k : Math.sqrt(g * d), n, th };
}

/**
 * Wavenumber from the exact linear dispersion relation, optionally Doppler-shifted by a uniform current:
 * (omega - k U)² = g k tanh(k d), omega = 2π/T.
 * Reference: Dean & Dalrymple, Water Wave Mechanics for Engineers and Scientists (1991) §3.4; current form
 * per DNV-RP-C205 (2019) §3.2.2 / Peregrine (1976).
 * Validity: T > 0, d > 0 (d = Infinity gives deep water). For an opposing current (U < 0) the wave is blocked
 * when no solution exists on the branch with relative group velocity > |U|; this throws RangeError.
 * @param {object} o
 * @param {number} o.T wave period in the fixed frame, s
 * @param {number} o.d water depth, m (Infinity allowed)
 * @param {number} [o.U=0] uniform current component in the direction of wave propagation, m/s (negative opposes)
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @returns {number} wavenumber k, rad/m
 */
export function wavenumber({ T, d, U = 0, g = G }) {
  reqFinitePos('T', T); reqPos('d', d); reqFinite('U', U); reqFinitePos('g', g);
  const om = TWO_PI / T;
  const f = (k) => sigmaOf(k, d, g).sig + k * U - om;
  const fp = (k) => sigmaOf(k, d, g).cg + U;
  // Upper bound of the admissible branch
  let kHi = Infinity;
  if (U < 0) {
    if (Math.sqrt(g * d) + U <= 0) throw new RangeError('U: opposing current exceeds the long-wave speed, wave is blocked');
    // relative group velocity decreases monotonically with k: locate cg(k) = -U
    let a = 0;
    let b = (g / (4 * U * U)) * 1.0000001 + 1e-12; // deep-water stopping wavenumber g/(4U²) bounds it from above
    while (fp(b) > 0) b *= 2;
    for (let i = 0; i < 200; i++) {
      const m = 0.5 * (a + b);
      if (fp(m) > 0) a = m; else b = m;
      if (b - a <= 1e-15 * b) break;
    }
    kHi = 0.5 * (a + b);
    if (f(kHi) < 0) throw new RangeError('U: wave is blocked by the opposing current (no solution of the dispersion relation)');
  }
  // Starting guess: Eckart-type explicit approximation on the intrinsic frequency, then safeguarded Newton.
  const k0 = (om * om) / g;
  const kdeep = k0;
  let k = Number.isFinite(d) ? k0 / Math.sqrt(Math.tanh(k0 * d)) : kdeep;
  if (U !== 0) {
    // one fixed-point correction for the Doppler shift
    const c = om / k;
    const fac = c / (c + U);
    if (fac > 0 && Number.isFinite(fac)) k *= fac;
  }
  let lo = 0;
  let hi = kHi;
  if (!(k < hi)) k = 0.5 * hi;
  for (let it = 0; it < 100; it++) {
    const fk = f(k);
    if (fk > 0) hi = Math.min(hi, k); else lo = Math.max(lo, k);
    let kn = k - fk / fp(k);
    if (!(kn > lo) || !(kn < hi)) kn = Number.isFinite(hi) ? 0.5 * (lo + hi) : 2 * k;
    if (Math.abs(kn - k) <= 4e-16 * kn) { k = kn; break; }
    k = kn;
  }
  // Polish on the squared form, which is what callers check residuals against.
  for (let it = 0; it < 3; it++) {
    const kd = k * d;
    const th = Number.isFinite(kd) ? Math.tanh(kd) : 1;
    const sech2 = Number.isFinite(kd) && kd < 350 ? 1 / (Math.cosh(kd) * Math.cosh(kd)) : 0;
    const r = g * k * th - (om - k * U) * (om - k * U);
    const dr = g * th + (Number.isFinite(d) ? g * kd * sech2 : 0) + 2 * U * (om - k * U);
    const kn = k - r / dr;
    if (!(kn > 0) || !(kn <= kHi * (1 + 1e-12))) break;
    k = kn;
  }
  return k;
}

/**
 * Linear wave properties: wavelength, celerity, group velocity, energy density and energy flux.
 * Reference: Dean & Dalrymple (1991) ch. 3–4; DNV-RP-C205 (2019) §3.2.
 * Validity: small-amplitude theory. With a current U the celerity and group velocity are absolute (fixed frame)
 * and `cRel`, `cgRel` are relative to the water; the energy flux uses the absolute group velocity.
 * @param {object} o
 * @param {number} o.T period, s
 * @param {number} o.d depth, m (Infinity allowed)
 * @param {number} [o.H] wave height, m (needed for energy, steepness and Ursell fields; otherwise they are NaN)
 * @param {number} [o.U=0] following current, m/s
 * @param {number} [o.rho=1025] water density, kg/m³
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @returns {{k:number,L:number,omega:number,omegaRel:number,c:number,cRel:number,cg:number,cgRel:number,n:number,
 *   kd:number,L0:number,E:number,P:number,steepness:number,ursell:number,regime:string}}
 *   k rad/m; L, L0 (deep-water length gT²/2π) m; omega, omegaRel rad/s; c, cRel, cg, cgRel m/s; n = cgRel/cRel (-);
 *   E energy per unit sea-surface area J/m²; P energy flux per unit crest width W/m; steepness H/L (-);
 *   ursell H L²/d³ (-); regime 'deep' | 'intermediate' | 'shallow'.
 */
export function waveProperties({ T, d, H = NaN, U = 0, rho = RHO_SEA, g = G }) {
  const k = wavenumber({ T, d, U, g });
  const omega = TWO_PI / T;
  const s = sigmaOf(k, d, g);
  const L = TWO_PI / k;
  const E = (rho * g * H * H) / 8;
  const cg = s.cg + U;
  return {
    k, L, omega, omegaRel: s.sig, c: omega / k, cRel: s.sig / k, cg, cgRel: s.cg, n: s.n, kd: k * d,
    L0: (g * T * T) / TWO_PI, E, P: E * cg, steepness: H / L,
    ursell: Number.isFinite(d) ? (H * L * L) / (d * d * d) : 0, regime: depthRegime(d / L),
  };
}

/**
 * Wave length from period and depth (linear dispersion). Reference and validity as `wavenumber`.
 * @param {number} T period, s
 * @param {number} d depth, m
 * @param {number} [g=9.80665] gravity, m/s²
 * @returns {number} wavelength, m
 */
export function wavelength(T, d, g = G) { return TWO_PI / wavenumber({ T, d, g }); }

/**
 * Phase speed (celerity) from period and depth (linear dispersion, no current).
 * @param {number} T period, s
 * @param {number} d depth, m
 * @param {number} [g=9.80665] gravity, m/s²
 * @returns {number} celerity, m/s
 */
export function celerity(T, d, g = G) { return TWO_PI / (T * wavenumber({ T, d, g })); }

/**
 * Group velocity cg = n c, n = ½ (1 + 2kd / sinh 2kd) (linear theory, no current).
 * Reference: Dean & Dalrymple (1991) Eq. 4.82.
 * @param {number} T period, s
 * @param {number} d depth, m
 * @param {number} [g=9.80665] gravity, m/s²
 * @returns {number} group velocity, m/s
 */
export function groupVelocity(T, d, g = G) { return sigmaOf(wavenumber({ T, d, g }), d, g).cg; }

/**
 * Mean wave energy per unit sea-surface area, E = ρ g H² / 8 (kinetic + potential, linear theory).
 * @param {number} H wave height, m
 * @param {number} [rho=1025] density, kg/m³
 * @param {number} [g=9.80665] gravity, m/s²
 * @returns {number} energy density, J/m²
 */
export function energyDensity(H, rho = RHO_SEA, g = G) { return (rho * g * H * H) / 8; }

/**
 * Mean wave energy flux (wave power) per unit crest width, P = E cg.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.T period, s
 * @param {number} o.d depth, m
 * @param {number} [o.rho=1025] density, kg/m³
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @returns {number} energy flux, W/m
 */
export function energyFlux({ H, T, d, rho = RHO_SEA, g = G }) {
  return energyDensity(H, rho, g) * groupVelocity(T, d, g);
}

/**
 * Relative-depth regime: shallow d/L < 1/20, deep d/L > 1/2, otherwise intermediate.
 * Reference: Shore Protection Manual (1984) Table 2-1.
 * @param {number} dOverL depth / wavelength, -
 * @returns {'deep'|'intermediate'|'shallow'}
 */
export function depthRegime(dOverL) {
  if (!(dOverL > 0)) throw new RangeError('dOverL must be positive');
  if (dOverL > 0.5) return 'deep';
  if (dOverL < 0.05) return 'shallow';
  return 'intermediate';
}

/**
 * Ursell number Ur = H L² / d³.
 * Reference: Ursell (1953); DNV-RP-C205 (2019) §3.2.1.
 * @param {number} H wave height, m
 * @param {number} L wavelength, m
 * @param {number} d depth, m
 * @returns {number} Ursell number, -
 */
export function ursellNumber(H, L, d) {
  reqFinitePos('H', H); reqFinitePos('L', L); reqFinitePos('d', d);
  return (H * L * L) / (d * d * d);
}

/**
 * Wave steepness H / L.
 * @param {number} H wave height, m
 * @param {number} L wavelength, m
 * @returns {number} steepness, -
 */
export function steepness(H, L) { reqFinitePos('L', L); return H / L; }

/* ------------------------------------------------------------------------------------------------------------ */
/* Depth-decay functions, written with decaying exponentials so that deep water and d = Infinity are exact       */
/* ------------------------------------------------------------------------------------------------------------ */

// Returns cosh k(z+d)/sinh kd, sinh k(z+d)/sinh kd, cosh k(z+d)/cosh kd.
function decay(k, z, d, out) {
  const e = Math.exp(k * z);
  const q = Number.isFinite(d) ? Math.exp(-2 * k * d) : 0;
  const r = Number.isFinite(d) ? Math.exp(-2 * k * (z + d)) : 0;
  out.cs = (e * (1 + r)) / (1 - q);
  out.ss = (e * (1 - r)) / (1 - q);
  out.cc = (e * (1 + r)) / (1 + q);
  return out;
}

function stretchedDecay(k, z, d, eta, stretching, out) {
  switch (stretching) {
    case 'none':
      return decay(k, z, d, out);
    case 'wheeler': {
      // Wheeler (1970): the SWL profile is mapped onto the instantaneous water column.
      if (!Number.isFinite(d)) return decay(k, z - eta, d, out);
      return decay(k, ((z - eta) * d) / (d + eta), d, out);
    }
    case 'vertical':
      return decay(k, Math.min(z, 0), d, out);
    case 'linear': {
      if (z <= 0) return decay(k, z, d, out);
      decay(k, 0, d, out);
      const c0 = out.cs; const s0 = out.ss; // s0 = 1
      const th = Number.isFinite(d) ? Math.tanh(k * d) : 1;
      out.cs = c0 + z * k * s0;
      out.ss = s0 + z * k * c0;
      out.cc = 1 + z * k * th;
      return out;
    }
    default:
      throw new RangeError("stretching must be 'none', 'wheeler', 'linear' or 'vertical'");
  }
}

/**
 * Airy (linear) wave surface elevation and kinematics at a point, with optional kinematic stretching to the
 * instantaneous free surface.
 *   eta = (H/2) cos θ,  u = (H/2) σ cosh k(z+d)/sinh kd · cos θ,  w = (H/2) σ sinh k(z+d)/sinh kd · sin θ,
 *   ax = ∂u/∂t, az = ∂w/∂t,  pDyn = ρ g (H/2) cosh k(z+d)/cosh kd · cos θ,   θ = k x − ω t + phase.
 * Stretching options (applied to the depth-decay functions):
 *   'none'     hyperbolic functions evaluated at z as they are (strictly valid for z ≤ 0);
 *   'wheeler'  Wheeler (1970) stretching, z' = (z − η) d / (d + η);
 *   'linear'   linear (Taylor) extrapolation above the still-water level, f(z) = f(0) + z f'(0), z > 0;
 *   'vertical' vertical extrapolation: values above the still-water level equal those at z = 0.
 * Reference: Dean & Dalrymple (1991) ch. 4; DNV-RP-C205 (2019) §3.2.2 and §3.5.2 (stretching);
 * Wheeler, J.D. (1970) J. Petroleum Technology 22(3).
 * Validity: small steepness and Ursell number; −d ≤ z. Points above the instantaneous surface return zero
 * kinematics with `wet = false` unless `clipAboveSurface` is false. With a current U, σ is the intrinsic
 * frequency and the returned velocities are wave-induced only (add the current yourself).
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.T period (fixed frame), s
 * @param {number} o.d depth, m (Infinity allowed)
 * @param {number} [o.x=0] horizontal position, m
 * @param {number} [o.z=0] elevation above still-water level, m
 * @param {number} [o.t=0] time, s
 * @param {number} [o.phase=0] additional phase, rad
 * @param {number} [o.U=0] uniform following current used in the dispersion relation, m/s
 * @param {'none'|'wheeler'|'linear'|'vertical'} [o.stretching='none']
 * @param {boolean} [o.clipAboveSurface=true] return zeros above the instantaneous surface
 * @param {number} [o.rho=1025] density, kg/m³
 * @param {number} [o.g=9.80665] gravity, m/s²
 * @param {number} [o.k] wavenumber if already known, rad/m (skips the dispersion solve)
 * @returns {{eta:number,u:number,w:number,ax:number,az:number,pDyn:number,k:number,omega:number,theta:number,wet:boolean}}
 *   eta m; u, w m/s; ax, az m/s²; pDyn dynamic pressure (excluding −ρ g z) Pa; k rad/m; omega rad/s; theta rad.
 */
export function airyKinematics({ H, T, d, x = 0, z = 0, t = 0, phase = 0, U = 0, stretching = 'none',
  clipAboveSurface = true, rho = RHO_SEA, g = G, k }) {
  if (!(H >= 0)) throw new RangeError('H must be non-negative');
  reqFinite('z', z);
  const kk = k ?? wavenumber({ T, d, U, g });
  const omega = TWO_PI / T;
  const sig = omega - kk * U;
  const a = 0.5 * H;
  const theta = kk * x - omega * t + phase;
  const c = Math.cos(theta); const s = Math.sin(theta);
  const eta = a * c;
  if (Number.isFinite(d) && z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  if (clipAboveSurface && z > eta + 1e-12 * Math.max(1, a)) {
    return { eta, u: 0, w: 0, ax: 0, az: 0, pDyn: 0, k: kk, omega, theta, wet: false };
  }
  const f = stretchedDecay(kk, z, d, eta, stretching, { cs: 0, ss: 0, cc: 0 });
  return {
    eta,
    u: a * sig * f.cs * c,
    w: a * sig * f.ss * s,
    ax: a * sig * omega * f.cs * s,
    az: -a * sig * omega * f.ss * c,
    pDyn: rho * g * a * f.cc * c,
    k: kk, omega, theta, wet: true,
  };
}

/**
 * Stokes second-order wave: surface elevation, kinematics and dynamic pressure at a point.
 *   eta = a cos θ + (k a²/4) cosh kd (2 + cosh 2kd)/sinh³ kd · cos 2θ
 *   u   = a ω cosh ks/sinh kd · cos θ + ¾ a² ω k cosh 2ks/sinh⁴ kd · cos 2θ          (s = z + d)
 *   w   = a ω sinh ks/sinh kd · sin θ + ¾ a² ω k sinh 2ks/sinh⁴ kd · sin 2θ
 *   p   = ρ g a cosh ks/cosh kd cos θ + ¾ ρ g (2 k a²)/sinh 2kd [cosh 2ks/sinh² kd − 1/3] cos 2θ
 *         − ¼ ρ g (2 k a²)/sinh 2kd [cosh 2ks − 1]
 * with a = H/2 and the linear dispersion relation (the frequency correction appears at third order).
 * Reference: Dean & Dalrymple (1991) §11.2; Shore Protection Manual (1984) Eqs. 2-45 to 2-52;
 * DNV-RP-C205 (2019) §3.2.3.
 * Validity: Ursell number below about 26 and steepness well below breaking; the second harmonic must stay
 * smaller than the first (`secondaryCrest` is true when a spurious bump appears in the trough).
 * Kinematics are evaluated at z for z ≤ eta (no stretching); above the surface zeros are returned.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.T period, s
 * @param {number} o.d depth, m (Infinity allowed)
 * @param {number} [o.x=0] m
 * @param {number} [o.z=0] m above still-water level
 * @param {number} [o.t=0] s
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{eta:number,eta1:number,eta2:number,u:number,w:number,ax:number,az:number,pDyn:number,k:number,
 *   omega:number,stokesDrift:number,secondaryCrest:boolean,wet:boolean}}
 *   eta, eta1, eta2 m; u, w m/s; ax, az m/s² (local); pDyn Pa; stokesDrift Lagrangian mean drift velocity at z, m/s.
 */
export function stokes2Kinematics({ H, T, d, x = 0, z = 0, t = 0, rho = RHO_SEA, g = G }) {
  if (!(H >= 0)) throw new RangeError('H must be non-negative');
  const k = wavenumber({ T, d, g });
  const omega = TWO_PI / T;
  const a = 0.5 * H;
  const theta = k * x - omega * t;
  const c1 = Math.cos(theta); const s1 = Math.sin(theta);
  const c2 = Math.cos(2 * theta); const s2 = Math.sin(2 * theta);
  const deep = !Number.isFinite(d) || k * d > 300;
  // second-order surface factor: cosh kd (2 + cosh 2kd) / sinh³ kd  -> 2 in deep water
  let surf;
  let invSinh4E; // helper: 1/sinh⁴ kd handled through exponentials below
  if (deep) surf = 2;
  else { const sh = Math.sinh(k * d); surf = (Math.cosh(k * d) * (2 + Math.cosh(2 * k * d))) / (sh * sh * sh); }
  const eta1 = a * c1;
  const eta2 = 0.25 * k * a * a * surf * c2;
  const eta = eta1 + eta2;
  const secondaryCrest = 0.25 * k * a * surf > 0.25;
  if (Number.isFinite(d) && z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  if (z > eta + 1e-12 * Math.max(1, a)) {
    return { eta, eta1, eta2, u: 0, w: 0, ax: 0, az: 0, pDyn: 0, k, omega, stokesDrift: 0, secondaryCrest, wet: false };
  }
  const f = decay(k, z, d, { cs: 0, ss: 0, cc: 0 });
  // cosh 2ks / sinh⁴ kd and sinh 2ks / sinh⁴ kd; both vanish in deep water (the second-order potential is zero)
  let c2s = 0; let s2s = 0; let cosh2ks = 0; let invSinh2 = 0; let invSinh2kd = 0;
  if (!deep) {
    const sh = Math.sinh(k * d);
    invSinh4E = 1 / (sh * sh * sh * sh);
    cosh2ks = Math.cosh(2 * k * (z + d));
    c2s = cosh2ks * invSinh4E;
    s2s = Math.sinh(2 * k * (z + d)) * invSinh4E;
    invSinh2 = 1 / (sh * sh);
    invSinh2kd = 1 / Math.sinh(2 * k * d);
  }
  const B = 0.75 * a * a * omega * k;
  const u = a * omega * f.cs * c1 + B * c2s * c2;
  const w = a * omega * f.ss * s1 + B * s2s * s2;
  const ax = a * omega * omega * f.cs * s1 + 2 * B * omega * c2s * s2;
  const az = -a * omega * omega * f.ss * c1 - 2 * B * omega * s2s * c2;
  let p2;
  if (deep) {
    // limits: cos 2θ term -> 0; steady term -> −½ ρ g k a² e^{2kz}
    p2 = -0.5 * rho * g * k * a * a * Math.exp(2 * k * z);
  } else {
    const pre = rho * g * 2 * k * a * a * invSinh2kd;
    p2 = 0.75 * pre * (cosh2ks * invSinh2 - 1 / 3) * c2 - 0.25 * pre * (cosh2ks - 1);
  }
  // Stokes drift: a² ω k cosh 2ks / (2 sinh² kd)
  const drift = deep ? a * a * omega * k * Math.exp(2 * k * z) : (a * a * omega * k * cosh2ks * invSinh2) / 2;
  return { eta, eta1, eta2, u, w, ax, az, pDyn: rho * g * a * f.cc * c1 + p2, k, omega, stokesDrift: drift, secondaryCrest, wet: true };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Stokes fifth order (Fenton 1985)                                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Coefficients of the fifth-order Stokes theory as functions of kd, with S = sech 2kd.
 * Reference: Fenton, J.D. (1985) "A fifth-order Stokes theory for steady waves", J. Waterway, Port, Coastal and
 * Ocean Eng. 111(2), 216–234, Table 1, with the C4 correction of Fenton (1990) "Nonlinear wave theories",
 * The Sea vol. 9, Table 1. The test suite verifies that the resulting solution satisfies both nonlinear
 * free-surface conditions with a residual that decays as the sixth power of the steepness.
 * Validity: kd > 0 (accuracy degrades for Ursell number above about 40, i.e. long waves in shallow water).
 * @param {number} kd wavenumber × depth, - (Infinity allowed)
 * @returns {{S:number,C0:number,C2:number,C4:number,A:number[][],B:number[][]}} dimensionless; A[i][j], B[i][j]
 *   are the coefficients A_ij, B_ij (unused entries are 0).
 */
export function stokes5Coefficients(kd) {
  reqPos('kd', kd);
  const fin = Number.isFinite(kd) && kd < 350;
  const S = fin ? 1 / Math.cosh(2 * kd) : 0;
  const sh = fin ? Math.sinh(kd) : Infinity;
  const th = fin ? Math.tanh(kd) : 1;
  const cth = 1 / th;
  const C0 = Math.sqrt(th);
  const m = 1 - S;
  const S2 = S * S; const S3 = S2 * S; const S4 = S3 * S; const S5 = S4 * S; const S6 = S5 * S; const S7 = S6 * S; const S8 = S7 * S;
  const A = [0, 1, 2, 3, 4, 5].map(() => [0, 0, 0, 0, 0, 0]);
  const B = [0, 1, 2, 3, 4, 5].map(() => [0, 0, 0, 0, 0, 0]);
  A[1][1] = 1 / sh;
  A[2][2] = (3 * S2) / (2 * m * m);
  A[3][1] = (-4 - 20 * S + 10 * S2 - 13 * S3) / (8 * sh * m * m * m);
  A[3][3] = (-2 * S2 + 11 * S3) / (8 * sh * m * m * m);
  A[4][2] = (12 * S - 14 * S2 - 264 * S3 - 45 * S4 - 13 * S5) / (24 * m ** 5);
  A[4][4] = (10 * S3 - 174 * S4 + 291 * S5 + 278 * S6) / (48 * (3 + 2 * S) * m ** 5);
  A[5][1] = (-1184 + 32 * S + 13232 * S2 + 21712 * S3 + 20940 * S4 + 12554 * S5 - 500 * S6 - 3341 * S7 - 670 * S8)
    / (64 * sh * (3 + 2 * S) * (4 + S) * m ** 6);
  A[5][3] = (4 * S + 105 * S2 + 198 * S3 - 1376 * S4 - 1302 * S5 - 117 * S6 + 58 * S7)
    / (32 * sh * (3 + 2 * S) * m ** 6);
  A[5][5] = (-6 * S3 + 272 * S4 - 1552 * S5 + 852 * S6 + 2029 * S7 + 430 * S8)
    / (64 * sh * (3 + 2 * S) * (4 + S) * m ** 6);
  B[1][1] = 1;
  B[2][2] = (cth * (1 + 2 * S)) / (2 * m);
  B[3][1] = (-3 * (1 + 3 * S + 3 * S2 + 2 * S3)) / (8 * m * m * m);
  B[3][3] = -B[3][1];
  B[4][2] = (cth * (6 - 26 * S - 182 * S2 - 204 * S3 - 25 * S4 + 26 * S5)) / (6 * (3 + 2 * S) * m ** 4);
  B[4][4] = (cth * (24 + 92 * S + 122 * S2 + 66 * S3 + 67 * S4 + 34 * S5)) / (24 * (3 + 2 * S) * m ** 4);
  B[5][3] = (9 * (132 + 17 * S - 2216 * S2 - 5897 * S3 - 6292 * S4 - 2687 * S5 + 194 * S6 + 467 * S7 + 82 * S8))
    / (128 * (3 + 2 * S) * (4 + S) * m ** 6);
  B[5][5] = (5 * (300 + 1579 * S + 3176 * S2 + 2949 * S3 + 1188 * S4 + 675 * S5 + 1326 * S6 + 827 * S7 + 130 * S8))
    / (384 * (3 + 2 * S) * (4 + S) * m ** 6);
  B[5][1] = -(B[5][3] + B[5][5]);
  const C2 = (C0 * (2 + 7 * S2)) / (4 * m * m);
  const C4 = (C0 * (4 + 32 * S - 116 * S2 - 400 * S3 - 71 * S4 + 146 * S5)) / (32 * m ** 5);
  return { S, C0, C2, C4, A, B };
}

// Assemble the dimensional description of a fifth-order wave for a trial wavenumber.
function stokes5Assemble(k, H, d, g) {
  // beyond kd = 25 every depth-dependent term differs from its deep-water limit by less than 1e-20, and
  // cosh(j kd) would eventually overflow, so the coefficients are evaluated at kd = 25
  const kd = Math.min(k * d, 25);
  const co = stokes5Coefficients(kd);
  const eps = 0.5 * k * H;
  const e2 = eps * eps;
  // mean fluid speed in the frame moving with the wave (= c − Eulerian current)
  const uBar = Math.sqrt(g / k) * (co.C0 + e2 * co.C2 + e2 * e2 * co.C4);
  // harmonic amplitudes: eta = Σ_j etaAmp[j] cos jθ ; velocity potential amplitudes phiAmp[j] (m/s per unit cosh ratio)
  const E = [eps, e2, e2 * eps, e2 * e2, e2 * e2 * eps];
  const etaAmp = new Array(6).fill(0);
  const velAmp = new Array(6).fill(0);
  const vs = co.C0 * Math.sqrt(g / k);
  for (let j = 1; j <= 5; j++) {
    let b = 0; let a = 0;
    for (let i = j; i <= 5; i += 1) { b += E[i - 1] * co.B[i][j]; a += E[i - 1] * co.A[i][j]; }
    etaAmp[j] = b / k;
    // A_ij multiplies cosh(jk y); rescaled to cosh(jk y)/cosh(jkd) so that deep water stays finite
    velAmp[j] = vs * j * a * Math.cosh(j * kd);
  }
  return { k, kd, eps, uBar, etaAmp, velAmp, co };
}

function stokes5Eta(w, theta) {
  let e = 0;
  for (let j = 1; j <= 5; j++) e += w.etaAmp[j] * Math.cos(j * theta);
  return e;
}

// Volume flux per unit width under the wave in the wave frame, from the stream function at the crest.
function stokes5Q(w, d) {
  const etaC = stokes5Eta(w, 0);
  let q = w.uBar * (d + etaC);
  for (let j = 1; j <= 5; j++) {
    const jk = j * w.k;
    // sinh(jk y)/cosh(jk d) at y = d + etaC
    const r = Math.exp(jk * etaC) * (1 - Math.exp(-2 * jk * (d + etaC))) / (1 + Math.exp(-2 * jk * d));
    q -= (w.velAmp[j] / jk) * r;
  }
  return q;
}

function rootNear(f, x0, lo0 = 0.6, hi0 = 1.5) {
  let a = x0 * lo0; let b = x0 * hi0;
  let fa = f(a); let fb = f(b);
  for (let i = 0; i < 40 && !(fa * fb <= 0); i++) {
    a *= 0.8; b *= 1.25; fa = f(a); fb = f(b);
  }
  if (!(fa * fb <= 0)) throw new RangeError('no wavenumber satisfies the specified H, T, d and current');
  return brent(f, a, b, { tol: 1e-14, maxIter: 200 });
}

/**
 * Solve a fifth-order Stokes wave for given height, depth and either period or wavelength.
 * The mean current may be given as the Eulerian mean (time-mean velocity at a point below the trough) or as the
 * Stokes (mass-transport, depth-averaged) drift velocity; exactly one applies, default zero Eulerian current.
 * The mass flux needed for the Stokes-drift closure is evaluated from the stream function at the crest.
 * Reference: Fenton (1985) J. WPCOE 111(2) — Eqs. 14 (surface), 12 (potential), 23/24 (period closure).
 * Validity: finite depth with Ursell number H L²/d³ below about 40 and H below about 90 % of the breaking height;
 * `valid` is false outside that range and the stream-function theory should be used instead.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.d depth, m (Infinity allowed when the current is Eulerian)
 * @param {number} [o.T] period, s (give T or L)
 * @param {number} [o.L] wavelength, m
 * @param {number} [o.currentEulerian=0] Eulerian mean current in the wave direction, m/s
 * @param {number} [o.currentStokes] depth-averaged mass-transport velocity, m/s (overrides currentEulerian)
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{H:number,d:number,g:number,k:number,L:number,T:number,c:number,omega:number,eps:number,uBar:number,
 *   currentEulerian:number,currentStokes:number,Q:number,R:number,etaAmp:number[],velAmp:number[],
 *   crest:number,trough:number,ursell:number,valid:boolean}}
 *   k rad/m; L m; T s; c phase speed m/s; eps = kH/2; uBar mean speed in the wave frame m/s; Q volume flux in the
 *   wave frame m²/s; R Bernoulli constant (bed datum) m²/s²; etaAmp[j] m and velAmp[j] m/s harmonic amplitudes;
 *   crest, trough elevations m.
 */
export function stokes5Wave({ H, d, T, L, currentEulerian = 0, currentStokes, g = G }) {
  reqFinitePos('H', H); reqPos('d', d);
  const useStokes = currentStokes !== undefined && currentStokes !== null;
  if (useStokes && !Number.isFinite(d)) throw new RangeError('currentStokes requires a finite depth');
  const cOf = (k) => {
    const w = stokes5Assemble(k, H, d, g);
    return useStokes ? currentStokes + stokes5Q(w, d) / d : currentEulerian + w.uBar;
  };
  let k;
  if (L !== undefined && L !== null) {
    reqFinitePos('L', L);
    k = TWO_PI / L;
  } else {
    reqFinitePos('T', T);
    const U0 = useStokes ? currentStokes : currentEulerian;
    const k0 = wavenumber({ T, d, U: U0, g });
    k = rootNear((kk) => kk * cOf(kk) * T - TWO_PI, k0);
  }
  const w = stokes5Assemble(k, H, d, g);
  const c = cOf(k);
  const Q = Number.isFinite(d) ? stokes5Q(w, d) : Infinity;
  const cE = c - w.uBar;
  const sol = {
    H, d, g, k, L: TWO_PI / k, T: TWO_PI / (k * c), c, omega: k * c, eps: w.eps, uBar: w.uBar,
    currentEulerian: cE, currentStokes: Number.isFinite(d) ? c - Q / d : cE, Q, R: 0,
    etaAmp: w.etaAmp, velAmp: w.velAmp, crest: stokes5Eta(w, 0), trough: stokes5Eta(w, Math.PI),
    ursell: Number.isFinite(d) ? (H * (TWO_PI / k) ** 2) / d ** 3 : 0, valid: true,
  };
  // Bernoulli constant from the crest (elevation datum at the still-water level so that deep water is finite).
  const kc = stokes5Raw(sol, 0, sol.crest);
  sol.R = 0.5 * ((kc.u - c) ** 2 + kc.w ** 2) + g * sol.crest;
  const Hb = Number.isFinite(d) ? Math.min(0.142 * sol.L * Math.tanh(k * d), 0.78 * d) : 0.142 * sol.L;
  sol.valid = sol.ursell <= 40 && H <= 0.9 * Hb;
  return sol;
}

// Fixed-frame velocity and gradients at phase theta and elevation z (no surface clipping).
function stokes5Raw(sol, theta, z) {
  const { k, d } = sol;
  let u = sol.currentEulerian; let w = 0; let ux = 0; let uz = 0;
  for (let j = 1; j <= 5; j++) {
    const amp = sol.velAmp[j];
    if (amp === 0) continue;
    const jk = j * k;
    const e = Math.exp(jk * z);
    const r = Number.isFinite(d) ? Math.exp(-2 * jk * (z + d)) : 0;
    const q = Number.isFinite(d) ? Math.exp(-2 * jk * d) : 0;
    const ch = (e * (1 + r)) / (1 + q);
    const sh = (e * (1 - r)) / (1 + q);
    const cj = Math.cos(j * theta); const sj = Math.sin(j * theta);
    u += amp * ch * cj;
    w += amp * sh * sj;
    ux -= amp * jk * ch * sj;
    uz += amp * jk * sh * cj;
  }
  return { u, w, ux, uz };
}

/**
 * Surface elevation, velocities, accelerations and pressure of a solved fifth-order Stokes wave.
 * Pressure follows from the steady Bernoulli equation in the frame moving with the wave.
 * Reference: Fenton (1985) Eqs. 12–16.
 * Validity: −d ≤ z ≤ eta; above the surface zeros are returned with wet = false.
 * @param {object} sol result of `stokes5Wave`
 * @param {object} p
 * @param {number} [p.x=0] m
 * @param {number} [p.z=0] m above still-water level
 * @param {number} [p.t=0] s
 * @param {number} [p.rho=1025] kg/m³
 * @returns {{eta:number,u:number,w:number,ax:number,az:number,axTotal:number,azTotal:number,pDyn:number,wet:boolean}}
 *   eta m; u, w m/s (u includes the Eulerian current); ax, az local accelerations m/s²; axTotal, azTotal particle
 *   (material) accelerations m/s²; pDyn = p + ρ g z, Pa.
 */
export function stokes5Kinematics(sol, { x = 0, z = 0, t = 0, rho = RHO_SEA } = {}) {
  const theta = sol.k * (x - sol.c * t);
  let eta = 0;
  for (let j = 1; j <= 5; j++) eta += sol.etaAmp[j] * Math.cos(j * theta);
  if (Number.isFinite(sol.d) && z < -sol.d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  if (z > eta + 1e-9 * sol.H) return { eta, u: 0, w: 0, ax: 0, az: 0, axTotal: 0, azTotal: 0, pDyn: 0, wet: false };
  const r = stokes5Raw(sol, theta, z);
  const ur = r.u - sol.c;
  // irrotational, incompressible: w_x = u_z, w_z = −u_x
  return {
    eta, u: r.u, w: r.w,
    ax: -sol.c * r.ux, az: -sol.c * r.uz,
    axTotal: ur * r.ux + r.w * r.uz, azTotal: ur * r.uz - r.w * r.ux,
    pDyn: rho * (sol.R - 0.5 * (ur * ur + r.w * r.w)),
    wet: true,
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Elliptic integrals and Jacobi elliptic functions                                                              */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Complete elliptic integrals of the first and second kind, K(m) and E(m), by the arithmetic–geometric mean.
 * Reference: Abramowitz & Stegun (1964) §17.6.
 * Validity: parameter 0 ≤ m < 1 (m = k², not the modulus). For m extremely close to 1 pass the complementary
 * parameter m1 = 1 − m explicitly to avoid cancellation.
 * @param {number} m parameter, -
 * @param {number} [m1=1-m] complementary parameter, -
 * @returns {{K:number,E:number}} dimensionless
 */
export function ellipticKE(m, m1 = 1 - m) {
  if (!(m >= 0) || !(m1 > 0) || m > 1) throw new RangeError('m must satisfy 0 <= m < 1');
  let a = 1; let b = Math.sqrt(m1); let c2 = m; let sum = 0.5 * c2; let p = 0.5;
  for (let i = 0; i < 60; i++) {
    const cn = 0.5 * (a - b);
    if (Math.abs(cn) < 1e-17 * a) break;
    const an = 0.5 * (a + b);
    b = Math.sqrt(a * b); a = an; p *= 2; sum += p * cn * cn;
  }
  const K = Math.PI / (2 * a);
  return { K, E: K * (1 - sum) };
}

/**
 * Jacobi elliptic functions sn(u|m), cn(u|m), dn(u|m) by descending Landen (AGM) transformation.
 * Reference: Abramowitz & Stegun (1964) §16.4.
 * Validity: real u, 0 ≤ m ≤ 1 (m = 1 returns tanh, sech, sech). Pass m1 = 1 − m when m is very close to 1.
 * @param {number} u argument, -
 * @param {number} m parameter, -
 * @param {number} [m1=1-m] complementary parameter, -
 * @returns {{sn:number,cn:number,dn:number}} dimensionless
 */
export function jacobiElliptic(u, m, m1 = 1 - m) {
  if (!(m >= 0) || m > 1 || !(m1 >= 0)) throw new RangeError('m must satisfy 0 <= m <= 1');
  if (m1 === 0) { const s = 1 / Math.cosh(u); return { sn: Math.tanh(u), cn: s, dn: s }; }
  if (m === 0) return { sn: Math.sin(u), cn: Math.cos(u), dn: 1 };
  const a = [1]; const c = [Math.sqrt(m)];
  let b = Math.sqrt(m1);
  let n = 0;
  while (Math.abs(c[n]) > 1e-16 * a[n] && n < 60) {
    const an = 0.5 * (a[n] + b);
    c.push(0.5 * (a[n] - b));
    b = Math.sqrt(a[n] * b);
    a.push(an);
    n++;
  }
  let phi = 2 ** n * a[n] * u;
  for (let i = n; i > 0; i--) phi = 0.5 * (phi + Math.asin((c[i] / a[i]) * Math.sin(phi)));
  const sn = Math.sin(phi);
  return { sn, cn: Math.cos(phi), dn: Math.sqrt(1 - m * sn * sn) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Shallow-water waves of permanent form: cnoidal and solitary (first order)                                     */
/* ------------------------------------------------------------------------------------------------------------ */

// Kinematics of a long wave of permanent form travelling at speed c with zero net mass flux:
// continuity in the moving frame gives the depth-uniform velocity u = c η/(d + η); w follows from continuity.
function longWaveKinematics(eta, ex, exx, c, d, z, g, rho) {
  const h = d + eta;
  const u = (c * eta) / h;
  const ux = (c * d * ex) / (h * h);
  const uxx = c * d * (exx / (h * h) - (2 * ex * ex) / (h * h * h));
  const s = z + d;
  const w = -s * ux;
  return {
    eta, u, w, ax: -c * ux, az: c * s * uxx,
    // hydrostatic under the instantaneous surface (first-order long-wave theory)
    pDyn: rho * g * eta,
  };
}

/**
 * First-order solitary wave (Boussinesq): eta = H sech²[ sqrt(3H/(4d³)) (x − c t) ], c = sqrt(g (d + H)).
 * Kinematics are depth-uniform horizontal velocity from mass conservation, u = c η/(d + η), with w from continuity
 * and hydrostatic pressure.
 * Reference: Boussinesq (1872); Munk, W.H. (1949) "The solitary wave theory and its application to surf problems";
 * Dean & Dalrymple (1991) §11.5.
 * Validity: H/d below about 0.3 for the first-order profile; the highest solitary wave has H/d ≈ 0.78–0.83.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.d depth, m
 * @param {number} [o.x=0] m (crest at x = c t)
 * @param {number} [o.z=0] m above still-water level
 * @param {number} [o.t=0] s
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{eta:number,u:number,w:number,ax:number,az:number,pDyn:number,c:number,beta:number,volume:number,
 *   energy:number,wet:boolean}} eta m; u, w m/s; ax, az m/s² (local); pDyn Pa; c m/s; beta decay wavenumber 1/m;
 *   volume above still water per unit crest width m²; energy per unit crest width J/m.
 */
export function solitaryWave({ H, d, x = 0, z = 0, t = 0, rho = RHO_SEA, g = G }) {
  reqFinitePos('H', H); reqFinitePos('d', d);
  if (z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  const c = Math.sqrt(g * (d + H));
  const beta = Math.sqrt((3 * H) / (4 * d * d * d));
  const X = beta * (x - c * t);
  const sech = 1 / Math.cosh(X); const th = Math.tanh(X);
  const eta = H * sech * sech;
  const ex = -2 * H * beta * sech * sech * th;
  const exx = 2 * H * beta * beta * sech * sech * (3 * th * th - 1);
  const kin = longWaveKinematics(eta, ex, exx, c, d, Math.min(z, eta), g, rho);
  const wet = z <= eta + 1e-12 * H;
  const out = wet ? kin : { eta, u: 0, w: 0, ax: 0, az: 0, pDyn: 0 };
  return {
    ...out, c, beta,
    volume: 4 * d * Math.sqrt((d * H) / 3),
    energy: (8 / (3 * Math.sqrt(3))) * rho * g * Math.pow(H * d, 1.5),
    wet,
  };
}

// First-order cnoidal relations for mean depth d, height H and complementary parameter m1 = 1 − m.
function cnoidalFromM1(H, d, m1, g) {
  const m = 1 - m1;
  const { K, E } = ellipticKE(m, m1);
  const L = d * Math.sqrt((16 * m * d) / (3 * H)) * K;
  const c = Math.sqrt(g * d) * (1 + (H / (m * d)) * (1 - 0.5 * m - (1.5 * E) / K));
  return { m, m1, K, E, L, c, T: c > 0 ? L / c : Infinity, trough: (H / m) * (1 - m - E / K) };
}

/**
 * Solve the first-order (Korteweg–de Vries) cnoidal wave for the elliptic parameter m given H, d and either the
 * period or the wavelength:
 *   L = d sqrt(16 m d / (3 H)) K(m),   c = sqrt(g d) [1 + (H/(m d)) (1 − m/2 − 3E/(2K))],
 *   trough elevation = (H/m)(1 − m − E/K),   eta = trough + H cn²(2K (x − c t)/L | m).
 * Reference: Korteweg & de Vries (1895); Wiegel, R.L. (1960) J. Fluid Mech. 7, 273–286;
 * Dingemans (1997) Water Wave Propagation over Uneven Bottoms §6.2.
 * Validity: long waves, d/L below about 1/8 and Ursell number above about 26; H/d below about 0.6.
 * Throws RangeError when the period is too short for a cnoidal solution to exist.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.d mean depth, m
 * @param {number} [o.T] period, s (give T or L)
 * @param {number} [o.L] wavelength, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{H:number,d:number,g:number,m:number,m1:number,K:number,E:number,L:number,T:number,c:number,
 *   crest:number,trough:number,ursell:number}} m, m1 = 1 − m (-); K, E complete elliptic integrals (-); L m; T s;
 *   c m/s; crest and trough elevations relative to the mean level, m.
 */
export function cnoidalWave({ H, d, T, L, g = G }) {
  reqFinitePos('H', H); reqFinitePos('d', d);
  // unknown expressed as q = −ln(1 − m) so that the solitary limit (m -> 1) is resolved
  let fq;
  if (L !== undefined && L !== null) {
    reqFinitePos('L', L);
    fq = (q) => cnoidalFromM1(H, d, Math.exp(-q), g).L - L;
  } else {
    reqFinitePos('T', T);
    fq = (q) => cnoidalFromM1(H, d, Math.exp(-q), g).T - T;
  }
  // Scan downward from the solitary limit; the physical branch is the one with the larger m.
  const n = 600;
  const qMax = 690; const qMin = 1e-6;
  let qb = qMax;
  if (!(fq(qb) > 0)) throw new RangeError('T (or L) is too long to resolve with the cnoidal parameter in double precision');
  let qa = NaN;
  for (let i = 1; i <= n; i++) {
    const q = qMax * Math.pow(qMin / qMax, i / n);
    const v = fq(q);
    if (v <= 0) { qa = q; break; }
    if (!Number.isFinite(v)) break;
    qb = q;
  }
  if (Number.isNaN(qa)) throw new RangeError('T (or L) is too short for cnoidal theory at this H and d');
  for (let i = 0; i < 200; i++) {
    const qm = 0.5 * (qa + qb);
    if (fq(qm) > 0) qb = qm; else qa = qm;
    if (qb - qa <= 1e-15 * qb) break;
  }
  const r = cnoidalFromM1(H, d, Math.exp(-0.5 * (qa + qb)), g);
  return { H, d, g, m: r.m, m1: r.m1, K: r.K, E: r.E, L: r.L, T: r.T, c: r.c, crest: r.trough + H, trough: r.trough,
    ursell: (H * r.L * r.L) / (d * d * d) };
}

/**
 * Surface elevation and kinematics of a solved first-order cnoidal wave. Horizontal velocity is depth-uniform,
 * u = c η/(d + η) (zero net mass flux); w from continuity; pressure hydrostatic under the instantaneous surface.
 * Reference: as `cnoidalWave`.
 * Validity: as `cnoidalWave`; −d ≤ z ≤ eta.
 * @param {object} sol result of `cnoidalWave`
 * @param {object} p
 * @param {number} [p.x=0] m (crest at x = c t)
 * @param {number} [p.z=0] m above mean water level
 * @param {number} [p.t=0] s
 * @param {number} [p.rho=1025] kg/m³
 * @returns {{eta:number,u:number,w:number,ax:number,az:number,pDyn:number,wet:boolean}} eta m; u, w m/s;
 *   ax, az m/s² (local); pDyn Pa.
 */
export function cnoidalKinematics(sol, { x = 0, z = 0, t = 0, rho = RHO_SEA } = {}) {
  const { H, d, m, m1, K, L, c, g } = sol;
  if (z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  const beta = (2 * K) / L;
  // reduce the phase to one wavelength about the crest to keep the elliptic argument small
  let s = (x - c * t) / L;
  s -= Math.round(s);
  const j = jacobiElliptic(beta * L * s, m, m1);
  const cn2 = j.cn * j.cn;
  const eta = sol.trough + H * cn2;
  if (z > eta + 1e-12 * H) return { eta, u: 0, w: 0, ax: 0, az: 0, pDyn: 0, wet: false };
  const ex = -2 * H * beta * j.cn * j.sn * j.dn;
  const sn2 = j.sn * j.sn; const dn2 = j.dn * j.dn;
  const exx = 2 * H * beta * beta * (sn2 * dn2 - cn2 * dn2 + m * sn2 * cn2);
  return { ...longWaveKinematics(eta, ex, exx, c, d, z, g, rho), wet: true };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Breaking limits and theory selection                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Miche steepness-limited breaking height, Hb = 0.142 L tanh(2π d / L).
 * Reference: Miche (1944); DNV-RP-C205 (2019) §3.4.6.1.
 * Validity: regular waves over a horizontal bed.
 * @param {number} L wavelength at depth d, m
 * @param {number} d depth, m (Infinity allowed)
 * @returns {number} breaking wave height, m
 */
export function breakingHeightMiche(L, d) {
  reqFinitePos('L', L); reqPos('d', d);
  return 0.142 * L * (Number.isFinite(d) ? Math.tanh((TWO_PI * d) / L) : 1);
}

/**
 * Depth-limited breaking height Hb = κ d with McCowan's κ = 0.78.
 * Reference: McCowan (1894); Shore Protection Manual (1984) Eq. 2-91; DNV-RP-C205 (2019) §3.4.6.2.
 * Validity: shallow water over a near-horizontal bed; κ varies from about 0.6 to 1.3 with beach slope.
 * @param {number} d depth, m
 * @param {number} [kappa=0.78] breaker index, -
 * @returns {number} breaking wave height, m
 */
export function breakingHeightMcCowan(d, kappa = 0.78) { reqFinitePos('d', d); reqFinitePos('kappa', kappa); return kappa * d; }

/**
 * Goda's breaker-height formula for regular waves on a slope,
 * Hb = A L0 { 1 − exp[ −1.5 (π d / L0) (1 + K tan^{4/3} β) ] },  L0 = g T²/(2π).
 * Reference: Goda, Y. (1974) Proc. 14th ICCE (A = 0.17, K = 15); Goda (2010) Coastal Eng. 57 revised K to 11.
 * Validity: bed slopes 1/100 to 1/10; A between 0.12 (lower bound for irregular waves) and 0.18.
 * @param {object} o
 * @param {number} o.T period, s
 * @param {number} o.d depth, m
 * @param {number} [o.slope=0] bed slope tan β, -
 * @param {number} [o.A=0.17] coefficient, -
 * @param {number} [o.K=15] slope coefficient, -
 * @param {number} [o.g=9.80665] m/s²
 * @returns {number} breaking wave height, m
 */
export function breakingHeightGoda({ T, d, slope = 0, A = 0.17, K = 15, g = G }) {
  reqFinitePos('T', T); reqFinitePos('d', d);
  if (!(slope >= 0)) throw new RangeError('slope must be non-negative');
  const L0 = (g * T * T) / TWO_PI;
  return A * L0 * (1 - Math.exp(-1.5 * ((Math.PI * d) / L0) * (1 + K * Math.pow(slope, 4 / 3))));
}

/**
 * Governing breaking height for a regular wave: the smaller of the Miche and McCowan limits, together with Goda.
 * Reference: as the individual functions.
 * @param {object} o
 * @param {number} o.T period, s
 * @param {number} o.d depth, m
 * @param {number} [o.slope=0] bed slope, -
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{miche:number,mccowan:number,goda:number,governing:number,L:number}} heights m; L linear wavelength m.
 */
export function breakingLimits({ T, d, slope = 0, g = G }) {
  const L = wavelength(T, d, g);
  const miche = breakingHeightMiche(L, d);
  const mccowan = breakingHeightMcCowan(d);
  return { miche, mccowan, goda: breakingHeightGoda({ T, d, slope, g }), governing: Math.min(miche, mccowan), L };
}

/**
 * Recommend a regular-wave theory from the dimensionless height H/(g T²) and depth d/(g T²), following the
 * layout of the applicability chart in API RP 2A-WSD (2014) Fig. 5.3 / ISO 19901-1 (2015) Fig. A.4
 * (after Atkins 1990 and Le Méhauté 1976).
 * The region boundaries of that chart are curves given only graphically; they are APPROXIMATED here by:
 *   breaking           H > min(0.142 L tanh kd, 0.78 d)            (Miche / McCowan with the linear wavelength);
 *   linear (Airy)      H/L < 0.006 (≈ H/gT² < 0.001 in deep water) and Ursell number < 5;
 *   cnoidal            Ursell number > 40 and d/L < 0.125;
 *   Stokes 5th order   Ursell number ≤ 40 and H ≤ 0.9 Hb;
 *   stream function    everything else below breaking (always admissible; order rises towards breaking:
 *                      3 below 0.25 Hb, 5 below 0.5 Hb, 7 below 0.75 Hb, 9 below 0.9 Hb, else 11).
 * Use the result as a screening recommendation, not as a substitute for the chart.
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.T apparent period, s
 * @param {number} o.d depth, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{theory:'airy'|'stokes5'|'cnoidal'|'streamFunction'|'breaking',streamFunctionOrder:number,
 *   HgT2:number,dgT2:number,ursell:number,steepness:number,dOverL:number,HoverHb:number,regime:string,
 *   approximate:true}} dimensionless ratios; streamFunctionOrder is the suggested Fourier order.
 */
export function selectWaveTheory({ H, T, d, g = G }) {
  reqFinitePos('H', H); reqFinitePos('T', T); reqFinitePos('d', d);
  const L = wavelength(T, d, g);
  const Hb = Math.min(breakingHeightMiche(L, d), 0.78 * d);
  const ur = (H * L * L) / (d * d * d);
  const st = H / L;
  const r = H / Hb;
  let theory;
  if (r > 1) theory = 'breaking';
  else if (st < 0.006 && ur < 5) theory = 'airy';
  else if (ur > 40 && d / L < 0.125) theory = r <= 0.75 ? 'cnoidal' : 'streamFunction';
  else if (ur <= 40 && r <= 0.9) theory = 'stokes5';
  else theory = 'streamFunction';
  const order = r < 0.25 ? 3 : r < 0.5 ? 5 : r < 0.75 ? 7 : r < 0.9 ? 9 : 11;
  return { theory, streamFunctionOrder: order, HgT2: H / (g * T * T), dgT2: d / (g * T * T), ursell: ur,
    steepness: st, dOverL: d / L, HoverHb: r, regime: depthRegime(d / L), approximate: true };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Shoaling and refraction                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Linear shoaling and refraction over straight parallel depth contours (Snell's law, energy-flux conservation):
 *   sin θ / c = sin θ0 / c0,   Ks = sqrt(cg0 / cg),   Kr = sqrt(cos θ0 / cos θ),   H = H0 Ks Kr.
 * Reference: Dean & Dalrymple (1991) §4.8; Shore Protection Manual (1984) §2.3.
 * Validity: linear waves, no breaking, dissipation, current or diffraction; |θ0| < 90°.
 * @param {object} o
 * @param {number} o.T period, s
 * @param {number} o.d local depth, m
 * @param {number} [o.d0=Infinity] depth at which H0 and theta0 are specified, m (default deep water)
 * @param {number} [o.theta0=0] angle between the wave direction and the contour normal at d0, rad
 * @param {number} [o.H0=1] wave height at d0, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{Ks:number,Kr:number,theta:number,H:number,k:number,L:number,c:number,cg:number}} Ks, Kr (-);
 *   theta rad; H m; k rad/m; L m; c, cg m/s at depth d.
 */
export function shoalingRefraction({ T, d, d0 = Infinity, theta0 = 0, H0 = 1, g = G }) {
  if (!(Math.abs(theta0) < Math.PI / 2)) throw new RangeError('theta0 must lie within (-π/2, π/2)');
  const k0 = wavenumber({ T, d: d0, g });
  const k = wavenumber({ T, d, g });
  const cg0 = sigmaOf(k0, d0, g).cg;
  const cg = sigmaOf(k, d, g).cg;
  const sinT = (Math.sin(theta0) * k0) / k;
  if (Math.abs(sinT) >= 1) throw new RangeError('theta0: wave is totally reflected (caustic) before reaching depth d');
  const theta = Math.asin(sinT);
  const Ks = Math.sqrt(cg0 / cg);
  const Kr = Math.sqrt(Math.cos(theta0) / Math.cos(theta));
  return { Ks, Kr, theta, H: H0 * Ks * Kr, k, L: TWO_PI / k, c: TWO_PI / (T * k), cg };
}
