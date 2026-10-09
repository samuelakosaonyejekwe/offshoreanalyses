// Current profiles, wind-driven current, stretching of the profile under waves, depth averaging and combined
// wave + current near-bed velocity.
//
// Conventions: z is elevation above still-water level (m, negative below the surface, −d at the seabed) unless an
// argument is named `zBed` (height above the seabed, m). Speeds in m/s.

import { wavenumber, G } from './waves.js';

function reqPos(name, v) {
  if (typeof v !== 'number' || !(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a positive finite number`);
}
function reqFinite(name, v) {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new RangeError(`${name} must be a finite number`);
}
function inColumn(z, d) {
  reqFinite('z', z);
  if (z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
}

/**
 * Seabed roughness parameter z0 by soil type, m.
 * Reference: DNV-RP-F109 (2021) Table 3-1.
 */
export const SEABED_ROUGHNESS = Object.freeze({
  silt: 5e-6, fineSand: 1e-5, mediumSand: 4e-5, coarseSand: 1e-4, gravel: 3e-4, pebble: 2e-3, cobble: 1e-2, boulder: 4e-2,
});

/**
 * Power-law (tidal) current profile, v(z) = v0 [(d + z)/d]^α for −d ≤ z ≤ 0, with α = 1/7 by default.
 * Above the still-water level the surface value is returned (apply stretching for wave crests).
 * Reference: DNV-RP-C205 (2019) §4.1.4.2; ISO 19901-1 (2015) §A.9.3.
 * Validity: tidal current in shallow, well-mixed water without stratification.
 * @param {number} z elevation above still-water level, m (≤ 0)
 * @param {object} o
 * @param {number} o.v0 current speed at the still-water level, m/s
 * @param {number} o.d water depth, m
 * @param {number} [o.alpha=1/7] exponent, -
 * @returns {number} current speed at z, m/s
 */
export function currentPowerLaw(z, { v0, d, alpha = 1 / 7 }) {
  reqFinite('v0', v0); reqPos('d', d); inColumn(z, d);
  const s = Math.min(1, Math.max(0, (d + z) / d));
  return v0 * Math.pow(s, alpha);
}

/**
 * Linear (wind-generated) current profile, v(z) = v0 (d0 + z)/d0 for −d0 ≤ z ≤ 0 and zero below, where d0 is the
 * reference depth of the wind-generated current (50 m by default, limited to the water depth).
 * Reference: DNV-RP-C205 (2019) §4.1.4.2.
 * Validity: local wind-driven current in open water.
 * @param {number} z elevation above still-water level, m
 * @param {object} o
 * @param {number} o.v0 wind-generated current speed at the still-water level, m/s
 * @param {number} [o.d0=50] reference depth, m
 * @param {number} [o.d=Infinity] water depth, m
 * @returns {number} current speed at z, m/s
 */
export function currentLinear(z, { v0, d0 = 50, d = Infinity }) {
  reqFinite('v0', v0); reqPos('d0', d0); reqFinite('z', z);
  if (!(d > 0)) throw new RangeError('d must be positive');
  const dd = Math.min(d0, d);
  if (z >= 0) return v0;
  if (z <= -dd) return 0;
  return (v0 * (dd + z)) / dd;
}

/**
 * Bilinear current profile through three points: surface (z = 0), a break point (z = zBreak) and the seabed.
 * Reference: piecewise-linear design profiles as commonly specified in metocean criteria (ISO 19901-1 §A.9.3).
 * @param {number} z elevation above still-water level, m
 * @param {object} o
 * @param {number} o.vSurface speed at z = 0, m/s
 * @param {number} o.vBreak speed at the break point, m/s
 * @param {number} o.vBed speed at the seabed, m/s
 * @param {number} o.zBreak elevation of the break point, m (−d < zBreak < 0)
 * @param {number} o.d water depth, m
 * @returns {number} current speed at z, m/s
 */
export function currentBilinear(z, { vSurface, vBreak, vBed, zBreak, d }) {
  reqPos('d', d); inColumn(z, d);
  if (!(zBreak < 0 && zBreak > -d)) throw new RangeError('zBreak must lie between -d and 0');
  if (z >= 0) return vSurface;
  if (z >= zBreak) return vBreak + ((vSurface - vBreak) * (z - zBreak)) / -zBreak;
  return vBed + ((vBreak - vBed) * (z + d)) / (zBreak + d);
}

/**
 * Logarithmic near-bed current profile, V(zBed) = V(zr) [ln(zBed + z0) − ln z0] / [ln(zr + z0) − ln z0].
 * Reference: DNV-RP-F109 (2021) §3.4.2, Eq. 3.2 (z0 from `SEABED_ROUGHNESS`).
 * Validity: turbulent boundary layer over a flat seabed without stratification; zBed ≥ 0.
 * @param {number} zBed height above the seabed, m
 * @param {object} o
 * @param {number} o.vRef current speed at the reference height, m/s
 * @param {number} o.zRef reference height above the seabed, m
 * @param {number} o.z0 bottom roughness parameter, m
 * @returns {number} current speed at zBed, m/s
 */
export function currentLogNearBed(zBed, { vRef, zRef, z0 }) {
  reqFinite('vRef', vRef); reqPos('zRef', zRef); reqPos('z0', z0);
  if (!(zBed >= 0)) throw new RangeError('zBed must be non-negative');
  return (vRef * Math.log(zBed / z0 + 1)) / Math.log(zRef / z0 + 1);
}

/**
 * Mean perpendicular current over the height of a pipe (or other object of height D) resting on the seabed in a
 * logarithmic boundary layer:
 *   Vc = V(zr) { [(1 + z0/D) ln(D/z0 + 1) − 1] / ln(zr/z0 + 1) } sin θ.
 * Reference: DNV-RP-F109 (2021) §3.4.2, Eq. 3.3 (the bracket is the exact average of the log profile over 0..D).
 * @param {object} o
 * @param {number} o.vRef current speed at the reference height, m/s
 * @param {number} o.zRef reference height above the seabed, m
 * @param {number} o.z0 bottom roughness parameter, m
 * @param {number} o.D object height (pipe outer diameter), m
 * @param {number} [o.theta=π/2] angle between current direction and the pipe axis, rad
 * @returns {number} mean perpendicular current over the object height, m/s
 */
export function meanCurrentOverHeight({ vRef, zRef, z0, D, theta = Math.PI / 2 }) {
  reqFinite('vRef', vRef); reqPos('zRef', zRef); reqPos('z0', z0); reqPos('D', D);
  return (vRef * ((1 + z0 / D) * Math.log(D / z0 + 1) - 1) * Math.sin(theta)) / Math.log(zRef / z0 + 1);
}

/**
 * Tabulated current profile with linear interpolation between points; values are held constant outside the table.
 * Reference: site-specific measured or hindcast profiles (ISO 19901-1 (2015) §9.3).
 * @param {number} z elevation above still-water level, m
 * @param {object} o
 * @param {ArrayLike<number>} o.zs elevations in ascending order (e.g. −d … 0), m
 * @param {ArrayLike<number>} o.vs speeds at zs, m/s
 * @returns {number} current speed at z, m/s
 */
export function currentTabulated(z, { zs, vs }) {
  const n = zs.length;
  if (n < 1 || vs.length !== n) throw new RangeError('zs and vs must have the same non-zero length');
  reqFinite('z', z);
  if (z <= zs[0]) return vs[0];
  if (z >= zs[n - 1]) return vs[n - 1];
  let lo = 0; let hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (zs[m] <= z) lo = m; else hi = m; }
  const dz = zs[hi] - zs[lo];
  if (!(dz > 0)) throw new RangeError('zs must be strictly ascending');
  return vs[lo] + ((vs[hi] - vs[lo]) * (z - zs[lo])) / dz;
}

/**
 * Wind-generated surface current by the percentage rule, v = factor × U(1 hour, 10 m).
 * Reference: DNV-RP-C205 (2019) §4.1.4.4 gives factor = 0.015 to 0.03 for deep water without site data; the
 * default here is the upper bound 0.03.
 * Validity: open water, steady wind, no coastal or density-driven circulation.
 * @param {number} U1h10 1-hour mean wind speed at 10 m, m/s
 * @param {number} [factor=0.03] ratio of surface current to wind speed (0.015–0.03), -
 * @returns {number} wind-generated current at the still-water level, m/s
 */
export function windDrivenSurfaceCurrent(U1h10, factor = 0.03) {
  if (!(U1h10 >= 0)) throw new RangeError('U1h10 must be non-negative');
  if (!(factor >= 0 && factor <= 0.1)) throw new RangeError('factor must lie in [0, 0.1]');
  return factor * U1h10;
}

/**
 * Total current as the sum of a tidal power-law component and a wind-generated linear component.
 * Reference: DNV-RP-C205 (2019) §4.1.4.1–4.1.4.2.
 * @param {number} z elevation above still-water level, m
 * @param {object} o
 * @param {number} o.vTide tidal current at the still-water level, m/s
 * @param {number} o.vWind wind-generated current at the still-water level, m/s
 * @param {number} o.d water depth, m
 * @param {number} [o.alpha=1/7] tidal exponent, -
 * @param {number} [o.d0=50] reference depth of the wind-generated current, m
 * @returns {number} current speed at z, m/s
 */
export function currentTideAndWind(z, { vTide, vWind, d, alpha = 1 / 7, d0 = 50 }) {
  return currentPowerLaw(Math.min(z, 0), { v0: vTide, d, alpha }) + currentLinear(z, { v0: vWind, d0, d });
}

/**
 * Elevation z' in the still-water current profile that corresponds to elevation z under a wave of local surface
 * elevation η (current stretching):
 *   'linear'     z' = (z − η) d/(d + η)                              (stretches the whole column uniformly);
 *   'nonlinear'  z = z' + η sinh[k(z' + d)]/sinh(k d), solved for z'  (stretching weighted towards the surface;
 *                k is the nonlinear wavenumber of the wave);
 *   'none'       z' = min(z, 0)                                       (profile unchanged, extended uniformly
 *                above the still-water level).
 * The stretched current at z is the still-water profile evaluated at z'.
 * Reference: API RP 2A-WSD (2014) §5.3.1.2 item 5 (linear and nonlinear stretching); ISO 19901-1 (2015) §A.9.4;
 * DNV-RP-C205 (2019) §4.1.5.
 * Validity: −d ≤ z ≤ η. Linear stretching is an acceptable approximation for slab or gently sheared profiles.
 * @param {object} o
 * @param {number} o.z elevation above still-water level, m
 * @param {number} o.eta instantaneous surface elevation at the location, m
 * @param {number} o.d still-water depth, m
 * @param {'linear'|'nonlinear'|'none'} [o.method='linear']
 * @param {number} [o.k] wavenumber, rad/m (required for 'nonlinear')
 * @returns {number} z', m (between −d and 0)
 */
export function stretchedElevation({ z, eta, d, method = 'linear', k }) {
  reqFinite('z', z); reqFinite('eta', eta); reqPos('d', d);
  if (!(d + eta > 0)) throw new RangeError('eta must be above the seabed');
  if (z < -d * (1 + 1e-12)) throw new RangeError('z is below the seabed');
  if (z > eta + 1e-9 * d) throw new RangeError('z is above the instantaneous free surface');
  if (method === 'none') return Math.min(z, 0);
  if (method === 'linear') return ((z - eta) * d) / (d + eta);
  if (method !== 'nonlinear') throw new RangeError("method must be 'linear', 'nonlinear' or 'none'");
  reqPos('k', k);
  const kd = k * d;
  // sinh k(z'+d)/sinh kd in decaying exponentials; f(z') = z' + η·ratio − z is monotonic for η > −d
  const ratio = (zp) => (Math.exp(k * zp) * (1 - Math.exp(-2 * k * (zp + d)))) / (1 - Math.exp(-2 * kd));
  let lo = -d; let hi = 0;
  let zp = ((z - eta) * d) / (d + eta);
  for (let i = 0; i < 100; i++) {
    const f = zp + eta * ratio(zp) - z;
    if (f > 0) hi = zp; else lo = zp;
    const cosr = (Math.exp(k * zp) * (1 + Math.exp(-2 * k * (zp + d)))) / (1 - Math.exp(-2 * kd));
    const df = 1 + eta * k * cosr;
    let zn = df > 0 ? zp - f / df : NaN;
    if (!(zn > lo && zn < hi)) zn = 0.5 * (lo + hi);
    if (Math.abs(zn - zp) < 1e-13 * d) { zp = zn; break; }
    zp = zn;
  }
  return Math.min(0, Math.max(-d, zp));
}

/**
 * Current speed at elevation z under a wave, from a still-water profile and a stretching method.
 * Reference and validity: see `stretchedElevation`.
 * @param {(z:number)=>number} profile still-water current profile v(z), m/s for −d ≤ z ≤ 0
 * @param {object} o arguments of `stretchedElevation`
 * @returns {number} current speed, m/s
 */
export function stretchedCurrent(profile, o) { return profile(stretchedElevation(o)); }

/**
 * Depth-averaged current, v̄ = (1/d) ∫_{−d}^{0} v(z) dz, by composite Simpson integration in a coordinate
 * clustered towards the seabed (handles the infinite slope of power-law profiles there).
 * Reference: definition; for the 1/7 power law v̄ = v0/(1 + α) = 0.875 v0.
 * @param {(z:number)=>number} profile current profile v(z), m/s
 * @param {object} o
 * @param {number} o.d water depth, m
 * @param {number} [o.n=400] number of intervals
 * @returns {number} depth-averaged speed, m/s
 */
export function depthAveragedCurrent(profile, { d, n = 400 }) {
  reqPos('d', d);
  const N = 2 * Math.ceil(n / 2);
  // substitution (z + d)/d = s⁴ removes the bed singularity of s^α-type profiles
  let sum = 0;
  for (let i = 0; i <= N; i++) {
    const s = i / N;
    const w = i === 0 || i === N ? 1 : i % 2 ? 4 : 2;
    sum += w * profile(-d + d * s ** 4) * 4 * s ** 3;
  }
  return sum / (3 * N);
}

/**
 * Combined wave and current velocity near the seabed by linear theory. The wave orbital velocity amplitude at
 * height zBed is Uw = (π H/T_r) cosh(k zBed)/sinh(k d), with k from the Doppler-shifted dispersion relation for
 * the current component along the wave direction (T_r = intrinsic period 2π/(ω − k Uc cos φ)).
 * The extreme combined speed is the vector sum of the steady current and the oscillatory velocity.
 * Reference: DNV-RP-C205 (2019) §3.2.2 and §4.1.5; DNV-RP-F109 (2021) §3.4 (regular-wave form); Soulsby (1997)
 * Dynamics of Marine Sands §4–5 for the definitions.
 * Validity: linear wave, uniform current over the wave-affected depth, no wave–current boundary-layer interaction
 * (bed shear stress enhancement is not modelled).
 * @param {object} o
 * @param {number} o.H wave height, m
 * @param {number} o.T wave period (fixed frame), s
 * @param {number} o.d water depth, m
 * @param {number} [o.Uc=0] steady current speed near the bed, m/s
 * @param {number} [o.angle=0] angle between current and wave propagation directions, rad
 * @param {number} [o.zBed=0] height above the seabed, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{Uw:number,Aw:number,k:number,Tr:number,Umax:number,Umin:number,ratio:number}} Uw orbital velocity
 *   amplitude m/s; Aw orbital semi-excursion m; k rad/m; Tr intrinsic period s; Umax largest and Umin smallest
 *   combined speed magnitude during a cycle, m/s; ratio = Uc/(Uc + Uw) (-)
 */
export function nearBedWaveCurrent({ H, T, d, Uc = 0, angle = 0, zBed = 0, g = G }) {
  if (!(H >= 0)) throw new RangeError('H must be non-negative');
  reqPos('d', d);
  if (!(zBed >= 0 && zBed <= d)) throw new RangeError('zBed must lie between 0 and d');
  const Ualong = Uc * Math.cos(angle);
  const k = wavenumber({ T, d, U: Ualong, g });
  const sig = (2 * Math.PI) / T - k * Ualong;
  const kd = k * d;
  const cs = (Math.exp(k * (zBed - d)) * (1 + Math.exp(-2 * k * zBed))) / (1 - Math.exp(-2 * kd));
  const Uw = 0.5 * H * sig * cs;
  const c = Math.cos(angle);
  const sp = Math.sqrt(Uc * Uc + Uw * Uw + 2 * Uc * Uw * Math.abs(c));
  // the smallest magnitude along the oscillation u = Uc + s·Uw (s in [−1, 1]) is the perpendicular offset if the
  // projection reverses, otherwise the end point
  const sMin = Uw > 0 ? Math.max(-1, Math.min(1, (-Uc * c) / Uw)) : 0;
  const mn = Math.sqrt(Math.max(0, Uc * Uc + sMin * sMin * Uw * Uw + 2 * Uc * Uw * sMin * c));
  return { Uw, Aw: sig > 0 ? Uw / sig : NaN, k, Tr: (2 * Math.PI) / sig, Umax: sp, Umin: mn, ratio: Uc + Uw > 0 ? Math.abs(Uc) / (Math.abs(Uc) + Uw) : NaN };
}
