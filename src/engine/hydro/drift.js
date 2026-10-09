// Mean and slowly varying wave drift, low-frequency response of a moored floater, and wind / current loads
// from projected areas and coefficients.
// Heading/direction convention: direction the environment travels towards, measured from the vessel +x (bow)
// axis towards +y (port), deg. Forces are in the vessel frame (Fx along +x, Fy along +y, Mz about +z).

import { G, dispersion, needPos } from './_wavehelp.js';

const DEG = Math.PI / 180;

/**
 * Mean wave drift force per unit length on a long two-dimensional body in beam regular waves:
 *   F = ½ rho g a² R² (1 + 2kh/sinh 2kh),
 * with R the amplitude reflection coefficient (deep water: ½ rho g a² R², Maruo's formula).
 * Reference: Maruo (1960) The drift of a body floating on waves, J. Ship Research 4(3); Longuet-Higgins (1977)
 * The mean forces exerted by waves on floating or submerged bodies, Proc. R. Soc. A 352 (finite depth,
 * F = ¼ rho g a² (1 + R² - T²)(1 + 2kh/sinh 2kh) with R² + T² = 1).
 * Valid for 2-D flow without energy dissipation. R is an input (1 = full reflection, the short-wave limit of
 * a wall-sided body).
 * @param {{a:number, R:number, T?:number, depth?:number, rho?:number, g?:number}} o wave amplitude m, reflection
 *   coefficient 0..1, wave period s (needed only for finite depth), depth m, kg/m³, m/s²
 * @returns {number} mean drift force per unit length in the wave direction, N/m
 */
export function maruoDrift({ a, R, T, depth = Infinity, rho = 1025, g = G }) {
  if (!(a >= 0)) throw new RangeError('a must be >= 0');
  if (!(R >= 0 && R <= 1)) throw new RangeError('R must be within 0..1');
  let fac = 1;
  if (Number.isFinite(depth)) {
    needPos(T, 'T');
    const kh = dispersion((2 * Math.PI) / T, depth, g) * depth;
    fac = kh > 25 ? 1 : 1 + (2 * kh) / Math.sinh(2 * kh);
  }
  return 0.5 * rho * g * a * a * R * R * fac;
}

/**
 * Mean drift force on a fixed, surface-piercing vertical circular cylinder in the short-wave limit, where the
 * incident wave is fully reflected by the up-wave half of the wall-sided waterline:
 *   F = ½ rho g a² ∫ cos²(theta) cos(theta) r d theta over -pi/2..pi/2 = (2/3) rho g a² r   (deep water).
 * Reference: Faltinsen (1990) Sea Loads on Ships and Offshore Structures ch. 5 (asymptotic drift-force
 * formula for small wavelengths, applied here to a circular waterline).
 * Documented approximation: an upper bound that the exact (diffraction) drift force approaches for
 * k r > about 2–3; it overestimates strongly for long waves (k r < 1), where the drift force tends to zero.
 * @param {{a:number, D:number, rho?:number, g?:number}} o wave amplitude m, diameter m
 * @returns {{force:number, coefficient:number}} N, and the drift coefficient F/a² in N/m²
 */
export function cylinderDriftShortWave({ a, D, rho = 1025, g = G }) {
  if (!(a >= 0)) throw new RangeError('a must be >= 0');
  needPos(D, 'D');
  const coefficient = (2 / 3) * rho * g * (D / 2);
  return { force: coefficient * a * a, coefficient };
}

/**
 * Mean drift force in an irregular sea from mean drift coefficients D(omega) = F_mean / a²:
 *   F_mean = 2 ∫ S(omega) D(omega) d omega.
 * Reference: Pinkster (1975); Faltinsen (1990) ch. 5; DNV-RP-C205 (2019) §7.4.
 * @param {{omega:ArrayLike<number>, S:ArrayLike<number>, D:ArrayLike<number>}} o rad/s, m²·s/rad, N/m²
 * @returns {number} mean drift force, N
 */
export function meanDriftForce({ omega, S, D }) {
  const n = omega.length;
  if (S.length !== n || D.length !== n || n < 2) throw new RangeError('omega, S and D must have equal length >= 2');
  let s = 0;
  for (let i = 1; i < n; i++) s += 0.5 * (S[i] * D[i] + S[i - 1] * D[i - 1]) * (omega[i] - omega[i - 1]);
  return 2 * s;
}

function lin(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0; let hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + ((ys[hi] - ys[lo]) * (x - xs[lo])) / (xs[hi] - xs[lo]);
}

/**
 * Spectrum of the slowly varying (difference-frequency) drift force by Newman's approximation, in which the
 * off-diagonal quadratic transfer function is replaced by the mean drift coefficient at the mean frequency:
 *   S_F(mu) = 8 ∫₀^∞ S(omega) S(omega + mu) [D(omega + mu/2)]² d omega.
 * Reference: Newman (1974) Second-order, slowly varying forces on vessels in irregular waves; Pinkster (1975);
 * Faltinsen (1990) ch. 5; DNV-RP-C205 (2019) §7.4.
 * Valid for small difference frequencies (natural periods long compared with the wave periods), deep water.
 * @param {object} o
 * @param {ArrayLike<number>} o.omega wave frequencies, rad/s, ascending
 * @param {ArrayLike<number>} o.S wave spectrum, m²·s/rad
 * @param {ArrayLike<number>} o.D mean drift coefficient F/a², N/m²
 * @param {ArrayLike<number>} o.mu difference frequencies, rad/s
 * @returns {Float64Array} S_F(mu), N²·s/rad (one-sided)
 */
export function newmanDriftSpectrum({ omega, S, D, mu }) {
  const n = omega.length;
  if (S.length !== n || D.length !== n || n < 2) throw new RangeError('omega, S and D must have equal length >= 2');
  const out = new Float64Array(mu.length);
  const wEnd = omega[n - 1];
  for (let im = 0; im < mu.length; im++) {
    const m = mu[im];
    if (!(m >= 0)) throw new RangeError('mu must be >= 0');
    let s = 0; let prev = 0;
    for (let i = 0; i < n; i++) {
      const w = omega[i];
      const s2 = w + m <= wEnd ? lin(omega, S, w + m) : 0;
      const d = lin(omega, D, w + 0.5 * m);
      const f = S[i] * s2 * d * d;
      if (i > 0) s += 0.5 * (f + prev) * (w - omega[i - 1]);
      prev = f;
    }
    out[im] = 8 * s;
  }
  return out;
}

/**
 * Low-frequency (slow-drift) response of a moored floater as a single degree of freedom with linearised
 * damping, m x'' + c x' + k x = F_mean + F_slow(t):
 *   sigma_x² = ∫ S_F(mu) / [(k - m mu²)² + (c mu)²] d mu      (numerical, when `mu`/`SF` are given)
 *   sigma_x² ≈ pi S_F(mu_n) / (2 c k)                          (lightly damped, white-noise approximation)
 * mean offset F_mean/k and most probable extreme x_mean + sigma sqrt(2 ln N), N = duration / Tn.
 * Reference: Faltinsen (1990) ch. 5 (slow-drift motions); DNV-OS-E301 (2021) ch. 2
 * (combination of mean and low-frequency offsets). Linear mooring stiffness about the mean offset; Gaussian
 * response assumed for the extreme (slow-drift maxima are in fact more skewed; treat as an estimate).
 * @param {object} o
 * @param {number} o.mass total mass including added mass, kg
 * @param {number} o.stiffness mooring stiffness, N/m
 * @param {number} [o.damping] linear damping, N·s/m (or give dampingRatio)
 * @param {number} [o.dampingRatio] fraction of critical
 * @param {ArrayLike<number>} o.mu difference frequencies, rad/s (must resolve the resonance peak)
 * @param {ArrayLike<number>} o.SF slow-drift force spectrum on mu, N²·s/rad
 * @param {number} [o.meanForce=0] N
 * @param {number} [o.duration=10800] s
 * @returns {{naturalPeriod:number, meanOffset:number, sigma:number, sigmaWhiteNoise:number, extreme:number,
 *   nCycles:number, damping:number}} s, m, m, m, m, -, N·s/m
 */
export function slowDriftResponse({ mass, stiffness, damping, dampingRatio, mu, SF, meanForce = 0, duration = 10800 }) {
  needPos(mass, 'mass'); needPos(stiffness, 'stiffness'); needPos(duration, 'duration');
  const wn = Math.sqrt(stiffness / mass);
  const c = damping !== undefined ? damping : 2 * dampingRatio * Math.sqrt(stiffness * mass);
  if (!(c > 0)) throw new RangeError('damping (or dampingRatio) must be > 0');
  if (!mu || mu.length < 2 || SF.length !== mu.length) throw new RangeError('mu and SF must have equal length >= 2');
  // integrate on a grid refined around the resonance so that the peak is resolved for light damping
  const grid = new Set(Array.from(mu));
  const half = c / (2 * mass);
  for (let j = -40; j <= 40; j++) {
    const w = wn + (j / 8) * half;
    if (w > mu[0] && w < mu[mu.length - 1]) grid.add(w);
  }
  const ws = Array.from(grid).sort((p, q) => p - q);
  let s = 0; let prev = 0;
  for (let i = 0; i < ws.length; i++) {
    const w = ws[i];
    const f = lin(mu, SF, w) / ((stiffness - mass * w * w) ** 2 + (c * w) ** 2);
    if (i > 0) s += 0.5 * (f + prev) * (w - ws[i - 1]);
    prev = f;
  }
  const sigma = Math.sqrt(s);
  const Tn = (2 * Math.PI) / wn;
  const nCycles = duration / Tn;
  const meanOffset = meanForce / stiffness;
  return {
    naturalPeriod: Tn, meanOffset, sigma, sigmaWhiteNoise: Math.sqrt((Math.PI * lin(mu, SF, wn)) / (2 * c * stiffness)),
    extreme: meanOffset + Math.sign(meanOffset || 1) * sigma * Math.sqrt(2 * Math.log(Math.max(nCycles, Math.E))), nCycles, damping: c,
  };
}

function coef(c, angleDeg) {
  if (typeof c === 'function') return c(angleDeg);
  if (typeof c === 'number') return c;
  if (c && c.angleDeg && c.value) return lin(c.angleDeg, c.value, angleDeg);
  throw new RangeError('coefficient must be a number, a function of angle (deg) or { angleDeg: [], value: [] }');
}

/**
 * Wind or current force and yaw moment on a vessel in the OCIMF form, with the coefficient curves supplied by
 * the caller (numbers, functions of the relative angle in deg, or tables { angleDeg, value }):
 *   wind:    Fx = ½ rho Cx(theta) V² A_T,  Fy = ½ rho Cy(theta) V² A_L,  Mz = ½ rho Cxy(theta) V² A_L L
 *   current: Fx = ½ rho Cx(theta) V² L T,  Fy = ½ rho Cy(theta) V² L T,  Mz = ½ rho Cxy(theta) V² L² T
 * theta is passed to the coefficient functions unchanged; use the angle and sign convention of the coefficient
 * source.
 * Reference: OCIMF (1994) Prediction of Wind and Current Loads on VLCCs, 2nd ed. (force definitions);
 * No coefficient data are embedded.
 * @param {object} o
 * @param {string} o.kind 'wind' | 'current'
 * @param {number} o.V speed, m/s
 * @param {number} o.angleDeg relative angle in the convention of the coefficients, deg
 * @param {number|function|object} o.Cx longitudinal coefficient
 * @param {number|function|object} o.Cy lateral coefficient
 * @param {number|function|object} [o.Cxy=0] yaw moment coefficient
 * @param {number} [o.AT] transverse (head-on) windage area, m² (wind)
 * @param {number} [o.AL] longitudinal (broadside) windage area, m² (wind)
 * @param {number} o.L length between perpendiculars, m
 * @param {number} [o.T] draft, m (current)
 * @param {number} [o.rho] density, kg/m³ (default 1.225 wind, 1025 current)
 * @returns {{Fx:number, Fy:number, Mz:number}} N, N, N·m
 */
export function ocimfForce({ kind, V, angleDeg, Cx, Cy, Cxy = 0, AT, AL, L, T, rho }) {
  if (!(V >= 0)) throw new RangeError('V must be >= 0');
  needPos(L, 'L');
  if (kind === 'wind') {
    needPos(AT, 'AT'); needPos(AL, 'AL');
    const q = 0.5 * (rho === undefined ? 1.225 : rho) * V * V;
    return { Fx: q * coef(Cx, angleDeg) * AT, Fy: q * coef(Cy, angleDeg) * AL, Mz: q * coef(Cxy, angleDeg) * AL * L };
  }
  if (kind === 'current') {
    needPos(T, 'T');
    const q = 0.5 * (rho === undefined ? 1025 : rho) * V * V;
    return { Fx: q * coef(Cx, angleDeg) * L * T, Fy: q * coef(Cy, angleDeg) * L * T, Mz: q * coef(Cxy, angleDeg) * L * L * T };
  }
  throw new RangeError("kind must be 'wind' or 'current'");
}

/**
 * Wind or current force on a floater from projected areas and drag coefficients, resolving the flow into the
 * vessel axes (cross-flow principle): with the flow travelling towards direction theta,
 *   Fx = ½ rho Cfront Afront V² cos(theta)|cos(theta)|,  Fy = ½ rho Cside Aside V² sin(theta)|sin(theta)|,
 *   Mz = Fy·xSide - Fx·yFront  (centres of the projected areas).
 * Reference: DNV-RP-C205 (2019) §5 (wind force on projected area), §6 (cross-flow principle).
 * Approximation: no lift, no shielding or interaction beyond what the coefficients include. Coefficients are
 * inputs (typical bluff-body values 0.7–1.3 on projected area).
 * @param {object} o
 * @param {number} o.V speed, m/s
 * @param {number} o.directionDeg direction the flow travels towards, deg from +x towards +y
 * @param {number} o.Afront projected area seen along x, m²
 * @param {number} o.Aside projected area seen along y, m²
 * @param {number} o.Cfront drag coefficient, -
 * @param {number} o.Cside drag coefficient, -
 * @param {number} [o.xSide=0] x of the centre of the side area, m
 * @param {number} [o.yFront=0] y of the centre of the front area, m
 * @param {number} o.rho fluid density, kg/m³ (air ≈ 1.226 at 15 °C, sea water ≈ 1025)
 * @returns {{Fx:number, Fy:number, Mz:number}} N, N, N·m
 */
export function projectedAreaForce({ V, directionDeg, Afront, Aside, Cfront, Cside, xSide = 0, yFront = 0, rho }) {
  if (!(V >= 0)) throw new RangeError('V must be >= 0');
  needPos(rho, 'rho');
  if (!(Afront >= 0 && Aside >= 0)) throw new RangeError('Afront and Aside must be >= 0');
  const c = Math.cos(directionDeg * DEG); const s = Math.sin(directionDeg * DEG);
  const q = 0.5 * rho * V * V;
  const Fx = q * Cfront * Afront * c * Math.abs(c);
  const Fy = q * Cside * Aside * s * Math.abs(s);
  return { Fx, Fy, Mz: Fy * xSide - Fx * yFront };
}

/**
 * Total mean environmental force on a floater versus heading: sum of the wind, current and mean wave-drift
 * contributions returned by caller-supplied functions of the relative direction.
 * Reference: superposition of mean loads, DNV-OS-E301 (2021) ch. 2; API RP 2SK (2005).
 * @param {object} o
 * @param {number[]} o.headingsDeg vessel headings to evaluate, deg
 * @param {{name:string, directionDeg:number, force:function(number): {Fx:number, Fy:number, Mz:number}}[]} o.loads
 *   each load has a global direction (travelling towards, deg) and a function of the direction relative to the
 *   vessel (global direction minus heading, deg) returning forces in the vessel frame, N and N·m
 * @returns {{headingDeg:number[], Fx:number[], Fy:number[], Mz:number[], F:number[], worst:{headingDeg:number,
 *   F:number}, best:{headingDeg:number, F:number}}} vessel-frame forces N, moments N·m, resultant N
 */
export function meanEnvironmentalForce({ headingsDeg, loads }) {
  if (!Array.isArray(headingsDeg) || headingsDeg.length === 0) throw new RangeError('headingsDeg must be a non-empty array');
  if (!Array.isArray(loads) || loads.length === 0) throw new RangeError('loads must be a non-empty array');
  const out = { headingDeg: [], Fx: [], Fy: [], Mz: [], F: [], worst: null, best: null };
  for (const h of headingsDeg) {
    let Fx = 0; let Fy = 0; let Mz = 0;
    for (const ld of loads) {
      const r = ld.force(ld.directionDeg - h);
      Fx += r.Fx; Fy += r.Fy; Mz += r.Mz || 0;
    }
    const F = Math.hypot(Fx, Fy);
    out.headingDeg.push(h); out.Fx.push(Fx); out.Fy.push(Fy); out.Mz.push(Mz); out.F.push(F);
    if (!out.worst || F > out.worst.F) out.worst = { headingDeg: h, F };
    if (!out.best || F < out.best.F) out.best = { headingDeg: h, F };
  }
  return out;
}
