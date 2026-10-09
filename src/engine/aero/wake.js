// Wind-turbine (and tidal-turbine) wake models: single-wake deficits, added turbulence, superposition,
// rotor-disc averaging, a steady wind-farm flow solver, induction/blockage and the deep-array asymptote.
//
// Conventions: SI units; x is distance downstream of the wake-generating rotor (m), r the radial distance from
// the wake centreline (m); "deficit" is the fractional velocity deficit ΔU/U∞ (-). Wind direction `wdDeg` is the
// meteorological direction the wind comes FROM, clockwise from north (+y = north, +x = east).

const DEG = Math.PI / 180;
const CT_MAX = 0.999;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}
function clampCt(Ct) {
  need(Number.isFinite(Ct) && Ct >= 0, 'Ct must be a non-negative finite number');
  return Math.min(Ct, CT_MAX);
}

/**
 * Axial induction factor from the thrust coefficient by 1-D momentum theory, a = ½ (1 − √(1 − Ct)).
 * Reference: Rankine–Froude actuator-disc theory. Valid for Ct < 1 (Ct is limited to 0.999).
 * @param {number} Ct thrust coefficient (-)
 * @returns {number} axial induction factor a (-)
 */
export function inductionFromCt(Ct) {
  return 0.5 * (1 - Math.sqrt(1 - clampCt(Ct)));
}

/**
 * Jensen/Park top-hat wake: δ = (1 − √(1 − Ct)) / (1 + k x / R)², wake diameter Dw = D + 2 k x.
 * Reference: Jensen (1983) Risø-M-2411; Katić, Højstrup & Jensen (1986). k ≈ 0.04–0.05 offshore, 0.075 onshore.
 * The model is the mass balance π R² u_r + π (Rw² − R²) U∞ = π Rw² u with u_r = U∞ √(1 − Ct).
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.x downstream distance (m), x >= 0
 * @param {number} o.D rotor diameter (m)
 * @param {number} [o.k=0.05] wake decay constant (-)
 * @returns {{deficit:number, wakeDiameter:number}} deficit ΔU/U∞ (-) inside the wake, wake diameter (m)
 */
export function jensenDeficit({ Ct, x, D, k = 0.05 }) {
  finitePos(D, 'D');
  need(Number.isFinite(x) && x >= 0, 'x must be >= 0');
  need(Number.isFinite(k) && k >= 0, 'k must be >= 0');
  const c = clampCt(Ct);
  const R = D / 2;
  const f = 1 + (k * x) / R;
  return { deficit: (1 - Math.sqrt(1 - c)) / (f * f), wakeDiameter: D + 2 * k * x };
}

/**
 * Frandsen et al. (2006) storage-space-area top-hat wake.
 * Deficit from the momentum balance over a control volume around the rotor,
 *   Ct A0 U∞²/2 = A u (U∞ − u)  →  u/U∞ = ½ + ½ √(1 − 2 Ct A0/A),
 * with the wake area from Dw = D (β^(n/2) + α x/D)^(1/n), β = ½ (1 + √(1 − Ct))/√(1 − Ct).
 * Reference: Frandsen et al. (2006) Wind Energy 9:39–53. The expansion exponent n and the decay constant α are
 * arguments (the paper discusses n = 2 and n = 3; α of order 0.5–0.9 for n = 2 — calibrate, not a fixed constant).
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.x downstream distance (m)
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.alpha wake expansion constant α (-)
 * @param {number} [o.n=2] expansion shape exponent (-)
 * @returns {{deficit:number, wakeDiameter:number, beta:number}}
 */
export function frandsenDeficit({ Ct, x, D, alpha, n = 2 }) {
  finitePos(D, 'D');
  finitePos(alpha, 'alpha');
  finitePos(n, 'n');
  need(Number.isFinite(x) && x >= 0, 'x must be >= 0');
  const c = clampCt(Ct);
  const s = Math.sqrt(1 - c);
  const beta = (0.5 * (1 + s)) / s;
  const Dw = D * Math.pow(Math.pow(beta, n / 2) + (alpha * x) / D, 1 / n);
  const rad = Math.max(0, 1 - (2 * c * D * D) / (Dw * Dw));
  return { deficit: 0.5 - 0.5 * Math.sqrt(rad), wakeDiameter: Dw, beta };
}

/**
 * Wake growth rate of the Gaussian model from local streamwise turbulence intensity,
 * k* = 0.3837 TI + 0.003678. Reference: Niayifar & Porté-Agel (2016) Energies 9:741
 * (fitted for 0.065 < TI < 0.15).
 * @param {number} TI turbulence intensity σu/U (-)
 * @returns {number} k* (-)
 */
export function wakeGrowthFromTI(TI) {
  need(Number.isFinite(TI) && TI >= 0, 'TI must be >= 0');
  return 0.3837 * TI + 0.003678;
}

/**
 * Bastankhah & Porté-Agel (2014) Gaussian wake:
 *   ΔU/U∞ = (1 − √(1 − Ct / (8 (σ/D)²))) exp(−r² / (2 σ²)),  σ/D = k* x/D + ε,  ε = 0.2 √β,
 *   β = ½ (1 + √(1 − Ct)) / √(1 − Ct).
 * Reference: Bastankhah & Porté-Agel (2014) Renewable Energy 70:116–123.
 * The model conserves the momentum deficit: ∫ ΔU (U∞ − ΔU) dA = ½ Ct A0 U∞² wherever the radicand is real.
 * Near wake: where Ct > 8 (σ/D)² the radicand is negative; it is then set to zero (centre deficit = 1) and
 * `nearWake` is returned true — results there are outside the model's validity.
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.x downstream distance (m), > 0
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.kStar wake growth rate k* (-)
 * @param {number} [o.r=0] radial distance from the wake centreline (m)
 * @returns {{deficit:number, centreDeficit:number, sigma:number, epsilon:number, nearWake:boolean}}
 *   sigma in m, deficits (-)
 */
export function gaussianDeficit({ Ct, x, D, kStar, r = 0 }) {
  finitePos(D, 'D');
  need(Number.isFinite(x) && x >= 0, 'x must be >= 0');
  need(Number.isFinite(kStar) && kStar >= 0, 'kStar must be >= 0');
  const c = clampCt(Ct);
  const s = Math.sqrt(1 - c);
  const epsilon = 0.2 * Math.sqrt((0.5 * (1 + s)) / s);
  const sigD = (kStar * x) / D + epsilon;
  const rad = 1 - c / (8 * sigD * sigD);
  const nearWake = rad < 0;
  const centre = 1 - Math.sqrt(Math.max(0, rad));
  const sigma = sigD * D;
  return {
    deficit: centre * Math.exp((-r * r) / (2 * sigma * sigma)),
    centreDeficit: centre,
    sigma,
    epsilon,
    nearWake,
  };
}

/**
 * Wake-centre lateral deflection behind a yawed rotor, Bastankhah & Porté-Agel (2016) J. Fluid Mech. 806:506–541
 * (their §6–7): initial skew θc0 = (0.3 γ / cos γ)(1 − √(1 − Ct cos γ)); potential-core length
 * x0/D = cos γ (1 + √(1 − Ct)) / (√2 (4 α I + 2 β (1 − √(1 − Ct)))), α = 2.32, β = 0.154; far-wake deflection
 *   δ = θc0 x0 + (θc0 E0 / 5.2) √(σy0 σz0 / (ky kz M0)) ln[ (1.6 + √M0)(1.6 s − √M0) / ((1.6 − √M0)(1.6 s + √M0)) ],
 * s = √(σy σz / (σy0 σz0)), σy0 = D cos γ/√8, σz0 = D/√8, C0 = 1 − √(1 − Ct), M0 = C0 (2 − C0),
 * E0 = C0² − 3 e^(1/12) C0 + 3 e^(1/3). Linear growth of the deflection is used inside the potential core.
 * Validity: |γ| ≲ 30°, flat terrain, neutral boundary layer.
 * @param {object} o
 * @param {number} o.Ct thrust coefficient of the yawed rotor (-)
 * @param {number} o.yawDeg yaw misalignment γ (deg)
 * @param {number} o.x downstream distance (m)
 * @param {number} o.D rotor diameter (m)
 * @param {number} o.TI ambient streamwise turbulence intensity (-)
 * @param {number} o.ky lateral wake growth rate (-)
 * @param {number} [o.kz=ky] vertical wake growth rate (-)
 * @param {number} [o.alpha=2.32] potential-core constant α (-)
 * @param {number} [o.beta=0.154] potential-core constant β (-)
 * @returns {{deflection:number, theta0:number, x0:number, sigmaY:number, sigmaZ:number}} deflection (m, sign of
 *   the yaw angle), initial skew angle (rad), core length (m), wake widths (m)
 */
export function yawDeflection({ Ct, yawDeg, x, D, TI, ky, kz = ky, alpha = 2.32, beta = 0.154 }) {
  finitePos(D, 'D');
  finitePos(ky, 'ky');
  finitePos(kz, 'kz');
  need(Number.isFinite(x) && x >= 0, 'x must be >= 0');
  need(Number.isFinite(TI) && TI >= 0, 'TI must be >= 0');
  need(Math.abs(yawDeg) < 60, 'yawDeg must be within ±60 deg');
  const c = clampCt(Ct);
  const g = yawDeg * DEG;
  const cg = Math.cos(g);
  const s = Math.sqrt(1 - c);
  const theta0 = ((0.3 * g) / cg) * (1 - Math.sqrt(1 - c * cg));
  const x0 = (D * cg * (1 + s)) / (Math.SQRT2 * (4 * alpha * TI + 2 * beta * (1 - s)));
  const sy0 = (D * cg) / Math.sqrt(8);
  const sz0 = D / Math.sqrt(8);
  if (x <= x0) return { deflection: theta0 * x, theta0, x0, sigmaY: sy0, sigmaZ: sz0 };
  const sy = ky * (x - x0) + sy0;
  const sz = kz * (x - x0) + sz0;
  const C0 = 1 - s;
  const M0 = C0 * (2 - C0);
  const E0 = C0 * C0 - 3 * Math.exp(1 / 12) * C0 + 3 * Math.exp(1 / 3);
  const q = Math.sqrt((sy * sz) / (sy0 * sz0));
  const rM = Math.sqrt(M0);
  const lnTerm = Math.log(((1.6 + rM) * (1.6 * q - rM)) / ((1.6 - rM) * (1.6 * q + rM)));
  const deflection = theta0 * x0 + ((theta0 * E0) / 5.2) * Math.sqrt((sy0 * sz0) / (ky * kz * M0)) * lnTerm;
  return { deflection, theta0, x0, sigmaY: sy, sigmaZ: sz };
}

/**
 * Time-averaged Gaussian wake seen in a fixed frame when the wake meanders as a passive tracer with Gaussian
 * centre-position statistics (simplified dynamic-wake-meandering picture, after Larsen et al. 2008):
 * convolution of two Gaussians gives σ_eff² = σ² + σ_m² in each direction and a centre deficit scaled by
 * σ² / (σ_eff,y σ_eff,z). The meander standard deviation can be estimated from Taylor's frozen-turbulence
 * hypothesis σ_m ≈ I_v x (lateral) and I_w x (vertical) for large-scale turbulence; pass it explicitly.
 * @param {object} o
 * @param {number} o.centreDeficit instantaneous (meandering-frame) centre deficit (-)
 * @param {number} o.sigma meandering-frame wake width (m)
 * @param {number} o.sigmaMeanderY lateral meander standard deviation (m)
 * @param {number} [o.sigmaMeanderZ=0] vertical meander standard deviation (m)
 * @returns {{centreDeficit:number, sigmaY:number, sigmaZ:number}}
 */
export function meanderedGaussian({ centreDeficit, sigma, sigmaMeanderY, sigmaMeanderZ = 0 }) {
  finitePos(sigma, 'sigma');
  need(sigmaMeanderY >= 0 && sigmaMeanderZ >= 0, 'meander standard deviations must be >= 0');
  const sy = Math.hypot(sigma, sigmaMeanderY);
  const sz = Math.hypot(sigma, sigmaMeanderZ);
  return { centreDeficit: (centreDeficit * sigma * sigma) / (sy * sz), sigmaY: sy, sigmaZ: sz };
}

/**
 * Wake-added turbulence intensity, Frandsen: I_add = 1 / (1.5 + 0.8 (x/D)/√Ct).
 * Reference: Frandsen (2007) Risø-R-1188; IEC 61400-1:2019 Annex E (ed. 3 Amd 1 Annex D), written as σ/V.
 * Combine with ambient as I = √(I0² + I_add²).
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.xOverD spacing in rotor diameters (-), > 0
 * @returns {number} added turbulence intensity (-)
 */
export function addedTurbulenceFrandsen({ Ct, xOverD }) {
  finitePos(xOverD, 'xOverD');
  finitePos(Ct, 'Ct');
  return 1 / (1.5 + (0.8 * xOverD) / Math.sqrt(Ct));
}

/**
 * Wake-added turbulence intensity, Crespo & Hernández (1996): I_add = 0.73 a^0.8325 I0^0.0325 (x/D)^−0.32,
 * a = ½ (1 − √(1 − Ct)). Reference: Crespo & Hernández (1996) J. Wind Eng. Ind. Aerodyn. 61:71–85.
 * Fitted for 5 < x/D < 15, 0.07 < I0 < 0.14, 0.1 < a < 0.4; `inRange` reports whether the call is inside.
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.xOverD downstream distance in diameters (-)
 * @param {number} o.I0 ambient turbulence intensity (-)
 * @returns {{added:number, inRange:boolean}}
 */
export function addedTurbulenceCrespo({ Ct, xOverD, I0 }) {
  finitePos(xOverD, 'xOverD');
  finitePos(I0, 'I0');
  const a = inductionFromCt(Ct);
  const added = 0.73 * Math.pow(a, 0.8325) * Math.pow(I0, 0.0325) * Math.pow(xOverD, -0.32);
  const inRange = xOverD >= 5 && xOverD <= 15 && I0 >= 0.07 && I0 <= 0.14 && a >= 0.1 && a <= 0.4;
  return { added, inRange };
}

/**
 * Effective turbulence intensity for fatigue in a wind farm, IEC 61400-1:2019 Annex E (Frandsen):
 *   I_eff = [ (1 − N p_w) I0^m + p_w Σ_i I_T,i^m ]^(1/m),  I_T,i = √(I_add,i² + I0²),
 *   I_add,i = 1 / (1.5 + 0.8 d_i/√Ct),  p_w = 0.06 (uniform wind-direction distribution).
 * @param {object} o
 * @param {number} o.I0 ambient (representative) turbulence intensity σ̂/V_hub (-)
 * @param {{distanceD:number, Ct:number}[]} o.neighbours neighbouring turbines: spacing (rotor diameters), Ct (-)
 * @param {number} o.m Wöhler exponent of the component material (-), e.g. 3–5 steel, 10–12 composites
 * @param {number} [o.pw=0.06] probability of a wake condition per neighbour (-)
 * @returns {{Ieff:number, wakeTI:number[]}} effective TI (-), in-wake TI per neighbour (-)
 */
export function effectiveTurbulence({ I0, neighbours, m, pw = 0.06 }) {
  finitePos(I0, 'I0');
  finitePos(m, 'm');
  need(Array.isArray(neighbours), 'neighbours must be an array');
  const N = neighbours.length;
  need(N * pw <= 1, 'neighbours.length * pw must not exceed 1');
  let sum = (1 - N * pw) * Math.pow(I0, m);
  const wakeTI = [];
  for (const nb of neighbours) {
    const Ia = addedTurbulenceFrandsen({ Ct: nb.Ct, xOverD: nb.distanceD });
    const It = Math.hypot(Ia, I0);
    wakeTI.push(It);
    sum += pw * Math.pow(It, m);
  }
  return { Ieff: Math.pow(sum, 1 / m), wakeTI };
}

/**
 * Combine fractional velocity deficits of several wakes at one point.
 *  - 'linear':  δ = Σ δ_i (Lissaman 1979)
 *  - 'rss':     δ = √(Σ δ_i²) (Katić et al. 1986)
 *  - 'energy':  1 − (1 − δ)² = Σ [1 − (1 − δ_i)²] (sum of kinetic-energy deficits)
 *  - 'max':     δ = max δ_i
 * The result is limited to 1 (no flow reversal).
 * @param {number[]|Float64Array} deficits individual deficits (-)
 * @param {'linear'|'rss'|'energy'|'max'} [method='rss']
 * @returns {number} combined deficit (-)
 */
export function superpose(deficits, method = 'rss') {
  let acc = 0;
  for (let i = 0; i < deficits.length; i++) {
    const d = deficits[i];
    if (method === 'linear') acc += d;
    else if (method === 'rss') acc += d * d;
    else if (method === 'energy') acc += 1 - (1 - d) * (1 - d);
    else if (method === 'max') acc = Math.max(acc, d);
    else throw new RangeError(`method '${method}' is not one of linear, rss, energy, max`);
  }
  return finishSuperposition(acc, method);
}

function finishSuperposition(acc, method) {
  let d = acc;
  if (method === 'rss') d = Math.sqrt(acc);
  else if (method === 'energy') d = 1 - Math.sqrt(Math.max(0, 1 - acc));
  return Math.min(1, Math.max(0, d));
}

/**
 * Area of intersection of two circles (lens formula).
 * @param {number} R1 radius 1 (m)
 * @param {number} R2 radius 2 (m)
 * @param {number} d centre distance (m)
 * @returns {number} overlap area (m²)
 */
export function circleOverlapArea(R1, R2, d) {
  need(R1 >= 0 && R2 >= 0 && d >= 0, 'R1, R2 and d must be >= 0');
  if (d >= R1 + R2) return 0;
  const rMin = Math.min(R1, R2);
  if (d <= Math.abs(R1 - R2)) return Math.PI * rMin * rMin;
  const a1 = Math.acos((d * d + R1 * R1 - R2 * R2) / (2 * d * R1));
  const a2 = Math.acos((d * d + R2 * R2 - R1 * R1) / (2 * d * R2));
  const k = 0.5 * Math.sqrt((-d + R1 + R2) * (d + R1 - R2) * (d - R1 + R2) * (d + R1 + R2));
  return R1 * R1 * a1 + R2 * R2 * a2 - k;
}

/**
 * Rotor-disc average of a top-hat wake deficit: δ̄ = δ · A_overlap / A_rotor (exact for a top-hat profile).
 * @param {object} o
 * @param {number} o.deficit in-wake deficit (-)
 * @param {number} o.wakeRadius wake radius at the rotor plane (m)
 * @param {number} o.rotorRadius downstream rotor radius (m)
 * @param {number} o.offset distance between wake centreline and rotor centre (m)
 * @returns {{deficit:number, overlapFraction:number}}
 */
export function tophatRotorAverage({ deficit, wakeRadius, rotorRadius, offset }) {
  finitePos(rotorRadius, 'rotorRadius');
  const f = circleOverlapArea(wakeRadius, rotorRadius, Math.abs(offset)) / (Math.PI * rotorRadius * rotorRadius);
  return { deficit: deficit * f, overlapFraction: f };
}

/**
 * Equal-weight quadrature points on the unit disc: nR equal-area rings (radius √((i+½)/nR)) × nTheta points,
 * successive rings staggered by half a step. Exact for constants; error O(1/nR²) for smooth integrands.
 * @param {number} [nR=3] rings
 * @param {number} [nTheta=8] points per ring
 * @returns {{y:Float64Array, z:Float64Array, n:number}} coordinates on the unit disc (-)
 */
export function rotorPoints(nR = 3, nTheta = 8) {
  need(Number.isInteger(nR) && nR >= 1 && Number.isInteger(nTheta) && nTheta >= 1, 'nR and nTheta must be integers >= 1');
  const n = nR * nTheta;
  const y = new Float64Array(n);
  const z = new Float64Array(n);
  let p = 0;
  for (let i = 0; i < nR; i++) {
    const r = Math.sqrt((i + 0.5) / nR);
    for (let j = 0; j < nTheta; j++) {
      const th = (2 * Math.PI * (j + 0.5 * (i % 2))) / nTheta;
      y[p] = r * Math.cos(th);
      z[p] = r * Math.sin(th);
      p++;
    }
  }
  return { y, z, n };
}

/**
 * Rotor-disc average of a Gaussian wake deficit by quadrature over the rotor.
 * @param {object} o
 * @param {number} o.centreDeficit deficit on the wake centreline (-)
 * @param {number} o.sigma wake width σ at the rotor plane (m)
 * @param {number} o.dy lateral offset of the rotor centre from the wake centreline (m)
 * @param {number} [o.dz=0] vertical offset (m)
 * @param {number} o.rotorRadius downstream rotor radius (m)
 * @param {{y:Float64Array, z:Float64Array, n:number}} [o.points] unit-disc points from `rotorPoints`
 * @returns {number} rotor-averaged deficit (-)
 */
export function gaussianRotorAverage({ centreDeficit, sigma, dy, dz = 0, rotorRadius, points = DEFAULT_POINTS }) {
  finitePos(sigma, 'sigma');
  finitePos(rotorRadius, 'rotorRadius');
  const inv = 1 / (2 * sigma * sigma);
  let s = 0;
  for (let p = 0; p < points.n; p++) {
    const yy = dy + rotorRadius * points.y[p];
    const zz = dz + rotorRadius * points.z[p];
    s += Math.exp(-(yy * yy + zz * zz) * inv);
  }
  return (centreDeficit * s) / points.n;
}

const DEFAULT_POINTS = rotorPoints(3, 8);

/**
 * Rotor-equivalent wind speed with vertical shear, IEC 61400-12-1:2017 (rotor equivalent wind speed):
 *   v_eq = (Σ v_i³ A_i / A)^(1/3), with A_i the area of the horizontal rotor strip at height z_i.
 * @param {object} o
 * @param {(z:number)=>number} o.profile wind speed (m/s) as a function of height above surface z (m)
 * @param {number} o.hubHeight hub height (m)
 * @param {number} o.D rotor diameter (m)
 * @param {number} [o.nSeg=40] number of horizontal strips
 * @returns {{rews:number, hubSpeed:number, ratio:number}} v_eq (m/s), hub-height speed (m/s), v_eq/v_hub (-)
 */
export function rotorEquivalentWindSpeed({ profile, hubHeight, D, nSeg = 40 }) {
  finitePos(D, 'D');
  finitePos(hubHeight, 'hubHeight');
  need(hubHeight > D / 2, 'hubHeight must exceed the rotor radius');
  need(typeof profile === 'function', 'profile must be a function of height');
  const R = D / 2;
  // area of the disc below normalised height t = (z − hub)/R
  const below = (t) => R * R * (Math.PI / 2 + Math.asin(t) + t * Math.sqrt(Math.max(0, 1 - t * t)));
  let sum = 0;
  for (let i = 0; i < nSeg; i++) {
    const t0 = -1 + (2 * i) / nSeg;
    const t1 = -1 + (2 * (i + 1)) / nSeg;
    const v = profile(hubHeight + R * 0.5 * (t0 + t1));
    sum += v * v * v * (below(t1) - below(t0));
  }
  const rews = Math.cbrt(sum / (Math.PI * R * R));
  const hubSpeed = profile(hubHeight);
  return { rews, hubSpeed, ratio: rews / hubSpeed };
}

/**
 * Axial velocity on the rotor axis induced by a rotor, from the semi-infinite vortex-cylinder model:
 *   u/U∞ = 1 − a (1 + ξ / √(1 + ξ²)),  ξ = x/R  (x < 0 upstream, x > 0 downstream).
 * Gives 1 at ξ → −∞, 1 − a at the rotor and 1 − 2a far downstream. Reference: Branlard & Gaunaa (2015)
 * Wind Energy 18:1973–1987; Medici et al. (2011). Valid on the axis for a uniformly loaded disc, a < 0.5.
 * @param {object} o
 * @param {number} o.a axial induction factor (-)
 * @param {number} o.x axial position relative to the rotor (m), negative upstream
 * @param {number} o.R rotor radius (m)
 * @returns {number} u/U∞ (-)
 */
export function inductionOnAxis({ a, x, R }) {
  finitePos(R, 'R');
  need(Number.isFinite(a) && a >= 0 && a < 0.5, 'a must be in [0, 0.5)');
  const xi = x / R;
  return 1 - a * (1 + xi / Math.sqrt(1 + xi * xi));
}

/**
 * Far-field axial velocity perturbation upstream of a rotor at an off-axis point. The vortex cylinder is replaced
 * by its far-field equivalent, a point source of strength Q = 2π a U∞ R² (the on-axis expansion of the cylinder
 * formula for ξ → −∞ is −a U∞ R²/(2 x²), which is the axial velocity of that source):
 *   Δu/U∞ = (a R²/2) x / (x² + r²)^(3/2)    (negative, i.e. a slow-down, for x < 0).
 * The magnitude is limited to the on-axis vortex-cylinder value at the same x so that the near field stays bounded.
 * Used for global-blockage estimates (cf. the Rankine half-body approach of Gribben & Hawkes 2019).
 * Returns 0 for x > 0 (points behind the rotor plane are handled by the wake models).
 * @param {object} o
 * @param {number} o.a axial induction factor (-)
 * @param {number} o.x axial position relative to the rotor (m), negative upstream
 * @param {number} o.r radial distance from the rotor axis (m)
 * @param {number} o.R rotor radius (m)
 * @returns {number} Δu/U∞ (-), <= 0
 */
export function inductionUpstream({ a, x, r, R }) {
  finitePos(R, 'R');
  need(Number.isFinite(a) && a >= 0 && a < 0.5, 'a must be in [0, 0.5)');
  if (x > 0) return 0;
  const d2 = x * x + r * r;
  if (d2 === 0) return -a;
  const src = (0.5 * a * R * R * x) / (d2 * Math.sqrt(d2));
  const cap = -(1 - inductionOnAxis({ a, x, R }));
  return Math.max(src, cap);
}

/**
 * Infinite wind-farm (deep-array) asymptote by Frandsen's two-layer roughness model.
 * Below hub height the profile is logarithmic on the surface roughness z0, above on an effective roughness z00;
 * the turbines add a distributed stress ρ c_t U_h² with c_t = π Ct / (8 s_x s_y). Matching velocity at hub height:
 *   z00 = h exp(−κ / √(c_t + (κ / ln(h/z0))²)).
 * Reference: Frandsen (1992) J. Wind Eng. Ind. Aerodyn. 39:251–265; Frandsen et al. (2006).
 * The hub-height speed ratio additionally assumes the same wind speed with and without the farm at a reference
 * height H (boundary-layer or blending height, an argument):
 *   U_h / U_h0 = [ln(h/z00) / ln(H/z00)] / [ln(h/z0) / ln(H/z0)].
 * @param {object} o
 * @param {number} o.Ct thrust coefficient (-)
 * @param {number} o.sx streamwise spacing (rotor diameters)
 * @param {number} o.sy cross-stream spacing (rotor diameters)
 * @param {number} o.hubHeight hub height h (m)
 * @param {number} o.z0 surface roughness length (m); open sea ≈ 1e-4–2e-4
 * @param {number} o.H reference height where the flow is unaffected by the farm (m), > hubHeight
 * @param {number} [o.kappa=0.4] von Kármán constant (-)
 * @returns {{z00:number, ct:number, speedRatio:number, powerRatio:number}} effective roughness (m), distributed
 *   thrust coefficient (-), U_h/U_h0 (-), its cube (-)
 */
export function deepArrayFrandsen({ Ct, sx, sy, hubHeight, z0, H, kappa = 0.4 }) {
  finitePos(sx, 'sx');
  finitePos(sy, 'sy');
  finitePos(hubHeight, 'hubHeight');
  finitePos(z0, 'z0');
  need(z0 < hubHeight, 'z0 must be smaller than hubHeight');
  need(H > hubHeight, 'H must exceed hubHeight');
  need(Ct >= 0, 'Ct must be >= 0');
  const ct = (Math.PI * Ct) / (8 * sx * sy);
  const k0 = kappa / Math.log(hubHeight / z0);
  const z00 = hubHeight * Math.exp(-kappa / Math.sqrt(ct + k0 * k0));
  const speedRatio =
    Math.log(hubHeight / z00) / Math.log(H / z00) / (Math.log(hubHeight / z0) / Math.log(H / z0));
  return { z00, ct, speedRatio, powerRatio: speedRatio ** 3 };
}

/**
 * Evaluate a curve given either as a function of wind speed or as a table { u: [], v: [] } (linear interpolation,
 * zero outside the table).
 */
function evalCurve(curve, U, name) {
  if (typeof curve === 'function') return curve(U);
  need(curve && curve.u && curve.v && curve.u.length === curve.v.length && curve.u.length > 1,
    `${name} must be a function or a table { u: [], v: [] }`);
  const { u, v } = curve;
  const n = u.length;
  if (U < u[0] || U > u[n - 1]) return 0;
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (u[mid] <= U) lo = mid;
    else hi = mid;
  }
  const t = (U - u[lo]) / (u[hi] - u[lo]);
  return v[lo] + t * (v[hi] - v[lo]);
}

/**
 * Steady wind-farm flow for one wind direction and one free-stream speed.
 * Turbines are processed upstream → downstream; each rotor sees the free stream reduced by the superposed
 * rotor-averaged deficits of all upstream wakes; its thrust coefficient and power follow from the supplied
 * curves at the rotor-averaged speed; local turbulence intensity (ambient + wake-added) sets the wake growth of
 * the Gaussian model. Without blockage the sweep is explicit (information only travels downstream); with
 * blockage the upstream induction of every rotor is added and the sweep is repeated until the speeds converge.
 *
 * Options:
 *  - model 'gaussian' (Bastankhah & Porté-Agel 2014, rotor quadrature) or 'jensen' (top-hat, exact overlap);
 *  - superposition 'linear' | 'rss' | 'energy' | 'max';
 *  - deficitReference 'freestream' (deficit × U∞, Katić) or 'local' (deficit × the wake-generating rotor's own
 *    inflow speed, Niayifar & Porté-Agel 2016);
 *  - turbulenceModel 'crespo' | 'frandsen' | 'none': the largest added TI among upstream wakes, weighted by the
 *    fraction of the rotor inside that wake (wake diameter 4σ for Gaussian, D + 2kx for Jensen);
 *  - reflect: image turbines mirrored in the ground/sea surface (Lissaman 1979) contribute like real wakes;
 *  - blockage: global-blockage estimate by summed upstream induction (`inductionUpstream`).
 *
 * @param {object} o
 * @param {{x:number, y:number, D:number, hubHeight:number, ct:Function|{u:number[],v:number[]},
 *          power:Function|{u:number[],v:number[]}}[]} o.turbines positions (m), rotor diameter (m), hub height
 *          (m), Ct(U) (-) and electrical power P(U) (W) as functions of rotor-averaged wind speed (m/s)
 * @param {number} o.U free-stream wind speed at hub height (m/s)
 * @param {number} o.wdDeg wind direction, from, clockwise from north (deg)
 * @param {number} [o.TI=0.06] ambient turbulence intensity (-)
 * @param {'gaussian'|'jensen'} [o.model='gaussian']
 * @param {number} [o.k=0.05] Jensen wake decay constant (-)
 * @param {number} [o.kStar] fixed Gaussian growth rate; default: from local TI (Niayifar & Porté-Agel)
 * @param {'linear'|'rss'|'energy'|'max'} [o.superposition='rss']
 * @param {'freestream'|'local'} [o.deficitReference='freestream']
 * @param {'crespo'|'frandsen'|'none'} [o.turbulenceModel='crespo']
 * @param {boolean} [o.reflect=false]
 * @param {boolean} [o.blockage=false]
 * @param {number} [o.maxIter=20] blockage iterations
 * @param {number} [o.tol=1e-6] convergence tolerance on speed (relative)
 * @param {{y:Float64Array,z:Float64Array,n:number}} [o.points] rotor quadrature points (`rotorPoints`)
 * @returns {{U:Float64Array, TI:Float64Array, Ct:Float64Array, power:Float64Array, deficit:Float64Array,
 *   blockage:Float64Array, order:Int32Array, totalPower:number, grossPower:number, efficiency:number,
 *   iterations:number}} rotor-averaged speed (m/s), local TI (-), Ct (-), power (W), combined wake deficit (-),
 *   blockage deficit (-) per turbine (input order); upstream→downstream order; farm power (W); power without
 *   interaction (W); their ratio (-)
 */
export function farmFlow({
  turbines,
  U,
  wdDeg,
  TI = 0.06,
  model = 'gaussian',
  k = 0.05,
  kStar,
  superposition = 'rss',
  deficitReference = 'freestream',
  turbulenceModel = 'crespo',
  reflect = false,
  blockage = false,
  maxIter = 20,
  tol = 1e-6,
  points = DEFAULT_POINTS,
}) {
  need(Array.isArray(turbines) && turbines.length > 0, 'turbines must be a non-empty array');
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  need(Number.isFinite(wdDeg), 'wdDeg must be finite');
  need(TI >= 0, 'TI must be >= 0');
  need(model === 'gaussian' || model === 'jensen', "model must be 'gaussian' or 'jensen'");
  need(['linear', 'rss', 'energy', 'max'].includes(superposition), 'superposition is not a known method');
  need(['crespo', 'frandsen', 'none'].includes(turbulenceModel), 'turbulenceModel is not a known model');
  const n = turbines.length;
  const th = wdDeg * DEG;
  const sn = Math.sin(th);
  const cs = Math.cos(th);
  const s = new Float64Array(n); // downstream coordinate
  const c = new Float64Array(n); // crosswind coordinate
  for (let i = 0; i < n; i++) {
    const t = turbines[i];
    finitePos(t.D, `turbines[${i}].D`);
    need(Number.isFinite(t.x) && Number.isFinite(t.y), `turbines[${i}] needs finite x and y`);
    s[i] = -(t.x * sn + t.y * cs);
    c[i] = t.x * cs - t.y * sn;
  }
  const order = Int32Array.from(Array.from({ length: n }, (_, i) => i).sort((p, q) => s[p] - s[q]));
  const Urot = new Float64Array(n).fill(U);
  const TIloc = new Float64Array(n).fill(TI);
  const Ct = new Float64Array(n);
  const P = new Float64Array(n);
  const def = new Float64Array(n);
  const blk = new Float64Array(n);
  let gross = 0;
  for (let i = 0; i < n; i++) gross += evalCurve(turbines[i].power, U, 'power');

  const sweep = () => {
    for (let oi = 0; oi < n; oi++) {
      const i = order[oi];
      const ti = turbines[i];
      const Ri = ti.D / 2;
      const hi = ti.hubHeight ?? 0;
      let acc = 0;
      let addMax = 0;
      for (let oj = 0; oj < oi; oj++) {
        const j = order[oj];
        const dx = s[i] - s[j];
        if (dx <= 1e-9 || Ct[j] <= 0) continue;
        const tj = turbines[j];
        const dy = c[i] - c[j];
        const hj = tj.hubHeight ?? 0;
        const scale = deficitReference === 'local' ? Urot[j] / (U || 1) : 1;
        const ctj = Math.min(Ct[j], CT_MAX);
        const sq = Math.sqrt(1 - ctj);
        const areaI = Math.PI * Ri * Ri;
        // single-wake formulas are written out here (same expressions as jensenDeficit / gaussianDeficit)
        // so that the pair loop does not allocate
        let wakeR;
        let centre = 0;
        let inv2s2 = 0;
        if (model === 'jensen') {
          const f = 1 + (2 * k * dx) / tj.D;
          wakeR = tj.D / 2 + k * dx;
          centre = (1 - sq) / (f * f);
        } else {
          const ks = kStar ?? 0.3837 * TIloc[j] + 0.003678;
          const sigD = (ks * dx) / tj.D + 0.2 * Math.sqrt((0.5 * (1 + sq)) / sq);
          const sigma = sigD * tj.D;
          centre = 1 - Math.sqrt(Math.max(0, 1 - ctj / (8 * sigD * sigD)));
          inv2s2 = 1 / (2 * sigma * sigma);
          wakeR = 2 * sigma;
        }
        const nImg = reflect && hj > 0 ? 2 : 1;
        for (let im = 0; im < nImg; im++) {
          const dz = im === 0 ? hi - hj : hi + hj;
          let d;
          if (model === 'jensen') d = (centre * circleOverlapArea(wakeR, Ri, Math.hypot(dy, dz))) / areaI;
          else {
            let sum = 0;
            for (let p = 0; p < points.n; p++) {
              const yy = dy + Ri * points.y[p];
              const zz = dz + Ri * points.z[p];
              sum += Math.exp(-(yy * yy + zz * zz) * inv2s2);
            }
            d = (centre * sum) / points.n;
          }
          d *= scale;
          if (superposition === 'linear') acc += d;
          else if (superposition === 'rss') acc += d * d;
          else if (superposition === 'energy') acc += 1 - (1 - d) * (1 - d);
          else acc = Math.max(acc, d);
        }
        if (turbulenceModel !== 'none') {
          const frac = circleOverlapArea(wakeR, Ri, Math.hypot(dy, hi - hj)) / areaI;
          if (frac > 0) {
            const xOverD = dx / tj.D;
            const add = turbulenceModel === 'frandsen'
              ? 1 / (1.5 + (0.8 * xOverD) / Math.sqrt(ctj))
              : 0.73 * Math.pow(0.5 * (1 - sq), 0.8325) * Math.pow(Math.max(TI, 1e-6), 0.0325) * Math.pow(xOverD, -0.32);
            addMax = Math.max(addMax, add * frac);
          }
        }
      }
      def[i] = finishSuperposition(acc, superposition);
      TIloc[i] = Math.hypot(TI, addMax);
      Urot[i] = U * Math.max(0, 1 - def[i] - blk[i]);
      Ct[i] = Math.max(0, evalCurve(ti.ct, Urot[i], 'ct'));
      P[i] = evalCurve(ti.power, Urot[i], 'power');
    }
  };

  sweep();
  let iterations = 1;
  if (blockage) {
    const prev = new Float64Array(n);
    for (; iterations < maxIter; iterations++) {
      prev.set(Urot);
      for (let i = 0; i < n; i++) {
        let b = 0;
        for (let j = 0; j < n; j++) {
          if (j === i || Ct[j] <= 0) continue;
          const dx = s[i] - s[j];
          if (dx > 0) continue;
          const tj = turbines[j];
          const r = Math.hypot(c[i] - c[j], (turbines[i].hubHeight ?? 0) - (tj.hubHeight ?? 0));
          b -= inductionUpstream({ a: Math.min(0.499, inductionFromCt(Ct[j])), x: dx, r, R: tj.D / 2 });
        }
        blk[i] = b;
      }
      sweep();
      let err = 0;
      for (let i = 0; i < n; i++) err = Math.max(err, Math.abs(Urot[i] - prev[i]));
      if (err <= tol * Math.max(U, 1e-9)) {
        iterations++;
        break;
      }
    }
  }
  let total = 0;
  for (let i = 0; i < n; i++) total += P[i];
  return {
    U: Urot,
    TI: TIloc,
    Ct,
    power: P,
    deficit: def,
    blockage: blk,
    order,
    totalPower: total,
    grossPower: gross,
    efficiency: gross > 0 ? total / gross : 1,
    iterations,
  };
}
