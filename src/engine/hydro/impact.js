// Fluid impact and related screening models: water entry (von Kármán, Wagner), slamming coefficients,
// slamming and deck-wetness probability, dam-break green-water flow, sloshing periods, wave-in-deck force
// and vortex-induced-motion screening.

import { G, needPos } from './_wavehelp.js';

const DEG = Math.PI / 180;

/**
 * Constant-velocity water entry of a two-dimensional wedge of deadrise angle beta (per unit length):
 *   von Kármán (no pile-up): wetted half-width c = V t / tan(beta),
 *   Wagner (pile-up factor pi/2): c = (pi/2) V t / tan(beta),
 *   added mass m_a = ½ rho pi c², force F = d(m_a V)/dt = rho pi c V dc/dt,
 * so F_vK = pi rho V³ t / tan²(beta) and F_Wagner = (pi/2)² F_vK. Wagner peak pressure at the spray root:
 *   p_max = ½ rho (dc/dt)² = ½ rho V² pi² / (4 tan²(beta)), i.e. Cp = pi²/(4 tan² beta) (stagnation-pressure
 *   contribution of order 1 neglected, as in the outer Wagner solution).
 * Reference: von Kármán (1929) NACA TN 321; Wagner (1932) ZAMM 12; Faltinsen (1990) ch. 9 (water impact and entry).
 * Valid for small deadrise (about 4° < beta < 30–40°), rigid body, no air cushion, before chine wetting
 * (c < half-beam); gravity and hydroelasticity neglected.
 * @param {object} o
 * @param {number} o.deadriseDeg deadrise angle, deg (0 < beta < 90)
 * @param {number} o.V entry velocity, m/s
 * @param {number} o.t time after first contact, s
 * @param {number} [o.halfBeam] half-beam of the wedge, m (limits the wetted width)
 * @param {number} [o.rho=1025] kg/m³
 * @returns {{vonKarman:{c:number, force:number, addedMass:number}, wagner:{c:number, force:number, addedMass:number,
 *   peakPressure:number, Cp:number}, ratio:number, chineWetted:boolean}} half-widths m, forces N/m, added mass kg/m,
 *   pressure Pa, ratio F_Wagner/F_vK (-)
 */
export function wedgeWaterEntry({ deadriseDeg, V, t, halfBeam, rho = 1025 }) {
  if (!(deadriseDeg > 0 && deadriseDeg < 90)) throw new RangeError('deadriseDeg must be within (0, 90)');
  if (!(V >= 0)) throw new RangeError('V must be >= 0');
  if (!(t >= 0)) throw new RangeError('t must be >= 0');
  const tb = Math.tan(deadriseDeg * DEG);
  const one = (factor) => {
    const c = (factor * V * t) / tb;
    const dc = (factor * V) / tb;
    return { c, force: rho * Math.PI * c * V * dc, addedMass: 0.5 * rho * Math.PI * c * c, dc };
  };
  const vk = one(1); const wg = one(Math.PI / 2);
  const Cp = (Math.PI * Math.PI) / (4 * tb * tb);
  return {
    vonKarman: { c: vk.c, force: vk.force, addedMass: vk.addedMass },
    wagner: { c: wg.c, force: wg.force, addedMass: wg.addedMass, peakPressure: 0.5 * rho * V * V * Cp, Cp },
    ratio: (Math.PI / 2) ** 2, chineWetted: halfBeam !== undefined && wg.c >= halfBeam,
  };
}

/**
 * Space-average slamming pressure coefficient C_Pa in p_s = ½ rho C_Pa v²:
 *   smooth circular cylinder: C_Pa >= 5.15;
 *   flat bottom (deadrise below 4°), allowing for cushioning and three-dimensional effects: C_Pa >= 2 pi;
 *   wedge with deadrise beta above 15°, allowing for pile-up: C_Pa >= 2.5 / (tan beta)^1.1.
 * Reference: DNV-RP-C205 (2019) §8 (slamming pressures). Between 4° and 15° the standard gives no
 * expression and this function throws: use wedgeWaterEntry or a dedicated analysis.
 * @param {{shape:string, deadriseDeg?:number}} o shape 'cylinder' | 'flat' | 'wedge'
 * @returns {number} C_Pa, -
 */
export function slammingPressureCoefficient({ shape, deadriseDeg }) {
  if (shape === 'cylinder') return 5.15;
  if (shape === 'flat') {
    if (deadriseDeg !== undefined && !(deadriseDeg >= 0 && deadriseDeg < 4)) throw new RangeError('deadriseDeg must be < 4 for flat-bottom slamming');
    return 2 * Math.PI;
  }
  if (shape === 'wedge') {
    if (!(deadriseDeg >= 15 && deadriseDeg < 90)) throw new RangeError('deadriseDeg must be within [15, 90) for the wedge expression');
    return 2.5 / Math.tan(deadriseDeg * DEG) ** 1.1;
  }
  throw new RangeError("shape must be 'cylinder', 'flat' or 'wedge'");
}

/**
 * Slamming coefficient of a horizontal circular cylinder versus submergence s (Campbell & Weynberg):
 *   C_s = 5.15 [ D/(D + 19 s) + 0.107 s/D ],  0 <= s <= D,
 * for use in F = ½ rho C_s D v² per unit length. C_s(0) = 5.15; pi is the von Kármán value at first contact.
 * Reference: Campbell & Weynberg (1980) Measurement of parameters affecting slamming, Southampton Univ. Rep.
 * 440; DNV-RP-C205 (2019) §8 (slamming on horizontal cylinders).
 * @param {{s:number, D:number}} o submergence m, diameter m
 * @returns {number} C_s, -
 */
export function cylinderSlamCoefficient({ s, D }) {
  needPos(D, 'D');
  if (!(s >= 0 && s <= D)) throw new RangeError('s must be within 0..D');
  return 5.15 * (D / (D + 19 * s) + (0.107 * s) / D);
}

/**
 * Probability and frequency of bottom / bow-flare slamming from the spectral moments of the relative vertical
 * motion and velocity at the section (Ochi's criterion: the section emerges and re-enters with a relative
 * velocity above a threshold):
 *   P = exp( -( d²/(2 m0r) + v_cr²/(2 m0v) ) ),   N per second = (1/2 pi) sqrt(m0v/m0r) · P.
 * The threshold defaults to Ochi's v_cr = 0.093 sqrt(g L) (12 ft/s for a 520 ft ship, Froude-scaled).
 * For bow-flare slamming use d = height of the flare above the still waterline.
 * Reference: Ochi (1964) Prediction of occurrence and severity of ship slamming at sea, 5th Symp. Naval
 * Hydrodynamics; Ochi & Motter (1973) SNAME Trans. 81. Narrow-band Gaussian relative motion.
 * @param {object} o
 * @param {number} o.m0r variance of relative motion, m²
 * @param {number} o.m0v variance of relative velocity, m²/s²
 * @param {number} o.d local draft (bottom slamming) or flare height (bow flare), m
 * @param {number} [o.L] ship length for the default threshold, m
 * @param {number} [o.vThreshold] threshold velocity, m/s (overrides the default)
 * @param {number} [o.g=9.80665]
 * @returns {{probability:number, ratePerHour:number, vThreshold:number, emergenceProbability:number}}
 */
export function slammingProbability({ m0r, m0v, d, L, vThreshold, g = G }) {
  needPos(m0r, 'm0r'); needPos(m0v, 'm0v');
  if (!(d >= 0)) throw new RangeError('d must be >= 0');
  let vc = vThreshold;
  if (vc === undefined) { needPos(L, 'L'); vc = 0.093 * Math.sqrt(g * L); }
  const emergence = Math.exp((-d * d) / (2 * m0r));
  const probability = emergence * Math.exp((-vc * vc) / (2 * m0v));
  return { probability, ratePerHour: (3600 / (2 * Math.PI)) * Math.sqrt(m0v / m0r) * probability, vThreshold: vc, emergenceProbability: emergence };
}

/**
 * Probability and frequency of deck wetness (green water): the relative motion exceeds the freeboard,
 *   P = exp(-f²/(2 m0r)),  N per second = (1/2 pi) sqrt(m2r/m0r) · P.
 * Reference: Ochi (1964), as slammingProbability. Rayleigh-distributed maxima of the relative motion; the
 * dynamic swell-up of the bow wave is not included (increase the relative motion or reduce the freeboard).
 * @param {{m0r:number, m2r?:number, freeboard:number}} o m², m²/s² (variance of relative velocity), m
 * @returns {{probability:number, ratePerHour:number|null}}
 */
export function deckWetnessProbability({ m0r, m2r, freeboard }) {
  needPos(m0r, 'm0r');
  if (!(freeboard >= 0)) throw new RangeError('freeboard must be >= 0');
  const probability = Math.exp((-freeboard * freeboard) / (2 * m0r));
  return { probability, ratePerHour: m2r > 0 ? (3600 / (2 * Math.PI)) * Math.sqrt(m2r / m0r) * probability : null };
}

/**
 * Green-water flow on deck by the dam-break analogy (Ritter's solution of the shallow-water equations for the
 * sudden release of a reservoir of depth h0 over a dry horizontal bed):
 *   h(x, t) = (2 sqrt(g h0) - x/t)² / (9 g),  u(x, t) = (2/3)(sqrt(g h0) + x/t),  -sqrt(g h0) <= x/t <= 2 sqrt(g h0);
 *   front speed 2 sqrt(g h0); at the dam (deck edge) h = (4/9) h0 and u = (2/3) sqrt(g h0).
 * h0 is the height of the water above the deck at the bulwark/deck edge (freeboard exceedance).
 * Loads: stagnation pressure ½ rho u², momentum flux per unit width rho h u² (force on a wall that stops the
 * flow), and hydrostatic deck pressure rho g h.
 * Reference: Ritter (1892) Die Fortpflanzung der Wasserwellen, Z. VDI 36; Buchner (2002) Green water on
 * ship-type offshore structures (dam-break analogy). Frictionless, dry horizontal deck, no ship motion;
 * impact peaks on structures can exceed the quasi-steady values.
 * @param {{h0:number, x:number, t:number, rho?:number, g?:number}} o height above deck m, distance from the deck
 *   edge m (>= 0), time after release s (> 0)
 * @returns {{h:number, u:number, frontSpeed:number, frontPosition:number, arrived:boolean, stagnationPressure:number,
 *   momentumFlux:number, hydrostaticPressure:number}} m, m/s, m/s, m, -, Pa, N/m, Pa
 */
export function damBreakDeckFlow({ h0, x, t, rho = 1025, g = G }) {
  needPos(h0, 'h0'); needPos(t, 't');
  if (!(x >= 0)) throw new RangeError('x must be >= 0');
  const c0 = Math.sqrt(g * h0);
  const xi = x / t;
  const arrived = xi <= 2 * c0;
  const h = arrived ? (2 * c0 - xi) ** 2 / (9 * g) : 0;
  const u = arrived ? (2 / 3) * (c0 + xi) : 0;
  return {
    h, u, frontSpeed: 2 * c0, frontPosition: 2 * c0 * t, arrived, stagnationPressure: 0.5 * rho * u * u,
    momentumFlux: rho * h * u * u, hydrostaticPressure: rho * g * h,
  };
}

// Zeros of J1'(x): radial wave numbers of the antisymmetric (m = 1) sloshing modes of an upright cylinder.
const J1P_ZEROS = [1.841183781, 5.331442774, 8.536316366];

/**
 * Natural sloshing periods by linear potential theory.
 *   Rectangular tank of length L (in the direction of motion), fill depth h:
 *     omega_n² = g (n pi / L) tanh(n pi h / L), n = 1, 2, 3
 *   Upright circular cylindrical tank of radius R: omega_n² = g (e_n / R) tanh(e_n h / R),
 *     e_n = 1.8412, 5.3314, 8.5363 (zeros of J1').
 * Reference: Lamb (1932) Hydrodynamics; Faltinsen & Timokha (2009) Sloshing; Abramson (1966) NASA SP-106.
 * Linear, inviscid, rigid smooth tank without internal structure.
 * @param {{shape:string, L?:number, R?:number, h:number, g?:number}} o shape 'rectangular' | 'cylindrical'; m
 * @returns {{periods:number[], omegas:number[]}} first three natural periods s and frequencies rad/s
 */
export function sloshingPeriods({ shape, L, R, h, g = G }) {
  needPos(h, 'h');
  let ks;
  if (shape === 'rectangular') { needPos(L, 'L'); ks = [1, 2, 3].map((n) => (n * Math.PI) / L); } else if (shape === 'cylindrical') { needPos(R, 'R'); ks = J1P_ZEROS.map((e) => e / R); } else throw new RangeError("shape must be 'rectangular' or 'cylindrical'");
  const omegas = ks.map((k) => Math.sqrt(g * k * Math.tanh(k * h)));
  return { periods: omegas.map((w) => (2 * Math.PI) / w), omegas };
}

/**
 * Sloshing resonance screening: compares the first sloshing period with excitation (vessel motion) periods
 * and flags those within a relative band.
 * Reference: resonance condition of linear sloshing theory (see sloshingPeriods). The band is a screening
 * input, not a code value (default ±20 %); class rules define their own separation requirements.
 * @param {{sloshingPeriod:number, motionPeriods:number[], band?:number}} o s, s, -
 * @returns {{ratios:number[], resonant:boolean[], any:boolean}} T_motion / T_slosh
 */
export function sloshingResonanceCheck({ sloshingPeriod, motionPeriods, band = 0.2 }) {
  needPos(sloshingPeriod, 'sloshingPeriod');
  if (!(band > 0 && band < 1)) throw new RangeError('band must be within (0, 1)');
  const ratios = motionPeriods.map((T) => T / sloshingPeriod);
  const resonant = ratios.map((r) => Math.abs(r - 1) <= band);
  return { ratios, resonant, any: resonant.some(Boolean) };
}

/**
 * Horizontal wave-in-deck force by the silhouette method:
 *   F = ½ rho Cd (alpha_wkf u + alpha_cbf U)² A,  A = width × (crest elevation - underside of deck),
 * applied at the mid-height of the wetted silhouette.
 * Drag coefficients of the method (inputs here): heavily equipped (solid) deck 2.5 end-on/broadside and 1.9
 * diagonal; moderately equipped 2.0 and 1.5; bare (no equipment) 1.6 and 1.2.
 * Reference: API RP 2A-WSD 21st ed. §C17.6.2 (wave/current deck force); API RP 2SIM (2014).
 * u is the crest-elevation fluid velocity from the wave theory in use; the wave kinematics factor
 * (API RP 2A: 0.85–0.95 for tropical storms, 0.95–1.00 for extra-tropical storms) and the current blockage factor are inputs.
 * Vertical (uplift) forces and impact dynamics are not covered.
 * @param {object} o
 * @param {number} o.crestElevation crest elevation, m (same datum as deckBottom)
 * @param {number} o.deckBottom elevation of the underside of the lowest deck beam, m
 * @param {number} [o.deckTop] top of the solid deck silhouette, m (limits the wetted height)
 * @param {number} o.width deck width normal to the wave, m
 * @param {number} o.u wave-induced horizontal velocity at the crest, m/s
 * @param {number} [o.U=0] in-line current speed, m/s
 * @param {number} o.Cd silhouette drag coefficient, -
 * @param {number} [o.waveKinematicsFactor=1]
 * @param {number} [o.currentBlockageFactor=1]
 * @param {number} [o.rho=1025]
 * @returns {{force:number, area:number, elevation:number, inundation:number}} N, m², m (point of application), m
 */
export function waveInDeckForce({ crestElevation, deckBottom, deckTop, width, u, U = 0, Cd, waveKinematicsFactor = 1, currentBlockageFactor = 1, rho = 1025 }) {
  needPos(width, 'width'); needPos(Cd, 'Cd');
  if (!Number.isFinite(crestElevation) || !Number.isFinite(deckBottom)) throw new RangeError('crestElevation and deckBottom must be finite');
  const top = deckTop === undefined ? crestElevation : Math.min(crestElevation, deckTop);
  const inundation = Math.max(top - deckBottom, 0);
  const area = width * inundation;
  const v = waveKinematicsFactor * u + currentBlockageFactor * U;
  return { force: 0.5 * rho * Cd * v * v * area, area, elevation: deckBottom + 0.5 * inundation, inundation };
}

/**
 * Vortex-induced-motion screening of a deep-draft floater (spar, deep-draft semi column):
 *   reduced velocity Ur = U Tn / D; lock-in when Ur lies within a band; resonance of the shedding frequency
 *   (St U / D) with the natural frequency occurs at Ur = 1/St.
 * The lock-in band and the response amplitude ratio A/D are inputs (model-test data govern); the default band
 * 4–10 is a screening range only, commonly quoted for transverse VIM of bare circular hulls.
 * Reference: DNV-RP-C205 (2019) §9 (reduced velocity, Strouhal relation); DNV-OS-E301 (2021) ch. 2 (VIM to
 * be included in mooring analysis). No prediction of amplitude is made here.
 * @param {object} o
 * @param {number} o.U current speed, m/s
 * @param {number} o.Tn natural period of the transverse (sway) or in-line (surge) mode, s
 * @param {number} o.D hull / column diameter, m
 * @param {number[]} [o.lockInRange=[4,10]] reduced-velocity band, -
 * @param {number} [o.AoverD] amplitude ratio inside lock-in (input), -
 * @param {number} [o.St=0.2] Strouhal number, -
 * @returns {{Ur:number, lockIn:boolean, sheddingPeriod:number, UrResonance:number, amplitude:number|null,
 *   criticalSpeeds:number[]}} -, -, s, -, motion amplitude m (when AoverD given and locked in, else 0 / null),
 *   current speeds at the band limits m/s
 */
export function vimScreening({ U, Tn, D, lockInRange = [4, 10], AoverD, St = 0.2 }) {
  if (!(U >= 0)) throw new RangeError('U must be >= 0');
  needPos(Tn, 'Tn'); needPos(D, 'D'); needPos(St, 'St');
  if (!(lockInRange[1] > lockInRange[0])) throw new RangeError('lockInRange must be [low, high]');
  const Ur = (U * Tn) / D;
  const lockIn = Ur >= lockInRange[0] && Ur <= lockInRange[1];
  return {
    Ur, lockIn, sheddingPeriod: U > 0 ? D / (St * U) : Infinity, UrResonance: 1 / St,
    amplitude: AoverD === undefined ? null : lockIn ? AoverD * D : 0,
    criticalSpeeds: [(lockInRange[0] * D) / Tn, (lockInRange[1] * D) / Tn],
  };
}
