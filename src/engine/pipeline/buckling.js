/**
 * Global buckling of pipelines: Hobbs lateral and vertical thermal-buckling modes, critical force and
 * temperature, infinite-mode screening force, upheaval-buckling download (Palmer et al.), and buckle
 * spacing / feed-in from axial friction.
 *
 * Sign convention: axial forces in this file are compressive-positive (P, P0). All quantities SI.
 */

import { requiredCover } from './pipesoil.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

const FOUR_PI2 = 4 * Math.PI * Math.PI;

/**
 * Constants of the Hobbs lateral buckling modes: P = k1 EI/L²,
 * P0 = P + k3 φ w L [ (1 + k2 E A φ w L⁵/(EI)²)^½ − 1 ], ŷ = k4 φ w L⁴/EI, M̂ = k5 φ w L².
 * Mode ∞: P0 = P + k2 E A (φ w/EI)² L⁶.
 *
 * Reference: Hobbs, R.E. (1984), "In-service buckling of heated pipelines", J. Transp. Eng. ASCE 110(2),
 * Table 1. The mode-1 and mode-∞ shape constants (k1, k4, k5, and k2 of mode ∞) are re-derived from the
 * beam-column solution in the test suite.
 */
export const HOBBS_CONSTANTS = Object.freeze({
  1: Object.freeze({ k1: 80.76, k2: 6.391e-5, k3: 0.5, k4: 2.407e-3, k5: 0.06938 }),
  2: Object.freeze({ k1: FOUR_PI2, k2: 1.743e-4, k3: 1.0, k4: 5.532e-3, k5: 0.1088 }),
  3: Object.freeze({ k1: 34.06, k2: 1.668e-4, k3: 1.294, k4: 1.032e-2, k5: 0.1434 }),
  4: Object.freeze({ k1: 28.20, k2: 2.144e-4, k3: 1.608, k4: 1.047e-2, k5: 0.1483 }),
  inf: Object.freeze({ k1: FOUR_PI2, k2: 4.7050e-5, k3: null, k4: 4.4495e-3, k5: 0.05066 }),
});

function modeKey(mode) {
  const k = mode === Infinity || mode === 'inf' || mode === 'infinity' || mode === '∞' ? 'inf' : String(mode);
  if (!HOBBS_CONSTANTS[k]) throw new RangeError('mode must be 1, 2, 3, 4 or "inf"');
  return k;
}

/**
 * Hobbs lateral buckling mode: buckle force, pre-buckle (far-field) force, amplitude, maximum moment and
 * bending strain for a buckle length L.
 *
 * With distinct axial and lateral friction coefficients the compatibility term becomes
 * k2 E A φ_L² w L⁵/(φ_A (EI)²) and the slip term k3 φ_A w L (identical to Hobbs when φ_A = φ_L).
 *
 * Reference: Hobbs (1984), J. Transp. Eng. 110(2), Eq. 8–14 and Table 1.
 * Validity: straight, perfectly elastic pipe on a rigid frictional seabed, small slopes, fully mobilised
 * lateral friction; no initial imperfection (the curves give post-buckling equilibrium, not snap loads).
 *
 * @param {object} o
 * @param {1|2|3|4|'inf'} o.mode
 * @param {number} o.L buckle length (of the central lobe; half-wave for mode ∞) [m]
 * @param {number} o.EI bending stiffness [N·m²]
 * @param {number} o.EA axial stiffness [N]
 * @param {number} o.w submerged weight [N/m]
 * @param {number} o.muLateral lateral friction coefficient φ_L [-]
 * @param {number} [o.muAxial=o.muLateral] axial friction coefficient φ_A [-]
 * @param {number} [o.D] outer diameter [m] (for the bending strain)
 * @returns {{P:number, P0:number, amplitude:number, moment:number, strain:(number|null), slipLength:(number|null)}}
 *   P force in the buckle [N]; P0 force away from the buckle [N]; amplitude [m]; moment [N·m]; strain [-];
 *   slipLength axial slip length each side [m] (null for mode ∞)
 */
export function hobbsLateral({ mode, L, EI, EA, w, muLateral, muAxial, D }) {
  const c = HOBBS_CONSTANTS[modeKey(mode)];
  pos('L', L); pos('EI', EI); pos('EA', EA); pos('w', w); pos('muLateral', muLateral);
  const fa = muAxial === undefined ? muLateral : pos('muAxial', muAxial);
  const P = (c.k1 * EI) / (L * L);
  let P0;
  let slip = null;
  if (c.k3 === null) P0 = P + c.k2 * EA * ((muLateral * w) / EI) ** 2 * L ** 6;
  else {
    P0 = P + c.k3 * fa * w * L * (Math.sqrt(1 + (c.k2 * EA * muLateral * muLateral * w * L ** 5) / (fa * EI * EI)) - 1);
    slip = (P0 - P) / (fa * w);
  }
  const moment = c.k5 * muLateral * w * L * L;
  return { P, P0, amplitude: (c.k4 * muLateral * w * L ** 4) / EI, moment, strain: D === undefined ? null : (moment * pos('D', D)) / (2 * EI), slipLength: slip };
}

/**
 * Hobbs vertical (upheaval) mode on a rigid seabed: P = 80.76 EI/L²,
 * P0 = P + (w L/EI) [ 1.597e-5 E A φ w L⁵ − 0.25 (φ EI)² ]^½, ŷ = 2.408e-3 w L⁴/EI, M̂ = 0.06938 w L²,
 * maximum slope 8.657e-3 w L³/EI.
 *
 * Reference: Hobbs (1984), J. Transp. Eng. 110(2), Eq. 1–7.
 * Validity: as hobbsLateral; requires the bracket to be non-negative (slip length ≥ 0), otherwise throws.
 *
 * @param {object} o
 * @param {number} o.L lifted length [m]
 * @param {number} o.EI [N·m²]
 * @param {number} o.EA [N]
 * @param {number} o.w submerged weight (plus any cover download) [N/m]
 * @param {number} o.muAxial axial friction coefficient φ [-]
 * @param {number} [o.D] outer diameter [m]
 * @returns {{P:number, P0:number, amplitude:number, moment:number, slope:number, strain:(number|null)}}
 */
export function hobbsVertical({ L, EI, EA, w, muAxial, D }) {
  pos('L', L); pos('EI', EI); pos('EA', EA); pos('w', w); pos('muAxial', muAxial);
  const P = (80.76 * EI) / (L * L);
  const br = 1.597e-5 * EA * muAxial * w * L ** 5 - 0.25 * (muAxial * EI) ** 2;
  if (br < 0) throw new RangeError('L is too short for the vertical mode solution (negative slip length)');
  const moment = 0.06938 * w * L * L;
  return { P, P0: P + ((w * L) / EI) * Math.sqrt(br), amplitude: (2.408e-3 * w * L ** 4) / EI, moment, slope: (8.657e-3 * w * L ** 3) / EI,
    strain: D === undefined ? null : (moment * pos('D', D)) / (2 * EI) };
}

function minimise(f, a, b) {
  // coarse log scan, then golden-section refinement in the bracketing interval
  const N = 400;
  let best = Infinity;
  let ib = 0;
  const xs = new Float64Array(N + 1);
  for (let i = 0; i <= N; i++) {
    xs[i] = a * (b / a) ** (i / N);
    let v;
    try { v = f(xs[i]); } catch { v = Infinity; }
    if (v < best) { best = v; ib = i; }
  }
  let lo = xs[Math.max(0, ib - 1)];
  let hi = xs[Math.min(N, ib + 1)];
  const g = (Math.sqrt(5) - 1) / 2;
  const safe = (x) => { try { return f(x); } catch { return Infinity; } };
  let c = hi - g * (hi - lo);
  let d = lo + g * (hi - lo);
  let fc = safe(c);
  let fd = safe(d);
  for (let i = 0; i < 120; i++) {
    if (fc < fd) { hi = d; d = c; fd = fc; c = hi - g * (hi - lo); fc = safe(c); }
    else { lo = c; c = d; fc = fd; d = lo + g * (hi - lo); fd = safe(d); }
  }
  return 0.5 * (lo + hi);
}

/**
 * Minimum ("safe") pre-buckle force of a Hobbs mode — the lowest P0 on the P0–L curve — with the
 * corresponding buckle length, amplitude, moment and the temperature rise that produces it in a fully
 * restrained line: ΔT = [P0 + H − Δp_i A_i (1 − 2ν)]/(E A α).
 *
 * Reference: Hobbs (1984); restrained force per DNV-RP-F110 (2018) Eq. 6.
 * Validity: as hobbsLateral / hobbsVertical. Below this force no buckled equilibrium of that mode exists.
 *
 * @param {object} o arguments of hobbsLateral without L, or with mode 'vertical' those of hobbsVertical
 * @param {1|2|3|4|'inf'|'vertical'} o.mode
 * @param {number} [o.alpha=11.7e-6] thermal expansion coefficient [1/K]
 * @param {number} [o.deltaPi=0] internal pressure rise [Pa]
 * @param {number} [o.Ai=0] bore area [m²]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.H=0] residual lay tension [N]
 * @param {number} [o.Lmin=1] search range [m]
 * @param {number} [o.Lmax=1e4] search range [m]
 * @returns {{L:number, P0min:number, P:number, amplitude:number, moment:number, strain:(number|null), deltaT:number}}
 */
export function hobbsMinimum(o) {
  const { mode, alpha = 11.7e-6, deltaPi = 0, Ai = 0, nu = 0.3, H = 0, Lmin = 1, Lmax = 1e4 } = o;
  const f = mode === 'vertical' ? (L) => hobbsVertical({ ...o, L }) : (L) => hobbsLateral({ ...o, L });
  const L = minimise((x) => f(x).P0, Lmin, Lmax);
  const r = f(L);
  return { L, P0min: r.P0, P: r.P, amplitude: r.amplitude, moment: r.moment, strain: r.strain,
    deltaT: (r.P0 + H - deltaPi * Ai * (1 - 2 * nu)) / (o.EA * alpha) };
}

/**
 * P0–L and P–L curves of a Hobbs mode for plotting.
 *
 * Reference: Hobbs (1984).
 *
 * @param {object} o arguments of hobbsLateral / hobbsVertical without L
 * @param {number} o.Lmin [m]
 * @param {number} o.Lmax [m]
 * @param {number} [o.n=100] number of points [-]
 * @returns {{L:Float64Array, P0:Float64Array, P:Float64Array, amplitude:Float64Array}} forces [N]; NaN where no solution
 */
export function hobbsCurve(o) {
  const { Lmin, Lmax, n = 100, mode } = o;
  pos('Lmin', Lmin); pos('Lmax', Lmax);
  const L = new Float64Array(n);
  const P0 = new Float64Array(n);
  const P = new Float64Array(n);
  const amplitude = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    L[i] = Lmin + ((Lmax - Lmin) * i) / (n - 1);
    try {
      const r = mode === 'vertical' ? hobbsVertical({ ...o, L: L[i] }) : hobbsLateral({ ...o, L: L[i] });
      P0[i] = r.P0; P[i] = r.P; amplitude[i] = r.amplitude;
    } catch { P0[i] = NaN; P[i] = NaN; amplitude[i] = NaN; }
  }
  return { L, P0, P, amplitude };
}

/**
 * Infinite-mode lateral buckling force in closed form (the minimum of the Hobbs mode-∞ curve):
 * S_∞ = c (EA)^¼ (EI)^¼ f_L^½ with c = (16π²/3)(3 k2/(4π²))^¼ ≈ 2.29, f_L = lateral soil resistance per
 * unit length, at half-wavelength L = [4π² (EI)³/(3 k2 EA f_L²)]^(1/8).
 * Used as a screening value: a line whose restrained compressive force stays below S_∞ cannot sustain a
 * lateral buckle on an even seabed.
 *
 * Reference: Hobbs (1984) mode ∞, minimised analytically; the same expression is the even-seabed
 * screening force of DNV-RP-F110 (2018) §6 (given there with the rounded constant).
 *
 * @param {object} o
 * @param {number} o.EA [N]
 * @param {number} o.EI [N·m²]
 * @param {number} o.fL lateral soil resistance [N/m]
 * @param {number} [o.S0] restrained compressive force to screen [N]
 * @returns {{Sinf:number, L:number, constant:number, susceptible:(boolean|null), uc:(number|null)}} Sinf [N]; L [m]
 */
export function infiniteModeForce({ EA, EI, fL, S0 }) {
  pos('EA', EA); pos('EI', EI); pos('fL', fL);
  const k2 = HOBBS_CONSTANTS.inf.k2;
  const constant = ((4 * FOUR_PI2) / 3) * ((3 * k2) / FOUR_PI2) ** 0.25;
  const Sinf = constant * EA ** 0.25 * EI ** 0.25 * Math.sqrt(fL);
  const L = ((FOUR_PI2 * EI ** 3) / (3 * k2 * EA * fL * fL)) ** 0.125;
  return { Sinf, L, constant, susceptible: S0 === undefined ? null : S0 > Sinf, uc: S0 === undefined ? null : S0 / Sinf };
}

/**
 * Dimensionless download required to hold a pipe down over an imperfection (universal design curve):
 * Φ_w = 0.0646 for Φ_L < 4.49; 5.68/Φ_L² − 88.35/Φ_L⁴ for 4.49 ≤ Φ_L ≤ 8.06; 9.6/Φ_L² − 343/Φ_L⁴ above,
 * with Φ_w = w EI/(δ P²) and Φ_L = L √(P/EI).
 *
 * Reference: Palmer, Ellinas, Richards & Guijt (1990), "Design of submarine pipelines against upheaval
 * buckling", OTC 6335.
 *
 * @param {number} phiL dimensionless imperfection length [-]
 * @returns {number} Φ_w [-]
 */
export function palmerPhiW(phiL) {
  pos('phiL', phiL);
  if (phiL < 4.49) return 0.0646;
  if (phiL <= 8.06) return 5.68 / phiL ** 2 - 88.35 / phiL ** 4;
  return 9.6 / phiL ** 2 - 343 / phiL ** 4;
}

/**
 * Length of a prop imperfection: a pipe of submerged weight w draped over a point support of height δ on a
 * rigid seabed touches down a distance ℓ = (72 EI δ/w)^¼ either side, so L = 2ℓ.
 *
 * Reference: beam statics (zero moment and slope at touchdown, zero slope at the prop); used with the
 * installation (empty or flooded) weight, without axial force.
 *
 * @param {object} o
 * @param {number} o.EI [N·m²]
 * @param {number} o.delta prop height [m]
 * @param {number} o.w submerged weight when the profile formed [N/m]
 * @returns {number} total imperfection length L [m]
 */
export function propImperfectionLength({ EI, delta, w }) {
  pos('EI', EI); pos('delta', delta); pos('w', w);
  return 2 * ((72 * EI * delta) / w) ** 0.25;
}

/**
 * Upheaval buckling: required total download, the part to be supplied by cover, and the cover height.
 * w_req = γ · Φ_w(Φ_L) δ P²/EI; cover resistance R = max(0, w_req − w_pipe); cover from requiredCover
 * (Schaminée et al. uplift model for sand / rock, vertical-slip model for clay).
 *
 * Reference: Palmer et al. (1990), OTC 6335; uplift per Schaminée et al. (1990), OTC 6486.
 * The safety factor γ on download is project-specific (Palmer et al. discuss the sensitivity to δ; DNV-RP-F110
 * uses partial factors on load and resistance instead) — argument.
 * Validity: small-slope beam-column, P the effective compressive force in operation.
 *
 * @param {object} o
 * @param {number} o.P effective compressive force [N]
 * @param {number} o.EI [N·m²]
 * @param {number} o.delta imperfection height [m]
 * @param {number} [o.L] imperfection length [m] (default: prop length from wInstall)
 * @param {number} [o.wInstall] submerged weight when laid [N/m] (needed when L is not given)
 * @param {number} o.wPipe submerged weight in operation [N/m]
 * @param {number} [o.safetyFactor=1] factor on the required download [-]
 * @param {{soil:'sand'|'clay', D:number, gammaSub?:number, f?:number, su?:number, Nc?:number}} [o.cover] uplift model data
 * @returns {{phiL:number, phiW:number, L:number, wRequired:number, coverResistance:number, coverHeight:(number|null), stableWithoutCover:boolean}}
 *   wRequired, coverResistance [N/m]; coverHeight above the crown [m]
 */
export function upheavalDownload({ P, EI, delta, L, wInstall, wPipe, safetyFactor = 1, cover }) {
  pos('P', P); pos('EI', EI); pos('delta', delta); nonNeg('wPipe', wPipe); pos('safetyFactor', safetyFactor);
  const len = L === undefined ? propImperfectionLength({ EI, delta, w: wInstall }) : pos('L', L);
  const phiL = len * Math.sqrt(P / EI);
  const phiW = palmerPhiW(phiL);
  const wRequired = (safetyFactor * phiW * delta * P * P) / EI;
  const coverResistance = Math.max(0, wRequired - wPipe);
  return { phiL, phiW, L: len, wRequired, coverResistance,
    coverHeight: cover ? requiredCover({ ...cover, R: coverResistance }) : null, stableWithoutCover: coverResistance === 0 };
}

/**
 * Buckle interaction along a line with planned triggers (snake lay, sleepers, buoyancy): slip length,
 * feed-in to each buckle and admissible trigger spacing from the axial friction build-up.
 *
 * After a buckle forms the force at its crown drops from P0 to Ppost and recovers at f = μ_a w per metre.
 *   isolated buckle: slip length each side (P0 − Ppost)/f, feed-in (P0 − Ppost)²/(E A f);
 *   buckles at spacing s < 2(P0 − Ppost)/f: virtual anchors midway, feed-in [(P0 − Ppost) s − f s²/4]/(E A),
 *     force at the virtual anchor Ppost + f s/2.
 *   Spacing limits: s ≤ 2 (Pcrit,unplanned − Ppost)/f keeps the force between triggers below the capacity
 *     of the untriggered pipe; s ≥ (Pcrit,trigger − Ppost)/f lets the force rebuild enough to fire the
 *     adjacent trigger.
 *
 * Reference: axial equilibrium with Coulomb friction (virtual anchor spacing concept of DNV-RP-F110 §7 and
 * the SAFEBUCK guideline); first-principles, no empirical coefficients.
 *
 * @param {object} o
 * @param {number} o.P0 fully restrained compressive force [N]
 * @param {number} o.Ppost compressive force remaining in the buckle [N]
 * @param {number} o.f axial friction [N/m]
 * @param {number} o.EA [N]
 * @param {number} [o.spacing] trigger spacing [m]
 * @param {number} [o.PcritUnplanned] critical buckling force of the pipe between triggers [N]
 * @param {number} [o.PcritTrigger] critical buckling force at a trigger [N]
 * @returns {{slipLength:number, feedInIsolated:number, feedIn:number, forceAtAnchor:number, interacting:boolean,
 *   maxSpacing:(number|null), minSpacing:(number|null)}} lengths [m]; forces [N]
 */
export function buckleSpacing({ P0, Ppost, f, EA, spacing, PcritUnplanned, PcritTrigger }) {
  pos('P0', P0); nonNeg('Ppost', Ppost); pos('f', f); pos('EA', EA);
  if (Ppost >= P0) throw new RangeError('Ppost must be below P0');
  const dP = P0 - Ppost;
  const slipLength = dP / f;
  const feedInIsolated = (dP * dP) / (EA * f);
  let feedIn = feedInIsolated;
  let forceAtAnchor = P0;
  let interacting = false;
  if (spacing !== undefined && pos('spacing', spacing) < 2 * slipLength) {
    interacting = true;
    feedIn = (dP * spacing - (f * spacing * spacing) / 4) / EA;
    forceAtAnchor = Ppost + (f * spacing) / 2;
  }
  return { slipLength, feedInIsolated, feedIn, forceAtAnchor, interacting,
    maxSpacing: PcritUnplanned === undefined ? null : Math.max(0, (2 * (PcritUnplanned - Ppost)) / f),
    minSpacing: PcritTrigger === undefined ? null : Math.max(0, (PcritTrigger - Ppost) / f) };
}
