/**
 * Vortex-induced vibration of slender members: shedding relations, modal analysis of tensioned beams, modal
 * screening in sheared current (power-in length, damping, empirical response amplitude, drag amplification),
 * curvature → stress → fatigue, suppression coverage, and the Facchinetti–de Langre–Biolley wake-oscillator
 * model of an elastically mounted rigid cylinder. Vortex-induced motion of floaters is handled elsewhere.
 *
 * Mode sets are exchanged in the shape `{ s: Float64Array [m], modes: [{ fHz, phi, curvature? }] }` (phi: transverse
 * shape at s, curvature: d²phi/ds² per unit phi [1/m²]) — produced by tensionedBeamModes here and by rodModes in
 * rod.js.
 */

import { eigSym } from '../math/eigen.js';
import { rk4 } from '../math/ode.js';
import { bisect } from '../math/roots.js';

const SEC_PER_YEAR = 365.25 * 24 * 3600;

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`viv: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}
function field(v, s, name) {
  const n = s.length, out = new Float64Array(n);
  if (typeof v === 'function') for (let i = 0; i < n; i++) out[i] = v(s[i]);
  else if (typeof v === 'number') out.fill(v);
  else if (v && v.length === n) out.set(v);
  else throw new RangeError(`viv: argument ${name} must be a number, a function of s or an array matching s`);
  return out;
}
function trap(f, s, w) {
  let a = 0;
  for (let i = 1; i < s.length; i++) a += 0.5 * (f[i] * (w ? w[i] : 1) + f[i - 1] * (w ? w[i - 1] : 1)) * (s[i] - s[i - 1]);
  return a;
}

/* ------------------------------------------------------------------------------------------------------- */
/* Shedding relations                                                                                        */

/**
 * Strouhal number of a smooth circular cylinder versus Reynolds number.
 *
 * Curve: St = 0.212·(1 − 21.2/Re) for 50 ≤ Re < 150 (Roshko 1954, NACA Rep. 1191, laminar range);
 * St = 0.198·(1 − 19.7/Re) for 250 ≤ Re ≤ 2×10⁵ (Roshko's subcritical fit as quoted by Blevins,
 * "Flow-Induced Vibration", §3.2); linear blend between the two for 150 ≤ Re < 250. Above 2×10⁵ (critical and
 * supercritical regimes) the Strouhal number depends strongly on roughness and turbulence (DNV-RP-C205 Fig. 9-1)
 * and the caller's `stHighRe` is returned.
 * Validity: Re ≥ 50 (no regular shedding below about Re = 47, returns 0), stationary smooth cylinder.
 * @param {number} Re Reynolds number U·D/ν [-]
 * @param {object} [o] `{ stHighRe? = 0.2 [-] value used above Re = 2×10⁵ }`
 * @returns {number} Strouhal number [-]
 */
export function strouhalNumber(Re, { stHighRe = 0.2 } = {}) {
  num(Re, 'Re', { min: 0 });
  const lo = (r) => 0.212 * (1 - 21.2 / r), hi = (r) => 0.198 * (1 - 19.7 / r);
  if (Re < 50) return 0;
  if (Re < 150) return lo(Re);
  if (Re < 250) { const f = (Re - 150) / 100; return (1 - f) * lo(Re) + f * hi(Re); }
  if (Re <= 2e5) return hi(Re);
  return stHighRe;
}

/**
 * Vortex shedding frequency fs = St·U/D (DNV-RP-C205 §9.1.1).
 * @param {object} o `{ U [m/s], D [m], St? = 0.2 [-] }` @returns {number} [Hz]
 */
export function sheddingFrequency({ U, D, St = 0.2 }) {
  num(U, 'U', { min: 0 }); num(D, 'D', { min: 0, strictMin: true }); num(St, 'St', { min: 0, strictMin: true });
  return St * U / D;
}

/**
 * Reduced velocity VR = U/(fn·D) (DNV-RP-C205 §9.1.3).
 * @param {object} o `{ U [m/s], fn [Hz], D [m] }` @returns {number} [-]
 */
export function reducedVelocity({ U, fn, D }) {
  num(U, 'U', { min: 0 }); num(fn, 'fn', { min: 0, strictMin: true }); num(D, 'D', { min: 0, strictMin: true });
  return U / (fn * D);
}

/**
 * Stability (mass-damping) parameters: KS = 2·me·δ/(ρ·D²) with δ = 2πζ (DNV-RP-C205 §9.1.6; identical to
 * 4π·me·ζ/(ρD²) of DNV-RP-F105) and the Skop–Griffin parameter SG = 2π·St²·KS.
 * @param {object} o `{ me [kg/m effective mass incl. added mass], zeta [-], rho? = 1025 [kg/m³], D [m], St? = 0.2 }`
 * @returns {{ Ks: number, SG: number, delta: number }} [-]
 */
export function stabilityParameter({ me, zeta, rho = 1025, D, St = 0.2 }) {
  num(me, 'me', { min: 0, strictMin: true }); num(zeta, 'zeta', { min: 0 }); num(D, 'D', { min: 0, strictMin: true });
  const delta = 2 * Math.PI * zeta, Ks = 2 * me * delta / (rho * D * D);
  return { Ks, SG: 2 * Math.PI * St * St * Ks, delta };
}

/**
 * Lock-in classification from the reduced velocity. Default windows are the DNV-RP-C205 §9 guidance for steady
 * current: in-line excitation 1.0 ≤ VR ≤ 4.5 (possible only for KS ≤ 1.8), cross-flow 3 ≤ VR ≤ 16 with the largest
 * cross-flow response typically at 5 ≤ VR ≤ 9. All limits are arguments.
 * @param {object} o `{ Vr [-], Ks? [-], inline? = [1.0, 4.5], crossflow? = [3, 16], crossflowPeak? = [5, 9],
 *   ksInlineMax? = 1.8 }`
 * @returns {{ inline: boolean, crossflow: boolean, crossflowPeak: boolean }}
 */
export function lockIn({ Vr, Ks, inline = [1.0, 4.5], crossflow = [3, 16], crossflowPeak = [5, 9], ksInlineMax = 1.8 }) {
  num(Vr, 'Vr', { min: 0 });
  const inR = (r) => Vr >= r[0] && Vr <= r[1];
  return { inline: inR(inline) && (Ks === undefined || Ks <= ksInlineMax), crossflow: inR(crossflow), crossflowPeak: inR(crossflowPeak) };
}

/**
 * In-line VIV onset reduced velocity, DNV-RP-F105 (2006) §4.3.5: VR,onset = 1.0/γon for Ksd < 0.4,
 * (0.6 + Ksd)/γon for 0.4 ≤ Ksd < 1.6, 2.2/γon for Ksd ≥ 1.6, with Ksd = Ks/γk.
 * @param {object} o `{ Ks [-], gammaOn? = 1.1, gammaK? = 1.0 }` (γk = 1.0 / 1.15 / 1.30 by safety class in F105)
 * @returns {number} onset reduced velocity [-]
 */
export function inlineOnsetF105({ Ks, gammaOn = 1.1, gammaK = 1.0 }) {
  num(Ks, 'Ks', { min: 0 }); num(gammaOn, 'gammaOn', { min: 0, strictMin: true }); num(gammaK, 'gammaK', { min: 0, strictMin: true });
  const k = Ks / gammaK;
  return (k < 0.4 ? 1.0 : k < 1.6 ? 0.6 + k : 2.2) / gammaOn;
}

/**
 * Keulegan–Carpenter number KC = Um·T/D for wave-induced (oscillatory-flow) shedding. In oscillatory flow the
 * shedding pattern and lift frequency depend on KC (Sumer & Fredsøe, "Hydrodynamics Around Cylindrical
 * Structures", Ch. 3; DNV-RP-C205 §9.4): the steady-current relations of this file apply only for large KC
 * (KC > ~40); a KC-dependent response model is not implemented.
 * @param {object} o `{ Um [m/s orbital velocity amplitude], T [s], D [m] }` @returns {number} KC [-]
 */
export function keuleganCarpenter({ Um, T, D }) {
  num(Um, 'Um', { min: 0 }); num(T, 'T', { min: 0, strictMin: true }); num(D, 'D', { min: 0, strictMin: true });
  return Um * T / D;
}

/**
 * Empirical maximum cross-flow response amplitude versus reduced damping.
 *
 * 'griffin' (Griffin & Ramberg 1982, as tabulated by Blevins, "Flow-Induced Vibration" 2nd ed., Table 3-2):
 *   A/D = 1.29·γ / [1 + 0.43·SG]^3.35
 * 'sarpkaya' (Sarpkaya 1979, same table): A/D = 0.32·γ / sqrt(0.06 + SG²)
 * with SG = 2π·St²·(2m·(2πζ)/(ρD²)) and mode factor γ = φmax·sqrt(∫φ²ds/∫φ⁴ds) (1.0 rigid cylinder, 1.155 sine
 * mode). The coefficients are overridable; guidance values from DNV-RP-C205 §9 / DNV-RP-F204 may be imposed with
 * `override`.
 * @param {object} o `{ SG [-], gamma? = 1, model? = 'griffin' | 'sarpkaya', coefficients?: number[], override?: number }`
 *   coefficients = [1.29, 0.43, 3.35] (griffin) or [0.32, 0.06] (sarpkaya)
 * @returns {number} single-amplitude A/D [-]
 */
export function responseAmplitude({ SG, gamma = 1, model = 'griffin', coefficients, override }) {
  num(SG, 'SG', { min: 0 }); num(gamma, 'gamma', { min: 0, strictMin: true });
  if (override !== undefined) return num(override, 'override', { min: 0 });
  if (model === 'griffin') { const c = coefficients ?? [1.29, 0.43, 3.35]; return c[0] * gamma / (1 + c[1] * SG) ** c[2]; }
  if (model === 'sarpkaya') { const c = coefficients ?? [0.32, 0.06]; return c[0] * gamma / Math.sqrt(c[1] + SG * SG); }
  throw new RangeError("viv: argument model must be 'griffin' or 'sarpkaya'");
}

/**
 * Drag amplification under cross-flow VIV: 'vandiver' Cd/Cd0 = 1 + 1.043·(2·Arms/D)^0.65 (Vandiver 1983, OTC 4490);
 * 'blevins' Cd/Cd0 = 1 + 2.1·(A/D) (Blevins 1990, §3.5).
 * @param {object} o `{ Cd0 [-], aRmsOverD? [-] (vandiver), aOverD? [-] (blevins; also used as Arms·√2 for vandiver),
 *   model? = 'vandiver' }`
 * @returns {{ Cd: number, factor: number }} amplified drag coefficient and amplification factor [-]
 */
export function dragAmplification({ Cd0, aRmsOverD, aOverD, model = 'vandiver' }) {
  num(Cd0, 'Cd0', { min: 0 });
  let f;
  if (model === 'vandiver') {
    const r = aRmsOverD ?? num(aOverD, 'aOverD', { min: 0 }) / Math.SQRT2;
    f = 1 + 1.043 * (2 * num(r, 'aRmsOverD', { min: 0 })) ** 0.65;
  } else if (model === 'blevins') f = 1 + 2.1 * num(aOverD ?? aRmsOverD * Math.SQRT2, 'aOverD', { min: 0 });
  else throw new RangeError("viv: argument model must be 'vandiver' or 'blevins'");
  return { Cd: Cd0 * f, factor: f };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Modes of a tensioned beam                                                                                 */

/**
 * Transverse modes of a pinned–pinned tensioned beam (riser, span, jumper leg) with uniform or linearly varying
 * effective tension.
 *
 * Uniform tension: closed form fn = (n/2L)·sqrt(T/m)·sqrt(1 + n²π²EI/(T·L²)), φn = sin(nπs/L)
 * (Blevins, "Formulas for Natural Frequency and Mode Shape", Table 8-? tensioned beam).
 * Varying tension T(s) = T0 + (T1 − T0)·s/L: Rayleigh–Ritz with a sine basis (exact integrals), which converges
 * to the exact Bessel-function solution of the hanging string for EI → 0. `wkbHz` gives the WKB string estimate
 * fn = n / (2·∫ds/c), c = sqrt(T(s)/m), for comparison.
 * Validity: uniform mass and EI, pinned ends, T > 0 everywhere (use rodModes for anything else).
 *
 * @param {object} o
 * @param {number} o.L length [m] @param {number} o.EI [N·m²] @param {number} o.m mass per length incl. added mass [kg/m]
 * @param {number} [o.T] uniform tension [N], or @param {number} [o.T0] tension at s = 0 and @param {number} [o.T1]
 *   at s = L [N]
 * @param {number} [o.nModes=5] @param {number} [o.nPoints=201] @param {number} [o.nBasis] sine terms (default
 *   3·nModes + 12) @param {boolean} [o.forceRitz=false]
 * @returns {{ s: Float64Array, modes: Array<{ n: number, fHz: number, wkbHz: number, phi: Float64Array,
 *   curvature: Float64Array }> }} s [m]; fHz [Hz]; phi [-] scaled to max|phi| = 1; curvature [1/m²]
 */
export function tensionedBeamModes(o) {
  const { L, EI, m, nModes = 5, nPoints = 201, forceRitz = false } = o ?? {};
  num(L, 'L', { min: 0, strictMin: true }); num(EI, 'EI', { min: 0 }); num(m, 'm', { min: 0, strictMin: true });
  const T0 = num(o.T0 ?? o.T, 'T0', { min: 0 }), T1 = num(o.T1 ?? o.T, 'T1', { min: 0 });
  if (!(T0 + T1 > 0) && !(EI > 0)) throw new RangeError('viv: arguments T and EI cannot both be zero');
  num(nModes, 'nModes', { min: 1, max: 200 });
  const s = new Float64Array(nPoints);
  for (let i = 0; i < nPoints; i++) s[i] = L * i / (nPoints - 1);
  const wkb = (n) => {
    if (!(T0 > 0 && T1 > 0)) return NaN;
    const tau = T0 === T1 ? L * Math.sqrt(m / T0) : 2 * L * Math.sqrt(m) * (Math.sqrt(T1) - Math.sqrt(T0)) / (T1 - T0);
    return n / (2 * tau);
  };
  const modes = [];
  if (T0 === T1 && !forceRitz) {
    for (let n = 1; n <= nModes; n++) {
      const k = n * Math.PI / L;
      const f = Math.sqrt((T0 * k * k + EI * k ** 4) / m) / (2 * Math.PI);
      const phi = new Float64Array(nPoints), curvature = new Float64Array(nPoints);
      for (let i = 0; i < nPoints; i++) { phi[i] = Math.sin(k * s[i]); curvature[i] = -k * k * phi[i]; }
      modes.push({ n, fHz: f, wkbHz: wkb(n), phi, curvature });
    }
    return { s, modes };
  }
  const nb = Math.round(o.nBasis ?? 3 * nModes + 12);
  const dT = T1 - T0, A = [];
  for (let j = 1; j <= nb; j++) {
    const row = new Array(nb);
    for (let k = 1; k <= nb; k++) {
      // ∫ T(s) φj' φk' ds with φ = sin(jπs/L); mass ∫ m φj φk ds = m L/2 δjk
      let kt;
      const pj = j * Math.PI / L, pk = k * Math.PI / L;
      if (j === k) kt = pj * pj * (T0 * L / 2 + dT * L / 4);
      else if ((j + k) % 2 === 1) {
        // ∫0^L s cos(pj s) cos(pk s) ds = ½[((−1)^(j−k) − 1)/(pj − pk)² + ((−1)^(j+k) − 1)/(pj + pk)²]
        kt = pj * pk * (dT / L) * 0.5 * (-2 / (pj - pk) ** 2 - 2 / (pj + pk) ** 2);
      } else kt = 0;
      const kb = j === k ? EI * pj ** 4 * L / 2 : 0;
      row[k - 1] = (kt + kb) / (m * L / 2);
    }
    A.push(row);
  }
  const eig = eigSym(A);
  for (let n = 1; n <= Math.min(nModes, nb); n++) {
    const lam = eig.values[n - 1], v = eig.vectors[n - 1];
    const phi = new Float64Array(nPoints), curvature = new Float64Array(nPoints);
    let big = 0;
    for (let i = 0; i < nPoints; i++) {
      let p = 0, c = 0;
      for (let j = 1; j <= nb; j++) { const kj = j * Math.PI / L, sn = Math.sin(kj * s[i]); p += v[j - 1] * sn; c -= v[j - 1] * kj * kj * sn; }
      phi[i] = p; curvature[i] = c;
      if (Math.abs(p) > Math.abs(big)) big = p;
    }
    for (let i = 0; i < nPoints; i++) { phi[i] /= big; curvature[i] /= big; }
    modes.push({ n, fHz: Math.sqrt(Math.max(lam, 0)) / (2 * Math.PI), wkbHz: wkb(n), phi, curvature });
  }
  return { s, modes };
}

function modeCurvature(s, mode) {
  if (mode.curvature && mode.curvature.length === s.length) return mode.curvature;
  const n = s.length, c = new Float64Array(n), p = mode.phi;
  for (let i = 1; i < n - 1; i++) {
    const h1 = s[i] - s[i - 1], h2 = s[i + 1] - s[i];
    c[i] = 2 * (p[i - 1] / (h1 * (h1 + h2)) - p[i] / (h1 * h2) + p[i + 1] / (h2 * (h1 + h2)));
  }
  c[0] = c[1]; c[n - 1] = c[n - 2];
  return c;
}

/* ------------------------------------------------------------------------------------------------------- */
/* Hydrodynamic damping and modal screening                                                                  */

/**
 * Hydrodynamic damping per unit length outside the power-in region, in the functional form of Venugopal (1996,
 * MIT PhD thesis) as used in modal VIV programs:
 *   still water       r_sw = (ω·π·ρ·D²/2)·[2√2/sqrt(Reω) + Csw·(A/D)²],  Reω = ω·D²/ν
 *   low reduced vel.  r_vl = r_sw + Cvl·ρ·D·V
 *   high reduced vel. r_vh = Cvh·ρ·V²/ω
 * The coefficients are arguments; the defaults Csw = 0.25, Cvl = 0.18, Cvh = 0.2 are the values commonly quoted
 * for this model and should be confirmed against the source before design use.
 * @param {object} o `{ zone: 'still' | 'low' | 'high', omega [rad/s], D [m], V? [m/s local current], aOverD? = 0 [-],
 *   rho? = 1025, nu? = 1.19e-6 [m²/s], Csw? = 0.25, Cvl? = 0.18, Cvh? = 0.2 }`
 * @returns {number} damping coefficient per unit length [N·s/m²]
 */
export function hydrodynamicDamping({ zone, omega, D, V = 0, aOverD = 0, rho = 1025, nu = 1.19e-6, Csw = 0.25, Cvl = 0.18, Cvh = 0.2 }) {
  num(omega, 'omega', { min: 0, strictMin: true }); num(D, 'D', { min: 0, strictMin: true }); num(V, 'V', { min: 0 });
  const rsw = (omega * Math.PI * rho * D * D / 2) * (2 * Math.SQRT2 / Math.sqrt(omega * D * D / nu) + Csw * aOverD * aOverD);
  if (zone === 'still') return rsw;
  if (zone === 'low') return rsw + Cvl * rho * D * V;
  if (zone === 'high') return Cvh * rho * V * V / omega;
  throw new RangeError("viv: argument zone must be 'still', 'low' or 'high'");
}

/**
 * Modal cross-flow VIV screening of a slender member in a sheared current.
 *
 * For each mode: power-in region = where the local reduced velocity U(s)/(fn·D) lies inside `vrWindow`
 * (weighted by 1 − suppression efficiency inside covered zones); effective mass
 * me = ∫m·φ²ds / ∫power-in φ²ds (DNV-RP-C205 §9.1.4 form); effective damping = structural + hydrodynamic damping
 * outside the power-in region (hydrodynamicDamping, iterated on amplitude); amplitude from responseAmplitude at
 * SG = 2π·St²·KS with the mode factor γ = sqrt(∫φ²/∫φ⁴) evaluated over the power-in region; an energy-balance
 * amplitude with constant lift coefficient Cl is reported alongside,
 * A = ∫in ½ρ·D·U²·Cl·|φ|ds / (ω·∫r·φ²ds), and the modal input power Π = ½·ω·A·∫in ½ρ·D·U²·Cl·|φ|ds ranks the modes.
 * Reference: Vandiver (1993), "Dimensionless parameters important to the prediction of VIV of long flexible
 * cylinders in ocean currents", J. Fluids Struct. 7; DNV-RP-C205 §9; DNV-RP-F204 §4.
 * Validity: screening level — single-mode lock-in per mode, cross-flow only, no mode competition in time, no
 * higher harmonics, constant St and Cl.
 *
 * @param {object} o
 * @param {Float64Array|number[]} o.s arc length [m] @param {Array<object>} o.modes `{ fHz, phi, curvature? }`
 * @param {number|Function|ArrayLike<number>} o.U current speed normal to the member along s [m/s]
 * @param {number|Function|ArrayLike<number>} o.D hydrodynamic diameter [m]
 * @param {number|Function|ArrayLike<number>} o.mass mass per length incl. added mass [kg/m]
 * @param {number} [o.zeta=0.003] structural damping ratio [-] @param {number} [o.rho=1025] @param {number} [o.St=0.2]
 * @param {[number, number]} [o.vrWindow=[4, 8]] power-in reduced-velocity window [-] (typical lock-in band; argument)
 * @param {number} [o.Cl=0.7] lift coefficient for the energy balance [-] (typical 0.3–0.9; argument)
 * @param {'griffin'|'sarpkaya'} [o.responseModel='griffin'] @param {number} [o.aOverDOverride] guidance amplitude
 * @param {{ zones: Array<[number, number]>, efficiency: number }} [o.suppression] covered arc-length intervals [m]
 *   and efficiency (fraction of excitation removed, 0–1) [-]
 * @param {object} [o.damping] overrides for hydrodynamicDamping `{ nu, Csw, Cvl, Cvh }`
 * @returns {{ modes: Array<object>, dominant: number, excited: number[] }} per mode `{ index, fHz, powerInLength [m],
 *   powerInFraction [-], vrMin, vrMax [-], me [kg/m], zetaEff [-], Ks, SG, gamma [-], aOverD [-] (single amplitude at
 *   the antinode), aRmsOverD [-], aOverDEnergy [-], power [W], dragFactor [-], excited }`; `dominant` = index of the
 *   mode with the largest input power (−1 if none); `excited` = indices sorted by power
 */
export function vivModalScreening(o) {
  const { s, modes, zeta = 0.003, rho = 1025, St = 0.2, vrWindow = [4, 8], Cl = 0.7, responseModel = 'griffin', suppression, damping = {} } = o ?? {};
  if (!s || s.length < 3 || !Array.isArray(modes) || modes.length === 0) throw new RangeError('viv: arguments s and modes are required');
  const n = s.length;
  const U = field(o.U, s, 'U'), D = field(o.D, s, 'D'), mass = field(o.mass, s, 'mass');
  num(zeta, 'zeta', { min: 0 }); num(vrWindow[0], 'vrWindow[0]', { min: 0 }); num(vrWindow[1], 'vrWindow[1]', { min: vrWindow[0], strictMin: true });
  const wt = new Float64Array(n).fill(1);
  if (suppression) {
    const eff = num(suppression.efficiency, 'suppression.efficiency', { min: 0, max: 1 });
    for (const [a, b] of suppression.zones ?? []) for (let i = 0; i < n; i++) if (s[i] >= a && s[i] <= b) wt[i] = 1 - eff;
  }
  const out = [];
  const inW = new Float64Array(n), f1 = new Float64Array(n), f2 = new Float64Array(n);
  modes.forEach((md, index) => {
    const fn = num(md.fHz, `modes[${index}].fHz`, { min: 0, strictMin: true }), om = 2 * Math.PI * fn, phi = md.phi;
    if (!phi || phi.length !== n) throw new RangeError(`viv: argument modes[${index}].phi must match s`);
    let pmax = 0, vrLo = Infinity, vrHi = -Infinity, Lin = 0;
    for (let i = 0; i < n; i++) {
      pmax = Math.max(pmax, Math.abs(phi[i]));
      const vr = Math.abs(U[i]) / (fn * D[i]);
      vrLo = Math.min(vrLo, vr); vrHi = Math.max(vrHi, vr);
      inW[i] = vr >= vrWindow[0] && vr <= vrWindow[1] ? wt[i] : 0;
    }
    for (let i = 1; i < n; i++) Lin += 0.5 * ((inW[i] > 0 ? 1 : 0) + (inW[i - 1] > 0 ? 1 : 0)) * (s[i] - s[i - 1]);
    for (let i = 0; i < n; i++) { const p = phi[i] / pmax; f1[i] = p * p; f2[i] = p ** 4; }
    const Mm = trap(f1, s, mass), I2in = trap(f1, s, inW), I4in = trap(f2, s, inW);
    const rec = { index, fHz: fn, powerInLength: Lin, powerInFraction: Lin / (s[n - 1] - s[0]), vrMin: vrLo, vrMax: vrHi, excited: I2in > 0 };
    if (!(I2in > 0)) { out.push({ ...rec, me: Infinity, zetaEff: NaN, Ks: Infinity, SG: Infinity, gamma: NaN, aOverD: 0, aRmsOverD: 0, aOverDEnergy: 0, power: 0, dragFactor: 1 }); return; }
    const me = Mm / I2in, gamma = Math.sqrt(I2in / I4in);
    let Dref = 0, wsum = 0;
    for (let i = 0; i < n; i++) { Dref += D[i] * inW[i] * f1[i]; wsum += inW[i] * f1[i]; }
    Dref /= wsum;
    // excitation integral ∫in ½ρ D U² Cl |φ| ds
    for (let i = 0; i < n; i++) f2[i] = 0.5 * rho * D[i] * U[i] * U[i] * Cl * Math.sqrt(f1[i]) * inW[i];
    const Fex = trap(f2, s);
    let aD = 0.5, zetaEff = zeta, Ks = 0, SG = 0, Rdamp = 0;
    for (let it = 0; it < 40; it++) {
      for (let i = 0; i < n; i++) {
        const vr = Math.abs(U[i]) / (fn * D[i]);
        const zone = vr < vrWindow[0] ? 'low' : 'high';
        const rh = hydrodynamicDamping({ zone, omega: om, D: D[i], V: Math.abs(U[i]), aOverD: aD * Math.sqrt(f1[i]), rho, ...damping });
        f2[i] = (1 - inW[i]) * rh * f1[i];
      }
      const Rh = trap(f2, s);
      zetaEff = zeta + Rh / (2 * om * Mm);
      Rdamp = 2 * zetaEff * om * Mm;
      Ks = 4 * Math.PI * me * zetaEff / (rho * Dref * Dref);
      SG = 2 * Math.PI * St * St * Ks;
      const aNew = responseAmplitude({ SG, gamma, model: responseModel, override: o.aOverDOverride });
      if (Math.abs(aNew - aD) < 1e-10) { aD = aNew; break; }
      aD = 0.5 * (aD + aNew);
    }
    const aE = Fex / (om * Rdamp) / Dref;
    const power = 0.5 * om * (aD * Dref) * Fex;
    out.push({ ...rec, me, zetaEff, Ks, SG, gamma, aOverD: aD, aRmsOverD: aD / Math.SQRT2, aOverDEnergy: aE, power,
      dragFactor: dragAmplification({ Cd0: 1, aRmsOverD: aD / Math.SQRT2 }).factor, Dref });
  });
  const excited = out.filter((m) => m.excited).sort((a, b) => b.power - a.power).map((m) => m.index);
  return { modes: out, dominant: excited.length ? excited[0] : -1, excited };
}

/**
 * VIV stress amplitude and fatigue damage rate along the member for one screened current case.
 * σ(s) = SCF·E·(Ds/2)·|κn(s)|·An (modal curvature × amplitude), stress range 2σ, damage per year =
 * fn·(seconds per year)/N(2σ) with the caller's S–N function. `modesUsed = 'dominant'` (default) takes the
 * dominant mode; 'all' sums the damage of every excited mode (conservative, no time sharing).
 * Reference: DNV-RP-F204 §4 (VIV fatigue), DNV-RP-C203 §2.
 * @param {object} o `{ s, modes, screening (vivModalScreening result), E? = 207e9 [Pa], stressDiameter [m outer
 *   diameter of the load-bearing wall], scf? = 1, snCyclesToFailure: (stressRangePa) => N, modesUsed? }`
 * @returns {{ s: Float64Array, stressAmplitude: Float64Array, damagePerYear: Float64Array, maxDamagePerYear: number,
 *   sMax: number, lifeYears: number }} stress [Pa] (dominant mode), damage [1/yr], location [m], life [yr]
 */
export function vivFatigue({ s, modes, screening, E = 207e9, stressDiameter, scf = 1, snCyclesToFailure, modesUsed = 'dominant' }) {
  num(stressDiameter, 'stressDiameter', { min: 0, strictMin: true }); num(scf, 'scf', { min: 0, strictMin: true });
  if (typeof snCyclesToFailure !== 'function') throw new RangeError('viv: argument snCyclesToFailure must be a function');
  const n = s.length, dmg = new Float64Array(n), sig = new Float64Array(n);
  const list = modesUsed === 'all' ? screening.excited : (screening.dominant >= 0 ? [screening.dominant] : []);
  list.forEach((idx, rank) => {
    const sm = screening.modes[idx], md = modes[idx];
    let pmax = 0; for (let i = 0; i < n; i++) pmax = Math.max(pmax, Math.abs(md.phi[i]));
    const kap = modeCurvature(s, md), A = sm.aOverD * sm.Dref / pmax;
    for (let i = 0; i < n; i++) {
      const sa = scf * E * 0.5 * stressDiameter * Math.abs(kap[i]) * A;
      if (rank === 0) sig[i] = sa;
      const N = snCyclesToFailure(2 * sa);
      if (N > 0 && Number.isFinite(N)) dmg[i] += sm.fHz * SEC_PER_YEAR / N;
    }
  });
  let mx = 0, sMax = s[0];
  for (let i = 0; i < n; i++) if (dmg[i] > mx) { mx = dmg[i]; sMax = s[i]; }
  return { s: Float64Array.from(s), stressAmplitude: sig, damagePerYear: dmg, maxDamagePerYear: mx, sMax, lifeYears: mx > 0 ? 1 / mx : Infinity };
}

/**
 * Long-term VIV fatigue from a scatter of current profiles: damage(s) = Σ pi·damage_i(s), life = 1/max damage.
 * Reference: DNV-RP-F204 §4.2 (weighting by probability of occurrence of each current profile).
 * @param {object} o vivModalScreening and vivFatigue arguments (without U) plus
 *   `cases: Array<{ probability [-], U }>` (probabilities should sum to ≤ 1; the remainder is treated as no VIV)
 * @returns {{ s: Float64Array, damagePerYear: Float64Array, lifeYears: number, sMax: number,
 *   cases: Array<{ probability: number, dominantMode: number, aOverD: number, maxDamagePerYear: number }> }}
 */
export function vivFatigueLife(o) {
  const { cases } = o ?? {};
  if (!Array.isArray(cases) || cases.length === 0) throw new RangeError('viv: argument cases must be a non-empty array');
  const n = o.s.length, dmg = new Float64Array(n), rows = [];
  let ptot = 0;
  for (const c of cases) {
    const p = num(c.probability, 'cases[].probability', { min: 0, max: 1 }); ptot += p;
    const scr = vivModalScreening({ ...o, U: c.U });
    const fat = vivFatigue({ ...o, screening: scr });
    for (let i = 0; i < n; i++) dmg[i] += p * fat.damagePerYear[i];
    rows.push({ probability: p, dominantMode: scr.dominant, aOverD: scr.dominant >= 0 ? scr.modes[scr.dominant].aOverD : 0, maxDamagePerYear: fat.maxDamagePerYear });
  }
  if (ptot > 1 + 1e-9) throw new RangeError('viv: argument cases[].probability sums to more than 1');
  let mx = 0, sMax = o.s[0];
  for (let i = 0; i < n; i++) if (dmg[i] > mx) { mx = dmg[i]; sMax = o.s[i]; }
  return { s: Float64Array.from(o.s), damagePerYear: dmg, lifeYears: mx > 0 ? 1 / mx : Infinity, sMax, cases: rows };
}

/**
 * Suppression coverage (strakes / fairings) needed to bring the dominant-mode amplitude down to a target:
 * coverage is grown from one end of the member and the screening repeated (bisection on covered length).
 * @param {object} o vivModalScreening arguments plus `{ efficiency [-], targetAOverD [-], from? = 'start' | 'end',
 *   tol? = 1e-3 [fraction of length] }`
 * @returns {{ coveredLength: number, coverageFraction: number, zone: [number, number], aOverD: number,
 *   achievable: boolean }} length [m]; amplitude with that coverage [-]; achievable = false when even full
 *   coverage does not reach the target
 */
export function suppressionCoverage(o) {
  const { efficiency, targetAOverD, from = 'start', tol = 1e-3 } = o ?? {};
  num(efficiency, 'efficiency', { min: 0, max: 1 }); num(targetAOverD, 'targetAOverD', { min: 0 });
  const s0 = o.s[0], s1 = o.s[o.s.length - 1], L = s1 - s0;
  const zoneOf = (f) => (from === 'end' ? [s1 - f * L, s1] : [s0, s0 + f * L]);
  const amp = (f) => {
    const scr = vivModalScreening({ ...o, suppression: f > 0 ? { zones: [zoneOf(f)], efficiency } : undefined });
    let a = 0; for (const m of scr.modes) if (m.excited) a = Math.max(a, m.aOverD);
    return a;
  };
  if (amp(0) <= targetAOverD) return { coveredLength: 0, coverageFraction: 0, zone: zoneOf(0), aOverD: amp(0), achievable: true };
  if (amp(1) > targetAOverD) return { coveredLength: L, coverageFraction: 1, zone: zoneOf(1), aOverD: amp(1), achievable: false };
  let lo = 0, hi = 1;
  while (hi - lo > tol) { const mid = 0.5 * (lo + hi); if (amp(mid) > targetAOverD) lo = mid; else hi = mid; }
  return { coveredLength: hi * L, coverageFraction: hi, zone: zoneOf(hi), aOverD: amp(hi), achievable: true };
}

/**
 * Reduced-velocity screening for members whose natural frequencies are supplied by the caller (free spans,
 * jumpers, conductors): reduced velocity, shedding frequency and in-line / cross-flow flags per frequency.
 * @param {object} o `{ frequenciesHz: number[], U [m/s], D [m], St? = 0.2, Ks? [-], ranges? (lockIn limit overrides) }`
 * @returns {Array<{ fHz: number, Vr: number, fsOverFn: number, inline: boolean, crossflow: boolean,
 *   crossflowPeak: boolean }>}
 */
export function spanVivScreening({ frequenciesHz, U, D, St = 0.2, Ks, ranges = {} }) {
  if (!Array.isArray(frequenciesHz)) throw new RangeError('viv: argument frequenciesHz must be an array');
  const fs = sheddingFrequency({ U, D, St });
  return frequenciesHz.map((fHz) => {
    const Vr = reducedVelocity({ U, fn: fHz, D });
    return { fHz, Vr, fsOverFn: fs / fHz, ...lockIn({ Vr, Ks, ...ranges }) };
  });
}

/* ------------------------------------------------------------------------------------------------------- */
/* Wake oscillator                                                                                           */

/**
 * Wake-oscillator model of cross-flow VIV of an elastically mounted rigid cylinder (Facchinetti, de Langre &
 * Biolley 2004, "Coupling of structure and wake oscillators in vortex-induced vibrations", J. Fluids Struct. 19),
 * acceleration coupling. Dimensionless equations (time t = Ωf·T, Ωf = 2π·St·U/D, y = Y/D, q = 2·CL/CL0):
 *   ÿ + (2ξδ + γ/μ)·ẏ + δ²·y = M·q,       q̈ + ε·(q² − 1)·q̇ + q = A·ÿ
 * with δ = 1/(St·Ur), Ur = U/(fs·D) (fs = structural frequency in still water), μ = (ms + mf)/(ρD²),
 * mf = CM·ρ·π·D²/4, γ = CD/(4π·St) (fluid-added damping), M = CL0/(16π²·St²·μ). Published values A = 12, ε = 0.3,
 * CL0 = 0.3, CD = 2.0 (γ = 0.8 at St = 0.2), CM = 1 are the defaults.
 * At lock-in (δ = 1) harmonic balance of these equations gives
 *   y0 = (CL0/2)/(SG + 4π²St²γ) · sqrt(1 + (A/ε)·(CL0/4)/(SG + 4π²St²γ)),   SG = 8π²·St²·μ·ξ,
 * returned as `lockInAmplitude` for cross-checking. Validity: one-degree-of-freedom rigid cylinder, uniform flow.
 *
 * @param {object} o
 * @param {number} o.massRatio structural mass / displaced fluid mass m* = ms/(ρπD²/4) [-]
 * @param {number} [o.zeta=0.001] structural damping ratio ξ [-]
 * @param {number[]} [o.Ur] reduced velocities (default 2 … 12 in steps of 0.25) [-]
 * @param {number} [o.St=0.2] @param {number} [o.A=12] @param {number} [o.eps=0.3] @param {number} [o.CL0=0.3]
 * @param {number} [o.CD=2.0] @param {number} [o.CM=1.0]
 * @param {number} [o.cycles=250] simulated wake cycles per point [-] @param {number} [o.stepsPerCycle=40]
 * @returns {{ Ur: number[], aOverD: number[], liftMagnification: number[], frequencyRatio: number[], peakUr: number,
 *   peakAOverD: number, lockInAmplitude: number, mu: number, gamma: number, M: number, SG: number }} amplitude A/D
 *   [-], q amplitude / 2 [-], response frequency / structural frequency [-] at each Ur; peak values; parameters
 */
export function wakeOscillator(o) {
  const { massRatio, zeta = 0.001, St = 0.2, A = 12, eps = 0.3, CL0 = 0.3, CD = 2.0, CM = 1.0, cycles = 250, stepsPerCycle = 40 } = o ?? {};
  num(massRatio, 'massRatio', { min: 0, strictMin: true }); num(zeta, 'zeta', { min: 0 });
  num(St, 'St', { min: 0, strictMin: true }); num(A, 'A', { min: 0 }); num(eps, 'eps', { min: 0, strictMin: true });
  let Ur = o.Ur;
  if (!Ur) { Ur = []; for (let v = 2; v <= 12 + 1e-9; v += 0.25) Ur.push(v); }
  const mu = (massRatio + CM) * Math.PI / 4, gamma = CD / (4 * Math.PI * St), M = CL0 / (16 * Math.PI ** 2 * St * St * mu);
  const SG = 8 * Math.PI ** 2 * St * St * mu * zeta;
  const aOverD = [], lift = [], fr = [];
  let state = [0, 0, 0.1, 0];
  const tEnd = 2 * Math.PI * cycles, nSteps = Math.round(cycles * stepsPerCycle);
  for (const ur of Ur) {
    num(ur, 'Ur[]', { min: 0, strictMin: true });
    const d = 1 / (St * ur), c = 2 * zeta * d + gamma / mu;
    const rhs = (t, y) => {
      const ydd = M * y[2] - c * y[1] - d * d * y[0];
      return [y[1], ydd, y[3], A * ydd - eps * (y[2] * y[2] - 1) * y[3] - y[2]];
    };
    let ymax = 0, qmax = 0, cross = 0, tFirst = 0, tLast = 0, prev = 0;
    const start = 0.6 * tEnd;
    const sol = rk4(rhs, 0, state, tEnd, nSteps, (t, y) => {
      if (t >= start) {
        ymax = Math.max(ymax, Math.abs(y[0])); qmax = Math.max(qmax, Math.abs(y[2]));
        if (prev < 0 && y[0] >= 0) { if (cross === 0) tFirst = t; tLast = t; cross++; }
      }
      prev = y[0];
    });
    const last = sol.y[sol.y.length - 1];
    state = [last[0], last[1], Math.abs(last[2]) > 1e-3 ? last[2] : 0.1, last[3]];   // continuation in Ur
    aOverD.push(ymax); lift.push(qmax / 2);
    fr.push(cross > 1 ? ((cross - 1) / (tLast - tFirst)) * 2 * Math.PI / d : NaN);
  }
  let k = 0; for (let i = 1; i < Ur.length; i++) if (aOverD[i] > aOverD[k]) k = i;
  const den = SG + 4 * Math.PI ** 2 * St * St * gamma;
  return {
    Ur: Array.from(Ur), aOverD, liftMagnification: lift, frequencyRatio: fr, peakUr: Ur[k], peakAOverD: aOverD[k],
    lockInAmplitude: (CL0 / 2) / den * Math.sqrt(1 + (A / eps) * (CL0 / 4) / den), mu, gamma, M, SG,
  };
}

/**
 * Exact natural frequency of mode n of a string with linearly varying tension (hanging riser without bending
 * stiffness): roots of J0(a)·Y0(b) − J0(b)·Y0(a) = 0, a = 2ω·sqrt(m·T0)/q, b = 2ω·sqrt(m·T1)/q, q = (T1 − T0)/L
 * (Irvine, "Cable Structures", §3; classical heavy-string solution). The Bessel functions are supplied by the
 * caller (math/special.js besselJ0, besselY0) to keep this module free of that dependency.
 * @param {object} o `{ n [-], L [m], m [kg/m], T0 [N], T1 [N] (T1 > T0 > 0), besselJ0, besselY0 }`
 * @returns {number} frequency [Hz]
 */
export function hangingStringFrequency({ n, L, m, T0, T1, besselJ0, besselY0 }) {
  num(n, 'n', { min: 1 }); num(L, 'L', { min: 0, strictMin: true }); num(m, 'm', { min: 0, strictMin: true });
  num(T0, 'T0', { min: 0, strictMin: true }); num(T1, 'T1', { min: T0, strictMin: true });
  const q = (T1 - T0) / L;
  const F = (w) => { const a = 2 * w * Math.sqrt(m * T0) / q, b = 2 * w * Math.sqrt(m * T1) / q; return besselJ0(a) * besselY0(b) - besselJ0(b) * besselY0(a); };
  // roots interlace close to the WKB values ωn = nπ / ∫ds/c
  const tau = 2 * Math.sqrt(m) * (Math.sqrt(T1) - Math.sqrt(T0)) / q;
  const w0 = (n - 0.5) * Math.PI / tau, w1 = (n + 0.5) * Math.PI / tau;
  return bisect(F, w0, w1, { tol: 1e-12 }) / (2 * Math.PI);
}
