// Wave-energy converters: resource, heaving point absorber (frequency domain, PTO control, limits), irregular
// waves and power matrices, Cummins time-domain model with nonlinear PTO, two-body absorber, bottom-hinged flap,
// simple parameterised devices, array factor and survival statistics.
//
// Units: SI. Frequencies are angular (rad/s) unless named otherwise; spectra S(ω) in m²·s/rad (one-sided, so
// that ∫ S dω = Hs²/16); excitation coefficients are force per unit wave amplitude (N/m, or N·m/m for flaps).

import { besselJ1 } from '../math/special.js';

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Linear dispersion relation ω² = g k tanh(k h): wavenumber, phase and group velocity.
 * @param {number} omega angular frequency (rad/s)
 * @param {number} [depth=Infinity] water depth (m)
 * @returns {{k:number, wavelength:number, c:number, cg:number}} k (rad/m), λ (m), phase and group velocity (m/s)
 */
export function dispersion(omega, depth = Infinity) {
  finitePos(omega, 'omega');
  need(depth > 0, 'depth must be > 0');
  const k0 = (omega * omega) / G;
  if (!Number.isFinite(depth) || k0 * depth > 40) {
    return { k: k0, wavelength: (2 * Math.PI) / k0, c: G / omega, cg: G / (2 * omega) };
  }
  let k = k0 / Math.sqrt(Math.tanh(k0 * depth));
  for (let i = 0; i < 60; i++) {
    const th = Math.tanh(k * depth);
    const dk = (G * k * th - omega * omega) / (G * th + G * k * depth * (1 - th * th));
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  const c = omega / k;
  const cg = 0.5 * c * (1 + (2 * k * depth) / Math.sinh(2 * k * depth));
  return { k, wavelength: (2 * Math.PI) / k, c, cg };
}

/**
 * Wave energy flux per unit crest length of an irregular sea state.
 * Deep water: J = ρ g² Hm0² Te / (64 π). Finite depth: J = ρ g (Hm0²/16) c_g(Te, h), the group velocity being
 * evaluated at the energy period (exact for a narrow spectrum; IEC TS 62600-101 integrates over the spectrum).
 * Regular wave of height H and period T: J = ρ g H² c_g / 8.
 * @param {object} o
 * @param {number} o.Hm0 significant wave height (m)
 * @param {number} o.Te energy period m−1/m0 (s)
 * @param {number} [o.depth=Infinity] water depth (m)
 * @param {number} [o.rho=1025] sea-water density (kg/m³)
 * @returns {number} J (W/m)
 */
export function wavePower({ Hm0, Te, depth = Infinity, rho = 1025 }) {
  need(Number.isFinite(Hm0) && Hm0 >= 0, 'Hm0 must be >= 0');
  finitePos(Te, 'Te');
  finitePos(rho, 'rho');
  if (!Number.isFinite(depth)) return (rho * G * G * Hm0 * Hm0 * Te) / (64 * Math.PI);
  return ((rho * G * Hm0 * Hm0) / 16) * dispersion((2 * Math.PI) / Te, depth).cg;
}

/**
 * Energy flux of a regular wave, J = ½ ρ g a² c_g (W/m).
 * @param {object} o
 * @param {number} o.amplitude wave amplitude a = H/2 (m)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} [o.depth=Infinity]
 * @param {number} [o.rho=1025]
 * @returns {number} J (W/m)
 */
export function regularWavePower({ amplitude, omega, depth = Infinity, rho = 1025 }) {
  need(amplitude >= 0, 'amplitude must be >= 0');
  return 0.5 * rho * G * amplitude * amplitude * dispersion(omega, depth).cg;
}

/**
 * Bretschneider / JONSWAP wave spectrum in angular frequency,
 *   S(ω) = A_γ (5/16) Hs² ω_p⁴ ω^−5 exp(−1.25 (ω_p/ω)⁴) γ^exp(−(ω − ω_p)²/(2 σ² ω_p²)),
 * σ = 0.07 (ω ≤ ω_p), 0.09 (ω > ω_p), A_γ = 1 − 0.287 ln γ. γ = 1 is the Bretschneider (two-parameter
 * Pierson–Moskowitz) form, for which m0 = Hs²/16 exactly. Reference: DNV-RP-C205 (2021) §3.5.5.
 * @param {object} o
 * @param {number} o.Hs significant wave height (m)
 * @param {number} o.Tp peak period (s)
 * @param {number} [o.gamma=1] peak-enhancement factor (-), 1–7
 * @param {number[]|Float64Array} o.omega angular frequencies (rad/s)
 * @returns {Float64Array} S(ω) (m²·s/rad)
 */
export function waveSpectrum({ Hs, Tp, gamma = 1, omega }) {
  need(Hs >= 0, 'Hs must be >= 0');
  finitePos(Tp, 'Tp');
  need(gamma >= 1 && gamma <= 10, 'gamma must be within [1, 10]');
  const wp = (2 * Math.PI) / Tp;
  const Ag = 1 - 0.287 * Math.log(gamma);
  const S = new Float64Array(omega.length);
  for (let i = 0; i < omega.length; i++) {
    const w = omega[i];
    if (!(w > 0)) continue;
    const sig = w <= wp ? 0.07 : 0.09;
    const r = Math.exp(-((w - wp) ** 2) / (2 * sig * sig * wp * wp));
    S[i] = Ag * (5 / 16) * Hs * Hs * wp ** 4 * w ** -5 * Math.exp(-1.25 * (wp / w) ** 4) * gamma ** r;
  }
  return S;
}

function trapz(y, x) {
  let s = 0;
  for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return s;
}

/**
 * Spectral periods and significant height from a spectrum: Hm0 = 4 √m0, Te = m−1/m0, Tz = 2π √(m0/m2)
 * (moments in angular frequency: T = 2π m_(n−1)/m_n).
 * @param {object} o
 * @param {number[]|Float64Array} o.omega angular frequencies (rad/s), ascending
 * @param {number[]|Float64Array} o.S spectrum (m²·s/rad)
 * @returns {{Hm0:number, Te:number, Tz:number, m0:number}}
 */
export function spectralPeriods({ omega, S }) {
  need(omega && S && omega.length === S.length && omega.length > 2, 'omega and S must be arrays of equal length > 2');
  const mom = (n) => trapz(Array.from(S, (s, i) => (omega[i] > 0 ? s * omega[i] ** n : 0)), omega);
  const m0 = mom(0);
  return { Hm0: 4 * Math.sqrt(m0), Te: (2 * Math.PI * mom(-1)) / m0, Tz: 2 * Math.PI * Math.sqrt(m0 / mom(2)), m0 };
}

/**
 * Mean wave power of a scatter diagram (joint occurrence of Hm0 and Te).
 * @param {object} o
 * @param {number[]} o.hs bin-centre significant wave heights (m)
 * @param {number[]} o.te bin-centre energy periods (s)
 * @param {number[][]} o.prob occurrence prob[iHs][iTe] (any scale; normalised internally)
 * @param {number} [o.depth=Infinity]
 * @param {number} [o.rho=1025]
 * @returns {{meanJ:number, annualWhPerM:number, J:number[][], contribution:number[][]}} mean flux (W/m), annual
 *   energy per metre crest (Wh/m), flux per cell (W/m) and each cell's share of the mean (-)
 */
export function scatterResource({ hs, te, prob, depth = Infinity, rho = 1025 }) {
  need(Array.isArray(hs) && Array.isArray(te) && Array.isArray(prob) && prob.length === hs.length, 'prob must have one row per hs value');
  let tot = 0;
  for (const row of prob) {
    need(row.length === te.length, 'each prob row must have one entry per te value');
    for (const p of row) {
      need(p >= 0, 'prob entries must be >= 0');
      tot += p;
    }
  }
  finitePos(tot, 'sum of prob');
  let mean = 0;
  const J = hs.map((h) => te.map((t) => wavePower({ Hm0: h, Te: t, depth, rho })));
  for (let i = 0; i < hs.length; i++) for (let j = 0; j < te.length; j++) mean += (J[i][j] * prob[i][j]) / tot;
  const contribution = J.map((row, i) => row.map((v, j) => (mean > 0 ? (v * prob[i][j]) / tot / mean : 0)));
  return { meanJ: mean, annualWhPerM: mean * 8760, J, contribution };
}

/**
 * Monthly mean wave power from one scatter diagram per month.
 * @param {object} o
 * @param {number[]} o.hs (m)
 * @param {number[]} o.te (s)
 * @param {number[][][]} o.monthly scatter diagrams prob[month][iHs][iTe]
 * @param {number} [o.depth=Infinity]
 * @param {number} [o.rho=1025]
 * @returns {number[]} mean J per month (W/m)
 */
export function monthlyResource({ hs, te, monthly, depth = Infinity, rho = 1025 }) {
  need(Array.isArray(monthly) && monthly.length > 0, 'monthly must be a non-empty array of scatter diagrams');
  return monthly.map((prob) => scatterResource({ hs, te, prob, depth, rho }).meanJ);
}

/**
 * Heave hydrodynamic coefficients of a truncated vertical circular cylinder (radius a, draft d) by documented
 * approximations — for screening when no diffraction/radiation table is available.
 *  - Hydrostatic stiffness K = ρ g π a².
 *  - Added mass A = c_a ρ D³ with c_a = 1/6: half the added mass of a thin circular disc in unbounded fluid,
 *    (8/3) ρ a³ / 2 (Lamb, Hydrodynamics §108; quoted as 0.167 ρ D³ by McCormick 1981). Frequency-independent
 *    here; c_a = π/24 gives the hemisphere high-frequency value.
 *  - Excitation: Froude–Krylov force on the flat bottom, integrating the undisturbed dynamic pressure over the
 *    disc, |F_e|/ζ = ρ g π a² P(d) · 2 J1(k a)/(k a) · f_d, with P(d) = e^(−k d) (deep) or
 *    cosh(k (h − d))/cosh(k h) (finite depth) and f_d an optional diffraction factor (default 1).
 *  - Radiation damping from the Haskind relation for a vertical-axis body in heave (Newman 1962; Falnes 2002,
 *    Ocean Waves and Oscillating Systems, ch. 5): B = k |F_e/ζ|² / (4 ρ g c_g), i.e. ω k |F_e/ζ|² / (2 ρ g²) in deep water.
 * With this B the maximum absorbed power |F_e|²/(8B) equals J/k, the axisymmetric heave limit.
 * Valid for small bodies (k a ≲ 1); errors grow as diffraction becomes important.
 * @param {object} o
 * @param {number} o.radius cylinder radius a (m)
 * @param {number} o.draft draft d (m)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} [o.depth=Infinity] water depth (m)
 * @param {number} [o.rho=1025]
 * @param {number} [o.addedMassCoeff=1/6] c_a in A = c_a ρ D³ (-)
 * @param {number} [o.diffractionFactor=1] multiplier f_d on the Froude–Krylov force (-)
 * @returns {{A:number, B:number, K:number, Fe:number, k:number, cg:number, displacedMass:number}} added mass
 *   (kg), radiation damping (N·s/m), stiffness (N/m), excitation per unit amplitude (N/m), k (rad/m), c_g (m/s)
 */
export function cylinderHeave({ radius, draft, omega, depth = Infinity, rho = 1025, addedMassCoeff = 1 / 6, diffractionFactor = 1 }) {
  finitePos(radius, 'radius');
  finitePos(draft, 'draft');
  finitePos(rho, 'rho');
  need(depth > draft, 'depth must exceed draft');
  const { k, cg } = dispersion(omega, depth);
  const Aw = Math.PI * radius * radius;
  const ka = k * radius;
  const jf = ka < 1e-6 ? 1 : (2 * besselJ1(ka)) / ka;
  const P = Number.isFinite(depth) && k * depth < 40
    ? Math.cosh(k * (depth - draft)) / Math.cosh(k * depth)
    : Math.exp(-k * draft);
  const Fe = rho * G * Aw * P * jf * diffractionFactor;
  return {
    A: addedMassCoeff * rho * (2 * radius) ** 3,
    B: (k * Fe * Fe) / (4 * rho * G * cg),
    K: rho * G * Aw,
    Fe: Math.abs(Fe),
    k,
    cg,
    displacedMass: rho * Aw * draft,
  };
}

function respond({ M, B, K, Fe, omega, amplitude, Bpto, Kpto }) {
  const re = K + Kpto - omega * omega * M;
  const im = omega * (B + Bpto);
  const x = (Fe * amplitude) / Math.hypot(re, im);
  const v = omega * x;
  return { x, v, F: x * Math.hypot(omega * Bpto, Kpto), P: 0.5 * Bpto * v * v };
}

/**
 * Single-degree-of-freedom wave absorber in a regular wave (frequency domain):
 *   [−ω² (m + A) + i ω (B + B_v + B_pto) + K + K_pto] X = F_e ζ,   P = ½ B_pto ω² |X|².
 * PTO strategies:
 *  - 'passive':  optimal linear damper, B_pto = √((B + B_v)² + (ω (m + A) − K/ω)²)  (Falnes 2002, ch. 3);
 *  - 'reactive': complex-conjugate control, B_pto = B + B_v, K_pto = ω² (m + A) − K, giving
 *                P_max = |F_e ζ|² / (8 (B + B_v));
 *  - 'fixed':    the given Bpto and Kpto.
 * Limits on stroke, velocity and PTO force amplitude are enforced by re-tuning the PTO damping (the stiffness
 * setting is kept): the damping maximising power among the settings that satisfy all limits is returned;
 * `feasible` is false when no damping satisfies them.
 * @param {object} o
 * @param {number} o.mass body mass m (kg)
 * @param {number} o.A added mass (kg)
 * @param {number} o.B radiation damping (N·s/m)
 * @param {number} o.K hydrostatic (plus mooring) stiffness (N/m)
 * @param {number} o.Fe excitation force per unit wave amplitude (N/m)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.amplitude wave amplitude ζ (m)
 * @param {'passive'|'reactive'|'fixed'} [o.control='passive']
 * @param {number} [o.Bpto] PTO damping (N·s/m) for 'fixed'
 * @param {number} [o.Kpto=0] PTO stiffness (N/m) for 'fixed' (may be negative)
 * @param {number} [o.Bvisc=0] linearised viscous damping (N·s/m); for quadratic drag ½ ρ Cd A_d v|v| the
 *   harmonic (Lorentz) equivalent is (8/(3π)) ½ ρ Cd A_d |v|
 * @param {{stroke?:number, velocity?:number, force?:number}} [o.limits] amplitude limits (m, m/s, N)
 * @returns {{power:number, stroke:number, velocity:number, ptoForce:number, Bpto:number, Kpto:number,
 *   constrained:boolean, feasible:boolean, powerUnconstrained:number, powerMax:number}} absorbed mean power
 *   (W), motion amplitudes, PTO force amplitude (N); powerMax = |F_e ζ|²/(8 (B + B_v))
 */
export function pointAbsorber({ mass, A, B, K, Fe, omega, amplitude, control = 'passive', Bpto, Kpto = 0, Bvisc = 0, limits = {} }) {
  finitePos(mass, 'mass');
  need(A >= 0, 'A must be >= 0');
  finitePos(B + Bvisc, 'B + Bvisc');
  need(K >= 0, 'K must be >= 0');
  finitePos(omega, 'omega');
  need(amplitude >= 0, 'amplitude must be >= 0');
  const M = mass + A;
  const Bt = B + Bvisc;
  let bp;
  let kp = Kpto;
  if (control === 'passive') {
    kp = 0;
    bp = Math.hypot(Bt, omega * M - K / omega);
  } else if (control === 'reactive') {
    kp = omega * omega * M - K;
    bp = Bt;
  } else {
    need(control === 'fixed', "control must be 'passive', 'reactive' or 'fixed'");
    need(Number.isFinite(Bpto) && Bpto >= 0, 'Bpto must be >= 0 for fixed control');
    bp = Bpto;
  }
  const base = { M, B: Bt, K, Fe, omega, amplitude, Kpto: kp };
  const free = respond({ ...base, Bpto: bp });
  const ok = (r) =>
    (limits.stroke === undefined || r.x <= limits.stroke * (1 + 1e-9)) &&
    (limits.velocity === undefined || r.v <= limits.velocity * (1 + 1e-9)) &&
    (limits.force === undefined || r.F <= limits.force * (1 + 1e-9));
  let best = free;
  let bBest = bp;
  let constrained = false;
  let feasible = true;
  if (!ok(free)) {
    constrained = true;
    feasible = false;
    const n = 600;
    const lo = Math.log(bp * 1e-4);
    const hi = Math.log(bp * 1e5);
    let bi = -1;
    let pBest = -1;
    const at = (i) => Math.exp(lo + ((hi - lo) * i) / n);
    for (let i = 0; i <= n; i++) {
      const r = respond({ ...base, Bpto: at(i) });
      if (ok(r) && r.P > pBest) {
        pBest = r.P;
        bi = i;
      }
    }
    if (bi >= 0) {
      feasible = true;
      bBest = at(bi);
      // refine toward the neighbouring infeasible setting with higher unconstrained power
      for (const nb of [bi - 1, bi + 1]) {
        if (nb < 0 || nb > n) continue;
        const rn = respond({ ...base, Bpto: at(nb) });
        if (ok(rn) || rn.P <= pBest) continue;
        let a = bBest;
        let b = at(nb);
        for (let it = 0; it < 60; it++) {
          const mid = Math.sqrt(a * b);
          if (ok(respond({ ...base, Bpto: mid }))) a = mid;
          else b = mid;
        }
        bBest = a;
      }
      best = respond({ ...base, Bpto: bBest });
    }
  }
  return {
    power: feasible ? best.P : 0,
    stroke: best.x,
    velocity: best.v,
    ptoForce: best.F,
    Bpto: bBest,
    Kpto: kp,
    constrained,
    feasible,
    powerUnconstrained: free.P,
    powerMax: (Fe * amplitude) ** 2 / (8 * Bt),
  };
}

/**
 * Capture width and capture-width ratio: CW = P / J (m), CWR = CW / characteristic width.
 * @param {object} o
 * @param {number} o.power absorbed power (W)
 * @param {number} o.J incident wave power per unit crest (W/m)
 * @param {number} [o.width] characteristic device width (m)
 * @returns {{captureWidth:number, ratio:number}} CW (m); CWR (-), NaN without a width
 */
export function captureWidth({ power, J, width }) {
  finitePos(J, 'J');
  return { captureWidth: power / J, ratio: width > 0 ? power / J / width : NaN };
}

/**
 * Theoretical limits for a heaving axisymmetric absorber.
 *  - Capture width ≤ λ/(2π) = 1/k (Budal & Falnes 1975; Evans 1976; Newman 1976).
 *  - Budal's upper bound for a body of volume V: with the excitation force bounded by the hydrostatic value
 *    ρ g A_w H/2 and the stroke amplitude by V/(2 A_w), P = ½ |F_e| |u| ≤ (π/4) ρ g V H / T
 *    (Budal & Falnes 1980; Falnes 2002, ch. 6).
 * @param {object} o
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.H wave height (m)
 * @param {number} o.volume body volume V (m³)
 * @param {number} [o.depth=Infinity]
 * @param {number} [o.rho=1025]
 * @returns {{maxCaptureWidth:number, radiationLimit:number, budalLimit:number, limit:number}} m, W, W, W
 */
export function heaveLimits({ omega, H, volume, depth = Infinity, rho = 1025 }) {
  need(H >= 0, 'H must be >= 0');
  need(volume >= 0, 'volume must be >= 0');
  const { k } = dispersion(omega, depth);
  const J = regularWavePower({ amplitude: H / 2, omega, depth, rho });
  const T = (2 * Math.PI) / omega;
  const radiationLimit = J / k;
  const budalLimit = ((Math.PI / 4) * rho * G * volume * H) / T;
  return { maxCaptureWidth: 1 / k, radiationLimit, budalLimit, limit: Math.min(radiationLimit, budalLimit) };
}

/**
 * Body description for the spectral and matrix functions from the truncated-cylinder approximations.
 * @param {object} o as for `cylinderHeave` without omega, plus
 * @param {number} [o.mass] body mass (kg); default the displaced mass ρ π a² d
 * @param {number} [o.Bvisc=0] linear viscous damping (N·s/m)
 * @returns {(omega:number)=>{M:number, B:number, K:number, Fe:number}} total inertia (kg), damping (N·s/m),
 *   stiffness (N/m), excitation per unit amplitude (N/m)
 */
export function cylinderBody({ mass, Bvisc = 0, ...geom }) {
  const probe = cylinderHeave({ ...geom, omega: 1 });
  const m = mass ?? probe.displacedMass;
  finitePos(m, 'mass');
  return (omega) => {
    const c = cylinderHeave({ ...geom, omega });
    return { M: m + c.A, B: c.B + Bvisc, K: c.K, Fe: c.Fe };
  };
}

/**
 * Response statistics and mean absorbed power of a linear absorber with a linear PTO in an irregular sea, by
 * superposition: σ_v² = ∫ ω² |H(ω)|² S(ω) dω with H = F_e / (−ω² M + i ω (B + B_pto) + K + K_pto), and
 * P̄ = B_pto σ_v² (Falnes 2002, ch. 6).
 * @param {object} o
 * @param {number[]|Float64Array} o.omega angular frequencies (rad/s), ascending
 * @param {number[]|Float64Array} o.S wave spectrum (m²·s/rad)
 * @param {(omega:number)=>{M:number, B:number, K:number, Fe:number}} o.body body coefficients
 * @param {number} o.Bpto PTO damping (N·s/m)
 * @param {number} [o.Kpto=0] PTO stiffness (N/m)
 * @returns {{power:number, sigmaX:number, sigmaV:number, sigmaF:number, Tz:number}} mean power (W), standard
 *   deviations of stroke (m), velocity (m/s) and PTO force (N), zero-crossing period of the motion (s)
 */
export function spectralResponse({ omega, S, body, Bpto, Kpto = 0 }) {
  need(omega && S && omega.length === S.length && omega.length > 2, 'omega and S must be arrays of equal length > 2');
  need(typeof body === 'function', 'body must be a function of omega');
  need(Bpto >= 0, 'Bpto must be >= 0');
  const n = omega.length;
  const sx = new Float64Array(n);
  const sv = new Float64Array(n);
  const sf = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const w = omega[i];
    if (!(w > 0) || S[i] === 0) continue;
    const b = body(w);
    const re = b.K + Kpto - w * w * b.M;
    const im = w * (b.B + Bpto);
    const h2 = (b.Fe * b.Fe) / (re * re + im * im);
    sx[i] = h2 * S[i];
    sv[i] = w * w * sx[i];
    sf[i] = (w * w * Bpto * Bpto + Kpto * Kpto) * sx[i];
  }
  const m0 = trapz(sx, omega);
  const m2 = trapz(sv, omega);
  return {
    power: Bpto * m2,
    sigmaX: Math.sqrt(m0),
    sigmaV: Math.sqrt(m2),
    sigmaF: Math.sqrt(trapz(sf, omega)),
    Tz: m2 > 0 ? 2 * Math.PI * Math.sqrt(m0 / m2) : NaN,
  };
}

function defaultOmega() {
  const w = new Float64Array(240);
  for (let i = 0; i < w.length; i++) w[i] = 0.15 + (i * (3.0 - 0.15)) / (w.length - 1);
  return w;
}

/**
 * Power matrix of a linear absorber over (Hs, Tp) with a rating cap.
 * PTO damping per sea state: a fixed value, or 'optimal' — the passive optimum at the peak frequency,
 * B_pto = √(B² + (ω_p M − K/ω_p)²), a common sub-optimal sea-state tuning.
 * @param {object} o
 * @param {number[]} o.hs significant wave heights (m)
 * @param {number[]} o.tp peak periods (s)
 * @param {(omega:number)=>{M:number, B:number, K:number, Fe:number}} o.body
 * @param {number|'optimal'} [o.Bpto='optimal'] PTO damping (N·s/m)
 * @param {number} [o.Kpto=0] PTO stiffness (N/m)
 * @param {number} [o.gamma=1] spectral peak-enhancement factor (-)
 * @param {number} [o.ratedPower=Infinity] cap on mean power per sea state (W)
 * @param {number} [o.efficiency=1] PTO conversion efficiency (-)
 * @param {number[]|Float64Array} [o.omega] frequency grid (rad/s); default 0.15–3.0 rad/s, 240 points
 * @returns {{hs:number[], tp:number[], power:number[][], Bpto:number[]}} power[iHs][iTp] (W), damping per Tp
 */
export function powerMatrix({ hs, tp, body, Bpto = 'optimal', Kpto = 0, gamma = 1, ratedPower = Infinity, efficiency = 1, omega }) {
  need(Array.isArray(hs) && Array.isArray(tp) && hs.length > 0 && tp.length > 0, 'hs and tp must be non-empty arrays');
  need(efficiency > 0 && efficiency <= 1, 'efficiency must be within (0, 1]');
  const w = omega ?? defaultOmega();
  const bUsed = tp.map((T) => {
    if (Bpto !== 'optimal') return Bpto;
    const wp = (2 * Math.PI) / T;
    const b = body(wp);
    return Math.hypot(b.B, wp * b.M - (b.K + Kpto) / wp);
  });
  // response is linear in Hs: compute once per Tp for Hs = 1 and scale by Hs²
  const unit = tp.map((T, j) => spectralResponse({ omega: w, S: waveSpectrum({ Hs: 1, Tp: T, gamma, omega: w }), body, Bpto: bUsed[j], Kpto }).power);
  const power = hs.map((H) => tp.map((_, j) => Math.min(ratedPower, efficiency * unit[j] * H * H)));
  return { hs: hs.slice(), tp: tp.slice(), power, Bpto: bUsed };
}

/**
 * Annual energy from a power matrix and a scatter diagram on the same (Hs, T) grid.
 * @param {object} o
 * @param {number[][]} o.power mean power per sea state power[iHs][iT] (W)
 * @param {number[][]} o.prob occurrence prob[iHs][iT] (normalised internally)
 * @param {number} [o.ratedPower] rated power (W) for the capacity factor; default the matrix maximum
 * @param {number} [o.availability=1]
 * @param {number} [o.arrayFactor=1] interaction factor q = P_array / (N P_isolated) (-)
 * @param {number} [o.devices=1] number of devices N
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{aepWh:number, meanPower:number, capacityFactor:number}} for the whole array
 */
export function matrixEnergy({ power, prob, ratedPower, availability = 1, arrayFactor = 1, devices = 1, hoursPerYear = 8760 }) {
  need(Array.isArray(power) && Array.isArray(prob) && power.length === prob.length, 'power and prob must have the same shape');
  need(availability >= 0 && availability <= 1, 'availability must be within [0, 1]');
  finitePos(arrayFactor, 'arrayFactor');
  let tot = 0;
  let mean = 0;
  let pMax = 0;
  for (let i = 0; i < power.length; i++) {
    need(power[i].length === prob[i].length, 'power and prob must have the same shape');
    for (let j = 0; j < power[i].length; j++) {
      tot += prob[i][j];
      mean += power[i][j] * prob[i][j];
      pMax = Math.max(pMax, power[i][j]);
    }
  }
  finitePos(tot, 'sum of prob');
  const one = (mean / tot) * availability;
  const rated = ratedPower ?? pMax;
  const meanPower = one * arrayFactor * devices;
  return { aepWh: meanPower * hoursPerYear, meanPower, capacityFactor: rated > 0 ? (one * arrayFactor) / rated : NaN };
}

/**
 * Radiation impulse-response kernel from the radiation damping by the cosine transform
 *   K(t) = (2/π) ∫₀^∞ B(ω) cos(ω t) dω   (Cummins 1962; Ogilvie 1964),
 * evaluated by the trapezoidal rule over the supplied frequencies (B is taken as zero beyond the last one, so
 * the table should extend until B has decayed).
 * @param {object} o
 * @param {number[]|Float64Array} o.omega angular frequencies (rad/s), ascending
 * @param {number[]|Float64Array} o.B radiation damping (N·s/m)
 * @param {number} o.dt kernel time step (s)
 * @param {number} o.tMax kernel length (s)
 * @returns {{dt:number, K:Float64Array}} kernel samples K(i dt) (N/m)
 */
export function radiationKernel({ omega, B, dt, tMax }) {
  need(omega && B && omega.length === B.length && omega.length > 2, 'omega and B must be arrays of equal length > 2');
  finitePos(dt, 'dt');
  finitePos(tMax, 'tMax');
  const n = Math.floor(tMax / dt) + 1;
  const K = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const t = i * dt;
    let s = 0;
    for (let j = 1; j < omega.length; j++) {
      s += 0.5 * (B[j] * Math.cos(omega[j] * t) + B[j - 1] * Math.cos(omega[j - 1] * t)) * (omega[j] - omega[j - 1]);
    }
    K[i] = (2 / Math.PI) * s;
  }
  return { dt, K };
}

/**
 * Infinite-frequency added mass from Ogilvie's relation A∞ = A(ω) + (1/ω) ∫₀^∞ K(t) sin(ω t) dt.
 * @param {object} o
 * @param {number} o.A added mass at frequency omega (kg)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {{dt:number, K:Float64Array}} o.kernel radiation kernel
 * @returns {number} A∞ (kg)
 */
export function infiniteAddedMass({ A, omega, kernel }) {
  finitePos(omega, 'omega');
  const { dt, K } = kernel;
  let s = 0;
  for (let i = 1; i < K.length; i++) s += 0.5 * (K[i] * Math.sin(omega * i * dt) + K[i - 1] * Math.sin(omega * (i - 1) * dt)) * dt;
  return A + s / omega;
}

/**
 * Time-domain single-degree-of-freedom absorber with the Cummins equation
 *   (m + A∞) ẍ + ∫₀ᵗ K(t − τ) ẋ(τ) dτ + (K_h + K_m) x = F_exc(t) + F_pto + F_drag + F_endstop + F_extra,
 * integrated with fixed-step RK4; the convolution uses the stored velocity history (trapezoidal rule, the
 * current-velocity term being updated within the RK stages).
 *
 * PTO models (`pto.type`):
 *  - 'linear':   F = −B v − K x.
 *  - 'coulomb':  F = −F_c tanh(v / v_s) (constant-force/rectified hydraulic PTO; v_s smooths the switch).
 *  - 'hydraulic': piston of area A_p pumping through rectifying valves into a gas accumulator feeding a
 *                motor: F = −(p − p_low) A_p tanh(v / v_s); gas volume dV_g/dt = −(A_p |v| − Q_m),
 *                p = p0 (V0/V_g)^γ (polytropic), motor flow Q_m = C_m (p − p_low); generated power
 *                P_gen = η_m Q_m (p − p_low). Valve losses, oil compressibility and low-pressure-side
 *                dynamics are neglected.
 * Control: `latching` holds the body for `latching.duration` each time its velocity changes sign (Budal &
 * Falnes 1980); `declutch.threshold` disconnects the PTO while |v| is below the threshold (Babarit et al. 2009).
 * End stops act as stiff springs beyond ±`endStop.stroke`. Quadratic drag F = −c_d v|v| with c_d = ½ ρ Cd A.
 *
 * @param {object} o
 * @param {number} o.mass body mass (kg)
 * @param {number} o.Ainf infinite-frequency added mass (kg)
 * @param {number} o.stiffness hydrostatic stiffness K_h (N/m)
 * @param {{dt:number, K:Float64Array}|((t:number)=>number)|null} [o.kernel] radiation kernel (N/m); null = none
 * @param {((t:number)=>number)|number[]|Float64Array} [o.excitation] excitation force (N) as a function of time
 *   or samples at the simulation step
 * @param {number} o.dt time step (s)
 * @param {number} o.duration simulated time (s)
 * @param {number} [o.x0=0] initial displacement (m)
 * @param {number} [o.v0=0] initial velocity (m/s)
 * @param {object} [o.pto] PTO definition (see above); default none
 * @param {number} [o.mooringStiffness=0] K_m (N/m)
 * @param {number} [o.dragCoeff=0] c_d = ½ ρ Cd A (kg/m)
 * @param {{stroke:number, stiffness:number}} [o.endStop]
 * @param {{duration:number}} [o.latching]
 * @param {{threshold:number}} [o.declutch]
 * @param {(x:number, v:number, t:number)=>number} [o.extraForce] additional force (N), e.g. nonlinear
 *   Froude–Krylov/hydrostatic correction for a varying waterplane
 * @param {number} [o.kernelLength] truncate the convolution memory to this duration (s)
 * @returns {{t:Float64Array, x:Float64Array, v:Float64Array, ptoForce:Float64Array, pressure:Float64Array|null,
 *   absorbedEnergy:number, meanAbsorbedPower:number, meanGeneratedPower:number, maxStroke:number,
 *   maxVelocity:number, maxPtoForce:number}} energies in J, powers in W (mean over the run)
 */
export function simulateCummins({
  mass,
  Ainf,
  stiffness,
  kernel = null,
  excitation,
  dt,
  duration,
  x0 = 0,
  v0 = 0,
  pto,
  mooringStiffness = 0,
  dragCoeff = 0,
  endStop,
  latching,
  declutch,
  extraForce,
  kernelLength,
}) {
  finitePos(mass, 'mass');
  need(Ainf >= 0, 'Ainf must be >= 0');
  need(stiffness >= 0, 'stiffness must be >= 0');
  finitePos(dt, 'dt');
  finitePos(duration, 'duration');
  const n = Math.round(duration / dt);
  need(n >= 2 && n <= 2e6, 'duration/dt must be between 2 and 2e6 steps');
  const M = mass + Ainf;
  const Kt = stiffness + mooringStiffness;
  // kernel resampled on the simulation step
  let Kk = null;
  if (kernel) {
    let nk;
    if (typeof kernel === 'function') {
      nk = Math.min(n + 1, kernelLength ? Math.floor(kernelLength / dt) + 1 : n + 1);
      Kk = new Float64Array(nk);
      for (let i = 0; i < nk; i++) Kk[i] = kernel(i * dt);
    } else {
      const tEnd = Math.min((kernel.K.length - 1) * kernel.dt, kernelLength ?? Infinity);
      nk = Math.floor(tEnd / dt + 1e-9) + 1;
      Kk = new Float64Array(nk);
      for (let i = 0; i < nk; i++) {
        const s = (i * dt) / kernel.dt;
        const j = Math.min(kernel.K.length - 2, Math.floor(s));
        Kk[i] = kernel.K[j] + (s - j) * (kernel.K[j + 1] - kernel.K[j]);
      }
    }
  }
  const excFn = typeof excitation === 'function';
  const excAt = (i, f) => {
    if (!excitation) return 0;
    if (excFn) return excitation((i + f) * dt);
    const a = excitation[Math.min(i, excitation.length - 1)] ?? 0;
    const b = excitation[Math.min(i + 1, excitation.length - 1)] ?? 0;
    return a + f * (b - a);
  };
  const type = pto?.type ?? 'none';
  need(['none', 'linear', 'coulomb', 'hydraulic'].includes(type), 'pto.type is not a known PTO model');
  const vs = pto?.vSmooth ?? 0.01;
  const hyd = type === 'hydraulic';
  if (hyd) {
    finitePos(pto.pistonArea, 'pto.pistonArea');
    finitePos(pto.p0, 'pto.p0');
    finitePos(pto.V0, 'pto.V0');
    need(pto.motorCoeff >= 0, 'pto.motorCoeff must be >= 0');
  }
  const gam = pto?.gamma ?? 1.4;
  const pLow = pto?.pLow ?? 0;
  const pressureOf = (Vg) => pto.p0 * Math.pow(pto.V0 / Vg, gam);
  const ptoForce = (x, v, Vg) => {
    if (type === 'none') return 0;
    if (declutch && Math.abs(v) < declutch.threshold) return 0;
    if (type === 'linear') return -(pto.B ?? 0) * v - (pto.K ?? 0) * x;
    if (type === 'coulomb') return -pto.force * Math.tanh(v / vs);
    return -(pressureOf(Vg) - pLow) * pto.pistonArea * Math.tanh(v / vs);
  };
  const t = new Float64Array(n + 1);
  const xs = new Float64Array(n + 1);
  const vsr = new Float64Array(n + 1);
  const fp = new Float64Array(n + 1);
  const pr = hyd ? new Float64Array(n + 1) : null;
  let x = x0;
  let v = v0;
  let Vg = hyd ? pto.V0 : 1;
  let latchedUntil = -1;
  let eAbs = 0;
  let eGen = 0;
  let maxX = Math.abs(x0);
  let maxV = Math.abs(v0);
  let maxF = 0;
  xs[0] = x;
  vsr[0] = v;
  fp[0] = ptoForce(x, v, Vg);
  if (hyd) pr[0] = pressureOf(Vg);
  const half0 = Kk ? 0.5 * Kk[0] * dt : 0;
  const acc = (xx, vv, VV, tt, fexc, mem) => {
    let f = fexc - mem - half0 * vv - Kt * xx + ptoForce(xx, vv, VV);
    if (dragCoeff) f -= dragCoeff * vv * Math.abs(vv);
    if (endStop && Math.abs(xx) > endStop.stroke) f -= endStop.stiffness * (Math.abs(xx) - endStop.stroke) * Math.sign(xx);
    if (extraForce) f += extraForce(xx, vv, tt);
    return f / M;
  };
  const dVg = (vv, VV) => (hyd ? -(pto.pistonArea * Math.abs(vv) - pto.motorCoeff * (pressureOf(VV) - pLow)) : 0);
  for (let i = 0; i < n; i++) {
    const time = i * dt;
    if (time < latchedUntil) {
      // body held: zero velocity, position frozen; the accumulator keeps discharging through the motor
      if (hyd) Vg += dVg(0, Vg) * dt;
      v = 0;
    } else {
      // memory term from past velocities (excluding the current one, handled inside the stages)
      let mem = 0;
      if (Kk) {
        const m = Math.min(i, Kk.length - 1);
        for (let j = 1; j <= m; j++) mem += Kk[j] * vsr[i - j];
        if (m === i && i > 0) mem -= 0.5 * Kk[m] * vsr[0];
        mem *= dt;
      }
      const f0 = excAt(i, 0);
      const f1 = excAt(i, 0.5);
      const f2 = excAt(i, 1);
      const k1x = v;
      const k1v = acc(x, v, Vg, time, f0, mem);
      const k1g = dVg(v, Vg);
      const k2x = v + 0.5 * dt * k1v;
      const k2v = acc(x + 0.5 * dt * k1x, k2x, Vg + 0.5 * dt * k1g, time + 0.5 * dt, f1, mem);
      const k2g = dVg(k2x, Vg + 0.5 * dt * k1g);
      const k3x = v + 0.5 * dt * k2v;
      const k3v = acc(x + 0.5 * dt * k2x, k3x, Vg + 0.5 * dt * k2g, time + 0.5 * dt, f1, mem);
      const k3g = dVg(k3x, Vg + 0.5 * dt * k2g);
      const k4x = v + dt * k3v;
      const k4v = acc(x + dt * k3x, k4x, Vg + dt * k3g, time + dt, f2, mem);
      const k4g = dVg(k4x, Vg + dt * k3g);
      const vOld = v;
      x += (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x);
      v += (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v);
      if (hyd) Vg = Math.max(1e-6 * pto.V0, Vg + (dt / 6) * (k1g + 2 * k2g + 2 * k3g + k4g));
      if (latching && vOld !== 0 && vOld * v < 0) {
        v = 0;
        latchedUntil = time + dt + latching.duration;
      }
    }
    t[i + 1] = (i + 1) * dt;
    xs[i + 1] = x;
    vsr[i + 1] = v;
    const F = ptoForce(x, v, Vg);
    fp[i + 1] = F;
    if (hyd) {
      const p = pressureOf(Vg);
      pr[i + 1] = p;
      eGen += (pto.motorEfficiency ?? 1) * pto.motorCoeff * (p - pLow) ** 2 * dt;
    }
    eAbs += -0.5 * (fp[i] * vsr[i] + F * v) * dt;
    maxX = Math.max(maxX, Math.abs(x));
    maxV = Math.max(maxV, Math.abs(v));
    maxF = Math.max(maxF, Math.abs(F));
  }
  const T = n * dt;
  return {
    t,
    x: xs,
    v: vsr,
    ptoForce: fp,
    pressure: pr,
    absorbedEnergy: eAbs,
    meanAbsorbedPower: eAbs / T,
    meanGeneratedPower: hyd ? eGen / T : (eAbs / T) * (pto?.efficiency ?? 1),
    maxStroke: maxX,
    maxVelocity: maxV,
    maxPtoForce: maxF,
  };
}

function cdiv(ar, ai, br, bi) {
  const d = br * br + bi * bi;
  return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d];
}

/**
 * Two-body heaving absorber (float 1 reacting against a spar/plate 2) with a PTO acting on the relative motion,
 * frequency domain (Falnes 1999, "Wave-energy conversion through relative motion between two single-mode
 * oscillating bodies", J. Offshore Mech. Arct. Eng. 121:32–38):
 *   [Z1 + Zp] X1 + [Zc − Zp] X2 = F1 ζ
 *   [Zc − Zp] X1 + [Z2 + Zp] X2 = F2 ζ
 * Z_j = −ω² (m_j + A_jj) + i ω B_jj + K_j,  Zc = −ω² A12 + i ω B12,  Zp = i ω B_pto + K_pto,
 * P = ½ B_pto ω² |X1 − X2|². Excitation coefficients may be complex ({ re, im }) to carry the phase between
 * the two bodies. With `optimise` the PTO damping maximising power is found by a bounded search.
 * @param {object} o
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.amplitude wave amplitude (m)
 * @param {{mass:number, A:number, B:number, K:number, Fe:number|{re:number, im:number}}} o.body1 float
 * @param {{mass:number, A:number, B:number, K:number, Fe:number|{re:number, im:number}}} o.body2 reaction body
 * @param {number} [o.A12=0] coupling added mass (kg)
 * @param {number} [o.B12=0] coupling radiation damping (N·s/m)
 * @param {number} [o.Bpto] PTO damping (N·s/m); required unless `optimise`
 * @param {number} [o.Kpto=0] PTO stiffness (N/m)
 * @param {boolean} [o.optimise=false]
 * @returns {{power:number, Bpto:number, x1:number, x2:number, relative:number, ptoForce:number}} amplitudes (m),
 *   PTO force amplitude (N)
 */
export function twoBodyAbsorber({ omega, amplitude, body1, body2, A12 = 0, B12 = 0, Bpto, Kpto = 0, optimise = false }) {
  finitePos(omega, 'omega');
  need(amplitude >= 0, 'amplitude must be >= 0');
  const cx = (f) => (typeof f === 'number' ? [f, 0] : [f.re, f.im]);
  const w2 = omega * omega;
  const z1 = [body1.K - w2 * (body1.mass + body1.A), omega * body1.B];
  const z2 = [body2.K - w2 * (body2.mass + body2.A), omega * body2.B];
  const zc = [-w2 * A12, omega * B12];
  const f1 = cx(body1.Fe).map((v) => v * amplitude);
  const f2 = cx(body2.Fe).map((v) => v * amplitude);
  const solve = (bp) => {
    const zp = [Kpto, omega * bp];
    const a = [z1[0] + zp[0], z1[1] + zp[1]];
    const b = [zc[0] - zp[0], zc[1] - zp[1]];
    const d = [z2[0] + zp[0], z2[1] + zp[1]];
    // determinant a d − b²
    const det = [a[0] * d[0] - a[1] * d[1] - (b[0] * b[0] - b[1] * b[1]), a[0] * d[1] + a[1] * d[0] - 2 * b[0] * b[1]];
    const n1 = [f1[0] * d[0] - f1[1] * d[1] - (b[0] * f2[0] - b[1] * f2[1]), f1[0] * d[1] + f1[1] * d[0] - (b[0] * f2[1] + b[1] * f2[0])];
    const n2 = [a[0] * f2[0] - a[1] * f2[1] - (b[0] * f1[0] - b[1] * f1[1]), a[0] * f2[1] + a[1] * f2[0] - (b[0] * f1[1] + b[1] * f1[0])];
    const X1 = cdiv(n1[0], n1[1], det[0], det[1]);
    const X2 = cdiv(n2[0], n2[1], det[0], det[1]);
    const rel = Math.hypot(X1[0] - X2[0], X1[1] - X2[1]);
    return { power: 0.5 * bp * w2 * rel * rel, Bpto: bp, x1: Math.hypot(X1[0], X1[1]), x2: Math.hypot(X2[0], X2[1]), relative: rel, ptoForce: rel * Math.hypot(omega * bp, Kpto) };
  };
  if (!optimise) {
    need(Number.isFinite(Bpto) && Bpto >= 0, 'Bpto must be >= 0 (or set optimise)');
    return solve(Bpto);
  }
  const ref = Math.max(body1.B + body2.B, omega * (body1.mass + body1.A) * 1e-3, 1);
  let best = null;
  let li = 0;
  const nS = 400;
  const at = (i) => ref * Math.exp(Math.log(1e-4) + (Math.log(1e8) * i) / nS);
  for (let i = 0; i <= nS; i++) {
    const r = solve(at(i));
    if (!best || r.power > best.power) {
      best = r;
      li = i;
    }
  }
  let a = Math.log(at(Math.max(0, li - 1)));
  let b = Math.log(at(Math.min(nS, li + 1)));
  const gr = (Math.sqrt(5) - 1) / 2;
  for (let i = 0; i < 60; i++) {
    const c = b - gr * (b - a);
    const d = a + gr * (b - a);
    if (solve(Math.exp(c)).power > solve(Math.exp(d)).power) b = d;
    else a = c;
  }
  const r = solve(Math.exp(0.5 * (a + b)));
  return r.power > best.power ? r : best;
}

/**
 * Bottom-hinged flap (oscillating wave surge converter), two-dimensional thin-flap hydrodynamics per unit width.
 * A flap hinged at the seabed and spanning the water column is a flap-type wavemaker; its far-field wave
 * amplitude per unit surface stroke amplitude is (Dean & Dalrymple 1991, Water Wave Mechanics for Engineers and Scientists, ch. 6)
 *   R_f = 4 (sinh kh / kh) (kh sinh kh − cosh kh + 1) / (sinh 2kh + 2kh).
 * Energy radiated to both sides gives the pitch radiation damping per unit width
 *   b = 2 ρ g c_g (h R_f)² / ω²,
 * and the two-dimensional Haskind relation b = (|X⁺|² + |X⁻|²)/(4 ρ g c_g) with |X⁺| = |X⁻| gives the
 * excitation moment per unit width and unit wave amplitude
 *   |X| = 2 ρ g c_g h R_f / ω.
 * The corresponding maximum absorbed power is |X ζ|²/(8 b) = J/2 per unit width (the antisymmetric-mode limit).
 * Valid for a thin, surface-piercing flap hinged at the bed in linear waves; finite-width (3-D) effects enter
 * through the caller's width factor.
 * @param {object} o
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.depth water depth h (m)
 * @param {number} [o.rho=1025]
 * @returns {{excitation:number, damping:number, waveRatio:number, k:number, cg:number}} |X| (N·m per m width
 *   per m amplitude = N/m), b (N·m·s/rad per m width), R_f (-)
 */
export function flapHydrodynamics2D({ omega, depth, rho = 1025 }) {
  finitePos(depth, 'depth');
  need(Number.isFinite(depth), 'depth must be finite for a bottom-hinged flap');
  const { k, cg } = dispersion(omega, depth);
  const kh = k * depth;
  need(kh < 300, 'k·depth is too large for the flap model');
  const Rf = (4 * (Math.sinh(kh) / kh) * (kh * Math.sinh(kh) - Math.cosh(kh) + 1)) / (Math.sinh(2 * kh) + 2 * kh);
  return {
    excitation: (2 * rho * G * cg * depth * Rf) / omega,
    damping: (2 * rho * G * cg * (depth * Rf) ** 2) / (omega * omega),
    waveRatio: Rf,
    k,
    cg,
  };
}

/**
 * Hydrostatic pitch stiffness of a buoyant flap about its hinge, C = g (ρ V z_b − m z_g): buoyancy moment minus
 * weight moment for small rotations (thin flap, waterplane contribution neglected).
 * @param {object} o
 * @param {number} o.volume displaced volume V (m³)
 * @param {number} o.buoyancyArm height of the centre of buoyancy above the hinge z_b (m)
 * @param {number} o.mass flap mass m (kg)
 * @param {number} o.gravityArm height of the centre of gravity above the hinge z_g (m)
 * @param {number} [o.rho=1025]
 * @returns {number} C (N·m/rad)
 */
export function flapStiffness({ volume, buoyancyArm, mass, gravityArm, rho = 1025 }) {
  need(volume >= 0 && mass >= 0, 'volume and mass must be >= 0');
  return G * (rho * volume * buoyancyArm - mass * gravityArm);
}

/**
 * Bottom-hinged flap converter, one-degree-of-freedom pitch model in a regular wave:
 *   [−ω² (I + I_a) + i ω (B + B_v + B_pto) + C] Θ = M_e ζ,  P = ½ B_pto ω² |Θ|²,
 * with M_e = w f_w |X| and B = w f_w² b from `flapHydrodynamics2D` (w = flap width, f_w a 3-D width factor).
 * The added inertia I_a is an argument (it depends on the evanescent near field and the flap thickness).
 * @param {object} o
 * @param {number} o.width flap width w (m)
 * @param {number} o.depth water depth (m)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.amplitude wave amplitude (m)
 * @param {number} o.inertia flap pitch inertia about the hinge I (kg·m²)
 * @param {number} o.addedInertia added pitch inertia I_a (kg·m²)
 * @param {number} o.stiffness hydrostatic pitch stiffness C (N·m/rad)
 * @param {number|'optimal'} [o.Bpto='optimal'] PTO damping (N·m·s/rad); 'optimal' = passive optimum
 * @param {number} [o.Bvisc=0] linearised viscous damping (N·m·s/rad)
 * @param {number} [o.widthFactor=1] 3-D correction on the excitation (-)
 * @param {number} [o.maxAngle] pitch amplitude limit (rad)
 * @param {number} [o.rho=1025]
 * @returns {{power:number, angle:number, Bpto:number, ptoTorque:number, captureWidth:number,
 *   captureWidthRatio:number, excitationMoment:number, radiationDamping:number, constrained:boolean}}
 */
export function flapConverter({ width, depth, omega, amplitude, inertia, addedInertia, stiffness, Bpto = 'optimal', Bvisc = 0, widthFactor = 1, maxAngle, rho = 1025 }) {
  finitePos(width, 'width');
  finitePos(inertia, 'inertia');
  need(addedInertia >= 0, 'addedInertia must be >= 0');
  const h2 = flapHydrodynamics2D({ omega, depth, rho });
  const Me = width * widthFactor * h2.excitation;
  const B = width * widthFactor * widthFactor * h2.damping;
  const r = pointAbsorber({
    mass: inertia,
    A: addedInertia,
    B,
    K: stiffness,
    Fe: Me,
    omega,
    amplitude,
    control: Bpto === 'optimal' ? 'passive' : 'fixed',
    Bpto: Bpto === 'optimal' ? undefined : Bpto,
    Bvisc,
    limits: maxAngle !== undefined ? { stroke: maxAngle } : {},
  });
  const J = regularWavePower({ amplitude, omega, depth, rho });
  return {
    power: r.power,
    angle: r.stroke,
    Bpto: r.Bpto,
    ptoTorque: r.ptoForce,
    captureWidth: J > 0 ? r.power / J : NaN,
    captureWidthRatio: J > 0 ? r.power / J / width : NaN,
    excitationMoment: Me,
    radiationDamping: B,
    constrained: r.constrained,
  };
}

/**
 * Parameterised device (attenuator, terminator or any device described by a capture-width ratio):
 * P = min(P_rated, η · CWR · L · J), with the capture-width ratio CWR and characteristic length L supplied by
 * the user (e.g. from tank tests; Babarit 2015 reviews measured CWR by device class). Not a hydrodynamic model.
 * @param {object} o
 * @param {number} o.J wave power per unit crest (W/m)
 * @param {number} o.captureWidthRatio CWR (-)
 * @param {number} o.length characteristic length L (m)
 * @param {number} [o.efficiency=1] PTO/conversion efficiency (-)
 * @param {number} [o.ratedPower=Infinity] rating (W)
 * @returns {number} power (W)
 */
export function captureWidthDevice({ J, captureWidthRatio, length, efficiency = 1, ratedPower = Infinity }) {
  need(J >= 0, 'J must be >= 0');
  need(captureWidthRatio >= 0, 'captureWidthRatio must be >= 0');
  finitePos(length, 'length');
  return Math.min(ratedPower, efficiency * captureWidthRatio * length * J);
}

/**
 * Overtopping device: hydraulic power of the overtopping discharge stored at the crest freeboard,
 * P = η ρ g q R_c w, with the mean overtopping rate q per unit width supplied by the user (from an overtopping
 * formula such as EurOtop or from tests). Not a wave model.
 * @param {object} o
 * @param {number} o.q mean overtopping discharge per metre width (m³/s/m)
 * @param {number} o.crestFreeboard reservoir head above mean sea level R_c (m)
 * @param {number} o.width ramp width w (m)
 * @param {number} [o.efficiency=1] turbine–generator efficiency (-)
 * @param {number} [o.rho=1025]
 * @returns {number} power (W)
 */
export function overtoppingPower({ q, crestFreeboard, width, efficiency = 1, rho = 1025 }) {
  need(q >= 0, 'q must be >= 0');
  need(crestFreeboard >= 0, 'crestFreeboard must be >= 0');
  finitePos(width, 'width');
  return efficiency * rho * G * q * crestFreeboard * width;
}

/**
 * Array power with an interaction factor q = P_array / (N P_isolated) supplied by the user (q < 1 destructive,
 * q > 1 constructive; for point-absorber arrays averaged over directions q ≈ 1, Babarit 2013).
 * @param {object} o
 * @param {number} o.devices number of devices N
 * @param {number} o.power isolated-device power (W)
 * @param {number} [o.q=1] interaction factor (-)
 * @returns {number} array power (W)
 */
export function arrayPower({ devices, power, q = 1 }) {
  need(Number.isInteger(devices) && devices >= 1, 'devices must be an integer >= 1');
  finitePos(q, 'q');
  return devices * power * q;
}

/**
 * Survival screening in an extreme sea state from the linear spectral response: most probable maximum of
 * stroke, velocity and PTO force in a storm of given duration, x_mpm = σ √(2 ln N), N = duration / T_z
 * (Rayleigh-distributed maxima of a narrow-band Gaussian process). Linear theory overstates resonant response in
 * extreme seas; treat as an upper screening value and check end-stop demand.
 * @param {object} o
 * @param {number[]|Float64Array} o.omega (rad/s)
 * @param {number[]|Float64Array} o.S extreme sea-state spectrum (m²·s/rad)
 * @param {(omega:number)=>{M:number, B:number, K:number, Fe:number}} o.body
 * @param {number} o.Bpto PTO damping in survival mode (N·s/m)
 * @param {number} [o.Kpto=0]
 * @param {number} [o.duration=10800] storm duration (s)
 * @param {number} [o.strokeLimit] available stroke (m)
 * @returns {{maxStroke:number, maxVelocity:number, maxPtoForce:number, sigmaX:number, cycles:number,
 *   strokeUtilisation:number}}
 */
export function survivalResponse({ omega, S, body, Bpto, Kpto = 0, duration = 10800, strokeLimit }) {
  finitePos(duration, 'duration');
  const r = spectralResponse({ omega, S, body, Bpto, Kpto });
  const N = Math.max(2, duration / r.Tz);
  const f = Math.sqrt(2 * Math.log(N));
  return {
    maxStroke: r.sigmaX * f,
    maxVelocity: r.sigmaV * f,
    maxPtoForce: r.sigmaF * f,
    sigmaX: r.sigmaX,
    cycles: N,
    strokeUtilisation: strokeLimit > 0 ? (r.sigmaX * f) / strokeLimit : NaN,
  };
}

/**
 * Slamming force hook, F = ½ ρ C_s A v² with the slamming coefficient supplied by the user
 * (DNV-RP-C205 §8.6 gives guidance by geometry).
 * @param {object} o
 * @param {number} o.Cs slamming coefficient (-)
 * @param {number} o.area projected area exposed to the slam (m²)
 * @param {number} o.velocity relative impact velocity (m/s)
 * @param {number} [o.rho=1025]
 * @returns {number} force (N)
 */
export function slammingForce({ Cs, area, velocity, rho = 1025 }) {
  need(Cs >= 0 && area >= 0, 'Cs and area must be >= 0');
  return 0.5 * rho * Cs * area * velocity * velocity;
}
