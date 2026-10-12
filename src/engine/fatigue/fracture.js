// Linear-elastic fracture mechanics, fatigue crack growth and failure-assessment-diagram checks.
//
// SI throughout: lengths in m, stresses in Pa, stress intensity in Pa·√m, growth rates in m/cycle.
// Unit trap: BS 7910 publishes Paris constants for ΔK in N/mm^(3/2) and da/dN in mm/cycle; they are converted by
// parisConstantToSI (1 MPa·√m = 31.623 N/mm^(3/2)).

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function positive(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

const N_MM32_PER_PA_SQRT_M = 1e-6 * Math.sqrt(1000);

// ---------------------------------------------------------------------------------------------------------------
// Stress-intensity solutions
// ---------------------------------------------------------------------------------------------------------------

/**
 * Centre through crack of length 2a in a plate of width W under remote tension:
 * K = σ·√(πa)·F, F = [1 − 0.025·α² + 0.06·α⁴]·√(sec(πα/2)), α = 2a/W.
 * Reference: Tada, Paris & Irwin, The Stress Analysis of Cracks Handbook (3rd ed.), centre-cracked strip (the
 * handbook quotes 0.1 % accuracy for any α).
 * @param {object} o
 * @param {number} o.a half crack length [m]
 * @param {number} o.W plate width [m]
 * @param {number} [o.stress=1] remote membrane stress [Pa]
 * @returns {{K:number, F:number}} K [Pa·√m], geometry factor F [-]
 */
export function centreCrackK({ a, W, stress = 1 }) {
  positive(a, 'a');
  positive(W, 'W');
  const al = (2 * a) / W;
  need(al < 1, 'a must be less than W/2');
  const F = (1 - 0.025 * al * al + 0.06 * al ** 4) / Math.sqrt(Math.cos((Math.PI * al) / 2));
  return { K: stress * Math.sqrt(Math.PI * a) * F, F };
}

/**
 * Single edge crack of depth a in a plate of width W under tension and/or bending:
 *   tension: F = √(2/(πα)·tan(πα/2))·[0.752 + 2.02·α + 0.37·(1 − sin(πα/2))³]/cos(πα/2)
 *   bending: F = √(2/(πα)·tan(πα/2))·[0.923 + 0.199·(1 − sin(πα/2))⁴]/cos(πα/2),  α = a/W
 * Reference: Tada, Paris & Irwin, The Stress Analysis of Cracks Handbook (3rd ed.), single edge notch strip. The
 * handbook quotes 0.5 % accuracy for any α; against the Brown & Srawley polynomial fits (a/W <= 0.6) the tension
 * factor agrees within 1 % and the bending factor within 2.5 %.
 * @param {object} o
 * @param {number} o.a crack depth [m]
 * @param {number} o.W plate width (ligament direction) [m]
 * @param {number} [o.tension=0] remote membrane stress [Pa]
 * @param {number} [o.bending=0] outer-fibre bending stress on the cracked side [Pa]
 * @returns {{K:number, Ft:number, Fb:number}} K [Pa·√m], tension and bending geometry factors [-]
 */
export function edgeCrackK({ a, W, tension = 0, bending = 0 }) {
  positive(a, 'a');
  positive(W, 'W');
  const al = a / W;
  need(al < 1, 'a must be less than W');
  const th = (Math.PI * al) / 2;
  const pre = Math.sqrt(Math.tan(th) / th) / Math.cos(th);
  const Ft = pre * (0.752 + 2.02 * al + 0.37 * (1 - Math.sin(th)) ** 3);
  const Fb = pre * (0.923 + 0.199 * (1 - Math.sin(th)) ** 4);
  return { K: (tension * Ft + bending * Fb) * Math.sqrt(Math.PI * a), Ft, Fb };
}

/**
 * Double edge cracks of depth a each in a plate of width W under tension:
 * F = (1.122 − 0.561·α − 0.205·α² + 0.471·α³ − 0.190·α⁴)/√(1 − α), α = 2a/W.
 * Reference: Tada, Paris & Irwin, double edge notch strip (the handbook quotes 0.5 % accuracy for any α).
 * @param {object} o
 * @param {number} o.a depth of each crack [m]
 * @param {number} o.W plate width [m]
 * @param {number} [o.stress=1] [Pa]
 * @returns {{K:number, F:number}}
 */
export function doubleEdgeCrackK({ a, W, stress = 1 }) {
  positive(a, 'a');
  positive(W, 'W');
  const al = (2 * a) / W;
  need(al < 1, 'a must be less than W/2');
  const F = (1.122 - 0.561 * al - 0.205 * al * al + 0.471 * al ** 3 - 0.190 * al ** 4) / Math.sqrt(1 - al);
  return { K: stress * Math.sqrt(Math.PI * a) * F, F };
}

/**
 * Complete elliptic integral of the second kind E(k) by the arithmetic–geometric mean.
 * @param {number} k modulus, 0 <= k <= 1 [-]
 * @returns {number} E(k) [-]
 */
export function ellipticE(k) {
  need(k >= 0 && k <= 1, 'k must be in [0, 1]');
  if (k === 1) return 1;
  let a = 1;
  let b = Math.sqrt(1 - k * k);
  let c = k;
  let sum = 0.5 * c * c;
  let pow = 0.5;
  for (let i = 0; i < 40 && Math.abs(c) > 1e-17; i++) {
    const an = 0.5 * (a + b);
    c = 0.5 * (a - b);
    b = Math.sqrt(a * b);
    a = an;
    pow *= 2;
    sum += pow * c * c;
  }
  return (Math.PI / (2 * a)) * (1 - sum);
}

/**
 * Embedded elliptical crack (semi-axes a <= c) in an infinite body under uniform tension normal to the crack:
 * K(φ) = σ·√(πa)/E(k)·[sin²φ + (a/c)²·cos²φ]^(1/4), k² = 1 − (a/c)²; φ = π/2 is the end of the minor axis.
 * Reference: Irwin (1962) J. Appl. Mech. 29, 651–654 (exact). No finite-thickness correction.
 * @param {object} o
 * @param {number} o.a semi-minor axis [m]
 * @param {number} o.c semi-major axis [m]
 * @param {number} [o.phi=Math.PI/2] parametric angle [rad]
 * @param {number} [o.stress=1] [Pa]
 * @returns {number} K [Pa·√m]
 */
export function embeddedEllipticalCrackK({ a, c, phi = Math.PI / 2, stress = 1 }) {
  positive(a, 'a');
  positive(c, 'c');
  need(a <= c, 'a must not exceed c');
  const r = a / c;
  const Ek = ellipticE(Math.sqrt(1 - r * r));
  return ((stress * Math.sqrt(Math.PI * a)) / Ek) * (Math.sin(phi) ** 2 + r * r * Math.cos(phi) ** 2) ** 0.25;
}

/**
 * Newman–Raju boundary-correction factors for a semi-elliptical surface crack in a finite plate under membrane
 * and through-thickness bending stress: K = (σm + H·σb)·√(πa/Q)·F.
 * For a/c <= 1:
 *   Q = 1 + 1.464·(a/c)^1.65,  F = [M1 + M2·(a/t)² + M3·(a/t)⁴]·fφ·g·fw
 *   M1 = 1.13 − 0.09·(a/c),  M2 = −0.54 + 0.89/(0.2 + a/c),  M3 = 0.5 − 1/(0.65 + a/c) + 14·(1 − a/c)^24
 *   g = 1 + [0.1 + 0.35·(a/t)²]·(1 − sin φ)²,  fφ = [(a/c)²·cos²φ + sin²φ]^(1/4)
 *   fw = [sec(π·c/(2b)·√(a/t))]^(1/2)
 *   H = H1 + (H2 − H1)·sin^p φ,  p = 0.2 + a/c + 0.6·a/t
 *   H1 = 1 − 0.34·a/t − 0.11·(a/c)·(a/t),  H2 = 1 + G21·(a/t) + G22·(a/t)²
 *   G21 = −1.22 − 0.12·(a/c),  G22 = 0.55 − 1.05·(a/c)^0.75 + 0.47·(a/c)^1.5
 * For 1 < a/c <= 2 the Newman–Raju (1984) extension is used:
 *   Q = 1 + 1.464·(c/a)^1.65, M1 = √(c/a)·(1 + 0.04·c/a), M2 = 0.2·(c/a)⁴, M3 = −0.11·(c/a)⁴,
 *   g = 1 + [0.1 + 0.35·(c/a)·(a/t)²]·(1 − sin φ)²,  fφ = [(c/a)²·sin²φ + cos²φ]^(1/4),
 *   p = 0.2 + c/a + 0.6·a/t,  H1 = 1 − (0.04 + 0.41·c/a)·a/t + (0.55 − 1.93·(c/a)^0.75 + 1.38·(c/a)^1.5)·(a/t)²,
 *   G21 = −2.11 + 0.77·c/a,  G22 = 0.55 − 0.72·(c/a)^0.75 + 0.14·(c/a)^1.5.
 * Reference: Newman & Raju (1981) Engineering Fracture Mechanics 15(1–2), 185–192; Newman & Raju (1984) NASA
 * TM-85793 (a/c > 1); the same expressions are reproduced in BS 7910 Annex M.
 * Validity: 0 < a/c <= 2, 0 <= a/t < 1 (bending accuracy degrades beyond a/t = 0.8), c/b < 0.5, 0 <= φ <= π.
 * φ = π/2 is the deepest point, φ = 0 the surface point.
 * @param {object} o
 * @param {number} o.a crack depth [m]
 * @param {number} o.c half surface length [m]
 * @param {number} o.t plate thickness [m]
 * @param {number} [o.b=Infinity] plate half-width [m]
 * @param {number} [o.phi=Math.PI/2] parametric angle [rad]
 * @returns {{Q:number, F:number, H:number, Mm:number, Mb:number, fw:number}} Mm = F/√Q and Mb = H·F/√Q are the
 *   BS 7910-style factors on σ·√(πa); all dimensionless
 */
export function newmanRaju({ a, c, t, b = Infinity, phi = Math.PI / 2 }) {
  positive(a, 'a');
  positive(c, 'c');
  positive(t, 't');
  need(b > 0, 'b must be positive');
  need(phi >= 0 && phi <= Math.PI, 'phi must be in [0, pi]');
  const at = a / t;
  const ac = a / c;
  need(at <= 1, 'a must not exceed t');
  need(ac <= 2, 'a/c must not exceed 2');
  need(c < b, 'c must be less than the half-width b');
  const s = Math.sin(phi);
  const co = Math.cos(phi);
  let Q, M1, M2, M3, g, fphi, p, H1, G21, G22;
  if (ac <= 1) {
    Q = 1 + 1.464 * ac ** 1.65;
    M1 = 1.13 - 0.09 * ac;
    M2 = -0.54 + 0.89 / (0.2 + ac);
    M3 = 0.5 - 1 / (0.65 + ac) + 14 * (1 - ac) ** 24;
    g = 1 + (0.1 + 0.35 * at * at) * (1 - s) ** 2;
    fphi = (ac * ac * co * co + s * s) ** 0.25;
    p = 0.2 + ac + 0.6 * at;
    H1 = 1 - 0.34 * at - 0.11 * ac * at;
    G21 = -1.22 - 0.12 * ac;
    G22 = 0.55 - 1.05 * ac ** 0.75 + 0.47 * ac ** 1.5;
  } else {
    const ca = c / a;
    Q = 1 + 1.464 * ca ** 1.65;
    M1 = Math.sqrt(ca) * (1 + 0.04 * ca);
    M2 = 0.2 * ca ** 4;
    M3 = -0.11 * ca ** 4;
    g = 1 + (0.1 + 0.35 * ca * at * at) * (1 - s) ** 2;
    fphi = (ca * ca * s * s + co * co) ** 0.25;
    p = 0.2 + ca + 0.6 * at;
    H1 = 1 - (0.04 + 0.41 * ca) * at + (0.55 - 1.93 * ca ** 0.75 + 1.38 * ca ** 1.5) * at * at;
    G21 = -2.11 + 0.77 * ca;
    G22 = 0.55 - 0.72 * ca ** 0.75 + 0.14 * ca ** 1.5;
  }
  const arg = ((Math.PI * c) / (2 * b)) * Math.sqrt(at);
  const fw = Number.isFinite(b) ? 1 / Math.sqrt(Math.cos(arg)) : 1;
  const F = (M1 + M2 * at * at + M3 * at ** 4) * fphi * g * fw;
  const H2 = 1 + G21 * at + G22 * at * at;
  const H = H1 + (H2 - H1) * s ** p;
  const rq = Math.sqrt(Q);
  return { Q, F, H, Mm: F / rq, Mb: (H * F) / rq, fw };
}

/**
 * Stress intensity at the deepest and surface points of a semi-elliptical surface crack in a plate, with
 * optional weld-toe magnification factors: K = (Mkm·Mm·σm + Mkb·Mb·σb)·√(πa).
 * Reference: Newman & Raju (1981); BS 7910:2019 Annex M (Mk supplied by the caller).
 * @param {object} o
 * @param {number} o.a crack depth [m]
 * @param {number} o.c half surface length [m]
 * @param {number} o.t thickness [m]
 * @param {number} [o.W=Infinity] full plate width [m]
 * @param {number} [o.membrane=0] membrane stress [Pa]
 * @param {number} [o.bending=0] bending stress at the cracked surface [Pa]
 * @param {{deepMembrane?:number, deepBending?:number, surfaceMembrane?:number, surfaceBending?:number}} [o.Mk]
 *   weld-toe magnification factors [-], default 1
 * @returns {{Kdeep:number, Ksurface:number, deep:object, surface:object}} K [Pa·√m]; factor sets from newmanRaju
 */
export function surfaceCrackK({ a, c, t, W = Infinity, membrane = 0, bending = 0, Mk = {} }) {
  const b = W / 2;
  const deep = newmanRaju({ a, c, t, b, phi: Math.PI / 2 });
  const surface = newmanRaju({ a, c, t, b, phi: 0 });
  const root = Math.sqrt(Math.PI * a);
  const { deepMembrane = 1, deepBending = 1, surfaceMembrane = 1, surfaceBending = 1 } = Mk;
  return {
    Kdeep: (deepMembrane * deep.Mm * membrane + deepBending * deep.Mb * bending) * root,
    Ksurface: (surfaceMembrane * surface.Mm * membrane + surfaceBending * surface.Mb * bending) * root,
    deep, surface,
  };
}

// ---------------------------------------------------------------------------------------------------------------
// Fatigue crack growth
// ---------------------------------------------------------------------------------------------------------------

/**
 * Convert a Paris coefficient published for ΔK in N/mm^(3/2) and da/dN in mm/cycle to SI (Pa·√m, m/cycle).
 * @param {number} A coefficient in N–mm units
 * @param {number} m exponent [-]
 * @returns {number} coefficient in SI units [m/cycle/(Pa·√m)^m]
 */
export function parisConstantToSI(A, m) {
  positive(A, 'A');
  positive(m, 'm');
  return 1e-3 * A * N_MM32_PER_PA_SQRT_M ** m;
}

/**
 * Convert a stress-intensity value between N/mm^(3/2) and Pa·√m.
 * @param {number} k value in N/mm^(3/2)
 * @returns {number} value in Pa·√m
 */
export function kFromNmm32(k) {
  return k / N_MM32_PER_PA_SQRT_M;
}

// [A stage A, m stage A, A stage B, m stage B] in N–mm units; the stage A/B transition is the intersection.
const BS7910_TWO_STAGE = {
  'air|low|mean': [1.21e-26, 8.16, 3.98e-13, 2.88],
  'air|low|mean+2sd': [4.37e-26, 8.16, 6.77e-13, 2.88],
  'air|high|mean': [4.80e-18, 5.10, 5.86e-13, 2.88],
  'air|high|mean+2sd': [2.10e-17, 5.10, 1.29e-12, 2.88],
  'marine-free|low|mean': [3.0e-14, 3.42, 1.27e-7, 1.30],
  'marine-free|low|mean+2sd': [8.55e-14, 3.42, 1.93e-7, 1.30],
  'marine-free|high|mean': [5.37e-14, 3.42, 5.67e-7, 1.11],
  'marine-free|high|mean+2sd': [1.72e-13, 3.42, 7.48e-7, 1.11],
  'marine-cp850|low|mean': [1.21e-26, 8.16, 5.16e-12, 2.67],
  'marine-cp850|low|mean+2sd': [4.37e-26, 8.16, 1.32e-11, 2.67],
  'marine-cp850|high|mean': [4.80e-18, 5.10, 6.00e-12, 2.67],
  'marine-cp850|high|mean+2sd': [2.10e-17, 5.10, 2.02e-11, 2.67],
};

/**
 * BS 7910 recommended fatigue crack growth laws for steels (yield <= 600 MPa in marine environments).
 *  curve 'simplified': single stage, A = 5.21e-13, m = 3 in air; A = 2.3e-12, m = 3 in a marine environment
 *    (N–mm units) — BS 7910:2019 clause 8 (simplified laws).
 *  curve 'mean' | 'mean+2sd': two-stage law of the BS 7910:2019 clause 8 tables of recommended crack growth laws for R < 0.5 ('low') or R >= 0.5
 *    ('high', recommended for as-welded joints); environments 'air', 'marine-free' (free corrosion) and
 *    'marine-cp850' (cathodic protection at −850 mV Ag/AgCl).
 * Threshold ΔK0: 63 N/mm^(3/2) for welded joints in air and with cathodic protection; 0 for free corrosion
 * (BS 7910:2019 clause 8). Override with `threshold`.
 *
 * Threshold trap. The default threshold describes constant-amplitude loading: every cycle with ΔK <= ΔK0 gives
 * no growth at all. A variable-amplitude history is often condensed to an equivalent constant-amplitude range
 * (e.g. the m-th root mean stress range); that range commonly sits below ΔK0 for a small flaw — 63 N/mm^(3/2) =
 * 2.0 MPa·√m is reached by Δσ·Y·√(πa) = 36 MPa at a = 1 mm (Y = 1) — and the calculation then reports an
 * arrested crack although the large cycles of the real history do propagate it. For variable-amplitude loading
 * BS 7910 advises taking ΔK0 = 0 unless every cycle of the spectrum is below the threshold: pass
 * `variableAmplitude: true` (an explicit `threshold` still takes precedence).
 * @param {object} [o]
 * @param {string} [o.environment='air'] 'air' | 'marine' (simplified only) | 'marine-free' | 'marine-cp850'
 * @param {string} [o.curve='simplified'] 'simplified' | 'mean' | 'mean+2sd'
 * @param {string} [o.rRatio='high'] 'low' (R < 0.5) | 'high' (R >= 0.5)
 * @param {number} [o.threshold] ΔK0 [Pa·√m]
 * @param {boolean} [o.variableAmplitude=false] the law will be used with an equivalent stress range standing for a
 *   variable-amplitude history: ΔK0 = 0 (unless `threshold` is given)
 * @returns {{stages: {A:number, m:number}[], threshold:number, published: {A:number, m:number}[],
 *   transition:(number|null), label:string, variableAmplitude:boolean}} stages in SI (m/cycle, Pa·√m); threshold
 *   [Pa·√m]; `published` in N–mm units; transition ΔK between the stages [Pa·√m]
 */
export function parisLawBs7910({ environment = 'air', curve = 'simplified', rRatio = 'high', threshold, variableAmplitude = false } = {}) {
  let published;
  if (curve === 'simplified') {
    if (environment === 'air') published = [{ A: 5.21e-13, m: 3 }];
    else if (environment.startsWith('marine')) published = [{ A: 2.3e-12, m: 3 }];
    else throw new RangeError(`environment: unknown environment "${environment}"`);
  } else {
    const row = BS7910_TWO_STAGE[`${environment}|${rRatio}|${curve}`];
    need(row !== undefined, `no BS 7910 two-stage law for environment "${environment}", rRatio "${rRatio}", curve "${curve}"`);
    published = [{ A: row[0], m: row[1] }, { A: row[2], m: row[3] }];
  }
  const stages = published.map((s) => ({ A: parisConstantToSI(s.A, s.m), m: s.m }));
  const th = threshold !== undefined ? threshold
    : variableAmplitude || environment === 'marine-free' || environment === 'marine' ? 0 : kFromNmm32(63);
  need(th >= 0, 'threshold must be non-negative');
  let transition = null;
  if (published.length === 2) {
    transition = kFromNmm32((published[1].A / published[0].A) ** (1 / (published[0].m - published[1].m)));
  }
  return { stages, threshold: th, published, transition, label: `BS 7910 ${environment} ${curve}${curve === 'simplified' ? '' : ` R ${rRatio}`}`,
    variableAmplitude: !!variableAmplitude };
}

/**
 * Crack growth rate da/dN for a one- or two-stage Paris–Erdogan law with threshold and optional Walker
 * mean-stress correction ΔK_eff = ΔK/(1 − R)^(1 − γ).
 * Two stages are combined as the lower of the two lines (stage A governs below their intersection), which is the
 * BS 7910 bilinear relationship.
 * Reference: Paris & Erdogan (1963) J. Basic Eng. 85, 528–534; Walker (1970) ASTM STP 462; BS 7910:2019 clause 8.
 * @param {number} deltaK stress-intensity range [Pa·√m]
 * @param {{stages: {A:number, m:number}[], threshold?:number}} law SI constants
 * @param {object} [o]
 * @param {number} [o.R=0] stress ratio [-], < 1 (used only with walkerGamma)
 * @param {number} [o.walkerGamma] Walker exponent γ [-], typically 0.3–1 (1 = no effect)
 * @returns {number} da/dN [m/cycle]
 */
export function crackGrowthRate(deltaK, law, { R = 0, walkerGamma } = {}) {
  need(deltaK >= 0, 'deltaK must be non-negative');
  need(law && Array.isArray(law.stages) && law.stages.length >= 1, 'law.stages must be a non-empty array');
  let dk = deltaK;
  if (walkerGamma !== undefined) {
    need(R < 1, 'R must be less than 1');
    dk = deltaK / (1 - R) ** (1 - walkerGamma);
  }
  if (dk <= (law.threshold || 0)) return 0;
  let rate = Infinity;
  for (const s of law.stages) rate = Math.min(rate, s.A * dk ** s.m);
  return rate;
}

/**
 * Closed-form Paris life for a constant geometry factor Y (K = Y·Δσ·√(πa)), single stage, no threshold:
 *   N = [a_f^(1−m/2) − a_0^(1−m/2)] / [(1 − m/2)·A·(Y·Δσ·√π)^m]   (m ≠ 2),  N = ln(a_f/a_0)/(A·(Y·Δσ·√π)²) (m = 2).
 * @param {object} o
 * @param {number} o.a0 initial crack size [m]
 * @param {number} o.af final crack size [m]
 * @param {number} o.stressRange Δσ [Pa]
 * @param {number} o.Y geometry factor [-]
 * @param {number} o.A Paris coefficient (SI) [m/cycle/(Pa·√m)^m]
 * @param {number} o.m Paris exponent [-]
 * @returns {number} cycles [cycles]
 */
export function parisLifeClosedForm({ a0, af, stressRange, Y, A, m }) {
  positive(a0, 'a0');
  need(af >= a0, 'af must not be less than a0');
  positive(stressRange, 'stressRange');
  positive(Y, 'Y');
  positive(A, 'A');
  positive(m, 'm');
  const k = A * (Y * stressRange * Math.sqrt(Math.PI)) ** m;
  if (Math.abs(m - 2) < 1e-12) return Math.log(af / a0) / k;
  const e = 1 - m / 2;
  return (af ** e - a0 ** e) / (e * k);
}

/**
 * Crack size after N cycles for a constant geometry factor (inverse of {@link parisLifeClosedForm}).
 * @param {object} o
 * @param {number} o.a0 [m]
 * @param {number} o.cycles [cycles]
 * @param {number} o.stressRange [Pa]
 * @param {number} o.Y [-]
 * @param {number} o.A SI Paris coefficient
 * @param {number} o.m [-]
 * @returns {number} crack size [m]; Infinity once the closed form diverges (unstable growth)
 */
export function parisCrackSizeClosedForm({ a0, cycles, stressRange, Y, A, m }) {
  positive(a0, 'a0');
  need(cycles >= 0, 'cycles must be non-negative');
  const k = A * (Y * stressRange * Math.sqrt(Math.PI)) ** m;
  if (Math.abs(m - 2) < 1e-12) return a0 * Math.exp(k * cycles);
  const e = 1 - m / 2;
  const base = a0 ** e + e * k * cycles;
  if (base <= 0) return Infinity;
  return base ** (1 / e);
}

/**
 * Numerical Paris integration of a one-parameter crack with an arbitrary stress-intensity solution.
 * Mid-point (second-order) stepping with the increment limited to a fraction of the current crack size.
 * The law's threshold applies to the range given: with an equivalent constant-amplitude range of a
 * variable-amplitude history pass `variableAmplitude: true` (threshold taken as zero, BS 7910:2019 clause 8),
 * otherwise a range below ΔK0 returns `arrested` with no growth at all (see {@link parisLawBs7910}).
 * @param {object} o
 * @param {number} o.a0 initial crack size [m]
 * @param {number} o.af final crack size [m]
 * @param {function(number):number} o.deltaK stress-intensity range as a function of crack size [Pa·√m]
 * @param {object} o.law SI Paris law, see {@link parisLawBs7910}
 * @param {number} [o.maxCycles=Infinity] stop after this many cycles [cycles]
 * @param {number} [o.stepFraction=0.005] largest relative growth per step [-]
 * @param {object} [o.rateOptions] passed to {@link crackGrowthRate}
 * @param {boolean} [o.variableAmplitude=false] ignore the law's threshold (ΔK0 = 0)
 * @returns {{cycles:number, a:number, reached:boolean, arrested:boolean, history: {cycles:number[], a:number[]}}}
 */
export function parisIntegrate({ a0, af, deltaK, law: lawIn, maxCycles = Infinity, stepFraction = 0.005, rateOptions, variableAmplitude = false }) {
  const law = variableAmplitude ? { ...lawIn, threshold: 0 } : lawIn;
  positive(a0, 'a0');
  need(af > a0, 'af must exceed a0');
  need(stepFraction > 0 && stepFraction <= 0.1, 'stepFraction must be in (0, 0.1]');
  let a = a0;
  let n = 0;
  const hist = { cycles: [0], a: [a0] };
  let k = 0;
  while (a < af && n < maxCycles) {
    const r1 = crackGrowthRate(deltaK(a), law, rateOptions);
    if (!(r1 > 0)) return { cycles: Infinity, a, reached: false, arrested: true, history: hist };
    let da = Math.min(stepFraction * a, af - a);
    const rMid = crackGrowthRate(deltaK(a + 0.5 * da), law, rateOptions);
    let dn = da / rMid;
    if (n + dn > maxCycles) {
      dn = maxCycles - n;
      da = rMid * dn;
    }
    a += da;
    n += dn;
    if (++k % 8 === 0 || a >= af) { hist.cycles.push(n); hist.a.push(a); }
  }
  return { cycles: n, a, reached: a >= af * (1 - 1e-12), arrested: false, history: hist };
}

/**
 * Fatigue growth of a semi-elliptical surface crack in a plate: depth a and half-length c are advanced
 * simultaneously with the Newman–Raju stress intensities at the deepest and surface points.
 * Stops at wall break-through (a = throughFraction·t), at fracture (K_max >= toughness), when 2c reaches
 * widthFraction·W, or at maxCycles.
 * Reference: BS 7910:2019 clause 8 (general procedure); Newman & Raju (1981).
 * Threshold: membraneRange and bendingRange are treated as constant-amplitude ranges, so with the default
 * threshold of {@link parisLawBs7910} (2.0 MPa·√m in air or with cathodic protection) a small flaw under a
 * modest range does not grow at all and the result is 'arrested' with cycles = Infinity. When the ranges are an
 * equivalent constant-amplitude representation of a variable-amplitude history that is wrong — the larger cycles
 * of the real history do grow the crack — and BS 7910:2019 clause 8 advises ΔK0 = 0: pass
 * `variableAmplitude: true`, which integrates with the law's threshold set to zero.
 * @param {object} o
 * @param {number} o.a0 initial depth [m]
 * @param {number} o.c0 initial half-length [m]
 * @param {number} o.t thickness [m]
 * @param {number} [o.W=Infinity] plate width [m]
 * @param {number} [o.membraneRange=0] membrane stress range [Pa]
 * @param {number} [o.bendingRange=0] bending stress range [Pa]
 * @param {object} o.law SI Paris law
 * @param {object} [o.Mk] weld-toe magnification factors, constants or function(a, c) -> factors as in surfaceCrackK
 * @param {number} [o.maxCycles=Infinity] [cycles]
 * @param {number} [o.throughFraction=1] a/t at which the crack is taken as through-thickness [-]
 * @param {number} [o.widthFraction=0.8] 2c/W limit [-]
 * @param {number} [o.toughness] fracture toughness K_mat [Pa·√m]; with `peakFactor` enables the fracture stop
 * @param {number} [o.peakFactor=1] K_max/ΔK ratio = 1/(1 − R) for the fracture check [-]
 * @param {number} [o.stepFraction=0.005] largest relative growth per step [-]
 * @param {object} [o.rateOptions] passed to {@link crackGrowthRate}
 * @param {boolean} [o.variableAmplitude=false] the stress ranges are equivalent ranges of a variable-amplitude
 *   history: ignore the law's threshold (ΔK0 = 0)
 * @returns {{cycles:number, a:number, c:number, status:string, history: {cycles:number[], a:number[], c:number[]}}}
 *   status: 'through-thickness' | 'fracture' | 'width' | 'max-cycles' | 'arrested'; cycles = Infinity if arrested
 */
export function surfaceCrackGrowth({ a0, c0, t, W = Infinity, membraneRange = 0, bendingRange = 0, law: lawIn, Mk, maxCycles = Infinity,
  throughFraction = 1, widthFraction = 0.8, toughness, peakFactor = 1, stepFraction = 0.005, rateOptions, variableAmplitude = false }) {
  const law = variableAmplitude ? { ...lawIn, threshold: 0 } : lawIn;
  positive(a0, 'a0');
  positive(c0, 'c0');
  positive(t, 't');
  need(a0 < t, 'a0 must be less than t');
  need(throughFraction > 0 && throughFraction <= 1, 'throughFraction must be in (0, 1]');
  need(stepFraction > 0 && stepFraction <= 0.1, 'stepFraction must be in (0, 0.1]');
  const aEnd = throughFraction * t;
  const cEnd = Number.isFinite(W) ? 0.5 * widthFraction * W : Infinity;
  const mkOf = (a, c) => (typeof Mk === 'function' ? Mk(a, c) : Mk || {});
  // Growth of an initially deep crack (a/c > 2) is outside the solution; the aspect ratio is held at the limit.
  const kOf = (a, c) => surfaceCrackK({ a: Math.min(a, t), c: Math.max(c, a / 2), t, W, membrane: membraneRange, bending: bendingRange, Mk: mkOf(a, c) });
  let a = a0;
  let c = c0;
  let n = 0;
  let step = 0;
  const hist = { cycles: [0], a: [a0], c: [c0] };
  const done = (status, cycles = n) => {
    hist.cycles.push(cycles); hist.a.push(a); hist.c.push(c);
    return { cycles, a, c, status, history: hist };
  };
  for (;;) {
    if (a >= aEnd * (1 - 1e-12)) return done('through-thickness');
    if (c >= cEnd) return done('width');
    if (n >= maxCycles) return done('max-cycles');
    const k1 = kOf(a, c);
    if (toughness !== undefined && peakFactor * Math.max(k1.Kdeep, k1.Ksurface) >= toughness) return done('fracture');
    const ra = crackGrowthRate(Math.max(k1.Kdeep, 0), law, rateOptions);
    const rc = crackGrowthRate(Math.max(k1.Ksurface, 0), law, rateOptions);
    if (!(ra > 0) && !(rc > 0)) return done('arrested', Infinity);
    let dn = Math.min(ra > 0 ? (stepFraction * a) / ra : Infinity, rc > 0 ? (stepFraction * c) / rc : Infinity);
    const k2 = kOf(Math.min(a + 0.5 * ra * dn, t), c + 0.5 * rc * dn);
    const ra2 = crackGrowthRate(Math.max(k2.Kdeep, 0), law, rateOptions);
    const rc2 = crackGrowthRate(Math.max(k2.Ksurface, 0), law, rateOptions);
    if (ra2 > 0 && a + ra2 * dn > aEnd) dn = (aEnd - a) / ra2;
    if (n + dn > maxCycles) dn = maxCycles - n;
    a += ra2 * dn;
    c += rc2 * dn;
    n += dn;
    if (++step % 8 === 0) { hist.cycles.push(n); hist.a.push(a); hist.c.push(c); }
    need(step < 2e6, 'stepFraction: crack growth integration did not finish in 2e6 steps');
  }
}

// ---------------------------------------------------------------------------------------------------------------
// Failure assessment diagram
// ---------------------------------------------------------------------------------------------------------------

/**
 * BS 7910:2019 Option 1 failure assessment curve for materials without a yield discontinuity:
 *   f(Lr) = (1 + 0.5·Lr²)^(−1/2)·[0.3 + 0.7·exp(−μ·Lr⁶)]            for Lr <= 1,  μ = min(0.001·E/σy, 0.6)
 *   f(Lr) = f(1)·Lr^((N − 1)/(2N))                                    for 1 < Lr <= Lr,max,  N = 0.3·(1 − σy/σu)
 *   f(Lr) = 0                                                         for Lr > Lr,max = (σy + σu)/(2σy)
 * Reference: BS 7910:2019 clause 7 (Option 1 assessment line); identical to R6 Rev. 4 Option 1 and FITNET.
 * Validity: continuous-yielding steels; materials with a Lüders plateau need the yield-plateau variant.
 * @param {number} Lr load ratio σref/σy [-]
 * @param {object} o
 * @param {number} o.yieldStrength σy (0.2 % proof) [Pa]
 * @param {number} o.tensileStrength σu [Pa]
 * @param {number} [o.E=207e9] Young's modulus [Pa]
 * @returns {{Kr:number, LrMax:number, mu:number, N:number}} allowable fracture ratio [-] and curve parameters
 */
export function fadOption1(Lr, { yieldStrength, tensileStrength, E = 207e9 }) {
  need(Lr >= 0, 'Lr must be non-negative');
  positive(yieldStrength, 'yieldStrength');
  positive(E, 'E');
  need(tensileStrength >= yieldStrength, 'tensileStrength must not be less than yieldStrength');
  const mu = Math.min((0.001 * E) / yieldStrength, 0.6);
  const N = 0.3 * (1 - yieldStrength / tensileStrength);
  const LrMax = (yieldStrength + tensileStrength) / (2 * yieldStrength);
  const f = (x) => (0.3 + 0.7 * Math.exp(-mu * x ** 6)) / Math.sqrt(1 + 0.5 * x * x);
  let Kr;
  if (Lr <= 1) Kr = f(Lr);
  else if (Lr <= LrMax && N > 0) Kr = f(1) * Lr ** ((N - 1) / (2 * N));
  else Kr = 0;
  return { Kr, LrMax, mu, N };
}

/**
 * Reference stress for a surface flaw in a plate under membrane and bending stress:
 *   σref = [Pb + √(Pb² + 9·Pm²·(1 − α″)²)] / [3·(1 − α″)²]
 *   α″ = (a/t)/(1 + t/c) for W >= 2(c + t), otherwise α″ = (2a/t)·(c/W).
 * Reference: BS 7910:2019 Annex P (flat plate, surface flaw, negligible bending restraint).
 * @param {object} o
 * @param {number} o.a depth [m]
 * @param {number} o.c half-length [m]
 * @param {number} o.t thickness [m]
 * @param {number} [o.W=Infinity] plate width [m]
 * @param {number} o.Pm primary membrane stress [Pa]
 * @param {number} [o.Pb=0] primary bending stress [Pa]
 * @returns {number} reference stress [Pa]
 */
export function referenceStressSurfaceCrack({ a, c, t, W = Infinity, Pm, Pb = 0 }) {
  positive(a, 'a');
  positive(c, 'c');
  positive(t, 't');
  need(a <= t, 'a must not exceed t');
  const alpha = W >= 2 * (c + t) ? a / t / (1 + t / c) : ((2 * a) / t) * (c / W);
  need(alpha < 1, 'crack occupies the full section');
  const k = (1 - alpha) ** 2;
  return (Pb + Math.sqrt(Pb * Pb + 9 * Pm * Pm * k)) / (3 * k);
}

/**
 * Reference stress for a through-thickness flaw of length 2a in a plate:
 * σref = [Pb + √(Pb² + 9·Pm²)] / [3·(1 − 2a/W)].
 * Reference: BS 7910:2019 Annex P (flat plate, through-thickness flaw).
 * @param {object} o
 * @param {number} o.a half-length [m]
 * @param {number} o.W plate width [m]
 * @param {number} o.Pm [Pa]
 * @param {number} [o.Pb=0] [Pa]
 * @returns {number} reference stress [Pa]
 */
export function referenceStressThroughCrack({ a, W, Pm, Pb = 0 }) {
  positive(a, 'a');
  positive(W, 'W');
  need(2 * a < W, '2a must be less than W');
  return (Pb + Math.sqrt(Pb * Pb + 9 * Pm * Pm)) / (3 * (1 - (2 * a) / W));
}

/**
 * Plasticity interaction term ρ for combined primary and secondary stresses (simplified procedure):
 *   χ = K_I^s·Lr/K_I^p,  ρ1 = 0.1·χ^0.714 − 0.007·χ² + 0.00003·χ⁵ (χ <= 4), ρ1 = 0.25 (χ > 4)
 *   ρ = ρ1 for Lr <= 0.8;  ρ = 4·ρ1·(1.05 − Lr) for 0.8 < Lr < 1.05;  ρ = 0 for Lr >= 1.05.
 * Kr = (K_I^p + K_I^s)/K_mat + ρ.
 * Reference: BS 7910:2005 Annex R (simplified ρ method, from R6); retained here in place of the tabulated V
 * factor of the 2013/2019 editions.
 * @param {object} o
 * @param {number} o.Kp stress intensity from primary stresses [Pa·√m]
 * @param {number} o.Ks stress intensity from secondary stresses [Pa·√m]
 * @param {number} o.Lr load ratio [-]
 * @returns {number} ρ [-]
 */
export function plasticityInteractionRho({ Kp, Ks, Lr }) {
  need(Lr >= 0, 'Lr must be non-negative');
  if (!(Ks > 0) || !(Kp > 0) || Lr >= 1.05) return 0;
  const chi = (Ks * Lr) / Kp;
  const rho1 = chi > 4 ? 0.25 : 0.1 * chi ** 0.714 - 0.007 * chi * chi + 0.00003 * chi ** 5;
  return Lr <= 0.8 ? rho1 : 4 * rho1 * (1.05 - Lr);
}

function bisectRoot(f, lo, hi, iters = 100) {
  let flo = f(lo);
  for (let i = 0; i < iters; i++) {
    const mid = 0.5 * (lo + hi);
    const fm = f(mid);
    if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else hi = mid;
    if (hi - lo < 1e-12 * Math.max(1, Math.abs(hi))) break;
  }
  return 0.5 * (lo + hi);
}

/**
 * Engineering critical assessment of a surface crack in a plate by the Option 1 FAD.
 * Kr = (K_I^p + K_I^s)/K_mat + ρ at the more onerous of the deepest and surface points; Lr = σref/σy.
 * Acceptable when the point lies inside the curve and Lr <= Lr,max. Partial safety factors are applied by the
 * caller to the inputs.
 * Reference: BS 7910:2019 clause 7 (Option 1), Annex M (Newman–Raju K solution), Annex P (reference stress);
 * ρ interaction per BS 7910:2005 Annex R.
 * @param {object} o
 * @param {number} o.a depth [m]
 * @param {number} o.c half-length [m]
 * @param {number} o.t thickness [m]
 * @param {number} [o.W=Infinity] plate width [m]
 * @param {number} o.Pm primary membrane stress [Pa]
 * @param {number} [o.Pb=0] primary bending stress [Pa]
 * @param {number} [o.Qm=0] secondary (residual/thermal) membrane stress [Pa]
 * @param {number} [o.Qb=0] secondary bending stress [Pa]
 * @param {number} o.yieldStrength [Pa]
 * @param {number} o.tensileStrength [Pa]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} o.toughness K_mat [Pa·√m]
 * @param {object} [o.Mk] weld-toe magnification factors (applied to primary and secondary alike)
 * @returns {{Kr:number, Lr:number, KrAllowable:number, LrMax:number, Kp:number, Ks:number, rho:number,
 *   referenceStress:number, acceptable:boolean, margin:number, reserveFactor:number, governingPoint:string}}
 *   margin = KrAllowable/Kr at the assessed Lr (>= 1 acceptable); reserveFactor = factor on the primary stresses
 *   that brings the point onto the assessment line with the secondary stresses unchanged
 */
export function assessSurfaceCrack(o) {
  const { a, c, t, W = Infinity, Pm, Pb = 0, Qm = 0, Qb = 0, yieldStrength, tensileStrength, E = 207e9, toughness, Mk } = o;
  positive(toughness, 'toughness');
  need(Pm >= 0, 'Pm must be non-negative');
  const mat = { yieldStrength, tensileStrength, E };
  const unitM = surfaceCrackK({ a, c, t, W, membrane: 1, Mk });
  const unitB = surfaceCrackK({ a, c, t, W, bending: 1, Mk });
  const refUnit = (f) => referenceStressSurfaceCrack({ a, c, t, W, Pm: f * Pm, Pb: f * Pb });
  const point = (f) => {
    const Lr = refUnit(f) / yieldStrength;
    let best = null;
    for (const loc of ['deep', 'surface']) {
      const km = loc === 'deep' ? unitM.Kdeep : unitM.Ksurface;
      const kb = loc === 'deep' ? unitB.Kdeep : unitB.Ksurface;
      const Kp = Math.max(0, f * (km * Pm + kb * Pb));
      const Ks = Math.max(0, km * Qm + kb * Qb);
      const rho = plasticityInteractionRho({ Kp, Ks, Lr });
      const Kr = (Kp + Ks) / toughness + rho;
      if (!best || Kr > best.Kr) best = { Kr, Lr, Kp, Ks, rho, governingPoint: loc };
    }
    return best;
  };
  const base = point(1);
  const fad = fadOption1(base.Lr, mat);
  const inside = (f) => {
    const p = point(f);
    const lim = fadOption1(p.Lr, mat);
    return p.Lr <= lim.LrMax && p.Kr <= lim.Kr;
  };
  let reserveFactor;
  if (!inside(0)) reserveFactor = 0;
  else if (Pm === 0 && Pb === 0) reserveFactor = Infinity;
  else {
    let hi = 1;
    while (inside(hi) && hi < 1e6) hi *= 2;
    reserveFactor = hi >= 1e6 ? Infinity : bisectRoot((f) => (inside(f) ? 1 : -1), 0, hi, 60);
  }
  const acceptable = base.Lr <= fad.LrMax && base.Kr <= fad.Kr;
  return {
    Kr: base.Kr, Lr: base.Lr, KrAllowable: fad.Kr, LrMax: fad.LrMax, Kp: base.Kp, Ks: base.Ks, rho: base.rho,
    referenceStress: base.Lr * yieldStrength, acceptable, margin: base.Kr > 0 ? fad.Kr / base.Kr : Infinity,
    reserveFactor, governingPoint: base.governingPoint,
  };
}

/**
 * Engineering critical assessment of a through-thickness crack (length 2a) in a plate, Option 1 FAD.
 * K from the Tada centre-crack solution with membrane stress Pm + Qm (bending ignored for K, kept in σref).
 * Reference: BS 7910:2019 clause 7, Annexes M and P.
 * @param {object} o
 * @param {number} o.a half-length [m]
 * @param {number} o.W plate width [m]
 * @param {number} o.Pm [Pa]
 * @param {number} [o.Pb=0] [Pa]
 * @param {number} [o.Qm=0] [Pa]
 * @param {number} o.yieldStrength [Pa]
 * @param {number} o.tensileStrength [Pa]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} o.toughness [Pa·√m]
 * @returns {{Kr:number, Lr:number, KrAllowable:number, LrMax:number, acceptable:boolean, margin:number}}
 */
export function assessThroughCrack({ a, W, Pm, Pb = 0, Qm = 0, yieldStrength, tensileStrength, E = 207e9, toughness }) {
  positive(toughness, 'toughness');
  const unit = centreCrackK({ a, W, stress: 1 }).K;
  const Lr = referenceStressThroughCrack({ a, W, Pm, Pb }) / yieldStrength;
  const Kp = unit * Pm;
  const Ks = unit * Qm;
  const rho = plasticityInteractionRho({ Kp, Ks, Lr });
  const Kr = (Kp + Ks) / toughness + rho;
  const fad = fadOption1(Lr, { yieldStrength, tensileStrength, E });
  return { Kr, Lr, KrAllowable: fad.Kr, LrMax: fad.LrMax, acceptable: Lr <= fad.LrMax && Kr <= fad.Kr, margin: Kr > 0 ? fad.Kr / Kr : Infinity };
}

/**
 * Critical (limiting) depth of a surface crack of fixed aspect ratio a/c: the largest depth for which the
 * assessment point stays inside the Option 1 FAD, by bisection.
 * @param {object} o all arguments of {@link assessSurfaceCrack} except a and c, plus:
 * @param {number} o.aspectRatio a/c held constant [-], 0 < a/c <= 2
 * @param {number} [o.aMin] smallest depth considered [m], default t/1000
 * @returns {{a:number, c:number, limitedBy:string}} critical depth and half-length [m]; limitedBy 'fad', 'wall'
 *   (acceptable up to full thickness) or 'none-acceptable' (a = 0)
 */
export function criticalCrackSize(o) {
  const { aspectRatio, t, W = Infinity } = o;
  need(aspectRatio > 0 && aspectRatio <= 2, 'aspectRatio must be in (0, 2]');
  positive(t, 't');
  const aMin = o.aMin === undefined ? t / 1000 : o.aMin;
  let aMax = 0.999 * t;
  if (Number.isFinite(W)) aMax = Math.min(aMax, 0.49 * W * aspectRatio * 0.999);
  const ok = (a) => {
    try { return assessSurfaceCrack({ ...o, a, c: a / aspectRatio }).acceptable; } catch (e) { return false; }
  };
  if (!ok(aMin)) return { a: 0, c: 0, limitedBy: 'none-acceptable' };
  if (ok(aMax)) return { a: aMax, c: aMax / aspectRatio, limitedBy: 'wall' };
  const a = bisectRoot((x) => (ok(x) ? 1 : -1), aMin, aMax, 80);
  return { a, c: a / aspectRatio, limitedBy: 'fad' };
}

/**
 * Limiting primary stress for a surface crack: the factor on (Pm, Pb) that brings the assessment point onto the
 * FAD, returned as stresses.
 * @param {object} o arguments of {@link assessSurfaceCrack}
 * @returns {{factor:number, Pm:number, Pb:number}} limiting stresses [Pa]
 */
export function limitingStress(o) {
  const r = assessSurfaceCrack(o);
  return { factor: r.reserveFactor, Pm: r.reserveFactor * o.Pm, Pb: r.reserveFactor * (o.Pb || 0) };
}

/**
 * Fracture toughness from CTOD: K_mat = √(m·σy·δ·E/(1 − ν²)), m = 1.517·(σy/σu)^(−0.3188) for
 * 0.3 < σy/σu < 0.98, otherwise m = 1.5.
 * Reference: BS 7910:2019 clause 7 (fracture toughness from CTOD).
 * @param {object} o
 * @param {number} o.ctod δ [m]
 * @param {number} o.yieldStrength [Pa]
 * @param {number} [o.tensileStrength] [Pa]; needed for the strength-ratio form of m
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] Poisson's ratio [-]
 * @param {number} [o.m] override of the constraint factor [-], typically 1–2
 * @returns {{K:number, m:number}} K_mat [Pa·√m]
 */
export function ctodToK({ ctod, yieldStrength, tensileStrength, E = 207e9, nu = 0.3, m }) {
  positive(ctod, 'ctod');
  positive(yieldStrength, 'yieldStrength');
  let mm = m;
  if (mm === undefined) {
    const r = tensileStrength ? yieldStrength / tensileStrength : NaN;
    mm = r > 0.3 && r < 0.98 ? 1.517 * r ** -0.3188 : 1.5;
  }
  positive(mm, 'm');
  return { K: Math.sqrt((mm * yieldStrength * ctod * E) / (1 - nu * nu)), m: mm };
}

/**
 * CTOD equivalent to a stress-intensity value (inverse of {@link ctodToK}).
 * @param {object} o
 * @param {number} o.K [Pa·√m]
 * @param {number} o.yieldStrength [Pa]
 * @param {number} [o.tensileStrength] [Pa]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.m] [-]
 * @returns {number} δ [m]
 */
export function kToCtod({ K, yieldStrength, tensileStrength, E = 207e9, nu = 0.3, m }) {
  positive(K, 'K');
  const mm = ctodToK({ ctod: 1e-4, yieldStrength, tensileStrength, E, nu, m }).m;
  return (K * K * (1 - nu * nu)) / (mm * yieldStrength * E);
}

/**
 * Lower-bound fracture toughness from Charpy V-notch energy.
 *  'lower-bound' (lower shelf and transition): K_mat = [(12·√Cv − 20)·(25/B)^0.25] + 20, K in MPa·√m, Cv in J,
 *    B in mm — BS 7910:2019 Annex J.
 *  'upper-shelf': K_mat = 0.54·Cv + 55 (MPa·√m) — BS 7910:2019 Annex J (fully ductile behaviour only).
 * Validity: ferritic steels, Cv from standard 10 mm × 10 mm specimens at the assessment temperature.
 * @param {object} o
 * @param {number} o.charpyEnergy Cv [J]
 * @param {number} [o.thickness=0.025] section thickness B [m]
 * @param {string} [o.regime='lower-bound'] 'lower-bound' | 'upper-shelf'
 * @returns {number} K_mat [Pa·√m]
 */
export function charpyToK({ charpyEnergy, thickness = 0.025, regime = 'lower-bound' }) {
  positive(charpyEnergy, 'charpyEnergy');
  positive(thickness, 'thickness');
  if (regime === 'upper-shelf') return (0.54 * charpyEnergy + 55) * 1e6;
  need(regime === 'lower-bound', `regime: unknown regime "${regime}"`);
  const B = thickness * 1000;
  return ((12 * Math.sqrt(charpyEnergy) - 20) * (25 / B) ** 0.25 + 20) * 1e6;
}

/**
 * Leak-before-break check for a surface crack in a plate or shell wall treated as a plate: the crack is assumed
 * to break through the wall and is re-characterised as a through-thickness crack of length 2c + t; the case is
 * leak-before-break when that through crack is acceptable on the FAD.
 * Reference: BS 7910:2019 Annex F (leak-before-break procedure, detectable-leakage route not included) and
 * Annex E (re-characterisation, length 2c + t). `recharacterisedLength` overrides the default length.
 * @param {object} o
 * @param {number} o.c half-length of the surface crack at break-through [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.W plate width (or developed shell width) [m]
 * @param {number} o.Pm [Pa]
 * @param {number} [o.Pb=0] [Pa]
 * @param {number} [o.Qm=0] [Pa]
 * @param {number} o.yieldStrength [Pa]
 * @param {number} o.tensileStrength [Pa]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} o.toughness [Pa·√m]
 * @param {number} [o.recharacterisedLength] total through-crack length [m], default 2c + t
 * @returns {{leakBeforeBreak:boolean, throughLength:number, criticalThroughLength:number, margin:number,
 *   assessment:object}} lengths [m]; margin = critical / re-characterised length
 */
export function leakBeforeBreak(o) {
  const { c, t, W } = o;
  positive(c, 'c');
  positive(t, 't');
  positive(W, 'W');
  const len = o.recharacterisedLength === undefined ? 2 * c + t : o.recharacterisedLength;
  need(len < W, 'the re-characterised crack is longer than the plate width');
  const assess = (half) => assessThroughCrack({ ...o, a: half });
  const assessment = assess(len / 2);
  const hiMax = 0.4999 * W;
  let crit;
  const lo = Math.min(len / 2, W * 1e-6);
  if (!assess(lo).acceptable) crit = 0;
  else if (assess(hiMax).acceptable) crit = 2 * hiMax;
  else crit = 2 * bisectRoot((x) => (assess(x).acceptable ? 1 : -1), lo, hiMax, 80);
  return { leakBeforeBreak: assessment.acceptable, throughLength: len, criticalThroughLength: crit, margin: crit / len, assessment };
}
