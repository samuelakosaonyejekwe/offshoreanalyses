/**
 * Special functions: error function family, gamma and beta functions, Bessel functions of integer order.
 *
 * All arguments and results are dimensionless. Every function is computed from convergent series, continued
 * fractions or exponentially convergent quadrature evaluated to working precision (no fitted approximations),
 * so the accuracy is close to machine precision (relative ~1e-14 or better unless stated).
 *
 * References: Abramowitz & Stegun, "Handbook of Mathematical Functions" (A&S), chs. 6, 7, 9, 26;
 * NIST Digital Library of Mathematical Functions (DLMF) §5, §7, §8, §10;
 * Matta & Reichel, Math. Comp. 25 (1971) 339–344 (trapezoidal evaluation of erfc);
 * Lanczos, SIAM J. Numer. Anal. B1 (1964) 86–96 (gamma function, g = 7, 9 terms);
 * Press et al., "Numerical Recipes" §6.2, §6.4 (incomplete gamma/beta by modified Lentz continued fractions).
 */

const EPS = 2.220446049250313e-16;
const FPMIN = 1e-300;
const SQRT_PI = 1.7724538509055160273;
const TWO_OVER_PI = 0.63661977236758134308;
const EULER = 0.57721566490153286061;

function checkNum(x, name) {
  if (typeof x !== 'number' || Number.isNaN(x)) throw new RangeError(`${name}: argument must be a number`);
}

// ---------------------------------------------------------------------------------------------------------------
// Error function family
// ---------------------------------------------------------------------------------------------------------------

/** erf for small |x| from the all-positive series erf x = (2/√π) e^(−x²) Σ 2ⁿ x^(2n+1) / (2n+1)!! (A&S 7.1.6). */
function erfSeries(x) {
  const x2 = 2 * x * x;
  let term = x;
  let sum = x;
  for (let n = 1; n < 200; n++) {
    term *= x2 / (2 * n + 1);
    sum += term;
    if (term <= EPS * sum && term >= -EPS * Math.abs(sum)) break;
  }
  return (2 / SQRT_PI) * Math.exp(-x * x) * sum;
}

const ERFC_H = 0.5;
const ERFC_TERMS = 14;
const ERFC_E = new Float64Array(ERFC_TERMS + 1);
const ERFC_K2 = new Float64Array(ERFC_TERMS + 1);
for (let k = 1; k <= ERFC_TERMS; k++) {
  ERFC_K2[k] = k * k * ERFC_H * ERFC_H;
  ERFC_E[k] = Math.exp(-ERFC_K2[k]);
}

/**
 * erfc for x >= 0.5 from the trapezoidal rule applied to e^(x²) erfc x = (x/π) ∫ e^(−t²)/(t² + x²) dt with the pole
 * correction (Matta & Reichel 1971); truncation error ~ e^(−π²/h²) = 7e-18 for h = 1/2.
 */
function erfcTrap(x) {
  const x2 = x * x;
  let s = 0.5 / x2;
  for (let k = 1; k <= ERFC_TERMS; k++) s += ERFC_E[k] / (ERFC_K2[k] + x2);
  let r = ((2 * ERFC_H * x) / Math.PI) * Math.exp(-x2) * s;
  if (x < 6.5) r -= 2 / Math.expm1((2 * Math.PI * x) / ERFC_H);
  return r;
}

/**
 * Error function erf(x) = (2/√π) ∫₀ˣ e^(−t²) dt (A&S 7.1.1).
 * Validity: all real x; relative accuracy ~1e-15.
 * @param {number} x argument [-]
 * @returns {number} erf(x) in [−1, 1]
 */
export function erf(x) {
  checkNum(x, 'erf');
  const ax = Math.abs(x);
  if (ax < 0.5) return erfSeries(x);
  if (ax > 6) return x > 0 ? 1 : -1;
  const r = 1 - erfcTrap(ax);
  return x > 0 ? r : -r;
}

/**
 * Complementary error function erfc(x) = 1 − erf(x), accurate in the far tail (A&S 7.1.2).
 * Validity: all real x; relative accuracy ~1e-14 down to the underflow limit (x ≈ 26.5).
 * @param {number} x argument [-]
 * @returns {number} erfc(x) in [0, 2]
 */
export function erfc(x) {
  checkNum(x, 'erfc');
  if (x < 0.5) {
    if (x > -0.5) return 1 - erfSeries(x);
    if (x < -6) return 2;
    return 2 - erfcTrap(-x);
  }
  if (x > 27) return 0;
  return erfcTrap(x);
}

/**
 * Inverse error function: y such that erf(y) = p.
 * Start value from Winitzki's closed-form approximation, refined by Halley iterations on erf (on erfc in the
 * tails so that precision is kept as |p| → 1).
 * Validity: −1 <= p <= 1 (±1 give ±Infinity); throws RangeError outside.
 * @param {number} p probability-like argument in [−1, 1]
 * @returns {number} erf⁻¹(p)
 */
export function erfinv(p) {
  checkNum(p, 'erfinv');
  if (p < -1 || p > 1) throw new RangeError('erfinv: p must be in [-1, 1]');
  if (p === 0) return 0;
  if (p === 1) return Infinity;
  if (p === -1) return -Infinity;
  const ap = Math.abs(p);
  const q = 1 - ap; // exact for ap >= 0.5
  const a = 0.147;
  const ln = Math.log(q * (1 + ap));
  const t = 2 / (Math.PI * a) + 0.5 * ln;
  let y = Math.sqrt(Math.sqrt(t * t - ln / a) - t);
  if (!(y > 0)) y = (ap * SQRT_PI) / 2;
  for (let it = 0; it < 12; it++) {
    const f = ap > 0.7 ? q - erfc(y) : erf(y) - ap;
    const fp = (2 / SQRT_PI) * Math.exp(-y * y);
    if (fp === 0) break;
    const u = f / fp;
    const dy = u / (1 + y * u); // Halley: f''/f' = −2y
    y -= dy;
    if (Math.abs(dy) <= 2 * EPS * Math.abs(y)) break;
  }
  return p > 0 ? y : -y;
}

// ---------------------------------------------------------------------------------------------------------------
// Gamma and beta functions
// ---------------------------------------------------------------------------------------------------------------

const LANCZOS_G = 7;
const LANCZOS = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
];
const SQRT_2PI = 2.5066282746310005024;

const FACT = new Float64Array(171);
FACT[0] = 1;
for (let i = 1; i < 171; i++) FACT[i] = FACT[i - 1] * i;

function lanczosSum(xm1) {
  let a = LANCZOS[0];
  for (let i = 1; i < 9; i++) a += LANCZOS[i] / (xm1 + i);
  return a;
}

/**
 * Gamma function Γ(x) (A&S 6.1.1). Exact factorials for integer x <= 171, Lanczos approximation (g = 7, n = 9,
 * relative error < 1e-14) otherwise, reflection formula Γ(x)Γ(1−x) = π / sin(πx) for x < 0.5.
 * Validity: real x not a non-positive integer (poles: throws RangeError); overflows to Infinity for x > 171.62.
 * @param {number} x argument [-]
 * @returns {number} Γ(x)
 */
export function gamma(x) {
  checkNum(x, 'gamma');
  if (x === Infinity) return Infinity;
  if (Number.isInteger(x)) {
    if (x <= 0) throw new RangeError('gamma: x must not be zero or a negative integer (pole)');
    return x <= 171 ? FACT[x - 1] : Infinity;
  }
  if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
  const xm1 = x - 1;
  const t = xm1 + LANCZOS_G + 0.5;
  // split the power so that t^(x−½) does not overflow before the exponential brings it back
  const half = Math.pow(t, 0.5 * (xm1 + 0.5));
  return SQRT_2PI * half * Math.exp(-t) * half * lanczosSum(xm1);
}

/**
 * Natural logarithm of |Γ(x)| (A&S 6.1.40 region; Lanczos g = 7 with reflection for x < 0.5).
 * Validity: real x not a non-positive integer (throws RangeError at the poles); absolute accuracy ~1e-14·max(1, |ln Γ|).
 * @param {number} x argument [-]
 * @returns {number} ln|Γ(x)|
 */
export function lgamma(x) {
  checkNum(x, 'lgamma');
  if (x === Infinity) return Infinity;
  if (x <= 0 && Number.isInteger(x)) throw new RangeError('lgamma: x must not be zero or a negative integer (pole)');
  if (x === 1 || x === 2) return 0;
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
  const xm1 = x - 1;
  const t = xm1 + LANCZOS_G + 0.5;
  return (xm1 + 0.5) * Math.log(t) - t + Math.log(SQRT_2PI * lanczosSum(xm1));
}

function gammaSeriesP(a, x) {
  let ap = a;
  let del = 1 / a;
  let sum = del;
  for (let n = 0; n < 100000; n++) {
    ap += 1;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * EPS) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
}

function gammaContFracQ(a, x) {
  let b = x + 1 - a;
  let c = 1 / FPMIN;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 100000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
}

/**
 * Regularised lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a) (A&S 6.5.1): series for x < a + 1
 * (A&S 6.5.29), continued fraction for the complement otherwise (A&S 6.5.31, modified Lentz).
 * Validity: a > 0, x >= 0 (throws RangeError otherwise).
 * @param {number} a shape parameter (> 0)
 * @param {number} x upper limit (>= 0)
 * @returns {number} P(a, x) in [0, 1]
 */
export function gammaP(a, x) {
  if (!(a > 0)) throw new RangeError('gammaP: a must be > 0');
  if (!(x >= 0)) throw new RangeError('gammaP: x must be >= 0');
  if (x === 0) return 0;
  if (x === Infinity) return 1;
  if (x < a + 1) return Math.min(1, gammaSeriesP(a, x));
  return 1 - gammaContFracQ(a, x);
}

/**
 * Regularised upper incomplete gamma function Q(a, x) = 1 − P(a, x) = Γ(a, x) / Γ(a) (A&S 6.5.3), accurate in
 * the upper tail.
 * Validity: a > 0, x >= 0 (throws RangeError otherwise).
 * @param {number} a shape parameter (> 0)
 * @param {number} x lower limit (>= 0)
 * @returns {number} Q(a, x) in [0, 1]
 */
export function gammaQ(a, x) {
  if (!(a > 0)) throw new RangeError('gammaQ: a must be > 0');
  if (!(x >= 0)) throw new RangeError('gammaQ: x must be >= 0');
  if (x === 0) return 1;
  if (x === Infinity) return 0;
  if (x < a + 1) return Math.max(0, 1 - gammaSeriesP(a, x));
  return gammaContFracQ(a, x);
}

function betaContFrac(a, b, x) {
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m < 100000; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

/**
 * Regularised incomplete beta function I_x(a, b) = B(x; a, b) / B(a, b) (A&S 26.5.1), by the continued fraction
 * A&S 26.5.8 (modified Lentz) with the symmetry I_x(a, b) = 1 − I_(1−x)(b, a) for fast convergence.
 * Validity: a > 0, b > 0, 0 <= x <= 1 (throws RangeError otherwise).
 * @param {number} a first shape parameter (> 0)
 * @param {number} b second shape parameter (> 0)
 * @param {number} x argument in [0, 1]
 * @returns {number} I_x(a, b) in [0, 1]
 */
export function betaI(a, b, x) {
  if (!(a > 0)) throw new RangeError('betaI: a must be > 0');
  if (!(b > 0)) throw new RangeError('betaI: b must be > 0');
  if (!(x >= 0 && x <= 1)) throw new RangeError('betaI: x must be in [0, 1]');
  if (x === 0) return 0;
  if (x === 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log1p(-x));
  if (x < (a + 1) / (a + b + 2)) return (bt * betaContFrac(a, b, x)) / a;
  return 1 - (bt * betaContFrac(b, a, 1 - x)) / b;
}

// ---------------------------------------------------------------------------------------------------------------
// Bessel functions of the first and second kind, orders 0 and 1 (and integer order n for J)
// ---------------------------------------------------------------------------------------------------------------

const J_SERIES_MAX = 2; // ascending series below this argument
const J_ASYMP_MIN = 25; // Hankel asymptotic expansion from this argument (smallest term ~ e^(−2x) < 1e-21)
const jy = { j0: 0, j1: 0, y0: 0, y1: 0 };

/** Ascending series A&S 9.1.10–9.1.13 for 0 < x < 2 (no significant cancellation in that range). */
function jySeries(x) {
  const q = 0.25 * x * x;
  let t0 = 1; // (−q)^k / (k!)²
  let t1 = 1; // (−q)^k / (k!(k+1)!)
  let j0 = 1;
  let j1 = 1;
  let hk = 0; // harmonic number H_k
  let s0 = 0;
  let s1 = 1 - 2 * EULER; // k = 0 term: ψ(1) + ψ(2)
  for (let k = 1; k < 60; k++) {
    t0 *= -q / (k * k);
    t1 *= -q / (k * (k + 1));
    hk += 1 / k;
    j0 += t0;
    j1 += t1;
    s0 -= hk * t0;
    s1 += (2 * hk + 1 / (k + 1) - 2 * EULER) * t1;
    if (Math.abs(t0) < 1e-18) break;
  }
  j1 *= 0.5 * x;
  const lg = Math.log(0.5 * x);
  jy.j0 = j0;
  jy.j1 = j1;
  jy.y0 = TWO_OVER_PI * ((lg + EULER) * j0 + s0);
  jy.y1 = -TWO_OVER_PI / x + TWO_OVER_PI * lg * j1 - (x / (2 * Math.PI)) * s1;
}

/**
 * Miller backward recurrence for 2 <= x < 25, normalised with 1 = J0 + 2 Σ J_2k (A&S 9.1.46); Y0 from the Neumann
 * series A&S 9.1.89 and Y1 from its derivative:
 *   Y1 = −(2/(πx)) J0 + (2/π)(ln(x/2) + γ) J1 + (2/π) Σ_{k>=1} (−1)^k (J_{2k−1} − J_{2k+1}) / k.
 */
function jyMiller(x) {
  const N = 2 * Math.ceil(0.7 * x + 16);
  const tox = 2 / x;
  let bjp = 0;
  let bj = 1;
  let norm = 0;
  let s0 = 0;
  let s1 = 0;
  let v0 = 0;
  let v1 = 0;
  for (let m = N; m >= 0; m--) {
    // bj holds the unnormalised J_m
    if (m === 0) {
      norm += bj;
      v0 = bj;
    } else if ((m & 1) === 0) {
      norm += 2 * bj;
      const k = m >> 1;
      s0 += ((k & 1) === 0 ? bj : -bj) / k;
    } else {
      const ka = (m + 1) >> 1;
      let c = ((ka & 1) === 0 ? 1 : -1) / ka;
      if (m >= 3) {
        const kb = (m - 1) >> 1;
        c -= ((kb & 1) === 0 ? 1 : -1) / kb;
      }
      s1 += c * bj;
      if (m === 1) v1 = bj;
    }
    const bjm = m * tox * bj - bjp;
    bjp = bj;
    bj = bjm;
  }
  const j0 = v0 / norm;
  const j1 = v1 / norm;
  const lg = Math.log(0.5 * x) + EULER;
  jy.j0 = j0;
  jy.j1 = j1;
  jy.y0 = TWO_OVER_PI * lg * j0 - 2 * TWO_OVER_PI * (s0 / norm);
  jy.y1 = (-TWO_OVER_PI / x) * j0 + TWO_OVER_PI * lg * j1 + TWO_OVER_PI * (s1 / norm);
}

/** Hankel asymptotic expansions A&S 9.2.5–9.2.10 for x >= 25, summed to the smallest term. */
function jyAsymp(x) {
  const sx = Math.sin(x);
  const cx = Math.cos(x);
  const amp = Math.sqrt(TWO_OVER_PI / x);
  const e8 = 8 * x;
  for (let n = 0; n <= 1; n++) {
    const mu = 4 * n * n;
    let P = 1;
    let Q = 0;
    let term = 1;
    let prev = Infinity;
    for (let k = 1; k < 200; k++) {
      const odd = 2 * k - 1;
      term *= (mu - odd * odd) / (k * e8);
      const at = Math.abs(term);
      if (at > prev) break;
      prev = at;
      // k = 1, 5, 9.. -> +Q; k = 2, 6.. -> −P; k = 3, 7.. -> −Q; k = 4, 8.. -> +P
      const r = k & 3;
      if (r === 1) Q += term;
      else if (r === 2) P -= term;
      else if (r === 3) Q -= term;
      else P += term;
      if (at < 1e-17) break;
    }
    // χ = x − (n/2 + 1/4)π, expanded to avoid subtracting a large multiple of π
    let cchi;
    let schi;
    if (n === 0) {
      cchi = (cx + sx) * Math.SQRT1_2;
      schi = (sx - cx) * Math.SQRT1_2;
    } else {
      cchi = (sx - cx) * Math.SQRT1_2;
      schi = -(sx + cx) * Math.SQRT1_2;
    }
    const j = amp * (P * cchi - Q * schi);
    const y = amp * (P * schi + Q * cchi);
    if (n === 0) {
      jy.j0 = j;
      jy.y0 = y;
    } else {
      jy.j1 = j;
      jy.y1 = y;
    }
  }
}

function jyEval(x) {
  if (x < J_SERIES_MAX) jySeries(x);
  else if (x < J_ASYMP_MIN) jyMiller(x);
  else jyAsymp(x);
}

/**
 * Bessel function of the first kind, order 0 (A&S 9.1).
 * Validity: all real x; absolute accuracy ~1e-15.
 * @param {number} x argument [-]
 * @returns {number} J0(x)
 */
export function besselJ0(x) {
  checkNum(x, 'besselJ0');
  const ax = Math.abs(x);
  if (ax === 0) return 1;
  if (ax === Infinity) return 0;
  jyEval(ax);
  return jy.j0;
}

/**
 * Bessel function of the first kind, order 1 (odd in x).
 * Validity: all real x; absolute accuracy ~1e-15.
 * @param {number} x argument [-]
 * @returns {number} J1(x)
 */
export function besselJ1(x) {
  checkNum(x, 'besselJ1');
  const ax = Math.abs(x);
  if (ax === 0) return 0;
  if (ax === Infinity) return 0;
  jyEval(ax);
  return x < 0 ? -jy.j1 : jy.j1;
}

/**
 * Bessel function of the second kind (Weber/Neumann), order 0.
 * Validity: x > 0 (throws RangeError otherwise).
 * @param {number} x argument (> 0)
 * @returns {number} Y0(x)
 */
export function besselY0(x) {
  if (!(x > 0)) throw new RangeError('besselY0: x must be > 0');
  if (x === Infinity) return 0;
  jyEval(x);
  return jy.y0;
}

/**
 * Bessel function of the second kind, order 1.
 * Validity: x > 0 (throws RangeError otherwise).
 * @param {number} x argument (> 0)
 * @returns {number} Y1(x)
 */
export function besselY1(x) {
  if (!(x > 0)) throw new RangeError('besselY1: x must be > 0');
  if (x === Infinity) return 0;
  jyEval(x);
  return jy.y1;
}

/**
 * Bessel function of the first kind of integer order n (A&S 9.1.27 recurrences): forward recurrence from J0, J1
 * when x > n (stable), Miller's backward recurrence with the normalisation 1 = J0 + 2 Σ J_2k otherwise.
 * Negative orders use J_(−n) = (−1)ⁿ J_n; negative arguments J_n(−x) = (−1)ⁿ J_n(x).
 * Validity: integer n (throws RangeError otherwise), all real x.
 * @param {number} n integer order
 * @param {number} x argument [-]
 * @returns {number} J_n(x)
 */
export function besselJn(n, x) {
  if (!Number.isInteger(n)) throw new RangeError('besselJn: n must be an integer');
  checkNum(x, 'besselJn');
  let sign = 1;
  if (n < 0) {
    n = -n;
    if (n & 1) sign = -sign;
  }
  if (x < 0) {
    x = -x;
    if (n & 1) sign = -sign;
  }
  if (n === 0) return sign * besselJ0(x);
  if (n === 1) return sign * besselJ1(x);
  if (x === 0 || x === Infinity) return 0;
  const tox = 2 / x;
  if (x > n) {
    jyEval(x);
    let bjm = jy.j0;
    let bj = jy.j1;
    for (let j = 1; j < n; j++) {
      const bjp = j * tox * bj - bjm;
      bjm = bj;
      bj = bjp;
    }
    return sign * bj;
  }
  const m = 2 * Math.ceil(0.5 * (n + Math.sqrt(160 * n) + 10));
  let bjp = 0;
  let bj = 1;
  let sum = 0;
  let ans = 0;
  for (let j = m; j >= 1; j--) {
    // bj = unnormalised J_j
    if (j === n) ans = bj;
    if ((j & 1) === 0) sum += 2 * bj;
    const bjm = j * tox * bj - bjp;
    bjp = bj;
    bj = bjm;
    if (Math.abs(bj) > 1e250) {
      bj *= 1e-250;
      bjp *= 1e-250;
      ans *= 1e-250;
      sum *= 1e-250;
    }
  }
  sum += bj; // J0 term
  return (sign * ans) / sum;
}

// ---------------------------------------------------------------------------------------------------------------
// Modified Bessel functions, orders 0 and 1
// ---------------------------------------------------------------------------------------------------------------

/** Ascending series A&S 9.6.10 (all terms positive): returns Σ q^k / (k!(k+n)!), q = x²/4. */
function iSeries(n, q) {
  let term = 1;
  let sum = 1;
  for (let k = 1; k < 500; k++) {
    term *= q / (k * (k + n));
    sum += term;
    if (term < EPS * sum) break;
  }
  return sum;
}

/** Asymptotic expansion A&S 9.7.1 for large x, without the e^x / sqrt(2πx) prefactor. */
function iAsympSum(n, x) {
  const mu = 4 * n * n;
  let term = 1;
  let sum = 1;
  let prev = Infinity;
  for (let k = 1; k < 200; k++) {
    const odd = 2 * k - 1;
    term *= -(mu - odd * odd) / (k * 8 * x);
    const at = Math.abs(term);
    if (at > prev) break;
    prev = at;
    sum += term;
    if (at < 1e-17 * Math.abs(sum)) break;
  }
  return sum;
}

/**
 * Modified Bessel function of the first kind, order 0 (A&S 9.6.12 series; A&S 9.7.1 asymptotic for |x| > 30).
 * Validity: all real x (overflows to Infinity for |x| > ~713).
 * @param {number} x argument [-]
 * @returns {number} I0(x) >= 1
 */
export function besselI0(x) {
  checkNum(x, 'besselI0');
  const ax = Math.abs(x);
  if (ax <= 30) return iSeries(0, 0.25 * ax * ax);
  // e^x / sqrt(2πx) split so that the result only overflows when I0 itself does
  const h = Math.exp(0.5 * ax);
  return (h / Math.sqrt(2 * Math.PI * ax)) * iAsympSum(0, ax) * h;
}

/**
 * Modified Bessel function of the first kind, order 1 (odd in x).
 * Validity: all real x (overflows to ±Infinity for |x| > ~713).
 * @param {number} x argument [-]
 * @returns {number} I1(x)
 */
export function besselI1(x) {
  checkNum(x, 'besselI1');
  const ax = Math.abs(x);
  let r;
  if (ax <= 30) {
    r = 0.5 * ax * iSeries(1, 0.25 * ax * ax);
  } else {
    const h = Math.exp(0.5 * ax);
    r = (h / Math.sqrt(2 * Math.PI * ax)) * iAsympSum(1, ax) * h;
  }
  return x < 0 ? -r : r;
}

/**
 * K_n(x) for x > 1 from the integral representation K_ν(x) = ∫₀^∞ exp(−x cosh t) cosh(νt) dt (A&S 9.6.24),
 * evaluated by the trapezoidal rule, which converges exponentially for this analytic, doubly exponentially
 * decaying integrand. The step is chosen so that the discretisation error is below ~1e-17 relative.
 */
function kIntegral(n, x) {
  const h = Math.min((Math.PI * Math.PI) / (40 + x), 0.68 / Math.sqrt(x));
  let sum = 0.5;
  for (let k = 1; k < 2000; k++) {
    const t = k * h;
    const sh = Math.sinh(0.5 * t);
    const arg = 2 * x * sh * sh; // x (cosh t − 1)
    if (arg > 45) break;
    const e = Math.exp(-arg);
    sum += n === 0 ? e : e * Math.cosh(t);
  }
  return h * sum * Math.exp(-x);
}

/**
 * Modified Bessel function of the second kind, order 0: ascending series A&S 9.6.13 for x <= 1, trapezoidal
 * evaluation of the integral representation A&S 9.6.24 otherwise.
 * Validity: x > 0 (throws RangeError otherwise); underflows to 0 for x > ~745.
 * @param {number} x argument (> 0)
 * @returns {number} K0(x) > 0
 */
export function besselK0(x) {
  if (!(x > 0)) throw new RangeError('besselK0: x must be > 0');
  if (x === Infinity) return 0;
  if (x > 1) return kIntegral(0, x);
  const q = 0.25 * x * x;
  let term = 1;
  let i0 = 1;
  let hk = 0;
  let s = 0;
  for (let k = 1; k < 60; k++) {
    term *= q / (k * k);
    hk += 1 / k;
    i0 += term;
    s += hk * term;
    if (term < 1e-18) break;
  }
  return -(Math.log(0.5 * x) + EULER) * i0 + s;
}

/**
 * Modified Bessel function of the second kind, order 1: ascending series A&S 9.6.11 for x <= 1, trapezoidal
 * evaluation of the integral representation A&S 9.6.24 otherwise.
 * Validity: x > 0 (throws RangeError otherwise).
 * @param {number} x argument (> 0)
 * @returns {number} K1(x) > 0
 */
export function besselK1(x) {
  if (!(x > 0)) throw new RangeError('besselK1: x must be > 0');
  if (x === Infinity) return 0;
  if (x > 1) return kIntegral(1, x);
  const q = 0.25 * x * x;
  let term = 1; // q^k / (k!(k+1)!)
  let i1 = 1;
  let hk = 0;
  let s = 1 - 2 * EULER; // ψ(1) + ψ(2)
  for (let k = 1; k < 60; k++) {
    term *= q / (k * (k + 1));
    hk += 1 / k;
    i1 += term;
    s += (2 * hk + 1 / (k + 1) - 2 * EULER) * term;
    if (term < 1e-18) break;
  }
  return 1 / x + Math.log(0.5 * x) * (0.5 * x * i1) - 0.25 * x * s;
}
