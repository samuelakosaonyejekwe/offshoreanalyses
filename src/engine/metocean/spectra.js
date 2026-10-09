// Wave spectra, directional spreading, spectral moments, short-term wave and crest statistics, near-bed
// wave-induced velocity and random sea-state synthesis.
//
// Conventions: spectra are one-sided variance densities S(ω) in m²·s/rad as functions of angular frequency
// ω (rad/s), so that m0 = ∫ S dω = Hs²/16. Convert to S(f) in m²/Hz with S(f) = 2π S(ω), f = ω/2π.

import { gamma as gammaFn } from '../math/special.js';
import { makeRng } from '../math/rng.js';

const G = 9.80665;
const TWO_PI = 2 * Math.PI;
const EULER = 0.5772156649015329;

function reqPos(name, v) {
  if (typeof v !== 'number' || !(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name} must be a positive finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Frequency spectra                                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Pierson–Moskowitz spectrum in the Hs–Tp form, S(ω) = (5/16) Hs² ωp⁴ ω⁻⁵ exp[−(5/4)(ω/ωp)⁻⁴].
 * Reference: DNV-RP-C205 (2019) §3.5.5.1; Pierson & Moskowitz (1964).
 * Validity: fully developed wind sea; ω > 0 (returns 0 for ω ≤ 0).
 * @param {number} omega angular frequency, rad/s
 * @param {number} Hs significant wave height, m
 * @param {number} Tp peak period, s
 * @returns {number} S(ω), m²·s/rad
 */
export function piersonMoskowitz(omega, Hs, Tp) {
  if (!(Hs >= 0)) throw new RangeError('Hs must be non-negative');
  reqPos('Tp', Tp);
  if (!(omega > 0)) return 0;
  const wp = TWO_PI / Tp;
  const r = wp / omega;
  const r4 = r * r * r * r;
  return ((5 / 16) * Hs * Hs * r4 * Math.exp(-1.25 * r4)) / omega;
}

/**
 * Peak-enhancement factor γ recommended when it is not given:
 * γ = 5 for Tp/√Hs ≤ 3.6;  γ = exp(5.75 − 1.15 Tp/√Hs) for 3.6 < Tp/√Hs < 5;  γ = 1 for Tp/√Hs ≥ 5
 * (Tp in s, Hs in m).
 * Reference: DNV-RP-C205 (2019) §3.5.5.5.
 * @param {number} Hs significant wave height, m
 * @param {number} Tp peak period, s
 * @returns {number} γ, -
 */
export function jonswapGamma(Hs, Tp) {
  reqPos('Hs', Hs); reqPos('Tp', Tp);
  const r = Tp / Math.sqrt(Hs);
  if (r <= 3.6) return 5;
  if (r >= 5) return 1;
  return Math.exp(5.75 - 1.15 * r);
}

/**
 * JONSWAP spectrum, S_J(ω) = Aγ S_PM(ω) γ^exp[−½((ω − ωp)/(σ ωp))²], Aγ = 1 − 0.287 ln γ,
 * σ = σa for ω ≤ ωp and σb for ω > ωp.
 * Reference: DNV-RP-C205 (2019) §3.5.5.2–3.5.5.3; Hasselmann et al. (1973).
 * Validity: 3.6 < Tp/√Hs < 5 for a reasonable model; the normalising factor Aγ is an approximation: m0 is
 * within 0.3 % of Hs²/16 for 1 ≤ γ ≤ 5 and about 2 % low at γ = 7.
 * @param {number} omega angular frequency, rad/s
 * @param {object} o
 * @param {number} o.Hs significant wave height, m
 * @param {number} o.Tp peak period, s
 * @param {number} [o.gamma] peak-enhancement factor (default from `jonswapGamma`), -
 * @param {number} [o.sigmaA=0.07] spectral width below the peak, -
 * @param {number} [o.sigmaB=0.09] spectral width above the peak, -
 * @returns {number} S(ω), m²·s/rad
 */
export function jonswap(omega, { Hs, Tp, gamma, sigmaA = 0.07, sigmaB = 0.09 }) {
  const gam = gamma ?? jonswapGamma(Hs, Tp);
  if (!(gam >= 1)) throw new RangeError('gamma must be >= 1');
  if (!(omega > 0)) return 0;
  const wp = TWO_PI / Tp;
  const sig = omega <= wp ? sigmaA : sigmaB;
  const x = (omega - wp) / (sig * wp);
  return (1 - 0.287 * Math.log(gam)) * piersonMoskowitz(omega, Hs, Tp) * Math.pow(gam, Math.exp(-0.5 * x * x));
}

/**
 * Period ratios of the JONSWAP spectrum as functions of γ (curve fits):
 *   Tm01/Tp = 0.7303 + 0.04936 γ − 0.006556 γ² + 0.0003610 γ³
 *   Tm02/Tp = 0.6673 + 0.05037 γ − 0.006230 γ² + 0.0003341 γ³      (Tm02 = Tz)
 * Reference: DNV-RP-C205 (2019) §3.5.5.4.
 * Validity: 1 ≤ γ ≤ 7.
 * @param {number} gamma peak-enhancement factor, -
 * @returns {{Tm01OverTp:number,TzOverTp:number}} dimensionless
 */
export function jonswapPeriodRatios(gamma) {
  if (!(gamma >= 1) || gamma > 7) throw new RangeError('gamma must lie in [1, 7]');
  const g2 = gamma * gamma; const g3 = g2 * gamma;
  return {
    Tm01OverTp: 0.7303 + 0.04936 * gamma - 0.006556 * g2 + 0.0003610 * g3,
    TzOverTp: 0.6673 + 0.05037 * gamma - 0.006230 * g2 + 0.0003341 * g3,
  };
}

// Linear wavenumber for angular frequency w at depth d (Newton on k tanh kd = w²/g).
function kOf(w, d, g) {
  const k0 = (w * w) / g;
  if (!Number.isFinite(d)) return k0;
  let k = k0 / Math.sqrt(Math.tanh(k0 * d));
  for (let i = 0; i < 50; i++) {
    const kd = k * d;
    const th = Math.tanh(kd);
    const f = k * th - k0;
    const df = th + kd * (1 - th * th);
    const kn = k - f / df;
    if (Math.abs(kn - k) <= 1e-15 * kn) { k = kn; break; }
    k = kn;
  }
  return k;
}

/**
 * Finite-depth transformation factor of the TMA spectrum (Kitaigorodskii similarity),
 * φ(ω, d) = [k⁻³ ∂k/∂ω]_d / [k⁻³ ∂k/∂ω]_∞ = tanh²(kd) / (1 + 2kd / sinh 2kd), with k from linear dispersion.
 * Limits: φ → 1 in deep water and φ → ½ ω² d / g in shallow water.
 * Reference: Kitaigorodskii, Krasitskii & Zaslavskii (1975) J. Phys. Oceanogr. 5; Bouws, Günther, Rosenthal &
 * Vincent (1985) J. Geophys. Res. 90(C1), Eq. 5; used by DNV-RP-C205 (2019) §3.5.6.
 * Validity: wind seas in finite depth over a gently sloping bed, non-breaking.
 * @param {number} omega angular frequency, rad/s
 * @param {number} d water depth, m
 * @param {number} [g=9.80665] m/s²
 * @returns {number} φ, -
 */
export function tmaFactor(omega, d, g = G) {
  reqPos('d', d);
  if (!(omega > 0)) return 0;
  const kd = kOf(omega, d, g) * d;
  if (kd > 350) return 1;
  const th = Math.tanh(kd);
  return (th * th) / (1 + (2 * kd) / Math.sinh(2 * kd));
}

/**
 * TMA finite-depth spectrum, S_TMA(ω) = S_JONSWAP(ω) φ(ω, d).
 * Note that Hs of the result is lower than the JONSWAP input Hs (the JONSWAP parameters describe the
 * equivalent deep-water spectrum); integrate with `spectrumMoments` for the actual m0.
 * Reference: Bouws et al. (1985); DNV-RP-C205 (2019) §3.5.6.
 * Validity: as `tmaFactor`.
 * @param {number} omega angular frequency, rad/s
 * @param {object} o JONSWAP parameters as in `jonswap` plus
 * @param {number} o.d water depth, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {number} S(ω), m²·s/rad
 */
export function tma(omega, o) {
  return jonswap(omega, o) * tmaFactor(omega, o.d, o.g ?? G);
}

/**
 * Bretschneider / ITTC / ISSC two-parameter spectrum S(ω) = A ω⁻⁵ exp(−B ω⁻⁴) with A = 4 B m0 = B Hs²/4 and B
 * fixed by the stated period, derived from the closed-form moments m_n = (A/4) B^{n/4−1} Γ(1 − n/4):
 *   periodType 'Tp' : B = (5/4) ωp⁴                (identical to the Pierson–Moskowitz Hs–Tp form)
 *   periodType 'Tz' : B = (1/π)(2π/Tz)⁴            (ITTC 1978: A = 4π³Hs²/Tz⁴, B = 16π³/Tz⁴)
 *   periodType 'T1' : B = [2π/(T1 Γ(3/4))]⁴        (ISSC 1964 mean-period form, B ≈ 691/T1⁴)
 * Reference: Bretschneider (1959); ITTC (1978) 15th Proc.; DNV-RP-C205 (2019) §3.5.5.1.
 * Validity: open-ocean single-peaked sea states.
 * @param {number} omega angular frequency, rad/s
 * @param {object} o
 * @param {number} o.Hs significant wave height, m
 * @param {number} o.T period of the stated type, s
 * @param {'Tp'|'Tz'|'T1'} [o.periodType='Tp']
 * @returns {number} S(ω), m²·s/rad
 */
export function bretschneider(omega, { Hs, T, periodType = 'Tp' }) {
  if (!(Hs >= 0)) throw new RangeError('Hs must be non-negative');
  reqPos('T', T);
  let B;
  const w = TWO_PI / T;
  if (periodType === 'Tp') B = 1.25 * w ** 4;
  else if (periodType === 'Tz') B = w ** 4 / Math.PI;
  else if (periodType === 'T1') B = (w / gammaFn(0.75)) ** 4;
  else throw new RangeError("periodType must be 'Tp', 'Tz' or 'T1'");
  if (!(omega > 0)) return 0;
  const w4 = omega ** 4;
  return ((0.25 * B * Hs * Hs) / (w4 * omega)) * Math.exp(-B / w4);
}

/**
 * Ochi–Hubble six-parameter double-peak spectrum (sum of two three-parameter components j = 1 swell, 2 wind sea):
 *   S(ω) = ¼ Σ_j [((4λ_j + 1)/4) ωp_j⁴]^λ_j / Γ(λ_j) · Hs_j² / ω^(4λ_j + 1) · exp[−((4λ_j + 1)/4)(ωp_j/ω)⁴].
 * Each component integrates to Hs_j²/16 and peaks at ωp_j; total Hs = sqrt(Hs1² + Hs2²). λ = 1 gives the
 * Pierson–Moskowitz shape.
 * Reference: Ochi, M.K. & Hubble, E.N. (1976) "Six-parameter wave spectra", Proc. 15th ICCE, 301–328;
 * DNV-RP-C205 (2019) §3.5.7.
 * Validity: λ_j > 0; combined swell and wind sea.
 * @param {number} omega angular frequency, rad/s
 * @param {object} o
 * @param {number} o.Hs1 significant height of the low-frequency component, m
 * @param {number} o.Tp1 peak period of the low-frequency component, s
 * @param {number} o.lambda1 shape parameter of the low-frequency component, -
 * @param {number} o.Hs2 significant height of the high-frequency component, m
 * @param {number} o.Tp2 peak period of the high-frequency component, s
 * @param {number} o.lambda2 shape parameter of the high-frequency component, -
 * @returns {number} S(ω), m²·s/rad
 */
export function ochiHubble(omega, { Hs1, Tp1, lambda1, Hs2, Tp2, lambda2 }) {
  reqPos('Tp1', Tp1); reqPos('Tp2', Tp2); reqPos('lambda1', lambda1); reqPos('lambda2', lambda2);
  if (!(Hs1 >= 0) || !(Hs2 >= 0)) throw new RangeError('Hs1 and Hs2 must be non-negative');
  if (!(omega > 0)) return 0;
  const part = (Hs, Tp, lam) => {
    if (Hs === 0) return 0;
    const wp = TWO_PI / Tp;
    const a = (4 * lam + 1) / 4;
    const r4 = (wp / omega) ** 4;
    // [a ωp⁴]^λ / ω^(4λ+1) = (a r⁴)^λ / ω, evaluated in logs to avoid overflow for large λ
    return (0.25 * Hs * Hs * Math.exp(lam * Math.log(a * r4) - a * r4)) / (gammaFn(lam) * omega);
  };
  return part(Hs1, Tp1, lambda1) + part(Hs2, Tp2, lambda2);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Directional spreading                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

function wrapPi(a) {
  let x = a % TWO_PI;
  if (x > Math.PI) x -= TWO_PI;
  if (x < -Math.PI) x += TWO_PI;
  return x;
}

/**
 * cosⁿ directional spreading function, D(θ) = Γ(1 + n/2) / [√π Γ(½ + n/2)] cosⁿ(θ − θp) for |θ − θp| ≤ π/2,
 * zero otherwise; ∫ D dθ = 1.
 * Reference: DNV-RP-C205 (2019) §3.5.8.4 (typically n = 2–4 for wind sea, n ≥ 6 for swell).
 * @param {number} theta direction, rad
 * @param {number} thetaP main direction, rad
 * @param {number} n spreading exponent (≥ 0), -
 * @returns {number} D(θ), 1/rad
 */
export function spreadingCosN(theta, thetaP, n) {
  if (!(n >= 0)) throw new RangeError('n must be non-negative');
  const a = wrapPi(theta - thetaP);
  if (Math.abs(a) > Math.PI / 2) return 0;
  return (gammaFn(1 + n / 2) / (Math.sqrt(Math.PI) * gammaFn(0.5 + n / 2))) * Math.pow(Math.max(Math.cos(a), 0), n);
}

/**
 * cos-2s directional spreading function, D(θ) = Γ(s + 1) / [2√π Γ(s + ½)] cos^{2s}((θ − θp)/2), |θ − θp| ≤ π;
 * ∫ D dθ = 1.
 * Reference: Longuet-Higgins, Cartwright & Smith (1963); Mitsuyasu et al. (1975); DNV-RP-C205 (2019) §3.5.8.6
 * (comparable to cosⁿ for s ≈ 2n + 1).
 * @param {number} theta direction, rad
 * @param {number} thetaP main direction, rad
 * @param {number} s spreading parameter (≥ 0; keep below about 80 to avoid overflow of the gamma function), -
 * @returns {number} D(θ), 1/rad
 */
export function spreadingCos2s(theta, thetaP, s) {
  if (!(s >= 0) || s > 85) throw new RangeError('s must lie in [0, 85]');
  const a = wrapPi(theta - thetaP);
  return (gammaFn(s + 1) / (2 * Math.sqrt(Math.PI) * gammaFn(s + 0.5))) * Math.pow(Math.abs(Math.cos(a / 2)), 2 * s);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Moments and derived parameters                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Spectral moments m_n = ∫ ωⁿ S(ω) dω, n = −1, 0, 1, 2, 4, of a tabulated spectrum (trapezoidal rule).
 * Reference: DNV-RP-C205 (2019) §3.5.3.
 * Validity: omega ascending and positive; the table must cover the energetic range (no tail is added).
 * @param {ArrayLike<number>} omega angular frequencies, rad/s
 * @param {ArrayLike<number>} S spectral ordinates, m²·s/rad
 * @returns {{mMinus1:number,m0:number,m1:number,m2:number,m4:number}} m_n in m²·(rad/s)ⁿ
 */
export function tabulatedMoments(omega, S) {
  const n = omega.length;
  if (n < 2 || S.length !== n) throw new RangeError('omega and S must have the same length >= 2');
  let a = 0; let b = 0; let c = 0; let d = 0; let e = 0;
  for (let i = 1; i < n; i++) {
    const w0 = omega[i - 1]; const w1 = omega[i];
    if (!(w1 > w0) || !(w0 > 0)) throw new RangeError('omega must be positive and strictly ascending');
    const h = 0.5 * (w1 - w0);
    const s0 = S[i - 1]; const s1 = S[i];
    a += h * (s0 / w0 + s1 / w1);
    b += h * (s0 + s1);
    c += h * (s0 * w0 + s1 * w1);
    d += h * (s0 * w0 * w0 + s1 * w1 * w1);
    e += h * (s0 * w0 ** 4 + s1 * w1 ** 4);
  }
  return { mMinus1: a, m0: b, m1: c, m2: d, m4: e };
}

/**
 * Spectral moments of a spectrum given as a function, by composite Simpson integration on [wMin, wMax], with an
 * optional analytical ω⁻⁵ tail beyond wMax (exact asymptote of the Pierson–Moskowitz family) for n < 4.
 * m4 of an ω⁻⁵ spectrum diverges logarithmically and is always reported truncated at wMax.
 * Reference: DNV-RP-C205 (2019) §3.5.3.
 * @param {(omega:number)=>number} spectrum S(ω) in m²·s/rad
 * @param {object} o
 * @param {number} o.wMin lower limit, rad/s (> 0)
 * @param {number} o.wMax upper limit, rad/s
 * @param {number} [o.n=4000] number of intervals (rounded up to even)
 * @param {boolean} [o.tail=false] add the ω⁻⁵ tail beyond wMax
 * @returns {{mMinus1:number,m0:number,m1:number,m2:number,m4:number}} m_n in m²·(rad/s)ⁿ
 */
export function spectrumMoments(spectrum, { wMin, wMax, n = 4000, tail = false }) {
  reqPos('wMin', wMin); reqPos('wMax', wMax);
  if (!(wMax > wMin)) throw new RangeError('wMax must exceed wMin');
  const N = 2 * Math.ceil(n / 2);
  const h = (wMax - wMin) / N;
  let a = 0; let b = 0; let c = 0; let d = 0; let e = 0;
  let sEnd = 0;
  for (let i = 0; i <= N; i++) {
    const w = wMin + i * h;
    const s = spectrum(w);
    const wt = (i === 0 || i === N ? 1 : i % 2 ? 4 : 2) * (h / 3);
    const sw = s * wt;
    a += sw / w; b += sw; c += sw * w; d += sw * w * w; e += sw * w ** 4;
    if (i === N) sEnd = s;
  }
  if (tail) {
    // ∫_{wMax}^∞ ωⁿ S(wMax)(ω/wMax)⁻⁵ dω = S(wMax) wMax^{n+1} / (4 − n)
    a += sEnd / 5; b += (sEnd * wMax) / 4; c += (sEnd * wMax * wMax) / 3; d += (sEnd * wMax ** 3) / 2;
  }
  return { mMinus1: a, m0: b, m1: c, m2: d, m4: e };
}

/**
 * Sea-state parameters from spectral moments:
 *   Hm0 = 4 √m0,  Tz = Tm02 = 2π √(m0/m2),  Tm01 = 2π m0/m1,  Te = Tm-10 = 2π m₋₁/m0,  Tc = 2π √(m2/m4),
 *   bandwidth ε = √(1 − m2²/(m0 m4)),  narrowness ν = √(m0 m2/m1² − 1).
 * Reference: DNV-RP-C205 (2019) §3.5.3; Cartwright & Longuet-Higgins (1956); Longuet-Higgins (1975).
 * @param {{mMinus1:number,m0:number,m1:number,m2:number,m4:number}} m moments in ω units
 * @returns {{Hm0:number,Tz:number,Tm01:number,Te:number,Tc:number,epsilon:number,nu:number}} Hm0 m; periods s;
 *   epsilon, nu dimensionless
 */
export function spectralParameters(m) {
  if (!(m.m0 > 0)) throw new RangeError('m0 must be positive');
  return {
    Hm0: 4 * Math.sqrt(m.m0),
    Tz: TWO_PI * Math.sqrt(m.m0 / m.m2),
    Tm01: (TWO_PI * m.m0) / m.m1,
    Te: (TWO_PI * m.mMinus1) / m.m0,
    Tc: TWO_PI * Math.sqrt(m.m2 / m.m4),
    epsilon: Math.sqrt(Math.max(0, 1 - (m.m2 * m.m2) / (m.m0 * m.m4))),
    nu: Math.sqrt(Math.max(0, (m.m0 * m.m2) / (m.m1 * m.m1) - 1)),
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Short-term statistics                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Rayleigh probability that an individual wave height exceeds h, P = exp(−2 h²/Hs²) (Hs = 4 √m0).
 * Reference: Longuet-Higgins (1952); DNV-RP-C205 (2019) §3.5.9.
 * Validity: narrow-banded Gaussian sea in deep water; conservative (by roughly 5–10 % in height) for real seas.
 * @param {number} h wave height, m
 * @param {number} Hs significant wave height, m
 * @returns {number} exceedance probability, -
 */
export function rayleighExceedance(h, Hs) {
  reqPos('Hs', Hs);
  if (!(h > 0)) return 1;
  return Math.exp((-2 * h * h) / (Hs * Hs));
}

/**
 * Most probable (modal) maximum wave height in N Rayleigh-distributed waves, Hmax = Hs √(½ ln N).
 * Reference: Longuet-Higgins (1952); DNV-RP-C205 (2019) §3.5.11.
 * Validity: N large (≳ 100), narrow-banded sea.
 * @param {number} Hs significant wave height, m
 * @param {number} N number of waves, -
 * @returns {number} most probable maximum height, m
 */
export function mostProbableMaxHeight(Hs, N) {
  reqPos('Hs', Hs);
  if (!(N > 1)) throw new RangeError('N must exceed 1');
  return Hs * Math.sqrt(0.5 * Math.log(N));
}

/**
 * Expected maximum wave height in N Rayleigh-distributed waves.
 * 'asymptotic': E[Hmax] = Hs √(½ ln N) [1 + γE/(2 ln N)], γE = 0.5772 (Euler's constant);
 * 'exact': E[Hmax] = ∫₀^∞ {1 − [1 − exp(−2h²/Hs²)]^N} dh by numerical quadrature.
 * Reference: Longuet-Higgins (1952) Eq. 67; DNV-RP-C205 (2019) §3.5.11.
 * Validity: independent Rayleigh heights; the asymptotic form needs N ≳ 50.
 * @param {number} Hs significant wave height, m
 * @param {number} N number of waves, -
 * @param {'asymptotic'|'exact'} [method='asymptotic']
 * @returns {number} expected maximum height, m
 */
export function expectedMaxHeight(Hs, N, method = 'asymptotic') {
  reqPos('Hs', Hs);
  if (!(N >= 1)) throw new RangeError('N must be >= 1');
  if (method === 'asymptotic') {
    if (!(N > 1)) throw new RangeError('N must exceed 1 for the asymptotic form');
    const ln = Math.log(N);
    return Hs * Math.sqrt(0.5 * ln) * (1 + EULER / (2 * ln));
  }
  if (method !== 'exact') throw new RangeError("method must be 'asymptotic' or 'exact'");
  // integrate in x = h/Hs up to where the exceedance of the maximum is below 1e-18
  const xMax = Math.sqrt(0.5 * (Math.log(N) + 42));
  const n = 4000;
  const h = xMax / n;
  let s = 0;
  for (let i = 0; i <= n; i++) {
    const x = i * h;
    const f = 1 - Math.exp(N * Math.log1p(-Math.exp(-2 * x * x)));
    s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * f;
  }
  return (Hs * s * h) / 3;
}

/**
 * Maximum wave height in a stationary sea state of given duration: N = duration / Tz waves.
 * Reference: DNV-RP-C205 (2019) §3.5.11 (Rayleigh heights).
 * @param {object} o
 * @param {number} o.Hs significant wave height, m
 * @param {number} o.Tz mean zero-up-crossing period, s
 * @param {number} o.duration sea-state duration, s (e.g. 10800 for 3 hours)
 * @returns {{N:number,mostProbable:number,expected:number,mostProbableOverHs:number,expectedOverHs:number}}
 *   N number of waves; heights m; ratios dimensionless
 */
export function maxWaveHeight({ Hs, Tz, duration }) {
  reqPos('Tz', Tz); reqPos('duration', duration);
  const N = duration / Tz;
  const mp = mostProbableMaxHeight(Hs, N);
  const ex = expectedMaxHeight(Hs, N, 'asymptotic');
  return { N, mostProbable: mp, expected: ex, mostProbableOverHs: mp / Hs, expectedOverHs: ex / Hs };
}

/**
 * Forristall second-order crest-height distribution: P(crest > x) = exp[−(x/(α Hs))^β], with the Weibull
 * parameters fitted to second-order simulations in terms of the steepness S1 = 2π Hs/(g T1²) and the Ursell number
 * Ur = Hs/(k1² d³) (k1 = linear wavenumber for the mean period T1 = Tm01):
 *   long-crested (2D):   α = 0.3536 + 0.2892 S1 + 0.1060 Ur,   β = 2 − 2.1597 S1 + 0.0968 Ur²
 *   short-crested (3D):  α = 0.3536 + 0.2568 S1 + 0.0800 Ur,   β = 2 − 1.7912 S1 − 0.5302 Ur + c Ur²,  c = 0.284
 * Reference: Forristall, G.Z. (2000) "Wave crest distributions: observations and second-order theory",
 * J. Phys. Oceanogr. 30(8), 1931–1943; DNV-RP-C205 (2019) §3.5.10.4. The last 3D coefficient is printed as 0.284
 * in DNV-RP-C205 and appears as 0.2824 in some transcriptions of the paper; it is exposed as `beta3Ur2`.
 * Validity: S1 < about 0.15 and Ur < about 1 (range of the simulations); crests measured from the mean level.
 * @param {object} o
 * @param {number} o.Hs significant wave height, m
 * @param {number} o.T1 mean wave period Tm01, s
 * @param {number} o.d water depth, m (Infinity allowed)
 * @param {'2D'|'3D'} [o.type='3D'] long- or short-crested parameterisation
 * @param {number} [o.N] number of waves for the extreme-crest fields
 * @param {number} [o.beta3Ur2=0.284] Ur² coefficient of the 3D β
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{alpha:number,beta:number,S1:number,Ur:number,k1:number,scale:number,mostProbableMax:number,
 *   expectedMax:number}} alpha, beta, S1, Ur dimensionless; k1 rad/m; scale = α Hs m; mostProbableMax =
 *   α Hs (ln N)^{1/β} and expectedMax ≈ mostProbableMax [1 + γE/(β ln N)] in m (NaN when N is not given).
 */
export function forristallCrest({ Hs, T1, d, type = '3D', N, beta3Ur2 = 0.284, g = G }) {
  reqPos('Hs', Hs); reqPos('T1', T1);
  if (!(d > 0)) throw new RangeError('d must be positive');
  const k1 = kOf(TWO_PI / T1, d, g);
  const S1 = (TWO_PI * Hs) / (g * T1 * T1);
  const Ur = Number.isFinite(d) ? Hs / (k1 * k1 * d * d * d) : 0;
  let alpha; let beta;
  if (type === '2D') {
    alpha = 0.3536 + 0.2892 * S1 + 0.1060 * Ur;
    beta = 2 - 2.1597 * S1 + 0.0968 * Ur * Ur;
  } else if (type === '3D') {
    alpha = 0.3536 + 0.2568 * S1 + 0.0800 * Ur;
    beta = 2 - 1.7912 * S1 - 0.5302 * Ur + beta3Ur2 * Ur * Ur;
  } else throw new RangeError("type must be '2D' or '3D'");
  let mp = NaN; let ex = NaN;
  if (N !== undefined && N !== null) {
    if (!(N > 1)) throw new RangeError('N must exceed 1');
    const ln = Math.log(N);
    mp = alpha * Hs * Math.pow(ln, 1 / beta);
    ex = mp * (1 + EULER / (beta * ln));
  }
  return { alpha, beta, S1, Ur, k1, scale: alpha * Hs, mostProbableMax: mp, expectedMax: ex };
}

/**
 * Exceedance probability of a crest height under the Forristall distribution.
 * Reference and validity: see `forristallCrest`.
 * @param {number} x crest height above mean level, m
 * @param {{alpha:number,beta:number,scale:number}} params result of `forristallCrest`
 * @returns {number} P(crest > x), -
 */
export function forristallExceedance(x, params) {
  if (!(x > 0)) return 1;
  return Math.exp(-Math.pow(x / params.scale, params.beta));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Near-bed wave-induced velocity                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Wave-induced near-bed velocity from a surface spectrum by first-order theory:
 *   S_UU(ω) = G²(ω) S_ηη(ω),  G(ω) = ω / sinh(k d),  M_n = ∫ ωⁿ S_UU dω,
 *   Us = 2 √M0 (significant velocity amplitude),  Tu = 2π √(M0/M2) (mean zero-up-crossing period).
 * The integrals are evaluated numerically, so any spectrum may be used (the DNV design charts are for JONSWAP).
 * Reference: DNV-RP-F109 (2021) §3.4.3, Eqs. 3.8–3.13.
 * Validity: linear wave theory, no spreading reduction (apply the directional reduction factor separately).
 * @param {(omega:number)=>number} spectrum surface-elevation spectrum S(ω), m²·s/rad
 * @param {object} o
 * @param {number} o.d water depth, m
 * @param {number} o.wMin lower integration limit, rad/s
 * @param {number} o.wMax upper integration limit, rad/s
 * @param {number} [o.n=2000] number of Simpson intervals
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{Us:number,Tu:number,M0:number,M2:number}} Us m/s; Tu s; M0 (m/s)²; M2 (m/s)²·(rad/s)²
 */
export function nearBedVelocity(spectrum, { d, wMin, wMax, n = 2000, g = G }) {
  reqPos('d', d); reqPos('wMin', wMin); reqPos('wMax', wMax);
  if (!(wMax > wMin)) throw new RangeError('wMax must exceed wMin');
  const N = 2 * Math.ceil(n / 2);
  const h = (wMax - wMin) / N;
  let M0 = 0; let M2 = 0;
  for (let i = 0; i <= N; i++) {
    const w = wMin + i * h;
    const kd = kOf(w, d, g) * d;
    const Gt = kd > 350 ? 0 : w / Math.sinh(kd);
    const v = (i === 0 || i === N ? 1 : i % 2 ? 4 : 2) * (h / 3) * Gt * Gt * spectrum(w);
    M0 += v; M2 += v * w * w;
  }
  return { Us: 2 * Math.sqrt(M0), Tu: M0 > 0 ? TWO_PI * Math.sqrt(M0 / M2) : NaN, M0, M2 };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Time-domain synthesis and regular-wave equivalents                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Random long-crested sea-surface realisation by superposition of harmonic components with deterministic
 * amplitudes a_i = √(2 S(ω_i) Δω) and uniformly random phases; one frequency is drawn at random inside each
 * equal-width band so that the record does not repeat after 2π/Δω.
 *   η(x, t) = Σ a_i cos(ω_i t − k_i x + φ_i)
 * Reference: DNV-RP-C205 (2019) §3.3.2.1; Tucker, Challenor & Carter (1984) for the limitations of the
 * deterministic-amplitude method.
 * Validity: linear Gaussian sea; use at least ~200 components for a 3-hour record. Deterministic for a seed.
 * @param {(omega:number)=>number} spectrum S(ω), m²·s/rad
 * @param {object} o
 * @param {number} o.wMin lowest frequency, rad/s
 * @param {number} o.wMax highest frequency, rad/s
 * @param {number} o.duration record length, s
 * @param {number} o.dt time step, s
 * @param {number} [o.nFreq=200] number of components
 * @param {number} [o.seed=1] seed of the random generator
 * @param {number} [o.x=0] position along the propagation direction, m
 * @param {number} [o.d=Infinity] depth used for the wavenumbers when x ≠ 0, m
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{t:Float64Array,eta:Float64Array,components:{omega:Float64Array,amp:Float64Array,phase:Float64Array,
 *   k:Float64Array}}} t s; eta m; component omega rad/s, amp m, phase rad, k rad/m
 */
export function synthesizeSeaState(spectrum, { wMin, wMax, duration, dt, nFreq = 200, seed = 1, x = 0, d = Infinity, g = G }) {
  reqPos('wMin', wMin); reqPos('wMax', wMax); reqPos('duration', duration); reqPos('dt', dt);
  if (!(wMax > wMin)) throw new RangeError('wMax must exceed wMin');
  if (!Number.isInteger(nFreq) || nFreq < 1) throw new RangeError('nFreq must be a positive integer');
  const rng = makeRng(seed);
  const dw = (wMax - wMin) / nFreq;
  const omega = new Float64Array(nFreq); const amp = new Float64Array(nFreq);
  const phase = new Float64Array(nFreq); const kk = new Float64Array(nFreq);
  for (let i = 0; i < nFreq; i++) {
    const wc = wMin + (i + 0.5) * dw;
    omega[i] = wMin + (i + rng.next()) * dw;
    amp[i] = Math.sqrt(2 * spectrum(wc) * dw);
    phase[i] = TWO_PI * rng.next();
    kk[i] = kOf(omega[i], d, g);
  }
  const n = Math.floor(duration / dt + 1e-9) + 1;
  const t = new Float64Array(n); const eta = new Float64Array(n);
  // per-component rotation recurrence avoids one cosine evaluation per sample and component
  const cr = new Float64Array(nFreq); const ci = new Float64Array(nFreq);
  const rr = new Float64Array(nFreq); const ri = new Float64Array(nFreq);
  for (let i = 0; i < nFreq; i++) {
    const p0 = phase[i] - kk[i] * x;
    cr[i] = amp[i] * Math.cos(p0); ci[i] = amp[i] * Math.sin(p0);
    rr[i] = Math.cos(omega[i] * dt); ri[i] = Math.sin(omega[i] * dt);
  }
  for (let j = 0; j < n; j++) {
    t[j] = j * dt;
    let s = 0;
    for (let i = 0; i < nFreq; i++) {
      s += cr[i];
      const a = cr[i] * rr[i] - ci[i] * ri[i];
      ci[i] = cr[i] * ri[i] + ci[i] * rr[i];
      cr[i] = a;
    }
    eta[j] = s;
    // re-seed the recurrence periodically to stop round-off drift on long records
    if ((j & 1023) === 1023) {
      const tn = (j + 1) * dt;
      for (let i = 0; i < nFreq; i++) {
        const p = omega[i] * tn + phase[i] - kk[i] * x;
        cr[i] = amp[i] * Math.cos(p); ci[i] = amp[i] * Math.sin(p);
      }
    }
  }
  return { t, eta, components: { omega, amp, phase, k: kk } };
}

/**
 * Regular-wave equivalents of an irregular sea state:
 *   energy-equivalent wave  H = Hs/√2 (same mean energy ρ g Hs²/16), T = Te if given, else Tp;
 *   significant wave        H = Hs, T = Tp;
 *   design (maximum) wave   H = most probable maximum in `duration` (Rayleigh), with the associated-period range
 *                           √(cLow Hmax) ≤ T ≤ √(cHigh Hmax), cLow = 6.5 s²/m, cHigh = 11 s²/m.
 * Reference: energy equivalence from E = ρ g H²/8; maximum wave per DNV-RP-C205 (2019) §3.5.11; period range per
 * NORSOK N-003 (2017) §6.2.2.3 (stated there for the 100-year wave; coefficients exposed as arguments).
 * @param {object} o
 * @param {number} o.Hs significant wave height, m
 * @param {number} o.Tp peak period, s
 * @param {number} [o.Tz] zero-up-crossing period, s (default from the JONSWAP fit with the DNV γ rule)
 * @param {number} [o.Te] energy period, s
 * @param {number} [o.duration=10800] sea-state duration, s
 * @param {number} [o.cLow=6.5] lower period coefficient, s²/m
 * @param {number} [o.cHigh=11] upper period coefficient, s²/m
 * @returns {{energy:{H:number,T:number},significant:{H:number,T:number},design:{H:number,TLow:number,THigh:number,N:number}}}
 *   heights m; periods s; N number of waves
 */
export function regularWaveEquivalents({ Hs, Tp, Tz, Te, duration = 10800, cLow = 6.5, cHigh = 11 }) {
  reqPos('Hs', Hs); reqPos('Tp', Tp);
  const tz = Tz ?? Tp * jonswapPeriodRatios(jonswapGamma(Hs, Tp)).TzOverTp;
  const m = maxWaveHeight({ Hs, Tz: tz, duration });
  return {
    energy: { H: Hs / Math.SQRT2, T: Te ?? Tp },
    significant: { H: Hs, T: Tp },
    design: { H: m.mostProbable, TLow: Math.sqrt(cLow * m.mostProbable), THigh: Math.sqrt(cHigh * m.mostProbable), N: m.N },
  };
}
