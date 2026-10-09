// Frequency-domain fatigue from a one-sided stress power spectral density S(f) [Pa²/Hz], f in Hz.
// S–N constants are passed in the neutral form { m1, a1, m2, a2, sKnee } of fatigue/sn.js (stress in Pa).

import { gamma, gammaP, gammaQ } from '../math/special.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function isSeries(x) {
  return Array.isArray(x) || ArrayBuffer.isView(x);
}
function checkSn(sn) {
  need(sn && sn.m1 > 0 && sn.a1 > 0, 'sn must provide positive m1 and a1');
}
function bilinear(sn) {
  return sn.m2 !== null && sn.m2 !== undefined && sn.sKnee > 0;
}

/**
 * Spectral moment m_n = ∫ fⁿ·S(f) df by the trapezoidal rule (any real order n).
 * @param {ArrayLike<number>} f frequencies, ascending [Hz]
 * @param {ArrayLike<number>} S one-sided PSD [unit²/Hz]
 * @param {number} order n [-]
 * @returns {number} moment [unit²·Hzⁿ]
 */
export function spectralMoment(f, S, order) {
  need(isSeries(f) && isSeries(S) && f.length === S.length && f.length >= 2, 'f and S must be arrays of equal length >= 2');
  const pw = (x) => (order === 0 ? 1 : x > 0 ? x ** order : 0);
  let sum = 0;
  let prev = pw(f[0]) * S[0];
  for (let i = 1; i < f.length; i++) {
    need(f[i] > f[i - 1], 'f must be strictly ascending');
    const cur = pw(f[i]) * S[i];
    sum += 0.5 * (cur + prev) * (f[i] - f[i - 1]);
    prev = cur;
  }
  return sum;
}

/**
 * Spectral moments, crossing/peak rates and bandwidth parameters of a stationary Gaussian process.
 * ν0 = √(m2/m0) mean zero up-crossing rate, νp = √(m4/m2) mean peak rate, α1 = m1/√(m0·m2),
 * α2 = m2/√(m0·m4) (irregularity factor), ε = √(1 − α2²), Vanmarcke δ = √(1 − α1²).
 * Reference: Rice (1944); Vanmarcke (1972); DNV-RP-C205 (spectral moments).
 * @param {ArrayLike<number>} f [Hz]
 * @param {ArrayLike<number>} S [unit²/Hz]
 * @returns {{m0:number,m1:number,m2:number,m4:number,sigma:number,nu0:number,nuP:number,alpha1:number,alpha2:number,
 *   epsilon:number,delta:number}} sigma in signal units, rates in Hz, bandwidth parameters dimensionless
 */
export function spectralParameters(f, S) {
  const m0 = spectralMoment(f, S, 0);
  const m1 = spectralMoment(f, S, 1);
  const m2 = spectralMoment(f, S, 2);
  const m4 = spectralMoment(f, S, 4);
  need(m0 > 0 && m2 > 0 && m4 > 0, 'S must have positive area');
  const alpha1 = Math.min(1, m1 / Math.sqrt(m0 * m2));
  const alpha2 = Math.min(1, m2 / Math.sqrt(m0 * m4));
  return {
    m0, m1, m2, m4, sigma: Math.sqrt(m0), nu0: Math.sqrt(m2 / m0), nuP: Math.sqrt(m4 / m2), alpha1, alpha2,
    epsilon: Math.sqrt(Math.max(0, 1 - alpha2 * alpha2)), delta: Math.sqrt(Math.max(0, 1 - alpha1 * alpha1)),
  };
}

/**
 * Narrow-band (Rayleigh) fatigue damage rate of a Gaussian process with standard deviation σ and cycle rate ν:
 *   single slope:  d = ν/a·(2√2σ)^m·Γ(1 + m/2)
 *   two slopes:    d = ν·[ (2√2σ)^m1/a1·Γ(1 + m1/2; x) + (2√2σ)^m2/a2·γ(1 + m2/2; x) ],  x = (S1/(2√2σ))²
 * Reference: Bendat (1964); DNV-RP-C203 (closed-form damage of narrow-band response) / DNV-OS-E301 Ch.2 Sec.2 [6.3].
 * Validity: narrow-banded Gaussian response (conservative for wide-band).
 * @param {object} o
 * @param {number} o.sigma standard deviation of stress [Pa]
 * @param {number} o.rate cycle rate ν (normally ν0) [Hz]
 * @param {object} o.sn S–N constants
 * @returns {number} damage per second [1/s]
 */
export function narrowBandDamageRate({ sigma, rate, sn }) {
  need(sigma >= 0, 'sigma must be non-negative');
  need(rate >= 0, 'rate must be non-negative');
  checkSn(sn);
  if (sigma === 0 || rate === 0) return 0;
  const q = 2 * Math.SQRT2 * sigma;
  const g1 = 1 + sn.m1 / 2;
  if (!bilinear(sn)) return (rate * q ** sn.m1 * gamma(g1)) / sn.a1;
  const x = (sn.sKnee / q) ** 2;
  const g2 = 1 + sn.m2 / 2;
  return rate * ((q ** sn.m1 * gamma(g1) * gammaQ(g1, x)) / sn.a1 + (q ** sn.m2 * gamma(g2) * gammaP(g2, x)) / sn.a2);
}

/**
 * Wirsching–Light rainflow correction factor on narrow-band damage:
 *   λ = a(m) + [1 − a(m)]·(1 − ε)^b(m),  a = 0.926 − 0.033·m,  b = 1.587·m − 2.323.
 * Reference: Wirsching & Light (1980) J. Struct. Div. ASCE 106(ST7), 1593–1607. Validity: single-slope S–N,
 * m ≈ 3–6.
 * @param {number} epsilon bandwidth parameter ε = √(1 − α2²) [-]
 * @param {number} m S–N inverse slope [-]
 * @returns {number} λ [-]
 */
export function wirschingLightFactor(epsilon, m) {
  need(epsilon >= 0 && epsilon <= 1, 'epsilon must be in [0, 1]');
  need(m > 0, 'm must be positive');
  const a = 0.926 - 0.033 * m;
  const b = 1.587 * m - 2.323;
  return a + (1 - a) * (1 - epsilon) ** b;
}

/**
 * Dirlik parameters and rainflow-range probability density.
 *   xm = (m1/m0)·√(m2/m4), γ = α2, D1 = 2(xm − γ²)/(1 + γ²), R = (γ − xm − D1²)/(1 − γ − D1 + D1²),
 *   D2 = (1 − γ − D1 + D1²)/(1 − R), D3 = 1 − D1 − D2, Q = 1.25(γ − D3 − D2·R)/D1,
 *   p(S) = [D1/Q·e^(−Z/Q) + D2·Z/R²·e^(−Z²/(2R²)) + D3·Z·e^(−Z²/2)]/(2√m0),  Z = S/(2√m0).
 * Reference: Dirlik (1985) PhD thesis, University of Warwick.
 * @param {{m0:number,m1:number,m2:number,m4:number}} mom spectral moments
 * @returns {{D1:number,D2:number,D3:number,Q:number,R:number,xm:number,gamma:number,pdf:(function(number):number)}}
 *   pdf(S) in 1/(signal unit), S = stress range
 */
export function dirlikParameters(mom) {
  const { m0, m1, m2, m4 } = mom;
  need(m0 > 0 && m2 > 0 && m4 > 0, 'moments must be positive');
  const g = Math.min(m2 / Math.sqrt(m0 * m4), 1 - 1e-9);
  const xm = (m1 / m0) * Math.sqrt(m2 / m4);
  const D1 = (2 * (xm - g * g)) / (1 + g * g);
  const R = (g - xm - D1 * D1) / (1 - g - D1 + D1 * D1);
  const D2 = (1 - g - D1 + D1 * D1) / (1 - R);
  const D3 = 1 - D1 - D2;
  const Q = (1.25 * (g - D3 - D2 * R)) / D1;
  const s0 = 2 * Math.sqrt(m0);
  // Degenerate (very narrow) spectra drive D1 -> 0 and Q -> 0/0; the first term then vanishes.
  const useExp = D1 > 1e-12 && Q > 0;
  const pdf = (S) => {
    const Z = S / s0;
    let p = D3 * Z * Math.exp(-0.5 * Z * Z);
    if (useExp) p += (D1 / Q) * Math.exp(-Z / Q);
    if (D2 !== 0 && R !== 0) p += ((D2 * Z) / (R * R)) * Math.exp((-Z * Z) / (2 * R * R));
    return p / s0;
  };
  return { D1, D2, D3, Q, R, xm, gamma: g, pdf };
}

// E[N(S)^-1] for a range PDF against a (possibly bilinear) S–N curve, by composite Simpson on [0, sMax].
function expectedDamagePerCycle(pdf, sn, sMax, n = 4000) {
  const h = sMax / n;
  let sum = 0;
  for (let i = 1; i <= n; i++) {
    const s = i * h;
    const N = bilinear(sn) && s < sn.sKnee ? sn.a2 * s ** -sn.m2 : sn.a1 * s ** -sn.m1;
    const w = i === n ? 1 : i % 2 === 1 ? 4 : 2;
    sum += (w * pdf(s)) / N;
  }
  return (sum * h) / 3;
}

/**
 * Dirlik fatigue damage rate. Single slope uses the closed form
 *   d = νp/a·(2√m0)^m·[ D1·Q^m·Γ(1 + m) + (√2)^m·Γ(1 + m/2)·(D2·|R|^m + D3) ];
 * a two-slope curve is integrated numerically over the Dirlik range density.
 * Reference: Dirlik (1985). Validity: stationary Gaussian stress, any bandwidth.
 * @param {object} o
 * @param {{m0:number,m1:number,m2:number,m4:number}} o.moments spectral moments of stress [Pa², Hz]
 * @param {object} o.sn S–N constants
 * @returns {number} damage per second [1/s]
 */
export function dirlikDamageRate({ moments, sn }) {
  checkSn(sn);
  const p = dirlikParameters(moments);
  const nuP = Math.sqrt(moments.m4 / moments.m2);
  const s0 = 2 * Math.sqrt(moments.m0);
  if (bilinear(sn)) return nuP * expectedDamagePerCycle(p.pdf, sn, 12 * s0 * Math.max(1, p.Q > 0 ? p.Q : 1));
  const m = sn.m1;
  let t = Math.SQRT2 ** m * gamma(1 + m / 2) * (p.D2 * Math.abs(p.R) ** m + p.D3);
  if (p.D1 > 1e-12 && p.Q > 0) t += p.D1 * p.Q ** m * gamma(1 + m);
  return (nuP * s0 ** m * t) / sn.a1;
}

/**
 * Tovo–Benasciutti damage rate: d = [b + (1 − b)·α2^(m−1)]·d_NB with d_NB the narrow-band rate at ν0 and
 *   b = (α1 − α2)·[1.112·(1 + α1·α2 − (α1 + α2))·e^(2.11·α2) + (α1 − α2)]/(α2 − 1)².
 * Reference: Tovo (2002) Int. J. Fatigue 24, 1137–1147; Benasciutti & Tovo (2005) Int. J. Fatigue 27, 867–877.
 * Validity: single-slope S–N (for a two-slope curve the weighting uses m1 and the narrow-band term is exact).
 * @param {object} o
 * @param {ArrayLike<number>} o.f [Hz]
 * @param {ArrayLike<number>} o.S stress PSD [Pa²/Hz]
 * @param {object} o.sn S–N constants
 * @returns {{damageRate:number, b:number, narrowBandRate:number}} rates in 1/s
 */
export function tovoBenasciuttiDamageRate({ f, S, sn }) {
  checkSn(sn);
  const p = spectralParameters(f, S);
  const nb = narrowBandDamageRate({ sigma: p.sigma, rate: p.nu0, sn });
  const { alpha1: a1, alpha2: a2 } = p;
  const b = a2 >= 1 - 1e-9 ? 1
    : Math.min(1, ((a1 - a2) * (1.112 * (1 + a1 * a2 - (a1 + a2)) * Math.exp(2.11 * a2) + (a1 - a2))) / (a2 - 1) ** 2);
  return { damageRate: (b + (1 - b) * a2 ** (sn.m1 - 1)) * nb, b, narrowBandRate: nb };
}

/**
 * α0.75 method: d = α0.75²·d_NB with α0.75 = m0.75/√(m0·m1.5).
 * Reference: Benasciutti & Tovo (2006) Probabilistic Engineering Mechanics 21, 287–299 (after Lutes & Larsen).
 * @param {object} o
 * @param {ArrayLike<number>} o.f [Hz]
 * @param {ArrayLike<number>} o.S [Pa²/Hz]
 * @param {object} o.sn
 * @returns {number} damage per second [1/s]
 */
export function alpha075DamageRate({ f, S, sn }) {
  const p = spectralParameters(f, S);
  const a = spectralMoment(f, S, 0.75) / Math.sqrt(p.m0 * spectralMoment(f, S, 1.5));
  return a * a * narrowBandDamageRate({ sigma: p.sigma, rate: p.nu0, sn });
}

/**
 * Fatigue damage from a stress PSD by a chosen spectral method.
 * @param {object} o
 * @param {ArrayLike<number>} o.f frequencies [Hz]
 * @param {ArrayLike<number>} o.S one-sided stress PSD [Pa²/Hz]
 * @param {object} o.sn S–N constants (stress in Pa)
 * @param {number} [o.duration=1] exposure time [s]
 * @param {string} [o.method='dirlik'] 'narrow-band' | 'wirsching-light' | 'dirlik' | 'tovo-benasciutti' | 'alpha075'
 * @returns {{damage:number, damageRate:number, cycles:number, parameters:object}} damage [-], rate [1/s],
 *   cycles = ν0·duration (νp·duration for Dirlik) [cycles]
 */
export function spectralDamage({ f, S, sn, duration = 1, method = 'dirlik' }) {
  need(duration >= 0, 'duration must be non-negative');
  checkSn(sn);
  const p = spectralParameters(f, S);
  let rate;
  let cyc = p.nu0;
  switch (method) {
    case 'narrow-band': rate = narrowBandDamageRate({ sigma: p.sigma, rate: p.nu0, sn }); break;
    case 'wirsching-light':
      rate = wirschingLightFactor(p.epsilon, sn.m1) * narrowBandDamageRate({ sigma: p.sigma, rate: p.nu0, sn });
      break;
    case 'dirlik': rate = dirlikDamageRate({ moments: p, sn }); cyc = p.nuP; break;
    case 'tovo-benasciutti': rate = tovoBenasciuttiDamageRate({ f, S, sn }).damageRate; break;
    case 'alpha075': rate = alpha075DamageRate({ f, S, sn }); break;
    default: throw new RangeError(`method: unknown spectral fatigue method "${method}"`);
  }
  return { damage: rate * duration, damageRate: rate, cycles: cyc * duration, parameters: p };
}

/**
 * Bimodal (low-frequency + wave-frequency) fatigue damage of two narrow-band Gaussian components.
 * Methods:
 *  'simple-summation'   d = d_L + d_W (DNV-OS-E301: acceptable when one component dominates);
 *  'combined-spectrum'  narrow-band with σY = √(σL² + σW²), νY = √(λL·νL² + λW·νW²) (conservative);
 *  'dual-narrow-band'   combined-spectrum damage times the Jiao–Moan correction
 *      ρ = (νe/νY)·[λL^(m/2+2)·(1 − √(λW/λL)) + √(π·λL·λW)·m·Γ((1+m)/2)/Γ(1+m/2)] + (νW/νY)·λW^(m/2),
 *      νe = √(λL²·νL² + λL·λW·νW²·δW²), λ = σ²/σY².
 * Reference: Jiao & Moan (1990) Probabilistic Engineering Mechanics 5(2), 76–83; DNV-OS-E301 Ch.2 Sec.2 [6.3].
 * Validity: single-slope S–N; the Jiao–Moan envelope formula assumes the low-frequency part carries the larger
 * share (λL >= λW) — outside that range the function falls back to the combined-spectrum value (ρ capped at 1).
 * @param {object} o
 * @param {number} o.sigmaLF standard deviation of the low-frequency stress [Pa]
 * @param {number} o.nuLF zero up-crossing rate of the low-frequency part [Hz]
 * @param {number} o.sigmaWF standard deviation of the wave-frequency stress [Pa]
 * @param {number} o.nuWF zero up-crossing rate of the wave-frequency part [Hz]
 * @param {object} o.sn S–N constants (single slope used: m1, a1)
 * @param {number} [o.deltaWF=0.1] Vanmarcke bandwidth of the wave-frequency part [-] (0.1 in DNV-OS-E301)
 * @param {string} [o.method='dual-narrow-band']
 * @returns {{damageRate:number, rho:number, sigmaY:number, nuY:number, rateLF:number, rateWF:number,
 *   combinedSpectrumRate:number}} rates in 1/s, sigmaY in Pa, nuY in Hz
 */
export function bimodalDamageRate({ sigmaLF, nuLF, sigmaWF, nuWF, sn, deltaWF = 0.1, method = 'dual-narrow-band' }) {
  need(sigmaLF >= 0 && sigmaWF >= 0, 'sigmaLF and sigmaWF must be non-negative');
  need(nuLF > 0 && nuWF > 0, 'nuLF and nuWF must be positive');
  need(deltaWF >= 0 && deltaWF <= 1, 'deltaWF must be in [0, 1]');
  checkSn(sn);
  const single = { m1: sn.m1, a1: sn.a1, m2: null, a2: null, sKnee: 0 };
  const m = sn.m1;
  const rateLF = narrowBandDamageRate({ sigma: sigmaLF, rate: nuLF, sn: single });
  const rateWF = narrowBandDamageRate({ sigma: sigmaWF, rate: nuWF, sn: single });
  const varY = sigmaLF * sigmaLF + sigmaWF * sigmaWF;
  const out = { damageRate: 0, rho: 1, sigmaY: Math.sqrt(varY), nuY: 0, rateLF, rateWF, combinedSpectrumRate: 0 };
  if (varY === 0) return out;
  const lL = (sigmaLF * sigmaLF) / varY;
  const lW = 1 - lL;
  const nuY = Math.sqrt(lL * nuLF * nuLF + lW * nuWF * nuWF);
  const cs = narrowBandDamageRate({ sigma: out.sigmaY, rate: nuY, sn: single });
  out.nuY = nuY;
  out.combinedSpectrumRate = cs;
  if (method === 'simple-summation') { out.damageRate = rateLF + rateWF; return out; }
  if (method === 'combined-spectrum') { out.damageRate = cs; return out; }
  need(method === 'dual-narrow-band', `method: unknown bimodal method "${method}"`);
  let rho = 1;
  if (lL > 0) {
    const nuE = Math.sqrt(lL * lL * nuLF * nuLF + lL * lW * nuWF * nuWF * deltaWF * deltaWF);
    rho = (nuE / nuY) * (lL ** (m / 2 + 2) * (1 - Math.sqrt(lW / lL))
      + (Math.sqrt(Math.PI * lL * lW) * m * gamma((1 + m) / 2)) / gamma(1 + m / 2)) + (nuWF / nuY) * lW ** (m / 2);
    if (!(rho > 0) || rho > 1) rho = Math.min(1, Math.max(rho, 0));
  }
  out.rho = rho;
  out.damageRate = rho * cs;
  return out;
}

/**
 * Stress PSD from a wave (or other load) spectrum and a stress transfer function: Sσ(f) = |H(f)|²·Sη(f).
 * @param {ArrayLike<number>} waveSpectrum Sη(f) [m²/Hz] (or load unit²/Hz)
 * @param {ArrayLike<number>} transfer |H(f)| [Pa/m] (stress per unit load amplitude)
 * @returns {Float64Array} stress PSD [Pa²/Hz]
 */
export function stressSpectrum(waveSpectrum, transfer) {
  need(isSeries(waveSpectrum) && isSeries(transfer) && waveSpectrum.length === transfer.length,
    'waveSpectrum and transfer must be arrays of equal length');
  const out = new Float64Array(waveSpectrum.length);
  for (let i = 0; i < out.length; i++) out[i] = transfer[i] * transfer[i] * waveSpectrum[i];
  return out;
}

/**
 * Long-term spectral fatigue over a scatter diagram: D = T·Σ pᵢ·dᵢ, where dᵢ is the short-term damage rate of
 * sea state (and heading) i from its stress PSD.
 * Reference: DNV-RP-C206 (spectral fatigue methodology); DNV-RP-C203 (frequency-domain fatigue).
 * @param {object} o
 * @param {object[]} o.states each { probability [-], f [Hz], S: stress PSD [Pa²/Hz] } or
 *   { probability, f, waveSpectrum [m²/Hz], transfer |H(f)| [Pa/m] }; optional `label`
 * @param {object} o.sn S–N constants
 * @param {number} o.duration exposure time [s] (e.g. one year = 31 557 600 s)
 * @param {string} [o.method='dirlik'] see {@link spectralDamage}
 * @returns {{damage:number, life:number, cycles:number, contributions: {label:(string|undefined),
 *   probability:number, damage:number, sigma:number, nu0:number}[]}} damage over `duration` [-], life = duration /
 *   damage [s], cycles [cycles]; per-state sigma [Pa] and nu0 [Hz]
 */
export function longTermSpectralDamage({ states, sn, duration, method = 'dirlik' }) {
  need(Array.isArray(states) && states.length > 0, 'states must be a non-empty array');
  need(duration > 0, 'duration must be positive');
  let total = 0;
  let cycles = 0;
  const contributions = [];
  for (let i = 0; i < states.length; i++) {
    const st = states[i];
    need(st.probability >= 0, `states[${i}].probability must be non-negative`);
    const S = st.S || stressSpectrum(st.waveSpectrum, st.transfer);
    let d = 0;
    let sigma = 0;
    let nu0 = 0;
    if (st.probability > 0 && spectralMoment(st.f, S, 0) > 0) {
      const r = spectralDamage({ f: st.f, S, sn, duration: duration * st.probability, method });
      d = r.damage;
      sigma = r.parameters.sigma;
      nu0 = r.parameters.nu0;
      cycles += r.cycles;
    }
    total += d;
    contributions.push({ label: st.label, probability: st.probability, damage: d, sigma, nu0 });
  }
  return { damage: total, life: total > 0 ? duration / total : Infinity, cycles, contributions };
}
