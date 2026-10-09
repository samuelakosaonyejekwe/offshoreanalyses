// Axial load-transfer curves for piles: shaft t–z and tip Q–z.
// Curves share the interface of py.js: { type, pu, p(z), k(z), meta } where p is the mobilised unit shaft
// friction t [Pa] (or tip load Q [N]) at local pile displacement z [m], and k = dp/dz.

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

function tableCurve(type, xs, fs, pu, odd, meta) {
  const n = xs.length;
  const f = (x) => {
    if (x >= xs[n - 1]) return fs[n - 1];
    let i = 1;
    while (x > xs[i]) i++;
    return fs[i - 1] + ((fs[i] - fs[i - 1]) * (x - xs[i - 1])) / (xs[i] - xs[i - 1]);
  };
  const df = (x) => {
    if (x >= xs[n - 1]) return 0;
    let i = 1;
    while (x > xs[i]) i++;
    return (fs[i] - fs[i - 1]) / (xs[i] - xs[i - 1]);
  };
  return {
    type, pu, meta,
    p: (x) => (x >= 0 ? f(x) : odd ? -f(-x) : 0),
    k: (x) => (x >= 0 ? df(x) : odd ? df(-x) : 0),
  };
}

/** API clay t–z table: local displacement / pile diameter and mobilised fraction t/tmax (peak at z/D = 0.01). */
export const API_TZ_CLAY = Object.freeze({ zOverD: [0, 0.0016, 0.0031, 0.0057, 0.0080, 0.0100, 0.0200], tOverTmax: [0, 0.30, 0.50, 0.75, 0.90, 1.00, null] });

/** API Q–z table: tip displacement / pile diameter and mobilised fraction Q/Qp (full mobilisation at z/D = 0.1). */
export const API_QZ = Object.freeze({ zOverD: [0, 0.002, 0.013, 0.042, 0.073, 0.100], qOverQp: [0, 0.25, 0.50, 0.75, 0.90, 1.00] });

/**
 * Clay t–z curve (tabulated): z/D = 0.0016, 0.0031, 0.0057, 0.0080, 0.0100 → t/tmax = 0.30, 0.50, 0.75, 0.90,
 * 1.00; at z/D = 0.0200 and beyond t = residualRatio · tmax. Linear between points, odd-symmetric (tension).
 * Reference: API RP 2GEO (2011) §8.4 (API RP 2A-WSD §6.7.2). Validity: non-carbonate
 * clays; residual ratio 0.70–0.90 depending on soil stress–strain behaviour, stress history, installation and
 * loading sequence.
 * @param {{ tmax: number, D: number, residualRatio?: number }} o tmax [Pa] peak unit shaft friction; D [m] pile
 *   diameter; residualRatio [-] tres/tmax in [0.7, 0.9] by the standard, any value in (0, 1] accepted (default 0.9)
 * @returns {object} curve { type: 'tz-clay', pu: tmax [Pa], p(z) [Pa], k(z) [Pa/m], meta: { zPeak [m] } }
 */
export function tzClay({ tmax, D, residualRatio = 0.9 }) {
  pos('tmax', tmax, true); pos('D', D);
  if (!(residualRatio > 0 && residualRatio <= 1)) throw new RangeError('residualRatio: must be in (0, 1]');
  const xs = API_TZ_CLAY.zOverD.map((r) => r * D);
  const fs = API_TZ_CLAY.tOverTmax.map((r) => (r === null ? residualRatio : r) * tmax);
  return tableCurve('tz-clay', xs, fs, tmax, true, { zPeak: 0.01 * D });
}

/**
 * Sand t–z curve: linear from the origin to tmax at z = zPeak, constant beyond; odd-symmetric.
 * Reference: API RP 2GEO (2011) §8.4 (API RP 2A-WSD §6.7.2: 0.1 in). Validity: non-carbonate sands.
 * @param {{ tmax: number, zPeak?: number }} o tmax [Pa]; zPeak [m] displacement at full mobilisation (default
 *   0.0025 m; the US-customary edition gives 0.1 in = 0.00254 m)
 * @returns {object} curve { type: 'tz-sand', pu: tmax [Pa], p, k, meta: { zPeak [m] } }
 */
export function tzSand({ tmax, zPeak = 0.0025 }) {
  pos('tmax', tmax, true); pos('zPeak', zPeak);
  return tableCurve('tz-sand', [0, zPeak], [0, tmax], tmax, true, { zPeak });
}

/**
 * Tip load–displacement (Q–z) curve: z/D = 0.002, 0.013, 0.042, 0.073, 0.100 → Q/Qp = 0.25, 0.50, 0.75, 0.90,
 * 1.00; constant beyond. No resistance in uplift (z < 0).
 * Reference: API RP 2GEO (2011) §8.4 (API RP 2A-WSD §6.7.3). Validity: sands and clays.
 * @param {{ Qp: number, D: number }} o Qp [N] total end bearing; D [m] pile diameter
 * @returns {object} curve { type: 'qz', pu: Qp [N], p(z) [N], k(z) [N/m], meta: { zFull [m] } }
 */
export function qz({ Qp, D }) {
  pos('Qp', Qp, true); pos('D', D);
  return tableCurve('qz', API_QZ.zOverD.map((r) => r * D), API_QZ.qOverQp.map((r) => r * Qp), Qp, false, { zFull: 0.1 * D });
}

/**
 * Hyperbolic load-transfer curve: p = z / (1/ki + Rf z / pmax), capped at pmax (the cap is active only for
 * Rf < 1). With Rf = 1 the curve approaches pmax asymptotically.
 * Reference: hyperbolic stress–strain form of Kondner (1963), commonly used for t–z (e.g. Kraft, Ray & Kagawa
 * 1981). Validity: generic alternative; ki and Rf are user inputs.
 * @param {{ pmax: number, ki: number, Rf?: number, odd?: boolean }} o pmax [Pa or N] ultimate value; ki [Pa/m or
 *   N/m] initial stiffness; Rf [-] failure ratio in (0, 1] (default 1); odd [bool] true for shaft friction
 *   (default), false for tip resistance in compression only
 * @returns {object} curve { type: 'hyperbolic', pu: pmax, p, k, meta: { ki, Rf } }
 */
export function hyperbolicCurve({ pmax, ki, Rf = 1, odd = true }) {
  pos('pmax', pmax, true); pos('ki', ki);
  if (!(Rf > 0 && Rf <= 1)) throw new RangeError('Rf: must be in (0, 1]');
  const f = (z) => (pmax === 0 ? 0 : Math.min(z / (1 / ki + (Rf * z) / pmax), pmax));
  const df = (z) => {
    if (pmax === 0) return 0;
    const d = 1 / ki + (Rf * z) / pmax;
    return z / d >= pmax ? 0 : 1 / (ki * d * d);
  };
  return { type: 'hyperbolic', pu: pmax, meta: { ki, Rf }, p: (z) => (z >= 0 ? f(z) : odd ? -f(-z) : 0), k: (z) => (z >= 0 ? df(z) : odd ? df(-z) : 0) };
}

/**
 * Scale the ordinate of a load-transfer curve, e.g. cyclic degradation of shaft friction or a set-up factor.
 * Reference: none (generic factor; the value is an engineering input). Validity: factor >= 0.
 * @param {object} curve curve object
 * @param {number} factor multiplier on resistance [-]
 * @returns {object} scaled curve
 */
export function scaleCurve(curve, factor) {
  pos('factor', factor, true);
  return { type: curve.type, pu: curve.pu * factor, meta: { ...curve.meta, factor }, p: (z) => factor * curve.p(z), k: (z) => factor * curve.k(z) };
}

/**
 * Sample a load-transfer curve for plotting.
 * @param {object} curve curve object
 * @param {{ zMax: number, n?: number, zMin?: number }} o zMax [m]; n [-] points (default 61); zMin [m] (default 0)
 * @returns {{ z: number[], p: number[] }} z [m], p [Pa or N]
 */
export function sampleCurve(curve, { zMax, n = 61, zMin = 0 }) {
  if (!(zMax > zMin)) throw new RangeError('zMax: must exceed zMin');
  if (!(n >= 2)) throw new RangeError('n: must be >= 2');
  const z = [], p = [];
  for (let i = 0; i < n; i++) { const zz = zMin + ((zMax - zMin) * i) / (n - 1); z.push(zz); p.push(curve.p(zz)); }
  return { z, p };
}
