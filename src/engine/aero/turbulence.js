/**
 * IEC 61400-1 wind conditions: turbulence models, extreme deterministic events, Kaimal and von Karman spectra,
 * spatial coherence, and Veers (Sandia) synthesis of a correlated turbulent field on a rotor-plane grid with
 * reduction to a rotor-effective wind speed.
 *
 * References: IEC 61400-1:2005+A1:2010 (ed. 3) and IEC 61400-1:2019 (ed. 4), clause 6.3 and Annex B/C (Annex C in
 * ed. 4); IEC 61400-3:2009 / IEC 61400-3-1:2019 clause 6 for offshore; Veers, P.S. (1988), "Three-dimensional wind
 * simulation", SAND88-0152.
 */

import { ifft, fft } from '../math/fft.js';
import { makeRng } from '../math/rng.js';

const TWO_PI = 2 * Math.PI;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/** Reference turbulence intensity Iref at 15 m/s by turbulence category (IEC 61400-1 Table 1; A+ added in ed. 4) [-]. */
export const TURBULENCE_CLASSES = { 'A+': 0.18, A: 0.16, B: 0.14, C: 0.12 };

/** Reference wind speed Vref (10-min mean, 50-year) by wind turbine class (IEC 61400-1 Table 1) [m/s]; Vave = 0.2 Vref. */
export const WIND_CLASSES = { I: 50, II: 42.5, III: 37.5 };

/** Normal wind profile power-law exponent: 0.2 (IEC 61400-1, 6.3.1.2), 0.14 offshore (IEC 61400-3, 6.3) [-]. */
export const NWP_EXPONENT = { onshore: 0.2, offshore: 0.14 };

/**
 * Longitudinal turbulence scale parameter: Lambda1 = 0.7 z for z <= 60 m, 42 m above (IEC 61400-1 Eq. 5).
 * @param {number} z hub height [m]
 * @returns {number} Lambda1 [m]
 */
export function turbulenceScale(z) {
  need(z > 0, 'z', 'must be > 0');
  return z <= 60 ? 0.7 * z : 42;
}

/**
 * Normal turbulence model: sigma1 = Iref (0.75 Vhub + b), b = 5.6 m/s (IEC 61400-1 Eq. 11), the representative
 * (90 % quantile) longitudinal standard deviation.
 * @param {object} o
 * @param {number} o.Vhub hub-height 10-min mean wind speed [m/s]
 * @param {number} o.Iref reference turbulence intensity [-]
 * @returns {number} sigma1 [m/s]
 */
export function ntmSigma(o) {
  need(o && o.Vhub >= 0, 'Vhub', 'must be >= 0');
  need(o.Iref > 0 && o.Iref < 1, 'Iref', 'must be in (0, 1)');
  return o.Iref * (0.75 * o.Vhub + 5.6);
}

/**
 * Extreme turbulence model: sigma1 = c Iref [0.072 (Vave/c + 3)(Vhub/c - 4) + 10], c = 2 m/s (IEC 61400-1 Eq. 19).
 * @param {object} o
 * @param {number} o.Vhub [m/s]
 * @param {number} o.Iref [-]
 * @param {number} o.Vave annual mean hub-height wind speed [m/s] (0.2 Vref for standard classes)
 * @returns {number} sigma1 [m/s]
 */
export function etmSigma(o) {
  need(o && o.Vhub >= 0, 'Vhub', 'must be >= 0');
  need(o.Iref > 0 && o.Iref < 1, 'Iref', 'must be in (0, 1)');
  need(o.Vave > 0, 'Vave', 'must be > 0');
  const c = 2;
  return c * o.Iref * (0.072 * (o.Vave / c + 3) * (o.Vhub / c - 4) + 10);
}

/**
 * Extreme wind speed model (IEC 61400-1, 6.3.2.1).
 * Steady: Ve50(z) = 1.4 Vref (z/zhub)^0.11, Ve1 = 0.8 Ve50 (3-s gusts). Turbulent: V50(z) = Vref (z/zhub)^0.11,
 * V1 = 0.8 V50 (10-min means) with sigma1 = 0.11 Vhub.
 * @param {object} o
 * @param {number} o.Vref reference wind speed [m/s]
 * @param {number} o.z height [m]
 * @param {number} o.hubHeight [m]
 * @returns {{Ve50:number, Ve1:number, V50:number, V1:number, sigma50:number, sigma1yr:number}} speeds [m/s];
 *   turbulence standard deviations for the 50-year and 1-year turbulent models (0.11 x hub-height mean) [m/s]
 */
export function ewm(o) {
  need(o && o.Vref > 0, 'Vref', 'must be > 0');
  need(o.z > 0 && o.hubHeight > 0, 'z/hubHeight', 'must be > 0');
  const s = Math.pow(o.z / o.hubHeight, 0.11);
  return { Ve50: 1.4 * o.Vref * s, Ve1: 0.8 * 1.4 * o.Vref * s, V50: o.Vref * s, V1: 0.8 * o.Vref * s, sigma50: 0.11 * o.Vref, sigma1yr: 0.11 * 0.8 * o.Vref };
}

/**
 * Extreme operating gust magnitude: Vgust = min{ 1.35 (Ve1 - Vhub), 3.3 sigma1 / (1 + 0.1 D/Lambda1) }
 * (IEC 61400-1 Eq. 17), Ve1 the 1-year extreme (3-s) speed at hub height.
 * @param {object} o
 * @param {number} o.Vhub [m/s] @param {number} o.sigma1 NTM standard deviation [m/s]
 * @param {number} o.D rotor diameter [m] @param {number} o.lambda1 turbulence scale [m]
 * @param {number} o.Ve1 1-year extreme wind speed at hub height [m/s]
 * @returns {number} Vgust [m/s]
 */
export function eogGust(o) {
  need(o && o.D > 0 && o.lambda1 > 0, 'D/lambda1', 'must be > 0');
  need(o.sigma1 >= 0, 'sigma1', 'must be >= 0');
  return Math.min(1.35 * (o.Ve1 - o.Vhub), 3.3 * o.sigma1 / (1 + 0.1 * o.D / o.lambda1));
}

/**
 * Extreme operating gust time shape ("Mexican hat"):
 * V(t) = Vhub - 0.37 Vgust sin(3 pi t/T) (1 - cos(2 pi t/T)) for 0 <= t <= T, T = 10.5 s; Vhub otherwise
 * (IEC 61400-1 Eq. 18). The peak, at t = T/2, is Vhub + 0.74 Vgust.
 * @param {number} t time since gust start [s]
 * @param {object} o
 * @param {number} o.Vhub [m/s] @param {number} o.Vgust [m/s] @param {number} [o.T=10.5] [s]
 * @returns {number} wind speed [m/s]
 */
export function eogSpeed(t, o) {
  const T = o.T ?? 10.5;
  if (t < 0 || t > T) return o.Vhub;
  return o.Vhub - 0.37 * o.Vgust * Math.sin(3 * Math.PI * t / T) * (1 - Math.cos(TWO_PI * t / T));
}

/**
 * Extreme direction change: theta_e = +-4 atan( sigma1 / (Vhub (1 + 0.1 D/Lambda1)) ), transient
 * theta(t) = 0.5 theta_e (1 - cos(pi t/T)), T = 6 s (IEC 61400-1 Eqs. 21-22).
 * @param {object} o
 * @param {number} o.Vhub [m/s] (> 0) @param {number} o.sigma1 [m/s] @param {number} o.D [m] @param {number} o.lambda1 [m]
 * @returns {{thetaEDeg:number, T:number, directionDeg:function(number):number}} magnitude [deg] (apply both signs),
 *   duration [s], direction history t [s] -> [deg]
 */
export function edc(o) {
  need(o && o.Vhub > 0, 'Vhub', 'must be > 0');
  const th = 4 * Math.atan(o.sigma1 / (o.Vhub * (1 + 0.1 * o.D / o.lambda1))) * 180 / Math.PI;
  const T = 6;
  return { thetaEDeg: th, T, directionDeg: (t) => (t <= 0 ? 0 : (t >= T ? th : 0.5 * th * (1 - Math.cos(Math.PI * t / T)))) };
}

/**
 * Extreme coherent gust with direction change: Vcg = 15 m/s, V(t) = Vhub + 0.5 Vcg (1 - cos(pi t/T)), T = 10 s;
 * theta_cg = 180 deg for Vhub < 4 m/s, 720 deg / Vhub for 4 m/s <= Vhub <= Vref; theta(t) = 0.5 theta_cg
 * (1 - cos(pi t/T)) (IEC 61400-1 Eqs. 23-26).
 * @param {object} o
 * @param {number} o.Vhub [m/s] (> 0)
 * @returns {{Vcg:number, thetaCgDeg:number, T:number, speed:function(number):number, directionDeg:function(number):number}}
 *   gust [m/s], direction change [deg], rise time [s], histories t [s] (since gust start) -> speed [m/s], direction [deg]
 */
export function ecd(o) {
  need(o && o.Vhub > 0, 'Vhub', 'must be > 0');
  const Vcg = 15; const T = 10;
  const th = o.Vhub < 4 ? 180 : 720 / o.Vhub;
  const ramp = (t) => (t <= 0 ? 0 : (t >= T ? 1 : 0.5 * (1 - Math.cos(Math.PI * t / T))));
  return { Vcg, thetaCgDeg: th, T, speed: (t) => o.Vhub + Vcg * ramp(t), directionDeg: (t) => th * ramp(t) };
}

/**
 * Extreme wind shear, transient increment at an offset from the hub (IEC 61400-1 Eqs. 27-28):
 * dV(t) = +-(offset/D) (2.5 + 0.2 beta sigma1 (D/Lambda1)^(1/4)) (1 - cos(2 pi t/T)), beta = 6.4, T = 12 s,
 * with offset = z - zhub for the vertical and y for the horizontal shear; added to the normal wind profile.
 * @param {object} o
 * @param {number} o.offset vertical (z - zhub) or lateral (y) distance from the hub [m]
 * @param {number} o.t time since event start [s]
 * @param {number} o.sigma1 [m/s] @param {number} o.D [m] @param {number} o.lambda1 [m]
 * @returns {number} speed increment [m/s] (positive-sign event; negate for the opposite sign)
 */
export function ewsIncrement(o) {
  const T = 12; const beta = 6.4;
  if (o.t < 0 || o.t > T) return 0;
  return (o.offset / o.D) * (2.5 + 0.2 * beta * o.sigma1 * Math.pow(o.D / o.lambda1, 0.25)) * (1 - Math.cos(TWO_PI * o.t / T));
}

/**
 * Sample a deterministic wind event into a time series (event inserted into steady wind).
 * @param {object} o
 * @param {function(number):number} o.fn value as a function of time since event start [s]
 * @param {number} o.dt [s] @param {number} o.duration [s] @param {number} [o.tStart=0] event start [s]
 * @returns {{t:Float64Array, v:Float64Array}} time [s], values (units of fn)
 */
export function eventSeries(o) {
  need(o && o.dt > 0 && o.duration > 0, 'dt/duration', 'must be > 0');
  const n = Math.round(o.duration / o.dt) + 1;
  const t = new Float64Array(n); const v = new Float64Array(n);
  for (let i = 0; i < n; i++) { t[i] = i * o.dt; v[i] = o.fn(t[i] - (o.tStart ?? 0)); }
  return { t, v };
}

/**
 * Kaimal spectral parameters at hub height (IEC 61400-1 Table B.1 / C.1): sigma_k = (1, 0.8, 0.5) sigma1 and
 * integral scales L_k = (8.1, 2.7, 0.66) Lambda1 for the longitudinal, lateral and upward components.
 * @param {object} o
 * @param {number} o.hubHeight [m] @param {number} o.sigma1 [m/s]
 * @returns {{lambda1:number, u:{sigma:number,L:number}, v:{sigma:number,L:number}, w:{sigma:number,L:number}, Lc:number}}
 *   scales [m], standard deviations [m/s], coherence scale Lc = 8.1 Lambda1 [m]
 */
export function kaimalParams(o) {
  const l1 = turbulenceScale(o.hubHeight);
  return { lambda1: l1, u: { sigma: o.sigma1, L: 8.1 * l1 }, v: { sigma: 0.8 * o.sigma1, L: 2.7 * l1 }, w: { sigma: 0.5 * o.sigma1, L: 0.66 * l1 }, Lc: 8.1 * l1 };
}

/**
 * Kaimal one-sided power spectral density: S(f) = sigma^2 (4 L/V) / (1 + 6 f L/V)^(5/3) (IEC 61400-1 Eq. B.14 /
 * C.14). Integrates to sigma^2 over 0..infinity exactly.
 * @param {number} f frequency [Hz] (>= 0)
 * @param {object} o
 * @param {number} o.Vhub mean wind speed [m/s] (> 0) @param {number} o.sigma standard deviation [m/s] @param {number} o.L integral scale [m]
 * @returns {number} S [m2/s2/Hz]
 */
export function kaimalSpectrum(f, o) {
  return o.sigma * o.sigma * 4 * o.L / o.Vhub / Math.pow(1 + 6 * f * o.L / o.Vhub, 5 / 3);
}

/**
 * von Karman isotropic-turbulence spectra (one-sided), x = f L / V with L the longitudinal integral scale:
 *   longitudinal  S_u = sigma^2 (4 L/V) / (1 + 70.8 x^2)^(5/6)
 *   lateral/vertical  S_v = S_w = sigma^2 (2 L/V) (1 + 188.8 x^2) / (1 + 70.8 x^2)^(11/6)
 * Reference: von Karman (1948); ESDU 85020; IEC 61400-1 ed. 2 Annex B (rounded constants 71 and 189). Each
 * integrates to sigma^2. Isotropic: all components share sigma.
 * @param {number} f frequency [Hz]
 * @param {object} o
 * @param {number} o.V mean speed [m/s] @param {number} o.sigma [m/s] @param {number} o.L longitudinal integral scale [m]
 * @param {'u'|'v'|'w'} [o.component='u']
 * @returns {number} S [m2/s2/Hz]
 */
export function vonKarmanSpectrum(f, o) {
  const x = f * o.L / o.V;
  const b = 1 + 70.8 * x * x;
  if ((o.component ?? 'u') === 'u') return o.sigma * o.sigma * 4 * o.L / o.V / Math.pow(b, 5 / 6);
  return o.sigma * o.sigma * 2 * o.L / o.V * (1 + 188.8 * x * x) / Math.pow(b, 11 / 6);
}

/**
 * IEC exponential coherence of the longitudinal component between two points separated by r in the rotor plane:
 * Coh(r, f) = exp( -a sqrt( (f r / Vhub)^2 + (0.12 r / Lc)^2 ) ), a = 12, Lc = 8.1 Lambda1 (IEC 61400-1 Eq. B.16 /
 * C.16). This is the magnitude of the (real, in-phase) co-coherence, not its square.
 * @param {number} f frequency [Hz] @param {number} r separation [m]
 * @param {object} o
 * @param {number} o.Vhub [m/s] @param {number} o.Lc coherence scale [m] @param {number} [o.a=12] decrement [-]
 * @returns {number} coherence [-]
 */
export function coherenceIec(f, r, o) {
  const a = o.a ?? 12;
  return Math.exp(-a * Math.sqrt((f * r / o.Vhub) ** 2 + (0.12 * r / o.Lc) ** 2));
}

function evenCount(duration, dt) {
  let n = Math.round(duration / dt);
  if (n % 2) n += 1;
  return n;
}

/**
 * Single-point turbulence time series from a one-sided spectrum by random-phase inverse FFT:
 * x(t) = sum_m sqrt(2 S(f_m) df) cos(2 pi f_m t + phi_m), f_m = m/T. Periodic with period T; zero mean; the sample
 * variance equals sum S df exactly (the band above the Nyquist frequency and below 1/T is not represented).
 * @param {object} o
 * @param {function(number):number} o.spectrum S(f) [unit^2/Hz], f in Hz
 * @param {number} o.duration [s] @param {number} o.dt [s] @param {number} o.seed integer seed
 * @returns {{t:Float64Array, x:Float64Array, variance:number}} series and its (exact) variance
 */
export function synthesizeSeries(o) {
  need(o && o.duration > 0 && o.dt > 0, 'duration/dt', 'must be > 0');
  const n = evenCount(o.duration, o.dt);
  const T = n * o.dt;
  const df = 1 / T;
  const rng = makeRng(o.seed ?? 1);
  const re = new Float64Array(n); const im = new Float64Array(n);
  let variance = 0;
  for (let m = 1; m < n / 2; m++) {
    const S = o.spectrum(m * df);
    const amp = Math.sqrt(0.5 * S * df) * n;
    const ph = TWO_PI * rng.next();
    re[m] = amp * Math.cos(ph); im[m] = amp * Math.sin(ph);
    re[n - m] = re[m]; im[n - m] = -im[m];
    variance += S * df;
  }
  const r = ifft(re, im);
  const t = new Float64Array(n); const x = new Float64Array(n);
  for (let i = 0; i < n; i++) { t[i] = i * o.dt; x[i] = r.re[i]; }
  return { t, x, variance };
}

/**
 * Veers (Sandia) method: longitudinal turbulence at the points of a rectangular rotor-plane grid with the Kaimal
 * spectrum (same at every point, hub-height parameters) and the IEC exponential coherence. For each frequency the
 * coherence matrix is factorised (Cholesky, C = H H^T) and the Fourier coefficient at point j is
 * sqrt(S df / 2) sum_k H_jk exp(i phi_k) with independent uniform phases; an inverse FFT gives the series.
 * Taylor's frozen-turbulence hypothesis is implied (no along-wind evolution). Cost O(nf * np^3 / 6).
 * @param {object} o
 * @param {number} o.Vhub mean hub-height speed [m/s] (> 0)
 * @param {number} o.sigma1 longitudinal standard deviation [m/s]
 * @param {number} o.hubHeight [m]
 * @param {number} o.width grid width [m] @param {number} o.height grid height [m] (centred on the hub)
 * @param {number} [o.ny=5] @param {number} [o.nz=5] grid points per side (ny*nz <= 121)
 * @param {number} o.duration [s] @param {number} o.dt [s] @param {number} o.seed integer seed
 * @param {number} [o.a=12] coherence decrement [-]
 * @param {number} [o.fMax] highest synthesised frequency [Hz]; default Nyquist
 * @param {boolean} [o.scaleToSigma=true] multiply all series by the common factor sigma1 / sqrt(sum S df), which
 *   restores the target standard deviation lost outside the represented band 1/T..Nyquist (for the Kaimal spectrum
 *   the band below 1/T holds roughly 2 L/(V T) of the variance); false keeps the spectral level exact
 * @returns {{t:Float64Array, y:number[], z:number[], points:Array<{y:number,z:number}>, u:Float64Array[], Vhub:number,
 *   hubHeight:number, dt:number, bandVariance:number, scale:number}} time [s]; lateral and vertical (relative to hub) grid coordinates [m]; fluctuation
 *   series u[j] [m/s] at points[j] (j = iz*ny + iy), zero mean; variance of the represented band before scaling
 *   [m2/s2] and the scale factor applied [-]
 */
export function synthesizeWindField(o) {
  need(o && o.Vhub > 0, 'Vhub', 'must be > 0');
  need(o.sigma1 >= 0, 'sigma1', 'must be >= 0');
  need(o.duration > 0 && o.dt > 0, 'duration/dt', 'must be > 0');
  const ny = o.ny ?? 5; const nz = o.nz ?? 5;
  need(Number.isInteger(ny) && Number.isInteger(nz) && ny >= 1 && nz >= 1 && ny * nz <= 121, 'ny/nz', 'integers with ny*nz <= 121');
  const np = ny * nz;
  const kp = kaimalParams({ hubHeight: o.hubHeight, sigma1: o.sigma1 });
  const a = o.a ?? 12;
  const ys = []; const zs = []; const points = [];
  for (let i = 0; i < ny; i++) ys.push(ny === 1 ? 0 : -0.5 * o.width + o.width * i / (ny - 1));
  for (let i = 0; i < nz; i++) zs.push(nz === 1 ? 0 : -0.5 * o.height + o.height * i / (nz - 1));
  for (let iz = 0; iz < nz; iz++) for (let iy = 0; iy < ny; iy++) points.push({ y: ys[iy], z: zs[iz] });
  const dist = new Float64Array(np * np);
  for (let i = 0; i < np; i++) for (let j = 0; j < np; j++) dist[i * np + j] = Math.hypot(points[i].y - points[j].y, points[i].z - points[j].z);
  const n = evenCount(o.duration, o.dt);
  const T = n * o.dt; const df = 1 / T;
  const mMax = Math.min(n / 2 - 1, o.fMax ? Math.floor(o.fMax / df) : n / 2 - 1);
  const rng = makeRng(o.seed ?? 1);
  const re = []; const im = [];
  for (let j = 0; j < np; j++) { re.push(new Float64Array(n)); im.push(new Float64Array(n)); }
  const H = new Float64Array(np * np);
  const cph = new Float64Array(np); const sph = new Float64Array(np);
  const c2 = (0.12 / kp.Lc) ** 2;
  let band = 0;
  for (let m = 1; m <= mMax; m++) band += kaimalSpectrum(m * df, { Vhub: o.Vhub, sigma: o.sigma1, L: kp.u.L }) * df;
  const scale = (o.scaleToSigma ?? true) && band > 0 ? o.sigma1 / Math.sqrt(band) : 1;
  for (let m = 1; m <= mMax; m++) {
    const f = m * df;
    const amp = scale * Math.sqrt(0.5 * kaimalSpectrum(f, { Vhub: o.Vhub, sigma: o.sigma1, L: kp.u.L }) * df) * n;
    const c1 = (f / o.Vhub) ** 2;
    const kf = a * Math.sqrt(c1 + c2);
    // in-place Cholesky of the coherence matrix (lower triangle in H)
    for (let i = 0; i < np; i++) {
      for (let j = 0; j <= i; j++) {
        let s = Math.exp(-kf * dist[i * np + j]);
        for (let k = 0; k < j; k++) s -= H[i * np + k] * H[j * np + k];
        if (i === j) H[i * np + i] = Math.sqrt(s > 1e-14 ? s : 1e-14);
        else H[i * np + j] = s / H[j * np + j];
      }
    }
    for (let k = 0; k < np; k++) { const ph = TWO_PI * rng.next(); cph[k] = Math.cos(ph); sph[k] = Math.sin(ph); }
    for (let j = 0; j < np; j++) {
      let xr = 0; let xi = 0;
      for (let k = 0; k <= j; k++) { const h = H[j * np + k]; xr += h * cph[k]; xi += h * sph[k]; }
      re[j][m] = amp * xr; im[j][m] = amp * xi;
      re[j][n - m] = amp * xr; im[j][n - m] = -amp * xi;
    }
  }
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) t[i] = i * o.dt;
  const u = [];
  for (let j = 0; j < np; j++) {
    const r = ifft(re[j], im[j]);
    u.push(Float64Array.from(r.re));
  }
  return { t, y: ys, z: zs, points, u, Vhub: o.Vhub, hubHeight: o.hubHeight, dt: o.dt, bandVariance: band, scale };
}

/**
 * Rotor-effective wind speed from a synthesised field: average over the grid points inside the rotor disc of
 * (mean profile + fluctuation); the mean profile is Vhub (1 + z/hubHeight)^shearExp.
 * `power` selects the averaging norm: 1 (thrust-like, linear), 2 or 3 (energy-flux equivalent,
 * (mean V^p)^(1/p)).
 * @param {object} field synthesizeWindField result
 * @param {object} o
 * @param {number} o.R rotor radius [m]
 * @param {number} [o.shearExp=0] @param {number} [o.power=1]
 * @returns {{t:Float64Array, u:Float64Array, nPoints:number}} time [s], rotor-effective speed [m/s], points used
 */
export function rotorEffectiveWind(field, o) {
  need(o && o.R > 0, 'R', 'must be > 0');
  const p = o.power ?? 1;
  const idx = []; const mean = [];
  field.points.forEach((pt, j) => {
    if (Math.hypot(pt.y, pt.z) <= o.R * (1 + 1e-9)) {
      idx.push(j);
      mean.push(field.Vhub * Math.pow(Math.max(1e-3, 1 + pt.z / field.hubHeight), o.shearExp ?? 0));
    }
  });
  need(idx.length >= 1, 'R', 'no grid point lies inside the rotor disc');
  const n = field.t.length;
  const u = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let k = 0; k < idx.length; k++) { const v = mean[k] + field.u[idx[k]][i]; s += p === 1 ? v : Math.pow(Math.max(0, v), p); }
    s /= idx.length;
    u[i] = p === 1 ? s : Math.pow(s, 1 / p);
  }
  return { t: field.t, u, nPoints: idx.length };
}

/**
 * Rotor-averaging filter: ratio of the spectrum of the disc-averaged longitudinal speed to the point spectrum,
 * |H(f)|^2 = (1/A^2) integral integral Coh(|x1 - x2|, f) dA1 dA2, evaluated numerically for the IEC exponential
 * coherence on a sunflower pattern of `nPts` equal-area points. Tends to its low-frequency plateau below one
 * (finite Lc) and falls roughly as (V / (f R))^2 at high frequency.
 * @param {number[]|Float64Array} f frequencies [Hz]
 * @param {object} o
 * @param {number} o.R rotor radius [m] @param {number} o.Vhub [m/s] @param {number} o.Lc coherence scale [m]
 * @param {number} [o.a=12] @param {number} [o.nPts=64]
 * @returns {Float64Array} |H|^2 [-] at each frequency
 */
export function rotorAveragingFilter(f, o) {
  need(o && o.R >= 0 && o.Vhub > 0 && o.Lc > 0, 'R/Vhub/Lc', 'must be > 0');
  const np = o.nPts ?? 64;
  const a = o.a ?? 12;
  const px = new Float64Array(np); const py = new Float64Array(np);
  const ga = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < np; i++) { const r = o.R * Math.sqrt((i + 0.5) / np); px[i] = r * Math.cos(i * ga); py[i] = r * Math.sin(i * ga); }
  const nd = np * (np - 1) / 2;
  const d = new Float64Array(nd);
  let q = 0;
  for (let i = 0; i < np; i++) for (let j = 0; j < i; j++) d[q++] = Math.hypot(px[i] - px[j], py[i] - py[j]);
  const c2 = (0.12 / o.Lc) ** 2;
  const one = (fr) => {
    const kf = a * Math.sqrt((fr / o.Vhub) ** 2 + c2);
    let s = 0;
    for (let i = 0; i < nd; i++) s += Math.exp(-kf * d[i]);
    return (np + 2 * s) / (np * np);
  };
  // tabulate on a log grid of frequency and interpolate (the filter is smooth)
  const nt = 48;
  const fLo = 1e-4; const fHi = Math.max(10, 50 * o.Vhub / Math.max(o.R, 1));
  const lf = new Float64Array(nt); const hv = new Float64Array(nt);
  for (let i = 0; i < nt; i++) { lf[i] = Math.log(fLo) + (Math.log(fHi) - Math.log(fLo)) * i / (nt - 1); hv[i] = one(Math.exp(lf[i])); }
  const h0 = one(0);
  const out = new Float64Array(f.length);
  for (let i = 0; i < f.length; i++) {
    const fr = f[i];
    if (fr <= fLo) { out[i] = h0 + (hv[0] - h0) * fr / fLo; continue; }
    if (fr >= fHi) { out[i] = hv[nt - 1]; continue; }
    const u = (Math.log(fr) - lf[0]) / (lf[1] - lf[0]);
    const k = Math.min(nt - 2, Math.floor(u));
    out[i] = hv[k] + (u - k) * (hv[k + 1] - hv[k]);
  }
  return out;
}

// Synthesis step of rotorEffectiveTurbulence [s]: Nyquist 5 Hz, far above the band the rotor average passes.
const RET_DT0 = 0.1;

/**
 * Rotor-effective turbulent wind speed generated directly in the frequency domain: single-point Kaimal spectrum
 * multiplied by the rotor-averaging filter, random phases. A fast substitute for the grid synthesis when only the
 * rotor-effective speed is needed (no information on shear fluctuation or blade-resolved turbulence).
 *
 * The realisation depends on the seed and the duration only, not on the output step. The series is synthesised
 * on a fixed grid — step 0.1 s, length the next power of two covering the duration (so the inverse FFT is a
 * radix-2 transform whatever the request), frequencies m/T0 with T0 that padded length — and then resampled to
 * `dt`: steps of 0.1/2^k s by zero-padding the spectrum (exact band-limited interpolation), multiples of the
 * synthesis step by picking samples, anything else by linear interpolation of the next finer 0.1/2^k s series.
 * Two runs that differ only in dt therefore see the same wind, which is what a time-step convergence check needs.
 * The returned window is recentred on Vhub and rescaled to the variance of the whole synthesised period (the
 * spectral sum), both evaluated on the 0.1 s synthesis grid; when duration/0.1 s is a power of two the window is
 * the whole period and neither correction does anything.
 * @param {object} o
 * @param {number} o.Vhub [m/s] @param {number} o.sigma1 [m/s] @param {number} o.hubHeight [m] @param {number} o.R [m]
 * @param {number} o.duration [s] @param {number} o.dt [s] @param {number} o.seed
 * @param {boolean} [o.scaleToSigma=true] scale so that the underlying point spectrum carries sigma1^2 within the
 *   represented band (see synthesizeWindField)
 * @returns {{t:Float64Array, u:Float64Array, sigmaEff:number}} speed (mean Vhub) [m/s] at t = 0, dt, 2 dt ...
 *   (an even number of samples, round(duration/dt)); sample standard deviation of the returned series [m/s]
 *   (equal to the root of the spectral sum on the 0.1 s grid, and within a few 1e-4 of it at other steps)
 */
export function rotorEffectiveTurbulence(o) {
  need(o && o.duration > 0 && o.dt > 0, 'duration/dt', 'must be > 0');
  const kp = kaimalParams({ hubHeight: o.hubHeight, sigma1: o.sigma1 });
  const n = evenCount(o.duration, o.dt);
  // fixed synthesis grid
  let n0 = 2;
  while (n0 * RET_DT0 < o.duration * (1 - 1e-12)) n0 *= 2;
  const df = 1 / (n0 * RET_DT0);
  const half = n0 / 2;
  const fr = new Float64Array(half);
  for (let m = 0; m < half; m++) fr[m] = m * df;
  const H2 = rotorAveragingFilter(fr, { R: o.R, Vhub: o.Vhub, Lc: kp.Lc });
  const Sk = new Float64Array(half);
  let band = 0;
  for (let m = 1; m < half; m++) { Sk[m] = kaimalSpectrum(fr[m], { Vhub: o.Vhub, sigma: o.sigma1, L: kp.u.L }); band += Sk[m] * df; }
  const sc2 = (o.scaleToSigma ?? true) && band > 0 ? o.sigma1 * o.sigma1 / band : 1;
  // finer output grid by zero-padding: step RET_DT0 / pad <= dt
  let pad = 1;
  while (RET_DT0 / pad > o.dt * (1 + 1e-9)) pad *= 2;
  const nf = n0 * pad;
  const dtf = RET_DT0 / pad;
  const rng = makeRng(o.seed ?? 1);
  const re = new Float64Array(nf); const im = new Float64Array(nf);
  for (let m = 1; m < half; m++) {
    const amp = Math.sqrt(0.5 * sc2 * Sk[m] * H2[m] * df) * nf;
    const ph = TWO_PI * rng.next();
    re[m] = amp * Math.cos(ph); im[m] = amp * Math.sin(ph);
    re[nf - m] = re[m]; im[nf - m] = -im[m];
  }
  const xf = ifft(re, im).re;
  // Mean and variance of the window, taken on the synthesis grid so that they are the same for every output step.
  // The window is recentred on Vhub and rescaled to the variance of the full period (the spectral sum), which a
  // window shorter than the period would otherwise miss in part.
  const nWin = Math.max(2, Math.min(n0, Math.round(o.duration / RET_DT0)));
  let mean = 0; let target = 0;
  for (let i = 0; i < nWin; i++) mean += xf[i * pad];
  mean /= nWin;
  let varW = 0;
  for (let i = 0; i < nWin; i++) { const d = xf[i * pad] - mean; varW += d * d; }
  varW /= nWin;
  for (let m = 1; m < half; m++) target += sc2 * Sk[m] * H2[m] * df;
  const gain = varW > 0 ? Math.sqrt(target / varW) : 1;
  const u = new Float64Array(n); const t = new Float64Array(n);
  const ratio = o.dt / dtf;
  const step = Math.round(ratio);
  const whole = Math.abs(ratio - step) < 1e-9 * ratio;
  for (let i = 0; i < n; i++) {
    t[i] = i * o.dt;
    let v;
    if (whole) v = xf[(i * step) % nf];
    else {
      const x = i * ratio; const j = Math.floor(x); const w = x - j;
      v = xf[j % nf] + w * (xf[(j + 1) % nf] - xf[j % nf]);
    }
    u[i] = o.Vhub + gain * (v - mean);
  }
  let mu = 0;
  for (let i = 0; i < n; i++) mu += u[i];
  mu /= n;
  let vs = 0;
  for (let i = 0; i < n; i++) vs += (u[i] - mu) ** 2;
  return { t, u, sigmaEff: Math.sqrt(vs / n) };
}

/**
 * Apply the rotor-averaging filter to a measured or synthesised single-point speed series (zero-phase, in the
 * frequency domain; the mean is preserved).
 * @param {object} o
 * @param {Float64Array|number[]} o.u point wind speed series [m/s] @param {number} o.dt [s]
 * @param {number} o.R [m] @param {number} o.Vhub [m/s] @param {number} o.Lc [m]
 * @returns {Float64Array} rotor-effective series [m/s]
 */
export function rotorEffectiveFromPoint(o) {
  const n = o.u.length;
  need(n >= 4, 'u', 'at least 4 samples required');
  const F = fft(Float64Array.from(o.u));
  const df = 1 / (n * o.dt);
  const fr = new Float64Array(n);
  for (let k = 0; k < n; k++) fr[k] = Math.min(k, n - k) * df;
  const H2 = rotorAveragingFilter(fr, { R: o.R, Vhub: o.Vhub, Lc: o.Lc });
  const re = new Float64Array(n); const im = new Float64Array(n);
  for (let k = 0; k < n; k++) { const h = k === 0 ? 1 : Math.sqrt(H2[k]); re[k] = F.re[k] * h; im[k] = F.im[k] * h; }
  return Float64Array.from(ifft(re, im).re);
}

/**
 * Offshore turbulence standard deviation from the Charnock sea-surface roughness (IEC 61400-3:2009, 6.3.1.3):
 *   z0 = (Ac / g) [ kappa Vhub / ln(zhub / z0) ]^2      (implicit, solved by fixed-point iteration)
 *   sigma1 = Vhub / ln(zhub / z0) + 1.28 * 1.44 * I15   (90 % quantile; I15 = average turbulence intensity at 15 m/s)
 * with Ac = 0.011 for open sea (0.034 near-coastal) and kappa = 0.4. Site data should be used where available.
 * @param {object} o
 * @param {number} o.Vhub [m/s] (> 0) @param {number} o.hubHeight [m]
 * @param {number} o.I15 average hub-height turbulence intensity at 15 m/s [-] (site value)
 * @param {number} [o.Ac=0.011] Charnock constant [-] @param {number} [o.kappa=0.4] @param {number} [o.g=9.80665]
 * @returns {{z0:number, sigma1:number, TI:number, uStar:number}} roughness [m], standard deviation [m/s], intensity [-],
 *   friction velocity [m/s]
 */
export function offshoreTurbulence(o) {
  need(o && o.Vhub > 0 && o.hubHeight > 0, 'Vhub/hubHeight', 'must be > 0');
  need(o.I15 >= 0 && o.I15 < 1, 'I15', 'must be in [0, 1)');
  const Ac = o.Ac ?? 0.011; const kappa = o.kappa ?? 0.4; const g = o.g ?? 9.80665;
  let z0 = 2e-4;
  for (let i = 0; i < 200; i++) {
    const zn = Ac / g * (kappa * o.Vhub / Math.log(o.hubHeight / z0)) ** 2;
    if (Math.abs(zn - z0) < 1e-14) { z0 = zn; break; }
    z0 = zn;
  }
  const ln = Math.log(o.hubHeight / z0);
  const sigma1 = o.Vhub / ln + 1.28 * 1.44 * o.I15;
  return { z0, sigma1, TI: sigma1 / o.Vhub, uStar: kappa * o.Vhub / ln };
}

/**
 * Peak factor of a stationary Gaussian process (Davenport 1964): g = sqrt(2 ln(nu T)) + 0.5772 / sqrt(2 ln(nu T)),
 * the expected maximum in standard deviations over duration T for mean zero-upcrossing rate nu.
 * @param {number} nu zero-upcrossing frequency [Hz] @param {number} T duration [s]
 * @returns {number} g [-]
 */
export function peakFactor(nu, T) {
  need(nu * T > 1, 'nu*T', 'must exceed 1');
  const b = Math.sqrt(2 * Math.log(nu * T));
  return b + 0.5772 / b;
}

/**
 * Gust factor G = Vgust / Vmean = 1 + g I for turbulence intensity I and peak factor g (statistical gust model;
 * g about 3 - 3.7 for a 3-s gust in 10 minutes, see peakFactor).
 * @param {number} intensity turbulence intensity sigma/V [-] @param {number} [g=3.5] peak factor [-]
 * @returns {number} gust factor [-]
 */
export function gustFactor(intensity, g = 3.5) {
  need(intensity >= 0, 'intensity', 'must be >= 0');
  return 1 + g * intensity;
}
