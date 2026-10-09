// Energy-yield assessment: power/thrust curves, wind statistics, gross and wake-affected farm AEP, loss cascade,
// uncertainty and exceedance levels, measure–correlate–predict, time-series production and extreme wind.
//
// Units: speeds m/s, power W, energy Wh (field names end in `Wh`), fractions (-). Hours per year default 8760.

import { gamma, erfinv } from '../math/special.js';
import { fitGumbel } from '../reliability/extreme.js';
import { farmFlow } from './wake.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}
function frac(v, name) {
  need(Number.isFinite(v) && v >= 0 && v <= 1, `${name} must be within [0, 1]`);
}

function interpTable(xs, ys, x) {
  const n = xs.length;
  if (x < xs[0] || x > xs[n - 1]) return 0;
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid;
    else hi = mid;
  }
  return ys[lo] + ((x - xs[lo]) / (xs[hi] - xs[lo])) * (ys[hi] - ys[lo]);
}

/**
 * Power of a tabulated curve at wind speed U (linear interpolation, zero outside the table).
 * @param {{u:number[], p:number[]}} curve wind speed (m/s), power (W)
 * @param {number} U wind speed (m/s)
 * @returns {number} power (W)
 */
export function curvePower(curve, U) {
  need(curve && curve.u && curve.p && curve.u.length === curve.p.length && curve.u.length > 1,
    'curve must be a table { u: [], p: [] } with at least two points');
  return interpTable(curve.u, curve.p, U);
}

/** Lower root a ∈ [0, 1/3] of Cp = 4 a (1 − a)² (ideal actuator disc). */
function inductionFromCp(Cp) {
  if (Cp <= 0) return 0;
  if (Cp >= 16 / 27) return 1 / 3;
  let lo = 0;
  let hi = 1 / 3;
  for (let i = 0; i < 60; i++) {
    const a = 0.5 * (lo + hi);
    if (4 * a * (1 - a) * (1 - a) < Cp) lo = a;
    else hi = a;
  }
  return 0.5 * (lo + hi);
}

/**
 * Generic pitch-regulated power and thrust curve: P = min(P_rated, ½ ρ A Cp U³) between cut-in and cut-out,
 * optionally blended near rated as P = (P_c^−n + P_r^−n)^(−1/n) (n = `smoothing`; Infinity gives the sharp knee).
 * Thrust coefficient from ideal actuator-disc momentum theory: the induction a solves
 * 4 a (1 − a)² = Cp_el / efficiency and Ct = 4 a (1 − a); above rated the rotor sheds power, so Cp, a and Ct fall.
 * Reference: Rankine–Froude momentum theory; Betz limit Cp ≤ 16/27. A screening curve, not a turbine model.
 * @param {object} o
 * @param {number} o.ratedPower rated electrical power (W)
 * @param {number} o.D rotor diameter (m)
 * @param {number} [o.Cp=0.45] electrical power coefficient below rated (-), ≤ 16/27 × efficiency
 * @param {number} [o.efficiency=0.94] drivetrain/electrical efficiency used to recover the aerodynamic Cp (-)
 * @param {number} [o.rho=1.225] air density (kg/m³)
 * @param {number} [o.cutIn=3] cut-in wind speed (m/s)
 * @param {number} [o.cutOut=25] cut-out wind speed (m/s)
 * @param {number} [o.smoothing=Infinity] knee blending exponent n (-); 8–12 mimics turbulence rounding
 * @param {number} [o.step=0.5] table step (m/s)
 * @returns {{u:number[], p:number[], cp:number[], ct:number[], ratedSpeed:number, ratedPower:number, D:number,
 *   rho:number}} table (m/s, W, -, -), sharp-knee rated speed (m/s)
 */
export function cubicPowerCurve({
  ratedPower,
  D,
  Cp = 0.45,
  efficiency = 0.94,
  rho = 1.225,
  cutIn = 3,
  cutOut = 25,
  smoothing = Infinity,
  step = 0.5,
}) {
  finitePos(ratedPower, 'ratedPower');
  finitePos(D, 'D');
  finitePos(rho, 'rho');
  finitePos(step, 'step');
  need(efficiency > 0 && efficiency <= 1, 'efficiency must be within (0, 1]');
  need(Cp > 0 && Cp <= (16 / 27) * efficiency + 1e-12, 'Cp must be positive and not exceed 16/27 × efficiency');
  need(cutIn >= 0 && cutOut > cutIn, 'cutOut must exceed cutIn >= 0');
  need(smoothing > 0, 'smoothing must be > 0');
  const A = (Math.PI * D * D) / 4;
  const q = 0.5 * rho * A;
  const ratedSpeed = Math.cbrt(ratedPower / (q * Cp));
  const u = [];
  const p = [];
  const cp = [];
  const ct = [];
  const nSteps = Math.round((cutOut - cutIn) / step);
  for (let i = 0; i <= nSteps; i++) {
    const U = i === nSteps ? cutOut : cutIn + i * step;
    const pc = q * Cp * U * U * U;
    let P;
    if (U <= 0) P = 0;
    else if (Number.isFinite(smoothing)) P = Math.pow(Math.pow(pc, -smoothing) + Math.pow(ratedPower, -smoothing), -1 / smoothing);
    else P = Math.min(pc, ratedPower);
    const cpi = U > 0 ? P / (q * U * U * U) : 0;
    const a = inductionFromCp(cpi / efficiency);
    u.push(U);
    p.push(P);
    cp.push(cpi);
    ct.push(4 * a * (1 - a));
  }
  return { u, p, cp, ct, ratedSpeed, ratedPower, D, rho };
}

/**
 * Power and thrust curve from a rotor-performance callback (for example a blade-element-momentum solver with its
 * own control schedule). The callback receives the wind speed and returns power and thrust.
 * @param {object} o
 * @param {(U:number)=>{power:number, thrust:number}} o.rotor returns shaft/electrical power (W) and thrust (N)
 * @param {number[]} o.speeds wind speeds (m/s), ascending
 * @param {number} o.D rotor diameter (m)
 * @param {number} [o.rho=1.225] air density (kg/m³)
 * @param {number} [o.ratedPower=Infinity] cap applied to the returned power (W)
 * @param {number} [o.efficiency=1] conversion efficiency applied to the callback power (-)
 * @returns {{u:number[], p:number[], cp:number[], ct:number[], thrust:number[], D:number, rho:number}}
 */
export function powerCurveFromRotor({ rotor, speeds, D, rho = 1.225, ratedPower = Infinity, efficiency = 1 }) {
  need(typeof rotor === 'function', 'rotor must be a function of wind speed');
  need(Array.isArray(speeds) && speeds.length > 1, 'speeds must hold at least two values');
  finitePos(D, 'D');
  const A = (Math.PI * D * D) / 4;
  const out = { u: [], p: [], cp: [], ct: [], thrust: [], D, rho };
  for (const U of speeds) {
    const r = rotor(U);
    need(r && Number.isFinite(r.power) && Number.isFinite(r.thrust), 'rotor callback must return finite power and thrust');
    const P = Math.min(ratedPower, Math.max(0, r.power * efficiency));
    out.u.push(U);
    out.p.push(P);
    out.thrust.push(r.thrust);
    out.cp.push(U > 0 ? P / (0.5 * rho * A * U ** 3) : 0);
    out.ct.push(U > 0 ? r.thrust / (0.5 * rho * A * U * U) : 0);
  }
  return out;
}

/**
 * Parse a turbine performance table in the NREL turbine-models CSV layout: a header row naming the columns
 * (wind speed, power, Cp, thrust, Ct — units in square brackets, e.g. "Wind Speed [m/s],Power [kW],Cp [-],
 * Thrust [kN],Ct [-]") followed by numeric rows. Missing Ct is derived from thrust when D and rho are given.
 * @param {string} text CSV text
 * @param {{D?:number, rho?:number}} [opts] rotor diameter (m) and air density (kg/m³) used to fill missing Cp/Ct
 * @returns {{u:number[], p:number[], cp:number[], ct:number[], thrust:number[]}} m/s, W, -, -, N
 */
export function parseTurbineCsv(text, opts = {}) {
  need(typeof text === 'string' && text.trim().length > 0, 'text must be a non-empty string');
  const lines = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && !l.startsWith('#'));
  const split = (l) => l.split(/[,;\t]/).map((c) => c.trim());
  const head = split(lines[0]).map((h) => h.toLowerCase());
  const unitScale = (h) => {
    const m = /\[([^\]]*)\]|\(([^)]*)\)/.exec(h);
    const unit = m ? (m[1] ?? m[2]).trim() : '';
    if (/^mw$/.test(unit) || /^mn$/.test(unit)) return 1e6;
    if (/^kw$/.test(unit) || /^kn$/.test(unit)) return 1e3;
    return 1;
  };
  const find = (re) => head.findIndex((h) => re.test(h.replace(/\[.*?\]|\(.*?\)/g, '').trim()));
  const iU = find(/wind\s*speed|^ws$|^u$|^v$/);
  const iP = find(/^power$|power(?!\s*coef)/);
  const iCp = find(/^cp$|power\s*coef/);
  const iCt = find(/^ct$|thrust\s*coef/);
  const iT = find(/^thrust$/);
  need(iU >= 0 && iP >= 0, 'CSV header must contain wind speed and power columns');
  const pScale = unitScale(head[iP]);
  const tScale = iT >= 0 ? unitScale(head[iT]) : 1;
  const out = { u: [], p: [], cp: [], ct: [], thrust: [] };
  const A = opts.D ? (Math.PI * opts.D * opts.D) / 4 : NaN;
  const rho = opts.rho ?? 1.225;
  for (let r = 1; r < lines.length; r++) {
    const c = split(lines[r]).map(Number);
    if (!Number.isFinite(c[iU]) || !Number.isFinite(c[iP])) continue;
    const U = c[iU];
    const P = c[iP] * pScale;
    const T = iT >= 0 ? c[iT] * tScale : NaN;
    let ct = iCt >= 0 ? c[iCt] : NaN;
    let cp = iCp >= 0 ? c[iCp] : NaN;
    if (!Number.isFinite(ct) && Number.isFinite(T) && Number.isFinite(A) && U > 0) ct = T / (0.5 * rho * A * U * U);
    if (!Number.isFinite(cp) && Number.isFinite(A) && U > 0) cp = P / (0.5 * rho * A * U ** 3);
    out.u.push(U);
    out.p.push(P);
    out.cp.push(cp);
    out.ct.push(ct);
    out.thrust.push(T);
  }
  need(out.u.length > 1, 'CSV must contain at least two numeric rows');
  for (let i = 1; i < out.u.length; i++) need(out.u[i] > out.u[i - 1], 'CSV wind speeds must be strictly ascending');
  return out;
}

/**
 * Site air-density correction of a power curve, IEC 61400-12-1:2017 (air-density normalisation).
 * Pitch-regulated (active power control): the curve is shifted in wind speed, P_site(V) = P_ref(V (ρ/ρ0)^(1/3)).
 * Stall-regulated: the power is scaled, P_site(V) = P_ref(V) ρ/ρ0.
 * @param {{u:number[], p:number[], ct?:number[]}} curve reference curve at density rho0
 * @param {object} o
 * @param {number} o.rho site air density (kg/m³)
 * @param {number} [o.rho0=1.225] reference density of the curve (kg/m³)
 * @param {'pitch'|'stall'} [o.regulation='pitch']
 * @returns {{u:number[], p:number[], ct?:number[], rho:number}} corrected curve
 */
export function airDensityCorrection(curve, { rho, rho0 = 1.225, regulation = 'pitch' }) {
  finitePos(rho, 'rho');
  finitePos(rho0, 'rho0');
  need(regulation === 'pitch' || regulation === 'stall', "regulation must be 'pitch' or 'stall'");
  const out = { ...curve, rho };
  if (regulation === 'stall') {
    out.p = curve.p.map((v) => (v * rho) / rho0);
    out.u = curve.u.slice();
  } else {
    const f = Math.cbrt(rho0 / rho);
    out.u = curve.u.map((v) => v * f);
    out.p = curve.p.slice();
  }
  return out;
}

/**
 * Moist-air density from pressure, temperature and relative humidity (ideal-gas mixture; saturation pressure by
 * the exponential fit used in IEC 61400-12-1: P_w = 0.0000205 exp(0.0631846 T), T in K;
 * gas constants R0 = 287.05 and Rw = 461.5 J/(kg·K)).
 * @param {object} o
 * @param {number} o.pressure barometric pressure (Pa)
 * @param {number} o.temperatureC air temperature (°C)
 * @param {number} [o.relativeHumidity=0] relative humidity (0–1)
 * @returns {number} density (kg/m³)
 */
export function airDensity({ pressure, temperatureC, relativeHumidity = 0 }) {
  finitePos(pressure, 'pressure');
  frac(relativeHumidity, 'relativeHumidity');
  const T = temperatureC + 273.15;
  finitePos(T, 'temperatureC + 273.15');
  const R0 = 287.05;
  const Rw = 461.5;
  const Pw = 0.0000205 * Math.exp(0.0631846 * T);
  return (1 / T) * (pressure / R0 - relativeHumidity * Pw * (1 / R0 - 1 / Rw));
}

/**
 * Weibull distribution of wind speed: pdf, cdf and moments. Mean = A Γ(1 + 1/k); k = 2 is the Rayleigh case.
 * @param {number} A scale parameter (m/s)
 * @param {number} k shape parameter (-)
 * @returns {{pdf:(u:number)=>number, cdf:(u:number)=>number, mean:number, std:number, cubeMean:number,
 *   moment:(n:number)=>number}} mean, std (m/s); cubeMean = E[U³] (m³/s³)
 */
export function weibull(A, k) {
  finitePos(A, 'A');
  finitePos(k, 'k');
  const moment = (n) => Math.pow(A, n) * gamma(1 + n / k);
  const mean = moment(1);
  return {
    pdf: (u) => (u < 0 ? 0 : u === 0 ? (k === 1 ? 1 / A : k < 1 ? Infinity : 0) : (k / A) * Math.pow(u / A, k - 1) * Math.exp(-Math.pow(u / A, k))),
    cdf: (u) => (u <= 0 ? 0 : 1 - Math.exp(-Math.pow(u / A, k))),
    mean,
    std: Math.sqrt(Math.max(0, moment(2) - mean * mean)),
    cubeMean: moment(3),
    moment,
  };
}

/**
 * Weibull scale parameter from the mean wind speed: A = Ū / Γ(1 + 1/k).
 * @param {number} mean mean wind speed (m/s)
 * @param {number} [k=2] shape (-)
 * @returns {number} A (m/s)
 */
export function weibullScaleFromMean(mean, k = 2) {
  finitePos(mean, 'mean');
  finitePos(k, 'k');
  return mean / gamma(1 + 1 / k);
}

/**
 * Vertical extrapolation of wind speed (or of a Weibull scale parameter; the shape k is kept).
 * Power law U2 = U1 (z2/z1)^α; logarithmic law U2 = U1 ln(z2/z0)/ln(z1/z0) (neutral stratification).
 * Reference: IEC 61400-1 §6.3 (α = 0.14 offshore normal wind profile in IEC 61400-3-1; 0.2 onshore).
 * @param {object} o
 * @param {number} o.U speed or Weibull A at height z1 (m/s)
 * @param {number} o.z1 measurement height (m)
 * @param {number} o.z2 target height (m)
 * @param {number} [o.alpha] power-law exponent (-); give either alpha or z0
 * @param {number} [o.z0] roughness length (m)
 * @returns {number} speed at z2 (m/s)
 */
export function extrapolateWind({ U, z1, z2, alpha, z0 }) {
  finitePos(z1, 'z1');
  finitePos(z2, 'z2');
  need(Number.isFinite(U) && U >= 0, 'U must be >= 0');
  if (z0 !== undefined) {
    finitePos(z0, 'z0');
    need(z0 < Math.min(z1, z2), 'z0 must be below both heights');
    return (U * Math.log(z2 / z0)) / Math.log(z1 / z0);
  }
  need(Number.isFinite(alpha), 'either alpha or z0 must be given');
  return U * Math.pow(z2 / z1, alpha);
}

/**
 * Gross annual energy of one turbine under a Weibull wind climate by numerical integration of P(u) f(u)
 * (composite Simpson on a fine uniform grid over the span of the power curve).
 * @param {object} o
 * @param {{u:number[], p:number[]}|((U:number)=>number)} o.curve power curve table or function (W)
 * @param {number} o.A Weibull scale at hub height (m/s)
 * @param {number} [o.k=2] Weibull shape (-)
 * @param {number} [o.uMin] lower integration limit (m/s); default first table speed (0 for a function)
 * @param {number} [o.uMax] upper integration limit (m/s); default last table speed (40 for a function)
 * @param {number} [o.n=4000] number of Simpson intervals (even)
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{aepWh:number, meanPower:number}} annual energy (Wh), mean power (W)
 */
export function aepWeibull({ curve, A, k = 2, uMin, uMax, n = 4000, hoursPerYear = 8760 }) {
  const w = weibull(A, k);
  const isFn = typeof curve === 'function';
  const P = isFn ? curve : (U) => curvePower(curve, U);
  const a = uMin ?? (isFn ? 0 : curve.u[0]);
  const b = uMax ?? (isFn ? 40 : curve.u[curve.u.length - 1]);
  need(b > a, 'uMax must exceed uMin');
  const m = n + (n % 2);
  const h = (b - a) / m;
  let s = 0;
  for (let i = 0; i <= m; i++) {
    const u = a + i * h;
    const f = u === 0 && k < 1 ? 0 : P(u) * w.pdf(u);
    s += f * (i === 0 || i === m ? 1 : i % 2 ? 4 : 2);
  }
  const meanPower = (s * h) / 3;
  return { aepWh: meanPower * hoursPerYear, meanPower };
}

/**
 * Gross annual energy by the method of bins, IEC 61400-12-1:2017 (annual energy production):
 *   AEP = N_h Σ_i [F(V_i) − F(V_(i−1))] (P_(i−1) + P_i)/2, over the points of the tabulated power curve.
 * @param {object} o
 * @param {{u:number[], p:number[]}} o.curve power curve (m/s, W)
 * @param {number} o.A Weibull scale (m/s)
 * @param {number} [o.k=2] Weibull shape (-)
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{aepWh:number, meanPower:number}}
 */
export function aepBins({ curve, A, k = 2, hoursPerYear = 8760 }) {
  const w = weibull(A, k);
  need(curve && curve.u && curve.u.length > 1, 'curve must be a table with at least two points');
  let s = 0;
  for (let i = 1; i < curve.u.length; i++) {
    s += (w.cdf(curve.u[i]) - w.cdf(curve.u[i - 1])) * 0.5 * (curve.p[i] + curve.p[i - 1]);
  }
  return { aepWh: s * hoursPerYear, meanPower: s };
}

/**
 * Normalise a wind rose: sector frequencies are scaled to sum to 1; sectors are validated.
 * @param {{dirDeg:number, frequency:number, A:number, k:number}[]} sectors direction from (deg), frequency (-),
 *   Weibull scale (m/s) and shape (-) at hub height per sector
 * @returns {{dirDeg:number, frequency:number, A:number, k:number}[]} normalised copy
 */
export function normaliseRose(sectors) {
  need(Array.isArray(sectors) && sectors.length > 0, 'sectors must be a non-empty array');
  let sum = 0;
  sectors.forEach((s, i) => {
    need(Number.isFinite(s.dirDeg), `sectors[${i}].dirDeg must be finite`);
    need(s.frequency >= 0, `sectors[${i}].frequency must be >= 0`);
    finitePos(s.A, `sectors[${i}].A`);
    finitePos(s.k, `sectors[${i}].k`);
    sum += s.frequency;
  });
  finitePos(sum, 'sum of sector frequencies');
  return sectors.map((s) => ({ ...s, frequency: s.frequency / sum }));
}

/**
 * Wind-farm annual energy: loop over direction sectors and wind-speed bins, solving the farm flow (wake.js) in
 * each state and weighting by the Weibull probability of the bin, P(bin) = F(u + Δ/2) − F(u − Δ/2).
 * @param {object} o
 * @param {{x:number, y:number, D:number, hubHeight:number, curve?:{u:number[],p:number[],ct:number[]},
 *   power?:Function|object, ct?:Function|object}[]} o.turbines positions (m) and either a `curve` table or
 *   `power`/`ct` as accepted by `farmFlow`
 * @param {{dirDeg:number, frequency:number, A:number, k:number}[]} o.rose wind rose at hub height
 * @param {number} [o.uMin=0.5] centre of the first speed bin (m/s)
 * @param {number} [o.uMax=30] upper limit of bin centres (m/s)
 * @param {number} [o.du=1] bin width (m/s)
 * @param {number} [o.subSectors=1] wind directions evaluated per sector (uniformly spread over the sector width)
 * @param {object} [o.wake] options passed to `farmFlow` (model, TI, superposition, …)
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{grossWh:number, aepWh:number, arrayEfficiency:number, wakeLoss:number, perTurbine:{grossWh:number,
 *   aepWh:number, efficiency:number}[], perSector:{dirDeg:number, frequency:number, grossWh:number,
 *   aepWh:number, efficiency:number}[], ratedPower:number, capacityFactor:number}} energies (Wh/yr),
 *   efficiencies (-); ratedPower = Σ max of each power curve (W)
 */
export function farmAep({ turbines, rose, uMin = 0.5, uMax = 30, du = 1, subSectors = 1, wake = {}, hoursPerYear = 8760 }) {
  need(Array.isArray(turbines) && turbines.length > 0, 'turbines must be a non-empty array');
  finitePos(du, 'du');
  need(uMax > uMin, 'uMax must exceed uMin');
  need(Number.isInteger(subSectors) && subSectors >= 1, 'subSectors must be an integer >= 1');
  const sectors = normaliseRose(rose);
  const width = 360 / sectors.length;
  const ts = turbines.map((t, i) => {
    if (t.curve) {
      need(t.curve.ct, `turbines[${i}].curve needs a ct column`);
      return { ...t, power: { u: t.curve.u, v: t.curve.p }, ct: { u: t.curve.u, v: t.curve.ct } };
    }
    need(t.power && t.ct, `turbines[${i}] needs either curve or power and ct`);
    return t;
  });
  const n = ts.length;
  const grossT = new Float64Array(n);
  const netT = new Float64Array(n);
  const perSector = [];
  const nBins = Math.floor((uMax - uMin) / du + 1e-9) + 1;
  for (const sec of sectors) {
    const w = weibull(sec.A, sec.k);
    let gS = 0;
    let nS = 0;
    for (let b = 0; b < nBins; b++) {
      const U = uMin + b * du;
      const prob = (w.cdf(U + du / 2) - w.cdf(Math.max(0, U - du / 2))) * sec.frequency;
      if (prob <= 0) continue;
      for (let d = 0; d < subSectors; d++) {
        const wd = sec.dirDeg + width * ((d + 0.5) / subSectors - 0.5);
        const res = farmFlow({ ...wake, turbines: ts, U, wdDeg: wd });
        if (res.grossPower <= 0) break;
        const wgt = (prob / subSectors) * hoursPerYear;
        for (let i = 0; i < n; i++) {
          netT[i] += res.power[i] * wgt;
        }
        // gross per turbine: its own curve at the undisturbed speed
        for (let i = 0; i < n; i++) {
          const g = typeof ts[i].power === 'function' ? ts[i].power(U) : interpTable(ts[i].power.u, ts[i].power.v, U);
          grossT[i] += g * wgt;
          gS += g * wgt;
        }
        nS += res.totalPower * wgt;
      }
    }
    perSector.push({ dirDeg: sec.dirDeg, frequency: sec.frequency, grossWh: gS, aepWh: nS, efficiency: gS > 0 ? nS / gS : 1 });
  }
  let gross = 0;
  let net = 0;
  let rated = 0;
  const perTurbine = [];
  for (let i = 0; i < n; i++) {
    gross += grossT[i];
    net += netT[i];
    perTurbine.push({ grossWh: grossT[i], aepWh: netT[i], efficiency: grossT[i] > 0 ? netT[i] / grossT[i] : 1 });
    const pw = ts[i].power;
    if (typeof pw !== 'function') rated += Math.max(...pw.v);
    else if (Number.isFinite(ts[i].ratedPower)) rated += ts[i].ratedPower;
    else rated = NaN;
  }
  const eff = gross > 0 ? net / gross : 1;
  return {
    grossWh: gross,
    aepWh: net,
    arrayEfficiency: eff,
    wakeLoss: 1 - eff,
    perTurbine,
    perSector,
    ratedPower: rated,
    capacityFactor: rated > 0 ? net / (rated * hoursPerYear) : NaN,
  };
}

/**
 * Loss cascade from gross to net energy. Losses are independent fractions applied multiplicatively in the order
 * given: E_net = E_gross Π (1 − l_i). Typical categories (DNV-style energy-loss taxonomy / IEC 61400-15 draft):
 * wake, blockage, electrical, availabilityTurbine, availabilityBop, availabilityGrid, icing, degradation,
 * hysteresis (high-wind), curtailment.
 * @param {object} o
 * @param {number} o.grossWh gross annual energy (Wh)
 * @param {Object<string, number>|{name:string, loss:number}[]} o.losses loss fractions (-), ordered
 * @param {number} [o.ratedPower] installed capacity (W) for capacity factor and full-load hours
 * @param {number} [o.hoursPerYear=8760]
 * @returns {{netWh:number, totalLoss:number, steps:{name:string, loss:number, energyWh:number,
 *   lostWh:number}[], capacityFactor:number, fullLoadHours:number}} full-load hours in h
 */
export function lossCascade({ grossWh, losses, ratedPower, hoursPerYear = 8760 }) {
  need(Number.isFinite(grossWh) && grossWh >= 0, 'grossWh must be >= 0');
  const list = Array.isArray(losses) ? losses : Object.entries(losses ?? {}).map(([name, loss]) => ({ name, loss }));
  let e = grossWh;
  const steps = [];
  for (const { name, loss } of list) {
    frac(loss, `loss '${name}'`);
    const lost = e * loss;
    e -= lost;
    steps.push({ name, loss, energyWh: e, lostWh: lost });
  }
  const hasRated = Number.isFinite(ratedPower) && ratedPower > 0;
  return {
    netWh: e,
    totalLoss: grossWh > 0 ? 1 - e / grossWh : 0,
    steps,
    capacityFactor: hasRated ? e / (ratedPower * hoursPerYear) : NaN,
    fullLoadHours: hasRated ? e / ratedPower : NaN,
  };
}

/**
 * Standard normal quantile z_p with P[Z ≤ z_p] = p.
 * @param {number} p probability (0–1, exclusive)
 * @returns {number}
 */
export function normalQuantile(p) {
  need(p > 0 && p < 1, 'p must be within (0, 1)');
  return Math.SQRT2 * erfinv(2 * p - 1);
}

/**
 * Combine independent uncertainty components into the standard uncertainty of annual energy and derive
 * exceedance levels. Wind-speed uncertainties (fractions of mean wind speed) are converted to energy with the
 * sensitivity s = (dE/E)/(dU/U); energy components are fractions of energy. Inter-annual variability of energy
 * (one-year standard deviation, fraction) is reduced by √n for an n-year averaging period.
 *   σ_E = √( Σ (s σ_U,i)² + Σ σ_E,j² + σ_IAV²/n ).
 * Exceedance convention: Pxx is the energy exceeded with probability xx %, assuming a normal distribution:
 *   Pxx = P50 (1 − z_xx σ_E),  z_75 = 0.6745, z_90 = 1.2816, z_99 = 2.3263.
 * @param {object} o
 * @param {number} o.p50Wh central (P50) net annual energy (Wh)
 * @param {Object<string, number>} [o.windSpeed] wind-speed uncertainty components (fractions of wind speed):
 *   e.g. measurement, longTermCorrection, verticalExtrapolation, horizontalExtrapolation
 * @param {number} [o.sensitivity] energy-to-wind-speed sensitivity s (-), typically 1.5–2.5 (required when
 *   `windSpeed` is given)
 * @param {Object<string, number>} [o.energy] energy uncertainty components (fractions of energy): e.g.
 *   wakeModel, powerCurve, losses
 * @param {number} [o.interAnnual=0] one-year inter-annual variability of energy (fraction)
 * @param {number} [o.years=1] averaging period n (years)
 * @param {number[]} [o.levels=[75, 90, 99]] exceedance levels (%)
 * @returns {{sigma:number, sigmaWh:number, components:{name:string, sigma:number}[], p50Wh:number,
 *   p75Wh:number, p90Wh:number, p99Wh:number, exceedance:Object<string, number>}} sigma as a fraction of P50
 */
export function aepUncertainty({ p50Wh, windSpeed = {}, sensitivity, energy = {}, interAnnual = 0, years = 1, levels = [75, 90, 99] }) {
  need(Number.isFinite(p50Wh), 'p50Wh must be finite');
  finitePos(years, 'years');
  need(interAnnual >= 0, 'interAnnual must be >= 0');
  const components = [];
  const wsNames = Object.keys(windSpeed);
  if (wsNames.length) finitePos(sensitivity, 'sensitivity');
  for (const name of wsNames) {
    need(windSpeed[name] >= 0, `windSpeed.${name} must be >= 0`);
    components.push({ name, sigma: sensitivity * windSpeed[name] });
  }
  for (const name of Object.keys(energy)) {
    need(energy[name] >= 0, `energy.${name} must be >= 0`);
    components.push({ name, sigma: energy[name] });
  }
  if (interAnnual > 0) components.push({ name: 'interAnnual', sigma: interAnnual / Math.sqrt(years) });
  let v = 0;
  for (const c of components) v += c.sigma * c.sigma;
  const sigma = Math.sqrt(v);
  const exceedance = {};
  for (const L of levels) {
    need(L > 0 && L < 100, 'levels must be within (0, 100)');
    exceedance[`p${L}`] = p50Wh * (1 - normalQuantile(L / 100) * sigma);
  }
  const at = (L) => p50Wh * (1 - normalQuantile(L / 100) * sigma);
  return { sigma, sigmaWh: sigma * p50Wh, components, p50Wh, p75Wh: at(75), p90Wh: at(90), p99Wh: at(99), exceedance };
}

function moments2(x, y) {
  const n = x.length;
  let mx = 0;
  let my = 0;
  for (let i = 0; i < n; i++) {
    mx += x[i];
    my += y[i];
  }
  mx /= n;
  my /= n;
  let sxx = 0;
  let syy = 0;
  let sxy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
    sxy += (x[i] - mx) * (y[i] - my);
  }
  return { mx, my, sxx, syy, sxy, n };
}

/**
 * Measure–correlate–predict on concurrent site and reference wind-speed series.
 *  - 'linear': ordinary least squares y = a + b x.
 *  - 'varianceRatio': y = μ_y + (σ_y/σ_x)(x − μ_x) (Rogers, Rogers & Manwell 2005), which preserves the variance
 *    of the predicted series.
 * The fitted relation is applied to the long-term reference series to give the long-term site estimate.
 * @param {object} o
 * @param {number[]} o.site concurrent site wind speeds (m/s)
 * @param {number[]} o.reference concurrent reference wind speeds (m/s)
 * @param {number[]} [o.longTermReference] long-term reference series (m/s)
 * @param {'linear'|'varianceRatio'} [o.method='linear']
 * @returns {{slope:number, intercept:number, r:number, r2:number, predict:(x:number)=>number,
 *   longTermMean:number, concurrentMean:number, correctionFactor:number, predicted:Float64Array|null}}
 *   longTermMean (m/s) is NaN without a long-term series; correctionFactor = longTermMean / concurrent site mean
 */
export function mcp({ site, reference, longTermReference, method = 'linear' }) {
  need(site && reference && site.length === reference.length && site.length >= 3,
    'site and reference must be concurrent series of equal length >= 3');
  need(method === 'linear' || method === 'varianceRatio', "method must be 'linear' or 'varianceRatio'");
  const m = moments2(reference, site);
  need(m.sxx > 0, 'reference series has zero variance');
  const r = m.sxy / Math.sqrt(m.sxx * m.syy || 1);
  const slope = method === 'linear' ? m.sxy / m.sxx : Math.sqrt(m.syy / m.sxx) * (r < 0 ? -1 : 1);
  const intercept = m.my - slope * m.mx;
  const predict = (x) => intercept + slope * x;
  let longTermMean = NaN;
  let predicted = null;
  if (longTermReference && longTermReference.length) {
    predicted = new Float64Array(longTermReference.length);
    let s = 0;
    for (let i = 0; i < predicted.length; i++) {
      predicted[i] = Math.max(0, predict(longTermReference[i]));
      s += predicted[i];
    }
    longTermMean = s / predicted.length;
  }
  return {
    slope,
    intercept,
    r,
    r2: r * r,
    predict,
    longTermMean,
    concurrentMean: m.my,
    correctionFactor: longTermMean / m.my,
    predicted,
  };
}

/**
 * Production from a wind-speed time series at hub height, with aggregation.
 * @param {object} o
 * @param {number[]|Float64Array} o.speeds hub-height wind speed per step (m/s)
 * @param {{u:number[], p:number[]}|((U:number)=>number)} o.curve power curve (W)
 * @param {number} [o.dtHours=1] time step (h)
 * @param {number[]} [o.times] time stamps (ms since epoch, UTC) for monthly/daily/diurnal aggregation
 * @param {number} [o.lossFactor=1] net-to-gross factor applied to every step (-)
 * @param {number} [o.ratedPower] rated power (W) for the capacity factor; default max of the curve table
 * @returns {{power:Float64Array, energyWh:number, meanPower:number, capacityFactor:number,
 *   monthly:{month:string, energyWh:number, hours:number}[], daily:{day:string, energyWh:number}[],
 *   diurnal:number[]}} power (W) per step; diurnal = mean power (W) by UTC hour (24 values, NaN without times)
 */
export function timeSeriesProduction({ speeds, curve, dtHours = 1, times, lossFactor = 1, ratedPower }) {
  need(speeds && speeds.length > 0, 'speeds must be a non-empty series');
  finitePos(dtHours, 'dtHours');
  need(!times || times.length === speeds.length, 'times must have the same length as speeds');
  const isFn = typeof curve === 'function';
  const n = speeds.length;
  const power = new Float64Array(n);
  let e = 0;
  const monthly = new Map();
  const daily = new Map();
  const dSum = new Float64Array(24);
  const dCnt = new Float64Array(24);
  for (let i = 0; i < n; i++) {
    const p = (isFn ? curve(speeds[i]) : curvePower(curve, speeds[i])) * lossFactor;
    power[i] = p;
    e += p * dtHours;
    if (times) {
      const d = new Date(times[i]);
      const iso = d.toISOString();
      const mk = iso.slice(0, 7);
      const dk = iso.slice(0, 10);
      const mo = monthly.get(mk) ?? { month: mk, energyWh: 0, hours: 0 };
      mo.energyWh += p * dtHours;
      mo.hours += dtHours;
      monthly.set(mk, mo);
      daily.set(dk, (daily.get(dk) ?? 0) + p * dtHours);
      dSum[d.getUTCHours()] += p;
      dCnt[d.getUTCHours()] += 1;
    }
  }
  const rated = ratedPower ?? (isFn ? NaN : Math.max(...curve.p));
  const meanPower = e / (n * dtHours);
  return {
    power,
    energyWh: e,
    meanPower,
    capacityFactor: rated > 0 ? meanPower / rated : NaN,
    monthly: [...monthly.values()],
    daily: [...daily.entries()].map(([day, energyWh]) => ({ day, energyWh })),
    diurnal: Array.from(dSum, (s, h) => (dCnt[h] ? s / dCnt[h] : NaN)),
  };
}

/**
 * Peak-period capacity factor: mean generation during the highest-demand fraction of time steps divided by the
 * rated power — a simple proxy for capacity value (capacity credit proper needs a loss-of-load calculation;
 * this is the "top demand hours" approximation, e.g. Milligan & Porter 2008).
 * @param {object} o
 * @param {number[]|Float64Array} o.generation generation per step (W)
 * @param {number[]|Float64Array} o.demand system demand per step (any consistent unit)
 * @param {number} o.ratedPower installed capacity (W)
 * @param {number} [o.topFraction=0.1] fraction of steps with the highest demand (-)
 * @returns {{capacityValue:number, steps:number}} capacity value (-), number of steps used
 */
export function capacityValue({ generation, demand, ratedPower, topFraction = 0.1 }) {
  need(generation && demand && generation.length === demand.length && generation.length > 0,
    'generation and demand must be series of equal non-zero length');
  finitePos(ratedPower, 'ratedPower');
  need(topFraction > 0 && topFraction <= 1, 'topFraction must be within (0, 1]');
  const idx = Array.from(demand, (_, i) => i).sort((a, b) => demand[b] - demand[a]);
  const m = Math.max(1, Math.round(topFraction * idx.length));
  let s = 0;
  for (let i = 0; i < m; i++) s += generation[idx[i]];
  return { capacityValue: s / m / ratedPower, steps: m };
}

/**
 * Turbulence-intensity summary from 10-minute statistics, binned by wind speed.
 * Representative TI per bin is the 90 % quantile assuming normality, mean + 1.28 σ (IEC 61400-1, normal turbulence model).
 * @param {object} o
 * @param {number[]} o.mean 10-min mean wind speed (m/s)
 * @param {number[]} o.std 10-min standard deviation (m/s)
 * @param {number} [o.binWidth=1] bin width (m/s)
 * @param {number} [o.minSpeed=3] lowest speed included (m/s)
 * @returns {{meanTI:number, bins:{u:number, n:number, meanTI:number, stdTI:number, representativeTI:number}[],
 *   ti15:number}} overall mean TI (-); per-bin statistics; representative TI of the 15 m/s bin (NaN if empty)
 */
export function turbulenceSummary({ mean, std, binWidth = 1, minSpeed = 3 }) {
  need(mean && std && mean.length === std.length && mean.length > 0, 'mean and std must be series of equal non-zero length');
  finitePos(binWidth, 'binWidth');
  const acc = new Map();
  let s = 0;
  let c = 0;
  for (let i = 0; i < mean.length; i++) {
    if (!(mean[i] >= minSpeed)) continue;
    const ti = std[i] / mean[i];
    const b = Math.round(mean[i] / binWidth);
    const a = acc.get(b) ?? { n: 0, s: 0, s2: 0 };
    a.n++;
    a.s += ti;
    a.s2 += ti * ti;
    acc.set(b, a);
    s += ti;
    c++;
  }
  const bins = [...acc.entries()].sort((p, q) => p[0] - q[0]).map(([b, a]) => {
    const m = a.s / a.n;
    const sd = a.n > 1 ? Math.sqrt(Math.max(0, (a.s2 - a.n * m * m) / (a.n - 1))) : 0;
    return { u: b * binWidth, n: a.n, meanTI: m, stdTI: sd, representativeTI: m + 1.28 * sd };
  });
  const b15 = bins.find((b) => Math.abs(b.u - 15) < binWidth / 2);
  return { meanTI: c ? s / c : NaN, bins, ti15: b15 ? b15.representativeTI : NaN };
}

/**
 * Wind-shear summary from concurrent speeds at two heights: power-law exponent α = ln(U2/U1)/ln(z2/z1) from
 * the mean speeds and as the mean of per-sample values, and the equivalent neutral roughness length from the
 * log law, ln z0 = (U2 ln z1 − U1 ln z2)/(U2 − U1).
 * @param {object} o
 * @param {number[]} o.u1 speeds at the lower height (m/s)
 * @param {number[]} o.u2 speeds at the upper height (m/s)
 * @param {number} o.z1 lower height (m)
 * @param {number} o.z2 upper height (m)
 * @param {number} [o.minSpeed=3] samples with u1 below this are ignored (m/s)
 * @returns {{alpha:number, alphaMeanOfSamples:number, z0:number, n:number}}
 */
export function shearSummary({ u1, u2, z1, z2, minSpeed = 3 }) {
  need(u1 && u2 && u1.length === u2.length && u1.length > 0, 'u1 and u2 must be series of equal non-zero length');
  finitePos(z1, 'z1');
  need(z2 > z1, 'z2 must exceed z1');
  const L = Math.log(z2 / z1);
  let s1 = 0;
  let s2 = 0;
  let sa = 0;
  let n = 0;
  for (let i = 0; i < u1.length; i++) {
    if (!(u1[i] >= minSpeed) || !(u2[i] > 0)) continue;
    s1 += u1[i];
    s2 += u2[i];
    sa += Math.log(u2[i] / u1[i]) / L;
    n++;
  }
  need(n > 0, 'no samples above minSpeed');
  const m1 = s1 / n;
  const m2 = s2 / n;
  const z0 = m2 > m1 ? Math.exp((m2 * Math.log(z1) - m1 * Math.log(z2)) / (m2 - m1)) : NaN;
  return { alpha: Math.log(m2 / m1) / L, alphaMeanOfSamples: sa / n, z0, n };
}

/**
 * Extreme wind speed for a return period by a Gumbel fit to annual maxima:
 *   V_T = μ − β ln(−ln(1 − 1/T)).
 * Reference: IEC 61400-1 (V_ref is the 50-year 10-min mean); Gumbel (1958). The optional gust factor converts
 * the averaging period (e.g. 10-min mean → 3-s gust; IEC 61400-1 uses V_e50 = 1.4 V_ref at hub height).
 * @param {object} o
 * @param {number[]} o.annualMaxima annual maximum wind speeds (m/s), at least 5 values recommended ≥ 10
 * @param {number} [o.returnPeriod=50] return period T (years), > 1
 * @param {number} [o.gustFactor=1] multiplier applied to the result (-)
 * @param {'mle'|'mom'|'lmom'} [o.method='mom'] fitting method
 * @returns {{value:number, loc:number, scale:number, returnPeriod:number}} extreme speed (m/s), Gumbel μ and β
 */
export function extremeWind({ annualMaxima, returnPeriod = 50, gustFactor = 1, method = 'mom' }) {
  need(annualMaxima && annualMaxima.length >= 5, 'annualMaxima must contain at least 5 values');
  need(returnPeriod > 1, 'returnPeriod must exceed 1 year');
  finitePos(gustFactor, 'gustFactor');
  const fit = fitGumbel(Array.from(annualMaxima), { method });
  const value = (fit.loc - fit.scale * Math.log(-Math.log(1 - 1 / returnPeriod))) * gustFactor;
  return { value, loc: fit.loc, scale: fit.scale, returnPeriod };
}

/**
 * Generation profile ready for revenue/dispatch studies: energy per step and, with a price series, revenue.
 * @param {object} o
 * @param {number[]|Float64Array} o.power generation per step (W)
 * @param {number} [o.dtHours=1] step (h)
 * @param {number[]} [o.prices] price per step (currency per MWh)
 * @param {number} [o.exportLimit=Infinity] export limit (W); the excess is reported as curtailed
 * @returns {{energyMWh:Float64Array, curtailedMWh:number, totalMWh:number, revenue:number,
 *   capturePrice:number}} revenue (currency), capture price (currency/MWh; NaN without prices)
 */
export function generationProfile({ power, dtHours = 1, prices, exportLimit = Infinity }) {
  need(power && power.length > 0, 'power must be a non-empty series');
  finitePos(dtHours, 'dtHours');
  need(!prices || prices.length === power.length, 'prices must have the same length as power');
  const energyMWh = new Float64Array(power.length);
  let total = 0;
  let curtailed = 0;
  let revenue = 0;
  for (let i = 0; i < power.length; i++) {
    const p = Math.min(power[i], exportLimit);
    curtailed += ((power[i] - p) * dtHours) / 1e6;
    energyMWh[i] = (p * dtHours) / 1e6;
    total += energyMWh[i];
    if (prices) revenue += energyMWh[i] * prices[i];
  }
  return { energyMWh, curtailedMWh: curtailed, totalMWh: total, revenue: prices ? revenue : NaN, capturePrice: prices && total > 0 ? revenue / total : NaN };
}
