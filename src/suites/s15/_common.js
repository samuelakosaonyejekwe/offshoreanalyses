// Shared helpers of the offshore installation and marine operations suite: unit constants, result builders,
// the recalled DNV-ST-N001 lifting factors (all overridable through module inputs), a seeded synthetic metocean
// series for weather-window work, and small structural checks used by several modules. Engineering units are
// converted to SI inside each module (docs/CONTRACT.md §2.1).
import { makeRng } from '../../engine/math/rng.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const RHO_STEEL = 7850;
// Screening emission factors (order of magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;     // t CO2e per t steel
export const CO2_CONCRETE = 0.3;  // t CO2e per m³ concrete
export const CO2_MGO = 3.206;     // t CO2 per t marine gas oil burnt

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible operation.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / Math.max(1, n - 1));
/** Decimate equally long arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 400) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
/** Piecewise-linear interpolation on ascending xs, clamped at the ends. */
export function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1; while (xs[i] < x) i++;
  return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
}

/* ------------------------------------------------------------------------------------------- result builders */

export const kpi = (key, label, value, unit, explain, limit, sense) => {
  const k = { key, label, value: finite(value, label), unit, explain };
  if (limit != null) { k.limit = limit; k.sense = sense || 'max'; }
  return k;
};
export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : (Math.abs(demand) > 0 ? 99 : 0);
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
/** Check written as required factor of safety over achieved factor of safety. */
export const fsCheck = (name, fs, fsReq, ref, explain) => (
  { name, demand: round(fsReq, 4), capacity: round(cap(fs), 4), unit: 'FS (required / achieved)', uc: fs > 0 ? Math.min(fsReq / fs, 99) : 99, ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc) => (uc > 1 ? 'fail' : uc > 0.9 ? 'warn' : 'ok');
export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const line = (title, xl, xu, yl, yu, series, explain, type = 'line') => (
  { type, title, x: { label: xl, unit: xu }, y: { label: yl, unit: yu }, series, explain });
export const bar = (title, categories, yl, yu, series, explain) => (
  { type: 'bar', title, categories, y: { label: yl, unit: yu }, series, explain });
export const tbl = (title, columns, rows, explain) => ({ title, columns: columns.map(([key, label, unit]) => (unit ? { key, label, unit } : { key, label })), rows, explain });

/* ------------------------------------------------------------------------------------------- standing text */

export const CONFIRM = 'Default recalled from the standard — confirm against the governing project edition.';
export const N001 = 'DNV-ST-N001 (2021, amended 2023)';
export const N103 = 'DNV-RP-N103 (2017, amended 2021)';
export const LIMIT_N001 = 'The DNV-ST-N001 factors used here (dynamic amplification by weight class, skew, weight and centre-of-gravity factors, consequence factors, alpha-factors) are recalled defaults exposed as criteria inputs: confirm each against the edition governing the project.';
export const LIMIT_N103 = 'The DNV-RP-N103 coefficients (slamming, drag, added mass, characteristic amplitudes) are recalled defaults exposed as inputs: confirm against the governing edition and the object geometry.';
export const LIMIT_ROM = 'Reduced-order planar rigid-body model, not a coupled multibody hydrodynamic simulation: radiation/diffraction, wave excitation in six degrees of freedom and structural flexibility are outside the module.';

/* ------------------------------------------------------------------------------------------- sustainability */

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (screening factor — confirm with the mill's product declaration).`;
export const vesselCo2Text = (days, burn = 40) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * burn)} t marine gas oil ≈ ${fmt(days * burn * CO2_MGO)} t CO₂ (screening: ${burn} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

/* ------------------------------------------------------------------------------------------- lifting factors */

/**
 * Dynamic amplification factor for lifts in air by weight class, as recalled from DNV-ST-N001 Sec. 16
 * (Table "DAF in air"). W is the static hook load in tonnes. Offshore: 1 + 0.25·√(100/W) up to 100 t, then
 * 1.25 / 1.20 / 1.15 / 1.10 for ≤ 300 / ≤ 1000 / ≤ 2500 / above. Inshore: 1.07 + 0.05·√(100/W), then
 * 1.12 / 1.10 / 1.08 / 1.05. Onshore: 1.10, then 1.05 / 1.05 / 1.03 / 1.03. Always overridable by the caller.
 */
export function dafN001(Wt, where = 'offshore') {
  const W = Math.max(Wt, 3);
  const cls = W <= 100 ? 0 : W <= 300 ? 1 : W <= 1000 ? 2 : W <= 2500 ? 3 : 4;
  if (where === 'onshore') return [1.10, 1.05, 1.05, 1.03, 1.03][cls];
  if (where === 'inshore') return cls === 0 ? 1.07 + 0.05 * Math.sqrt(100 / W) : [0, 1.12, 1.10, 1.08, 1.05][cls];
  return cls === 0 ? 1 + 0.25 * Math.sqrt(100 / W) : [0, 1.25, 1.20, 1.15, 1.10][cls];
}
export const LIFT_SITES = [{ value: 'offshore', label: 'Offshore' }, { value: 'inshore', label: 'Inshore (sheltered)' }, { value: 'onshore', label: 'Onshore / quay' }];

/** Common lifted-weight inputs: weight, contingency, rigging weight, location, DAF override. */
export function liftWeightInputs({ W = 1200, rig = 60, group = 'Lifted object' } = {}) {
  return [
    { key: 'W', label: 'Net (weighed or calculated) weight', unit: 't', type: 'number', default: W, min: 1, max: 60000, group, role: 'load', bind: 'install.liftWeight', uncertain: { dist: 'normal', cov: 0.03 },
      help: 'Weight report value without contingency. Bound to the lifted weight on the project bus when an upstream module publishes it.' },
    { key: 'fWeight', label: 'Weight contingency factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group, role: 'criteria',
      help: `DNV-ST-N001: about 1.03 for a weighed object, 1.05–1.10 for a calculated weight, more for an early estimate. ${CONFIRM}` },
    { key: 'Wrig', label: 'Rigging weight (slings, shackles, spreader)', unit: 't', type: 'number', default: rig, min: 0, max: 3000, group, role: 'load', help: 'Everything hanging from the hook that is not the object; typically 3–7 % of the lifted weight.' },
    { key: 'where', label: 'Lift location', type: 'select', default: 'offshore', options: LIFT_SITES, group: 'Operation', role: 'operation', help: 'Selects the dynamic-amplification column of the recalled DNV-ST-N001 table.' },
    { key: 'dafUser', label: 'DAF override (0 = table by weight class)', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Criteria', role: 'criteria',
      help: `Enter a DAF from a lift dynamic analysis or from the governing edition's table to replace the recalled table value. ${CONFIRM}` },
  ];
}

/* ------------------------------------------------------------------------------------------- small checks */

/** EN 1993-1-1 §6.3.1.2 column reduction factor χ for relative slenderness λ̄ and imperfection factor α. */
export function chiEC3(lambdaBar, alpha = 0.49) {
  if (lambdaBar <= 0.2) return 1;
  const phi = 0.5 * (1 + alpha * (lambdaBar - 0.2) + lambdaBar * lambdaBar);
  return Math.min(1, 1 / (phi + Math.sqrt(phi * phi - lambdaBar * lambdaBar)));
}
/** Thin-walled tube section properties, SI. */
export function tube(D, t) {
  const d = D - 2 * t;
  const A = Math.PI / 4 * (D * D - d * d), I = Math.PI / 64 * (D ** 4 - d ** 4);
  return { A, I, W: I / (D / 2), r: Math.sqrt(I / A), Wp: (D ** 3 - d ** 3) / 6, mass: A * RHO_STEEL };
}
/**
 * Padeye plate at the pin hole by plain section statics: pin bearing F/(dpin·t), tear-out on two shear planes
 * F/(2·(R − dh/2)·t) and net tension F/((2R − dh)·t). Returned stresses in Pa; the allowables are the caller's.
 */
export function padeyeStress({ F, t, R, dHole, dPin }) {
  return { bearing: F / (dPin * t), tearOut: F / (2 * (R - dHole / 2) * t), netTension: F / ((2 * R - dHole) * t) };
}

/* ------------------------------------------------------------------------------------------- sea state */

export function seaInputs({ Hs = 1.5, Tp = 7, group = 'Sea state', role = 'environment', hsMax = 12 } = {}) {
  return [
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: Hs, min: 0, max: hsMax, group, role, design: true, uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Design (limiting) sea state of the operation, not the forecast limit: the forecast limit is the design value times the alpha-factor.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: Tp, min: 2, max: 25, group, role, uncertain: { dist: 'normal', cov: 0.1 }, help: 'Check the whole range of periods that can accompany the limiting wave height; resonance often governs.' },
  ];
}

/**
 * Seeded synthetic metocean series for weather-window work when no measured or hindcast series is supplied:
 * ln Hs follows a first-order autoregressive process (persistence time tauHours) about a seasonal mean with a
 * winter peak; wind follows Hs with independent scatter, Tp = tpFactor·√Hs with scatter. monthlyMean (12 values
 * of mean Hs, m) calibrates the seasonal cycle when given. This is a stand-in with realistic persistence, not a
 * hindcast: results obtained from it are indicative only.
 */
export function syntheticWeather({ years = 3, dtHours = 3, hsMean = 1.8, seasonality = 0.35, sigmaLn = 0.5, tauHours = 36, tpFactor = 4.4, windFactor = 5.5, monthlyMean = null, seed = 7 }) {
  const perYear = Math.round(8760 / dtHours), n = perYear * years;
  const rng = makeRng(seed);
  const phi = Math.exp(-dtHours / tauHours), s = Math.sqrt(1 - phi * phi);
  const hs = new Float64Array(n), tp = new Float64Array(n), wind = new Float64Array(n), timeMs = new Float64Array(n);
  const t0 = Date.UTC(2001, 0, 1);
  let z = rng.normal(), zw = rng.normal();
  for (let i = 0; i < n; i++) {
    z = phi * z + s * rng.normal(); zw = phi * zw + s * rng.normal();
    const frac = (i % perYear) / perYear; // 0 = 1 January
    let m;
    if (monthlyMean) { const k = frac * 12 - 0.5, a = Math.floor(k), f = k - a; m = monthlyMean[(a + 12) % 12] * (1 - f) + monthlyMean[(a + 13) % 12] * f; }
    else m = hsMean * (1 + seasonality * Math.cos(2 * Math.PI * (frac - 0.04)));
    // mean of a lognormal = exp(mu + σ²/2): keep the target mean
    hs[i] = Math.exp(Math.log(Math.max(m, 0.05)) - sigmaLn * sigmaLn / 2 + sigmaLn * z);
    tp[i] = Math.max(2.5, tpFactor * Math.sqrt(hs[i]) * (1 + 0.12 * rng.normal()));
    wind[i] = Math.max(0, windFactor * Math.sqrt(hs[i]) * (1 + 0.25 * zw));
    timeMs[i] = t0 + i * dtHours * 3.6e6;
  }
  return { hs, tp, wind, timeMs, dtHours };
}

/** JONSWAP spectral moments helper is in the engines; this is the Rayleigh most probable maximum factor √(2·ln N). */
export const mpmFactor = (durationS, Tz) => Math.sqrt(2 * Math.log(Math.max(durationS / Tz, 2)));
