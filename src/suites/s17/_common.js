// Shared helpers of the economics suite: validation, rounding, advice records, table coercion, the wording
// required on every price-like input, percentile conventions and the screening emission factors.
// Money is carried in plain project-currency units (shown as USD) so that bus values bind without scaling;
// headline results are reported in millions.
import { irr as irrRoots } from '../../engine/econ/cashflow.js';
import { dist } from '../../engine/math/stats.js';
import { makeRng } from '../../engine/math/rng.js';

export const CUR = 'USD';
export const MU = `M ${CUR}`;
export const PLACEHOLDER = 'Indicative placeholder — replace with current quotations or indices.';
export const CUR_NOTE = `Amounts are in the project currency, labelled ${CUR}; no default is a market quotation.`;

// Screening emission factors (cradle-to-gate order of magnitude, not project EPD values).
export const CO2_STEEL = 1.9;   // t CO2e per t steel
export const CO2_CONC = 0.3;    // t CO2e per m³ concrete
export const CO2_MGO = 3.206;   // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 25; // main installation vessel fuel burn, t/day (screening)
export const SCREEN = `screening factors: ${CO2_STEEL} t CO₂e/t steel, ${CO2_CONC} t CO₂e/m³ concrete, ${CO2_MGO} t CO₂ per t marine gas oil`;

export const PCT_CUM = 'Percentile convention: cumulative (non-exceedance) — Pxx is the value not exceeded with xx % probability, so P10 is the low value and P90 the high value.';
export const PCT_COST = `${PCT_CUM} For a cost, P90 is the high-cost (conservative) case.`;
export const PCT_VALUE = `${PCT_CUM} For NPV and energy the downside case is therefore P10, which equals the lender-style "P90 exceedance" figure.`;

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check the inputs for zero durations, zero production or rates at or below −100 %.`);
  return v;
}
export const round = (v, n = 7) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const mil = (v) => round(v / 1e6, 7);
export const fM = (v) => `${fmt(v / 1e6)} ${MU}`;
export const pct = (v, n = 3) => `${fmt(v * 100, n)} %`;
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
export const sum = (a) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s; };
export const range = (n, f = (i) => i) => Array.from({ length: n }, (_, i) => f(i));

/** Decimate a sample to at most n values for plotting (keeps order). */
export function thinSample(a, n = 1500) {
  if (a.length <= n) return Array.from(a);
  const out = []; const step = a.length / n;
  for (let i = 0; i < n; i++) out.push(a[Math.floor(i * step)]);
  return out;
}

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.max(demand, 0) / capacity : 99;
  return { name, demand: round(demand), capacity: round(capacity), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
/** Check of a quantity that must be at least a minimum (uc = required / achieved). */
export const checkMin = (name, achieved, required, unit, ref, explain) => {
  const uc = achieved > 0 ? required / achieved : 99;
  return { name, demand: round(required), capacity: round(achieved), unit, uc: Number.isFinite(uc) ? Math.min(Math.max(uc, 0), 99) : 99, ref, explain };
};

/** Input of role `cost` carrying the mandatory placeholder wording. */
export const cost = (o) => ({ type: 'number', role: 'cost', ...o, help: `${o.help ? `${o.help} ` : ''}${PLACEHOLDER}` });

/** Coerce and validate a table input: numeric columns must be finite, optional bounds per column. */
export function tableRows(rows, what, cols, { min = 1 } = {}) {
  need(Array.isArray(rows) && rows.length >= min, `${what}: at least ${min} row${min > 1 ? 's are' : ' is'} required.`);
  return rows.map((r, i) => {
    const o = { ...r };
    for (const [k, spec] of Object.entries(cols)) {
      const v = Number(r[k]);
      need(Number.isFinite(v), `${what}, row ${i + 1}: '${k}' must be a number.`);
      if (spec.min != null) need(v >= spec.min, `${what}, row ${i + 1}: '${k}' must be ≥ ${spec.min}.`);
      if (spec.max != null) need(v <= spec.max, `${what}, row ${i + 1}: '${k}' must be ≤ ${spec.max}.`);
      o[k] = v;
    }
    return o;
  });
}
/** Reject duplicated names in a table (duplicate cost lines silently double count). */
export function uniqueNames(rows, what, key = 'name') {
  const seen = new Set();
  for (const r of rows) {
    const n = String(r[key] ?? '').trim();
    need(n.length > 0, `${what}: every row needs a ${key}.`);
    need(!seen.has(n), `${what}: '${n}' appears twice — merge or rename the duplicate so it is not double counted.`);
    seen.add(n);
  }
}
/** Series input → validated numeric arrays. */
export function seriesXY(s, what, minPts = 3) {
  const x = (s?.x || []).map(Number); const y = (s?.y || []).map(Number);
  need(x.length === y.length && x.length >= minPts, `${what}: at least ${minPts} points with matching x and y are required.`);
  need(x.every(Number.isFinite) && y.every(Number.isFinite), `${what}: all values must be numbers.`);
  return { x, y };
}

/** Distribution spec from a three-point estimate and a shape name. */
export function threePoint(shape, low, mode, high, what = 'three-point estimate') {
  need(low <= mode && mode <= high && high > low, `${what}: need low ≤ most likely ≤ high and high > low (got ${low}, ${mode}, ${high}).`);
  if (shape === 'triangular') return { dist: 'triangular', min: low, mode, max: high };
  if (shape === 'pert') return { dist: 'pert', min: low, mode, max: high };
  if (shape === 'uniform') return { dist: 'uniform', min: low, max: high };
  // normal / lognormal: low and high read as the cumulative P10 and P90 (z = 1.2816)
  if (shape === 'normal') return { dist: 'normal', mean: 0.5 * (low + high), std: (high - low) / (2 * 1.2816) };
  if (shape === 'lognormal') {
    need(low > 0, `${what}: a lognormal needs a positive low value.`);
    return { dist: 'lognormal', mu: 0.5 * (Math.log(low) + Math.log(high)), sigma: (Math.log(high) - Math.log(low)) / (2 * 1.2816) };
  }
  throw new Error(`${what}: distribution '${shape}' is not one of triangular, pert, uniform, normal, lognormal.`);
}
export const SHAPES = [
  { value: 'pert', label: 'Beta-PERT (min, mode, max)' }, { value: 'triangular', label: 'Triangular (min, mode, max)' },
  { value: 'uniform', label: 'Uniform (min, max)' }, { value: 'normal', label: 'Normal (low/high = P10/P90)' },
  { value: 'lognormal', label: 'Lognormal (low/high = P10/P90)' },
];

/**
 * IRR with diagnostics. A flow with late outflows (decommissioning) has two sign changes, so a second, negative
 * root usually exists; when Norstrom's criterion guarantees one positive root that root is the economic IRR.
 * Returns value = null when no rate can be used for a decision.
 */
export function irrReport(cf) {
  const r = irrRoots(Array.from(cf));
  const pos = r.roots.filter((x) => x > 0);
  let value = r.irr; let text;
  if (r.status === 'unique') text = `Unique IRR of ${pct(r.irr)} (${r.descartesSignChanges} sign change${r.descartesSignChanges === 1 ? '' : 's'} in the flows).`;
  else if (r.status === 'none') text = 'No real IRR exists: the cash flow never changes sign in present-value terms at any rate. Use NPV and MIRR.';
  else if (r.norstromUnique && pos.length === 1) {
    value = pos[0];
    text = `${r.roots.length} mathematical IRRs (${r.roots.map((x) => pct(x)).join(', ')}) because the flows change sign ${r.descartesSignChanges} times; Norstrom's test shows a single positive root, ${pct(value)}, which is the economic IRR.`;
  } else text = `${r.roots.length} IRRs (${r.roots.map((x) => pct(x)).join(', ')}): the flows change sign ${r.descartesSignChanges} times and no single rate describes the return. Decide on NPV and MIRR.`;
  return { value, status: r.status, roots: r.roots, count: r.roots.length, signChanges: r.descartesSignChanges, text };
}

/** Plain NPV of a typed or plain array at a constant rate (end-of-period, index 0 = time zero). */
export function pv(cf, rate) {
  let s = 0; let d = 1; const q = 1 / (1 + rate);
  for (let t = 0; t < cf.length; t++) { s += cf[t] * d; d *= q; }
  return s;
}
/**
 * Economic IRR for Monte Carlo use: the highest rate in [−45 %, 100 %] at which NPV changes sign, by a coarse
 * downward scan and bisection (a late outflow such as decommissioning adds a second, strongly negative root that
 * is not the return on the investment). NaN when no sign change is found.
 */
const IRR_SCAN = [1, 0.6, 0.4, 0.3, 0.22, 0.16, 0.12, 0.09, 0.06, 0.03, 0, -0.04, -0.08, -0.13, -0.2, -0.3, -0.45];
export function irrFast(cf) {
  let hi = IRR_SCAN[0]; let fhi = pv(cf, hi);
  for (let j = 1; j < IRR_SCAN.length; j++) {
    let lo = IRR_SCAN[j]; let flo = pv(cf, lo);
    if (flo === 0) return lo;
    if (flo * fhi < 0) {
      for (let i = 0; i < 40; i++) {
        const mid = 0.5 * (lo + hi); const fm = pv(cf, mid);
        if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; }
      }
      return 0.5 * (lo + hi);
    }
    hi = lo; fhi = flo;
  }
  return NaN;
}
export const crf = (i, n) => (Math.abs(i) < 1e-12 ? 1 / n : i / (1 - (1 + i) ** -n));

/* ------------------------------------------------------------------------------------------ fast sampling */

/**
 * Inverse CDF u → x of a distribution spec. Families without a closed-form quantile (PERT, beta, gamma) are
 * inverted from their CDF tabulated on 400 intervals, which is far cheaper than a Newton solve per sample.
 */
export function quantileFn(spec) {
  const d = dist(spec);
  if (!['pert', 'beta', 'gamma'].includes(spec.dist)) return (u) => d.icdf(u);
  const m = 400;
  const lo = Number.isFinite(d.support[0]) ? d.support[0] : d.icdf(1e-7); const hi = Number.isFinite(d.support[1]) ? d.support[1] : d.icdf(1 - 1e-7);
  const xs = new Float64Array(m + 1); const F = new Float64Array(m + 1);
  for (let k = 0; k <= m; k++) { xs[k] = lo + (hi - lo) * k / m; F[k] = d.cdf(xs[k]); }
  F[0] = Math.min(F[0], 0); F[m] = Math.max(F[m], 1);
  return (u) => {
    let a = 0; let b = m;
    while (b - a > 1) { const c = (a + b) >> 1; if (F[c] <= u) a = c; else b = c; }
    const dF = F[b] - F[a];
    return dF > 0 ? xs[a] + (xs[b] - xs[a]) * (u - F[a]) / dF : xs[a];
  };
}
// Sampling-grade normal CDF (Abramowitz & Stegun 26.2.17, |error| < 7.5e-8) and inverse (Acklam, relative error
// < 1.2e-9): accurate far beyond Monte Carlo noise and several times faster than the refined library versions.
function phi(z) {
  const t = 1 / (1 + 0.2316419 * Math.abs(z));
  const p = 0.3989422804014327 * Math.exp(-0.5 * z * z) * t * (0.31938153 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
  return z >= 0 ? 1 - p : p;
}
function phiInv(p) {
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.3577518672690, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  if (p < 0.02425) { const q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p > 1 - 0.02425) { const q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  const q = p - 0.5; const r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}
/** n stratified uniforms (one per stratum, Latin hypercube), shuffled with the seeded generator. */
export function stratified(n, rng) {
  const a = new Float64Array(n);
  for (let k = 0; k < n; k++) a[k] = (k + rng.next()) / n;
  for (let k = n - 1; k > 0; k--) { const j = Math.floor(rng.next() * (k + 1)); const t = a[k]; a[k] = a[j]; a[j] = t; }
  return a;
}
/**
 * Latin-hypercube sample of named distributions with one common rank correlation, imposed through a one-factor
 * Gaussian copula: z_i = √ρ w + √(1 − ρ) e_i with ρ = 2 sin(π ρ_s / 6) for a target Spearman correlation ρ_s
 * (Pearson 1907). Seeded and reproducible. Returns { name: Float64Array(n) }.
 */
export function sampleLHS(vars, n, seed, rhoRank = 0) {
  const rng = makeRng(seed);
  const names = Object.keys(vars);
  const r = rhoRank > 0 ? 2 * Math.sin(Math.PI * rhoRank / 6) : 0;
  const a = Math.sqrt(r); const b = Math.sqrt(1 - r);
  const w = r > 0 ? stratified(n, rng).map((u) => phiInv(Math.min(Math.max(u, 1e-12), 1 - 1e-12))) : null;
  const out = {};
  for (const nm of names) {
    const q = quantileFn(vars[nm]); const e = stratified(n, rng); const col = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      let u = e[k];
      if (w) u = phi(a * w[k] + b * phiInv(Math.min(Math.max(u, 1e-12), 1 - 1e-12)));
      col[k] = q(Math.min(Math.max(u, 1e-12), 1 - 1e-12));
    }
    out[nm] = col;
  }
  return out;
}
/** Sorted copy and percentile reader (cumulative, p in percent) of a sample. */
export function percentiles(sample) {
  const sorted = Float64Array.from(sample).sort();
  const n = sorted.length;
  return (p) => { const h = (n - 1) * p / 100; const i = Math.floor(h); const j = Math.min(i + 1, n - 1); return sorted[i] + (h - i) * (sorted[j] - sorted[i]); };
}
export function meanStd(a) {
  let m = 0; for (let k = 0; k < a.length; k++) m += a[k]; m /= a.length;
  let v = 0; for (let k = 0; k < a.length; k++) v += (a[k] - m) ** 2;
  return { mean: m, std: Math.sqrt(v / Math.max(1, a.length - 1)) };
}
/** Spearman rank correlation of two samples. */
export function rankCorr(a, b) {
  const rk = (x) => { const idx = Array.from(x.keys()).sort((i, j) => x[i] - x[j]); const r = new Float64Array(x.length); idx.forEach((i, k) => { r[i] = k; }); return r; };
  const ra = rk(a); const rb = rk(b); const n = a.length; const mu = (n - 1) / 2;
  let sab = 0; let saa = 0; let sbb = 0;
  for (let k = 0; k < n; k++) { sab += (ra[k] - mu) * (rb[k] - mu); saa += (ra[k] - mu) ** 2; sbb += (rb[k] - mu) ** 2; }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : 0;
}
