// Generic studies that work on any analysis module through the module contract (docs/CONTRACT.md §2):
// run, verification, sensitivity, parameter sweeps, uncertainty propagation, discretisation convergence,
// calibration against measurements, validation metrics and the design advisor.
// Pure: runs in the page, in a Worker and under Node.
import { makeRng } from '../engine/math/rng.js';
import { dist, metrics as errorMetrics, percentile } from '../engine/math/stats.js';

export const ROLE_LABEL = {
  geometry: 'Geometry', material: 'Material', load: 'Load', environment: 'Environment', soil: 'Soil', ic: 'Initial condition',
  bc: 'Boundary condition', operation: 'Operation', control: 'Analysis control', criteria: 'Acceptance criterion', cost: 'Cost',
};

/** Default input set: module defaults, overridden by bound project-bus values, then by saved user values. */
export function defaults(mod, bus = {}, saved = {}, factors = {}) {
  const out = {};
  const linked = {};
  for (const inp of mod.inputs || []) {
    let v = inp.default;
    if (inp.type === 'table') v = (inp.default || []).map((r) => ({ ...r }));
    if (inp.type === 'series') v = { x: [...(inp.default?.x || [])], y: [...(inp.default?.y || [])] };
    if (inp.bind && Number.isFinite(bus[inp.bind]) && (inp.type === 'number' || inp.type === 'int' || !inp.type)) {
      const conv = fromSI(bus[inp.bind], inp.unit);
      if (Number.isFinite(conv) && (inp.min == null || conv >= inp.min) && (inp.max == null || conv <= inp.max)) { v = conv; linked[inp.key] = inp.bind; }
    }
    // a factor confirmed in the project's code-factor register replaces the built-in default everywhere
    if (inp.role === 'criteria' && Number.isFinite(factors[inp.key])) v = factors[inp.key];
    if (saved && Object.prototype.hasOwnProperty.call(saved, inp.key)) { v = saved[inp.key]; delete linked[inp.key]; }
    out[inp.key] = v;
  }
  return { values: out, linked };
}

// Conversion from the SI bus value to the engineering unit shown in a form field.
const UNIT_FACTORS = {
  mm: 1e-3, cm: 1e-2, km: 1e3, in: 0.0254, ft: 0.3048, 'm²': 1, 'mm²': 1e-6, kN: 1e3, MN: 1e6, t: 1e3, te: 1e3, tonne: 1e3, 'kN·m': 1e3, 'MN·m': 1e6, kNm: 1e3, MNm: 1e6,
  'kPa/m': 1e3, 'MPa/m': 1e6, 'N/mm²': 1e6, 'kg/m': 1, 'N/m': 1, 'kN·m/rad': 1e3, 'MN·m/rad': 1e6, 'GN·m/rad': 1e9, 'GN/m': 1e9, 'kt': 1e6, 'µm': 1e-6, 'kV': 1e3, 'kA': 1e3, 'MVA': 1e6, yr: 1, years: 1,
  kPa: 1e3, MPa: 1e6, GPa: 1e9, bar: 1e5, barg: 1e5, psi: 6894.757, 'kN/m': 1e3, 'MN/m': 1e6, 'kN/m³': 1e3, 'kN/m²': 1e3, kW: 1e3, MW: 1e6, GW: 1e9, 'MWh': 3.6e9, 'GWh': 3.6e12, 'kWh': 3.6e6,
  'GWh/yr': 1e9, 'MWh/yr': 1e6, min: 60, h: 3600, hr: 3600, day: 86400, days: 86400, kn: 0.514444, knots: 0.514444, 'km/h': 1 / 3.6, '%': 0.01, g: 9.80665, 'mm/yr': 1e-3,
};
export function fromSI(v, unit) {
  if (!unit) return v;
  if (unit === 'deg' || unit === '°') return v; // bus angles named *Deg are already degrees
  if (unit === '°C' || unit === 'K') return v;  // bus temperatures are °C by convention
  const f = UNIT_FACTORS[unit];
  return f ? v / f : v;
}

export function validate(mod, values) {
  const errors = {};
  for (const inp of mod.inputs || []) {
    const v = values[inp.key];
    if (hidden(inp, values)) continue;
    if (inp.type === 'number' || inp.type === 'int' || !inp.type) {
      if (typeof v !== 'number' || !Number.isFinite(v)) { errors[inp.key] = 'Enter a number'; continue; }
      if (inp.min != null && v < inp.min) errors[inp.key] = `Must be ≥ ${inp.min}`;
      if (inp.max != null && v > inp.max) errors[inp.key] = `Must be ≤ ${inp.max}`;
      if (inp.type === 'int' && !Number.isInteger(v)) errors[inp.key] = 'Whole number required';
    }
  }
  return errors;
}
export const hidden = (inp, values) => !!inp.showIf && !([].concat(inp.showIf.in ?? inp.showIf.eq)).includes(values[inp.showIf.key]);

/** The exact initial and boundary conditions (and every other input, by role) used for a run. */
export function conditionReport(mod, values) {
  const rows = [];
  for (const inp of mod.inputs || []) {
    if (hidden(inp, values)) continue;
    let v = values[inp.key];
    if (inp.type === 'select') v = inp.options?.find((o) => o.value === v)?.label ?? v;
    else if (inp.type === 'table') v = `${(v || []).length} rows`;
    else if (inp.type === 'series') v = `${(v?.x || []).length} points`;
    else if (inp.type === 'geometry') v = v ? (v.name || 'imported model') : 'none (parametric)';
    else if (inp.type === 'bool') v = v ? 'yes' : 'no';
    rows.push({ role: inp.role || 'operation', label: inp.label, value: v, unit: inp.unit || '', key: inp.key });
  }
  return rows;
}

/**
 * Status of the standards-derived factors (role 'criteria') used in a run: confirmed in the register,
 * changed by the user without a sign-off, or still the built-in default recalled from the standard.
 */
export function factorReport(mod, values, confirmed = {}) {
  const rows = [];
  for (const inp of mod.inputs || []) {
    if (inp.role !== 'criteria' || (inp.type || 'number') !== 'number' || hidden(inp, values)) continue;
    const c = confirmed[`${mod.id}.${inp.key}`];
    const v = values[inp.key];
    rows.push({ key: inp.key, label: inp.label, unit: inp.unit || '', value: v, default: inp.default,
      status: c && c.value === v ? 'confirmed' : v !== inp.default ? 'edited' : 'default', edition: c?.edition || '', by: c?.by || '' });
  }
  return { rows, unconfirmed: rows.filter((r) => r.status !== 'confirmed').length, total: rows.length };
}

export function summarise(result) {
  const kpis = {};
  let fails = 0, maxUC = 0, governing = '';
  for (const k of result.kpis || []) {
    kpis[k.key] = k.value;
    const ok = kpiOk(k);
    if (ok === false) fails++;
  }
  for (const c of result.checks || []) {
    if (Number.isFinite(c.uc) && c.uc > maxUC) { maxUC = c.uc; governing = c.name; }
    if (c.uc > 1) fails++;
  }
  return { status: fails ? 'fail' : 'ok', fails, maxUC, governing, kpis, quantities: result.quantities || null };
}
export function kpiOk(k) {
  if (k.limit == null || !Number.isFinite(k.value)) return null;
  return k.sense === 'min' ? k.value >= k.limit : k.value <= k.limit;
}

async function evaluate(mod, values, ctx) {
  const res = await mod.compute(values, ctx);
  const s = summarise(res);
  return { res, vals: { ...s.kpis, maxUC: s.maxUC }, status: s.status };
}

/** Run the module's verification cases against their independent reference values. */
export async function verify(mod, ctx) {
  const out = [];
  for (const vc of mod.verify || []) {
    if (vc.pending) { out.push({ name: vc.name, status: 'pending', rows: [] }); continue; }
    try {
      const { values } = defaults(mod);
      const { res } = await evaluate(mod, { ...values, ...(vc.inputs || {}) }, ctx);
      const rows = (vc.expect || []).map((e) => {
        const k = (res.kpis || []).find((x) => x.key === e.kpi);
        const got = k?.value;
        const err = e.abs != null ? Math.abs(got - e.value) : Math.abs(got - e.value) / Math.max(Math.abs(e.value), 1e-300);
        const tol = e.abs != null ? e.abs : (e.tol ?? 0.01);
        return { kpi: e.kpi, label: k?.label || e.kpi, unit: k?.unit || '', expected: e.value, got, error: err, tol, absolute: e.abs != null, pass: Number.isFinite(got) && err <= tol };
      });
      out.push({ name: vc.name, source: vc.source || '', status: rows.length && rows.every((r) => r.pass) ? 'pass' : 'fail', rows });
    } catch (err) {
      out.push({ name: vc.name, source: vc.source || '', status: 'fail', error: err.message, rows: [] });
    }
  }
  return out;
}

export const numericInputs = (mod, values) => (mod.inputs || []).filter((i) => (i.type === 'number' || !i.type) && !hidden(i, values) && Number.isFinite(values[i.key]));

/** One-at-a-time sensitivity: each numeric input moved ±pct → swing of every KPI (tornado data). */
export async function sensitivity(mod, values, ctx, { pct = 10, keys, onProgress } = {}) {
  const base = await evaluate(mod, values, ctx);
  const inputs = numericInputs(mod, values).filter((i) => !keys || keys.includes(i.key)).filter((i) => i.role !== 'control');
  const rows = [];
  let n = 0;
  for (const inp of inputs) {
    const v0 = values[inp.key];
    const d = v0 !== 0 ? Math.abs(v0) * pct / 100 : (inp.max != null && inp.min != null ? (inp.max - inp.min) * pct / 100 : pct / 100);
    const lo = clamp(v0 - d, inp), hi = clamp(v0 + d, inp);
    const row = { key: inp.key, name: inp.label, unit: inp.unit, low: lo, high: hi, out: {} };
    for (const [side, v] of [['low', lo], ['high', hi]]) {
      try { row.out[side] = (await evaluate(mod, { ...values, [inp.key]: v }, ctx)).vals; } catch { row.out[side] = null; }
    }
    rows.push(row);
    onProgress?.(++n / inputs.length, inp.label);
  }
  return { base: base.vals, pct, rows, kpis: (base.res.kpis || []).map((k) => ({ key: k.key, label: k.label, unit: k.unit })) };
}
const clamp = (v, inp) => Math.min(inp.max ?? Infinity, Math.max(inp.min ?? -Infinity, v));

/** Parameter sweep over one or two inputs. */
export async function sweep(mod, values, ctx, { key, from, to, n = 21, key2, from2, to2, n2 = 9, onProgress } = {}) {
  const xs = Array.from({ length: n }, (_, i) => from + (to - from) * i / Math.max(1, n - 1));
  const ys = key2 ? Array.from({ length: n2 }, (_, i) => from2 + (to2 - from2) * i / Math.max(1, n2 - 1)) : [null];
  const grid = [];
  let c = 0;
  for (const y of ys) {
    const row = [];
    for (const x of xs) {
      const v = { ...values, [key]: x };
      if (key2) v[key2] = y;
      try { row.push((await evaluate(mod, v, ctx)).vals); } catch { row.push(null); }
      onProgress?.(++c / (xs.length * ys.length));
    }
    grid.push(row);
  }
  return { xs, ys: key2 ? ys : null, grid };
}

/** Uncertainty propagation by Latin hypercube sampling of the selected inputs. */
export async function monteCarlo(mod, values, ctx, { n = 300, seed = 1, vars = {}, onProgress } = {}) {
  const names = Object.keys(vars);
  if (!names.length) throw new Error('Select at least one uncertain input.');
  const rng = makeRng(seed);
  const dists = names.map((k) => dist(toDistSpec(vars[k], values[k])));
  // Latin hypercube: one stratified, independently shuffled column per variable
  const cols = names.map(() => { const p = Array.from({ length: n }, (_, i) => (i + rng.next()) / n); rng.shuffle(p); return p; });
  const samples = { inputs: Object.fromEntries(names.map((k) => [k, new Float64Array(n)])), outputs: {} };
  let failed = 0, nFail = 0;
  for (let i = 0; i < n; i++) {
    const v = { ...values };
    names.forEach((k, j) => {
      const inp = mod.inputs.find((q) => q.key === k);
      v[k] = clamp(dists[j].icdf(Math.min(1 - 1e-12, Math.max(1e-12, cols[j][i]))), inp || {});
      samples.inputs[k][i] = v[k];
    });
    try {
      const e = await evaluate(mod, v, ctx);
      for (const [k, val] of Object.entries(e.vals)) (samples.outputs[k] ||= new Float64Array(n).fill(NaN))[i] = val;
      if (e.status === 'fail') nFail++;
    } catch { failed++; }
    if (i % 5 === 0) onProgress?.((i + 1) / n);
  }
  const stats = {};
  for (const [k, arr] of Object.entries(samples.outputs)) {
    const a = Array.from(arr).filter(Number.isFinite);
    if (!a.length) continue;
    const mean = a.reduce((s, x) => s + x, 0) / a.length;
    const sd = Math.sqrt(a.reduce((s, x) => s + (x - mean) ** 2, 0) / Math.max(1, a.length - 1));
    stats[k] = { mean, std: sd, p05: percentile(a, 5), p10: percentile(a, 10), p50: percentile(a, 50), p90: percentile(a, 90), p95: percentile(a, 95), min: a.reduce((m, x) => Math.min(m, x), Infinity), max: a.reduce((m, x) => Math.max(m, x), -Infinity), n: a.length };
  }
  // rank correlation of each input with each output: which uncertainty matters
  const influence = {};
  for (const out of Object.keys(stats)) influence[out] = names.map((k) => ({ key: k, rho: spearman(samples.inputs[k], samples.outputs[out]) })).sort((a, b) => Math.abs(b.rho) - Math.abs(a.rho));
  const pFail = nFail / Math.max(1, n - failed);
  return { n, seed, failedRuns: failed, pFail, pFailSE: Math.sqrt(pFail * (1 - pFail) / Math.max(1, n - failed)), stats, influence, samples };
}
function toDistSpec(v, base) {
  if (v.dist === 'uniform' || v.dist === 'triangular' || v.dist === 'pert') return { dist: v.dist, min: v.min, max: v.max, mode: v.mode ?? base };
  const mean = v.mean ?? base;
  return { dist: v.dist || 'normal', mean, std: v.std ?? Math.abs(mean) * (v.cov ?? 0.1) };
}
function spearman(a, b) {
  const idx = [];
  for (let i = 0; i < a.length; i++) if (Number.isFinite(a[i]) && Number.isFinite(b[i])) idx.push(i);
  if (idx.length < 3) return 0;
  const rank = (arr) => { const o = idx.slice().sort((p, q) => arr[p] - arr[q]); const r = new Map(); o.forEach((id, k) => r.set(id, k)); return idx.map((id) => r.get(id)); };
  const ra = rank(a), rb = rank(b), m = (idx.length - 1) / 2;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < idx.length; i++) { sab += (ra[i] - m) * (rb[i] - m); saa += (ra[i] - m) ** 2; sbb += (rb[i] - m) ** 2; }
  return saa && sbb ? sab / Math.sqrt(saa * sbb) : 0;
}

/**
 * Grid Convergence Index on three successive levels (Roache; Celik et al. 2008, J. Fluids Eng. 130, 078001).
 * f1 is the finest solution; h1 < h2 < h3 are representative sizes.
 */
export function gci3(h1, h2, h3, f1, f2, f3, Fs = 1.25) {
  const r21 = h2 / h1, r32 = h3 / h2;
  const e21 = f2 - f1, e32 = f3 - f2;
  if (Math.abs(e21) < 1e-14 * Math.max(1, Math.abs(f1))) return { p: Infinity, fExt: f1, gciFine: 0, ea21: 0, behaviour: 'converged', asymptotic: 1 };
  const ratio = e32 / e21;
  const behaviour = ratio < 0 ? 'oscillatory' : Math.abs(ratio) < 1 ? 'divergent' : 'monotone';
  const s = Math.sign(ratio) || 1;
  let p = Math.abs(Math.log(Math.abs(ratio))) / Math.log(r21);
  for (let it = 0; it < 60; it++) {
    const q = Math.log((r21 ** p - s) / (r32 ** p - s));
    const pn = Math.abs(Math.log(Math.abs(ratio)) + q) / Math.log(r21);
    if (!Number.isFinite(pn)) break;
    if (Math.abs(pn - p) < 1e-10) { p = pn; break; }
    p = pn;
  }
  const pUse = Math.min(Math.max(p, 0.5), 8); // guard the extrapolation against a noisy observed order
  const fExt = (r21 ** pUse * f1 - f2) / (r21 ** pUse - 1);
  const ea21 = Math.abs((f1 - f2) / (f1 || 1e-300));
  const ea32 = Math.abs((f2 - f3) / (f2 || 1e-300));
  const gciFine = Fs * ea21 / (r21 ** pUse - 1);
  const gciCoarse = Fs * ea32 / (r32 ** pUse - 1);
  return { p, pUsed: pUse, fExt, ea21, gciFine, gciCoarse, behaviour, asymptotic: gciCoarse / (r21 ** pUse * gciFine || 1e-300), r21, r32 };
}

/**
 * Discretisation convergence study on the module's declared control (mesh or meshTime).
 * Runs every level, reports observed order, Richardson extrapolation and GCI for each consecutive triplet,
 * and selects the coarsest level whose estimated discretisation error meets the tolerance.
 */
export async function convergence(mod, values, ctx, { which = 'mesh', levels, qoi, tol = 0.02, onProgress } = {}) {
  const spec = mod[which];
  if (!spec) throw new Error('This module has no discretisation control: its result is closed-form.');
  const lv = (levels && levels.length ? levels : spec.levels).slice().sort((a, b) => a - b);
  if (lv.length < 3) throw new Error('At least three refinement levels are needed.');
  const key = qoi || spec.qoi;
  const runs = [];
  for (let i = 0; i < lv.length; i++) {
    const t0 = Date.now();
    const e = await evaluate(mod, { ...values, [spec.key]: lv[i] }, ctx);
    const kpi = (e.res.kpis || []).find((k) => k.key === key);
    if (!kpi) throw new Error(`Quantity of interest '${key}' was not returned by the module.`);
    // control value N (elements, panels, steps): representative size h ∝ 1/N; a control that is itself a size uses h ∝ value
    runs.push({ level: lv[i], h: spec.isSize ? lv[i] : 1 / lv[i], value: kpi.value, label: kpi.label, unit: kpi.unit, ms: Date.now() - t0, all: e.vals });
    onProgress?.((i + 1) / lv.length, `${spec.label || spec.key} = ${lv[i]}`);
  }
  runs.sort((a, b) => b.h - a.h); // coarse → fine
  const triplets = [];
  for (let i = 0; i + 2 < runs.length; i++) {
    const [c, m, f] = [runs[i], runs[i + 1], runs[i + 2]];
    triplets.push({ levels: [c.level, m.level, f.level], ...gci3(f.h, m.h, c.h, f.value, m.value, c.value) });
  }
  const last = triplets[triplets.length - 1];
  const fExt = last.fExt;
  for (const r of runs) r.errEst = Math.abs((r.value - fExt) / (fExt || 1e-300));
  // selection: coarsest level within tolerance of the extrapolated value, provided the sequence behaves
  const trustworthy = last.behaviour === 'monotone' || last.behaviour === 'converged';
  let selected = runs[runs.length - 1];
  if (trustworthy) for (const r of runs) if (r.errEst <= tol) { selected = r; break; }
  return { which, key: spec.key, controlLabel: spec.label || spec.key, qoi: key, qoiLabel: runs[0].label, unit: runs[0].unit, tol, runs, triplets, fExt, observedOrder: last.p, gciFine: last.gciFine, behaviour: last.behaviour, asymptotic: last.asymptotic, selected: selected.level, selectedErr: selected.errEst, trustworthy, met: trustworthy && selected.errEst <= tol };
}

/**
 * Calibration: adjust chosen inputs within bounds so that predicted KPIs match measurements.
 * data rows: { set: { inputKey: value, ... }, kpi, value, sigma? }. Weighted least squares by Nelder–Mead.
 */
export async function calibrate(mod, values, ctx, { params, data, onProgress } = {}) {
  if (!params?.length) throw new Error('Choose at least one parameter to calibrate.');
  if (!data?.length) throw new Error('Provide measured data.');
  if (data.length < params.length) throw new Error(`Not identifiable: ${params.length} parameters but only ${data.length} measurements.`);
  const predict = async (x) => {
    const v = { ...values };
    params.forEach((p, i) => { v[p.key] = x[i]; });
    const out = [];
    for (const d of data) {
      try { out.push((await evaluate(mod, { ...v, ...(d.set || {}) }, ctx)).vals[d.kpi]); } catch { out.push(NaN); }
    }
    return out;
  };
  const cost = (pred) => pred.reduce((s, pv, i) => s + (Number.isFinite(pv) ? ((pv - data[i].value) / (data[i].sigma || Math.abs(data[i].value) || 1)) ** 2 : 1e6), 0);
  const x0 = params.map((p) => values[p.key]);
  const before = await predict(x0);
  // module evaluations may be async, so the simplex search below awaits its objective
  const res = await asyncNelderMead(async (x) => cost(await predict(x)), x0, params.map((p) => [p.min, p.max]), { maxIter: 60 * params.length, onProgress });
  const after = await predict(res.x);
  // local identifiability: finite-difference sensitivities → parameter correlation matrix
  const J = [];
  for (let j = 0; j < params.length; j++) {
    const dx = Math.max(1e-6, Math.abs(res.x[j]) * 1e-3);
    const xp = res.x.slice(); xp[j] += dx;
    const pp = await predict(xp);
    J.push(pp.map((v, i) => (v - after[i]) / dx / (data[i].sigma || Math.abs(data[i].value) || 1)));
  }
  const JTJ = params.map((_, a) => params.map((__, b) => J[a].reduce((s, v, i) => s + v * J[b][i], 0)));
  const corr = params.map((_, a) => params.map((__, b) => JTJ[a][b] / Math.sqrt((JTJ[a][a] || 1e-300) * (JTJ[b][b] || 1e-300))));
  const atBound = params.map((p, i) => Math.abs(res.x[i] - p.min) < 1e-9 * Math.max(1, Math.abs(p.min)) || Math.abs(res.x[i] - p.max) < 1e-9 * Math.max(1, Math.abs(p.max)));
  const obs = data.map((d) => d.value);
  return {
    params: params.map((p, i) => ({ key: p.key, label: mod.inputs.find((q) => q.key === p.key)?.label || p.key, initial: x0[i], fitted: res.x[i], min: p.min, max: p.max, atBound: atBound[i] })),
    iterations: res.iterations, costBefore: cost(before), costAfter: res.f, before, after, observed: obs,
    metricsBefore: errorMetrics(before, obs), metricsAfter: errorMetrics(after, obs), correlation: corr,
    warnings: [
      ...(atBound.some(Boolean) ? ['A fitted parameter sits on its bound: the measurements are pulling it outside the physically justified range. Check the model, geometry and boundary conditions before trusting this fit.'] : []),
      ...(corr.some((r, a) => r.some((c, b) => a !== b && Math.abs(c) > 0.95)) ? ['Two parameters are almost perfectly correlated (|ρ| > 0.95): the data cannot separate them. Fix one of them or add measurements of a different kind.'] : []),
    ],
  };
}

async function asyncNelderMead(f, x0, bounds, { maxIter = 200, tol = 1e-9, onProgress } = {}) {
  const n = x0.length;
  const clip = (x) => x.map((v, i) => Math.min(bounds[i][1] ?? Infinity, Math.max(bounds[i][0] ?? -Infinity, v)));
  let simplex = [clip(x0)];
  for (let i = 0; i < n; i++) {
    const x = x0.slice();
    const span = Number.isFinite(bounds[i][1] - bounds[i][0]) ? (bounds[i][1] - bounds[i][0]) * 0.1 : Math.abs(x[i]) * 0.1 || 0.1;
    x[i] += x[i] + span > (bounds[i][1] ?? Infinity) ? -span : span;
    simplex.push(clip(x));
  }
  let fs = [];
  for (const s of simplex) fs.push(await f(s));
  let it = 0;
  for (; it < maxIter; it++) {
    const order = fs.map((v, i) => i).sort((a, b) => fs[a] - fs[b]);
    simplex = order.map((i) => simplex[i]); fs = order.map((i) => fs[i]);
    if (Math.abs(fs[n] - fs[0]) <= tol * (Math.abs(fs[0]) + tol)) break;
    const c = Array.from({ length: n }, (_, j) => simplex.slice(0, n).reduce((s, x) => s + x[j], 0) / n);
    const along = (t) => clip(c.map((cj, j) => cj + t * (simplex[n][j] - cj)));
    const xr = along(-1), fr = await f(xr);
    if (fr < fs[0]) { const xe = along(-2), fe = await f(xe); if (fe < fr) { simplex[n] = xe; fs[n] = fe; } else { simplex[n] = xr; fs[n] = fr; } }
    else if (fr < fs[n - 1]) { simplex[n] = xr; fs[n] = fr; }
    else {
      const xc = along(fr < fs[n] ? -0.5 : 0.5), fc = await f(xc);
      if (fc < Math.min(fr, fs[n])) { simplex[n] = xc; fs[n] = fc; }
      else for (let i = 1; i <= n; i++) { simplex[i] = clip(simplex[i].map((v, j) => simplex[0][j] + 0.5 * (v - simplex[0][j]))); fs[i] = await f(simplex[i]); }
    }
    onProgress?.(Math.min(0.99, it / maxIter));
  }
  const best = fs.indexOf(Math.min(...fs));
  return { x: simplex[best], f: fs[best], iterations: it };
}

/** Validation: compare predictions with independent measurements that were not used for calibration. */
export async function validation(mod, values, ctx, { data, acceptance = 0.1 } = {}) {
  if (!data?.length) throw new Error('Provide independent measured data.');
  const pred = [];
  for (const d of data) {
    try { pred.push((await evaluate(mod, { ...values, ...(d.set || {}) }, ctx)).vals[d.kpi]); } catch { pred.push(NaN); }
  }
  const obs = data.map((d) => d.value);
  const m = errorMetrics(pred, obs);
  const rel = pred.map((p, i) => (p - obs[i]) / (Math.abs(obs[i]) || 1e-300));
  const within = rel.filter((r) => Math.abs(r) <= acceptance).length / rel.length;
  // under-prediction of the largest measured values is what hurts a design, whatever the average error says
  const top = obs.map((o, i) => [o, i]).sort((a, b) => Math.abs(b[0]) - Math.abs(a[0])).slice(0, Math.max(1, Math.ceil(obs.length * 0.2)));
  const peakBias = top.reduce((s, [, i]) => s + rel[i], 0) / top.length;
  return { predicted: pred, observed: obs, relative: rel, metrics: m, within, acceptance, peakBias, pass: within >= 0.9 && Math.abs(m.bias) <= acceptance * Math.max(...obs.map(Math.abs)) && peakBias > -acceptance };
}

/**
 * Design advisor: for a run that fails, find for each designated design variable the smallest change that
 * brings every check and limit back within its criterion (bisection on a monotone search along the variable).
 */
export async function advise(mod, values, ctx, { onProgress } = {}) {
  const base = await evaluate(mod, values, ctx);
  const designVars = (mod.inputs || []).filter((i) => i.design && (i.type === 'number' || !i.type) && !hidden(i, values));
  const options = [];
  const ok = async (v) => { try { return (await evaluate(mod, v, ctx)).status === 'ok'; } catch { return false; } };
  let n = 0;
  if (base.status === 'fail') {
    for (const inp of designVars) {
      const v0 = values[inp.key];
      let found = null;
      for (const dir of [1, -1]) {
        const end = dir > 0 ? (inp.max ?? v0 * 4) : (inp.min ?? v0 / 4);
        if (end === v0 || !(await ok({ ...values, [inp.key]: end }))) continue;
        let a = v0, b = end;
        for (let it = 0; it < 26; it++) { const mid = 0.5 * (a + b); if (await ok({ ...values, [inp.key]: mid })) b = mid; else a = mid; if (Math.abs(b - a) <= 1e-4 * Math.abs(end - v0)) break; }
        if (!found || Math.abs(b - v0) < Math.abs(found.value - v0)) found = { value: b, dir };
      }
      if (found) {
        const after = await evaluate(mod, { ...values, [inp.key]: found.value }, ctx);
        options.push({ key: inp.key, label: inp.label, unit: inp.unit, from: v0, to: found.value, change: (found.value - v0) / (Math.abs(v0) || 1), maxUC: after.vals.maxUC, quantities: after.res.quantities || null });
      }
      onProgress?.(++n / designVars.length, inp.label);
    }
    options.sort((a, b) => Math.abs(a.change) - Math.abs(b.change));
  } else {
    // passing design: how much margin can be given back on each design variable before the first criterion binds
    for (const inp of designVars) {
      const v0 = values[inp.key];
      let best = null;
      for (const dir of [-1, 1]) {
        const end = dir > 0 ? (inp.max ?? v0 * 4) : (inp.min ?? v0 / 4);
        if (end === v0 || (await ok({ ...values, [inp.key]: end }))) continue;
        let a = v0, b = end;
        for (let it = 0; it < 26; it++) { const mid = 0.5 * (a + b); if (await ok({ ...values, [inp.key]: mid })) a = mid; else b = mid; if (Math.abs(b - a) <= 1e-4 * Math.abs(end - v0)) break; }
        if (!best || Math.abs(a - v0) > Math.abs(best.value - v0)) best = { value: a };
      }
      if (best && Math.abs(best.value - v0) > 1e-3 * Math.abs(v0 || 1)) options.push({ key: inp.key, label: inp.label, unit: inp.unit, from: v0, to: best.value, change: (best.value - v0) / (Math.abs(v0) || 1), limitOf: true });
      onProgress?.(++n / designVars.length, inp.label);
    }
  }
  return { status: base.status, maxUC: base.vals.maxUC, designVars: designVars.length, options };
}
