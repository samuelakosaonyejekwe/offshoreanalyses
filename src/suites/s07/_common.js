// Shared helpers of the marine-structures suite: constants, unit handling, sea-level-rise inputs, check and
// advice records. Engineering units are converted to SI inside each module (docs/CONTRACT.md §2.1).

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;     // t CO2e per t steel
export const CO2_CONC = 0.3;      // t CO2e per m³ structural concrete
export const CO2_MGO = 3.206;     // t CO2 per t marine gas oil burnt
export const CO2_ROCK = 0.02;     // t CO2e per t quarried and placed rock; 0.01–0.05 depending on haul (screening)
export const MGO_T_PER_DAY = 12;  // typical rock-placing / crane barge spread fuel burn, t/day (screening)

export const EDITION_NOTE = 'Coefficients are carried as recalled from the manuals and are overridable: confirm them against the current edition of EurOtop and the CIRIA Rock Manual before design use.';

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible structure and sea state.`);
  return v;
}
export const r = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const f = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

/** Piecewise-linear interpolation with end clamping (xs ascending). */
export function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1; while (xs[i] < x) i++;
  const t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
  return ys[i - 1] + t * (ys[i] - ys[i - 1]);
}

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}

/** Code-check record; utilisation is demand / capacity, kept finite and non-negative. */
export function check(name, demand, capacity, unit, ref, explain) {
  const d = Math.max(0, demand);
  const uc = capacity > 0 ? d / capacity : (d > 0 ? 99 : 0);
  return { name, demand: r(d), capacity: r(capacity), unit, uc: r(Math.min(uc, 99)), ref, explain };
}
/** Check expressed as a required safety factor: uc = required / achieved. */
export function sfCheck(name, sf, required, ref, explain) {
  const s = Number.isFinite(sf) ? Math.max(sf, 1e-6) : 1e6;
  return { name, demand: r(required), capacity: r(Math.min(s, 999)), unit: '-', uc: r(Math.min(required / s, 99)), ref, explain };
}

export const eng = (severity, text, action, basis, sustainability) => ({ kind: 'engineering', severity, text, action, basis, ...(sustainability ? { sustainability } : {}) });
export const eco = (severity, text, action, basis, sustainability) => ({ kind: 'economic', severity, text, action, basis, ...(sustainability ? { sustainability } : {}) });

/* ------------------------------------------------------------------------------------- sea-level rise */

/** Design sea-level-rise allowance and the scenario table used for the adaptation pathway. */
export function slrInputs({ group = 'Water level', slr = 0.3 } = {}) {
  return [
    { key: 'slr', label: 'Sea-level-rise allowance in the design case', unit: 'm', type: 'number', default: slr, min: 0, max: 3, group, role: 'bc',
      uncertain: { dist: 'lognormal', cov: 0.4 },
      help: 'Added to the still-water level for the checks. Use the rise expected at the end of the design life for the adopted climate scenario.' },
    { key: 'slrScen', label: 'Sea-level-rise scenarios', type: 'table', group, role: 'bc',
      columns: [{ key: 'name', label: 'Scenario', type: 'text' }, { key: 'rise', label: 'Rise', unit: 'm', type: 'number' }],
      default: [{ name: 'Present day', rise: 0 }, { name: '2050 central', rise: 0.25 }, { name: '2100 intermediate', rise: 0.6 }, { name: '2100 high', rise: 1.0 }, { name: '2100 high-end (H++)', rise: 1.6 }],
      help: 'Placeholders: replace with the regional projections of the current IPCC assessment or national guidance. Each row is re-analysed to build the adaptation pathway.' },
  ];
}
/** Validated scenario rows [{ name, rise }] (always includes at least the design allowance). */
export function slrRows(inp) {
  const rows = (inp.slrScen || []).filter((s) => s && Number.isFinite(Number(s.rise))).map((s, i) => ({ name: String(s.name || `Scenario ${i + 1}`), rise: Number(s.rise) }));
  for (const s of rows) need(s.rise >= 0 && s.rise <= 5, `Sea-level-rise scenario '${s.name}': rise must be between 0 and 5 m.`);
  if (!rows.length) rows.push({ name: 'Design allowance', rise: inp.slr || 0 });
  return rows.slice(0, 12);
}
export const SLR_BASIS = 'IPCC AR6 WG1 (2021) ch. 9 sea-level projections; adaptive pathways per ISO 14090:2019 and national coastal-allowance guidance (confirm the current edition)';
