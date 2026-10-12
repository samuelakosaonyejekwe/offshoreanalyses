// Shared helpers of the fixed-platform suite: unit handling, advice / check records, screening sustainability
// factors and the caveat texts repeated in module limits and input help. Engineering units are converted to SI
// inside each module (docs/CONTRACT.md §2.1); the helpers here are unit-neutral unless stated.

export const G = 9.80665;
export const DEG = Math.PI / 180;
// Screening emission factors (order of magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 35;   // heavy-lift / construction vessel fuel burn, t/day (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
/** Run an engine call and turn range / convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible structure.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(Math.max(v, 0), hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : 99;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: cap(uc), ref, explain };
};
/** Check written as required factor over achieved factor (reserve-type criteria). */
export const fsCheck = (name, achieved, required, unit, ref, explain) => (
  { name, demand: round(required, 4), capacity: round(cap(achieved, 999), 4), unit, uc: achieved > 0 ? cap(required / achieved) : 99, ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

export const steelCo2 = (t) => `${fmt(Math.abs(t))} t of steel ≈ ${fmt(Math.abs(t) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (a screening factor — confirm with the mill's product declaration).`;
export const vesselCo2 = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

// Caveats required by the engine authors, repeated wherever the corresponding engine is used.
export const CONFIRM = 'Value recalled from the standard — confirm against the governing project edition.';
export const CODE_CAVEAT = 'Tubular member and joint equations and partial factors are implemented as recalled from ISO 19902 / API RP 2A and are overridable criteria: confirm them against the governing edition before relying on a utilisation.';
export const MORISON_CAVEAT = 'Morison loading is valid for slender members (D/λ < 0.2); diffraction on large caissons and CFD are outside the suite.';
export const PUSHOVER_CAVEAT = 'The pushover uses plastic hinges at member ends with a simple constant post-buckling force, not full nonlinear shell finite elements; explicit nonlinear FE is outside the suite.';
export const BASIS_MEMBER = 'ISO 19902:2020 §13 (tubular members); API RP 2A-WSD 22nd ed. §6';
export const BASIS_COST = 'AACE RP 18R-97 estimate classes; whole-life costing per ISO 15686-5';
