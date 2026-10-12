// Shared helpers of the riser suite: unit handling, pipe-section inputs, current profiles, S–N selection,
// screening sustainability factors and advice records. Engineering units are converted to SI here so that every
// module calls the engines in SI (docs/CONTRACT.md §2.1).
import { pipeProperties } from '../../engine/lines/riser.js';
import { snCurve, snConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 25;   // typical installation / drilling vessel fuel burn, t/day (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible riser.`);
  return v;
}
export const round = (v, n = 3) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
export const scaled = (a, f) => Array.from(a, (v) => v * f);

/* ------------------------------------------------------------------------------------------------ inputs */

/** Pipe cross-section, material and contents inputs shared by the modules. */
export function pipeInputs({ D = 273.1, t = 20.6, mat = 'api5l-x65', rhoInt = 700, coat = 3, group = 'Riser pipe', design = true } = {}) {
  return [
    { key: 'D', label: 'Steel outer diameter', unit: 'mm', type: 'number', default: D, min: 50, max: 1500, group, role: 'geometry', bind: 'pipeline.D',
      help: 'Outside diameter of the steel pipe (10¾ in = 273.1 mm).' },
    { key: 't', label: 'Wall thickness', unit: 'mm', type: 'number', default: t, min: 3, max: 80, group, role: 'geometry', design,
      uncertain: { dist: 'normal', cov: 0.03 }, help: 'Nominal wall. Inside diameter = D − 2t.' },
    { key: 'mat', label: 'Line-pipe steel', type: 'select', source: 'materials:steel', default: mat, group, role: 'material',
      help: 'Gives Young’s modulus, Poisson’s ratio, density, yield and tensile strength.' },
    { key: 'rhoInt', label: 'Contents density', unit: 'kg/m³', type: 'number', default: rhoInt, min: 0, max: 2200, group: 'Internal fluid', role: 'load',
      uncertain: { dist: 'normal', cov: 0.08 }, help: 'Oil ≈ 700–850, gas ≈ 100–250, sea water 1025, drilling mud up to 2000.' },
    { key: 'coatT', label: 'External coating thickness', unit: 'mm', type: 'number', default: coat, min: 0, max: 150, group, role: 'geometry',
      help: 'Anti-corrosion or thermal coating; adds weight, buoyancy and drag diameter, no stiffness.' },
    { key: 'coatRho', label: 'Coating density', unit: 'kg/m³', type: 'number', default: 900, min: 100, max: 3500, group, role: 'material' },
  ];
}

export function waterInputs({ depth = 1000 } = {}) {
  return [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: depth, min: 50, max: 3500, group: 'Site', role: 'environment', bind: 'site.depth',
      help: 'Still-water depth at the riser base.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
  ];
}

export function hydroInputs({ Cdn = 1.2, Ca = 1.0 } = {}) {
  return [
    { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: Cdn, min: 0.3, max: 3, group: 'Hydrodynamic coefficients', role: 'load',
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'Bare pipe 0.7–1.2; straked pipe 1.4–2.0; raise it when VIV amplifies drag.' },
    { key: 'Cdt', label: 'Tangential drag coefficient', unit: '-', type: 'number', default: 0.01, min: 0, max: 0.2, group: 'Hydrodynamic coefficients', role: 'load' },
    { key: 'Ca', label: 'Added-mass coefficient', unit: '-', type: 'number', default: Ca, min: 0, max: 2.5, group: 'Hydrodynamic coefficients', role: 'load' },
  ];
}

/**
 * Section record for the rod / riser engines from the form values (SI).
 * Bending and axial stiffness come from the steel wall only; coating adds mass, buoyancy and drag diameter.
 */
export function pipeFrom(inp, ctx, { rhoInt, tOverride, buoyancy } = {}) {
  const mat = ctx.material(inp.mat);
  need(mat, `Steel '${inp.mat}' is not in the materials library.`);
  const D = inp.D / 1e3, t = (tOverride ?? inp.t) / 1e3;
  need(t < D / 2, `Wall thickness ${fmt(t * 1e3)} mm is not smaller than the pipe radius ${fmt(D * 500)} mm.`);
  const rhoW = inp.rhoW ?? 1025;
  const coatings = inp.coatT > 0 ? [{ thickness: inp.coatT / 1e3, rho: inp.coatRho ?? 900 }] : [];
  const props = pipeProperties({ D, t, E: mat.E, rhoSteel: mat.rho, rhoInt: rhoInt ?? inp.rhoInt ?? 0, rhoExt: rhoW, coatings, buoyancy, Ca: inp.Ca ?? 1 });
  const section = { EA: props.EA, EI: props.EI, w: props.w, D: props.Dhydro, Cdn: inp.Cdn ?? 1.2, Cdt: inp.Cdt ?? 0.01, Ai: props.Ai, Ae: props.Ae, m: props.m, ma: props.ma };
  return { D, t, mat, props, section, rhoW };
}

/** Steel mass of a pipe string in tonnes. */
export const steelTonnes = (D, t, L, rho = 7850) => Math.PI * (D - t) * t * L * rho / 1e3;

/** Wall thickness that would bring a stress-type utilisation to the target, assuming stress ∝ 1/t (thin wall). */
export const thicknessForUc = (t, uc, target = 0.9) => t * uc / target;

/** Current profile inputs: surface and near-bed speed with a power-law in between. */
export function currentInputs({ Us = 0.8, Ub = 0.1, bind = 'site.Uc100', group = 'Current' } = {}) {
  return [
    { key: 'Us', label: 'Surface current speed', unit: 'm/s', type: 'number', default: Us, min: 0, max: 3.5, group, role: 'bc', bind,
      uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'In-plane current at the surface (the rod model is two-dimensional: out-of-plane current is not modelled).' },
    { key: 'Ub', label: 'Near-bed current speed', unit: 'm/s', type: 'number', default: Ub, min: 0, max: 2, group, role: 'bc', bind: 'site.UcBed100' },
    { key: 'curExp', label: 'Profile exponent', unit: '-', type: 'number', default: 1.5, min: 0.1, max: 8, group, role: 'bc',
      help: 'U(z) = Ub + (Us − Ub)·(height/depth)^n. n = 1 linear shear; larger n concentrates the current near the surface.' },
  ];
}
/** U as a function of height above the seabed [m] → m/s. */
export function currentFn(inp, depth) {
  const { Us, Ub, curExp } = inp;
  return (h) => Ub + (Us - Ub) * Math.min(1, Math.max(0, h / depth)) ** curExp;
}

/** S–N curve selection (DNV-RP-C203) shared by the fatigue modules. */
export function snInputs({ curve = 'F1', env = 'seawater-cp', scf = 1.2, group = 'Fatigue resistance' } = {}) {
  return [
    { key: 'snCurve', label: 'S–N curve (weld class)', type: 'select', default: curve, options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group, role: 'material',
      help: 'Girth welds: D or E for ground/high-quality welds at the cap, F1 for single-sided root, C1/B1 for machined base material. Confirm the class in the project welding specification.' },
    { key: 'snEnv', label: 'Environment', type: 'select', default: env, options: [{ value: 'air', label: 'Air (internal dry service)' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion / sour' }], group, role: 'material' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: scf, min: 1, max: 5, group, role: 'material', uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'Girth-weld misalignment SCF (DNV-RP-C203 §3.3.7), typically 1.1–1.4 for counter-bored riser joints.' },
  ];
}
export function snFrom(inp, thickness) {
  const curve = run(() => snCurve(inp.snCurve, inp.snEnv), 'S–N curve');
  return { curve, sn: snConstants(curve, { thickness, scf: inp.scf }), label: `DNV-RP-C203 ${inp.snCurve} (${inp.snEnv})` };
}

/* ------------------------------------------------------------------------------------------------ advice */

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? uc : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

/** Sustainability sentence for a change in steel mass [t]. */
export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (confirm with the mill's product declaration).`;
export const vesselCo2Text = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const LIMIT_2D = 'The rod model is two-dimensional (in-plane): out-of-plane current, cross offsets, torsion and lateral seabed friction are not modelled, and the seabed is a row of linear lift-off springs.';
export const CONFIRM = 'Default recalled from the standard — confirm against the governing project edition.';

/* ------------------------------------------------------------------------------------------------ results */

export const kpi = (key, label, value, unit, explain, limit, sense) => (
  limit == null ? { key, label, value, unit, explain } : { key, label, value, unit, limit, sense, explain });
export const plot = (title, xl, xu, yl, yu, series, explain, type = 'line') => (
  { type, title, x: { label: xl, unit: xu }, y: { label: yl, unit: yu }, series, explain });
export const meshInput = (def = 160, levelsHint = '80 / 160 / 320') => ({
  key: 'nElem', label: 'Beam elements along the riser', type: 'int', default: def, min: 20, max: 1280, group: 'Solver', role: 'control',
  help: `Spatial discretisation of the rod model (${levelsHint} in the convergence study). The mesh is graded towards the hang-off and the touchdown zone.` });

/** Reject a result that carries a non-finite headline number, with a message the user can act on. */
export function finish(res) {
  for (const k of res.kpis) finite(k.value, k.label);
  for (const c of res.checks ?? []) if (!Number.isFinite(c.uc) || c.uc < 0) throw new Error(`Check '${c.name}' could not be evaluated. Review the inputs it depends on.`);
  return res;
}

/**
 * Independent hand formula for apparent (submerged) weight per metre [N/m] used by the verification cases:
 * (steel + coating + contents mass − displaced water mass)·g, lengths in mm.
 */
export function handWeight({ D = 273.1, t = 20.6, coatT = 3, coatRho = 900, rhoInt = 700, rhoSteel = 7850, rhoW = 1025 } = {}) {
  const d = D / 1e3, w = t / 1e3, c = coatT / 1e3, q = Math.PI / 4;
  const steel = q * (d * d - (d - 2 * w) ** 2) * rhoSteel, coat = q * ((d + 2 * c) ** 2 - d * d) * coatRho;
  const cont = q * (d - 2 * w) ** 2 * rhoInt, disp = q * (d + 2 * c) ** 2 * rhoW;
  return (steel + coat + cont - disp) * 9.80665;
}

/** Stations table rows from parallel arrays, decimated to about n rows. */
export function stations(cols, n = 40) {
  const keys = Object.keys(cols), len = cols[keys[0]].length, rows = [];
  const step = Math.max(1, Math.floor((len - 1) / n));
  for (let i = 0; i < len; i += step) { const r = {}; for (const k of keys) r[k] = round(cols[k][i], 5); rows.push(r); }
  return rows;
}
export const maxAbs = (a) => { let m = 0; for (const v of a) if (Math.abs(v) > m) m = Math.abs(v); return m; };
export const argMaxAbs = (a) => { let k = 0; for (let i = 1; i < a.length; i++) if (Math.abs(a[i]) > Math.abs(a[k])) k = i; return k; };
