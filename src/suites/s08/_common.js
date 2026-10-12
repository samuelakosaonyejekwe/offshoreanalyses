// Shared helpers of the subsea pipeline suite: the common pipe build-up inputs and their conversion to SI section
// and weight records, the route (chainage / elevation) input, result builders and the standing caveats of the
// engines. Engineering units are converted to SI here so every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { pipeSection, pipeWeights } from '../../engine/pipeline/pipeprops.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const BAR = 1e5;
export const MPA = 1e6;
export const KELVIN = 273.15;
export const YEAR = 365.25 * 86400;
// Screening emission factors (cradle-to-gate order of magnitude, not project EPD values).
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m³ concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil
export const MGO_T_PER_DAY = 35;   // pipelay / construction vessel fuel burn, t/day (screening)
export const KM_PER_DAY = 3;       // S-lay progress for a 12-inch line, km/day (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible pipe and load case.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / Math.max(1, n - 1));
export const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
export const arr = (a, f = 1, n = 4) => Array.from(a, (v) => round(v * f, n));

/* ------------------------------------------------------------------------------------------- standing text */

export const CONFIRM = 'Recalled default — confirm against the edition of the standard governing the project.';
export const LIMIT_F101 = 'DNV-ST-F101 partial factors, test-pressure factors, strain factors and de-rating values are recalled defaults exposed as criteria inputs; confirm them against the governing edition.';
export const LIMIT_NO_FE = 'Closed-form and one-dimensional models only: explicit non-linear finite-element analysis (3-D route, plasticity, contact) is outside this module.';
export const LIMIT_NO_CFD = '3-D CFD, CFD–DEM deposition and full two-fluid six/seven-equation models are outside this module.';
export const STD_F101 = 'DNV-ST-F101 (2021)';
export const ECON_BASIS = 'Quantity-based estimating per AACE RP 18R-97; whole-life costing per ISO 15686-5';
export const RBI_BASIS = 'Risk-based inspection per DNV-RP-G101 / API RP 580; whole-life costing per ISO 15686-5';

/* ------------------------------------------------------------------------------------------- result builders */

export const kpi = (key, label, value, unit, explain, limit, sense) => {
  const k = { key, label, value: finite(value, label), unit, explain };
  if (limit != null) { k.limit = limit; k.sense = sense || 'max'; }
  return k;
};
export const chk = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : (Math.abs(demand) > 0 ? 99 : 0);
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
/** Check whose utilisation is supplied directly (interaction equations). */
export const chkUc = (name, uc, ref, explain) => ({ name, demand: round(uc, 5), capacity: 1, unit: '-', uc: Number.isFinite(uc) ? clamp(uc, 0, 99) : 99, ref, explain });
export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const plot = (type, title, xl, xu, yl, yu, series, explain, extra = {}) => (
  { type, title, x: { label: xl, unit: xu }, y: { label: yl, unit: yu }, series, explain, ...extra });
export const bar = (title, categories, yl, yu, series, explain) => (
  { type: 'bar', title, categories, y: { label: yl, unit: yu }, series, explain });
export const tbl = (title, columns, rows, explain) => ({ title, columns, rows, explain });
export const col = (key, label, unit) => (unit ? { key, label, unit } : { key, label });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc) => (uc > 1 ? 'fail' : uc > 0.9 ? 'warn' : 'ok');
export const num = (key, label, unit, def, min, max, group, role, help, extra = {}) => (
  { key, label, unit, type: 'number', default: def, min, max, group, role, help, ...extra });
export const sel = (key, label, def, options, group, role, help, extra = {}) => (
  { key, label, type: 'select', default: def, options: options.map((o) => (Array.isArray(o) ? { value: o[0], label: o[1] } : { value: o, label: String(o) })), group, role, help, ...extra });

/* ------------------------------------------------------------------------------------------- pipe build-up */

const GROUP = 'Pipe';
/**
 * Common pipe inputs. `opts.omit` lists keys to leave out; `opts.defaults` overrides defaults by key.
 * Default asset: 12-inch (323.9 mm) X65 flowline, 14.3 mm wall, 3 mm 3LPP, 40 mm concrete weight coating.
 */
export function pipeInputs({ omit = [], defaults = {} } = {}) {
  const list = [
    num('D', 'Steel outer diameter', 'mm', 323.9, 50, 1600, GROUP, 'geometry', 'Outside diameter of the steel pipe (12 in = 323.9 mm).', { bind: 'pipeline.D' }),
    num('t', 'Nominal wall thickness', 'mm', 14.3, 3, 80, GROUP, 'geometry', 'Nominal steel wall thickness before tolerances and corrosion allowance.', { bind: 'pipeline.t', design: true, uncertain: { dist: 'normal', cov: 0.02 } }),
    { key: 'steel', label: 'Line-pipe steel', type: 'select', source: 'materials:steel', default: 'api5l-x65', group: GROUP, role: 'material',
      help: 'Supplies SMYS, SMTS, Young\'s modulus, Poisson\'s ratio, density and thermal expansion coefficient.' },
    num('tCoat', 'Anti-corrosion coating thickness', 'mm', 3, 0, 120, GROUP, 'geometry', 'Corrosion or thermal insulation coating directly on the steel (3LPP ≈ 3 mm; wet insulation 30–80 mm).'),
    num('rhoCoat', 'Coating density', 'kg/m³', 900, 300, 2000, GROUP, 'material', 'Polypropylene ≈ 900, FBE ≈ 1400, syntactic PP foam ≈ 700 kg/m³.'),
    num('tConc', 'Concrete weight-coating thickness', 'mm', 40, 0, 200, GROUP, 'geometry', 'Zero for a line without weight coating.', { design: true }),
    { key: 'concrete', label: 'Concrete weight coating', type: 'select', source: 'materials:concrete', default: 'cwc-3040', group: GROUP, role: 'material',
      help: 'Supplies the density and water absorption of the weight coating.' },
    num('rhoCont', 'Contents density', 'kg/m³', 800, 0, 1300, GROUP, 'operation', 'Density of the product in operation (0 = empty, 1025 = sea-water flooded).', { uncertain: { dist: 'normal', cov: 0.05 } }),
    num('rhoW', 'Sea-water density', 'kg/m³', 1025, 990, 1050, GROUP, 'environment', 'Used for buoyancy and hydrodynamic loads.', { bind: 'site.rhoW' }),
    num('wsGiven', 'Submerged weight override', 'N/m', 0, 0, 50000, GROUP, 'load', 'As-built submerged weight in operation. Zero = calculate from the build-up above.', { bind: 'pipeline.submergedWeight' }),
  ];
  return list.filter((i) => !omit.includes(i.key)).map((i) => (defaults[i.key] !== undefined ? { ...i, default: defaults[i.key] } : i));
}

const STEEL_FALLBACK = { E: 207e9, nu: 0.3, rho: 7850, alpha: 1.17e-5, fy: 450e6, fu: 535e6, name: 'API 5L X65' };

/**
 * SI pipe record from the common inputs: section, material, weights.
 * @returns {{ D, t, Di, As, Ai, Ae, I, Z, E, nu, alpha, smys, smts, rhoS, EI, EA, Do, W, ws, wsEmpty, wsFlood, buoy,
 *   massOp, steelKgM, concM3M, grade }} lengths m, areas m², forces N/m, masses kg/m
 */
export function pipeFrom(inp, ctx) {
  const D = inp.D / 1e3, t = inp.t / 1e3;
  need(t < D / 2.2, `Wall thickness ${inp.t} mm is not possible for a ${inp.D} mm pipe: it must be well below the radius.`);
  const m = { ...STEEL_FALLBACK, ...((ctx?.material && ctx.material(inp.steel)) || {}) };
  const conc = (ctx?.material && ctx.material(inp.concrete)) || { rho: 3040, waterAbsorption: 0.03 };
  const sec = run(() => pipeSection({ D, t }), 'Pipe section');
  const coatings = [];
  if ((inp.tCoat ?? 0) > 0) coatings.push({ name: 'coating', thickness: inp.tCoat / 1e3, density: inp.rhoCoat ?? 900 });
  if ((inp.tConc ?? 0) > 0) coatings.push({ name: 'concrete', thickness: inp.tConc / 1e3, density: conc.rho, absorption: conc.waterAbsorption ?? 0 });
  const rhoW = inp.rhoW ?? 1025;
  const W = run(() => pipeWeights({ D, t, rhoSteel: m.rho, coatings, rhoContents: inp.rhoCont ?? 0, rhoWater: rhoW }), 'Weight build-up');
  const ws = (inp.wsGiven ?? 0) > 0 ? inp.wsGiven : W.submergedOperating;
  const Do = W.Dtotal;
  const rConcIn = D / 2 + (inp.tCoat ?? 0) / 1e3;
  return {
    D, t, Di: sec.Di, As: sec.As, Ai: sec.Ai, Ae: sec.Ae, I: sec.I, Z: sec.Z, sec,
    E: m.E, nu: m.nu, alpha: m.alpha, smys: m.fy, smts: m.fu, rhoS: m.rho, grade: m.name || inp.steel,
    EI: m.E * sec.I, EA: m.E * sec.As, Do, W, ws, wsEmpty: W.submergedEmpty, wsFlood: W.submergedFlooded, buoy: W.buoyancy,
    massOp: W.massOperating, steelKgM: W.massSteel, concM3M: Math.PI * (Do * Do / 4 - rConcIn * rConcIn), coatings, concRho: conc.rho, concAbs: conc.waterAbsorption ?? 0, rhoW,
  };
}

/** Steel mass change per metre for a wall-thickness change dt [m] → kg/m. */
export const steelPerMetre = (p, dt) => p.rhoS * Math.PI * (p.D - p.t) * dt;
export const steelCo2 = (tonnes) => `${fmt(Math.abs(tonnes))} t of line-pipe steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (replace with the mill's product declaration).`;
export const concreteCo2 = (m3) => `${fmt(Math.abs(m3))} m³ of concrete weight coating ≈ ${fmt(Math.abs(m3) * CO2_CONCRETE)} t CO₂e at the screening factor of ${CO2_CONCRETE} t CO₂e per m³.`;
export const vesselCo2 = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const lengthInput = (def = 10, group = 'Route') => num('Lkm', 'Pipeline length', 'km', def, 0.05, 1500, group, 'geometry', 'Length of the pipeline section assessed; used for quantities.', { bind: 'pipeline.length' });
export const depthInput = (def = 120, group = 'Environment') => num('depth', 'Water depth', 'm', def, 0, 3500, group, 'environment', 'Still-water depth at the section; sets the external hydrostatic pressure.', { bind: 'site.depth' });
export const safetyClassInput = (group = 'Acceptance criteria') => sel('safetyClass', 'Safety class', 'medium', [['low', 'Low'], ['medium', 'Medium'], ['high', 'High']], group, 'criteria',
  'DNV-ST-F101 safety class of the section and phase (temporary phases are usually low; operation outside platform zones medium; near platforms high).');

/* ------------------------------------------------------------------------------------------- route */

export const DEFAULT_ROUTE = [
  { kp: 0, z: -320 }, { kp: 2, z: -335 }, { kp: 4, z: -310 }, { kp: 6, z: -325 }, { kp: 8, z: -300 }, { kp: 10, z: -305 },
];
/** Route as an imported 3-D polyline, with a chainage / elevation table as the parametric fallback. */
export function routeInputs({ group = 'Route', table = DEFAULT_ROUTE } = {}) {
  return [
    { key: 'route', label: 'Imported route centreline (optional)', type: 'geometry', kinds: ['polyline'], default: null, group, role: 'geometry',
      help: 'A 3-D polyline (x, y, z with z up, metres). Chainage is the horizontal distance along it and z the pipe elevation. Leave empty to use the table.' },
    { key: 'profile', label: 'Route profile', type: 'table', group, role: 'geometry',
      columns: [{ key: 'kp', label: 'Chainage (KP)', unit: 'km', type: 'number' }, { key: 'z', label: 'Pipe elevation (negative below sea level)', unit: 'm', type: 'number' }],
      default: table, help: 'Chainage from the inlet and elevation of the pipe axis; at least two rows, chainage increasing.' },
  ];
}
/** @returns {{x:number,z:number}[]} horizontal distance x [m] and elevation z [m], inlet first */
export function routeFrom(inp) {
  const g = inp.route;
  let pts = [];
  if (g && g.points && g.points.length >= 6) {
    const P = g.points; let s = 0;
    for (let i = 0; i < P.length / 3; i++) {
      if (i > 0) s += Math.hypot(P[3 * i] - P[3 * i - 3], P[3 * i + 1] - P[3 * i - 2]);
      pts.push({ x: s, z: P[3 * i + 2] });
    }
    // vertical legs (risers) have no horizontal run: give them a nominal one so the chainage stays increasing
    pts = pts.filter((p, i) => i === 0 || p.x > pts[i - 1].x + 1e-6 || Math.abs(p.z - pts[i - 1].z) > 1e-6);
    for (let i = 1; i < pts.length; i++) if (pts[i].x <= pts[i - 1].x) pts[i].x = pts[i - 1].x + 0.01 * Math.abs(pts[i].z - pts[i - 1].z) + 1e-3;
  } else {
    const rows = (inp.profile || []).map((r) => ({ x: Number(r.kp) * 1e3, z: Number(r.z) })).filter((r) => Number.isFinite(r.x) && Number.isFinite(r.z));
    pts = rows;
  }
  need(pts.length >= 2, 'The route needs at least two points (chainage and elevation).');
  for (let i = 1; i < pts.length; i++) need(pts[i].x > pts[i - 1].x, `Route chainage must increase: row ${i + 1} (${fmt(pts[i].x / 1e3)} km) does not follow row ${i} (${fmt(pts[i - 1].x / 1e3)} km).`);
  return pts;
}
/** Along-pipe length of a route [m]. */
export const routeLength = (pts) => pts.slice(1).reduce((s, p, i) => s + Math.hypot(p.x - pts[i].x, p.z - pts[i].z), 0);
/** Linear interpolation in a sorted table. */
export function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1; while (xs[i] < x) i++;
  const w = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
  return ys[i - 1] * (1 - w) + ys[i] * w;
}
