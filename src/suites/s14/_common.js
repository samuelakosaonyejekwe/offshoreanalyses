// Shared helpers of the flexible-pipe suite: the layer-stack table and its conversion to the engine stack, the
// common load inputs, result builders and the standing limitations of the cross-section models. Engineering units
// are converted to SI here so every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { buildFlexibleStack, massAndWeight } from '../../engine/xsection/flexible.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const BAR = 1e5;
// Screening emission factors (cradle-to-gate order of magnitude, not project EPD values).
export const CO2_STEEL = 1.9;     // t CO2e per t steel
export const CO2_MGO = 3.206;     // t CO2 per t marine gas oil
export const MGO_T_PER_DAY = 30;  // flex-lay / construction vessel fuel burn, t/day (screening)

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
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}

/* ------------------------------------------------------------------------------------------- standing text */

export const UNBONDED = 'Unbonded flexible pipe only: bonded (vulcanised) hose constructions are outside these models.';
export const LIMIT_AXI = 'Axisymmetric elastic layer model: carcass and pressure armour are helical strips of equivalent area (no interlock-profile contact mechanics), wires are radially rigid, there is no plasticity and no 3-D finite-element contact. ' + UNBONDED;
export const GENERIC = 'The default layer stack is a generic, indicative 6-inch riser — replace it with the manufacturer\'s pipe data sheet before using any result.';
export const CONFIRM = 'Default recalled from the standard — confirm against the edition governing the project.';
export const LIMIT_2D = 'The rod solver is two-dimensional (in-plane): out-of-plane current, cross offsets and torsion of the global line are not modelled.';
export const STD_17J = 'API Spec 17J (4th ed., 2014) §5.3 and Table 6';
export const STD_17B = 'API RP 17B (5th ed., 2014) §8';

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

/* ------------------------------------------------------------------------------------------- layer stack */

const TYPES = ['carcass', 'pressureSheath', 'pressureArmour', 'backupArmour', 'antiWear', 'tensileArmour', 'insulation', 'outerSheath', 'tape'];
const TYPE_ALIAS = Object.fromEntries(TYPES.map((t) => [t.toLowerCase(), t]));
Object.assign(TYPE_ALIAS, { 'pressure sheath': 'pressureSheath', 'pressure armour': 'pressureArmour', 'tensile armour': 'tensileArmour', 'anti-wear': 'antiWear', antiwear: 'antiWear', 'outer sheath': 'outerSheath', 'backup armour': 'backupArmour', liner: 'pressureSheath' });

// Generic 6-inch rough-bore dynamic riser (same construction as the engine demonstration design). Indicative only.
export const DEFAULT_LAYERS = [
  { name: 'Carcass', type: 'carcass', t: 7, alphaDeg: 87.5, lay: 1, fill: 0.55, wireW: 0, E: 193, yield: 500, uts: 750, lock: 9 },
  { name: 'Pressure sheath', type: 'pressureSheath', t: 8, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 },
  { name: 'Pressure armour', type: 'pressureArmour', t: 8, alphaDeg: 87, lay: 1, fill: 0.9, wireW: 0, E: 207, yield: 800, uts: 1000, lock: 7 },
  { name: 'Anti-wear tape 1', type: 'antiWear', t: 1.5, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 },
  { name: 'Tensile armour 1', type: 'tensileArmour', t: 4, alphaDeg: 35, lay: 1, fill: 0.9, wireW: 12, E: 207, yield: 1100, uts: 1400, lock: 0 },
  { name: 'Anti-wear tape 2', type: 'antiWear', t: 1.5, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 },
  { name: 'Tensile armour 2', type: 'tensileArmour', t: 4, alphaDeg: 35, lay: -1, fill: 0.9, wireW: 12, E: 207, yield: 1100, uts: 1400, lock: 0 },
  { name: 'Outer sheath', type: 'outerSheath', t: 8, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 },
];

const LAYER_COLUMNS = [
  { key: 'name', label: 'Layer', type: 'text' },
  { key: 'type', label: 'Type (carcass, pressureSheath, pressureArmour, backupArmour, antiWear, tensileArmour, insulation, outerSheath, tape)', type: 'text' },
  { key: 't', label: 'Thickness', unit: 'mm', type: 'number' },
  { key: 'alphaDeg', label: 'Lay angle from pipe axis (0 = polymer layer)', unit: 'deg', type: 'number' },
  { key: 'lay', label: 'Lay direction (+1 / −1)', type: 'number' },
  { key: 'fill', label: 'Metal fill fraction', unit: '-', type: 'number' },
  { key: 'wireW', label: 'Wire width (0 = unknown)', unit: 'mm', type: 'number' },
  { key: 'E', label: 'Young\'s modulus (0 = from material library)', unit: 'GPa', type: 'number' },
  { key: 'yield', label: 'Yield strength', unit: 'MPa', type: 'number' },
  { key: 'uts', label: 'Tensile strength', unit: 'MPa', type: 'number' },
  { key: 'lock', label: 'Locking strain of interlocked profile', unit: '%', type: 'number' },
];

/** Bore, layer table, polymer materials and optional imported cross-section shared by every module. */
export function stackInputs() {
  const group = 'Pipe construction';
  return [
    { key: 'boreD', label: 'Bore (inner) diameter', unit: 'mm', type: 'number', default: 152.4, min: 25, max: 600, group, role: 'geometry',
      help: 'Inside diameter of the innermost layer (6 in = 152.4 mm). ' + GENERIC },
    { key: 'layers', label: 'Layer stack, bore outwards', type: 'table', columns: LAYER_COLUMNS, default: DEFAULT_LAYERS, group, role: 'geometry',
      help: 'One row per layer from the bore outwards. Metallic layers carry a lay angle, fill fraction and strengths; polymer layers (lay angle 0) take their properties from the material selections below unless a modulus is typed. Carcass and pressure armour are modelled as helical strips of equivalent area — the interlocked profile itself is not resolved. Unbonded construction only.' },
    { key: 'matSheath', label: 'Pressure sheath / tape polymer', type: 'select', source: 'materials:polymer', default: 'pa11', group, role: 'material',
      help: 'Supplies modulus, Poisson\'s ratio, density, conductivity and allowable strain of the pressure sheath and anti-wear tapes.' },
    { key: 'matOuter', label: 'Outer sheath polymer', type: 'select', source: 'materials:polymer', default: 'hdpe', group, role: 'material',
      help: 'Supplies the properties of the outer sheath.' },
    { key: 'xsec', label: 'Imported cross-section (optional)', type: 'geometry', kinds: ['section'], default: null, group, role: 'geometry',
      help: 'A section model whose "layers" list (same fields as the table) replaces the table, and whose innermost loop sets the bore. Leave empty to use the form fields.' },
  ];
}

const POLY_ALPHA = { pa11: 1.1e-4, pa12: 1.1e-4, hdpe: 1.8e-4, pvdf: 1.4e-4 };
const POLY_DEFAULT = { E: 0.35e9, nu: 0.42, rho: 1040, k: 0.25, cp: 2000, fy: 25e6, strainLimit: 0.077 };

function polymerProps(id, ctx) {
  const m = (ctx?.material && ctx.material(id)) || null;
  const p = { ...POLY_DEFAULT, ...(m || {}) };
  return {
    E: p.E || POLY_DEFAULT.E, nu: Math.min(p.nu ?? 0.42, 0.45), rho: p.rho || 1000, k: p.k || 0.25, cp: p.cp || 2000, fy: p.fy || 25e6,
    alphaT: POLY_ALPHA[id] ?? 1.5e-4, strainLimit: p.strainLimit || 0.077, strainLimitDynamic: id === 'pvdf' ? 0.035 : (p.strainLimit || 0.077), tMax: p.tMax ?? (id === 'hdpe' ? 60 : 90), id,
  };
}

/**
 * Engine layer stack from the form values.
 * Metallic layers: carcass taken as austenitic stainless (density 7980 kg/m³, k 15 W/m/K, α 16e-6 /K), the other
 * helical layers as carbon steel (7850, 45, 11.7e-6). Polymer layers from the materials library.
 * @returns {{ stack, defs, idx: { carcass, sheath, pa, ta: number[], outer }, od, bore, sheathMat, outerMat }}
 */
export function stackFrom(inp, ctx, { polymerEFactor = 1, strain } = {}) {
  let rows = inp.layers;
  let bore = inp.boreD / 1e3;
  const sec = inp.xsec;
  if (sec && Array.isArray(sec.layers) && sec.layers.length) rows = sec.layers;
  if (sec && Array.isArray(sec.loops) && sec.loops[0]?.length >= 6) {
    // innermost loop: mean distance of its vertices from their centroid is the bore radius
    let rMin = Infinity;
    for (const lp of sec.loops) {
      const n = lp.length / 2; let cx = 0, cy = 0;
      for (let i = 0; i < n; i++) { cx += lp[2 * i]; cy += lp[2 * i + 1]; }
      cx /= n; cy /= n; let r = 0;
      for (let i = 0; i < n; i++) r += Math.hypot(lp[2 * i] - cx, lp[2 * i + 1] - cy);
      rMin = Math.min(rMin, r / n);
    }
    if (Number.isFinite(rMin) && rMin > 0) bore = 2 * rMin;
  }
  need(Array.isArray(rows) && rows.length > 0, 'The layer stack is empty: enter at least a pressure sheath.');
  const sheathMat = polymerProps(inp.matSheath, ctx), outerMat = polymerProps(inp.matOuter, ctx);
  const tapeMat = polymerProps('pa11', ctx), foam = { ...POLY_DEFAULT, E: 0.2e9, nu: 0.4, rho: 700, k: 0.17, cp: 2000, alphaT: 1.5e-4 };
  const defs = rows.map((r, i) => {
    const type = TYPE_ALIAS[String(r.type ?? '').trim().toLowerCase()];
    need(type, `Layer ${i + 1} ('${r.name ?? ''}'): type '${r.type}' is not recognised. Use one of ${TYPES.join(', ')}.`);
    const t = Number(r.t ?? (r.thickness != null ? r.thickness * 1e3 : NaN)) / 1e3;
    need(t > 0, `Layer ${i + 1} ('${r.name ?? type}'): thickness must be positive.`);
    const name = String(r.name || `${type} ${i + 1}`);
    const a = Number(r.alphaDeg) || 0;
    if (a > 0) {
      need(a < 90, `${name}: lay angle must be below 90°.`);
      need(!['pressureSheath', 'outerSheath'].includes(type), `${name}: a sealing sheath cannot have a lay angle — set it to 0.`);
      const fill = Number(r.fill);
      need(fill > 0 && fill <= 1, `${name}: fill fraction must be in (0, 1].`);
      const ss = type === 'carcass';
      const L = { name, type, thickness: t, alphaDeg: a, lay: Number(r.lay) < 0 ? -1 : 1, fill, E: (Number(r.E) || (ss ? 193 : 207)) * 1e9, nu: 0.3,
        rho: ss ? 7980 : 7850, k: ss ? 15 : 45, cp: 500, alphaT: ss ? 1.6e-5 : 1.17e-5 };
      if (Number(r.yield) > 0) L.yield = r.yield * 1e6;
      if (Number(r.uts) > 0) L.uts = r.uts * 1e6;
      if (Number(r.lock) > 0) L.lockStrain = r.lock / 100;
      if (Number(r.wireW) > 0) { L.wireWidth = r.wireW / 1e3; L.wireThickness = t; }
      return L;
    }
    need(!['carcass', 'pressureArmour', 'backupArmour', 'tensileArmour'].includes(type), `${name}: a metallic ${type} layer needs a lay angle between 0 and 90°.`);
    const m = type === 'pressureSheath' ? sheathMat : type === 'outerSheath' ? outerMat : type === 'insulation' ? foam : tapeMat;
    const L = { name, type, thickness: t, E: (Number(r.E) > 0 ? r.E * 1e9 : m.E) * polymerEFactor, nu: m.nu, rho: m.rho, k: m.k, cp: m.cp, alphaT: m.alphaT, fy: m.fy };
    if (type === 'pressureSheath' || type === 'outerSheath') {
      L.strainLimit = strain?.static ?? m.strainLimit; L.strainLimitDynamic = strain?.dynamic ?? m.strainLimitDynamic;
    }
    return L;
  });
  need(defs.some((d) => d.type === 'pressureSheath' || d.type === 'outerSheath'), 'The stack needs at least one sealing layer (pressureSheath or outerSheath).');
  const stack = run(() => buildFlexibleStack({ innerDiameter: bore, layers: defs }), 'Layer stack');
  const find = (ty) => stack.layers.findIndex((L) => L.type === ty);
  const idx = { carcass: find('carcass'), sheath: find('pressureSheath'), pa: find('pressureArmour'), outer: stack.sealOuter,
    ta: stack.layers.filter((L) => L.type === 'tensileArmour').map((L) => L.index) };
  return { stack, defs, idx, od: 2 * stack.outerRadius, bore, sheathMat, outerMat };
}

/** Steel and polymer mass per metre of pipe [kg/m]. */
export function massSplit(stack) {
  let steel = 0, polymer = 0;
  for (const L of stack.layers) {
    if (L.kind === 'helical') steel += L.rho * L.nA / Math.cos(L.alpha);
    else polymer += L.rho * Math.PI * (L.rOut ** 2 - L.rIn ** 2);
  }
  return { steel, polymer };
}
export const weights = (stack, inp, flooded = false) => run(() => massAndWeight(stack, { rhoContent: inp.rhoContent ?? 0, rhoSeawater: inp.rhoW ?? 1025, annulusFlooded: flooded }), 'Mass and weight');
export const steelCo2 = (t) => `${fmt(Math.abs(t))} t of armour and carcass steel ≈ ${fmt(Math.abs(t) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (replace with the manufacturer's product declaration).`;
export const vesselCo2 = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

/* ------------------------------------------------------------------------------------------- common loads */

export function loadInputs({ Te = 500, pInt = 200, depth = 1000, torque = false, group = 'Loads and pressure boundaries' } = {}) {
  const list = [
    { key: 'Te', label: 'Effective tension', unit: 'kN', type: 'number', default: Te, min: -2000, max: 20000, group, role: 'bc', bind: 'flexible.topTension',
      uncertain: { dist: 'normal', cov: 0.1 }, help: 'Effective tension at the section (true wall force minus the pressure end-cap terms). Negative = effective compression.' },
    { key: 'pInt', label: 'Internal (bore) pressure', unit: 'bar', type: 'number', default: pInt, min: 0, max: 1500, group, role: 'bc',
      uncertain: { dist: 'normal', cov: 0.05 }, help: 'Absolute bore pressure acting on the inner face of the pressure sheath (the carcass is not leak-tight).' },
    { key: 'depth', label: 'Water depth at the section', unit: 'm', type: 'number', default: depth, min: 0, max: 3500, group, role: 'bc', bind: 'site.depth',
      help: 'Sets the external hydrostatic pressure ρ·g·depth on the outer sheath.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group, role: 'environment', bind: 'site.rhoW' },
    { key: 'annulus', label: 'Annulus condition', type: 'select', default: 'dry', options: [{ value: 'dry', label: 'Dry (intact outer sheath, vented)' }, { value: 'flooded', label: 'Flooded with sea water' }], group, role: 'ic',
      help: 'Dry: the annulus is at the vent pressure and external pressure is carried through the armour. Flooded: the annulus is at external pressure and the carcass alone resists collapse.' },
    { key: 'pAnn', label: 'Annulus gas pressure (dry)', unit: 'bar', type: 'number', default: 1, min: 0, max: 50, group, role: 'bc', showIf: { key: 'annulus', in: ['dry'] },
      help: 'Absolute pressure of the vented annulus; a few bar at most with a working vent system.' },
  ];
  if (torque) list.push({ key: 'torque', label: 'Applied torque', unit: 'kN·m', type: 'number', default: 0, min: -500, max: 500, group, role: 'bc',
    help: 'Torque about the pipe axis. Positive tightens layers of lay +1.' });
  return list;
}
/** Engine load record (SI) from the common load inputs. */
export function loadsFrom(inp) {
  const pExt = (inp.rhoW ?? 1025) * G * (inp.depth ?? 0);
  const o = { tension: (inp.Te ?? 0) * 1e3, pInt: (inp.pInt ?? 0) * BAR, pExt, annulus: inp.annulus ?? 'dry', pAnnulus: (inp.pAnn ?? 0) * BAR };
  if (inp.torque != null) o.torque = inp.torque * 1e3;
  return o;
}

/** API Spec 17J Table 6 utilisation factors as overridable criteria. */
export function factorInputs({ group = 'Acceptance criteria (API 17J Table 6)' } = {}) {
  const f = (key, label, def) => ({ key, label, unit: '-', type: 'number', default: def, min: 0.3, max: 1, group, role: 'criteria', help: CONFIRM });
  return [
    { key: 'condition', label: 'Load condition', type: 'select', default: 'normal', options: [{ value: 'normal', label: 'Normal operation (recurrent)' }, { value: 'abnormal', label: 'Abnormal / extreme operation' }, { value: 'hydrotest', label: 'Hydrostatic test (FAT)' }], group, role: 'criteria',
      help: 'Selects which permissible utilisation applies to the armour layers.' },
    f('fTaN', 'Tensile armour factor — normal', 0.67), f('fTaA', 'Tensile armour factor — abnormal', 0.85), f('fTaH', 'Tensile armour factor — hydrotest', 0.91),
    f('fPaN', 'Pressure armour factor — normal', 0.55), f('fPaA', 'Pressure armour factor — abnormal', 0.85), f('fPaH', 'Pressure armour factor — hydrotest', 0.91),
  ];
}
export const factorsFrom = (inp) => ({
  tensileArmour: { normal: inp.fTaN, abnormal: inp.fTaA, hydrotest: inp.fTaH },
  pressureArmour: { normal: inp.fPaN, abnormal: inp.fPaA, hydrotest: inp.fPaH },
});
export const muInput = (group = 'Interlayer contact') => ({ key: 'mu', label: 'Interlayer friction coefficient', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group, role: 'material',
  uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Steel wire on anti-wear tape: typically 0.05–0.2, higher for steel on steel or a dry, aged annulus. No normative value — calibrate against bending-hysteresis tests.' });

/** Structural capacity of a layer as API 17J defines it: 0.9·UTS when known, else yield [Pa]. */
export const capacityOf = (L) => (L.uts ? 0.9 * L.uts : L.yield);
