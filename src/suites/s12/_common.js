// Shared helpers of the subsea-umbilical suite: unit handling, the cross-section definition (component and layer
// tables → engine assembler), the global line properties and configuration (free-hanging catenary or lazy wave on
// the 2-D rod solver), acceptance-criteria inputs, screening sustainability factors and advice records.
// Engineering units are converted to SI here so that every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { assembleUmbilical } from '../../engine/xsection/umbilical.js';
import { solveRiser, solveWaveRiser } from '../../engine/lines/riser.js';
import { buoyancyModules } from '../../engine/xsection/ancillaries.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 20;   // typical umbilical lay / construction vessel fuel burn, t/day (screening)
export const CO2_POLYMER = 2.5;    // t CO2e per t thermoplastic / polyurethane (screening, typical PE–PU range 2–4)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible umbilical.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
export const scaled = (a, f) => Array.from(a, (v) => v * f);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

/* ------------------------------------------------------------------------------------------------ advice */

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? uc : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc, warnAt = 0.9) => (uc > 1 ? 'fail' : uc > warnAt ? 'warn' : 'ok');

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (super-duplex tube is several times higher per tonne — confirm with the mill's product declaration).`;
export const vesselCo2Text = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const LIMIT_2D = 'The rod model is two-dimensional (in-plane): out-of-plane current, cross offsets, torsion and lateral seabed friction are not modelled, and the seabed is a row of linear lift-off springs.';
export const LIMIT_RIGID = 'The cross-section assembler treats the bundle as radially rigid, so the axial stiffness is an upper estimate for filler-rich bundles (calibrate against a measured tensile test); generic component data is indicative only.';
export const CONFIRM = 'Default in the ISO 13628-5 / API 17E pattern — confirm against the governing edition and the purchaser specification.';
export const STD_UMB = 'ISO 13628-5:2009 / API Spec 17E (5th ed.)';

/* --------------------------------------------------------------------------------- cross-section definition */

const KINDS = ['tube', 'hose', 'cable', 'fibre', 'filler', 'armour', 'strength'];

// Generic 145 mm steel-tube dynamic umbilical with two contra-wound armour layers: demonstration data, not a
// manufacturer design. Cables and fibre units are entered as solid rods with an effective modulus (measured EA / area).
export const DEFAULT_COMPONENTS = [
  { name: 'Centre tube 1 in (methanol)', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 345, rhoF: 792 },
  { name: 'HP hydraulic tube 1/2 in', kind: 'tube', layer: 1, count: 5, od: 16.5, wall: 1.9, E: 200, rho: 7800, fy: 550, p: 690, rhoF: 1060 },
  { name: 'Filler (inner layer)', kind: 'filler', layer: 1, count: 3, od: 16.5, wall: 0, E: 0.8, rho: 950, fy: 0, p: 0, rhoF: 0 },
  { name: 'LP hydraulic tube 5/8 in', kind: 'tube', layer: 2, count: 4, od: 20, wall: 1.8, E: 200, rho: 7800, fy: 550, p: 345, rhoF: 1060 },
  { name: 'Electrical quad cable', kind: 'cable', layer: 2, count: 3, od: 20, wall: 0, E: 12, rho: 2600, fy: 0, p: 0, rhoF: 0 },
  { name: 'Fibre-optic unit', kind: 'fibre', layer: 2, count: 1, od: 20, wall: 0, E: 5, rho: 1800, fy: 0, p: 0, rhoF: 0 },
  { name: 'Filler (outer layer)', kind: 'filler', layer: 2, count: 5, od: 20, wall: 0, E: 0.8, rho: 950, fy: 0, p: 0, rhoF: 0 },
  { name: 'Armour wire (inner layer)', kind: 'armour', layer: 3, count: 68, od: 5, wall: 0, E: 200, rho: 7850, fy: 1100, p: 0, rhoF: 0 },
  { name: 'Armour wire (outer layer)', kind: 'armour', layer: 4, count: 74, od: 5, wall: 0, E: 200, rho: 7850, fy: 1100, p: 0, rhoF: 0 },
];
export const DEFAULT_LAYERS = [
  { layer: 1, R: 25.5, alpha: 8, lay: 1 },
  { layer: 2, R: 45, alpha: 10, lay: -1 },
  { layer: 3, R: 59, alpha: 15, lay: 1 },
  { layer: 4, R: 64, alpha: 11, lay: -1 },   // shorter lay than the inner armour to balance torque (R·sinα·cos²α weighting)
];

/** Cross-section definition inputs: component table, layer table, outer sheath, optional imported section. */
export function sectionInputs({ pressures = true } = {}) {
  return [
    { key: 'sectionGeom', label: 'Imported cross-section outline (optional)', type: 'geometry', kinds: ['section'], default: null, group: 'Cross-section', role: 'geometry',
      help: 'DXF / SVG / JSON section. Only the outer boundary is used (overall diameter); the internal lay-up always comes from the tables below, which are the parametric fallback.' },
    { key: 'components', label: 'Internal components', type: 'table', group: 'Cross-section', role: 'geometry',
      columns: [
        { key: 'name', label: 'Component', type: 'text' },
        { key: 'kind', label: 'Kind (tube / hose / cable / fibre / filler / armour / strength)', type: 'text' },
        { key: 'layer', label: 'Layer (0 = centre)', type: 'int' },
        { key: 'count', label: 'Number', type: 'int' },
        { key: 'od', label: 'Outside diameter', unit: 'mm', type: 'number' },
        { key: 'wall', label: 'Wall (0 = solid)', unit: 'mm', type: 'number' },
        { key: 'E', label: 'Modulus (effective)', unit: 'GPa', type: 'number' },
        { key: 'rho', label: 'Density', unit: 'kg/m³', type: 'number' },
        { key: 'fy', label: 'Yield (0 = not stress-checked)', unit: 'MPa', type: 'number' },
        ...(pressures ? [{ key: 'p', label: 'Internal pressure', unit: 'bar', type: 'number' }, { key: 'rhoF', label: 'Fluid density', unit: 'kg/m³', type: 'number' }] : []),
      ],
      default: DEFAULT_COMPONENTS.map((r) => ({ ...r })),
      help: 'One row per component type. Cables, hoses and fibre units are solid rods with an effective modulus = measured EA / gross area. Generic demonstration data — replace with the manufacturer data sheet.' },
    { key: 'layers', label: 'Helical layers', type: 'table', group: 'Cross-section', role: 'geometry',
      columns: [
        { key: 'layer', label: 'Layer', type: 'int' },
        { key: 'R', label: 'Pitch-circle radius', unit: 'mm', type: 'number' },
        { key: 'alpha', label: 'Lay angle from the axis', unit: 'deg', type: 'number' },
        { key: 'lay', label: 'Hand (+1 / −1)', type: 'int' },
      ],
      default: DEFAULT_LAYERS.map((r) => ({ ...r })),
      help: 'Contra-wound layers (alternating hand) balance torque. Lay length = 2πR / tan(lay angle).' },
    { key: 'sheathT', label: 'Outer sheath thickness', unit: 'mm', type: 'number', default: 6, min: 0, max: 30, group: 'Cross-section', role: 'geometry', design: true,
      help: 'Extruded over the outermost layer. 0 = no sheath (roving only).' },
    { key: 'sheathMat', label: 'Outer sheath material', type: 'select', source: 'materials:polymer', default: 'hdpe', group: 'Cross-section', role: 'material' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
  ];
}

/** Overall diameter from an imported section model (bounding circle of all loops) or null. */
function importedDiameter(geom) {
  if (!geom || !Array.isArray(geom.loops) || !geom.loops.length) return null;
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const loop of geom.loops) for (let i = 0; i + 1 < loop.length; i += 2) {
    xmin = Math.min(xmin, loop[i]); xmax = Math.max(xmax, loop[i]); ymin = Math.min(ymin, loop[i + 1]); ymax = Math.max(ymax, loop[i + 1]);
  }
  const d = Math.max(xmax - xmin, ymax - ymin);
  return Number.isFinite(d) && d > 0 ? d : null;
}

/**
 * Assemble the cross-section from the form tables (SI out).
 * Returns { sec (engine section), def, OD, rows (components in SI), sheath, tubes (stress-checked tubes) }.
 */
export function buildSection(inp, ctx, { pressureScale = 1 } = {}) {
  const comps = inp.components || [];
  need(comps.length > 0, 'The component table is empty: add at least one tube, cable or wire.');
  const layerRows = new Map((inp.layers || []).map((l) => [Math.round(l.layer), l]));
  const mk = (r, i) => {
    const kind = String(r.kind || 'filler').trim().toLowerCase();
    need(KINDS.includes(kind), `Component row ${i + 1}: kind '${r.kind}' is not one of ${KINDS.join(', ')}.`);
    need(r.od > 0, `Component '${r.name}': outside diameter must be positive.`);
    need(r.wall >= 0 && r.wall < r.od / 2, `Component '${r.name}': wall ${r.wall} mm must be smaller than the radius ${r.od / 2} mm (0 = solid).`);
    need(r.count >= 1, `Component '${r.name}': number must be at least 1.`);
    need(r.E >= 0 && r.rho >= 0, `Component '${r.name}': modulus and density cannot be negative.`);
    const od = r.od / 1e3, id = r.wall > 0 ? (r.od - 2 * r.wall) / 1e3 : 0;
    return { name: r.name || `component ${i + 1}`, kind, od, id, E: r.E * 1e9, rho: r.rho, yieldStress: r.fy > 0 ? r.fy * 1e6 : undefined,
      pressure: (r.p || 0) * 1e5 * pressureScale, fluidDensity: r.rhoF || 0 };
  };
  const rows = comps.map((r, i) => ({ ...mk(r, i), layer: Math.round(r.layer), count: Math.round(r.count) }));
  const centre = rows.filter((r) => r.layer === 0);
  need(centre.length <= 1 && (centre[0]?.count ?? 1) === 1, 'Layer 0 is the straight centre element: give at most one component with number = 1.');
  const ids = [...new Set(rows.filter((r) => r.layer > 0).map((r) => r.layer))].sort((a, b) => a - b);
  need(ids.length > 0 || centre.length === 1, 'No helical layer and no centre element are defined.');
  let prevR = centre.length ? centre[0].od / 2 : 0;
  const layers = ids.map((id) => {
    const L = layerRows.get(id);
    need(L, `Layer ${id} is used in the component table but has no row in the helical-layer table (under-defined lay-up).`);
    need(L.R > 0 && L.alpha > 0 && L.alpha < 60, `Layer ${id}: pitch radius must be positive and the lay angle between 0 and 60°.`);
    need(L.lay === 1 || L.lay === -1, `Layer ${id}: hand must be +1 or −1.`);
    const els = rows.filter((r) => r.layer === id);
    const R = L.R / 1e3, dMax = Math.max(...els.map((e) => e.od));
    need(R - dMax / 2 >= prevR - 1e-6, `Layer ${id}: pitch radius ${L.R} mm puts its elements inside the layer below (needs ≥ ${fmt((prevR + dMax / 2) * 1e3, 4)} mm).`);
    // elements must fit round the pitch circle: Σ n·d / cos α ≤ 2πR
    const occupied = els.reduce((a, e) => a + e.count * e.od / Math.cos(L.alpha * DEG), 0);
    need(occupied <= 2 * Math.PI * R * 1.001, `Layer ${id}: ${fmt(occupied * 1e3, 4)} mm of elements do not fit on a pitch circle of ${fmt(2 * Math.PI * R * 1e3, 4)} mm circumference (over-filled layer).`);
    prevR = R + dMax / 2;
    return { R, alphaDeg: L.alpha, lay: L.lay, elements: els.map((e) => ({ component: e, count: e.count })) };
  });
  const sm = ctx.material(inp.sheathMat);
  need(sm, `Sheath material '${inp.sheathMat}' is not in the materials library.`);
  const tS = (inp.sheathT || 0) / 1e3;
  const sheath = tS > 0 ? { ri: prevR, ro: prevR + tS, E: sm.E ?? 0.8e9, nu: sm.nu ?? 0.4, rho: sm.rho ?? 950, k: sm.k ?? 0.4 } : null;
  const def = { centre: centre[0], layers, sheaths: sheath ? [sheath] : [] };
  const imp = importedDiameter(inp.sectionGeom);
  if (imp) { need(imp >= 2 * prevR, `The imported outline (${fmt(imp * 1e3, 4)} mm) is smaller than the lay-up envelope (${fmt(2 * prevR * 1e3, 4)} mm): check the drawing scale.`); def.outerDiameter = imp; }
  const sec = run(() => assembleUmbilical(def, { rhoSeawater: inp.rhoW ?? 1025 }), 'The cross-section could not be assembled');
  const tubes = sec.elements.filter((e) => e.component.id > 0 && e.component.yieldStress);
  return { sec, def, OD: sec.outerDiameter, rows, sheath, sheathMat: sm, tubes, layers };
}

/** Steel and polymer mass per metre of the assembled section [kg/m]. */
export function sectionMasses(built) {
  let steel = 0, polymer = 0, other = 0;
  for (const e of built.sec.elements) {
    const m = e.count * e.component.mass / Math.cos(e.alpha);
    if (['tube', 'armour', 'strength'].includes(e.kind)) steel += m; else if (['filler', 'hose'].includes(e.kind)) polymer += m; else other += m;
  }
  if (built.sheath) polymer += built.sheath.rho * Math.PI * (built.sheath.ro ** 2 - built.sheath.ri ** 2);
  return { steel, polymer, other };
}

/** Acceptance-criteria inputs shared by the capacity modules (utilisation factors by load condition). */
export function criteriaInputs({ condition = 'normal' } = {}) {
  return [
    { key: 'condition', label: 'Load condition', type: 'select', default: condition, group: 'Acceptance criteria', role: 'criteria',
      options: [{ value: 'normal', label: 'Normal operation' }, { value: 'abnormal', label: 'Abnormal operation / installation' }, { value: 'accidental', label: 'Accidental / survival' }] },
    { key: 'etaNormal', label: 'Utilisation factor — normal', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Allowable von Mises stress / yield. ${CONFIRM}` },
    { key: 'etaAbnormal', label: 'Utilisation factor — abnormal', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'etaAccidental', label: 'Utilisation factor — accidental', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'muInt', label: 'Friction coefficient between elements', unit: '-', type: 'number', default: 0.2, min: 0, max: 0.6, group: 'Acceptance criteria', role: 'material',
      uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'No normative value: 0.1–0.3 for steel on polymer. Calibrate against cyclic bending tests.' },
  ];
}
export const etaOf = (inp) => ({ normal: inp.etaNormal, abnormal: inp.etaAbnormal, accidental: inp.etaAccidental })[inp.condition];

/* ------------------------------------------------------------------------------------ global line properties */

/** Equivalent-beam properties of the umbilical for the global (rod) analyses; defaults are close to the generic section (the bus values replace them once the cross-section module has run). */
export function lineInputs({ group = 'Umbilical (global properties)' } = {}) {
  return [
    { key: 'OD', label: 'Overall diameter', unit: 'mm', type: 'number', default: 145, min: 30, max: 400, group, role: 'geometry', bind: 'umbilical.OD',
      help: 'Outside diameter over the sheath (buoyancy and drag).' },
    { key: 'mass', label: 'Mass in air, tubes filled', unit: 'kg/m', type: 'number', default: 40.9, min: 1, max: 300, group, role: 'geometry', bind: 'umbilical.massPerM',
      uncertain: { dist: 'normal', cov: 0.03 }, help: 'From the cross-section module or the manufacturer data sheet.' },
    { key: 'EA', label: 'Axial stiffness EA', unit: 'MN', type: 'number', default: 720, min: 1, max: 5000, group, role: 'material', bind: 'umbilical.EA',
      uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'The assembler value is an upper estimate (radially rigid bundle); a measured value from the tensile test is preferred.' },
    { key: 'EI', label: 'Bending stiffness EI', unit: 'N·m²', type: 'number', default: 17100, min: 100, max: 5e6, group, role: 'material', bind: 'umbilical.EI', design: true,
      uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Slip (large-curvature) bending stiffness for global shape; the stick value applies only below the slip curvature.' },
    { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Hydrodynamic coefficients', role: 'load',
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'Bare sheathed umbilical 1.0–1.2; raise it when VIV amplifies drag. Calibrate from tank tests where available.' },
    { key: 'Cdt', label: 'Tangential drag coefficient', unit: '-', type: 'number', default: 0.01, min: 0, max: 0.2, group: 'Hydrodynamic coefficients', role: 'load' },
    { key: 'Ca', label: 'Added-mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 2.5, group: 'Hydrodynamic coefficients', role: 'load' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
  ];
}

/** Rod-section record (SI) from the global line inputs. */
export function lineFrom(inp) {
  const D = inp.OD / 1e3, rhoW = inp.rhoW ?? 1025;
  const Ae = Math.PI * D * D / 4;
  const w = (inp.mass - rhoW * Ae) * G;
  const ma = (inp.Ca ?? 1) * rhoW * Ae;
  return { D, Ae, w, m: inp.mass, ma, EA: inp.EA * 1e6, EI: inp.EI, rhoW,
    section: { EA: inp.EA * 1e6, EI: inp.EI, w, D, Cdn: inp.Cdn ?? 1.2, Cdt: inp.Cdt ?? 0.01, m: inp.mass, ma } };
}

/** Configuration inputs: free-hanging catenary or lazy wave between a floater hang-off and a seabed termination. */
export function configInputs({ config = 'lazy-wave' } = {}) {
  const wave = { key: 'config', in: ['lazy-wave'] };
  return [
    { key: 'config', label: 'Configuration', type: 'select', default: config, group: 'Configuration', role: 'ic',
      options: [{ value: 'lazy-wave', label: 'Lazy wave (distributed buoyancy)' }, { value: 'catenary', label: 'Free-hanging catenary' }],
      help: 'Installed geometry that the analysis starts from.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 350, min: 30, max: 3000, group: 'Configuration', role: 'environment', bind: 'site.depth' },
    { key: 'hangoffDepth', label: 'Hang-off depth below still water', unit: 'm', type: 'number', default: 10, min: -40, max: 200, group: 'Configuration', role: 'bc',
      help: 'Elevation of the I-tube exit / hang-off; negative = above the water line.' },
    { key: 'length', label: 'Total length hang-off to termination', unit: 'm', type: 'number', default: 740, min: 40, max: 12000, group: 'Configuration', role: 'geometry', design: true },
    { key: 'span', label: 'Horizontal distance hang-off to termination', unit: 'm', type: 'number', default: 480, min: 5, max: 11000, group: 'Configuration', role: 'bc', design: true,
      help: 'Position of the subsea termination (UTA / SDU) relative to the hang-off. The termination end is pinned on the seabed.' },
    { key: 'upperLength', label: 'Upper catenary length (hang-off to first buoyancy module)', unit: 'm', type: 'number', default: 370, min: 10, max: 6000, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'buoyLength', label: 'Buoyancy-section length', unit: 'm', type: 'number', default: 130, min: 5, max: 3000, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'modOD', label: 'Buoyancy-module outside diameter', unit: 'm', type: 'number', default: 0.55, min: 0.1, max: 3, group: 'Buoyancy section', role: 'geometry', showIf: wave },
    { key: 'modL', label: 'Buoyancy-module length', unit: 'm', type: 'number', default: 0.9, min: 0.2, max: 4, group: 'Buoyancy section', role: 'geometry', showIf: wave },
    { key: 'modPitch', label: 'Buoyancy-module pitch', unit: 'm', type: 'number', default: 3, min: 0.3, max: 30, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'foamRho', label: 'Syntactic-foam density', unit: 'kg/m³', type: 'number', default: 450, min: 250, max: 800, group: 'Buoyancy section', role: 'material', showIf: wave,
      uncertain: { dist: 'normal', cov: 0.04 }, help: 'Rises with depth rating (≈ 350 at 500 m to ≈ 600 at 3000 m, indicative). Water absorption over life reduces uplift.' },
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 200, min: 5, max: 1e5, group: 'Seabed', role: 'soil',
      uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Force per metre of umbilical per metre of penetration (kN/m/m). Soft clay 20–200, sand 1000+.' },
    { key: 'muAx', label: 'Axial seabed friction coefficient', unit: '-', type: 'number', default: 0.4, min: 0, max: 1.2, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 } },
  ];
}

/** Current profile inputs: surface and near-bed speed with a power law in between (in-plane). */
export function currentInputs({ Us = 0.8, Ub = 0.15, bind = 'site.Uc100', group = 'Current' } = {}) {
  return [
    { key: 'Us', label: 'Surface current speed', unit: 'm/s', type: 'number', default: Us, min: -3, max: 3, group, role: 'bc', bind,
      uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'In-plane current at the surface, positive from the hang-off towards the termination (far direction); negative = near.' },
    { key: 'Ub', label: 'Near-bed current speed', unit: 'm/s', type: 'number', default: Ub, min: -2, max: 2, group, role: 'bc', bind: 'site.UcBed100' },
    { key: 'curExp', label: 'Profile exponent', unit: '-', type: 'number', default: 1.5, min: 0.1, max: 8, group, role: 'bc',
      help: 'U(z) = Ub + (Us − Ub)·(height/depth)^n. n = 1 linear shear; larger n concentrates the current near the surface.' },
  ];
}
/** U as a function of elevation z above the seabed [m] → m/s. */
export function currentFn(inp, depth) {
  const { Us, Ub, curExp } = inp;
  return (z) => Ub + (Us - Ub) * Math.min(1, Math.max(0, z / depth)) ** curExp;
}

/**
 * Build and solve the static configuration with the 2-D rod solver. Seabed is z = 0; the hang-off is at
 * (offset, depth − hangoffDepth) and the termination at (span, 0), pinned.
 * Returns { res, line, h, buoy (engine buoyancy record or null), buoySec, s0/s1 of the buoyancy section }.
 */
export function solveConfig(inp, { nElem = 160, current, offset = 0, heave = 0, initial } = {}) {
  const line = lineFrom(inp);
  const h = inp.depth - inp.hangoffDepth;
  need(h > 10, 'Hang-off depth must leave at least 10 m between the hang-off and the seabed.');
  need(line.w > 0, `The umbilical is buoyant in sea water (apparent weight ${fmt(line.w)} N/m): it cannot hang as a catenary. Increase mass or reduce diameter.`);
  const span = inp.span - offset, top = { x: 0, z: h + heave }, bottom = { x: span, z: 0 };
  need(span > 1, 'The vessel offset moves the hang-off beyond the termination.');
  need(inp.length > Math.hypot(span, h) * 1.002, `Under-constrained geometry: the length (${fmt(inp.length, 5)} m) is shorter than the straight distance hang-off to termination (${fmt(Math.hypot(span, h), 5)} m).`);
  need(inp.length < span + h + (inp.config === 'lazy-wave' ? 3 * inp.buoyLength : 0) + 0.6 * span, `Over-length: ${fmt(inp.length, 5)} m cannot be laid out between a hang-off ${fmt(h, 4)} m above the seabed and a termination ${fmt(span, 4)} m away without piling on the seabed. Shorten the line or move the termination.`);
  const seabed = { z: 0, kn: inp.kn * 1e3, mu: inp.muAx };
  let res, buoy = null, buoySec = null, s0 = null, s1 = null;
  if (inp.config === 'lazy-wave') {
    const lower = inp.length - inp.upperLength - inp.buoyLength;
    need(lower > 5, `Over-constrained lengths: upper (${inp.upperLength} m) + buoyancy section (${inp.buoyLength} m) leave no lower catenary within the total length of ${inp.length} m.`);
    need(inp.modOD > line.D * 1.2, 'The buoyancy-module diameter must exceed the umbilical diameter by at least 20 %.');
    need(inp.modPitch >= inp.modL, 'Buoyancy-module pitch cannot be smaller than the module length.');
    buoy = run(() => buoyancyModules({ outerDiameter: inp.modOD, innerDiameter: line.D, moduleLength: inp.modL, foamDensity: inp.foamRho, rhoW: line.rhoW }), 'Buoyancy module');
    const up = buoy.netBuoyancy / inp.modPitch;                     // N/m net uplift smeared along the section
    const foamMass = inp.foamRho * buoy.foamVolume / inp.modPitch;  // kg/m
    // drag on the projected area of modules + bare gaps (equivalent diameter); added mass on the smeared volume
    const Deq = (inp.modOD * inp.modL + line.D * (inp.modPitch - inp.modL)) / inp.modPitch;
    const AeB = line.Ae + buoy.foamVolume / inp.modPitch;
    buoySec = { ...line.section, w: line.w - up, D: Deq, m: line.m + foamMass, ma: (inp.Ca ?? 1) * line.rhoW * AeB };
    s0 = inp.upperLength; s1 = inp.upperLength + inp.buoyLength;
    res = run(() => solveWaveRiser({ pipe: line.section, buoyant: buoySec, upperLength: inp.upperLength, buoyancyLength: inp.buoyLength, lowerLength: lower,
      top, bottom, seabed, current, nElem, initial }), 'The lazy-wave configuration did not solve — the shape is probably not feasible or the current is too strong for it (reduce the current, lengthen the line or adjust the buoyancy section)');
  } else {
    res = run(() => solveRiser({ sections: [{ ...line.section, length: inp.length }], top, bottom, seabed, current, nElem, initial }), 'The catenary configuration did not solve — the shape is probably not feasible or the current is too strong for it (reduce the current or lengthen the line)');
  }
  need(res.converged, 'The static configuration did not converge. Check that the length, span and buoyancy describe a feasible shape (try a longer line, a smaller span or less buoyancy).');
  return { res, line, h, buoy, buoySec, s0, s1 };
}

/**
 * Apply an imported route / centreline polyline (z up): the first point is the hang-off, the last the termination.
 * Only the end positions are used — they replace the horizontal distance and the hang-off height above the seabed.
 */
export function applyRoute(inp) {
  const p = inp.route?.points;
  if (!p || p.length < 6) return inp;
  const n = Math.floor(p.length / 3), a = 3 * (n - 1);
  const span = Math.hypot(p[a] - p[0], p[a + 1] - p[1]), h = p[2] - p[a + 2];
  need(Number.isFinite(span) && Number.isFinite(h), 'The imported route contains non-numeric coordinates.');
  need(span > 5 && h > 10, `The imported route does not run from a hang-off down to a seabed termination (horizontal distance ${fmt(span)} m, height difference ${fmt(h)} m): the first point must be the hang-off and the last the termination, with z upwards.`);
  need(h < inp.depth + 40, `The imported route spans ${fmt(h)} m vertically, more than the water depth of ${inp.depth} m plus hang-off elevation: check units or the depth input.`);
  return { ...inp, span, hangoffDepth: inp.depth - h };
}

/** Rod current record for a profile function of elevation. */
export const rodCurrent = (inp, rhoW) => ((inp.Us || inp.Ub) ? { rho: rhoW, profile: currentFn(inp, inp.depth) } : undefined);

/** 3-D line plot (z up) of one or more planar shapes with nodal values. */
export function shape3d(title, shapes, { valueLabel, valueUnit, explain, depth, y = 0 } = {}) {
  const nodes = [], lines = [], nodeValues = [];
  for (const sh of shapes) {
    const [x, z, v] = thin([sh.x, sh.z, sh.v ?? sh.x.map(() => 0)], 200);
    const base = nodes.length;
    for (let i = 0; i < x.length; i++) { nodes.push([x[i], sh.y ?? y, z[i] - (depth ?? 0)]); nodeValues.push(v[i]); if (i) lines.push([base + i - 1, base + i]); }
  }
  return { type: 'model3d', title, nodes, lines, tris: [], nodeValues, valueLabel, valueUnit, explain };
}
