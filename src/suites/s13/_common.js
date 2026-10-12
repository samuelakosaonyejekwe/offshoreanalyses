// Shared helpers of the submarine power and communication cable suite: unit handling, the parametric cable
// construction (form fields → IEC 60287 geometry and → mechanical cross-section), material take-off per metre,
// the global line properties and dynamic configuration (free-hanging catenary or lazy wave on the 2-D rod solver),
// screening sustainability factors and advice records.
// Engineering units are converted to SI here so that every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { cableGeometry } from '../../engine/electrical/ampacity.js';
import { buildCable } from '../../engine/xsection/cablemech.js';
import { solveRiser, solveWaveRiser } from '../../engine/lines/riser.js';
import { buoyancyModules } from '../../engine/xsection/ancillaries.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const MU0 = 4e-7 * Math.PI;
export const HOURS_YEAR = 8766;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 25;   // typical cable-lay / trenching-support vessel fuel burn, t/day (screening)
export const CO2_COPPER = 4.0;     // t CO2e per t refined copper (screening, typical 3–5)
export const CO2_ALUMINIUM = 12;   // t CO2e per t primary aluminium (screening, 8–17 by smelter power mix)
export const CO2_LEAD = 1.5;       // t CO2e per t lead (screening)
export const CO2_POLYMER = 2.5;    // t CO2e per t polyethylene / polypropylene (screening)
export const CO2_ROCK = 0.01;      // t CO2e per t quarried rock at the quay (screening, excludes placement vessel)
export const CO2_GRID = 0.4;       // t CO2 per MWh replaced when lost energy is made up by a gas-dominated grid (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible cable.`);
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
export const logspace = (a, b, n) => linspace(Math.log(a), Math.log(b), n).map(Math.exp);

/* ------------------------------------------------------------------------------------------------ advice */

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? uc : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc, warnAt = 0.9) => (uc > 1 ? 'fail' : uc > warnAt ? 'warn' : 'ok');

export const vesselCo2Text = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;
export const lossCo2Text = (mwhYr) => `${fmt(mwhYr)} MWh/yr of losses ≈ ${fmt(mwhYr * CO2_GRID)} t CO₂ a year if made up by gas-fired generation (screening factor ${CO2_GRID} t CO₂/MWh — use the grid factor of the project).`;

export const LIMIT_2D = 'The rod model is two-dimensional (in-plane): out-of-plane current, cross offsets, torsion and lateral seabed friction are not modelled, and the seabed is a row of linear lift-off springs.';
export const LIMIT_GENERIC = 'The parametric cable construction is generic and indicative: replace it with the manufacturer data sheet and type-test results before design use.';
export const LIMIT_FE = 'Two- and three-dimensional electromagnetic finite elements are outside the suite.';
export const CONFIRM = 'Overridable default recalled from the document — confirm against the governing edition and the cable manufacturer.';
export const STD_CABLE = 'DNV-ST-0359 (2021) / DNV-RP-0360 (2016)';
export const TB880 = 'IEC 60287-1-1 is known to overestimate the armour loss of three-core submarine cables (CIGRE TB 880, 2022).';

/* ----------------------------------------------------------------------------------- cable construction */

// Library conductor ids → metal keys of the electrical engine.
const METAL_KEY = { cu: 'Cu', al: 'Al', 'lead-sheath': 'lead', 'steel-armour': 'steel' };
const RHO_T_PE = 3.5;        // K·m/W, polyethylene coverings (IEC 60287-2-1 Table 1)
const RHO_T_BEDDING = 6.0;   // K·m/W, fibrous / polypropylene-yarn bedding and serving (IEC 60287-2-1 Table 1)

/**
 * Parametric power-cable construction inputs (conductor outwards). Defaults: generic 66 kV three-core
 * 630 mm² copper / XLPE / lead-sheath / single-wire-armoured cable of about 150 mm diameter (indicative data).
 */
export function cableInputs({ mech = false, group = 'Cable construction' } = {}) {
  const g = group;
  return [
    { key: 'sectionGeom', label: 'Imported cross-section outline (optional)', type: 'geometry', kinds: ['section'], default: null, group: g, role: 'geometry',
      help: 'DXF / SVG / JSON section. Only the outer boundary is used (it must enclose the parametric build-up); the layer build-up always comes from the fields below, which are the parametric fallback.' },
    { key: 'cores', label: 'Cores', type: 'select', default: '3', group: g, role: 'geometry', options: [{ value: '3', label: 'Three-core (separately sheathed cores, AC)' }, { value: '1', label: 'Single-core' }] },
    { key: 'condMat', label: 'Conductor metal', type: 'select', source: 'materials:conductor', default: 'cu', group: g, role: 'material', help: 'Copper or aluminium (the sheath and armour entries of the library are not conductor metals).' },
    { key: 'area', label: 'Conductor cross-section', unit: 'mm²', type: 'number', default: 630, min: 35, max: 3000, group: g, role: 'geometry', design: true, help: 'Nominal IEC 60228 class-2 section per core.' },
    { key: 'insMat', label: 'Insulation', type: 'select', source: 'materials:insulation', default: 'xlpe', group: g, role: 'material' },
    { key: 'tScreenC', label: 'Conductor screen thickness', unit: 'mm', type: 'number', default: 1.0, min: 0, max: 4, group: g, role: 'geometry' },
    { key: 'tIns', label: 'Insulation thickness', unit: 'mm', type: 'number', default: 9, min: 2, max: 35, group: g, role: 'geometry', design: true, help: 'About 5.5 mm at 33 kV, 9–11 mm at 66 kV, 15–18 mm at 132 kV, 23–27 mm at 220 kV (indicative).' },
    { key: 'tScreenI', label: 'Insulation screen and swelling tapes', unit: 'mm', type: 'number', default: 1.0, min: 0, max: 5, group: g, role: 'geometry' },
    { key: 'sheathMat', label: 'Metallic sheath / screen metal', type: 'select', source: 'materials:conductor', default: 'lead-sheath', group: g, role: 'material' },
    { key: 'tSheath', label: 'Metallic sheath thickness', unit: 'mm', type: 'number', default: 1.8, min: 0, max: 6, group: g, role: 'geometry', help: '0 = no metallic sheath (wet design without a radial water barrier).' },
    { key: 'tJacket', label: 'Core jacket (polyethylene) thickness', unit: 'mm', type: 'number', default: 2.0, min: 0, max: 8, group: g, role: 'geometry' },
    { key: 'tBed', label: 'Armour bedding thickness', unit: 'mm', type: 'number', default: 2.0, min: 0, max: 10, group: g, role: 'geometry' },
    { key: 'armourMat', label: 'Armour wire metal', type: 'select', source: 'materials:conductor', default: 'steel-armour', group: g, role: 'material' },
    { key: 'dWire', label: 'Armour wire diameter', unit: 'mm', type: 'number', default: 5, min: 0, max: 9, group: g, role: 'geometry', design: true, help: '0 = unarmoured. Galvanised steel 4–8 mm is usual.' },
    ...(mech ? [
      { key: 'nArmour', label: 'Armour layers', type: 'select', default: '2', group: g, role: 'geometry', options: [{ value: '1', label: 'Single layer' }, { value: '2', label: 'Two contra-wound layers (torque-balanced, dynamic)' }] },
      { key: 'armAngle', label: 'Armour lay angle from the axis', unit: 'deg', type: 'number', default: 12, min: 4, max: 40, group: g, role: 'geometry', help: 'Lay length = 2πR / tan(angle). Long lays (small angle) give axial stiffness; short lays give flexibility.' },
      { key: 'coreLay', label: 'Core lay length (three-core)', unit: 'm', type: 'number', default: 2.2, min: 0.5, max: 10, group: g, role: 'geometry', help: 'Pitch of the laid-up cores, typically 10–20 laid-up diameters.' },
      { key: 'armUts', label: 'Armour wire tensile strength', unit: 'MPa', type: 'number', default: 500, min: 300, max: 1900, group: g, role: 'material', uncertain: { dist: 'lognormal', cov: 0.05 }, help: 'Galvanised mild steel grade 34/65 ≈ 340–700 MPa; high-tensile wires for deep water.' },
    ] : []),
    { key: 'tServ', label: 'Outer serving / sheath thickness', unit: 'mm', type: 'number', default: 4, min: 0.5, max: 12, group: g, role: 'geometry' },
    ...(mech ? [] : [
      { key: 'ks', label: 'Skin-effect coefficient ks', unit: '-', type: 'number', default: 1, min: 0.2, max: 1, group: 'IEC 60287 constants', role: 'criteria', help: `IEC 60287-1-1 Table 2: 1 for round stranded conductors, 0.435 for Milliken (0.35–0.80 by wire treatment in the 2014 amendment). ${CONFIRM}` },
      { key: 'kp', label: 'Proximity-effect coefficient kp', unit: '-', type: 'number', default: 1, min: 0.2, max: 1, group: 'IEC 60287 constants', role: 'criteria', help: `IEC 60287-1-1 Table 2: 1 for round stranded (0.8 dried and impregnated), 0.37 for Milliken. ${CONFIRM}` },
      { key: 'rhoTCover', label: 'Thermal resistivity of polyethylene coverings', unit: 'K·m/W', type: 'number', default: RHO_T_PE, min: 1, max: 8, group: 'IEC 60287 constants', role: 'criteria', help: `IEC 60287-2-1 Table 1: polyethylene 3.5, PVC 5.0–6.0. Applies to the core jackets and the serving. ${CONFIRM}` },
      { key: 'rhoTBed', label: 'Thermal resistivity of bedding and fillers', unit: 'K·m/W', type: 'number', default: RHO_T_BEDDING, min: 1, max: 10, group: 'IEC 60287 constants', role: 'criteria', help: `IEC 60287-2-1 Table 1: fibrous / polypropylene-yarn bedding 6.0; water-filled interstices behave better. ${CONFIRM}` },
    ]),
  ];
}

function importedDiameter(geom) {
  if (!geom || !Array.isArray(geom.loops) || !geom.loops.length) return null;
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const loop of geom.loops) for (let i = 0; i + 1 < loop.length; i += 2) {
    xmin = Math.min(xmin, loop[i]); xmax = Math.max(xmax, loop[i]); ymin = Math.min(ymin, loop[i + 1]); ymax = Math.max(ymax, loop[i + 1]);
  }
  const d = Math.max(xmax - xmin, ymax - ymin);
  return Number.isFinite(d) && d > 0 ? d : null;
}

function mats(inp, ctx) {
  const cond = ctx.material(inp.condMat), ins = ctx.material(inp.insMat), sh = ctx.material(inp.sheathMat), ar = ctx.material(inp.armourMat);
  need(cond && ['cu', 'al'].includes(cond.id), `Conductor metal '${inp.condMat}' must be copper or aluminium.`);
  need(ins && ins.rhoT > 0 && ins.epsR > 0, `Insulation '${inp.insMat}' lacks thermal resistivity or permittivity in the materials library.`);
  need(sh && METAL_KEY[sh.id], `Sheath metal '${inp.sheathMat}' is not a cable metal of the library.`);
  need(ar && METAL_KEY[ar.id], `Armour metal '${inp.armourMat}' is not a cable metal of the library.`);
  return { cond, ins, sh, ar };
}

/**
 * Electrical / thermal description for the IEC 60287 engine and its geometry record.
 * Returns { desc, geom, m (library records), cores, De (m) }.
 */
export function buildElec(inp, ctx, { T1, T2, T3 } = {}) {
  const m = mats(inp, ctx), cores = Number(inp.cores);
  const mm = (v) => (v || 0) / 1e3;
  const desc = {
    cores,
    conductor: { material: METAL_KEY[m.cond.id], area: inp.area * 1e-6, ...(inp.ks != null ? { ks: inp.ks } : {}), ...(inp.kp != null ? { kp: inp.kp } : {}) },
    conductorScreenThickness: mm(inp.tScreenC),
    insulation: { thickness: mm(inp.tIns), rhoT: m.ins.rhoT, epsR: m.ins.epsR, tanDelta: m.ins.tanDelta ?? 0.001, heatCapacity: m.ins.cp && m.ins.rho ? m.ins.cp * m.ins.rho : 2.4e6 },
    insulationScreenThickness: mm(inp.tScreenI),
    ...(inp.tSheath > 0 ? { sheath: { material: METAL_KEY[m.sh.id], thickness: mm(inp.tSheath) } } : {}),
    ...(inp.tJacket > 0 ? { coreJacket: { thickness: mm(inp.tJacket), rhoT: inp.rhoTCover ?? RHO_T_PE } } : {}),
    ...(inp.tBed > 0 ? { bedding: { thickness: mm(inp.tBed), rhoT: inp.rhoTBed ?? RHO_T_BEDDING } } : {}),
    ...(inp.dWire > 0 ? { armour: { material: METAL_KEY[m.ar.id], wireDiameter: mm(inp.dWire), magnetic: m.ar.id === 'steel-armour' && cores === 3 } } : {}),
    serving: { thickness: mm(inp.tServ), rhoT: inp.rhoTCover ?? RHO_T_PE },
    ...(T1 != null ? { T1 } : {}), ...(T2 != null ? { T2 } : {}), ...(T3 != null ? { T3 } : {}),
  };
  let geom;
  try { geom = cableGeometry(desc); } catch (e) {
    if (/thickness ratio/.test(e.message)) {
      // IEC 60287-2-1 Fig. 6 (geometric factor of three-core SL cables) covers bedding thickness ≤ 0.15 × core diameter
      const { bedding, ...rest } = desc, dCore = cableGeometry(rest).coreDiameter;
      throw new Error(`The armour bedding of ${inp.tBed} mm is outside the IEC 60287-2-1 geometric factor for three-core cables, which covers bedding up to 15 % of the core diameter (${fmt(0.15 * dCore * 1e3)} mm for this core of ${fmt(dCore * 1e3)} mm): reduce the bedding thickness.`);
    }
    throw new Error(`The cable construction could not be built: ${e.message}`);
  }
  const imp = importedDiameter(inp.sectionGeom);
  if (imp) need(imp >= geom.De * 0.98, `The imported outline (${fmt(imp * 1e3, 4)} mm) is smaller than the parametric build-up (${fmt(geom.De * 1e3, 4)} mm): check the drawing scale or the layer thicknesses.`);
  return { desc, geom, m, cores, De: geom.De };
}

/**
 * Mass per metre by material [kg/m] from the IEC geometry record: conductor metal, lead (or other sheath metal),
 * steel armour, polymers (insulation, screens, jackets, bedding with fillers, serving).
 */
export function massesPerMetre(built, { fillerRho = 600 } = {}) {
  const { geom, m, cores } = built;
  let conductor = 0, sheath = 0, armour = 0, polymer = 0;
  for (const L of geom.layers) {
    const n = ['conductor', 'insulation', 'sheath', 'coreJacket'].includes(L.name) ? cores : 1;
    if (L.name === 'conductor') conductor += n * L.area * m.cond.rho;
    else if (L.name === 'sheath') sheath += n * L.area * m.sh.rho;
    else if (L.name === 'armour') armour += L.area * m.ar.rho;
    else if (L.name === 'insulation') polymer += n * L.area * (m.ins.rho ?? 920);
    else if (L.name === 'bedding') polymer += L.area * fillerRho;       // bedding region includes the interstitial fillers
    else polymer += n * L.area * 950;
  }
  const al = m.cond.id === 'al';
  const leadSheath = m.sh.id === 'lead-sheath';
  return {
    copper: (al ? 0 : conductor) + (m.sh.id === 'cu' ? sheath : 0), aluminium: (al ? conductor : 0) + (m.sh.id === 'al' ? sheath : 0),
    lead: leadSheath ? sheath : 0, steel: (m.ar.id === 'steel-armour' ? armour : 0) + (m.sh.id === 'steel-armour' ? sheath : 0),
    polymer, total: conductor + sheath + armour + polymer,
  };
}

/** Quantities record (tonnes, km) for a route length [m]. Lead has no canonical key and is reported as `lead_t`. */
export function cableQuantities(mass, lengthM, extra = {}) {
  const t = (kgm) => round(kgm * lengthM / 1e3, 5);
  return { cable_km: round(lengthM / 1e3, 6), copper_t: t(mass.copper), aluminium_t: t(mass.aluminium), lead_t: t(mass.lead), steel_t: t(mass.steel), polymer_t: t(mass.polymer), ...extra };
}
export const embodiedCo2 = (mass, lengthM) => (mass.copper * CO2_COPPER + mass.aluminium * CO2_ALUMINIUM + mass.lead * CO2_LEAD + mass.steel * CO2_STEEL + mass.polymer * CO2_POLYMER) * lengthM / 1e3;

/**
 * Mechanical cross-section (engine buildCable) from the same form fields.
 * Polymer moduli are room-temperature handbook values; the result brackets the bending stiffness (EImin / EImax).
 */
export function buildMech(inp, ctx, { rhoW = 1025 } = {}) {
  const m = mats(inp, ctx), cores = Number(inp.cores);
  const mm = (v) => (v || 0) / 1e3;
  const sol = (name, t, E, rho) => (t > 0 ? [{ type: 'solid', name, thickness: mm(t), E, rho, nu: 0.42 }] : []);
  const insE = m.ins.E ?? 0.2e9, insRho = m.ins.rho ?? 920;
  const core = [
    { type: 'conductor', name: 'Conductor', area: inp.area * 1e-6, fillFactor: 0.9, E: m.cond.E ?? 110e9, rho: m.cond.rho, strandDiameter: 0.003 },
    ...sol('Conductor screen', inp.tScreenC, insE, insRho), ...sol('Insulation', inp.tIns, insE, insRho), ...sol('Insulation screen', inp.tScreenI, insE, insRho),
    ...(inp.tSheath > 0 ? [{ type: 'sheath', name: 'Metallic sheath', thickness: mm(inp.tSheath), E: m.sh.E ?? 16e9, rho: m.sh.rho, nu: 0.42 }] : []),
    ...sol('Core jacket', inp.tJacket, 0.8e9, 950),
  ];
  const nArm = inp.dWire > 0 ? Number(inp.nArmour ?? 1) : 0;
  const arm = (lay, i) => ({ type: 'armour', name: `Armour layer ${i}`, wireDiameter: mm(inp.dWire), alphaDeg: inp.armAngle ?? 12, lay, E: m.ar.E ?? 200e9, rho: m.ar.rho, uts: (inp.armUts ?? 500) * 1e6 });
  const outer = [
    ...sol('Bedding', inp.tBed, 0.05e9, 600),
    ...(nArm >= 1 ? [arm(1, 1)] : []), ...(nArm >= 2 ? [arm(-1, 2)] : []),
    ...sol('Serving', inp.tServ, 0.8e9, 950),
  ];
  const def = { cores, core, coreLayLength: inp.coreLay ?? 2.2, coreLay: -1, fillerDensity: 600, outer };
  const cab = run(() => buildCable(def, { rhoSeawater: rhoW }), 'The mechanical cross-section could not be built');
  return { cab, def, m, cores, nArm, armourLayers: cab.layers.filter((L) => L.type === 'armour') };
}

/* ------------------------------------------------------------------------------------ global line properties */

/** Equivalent-beam properties of the cable for the global (rod) analyses; defaults match the generic dynamic cable of the cross-section module. */
export function lineInputs({ group = 'Cable (global properties)', hydro = true } = {}) {
  return [
    { key: 'OD', label: 'Overall diameter', unit: 'mm', type: 'number', default: 160, min: 15, max: 400, group, role: 'geometry', bind: 'cable.OD',
      help: 'Outside diameter over the serving (buoyancy and drag).' },
    { key: 'mass', label: 'Mass in air', unit: 'kg/m', type: 'number', default: 63, min: 0.3, max: 300, group, role: 'geometry', bind: 'cable.massPerM',
      uncertain: { dist: 'normal', cov: 0.03 }, help: 'From the tension–bending–torsion module or the manufacturer data sheet.' },
    { key: 'EA', label: 'Axial stiffness EA', unit: 'MN', type: 'number', default: 800, min: 1, max: 5000, group, role: 'material', bind: 'cable.EA',
      uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Measured in the tensile test where available; the calculated value is an upper estimate.' },
    { key: 'EI', label: 'Bending stiffness EI', unit: 'N·m²', type: 'number', default: 12000, min: 50, max: 5e6, group, role: 'material', bind: 'cable.EI', design: true,
      uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Slip (large-curvature) bending stiffness; temperature and curvature dependent — use the measured value for dynamic analysis.' },
    ...(hydro ? [
      { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Hydrodynamic coefficients', role: 'load',
        uncertain: { dist: 'normal', cov: 0.15 }, help: 'Bare served cable 1.0–1.2 (DNV-RP-C205 §6.7); raise it for marine growth or VIV. Calibrate from oscillatory-flow tests where available.' },
      { key: 'Cdt', label: 'Tangential drag coefficient', unit: '-', type: 'number', default: 0.01, min: 0, max: 0.2, group: 'Hydrodynamic coefficients', role: 'load' },
      { key: 'Ca', label: 'Added-mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 2.5, group: 'Hydrodynamic coefficients', role: 'load' },
    ] : []),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
  ];
}

/** Rod-section record (SI) from the global line inputs. */
export function lineFrom(inp) {
  const D = inp.OD / 1e3, rhoW = inp.rhoW ?? 1025;
  const Ae = Math.PI * D * D / 4;
  const w = (inp.mass - rhoW * Ae) * G;
  const ma = (inp.Ca ?? 1) * rhoW * Ae;
  const EA = inp.EA * 1e6, EI = inp.EI;
  return { D, Ae, w, m: inp.mass, ma, EA, EI, rhoW,
    section: { EA, EI, w, D, Cdn: inp.Cdn ?? 1.2, Cdt: inp.Cdt ?? 0.01, m: inp.mass, ma } };
}

/** Configuration inputs: free-hanging catenary or lazy wave between a floater hang-off and a seabed anchor point. */
export function configInputs({ config = 'lazy-wave' } = {}) {
  const wave = { key: 'config', in: ['lazy-wave'] };
  return [
    { key: 'config', label: 'Configuration', type: 'select', default: config, group: 'Configuration', role: 'ic',
      options: [{ value: 'lazy-wave', label: 'Lazy wave (distributed buoyancy)' }, { value: 'catenary', label: 'Free-hanging catenary' }],
      help: 'Installed geometry that the analysis starts from (initial equilibrium under apparent weight and seabed support).' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 150, min: 20, max: 2000, group: 'Configuration', role: 'environment', bind: 'site.depth' },
    { key: 'hangoffDepth', label: 'Hang-off depth below still water', unit: 'm', type: 'number', default: 12, min: -40, max: 200, group: 'Configuration', role: 'bc',
      help: 'Elevation of the I-tube exit / hang-off on the floater; negative = above the water line. The hang-off follows the platform.' },
    { key: 'length', label: 'Dynamic length hang-off to seabed anchor point', unit: 'm', type: 'number', default: 330, min: 30, max: 8000, group: 'Configuration', role: 'geometry', design: true },
    { key: 'span', label: 'Horizontal distance hang-off to anchor point', unit: 'm', type: 'number', default: 215, min: 5, max: 7500, group: 'Configuration', role: 'bc', design: true,
      help: 'Position of the seabed tether / transition joint relative to the hang-off; that end is pinned on the seabed.' },
    { key: 'upperLength', label: 'Upper catenary length (hang-off to first buoyancy module)', unit: 'm', type: 'number', default: 160, min: 10, max: 5000, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'buoyLength', label: 'Buoyancy-section length', unit: 'm', type: 'number', default: 80, min: 5, max: 2000, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'modOD', label: 'Buoyancy-module outside diameter', unit: 'm', type: 'number', default: 0.6, min: 0.1, max: 3, group: 'Buoyancy section', role: 'geometry', showIf: wave },
    { key: 'modL', label: 'Buoyancy-module length', unit: 'm', type: 'number', default: 0.9, min: 0.2, max: 4, group: 'Buoyancy section', role: 'geometry', showIf: wave },
    { key: 'modPitch', label: 'Buoyancy-module pitch', unit: 'm', type: 'number', default: 2, min: 0.3, max: 30, group: 'Buoyancy section', role: 'geometry', showIf: wave, design: true },
    { key: 'foamRho', label: 'Syntactic-foam density', unit: 'kg/m³', type: 'number', default: 420, min: 250, max: 800, group: 'Buoyancy section', role: 'material', showIf: wave,
      uncertain: { dist: 'normal', cov: 0.04 }, help: 'Rises with depth rating (indicative). Water absorption over life reduces uplift.' },
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 200, min: 5, max: 1e5, group: 'Seabed', role: 'soil',
      uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Force per metre of cable per metre of penetration (kN/m/m). Soft clay 20–200, sand 1000+.' },
    { key: 'muAx', label: 'Axial seabed friction coefficient', unit: '-', type: 'number', default: 0.4, min: 0, max: 1.2, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 } },
  ];
}

/** Current profile inputs: surface and near-bed speed with a power law in between (in-plane). */
export function currentInputs({ Us = 0.8, Ub = 0.2, bind = 'site.Uc100', group = 'Current' } = {}) {
  return [
    { key: 'Us', label: 'Surface current speed', unit: 'm/s', type: 'number', default: Us, min: -3, max: 3, group, role: 'bc', bind,
      uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'In-plane current at the surface, positive from the hang-off towards the anchor point (far direction); negative = near.' },
    { key: 'Ub', label: 'Near-bed current speed', unit: 'm/s', type: 'number', default: Ub, min: -2, max: 2, group, role: 'bc', bind: 'site.UcBed100' },
    { key: 'curExp', label: 'Profile exponent', unit: '-', type: 'number', default: 1.5, min: 0.1, max: 8, group, role: 'bc',
      help: 'U(z) = Ub + (Us − Ub)·(height/depth)^n. n = 1 linear shear; larger n concentrates the current near the surface.' },
  ];
}
export function currentFn(inp, depth) {
  const { Us, Ub, curExp } = inp;
  return (z) => Ub + (Us - Ub) * Math.min(1, Math.max(0, z / depth)) ** curExp;
}
export const rodCurrent = (inp, rhoW) => ((inp.Us || inp.Ub) ? { rho: rhoW, profile: currentFn(inp, inp.depth) } : undefined);

/**
 * Build and solve the static configuration with the 2-D rod solver. Seabed is z = 0; the hang-off is at
 * (offset, depth − hangoffDepth) and the anchor point at (span, 0), pinned.
 */
export function solveConfig(inp, { nElem = 160, current, offset = 0, heave = 0, initial } = {}) {
  const line = lineFrom(inp);
  const h = inp.depth - inp.hangoffDepth;
  need(h > 10, 'Hang-off depth must leave at least 10 m between the hang-off and the seabed.');
  need(line.w > 0, `The cable is buoyant in sea water (apparent weight ${fmt(line.w)} N/m): it cannot hang as a catenary. Increase mass or reduce diameter.`);
  const span = inp.span - offset, top = { x: 0, z: h + heave }, bottom = { x: span, z: 0 };
  need(span > 1, 'The vessel offset moves the hang-off beyond the anchor point.');
  need(inp.length > Math.hypot(span, h) * 1.002, `Under-constrained geometry: the length (${fmt(inp.length, 5)} m) is shorter than the straight distance hang-off to anchor point (${fmt(Math.hypot(span, h), 5)} m).`);
  need(inp.length < span + h + (inp.config === 'lazy-wave' ? 3 * inp.buoyLength : 0) + 0.6 * span, `Over-length: ${fmt(inp.length, 5)} m cannot be laid out between a hang-off ${fmt(h, 4)} m above the seabed and an anchor point ${fmt(span, 4)} m away without piling on the seabed. Shorten the line or move the anchor point.`);
  const seabed = { z: 0, kn: inp.kn * 1e3, mu: inp.muAx };
  let res, buoy = null, buoySec = null, s0 = null, s1 = null;
  if (inp.config === 'lazy-wave') {
    const lower = inp.length - inp.upperLength - inp.buoyLength;
    need(lower > 5, `Over-constrained lengths: upper (${inp.upperLength} m) + buoyancy section (${inp.buoyLength} m) leave no lower catenary within the total length of ${inp.length} m.`);
    need(inp.modOD > line.D * 1.2, 'The buoyancy-module diameter must exceed the cable diameter by at least 20 %.');
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
 * Apply an imported route / centreline polyline (z up): the first point is the hang-off, the last the anchor point.
 * Only the end positions are used — they replace the horizontal distance and the hang-off height above the seabed.
 */
export function applyRoute(inp) {
  const p = inp.route?.points;
  if (!p || p.length < 6) return inp;
  const n = Math.floor(p.length / 3), a = 3 * (n - 1);
  const span = Math.hypot(p[a] - p[0], p[a + 1] - p[1]), h = p[2] - p[a + 2];
  need(Number.isFinite(span) && Number.isFinite(h), 'The imported route contains non-numeric coordinates.');
  need(span > 5 && h > 10, `The imported route does not run from a hang-off down to a seabed anchor point (horizontal distance ${fmt(span)} m, height difference ${fmt(h)} m): the first point must be the hang-off and the last the anchor point, with z upwards.`);
  need(h < inp.depth + 40, `The imported route spans ${fmt(h)} m vertically, more than the water depth of ${inp.depth} m plus hang-off elevation: check units or the depth input.`);
  return { ...inp, span, hangoffDepth: inp.depth - h };
}
export const routeInput = (group = 'Configuration') => ({ key: 'route', label: 'Imported centreline / route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group, role: 'geometry',
  help: 'CSV / XYZ / DXF / GeoJSON polyline, z up, first point at the hang-off and last at the seabed anchor point. Only the end positions are used; the shape between them is solved. Leave empty to use the form values.' });

/** Total 3-D length [m] of an imported polyline, or null. */
export function polylineLength(geom) {
  const p = geom?.points;
  if (!p || p.length < 6) return null;
  let L = 0;
  for (let i = 3; i + 2 < p.length; i += 3) L += Math.hypot(p[i] - p[i - 3], p[i + 1] - p[i - 2], p[i + 2] - p[i - 1]);
  return Number.isFinite(L) && L > 0 ? L : null;
}

/** 3-D line plot (z up) of one or more planar shapes with nodal values. */
export function shape3d(title, shapes, { valueLabel, valueUnit, explain, depth, y = 0 } = {}) {
  const nodes = [], lines = [], nodeValues = [];
  for (const sh of shapes) {
    const [x, z, v] = thin([sh.x, sh.z, sh.v ?? Array.from(sh.x, () => 0)], 200);
    const base = nodes.length;
    for (let i = 0; i < x.length; i++) { nodes.push([x[i], sh.y ?? y, z[i] - (depth ?? 0)]); nodeValues.push(v[i]); if (i) lines.push([base + i - 1, base + i]); }
  }
  return { type: 'model3d', title, nodes, lines, tris: [], nodeValues, valueLabel, valueUnit, explain };
}
