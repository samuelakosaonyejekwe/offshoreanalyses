// Shared helpers of the subsea support structure suite: unit handling, common input groups, screening
// sustainability factors and advice/check records. Engineering units are converted to SI inside each module
// (docs/CONTRACT.md §2.1); the helpers here are unit-neutral unless stated.

export const G = 9.80665;
export const DEG = Math.PI / 180;
// Screening emission factors (order of magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m³ concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 20;   // construction / survey vessel fuel burn, t/day (screening)
export const NU_SEA = 1.36e-6;     // kinematic viscosity of sea water at about 10 °C, m²/s

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible structure.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
/** Safety factors can be infinite when the driving action vanishes: report a finite cap instead. */
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
/** Check written as required factor of safety over achieved factor of safety. */
export const fsCheck = (name, fs, fsReq, ref, explain) => (
  { name, demand: round(fsReq, 4), capacity: round(cap(fs), 4), unit: 'FS (required / achieved)', uc: fs > 0 ? Math.min(fsReq / fs, 99) : 99, ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

/** Sustainability sentences built on the screening factors above. */
export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (screening factor — confirm with the mill's product declaration).`;
export const concreteCo2Text = (m3) => `${fmt(Math.abs(m3))} m³ of concrete ≈ ${fmt(Math.abs(m3) * CO2_CONCRETE)} t CO₂e at the screening factor of ${CO2_CONCRETE} t CO₂e per m³.`;
export const vesselCo2Text = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const CONFIRM = 'Default recalled from the standard — confirm against the governing project edition.';
export const NOT_FE = 'Large-deformation soil finite elements, coupled Eulerian–Lagrangian and particle methods, and coupled effective-stress continuum analysis are outside this module.';

/* ------------------------------------------------------------------------------------------------ inputs */

/** Mudmat plan dimensions. B is the side in the direction of the governing horizontal load. */
export function matInputs({ B = 10, L = 14, group = 'Mudmat geometry', design = true } = {}) {
  return [
    { key: 'B', label: 'Mudmat width B', unit: 'm', type: 'number', default: B, min: 0.5, max: 60, group, role: 'geometry', design, bind: 'subsea.matB',
      help: 'Plan dimension along x, the direction of the governing horizontal load.' },
    { key: 'L', label: 'Mudmat length L', unit: 'm', type: 'number', default: L, min: 0.5, max: 80, group, role: 'geometry', design, bind: 'subsea.matL',
      help: 'Plan dimension along y.' },
  ];
}

/** Soil strength inputs shared by the geotechnical modules. */
export function soilInputs({ soil = 'clay', su0 = 5, suGrad = 1.5, phi = 32, gammaSub = 6, group = 'Soil' } = {}) {
  return [
    { key: 'soil', label: 'Soil response', type: 'select', default: soil, group, role: 'soil',
      options: [{ value: 'clay', label: 'Clay — undrained (su)' }, { value: 'sand', label: 'Sand — drained (c′, φ′)' }],
      help: 'Undrained for clay under short-term loading; drained for sand and for long-term loading of clay.' },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: su0, min: 0.2, max: 500, group, role: 'soil', bind: 'site.su0',
      uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'soil', in: ['clay'] }, help: 'From CPT / T-bar / laboratory tests. Soft deep-water clay is typically 2–10 kPa at the mudline.' },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: suGrad, min: 0, max: 50, group, role: 'soil',
      uncertain: { dist: 'normal', cov: 0.2 }, showIf: { key: 'soil', in: ['clay'] }, help: 'Increase of su with depth; 1–2 kPa/m for normally consolidated clay.' },
    { key: 'phiDeg', label: 'Effective friction angle φ′', unit: 'deg', type: 'number', default: phi, min: 15, max: 45, group, role: 'soil', bind: 'site.phiDeg',
      uncertain: { dist: 'normal', cov: 0.07 }, showIf: { key: 'soil', in: ['sand'] } },
    { key: 'cEff', label: 'Effective cohesion c′', unit: 'kPa', type: 'number', default: 0, min: 0, max: 100, group, role: 'soil', showIf: { key: 'soil', in: ['sand'] } },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: gammaSub, min: 2, max: 13, group, role: 'soil', bind: 'site.gammaSub',
      uncertain: { dist: 'normal', cov: 0.05 } },
  ];
}

export function waterInputs({ depth = 150, group = 'Site' } = {}) {
  return [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: depth, min: 5, max: 3500, group, role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group, role: 'environment', bind: 'site.rhoW' },
  ];
}

/** Resultant foundation loads at mudmat centre, mudline level (published by the interface-reaction module). */
export function loadInputs({ V = 1000, H = 120, M = 600, T = 0, group = 'Design loads at mudline', torsion = false } = {}) {
  const out = [
    { key: 'V', label: 'Vertical load V', unit: 'kN', type: 'number', default: V, min: 1, max: 2e5, group, role: 'bc', bind: 'subsea.designV',
      uncertain: { dist: 'normal', cov: 0.05 }, help: 'Submerged weight of structure and equipment plus vertical interface reactions, compression positive.' },
    { key: 'H', label: 'Horizontal load H (along B)', unit: 'kN', type: 'number', default: H, min: 0, max: 1e5, group, role: 'bc', bind: 'subsea.designH',
      uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Resultant of hydrodynamic load and pipeline / jumper reactions.' },
    { key: 'M', label: 'Overturning moment M (eccentricity along B)', unit: 'kN·m', type: 'number', default: M, min: 0, max: 1e6, group, role: 'bc', bind: 'subsea.designM',
      uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Moment at mudline about the mat centre-line, including H × lever arm.' },
  ];
  if (torsion) out.push({ key: 'T', label: 'Torsion T about the vertical axis', unit: 'kN·m', type: 'number', default: T, min: 0, max: 1e6, group, role: 'bc', bind: 'subsea.designT',
    help: 'From eccentric pipeline / jumper pulls.' });
  return out;
}

/** Soil parameters in SI from the shared soil inputs. suAt(z) is the undrained strength profile. */
export function soilFrom(inp) {
  const clay = inp.soil === 'clay';
  const su0 = (inp.su0 ?? 0) * 1e3, k = (inp.suGrad ?? 0) * 1e3;
  return { clay, su0, k, suAt: (z) => su0 + k * z, phi: clay ? 0 : inp.phiDeg * DEG, c: clay ? 0 : (inp.cEff ?? 0) * 1e3, gamma: inp.gammaSub * 1e3 };
}
