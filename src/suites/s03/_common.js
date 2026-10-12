// Shared helpers of the offshore wind turbine foundation suite: unit handling, the common turbine / support
// structure / soil input groups, the reference-turbine builder, soil-profile builder, advice and check records and
// the screening sustainability factors. Engineering units are converted to SI here so that every module calls the
// engines in SI (docs/CONTRACT.md §2.1). Elevations: z up, z = 0 at still-water level, mudline at z = −depth.
import { nrel5mw, scaleRotor } from '../../engine/aero/rotors.js';
import { cpSurface, bemSolve } from '../../engine/aero/bem.js';
import { operatingSchedule, controllerFromRotor } from '../../engine/aero/control.js';
import { createProfile } from '../../engine/geotech/soilprofile.js';
import { cptToLayers } from '../../engine/geotech/cpt.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order of magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m³ concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 35;   // jack-up / heavy-lift installation vessel fuel burn, t/day (screening)
export const NU_SEA = 1.36e-6;     // kinematic viscosity of sea water at about 10 °C, m²/s

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible foundation.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 600) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round((i * (len - 1)) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
export const scaled = (a, f) => Array.from(a, (v) => v * f);

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : 99;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc) => (uc > 1 ? 'fail' : uc > 0.9 ? 'warn' : 'ok');

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (screening factor — confirm with the mill's product declaration).`;
export const concreteCo2Text = (m3) => `${fmt(Math.abs(m3))} m³ of concrete ≈ ${fmt(Math.abs(m3) * CO2_CONCRETE)} t CO₂e at the screening factor of ${CO2_CONCRETE} t CO₂e per m³.`;
export const vesselCo2Text = (days) => `${fmt(days)} installation-vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

// Standing caveats carried into `limits`, input help and advice.
export const LIMIT_AERO = 'Airfoil polars are generic parametric polars keyed by relative thickness (not the measured DU/NACA tables) and rotors other than the NREL 5-MW are geometric scalings of it, not published reference machines: rotor loads are for preliminary foundation loading and must be replaced by the turbine supplier\'s load documents.';
export const LIMIT_COUPLED = 'The coupled simulator has rigid blades, a modal (linear, fore-aft) support structure, rotor-effective wind without blade-resolved turbulence or 1P/3P excitation, and linear (Airy) waves.';
export const LIMIT_PY = 'Soil reactions are API-type p–y, t–z and Q–z curves; large-diameter monopile soil reactions (PISA-type distributed moment, base shear and base moment) are not generated and enter only as user-supplied curves or multipliers. The sand initial modulus k(φ′) is a fitted curve to the API figure.';
export const LIMIT_CODE = 'Load, material and fatigue factors are editable criteria recalled from IEC 61400-3-1:2019 and DNV-ST-0126 (2021): confirm them against the edition governing the project.';
export const LIMIT_HIFI = 'Blade-resolved CFD, potential-flow panel methods and three-dimensional soil finite elements are outside the module.';
export const CONFIRM = 'Default recalled from the standard — confirm against the governing project edition.';
export const HELP_SUPPLIER = 'Preliminary value from the generic rotor model — replace with the turbine supplier\'s figure when available.';

/* ---------------------------------------------------------------------------------------------- inputs */

/** Turbine definition: rating selects a geometric scaling of the NREL 5-MW reference rotor. */
export function turbineInputs({ group = 'Turbine' } = {}) {
  return [
    { key: 'ratedMW', label: 'Rated power', unit: 'MW', type: 'number', default: 5, min: 2, max: 22, group, role: 'geometry',
      help: '5 MW uses the NREL 5-MW reference planform; other ratings are geometric scalings of it at constant tip speed (generic rotor, generic polars).' },
    { key: 'hubHeight', label: 'Hub height above still water', unit: 'm', type: 'number', default: 90, min: 50, max: 200, group, role: 'geometry',
      help: 'Must leave blade-tip clearance above the interface platform.' },
    { key: 'rnaT', label: 'Rotor–nacelle mass (0 = scaled reference)', unit: 't', type: 'number', default: 0, min: 0, max: 2500, group, role: 'geometry', bind: 'owt.rna',
      uncertain: { dist: 'normal', cov: 0.03 }, help: 'NREL 5-MW: 350 t. With 0 the reference mass is scaled with rating (exponent 2.5 on the length scale). ' + HELP_SUPPLIER },
  ];
}

/** Monopile, transition and tower geometry. */
export function supportInputs({ group = 'Support structure', embed = true, D = 6, t = 70, Lp = 32 } = {}) {
  const out = [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 25, min: 5, max: 70, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'D', label: 'Monopile outer diameter', unit: 'm', type: 'number', default: D, min: 2, max: 14, group, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.005 } },
    { key: 't', label: 'Monopile wall thickness', unit: 'mm', type: 'number', default: t, min: 20, max: 160, group, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.02 },
      help: 'API RP 2A minimum for driving is 6.35 mm + D/100.' },
  ];
  if (embed) out.push({ key: 'Lp', label: 'Embedded length', unit: 'm', type: 'number', default: Lp, min: 8, max: 90, group, role: 'ic', design: true,
    help: 'Installed penetration below the original mudline (initial condition of the in-place analysis).' });
  out.push(
    { key: 'zInt', label: 'Interface (tower base) elevation', unit: 'm', type: 'number', default: 12, min: 4, max: 30, group, role: 'geometry', help: 'Above still water; the monopile / transition piece runs from the mudline to this level.' },
    { key: 'towerDb', label: 'Tower base diameter', unit: 'm', type: 'number', default: 6, min: 2, max: 12, group, role: 'geometry', design: true },
    { key: 'towerTb', label: 'Tower base wall thickness', unit: 'mm', type: 'number', default: 35, min: 10, max: 120, group, role: 'geometry', design: true },
    { key: 'towerDt', label: 'Tower top diameter', unit: 'm', type: 'number', default: 3.87, min: 1.5, max: 10, group, role: 'geometry' },
    { key: 'towerTt', label: 'Tower top wall thickness', unit: 'mm', type: 'number', default: 25, min: 8, max: 100, group, role: 'geometry' },
    { key: 'towerMassFactor', label: 'Tower mass factor (flanges, internals)', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.3, group, role: 'geometry', help: 'Multiplier on the tower shell mass for flanges, platforms, ladders, cables and paint.' },
  );
  return out;
}

/** Foundation boundary condition of the stick model: three-spring matrix at the mudline or a clamp. */
export function foundationInputs({ group = 'Foundation boundary' } = {}) {
  return [
    { key: 'fnd', label: 'Support at the mudline', type: 'select', default: 'springs', group, role: 'bc',
      options: [{ value: 'springs', label: 'Coupled lateral / rocking springs' }, { value: 'fixed', label: 'Clamped at the mudline (rigid soil)' }],
      help: 'The clamp is the rigid-soil limit; real monopiles are 5–15 % softer.' },
    { key: 'Kh', label: 'Lateral stiffness K_LL', unit: 'GN/m', type: 'number', default: 2.0, min: 0.01, max: 500, group, role: 'bc', bind: 'foundation.Kh', showIf: { key: 'fnd', in: ['springs'] },
      uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Mudline stiffness from the pile–soil interaction module (secant at working load).' },
    { key: 'Kr', label: 'Rocking stiffness K_RR', unit: 'GN·m/rad', type: 'number', default: 208, min: 0.1, max: 1e5, group, role: 'bc', bind: 'foundation.Kr', showIf: { key: 'fnd', in: ['springs'] },
      uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Khr', label: 'Coupling stiffness K_LR', unit: 'MN', type: 'number', default: -16100, min: -5e6, max: 0, group, role: 'bc', bind: 'foundation.Khr', showIf: { key: 'fnd', in: ['springs'] },
      help: 'Negative in the sign convention [H; M] = K·[y; θ]. The matrix must stay positive definite: K_LR² < K_LL·K_RR.' },
  ];
}

const LAYER_COLS = [
  { key: 'type', label: 'Type (sand / clay)', type: 'text' }, { key: 'thickness', label: 'Thickness', unit: 'm', type: 'number' },
  { key: 'gammaSub', label: 'γ′', unit: 'kN/m³', type: 'number' }, { key: 'phiDeg', label: 'φ′ (sand)', unit: 'deg', type: 'number' },
  { key: 'Dr', label: 'Relative density (sand)', unit: '-', type: 'number' }, { key: 'suTop', label: 'su at layer top (clay)', unit: 'kPa', type: 'number' },
  { key: 'suGrad', label: 'su gradient (clay)', unit: 'kPa/m', type: 'number' }, { key: 'eps50', label: 'ε50 (clay)', unit: '-', type: 'number' },
];

/** Soil description: uniform sand, uniform clay or a layer table from the borehole / CPT interpretation. */
export function soilInputs({ group = 'Soil', soil = 'sand' } = {}) {
  return [
    { key: 'soil', label: 'Soil profile', type: 'select', default: soil, group, role: 'soil',
      options: [{ value: 'sand', label: 'Uniform sand (φ′, relative density)' }, { value: 'clay', label: 'Uniform clay (su increasing with depth)' }, { value: 'layers', label: 'Layer table (borehole interpretation)' }, { value: 'cpt', label: 'CPT record (cone resistance and sleeve friction)' }],
      help: 'Initial geotechnical state: in-situ effective stress follows from the submerged unit weights.' },
    { key: 'phiDeg', label: 'Effective friction angle φ′', unit: 'deg', type: 'number', default: 36, min: 29, max: 44, group, role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] },
      uncertain: { dist: 'normal', cov: 0.06 }, help: 'Sets the API sand p–y curves; the initial modulus k(φ′) is a fitted curve to the API figure.' },
    { key: 'Dr', label: 'Relative density', unit: '-', type: 'number', default: 0.7, min: 0.2, max: 1, group, role: 'soil', showIf: { key: 'soil', in: ['sand'] },
      help: 'Selects the API axial friction and end-bearing class (medium dense ≈ 0.5, dense ≈ 0.75, very dense ≥ 0.85). Below 0.35 the API table does not apply; the lowest tabulated class is used as a screening value and a CPT-based method is needed.' },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 40, min: 1, max: 600, group, role: 'soil', bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] },
      uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 3, min: 0, max: 30, group, role: 'soil', bind: 'site.suGrad', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'eps50', label: 'Strain at 50 % strength ε50', unit: '-', type: 'number', default: 0.008, min: 0.002, max: 0.03, group, role: 'soil', showIf: { key: 'soil', in: ['clay'] },
      help: '0.02 soft, 0.01 firm, 0.005 stiff clay (Matlock).' },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 10, min: 4, max: 13, group, role: 'soil', bind: 'site.gammaSub', showIf: { key: 'soil', in: ['sand', 'clay'] },
      uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'layers', label: 'Soil layers (top down)', type: 'table', group, role: 'soil', showIf: { key: 'soil', in: ['layers'] }, columns: LAYER_COLS,
      default: [
        { type: 'sand', thickness: 8, gammaSub: 9.5, phiDeg: 33, Dr: 0.55, suTop: 0, suGrad: 0, eps50: 0 },
        { type: 'clay', thickness: 12, gammaSub: 9, phiDeg: 0, Dr: 0, suTop: 80, suGrad: 3, eps50: 0.007 },
        { type: 'sand', thickness: 80, gammaSub: 10, phiDeg: 38, Dr: 0.8, suTop: 0, suGrad: 0, eps50: 0 },
      ],
      help: 'From the borehole log and CPT interpretation; the last layer should extend well below the pile toe.' },
    { key: 'cpt', label: 'CPT record', type: 'table', group, role: 'soil', showIf: { key: 'soil', in: ['cpt'] },
      columns: [{ key: 'z', label: 'Depth below mudline', unit: 'm', type: 'number' }, { key: 'qc', label: 'Cone resistance qc', unit: 'MPa', type: 'number' }, { key: 'fs', label: 'Sleeve friction fs', unit: 'kPa', type: 'number' }],
      default: [[1, 4, 25], [3, 7, 40], [5, 10, 60], [7, 12, 75], [9, 1.4, 55], [12, 1.6, 62], [15, 1.8, 70], [18, 2.0, 78], [21, 22, 130], [25, 26, 150], [30, 31, 180], [36, 36, 210], [45, 42, 250]].map(([z, qc, fs]) => ({ z, qc, fs })),
      help: 'Depth-averaged cone readings (about one row per metre or per layer). Layers, strength, friction angle and relative density are derived by soil-behaviour-type index; the deepest layer is extended below the pile toe.' },
    { key: 'Nkt', label: 'Cone factor Nkt', unit: '-', type: 'number', default: 15, min: 8, max: 25, group, role: 'soil', showIf: { key: 'soil', in: ['cpt'] }, help: 'su = (qt − σv)/Nkt; 10–20 depending on clay type, calibrated against laboratory strength.' },
  ];
}

export function steelInputs({ group = 'Material', id = 's355' } = {}) {
  return [
    { key: 'steel', label: 'Steel grade', type: 'select', source: 'materials:steel', default: id, group, role: 'material', help: 'Yield strength is the nominal value; reduce for plate thickness above 16 mm per EN 10025 / EN 10225.' },
  ];
}
export function steelOf(inp, ctx) {
  const m = ctx && typeof ctx.material === 'function' ? ctx.material(inp.steel) : null;
  return { E: m?.E ?? 210e9, fy: m?.fy ?? 355e6, fu: m?.fu ?? 490e6, rho: m?.rho ?? 7850, nu: m?.nu ?? 0.3, name: m?.name ?? 'S355' };
}

/* ------------------------------------------------------------------------------------------ builders */

/** Thin-walled tube section properties, D and t in metres. */
export function tube(D, t) {
  need(t > 0 && 2 * t < D, 'Wall thickness must be positive and less than half the diameter.');
  const Di = D - 2 * t;
  const A = (Math.PI * (D * D - Di * Di)) / 4, I = (Math.PI * (D ** 4 - Di ** 4)) / 64;
  return { D, t, A, I, W: I / (D / 2), Zp: (D ** 3 - Di ** 3) / 6 };
}

// API RP 2GEO Table 1 gives no axial parameters for loose sand (CPT-based methods are required there): the lowest
// tabulated class is used below a relative density of 0.35 as a screening assumption, stated in the input help.
const API_SAND = (Dr) => (Dr >= 0.85 ? 'very dense sand' : Dr >= 0.65 ? 'dense sand' : Dr >= 0.35 ? 'medium-dense sand' : 'medium-dense sand-silt');

/** Soil profile (engine object) from the shared soil inputs. Depth z downwards from the original mudline. */
export function soilProfile(inp) {
  let layers;
  if (inp.soil === 'layers') {
    const rows = (inp.layers || []).filter((r) => Number(r.thickness) > 0);
    need(rows.length >= 1, 'The soil layer table is empty: enter at least one layer with a positive thickness.');
    let top = 0;
    layers = rows.map((r, i) => {
      const type = String(r.type || '').trim().toLowerCase();
      need(type === 'sand' || type === 'clay', `Soil layer ${i + 1}: type must be "sand" or "clay".`);
      need(Number(r.gammaSub) > 0, `Soil layer ${i + 1}: submerged unit weight must be positive.`);
      const base = { type, top, thickness: Number(r.thickness), gammaSub: Number(r.gammaSub) * 1e3 };
      top += Number(r.thickness);
      if (type === 'sand') {
        need(r.phiDeg >= 29 && r.phiDeg <= 45, `Soil layer ${i + 1}: friction angle must be between 29° and 45° for the API sand curves.`);
        const Dr = Number(r.Dr) > 0 ? Number(r.Dr) : 0.6;
        return { ...base, phiDeg: Number(r.phiDeg), Dr, apiSand: API_SAND(Dr) };
      }
      need(Number(r.suTop) > 0, `Soil layer ${i + 1}: undrained strength at the layer top must be positive.`);
      return { ...base, suTop: Number(r.suTop) * 1e3, suGrad: (Number(r.suGrad) || 0) * 1e3, eps50: Number(r.eps50) > 0 ? Number(r.eps50) : 0.01, OCR: 1 };
    });
  } else if (inp.soil === 'cpt') {
    const rows = (inp.cpt || []).map((r) => ({ z: Number(r.z), qc: Number(r.qc) * 1e6, fs: Number(r.fs) * 1e3 })).filter((r) => r.z > 0 && r.qc > 0 && r.fs >= 0).sort((a, b) => a.z - b.z);
    need(rows.length >= 3, 'The CPT table needs at least three depths with cone resistance and sleeve friction.');
    for (let i = 1; i < rows.length; i++) need(rows[i].z > rows[i - 1].z, 'The CPT table has repeated depths.');
    const res = run(() => cptToLayers({ z: rows.map((r) => r.z), qc: rows.map((r) => r.qc), fs: rows.map((r) => r.fs), Nkt: inp.Nkt ?? 15, eps50: 0.01 }), 'CPT interpretation');
    layers = res.layers.map((l, i, all) => {
      const thickness = i === all.length - 1 ? Math.max(l.thickness, 300 - l.top) : l.thickness;
      const cohesive = l.suTop > 0 && !(l.phiDeg > 0);
      if (cohesive) return { type: 'clay', top: l.top, thickness, gammaSub: l.gammaSub, suTop: l.suTop, suGrad: Math.max(l.suGrad || 0, 0), eps50: l.eps50 ?? 0.01, OCR: l.OCR ?? 1 };
      // the API sand curves are defined for 29°–45°: correlations outside that range are clamped
      const Dr = Math.min(Math.max(l.Dr ?? 0.6, 0.05), 1);
      return { type: 'sand', top: l.top, thickness, gammaSub: l.gammaSub, phiDeg: Math.min(Math.max(l.phiDeg ?? 33, 29), 44), Dr, apiSand: API_SAND(Dr) };
    });
  } else if (inp.soil === 'clay') {
    layers = [{ type: 'clay', top: 0, thickness: 300, gammaSub: inp.gammaSub * 1e3, suTop: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, eps50: inp.eps50, OCR: 1 }];
  } else {
    layers = [{ type: 'sand', top: 0, thickness: 300, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phiDeg, Dr: inp.Dr, apiSand: API_SAND(inp.Dr) }];
  }
  return run(() => createProfile(layers), 'Soil profile');
}
export const soilLabel = (inp) => (inp.soil === 'cpt' ? `profile interpreted from ${(inp.cpt || []).length} CPT readings` : inp.soil === 'layers' ? `${(inp.layers || []).length}-layer profile` : inp.soil === 'clay' ? `clay, su = ${fmt(inp.su0)} + ${fmt(inp.suGrad)}·z kPa` : `sand, φ′ = ${fmt(inp.phiDeg)}°`);

// Cp/Ct/Cq surface of the reference planform. With Reynolds-independent polars the coefficients depend on tip-speed
// ratio and pitch only, so one surface serves every geometric scaling; it is cached for the session.
let SURF = null;
const TSR = [0.5, 1, 2, 3, 4, 5, 6, 7, 7.5, 8, 8.5, 9, 10, 11, 12, 14, 16];
const PITCH = [-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 23, 26, 30, 40, 60, 90];
const TCACHE = new Map();

/**
 * Turbine record from the shared inputs: rotor, coefficient surface, steady operating schedule and baseline
 * controller. Results are cached by (rating, hub height, RNA mass, air density) because the generic studies call
 * compute hundreds of times.
 */
export function turbine(inp, { rhoAir = 1.225 } = {}) {
  const key = [inp.ratedMW, inp.hubHeight, inp.rnaT || 0, rhoAir].join('|');
  if (TCACHE.has(key)) return TCACHE.get(key);
  const ref = nrel5mw();
  if (!SURF) SURF = run(() => cpSurface({ rotor: ref, tsr: TSR, pitchDeg: PITCH }), 'Rotor coefficient surface');
  const P = inp.ratedMW * 1e6;
  const scaledRotor = Math.abs(inp.ratedMW - 5) > 1e-9;
  const rotor = scaledRotor ? run(() => scaleRotor(ref, { ratedPower: P, hubHeight: inp.hubHeight, massExponent: 2.5 }), 'Rotor scaling') : { ...ref, hubHeight: inp.hubHeight };
  const s = scaledRotor ? rotor.scaleFactor : 1;
  need(inp.hubHeight > rotor.tipRadius + 15, `Hub height ${inp.hubHeight} m leaves less than 15 m between the blade tip (rotor radius ${fmt(rotor.tipRadius)} m) and the sea: raise the hub.`);
  const surface = { ...SURF, R: rotor.rotorRadius, area: rotor.sweptArea };
  const schedule = run(() => operatingSchedule({ surface, ratedPower: P, ratedRotorSpeed: rotor.ratedRotorSpeed, minRotorSpeed: rotor.minRotorSpeed,
    efficiency: rotor.drivetrainEfficiency, gearboxRatio: rotor.gearboxRatio, rho: rhoAir, cutIn: rotor.cutIn, cutOut: rotor.cutOut }), 'Operating schedule');
  const controller = run(() => controllerFromRotor({ rotor, surface, rho: rhoAir }), 'Controller design');
  // Parked, feathered rotor: thrust coefficient referred to the swept area from a BEM run at standstill.
  const parked = run(() => bemSolve({ rotor: ref, U: 50, omega: 0, pitchDeg: 90, rho: rhoAir }), 'Parked rotor');
  const iMax = schedule.thrust.indexOf(Math.max(...schedule.thrust));
  const out = {
    rotor, surface, schedule, controller, scale: s, R: rotor.rotorRadius, area: rotor.sweptArea, rhoAir,
    rnaMass: inp.rnaT > 0 ? inp.rnaT * 1e3 : rotor.rnaMass, hubOffset: 2.4 * s, zTowerTop: inp.hubHeight - 2.4 * s,
    parkedCt: Math.max(parked.Ct, 0.002), thrustRated: schedule.thrust[iMax], windAtMaxThrust: schedule.U[iMax],
    rpmRated: (rotor.ratedRotorSpeed * 30) / Math.PI, rpmMin: (rotor.minRotorSpeed * 30) / Math.PI,
    name: scaledRotor ? `generic ${fmt(inp.ratedMW)}-MW rotor scaled from the NREL 5-MW reference` : 'NREL 5-MW reference rotor (generic polars)',
  };
  TCACHE.set(key, out);
  return out;
}

/** Steady operating point (thrust N, rotor speed rad/s, pitch deg, Ct, power W) interpolated from the schedule. */
export function operatingPoint(tb, V) {
  const s = tb.schedule;
  if (!(V >= s.U[0] && V <= s.U[s.U.length - 1])) return { thrust: 0, omega: 0, pitchDeg: 90, ct: 0, power: 0, parked: true };
  let i = 0;
  while (i < s.U.length - 2 && s.U[i + 1] < V) i++;
  const w = (V - s.U[i]) / (s.U[i + 1] - s.U[i]);
  const L = (a) => a[i] + w * (a[i + 1] - a[i]);
  return { thrust: L(s.thrust), omega: L(s.omega), pitchDeg: L(s.pitchDeg), ct: L(s.ct), power: L(s.power), parked: false };
}

/** Tubular segments (monopile to the interface, tapered tower above) for the beam models; SI, z = 0 at still water. */
export function segments(inp, tb, steel = { E: 210e9, rho: 7850 }) {
  need(tb.zTowerTop > inp.zInt + 10, 'The tower top must be at least 10 m above the interface level: raise the hub or lower the interface.');
  need(inp.towerDt <= inp.towerDb, 'The tower top diameter cannot exceed the base diameter.');
  return [
    { zBase: -inp.depth, zTop: inp.zInt, dBase: inp.D, dTop: inp.D, tBase: inp.t / 1e3, tTop: inp.t / 1e3, density: steel.rho, E: steel.E },
    { zBase: inp.zInt, zTop: tb.zTowerTop, dBase: inp.towerDb, dTop: inp.towerDt, tBase: inp.towerTb / 1e3, tTop: inp.towerTt / 1e3, density: steel.rho * (inp.towerMassFactor ?? 1.05), E: steel.E },
  ];
}

/** Foundation argument of the beam models from the shared boundary inputs (SI). */
export function foundationOf(inp) {
  if (inp.fnd === 'fixed') return 'fixed';
  const Kh = inp.Kh * 1e9, Kr = inp.Kr * 1e9, Khr = inp.Khr * 1e6;
  need(Kh * Kr - Khr * Khr > 0, `The foundation stiffness matrix is not positive definite: |K_LR| = ${fmt(Math.abs(inp.Khr))} MN must be below √(K_LL·K_RR) = ${fmt(Math.sqrt(Kh * Kr) / 1e6)} MN. Reduce the coupling term or raise K_LL / K_RR.`);
  return { Kh, Kr, Khr };
}

/** Steel mass [kg] of the primary structure: embedded pile + monopile to the interface + tower. */
export function steelMass(inp, tb, rho = 7850) {
  const mp = tube(inp.D, inp.t / 1e3).A * rho * ((inp.Lp ?? 0) + inp.depth + inp.zInt);
  const H = tb.zTowerTop - inp.zInt;
  const a = (s) => tube(inp.towerDb + (inp.towerDt - inp.towerDb) * s, (inp.towerTb + (inp.towerTt - inp.towerTb) * s) / 1e3).A;
  const mt = (rho * H * (a(0) + 4 * a(0.5) + a(1))) / 6; // Simpson: A(s) is quadratic for a linear taper
  return { pile: mp, tower: mt, total: mp + mt };
}
