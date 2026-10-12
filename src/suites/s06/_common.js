// Shared helpers of the dolphin suite: unit handling, common input groups (pile, soil, pile layout), the
// rigid-cap pile-group model on nonlinear pile-head stiffnesses, the equivalent-fixity frame model used by the
// buckling and seismic modules, and advice / check records. Engineering units are converted to SI inside each
// module (docs/CONTRACT.md §2.1); functions here take SI unless the argument is the raw input record `inp`.
import { createProfile, applyScour } from '../../engine/geotech/soilprofile.js';
import { solveLateralPile, pileHeadStiffness, equivalentCantilever, axialLoadSettlement } from '../../engine/geotech/pile.js';
import { analysePileGroup } from '../../engine/geotech/pilegroup.js';
import { isoTubularCheck } from '../../engine/structural/tubular.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
// Screening emission factors (cradle-to-gate order of magnitude), not project-specific declarations.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m³ concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 15;   // piling barge / jack-up spread fuel burn, t/day (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible dolphin.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
export const arr = (a) => Array.from(a, (v) => round(v, 5));

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
export const ucCheck = (name, uc, ref, explain) => ({ name, demand: round(uc, 4), capacity: 1, unit: '-', uc: Number.isFinite(uc) ? Math.min(Math.max(uc, 0), 99) : 99, ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (screening factor — confirm with the mill's product declaration).`;
export const vesselCo2Text = (days) => `${fmt(days)} marine-spread day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const CONFIRM = 'Overridable default — confirm against the edition of the standard in force on the project.';
export const NOT_HIFI = 'Explicit nonlinear contact finite elements, three-dimensional soil continuum finite elements and CFD-derived vessel loading are outside this module.';
export const ECON_BASIS = 'AACE RP 18R-97 estimate classes; whole-life costing per ISO 15686-5';

/* ------------------------------------------------------------------------------------------------ inputs */

/** Tubular steel pile: section, material, corrosion state, free and embedded lengths. */
export function pileInputs({ D = 1.4, t = 28, Lfree = 20, Lemb = 30, group = 'Pile' } = {}) {
  return [
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: D, min: 0.2, max: 6, group, role: 'geometry', design: true, bind: 'dolphin.pileD' },
    { key: 't', label: 'Pile wall thickness (as built)', unit: 'mm', type: 'number', default: t, min: 6, max: 150, group, role: 'geometry', design: true, bind: 'dolphin.pileT',
      uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'corr', label: 'Corrosion allowance / measured wall loss', unit: 'mm', type: 'number', default: 3, min: 0, max: 40, group, role: 'ic',
      help: 'Deducted from the wall for stiffness and strength. Use the design allowance for a new dolphin and the measured loss for an existing one.' },
    { key: 'steel', label: 'Pile steel grade', type: 'select', source: 'materials:steel', default: 's355', group, role: 'material' },
    { key: 'Lfree', label: 'Free length, mudline to pile head', unit: 'm', type: 'number', default: Lfree, min: 0, max: 80, group, role: 'geometry', bind: 'dolphin.Lfree',
      help: 'Water depth plus the height of the pile-cap soffit above the water line. Scour is added separately.' },
    { key: 'Lemb', label: 'Embedded length', unit: 'm', type: 'number', default: Lemb, min: 3, max: 120, group, role: 'geometry', design: true, bind: 'dolphin.Lemb' },
  ];
}

export const SOIL_DEFAULT = [
  { type: 'clay', thick: 6, gamma: 7, su: 30, suGrad: 2, phi: 0, Dr: 0, eps50: 0.015 },
  { type: 'sand', thick: 12, gamma: 10, su: 0, suGrad: 0, phi: 36, Dr: 0.7, eps50: 0 },
  { type: 'clay', thick: 40, gamma: 9, su: 120, suGrad: 3, phi: 0, Dr: 0, eps50: 0.007 },
];

/** Layered soil with p–y loading type, group p-multiplier and scour. */
export function soilInputs({ group = 'Soil', layers = SOIL_DEFAULT } = {}) {
  return [
    { key: 'layers', label: 'Soil layers (from mudline down)', type: 'table', group, role: 'soil',
      columns: [{ key: 'type', label: 'Type (clay / sand)', type: 'text' }, { key: 'thick', label: 'Thickness', unit: 'm', type: 'number' },
        { key: 'gamma', label: "γ′", unit: 'kN/m³', type: 'number' }, { key: 'su', label: 'su at top', unit: 'kPa', type: 'number' },
        { key: 'suGrad', label: 'su gradient', unit: 'kPa/m', type: 'number' }, { key: 'phi', label: "φ′", unit: 'deg', type: 'number' },
        { key: 'Dr', label: 'Relative density', unit: '-', type: 'number' }, { key: 'eps50', label: 'ε50', unit: '-', type: 'number' }],
      default: layers,
      help: 'Clay rows use su, its gradient and ε50 (Matlock p–y, API α method); sand rows use φ′ and relative density (API sand p–y, API β method). The last layer is extended to the pile toe if needed.' },
    { key: 'soilModel', label: 'Lateral soil model', type: 'select', default: 'py', group, role: 'bc',
      options: [{ value: 'py', label: 'Nonlinear p–y curves from the layers' }, { value: 'linear', label: 'Linear Winkler springs (uniform modulus)' }, { value: 'fixed', label: 'Built-in at an equivalent fixity depth' }],
      help: 'p–y curves for design. The linear model is for calibration against a lateral load test; the fixity model is the classical equivalent cantilever for quick checks.' },
    { key: 'kLin', label: 'Modulus of subgrade reaction', unit: 'kN/m²', type: 'number', default: 20000, min: 100, max: 5e6, group, role: 'soil', showIf: { key: 'soilModel', in: ['linear'] }, uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Soil reaction per unit pile length per unit deflection (already includes the pile width).' },
    { key: 'fixDepth', label: 'Depth of fixity below the mudline', unit: 'm', type: 'number', default: 6, min: 0.5, max: 40, group, role: 'bc', showIf: { key: 'soilModel', in: ['fixed'] }, help: 'About 4–6 D in stiff soil, 6–10 D in soft soil.' },
    { key: 'loading', label: 'Lateral soil response', type: 'select', default: 'static', group, role: 'soil', showIf: { key: 'soilModel', in: ['py'] },
      options: [{ value: 'static', label: 'Static p–y curves' }, { value: 'cyclic', label: 'Cyclic p–y curves (degraded)' }],
      help: 'Cyclic curves represent stiffness and strength degradation under repeated berthing and storm loading (API RP 2GEO §8.5).' },
    { key: 'pMult', label: 'Group p-multiplier', unit: '-', type: 'number', default: 1, min: 0.2, max: 1, group, role: 'soil', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Shadowing reduction of lateral soil resistance for closely spaced piles (about 0.8 / 0.4 / 0.3 by row at 3 D, 1.0 / 0.85 / 0.7 at 5 D in AASHTO LRFD Table 10.7.2.4-1). Overridable; 1 for an isolated pile.' },
    { key: 'scour', label: 'Local scour depth at the pile', unit: 'm', type: 'number', default: 1.5, min: 0, max: 15, group, role: 'ic', bind: 'dolphin.scourDepth',
      help: 'Soil support removed around the pile; overburden is reduced to 6 D below (API RP 2GEO scour guidance). About 1.5 D for an unprotected isolated pile.' },
  ];
}

export const PILES_DEFAULT = [
  { x: -2, y: -2, batter: 0, azimuth: 0 }, { x: 2, y: -2, batter: 0, azimuth: 0 },
  { x: 2, y: 2, batter: 0, azimuth: 0 }, { x: -2, y: 2, batter: 0, azimuth: 0 },
];

/** Pile-head layout under the cap; an imported point set replaces the table. */
export function layoutInputs({ group = 'Pile layout', piles = PILES_DEFAULT } = {}) {
  return [
    { key: 'layoutGeom', label: 'Imported pile-head layout (optional)', type: 'geometry', kinds: ['points', 'frame'], default: null, group, role: 'geometry',
      help: 'Point set (CSV / XYZ) or frame model (SDNF, CIS/2, node–member JSON): pile heads are the points, or the frame nodes at the highest elevation. When given it replaces the x, y columns of the table; rake is taken from the table row by row (0 if absent).' },
    { key: 'piles', label: 'Pile heads (plan position and rake)', type: 'table', group, role: 'geometry',
      columns: [{ key: 'x', label: 'x (along berthing line)', unit: 'm', type: 'number' }, { key: 'y', label: 'y (towards the vessel)', unit: 'm', type: 'number' },
        { key: 'batter', label: 'Rake from vertical', unit: 'deg', type: 'number' }, { key: 'azimuth', label: 'Rake azimuth from +x', unit: 'deg', type: 'number' }],
      default: piles, help: 'Origin at the cap centre, z up. A single row is a mono-pile dolphin.' },
    { key: 'headFix', label: 'Pile-head connection', type: 'select', default: 'fixed', group, role: 'bc',
      options: [{ value: 'fixed', label: 'Moment-fixed into the cap' }, { value: 'pinned', label: 'Pinned (shear and axial only)' }],
      help: 'Structural continuity between pile and cap. A pinned multi-pile group of vertical piles has no rotational restraint other than push–pull.' },
  ];
}

/** Pile-cap dimensions and sustained weight. */
export function capInputs({ group = 'Pile cap' } = {}) {
  return [
    { key: 'capB', label: 'Cap width (along berth)', unit: 'm', type: 'number', default: 6, min: 0.5, max: 40, group, role: 'geometry' },
    { key: 'capL', label: 'Cap length (normal to berth)', unit: 'm', type: 'number', default: 6, min: 0.5, max: 40, group, role: 'geometry' },
    { key: 'capT', label: 'Cap thickness', unit: 'm', type: 'number', default: 1.5, min: 0.2, max: 6, group, role: 'geometry', design: true },
    { key: 'gammaConc', label: 'Cap unit weight', unit: 'kN/m³', type: 'number', default: 25, min: 10, max: 80, group, role: 'material', help: '25 kN/m³ for reinforced concrete; use the steel deck weight divided by its envelope volume for a steel cap.' },
    { key: 'Wequip', label: 'Equipment and fender-panel weight', unit: 'kN', type: 'number', default: 250, min: 0, max: 1e5, group, role: 'ic', help: 'Sustained load on the cap: fender units and panel, hooks, capstans, walkway reactions.' },
  ];
}

/** Horizontal and vertical loads applied to the cap, with their partial factors. */
export function capLoadInputs({ group = 'Applied loads', Fy = 1575, Fx = 315 } = {}) {
  return [
    { key: 'Fy', label: 'Horizontal load normal to the berth', unit: 'kN', type: 'number', default: Fy, min: -1e5, max: 1e5, group, role: 'bc', bind: 'dolphin.fenderReaction', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Fender reaction (breasting dolphin) or line pull (mooring dolphin), unfactored.' },
    { key: 'Fx', label: 'Horizontal load along the berth', unit: 'kN', type: 'number', default: Fx, min: -1e5, max: 1e5, group, role: 'bc', help: 'Fender friction (μ·R) or the along-berth component of the line pull.' },
    { key: 'zLoad', label: 'Load height above the pile heads', unit: 'm', type: 'number', default: 1.5, min: -5, max: 20, group, role: 'bc', help: 'Fender centreline or hook elevation above pile-head level.' },
    { key: 'ecc', label: 'Plan eccentricity of the normal load', unit: 'm', type: 'number', default: 0.5, min: -20, max: 20, group, role: 'bc', help: 'Offset of the line of action from the cap centre along the berth; produces torsion.' },
    { key: 'Mz', label: 'Additional torsion', unit: 'kN·m', type: 'number', default: 0, min: -1e6, max: 1e6, group, role: 'bc' },
    { key: 'Fup', label: 'Vertical applied load (up +)', unit: 'kN', type: 'number', default: 0, min: -1e5, max: 1e5, group, role: 'bc', bind: 'dolphin.mooringUplift', help: 'Uplift component of mooring lines (positive) or live load on the cap (negative).' },
    { key: 'gammaF', label: 'Partial factor on applied loads', unit: '-', type: 'number', default: 1.2, min: 0.8, max: 2, group: 'Design factors', role: 'criteria', help: `Applied to berthing or mooring loads; member resistance factors are those of ISO 19902 §13. ${CONFIRM}` },
    { key: 'gammaG', label: 'Partial factor on self-weight', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.5, group: 'Design factors', role: 'criteria', help: `Use 0.9 when the weight relieves pile tension. ${CONFIRM}` },
  ];
}

/** Factored loads at the cap centre, pile-head level (kN → N). z up; Fz negative = down. W is the unfactored weight. */
export function capLoad(inp) {
  const g = inp.gammaF ?? 1, Fx = g * (inp.Fx ?? 0) * 1e3, Fy = g * (inp.Fy ?? 0) * 1e3;
  const W = (inp.capB * inp.capL * inp.capT * inp.gammaConc + inp.Wequip) * 1e3;
  return { Fx, Fy, Fz: -(inp.gammaG ?? 1) * W + g * (inp.Fup ?? 0) * 1e3, Mx: -Fy * (inp.zLoad ?? 0), My: Fx * (inp.zLoad ?? 0), Mz: Fy * (inp.ecc ?? 0) + g * (inp.Mz ?? 0) * 1e3, W };
}

/** Inputs of the equivalent-fixity frame model shared by the buckling and seismic modules. */
export function frameInputs({ group = 'Frame model' } = {}) {
  return [
    { key: 'LeUser', label: 'Pile length from head to fixity (0 = from the soil model)', unit: 'm', type: 'number', default: 0, min: 0, max: 150, group, role: 'bc', help: 'When zero, the length of the cantilever (fixed head: fixed–guided column) with the same lateral head stiffness as the pile in soil is used.' },
    { key: 'Href', label: 'Pile shear at which the soil stiffness is taken', unit: 'kN', type: 'number', default: 480, min: 1, max: 1e5, group, role: 'bc', bind: 'dolphin.pileShearMax', help: 'The secant soil stiffness, and with it the depth of fixity, depends on the load level.' },
    { key: 'shear', label: 'Include shear deformation of the piles', type: 'bool', default: true, group, role: 'control' },
    { key: 'nSub', label: 'Elements per pile', unit: '-', type: 'int', default: 4, min: 1, max: 32, group, role: 'control' },
  ];
}

/** Length head-to-fixity of the equivalent column, from the secant pile-head stiffness or the user value. */
export function fixityLength(inp, P, profile, fixed) {
  if (inp.LeUser > 0) return inp.LeUser;
  const H = inp.Href * 1e3;
  const hint = `The pile-head stiffness at the reference shear of ${fmt(inp.Href)} kN is not usable: the load is at or beyond the lateral capacity of the pile, or the embedment is too short for fixity. Reduce the reference shear, lengthen the pile or enter the length to fixity directly.`;
  let Le;
  try {
    const sec = pileHeadStiffness({ ...latOpts(inp, P, profile), head: fixed ? { H, fixedRotation: true } : { H, M: 0 } }).secant;
    const eq = equivalentCantilever({ KHH: sec.KHH, KHM: sec.KHM, KMM: sec.KMM, EI: P.EI });
    Le = fixed ? eq.fixedHead : eq.lateral;
  } catch { throw new Error(hint); }
  need(Number.isFinite(Le) && Le > 0.8 * inp.Lfree && Le < inp.Lfree + 1.2 * inp.Lemb, hint);
  return Le;
}

/* ------------------------------------------------------------------------------------------- data builders */

/** Section and material of the pile in SI, net of the corrosion allowance. */
export function pileProps(inp, ctx) {
  const m = ctx?.material ? ctx.material(inp.steel) : null;
  need(m && m.fy > 0, `Steel grade '${inp.steel}' is not in the materials library.`);
  const D = inp.D, tNom = inp.t / 1e3, t = tNom - (inp.corr ?? 0) / 1e3;
  need(t >= 0.004, `The wall remaining after ${fmt(inp.corr)} mm of corrosion is ${fmt(t * 1e3)} mm: below 4 mm the tube cannot be assessed as a structural member. Increase the wall or reduce the loss.`);
  need(t < D / 2, 'The wall thickness must be less than half the diameter.');
  const Di = D - 2 * t;
  const A = Math.PI / 4 * (D * D - Di * Di), I = Math.PI / 64 * (D ** 4 - Di ** 4);
  const Zp = (D ** 3 - Di ** 3) / 6;
  return { D, t, tNom, E: m.E, fy: m.fy, rho: m.rho ?? 7850, A, I, EI: m.E * I, Z: 2 * I / D, Zp, Mp: Zp * m.fy, My: 2 * I / D * m.fy, name: m.name };
}

/** Steel mass of the piles, t (as-built wall). */
export function pileSteelT(inp, n = 1) {
  const tN = inp.t / 1e3;
  return n * Math.PI * (inp.D - tN) * tN * (inp.Lfree + inp.Lemb) * 7.85;
}

/** Soil profile (engine record, SI) from the layer table, extended to the pile toe, with local scour applied. */
export function profileFrom(inp, D) {
  const rows = inp.layers || [];
  need(rows.length >= 1, 'Give at least one soil layer.');
  let top = 0;
  const layers = rows.map((r, i) => {
    const type = String(r.type || '').trim().toLowerCase();
    need(type === 'clay' || type === 'sand', `Soil layer ${i + 1}: type must be "clay" or "sand".`);
    need(r.thick > 0 && r.gamma > 0, `Soil layer ${i + 1}: thickness and γ′ must be positive.`);
    const l = { type, top, thickness: r.thick, gammaSub: r.gamma * 1e3 };
    if (type === 'clay') {
      need(r.su > 0, `Soil layer ${i + 1}: clay needs su > 0.`);
      l.suTop = r.su * 1e3; l.suGrad = (r.suGrad || 0) * 1e3; l.eps50 = r.eps50 > 0 ? r.eps50 : 0.01;
    } else {
      need(r.phi >= 29 && r.phi <= 45, `Soil layer ${i + 1}: sand needs φ′ between 29° and 45° (range of the API p–y stiffness correlation).`);
      need(r.Dr >= 0.35 && r.Dr <= 1, `Soil layer ${i + 1}: sand needs relative density between 0.35 and 1 (API RP 2GEO Table 1 has no design parameters for loose sand).`);
      l.phiDeg = r.phi; l.Dr = r.Dr;
    }
    top += r.thick;
    return l;
  });
  const needDepth = inp.Lemb + 2;
  if (top < needDepth) layers[layers.length - 1].thickness += needDepth - top;
  let profile = run(() => createProfile(layers), 'Soil profile');
  const sc = inp.scour ?? 0;
  need(sc < inp.Lemb * 0.6, `Scour of ${fmt(sc)} m removes more than 60 % of the ${fmt(inp.Lemb)} m embedment: the pile is not a credible foundation.`);
  if (sc > 0) profile = run(() => applyScour(profile, { general: 0, local: sc, D, reductionDepth: Math.max(6 * D, sc) }), 'Scour');
  return profile;
}

/** Options of the lateral pile solver for the pile of this dolphin (head loads added by the caller). */
export function latOpts(inp, P, profile, extra = {}) {
  const o = { L: inp.Lemb, stickUp: inp.Lfree, nElem: inp.nElem ?? 40, D: P.D, t: P.t, E: P.E, profile, loading: inp.loading ?? 'static', pMult: inp.pMult ?? 1, ...extra };
  const model = inp.soilModel ?? 'py';
  if (model === 'linear') { const k = inp.kLin * 1e3 * (inp.pMult ?? 1), sc = inp.scour ?? 0; o.pyAt = (z) => (z >= sc ? k : 0); }
  if (model === 'fixed') {
    need(inp.fixDepth < inp.Lemb, `The depth of fixity (${fmt(inp.fixDepth)} m) must be less than the embedded length (${fmt(inp.Lemb)} m).`);
    // very stiff lateral springs at every node below the fixity depth hold the pile there (the clamp sharpens as the
    // element size is reduced); no support above it
    const kBig = 1e5 * P.EI, zf = inp.fixDepth + (inp.scour ?? 0);
    o.pyAt = (z) => (z >= zf - 1e-9 ? kBig : 0);
  }
  return o;
}

/** Pile heads from the table or an imported point / frame model. */
export function pilesFrom(inp) {
  const rows = (inp.piles || []).map((r) => ({ x: +r.x || 0, y: +r.y || 0, batter: +r.batter || 0, azimuth: +r.azimuth || 0 }));
  const g = inp.layoutGeom;
  let out = rows;
  if (g) {
    let pts = [];
    if (g.kind === 'points' && g.positions) for (let i = 0; i + 2 < g.positions.length; i += 3) pts.push([g.positions[i], g.positions[i + 1], g.positions[i + 2]]);
    else if (g.kind === 'frame' && g.nodes) {
      const all = []; for (let i = 0; i + 2 < g.nodes.length; i += 3) all.push([g.nodes[i], g.nodes[i + 1], g.nodes[i + 2]]);
      const zTop = Math.max(...all.map((p) => p[2]));
      pts = all.filter((p) => p[2] > zTop - 0.05);
    }
    need(pts.length >= 1 && pts.length <= 60, 'The imported layout must hold between 1 and 60 pile heads (points, or frame nodes at the top elevation).');
    const cx = pts.reduce((s, p) => s + p[0], 0) / pts.length, cy = pts.reduce((s, p) => s + p[1], 0) / pts.length;
    out = pts.map((p, i) => ({ x: p[0] - cx, y: p[1] - cy, batter: rows[i]?.batter ?? 0, azimuth: rows[i]?.azimuth ?? 0 }));
  }
  need(out.length >= 1 && out.length <= 60, 'Give between 1 and 60 piles.');
  for (const [i, p] of out.entries()) need(p.batter >= 0 && p.batter <= 30, `Pile ${i + 1}: rake must be between 0° and 30° from vertical.`);
  for (let i = 0; i < out.length; i++) for (let j = i + 1; j < out.length; j++) need(Math.hypot(out[i].x - out[j].x, out[i].y - out[j].y) > 1e-6, `Piles ${i + 1} and ${j + 1} share the same head position.`);
  return out;
}

/** Smallest centre-to-centre spacing of the group, m (Infinity for a single pile). */
export function minSpacing(piles) {
  let s = Infinity;
  for (let i = 0; i < piles.length; i++) for (let j = i + 1; j < piles.length; j++) s = Math.min(s, Math.hypot(piles[i].x - piles[j].x, piles[i].y - piles[j].y));
  return s;
}

/**
 * Rigid-cap pile group on secant pile-head stiffnesses from the nonlinear BNWF pile.
 * The secant stiffness is taken at the largest pile shear and updated twice, so the softening of the soil under
 * the working load is carried into the group (Poulos & Davis 1980 ch. 9 stiffness method with nonlinear heads).
 * load: { Fx, Fy, Fz (down negative), Mx, My, Mz } at the cap centre at pile-head level, SI.
 */
export function analyseDolphin(inp, P, profile, load, opt = {}) {
  const piles = opt.piles ?? pilesFrom(inp);
  const n = piles.length, fixed = (inp.headFix ?? 'fixed') === 'fixed' && !opt.forceFree;
  const Hh = Math.hypot(load.Fx || 0, load.Fy || 0);
  // A single pile carries the cap moment itself: head is free, with the moment applied.
  const single = n === 1;
  const base = latOpts(inp, P, profile, opt.lat);
  // axial head stiffness: soil t–z / Q–z secant in series with the free length
  const ax = run(() => axialLoadSettlement({ profile, L: inp.Lemb, D: P.D, t: P.t, E: P.E, nElem: 30 }, { steps: 30 }), 'Axial pile response');
  const Kax = 1 / (1 / ax.headStiffness + inp.Lfree / (P.E * P.A));
  let Href = Math.max(Hh / n, 1e3), grp = null, sec = null;
  if (single) {
    const Mh = Math.hypot(load.Mx || 0, load.My || 0);
    const sol = run(() => solveLateralPile({ ...base, axialLoad: Math.max(-(load.Fz || 0), 0), head: { H: Math.max(Hh, 1), M: Mh } }), 'Lateral pile analysis');
    sec = run(() => pileHeadStiffness({ ...base, head: { H: Math.max(Hh, 1e3), M: Mh } }), 'Pile-head stiffness').secant;
    const eq = equivalentCantilever({ KHH: sec.KHH, KHM: sec.KHM, KMM: sec.KMM, EI: P.EI });
    const rows = [{ id: 1, axial: -(load.Fz || 0), shear: Hh, moment: Math.abs(sol.maxMoment), headMoment: Mh, torsion: Math.abs(load.Mz || 0) }];
    return { n, piles, single, fixed: false, Kax, sec, eq, sol, rows, disp: { u: sol.head.y, rot: sol.head.theta, uz: (load.Fz || 0) / Kax, rz: 0 }, Kgroup: Hh > 0 ? Hh / Math.max(sol.head.y, 1e-12) : sec.KHH - sec.KHM * sec.KHM / sec.KMM, Le: eq.lateral, Kfac: 2, axCap: ax.capacity, mobil: sol.maxMobilisation };
  }
  for (let it = 0; it < 3; it++) {
    sec = run(() => pileHeadStiffness({ ...base, head: fixed ? { H: Href, fixedRotation: true } : { H: Href, M: 0 } }), 'Pile-head stiffness').secant;
    const Ktor = P.E / 2.6 * 2 * P.I / (inp.Lfree + 5 * P.D); // St-Venant torsion of the tube over the free length plus about five diameters
    const stiffness = { KHH: sec.KHH, KHM: Math.min(sec.KHM, 0), KMM: sec.KMM, Kaxial: Kax, Ktorsion: Ktor, head: fixed ? 'fixed' : 'pinned' };
    grp = run(() => analysePileGroup({ piles: piles.map((p, i) => ({ id: i + 1, x: p.x, y: p.y, z: 0, batterDeg: p.batter, azimuthDeg: p.azimuth, stiffness })), load }), 'Pile-group analysis');
    const Hmax = Math.max(...grp.piles.map((p) => p.shear));
    if (Math.abs(Hmax - Href) < 0.02 * Href || Hmax < 1e3) break;
    Href = Hmax;
  }
  const Hmax = Math.max(...grp.piles.map((p) => p.shear), 1);
  // Bending along the governing pile: fixed head → largest moment at the cap; pinned → below the mudline.
  const sol = run(() => solveLateralPile({ ...base, head: fixed ? { H: Hmax, fixedRotation: true } : { H: Hmax, M: 0 } }), 'Lateral pile analysis');
  const Mgov = Math.abs(sol.maxMoment);
  const rows = grp.piles.map((p) => ({ id: p.id, axial: p.axial, shear: p.shear, headMoment: p.moment, moment: Math.max(p.moment, Mgov * p.shear / Hmax), torsion: Math.abs(p.torsion) }));
  const eq = equivalentCantilever({ KHH: sec.KHH, KHM: sec.KHM, KMM: sec.KMM, EI: P.EI });
  const d = grp.displacement, u = Math.hypot(d.ux, d.uy);
  return { n, piles, single, fixed, Kax, sec, eq, sol, rows, grp, disp: { u, ux: d.ux, uy: d.uy, rot: Math.hypot(d.rx, d.ry), uz: d.uz, rz: d.rz }, Kgroup: Hh > 0 && u > 0 ? Hh / u : grp.K[0][0],
    Le: fixed ? eq.fixedHead : eq.lateral, Kfac: fixed ? 1 : 2, axCap: ax.capacity, mobil: sol.maxMobilisation };
}

/** Governing utilisation of an ISO tubular check, capped so that an unstable member (infinite amplification) reports 99. */
export function govUc(chk) {
  let g = chk.governing;
  if (!g || !Number.isFinite(g.uc)) g = chk.checks.reduce((a, c) => ((Number.isFinite(c.uc) ? c.uc : Infinity) > (Number.isFinite(a.uc) ? a.uc : Infinity) ? c : a), chk.checks[0]);
  return { uc: Number.isFinite(g.uc) ? Math.min(Math.max(g.uc, 0), 99) : 99, name: g.name ?? '', ref: g.ref ?? 'ISO 19902:2020 §13' };
}

/** ISO 19902 member checks of every pile of an analysed dolphin. Returns rows with uc and the governing check. */
export function pileChecks(res, P, { gammaScale = 1 } = {}) {
  return res.rows.map((r) => {
    const c = run(() => isoTubularCheck({ D: P.D, t: P.t, fy: P.fy, E: P.E, L: res.Le, Ky: res.Kfac, Kz: res.Kfac, Cmy: 0.85, Cmz: 0.85,
      N: -r.axial, My: r.moment * gammaScale, V: r.shear * gammaScale, T: r.torsion * gammaScale }), 'Tubular member check');
    const g = govUc(c);
    return { ...r, uc: g.uc, gov: g.name, ref: g.ref, sigma: Math.abs(r.axial) / P.A + r.moment / P.Z, warnings: c.warnings };
  });
}

/**
 * Frame model of the dolphin on equivalent-fixity cantilevers: each pile runs from a fixed base at the apparent
 * fixity depth to its head, the heads are tied by stiff cap beams to a centre node carrying the cap mass.
 * Le [m] is the pile length head-to-fixity; nSub elements per pile.
 */
export function frameModel({ piles, P, Le, nSub, capMass = 0, pileRho, fixedHead = true, fy, extraMass = 0 }) {
  const nodes = [], members = [], supports = [], heads = [];
  for (const p of piles) {
    const b = p.batter * DEG, a = p.azimuth * DEG;
    const dir = [Math.sin(b) * Math.cos(a), Math.sin(b) * Math.sin(a), -Math.cos(b)];
    const i0 = nodes.length;
    for (let k = 0; k <= nSub; k++) { const s = Le * k / nSub; nodes.push([p.x + dir[0] * s, p.y + dir[1] * s, dir[2] * s]); }
    for (let k = 0; k < nSub; k++) members.push({ a: i0 + k, b: i0 + k + 1, section: 'pile', material: 'steel', group: 'pile', ...(extraMass > 0 ? { extraMass } : {}) });
    supports.push({ node: i0 + nSub, fix: 'fixed' });
    heads.push(i0);
  }
  let centre = heads[0];
  if (piles.length > 1) {
    centre = nodes.length; nodes.push([0, 0, 0]);
    for (const h of heads) members.push({ a: centre, b: h, section: 'cap', material: 'rigid', group: 'cap', ...(fixedHead ? {} : { releases: { b: 'ry,rz' } }) });
  }
  // Cap beams: stiff relative to the piles without wrecking the conditioning of the stiffness matrix.
  const big = 1e3;
  const model = {
    nodes, members, supports,
    materials: { steel: { E: P.E, nu: 0.3, rho: pileRho ?? P.rho, fy: fy ?? P.fy }, rigid: { E: P.E, nu: 0.3, rho: 1e-6, fy: fy ?? P.fy } },
    sections: { pile: { type: 'tube', D: P.D, t: P.t }, cap: { A: P.A * big, Iy: P.I * big, Iz: P.I * big, J: 2 * P.I * big, Asy: P.A * big, Asz: P.A * big } },
    nodeMasses: capMass > 0 ? [{ node: centre, m: capMass }] : [],
  };
  return { model, heads, centre };
}

/** Lines of a frame model for a model3d plot, with an optional deformed shape from a dof vector. */
export function framePlot(model, { disp, scale = 1, title, explain, values, valueLabel, valueUnit }) {
  const nodes = model.nodes.map((p) => p.map((v) => round(v, 5)));
  const lines = model.members.map((m) => [m.a, m.b]);
  const out = { type: 'model3d', title, nodes, lines, explain };
  if (disp) out.deformed = nodes.map((p, i) => [0, 1, 2].map((k) => round(p[k] + scale * disp[6 * i + k], 5)));
  if (values) { out.nodeValues = values.map((v) => round(v, 5)); out.valueLabel = valueLabel; out.valueUnit = valueUnit; }
  return out;
}
