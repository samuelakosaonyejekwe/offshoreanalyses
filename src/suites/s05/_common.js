// Shared helpers of the wharf, jetty and quay suite: unit handling, check and advice records, the common input
// groups (pile bent, piles, soil, reinforced-concrete section), the transverse bent frame model with its load
// cases, and plotting helpers. Engineering units are converted to SI inside each module (docs/CONTRACT.md §2.1);
// functions here take SI unless the argument is the raw input record `inp`.
import { wharfBentModel } from '../../engine/structural/builders.js';
import { analyzeFrame, memberDiagram } from '../../engine/structural/frame3d.js';
import { isoTubularCheck } from '../../engine/structural/tubular.js';
import { rcSection, concreteProperties } from '../../engine/structural/concrete.js';
import { createProfile, applyScour } from '../../engine/geotech/soilprofile.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (cradle-to-gate order of magnitude), not project-specific declarations.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONC = 0.3;       // t CO2e per m³ structural concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const CO2_ROCK = 0.02;      // t CO2e per t quarried and placed rock (0.01–0.05 with haul distance)
export const MGO_T_PER_DAY = 12;   // piling / crane barge spread fuel burn, t/day (screening)

export const CONFIRM = 'Overridable default: confirm against the edition of the standard in force on the project.';
export const NOT_HIFI = 'CFD wave loading, explicit nonlinear contact finite elements and three-dimensional soil continuum finite elements are outside this module.';
export const ECON_BASIS = 'AACE RP 18R-97 estimate classes; whole-life costing per ISO 15686-5';
export const RC_NOTE = 'Reinforced-concrete sections are analysed in uniaxial bending with the EN 1992-1-1 material models; the partial factors and αcc are inputs to be set to the governing national annex or code.';

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible structure.`);
  return v;
}
export const r = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const f = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
export const arr = (a, n = 4) => Array.from(a, (v) => r(v, n));

/** Piecewise-linear interpolation with end clamping (xs ascending). */
export function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + ((ys[hi] - ys[lo]) * (x - xs[lo])) / (xs[hi] - xs[lo]);
}

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 300) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round((i * (len - 1)) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}

/** Largest x in [lo, hi] for which ok(x) holds, assuming ok is true at lo and switches once (bisection). */
export function maxWhere(ok, lo, hi, iters = 40) {
  if (!ok(lo)) return lo;
  if (ok(hi)) return hi;
  let a = lo, b = hi;
  for (let i = 0; i < iters; i++) { const m = 0.5 * (a + b); if (ok(m)) a = m; else b = m; }
  return a;
}

/** Code-check record; utilisation is demand / capacity, kept finite and non-negative. */
export function check(name, demand, capacity, unit, ref, explain) {
  const d = Math.abs(demand);
  const uc = capacity > 0 ? d / capacity : (d > 0 ? 99 : 0);
  return { name, demand: r(d), capacity: r(capacity), unit, uc: r(cap(uc)), ref, explain };
}
/** Check expressed as a required safety factor: uc = required / achieved. */
export function sfCheck(name, sf, required, ref, explain) {
  const s = Number.isFinite(sf) ? Math.max(sf, 1e-6) : 1e6;
  return { name, demand: r(required), capacity: r(Math.min(s, 999)), unit: '-', uc: r(cap(required / s)), ref, explain };
}
export const ucCheck = (name, uc, ref, explain) => ({ name, demand: r(cap(Math.max(uc, 0))), capacity: 1, unit: '-', uc: r(cap(Math.max(uc, 0))), ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

export const eng = (severity, text, action, basis, sustainability) => ({ kind: 'engineering', severity, text, action, basis, ...(sustainability ? { sustainability } : {}) });
export const eco = (severity, text, action, basis, sustainability) => ({ kind: 'economic', severity, text, action, basis, ...(sustainability ? { sustainability } : {}) });
export const steelCo2 = (t) => `${f(Math.abs(t))} t of steel ≈ ${f(Math.abs(t) * CO2_STEEL)} t CO₂e at ${CO2_STEEL} t CO₂e/t (screening factor, not a product declaration)`;
export const concCo2 = (m3) => `${f(Math.abs(m3))} m³ of concrete ≈ ${f(Math.abs(m3) * CO2_CONC)} t CO₂e at ${CO2_CONC} t CO₂e/m³ (screening factor, not a product declaration)`;
export const spreadCo2 = (days) => `${f(days)} marine-spread day(s) ≈ ${f(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel at ${CO2_MGO} t CO₂ per t)`;

/* ------------------------------------------------------------------------------------------------ inputs */

/** Tubular steel pile of the wharf. */
export function pileInputs({ D = 1.067, t = 25, group = 'Piles', corr = 3 } = {}) {
  return [
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: D, min: 0.3, max: 3, group, role: 'geometry', design: true, bind: 'wharf.pileD' },
    { key: 't', label: 'Pile wall thickness (as built)', unit: 'mm', type: 'number', default: t, min: 8, max: 80, group, role: 'geometry', design: true, bind: 'wharf.pileT', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'corr', label: 'Corrosion allowance / measured wall loss', unit: 'mm', type: 'number', default: corr, min: 0, max: 30, group, role: 'ic',
      help: 'Deducted from the wall for stiffness and strength: the design allowance for a new wharf, the measured loss for an existing one.' },
    { key: 'steel', label: 'Pile steel grade', type: 'select', source: 'materials:steel', default: 's355', group, role: 'material' },
  ];
}

/** Section and material of the pile in SI, net of the corrosion allowance. */
export function pileProps(inp, ctx) {
  const m = ctx?.material ? ctx.material(inp.steel) : null;
  need(m && m.fy > 0, `Steel grade '${inp.steel}' is not in the materials library.`);
  const D = inp.D, tNom = inp.t / 1e3, t = tNom - (inp.corr ?? 0) / 1e3;
  need(t >= 0.004, `The wall remaining after ${f(inp.corr)} mm of corrosion is ${f(t * 1e3)} mm: below 4 mm the tube cannot be assessed as a structural member. Increase the wall or reduce the loss.`);
  need(t < D / 2, 'The pile wall thickness must be less than half the diameter.');
  const Di = D - 2 * t;
  const A = (Math.PI / 4) * (D * D - Di * Di), I = (Math.PI / 64) * (D ** 4 - Di ** 4);
  const Zp = (D ** 3 - Di ** 3) / 6;
  return { D, t, tNom, E: m.E, fy: m.fy, rho: m.rho ?? 7850, A, I, EI: m.E * I, Z: (2 * I) / D, Zp, Mp: Zp * m.fy, My: ((2 * I) / D) * m.fy, name: m.name, massPerM: Math.PI * (D - tNom) * tNom * 7850 };
}

/** Geometry of one transverse bent (deck beam on a row of piles) and of the deck it carries. */
export function bentInputs({ group = 'Wharf geometry' } = {}) {
  return [
    { key: 'frameGeom', label: 'Imported bent frame (optional)', type: 'geometry', kinds: ['frame'], default: null, group, role: 'geometry',
      help: 'Frame model of one transverse bent (SDNF, CIS/2, IFC or node–member JSON) in the x–z plane, x from the berthing face landward, z up from the still-water datum. Nodes at the top elevation are the deck; the lowest nodes are fixed. When null the bent is generated from the fields below.' },
    { key: 'nPiles', label: 'Piles per bent', unit: '-', type: 'int', default: 4, min: 1, max: 12, group, role: 'geometry', design: true },
    { key: 'sx', label: 'Pile spacing across the wharf', unit: 'm', type: 'number', default: 6, min: 1.5, max: 15, group, role: 'geometry', design: true },
    { key: 'sy', label: 'Bent spacing along the wharf', unit: 'm', type: 'number', default: 6, min: 1.5, max: 15, group, role: 'geometry', design: true, help: 'Tributary length of quay carried by one bent.' },
    { key: 'overhang', label: 'Deck overhang beyond the end piles', unit: 'm', type: 'number', default: 1.5, min: 0, max: 5, group, role: 'geometry' },
    { key: 'deckLevel', label: 'Deck beam axis level above datum', unit: 'm', type: 'number', default: 4.0, min: 0.5, max: 15, group, role: 'geometry', help: 'Datum is the still-water level used for buoyancy (z = 0).' },
    { key: 'depth', label: 'Water depth at the berthing face', unit: 'm', type: 'number', default: 12, min: 2, max: 30, group, role: 'bc', bind: 'site.depth', help: 'Dredged depth below datum at the seaward pile row.' },
    { key: 'bedRise', label: 'Rise of the seabed to the landward pile row', unit: 'm', type: 'number', default: 6, min: 0, max: 25, group, role: 'bc', help: 'Under-wharf slope: seabed level rises linearly from the berthing face to the landward row. 0 for a level bed.' },
    { key: 'fixDepth', label: 'Depth of apparent fixity below the seabed (0 = 5 D)', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group, role: 'bc', bind: 'wharf.fixityDepth',
      help: 'The piles are built-in at this depth. Take it from the soil–structure interaction module (equivalent cantilever of the pile in soil); about 4–6 D in stiff soil, 6–10 D in soft soil.' },
    { key: 'rakerBatter', label: 'Raker pile batter at the landward row (0 = none)', unit: '-', type: 'number', default: 0, min: 0, max: 12, group, role: 'geometry', design: true, help: '1 horizontal to this many vertical; a pair of opposed rakers is added at the landward pile head. 0 leaves a frame of vertical piles.' },
    { key: 'headFix', label: 'Pile-head connection', type: 'select', default: 'fixed', group, role: 'bc', options: [{ value: 'fixed', label: 'Moment connection into the deck beam' }, { value: 'pinned', label: 'Pinned (shear and axial only)' }],
      help: 'Structural continuity between pile and deck. A frame of pinned vertical piles relies on pile cantilever action alone.' },
    { key: 'capB', label: 'Deck (cap) beam width', unit: 'm', type: 'number', default: 1.2, min: 0.3, max: 4, group: 'Deck', role: 'geometry' },
    { key: 'capH', label: 'Deck (cap) beam depth', unit: 'm', type: 'number', default: 1.4, min: 0.3, max: 4, group: 'Deck', role: 'geometry', design: true },
    { key: 'slabT', label: 'Deck slab thickness', unit: 'm', type: 'number', default: 0.45, min: 0.15, max: 1.5, group: 'Deck', role: 'geometry', design: true },
    { key: 'conc', label: 'Deck concrete class', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Deck', role: 'material' },
    { key: 'kCr', label: 'Cracked stiffness of the deck beam (fraction of gross)', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Deck', role: 'material', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Flexural stiffness used in the frame analysis; 0.35–0.7 of the gross section for cracked beams. Calibrate against measured deck deflection where available.' },
  ];
}

/** Characteristic actions on one bent. */
export function bentLoadInputs({ group = 'Actions on the bent', fender = 600, moor = 400 } = {}) {
  return [
    { key: 'gSdl', label: 'Superimposed dead load (surfacing, services)', unit: 'kPa', type: 'number', default: 2.5, min: 0, max: 30, group, role: 'ic', help: 'Permanent load present before any variable action.' },
    { key: 'qLive', label: 'Uniform cargo / storage load on the deck', unit: 'kPa', type: 'number', default: 30, min: 0, max: 200, group, role: 'bc', design: true, uncertain: { dist: 'gumbel', cov: 0.2 },
      help: 'For example 30 kPa general cargo, 40–60 kPa for containers stacked 3–4 high, 100+ kPa for heavy bulk.' },
    { key: 'crane', label: 'Crane rail and equipment loads on the bent', type: 'table', group, role: 'bc',
      columns: [{ key: 'name', label: 'Load', type: 'text' }, { key: 'x', label: 'x from the seaward pile', unit: 'm', type: 'number' }, { key: 'V', label: 'Vertical', unit: 'kN', type: 'number' }, { key: 'H', label: 'Horizontal (landward +)', unit: 'kN', type: 'number' }],
      default: [{ name: 'Seaside crane rail', x: 0, V: 1800, H: 90 }, { name: 'Landside crane rail', x: 18, V: 1400, H: 70 }],
      help: 'Rail loads per bent: wheel load × wheels carried by one bent (corner load spread by the rail beam). Take the values from the crane supplier for the operating and stowed conditions.' },
    { key: 'Ffender', label: 'Fender reaction on the bent (landward)', unit: 'kN', type: 'number', default: fender, min: 0, max: 20000, group, role: 'bc', bind: 'wharf.fenderReactionBent', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Share of the fender reaction carried by one bent after spreading through the deck diaphragm.' },
    { key: 'Fmoor', label: 'Mooring pull on the bent (seaward)', unit: 'kN', type: 'number', default: moor, min: 0, max: 20000, group, role: 'bc', bind: 'wharf.bollardLoadBent', help: 'Horizontal bollard load shared by one bent.' },
    { key: 'Fenv', label: 'Wave and current force on the piles of the bent', unit: 'kN', type: 'number', default: 120, min: 0, max: 20000, group, role: 'bc', bind: 'wharf.envForceBent', help: 'Resultant from the environmental loading module, applied as a uniform line load over the wetted length of the vertical piles (landward).' },
  ];
}

/** Partial factors of the limit-state combinations. */
export function factorInputs({ group = 'Partial factors' } = {}) {
  return [
    { key: 'gG', label: 'γG permanent actions', unit: '-', type: 'number', default: 1.35, min: 0.9, max: 1.6, group, role: 'criteria', help: `EN 1990 Table A1.2(B) STR set; BS 6349-2 uses the same framework. ${CONFIRM}` },
    { key: 'gQ', label: 'γQ cargo and crane actions', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group, role: 'criteria', help: CONFIRM },
    { key: 'gB', label: 'γ berthing and mooring actions', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group, role: 'criteria', help: `Applied to the characteristic (normal) fender reaction and mooring pull; use 1.0 with the abnormal-berthing reaction. ${CONFIRM}` },
    { key: 'psi', label: 'ψ0 combination factor for accompanying actions', unit: '-', type: 'number', default: 0.7, min: 0, max: 1, group, role: 'criteria', help: CONFIRM },
  ];
}

export const SOIL_DEFAULT = [
  { type: 'sand', thick: 6, gamma: 9, su: 0, suGrad: 0, phi: 32, Dr: 0.5, eps50: 0 },
  { type: 'clay', thick: 8, gamma: 8, su: 60, suGrad: 3, phi: 0, Dr: 0, eps50: 0.01 },
  { type: 'sand', thick: 30, gamma: 10, su: 0, suGrad: 0, phi: 37, Dr: 0.75, eps50: 0 },
];

/** Layered soil below the seabed at the pile. */
export function soilInputs({ group = 'Soil', layers = SOIL_DEFAULT, Lemb = 28, scour = 1 } = {}) {
  return [
    { key: 'layers', label: 'Soil layers (from the seabed down)', type: 'table', group, role: 'soil',
      columns: [{ key: 'type', label: 'Type (clay / sand)', type: 'text' }, { key: 'thick', label: 'Thickness', unit: 'm', type: 'number' },
        { key: 'gamma', label: 'γ′', unit: 'kN/m³', type: 'number' }, { key: 'su', label: 'su at top', unit: 'kPa', type: 'number' },
        { key: 'suGrad', label: 'su gradient', unit: 'kPa/m', type: 'number' }, { key: 'phi', label: 'φ′', unit: 'deg', type: 'number' },
        { key: 'Dr', label: 'Relative density', unit: '-', type: 'number' }, { key: 'eps50', label: 'ε50', unit: '-', type: 'number' }],
      default: layers,
      help: 'From the borehole logs. Clay rows use su, its gradient and ε50 (Matlock p–y, API α method); sand rows use φ′ and relative density (API sand p–y and β method). The last layer is extended to the pile toe if needed.' },
    { key: 'Lemb', label: 'Pile penetration below the seabed', unit: 'm', type: 'number', default: Lemb, min: 3, max: 90, group, role: 'geometry', design: true, bind: 'wharf.pileLemb' },
    { key: 'scour', label: 'Scour depth at the pile', unit: 'm', type: 'number', default: scour, min: 0, max: 12, group, role: 'ic', bind: 'wharf.scourDepth',
      help: 'Soil support removed around the pile (propeller wash, currents); from the scour module or the latest bathymetric survey.' },
  ];
}

/** Soil profile (engine record, SI) from the layer table, extended to the pile toe, with local scour applied. */
export function profileFrom(inp, D) {
  const rows = inp.layers || [];
  need(rows.length >= 1, 'Give at least one soil layer.');
  let top = 0;
  const layers = rows.map((q, i) => {
    const type = String(q.type || '').trim().toLowerCase();
    need(type === 'clay' || type === 'sand', `Soil layer ${i + 1}: type must be "clay" or "sand".`);
    need(q.thick > 0 && q.gamma > 0, `Soil layer ${i + 1}: thickness and γ′ must be positive.`);
    const l = { type, top, thickness: Number(q.thick), gammaSub: q.gamma * 1e3 };
    if (type === 'clay') {
      need(q.su > 0, `Soil layer ${i + 1}: clay needs su > 0.`);
      l.suTop = q.su * 1e3; l.suGrad = (q.suGrad || 0) * 1e3; l.eps50 = q.eps50 > 0 ? q.eps50 : 0.01;
    } else {
      need(q.phi >= 29 && q.phi <= 45, `Soil layer ${i + 1}: sand needs φ′ between 29° and 45° (range of the API p–y stiffness correlation).`);
      need(q.Dr >= 0.35 && q.Dr <= 1, `Soil layer ${i + 1}: sand needs relative density between 0.35 and 1 (API RP 2GEO has no design parameters for loose sand).`);
      l.phiDeg = Number(q.phi); l.Dr = Number(q.Dr);
    }
    top += Number(q.thick);
    return l;
  });
  const needDepth = inp.Lemb + 2;
  if (top < needDepth) layers[layers.length - 1].thickness += needDepth - top;
  let profile = run(() => createProfile(layers), 'Soil profile');
  const sc = inp.scour ?? 0;
  need(sc < inp.Lemb * 0.6, `Scour of ${f(sc)} m removes more than 60 % of the ${f(inp.Lemb)} m penetration: the pile is not a credible foundation.`);
  if (sc > 0) profile = run(() => applyScour(profile, { general: 0, local: sc, D, reductionDepth: Math.max(6 * D, sc) }), 'Scour');
  return profile;
}

/** Reinforced-concrete rectangular member (deck beam or slab strip): bars and design factors. */
export function rcInputs({ group = 'Reinforcement', nBot = 10, dBot = 32, nTop = 12, dTop = 32, cover = 60, links = 16, linkLegs = 4, linkS = 200 } = {}) {
  return [
    { key: 'nBot', label: 'Bottom bars', unit: '-', type: 'int', default: nBot, min: 2, max: 60, group, role: 'geometry', design: true },
    { key: 'dBot', label: 'Bottom bar diameter', unit: 'mm', type: 'number', default: dBot, min: 8, max: 50, group, role: 'geometry', design: true },
    { key: 'nTop', label: 'Top bars', unit: '-', type: 'int', default: nTop, min: 2, max: 60, group, role: 'geometry' },
    { key: 'dTop', label: 'Top bar diameter', unit: 'mm', type: 'number', default: dTop, min: 8, max: 50, group, role: 'geometry' },
    { key: 'cover', label: 'Nominal cover to the links', unit: 'mm', type: 'number', default: cover, min: 20, max: 120, group, role: 'geometry', help: 'Marine exposure XS3 typically needs 50–75 mm depending on concrete and design life.' },
    { key: 'dLink', label: 'Link diameter', unit: 'mm', type: 'number', default: links, min: 0, max: 25, group, role: 'geometry', help: '0 for a slab without shear reinforcement.' },
    { key: 'linkLegs', label: 'Link legs across the width', unit: '-', type: 'int', default: linkLegs, min: 0, max: 20, group, role: 'geometry' },
    { key: 'linkS', label: 'Link spacing', unit: 'mm', type: 'number', default: linkS, min: 50, max: 600, group, role: 'geometry' },
    { key: 'rebar', label: 'Reinforcing steel', type: 'select', source: 'materials:steel', default: 'reinf-b500', group, role: 'material' },
    { key: 'gammaC', label: 'γc concrete', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Resistance factors', role: 'criteria', help: `EN 1992-1-1 §2.4.2.4. ${CONFIRM}` },
    { key: 'gammaS', label: 'γs reinforcement', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Resistance factors', role: 'criteria', help: CONFIRM },
    { key: 'alphaCC', label: 'αcc long-term strength coefficient', unit: '-', type: 'number', default: 0.85, min: 0.7, max: 1, group: 'Resistance factors', role: 'criteria', help: 'EN 1992-1-1 §3.1.6 recommends 1.0; several national annexes (and EN 1992-2) use 0.85. Confirm for the governing code.' },
  ];
}

/** RC rectangular section b × h [m] with the bar layers of rcInputs; returns section, material and derived data (SI). */
export function rcMember(inp, ctx, b, h) {
  const cm = ctx?.material ? ctx.material(inp.conc) : null, sm = ctx?.material ? ctx.material(inp.rebar) : null;
  need(cm && cm.fck > 0, `Concrete '${inp.conc}' is not in the materials library.`);
  need(sm && sm.fy > 0, `Reinforcing steel '${inp.rebar}' is not in the materials library.`);
  const c = inp.cover / 1e3, dl = (inp.dLink || 0) / 1e3, db = inp.dBot / 1e3, dt = inp.dTop / 1e3;
  const dEff = h - c - dl - db / 2, dTopDepth = c + dl + dt / 2;
  need(dEff > dTopDepth + 0.02, `The section depth of ${f(h * 1e3)} mm is too small for the cover and bars given.`);
  const AsB = inp.nBot * Math.PI * db * db / 4, AsT = inp.nTop * Math.PI * dt * dt / 4;
  const clear = (b - 2 * (c + dl) - inp.nBot * db) / Math.max(inp.nBot - 1, 1);
  need(clear >= Math.min(db, 0.02) - 1e-9 || inp.nBot > 30, `The ${inp.nBot} bottom bars of ${f(inp.dBot)} mm do not fit in a width of ${f(b * 1e3)} mm in one layer (clear spacing ${f(clear * 1e3)} mm). Reduce the number, use a larger bar or widen the member.`);
  const props = concreteProperties({ fck: cm.fck });
  const material = { fck: cm.fck, fyk: sm.fy, gammaC: inp.gammaC ?? 1.5, gammaS: inp.gammaS ?? 1.15, alphaCC: inp.alphaCC ?? 0.85 };
  // sagging section: bottom bars in tension; hogging section: the same member turned over
  const sag = rcSection({ shape: 'rect', b, h, bars: [{ depth: dTopDepth, area: AsT, diameter: dt }, { depth: dEff, area: AsB, diameter: db }] });
  const hog = rcSection({ shape: 'rect', b, h, bars: [{ depth: h - dEff, area: AsB, diameter: db }, { depth: h - dTopDepth, area: AsT, diameter: dt }] });
  const AswS = inp.dLink > 0 && inp.linkLegs > 0 ? inp.linkLegs * Math.PI * dl * dl / 4 / (inp.linkS / 1e3) : 0;
  return { b, h, d: dEff, dHog: h - dTopDepth, AsB, AsT, db, dt, cover: c + dl, material, props, sag, hog, AswS, fyk: sm.fy, fck: cm.fck, Es: sm.E ?? 200e9, spacing: (b - 2 * (c + dl) - db) / Math.max(inp.nBot - 1, 1), rhoC: cm.rho ?? 2500 };
}

/* ------------------------------------------------------------------------------------------- bent model */

/**
 * Frame model of one transverse bent in the x–z plane (x landward from the seaward pile row, z up from the
 * still-water datum): a cracked reinforced-concrete deck beam on steel tubular piles built-in at the apparent
 * fixity depth. Load cases (characteristic, per bent): G permanent, Q uniform cargo, C crane and equipment,
 * B fender reaction, M mooring pull, E wave and current on the piles.
 */
export function bentModel(inp, ctx, opt = {}) {
  const P = pileProps(inp, ctx);
  const cm = ctx?.material ? ctx.material(inp.conc) : null;
  need(cm && cm.fck > 0, `Concrete '${inp.conc}' is not in the materials library.`);
  const Ec = concreteProperties({ fck: cm.fck }).Ecm, b = inp.capB, h = inp.capH, kc = inp.kCr ?? 0.5;
  const concMat = { E: Ec, nu: 0.2, rho: 2500, fy: cm.fck };
  const steelMat = { E: P.E, nu: 0.3, rho: P.rho, fy: P.fy };
  const deckSection = { A: b * h, Iy: (kc * b * h ** 3) / 12, Iz: (h * b ** 3) / 12, J: (kc * Math.min(b, h) ** 3 * Math.max(b, h)) / 4.5, Asy: (5 / 6) * b * h, Asz: (5 / 6) * b * h };
  let model, meta;
  const g = inp.frameGeom;
  if (g && g.kind === 'frame' && g.nodes && g.members?.length) {
    const nodes = []; for (let i = 0; i + 2 < g.nodes.length; i += 3) nodes.push([g.nodes[i], 0, g.nodes[i + 2]]);
    need(nodes.length >= 2 && nodes.length <= 400, 'The imported bent must have between 2 and 400 nodes.');
    const zTop = Math.max(...nodes.map((p) => p[2])), zBot = Math.min(...nodes.map((p) => p[2]));
    need(zTop - zBot > 1, 'The imported bent has no height: give it in the x–z plane with z up.');
    const isDeck = (m) => nodes[m.a][2] > zTop - 0.05 && nodes[m.b][2] > zTop - 0.05;
    const members = g.members.map((m) => (isDeck(m)
      ? { a: m.a, b: m.b, section: deckSection, material: concMat, group: 'deck' }
      : { a: nodes[m.a][2] < nodes[m.b][2] ? m.a : m.b, b: nodes[m.a][2] < nodes[m.b][2] ? m.b : m.a, section: { type: 'tube', D: m.section?.D > 0 ? m.section.D : P.D, t: m.section?.t > 0 ? m.section.t : P.t }, material: steelMat, group: 'pile', flooded: true }));
    const deckNodes = nodes.map((p, i) => i).filter((i) => nodes[i][2] > zTop - 0.05).sort((p, q) => nodes[p][0] - nodes[q][0]);
    const toeNodes = nodes.map((p, i) => i).filter((i) => nodes[i][2] < zBot + 0.05);
    need(deckNodes.length >= 1 && members.some((m) => m.group === 'pile'), 'The imported bent needs deck nodes at the top elevation and at least one pile member below them.');
    const supports = [...toeNodes.map((n) => ({ node: n, fix: 'fixed' })), ...nodes.map((p, i) => ({ node: i, fix: [0, 1, 0, 1, 0, 1] })).filter((s) => !toeNodes.includes(s.node))];
    const used = new Set(members.flatMap((m) => [m.a, m.b]));
    need(nodes.every((p, i) => used.has(i)), 'The imported bent has nodes that no member connects.');
    model = { nodes, members, supports };
    const idx = (grp) => members.map((m, i) => (m.group === grp ? i : -1)).filter((i) => i >= 0);
    meta = { headNodes: deckNodes.filter((n) => members.some((m) => m.group === 'pile' && (m.a === n || m.b === n))), toeNodes, deckNodes, deckMembers: idx('deck'), pileMembers: idx('pile'), imported: true };
    need(meta.headNodes.length >= 1, 'No pile of the imported bent reaches the deck elevation.');
    inp = { ...inp, deckLevel: zTop };
  } else {
    const n = Math.round(inp.nPiles);
    need(!(inp.headFix === 'pinned' && n < 2 && !(inp.rakerBatter > 0)), 'A single pile with a pinned head is a mechanism: use a moment connection or add piles.');
    const pileX = Array.from({ length: n }, (_, i) => i * inp.sx);
    const seabedZ = pileX.map((x, i) => -inp.depth + (n > 1 ? (inp.bedRise * i) / (n - 1) : 0));
    need(Math.max(...seabedZ) < inp.deckLevel - 0.5, 'The seabed at the landward row is at or above the deck: reduce the rise of the seabed.');
    const fix = (inp.fixDepth > 0 ? inp.fixDepth : 5 * P.D) + (opt.scour ?? 0);
    const rakers = inp.rakerBatter > 0 ? [{ pile: n - 1, batter: inp.rakerBatter, direction: 1 }, { pile: n - 1, batter: inp.rakerBatter, direction: -1 }] : [];
    const m0 = run(() => wharfBentModel({ pileX, deckElevation: inp.deckLevel, seabedZ, fixityDepth: fix, rakers, deckSection, pileSection: { type: 'tube', D: P.D, t: P.t }, overhang: inp.overhang, pinnedHeads: inp.headFix === 'pinned' && n >= 2, floodedPiles: true, material: steelMat }), 'Bent geometry');
    model = { nodes: m0.nodes, members: m0.members, supports: m0.supports };
    meta = { ...m0.meta, seabedZ, fix };
    for (const i of meta.deckMembers) model.members[i].material = concMat;
  }
  // deck stations sorted along x for placing loads
  const X = (nd) => model.nodes[nd][0];
  const deckM = meta.deckMembers.map((i) => { const m = model.members[i]; const xa = X(m.a), xb = X(m.b); return { i, x0: Math.min(xa, xb), x1: Math.max(xa, xb), flip: xa > xb }; }).sort((p, q) => p.x0 - q.x0);
  const xMin = deckM.length ? deckM[0].x0 : X(meta.headNodes[0]), xMax = deckM.length ? deckM[deckM.length - 1].x1 : xMin, x0 = X(meta.headNodes[0]);
  const sy = inp.sy;
  const udl = (q) => deckM.map((d) => ({ member: d.i, type: 'uniform', q: [0, 0, -q], axes: 'global' }));
  const points = (rows) => {
    const nodal = [], member = [];
    for (const [k, c] of rows.entries()) {
      const x = x0 + Number(c.x || 0), V = Number(c.V || 0) * 1e3, H = Number(c.H || 0) * 1e3;
      need(Number.isFinite(x) && Number.isFinite(V) && Number.isFinite(H), `Crane load row ${k + 1}: position and loads must be numbers.`);
      need(x >= xMin - 1e-6 && x <= xMax + 1e-6, `Crane load '${c.name || k + 1}' at x = ${f(c.x)} m is off the deck beam (it spans ${f(xMin - x0)} m to ${f(xMax - x0)} m from the seaward pile).`);
      const d = deckM.find((q) => x >= q.x0 - 1e-9 && x <= q.x1 + 1e-9);
      if (!d) { nodal.push({ node: meta.headNodes[0], F: [H, 0, -V, 0, 0, 0] }); continue; }
      const s = d.flip ? d.x1 - x : x - d.x0;
      member.push({ member: d.i, type: 'point', s: clamp(s, 0, d.x1 - d.x0), P: [H, 0, -V], axes: 'global' });
    }
    return { nodal, member };
  };
  const cr = points(inp.crane || []);
  // wave and current: uniform line load on each vertical pile between seabed and datum
  const envMember = [];
  const piles = meta.pileMembers.map((i) => model.members[i]).map((m, k) => ({ i: meta.pileMembers[k], m, vertical: Math.abs(model.nodes[m.a][0] - model.nodes[m.b][0]) < 1e-6 }));
  const vert = piles.filter((p) => p.vertical && Math.min(model.nodes[p.m.a][2], model.nodes[p.m.b][2]) < 0);
  const Fenv = (inp.Fenv || 0) * 1e3;
  if (Fenv > 0 && vert.length) {
    for (const p of vert) {
      const za = model.nodes[p.m.a][2], zb = model.nodes[p.m.b][2], zLo = Math.min(za, zb), zHi = Math.max(za, zb);
      const bed = meta.seabedZ ? Math.max(zLo + (meta.fix ?? 0), zLo) : zLo, top = Math.min(0, zHi);
      const wet = Math.max(top - bed, 0.1);
      const s1 = za < zb ? bed - za : za - top, s2 = s1 + wet;
      envMember.push({ member: p.i, type: 'uniform', q: [Fenv / vert.length / wet, 0, 0], axes: 'global', s1: Math.max(s1, 0), s2: Math.min(s2, zHi - zLo) });
    }
  }
  const seaNode = meta.deckNodes[0], headSea = meta.headNodes[0];
  const wSlab = (inp.slabT * 25e3 + (inp.gSdl || 0) * 1e3) * sy - 0; // slab and surfacing carried by the bent, N/m
  model.loadCases = [
    { name: 'G', selfWeight: { buoyancy: { waterline: 0 } }, member: udl(wSlab) },
    { name: 'Q', member: udl((inp.qLive || 0) * 1e3 * sy) },
    { name: 'C', nodal: cr.nodal, member: cr.member },
    { name: 'B', nodal: [{ node: headSea, F: [(inp.Ffender || 0) * 1e3, 0, 0, 0, 0, 0] }] },
    { name: 'M', nodal: [{ node: headSea, F: [-(inp.Fmoor || 0) * 1e3, 0, 0, 0, 0, 0] }] },
    { name: 'E', nodal: envMember.length ? [] : [{ node: headSea, F: [Fenv, 0, 0, 0, 0, 0] }], member: envMember },
  ];
  const deckLen = xMax - xMin;
  const pileLen = piles.reduce((s, p) => s + Math.hypot(model.nodes[p.m.a][0] - model.nodes[p.m.b][0], model.nodes[p.m.a][2] - model.nodes[p.m.b][2]), 0);
  return { model, meta, P, Ec, concMat, deckSection, deckLen, pileLen, nPileMembers: piles.length, seaNode, headSea, sy, x0 };
}

/** Standard limit-state combinations of the bent load cases. */
export function bentCombinations(inp) {
  const gG = inp.gG ?? 1.35, gQ = inp.gQ ?? 1.5, gB = inp.gB ?? 1.5, psi = inp.psi ?? 0.7;
  return [
    { name: 'ULS operating (cargo + crane)', factors: { G: gG, Q: gQ, C: gQ, E: psi * gQ }, ls: 'ULS' },
    { name: 'ULS berthing', factors: { G: gG, Q: psi * gQ, C: psi * gQ, B: gB }, ls: 'ULS' },
    { name: 'ULS mooring + environment', factors: { G: gG, Q: psi * gQ, C: psi * gQ, M: gB, E: gQ }, ls: 'ULS' },
    { name: 'ULS berthing, deck unloaded', factors: { G: 1.0, B: gB }, ls: 'ULS' },
    { name: 'SLS characteristic', factors: { G: 1, Q: 1, C: 1, B: psi, E: psi }, ls: 'SLS' },
    { name: 'SLS berthing', factors: { G: 1, Q: psi, C: psi, B: 1 }, ls: 'SLS' },
  ];
}

/**
 * Solve the bent once for its six characteristic load cases; `combo(factors)` then superposes them (linear
 * analysis) into pile forces with their ISO 19902 utilisation, the deck-beam moment and shear envelope, the
 * displacement vector and the vertical reaction. Frame convention for a beam with z up: My > 0 is hogging.
 */
export function solveBent(B, inp) {
  const out = run(() => analyzeFrame(B.model), 'Frame analysis');
  const K = inp.Kpile > 0 ? inp.Kpile : (inp.headFix === 'pinned' ? 2.1 : 1.2);
  const base = {};
  for (const c of B.model.loadCases) {
    const res = out.results[c.name];
    base[c.name] = { disp: res.displacements, piles: B.meta.pileMembers.map((i) => [res.memberForces[i].endA, res.memberForces[i].endB]),
      deck: B.meta.deckMembers.map((i) => { const d = memberDiagram(res, i, 13); return { s: d.s, My: d.My, Vz: d.Vz }; }), Rz: res.reactionSum[2], Rx: res.reactionSum[0] };
  }
  const geo = B.meta.pileMembers.map((i) => { const m = B.model.members[i]; return { L: Math.hypot(...[0, 1, 2].map((q) => B.model.nodes[m.a][q] - B.model.nodes[m.b][q])), D: m.section.D ?? B.P.D, t: m.section.t ?? B.P.t, group: m.group, xHead: B.model.nodes[m.b][0] }; });
  const combo = (factors) => {
    const fs = Object.entries(factors).filter(([c, v]) => v !== 0 && base[c]);
    const piles = geo.map((g, k) => {
      const end = [0, 1].map((e) => { const q = { N: 0, Vy: 0, Vz: 0, My: 0, Mz: 0 }; for (const [c, v] of fs) { const b = base[c].piles[k][e]; q.N += v * b.N; q.Vy += v * b.Vy; q.Vz += v * b.Vz; q.My += v * b.My; q.Mz += v * b.Mz; } return q; });
      const M = Math.max(...end.map((q) => Math.hypot(q.My, q.Mz))), V = Math.max(...end.map((q) => Math.hypot(q.Vy, q.Vz)));
      const N = Math.abs(end[0].N) > Math.abs(end[1].N) ? end[0].N : end[1].N; // tension positive
      const c = run(() => isoTubularCheck({ D: g.D, t: g.t, fy: B.P.fy, E: B.P.E, L: g.L, Ky: K, Kz: K, N, My: M, V }), 'Tubular member check');
      // a member beyond its Euler load returns an unbounded utilisation: report it at the cap so that results stay finite
      const uc = cap(c.governing?.uc ?? Math.max(...c.checks.map((q) => q.uc)));
      return { id: k + 1, group: g.group, x: g.xHead, N, M, V, Mhead: Math.hypot(end[1].My, end[1].Mz), L: g.L, uc, gov: c.governing?.name ?? '', ref: c.governing?.ref ?? 'ISO 19902:2020 §13' };
    });
    let Msag = 0, Mhog = 0, Vd = 0; const diag = [];
    for (let k = 0; k < B.meta.deckMembers.length; k++) {
      const n = base.G.deck[k].s.length, My = new Float64Array(n), Vz = new Float64Array(n);
      for (const [c, v] of fs) { const b = base[c].deck[k]; for (let i = 0; i < n; i++) { My[i] += v * b.My[i]; Vz[i] += v * b.Vz[i]; } }
      for (let i = 0; i < n; i++) { Msag = Math.max(Msag, -My[i]); Mhog = Math.max(Mhog, My[i]); Vd = Math.max(Vd, Math.abs(Vz[i])); }
      diag.push({ member: B.meta.deckMembers[k], s: base.G.deck[k].s, My, Vz });
    }
    const disp = new Float64Array(base.G.disp.length);
    for (const [c, v] of fs) { const d = base[c].disp; for (let i = 0; i < disp.length; i++) disp[i] += v * d[i]; }
    let Rz = 0, Rx = 0; for (const [c, v] of fs) { Rz += v * base[c].Rz; Rx += v * base[c].Rx; }
    return { piles, Msag, Mhog, Vd, diag, disp, Rz, Rx, ucMax: Math.max(...piles.map((p) => p.uc)), sway: Math.abs(disp[6 * B.headSea]) };
  };
  return { out, combo, K, base };
}

/** Displacement of a node in a result: [ux, uz] in metres. */
export const nodeDisp = (res, node) => [res.displacements[6 * node], res.displacements[6 * node + 2]];

/** Lines of a frame model for a model3d plot, with an optional deformed shape from a dof vector. */
export function framePlot(model, { disp, scale = 1, title, explain, values, valueLabel, valueUnit }) {
  const nodes = model.nodes.map((p) => p.map((v) => r(v, 5)));
  const lines = model.members.map((m) => [m.a, m.b]);
  const out = { type: 'model3d', title, nodes, lines, explain };
  if (disp) out.deformed = nodes.map((p, i) => [0, 1, 2].map((k) => r(p[k] + scale * disp[6 * i + k], 5)));
  if (values) { out.nodeValues = values.map((v) => r(v, 5)); out.valueLabel = valueLabel; out.valueUnit = valueUnit; }
  return out;
}

/** Scale that makes the largest displacement about 5 % of the model size. */
export function deformScale(model, disp) {
  let dmax = 0; for (let i = 0; i < model.nodes.length; i++) dmax = Math.max(dmax, Math.hypot(disp[6 * i], disp[6 * i + 1], disp[6 * i + 2]));
  const xs = model.nodes.map((p) => p[0]), zs = model.nodes.map((p) => p[2]);
  const size = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...zs) - Math.min(...zs));
  return dmax > 0 ? r((0.05 * size) / dmax, 2) : 1;
}

/** Steel and concrete quantities of one bent with its tributary deck, per bent. */
export function bentQuantities(B, inp) {
  const steel = (B.pileLen + B.nPileMembers * (inp.Lemb ?? 0)) * B.P.massPerM / 1e3;
  const concrete = B.deckLen * inp.capB * inp.capH + B.deckLen * inp.sy * inp.slabT;
  return { steel_t: r(steel), concrete_m3: r(concrete) };
}

/* ------------------------------------------------------------------------------------- single pile in soil */

/** Inputs of the single-pile lateral model shared by the foundation modules. */
export function pileSoilModelInputs({ group = 'Pile model', Lfree = 16 } = {}) {
  return [
    { key: 'Lfree', label: 'Free length, seabed to pile head', unit: 'm', type: 'number', default: Lfree, min: 0, max: 50, group, role: 'geometry', help: 'Water depth at the pile plus the height of the deck beam axis above the water line. Scour is added separately.' },
    { key: 'headFix', label: 'Pile-head connection', type: 'select', default: 'fixed', group, role: 'bc', options: [{ value: 'fixed', label: 'Moment connection into the deck (no head rotation)' }, { value: 'free', label: 'Pinned head (free to rotate)' }] },
    { key: 'loading', label: 'p–y loading type', type: 'select', default: 'cyclic', group, role: 'bc', options: [{ value: 'cyclic', label: 'Cyclic (berthing, waves)' }, { value: 'static', label: 'Static (sustained)' }] },
    { key: 'pMult', label: 'p-multiplier (group shadowing)', unit: '-', type: 'number', default: 1, min: 0.2, max: 1, group, role: 'soil', help: 'Reduction of lateral soil resistance for trailing rows at spacings below about 5 D (for example 0.8 / 0.4 / 0.3 at 3 D, 1.0 / 0.85 / 0.7 at 5 D for leading, second and third rows). 1 for rows more than 5 D apart.' },
    { key: 'nElem', label: 'Pile elements over the embedded length', unit: '-', type: 'int', default: 40, min: 10, max: 320, group: 'Numerics', role: 'control' },
  ];
}

/** Options of the engine's lateral pile solver from the module inputs (SI). */
export function latOpts(inp, P, profile, extra = {}) {
  return { L: inp.Lemb, stickUp: inp.Lfree ?? 0, nElem: Math.round(inp.nElem ?? 40), D: P.D, t: P.t, E: P.E, profile, loading: inp.loading ?? 'cyclic', pMult: inp.pMult ?? 1, ...extra };
}

/**
 * Lateral load at which the first plastic hinge forms in the pile, or at which the soil no longer provides
 * equilibrium, by load stepping on the nonlinear beam on p–y springs; the last step is interpolated on moment.
 * `solve` is the engine's lateral pile solver (passed in so that this file stays free of that import cycle).
 */
export function lateralLimit(solve, base, P, fixedHead) {
  const arm0 = (base.stickUp || 0) + 3 * P.D;
  const dH = ((fixedHead ? 2 : 1) * P.Mp) / arm0 / 12;
  const H = [0], y = [0], M = [0];
  let mode = 'plastic hinge in the pile', Hult = null;
  for (let k = 1; k <= 60; k++) {
    const h = k * dH;
    let s;
    try { s = solve({ ...base, head: fixedHead ? { H: h, fixedRotation: true } : { H: h, M: 0 } }); } catch { mode = 'soil failure (no lateral equilibrium)'; Hult = H[H.length - 1]; break; }
    const m = Math.abs(s.maxMoment);
    if (!(s.head.y > y[y.length - 1]) || !Number.isFinite(m)) { mode = 'soil failure (no lateral equilibrium)'; Hult = H[H.length - 1]; break; }
    if (m >= P.Mp) { const q = (P.Mp - M[M.length - 1]) / (m - M[M.length - 1]); Hult = H[H.length - 1] + q * (h - H[H.length - 1]); H.push(Hult); y.push(y[y.length - 1] + q * (s.head.y - y[y.length - 1])); M.push(P.Mp); break; }
    H.push(h); y.push(s.head.y); M.push(m);
  }
  if (Hult === null) { Hult = H[H.length - 1]; mode = 'load stepping limit reached without failure'; }
  return { Hult, mode, H, y, M };
}

/**
 * Copy of a bent model with every pile split into nSub elements for dynamic analysis (mass distribution along the
 * piles). Node numbers of the original model are kept; `pileSubs[k]` lists the element indices of original pile k.
 * A pinned head stays a release on the top element only; new nodes are held out of the plane of the bent.
 */
export function subdivideBent(model, meta, nSub) {
  const n = Math.max(1, Math.round(nSub));
  const nodes = model.nodes.map((p) => p.slice()), members = [], supports = model.supports.slice(), map = new Map();
  const pileSet = new Set(meta.pileMembers);
  model.members.forEach((m, i) => {
    if (!pileSet.has(i) || n === 1) { map.set(i, [members.length]); members.push({ ...m }); return; }
    const a = model.nodes[m.a], b = model.nodes[m.b], ids = [m.a];
    for (let k = 1; k < n; k++) { ids.push(nodes.length); nodes.push([0, 1, 2].map((q) => a[q] + ((b[q] - a[q]) * k) / n)); supports.push({ node: nodes.length - 1, fix: [0, 1, 0, 1, 0, 1] }); }
    ids.push(m.b);
    const list = [];
    for (let k = 0; k < n; k++) { const { releases, ...rest } = m; list.push(members.length); members.push({ ...rest, a: ids[k], b: ids[k + 1], ...(releases && k === n - 1 ? { releases } : {}) }); }
    map.set(i, list);
  });
  return { model: { ...model, nodes, members, supports }, pileSubs: meta.pileMembers.map((i) => map.get(i)), deckMembers: meta.deckMembers.map((i) => map.get(i)[0]) };
}
