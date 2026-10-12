// Shared helpers of the floating-structures suite: unit handling, parametric hull generators (closed surface
// meshes with the keel at z = 0, x forward, y to port, origin amidships on the centreline), mass properties,
// the motion model (slender-member, box-barge or imported coefficients), environment forces, advice records
// and screening sustainability factors. Engineering units are converted to SI here (docs/CONTRACT.md §2.1).
import { boxMesh, cylinderMesh, mergeMeshes, meshHydrostatics } from '../../engine/hydro/hydrostatics.js';
import { buildFloater, floaterRAO, boxBargeRAO, naturalPeriods, waveSpectrum } from '../../engine/hydro/floater.js';
import { tendonStiffness } from '../../engine/lines/mooringsystem.js';
import { projectedAreaForce, cylinderDriftShortWave, maruoDrift } from '../../engine/hydro/drift.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
/** Runs an engine call and turns RangeErrors into user-readable Errors. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check the inputs for an inconsistent combination.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
export const geomspace = (a, b, n) => Array.from({ length: n }, (_, i) => a * (b / a) ** (i / (n - 1)));
export const rnd = (arr, n = 4) => Array.from(arr, (v) => round(v, n));
/** Decimates parallel arrays to about n points (keeps the last point). */
export function thin(arrs, n = 400) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const step = Math.ceil(len / n);
  return arrs.map((a) => { const o = []; for (let i = 0; i < len; i += step) o.push(a[i]); if ((len - 1) % step) o.push(a[len - 1]); return o; });
}
export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });
export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : 99;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? Math.min(uc, 99) : 99, ref, explain };
};
/** Check with a lower limit (value must be at least `required`): utilisation = required / value. */
export const checkMin = (name, value, required, unit, ref, explain) => {
  const uc = value > 0 ? required / value : 99;
  return { name, demand: round(required, 5), capacity: round(value, 5), unit, uc: Number.isFinite(uc) ? Math.min(Math.max(uc, 0), 99) : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc) => (uc > 1 ? 'fail' : uc > 0.9 ? 'warn' : 'ok');
export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (confirm with the mill's product declaration).`;
export const fuelCo2Text = (tonnesFuel) => `${fmt(tonnesFuel)} t of marine gas oil ≈ ${fmt(tonnesFuel * CO2_MGO)} t CO₂ at the screening factor of ${CO2_MGO} t CO₂ per tonne of fuel.`;

export const LIMIT_HYDRO = 'There is no potential-flow panel (BEM) solver in the suite: motions come from long-wave Froude–Krylov and Morison loading on slender members (valid for member diameter / wavelength < 0.2, no diffraction, frequency-independent strip added mass), from a closed-form box-barge model (heave, pitch and a one-degree-of-freedom roll with user gyradius and damping; no surge, sway or yaw) or from imported coefficients; radiation damping is a user fraction of critical.';
export const LIMIT_CFD = 'CFD and hydroelastic finite-element analysis are outside the suite.';

/* --------------------------------------------------------------------------------------------- hull geometry */

/** Closed, outward-oriented vertical prism over a convex plan polygon given counter-clockwise ([x, y] pairs). */
export function prismMesh(poly, zMin, zMax) {
  const n = poly.length; const pos = []; const idx = [];
  for (const z of [zMin, zMax]) for (const [x, y] of poly) pos.push(x, y, z);
  let cx = 0; let cy = 0;
  for (const [x, y] of poly) { cx += x / n; cy += y / n; }
  pos.push(cx, cy, zMin, cx, cy, zMax);
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    idx.push(i, j, n + j, i, n + j, n + i);
    idx.push(2 * n, j, i);
    idx.push(2 * n + 1, n + i, n + j);
  }
  return { kind: 'surface', name: 'prism', positions: new Float32Array(pos), indices: new Uint32Array(idx) };
}

/** Wall-sided ship-like plan: parallel mid-body with parabolic waterline ends over `taper`·L at bow and stern. */
function shipPlan(L, B, taper, nEnd = 10) {
  const a = Math.max(taper, 0.01) * L; const pts = [];
  const half = (x) => { const d = Math.abs(x) - (L / 2 - a); return d <= 0 ? B / 2 : (B / 2) * (1 - (d / a) ** 2); };
  const xs = [];
  for (let i = 0; i <= nEnd; i++) xs.push(-L / 2 + a * i / nEnd);
  for (let i = 1; i <= nEnd; i++) xs.push(L / 2 - a + a * i / nEnd);
  // starboard side from stern to bow, then port side back (counter-clockwise seen from above)
  for (const x of xs) { const y = half(x); if (y > 1e-9 || x === xs[0] || x === xs[xs.length - 1]) pts.push([x, -y]); }
  for (let i = xs.length - 2; i >= 1; i--) pts.push([xs[i], half(xs[i])]);
  return pts;
}

export const HULL_TYPES = [
  { value: 'semi', label: 'Semi-submersible (twin pontoon, four columns)' },
  { value: 'ship', label: 'Ship-shaped hull (FPSO / FSO)' },
  { value: 'barge', label: 'Box barge' },
  { value: 'spar', label: 'Spar (deep-draft cylinder)' },
  { value: 'tlp', label: 'Tension-leg platform (four columns, ring pontoon)' },
];
const DRAFT0 = { semi: 22, ship: 16, barge: 6, spar: 90, tlp: 28 };
const KG0 = { semi: 21, ship: 13, barge: 7, spar: 38, tlp: 36 };

/**
 * Hull definition inputs. `type` sets the default hull type; draft and KG defaults follow it. `only` restricts the
 * selectable types (e.g. ship-shaped modules).
 */
export function hullInputs({ type = 'semi', only, withMesh = true, draft, kg } = {}) {
  const opts = only ? HULL_TYPES.filter((o) => only.includes(o.value)) : HULL_TYPES;
  const si = (types) => ({ key: 'hullType', in: types });
  const g = 'Hull geometry';
  const out = [
    { key: 'hullType', label: 'Hull type', type: 'select', options: opts, default: type, group: g, role: 'geometry',
      help: 'Parametric hull used when no surface is imported. Floating wind platforms are represented by the semi-submersible or spar with the turbine in the mass and windage data.' },
  ];
  if (withMesh) out.push({ key: 'hull', label: 'Imported hull surface (closed)', type: 'geometry', kinds: ['surface'], default: null, group: g, role: 'geometry',
    help: 'Closed triangulated hull (STL, OBJ, GDF …): x forward, y to port, z up. The lowest point is moved to z = 0 (keel) and draft is measured from there. When empty the parametric hull below is meshed.' });
  out.push(
    { key: 'L', label: 'Length', unit: 'm', type: 'number', default: 260, min: 20, max: 500, group: g, role: 'geometry', showIf: si(['ship', 'barge']), design: true, help: 'Overall length of the ship-shaped hull or barge.' },
    { key: 'B', label: 'Beam', unit: 'm', type: 'number', default: 46, min: 5, max: 100, group: g, role: 'geometry', showIf: si(['ship', 'barge']), design: true, help: 'Moulded breadth.' },
    { key: 'Dhull', label: 'Hull depth (keel to deck)', unit: 'm', type: 'number', default: 27, min: 2, max: 60, group: g, role: 'geometry', showIf: si(['ship', 'barge']), help: 'Moulded depth to the main deck.' },
    { key: 'endTaper', label: 'Bow / stern taper length (fraction of L)', unit: '-', type: 'number', default: 0.12, min: 0.02, max: 0.45, group: g, role: 'geometry', showIf: si(['ship']), help: 'Each end tapers parabolically in plan over this share of the length; the sides are vertical (wall-sided).' },
    { key: 'Lpon', label: 'Pontoon length', unit: 'm', type: 'number', default: 110, min: 30, max: 200, group: g, role: 'geometry', showIf: si(['semi']), help: 'Length of each of the two pontoons.' },
    { key: 'bPon', label: 'Pontoon width', unit: 'm', type: 'number', default: 17, min: 4, max: 40, group: g, role: 'geometry', showIf: si(['semi']), design: true, help: 'Breadth of each rectangular pontoon.' },
    { key: 'hPon', label: 'Pontoon height', unit: 'm', type: 'number', default: 9, min: 2, max: 20, group: g, role: 'geometry', showIf: si(['semi']), design: true, help: 'Depth of the pontoon section.' },
    { key: 'Dcol', label: 'Column diameter', unit: 'm', type: 'number', default: 15, min: 4, max: 35, group: g, role: 'geometry', showIf: si(['semi']), design: true, help: 'Diameter of each of the four circular columns standing on the pontoons.' },
    { key: 'sColX', label: 'Column spacing, longitudinal', unit: 'm', type: 'number', default: 70, min: 20, max: 150, group: g, role: 'geometry', showIf: si(['semi']), design: true, help: 'Centre-to-centre distance of the columns along a pontoon.' },
    { key: 'sColY', label: 'Column spacing, transverse', unit: 'm', type: 'number', default: 70, min: 20, max: 120, group: g, role: 'geometry', showIf: si(['semi']), design: true, help: 'Centre-to-centre distance between the two pontoons.' },
    { key: 'Hsemi', label: 'Keel to top of columns', unit: 'm', type: 'number', default: 40, min: 15, max: 80, group: g, role: 'geometry', showIf: si(['semi']), help: 'Height of the watertight hull from pontoon bottom to the deck-box underside.' },
    { key: 'Dspar', label: 'Spar diameter', unit: 'm', type: 'number', default: 30, min: 5, max: 50, group: g, role: 'geometry', showIf: si(['spar']), design: true, help: 'Diameter of the cylindrical hull.' },
    { key: 'Hspar', label: 'Spar hull height', unit: 'm', type: 'number', default: 105, min: 30, max: 250, group: g, role: 'geometry', showIf: si(['spar']), help: 'Keel to top of the hard tank.' },
    { key: 'DcolT', label: 'Column diameter', unit: 'm', type: 'number', default: 18, min: 5, max: 35, group: g, role: 'geometry', showIf: si(['tlp']), design: true, help: 'Diameter of each of the four corner columns.' },
    { key: 'sColT', label: 'Column spacing (square)', unit: 'm', type: 'number', default: 60, min: 25, max: 120, group: g, role: 'geometry', showIf: si(['tlp']), design: true, help: 'Centre-to-centre spacing of the columns in both directions.' },
    { key: 'bPonT', label: 'Ring-pontoon width', unit: 'm', type: 'number', default: 9, min: 3, max: 20, group: g, role: 'geometry', showIf: si(['tlp']), help: 'Breadth of the rectangular pontoons joining the columns at the keel.' },
    { key: 'hPonT', label: 'Ring-pontoon height', unit: 'm', type: 'number', default: 8, min: 3, max: 20, group: g, role: 'geometry', showIf: si(['tlp']), help: 'Depth of the pontoon section.' },
    { key: 'HtlpHull', label: 'Keel to top of columns', unit: 'm', type: 'number', default: 52, min: 20, max: 90, group: g, role: 'geometry', showIf: si(['tlp']), help: 'Height of the watertight columns.' },
    { key: 'tendonPretensionFrac', label: 'Total tendon pretension (fraction of buoyancy)', unit: '-', type: 'number', default: 0.25, min: 0.05, max: 0.6, group: g, role: 'bc', showIf: si(['tlp']),
      help: 'Excess buoyancy carried by the tendons at the design draft: weight = buoyancy × (1 − this fraction). Typical 0.2–0.35.' },
    { key: 'draft', label: 'Draft', unit: 'm', type: 'number', default: draft ?? DRAFT0[type], min: 0.5, max: 240, group: 'Loading condition', role: 'ic', design: true,
      help: 'Draft at the keel on even keel. The displacement follows from the hull volume at this draft.' },
    { key: 'KG', label: 'Centre of gravity above keel (KG)', unit: 'm', type: 'number', default: kg ?? KG0[type], min: 0.2, max: 200, group: 'Loading condition', role: 'ic',
      uncertain: { dist: 'normal', cov: 0.03 }, help: 'Vertical centre of gravity of the loaded unit, from the weight report or inclining test, before free-surface correction.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1040, group: 'Loading condition', role: 'environment', bind: 'site.rhoW', help: 'Density of the water the unit floats in.' },
  );
  return out;
}

function parametricHull(inp, nTheta = 48) {
  const t = inp.hullType;
  if (t === 'barge') {
    return { mesh: boxMesh({ xMin: -inp.L / 2, xMax: inp.L / 2, yMin: -inp.B / 2, yMax: inp.B / 2, zMin: 0, zMax: inp.Dhull }), depth: inp.Dhull, Lc: inp.L, Bc: inp.B, Cwp: 1 };
  }
  if (t === 'ship') {
    const plan = shipPlan(inp.L, inp.B, inp.endTaper);
    return { mesh: prismMesh(plan, 0, inp.Dhull), depth: inp.Dhull, Lc: inp.L, Bc: inp.B, Cwp: 1 - (2 / 3) * inp.endTaper };
  }
  if (t === 'spar') return { mesh: cylinderMesh({ D: inp.Dspar, zMin: 0, zMax: inp.Hspar, nTheta }), depth: inp.Hspar, Lc: inp.Dspar, Bc: inp.Dspar };
  if (t === 'semi') {
    need(inp.Hsemi > inp.hPon + 1, 'The semi-submersible column top must be above the pontoon deck: increase the keel-to-top height or reduce the pontoon height.');
    need(inp.sColX + inp.Dcol <= inp.Lpon + 1e-9, `Columns overhang the pontoon ends: longitudinal spacing + column diameter (${fmt(inp.sColX + inp.Dcol)} m) exceeds the pontoon length (${fmt(inp.Lpon)} m).`);
    need(inp.sColY > Math.max(inp.bPon, inp.Dcol), 'The two pontoons overlap: increase the transverse column spacing or reduce the pontoon width.');
    const parts = [];
    for (const y of [inp.sColY / 2, -inp.sColY / 2]) {
      parts.push(boxMesh({ xMin: -inp.Lpon / 2, xMax: inp.Lpon / 2, yMin: y - inp.bPon / 2, yMax: y + inp.bPon / 2, zMin: 0, zMax: inp.hPon }));
      for (const x of [inp.sColX / 2, -inp.sColX / 2]) parts.push(cylinderMesh({ D: inp.Dcol, zMin: inp.hPon, zMax: inp.Hsemi, x0: x, y0: y, nTheta }));
    }
    return { mesh: mergeMeshes(parts), depth: inp.Hsemi, Lc: inp.Lpon, Bc: inp.sColY + Math.max(inp.bPon, inp.Dcol) };
  }
  if (t === 'tlp') {
    const s = inp.sColT; const D = inp.DcolT;
    need(s > D + 2, 'TLP columns overlap: increase the column spacing or reduce the diameter.');
    need(inp.HtlpHull > inp.hPonT + 1, 'The TLP column top must be above the pontoons.');
    const parts = [];
    for (const x of [s / 2, -s / 2]) for (const y of [s / 2, -s / 2]) parts.push(cylinderMesh({ D, zMin: 0, zMax: inp.HtlpHull, x0: x, y0: y, nTheta }));
    // pontoons span between column faces on the centrelines (the small crescent at each curved face is neglected)
    for (const c of [s / 2, -s / 2]) {
      parts.push(boxMesh({ xMin: -s / 2 + D / 2, xMax: s / 2 - D / 2, yMin: c - inp.bPonT / 2, yMax: c + inp.bPonT / 2, zMin: 0, zMax: inp.hPonT }));
      parts.push(boxMesh({ xMin: c - inp.bPonT / 2, xMax: c + inp.bPonT / 2, yMin: -s / 2 + D / 2, yMax: s / 2 - D / 2, zMin: 0, zMax: inp.hPonT }));
    }
    return { mesh: mergeMeshes(parts), depth: inp.HtlpHull, Lc: s + D, Bc: s + D };
  }
  throw new Error(`Unknown hull type '${t}'.`);
}

function importedHull(model) {
  need(model && model.positions && model.indices && model.indices.length >= 12, 'The imported hull is not a closed triangulated surface (positions and indices are required).');
  const p = Float32Array.from(model.positions);
  let zMin = Infinity; let zMax = -Infinity; let xMin = Infinity; let xMax = -Infinity; let yMin = Infinity; let yMax = -Infinity;
  for (let i = 0; i < p.length; i += 3) {
    need(Number.isFinite(p[i]) && Number.isFinite(p[i + 1]) && Number.isFinite(p[i + 2]), 'The imported hull contains non-finite coordinates.');
    xMin = Math.min(xMin, p[i]); xMax = Math.max(xMax, p[i]); yMin = Math.min(yMin, p[i + 1]); yMax = Math.max(yMax, p[i + 1]); zMin = Math.min(zMin, p[i + 2]); zMax = Math.max(zMax, p[i + 2]);
  }
  for (let i = 2; i < p.length; i += 3) p[i] -= zMin;
  return { mesh: { kind: 'surface', name: model.name || 'imported hull', positions: p, indices: Uint32Array.from(model.indices) }, depth: zMax - zMin, Lc: xMax - xMin, Bc: yMax - yMin, imported: true };
}

/**
 * Builds the hull (mesh, displacement at the design draft, mass and centre of gravity) from the form inputs.
 * @returns {{type, imported, mesh, depth, Lc, Bc, T, rho, hs, volume, mass (kg), cog, KG, pretension (N)}}
 */
export function hullFrom(inp, { nTheta = 48 } = {}) {
  const h = inp.hull ? importedHull(inp.hull) : parametricHull(inp, nTheta);
  const T = inp.draft; const rho = inp.rhoW ?? 1025;
  need(T < h.depth, `Draft ${fmt(T)} m is not below the hull depth ${fmt(h.depth)} m: the unit would be submerged. Reduce the draft.`);
  const type = h.imported ? 'imported' : inp.hullType;
  const pre = type === 'tlp' ? inp.tendonPretensionFrac : 0;
  let hs = run(() => meshHydrostatics({ mesh: h.mesh, draft: T, rho }), 'Hydrostatics of the hull');
  need(hs.volume > 0, 'The hull has no displaced volume at this draft. A negative volume means the imported surface is oriented inwards: flip the normals.');
  const mass = hs.displacement * (1 - pre);
  const cog = [hs.cobBody ? hs.cobBody[0] : 0, 0, inp.KG];
  // LCG is placed over the centre of buoyancy and TCG on the centreline: even keel by construction
  cog[1] = hs.cobBody ? hs.cobBody[1] : 0;
  hs = meshHydrostatics({ mesh: h.mesh, draft: T, rho, cog, mass });
  return { type, imported: !!h.imported, mesh: h.mesh, depth: h.depth, Lc: h.Lc, Bc: h.Bc, Cwp: h.Cwp, T, rho, hs, volume: hs.volume, mass, cog, KG: inp.KG, pretension: hs.displacement * pre * G };
}

/** model3d plot of a hull mesh; node values = hydrostatic pressure head at the given draft (m of water). */
export function hullPlot(mesh, draft, title = 'Hull model', explain) {
  const nodes = []; const vals = [];
  for (let i = 0; i < mesh.positions.length; i += 3) {
    nodes.push([round(mesh.positions[i], 5), round(mesh.positions[i + 1], 5), round(mesh.positions[i + 2], 5)]);
    vals.push(round(Math.max(0, draft - mesh.positions[i + 2]), 4));
  }
  const tris = [];
  for (let i = 0; i < mesh.indices.length; i += 3) tris.push([mesh.indices[i], mesh.indices[i + 1], mesh.indices[i + 2]]);
  return { type: 'model3d', title, nodes, tris, lines: [], nodeValues: vals, valueLabel: 'Hydrostatic head', valueUnit: 'm',
    explain: explain || 'The closed hull surface used for the hydrostatic integrals, coloured by the still-water pressure head (zero above the waterline). Check that the shape and the waterline are what you intended before trusting the numbers.' };
}

/* --------------------------------------------------------------------------------------------- environment */

export function seaInputs({ Hs = 6, Tp = 11, ret = '1', heading = 0, withDuration = true } = {}) {
  const g = 'Sea state';
  const out = [
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: Hs, min: 0, max: 25, group: g, role: 'environment', bind: `site.Hs${ret}`, uncertain: { dist: 'normal', cov: 0.1 }, help: 'Hs of the design or operating sea state.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: Tp, min: 3, max: 30, group: g, role: 'environment', bind: `site.Tp${ret}`, uncertain: { dist: 'normal', cov: 0.08 }, help: 'Spectral peak period.' },
    { key: 'gamma', label: 'JONSWAP peak factor γ', unit: '-', type: 'number', default: 2.5, min: 1, max: 7, group: g, role: 'environment', help: '1.0 gives a Pierson–Moskowitz spectrum; 3.3 is the North Sea mean; storms in fetch-limited seas are peakier.' },
    { key: 'headingDeg', label: 'Wave heading (propagation direction)', unit: 'deg', type: 'number', default: heading, min: 0, max: 360, group: g, role: 'environment', help: '0° = waves travel towards +x (following sea), 90° = beam sea towards port, 180° = head sea.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 350, min: 20, max: 4000, group: g, role: 'environment', bind: 'site.depth', help: 'Still-water depth at the site.' },
  ];
  if (withDuration) out.push({ key: 'durationH', label: 'Sea-state duration', unit: 'h', type: 'number', default: 3, min: 0.5, max: 12, group: g, role: 'environment', help: 'Exposure used for the most probable maximum (3 h is the customary storm peak).' });
  return out;
}

/* --------------------------------------------------------------------------------------------- motion model */

export function motionInputs({ nFreq = true, zeta = 0.04 } = {}) {
  const g = 'Mass distribution and damping';
  const out = [
    { key: 'kRoll', label: 'Roll gyradius / reference dimension', unit: '-', type: 'number', default: 0.42, min: 0.15, max: 0.8, group: g, role: 'ic', uncertain: { dist: 'normal', cov: 0.05 },
      help: 'Roll radius of gyration about the centre of gravity as a fraction of the overall beam (ship, barge, semi, TLP) or of the hull height (spar). For the box-barge model it includes added inertia (typical 0.35–0.45 B).' },
    { key: 'kPitch', label: 'Pitch gyradius / reference dimension', unit: '-', type: 'number', default: 0.3, min: 0.15, max: 0.8, group: g, role: 'ic',
      help: 'Pitch radius of gyration as a fraction of the overall length (ship, barge, semi, TLP) or of the hull height (spar). Typical 0.25 L for ships, 0.3–0.4 for column-stabilised units.' },
    { key: 'kYaw', label: 'Yaw gyradius / reference dimension', unit: '-', type: 'number', default: 0.35, min: 0.15, max: 0.8, group: g, role: 'ic', help: 'Yaw radius of gyration as a fraction of the overall length (diameter for a spar).' },
    { key: 'zetaRad', label: 'Radiation damping (fraction of critical)', unit: '-', type: 'number', default: zeta, min: 0.005, max: 0.3, group: g, role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'There is no radiation solver: wave-making damping is this fraction of critical in every restored degree of freedom (typical 0.01–0.10). Calibrate against a decay test or a diffraction analysis.' },
    { key: 'zetaRoll', label: 'Roll damping of the box model (fraction of critical)', unit: '-', type: 'number', default: 0.06, min: 0.01, max: 0.3, group: g, role: 'bc', showIf: { key: 'hullType', in: ['ship', 'barge'] },
      help: 'Total roll damping (wave-making plus bilge-keel and eddy damping) for the one-degree-of-freedom roll model; 0.02–0.10 is typical.' },
    { key: 'CdCol', label: 'Column drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group: g, role: 'bc', showIf: { key: 'hullType', in: ['semi', 'spar', 'tlp'] }, help: 'Morison drag coefficient of the circular columns (0.6–1.2).' },
    { key: 'CdPon', label: 'Pontoon drag coefficient', unit: '-', type: 'number', default: 2.0, min: 0.5, max: 3, group: g, role: 'bc', showIf: { key: 'hullType', in: ['semi', 'tlp'] }, help: 'Morison drag coefficient of the sharp-edged rectangular pontoons (1.5–2.5).' },
    { key: 'Kmoor', label: 'Horizontal station-keeping stiffness', unit: 'kN/m', type: 'number', default: 250, min: 0, max: 50000, group: 'Station keeping', role: 'bc', bind: 'mooring.Kxx',
      help: 'Linearised mooring (or DP) restoring stiffness in surge and sway. Yaw stiffness is taken as this value times the squared half-breadth. For a TLP the tendon stiffness is added automatically.' },
    { key: 'tendonEA', label: 'Tendon axial stiffness EA (each corner)', unit: 'MN', type: 'number', default: 30000, min: 1000, max: 300000, group: 'Station keeping', role: 'bc', showIf: { key: 'hullType', in: ['tlp'] },
      help: 'Summed EA of the tendons at one column. Tendons run from the keel to the seabed.' },
  ];
  if (nFreq) out.push({ key: 'nFreq', label: 'Frequency points', unit: '-', type: 'int', default: 100, min: 40, max: 1200, group: 'Numerical controls', role: 'control',
    help: 'Number of frequencies between 0.08 and 2.5 rad/s (geometric spacing). Narrow resonance peaks need a fine grid: run the refinement study.' });
  return out;
}

export const omegaGrid = (n) => geomspace(0.08, 2.5, Math.max(8, Math.round(n)));

/** Floater members in the floating position (origin on the still-water plane above the keel origin). */
function membersFrom(inp, hull) {
  const T = hull.T; const t = hull.type; const m = [];
  if (t === 'spar') m.push({ type: 'cylinder', p1: [0, 0, -T], p2: [0, 0, hull.depth - T], D: inp.Dspar, Cd: inp.CdCol, end1: 'half-disc' });
  else if (t === 'semi') {
    for (const y of [inp.sColY / 2, -inp.sColY / 2]) {
      m.push({ type: 'rect', p1: [-inp.Lpon / 2, y, -T + inp.hPon / 2], p2: [inp.Lpon / 2, y, -T + inp.hPon / 2], b: inp.bPon, h: inp.hPon, Cd: inp.CdPon });
      for (const x of [inp.sColX / 2, -inp.sColX / 2]) m.push({ type: 'cylinder', p1: [x, y, -T + inp.hPon], p2: [x, y, hull.depth - T], D: inp.Dcol, Cd: inp.CdCol });
    }
  } else if (t === 'tlp') {
    const s = inp.sColT; const D = inp.DcolT;
    for (const x of [s / 2, -s / 2]) for (const y of [s / 2, -s / 2]) m.push({ type: 'cylinder', p1: [x, y, -T], p2: [x, y, hull.depth - T], D, Cd: inp.CdCol, end1: 'half-disc' });
    for (const c of [s / 2, -s / 2]) {
      m.push({ type: 'rect', p1: [-s / 2 + D / 2, c, -T + inp.hPonT / 2], p2: [s / 2 - D / 2, c, -T + inp.hPonT / 2], b: inp.bPonT, h: inp.hPonT, Cd: inp.CdPon });
      m.push({ type: 'rect', p1: [c, -s / 2 + D / 2, -T + inp.hPonT / 2], p2: [c, s / 2 - D / 2, -T + inp.hPonT / 2], b: inp.bPonT, h: inp.hPonT, Cd: inp.CdPon });
    }
  }
  return m;
}

/**
 * Motion model of the unit. Slender-member model for semi, spar and TLP; closed-form box model for ship,
 * barge and imported hulls. Returns helpers that evaluate six-degree-of-freedom RAOs on a frequency grid.
 */
export function motionFrom(inp, hull) {
  const t = hull.type; const rho = hull.rho;
  const slender = t === 'semi' || t === 'spar' || t === 'tlp';
  const refRoll = t === 'spar' ? hull.depth : hull.Bc; const refPitch = t === 'spar' ? hull.depth : hull.Lc; const refYaw = hull.Lc;
  const gyr = [inp.kRoll * refRoll, inp.kPitch * refPitch, inp.kYaw * refYaw];
  const depth = inp.depth ?? 350;
  const Kh = (inp.Kmoor ?? 0) * 1e3;
  const Kext = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  Kext[0][0] = Kh; Kext[1][1] = Kh; Kext[5][5] = Kh * (Math.max(hull.Bc, hull.Lc) / 2) ** 2;
  let tendons = null;
  if (t === 'tlp') {
    const Lt = depth - hull.T;
    need(Lt > 10, `Water depth ${fmt(depth)} m leaves no room for tendons below a ${fmt(hull.T)} m draft.`);
    const s = inp.sColT;
    tendons = [];
    for (const x of [s / 2, -s / 2]) for (const y of [s / 2, -s / 2]) tendons.push({ x, y, L: Lt, EA: inp.tendonEA * 1e6, T0: hull.pretension / 4 });
    const ks = tendonStiffness(tendons);
    // tendons attach at the keel: the lateral geometric stiffness acts at z = −T, coupling surge–pitch and sway–roll
    Kext[0][0] += ks.surge; Kext[1][1] += ks.sway; Kext[2][2] += ks.heave; Kext[5][5] += ks.yaw;
    Kext[3][3] += ks.roll + ks.sway * hull.T * hull.T; Kext[4][4] += ks.pitch + ks.surge * hull.T * hull.T;
    Kext[0][4] = Kext[4][0] = -ks.surge * hull.T; Kext[1][3] = Kext[3][1] = ks.sway * hull.T;
    // pretension acts as a weight at the keel in the roll/pitch gravity stiffness
    Kext[3][3] += hull.pretension * hull.T; Kext[4][4] += hull.pretension * hull.T;
  }
  const cogW = [hull.cog[0], hull.cog[1], hull.KG - hull.T];
  if (slender) {
    const model = run(() => buildFloater({ members: membersFrom(inp, hull), mass: hull.mass, cog: cogW, gyradii: gyr, rho, depth, Kext, nStrips: 12 }), 'Slender-member floater model');
    const np = naturalPeriods(model);
    need(t === 'tlp' || (model.C[3][3] > 0 && model.C[4][4] > 0), `The unit is initially unstable (negative roll or pitch restoring): lower KG below ${fmt(hull.KG + Math.min(hull.hs.GMt, hull.hs.GMl))} m or change the draft.`);
    return {
      kind: 'slender', model, tendons, gyr, Kext, depth, validOmegaMax: model.validOmegaMax,
      Tn: np.uncoupled.map((v) => (Number.isFinite(v) && v > 0 ? v : null)),
      rao(omega, headingDeg, sea) {
        return run(() => floaterRAO({ model, omega, headingDeg, seaState: sea && sea.Hs > 0 ? sea : undefined, radiationDampingRatio: inp.zetaRad, viscous: !!(sea && sea.Hs > 0), tol: 2e-3, maxIter: 40 }), 'RAO solution');
      },
      // single-pass solve re-using the linearised drag damping of a previous sea-state solve (drag excitation omitted)
      raoLin(omega, headingDeg, viscousDamping) {
        return run(() => floaterRAO({ model, omega, headingDeg, radiationDampingRatio: inp.zetaRad, viscous: false, extraDamping: viscousDamping || undefined }), 'RAO solution');
      },
    };
  }
  const GMt = hull.hs.GMt;
  need(GMt > 0, `The unit is initially unstable (GM = ${fmt(GMt)} m): lower KG or increase the beam.`);
  const L = hull.Lc; const B = hull.Bc; const T = hull.T;
  const Cb = Math.min(1, hull.volume / (L * B * T));
  const kxx = inp.kRoll * B;
  const box0 = boxBargeRAO({ L, B, T, omega: [0.5], headingDeg: 90, Cb, GMt, kxx, rollDampingRatio: inp.zetaRoll ?? 0.06, rho });
  // pitch natural period from the same strip assumptions: added inertia taken equal to the rigid inertia (screening)
  const TnPitch = 2 * Math.PI * Math.sqrt(2 * (inp.kPitch * L) ** 2 / (G * Math.max(hull.hs.GMl, 1e-6)));
  const mTot = hull.mass;
  const TnSurge = Kh > 0 ? 2 * Math.PI * Math.sqrt(1.1 * mTot / Kh) : null;
  return {
    kind: 'box', model: null, tendons, gyr, Kext, depth, validOmegaMax: Infinity, Cb,
    Tn: [TnSurge, Kh > 0 ? 2 * Math.PI * Math.sqrt(1.8 * mTot / Kh) : null, box0.heaveNaturalPeriod, box0.rollNaturalPeriod, TnPitch, null],
    raoLin(omega, headingDeg) { return this.rao(omega, headingDeg); },
    rao(omega, headingDeg) {
      const r = run(() => boxBargeRAO({ L, B, T, omega, headingDeg, Cb, GMt, kxx, rollDampingRatio: inp.zetaRoll ?? 0.06, rho }), 'Box-barge RAO');
      const n = omega.length; const z = () => new Float64Array(n);
      const re = [z(), z(), Float64Array.from(r.heave), z(), Float64Array.from(r.pitch), z()]; const im = [z(), z(), z(), z(), z(), z()];
      for (let i = 0; i < n; i++) { const a = r.roll ? r.roll[i] : 0; const ph = r.rollPhase ? r.rollPhase[i] : 0; re[3][i] = a * Math.cos(ph); im[3][i] = a * Math.sin(ph); }
      const amp = re.map((c, d) => Float64Array.from(c, (v, i) => Math.hypot(v, im[d][i])));
      const phase = re.map((c, d) => Float64Array.from(c, (v, i) => Math.atan2(im[d][i], v)));
      return { omega: Float64Array.from(omega), period: Float64Array.from(omega, (w) => 2 * Math.PI / w), headingDeg, depth, g: G, re, im, amp, phase, force: null, viscousDamping: null, converged: true, iterations: 0, validOmegaMax: Infinity };
    },
  };
}

export const spectrum = (omega, sea) => (sea.Hs > 0 ? waveSpectrum(omega, sea) : new Float64Array(omega.length));
export const seaFrom = (inp) => ({ Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma ?? 2.5 });
export const DOF = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'];
export const DOF_UNIT = ['m', 'm', 'm', 'deg', 'deg', 'deg'];
export const dofScale = (d) => (d < 3 ? 1 : 1 / DEG);
export const trapz = (x, y) => { let s = 0; for (let i = 1; i < x.length; i++) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]); return s; };

/* --------------------------------------------------------------------------------------------- mean forces */

export function windCurrentInputs({ Uw = 36, Uc = 1.0, ret = '100' } = {}) {
  const g = 'Wind and current';
  return [
    { key: 'Uw', label: 'Wind speed (1-hour mean at 10 m)', unit: 'm/s', type: 'number', default: Uw, min: 0, max: 80, group: g, role: 'environment', bind: `site.Uw${ret}`, uncertain: { dist: 'normal', cov: 0.1 }, help: 'Mean wind used for the steady wind force.' },
    { key: 'Uc', label: 'Surface current speed', unit: 'm/s', type: 'number', default: Uc, min: 0, max: 4, group: g, role: 'environment', bind: `site.Uc${ret}`, uncertain: { dist: 'normal', cov: 0.15 }, help: 'Current speed taken as uniform over the draft.' },
    { key: 'AwFront', label: 'Windage area, head-on', unit: 'm²', type: 'number', default: 2400, min: 0, max: 30000, group: g, role: 'load', help: 'Projected area above water seen along the x-axis (hull freeboard, deck, topsides, turbine).' },
    { key: 'AwSide', label: 'Windage area, beam-on', unit: 'm²', type: 'number', default: 3200, min: 0, max: 40000, group: g, role: 'load', help: 'Projected area above water seen along the y-axis.' },
    { key: 'Cwind', label: 'Wind force coefficient', unit: '-', type: 'number', default: 1.0, min: 0.4, max: 2, group: g, role: 'load', help: 'Combined shape × height coefficient on the projected area. There is no built-in table: take it from wind-tunnel data or the MODU Code tables.' },
    { key: 'Ccur', label: 'Current drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 2.5, group: g, role: 'load', help: 'Drag coefficient on the projected underwater area. User input (model tests, OCIMF curves).' },
    { key: 'driftR', label: 'Wave-drift reflection coefficient', unit: '-', type: 'number', default: 0.7, min: 0, max: 1, group: g, role: 'load',
      help: 'Short-wave reflection coefficient for the mean wave-drift force on wall-sided hulls (Maruo). Columns use the short-wave cylinder limit. Screening only; replace with QTFs from a diffraction analysis.' },
  ];
}

/** Projected underwater areas (head-on, beam-on) in m². */
export function underwaterAreas(inp, hull) {
  const T = hull.T; const t = hull.type;
  if (t === 'semi') return { front: 2 * inp.bPon * inp.hPon + 2 * inp.Dcol * (T - inp.hPon), side: inp.Lpon * inp.hPon + 2 * inp.Dcol * Math.max(T - inp.hPon, 0) };
  if (t === 'tlp') return { front: 2 * inp.DcolT * T + (inp.sColT - inp.DcolT) * inp.hPonT, side: 2 * inp.DcolT * T + (inp.sColT - inp.DcolT) * inp.hPonT };
  if (t === 'spar') return { front: inp.Dspar * T, side: inp.Dspar * T };
  return { front: hull.Bc * T, side: hull.Lc * T };
}

/** Mean wave-drift coefficient D(ω) = F/a² (N/m²) seen from the given heading; screening model, see windCurrentInputs. */
export function driftCoefficient(inp, hull, omega, headingDeg) {
  const t = hull.type; const rho = hull.rho; const T = hull.T;
  const c = Math.abs(Math.cos(headingDeg * DEG)); const s = Math.abs(Math.sin(headingDeg * DEG));
  return Array.from(omega, (w) => {
    const k = w * w / G;
    // share of the deep-water wave energy flux above the keel: the part a wall-sided body can reflect
    const flux = 1 - Math.exp(-2 * k * T);
    if (t === 'semi') return 4 * cylinderDriftShortWave({ a: 1, D: inp.Dcol, rho }).force * (1 - Math.exp(-2 * k * Math.max(T - inp.hPon, 0)));
    if (t === 'tlp') return 4 * cylinderDriftShortWave({ a: 1, D: inp.DcolT, rho }).force * flux;
    if (t === 'spar') return cylinderDriftShortWave({ a: 1, D: inp.Dspar, rho }).force * flux;
    const width = hull.Bc * c + hull.Lc * s;
    return maruoDrift({ a: 1, R: (inp.driftR ?? 0.7) * Math.sqrt(flux), rho }) * width;
  });
}

/** Steady wind and current forces (N) travelling towards `dirDeg` (vessel frame). */
export function windCurrentForce(inp, hull, dirDeg, Uw = inp.Uw, Uc = inp.Uc) {
  const ua = underwaterAreas(inp, hull);
  const wind = projectedAreaForce({ V: Uw, directionDeg: dirDeg, Afront: inp.AwFront, Aside: inp.AwSide, Cfront: inp.Cwind, Cside: inp.Cwind, rho: 1.226 });
  const cur = projectedAreaForce({ V: Uc, directionDeg: dirDeg, Afront: ua.front, Aside: ua.side, Cfront: inp.Ccur, Cside: inp.Ccur, rho: hull.rho });
  return { wind, cur, areas: ua };
}

/** Hull steel mass estimate (t) from the displaced volume: structural density input in kg per m³ of enclosed hull volume. */
export function hullEnclosedVolume(hull) {
  const hs = meshHydrostatics({ mesh: hull.mesh, draft: hull.depth * 0.9999, rho: hull.rho });
  return hs.volume / 0.9999;
}

/* --------------------------------------------------------------------------------------------- statistics */

/** Short-term statistics of a linear response (amplitude per unit wave amplitude) in a spectrum; zeros in a flat calm. */
export function stats(omega, amp, S, duration = 10800) {
  const y = Array.from(amp, (a, i) => a * a * S[i]);
  const m0 = trapz(omega, y);
  if (!(m0 > 1e-300)) return { m0: 0, m2: 0, sigma: 0, significant: 0, Tz: 0, mpm: 0, nCycles: 0 };
  const m2 = trapz(omega, y.map((v, i) => v * omega[i] * omega[i]));
  const Tz = 2 * Math.PI * Math.sqrt(m0 / m2); const n = Math.max(duration / Tz, 2);
  const sigma = Math.sqrt(m0);
  return { m0, m2, sigma, significant: 2 * sigma, Tz, mpm: sigma * Math.sqrt(2 * Math.log(n)), nCycles: n };
}

/** Point fixed to the unit, given relative to the keel origin, expressed about the waterline motion origin. */
export const pointInputs = ({ x = 30, y = 0, z = 55, group = 'Point of interest', what = 'helideck, crane tip or process equipment' } = {}) => [
  { key: 'px', label: 'Point x (forward of midship)', unit: 'm', type: 'number', default: x, min: -300, max: 300, group, role: 'geometry', help: `Longitudinal position of the ${what}.` },
  { key: 'py', label: 'Point y (to port of centreline)', unit: 'm', type: 'number', default: y, min: -100, max: 100, group, role: 'geometry', help: 'Transverse position.' },
  { key: 'pz', label: 'Point height above keel', unit: 'm', type: 'number', default: z, min: 0, max: 300, group, role: 'geometry', help: 'Vertical position measured from the keel.' },
];

/** Rigid-body mass matrix 6×6 about the waterline origin from mass, centre of gravity and gyradii about G. */
export function massMatrix(m, cog, gyr) {
  const M = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  const [xg, yg, zg] = cog; const off = [[0, zg, -yg], [-zg, 0, xg], [yg, -xg, 0]]; const r2 = xg * xg + yg * yg + zg * zg;
  for (let i = 0; i < 3; i++) {
    M[i][i] = m;
    for (let j = 0; j < 3; j++) { M[i][3 + j] = m * off[i][j]; M[3 + j][i] = m * off[i][j]; M[3 + i][3 + j] = m * ((i === j ? r2 + gyr[i] * gyr[i] : 0) - cog[i] * cog[j]); }
  }
  return M;
}
