// Shared helpers of the topsides and deck structures suite: the parametric multi-level deck frame (or an imported
// frame), its load cases, member code checks, plots and the advice/check records. Engineering units are converted
// to SI inside each module (docs/CONTRACT.md §2.1); everything returned by the builders here is SI.
import { analyzeFrame, memberDiagram, sectionStress, refineModel } from '../../engine/structural/frame3d.js';
import { iSection, tubeSection, resolveSection } from '../../engine/structural/sections.js';
import { ec3MemberCheck, ec3Chi } from '../../engine/structural/steelcheck.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { fromImportedFrame } from '../../engine/structural/builders.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
// Screening emission factors (order of magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m³ concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 45;   // heavy-lift / transport spread fuel burn, t/day (screening)

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
export const cap = (v, hi = 99) => (Number.isFinite(v) ? Math.min(v, hi) : hi);
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
export const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 400) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const d = Math.abs(demand);
  const uc = capacity > 0 ? d / capacity : (d > 0 ? 99 : 0);
  return { name, demand: round(d, 5), capacity: round(capacity, 5), unit, uc: round(Math.min(Number.isFinite(uc) ? uc : 99, 99), 5), ref, explain };
};
/** Check written as required factor over achieved factor (buckling factors, safety factors, margins). */
export const fsCheck = (name, fs, fsReq, ref, explain) => (
  { name, demand: round(fsReq, 4), capacity: round(cap(fs, 999), 4), unit: 'factor (required / achieved)', uc: round(fs > 0 ? Math.min(fsReq / fs, 99) : 99, 5), ref, explain });
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (a screening factor — confirm with the mill's product declaration).`;
export const vesselCo2Text = (days) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const CONFIRM = 'Default recalled from the standard — confirm against the governing project edition.';
export const NOT_FE = 'Explicit nonlinear finite elements, CFD-generated loads and weld-level solid models are outside this module.';

/** Steel record in SI from the materials library. */
export function steelOf(ctx, id) {
  const m = ctx.material(id);
  need(m && m.fy, `Steel '${id}' is not in the materials library.`);
  return { E: m.E, nu: m.nu ?? 0.3, rho: m.rho ?? 7850, fy: m.fy, fu: m.fu ?? 1.3 * m.fy, alpha: m.alpha ?? 1.2e-5, G: m.E / (2 * (1 + (m.nu ?? 0.3))), name: m.name };
}

/* ------------------------------------------------------------------------------------------------ inputs */

export const DEFAULT_EQUIPMENT = [
  { name: 'Separator train', W: 180, x: 9, y: 7, level: 1, a: 10, b: 4, zc: 2.5 },
  { name: 'Gas compressor skid', W: 140, x: 21, y: 12, level: 1, a: 8, b: 5, zc: 2.2 },
  { name: 'Power generation', W: 120, x: 8, y: 12, level: 2, a: 9, b: 5, zc: 2.5 },
  { name: 'Living quarters / utilities', W: 220, x: 22, y: 8, level: 2, a: 10, b: 8, zc: 4.0 },
  { name: 'Pipe rack and piping', W: 90, x: 15, y: 3, level: 1, a: 20, b: 3, zc: 3.0 },
];

/**
 * Inputs describing the deck: an imported frame or the parametric multi-level grid, its sections, material and
 * equipment table. `over` overrides defaults by key; `omit` removes keys a module does not use.
 */
export function deckInputs(over = {}, omit = []) {
  const g = 'Deck geometry', m = 'Members', l = 'Permanent and equipment loads';
  const list = [
    { key: 'geom', label: 'Imported frame model (optional)', type: 'geometry', kinds: ['frame'], default: null, group: g, role: 'geometry',
      help: 'SDNF, CIS/2, IFC, structured JSON/CSV or FE beam model. When empty the parametric deck below is analysed. Imported members keep their own sections; members without one take the column tube. The lowest nodes become the jacket interface.' },
    { key: 'Lx', label: 'Deck length (x)', unit: 'm', type: 'number', default: 30, min: 4, max: 150, group: g, role: 'geometry', design: true },
    { key: 'Ly', label: 'Deck width (y)', unit: 'm', type: 'number', default: 20, min: 4, max: 100, group: g, role: 'geometry', design: true },
    { key: 'nx', label: 'Bays along x', type: 'int', default: 3, min: 1, max: 8, group: g, role: 'geometry', help: 'Grid lines of transverse beams; the girder span is Lx / bays.' },
    { key: 'ny', label: 'Bays along y', type: 'int', default: 2, min: 1, max: 6, group: g, role: 'geometry' },
    { key: 'nLevels', label: 'Deck levels', type: 'int', default: 2, min: 1, max: 4, group: g, role: 'geometry' },
    { key: 'zLow', label: 'Lowest deck above the jacket interface', unit: 'm', type: 'number', default: 4, min: 0.5, max: 20, group: g, role: 'geometry', help: 'Length of the deck-leg stubs between the stab-in / interface elevation and the cellar deck.' },
    { key: 'hStorey', label: 'Height between deck levels', unit: 'm', type: 'number', default: 7, min: 2.5, max: 15, group: g, role: 'geometry' },
    { key: 'legs', label: 'Deck legs to the jacket at', type: 'select', default: 'corners', group: g, role: 'bc',
      options: [{ value: 'corners', label: 'Four corners (4-leg jacket)' }, { value: 'perimeter', label: 'Every perimeter grid point' }, { value: 'all', label: 'Every grid point' }] },
    { key: 'columns', label: 'Columns between deck levels at', type: 'select', default: 'perimeter', group: g, role: 'geometry',
      options: [{ value: 'corners', label: 'Four corners' }, { value: 'perimeter', label: 'Every perimeter grid point' }, { value: 'all', label: 'Every grid point' }] },
    { key: 'bracing', label: 'Vertical bracing between levels', type: 'select', default: 'perimeter', group: g, role: 'geometry',
      options: [{ value: 'none', label: 'None (moment frame)' }, { value: 'perimeter', label: 'Diagonals in the perimeter faces (truss rows)' }],
      help: 'Diagonals turn the perimeter frames into storey-deep trusses. They need columns at the perimeter grid points.' },
    { key: 'support', label: 'Restraint at the jacket interface', type: 'select', default: 'fixed', group: g, role: 'bc',
      options: [{ value: 'fixed', label: 'Fixed (welded stab-in, stiff jacket)' }, { value: 'pinned', label: 'Pinned (translations only)' }],
      help: 'Pinned legs need bracing or moment-resisting column heads to be stable.' },
    { key: 'gh', label: 'Girder depth (along x)', unit: 'mm', type: 'number', default: 1500, min: 150, max: 4000, group: m, role: 'geometry', design: true },
    { key: 'gb', label: 'Girder flange width', unit: 'mm', type: 'number', default: 500, min: 80, max: 1500, group: m, role: 'geometry' },
    { key: 'gtw', label: 'Girder web thickness', unit: 'mm', type: 'number', default: 20, min: 4, max: 80, group: m, role: 'geometry' },
    { key: 'gtf', label: 'Girder flange thickness', unit: 'mm', type: 'number', default: 40, min: 6, max: 120, group: m, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'bh', label: 'Beam depth (along y)', unit: 'mm', type: 'number', default: 1200, min: 150, max: 4000, group: m, role: 'geometry', design: true },
    { key: 'bb', label: 'Beam flange width', unit: 'mm', type: 'number', default: 450, min: 80, max: 1500, group: m, role: 'geometry' },
    { key: 'btw', label: 'Beam web thickness', unit: 'mm', type: 'number', default: 18, min: 4, max: 80, group: m, role: 'geometry' },
    { key: 'btf', label: 'Beam flange thickness', unit: 'mm', type: 'number', default: 35, min: 6, max: 120, group: m, role: 'geometry', design: true },
    { key: 'cD', label: 'Leg / column diameter', unit: 'mm', type: 'number', default: 914, min: 100, max: 3000, group: m, role: 'geometry', design: true },
    { key: 'ct', label: 'Leg / column wall thickness', unit: 'mm', type: 'number', default: 32, min: 5, max: 120, group: m, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'brD', label: 'Brace diameter', unit: 'mm', type: 'number', default: 508, min: 60, max: 2000, group: m, role: 'geometry', design: true, showIf: { key: 'bracing', in: ['perimeter'] } },
    { key: 'brt', label: 'Brace wall thickness', unit: 'mm', type: 'number', default: 22, min: 4, max: 80, group: m, role: 'geometry', showIf: { key: 'bracing', in: ['perimeter'] } },
    { key: 'tPlate', label: 'Deck plate thickness', unit: 'mm', type: 'number', default: 10, min: 0, max: 50, group: m, role: 'geometry', design: true, help: 'Counted as weight and mass on every level; 0 for grated decks.' },
    { key: 'secSteel', label: 'Secondary steel allowance', unit: 'kg/m²', type: 'number', default: 70, min: 0, max: 400, group: m, role: 'load', help: 'Stiffeners, stringers, grating, handrails and outfitting steel per unit deck area.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: m, role: 'material' },
    { key: 'equip', label: 'Equipment', type: 'table', group: l, role: 'load',
      help: 'Operating weight, plan position of the centre, deck level (1 = lowest), footprint a × b and height of the centre of gravity above the deck. Each footprint is spread to the surrounding grid nodes.',
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'W', label: 'Weight', unit: 't', type: 'number' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' },
        { key: 'level', label: 'Level', type: 'number' }, { key: 'a', label: 'Footprint a (x)', unit: 'm', type: 'number' }, { key: 'b', label: 'Footprint b (y)', unit: 'm', type: 'number' }, { key: 'zc', label: 'CoG above deck', unit: 'm', type: 'number' }],
      default: DEFAULT_EQUIPMENT },
    { key: 'live', label: 'Variable (live) area load', unit: 'kPa', type: 'number', default: 5, min: 0, max: 50, group: l, role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'Laydown, storage, personnel and maintenance loads on every deck level (ISO 19901-3 Table 1 gives 5–15 kPa by area type).' },
  ];
  return list.filter((i) => !omit.includes(i.key)).map((i) => (over[i.key] !== undefined ? { ...i, default: over[i.key] } : i));
}

export const solverInput = (def = 2) => ({ key: 'nSub', label: 'Elements per member', type: 'int', default: def, min: 1, max: 16, group: 'Solver', role: 'control',
  help: 'Beam elements are exact at their nodes; subdividing adds sampling points (mid-span deflection, mode shapes, buckling shapes).' });

/* ------------------------------------------------------------------------------------------------ model */

const takes = (rule, i, j, nx, ny) => { const ei = i === 0 || i === nx, ej = j === 0 || j === ny; return rule === 'all' || (rule === 'perimeter' && (ei || ej)) || (rule === 'corners' && ei && ej); };

/**
 * Deck frame (SI) from the shared inputs. Girders run along x, beams along y at every level, tubular legs go
 * down to the interface and tubular columns and diagonals connect the levels. Deck plate and secondary steel
 * are carried as non-structural mass on the grid members by tributary area (two-way split).
 */
export function buildDeck(inp, ctx, { wall = 1 } = {}) {
  const steel = steelOf(ctx, inp.mat);
  const tube = (D, t) => { need(t * wall < D / 2, 'A tube wall thickness must be smaller than its radius.'); return tubeSection({ D: D / 1e3, t: (t * wall) / 1e3 }); };
  const isec = (h, b, tw, tf, what) => {
    need(2 * tf < h, `${what}: flange thickness must be less than half the depth.`); need(tw < b, `${what}: web thickness must be less than the flange width.`);
    return run(() => iSection({ h: h / 1e3, b: b / 1e3, tw: (tw * wall) / 1e3, tf: (tf * wall) / 1e3, welded: true }), what);
  };
  const sections = { girder: isec(inp.gh, inp.gb, inp.gtw, inp.gtf, 'Girder'), beam: isec(inp.bh, inp.bb, inp.btw, inp.btf, 'Beam'), column: tube(inp.cD, inp.ct), brace: tube(inp.brD ?? 508, inp.brt ?? 22) };
  const areaMass = (inp.tPlate / 1e3) * steel.rho + inp.secSteel;   // kg/m² of plate + secondary steel
  let model, levelsZ, levelNodes, supNodes, grid = null, trib, planArea, span, x0 = 0, y0 = 0, Lx = inp.Lx, Ly = inp.Ly;
  if (inp.geom) {
    need(inp.geom.kind === 'frame', 'The imported model must be a frame (nodes and members).');
    model = run(() => fromImportedFrame(inp.geom, { defaultSection: { type: 'tube', D: inp.cD / 1e3, t: inp.ct / 1e3 }, supportType: inp.support === 'pinned' ? 'pinned' : 'fixed' }), 'Imported frame');
    model.members.forEach((mm) => { delete mm.material; });
    supNodes = model.supports.map((s) => s.node);
    need(supNodes.length >= 1, 'The imported frame has no nodes at its lowest level to act as the jacket interface.');
    const xs = model.nodes.map((n) => n[0]), ys = model.nodes.map((n) => n[1]), zs = model.nodes.map((n) => n[2]);
    const zmin = Math.min(...zs), zmax = Math.max(...zs);
    need(zmax - zmin > 1e-6, 'The imported frame is flat: it needs members above the support level.');
    x0 = Math.min(...xs); y0 = Math.min(...ys); Lx = Math.max(...xs) - x0; Ly = Math.max(...ys) - y0;
    need(Lx > 1e-6 && Ly > 1e-6, 'The imported frame has no plan extent: a deck needs a plan area in x and y.');
    // Deck levels: elevations that carry horizontal members.
    const tol = 1e-3 * (zmax - zmin), lev = [];
    model.members.forEach((mm) => { const za = zs[mm.a], zb = zs[mm.b]; if (Math.abs(za - zb) < tol && za > zmin + tol && !lev.some((z) => Math.abs(z - za) < 10 * tol)) lev.push(za); });
    need(lev.length, 'The imported frame has no horizontal members above its supports: no deck level was found.');
    lev.sort((a, b) => a - b);
    levelsZ = lev.map((z) => z - zmin);
    model.nodes = model.nodes.map((n) => [n[0], n[1], n[2] - zmin]);
    levelNodes = lev.map((z) => model.nodes.map((n, i) => i).filter((i) => Math.abs(zs[i] - z) < 10 * tol));
    planArea = Lx * Ly;
    // Area mass smeared over the horizontal members of each level.
    trib = new Array(model.members.length).fill(0);
    lev.forEach((z) => {
      const ids = model.members.map((mm, i) => i).filter((i) => Math.abs(zs[model.members[i].a] - z) < 10 * tol && Math.abs(zs[model.members[i].b] - z) < 10 * tol);
      const len = ids.reduce((s, i) => { const a = model.nodes[model.members[i].a], b = model.nodes[model.members[i].b]; return s + Math.hypot(b[0] - a[0], b[1] - a[1]); }, 0);
      ids.forEach((i) => { trib[i] = planArea / len; });
    });
    const lens = model.members.map((mm) => { const a = model.nodes[mm.a], b = model.nodes[mm.b]; return Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]); });
    span = Math.max(...lens.filter((v, i) => trib[i] > 0));
  } else {
    const { nx, ny, nLevels } = inp;
    need(!(inp.bracing === 'perimeter' && nLevels > 1 && inp.columns === 'corners' && (nx > 1 || ny > 1)), 'Perimeter bracing needs columns at the perimeter grid points: set the inter-deck columns to "perimeter" or "all", or switch the bracing off.');
    const dx = Lx / nx, dy = Ly / ny, nodes = [], members = [];
    levelsZ = Array.from({ length: nLevels }, (_, k) => inp.zLow + k * inp.hStorey);
    grid = levelsZ.map((z) => { const gk = []; for (let i = 0; i <= nx; i++) { gk.push([]); for (let j = 0; j <= ny; j++) { nodes.push([i * dx, j * dy, z]); gk[i].push(nodes.length - 1); } } return gk; });
    trib = [];
    // Two-way tributary area of a dx × dy cell: trapezoids to the long edges, triangles to the short edges.
    const s = Math.min(dx, dy), lg = Math.max(dx, dy), aLong = (s / 2) * (lg - s / 2), aShort = (s * s) / 4;
    const ax = dx >= dy ? aLong : aShort, ay = dx >= dy ? aShort : aLong;
    grid.forEach((gk) => {
      for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) { members.push({ a: gk[i][j], b: gk[i + 1][j], section: 'girder', group: 'Girder' }); trib.push(((j === 0 || j === ny ? 1 : 2) * ax) / dx); }
      for (let i = 0; i <= nx; i++) for (let j = 0; j < ny; j++) { members.push({ a: gk[i][j], b: gk[i][j + 1], section: 'beam', group: 'Beam' }); trib.push(((i === 0 || i === nx ? 1 : 2) * ay) / dy); }
    });
    supNodes = [];
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) if (takes(inp.legs, i, j, nx, ny)) { nodes.push([i * dx, j * dy, 0]); supNodes.push(nodes.length - 1); members.push({ a: nodes.length - 1, b: grid[0][i][j], section: 'column', group: 'Deck leg' }); trib.push(0); }
    for (let k = 0; k + 1 < nLevels; k++) {
      for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) if (takes(inp.columns, i, j, nx, ny)) { members.push({ a: grid[k][i][j], b: grid[k + 1][i][j], section: 'column', group: 'Column' }); trib.push(0); }
      if (inp.bracing === 'perimeter') {
        // Diagonals rise toward mid-span so that they work in compression under gravity, as in an N-truss.
        for (const j of [0, ny]) for (let i = 0; i < nx; i++) { const up = i < nx / 2; members.push({ a: grid[k][up ? i : i + 1][j], b: grid[k + 1][up ? i + 1 : i][j], section: 'brace', group: 'Brace' }); trib.push(0); }
        for (const i of [0, nx]) for (let j = 0; j < ny; j++) { const up = j < ny / 2; members.push({ a: grid[k][i][up ? j : j + 1], b: grid[k + 1][i][up ? j + 1 : j], section: 'brace', group: 'Brace' }); trib.push(0); }
      }
    }
    model = { nodes, members, supports: supNodes.map((n) => ({ node: n, fix: inp.support === 'pinned' ? 'pinned' : 'fixed' })) };
    levelNodes = grid.map((gk) => gk.flat());
    planArea = Lx * Ly; span = Math.max(dx, dy);
  }
  model.materials = { steel: { E: steel.E, nu: steel.nu, rho: steel.rho, alpha: steel.alpha, fy: steel.fy } };
  model.sections = sections;
  const orig = model.members.map((mm, i) => {
    mm.material = 'steel';
    if (trib[i] > 0) mm.extraMass = areaMass * trib[i];
    const sec = resolveSection(typeof mm.section === 'string' ? sections[mm.section] : mm.section);
    const a = model.nodes[mm.a], b = model.nodes[mm.b], L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    need(L > 1e-9, 'The frame contains a zero-length member.');
    const kind = sec.D && sec.t && !sec.h ? 'tube' : (sec.h && sec.tw && sec.tf && sec.Wply ? 'I' : 'other');
    const group = mm.group || (Math.abs(b[2] - a[2]) < 1e-6 ? 'Beam' : 'Column');
    return { sec, L, kind, group, mid: [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2], horizontal: Math.abs(b[2] - a[2]) < 1e-6 * L + 1e-9 };
  });
  const frameSteelT = orig.reduce((s, o) => s + o.sec.A * o.L * steel.rho, 0) / 1e3;
  const plateSteelT = (areaMass * planArea * levelsZ.length) / 1e3;
  return { model, steel, sections, orig, levelsZ, levelNodes, supNodes, grid, trib, planArea, span, x0, y0, Lx, Ly, areaMass, frameSteelT, plateSteelT, steelT: frameSteelT + plateSteelT, imported: !!inp.geom,
    dx: grid ? Lx / inp.nx : span, dy: grid ? Ly / inp.ny : span };
}

/** Equipment rows validated and converted to SI: W [N], mass [kg], position, level index (0-based). */
export function equipmentOf(deck, inp) {
  const rows = (inp.equip || []).filter((e) => e && Number(e.W) > 0);
  return rows.map((e, n) => {
    const lev = Math.round(Number(e.level) || 1);
    need(lev >= 1 && lev <= deck.levelsZ.length, `Equipment '${e.name || n + 1}' is on level ${lev}, but the deck has ${deck.levelsZ.length} level(s).`);
    const x = Number(e.x), y = Number(e.y), a = Math.max(Number(e.a) || 0, 0), b = Math.max(Number(e.b) || 0, 0), zc = Math.max(Number(e.zc) || 0, 0);
    need(Number.isFinite(x) && Number.isFinite(y), `Equipment '${e.name || n + 1}' needs a plan position x, y.`);
    need(x >= deck.x0 - 1e-6 && x <= deck.x0 + deck.Lx + 1e-6 && y >= deck.y0 - 1e-6 && y <= deck.y0 + deck.Ly + 1e-6, `Equipment '${e.name || n + 1}' at (${x}, ${y}) m is outside the deck plan (${fmt(deck.Lx)} × ${fmt(deck.Ly)} m).`);
    return { name: String(e.name || `Item ${n + 1}`), mass: Number(e.W) * 1e3, W: Number(e.W) * 1e3 * G, x, y, k: lev - 1, a, b, zc, z: deck.levelsZ[lev - 1] + zc };
  });
}

/** Shares of each equipment item on the frame nodes: [{ node, w (fraction), eq }] — the resultant is preserved. */
export function equipmentShares(deck, eqs) {
  const out = [];
  eqs.forEach((e, q) => {
    const acc = new Map();
    const add = (n, w) => acc.set(n, (acc.get(n) || 0) + w);
    if (deck.grid) {
      const gk = deck.grid[e.k], nx = gk.length - 1, ny = gk[0].length - 1, dx = deck.Lx / nx, dy = deck.Ly / ny;
      // Sample the footprint at 3 × 3 points and distribute each sample bilinearly to its grid cell.
      for (let p = 0; p < 3; p++) for (let r = 0; r < 3; r++) {
        const x = clamp(e.x + ((p - 1) * e.a) / 3, 0, deck.Lx), y = clamp(e.y + ((r - 1) * e.b) / 3, 0, deck.Ly);
        const i = Math.min(Math.floor(x / dx), nx - 1), j = Math.min(Math.floor(y / dy), ny - 1), u = x / dx - i, v = y / dy - j;
        add(gk[i][j], ((1 - u) * (1 - v)) / 9); add(gk[i + 1][j], (u * (1 - v)) / 9); add(gk[i][j + 1], ((1 - u) * v) / 9); add(gk[i + 1][j + 1], (u * v) / 9);
      }
    } else {
      let best = -1, bd = Infinity;
      for (const n of deck.levelNodes[e.k]) { const c = deck.model.nodes[n], d = Math.hypot(c[0] - e.x, c[1] - e.y); if (d < bd) { bd = d; best = n; } }
      add(best, 1);
    }
    for (const [node, w] of acc) if (w > 1e-12) out.push({ node, w, eq: q });
  });
  return out;
}

/**
 * Static load case. `dir` is the acceleration field in units of g acting on every mass ([0, 0, -1] is gravity;
 * [0.2, 0, -1] adds a 0.2 g inertia load along +x). Horizontal equipment inertia acts at the equipment centre of
 * gravity, so it also brings a couple onto the deck nodes.
 */
export function loadCase(deck, inp, { name = 'design', gDead = 1, gEquip = 1, gLive = 1, live = (inp.live ?? 0) * 1e3, dir = [0, 0, -1], extra = [], eqs } = {}) {
  const e = eqs || equipmentOf(deck, inp);
  const nodal = extra.map((x) => ({ node: x.node, F: x.F.slice() }));
  for (const s of equipmentShares(deck, e)) {
    const q = e[s.eq], f = q.W * gEquip * s.w;
    nodal.push({ node: s.node, F: [f * dir[0], f * dir[1], f * dir[2], -q.zc * f * dir[1], q.zc * f * dir[0], 0] });
  }
  const member = [];
  if (live > 0 && gLive > 0) deck.trib.forEach((t, i) => { if (t > 0) member.push({ member: i, type: 'uniform', q: [live * gLive * t * dir[0], live * gLive * t * dir[1], live * gLive * t * dir[2]] }); });
  const lc = { name, nodal, member };
  if (gDead > 0) lc.selfWeight = { g: [G * dir[0], G * dir[1], G * dir[2]], factor: gDead };
  return lc;
}

/** Weight breakdown [N] and centre of gravity [m] of the deck with its equipment (and optionally the live load). */
export function weights(deck, inp, { liveFrac = 0, eqs } = {}) {
  const e = eqs || equipmentOf(deck, inp);
  let W = 0, sx = 0, sy = 0, sz = 0;
  const add = (w, x, y, z) => { W += w; sx += w * x; sy += w * y; sz += w * z; };
  let frame = 0, plate = 0, live = 0;
  deck.orig.forEach((o, i) => {
    const wf = o.sec.A * o.L * deck.steel.rho * G, wp = (deck.model.members[i].extraMass || 0) * o.L * G, wl = liveFrac * (inp.live ?? 0) * 1e3 * deck.trib[i] * o.L;
    frame += wf; plate += wp; live += wl; add(wf + wp + wl, o.mid[0], o.mid[1], o.mid[2]);
  });
  let equip = 0;
  for (const q of e) { equip += q.W; add(q.W, q.x, q.y, q.z); }
  return { frame, plate, equip, live, total: W, cog: [sx / W, sy / W, sz / W] };
}

/** Solve one or more load cases on the refined model. Returns { ref, sol, r(name) }. */
export function solveDeck(deck, cases, nSub = 1, extraModel = {}) {
  const model = { ...deck.model, ...extraModel, loadCases: cases };
  const ref = run(() => refineModel(model, Math.max(1, Math.round(nSub))), 'Mesh refinement');
  const sol = run(() => analyzeFrame(ref.model, { cases: cases.map((c) => c.name) }), 'Frame analysis (check supports, bracing and connectivity)');
  for (const c of cases) {
    const r = sol.results[c.name];
    need(r && r.displacements, 'The frame analysis returned no result: check supports and connectivity.');
    for (let i = 0; i < r.displacements.length; i++) need(Number.isFinite(r.displacements[i]), 'The frame is unstable (a mechanism): fix the legs, add bracing or check that all members are connected.');
  }
  return { ref, sol, model };
}

/** Extreme section forces of every original member over its sub-elements (N signed at the governing end). */
export function memberForces(deck, ref, r) {
  return deck.orig.map((o, i) => {
    let N = 0, My = 0, Mz = 0, Vy = 0, Vz = 0, T = 0, defl = 0;
    for (const k of ref.memberMap[i]) {
      const d = memberDiagram(r, k, 9);
      for (const n of [d.N[0], d.N[d.N.length - 1]]) if (Math.abs(n) > Math.abs(N)) N = n;
      My = Math.max(My, d.maxAbs.My); Mz = Math.max(Mz, d.maxAbs.Mz); Vy = Math.max(Vy, d.maxAbs.Vy); Vz = Math.max(Vz, d.maxAbs.Vz); T = Math.max(T, d.maxAbs.T);
    }
    return { i, N, My, Mz, Vy, Vz, T, defl };
  });
}

/**
 * Elastic stresses of any section [Pa]. Imported generic sections { A, Iy, Iz, J } carry no outline, so the extreme
 * fibre is not known: it is placed at 1.4 radii of gyration (between an I-section, 1.2, and a thin tube, 1.41) and
 * the shear area is taken as half the area. Sections with a known shape use the exact section stress.
 */
export function stressOf(sec, f) {
  const ss = sectionStress({ section: sec, ...f });
  if (Number.isFinite(ss.vonMises)) return ss;
  const cy = 1.4 * Math.sqrt(sec.Iy / sec.A), cz = 1.4 * Math.sqrt(sec.Iz / sec.A);
  const sa = (f.N || 0) / sec.A, sb = (Math.abs(f.My || 0) * cy) / sec.Iy + (Math.abs(f.Mz || 0) * cz) / sec.Iz;
  const tau = Math.hypot(f.Vy || 0, f.Vz || 0) / (0.5 * sec.A) + (Math.abs(f.T || 0) * Math.max(cy, cz)) / Math.max(sec.J, 1e-12);
  const s = Math.max(Math.abs(sa + sb), Math.abs(sa - sb));
  return { sigmaAxial: sa, sigmaBending: sb, sigmaMax: sa + sb, sigmaMin: sa - sb, tau, vonMises: Math.hypot(s, Math.sqrt(3) * tau) };
}

/**
 * Code check of every member: EN 1993-1-1 for I-sections, ISO 19902 §13 for tubulars, von Mises for anything else.
 * fyFac scales the yield strength (elevated temperature, accidental limit state); K is the effective length factor;
 * latFrac is the laterally unrestrained fraction of a deck beam bay (stringers and plate restrain the compression flange).
 */
export function memberChecks(deck, forces, { K = 1, gammaM0 = 1.0, gammaM1 = 1.0, fyFac = 1, EFac = 1, isoGamma, latFrac = 0.25 } = {}) {
  const st = deck.steel, fy = st.fy * fyFac, E = st.E * EFac;
  return forces.map((f) => {
    const o = deck.orig[f.i];
    const ss = stressOf(o.sec, { N: f.N, Vy: f.Vy, Vz: f.Vz, T: f.T, My: f.My, Mz: f.Mz });
    let uc = ss.vonMises / (fy / gammaM0), ref = 'von Mises yield';
    if (o.kind === 'other' && f.N < 0) {
      // Generic section in compression: flexural buckling on curve c about the weaker axis, bending stress added linearly.
      const Imin = Math.min(o.sec.Iy, o.sec.Iz), lam = Math.sqrt((o.sec.A * fy) / ((Math.PI ** 2 * E * Imin) / (K * o.L) ** 2)), ub = -f.N / ((ec3Chi(lam, 'c') * o.sec.A * fy) / gammaM1) + ss.sigmaBending / (fy / gammaM0);
      if (ub > uc) { uc = ub; ref = 'EN 1993-1-1 §6.3.1 flexural buckling (curve c) + bending, generic section'; }
    }
    try {
      if (o.kind === 'I') {
        const c = ec3MemberCheck({ section: o.sec, fy, E, G: E / 2.6, NEd: f.N, MyEd: f.My, MzEd: f.Mz, VzEd: f.Vz, Lcry: K * o.L, Lcrz: (o.horizontal ? latFrac : K) * o.L, LLT: (o.horizontal ? latFrac : 1) * o.L, gammaM0, gammaM1 });
        uc = c.governing.uc; ref = `EN 1993-1-1 ${c.governing.ref || ''} ${c.governing.name}`.replace(/\s+/g, ' ');
      } else if (o.kind === 'tube') {
        const c = tubularMemberCheck({ code: 'ISO', section: o.sec, fy, L: o.L, forces: { N: f.N, My: f.My, Mz: f.Mz, Vy: Math.hypot(f.Vy, f.Vz), T: f.T }, options: { Ky: K, Kz: K, E, ...(isoGamma ? { gamma: isoGamma } : {}) } });
        uc = c.governing.uc; ref = `ISO 19902 §13 ${c.governing.name}`;
      }
    } catch (e) { ref = `von Mises yield (code check not applicable: ${e.message})`; }
    uc = Number.isFinite(uc) ? Math.min(Math.max(uc, 0), 99) : 99;
    return { ...f, group: o.group, L: o.L, kind: o.kind, vm: ss.vonMises, uc, ref };
  });
}

export function groupSummary(rows) {
  const by = {};
  for (const x of rows) { const g = by[x.group] || (by[x.group] = { group: x.group, uc: -1, n: 0 }); g.n++; if (x.uc > g.uc) Object.assign(g, { uc: x.uc, N: x.N, My: x.My, Mz: x.Mz, V: Math.hypot(x.Vy, x.Vz), T: x.T, vm: x.vm, ref: x.ref, member: x.i }); }
  return Object.values(by);
}

/** Nodal values on the refined mesh from per-original-member values (largest of the members meeting at a node). */
export function nodeValuesFromMembers(ref, values) {
  const out = new Array(ref.model.nodes.length).fill(0);
  values.forEach((v, i) => { for (const n of ref.nodeMap[i]) out[n] = Math.max(out[n], v); });
  return out;
}

/** 3-D frame plot with nodal values and an exaggerated deformed shape. u: 6 dofs per node [m, rad]. */
export function framePlot({ title, ref, u, values, valueLabel, valueUnit, explain, scaleTo = 0.06 }) {
  const P = ref.model.nodes, n = P.length;
  const plot = { type: 'model3d', title, nodes: P.map((p) => [round(p[0], 6), round(p[1], 6), round(p[2], 6)]), lines: ref.model.members.map((m) => [m.a, m.b]), nodeValues: values.map((v) => round(v, 4)), valueLabel, valueUnit, explain };
  if (u) {
    let umax = 0, size = 0;
    for (let i = 0; i < n; i++) { umax = Math.max(umax, Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2])); size = Math.max(size, Math.abs(P[i][0]), Math.abs(P[i][1]), Math.abs(P[i][2])); }
    const sc = umax > 0 ? (scaleTo * size) / umax : 0;
    plot.deformed = P.map((p, i) => [round(p[0] + sc * u[6 * i], 6), round(p[1] + sc * u[6 * i + 1], 6), round(p[2] + sc * u[6 * i + 2], 6)]);
    plot.explain = `${explain} The deformed shape is exaggerated ${fmt(sc, 2)} times.`;
  }
  return plot;
}

/** Translation magnitudes [m] per node. */
export function nodeDisp(u, n) { const d = new Array(n); for (let i = 0; i < n; i++) d[i] = Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2]); return d; }

/** Sum of support reactions about the plan centre of the interface: forces [N] and moments [N·m]. */
export function interfaceResultant(deck, r) {
  const xc = deck.x0 + deck.Lx / 2, yc = deck.y0 + deck.Ly / 2;
  const R = [0, 0, 0, 0, 0, 0], rows = [];
  for (const s of r.reactions) {
    const p = deck.model.nodes[s.node]; if (!p) continue;
    const f = s.R; R[0] += f[0]; R[1] += f[1]; R[2] += f[2];
    R[3] += f[3] + (p[1] - yc) * f[2] - p[2] * f[1]; R[4] += f[4] + p[2] * f[0] - (p[0] - xc) * f[2]; R[5] += f[5] + (p[0] - xc) * f[1] - (p[1] - yc) * f[0];
    rows.push({ node: s.node, x: p[0], y: p[1], Rx: f[0], Ry: f[1], Rz: f[2], Mx: f[3], My: f[4] });
  }
  return { R, rows };
}

export const DECK_LIMITS = 'Linear elastic space frame with rigid joints: girders and beams are welded I-sections, legs, columns and braces are tubes; the deck plate is weight and mass only (no diaphragm or composite action). Joint, weld and local plate checks are separate modules.';

/** model3d plot of a 2-D Q4/T3 mesh (in the x–y plane, or lifted by `z`) with nodal values. */
export function meshPlot({ title, mesh, values, valueLabel, valueUnit, explain, z, deformed }) {
  const nn = mesh.nn, nodes = [], tris = [];
  for (let i = 0; i < nn; i++) nodes.push([round(mesh.nodes[2 * i], 6), round(mesh.nodes[2 * i + 1], 6), z ? round(z[i], 6) : 0]);
  for (let e = 0; e < mesh.ne; e++) {
    const c = mesh.conn, o = e * mesh.nen;
    tris.push([c[o], c[o + 1], c[o + 2]]);
    if (mesh.nen >= 4 && mesh.type !== 'T6') tris.push([c[o], c[o + 2], c[o + 3]]);
  }
  const plot = { type: 'model3d', title, nodes, tris, lines: [], nodeValues: Array.from(values, (v) => round(v, 4)), valueLabel, valueUnit, explain };
  if (deformed) plot.deformed = deformed;
  return plot;
}

/**
 * Members for eigenvalue (modal and buckling) analyses. Welded plate and stringers stop a deck girder twisting
 * about its own axis; a bare open I-section (no warping stiffness in a beam element) would otherwise show clusters
 * of fictitious twist modes and a fictitious torsional buckling load. The St-Venant constant of the deck girders
 * and beams is therefore raised to a nominal 5 % of the strong-axis second moment of area.
 */
export function twistRestrained(deck) {
  return deck.model.members.map((m, i) => { const o = deck.orig[i]; return o.horizontal && o.kind === 'I' && o.sec.J < 0.05 * o.sec.Iy ? { ...m, section: { ...o.sec, J: 0.05 * o.sec.Iy } } : { ...m }; });
}

/**
 * Model for dynamic analysis: equipment as lumped node masses, plate and secondary steel as member mass, plus a
 * fraction of the live load as mass on the grid members, with twist of the deck girders restrained by the plating.
 * Returns a refined copy ({ model, memberMap, nodeMap }).
 */
export function dynamicModel(deck, inp, { liveFrac = 0, nSub = 1, eqs } = {}) {
  const e = eqs || equipmentOf(deck, inp);
  const liveMass = (liveFrac * (inp.live ?? 0) * 1e3) / G;
  const members = twistRestrained(deck).map((m, i) => (deck.trib[i] > 0 ? { ...m, extraMass: (m.extraMass || 0) + liveMass * deck.trib[i] } : m));
  const acc = new Map();
  for (const s of equipmentShares(deck, e)) acc.set(s.node, (acc.get(s.node) || 0) + e[s.eq].mass * s.w);
  const nodeMasses = [...acc].map(([node, m]) => ({ node, m }));
  const model = { ...deck.model, members, nodeMasses, loadCases: [{ name: 'none', nodal: [] }] };
  return run(() => refineModel(model, Math.max(1, Math.round(nSub))), 'Mesh refinement');
}
