// Shared helpers of the mooring suite: unit handling, line-component properties, mooring-system inputs and
// builder, safety-factor criteria, advice records and screening sustainability factors. Engineering units are
// converted to SI here so every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { chainProperties, wireRopeProperties, fibreRopeProperties } from '../../engine/lines/lineprops.js';
import { spreadLayout, clusterLayout, turretLayout, lineForces } from '../../engine/lines/mooringsystem.js';
import { solveMultiSegment, compileLine } from '../../engine/lines/multiseg.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const YEAR_S = 31557600;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_POLY = 3.0;       // t CO2e per t polyester / HMPE fibre (screening)
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 30;   // anchor-handling vessel fuel burn, t/day (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };
/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible mooring.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 400) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
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

export const steelCo2Text = (tonnes) => `${fmt(Math.abs(tonnes))} t of steel ≈ ${fmt(Math.abs(tonnes) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e per tonne (confirm with the mill's product declaration).`;
export const vesselCo2Text = (days) => `${fmt(days)} anchor-handling vessel day(s) ≈ ${fmt(days * MGO_T_PER_DAY)} t marine gas oil ≈ ${fmt(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;

export const CONFIRM = 'Default recalled from the standard; the 2021 revision of DNV-OS-E301 changed several factors — confirm against the edition governing the project.';
export const LIMIT_QS = 'The system solver is quasi-static: each line is a static elastic catenary in its own vertical plane, with no current or wave load on the lines, no line inertia or drag and a flat seabed at the anchor level.';
export const LIMIT_PROPS = 'Chain, wire and fibre-rope mass, stiffness and break load come from indicative estimators and must be replaced by certified manufacturer data for design.';

/* ------------------------------------------------------------------------------------- component properties */

export const LINE_TYPES = [
  { value: 'chain-studless', label: 'Studless chain' },
  { value: 'chain-studlink', label: 'Stud-link chain' },
  { value: 'wire-six', label: 'Six-strand wire rope (IWRC)' },
  { value: 'wire-spiral', label: 'Spiral-strand wire rope' },
  { value: 'polyester', label: 'Polyester rope' },
  { value: 'hmpe', label: 'HMPE rope' },
  { value: 'nylon', label: 'Nylon rope' },
  { value: 'custom', label: 'Custom (enter properties)' },
];
const GRADES = ['R3', 'R3S', 'R4', 'R4S', 'R5'];

/**
 * Indicative properties of a line component (SI): mass [kg/m], submerged weight w [N/m], EA [N], MBL [N],
 * volume-equivalent diameter dVol [m]. `over` = { w, EA, MBL, mass } in SI replaces the estimator where > 0.
 */
export function componentProps({ type, d, grade = 'R4', rhoW = 1025, over = {} }) {
  const t = String(type || '').trim().toLowerCase();
  let p, family;
  if (t === 'chain-studless' || t === 'chain-studlink' || t === 'chain') {
    const g = String(grade || 'R4').toUpperCase().replace(/^CHAIN-/, '');
    need(GRADES.includes(g), `Chain grade '${grade}' is not one of ${GRADES.join(', ')}.`);
    p = run(() => chainProperties({ d, grade: g, type: t === 'chain-studlink' ? 'studlink' : 'studless', rhoWater: rhoW }), 'Chain properties');
    family = 'chain';
  } else if (t === 'wire-six' || t === 'wire-spiral') {
    p = run(() => wireRopeProperties({ d, construction: t === 'wire-six' ? 'six-strand-iwrc' : 'spiral-strand', rhoWater: rhoW }), 'Wire-rope properties');
    family = 'wire';
  } else if (t === 'polyester' || t === 'hmpe' || t === 'nylon') {
    p = run(() => fibreRopeProperties({ d, material: t, rhoWater: rhoW }), 'Fibre-rope properties');
    family = 'fibre';
  } else if (t === 'custom') {
    need(over.w !== undefined && over.EA > 0 && over.MBL > 0, 'A custom line component needs submerged weight, axial stiffness and breaking load.');
    const mass = over.mass > 0 ? over.mass : Math.max(Math.abs(over.w) / G / 0.87, 1);
    // volume-equivalent diameter from the displaced volume implied by mass and submerged weight
    const vol = Math.max((mass * G - over.w) / (rhoW * G), 1e-6);
    p = { d: d > 0 ? d : Math.sqrt(4 * vol / Math.PI), mass, w: over.w, EA: over.EA, MBL: over.MBL, dVol: Math.sqrt(4 * vol / Math.PI) };
    family = 'custom';
  } else {
    throw new Error(`Line component type '${type}' is not recognised. Use one of: ${LINE_TYPES.map((o) => o.value).join(', ')}.`);
  }
  const out = { type: t, family, d: p.d, mass: p.mass, w: p.w, EA: p.EA, MBL: p.MBL, dVol: p.dVol };
  if (t !== 'custom') {
    if (Number.isFinite(over.w) && over.w !== 0) out.w = over.w;
    if (over.EA > 0) out.EA = over.EA;
    if (over.MBL > 0) out.MBL = over.MBL;
    if (over.mass > 0) out.mass = over.mass;
  }
  return out;
}

/* ------------------------------------------------------------------------------------- single-line inputs */

/** Inputs describing one uniform line between an anchor and a fairlead (single-line modules). */
export function lineInputs({ type = 'chain-studless', d = 95, L = 1000, X = 900, h = 280, group = 'Line', withGeometry = true } = {}) {
  const custom = { key: 'lineType', in: ['custom'] };
  const out = [
    { key: 'lineType', label: 'Line component', type: 'select', default: type, options: LINE_TYPES, group, role: 'material',
      help: 'Mass, stiffness and break load are estimated from the diameter (indicative, catalogue-style estimators). Choose Custom to enter certified values.' },
    { key: 'chainGrade', label: 'Chain grade', type: 'select', source: 'materials:chain', default: 'chain-r4', group, role: 'material', showIf: { key: 'lineType', in: ['chain-studless', 'chain-studlink'] },
      help: 'Grade sets the break-load constant of IACS UR W22 (R3 … R5).' },
    { key: 'd', label: 'Nominal diameter', unit: 'mm', type: 'number', default: d, min: 20, max: 300, group, role: 'geometry', design: true,
      help: 'Chain bar diameter or rope diameter.' },
    { key: 'wCustom', label: 'Submerged weight (custom)', unit: 'N/m', type: 'number', default: 1500, min: 0.001, max: 20000, group, role: 'material', showIf: custom },
    { key: 'EAcustom', label: 'Axial stiffness EA (custom)', unit: 'MN', type: 'number', default: 800, min: 1, max: 1e8, group, role: 'material', showIf: custom },
    { key: 'MBLcustom', label: 'Minimum breaking load (custom)', unit: 'kN', type: 'number', default: 9000, min: 10, max: 60000, group, role: 'material', showIf: custom },
    { key: 'massCustom', label: 'Mass in air (custom)', unit: 'kg/m', type: 'number', default: 176, min: 0.1, max: 3000, group, role: 'material', showIf: custom },
  ];
  if (withGeometry) out.push(
    { key: 'L', label: 'Unstretched line length', unit: 'm', type: 'number', default: L, min: 10, max: 10000, group: 'Geometry', role: 'geometry', design: true,
      help: 'Paid-out length between anchor shackle and fairlead.' },
    { key: 'X', label: 'Horizontal anchor–fairlead distance', unit: 'm', type: 'number', default: X, min: 1, max: 10000, group: 'Geometry', role: 'bc',
      help: 'Boundary condition: the fairlead position relative to the anchor. Moving the fairlead away raises tension.' },
    { key: 'h', label: 'Fairlead height above seabed', unit: 'm', type: 'number', default: h, min: 5, max: 3500, group: 'Geometry', role: 'bc', bind: 'site.depth',
      help: 'Water depth less fairlead submergence. Bound to the site depth when available.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Geometry', role: 'environment', bind: 'site.rhoW' },
  );
  return out;
}

/** Line property record (SI) from the single-line form. */
export function lineFrom(inp, ctx) {
  const d = inp.d / 1e3;
  let grade = 'R4';
  if (inp.lineType.startsWith('chain')) {
    const m = ctx.material(inp.chainGrade);
    need(m && m.grade, `Chain grade '${inp.chainGrade}' is not in the materials library.`);
    grade = m.grade;
  }
  const over = inp.lineType === 'custom' ? { w: inp.wCustom, EA: inp.EAcustom * 1e6, MBL: inp.MBLcustom * 1e3, mass: inp.massCustom } : {};
  const p = componentProps({ type: inp.lineType, d, grade, rhoW: inp.rhoW ?? 1025, over });
  need(p.w > 0, 'The line is neutrally buoyant or buoyant: a catenary on the seabed needs a positive submerged weight.');
  return { ...p, grade };
}

/** Check that a single line can span anchor to fairlead (not over-constrained). */
export function checkSpan(L, X, h, EA, MBL) {
  const straight = Math.hypot(X, h);
  const maxStretch = L * (1 + MBL / EA);
  need(maxStretch > straight, `Over-constrained: the anchor–fairlead distance (${fmt(straight, 5)} m straight) exceeds the line length even at breaking strain (${fmt(maxStretch, 5)} m). Lengthen the line or move the anchor in.`);
}

/* ------------------------------------------------------------------------------------- system inputs */

const SEG_DEFAULT = [
  { name: 'Bottom chain', type: 'chain-studless', d: 95, grade: 'R4', L: 650, w: 0, EA: 0, MBL: 0 },
  { name: 'Spiral-strand wire', type: 'wire-spiral', d: 100, grade: '', L: 400, w: 0, EA: 0, MBL: 0 },
  { name: 'Top chain', type: 'chain-studless', d: 95, grade: 'R4', L: 100, w: 0, EA: 0, MBL: 0 },
];

/** Mooring-system inputs shared by every system-level module. */
export function systemInputs({ anchorRadius = 1130, nLines = 8, segments = SEG_DEFAULT, depth = 300 } = {}) {
  return [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: depth, min: 20, max: 3500, group: 'Site', role: 'environment', bind: 'site.depth',
      help: 'Still-water depth; the seabed is flat at the anchor level in the quasi-static solver.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'layout', label: 'Line arrangement', type: 'select', default: 'spread', group: 'Mooring configuration', role: 'geometry',
      options: [{ value: 'spread', label: 'Equally spread' }, { value: 'clustered', label: 'Clustered groups' }, { value: 'turret', label: 'Turret (weathervaning)' }],
      help: 'Spread: lines equally spaced in azimuth. Clustered: groups of lines (semi-submersible corners). Turret: lines hang from a turret that does not transmit yaw.' },
    { key: 'nLines', label: 'Number of lines', type: 'int', default: nLines, min: 3, max: 24, group: 'Mooring configuration', role: 'geometry', design: true,
      help: 'For clustered layouts this is the total; it must divide evenly into the groups.' },
    { key: 'nGroups', label: 'Number of groups', type: 'int', default: 4, min: 2, max: 8, group: 'Mooring configuration', role: 'geometry', showIf: { key: 'layout', in: ['clustered'] } },
    { key: 'groupSpreadDeg', label: 'Angle between lines in a group', unit: 'deg', type: 'number', default: 5, min: 0.5, max: 30, group: 'Mooring configuration', role: 'geometry', showIf: { key: 'layout', in: ['clustered'] } },
    { key: 'headingDeg', label: 'Heading of line 1', unit: 'deg', type: 'number', default: 22.5, min: -180, max: 360, group: 'Mooring configuration', role: 'geometry',
      help: 'Azimuth of the first line (or first group axis) from the global x-axis, anticlockwise.' },
    { key: 'anchorRadius', label: 'Anchor radius', unit: 'm', type: 'number', default: anchorRadius, min: 30, max: 12000, group: 'Mooring configuration', role: 'bc', design: true,
      help: 'Horizontal distance from the mooring centre to each anchor. Boundary condition at the anchor end: anchors are fixed points. With the line length it sets the pretension.' },
    { key: 'fairleadRadius', label: 'Fairlead radius', unit: 'm', type: 'number', default: 40, min: 0.5, max: 250, group: 'Mooring configuration', role: 'geometry',
      help: 'Horizontal distance of the fairleads from the vessel (or turret) centre.' },
    { key: 'fairleadDraft', label: 'Fairlead depth below water line', unit: 'm', type: 'number', default: 20, min: -30, max: 200, group: 'Mooring configuration', role: 'geometry' },
    { key: 'turretX', label: 'Turret position forward of midship', unit: 'm', type: 'number', default: 90, min: -200, max: 200, group: 'Mooring configuration', role: 'geometry', showIf: { key: 'layout', in: ['turret'] } },
    { key: 'layoutGeom', label: 'Anchor coordinates (imported)', type: 'geometry', kinds: ['points', 'polyline'], default: null, group: 'Mooring configuration', role: 'geometry',
      help: 'Optional as-surveyed anchor positions (x, y, z with z up, one point per line in line order). When empty the anchors are generated from the radius and headings above.' },
    { key: 'segments', label: 'Line make-up (anchor → fairlead)', type: 'table', group: 'Line properties', role: 'material',
      columns: [
        { key: 'name', label: 'Segment', type: 'text' },
        { key: 'type', label: 'Type', type: 'text' },
        { key: 'd', label: 'Diameter', unit: 'mm', type: 'number' },
        { key: 'grade', label: 'Chain grade', type: 'text' },
        { key: 'L', label: 'Length', unit: 'm', type: 'number' },
        { key: 'w', label: 'Submerged weight (0 = estimate)', unit: 'N/m', type: 'number' },
        { key: 'EA', label: 'EA (0 = estimate)', unit: 'MN', type: 'number' },
        { key: 'MBL', label: 'MBL (0 = estimate)', unit: 'kN', type: 'number' },
      ],
      default: segments,
      help: `Type is one of: ${LINE_TYPES.map((o) => o.value).join(', ')}. Grade R3, R3S, R4, R4S or R5 for chain. Leave weight, EA and MBL at 0 to use the indicative estimators; enter certified values to override. A negative weight is a buoyant section.` },
    { key: 'clumps', label: 'Clump weights and buoys', type: 'table', group: 'Line properties', role: 'geometry',
      columns: [{ key: 'junction', label: 'At junction no.', type: 'number' }, { key: 'weight', label: 'Net submerged weight (− = buoy)', unit: 'kN', type: 'number' }],
      default: [],
      help: 'Point loads at segment junctions: junction 1 is between segments 1 and 2. Negative weight is net buoyancy.' },
    { key: 'CB', label: 'Seabed friction coefficient', unit: '-', type: 'number', default: 0.6, min: 0, max: 2, group: 'Line properties', role: 'bc', uncertain: { dist: 'uniform', min: 0.4, max: 1.0 },
      help: 'Axial friction of the grounded line (chain on clay ≈ 0.6–1.0, wire ≈ 0.3–0.5). It lowers the anchor tension, not the fairlead tension.' },
  ];
}

/** Mean environmental load on the floater (boundary condition at the vessel end). */
export function loadInputs({ F = 3000, dir = 0 } = {}) {
  return [
    { key: 'Fmean', label: 'Mean environmental force', unit: 'kN', type: 'number', default: F, min: 0, max: 60000, group: 'Environmental load', role: 'bc', bind: 'vessel.Fmean',
      uncertain: { dist: 'gumbel', cov: 0.2 }, help: 'Steady wind + current + mean wave-drift force on the hull. The solver applies no current or wave load to the lines themselves.' },
    { key: 'loadDirDeg', label: 'Direction the load acts towards', unit: 'deg', type: 'number', default: dir, min: -180, max: 360, group: 'Environmental load', role: 'bc',
      help: 'From the global x-axis, anticlockwise. A load between two lines is usually the worst intact case; a load in line with one line the worst damaged case.' },
    { key: 'Mz', label: 'Mean yaw moment', unit: 'kN·m', type: 'number', default: 0, min: -5e6, max: 5e6, group: 'Environmental load', role: 'bc' },
  ];
}
export function fextFrom(inp) {
  const F = (inp.Fmean ?? 0) * 1e3, a = (inp.loadDirDeg ?? 0) * DEG;
  return [F * Math.cos(a), F * Math.sin(a), 0, 0, 0, (inp.Mz ?? 0) * 1e3];
}

/** Safety-factor criteria, overridable (see CONFIRM). */
export function criteriaInputs({ damaged = false, offset = true, offsetDefault = 8 } = {}) {
  const g = 'Acceptance criteria';
  const out = [
    { key: 'sfIntact', label: 'API RP 2SK safety factor, intact', unit: '-', type: 'number', default: 1.67, min: 1, max: 4, group: g, role: 'criteria',
      help: `Required MBL / maximum tension, dynamic analysis (2.0 for quasi-static). ${CONFIRM}` },
    { key: 'gammaMean', label: 'DNV partial factor on mean tension', unit: '-', type: 'number', default: damaged ? 1.0 : 1.1, min: 0.8, max: 3, group: g, role: 'criteria',
      help: `DNV-OS-E301 editions up to 2013, consequence class 1: ULS 1.10, ALS 1.00; class 2: ULS 1.40, ALS 1.00. From the 2018 edition the format changed to a pretension factor (1.20) with an environmental-tension factor — if that edition governs, enter its factors here and in the dynamic factor. ${CONFIRM}` },
    { key: 'gammaDyn', label: 'DNV partial factor on dynamic tension', unit: '-', type: 'number', default: damaged ? 1.1 : 1.5, min: 0.8, max: 3.5, group: g, role: 'criteria',
      help: `DNV-OS-E301 editions up to 2013, consequence class 1: ULS 1.50, ALS 1.10; class 2: ULS 2.10, ALS 1.25. The 2018 edition applies 1.45 (class 1) or 1.90 (class 2) to the environmental tension in time-domain analysis. ${CONFIRM}` },
    { key: 'strengthFactor', label: 'Characteristic strength / MBL', unit: '-', type: 'number', default: 0.95, min: 0.7, max: 1, group: g, role: 'criteria',
      help: 'DNV-OS-E301: Sc = 0.95·MBS when strength statistics of the component are not available.' },
  ];
  if (damaged) out.splice(1, 0, { key: 'sfDamaged', label: 'API RP 2SK safety factor, one line broken', unit: '-', type: 'number', default: 1.25, min: 1, max: 3, group: g, role: 'criteria',
    help: `Required MBL / maximum tension in the damaged condition, dynamic analysis (1.43 quasi-static). ${CONFIRM}` });
  if (offset) out.push({ key: 'offsetLimitPct', label: 'Allowable offset', unit: '% depth', type: 'number', default: offsetDefault, min: 0.5, max: 60, group: g, role: 'criteria',
    help: 'Project limit set by risers, umbilicals or gangway (typically 5–8 % of depth intact with rigid risers, 10–12 % damaged, more for flexibles).' });
  return out;
}

/* ------------------------------------------------------------------------------------- system builder */

function importedAnchors(model, n) {
  const pos = model.positions || model.points;
  need(pos && pos.length >= 3, 'The imported anchor layout has no coordinates.');
  const m = Math.floor(pos.length / 3);
  need(m === n, `The imported layout has ${m} point(s) but the system has ${n} line(s): supply one anchor per line, in line order.`);
  const out = [];
  for (let i = 0; i < m; i++) out.push([pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]]);
  return out;
}

/**
 * Build the engine system record from the form. Returns the system plus per-segment property records,
 * the weakest-segment MBL [N], total line mass [kg] split by family and the nominal line length [m].
 */
export function buildSystem(inp) {
  const rhoW = inp.rhoW ?? 1025;
  need(Array.isArray(inp.segments) && inp.segments.length > 0, 'The line make-up table is empty: add at least one segment.');
  const segs = inp.segments.map((r, i) => {
    need(r.L > 0, `Segment ${i + 1}: length must be positive.`);
    const type = String(r.type || '').trim().toLowerCase();
    if (type !== 'custom') need(r.d > 0, `Segment ${i + 1}: diameter must be positive.`);
    const p = componentProps({ type, d: (r.d || 0) / 1e3, grade: r.grade || 'R4', rhoW,
      over: { w: Number(r.w) || (type === 'custom' ? r.w : undefined), EA: (Number(r.EA) || 0) * 1e6, MBL: (Number(r.MBL) || 0) * 1e3 } });
    return { ...p, L: r.L, name: r.name || `segment ${i + 1}`, grade: r.grade };
  });
  const n = segs.length;
  const points = (inp.clumps || []).filter((c) => Number.isFinite(c.junction) && Number.isFinite(c.weight) && c.weight !== 0).map((c) => {
    need(Number.isInteger(c.junction) && c.junction >= 1 && c.junction <= n - 1, `Clump/buoy junction ${c.junction} is outside 1…${n - 1} (junction k lies between segments k and k+1).`);
    return { junction: c.junction, weight: c.weight * 1e3 };
  });
  const Ltot = segs.reduce((a, s) => a + s.L, 0);
  const MBL = Math.min(...segs.map((s) => s.MBL));
  const depth = inp.depth, fz = -inp.fairleadDraft;
  need(inp.fairleadDraft < depth, 'The fairlead is at or below the seabed: reduce the fairlead depth or increase the water depth.');
  need(inp.anchorRadius > inp.fairleadRadius, 'Anchor radius must exceed the fairlead radius.');
  const span = inp.anchorRadius - inp.fairleadRadius, hF = depth - inp.fairleadDraft;
  need(Ltot > hF, `Under-length line: total length ${fmt(Ltot, 5)} m is shorter than the fairlead height above seabed ${fmt(hF, 5)} m.`);
  const stretchMax = segs.reduce((a, s) => a + s.L * (1 + MBL / s.EA), 0);
  need(stretchMax > Math.hypot(span, hF), `Over-constrained: the anchors are ${fmt(Math.hypot(span, hF), 5)} m from the fairleads but the line is only ${fmt(Ltot, 5)} m long (${fmt(stretchMax, 5)} m at breaking strain). Reduce the anchor radius or lengthen a segment.`);
  const base = { anchorRadius: inp.anchorRadius, fairleadRadius: inp.fairleadRadius, fairleadZ: fz, depth, heading: (inp.headingDeg ?? 0) * DEG,
    segments: segs.map((s) => ({ L: s.L, w: s.w, EA: s.EA, name: s.name })), points: points.length ? points : undefined, MBL, CB: inp.CB ?? 0 };
  let system;
  if (inp.layout === 'clustered') {
    need(inp.nLines % inp.nGroups === 0, `${inp.nLines} lines cannot be divided evenly into ${inp.nGroups} groups.`);
    system = run(() => clusterLayout({ ...base, nGroups: inp.nGroups, linesPerGroup: inp.nLines / inp.nGroups, groupSpread: inp.groupSpreadDeg * DEG }), 'Mooring layout');
  } else if (inp.layout === 'turret') {
    system = run(() => turretLayout({ ...base, nLines: inp.nLines, turretCenter: [inp.turretX ?? 0, 0, 0] }), 'Mooring layout');
  } else {
    system = run(() => spreadLayout({ ...base, nLines: inp.nLines }), 'Mooring layout');
  }
  if (inp.layoutGeom) {
    const anchors = importedAnchors(inp.layoutGeom, system.lines.length);
    system.lines.forEach((ln, i) => { ln.anchor = anchors[i]; });
  }
  const massBy = { chain: 0, wire: 0, fibre: 0, custom: 0 };
  for (const s of segs) massBy[s.family] += s.mass * s.L * system.lines.length;
  return { system, segs, MBL, Ltot, depth, hF, nLines: system.lines.length, massBy,
    steel_t: (massBy.chain + massBy.wire + massBy.custom) / 1e3, polymer_t: massBy.fibre / 1e3 };
}

/** Static line forces at a pose with a readable error. */
export const forcesAt = (system, pose, opts) => run(() => lineForces(system, pose, opts), 'Line solve');

/** Bearing of each line (fairlead → anchor), degrees. */
export const lineBearingDeg = (ln) => Math.atan2(ln.anchor[1] - ln.fairlead[1], ln.anchor[0] - ln.fairlead[0]) / DEG;

/** Plan view of the system at a pose: one polyline per line plus the vessel position. */
export function planPlot(system, pose, title = 'Plan view of the mooring') {
  const r = forcesAt(system, pose, { nProfile: 12 });
  const series = r.lines.filter((l) => !l.broken).map((l) => ({ name: l.name, x: l.points.map((p) => round(p[0], 5)), y: l.points.map((p) => round(p[1], 5)) }));
  return { type: 'line', title, x: { label: 'x (east)', unit: 'm' }, y: { label: 'y (north)', unit: 'm' }, series,
    explain: 'Each line from its anchor (outer end) to its fairlead, seen from above at the solved vessel position. Uneven spacing or a line running close to another shows where interference and load concentration will occur.' };
}

/** Elevation of the given lines (distance from anchor vs height above seabed), with the tension along one of them. */
export function elevationPlot(system, pose, indices, title = 'Line profiles (elevation)') {
  const r = forcesAt(system, pose, { nProfile: 40 });
  const series = indices.map((i) => r.lines[i]).filter((l) => l && !l.broken).map((l) => {
    const pr = l.solution.profile;
    return { name: l.name, x: Array.from(pr.x, (v) => round(v, 5)), y: Array.from(pr.z, (v) => round(v, 5)) };
  });
  return { type: 'line', title, x: { label: 'Horizontal distance from anchor', unit: 'm' }, y: { label: 'Height above seabed', unit: 'm' }, series,
    explain: 'Catenary shape of the selected lines. The flat part at zero height is line resting on the seabed; the point where it lifts off is the touchdown. A line with no flat part pulls its anchor upwards.' };
}

/** Table of per-line static results from an equilibrium / lineForces result. */
export function lineTable(lines, MBL, title = 'Line results') {
  return { title,
    columns: [{ key: 'line', label: 'Line' }, { key: 'T', label: 'Fairlead tension', unit: 'kN' }, { key: 'util', label: 'T / MBL', unit: '%' },
      { key: 'aH', label: 'Anchor horizontal load', unit: 'kN' }, { key: 'aV', label: 'Anchor vertical load', unit: 'kN' },
      { key: 'up', label: 'Anchor uplift angle', unit: 'deg' }, { key: 'gr', label: 'Grounded length', unit: 'm' }, { key: 'ang', label: 'Fairlead angle from horizontal', unit: 'deg' }],
    rows: lines.filter((l) => !l.broken).map((l) => ({ line: l.name, T: round(l.T / 1e3), util: round(100 * l.T / MBL, 3), aH: round(l.anchorH / 1e3), aV: round(l.anchorV / 1e3),
      up: round(l.upliftAngle / DEG, 3), gr: round(l.groundedLength), ang: round(Math.atan2(l.V, l.H) / DEG, 3) })),
    explain: 'One row per line: tension at the fairlead and its share of the breaking load, the horizontal and vertical load delivered to the anchor, the uplift angle at the anchor (must be zero for a drag anchor) and the length resting on the seabed.' };
}

/** Standard economic + sustainability advice on line steel for system modules. */
export function steelAdvice(b, ucGov, label = 'line tension') {
  const t = b.steel_t;
  if (ucGov <= 1) {
    // tension capacity ∝ d² (MBL ≈ c·d²): the diameter could shrink by sqrt(uc/0.95) and the mass by uc/0.95
    const save = t * Math.max(0, 1 - ucGov / 0.95);
    return adv('economic', 'info', `The system carries ${fmt(t)} t of chain and wire; the governing ${label} utilisation is ${fmt(ucGov)}.`,
      save > 0.02 * t ? `If nothing else governs (fatigue, corrosion allowance, anchor uplift), a utilisation of 0.95 would be reached with about ${fmt(save)} t less line steel. Re-run with the smaller diameter before committing.` : 'The lines are sized close to their limit; hold the diameter and spend optimisation effort on length and pretension.',
      'Whole-life costing per ISO 15686-5; quantity-based Class 4 estimate per AACE RP 18R-97 (line procurement scales with tonnage).',
      `${steelCo2Text(save > 0 ? save : t)} Mass scales with diameter squared, so a higher grade reaches the same strength with less steel.`);
  }
  const add = t * (ucGov / 0.95 - 1);
  return adv('economic', 'warn', `Bringing the governing ${label} utilisation from ${fmt(ucGov)} to 0.95 by diameter alone needs about ${fmt(add)} t more line steel (${fmt(100 * (ucGov / 0.95 - 1))} %).`,
    'Compare against a higher chain grade (R4S/R5: +11 % / +17 % break load at the same mass), more lines, or a polyester insert before buying heavier chain.',
    'Whole-life costing per ISO 15686-5; quantity-based Class 4 estimate per AACE RP 18R-97.',
    `${steelCo2Text(add)} A grade change adds strength without the extra tonnage.`);
}

/* ------------------------------------------------------------------------------------- tabulated system */


/**
 * Tabulated quasi-static system for horizontal vessel translations at fixed heading and draught.
 * Each line's fairlead tension components H(X), V(X) are solved once with the engine's multi-segment catenary
 * on a grid of spans X and interpolated with cubic Hermite polynomials using the engine's tangent stiffness
 * dH/dX, dV/dX; queries outside the grid go straight to the engine. Same physics as mooringsystem.lineForces,
 * but cheap enough for direction sweeps, time stepping and reliability sampling.
 */
export function fastSystem(system, { maxOffset = 50, yaw = 0, nGrid = 40 } = {}) {
  const c = Math.cos(yaw), s = Math.sin(yaw);
  const tc = system.turretCenter;
  const geo = system.lines.map((ln) => {
    let ax, ay;
    if (ln.turret && tc) { ax = tc[0] * c - tc[1] * s + ln.fairlead[0] - tc[0]; ay = tc[0] * s + tc[1] * c + ln.fairlead[1] - tc[1]; }
    else { ax = ln.fairlead[0] * c - ln.fairlead[1] * s; ay = ln.fairlead[0] * s + ln.fairlead[1] * c; }
    return { ax, ay, Z: ln.fairlead[2] - ln.anchor[2], X0: Math.hypot(ln.anchor[0] - ax, ln.anchor[1] - ay) };
  });
  const tables = new Map();
  const compiled = compileLine(system.lines[0]);
  const direct = (X, Z, guess) => solveMultiSegment({ X, Z, line: compiled, guess });
  const tableFor = (Z) => {
    const key = Z.toFixed(4);
    if (tables.has(key)) return tables.get(key);
    const mine = geo.filter((g) => g.Z.toFixed(4) === key);
    const lo = Math.max(1, Math.min(...mine.map((g) => g.X0)) - maxOffset - 1), hi = Math.max(...mine.map((g) => g.X0)) + maxOffset + 1;
    const t = { lo, hi, dx: (hi - lo) / (nGrid - 1), n: 0, H: [], V: [], dH: [], dV: [], aH: [], aV: [], gr: [] };
    let guess;
    for (let i = 0; i < nGrid; i++) {
      let sol;
      try { sol = direct(lo + i * t.dx, Z, guess); } catch (e) { break; }   // beyond the elastic reach of the line
      if (!sol.slack) guess = { H: sol.H, V: sol.V };
      t.H.push(sol.H); t.V.push(sol.V); t.dH.push(sol.stiffness[0][0]); t.dV.push(sol.stiffness[1][0]);
      t.aH.push(sol.anchorH); t.aV.push(sol.anchorV); t.gr.push(sol.groundedLength); t.n++;
    }
    need(t.n >= 2, 'The line cannot span the anchor–fairlead distance within its elastic range: reduce the anchor radius or lengthen the line.');
    tables.set(key, t);
    return t;
  };
  const lineAt = (g, X) => {
    const t = tableFor(g.Z);
    const u = (X - t.lo) / t.dx;
    if (!(u >= 0 && u <= t.n - 1)) {
      const sol = run(() => direct(X, g.Z), 'Line solve (a line is stretched beyond the tabulated range — it may be near its elastic limit)');
      return { H: sol.H, V: sol.V, anchorH: sol.anchorH, anchorV: sol.anchorV, groundedLength: sol.groundedLength };
    }
    const i = Math.min(Math.floor(u), t.n - 2), p = u - i, p2 = p * p, p3 = p2 * p;
    const h00 = 2 * p3 - 3 * p2 + 1, h10 = p3 - 2 * p2 + p, h01 = -2 * p3 + 3 * p2, h11 = p3 - p2;
    const lin = (a) => a[i] + p * (a[i + 1] - a[i]);
    return { H: Math.max(0, h00 * t.H[i] + h10 * t.dx * t.dH[i] + h01 * t.H[i + 1] + h11 * t.dx * t.dH[i + 1]),
      V: h00 * t.V[i] + h10 * t.dx * t.dV[i] + h01 * t.V[i + 1] + h11 * t.dx * t.dV[i + 1],
      anchorH: lin(t.aH), anchorV: lin(t.aV), groundedLength: lin(t.gr) };
  };
  const forces = (x, y, exclude = []) => {
    const F = [0, 0, 0, 0, 0, 0], lines = [];
    let Tmax = 0, iMax = -1;
    system.lines.forEach((ln, i) => {
      const name = ln.name ?? `line ${i + 1}`;
      if (exclude.includes(i)) { lines.push({ index: i, name, broken: true, H: 0, V: 0, T: 0, anchorH: 0, anchorV: 0, upliftAngle: 0, groundedLength: 0 }); return; }
      const g = geo[i];
      const dx = ln.anchor[0] - (x + g.ax), dy = ln.anchor[1] - (y + g.ay), X = Math.hypot(dx, dy);
      const r = lineAt(g, X);
      const ux = dx / X, uy = dy / X, fx = r.H * ux, fy = r.H * uy;
      F[0] += fx; F[1] += fy; F[2] -= r.V; F[5] += g.ax * fy - g.ay * fx;
      const T = Math.hypot(r.H, r.V);
      lines.push({ index: i, name, broken: false, H: r.H, V: r.V, T, anchorH: r.anchorH, anchorV: r.anchorV, anchorT: Math.hypot(r.anchorH, r.anchorV),
        upliftAngle: r.anchorV > 0 ? Math.atan2(r.anchorV, Math.max(r.anchorH, 1e-9)) : 0, groundedLength: r.groundedLength, X });
      if (T > Tmax) { Tmax = T; iMax = i; }
    });
    return { F, lines, maxTension: Tmax, maxTensionLine: iMax };
  };
  /** Surge–sway equilibrium under Fext = [Fx, Fy, …] from the starting point (x0, y0). */
  const solve = (Fext, start = [0, 0], exclude = []) => {
    let x = start[0], y = start[1];
    let cur = forces(x, y, exclude);
    const scale = cur.lines.reduce((a, l) => a + l.H, 0) + Math.hypot(Fext[0], Fext[1]) + 1;
    const h = 0.01, cap = 0.03 * Math.abs(system.lines[0].anchor[2]) + 3;
    let f = Math.hypot(cur.F[0] + Fext[0], cur.F[1] + Fext[1]);
    for (let it = 0; it < 80 && f > 1e-8 * scale; it++) {
      const ax = forces(x + h, y, exclude), ay = forces(x, y + h, exclude);
      const r1 = cur.F[0] + Fext[0], r2 = cur.F[1] + Fext[1];
      const a = (ax.F[0] - cur.F[0]) / h, b = (ay.F[0] - cur.F[0]) / h, cc = (ax.F[1] - cur.F[1]) / h, d = (ay.F[1] - cur.F[1]) / h;
      const det = a * d - b * cc;
      need(Math.abs(det) > 0, 'The mooring has no horizontal stiffness at this position (all remaining lines slack).');
      let dx = -(d * r1 - b * r2) / det, dy = -(-cc * r1 + a * r2) / det;
      const st = Math.hypot(dx, dy);
      if (st > cap) { dx *= cap / st; dy *= cap / st; }
      let alpha = 1, moved = false;
      for (let ls = 0; ls < 14; ls++) {
        let nx = null;
        try { nx = forces(x + alpha * dx, y + alpha * dy, exclude); } catch (e) { nx = null; }
        const fn = nx ? Math.hypot(nx.F[0] + Fext[0], nx.F[1] + Fext[1]) : Infinity;
        if (fn < f) { x += alpha * dx; y += alpha * dy; cur = nx; f = fn; moved = true; break; }
        alpha *= 0.5;
      }
      if (!moved) break;
    }
    need(f <= 1e-5 * scale, 'Mooring equilibrium did not converge: the applied load probably exceeds what the remaining lines can restore within their elastic range.');
    return { x, y, pose: [x, y, 0, 0, 0, yaw], offset: Math.hypot(x - start[0], y - start[1]), ...cur };
  };
  /** Restoring force (opposing the offset) and line tensions along a direction from a base point. */
  const curve = (direction, offsets, base = [0, 0], exclude = []) => {
    const cd = Math.cos(direction), sd = Math.sin(direction);
    const out = { offset: [], restoring: [], maxTension: [], lineTensions: [], results: [] };
    for (const d of offsets) {
      const r = forces(base[0] + d * cd, base[1] + d * sd, exclude);
      out.offset.push(d); out.restoring.push(-(r.F[0] * cd + r.F[1] * sd)); out.maxTension.push(r.maxTension); out.lineTensions.push(r.lines.map((l) => l.T)); out.results.push(r);
    }
    return out;
  };
  return { forces, solve, curve, geo };
}

/** Linear interpolation on an ascending grid with end clamping. */
export function interp1(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}

/**
 * Build the system and solve the calm-water and mean-load states with the heading held (surge and sway free).
 * `maxOffset` sizes the tabulated range; larger excursions are still solved, directly by the engine.
 */
export function meanState(inp, { maxOffset } = {}) {
  const b = buildSystem(inp);
  const Fext = fextFrom(inp);
  const fs = fastSystem(b.system, { maxOffset: maxOffset ?? Math.max(0.15 * b.depth, 30) });
  const calm = fs.solve([0, 0, 0, 0, 0, 0]);
  const eq = fs.solve(Fext, [calm.x, calm.y]);
  return { b, sys: b.system, fs, Fext, calm, eq };
}
