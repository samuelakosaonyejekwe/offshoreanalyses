// Tabular geometry: delimited text (CSV, TSV, semicolon, pipe or whitespace separated; XYZ lists) and the
// table → model logic shared with the JSON and XML readers.
//
// These files have no standard schema, so columns are mapped to roles by header name (see ROLE_NAMES) and the
// mapping that was used is always returned in `meta.mapping`. Pass `opts.mapping = { role: columnNameOrIndex }`
// to override any of it, and `opts.schema` to force the interpretation:
//
//   'points'  x y z soundings / point cloud          → points, or grid when the points form a regular lattice
//   'route'   KP / easting / northing / depth rows   → polyline with attributes (one per `line` value)
//   'frame'   node table + member table              → frame
//   'hull'    station / waterline / half-breadth     → lofted surface mirrored about the centre plane
//   'blade'   radius / chord / twist stations        → polyline along the span with attributes
//   'layout'  named positions (turbines, piles …)    → points with names and attributes
//   'section' circular components or layers          → section
//   'panels'  x1 y1 z1 … x4 y4 z4 per row            → surface

import {
  lineCursor, baseName, Grow, Welder, finish, surfaceModel, pointsModel, unitAdvice, cleanName, extOf, tally,
} from './common.js';
import { resolveUnits, suggestUnits, normaliseUnit, UNIT_SCALE, verticalConvention, geoProjector, looksGeographic, epsgInfo } from '../units.js';
import { resolveSection } from './sections.js';

// ---------------------------------------------------------------------------------------------------------
// Reading delimited text into column tables
// ---------------------------------------------------------------------------------------------------------

const keyOf = (s) => String(s).toLowerCase().replace(/[^a-z0-9]/g, '');

/** Split "Depth (m)" / "x[mm]" / "KP_km" into a name and a unit. */
export function splitHeader(h) {
  const raw = String(h).trim().replace(/^["']|["']$/g, '');
  let m = /^(.*?)[\s_]*[([]\s*([^)\]]*)\s*[)\]]\s*$/.exec(raw);
  if (m && m[1]) return { name: m[1].trim(), unit: m[2].trim() || null };
  m = /^(.+?)[_ ](mm|cm|m|km|ft|in|deg|rad)$/i.exec(raw);
  if (m && m[1].length > 0) return { name: m[1].trim(), unit: m[2] };
  return { name: raw, unit: null };
}

function splitLine(line, delim) {
  if (delim === ' ') return line.trim().split(/\s+/);
  if (line.indexOf('"') < 0) { const p = line.split(delim); for (let i = 0; i < p.length; i++) p[i] = p[i].trim(); return p; }
  const out = [];
  let cur = '', q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; } else if (c === '"') q = true; else if (c === delim) { out.push(cur.trim()); cur = ''; } else cur += c;
  }
  out.push(cur.trim());
  return out;
}

function sniffDelimiter(line) {
  let best = ' ', n = 0;
  for (const d of [',', ';', '\t', '|']) {
    let c = 0, q = false;
    for (let i = 0; i < line.length; i++) { const ch = line[i]; if (ch === '"') q = !q; else if (!q && ch === d) c++; }
    if (c > n) { n = c; best = d; }
  }
  return best;
}

class Column {
  constructor(name, unit) {
    this.name = name; this.key = keyOf(name); this.unit = unit || null;
    this.g = new Grow(Float64Array, 256); this.str = null; this.num = null;
  }
  push(tok) {
    if (tok === undefined || tok === '') this.pushParsed('', NaN); else this.pushParsed(tok, Number(tok));
  }
  /** Append a cell whose numeric value `v` (NaN when not a number) has already been parsed from `tok`. */
  pushParsed(tok, v) {
    if (v !== v && tok !== '') {
      if (tok.length < 5 && /^(nan|na|n\/a|null|-|--)$/i.test(tok)) tok = '';
      else if (!this.str) { this.str = new Array(this.g.n); for (let i = 0; i < this.g.n; i++) { const x = this.g.a[i]; this.str[i] = x === x ? String(x) : ''; } }
    }
    if (this.str) this.str.push(tok);
    this.g.push(v);
  }
  done() { this.num = this.g.finish(); this.g = null; return this; }
  get n() { return this.num.length; }
  /** Text of row i whatever the column type. */
  text(i) { return this.str ? this.str[i] : (this.num[i] === this.num[i] ? String(this.num[i]) : ''); }
}

/** Make a table from column arrays (used by the JSON / XML readers). */
export function makeTable(name, headers, columns) {
  const cols = headers.map((h, k) => {
    const { name: nm, unit } = splitHeader(h);
    const c = new Column(nm, unit);
    for (const v of columns[k]) c.push(v === null || v === undefined ? '' : typeof v === 'number' ? (v === v ? String(v) : '') : String(v).trim());
    return c.done();
  });
  return { name, headerless: false, cols, n: cols.length ? cols[0].n : 0 };
}

const COMMENT = /^(#|\/\/|!|%)/;
const UNIT_TOKEN = /^[([]?\s*(m|mm|cm|km|ft|in|deg|rad|degrees?|met(?:er|re)s?|millimet(?:er|re)s?|feet|inch(?:es)?|-|)\s*[)\]]?$/i;

/**
 * Read delimited text into one or more tables. Blocks introduced by a `[name]` or `*NAME` line, or by a new
 * header row after data, become separate tables (node table + member table in one file).
 * @returns {{ tables: object[], delimiter: string, comments: string[] }}
 */
export function readDelimited(input) {
  const cur = lineCursor(input);
  const tables = [], comments = [];
  let delim = null, t = null, pendingName = null, decimalComma = false;
  const isNum = (s) => s !== '' && Number(s) === Number(s);
  const start = (headers, headerless) => {
    t = { name: pendingName || (tables.length ? `table ${tables.length + 1}` : baseName(input.name)), headerless, cols: headers.map((h) => { const s = splitHeader(h); return new Column(s.name, s.unit); }), n: 0, ragged: 0, wide: null };
    pendingName = null;
    tables.push(t);
  };
  // Bulk path for the rows of a large file: one pass over the characters, each cell parsed once, no token
  // array. Returns false (nothing appended) when the row needs the general path below.
  const tokS = [], tokV = [];
  const fastRow = (s) => {
    const c0 = s.charCodeAt(0);
    if (!((c0 >= 48 && c0 <= 57) || c0 === 45 || c0 === 46 || c0 === 43) || decimalComma || s.indexOf('"') >= 0) return false;
    const cols = t.cols, nc = cols.length, L = s.length;
    const ws = delim === ' ', dc = ws ? 32 : delim.charCodeAt(0);
    let i = 0, k = 0;
    while (i <= L) {
      if (ws) { while (i < L) { const c = s.charCodeAt(i); if (c === 32 || c === 9) i++; else break; } if (i >= L) break; }
      let e = i;
      if (ws) { while (e < L) { const c = s.charCodeAt(e); if (c === 32 || c === 9) break; e++; } } else { while (e < L && s.charCodeAt(e) !== dc) e++; }
      if (k >= nc) return false;
      let a = i, b = e;
      if (!ws) { while (a < b && s.charCodeAt(a) <= 32) a++; while (b > a && s.charCodeAt(b - 1) <= 32) b--; }
      const tk = a === i && b === e ? s.slice(i, e) : s.slice(a, b);
      tokS[k] = tk;
      tokV[k] = tk === '' ? NaN : +tk;
      k++;
      i = e + 1;
      if (!ws && e >= L) break;
    }
    // a trailing delimiter leaves one empty cell too many
    if (!ws && k === nc + 1 && tokS[k - 1] === '') k--;
    if (k !== nc || tokV[0] !== tokV[0]) return false;
    for (let q = 0; q < nc; q++) cols[q].pushParsed(tokS[q], tokV[q]);
    t.n++;
    return true;
  };
  for (let line = cur.next(); line !== null; line = cur.next()) {
    const s = line.trim();
    if (!s) continue;
    if (t !== null && t.n > 0 && delim !== null && fastRow(s)) continue;
    if (COMMENT.test(s)) { if (comments.length < 50) comments.push(s.replace(COMMENT, '').trim()); continue; }
    const sec = /^\[([^\]]+)\]$/.exec(s) || /^\*\s*([A-Za-z][\w \-]*)$/.exec(s);
    if (sec) { pendingName = sec[1].trim(); t = null; continue; }
    if (delim === null) delim = sniffDelimiter(s);
    let tok = splitLine(delim === ' ' ? s : line, delim);
    if (delim !== ' ' && tok.length > 1 && tok[tok.length - 1] === '') tok.pop();
    // whitespace-separated headers: "r (m)  chord (mm)" keeps each bracketed unit with its name
    if (delim === ' ' && /[([]/.test(s) && !isNum(tok[0])) { const merged = []; for (const v of tok) { if (merged.length && /^[([][^)\]]*[)\]]$/.test(v)) merged[merged.length - 1] += ` ${v}`; else merged.push(v); } tok = merged; }
    if (delim === ';' && !decimalComma && tok.some((v) => /^-?\d+,\d+$/.test(v))) decimalComma = true;
    if (decimalComma) tok = tok.map((v) => (/^-?\d*,\d+([eE][-+]?\d+)?$/.test(v) ? v.replace(',', '.') : v));
    const nNum = tok.reduce((c, v) => c + (isNum(v) ? 1 : 0), 0);
    const textRow = nNum === 0 && tok.some((v) => v !== '');
    // "station, 0, 1, 2, 4": a text corner cell followed by numbers is the header of an offset matrix
    const wideHeader = !textRow && tok.length >= 3 && !isNum(tok[0]) && tok[0] !== '' && nNum === tok.length - 1 && (!t || t.n === 0 || pendingName) && /^(station|stn|sta|frame|x|waterline|wl|z|height)/i.test(keyOf(tok[0]));
    if (!t || (textRow && t.n > 0 && tok.length >= 2) || (wideHeader && (!t || t.n > 0))) {
      if (textRow || wideHeader) { start(tok, false); if (wideHeader) t.wide = { corner: tok[0], values: tok.slice(1).map(Number) }; continue; }
      start(tok.map((_, i) => `c${i + 1}`), true);
    } else if (textRow && t.n === 0 && !t.headerless) {
      // a second text row straight after the header: a units row
      if (tok.every((v) => UNIT_TOKEN.test(v))) { tok.forEach((v, i) => { const u = v.replace(/[()[\]\s]/g, ''); if (t.cols[i] && u && u !== '-') t.cols[i].unit = u; }); continue; }
    }
    if (tok.length !== t.cols.length) {
      t.ragged++;
      while (t.cols.length < tok.length) { const c = new Column(`c${t.cols.length + 1}`, null); for (let i = 0; i < t.n; i++) c.push(''); t.cols.push(c); }
    }
    for (let i = 0; i < t.cols.length; i++) t.cols[i].push(tok[i]);
    t.n++;
  }
  for (const tb of tables) tb.cols.forEach((c) => c.done());
  return { tables: tables.filter((tb) => tb.n > 0), delimiter: delim === ' ' ? 'whitespace' : delim === '\t' ? 'tab' : delim, comments };
}

// ---------------------------------------------------------------------------------------------------------
// Column roles
// ---------------------------------------------------------------------------------------------------------

/** Header names (lower case, punctuation removed) recognised for each role. */
export const ROLE_NAMES = {
  x: ['x', 'xm', 'easting', 'eastings', 'east', 'xcoord', 'xcoordinate', 'utmx', 'utme', 'eastingm', 'xpos', 'px'],
  y: ['y', 'ym', 'northing', 'northings', 'north', 'ycoord', 'ycoordinate', 'utmy', 'utmn', 'northingm', 'ypos', 'py'],
  z: ['z', 'zm', 'elevation', 'elev', 'height', 'level', 'alt', 'altitude', 'zcoord', 'zcoordinate', 'rl', 'ele', 'zpos', 'pz'],
  depth: ['depth', 'waterdepth', 'wd', 'sounding', 'soundings', 'dep', 'seabeddepth', 'depthm', 'tvd'],
  lon: ['lon', 'long', 'longitude', 'lng', 'londeg', 'longitudedeg'],
  lat: ['lat', 'latitude', 'latdeg', 'latitudedeg'],
  kp: ['kp', 'chainage', 'ch', 'station', 'sta', 'stationing', 'distance', 'dist', 'md', 'arclength', 'along', 'range'],
  line: ['line', 'lineid', 'linename', 'route', 'routeid', 'branch', 'cable', 'string', 'riser', 'pipeline', 'leg'],
  id: ['id', 'node', 'nodeid', 'nid', 'no', 'nr', 'number', 'name', 'label', 'tag', 'point', 'pointid', 'joint', 'jointid', 'member', 'memberid', 'element', 'elementid', 'eid', 'mid'],
  a: ['nodea', 'n1', 'startnode', 'from', 'nodei', 'end1', 'na', 'node1', 'inode', 'start', 'ja', 'jointa', 'refstart', 'a', 'i'],
  b: ['nodeb', 'n2', 'endnode', 'to', 'nodej', 'end2', 'nb', 'node2', 'jnode', 'end', 'jb', 'jointb', 'refend', 'b', 'j'],
  x1: ['x1', 'xa', 'xstart', 'xi', 'startx', 'ax'], y1: ['y1', 'ya', 'ystart', 'yi', 'starty', 'ay'], z1: ['z1', 'za', 'zstart', 'zi', 'startz', 'az'],
  x2: ['x2', 'xb', 'xend', 'xj', 'endx', 'bx'], y2: ['y2', 'yb', 'yend', 'yj', 'endy', 'by'], z2: ['z2', 'zb', 'zend', 'zj', 'endz', 'bz'],
  D: ['d', 'od', 'diameter', 'outerdiameter', 'dia', 'do', 'outsidediameter', 'dout', 'diam'],
  t: ['t', 'wt', 'thk', 'wall', 'wallthickness', 'tw', 'thick'],
  ID: ['innerdiameter', 'insidediameter', 'di', 'din', 'bore'],
  section: ['section', 'profile', 'sec', 'size', 'sectionname', 'shape', 'designation'],
  material: ['material', 'mat', 'grade', 'steel'],
  group: ['group', 'type', 'layer', 'set', 'class', 'category', 'membertype', 'kind'],
  A: ['area', 'ax', 'csa'], Iy: ['iy', 'iyy'], Iz: ['iz', 'izz'], J: ['ixx', 'it', 'torsion', 'jt', 'j'],
  station: ['station', 'stn', 'sta', 'frame', 'section', 'x', 'xstation', 'stationx'],
  halfBreadth: ['halfbreadth', 'hb', 'offset', 'y', 'halfbeam', 'breadth', 'halfwidth'],
  waterline: ['waterline', 'wl', 'z', 'height', 'heightabovebase', 'hab', 'draft', 'draught'],
  r: ['r', 'radius', 'span', 'rm', 'spanwise', 'spanposition', 'blspn', 'rnodes', 'z', 'station'],
  chord: ['chord', 'c', 'chordm', 'blchord'],
  twist: ['twist', 'twistdeg', 'beta', 'pitch', 'bltwist', 'aerotwst', 'twistangle'],
  relThickness: ['thickness', 'tc', 'relthickness', 'relativethickness', 'toverc', 'thicknessratio'],
  airfoil: ['airfoil', 'aerofoil', 'foil', 'profile', 'afname', 'polar'],
  thickness: ['thickness', 'thk', 't'],
};

const HINT_ROUTE = /route|line|centre|center|riser|cable|pipe|mooring|umbilical|profile|path|track|align|catenary|flowline|tow|trajectory/i;

function findRole(t, role, used, mapping) {
  if (mapping && mapping[role] !== undefined && mapping[role] !== null) {
    const m = mapping[role];
    const i = typeof m === 'number' ? m : t.cols.findIndex((c) => c.name === m || c.key === keyOf(m));
    if (i >= 0 && i < t.cols.length) return i;
    throw new RangeError(`mapping: column "${m}" given for role "${role}" does not exist in table "${t.name}" (columns: ${t.cols.map((c) => c.name).join(', ')}).`);
  }
  if (mapping && mapping[role] === null) return -1;
  for (const nm of ROLE_NAMES[role]) {
    const i = t.cols.findIndex((c, k) => c.key === nm && !used.has(k));
    if (i >= 0) return i;
  }
  return -1;
}

function pick(t, roles, mapping) {
  const used = new Set(), out = {};
  for (const r of roles) { const i = findRole(t, r, used, mapping); out[r] = i; if (i >= 0) used.add(i); }
  return out;
}

const isNumeric = (c) => c && !c.str;
const isIntCol = (c) => { if (!isNumeric(c)) return false; for (let i = 0; i < c.n; i++) if (!Number.isInteger(c.num[i])) return false; return true; };

/**
 * Decide what a table is. Returns `{ schema, basis, roles }` where `basis` says whether the caller forced it
 * ('option'), the column names decided ('headers'), the table/file name did ('name') or it is the default
 * for bare numbers ('default').
 */
export function classify(t, opts, hasMembersElsewhere, hasNodesElsewhere) {
  const mapping = opts.mapping || null;
  const forced = opts.schema || null;
  const keys = new Set(t.cols.map((c) => c.key));
  const has = (...names) => names.some((n) => keys.has(n));
  const hint = `${t.name} ${opts.nameHint || ''}`;
  if (forced && !['points', 'grid', 'route', 'frame', 'hull', 'blade', 'layout', 'section', 'panels'].includes(forced)) throw new RangeError(`schema: "${forced}" is not one of points, grid, route, frame, hull, blade, layout, section, panels.`);
  if (t.headerless && !mapping && !forced) {
    const n = t.cols.length;
    if (n === 12 && t.cols.every(isNumeric)) return { schema: 'panels', basis: 'default' };
    if (n === 1) return { schema: null, basis: 'default' };
    return { schema: 'points', basis: 'default', headerless: true };
  }
  if (forced === 'panels' || (!forced && has('x1') && has('z4') && has('y3'))) return { schema: 'panels', basis: forced ? 'option' : 'headers' };
  if (t.wide || forced === 'hull' || (!forced && has('halfbreadth', 'hb', 'halfbeam', 'halfwidth') ) || (!forced && has('offset') && has('station', 'stn', 'frame') && has('waterline', 'wl', 'z'))) return { schema: 'hull', basis: forced ? 'option' : 'headers' };
  if (forced === 'blade' || (!forced && has('chord', 'blchord') && has('r', 'radius', 'span', 'rm', 'blspn', 'rnodes', 'spanwise', 'spanposition', 'z', 'station'))) return { schema: 'blade', basis: forced ? 'option' : 'headers' };
  const inline = has('x1', 'xa', 'xstart', 'xi', 'startx') && has('x2', 'xb', 'xend', 'xj', 'endx');
  const pairs = [['nodea', 'nodeb'], ['n1', 'n2'], ['startnode', 'endnode'], ['from', 'to'], ['nodei', 'nodej'], ['end1', 'end2'], ['na', 'nb'], ['node1', 'node2'], ['inode', 'jnode'], ['start', 'end'], ['ja', 'jb'], ['jointa', 'jointb'], ['refstart', 'refend'], ['a', 'b'], ['i', 'j']];
  const pair = pairs.find(([a, b]) => keys.has(a) && keys.has(b));
  const mapPair = mapping && mapping.a !== undefined && mapping.b !== undefined;
  const xy = (has(...ROLE_NAMES.x) && has(...ROLE_NAMES.y)) || (has(...ROLE_NAMES.lon) && has(...ROLE_NAMES.lat)) || (mapping && mapping.x !== undefined);
  if (forced === 'frame' || !forced) {
    if (inline) return { schema: 'frame', part: 'inline', basis: forced ? 'option' : 'headers' };
    if ((pair || mapPair) && !xy) return { schema: 'frame', part: 'members', basis: forced ? 'option' : 'headers', pair };
    if (forced === 'frame' || hasMembersElsewhere) { if (xy) return { schema: 'frame', part: 'nodes', basis: forced ? 'option' : 'headers' }; }
  }
  if (forced === 'section' || (!forced && !xy && has(...ROLE_NAMES.D, ...ROLE_NAMES.ID) && !has(...ROLE_NAMES.kp) && !hasNodesElsewhere) || (!forced && /cross.?section|layers?|components?/i.test(t.name) && has(...ROLE_NAMES.D))) return { schema: 'section', basis: forced ? 'option' : 'headers' };
  if (forced === 'route' || (!forced && (has(...ROLE_NAMES.kp.filter((k) => k !== 'station' && k !== 'sta')) || (has('station', 'sta') && xy)))) return { schema: 'route', basis: forced ? 'option' : 'headers' };
  if (!forced && !xy && has(...ROLE_NAMES.kp) && has(...ROLE_NAMES.depth, ...ROLE_NAMES.z)) return { schema: 'route', basis: 'headers' };
  if (forced === 'layout') return { schema: 'layout', basis: 'option' };
  if (forced === 'points' || forced === 'grid') return { schema: 'points', basis: 'option', grid: forced === 'grid' ? true : undefined };
  if (!xy) {
    if (has('x') && has('z')) return { schema: 'route', basis: 'headers' };
    return { schema: null, basis: 'headers' };
  }
  if (has(...ROLE_NAMES.line)) return { schema: 'route', basis: 'headers' };
  if (prefers(opts, 'polyline')) return { schema: 'route', basis: 'option' };
  if (HINT_ROUTE.test(hint)) return { schema: 'route', basis: 'name' };
  const idc = t.cols.find((c) => ROLE_NAMES.id.includes(c.key) || ROLE_NAMES.group.includes(c.key));
  if (idc && idc.str && t.n <= 20000) return { schema: 'layout', basis: 'headers' };
  return { schema: 'points', basis: 'headers' };
}

const prefers = (opts, kind) => Array.isArray(opts.kinds) && opts.kinds.includes(kind) && !opts.kinds.includes('points') && !opts.kinds.includes('grid');

// ---------------------------------------------------------------------------------------------------------
// Builders
// ---------------------------------------------------------------------------------------------------------

function colScale(c, units) {
  const u = c && c.unit ? normaliseUnit(c.unit) : null;
  return u ? UNIT_SCALE[u] : units.scale;
}

function describe(t, roles, units, extra) {
  const m = {};
  for (const r in roles) {
    const i = roles[r];
    if (i < 0 || i === undefined) continue;
    const c = t.cols[i];
    m[r] = { column: c.name, index: i, unit: c.unit || null };
  }
  return Object.assign(m, extra || {});
}

/** Columns not consumed by a role become attributes (numbers as Float64Array, text as string arrays). */
function restAttrs(t, roles, rows) {
  const used = new Set(Object.values(roles).filter((i) => i >= 0));
  const attrs = {}, units = {};
  t.cols.forEach((c, i) => {
    if (used.has(i)) return;
    let key = c.name.replace(/[^\w]+/g, '_').replace(/^_+|_+$/g, '') || `c${i + 1}`;
    while (key in attrs) key += '_';
    if (c.str) attrs[key] = rows ? rows.map((r) => c.str[r]) : c.str.slice();
    else if (rows) { const a = new Float64Array(rows.length); rows.forEach((r, k) => { a[k] = c.num[r]; }); attrs[key] = a; } else attrs[key] = c.num.slice();
    if (c.unit) units[key] = c.unit;
  });
  return { attrs, units };
}

function horizontal(t, opts, units, warnings, ctx) {
  const mapping = opts.mapping || null;
  const r = pick(t, ['lon', 'lat', 'x', 'y', 'z', 'depth'], mapping);
  let geographic = r.lon >= 0 && r.lat >= 0 && !(mapping && mapping.x !== undefined);
  const xi = geographic ? r.lon : r.x, yi = geographic ? r.lat : r.y;
  if (geographic) { r.x = -1; r.y = -1; } else { r.lon = -1; r.lat = -1; }
  if (xi < 0 || yi < 0) return null;
  const xc = t.cols[xi], yc = t.cols[yi];
  if (xc.str || yc.str) throw new Error(`Table "${t.name}": the coordinate columns "${xc.name}" and "${yc.name}" must be numeric.`);
  const crsOpt = opts.crs || ctx.crsHint || null;
  if (!geographic && crsOpt) { const e = /(\d{4,5})/.exec(String(crsOpt)); const info = e ? epsgInfo(Number(e[1])) : null; if (/^geographic$/i.test(crsOpt) || (info && info.kind === 'geographic')) geographic = true; }
  return { r, xi, yi, geographic, crsOpt };
}

function coordsOf(t, opts, units, warnings, ctx, rows) {
  const h = horizontal(t, opts, units, warnings, ctx);
  if (!h) return null;
  const { r, xi, yi } = h;
  const n = rows ? rows.length : t.n;
  const at = (k) => (rows ? rows[k] : k);
  const xs = t.cols[xi].num, ys = t.cols[yi].num;
  const zi = r.z >= 0 ? r.z : r.depth;
  const zc = zi >= 0 ? t.cols[zi] : null;
  if (zc && zc.str) throw new Error(`Table "${t.name}": the vertical column "${zc.name}" must be numeric.`);
  const vert = zc ? verticalConvention(opts, t.headerless ? '' : r.z >= 0 ? 'z' : zc.key) : null;
  const out = new Float64Array(n * 3);
  let crs = null;
  if (h.geographic) {
    let lo = Infinity, hi = -Infinity, la = Infinity, lb = -Infinity;
    for (let k = 0; k < n; k++) { const x = xs[at(k)], y = ys[at(k)]; if (x === x && y === y) { if (x < lo) lo = x; if (x > hi) hi = x; if (y < la) la = y; if (y > lb) lb = y; } }
    if (!looksGeographic(lo, hi, la, lb)) throw new Error(`Table "${t.name}": columns "${t.cols[xi].name}" / "${t.cols[yi].name}" are taken as longitude / latitude but hold values outside ±180 / ±90. Map them as x / y instead (opts.mapping) if they are projected coordinates.`);
    const proj = ctx.projector || (ctx.projector = geoProjector(opts, (lo + hi) / 2, (la + lb) / 2));
    for (let k = 0; k < n; k++) { const p = proj.forward(xs[at(k)], ys[at(k)]); out[3 * k] = p[0]; out[3 * k + 1] = p[1]; }
    crs = { source: 'EPSG:4326', basis: r.lon >= 0 ? 'headers' : 'option', ...proj.crs };
  } else {
    const sx = colScale(t.cols[xi], units), sy = colScale(t.cols[yi], units);
    for (let k = 0; k < n; k++) { out[3 * k] = xs[at(k)] * sx; out[3 * k + 1] = ys[at(k)] * sy; }
    if (h.crsOpt) crs = { source: String(h.crsOpt), basis: opts.crs ? 'option' : 'file', units: 'm' };
  }
  if (zc) { const sz = colScale(zc, units) * vert.sign; for (let k = 0; k < n; k++) out[3 * k + 2] = zc.num[at(k)] * sz; }
  const roles = { ...r };
  return { coords: out, roles, crs, vertical: vert ? { ...vert, column: zc.name } : null, zi };
}

/** Try to interpret scattered points as a regular lattice. Returns a grid model or null. */
export function gridFromPoints(name, c, tolRel = 1e-4) {
  const n = c.length / 3;
  if (n < 4) return null;
  const axis = (o) => {
    const v = new Float64Array(n);
    for (let i = 0; i < n; i++) v[i] = c[3 * i + o];
    v.sort();
    const range = v[n - 1] - v[0];
    if (!(range > 0)) return null;
    // the pitch is the smallest gap that is not rounding noise
    let pitch = Infinity;
    for (let i = 1; i < n; i++) { const d = v[i] - v[i - 1]; if (d > range * 1e-9 && d < pitch) pitch = d; }
    const cnt = Math.round(range / pitch) + 1;
    if (!(cnt >= 2) || cnt > n || Math.abs(range / pitch - (cnt - 1)) > tolRel * cnt) return null;
    return { min: v[0], pitch, cnt };
  };
  const ax = axis(0), ay = axis(1);
  if (!ax || !ay || ax.cnt * ay.cnt > 2 * n || ax.cnt * ay.cnt < n) return null;
  const z = new Float32Array(ax.cnt * ay.cnt).fill(NaN);
  const seen = new Uint8Array(z.length);
  for (let i = 0; i < n; i++) {
    const fx = (c[3 * i] - ax.min) / ax.pitch, fy = (c[3 * i + 1] - ay.min) / ay.pitch;
    const ix = Math.round(fx), iy = Math.round(fy);
    if (Math.abs(fx - ix) > 0.01 || Math.abs(fy - iy) > 0.01) return null;
    const k = iy * ax.cnt + ix;
    if (seen[k]) return null;
    seen[k] = 1;
    z[k] = c[3 * i + 2];
  }
  return { kind: 'grid', name, nx: ax.cnt, ny: ay.cnt, x0: ax.min, y0: ay.min, dx: ax.pitch, dy: ay.pitch, z };
}

/** Column roles of a table without a header row: [id] x y [z], by position. */
function positionalMapping(t) {
  const n = t.cols.length;
  const firstIsId = n >= 4 && isIntCol(t.cols[0]) && !isIntCol(t.cols[1]) && new Set(t.cols[0].num).size === t.n;
  const o = firstIsId ? 1 : 0;
  const m = { x: o, y: o + 1 };
  if (n - o >= 3) m.z = o + 2;
  if (firstIsId) m.id = 0;
  return m;
}

function buildPoints(t, cls, opts, units, warnings, ctx, out) {
  const res = coordsOf(t, opts, units, warnings, ctx);
  if (!res) return false;
  const roles = { ...res.roles, ...pick(t, ['id'], opts.mapping) };
  const { attrs, units: attrUnits } = restAttrs(t, roles);
  if (roles.id >= 0) attrs.names = Array.from({ length: t.n }, (_, i) => t.cols[roles.id].text(i));
  const wantGrid = cls.schema === 'points' && cls.grid !== false && opts.grid !== false && res.zi >= 0 && !Object.keys(attrs).length;
  const g = wantGrid ? gridFromPoints(t.name, res.coords) : null;
  if (cls.grid === true && !g) warnings.push(`Table "${t.name}": the points do not form a regular lattice, so they are returned as points rather than a grid.`);
  if (g) {
    if (res.crs) g.crs = res.crs.source;
    out.models.push(g);
    out.schema.push({ table: t.name, type: 'grid', basis: cls.basis });
    tally(out.hist, 'gridCell', g.nx * g.ny);
  } else {
    out.models.push(pointsModel(t.name, res.coords, attrs));
    out.schema.push({ table: t.name, type: cls.schema === 'layout' ? 'layout' : 'points', basis: cls.basis });
    tally(out.hist, 'point1', t.n);
  }
  out.mapping[t.name] = describe(t, roles, units, Object.keys(attrUnits).length ? { attrUnits } : null);
  if (res.crs) out.crs = res.crs;
  if (res.vertical) out.vertical = res.vertical;
  return true;
}

function buildRoute(t, cls, opts, units, warnings, ctx, out) {
  const extra = pick(t, ['kp', 'line'], opts.mapping);
  let res = coordsOf(t, opts, units, warnings, ctx);
  let profile = false;
  if (!res) {
    // profile only: chainage + depth/elevation → a line in the x–z plane; or x + z columns
    const xz = t.cols.some((c) => c.key === 'x') && t.cols.some((c) => c.key === 'z');
    if (extra.kp < 0 && !xz) return false;
    const v = pick(t, ['z', 'depth'], opts.mapping);
    const zi = v.z >= 0 ? v.z : v.depth;
    if (zi < 0) return false;
    profile = true;
    const xcol = extra.kp >= 0 ? t.cols[extra.kp] : t.cols.find((c) => c.key === 'x');
    const s = extra.kp >= 0 ? kpScale(xcol, units) : colScale(xcol, units);
    const vert = verticalConvention(opts, v.z >= 0 ? 'z' : t.cols[zi].key);
    const c = new Float64Array(t.n * 3);
    const sz = colScale(t.cols[zi], units) * vert.sign;
    for (let i = 0; i < t.n; i++) { c[3 * i] = xcol.num[i] * s; c[3 * i + 2] = t.cols[zi].num[i] * sz; }
    res = { coords: c, roles: { ...v, x: extra.kp >= 0 ? -1 : t.cols.indexOf(xcol) }, crs: null, vertical: { ...vert, column: t.cols[zi].name } };
  }
  const roles = { ...res.roles, ...extra };
  const lineCol = roles.line >= 0 ? t.cols[roles.line] : null;
  const groups = new Map();
  if (lineCol) for (let i = 0; i < t.n; i++) { const k = lineCol.text(i); let g = groups.get(k); if (!g) groups.set(k, (g = [])); g.push(i); } else groups.set(t.name, null);
  for (const [gname, rows] of groups) {
    const n = rows ? rows.length : t.n;
    const pts = new Float64Array(n * 3);
    let dropped = 0, k = 0;
    const keep = [];
    for (let j = 0; j < n; j++) {
      const i = rows ? rows[j] : j;
      const x = res.coords[3 * i], y = res.coords[3 * i + 1], z = res.coords[3 * i + 2];
      if (!(x === x && y === y && z === z)) { dropped++; continue; }
      pts[3 * k] = x; pts[3 * k + 1] = y; pts[3 * k + 2] = z; k++; keep.push(i);
    }
    if (dropped) warnings.push(`Table "${t.name}": ${dropped} row(s) with missing coordinates were dropped from "${gname}".`);
    if (k < 2) { warnings.push(`Table "${t.name}": "${gname}" has fewer than two valid points and was skipped.`); continue; }
    const { attrs, units: attrUnits } = restAttrs(t, roles, keep);
    const points = k === n ? pts : pts.slice(0, 3 * k);
    if (roles.kp >= 0) {
      const s = kpScale(t.cols[roles.kp], units);
      attrs.kp = Float64Array.from(keep, (i) => t.cols[roles.kp].num[i] * s);
    } else {
      // no stationing in the file: cumulative 3-D length along the line
      const kp = new Float64Array(k);
      for (let i = 1; i < k; i++) kp[i] = kp[i - 1] + Math.hypot(points[3 * i] - points[3 * i - 3], points[3 * i + 1] - points[3 * i - 2], points[3 * i + 2] - points[3 * i - 1]);
      attrs.kp = kp;
    }
    const model = { kind: 'polyline', name: cleanName(gname, t.name), points, attrs };
    if (Object.keys(attrUnits).length) model.attrUnits = attrUnits;
    out.models.push(model);
    tally(out.hist, 'polyline');
  }
  out.schema.push({ table: t.name, type: profile ? 'route-profile' : 'route', basis: cls.basis });
  const m = describe(t, roles, units);
  if (m.kp) { m.kp.scale = kpScale(t.cols[roles.kp], units); m.kp.note = 'attrs.kp is in metres along the line'; } else m.kp = { column: null, note: 'attrs.kp computed as cumulative 3-D length, metres' };
  out.mapping[t.name] = m;
  if (res.crs) out.crs = res.crs;
  if (res.vertical) out.vertical = res.vertical;
  return true;
}

// "KP" is kilometre point by definition; chainage, station and distance columns follow the file's length unit.
function kpScale(c, units) {
  const u = c.unit ? normaliseUnit(c.unit) : null;
  if (u) return UNIT_SCALE[u];
  return c.key === 'kp' ? 1000 : units.scale;
}

function sectionOf(t, roles, i, units, opts, stats) {
  const dimU = opts.sectionUnits ? normaliseUnit(opts.sectionUnits) : null;
  if (opts.sectionUnits && !dimU) throw new RangeError(`sectionUnits: unknown length unit "${opts.sectionUnits}".`);
  const sc = (ci) => { const c = t.cols[ci]; const u = c.unit ? normaliseUnit(c.unit) : null; return u ? UNIT_SCALE[u] : dimU ? UNIT_SCALE[dimU] : units.scale; };
  if (roles.D >= 0) {
    const D = t.cols[roles.D].num[i] * sc(roles.D);
    let tt = roles.t >= 0 ? t.cols[roles.t].num[i] * sc(roles.t) : NaN;
    if (!(tt === tt) && roles.ID >= 0) tt = (D - t.cols[roles.ID].num[i] * sc(roles.ID)) / 2;
    if (D > 0 && tt > 0 && 2 * tt <= D) { if (D > 20 && !t.cols[roles.D].unit && !dimU) stats.bigD++; return { type: 'tube', D, t: tt }; }
    if (D > 0 && !(tt === tt)) return { type: 'circle', D, A: (Math.PI * D * D) / 4, Iy: (Math.PI * D ** 4) / 64, Iz: (Math.PI * D ** 4) / 64, J: (Math.PI * D ** 4) / 32 };
  }
  if (roles.A >= 0 && t.cols[roles.A].num[i] > 0) {
    const s = sc(roles.A), v = (r, p) => (roles[r] >= 0 ? t.cols[roles[r]].num[i] * s ** p : 0);
    return { A: v('A', 2), Iy: v('Iy', 4), Iz: v('Iz', 4), J: v('J', 4) };
  }
  if (roles.section >= 0) {
    const nm = t.cols[roles.section].text(i);
    if (nm) {
      const r = resolveSection(nm, { dimUnit: opts.sectionUnits || 'mm' });
      if (r.matched) return { ...r.section, name: nm };
      stats.unresolved.add(nm);
      return { name: nm };
    }
  }
  return null;
}

function buildFrame(tables, classes, opts, units, warnings, ctx, out) {
  const nodes = new Grow(Float64Array, 256), members = [];
  const nodeIndex = new Map(), nodeNames = [];
  const stats = { bigD: 0, unresolved: new Set() };
  const weld = new Welder(0);
  let missing = 0, zero = 0;
  const norm = (s) => { const v = Number(s); return s !== '' && v === v ? String(v) : String(s).trim(); };
  tables.forEach((t, k) => {
    if (classes[k].schema !== 'frame' || classes[k].part !== 'nodes') return;
    const res = coordsOf(t, opts, units, warnings, ctx);
    if (!res) throw new Error(`Table "${t.name}": a node table needs x and y (and usually z) columns; found ${t.cols.map((c) => c.name).join(', ')}.`);
    const id = pick(t, ['id'], opts.mapping).id;
    for (let i = 0; i < t.n; i++) {
      const key = id >= 0 ? norm(t.cols[id].text(i)) : String(i + 1);
      if (nodeIndex.has(key)) { warnings.push(`Table "${t.name}": node id "${key}" is defined more than once; the first definition is kept.`); continue; }
      nodeIndex.set(key, nodes.n / 3);
      nodeNames.push(key);
      nodes.push3(res.coords[3 * i], res.coords[3 * i + 1], res.coords[3 * i + 2]);
    }
    out.mapping[t.name] = describe(t, { ...res.roles, id }, units);
    out.schema.push({ table: t.name, type: 'frame-nodes', basis: classes[k].basis });
    if (res.crs) out.crs = res.crs;
  });
  const hadNodeTable = nodes.n > 0;
  tables.forEach((t, k) => {
    const cls = classes[k];
    if (cls.schema !== 'frame' || cls.part === 'nodes') return;
    const common = pick(t, ['D', 't', 'ID', 'section', 'material', 'group', 'A', 'Iy', 'Iz', 'J'], opts.mapping);
    let roles;
    if (cls.part === 'inline') {
      const e = pick(t, ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], opts.mapping);
      if (e.x1 < 0 || e.y1 < 0 || e.x2 < 0 || e.y2 < 0) throw new Error(`Table "${t.name}": member end coordinates need columns x1, y1, z1, x2, y2, z2.`);
      const idc = pick(t, ['id'], opts.mapping).id;
      roles = { ...e, ...common, id: idc };
      const s = (ci) => (ci >= 0 ? colScale(t.cols[ci], units) : 0);
      const v = (ci, i) => (ci >= 0 ? t.cols[ci].num[i] * s(ci) : 0);
      for (let i = 0; i < t.n; i++) {
        const base = hadNodeTable ? nodes.n / 3 : 0;
        let a, b;
        if (hadNodeTable) { a = base; nodes.push3(v(e.x1, i), v(e.y1, i), v(e.z1, i)); b = base + 1; nodes.push3(v(e.x2, i), v(e.y2, i), v(e.z2, i)); } else { a = weld.add(v(e.x1, i), v(e.y1, i), v(e.z1, i)); b = weld.add(v(e.x2, i), v(e.y2, i), v(e.z2, i)); }
        if (a === b) { zero++; continue; }
        members.push(memberRec(t, roles, i, a, b, units, opts, stats));
      }
    } else {
      const used = new Set();
      let ai, bi;
      if (cls.pair && !(opts.mapping && (opts.mapping.a !== undefined || opts.mapping.b !== undefined))) { ai = t.cols.findIndex((c) => c.key === cls.pair[0]); bi = t.cols.findIndex((c) => c.key === cls.pair[1]); } else { ai = findRole(t, 'a', used, opts.mapping); used.add(ai); bi = findRole(t, 'b', used, opts.mapping); }
      if (ai < 0 || bi < 0) throw new Error(`Table "${t.name}": a member table needs two node-reference columns (for example node_a / node_b, n1 / n2, start / end).`);
      const skip = new Set([ai, bi]);
      const idc = ROLE_NAMES.id.map((nm) => t.cols.findIndex((c, q) => c.key === nm && !skip.has(q))).find((q) => q >= 0) ?? -1;
      // D / t etc. must not grab the end-node columns ('t' vs 'to', 'a' vs area)
      for (const r in common) if (skip.has(common[r])) common[r] = -1;
      roles = { a: ai, b: bi, ...common, id: idc };
      if (!nodeIndex.size) throw new Error(`Table "${t.name}" lists members by node id but no node table (id, x, y, z) was found. Put both tables in the file (separated by [nodes] / [members] lines) or upload the node file with it.`);
      for (let i = 0; i < t.n; i++) {
        const a = nodeIndex.get(norm(t.cols[ai].text(i))), b = nodeIndex.get(norm(t.cols[bi].text(i)));
        if (a === undefined || b === undefined) { missing++; continue; }
        if (a === b) { zero++; continue; }
        members.push(memberRec(t, roles, i, a, b, units, opts, stats));
      }
    }
    out.mapping[t.name] = describe(t, roles, units);
    out.schema.push({ table: t.name, type: cls.part === 'inline' ? 'frame-members-inline' : 'frame-members', basis: cls.basis });
  });
  const xyz = hadNodeTable ? nodes.finish() : weld.finish();
  if (!members.length) {
    if (!xyz.length) return false;
    warnings.push('The node table was read but no members could be built; the nodes are returned as points.');
    out.models.push(pointsModel(ctx.name, xyz, { names: nodeNames }));
    return true;
  }
  if (missing) warnings.push(`${missing} member(s) reference node ids that are not in the node table and were skipped.`);
  if (zero) warnings.push(`${zero} zero-length member(s) (both ends on the same node) were skipped.`);
  if (stats.bigD) warnings.push(`${stats.bigD} member diameter(s) exceed 20 m as read: the D / t columns carry no unit and were taken in the file's length unit. Put the unit in the header (for example "D (mm)") or set sectionUnits.`);
  if (stats.unresolved.size) warnings.push(`Section name(s) not resolved to properties (kept by name only): ${[...stats.unresolved].slice(0, 12).join(', ')}${stats.unresolved.size > 12 ? ' …' : ''}.`);
  const model = { kind: 'frame', name: ctx.name, nodes: xyz, members };
  if (hadNodeTable) model.nodeNames = nodeNames.concat(Array.from({ length: xyz.length / 3 - nodeNames.length }, () => ''));
  out.models.push(model);
  tally(out.hist, 'line2', members.length);
  if (missing) out.skipped.membersWithMissingNodes = missing;
  if (zero) out.skipped.zeroLengthMembers = zero;
  return true;
}

function memberRec(t, roles, i, a, b, units, opts, stats) {
  const m = { a, b };
  const sec = sectionOf(t, roles, i, units, opts, stats);
  if (sec) { if (sec.type === 'tube' || sec.A !== undefined) m.section = sec; else m.sectionName = sec.name; }
  if (roles.material >= 0 && t.cols[roles.material].text(i)) m.material = t.cols[roles.material].text(i);
  if (roles.group >= 0 && t.cols[roles.group].text(i)) m.group = t.cols[roles.group].text(i);
  if (roles.id >= 0) m.id = t.cols[roles.id].text(i);
  return m;
}

/**
 * Hull offsets → surface. Each station is a transverse section (half-breadth y against height z). Stations
 * with the full number of offsets are used as given (knuckles survive); shorter ones are resampled to that
 * number by arc length. Neighbouring stations are joined by quad strips. The starboard half
 * (y ≥ 0) is mirrored about the centre plane y = 0; points on the centre plane are shared by both halves.
 * Deck, transom and bow are left open (the table does not define them).
 */
export function loftHull(name, stations) {
  const st = stations.filter((s) => s.y.length >= 1).sort((p, q) => p.x - q.x);
  if (st.length < 2) return null;
  const N = Math.max(2, ...st.map((s) => s.y.length));
  const ns = st.length;
  const half = new Float64Array(ns * N * 3);
  st.forEach((s, i) => {
    const n = s.y.length;
    const order = Array.from({ length: n }, (_, k) => k).sort((p, q) => s.z[p] - s.z[q] || s.y[p] - s.y[q]);
    const cum = new Float64Array(n);
    for (let k = 1; k < n; k++) cum[k] = cum[k - 1] + Math.hypot(s.y[order[k]] - s.y[order[k - 1]], s.z[order[k]] - s.z[order[k - 1]]);
    const L = cum[n - 1];
    let seg = 0;
    for (let j = 0; j < N; j++) {
      let y, z;
      if (n === N) { y = s.y[order[j]]; z = s.z[order[j]]; } else if (n === 1 || !(L > 0)) { y = s.y[order[0]]; z = s.z[order[0]]; } else {
        const target = (L * j) / (N - 1);
        while (seg < n - 2 && cum[seg + 1] < target) seg++;
        const d = cum[seg + 1] - cum[seg];
        const f = d > 0 ? Math.min(1, Math.max(0, (target - cum[seg]) / d)) : 0;
        y = s.y[order[seg]] + f * (s.y[order[seg + 1]] - s.y[order[seg]]);
        z = s.z[order[seg]] + f * (s.z[order[seg + 1]] - s.z[order[seg]]);
      }
      const o = 3 * (i * N + j);
      half[o] = s.x; half[o + 1] = Math.max(0, y); half[o + 2] = z;
    }
  });
  const weld = new Welder(0, ns * N * 2);
  const sb = new Uint32Array(ns * N), pt = new Uint32Array(ns * N);
  for (let k = 0; k < ns * N; k++) sb[k] = weld.add(half[3 * k], half[3 * k + 1], half[3 * k + 2]);
  for (let k = 0; k < ns * N; k++) pt[k] = weld.add(half[3 * k], -half[3 * k + 1], half[3 * k + 2]);
  const s1 = new Grow(Uint32Array, ns * N * 6), s2 = new Grow(Uint32Array, ns * N * 6);
  const tri = (g, a, b, c) => { if (a !== b && b !== c && a !== c) { g.push(a); g.push(b); g.push(c); } };
  for (let i = 0; i + 1 < ns; i++) for (let j = 0; j + 1 < N; j++) {
    const a = i * N + j, b = a + 1, c = a + N + 1, d = a + N;
    // starboard: outward normal +y; port: mirrored, so the winding is reversed
    tri(s1, sb[a], sb[b], sb[c]); tri(s1, sb[a], sb[c], sb[d]);
    tri(s2, pt[a], pt[c], pt[b]); tri(s2, pt[a], pt[d], pt[c]);
  }
  const indices = new Uint32Array(s1.n + s2.n);
  indices.set(s1.view()); indices.set(s2.view(), s1.n);
  return surfaceModel(name, weld.finish(), indices, [{ name: 'starboard', start: 0, count: s1.n }, { name: 'port', start: s1.n, count: s2.n }], {
    stations: st.map((s) => ({ x: s.x, y: Float64Array.from(s.y), z: Float64Array.from(s.z) })),
  });
}

function buildHull(t, cls, opts, units, warnings, ctx, out) {
  const stations = new Map();
  const add = (x, y, z) => { if (!(x === x && y === y && z === z)) return; let s = stations.get(x); if (!s) stations.set(x, (s = { x, y: [], z: [] })); s.y.push(y); s.z.push(z); };
  const s = units.scale;
  let roles = {};
  if (t.wide) {
    const cornerIsStation = /^(station|stn|sta|frame|x)/.test(keyOf(t.wide.corner));
    for (let i = 0; i < t.n; i++) for (let j = 1; j < t.cols.length; j++) {
      const a = t.cols[0].num[i], b = t.wide.values[j - 1], v = t.cols[j].num[i];
      if (cornerIsStation) add(a * s, v * s, b * s); else add(b * s, v * s, a * s);
    }
    out.mapping[t.name] = { layout: 'matrix', rows: cornerIsStation ? 'station (x)' : 'waterline (z)', columns: cornerIsStation ? 'waterline (z)' : 'station (x)', cells: 'half-breadth (y)' };
  } else {
    roles = pick(t, ['halfBreadth', 'station', 'waterline'], opts.mapping);
    if (roles.station < 0 || roles.halfBreadth < 0 || roles.waterline < 0) throw new Error(`Table "${t.name}": a hull offset table needs station (x), waterline / height (z) and half-breadth (y) columns; found ${t.cols.map((c) => c.name).join(', ')}.`);
    const sx = colScale(t.cols[roles.station], units), sy = colScale(t.cols[roles.halfBreadth], units), sz = colScale(t.cols[roles.waterline], units);
    for (let i = 0; i < t.n; i++) add(t.cols[roles.station].num[i] * sx, t.cols[roles.halfBreadth].num[i] * sy, t.cols[roles.waterline].num[i] * sz);
    out.mapping[t.name] = describe(t, roles, units, { layout: 'long' });
  }
  const model = loftHull(t.name, [...stations.values()]);
  if (!model) throw new Error(`Table "${t.name}": a hull offset table needs at least two stations with offsets.`);
  warnings.push('Hull lofted from offsets: sections are joined linearly between stations and mirrored about y = 0; deck, transom and bow closures are not generated, so the surface is open.');
  out.models.push(model);
  out.schema.push({ table: t.name, type: 'hull', basis: cls.basis, stations: stations.size });
  tally(out.hist, 'tri3', model.indices.length / 3);
  return true;
}

function buildBlade(t, cls, opts, units, warnings, ctx, out) {
  const roles = pick(t, ['chord', 'twist', 'r', 'relThickness', 'airfoil'], opts.mapping);
  if (roles.r < 0 || roles.chord < 0) throw new Error(`Table "${t.name}": a blade table needs a radius / span column and a chord column; found ${t.cols.map((c) => c.name).join(', ')}.`);
  const rs = colScale(t.cols[roles.r], units), cs = colScale(t.cols[roles.chord], units);
  const order = Array.from({ length: t.n }, (_, i) => i).filter((i) => t.cols[roles.r].num[i] === t.cols[roles.r].num[i]).sort((a, b) => t.cols[roles.r].num[a] - t.cols[roles.r].num[b]);
  const n = order.length;
  if (n < 2) throw new Error(`Table "${t.name}": a blade table needs at least two stations.`);
  const pts = new Float64Array(3 * n);
  const { attrs, units: attrUnits } = restAttrs(t, roles, order);
  attrs.r = new Float64Array(n); attrs.chord = new Float64Array(n);
  order.forEach((i, k) => { attrs.r[k] = t.cols[roles.r].num[i] * rs; attrs.chord[k] = t.cols[roles.chord].num[i] * cs; pts[3 * k + 2] = attrs.r[k]; });
  if (roles.twist >= 0) {
    const u = (t.cols[roles.twist].unit || '').toLowerCase();
    const f = u === 'rad' ? 180 / Math.PI : 1;
    attrs.twistDeg = Float64Array.from(order, (i) => t.cols[roles.twist].num[i] * f);
  }
  if (roles.relThickness >= 0 && isNumeric(t.cols[roles.relThickness])) attrs.relThickness = Float64Array.from(order, (i) => t.cols[roles.relThickness].num[i]);
  if (roles.airfoil >= 0) attrs.airfoil = order.map((i) => t.cols[roles.airfoil].text(i));
  const model = { kind: 'polyline', name: t.name, points: pts, attrs };
  if (Object.keys(attrUnits).length) model.attrUnits = attrUnits;
  out.models.push(model);
  out.schema.push({ table: t.name, type: 'blade', basis: cls.basis, note: 'span axis along +z; chord in metres, twist in degrees' });
  out.mapping[t.name] = describe(t, roles, units);
  tally(out.hist, 'polyline');
  return true;
}

function circleLoop(cx, cy, r, n = 72) {
  const l = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) { const a = (2 * Math.PI * i) / n; l[2 * i] = cx + r * Math.cos(a); l[2 * i + 1] = cy + r * Math.sin(a); }
  return l;
}

function buildSection(t, cls, opts, units, warnings, ctx, out) {
  const roles = pick(t, ['x', 'y', 'D', 'ID', 'thickness', 'material', 'id', 'group'], opts.mapping);
  // beside an "od" column, a numeric "id" column is the inner diameter, not an identifier
  if (roles.ID < 0 && roles.id >= 0 && t.cols[roles.id].key === 'id' && isNumeric(t.cols[roles.id]) && roles.D >= 0 && t.cols[roles.D].key === 'od') {
    roles.ID = roles.id;
    roles.id = t.cols.findIndex((c, k) => k !== roles.ID && ['name', 'label', 'tag', 'component', 'layer', 'element', 'item'].includes(c.key));
  }
  if (roles.D < 0 && roles.ID < 0) throw new Error(`Table "${t.name}": a cross-section table needs an outer-diameter column (od / d / diameter), optionally with id / thickness and x, y positions.`);
  const dimU = opts.sectionUnits ? normaliseUnit(opts.sectionUnits) : null;
  const sc = (ci) => { const c = t.cols[ci]; const u = c.unit ? normaliseUnit(c.unit) : null; return u ? UNIT_SCALE[u] : dimU ? UNIT_SCALE[dimU] : units.scale; };
  const loops = [], layers = [];
  for (let i = 0; i < t.n; i++) {
    const x = roles.x >= 0 ? t.cols[roles.x].num[i] * sc(roles.x) : 0, y = roles.y >= 0 ? t.cols[roles.y].num[i] * sc(roles.y) : 0;
    let od = roles.D >= 0 ? t.cols[roles.D].num[i] * sc(roles.D) : NaN;
    let id = roles.ID >= 0 ? t.cols[roles.ID].num[i] * sc(roles.ID) : NaN;
    const th = roles.thickness >= 0 ? t.cols[roles.thickness].num[i] * sc(roles.thickness) : NaN;
    if (!(od === od) && id === id && th === th) od = id + 2 * th;
    if (!(id === id) && th === th) id = od - 2 * th;
    if (!(od > 0)) continue;
    const layer = { name: roles.id >= 0 ? t.cols[roles.id].text(i) : `component ${i + 1}`, x: x || 0, y: y || 0, od, id: id > 0 ? id : 0, loops: [loops.length] };
    loops.push(circleLoop(layer.x, layer.y, od / 2));
    if (id > 0 && id < od) { layer.loops.push(loops.length); loops.push(circleLoop(layer.x, layer.y, id / 2)); }
    if (roles.material >= 0) layer.material = t.cols[roles.material].text(i);
    if (roles.group >= 0) layer.group = t.cols[roles.group].text(i);
    layers.push(layer);
  }
  if (!loops.length) return false;
  out.models.push({ kind: 'section', name: t.name, loops, layers });
  out.schema.push({ table: t.name, type: 'section', basis: cls.basis, note: 'circular components; each layer lists its outer (and inner) loop' });
  out.mapping[t.name] = describe(t, roles, units);
  tally(out.hist, 'loop', loops.length);
  return true;
}

function buildPanels(t, cls, opts, units, warnings, ctx, out) {
  let idx;
  if (t.headerless) {
    if (t.cols.length < 12) throw new Error(`Table "${t.name}": a panel table needs twelve numbers per row (x1 y1 z1 … x4 y4 z4); found ${t.cols.length} columns.`);
    idx = Array.from({ length: 12 }, (_, i) => i);
  } else {
    idx = [];
    for (let k = 1; k <= 4; k++) for (const a of ['x', 'y', 'z']) idx.push(t.cols.findIndex((c) => c.key === `${a}${k}`));
    if (idx.some((i) => i < 0)) throw new Error(`Table "${t.name}": a panel table needs the twelve columns x1 y1 z1 … x4 y4 z4.`);
  }
  const weld = new Welder(0, t.n * 4), tris = new Grow(Uint32Array, t.n * 6), quads = new Grow(Uint32Array, t.n * 4);
  const s = units.scale;
  let nq = 0, nt = 0, bad = 0;
  for (let i = 0; i < t.n; i++) {
    const v = [0, 0, 0, 0];
    let ok = true;
    for (let k = 0; k < 4; k++) { const x = t.cols[idx[3 * k]].num[i] * s, y = t.cols[idx[3 * k + 1]].num[i] * s, z = t.cols[idx[3 * k + 2]].num[i] * s; if (!(x === x && y === y && z === z)) ok = false; v[k] = weld.add(x, y, z); }
    const u = v.filter((p, k) => v.indexOf(p) === k);
    if (!ok || u.length < 3) { bad++; continue; }
    quads.push(v[0]); quads.push(v[1]); quads.push(v[2]); quads.push(v[3]);
    if (u.length === 3) { tris.push(u[0]); tris.push(u[1]); tris.push(u[2]); nt++; } else { tris.push(v[0]); tris.push(v[1]); tris.push(v[2]); tris.push(v[0]); tris.push(v[2]); tris.push(v[3]); nq++; }
  }
  if (!tris.n) return false;
  if (bad) { warnings.push(`${bad} degenerate or incomplete panel row(s) were skipped.`); out.skipped.degeneratePanels = bad; }
  out.models.push(surfaceModel(t.name, weld.finish(), tris.finish(), null, { quads: quads.finish() }));
  out.schema.push({ table: t.name, type: 'panels', basis: cls.basis });
  if (nq) tally(out.hist, 'quad4', nq);
  if (nt) tally(out.hist, 'tri3', nt);
  return true;
}

/**
 * Turn tables into models.
 * @param {object[]} tables from readDelimited() / makeTable()
 * @param {object} opts importer options (`mapping`, `schema`, `units`, `sectionUnits`, `vertical`, `crs`,
 *   `projection`, `origin`, `kinds`, `grid`)
 * @param {object} ctx `{ name, format, unitHint?, crsHint?, extra? }`
 */
export function modelsFromTables(tables, opts, ctx) {
  const warnings = [];
  if (!tables.length) throw new Error(`${ctx.name}: no data rows found (expected rows of delimited values, optionally under a header line).`);
  const headerUnit = (() => {
    for (const t of tables) for (const c of t.cols) { if (c.unit && ['x', 'y', 'z', 'easting', 'northing', 'x1', 'station'].includes(c.key) && normaliseUnit(c.unit)) return c.unit; }
    return null;
  })();
  const units = resolveUnits(opts, ctx.unitHint || headerUnit, 'm');
  const base = { ...opts, nameHint: ctx.name };
  let classes = tables.map((t) => classify(t, base, false, false));
  if (classes.some((c) => c.schema === 'frame' && c.part === 'members')) classes = tables.map((t, k) => (classes[k].schema === 'frame' ? classes[k] : classify(t, base, true, false)));
  if (classes.some((c) => c.schema === 'frame' && c.part === 'nodes')) classes = tables.map((t, k) => (classes[k].schema === 'frame' ? classes[k] : classify(t, base, true, true)));
  const out = { models: [], schema: [], mapping: {}, hist: {}, skipped: {}, crs: null, vertical: null };
  if (classes.some((c) => c.schema === 'frame')) buildFrame(tables, classes, base, units, warnings, ctx, out);
  tables.forEach((t, k) => {
    const cls = classes[k];
    // without a header row the columns are taken by position, whatever interpretation is asked for
    const o = t.headerless && !opts.mapping && t.cols.length >= 2 && cls.schema !== 'panels' ? { ...base, mapping: positionalMapping(t) } : base;
    if (t.ragged) warnings.push(`Table "${t.name}": ${t.ragged} row(s) have a different number of fields from the header; missing cells are treated as empty.`);
    let ok = true;
    switch (cls.schema) {
      case 'frame': break;
      case 'points': case 'layout': ok = buildPoints(t, cls, o, units, warnings, ctx, out); break;
      case 'route': ok = buildRoute(t, cls, o, units, warnings, ctx, out); break;
      case 'hull': ok = buildHull(t, cls, o, units, warnings, ctx, out); break;
      case 'blade': ok = buildBlade(t, cls, o, units, warnings, ctx, out); break;
      case 'section': ok = buildSection(t, cls, o, units, warnings, ctx, out); break;
      case 'panels': ok = buildPanels(t, cls, o, units, warnings, ctx, out); break;
      default: ok = false;
    }
    if (!ok) { out.skipped[`table:${t.name}`] = t.n; warnings.push(`Table "${t.name}" (${t.cols.map((c) => c.name).join(', ')}) was not recognised as geometry and was skipped. Map its columns with opts.mapping, for example { x: '…', y: '…', z: '…' }.`); }
  });
  if (!out.models.length) {
    const t = tables[0];
    throw new Error(`${ctx.name}: could not find coordinate columns. Columns present: ${t.cols.map((c) => c.name).join(', ')}. Expected headers such as x, y, z / easting, northing, depth / KP, or pass opts.mapping = { x: 'column', y: 'column', z: 'column' }.`);
  }
  const single = tables.length === 1 && out.schema.length === 1;
  const meta = {
    format: ctx.format, units, crs: out.crs, vertical: out.vertical, elementTypes: out.hist, skipped: out.skipped,
    schema: single ? out.schema[0] : out.schema, mapping: single ? out.mapping[tables[0].name] : out.mapping,
    tables: tables.map((t) => ({ name: t.name, rows: t.n, columns: t.cols.map((c) => c.name) })),
    ...(ctx.extra || {}),
  };
  const r = finish(out.models, warnings, meta);
  if (!out.crs) {
    // millimetre hints only make sense for structures; survey tables are routinely kilometres across
    const structural = out.schema.some((sc) => /frame|hull|panels|blade|section/.test(sc.type));
    unitAdvice(units, r.meta.bbox, structural ? suggestUnits : () => null, r.warnings);
  }
  return r;
}

/** Delimited text table (CSV / TSV / XYZ / whitespace). */
export function parseTable(input, opts = {}) {
  let tables = [], delimiter = null, comments = [];
  const sources = input.bytes || typeof input.text === 'string' ? [input] : [...(input.files || [])].filter(([k]) => ['csv', 'tsv', 'txt', 'xyz', 'dat', 'tab'].includes(extOf(k))).map(([k, v]) => ({ name: k, bytes: v }));
  if (!sources.length) throw new Error(`${input.name || 'table'}: no content supplied (expected delimited text).`);
  for (const src of sources) {
    const r = readDelimited(src);
    if (sources.length > 1) r.tables.forEach((t) => { if (/^table \d+$/.test(t.name) || t.name === baseName(src.name)) t.name = baseName(src.name); });
    tables = tables.concat(r.tables);
    delimiter = delimiter || r.delimiter;
    comments = comments.concat(r.comments);
  }
  let unitHint = null, crsHint = null;
  for (const c of comments) {
    const u = /\bunits?\s*[:=]\s*([A-Za-z"'-]+)/i.exec(c);
    if (u && !unitHint && normaliseUnit(u[1])) unitHint = u[1];
    const e = /\b(EPSG\s*:\s*\d{4,5})\b/i.exec(c);
    if (e && !crsHint) crsHint = e[1].replace(/\s/g, '').toUpperCase();
  }
  return modelsFromTables(tables, opts, { name: baseName(input.name), format: 'table', unitHint, crsHint, extra: { delimiter, comments: comments.slice(0, 10) } });
}
