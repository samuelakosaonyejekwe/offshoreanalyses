// ISO 10303-21 ("STEP physical file", Part 21) tokenizer, the kernel-free STEP geometry subsets built on it,
// and the CIMsteel Integration Standards (CIS/2) member reader. The IFC reader (ifc.js) uses the same
// tokenizer.
//
// STEP (AP203 / AP214 / AP242), read without a geometry kernel:
//   - tessellated geometry: TRIANGULATED_FACE, COMPLEX_TRIANGULATED_FACE, TRIANGULATED_SURFACE_SET,
//     COMPLEX_TRIANGULATED_SURFACE_SET (triangles, strips and fans over a COORDINATES_LIST), TESSELLATED_CURVE_SET
//   - faceted B-rep: faces bounded by POLY_LOOPs
//   - wireframe: POLYLINE, TRIMMED_CURVE on LINE / CIRCLE, B_SPLINE_CURVE_WITH_KNOTS (rational too),
//     COMPOSITE_CURVE, and EDGE_CURVEs of files that contain no faces; CARTESIAN_POINTs of geometric sets
//   - the length unit of the file
// Exact B-rep solids and surfaces (ADVANCED_FACE on analytic or spline surfaces) need a geometry kernel: the
// reader counts them and sets `meta.needsAdapter = 'occt'`. Assembly placements (NEXT_ASSEMBLY_USAGE_OCCURRENCE
// with ITEM_DEFINED_TRANSFORMATION) are not applied by this reader; a warning says so when they are present.
//
// CIS/2 (STRUCTURAL_FRAME_SCHEMA): analysis-model elements (ELEMENT_CURVE_SIMPLE joined to NODEs through
// ELEMENT_NODE_CONNECTIVITY) and design-model parts (LOCATED_PART placing a PART_PRISMATIC_SIMPLE of a given
// length along the local x axis) → frame, with SECTION_PROFILE names resolved through sections.js.
// Connections, bolts, welds, plates (PART_SHEET), loads and results are not read.

import { inputBytes, baseName, Grow, Welder, finish, surfaceModel, pointsModel, tally, V3, cleanName, groupsFromBuckets, triangulate3D, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits, normaliseUnit, UNIT_SCALE } from '../units.js';
import { sampleCurve, arcPoints } from './nurbs.js';
import { resolveSection } from './sections.js';

/** Reference to another entity instance (#id). */
export class Ref { constructor(id) { this.id = id; } }
/** Typed parameter such as IFCLABEL('x') or LENGTH_MEASURE(25.4). */
export class Typed { constructor(type, args) { this.type = type; this.args = args; } get value() { return this.args[0]; } }

function decodeP21String(s) {
  if (s.indexOf('\\') < 0) return s;
  return s
    .replace(/\\X2\\((?:[0-9A-Fa-f]{4})*)\\X0\\/g, (m, h) => { let o = ''; for (let i = 0; i < h.length; i += 4) o += String.fromCharCode(parseInt(h.slice(i, i + 4), 16)); return o; })
    .replace(/\\X4\\((?:[0-9A-Fa-f]{8})*)\\X0\\/g, (m, h) => { let o = ''; for (let i = 0; i < h.length; i += 8) o += String.fromCodePoint(parseInt(h.slice(i, i + 8), 16)); return o; })
    .replace(/\\X\\([0-9A-Fa-f]{2})/g, (m, h) => String.fromCharCode(parseInt(h, 16)))
    .replace(/\\S\\(.)/g, (m, c) => String.fromCharCode(c.charCodeAt(0) + 128))
    .replace(/\\\\/g, '\\');
}

/** Parse the text of a parameter list "( … )" starting at `pos` (which must be the opening parenthesis). */
function parseList(s, pos) {
  const out = [];
  let i = pos + 1;
  const n = s.length;
  for (;;) {
    while (i < n) { const c = s.charCodeAt(i); if (c === 32 || c === 44 || c === 10 || c === 13 || c === 9) i++; else break; }
    if (i >= n) return [out, i];
    const ch = s[i];
    if (ch === ')') return [out, i + 1];
    if (ch === '(') { const [l, e] = parseList(s, i); out.push(l); i = e; } else if (ch === '#') {
      let e = i + 1;
      while (e < n && s.charCodeAt(e) >= 48 && s.charCodeAt(e) <= 57) e++;
      out.push(new Ref(Number(s.slice(i + 1, e)))); i = e;
    } else if (ch === "'") {
      let e = i + 1, str = '';
      for (;;) { const q = s.indexOf("'", e); if (q < 0) { str += s.slice(e); e = n; break; } str += s.slice(e, q); if (s[q + 1] === "'") { str += "'"; e = q + 2; } else { e = q + 1; break; } }
      out.push(decodeP21String(str)); i = e;
    } else if (ch === '$') { out.push(null); i++; } else if (ch === '*') { out.push(undefined); i++; } else if (ch === '.') {
      const e = s.indexOf('.', i + 1);
      out.push(s.slice(i, e < 0 ? n : e + 1)); i = e < 0 ? n : e + 1;
    } else if (ch === '"') { const e = s.indexOf('"', i + 1); out.push(s.slice(i + 1, e < 0 ? n : e)); i = e < 0 ? n : e + 1; } else if ((ch >= '0' && ch <= '9') || ch === '-' || ch === '+') {
      let e = i + 1;
      while (e < n) { const c = s.charCodeAt(e); if ((c >= 48 && c <= 57) || c === 46 || c === 69 || c === 101 || c === 43 || c === 45) e++; else break; }
      out.push(Number(s.slice(i, e))); i = e;
    } else if (/[A-Za-z_]/.test(ch)) {
      let e = i + 1;
      while (e < n && /[A-Za-z0-9_]/.test(s[e])) e++;
      const type = s.slice(i, e).toUpperCase();
      while (e < n && /\s/.test(s[e])) e++;
      if (s[e] === '(') { const [l, e2] = parseList(s, e); out.push(new Typed(type, l)); i = e2; } else { out.push(type); i = e; }
    } else i++;
  }
}

/**
 * Index of the DATA section of a Part 21 file. Statements are located once (id, type, byte range) and their
 * parameter lists are parsed on first use, so files much larger than what is actually needed stay cheap.
 */
export class P21 {
  constructor(bytes, fileName = 'file') {
    this.bytes = bytes;
    this.fileName = fileName;
    const head = String.fromCharCode.apply(null, bytes.subarray(0, Math.min(bytes.length, 24))).replace(/^﻿|^ï»¿/, '');
    if (!/^\s*ISO-10303-21\s*;/.test(head)) throw new Error(`${fileName}: not an ISO 10303-21 file (expected the first line "ISO-10303-21;").`);
    this.ids = new Grow(Uint32Array, 4096); this.start = new Grow(Uint32Array, 4096); this.end = new Grow(Uint32Array, 4096); this.tcode = new Grow(Uint16Array, 4096);
    this.types = []; this.typeIndex = new Map(); this.byId = new Map(); this.cache = new Map(); this.typeLists = null;
    this.header = { description: '', name: '', schema: [], originatingSystem: '', timeStamp: '' };
    this.index();
  }
  slice(a, b) {
    const by = this.bytes;
    if (b - a < 8192) return String.fromCharCode.apply(null, by.subarray(a, b));
    let s = '';
    for (let p = a; p < b; p += 8192) s += String.fromCharCode.apply(null, by.subarray(p, Math.min(b, p + 8192)));
    return s;
  }
  index() {
    const b = this.bytes, n = b.length;
    // statement end: ';' outside strings and comments
    const stmtEnd = (p) => {
      for (; p < n; p++) {
        const c = b[p];
        if (c === 59) return p;
        if (c === 39) { for (p++; p < n; p++) { if (b[p] === 39) { if (b[p + 1] === 39) p++; else break; } } } else if (c === 47 && b[p + 1] === 42) { for (p += 2; p < n - 1; p++) if (b[p] === 42 && b[p + 1] === 47) { p++; break; } }
      }
      return n;
    };
    let p = 0, inData = false;
    while (p < n) {
      const c = b[p];
      if (c <= 32) { p++; continue; }
      if (c === 47 && b[p + 1] === 42) { for (p += 2; p < n - 1; p++) if (b[p] === 42 && b[p + 1] === 47) { p += 2; break; } continue; }
      const e = stmtEnd(p);
      if (c === 35 && inData) { // '#'
        let q = p + 1, id = 0;
        while (q < e && b[q] >= 48 && b[q] <= 57) { id = id * 10 + (b[q] - 48); q++; }
        while (q < e && (b[q] <= 32 || b[q] === 61)) q++;
        let t = q;
        while (t < e && b[t] !== 40 && b[t] > 32) t++;
        const type = b[q] === 40 ? '(' : this.slice(q, t).toUpperCase();
        let code = this.typeIndex.get(type);
        if (code === undefined) { code = this.types.length; this.types.push(type); this.typeIndex.set(type, code); }
        this.byId.set(id, this.ids.n);
        this.ids.push(id); this.start.push(b[q] === 40 ? q : t); this.end.push(e); this.tcode.push(code);
      } else {
        const kw = this.slice(p, Math.min(e, p + 24)).toUpperCase();
        if (kw.startsWith('DATA')) inData = true;
        else if (kw.startsWith('ENDSEC')) inData = false;
        else if (!inData && /^FILE_(SCHEMA|NAME|DESCRIPTION)/.test(kw)) {
          const s = this.slice(p, e);
          const [args] = parseList(s, s.indexOf('('));
          if (kw.startsWith('FILE_SCHEMA')) this.header.schema = (args[0] || []).map(String);
          else if (kw.startsWith('FILE_NAME')) { this.header.name = String(args[0] ?? ''); this.header.timeStamp = String(args[1] ?? ''); this.header.originatingSystem = String(args[5] ?? ''); } else this.header.description = (args[0] || []).join(' ');
        }
      }
      p = e + 1;
    }
  }
  get count() { return this.ids.n; }
  /** Type name of #id ('(' for a complex instance), or undefined. */
  typeOf(id) { const k = this.byId.get(id); return k === undefined ? undefined : this.types[this.tcode.a[k]]; }
  /** Parsed instance `{ id, type, args }`; complex instances have `parts: [{ type, args }]` and `args` of the first part. */
  get(id) {
    if (id instanceof Ref) id = id.id;
    let v = this.cache.get(id);
    if (v !== undefined) return v;
    const k = this.byId.get(id);
    if (k === undefined) return null;
    const type = this.types[this.tcode.a[k]];
    const s = this.slice(this.start.a[k], this.end.a[k]);
    if (type === '(') {
      const [items] = parseList(s, 0);
      const parts = items.filter((x) => x instanceof Typed).map((x) => ({ type: x.type, args: x.args }));
      v = { id, type: '(', parts, args: parts.length ? parts[0].args : [] };
    } else v = { id, type, args: s[0] === '(' ? parseList(s, 0)[0] : [] };
    if (this.cache.size > 400000) this.cache.clear();
    this.cache.set(id, v);
    return v;
  }
  /** Ids of all instances of the given type name(s). */
  idsOf(...names) {
    if (!this.typeLists) { this.typeLists = new Map(); for (let k = 0; k < this.ids.n; k++) { const t = this.tcode.a[k]; let l = this.typeLists.get(t); if (!l) this.typeLists.set(t, (l = [])); l.push(this.ids.a[k]); } }
    let out = [];
    for (const nm of names) { const c = this.typeIndex.get(nm); if (c !== undefined) out = out.concat(this.typeLists.get(c) || []); }
    return out;
  }
  /** Number of instances per type name. */
  histogram() { const h = {}; for (let k = 0; k < this.ids.n; k++) { const t = this.types[this.tcode.a[k]]; h[t] = (h[t] || 0) + 1; } return h; }
  /** Part of a (possibly complex) instance with the given type, or null. */
  part(inst, type) { if (!inst) return null; if (inst.type === type) return inst; if (inst.parts) return inst.parts.find((x) => x.type === type) || null; return null; }
  /** [x, y, z] of a CARTESIAN_POINT / DIRECTION-like instance reference. */
  xyz(ref) { const e = ref instanceof Ref ? this.get(ref.id) : ref; if (!e || !Array.isArray(e.args[e.args.length - 1])) return null; const c = e.args[e.args.length - 1]; return [Number(c[0]) || 0, Number(c[1]) || 0, Number(c[2]) || 0]; }
  /** Right-handed frame of an AXIS2_PLACEMENT_3D: `{ o, x, y, z }`. */
  placement(ref) {
    const e = ref instanceof Ref ? this.get(ref.id) : ref;
    const I = { o: [0, 0, 0], x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
    if (!e) return I;
    const refs = e.args.filter((a) => a instanceof Ref);
    const o = (refs[0] && this.xyz(refs[0])) || [0, 0, 0];
    // attributes: (name?,) location, axis, ref_direction — axis and ref_direction optional
    const tail = e.args.slice(e.args.findIndex((a) => a instanceof Ref) + 1);
    let z = tail[0] instanceof Ref ? this.xyz(tail[0]) : null, x = tail[1] instanceof Ref ? this.xyz(tail[1]) : null;
    z = z ? V3.norm(z) : [0, 0, 1];
    if (!x) x = Math.abs(z[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    const d = V3.dot(x, z);
    x = V3.norm([x[0] - d * z[0], x[1] - d * z[1], x[2] - d * z[2]]);
    if (V3.len(x) === 0) x = V3.norm(V3.cross([0, 1, 0], z));
    return { o, x, y: V3.cross(z, x), z };
  }
}

/** Length unit of a STEP / CIS/2 file: { name, scale (metres per unit) } or null. */
export function p21LengthUnit(p) {
  const SI = { '.EXA.': 1e18, '.PETA.': 1e15, '.TERA.': 1e12, '.GIGA.': 1e9, '.MEGA.': 1e6, '.KILO.': 1e3, '.HECTO.': 100, '.DECA.': 10, '.DECI.': 0.1, '.CENTI.': 0.01, '.MILLI.': 1e-3, '.MICRO.': 1e-6, '.NANO.': 1e-9 };
  const NAMES = { 1: 'm', 0.001: 'mm', 0.01: 'cm', 1000: 'km', 0.1: 'dm', 0.000001: 'um' };
  const unitOf = (e, depth = 0) => {
    if (!e || depth > 4) return null;
    const si = p.part(e, 'SI_UNIT');
    const isLen = e.type === '(' ? !!p.part(e, 'LENGTH_UNIT') : true;
    if (si && isLen) {
      const a = si.args.filter((x) => typeof x === 'string' || x === null || x === undefined);
      const nm = a.find((x) => typeof x === 'string' && /METRE/.test(x));
      if (!nm) return null;
      const pre = a.find((x) => typeof x === 'string' && SI[x]);
      const scale = pre ? SI[pre] : 1;
      return { name: NAMES[scale] || `${scale} m`, scale };
    }
    const cb = p.part(e, 'CONVERSION_BASED_UNIT');
    if (cb && isLen) {
      const nm = cb.args.find((x) => typeof x === 'string' && !x.startsWith('.'));
      const u = normaliseUnit(nm);
      if (u) return { name: u, scale: UNIT_SCALE[u] };
      const mref = cb.args.find((x) => x instanceof Ref);
      const m = mref && p.get(mref.id);
      if (m) { const val = m.args.find((x) => x instanceof Typed); const base = unitOf(p.get(m.args.find((x) => x instanceof Ref)), depth + 1); if (val && base) return { name: nm || 'custom', scale: Number(val.value) * base.scale }; }
    }
    return null;
  };
  for (const id of p.idsOf('(')) {
    const e = p.get(id);
    const g = p.part(e, 'GLOBAL_UNIT_ASSIGNED_CONTEXT');
    if (!g) continue;
    for (const r of (g.args[0] || [])) { if (!(r instanceof Ref)) continue; const u = p.get(r.id); if (u && (u.type !== '(' || p.part(u, 'LENGTH_UNIT'))) { const res = unitOf(u); if (res) return res; } }
  }
  for (const id of p.idsOf('(')) { const e = p.get(id); if (p.part(e, 'LENGTH_UNIT')) { const r = unitOf(e); if (r) return r; } }
  // files without a unit context (CIS/2): a plain SI length unit instance
  for (const id of p.idsOf('SI_UNIT')) { const r = unitOf(p.get(id)); if (r) return r; }
  return null;
}

// ---------------------------------------------------------------------------------------------------------
// STEP geometry subsets
// ---------------------------------------------------------------------------------------------------------

const BREP_TYPES = ['MANIFOLD_SOLID_BREP', 'BREP_WITH_VOIDS', 'SHELL_BASED_SURFACE_MODEL', 'ADVANCED_FACE', 'FACETED_BREP', 'CLOSED_SHELL', 'OPEN_SHELL'];

/**
 * STEP AP203 / AP214 / AP242 without a geometry kernel.
 * With `opts.allowEmpty` the function returns an empty model list (and `meta.needsAdapter`) instead of
 * throwing when the file holds only exact B-rep geometry.
 */
export function parseStep(input, opts = {}) {
  const fname = input.name || 'STEP';
  const p = new P21(inputBytes(input), fname);
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};
  const lu = p21LengthUnit(p);
  const units = resolveUnits(opts, lu ? lu.name : null, 'mm');
  const s = opts.units || opts.scale || !lu ? units.scale : lu.scale;
  if (lu && !opts.units && !opts.scale) { units.scale = lu.scale; units.unit = normaliseUnit(lu.name) || lu.name; units.basis = 'file'; }
  const models = [];

  // --- tessellated ---
  const coordCache = new Map();
  const coords = new Grow(Float64Array, 3072);
  const listBase = (ref) => {
    if (!(ref instanceof Ref)) return -1;
    let b = coordCache.get(ref.id);
    if (b !== undefined) return b;
    const e = p.get(ref.id);
    const pts = e && e.args.find((a) => Array.isArray(a) && Array.isArray(a[0]));
    if (!pts) { coordCache.set(ref.id, -1); return -1; }
    b = coords.n / 3;
    for (const q of pts) coords.push3(q[0] * s, q[1] * s, (q[2] || 0) * s);
    coordCache.set(ref.id, b);
    return b;
  };
  const owner = new Map();
  for (const id of p.idsOf('TESSELLATED_SOLID', 'TESSELLATED_SHELL')) { const e = p.get(id); const nm = cleanName(e.args[0], `${e.type.toLowerCase().replace('_', ' ')} #${id}`); for (const r of e.args[1] || []) if (r instanceof Ref) owner.set(r.id, nm); }
  const buckets = new Map();
  let badTri = 0;
  for (const id of p.idsOf('TRIANGULATED_FACE', 'COMPLEX_TRIANGULATED_FACE', 'TRIANGULATED_SURFACE_SET', 'COMPLEX_TRIANGULATED_SURFACE_SET')) {
    const e = p.get(id);
    const complex = e.type.startsWith('COMPLEX');
    const base = listBase(e.args[1]);
    if (base < 0) { tally(skipped, `${e.type} without a coordinates list`); continue; }
    const lists = e.args.filter(Array.isArray);
    // … pnindex, triangles   |   … pnindex, triangle_strips, triangle_fans   (normals precede and are lists of lists of reals)
    const tail = complex ? lists.slice(-3) : lists.slice(-2);
    const pn = tail.length === (complex ? 3 : 2) && tail[0].every((v) => typeof v === 'number') ? tail[0] : [];
    const at = (k) => base + (pn.length ? pn[k - 1] : k) - 1;
    const key = owner.get(id) || cleanName(e.args[0], name);
    let g = buckets.get(key);
    if (!g) buckets.set(key, (g = new Grow(Uint32Array, 3072)));
    const nV = coords.n / 3;
    const tri = (a, b, c) => { const i = at(a), j = at(b), k = at(c); if (!(i >= 0 && j >= 0 && k >= 0 && i < nV && j < nV && k < nV)) { badTri++; return; } if (i === j || j === k || i === k) return; g.push(i); g.push(j); g.push(k); };
    if (!complex) { for (const t of lists[lists.length - 1] || []) if (Array.isArray(t) && t.length >= 3) tri(t[0], t[1], t[2]); } else {
      for (const st of lists[lists.length - 2] || []) if (Array.isArray(st)) for (let i = 0; i + 2 < st.length; i++) { if (i & 1) tri(st[i + 1], st[i], st[i + 2]); else tri(st[i], st[i + 1], st[i + 2]); }
      for (const fan of lists[lists.length - 1] || []) if (Array.isArray(fan)) for (let i = 1; i + 1 < fan.length; i++) tri(fan[0], fan[i], fan[i + 1]);
    }
    tally(hist, e.type);
  }
  if (badTri) { skipped.trianglesWithInvalidIndices = badTri; warnings.push(`${badTri} tessellated triangle(s) index points outside their coordinates list and were skipped.`); }

  // --- faceted B-rep (poly loops) ---
  let polyFaces = 0, edgeFaces = 0;
  {
    const faceTypes = p.idsOf('FACE_SURFACE', 'ADVANCED_FACE', 'FACE');
    const fc = new Grow(Float64Array, 3072), fi = new Grow(Uint32Array, 3072);
    for (const id of faceTypes) {
      const e = p.get(id);
      const bounds = (e.args.find(Array.isArray) || []).filter((r) => r instanceof Ref).map((r) => p.get(r.id)).filter(Boolean);
      const rings = [];
      let poly = bounds.length > 0;
      for (const bnd of bounds.sort((a, b) => (a.type === 'FACE_OUTER_BOUND' ? -1 : 0) - (b.type === 'FACE_OUTER_BOUND' ? -1 : 0))) {
        const loop = p.get(bnd.args.find((a) => a instanceof Ref));
        if (!loop || loop.type !== 'POLY_LOOP') { poly = false; break; }
        const pts = (loop.args.find(Array.isArray) || []).map((r) => p.xyz(r)).filter(Boolean);
        if (pts.length >= 3) rings.push(bnd.args.includes('.F.') ? pts.reverse() : pts);
      }
      if (!poly) { edgeFaces++; continue; }
      if (!rings.length) continue;
      const total = rings.reduce((a, r) => a + r.length, 0);
      const flat = new Float64Array(total * 3), holes = [];
      let o = 0;
      rings.forEach((r, k) => { if (k) holes.push(o); for (const q of r) { flat[3 * o] = q[0] * s; flat[3 * o + 1] = q[1] * s; flat[3 * o + 2] = q[2] * s; o++; } });
      const t = triangulate3D(flat, holes.length ? holes : undefined);
      const base = fc.n / 3;
      for (const v of flat) fc.push(v);
      for (const k of t) fi.push(base + k);
      polyFaces++;
    }
    if (fi.n) {
      // weld the per-face copies of the shared corners
      const w = new Welder(0, fc.n / 3);
      const map = new Uint32Array(fc.n / 3);
      for (let i = 0; i < map.length; i++) map[i] = w.add(fc.a[3 * i], fc.a[3 * i + 1], fc.a[3 * i + 2]);
      const idx = fi.finish();
      for (let i = 0; i < idx.length; i++) idx[i] = map[idx[i]];
      models.push(surfaceModel(buckets.size ? `${name} (faceted)` : name, w.finish(), idx));
      tally(hist, 'POLY_LOOP face', polyFaces);
    }
  }
  if (coords.n && [...buckets.values()].some((g) => g.n)) {
    const { indices, groups } = groupsFromBuckets(buckets);
    models.push(surfaceModel(name, coords.finish(), indices, groups.length > 1 ? groups : null));
  }

  // --- curves ---
  const angleFactor = (() => { for (const id of p.idsOf('(')) { const e = p.get(id); if (p.part(e, 'PLANE_ANGLE_UNIT') && p.part(e, 'CONVERSION_BASED_UNIT')) { const nm = p.part(e, 'CONVERSION_BASED_UNIT').args.find((x) => typeof x === 'string' && !x.startsWith('.')); if (/DEG/i.test(nm || '')) return Math.PI / 180; } } return 1; })();
  const curve = (ref, depth = 0) => {
    const e = ref instanceof Ref ? p.get(ref.id) : ref;
    if (!e || depth > 8) return null;
    const bs = p.part(e, 'B_SPLINE_CURVE_WITH_KNOTS');
    if (bs) {
      const base = e.type === '(' ? p.part(e, 'B_SPLINE_CURVE') : e;
      const o = e.type === '(' ? 0 : 1;
      const deg = base.args[o], cps = (base.args[o + 1] || []).map((r) => p.xyz(r)).filter(Boolean);
      const kl = bs.args.filter(Array.isArray).filter((a) => a.every((v) => typeof v === 'number'));
      const mults = kl[kl.length - 2] || [], kn = kl[kl.length - 1] || [];
      const U = [];
      kn.forEach((k, i) => { for (let m = 0; m < (mults[i] || 1); m++) U.push(k); });
      const rat = p.part(e, 'RATIONAL_B_SPLINE_CURVE');
      const W = rat ? rat.args.find(Array.isArray) : null;
      const pts = sampleCurve(deg, U, Float64Array.from(cps.flat()), W && W.length === cps.length ? W : null, 8);
      if (pts) tally(hist, 'B_SPLINE_CURVE');
      return pts;
    }
    switch (e.type) {
      case 'POLYLINE': { const pts = (e.args.find(Array.isArray) || []).map((r) => p.xyz(r)).filter(Boolean); tally(hist, 'POLYLINE'); return pts.length >= 2 ? Float64Array.from(pts.flat()) : null; }
      case 'TRIMMED_CURVE': {
        const basis = p.get(e.args.find((a) => a instanceof Ref));
        const trims = e.args.filter(Array.isArray);
        const sense = !e.args.includes('.F.');
        const pick = (t) => { const out = { pt: null, par: null }; for (const v of t || []) { if (v instanceof Ref) out.pt = p.xyz(v); else if (v instanceof Typed) out.par = Number(v.value); else if (typeof v === 'number') out.par = v; } return out; };
        const t1 = pick(trims[0]), t2 = pick(trims[1]);
        if (!basis) return null;
        if (basis.type === 'LINE') {
          const o = p.xyz(basis.args.find((a) => a instanceof Ref));
          const vec = p.get(basis.args.filter((a) => a instanceof Ref)[1]);
          const d = vec ? p.xyz(vec.args.find((a) => a instanceof Ref)) : null, mag = vec ? vec.args.find((a) => typeof a === 'number') : 1;
          const at = (t) => (t.pt ? t.pt : d && t.par !== null ? [o[0] + d[0] * mag * t.par, o[1] + d[1] * mag * t.par, o[2] + d[2] * mag * t.par] : null);
          const a = at(t1), b = at(t2);
          tally(hist, 'LINE');
          return a && b ? Float64Array.from([...a, ...b]) : null;
        }
        if (basis.type === 'CIRCLE') {
          const fr = p.placement(basis.args.find((a) => a instanceof Ref)), r = basis.args.find((a) => typeof a === 'number');
          const ang = (t) => (t.par !== null ? t.par * angleFactor : t.pt ? Math.atan2(V3.dot(V3.sub(t.pt, fr.o), fr.y), V3.dot(V3.sub(t.pt, fr.o), fr.x)) : null);
          let a0 = ang(t1), a1 = ang(t2);
          if (a0 === null || a1 === null) { a0 = 0; a1 = 2 * Math.PI; }
          let sw = a1 - a0;
          if (sense) { while (sw <= 1e-12) sw += 2 * Math.PI; } else while (sw >= -1e-12) sw -= 2 * Math.PI;
          tally(hist, 'CIRCLE');
          return arcPoints(fr.o, fr.x, fr.y, r, a0, sw);
        }
        const full = curve(basis, depth + 1);
        if (full) tally(skipped, `TRIMMED_CURVE on ${basis.type} drawn untrimmed`);
        return full;
      }
      case 'COMPOSITE_CURVE': case 'COMPOSITE_CURVE_ON_SURFACE': {
        const out = [];
        for (const sr of e.args.find(Array.isArray) || []) {
          const seg = p.get(sr);
          if (!seg) continue;
          let c = curve(seg.args.find((a) => a instanceof Ref), depth + 1);
          if (!c) continue;
          if (seg.args.includes('.F.')) { const r = new Float64Array(c.length); for (let i = 0; i < c.length; i += 3) { r[i] = c[c.length - 3 - i]; r[i + 1] = c[c.length - 2 - i]; r[i + 2] = c[c.length - 1 - i]; } c = r; }
          for (let i = 0; i < c.length; i += 3) { const n = out.length; if (n && Math.hypot(out[n - 3] - c[i], out[n - 2] - c[i + 1], out[n - 1] - c[i + 2]) < 1e-9) continue; out.push(c[i], c[i + 1], c[i + 2]); }
        }
        return out.length >= 6 ? Float64Array.from(out) : null;
      }
      case 'LINE': { const o = p.xyz(e.args.find((a) => a instanceof Ref)); const vec = p.get(e.args.filter((a) => a instanceof Ref)[1]); const d = vec ? p.xyz(vec.args.find((a) => a instanceof Ref)) : null; const mag = vec ? vec.args.find((a) => typeof a === 'number') : 1; tally(hist, 'LINE'); return o && d ? Float64Array.from([...o, o[0] + d[0] * mag, o[1] + d[1] * mag, o[2] + d[2] * mag]) : null; }
      case 'CIRCLE': { const fr = p.placement(e.args.find((a) => a instanceof Ref)); tally(hist, 'CIRCLE'); return arcPoints(fr.o, fr.x, fr.y, e.args.find((a) => typeof a === 'number'), 0, 2 * Math.PI); }
      case 'EDGE_CURVE': {
        const refs = e.args.filter((a) => a instanceof Ref);
        const vp = (r) => { const v = p.get(r); return v ? p.xyz(v.args.find((a) => a instanceof Ref)) : null; };
        const a = vp(refs[0]), b = vp(refs[1]), geom = p.get(refs[2]);
        const same = !e.args.includes('.F.');
        if (!a || !b) return null;
        if (geom && geom.type === 'CIRCLE') {
          const fr = p.placement(geom.args.find((x) => x instanceof Ref)), r = geom.args.find((x) => typeof x === 'number');
          const ang = (q) => Math.atan2(V3.dot(V3.sub(q, fr.o), fr.y), V3.dot(V3.sub(q, fr.o), fr.x));
          const a0 = ang(a);
          let sw = ang(b) - a0;
          if (same) { while (sw <= 1e-12) sw += 2 * Math.PI; } else while (sw >= -1e-12) sw -= 2 * Math.PI;
          tally(hist, 'CIRCLE');
          return arcPoints(fr.o, fr.x, fr.y, r, a0, sw);
        }
        if (geom && p.part(geom, 'B_SPLINE_CURVE_WITH_KNOTS')) { const c = curve(geom, depth + 1); if (c) return c; }
        if (geom && geom.type !== 'LINE') tally(skipped, `EDGE_CURVE on ${geom.type === '(' ? 'a complex curve' : geom.type} drawn as its chord`);
        else tally(hist, 'LINE');
        return Float64Array.from([...a, ...b]);
      }
      default: tally(skipped, `curve type ${e.type === '(' ? 'complex instance' : e.type}`); return null;
    }
  };
  const scaleC = (c) => { if (s !== 1) for (let i = 0; i < c.length; i++) c[i] *= s; return c; };
  const polylines = [], pts = [];
  const seen = new Set();
  for (const id of p.idsOf('GEOMETRIC_CURVE_SET', 'GEOMETRIC_SET')) {
    const e = p.get(id);
    for (const r of e.args.find(Array.isArray) || []) {
      if (!(r instanceof Ref) || seen.has(r.id)) continue;
      seen.add(r.id);
      const t = p.typeOf(r.id);
      if (t === 'CARTESIAN_POINT') { const q = p.xyz(r); if (q) pts.push(q[0] * s, q[1] * s, q[2] * s); continue; }
      const c = curve(r);
      if (c && c.length >= 6) polylines.push({ c: scaleC(c), name: cleanName(p.get(r.id).args[0], `${name} curve ${polylines.length + 1}`) });
    }
  }
  for (const id of p.idsOf('TESSELLATED_CURVE_SET')) {
    const e = p.get(id);
    const base = listBase(e.args[1]);
    if (base < 0) continue;
    const all = coords.a;
    for (const strip of e.args[e.args.length - 1] || []) { if (!Array.isArray(strip) || strip.length < 2) continue; const c = new Float64Array(strip.length * 3); strip.forEach((k, i) => { const o = 3 * (base + k - 1); c[3 * i] = all[o]; c[3 * i + 1] = all[o + 1]; c[3 * i + 2] = all[o + 2]; }); polylines.push({ c, name: cleanName(e.args[0], `${name} curve ${polylines.length + 1}`) }); tally(hist, 'TESSELLATED_CURVE_SET strip'); }
  }
  const brep = {};
  const h = p.histogram();
  for (const t of BREP_TYPES) if (h[t]) brep[t] = h[t];
  const hasFaces = edgeFaces > 0 || polyFaces > 0 || buckets.size > 0;
  if (!hasFaces && !polylines.length) {
    // a pure wireframe file: its edges are the geometry
    for (const id of p.idsOf('EDGE_CURVE')) { const c = curve(new Ref(id)); if (c && c.length >= 6) polylines.push({ c: scaleC(c), name: cleanName(p.get(id).args[0], `${name} edge ${polylines.length + 1}`) }); }
  }
  for (const l of polylines) models.push({ kind: 'polyline', name: l.name, points: l.c });
  if (pts.length) { models.push(pointsModel(`${name} (points)`, Float64Array.from(pts))); tally(hist, 'CARTESIAN_POINT', pts.length / 3); }

  const needs = edgeFaces > 0;
  if (needs) { skipped['B-rep faces (exact geometry)'] = edgeFaces; }
  if (h.NEXT_ASSEMBLY_USAGE_OCCURRENCE && models.length) warnings.push(`The file is an assembly (${h.NEXT_ASSEMBLY_USAGE_OCCURRENCE} component placements). Placements are not applied by the built-in reader: parts appear at their own coordinates.`);
  if (needs && models.length) warnings.push(`${edgeFaces} exact B-rep face(s) were not read: analytic and spline surfaces need the geometry-kernel reader (adapter "occt"). Only the tessellated / faceted / wireframe content is returned.`);
  for (const k in skipped) if (!/B-rep/.test(k)) warnings.push(`Skipped ${skipped[k]} × ${k}.`);
  const meta = {
    format: 'step', units, crs: null, elementTypes: hist, skipped, schema: p.header.schema, originatingSystem: p.header.originatingSystem, instances: p.count,
    products: p.idsOf('PRODUCT').slice(0, 200).map((id) => String(p.get(id).args[1] ?? p.get(id).args[0] ?? '')), brep, needsAdapter: needs ? 'occt' : null,
  };
  if (!models.length) {
    if (opts.allowEmpty) return finish([], warnings, meta);
    if (needs) throw new Error(`${fname}: the file holds exact B-rep geometry only (${Object.entries(brep).map(([k, v]) => `${v} ${k}`).join(', ')}), which needs the geometry-kernel reader (adapter "occt"). Without it, export the model as tessellated STEP (AP242 with tessellation), STL, OBJ or glTF.`);
    throw new Error(`${fname}: no geometry this reader understands was found (schema ${p.header.schema.join(', ') || 'unknown'}; expected tessellated faces, poly-loop faces, curves or points).`);
  }
  const r = finish(models, warnings, meta);
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}

// ---------------------------------------------------------------------------------------------------------
// CIS/2
// ---------------------------------------------------------------------------------------------------------

/** CIS/2 (CIMsteel) structural frame exchange file. */
export function parseCis2(input, opts = {}) {
  const fname = input.name || 'CIS/2';
  const p = new P21(inputBytes(input), fname);
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};
  const lu = p21LengthUnit(p);
  const units = resolveUnits(opts, lu ? lu.name : null, 'mm');
  if (lu && !opts.units && !opts.scale) { units.scale = lu.scale; units.unit = normaliseUnit(lu.name) || lu.name; units.basis = 'file'; }
  const s = units.scale;
  const dimUnit = ['in', 'ft', 'us-ft'].includes(units.unit) ? 'in' : 'mm';
  const refsOf = (e) => e.args.flat(3).filter((a) => a instanceof Ref);
  const firstOfType = (e, test) => { for (const r of refsOf(e)) { const t = p.typeOf(r.id); if (t && test(t)) return p.get(r.id); } return null; };
  const unresolved = new Set();
  const sectionFor = (e) => {
    const prof = firstOfType(e, (t) => /^SECTION_PROFILE/.test(t));
    if (!prof) return null;
    const nm = prof.args.find((a) => typeof a === 'string' && !a.startsWith('.') && a.trim() !== '' && !/^\d+$/.test(a.trim()));
    if (!nm) return null;
    const r = resolveSection(nm, { dimUnit });
    if (r.matched) return { section: { ...r.section, name: nm } };
    unresolved.add(nm);
    return { sectionName: nm };
  };
  const lengthOf = (e) => {
    for (const r of refsOf(e)) {
      const t = p.typeOf(r.id);
      if (!t || !/LENGTH_MEASURE_WITH_UNIT|MEASURE_WITH_UNIT/.test(t)) continue;
      const m = p.get(r.id);
      const v = m.args.find((a) => a instanceof Typed);
      if (v && Number(v.value) > 0) return Number(v.value);
    }
    return null;
  };
  const weld = new Welder(1e-6, 256);
  const members = [];
  const nodeNames = [];
  let zero = 0;

  // analysis model: elements joined to nodes
  const nodeOf = new Map();
  for (const id of p.idsOf('NODE')) {
    const e = p.get(id);
    const q = firstOfType(e, (t) => t === 'CARTESIAN_POINT');
    if (!q) continue;
    const c = p.xyz(q);
    const k = weld.add(c[0] * s, c[1] * s, c[2] * s);
    nodeOf.set(id, k);
    nodeNames[k] = String(e.args.find((a) => typeof a === 'string') ?? '');
  }
  const ends = new Map();
  for (const id of p.idsOf('ELEMENT_NODE_CONNECTIVITY')) {
    const e = p.get(id);
    let node = null, elem = null;
    for (const r of refsOf(e)) { const t = p.typeOf(r.id); if (t === 'NODE') node = r.id; else if (t && /^ELEMENT_CURVE/.test(t)) elem = r.id; }
    if (node === null || elem === null || !nodeOf.has(node)) continue;
    const order = e.args.find((a) => typeof a === 'number') ?? 0;
    let l = ends.get(elem);
    if (!l) ends.set(elem, (l = []));
    l.push({ order, node: nodeOf.get(node) });
  }
  for (const [elem, l] of ends) {
    l.sort((a, b) => a.order - b.order);
    if (l.length < 2 || l[0].node === l[l.length - 1].node) { zero++; continue; }
    const e = p.get(elem);
    const m = { a: l[0].node, b: l[l.length - 1].node, id: elem, group: 'element', name: String(e.args.find((a) => typeof a === 'string') ?? '') };
    Object.assign(m, sectionFor(e) || {});
    members.push(m);
    tally(hist, e.type);
  }

  // design model: located prismatic parts
  for (const id of p.idsOf('LOCATED_PART')) {
    const e = p.get(id);
    const part = firstOfType(e, (t) => /^PART_PRISMATIC/.test(t));
    const cs = firstOfType(e, (t) => /^COORD_SYSTEM/.test(t));
    if (!part) { tally(skipped, 'LOCATED_PART that is not a prismatic part'); continue; }
    const L = lengthOf(part);
    const pl = cs ? firstOfType(cs, (t) => t === 'AXIS2_PLACEMENT_3D') : null;
    if (!L || !pl) { tally(skipped, 'LOCATED_PART without a length or placement'); continue; }
    const fr = p.placement(pl);
    // CIS/2: the member's longitudinal axis is the local x axis of its coordinate system
    const a = weld.add(fr.o[0] * s, fr.o[1] * s, fr.o[2] * s), b = weld.add((fr.o[0] + fr.x[0] * L) * s, (fr.o[1] + fr.x[1] * L) * s, (fr.o[2] + fr.x[2] * L) * s);
    if (a === b) { zero++; continue; }
    const m = { a, b, id, group: 'part', name: String(e.args.find((x) => typeof x === 'string') ?? ''), orient: fr.z };
    Object.assign(m, sectionFor(part) || {});
    const mat = firstOfType(part, (t) => /^MATERIAL/.test(t));
    if (mat) { const nm = mat.args.find((x) => typeof x === 'string' && !x.startsWith('.')); if (nm) m.material = nm; }
    members.push(m);
    tally(hist, part.type);
  }
  const h = p.histogram();
  for (const t of ['PART_SHEET_BOUNDED_SIMPLE', 'PART_SHEET_BOUNDED_COMPLEX', 'JOINT_SYSTEM_WELDED', 'JOINT_SYSTEM_MECHANICAL', 'FASTENER_SIMPLE_BOLT', 'LOAD_CASE', 'ELEMENT_SURFACE_SIMPLE']) if (h[t]) skipped[t] = h[t];
  if (!members.length) throw new Error(`${fname}: no members were found (schema ${p.header.schema.join(', ') || 'unknown'}; expected ELEMENT_CURVE_SIMPLE with ELEMENT_NODE_CONNECTIVITY, or LOCATED_PART with PART_PRISMATIC_SIMPLE).`);
  if (zero) { skipped.zeroLengthMembers = zero; warnings.push(`${zero} zero-length or incompletely connected member(s) were skipped.`); }
  if (unresolved.size) warnings.push(`Section name(s) not resolved to properties (kept by name only): ${[...unresolved].slice(0, 12).join(', ')}${unresolved.size > 12 ? ' …' : ''}.`);
  for (const k in skipped) if (k !== 'zeroLengthMembers') warnings.push(`${skipped[k]} × ${k} not read.`);
  const model = { kind: 'frame', name, nodes: weld.finish(), members };
  if (nodeNames.some(Boolean)) model.nodeNames = Array.from({ length: model.nodes.length / 3 }, (_, i) => nodeNames[i] || '');
  const r = finish([model], warnings, { format: 'cis2', units, crs: null, elementTypes: hist, skipped, schema: p.header.schema, originatingSystem: p.header.originatingSystem, instances: p.count });
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}
