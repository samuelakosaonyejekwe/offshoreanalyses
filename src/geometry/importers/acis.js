// ACIS SAT (text) — boundary representation with analytic and B-spline geometry. The binary form (SAB) is
// decoded by sab.js into the same records and read by the code below.
//
// File structure (as published for the SAT save format):
//   line 1   version (major × 100 + minor), record count, body count, history flag
//   line 2   product, ACIS version and date as counted strings          (version 2.0 and later)
//   line 3   millimetres per model unit, resabs, resnor                 (version 2.0 and later)
//   records  "<type> <fields …> #", optionally prefixed by a sequence number "-n". Pointers are "$n" (record
//            index, −1 = null). Up to 6.x a record starts with its attribute pointer; from 7.0 on that is
//            followed by an integer and a history pointer ("body $-1 -1 $-1 …"; two integers from R27 on) and
//            strings are "@len text". The data ends at "End-of-ACIS-data" or where a history section begins.
//
// Read: body → lump → shell → face → loop → coedge → edge → vertex → point, with the body transform; tolerant
// topology (tcoedge, tedge, tvertex) is read like the exact one.
//   surfaces  plane-surface, cone-surface (cylinders and cones, elliptical ratio), sphere-surface, torus-surface,
//             spline-surface when the record stores a B-spline (exactsur, or the approximation kept with a
//             procedural definition: rbblnsur, sweepsur, sumsur, skinsur, offsur, rot_spl_sur …), trimmed by the
//             face boundary
//   curves    straight-curve, ellipse-curve (circles, ellipses and their arcs), intcurve-curve when the record
//             stores its B-spline (exactcur, or the stored approximation of a procedural curve)
// Sub-type blocks "{ ref n }" refer to the n-th block of the file (counted where each block opens).
// Each face is triangulated in the parameter space of its surface between its boundary loops (brep.js), the
// boundary points being those of the shared edges, so neighbouring faces meet without gaps.
//
// Conventions relied on (stated here because they decide the geometry; checked on files written by Autodesk
// Inventor, AutoCAD-family writers, Revit, Fusion, SESAM and ACIS 1.x – 28 — see docs/geometry.md):
//   * cone-surface: radius at axial distance z from the base ellipse is R + z·sin/cos; the surface normal points
//     away from the axis when cos > 0 and towards it when cos < 0. sphere / torus: a negative (minor) radius
//     turns the normal inwards.
//   * a face whose sense is "reversed" has its normal opposite to the surface normal; loops run counter-clockwise
//     about the face normal (material on the left). Before 2.0 coedge and edge senses are written as 0 / 1.
//   * an edge stores the curve parameters of its end vertices (5.0 and later); a reversed edge runs through the
//     curve at the negated parameters, and so does a reversed intcurve through its B-spline.
//   * transform: row-vector convention, p′ = scale·(p·A) + t.
//   * B-spline knot lists hold (value, multiplicity) pairs with one end knot fewer than the usual clamped
//     vector at each end; rational control points are stored as x y z w, not premultiplied.
//
// Not read: spline-surface faces without a stored B-spline, or whose boundary cannot be placed on it — counted
// and reported per reason; spindle tori; pcurves, attributes, subshells, history. Record types that are not
// recognised are reported, never interpreted.

import { inputText, baseName, finish, surfaceModel, Welder, V3, tally, cleanName, Warnings } from './common.js';
import { resolveUnits } from '../units.js';
import { sampleCurve } from './nurbs.js';
import { stepper, planeSurface, coneSurface, sphereSurface, torusSurface, splineSurface, tessellateFace, BodyMesh, CROSSING } from './brep.js';

const TAU = 2 * Math.PI, D2R = Math.PI / 180;
const NUM = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const UNIT_BY_MM = [[1, 'mm'], [10, 'cm'], [100, 'dm'], [1000, 'm'], [1e6, 'km'], [25.4, 'in'], [304.8, 'ft'], [914.4, 'yd'], [0.001, 'um']];
// record kinds that carry no geometry of their own and are passed over without comment
const QUIET = new Set(['attrib', 'pcurve', 'asmheader', 'subshell', 'wire', 'lump', 'shell', 'body', 'loop', 'coedge', 'edge', 'vertex', 'point', 'face', 'transform', 'curve', 'surface', 'eye_refinement', 'vertex_template', 'history_stream', 'delta_state', 'annotation', 'pattern', 'law_data', 'entity']);

/** Record names that end the entity data: the end marker, or the start of a history section (not read). */
export const END_MARK = /^(End-of-(ACIS|ASM)-data|Begin-of-(ACIS|ASM)-History-Data)$/i;

const isWs = (c) => c === 32 || (c >= 9 && c <= 13);
const lerp3 = (a, b, f) => [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
const dist3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

// ---------------------------------------------------------------------------------------------------------
// Text → header and records
// ---------------------------------------------------------------------------------------------------------

/**
 * Split a SAT text into its header and records.
 * @returns {{ header: object, records: Array<{ type: string, kind: string, tok: object[], idx: number }>, numbered: boolean }}
 */
export function readSat(text, fname = 'SAT') {
  const n = text.length;
  let i = 0;
  const next = () => {
    while (i < n && isWs(text.charCodeAt(i))) i++;
    if (i >= n) return null;
    if (text[i] === '@') {
      let j = i + 1;
      while (j < n && text.charCodeAt(j) >= 48 && text.charCodeAt(j) <= 57) j++;
      if (j > i + 1) { const len = Number(text.slice(i + 1, j)); const s = text.substr(j + 1, len); i = j + 1 + len; return { t: 's', v: s }; }
    }
    let j = i;
    while (j < n && !isWs(text.charCodeAt(j))) j++;
    if (j - i > 1 && text[j - 1] === '#') j--;
    const w = text.slice(i, j);
    i = j;
    if (w === '#' || w === '{' || w === '}') return { t: w };
    if (w[0] === '$' && /^\$-?\d+$/.test(w)) return { t: 'p', v: Number(w.slice(1)) };
    if (NUM.test(w)) return { t: 'n', v: Number(w), raw: w };
    return { t: 'w', v: w };
  };
  const notSat = () => new Error(`${fname}: not an ACIS SAT text file (expected a first line of four integers: version, record count, body count, flag).`);
  const h = [next(), next(), next(), next()];
  if (h.some((t) => !t || t.t !== 'n' || !Number.isInteger(t.v))) throw notSat();
  const header = { version: h[0].v, records: h[1].v, bodies: h[2].v, flag: h[3].v, product: null, acis: null, date: null, mmPerUnit: null, resabs: null, resnor: null };
  if (header.version < 100 || header.version > 99999) throw notSat();
  if (header.version >= 200) {
    // three counted strings ("@len text" from 7.0, "len text" before), then three numbers
    const save = i;
    const strs = [];
    for (let k = 0; k < 3; k++) {
      const t = next();
      if (t && t.t === 's') strs.push(t.v);
      else if (t && t.t === 'n' && /^\d+$/.test(t.raw)) { strs.push(text.substr(i + 1, t.v)); i += 1 + t.v; } else break;
    }
    const nums = strs.length === 3 ? [next(), next(), next()] : [];
    if (nums.length === 3 && nums.every((t) => t && t.t === 'n')) {
      [header.product, header.acis, header.date] = strs.map((s) => s.trim());
      header.mmPerUnit = nums[0].v; header.resabs = nums[1].v; header.resnor = nums[2].v;
    } else i = save;
    // from R27 on a fourth header line may follow: "T @len <stream identifier>" or "F"
    const mark = i, flag = next();
    if (flag && flag.t === 'w' && flag.v === 'T') { const id = next(); if (!id || id.t !== 's') i = mark; } else if (!(flag && flag.t === 'w' && flag.v === 'F' && (() => { const k = i, q = next(); i = k; return q && q.t === 'w'; })())) i = mark;
  }
  const records = [];
  let numbered = false;
  for (;;) {
    let t = next();
    if (!t) break;
    let idx = null;
    if (t.t === 'n') { idx = Math.abs(t.v); numbered = true; t = next(); }
    if (!t || t.t !== 'w') {
      if (!records.length) throw new Error(`${fname}: not an ACIS SAT text file (no entity record follows the header).`);
      throw new Error(`${fname}: record ${records.length} does not start with an entity name${t ? '' : ' — the file is truncated'}.`);
    }
    const type = t.v;
    if (END_MARK.test(type)) break;
    const tok = [];
    for (;;) {
      const q = next();
      if (!q) throw new Error(`${fname}: the SAT file is truncated: record ${records.length} (${type}) is not closed by "#".`);
      if (q.t === '#') break;
      tok.push(q);
    }
    records.push({ type, kind: type.slice(type.lastIndexOf('-') + 1), tok, idx: idx === null ? records.length : idx });
  }
  if (!records.length) throw new Error(`${fname}: the SAT file holds a header but no entity records.`);
  return { header, records, numbered };
}

// ---------------------------------------------------------------------------------------------------------
// B-spline data inside "{ … }" sub-type blocks
// ---------------------------------------------------------------------------------------------------------

const isN = (t) => t && t.t === 'n';
const isW = (t, ...v) => t && t.t === 'w' && v.includes(t.v);
// a keyword of the text form, or the enumeration value the binary form stores in its place
const isTag = (t) => t && (t.t === 'w' || t.t === 'e');

function knotList(tok, j, count) {
  const vals = [];
  let sum = 0;
  for (let k = 0; k < count; k++, j += 2) {
    if (!isN(tok[j]) || !isN(tok[j + 1]) || !(tok[j + 1].v >= 1) || !Number.isInteger(tok[j + 1].v)) return null;
    if (vals.length && tok[j].v < vals[vals.length - 1]) return null;
    for (let m = 0; m < tok[j + 1].v; m++) vals.push(tok[j].v);
    sum += tok[j + 1].v;
  }
  if (!vals.length) return null;
  return { U: [vals[0], ...vals, vals[vals.length - 1]], sum, j };
}

/** B-spline curve stored at tok[i] = nubs | nurbs. */
function bs3Curve(tok, i) {
  const rational = tok[i].v === 'nurbs';
  if (!isN(tok[i + 1]) || !isTag(tok[i + 2]) || !isN(tok[i + 3])) return null;
  const deg = tok[i + 1].v;
  if (!(deg >= 1 && deg <= 25)) return null;
  const k = knotList(tok, i + 4, tok[i + 3].v);
  if (!k) return null;
  const ncp = k.sum - deg + 1, dim = rational ? 4 : 3;
  if (ncp < deg + 1) return null;
  const P = new Float64Array(ncp * 3), W = rational ? new Float64Array(ncp) : null;
  let j = k.j;
  for (let c = 0; c < ncp; c++, j += dim) {
    for (let d = 0; d < dim; d++) if (!isN(tok[j + d])) return null;
    P[3 * c] = tok[j].v; P[3 * c + 1] = tok[j + 1].v; P[3 * c + 2] = tok[j + 2].v;
    if (rational) { if (!(tok[j + 3].v > 0)) return null; W[c] = tok[j + 3].v; }
  }
  return { deg, U: k.U, P, W, end: j };
}

/** B-spline surface stored at tok[i] = nubs | nurbs; control points kept in file order with both counts. */
function bs3Surface(tok, i) {
  const rational = tok[i].v === 'nurbs';
  if (!isN(tok[i + 1]) || !isN(tok[i + 2])) return null;
  const pu = tok[i + 1].v, pv = tok[i + 2].v;
  let j = i + 3;
  // newer files name the rational direction ("u", "v", "both") before the closure and singularity keywords
  if (isW(tok[j], 'u', 'v', 'both') && isTag(tok[j + 4])) j++;
  const words = [];
  while (isTag(tok[j]) && words.length < 4) words.push(tok[j++].v);
  if (words.length !== 4) return null;
  if (!isN(tok[j]) || !isN(tok[j + 1]) || !(pu >= 1 && pv >= 1 && pu <= 25 && pv <= 25)) return null;
  const ku = knotList(tok, j + 2, tok[j].v);
  if (!ku) return null;
  const kv = knotList(tok, ku.j, tok[j + 1].v);
  if (!kv) return null;
  const nu = ku.sum - pu + 1, nv = kv.sum - pv + 1, dim = rational ? 4 : 3;
  if (nu < pu + 1 || nv < pv + 1) return null;
  const P = new Float64Array(nu * nv * 3), W = rational ? new Float64Array(nu * nv) : null;
  j = kv.j;
  for (let c = 0; c < nu * nv; c++, j += dim) {
    for (let d = 0; d < dim; d++) if (!isN(tok[j + d])) return null;
    P[3 * c] = tok[j].v; P[3 * c + 1] = tok[j + 1].v; P[3 * c + 2] = tok[j + 2].v;
    if (rational) { if (!(tok[j + 3].v > 0)) return null; W[c] = tok[j + 3].v; }
  }
  return { pu, pv, U: ku.U, V: kv.U, nu, nv, P, W };
}

// ---------------------------------------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------------------------------------

const ANALYTIC = /^(plane|cone|sphere|torus)-surface$/;
// names the binary form (and the text form of recent releases) uses for the sub-types, mapped to the classic ones
const SUBTYPE_ALIAS = {
  exact_int_cur: 'exactcur', surf_int_cur: 'surfintcur', int_int_cur: 'intcur', par_int_cur: 'parcur', off_int_cur: 'offintcur', law_int_cur: 'lawintcur', proj_int_cur: 'projcur', blend_int_cur: 'bldcur', spring_int_cur: 'springcur', offset_int_cur: 'offsetintcur', sil_int_cur: 'silcur', helix_int_cur: 'helixcur',
  exact_spl_sur: 'exactsur', sum_spl_sur: 'sumsur', cyl_spl_sur: 'cylsur', rot_spl_sur: 'rotsur', off_spl_sur: 'offsur', sweep_spl_sur: 'sweepsur', skin_spl_sur: 'skinsur', loft_spl_sur: 'loftsur', net_spl_sur: 'netsur', rb_blend_spl_sur: 'rbblnsur', law_spl_sur: 'lawsur', pipe_spl_sur: 'pipesur', ruled_spl_sur: 'ruledsur', tube_spl_sur: 'tubesur', vertex_blend_spl_sur: 'vertexblendsur', g2_blend_spl_sur: 'g2blnsur', sub_spl_sur: 'subsur', var_blend_spl_sur: 'varblendsplsur', exp_par_cur: 'exppc',
};

/**
 * Tessellate the bodies of a decoded ACIS file. Coordinates stay in model units; the body transforms are applied.
 * @param {{ header: object, records: object[] }} sat result of readSat (or of sab.js readSab)
 * @param {{ deflection?: number, angularStep?: number, name?: string }} [opts] `deflection` in model units
 * @returns {{ bodies: Array<{ name, coords: Float64Array, indices: Uint32Array, groups: object[], closed: boolean, oriented: boolean, volume: number|null }>,
 *   wires: Array<{ a: number[], b: number[], group: string }>, warnings: string[], hist: object, skipped: object, stats: object }}
 */
export function acisBodies(sat, opts = {}) {
  const { header, records } = sat;
  const name = opts.name || 'body';
  const W = new Warnings();
  const extra = [];
  const skipped = {}, hist = {}, faceSkips = {};
  const modern = header.version >= 700;

  // ---- record access ----
  const byIdx = new Map();
  for (const r of records) if (!byIdx.has(r.idx)) byIdx.set(r.idx, r);
  if (header.records > 0 && header.records !== records.length) extra.push(`The header announces ${header.records} records but ${records.length} were found.`);
  const deref = (k) => {
    if (!(k >= 0)) return null;
    const r = byIdx.get(k);
    if (!r) W.add('A record points to an entity that is not in the file; the reference was ignored');
    return r || null;
  };
  const dataOf = (r) => {
    if (r.d) return r.d;
    const t = r.tok;
    let k = t[0] && t[0].t === 'p' ? 1 : 0;
    if (modern && k) {
      let j = k;
      while (j < k + 2 && isN(t[j]) && Number.isInteger(t[j].v)) j++;
      if (j > k && t[j] && t[j].t === 'p') k = j + 1;
    }
    return (r.d = t.slice(k));
  };
  const firstOf = (r, kind) => { for (const t of dataOf(r)) if (t.t === 'p') { const q = t.v >= 0 ? byIdx.get(t.v) : null; if (q && q.kind === kind) return q; } return null; };
  /** Sense of a face, coedge, edge or curve: true = reversed. Before 2.0 the sense is written as 0 / 1 after the last pointer of the record's fixed part. */
  const reversedOf = (r) => {
    const d = dataOf(r);
    let np = 0;
    for (let i = 0; i < d.length; i++) {
      const t = d[i];
      if (t.t === '{') break;
      if (t.t === 'b') return t.v;
      if (t.t === 'w' && (t.v === 'forward' || t.v === 'reversed')) return t.v === 'reversed';
      if (t.t === 'p') np++;
      else if (t.t === 'n' && d[i - 1].t === 'p' && (t.v === 0 || t.v === 1) && (((r.kind === 'coedge' || r.kind === 'edge') && np === 4) || (r.kind === 'face' && np === 5))) return t.v === 1;
    }
    return false;
  };
  // numbers outside sub-type blocks; the bounds of finite intervals ("F <value>", binary: true + value) are not data
  const numsOf = (r) => {
    const o = [], d = dataOf(r);
    let depth = 0;
    for (let i = 0; i < d.length; i++) {
      if (d[i].t === '{') depth++; else if (d[i].t === '}') depth--;
      else if (!depth && d[i].t === 'n' && !(i && (isW(d[i - 1], 'F') || (d[i - 1].t === 'b' && d[i - 1].v)))) o.push(d[i].v);
    }
    return o;
  };

  // ---- sub-type blocks, numbered where they open, for "{ ref n }" ----
  let blocks = null;
  const blockTable = () => {
    if (blocks) return blocks;
    blocks = [];
    for (const r of records) for (let i = 0; i < r.tok.length; i++) if (r.tok[i].t === '{' && !(isW(r.tok[i + 1], 'ref') && isN(r.tok[i + 2]))) blocks.push({ tok: r.tok, at: i });
    return blocks;
  };
  /** Follow "{ ref n }" at tok[at] to the block that holds the data. */
  const resolveBlock = (tok, at) => {
    for (let guard = 0; guard < 8; guard++) {
      if (!(isW(tok[at + 1], 'ref') && isN(tok[at + 2]))) return { tok, at };
      const b = blockTable()[tok[at + 2].v];
      if (!b) return null;
      tok = b.tok; at = b.at;
    }
    return null;
  };
  /** Name of a block and, for curves, the position of the B-spline that follows the name. */
  const blockHead = (b) => {
    const nm = b.tok[b.at + 1];
    const name = nm && (nm.t === 'w' || nm.t === 's') ? nm.v : '';
    let j = b.at + 2;
    // from R21 on the name is followed by the release the definition was saved for
    if (isN(b.tok[j]) && b.tok[j + 1] && (isW(b.tok[j + 1], 'full', 'none', 'summary', 'nubs', 'nurbs', 'nullbs') || b.tok[j + 1].t === 'e')) j++;
    if (isW(b.tok[j], 'full') || (b.tok[j] && b.tok[j].t === 'e')) j++;
    return { name: SUBTYPE_ALIAS[name] || name, at: isW(b.tok[j], 'nubs', 'nurbs') ? j : -1 };
  };
  /** B-spline surface stored at the top level of a spline-surface block: "… full nubs|nurbs pu pv …". */
  const storedSurface = (b) => {
    let depth = 0;
    for (let i = b.at; i < b.tok.length; i++) {
      const t = b.tok[i];
      if (t.t === '{') depth++;
      else if (t.t === '}') { if (--depth === 0) break; } else if (depth === 1 && isW(t, 'nubs', 'nurbs') && isN(b.tok[i + 1]) && isN(b.tok[i + 2]) && isTag(b.tok[i + 3])) { const s = bs3Surface(b.tok, i); if (s) return s; }
    }
    return null;
  };

  // ---- sum surfaces kept as their definition: S(u, v) = c1(u) + c2(v) − P ----
  const isBool = (t) => t && (t.t === 'b' || isW(t, 'F', 'I'));
  const finite = (t) => (t.t === 'b' ? t.v : t.v === 'F');
  /** Interval end at tok[j]: "I", or "F value" (binary: false, or true and the value). Returns [value|null, next]. */
  const bound = (tok, j) => (isBool(tok[j]) ? (finite(tok[j]) ? (isN(tok[j + 1]) ? [tok[j + 1].v, j + 2] : null) : [null, j + 1]) : null);
  /** Curve written inside a surface definition: "intcurve sense { … } range" or "straight P D range". */
  const inlineCurve = (tok, j) => {
    if (isW(tok[j], 'straight')) {
      const n = tok.slice(j + 1, j + 7);
      if (n.length < 6 || !n.every(isN)) return null;
      const a = bound(tok, j + 7), b = a && bound(tok, a[1]);
      return b ? { line: { O: n.slice(0, 3).map((t) => t.v), D: n.slice(3, 6).map((t) => t.v) }, range: [a[0], b[0]], next: b[1] } : null;
    }
    if (isW(tok[j], 'intcurve') && tok[j + 2] && tok[j + 2].t === '{') {
      const rev = tok[j + 1].t === 'b' ? tok[j + 1].v : isW(tok[j + 1], 'reversed');
      const blk = resolveBlock(tok, j + 2);
      const head = blk && blockHead(blk);
      const bs = head && head.at >= 0 ? bs3Curve(blk.tok, head.at) : null;
      let k = j + 2, depth = 0;
      for (; k < tok.length; k++) { if (tok[k].t === '{') depth++; else if (tok[k].t === '}' && --depth === 0) break; }
      const a = bound(tok, k + 1), b = a && bound(tok, a[1]);
      if (!bs || !b) return null;
      if (rev) {
        // the same curve run backwards: parameter t ↦ −t
        const m = bs.P.length / 3, P = new Float64Array(bs.P.length), Wt = bs.W ? new Float64Array(m) : null;
        for (let i = 0; i < m; i++) { const q = m - 1 - i; P[3 * i] = bs.P[3 * q]; P[3 * i + 1] = bs.P[3 * q + 1]; P[3 * i + 2] = bs.P[3 * q + 2]; if (Wt) Wt[i] = bs.W[q]; }
        return { bs: { deg: bs.deg, U: bs.U.map((u) => -u).reverse(), P, W: Wt }, range: [a[0], b[0]], next: b[1] };
      }
      return { bs, range: [a[0], b[0]], next: b[1] };
    }
    return null;
  };
  /** Control net of the sum surface of a block "{ sumsur [release] curve1 curve2 P … }", or null. */
  const sumSurface = (b) => {
    let j = b.at + 2;
    if (isN(b.tok[j])) j++;
    const c1 = inlineCurve(b.tok, j), c2 = c1 && inlineCurve(b.tok, c1.next);
    if (!c2) return null;
    const pn = b.tok.slice(c2.next, c2.next + 3);
    if (pn.length < 3 || !pn.every(isN)) return null;
    const P0 = pn.map((t) => t.v);
    // the parameter ranges of the surface follow the approximation keyword; they bound the straight lines
    let k = c2.next + 3;
    if (isTag(b.tok[k])) k++;
    const ranges = [];
    for (let q = 0; q < 4; q++) { const r = bound(b.tok, k); if (!r) break; ranges.push(r[0]); k = r[1]; }
    const asSpline = (c, lo, hi) => {
      if (c.bs) return c.bs;
      const t0 = c.range[0] ?? lo, t1 = c.range[1] ?? hi;
      if (!(Number.isFinite(t0) && Number.isFinite(t1) && t1 > t0)) return null;
      const { O, D } = c.line;
      return { deg: 1, U: [t0, t0, t1, t1], P: Float64Array.from([O[0] + t0 * D[0], O[1] + t0 * D[1], O[2] + t0 * D[2], O[0] + t1 * D[0], O[1] + t1 * D[1], O[2] + t1 * D[2]]), W: null };
    };
    const A = asSpline(c1, ranges[0], ranges[1]), B = asSpline(c2, ranges[2], ranges[3]);
    if (!A || !B) return null;
    const nu = A.P.length / 3, nv = B.P.length / 3;
    const P = new Float64Array(nu * nv * 3), Wt = A.W || B.W ? new Float64Array(nu * nv) : null;
    for (let iv = 0; iv < nv; iv++) for (let iu = 0; iu < nu; iu++) {
      const o = iu + nu * iv;
      for (let d = 0; d < 3; d++) P[3 * o + d] = A.P[3 * iu + d] + B.P[3 * iv + d] - P0[d];
      if (Wt) Wt[o] = (A.W ? A.W[iu] : 1) * (B.W ? B.W[iv] : 1);
    }
    return { pu: A.deg, pv: B.deg, U: A.U, V: B.U, nu, nv, P, W: Wt, exact: true };
  };

  for (const r of records) if (!QUIET.has(r.kind)) tally(skipped, `record type "${r.type}"`);

  const pointOf = (v) => {
    if (!v) return null;
    if (v.pos !== undefined) return v.pos;
    const pt = firstOf(v, 'point');
    const n = pt ? numsOf(pt) : [];
    return (v.pos = n.length >= 3 && n.slice(-3).every(Number.isFinite) ? n.slice(-3) : null);
  };

  const transformOf = (body) => {
    const tr = body ? firstOf(body, 'transform') : null;
    if (!tr) return null;
    const n = numsOf(tr);
    if (n.length < 13 || !n.slice(-13).every(Number.isFinite)) { W.add('A transform record could not be read; its body is placed without it'); return null; }
    const a = n.slice(-13);
    const A = a.slice(0, 9), t = a.slice(9, 12), s = a[12];
    const det = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
    return { apply: (p) => [s * (p[0] * A[0] + p[1] * A[3] + p[2] * A[6]) + t[0], s * (p[0] * A[1] + p[1] * A[4] + p[2] * A[7]) + t[1], s * (p[0] * A[2] + p[1] * A[5] + p[2] * A[8]) + t[2]], mirror: det * s < 0 };
  };

  // everything a body owns, by following pointers (robust against differences in field order between versions)
  const owned = (root) => {
    const seen = new Set([root]);
    const stack = [root];
    while (stack.length) {
      const r = stack.pop();
      for (const t of r.tok) {
        if (t.t !== 'p' || t.v < 0) continue;
        const q = byIdx.get(t.v);
        if (!q || seen.has(q) || q.kind === 'body' || q.kind === 'attrib' || q.kind === 'transform') continue;
        seen.add(q);
        stack.push(q);
      }
    }
    return [...seen].sort((a, b) => a.idx - b.idx);
  };

  const bodyRecs = records.filter((r) => r.kind === 'body');
  const roots = bodyRecs.length ? bodyRecs.map((b) => ({ body: b, recs: owned(b) })) : [{ body: null, recs: records }];
  const bodies = [], wires = [];
  const anyFace = records.some((r) => r.kind === 'face');
  let facesRead = 0, faceTotal = 0, wireEdges = 0, chord = 0, approxCurves = 0;
  const angStep = (opts.angularStep ?? 15) * D2R;
  const skipFace = (why) => { tally(faceSkips, why); };

  roots.forEach((root, bi) => {
    const bodyName = roots.length > 1 ? `${name} body ${bi + 1}` : name;
    const tf = transformOf(root.body);
    const faces = root.recs.filter((r) => r.kind === 'face');
    faceTotal += faces.length;
    // extent of the body for the default deflection: its vertices and the spheres, tori and full circles it holds
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    const see = (p, r = 0) => { for (let k = 0; k < 3; k++) { if (p[k] - r < lo[k]) lo[k] = p[k] - r; if (p[k] + r > hi[k]) hi[k] = p[k] + r; } };
    for (const r of root.recs) {
      if (r.kind === 'vertex') { const p = pointOf(r); if (p) see(p); } else if (r.type === 'sphere-surface') { const n = numsOf(r); if (n.length >= 4) see(n, Math.abs(n[3])); } else if (r.type === 'torus-surface') { const n = numsOf(r); if (n.length >= 8) see(n, Math.abs(n[6]) + Math.abs(n[7])); } else if (r.kind === 'edge') {
        // a full circle or ellipse reaches beyond its single vertex; an arc is covered well enough by its ends
        const ps = dataOf(r).filter((t) => t.t === 'p');
        const cv = ps.length >= 2 && ps[0].v === ps[1].v ? firstOf(r, 'curve') : null;
        if (cv && cv.type === 'ellipse-curve') { const n = numsOf(cv); if (n.length >= 9) see(n, Math.hypot(n[6], n[7], n[8])); }
      }
    }
    const diag = lo[0] === Infinity ? 0 : Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    if (!(diag > 0)) {
      // nothing with an extent: only faces on surfaces that are not read, or records without coordinates
      for (const r of faces) {
        const sr = firstOf(r, 'surface');
        if (sr && !ANALYTIC.test(sr.type)) skipFace(sr.type === 'spline-surface' ? 'spline surface on a body without any vertex' : `surface type "${sr.type}"`); else skipFace('face without any vertex or extent');
      }
      return;
    }
    let defl = opts.deflection !== undefined ? opts.deflection : 0.01 * diag;
    // One pass over the edges and faces of the body at a given facet size. Where the facets of two boundary
    // loops of a face cross (a hole close to a coarsely faceted outline), the body is done again more finely.
    const attempt = (defl) => {
      const loc = { W: new Warnings(), hist: {}, skips: {}, read: 0, approx: 0, crossing: 0, skip(why) { tally(this.skips, why); } };
      for (const r of root.recs) { delete r.pts; delete r.surf; delete r.spl; }
      const stepFor = stepper(defl, angStep);
      const tol = diag * 1e-9;
      // an end point computed from the stored curve parameters must meet its vertex this closely to be believed
      const meet = Math.max(1e-4 * diag, 100 * (header.resabs > 0 ? header.resabs : 0));

      // ---- edges: one polyline per edge, shared by every face that uses it ----
      const edgePts = (e) => {
        if (e.pts) return e.pts;
        const d = dataOf(e);
        const ends = [];
        for (let i = 0; i < d.length && ends.length < 2; i++) if (d[i].t === 'p') { const q = d[i].v >= 0 ? byIdx.get(d[i].v) : null; ends.push({ v: q && q.kind === 'vertex' ? q : null, t: isN(d[i + 1]) ? d[i + 1].v : null }); }
        const Sp = pointOf(ends[0] && ends[0].v), Ep = pointOf(ends[1] && ends[1].v ? ends[1].v : ends[0] && ends[0].v);
        const cv = firstOf(e, 'curve');
        const rev = reversedOf(e);
        const chordPts = (why) => { if (why) loc.W.add(why); return Sp && Ep && Sp !== Ep ? [Sp, Ep] : Sp ? [Sp] : []; };
        if (!cv) return (e.pts = Sp ? [Sp] : []);
        const same = !Sp || !Ep || Sp === Ep || dist3(Sp, Ep) <= tol;
        // curve parameters of the two ends, when the record gives them
        let c0 = null, c1 = null;
        if (ends.length === 2 && ends[0].t !== null && ends[1].t !== null && ends[1].t > ends[0].t) { c0 = rev ? -ends[0].t : ends[0].t; c1 = rev ? -ends[1].t : ends[1].t; }
        if (cv.type === 'straight-curve') return (e.pts = chordPts());
        if (cv.type === 'ellipse-curve') {
          const n = numsOf(cv);
          if (n.length < 9 || !n.slice(0, 9).every(Number.isFinite)) return (e.pts = chordPts('An ellipse-curve record could not be read; its edge was replaced by a straight chord'));
          const C = n.slice(0, 3), N = V3.norm(n.slice(3, 6)), M = n.slice(6, 9), ratio = n.length >= 10 && n[9] > 0 ? n[9] : 1;
          const Y = V3.scale(V3.cross(N, M), ratio), m2 = V3.dot(M, M), y2 = V3.dot(Y, Y);
          if (!(m2 > 0) || !(y2 > 0)) return (e.pts = chordPts('An ellipse-curve record has no radius; its edge was replaced by a straight chord'));
          const at = (t) => { const c = Math.cos(t), s = Math.sin(t); return [C[0] + c * M[0] + s * Y[0], C[1] + c * M[1] + s * Y[1], C[2] + c * M[2] + s * Y[2]]; };
          const ang = (p) => { const dd = V3.sub(p, C); return Math.atan2(V3.dot(dd, Y) / y2, V3.dot(dd, M) / m2); };
          let a0, sweep;
          if (c0 !== null && Sp && Ep && Math.abs(c1 - c0) <= TAU * (1 + 1e-9) && dist3(at(c0), Sp) <= meet && dist3(at(c1), Ep) <= meet) { a0 = c0; sweep = c1 - c0; } else {
            const dir = rev ? -1 : 1;
            a0 = Sp ? ang(Sp) : 0;
            sweep = dir * TAU;
            if (!same) { sweep = ang(Ep) - a0; if (dir > 0) { while (sweep <= 1e-12) sweep += TAU; } else { while (sweep >= -1e-12) sweep -= TAU; } }
          }
          const closedArc = same || !Ep;
          const k = Math.max(closedArc ? 3 : 1, Math.ceil(Math.abs(sweep) / stepFor(Math.sqrt(Math.max(m2, y2))) - 1e-9));
          const pts = [];
          for (let i = 0; i <= k; i++) {
            if (i === 0 && Sp) pts.push(Sp);
            else if (i === k && (closedArc ? Sp : Ep)) pts.push(closedArc ? Sp : Ep);
            else pts.push(at(a0 + (sweep * i) / k));
          }
          if (closedArc && !Sp) pts[k] = pts[0];
          return (e.pts = pts);
        }
        if (cv.type === 'intcurve-curve') {
          const at0 = cv.tok.findIndex((t) => t.t === '{');
          const blk = at0 >= 0 ? resolveBlock(cv.tok, at0) : null;
          const head = blk ? blockHead(blk) : null;
          const bs = head && head.at >= 0 ? bs3Curve(blk.tok, head.at) : null;
          if (!bs) return (e.pts = chordPts(`An intcurve-curve without stored B-spline data${head && head.name ? ` ("${head.name}")` : ''} was replaced by a straight chord between its end points`));
          const curveRev = reversedOf(cv);
          const kLo = bs.U[bs.deg], kHi = bs.U[bs.U.length - 1 - bs.deg];
          const toPts = (flat) => { const o = []; for (let i = 0; i < flat.length; i += 3) o.push([flat[i], flat[i + 1], flat[i + 2]]); return o; };
          const ends3 = sampleCurve(bs.deg, bs.U, bs.P, bs.W, 1, kLo, kHi);
          if (!ends3) return (e.pts = chordPts('An intcurve-curve whose B-spline data are inconsistent was replaced by a straight chord'));
          const closedCurve = Math.hypot(ends3[0] - ends3[ends3.length - 3], ends3[1] - ends3[ends3.length - 2], ends3[2] - ends3[ends3.length - 1]) <= meet;
          let pts = null;
          if (c0 !== null && Sp && Ep) {
            // parameters on the B-spline, increasing; `back` when the edge runs against it
            let b0 = curveRev ? -c0 : c0, b1 = curveRev ? -c1 : c1;
            const back = b1 < b0;
            if (back) { const x = b0; b0 = b1; b1 = x; }
            const T = kHi - kLo, slack = 1e-6 * T;
            let piece = null;
            if (closedCurve && b1 - b0 <= T + slack) {
              const sh = T * Math.floor((b0 - kLo + slack) / T);
              b0 -= sh; b1 -= sh;
              if (b0 < kLo) b0 = kLo;
              if (b1 <= kHi + slack) piece = toPts(sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8, b0, Math.min(b1, kHi)));
              else piece = toPts(sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8, b0, kHi)).concat(toPts(sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8, kLo, b1 - T)).slice(1));
            } else if (b0 >= kLo - slack && b1 <= kHi + slack) piece = toPts(sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8, Math.max(b0, kLo), Math.min(b1, kHi)));
            if (piece && piece.length >= 2) {
              if (back) piece.reverse();
              if (dist3(piece[0], Sp) <= meet && dist3(piece[piece.length - 1], Ep) <= meet) pts = piece;
            }
          }
          if (!pts) {
            // no usable parameters: cut the sampled curve where the vertices lie on it
            const sm = toPts(sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8));
            if (!Sp || !Ep) return (e.pts = sm);
            const offSeg = (p, a, b) => { const ee = V3.sub(b, a), l2 = V3.dot(ee, ee); const f = l2 > 0 ? Math.max(0, Math.min(1, V3.dot(V3.sub(p, a), ee) / l2)) : 0; return dist3(p, lerp3(a, b, f)); };
            const near = (p) => {
              let b = 0, dd = Infinity;
              sm.forEach((q, i) => { const x = dist3(p, q); if (x < dd) { dd = x; b = i; } });
              const sp = Math.max(b > 0 ? dist3(sm[b], sm[b - 1]) : 0, b + 1 < sm.length ? dist3(sm[b], sm[b + 1]) : 0);
              const off = Math.min(b > 0 ? offSeg(p, sm[b - 1], sm[b]) : Infinity, b + 1 < sm.length ? offSeg(p, sm[b], sm[b + 1]) : Infinity);
              return { i: b, ok: off <= 0.1 * sp + tol };
            };
            if (same) { if (closedCurve && dist3(Sp, sm[0]) <= tol * 1e3 + defl * 1e-3) pts = (rev !== curveRev) ? sm.slice().reverse() : sm.slice(); } else if (!closedCurve) {
              const a = near(Sp), b = near(Ep);
              if (a.ok && b.ok && a.i !== b.i) pts = a.i < b.i ? sm.slice(a.i, b.i + 1) : sm.slice(b.i, a.i + 1).reverse();
            }
          }
          if (!pts) return (e.pts = chordPts('An intcurve-curve whose B-spline does not pass through its edge vertices in an unambiguous way was replaced by a straight chord'));
          if (head.name !== 'exactcur') loc.approx++;
          pts[0] = Sp; pts[pts.length - 1] = same ? Sp : Ep;
          return (e.pts = pts);
        }
        return (e.pts = chordPts(`A curve of type "${cv.type}" is not read; its edge was replaced by a straight chord`));
      };

      // how many coedges of the faces use each edge: 2 on a closed manifold solid
      const usedEdges = new Map();
      const loopsOf = (face) => {
        const loops = [];
        const seenLoops = new Set();
        for (let lp = firstOf(face, 'loop'); lp && !seenLoops.has(lp); lp = firstOf(lp, 'loop')) {
          seenLoops.add(lp);
          const c0 = firstOf(lp, 'coedge');
          const pts = [];
          const seen = new Set();
          for (let c = c0; c && !seen.has(c);) {
            seen.add(c);
            const e = firstOf(c, 'edge');
            if (e) {
              usedEdges.set(e, (usedEdges.get(e) || 0) + 1);
              const ep = edgePts(e);
              const seq = reversedOf(c) ? ep.slice().reverse() : ep;
              for (let i = 0; i + 1 < seq.length; i++) if (!pts.length || pts[pts.length - 1] !== seq[i]) pts.push(seq[i]);
              if (seq.length === 1 && (!pts.length || pts[pts.length - 1] !== seq[0])) pts.push(seq[0]);
            }
            const nx = dataOf(c).find((t) => t.t === 'p');
            c = nx ? deref(nx.v) : null;
            if (c && c.kind !== 'coedge') c = null;
          }
          while (pts.length > 1 && pts[0] === pts[pts.length - 1]) pts.pop();
          if (pts.length) loops.push(pts);
        }
        return loops;
      };

      const surfaceOf = (sr) => {
        if (sr.surf !== undefined) return sr.surf;
        const n = numsOf(sr);
        const ok = (k) => n.length >= k && n.slice(0, k).every(Number.isFinite);
        let s = null;
        if (sr.type === 'plane-surface') {
          if (ok(6) && V3.len(n.slice(3, 6)) > 0) s = planeSurface(n.slice(0, 3), V3.norm(n.slice(3, 6)), n.length >= 9 ? n.slice(6, 9) : null);
        } else if (sr.type === 'cone-surface') {
          if (ok(12)) {
            const A = n.slice(3, 6), M = n.slice(6, 9), ratio = n[9], sin = n[10], cos = n[11], R0 = V3.len(M);
            if (V3.len(A) > 0 && R0 > 0 && ratio > 0 && Math.abs(cos) > 1e-12 && Math.abs(Math.hypot(sin, cos) - 1) <= 1e-3) s = coneSurface(n.slice(0, 3), V3.norm(A), M, R0, ratio, sin / cos, cos < 0, stepFor);
          }
        } else if (sr.type === 'sphere-surface') {
          if (ok(4) && Math.abs(n[3]) > 0) {
            const Z = n.length >= 10 && V3.len(n.slice(7, 10)) > 0 ? V3.norm(n.slice(7, 10)) : [0, 0, 1];
            s = sphereSurface(n.slice(0, 3), Math.abs(n[3]), Z, n.length >= 7 ? n.slice(4, 7) : null, n[3] < 0, stepFor);
          }
        } else if (sr.type === 'torus-surface') {
            if (ok(8) && V3.len(n.slice(3, 6)) > 0) s = torusSurface(n.slice(0, 3), V3.norm(n.slice(3, 6)), n[6], Math.abs(n[7]), n.length >= 11 ? n.slice(8, 11) : null, n[7] < 0, stepFor);
        }
        return (sr.surf = s);
      };
      /** Stored B-spline of a spline-surface record: { nets: both readings of the control-net order, name } or { why }. */
      const splineOf = (sr) => {
        if (sr.spl !== undefined) return sr.spl;
        const at0 = sr.tok.findIndex((t) => t.t === '{');
        const blk = at0 >= 0 ? resolveBlock(sr.tok, at0) : null;
        if (!blk) return (sr.spl = { why: 'spline surface whose definition could not be found' });
        const nm = blockHead(blk).name || 'unnamed';
        let bs = storedSurface(blk), how = nm === 'exactsur' ? 'spline' : `spline (stored approximation of ${nm})`;
        if (!bs && nm === 'sumsur') { bs = sumSurface(blk); how = 'spline (sum of two curves)'; }
        if (!bs) return (sr.spl = { why: `procedural spline surface without a stored B-spline ("${nm}")` });
        // the net is listed row by row with u running fastest (seen in every file checked); the other reading is
      // tried only when the boundary of the face does not lie on the first
        const vFast = new Float64Array(bs.P.length), vFastW = bs.W ? new Float64Array(bs.W.length) : null;
        for (let iu = 0; iu < bs.nu; iu++) for (let iv = 0; iv < bs.nv; iv++) {
          const from = iv + bs.nv * iu, to = iu + bs.nu * iv;
          vFast[3 * to] = bs.P[3 * from]; vFast[3 * to + 1] = bs.P[3 * from + 1]; vFast[3 * to + 2] = bs.P[3 * from + 2];
          if (vFastW) vFastW[to] = bs.W[from];
        }
        const mk = (P, Wt) => splineSurface({ pu: bs.pu, pv: bs.pv, U: bs.U, V: bs.V, nu: bs.nu, nv: bs.nv, P, W: Wt }, defl, angStep);
        const nets = (bs.exact ? [mk(bs.P, bs.W)] : [mk(bs.P, bs.W), mk(vFast, vFastW)]).filter(Boolean);
        return (sr.spl = nets.length ? { nets, name: nm, how } : { why: `spline surface whose stored B-spline is collapsed or inconsistent ("${nm}")` });
      };

      // ---- faces ----
      const mesh = new BodyMesh(diag);
      let faceNo = 0;
      for (const face of faces) {
        const sr = firstOf(face, 'surface');
        const flipFace = reversedOf(face);
        const label = `face ${++faceNo}`;
        if (!sr) { loc.skip('face without a surface record'); continue; }
        let loops;
        try { loops = loopsOf(face); } catch (e) { loc.skip(`boundary could not be followed (${e.message})`); continue; }

        if (sr.type === 'spline-surface') {
          const sp = splineOf(sr);
          if (sp.why) { loc.skip(sp.why); continue; }
          const srRev = reversedOf(sr);
          let done = false, why = '';
          for (const surf of sp.nets) {
            let res;
            try { res = tessellateFace(surf, loops, (surf.inward !== srRev) !== flipFace, Math.max(2 * defl, 1e-3 * diag)); } catch (e) { res = `triangulation failed (${e.message})`; }
            if (typeof res === 'string') { if (!why) why = res; continue; }
            if (mesh.addFace(res.pts, res.indices, (surf.inward !== srRev) !== flipFace, label)) { loc.read++; tally(loc.hist, sp.how); done = true; break; }
          }
          if (!done) { if (why === CROSSING) loc.crossing++; loc.skip(`spline surface: ${why || 'its boundary is degenerate'}`); }
          continue;
        }

        const surf = surfaceOf(sr);
        if (!surf) { loc.skip(ANALYTIC.test(sr.type) ? `${sr.type} record with missing or unsupported values` : `surface type "${sr.type}"`); continue; }
        const flip = surf.inward !== flipFace;
        let res;
        try { res = tessellateFace(surf, loops, flip, Math.max(2 * defl, 1e-3 * diag)); } catch (e) { res = `triangulation failed (${e.message})`; }
        if (typeof res === 'string') { if (res === CROSSING) loc.crossing++; loc.skip(`${surf.type}: ${res}`); continue; }
        if (!mesh.addFace(res.pts, res.indices, flip, label)) { loc.skip(`${surf.type}: its boundary is degenerate`); continue; }
        loc.read++;
        tally(loc.hist, surf.type);
      }

      return { loc, mesh, usedEdges, edgePts };
    };
    let pass = attempt(defl);
    for (let k = 0; k < 3 && pass.loc.crossing; k++) { defl /= 4; pass = attempt(defl); }
    if (defl !== (opts.deflection !== undefined ? opts.deflection : 0.01 * diag)) extra.push(`${bodyName}: tessellated with a chord error of ${defl.toPrecision(3)} model units instead of ${(opts.deflection !== undefined ? opts.deflection : 0.01 * diag).toPrecision(3)}, because the coarser facets of neighbouring boundary loops crossed.`);
    chord = Math.max(chord, defl);
    const { loc, mesh, usedEdges, edgePts } = pass;
    facesRead += loc.read;
    for (const k in loc.hist) tally(hist, k, loc.hist[k]);
    for (const k in loc.skips) tally(faceSkips, k, loc.skips[k]);
    if (!mesh.empty) {
      const m = mesh.finish(tf);
      if (m.reoriented) extra.push(`${bodyName}: the faces enclose a negative volume as stored, so their orientation was reversed to point outwards.`);
      let free = 0, shared = 0;
      for (const [e, n] of usedEdges) { if (firstOf(e, 'curve') === null && edgePts(e).length < 2) continue; if (n === 1) free++; else if (n > 2) shared++; }
      bodies.push({ name: cleanName(bodyName), coords: m.coords, indices: m.indices, groups: m.groups, closed: m.closed, oriented: m.oriented, volume: m.volume, faces: faces.length, facesRead: m.groups.length, freeEdges: free, nonManifoldEdges: shared });
    }

    // ---- wire edges (edges no face uses) ----
    const loose = root.recs.filter((r) => r.kind === 'edge' && !usedEdges.has(r));
    if (loose.length && anyFace) wireEdges += loose.length;
    else for (const e of loose) {
      const pts = edgePts(e);
      for (let i = 0; i + 1 < pts.length; i++) wires.push({ a: tf ? tf.apply(pts[i]) : pts[i], b: tf ? tf.apply(pts[i + 1]) : pts[i + 1], group: bodyName });
    }
    approxCurves += loc.approx;
    for (const [t, n] of loc.W.m) for (let q = 0; q < n; q++) W.add(t);
  });

  const warnings = [];
  const nSkip = Object.values(faceSkips).reduce((a, b) => a + b, 0);
  if (nSkip) warnings.push(`${nSkip} of ${faceTotal} face${faceTotal > 1 ? 's' : ''} could not be read (${Object.entries(faceSkips).map(([k, v]) => `${v} × ${k}`).join('; ')}) — export those bodies as STEP to get them complete.`);
  if (approxCurves) warnings.push(`${approxCurves} edge${approxCurves > 1 ? 's follow' : ' follows'} the B-spline approximation stored with a procedural curve.`);
  if (wireEdges) warnings.push(`${wireEdges} wire edge${wireEdges > 1 ? 's' : ''} outside any face ${wireEdges > 1 ? 'were' : 'was'} not imported (wire bodies are read only from files without faces).`);
  for (const k in skipped) warnings.push(`${skipped[k]} × ${k} not recognised and ignored.`);
  for (const w of W.list()) warnings.push(`${w}.`);
  warnings.push(...extra);
  for (const k in faceSkips) skipped[`face: ${k}`] = faceSkips[k];
  return {
    bodies, wires, warnings, hist, skipped,
    stats: { bodies: bodyRecs.length, faces: faceTotal, facesRead, facesSkipped: nSkip, closedBodies: bodies.filter((b) => b.closed).length, chord, angularStepDeg: opts.angularStep ?? 15 },
  };
}

/** Importer result from a decoded file; shared by the text and the binary reader. */
export function acisModels(sat, input, opts = {}, format = 'acis') {
  const fname = input.name || (format === 'acis' ? 'SAT' : 'SAB');
  const label = format === 'acis' ? 'SAT' : 'SAB';
  if (opts.deflection !== undefined && !(opts.deflection > 0)) throw new RangeError('deflection: must be a positive length in metres.');
  if (opts.angularStep !== undefined && !(opts.angularStep > 0 && opts.angularStep <= 90)) throw new RangeError('angularStep: must be within 0…90 degrees.');
  const { header } = sat;
  const name = baseName(input.name);
  const lead = [];

  let units;
  const mm = header.mmPerUnit;
  if (mm > 0) {
    const known = UNIT_BY_MM.find(([v]) => Math.abs(v - mm) <= 1e-9 * v);
    units = known || opts.units || (Number.isFinite(opts.scale) && opts.scale > 0) ? resolveUnits(opts, known ? known[1] : null, 'mm') : { unit: 'custom', scale: mm / 1000, basis: 'file' };
  } else {
    units = resolveUnits(opts, null, 'mm');
    if (units.basis === 'assumed') lead.push('The file declares no length unit; millimetres were assumed. Re-import with units set if the model was drawn in another unit.');
  }
  const S = units.scale;

  const res = acisBodies(sat, { deflection: opts.deflection !== undefined ? opts.deflection / S : undefined, angularStep: opts.angularStep, name });
  const models = [];
  for (const b of res.bodies) {
    const c = b.coords;
    for (let i = 0; i < c.length; i++) c[i] *= S;
    models.push(surfaceModel(b.name, c, b.indices, b.groups));
  }
  if (res.wires.length) {
    let lo = Infinity, hi = -Infinity;
    for (const w of res.wires) for (const p of [w.a, w.b]) for (const v of p) { if (v * S < lo) lo = v * S; if (v * S > hi) hi = v * S; }
    const weld = new Welder(Math.max(hi - lo, 1e-12) * 1e-9, res.wires.length * 2);
    const members = [];
    for (const w of res.wires) {
      const a = weld.add(w.a[0] * S, w.a[1] * S, w.a[2] * S), b = weld.add(w.b[0] * S, w.b[1] * S, w.b[2] * S);
      if (a !== b) members.push({ a, b, group: w.group });
    }
    if (members.length) { models.push({ kind: 'frame', name: cleanName(models.length ? `${name} (wires)` : name), nodes: weld.finish(), members }); tally(res.hist, 'wire edge', res.wires.length); }
  }
  // the unit note goes after the notes on what was read, as before
  const warnings = res.warnings.slice();
  const at = warnings.findIndex((w) => /^The header announces|the faces enclose a negative volume/.test(w));
  warnings.splice(at < 0 ? warnings.length : at, 0, ...lead);
  if (!models.length) throw new Error(`${fname}: no geometry could be read from the ${label} file.${warnings.length ? ` ${warnings.join(' ')}` : ''}`);
  return finish(models, warnings, {
    format, units, crs: null, elementTypes: res.hist, skipped: res.skipped,
    acisVersion: header.version, product: header.product, acisRelease: header.acis, saved: header.date,
    bodies: res.stats.bodies, closedBodies: res.stats.closedBodies, orientedBodies: res.bodies.filter((b) => b.closed && b.oriented).length,
    solidBodies: res.bodies.filter((b) => !b.freeEdges && !b.nonManifoldEdges).length,
    faces: res.stats.facesRead, faceTotal: res.stats.faces, facesSkipped: res.stats.facesSkipped,
    volumes: res.bodies.map((b) => (b.volume === null ? null : b.volume * S * S * S)), resabs: header.resabs,
    tessellation: { chord: res.stats.chord * S, angularStepDeg: res.stats.angularStepDeg },
  });
}

/**
 * ACIS SAT (text).
 * @param {{ name?: string, text?: string, bytes?: Uint8Array }} input
 * @param {{ units?: string, scale?: number, deflection?: number, angularStep?: number }} [opts] `deflection`:
 *   largest distance between a curved edge or surface and its facets [m] (default 1 % of the bounding-box
 *   diagonal of each body); `angularStep`: largest angle one facet spans [deg], default 15
 */
export function parseAcis(input, opts = {}) {
  const fname = input.name || 'SAT';
  if (input.bytes && input.bytes.length >= 15) {
    const sig = String.fromCharCode(...input.bytes.subarray(0, 15));
    if (sig === 'ACIS BinaryFile' || sig.startsWith('ASM BinaryFile')) throw new Error(`${fname}: this is a binary ACIS file (SAB), which the SAT (text) reader does not take. Open it as ACIS SAB (importers/sab.js reads it).`);
  }
  return acisModels(readSat(inputText(input), fname), input, opts, 'acis');
}
