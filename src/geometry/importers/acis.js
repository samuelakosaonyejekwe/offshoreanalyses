// ACIS SAT (text) — boundary representation with analytic geometry.
//
// File structure (as published for the SAT save format):
//   line 1   version (major × 100 + minor), record count, body count, history flag
//   line 2   product, ACIS version and date as counted strings          (version 2.0 and later)
//   line 3   millimetres per model unit, resabs, resnor                 (version 2.0 and later)
//   records  "<type> <fields …> #", optionally prefixed by a sequence number "-n". Pointers are "$n" (record
//            index, −1 = null). Up to 6.x a record starts with its attribute pointer; from 7.0 on that is
//            followed by an integer and a history pointer ("body $-1 -1 $-1 …") and strings are "@len text".
//
// Read: body → lump → shell → face → loop → coedge → edge → vertex → point, with the body transform.
//   surfaces  plane-surface, cone-surface (cylinders and cones, elliptical ratio), sphere-surface, torus-surface
//   curves    straight-curve, ellipse-curve (circles, ellipses and their arcs), intcurve-curve when the record
//             stores its B-spline (exactcur, or the stored approximation of a procedural curve)
// Each face is triangulated in the parameter space of its surface between its boundary loops, the boundary
// points being those of the shared edges, so neighbouring faces meet without gaps.
//
// Conventions relied on (stated here because they decide the geometry):
//   * cone-surface: radius at axial distance z from the base ellipse is R + z·sin/cos; the surface normal points
//     away from the axis when cos > 0 and towards it when cos < 0. sphere / torus: a negative (minor) radius
//     turns the normal inwards.
//   * a face whose sense is "reversed" has its normal opposite to the surface normal; loops run counter-clockwise
//     about the face normal (material on the left).
//   * transform: row-vector convention, p′ = scale·(p·A) + t.
//   * B-spline knot lists hold (value, multiplicity) pairs with one end knot fewer than the usual clamped
//     vector at each end; rational control points are stored as x y z w, not premultiplied.
//
// Not read: spline-surface faces without stored control points (and other procedural surfaces) — counted and
// reported; trimming of explicit B-spline surfaces; pcurves, attributes, subshells, history. SAB (binary) is
// not read. Record types that are not recognised are reported, never interpreted.

import { inputText, baseName, finish, surfaceModel, Welder, V3, tally, cleanName, Warnings } from './common.js';
import { resolveUnits } from '../units.js';
import { sampleCurve, tessellateSurface } from './nurbs.js';
import { triangulateSection } from '../../mesh/generate.js';

const TAU = 2 * Math.PI, D2R = Math.PI / 180;
const NUM = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const UNIT_BY_MM = [[1, 'mm'], [10, 'cm'], [100, 'dm'], [1000, 'm'], [1e6, 'km'], [25.4, 'in'], [304.8, 'ft'], [914.4, 'yd'], [0.001, 'um']];
// record kinds that carry no geometry of their own and are passed over without comment
const QUIET = new Set(['attrib', 'pcurve', 'asmheader', 'subshell', 'wire', 'lump', 'shell', 'body', 'loop', 'coedge', 'edge', 'vertex', 'point', 'face', 'transform', 'curve', 'surface', 'eye_refinement', 'vertex_template', 'history_stream', 'delta_state', 'annotation']);

const isWs = (c) => c === 32 || (c >= 9 && c <= 13);
const lerp3 = (a, b, f) => [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
const dist3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
const nearest = (u, ref) => u + TAU * Math.round((ref - u) / TAU);

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
    if (/^(End-of-(ACIS|ASM)-data|Begin-of-(ACIS|ASM)-History-Data)$/i.test(type)) break;
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
  if (!isN(tok[i + 1]) || !(tok[i + 2] && tok[i + 2].t === 'w') || !isN(tok[i + 3])) return null;
  const deg = tok[i + 1].v, closure = tok[i + 2].v;
  if (!(deg >= 1 && deg <= 25) || closure === 'periodic') return null;
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
  return { deg, U: k.U, P, W };
}

/** B-spline surface stored at tok[i] = nubs | nurbs; control points kept in file order with both counts. */
function bs3Surface(tok, i) {
  const rational = tok[i].v === 'nurbs';
  if (!isN(tok[i + 1]) || !isN(tok[i + 2])) return null;
  const pu = tok[i + 1].v, pv = tok[i + 2].v;
  let j = i + 3;
  const words = [];
  while (tok[j] && tok[j].t === 'w' && words.length < 4) words.push(tok[j++].v);
  if (words.length !== 4 || words[0] === 'periodic' || words[1] === 'periodic') return null;
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

/** Sub-type block of a record: its name and the position of a stored B-spline, if one follows the name. */
function subtype(tok) {
  const a = tok.findIndex((t) => t.t === '{');
  if (a < 0) return null;
  const name = tok[a + 1] && tok[a + 1].t === 'w' ? tok[a + 1].v : '';
  let j = a + 2;
  if (isW(tok[j], 'full')) j++;
  return { name, at: isW(tok[j], 'nubs', 'nurbs') ? j : -1, reversed: tok.slice(0, a).some((t) => isW(t, 'reversed')) };
}

// ---------------------------------------------------------------------------------------------------------
// Analytic surfaces: 3-D ↔ parameter space. Parameters are chosen so that ∂u × ∂v points away from the axis /
// centre (plane: along its normal); `inward` records that the ACIS surface normal is the opposite.
// ---------------------------------------------------------------------------------------------------------

function anyPerp(n) {
  const s = Math.abs(n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  return V3.norm(V3.sub(s, V3.scale(n, V3.dot(s, n))));
}
function perpTo(x, n) {
  if (!x || !(V3.len(x) > 0)) return anyPerp(n);
  const p = V3.sub(x, V3.scale(n, V3.dot(x, n)));
  return V3.len(p) > 1e-9 * V3.len(x) ? V3.norm(p) : anyPerp(n);
}

function makeSurface(type, n, stepFor) {
  const ok = (k) => n.length >= k && n.slice(0, k).every(Number.isFinite);
  if (type === 'plane-surface') {
    if (!ok(6)) return null;
    const O = n.slice(0, 3), N = V3.norm(n.slice(3, 6));
    if (!(V3.len(N) > 0)) return null;
    const e1 = perpTo(n.length >= 9 ? n.slice(6, 9) : null, N), e2 = V3.cross(N, e1);
    return {
      type: 'plane', inward: false, perU: false, perV: false,
      toParam: (p) => { const d = V3.sub(p, O); return [V3.dot(d, e1), V3.dot(d, e2)]; },
      fromParam: (u, v) => [O[0] + u * e1[0] + v * e2[0], O[1] + u * e1[1] + v * e2[1], O[2] + u * e1[2] + v * e2[2]],
      scales: () => ({ su: 1, sv: 1, size: undefined }),
    };
  }
  if (type === 'cone-surface') {
    if (!ok(12)) return null;
    const C = n.slice(0, 3), A = V3.norm(n.slice(3, 6)), M = n.slice(6, 9), ratio = n[9], sin = n[10], cos = n[11];
    const R0 = V3.len(M);
    if (!(V3.len(A) > 0) || !(R0 > 0) || !(ratio > 0) || !(Math.abs(cos) > 1e-12) || Math.abs(Math.hypot(sin, cos) - 1) > 1e-3) return null;
    const X = perpTo(M, A), Y = V3.cross(A, X), tan = sin / cos;
    const eps = 1e-9 * R0;
    const surf = {
      type: Math.abs(sin) < 1e-12 ? 'cylinder' : 'cone', inward: cos < 0, perU: true, perV: false, straightV: true,
      apex: Math.abs(tan) > 1e-12 ? -R0 / tan : null,
      toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y) / ratio; return [Math.hypot(x, y) > eps ? Math.atan2(y, x) : NaN, V3.dot(d, A)]; },
      fromParam: (u, v) => { const r = R0 + v * tan, cx = r * Math.cos(u), cy = r * ratio * Math.sin(u); return [C[0] + cx * X[0] + cy * Y[0] + v * A[0], C[1] + cx * X[1] + cy * Y[1] + v * A[1], C[2] + cx * X[2] + cy * Y[2] + v * A[2]]; },
      radiusAt: (v) => Math.abs(R0 + v * tan) * Math.max(1, ratio),
      // generators are straight: only the angle a triangle spans matters, so the axial extent of a face is
      // squeezed below one target size (see `squeeze`) and triangles run from one end of the face to the other
      scales: (rmax) => { const st = stepFor(rmax); return { su: 1 / st, sv: 1 / (8 * rmax * st), size: 1, stepU: st, stepV: Infinity, squeeze: 0.6 }; },
    };
    return surf;
  }
  if (type === 'sphere-surface') {
    if (!ok(4)) return null;
    const C = n.slice(0, 3), R = Math.abs(n[3]);
    if (!(R > 0)) return null;
    const Z = n.length >= 10 && V3.len(n.slice(7, 10)) > 0 ? V3.norm(n.slice(7, 10)) : [0, 0, 1];
    const X = perpTo(n.length >= 7 ? n.slice(4, 7) : null, Z), Y = V3.cross(Z, X);
    return {
      type: 'sphere', inward: n[3] < 0, perU: true, perV: false, poles: true, C, R,
      toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y), z = V3.dot(d, Z); const h = Math.hypot(x, y); return [h > 1e-9 * R ? Math.atan2(y, x) : NaN, Math.atan2(z, h)]; },
      fromParam: (u, v) => { const c = R * Math.cos(v), cx = c * Math.cos(u), cy = c * Math.sin(u), cz = R * Math.sin(v); return [C[0] + cx * X[0] + cy * Y[0] + cz * Z[0], C[1] + cx * X[1] + cy * Y[1] + cz * Z[1], C[2] + cx * X[2] + cy * Y[2] + cz * Z[2]]; },
      scales: () => { const st = stepFor(R); return { su: 1 / st, sv: 1 / st, size: 1, stepU: st, stepV: st }; },
    };
  }
  if (type === 'torus-surface') {
    if (!ok(8)) return null;
    const C = n.slice(0, 3), N = V3.norm(n.slice(3, 6)), R = n[6], r = Math.abs(n[7]);
    // spindle ("apple" / "lemon") tori, where the tube crosses the axis, are not handled
    if (!(V3.len(N) > 0) || !(r > 0) || !(R > r)) return null;
    const X = perpTo(n.length >= 11 ? n.slice(8, 11) : null, N), Y = V3.cross(N, X);
    return {
      type: 'torus', inward: n[7] < 0, perU: true, perV: true, C, R, r,
      toParam: (p) => { const d = V3.sub(p, C); const x = V3.dot(d, X), y = V3.dot(d, Y); return [Math.atan2(y, x), Math.atan2(V3.dot(d, N), Math.hypot(x, y) - R)]; },
      fromParam: (u, v) => { const c = R + r * Math.cos(v), cx = c * Math.cos(u), cy = c * Math.sin(u), cz = r * Math.sin(v); return [C[0] + cx * X[0] + cy * Y[0] + cz * N[0], C[1] + cx * X[1] + cy * Y[1] + cz * N[1], C[2] + cx * X[2] + cy * Y[2] + cz * N[2]]; },
      // around the axis the facets may not sag more than they do around the tube, or the chords of a slender
      // torus would cut through it
      scales: () => { const b = stepFor(r), a = Math.min(stepFor(R + r), 2 * Math.acos(1 - (r * (1 - Math.cos(b / 2))) / (R + r))); return { su: 1 / a, sv: 1 / b, size: 1, stepU: a, stepV: b }; },
    };
  }
  return null;
}

// ---------------------------------------------------------------------------------------------------------
// Face domain in parameter space
// ---------------------------------------------------------------------------------------------------------

const areaOf = (L) => { let a = 0; for (let i = 0, j = L.length - 1; i < L.length; j = i++) a += L[j].x * L[i].y - L[i].x * L[j].y; return a / 2; };
const meanOf = (L, k) => { let s = 0; for (const q of L) s += q[k]; return s / L.length; };
const shifted = (L, dx, dy = 0) => (dx || dy ? L.map((q) => ({ x: q.x + dx, y: q.y + dy, p: q.p })) : L);

/**
 * Boundary loops of a face in the parameter plane of a periodic surface, unwrapped. A vertex at a pole or apex
 * (where u is undefined) becomes the stretch of the degenerate row between the u of its neighbours.
 * @returns {{ pts: Array<{x, y, p}>, wu: number, wv: number }|null}
 */
function unwrapLoop(surf, pts3) {
  const raw = pts3.map((p) => { const [u, v] = surf.toParam(p); return { u, v, p, sing: u !== u }; });
  const m = raw.length;
  if (raw.every((q) => q.sing)) return null;
  const out = [];
  const firstSing = raw.findIndex((q, i) => q.sing && !raw[(i + 1) % m].sing);
  const fixV = (v, ref) => (surf.perV ? nearest(v, ref) : v);
  if (firstSing < 0) {
    let lu = raw[0].u, lv = raw[0].v;
    out.push({ x: lu, y: lv, p: raw[0].p });
    for (let i = 1; i < m; i++) { lu = nearest(raw[i].u, lu); lv = fixV(raw[i].v, lv); out.push({ x: lu, y: lv, p: raw[i].p }); }
    return { pts: out, wu: Math.round((nearest(raw[0].u, lu) - raw[0].u) / TAU), wv: surf.perV ? Math.round((nearest(raw[0].v, lv) - raw[0].v) / TAU) : 0 };
  }
  // start just after a degenerate vertex so that the loop closes along the degenerate row without winding
  const seq = [];
  for (let k = 1; k <= m; k++) seq.push(raw[(firstSing + k) % m]);
  let lu = seq[0].u, lv = seq[0].v;
  out.push({ x: lu, y: lv, p: seq[0].p });
  for (let i = 1; i < m; i++) {
    const q = seq[i];
    if (!q.sing) { lu = nearest(q.u, lu); lv = fixV(q.v, lv); out.push({ x: lu, y: lv, p: q.p }); continue; }
    out.push({ x: lu, y: q.v, p: q.p });
    let nx = null;
    for (let k = i + 1; k < m && !nx; k++) if (!seq[k].sing) nx = seq[k];
    const un = nx ? nearest(nx.u, lu) : out[0].x;
    if (un !== lu) out.push({ x: un, y: q.v, p: q.p });
    lu = un;
  }
  return { pts: out, wu: 0, wv: 0 };
}

/**
 * Polygons (even–odd rule) bounding a face on a periodic surface.
 * @param {boolean} flip face normal opposes ∂u × ∂v, so loops run clockwise in the parameter plane
 * @returns {Array<Array<{x, y, p}>>|string} loops, or the reason the domain could not be built
 */
function periodicDomain(surf, loops, flip, sc) {
  const orient = flip ? -1 : 1;
  const mk = (x, y) => ({ x, y, p: surf.fromParam(x, y) });
  const copy = (q, dx, dy) => ({ x: q.x + dx, y: q.y + dy, p: q.p });
  const divs = (a, b) => Math.max(1, Math.ceil(Math.max(Math.abs(b.x - a.x) / sc.stepU, Number.isFinite(sc.stepV) ? Math.abs(b.y - a.y) / sc.stepV : 0) - 1e-9));
  // interior points of the straight parameter-space line a → b (shared by both sides of a seam)
  const between = (a, b) => { const k = divs(a, b), o = []; for (let i = 1; i < k; i++) o.push(mk(a.x + ((b.x - a.x) * i) / k, a.y + ((b.y - a.y) * i) / k)); return o; };
  const wind = loops.filter((l) => l.wu || l.wv), flat = loops.filter((l) => !l.wu && !l.wv);
  if (wind.some((l) => l.wu && l.wv)) return 'a boundary loop winds around both directions of the surface';

  // full parameter rectangle of a closed surface, cut where no hole lies
  const rectangle = (holes) => {
    const place = (list, k, lo) => list.map((h) => { const d = -TAU * Math.floor((meanOf(h, k) - lo) / TAU); return k === 'x' ? shifted(h, d, 0) : shifted(h, 0, d); });
    const clear = (hs, k, lo) => hs.every((h) => h.every((q) => q[k] > lo + 1e-6 && q[k] < lo + TAU - 1e-6));
    let u0 = null, hs = holes;
    for (let c = 0; c < 16 && u0 === null; c++) { const lo = -Math.PI + 0.0137 + (c * Math.PI) / 8; const t = place(holes, 'x', lo); if (clear(t, 'x', lo)) { u0 = lo; hs = t; } }
    if (u0 === null) return 'its holes leave no meridian free for the seam';
    if (surf.poles) {
      const s = mk(u0, -Math.PI / 2), nth = mk(u0, Math.PI / 2);
      const side = between(s, nth);
      return [[s, copy(s, TAU, 0), ...side.map((q) => copy(q, TAU, 0)), copy(nth, TAU, 0), nth, ...side.slice().reverse()], ...hs];
    }
    let v0 = null;
    for (let c = 0; c < 16 && v0 === null; c++) { const lo = -Math.PI + 0.0137 + (c * Math.PI) / 8; const t = place(hs, 'y', lo); if (clear(t, 'y', lo)) { v0 = lo; hs = t; } }
    if (v0 === null) return 'its holes leave no parallel free for the seam';
    const o = mk(u0, v0), bu = between(o, { x: u0 + TAU, y: v0 }), bv = between(o, { x: u0, y: v0 + TAU });
    return [[o, ...bu, copy(o, TAU, 0), ...bv.map((q) => copy(q, TAU, 0)), copy(o, TAU, TAU), ...bu.slice().reverse().map((q) => copy(q, 0, TAU)), copy(o, 0, TAU), ...bv.slice().reverse()], ...hs];
  };

  if (!wind.length) {
    if (!flat.length) return surf.poles || surf.perV ? rectangle([]) : 'the face has no boundary on an unbounded surface';
    const withArea = flat.map((l) => ({ l, a: areaOf(l.pts) }));
    if (withArea.some((q) => q.a * orient > 0) || !(surf.poles || surf.perV)) {
      withArea.sort((p, q) => Math.abs(q.a) - Math.abs(p.a));
      const ref = withArea[0].l.pts, mu = meanOf(ref, 'x'), mv = meanOf(ref, 'y');
      return withArea.map((q) => shifted(q.l.pts, TAU * Math.round((mu - meanOf(q.l.pts, 'x')) / TAU), surf.perV ? TAU * Math.round((mv - meanOf(q.l.pts, 'y')) / TAU) : 0));
    }
    return rectangle(flat.map((l) => l.pts));
  }

  // work with "along" = the winding direction; swap coordinates for loops that wind in v
  const alongV = !!wind[0].wv;
  if (wind.some((l) => !!l.wv !== alongV)) return 'its boundary loops wind around different directions of the surface';
  const sw = (L) => (alongV ? L.map((q) => ({ x: q.y, y: q.x, p: q.p })) : L);
  const back = (polys) => (alongV ? polys.map(sw) : polys);
  const mkA = (x, y) => (alongV ? { x, y, p: surf.fromParam(y, x) } : mk(x, y));
  const betweenA = (a, b) => (alongV ? between({ x: a.y, y: a.x }, { x: b.y, y: b.x }).map((q) => ({ x: q.y, y: q.x, p: q.p })) : between(a, b));
  const A = sw(wind[0].pts), wA = alongV ? wind[0].wv : wind[0].wu;
  const crossPeriodic = alongV ? surf.perU : surf.perV;
  // side of loop A on which the face lies, as a sign on the cross coordinate
  const side = Math.sign(wA) * (alongV ? -1 : 1) * orient;
  const lo = Math.min(A[0].x, A[0].x + TAU * wA);
  const holes = flat.map((l) => { const h = sw(l.pts); return shifted(h, -TAU * Math.floor((meanOf(h, 'x') - lo) / TAU), 0); });
  if (holes.some((h) => h.some((q) => !(q.x > lo + 1e-6 && q.x < lo + TAU - 1e-6)))) return 'an inner loop crosses the seam chosen for the face';
  const a0 = A[0], aEnd = copy(a0, TAU * wA, 0);

  if (wind.length === 1) {
    if (alongV) return 'a single boundary loop winds around the tube of the surface';
    let vp;
    if (surf.poles) vp = side * Math.PI / 2;
    else if (surf.apex !== null && surf.apex !== undefined) vp = surf.apex;
    else return 'a single closed boundary leaves the face unbounded along the axis';
    const pole = mkA(a0.x, vp);
    const seam = surf.straightV ? [] : betweenA(a0, pole);
    return back([[...A, aEnd, ...seam.map((q) => copy(q, TAU * wA, 0)), copy(pole, TAU * wA, 0), pole, ...seam.slice().reverse()], ...holes]);
  }
  if (wind.length !== 2) return `${wind.length} boundary loops wind around the surface`;
  let B = sw(wind[1].pts);
  const wB = alongV ? wind[1].wv : wind[1].wu;
  if (Math.sign(wB) === Math.sign(wA)) B = B.slice().reverse();
  // start B at its vertex nearest (in angle) to the start of A, so that the seam is as short as possible
  let best = 0, bd = Infinity;
  for (let k = 0; k < B.length; k++) { const d = Math.abs(nearest(B[k].x, a0.x) - a0.x); if (d < bd) { bd = d; best = k; } }
  const base = B[best].x;
  // B runs against A: rebuild it continuously from its new start
  const Br = [];
  for (let k = 0; k < B.length; k++) { const q = B[(best + k) % B.length]; let x = q.x; if (k) { x = nearest(q.x, Br[k - 1].x); } Br.push({ x, y: q.y, p: q.p }); }
  let dx = nearest(base, a0.x) - base, dy = 0;
  if (crossPeriodic) { const ca = meanOf(A, 'y'), cb = meanOf(Br, 'y'); let d = side * (cb - ca); d -= TAU * Math.floor(d / TAU); dy = ca + side * d - cb; }
  const Bs = shifted(Br, dx, dy);
  const b0 = Bs[0];
  const seam = betweenA(a0, b0);
  return back([[...A, aEnd, ...seam.map((q) => copy(q, TAU * wA, 0)), ...Bs.map((q) => copy(q, TAU * wA, 0)), b0, ...seam.slice().reverse()], ...holes]);
}

/**
 * Triangulate a face domain and lift it to 3-D. Vertices the triangulation adds on a boundary segment are put on
 * the chord between the 3-D end points of that segment, so that a neighbouring face can be matched to them.
 */
function liftDomain(surf, polys, sc) {
  const loops = [];
  const input = [];
  for (const L of polys) {
    const flat = [];
    const kept = [];
    for (const q of L) {
      const X = q.x * sc.su, Y = q.y * sc.sv;
      if (kept.length && kept[kept.length - 1].X === X && kept[kept.length - 1].Y === Y) continue;
      kept.push({ X, Y, p: q.p });
    }
    if (kept.length > 1 && kept[0].X === kept[kept.length - 1].X && kept[0].Y === kept[kept.length - 1].Y) kept.pop();
    if (kept.length < 3) continue;
    for (const k of kept) flat.push(k.X, k.Y);
    loops.push(Float64Array.from(flat));
    input.push(kept);
  }
  if (!loops.length) return null;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const L of input) for (const q of L) { if (q.X < x0) x0 = q.X; if (q.X > x1) x1 = q.X; if (q.Y < y0) y0 = q.Y; if (q.Y > y1) y1 = q.Y; }
  const span = Math.max(x1 - x0, y1 - y0);
  if (!(span > 0)) return null;
  const tri = triangulateSection({ loops }, { minAngle: 0, size: sc.size });
  if (!tri.indices.length) return null;
  const cell = span * 1e-7;
  const grid = new Map();
  for (const L of input) for (const q of L) grid.set(`${Math.round(q.X / cell)},${Math.round(q.Y / cell)}`, q);
  const onBoundary = new Uint8Array(tri.xy.length / 2);
  for (let i = 0; i < tri.boundary.length; i++) onBoundary[tri.boundary[i]] = 1;
  const pts = [];
  for (let i = 0; i < tri.xy.length / 2; i++) {
    const X = tri.xy[2 * i], Y = tri.xy[2 * i + 1];
    const gx = Math.round(X / cell), gy = Math.round(Y / cell);
    let hit = null;
    for (let a = -1; a <= 1 && !hit; a++) for (let b = -1; b <= 1 && !hit; b++) { const q = grid.get(`${gx + a},${gy + b}`); if (q && Math.abs(q.X - X) <= cell && Math.abs(q.Y - Y) <= cell) hit = q; }
    if (hit) { pts.push(hit.p); continue; }
    if (onBoundary[i]) {
      let bd = Infinity, bp = null;
      for (const L of input) for (let k = 0, j = L.length - 1; k < L.length; j = k++) {
        const a = L[j], b = L[k], ex = b.X - a.X, ey = b.Y - a.Y, l2 = ex * ex + ey * ey;
        const f = l2 > 0 ? Math.max(0, Math.min(1, ((X - a.X) * ex + (Y - a.Y) * ey) / l2)) : 0;
        const d = Math.hypot(X - a.X - f * ex, Y - a.Y - f * ey);
        if (d < bd) { bd = d; bp = a.p === b.p ? a.p : lerp3(a.p, b.p, f); }
      }
      pts.push(bp);
      continue;
    }
    pts.push(surf.fromParam(X / sc.su, Y / sc.sv));
  }
  return { pts, indices: tri.indices };
}

// ---------------------------------------------------------------------------------------------------------
// Mesh clean-up shared by all faces of a body
// ---------------------------------------------------------------------------------------------------------

/**
 * Close T-junctions: wherever a vertex lies on an edge that only one triangle uses, that triangle is split there.
 * They arise where the triangulation of one face subdivided a shared boundary segment and its neighbour did not.
 * @param {Float64Array} P welded coordinates @param {number[]} T triangles (flat) @param {number[]} F face id per triangle
 */
function closeTJunctions(P, T, F, tol) {
  for (let pass = 0; pass < 12; pass++) {
    const nV = P.length / 3;
    const use = new Map();
    for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) {
      const a = T[t + k], b = T[t + (k + 1) % 3];
      const key = a < b ? a * nV + b : b * nV + a;
      const e = use.get(key);
      if (e) e.n++; else use.set(key, { n: 1, t, k });
    }
    const open = [];
    const bv = new Set();
    for (const e of use.values()) if (e.n === 1) { open.push(e); bv.add(T[e.t + e.k]); bv.add(T[e.t + (e.k + 1) % 3]); }
    if (!open.length) return;
    const cand = [...bv];
    const split = new Map(); // triangle offset → { k, mids }
    for (const e of open) {
      if (split.has(e.t)) continue;
      const a = T[e.t + e.k], b = T[e.t + (e.k + 1) % 3];
      const ax = P[3 * a], ay = P[3 * a + 1], az = P[3 * a + 2];
      const ex = P[3 * b] - ax, ey = P[3 * b + 1] - ay, ez = P[3 * b + 2] - az;
      const l2 = ex * ex + ey * ey + ez * ez;
      if (!(l2 > 0)) continue;
      const mids = [];
      for (const c of cand) {
        if (c === a || c === b) continue;
        const dx = P[3 * c] - ax, dy = P[3 * c + 1] - ay, dz = P[3 * c + 2] - az;
        const f = (dx * ex + dy * ey + dz * ez) / l2;
        if (!(f > 0 && f < 1)) continue;
        if (Math.hypot(dx - f * ex, dy - f * ey, dz - f * ez) <= tol && f * Math.sqrt(l2) > tol && (1 - f) * Math.sqrt(l2) > tol) mids.push([f, c]);
      }
      if (mids.length) { mids.sort((p, q) => p[0] - q[0]); split.set(e.t, { k: e.k, mids: mids.map((q) => q[1]) }); }
    }
    if (!split.size) return;
    const nT = [], nF = [];
    for (let t = 0; t < T.length; t += 3) {
      const s = split.get(t);
      if (!s) { nT.push(T[t], T[t + 1], T[t + 2]); nF.push(F[t / 3]); continue; }
      const a = T[t + s.k], b = T[t + (s.k + 1) % 3], o = T[t + (s.k + 2) % 3];
      const chain = [a, ...s.mids, b];
      for (let i = 0; i + 1 < chain.length; i++) { nT.push(chain[i], chain[i + 1], o); nF.push(F[t / 3]); }
    }
    T.length = 0; F.length = 0;
    for (const v of nT) T.push(v);
    for (const v of nF) F.push(v);
  }
}

/** True when every edge is used by exactly two triangles. */
function isClosed(T, nV) {
  const use = new Map();
  for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) { const a = T[t + k], b = T[t + (k + 1) % 3]; const key = a < b ? a * nV + b : b * nV + a; use.set(key, (use.get(key) || 0) + 1); }
  for (const n of use.values()) if (n !== 2) return false;
  return use.size > 0;
}

function volumeOf(P, T) {
  let s = 0;
  for (let t = 0; t < T.length; t += 3) {
    const a = 3 * T[t], b = 3 * T[t + 1], c = 3 * T[t + 2];
    s += P[a] * (P[b + 1] * P[c + 2] - P[b + 2] * P[c + 1]) - P[a + 1] * (P[b] * P[c + 2] - P[b + 2] * P[c]) + P[a + 2] * (P[b] * P[c + 1] - P[b + 1] * P[c]);
  }
  return s / 6;
}

// ---------------------------------------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------------------------------------

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
    if (sig === 'ACIS BinaryFile' || sig.startsWith('ASM BinaryFile')) throw new Error(`${fname}: this is a binary ACIS file (SAB), which is not read. Save it as SAT (text) or STEP in the originating program.`);
  }
  if (opts.deflection !== undefined && !(opts.deflection > 0)) throw new RangeError('deflection: must be a positive length in metres.');
  if (opts.angularStep !== undefined && !(opts.angularStep > 0 && opts.angularStep <= 90)) throw new RangeError('angularStep: must be within 0…90 degrees.');
  const { header, records } = readSat(inputText(input), fname);
  const name = baseName(input.name);
  const W = new Warnings();
  const extra = [];
  const skipped = {}, hist = {};
  const modern = header.version >= 700;

  // ---- units ----
  let units;
  const mm = header.mmPerUnit;
  if (mm > 0) {
    const known = UNIT_BY_MM.find(([v]) => Math.abs(v - mm) <= 1e-9 * v);
    units = known || opts.units || (Number.isFinite(opts.scale) && opts.scale > 0) ? resolveUnits(opts, known ? known[1] : null, 'mm') : { unit: 'custom', scale: mm / 1000, basis: 'file' };
  } else {
    units = resolveUnits(opts, null, 'mm');
    if (units.basis === 'assumed') extra.push('The file declares no length unit; millimetres were assumed. Re-import with units set if the model was drawn in another unit.');
  }
  const S = units.scale;

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
  const allOf = (r, kind) => { const o = []; for (const t of dataOf(r)) if (t.t === 'p') { const q = t.v >= 0 ? byIdx.get(t.v) : null; if (q && q.kind === kind) o.push(q); } return o; };
  const reversedWord = (r) => { for (const t of dataOf(r)) { if (t.t === '{') break; if (t.t === 'w' && (t.v === 'forward' || t.v === 'reversed')) return t.v === 'reversed'; } return false; };
  // numbers outside sub-type blocks; "F <value>" interval bounds are not data
  const numsOf = (r) => {
    const o = [], d = dataOf(r);
    let depth = 0;
    for (let i = 0; i < d.length; i++) {
      if (d[i].t === '{') depth++; else if (d[i].t === '}') depth--;
      else if (!depth && d[i].t === 'n' && !(i && isW(d[i - 1], 'F'))) o.push(d[i].v);
    }
    return o;
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

  const bodies = records.filter((r) => r.kind === 'body');
  const roots = bodies.length ? bodies.map((b) => ({ body: b, recs: owned(b) })) : [{ body: null, recs: records }];
  const models = [];
  const anyFace = records.some((r) => r.kind === 'face');
  let freeForm = 0, untrimmed = 0, facesRead = 0, wireEdges = 0, closedBodies = 0;
  const frameNodes = [], frameMembers = [];
  let chordM = 0;
  const angStep = (opts.angularStep ?? 15) * D2R;

  roots.forEach((root, bi) => {
    const bodyName = roots.length > 1 ? `${name} body ${bi + 1}` : name;
    const tf = transformOf(root.body);
    // extent of the body for the default deflection: its vertices and the spheres, tori and circles it holds
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    const see = (p, r = 0) => { for (let k = 0; k < 3; k++) { if (p[k] - r < lo[k]) lo[k] = p[k] - r; if (p[k] + r > hi[k]) hi[k] = p[k] + r; } };
    for (const r of root.recs) {
      if (r.kind === 'vertex') { const p = pointOf(r); if (p) see(p); } else if (r.type === 'sphere-surface') { const n = numsOf(r); if (n.length >= 4) see(n, Math.abs(n[3])); } else if (r.type === 'torus-surface') { const n = numsOf(r); if (n.length >= 8) see(n, Math.abs(n[6]) + Math.abs(n[7])); } else if (r.type === 'ellipse-curve') { const n = numsOf(r); if (n.length >= 9) see(n, Math.hypot(n[6], n[7], n[8])); }
    }
    const diag = lo[0] === Infinity ? 0 : Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    if (!(diag > 0)) {
      // nothing with an extent: only faces on surfaces that are not read, or records without coordinates
      for (const r of root.recs) {
        if (r.kind !== 'face') continue;
        const sr = firstOf(r, 'surface');
        if (sr && !/^(plane|cone|sphere|torus)-surface$/.test(sr.type)) freeForm++; else W.add('A face without any vertex or extent was skipped');
      }
      return;
    }
    const defl = opts.deflection !== undefined ? opts.deflection / S : 0.01 * diag;
    chordM = Math.max(chordM, defl * S);
    const stepFor = (r) => {
      let s = angStep;
      if (r > defl) s = Math.min(s, 2 * Math.acos(1 - defl / r));
      return Math.max(s, 0.25 * D2R);
    };
    const tol = diag * 1e-9;

    // ---- edges: one polyline per edge, shared by every face that uses it ----
    const edgePts = (e) => {
      if (e.pts) return e.pts;
      const vs = allOf(e, 'vertex');
      const Sp = pointOf(vs[0]), Ep = pointOf(vs.length > 1 ? vs[1] : vs[0]);
      const cv = firstOf(e, 'curve');
      const rev = reversedWord(e);
      const chord = (why) => { if (why) W.add(why); return Sp && Ep && Sp !== Ep ? [Sp, Ep] : Sp ? [Sp] : []; };
      if (!cv) return (e.pts = Sp ? [Sp] : []);
      const same = !Sp || !Ep || Sp === Ep || dist3(Sp, Ep) <= tol;
      if (cv.type === 'straight-curve') return (e.pts = chord());
      if (cv.type === 'ellipse-curve') {
        const n = numsOf(cv);
        if (n.length < 9 || !n.slice(0, 9).every(Number.isFinite)) return (e.pts = chord('An ellipse-curve record could not be read; its edge was replaced by a straight chord'));
        const C = n.slice(0, 3), N = V3.norm(n.slice(3, 6)), M = n.slice(6, 9), ratio = n.length >= 10 && n[9] > 0 ? n[9] : 1;
        const Y = V3.scale(V3.cross(N, M), ratio), m2 = V3.dot(M, M), y2 = V3.dot(Y, Y);
        if (!(m2 > 0) || !(y2 > 0)) return (e.pts = chord('An ellipse-curve record has no radius; its edge was replaced by a straight chord'));
        const ang = (p) => { const d = V3.sub(p, C); return Math.atan2(V3.dot(d, Y) / y2, V3.dot(d, M) / m2); };
        const dir = rev ? -1 : 1;
        const a0 = Sp ? ang(Sp) : 0;
        let sweep = dir * TAU;
        if (!same) { sweep = ang(Ep) - a0; if (dir > 0) { while (sweep <= 1e-12) sweep += TAU; } else { while (sweep >= -1e-12) sweep -= TAU; } }
        const k = Math.max(same ? 3 : 1, Math.ceil(Math.abs(sweep) / stepFor(Math.sqrt(Math.max(m2, y2))) - 1e-9));
        const pts = [];
        for (let i = 0; i <= k; i++) {
          if (i === 0 && Sp) pts.push(Sp);
          else if (i === k && (same ? Sp : Ep)) pts.push(same ? Sp : Ep);
          else { const t = a0 + (sweep * i) / k, c = Math.cos(t), s = Math.sin(t); pts.push([C[0] + c * M[0] + s * Y[0], C[1] + c * M[1] + s * Y[1], C[2] + c * M[2] + s * Y[2]]); }
        }
        if (same && !Sp) pts[k] = pts[0];
        return (e.pts = pts);
      }
      if (cv.type === 'intcurve-curve') {
        const st = subtype(dataOf(cv));
        const bs = st && st.at >= 0 ? bs3Curve(dataOf(cv), st.at) : null;
        const flat = bs ? sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8) : null;
        if (!flat) return (e.pts = chord('An intcurve-curve without stored B-spline data was replaced by a straight chord between its end points'));
        const sm = [];
        for (let i = 0; i < flat.length; i += 3) sm.push([flat[i], flat[i + 1], flat[i + 2]]);
        if (!Sp || !Ep) return (e.pts = sm);
        const closedCurve = dist3(sm[0], sm[sm.length - 1]) <= tol * 1e3;
        // a vertex belongs to the curve when it lies on the sampled polyline to within a tenth of the local spacing
        const offSeg = (p, a, b) => { const e = V3.sub(b, a), l2 = V3.dot(e, e); const f = l2 > 0 ? Math.max(0, Math.min(1, V3.dot(V3.sub(p, a), e) / l2)) : 0; return dist3(p, lerp3(a, b, f)); };
        const near = (p) => {
          let b = 0, d = Infinity;
          sm.forEach((q, i) => { const x = dist3(p, q); if (x < d) { d = x; b = i; } });
          const sp = Math.max(b > 0 ? dist3(sm[b], sm[b - 1]) : 0, b + 1 < sm.length ? dist3(sm[b], sm[b + 1]) : 0);
          const off = Math.min(b > 0 ? offSeg(p, sm[b - 1], sm[b]) : Infinity, b + 1 < sm.length ? offSeg(p, sm[b], sm[b + 1]) : Infinity);
          return { i: b, ok: off <= 0.1 * sp + tol };
        };
        let pts = null;
        if (same) { if (closedCurve && dist3(Sp, sm[0]) <= tol * 1e3 + defl * 1e-3) pts = sm.slice(); } else if (!closedCurve) {
          const a = near(Sp), b = near(Ep);
          if (a.ok && b.ok && a.i !== b.i) pts = a.i < b.i ? sm.slice(a.i, b.i + 1) : sm.slice(b.i, a.i + 1).reverse();
        }
        if (!pts) return (e.pts = chord('An intcurve-curve whose B-spline does not pass through its edge vertices in an unambiguous way was replaced by a straight chord'));
        if (st.name !== 'exactcur') W.add(`The stored approximating B-spline of a procedural curve ("${st.name}") was used for its edge`);
        pts[0] = Sp; pts[pts.length - 1] = same ? Sp : Ep;
        return (e.pts = pts);
      }
      return (e.pts = chord(`A curve of type "${cv.type}" is not read; its edge was replaced by a straight chord`));
    };

    const usedEdges = new Set();
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
            usedEdges.add(e);
            const ep = edgePts(e);
            const seq = reversedWord(c) ? ep.slice().reverse() : ep;
            for (let i = 0; i + 1 < seq.length; i++) if (!pts.length || pts[pts.length - 1] !== seq[i]) pts.push(seq[i]);
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

    // ---- faces ----
    const weld = new Welder(tol, 4096);
    const T = [], F = [];
    const faceNames = [];
    const faces = root.recs.filter((r) => r.kind === 'face');
    for (const face of faces) {
      const sr = firstOf(face, 'surface');
      const flipFace = reversedWord(face);
      const label = `face ${faceNames.length + 1}`;
      const emit = (pts, idx, flip) => {
        const ids = pts.map((p) => weld.add(p[0], p[1], p[2]));
        let n = 0;
        for (let i = 0; i < idx.length; i += 3) {
          const a = ids[idx[i]], b = ids[idx[i + (flip ? 2 : 1)]], c = ids[idx[i + (flip ? 1 : 2)]];
          if (a === b || b === c || a === c) continue;
          T.push(a, b, c); F.push(faceNames.length); n++;
        }
        if (n) { faceNames.push(label); facesRead++; }
        return n;
      };
      if (!sr) { W.add('A face without a surface record was skipped'); continue; }
      let loops;
      try { loops = loopsOf(face); } catch (e) { W.add(`A face was skipped because its boundary could not be followed (${e.message})`); continue; }

      if (sr.type === 'spline-surface') {
        const st = subtype(dataOf(sr));
        const bs = st && st.name === 'exactsur' && st.at >= 0 ? bs3Surface(dataOf(sr), st.at) : null;
        const corners = [];
        for (const L of loops) for (const p of L) corners.push(p);
        let done = false;
        if (bs && corners.length) {
          // the order of the control net in the record is settled by the face's own boundary: the net must pass by it
          const alt = new Float64Array(bs.P.length), altW = bs.W ? new Float64Array(bs.W.length) : null;
          for (let iu = 0; iu < bs.nu; iu++) for (let iv = 0; iv < bs.nv; iv++) {
            const from = iv + bs.nv * iu, to = iu + bs.nu * iv;
            alt[3 * to] = bs.P[3 * from]; alt[3 * to + 1] = bs.P[3 * from + 1]; alt[3 * to + 2] = bs.P[3 * from + 2];
            if (altW) altW[to] = bs.W[from];
          }
          let best = null;
          for (const [P, Wt] of [[bs.P, bs.W], [alt, altW]]) {
            const m = tessellateSurface(bs.pu, bs.pv, bs.U, bs.V, P, Wt, bs.nu - 1, bs.nv - 1, 6);
            if (!m) continue;
            let miss = 0, edge = 0;
            for (const p of corners) { let d = Infinity; for (let i = 0; i < m.coords.length; i += 3) { const x = Math.hypot(m.coords[i] - p[0], m.coords[i + 1] - p[1], m.coords[i + 2] - p[2]); if (x < d) d = x; } if (d > miss) miss = d; }
            for (let i = 0; i + 3 < m.coords.length; i += 3) { const x = Math.hypot(m.coords[i + 3] - m.coords[i], m.coords[i + 4] - m.coords[i + 1], m.coords[i + 5] - m.coords[i + 2]); if ((i / 3 + 1) % m.nu && x > edge) edge = x; }
            for (let i = 0; i + 3 * m.nu < m.coords.length; i += 3) { const j = i + 3 * m.nu; const x = Math.hypot(m.coords[j] - m.coords[i], m.coords[j + 1] - m.coords[i + 1], m.coords[j + 2] - m.coords[i + 2]); if (x > edge) edge = x; }
            if (miss <= edge + tol && (!best || miss < best.miss - tol)) best = { m, miss };
          }
          if (best) {
            const pts = [];
            for (let i = 0; i < best.m.coords.length; i += 3) pts.push([best.m.coords[i], best.m.coords[i + 1], best.m.coords[i + 2]]);
            if (emit(pts, best.m.indices, st.reversed !== flipFace)) { untrimmed++; tally(hist, 'spline (untrimmed)'); done = true; }
          }
        }
        if (!done) freeForm++;
        continue;
      }

      const surf = makeSurface(sr.type, numsOf(sr), stepFor);
      if (!surf) {
        if (/^(plane|cone|sphere|torus)-surface$/.test(sr.type)) W.add(`A ${sr.type} record could not be read (missing or unsupported values); its face was skipped`);
        else freeForm++;
        continue;
      }
      const flip = surf.inward !== flipFace;
      try {
        let polys, sc;
        if (surf.type === 'plane') {
          sc = surf.scales();
          polys = loops.map((L) => L.map((p) => { const [x, y] = surf.toParam(p); return { x, y, p }; }));
          if (!polys.length) { W.add('A planar face without boundary loops was skipped'); continue; }
        } else {
          let rmax = 0;
          if (surf.radiusAt) { for (const L of loops) for (const p of L) rmax = Math.max(rmax, surf.radiusAt(surf.toParam(p)[1])); if (!(rmax > 0)) { W.add('A conical face without boundary was skipped'); continue; } }
          sc = surf.scales(rmax);
          const un = loops.map((L) => unwrapLoop(surf, L)).filter(Boolean);
          polys = periodicDomain(surf, un, flip, sc);
          if (typeof polys === 'string') { W.add(`A ${surf.type} face was skipped: ${polys}`); continue; }
        }
        if (sc.squeeze) {
          let v0 = Infinity, v1 = -Infinity;
          for (const L of polys) for (const q of L) { if (q.y < v0) v0 = q.y; if (q.y > v1) v1 = q.y; }
          if (v1 > v0) sc.sv = sc.squeeze / (v1 - v0);
        }
        const lifted = liftDomain(surf, polys, sc);
        if (!lifted || !emit(lifted.pts, lifted.indices, flip)) { W.add(`A ${surf.type} face with a degenerate boundary was skipped`); continue; }
        tally(hist, surf.type);
      } catch (e) {
        W.add(`A ${surf.type} face could not be triangulated and was skipped (${e.message})`);
      }
    }

    if (T.length) {
      const P = weld.finish();
      closeTJunctions(P, T, F, diag * 1e-7);
      const closed = isClosed(T, P.length / 3);
      if (closed && volumeOf(P, T) < 0) {
        for (let t = 0; t < T.length; t += 3) { const x = T[t + 1]; T[t + 1] = T[t + 2]; T[t + 2] = x; }
        extra.push(`${bodyName}: the faces enclose a negative volume as stored, so their orientation was reversed to point outwards.`);
      }
      const mirror = tf && tf.mirror;
      const coords = new Float64Array(P.length);
      for (let i = 0; i < P.length; i += 3) { const p = tf ? tf.apply([P[i], P[i + 1], P[i + 2]]) : [P[i], P[i + 1], P[i + 2]]; coords[i] = p[0] * S; coords[i + 1] = p[1] * S; coords[i + 2] = p[2] * S; }
      // triangles grouped by face, in face order
      const order = Array.from({ length: F.length }, (_, i) => i).sort((a, b) => F[a] - F[b] || a - b);
      const indices = new Uint32Array(T.length);
      const groups = [];
      order.forEach((t, k) => {
        indices[3 * k] = T[3 * t]; indices[3 * k + 1] = T[3 * t + (mirror ? 2 : 1)]; indices[3 * k + 2] = T[3 * t + (mirror ? 1 : 2)];
        const g = groups[groups.length - 1];
        if (g && g.face === F[t]) g.count += 3; else groups.push({ name: faceNames[F[t]], start: 3 * k, count: 3, face: F[t] });
      });
      for (const g of groups) delete g.face;
      models.push(surfaceModel(cleanName(bodyName), coords, indices, groups));
      if (closed) closedBodies++;
    }

    // ---- wire edges (edges no face uses) ----
    const wires = root.recs.filter((r) => r.kind === 'edge' && !usedEdges.has(r));
    if (wires.length && anyFace) wireEdges += wires.length;
    else for (const e of wires) {
      const pts = edgePts(e);
      for (let i = 0; i + 1 < pts.length; i++) {
        const a = tf ? tf.apply(pts[i]) : pts[i], b = tf ? tf.apply(pts[i + 1]) : pts[i + 1];
        frameNodes.push([a[0] * S, a[1] * S, a[2] * S], [b[0] * S, b[1] * S, b[2] * S]);
        frameMembers.push(bodyName);
      }
    }
  });

  if (frameMembers.length) {
    let lo = Infinity, hi = -Infinity;
    for (const p of frameNodes) for (const v of p) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const weld = new Welder(Math.max(hi - lo, 1e-12) * 1e-9, frameNodes.length);
    const members = [];
    for (let i = 0; i < frameMembers.length; i++) {
      const a = weld.add(...frameNodes[2 * i]), b = weld.add(...frameNodes[2 * i + 1]);
      if (a !== b) members.push({ a, b, group: frameMembers[i] });
    }
    if (members.length) { models.push({ kind: 'frame', name: cleanName(models.length ? `${name} (wires)` : name), nodes: weld.finish(), members }); tally(hist, 'wire edge', frameMembers.length); }
  }

  const warnings = [];
  if (freeForm) warnings.push(`${freeForm} free-form face${freeForm > 1 ? 's' : ''} could not be read — export those bodies as STEP.`);
  if (untrimmed) warnings.push(`${untrimmed} B-spline surface face${untrimmed > 1 ? 's were' : ' was'} tessellated over the whole surface: trimming by the face boundary was not applied, so the patch may extend beyond the face and does not join its neighbours.`);
  if (wireEdges) warnings.push(`${wireEdges} wire edge${wireEdges > 1 ? 's' : ''} outside any face ${wireEdges > 1 ? 'were' : 'was'} not imported (wire bodies are read only from files without faces).`);
  for (const k in skipped) warnings.push(`${skipped[k]} × ${k} not recognised and ignored.`);
  for (const w of W.list()) warnings.push(`${w}.`);
  warnings.push(...extra);
  if (!models.length) throw new Error(`${fname}: no geometry could be read from the SAT file.${warnings.length ? ` ${warnings.join(' ')}` : ''}`);
  if (freeForm) skipped['free-form face'] = freeForm;
  return finish(models, warnings, {
    format: 'acis', units, crs: null, elementTypes: hist, skipped,
    acisVersion: header.version, product: header.product, acisRelease: header.acis, saved: header.date,
    bodies: bodies.length, closedBodies, faces: facesRead, resabs: header.resabs,
    tessellation: { chord: chordM, angularStepDeg: opts.angularStep ?? 15 },
  });
}
