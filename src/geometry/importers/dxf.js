// AutoCAD DXF, ASCII form (R12 – 2018).
//
// Entities read: LINE, LWPOLYLINE and 2-D/3-D POLYLINE (bulge arcs included), polyface and polygon-mesh
// POLYLINE, 3DFACE, SOLID/TRACE, CIRCLE, ARC, ELLIPSE, SPLINE (NURBS evaluation; fit points when no control
// points), POINT and INSERT (nested blocks, scale, rotation, arrays). Extrusion directions (group 210) are
// applied with the arbitrary-axis algorithm. $INSUNITS gives the length unit; layers are kept as groups.
//
// Mapping to model kinds:
//   3DFACE / polyface / mesh / SOLID   → surface, one group per layer
//   closed planar loops                → section (loops grouped by layer in `layers`)
//   LINE                               → frame (end points welded; member.group = layer) — unless the
//                                        drawing is flat and every line/arc closes into loops, in which case
//                                        the loops are section outlines (opts.dxfLines = 'frame' | 'section'
//                                        overrides the choice)
//   open polylines, arcs, splines      → polyline
//   POINT                              → points
//
//   3DSOLID / REGION / BODY / SURFACE  → surface, one model per ACIS body with one group per B-rep face
//
// Solids are stored as ACIS data: as SAT text in groups 1 and 3 of the entity (R13 – R2010, each character
// replaced by a fixed substitution, see decodeAcisText), or as SAB bytes in the ACDSDATA section keyed by the
// entity handle (R2013 and later). Both are decoded here and tessellated by acis.js / sab.js; the coordinates
// are drawing units in the space of the entity, so block transforms and $INSUNITS apply as for everything else.
//
// Not read: MESH, HATCH, text, dimensions, leaders, images, viewports, paper space layout. These are counted
// in `meta.skipped`, as are solids whose ACIS data are missing or cannot be decoded. Binary DXF is decoded by
// dxfbinary.js and enters at dxfModels().

import { lineCursor, baseName, Grow, Welder, finish, surfaceModel, pointsModel, M4, V3, tally, cleanName, groupsFromBuckets, unitAdvice, closeLoop } from './common.js';
import { resolveUnits, suggestUnits } from '../units.js';
import { sampleCurve, arcPoints } from './nurbs.js';
import { readSat, acisBodies } from './acis.js';
import { readSab, isSab } from './sab.js';

const INSUNITS = { 1: 'in', 2: 'ft', 3: 'mi', 4: 'mm', 5: 'cm', 6: 'm', 7: 'km', 9: 'mil', 10: 'yd', 13: 'um', 14: 'dm' };
const D2R = Math.PI / 180;
const ACIS_ENTITY = new Set(['3DSOLID', 'REGION', 'BODY', 'SURFACE', 'EXTRUDEDSURFACE', 'REVOLVEDSURFACE', 'LOFTEDSURFACE', 'SWEPTSURFACE', 'PLANESURFACE', 'NURBSSURFACE']);

/**
 * Undo the character substitution AutoCAD applies to the SAT text it stores in groups 1 and 3: a space stays a
 * space, every other character c becomes 159 − c (which exchanges "_" and "@", digits and lower-case letters,
 * and so on), and the space written after a "^" (the substitute of "A") is dropped.
 */
export function decodeAcisText(line) {
  let out = '';
  for (let i = 0; i < line.length; i++) {
    const c = line.charCodeAt(i);
    if (c === 32) out += ' ';
    else if (c < 33 || c > 126) out += line[i];
    else { out += String.fromCharCode(159 - c); if (c === 0x5e && line.charCodeAt(i + 1) === 32) i++; }
  }
  return out;
}

const hexBytes = (hex) => { const n = hex.length >> 1, b = new Uint8Array(n); for (let i = 0; i < n; i++) b[i] = parseInt(hex.substr(2 * i, 2), 16); return b; };

/** Object-coordinate-system frame for an extrusion direction (arbitrary axis algorithm, DXF reference). */
export function ocsMatrix(nx, ny, nz) {
  const n = V3.norm([nx, ny, nz]);
  if (n[0] === 0 && n[1] === 0 && n[2] > 0) return null;
  const ax = V3.norm(Math.abs(n[0]) < 1 / 64 && Math.abs(n[1]) < 1 / 64 ? V3.cross([0, 1, 0], n) : V3.cross([0, 0, 1], n));
  return M4.fromAxes([0, 0, 0], ax, V3.cross(n, ax), n);
}

/** Points of a bulge arc from (x1, y1) to (x2, y2), excluding the start point. bulge = tan(θ / 4). */
export function bulgeArc(x1, y1, x2, y2, bulge, out) {
  const dx = x2 - x1, dy = y2 - y1, c = Math.hypot(dx, dy);
  if (!(c > 0) || !bulge) { out.push(x2, y2); return; }
  const d = c / 2, s = bulge * d, r = (d * d + s * s) / (2 * s);
  // the centre lies to the left of the chord for a counter-clockwise (positive) bulge
  const cx = (x1 + x2) / 2 - (dy / c) * (r - s), cy = (y1 + y2) / 2 + (dx / c) * (r - s);
  const a0 = Math.atan2(y1 - cy, x1 - cx), sweep = 4 * Math.atan(bulge);
  const n = Math.max(2, Math.ceil(Math.abs(sweep) / (5 * D2R)));
  const R = Math.abs(r);
  for (let k = 1; k < n; k++) { const t = a0 + (sweep * k) / n; out.push(cx + R * Math.cos(t), cy + R * Math.sin(t)); }
  out.push(x2, y2);
}

function readPairs(input) {
  const cur = lineCursor(input);
  const sections = { HEADER: [], TABLES: [], BLOCKS: [], ENTITIES: [], ACDSDATA: [] };
  let target = null, n = 0, bad = 0;
  for (;;) {
    const a = cur.next();
    if (a === null) break;
    const b = cur.next();
    if (b === null) break;
    const code = parseInt(a, 10);
    if (code !== code) { bad++; if (bad > 20 && n < 5) return null; continue; }
    n++;
    // the text of groups 1 and 3 is kept as written: embedded ACIS lines are split across them at any character
    const v = (code === 1 || code === 3) && (target === sections.ENTITIES || target === sections.BLOCKS) ? b : b.trim();
    if (code === 0 && v === 'SECTION') { const c2 = cur.next(), nm = cur.next(); if (c2 === null || nm === null) break; target = sections[nm.trim().toUpperCase()] || null; continue; }
    if (code === 0 && v === 'ENDSEC') { target = null; continue; }
    if (code === 0 && v === 'EOF') break;
    if (target) target.push(code, v);
  }
  return n ? sections : null;
}

/** Split a flat [code, value, …] list into records that start at group 0. */
function records(flat) {
  const out = [];
  let cur = null;
  for (let i = 0; i < flat.length; i += 2) {
    if (flat[i] === 0) { cur = { type: String(flat[i + 1]).trim().toUpperCase(), g: [] }; out.push(cur); } else if (cur) cur.g.push(flat[i], flat[i + 1]);
  }
  return out;
}

const first = (g, code, d) => { for (let i = 0; i < g.length; i += 2) if (g[i] === code) return code === 1 || code === 3 ? String(g[i + 1]).trim() : g[i + 1]; return d; };
const num = (g, code, d = 0) => { const v = first(g, code); if (v === undefined) return d; const x = Number(v); return x === x ? x : d; };
const all = (g, code) => { const o = []; for (let i = 0; i < g.length; i += 2) if (g[i] === code) o.push(Number(g[i + 1])); return o; };

/** AutoCAD DXF (ASCII). */
export function parseDxf(input, opts = {}) {
  const fname = input.name || 'DXF';
  if (input.bytes && input.bytes.length > 22 && String.fromCharCode(...input.bytes.subarray(0, 18)) === 'AutoCAD Binary DXF') throw new Error(`${fname}: binary DXF is not read. Save the drawing as ASCII DXF (SAVEAS → "AutoCAD DXF", not "Binary DXF").`);
  const sec = readPairs(input);
  if (!sec || (!sec.ENTITIES.length && !sec.BLOCKS.length)) throw new Error(`${fname}: not an ASCII DXF drawing with entities (expected group-code / value line pairs with "0 / SECTION … 2 / ENTITIES").`);
  return dxfModels(sec, input, opts);
}

/**
 * Build the models from the group pairs of a drawing, whatever encoding they were stored in.
 * @param {{ HEADER: Array, TABLES: Array, BLOCKS: Array, ENTITIES: Array, ACDSDATA?: Array }} sec flat
 *   [code, value, …] lists per section, codes as numbers and values as trimmed strings (exactly what the ASCII
 *   form holds; groups 1 and 3 of entities untrimmed)
 */
export function dxfModels(sec, input, opts = {}) {
  const fname = input.name || 'DXF';
  const name = baseName(input.name);
  const warnings = [];
  const hist = {}, skipped = {};

  const header = {};
  for (let i = 0; i < sec.HEADER.length; i += 2) if (sec.HEADER[i] === 9 && i + 3 < sec.HEADER.length) header[String(sec.HEADER[i + 1]).trim()] = String(sec.HEADER[i + 3]).trim();
  const insunits = Number(header.$INSUNITS);
  const layers = [];
  for (const r of records(sec.TABLES)) if (r.type === 'LAYER') layers.push(first(r.g, 2, ''));

  // fold VERTEX records into their POLYLINE
  const fold = (recs) => {
    const out = [];
    let pl = null;
    for (const r of recs) {
      if (r.type === 'POLYLINE') { pl = r; pl.verts = []; out.push(pl); } else if (r.type === 'VERTEX' && pl) pl.verts.push(r.g);
      else if (r.type === 'SEQEND') pl = null;
      else { pl = null; out.push(r); }
    }
    return out;
  };
  const blocks = new Map();
  {
    let cur = null;
    for (const r of fold(records(sec.BLOCKS))) {
      if (r.type === 'BLOCK') { cur = { name: first(r.g, 2, ''), base: [num(r.g, 10), num(r.g, 20), num(r.g, 30)], ents: [] }; blocks.set(cur.name.toUpperCase(), cur); } else if (r.type === 'ENDBLK') cur = null;
      else if (cur) cur.ents.push(r);
    }
  }
  const entities = fold(records(sec.ENTITIES)).filter((r) => num(r.g, 67, 0) !== 1);

  const units = resolveUnits(opts, INSUNITS[insunits] || null, 'm');
  const s = units.scale;

  // SAB data of the solids of R2013 and later drawings: ACDSRECORD entries with the owner handle and the bytes
  const sabByHandle = new Map();
  for (const r of records(sec.ACDSDATA || [])) {
    if (r.type !== 'ACDSRECORD') continue;
    let part = '', handle = '', hex = [];
    for (let i = 0; i < r.g.length; i += 2) {
      if (r.g[i] === 2) part = r.g[i + 1];
      else if (r.g[i] === 320 && part === 'AcDbDs::ID') handle = String(r.g[i + 1]).toUpperCase();
      else if (r.g[i] === 310 && part === 'ASM_Data') hex.push(r.g[i + 1]);
    }
    if (handle && hex.length) sabByHandle.set(handle, hex);
  }

  // primitives in world coordinates
  const segs = [], polys = [], faces = [], points = [], solids = [];
  // `solidBodies`: bodies whose topology is that of a closed solid (every edge shared by two faces);
  // `closedBodies`: bodies whose triangles came out watertight and consistently oriented
  const acis = { entities: 0, text: 0, binary: 0, bodies: 0, solidBodies: 0, closedBodies: 0, faces: 0, facesRead: 0 };
  const acisNotes = new Map();
  if (opts.deflection !== undefined && !(opts.deflection > 0)) throw new RangeError('deflection: must be a positive length in metres.');
  if (opts.angularStep !== undefined && !(opts.angularStep > 0 && opts.angularStep <= 90)) throw new RangeError('angularStep: must be within 0…90 degrees.');
  /** Decode and tessellate the ACIS data of one entity; false (with the reason tallied) when that is not possible. */
  const emitAcis = (r, M, layer) => {
    const g = r.g;
    const lines = [];
    for (let i = 0; i < g.length; i += 2) { if (g[i] === 1) lines.push(String(g[i + 1])); else if (g[i] === 3 && lines.length) lines[lines.length - 1] += String(g[i + 1]); }
    let sat = null, why = '';
    try {
      if (lines.length) {
        const plain = /^\s*\d{3,5}\s+\d+\s+\d+\s+\d+\s*$/.test(lines[0]);
        sat = readSat((plain ? lines : lines.map(decodeAcisText)).join('\n'), `${fname}: ${r.type}`);
        acis.text++;
      } else {
        const hex = sabByHandle.get(String(first(g, 5, '')).toUpperCase());
        if (!hex) why = 'no ACIS data stored with the entity';
        else {
          const bytes = hexBytes(hex.join(''));
          if (!isSab(bytes)) why = 'the data stored for the entity are not ACIS binary data';
          else { sat = readSab(bytes, `${fname}: ${r.type}`); acis.binary++; }
        }
      }
    } catch (e) { why = e.message.replace(/^.*?: /, '').replace(/\.$/, ''); }
    let res = null;
    if (sat) { try { res = acisBodies(sat, { deflection: opts.deflection !== undefined ? opts.deflection / s : undefined, angularStep: opts.angularStep, name: `${layer} ${r.type}` }); } catch (e) { why = e.message; } }
    if (!res || (!res.bodies.length && !res.wires.length)) {
      tally(skipped, `${r.type} (embedded ACIS geometry: ${why || (res && res.warnings[0]) || 'no faces could be read'})`);
      return false;
    }
    acis.entities++;
    acis.bodies += res.bodies.length; acis.faces += res.stats.faces; acis.facesRead += res.stats.facesRead;
    for (const b of res.bodies) { if (!b.freeEdges && !b.nonManifoldEdges) acis.solidBodies++; if (b.closed && b.oriented) acis.closedBodies++; }
    for (const w of res.warnings) acisNotes.set(w, (acisNotes.get(w) || 0) + 1);
    const mirror = M ? M[0] * (M[5] * M[10] - M[6] * M[9]) - M[4] * (M[1] * M[10] - M[2] * M[9]) + M[8] * (M[1] * M[6] - M[2] * M[5]) < 0 : false;
    for (const b of res.bodies) {
      const c = b.coords;
      if (M) for (let i = 0; i < c.length; i += 3) { const p = M4.apply(M, c[i], c[i + 1], c[i + 2]); c[i] = p[0]; c[i + 1] = p[1]; c[i + 2] = p[2]; }
      if (mirror) for (let i = 0; i < b.indices.length; i += 3) { const x = b.indices[i + 1]; b.indices[i + 1] = b.indices[i + 2]; b.indices[i + 2] = x; }
      solids.push({ c, indices: b.indices, groups: b.groups, layer, type: r.type });
    }
    for (const w of res.wires) segs.push({ a: P(M, w.a[0], w.a[1], w.a[2]), b: P(M, w.b[0], w.b[1], w.b[2]), layer });
    tally(hist, r.type);
    return true;
  };
  const P = (M, x, y, z) => (M ? M4.apply(M, x, y, z) : [x, y, z]);
  const emitPoly = (M, flat2or3, dim, z, closed, layer, type) => {
    const n = flat2or3.length / dim;
    if (n < 2) return;
    const c = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) { const p = P(M, flat2or3[dim * i], flat2or3[dim * i + 1], dim === 3 ? flat2or3[dim * i + 2] : z); c[3 * i] = p[0]; c[3 * i + 1] = p[1]; c[3 * i + 2] = p[2]; }
    polys.push({ c, closed, layer, type });
  };
  const frame = (M, g) => { const o = ocsMatrix(num(g, 210, 0), num(g, 220, 0), num(g, 230, 1)); return o ? (M ? M4.mul(M, o) : o) : M; };

  const emit = (r, M, layerOverride, depth) => {
    const g = r.g;
    let layer = first(g, 8, '0');
    if (layer === '0' && layerOverride) layer = layerOverride;
    switch (r.type) {
      case 'LINE': segs.push({ a: P(M, num(g, 10), num(g, 20), num(g, 30)), b: P(M, num(g, 11), num(g, 21), num(g, 31)), layer }); tally(hist, 'LINE'); break;
      case 'POINT': points.push(P(M, num(g, 10), num(g, 20), num(g, 30))); tally(hist, 'POINT'); break;
      case '3DFACE': case 'SOLID': case 'TRACE': {
        const T = r.type === '3DFACE' ? M : frame(M, g);
        const v = [0, 1, 2, 3].map((k) => P(T, num(g, 10 + k), num(g, 20 + k), num(g, 30 + k)));
        // SOLID / TRACE list their corners in a Z pattern: 1, 2, 4, 3
        faces.push({ v: r.type === '3DFACE' ? v : [v[0], v[1], v[3], v[2]], layer });
        tally(hist, r.type);
        break;
      }
      case 'CIRCLE': case 'ARC': {
        const T = frame(M, g);
        const rad = num(g, 40), cx = num(g, 10), cy = num(g, 20), cz = num(g, 30);
        const full = r.type === 'CIRCLE';
        const a0 = full ? 0 : num(g, 50) * D2R;
        let sw = full ? 2 * Math.PI : num(g, 51) * D2R - a0;
        if (!full) { while (sw <= 0) sw += 2 * Math.PI; }
        let pts = arcPoints([cx, cy, cz], [1, 0, 0], [0, 1, 0], rad, a0, sw);
        if (full) pts = pts.subarray(0, pts.length - 3);
        emitPoly(T, pts, 3, 0, full, layer, r.type);
        tally(hist, r.type);
        break;
      }
      case 'ELLIPSE': {
        // centre and major axis are in world coordinates; the minor axis is N × major scaled by the ratio
        const c = [num(g, 10), num(g, 20), num(g, 30)], mj = [num(g, 11), num(g, 21), num(g, 31)], ratio = num(g, 40, 1);
        const nrm = V3.norm([num(g, 210, 0), num(g, 220, 0), num(g, 230, 1)]);
        const mn = V3.scale(V3.cross(nrm, mj), ratio);
        const t0 = num(g, 41, 0), t1 = num(g, 42, 2 * Math.PI);
        let sw = t1 - t0;
        while (sw <= 0) sw += 2 * Math.PI;
        const full = Math.abs(sw - 2 * Math.PI) < 1e-9;
        let pts = arcPoints(c, mj, mn, 1, t0, sw);
        if (full) pts = pts.subarray(0, pts.length - 3);
        emitPoly(M, pts, 3, 0, full, layer, 'ELLIPSE');
        tally(hist, 'ELLIPSE');
        break;
      }
      case 'LWPOLYLINE': {
        const T = frame(M, g);
        const closed = (num(g, 70, 0) & 1) === 1, elev = num(g, 38, 0);
        const xs = [], ys = [], bl = [];
        for (let i = 0; i < g.length; i += 2) { if (g[i] === 10) { xs.push(Number(g[i + 1])); bl.push(0); } else if (g[i] === 20) ys.push(Number(g[i + 1])); else if (g[i] === 42 && bl.length) bl[bl.length - 1] = Number(g[i + 1]); }
        const n = Math.min(xs.length, ys.length);
        if (n < 2) { tally(skipped, 'LWPOLYLINE with fewer than two vertices'); break; }
        const flat = [xs[0], ys[0]];
        for (let i = 1; i < n; i++) bulgeArc(xs[i - 1], ys[i - 1], xs[i], ys[i], bl[i - 1], flat);
        if (closed && bl[n - 1]) { bulgeArc(xs[n - 1], ys[n - 1], xs[0], ys[0], bl[n - 1], flat); flat.length -= 2; }
        emitPoly(T, flat, 2, elev, closed, layer, 'LWPOLYLINE');
        tally(hist, 'LWPOLYLINE');
        break;
      }
      case 'POLYLINE': {
        const flags = num(g, 70, 0);
        const vs = r.verts || [];
        if (flags & 64) {
          // polyface mesh: vertices with flag 128 are coordinates, records with 71–74 are faces (1-based, sign = edge visibility)
          const vv = [], ff = [];
          for (const v of vs) { const vf = num(v, 70, 0); if ((vf & 128) && (vf & 64)) vv.push(P(M, num(v, 10), num(v, 20), num(v, 30))); else if (vf & 128) ff.push([71, 72, 73, 74].map((c) => Math.abs(num(v, c, 0))).filter((k) => k > 0)); }
          for (const f of ff) { if (f.length >= 3 && f.every((k) => k <= vv.length)) faces.push({ v: f.map((k) => vv[k - 1]), layer: first(g, 8, layer) === '0' ? layer : first(g, 8, layer) }); else tally(skipped, 'polyface record with invalid vertex index'); }
          tally(hist, 'POLYFACE');
        } else if (flags & 16) {
          const m = num(g, 71, 0), n = num(g, 72, 0);
          const vv = vs.map((v) => P(M, num(v, 10), num(v, 20), num(v, 30)));
          if (m * n !== vv.length || m < 2 || n < 2) { tally(skipped, 'polygon mesh with inconsistent vertex count'); break; }
          const cm = flags & 1 ? m : m - 1, cn = flags & 32 ? n : n - 1;
          for (let i = 0; i < cm; i++) for (let j = 0; j < cn; j++) faces.push({ v: [vv[i * n + j], vv[((i + 1) % m) * n + j], vv[((i + 1) % m) * n + ((j + 1) % n)], vv[i * n + ((j + 1) % n)]], layer });
          tally(hist, 'POLYGON MESH');
        } else if (flags & 8) {
          const flat = [];
          for (const v of vs) flat.push(num(v, 10), num(v, 20), num(v, 30));
          emitPoly(M, flat, 3, 0, (flags & 1) === 1, layer, 'POLYLINE 3D');
          tally(hist, 'POLYLINE');
        } else {
          const T = frame(M, g), elev = num(g, 30, 0), closed = (flags & 1) === 1;
          // spline-fit / curve-fit frames carry extra control vertices (flag 16): skip those
          const use = vs.filter((v) => !(num(v, 70, 0) & 16));
          if (use.length < 2) { tally(skipped, 'POLYLINE with fewer than two vertices'); break; }
          const flat = [num(use[0], 10), num(use[0], 20)];
          for (let i = 1; i < use.length; i++) bulgeArc(num(use[i - 1], 10), num(use[i - 1], 20), num(use[i], 10), num(use[i], 20), num(use[i - 1], 42, 0), flat);
          const lb = num(use[use.length - 1], 42, 0);
          if (closed && lb) { bulgeArc(num(use[use.length - 1], 10), num(use[use.length - 1], 20), num(use[0], 10), num(use[0], 20), lb, flat); flat.length -= 2; }
          emitPoly(T, flat, 2, elev, closed, layer, 'POLYLINE');
          tally(hist, 'POLYLINE');
        }
        break;
      }
      case 'SPLINE': {
        const flags = num(g, 70, 0), deg = num(g, 71, 3);
        const knots = all(g, 40), w = all(g, 41);
        const cx = all(g, 10), cy = all(g, 20), cz = all(g, 30);
        let pts = null;
        if (cx.length > deg) {
          const ctrl = new Float64Array(cx.length * 3);
          for (let i = 0; i < cx.length; i++) { ctrl[3 * i] = cx[i]; ctrl[3 * i + 1] = cy[i] || 0; ctrl[3 * i + 2] = cz[i] || 0; }
          pts = sampleCurve(deg, knots, ctrl, w.length === cx.length ? w : null, 8);
        }
        if (!pts) {
          const fx = all(g, 11), fy = all(g, 21), fz = all(g, 31);
          if (fx.length >= 2) { pts = new Float64Array(fx.length * 3); for (let i = 0; i < fx.length; i++) { pts[3 * i] = fx[i]; pts[3 * i + 1] = fy[i] || 0; pts[3 * i + 2] = fz[i] || 0; } tally(skipped, 'SPLINE drawn through its fit points'); }
        }
        if (!pts) { tally(skipped, 'SPLINE with an invalid definition'); break; }
        const n = pts.length;
        const closedGeom = Math.hypot(pts[0] - pts[n - 3], pts[1] - pts[n - 2], pts[2] - pts[n - 1]) < 1e-9 * (1 + Math.abs(pts[0]) + Math.abs(pts[1]));
        const closed = (flags & 1) === 1 || closedGeom;
        emitPoly(M, closedGeom ? pts.subarray(0, n - 3) : pts, 3, 0, closed, layer, 'SPLINE');
        tally(hist, 'SPLINE');
        break;
      }
      case 'INSERT': {
        const blk = blocks.get(String(first(g, 2, '')).toUpperCase());
        if (!blk) { tally(skipped, 'INSERT of an undefined block'); break; }
        if (depth > 16) { tally(skipped, 'INSERT nested deeper than 16 levels'); break; }
        const T0 = frame(M, g) || M4.identity();
        const cols = Math.max(1, num(g, 70, 1)), rows = Math.max(1, num(g, 71, 1));
        const R = M4.rotation(0, 0, 1, num(g, 50, 0) * D2R), S = M4.scaling(num(g, 41, 1), num(g, 42, 1), num(g, 43, 1));
        const B = M4.translation(-blk.base[0], -blk.base[1], -blk.base[2]);
        for (let i = 0; i < cols; i++) for (let j = 0; j < rows; j++) {
          const T = M4.mul(T0, M4.mul(M4.translation(num(g, 10), num(g, 20), num(g, 30)), M4.mul(R, M4.mul(M4.translation(i * num(g, 44, 0), j * num(g, 45, 0), 0), M4.mul(S, B)))));
          for (const e of blk.ents) emit(e, T, layer, depth + 1);
        }
        tally(hist, 'INSERT');
        break;
      }
      default: if (ACIS_ENTITY.has(r.type)) emitAcis(r, M, layer); else tally(skipped, r.type);
    }
  };
  for (const r of entities) emit(r, null, null, 0);
  const lost = Object.keys(skipped).filter((k) => /ACIS/.test(k));
  if (lost.length) warnings.push(`${lost.reduce((n, k) => n + skipped[k], 0)} solid / surface entit${lost.reduce((n, k) => n + skipped[k], 0) > 1 ? 'ies' : 'y'} stored as ACIS data could not be read (${lost.map((k) => `${skipped[k]} × ${k}`).join('; ')}). Export them as STEP or IGES, or mesh them (3DFACE / polyface), to bring them in.`);
  if (acis.entities) {
    for (const [w, n] of acisNotes) warnings.push(`ACIS solids: ${w}${n > 1 ? ` (in ${n} entities)` : ''}`);
  }
  for (const k in skipped) if (!/ACIS/.test(k)) warnings.push(`Skipped ${skipped[k]} × ${k}.`);

  // extent, for tolerances
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  const see = (x, y, z) => { if (x < lo[0]) lo[0] = x; if (y < lo[1]) lo[1] = y; if (z < lo[2]) lo[2] = z; if (x > hi[0]) hi[0] = x; if (y > hi[1]) hi[1] = y; if (z > hi[2]) hi[2] = z; };
  for (const sg of segs) { see(...sg.a); see(...sg.b); }
  for (const p of polys) for (let i = 0; i < p.c.length; i += 3) see(p.c[i], p.c[i + 1], p.c[i + 2]);
  for (const f of faces) for (const v of f.v) see(...v);
  for (const p of points) see(...p);
  for (const q of solids) for (let i = 0; i < q.c.length; i += 3) see(q.c[i], q.c[i + 1], q.c[i + 2]);
  if (lo[0] === Infinity) throw new Error(`${fname}: the drawing contains no geometry this reader keeps${Object.keys(skipped).length ? ` (found only: ${Object.keys(skipped).join(', ')})` : ''}.`);
  const ext = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-12);
  const tol = ext * 1e-7;
  const flat = hi[2] - lo[2] <= tol;
  const models = [];

  if (faces.length) {
    const weld = new Welder(tol, faces.length * 3);
    const buckets = new Map();
    for (const f of faces) {
      let g = buckets.get(f.layer);
      if (!g) buckets.set(f.layer, (g = new Grow(Uint32Array, 3072)));
      const ids = f.v.map((v) => weld.add(v[0] * s, v[1] * s, v[2] * s)).filter((k, i, a) => a.indexOf(k) === i);
      if (ids.length === 3) { g.push(ids[0]); g.push(ids[1]); g.push(ids[2]); } else if (ids.length === 4) { g.push(ids[0]); g.push(ids[1]); g.push(ids[2]); g.push(ids[0]); g.push(ids[2]); g.push(ids[3]); }
    }
    const { indices, groups } = groupsFromBuckets(buckets);
    if (indices.length) models.push(surfaceModel(name, weld.finish(), indices, groups));
  }

  if (solids.length) {
    const counts = {};
    for (const q of solids) {
      const key = `${q.layer} ${q.type}`;
      counts[key] = (counts[key] || 0) + 1;
      const c = q.c;
      if (s !== 1) for (let i = 0; i < c.length; i++) c[i] *= s;
      models.push(surfaceModel(cleanName(`${key} ${counts[key]}`), c, q.indices, q.groups, { layer: q.layer, entity: q.type }));
    }
  }

  // closed planar loops → section; optionally chain loose LINE/ARC pieces into loops first
  const loops = [], loopLayer = [];
  const open = [];
  const isPlanar = (c) => { for (let i = 5; i < c.length; i += 3) if (Math.abs(c[i] - c[2]) > tol) return false; return true; };
  for (const p of polys) {
    if (p.closed && isPlanar(p.c) && p.c.length >= 9) { const l = new Float64Array((p.c.length / 3) * 2); for (let i = 0; i < p.c.length / 3; i++) { l[2 * i] = p.c[3 * i] * s; l[2 * i + 1] = p.c[3 * i + 1] * s; } const cl = closeLoop(l); if (cl) { loops.push(cl); loopLayer.push(p.layer); continue; } }
    open.push(p);
  }
  const mode = opts.dxfLines || 'auto';
  if (!['auto', 'frame', 'section'].includes(mode)) throw new RangeError(`dxfLines: "${mode}" is not one of auto, frame, section.`);
  let frameSegs = segs, rest = open, joined = 0;
  if (flat && mode !== 'frame' && (segs.length || open.some((p) => !p.closed))) {
    const pieces = segs.map((sg) => ({ c: Float64Array.from([...sg.a, ...sg.b]), layer: sg.layer, seg: true })).concat(open.filter((p) => !p.closed).map((p) => ({ c: p.c, layer: p.layer, src: p })));
    const weld = new Welder(tol * 10, pieces.length * 2);
    const ends = pieces.map((pc) => [weld.add(pc.c[0], pc.c[1], 0), weld.add(pc.c[pc.c.length - 3], pc.c[pc.c.length - 2], 0)]);
    const at = new Map();
    ends.forEach((e, i) => { for (const k of e) { let a = at.get(k); if (!a) at.set(k, (a = [])); a.push(i); } });
    const used = new Uint8Array(pieces.length);
    const found = [];
    for (let i = 0; i < pieces.length; i++) {
      if (used[i] || ends[i][0] === ends[i][1]) continue;
      const chain = [[i, false]];
      const seen = new Set([i]);
      let node = ends[i][1], ok = false;
      for (let guard = 0; guard <= pieces.length; guard++) {
        if (node === ends[i][0]) { ok = true; break; }
        const inc = at.get(node);
        if (!inc || inc.length !== 2) break;
        const nx = inc[0] === chain[chain.length - 1][0] ? inc[1] : inc[0];
        if (seen.has(nx) || used[nx]) break;
        seen.add(nx);
        const rev = ends[nx][0] !== node;
        chain.push([nx, rev]);
        node = rev ? ends[nx][0] : ends[nx][1];
      }
      if (ok && chain.length >= 2) { for (const [k] of chain) used[k] = 1; found.push(chain); }
    }
    const allClose = found.length > 0 && used.every((u) => u === 1);
    if (found.length && (mode === 'section' || allClose)) {
      for (const chain of found) {
        const l = [];
        for (const [k, rev] of chain) { const c = pieces[k].c, n = c.length / 3; for (let q = 0; q < n - 1; q++) { const i = rev ? n - 1 - q : q; l.push(c[3 * i] * s, c[3 * i + 1] * s); } }
        const cl = closeLoop(l);
        if (cl) { loops.push(cl); loopLayer.push(pieces[chain[0][0]].layer); joined++; }
      }
      frameSegs = segs.filter((_, i) => !used[i]);
      const gone = new Set(pieces.filter((pc, i) => used[i] && pc.src).map((pc) => pc.src));
      rest = open.filter((p) => !gone.has(p));
    }
  }
  if (loops.length) {
    const byLayer = new Map();
    loopLayer.forEach((l, i) => { let a = byLayer.get(l); if (!a) byLayer.set(l, (a = [])); a.push(i); });
    models.push({ kind: 'section', name: models.length || frameSegs.length || rest.length ? `${name} (profiles)` : name, loops, layers: [...byLayer].map(([nm, idx]) => ({ name: nm, loops: idx })) });
  }
  if (frameSegs.length) {
    const weld = new Welder(tol * 10, frameSegs.length * 2);
    const members = [];
    let zero = 0;
    for (const sg of frameSegs) { const a = weld.add(sg.a[0] * s, sg.a[1] * s, sg.a[2] * s), b = weld.add(sg.b[0] * s, sg.b[1] * s, sg.b[2] * s); if (a === b) { zero++; continue; } members.push({ a, b, group: sg.layer }); }
    if (zero) { skipped['zero-length LINE'] = zero; warnings.push(`${zero} zero-length LINE(s) were dropped.`); }
    if (members.length) models.push({ kind: 'frame', name: models.length ? `${name} (lines)` : name, nodes: weld.finish(), members });
  }
  const counts = {};
  for (const p of rest) {
    const c = Float64Array.from(p.c, (v) => v * s);
    counts[p.layer] = (counts[p.layer] || 0) + 1;
    const m = { kind: 'polyline', name: cleanName(`${p.layer} ${counts[p.layer]}`), points: c, layer: p.layer, entity: p.type };
    if (p.closed) m.closed = true;
    models.push(m);
  }
  if (points.length) { const c = new Float64Array(points.length * 3); points.forEach((p, i) => { c[3 * i] = p[0] * s; c[3 * i + 1] = p[1] * s; c[3 * i + 2] = p[2] * s; }); models.push(pointsModel(`${name} (points)`, c)); }
  if (!models.length) throw new Error(`${fname}: no usable geometry was found in the drawing.`);
  if (joined) warnings.push(`${joined} closed outline(s) were assembled from separate LINE / ARC / open polyline entities and returned as section loops (set dxfLines = 'frame' to keep lines as members).`);
  const r = finish(models, warnings, {
    format: 'dxf', units, crs: null, elementTypes: hist, skipped, layers: layers.length ? layers : [...new Set([...segs.map((x) => x.layer), ...polys.map((x) => x.layer), ...faces.map((x) => x.layer)])],
    blocks: blocks.size, insunits: Number.isFinite(insunits) ? insunits : null, acadVersion: header.$ACADVER || null, planar: flat,
    ...(acis.entities || lost.length ? { acis: { ...acis, encoding: acis.text && acis.binary ? 'SAT text and SAB binary' : acis.binary ? 'SAB binary (ACDSDATA)' : 'SAT text (groups 1 / 3)' } } : {}),
  });
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}
