// Triangulated-surface and scan formats: STL (ASCII and binary), Wavefront OBJ, Stanford PLY (ASCII and both
// binary byte orders) and Geomview OFF.
//
// None of these formats records a length unit. Metres are assumed unless `opts.units` says otherwise, and the
// assumption is reported in `meta.units` together with a plausibility hint. None records an up axis either:
// z up is assumed (the engineering convention); pass `opts.upAxis = 'y'` for y-up models.

import {
  inputBytes, lineCursor, baseName, Grow, Welder, finish, surfaceModel, pointsModel, bboxOf, unitAdvice,
  cleanName, groupsFromBuckets, pushPolygon, decodeText,
} from './common.js';
import { resolveUnits, suggestUnits, scaleArray, yUpToZUp } from '../units.js';

function orient(coords, opts, scale) {
  scaleArray(coords, scale);
  if (opts.upAxis === 'y') yUpToZUp(coords);
  else if (opts.upAxis && opts.upAxis !== 'z') throw new RangeError(`upAxis: "${opts.upAxis}" is not one of z, y.`);
  return coords;
}

function upAxisMeta(opts) {
  return { source: opts.upAxis === 'y' ? 'y' : 'z', basis: opts.upAxis ? 'option' : 'assumed' };
}

// ---------------------------------------------------------------------------------------------------------
// STL
// ---------------------------------------------------------------------------------------------------------

/**
 * STL, ASCII or binary. Coincident facet corners are welded into shared vertices (bit-identical coordinates),
 * so the result is an indexed mesh on which adjacency and watertightness can be checked.
 * Each `solid … endsolid` block of an ASCII file becomes a group.
 */
export function parseStl(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  const warnings = [];
  if (bytes.length < 15) throw new Error(`${input.name || 'STL'}: file is too short to be an STL (expected "solid …" text or an 80-byte header followed by a triangle count).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const nBin = bytes.length >= 84 ? dv.getUint32(80, true) : -1;
  const startsSolid = decodeText(bytes.subarray(0, 80)).trimStart().toLowerCase().startsWith('solid');
  let binary = nBin >= 0 && 84 + 50 * nBin === bytes.length;
  if (!binary && !startsSolid) {
    if (nBin < 0) throw new Error(`${input.name || 'STL'}: not an STL file (expected "solid …" text or a binary file of 84 + 50·n bytes).`);
    const fit = Math.floor((bytes.length - 84) / 50);
    warnings.push(`Binary STL header declares ${nBin} triangles but the file holds ${fit}; read ${Math.min(fit, nBin)}.`);
    binary = true;
  }
  // size the vertex table from what the file can actually hold, not from the count a header claims
  const weld = new Welder(0, binary ? Math.max(1024, Math.min(nBin, Math.floor((bytes.length - 84) / 50)) >> 1) : 4096);
  const buckets = new Map();
  let degenerate = 0, facets = 0, unitHint = null, skipped = 0, origin = null;
  const addTri = (g, a, b, c) => {
    facets++;
    if (a === b || b === c || a === c) { degenerate++; return; }
    g.push(a); g.push(b); g.push(c);
  };
  if (binary) {
    const head = decodeText(bytes.subarray(0, 80));
    const m = /UNITS\s*=\s*([A-Za-z]+)/i.exec(head);
    if (m) unitHint = m[1];
    // written by this application's exporter for models stored relative to an origin (UTM-sized coordinates)
    const og = /ORIGIN\s*=\s*(-?[\d.]+(?:e[-+]?\d+)?)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s+(-?[\d.]+(?:e[-+]?\d+)?)/i.exec(head);
    if (og) origin = [Number(og[1]), Number(og[2]), Number(og[3])];
    const n = Math.min(nBin, Math.floor((bytes.length - 84) / 50));
    const g = new Grow(Uint32Array, Math.max(3, n * 3));
    buckets.set(name, g);
    let p = 84;
    for (let i = 0; i < n; i++, p += 50) {
      const a = weld.add(dv.getFloat32(p + 12, true), dv.getFloat32(p + 16, true), dv.getFloat32(p + 20, true));
      const b = weld.add(dv.getFloat32(p + 24, true), dv.getFloat32(p + 28, true), dv.getFloat32(p + 32, true));
      const c = weld.add(dv.getFloat32(p + 36, true), dv.getFloat32(p + 40, true), dv.getFloat32(p + 44, true));
      addTri(g, a, b, c);
    }
  } else {
    const cur = lineCursor({ bytes });
    let g = null, solidNo = 0;
    const tri = [0, 0, 0];
    let k = 0;
    for (let line = cur.next(); line !== null; line = cur.next()) {
      const t = line.trim();
      if (!t) continue;
      const c0 = t.charCodeAt(0) | 32;
      if (c0 === 118 && /^vertex\s/i.test(t)) { // v
        const p = t.split(/\s+/);
        const x = Number(p[1]), y = Number(p[2]), z = Number(p[3]);
        if (!(x === x && y === y && z === z)) throw new Error(`${input.name || 'STL'}: line ${cur.lineNo}: expected "vertex x y z" with three numbers, found "${t.slice(0, 60)}".`);
        if (!g) { g = new Grow(Uint32Array, 3072); buckets.set(name, g); }
        tri[k++] = weld.add(x, y, z);
        if (k === 3) { addTri(g, tri[0], tri[1], tri[2]); k = 0; }
      } else if (c0 === 115 && /^solid\b/i.test(t)) { // s
        solidNo++;
        let nm = cleanName(t.slice(5), `solid ${solidNo}`);
        while (buckets.has(nm)) nm += ` (${solidNo})`;
        g = new Grow(Uint32Array, 3072);
        buckets.set(nm, g);
        k = 0;
      } else if (c0 === 101 && /^endloop\b/i.test(t)) {
        if (k !== 0) { skipped++; k = 0; }
      } else if (!/^(facet|outer|endfacet|endsolid|endloop)\b/i.test(t)) skipped++;
    }
    if (!weld.count) throw new Error(`${input.name || 'STL'}: no "vertex x y z" records found in an ASCII STL (expected facet … outer loop … vertex … endloop … endfacet).`);
  }
  const units = resolveUnits(opts, unitHint, 'm');
  const coords = orient(weld.finish(), opts, units.scale);
  const { indices, groups } = groupsFromBuckets(buckets);
  if (degenerate) warnings.push(`${degenerate} degenerate facet(s) (two or more coincident corners) were dropped.`);
  if (skipped) warnings.push(`${skipped} unrecognised or incomplete record(s) were skipped.`);
  const model = surfaceModel(name, coords, indices, groups.length > 1 ? groups : null);
  if (origin && !model.origin && origin.every(Number.isFinite)) model.origin = origin.map((v) => v * units.scale);
  unitAdvice(units, bboxOf(coords), suggestUnits, warnings);
  return finish([model], warnings, {
    format: 'stl', encoding: binary ? 'binary' : 'ascii', units, upAxis: upAxisMeta(opts), crs: null,
    elementTypes: { tri3: indices.length / 3 }, skipped: { degenerateFacets: degenerate, unknownRecords: skipped }, facets,
  });
}

// ---------------------------------------------------------------------------------------------------------
// OBJ
// ---------------------------------------------------------------------------------------------------------

/**
 * Wavefront OBJ: `v`, `f` (triangles, quads and n-gons; negative indices; v/vt/vn forms), `l` polylines and
 * `p` points. Groups follow `g` (within the current object), else `o`, else `usemtl`. Texture coordinates, normals,
 * free-form curves and surfaces (`curv`, `surf`, `vp` …) are not kept and are counted as skipped.
 */
export function parseObj(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const warnings = [];
  const v = new Grow(Float64Array, 3 * 4096);
  const buckets = new Map();
  const lines = [], pts = new Grow(Uint32Array, 16);
  const hist = {}, skipped = {};
  let object = '', group = '', mtl = '', g = null, gkey = null, bad = 0, unitHint = null;
  const bucket = () => {
    const key = group || object || mtl || name;
    if (key !== gkey) { gkey = key; g = buckets.get(key); if (!g) { g = new Grow(Uint32Array, 3072); buckets.set(key, g); } }
    return g;
  };
  const idx = (tok) => {
    const s = tok.indexOf('/');
    const i = parseInt(s < 0 ? tok : tok.slice(0, s), 10);
    const n = v.n / 3;
    const r = i < 0 ? n + i : i - 1;
    return r >= 0 && r < n ? r : -1;
  };
  let pending = '';
  for (let raw = cur.next(); raw !== null; raw = cur.next()) {
    let line = pending + raw;
    if (line.endsWith('\\')) { pending = line.slice(0, -1) + ' '; continue; }
    pending = '';
    const t = line.trim();
    if (!t) continue;
    if (t[0] === '#') { const m = /units?\s*[:=]?\s*(mm|millimet\w+|cm|centimet\w+|m|met\w+|in|inch\w*|ft|feet|foot)\b/i.exec(t); if (m && !unitHint) unitHint = m[1]; continue; }
    const p = t.split(/\s+/);
    switch (p[0]) {
      case 'v': v.push3(Number(p[1]), Number(p[2]), Number(p[3] ?? 0)); break;
      case 'f': {
        const ids = [];
        let ok = true;
        for (let i = 1; i < p.length; i++) { const k = idx(p[i]); if (k < 0) ok = false; else if (!ids.length || (k !== ids[ids.length - 1] && !(i === p.length - 1 && k === ids[0]))) ids.push(k); }
        if (!ok || ids.length < 3) { bad++; break; }
        hist[ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon'] = (hist[ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon'] || 0) + 1;
        pushPolygon(bucket(), ids, v.a);
        break;
      }
      case 'l': {
        const ids = [];
        for (let i = 1; i < p.length; i++) { const k = idx(p[i]); if (k >= 0) ids.push(k); }
        if (ids.length >= 2) { lines.push({ ids, name: object || group || name }); hist.polyline = (hist.polyline || 0) + 1; } else bad++;
        break;
      }
      case 'p': for (let i = 1; i < p.length; i++) { const k = idx(p[i]); if (k >= 0) { pts.push(k); hist.point1 = (hist.point1 || 0) + 1; } } break;
      case 'o': object = cleanName(t.slice(1), ''); group = ''; break;
      case 'g': group = cleanName(t.slice(1), ''); break;
      case 'usemtl': mtl = cleanName(t.slice(6), ''); break;
      case 'vn': case 'vt': case 's': case 'mtllib': break;
      default: skipped[p[0]] = (skipped[p[0]] || 0) + 1;
    }
  }
  if (!v.n) throw new Error(`${input.name || 'OBJ'}: no vertex records found (expected lines of the form "v x y z").`);
  const units = resolveUnits(opts, unitHint, 'm');
  const all = orient(v.finish(), opts, units.scale);
  const models = [];
  const { indices, groups } = groupsFromBuckets(buckets);
  if (indices.length) {
    // keep only the vertices the faces use so that curve/point vertices do not appear as stray surface nodes
    const map = new Int32Array(all.length / 3).fill(-1);
    const used = new Grow(Float64Array, all.length);
    for (let i = 0; i < indices.length; i++) {
      const k = indices[i];
      if (map[k] < 0) { map[k] = used.n / 3; used.push3(all[3 * k], all[3 * k + 1], all[3 * k + 2]); }
      indices[i] = map[k];
    }
    models.push(surfaceModel(name, used.finish(), indices, groups.length > 1 ? groups : null));
  }
  lines.forEach((l, i) => {
    const pp = new Float64Array(l.ids.length * 3);
    l.ids.forEach((k, j) => { pp[3 * j] = all[3 * k]; pp[3 * j + 1] = all[3 * k + 1]; pp[3 * j + 2] = all[3 * k + 2]; });
    models.push({ kind: 'polyline', name: lines.length > 1 ? `${l.name} ${i + 1}` : l.name, points: pp });
  });
  if (pts.n || !models.length) {
    const ids = pts.n ? pts.view() : null;
    const n = ids ? ids.length : all.length / 3;
    const pp = new Float64Array(n * 3);
    for (let j = 0; j < n; j++) { const k = ids ? ids[j] : j; pp[3 * j] = all[3 * k]; pp[3 * j + 1] = all[3 * k + 1]; pp[3 * j + 2] = all[3 * k + 2]; }
    models.push(pointsModel(models.length ? `${name} (points)` : name, pp));
    if (!ids) warnings.push('The file has vertices but no faces, lines or points; the vertices are returned as a point set.');
  }
  if (bad) { skipped.invalidElements = bad; warnings.push(`${bad} face/line record(s) with missing or out-of-range vertex indices were skipped.`); }
  for (const k in skipped) if (k !== 'invalidElements') warnings.push(`Skipped ${skipped[k]} "${k}" record(s) (not geometry this reader keeps).`);
  unitAdvice(units, bboxOf(all), suggestUnits, warnings);
  return finish(models, warnings, { format: 'obj', units, upAxis: upAxisMeta(opts), crs: null, elementTypes: hist, skipped });
}

// ---------------------------------------------------------------------------------------------------------
// PLY
// ---------------------------------------------------------------------------------------------------------

const PLY_TYPES = {
  char: [1, 'Int8'], int8: [1, 'Int8'], uchar: [1, 'Uint8'], uint8: [1, 'Uint8'],
  short: [2, 'Int16'], int16: [2, 'Int16'], ushort: [2, 'Uint16'], uint16: [2, 'Uint16'],
  int: [4, 'Int32'], int32: [4, 'Int32'], uint: [4, 'Uint32'], uint32: [4, 'Uint32'],
  float: [4, 'Float32'], float32: [4, 'Float32'], double: [8, 'Float64'], float64: [8, 'Float64'],
};

/**
 * Stanford PLY, ASCII or binary (little/big endian). Reads the `vertex` element (x, y, z plus every other
 * scalar property as a per-vertex attribute: colour, normals, intensity …) and the `face` element
 * (`vertex_indices` lists; polygons are triangulated). A file without faces is a point cloud.
 * Other elements (edge, material, tristrips …) are read past and reported as skipped.
 */
export function parsePly(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  const warnings = [];
  const cur = lineCursor({ bytes });
  const first = cur.next();
  if (first === null || first.trim() !== 'ply') throw new Error(`${input.name || 'PLY'}: not a PLY file (expected the first line to be "ply").`);
  let fmt = null, unitHint = null;
  const elements = [];
  for (;;) {
    const l = cur.next();
    if (l === null) throw new Error(`${input.name || 'PLY'}: header is not terminated (expected a line "end_header").`);
    const p = l.trim().split(/\s+/);
    if (p[0] === 'end_header') break;
    if (p[0] === 'format') fmt = p[1];
    else if (p[0] === 'comment' || p[0] === 'obj_info') { const m = /units?\s*[:=]?\s*([A-Za-z]+)/i.exec(l); if (m && !unitHint) unitHint = m[1]; } else if (p[0] === 'element') elements.push({ name: p[1], count: Number(p[2]), props: [] });
    else if (p[0] === 'property' && elements.length) {
      const e = elements[elements.length - 1];
      if (p[1] === 'list') {
        if (!PLY_TYPES[p[2]] || !PLY_TYPES[p[3]]) throw new Error(`${input.name || 'PLY'}: unknown list property types "${p[2]} ${p[3]}".`);
        e.props.push({ name: p[4], list: true, ct: PLY_TYPES[p[2]], it: PLY_TYPES[p[3]] });
      } else {
        if (!PLY_TYPES[p[1]]) throw new Error(`${input.name || 'PLY'}: unknown property type "${p[1]}".`);
        e.props.push({ name: p[2], list: false, t: PLY_TYPES[p[1]], type: p[1] });
      }
    }
  }
  if (!['ascii', 'binary_little_endian', 'binary_big_endian'].includes(fmt)) throw new Error(`${input.name || 'PLY'}: unsupported format "${fmt}" (expected ascii, binary_little_endian or binary_big_endian).`);
  const ascii = fmt === 'ascii', le = fmt !== 'binary_big_endian';
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let off = cur.offset;
  let coords = null;
  const attrs = {};
  const tris = new Grow(Uint32Array, 3072);
  const hist = {}, skipped = {};
  let bad = 0;
  const truncated = (e) => new Error(`${input.name || 'PLY'}: file ends inside element "${e.name}" (header declares ${e.count} entries).`);

  for (const e of elements) {
    const isVertex = e.name === 'vertex', isFace = e.name === 'face';
    if (!(e.count >= 0) || e.count > bytes.length) throw new Error(`${input.name || 'PLY'}: element "${e.name}" declares ${e.count} entries, more than the file can hold (corrupt header).`);
    if (isVertex) {
      const ix = e.props.findIndex((p) => p.name === 'x'), iy = e.props.findIndex((p) => p.name === 'y'), iz = e.props.findIndex((p) => p.name === 'z');
      if (ix < 0 || iy < 0) throw new Error(`${input.name || 'PLY'}: the vertex element has no x/y properties.`);
      coords = new Float64Array(e.count * 3);
      const cols = e.props.map((p, k) => {
        if (p.list || k === ix || k === iy || k === iz) return null;
        const T = p.t[1] === 'Float64' ? Float64Array : p.t[1] === 'Float32' ? Float32Array : globalThis[`${p.t[1]}Array`];
        return (attrs[p.name] = new T(e.count));
      });
      for (let i = 0; i < e.count; i++) {
        let tok = null;
        if (ascii) { const l = cur.nextNonBlank(); if (l === null) throw truncated(e); tok = l.split(/\s+/); }
        let t = 0;
        for (let k = 0; k < e.props.length; k++) {
          const p = e.props[k];
          let val;
          if (p.list) {
            if (ascii) { const n = Number(tok[t++]); t += n; } else { const n = dv[`get${p.ct[1]}`](off, le); off += p.ct[0] + n * p.it[0]; }
            continue;
          }
          if (ascii) val = Number(tok[t++]);
          else { if (off + p.t[0] > bytes.length) throw truncated(e); val = dv[`get${p.t[1]}`](off, le); off += p.t[0]; }
          if (k === ix) coords[3 * i] = val; else if (k === iy) coords[3 * i + 1] = val; else if (k === iz) coords[3 * i + 2] = val; else cols[k][i] = val;
        }
      }
    } else if (isFace) {
      const li = e.props.findIndex((p) => p.list && /^vertex_ind/.test(p.name));
      if (li < 0) { skipped.face = e.count; warnings.push('The face element has no vertex_indices list; faces were skipped.'); }
      const nV = coords ? coords.length / 3 : 0;
      for (let i = 0; i < e.count; i++) {
        let tok = null;
        if (ascii) { const l = cur.nextNonBlank(); if (l === null) throw truncated(e); tok = l.split(/\s+/); }
        let t = 0;
        for (let k = 0; k < e.props.length; k++) {
          const p = e.props[k];
          if (!p.list) { if (ascii) t++; else off += p.t[0]; continue; }
          let n;
          if (ascii) n = Number(tok[t++]); else { if (off + p.ct[0] > bytes.length) throw truncated(e); n = dv[`get${p.ct[1]}`](off, le); off += p.ct[0]; }
          if (k !== li) { if (ascii) t += n; else off += n * p.it[0]; continue; }
          const ids = new Array(n);
          let ok = n >= 3;
          if (!ascii && off + n * p.it[0] > bytes.length) throw truncated(e);
          for (let j = 0; j < n; j++) {
            ids[j] = ascii ? Number(tok[t++]) : dv[`get${p.it[1]}`](off + j * p.it[0], le);
            if (!(ids[j] >= 0 && ids[j] < nV)) ok = false;
          }
          if (!ascii) off += n * p.it[0];
          if (!ok) { bad++; continue; }
          const key = n === 3 ? 'tri3' : n === 4 ? 'quad4' : 'polygon';
          hist[key] = (hist[key] || 0) + 1;
          pushPolygon(tris, ids, coords);
        }
      }
    } else {
      skipped[e.name] = e.count;
      for (let i = 0; i < e.count; i++) {
        if (ascii) { if (cur.nextNonBlank() === null) break; continue; }
        for (const p of e.props) {
          if (p.list) { const n = dv[`get${p.ct[1]}`](off, le); off += p.ct[0] + n * p.it[0]; } else off += p.t[0];
        }
      }
      if (e.count) warnings.push(`Skipped ${e.count} "${e.name}" element(s) (only vertex and face elements are read).`);
    }
  }
  if (!coords || !coords.length) throw new Error(`${input.name || 'PLY'}: no vertex element with data found (expected "element vertex N" with x, y, z properties).`);
  if (bad) { skipped.invalidFaces = bad; warnings.push(`${bad} face(s) with fewer than three or out-of-range vertex indices were skipped.`); }
  const units = resolveUnits(opts, unitHint, 'm');
  orient(coords, opts, units.scale);
  const keep = {};
  if (attrs.red && attrs.green && attrs.blue) {
    const n = coords.length / 3, rgb = new Uint8Array(3 * n);
    const s = attrs.red instanceof Uint8Array ? 1 : attrs.red instanceof Uint16Array ? 1 / 257 : 255;
    for (let i = 0; i < n; i++) { rgb[3 * i] = attrs.red[i] * s; rgb[3 * i + 1] = attrs.green[i] * s; rgb[3 * i + 2] = attrs.blue[i] * s; }
    keep.rgb = rgb;
  }
  for (const k in attrs) if (!['red', 'green', 'blue'].includes(k) || !keep.rgb) keep[k] = attrs[k];
  let model;
  if (tris.n) {
    model = surfaceModel(name, coords, tris.finish(), null, Object.keys(keep).length ? { attrs: keep } : null);
  } else {
    hist.point1 = coords.length / 3;
    model = pointsModel(name, coords, keep);
  }
  unitAdvice(units, bboxOf(coords), suggestUnits, warnings);
  return finish([model], warnings, { format: 'ply', encoding: fmt, units, upAxis: upAxisMeta(opts), crs: null, elementTypes: hist, skipped, vertexProperties: Object.keys(attrs) });
}

// ---------------------------------------------------------------------------------------------------------
// OFF
// ---------------------------------------------------------------------------------------------------------

/**
 * Geomview OFF (also COFF / NOFF / CNOFF / STOFF prefixes; per-vertex colours and normals are read past).
 * Faces with more than three vertices are triangulated. 4-D and n-D variants (4OFF, nOFF) are rejected.
 */
export function parseOff(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const warnings = [];
  const queue = [];
  let qi = 0;
  const tokenLine = () => {
    for (;;) {
      const l = cur.next();
      if (l === null) return null;
      const h = l.indexOf('#');
      const t = (h >= 0 ? l.slice(0, h) : l).trim();
      if (t) return t.split(/\s+/);
    }
  };
  const tok = () => {
    while (qi >= queue.length) { const l = tokenLine(); if (l === null) return null; queue.length = 0; qi = 0; for (const s of l) queue.push(s); }
    return queue[qi++];
  };
  const head = tok();
  const m = head === null ? null : /^(ST)?(C)?(N)?(4)?(n)?OFF$/.exec(head);
  if (!m) throw new Error(`${input.name || 'OFF'}: not an OFF file (expected the keyword OFF, COFF, NOFF or CNOFF on the first line).`);
  if (m[4] || m[5]) throw new Error(`${input.name || 'OFF'}: ${head} (4-D / n-D homogeneous coordinates) is not supported; export plain 3-D OFF.`);
  const perVertex = 3 + (m[3] ? 3 : 0) + (m[1] ? 2 : 0);
  const nV = Number(tok()), nF = Number(tok());
  tok(); // edge count, unused by the format itself
  if (!(nV > 0) || !(nF >= 0) || nV > cur.len || nF > cur.len) throw new Error(`${input.name || 'OFF'}: expected "vertices faces edges" counts after the OFF keyword (found ${nV}, ${nF}).`);
  const coords = new Float64Array(nV * 3);
  // vertex rows may carry colours (C) of unknown arity, so they are read line by line
  queue.length = 0; qi = 0;
  for (let i = 0; i < nV; i++) {
    const l = tokenLine();
    if (l === null || l.length < 3) throw new Error(`${input.name || 'OFF'}: vertex ${i + 1} of ${nV}: expected at least ${perVertex} numbers per line.`);
    coords[3 * i] = Number(l[0]); coords[3 * i + 1] = Number(l[1]); coords[3 * i + 2] = Number(l[2]);
  }
  const tris = new Grow(Uint32Array, Math.max(3, nF * 3));
  const hist = {};
  let bad = 0, read = 0;
  for (let i = 0; i < nF; i++) {
    const l = tokenLine();
    if (l === null) break;
    read++;
    const n = Number(l[0]);
    const ids = [];
    let ok = n >= 3 && l.length > n;
    for (let j = 1; j <= n && ok; j++) { const k = Number(l[j]); if (!(k >= 0 && k < nV)) ok = false; ids.push(k); }
    if (!ok) { bad++; continue; }
    const key = n === 3 ? 'tri3' : n === 4 ? 'quad4' : 'polygon';
    hist[key] = (hist[key] || 0) + 1;
    pushPolygon(tris, ids, coords);
  }
  if (read < nF) warnings.push(`The header declares ${nF} faces but the file ends after ${read}.`);
  if (bad) warnings.push(`${bad} face(s) with fewer than three or out-of-range vertex indices were skipped.`);
  const units = resolveUnits(opts, null, 'm');
  orient(coords, opts, units.scale);
  const model = tris.n ? surfaceModel(name, coords, tris.finish()) : pointsModel(name, coords);
  if (!tris.n) hist.point1 = nV;
  unitAdvice(units, bboxOf(coords), suggestUnits, warnings);
  return finish([model], warnings, { format: 'off', units, upAxis: upAxisMeta(opts), crs: null, elementTypes: hist, skipped: bad ? { invalidFaces: bad } : {} });
}
