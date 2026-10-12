// Exporters for the geometry model kinds of docs/CONTRACT.md §4.
//
// Every exporter writes world coordinates in metres: models that store single-precision positions relative to
// an integer-metre `origin` (see packPositions in importers/common.js) get the origin added back — except
// binary STL, whose 32-bit floats cannot hold UTM-sized coordinates: there the positions are written as stored
// and the origin is recorded in the 80-byte header ("ORIGIN=x y z"), which the STL importer restores.
//
//   format      kinds                                         output
//   stl         surface, grid                                 Uint8Array (binary STL)
//   obj         surface, grid, frame, polyline, points        string
//   ply         surface, grid, points                         Uint8Array (binary little endian) or string (ascii)
//   vtk         surface, volume, frame, polyline, points, grid string (legacy ASCII, nodal scalars / vectors)
//   csv         points, polyline, frame, grid, section        string
//   geojson     points, polyline, frame, section              string
//   esri-ascii  grid                                          string
//   json        all kinds (lossless application format)       string

import { base64Encode, ELEM } from './importers/common.js';
import { localToLonLat } from './units.js';
import { OA_FORMAT, OA_VERSION } from './importers/json.js';

const enc = (s) => new TextEncoder().encode(s);
const fmt = (v) => (Number.isFinite(v) ? String(Number(v.toPrecision(15))) : 'NaN');

function unsupported(format, model, kinds) {
  return new Error(`${format} export does not handle a "${model && model.kind}" model (supported kinds: ${kinds.join(', ')}).`);
}

/** World coordinates (Float64Array, flat xyz) of a model's nodes. */
export function worldCoords(model) {
  const src = model.kind === 'frame' ? model.nodes : model.kind === 'polyline' ? model.points : model.positions;
  const o = model.origin || [0, 0, 0];
  const out = new Float64Array(src.length);
  for (let i = 0; i < src.length; i += 3) { out[i] = src[i] + o[0]; out[i + 1] = src[i + 1] + o[1]; out[i + 2] = src[i + 2] + o[2]; }
  return out;
}

/** Triangulate a grid model (cells with a no-data corner are left out). Returns a surface-like object. */
export function gridToSurface(g) {
  const { nx, ny } = g;
  const pos = new Float64Array(nx * ny * 3);
  for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) { const k = 3 * (iy * nx + ix); pos[k] = g.x0 + ix * g.dx; pos[k + 1] = g.y0 + iy * g.dy; pos[k + 2] = g.z[iy * nx + ix]; }
  const idx = [];
  for (let iy = 0; iy + 1 < ny; iy++) for (let ix = 0; ix + 1 < nx; ix++) {
    const a = iy * nx + ix, b = a + 1, c = a + nx + 1, d = a + nx;
    const ok = (k) => g.z[k] === g.z[k];
    if (ok(a) && ok(b) && ok(c)) idx.push(a, b, c);
    if (ok(a) && ok(c) && ok(d)) idx.push(a, c, d);
  }
  return { kind: 'surface', name: g.name, positions: pos, indices: Uint32Array.from(idx) };
}

const asSurface = (m, format) => { if (m.kind === 'surface') return m; if (m.kind === 'grid') return gridToSurface(m); throw unsupported(format, m, ['surface', 'grid']); };

/** Binary STL. */
export function exportStl(model) {
  const s = asSurface(model, 'STL');
  const p = s.positions, t = s.indices, n = t.length / 3;
  const out = new Uint8Array(84 + 50 * n);
  const o = s.origin;
  out.set(enc(`UNITS=m${o ? ` ORIGIN=${o[0]} ${o[1]} ${o[2]}` : ''}`).subarray(0, 80));
  const dv = new DataView(out.buffer);
  dv.setUint32(80, n, true);
  for (let f = 0, q = 84; f < n; f++, q += 50) {
    const a = 3 * t[3 * f], b = 3 * t[3 * f + 1], c = 3 * t[3 * f + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2], vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const l = Math.hypot(nx, ny, nz) || 1;
    nx /= l; ny /= l; nz /= l;
    dv.setFloat32(q, nx, true); dv.setFloat32(q + 4, ny, true); dv.setFloat32(q + 8, nz, true);
    for (let k = 0; k < 3; k++) { const v = 3 * t[3 * f + k]; dv.setFloat32(q + 12 + 12 * k, p[v], true); dv.setFloat32(q + 16 + 12 * k, p[v + 1], true); dv.setFloat32(q + 20 + 12 * k, p[v + 2], true); }
  }
  return out;
}

/** Wavefront OBJ (one or several models in one file). */
export function exportObj(models) {
  const list = Array.isArray(models) ? models : [models];
  const lines = ['# units: m'];
  let base = 0;
  for (const src of list) {
    const m = src.kind === 'grid' ? gridToSurface(src) : src;
    if (!['surface', 'frame', 'polyline', 'points'].includes(m.kind)) throw unsupported('OBJ', src, ['surface', 'grid', 'frame', 'polyline', 'points']);
    const c = m.kind === 'surface' && !m.origin ? m.positions : worldCoords(m);
    lines.push(`o ${String(m.name || m.kind).replace(/\s+/g, '_')}`);
    for (let i = 0; i < c.length; i += 3) lines.push(`v ${fmt(c[i])} ${fmt(c[i + 1])} ${fmt(c[i + 2])}`);
    if (m.kind === 'surface') {
      const groups = m.groups && m.groups.length ? m.groups : [{ name: null, start: 0, count: m.indices.length }];
      for (const g of groups) { if (g.name !== null) lines.push(`g ${String(g.name).replace(/\s+/g, '_')}`); for (let i = g.start; i < g.start + g.count; i += 3) lines.push(`f ${m.indices[i] + 1 + base} ${m.indices[i + 1] + 1 + base} ${m.indices[i + 2] + 1 + base}`); }
    } else if (m.kind === 'frame') for (const mb of m.members) lines.push(`l ${mb.a + 1 + base} ${mb.b + 1 + base}`);
    else if (m.kind === 'polyline') { const n = c.length / 3; let s = 'l'; for (let i = 0; i < n; i++) s += ` ${i + 1 + base}`; if (m.closed) s += ` ${1 + base}`; lines.push(s); } else { const n = c.length / 3; for (let i = 0; i < n; i += 64) { let s = 'p'; for (let k = i; k < Math.min(n, i + 64); k++) s += ` ${k + 1 + base}`; lines.push(s); } }
    base += c.length / 3;
  }
  return lines.join('\n') + '\n';
}

/** Stanford PLY: binary little endian (default) or ASCII. Vertex colours (`attrs.rgb`) are written when present. */
export function exportPly(model, opts = {}) {
  const m = model.kind === 'grid' ? gridToSurface(model) : model;
  if (m.kind !== 'surface' && m.kind !== 'points') throw unsupported('PLY', model, ['surface', 'grid', 'points']);
  const c = worldCoords(m);
  const n = c.length / 3, nf = m.kind === 'surface' ? m.indices.length / 3 : 0;
  const rgb = m.attrs && m.attrs.rgb && m.attrs.rgb.length === 3 * n ? m.attrs.rgb : null;
  const ascii = opts.binary === false;
  const head = `ply\nformat ${ascii ? 'ascii' : 'binary_little_endian'} 1.0\ncomment units: m\nelement vertex ${n}\nproperty double x\nproperty double y\nproperty double z\n${rgb ? 'property uchar red\nproperty uchar green\nproperty uchar blue\n' : ''}${nf ? `element face ${nf}\nproperty list uchar int vertex_indices\n` : ''}end_header\n`;
  if (ascii) {
    const lines = [];
    for (let i = 0; i < n; i++) lines.push(`${fmt(c[3 * i])} ${fmt(c[3 * i + 1])} ${fmt(c[3 * i + 2])}${rgb ? ` ${rgb[3 * i]} ${rgb[3 * i + 1]} ${rgb[3 * i + 2]}` : ''}`);
    for (let f = 0; f < nf; f++) lines.push(`3 ${m.indices[3 * f]} ${m.indices[3 * f + 1]} ${m.indices[3 * f + 2]}`);
    return head + lines.join('\n') + '\n';
  }
  const h = enc(head);
  const vs = 24 + (rgb ? 3 : 0);
  const out = new Uint8Array(h.length + n * vs + nf * 13);
  out.set(h);
  const dv = new DataView(out.buffer);
  let p = h.length;
  for (let i = 0; i < n; i++) { dv.setFloat64(p, c[3 * i], true); dv.setFloat64(p + 8, c[3 * i + 1], true); dv.setFloat64(p + 16, c[3 * i + 2], true); p += 24; if (rgb) { out[p++] = rgb[3 * i]; out[p++] = rgb[3 * i + 1]; out[p++] = rgb[3 * i + 2]; } }
  for (let f = 0; f < nf; f++) { out[p++] = 3; for (let k = 0; k < 3; k++) { dv.setInt32(p, m.indices[3 * f + k], true); p += 4; } }
  return out;
}

const VTK_TYPE = { tet4: [10], tet10: [24], hex8: [12], hex20: [25], wedge6: [13, [0, 2, 1, 3, 5, 4]], pyr5: [14], tri3: [5], tri6: [22], quad4: [9], quad8: [23], quad9: [28] };

/**
 * Legacy VTK, ASCII. Nodal fields: `opts.scalars = { name: ArrayLike }` with one value (SCALARS) or three
 * values (VECTORS) per node; a model's own `nodeData` (as attached by the VTK importer) is written too.
 */
export function exportVtk(model, opts = {}) {
  const L = ['# vtk DataFile Version 3.0', String(opts.title || model.name || 'geometry').replace(/[\r\n]+/g, ' ').slice(0, 255), 'ASCII'];
  let n;
  if (model.kind === 'grid') {
    n = model.nx * model.ny;
    L.push('DATASET STRUCTURED_POINTS', `DIMENSIONS ${model.nx} ${model.ny} 1`, `ORIGIN ${fmt(model.x0)} ${fmt(model.y0)} 0`, `SPACING ${fmt(model.dx)} ${fmt(model.dy)} 1`);
  } else {
    if (!['surface', 'volume', 'frame', 'polyline', 'points'].includes(model.kind)) throw unsupported('VTK', model, ['surface', 'volume', 'frame', 'polyline', 'points', 'grid']);
    const c = worldCoords(model);
    n = c.length / 3;
    const pts = [];
    for (let i = 0; i < c.length; i += 3) pts.push(`${fmt(c[i])} ${fmt(c[i + 1])} ${fmt(c[i + 2])}`);
    const header = (ds) => L.push(`DATASET ${ds}`, `POINTS ${n} double`, ...pts);
    if (model.kind === 'volume' || (model.kind === 'surface' && model.elements && opts.elements)) {
      header('UNSTRUCTURED_GRID');
      const cells = [], types = [];
      let size = 0;
      for (const b of model.kind === 'volume' ? model.cells : model.elements) {
        const def = VTK_TYPE[b.type], k = ELEM[b.type].n;
        if (!def) throw new Error(`VTK export: cell type "${b.type}" has no VTK equivalent.`);
        for (let e = 0; e < b.conn.length / k; e++) { const ids = Array.from(b.conn.subarray(e * k, (e + 1) * k)); cells.push(`${k} ${(def[1] ? def[1].map((q) => ids[q]) : ids).join(' ')}`); types.push(def[0]); size += k + 1; }
      }
      L.push(`CELLS ${cells.length} ${size}`, ...cells, `CELL_TYPES ${types.length}`, ...types);
    } else {
      header('POLYDATA');
      if (model.kind === 'surface') { const t = model.indices; L.push(`POLYGONS ${t.length / 3} ${(t.length / 3) * 4}`); for (let i = 0; i < t.length; i += 3) L.push(`3 ${t[i]} ${t[i + 1]} ${t[i + 2]}`); } else if (model.kind === 'frame') { L.push(`LINES ${model.members.length} ${model.members.length * 3}`); for (const m of model.members) L.push(`2 ${m.a} ${m.b}`); } else if (model.kind === 'polyline') { const k = n + (model.closed ? 1 : 0); L.push(`LINES 1 ${k + 1}`, `${k} ${Array.from({ length: n }, (_, i) => i).join(' ')}${model.closed ? ' 0' : ''}`); } else { L.push(`VERTICES ${n} ${2 * n}`); for (let i = 0; i < n; i++) L.push(`1 ${i}`); }
    }
  }
  const fields = {};
  if (model.kind === 'grid') fields.elevation = model.z;
  if (model.nodeData) for (const k in model.nodeData) fields[k] = model.nodeData[k];
  for (const k in opts.scalars || {}) fields[k] = opts.scalars[k];
  const names = Object.keys(fields);
  if (names.length) {
    L.push(`POINT_DATA ${n}`);
    for (const nm of names) {
      const a = fields[nm];
      const safe = nm.replace(/\s+/g, '_');
      if (a.length === n) { L.push(`SCALARS ${safe} double 1`, 'LOOKUP_TABLE default'); for (let i = 0; i < n; i += 9) L.push(Array.from(a.subarray ? a.subarray(i, i + 9) : a.slice(i, i + 9), fmt).join(' ')); } else if (a.length === 3 * n) { L.push(`VECTORS ${safe} double`); for (let i = 0; i < n; i++) L.push(`${fmt(a[3 * i])} ${fmt(a[3 * i + 1])} ${fmt(a[3 * i + 2])}`); } else throw new RangeError(`VTK export: nodal field "${nm}" has ${a.length} values for ${n} nodes (expected ${n} or ${3 * n}).`);
    }
  }
  return L.join('\n') + '\n';
}

const csvCell = (v) => { const s = v === undefined || v === null ? '' : typeof v === 'number' ? (Number.isFinite(v) ? fmt(v) : '') : String(v); return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };

/** CSV. Frames are written as a [nodes] table and a [members] table, which the table importer reads back. */
export function exportCsv(model) {
  const rows = [];
  const push = (...c) => rows.push(c.map(csvCell).join(','));
  switch (model.kind) {
    case 'points': case 'polyline': {
      const c = worldCoords(model), n = c.length / 3;
      const attrs = model.attrs || {};
      const keys = Object.keys(attrs).filter((k) => attrs[k] && attrs[k].length === n && !(k === 'kp' && model.kind === 'polyline'));
      const kp = model.kind === 'polyline' && attrs.kp && attrs.kp.length === n ? attrs.kp : null;
      push(...(kp ? ['chainage (m)'] : []), 'x (m)', 'y (m)', 'z (m)', ...keys.map((k) => (k === 'names' ? 'name' : k)));
      for (let i = 0; i < n; i++) push(...(kp ? [kp[i]] : []), c[3 * i], c[3 * i + 1], c[3 * i + 2], ...keys.map((k) => attrs[k][i]));
      break;
    }
    case 'frame': {
      const c = worldCoords(model), n = c.length / 3;
      rows.push('[nodes]');
      push('id', 'x (m)', 'y (m)', 'z (m)');
      for (let i = 0; i < n; i++) push(model.nodeNames && model.nodeNames[i] ? model.nodeNames[i] : i + 1, c[3 * i], c[3 * i + 1], c[3 * i + 2]);
      rows.push('[members]');
      push('id', 'node_a', 'node_b', 'D (m)', 't (m)', 'area (m2)', 'Iy (m4)', 'Iz (m4)', 'J (m4)', 'section', 'material', 'group');
      const nid = (i) => (model.nodeNames && model.nodeNames[i] ? model.nodeNames[i] : i + 1);
      model.members.forEach((m, i) => { const s = m.section || {}; push(m.id ?? i + 1, nid(m.a), nid(m.b), s.type === 'tube' ? s.D : '', s.type === 'tube' ? s.t : '', s.A ?? '', s.Iy ?? '', s.Iz ?? '', s.J ?? '', s.name ?? s.designation ?? m.sectionName ?? '', m.material ?? '', m.group ?? ''); });
      break;
    }
    case 'grid': {
      push('x (m)', 'y (m)', 'z (m)');
      for (let iy = 0; iy < model.ny; iy++) for (let ix = 0; ix < model.nx; ix++) { const z = model.z[iy * model.nx + ix]; if (z === z) push(model.x0 + ix * model.dx, model.y0 + iy * model.dy, z); }
      break;
    }
    case 'section': {
      push('loop', 'x (m)', 'y (m)');
      model.loops.forEach((l, k) => { for (let i = 0; i < l.length; i += 2) push(k + 1, l[i], l[i + 1]); });
      break;
    }
    default: throw unsupported('CSV', model, ['points', 'polyline', 'frame', 'grid', 'section']);
  }
  return rows.join('\n') + '\n';
}

/**
 * GeoJSON FeatureCollection.
 * With `opts.origin = { lon, lat }` the metric coordinates are taken as a local tangent plane about that
 * origin (the importers' default projection) and written as WGS 84 longitude/latitude, as RFC 7946 requires.
 * Without it the projected coordinates are written unchanged with a legacy named `crs` member when `opts.crs`
 * or the model's own `crs` gives an EPSG code — such a file is outside RFC 7946 but is read by GIS software
 * and by this application.
 */
export function exportGeoJson(models, opts = {}) {
  const list = Array.isArray(models) ? models : [models];
  const o = opts.origin && Number.isFinite(opts.origin.lon) && Number.isFinite(opts.origin.lat) ? opts.origin : null;
  // lines always carry their elevation; points and polygon rings only when it is not zero
  const P = (x, y, z, three) => { const q = o ? localToLonLat(x, y, o.lon, o.lat) : [x, y]; return three || z ? [q[0], q[1], z] : [q[0], q[1]]; };
  const features = [];
  let crs = opts.crs || null;
  for (const m of list) {
    if (m.crs && !crs) crs = m.crs;
    if (m.kind === 'points') {
      const c = worldCoords(m), n = c.length / 3, a = m.attrs || {};
      for (let i = 0; i < n; i++) { const props = {}; for (const k in a) if (a[k] && a[k].length === n) props[k === 'names' ? 'name' : k] = typeof a[k][i] === 'number' && !Number.isFinite(a[k][i]) ? null : a[k][i]; features.push({ type: 'Feature', properties: props, geometry: { type: 'Point', coordinates: P(c[3 * i], c[3 * i + 1], c[3 * i + 2]) } }); }
    } else if (m.kind === 'polyline') {
      const c = worldCoords(m), coords = [];
      for (let i = 0; i < c.length; i += 3) coords.push(P(c[i], c[i + 1], c[i + 2], true));
      if (m.closed) coords.push(coords[0]);
      features.push({ type: 'Feature', properties: { name: m.name, ...(m.properties || {}) }, geometry: { type: 'LineString', coordinates: coords } });
    } else if (m.kind === 'frame') {
      const c = worldCoords(m);
      for (const mb of m.members) { const props = { group: mb.group ?? null, material: mb.material ?? null }; if (mb.id !== undefined) props.id = mb.id; if (mb.section && mb.section.type === 'tube') { props.D = mb.section.D; props.t = mb.section.t; } features.push({ type: 'Feature', properties: props, geometry: { type: 'LineString', coordinates: [P(c[3 * mb.a], c[3 * mb.a + 1], c[3 * mb.a + 2], true), P(c[3 * mb.b], c[3 * mb.b + 1], c[3 * mb.b + 2], true)] } }); }
    } else if (m.kind === 'section') {
      // the first loop is the outer boundary, the others are holes
      const rings = m.loops.map((l) => { const r = []; for (let i = 0; i < l.length; i += 2) r.push(P(l[i], l[i + 1], 0)); r.push(r[0]); return r; });
      features.push({ type: 'Feature', properties: { name: m.name, ...(m.properties || {}) }, geometry: { type: 'Polygon', coordinates: rings } });
    } else throw unsupported('GeoJSON', m, ['points', 'polyline', 'frame', 'section']);
  }
  const fc = { type: 'FeatureCollection', features };
  if (!o && crs) { const e = /(\d{4,6})/.exec(String(crs)); if (e) fc.crs = { type: 'name', properties: { name: `urn:ogc:def:crs:EPSG::${e[1]}` } }; }
  return JSON.stringify(fc);
}

/** ESRI ASCII grid. Rows are written north to south; no-data cells as `opts.nodata` (default −9999). */
export function exportEsriAscii(grid, opts = {}) {
  if (!grid || grid.kind !== 'grid') throw unsupported('ESRI ASCII grid', grid, ['grid']);
  const nodata = opts.nodata ?? -9999;
  const L = [`ncols ${grid.nx}`, `nrows ${grid.ny}`, `xllcorner ${fmt(grid.x0 - grid.dx / 2)}`, `yllcorner ${fmt(grid.y0 - grid.dy / 2)}`];
  // square cells use the classic keyword; otherwise the dx / dy form understood by GDAL
  if (Math.abs(grid.dx - grid.dy) <= 1e-12 * Math.abs(grid.dx)) L.push(`cellsize ${fmt(grid.dx)}`); else L.push(`dx ${fmt(grid.dx)}`, `dy ${fmt(grid.dy)}`);
  L.push(`NODATA_value ${nodata}`);
  for (let iy = grid.ny - 1; iy >= 0; iy--) {
    const row = new Array(grid.nx);
    for (let ix = 0; ix < grid.nx; ix++) { const v = grid.z[iy * grid.nx + ix]; row[ix] = v === v ? String(Number(v.toPrecision(7))) : String(nodata); }
    L.push(row.join(' '));
  }
  return L.join('\n') + '\n';
}

function pack(v) {
  if (ArrayBuffer.isView(v) && !(v instanceof DataView)) {
    // NaN cannot be written in JSON, so typed arrays always travel as base64 of their little-endian bytes
    return { $typed: v.constructor.name, base64: base64Encode(new Uint8Array(v.buffer, v.byteOffset, v.byteLength)) };
  }
  if (Array.isArray(v)) return v.map(pack);
  if (v && typeof v === 'object') { const o = {}; for (const k in v) if (v[k] !== undefined && typeof v[k] !== 'function') o[k] = pack(v[k]); return o; }
  return typeof v === 'number' && !Number.isFinite(v) ? null : v;
}

/** The application's own lossless JSON format (read back by importers/json.js). */
export function exportJson(models, opts = {}) {
  const list = Array.isArray(models) ? models : [models];
  return JSON.stringify({ format: OA_FORMAT, version: OA_VERSION, units: 'm', crs: opts.crs ?? null, source: opts.source ?? null, models: list.map(pack) });
}

/** Export formats offered to the user. */
export const EXPORTERS = Object.freeze([
  { id: 'stl', name: 'STL (binary)', extension: 'stl', mime: 'model/stl', kinds: ['surface', 'grid'], multi: false, write: (m) => exportStl(m) },
  { id: 'obj', name: 'Wavefront OBJ', extension: 'obj', mime: 'model/obj', kinds: ['surface', 'grid', 'frame', 'polyline', 'points'], multi: true, write: (m) => exportObj(m) },
  { id: 'ply', name: 'Stanford PLY', extension: 'ply', mime: 'application/octet-stream', kinds: ['surface', 'grid', 'points'], multi: false, write: (m, o) => exportPly(m, o) },
  { id: 'vtk', name: 'VTK legacy (ASCII)', extension: 'vtk', mime: 'text/plain', kinds: ['surface', 'volume', 'frame', 'polyline', 'points', 'grid'], multi: false, write: (m, o) => exportVtk(m, o) },
  { id: 'csv', name: 'CSV table', extension: 'csv', mime: 'text/csv', kinds: ['points', 'polyline', 'frame', 'grid', 'section'], multi: false, write: (m) => exportCsv(m) },
  { id: 'geojson', name: 'GeoJSON', extension: 'geojson', mime: 'application/geo+json', kinds: ['points', 'polyline', 'frame', 'section'], multi: true, write: (m, o) => exportGeoJson(m, o) },
  { id: 'esri-ascii', name: 'ESRI ASCII grid', extension: 'asc', mime: 'text/plain', kinds: ['grid'], multi: false, write: (m, o) => exportEsriAscii(m, o) },
  { id: 'json', name: 'Application geometry JSON', extension: 'json', mime: 'application/json', kinds: ['surface', 'volume', 'frame', 'polyline', 'grid', 'points', 'section'], multi: true, write: (m, o) => exportJson(m, o) },
]);

/**
 * Export one model (or, for formats that hold several, an array of models).
 * @returns {{ data: Uint8Array|string, fileName: string, mime: string }}
 */
export function exportModels(formatId, models, opts = {}) {
  const ex = EXPORTERS.find((e) => e.id === formatId);
  if (!ex) throw new RangeError(`export: unknown format "${formatId}" (available: ${EXPORTERS.map((e) => e.id).join(', ')}).`);
  const list = Array.isArray(models) ? models : [models];
  if (!list.length) throw new RangeError('export: no model given.');
  if (!ex.multi && list.length > 1) throw new RangeError(`export: ${ex.name} holds one model per file; ${list.length} were given.`);
  for (const m of list) if (!ex.kinds.includes(m.kind)) throw unsupported(ex.name, m, ex.kinds);
  const base = String(opts.fileName || list[0].name || 'geometry').replace(/[\\/:*?"<>|]+/g, '_').replace(/\.[A-Za-z0-9]+$/, '');
  return { data: ex.write(ex.multi ? list : list[0], opts), fileName: `${base}.${ex.extension}`, mime: ex.mime };
}
