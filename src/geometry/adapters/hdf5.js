// HDF5 containers through the vendored h5wasm build of the HDF5 library (see docs/vendor.md).
//
// HDF5 is a file system in a file; what the datasets mean is decided by the convention the writer followed.
// The convention is recognised from the object tree, never from the file extension:
//
//   root group /ENS_MAA                        MED (SALOME, Code_Aster, Gmsh)      → frame / surface / volume
//   root label "Root Node of HDF5 File",
//     a CGNSLibraryVersion or CGNSBase_t node  CGNS on HDF5                        → frame / surface / volume
//   datasets coordx | coord and connect<n>     Exodus II on NetCDF-4               → frame / surface / volume
//   group /meshes/<n> with nodes and faces     Fluent mesh / case (.msh.h5, .cas.h5) → surface (2-D: frame)
//   a variable with cf_role = mesh_topology    UGRID 2-D mesh on NetCDF-4          → surface
//   a 2-D variable on x/y or lon/lat axes      CF grid (GEBCO, ETOPO, EMODnet …)   → grid (points when uneven)
//   anything else                              generic: `opts.datasets` names what to read; without it the
//                                              tree is listed in `meta.tree` and nothing is imported
//
// CF grids. The variable is `opts.variable`, else the first of elevation, z, Band1, bathymetry, topo, depth,
// DEPTH, depth_average, altitude, height … that lies on recognisable axes (lon/lat, x/y, COLUMNS/LINES), else
// one with a height/depth standard name, else the first gridded variable. Leading dimensions (time, level) are
// taken at index 0. scale_factor, add_offset, _FillValue, missing_value and valid_range / valid_min /
// valid_max are applied. The variable is read in strips of rows and block-averaged on the fly to at most
// `opts.maxCells` nodes (default 1 500 000), so a 15-arc-second tile never has to exist uncompressed in
// memory. `positive = "down"` or a depth name flips the sign so that z is up.
//
// MED. Coordinates and connectivities are stored component by component ("non-interlaced"). The first
// computation step of the mesh is read. MED numbers the nodes of solid cells with the opposite handedness of
// the VTK / Abaqus order used in this application; the corner nodes are mapped as
//   TE4  [0, 2, 1, 3]              HE8  [0, 3, 2, 1, 4, 7, 6, 5]
//   PE6  [0, 2, 1, 3, 5, 4]        PY5  [0, 3, 2, 1, 4]
// (canonical node k = MED node map[k]; the same tables as the Code_Aster and Gmsh MED readers). SE2, TR3 and
// QU4 are identical in both conventions. Quadratic cells (SE3, TR6, TR7, QU8, QU9, T10, H20, H27, P15, P18,
// P13) are read on their corner nodes. Families become sets: every group a family belongs to is an element
// set (families of cells) or a node set (families of nodes). Polygons, polyhedra and structured MED meshes are
// reported, not read.
//
// CGNS. Every zone of every base is read. Unstructured zones: all element sections including MIXED; the linear
// and quadratic CGNS orders are those of this application except HEXA_20 / HEXA_27 (mid-edge nodes of the top
// face and of the vertical edges are exchanged). Structured zones become hexahedra (3-D) or quadrilaterals
// (2-D) on the vertex lattice. Sections become element sets and surface groups; boundary conditions given as
// element or vertex lists / ranges become element or node sets. Zones are placed side by side: nodes on zone
// interfaces are not merged. NGON_n / NFACE_n polyhedral sections are reported, not read.
//
// Fluent (common fluids format, written by Fluent 2020 R1 and later by default). The file holds the tables of
// the legacy .msh / .cas format: /meshes/1/nodes/coords/<k> (n × dimension reals, first node number in the
// minId attribute), /meshes/1/faces/nodes/<k>/nnodes and …/nodes (node count per face, node numbers),
// /meshes/1/faces/c0/<k> and …/c1/<k> (adjacent cells, each over its own minId – maxId range of face numbers;
// faces outside every c1 range have no second cell) and faces/zoneTopology (id, zoneType, minId, maxId and
// the names as one ";"-separated string). As for the legacy files (importers/fluent.js, which builds the
// result) only the boundary faces are returned, grouped by face zone; the first mesh of the file is read.
//
// The runtime is 4.2 MB (the WebAssembly binary is embedded in its JavaScript module) and is loaded on first
// use only, relative to this module. The bytes of the file are mounted in the library's in-memory file system
// without being copied; the library then reads the parts it needs into its own memory.

import { inputBytes, baseName, FeMesh, triangulate3D, finish } from '../importers/common.js';
import { parsePrj, epsgInfo, normaliseUnit } from '../units.js';
import { finishGrid, decimationFactor } from '../importers/grids.js';
import { feResult } from '../importers/fe.js';
import { Features } from '../importers/vector.js';
import { faceZone, addFace, fluentBoundary } from '../importers/fluent.js';

export const id = 'hdf5';
export const title = 'HDF5 containers: NetCDF-4 grids and UGRID meshes, MED, CGNS, Exodus II, Fluent';
export const formats = ['hdf5', 'netcdf4', 'med', 'cgns', 'exodus-hdf5'];
export const extensions = ['h5', 'hdf5', 'he5', 'nc', 'nc4', 'med', 'rmed', 'cgns', 'e', 'exo', 'ex2', 'exii', 'g'];
export const library = { name: 'h5wasm', version: '0.10.3', licence: 'NIST software licence; HDF5 library BSD-3-Clause' };
/** Download size of the runtime (two JavaScript modules, WebAssembly embedded) [bytes]. */
export const sizeBytes = 48561 + 4150134;
/** What to convert a file to so that it can be imported when this adapter cannot run. */
export const convertTo = 'NetCDF classic (nccopy -k classic in.nc out.nc), or for meshes Gmsh .msh, VTK .vtu or a Nastran / Abaqus deck';

const DEFAULT_MAX_CELLS = 1500000;
const TREE_LIMIT = 2000;
const BASE = new URL('../../../vendor/h5wasm/', import.meta.url);
const isNode = typeof process !== 'undefined' && !!(process.versions && process.versions.node) && typeof window === 'undefined';

/**
 * True when the library can run here: WebAssembly is supported and, under Node (where it can be checked
 * without a request), the vendored files are on disk.
 */
export function available() {
  if (typeof WebAssembly !== 'object' || typeof WebAssembly.instantiate !== 'function') return false;
  if (!isNode || typeof process.getBuiltinModule !== 'function') return true;
  try { const fs = process.getBuiltinModule('node:fs'); return ['hdf5_hl.js', 'hdf5_util.js'].every((f) => fs.existsSync(new URL(f, BASE))); } catch { return true; }
}

let loading = null;
let serial = 0;

/**
 * Load and initialise the library (once; later calls return the same instance).
 * @returns {Promise<{ h5: object, FS: object }>} the h5wasm module (File, Group, Dataset) and its in-memory
 *   file system
 */
export function load() {
  if (!loading) {
    loading = (async () => {
      const h5 = await import(new URL('hdf5_hl.js', BASE).href);
      const mod = await h5.ready;
      // library errors become exceptions carrying the HDF5 error stack instead of text on the console
      mod.activate_throwing_error_handler();
      return { h5, FS: mod.FS };
    })().catch((err) => { loading = null; throw new Error(`The HDF5 library could not be started: ${err && err.message ? err.message : err}. HDF5 import needs the files in vendor/h5wasm/ and WebAssembly support; otherwise convert the file to ${convertTo}.`); });
  }
  return loading;
}

const SIGNATURE = [0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a];

/** True when the content carries the HDF5 signature (at 0 or after a user block of 512, 1024, 2048 … bytes). */
export function sniff(input) {
  const b = inputBytes(input);
  for (let off = 0; off + 8 <= b.length; off = off ? off * 2 : 512) {
    let ok = true;
    for (let i = 0; i < 8; i++) if (b[off + i] !== SIGNATURE[i]) { ok = false; break; }
    if (ok) return true;
    if (off >= 1 << 24) break;
  }
  return false;
}

// ---------------------------------------------------------------------------------------------------------
// Access helpers
// ---------------------------------------------------------------------------------------------------------

const isBig = (v) => v instanceof BigInt64Array || v instanceof BigUint64Array;
const isGroup = (o) => !!o && o.type === 'Group';
const isDataset = (o) => !!o && o.type === 'Dataset';
const trimText = (s) => String(s).replace(/\u0000[\s\S]*$/, '').trim();
const charsToText = (a) => { let s = ''; for (let i = 0; i < a.length; i++) { const c = a[i] & 0xff; if (!c) break; s += String.fromCharCode(c); } return s.trim(); };

/** The innermost line of an HDF5 error stack, which names the actual cause. */
function reason(err) {
  const msg = err && err.message ? err.message : String(err);
  const lines = msg.split('\n').map((l) => l.trim()).filter((l) => /^#\d+:/.test(l));
  if (!lines.length) return msg.split('\n')[0].slice(0, 300);
  return lines[lines.length - 1].replace(/^#\d+:\s*\S+\s+line\s+\d+\s+in\s+/, '').slice(0, 300);
}

/** Attribute value as plain JavaScript: number, string, number[] or string[]; undefined for references. */
function plain(v) {
  if (typeof v === 'bigint') return Number(v);
  if (typeof v === 'string') return trimText(v);
  if (typeof v === 'number' || typeof v === 'boolean') return v;
  if (ArrayBuffer.isView(v)) { const a = Array.from(v, Number); return a.length === 1 ? a[0] : a; }
  if (Array.isArray(v)) {
    if (!v.every((s) => typeof s === 'string')) return undefined;
    return v.length === 1 ? trimText(v[0]) : v.map(trimText);
  }
  return undefined;
}

const SKIP_ATTRS = new Set(['DIMENSION_LIST', 'REFERENCE_LIST', '_Netcdf4Coordinates', '_Netcdf4Dimid', '_NCProperties', '_nc3_strict']);

function attrsOf(obj) {
  const out = {};
  let all;
  try { all = obj.attrs; } catch { return out; }
  for (const k in all) {
    if (SKIP_ATTRS.has(k)) continue;
    try { const v = plain(all[k].value); if (v !== undefined) out[k] = v; } catch { /* an attribute of a type the library cannot convert */ }
  }
  return out;
}

function child(group, name) {
  try { return group.get(name); } catch { return null; }
}

function keysOf(group) {
  try { return group.keys(); } catch { return []; }
}

/** Numeric content of a dataset (or of a hyperslab of it) as a typed array; 64-bit integers become doubles. */
function numbers(ds, ranges) {
  let v = ranges ? ds.slice(ranges) : ds.value;
  if (v === null || v === undefined) return new Float64Array(0);
  if (typeof v === 'number') return Float64Array.of(v);
  if (typeof v === 'bigint') return Float64Array.of(Number(v));
  if (Array.isArray(v)) v = v.flat(Infinity);
  if (isBig(v) || Array.isArray(v)) { const o = new Float64Array(v.length); for (let i = 0; i < v.length; i++) o[i] = Number(v[i]); return o; }
  if (ArrayBuffer.isView(v)) return v;
  return new Float64Array(0);
}

/** Text of a character dataset: C1 / int8 arrays, fixed strings, or rows of single characters. */
function textOf(ds) {
  let v;
  try { v = ds.value; } catch { return ''; }
  if (typeof v === 'string') return trimText(v);
  if (ArrayBuffer.isView(v)) return charsToText(v);
  if (Array.isArray(v)) return trimText(v.flat(Infinity).map((c) => (typeof c === 'number' ? String.fromCharCode(c & 0xff) : String(c))).join(''));
  return '';
}

/** Rows of a 2-D character dataset (NetCDF char arrays) or of a 1-D string dataset, trimmed. */
function textRows(ds) {
  let v;
  try { v = ds.value; } catch { return []; }
  const shape = ds.shape || [];
  if (typeof v === 'string') return [trimText(v)];
  if (ArrayBuffer.isView(v) && shape.length === 2) { const out = []; for (let i = 0; i < shape[0]; i++) out.push(charsToText(v.subarray(i * shape[1], (i + 1) * shape[1]))); return out; }
  if (!Array.isArray(v)) return [];
  if (shape.length === 2 && v.length === shape[0] * shape[1]) {
    const out = [];
    for (let i = 0; i < shape[0]; i++) { let s = ''; for (let k = 0; k < shape[1]; k++) { const c = v[i * shape[1] + k]; if (c === '' || c === '\u0000') break; s += c; } out.push(s.trim()); }
    return out;
  }
  return v.map((s) => (Array.isArray(s) ? charsToText(s) : trimText(s)));
}

const dtypeLabel = (d) => (typeof d === 'string' ? d : d && d.compound_type ? 'compound' : d && d.array_type ? 'array' : 'other');
const isNumericType = (d) => typeof d === 'string' && /^[<>|=]?[bBhHiIlLqQefd]\d*$/.test(d);

/**
 * List the objects of the file: `[{ path, type: 'group' | 'dataset' | 'link' | 'datatype', shape?, dtype? }]`,
 * depth first, at most TREE_LIMIT entries.
 */
function listTree(file) {
  const out = [];
  let truncated = false;
  const visit = (group, prefix, depth) => {
    for (const k of keysOf(group)) {
      if (out.length >= TREE_LIMIT) { truncated = true; return; }
      const path = `${prefix}/${k}`;
      const o = child(group, k);
      if (isDataset(o)) {
        let shape = null, dtype = 'unknown';
        try { shape = o.shape; dtype = dtypeLabel(o.dtype); } catch { /* listed without its description */ }
        out.push({ path, type: 'dataset', shape: shape || [], dtype });
      } else if (isGroup(o)) {
        out.push({ path, type: 'group' });
        // soft links can point back up the tree: the depth limit ends such loops
        if (depth < 12) visit(o, path, depth + 1); else truncated = true;
      } else out.push({ path, type: o && o.type === 'Datatype' ? 'datatype' : 'link' });
    }
  };
  visit(file, '', 0);
  return { entries: out, truncated };
}

// ---------------------------------------------------------------------------------------------------------
// NetCDF-4: variables, CF grids, UGRID
// ---------------------------------------------------------------------------------------------------------

/** Every dataset of the file as a NetCDF-style variable `{ name, path, group, ds, shape, dtype, atts, dimOnly }`. */
function collectVariables(file, tree) {
  const vars = [];
  for (const e of tree.entries) {
    if (e.type !== 'dataset') continue;
    const ds = child(file, e.path);
    if (!isDataset(ds)) continue;
    const atts = attrsOf(ds);
    const cut = e.path.lastIndexOf('/');
    vars.push({
      name: e.path.slice(cut + 1), path: e.path, group: e.path.slice(0, cut), ds, shape: e.shape, dtype: e.dtype, atts,
      // a dimension that has no coordinate variable is stored as an empty scale carrying this marker
      dimOnly: atts.CLASS === 'DIMENSION_SCALE' && /^This is a netCDF dimension but not a netCDF variable/.test(String(atts.NAME || '')),
    });
  }
  return vars;
}

const isX = (v) => !!v && (v.atts.axis === 'X' || /^(longitude|projection_x_coordinate|grid_longitude)$/.test(v.atts.standard_name || '') || /degrees?_?e/i.test(v.atts.units || '') || /^(x|lon|long|longitude|easting|xc|eastings?|columns?)$/i.test(v.name));
const isY = (v) => !!v && (v.atts.axis === 'Y' || /^(latitude|projection_y_coordinate|grid_latitude)$/.test(v.atts.standard_name || '') || /degrees?_?n/i.test(v.atts.units || '') || /^(y|lat|latitude|northing|yc|northings?|lines?|rows?)$/i.test(v.name));
const isDegrees = (v) => !!v && (/degree/i.test(v.atts.units || '') || /^(longitude|latitude)$/.test(v.atts.standard_name || '') || /^(lon|lat|longitude|latitude)$/i.test(v.name));
// in order of preference: GEBCO / EMODnet "elevation", ETOPO and GMT "z", GDAL "Band1", then the rest
const Z_NAMES = ['elevation', 'z', 'band1', 'bathymetry', 'bathy', 'topo', 'topography', 'depth', 'depth_average', 'depth_smooth', 'deptho', 'altitude', 'height', 'elev', 'dem', 'h'];
const DEPTH_NAME = /^(depth|deptho|depth_average|depth_smooth|depth_min|depth_max)$/i;

function gridCrs(vars, v, x, y) {
  const gm = v.atts.grid_mapping ? vars.find((q) => q.name === v.atts.grid_mapping) : vars.find((q) => q.atts.grid_mapping_name);
  if (gm) {
    const wkt = gm.atts.crs_wkt || gm.atts.spatial_ref;
    if (typeof wkt === 'string') { const info = parsePrj(wkt); if (info.kind !== 'unknown') return info; }
    const code = gm.atts.epsg_code || gm.atts.EPSG_code;
    const m = code ? /(\d{4,5})/.exec(String(code)) : null;
    const info = m ? epsgInfo(Number(m[1])) : null;
    if (info) return { ...info, unitScale: 1 };
    if (gm.atts.grid_mapping_name === 'latitude_longitude') return { kind: 'geographic', epsg: 4326, name: 'WGS 84', unit: 'deg' };
  }
  if (isDegrees(x) && isDegrees(y)) return { kind: 'geographic', epsg: 4326, name: 'WGS 84 (assumed from degree units)', unit: 'deg' };
  const u = x && normaliseUnit(x.atts.units);
  if (u) return { kind: 'projected', epsg: null, name: '', unit: u, unitScale: { m: 1, km: 1000, ft: 0.3048, 'us-ft': 1200 / 3937, mm: 0.001, cm: 0.01 }[u] || 1 };
  return null;
}

function regular(c) {
  if (c.length < 2) return null;
  const d = (c[c.length - 1] - c[0]) / (c.length - 1);
  if (d === 0) return null;
  for (let i = 1; i < c.length; i++) if (Math.abs(c[i] - c[i - 1] - d) > 1e-4 * Math.abs(d)) return null;
  return d;
}

/** Coordinate variable of dimension `i` of `v`: the attached dimension scale, else a 1-D variable that fits. */
function axisOf(vars, byPath, v, i, test) {
  try {
    const s = v.ds.get_attached_scales(i);
    if (s && s.length) { const a = byPath.get(s[0]); if (a && !a.dimOnly && a.shape.length === 1) return a; if (a && a.dimOnly) return null; }
  } catch { /* no dimension scales: plain HDF5 written with CF attributes */ }
  const fits = (q) => q !== v && !q.dimOnly && q.shape.length === 1 && q.shape[0] === v.shape[i] && isNumericType(q.dtype) && test(q);
  return vars.find((q) => q.group === v.group && fits(q)) || vars.find(fits) || null;
}

function gridCandidates(vars) {
  const byPath = new Map(vars.map((v) => [v.path, v]));
  const cands = [];
  for (const v of vars) {
    const r = v.shape.length;
    if (r < 2 || v.dimOnly || !isNumericType(v.dtype)) continue;
    const any = () => true;
    let a = axisOf(vars, byPath, v, r - 2, (q) => isY(q) || isX(q)), b = axisOf(vars, byPath, v, r - 1, (q) => isX(q) || isY(q));
    let yx = null;
    if (isY(a) && isX(b)) yx = true; else if (isX(a) && isY(b)) yx = false;
    else {
      // attached scales with names that say nothing (dimension "n", "m"): rows are taken as y
      a = axisOf(vars, byPath, v, r - 2, any); b = axisOf(vars, byPath, v, r - 1, any);
      let attached = false;
      try { attached = v.ds.get_attached_scales(r - 2).length > 0 && v.ds.get_attached_scales(r - 1).length > 0; } catch { attached = false; }
      if (attached && a && b && a !== b && r === 2 && !isX(a) && !isY(b)) yx = true;
    }
    if (yx === null || !a || !b || a === b) continue;
    cands.push({ v, y: yx ? a : b, x: yx ? b : a, yx });
  }
  return cands;
}

function cfGrid(ctx, vars, cands) {
  const { opts, fname, name, warnings } = ctx;
  let pick = null;
  if (opts.variable) {
    pick = cands.find((c) => c.v.name === opts.variable || c.v.path === opts.variable);
    if (!pick) throw new Error(`${fname}: variable "${opts.variable}" is not a gridded variable on x/y axes. Gridded variables: ${cands.map((c) => c.v.name).join(', ') || 'none'}.`);
  } else {
    for (const zn of Z_NAMES) { pick = cands.find((c) => c.v.name.toLowerCase() === zn); if (pick) break; }
    if (!pick) pick = cands.find((c) => /altitude|depth|height|elevation/.test(c.v.atts.standard_name || '')) || cands[0];
  }
  const { v, x, y, yx } = pick;
  if (cands.length > 1 && !opts.variable) warnings.push(`Several gridded variables are present (${cands.map((c) => c.v.name).join(', ')}); "${v.name}" was used. Choose another with opts.variable.`);
  const r = v.shape.length;
  const n1 = v.shape[r - 2], n2 = v.shape[r - 1];
  const lead = v.shape.slice(0, -2);
  if (lead.some((l) => l === 0)) throw new Error(`${fname}: variable "${v.name}" has no records.`);
  if (lead.some((l) => l > 1)) warnings.push(`"${v.name}" has leading dimension(s) of size ${lead.join(' × ')}; the first index of each was used.`);
  const xs = numbers(x.ds), ys = numbers(y.ds);
  const nx = xs.length, ny = ys.length;
  if ((yx ? n2 : n1) !== nx || (yx ? n1 : n2) !== ny) throw new Error(`${fname}: variable "${v.name}" is ${v.shape.join(' × ')} but its axes "${y.name}" and "${x.name}" have ${ny} and ${nx} values.`);
  const dx = regular(xs), dy = regular(ys);
  const crs = gridCrs(vars, v, x, y);
  const depthLike = v.atts.positive === 'down' || /depth/.test(v.atts.standard_name || '') || DEPTH_NAME.test(v.name);
  const fill = v.atts._FillValue, miss = v.atts.missing_value;
  const vr = Array.isArray(v.atts.valid_range) ? v.atts.valid_range : [v.atts.valid_min ?? -Infinity, v.atts.valid_max ?? Infinity];
  const sf = typeof v.atts.scale_factor === 'number' ? v.atts.scale_factor : 1, ao = typeof v.atts.add_offset === 'number' ? v.atts.add_offset : 0;
  const maxCells = opts.maxCells ?? DEFAULT_MAX_CELLS;
  const extra = { flavour: 'netcdf4-cf-grid', convention: ctx.rootAtts.Conventions || null, variable: v.name, longName: v.atts.long_name || null, standardName: v.atts.standard_name || null, packing: sf !== 1 || ao !== 0 ? { scaleFactor: sf, addOffset: ao } : null, axes: { x: x.name, y: y.name }, tree: ctx.tree };
  const leadSel = lead.map(() => [0, 1]);
  const rowsPerStrip = Math.max(1, Math.floor(2000000 / Math.max(1, n2)));
  /** Walk the variable strip by strip: cb(value in physical units, ix, iy) for every valid cell kept. */
  const scan = (stepA, stepB, cb) => {
    for (let a0 = 0; a0 < n1; a0 += rowsPerStrip * stepA) {
      const a1 = Math.min(n1, a0 + rowsPerStrip * stepA);
      let raw;
      try { raw = numbers(v.ds, [...leadSel, [a0, a1, stepA], [0, n2, stepB]]); } catch (err) { throw new Error(`${fname}: variable "${v.name}" could not be read (${reason(err)}). The file is probably truncated, or the variable is compressed with a filter that needs a plug-in (deflate and shuffle are built in; Blosc, Zstandard and the like are not).`); }
      const nb = Math.ceil(n2 / stepB);
      if (raw.length !== Math.ceil((a1 - a0) / stepA) * nb) throw new Error(`${fname}: variable "${v.name}" returned ${raw.length} values for a strip of ${Math.ceil((a1 - a0) / stepA)} × ${nb}.`);
      for (let i = 0, a = a0; a < a1; a += stepA) {
        for (let b = 0; b < n2; b += stepB, i++) {
          const q = raw[i];
          if (q === fill || q === miss || q < vr[0] || q > vr[1] || q !== q) continue;
          cb(q * sf + ao, yx ? b : a, yx ? a : b);
        }
      }
    }
  };
  if (dx === null || dy === null) {
    // unevenly spaced axes cannot be a grid model: return the nodes as soundings, thinned to maxCells
    const s = decimationFactor(nx, ny, maxCells);
    const F = new Features(name);
    scan(s, s, (q, ix, iy) => F.addPoint(xs[ix], ys[iy], depthLike && !opts.vertical ? -q : q, null));
    warnings.push(`The coordinate axes are not evenly spaced, so the grid nodes are returned as points${s > 1 ? ` (every ${s}th row and column, to stay under maxCells)` : ''}.`);
    return F.build(opts, { format: 'netcdf4', fileName: fname, crs: crs ? { ...crs, source: crs.epsg ? `EPSG:${crs.epsg}` : crs.name } : null, zScale: 1, warnings, expected: 'expected valid cells', extra });
  }
  const f = decimationFactor(nx, ny, maxCells);
  const mx = Math.ceil(nx / f), my = Math.ceil(ny / f);
  const sum = new Float64Array(mx * my), cnt = new Uint32Array(mx * my);
  // output node (ix, iy) counts from the west and from the south whatever the storage order of the file
  const colOf = new Uint32Array(nx), rowOf = new Uint32Array(ny);
  for (let i = 0; i < nx; i++) colOf[i] = Math.floor((dx > 0 ? i : nx - 1 - i) / f);
  for (let i = 0; i < ny; i++) rowOf[i] = Math.floor((dy > 0 ? i : ny - 1 - i) / f) * mx;
  scan(1, 1, (q, ix, iy) => { const k = rowOf[iy] + colOf[ix]; sum[k] += q; cnt[k]++; });
  const z = new Float32Array(mx * my);
  for (let k = 0; k < z.length; k++) z[k] = cnt[k] ? sum[k] / cnt[k] : NaN;
  const ax = Math.abs(dx), ay = Math.abs(dy);
  return finishGrid({ name, nx: mx, ny: my, x0: Math.min(xs[0], xs[nx - 1]) + ((f - 1) / 2) * ax, y0: Math.min(ys[0], ys[ny - 1]) + ((f - 1) / 2) * ay, dx: ax * f, dy: ay * f, z }, { ...opts, maxCells: 0 }, {
    format: 'netcdf4', crs, zName: depthLike ? 'depth' : 'elevation', zUnit: normaliseUnit(v.atts.units) ? v.atts.units : (crs && crs.kind === 'geographic' ? 'm' : null), warnings, decimated: f, sourceSize: [nx, ny], extra,
  });
}

function ugrid(ctx, vars, topo) {
  const { opts, fname, name } = ctx;
  const V = (nm) => vars.find((v) => v.name === nm && v.group === topo.group) || vars.find((v) => v.name === nm);
  const [xn, yn] = String(topo.atts.node_coordinates || '').trim().split(/\s+/);
  const xv = V(xn), yv = V(yn), fv = V(topo.atts.face_node_connectivity);
  if (!xv || !yv || !fv || fv.shape.length !== 2) throw new Error(`${fname}: the UGRID mesh topology "${topo.name}" does not name node coordinates and a face_node_connectivity variable that exist in the file.`);
  const x = numbers(xv.ds), y = numbers(yv.ds);
  const n = x.length;
  const zv = vars.find((v) => v.shape.length === 1 && v.shape[0] === n && v !== xv && v !== yv && !v.dimOnly && isNumericType(v.dtype) && (/altitude|depth|elevation|height/.test(v.atts.standard_name || '') || /node_z|^z$|bath|depth|elev/i.test(v.name)));
  const z = zv ? numbers(zv.ds) : new Float64Array(n);
  const down = zv && (zv.atts.positive === 'down' || /depth/.test(zv.atts.standard_name || ''));
  const fill = zv ? zv.atts._FillValue : undefined;
  const coords = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) { coords[3 * i] = x[i]; coords[3 * i + 1] = y[i]; coords[3 * i + 2] = z[i] === fill || z[i] !== z[i] ? 0 : (down ? -z[i] : z[i]); }
  const conn = numbers(fv.ds);
  const [nf, mxn] = fv.shape;
  const start = typeof fv.atts.start_index === 'number' ? fv.atts.start_index : 0, ffill = fv.atts._FillValue;
  const tris = [];
  const hist = {};
  for (let f = 0; f < nf; f++) {
    const ids = [];
    for (let k = 0; k < mxn; k++) { const q = conn[f * mxn + k]; if (q === ffill || q < start || q !== q) break; ids.push(q - start); }
    if (ids.length < 3 || ids.some((k) => k >= n)) continue;
    const kind = ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon';
    hist[kind] = (hist[kind] || 0) + 1;
    if (ids.length === 3) tris.push(ids[0], ids[1], ids[2]);
    else if (ids.length === 4) tris.push(ids[0], ids[1], ids[2], ids[0], ids[2], ids[3]);
    else { const flat = new Float64Array(ids.length * 3); ids.forEach((k, i) => { flat[3 * i] = coords[3 * k]; flat[3 * i + 1] = coords[3 * k + 1]; flat[3 * i + 2] = 0; }); for (const k of triangulate3D(flat)) tris.push(ids[k]); }
  }
  if (!tris.length) throw new Error(`${fname}: the UGRID mesh "${topo.name}" has no faces with three or more nodes.`);
  const F = new Features(name);
  F.addMesh(coords, Uint32Array.from(tris), topo.name);
  const crs = gridCrs(vars, zv || xv, xv, yv);
  const r = F.build(opts, { format: 'netcdf4', fileName: fname, crs: crs ? { ...crs, source: crs.epsg ? `EPSG:${crs.epsg}` : crs.name } : null, zScale: 1, warnings: ctx.warnings, expected: '', extra: { flavour: 'netcdf4-ugrid', convention: 'UGRID', mesh: topo.name, nodeElevation: zv ? zv.name : null, tree: ctx.tree } });
  r.meta.elementTypes = hist;
  if (!zv) r.warnings.push('The UGRID mesh has no node elevation variable; z is 0.');
  return r;
}

// ---------------------------------------------------------------------------------------------------------
// Exodus II on NetCDF-4
// ---------------------------------------------------------------------------------------------------------

const HEX20_EDGES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15];

function exodus(ctx) {
  const { file, opts, fname, name } = ctx;
  const fe = new FeMesh(name);
  const D = (nm) => { const o = child(file, nm); return isDataset(o) ? o : null; };
  const cx = D('coordx'), cy = D('coordy'), cz = D('coordz'), c = D('coord');
  let X, Y, Z, nNodes;
  if (!cx && !(c && (c.shape || []).length === 2)) throw new Error(`${fname}: the Exodus file has no coordx / coord dataset.`);
  if (cx) { X = numbers(cx); nNodes = X.length; Y = cy ? numbers(cy) : new Float64Array(nNodes); Z = cz ? numbers(cz) : new Float64Array(nNodes); } else {
    const all = numbers(c);
    const nd = c.shape[0];
    nNodes = c.shape[1];
    X = all.subarray(0, nNodes); Y = nd > 1 ? all.subarray(nNodes, 2 * nNodes) : new Float64Array(nNodes); Z = nd > 2 ? all.subarray(2 * nNodes, 3 * nNodes) : new Float64Array(nNodes);
  }
  for (let i = 0; i < nNodes; i++) fe.addNode(i + 1, X[i], Y[i], Z[i]);
  const ids = D('eb_prop1') ? numbers(D('eb_prop1')) : [];
  const names = D('eb_names') ? textRows(D('eb_names')) : [];
  let eid = 0, blk = 0;
  const usedTags = new Set();
  const blocks = keysOf(file).map((k) => /^connect(\d+)$/.exec(k)).filter(Boolean).sort((a, b) => a[1] - b[1]);
  for (const m of blocks) {
    const ds = D(m[0]);
    if (!ds || (ds.shape || []).length !== 2) continue;
    const k = Number(m[1]);
    const [ne, npe] = ds.shape;
    const et = String(attrsOf(ds).elem_type || '').toUpperCase();
    let type = null, order = null;
    if (/^(BAR|BEAM|TRUSS|EDGE|ROD)/.test(et)) type = npe === 3 ? 'line3' : npe === 2 ? 'line2' : null;
    else if (/^TRI/.test(et)) type = npe === 6 ? 'tri6' : npe === 3 ? 'tri3' : null;
    else if (/^(QUAD|SHELL)/.test(et)) type = npe === 9 ? 'quad9' : npe === 8 ? 'quad8' : npe === 4 ? 'quad4' : null;
    else if (/^TET/.test(et)) type = npe === 10 ? 'tet10' : npe === 4 ? 'tet4' : null;
    else if (/^HEX/.test(et)) { type = npe === 27 ? 'hex27' : npe === 20 ? 'hex20' : npe === 8 ? 'hex8' : null; if (npe === 20) order = HEX20_EDGES; if (npe === 27) order = HEX20_EDGES.concat([20, 21, 22, 23, 24, 25, 26]); } else if (/^WEDGE/.test(et)) type = npe === 15 ? 'wedge15' : npe === 6 ? 'wedge6' : null;
    else if (/^PYR/.test(et)) type = npe === 13 ? 'pyr13' : npe === 5 ? 'pyr5' : null;
    // block ids are positive and unique in a valid file; a writer that numbers from 0 would merge two blocks
    let tag = ids[k - 1] > 0 && !usedTags.has(ids[k - 1]) ? ids[k - 1] : k;
    while (usedTags.has(tag)) tag++;
    usedTags.add(tag);
    const label = names[k - 1] || `block ${tag}`;
    fe.setProp(tag, { name: label });
    blk++;
    if (!type) { fe.skip(`element type ${et || '?'} with ${npe} nodes`, ne); eid += ne; continue; }
    const conn = numbers(ds);
    const nodes = new Array(npe);
    const set = [];
    for (let e = 0; e < ne; e++) { for (let q = 0; q < npe; q++) nodes[q] = conn[e * npe + (order ? order[q] : q)]; fe.addElement(type, nodes, ++eid, tag); set.push(eid); }
    fe.addElemSet(label, set);
  }
  const nsNames = D('ns_names') ? textRows(D('ns_names')) : [], nsIds = D('ns_prop1') ? numbers(D('ns_prop1')) : [];
  let sideSets = 0;
  for (const k of keysOf(file)) {
    const m = /^node_ns(\d+)$/.exec(k);
    if (m && D(k)) fe.addNodeSet(nsNames[m[1] - 1] || `nodeset ${nsIds[m[1] - 1] || m[1]}`, numbers(D(k)));
    if (/^elem_ss\d+$/.test(k)) sideSets++;
  }
  const warnings = ctx.warnings;
  if (sideSets) warnings.push(`${sideSets} side set(s) were not read (node sets and element blocks are).`);
  return feResult(fe, 'exodus-hdf5', opts, { name: fname, expected: 'expected a coordx or coord dataset with nodes', warnings, skipped: sideSets ? { sideSets } : {}, extra: { flavour: 'exodus-netcdf4', title: ctx.rootAtts.title || null, blocks: blk, tree: ctx.tree } });
}

// ---------------------------------------------------------------------------------------------------------
// MED
// ---------------------------------------------------------------------------------------------------------

/** MED geometric type → [canonical type, canonical node k = MED node map[k]]; quadratic types on their corners. */
const MED_TYPES = {
  SE2: ['line2', [0, 1]], SE3: ['line2', [0, 1]],
  TR3: ['tri3', [0, 1, 2]], TR6: ['tri3', [0, 1, 2]], TR7: ['tri3', [0, 1, 2]],
  QU4: ['quad4', [0, 1, 2, 3]], QU8: ['quad4', [0, 1, 2, 3]], QU9: ['quad4', [0, 1, 2, 3]],
  TE4: ['tet4', [0, 2, 1, 3]], T10: ['tet4', [0, 2, 1, 3]],
  HE8: ['hex8', [0, 3, 2, 1, 4, 7, 6, 5]], H20: ['hex8', [0, 3, 2, 1, 4, 7, 6, 5]], H27: ['hex8', [0, 3, 2, 1, 4, 7, 6, 5]],
  PE6: ['wedge6', [0, 2, 1, 3, 5, 4]], P15: ['wedge6', [0, 2, 1, 3, 5, 4]], P18: ['wedge6', [0, 2, 1, 3, 5, 4]],
  PY5: ['pyr5', [0, 3, 2, 1, 4]], P13: ['pyr5', [0, 3, 2, 1, 4]],
};
const MED_NODES = { SE2: 2, SE3: 3, TR3: 3, TR6: 6, TR7: 7, QU4: 4, QU8: 8, QU9: 9, TE4: 4, T10: 10, HE8: 8, H20: 20, H27: 27, PE6: 6, P15: 15, P18: 18, PY5: 5, P13: 13, PO1: 1 };
const MED_QUADRATIC = new Set(['SE3', 'TR6', 'TR7', 'QU8', 'QU9', 'T10', 'H20', 'H27', 'P15', 'P18', 'P13']);

/** Family number → { name, groups } from /FAS/<mesh>/{ELEME,NOEUD}. */
function medFamilies(file, meshName) {
  const out = new Map();
  const fas = child(file, `FAS/${meshName}`);
  if (!isGroup(fas)) return out;
  for (const kind of ['ELEME', 'NOEUD']) {
    const g = child(fas, kind);
    if (!isGroup(g)) continue;
    for (const fam of keysOf(g)) {
      const fg = child(g, fam);
      if (!isGroup(fg)) continue;
      const num = attrsOf(fg).NUM;
      if (typeof num !== 'number') continue;
      const groups = [];
      const nom = child(fg, 'GRO/NOM');
      if (isDataset(nom)) {
        let v = null;
        try { v = nom.value; } catch { v = null; }
        // group names are 80-character fields, stored as int8[80] items or as fixed strings
        if (ArrayBuffer.isView(v)) for (let i = 0; i + 80 <= v.length; i += 80) groups.push(charsToText(v.subarray(i, i + 80)));
        else if (Array.isArray(v)) for (const row of v) groups.push(Array.isArray(row) || ArrayBuffer.isView(row) ? charsToText(row) : trimText(row));
        else if (typeof v === 'string') for (let i = 0; i < v.length; i += 80) groups.push(v.slice(i, i + 80).trim());
      }
      out.set(num, { name: fam, groups: groups.filter(Boolean) });
    }
  }
  return out;
}

function med(ctx) {
  const { file, opts, fname, warnings } = ctx;
  const ens = child(file, 'ENS_MAA');
  const meshNames = keysOf(ens).filter((k) => isGroup(child(ens, k)));
  if (!meshNames.length) throw new Error(`${fname}: the MED file has no mesh under /ENS_MAA (it may hold fields or profiles only).`);
  let meshName = meshNames[0];
  if (opts.mesh) { if (!meshNames.includes(opts.mesh)) throw new Error(`${fname}: MED mesh "${opts.mesh}" is not in the file. Meshes: ${meshNames.join(', ')}.`); meshName = opts.mesh; } else if (meshNames.length > 1) warnings.push(`The file holds ${meshNames.length} meshes (${meshNames.join(', ')}); "${meshName}" was read. Choose another with opts.mesh.`);
  const mg = child(ens, meshName);
  const ma = attrsOf(mg);
  if (ma.TYP === 1) throw new Error(`${fname}: MED mesh "${meshName}" is a structured grid, which is not read. Convert it to an unstructured mesh in SALOME (or export VTK / UNV).`);
  // MED 3 and 4 keep each computation step in its own group; older files have NOE and MAI directly under the mesh
  let step = mg, stepName = null;
  if (!isGroup(child(mg, 'NOE'))) {
    const steps = keysOf(mg).filter((k) => isGroup(child(mg, `${k}/NOE`)));
    if (!steps.length) throw new Error(`${fname}: MED mesh "${meshName}" has no node group (NOE).`);
    stepName = steps[0];
    step = child(mg, stepName);
    if (steps.length > 1) warnings.push(`MED mesh "${meshName}" has ${steps.length} computation steps; the first was read.`);
  }
  const cooDs = child(step, 'NOE/COO');
  if (!isDataset(cooDs)) throw new Error(`${fname}: MED mesh "${meshName}" has no node coordinates (NOE/COO).`);
  const coo = numbers(cooDs);
  const nbr = attrsOf(cooDs).NBR;
  const dim = typeof ma.ESP === 'number' && ma.ESP >= 1 && ma.ESP <= 3 ? ma.ESP : typeof nbr === 'number' && nbr > 0 ? Math.round(coo.length / nbr) : 3;
  const nNodes = Math.floor(coo.length / dim);
  if (!(dim >= 1 && dim <= 3) || nNodes * dim !== coo.length) throw new Error(`${fname}: the MED coordinate array has ${coo.length} values, which is not a whole number of ${dim}-D nodes.`);
  const fe = new FeMesh(meshName || ctx.name);
  // all x, then all y, then all z; connectivities refer to nodes by position (1-based), not by their NUM label
  for (let i = 0; i < nNodes; i++) fe.addNode(i + 1, coo[i], dim > 1 ? coo[nNodes + i] : 0, dim > 2 ? coo[2 * nNodes + i] : 0);
  const families = medFamilies(file, meshName);
  const famLabel = (num) => { const f = families.get(num); return f ? (f.groups.length ? f.groups.join(' + ') : f.name) : `family ${num}`; };
  // families that belong to the same groups share one tag, so that a group is one surface group / cell tag
  const tagOf = new Map(), tagLabels = new Map();
  const tagFor = (fnum) => {
    if (!fnum) return 0;
    let t = tagOf.get(fnum);
    if (t === undefined) {
      const lab = famLabel(fnum);
      t = tagLabels.get(lab);
      if (t === undefined) { t = tagLabels.size + 1; tagLabels.set(lab, t); fe.setProp(t, { name: lab }); }
      tagOf.set(fnum, t);
    }
    return t;
  };
  const nodeFam = child(step, 'NOE/FAM');
  if (isDataset(nodeFam)) {
    const nf = numbers(nodeFam);
    const sets = new Map();
    for (let i = 0; i < Math.min(nf.length, nNodes); i++) {
      const fam = families.get(nf[i]);
      if (!nf[i] || !fam) continue;
      for (const g of fam.groups) { let s = sets.get(g); if (!s) { s = []; sets.set(g, s); } s.push(i + 1); }
    }
    for (const [g, s] of sets) fe.addNodeSet(g, s);
  }
  const mai = child(step, 'MAI');
  const hist = {}, quadratic = {};
  let eid = 0;
  const elemSets = new Map();
  for (const t of isGroup(mai) ? keysOf(mai) : []) {
    const tg = child(mai, t);
    const nodDs = isGroup(tg) ? child(tg, 'NOD') : null;
    if (!isDataset(nodDs)) continue;
    const def = MED_TYPES[t];
    const npe = MED_NODES[t];
    const total = (nodDs.shape || []).reduce((a, b) => a * b, 1);
    if (!def) {
      // polygons (POG) and polyhedra (POE) index their connectivity through INN
      const inn = child(tg, 'INN');
      const ne = npe ? Math.floor(total / npe) : isDataset(inn) ? Math.max(1, ((inn.shape || [1])[0] || 1) - 1) : 1;
      if (t !== 'PO1') fe.skip(`MED cell type ${t === 'POG' ? 'POG (polygons)' : t === 'POE' ? 'POE (polyhedra)' : t}`, ne);
      eid += ne;
      continue;
    }
    const nod = numbers(nodDs);
    const ne = Math.floor(nod.length / npe);
    if (ne * npe !== nod.length) throw new Error(`${fname}: the MED connectivity of ${t} has ${nod.length} entries, which is not a multiple of ${npe} nodes per cell.`);
    const famDs = child(tg, 'FAM');
    const fam = isDataset(famDs) ? numbers(famDs) : null;
    const [type, map] = def;
    const nodes = new Array(map.length);
    for (let e = 0; e < ne; e++) {
      // component k of every cell is stored before component k + 1
      for (let k = 0; k < map.length; k++) nodes[k] = nod[map[k] * ne + e];
      const fnum = fam && e < fam.length ? fam[e] : 0;
      fe.addElement(type, nodes, ++eid, tagFor(fnum));
      const f = fnum ? families.get(fnum) : null;
      if (f) for (const g of f.groups) { let s = elemSets.get(g); if (!s) { s = []; elemSets.set(g, s); } s.push(eid); }
    }
    hist[t] = ne;
    if (MED_QUADRATIC.has(t)) quadratic[t] = ne;
  }
  for (const [g, s] of elemSets) fe.addElemSet(g, s);
  const qk = Object.keys(quadratic);
  if (qk.length) warnings.push(`Quadratic MED cells (${qk.map((k) => `${quadratic[k]} ${k}`).join(', ')}) were read on their corner nodes; the mid-side nodes are not kept.`);
  const info = attrsOf(child(file, 'INFOS_GENERALES') || {});
  const version = typeof info.MAJ === 'number' ? `${info.MAJ}.${info.MIN ?? 0}.${info.REL ?? 0}` : null;
  // axis units are 16-character fields; a mesh in one length unit repeats it
  const uni = String(ma.UNI || '').match(/.{1,16}/g) || [];
  const unitHint = uni.map((u) => normaliseUnit(u.trim())).find(Boolean) || null;
  return feResult(fe, 'med', opts, {
    name: fname, unitHint, expected: 'expected node coordinates in NOE/COO', warnings,
    extra: { flavour: 'med', medVersion: version, mesh: meshName, meshes: meshNames, description: ma.DES || null, spaceDimension: dim, medTypes: hist, families: [...families].map(([num, f]) => ({ number: num, name: f.name, groups: f.groups })), tree: ctx.tree },
  });
}

// ---------------------------------------------------------------------------------------------------------
// CGNS
// ---------------------------------------------------------------------------------------------------------

/** CGNS ElementType_t → [canonical type, nodes per element, node order]. */
const CGNS_TYPES = {
  3: ['line2', 2], 4: ['line3', 3], 5: ['tri3', 3], 6: ['tri6', 6], 7: ['quad4', 4], 8: ['quad8', 8], 9: ['quad9', 9],
  10: ['tet4', 4], 11: ['tet10', 10], 12: ['pyr5', 5], 13: ['pyr14', 14], 14: ['wedge6', 6], 15: ['wedge15', 15], 16: ['wedge18', 18],
  17: ['hex8', 8], 18: ['hex20', 20, HEX20_EDGES], 19: ['hex27', 27, HEX20_EDGES.concat([20, 21, 22, 23, 24, 25, 26])], 21: ['pyr13', 13],
};
const CGNS_NAMES = { 2: 'NODE', 20: 'MIXED', 22: 'NGON_n', 23: 'NFACE_n' };
const CGNS_UNITS = { meter: 'm', centimeter: 'cm', millimeter: 'mm', foot: 'ft', inch: 'in' };

function cgns(ctx) {
  const { file, opts, fname, warnings } = ctx;
  const label = (g) => attrsOf(g).label || '';
  const kids = (g, want) => keysOf(g).map((k) => [k, child(g, k)]).filter(([, o]) => isGroup(o) && (!want || label(o) === want));
  const data = (g) => { const d = isGroup(g) ? child(g, ' data') : null; return isDataset(d) ? d : null; };
  const nums = (g) => { const d = data(g); return d ? numbers(d) : new Float64Array(0); };
  const text = (g) => { const d = data(g); return d ? textOf(d) : ''; };
  const bases = kids(file, 'CGNSBase_t');
  if (!bases.length) throw new Error(`${fname}: the CGNS file has no base (CGNSBase_t node).`);
  const fe = new FeMesh(ctx.name);
  let version = null;
  const lv = child(file, 'CGNSLibraryVersion');
  if (isGroup(lv)) { const q = nums(lv); if (q.length) version = Number(q[0].toFixed(2)); }
  let nodeBase = 0, elemBase = 0, top = 0, tag = 0, structured = 0, zonesRead = 0, bcSkipped = 0, unitHint = null;
  const zoneInfo = [];
  const sectionTypes = {};
  for (const [baseName_, base] of bases) {
    const bd = nums(base);
    // a family may list the names it is known by (FamilyName_t children); sets are filed under those too
    const alias = new Map();
    for (const [famName, famNode] of kids(base, 'Family_t')) { const names = kids(famNode, 'FamilyName_t').map(([k, o]) => text(o) || k).filter((n) => n && n !== famName); if (names.length) alias.set(famName, names); }
    const physDim = bd.length > 1 ? bd[1] : 3;
    const zones = kids(base, 'Zone_t');
    for (const [, du] of kids(base, 'DimensionalUnits_t')) {
      // five 32-character fields: mass, length, time, temperature, angle
      const q = nums(du);
      let s = '';
      for (let i = 32; i < Math.min(64, q.length); i++) s += String.fromCharCode(q[i] & 0xff);
      const u = CGNS_UNITS[trimText(s).toLowerCase()];
      if (u) unitHint = u;
    }
    for (const [zoneName, zone] of zones) {
      const zt = (kids(zone, 'ZoneType_t')[0] || [])[1];
      const kind = zt ? text(zt) : 'Structured';
      const size = nums(zone);
      const gc = (kids(zone, 'GridCoordinates_t').find(([k]) => k === 'GridCoordinates') || kids(zone, 'GridCoordinates_t')[0] || [])[1];
      if (!gc) { warnings.push(`CGNS zone "${zoneName}" has no GridCoordinates and was left out.`); continue; }
      const arrays = kids(gc, 'DataArray_t');
      const coord = (axis, alt) => { const a = arrays.find(([k]) => k === `Coordinate${axis}`) || arrays.find(([k]) => k === alt); return a ? nums(a[1]) : null; };
      const X = coord('X', 'CoordinateR'), Y = coord('Y', 'CoordinateTheta'), Z = coord('Z', 'CoordinatePhi');
      if (arrays.some(([k]) => /Coordinate(R|Theta|Phi)$/.test(k)) && !arrays.some(([k]) => k === 'CoordinateX')) throw new Error(`${fname}: CGNS zone "${zoneName}" stores cylindrical or spherical coordinates, which are not read. Export the grid in Cartesian coordinates.`);
      if (!X) { warnings.push(`CGNS zone "${zoneName}" has no CoordinateX array and was left out.`); continue; }
      const nv = X.length;
      for (let i = 0; i < nv; i++) fe.addNode(nodeBase + i + 1, X[i], Y ? Y[i] : 0, Z && physDim > 2 ? Z[i] : 0);
      const sets = [];
      if (/^Unstructured/i.test(kind)) {
        if (size.length && size[0] !== nv) warnings.push(`CGNS zone "${zoneName}" declares ${size[0]} vertices but stores ${nv} coordinates; the stored ones were used.`);
        for (const [secName, sec] of kids(zone, 'Elements_t')) {
          const sd = nums(sec);
          const et = sd[0];
          const range = nums((kids(sec, 'IndexRange_t').find(([k]) => k === 'ElementRange') || [])[1]);
          const connNode = kids(sec, 'DataArray_t').find(([k]) => k === 'ElementConnectivity');
          const first = range.length ? range[0] : 1, count = range.length > 1 ? range[1] - range[0] + 1 : 0;
          const tname = CGNS_TYPES[et] ? CGNS_TYPES[et][0] : CGNS_NAMES[et] || `type ${et}`;
          sectionTypes[tname] = (sectionTypes[tname] || 0) + count;
          if (!connNode || et === 22 || et === 23 || et === 2 || (!CGNS_TYPES[et] && et !== 20)) {
            if (et !== 2) fe.skip(`CGNS ${et === 22 || et === 23 ? 'polyhedral ' : ''}section of ${tname} elements`, count || 1);
            continue;
          }
          const conn = nums(connNode[1]);
          const secTag = ++tag;
          fe.setProp(secTag, { name: secName });
          const ids = [];
          const emit = (def, at, k) => {
            const [type, npe, order] = def;
            const nodes = new Array(npe);
            for (let q = 0; q < npe; q++) nodes[q] = nodeBase + conn[at + (order ? order[q] : q)];
            const eid = elemBase + first + k;
            fe.addElement(type, nodes, eid, secTag);
            ids.push(eid);
            if (eid > top) top = eid;
          };
          if (et === 20) {
            // MIXED: each element is its type code followed by its nodes
            let p = 0, k = 0;
            while (p < conn.length) {
              const def = CGNS_TYPES[conn[p]];
              if (!def) { fe.skip(`CGNS element type ${CGNS_NAMES[conn[p]] || conn[p]} in MIXED section "${secName}" (the rest of the section was left out)`, Math.max(1, count - k)); break; }
              if (p + 1 + def[1] > conn.length) break;
              emit(def, p + 1, k);
              p += 1 + def[1]; k++;
            }
          } else {
            const def = CGNS_TYPES[et];
            const ne = Math.floor(conn.length / def[1]);
            for (let k = 0; k < ne; k++) emit(def, k * def[1], k);
          }
          fe.addElemSet(secName, ids);
          sets.push(secName);
        }
      } else {
        // vertex sizes come first in the zone size array; coordinates run with i fastest
        const idim = Math.max(1, Math.round(size.length / 3));
        const ni = size[0] || 1, nj = idim > 1 ? size[1] : 1, nk = idim > 2 ? size[2] : 1;
        if (ni * nj * nk !== nv) throw new Error(`${fname}: structured CGNS zone "${zoneName}" is ${ni} × ${nj} × ${nk} vertices but stores ${nv} coordinates (rind planes and partial coordinate arrays are not read).`);
        const secTag = ++tag;
        fe.setProp(secTag, { name: zoneName });
        const at = (i, j, k) => nodeBase + 1 + i + ni * (j + nj * k);
        const ids = [];
        let e = 0;
        if (idim >= 3 && nk > 1) {
          for (let k = 0; k < nk - 1; k++) for (let j = 0; j < nj - 1; j++) for (let i = 0; i < ni - 1; i++) {
            const eid = elemBase + (++e);
            fe.addElement('hex8', [at(i, j, k), at(i + 1, j, k), at(i + 1, j + 1, k), at(i, j + 1, k), at(i, j, k + 1), at(i + 1, j, k + 1), at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1)], eid, secTag);
            ids.push(eid);
          }
        } else if (idim >= 2 && nj > 1) {
          for (let j = 0; j < nj - 1; j++) for (let i = 0; i < ni - 1; i++) { const eid = elemBase + (++e); fe.addElement('quad4', [at(i, j, 0), at(i + 1, j, 0), at(i + 1, j + 1, 0), at(i, j + 1, 0)], eid, secTag); ids.push(eid); }
        } else for (let i = 0; i < ni - 1; i++) { const eid = elemBase + (++e); fe.addElement('line2', [at(i, 0, 0), at(i + 1, 0, 0)], eid, secTag); ids.push(eid); }
        fe.addElemSet(zoneName, ids);
        top = Math.max(top, elemBase + e);
        structured++;
      }
      // boundary conditions given on elements or vertices become sets (unstructured zones)
      for (const [, zbc] of kids(zone, 'ZoneBC_t')) {
        for (const [bcName, bc] of kids(zbc, 'BC_t')) {
          if (!/^Unstructured/i.test(kind)) { bcSkipped++; continue; }
          const find = (nm) => keysOf(bc).includes(nm) ? child(bc, nm) : null;
          const loc = find('GridLocation') ? text(find('GridLocation')) : 'Vertex';
          const onElems = !!find('ElementRange') || !!find('ElementList') || /Center/i.test(loc);
          const list = find('PointList') || find('ElementList'), rng = find('PointRange') || find('ElementRange');
          let members = null;
          if (list) members = Array.from(nums(list)); else if (rng) { const q = nums(rng); if (q.length >= 2) { members = []; for (let i = q[0]; i <= q[q.length - 1]; i++) members.push(i); } }
          if (!members || !members.length) { bcSkipped++; continue; }
          const fam = find('FamilyName') ? text(find('FamilyName')) : '';
          for (const nm of new Set([bcName, ...(fam ? [fam, ...(alias.get(fam) || [])] : [])])) {
            if (onElems) fe.addElemSet(nm, members.map((i) => elemBase + i)); else fe.addNodeSet(nm, members.map((i) => nodeBase + i));
          }
        }
      }
      zoneInfo.push({ base: baseName_, zone: zoneName, type: /^Unstructured/i.test(kind) ? 'Unstructured' : 'Structured', vertices: nv, sections: sets });
      zonesRead++;
      // node and element numbers restart in every zone: keep the next zone clear of this one
      nodeBase += nv;
      elemBase = top;
    }
  }
  if (zonesRead > 1) warnings.push(`${zonesRead} zones were read side by side; nodes on zone interfaces are not merged (weld them with the repair tool if one connected mesh is needed).`);
  if (bcSkipped) warnings.push(`${bcSkipped} boundary condition(s) on structured index ranges or without a point list were not turned into sets.`);
  return feResult(fe, 'cgns', opts, { name: fname, unitHint, expected: 'expected zones with GridCoordinates', warnings, extra: { flavour: 'cgns', cgnsVersion: version, bases: bases.map(([k]) => k), zones: zoneInfo, structuredZones: structured, sectionTypes, tree: ctx.tree } });
}

// ---------------------------------------------------------------------------------------------------------
// Fluent common fluids format
// ---------------------------------------------------------------------------------------------------------

function fluentCff(ctx) {
  const { file, opts, fname, warnings } = ctx;
  const meshes = child(file, 'meshes');
  const names = keysOf(meshes).filter((k) => isGroup(child(meshes, k)));
  const mesh = child(meshes, names[0]);
  if (names.length > 1) warnings.push(`The file holds ${names.length} meshes; the first ("${names[0]}") was read.`);
  const at = (g, ...path) => { let o = g; for (const k of path) o = isGroup(o) ? child(o, k) : null; return o; };
  const sets = (g) => keysOf(g).map((k) => child(g, k)).filter(isDataset);
  const first = (o) => { const v = attrsOf(o).minId; return Number.isFinite(v) && v >= 1 ? v : 1; };
  const topo = (g) => {
    const col = (k) => { const d = at(g, k); return isDataset(d) ? numbers(d) : []; };
    const nm = at(g, 'name');
    return { id: col('id'), type: col('zoneType'), min: col('minId'), max: col('maxId'), shape: col('cellType'), names: isDataset(nm) ? textOf(nm).split(';') : [] };
  };
  const ma = attrsOf(mesh);

  const coordSets = sets(at(mesh, 'nodes', 'coords'));
  if (!coordSets.length) throw new Error(`${fname}: the Fluent mesh "${names[0]}" has no node coordinates (meshes/${names[0]}/nodes/coords).`);
  let nNodes = 0, nodesRead = 0, coordDim = 0;
  const parts = coordSets.map((d) => { const shape = d.shape || []; const n = shape[0] || 0, nd = shape[1] || 1, f = first(d); nNodes = Math.max(nNodes, f + n - 1); return { d, n, nd, f }; });
  const coords = new Float64Array(nNodes * 3).fill(NaN);
  for (const { d, n, nd, f } of parts) {
    if (nd !== 2 && nd !== 3) throw new Error(`${fname}: the node coordinates of the Fluent mesh have ${nd} component(s) per node; expected 2 or 3.`);
    const v = numbers(d);
    for (let i = 0, o = (f - 1) * 3; i < n; i++, o += 3) { coords[o] = v[i * nd]; coords[o + 1] = v[i * nd + 1]; coords[o + 2] = nd === 3 ? v[i * nd + 2] : 0; }
    nodesRead += n; coordDim = Math.max(coordDim, nd);
  }

  // the two adjacent cells by face number
  const faces = at(mesh, 'faces');
  const sections = keysOf(at(faces, 'nodes')).map((k) => at(faces, 'nodes', k)).filter((g) => isDataset(at(g, 'nnodes')) && isDataset(at(g, 'nodes')));
  if (!sections.length) throw new Error(`${fname}: the Fluent mesh "${names[0]}" has no faces (meshes/${names[0]}/faces/nodes).`);
  let nFaces = 0;
  for (const g of sections) nFaces = Math.max(nFaces, first(g) + (at(g, 'nnodes').shape || [0])[0] - 1);
  const cell = [new Float64Array(nFaces + 1), new Float64Array(nFaces + 1)];
  ['c0', 'c1'].forEach((k, side) => { for (const d of sets(at(faces, k))) { const v = numbers(d), f = first(d); for (let i = 0; i < v.length && f + i <= nFaces; i++) cell[side][f + i] = v[i]; } });

  const ft = topo(at(faces, 'zoneTopology'));
  const zoneNames = new Map();
  const order = Array.from(ft.id, (_, i) => i).sort((a, b) => ft.min[a] - ft.min[b]);
  const faceZones = order.map((i) => { zoneNames.set(ft.id[i], { name: (ft.names[i] || '').trim() || `zone ${ft.id[i]}` }); return faceZone(ft.id[i], ft.type[i] ?? 0); });
  const stray = faceZone(0, 0);
  const ids = new Float64Array(4096);
  let zi = 0, bad = 0;
  for (const g of sections) {
    const nn = numbers(at(g, 'nnodes')), nodes = numbers(at(g, 'nodes'));
    let o = 0;
    for (let i = 0, f = first(g); i < nn.length; i++, f++) {
      const k = nn[i];
      if (!(k >= 0) || k > ids.length || o + k > nodes.length) throw new Error(`${fname}: the face table of the Fluent mesh is inconsistent at face ${f} (${k} nodes, ${nodes.length - o} node numbers left). The file is damaged.`);
      for (let j = 0; j < k; j++) ids[j] = nodes[o + j];
      o += k;
      // zones cover ascending ranges of face numbers
      if (zi >= order.length || f < ft.min[order[zi]] || f > ft.max[order[zi]]) { zi = order.findIndex((z) => f >= ft.min[z] && f <= ft.max[z]); if (zi < 0) { zi = order.length; bad++; } }
      addFace(zi < order.length ? faceZones[zi] : stray, ids, k, cell[0][f], cell[1][f], !!opts.twoSided);
    }
  }
  if (bad) { faceZones.push(stray); warnings.push(`${bad} face(s) lie outside every face zone of the file; those on the boundary are returned as "zone 0".`); }
  const ct = topo(at(mesh, 'cells', 'zoneTopology'));
  const cellZones = Array.from(ct.id, (zid, i) => { zoneNames.set(zid, { name: (ct.names[i] || '').trim() || `zone ${zid}` }); return { id: zid, count: Math.max(0, ct.max[i] - ct.min[i] + 1), typeCode: ct.type[i] ?? 1, shapeCode: ct.shape[i] ?? 0 }; });
  const text = (k) => { const d = at(file, 'settings', k); return isDataset(d) ? textOf(d) : null; };
  return fluentBoundary({ coords, nNodes, nodesRead, faceZones, zoneNames, cellZones, dim: ma.dimension === 2 || ma.dimension === 3 ? ma.dimension : coordDim, declNodes: ma.nodeCount || 0, declFaces: ma.faceCount || 0, nCells: ma.cellCount || cellZones.reduce((s, c) => s + c.count, 0) || null, notes: warnings }, { name: fname.replace(/\.h5$/i, '') }, opts, {
    flavour: 'fluent', encoding: 'hdf5', writer: text('Origin'), version: text('Version'), mesh: names[0], meshes: names, tree: ctx.tree,
  });
}

// ---------------------------------------------------------------------------------------------------------
// Generic HDF5: datasets named by the caller
// ---------------------------------------------------------------------------------------------------------

function generic(ctx) {
  const { file, opts, fname, name, warnings } = ctx;
  const spec = opts.datasets;
  const D = (path, what) => {
    const o = child(file, String(path).replace(/^\/+/, ''));
    if (!isDataset(o)) throw new Error(`${fname}: opts.datasets.${what} = "${path}" is not a dataset of the file. Datasets: ${ctx.tree.entries.filter((e) => e.type === 'dataset').slice(0, 20).map((e) => e.path).join(', ') || 'none'}.`);
    return o;
  };
  const extra = { flavour: 'generic', datasets: spec, tree: ctx.tree };
  if (spec.grid) {
    const ds = D(spec.grid, 'grid');
    const shape = ds.shape || [];
    if (shape.length !== 2) throw new Error(`${fname}: the grid dataset "${spec.grid}" is ${shape.length}-dimensional (${shape.join(' × ')}); a [ny, nx] dataset is needed.`);
    const [ny, nx] = shape;
    const raw = numbers(ds);
    let x0 = Number(spec.x0) || 0, y0 = Number(spec.y0) || 0, dx = Number(spec.dx) || 1, dy = Number(spec.dy) || 1;
    if (spec.x && spec.y) {
      const xs = numbers(D(spec.x, 'x')), ys = numbers(D(spec.y, 'y'));
      const rx = regular(xs), ry = regular(ys);
      if (xs.length !== nx || ys.length !== ny || rx === null || ry === null) throw new Error(`${fname}: the axis datasets must be evenly spaced and have ${nx} (x) and ${ny} (y) values; found ${xs.length} and ${ys.length}.`);
      x0 = Math.min(xs[0], xs[nx - 1]); y0 = Math.min(ys[0], ys[ny - 1]); dx = rx; dy = ry;
    } else if (!(spec.dx && spec.dy)) warnings.push('No axis datasets (opts.datasets.x, .y) or spacing (.dx, .dy) were given: the grid has unit spacing from the origin.');
    const z = new Float32Array(nx * ny);
    const nodata = spec.nodata;
    for (let iy = 0; iy < ny; iy++) { const oy = dy > 0 ? iy : ny - 1 - iy; for (let ix = 0; ix < nx; ix++) { const q = raw[iy * nx + ix]; z[oy * nx + (dx > 0 ? ix : nx - 1 - ix)] = q === nodata ? NaN : q; } }
    return finishGrid({ name, nx, ny, x0, y0, dx: Math.abs(dx), dy: Math.abs(dy), z }, { ...opts, maxCells: opts.maxCells ?? DEFAULT_MAX_CELLS }, { format: 'hdf5', crs: null, zName: 'elevation', warnings, extra });
  }
  if (spec.positions && spec.elements) {
    const pos = D(spec.positions, 'positions');
    const cols = (pos.shape || [])[1];
    if ((pos.shape || []).length !== 2 || !(cols === 2 || cols === 3)) throw new Error(`${fname}: the positions dataset "${spec.positions}" must be [n, 3] or [n, 2]; it is ${(pos.shape || []).join(' × ')}.`);
    const p = numbers(pos);
    const n = pos.shape[0];
    const fe = new FeMesh(name);
    for (let i = 0; i < n; i++) fe.addNode(i, p[cols * i], p[cols * i + 1], cols > 2 ? p[cols * i + 2] : 0);
    const npe = { line2: 2, tri3: 3, quad4: 4, tet4: 4, hex8: 8, wedge6: 6, pyr5: 5 };
    const blocks = [];
    let lo = Infinity, hi = -Infinity;
    for (const type in spec.elements) {
      if (!npe[type]) throw new Error(`${fname}: opts.datasets.elements has the unknown element type "${type}" (use ${Object.keys(npe).join(', ')}).`);
      const ds = D(spec.elements[type], `elements.${type}`);
      const c = numbers(ds);
      if (c.length % npe[type]) throw new Error(`${fname}: the ${type} connectivity "${spec.elements[type]}" has ${c.length} entries, not a multiple of ${npe[type]}.`);
      for (let i = 0; i < c.length; i++) { if (c[i] < lo) lo = c[i]; if (c[i] > hi) hi = c[i]; }
      blocks.push([type, c]);
    }
    // 1-based when the caller says so, or when the indices can only be 1-based
    const base = Number.isInteger(spec.indexBase) ? spec.indexBase : (lo >= 1 && hi === n ? 1 : 0);
    let eid = 0;
    for (const [type, c] of blocks) { const k = npe[type]; const nodes = new Array(k); for (let e = 0; e < c.length / k; e++) { for (let q = 0; q < k; q++) nodes[q] = c[e * k + q] - base; fe.addElement(type, nodes, ++eid, 0); } }
    return feResult(fe, 'hdf5', opts, { name: fname, expected: 'expected rows in the positions dataset', warnings, extra: { ...extra, indexBase: base } });
  }
  if (spec.points || (spec.x && spec.y && spec.z)) {
    let coords;
    if (spec.points) {
      const ds = D(spec.points, 'points');
      const cols = (ds.shape || [])[1];
      if ((ds.shape || []).length !== 2 || !(cols >= 2)) throw new Error(`${fname}: the points dataset "${spec.points}" must be [n, 3] (or [n, 2]); it is ${(ds.shape || []).join(' × ')}.`);
      const p = numbers(ds);
      const n = ds.shape[0];
      coords = new Float64Array(3 * n);
      for (let i = 0; i < n; i++) { coords[3 * i] = p[cols * i]; coords[3 * i + 1] = p[cols * i + 1]; coords[3 * i + 2] = cols > 2 ? p[cols * i + 2] : 0; }
    } else {
      const x = numbers(D(spec.x, 'x')), y = numbers(D(spec.y, 'y')), z = numbers(D(spec.z, 'z'));
      if (x.length !== y.length || x.length !== z.length) throw new Error(`${fname}: the x, y and z datasets have ${x.length}, ${y.length} and ${z.length} values; they must be equally long.`);
      coords = new Float64Array(3 * x.length);
      for (let i = 0; i < x.length; i++) { coords[3 * i] = x[i]; coords[3 * i + 1] = y[i]; coords[3 * i + 2] = z[i]; }
    }
    if (!coords.length) throw new Error(`${fname}: the points dataset is empty.`);
    const F = new Features(name);
    for (let i = 0; i < coords.length; i += 3) F.addPoint(coords[i], coords[i + 1], coords[i + 2], null);
    return F.build(opts, { format: 'hdf5', fileName: fname, crs: null, warnings, expected: 'expected rows of x, y, z', extra });
  }
  throw new Error(`${fname}: opts.datasets must name what to read: { grid, x?, y? } for a [ny, nx] raster, { points } or { x, y, z } for a point set, or { positions, elements: { tri3 | quad4 | tet4 | hex8 | wedge6 | pyr5 | line2: path } } for a mesh.`);
}

// ---------------------------------------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------------------------------------

/**
 * Decide the convention of an open file from its object tree.
 * @returns {'med'|'cgns'|'exodus'|'fluent'|null} null: NetCDF-style variables decide (UGRID, CF grid or generic)
 */
function meshFlavour(file, rootAtts) {
  const keys = keysOf(file);
  if (keys.includes('ENS_MAA') && isGroup(child(file, 'ENS_MAA'))) return 'med';
  if (rootAtts.label === 'Root Node of HDF5 File' || keys.includes('CGNSLibraryVersion')) return 'cgns';
  for (const k of keys.slice(0, 32)) { const o = child(file, k); if (isGroup(o) && attrsOf(o).label === 'CGNSBase_t') return 'cgns'; }
  if ((keys.includes('coordx') || keys.includes('coord')) && keys.some((k) => /^connect\d+$/.test(k))) return 'exodus';
  if (keys.includes('meshes')) { const ms = child(file, 'meshes'); if (isGroup(ms) && keysOf(ms).some((k) => { const m = child(ms, k); return isGroup(m) && isGroup(child(m, 'faces')) && isGroup(child(m, 'nodes')); })) return 'fluent'; }
  return null;
}

/**
 * Import an HDF5 file: NetCDF-4 CF grid or UGRID mesh, MED, CGNS, Exodus II, Fluent mesh / case, or datasets
 * named by the caller.
 * @param {{ name?: string, bytes?: Uint8Array }} input file content
 * @param {object} [opts] options of the grid, vector and mesh importers (`units`, `scale`, `crs`, `vertical`,
 *   `zUnits`, `projection`, `origin`) and
 * @param {number} [opts.maxCells] grids: block-average to at most this many nodes (default 1 500 000; 0 keeps all)
 * @param {string} [opts.variable] grids: name of the variable to read
 * @param {string} [opts.mesh] MED: name of the mesh to read (default: the first)
 * @param {object} [opts.datasets] generic HDF5: `{ grid, x?, y?, x0?, y0?, dx?, dy?, nodata? }`,
 *   `{ points }`, `{ x, y, z }` or `{ positions, elements: { tri3: path, … }, indexBase? }` (dataset paths)
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} meta.flavour is one of
 *   'netcdf4-cf-grid', 'netcdf4-ugrid', 'exodus-netcdf4', 'med', 'cgns', 'fluent', 'generic'; meta.tree is
 *   `{ entries: [{ path, type, shape?, dtype? }], truncated }`
 */
export async function parse(input, opts = {}) {
  const fname = (input && input.name) || 'HDF5';
  const bytes = inputBytes(input || {});
  if (!sniff({ bytes })) {
    const classic = bytes.length > 3 && bytes[0] === 0x43 && bytes[1] === 0x44 && bytes[2] === 0x46;
    const adf = bytes.length > 32 && /^@\(#\)ADF Database Version/.test(String.fromCharCode.apply(null, bytes.subarray(0, 28)));
    throw new Error(`${fname}: not an HDF5 file (the HDF5 signature was not found).${classic ? ' It is a NetCDF classic file, which is read directly (import it as NetCDF).' : adf ? ' It is a CGNS file in the older ADF format: convert it with "cgnsconvert -h in.cgns out.cgns".' : ''}`);
  }
  const { h5, FS } = await load();
  const path = `/import-${++serial}.h5`;
  let file = null;
  try {
    try {
      // the in-memory file system takes the caller's array as the file content instead of copying it; the
      // file is opened read-only, so the array is never written to
      FS.writeFile(path, bytes, { canOwn: true });
      file = new h5.File(path, 'r');
    } catch (err) {
      file = null;
      throw new Error(`${fname}: the HDF5 file could not be opened (${reason(err)}). The file is probably truncated or damaged; "h5dump -H" or "ncdump -h" shows whether it can be read at all.`);
    }
    const rootAtts = attrsOf(file);
    const tree = listTree(file);
    const warnings = [];
    if (tree.truncated) warnings.push(`The file has more than ${TREE_LIMIT} objects (or very deep groups); only the first ${tree.entries.length} are listed and searched.`);
    const ctx = { h5, file, opts, fname, name: baseName(fname), warnings, tree, rootAtts };
    const run = (what, fn) => {
      try { return fn(); } catch (err) {
        // errors raised by the library carry its error stack; anything else is already a message for the user
        if (!/HDF5-DIAG|H5\w*\.c line \d+/.test(String(err && err.message))) throw err;
        throw new Error(`${fname}: a dataset of the ${what} could not be read (${reason(err)}). The file is probably truncated or damaged.`);
      }
    };
    if (opts.datasets) return run('file', () => generic(ctx));
    const flavour = meshFlavour(file, rootAtts);
    if (flavour === 'med') return run('MED mesh', () => med(ctx));
    if (flavour === 'cgns') return run('CGNS file', () => cgns(ctx));
    if (flavour === 'exodus') return run('Exodus mesh', () => exodus(ctx));
    if (flavour === 'fluent') return run('Fluent mesh', () => fluentCff(ctx));
    return run('file', () => {
      const vars = collectVariables(file, tree);
      const topo = vars.find((v) => v.atts.cf_role === 'mesh_topology' && v.atts.face_node_connectivity);
      if (topo && !opts.variable) return ugrid(ctx, vars, topo);
      const cands = gridCandidates(vars);
      if (cands.length) return cfGrid(ctx, vars, cands);
      if (opts.variable) throw new Error(`${fname}: variable "${opts.variable}" is not a gridded variable on x/y axes; the file has none.`);
      const sets = tree.entries.filter((e) => e.type === 'dataset');
      warnings.push(`No geometry convention was recognised in this HDF5 file (not a CF grid, UGRID mesh, MED, CGNS or Exodus file), so nothing was imported. It holds ${sets.length} dataset(s)${sets.length ? `: ${sets.slice(0, 12).map((e) => `${e.path} [${e.shape.join(' × ')}]`).join(', ')}${sets.length > 12 ? ', …' : ''}` : ''}. Name the datasets to read with opts.datasets ({ grid, x, y }, { points } or { positions, elements }); the full listing is in meta.tree.`);
      return finish([], warnings, { format: 'hdf5', flavour: 'generic', units: null, crs: null, tree });
    });
  } finally {
    try { if (file) file.close(); } catch { /* the handle is gone either way */ }
    try { FS.unlink(path); } catch { /* nothing was written */ }
  }
}
