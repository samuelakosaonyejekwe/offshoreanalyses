// NetCDF classic files (CDF-1, CDF-2 64-bit offset, CDF-5 64-bit data) read directly from the bytes:
//   - CF-convention gridded elevation / bathymetry (and legacy GMT grids) → grid
//   - UGRID 2-D mesh topologies (face_node_connectivity) → surface
//   - Exodus II meshes written on the classic format → frame / surface / volume
// NetCDF-4 and Exodus files stored as HDF5 are not read here (see the hdf5 adapter in the registry).
//
// Grid reading: the variable is chosen by `opts.variable`, else by name (elevation, z, Band1, bathymetry,
// topo, depth, altitude, height), else the first 2-D numeric variable on recognised x/y axes. Leading
// dimensions (time, level) are taken at index 0. scale_factor, add_offset, _FillValue, missing_value and
// valid_range are applied. `positive = "down"` or a depth standard name flips the sign so that z is up.

import { inputBytes, baseName, FeMesh, decodeText, triangulate3D } from './common.js';
import { parsePrj, epsgInfo, normaliseUnit } from '../units.js';
import { finishGrid } from './grids.js';
import { feResult } from './fe.js';
import { Features } from './vector.js';

const NC_SIZE = [0, 1, 1, 2, 4, 4, 8, 1, 2, 4, 8, 8];
const NC_GET = [null, 'getInt8', 'getUint8', 'getInt16', 'getInt32', 'getFloat32', 'getFloat64', 'getUint8', 'getUint16', 'getUint32', 'getBigInt64', 'getBigUint64'];

/**
 * Parse the header of a classic NetCDF file.
 * @returns {{ version: number, numrecs: number, dims: {name, len}[], atts: object, vars: object[], read: Function, readText: Function }}
 */
export function openNetcdf(bytes, fileName = 'NetCDF') {
  if (bytes.length < 32 || bytes[0] !== 0x43 || bytes[1] !== 0x44 || bytes[2] !== 0x46) {
    if (bytes.length >= 8 && bytes[0] === 0x89 && bytes[1] === 0x48 && bytes[2] === 0x44 && bytes[3] === 0x46) throw new Error(`${fileName}: this is a NetCDF-4 / HDF5 file, which the classic NetCDF reader cannot open. Convert it with "nccopy -k classic in.nc out.nc" (or "ncks -3").`);
    throw new Error(`${fileName}: not a NetCDF classic file (expected the magic bytes "CDF" followed by 1, 2 or 5).`);
  }
  const version = bytes[3];
  if (![1, 2, 5].includes(version)) throw new Error(`${fileName}: NetCDF classic format version ${version} is not known (expected 1, 2 or 5).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const wide = version === 5;
  let p = 4;
  const u32 = () => { const v = dv.getUint32(p, false); p += 4; return v; };
  const nn = () => { if (!wide) return u32(); const v = Number(dv.getBigUint64(p, false)); p += 8; return v; };
  const str = () => { const n = nn(); const s = decodeText(bytes.subarray(p, p + n)); p += n + ((4 - (n % 4)) % 4); return s; };
  const values = (type, n) => {
    const size = NC_SIZE[type];
    if (!size) throw new Error(`${fileName}: unknown NetCDF data type ${type} in the header.`);
    let out;
    if (type === 2) out = decodeText(bytes.subarray(p, p + n)).replace(/\u0000+$/, '');
    else { out = new Array(n); for (let i = 0; i < n; i++) { const v = dv[NC_GET[type]](p + i * size, false); out[i] = typeof v === 'bigint' ? Number(v) : v; } }
    p += n * size + ((4 - ((n * size) % 4)) % 4);
    return out;
  };
  const attList = () => {
    const tag = u32(), n = nn();
    const o = {};
    if (tag === 0 && n === 0) return o;
    if (tag !== 0x0c) throw new Error(`${fileName}: corrupt NetCDF header (expected an attribute list).`);
    for (let i = 0; i < n; i++) { const nm = str(); const type = u32(); const cnt = nn(); const v = values(type, cnt); o[nm] = typeof v === 'string' ? v : v.length === 1 ? v[0] : v; }
    return o;
  };
  let numrecs, atts;
  const dims = [], vars = [];
  try {
  numrecs = nn();
  { const tag = u32(), n = nn(); if (tag !== 0x0a && !(tag === 0 && n === 0)) throw new Error(`${fileName}: corrupt NetCDF header (expected a dimension list).`); for (let i = 0; i < n; i++) dims.push({ name: str(), len: nn() }); }
  atts = attList();
  {
    const tag = u32(), n = nn();
    if (tag !== 0x0b && !(tag === 0 && n === 0)) throw new Error(`${fileName}: corrupt NetCDF header (expected a variable list).`);
    for (let i = 0; i < n; i++) {
      const name = str();
      const nd = nn();
      const dimids = [];
      for (let k = 0; k < nd; k++) dimids.push(nn());
      const a = attList();
      const type = u32();
      const vsize = nn();
      let begin;
      if (version === 1) begin = u32(); else { begin = Number(dv.getBigUint64(p, false)); p += 8; }
      const record = nd > 0 && dims[dimids[0]] && dims[dimids[0]].len === 0;
      vars.push({ name, dimids, dims: dimids.map((d) => (dims[d] ? dims[d].name : '?')), shape: dimids.map((d) => (dims[d] ? dims[d].len : 0)), atts: a, type, vsize, begin, record });
    }
  }
  } catch (e) {
    // reading past the end of the buffer or an absurd count: the header is cut short or damaged
    if (e instanceof RangeError) throw new Error(`${fileName}: the NetCDF header is truncated or corrupt (it ends after ${bytes.length} bytes while the dimension, attribute and variable lists are still being read).`);
    throw e;
  }
  const recVars = vars.filter((v) => v.record);
  let recSize = recVars.reduce((s, v) => s + v.vsize, 0);
  if (recVars.length === 1) { const v = recVars[0]; recSize = NC_SIZE[v.type] * v.shape.slice(1).reduce((a, b) => a * b, 1); }
  if (numrecs === 0xffffffff && recSize) numrecs = Math.floor((bytes.length - Math.min(...recVars.map((v) => v.begin))) / recSize);
  for (const v of vars) if (v.record) v.shape[0] = numrecs;
  /**
   * Read `count` values of variable `v` starting at flat index `start` (row-major over its shape).
   * A record variable is read across records transparently.
   */
  const read = (v, start = 0, count) => {
    const size = NC_SIZE[v.type], get = NC_GET[v.type];
    const total = v.shape.reduce((a, b) => a * b, 1);
    const n = count === undefined ? total - start : count;
    if (!size) throw new Error(`${fileName}: variable "${v.name}" has an unknown data type (${v.type}).`);
    if (!(n >= 0) || n * size > bytes.length) throw new Error(`${fileName}: variable "${v.name}" declares ${n} values, more than the file holds (truncated NetCDF).`);
    const out = new Float64Array(n);
    const perRec = v.record ? v.shape.slice(1).reduce((a, b) => a * b, 1) : total;
    for (let i = 0; i < n; i++) {
      const k = start + i;
      const o = v.record ? v.begin + Math.floor(k / perRec) * recSize + (k % perRec) * size : v.begin + k * size;
      if (o + size > bytes.length) throw new Error(`${fileName}: variable "${v.name}" extends past the end of the file (truncated NetCDF).`);
      const x = dv[get](o, false);
      out[i] = typeof x === 'bigint' ? Number(x) : x;
    }
    return out;
  };
  /** Rows of a 2-D char variable as trimmed strings. */
  const readText = (v) => {
    if (v.type !== 2 || v.shape.length !== 2) return [];
    const [n, w] = v.shape;
    const out = [];
    for (let i = 0; i < n; i++) out.push(decodeText(bytes.subarray(v.begin + i * w, v.begin + (i + 1) * w)).replace(/\u0000.*$/, '').trim());
    return out;
  };
  return { version, numrecs, dims, atts, vars, read, readText };
}

const isX = (v) => v && (v.atts.axis === 'X' || /^(longitude|projection_x_coordinate|grid_longitude)$/.test(v.atts.standard_name || '') || /degrees?_?e/i.test(v.atts.units || '') || /^(x|lon|long|longitude|easting|xc|x_range|eastings?)$/i.test(v.name));
const isY = (v) => v && (v.atts.axis === 'Y' || /^(latitude|projection_y_coordinate|grid_latitude)$/.test(v.atts.standard_name || '') || /degrees?_?n/i.test(v.atts.units || '') || /^(y|lat|latitude|northing|yc|y_range|northings?)$/i.test(v.name));
const isDegrees = (v) => v && (/degree/i.test(v.atts.units || '') || /^(longitude|latitude)$/.test(v.atts.standard_name || '') || /^(lon|lat|longitude|latitude)$/i.test(v.name));
const Z_NAMES = /^(elevation|z|band1|bathymetry|bathy|topo|topography|depth|altitude|height|deptho|elev|dem|h)$/i;

function gridCrs(nc, v, x, y) {
  const gm = v.atts.grid_mapping ? nc.vars.find((q) => q.name === v.atts.grid_mapping) : nc.vars.find((q) => q.atts.grid_mapping_name);
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

function exodus(nc, input, opts) {
  const name = baseName(input.name);
  const fe = new FeMesh(name);
  const V = (nm) => nc.vars.find((v) => v.name === nm);
  const nNodes = (nc.dims.find((d) => d.name === 'num_nodes') || { len: 0 }).len;
  const cx = V('coordx'), cy = V('coordy'), cz = V('coordz'), c = V('coord');
  let X, Y, Z;
  if (cx) { X = nc.read(cx); Y = cy ? nc.read(cy) : new Float64Array(nNodes); Z = cz ? nc.read(cz) : new Float64Array(nNodes); } else if (c) { const all = nc.read(c); const nd = c.shape[0]; X = all.subarray(0, nNodes); Y = nd > 1 ? all.subarray(nNodes, 2 * nNodes) : new Float64Array(nNodes); Z = nd > 2 ? all.subarray(2 * nNodes, 3 * nNodes) : new Float64Array(nNodes); } else throw new Error(`${input.name || 'Exodus'}: the Exodus file has no coordx / coord variable.`);
  for (let i = 0; i < nNodes; i++) fe.addNode(i + 1, X[i], Y[i], Z[i]);
  const ids = V('eb_prop1') ? nc.read(V('eb_prop1')) : [];
  const names = V('eb_names') ? nc.readText(V('eb_names')) : [];
  const HEX20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15];
  let eid = 0, blk = 0;
  for (const v of nc.vars) {
    const m = /^connect(\d+)$/.exec(v.name);
    if (!m || v.shape.length !== 2) continue;
    const k = Number(m[1]);
    const [ne, npe] = v.shape;
    const et = String(v.atts.elem_type || '').toUpperCase();
    let type = null, order = null;
    if (/^(BAR|BEAM|TRUSS|EDGE|ROD)/.test(et)) type = npe === 3 ? 'line3' : npe === 2 ? 'line2' : null;
    else if (/^TRI/.test(et)) type = npe === 6 ? 'tri6' : npe === 3 ? 'tri3' : null;
    else if (/^(QUAD|SHELL)/.test(et)) type = npe === 9 ? 'quad9' : npe === 8 ? 'quad8' : npe === 4 ? 'quad4' : null;
    else if (/^TET/.test(et)) type = npe === 10 ? 'tet10' : npe === 4 ? 'tet4' : null;
    else if (/^HEX/.test(et)) { type = npe === 27 ? 'hex27' : npe === 20 ? 'hex20' : npe === 8 ? 'hex8' : null; if (npe === 20) order = HEX20; if (npe === 27) order = HEX20.concat([20, 21, 22, 23, 24, 25, 26]); } else if (/^WEDGE/.test(et)) type = npe === 15 ? 'wedge15' : npe === 6 ? 'wedge6' : null;
    else if (/^PYR/.test(et)) type = npe === 13 ? 'pyr13' : npe === 5 ? 'pyr5' : null;
    const tag = ids[k - 1] || k;
    fe.setProp(tag, { name: names[k - 1] || `block ${tag}` });
    blk++;
    if (!type) { fe.skip(`element type ${et || '?'} with ${npe} nodes`, ne); eid += ne; continue; }
    const conn = nc.read(v);
    const nodes = new Array(npe);
    const set = [];
    for (let e = 0; e < ne; e++) { for (let q = 0; q < npe; q++) nodes[q] = conn[e * npe + (order ? order[q] : q)]; fe.addElement(type, nodes, ++eid, tag); set.push(eid); }
    fe.addElemSet(names[k - 1] || `block ${tag}`, set);
  }
  const nsNames = V('ns_names') ? nc.readText(V('ns_names')) : [], nsIds = V('ns_prop1') ? nc.read(V('ns_prop1')) : [];
  for (const v of nc.vars) { const m = /^node_ns(\d+)$/.exec(v.name); if (m) fe.addNodeSet(nsNames[m[1] - 1] || `nodeset ${nsIds[m[1] - 1] || m[1]}`, nc.read(v)); }
  const sideSets = nc.vars.filter((v) => /^elem_ss\d+$/.test(v.name)).length;
  const warnings = [];
  if (sideSets) warnings.push(`${sideSets} side set(s) were not read (node sets and element blocks are).`);
  return feResult(fe, 'exodus', opts, { name: input.name || 'Exodus', expected: 'expected num_nodes > 0', warnings, skipped: sideSets ? { sideSets } : {}, extra: { title: nc.atts.title || null, blocks: blk, netcdfVersion: nc.version } });
}

function ugrid(nc, topo, input, opts) {
  const V = (nm) => nc.vars.find((v) => v.name === nm);
  const [xn, yn] = String(topo.atts.node_coordinates || '').trim().split(/\s+/);
  const xv = V(xn), yv = V(yn), fv = V(topo.atts.face_node_connectivity);
  if (!xv || !yv || !fv) throw new Error(`${input.name || 'NetCDF'}: the UGRID mesh topology "${topo.name}" does not name node coordinates and a face_node_connectivity variable that exist in the file.`);
  const x = nc.read(xv), y = nc.read(yv);
  const n = x.length;
  const nodeDim = xv.dims[0];
  const zv = nc.vars.find((v) => v.dims.length === 1 && v.dims[0] === nodeDim && v !== xv && v !== yv && (/altitude|depth|elevation|height/.test(v.atts.standard_name || '') || /node_z|^z$|bath|depth|elev/i.test(v.name)));
  const z = zv ? nc.read(zv) : new Float64Array(n);
  const down = zv && (zv.atts.positive === 'down' || /depth/.test(zv.atts.standard_name || ''));
  const fill = zv ? zv.atts._FillValue : undefined;
  const coords = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) { coords[3 * i] = x[i]; coords[3 * i + 1] = y[i]; coords[3 * i + 2] = z[i] === fill || z[i] !== z[i] ? 0 : (down ? -z[i] : z[i]); }
  const conn = nc.read(fv);
  const [nf, mx] = fv.shape;
  const start = fv.atts.start_index ?? 0, ffill = fv.atts._FillValue;
  const tris = [];
  const hist = {};
  for (let f = 0; f < nf; f++) {
    const ids = [];
    for (let k = 0; k < mx; k++) { const v = conn[f * mx + k]; if (v === ffill || v < start || v !== v) break; ids.push(v - start); }
    if (ids.length < 3 || ids.some((k) => k >= n)) continue;
    hist[ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon'] = (hist[ids.length === 3 ? 'tri3' : ids.length === 4 ? 'quad4' : 'polygon'] || 0) + 1;
    if (ids.length === 3) tris.push(ids[0], ids[1], ids[2]);
    else if (ids.length === 4) tris.push(ids[0], ids[1], ids[2], ids[0], ids[2], ids[3]);
    else { const flat = new Float64Array(ids.length * 3); ids.forEach((k, i) => { flat[3 * i] = coords[3 * k]; flat[3 * i + 1] = coords[3 * k + 1]; flat[3 * i + 2] = 0; }); const t = triangulate3D(flat); for (const k of t) tris.push(ids[k]); }
  }
  if (!tris.length) throw new Error(`${input.name || 'NetCDF'}: the UGRID mesh "${topo.name}" has no faces with three or more nodes.`);
  const F = new Features(baseName(input.name));
  F.addMesh(coords, Uint32Array.from(tris), topo.name);
  const crs = gridCrs(nc, zv || xv, xv, yv);
  const r = F.build(opts, { format: 'netcdf', fileName: input.name, crs: crs ? { ...crs, source: crs.epsg ? `EPSG:${crs.epsg}` : crs.name } : null, zScale: 1, expected: '', extra: { convention: 'UGRID', mesh: topo.name, nodeElevation: zv ? zv.name : null, netcdfVersion: nc.version } });
  r.meta.elementTypes = hist;
  if (!zv) r.warnings.push('The UGRID mesh has no node elevation variable; z is 0.');
  return r;
}

/** NetCDF classic: CF grid, GMT grid, UGRID mesh or Exodus II mesh. */
export function parseNetcdf(input, opts = {}) {
  const bytes = inputBytes(input);
  const nc = openNetcdf(bytes, input.name || 'NetCDF');
  const name = baseName(input.name);
  if (nc.dims.some((d) => d.name === 'num_nodes') && nc.vars.some((v) => v.name === 'coordx' || v.name === 'coord')) return exodus(nc, input, opts);
  const topo = nc.vars.find((v) => v.atts.cf_role === 'mesh_topology' && v.atts.face_node_connectivity);
  if (topo && !opts.variable) return ugrid(nc, topo, input, opts);
  const V = (nm) => nc.vars.find((v) => v.name === nm);
  const warnings = [];
  const mask = (data, v) => {
    const fill = v.atts._FillValue, miss = v.atts.missing_value;
    const vr = Array.isArray(v.atts.valid_range) ? v.atts.valid_range : [v.atts.valid_min ?? -Infinity, v.atts.valid_max ?? Infinity];
    const sf = v.atts.scale_factor ?? 1, ao = v.atts.add_offset ?? 0;
    const out = new Float32Array(data.length);
    for (let i = 0; i < data.length; i++) { const q = data[i]; out[i] = q === fill || q === miss || q < vr[0] || q > vr[1] || q !== q ? NaN : q * sf + ao; }
    return out;
  };

  // legacy GMT grid: x_range, y_range, spacing, dimension, z (rows north → south)
  if (V('x_range') && V('y_range') && V('dimension') && V('z')) {
    const xr = nc.read(V('x_range')), yr = nc.read(V('y_range')), dim = nc.read(V('dimension'));
    const nx = dim[0], ny = dim[1];
    const zv = V('z');
    const raw = mask(nc.read(zv), zv);
    const pixel = zv.atts.node_offset === 1;
    const dx = (xr[1] - xr[0]) / (pixel ? nx : nx - 1), dy = (yr[1] - yr[0]) / (pixel ? ny : ny - 1);
    const z = new Float32Array(nx * ny);
    for (let r = 0; r < ny; r++) z.set(raw.subarray(r * nx, (r + 1) * nx), (ny - 1 - r) * nx);
    const xu = V('x_range').atts.units || '';
    return finishGrid({ name, nx, ny, x0: xr[0] + (pixel ? dx / 2 : 0), y0: yr[0] + (pixel ? dy / 2 : 0), dx, dy, z }, opts, { format: 'netcdf', crs: /degree/i.test(xu) ? { kind: 'geographic', epsg: 4326, name: 'WGS 84 (assumed from degree units)', unit: 'deg' } : null, zName: zv.atts.long_name || 'z', warnings, extra: { convention: 'GMT (legacy)', variable: 'z', netcdfVersion: nc.version } });
  }

  const coordOf = (dimName) => nc.vars.find((v) => v.name === dimName && v.dims.length === 1);
  const cands = [];
  for (const v of nc.vars) {
    if (v.dims.length < 2 || v.type === 2) continue;
    const a = coordOf(v.dims[v.dims.length - 2]), b = coordOf(v.dims[v.dims.length - 1]);
    let yx = null;
    if (isY(a) && isX(b)) yx = true; else if (isX(a) && isY(b)) yx = false; else if (a && b && v.dims.length === 2 && !isX(a) && !isY(b)) yx = true;
    if (yx === null) continue;
    cands.push({ v, y: yx ? a : b, x: yx ? b : a, yx });
  }
  let pick = null;
  if (opts.variable) { pick = cands.find((c) => c.v.name === opts.variable); if (!pick) throw new Error(`${input.name || 'NetCDF'}: variable "${opts.variable}" is not a gridded variable on x/y axes. Gridded variables: ${cands.map((c) => c.v.name).join(', ') || 'none'}.`); } else pick = cands.find((c) => Z_NAMES.test(c.v.name)) || cands.find((c) => /altitude|depth|height|elevation/.test(c.v.atts.standard_name || '')) || cands[0];
  if (!pick) throw new Error(`${input.name || 'NetCDF'}: no gridded variable on recognisable x/y (or lon/lat) coordinate axes was found. Variables: ${nc.vars.map((v) => `${v.name}(${v.dims.join(', ')})`).join(', ') || 'none'}. This reader handles CF grids, legacy GMT grids, UGRID meshes and Exodus meshes.`);
  const { v, x, y, yx } = pick;
  if (cands.length > 1 && !opts.variable) warnings.push(`Several gridded variables are present (${cands.map((c) => c.v.name).join(', ')}); "${v.name}" was used. Choose another with opts.variable.`);
  const n1 = v.shape[v.shape.length - 2], n2 = v.shape[v.shape.length - 1];
  const lead = v.shape.slice(0, -2);
  if (lead.some((l) => l > 1)) warnings.push(`"${v.name}" has leading dimension(s) ${v.dims.slice(0, -2).map((d, i) => `${d} = ${lead[i]}`).join(', ')}; the first index of each was used.`);
  if (lead.some((l) => l === 0)) throw new Error(`${input.name || 'NetCDF'}: variable "${v.name}" has no records.`);
  const raw = mask(nc.read(v, 0, n1 * n2), v);
  const xs = nc.read(x), ys = nc.read(y);
  const nx = xs.length, ny = ys.length;
  const dx = regular(xs), dy = regular(ys);
  const at = (ix, iy) => (yx ? raw[iy * nx + ix] : raw[ix * ny + iy]);
  const crs = gridCrs(nc, v, x, y);
  const depthLike = v.atts.positive === 'down' || /depth/.test(v.atts.standard_name || '') || /^(depth|deptho)$/i.test(v.name);
  const extra = { convention: nc.atts.Conventions || null, variable: v.name, longName: v.atts.long_name || null, standardName: v.atts.standard_name || null, netcdfVersion: nc.version, axes: { x: x.name, y: y.name } };
  if (dx === null || dy === null) {
    // unevenly spaced axes cannot be a grid model: return the nodes as soundings
    const F = new Features(name);
    for (let iy = 0; iy < ny; iy++) for (let ix = 0; ix < nx; ix++) { const q = at(ix, iy); if (q === q) F.addPoint(xs[ix], ys[iy], depthLike && !opts.vertical ? -q : q, null); }
    warnings.push('The coordinate axes are not evenly spaced, so the grid nodes are returned as points.');
    return F.build(opts, { format: 'netcdf', fileName: input.name, crs: crs ? { ...crs, source: crs.epsg ? `EPSG:${crs.epsg}` : crs.name } : null, zScale: 1, warnings, expected: 'expected valid cells', extra });
  }
  const z = new Float32Array(nx * ny);
  for (let iy = 0; iy < ny; iy++) { const oy = dy > 0 ? iy : ny - 1 - iy; for (let ix = 0; ix < nx; ix++) z[oy * nx + (dx > 0 ? ix : nx - 1 - ix)] = at(ix, iy); }
  return finishGrid({ name, nx, ny, x0: dx > 0 ? xs[0] : xs[nx - 1], y0: dy > 0 ? ys[0] : ys[ny - 1], dx: Math.abs(dx), dy: Math.abs(dy), z }, opts, {
    format: 'netcdf', crs, zName: depthLike ? 'depth' : 'elevation', zUnit: normaliseUnit(v.atts.units) ? v.atts.units : (crs && crs.kind === 'geographic' ? 'm' : null), warnings, extra,
  });
}
