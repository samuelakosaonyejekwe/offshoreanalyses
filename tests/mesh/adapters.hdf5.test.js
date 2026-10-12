import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import * as hdf5 from '../../src/geometry/adapters/hdf5.js';
import { importFile } from '../../src/geometry/registry.js';
import { inspect, cellVolume } from '../../src/mesh/topology.js';
import { close, rel } from './helpers.js';

// ---------------------------------------------------------------------------------------------------------
// Fixtures are written in memory with the writer of the vendored HDF5 library itself, following the object
// layouts of the NetCDF-4 (dimension scales + CF attributes), MED 4, CGNS / HDF5 and Exodus II conventions.
// The MED and CGNS trees and the MED node orders below were copied from files exported by Gmsh 4.15 (MED 4.1
// and CGNS 3.4 libraries) for the same four cells.
// ---------------------------------------------------------------------------------------------------------

const h5 = await import(new URL('../../vendor/h5wasm/hdf5_hl.js', import.meta.url).href);
const lib = await h5.ready;
let serial = 0;

/** Build an HDF5 file with `fill(file)` and return its bytes. */
function write(fill) {
  const path = `/fixture-${++serial}.h5`;
  const f = new h5.File(path, 'w');
  try { fill(f); } finally { f.close(); }
  const bytes = lib.FS.readFile(path).slice();
  lib.FS.unlink(path);
  return bytes;
}

const chars = (s, n = s.length) => { const a = new Int8Array(n).fill(32); for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i); return a; };
const i64 = (a) => BigInt64Array.from(a, (v) => BigInt(v));
const f64 = (a) => Float64Array.from(a);
// 64-bit integers (the integer type of the MED library) have to be declared; the rest is inferred
const setAttrs = (obj, attrs) => { for (const k in attrs) { if (typeof attrs[k] === 'bigint') obj.create_attribute(k, attrs[k], null, '<q'); else obj.create_attribute(k, attrs[k]); } };

/** A NetCDF-4 coordinate variable: a 1-D dataset that is also the dimension scale of its dimension. */
function axis(f, name, values, attrs = {}) {
  const d = f.create_dataset({ name, data: f64(values) });
  d.make_scale(name);
  setAttrs(d, attrs);
  return d;
}

/** A NetCDF-4 variable on named dimensions. */
function variable(f, name, data, shape, dims, attrs = {}, extra = {}) {
  const d = f.create_dataset({ name, data, shape, ...extra });
  dims.forEach((dim, i) => { if (dim) d.attach_scale(i, dim); });
  setAttrs(d, attrs);
  return d;
}

const LAT = { standard_name: 'latitude', units: 'degrees_north', axis: 'Y' };
const LON = { standard_name: 'longitude', units: 'degrees_east', axis: 'X' };
const wgs84 = (f) => setAttrs(f.create_dataset({ name: 'crs', data: new Int32Array([0]), shape: [] }), { grid_mapping_name: 'latitude_longitude', epsg_code: 'EPSG:4326' });

test('vendored HDF5 library matches the record in docs/vendor.md (size and SHA-256) and reports itself available', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  let total = 0;
  for (const f of ['hdf5_hl.js', 'hdf5_util.js', 'LICENSE.txt']) {
    const buf = await readFile(new URL(`../../vendor/h5wasm/${f}`, import.meta.url));
    if (f !== 'LICENSE.txt') total += buf.length;
    const sha = createHash('sha256').update(buf).digest('hex');
    assert.ok(doc.includes(sha), `${f}: SHA-256 ${sha} is recorded`);
    assert.ok(doc.includes(String(buf.length)), `${f}: size ${buf.length} is recorded`);
  }
  assert.equal(hdf5.sizeBytes, total);
  assert.equal(hdf5.id, 'hdf5');
  assert.equal(hdf5.library.name, 'h5wasm'); assert.equal(hdf5.library.version, '0.10.3');
  assert.ok(doc.includes('h5wasm') && doc.includes('0.10.3'));
  for (const f of ['hdf5', 'netcdf4', 'med', 'cgns', 'exodus-hdf5']) assert.ok(hdf5.formats.includes(f), f);
  assert.equal(hdf5.available(), true);
  const first = await hdf5.load();
  assert.equal(typeof first.h5.File, 'function');
  assert.equal(await hdf5.load(), first);
});

test('GEBCO layout: int16 elevation(lat, lon) with _FillValue and a WGS 84 grid mapping → grid in local metres', async () => {
  const nx = 12, ny = 8;
  const z = new Int16Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) z[j * nx + i] = -100 - 2 * j + i;
  z[3 * nx + 4] = -32767;
  const bytes = write((f) => {
    setAttrs(f, { Conventions: 'CF-1.6', title: 'The GEBCO_2023 Grid - a continuous terrain model for oceans and land' });
    axis(f, 'lat', Array.from({ length: ny }, (_, j) => 56 + (j + 0.5) / 240), LAT);
    axis(f, 'lon', Array.from({ length: nx }, (_, i) => 2 + (i + 0.5) / 240), LON);
    wgs84(f);
    variable(f, 'elevation', z, [ny, nx], ['lat', 'lon'], { standard_name: 'height_above_mean_sea_level', long_name: 'Elevation relative to sea level', units: 'm', grid_mapping: 'crs' }, { chunks: [4, 6], compression: 'gzip', compression_opts: 4 });
    f.get('elevation').create_attribute('_FillValue', new Int16Array([-32767]));
  });
  assert.equal(hdf5.sniff({ bytes }), true);
  const r = await hdf5.parse({ name: 'tiles/gebco_2023_n57_e002.nc', bytes });
  assert.deepEqual(r.warnings, []);
  assert.equal(r.models.length, 1);
  const g = r.models[0];
  assert.equal(g.kind, 'grid'); assert.equal(g.name, 'gebco_2023_n57_e002');
  assert.equal(g.nx, nx); assert.equal(g.ny, ny);
  assert.ok(g.z instanceof Float32Array);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { if (j === 3 && i === 4) assert.ok(Number.isNaN(g.z[j * nx + i])); else assert.equal(g.z[j * nx + i], -100 - 2 * j + i); }
  assert.equal(g.crs, 'EPSG:4326');
  // 15 arc-seconds: about 463.3 m of latitude; a parallel at 56.017° N is cos(lat) times shorter
  rel(g.dy, (Math.PI / 180 / 240) * 6371008.8, 5e-3);
  rel(g.dx / g.dy, Math.cos((56 + 4 / 240) * Math.PI / 180), 5e-3);
  close(g.x0, -((nx - 1) / 2) * g.dx, 1e-6); close(g.y0, -((ny - 1) / 2) * g.dy, 1e-6);
  assert.equal(r.meta.format, 'netcdf4');
  assert.equal(r.meta.flavour, 'netcdf4-cf-grid');
  assert.equal(r.meta.variable, 'elevation'); assert.equal(r.meta.convention, 'CF-1.6');
  assert.deepEqual(r.meta.axes, { x: 'lon', y: 'lat' });
  assert.equal(r.meta.crs.epsg, 4326);
  close(r.meta.crs.origin.lon, 2 + 6 / 240, 1e-9); close(r.meta.crs.origin.lat, 56 + 4 / 240, 1e-9);
  assert.deepEqual(r.meta.skipped, { noDataCells: 1 });
  assert.deepEqual(r.meta.decimation, { factor: 1, source: [nx, ny] });
  assert.equal(r.meta.packing, null);
  assert.deepEqual(r.meta.tree.entries.map((e) => e.path).sort(), ['/crs', '/elevation', '/lat', '/lon']);
  assert.deepEqual(r.meta.tree.entries.find((e) => e.path === '/elevation'), { path: '/elevation', type: 'dataset', shape: [ny, nx], dtype: '<h' });
  assert.ok(inspect(g).ok);
});

test('ETOPO layout: packed z(lat, lon), latitude stored north to south, block-averaged to maxCells', async () => {
  const nx = 8, ny = 6;
  const raw = new Int16Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) raw[j * nx + i] = 100 * j + i;
  const bytes = write((f) => {
    axis(f, 'lat', Array.from({ length: ny }, (_, j) => 10 - 0.1 * j), { long_name: 'latitude', units: 'degrees_north' });
    axis(f, 'lon', Array.from({ length: nx }, (_, i) => 20 + 0.1 * i), { long_name: 'longitude', units: 'degrees_east' });
    wgs84(f);
    variable(f, 'z', raw, [ny, nx], ['lat', 'lon'], { long_name: 'z', grid_mapping: 'crs', scale_factor: 0.5, add_offset: -1000 });
  });
  // every value is −1000 + 0.5 · (100 j + i); row 0 of the model is the southern row, j = 5 of the file
  const full = await hdf5.parse({ name: 'ETOPO_2022_v1_60s.nc', bytes }, { projection: 'none' });
  const g = full.models[0];
  assert.equal(g.nx, nx); assert.equal(g.ny, ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) assert.equal(g.z[(ny - 1 - j) * nx + i], -1000 + 0.5 * (100 * j + i));
  close(g.x0, 20, 1e-12); close(g.y0, 9.5, 1e-12); close(g.dx, 0.1, 1e-12); close(g.dy, 0.1, 1e-12);
  assert.deepEqual(full.meta.packing, { scaleFactor: 0.5, addOffset: -1000 });
  assert.equal(full.meta.variable, 'z');
  // 48 nodes into at most 12: blocks of 2 × 2, node positions at the block centres
  const r = await hdf5.parse({ name: 'ETOPO_2022_v1_60s.nc', bytes }, { projection: 'none', maxCells: 12 });
  const b = r.models[0];
  assert.equal(b.nx, 4); assert.equal(b.ny, 3);
  close(b.x0, 20.05, 1e-12); close(b.y0, 9.55, 1e-12); close(b.dx, 0.2, 1e-12); close(b.dy, 0.2, 1e-12);
  assert.equal(b.z[0], -1000 + 0.5 * (500 + 501 + 400 + 401) / 4);
  assert.equal(b.z[2 * 4 + 3], -1000 + 0.5 * (106 + 107 + 6 + 7) / 4);
  assert.deepEqual(r.meta.decimation, { factor: 2, source: [nx, ny] });
  assert.deepEqual(r.warnings, ['The raster (8 × 6 cells) was block-averaged by 2 in each direction to 4 × 3 to stay under maxCells.']);
});

test('default limit of 1.5 million nodes: a 1300 × 1200 compressed grid is read in strips and halved', async () => {
  const nx = 1300, ny = 1200;
  const z = new Float32Array(nx * ny);
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) z[j * nx + i] = i + 10 * j;
  const bytes = write((f) => {
    axis(f, 'y', Array.from({ length: ny }, (_, j) => 2000 + 5 * j), { standard_name: 'projection_y_coordinate', units: 'm' });
    axis(f, 'x', Array.from({ length: nx }, (_, i) => 1000 + 5 * i), { standard_name: 'projection_x_coordinate', units: 'm' });
    variable(f, 'elevation', z, [ny, nx], ['y', 'x'], { units: 'm' }, { chunks: [100, 130], compression: 'gzip', compression_opts: 1 });
  });
  const r = await hdf5.parse({ name: 'site.nc', bytes });
  const g = r.models[0];
  assert.equal(g.nx, 650); assert.equal(g.ny, 600);
  assert.ok(g.nx * g.ny <= 1500000);
  close(g.x0, 1002.5, 1e-9); close(g.y0, 2002.5, 1e-9); close(g.dx, 10, 1e-9); close(g.dy, 10, 1e-9);
  // mean of i + 10 j over a 2 × 2 block
  assert.equal(g.z[0], 5.5); assert.equal(g.z[1], 7.5); assert.equal(g.z[650], 25.5); assert.equal(g.z[650 * 600 - 1], 1298.5 + 10 * 1198.5);
  assert.equal(r.meta.decimation.factor, 2);
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  const all = await hdf5.parse({ name: 'site.nc', bytes }, { maxCells: 0 });
  assert.equal(all.models[0].nx, nx); assert.equal(all.models[0].ny, ny);
  assert.equal(all.models[0].z[3 * nx + 7], 37);
  close(all.models[0].x0, 1000, 1e-9);
});

test('EMODnet naming: elevation among flag layers; older COLUMNS / LINES tiles with DEPTH positive down', async () => {
  const nx = 5, ny = 4;
  const el = Float32Array.from({ length: nx * ny }, (_, k) => -50 - k);
  const tile = write((f) => {
    axis(f, 'lat', [52, 52.01, 52.02, 52.03], LAT);
    axis(f, 'lon', [3, 3.01, 3.02, 3.03, 3.04], LON);
    wgs84(f);
    variable(f, 'value_count', new Int32Array(nx * ny).fill(3), [ny, nx], ['lat', 'lon'], {});
    variable(f, 'interpolation_flag', new Int8Array(nx * ny), [ny, nx], ['lat', 'lon'], {});
    variable(f, 'elevation', el, [ny, nx], ['lat', 'lon'], { units: 'm', grid_mapping: 'crs', missing_value: -50 });
  });
  const r = await hdf5.parse({ name: 'E4_2022.nc', bytes: tile });
  assert.equal(r.meta.variable, 'elevation');
  assert.ok(r.warnings.some((w) => /Several gridded variables are present \(.*value_count.*\); "elevation" was used/.test(w)));
  assert.ok(Number.isNaN(r.models[0].z[0]), 'missing_value is no data');
  assert.equal(r.models[0].z[7], -57);
  const flags = await hdf5.parse({ name: 'E4_2022.nc', bytes: tile }, { variable: 'value_count' });
  assert.equal(flags.models[0].z[7], 3);
  await assert.rejects(() => hdf5.parse({ name: 'E4_2022.nc', bytes: tile }, { variable: 'nope' }), /variable "nope" is not a gridded variable.*Gridded variables: .*elevation/);

  const old = write((f) => {
    axis(f, 'LINES', [52, 52.01, 52.02, 52.03], { units: 'degrees_north' });
    axis(f, 'COLUMNS', [3, 3.01, 3.02, 3.03, 3.04], { units: 'degrees_east' });
    variable(f, 'DEPTH_SMOOTH', Float32Array.from(el, (v) => -v + 1), [ny, nx], ['LINES', 'COLUMNS'], {});
    variable(f, 'DEPTH', Float32Array.from(el, (v) => -v), [ny, nx], ['LINES', 'COLUMNS'], { units: 'm' });
  });
  const o = await hdf5.parse({ name: 'D5.dtm.nc', bytes: old });
  assert.equal(o.meta.variable, 'DEPTH');
  assert.deepEqual(o.meta.axes, { x: 'COLUMNS', y: 'LINES' });
  assert.equal(o.meta.vertical.sign, -1, 'depths are turned into elevations');
  assert.equal(o.models[0].z[7], -57);
  assert.equal(o.models[0].crs, 'EPSG:4326');
});

test('CF details: record dimension, valid_min, axes without dimension scales, uneven axes returned as points', async () => {
  const nx = 4, ny = 5;
  const raw = Int16Array.from({ length: 2 * nx * ny }, (_, k) => k);
  raw[2 * nx + 1] = -9999;
  const bytes = write((f) => {
    axis(f, 'time', [0, 1]);
    axis(f, 'y', [1000, 990, 980, 970, 960], { standard_name: 'projection_y_coordinate', units: 'm' });
    axis(f, 'x', [500, 525, 550, 575], { standard_name: 'projection_x_coordinate', units: 'm' });
    variable(f, 'depth', raw, [2, ny, nx], ['time', 'y', 'x'], { scale_factor: 0.1, add_offset: 50, positive: 'down', units: 'm', valid_min: 1 });
    f.get('depth').create_attribute('_FillValue', new Int16Array([-9999]));
  });
  const r = await hdf5.parse({ name: 'survey.nc4', bytes });
  const g = r.models[0];
  assert.equal(g.nx, nx); assert.equal(g.ny, ny);
  close(g.x0, 500, 1e-9); close(g.y0, 960, 1e-9); close(g.dx, 25, 1e-9); close(g.dy, 10, 1e-9);
  // first time step; y runs north to south in the file; depth 50 + 0.1 · raw becomes a negative elevation
  close(g.z[0], -(50 + 0.1 * 16), 1e-5); close(g.z[4 * nx + 3], -(50 + 0.1 * 3), 1e-5);
  assert.ok(Number.isNaN(g.z[4 * nx]), 'raw 0 is below valid_min');
  assert.ok(Number.isNaN(g.z[2 * nx + 1]), '_FillValue');
  assert.deepEqual(r.meta.skipped, { noDataCells: 2 });
  assert.ok(r.warnings.some((w) => /"depth" has leading dimension\(s\) of size 2; the first index of each was used/.test(w)));
  assert.equal(r.meta.vertical.sign, -1);

  // plain HDF5 with CF names but no dimension scales: the axes are found by name and length
  const plain = write((f) => {
    f.create_dataset({ name: 'northing', data: f64([0, 2, 4]) });
    f.create_dataset({ name: 'easting', data: f64([10, 11, 12, 13]) });
    f.create_dataset({ name: 'bathymetry', data: Float64Array.from({ length: 12 }, (_, k) => -k), shape: [3, 4] });
  });
  const p = await hdf5.parse({ name: 'grid.h5', bytes: plain });
  assert.deepEqual([p.models[0].nx, p.models[0].ny, p.models[0].dx, p.models[0].dy, p.models[0].x0, p.models[0].y0], [4, 3, 1, 2, 10, 0]);
  assert.equal(p.models[0].z[5], -5);
  assert.deepEqual(p.meta.axes, { x: 'easting', y: 'northing' });

  const uneven = write((f) => {
    axis(f, 'y', [0, 1, 3], { units: 'm' });
    axis(f, 'x', [0, 10, 20, 30], { units: 'm' });
    variable(f, 'elevation', Float32Array.from({ length: 12 }, (_, k) => k), [3, 4], ['y', 'x'], {});
    f.get('elevation').create_attribute('_FillValue', new Float32Array([5]));
  });
  const u = await hdf5.parse({ name: 'uneven.nc', bytes: uneven });
  assert.equal(u.models[0].kind, 'points');
  assert.equal(u.models[0].positions.length / 3, 11);
  assert.ok(u.warnings.some((w) => /not evenly spaced, so the grid nodes are returned as points/.test(w)));
  const q = u.models[0], o = q.origin || [0, 0, 0];
  // node 6 of the file (x = 20, y = 1, z = 6) is the sixth point kept, the fill value before it being dropped
  close(o[0] + q.positions[15], 20, 1e-6); close(o[1] + q.positions[16], 1, 1e-6); close(o[2] + q.positions[17], 6, 1e-6);
});

test('UGRID 2-D mesh on NetCDF-4: triangles and quadrilaterals → surface with node elevations', async () => {
  const bytes = write((f) => {
    setAttrs(f, { Conventions: 'CF-1.8 UGRID-1.0' });
    setAttrs(f.create_dataset({ name: 'mesh2d', data: new Int32Array([0]), shape: [] }), { cf_role: 'mesh_topology', topology_dimension: 2, node_coordinates: 'node_x node_y', face_node_connectivity: 'face_nodes' });
    setAttrs(f.create_dataset({ name: 'node_x', data: f64([0, 10, 10, 0, 20]) }), { units: 'm' });
    setAttrs(f.create_dataset({ name: 'node_y', data: f64([0, 0, 10, 10, 5]) }), { units: 'm' });
    setAttrs(f.create_dataset({ name: 'node_z', data: f64([5, 6, 7, 8, 9]) }), { standard_name: 'sea_floor_depth_below_geoid', positive: 'down' });
    const fn = f.create_dataset({ name: 'face_nodes', data: new Int32Array([1, 2, 3, 4, 2, 5, 3, -999]), shape: [2, 4] });
    setAttrs(fn, { cf_role: 'face_node_connectivity', start_index: 1 });
    fn.create_attribute('_FillValue', new Int32Array([-999]));
  });
  const r = await hdf5.parse({ name: 'flow.nc', bytes });
  assert.equal(r.meta.flavour, 'netcdf4-ugrid');
  assert.deepEqual(r.meta.elementTypes, { quad4: 1, tri3: 1 });
  const m = r.models[0];
  assert.equal(m.kind, 'surface'); assert.equal(m.name, 'mesh2d');
  assert.deepEqual(Array.from(m.indices), [0, 1, 2, 0, 2, 3, 1, 4, 2]);
  const o = m.origin || [0, 0, 0];
  assert.deepEqual([0, 1, 2, 3, 4].map((k) => o[2] + m.positions[3 * k + 2]), [-5, -6, -7, -8, -9]);
  close(o[0] + m.positions[12], 20, 1e-6);
  // plan area: the square of 10 × 10 and the triangle of base 10 and height 10 beside it
  let plan = 0;
  for (let t = 0; t < m.indices.length; t += 3) { const [a, b, c] = [0, 1, 2].map((k) => m.indices[t + k]); plan += 0.5 * ((m.positions[3 * b] - m.positions[3 * a]) * (m.positions[3 * c + 1] - m.positions[3 * a + 1]) - (m.positions[3 * c] - m.positions[3 * a]) * (m.positions[3 * b + 1] - m.positions[3 * a + 1])); }
  close(plan, 150, 1e-6);
  assert.equal(r.meta.nodeElevation, 'node_z');
});

// The four solid cells of the reference export: a unit cube (nodes 1–8), a pyramid on its top face (apex 9),
// a tetrahedron (10–13) and a wedge (14–19), each positively oriented in the VTK / Gmsh order.
const CELL_XYZ = [
  [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0.5, 0.5, 2],
  [2, 0, 0], [3, 0, 0], [2, 1, 0], [2, 0, 1],
  [4, 0, 0], [5, 0, 0], [4, 1, 0], [4, 0, 2], [5, 0, 2], [4, 1, 2],
];
const CANONICAL = { hex8: [1, 2, 3, 4, 5, 6, 7, 8], pyr5: [5, 6, 7, 8, 9], tet4: [10, 11, 12, 13], wedge6: [14, 15, 16, 17, 18, 19] };
const CELL_VOLUME = { hex8: 1, pyr5: 1 / 3, tet4: 1 / 6, wedge6: 1 };

/** Check a volume model against the four reference cells: source node order and positive volumes. */
function checkCells(vol) {
  for (const type of Object.keys(CANONICAL)) {
    const blk = vol.cells.find((c) => c.type === type);
    assert.ok(blk, `${type} block`);
    const n = CANONICAL[type].length;
    assert.deepEqual(Array.from(blk.conn.subarray(0, n), (k) => vol.nodeIds[k]), CANONICAL[type], `${type}: nodes in the order of this application`);
    rel(cellVolume(vol.positions, blk.conn, 0, type), CELL_VOLUME[type], 1e-9, `${type} volume`);
  }
}

/** MED 4 tree as written by the MED library: non-interlaced coordinates and connectivities, 64-bit integers. */
function medFile({ withStep = true } = {}) {
  return write((f) => {
    setAttrs(f.create_group('INFOS_GENERALES'), { MAJ: 4n, MIN: 1n, REL: 0n });
    const mesh = f.create_group('ENS_MAA').create_group('platform');
    setAttrs(mesh, { DES: 'four cells', DIM: 3n, ESP: 3n, NOM: '', TYP: 0n, UNI: '', UNT: '' });
    const step = withStep ? mesh.create_group('-0000000000000000001-0000000000000000001') : mesh;
    const noe = step.create_group('NOE');
    const n = CELL_XYZ.length;
    setAttrs(noe.create_dataset({ name: 'COO', data: f64([0, 1, 2].flatMap((k) => CELL_XYZ.map((p) => p[k]))) }), { NBR: BigInt(n) });
    noe.create_dataset({ name: 'FAM', data: i64(CELL_XYZ.map((_, i) => (i < 4 ? 1 : 0))) });
    const mai = step.create_group('MAI');
    // connectivities exactly as exported for these cells (MED node order), 1-based
    const cells = {
      HE8: [[[1, 4, 3, 2, 5, 8, 7, 6]], [-1]],
      PY5: [[[5, 8, 7, 6, 9]], [-1]],
      TE4: [[[10, 12, 11, 13]], [-2]],
      PE6: [[[14, 16, 15, 17, 19, 18]], [-2]],
      TR3: [[[10, 12, 11], [10, 11, 13]], [-3, -4]],
      QU4: [[[1, 4, 3, 2]], [-3]],
      SE2: [[[1, 2], [2, 3]], [0, -5]],
    };
    for (const t in cells) {
      const [conn, fam] = cells[t];
      const g = mai.create_group(t);
      const npe = conn[0].length;
      // node k of every cell before node k + 1
      const nod = [];
      for (let k = 0; k < npe; k++) for (const c of conn) nod.push(c[k]);
      setAttrs(g.create_dataset({ name: 'NOD', data: i64(nod) }), { NBR: BigInt(conn.length) });
      g.create_dataset({ name: 'FAM', data: i64(fam) });
    }
    const fas = f.create_group('FAS').create_group('platform');
    setAttrs(fas.create_group('FAMILLE_ZERO'), { NUM: 0n });
    const family = (kind, name, num, groups) => {
      const parent = fas.keys().includes(kind) ? fas.get(kind) : fas.create_group(kind);
      const g = parent.create_group(name);
      setAttrs(g, { NUM: BigInt(num) });
      if (!groups.length) return;
      const gro = g.create_group('GRO');
      setAttrs(gro, { NBR: BigInt(groups.length) });
      const nom = new Int8Array(80 * groups.length);
      groups.forEach((s, i) => nom.set(chars(s, 80), 80 * i));
      gro.create_dataset({ name: 'NOM', data: nom, shape: [groups.length, 80] });
    };
    family('ELEME', 'F_3D_1', -1, ['steel']);
    family('ELEME', 'F_3D_2', -2, ['concrete', 'ballast tanks']);
    family('ELEME', 'F_2D_3', -3, ['skin']);
    family('ELEME', 'F_2D_4', -4, ['skin']);
    family('ELEME', 'F_1D_5', -5, []);
    family('NOEUD', 'F_0D_1', 1, ['fixed']);
  });
}

test('MED mesh: SE2, TR3, QU4, TE4, HE8, PE6, PY5 with the documented node-order mapping; families → sets', async () => {
  const r = await hdf5.parse({ name: 'platform.med', bytes: medFile() });
  assert.equal(r.meta.flavour, 'med'); assert.equal(r.meta.format, 'med');
  assert.equal(r.meta.medVersion, '4.1.0'); assert.equal(r.meta.mesh, 'platform');
  assert.deepEqual(r.meta.medTypes, { HE8: 1, PE6: 1, PY5: 1, QU4: 1, SE2: 2, TE4: 1, TR3: 2 });
  assert.deepEqual(r.meta.elementTypes, { hex8: 1, wedge6: 1, pyr5: 1, quad4: 1, line2: 2, tet4: 1, tri3: 2 });
  assert.deepEqual(r.meta.elementSets.sort(), ['ballast tanks', 'concrete', 'skin', 'steel']);
  assert.deepEqual(r.meta.nodeSets, ['fixed']);
  assert.deepEqual(r.meta.families.find((q) => q.number === -2), { number: -2, name: 'F_3D_2', groups: ['concrete', 'ballast tanks'] });
  assert.deepEqual(r.models.map((m) => m.kind), ['frame', 'surface', 'volume']);
  const [frame, shell, vol] = r.models;

  // solids: every cell comes out in the order of this application and with a positive volume
  checkCells(vol);
  const i = inspect(vol);
  assert.ok(i.ok, JSON.stringify(i.issues));
  assert.ok(!i.issues.some((q) => q.code === 'inverted-cells'));
  rel(i.stats.volume, 1 + 1 / 3 + 1 / 6 + 1, 1e-9);
  assert.deepEqual(Object.fromEntries(Object.entries(vol.sets.cells).map(([k, v]) => [k, v.length])), { steel: 2, concrete: 2, 'ballast tanks': 2 });
  assert.deepEqual(Object.values(vol.tagNames).sort(), ['concrete + ballast tanks', 'steel']);
  assert.deepEqual(Array.from(vol.sets.nodes.fixed, (k) => vol.nodeIds[k]).sort(), [1, 2, 3, 4]);

  // shells: the two families of the group "skin" are one surface group
  assert.equal(shell.indices.length / 3, 4);
  assert.deepEqual(shell.groups.map((g) => [g.name, g.count / 3]), [['skin', 4]]);
  assert.deepEqual(shell.elements.map((e) => [e.type, e.ids.length]), [['quad4', 1], ['tri3', 2]]);
  // TR3 and QU4 keep the MED order: 10-12-11 lies in z = 0 seen from below
  const tri = shell.elements.find((e) => e.type === 'tri3');
  assert.deepEqual(Array.from(tri.conn.subarray(0, 3), (k) => shell.nodeIds[k]), [10, 12, 11]);
  assert.equal(shell.sets.elements.skin.length, 3);

  // beams: one member in family zero, one in a family without a group
  assert.equal(frame.members.length, 2);
  assert.deepEqual(frame.members.map((m) => [frame.nodeIds[m.a], frame.nodeIds[m.b], m.group]), [[1, 2, '0'], [2, 3, 'F_1D_5']]);
  assert.ok(r.warnings.every((w) => /no length unit/.test(w)), r.warnings.join(' | '));
});

test('MED variants: mesh without a computation-step group, 2-D space, quadratic cells on their corners, polygons reported', async () => {
  const old = await hdf5.parse({ name: 'old.med', bytes: medFile({ withStep: false }) });
  checkCells(old.models.find((m) => m.kind === 'volume'));

  const bytes = write((f) => {
    const mesh = f.create_group('ENS_MAA').create_group('plate');
    setAttrs(mesh, { DIM: 2n, ESP: 2n, TYP: 0n, UNI: `${'mm'.padEnd(16)}${'mm'.padEnd(16)}` });
    const step = mesh.create_group('-0000000000000000001-0000000000000000001');
    // six nodes of one TR6 (corners, then mid-sides) and four more of a T10's corners are not needed here
    step.create_group('NOE').create_dataset({ name: 'COO', data: f64([0, 2000, 0, 1000, 1000, 0, /* y */ 0, 0, 2000, 0, 1000, 1000]) });
    const mai = step.create_group('MAI');
    mai.create_group('TR6').create_dataset({ name: 'NOD', data: i64([1, 2, 3, 4, 5, 6]) });
    const pog = mai.create_group('POG');
    pog.create_dataset({ name: 'NOD', data: i64([1, 4, 2, 5, 3]) });
    pog.create_dataset({ name: 'INN', data: i64([1, 6]) });
  });
  const r = await hdf5.parse({ name: 'plate.med', bytes });
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.equal(r.meta.spaceDimension, 2);
  const m = r.models[0];
  assert.equal(m.kind, 'surface');
  assert.equal(m.indices.length, 3);
  rel(inspect(m).stats.area, 0.5 * 2 * 2, 1e-9);
  assert.ok(r.warnings.some((w) => /Quadratic MED cells \(1 TR6\) were read on their corner nodes/.test(w)));
  assert.ok(r.warnings.some((w) => /Skipped 1 record\(s\) of unsupported type MED cell type POG \(polygons\)/.test(w)), r.warnings.join(' | '));

  const two = write((f) => {
    const ens = f.create_group('ENS_MAA');
    for (const name of ['a', 'b']) {
      const noe = ens.create_group(name).create_group('NOE');
      noe.create_dataset({ name: 'COO', data: f64(name === 'a' ? [0, 1, 0, 0, 0, 1, 0, 0, 0] : [0, 5, 0, 0, 0, 5, 0, 0, 0]) });
      ens.get(name).create_group('MAI').create_group('TR3').create_dataset({ name: 'NOD', data: i64([1, 2, 3]) });
    }
  });
  const first = await hdf5.parse({ name: 'two.med', bytes: two });
  assert.ok(first.warnings.some((w) => /holds 2 meshes \(a, b\); "a" was read/.test(w)));
  const second = await hdf5.parse({ name: 'two.med', bytes: two }, { mesh: 'b' });
  rel(inspect(second.models[0]).stats.area, 12.5, 1e-9);
  await assert.rejects(() => hdf5.parse({ name: 'two.med', bytes: two }, { mesh: 'c' }), /MED mesh "c" is not in the file\. Meshes: a, b/);
});

/** CGNS node: a group with name / label / type attributes and its payload in the dataset " data". */
function cg(parent, name, label, type, data, shape) {
  const g = parent.create_group(name);
  setAttrs(g, { name, label, type });
  if (data) g.create_dataset(shape ? { name: ' data', data, shape } : { name: ' data', data });
  return g;
}

function cgnsRoot(f) {
  setAttrs(f, { name: 'HDF5 MotherNode', label: 'Root Node of HDF5 File', type: 'MT' });
  f.create_dataset({ name: ' format', data: chars('IEEE_LITTLE_32') });
  f.create_dataset({ name: ' hdf5version', data: chars('HDF5 Version 1.14.6') });
  cg(f, 'CGNSLibraryVersion', 'CGNSLibraryVersion_t', 'R4', new Float32Array([3.4]));
}

function coordinates(zone, xyz) {
  const gc = cg(zone, 'GridCoordinates', 'GridCoordinates_t', 'MT');
  ['X', 'Y', 'Z'].forEach((a, k) => cg(gc, `Coordinate${a}`, 'DataArray_t', 'R8', f64(xyz.map((p) => p[k]))));
}

function section(zone, name, type, range, conn) {
  const s = cg(zone, name, 'Elements_t', 'I4', new Int32Array([type, 0]));
  cg(s, 'ElementRange', 'IndexRange_t', 'I4', new Int32Array(range));
  cg(s, 'ElementConnectivity', 'DataArray_t', 'I4', new Int32Array(conn));
  return s;
}

test('CGNS: MIXED and typed sections, boundary conditions → sets, a structured zone → hexahedra, units from the file', async () => {
  const ni = 3, nj = 2, nk = 2;
  const lattice = [];
  for (let k = 0; k < nk; k++) for (let j = 0; j < nj; j++) for (let i = 0; i < ni; i++) lattice.push([10000 + 1000 * i, 2000 * j, 3000 * k]);
  const bytes = write((f) => {
    cgnsRoot(f);
    const base = cg(f, 'Base', 'CGNSBase_t', 'I4', new Int32Array([3, 3]));
    // mass, length, time, temperature, angle — 32 characters each
    cg(base, 'DimensionalUnits', 'DimensionalUnits_t', 'C1', Int8Array.from([...chars('Kilogram', 32), ...chars('Millimeter', 32), ...chars('Second', 32), ...chars('Kelvin', 32), ...chars('Degree', 32)]), [5, 32]);
    const fam = cg(base, 'F_inlet', 'Family_t', 'MT');
    cg(fam, 'inflow plane', 'FamilyName_t', 'C1', chars('inflow plane'));
    const zone = cg(base, 'zone1 cells', 'Zone_t', 'I4', new Int32Array([CELL_XYZ.length, 4, 0]), [3, 1]);
    cg(zone, 'ZoneType', 'ZoneType_t', 'C1', chars('Unstructured'));
    coordinates(zone, CELL_XYZ.map((p) => p.map((v) => 1000 * v)));
    // MIXED: HEXA_8 = 17, PYRA_5 = 12, TETRA_4 = 10, PENTA_6 = 14, each followed by its nodes
    section(zone, 'solids', 20, [1, 4], [17, ...CANONICAL.hex8, 12, ...CANONICAL.pyr5, 10, ...CANONICAL.tet4, 14, ...CANONICAL.wedge6]);
    section(zone, 'skin', 7, [5, 6], [1, 4, 3, 2, 5, 6, 7, 8]);
    section(zone, 'legs', 3, [7, 7], [1, 5]);
    section(zone, 'faces', 22, [8, 9], [1, 2, 3, 2, 3, 4]);
    const zbc = cg(zone, 'ZoneBC', 'ZoneBC_t', 'MT');
    const inlet = cg(zbc, 'inlet', 'BC_t', 'C1', chars('FamilySpecified'));
    cg(inlet, 'GridLocation', 'GridLocation_t', 'C1', chars('FaceCenter'));
    cg(inlet, 'PointList', 'IndexArray_t', 'I4', new Int32Array([6]), [1, 1]);
    cg(inlet, 'FamilyName', 'FamilyName_t', 'C1', chars('F_inlet'));
    const pinned = cg(zbc, 'pinned', 'BC_t', 'C1', chars('BCWall'));
    cg(pinned, 'PointRange', 'IndexRange_t', 'I4', new Int32Array([1, 4]), [2, 1]);

    const block = cg(base, 'zone2 block', 'Zone_t', 'I4', new Int32Array([ni, nj, nk, ni - 1, nj - 1, nk - 1, 0, 0, 0]), [3, 3]);
    cg(block, 'ZoneType', 'ZoneType_t', 'C1', chars('Structured'));
    coordinates(block, lattice);
  });
  const r = await hdf5.parse({ name: 'jacket.cgns', bytes });
  assert.equal(r.meta.flavour, 'cgns'); assert.equal(r.meta.cgnsVersion, 3.4);
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.deepEqual(r.meta.zones.map((z) => [z.zone, z.type, z.vertices]), [['zone1 cells', 'Unstructured', 19], ['zone2 block', 'Structured', 12]]);
  assert.equal(r.meta.structuredZones, 1);
  assert.deepEqual(r.meta.elementTypes, { hex8: 3, pyr5: 1, tet4: 1, wedge6: 1, quad4: 2, line2: 1 });
  assert.deepEqual(r.models.map((m) => m.kind), ['frame', 'surface', 'volume']);
  const [frame, shell, vol] = r.models;
  // coordinates are in millimetres in the file: the unit cube of 1000 is 1 m³
  checkCells(vol);
  const i = inspect(vol);
  assert.ok(i.ok, JSON.stringify(i.issues));
  // the structured zone: two right-handed hexahedra of 1 × 2 × 3 m beside the unstructured cells
  const hex = vol.cells.find((c) => c.type === 'hex8');
  assert.equal(hex.ids.length, 3);
  rel(cellVolume(vol.positions, hex.conn, 8, 'hex8'), 6, 1e-9); rel(cellVolume(vol.positions, hex.conn, 16, 'hex8'), 6, 1e-9);
  rel(i.stats.volume, 2.5 + 12, 1e-9);
  // node numbers restart in the second zone and are kept apart
  assert.deepEqual(Array.from(hex.conn.subarray(8, 16), (k) => vol.nodeIds[k] - 19), [1, 2, 5, 4, 7, 8, 11, 10]);
  assert.deepEqual(Object.fromEntries(Object.entries(vol.sets.cells).map(([k, v]) => [k, v.length])), { solids: 4, 'zone2 block': 2 });
  assert.deepEqual(Array.from(vol.sets.nodes.pinned, (k) => vol.nodeIds[k]).sort(), [1, 2, 3, 4]);
  assert.deepEqual(shell.groups.map((g) => [g.name, g.count / 3]), [['skin', 4]]);
  // the boundary condition names element 6, the second quadrilateral; its family is known as "inflow plane"
  for (const name of ['inlet', 'F_inlet', 'inflow plane']) assert.deepEqual(Array.from(shell.sets.elements[name]), [1], name);
  assert.equal(frame.members.length, 1);
  assert.equal(frame.members[0].group, 'legs');
  assert.ok(r.warnings.some((w) => /Skipped 2 record\(s\) of unsupported type CGNS polyhedral section of NGON_n elements/.test(w)), r.warnings.join(' | '));
  assert.ok(r.warnings.some((w) => /2 zones were read side by side; nodes on zone interfaces are not merged/.test(w)));
});

test('CGNS HEXA_20: top-face and vertical mid-edge nodes are exchanged into the order of this application', async () => {
  const c = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
  const mid = (a, b) => c[a].map((v, k) => (v + c[b][k]) / 2);
  // SIDS order: corners, mid-edges of the bottom face, of the vertical edges, of the top face
  const xyz = [...c, mid(0, 1), mid(1, 2), mid(2, 3), mid(3, 0), mid(0, 4), mid(1, 5), mid(2, 6), mid(3, 7), mid(4, 5), mid(5, 6), mid(6, 7), mid(7, 4)];
  const bytes = write((f) => {
    cgnsRoot(f);
    const zone = cg(cg(f, 'Base', 'CGNSBase_t', 'I4', new Int32Array([3, 3])), 'z', 'Zone_t', 'I4', new Int32Array([20, 1, 0]), [3, 1]);
    cg(zone, 'ZoneType', 'ZoneType_t', 'C1', chars('Unstructured'));
    coordinates(zone, xyz);
    section(zone, 'brick', 18, [1, 1], Array.from({ length: 20 }, (_, k) => k + 1));
  });
  const r = await hdf5.parse({ name: 'brick.cgns', bytes });
  const blk = r.models[0].cells[0];
  assert.equal(blk.type, 'hex20');
  const zOf = (k) => r.models[0].positions[3 * blk.conn[k] + 2];
  assert.deepEqual([8, 9, 10, 11].map(zOf), [0, 0, 0, 0], 'bottom edges');
  assert.deepEqual([12, 13, 14, 15].map(zOf), [1, 1, 1, 1], 'top edges');
  assert.deepEqual([16, 17, 18, 19].map(zOf), [0.5, 0.5, 0.5, 0.5], 'vertical edges');
  rel(cellVolume(r.models[0].positions, blk.conn, 0, 'hex20'), 1, 1e-9);
  assert.ok(r.warnings.every((w) => /no length unit/.test(w)));
});

test('Exodus II on NetCDF-4: element blocks with names, node sets; side sets reported', async () => {
  const names = (list, w = 33) => { const out = []; for (const s of list) for (let k = 0; k < w; k++) out.push(s[k] || ''); return out; };
  const bytes = write((f) => {
    setAttrs(f, { title: 'two blocks', api_version: new Float32Array([8.11]) });
    f.create_dataset({ name: 'coordx', data: f64(CELL_XYZ.map((p) => p[0])) });
    f.create_dataset({ name: 'coordy', data: f64(CELL_XYZ.map((p) => p[1])) });
    f.create_dataset({ name: 'coordz', data: f64(CELL_XYZ.map((p) => p[2])) });
    f.create_dataset({ name: 'connect1', data: i64(CANONICAL.hex8), shape: [1, 8], dtype: '<q' }).create_attribute('elem_type', 'HEX8');
    f.create_dataset({ name: 'connect2', data: new Int32Array(CANONICAL.tet4), shape: [1, 4] }).create_attribute('elem_type', 'TETRA4');
    f.create_dataset({ name: 'connect3', data: new Int32Array(CANONICAL.wedge6), shape: [1, 6] }).create_attribute('elem_type', 'WEDGE6');
    f.create_dataset({ name: 'connect4', data: new Int32Array(CANONICAL.pyr5), shape: [1, 5] }).create_attribute('elem_type', 'PYRAMID5');
    f.create_dataset({ name: 'connect5', data: new Int32Array([1, 2, 3, 4, 5, 6, 7]), shape: [1, 7] }).create_attribute('elem_type', 'NSIDED');
    f.create_dataset({ name: 'eb_prop1', data: new Int32Array([10, 20, 30, 40, 50]) });
    f.create_dataset({ name: 'eb_names', data: names(['deck', 'brace', '', 'cap', 'odd']), shape: [5, 33], dtype: 'S1' });
    f.create_dataset({ name: 'node_ns1', data: new Int32Array([1, 2, 3, 4]) });
    f.create_dataset({ name: 'ns_names', data: names(['mudline']), shape: [1, 33], dtype: 'S1' });
    f.create_dataset({ name: 'elem_ss1', data: new Int32Array([1]) });
    f.create_dataset({ name: 'side_ss1', data: new Int32Array([5]) });
  });
  const r = await hdf5.parse({ name: 'model.exo', bytes });
  assert.equal(r.meta.flavour, 'exodus-netcdf4'); assert.equal(r.meta.format, 'exodus-hdf5');
  assert.equal(r.meta.title, 'two blocks'); assert.equal(r.meta.blocks, 5);
  const vol = r.models[0];
  assert.equal(vol.kind, 'volume');
  checkCells(vol);
  assert.ok(inspect(vol).ok);
  assert.deepEqual(Object.keys(vol.sets.cells).sort(), ['block 30', 'brace', 'cap', 'deck']);
  assert.deepEqual(vol.tagNames, { 10: 'deck', 20: 'brace', 30: 'block 30', 40: 'cap' });
  assert.deepEqual(Array.from(vol.sets.nodes.mudline, (k) => vol.nodeIds[k]), [1, 2, 3, 4]);
  assert.deepEqual(r.meta.skipped, { 'element type NSIDED with 7 nodes': 1, sideSets: 1 });
  assert.ok(r.warnings.some((w) => /1 side set\(s\) were not read/.test(w)));
});

test('the convention is decided from the object tree, not from the file name', async () => {
  const med = medFile();
  const viaNc = await importFile({ name: 'mesh.nc', bytes: med });
  assert.equal(viaNc.meta.via, 'adapter:hdf5');
  assert.equal(viaNc.meta.flavour, 'med');
  assert.deepEqual(viaNc.models.map((m) => m.kind), ['frame', 'surface', 'volume']);
  assert.equal((await importFile({ name: 'mesh.cgns', bytes: med })).meta.flavour, 'med');
  assert.equal((await importFile({ name: 'mesh.med', bytes: med })).format, 'med');
  const grid = write((f) => {
    axis(f, 'y', [0, 1, 2], { units: 'm' }); axis(f, 'x', [0, 1], { units: 'm' });
    variable(f, 'elevation', new Float32Array([1, 2, 3, 4, 5, 6]), [3, 2], ['y', 'x'], {});
  });
  for (const name of ['bathy.med', 'bathy.cgns', 'bathy.h5', 'bathy.nc', 'bathy.exo']) {
    const r = await importFile({ name, bytes: grid });
    assert.equal(r.meta.flavour, 'netcdf4-cf-grid', name);
    assert.equal(r.models[0].kind, 'grid', name);
  }
});

test('generic HDF5: nothing is guessed — the tree is listed; opts.datasets names a grid, points or a mesh', async () => {
  const bytes = write((f) => {
    const run = f.create_group('run_01');
    run.create_group('seabed').create_dataset({ name: 'height', data: Float32Array.from({ length: 12 }, (_, k) => -k), shape: [3, 4] });
    run.get('seabed').create_dataset({ name: 'east', data: f64([100, 110, 120, 130]) });
    run.get('seabed').create_dataset({ name: 'north', data: f64([40, 20, 0]) });
    run.create_dataset({ name: 'soundings', data: f64([0, 0, -5, 1, 0, -6, 0, 1, -7]), shape: [3, 3] });
    const mesh = run.create_group('mesh');
    mesh.create_dataset({ name: 'xyz', data: f64([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]), shape: [4, 3] });
    mesh.create_dataset({ name: 'tets', data: new Int32Array([1, 2, 3, 4]), shape: [1, 4] });
    mesh.create_dataset({ name: 'tris', data: new Int32Array([1, 3, 2]), shape: [1, 3] });
    f.create_dataset({ name: 'note', data: 'calibration run' });
  });
  const r = await hdf5.parse({ name: 'results.h5', bytes });
  assert.deepEqual(r.models, []);
  assert.equal(r.meta.flavour, 'generic');
  assert.equal(r.warnings.length, 1);
  assert.match(r.warnings[0], /No geometry convention was recognised.*holds 8 dataset\(s\): .*\/run_01\/seabed\/height \[3 × 4\].*opts\.datasets/);
  assert.deepEqual(r.meta.tree.entries.map((e) => `${e.type} ${e.path}`), [
    'dataset /note', 'group /run_01', 'group /run_01/mesh', 'dataset /run_01/mesh/tets', 'dataset /run_01/mesh/tris', 'dataset /run_01/mesh/xyz',
    'group /run_01/seabed', 'dataset /run_01/seabed/east', 'dataset /run_01/seabed/height', 'dataset /run_01/seabed/north', 'dataset /run_01/soundings',
  ]);
  assert.deepEqual(r.meta.tree.entries.find((e) => e.path === '/run_01/soundings'), { path: '/run_01/soundings', type: 'dataset', shape: [3, 3], dtype: '<d' });
  assert.equal(r.meta.tree.truncated, false);

  const g = await hdf5.parse({ name: 'results.h5', bytes }, { datasets: { grid: '/run_01/seabed/height', x: '/run_01/seabed/east', y: '/run_01/seabed/north' } });
  const grid = g.models[0];
  assert.deepEqual([grid.kind, grid.nx, grid.ny, grid.x0, grid.y0, grid.dx, grid.dy], ['grid', 4, 3, 100, 0, 10, 20]);
  // the northern row is first in the dataset and last in the model
  assert.deepEqual(Array.from(grid.z.subarray(8, 12)), [-0, -1, -2, -3]);
  assert.equal(g.meta.flavour, 'generic');
  const bare = await hdf5.parse({ name: 'results.h5', bytes }, { datasets: { grid: 'run_01/seabed/height' } });
  assert.ok(bare.warnings.some((w) => /unit spacing from the origin/.test(w)));

  const p = await hdf5.parse({ name: 'results.h5', bytes }, { datasets: { points: '/run_01/soundings' } });
  assert.equal(p.models[0].kind, 'points');
  assert.equal(p.models[0].positions.length, 9);

  const m = await hdf5.parse({ name: 'results.h5', bytes }, { datasets: { positions: '/run_01/mesh/xyz', elements: { tet4: '/run_01/mesh/tets', tri3: '/run_01/mesh/tris' } } });
  assert.deepEqual(m.models.map((q) => q.kind), ['surface', 'volume']);
  assert.equal(m.meta.indexBase, 1, 'indices 1 … n can only be 1-based');
  rel(cellVolume(m.models[1].positions, m.models[1].cells[0].conn, 0, 'tet4'), 1 / 6, 1e-9);

  await assert.rejects(() => hdf5.parse({ name: 'results.h5', bytes }, { datasets: { grid: '/nope' } }), /opts\.datasets\.grid = "\/nope" is not a dataset of the file\. Datasets: \/note, /);
  await assert.rejects(() => hdf5.parse({ name: 'results.h5', bytes }, { datasets: { grid: '/run_01/soundings', x: '/run_01/seabed/east', y: '/run_01/seabed/north' } }), /axis datasets must be evenly spaced and have 3 \(x\) and 3 \(y\) values/);
  await assert.rejects(() => hdf5.parse({ name: 'results.h5', bytes }, { datasets: { surface: '/x' } }), /opts\.datasets must name what to read/);
});

test('files that are not readable HDF5 fail with the reason; the library keeps working afterwards', async () => {
  const good = write((f) => {
    axis(f, 'y', [0, 1, 2], { units: 'm' }); axis(f, 'x', [0, 1], { units: 'm' });
    variable(f, 'elevation', Float32Array.from({ length: 6000 }, (_, k) => k), [3, 2000], [null, null], {}, { chunks: [1, 500], compression: 'gzip' });
    variable(f, 'z', new Float32Array([1, 2, 3, 4, 5, 6]), [3, 2], ['y', 'x'], {});
  });
  await assert.rejects(() => hdf5.parse({ name: 'notes.h5', bytes: new TextEncoder().encode('not a container') }), /notes\.h5: not an HDF5 file \(the HDF5 signature was not found\)\.$/);
  await assert.rejects(() => hdf5.parse({ name: 'classic.nc', bytes: Uint8Array.from([0x43, 0x44, 0x46, 2, 0, 0, 0, 0, 0, 0]) }), /It is a NetCDF classic file, which is read directly/);
  const adf = new Uint8Array(200); adf.set(new TextEncoder().encode('@(#)ADF Database Version B02012>'));
  await assert.rejects(() => hdf5.parse({ name: 'old.cgns', bytes: adf }), /CGNS file in the older ADF format: convert it with "cgnsconvert -h in\.cgns out\.cgns"/);
  await assert.rejects(() => hdf5.parse({ name: 'cut.nc', bytes: good.slice(0, Math.floor(good.length / 2)) }), /cut\.nc: the HDF5 file could not be opened \(.*truncated file.*\)\. The file is probably truncated or damaged/);
  await assert.rejects(() => hdf5.parse({ name: 'empty.nc', bytes: new Uint8Array(0) }), /not an HDF5 file/);
  assert.equal(hdf5.sniff({ bytes: good }), true);
  // the signature may follow a user block of 512, 1024, 2048 … bytes
  const shifted = new Uint8Array(1024 + 16); shifted.set(good.subarray(0, 8), 1024);
  assert.equal(hdf5.sniff({ bytes: shifted }), true);
  assert.equal(hdf5.sniff({ bytes: new Uint8Array(4096) }), false);
  const r = await hdf5.parse({ name: 'ok.nc', bytes: good });
  assert.equal(r.meta.variable, 'z');
  assert.deepEqual(Array.from(r.models[0].z), [1, 2, 3, 4, 5, 6]);
});
