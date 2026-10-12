import test from 'node:test';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import * as sqlite from '../../src/geometry/adapters/sqlite.js';
import { lonLatToLocal, utmForward } from '../../src/geometry/units.js';
import { importFile } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { close, closeV } from './helpers.js';

// ---------------------------------------------------------------------------------------------------------
// A SQLite database written by hand (https://www.sqlite.org/fileformat2.html): one B-tree per table,
// records in the serial-type format, payloads that do not fit on a page continued on overflow pages.
// ---------------------------------------------------------------------------------------------------------

const enc = new TextEncoder();
const cat = (...parts) => { const out = new Uint8Array(parts.reduce((s, p) => s + p.length, 0)); let o = 0; for (const p of parts) { out.set(p, o); o += p.length; } return out; };

function varint(v) {
  const out = [v & 0x7f];
  v = Math.floor(v / 128);
  while (v > 0) { out.unshift((v & 0x7f) | 0x80); v = Math.floor(v / 128); }
  return Uint8Array.from(out);
}

function record(values) {
  const types = [], bodies = [];
  for (const v of values) {
    if (v === null) { types.push(0); bodies.push(new Uint8Array(0)); } else if (typeof v === 'string') { const b = enc.encode(v); types.push(13 + 2 * b.length); bodies.push(b); } else if (v instanceof Uint8Array) { types.push(12 + 2 * v.length); bodies.push(v); } else if (Number.isInteger(v)) {
      if (v === 0) { types.push(8); bodies.push(new Uint8Array(0)); } else if (v === 1) { types.push(9); bodies.push(new Uint8Array(0)); } else if (v >= -128 && v <= 127) { types.push(1); bodies.push(Uint8Array.from([v & 0xff])); } else if (v >= -32768 && v <= 32767) { const b = new Uint8Array(2); new DataView(b.buffer).setInt16(0, v, false); types.push(2); bodies.push(b); } else if (v >= -8388608 && v <= 8388607) { types.push(3); bodies.push(Uint8Array.from([(v >> 16) & 0xff, (v >> 8) & 0xff, v & 0xff])); } else { const b = new Uint8Array(4); new DataView(b.buffer).setInt32(0, v, false); types.push(4); bodies.push(b); }
    } else { const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, v, false); types.push(7); bodies.push(b); }
  }
  const th = cat(...types.map(varint));
  // the header-length varint counts itself; one byte is enough for the small headers used here
  assert.ok(th.length + 1 < 128);
  return cat(Uint8Array.from([th.length + 1]), th, ...bodies);
}

class DbWriter {
  constructor(pageSize) { this.ps = pageSize; this.pages = [null]; }
  alloc() { this.pages.push(new Uint8Array(this.ps)); return this.pages.length; }
  page(n) { return this.pages[n - 1]; }
  /** Table leaf page holding `rows` = [[rowid, values]]; long payloads spill to overflow pages. */
  leaf(pn, rows) {
    const ps = this.ps, pg = this.page(pn), dv = new DataView(pg.buffer);
    const h = pn === 1 ? 100 : 0;
    const X = ps - 35, M = Math.floor(((ps - 12) * 32) / 255) - 23;
    let end = ps;
    pg[h] = 0x0d;
    dv.setUint16(h + 3, rows.length, false);
    rows.forEach(([rowid, values], i) => {
      const payload = record(values);
      let cell;
      if (payload.length <= X) cell = cat(varint(payload.length), varint(rowid), payload);
      else {
        const K = M + ((payload.length - M) % (ps - 4));
        const local = K <= X ? K : M;
        let rest = payload.subarray(local);
        const first = this.alloc();
        let cur = first;
        for (;;) {
          const o = this.page(cur);
          o.set(rest.subarray(0, ps - 4), 4);
          rest = rest.subarray(Math.min(rest.length, ps - 4));
          if (!rest.length) break;
          const next = this.alloc();
          new DataView(o.buffer).setUint32(0, next, false);
          cur = next;
        }
        const ptr = new Uint8Array(4); new DataView(ptr.buffer).setUint32(0, first, false);
        cell = cat(varint(payload.length), varint(rowid), payload.subarray(0, local), ptr);
      }
      end -= cell.length;
      assert.ok(end >= h + 8 + 2 * rows.length, `page ${pn} overflows`);
      pg.set(cell, end);
      dv.setUint16(h + 8 + 2 * i, end, false);
    });
    dv.setUint16(h + 5, end, false);
  }
  /** Table interior page: `children` = [[leftChildPage, largestRowidInIt]], then the right-most child. */
  interior(pn, children, rightMost) {
    const pg = this.page(pn), dv = new DataView(pg.buffer);
    let end = this.ps;
    pg[0] = 0x05;
    dv.setUint16(3, children.length, false);
    dv.setUint32(8, rightMost, false);
    children.forEach(([child, key], i) => {
      const ptr = new Uint8Array(4); new DataView(ptr.buffer).setUint32(0, child, false);
      const cell = cat(ptr, varint(key));
      end -= cell.length;
      pg.set(cell, end);
      dv.setUint16(12 + 2 * i, end, false);
    });
    dv.setUint16(5, end, false);
  }
  finish({ encoding = 1, readVersion = 1 } = {}) {
    const p1 = this.page(1), dv = new DataView(p1.buffer);
    p1.set(enc.encode('SQLite format 3\0'));
    dv.setUint16(16, this.ps === 65536 ? 1 : this.ps, false);
    p1[18] = readVersion; p1[19] = readVersion; p1[20] = 0; p1[21] = 64; p1[22] = 32; p1[23] = 32;
    dv.setUint32(24, 1, false); dv.setUint32(28, this.pages.length, false);
    dv.setUint32(40, 1, false); dv.setUint32(44, 4, false); dv.setUint32(56, encoding, false);
    dv.setUint32(60, 10200, false); dv.setUint32(68, 0x47504b47, false);
    dv.setUint32(92, 1, false); dv.setUint32(96, 3045001, false);
    return cat(...this.pages);
  }
}

// ---- GeoPackageBinary and WKB builders ----
const f64 = (vals, le) => { const b = new Uint8Array(8 * vals.length); const dv = new DataView(b.buffer); vals.forEach((v, i) => dv.setFloat64(8 * i, v, le)); return b; };
const u32 = (v, le) => { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, v >>> 0, le); return b; };
/** dim: 'XY' | 'Z' | 'M' | 'ZM'; style: 'iso' | 'ewkb'. */
const wkbCode = (base, dim, style) => (style === 'ewkb' ? (base | (dim.includes('Z') ? 0x80000000 : 0) | (dim === 'M' || dim === 'ZM' ? 0x40000000 : 0)) >>> 0 : base + { XY: 0, Z: 1000, M: 2000, ZM: 3000 }[dim]);
const wkbHead = (base, { le = true, dim = 'XY', style = 'iso' } = {}) => cat(Uint8Array.from([le ? 1 : 0]), u32(wkbCode(base, dim, style), le));
const wkbPoint = (c, o = {}) => cat(wkbHead(1, o), f64(c, o.le ?? true));
const wkbLine = (pts, o = {}) => cat(wkbHead(2, o), u32(pts.length, o.le ?? true), f64(pts.flat(), o.le ?? true));
const wkbPolygon = (rings, o = {}) => cat(wkbHead(3, o), u32(rings.length, o.le ?? true), ...rings.map((r) => cat(u32(r.length, o.le ?? true), f64(r.flat(), o.le ?? true))));
const wkbMulti = (base, parts, o = {}) => cat(wkbHead(base, o), u32(parts.length, o.le ?? true), ...parts);
/** GeoPackageBinary header: magic, version 0, flags (byte order, envelope code, empty), srs_id, envelope. */
function gp(srsId, wkb, { le = true, envelope = null, empty = false } = {}) {
  const code = envelope ? { 4: 1, 6: 2, 8: 4 }[envelope.length] : 0;
  const head = new Uint8Array(8);
  head[0] = 0x47; head[1] = 0x50; head[2] = 0; head[3] = (le ? 1 : 0) | (code << 1) | (empty ? 0x10 : 0);
  new DataView(head.buffer).setInt32(4, srsId, le);
  return cat(head, envelope ? f64(envelope, le) : new Uint8Array(0), wkb);
}

const SQL = {
  srs: 'CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL, srs_id INTEGER PRIMARY KEY, organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL, definition TEXT NOT NULL, description TEXT)',
  contents: 'CREATE TABLE gpkg_contents (table_name TEXT NOT NULL, data_type TEXT NOT NULL, identifier TEXT, description TEXT DEFAULT \'\', last_change TEXT, min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER)',
  geom: 'CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL, column_name TEXT NOT NULL, geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL, z TINYINT NOT NULL, m TINYINT NOT NULL)',
  wells: 'CREATE TABLE wells (fid INTEGER PRIMARY KEY NOT NULL, geom POINT, name TEXT, "water depth" REAL, slots INTEGER)',
  route: 'CREATE TABLE route (fid INTEGER PRIMARY KEY, geom LINESTRING, name TEXT, kp_end DOUBLE)',
  zones: 'CREATE TABLE zones (fid INTEGER PRIMARY KEY, name TEXT, geom POLYGON)',
  tiles: 'CREATE TABLE bathy (id INTEGER PRIMARY KEY, zoom_level INTEGER, tile_column INTEGER, tile_row INTEGER, tile_data BLOB)',
};

const WELLS = [[2.90, 53.40, 'A-1', 31.5, 12], [2.95, 53.42, 'A-2', 33.25, 6], [3.00, 53.44, 'B-1', 40, 0], [3.05, 53.46, 'B-2', 41.75, 1], [3.10, 53.48, 'C-1', 300.5, 70000]];
const ROUTE = Array.from({ length: 120 }, (_, i) => [2.9 + 0.003 * i, 53.4 + 0.001 * i + 0.0002 * Math.sin(i), -30 - 0.25 * i]);

/** A five-layer GeoPackage on 2048-byte pages; the `wells` table has an interior root over two leaves. */
function handWrittenGeoPackage(o = {}) {
  const w = new DbWriter(2048);
  w.pages = [];
  const schemaPage = w.alloc(), srs = w.alloc(), contents = w.alloc(), geom = w.alloc(), wellsRoot = w.alloc(), wellsA = w.alloc(), wellsB = w.alloc(), route = w.alloc(), zones = w.alloc(), tiles = w.alloc();
  w.leaf(schemaPage, [
    [1, ['table', 'gpkg_spatial_ref_sys', 'gpkg_spatial_ref_sys', srs, SQL.srs]],
    [2, ['table', 'gpkg_contents', 'gpkg_contents', contents, SQL.contents]],
    [3, ['table', 'gpkg_geometry_columns', 'gpkg_geometry_columns', geom, SQL.geom]],
    [4, ['table', 'wells', 'wells', wellsRoot, SQL.wells]],
    [5, ['table', 'route', 'route', route, SQL.route]],
    [6, ['table', 'zones', 'zones', zones, SQL.zones]],
    [7, ['table', 'bathy', 'bathy', tiles, SQL.tiles]],
  ]);
  // INTEGER PRIMARY KEY columns are stored as NULL: the value is the rowid
  w.leaf(srs, [
    [-0 + 4326, ['WGS 84 geodetic', null, 'EPSG', 4326, 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]', 'longitude/latitude']],
    [32631, ['WGS 84 / UTM zone 31N', null, 'EPSG', 32631, 'undefined', null]],
  ]);
  w.leaf(contents, [
    [1, ['wells', 'features', 'Wells', '', '2024-01-01T00:00:00.000Z', 2.9, 53.4, 3.1, 53.48, 4326]],
    [2, ['route', 'features', 'Export route', '', '2024-01-01T00:00:00.000Z', null, null, null, null, 4326]],
    [3, ['zones', 'features', 'Exclusion zones', '', '2024-01-01T00:00:00.000Z', null, null, null, null, 32631]],
    [4, ['bathy', '2d-gridded-coverage', 'Bathymetry tiles', '', '2024-01-01T00:00:00.000Z', null, null, null, null, 32631]],
  ]);
  w.leaf(geom, [
    [1, ['wells', 'geom', 'POINT', 4326, 0, 0]],
    [2, ['route', 'geom', 'LINESTRING', 4326, 1, 0]],
    [3, ['zones', 'geom', 'POLYGON', 32631, 0, 0]],
  ]);
  const well = (i) => { const [lon, lat, name, wd, slots] = WELLS[i]; return [i + 1, [null, gp(4326, wkbPoint([lon, lat]), { envelope: i % 2 ? [lon, lon, lat, lat] : null }), name, wd, slots]]; };
  w.leaf(wellsA, [well(0), well(1), well(2)]);
  w.leaf(wellsB, [well(3), well(4)]);
  w.interior(wellsRoot, [[wellsA, 3]], wellsB);
  // 120 XYZ points = 2897 bytes of geometry: more than a page, so the record continues on an overflow page
  w.leaf(route, [[1, [null, gp(4326, wkbLine(ROUTE, { dim: 'Z' })), 'Export line', 12.5]]]);
  const outer = [[500000, 5900000], [500400, 5900000], [500400, 5900300], [500000, 5900300], [500000, 5900000]];
  const hole = [[500100, 5900100], [500100, 5900200], [500200, 5900200], [500200, 5900100], [500100, 5900100]];
  w.leaf(zones, [
    [1, [null, 'Anchor exclusion', gp(32631, wkbPolygon([outer, hole], { le: false }), { le: false, envelope: [500000, 500400, 5900000, 5900300] })]],
    [2, [null, 'no geometry', null]],
  ]);
  w.leaf(tiles, [[1, [null, 0, 0, 0, Uint8Array.from([0x89, 0x50, 0x4e, 0x47])]]]);
  return w.finish(o);
}

const python = (() => { try { execFileSync('python3', ['-I', '-c', 'import sqlite3'], { stdio: 'ignore' }); return true; } catch { return false; } })();

test('adapter interface: always available, no runtime to download', () => {
  assert.equal(sqlite.id, 'sqlite');
  assert.equal(sqlite.available(), true);
  assert.equal(sqlite.sizeBytes, 0);
  assert.deepEqual(sqlite.formats, ['geopackage']);
  assert.ok(sqlite.extensions.includes('gpkg'));
  assert.equal(typeof sqlite.library.name, 'string');
});

test('hand-written SQLite file: header, schema, interior + leaf pages, overflow chain, rowid alias', () => {
  const bytes = handWrittenGeoPackage();
  assert.equal(bytes.length, 11 * 2048); // ten table pages + one overflow page
  const db = sqlite.openSqlite(bytes, 'hand.gpkg');
  assert.equal(db.header.pageSize, 2048); assert.equal(db.header.encoding, 'utf-8'); assert.equal(db.header.pages, 11);
  assert.equal(db.header.wal, false); assert.equal(db.header.applicationId, 0x47504b47);
  assert.deepEqual(db.schema.map((s) => s.name), ['gpkg_spatial_ref_sys', 'gpkg_contents', 'gpkg_geometry_columns', 'wells', 'route', 'zones', 'bathy']);
  assert.deepEqual(db.columnsOf('wells').map((c) => [c.name, c.type, c.rowidAlias]), [['fid', 'INTEGER', true], ['geom', 'POINT', false], ['name', 'TEXT', false], ['water depth', 'REAL', false], ['slots', 'INTEGER', false]]);
  assert.equal(db.count('wells'), 5);
  const rows = [];
  db.scan('wells', (rowid, v) => { rows.push([rowid, v[0], v[2], v[3], v[4]]); });
  // serial types 8 (0), 9 (1), 1-byte, 3-byte and float64, in rowid order across the two leaves
  assert.deepEqual(rows, [[1, 1, 'A-1', 31.5, 12], [2, 2, 'A-2', 33.25, 6], [3, 3, 'B-1', 40, 0], [4, 4, 'B-2', 41.75, 1], [5, 5, 'C-1', 300.5, 70000]]);
  let blob = null;
  db.scan('route', (rowid, v) => { blob = v[1]; });
  assert.equal(blob.length, 8 + 9 + 120 * 24);
  const stop = [];
  db.scan('wells', (rowid) => { stop.push(rowid); return rowid < 2; });
  assert.deepEqual(stop, [1, 2]);
  assert.throws(() => db.scan('nothing', () => {}), /no table "nothing"/);
});

test('hand-written file is a database SQLite itself accepts', { skip: python ? false : 'python3 with the sqlite3 module is not installed' }, () => {
  const dir = mkdtempSync(join(tmpdir(), 'gpkg-'));
  try {
    const f = join(dir, 'hand.gpkg');
    writeFileSync(f, handWrittenGeoPackage());
    const out = execFileSync('python3', ['-I', '-c', 'import sqlite3,sys,json\nc=sqlite3.connect(sys.argv[1])\nprint(json.dumps([c.execute("PRAGMA integrity_check").fetchall(), c.execute("SELECT fid,name,slots FROM wells ORDER BY fid").fetchall(), c.execute("SELECT length(geom) FROM route").fetchall()]))', f], { encoding: 'utf8' });
    assert.deepEqual(JSON.parse(out), [[['ok']], [[1, 'A-1', 12], [2, 'A-2', 6], [3, 'B-1', 0], [4, 'B-2', 1], [5, 'C-1', 70000]], [[2897]]]);
  } finally { rmSync(dir, { recursive: true, force: true }); }
});

test('GeoPackage layers → points, polyline and section; lon/lat in local metres about one origin; tiles reported', async () => {
  const bytes = handWrittenGeoPackage();
  const r = await sqlite.parse({ name: 'field.gpkg', bytes }, { origin: { lon: 3, lat: 53.44 } });
  assert.deepEqual(r.meta.counts.byKind, { points: 1, polyline: 1, section: 1 });
  assert.equal(r.meta.format, 'geopackage');
  assert.deepEqual(r.meta.sqlite, { pageSize: 2048, encoding: 'utf-8', pages: 11, journal: 'rollback' });
  const [pts, line, zone] = r.models;
  // points with attribute columns
  assert.equal(pts.kind, 'points'); assert.equal(pts.name, 'wells'); assert.equal(pts.layer, 'wells');
  assert.equal(pts.positions.length, 15);
  WELLS.forEach(([lon, lat], i) => { const e = lonLatToLocal(lon, lat, 3, 53.44); close(pts.positions[3 * i], e[0], 0.01); close(pts.positions[3 * i + 1], e[1], 0.01); });
  close(pts.positions[0], -6648.6, 5); // 0.1° of longitude at 53.44° N
  assert.deepEqual(Array.from(pts.attrs['water depth']), [31.5, 33.25, 40, 41.75, 300.5]);
  assert.deepEqual(Array.from(pts.attrs.slots), [12, 6, 0, 1, 70000]);
  assert.deepEqual(Array.from(pts.attrs.fid), [1, 2, 3, 4, 5]);
  assert.deepEqual(pts.attrs.names, ['A-1', 'A-2', 'B-1', 'B-2', 'C-1']);
  // polyline with z, read through the overflow chain
  assert.equal(line.kind, 'polyline'); assert.equal(line.name, 'Export line');
  assert.equal(line.points.length, 360);
  ROUTE.forEach(([lon, lat, z], i) => { const e = lonLatToLocal(lon, lat, 3, 53.44); closeV(Array.from(line.points.subarray(3 * i, 3 * i + 3)), [e[0], e[1], z], 1e-6); });
  assert.deepEqual(line.properties, { fid: 1, name: 'Export line', kp_end: 12.5 });
  // polygon with a hole, big-endian, projected CRS: coordinates unchanged
  assert.equal(zone.kind, 'section'); assert.equal(zone.name, 'Anchor exclusion');
  assert.equal(zone.loops.length, 2);
  assert.deepEqual(Array.from(zone.loops[0]), [500000, 5900000, 500400, 5900000, 500400, 5900300, 500000, 5900300]);
  assert.ok(inspect(zone).ok);
  // CRS per layer
  const L = Object.fromEntries(r.meta.layers.map((l) => [l.name, l]));
  assert.equal(L.wells.crs.epsg, 4326); assert.equal(L.wells.crs.kind, 'geographic'); assert.equal(L.wells.crs.projection, 'local-tangent-plane');
  assert.deepEqual(L.wells.crs.origin, { lon: 3, lat: 53.44 });
  assert.equal(L.zones.crs.epsg, 32631); assert.equal(L.zones.crs.kind, 'projected'); assert.equal(L.zones.srsId, 32631);
  assert.equal(L.wells.rows, 5); assert.equal(L.wells.geometryType, 'POINT'); assert.deepEqual(L.wells.geometryTypes, { Point: 5 });
  assert.equal(L.route.models, 1);
  assert.equal(r.meta.crs.epsg, 4326);
  // what was not read is said
  assert.deepEqual(r.meta.tiles, [{ name: 'bathy', identifier: 'Bathymetry tiles', dataType: '2d-gridded-coverage', tiles: 1, srsId: 32631 }]);
  assert.ok(r.warnings.some((w) => /"bathy" is a tiled gridded coverage \(elevation tiles\) with 1 tile\(s\): raster tiles are not decoded.*gdal_translate/.test(w)));
  assert.ok(r.warnings.some((w) => /Layer "zones": 1 row\(s\) have no geometry/.test(w)));
  assert.ok(r.warnings.some((w) => /different coordinate systems/.test(w)));
});

test('without an origin the lon/lat layers share the centre of all of them; layer selection; UTM projection', async () => {
  const bytes = handWrittenGeoPackage();
  const r = await sqlite.parse({ name: 'field.gpkg', bytes });
  // wells span 2.9–3.1° E, the route 2.9–3.257° E: the union decides
  const o = r.meta.layers[0].crs.origin;
  close(o.lon, 3.0785, 1e-9);
  assert.deepEqual(r.meta.layers[1].crs.origin, o);
  assert.equal(r.meta.layers[0].crs.originBasis, 'data centre of all layers');
  const e = lonLatToLocal(2.9, 53.4, o.lon, o.lat);
  // the cloud is more than 10 km from the origin, so its single-precision positions are stored about model.origin
  close(r.models[0].positions[0] + r.models[0].origin[0], e[0], 0.01); close(r.models[1].points[0], e[0], 1e-6); close(r.models[1].points[1], e[1], 1e-6);

  const one = await sqlite.parse({ name: 'field.gpkg', bytes }, { layers: ['WELLS'], projection: 'utm' });
  assert.equal(one.models.length, 1);
  const u = utmForward(53.4, 2.9, 31);
  assert.equal(one.models[0].origin.length, 3); // UTM northings need the integer-metre origin for single precision
  close(one.models[0].positions[0] + one.models[0].origin[0], u.easting, 0.05);
  close(one.models[0].positions[1] + one.models[0].origin[1], u.northing, 0.5);
  assert.equal(one.meta.crs.projection, 'utm'); assert.equal(one.meta.crs.zone, 31);
  const none = await sqlite.parse({ name: 'field.gpkg', bytes }, { layers: 'pipes' });
  assert.equal(none.models.length, 0);
  assert.ok(none.warnings.some((w) => /Layer "pipes" was asked for but is not a feature table.*wells, route, zones/.test(w)));
  // through the registry: detected from the SQLite signature
  const viaRegistry = await importFile({ name: 'field.gpkg', bytes });
  assert.equal(viaRegistry.format, 'geopackage'); assert.equal(viaRegistry.meta.via, 'adapter:sqlite'); assert.equal(viaRegistry.models.length, 3);
});

test('GeoPackageBinary header and WKB: byte orders, envelopes, ISO and EWKB codes, Z / M / ZM, Multi*, collections', () => {
  for (const le of [true, false]) {
    for (const [envelope, off] of [[null, 8], [[1, 2, 3, 4], 40], [[1, 2, 3, 4, 5, 6], 56], [[1, 2, 3, 4, 5, 6, 7, 8], 72]]) {
      const h = sqlite.gpkgHeader(gp(-1, wkbPoint([7, 8], { le }), { le, envelope }));
      assert.equal(h.srsId, -1); assert.equal(h.wkbOffset, off); assert.deepEqual(h.envelope, envelope); assert.equal(h.empty, false);
    }
    for (const style of ['iso', 'ewkb']) {
      for (const [dim, c, z] of [['XY', [1, 2], 0], ['Z', [1, 2, 3], 3], ['M', [1, 2, 9], 0], ['ZM', [1, 2, 3, 9], 3]]) {
        const p = sqlite.parseWkb(wkbPoint(c, { le, dim, style }));
        assert.equal(p.type, 'Point'); assert.deepEqual(Array.from(p.coords), [1, 2, z]);
        assert.equal(p.hasZ, dim.includes('Z')); assert.equal(p.hasM, dim.includes('M'));
        assert.equal(p.end, 5 + 8 * c.length);
        const l = sqlite.parseWkb(wkbLine([c, c.map((v) => v + 10)], { le, dim, style }));
        assert.deepEqual(Array.from(l.coords), [1, 2, z, 11, 12, z ? z + 10 : 0]);
      }
    }
    const o = { le };
    const mp = sqlite.parseWkb(wkbMulti(4, [wkbPoint([1, 2], o), wkbPoint([3, 4], { le: !le })], o));
    assert.equal(mp.type, 'MultiPoint'); assert.deepEqual(mp.parts.map((q) => Array.from(q.coords)), [[1, 2, 0], [3, 4, 0]]);
    const ml = sqlite.parseWkb(wkbMulti(5, [wkbLine([[0, 0], [1, 1]], o), wkbLine([[2, 2], [3, 3], [4, 4]], o)], o));
    assert.deepEqual(ml.parts.map((q) => q.coords.length / 3), [2, 3]);
    const ring = [[0, 0], [4, 0], [4, 4], [0, 0]];
    const mpoly = sqlite.parseWkb(wkbMulti(6, [wkbPolygon([ring], o), wkbPolygon([ring, ring], o)], o));
    assert.deepEqual(mpoly.parts.map((q) => q.rings.length), [1, 2]);
    const gc = sqlite.parseWkb(wkbMulti(7, [wkbPoint([1, 2], o), wkbMulti(7, [wkbLine([[0, 0], [1, 1]], o)], o), wkbPolygon([ring], o)], o));
    assert.deepEqual(gc.parts.map((q) => q.type), ['Point', 'GeometryCollection', 'Polygon']);
    assert.equal(gc.parts[1].parts[0].type, 'LineString');
  }
  // EWKB with an embedded SRID
  const withSrid = cat(Uint8Array.from([1]), u32(0x20000001, true), u32(4326, true), f64([5, 6], true));
  const s = sqlite.parseWkb(withSrid);
  assert.equal(s.srid, 4326); assert.deepEqual(Array.from(s.coords), [5, 6, 0]);
  assert.throws(() => sqlite.parseWkb(wkbLine([[0, 0], [1, 1]]).subarray(0, 20)), /truncated/);
  assert.throws(() => sqlite.parseWkb(cat(Uint8Array.from([1]), u32(99, true))), /unknown WKB geometry type code 99/);
  assert.throws(() => sqlite.gpkgHeader(Uint8Array.from([1, 2, 3, 4, 5, 6, 7, 8])), /does not start with "GP"/);
  assert.equal(sqlite.gpkgHeader(gp(0, new Uint8Array(0), { empty: true })).empty, true);
});

test('CREATE TABLE column parser: quoting, constraints, table-level primary key, comments', () => {
  const c = sqlite.parseColumns('CREATE TABLE "my table" (\n  [id] INTEGER NOT NULL, -- key\n  "the geom" MULTIPOLYGON,\n  `na,me` VARCHAR(40) DEFAULT \'a,b\' CHECK (length(`na,me`) < 41),\n  depth DOUBLE PRECISION /* m */,\n  bare,\n  CONSTRAINT pk PRIMARY KEY (id),\n  UNIQUE (depth), FOREIGN KEY (bare) REFERENCES t(x)\n)');
  assert.deepEqual(c.map((q) => [q.name, q.type, q.rowidAlias]), [['id', 'INTEGER', true], ['the geom', 'MULTIPOLYGON', false], ['na,me', 'VARCHAR', false], ['depth', 'DOUBLE PRECISION', false], ['bare', '', false]]);
  assert.equal(sqlite.parseColumns('CREATE TABLE t (a INT PRIMARY KEY, b TEXT)')[0].rowidAlias, false); // INT is not INTEGER
  assert.equal(sqlite.parseColumns('CREATE TABLE t (a INTEGER PRIMARY KEY DESC, b TEXT)')[0].rowidAlias, false);
});

test('files that are not readable GeoPackages fail with what to do; WAL mode and odd lengths are flagged', async () => {
  await assert.rejects(() => sqlite.parse({ name: 'a.gpkg', bytes: enc.encode('hello world, this is not a database') }), /a\.gpkg: not a SQLite database/);
  await assert.rejects(() => sqlite.parse({ name: 'a.gpkg', text: 'x' }), /binary SQLite file/);
  await assert.rejects(() => sqlite.parse({ name: 'a.gpkg', bytes: handWrittenGeoPackage().subarray(0, 300) }), /truncated/);
  // a plain SQLite database
  const w = new DbWriter(512); w.pages = [];
  const p1 = w.alloc(), t = w.alloc();
  w.leaf(p1, [[1, ['table', 'notes', 'notes', t, 'CREATE TABLE notes (id INTEGER PRIMARY KEY, body TEXT)']]]);
  w.leaf(t, [[1, [null, 'x']]]);
  await assert.rejects(() => sqlite.parse({ name: 'notes.sqlite', bytes: w.finish() }), /not a GeoPackage \(it has no gpkg_contents table\)\. Tables: notes\./);
  // cut off behind the schema page: a table root that is not in the file
  await assert.rejects(() => sqlite.parse({ name: 'cut.gpkg', bytes: handWrittenGeoPackage().subarray(0, 4 * 2048) }), /refers to page \d+, but the file holds 4 page\(s\).*truncated or damaged/);
  const wal = await sqlite.parse({ name: 'wal.gpkg', bytes: cat(handWrittenGeoPackage({ readVersion: 2 }), new Uint8Array(100)) });
  assert.equal(wal.meta.sqlite.journal, 'wal');
  assert.ok(wal.warnings.some((s) => /WAL journal mode.*wal_checkpoint/.test(s)));
  assert.ok(wal.warnings.some((s) => /not a whole number of pages/.test(s)));
  assert.equal(wal.models.length, 3);
});

// ---------------------------------------------------------------------------------------------------------
// A larger GeoPackage written by SQLite itself: multi-level B-trees, overflow chains, indexes, a WITHOUT ROWID
// table, UTF-8 text, columns added after rows were written.
// ---------------------------------------------------------------------------------------------------------

const PY = String.raw`
import sqlite3, struct, sys, math
c = sqlite3.connect(sys.argv[1])
c.execute('PRAGMA page_size=4096')
c.execute('PRAGMA application_id=1196444487')
c.execute('PRAGMA user_version=10300')
if sys.argv[2] == 'wal': c.execute('PRAGMA journal_mode=WAL')
c.executescript('''
CREATE TABLE gpkg_spatial_ref_sys (srs_name TEXT NOT NULL, srs_id INTEGER PRIMARY KEY, organization TEXT NOT NULL, organization_coordsys_id INTEGER NOT NULL, definition TEXT NOT NULL, description TEXT);
CREATE TABLE gpkg_contents (table_name TEXT NOT NULL PRIMARY KEY, data_type TEXT NOT NULL, identifier TEXT UNIQUE, description TEXT DEFAULT '', last_change DATETIME NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')), min_x DOUBLE, min_y DOUBLE, max_x DOUBLE, max_y DOUBLE, srs_id INTEGER, CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id));
CREATE TABLE gpkg_geometry_columns (table_name TEXT NOT NULL, column_name TEXT NOT NULL, geometry_type_name TEXT NOT NULL, srs_id INTEGER NOT NULL, z TINYINT NOT NULL, m TINYINT NOT NULL, CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name));
CREATE TABLE soundings (fid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, geom POINT, depth REAL, line TEXT);
CREATE TABLE cables (fid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, name TEXT, geom MULTILINESTRING);
CREATE TABLE blocks (fid INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, geom MULTIPOLYGON, licence TEXT);
CREATE TABLE lookup (code TEXT PRIMARY KEY, meaning TEXT) WITHOUT ROWID;
CREATE TABLE tiles_dem (id INTEGER PRIMARY KEY AUTOINCREMENT, zoom_level INTEGER NOT NULL, tile_column INTEGER NOT NULL, tile_row INTEGER NOT NULL, tile_data BLOB NOT NULL);
CREATE INDEX soundings_line ON soundings(line);
''')
c.execute("INSERT INTO gpkg_spatial_ref_sys VALUES ('WGS 84 / UTM zone 31N', 32631, 'EPSG', 32631, 'PROJCS[\"WGS 84 / UTM zone 31N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",3],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1]]', NULL)")
c.execute("INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined Cartesian SRS', -1, 'NONE', -1, 'undefined', NULL)")
for t, k, g in (('soundings', 'features', 'POINT'), ('cables', 'features', 'MULTILINESTRING'), ('blocks', 'features', 'MULTIPOLYGON')):
    c.execute('INSERT INTO gpkg_contents (table_name, data_type, identifier, srs_id) VALUES (?,?,?,32631)', (t, k, t))
    c.execute('INSERT INTO gpkg_geometry_columns VALUES (?,?,?,32631,1,0)', (t, 'geom', g))
c.execute("INSERT INTO gpkg_contents (table_name, data_type, identifier, srs_id) VALUES ('tiles_dem','tiles','dem',32631)")
c.execute("INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('lookup','attributes','lookup')")
def gp(wkb, srs=32631): return b'GP' + bytes([0, 1]) + struct.pack('<i', srs) + wkb
def pt(x, y, z): return struct.pack('<BIddd', 1, 1001, x, y, z)
def ls(p): return struct.pack('<BII', 1, 1002, len(p)) + b''.join(struct.pack('<ddd', *q) for q in p)
def poly(rings): return struct.pack('>BII', 0, 3, len(rings)) + b''.join(struct.pack('>I', len(r)) + b''.join(struct.pack('>dd', *q) for q in r) for r in rings)
N = 20000
c.executemany('INSERT INTO soundings (geom, depth, line) VALUES (?,?,?)', ((gp(pt(500000 + (i % 200) * 5.0, 5900000 + (i // 200) * 5.0, -40 - 0.001 * i)), 40 + 0.001 * i, 'L%03d' % (i // 200)) for i in range(N)))
for k in range(3):
    parts = [[(500000 + 10.0 * i, 5900000 + 100.0 * k + 20.0 * j + math.sin(i / 50.0), -41.0 - j) for i in range(700 + 300 * j)] for j in range(2)]
    c.execute('INSERT INTO cables (name, geom) VALUES (?,?)', ('Kabel Øst %d – 北' % (k + 1), gp(struct.pack('<BII', 1, 1005, 2) + b''.join(ls(p) for p in parts))))
sq = lambda x, y, s: [(x, y), (x + s, y), (x + s, y + s), (x, y + s), (x, y)]
c.execute('INSERT INTO blocks (geom, licence) VALUES (?,?)', (gp(struct.pack('<BII', 1, 6, 2) + poly([sq(500000, 5900000, 1000), sq(500200, 5900200, 100)[::-1]]) + poly([sq(502000, 5900000, 500)])), 'PL-001'))
c.execute('INSERT INTO blocks (geom, licence) VALUES (?,?)', (None, 'PL-002'))
c.execute('ALTER TABLE blocks ADD COLUMN operator TEXT')
c.execute("INSERT INTO blocks (geom, licence, operator) VALUES (?,?,?)", (gp(struct.pack('<BII', 1, 6, 1) + poly([sq(504000, 5900000, 250)])), 'PL-003', 'Nordsee AS'))
c.execute("INSERT INTO lookup VALUES ('a','alpha')")
c.execute('INSERT INTO tiles_dem (zoom_level, tile_column, tile_row, tile_data) VALUES (0,0,0,?)', (b'\x89PNG' + bytes(9000),))
c.execute('DELETE FROM soundings WHERE fid % 1000 = 0')
c.commit()
c.close()
`;

function sqliteBuilt(mode) {
  const dir = mkdtempSync(join(tmpdir(), 'gpkg-'));
  try {
    const f = join(dir, 'survey.gpkg');
    execFileSync('python3', ['-I', '-c', PY, f, mode], { stdio: ['ignore', 'ignore', 'pipe'] });
    return new Uint8Array(readFileSync(f));
  } finally { rmSync(dir, { recursive: true, force: true }); }
}

test('GeoPackage written by SQLite: 19 980 soundings over a multi-level B-tree, long multi-line blobs, late columns', { skip: python ? false : 'python3 with the sqlite3 module is not installed, so no SQLite-written fixture can be generated' }, async () => {
  const bytes = sqliteBuilt('rollback');
  const db = sqlite.openSqlite(bytes, 'survey.gpkg');
  assert.equal(db.header.pageSize, 4096);
  assert.equal(db.count('soundings'), 19980);
  assert.equal(db.tables.get('lookup').withoutRowid, true);
  assert.throws(() => db.scan('lookup', () => {}), /WITHOUT ROWID/);
  const r = await sqlite.parse({ name: 'survey.gpkg', bytes });
  assert.equal(r.meta.application, 'GeoPackage (user_version 10300)');
  const L = Object.fromEntries(r.meta.layers.map((l) => [l.name, l]));
  assert.equal(L.soundings.rows, 19980); assert.equal(L.cables.rows, 3); assert.equal(L.blocks.rows, 3);
  assert.equal(L.soundings.crs.epsg, 32631); assert.equal(L.soundings.crs.kind, 'projected');
  const pts = r.models.find((m) => m.kind === 'points');
  assert.equal(pts.positions.length / 3, 19980);
  // rows deleted with "fid % 1000 = 0" are gone; fids keep their values and order
  assert.equal(pts.attrs.fid[0], 1); assert.equal(pts.attrs.fid[998], 999); assert.equal(pts.attrs.fid[999], 1001); assert.equal(pts.attrs.fid[19979], 19999);
  const o = pts.origin;
  for (const k of [0, 5000, 19979]) {
    const i = pts.attrs.fid[k] - 1;
    closeV([pts.positions[3 * k] + o[0], pts.positions[3 * k + 1] + o[1], pts.positions[3 * k + 2] + o[2]], [500000 + (i % 200) * 5, 5900000 + Math.floor(i / 200) * 5, -40 - 0.001 * i], 1e-2);
    close(pts.attrs.depth[k], 40 + 0.001 * i, 1e-12);
    assert.equal(pts.attrs.line[k], `L${String(Math.floor(i / 200)).padStart(3, '0')}`);
  }
  // six polylines (three cables × two parts), each long enough to need overflow pages
  const lines = r.models.filter((m) => m.kind === 'polyline');
  assert.deepEqual(lines.map((m) => m.points.length / 3), [700, 1000, 700, 1000, 700, 1000]);
  assert.equal(lines[2].name, 'Kabel Øst 2 – 北');
  closeV(Array.from(lines[3].points.subarray(3 * 999, 3 * 1000)), [500000 + 9990, 5900000 + 100 + 20 + Math.sin(999 / 50), -42], 1e-9);
  // polygons: big-endian WKB inside a little-endian collection; a hole; a column added later
  const secs = r.models.filter((m) => m.kind === 'section');
  assert.equal(secs.length, 3);
  assert.deepEqual(secs.map((s) => s.loops.length), [2, 1, 1]);
  assert.deepEqual(secs[0].properties, { fid: 1, licence: 'PL-001' });
  assert.deepEqual(secs[2].properties, { fid: 3, licence: 'PL-003', operator: 'Nordsee AS' });
  assert.deepEqual(Array.from(secs[2].loops[0]), [504000, 5900000, 504250, 5900000, 504250, 5900250, 504000, 5900250]);
  assert.deepEqual(r.meta.tiles, [{ name: 'tiles_dem', identifier: 'dem', dataType: 'tiles', tiles: 1, srsId: 32631 }]);
  assert.ok(r.warnings.some((w) => /"tiles_dem" is a raster tile pyramid with 1 tile\(s\): raster tiles are not decoded/.test(w)));
  assert.ok(r.warnings.some((w) => /Table "lookup" is a WITHOUT ROWID table; such tables are not read/.test(w)));
  assert.ok(r.warnings.some((w) => /Layer "blocks": 1 row\(s\) have no geometry/.test(w)));
  assert.deepEqual(r.meta.otherContents, [{ name: 'lookup', dataType: 'attributes' }]);
  assert.deepEqual(r.meta.bbox.min.map(Math.round), [500000, 5899999, -60]); // the first cable meanders 1 m south of the grid
  const few = await sqlite.parse({ name: 'survey.gpkg', bytes }, { layers: ['soundings'], maxFeatures: 250 });
  assert.equal(few.models[0].positions.length / 3, 250);
  assert.ok(few.warnings.some((w) => /stopped after maxFeatures = 250/.test(w)));
});

test('GeoPackage in WAL mode written by SQLite is read as of its checkpoint, with a warning', { skip: python ? false : 'python3 with the sqlite3 module is not installed, so no SQLite-written fixture can be generated' }, async () => {
  const r = await sqlite.parse({ name: 'survey.gpkg', bytes: sqliteBuilt('wal') });
  assert.equal(r.meta.sqlite.journal, 'wal');
  assert.ok(r.warnings.some((w) => /WAL journal mode/.test(w)));
  assert.equal(r.models.find((m) => m.kind === 'points').positions.length / 3, 19980);
});
