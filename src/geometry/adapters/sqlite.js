// OGC GeoPackage (vector feature tables) read straight from the SQLite file, in plain JavaScript.
//
// No SQL engine is needed to list the rows of a table: this module reads the SQLite 3 file format itself
// (https://www.sqlite.org/fileformat2.html) — the 100-byte database header, table B-tree interior and leaf
// pages, varints, the record format with its serial types and overflow page chains — and finds the tables
// through sqlite_schema on page 1. It is read-only and deliberately small:
//   - index B-trees are never needed for a full table scan and are ignored;
//   - WITHOUT ROWID tables are stored as index B-trees and are reported, not read;
//   - a database in WAL journal mode is read as of its last checkpoint (changes still in a side-car "-wal"
//     file are not part of the main file) and that is said in a warning.
//
// On top of it, GeoPackage 1.0–1.4 (OGC 12-128): gpkg_contents, gpkg_geometry_columns and gpkg_spatial_ref_sys
// are read, each feature table is scanned and its GeoPackageBinary geometries (header "GP" + WKB) become
//   Point / MultiPoint            → one `points` model per layer, attribute columns as per-point attributes
//   LineString / MultiLineString  → one `polyline` per line, attribute columns in `properties`
//   Polygon / MultiPolygon        → one `section` per flat polygon (outer loop, then holes); polygons with
//                                   varying z become closed polylines
// exactly as for GeoJSON and Shapefile (importers/vector.js). XY, XYZ, XYM and XYZM coordinates, both byte
// orders and both ISO (1001, 2001, 3001 …) and EWKB (high-bit flags) type codes are understood; M is dropped.
// Longitude/latitude layers are projected to metres with units.geoProjector — one common origin for all
// layers of the file so that they stay aligned.
//
// Tile pyramids (imagery and tiled gridded coverage, i.e. elevation stored as PNG / TIFF tiles) are listed but
// not decoded.

import { inputBytes, baseName, finish } from '../importers/common.js';
import { Features, crsFromName } from '../importers/vector.js';
import { parsePrj } from '../units.js';

export const id = 'sqlite';
export const title = 'GeoPackage vector layers (SQLite container)';
export const formats = ['geopackage'];
export const extensions = ['gpkg', 'sqlite'];
export const library = { name: 'built-in SQLite file-format reader', version: '1', licence: 'part of this application' };
/** Download size of a runtime [bytes]: none, the reader is this module. */
export const sizeBytes = 0;
/** What to convert the parts that are not read to. */
export const convertTo = 'GeoTIFF for raster tiles (gdal_translate in.gpkg out.tif); GeoJSON or Shapefile for anything else (ogr2ogr -f GeoJSON out.geojson in.gpkg layer)';

/** Always true: nothing has to be downloaded or compiled. */
export function available() { return true; }

// ---------------------------------------------------------------------------------------------------------
// SQLite file format
// ---------------------------------------------------------------------------------------------------------

const MAGIC = 'SQLite format 3\u0000';
const ENCODINGS = { 1: 'utf-8', 2: 'utf-16le', 3: 'utf-16be' };

/** True when the bytes start with the SQLite 3 signature. */
export function sniff(input) {
  const b = inputBytes(input);
  if (b.length < 16) return false;
  for (let i = 0; i < 16; i++) if (b[i] !== MAGIC.charCodeAt(i)) return false;
  return true;
}

/**
 * Open a SQLite 3 database image for reading.
 * @param {Uint8Array} bytes the whole database file (not copied)
 * @param {string} [fileName] used in error messages
 * @returns {{ header: object, schema: object[], tables: Map<string, object>, scan: Function, count: Function,
 *   columnsOf: Function }}
 *   `schema` rows are `{ type, name, tblName, rootPage, sql }`; `tables` maps lower-case table names to
 *   `{ name, rootPage, sql, columns: [{ name, type, rowidAlias }], withoutRowid, virtual }`;
 *   `scan(table, (rowid, values) => void|false)` walks a table in rowid order (`values` are numbers, strings,
 *   Uint8Array views for blobs, or null; returning false stops); `count(table)` counts rows without decoding.
 */
export function openSqlite(bytes, fileName = 'SQLite') {
  if (!sniff({ bytes })) throw new Error(`${fileName}: not a SQLite database (the file does not start with "SQLite format 3").`);
  if (bytes.length < 512) throw new Error(`${fileName}: the SQLite file is truncated (${bytes.length} bytes; the smallest database is one 512-byte page).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let pageSize = dv.getUint16(16, false);
  if (pageSize === 1) pageSize = 65536;
  if (pageSize < 512 || (pageSize & (pageSize - 1)) !== 0) throw new Error(`${fileName}: invalid SQLite page size ${pageSize} in the header (expected a power of two from 512 to 65536).`);
  const reserved = bytes[20];
  const usable = pageSize - reserved;
  if (usable < 480) throw new Error(`${fileName}: invalid SQLite header (usable page size ${usable}).`);
  const encCode = dv.getUint32(56, false) || 1;
  const encoding = ENCODINGS[encCode];
  if (!encoding) throw new Error(`${fileName}: unknown SQLite text encoding code ${encCode} (expected 1, 2 or 3).`);
  const pagesInFile = Math.floor(bytes.length / pageSize);
  const header = {
    pageSize, reserved, encoding, pages: pagesInFile, declaredPages: dv.getUint32(28, false),
    writeVersion: bytes[18], readVersion: bytes[19], wal: bytes[18] === 2 || bytes[19] === 2,
    schemaFormat: dv.getUint32(44, false), applicationId: dv.getUint32(68, false), userVersion: dv.getUint32(60, false),
    sqliteVersion: dv.getUint32(96, false),
  };
  if (bytes.length % pageSize) header.partialLastPage = true;
  const text = new TextDecoder(encoding);

  const pageOffset = (n) => {
    if (!(n >= 1 && n <= pagesInFile)) throw new Error(`${fileName}: the database refers to page ${n}, but the file holds ${pagesInFile} page(s) of ${pageSize} bytes — it is truncated or damaged.`);
    return (n - 1) * pageSize;
  };

  // varint: big-endian, 7 bits per byte, the ninth byte contributes all 8 bits. Returns the value; the
  // position after it is left in `vp`.
  let vp = 0;
  const varint = (p) => {
    let v = 0;
    for (let i = 0; i < 8; i++) {
      const b = bytes[p++];
      v = v * 128 + (b & 0x7f);
      if (b < 0x80) { vp = p; return v; }
    }
    // nine bytes: a full 64-bit value (negative rowids land here)
    let big = 0n;
    for (let i = 0; i < 8; i++) big = (big << 7n) | BigInt(bytes[p - 8 + i] & 0x7f);
    big = BigInt.asIntN(64, (big << 8n) | BigInt(bytes[p]));
    vp = p + 1;
    return Number(big);
  };

  // Payload of a table-leaf cell: local part on the page, the rest along a chain of overflow pages.
  const maxLocal = usable - 35;
  const minLocal = Math.floor(((usable - 12) * 32) / 255) - 23;
  const payload = (p, size) => {
    if (size <= maxLocal) return bytes.subarray(p, p + size);
    const k = minLocal + ((size - minLocal) % (usable - 4));
    const local = k <= maxLocal ? k : minLocal;
    const out = new Uint8Array(size);
    out.set(bytes.subarray(p, p + local));
    let filled = local;
    let next = dv.getUint32(p + local, false);
    let hops = 0;
    while (filled < size) {
      if (!next || ++hops > pagesInFile) throw new Error(`${fileName}: an overflow page chain ends early or loops (${filled} of ${size} payload bytes found) — the database is damaged.`);
      const o = pageOffset(next);
      const n = Math.min(usable - 4, size - filled);
      out.set(bytes.subarray(o + 4, o + 4 + n), filled);
      filled += n;
      next = dv.getUint32(o, false);
    }
    return out;
  };

  // Record: header length varint, one serial-type varint per column, then the column values.
  const record = (rec) => {
    const rdv = new DataView(rec.buffer, rec.byteOffset, rec.byteLength);
    // varints are decoded on `bytes`; records in overflow copies need their own reader
    let p = 0;
    const rv = () => {
      let v = 0;
      for (let i = 0; i < 8; i++) { const b = rec[p++]; v = v * 128 + (b & 0x7f); if (b < 0x80) return v; }
      return v * 256 + rec[p++];
    };
    const hdr = rv();
    const types = [];
    while (p < hdr && p < rec.length) types.push(rv());
    let o = hdr;
    const out = new Array(types.length);
    for (let i = 0; i < types.length; i++) {
      const t = types[i];
      switch (t) {
        case 0: out[i] = null; break;
        case 1: out[i] = rdv.getInt8(o); o += 1; break;
        case 2: out[i] = rdv.getInt16(o, false); o += 2; break;
        case 3: out[i] = (rdv.getInt8(o) * 65536) + rdv.getUint16(o + 1, false); o += 3; break;
        case 4: out[i] = rdv.getInt32(o, false); o += 4; break;
        case 5: out[i] = rdv.getInt16(o, false) * 4294967296 + rdv.getUint32(o + 2, false); o += 6; break;
        case 6: { const b = rdv.getBigInt64(o, false); out[i] = Number(b); o += 8; break; }
        case 7: out[i] = rdv.getFloat64(o, false); o += 8; break;
        case 8: out[i] = 0; break;
        case 9: out[i] = 1; break;
        case 10: case 11: throw new Error(`${fileName}: a record uses reserved serial type ${t} — the database is damaged or written by an unsupported SQLite build.`);
        default: {
          const len = (t - 12 - (t & 1)) / 2;
          if (o + len > rec.length) throw new Error(`${fileName}: a record is longer than its payload — the database is damaged.`);
          out[i] = t & 1 ? text.decode(rec.subarray(o, o + len)) : rec.subarray(o, o + len);
          o += len;
        }
      }
    }
    return out;
  };

  /** Depth-first walk of a table B-tree. `onLeafCell(cellOffset)` returning false stops the walk. */
  const walk = (root, onLeaf) => {
    const stack = [root];
    let visited = 0;
    while (stack.length) {
      const pn = stack.pop();
      if (++visited > pagesInFile + 1) throw new Error(`${fileName}: a table B-tree loops back on itself — the database is damaged.`);
      const base = pageOffset(pn);
      const h = base + (pn === 1 ? 100 : 0);
      const type = bytes[h];
      const nCells = dv.getUint16(h + 3, false);
      if (type === 0x0d) {
        if (onLeaf(base, h + 8, nCells) === false) return;
      } else if (type === 0x05) {
        // children are pushed right to left so that they pop in key order
        stack.push(dv.getUint32(h + 8, false));
        for (let i = nCells - 1; i >= 0; i--) stack.push(dv.getUint32(base + dv.getUint16(h + 12 + 2 * i, false), false));
      } else if (type === 0x02 || type === 0x0a) {
        const e = new Error(`${fileName}: page ${pn} is an index B-tree page, not a table page.`);
        e.code = 'INDEX_BTREE';
        throw e;
      } else throw new Error(`${fileName}: page ${pn} has unknown B-tree page type 0x${type.toString(16)} — the database is damaged.`);
    }
  };

  const scanRoot = (root, cb) => {
    walk(root, (base, ptrs, nCells) => {
      for (let i = 0; i < nCells; i++) {
        const c = base + dv.getUint16(ptrs + 2 * i, false);
        const size = varint(c);
        const rowid = varint(vp);
        if (cb(rowid, record(payload(vp, size))) === false) return false;
      }
      return true;
    });
  };
  const countRoot = (root) => { let n = 0; walk(root, (base, ptrs, nCells) => { n += nCells; }); return n; };

  // ---- schema ----
  const schema = [];
  scanRoot(1, (rowid, v) => { schema.push({ type: String(v[0]), name: String(v[1]), tblName: String(v[2]), rootPage: Number(v[3]) || 0, sql: v[4] === null || v[4] === undefined ? '' : String(v[4]) }); });
  const tables = new Map();
  for (const s of schema) {
    if (s.type !== 'table') continue;
    const virtual = /^\s*CREATE\s+VIRTUAL\s+TABLE/i.test(s.sql);
    tables.set(s.name.toLowerCase(), { name: s.name, rootPage: s.rootPage, sql: s.sql, columns: virtual ? [] : parseColumns(s.sql), withoutRowid: /\)\s*(STRICT\s*,\s*)?WITHOUT\s+ROWID\b[^)]*$/i.test(s.sql), virtual });
  }
  const resolve = (t) => {
    const tab = typeof t === 'string' ? tables.get(t.toLowerCase()) : t;
    if (!tab) throw new Error(`${fileName}: the database has no table "${t}".`);
    if (tab.virtual || !tab.rootPage) { const e = new Error(`${fileName}: table "${tab.name}" is a virtual table and has no rows of its own.`); e.code = 'VIRTUAL'; throw e; }
    if (tab.withoutRowid) { const e = new Error(`${fileName}: table "${tab.name}" is a WITHOUT ROWID table (stored as an index B-tree), which this reader does not read.`); e.code = 'WITHOUT_ROWID'; throw e; }
    return tab;
  };
  return {
    header, schema, tables,
    /** Walk the rows of a table; INTEGER PRIMARY KEY columns (stored as NULL) are filled with the rowid. */
    scan(t, cb) {
      const tab = resolve(t);
      const alias = tab.columns.findIndex((c) => c.rowidAlias);
      const n = tab.columns.length;
      scanRoot(tab.rootPage, (rowid, v) => {
        // rows written before ALTER TABLE ADD COLUMN are shorter than the table: the new columns read as NULL
        while (v.length < n) v.push(null);
        if (alias >= 0 && v[alias] === null) v[alias] = rowid;
        return cb(rowid, v);
      });
    },
    count(t) { return countRoot(resolve(t).rootPage); },
    columnsOf(t) { const tab = tables.get(String(t).toLowerCase()); return tab ? tab.columns : []; },
  };
}

/** Column list of a CREATE TABLE statement: `[{ name, type, rowidAlias }]`. */
export function parseColumns(sql) {
  const open = sql.indexOf('(');
  if (open < 0) return [];
  // body up to the matching parenthesis, split at top-level commas (quotes and nested parentheses respected)
  const parts = [];
  let depth = 0, cur = '', quote = null;
  for (let i = open + 1; i < sql.length; i++) {
    const ch = sql[i];
    if (quote) { cur += ch; if (ch === quote) { if (sql[i + 1] === quote && quote !== ']') { cur += sql[++i]; } else quote = null; } continue; }
    if (ch === '"' || ch === '\'' || ch === '`') { quote = ch; cur += ch; continue; }
    if (ch === '[') { quote = ']'; cur += ch; continue; }
    if (ch === '-' && sql[i + 1] === '-') { while (i < sql.length && sql[i] !== '\n') i++; cur += ' '; continue; }
    if (ch === '/' && sql[i + 1] === '*') { const e = sql.indexOf('*/', i + 2); i = e < 0 ? sql.length : e + 1; cur += ' '; continue; }
    if (ch === '(') depth++;
    if (ch === ')') { if (depth === 0) break; depth--; }
    if (ch === ',' && depth === 0) { parts.push(cur); cur = ''; continue; }
    cur += ch;
  }
  if (cur.trim()) parts.push(cur);
  const unq = (s) => (/^".*"$|^`.*`$|^'.*'$/s.test(s) ? s.slice(1, -1).split(s[0] + s[0]).join(s[0]) : /^\[.*\]$/s.test(s) ? s.slice(1, -1) : s);
  const cols = [];
  let tablePk = null;
  for (const raw of parts) {
    const p = raw.trim();
    if (!p) continue;
    const pk = /^(?:CONSTRAINT\s+\S+\s+)?PRIMARY\s+KEY\s*\(([^)]*)\)/i.exec(p);
    if (pk) { const names = pk[1].split(',').map((s) => unq(s.trim().replace(/\s+(ASC|DESC|COLLATE\s+\S+)\s*$/i, '').trim())); if (names.length === 1) tablePk = names[0].toLowerCase(); continue; }
    if (/^(CONSTRAINT|UNIQUE|CHECK|FOREIGN\s+KEY)\b/i.test(p)) continue;
    const m = /^("(?:[^"]|"")*"|`(?:[^`]|``)*`|\[[^\]]*\]|'(?:[^']|'')*'|[^\s(]+)\s*([\s\S]*)$/.exec(p);
    if (!m) continue;
    const rest = m[2] || '';
    const type = (/^([^\s(]+(?:\s+(?!PRIMARY|NOT|NULL|UNIQUE|CHECK|DEFAULT|COLLATE|REFERENCES|GENERATED|AS|CONSTRAINT)[^\s(]+)*)/i.exec(rest) || ['', ''])[1];
    const isType = type && !/^(PRIMARY|NOT|NULL|UNIQUE|CHECK|DEFAULT|COLLATE|REFERENCES|GENERATED|AS|CONSTRAINT)$/i.test(type);
    cols.push({ name: unq(m[1]), type: isType ? type.toUpperCase() : '', rowidAlias: /^INTEGER$/i.test(isType ? type : '') && /\bPRIMARY\s+KEY\b/i.test(rest) && !/\bPRIMARY\s+KEY\s+DESC\b/i.test(rest) });
  }
  if (tablePk) { const c = cols.find((q) => q.name.toLowerCase() === tablePk); if (c && c.type === 'INTEGER') c.rowidAlias = true; }
  return cols;
}

// ---------------------------------------------------------------------------------------------------------
// GeoPackageBinary and WKB
// ---------------------------------------------------------------------------------------------------------

const ENVELOPE_BYTES = [0, 32, 48, 48, 64];
const WKB_NAMES = { 1: 'Point', 2: 'LineString', 3: 'Polygon', 4: 'MultiPoint', 5: 'MultiLineString', 6: 'MultiPolygon', 7: 'GeometryCollection', 8: 'CircularString', 9: 'CompoundCurve', 10: 'CurvePolygon', 11: 'MultiCurve', 12: 'MultiSurface', 13: 'Curve', 14: 'Surface', 15: 'PolyhedralSurface', 16: 'TIN', 17: 'Triangle' };

/**
 * Header of a GeoPackageBinary blob.
 * @returns {{ version: number, srsId: number, empty: boolean, extended: boolean, envelope: number[]|null, wkbOffset: number }}
 */
export function gpkgHeader(blob) {
  if (blob.length < 8 || blob[0] !== 0x47 || blob[1] !== 0x50) throw new Error('not a GeoPackage geometry (the blob does not start with "GP")');
  const flags = blob[3];
  const le = (flags & 1) === 1;
  const env = (flags >> 1) & 7;
  if (env > 4) throw new Error(`invalid GeoPackage geometry envelope code ${env}`);
  const dv = new DataView(blob.buffer, blob.byteOffset, blob.byteLength);
  const nb = ENVELOPE_BYTES[env];
  if (blob.length < 8 + nb) throw new Error('truncated GeoPackage geometry header');
  let envelope = null;
  if (nb) { envelope = []; for (let i = 0; i < nb; i += 8) envelope.push(dv.getFloat64(8 + i, le)); }
  return { version: blob[2], srsId: dv.getInt32(4, le), empty: !!(flags & 0x10), extended: !!(flags & 0x20), envelope, wkbOffset: 8 + nb };
}

/**
 * Decode a WKB geometry (ISO 13249-3 / OGC 06-103 codes and PostGIS EWKB flags).
 * @param {Uint8Array} b bytes
 * @param {number} [offset]
 * @returns {{ type: string, hasZ: boolean, hasM: boolean, srid?: number, coords?: Float64Array, rings?: Float64Array[],
 *   parts?: object[], end: number }} coordinates flat x, y, z (z = 0 without Z; M dropped); `parts` for Multi*
 *   and collections; `end` is the offset after the geometry
 */
export function parseWkb(b, offset = 0) {
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength);
  let p = offset;
  const need = (n) => { if (p + n > b.length) throw new Error(`WKB geometry is truncated at byte ${p} of ${b.length}`); };
  const geom = (depth) => {
    if (depth > 32) throw new Error('WKB geometry is nested too deeply');
    need(5);
    const bo = b[p];
    if (bo > 1) throw new Error(`invalid WKB byte-order mark ${bo}`);
    const le = bo === 1;
    let code = dv.getUint32(p + 1, le);
    p += 5;
    let hasZ = false, hasM = false, srid;
    if (code & 0xe0000000) {
      // EWKB: dimensionality and SRID presence in the top bits
      hasZ = !!(code & 0x80000000); hasM = !!(code & 0x40000000);
      if (code & 0x20000000) { need(4); srid = dv.getInt32(p, le); p += 4; }
      code &= 0x1fffffff;
    }
    const dimCode = Math.floor(code / 1000);
    if (dimCode === 1 || dimCode === 3) hasZ = true;
    if (dimCode === 2 || dimCode === 3) hasM = true;
    const base = code % 1000;
    const type = WKB_NAMES[base];
    if (!type || dimCode > 3) throw new Error(`unknown WKB geometry type code ${code}`);
    const stride = 2 + (hasZ ? 1 : 0) + (hasM ? 1 : 0);
    const seq = (n) => {
      need(8 * stride * n);
      const out = new Float64Array(3 * n);
      for (let i = 0; i < n; i++) {
        out[3 * i] = dv.getFloat64(p, le); out[3 * i + 1] = dv.getFloat64(p + 8, le);
        out[3 * i + 2] = hasZ ? dv.getFloat64(p + 16, le) : 0;
        p += 8 * stride;
      }
      return out;
    };
    const count = () => { need(4); const n = dv.getUint32(p, le); p += 4; return n; };
    const g = { type, hasZ, hasM };
    if (srid !== undefined) g.srid = srid;
    switch (base) {
      case 1: g.coords = seq(1); break;
      case 2: case 8: g.coords = seq(count()); break;
      case 3: case 17: { const n = count(); g.rings = []; for (let i = 0; i < n; i++) g.rings.push(seq(count())); break; }
      case 4: case 5: case 6: case 7: case 9: case 10: case 11: case 12: case 15: case 16: { const n = count(); g.parts = []; for (let i = 0; i < n; i++) g.parts.push(geom(depth + 1)); break; }
      default: throw new Error(`WKB geometry type ${type} cannot be decoded`);
    }
    return g;
  };
  const g = geom(0);
  g.end = p;
  return g;
}

const dropClosing = (c) => { const n = c.length; return n >= 6 && c[0] === c[n - 3] && c[1] === c[n - 2] && c[2] === c[n - 1] ? c.subarray(0, n - 3) : c; };

function addGeometry(F, g, props, name, curved) {
  switch (g.type) {
    case 'Point': F.addPoint(g.coords[0], g.coords[1], g.coords[2], props, name); break;
    case 'LineString': F.addLine(g.coords, props, name); break;
    case 'CircularString': curved.n++; F.addLine(g.coords, props, name); break;
    case 'Polygon': case 'Triangle': if (g.rings.length) F.addPolygon(g.rings.map(dropClosing), props, name); break;
    default: for (const s of g.parts || []) addGeometry(F, s, props, name, curved);
  }
}

// ---------------------------------------------------------------------------------------------------------
// GeoPackage
// ---------------------------------------------------------------------------------------------------------

const rowsOf = (db, table) => {
  const cols = db.columnsOf(table).map((c) => c.name.toLowerCase());
  const out = [];
  db.scan(table, (rowid, v) => { const r = {}; cols.forEach((c, i) => { r[c] = v[i]; }); out.push(r); });
  return out;
};

function srsInfo(srs, srsId) {
  if (srsId === 0) return { epsg: null, name: 'undefined geographic coordinate system', kind: 'geographic', unit: 'deg', source: 'srs_id 0' };
  if (srsId === -1 || srsId === null || srsId === undefined) return null;
  const row = srs.get(srsId);
  let info = null;
  if (row && /^epsg$/i.test(String(row.organization || ''))) info = crsFromName(`EPSG:${row.organization_coordsys_id}`);
  else if (!row && srsId > 1023) info = crsFromName(`EPSG:${srsId}`);
  if ((!info || info.kind === 'unknown') && row && typeof row.definition === 'string' && row.definition.length > 12 && !/^undefined$/i.test(row.definition)) {
    const w = parsePrj(row.definition);
    if (w && w.kind !== 'unknown') info = { ...w, epsg: w.epsg ?? (info && info.epsg) ?? null, name: w.name || (row && row.srs_name) || '', source: (info && info.source) || (w.epsg ? `EPSG:${w.epsg}` : String(row.srs_name || `srs_id ${srsId}`)) };
  }
  if (info && row && row.srs_name && (!info.name || /^EPSG:\d+$/.test(info.name))) info = { ...info, name: String(row.srs_name) };
  return info || { epsg: null, name: row ? String(row.srs_name || `srs_id ${srsId}`) : `srs_id ${srsId}`, kind: 'unknown', unit: null, source: `srs_id ${srsId}` };
}

/**
 * Import the vector layers of a GeoPackage.
 * @param {{ name?: string, bytes: Uint8Array }} input the .gpkg file
 * @param {object} [opts] the vector-importer options (`projection`, `origin`, `utmZone`, `crs`, `units`,
 *   `vertical`, `polygons`) plus:
 * @param {string|string[]} [opts.layers] feature table name(s) to read (default: all feature tables)
 * @param {number} [opts.maxFeatures] stop reading a layer after this many rows (default: no limit)
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} meta: { format, application,
 *   sqlite: { pageSize, encoding, pages, journal }, layers: [{ name, identifier, dataType, geometryColumn,
 *   geometryType, srsId, crs, rows, features, models }], tiles: [...], crs, units, counts, bbox }
 */
export async function parse(input, opts = {}) {
  const fname = (input && input.name) || 'GeoPackage';
  if (!input || !input.bytes) throw new Error(`${fname}: a GeoPackage is a binary SQLite file — expected its bytes.`);
  const bytes = input.bytes;
  const db = openSqlite(bytes, fname);
  const warnings = [];
  if (db.header.wal) warnings.push('The database is in WAL journal mode: it was read as of its last checkpoint. Changes that are still in a side-car "-wal" file are not included — run "PRAGMA wal_checkpoint(TRUNCATE)" (or open and close it in a GIS) before importing.');
  if (db.header.partialLastPage) warnings.push('The file length is not a whole number of pages; the incomplete last page was ignored.');
  const has = (t) => db.tables.has(t);
  if (!has('gpkg_contents')) {
    const names = [...db.tables.values()].map((t) => t.name);
    if (has('geometry_columns') && has('spatial_ref_sys')) throw new Error(`${fname}: this SQLite database is a SpatiaLite database, not a GeoPackage (it has geometry_columns but no gpkg_contents). Convert it with "ogr2ogr -f GPKG out.gpkg in.sqlite" or export GeoJSON.`);
    throw new Error(`${fname}: this SQLite database is not a GeoPackage (it has no gpkg_contents table). Tables: ${names.slice(0, 20).join(', ') || 'none'}${names.length > 20 ? ', …' : ''}.`);
  }
  const contents = rowsOf(db, 'gpkg_contents');
  const geomCols = has('gpkg_geometry_columns') ? rowsOf(db, 'gpkg_geometry_columns') : [];
  const srs = new Map();
  if (has('gpkg_spatial_ref_sys')) for (const r of rowsOf(db, 'gpkg_spatial_ref_sys')) srs.set(Number(r.srs_id), r);

  const wanted = opts.layers ? new Set((Array.isArray(opts.layers) ? opts.layers : [opts.layers]).map((s) => String(s).toLowerCase())) : (opts.layer ? new Set([String(opts.layer).toLowerCase()]) : null);
  const layers = [], tiles = [], other = [];
  const collected = [];
  const unread = [];
  for (const c of contents) {
    const tname = String(c.table_name);
    const type = String(c.data_type || '').toLowerCase();
    if (type === 'tiles' || type === '2d-gridded-coverage') {
      let n = null;
      try { n = db.count(tname); } catch { /* the table is listed but missing */ }
      tiles.push({ name: tname, identifier: c.identifier ?? null, dataType: type, tiles: n, srsId: c.srs_id ?? null });
      continue;
    }
    if (type !== 'features') { other.push({ name: tname, dataType: type || null }); continue; }
    if (wanted && !wanted.has(tname.toLowerCase())) continue;
    const gc = geomCols.find((g) => String(g.table_name).toLowerCase() === tname.toLowerCase());
    const tab = db.tables.get(tname.toLowerCase());
    if (!gc || !tab) { warnings.push(`Feature layer "${tname}" is listed in gpkg_contents but ${tab ? 'has no entry in gpkg_geometry_columns' : 'its table is missing'}; it was skipped.`); continue; }
    if (tab.withoutRowid) { warnings.push(`Feature layer "${tname}" is a WITHOUT ROWID table, which this reader does not read; it was skipped. Export it with "ogr2ogr -f GeoJSON out.geojson in.gpkg ${tname}".`); unread.push(tname); continue; }
    const cols = tab.columns;
    const gi = cols.findIndex((q) => q.name.toLowerCase() === String(gc.column_name).toLowerCase());
    if (gi < 0) { warnings.push(`Feature layer "${tname}": geometry column "${gc.column_name}" is not a column of the table; the layer was skipped.`); continue; }
    const F = new Features(tname);
    const nameCol = cols.findIndex((q, i) => i !== gi && /^(name|label|title|id_name)$/i.test(q.name));
    let rows = 0, nullGeom = 0, bad = 0, firstBad = null, srsSeen = null, mixedSrs = false, mDropped = false, stopped = false;
    const curved = { n: 0 };
    const types = {};
    try {
      db.scan(tab, (rowid, v) => {
        if (opts.maxFeatures > 0 && rows >= opts.maxFeatures) { stopped = true; return false; }
        rows++;
        const blob = v[gi];
        if (!(blob instanceof Uint8Array) || !blob.length) { nullGeom++; return true; }
        let g;
        try {
          const h = gpkgHeader(blob);
          if (srsSeen === null) srsSeen = h.srsId; else if (srsSeen !== h.srsId) mixedSrs = true;
          if (h.extended) throw new Error('extended (non-standard) geometry encoding');
          if (h.empty) { nullGeom++; return true; }
          g = parseWkb(blob, h.wkbOffset);
        } catch (e) { bad++; if (!firstBad) firstBad = `row ${rowid}: ${e.message}`; return true; }
        if (g.hasM) mDropped = true;
        types[g.type] = (types[g.type] || 0) + 1;
        const props = {};
        for (let i = 0; i < cols.length; i++) {
          if (i === gi) continue;
          const q = v[i];
          if (q === null || q === undefined || q instanceof Uint8Array) continue;
          props[cols[i].name] = q;
        }
        const nm = nameCol >= 0 && v[nameCol] !== null && v[nameCol] !== undefined && v[nameCol] !== '' ? String(v[nameCol]) : undefined;
        addGeometry(F, g, props, g.type === 'Point' || g.type === 'MultiPoint' ? undefined : nm, curved);
        return true;
      });
    } catch (e) {
      if (e && e.code === 'INDEX_BTREE') { warnings.push(`Feature layer "${tname}" is stored as an index B-tree (WITHOUT ROWID), which this reader does not read; it was skipped.`); unread.push(tname); continue; }
      throw e;
    }
    const srsId = gc.srs_id !== null && gc.srs_id !== undefined ? Number(gc.srs_id) : srsSeen;
    const layer = { name: tname, identifier: c.identifier ?? null, description: c.description || null, dataType: type, geometryColumn: String(gc.column_name), geometryType: String(gc.geometry_type_name || ''), srsId, rows, features: 0, geometryTypes: types, columns: cols.filter((q, i) => i !== gi).map((q) => ({ name: q.name, type: q.type })) };
    if (nullGeom) warnings.push(`Layer "${tname}": ${nullGeom} row(s) have no geometry (NULL or empty) and were skipped.`);
    if (bad) warnings.push(`Layer "${tname}": ${bad} geometry blob(s) could not be decoded and were skipped (first: ${firstBad}).`);
    if (mixedSrs) warnings.push(`Layer "${tname}": the geometries carry different srs_id values; the one in gpkg_geometry_columns (${srsId}) was used for all.`);
    if (mDropped) warnings.push(`Layer "${tname}": measure (M) values were dropped.`);
    if (curved.n) warnings.push(`Layer "${tname}": ${curved.n} circular arc string(s) were read as straight segments through their control points.`);
    if (stopped) warnings.push(`Layer "${tname}": reading stopped after maxFeatures = ${opts.maxFeatures} rows.`);
    layers.push(layer);
    collected.push({ F, layer, crs: srsInfo(srs, srsId) });
  }
  if (wanted) for (const w of wanted) if (!collected.some((q) => q.layer.name.toLowerCase() === w) && !unread.includes(w)) warnings.push(`Layer "${w}" was asked for but is not a feature table of this GeoPackage (feature tables: ${contents.filter((c) => String(c.data_type).toLowerCase() === 'features').map((c) => c.table_name).join(', ') || 'none'}).`);
  for (const t of tiles) warnings.push(`"${t.name}" is a ${t.dataType === 'tiles' ? 'raster tile pyramid' : 'tiled gridded coverage (elevation tiles)'}${t.tiles !== null ? ` with ${t.tiles} tile(s)` : ''}: raster tiles are not decoded. Export it as GeoTIFF, for example "gdal_translate in.gpkg out.tif${t.dataType === 'tiles' ? '' : ` -oo TABLE=${t.name}`}".`);
  for (const s of db.schema) if (s.type === 'table' && /\bWITHOUT\s+ROWID\b/i.test(s.sql) && !unread.includes(s.name) && !/^(gpkg_|rtree_|sqlite_)/i.test(s.name)) warnings.push(`Table "${s.name}" is a WITHOUT ROWID table; such tables are not read (it is not a feature layer, so nothing is missing from the geometry).`);

  // One projection origin for every longitude/latitude layer, so that the layers stay aligned.
  let common = null;
  if (!(opts.origin && Number.isFinite(opts.origin.lon) && Number.isFinite(opts.origin.lat))) {
    const bb = [Infinity, Infinity, -Infinity, -Infinity];
    for (const q of collected) {
      const geo = opts.crs ? /4326|geographic|CRS84/i.test(String(opts.crs)) : q.crs && q.crs.kind === 'geographic';
      if (!geo || q.F.empty) continue;
      const b = q.F.b;
      if (b[0] < bb[0]) bb[0] = b[0]; if (b[1] < bb[1]) bb[1] = b[1]; if (b[2] > bb[2]) bb[2] = b[2]; if (b[3] > bb[3]) bb[3] = b[3];
    }
    if (bb[0] !== Infinity) common = { lon: (bb[0] + bb[2]) / 2, lat: (bb[1] + bb[3]) / 2 };
  }
  const models = [];
  let crsMeta = null, units = null, vertical = null;
  const elementTypes = {};
  for (const q of collected) {
    if (q.F.empty) { if (q.layer.rows) warnings.push(`Layer "${q.layer.name}" has ${q.layer.rows} row(s) but no readable geometry.`); q.layer.crs = q.crs ? (q.crs.source || q.crs.name) : null; q.layer.models = 0; continue; }
    const r = q.F.build(common ? { ...opts, origin: common } : opts, { format: 'geopackage', fileName: fname, crs: q.crs, warnings: [], expected: 'expected GeoPackageBinary geometries' });
    if (common && r.meta.crs && r.meta.crs.originBasis === 'option') r.meta.crs.originBasis = 'data centre of all layers';
    for (const m of r.models) { m.layer = q.layer.name; models.push(m); }
    for (const w of r.warnings) warnings.push(`Layer "${q.layer.name}": ${w}`);
    q.layer.crs = r.meta.crs;
    q.layer.features = r.meta.features;
    q.layer.models = r.models.length;
    for (const k in r.meta.elementTypes) elementTypes[k] = (elementTypes[k] || 0) + r.meta.elementTypes[k];
    if (!crsMeta) { crsMeta = r.meta.crs; units = r.meta.units; vertical = r.meta.vertical; } else if (JSON.stringify(crsMeta) !== JSON.stringify(r.meta.crs) && !warnings.some((w) => /different coordinate systems/.test(w))) warnings.push('The layers use different coordinate systems; each was converted on its own (see meta.layers[].crs), so they may not overlay.');
  }
  if (!models.length) {
    if (!tiles.length && !layers.length) throw new Error(`${fname}: the GeoPackage has no feature layers${wanted ? ' with the requested name(s)' : ''} (gpkg_contents lists: ${contents.map((c) => `${c.table_name} [${c.data_type}]`).join(', ') || 'nothing'}).`);
    if (!tiles.length && !warnings.length) warnings.push('The feature layers of the GeoPackage are empty.');
  }
  return finish(models, warnings, {
    format: 'geopackage', units, crs: crsMeta, vertical, elementTypes,
    application: db.header.applicationId === 0x47504b47 ? `GeoPackage (user_version ${db.header.userVersion})` : db.header.applicationId === 0x47503130 ? 'GeoPackage 1.0' : db.header.applicationId === 0x47503131 ? 'GeoPackage 1.1' : null,
    sqlite: { pageSize: db.header.pageSize, encoding: db.header.encoding, pages: db.header.pages, journal: db.header.wal ? 'wal' : 'rollback' },
    layers, tiles, otherContents: other, reader: 'built-in', name: baseName(fname),
  });
}
