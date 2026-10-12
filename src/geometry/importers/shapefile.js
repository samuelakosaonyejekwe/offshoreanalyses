// ESRI Shapefile: the .shp geometry file with its optional .dbf attribute table, .prj coordinate system and
// .cpg code page — supplied as a file set or zipped.
//
// Shape types read: Point, MultiPoint, PolyLine, Polygon and their Z and M variants. MultiPatch (31) records
// are counted as skipped. Polygon holes are attached to the outer ring that contains them (outer rings are
// clockwise in a shapefile). The .shx index is not needed. Attributes of type C, N, F, L and D are read;
// memo fields are ignored.

import { baseName, findFile, decodeText, extOf, tally } from './common.js';
import { parsePrj, crsLabel } from '../units.js';
import { Features } from './vector.js';

function readDbf(bytes, encoding) {
  if (!bytes || bytes.length < 33) return null;
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const n = dv.getUint32(4, true), headerSize = dv.getUint16(8, true), recSize = dv.getUint16(10, true);
  let dec;
  try { dec = new TextDecoder(encoding || 'windows-1252'); } catch { dec = new TextDecoder('utf-8'); }
  const fields = [];
  let off = 1;
  for (let p = 32; p + 32 <= headerSize && bytes[p] !== 0x0d; p += 32) {
    let nm = '';
    for (let i = 0; i < 11 && bytes[p + i] !== 0; i++) nm += String.fromCharCode(bytes[p + i]);
    const len = bytes[p + 16];
    fields.push({ name: nm.trim(), type: String.fromCharCode(bytes[p + 11]), off, len });
    off += len;
  }
  const record = (i) => {
    const base = headerSize + i * recSize;
    if (i >= n || base + recSize > bytes.length || bytes[base] === 0x2a) return null;
    const o = {};
    for (const f of fields) {
      const s = dec.decode(bytes.subarray(base + f.off, base + f.off + f.len)).replace(/\u0000/g, '').trim();
      if (s === '') continue;
      if (f.type === 'N' || f.type === 'F') { const v = Number(s); if (v === v) o[f.name] = v; } else if (f.type === 'L') o[f.name] = /^[YyTt]/.test(s); else o[f.name] = s;
    }
    return o;
  };
  return { count: n, fields: fields.map((f) => f.name), record };
}

function ringArea2(c) {
  let a = 0;
  for (let i = 0, n = c.length / 3, j = n - 1; i < n; j = i++) a += c[3 * j] * c[3 * i + 1] - c[3 * i] * c[3 * j + 1];
  return a;
}

function pointInRing(x, y, c) {
  let inside = false;
  for (let i = 0, n = c.length / 3, j = n - 1; i < n; j = i++) {
    const xi = c[3 * i], yi = c[3 * i + 1], xj = c[3 * j], yj = c[3 * j + 1];
    if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/** ESRI Shapefile (file set, or the bare .shp). */
export function parseShapefile(input, opts = {}) {
  let shp = input.bytes || null, files = input.files || null;
  let stem = baseName(input.name);
  if (!shp && files) {
    let key = input.name && files.has(input.name) ? input.name : null;
    if (!key) for (const k of files.keys()) if (extOf(k) === 'shp') { key = k; break; }
    if (key) { shp = files.get(key); stem = baseName(key); }
  }
  if (!shp) throw new Error(`${input.name || 'shapefile'}: no .shp file was supplied (a shapefile is a set: .shp geometry, .dbf attributes, .prj coordinate system — upload them together or zipped).`);
  if (shp.length < 100) throw new Error(`${input.name || 'shapefile'}: the .shp file is shorter than its 100-byte header.`);
  const dv = new DataView(shp.buffer, shp.byteOffset, shp.byteLength);
  if (dv.getInt32(0, false) !== 9994) throw new Error(`${input.name || 'shapefile'}: not a shapefile (expected file code 9994 at the start of the .shp).`);
  const side = (ext) => { if (!files) return null; const f = findFile(files, `${stem}.${ext}`); if (f) return f; for (const [k, v] of files) if (extOf(k) === ext) return v; return null; };
  const cpg = side('cpg');
  const enc = cpg ? decodeText(cpg).trim().toLowerCase().replace(/^(utf)(\d)/, '$1-$2') : null;
  const dbf = readDbf(side('dbf'), enc && enc !== 'system' ? enc : null);
  const prjBytes = side('prj');
  const prj = prjBytes ? parsePrj(decodeText(prjBytes)) : null;
  const warnings = [];
  const F = new Features(stem);
  const fileType = dv.getInt32(32, true);
  let p = 100, rec = 0;
  const flat = (o, n, zs) => { const c = new Float64Array(n * 3); for (let i = 0; i < n; i++) { c[3 * i] = dv.getFloat64(o + 16 * i, true); c[3 * i + 1] = dv.getFloat64(o + 16 * i + 8, true); c[3 * i + 2] = zs >= 0 ? dv.getFloat64(zs + 8 * i, true) : 0; } return c; };
  while (p + 12 <= shp.length) {
    const len = dv.getInt32(p + 4, false) * 2;
    const o = p + 8;
    if (len < 4 || o + len > shp.length) { warnings.push(`The .shp file ends inside record ${rec + 1}; reading stopped there.`); break; }
    const type = dv.getInt32(o, true);
    const props = dbf ? dbf.record(rec) : null;
    const base = type % 10, hasZ = type >= 10 && type < 20;
    if (type === 0) tally(F.skipped, 'nullShape');
    else if (type === 31) tally(F.skipped, 'multiPatch');
    else if (base === 1) F.addPoint(dv.getFloat64(o + 4, true), dv.getFloat64(o + 12, true), hasZ && len >= 28 ? dv.getFloat64(o + 20, true) : 0, props);
    else if (base === 8) {
      const n = dv.getInt32(o + 36, true);
      const zs = hasZ ? o + 40 + 16 * n + 16 : -1;
      const c = flat(o + 40, n, zs);
      for (let i = 0; i < n; i++) F.addPoint(c[3 * i], c[3 * i + 1], c[3 * i + 2], props);
    } else if (base === 3 || base === 5) {
      const nParts = dv.getInt32(o + 36, true), nPts = dv.getInt32(o + 40, true);
      const po = o + 44, xy = po + 4 * nParts;
      const zs = hasZ ? xy + 16 * nPts + 16 : -1;
      const parts = [];
      for (let k = 0; k < nParts; k++) {
        const s = dv.getInt32(po + 4 * k, true), e = k + 1 < nParts ? dv.getInt32(po + 4 * (k + 1), true) : nPts;
        if (e > s) parts.push(flat(xy + 16 * s, e - s, zs >= 0 ? zs + 8 * s : -1));
      }
      if (base === 3) for (const c of parts) F.addLine(c, props);
      else {
        const rings = parts.map((c) => { const n = c.length; return n >= 6 && c[0] === c[n - 3] && c[1] === c[n - 2] ? c.subarray(0, n - 3) : c; }).filter((c) => c.length >= 9);
        // clockwise rings (negative signed area) are outer boundaries, counter-clockwise ones are holes
        const outers = [], holes = [];
        for (const r of rings) (ringArea2(r) <= 0 ? outers : holes).push(r);
        if (!outers.length && holes.length) { outers.push(...holes); holes.length = 0; }
        const groups = outers.map((r) => [r]);
        for (const h of holes) { const g = groups.find((gr) => pointInRing(h[0], h[1], gr[0])) || groups[0]; g.push(h); }
        for (const g of groups) F.addPolygon(g, props);
      }
    } else tally(F.skipped, `shapeType${type}`);
    rec++;
    p = o + len;
  }
  if (!dbf) warnings.push('No .dbf attribute table was supplied; features carry no attributes.');
  if (!prj) warnings.push('No .prj file was supplied; coordinates are taken as projected metres unless crs is given.');
  if (F.skipped.multiPatch) warnings.push(`${F.skipped.multiPatch} MultiPatch record(s) were skipped (not read).`);
  const crs = prj && prj.kind !== 'unknown' ? { epsg: prj.epsg, name: prj.name, kind: prj.kind === 'local' ? 'projected' : prj.kind, unit: prj.unit, unitScale: prj.unitScale || 1, source: crsLabel(prj) || prj.name, utm: prj.utm } : null;
  return F.build(opts, {
    format: 'shapefile', fileName: input.name, crs, warnIfDegrees: !prj, warnings, expected: 'expected Point, PolyLine or Polygon records in the .shp',
    extra: { shapeType: fileType, records: rec, attributeFields: dbf ? dbf.fields : [], verticalCrs: prj ? prj.vertical : null },
  });
}
