import test from 'node:test';
import assert from 'node:assert/strict';
import { deflateSync } from 'node:zlib';
import { parseEsriAscii, parseSurferGrd, blockAverage } from '../../src/geometry/importers/grids.js';
import { parseGeoTiff, lzwDecode, packBitsDecode } from '../../src/geometry/importers/geotiff.js';
import { parseGeoJson, parseWkt, crsFromName } from '../../src/geometry/importers/vector.js';
import { parseShapefile } from '../../src/geometry/importers/shapefile.js';
import { parseLandXml, parseGml, parseKml, parseGpx, profileElevation } from '../../src/geometry/importers/gisxml.js';
import { parseNetcdf } from '../../src/geometry/importers/netcdf.js';
import { parseLas, parsePts, parsePtx, parseXyzRgb } from '../../src/geometry/importers/pointcloud.js';
import { makeZip } from '../../src/geometry/zip.js';
import { utmForward } from '../../src/geometry/units.js';
import { enc, close, closeArr, concat, assertContract, surfaceArea, checkElements, assertTruncationSafe } from './fixtures/helpers.js';

// ---- grids -----------------------------------------------------------------------------------------------

test('ESRI ASCII grid: corner registration, north-up rows flipped, NODATA → NaN, depth option', () => {
  const text = 'ncols 3\nnrows 2\nxllcorner 100\nyllcorner 200\ncellsize 10\nNODATA_value -9999\n1 2 3\n4 -9999 6\n';
  const r = parseEsriAscii({ name: 'dem.asc', bytes: enc(text) });
  assertContract(r, ['grid']);
  const g = r.models[0];
  assert.deepEqual([g.nx, g.ny, g.x0, g.y0, g.dx, g.dy], [3, 2, 105, 205, 10, 10]);
  // file row 2 is the southern row → grid row 0
  assert.deepEqual([...g.z].map((v) => (Number.isNaN(v) ? null : v)), [4, null, 6, 1, 2, 3]);
  assert.equal(r.meta.skipped.noDataCells, 1);
  assert.deepEqual(r.meta.bbox, { min: [105, 205, 1], max: [125, 215, 6] });
  const d = parseEsriAscii({ name: 'dem.asc', text }, { vertical: 'depth', units: 'ft', zUnits: 'm' });
  assert.equal(d.models[0].z[0], -4);
  close(d.models[0].dx, 3.048, 1e-12);
  assert.throws(() => parseEsriAscii({ name: 'x.asc', text: 'hello world\n1 2 3\n' }), /expected header lines ncols, nrows/);
});

test('ESRI ASCII grid with a geographic .prj: spacing converted to metres on the local tangent plane', () => {
  const text = 'ncols 2\nnrows 2\nxllcenter 3.0\nyllcenter 54.0\ncellsize 0.01\n-30 -31\n-32 -33\n';
  const prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]';
  const r = parseEsriAscii({ name: 'bathy.asc', text, files: new Map([['bathy.prj', enc(prj)]]) }, { origin: { lon: 3.0, lat: 54.0 } });
  const g = r.models[0];
  assert.equal(g.crs, 'EPSG:4326');
  close(g.x0, 0, 1e-6);
  // 0.01° of latitude at 54°N: M·Δφ, M = a(1 − e²)/(1 − e² sin²φ)^1.5
  const f = 1 / 298.257223563, e2 = f * (2 - f), sp = Math.sin((54 * Math.PI) / 180);
  close(g.dy, ((6378137 * (1 - e2)) / (1 - e2 * sp * sp) ** 1.5) * ((0.01 * Math.PI) / 180), 1e-6);
  close(g.dy, 1113.0, 0.5, 'hand value');
  assert.equal(r.meta.crs.projection, 'local-tangent-plane');
  const plain = parseEsriAscii({ name: 'bathy.asc', text });
  assert.ok(plain.warnings.some((w) => /would fit longitude\/latitude/.test(w)));
  assert.equal(plain.models[0].dx, 0.01);
});

test('block averaging to maxCells ignores no-data and keeps the lattice registered', () => {
  const z = Float32Array.from([1, 3, NaN, 7, 5, 7, 9, 11, 2, 2, 4, 4, 2, 2, 4, 8]);
  const b = blockAverage(z, 4, 4, 2);
  assert.deepEqual([b.nx, b.ny, ...b.z], [2, 2, 4, 9, 2, 5]);
  const rows = [];
  for (let r = 0; r < 4; r++) rows.push([0, 1, 2, 3].map((c) => z[(3 - r) * 4 + c]).map((v) => (Number.isNaN(v) ? -9999 : v)).join(' '));
  const g = parseEsriAscii({ name: 'd.asc', text: `ncols 4\nnrows 4\nxllcorner 0\nyllcorner 0\ncellsize 1\nnodata_value -9999\n${rows.join('\n')}\n` }, { maxCells: 4 });
  assert.deepEqual([g.models[0].nx, g.models[0].ny, g.models[0].x0, g.models[0].y0, g.models[0].dx], [2, 2, 1, 1, 2]);
  assert.deepEqual([...g.models[0].z], [4, 9, 2, 5]);
  assert.equal(g.meta.decimation.factor, 2);
  assert.ok(g.warnings.some((w) => /block-averaged by 2/.test(w)));
});

test('Surfer grids: DSAA text, DSBB and DSRR binary', () => {
  const vals = [1, 2, 3, 4, 1.70141e38, 6];
  const a = parseSurferGrd({ name: 's.grd', text: `DSAA\n3 2\n0 20\n100 110\n1 6\n1 2 3\n4 1.70141e38 6\n` });
  assertContract(a, ['grid']);
  assert.deepEqual([a.models[0].nx, a.models[0].ny, a.models[0].x0, a.models[0].y0, a.models[0].dx, a.models[0].dy], [3, 2, 0, 100, 10, 10]);
  assert.deepEqual([...a.models[0].z].map((v) => (Number.isNaN(v) ? null : v)), [1, 2, 3, 4, null, 6]);
  const b = new Uint8Array(56 + 24);
  const dv = new DataView(b.buffer);
  b.set(enc('DSBB'));
  dv.setInt16(4, 3, true); dv.setInt16(6, 2, true);
  [0, 20, 100, 110, 1, 6].forEach((v, i) => dv.setFloat64(8 + 8 * i, v, true));
  vals.forEach((v, i) => dv.setFloat32(56 + 4 * i, v, true));
  const rb = parseSurferGrd({ name: 's.grd', bytes: b });
  assert.deepEqual([...rb.models[0].z].map((v) => (Number.isNaN(v) ? null : v)), [1, 2, 3, 4, null, 6]);
  const c = new Uint8Array(12 + 8 + 72 + 8 + 48);
  const cv = new DataView(c.buffer);
  c.set(enc('DSRR')); cv.setInt32(4, 4, true); cv.setInt32(8, 2, true);
  c.set(enc('GRID'), 12); cv.setInt32(16, 72, true);
  cv.setInt32(20, 2, true); cv.setInt32(24, 3, true);
  [0, 100, 10, 10, 1, 6, 0, 1.70141e38].forEach((v, i) => cv.setFloat64(28 + 8 * i, v, true));
  c.set(enc('DATA'), 92); cv.setInt32(96, 48, true);
  vals.forEach((v, i) => cv.setFloat64(100 + 8 * i, v, true));
  const rc = parseSurferGrd({ name: 's.grd', bytes: c });
  assert.deepEqual([rc.models[0].dx, ...[...rc.models[0].z].map((v) => (Number.isNaN(v) ? null : v))], [10, 1, 2, 3, 4, null, 6]);
  assert.throws(() => parseSurferGrd({ name: 'x.grd', text: 'NOPE\n' }), /expected the file to start with DSAA, DSBB or DSRR/);
});

// ---- GeoTIFF ---------------------------------------------------------------------------------------------

function lzwEncode(data) {
  const out = [];
  let buf = 0, cnt = 0, width = 9, next = 258, cur = -1;
  const table = new Map();
  const put = (code) => { buf = (buf << width) | code; cnt += width; while (cnt >= 8) { out.push((buf >> (cnt - 8)) & 255); cnt -= 8; } buf &= (1 << cnt) - 1; };
  put(256);
  for (const b of data) {
    if (cur < 0) { cur = b; continue; }
    const key = cur * 256 + b;
    if (table.has(key)) { cur = table.get(key); continue; }
    put(cur);
    table.set(key, next++);
    if (next === 4094) { put(256); table.clear(); next = 258; width = 9; } else if (next > (1 << width) - 1 && width < 12) width++;
    cur = b;
  }
  if (cur >= 0) put(cur);
  put(257);
  if (cnt) out.push((buf << (8 - cnt)) & 255);
  return Uint8Array.from(out);
}

function packBitsEncode(data) {
  const out = [];
  for (let i = 0; i < data.length;) {
    let run = 1;
    while (i + run < data.length && data[i + run] === data[i] && run < 128) run++;
    if (run >= 3) { out.push(257 - run, data[i]); i += run; } else { const n = Math.min(128, data.length - i); let k = 0; while (k < n && !(k + 2 < n && data[i + k] === data[i + k + 1] && data[i + k] === data[i + k + 2])) k++; k = Math.max(1, k); out.push(k - 1, ...data.subarray(i, i + k)); i += k; }
  }
  return Uint8Array.from(out);
}

/** Minimal GeoTIFF writer for fixtures (single band, classic or BigTIFF, either byte order). */
function makeTiff({ width, height, values, bits = 32, format = 3, tile = null, rowsPerStrip = height, compression = 1, predictor = 1, le = true, big = false, geo = {}, nodata }) {
  const bps = bits / 8;
  const setter = format === 3 ? (bits === 32 ? 'setFloat32' : 'setFloat64') : format === 2 ? `setInt${bits}` : `setUint${bits}`;
  const bw = tile ? tile[0] : width, bh = tile ? tile[1] : rowsPerStrip;
  const across = Math.ceil(width / bw), down = Math.ceil(height / bh);
  const blocks = [];
  for (let by = 0; by < down; by++) for (let bx = 0; bx < across; bx++) {
    const rows = tile ? bh : Math.min(bh, height - by * bh);
    const raw = new Uint8Array(bw * rows * bps);
    const dv = new DataView(raw.buffer);
    for (let r = 0; r < rows; r++) for (let c = 0; c < bw; c++) {
      const x = bx * bw + c, y = by * bh + r;
      const v = x < width && y < height ? values[y * width + x] : 0;
      if (bps === 1) dv[setter](r * bw + c, v); else dv[setter]((r * bw + c) * bps, v, le);
    }
    if (predictor === 2) for (let r = 0; r < rows; r++) for (let c = bw - 1; c >= 1; c--) {
      const o = (r * bw + c) * bps;
      if (bps === 1) raw[o] = (raw[o] - raw[o - 1]) & 255; else if (bps === 2) dv.setUint16(o, (dv.getUint16(o, le) - dv.getUint16(o - 2, le)) & 0xffff, le); else dv.setUint32(o, (dv.getUint32(o, le) - dv.getUint32(o - 4, le)) >>> 0, le);
    }
    if (predictor === 3) for (let r = 0; r < rows; r++) {
      const row = raw.subarray(r * bw * bps, (r + 1) * bw * bps), tmp = new Uint8Array(row.length);
      for (let c = 0; c < bw; c++) for (let k = 0; k < bps; k++) tmp[k * bw + c] = row[c * bps + (le ? bps - 1 - k : k)];
      for (let i = tmp.length - 1; i >= 1; i--) tmp[i] = (tmp[i] - tmp[i - 1]) & 255;
      row.set(tmp);
    }
    blocks.push(compression === 5 ? lzwEncode(raw) : compression === 8 ? new Uint8Array(deflateSync(Buffer.from(raw))) : compression === 32773 ? packBitsEncode(raw) : raw);
  }
  const entries = [[256, 4, [width]], [257, 4, [height]], [258, 3, [bits]], [259, 3, [compression]], [262, 3, [1]], [277, 3, [1]], [339, 3, [format]]];
  if (predictor !== 1) entries.push([317, 3, [predictor]]);
  const head = big ? 16 : 8;
  let pos = head;
  const offs = blocks.map((b) => { const o = pos; pos += b.length; return o; });
  if (tile) entries.push([322, 4, [bw]], [323, 4, [bh]], [324, 4, offs], [325, 4, blocks.map((b) => b.length)]);
  else entries.push([278, 4, [bh]], [273, 4, offs], [279, 4, blocks.map((b) => b.length)]);
  if (geo.scale) entries.push([33550, 12, geo.scale]);
  if (geo.tie) entries.push([33922, 12, geo.tie]);
  if (geo.transform) entries.push([34264, 12, geo.transform]);
  if (geo.keys) entries.push([34735, 3, [1, 1, 0, geo.keys.length, ...geo.keys.flat()]]);
  if (nodata !== undefined) entries.push([42113, 2, [...enc(String(nodata)), 0]]);
  entries.sort((a, b) => a[0] - b[0]);
  const SZ = { 2: 1, 3: 2, 4: 4, 12: 8 };
  const esz = big ? 20 : 12, inline = big ? 8 : 4;
  const ifdOff = pos + (pos % 2);
  const ifdLen = (big ? 8 : 2) + entries.length * esz + (big ? 8 : 4);
  let extra = ifdOff + ifdLen;
  const total = extra + entries.reduce((s, e) => s + (SZ[e[1]] * e[2].length > inline ? SZ[e[1]] * e[2].length + 1 : 0), 0);
  const out = new Uint8Array(total + 2);
  const dv = new DataView(out.buffer);
  out.set(le ? [0x49, 0x49] : [0x4d, 0x4d]);
  dv.setUint16(2, big ? 43 : 42, le);
  if (big) { dv.setUint16(4, 8, le); dv.setBigUint64(8, BigInt(ifdOff), le); } else dv.setUint32(4, ifdOff, le);
  blocks.forEach((b, i) => out.set(b, offs[i]));
  const writeVals = (o, type, vals) => vals.forEach((v, i) => { if (type === 2) out[o + i] = v; else if (type === 3) dv.setUint16(o + 2 * i, v, le); else if (type === 4) dv.setUint32(o + 4 * i, v, le); else dv.setFloat64(o + 8 * i, v, le); });
  let p = ifdOff;
  if (big) { dv.setBigUint64(p, BigInt(entries.length), le); p += 8; } else { dv.setUint16(p, entries.length, le); p += 2; }
  for (const [tag, type, vals] of entries) {
    dv.setUint16(p, tag, le); dv.setUint16(p + 2, type, le);
    if (big) dv.setBigUint64(p + 4, BigInt(vals.length), le); else dv.setUint32(p + 4, vals.length, le);
    const vo = p + (big ? 12 : 8);
    if (SZ[type] * vals.length <= inline) writeVals(vo, type, vals);
    else { if (big) dv.setBigUint64(vo, BigInt(extra), le); else dv.setUint32(vo, extra, le); writeVals(extra, type, vals); extra += SZ[type] * vals.length; extra += extra % 2; }
    p += esz;
  }
  return out;
}

test('LZW and PackBits decoders invert the fixture encoders (table growth and reset included)', () => {
  let s = 12345;
  const rnd = () => { s = (Math.imul(s, 1664525) + 1013904223) >>> 0; return s; };
  const data = new Uint8Array(30000);
  for (let i = 0; i < data.length; i++) data[i] = i % 97 < 60 ? (rnd() >>> 24) & 7 : 200;
  assert.deepEqual(lzwDecode(lzwEncode(data), data.length), data);
  assert.deepEqual(packBitsDecode(packBitsEncode(data), data.length), data);
  assert.deepEqual([...packBitsDecode(Uint8Array.from([0xfe, 0xaa, 0x02, 0x80, 0x00, 0x2a, 0xfd, 0xaa]), 10)], [0xaa, 0xaa, 0xaa, 0x80, 0x00, 0x2a, 0xaa, 0xaa, 0xaa, 0xaa], 'TIFF 6.0 PackBits example prefix');
  assert.throws(() => lzwDecode(Uint8Array.from([0x80, 0x7f, 0xff]), 4), /Corrupt LZW/);
});

const W = 7, H = 5;
const ELEV = Float32Array.from({ length: W * H }, (_, i) => -20 - (i % W) * 1.5 - Math.floor(i / W) * 0.25);
const GEO = { scale: [2, 2, 0], tie: [0, 0, 0, 500000, 6000010, 0], keys: [[1024, 0, 1, 1], [1025, 0, 1, 1], [3072, 0, 1, 32631]] };

for (const [label, o] of Object.entries({
  'strips, uncompressed': { rowsPerStrip: 2 },
  'strips, LZW + floating-point predictor': { rowsPerStrip: 3, compression: 5, predictor: 3 },
  'tiles, Deflate, big-endian': { tile: [4, 4], compression: 8, le: false },
  'tiles, PackBits': { tile: [4, 4], compression: 32773 },
  'BigTIFF strips, LZW': { big: true, compression: 5 },
})) {
  test(`GeoTIFF float32: ${label}`, async () => {
    const r = await parseGeoTiff({ name: 'bathy.tif', bytes: makeTiff({ width: W, height: H, values: ELEV, geo: GEO, ...o }) });
    assertContract(r, ['grid']);
    const g = r.models[0];
    assert.deepEqual([g.nx, g.ny, g.dx, g.dy], [W, H, 2, 2]);
    // pixel-is-area: nodes at cell centres; upper-left corner (500000, 6000010)
    close(g.x0, 500001, 1e-9);
    close(g.y0, 6000010 - 2 * H + 1, 1e-9);
    for (let iy = 0; iy < H; iy++) for (let ix = 0; ix < W; ix++) close(g.z[iy * W + ix], ELEV[(H - 1 - iy) * W + ix], 1e-6, `${label} z(${ix},${iy})`);
    assert.equal(g.crs, 'EPSG:32631');
    assert.equal(r.meta.crs.kind, 'projected');
    assert.equal(r.meta.units.unit, 'm');
  });
}

test('GeoTIFF int16 with horizontal predictor, GDAL_NODATA, pixel-is-point and block averaging', async () => {
  const w = 8, h = 4;
  const v = Int16Array.from({ length: w * h }, (_, i) => (i === 9 ? -32768 : -100 - i));
  const geo = { scale: [5, 5, 0], tie: [0, 0, 0, 1000, 2000, 0], keys: [[1024, 0, 1, 1], [1025, 0, 1, 2], [3072, 0, 1, 25832]] };
  const bytes = makeTiff({ width: w, height: h, values: v, bits: 16, format: 2, compression: 5, predictor: 2, rowsPerStrip: 2, geo, nodata: -32768 });
  const r = await parseGeoTiff({ name: 'dem.tif', bytes });
  const g = r.models[0];
  assert.deepEqual([g.x0, g.y0], [1000, 2000 - 5 * (h - 1)], 'pixel-is-point: the tie point is the first node');
  assert.ok(Number.isNaN(g.z[(h - 1 - 1) * w + 1]));
  assert.equal(g.z[(h - 1) * w], -100);
  assert.equal(r.meta.skipped.noDataCells, 1);
  const d = await parseGeoTiff({ name: 'dem.tif', bytes }, { maxCells: 8 });
  assert.deepEqual([d.models[0].nx, d.models[0].ny, d.models[0].dx], [4, 2, 10]);
  // north-west 2 × 2 block: pixels 0, 1, 8 and (no-data) 9
  close(d.models[0].z[4], (-100 - 101 - 108) / 3, 1e-4);
  close(d.models[0].x0, 1000 - 2.5 + 5, 1e-9);
  assert.equal(d.meta.decimation.factor, 2);
  assert.deepEqual(d.meta.decimation.source, [w, h]);
});

test('GeoTIFF geographic CRS → metres; unsupported compression and non-TIFF input explain themselves', async () => {
  const geo = { scale: [0.01, 0.01, 0], tie: [0, 0, 0, 3.0, 54.05, 0], keys: [[1024, 0, 1, 2], [1025, 0, 1, 1], [2048, 0, 1, 4326]] };
  const r = await parseGeoTiff({ name: 'g.tif', bytes: makeTiff({ width: W, height: H, values: ELEV, geo }) });
  assert.equal(r.models[0].crs, 'EPSG:4326');
  close(r.models[0].dy, 1113.0, 0.6);
  assert.equal(r.meta.crs.projection, 'local-tangent-plane');
  const bad = makeTiff({ width: W, height: H, values: ELEV, geo: GEO });
  const dv = new DataView(bad.buffer);
  const ifd = dv.getUint32(4, true);
  for (let i = 0; i < dv.getUint16(ifd, true); i++) if (dv.getUint16(ifd + 2 + 12 * i, true) === 259) dv.setUint16(ifd + 2 + 12 * i + 8, 7, true);
  await assert.rejects(parseGeoTiff({ name: 'j.tif', bytes: bad }), /compression JPEG is not decoded here\. Re-write the file with LZW or Deflate/);
  await assert.rejects(parseGeoTiff({ name: 'x.tif', bytes: enc('this is not a tiff file at all') }), /expected the byte-order mark II or MM/);
  const rot = makeTiff({ width: W, height: H, values: ELEV, geo: { transform: [1, 0.5, 0, 0, 0.5, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1] } });
  await assert.rejects(parseGeoTiff({ name: 'r.tif', bytes: rot }), /rotated or sheared/);
});

// ---- GeoJSON / WKT ---------------------------------------------------------------------------------------

test('GeoJSON: points with properties, 3-D route, polygon with hole → section; lon/lat projected', () => {
  const fc = { type: 'FeatureCollection', features: [
    { type: 'Feature', properties: { name: 'WTG-01', rating: 15 }, geometry: { type: 'Point', coordinates: [3.0, 54.0] } },
    { type: 'Feature', properties: { name: 'WTG-02', rating: 15 }, geometry: { type: 'Point', coordinates: [3.02, 54.0, 5] } },
    { type: 'Feature', properties: { name: 'export cable' }, geometry: { type: 'LineString', coordinates: [[3.0, 54.0, -30], [3.0, 54.01, -32]] } },
    { type: 'Feature', properties: { name: 'lease' }, geometry: { type: 'Polygon', coordinates: [[[3, 54], [3.02, 54], [3.02, 54.01], [3, 54.01], [3, 54]], [[3.005, 54.002], [3.005, 54.004], [3.01, 54.004], [3.005, 54.002]]] } },
    { type: 'Feature', properties: null, geometry: null },
  ] };
  const r = parseGeoJson({ name: 'site.geojson', text: JSON.stringify(fc) }, { origin: { lon: 3, lat: 54 } });
  assertContract(r, ['points', 'polyline', 'section']);
  const pts = r.models.find((m) => m.kind === 'points');
  assert.deepEqual(pts.attrs.names, ['WTG-01', 'WTG-02']);
  assert.deepEqual([...pts.attrs.rating], [15, 15]);
  close(pts.positions[3], 1311.5, 0.05);
  close(pts.positions[5], 5, 1e-9);
  const line = r.models.find((m) => m.kind === 'polyline');
  assert.equal(line.name, 'export cable');
  close(line.points[4], 1113.0, 0.5);
  close(line.points[5], -32, 1e-12);
  const sec = r.models.find((m) => m.kind === 'section');
  assert.deepEqual(sec.loops.map((l) => l.length / 2), [4, 3], 'closing vertices removed');
  assert.equal(sec.properties.name, 'lease');
  assert.equal(r.meta.crs.source, 'EPSG:4326');
  assert.equal(r.meta.skipped.featureWithoutGeometry, 1);
  assert.deepEqual(r.meta.features, { points: 2, lines: 1, polygons: 1, polygonsAsSections: 1, polygonsAsPolylines: 0, meshes: 0 });
});

test('GeoJSON with a projected crs member, newline-delimited features, UTM projection option', () => {
  const g = { type: 'FeatureCollection', crs: { type: 'name', properties: { name: 'urn:ogc:def:crs:EPSG::32631' } }, features: [{ type: 'Feature', properties: {}, geometry: { type: 'LineString', coordinates: [[500000, 6000000, -30], [500100, 6000000, -31]] } }] };
  const r = parseGeoJson({ name: 'r.geojson', text: JSON.stringify(g) });
  closeArr([...r.models[0].points], [500000, 6000000, -30, 500100, 6000000, -31]);
  assert.equal(r.meta.crs.epsg, 32631);
  const nd = parseGeoJson({ name: 'p.geojsonl', text: '{"type":"Feature","properties":{"name":"a"},"geometry":{"type":"Point","coordinates":[3,54]}}\n{"type":"Feature","properties":{"name":"b"},"geometry":{"type":"Point","coordinates":[3.1,54]}}\n' }, { projection: 'utm' });
  const u = utmForward(54, 3);
  const m = nd.models[0];
  close(m.positions[0] + m.origin[0], u.easting, 0.05);
  close(m.positions[1] + m.origin[1], u.northing, 0.6);
  assert.equal(nd.meta.crs.zone, 31);
  assert.throws(() => parseGeoJson({ name: 'x.json', text: '{"type":"Topology","objects":{}}' }), /TopoJSON is not read/);
  assert.throws(() => parseGeoJson({ name: 'x.json', text: '{"a":1}' }), /not GeoJSON/);
  assert.equal(crsFromName('http://www.opengis.net/def/crs/EPSG/0/4326').latFirst, true);
  assert.equal(crsFromName('EPSG:4326').latFirst, false);
});

test('WKT: Z line, polygon with hole, multipoint forms, TIN, EWKT SRID', () => {
  const text = 'SRID=32631;LINESTRING Z (500000 6000000 -30, 500100 6000000 -31)\nPOLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 2 4, 4 4, 2 2))\nMULTIPOINT ((1 2), (3 4))\nMULTIPOINT Z (5 6 7, 8 9 10)\nTIN Z (((0 0 0, 1 0 0, 0 1 1, 0 0 0)), ((1 0 0, 1 1 1, 0 1 1, 1 0 0)))\nPOINT M (1 1 99)\n';
  const r = parseWkt({ name: 'geom.wkt', text });
  assertContract(r, ['points', 'polyline', 'section', 'surface']);
  assert.equal(r.models.find((m) => m.kind === 'polyline').points[5], -31);
  assert.deepEqual(r.models.find((m) => m.kind === 'section').loops.map((l) => l.length / 2), [4, 3]);
  const pts = r.models.find((m) => m.kind === 'points');
  assert.equal(pts.positions.length / 3, 5);
  assert.equal(pts.positions[14], 0, 'M ordinate is not an elevation');
  close(surfaceArea(r.models.find((m) => m.kind === 'surface')), Math.SQRT2, 1e-6);
  assert.equal(r.meta.crs.epsg, 32631);
  assert.throws(() => parseWkt({ name: 'x.wkt', text: 'hello there\n' }), /no WKT geometry found/);
});

// ---- Shapefile -------------------------------------------------------------------------------------------

function makeShp(type, shapes) {
  const recs = shapes.map((s) => {
    let b;
    if (type === 1 || type === 11) { b = new Uint8Array(type === 11 ? 36 : 20); const dv = new DataView(b.buffer); dv.setInt32(0, type, true); dv.setFloat64(4, s[0], true); dv.setFloat64(12, s[1], true); if (type === 11) dv.setFloat64(20, s[2], true); } else {
      const parts = s, n = parts.reduce((a, p) => a + p.length, 0), z = type > 10;
      b = new Uint8Array(44 + 4 * parts.length + 16 * n + (z ? 16 + 8 * n : 0));
      const dv = new DataView(b.buffer);
      dv.setInt32(0, type, true); dv.setInt32(36, parts.length, true); dv.setInt32(40, n, true);
      let o = 44, k = 0;
      for (const p of parts) { dv.setInt32(o, k, true); o += 4; k += p.length; }
      for (const p of parts) for (const q of p) { dv.setFloat64(o, q[0], true); dv.setFloat64(o + 8, q[1], true); o += 16; }
      if (z) { o += 16; for (const p of parts) for (const q of p) { dv.setFloat64(o, q[2], true); o += 8; } }
    }
    return b;
  });
  const total = 100 + recs.reduce((a, r) => a + 8 + r.length, 0);
  const out = new Uint8Array(total);
  const dv = new DataView(out.buffer);
  dv.setInt32(0, 9994, false); dv.setInt32(24, total / 2, false); dv.setInt32(28, 1000, true); dv.setInt32(32, type, true);
  let p = 100;
  recs.forEach((r, i) => { dv.setInt32(p, i + 1, false); dv.setInt32(p + 4, r.length / 2, false); out.set(r, p + 8); p += 8 + r.length; });
  return out;
}

function makeDbf(fields, rows) {
  const recSize = 1 + fields.reduce((a, f) => a + f.len, 0), headerSize = 32 + 32 * fields.length + 1;
  const out = new Uint8Array(headerSize + recSize * rows.length + 1).fill(0x20);
  const dv = new DataView(out.buffer);
  out.fill(0, 0, headerSize);
  out[0] = 3; dv.setUint32(4, rows.length, true); dv.setUint16(8, headerSize, true); dv.setUint16(10, recSize, true);
  fields.forEach((f, i) => { out.set(enc(f.name), 32 + 32 * i); out[32 + 32 * i + 11] = f.type.charCodeAt(0); out[32 + 32 * i + 16] = f.len; });
  out[headerSize - 1] = 0x0d;
  rows.forEach((r, i) => { let o = headerSize + i * recSize + 1; fields.forEach((f, k) => { const s = String(r[k]); out.set(enc(f.type === 'N' ? s.padStart(f.len) : s.padEnd(f.len)).subarray(0, f.len), o); o += f.len; }); });
  out[out.length - 1] = 0x1a;
  return out;
}

const PRJ_UTM31 = 'PROJCS["WGS_1984_UTM_Zone_31N",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",3.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]';

test('Shapefile: PolyLineZ with .dbf attributes and UTM .prj', () => {
  const files = new Map([
    ['routes.shp', makeShp(13, [[[[500000, 6000000, -30], [500100, 6000000, -31], [500200, 6000050, -33]]], [[[1, 1, 0], [2, 2, 0]], [[5, 5, -1], [6, 6, -2]]]])],
    ['routes.dbf', makeDbf([{ name: 'NAME', type: 'C', len: 12 }, { name: 'OD_MM', type: 'N', len: 8 }], [['export A', 273.1], ['spur', 114.3]])],
    ['routes.prj', enc(PRJ_UTM31)],
  ]);
  const r = parseShapefile({ name: 'routes.shp', files });
  assertContract(r, ['polyline']);
  assert.deepEqual(r.models.map((m) => [m.name, m.points.length / 3]), [['export A', 3], ['spur', 2], ['spur', 2]]);
  closeArr([...r.models[0].points], [500000, 6000000, -30, 500100, 6000000, -31, 500200, 6000050, -33]);
  assert.deepEqual(r.models[0].properties, { NAME: 'export A', OD_MM: 273.1 });
  assert.equal(r.meta.crs.source, 'EPSG:32631');
  assert.deepEqual(r.meta.attributeFields, ['NAME', 'OD_MM']);
  assert.equal(r.meta.shapeType, 13);
});

test('Shapefile: polygons with holes (ring direction), points, bare .shp, malformed input', () => {
  // outer rings clockwise, hole counter-clockwise
  const outer = [[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]], hole = [[2, 2], [4, 2], [4, 4], [2, 4], [2, 2]], other = [[20, 0], [20, 5], [25, 5], [25, 0], [20, 0]];
  const r = parseShapefile({ name: 'areas.shp', bytes: makeShp(5, [[outer, hole, other]]) });
  assertContract(r, ['section']);
  assert.deepEqual(r.models.map((m) => m.loops.map((l) => l.length / 2)), [[4, 4], [4]]);
  assert.ok(r.warnings.some((w) => /No \.dbf/.test(w)) && r.warnings.some((w) => /No \.prj/.test(w)));
  const p = parseShapefile({ name: 'pts.shp', bytes: makeShp(11, [[1, 2, 3], [4, 5, 6]]) });
  closeArr([...p.models[0].positions], [1, 2, 3, 4, 5, 6]);
  assert.throws(() => parseShapefile({ name: 'x.shp', bytes: new Uint8Array(120) }), /expected file code 9994/);
  assert.throws(() => parseShapefile({ name: 'x.shp', files: new Map([['a.dbf', new Uint8Array(4)]]) }), /no \.shp file was supplied/);
});

// ---- LandXML / GML / KML / GPX ---------------------------------------------------------------------------

test('LandXML: TIN (northing-easting order, invisible faces), CgPoints, alignment with arc and profile', () => {
  const xml = `<?xml version="1.0"?><LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2" version="1.2">
  <Units><Imperial linearUnit="foot" areaUnit="squareFoot"/></Units>
  <CoordinateSystem epsgCode="32631" name="UTM31"/>
  <Application name="survey"/>
  <CgPoints><CgPoint name="BM1" code="BM">100 200 10</CgPoint><CgPoint name="BM2">110 200 12</CgPoint></CgPoints>
  <Surfaces><Surface name="seabed"><Definition surfType="TIN">
    <Pnts><P id="1">0 0 -10</P><P id="2">0 100 -12</P><P id="3">100 100 -14</P><P id="4">100 0 -11</P></Pnts>
    <Faces><F>1 2 3</F><F>1 4 3</F><F i="1">1 2 4</F></Faces>
  </Definition></Surface></Surfaces>
  <Alignments><Alignment name="CL1" length="257.08" staStart="1000">
    <CoordGeom>
      <Line><Start>0 0</Start><End>0 100</End></Line>
      <Curve rot="ccw" radius="100"><Start>0 100</Start><Center>100 100</Center><End>100 200</End></Curve>
    </CoordGeom>
    <Profile><ProfAlign name="p"><PVI>1000 0</PVI><ParaCurve length="100">1100 10</ParaCurve><PVI>1300 10</PVI></ProfAlign></Profile>
  </Alignment></Alignments>
  <Parcels><Parcel name="lot"/></Parcels>
</LandXML>`;
  const ft = 0.3048;
  const r = parseLandXml({ name: 'site.xml', bytes: enc(xml) });
  assertContract(r, ['points', 'surface', 'polyline']);
  assert.equal(r.meta.units.unit, 'ft');
  const pts = r.models.find((m) => m.kind === 'points');
  closeArr([...pts.positions].slice(0, 3), [200 * ft, 100 * ft, 10 * ft], 1e-5);
  assert.deepEqual(pts.attrs.names, ['BM1', 'BM2']);
  const s = r.models.find((m) => m.kind === 'surface');
  assert.equal(s.indices.length, 6);
  close(s.positions[3], 100 * ft, 1e-5, 'second point: easting 100');
  // the two faces are wound oppositely in the file; both come out with the normal up
  const p = s.positions, t = s.indices;
  for (let f = 0; f < 6; f += 3) assert.ok((p[3 * t[f + 1]] - p[3 * t[f]]) * (p[3 * t[f + 2] + 1] - p[3 * t[f] + 1]) - (p[3 * t[f + 1] + 1] - p[3 * t[f] + 1]) * (p[3 * t[f + 2]] - p[3 * t[f]]) > 0);
  const al = r.models.find((m) => m.kind === 'polyline');
  const n = al.points.length / 3;
  close(al.attrs.kp[0], 1000 * ft, 1e-9);
  // 100 ft tangent + quarter circle of radius 100 ft
  close(al.attrs.kp[n - 1], (1000 + 100 + (Math.PI / 2) * 100) * ft, 1e-6);
  closeArr([al.points[3 * n - 3], al.points[3 * n - 2]], [200 * ft, 100 * ft], 1e-9);
  // parabolic crest curve: at the PVI station the curve sits (g1 − g2)·L/8 below the PVI
  const pv = [{ s: 1000, z: 0, L: 0 }, { s: 1100, z: 10, L: 100 }, { s: 1300, z: 10, L: 0 }];
  close(profileElevation(pv, 1100), 10 - (0.1 * 100) / 8, 1e-12);
  close(profileElevation(pv, 1050), 5, 1e-12);
  close(profileElevation(pv, 1200), 10, 1e-12);
  const iPvi = [...al.attrs.kp].findIndex((v) => Math.abs(v - 1100 * ft) < 1e-9);
  assert.ok(iPvi > 0, 'a vertex is inserted at the PVI station');
  close(al.points[3 * iPvi + 2], 8.75 * ft, 1e-9);
  assert.equal(r.meta.crs.source, 'EPSG:32631');
  assert.equal(r.meta.skipped['invisible TIN faces'], 1);
  assert.equal(r.meta.skipped['<Parcels>'], 1);
  assert.throws(() => parseLandXml({ name: 'x.xml', text: '<Other/>' }), /expected root element <LandXML>/);
  assert.throws(() => parseLandXml({ name: 'x.xml', text: '<LandXML><Units/></LandXML>' }), /no geometry was found/);
});

test('LandXML pipe network → frame with tube sections', () => {
  const xml = `<LandXML><Units><Metric linearUnit="meter"/></Units><PipeNetworks><PipeNetwork name="outfall">
   <Structs><Struct name="MH1" elevSump="-2"><Center>0 0</Center><Invert elev="-1.5" refPipe="P1"/></Struct><Struct name="MH2" elevSump="-4"><Center>0 50</Center></Struct></Structs>
   <Pipes><Pipe name="P1" refStart="MH1" refEnd="MH2"><CircPipe diameter="0.9" thickness="0.05"/></Pipe></Pipes></PipeNetwork></PipeNetworks></LandXML>`;
  const r = parseLandXml({ name: 'net.xml', text: xml });
  const f = r.models[0];
  assert.equal(f.kind, 'frame');
  closeArr([...f.nodes], [0, 0, -1.5, 50, 0, -4]);
  assert.deepEqual(f.members[0].section, { type: 'tube', D: 1.0, t: 0.05 });
});

test('GML: latitude-first URN CRS, posList dimension, polygon with hole, CityGML solid → surface', () => {
  const gml = `<?xml version="1.0"?><wfs:FeatureCollection xmlns:wfs="http://www.opengis.net/wfs" xmlns:gml="http://www.opengis.net/gml" xmlns:app="urn:x">
  <gml:boundedBy><gml:Envelope srsName="urn:ogc:def:crs:EPSG::4326"><gml:lowerCorner>54 3</gml:lowerCorner><gml:upperCorner>54.01 3.02</gml:upperCorner></gml:Envelope></gml:boundedBy>
  <gml:featureMember><app:Cable><gml:name>IAC-1</gml:name><app:geom><gml:LineString srsName="urn:ogc:def:crs:EPSG::4326" srsDimension="3"><gml:posList>54 3 -30 54.01 3 -32</gml:posList></gml:LineString></app:geom></app:Cable></gml:featureMember>
  <gml:featureMember><app:Zone><app:geom><gml:Polygon srsName="urn:ogc:def:crs:EPSG::4326"><gml:exterior><gml:LinearRing><gml:posList>54 3 54 3.02 54.01 3.02 54.01 3 54 3</gml:posList></gml:LinearRing></gml:exterior></gml:Polygon></app:geom></app:Zone></gml:featureMember>
  <gml:featureMember><app:Mast><app:geom><gml:Point srsName="urn:ogc:def:crs:EPSG::4326"><gml:pos>54 3.02</gml:pos></gml:Point></app:geom></app:Mast></gml:featureMember>
</wfs:FeatureCollection>`;
  const r = parseGml({ name: 'site.gml', text: gml }, { origin: { lon: 3, lat: 54 } });
  assertContract(r, ['points', 'polyline', 'section']);
  const line = r.models.find((m) => m.kind === 'polyline');
  assert.equal(line.name, 'IAC-1');
  close(line.points[0], 0, 1e-6);
  close(line.points[4], 1113.0, 0.5, 'latitude was the first ordinate');
  close(line.points[5], -32, 1e-12);
  close(r.models.find((m) => m.kind === 'points').positions[0], 1311.5, 0.05);
  assert.equal(r.models.find((m) => m.kind === 'section').loops[0].length, 8);
  assert.match(r.meta.axisOrder, /latitude\/northing first/);
  assert.equal(r.meta.features.points, 1, 'the envelope corners are not features');

  const ring = (pts) => `<gml:surfaceMember><gml:Polygon><gml:exterior><gml:LinearRing><gml:posList>${[...pts, pts[0]].map((p) => p.join(' ')).join(' ')}</gml:posList></gml:LinearRing></gml:exterior></gml:Polygon></gml:surfaceMember>`;
  const V = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [0, 0, 5], [10, 0, 5], [10, 10, 5], [0, 10, 5]];
  const faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
  const city = `<core:CityModel xmlns:core="http://www.opengis.net/citygml/2.0" xmlns:bldg="http://www.opengis.net/citygml/building/2.0" xmlns:gml="http://www.opengis.net/gml"><core:cityObjectMember><bldg:Building gml:id="B1"><bldg:lod1Solid><gml:Solid srsName="EPSG:25832" srsDimension="3"><gml:exterior><gml:CompositeSurface>${faces.map((f) => ring(f.map((k) => V[k]))).join('')}</gml:CompositeSurface></gml:exterior></gml:Solid></bldg:lod1Solid></bldg:Building></core:cityObjectMember></core:CityModel>`;
  const c = parseGml({ name: 'city.gml', text: city });
  assertContract(c, ['surface']);
  close(surfaceArea(c.models[0]), 2 * 100 + 4 * 50, 1e-6);
  assert.equal(c.models[0].name, 'Building');
  assert.equal(c.meta.crs.epsg, 25832);
  assert.throws(() => parseGml({ name: 'x.gml', text: '<a><b/></a>' }), /no geometry was found/);
});

const KML = `<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>
  <Placemark><name>OSS</name><ExtendedData><Data name="capacity"><value>800</value></Data></ExtendedData><Point><coordinates>3.0,54.0,25</coordinates></Point></Placemark>
  <Folder><Placemark><name>route</name><LineString><coordinates>
     3.0,54.0,0 3.0,54.01,0
     3.02,54.01
  </coordinates></LineString></Placemark></Folder>
  <Placemark><name>area</name><MultiGeometry><Polygon><outerBoundaryIs><LinearRing><coordinates>3,54 3.02,54 3.02,54.01 3,54</coordinates></LinearRing></outerBoundaryIs></Polygon></MultiGeometry></Placemark>
</Document></kml>`;

test('KML and KMZ: placemarks, extended data, optional altitude', async () => {
  const r = await parseKml({ name: 'site.kml', text: KML }, { origin: { lon: 3, lat: 54 } });
  assertContract(r, ['points', 'polyline', 'section']);
  const pts = r.models.find((m) => m.kind === 'points');
  assert.deepEqual(pts.attrs.names, ['OSS']);
  assert.deepEqual([...pts.attrs.capacity], [800]);
  close(pts.positions[2], 25, 1e-9);
  const line = r.models.find((m) => m.kind === 'polyline');
  assert.equal(line.points.length, 9);
  close(line.points[6], 1311.5, 0.6);
  assert.equal(r.models.find((m) => m.kind === 'section').loops[0].length, 6);
  const z = await parseKml({ name: 'site.kmz', bytes: makeZip([{ name: 'doc.kml', data: enc(KML) }, { name: 'files/icon.png', data: new Uint8Array(4) }]) });
  assert.equal(z.meta.counts.models, 3);
  await assert.rejects(parseKml({ name: 'x.kml', text: '<html/>' }), /not a KML document/);
});

test('GPX: waypoints, route and track segments', () => {
  const gpx = `<gpx version="1.1" creator="unit test" xmlns="http://www.topografix.com/GPX/1/1"><wpt lat="54.0" lon="3.0"><ele>12</ele><name>WP1</name><sym>Buoy</sym></wpt>
  <rte><name>R</name><rtept lat="54.0" lon="3.0"/><rtept lat="54.01" lon="3.0"/></rte>
  <trk><name>survey</name><trkseg><trkpt lat="54.0" lon="3.0"><ele>1</ele></trkpt><trkpt lat="54.0" lon="3.02"><ele>2</ele></trkpt></trkseg><trkseg><trkpt lat="54.01" lon="3.0"/><trkpt lat="54.01" lon="3.02"/></trkseg></trk></gpx>`;
  const r = parseGpx({ name: 'track.gpx', text: gpx }, { origin: { lon: 3, lat: 54 } });
  assertContract(r, ['points', 'polyline']);
  assert.deepEqual(r.models.map((m) => m.name), ['track (points)', 'R', 'survey (1)', 'survey (2)']);
  close(r.models[2].points[3], 1311.5, 0.05);
  close(r.models[2].points[5], 2, 1e-12);
  assert.deepEqual(r.models[0].attrs.sym, ['Buoy']);
  assert.equal(r.meta.creator, 'unit test');
  assert.throws(() => parseGpx({ name: 'x.gpx', text: '<kml/>' }), /expected root element <gpx>/);
});

// ---- NetCDF ----------------------------------------------------------------------------------------------

/** Minimal NetCDF classic (CDF-1) writer for fixtures: fixed-size variables only. */
function makeNetcdf({ dims, gatts = {}, vars }) {
  const parts = [];
  const u32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, v, false); parts.push(b); };
  const pad = (n) => { if (n % 4) parts.push(new Uint8Array(4 - (n % 4))); };
  const str = (s) => { const b = enc(s); u32(b.length); parts.push(b); pad(b.length); };
  const TYPES = { char: 2, short: 3, int: 4, float: 5, double: 6 };
  const SET = { 3: ['setInt16', 2], 4: ['setInt32', 4], 5: ['setFloat32', 4], 6: ['setFloat64', 8] };
  const encode = (type, vals) => {
    if (type === 2) return typeof vals === 'string' ? enc(vals) : Uint8Array.from(vals);
    const [fn, sz] = SET[type];
    const b = new Uint8Array(vals.length * sz), dv = new DataView(b.buffer);
    vals.forEach((v, i) => dv[fn](i * sz, v, false));
    return b;
  };
  const atts = (o) => {
    const k = Object.keys(o);
    if (!k.length) { u32(0); u32(0); return; }
    u32(0x0c); u32(k.length);
    for (const nm of k) {
      const v = o[nm];
      str(nm);
      const type = typeof v === 'string' ? 2 : v.type ? TYPES[v.type] : 6;
      const vals = typeof v === 'string' ? v : v.type ? v.values : [].concat(v);
      const b = encode(type, vals);
      u32(type); u32(typeof vals === 'string' ? b.length : vals.length); parts.push(b); pad(b.length);
    }
  };
  parts.push(Uint8Array.from([0x43, 0x44, 0x46, 1]));
  u32(0);
  const dn = Object.keys(dims);
  u32(0x0a); u32(dn.length);
  for (const d of dn) { str(d); u32(dims[d]); }
  atts(gatts);
  u32(0x0b); u32(vars.length);
  const slots = [];
  const datas = vars.map((v) => { const b = encode(TYPES[v.type], v.data); return b.length % 4 ? concat(b, new Uint8Array(4 - (b.length % 4))) : b; });
  vars.forEach((v, i) => {
    str(v.name); u32(v.dims.length);
    for (const d of v.dims) u32(dn.indexOf(d));
    atts(v.atts || {});
    u32(TYPES[v.type]); u32(datas[i].length);
    const slot = new Uint8Array(4); parts.push(slot); slots.push(slot);
  });
  let off = parts.reduce((s, p) => s + p.length, 0);
  datas.forEach((d, i) => { new DataView(slots[i].buffer).setUint32(0, off, false); off += d.length; });
  return concat(...parts, ...datas);
}

test('NetCDF classic CF grid: packed shorts, fill value, north-to-south latitude, positive down', () => {
  const lat = [54.02, 54.01, 54.0], lon = [3.0, 3.01];
  const bytes = makeNetcdf({
    dims: { time: 1, lat: 3, lon: 2 }, gatts: { Conventions: 'CF-1.8' },
    vars: [
      { name: 'lat', dims: ['lat'], type: 'double', data: lat, atts: { units: 'degrees_north', standard_name: 'latitude' } },
      { name: 'lon', dims: ['lon'], type: 'double', data: lon, atts: { units: 'degrees_east' } },
      { name: 'deptho', dims: ['time', 'lat', 'lon'], type: 'short', data: [300, 310, 320, -32768, 340, 350], atts: { scale_factor: { type: 'double', values: [0.1] }, add_offset: { type: 'double', values: [1.0] }, _FillValue: { type: 'short', values: [-32768] }, positive: 'down', units: 'm', standard_name: 'sea_floor_depth_below_sea_level' } },
      { name: 'mask', dims: ['lat', 'lon'], type: 'short', data: [1, 1, 1, 0, 1, 1] },
    ],
  });
  const r = parseNetcdf({ name: 'gebco.nc', bytes }, { origin: { lon: 3, lat: 54 } });
  assertContract(r, ['grid']);
  const g = r.models[0];
  assert.deepEqual([g.nx, g.ny], [2, 3]);
  close(g.x0, 0, 1e-6); close(g.y0, 0, 1e-6);
  close(g.dy, 1113.0, 0.5);
  // southern row first; depth 34 + 1 → z = −35
  closeArr([...g.z].map((v) => (Number.isNaN(v) ? 0 : v)), [-35, -36, -33, 0, -31, -32], 1e-5);
  assert.equal(r.meta.variable, 'deptho');
  assert.equal(r.meta.vertical.source, 'depth-positive-down');
  assert.equal(r.meta.crs.kind, 'geographic');
  assert.ok(r.warnings.some((w) => /Several gridded variables/.test(w)));
  const m = parseNetcdf({ name: 'gebco.nc', bytes }, { variable: 'mask', projection: 'none' });
  assert.deepEqual([...m.models[0].z], [1, 1, 1, 0, 1, 1], 'rows reversed: the file runs north to south');
  close(m.models[0].dx, 0.01, 1e-9);
  assert.throws(() => parseNetcdf({ name: 'gebco.nc', bytes }, { variable: 'nope' }), /Gridded variables: deptho, mask/);
  assert.throws(() => parseNetcdf({ name: 'x.nc', bytes: concat(Uint8Array.from([0x89, 0x48, 0x44, 0x46, 13, 10, 26, 10]), new Uint8Array(40)) }), /NetCDF-4 \/ HDF5 file.*nccopy -k classic/);
  assert.throws(() => parseNetcdf({ name: 'x.nc', bytes: new Uint8Array(64) }), /expected the magic bytes "CDF"/);
});

test('NetCDF: projected grid in km, UGRID mesh, Exodus II mesh on the classic format', () => {
  const g = parseNetcdf({ name: 'local.nc', bytes: makeNetcdf({ dims: { y: 2, x: 3 }, vars: [
    { name: 'x', dims: ['x'], type: 'double', data: [500, 501, 502], atts: { units: 'km', axis: 'X' } },
    { name: 'y', dims: ['y'], type: 'double', data: [6000, 6001], atts: { units: 'km', axis: 'Y' } },
    { name: 'elevation', dims: ['y', 'x'], type: 'float', data: [-1, -2, -3, -4, -5, -6], atts: { units: 'm' } },
  ] }) });
  assert.deepEqual([g.models[0].x0, g.models[0].dx, g.models[0].y0, ...g.models[0].z], [500000, 1000, 6000000, -1, -2, -3, -4, -5, -6]);
  const u = parseNetcdf({ name: 'mesh.nc', bytes: makeNetcdf({ dims: { node: 5, face: 2, max: 4 }, gatts: { Conventions: 'UGRID-1.0' }, vars: [
    { name: 'mesh2d', dims: [], type: 'int', data: [0], atts: { cf_role: 'mesh_topology', topology_dimension: { type: 'int', values: [2] }, node_coordinates: 'node_x node_y', face_node_connectivity: 'face_nodes' } },
    { name: 'node_x', dims: ['node'], type: 'double', data: [0, 10, 10, 0, 20], atts: { units: 'm' } },
    { name: 'node_y', dims: ['node'], type: 'double', data: [0, 0, 10, 10, 5], atts: { units: 'm' } },
    { name: 'node_z', dims: ['node'], type: 'double', data: [1, 2, 3, 4, 5], atts: { standard_name: 'sea_floor_depth_below_geoid', positive: 'down' } },
    { name: 'face_nodes', dims: ['face', 'max'], type: 'int', data: [1, 2, 3, 4, 2, 5, 3, -999], atts: { start_index: { type: 'int', values: [1] }, _FillValue: { type: 'int', values: [-999] } } },
  ] }) });
  assertContract(u, ['surface']);
  assert.deepEqual(u.meta.elementTypes, { quad4: 1, tri3: 1 });
  assert.equal(u.models[0].positions[14], -5);
  const HEXN = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
  const names = (list, w) => list.flatMap((s) => [...enc(s.padEnd(w, '\0'))]);
  const e = parseNetcdf({ name: 'block.e', bytes: makeNetcdf({ dims: { len_name: 8, num_dim: 3, num_nodes: 8, num_elem: 1, num_el_blk: 1, num_el_in_blk1: 1, num_nod_per_el1: 8, num_node_sets: 1, num_nod_ns1: 4 }, gatts: { title: 'cube' }, vars: [
    { name: 'coordx', dims: ['num_nodes'], type: 'double', data: HEXN.map((p) => p[0]) }, { name: 'coordy', dims: ['num_nodes'], type: 'double', data: HEXN.map((p) => p[1]) }, { name: 'coordz', dims: ['num_nodes'], type: 'double', data: HEXN.map((p) => p[2]) },
    { name: 'eb_prop1', dims: ['num_el_blk'], type: 'int', data: [100] }, { name: 'eb_names', dims: ['num_el_blk', 'len_name'], type: 'char', data: names(['grout'], 8) },
    { name: 'connect1', dims: ['num_el_in_blk1', 'num_nod_per_el1'], type: 'int', data: [1, 2, 3, 4, 5, 6, 7, 8], atts: { elem_type: 'HEX8' } },
    { name: 'ns_prop1', dims: ['num_node_sets'], type: 'int', data: [7] }, { name: 'node_ns1', dims: ['num_nod_ns1'], type: 'int', data: [1, 2, 3, 4] },
  ] }) });
  assertContract(e, ['volume']);
  checkElements(e.models[0]);
  assert.equal(e.meta.format, 'exodus');
  assert.deepEqual(e.models[0].tagNames, { 100: 'grout' });
  assert.deepEqual([...e.models[0].sets.nodes['nodeset 7']], [0, 1, 2, 3]);
});

// ---- point clouds ----------------------------------------------------------------------------------------

function makeLas(points, { format = 2, scale = 0.01, offset = [500000, 6000000, 0], vlrs = [] } = {}) {
  const recLen = { 0: 20, 1: 28, 2: 26, 3: 34, 6: 30, 7: 36 }[format];
  const vlrBytes = vlrs.reduce((s, v) => s + 54 + v.data.length, 0);
  const dataOff = 227 + vlrBytes;
  const out = new Uint8Array(dataOff + recLen * points.length);
  const dv = new DataView(out.buffer);
  out.set(enc('LASF')); out[24] = 1; out[25] = 2; out.set(enc('unit test'), 26);
  dv.setUint16(94, 227, true); dv.setUint32(96, dataOff, true); dv.setUint32(100, vlrs.length, true); out[104] = format; dv.setUint16(105, recLen, true); dv.setUint32(107, points.length, true);
  [scale, scale, scale, ...offset].forEach((v, i) => dv.setFloat64(131 + 8 * i, v, true));
  let p = 227;
  for (const v of vlrs) { out.set(enc(v.user), p + 2); dv.setUint16(p + 18, v.id, true); dv.setUint16(p + 20, v.data.length, true); out.set(v.data, p + 54); p += 54 + v.data.length; }
  points.forEach((pt, i) => {
    const o = dataOff + i * recLen;
    for (let k = 0; k < 3; k++) dv.setInt32(o + 4 * k, Math.round((pt[k] - offset[k]) / scale), true);
    dv.setUint16(o + 12, pt[3] ?? 0, true);
    out[o + (format >= 6 ? 16 : 15)] = pt[4] ?? 0;
    const ro = { 2: 20, 3: 28, 7: 30 }[format];
    if (ro) for (let k = 0; k < 3; k++) dv.setUint16(o + ro + 2 * k, (pt[5 + k] ?? 0) * 257, true);
  });
  return out;
}

test('LAS: scaled integer coordinates, intensity, classification, 16-bit colour, GeoKey CRS, thinning', () => {
  const pts = [[500000.12, 6000000.5, -30.25, 100, 2, 255, 0, 10], [500010.0, 6000005.0, -31.0, 200, 9, 0, 128, 20], [500020.5, 6000009.75, -29.5, 300, 2, 1, 2, 3], [500030, 6000010, -28, 400, 2, 4, 5, 6]];
  const keys = new Uint8Array(new Uint16Array([1, 1, 0, 2, 1024, 0, 1, 1, 3072, 0, 1, 32631]).buffer);
  const r = parseLas({ name: 'mbes.las', bytes: makeLas(pts, { vlrs: [{ user: 'LASF_Projection', id: 34735, data: keys }] }) });
  assertContract(r, ['points']);
  const m = r.models[0];
  for (let i = 0; i < 4; i++) for (let k = 0; k < 3; k++) close(m.positions[3 * i + k] + m.origin[k], pts[i][k], 1e-3);
  assert.deepEqual([...m.attrs.intensity], [100, 200, 300, 400]);
  assert.deepEqual([...m.attrs.classification], [2, 9, 2, 2]);
  assert.deepEqual([...m.attrs.rgb.slice(0, 6)], [255, 0, 10, 0, 128, 20]);
  assert.equal(r.meta.crs.source, 'EPSG:32631');
  assert.deepEqual(r.meta.classes, { 2: 3, 9: 1 });
  assert.equal(r.meta.pointFormat, 2);
  const t = parseLas({ name: 'mbes.las', bytes: makeLas(pts) }, { maxPoints: 2 });
  assert.equal(t.models[0].positions.length / 3, 2);
  assert.equal(t.meta.thinning, 2);
  assert.ok(t.warnings.some((w) => /every 2nd point was kept/.test(w)));
  const f7 = parseLas({ name: 'f7.las', bytes: makeLas(pts, { format: 7 }) });
  assert.deepEqual([...f7.models[0].attrs.classification], [2, 9, 2, 2]);
  assert.deepEqual([...f7.models[0].attrs.rgb.slice(3, 6)], [0, 128, 20]);
  const laz = makeLas(pts); laz[104] |= 0x80;
  assert.throws(() => parseLas({ name: 'c.laz', bytes: laz }), /LASzip-compressed \(LAZ\)\. Decompress first/);
  assert.throws(() => parseLas({ name: 'x.las', bytes: new Uint8Array(300) }), /expected the signature "LASF"/);
});

test('PTS, PTX (registration matrix) and XYZRGB text clouds', () => {
  const pts = parsePts({ name: 'scan.pts', text: '2\n1.0 2.0 3.0 -1200 255 128 0\n4.0 5.0 6.0 -800 10 20 30\n' });
  assertContract(pts, ['points']);
  assert.deepEqual([...pts.models[0].positions], [1, 2, 3, 4, 5, 6]);
  assert.deepEqual([...pts.models[0].attrs.rgb], [255, 128, 0, 10, 20, 30]);
  assert.deepEqual([...pts.models[0].attrs.intensity], [-1200, -800]);
  // one scan, rotated 90° about z and shifted by (10, 0, 0): row-vector convention
  const ptx = parsePtx({ name: 'scan.ptx', text: '2\n1\n10 0 0\n0 1 0\n-1 0 0\n0 0 1\n0 1 0 0\n-1 0 0 0\n0 0 1 0\n10 0 0 1\n1 0 0 0.5\n0 0 0 0.5\n' });
  closeArr([...ptx.models[0].positions], [10, 1, 0], 1e-9);
  assert.equal(ptx.meta.skipped.emptyReturns, 1);
  assert.equal(ptx.meta.scans, 1);
  const x = parseXyzRgb({ name: 'cloud.xyzrgb', text: '0 0 0 0.5 0.25 1.0\n1 0 0 1.0 0.0 0.0\nbad line\n' });
  assert.deepEqual([...x.models[0].attrs.rgb], [127, 63, 255, 255, 0, 0]);
  assert.equal(x.meta.skipped.invalidLines, 1);
  const n = parseXyzRgb({ name: 'cloud.xyzn', text: '0 0 0 0 0 1\n1 0 0 0 1 0\n' });
  assert.deepEqual([...n.models[0].attrs.normals], [0, 0, 1, 0, 1, 0]);
  assert.throws(() => parsePts({ name: 'e.pts', text: '0\n' }), /no points were found/);
  assert.throws(() => parsePtx({ name: 'e.ptx', text: 'abc\ndef\n' }), /not a PTX file/);
});

test('truncated binary geospatial files fail with a worded error, never an internal one', async () => {
  for (const o of [{ rowsPerStrip: 2 }, { tile: [4, 4], compression: 8 }, { compression: 5, predictor: 3 }, { big: true, compression: 32773 }]) await assertTruncationSafe(makeTiff({ width: W, height: H, values: ELEV, geo: GEO, ...o }), (b) => parseGeoTiff({ name: 't.tif', bytes: b }), `GeoTIFF ${JSON.stringify(o)}`);
  const keys = new Uint8Array(new Uint16Array([1, 1, 0, 2, 1024, 0, 1, 1, 3072, 0, 1, 32631]).buffer);
  await assertTruncationSafe(makeLas([[1, 2, 3, 4, 2, 9, 9, 9], [4, 5, 6, 7, 2, 1, 1, 1]], { offset: [0, 0, 0], vlrs: [{ user: 'LASF_Projection', id: 34735, data: keys }] }), (b) => parseLas({ name: 'p.las', bytes: b }), 'LAS');
  const shp = makeShp(13, [[[[0, 0, 1], [1, 1, 2], [2, 0, 3]]], [[[5, 5, 0], [6, 6, 0]]]]);
  const dbf = makeDbf([{ name: 'NAME', type: 'C', len: 8 }], [['a'], ['b']]);
  await assertTruncationSafe(shp, (b) => parseShapefile({ name: 'r.shp', files: new Map([['r.shp', b], ['r.dbf', dbf]]) }), 'SHP');
  await assertTruncationSafe(dbf, (b) => parseShapefile({ name: 'r.shp', files: new Map([['r.shp', shp], ['r.dbf', b]]) }), 'DBF');
  await assertTruncationSafe(makeShp(5, [[[[0, 0], [0, 10], [10, 10], [10, 0], [0, 0]]]]), (b) => parseShapefile({ name: 'a.shp', bytes: b }), 'SHP polygon');
  const nc = makeNetcdf({ dims: { y: 2, x: 3 }, gatts: { Conventions: 'CF-1.8' }, vars: [
    { name: 'x', dims: ['x'], type: 'double', data: [0, 1, 2], atts: { units: 'm', axis: 'X' } }, { name: 'y', dims: ['y'], type: 'double', data: [0, 1], atts: { units: 'm', axis: 'Y' } },
    { name: 'z', dims: ['y', 'x'], type: 'float', data: [1, 2, 3, 4, 5, 6], atts: { units: 'm' } }] });
  await assertTruncationSafe(nc, (b) => parseNetcdf({ name: 'g.nc', bytes: b }), 'NetCDF');
  const dsbb = new Uint8Array(56 + 24);
  dsbb.set(enc('DSBB'));
  const dv = new DataView(dsbb.buffer);
  dv.setInt16(4, 3, true); dv.setInt16(6, 2, true);
  [0, 20, 100, 110, 1, 6].forEach((v, i) => dv.setFloat64(8 + 8 * i, v, true));
  await assertTruncationSafe(dsbb, (b) => parseSurferGrd({ name: 's.grd', bytes: b }), 'Surfer DSBB');
  await assertTruncationSafe(makeZip([{ name: 'doc.kml', data: enc(KML) }]), (b) => parseKml({ name: 's.kmz', bytes: b }), 'KMZ');
});
