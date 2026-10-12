import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import * as laz from '../../src/geometry/adapters/laz.js';
import { importFile } from '../../src/geometry/registry.js';
import { close } from './helpers.js';

// ---------------------------------------------------------------------------------------------------------
// Fixtures. Three small LAZ files in tests/mesh/fixtures, compressed with the LASzip implementation of the
// laz-rs project (Python package lazrs 0.8.2) and checked by reading them back with laspy 2.7.0 — that is,
// written by a different LASzip implementation than the laz-perf decoder under test:
//   cloud-format1-fixed-chunks.laz   LAS 1.2, point format 1, 60 points in fixed chunks of 25, GeoTIFF keys
//                                    for EPSG:32631, pointwise compressor (2)
//   cloud-format7-layered.laz        LAS 1.4, point format 7 (with colour), 40 points in variable-size chunks
//                                    of 10, 17 and 13, layered compressor (3)
//   cloud-format3-one-chunk.laz      LAS 1.2, point format 3 (with colour), 30 points in one chunk
// Point i of every file is, in integer record units (scale 0.01, offset 500 000 / 5 900 000 / 0):
//   X = 1000 + 37 i + i² mod 11,  Y = −500 + 53 i − 7 i mod 13,  Z = −4000 + 29 i mod 300,
//   intensity = 977 i mod 65536,  classification = (2, 5, 9, 40)[i mod 4]  (40 → 7 in formats 1 and 3, whose
//   classification field has five bits),  red = 4 i mod 256,  green = 255 − i,  blue = 128  (× 256 in the file)
// ---------------------------------------------------------------------------------------------------------

const fixture = async (name) => new Uint8Array(await readFile(new URL(`./fixtures/${name}`, import.meta.url)));
const vendor = (name) => readFile(new URL(`../../vendor/laz-perf/${name}`, import.meta.url));

const expected = (i, wideClass) => ({
  x: 500000 + 0.01 * (1000 + 37 * i + ((i * i) % 11)),
  y: 5900000 + 0.01 * (-500 + 53 * i - ((i * 7) % 13)),
  z: 0.01 * (-4000 + ((i * 29) % 300)),
  intensity: (i * 977) % 65536,
  classification: [2, 5, 9, wideClass ? 40 : 7][i % 4],
  rgb: [(i * 4) % 256, 255 - i, 128],
});

function checkPoints(model, indices, wideClass, hasColour) {
  const p = model.positions, o = model.origin || [0, 0, 0];
  assert.equal(p.length, 3 * indices.length);
  indices.forEach((i, k) => {
    const e = expected(i, wideClass);
    close(o[0] + p[3 * k], e.x, 1e-5, `x of point ${i}`);
    close(o[1] + p[3 * k + 1], e.y, 1e-5, `y of point ${i}`);
    close(o[2] + p[3 * k + 2], e.z, 1e-5, `z of point ${i}`);
    assert.equal(model.attrs.intensity[k], e.intensity, `intensity of point ${i}`);
    assert.equal(model.attrs.classification[k], e.classification, `classification of point ${i}`);
    if (hasColour) assert.deepEqual(Array.from(model.attrs.rgb.subarray(3 * k, 3 * k + 3)), e.rgb, `colour of point ${i}`);
  });
}

const range = (n, step = 1) => Array.from({ length: Math.ceil(n / step) }, (_, k) => k * step);

/** LAS public header (+ optional VLRs) with the compression bit set, no point data. */
function lazHeader({ minor = 2, fmt = 1, recLen = 28, count = 10, vlrs = [], compressed = true } = {}) {
  const hs = minor >= 4 ? 375 : 227;
  const body = vlrs.reduce((s, v) => s + 54 + v.data.length, 0);
  const b = new Uint8Array(hs + body);
  const dv = new DataView(b.buffer);
  b.set(new TextEncoder().encode('LASF'));
  b[24] = 1; b[25] = minor;
  dv.setUint16(94, hs, true); dv.setUint32(96, hs + body, true); dv.setUint32(100, vlrs.length, true);
  b[104] = fmt | (compressed ? 0x80 : 0);
  dv.setUint16(105, recLen, true); dv.setUint32(107, count, true);
  for (let k = 0; k < 3; k++) dv.setFloat64(131 + 8 * k, 0.001, true);
  let p = hs;
  for (const v of vlrs) {
    b.set(new TextEncoder().encode(v.user), p + 2);
    dv.setUint16(p + 18, v.id, true); dv.setUint16(p + 20, v.data.length, true);
    b.set(v.data, p + 54);
    p += 54 + v.data.length;
  }
  return b;
}

/** Payload of the "laszip encoded" record. */
function laszipRecord({ compressor = 2, coder = 0, chunk = 50000, items = [[6, 20, 2], [7, 8, 2]] } = {}) {
  const d = new Uint8Array(34 + 6 * items.length);
  const dv = new DataView(d.buffer);
  dv.setUint16(0, compressor, true); dv.setUint16(2, coder, true);
  d[4] = 3; d[5] = 4; dv.setUint16(6, 3, true);
  dv.setUint32(12, chunk, true);
  dv.setBigInt64(16, -1n, true); dv.setBigInt64(24, -1n, true);
  dv.setUint16(32, items.length, true);
  items.forEach(([type, size, version], k) => { dv.setUint16(34 + 6 * k, type, true); dv.setUint16(36 + 6 * k, size, true); dv.setUint16(38 + 6 * k, version, true); });
  return d;
}

test('vendored LAZ decoder matches the record in docs/vendor.md (size and SHA-256) and reports itself available', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  let total = 0;
  for (const f of ['laz-perf.js', 'laz-perf.wasm']) {
    const buf = await vendor(f);
    total += buf.length;
    const sha = createHash('sha256').update(buf).digest('hex');
    assert.ok(doc.includes(sha), `${f}: SHA-256 ${sha} is recorded`);
    assert.ok(doc.includes(String(buf.length)), `${f}: size ${buf.length} is recorded`);
  }
  assert.equal(laz.sizeBytes, total);
  assert.ok(doc.includes(createHash('sha256').update(await vendor('LICENSE.txt')).digest('hex')));
  assert.equal(laz.id, 'laz');
  assert.deepEqual(laz.formats, ['laz']);
  assert.deepEqual(laz.library, { name: 'laz-perf', version: '0.0.7', licence: 'Apache-2.0' });
  assert.ok(doc.includes('laz-perf') && doc.includes('0.0.7') && doc.includes('Apache-2.0'));
  assert.equal(laz.available(), true);
  const mod = await laz.load();
  assert.equal(typeof mod.LASZip, 'function');
  assert.equal(await laz.load(), mod);
});

test('header and laszip record: version, point format, counts, chunking and items of the three fixtures', async () => {
  const a = laz.readLazHeader(await fixture('cloud-format1-fixed-chunks.laz'), 'a.laz');
  assert.equal(a.version, '1.2'); assert.equal(a.pointFormat, 1); assert.equal(a.compressed, true);
  assert.equal(a.recordLength, 28); assert.equal(a.pointCount, 60); assert.equal(a.headerSize, 227);
  assert.deepEqual(a.scale, [0.01, 0.01, 0.01]); assert.deepEqual(a.offset, [500000, 5900000, 0]);
  assert.equal(a.systemId, 'OFFSHORE SURVEY');
  assert.deepEqual(a.vlrs.map((v) => [v.userId, v.recordId, v.length]), [['LASF_Projection', 34735, 32], ['laszip encoded', 22204, 46]]);
  assert.equal(a.dataOffset, 227 + 54 + 32 + 54 + 46);
  assert.deepEqual({ ...a.laszip, items: a.laszip.items.map((i) => [i.name, i.size, i.version]) }, {
    compressor: 2, compressorName: 'pointwise, chunked', coder: 0, version: '2.2.0', options: 0, chunkSize: 25, variableChunks: false, recordLength: 28,
    items: [['POINT10', 20, 2], ['GPSTIME11', 8, 2]],
  });
  close(a.min[0], 500010, 1e-9); close(a.max[2], -37.01, 1e-9);

  const b = laz.readLazHeader(await fixture('cloud-format7-layered.laz'));
  assert.equal(b.version, '1.4'); assert.equal(b.pointFormat, 7); assert.equal(b.headerSize, 375);
  assert.equal(b.pointCount, 40, 'taken from the 64-bit count of LAS 1.4');
  assert.equal(b.laszip.compressor, 3); assert.equal(b.laszip.variableChunks, true); assert.equal(b.laszip.recordLength, 36);
  assert.deepEqual(b.laszip.items.map((i) => i.name), ['POINT14', 'RGB14']);

  const c = laz.readLazHeader(await fixture('cloud-format3-one-chunk.laz'));
  assert.equal(c.pointFormat, 3); assert.equal(c.pointCount, 30); assert.equal(c.laszip.chunkSize, 50000);
  assert.deepEqual(c.laszip.items.map((i) => i.name), ['POINT10', 'GPSTIME11', 'RGB12']);

  // a hand-built header: the record is found among other VLRs and described without any decompression
  const h = laz.readLazHeader(lazHeader({ fmt: 0, recLen: 22, count: 7, vlrs: [{ user: 'other', id: 1, data: new Uint8Array(5) }, { user: 'laszip encoded', id: 22204, data: laszipRecord({ chunk: 0xffffffff, items: [[6, 20, 2], [0, 2, 2]] }) }] }));
  assert.equal(h.pointCount, 7); assert.equal(h.laszip.variableChunks, true); assert.equal(h.laszip.recordLength, 22);
  assert.deepEqual(h.laszip.items.map((i) => i.name), ['POINT10', 'BYTE']);
  assert.equal(h.laszip.version, '3.4.3');
  assert.equal(laz.readLazHeader(lazHeader({ vlrs: [] })).laszip, null);
  assert.equal(laz.sniff({ bytes: lazHeader() }), true);
  assert.equal(laz.sniff({ bytes: lazHeader({ compressed: false }) }), false);
  assert.equal(laz.sniff({ bytes: new Uint8Array(300) }), false);
});

test('format 1, fixed chunks of 25: all 60 points decoded exactly; CRS from the GeoTIFF keys; five-bit classes', async () => {
  const bytes = await fixture('cloud-format1-fixed-chunks.laz');
  const r = await laz.parse({ name: 'survey/cloud.laz', bytes });
  assert.equal(r.models.length, 1);
  const m = r.models[0];
  assert.equal(m.kind, 'points'); assert.equal(m.name, 'cloud');
  assert.ok(m.positions instanceof Float32Array);
  checkPoints(m, range(60), false, false);
  assert.equal(m.attrs.rgb, undefined);
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.format, 'laz');
  assert.equal(r.meta.pointFormat, 1); assert.equal(r.meta.version, '1.2');
  assert.equal(r.meta.pointCount, 60); assert.equal(r.meta.thinning, 1);
  assert.deepEqual(r.meta.skipped, { thinnedPoints: 0 });
  assert.deepEqual(r.meta.classes, { 2: 15, 5: 15, 7: 15, 9: 15 });
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.equal(r.meta.crs.epsg, 32631);
  assert.deepEqual(r.meta.compression, { compressor: 'pointwise, chunked', version: '2.2.0', chunkSize: 25, items: ['POINT10 v2', 'GPSTIME11 v2'] });
  assert.equal(r.meta.decoder, 'laz-perf 0.0.7');
  assert.equal(r.meta.counts.points, 60);
  close(r.meta.bbox.min[0], 500010, 1e-4); close(r.meta.bbox.min[2], -40, 1e-4);
});

test('format 7, layered compressor with chunks of 10 + 17 + 13 points: 8-bit classes, colour, LAS 1.4 count', async () => {
  const r = await laz.parse({ name: 'b.laz', bytes: await fixture('cloud-format7-layered.laz') });
  const m = r.models[0];
  checkPoints(m, range(40), true, true);
  assert.deepEqual(r.meta.classes, { 2: 10, 5: 10, 9: 10, 40: 10 });
  assert.equal(r.meta.pointFormat, 7); assert.equal(r.meta.version, '1.4');
  assert.equal(r.meta.compression.compressor, 'layered, chunked');
  assert.equal(r.meta.compression.chunkSize, 'variable');
  assert.equal(r.meta.crs, null);
  assert.deepEqual(r.warnings, ['The file declares no length unit; m was assumed.']);
});

test('format 3 in a single chunk; options of the LAS importer apply (units, depth-positive vertical)', async () => {
  const bytes = await fixture('cloud-format3-one-chunk.laz');
  const r = await laz.parse({ name: 'c.laz', bytes });
  checkPoints(r.models[0], range(30), false, true);
  const ft = await laz.parse({ name: 'c.laz', bytes }, { units: 'ft', vertical: 'depth' });
  const p = ft.models[0].positions, o = ft.models[0].origin || [0, 0, 0];
  const e = expected(5, false);
  close(o[0] + p[15], e.x * 0.3048, 1e-4); close(o[1] + p[16], e.y * 0.3048, 1e-4); close(o[2] + p[17], -e.z * 0.3048, 1e-5);
  assert.deepEqual(ft.meta.units, { unit: 'ft', scale: 0.3048, basis: 'option' });
});

test('thinning: every n-th record of the decoded sequence is kept so that at most maxPoints remain', async () => {
  const bytes = await fixture('cloud-format1-fixed-chunks.laz');
  const r = await laz.parse({ name: 'a.laz', bytes }, { maxPoints: 20 });
  checkPoints(r.models[0], range(60, 3), false, false);
  assert.equal(r.meta.thinning, 3); assert.equal(r.meta.pointCount, 60);
  assert.deepEqual(r.meta.skipped, { thinnedPoints: 40 });
  assert.equal(r.warnings[0], 'The cloud has 60 points; every 3rd point was kept (20) to stay under maxPoints.');
  // 60 points into at most 25: every third (20 kept), not every second (30)
  assert.equal((await laz.parse({ name: 'a.laz', bytes }, { maxPoints: 25 })).models[0].positions.length / 3, 20);
  const r7 = await laz.parse({ name: 'b.laz', bytes: await fixture('cloud-format7-layered.laz') }, { maxPoints: 7 });
  checkPoints(r7.models[0], range(40, 6), true, true);
  assert.equal(r7.models[0].positions.length / 3, 7);
  // 0 switches thinning off; the default (2 000 000) leaves small clouds alone
  assert.equal((await laz.parse({ name: 'a.laz', bytes }, { maxPoints: 0 })).meta.thinning, 1);
});

test('through the registry: a .laz file is routed to the adapter and reported as format "laz"', async () => {
  const bytes = await fixture('cloud-format7-layered.laz');
  const r = await importFile({ name: 'scan.laz', bytes });
  assert.equal(r.format, 'laz');
  assert.equal(r.meta.via, 'adapter:laz');
  assert.equal(r.models[0].kind, 'points');
  checkPoints(r.models[0], range(40), true, true);
  // the content decides, not the name
  assert.equal((await importFile({ name: 'scan.las', bytes })).format, 'laz');
});

test('files that cannot be decoded fail with the reason and what to do; the decoder recovers afterwards', async () => {
  const good = await fixture('cloud-format7-layered.laz');
  await assert.rejects(() => laz.parse({ name: 'notes.laz', bytes: new TextEncoder().encode('x y z\n1 2 3\n') }), /notes\.laz: not a LAS \/ LAZ file.*"LASF"/);
  await assert.rejects(() => laz.parse({ name: 'short.laz', bytes: good.slice(0, 100) }), /short\.laz: the file is truncated \(100 bytes/);
  await assert.rejects(() => laz.parse({ name: 't.laz', bytes: good.slice(0, 300) }), /header size in the file \(375\) is not possible for a 300-byte/);
  await assert.rejects(() => laz.parse({ name: 'novlr.laz', bytes: lazHeader({ vlrs: [] }) }), /novlr\.laz: .*no "laszip encoded" record \(VLR 22204\).*uncompressed LAS/);
  const withRec = (o, h = {}) => lazHeader({ ...h, vlrs: [{ user: 'laszip encoded', id: 22204, data: laszipRecord(o) }] });
  await assert.rejects(() => laz.parse({ name: 'old.laz', bytes: withRec({ compressor: 1 }) }), /original pointwise LASzip stream without chunks.*laszip -i in\.laz -o out\.laz/);
  await assert.rejects(() => laz.parse({ name: 'odd.laz', bytes: withRec({ compressor: 9 }) }), /unknown LASzip compressor 9/);
  await assert.rejects(() => laz.parse({ name: 'wave.laz', bytes: withRec({ items: [[6, 20, 2], [7, 8, 2], [9, 29, 1]] }, { fmt: 4, recLen: 57 }) }), /format 4 carries full-waveform packets/);
  await assert.rejects(() => laz.parse({ name: 'mix.laz', bytes: withRec({ compressor: 2, items: [[10, 30, 3]] }, { fmt: 6, recLen: 30, minor: 4 }) }), /format 6 must use the layered LASzip compressor \(3\)/);
  await assert.rejects(() => laz.parse({ name: 'len.laz', bytes: withRec({}, { recLen: 30 }) }), /describes points of 28 bytes but the header says 30/);
  await assert.rejects(() => laz.parse({ name: 'end.laz', bytes: withRec({}) }), /ends before the compressed point data/);
  // compressed data cut off part-way
  await assert.rejects(() => laz.parse({ name: 'cut.laz', bytes: good.slice(0, 700) }), /cut\.laz: the point records could not be decompressed after \d+ of 40 points.*truncated or corrupt.*uncompressed LAS/);
  // … and the next file is read normally
  const again = await laz.parse({ name: 'b.laz', bytes: good });
  checkPoints(again.models[0], range(40), true, true);
});

test('damaged compressed bytes cannot be detected by the format; points outside the header box are flagged', async () => {
  const bytes = (await fixture('cloud-format7-layered.laz')).slice();
  for (let i = 900; i < 1000; i++) bytes[i] ^= 0x5a;
  let r = null, failure = null;
  try { r = await laz.parse({ name: 'bad.laz', bytes }); } catch (e) { failure = e; }
  // either the decoder runs off the stream, or it returns coordinates that the header contradicts
  if (failure) assert.match(failure.message, /could not be decompressed/);
  else assert.ok(r.warnings.some((w) => /lie outside the bounding box stated in the file header/.test(w)), r.warnings.join(' | '));
});

test('an uncompressed LAS file with a .laz name is read without the decoder', async () => {
  // format 0, three points, written by hand
  const h = lazHeader({ fmt: 0, recLen: 20, count: 3, compressed: false });
  const b = new Uint8Array(h.length + 60);
  b.set(h);
  const dv = new DataView(b.buffer);
  [[1000, 2000, -500], [1500, 2500, -600], [2000, 3000, -700]].forEach((p, i) => { p.forEach((v, k) => dv.setInt32(227 + 20 * i + 4 * k, v, true)); b[227 + 20 * i + 15] = 2; });
  const r = await laz.parse({ name: 'plain.laz', bytes: b });
  assert.equal(r.meta.format, 'las');
  const m = r.models[0], o = m.origin || [0, 0, 0];
  close(o[0] + m.positions[3], 1.5, 1e-6); close(o[2] + m.positions[8], -0.7, 1e-6);
  assert.deepEqual(Array.from(m.attrs.classification), [2, 2, 2]);
});
