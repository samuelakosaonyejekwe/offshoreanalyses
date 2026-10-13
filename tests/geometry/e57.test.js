import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseE57 } from '../../src/geometry/importers/e57.js';
import { enc, close, concat, assertContract, assertTruncationSafe } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const world = (m, i) => [0, 1, 2].map((k) => m.positions[3 * i + k] + (m.origin ? m.origin[k] : 0));
const INTERNAL = /Cannot read|is not a function|is not iterable|is not defined|Invalid (typed )?array length|Offset is outside|Invalid array buffer length|Invalid DataView length|Start offset|Maximum call stack|out of memory|Invalid string length|source is too large|BigInt/i;

// ---------------------------------------------------------------------------------------------------------
// An E57 writer for the tests: 48-byte header, binary sections of data packets (optionally an index and an
// empty packet in between), the XML section, all cut into 1024-byte pages with their CRC-32C. Each field's
// bytestream is split over the packets at byte positions of its own, so that values straddle packets and the
// packets hold different numbers of values per field.
// ---------------------------------------------------------------------------------------------------------

const CRC_TABLE = Int32Array.from({ length: 256 }, (_, i) => { let c = i; for (let k = 0; k < 8; k++) c = c & 1 ? (c >>> 1) ^ 0x82f63b78 : c >>> 1; return c; });
function crc32c(b) {
  let c = ~0;
  for (let i = 0; i < b.length; i++) c = (c >>> 8) ^ CRC_TABLE[(c ^ b[i]) & 255];
  return ~c >>> 0;
}
const phys = (L) => Math.floor(L / 1020) * 1024 + (L % 1020);
const F = (name, precision = 'double') => ({ name, type: 'Float', precision });
const I = (name, min, max) => ({ name, type: 'Integer', min: BigInt(min), max: BigInt(max) });
const S = (name, min, max, scale, offset = 0) => ({ name, type: 'ScaledInteger', min: BigInt(min), max: BigInt(max), scale, offset });

function bytestream(f, values) {
  if (f.xml || f.type === 'Float') return new Uint8Array((f.precision === 'single' ? Float32Array : Float64Array).from(values).buffer); // fields declared by hand: doubles
  const bits = f.max === f.min ? 0 : (f.max - f.min).toString(2).length;
  let acc = 0n, pos = 0n;
  for (const v of values) { acc |= (BigInt(v) - f.min) << pos; pos += BigInt(bits); }
  const out = new Uint8Array(Math.ceil(Number(pos) / 8));
  for (let i = 0; i < out.length; i++) out[i] = Number((acc >> BigInt(8 * i)) & 255n);
  return out;
}
const fieldXml = (f) => (f.xml ? f.xml : f.type === 'Float' ? `<${f.name} type="Float"${f.precision === 'single' ? ' precision="single"' : ''}/>` : `<${f.name} type="${f.type}" minimum="${f.min}" maximum="${f.max}"${f.type === 'ScaledInteger' ? ` scale="${f.scale}" offset="${f.offset}"` : ''}/>`);

/**
 * @param {object[]} scans { name, guid, fields, values (one array per field; raw integers for integer fields), packets,
 *   filler (index + empty packet after the first data packet), pose: { q: [w,x,y,z], t }, xml (extra scan
 *   children), codecs, recordCount (override) }
 * @param {object} [o] { images, version, pageSize, mutate(logical bytes, info) }
 */
function makeE57(scans, o = {}) {
  const chunks = [new Uint8Array(48)];
  let L = 48;
  const put = (u) => { chunks.push(u); L += u.length; };
  let xml = '<?xml version="1.0" encoding="UTF-8"?>\n<e57Root type="Structure" xmlns="http://www.astm.org/COMMIT/E57/2010-e57-v1.0" xmlns:nor="http://www.libe57.org/E57_EXT_surface_normals.txt"><formatName type="String"><![CDATA[ASTM E57 3D Imaging Data File]]></formatName><guid type="String"><![CDATA[{FILE-GUID}]]></guid><versionMajor type="Integer">1</versionMajor><versionMinor type="Integer"/><e57LibraryVersion type="String"><![CDATA[unit test]]></e57LibraryVersion>' + (o.rootXml || '') + '<data3D type="Vector" allowHeterogeneousChildren="1">';
  const sections = [];
  for (const s of scans) {
    const n = s.values.length ? s.values[0].length : 0, P = s.packets || 1;
    const streams = s.fields.map((f, k) => bytestream(f, s.values[k] || []));
    const pk = [];
    for (let p = 0; p < P; p++) {
      const parts = streams.map((b) => b.subarray(Math.floor((b.length * p) / P), Math.floor((b.length * (p + 1)) / P)));
      let len = 6 + 2 * parts.length + parts.reduce((a, b) => a + b.length, 0);
      len += (4 - (len % 4)) % 4;
      const u = new Uint8Array(len), d = new DataView(u.buffer);
      u[0] = 1; d.setUint16(2, len - 1, true); d.setUint16(4, parts.length, true);
      let q = 6 + 2 * parts.length;
      parts.forEach((b, k) => { d.setUint16(6 + 2 * k, b.length, true); u.set(b, q); q += b.length; });
      pk.push(u);
      if (s.filler && p === 0) {
        const ix = new Uint8Array(32); ix[0] = 0; ix[2] = 31; ix[4] = 1; // index packet with one entry
        const em = new Uint8Array(8); em[0] = 2; em[2] = 7;               // empty packet
        pk.push(ix, em);
      }
    }
    const body = concat(...pk), h = new Uint8Array(32), hd = new DataView(h.buffer);
    h[0] = 1; hd.setBigUint64(8, BigInt(32 + body.length), true); hd.setBigUint64(16, BigInt(phys(L + 32)), true);
    sections.push({ start: L, data: L + 32, end: L + 32 + body.length });
    const pose = s.pose ? `<pose type="Structure"><rotation type="Structure">${['w', 'x', 'y', 'z'].map((c, k) => `<${c} type="Float">${s.pose.q[k]}</${c}>`).join('')}</rotation><translation type="Structure">${['x', 'y', 'z'].map((c, k) => `<${c} type="Float">${s.pose.t[k]}</${c}>`).join('')}</translation></pose>` : '';
    xml += `<vectorChild type="Structure"><guid type="String"><![CDATA[${s.guid || ''}]]></guid>${s.name ? `<name type="String"><![CDATA[${s.name}]]></name>` : ''}${pose}${s.xml || ''}<points type="CompressedVector" fileOffset="${phys(L)}" recordCount="${s.recordCount ?? n}"><prototype type="Structure">${s.fields.map(fieldXml).join('')}</prototype><codecs type="Vector" allowHeterogeneousChildren="1">${s.codecs || ''}</codecs></points></vectorChild>`;
    put(h); put(body);
  }
  xml += `</data3D><images2D type="Vector" allowHeterogeneousChildren="1">${'<vectorChild type="Structure"><guid type="String">image</guid></vectorChild>'.repeat(o.images || 0)}</images2D></e57Root>\n`;
  const xb = enc(xml), xmlL = L;
  put(xb);
  const logical = concat(...chunks), nPages = Math.ceil(L / 1020), hv = new DataView(logical.buffer);
  logical.set(enc('ASTM-E57'));
  hv.setUint32(8, o.version ?? 1, true); hv.setBigUint64(16, BigInt(nPages * 1024), true); hv.setBigUint64(24, BigInt(phys(xmlL)), true); hv.setBigUint64(32, BigInt(xb.length), true); hv.setBigUint64(40, BigInt(o.pageSize ?? 1024), true);
  if (o.mutate) o.mutate(logical, { xmlL, sections });
  const out = new Uint8Array(nPages * 1024), ov = new DataView(out.buffer);
  for (let p = 0; p < nPages; p++) {
    out.set(logical.subarray(p * 1020, Math.min(L, (p + 1) * 1020)), p * 1024);
    ov.setUint32(p * 1024 + 1020, crc32c(out.subarray(p * 1024, p * 1024 + 1020)), false);
  }
  return out;
}

const parse = (bytes, opts) => parseE57({ name: 'scan.e57', bytes }, opts);
const XYZ = [[1, 2, 3], [-4.5, 5.25, 6.125], [7, -8, 9], [10.5, 11, -12], [0.001, 0.002, 0.003]];
const col = (k) => XYZ.map((p) => p[k]);
// 90° about z, then a shift: (x, y, z) → (−y + 100, x + 200, z + 300)
const POSE = { q: [Math.SQRT1_2, 0, 0, Math.SQRT1_2], t: [100, 200, 300] };

test('E57: float coordinates (single and double), pose, invalid state, intensity, colour, unused fields', () => {
  const bytes = makeE57([{
    name: 'station 1', guid: '{SCAN-1}', pose: POSE, packets: 3, filler: true,
    xml: '<sensorVendor type="String"><![CDATA[Test]]></sensorVendor><sensorModel type="String"><![CDATA[TLS-1]]></sensorModel><intensityLimits type="Structure"><intensityMinimum type="Float">0</intensityMinimum><intensityMaximum type="Float">1</intensityMaximum></intensityLimits><colorLimits type="Structure"><colorRedMinimum type="Integer">0</colorRedMinimum><colorRedMaximum type="Integer">255</colorRedMaximum><colorGreenMinimum type="Integer">0</colorGreenMinimum><colorGreenMaximum type="Integer">255</colorGreenMaximum><colorBlueMinimum type="Integer">0</colorBlueMinimum><colorBlueMaximum type="Integer">255</colorBlueMaximum></colorLimits>',
    fields: [F('cartesianX', 'single'), F('cartesianY'), F('cartesianZ', 'single'), I('cartesianInvalidState', 0, 2), F('intensity', 'single'), I('colorRed', 0, 255), I('colorGreen', 0, 255), I('colorBlue', 0, 255), I('rowIndex', 0, 1000), I('columnIndex', 0, 70000), F('timeStamp'), I('returnIndex', 0, 0), { name: 'nor:normalX', type: 'Float', precision: 'single' }],
    values: [col(0), col(1), col(2), [0, 0, 2, 0, 1], [0.1, 0.2, 0.3, 0.4, 0.5], [255, 0, 9, 128, 9], [0, 255, 9, 64, 9], [10, 20, 9, 32, 9], [0, 1, 2, 3, 1000], [69999, 5, 4, 3, 2], [1e9, 2e9, 3e9, 4e9, 5e9], [0, 0, 0, 0, 0], [0, 0, 1, 1, 0]],
  }], { images: 2 });
  assert.equal(bytes.length % 1024, 0);
  const r = parse(bytes);
  assertContract(r, ['points']);
  const m = r.models[0];
  assert.equal(m.name, 'scan');
  assert.equal(m.positions.length, 9);
  [0, 1, 3].forEach((src, i) => { const [x, y, z] = XYZ[src], w = world(m, i); close(w[0], -y + 100, 1e-4); close(w[1], x + 200, 1e-4); close(w[2], z + 300, 1e-4); });
  assert.ok(m.attrs.intensity instanceof Float32Array && m.attrs.rgb instanceof Uint8Array);
  assert.deepEqual([...m.attrs.intensity], [0.1, 0.2, 0.4].map(Math.fround));
  assert.deepEqual([...m.attrs.rgb], [255, 0, 10, 0, 255, 20, 128, 64, 32]);
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.equal(r.meta.format, 'e57');
  assert.equal(r.meta.version, '1.0');
  assert.equal(r.meta.pointCount, 3);
  assert.equal(r.meta.thinning, 1);
  assert.deepEqual(r.meta.skipped, { thinnedPoints: 0, invalidPoints: 2, images2D: 2, scans: 0 });
  assert.deepEqual(r.meta.elementTypes, { point1: 3 });
  const s = r.meta.scans[0];
  assert.deepEqual([s.name, s.guid, s.records, s.points, s.invalid, s.coordinates, s.sensor], ['station 1', '{SCAN-1}', 5, 3, 2, 'cartesian', 'Test TLS-1']);
  assert.deepEqual(s.intensityLimits, [0, 1]);
  close(s.pose.rotation[0], Math.SQRT1_2); assert.deepEqual(s.pose.translation, [100, 200, 300]);
  assert.deepEqual(r.meta.ignoredFields, { rowIndex: 1, columnIndex: 1, timeStamp: 1, returnIndex: 1, 'nor:normalX': 1 });
  assert.deepEqual(r.meta.checksums.failed, 0);
  assert.equal(r.meta.library, 'unit test');
  assert.ok(r.warnings.some((w) => /^2 point\(s\) flagged invalid/.test(w)));
  assert.ok(r.warnings.some((w) => /^2 image\(s\) \(images2D\) were not read/.test(w)));
  // metres by definition; the user may still override, and depth flips z
  const mm = parse(bytes, { units: 'mm', vertical: 'depth' });
  close(world(mm.models[0], 0)[2], -0.303, 1e-6);
});

test('E57: scaled integers and bit-packed integers of every width, values straddling packets', () => {
  // an LCG so that the fixture is the same on every run
  let seed = 12345n;
  const rnd = (range) => { seed = (seed * 6364136223846793005n + 1442695040888963407n) & 0xffffffffffffffffn; return range === 0n ? 0n : (seed >> 3n) % (range + 1n); };
  for (const w of [1, 2, 7, 8, 13, 19, 24, 25, 31, 32, 33, 40, 48, 53, 54, 63, 64]) {
    const min = -(2n ** BigInt(w - 1)) + 3n, max = min + 2n ** BigInt(w) - 1n, scale = 200 / 2 ** w, N = 41;
    const raw = Array.from({ length: N }, (_, i) => (i === 0 ? min : i === 1 ? max : min + rnd(max - min)));
    const ys = Array.from({ length: N }, () => 100000n + rnd(999n)), zs = Array.from({ length: N }, () => -5n + rnd(1023n)), ks = Array.from({ length: N }, () => rnd(4095n));
    const r = parse(makeE57([{ packets: 4, filler: true, fields: [S('cartesianX', min, max, scale, 0.5), I('pad', 0, 6), S('cartesianY', 100000, 100999, 0.001, -100), I('cartesianZ', -5, 1018), S('const', 7, 7, 2, 1), I('intensity', 0, 4095)], values: [raw, ks.map((v) => v % 7n), ys, zs, raw.map(() => 7n), ks] }]));
    const m = r.models[0];
    assert.equal(m.positions.length, 3 * N, `width ${w}`);
    for (let i = 0; i < N; i++) {
      const p = world(m, i);
      close(p[0], Number(raw[i]) * scale + 0.5, 2e-5, `width ${w} x[${i}]`);
      close(p[1], Number(ys[i]) * 0.001 - 100, 1e-5, `width ${w} y[${i}]`);
      close(p[2], Number(zs[i]), 1e-5, `width ${w} z[${i}]`);
      assert.equal(m.attrs.intensity[i], Number(ks[i]), `width ${w} intensity[${i}]`);
    }
    assert.deepEqual(r.meta.scans[0].intensityLimits, [0, 4095]);
    assert.deepEqual(r.warnings, []);
  }
  // 64-bit raw values keep every bit until the scale is applied: 2^62 + 1 and 2^62 + 3 differ after × 1 − 2^62
  const big = 2n ** 62n;
  const r = parse(makeE57([{ fields: [S('cartesianX', -(2n ** 63n), 2n ** 63n - 1n, 1, 0), F('cartesianY'), F('cartesianZ'), I('intensity', big, big + 15n)], values: [[1n, -2n, 3n], [0, 0, 0], [0, 0, 0], [big + 1n, big + 3n, big + 15n]] }]));
  assert.deepEqual([0, 1, 2].map((i) => world(r.models[0], i)[0]), [1, -2, 3]);
  assert.deepEqual(r.meta.scans[0].intensityLimits, [2 ** 62, 2 ** 62 + 15]);
});

test('E57: spherical scans, several scans merged, colour ranges, scans without colour', () => {
  const sph = [[2, 0, 0], [2, Math.PI / 2, 0], [4, Math.PI, Math.PI / 6], [1, 0, Math.PI / 2], [9, 1, 1]];
  const bytes = makeE57([
    { name: 'A', guid: 'a', pose: POSE, fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ'), I('colorRed', 0, 65535), I('colorGreen', 0, 65535), I('colorBlue', 0, 65535)], values: [[1, 2], [0, 0], [0, 0], [65535, 0], [0, 32768], [257, 65535]] },
    { name: 'B', guid: 'b', packets: 2, fields: [S('sphericalRange', 0, 100000, 0.001), F('sphericalAzimuth'), F('sphericalElevation', 'single'), I('sphericalInvalidState', 0, 1), F('intensity')], values: [sph.map((p) => Math.round(p[0] * 1000)), sph.map((p) => p[1]), sph.map((p) => p[2]), [0, 0, 0, 0, 1], [5, 6, 7, 8, 9]] },
    { name: 'C', guid: 'c', pose: { q: [1, 0, 0, 0], t: [0, 0, -50] }, xml: '<colorLimits type="Structure"><colorRedMinimum type="Float">0</colorRedMinimum><colorRedMaximum type="Float">1</colorRedMaximum><colorGreenMinimum type="Float">0</colorGreenMinimum><colorGreenMaximum type="Float">1</colorGreenMaximum><colorBlueMinimum type="Float">0</colorBlueMinimum><colorBlueMaximum type="Float">1</colorBlueMaximum></colorLimits>', fields: [F('cartesianX', 'single'), F('cartesianY', 'single'), F('cartesianZ', 'single'), F('colorRed', 'single'), F('colorGreen', 'single'), F('colorBlue', 'single'), I('cartesianInvalidState', 0, 0)], values: [[7], [8], [9], [1], [0.5], [0], [0]] },
    { name: 'empty', guid: 'd', fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ')], values: [] },
  ]);
  const r = parse(bytes);
  assertContract(r, ['points']);
  const m = r.models[0];
  assert.equal(r.models.length, 1);
  assert.equal(m.positions.length / 3, 7);
  const want = [[100, 201, 300], [100, 202, 300], [2, 0, 0], [0, 2, 0], [-4 * Math.cos(Math.PI / 6), 0, 2], [0, 0, 1], [7, 8, -41]];
  want.forEach((p, i) => world(m, i).forEach((v, k) => close(v, p[k], 1e-5, `point ${i}`)));
  assert.deepEqual([...m.attrs.rgb], [255, 0, 1, 0, 128, 255, ...new Array(12).fill(160), 255, 128, 0]);
  assert.deepEqual([...m.attrs.intensity], [0, 0, 5, 6, 7, 8, 0]);
  assert.deepEqual(r.meta.scans.map((s) => [s.name, s.guid, s.records, s.points, s.coordinates]), [['A', 'a', 2, 2, 'cartesian'], ['B', 'b', 5, 4, 'spherical'], ['C', 'c', 1, 1, 'cartesian'], ['empty', 'd', 0, 0, 'cartesian']]);
  assert.equal(r.meta.skipped.invalidPoints, 1);
  assert.deepEqual(r.meta.attributes, ['intensity', 'rgb']);
  assert.ok(r.warnings.some((w) => /^1 of the 4 scans carry no colour; their points are grey/.test(w)));
  assert.ok(r.warnings.some((w) => /^2 of the 4 scans carry no intensity/.test(w)));
  close(r.meta.bbox.min[2], -41, 1e-5); close(r.meta.bbox.max[1], 202, 1e-5);
});

test('E57: a projected system in coordinateMetadata becomes the CRS; lengths stay metres', () => {
  const wkt = 'PROJCS["NAD83 / Texas South Central (ftUS)",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["standard_parallel_1",30.28],PARAMETER["standard_parallel_2",28.38],PARAMETER["latitude_of_origin",27.83],PARAMETER["central_meridian",-99],PARAMETER["false_easting",1968500],PARAMETER["false_northing",13123333.333],UNIT["US survey foot",0.3048006096012192],AUTHORITY["EPSG","2278"]]';
  const scan = { fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ')], values: [[950000.25, 950010], [4200000.5, 4200020], [12.5, 13]] };
  const r = parse(makeE57([scan], { rootXml: `<coordinateMetadata type="String"><![CDATA[${wkt}]]></coordinateMetadata>` }));
  assertContract(r, ['points']);
  assert.deepEqual([r.meta.crs.kind, r.meta.crs.epsg, r.meta.crs.units], ['projected', 2278, 'm']);
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.equal(r.meta.coordinateMetadata, wkt);
  world(r.models[0], 0).forEach((v, k) => close(v, [950000.25, 4200000.5, 12.5][k], 1e-9));
  assert.deepEqual(r.models[0].origin, [950005, 4200010, 13]);
  const g = parse(makeE57([scan], { rootXml: '<coordinateMetadata type="String"><![CDATA[GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]]]></coordinateMetadata>' }));
  assert.equal(g.meta.crs, null);
});

test('E57: thinning to maxPoints keeps every n-th valid point, as the LAS reader does', () => {
  const n = 11, xs = Array.from({ length: n }, (_, i) => i);
  const bytes = makeE57([{ fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ'), I('cartesianInvalidState', 0, 2)], values: [xs.slice(0, 6), xs.slice(0, 6), xs.slice(0, 6), [0, 2, 0, 0, 0, 0]] }, { packets: 2, fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ')], values: [xs.slice(6), xs.slice(6), xs.slice(6)] }]);
  const all = parse(bytes);
  assert.deepEqual([...all.models[0].positions].filter((_, i) => i % 3 === 0), [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
  const t = parse(bytes, { maxPoints: 4 });
  // 11 records → stride 3 over the 10 valid points
  assert.deepEqual([...t.models[0].positions].filter((_, i) => i % 3 === 0), [0, 4, 7, 10]);
  assert.equal(t.meta.thinning, 3);
  assert.equal(t.meta.pointCount, 10);
  assert.deepEqual(t.meta.skipped, { thinnedPoints: 6, invalidPoints: 1, images2D: 0, scans: 0 });
  assert.ok(t.warnings.includes('The cloud has 10 points; every 3rd point was kept (4) to stay under maxPoints.'));
  assert.equal(parse(bytes, { maxPoints: 0 }).models[0].positions.length, 30);
  assert.equal(parse(bytes, { maxPoints: 5 }).meta.thinning, 3);
  assert.equal(parse(bytes, { maxPoints: 6 }).meta.thinning, 2);
});

test('E57: unsupported variants and damaged files are refused in words', () => {
  const xyz = [F('cartesianX'), F('cartesianY'), F('cartesianZ')], v = [[1, 2], [3, 4], [5, 6]];
  const one = (scan, o) => makeE57([{ fields: xyz, values: v, ...scan }], o);
  assert.throws(() => parse(enc('not an E57 file, but long enough to hold a header of 48 bytes')), /scan\.e57: not an E57 file \(expected the signature "ASTM-E57"/);
  assert.throws(() => parse(new Uint8Array(10)), /not an E57 file/);
  assert.throws(() => parse(one({}, { version: 2 })), /E57 format version 2\.0 is not known/);
  assert.throws(() => parse(one({}, { pageSize: 4096 })), /page size of 4096 bytes; E57 files use 1024/);
  // a writer that crashed leaves the XML offset at zero
  assert.throws(() => parse(one({}, { mutate: (b) => b.fill(0, 24, 40) })), /gives no XML section .* was not closed by the program that wrote it/);
  assert.throws(() => parse(one({}, { mutate: (b) => new DataView(b.buffer).setBigUint64(24, 1n << 40n, true) })), /XML section \(offset 1099511627776, \d+ bytes\) lies outside the file/);
  assert.throws(() => parse(one({}, { mutate: (b, i) => b.set(enc('<html>'), i.xmlL + 39) })), /starts with <html> where <e57Root> was expected/);
  assert.throws(() => parse(makeE57([])), /the file holds no scans \(data3D is empty\)/);
  assert.throws(() => parse(makeE57([], { images: 3 })), /holds no scans \(data3D is empty; it holds 3 image\(s\), which are not read\)/);
  assert.throws(() => parse(one({ values: [] })), /no points were found \(the 1 scan\(s\) hold no point records\)/);
  // field types and codecs
  assert.throws(() => parse(one({ fields: [{ name: 'cartesianX', type: 'String', xml: '<cartesianX type="String"/>' }, xyz[1], xyz[2]] })), /Scan 1 was not read: field cartesianX: element type "String" is not one this reader decodes \(Float, Integer, ScaledInteger\)/);
  assert.throws(() => parse(one({ fields: [{ name: 'cartesianX', type: 'Float', xml: '<cartesianX type="Float" precision="half"/>' }, xyz[1], xyz[2]] })), /Float precision "half" is not known/);
  assert.throws(() => parse(one({ fields: [{ name: 'cartesianX', type: 'Integer', xml: '<cartesianX type="Integer" minimum="-9223372036854775808" maximum="18446744073709551616"/>' }, xyz[1], xyz[2]] })), /field cartesianX: the range .* needs 65 bits per value; E57 integers have at most 64/);
  assert.throws(() => parse(one({ fields: [{ name: 'cartesianX', type: 'Integer', xml: '<cartesianX type="Integer" minimum="5" maximum="1"/>' }, xyz[1], xyz[2]] })), /maximum \(1\) is below minimum \(5\)/);
  assert.throws(() => parse(one({ fields: [{ name: 'cartesianX', type: 'Integer', xml: '<cartesianX type="Integer" minimum="abc" maximum="1"/>' }, xyz[1], xyz[2]] })), /minimum="abc" is not an integer/);
  assert.throws(() => parse(one({ codecs: '<vectorChild type="Structure"><inputs type="Vector"/><codec type="Structure"><cdc:zip type="Integer" xmlns:cdc="x"/></codec></vectorChild>' })), /its points name 1 codec\(s\) in the codecs list; only the standard encoding .* not E57 extension codecs/);
  assert.throws(() => parse(one({ fields: [xyz[0], xyz[1], F('intensity')] })), /Scan 1 was not read: its records have no complete coordinate triple \(fields: cartesianX, cartesianY, intensity\)/);
  assert.throws(() => parse(one({ fields: [xyz[0], xyz[1], xyz[2], I('cartesianInvalidState', 0, 2)], values: [...v, [2, 1]] })), /all 2 points of the 1 scan\(s\) are flagged invalid/);
  assert.throws(() => parse(one({}, { mutate: (b, i) => { b[i.sections[0].start] = 9; } })), /binary section at offset \d+ does not start with a CompressedVector header \(id 9\)/);
  assert.throws(() => parse(one({}, { mutate: (b, i) => { b[i.sections[0].data] = 7; } })), /a packet of unknown type 7 follows after 0 of its 2 records/);
  assert.throws(() => parse(one({}, { mutate: (b, i) => { b[i.sections[0].data + 4] = 5; } })), /a data packet holds 5 bytestreams where the prototype has 3 fields/);
  assert.throws(() => parse(one({}, { mutate: (b, i) => { b[i.sections[0].data + 7] = 200; } })), /bytestreams of a data packet are longer than the packet/);
  // a field of an unknown type that is not used does not matter: its bytestream is skipped
  const odd = parse(one({ fields: [...xyz, { name: 'ext:label', type: 'String', xml: '<ext:label type="String"/>' }, { name: 'timeStamp', type: 'Integer', xml: '<timeStamp type="Integer" minimum="0" maximum="99999999999999999999999"/>' }] }));
  assert.equal(odd.models[0].positions.length, 6);
  assert.deepEqual(odd.meta.ignoredFields, { 'ext:label': 1, timeStamp: 1 });
  // a scan that cannot be read does not take the others with it
  const two = parse(makeE57([{ name: 'bad', fields: xyz, values: v, codecs: '<vectorChild type="Structure"/>' }, { name: 'good', fields: xyz, values: v }, { fields: [F('rowIndex')], values: [[1]] }]));
  assert.equal(two.models[0].positions.length, 6);
  assert.equal(two.meta.skipped.scans, 2);
  assert.ok(two.warnings.some((w) => /^Scan 1 \("bad"\) was not read: its points name 1 codec/.test(w)));
  assert.ok(two.warnings.some((w) => /^Scan 3 was not read: its records have no complete coordinate triple/.test(w)));
  // a record count larger than the data: what is there is read
  const over = parse(one({ recordCount: 1000 }));
  assert.equal(over.models[0].positions.length, 6);
  assert.ok(over.warnings.some((w) => /^Scan 1 declares 1000 records but its binary section has room for 2; read 2\./.test(w)));
  const lost = parse(makeE57([{ fields: xyz, values: [[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]], packets: 2 }], { mutate: (b, i) => new DataView(b.buffer).setBigUint64(i.sections[0].start + 8, 32n + 64n, true) }));
  assert.equal(lost.models[0].positions.length, 6);
  assert.ok(lost.warnings.some((w) => /declares 4 records but its binary section has room for 2; read 2/.test(w)));
});

test('E57: page checksums — an XML page that fails is an error, data pages are counted', () => {
  const n = 600, xs = Array.from({ length: n }, (_, i) => i * 0.5);
  const good = makeE57([{ fields: [F('cartesianX'), F('cartesianY'), F('cartesianZ')], values: [xs, xs, xs] }]);
  assert.ok(good.length > 14 * 1024);
  assert.equal(parse(good).meta.checksums.failed, 0);
  const data = good.slice(); data[3 * 1024 + 500] ^= 0x10; data[7 * 1024 + 1021] ^= 0x01;
  const r = parse(data);
  assert.equal(r.models[0].positions.length, 3 * n);
  assert.equal(r.meta.checksums.failed, 2);
  assert.ok(r.warnings.some((w) => /^2 of the \d+ data pages read fail their CRC-32C checksum; the points stored in them may be wrong\./.test(w)));
  const xml = good.slice(); xml[xml.length - 600] ^= 0x20;
  assert.throws(() => parse(xml), /scan\.e57: 1 of the \d+ pages of the XML section fail their CRC-32C checksum — the file is corrupt\./);
});

test('E57: truncated and garbled files fail with a worded error, never an internal one', async () => {
  const xs = Array.from({ length: 300 }, (_, i) => i);
  const scans = [
    { name: 'A', pose: POSE, packets: 3, filler: true, fields: [S('cartesianX', 0, 100000, 0.001), F('cartesianY'), F('cartesianZ', 'single'), I('cartesianInvalidState', 0, 2), I('colorRed', 0, 255), I('colorGreen', 0, 255), I('colorBlue', 0, 255), S('intensity', 0, 2n ** 60n, 1e-12)], values: [xs, xs, xs, xs.map((i) => i % 3), xs.map((i) => i % 256), xs.map((i) => i % 256), xs.map((i) => i % 256), xs] },
    { name: 'B', packets: 2, fields: [F('sphericalRange'), F('sphericalAzimuth'), F('sphericalElevation'), I('rowIndex', 0, 5000)], values: [xs, xs, xs, xs] },
  ];
  const bytes = makeE57(scans, { images: 1 });
  assert.equal(parse(bytes).models[0].positions.length, 3 * 400);
  await assertTruncationSafe(bytes, (b) => parse(b), 'E57');
  await assertTruncationSafe(real('e57-pdal-two-scans.e57'), (b) => parse(b), 'E57 (PDAL)');
  // a file cut on a page boundary keeps valid pages but loses the XML section at its end
  assert.throws(() => parse(bytes.slice(0, 4096)), /XML section \(offset \d+, \d+ bytes\) lies outside the file \(4096 bytes; the header declares \d+\)/);
  // garbled content with correct checksums: flip bytes of the logical stream before the pages are written
  let seed = 99;
  const rnd = (n) => { seed = (Math.imul(seed, 1103515245) + 12345) >>> 0; return (seed >>> 8) % n; };
  let read = 0, refused = 0, slowest = 0;
  for (let it = 0; it < 300; it++) {
    const g = makeE57(scans, { images: 1, mutate: (b, info) => {
      const region = it % 4; // header, packets, XML, anywhere
      const [lo, hi] = region === 0 ? [8, 48] : region === 1 ? [48, info.xmlL] : region === 2 ? [info.xmlL, b.length] : [0, b.length];
      for (let k = 0, m = 1 + rnd(6); k < m; k++) b[lo + rnd(hi - lo)] = rnd(4) ? rnd(256) : 0xff;
    } });
    const tp = Date.now();
    try {
      const r = parse(g);
      assert.ok(Array.isArray(r.models) && r.models[0].positions instanceof Float32Array && Array.isArray(r.warnings));
      read++;
    } catch (e) {
      assert.ok(e instanceof Error && !(e instanceof TypeError) && !(e instanceof RangeError) && !(e instanceof ReferenceError) && !(e instanceof SyntaxError) && !INTERNAL.test(e.message) && e.message.length > 15, `garbled file ${it}: ${e && e.constructor.name}: ${e && e.message}\n${e && e.stack}`);
      assert.match(e.message, /^scan\.e57: /);
      refused++;
    }
    slowest = Math.max(slowest, Date.now() - tp);
  }
  assert.ok(read > 40 && refused > 40, `garbled files: ${read} read, ${refused} refused`);
  assert.ok(slowest < 5000, `no garbled file takes long (slowest ${slowest} ms)`);
});

// ---------------------------------------------------------------------------------------------------------
// Real files (sources and licences in fixtures/real/README.md). The expected numbers are those of pye57
// (libE57Format) on the same file, or of the "e57" Python package where pye57 does not read the file.
// ---------------------------------------------------------------------------------------------------------

test('real E57, CloudCompare (via the PDAL test data): two scans, float coordinates, intensity and colour', () => {
  const r = parseE57({ name: 'A_moved_B.e57', bytes: real('e57-pdal-two-scans.e57') });
  assertContract(r, ['points']);
  const m = r.models[0];
  assert.equal(m.positions.length, 18);
  assert.deepEqual(r.meta.scans.map((s) => [s.name, s.guid, s.records, s.points]), [['B2 - Cloud', '{7c970a20-799d-45a7-9147-10319a9d4fcc}', 2, 2], ['A4 - Cloud', '{13288ae4-79e2-414f-ae0c-dfd05eabb5da}', 4, 4]]);
  assert.match(r.meta.scans[0].description, /CloudCompare/);
  assert.deepEqual([...m.positions.slice(0, 9)], [-42.63949966430664, 1.8725999593734741, 0.666700005531311, -42.64229965209961, 2.057800054550171, 0.5101000070571899, -43.46733474731445, 0.922502875328064, 0.3357999920845032]);
  assert.deepEqual([...m.attrs.intensity.slice(0, 2)], [0, 0.873046875]);
  assert.deepEqual([...m.attrs.rgb.slice(0, 6)], [255, 255, 0, 255, 255, 0]);
  assert.deepEqual(r.meta.bbox, { min: [-44.11943817138672, 0.7053028345108032, -0.607699990272522], max: [-42.63949966430664, 2.057800054550171, 0.666700005531311] });
  assert.deepEqual(r.warnings, []);
});

test('real E57, LAStools las2e57: scaled-integer coordinates in 10 bits, 16-bit colour, LAS extension fields', () => {
  const r = parseE57({ name: 'ColourRepresentation.e57', bytes: real('e57-las2e57-colour.e57') });
  assertContract(r, ['points']);
  const m = r.models[0];
  assert.equal(m.positions.length / 3, 153);
  assert.equal(r.meta.scans[0].sensor, 'LAStools (c) by rapidlasso GmbH');
  [-0.5, -0.015, -0.432, -0.5, 0.169, 0.485, -0.5, -0.16, 0.026, -0.5, 0.433, 0.35].forEach((v, i) => close(m.positions[i], v, 1e-7));
  assert.deepEqual(r.meta.bbox, { min: [-0.5, -0.5, -0.5], max: [0.5, 0.5, 0.5] });
  // 65280 of 65535
  assert.deepEqual([...m.attrs.rgb.slice(0, 6)], [0, 0, 254, 0, 0, 254]);
  assert.deepEqual(r.meta.ignoredFields, { returnIndex: 1, returnCount: 1, 'las:pointSourceId': 1 });
});

test('real E57, libE57Format 3: spherical coordinates with invalid points; ScaledInteger intensity and limits', () => {
  const s = parseE57({ name: 'tiny_spherical.e57', bytes: real('e57-spherical.e57') });
  assertContract(s, ['points']);
  assert.equal(s.models[0].positions.length / 3, 180);
  assert.deepEqual([s.meta.scans[0].records, s.meta.scans[0].coordinates, s.meta.skipped.invalidPoints], [360, 'spherical', 180]);
  [1, 0, 0, 0.9996957222246887, 0.017440905663028396, 0.017443559708706372, 0.9987832592386133, 0.034860583806158564, 0.0348818113260568].forEach((v, i) => close(s.models[0].positions[i], v, 1e-7));
  close(s.meta.bbox.min[1], -0.499998573193859, 1e-7); close(s.meta.bbox.max[2], 0.9999996829318346, 1e-7);
  const i = parseE57({ name: 'scaled_integer_intensity.e57', bytes: real('e57-scaled-intensity.e57') });
  assert.deepEqual([...i.models[0].positions], [1.1, 2.2, 3.3, 4.4, 5.5, 6.6].map(Math.fround));
  assert.deepEqual([...i.models[0].attrs.intensity], [-66.6, 66.6].map(Math.fround));
  i.meta.scans[0].intensityLimits.forEach((v, k) => close(v, [-66.6, 66.6][k], 1e-9));
});

test('real E57: pose with a translation only (pye57), normals extension skipped, images counted, colour', () => {
  const t = parseE57({ name: 'translationWithoutRotation.e57', bytes: real('e57-pye57-translation.e57') });
  assert.deepEqual([...t.models[0].positions], [4, 5, 6, 5, 5, 6, 4, 6, 6, 5, 6, 6]);
  assert.deepEqual(t.meta.scans[0].pose, { rotation: [1, 0, 0, 0], translation: [4, 5, 6] });
  const n = parseE57({ name: 'testWithNormals.e57', bytes: real('e57-normals.e57') });
  assert.deepEqual([...n.models[0].positions], [1, 0, 0, 2, 0, 0]);
  assert.deepEqual(n.meta.ignoredFields, { 'nor:normalX': 1, 'nor:normalY': 1, 'nor:normalZ': 1 });
  const im = parseE57({ name: 'tiny_pc_and_images.e57', bytes: real('e57-images.e57') });
  assert.deepEqual([...im.models[0].positions], [0, 0, 0, 1, 1, 1]);
  assert.equal(im.meta.skipped.images2D, 4);
  assert.deepEqual(im.warnings, ['4 image(s) (images2D) were not read.']);
  const c = parseE57({ name: 'tinyCartesianFloatRgb.e57', bytes: real('e57-cloudcompare-rgb.e57') });
  assertContract(c, ['points']);
  assert.equal(c.models[0].positions.length / 3, 2090);
  assert.deepEqual([...c.models[0].positions.slice(0, 3)], [-8.673629760742188, 4.575418472290039, 295.5494689941406]);
  assert.deepEqual([...c.models[0].attrs.rgb.slice(0, 6)], [219, 218, 222, 231, 232, 232]);
  assert.deepEqual(c.meta.bbox, { min: [-9.779529571533203, 4.513879299163818, 295.5246887207031], max: [-6.774238109588623, 7.51546049118042, 298.53216552734375] });
});

test('real E57, libE57Format test data: a page with a wrong checksum is reported', () => {
  assert.throws(() => parseE57({ name: 'bad-crc.e57', bytes: real('e57-bad-crc.e57') }), /bad-crc\.e57: 1 of the 1 pages of the XML section fail their CRC-32C checksum — the file is corrupt\./);
});
