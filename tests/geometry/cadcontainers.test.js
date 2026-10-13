import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { deflateSync } from 'node:zlib';
import { parseSolidWorks, parseSolidEdge, parseInventor } from '../../src/geometry/importers/cadcontainers.js';
import { parseParasolid } from '../../src/geometry/importers/parasolid.js';
import { isCfb, readCfb, cfbProgram } from '../../src/geometry/cfb.js';
import { importFile, formatFor, GeometryFormatError } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { signedVolume } from '../../src/mesh/measure.js';
import { enc, close, concat, assertContract, assertTruncationSafe } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const vol = (m) => signedVolume(m.positions, m.indices);
const u32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, v, true); return b; };

/**
 * Minimal compound file (version 3, 512-byte sectors): the streams `{ path: bytes }`, at most one level of
 * storages, streams under 4096 bytes in the mini stream.
 */
function makeCfb(streams) {
  const SEC = 512, END = 0xfffffffe, FREE = 0xffffffff, NONE = 0xffffffff;
  const entries = [{ name: 'Root Entry', type: 5, child: NONE, left: NONE, right: NONE, start: END, size: 0 }];
  const dirs = new Map();
  const link = (parent, id) => { if (parent.child === NONE) parent.child = id; else { let e = entries[parent.child]; while (e.right !== NONE) e = entries[e.right]; e.right = id; } };
  const big = [], small = [];
  for (const [path, data] of Object.entries(streams)) {
    const parts = path.split('/');
    let parent = entries[0];
    if (parts.length === 2) {
      if (!dirs.has(parts[0])) { entries.push({ name: parts[0], type: 1, child: NONE, left: NONE, right: NONE, start: 0, size: 0 }); dirs.set(parts[0], entries.length - 1); link(entries[0], entries.length - 1); }
      parent = entries[dirs.get(parts[0])];
    }
    const e = { name: parts[parts.length - 1], type: 2, child: NONE, left: NONE, right: NONE, start: END, size: data.length, data };
    entries.push(e);
    link(parent, entries.length - 1);
    if (data.length >= 4096) big.push(e); else if (data.length) small.push(e);
  }
  // mini stream
  const miniFat = [];
  const miniParts = [];
  for (const e of small) { const n = Math.ceil(e.size / 64); e.start = miniFat.length; for (let i = 0; i < n; i++) miniFat.push(i === n - 1 ? END : miniFat.length + 1); const pad = new Uint8Array(n * 64); pad.set(e.data); miniParts.push(pad); }
  const mini = concat(...miniParts);
  // sector plan: 0 FAT, then directory, mini FAT, mini stream, big streams
  const nDir = Math.ceil(entries.length / 4), nMiniFat = Math.max(1, Math.ceil(miniFat.length / 128)), nMini = Math.ceil(mini.length / SEC);
  const fat = [0xfffffffd];
  const chainOf = (n) => { const s = fat.length; for (let i = 0; i < n; i++) fat.push(i === n - 1 ? END : fat.length + 1); return n ? s : END; };
  const dirStart = chainOf(nDir), miniFatStart = chainOf(nMiniFat);
  entries[0].start = chainOf(nMini); entries[0].size = mini.length;
  for (const e of big) e.start = chainOf(Math.ceil(e.size / SEC));
  assert.ok(fat.length <= 128, 'test container fits one FAT sector');
  const out = new Uint8Array((fat.length + 1) * SEC);
  const dv = new DataView(out.buffer);
  out.set([0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1]);
  dv.setUint16(0x18, 0x3e, true); dv.setUint16(0x1a, 3, true); dv.setUint16(0x1c, 0xfffe, true); dv.setUint16(0x1e, 9, true); dv.setUint16(0x20, 6, true);
  dv.setUint32(0x2c, 1, true); dv.setUint32(0x30, dirStart, true); dv.setUint32(0x38, 4096, true); dv.setUint32(0x3c, miniFatStart, true); dv.setUint32(0x40, nMiniFat, true);
  dv.setUint32(0x44, END, true); dv.setUint32(0x48, 0, true);
  for (let i = 0; i < 109; i++) dv.setUint32(0x4c + 4 * i, i === 0 ? 0 : FREE, true);
  const at = (sec) => (sec + 1) * SEC;
  for (let i = 0; i < 128; i++) dv.setUint32(at(0) + 4 * i, i < fat.length ? fat[i] : FREE, true);
  entries.forEach((e, i) => {
    const o = at(dirStart) + 128 * i;
    for (let k = 0; k < e.name.length; k++) dv.setUint16(o + 2 * k, e.name.charCodeAt(k), true);
    dv.setUint16(o + 0x40, 2 * e.name.length + 2, true);
    out[o + 0x42] = e.type; out[o + 0x43] = 1;
    dv.setUint32(o + 0x44, e.left, true); dv.setUint32(o + 0x48, e.right, true); dv.setUint32(o + 0x4c, e.child, true);
    dv.setUint32(o + 0x74, e.start, true); dv.setUint32(o + 0x78, e.size, true);
  });
  for (let i = entries.length; i < nDir * 4; i++) { const o = at(dirStart) + 128 * i; dv.setUint32(o + 0x44, NONE, true); dv.setUint32(o + 0x48, NONE, true); dv.setUint32(o + 0x4c, NONE, true); }
  for (let i = 0; i < nMiniFat * 128; i++) dv.setUint32(at(miniFatStart) + 4 * i, i < miniFat.length ? miniFat[i] : FREE, true);
  if (nMini) out.set(mini, at(entries[0].start));
  for (const e of big) out.set(e.data, at(e.start));
  return out;
}

test('compound files: directory tree, regular and mini streams of a real SolidWorks part (compound-file generation)', () => {
  const bytes = real('solidworks-ole-washer.sldprt');
  assert.equal(isCfb(bytes), true);
  assert.equal(cfbProgram(bytes), 'solidworks');
  const cfb = readCfb(bytes, 'bearing.sldprt');
  const sizes = Object.fromEntries(cfb.entries.filter((e) => !e.storage).map((e) => [e.path, e.size]));
  // sizes as listed by an independent reader of the same file (olefile)
  assert.equal(Object.keys(sizes).length, 34);
  assert.deepEqual([sizes['Contents/Config-0-Partition'], sizes['Contents/Definition'], sizes.Header2, sizes.PreviewPNG], [2169, 3165, 1956, 6161]);
  assert.equal(cfb.entries.filter((e) => e.storage).map((e) => e.path).includes('Contents'), true);
  // a PNG in a regular stream, a small stream in the mini stream
  assert.deepEqual(Array.from(cfb.read('PreviewPNG').subarray(0, 8)), [0x89, 0x50, 0x4e, 0x47, 13, 10, 26, 10]);
  assert.equal(cfb.read('Contents/Config-0-Partition').length, 2169);
  assert.equal(cfb.read('no/such/stream'), null);
  // the test writer and the reader agree
  const mine = readCfb(makeCfb({ 'A/big': new Uint8Array(5000).fill(7), 'A/small': enc('hello'), top: enc('x'.repeat(200)), empty: new Uint8Array(0) }));
  assert.deepEqual(mine.entries.map((e) => `${e.path}:${e.size}`), ['A:0', 'A/big:5000', 'A/small:5', 'empty:0', 'top:200']);
  assert.deepEqual([mine.read('A/big').length, mine.read('A/big')[4999], new TextDecoder().decode(mine.read('A/small')), mine.read('top').length, mine.read('empty').length], [5000, 7, 'hello', 200, 0]);
});

test('SolidWorks part, compound file (written with Parasolid 24): a thrust washer 26 × 12 × 1 mm', async () => {
  const bytes = real('solidworks-ole-washer.sldprt');
  assert.equal(formatFor('AS1226.SLDPRT', bytes).id, 'solidworks');
  assert.equal(formatFor('renamed.bin', bytes).id, 'solidworks');
  const r = await importFile({ name: 'AS1226.SLDPRT', bytes }, { angularStep: 2 });
  assertContract(r, ['surface']);
  assert.equal(r.format, 'solidworks');
  assert.deepEqual(r.warnings, []);
  assert.deepEqual([r.meta.container, r.meta.kernelStream, r.meta.configuration, r.meta.configurations], ['SolidWorks compound file (2014 or earlier)', 'Contents/Config-0-Partition', 0, [0]]);
  assert.deepEqual([r.meta.schema, r.meta.encoding, r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.faces, r.meta.faceTotal], ['SCH_2401264_20000_13006', 'neutral binary', 1, 1, 1, 4, 4]);
  assert.deepEqual(r.meta.elementTypes, { cylinder: 2, plane: 2 });
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  const m = r.models[0];
  assert.equal(m.name, 'AS1226');
  assert.ok(inspect(m).watertight);
  const exact = Math.PI * (0.013 ** 2 - 0.006 ** 2) * 0.001;
  assert.ok(vol(m) < exact && vol(m) > exact * 0.999, `${vol(m)} against ${exact}`);
  [0.026, 0.001, 0.026].forEach((v, k) => close(r.meta.bbox.max[k] - r.meta.bbox.min[k], v, 1e-8));
});

test('SolidWorks part, record file (written with Parasolid 34): a sleeve 7.4 × 6.2 × 15 mm, as its STEP export describes it', async () => {
  const bytes = real('solidworks-record-sleeve.sldprt');
  assert.equal(isCfb(bytes), false);
  assert.equal(formatFor('sleeve.sldprt', bytes).id, 'solidworks');
  assert.equal(formatFor('renamed.bin', bytes).id, 'solidworks');
  const r = await importFile({ name: 'sleeve.sldprt', bytes }, { angularStep: 2 });
  assertContract(r, ['surface']);
  assert.equal(r.format, 'solidworks');
  assert.deepEqual([r.meta.container, r.meta.kernelStream, r.meta.schema], ['SolidWorks record file (2015 or later)', 'Contents/Config-0-Partition', 'SCH_3401254_34101_13006']);
  assert.deepEqual([r.meta.bodies, r.meta.closedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 4, 4]);
  // the STEP file of the same part in the same repository: cylinders of radius 3.7 and 3.1 mm, 15 mm long
  const exact = Math.PI * (0.0037 ** 2 - 0.0031 ** 2) * 0.015;
  const m = r.models[0];
  assert.ok(inspect(m).watertight);
  assert.ok(vol(m) < exact && vol(m) > exact * 0.999, `${vol(m)} against ${exact}`);
  [0.0074, 0.0074, 0.015].forEach((v, k) => close(r.meta.bbox.max[k] - r.meta.bbox.min[k], v, 1e-8));
  await assert.rejects(parseSolidWorks({ name: 's.sldprt', bytes }, { configuration: 3 }), /configuration: the file has body data for configuration 0, not for 3/);
});

test('Solid Edge part: the Parasolid partition stream of the compound file is read (container built here around a real Parasolid body)', async () => {
  const xb = real('solidworks-bell-crank.x_b');
  const par = makeCfb({ 'PARASOLID/STREAM0.P_B': xb, 'PARASOLID/STREAM0.D_B': new Uint8Array(4200), JSitesList: enc('OLEM\0\0\0\0'), BuildVersions: new Uint8Array(8) });
  assert.equal(cfbProgram(par), 'solid-edge');
  assert.equal(formatFor('crank.par', par).id, 'solid-edge');
  const r = await importFile({ name: 'crank.par', bytes: par });
  assert.equal(r.format, 'solid-edge');
  assert.deepEqual([r.meta.container, r.meta.kernelStream], ['Solid Edge compound file', 'PARASOLID/STREAM0.P_B']);
  const direct = parseParasolid({ name: 'crank.x_b', bytes: xb });
  assert.deepEqual([r.meta.faces, r.meta.faceTotal, r.meta.volumes, r.meta.bbox], [direct.meta.faces, direct.meta.faceTotal, direct.meta.volumes, direct.meta.bbox]);
  assert.deepEqual(r.warnings, direct.warnings);
  // two partitions in one file: both read, the models side by side
  const two = parseSolidEdge({ name: 'two.par', bytes: makeCfb({ 'PARASOLID/STREAM0.P_B': xb, 'PARASOLID/STREAM1.P_B': xb, JSitesList: enc('OLEM') }) });
  assert.deepEqual([two.models.length, two.meta.bodies, two.meta.faces, two.meta.kernelStream], [2, 2, 2 * direct.meta.faces, ['PARASOLID/STREAM0.P_B', 'PARASOLID/STREAM1.P_B']]);
  assert.throws(() => parseSolidEdge({ name: 'a.asm', bytes: makeCfb({ JSitesList: enc('OLEM'), BuildVersions: new Uint8Array(8) }) }), /a\.asm: a Solid Edge assembly stores references to its part files and no bodies of its own/);
  assert.throws(() => parseSolidEdge({ name: 'd.par', bytes: makeCfb({ JSitesList: enc('OLEM') }) }), /d\.par: the file holds no body data \(no "PARASOLID\/STREAMn\.P_B" stream\)/);
  assert.throws(() => parseSolidEdge({ name: 'w.par', bytes: makeCfb({ WordDocument: enc('x') }) }), /w\.par: this compound file has none of the streams of a Solid Edge document/);
  assert.throws(() => parseSolidEdge({ name: 't.par', bytes: enc('plain text, long enough to be told apart') }), /t\.par: not a Solid Edge file/);
});

test('Inventor part: the binary ACIS stream of the B-rep segment (segment streams of a real Inventor part, build 180, put back into a compound file)', async () => {
  // the two streams of the segment as they are in the part file; the other segments (browser, graphics …) are left out
  const B = real('inventor-pin.B.bin'), M = real('inventor-pin.M.bin');
  const ipt = makeCfb({ 'RSeStorage/Bwjduekaj5efbff5eenfnwjbxfb': B, 'RSeStorage/Mwjduekaj5efbff5eenfnwjbxfb': M, 'RSeStorage/RSeSegInfo': new Uint8Array(300), UFRxDoc: new Uint8Array(100) });
  assert.equal(cfbProgram(ipt), 'inventor');
  assert.equal(formatFor('pin.ipt', ipt).id, 'inventor');
  const r = await importFile({ name: 'pin.ipt', bytes: ipt }, { angularStep: 2 });
  assertContract(r, ['surface']);
  assert.equal(r.format, 'inventor');
  assert.deepEqual(r.warnings, []);
  assert.deepEqual([r.meta.container, r.meta.kernelStream, r.meta.segments], ['Inventor compound file', 'RSeStorage/Bwjduekaj5efbff5eenfnwjbxfb (PmBRepSegment)', ['PmBRepSegment']]);
  assert.deepEqual([r.meta.bodies, r.meta.closedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 8, 8]);
  assert.deepEqual(r.meta.elementTypes, { plane: 4, cylinder: 3, cone: 1 });
  // Inventor keeps lengths in centimetres; the model is in metres
  assert.deepEqual(r.meta.units, { unit: 'cm', scale: 0.01, basis: 'file' });
  const m = r.models[0];
  assert.ok(inspect(m).watertight);
  // a stepped pin 5.2 mm across and 15.125 mm long
  [0.0052, 0.0052, 0.015125].forEach((v, k) => close(r.meta.bbox.max[k] - r.meta.bbox.min[k], v, 1e-7));
  close(vol(m), 1.2624e-7, 2e-10);
  // Zstandard-compressed segments (Inventor 2026), assemblies, foreign files
  const z = B.slice(); z[17] = 2; z[18] = 0x28; z[19] = 0xb5; z[20] = 0x2f; z[21] = 0xfd;
  await assert.rejects(parseInventor({ name: 'new.ipt', bytes: makeCfb({ 'RSeStorage/Bx': z, 'RSeStorage/Mx': M }) }), /new\.ipt: the segments of this Inventor file are compressed with Zstandard, as Inventor 2026 writes them, which this reader does not decode\. Use Save Copy As → SAT \(ACIS text\) or STEP in Inventor; both are read\./);
  const dc = M.slice(); const at = Buffer.from(M).indexOf(Buffer.from('PmBRepSegment', 'utf16le')); dc.set(Buffer.from('AmBxxxSegment', 'utf16le'), at);
  await assert.rejects(parseInventor({ name: 'rig.iam', bytes: makeCfb({ 'RSeStorage/Bx': B, 'RSeStorage/Mx': dc }) }), /rig\.iam: an Inventor assembly, presentation or drawing stores references to its part files and no bodies of its own/);
  await assert.rejects(parseInventor({ name: 'w.ipt', bytes: makeCfb({ WordDocument: enc('x') }) }), /w\.ipt: this compound file has none of the segment streams of an Inventor document/);
  await assert.rejects(parseInventor({ name: 't.ipt', bytes: enc('no container at all, just text') }), /t\.ipt: not an Inventor file/);
  const cut = B.slice(0, 600);
  await assert.rejects(parseInventor({ name: 'd.ipt', bytes: makeCfb({ 'RSeStorage/Bx': cut, 'RSeStorage/Mx': M }) }), /d\.ipt: /);
  await assertTruncationSafe(ipt, (b) => parseInventor({ name: 'cut.ipt', bytes: b }), 'Inventor');
});

test('SolidWorks: assemblies, empty parts, foreign and damaged files are refused with the reason', async () => {
  const MARK = [0x23, 0x1d, 0xd5, 0x71, 0xda, 0x81, 0x48, 0xa2, 0xa8, 0x58, 0x98, 0xb2, 0x1b, 0x89, 0xef, 0x99];
  const block = (data) => { const z = new Uint8Array(deflateSync(data)); return concat(u32(16 + 8 + z.length + 8), Uint8Array.from(MARK), u32(data.length), u32(z.length), z, new Uint8Array(8)); };
  const xb = real('solidworks-bell-crank.x_b');
  const body = xb.subarray(xb.indexOf(0x50, 600)); // from the flag sequence "PS", as the streams of a part file hold it
  assert.deepEqual(Array.from(body.subarray(0, 4)), [0x50, 0x53, 0, 0]);
  // a container built here: the partition block is found, inflated and read
  const made = await parseSolidWorks({ name: 'made.sldprt', bytes: makeCfb({ 'Contents/Config-0-Partition': block(body), 'Contents/Config-2-Partition': block(body), Header2: new Uint8Array(10) }) });
  assert.deepEqual([made.meta.faces, made.meta.configuration, made.meta.configurations], [25, 0, [0, 2]]);
  assert.match(made.warnings.join(' '), /The file stores bodies for 2 configurations \(0, 2\); configuration 0 was read/);
  const asm = makeCfb({ 'Contents/Config-0': new Uint8Array(40), Header2: new Uint8Array(10) });
  await assert.rejects(importFile({ name: 'rig.SLDASM', bytes: asm }), /rig\.SLDASM: a SolidWorks assembly stores references to its part files and no bodies of its own\. Import the part files/);
  await assert.rejects(parseSolidWorks({ name: 'p.sldprt', bytes: asm }), /p\.sldprt: the file holds no body data \(no "Contents\/Config-n-Partition" stream\)/);
  await assert.rejects(parseSolidWorks({ name: 'x.sldprt', bytes: makeCfb({ WordDocument: enc('x') }) }), /x\.sldprt: this compound file has none of the streams of a SolidWorks document/);
  await assert.rejects(parseSolidWorks({ name: 'n.sldprt', bytes: new Uint8Array(600) }), /n\.sldprt: not a SolidWorks file — neither a compound \(OLE2\) file nor the record stream of SolidWorks 2015 and later was found/);
  await assert.rejects(parseSolidWorks({ name: 'd.sldprt', bytes: makeCfb({ 'Contents/Config-0-Partition': concat(u32(40), Uint8Array.from(MARK), u32(500), u32(12), Uint8Array.from([0x78, 1, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]), new Uint8Array(8)), Header2: new Uint8Array(10) }) }), /d\.sldprt: a compressed block of the body data could not be inflated — the file is damaged/);
  await assert.rejects(parseSolidWorks({ name: 'e.sldprt', bytes: makeCfb({ 'Contents/Config-0-Partition': new Uint8Array(64), Header2: new Uint8Array(10) }) }), /e\.sldprt: the body data of configuration 0 hold no Parasolid partition/);
  // formats that stay unread say what to export, each for its own program
  await assert.rejects(importFile({ name: 'house.rvt', bytes: makeCfb({ BasicFileInfo: new Uint8Array(8) }) }), (e) => e instanceof GeometryFormatError && e.reason === 'convert' && e.formatId === 'revit');
  for (const name of ['solidworks-ole-washer.sldprt', 'solidworks-record-sleeve.sldprt']) await assertTruncationSafe(real(name), (b) => parseSolidWorks({ name: 'cut.sldprt', bytes: b }), name);
  const par = makeCfb({ 'PARASOLID/STREAM0.P_B': xb, JSitesList: enc('OLEM') });
  await assertTruncationSafe(par, async (b) => parseSolidEdge({ name: 'cut.par', bytes: b }), 'Solid Edge');
  // damaged directory and allocation table of a compound file
  let seed = 11;
  const rnd = () => (seed = (seed * 1103515245 + 12345) >>> 0) >>> 24;
  const good = real('solidworks-ole-washer.sldprt');
  for (let k = 0; k < 60; k++) {
    const bad = good.slice();
    const at = k < 20 ? 0x1e + rnd() % 0x60 : 512 + (rnd() * 256 + rnd()) % 4096;
    for (let i = 0; i < 6; i++) bad[at + i] = rnd();
    try { const r = await parseSolidWorks({ name: 'bad.sldprt', bytes: bad }); assert.ok(Array.isArray(r.models)); } catch (e) { assert.ok(e instanceof Error && !(e instanceof TypeError) && !(e instanceof ReferenceError) && e.message.length > 15, `damage ${k}: ${e && e.stack}`); }
  }
});
