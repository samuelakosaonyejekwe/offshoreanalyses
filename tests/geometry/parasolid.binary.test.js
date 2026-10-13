import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseParasolid, readXt, readXb, isBinaryXt } from '../../src/geometry/importers/parasolid.js';
import { importFile, formatFor } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { signedVolume } from '../../src/mesh/measure.js';
import { close, concat, assertContract, assertTruncationSafe } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const vol = (m) => signedVolume(m.positions, m.indices);

// ---- neutral-binary encoder for the nodes of a text file, as the format reference describes the encoding ----
const bytesOf = (s) => Uint8Array.from(s, (c) => c.charCodeAt(0));
let LE = false; // bare binary of a little-endian machine instead of neutral binary
const i16 = (v) => { const b = new Uint8Array(2); new DataView(b.buffer).setInt16(0, v, LE); return b; };
const i32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setInt32(0, v, LE); return b; };
const f64 = (v) => { const b = new Uint8Array(8); new DataView(b.buffer).setFloat64(0, v === null ? -3.14158e13 : v, LE); return b; };
const ptr = (i) => (i < 32767 ? i16(i + 1) : concat(i16(-((i % 32767) + 1)), i16(Math.floor(i / 32767))));
const short = (s) => concat(Uint8Array.of(s.length), bytesOf(s));
function value(type, v) {
  switch (type) {
    case 'd': case 't': return i32(v === null ? -32764 : v);
    case 'n': case 'w': return i16(v);
    case 'u': return Uint8Array.of(v);
    case 'p': return ptr(v);
    case 'f': return f64(v);
    case 'c': return bytesOf(v);
    case 'l': return Uint8Array.of(v ? 1 : 0);
    case 'v': case 'h': return v === null ? concat(f64(null), f64(null), f64(null)) : concat(...v.map(f64));
    default: return concat(...v.map(f64)); // i, b
  }
}
/** Binary form of a text transmit file: same nodes, same schema (every node type of a file with an embedded schema written as a described one). */
function toBinary(text, { header = true, userField = 0, bare = false } = {}) {
  LE = bare;
  try { return encode(text, header, userField, bare); } finally { LE = false; }
}
function encode(text, header, userField, bare) {
  const xt = readXt(text);
  const data = text.slice(text.indexOf('\n', text.indexOf('**END_OF_HEADER')) + 1).replace(/[\r\n]/g, '');
  const m = /^T(\d+) /.exec(data), mod = data.substr(m[0].length, Number(m[1]));
  const embedded = xt.header.baseSchema !== null;
  const out = [];
  if (header) out.push(bytesOf(text.slice(0, text.indexOf('\n', text.indexOf('**END_OF_HEADER')) + 1).replace('FORMAT=text', 'FORMAT=binary')));
  out.push(bare ? Uint8Array.of(0x42) : Uint8Array.of(0x50, 0x53, 0, 0), bare ? i32(mod.length) : i16(mod.length), bytesOf(mod), i32(xt.header.schema.length), bytesOf(xt.header.schema));
  if (embedded) out.push(i16(250));
  out.push(i32(userField));
  const seen = new Set();
  for (const node of xt.nodes.values()) {
    const fields = xt.schema.get(node.type);
    out.push(i16(node.type));
    if (embedded && !seen.has(node.type)) {
      seen.add(node.type);
      // described in full: drop every field of the base schema (if it has the type), then insert the fields in force
      out.push(Uint8Array.of(fields.length));
      const named = xt.typeNames.get(node.type);
      if (named !== undefined) out.push(short(named), short('described'));
      else out.push(bytesOf('D'.repeat(40)));
      for (const f of fields) {
        if (named === undefined) out.push(bytesOf('I'));
        out.push(short(f.name), i16(f.type === 'p' ? 1012 : 0), ptr(f.n));
        if (f.type !== 'p') out.push(short(f.type));
        if (f.n === 1) out.push(Uint8Array.of(f.skip ? 0 : 1));
      }
      if (named === undefined) out.push(bytesOf('Z'));
    }
    const vf = fields.find((f) => f.n === 1);
    const stored = (f) => node[f.name === 'type' ? 'kind' : f.name === 'index' || f.name === 'name' ? `${f.name}_` : f.name];
    if (vf) out.push(i32(vf.skip ? 0 : stored(vf).length));
    out.push(ptr(node.index));
    for (const f of fields) {
      if (f.skip) continue;
      const v = stored(f);
      if (f.n === 0) out.push(value(f.type, v));
      else if (f.type === 'c') out.push(bytesOf(v));
      else for (const x of v) out.push(value(f.type, x));
    }
    for (let k = 0; k < userField; k++) out.push(i32(0));
  }
  out.push(i16(1), ptr(0));
  return concat(...out);
}

test('Parasolid binary: file written by SolidWorks (schema 19008) — same body as the SolidWorks part it was exported from', async () => {
  const bytes = real('solidworks-bell-crank.x_b');
  assert.equal(isBinaryXt(bytes), true);
  assert.equal(formatFor('Bell_Crank.x_b', bytes).id, 'parasolid');
  const xb = readXb(bytes, 'Bell_Crank.x_b');
  assert.deepEqual(xb.header, { modeller: 1901315, schema: 'SCH_1901315_19008_13006', schemaNumber: 19008, baseSchema: 13006, application: 'SolidWorks 2009-2009077', date: 'Mon Sep 21 17:12:13 2009', format: 'binary', partition: false });
  assert.equal(xb.layout, 'described by the file');
  assert.equal(xb.end, bytes.length);
  assert.deepEqual([xb.types.BODY, xb.types.FACE, xb.types.PLANE, xb.types.CYLINDER, xb.types.BLENDED_EDGE], [1, 27, 11, 17, 4]);
  const r = await importFile({ name: 'Bell_Crank.x_b', bytes }, { angularStep: 5 });
  assertContract(r, ['surface']);
  assert.equal(r.format, 'parasolid');
  assert.equal(r.meta.encoding, 'neutral binary');
  assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 25, 27]);
  // the two rolling-ball blends of the part are kept by Parasolid as a definition, which is not read
  assert.deepEqual(r.meta.skipped, { 'face: BLENDED_EDGE surface, which is not read': 2 });
  assert.match(r.warnings[0], /^2 of 27 faces could not be read \(2 × BLENDED_EDGE surface, which is not read\)/);
  assert.deepEqual(r.meta.elementTypes, { torus: 1, plane: 11, cylinder: 13 });
  // the extent of the part: 41.2 mm × 23.4 mm × 8.4 mm, as the SolidWorks part file of the same crank gives it
  const size = [0, 1, 2].map((k) => r.meta.bbox.max[k] - r.meta.bbox.min[k]);
  close(size[0], 0.0412213, 2e-6); close(size[1], 0.0233959, 2e-6); close(size[2], 0.0083566, 2e-6);
});

test('Parasolid binary: the nodes of real text files, encoded as neutral binary, give the same geometry', () => {
  for (const name of ['onshape-sphere-in-box.x_t', 'solidworks-plate.x_t', 'solidworks-load-cell.x_t', 'solidworks-bars-v10.x_t']) {
    const text = new TextDecoder().decode(real(name));
    const a = parseParasolid({ name, text });
    for (const o of [{}, { header: false }, { userField: 2 }, { bare: true }]) {
      const bin = toBinary(text, o);
      const b = parseParasolid({ name, bytes: bin });
      assert.equal(b.meta.encoding, o.bare ? 'bare binary (little-endian)' : 'neutral binary');
      assert.deepEqual(b.meta.nodeTypes, a.meta.nodeTypes, name);
      assert.deepEqual([b.meta.schema, b.meta.bodies, b.meta.closedBodies, b.meta.faces, b.meta.elementTypes], [a.meta.schema, a.meta.bodies, a.meta.closedBodies, a.meta.faces, a.meta.elementTypes], name);
      // text holds the shortest decimal form of every double, so the round trip is exact
      assert.deepEqual(b.meta.volumes, a.meta.volumes, name);
      assert.deepEqual(b.meta.bbox, a.meta.bbox, name);
      for (const m of b.models) assert.ok(inspect(m).watertight, name);
      close(b.models.reduce((s, m) => s + vol(m), 0), a.models.reduce((s, m) => s + vol(m), 0), 1e-15);
    }
  }
  // the file of a release that does not describe its schema is read with the same trial of layouts
  const old = parseParasolid({ name: 'bars.x_b', bytes: toBinary(new TextDecoder().decode(real('solidworks-bars-v10.x_t'))) });
  assert.equal(old.meta.schemaLayout, 'old-attdef, bare-fin');
});

test('Parasolid binary: pointer indices from 32767 on are written as a pair of shorts', () => {
  // a point node far up the index range, referenced from nowhere: the stream must still read to its terminator
  const text = new TextDecoder().decode(real('solidworks-plate.x_t'));
  const bin = toBinary(text);
  const tail = concat(i16(29), ptr(70000), i32(1), ptr(0), ptr(0), ptr(0), ptr(0), f64(1), f64(2), f64(3), i16(1), ptr(0));
  const xb = readXb(concat(bin.subarray(0, bin.length - 4), tail), 'big.x_b');
  assert.deepEqual(xb.nodes.get(70000).pvec, [1, 2, 3]);
  assert.deepEqual(Array.from(ptr(70000)), [0xee, 0x8d, 0, 2]); // −(4466 + 1), 2
});

test('Parasolid binary: bare binary, damaged and truncated files are refused with the reason', async () => {
  const text = new TextDecoder().decode(real('solidworks-plate.x_t'));
  const bin = toBinary(text), bare = toBinary(text, { header: false });
  const P = (bytes, name = 'p.x_b') => parseParasolid({ name, bytes });
  const head = bin.length - bare.length;
  // bare binary of a big-endian machine: the flag "B", then the length of the modeller version as 0, 0, 0, n
  const b2 = concat(bin.subarray(0, head), Uint8Array.of(0x42, 0, 0, 0, 0x33), bin.subarray(head + 6));
  assert.throws(() => P(b2), /p\.x_b: this is a "bare binary" Parasolid file written in the number format of a machine that is not little-endian; only neutral binary \(flag "PS"\), little-endian bare binary and text files are read\. Save the part again as Parasolid text \(\.x_t\), or export STEP\./);
  assert.throws(() => readXb(Uint8Array.of(1, 2, 3, 4, 5, 6, 7, 8), 'n.x_b'), /n\.x_b: not a Parasolid binary transmit file \(the data do not start with the flag sequence "PS", 0, 0\)/);
  assert.throws(() => P(bare.subarray(0, bare.length - 4)), /p\.x_b: the file is truncated: the node stream ends after \d+ nodes without its terminator/);
  assert.throws(() => P(bare.subarray(0, Math.floor(bare.length / 2))), /p\.x_b: the node stream could not be interpreted \(schema SCH_1901315_19008_13006\): after \d+ nodes, the file ends inside /);
  // a node type the file does not describe
  const odd = bare.slice(); odd[bare.length - 4] = 0x7f; odd[bare.length - 3] = 0x7f;
  assert.throws(() => P(odd), /after \d+ nodes, (the file ends inside|.* is not a node type|an index was expected)/);
  // a pointer of zero (indices are stored plus one)
  const at = bare.length - 4 - 40;
  const zero = bare.slice(); zero.fill(0, at, at + 36);
  assert.throws(() => P(zero), /p\.x_b: the node stream could not be interpreted/);
  await assertTruncationSafe(real('solidworks-bell-crank.x_b'), (b) => P(b), 'x_b');
  await assertTruncationSafe(bare, (b) => P(b), 'x_b without header');
  // random bytes after a valid flag sequence
  let seed = 7;
  const rnd = () => (seed = (seed * 1103515245 + 12345) >>> 0) >>> 24;
  for (let k = 0; k < 40; k++) {
    const junk = Uint8Array.from({ length: 300 }, rnd);
    const g = concat(bare.subarray(0, 100 + 7 * k), junk);
    assert.throws(() => P(g), (e) => e instanceof Error && !(e instanceof TypeError) && !(e instanceof RangeError) && /^p\.x_b: /.test(e.message), `junk ${k}`);
  }
});
