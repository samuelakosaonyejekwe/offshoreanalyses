import test from 'node:test';
import assert from 'node:assert/strict';
import { parseDxf } from '../../src/geometry/importers/dxf.js';
import { parseDxfBinary, readBinaryDxfPairs, groupType } from '../../src/geometry/importers/dxfbinary.js';
import { importFile, formatFor } from '../../src/geometry/registry.js';
import { enc, assertContract } from './fixtures/helpers.js';

// ---- fixture: one list of (group code, value) pairs, written out in the three encodings -------------------

const SENTINEL = [...'AutoCAD Binary DXF\r\n\x1a\0'].map((c) => c.charCodeAt(0));
const ascii = (pairs) => { let s = ''; for (let i = 0; i < pairs.length; i += 2) s += `${String(pairs[i]).padStart(3)}\n${pairs[i + 1]}\n`; return s; };

/** Binary DXF from group pairs. `wide`: 16-bit group codes (R13 and later); otherwise one byte, 255 = escape. */
function binary(pairs, wide) {
  const out = [...SENTINEL];
  const dv = new DataView(new ArrayBuffer(8));
  const put = (n, set) => { set(); for (let k = 0; k < n; k++) out.push(dv.getUint8(k)); };
  for (let i = 0; i < pairs.length; i += 2) {
    const code = pairs[i], v = pairs[i + 1];
    if (wide) put(2, () => dv.setInt16(0, code, true));
    else if (code >= 255) { out.push(255); put(2, () => dv.setInt16(0, code, true)); } else out.push(code);
    const t = groupType(code);
    if (t === 0) { for (const b of enc(String(v))) out.push(b); out.push(0); }
    else if (t === 1) put(8, () => dv.setFloat64(0, Number(v), true));
    else if (t === 2) put(2, () => dv.setInt16(0, Number(v), true));
    else if (t === 3) put(4, () => dv.setInt32(0, Number(v), true));
    else if (t === 4) put(8, () => dv.setBigInt64(0, BigInt(v), true));
    else if (t === 5) out.push(Number(v) ? 1 : 0);
    else if (t === 6) { const hex = String(v); out.push(hex.length / 2); for (let k = 0; k < hex.length; k += 2) out.push(parseInt(hex.slice(k, k + 2), 16)); }
    else throw new Error(`fixture: group code ${code}`);
  }
  return Uint8Array.from(out);
}

const LINE = (layer, a, b) => [0, 'LINE', 8, layer, 10, a[0], 20, a[1], 30, a[2], 11, b[0], 21, b[1], 31, b[2]];
const face = (v, layer) => [0, '3DFACE', 8, layer, ...v.flatMap((p, k) => [10 + k, p[0], 20 + k, p[1], 30 + k, p[2]])];
const file = ({ header = [], tables = [], blocks = [], entities }) => [0, 'SECTION', 2, 'HEADER', ...header, 0, 'ENDSEC', 0, 'SECTION', 2, 'TABLES', ...tables, 0, 'ENDSEC', 0, 'SECTION', 2, 'BLOCKS', ...blocks, 0, 'ENDSEC', 0, 'SECTION', 2, 'ENTITIES', ...entities, 0, 'ENDSEC', 0, 'EOF'];

// every entity kind the ASCII reader handles, with layers, a block, units and extended data
const DRAWING = file({
  header: [9, '$ACADVER', 1, 'AC1009', 9, '$INSUNITS', 70, 4],
  tables: [0, 'TABLE', 2, 'LAYER', 70, 3, 0, 'LAYER', 2, 'LEGS', 70, 0, 62, 7, 6, 'CONTINUOUS', 0, 'LAYER', 2, 'SKIN', 70, 0, 62, 3, 6, 'CONTINUOUS', 0, 'ENDTAB'],
  blocks: [0, 'BLOCK', 8, '0', 2, 'PAD', 70, 0, 10, 0, 20, 0, 30, 0, ...face([[0, 0, 0], [1000, 0, 0], [1000, 1000, 0], [0, 1000, 0]], '0'), 0, 'ENDBLK'],
  entities: [
    ...LINE('LEGS', [0, 0, 0], [0, 0, 10000.125]), ...LINE('LEGS', [8000, 0, 0], [8000, 0, 10000.125]), ...LINE('BRACES', [0, 0, 10000.125], [8000, 0, 10000.125]),
    ...face([[0, 0, 0], [1000, 0, 0], [1000, 1000, 500], [1000, 1000, 500]], 'SKIN'),
    0, 'POLYLINE', 8, 'MESH', 66, 1, 70, 64, 71, 4, 72, 1,
    ...[[0, 0, 5000], [1000, 0, 5000], [1000, 1000, 5000], [0, 1000, 5000]].flatMap((p) => [0, 'VERTEX', 8, 'MESH', 10, p[0], 20, p[1], 30, p[2], 70, 192]),
    0, 'VERTEX', 8, 'MESH', 10, 0, 20, 0, 30, 0, 70, 128, 71, 1, 72, 2, 73, 3, 74, -4, 0, 'SEQEND',
    0, 'INSERT', 8, 'PADS', 2, 'PAD', 10, 10000, 20, 0, 30, 0, 41, 2, 42, 3, 43, 1, 50, 90,
    0, 'LWPOLYLINE', 8, 'PLATE', 90, 4, 70, 1, 38, 250, 10, 0, 20, 0, 42, 0.41421356237309503, 10, 400, 20, 0, 10, 400, 20, 300, 10, 0, 20, 300,
    0, 'POLYLINE', 8, 'ROUTE', 66, 1, 70, 8, 0, 'VERTEX', 8, 'ROUTE', 10, 0, 20, 0, 30, 0, 70, 32, 0, 'VERTEX', 8, 'ROUTE', 10, 500, 20, 100, 30, -20, 70, 32, 0, 'VERTEX', 8, 'ROUTE', 10, 900, 20, 400, 30, -35.5, 70, 32, 0, 'SEQEND',
    0, 'CIRCLE', 8, 'C', 10, 5000, 20, 0, 30, 0, 40, 1000, 210, 0, 220, 0, 230, -1,
    0, 'ARC', 8, 'A', 10, 0, 20, 7000, 30, 100, 40, 300, 50, 10, 51, 200,
    0, 'ELLIPSE', 8, 'E', 10, 3000, 20, 7000, 30, 0, 11, 600, 21, 0, 31, 0, 40, 0.5, 41, 0, 42, 6.283185307179586,
    0, 'SPLINE', 8, 'S', 70, 4, 71, 2, 72, 6, 73, 3, 40, 0, 40, 0, 40, 0, 40, 1, 40, 1, 40, 1, 10, 1000, 20, 0, 30, 9000, 10, 1000, 20, 1000, 30, 9000, 10, 0, 20, 1000, 30, 9000, 41, 1, 41, Math.SQRT1_2, 41, 1,
    0, 'POINT', 8, 'P', 10, 1.5, 20, -2.25, 30, 1e-3,
    // groups whose values are not strings or doubles: boolean, 32-bit, binary chunk, extended data
    0, 'TEXT', 8, 'NOTES', 10, 0, 20, 0, 30, 0, 40, 2.5, 1, 'north arrow', 290, 1, 90, 70000, 310, '00FF10A5', 1001, 'APP', 1000, 'tag', 1040, 1.25, 1070, 7, 1071, 123456, 1004, 'CAFE',
  ],
});

const plain = (r) => JSON.parse(JSON.stringify(r.models, (k, v) => (ArrayBuffer.isView(v) ? { type: v.constructor.name, data: Array.from(v) } : v)));

test('binary DXF: both group-code widths decode to the pairs of the ASCII form and give identical models', () => {
  const ref = parseDxf({ name: 'drawing.dxf', text: ascii(DRAWING) });
  assertContract(ref, ['surface', 'section', 'frame', 'polyline', 'points']);
  for (const wide of [true, false]) {
    const bytes = binary(DRAWING, wide);
    const { pairs, codeBytes } = readBinaryDxfPairs({ name: 'drawing.dxf', bytes });
    assert.equal(codeBytes, wide ? 2 : 1);
    assert.equal(pairs.length, DRAWING.length);
    for (let i = 0; i < pairs.length; i += 2) {
      assert.equal(pairs[i], DRAWING[i], `code of pair ${i / 2}`);
      assert.equal(pairs[i + 1], String(DRAWING[i + 1]), `value of pair ${i / 2} (group ${pairs[i]})`);
    }
    const r = parseDxfBinary({ name: 'drawing.dxf', bytes });
    assertContract(r, ['surface', 'section', 'frame', 'polyline', 'points']);
    assert.deepEqual(plain(r), plain(ref));
    assert.deepEqual(r.warnings, ref.warnings);
    assert.deepEqual(r.meta.elementTypes, ref.meta.elementTypes);
    assert.deepEqual(r.meta.skipped, ref.meta.skipped);
    assert.deepEqual(r.meta.layers, ref.meta.layers);
    assert.deepEqual(r.meta.bbox, ref.meta.bbox);
    assert.equal(r.meta.units.unit, 'mm');
    assert.equal(r.meta.units.basis, 'file');
    assert.equal(r.meta.format, 'dxf-binary');
    assert.match(r.meta.encoding, wide ? /16-bit group codes/ : /single-byte group codes/);
  }
  // what the comparison rests on: the fixture really exercises the entity kinds
  for (const k of ['LINE', '3DFACE', 'POLYFACE', 'INSERT', 'LWPOLYLINE', 'POLYLINE', 'CIRCLE', 'ARC', 'ELLIPSE', 'SPLINE', 'POINT']) assert.ok(ref.meta.elementTypes[k] >= 1, k);
  assert.equal(ref.meta.blocks, 1);
  assert.deepEqual(ref.meta.layers, ['LEGS', 'SKIN']);
});

test('binary DXF: read through the registry; an option such as dxfLines reaches the shared entity handling', async () => {
  const sq = file({ entities: [...LINE('P', [0, 0, 0], [2, 0, 0]), ...LINE('P', [2, 0, 0], [2, 1, 0]), ...LINE('P', [2, 1, 0], [0, 1, 0]), ...LINE('P', [0, 1, 0], [0, 0, 0])] });
  const bytes = binary(sq, true);
  const f = formatFor('outline.dxf', bytes);
  assert.equal(f.id, 'dxf-binary');
  assert.equal(f.status, 'native');
  const r = await importFile({ name: 'outline.dxf', bytes });
  assert.equal(r.format, 'dxf-binary');
  assert.deepEqual(r.models.map((m) => m.kind), ['section']);
  assert.deepEqual(Array.from(r.models[0].loops[0]), [0, 0, 2, 0, 2, 1, 0, 1]);
  const fr = await importFile({ name: 'outline.dxf', bytes: binary(sq, false) }, { dxfLines: 'frame' });
  assert.deepEqual(fr.models.map((m) => m.kind), ['frame']);
  assert.equal(fr.models[0].members.length, 4);
});

test('binary DXF: bad sentinel, truncated files and unknown group codes are rejected with the reason', () => {
  const good = binary(DRAWING, true);
  const bad = good.slice();
  bad[3] = 0x58;
  assert.throws(() => parseDxfBinary({ name: 'b.dxf', bytes: bad }), /b\.dxf: not a binary DXF file \(the 22-byte sentinel/);
  assert.throws(() => parseDxfBinary({ name: 'b.dxf', bytes: enc('  0\nSECTION\n  2\nENTITIES\n') }), /not a binary DXF file/);
  assert.throws(() => parseDxfBinary({ name: 'b.dxf', bytes: good.subarray(0, 22) }), /truncated: it ends right after the sentinel/);
  assert.throws(() => parseDxfBinary({ name: 'b.dxf', bytes: good.subarray(0, 27) }), /truncated: it ends inside the first groups/);
  // every cut short of the full file is reported as truncation, never read as a shorter drawing
  for (const wide of [true, false]) {
    const b = binary(DRAWING, wide);
    for (const cut of [60, 200, 333, b.length >> 1, b.length - 9, b.length - 4, b.length - 1]) assert.throws(() => parseDxfBinary({ name: 't.dxf', bytes: b.subarray(0, cut) }), /t\.dxf: the binary DXF file is truncated: /, `cut at ${cut}`);
  }
  // a cut on a pair boundary leaves well-formed pairs but no EOF group
  const noEof = binary(DRAWING.slice(0, -2), true);
  assert.throws(() => parseDxfBinary({ name: 't.dxf', bytes: noEof }), /truncated: it ends after \d+ groups without the closing "0 EOF" group/);
  const odd = binary(file({ entities: LINE('L', [0, 0, 0], [1, 0, 0]) }), true);
  // replace the group code of the "0 ENDSEC" that closes ENTITIES by a code outside the tables
  let hit = -1;
  for (let i = 22; i + 9 < odd.length; i++) if (odd[i] === 0 && odd[i + 1] === 0 && String.fromCharCode(...odd.subarray(i + 2, i + 8)) === 'ENDSEC') hit = i;
  assert.ok(hit > 22);
  new DataView(odd.buffer).setInt16(hit, 5000, true);
  assert.throws(() => parseDxfBinary({ name: 'u.dxf', bytes: odd }), /u\.dxf: binary DXF could not be decoded: group code 5000 at byte \d+ .* is not in the DXF group-code tables/);
  assert.throws(() => parseDxfBinary({ name: 'e.dxf', bytes: binary(file({ entities: [] }), true) }), /no ENTITIES or BLOCKS section with content/);
});
