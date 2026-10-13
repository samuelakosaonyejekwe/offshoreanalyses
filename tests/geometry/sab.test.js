import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseSab, readSab, isSab } from '../../src/geometry/importers/sab.js';
import { parseAcis, readSat } from '../../src/geometry/importers/acis.js';
import { importFile, formatFor } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { enc, close, assertContract } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const vol = (m) => signedVolume(m.positions, m.indices);

// ---------------------------------------------------------------------------------------------------------
// A SAB writer for the tests: the records of a SAT file, value by value, in the tagged binary form. It is used
// to check the decoder against the text reader on the records of real SAT files — no body was available to
// us in both forms as written by a CAD system.
// ---------------------------------------------------------------------------------------------------------

const BOOL = { forward: 0, reversed: 1, single: 0, double: 1, out: 0, in: 1, forward_v: 0, reverse_v: 1, I: 0, F: 1 };
function toSab(text, { wide = false, signature = wide ? 'ASM BinaryFile8' : 'ACIS BinaryFile' } = {}) {
  const sat = readSat(text);
  const out = [];
  const dv = new DataView(new ArrayBuffer(8));
  const raw = (n, set) => { set(); for (let k = 0; k < n; k++) out.push(dv.getUint8(k)); };
  const int = (v) => (wide ? raw(8, () => dv.setBigInt64(0, BigInt(v), true)) : raw(4, () => dv.setInt32(0, v, true)));
  const f64 = (v) => raw(8, () => dv.setFloat64(0, v, true));
  const str = (tag, s) => { const b = enc(s); out.push(tag, b.length, ...b); };
  for (const c of signature) out.push(c.charCodeAt(0));
  const h = sat.header;
  int(h.version); int(0); int(h.bodies); int(h.flag);
  str(0x07, h.product); str(0x07, h.acis); str(0x07, h.date);
  out.push(0x06); f64(h.mmPerUnit); out.push(0x06); f64(h.resabs); out.push(0x06); f64(h.resnor);
  for (const r of sat.records) {
    const parts = r.type.split('-');
    parts.forEach((p, i) => str(i + 1 < parts.length ? 0x0e : 0x0d, p));
    if (r.type === 'transform') {
      // attribute pointer, id, then the matrix as text
      r.tok.slice(0, 2).forEach((t) => (t.t === 'p' ? (out.push(0x0c), int(t.v)) : (out.push(0x04), int(t.v))));
      const lit = enc(r.tok.slice(2).map((t) => t.raw ?? t.v).join(' ') + ' ');
      out.push(0x12); raw(4, () => dv.setUint32(0, lit.length, true)); out.push(...lit);
    } else for (const t of r.tok) {
      if (t.t === 'p') { out.push(0x0c); int(t.v); } else if (t.t === 'n') { if (/^-?\d+$/.test(t.raw)) { out.push(0x04); int(t.v); } else { out.push(0x06); f64(t.v); } } else if (t.t === '{') out.push(0x0f);
      else if (t.t === '}') out.push(0x10);
      else if (t.t === 's') str(0x07, t.v);
      else if (t.v in BOOL) out.push(BOOL[t.v] ? 0x0a : 0x0b);
      else str(0x07, t.v);
    }
    out.push(0x11);
  }
  for (const [tag, s] of [[0x0e, 'End'], [0x0e, 'of'], [0x0e, 'ACIS'], [0x0d, 'data']]) str(tag, s);
  return Uint8Array.from(out);
}

test('real SAB, AutoCAD (ASM 218): a flange plate with a bore and eight bolt holes', async () => {
  const bytes = real('autocad-flange.sab');
  assert.ok(isSab(bytes));
  const sab = readSab(bytes, 'flange.sab');
  assert.deepEqual([sab.header.version, sab.header.product, sab.header.mmPerUnit, sab.wide], [21800, 'Autodesk AutoCAD', 25.4, false]);
  assert.match(sab.header.acis, /^ASM \d+\.\d+\.\d+\.\d+ NT$/);
  assert.equal(sab.records[0].type, 'asmheader');
  // values keep their kind: pointers, numbers, logicals — and the keyword pairs of the text form stay logicals
  const face = sab.records.find((q) => q.type === 'face');
  assert.deepEqual(face.tok.slice(-2).map((t) => t.t), ['b', 'b']);
  assert.ok(sab.records.some((q) => q.type === 'cone-surface') && sab.records.some((q) => q.type === 'ellipse-curve'));
  const r = await importFile({ name: 'flange.sab', bytes });
  assert.equal(r.format, 'acis-sab');
  assertContract(r, ['surface']);
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.encoding, 'binary, 32-bit integers');
  assert.deepEqual(r.meta.units, { unit: 'in', scale: 0.0254, basis: 'file' });
  assert.deepEqual(r.meta.elementTypes, { cylinder: 9, plane: 2 });
  assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 1, 1, 11, 11]);
  const m = r.models[0];
  assert.ok(inspect(m).watertight);
  // a plate: 5 in thick, 115 in across; its volume is the area of one flat side times the thickness
  close(r.meta.bbox.max[2] - r.meta.bbox.min[2], 5 * 0.0254, 1e-7);
  close(r.meta.bbox.max[0] - r.meta.bbox.min[0], 115 * 0.0254, 1e-6);
  const flat = m.groups.map((g) => ({ g, sub: { positions: m.positions, indices: m.indices.subarray(g.start, g.start + g.count) } })).filter(({ sub }) => { let z0 = Infinity, z1 = -Infinity; for (const i of sub.indices) { z0 = Math.min(z0, m.positions[3 * i + 2]); z1 = Math.max(z1, m.positions[3 * i + 2]); } return z1 - z0 < 1e-9; });
  assert.equal(flat.length, 2);
  const a0 = surfaceArea(flat[0].sub.positions, flat[0].sub.indices), a1 = surfaceArea(flat[1].sub.positions, flat[1].sub.indices);
  close(a0, a1, 1e-9 * a0);
  close(vol(m), a0 * 5 * 0.0254, 1e-6 * vol(m), 'prism');
  // against the dimensions: outer radius 57.5 in, facets inscribed, the holes take a few per cent
  const disc = Math.PI * (57.5 * 0.0254) ** 2 * 5 * 0.0254;
  assert.ok(vol(m) < disc && vol(m) > 0.9 * disc, `${vol(m)} of ${disc}`);
});

test('real SAB, AutoCAD: a collar with two cylindrical faces', () => {
  const r = parseSab({ name: 'collar.sab', bytes: real('autocad-collar.sab') }, { angularStep: 5 });
  assert.deepEqual(r.warnings, []);
  assert.deepEqual(r.meta.elementTypes, { plane: 4, cylinder: 2 });
  assert.deepEqual([r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 6, 6]);
  assert.ok(inspect(r.models[0]).watertight);
  assert.ok(r.meta.volumes[0] > 0);
  close(vol(r.models[0]), r.meta.volumes[0], 1e-6 * r.meta.volumes[0]);
});

test('SAB decoder against the SAT reader: the records of real SAT files, re-encoded value by value, give the same geometry', () => {
  for (const name of ['fusion-box.sat', 'acis28-rod.sat', 'revit-solids.sat']) {
    const text = new TextDecoder().decode(real(name));
    const a = parseAcis({ name, text });
    for (const wide of [false, true]) {
      const bytes = toSab(text, { wide });
      const sab = readSab(bytes);
      assert.equal(sab.wide, wide);
      assert.equal(sab.records.length, readSat(text).records.length, name);
      const b = parseSab({ name, bytes });
      assert.equal(b.meta.format, 'acis-sab');
      assert.equal(b.meta.encoding, `binary, ${wide ? 64 : 32}-bit integers`);
      assert.equal(b.models.length, a.models.length);
      a.models.forEach((m, k) => { assert.deepEqual(b.models[k].positions, m.positions, `${name} body ${k + 1}`); assert.deepEqual(b.models[k].indices, m.indices); assert.deepEqual(b.models[k].groups, m.groups); });
      assert.deepEqual(b.meta.elementTypes, a.meta.elementTypes);
      assert.deepEqual(b.meta.units, a.meta.units);
      assert.deepEqual(b.meta.volumes, a.meta.volumes);
      assert.deepEqual(b.warnings, a.warnings);
    }
  }
});

test('SAB: damaged, truncated and foreign data are rejected with the reason', async () => {
  const good = real('autocad-collar.sab');
  assert.throws(() => readSab(enc('700 0 1 0\nbody $-1 #\n'), 'x.sab'), /x\.sab: not an ACIS binary \(SAB\) file: the signature "ACIS BinaryFile" \/ "ASM BinaryFile4" is missing\./);
  assert.equal(isSab(enc('ACIS BinaryFil')), false);
  for (const cut of [15, 20, 60, 400, good.length - 40]) assert.throws(() => parseSab({ name: 't.sab', bytes: good.slice(0, cut) }), /t\.sab: (the SAB file is truncated: it ends inside|the SAB header|the SAB file holds a header but no entity records)/, `cut at ${cut}`);
  // a tag the format does not have: the length of its value is unknown, so reading stops there
  const bad = good.slice();
  const at = bad.indexOf(0x11, 300) + 1;
  // corrupt the first value of a record (behind its name)
  let p = at + 2 + bad[at + 1];
  bad[p] = 0x7f;
  assert.throws(() => parseSab({ name: 'b.sab', bytes: bad }), /b\.sab: unknown value tag 0x7f at byte \d+ in record \d+ \([\w-]+\) of the SAB file/);
  // a header that is not laid out as expected
  const hdr = good.slice();
  hdr[31] = 0x06;
  assert.throws(() => parseSab({ name: 'h.sab', bytes: hdr }), /h\.sab: the SAB header is not laid out as expected/);
  // the text reader names the right reader when it is handed binary data, and the registry routes by content
  assert.throws(() => parseAcis({ name: 'b.sat', bytes: good }), /b\.sat: this is a binary ACIS file \(SAB\), which the SAT \(text\) reader does not take/);
  assert.equal(formatFor('part.sab').status, 'native');
  assert.equal(formatFor('part.sab').via, 'importers/sab.js');
  assert.equal((await importFile({ name: 'renamed.sat', bytes: good })).format, 'acis-sab');
});
