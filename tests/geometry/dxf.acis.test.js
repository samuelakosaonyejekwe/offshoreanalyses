import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseDxf, decodeAcisText } from '../../src/geometry/importers/dxf.js';
import { parseDxfBinary, groupType } from '../../src/geometry/importers/dxfbinary.js';
import { importFile } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { enc, close, assertContract } from './fixtures/helpers.js';

// Solids stored as ACIS data inside DXF drawings: real drawings (fixtures/real/README.md) for the two ways the
// data are stored, and drawings assembled from their entities for blocks, mirroring and the binary encoding.

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const text = (name) => new TextDecoder().decode(real(name));
const vol = (m) => signedVolume(m.positions, m.indices);

test('ACIS text in DXF: the character substitution of groups 1 and 3', () => {
  // first lines of a real R2000 drawing: a space stays, every other character c becomes 159 − c
  assert.equal(decodeAcisText('koo o n o '), '400 0 1 0 ');
  assert.equal(decodeAcisText('=0;& {n {m {rn {l |'), 'body $1 $2 $-1 $3 #');
  assert.equal(decodeAcisText('+0-*,r,*-9><: {rn o o o o o n lnqffffffffffffffl no n o o 90-(>-;@) V V V V |'), 'torus-surface $-1 0 0 0 0 0 1 31.999999999999993 10 1 0 0 forward_v I I I I #');
  // "A" is written as "^" followed by a blank that is not part of the text
  assert.equal(decodeAcisText('ll P/:1 [:,681 ^ 336>1<: ^ \\VL ]*63;:- nm ^ \\VL lmqo QK'), '33 Open Design Alliance ACIS Builder 12 ACIS 32.0 NT');
  assert.equal(decodeAcisText('@'), '_');
  assert.equal(decodeAcisText('_'), '@');
  // the substitution is its own inverse apart from that blank
  const plain = 'cone-surface $-1 -1 $-1 0 0 0 0 0 1 10 0 0 1 I I 0 1 10 forward I I I I #';
  assert.equal(decodeAcisText(decodeAcisText(plain)), plain);
});

test('real DXF R2010 with 3DSOLID entities stored as SAT text: four solids, one surface each', async () => {
  const r = await importFile({ name: 'solids.dxf', bytes: real('ezdxf-3dsolids-r2010.dxf') }, { angularStep: 5 });
  assert.equal(r.format, 'dxf');
  assertContract(r, ['surface']);
  assert.equal(r.meta.acadVersion, 'AC1024');
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.deepEqual(r.meta.elementTypes, { '3DSOLID': 4 });
  assert.deepEqual(r.meta.acis, { entities: 4, text: 4, binary: 0, bodies: 4, solidBodies: 4, closedBodies: 4, faces: 32, facesRead: 32, encoding: 'SAT text (groups 1 / 3)' });
  assert.deepEqual(r.warnings, []);
  assert.deepEqual(r.models.map((m) => [m.name, m.kind, m.entity, m.layer, m.groups.length]), [['0 3DSOLID 1', 'surface', '3DSOLID', '0', 9], ['0 3DSOLID 2', 'surface', '3DSOLID', '0', 8], ['0 3DSOLID 3', 'surface', '3DSOLID', '0', 8], ['0 3DSOLID 4', 'surface', '3DSOLID', '0', 7]]);
  for (const m of r.models) { assert.ok(inspect(m).watertight, m.name); assert.ok(vol(m) > 0); }
  // the first is a 10 mm cube with a 5 mm cube cut from one corner: volume 875 mm³, and the surface of the whole cube
  close(vol(r.models[0]), 875e-9, 1e-13);
  close(surfaceArea(r.models[0].positions, r.models[0].indices), 600e-6, 1e-10);
  // the others are the cube with rounded or cut edges: between 85 % and 100 % of it, each in its own 10 mm cell
  r.models.forEach((m, k) => {
    assert.ok(vol(m) > 0.85e-6 && vol(m) < 1e-6 + 1e-12, `${m.name}: ${vol(m)}`);
    let x0 = Infinity, x1 = -Infinity;
    for (let i = 0; i < m.positions.length; i += 3) { x0 = Math.min(x0, m.positions[i]); x1 = Math.max(x1, m.positions[i]); }
    close(x0, 0.02 * k, 1e-8); close(x1, 0.02 * k + 0.01, 1e-8);
  });
});

test('real DXF R2018 with a 3DSOLID stored as SAB in ACDSDATA: a sphere of radius 5 m', async () => {
  const R = 5;
  let prev = Infinity;
  for (const angularStep of [15, 5]) {
    const r = await importFile({ name: 'sphere.dxf', bytes: real('ezdxf-sphere-r2018.dxf') }, { angularStep });
    assert.equal(r.meta.acadVersion, 'AC1032');
    assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
    assert.deepEqual(r.meta.acis, { entities: 1, text: 0, binary: 1, bodies: 1, solidBodies: 1, closedBodies: 1, faces: 1, facesRead: 1, encoding: 'SAB binary (ACDSDATA)' });
    const m = r.models.find((q) => q.entity === '3DSOLID');
    assert.ok(inspect(m).watertight);
    for (let i = 0; i < m.positions.length; i += 3) close(Math.hypot(m.positions[i], m.positions[i + 1], m.positions[i + 2]), R, 1e-5);
    const err = 1 - vol(m) / ((4 / 3) * Math.PI * R ** 3);
    assert.ok(err > 0 && err < prev, 'inscribed, converging');
    prev = err;
  }
  assert.ok(prev < 3e-3);
  // without the record that holds the data the solid is reported, not invented
  const noData = text('ezdxf-sphere-r2018.dxf').replace(/ASM_Data/g, 'ASM_Gone');
  await assert.rejects(importFile({ name: 'empty.dxf', bytes: enc(noData) }), /empty\.dxf: the drawing contains no geometry this reader keeps \(found only: 3DSOLID \(embedded ACIS geometry: no ACIS data stored with the entity\)\)/);
});

/** The group pairs of the 3DSOLID entities of a drawing text. */
function solids(dxf) {
  const lines = dxf.split(/\r?\n/);
  const pairs = [];
  for (let i = 0; i + 1 < lines.length; i += 2) pairs.push([parseInt(lines[i], 10), lines[i + 1]]);
  const out = [];
  pairs.forEach((p, a) => { if (p[0] === 0 && p[1].trim() === '3DSOLID') { let b = a + 1; while (pairs[b][0] !== 0) b++; out.push(pairs.slice(a, b)); } });
  return out;
}
/** The notched cube of the R2010 fixture: the solid with nine faces. */
function solidPairs(dxf) {
  for (const s of solids(dxf)) {
    const r = parseDxf({ name: 'one.dxf', text: `  0\nSECTION\n  2\nENTITIES\n${asText(s)}  0\nENDSEC\n  0\nEOF\n` });
    if (r.models[0].groups.length === 9) return s;
  }
  throw new Error('fixture: no solid with nine faces');
}
const asText = (pairs) => pairs.map(([c, v]) => `${String(c).padStart(3)}\n${v}\n`).join('');

test('3DSOLID inside a block: placed by its INSERT, mirrored copies keep an outward orientation', () => {
  const solid = solidPairs(text('ezdxf-3dsolids-r2010.dxf'));
  const dxf = (inserts) => `  0\nSECTION\n  2\nHEADER\n  9\n$INSUNITS\n 70\n4\n  0\nENDSEC\n  0\nSECTION\n  2\nBLOCKS\n  0\nBLOCK\n  2\nPART\n 10\n0\n 20\n0\n 30\n0\n${asText(solid)}  0\nENDBLK\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n${inserts}  0\nENDSEC\n  0\nEOF\n`;
  const ins = (x, sx) => `  0\nINSERT\n  8\nPARTS\n  2\nPART\n 10\n${x}\n 20\n0\n 30\n0\n 41\n${sx}\n 42\n1\n 43\n1\n`;
  const r = parseDxf({ name: 'blocks.dxf', text: dxf(ins(1000, 1) + ins(2000, -1) + ins(5000, 2)) });
  assert.deepEqual(r.meta.elementTypes, { '3DSOLID': 3, INSERT: 3 });
  assert.equal(r.meta.acis.bodies, 3);
  const [plain, mirrored, stretched] = r.models;
  assert.equal(plain.layer, 'PARTS');
  for (const m of r.models) assert.ok(inspect(m).watertight);
  // single-precision positions a few metres from the origin: the last digits of a 10 mm solid are lost
  close(vol(plain), 875e-9, 1e-11);
  close(vol(mirrored), 875e-9, 1e-11, 'mirrored copy still has its normals outwards');
  close(vol(stretched), 2 * 875e-9, 2e-11);
  const xr = (m) => { let a = Infinity, b = -Infinity; for (let i = 0; i < m.positions.length; i += 3) { a = Math.min(a, m.positions[i]); b = Math.max(b, m.positions[i]); } return [a, b]; };
  assert.deepEqual(xr(plain).map((v) => +v.toFixed(6)), [1, 1.01]);
  assert.deepEqual(xr(mirrored).map((v) => +v.toFixed(6)), [1.99, 2]);
  assert.deepEqual(xr(stretched).map((v) => +v.toFixed(6)), [5, 5.02]);
});

test('binary DXF: solids stored as SAT text and as SAB chunks are read as from the ASCII form', () => {
  const SENTINEL = [...'AutoCAD Binary DXF\r\n\x1a\0'].map((c) => c.charCodeAt(0));
  const toBinary = (dxf) => {
    const lines = dxf.split(/\r?\n/);
    const out = [...SENTINEL];
    const dv = new DataView(new ArrayBuffer(8));
    const put = (n, set) => { set(); for (let k = 0; k < n; k++) out.push(dv.getUint8(k)); };
    for (let i = 0; i + 1 < lines.length; i += 2) {
      const code = parseInt(lines[i], 10), raw = lines[i + 1], v = code === 1 || code === 3 ? raw : raw.trim();
      put(2, () => dv.setInt16(0, code, true));
      const t = groupType(code);
      if (t === 0) { for (const b of enc(v)) out.push(b); out.push(0); } else if (t === 1) put(8, () => dv.setFloat64(0, Number(v), true));
      else if (t === 2) put(2, () => dv.setInt16(0, Number(v), true));
      else if (t === 3) put(4, () => dv.setInt32(0, Number(v), true));
      else if (t === 4) put(8, () => dv.setBigInt64(0, BigInt(v), true));
      else if (t === 5) out.push(Number(v) ? 1 : 0);
      else if (t === 6) { out.push(v.length / 2); for (let k = 0; k < v.length; k += 2) out.push(parseInt(v.slice(k, k + 2), 16)); } else throw new Error(`group ${code}`);
    }
    return Uint8Array.from(out);
  };
  // keep the header, the entities and the solid data; tables, classes and objects use groups this writer need not know
  const keep = (dxf, names) => { const parts = dxf.split(/\n\s*0\r?\nSECTION\r?\n/).slice(1).map((p) => `  0\nSECTION\n${p}`); return parts.filter((p) => names.includes(p.split(/\r?\n/)[3].trim())).join('\n').replace(/\n\s*0\r?\nEOF[\s\S]*$/, '') + '\n  0\nEOF\n'; };
  for (const [name, sections] of [['ezdxf-3dsolids-r2010.dxf', ['HEADER', 'ENTITIES']], ['ezdxf-sphere-r2018.dxf', ['HEADER', 'ENTITIES', 'ACDSDATA']]]) {
    const dxf = keep(text(name), sections);
    const a = parseDxf({ name, text: dxf });
    const b = parseDxfBinary({ name, bytes: toBinary(dxf) });
    assert.equal(b.meta.format, 'dxf-binary');
    assert.deepEqual(b.meta.acis, a.meta.acis, name);
    assert.ok(a.meta.acis.bodies >= 1);
    assert.equal(b.models.length, a.models.length);
    a.models.forEach((m, k) => { assert.deepEqual(b.models[k].positions, m.positions); assert.deepEqual(b.models[k].indices, m.indices); });
  }
});

test('ACIS data that cannot be decoded are counted and named; the rest of the drawing is kept', () => {
  const solid = solidPairs(text('ezdxf-3dsolids-r2010.dxf'));
  const broken = solid.map(([c, v], k) => (c === 1 && k > 12 && k < 16 ? [c, 'zzzz'] : [c, v]));
  const cut = solid.filter(([c]) => c !== 1 && c !== 3).concat([[1, solid.find(([c]) => c === 1)[1]]]);
  const line = '  0\nLINE\n  8\nL\n 10\n0\n 20\n0\n 30\n0\n 11\n100\n 21\n0\n 31\n50\n';
  const dxf = (ents) => `  0\nSECTION\n  2\nENTITIES\n${ents}  0\nENDSEC\n  0\nEOF\n`;
  for (const bad of [broken, cut]) {
    const r = parseDxf({ name: 'mix.dxf', text: dxf(line + asText(bad) + asText(solid)) });
    assert.equal(r.meta.acis.entities, 1);
    assert.equal(r.models.filter((m) => m.kind === 'surface').length, 1);
    assert.equal(r.models.filter((m) => m.kind === 'frame').length, 1);
    const k = Object.keys(r.meta.skipped).filter((q) => /^3DSOLID \(embedded ACIS geometry: /.test(q));
    assert.equal(k.length, 1, JSON.stringify(r.meta.skipped));
    assert.ok(r.warnings.some((w) => /^1 solid \/ surface entity stored as ACIS data could not be read \(1 × 3DSOLID \(embedded ACIS geometry: .+\)\)\. Export them as STEP or IGES, or mesh them/.test(w)), r.warnings.join(' | '));
  }
  assert.throws(() => parseDxf({ name: 'a.dxf', text: dxf(asText(solid)) }, { deflection: -1 }), /deflection: must be a positive length/);
});
