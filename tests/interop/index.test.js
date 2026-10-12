import test from 'node:test';
import assert from 'node:assert/strict';
import { SOLVER_EXPORTS, exportCase, packCase, crc32 } from '../../src/interop/index.js';
import { isZip, zipEntries, readZip, crc32 as zipCrc32 } from '../../src/geometry/zip.js';
import { box } from '../../src/mesh/generate.js';

const HULL = box({ size: [4, 2, 3], center: [0, 0, -0.5], divisions: [4, 2, 3], name: 'hull' });
const FRAME = { kind: 'frame', name: 'beam', nodes: Float64Array.from([0, 0, 0, 3, 0, 0]), members: [{ a: 0, b: 1, section: { type: 'tube', D: 0.3, t: 0.01 } }] };

test('SOLVER_EXPORTS: every entry is complete and writes a case for the model kinds it lists', () => {
  assert.deepEqual(SOLVER_EXPORTS.map((e) => e.id), ['nemoh', 'capytaine', 'wamit-gdf', 'calculix', 'openfoam-wave-tank', 'openfoam-current']);
  const options = {
    nemoh: { omegas: [0.5, 1], headingsDeg: [0, 90] },
    capytaine: { omegas: [0.5, 1] },
    'wamit-gdf': {},
    calculix: { material: { E: 2.1e11, nu: 0.3, rho: 7850 }, supports: [{ nodes: [0] }], loads: [{ type: 'gravity' }] },
    'openfoam-wave-tank': { domain: { length: 60, width: 20, height: 8, waterDepth: 5 }, wave: { H: 1, T: 5 }, cellSize: 1, endTime: 10 },
    'openfoam-current': { domain: { length: 60, width: 30, height: 30 }, current: 1, cellSize: 2 },
  };
  for (const e of SOLVER_EXPORTS) {
    for (const k of ['id', 'solver', 'purpose', 'notes', 'extension']) assert.equal(typeof e[k], 'string', `${e.id}.${k}`);
    assert.ok(Array.isArray(e.kinds) && e.kinds.length > 0 && Array.isArray(e.reads) && e.reads.length > 0);
    assert.equal(typeof e.write, 'function');
    assert.equal(typeof e.read, 'function');
    assert.equal(typeof e.executed, 'boolean');
    const model = e.kinds.includes('surface') ? HULL : FRAME;
    const r = exportCase(e.id, e.id === 'calculix' ? FRAME : model, options[e.id]);
    assert.equal(r.solver, e.solver);
    assert.ok(r.files.length >= 1 && r.notes.length >= 1);
    for (const f of r.files) {
      assert.equal(typeof f.path, 'string');
      assert.ok(typeof f.text === 'string' || f.bytes instanceof Uint8Array, `${e.id}: ${f.path} has text or bytes`);
    }
    // a whole case always packs
    assert.ok(isZip(packCase(r.files)));
  }
  assert.throws(() => exportCase('abaqus', HULL), /unknown id "abaqus"/);
  assert.throws(() => exportCase('nemoh', FRAME, options.nemoh), /does not handle a "frame" model/);
  // the readers behind the registry
  const forces = SOLVER_EXPORTS.find((e) => e.id === 'openfoam-current').read('# Force\n# Time total_x\n1 1 2 3 1 2 3 0 0 0\n');
  assert.deepEqual([forces.Fx[0], forces.Fz[0]], [1, 3]);
  const dat = SOLVER_EXPORTS.find((e) => e.id === 'calculix').read(' displacements (vx,vy,vz) for set NALL and time  0.1000000E+01\n\n         2  1.000000E-03  0.000000E+00  0.000000E+00\n');
  assert.equal(dat.displacement[3], 1e-3);
});

test('packCase: stored ZIP readable by the app\'s ZIP reader, CRC-32 and contents intact', async () => {
  // CRC-32 check value of the ASCII string "123456789" (ISO 3309 / ITU-T V.42)
  assert.equal(crc32(new TextEncoder().encode('123456789')), 0xcbf43926);
  const bytes = Uint8Array.from({ length: 300 }, (_, i) => (i * 37 + 11) & 255);
  const files = [
    { path: 'Nemoh.cal', text: 'rho 1025 — ρ, ω\n' },
    { path: 'mesh/hull.dat', text: '2 0\n' },
    { path: 'results/', text: '' },
    { path: 'constant/triSurface/hull.stl', bytes },
    { path: 'Allrun', text: '#!/bin/sh\n' },
    { path: 'empty.txt', text: '' },
  ];
  const zip = packCase(files, { root: 'case' });
  assert.ok(zip instanceof Uint8Array && isZip(zip));
  const entries = zipEntries(zip);
  assert.deepEqual(entries.map((e) => e.name), ['case/Nemoh.cal', 'case/mesh/hull.dat', 'case/results/', 'case/constant/triSurface/hull.stl', 'case/Allrun', 'case/empty.txt']);
  for (const e of entries) { assert.equal(e.method, 0); assert.equal(e.compressedSize, e.size); assert.equal(e.encrypted, false); }
  assert.equal(entries[2].isDirectory, true);
  assert.equal(entries[3].crc, zipCrc32(bytes));
  assert.equal(entries[3].size, 300);
  const back = await readZip(zip);
  assert.equal(back.size, 5);
  assert.deepEqual(Array.from(back.get('case/constant/triSurface/hull.stl')), Array.from(bytes));
  assert.equal(new TextDecoder().decode(back.get('case/Nemoh.cal')), 'rho 1025 — ρ, ω\n');
  assert.equal(back.get('case/empty.txt').length, 0);
  // layout: local headers + data, central directory, 22-byte end record, nothing else
  const enc = new TextEncoder();
  let expected = 22;
  for (const f of files) { const n = enc.encode(`case/${f.path}`).length; expected += 30 + 46 + 2 * n + (f.bytes ? f.bytes.length : enc.encode(f.text).length); }
  assert.equal(zip.length, expected);
  // scripts are marked executable for unzip (Unix mode in the external attributes)
  const dv = new DataView(zip.buffer);
  let p = dv.getUint32(zip.length - 6, true);
  const modes = {};
  for (let i = 0; i < files.length; i++) { const nl = dv.getUint16(p + 28, true); modes[new TextDecoder().decode(zip.subarray(p + 46, p + 46 + nl))] = dv.getUint32(p + 38, true) >>> 16; p += 46 + nl; }
  assert.equal(modes['case/Allrun'] & 0o777, 0o755);
  assert.equal(modes['case/Nemoh.cal'] & 0o777, 0o644);
  // reproducible, and without a root the paths are kept
  assert.deepEqual(Array.from(packCase(files, { root: 'case' })), Array.from(zip));
  assert.deepEqual(zipEntries(packCase(files)).map((e) => e.name)[1], 'mesh/hull.dat');
  assert.throws(() => packCase([]), /non-empty/);
  assert.throws(() => packCase([{ path: '../x', text: '' }]), /invalid path/);
  assert.throws(() => packCase([{ path: 'a', text: '' }, { path: 'a', text: '' }]), /duplicate/);
});
