import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseAcis, readSat } from '../../src/geometry/importers/acis.js';
import { importFile } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { close, assertContract } from './fixtures/helpers.js';

// SAT files written by CAD systems (tests/geometry/fixtures/real/README.md gives source and licence of each).
// What is checked is what can be known without the originating program: sizes that follow from the stored
// dimensions, closure of the surface, and that nothing is passed off as read when it was not.

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const vol = (m) => signedVolume(m.positions, m.indices);
const area = (m) => surfaceArea(m.positions, m.indices);

test('real SAT, Fusion 360 (ASM 222, layout 7.0): a 100 × 50 × 100 mm box', async () => {
  const r = await importFile({ name: 'fusion-box.sat', bytes: real('fusion-box.sat') });
  assert.equal(r.format, 'acis');
  assertContract(r, ['surface']);
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.acisVersion, 700);
  assert.equal(r.meta.product, 'Autodesk Translation Framework');
  assert.match(r.meta.acisRelease, /^ASM 222/);
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 1, 1, 6, 6]);
  assert.deepEqual(r.meta.elementTypes, { plane: 6 });
  const m = r.models[0];
  assert.ok(inspect(m).watertight);
  // positions are single precision
  close(vol(m), 0.1 * 0.05 * 0.1, 1e-10);
  close(area(m), 2 * (0.1 * 0.05 + 0.1 * 0.1 + 0.05 * 0.1), 1e-8);
  close(r.meta.volumes[0], 5e-4, 1e-12);
  // the attribute records Fusion writes are passed over, not reported as unknown geometry
  const sat = readSat(new TextDecoder().decode(real('fusion-box.sat')));
  assert.ok(sat.records.some((q) => q.kind === 'attrib'));
  assert.deepEqual(Object.keys(r.meta.skipped), []);
});

test('real SAT, ACIS 28 (four-line header, bounding boxes in the records): a rod Ø 3.175 × 100 mm', () => {
  const bytes = real('acis28-rod.sat');
  const sat = readSat(new TextDecoder().decode(bytes));
  assert.equal(sat.header.version, 2800);
  assert.equal(sat.header.acis, 'ACIS 28.0 NT');
  // the stream identifier of the fourth header line is not taken for a record
  assert.equal(sat.records[0].type, 'body');
  const R = 1.5875e-3, L = 0.1;
  let prev = Infinity;
  for (const angularStep of [15, 5, 2]) {
    const r = parseAcis({ name: 'rod.sat', bytes }, { angularStep });
    assert.deepEqual(r.warnings, []);
    assert.deepEqual(r.meta.elementTypes, { cylinder: 1, plane: 2 });
    assert.deepEqual([r.meta.bodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 1, 3, 3]);
    const m = r.models[0];
    assert.ok(inspect(m).watertight);
    const n = 360 / angularStep;
    // an inscribed n-sided prism
    close(vol(m), (n / 2) * Math.sin((2 * Math.PI) / n) * R * R * L, 1e-12, `prism of ${n} sides`);
    const err = Math.PI * R * R * L - vol(m);
    assert.ok(err > 0 && err < prev);
    prev = err;
  }
  assert.ok(prev / (Math.PI * R * R * L) < 3e-4);
});

test('real SAT, Revit 2014: bodies with geometry the file itself marks as illegal are read as far as they are sound', () => {
  const bytes = real('revit-solids.sat');
  const sat = readSat(new TextDecoder().decode(bytes));
  // four edges carry the attribute Revit attaches to geometry that failed its own check
  assert.equal(sat.records.filter((q) => q.type === 'illegal_geometry-sys-attrib').length, 4);
  // tolerant topology and sub-type references are present in the file
  assert.ok(sat.records.some((q) => q.type === 'tedge-edge') && sat.records.some((q) => q.type === 'tcoedge-coedge') && sat.records.some((q) => q.type === 'tvertex-vertex'));
  assert.ok(sat.records.some((q) => q.type === 'intcurve-curve' && q.tok.some((t, i) => t.t === '{' && q.tok[i + 1].v === 'ref')));
  const r = parseAcis({ name: 'revit.sat', bytes });
  assertContract(r, ['surface']);
  assert.equal(r.meta.product, 'Autodesk Revit 2014');
  assert.equal(r.meta.bodies, 3);
  assert.equal(r.models.length, 3);
  assert.deepEqual([r.meta.faces, r.meta.faceTotal, r.meta.facesSkipped], [10, 15, 5]);
  // the curves of the marked edges lie 200 mm from their vertices: the faces that depend on them are refused, with the distance
  assert.ok(r.warnings.some((w) => /^5 of 15 faces could not be read \(.*its boundary lies 2\.0e\+2 model units off.*\) — export those bodies as STEP/.test(w)), r.warnings.join(' | '));
  assert.ok(r.warnings.some((w) => /An intcurve-curve whose B-spline does not pass through its edge vertices .* \(×4\)/.test(w)));
  assert.equal(r.meta.closedBodies, 0);
  // every face that was read is planar and lies where its vertices are
  assert.deepEqual(r.meta.elementTypes, { plane: 10 });
  for (const m of r.models) for (const v of m.positions) assert.ok(Number.isFinite(v));
});
