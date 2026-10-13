import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseHst } from '../../src/geometry/importers/panels.js';
import { exportHst, exportModels } from '../../src/geometry/exporters.js';
import { importFile, formatFor } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { enc, close, assertContract, assertTruncationSafe, surfaceArea, surfaceVolume } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));

// Wetted quarter of a box 2 × 2 × 1 below z = 0 (x ≥ 0, y ≥ 0): bottom and two sides, normals outwards.
// Node numbers start at 0 and skip values, as in files written by HydroStar's mesher.
const NODES = [[0, 0, 0, -1], [1, 1, 0, -1], [2, 1, 1, -1], [3, 0, 1, -1], [10, 1, 0, 0], [11, 1, 1, 0], [12, 0, 1, 0]];
const PANELS = [[0, 3, 2, 1], [1, 2, 11, 10], [2, 3, 12, 11]];
const hst = ({ sym = 'SYMMETRY_BODY 1 2', type = 1, coords = 'COORDINATES', panel = 'PANEL', end = 'ENDPANEL', tail = 'ENDFILE\n' } = {}) => `PROJECT quarter box
USER test

NBBODY 1
${sym}
RHO 1025.0
GRAVITY 9.81
NUMPANEL 1 1 3
${coords}
${NODES.map((n) => n.join(' ')).join('\n')}
ENDCOORDINATES
${panel} TYPE ${type}
${PANELS.map((p, i) => (type ? [i + 1, ...p] : p).join(' ')).join('\n')}
${end}
${tail}`;

test('HydroStar: quarter body with SYMMETRY_BODY 1 2 is mirrored about y = 0 and x = 0', () => {
  const r = parseHst({ name: 'box.hst', bytes: enc(hst()) });
  assertContract(r, ['surface']);
  const m = r.models[0];
  assert.equal(r.meta.format, 'hydrostar');
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.deepEqual([r.meta.symmetry.x, r.meta.symmetry.y, r.meta.symmetry.expanded], [true, true, true]);
  assert.deepEqual([r.meta.panelsStored, r.meta.nodes, r.meta.project, r.meta.user, r.meta.declaredBodies], [3, 7, 'quarter box', 'test', 1]);
  assert.deepEqual(r.meta.elementTypes, { quad4: 12 });
  assert.deepEqual(m.groups.map((g) => g.name), ['as stored', 'mirror x', 'mirror y', 'mirror xy']);
  assert.equal(m.quads.length, 48);
  assert.deepEqual([r.meta.bbox.min, r.meta.bbox.max], [[-1, -1, -1], [1, 1, 0]]);
  // bottom 2 × 2 and four sides 2 × 1; the open top at z = 0 closes the volume integral
  close(surfaceArea(m), 12, 1e-6);
  close(surfaceVolume(m), 4, 1e-6);
  assert.equal(inspect(m).stats.boundaryEdges, 8);
  // RHO and GRAVITY are not geometry: named, not silently dropped
  assert.deepEqual(r.meta.keywordsSkipped, { RHO: 1, GRAVITY: 1 });
  assert.match(r.warnings.join(' '), /not read: RHO, GRAVITY/);

  const stored = parseHst({ name: 'box.hst', text: hst() }, { expandSymmetry: false });
  close(surfaceArea(stored.models[0]), 3, 1e-6);
  assert.equal(stored.meta.symmetry.expanded, false);
});

test('HydroStar: SYMMETRY 1 (half body), panel types 0 and 1, keyword spellings of other writers', () => {
  const half = parseHst({ name: 'h.hst', text: hst({ sym: 'SYMMETRY 1' }) });
  close(surfaceArea(half.models[0]), 6, 1e-6);
  assert.deepEqual(half.models[0].groups.map((g) => g.name), ['as stored', 'mirror y']);
  assert.deepEqual([half.meta.bbox.min, half.meta.bbox.max], [[0, -1, -1], [1, 1, 0]]);
  const none = parseHst({ name: 'h.hst', text: hst({ sym: 'SYMMETRY_BODY 1 0', type: 0 }) });
  close(surfaceArea(none.models[0]), 3, 1e-6);
  assert.equal(none.models[0].groups, undefined);
  const other = parseHst({ name: 'h.hst', text: hst({ coords: 'COORDINATES TYPE 0', panel: 'PANELS', end: 'ENDPANELS', type: 0 }) });
  close(surfaceArea(other.models[0]), 12, 1e-6);
  // triangles (three node numbers) and a node table without numbers (numbered from 1 in order)
  const tri = parseHst({ name: 't.hst', text: 'COORDINATES\n0 0 -1\n1 0 -1\n1 1 -1\n0 1 -1\nENDCOORDINATES\nPANEL TYPE 0\n1 3 2\n1 4 3\nENDPANEL\nENDFILE\n' });
  assert.deepEqual(tri.meta.elementTypes, { tri3: 2 });
  close(surfaceArea(tri.models[0]), 1, 1e-6);
  const bad = parseHst({ name: 'h.hst', text: hst({ sym: 'SYMMETRY_BODY 1 7' }) });
  assert.match(bad.warnings.join(' '), /SYMMETRY_BODY = 7 is not a known symmetry code/);
  close(surfaceArea(bad.models[0]), 3, 1e-6);
});

test('HydroStar: bodies of a multi-body file are separated through NUMPANEL, each with its own symmetry', () => {
  const text = `NBBODY 2\nSYMMETRY_BODY 1 1\nSYMMETRY_BODY 2 0\nNUMPANEL 1 1 1\nNUMPANEL 2 2 2\nCOORDINATES\n1 0 0 -1\n2 1 0 -1\n3 1 1 -1\n4 0 1 -1\n5 5 2 -1\n6 6 2 -1\n7 6 3 -1\n8 5 3 -1\nENDCOORDINATES\nPANEL TYPE 1\n1 1 4 3 2\n2 5 8 7 6\nENDPANEL\nENDFILE\n`;
  const r = parseHst({ name: 'two.hst', text });
  assert.deepEqual(r.models[0].groups.map((g) => g.name), ['body 1', 'body 1 (mirror y)', 'body 2']);
  assert.deepEqual(r.meta.bodies, [{ body: 1, panels: 1, symmetry: 1 }, { body: 2, panels: 1, symmetry: 0 }]);
  close(surfaceArea(r.models[0]), 3, 1e-6);
  assert.throws(() => parseHst({ name: 'two.hst', text: text.replace('NUMPANEL 2 2 2', 'NUMPANEL 2 2 9') }), /two\.hst: the bodies carry different symmetry codes/);
});

test('HydroStar: file written by meshmagick, and the exporter read back', async () => {
  const bytes = real('meshmagick-box.hst');
  assert.equal(formatFor('box.hst', bytes).id, 'hydrostar');
  assert.equal(formatFor('hull.txt', bytes).id, 'hydrostar');
  const r = await importFile({ name: 'box.hst', bytes });
  assert.equal(r.format, 'hydrostar');
  const m = r.models[0];
  // a closed box 4 × 2 × 1 of 108 panels on 110 nodes
  assert.deepEqual([r.meta.panelsStored, r.meta.nodes, m.positions.length / 3], [108, 110, 110]);
  close(surfaceArea(m), 28, 1e-5);
  close(surfaceVolume(m), 8, 1e-5);
  assert.equal(inspect(m).watertight, true);

  const out = exportModels('hst', m);
  assert.equal(out.fileName, 'box.hst');
  assert.match(out.data, /^PROJECT box\nNBBODY 1\nSYMMETRY_BODY 1 0\nNUMPANEL 1 1 108\nCOORDINATES\n1 -2\.0+e\+0 /);
  const back = parseHst({ name: 'b.hst', text: out.data });
  assert.deepEqual(back.meta.elementTypes, r.meta.elementTypes);
  assert.deepEqual(Array.from(back.models[0].quads), Array.from(m.quads));
  assert.deepEqual(Array.from(back.models[0].positions), Array.from(m.positions));
  // a surface without panels is written triangle by triangle
  const t = parseHst({ name: 't.hst', text: exportHst({ kind: 'surface', name: 'tri', positions: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0]), indices: Uint32Array.from([0, 1, 2]), origin: [10, 20, -5] }) });
  assert.deepEqual(t.meta.elementTypes, { tri3: 1 });
  assert.deepEqual([t.meta.bbox.min, t.meta.bbox.max], [[10, 20, -5], [11, 21, -5]]);
  assert.throws(() => exportHst({ kind: 'points', positions: new Float32Array(3) }), /HydroStar export does not handle a "points" model/);
});

test('HydroStar: malformed and truncated files are refused or reported, never guessed', async () => {
  const P = (text) => parseHst({ name: 'x.hst', text });
  assert.throws(() => P('hello world\n1 2 3\n'), /x\.hst: not a HydroStar mesh \(expected keyword lines such as PROJECT/);
  assert.throws(() => P('PROJECT a\nNBBODY 1\n'), /no COORDINATES … ENDCOORDINATES node table was found/);
  assert.throws(() => P('COORDINATES\n1 0 0 0\n2 1 0 0\n'), /the file ends inside the node table \(no ENDCOORDINATES line\)/);
  assert.throws(() => P('COORDINATES\n1 0 0 0\nENDCOORDINATES\nENDFILE\n'), /holds 1 nodes but no panel table/);
  assert.throws(() => P('COORDINATES\n1 0 0 0\nENDCOORDINATES\nPANEL TYPE 2\n'), /only PANEL TYPE 0 .* and PANEL TYPE 1/);
  assert.throws(() => P('PANEL TYPE 0\n1 2 3 4\nENDPANEL\n'), /the panel table comes before the node table/);
  assert.throws(() => P(hst().replace('ENDPANEL\nENDFILE\n', '')), /the file ends inside the panel table/);
  assert.throws(() => P('COORDINATES\n1 0 0 0\n2 1 0 0\n3 1 1 0\nENDCOORDINATES\nPANEL TYPE 0\n1 2 9\nENDPANEL\n'), /no panel of the panel table could be read \(0 malformed row\(s\), 1 naming nodes that are not in the node table\)/);
  const r = P(hst({ tail: '' }).replace('1 0 3 2 1', '1 0 3 2 1\n9 0 3 77 1\nx y z').replace('3 0 1 -1', '3 0 1 -1\n4 a b c'));
  assert.match(r.warnings.join(' '), /no ENDFILE line/);
  assert.match(r.warnings.join(' '), /1 row\(s\) of the node table without three numeric coordinates were skipped/);
  assert.match(r.warnings.join(' '), /1 row\(s\) of the panel table that are not three or four node numbers were skipped/);
  assert.match(r.warnings.join(' '), /1 panel\(s\) reference node numbers that are not in the node table/);
  await assertTruncationSafe(enc(hst()), (b) => parseHst({ name: 'b.hst', bytes: b }), 'HST');
  await assertTruncationSafe(real('meshmagick-box.hst'), (b) => parseHst({ name: 'b.hst', bytes: b }), 'HST meshmagick');
});
