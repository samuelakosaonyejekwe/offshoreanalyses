import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseParasolid, readXt } from '../../src/geometry/importers/parasolid.js';
import { importFile, formatFor } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { surfaceArea, signedVolume } from '../../src/mesh/measure.js';
import { close, assertContract } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const text = (name) => new TextDecoder().decode(real(name));
const vol = (m) => signedVolume(m.positions, m.indices);
const area = (m) => surfaceArea(m.positions, m.indices);

// The worked example of the Parasolid XT Format Reference: a sheet body with one circular face of radius 1 and
// a colour attribute on its edge, schema 12006 (no schema description in the file).
const HEAD = (sch, usfld = 0) => `**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**************************
**PARASOLID !"#$%&'()*+,-./:;<=>?@[\\]^_\`{|}~0123456789**************************
**PART1;MC=osf65;MC_MODEL=alpha;OS=OSF1;APPL=unknown;FORMAT=text;GUISE=transmit;DATE=29-mar-2000;
**PART2;SCH=${sch};USFLD_SIZE=${usfld};
**PART3;
**END_OF_HEADER*****************************************************************
`;
const NODES = `12 1 12 0 2 0 0 0 0 1e3 1e-8 0 0 0 1 0 3 1 3 4 5 0 6 7 0 70 2 0 1 0 0 4 1 20 8 8 8 1 T13 3 3 0 1 0 9 0 0 6 9 50 4 11 0 9 0 0 0 +0 0 0 0 0 1 1 0 0 31 5 10 0 7 0 0 0 +0 0 0 0 0 1 1 0 0 1 19 6 5 0 1 0 0 3 V16 7 6 0 ?10 0 0 5 0 0 1 17 10 0 11 10 10 0 12 7 0 0 +15 11 7 0 10 9 0 17 12 0 0 0 0 0 10 7 0 0 -14 9 2 13 ?0 0 11 3 4 +0 0 0 0 3 81 1 13 12 14 9 0 0 0 0 15 80 1 14 0 16 8001 0 0 0 0 3 5 0 0 FFFFTFTFFFFFF2 83 3 15 1 2 3 79 15 16 SDL/TYSA_COLOUR74 20 8 1 0 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 `;
// the stream is cut into lines at arbitrary places, as files are
const wrap = (s) => s.replace(/(.{1,80})/g, '$1\n');
const example = (nodes) => HEAD('SCH_1200000_12006') + wrap(`T51 : TRANSMIT FILE created by modeller version 120000017 SCH_1200000_120060 ${nodes}`);
const EXAMPLE = example(NODES);

test('Parasolid text: header, node stream and schema of the reference example (schema 12006)', () => {
  const xt = readXt(EXAMPLE, 'circle.x_t');
  assert.deepEqual(xt.header, { modeller: 1200000, schema: 'SCH_1200000_12006', schemaNumber: 12006, baseSchema: null, application: 'unknown', date: '29-mar-2000', format: 'text', partition: false });
  assert.equal(xt.userField, 0);
  // an attribute definition of that release has 13 owner classes and no field names: the layout is found by trying
  assert.equal(xt.layout, 'old-attdef');
  assert.deepEqual(xt.types, { BODY: 1, LIST: 1, SHELL: 1, PLANE: 1, CIRCLE: 1, REGION: 1, EDGE: 1, FIN: 2, LOOP: 1, FACE: 1, ATTRIBUTE: 1, ATTRIB_DEF: 1, REAL_VALUES: 1, ATT_DEF_ID: 1, POINTER_LIS_BLOCK: 1 });
  assert.equal(xt.nodes.size, 16);
  const n = (i) => xt.nodes.get(i);
  assert.equal(n(1).name, 'BODY');
  assert.deepEqual([n(1).res_size, n(1).res_linear, n(1).body_type, n(1).shell, n(1).boundary_surface, n(1).boundary_curve, n(1).region, n(1).edge], [1e3, 1e-8, 3, 3, 4, 5, 6, 7]);
  assert.deepEqual([n(5).name, n(5).sense, n(5).centre, n(5).normal, n(5).x_axis, n(5).radius], ['CIRCLE', '+', [0, 0, 0], [0, 0, 1], [1, 0, 0], 1]);
  // unset values ("?") and characters that run into the next number
  assert.deepEqual([n(7).tolerance, n(7).fin, n(7).curve], [null, 10, 5]);
  assert.deepEqual([n(9).name, n(9).tolerance, n(9).loop, n(9).surface, n(9).sense], ['FACE', null, 11, 4, '+']);
  assert.deepEqual([n(6).name, n(6).kind], ['REGION', 'V']);
  assert.deepEqual([n(10).sense, n(12).sense, n(10).other, n(12).loop], ['+', '-', 12, 0]);
  // variable-length nodes and character arrays
  assert.equal(n(16).string, 'SDL/TYSA_COLOUR');
  assert.deepEqual(n(15).values, [1, 2, 3]);
  assert.equal(n(8).entries.length, 20);
  assert.deepEqual([n(14).type_id, n(14).legal_owners.length, n(14).fields], [8001, 13, [2]]);

  const r = parseParasolid({ name: 'circle.x_t', text: EXAMPLE });
  assertContract(r, ['surface']);
  assert.deepEqual(r.warnings, []);
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.deepEqual([r.meta.bodies, r.meta.instances, r.meta.closedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 0, 1, 1]);
  assert.deepEqual(r.meta.elementTypes, { plane: 1 });
  // a disc of radius 1 in the plane z = 0: 24 chords at the default 15°
  close(area(r.models[0]), 12 * Math.sin(Math.PI / 12), 1e-6);
  close(area(parseParasolid({ name: 'circle.x_t', text: EXAMPLE }, { angularStep: 1 }).models[0]), Math.PI, 2e-4);
  for (let i = 2; i < r.models[0].positions.length; i += 3) assert.equal(r.models[0].positions[i], 0);
});

test('Parasolid text: schema edits in front of the first node of a type (files of Parasolid 14 and later)', () => {
  // the same body as a newer release writes it: the base schema is named, and BODY, LIST, REGION and the pointer
  // block carry their differences from it; every other type is flagged "as the base" (255). User field of one word.
  const body = '12 27 CCCCCCCCCCCCCDI5 owner1040 0 CCCCCCCCCA16 index_map_offset0 0 1 dA9 index_map82 0 A17 node_id_index_map82 0 A20 schema_embedding_map82 0 Z1 12 0 2 0 0 0 0 1e3 1e-8 0 0 0 1 0 3 1 3 4 5 0 6 7 0 0 0 0 0 9 ';
  const list = '70 11 CI9 list_type0 0 1 uI10 notransmit0 0 1 lCCCDCCDI12 finger_index0 0 1 dI12 finger_block1012 0 CZ2 0 4 T1 0 0 1 20 1 8 8 9 ';
  const rest = '13 255 3 3 0 1 0 9 0 0 6 9 9 50 255 4 11 0 9 0 0 0 +0 0 0 0 0 1 1 0 0 9 31 255 5 10 0 7 0 0 0 +0 0 0 0 0 1 1 0 0 1 9 19 8 CCCCCCCA5 owner12 0 Z6 5 0 1 0 0 3 V0 9 16 255 7 6 0 ?10 0 0 5 0 0 1 9 17 255 10 0 11 10 10 0 12 7 0 0 +9 15 255 11 7 0 10 9 0 9 17 12 0 0 0 0 0 10 7 0 0 -9 14 255 9 2 0 ?0 0 11 3 4 +0 0 0 0 3 9 74 4 CI16 index_map_offset0 0 1 dCCZ20 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 1 0 ';
  const mk = (stream) => HEAD('SCH_2700142_26105', 1) + wrap(stream);
  const stream = `T51 : TRANSMIT FILE created by modeller version 270014223 SCH_2700142_26105_13006196 1 ${body}${list}${rest}`;
  const file = mk(stream);
  const xt = readXt(file, 'new.x_t');
  assert.deepEqual([xt.header.schemaNumber, xt.header.baseSchema, xt.userField, xt.layout], [26105, 13006, 1, 'described by the file']);
  assert.deepEqual(xt.schema.get(12).map((f) => f.name).slice(-4), ['index_map_offset', 'index_map', 'node_id_index_map', 'schema_embedding_map']);
  assert.equal(xt.schema.get(12).length, 27);
  assert.deepEqual(xt.schema.get(70).map((f) => `${f.name}:${f.type}`), ['node_id:d', 'list_type:u', 'notransmit:l', 'owner:p', 'next:p', 'previous:p', 'list_length:d', 'block_length:d', 'finger_index:d', 'finger_block:p', 'list_block:p']);
  assert.equal(xt.nodes.get(2).notransmit, true);
  assert.equal(xt.nodes.get(6).owner, 0);
  const r = parseParasolid({ name: 'new.x_t', text: file });
  close(area(r.models[0]), 12 * Math.sin(Math.PI / 12), 1e-6);
  // an edit that does not add up, a field type or a node type nobody described: the file is refused, not guessed at
  assert.throws(() => readXt(mk(stream.replace('12 27 CCCC', '12 26 CCCC')), 'bad.x_t'), /bad\.x_t: the node stream could not be interpreted \(schema SCH_2700142_26105_13006\): after 0 nodes, the schema edit of node type 12 gives 27 fields where 26 are announced\. Nothing was read from the file; export STEP/);
  assert.throws(() => readXt(mk(stream.replace('13 255 3 3 0', '299 255 3 3 0')), 'bad.x_t'), /after 2 nodes, node type 299 is declared equal to the base schema, which has no such type/);
  assert.throws(() => readXt(mk(stream.replace('_13006196 1 ', '_12006196 1 ')), 'bad.x_t'), /written against the base schema 12006; only 13006 is known to this reader/);
});

test('real Parasolid, Onshape (schema 32001, embedded): a 100 mm block with a spherical cavity of radius 10 mm', async () => {
  const exact = 1e-3 - (4 / 3) * Math.PI * 0.01 ** 3;
  let prev = Infinity;
  for (const angularStep of [15, 5, 2]) {
    const r = await importFile({ name: 'cavity.x_t', bytes: real('onshape-sphere-in-box.x_t') }, { angularStep });
    assert.equal(r.format, 'parasolid');
    assertContract(r, ['surface']);
    assert.deepEqual(r.warnings, []);
    assert.equal(r.meta.schema, 'SCH_3300229_32001_13006');
    assert.equal(r.meta.schemaLayout, 'described by the file');
    assert.equal(r.meta.application, 'Onshape');
    assert.deepEqual(r.meta.elementTypes, { plane: 6, sphere: 2 });
    assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 1, 1, 8, 8]);
    const m = r.models[0];
    assert.ok(inspect(m).watertight);
    // the cavity is faceted from inside the sphere, so the solid comes out slightly too large
    const err = vol(m) - exact;
    assert.ok(err > 0 && err < prev, `${angularStep}°: ${err}`);
    prev = err;
    close(r.meta.volumes[0], vol(m), 1e-9);
  }
  assert.ok(prev / exact < 2e-6);
});

test('real Parasolid, SolidWorks (schema 10004, not described by the file): four bars with exact volumes', async () => {
  const r = await importFile({ name: 'bars.x_t', bytes: real('solidworks-bars-v10.x_t') });
  assert.equal(r.meta.schema, 'SCH_1000230_10004');
  // read with the layout of that release, found by trying the known ones against the whole stream
  assert.equal(r.meta.schemaLayout, 'old-attdef, bare-fin');
  assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [4, 4, 4, 4, 24, 24]);
  assert.deepEqual(r.models.map((m) => m.name), ['bars body 1', 'bars body 2', 'bars body 3', 'bars body 4']);
  [12, 40, 20, 8].forEach((v, k) => { close(vol(r.models[k]), v, 1e-5 * v); assert.ok(inspect(r.models[k]).watertight); });
  assert.deepEqual(r.meta.bbox, { min: [-4.440892098500626e-16, -2.220446049250313e-16, 0], max: [2, 2, 20] });
});

test('real Parasolid, SolidWorks (schemas 19008 and 20000): machined parts with holes and countersinks', () => {
  for (const [name, hist, faces, size] of [['solidworks-load-cell.x_t', { cylinder: 16, plane: 27 }, 43, [0.0889, 0.0635, 0.03048]], ['solidworks-plate.x_t', { cone: 4, plane: 14, cylinder: 8 }, 26, [0.0889, 0.0889, 0.0254]]]) {
    const r = parseParasolid({ name, bytes: real(name) }, { angularStep: 5 });
    assert.deepEqual(r.warnings, [], name);
    assert.match(r.meta.application, /^SolidWorks/);
    assert.deepEqual(r.meta.elementTypes, hist);
    assert.deepEqual([r.meta.bodies, r.meta.solidBodies, r.meta.closedBodies, r.meta.orientedBodies, r.meta.faces, r.meta.faceTotal], [1, 1, 1, 1, faces, faces]);
    const m = r.models[0];
    assert.ok(inspect(m).watertight, name);
    assert.equal(m.groups.length, faces);
    [0, 1, 2].forEach((k) => close(r.meta.bbox.max[k] - r.meta.bbox.min[k], size[k], 1e-7, `${name} extent ${k}`));
    // a part with holes: less than its bounding block, more than half of it
    const block = size[0] * size[1] * size[2];
    assert.ok(vol(m) > 0.5 * block && vol(m) < block, `${name}: ${vol(m)} of ${block}`);
  }
});

test('Parasolid: truncated and damaged files are refused with the reason; nothing is guessed', async () => {
  const good = text('solidworks-plate.x_t');
  await assert.rejects(importFile({ name: 'part.x_b', bytes: Uint8Array.from([0x50, 0x53, 0, 0, 0, 51, 58, 32, 84, 82]) }), /part\.x_b: the start of the data could not be read: the file ends inside the modeller version\./);
  assert.throws(() => parseParasolid({ name: 'n.x_t', text: 'solid cube\n' }), /n\.x_t: not a Parasolid text transmit file \(the "\*\*…" header ending in "\*\*END_OF_HEADER" is missing\)/);
  const cut = good.slice(0, Math.floor(good.length * 0.6));
  assert.throws(() => parseParasolid({ name: 't.x_t', text: cut }), /t\.x_t: (the file is truncated: the node stream ends after \d+ nodes without its terminator|the node stream could not be interpreted .*the file ends inside)/);
  // a letter where a number belongs
  const at = good.indexOf(' 50 ', good.indexOf('END_OF_HEADER') + 400);
  assert.throws(() => parseParasolid({ name: 'd.x_t', text: `${good.slice(0, at + 9)}x${good.slice(at + 10)}` }), /d\.x_t: the node stream could not be interpreted \(schema SCH_1901315_19008_13006\): after \d+ nodes, .* Nothing was read from the file; export STEP from the originating program\./);
  // a file of a release whose layout is not known in full is refused as a whole
  const odd = example(NODES.replace('17 10 0 11 10 10 0 12 7 0 0 +', '17 10 0 11 10 10 0 12 7 0 0 0 +'));
  assert.throws(() => parseParasolid({ name: 'o.x_t', text: odd }), /o\.x_t: the node stream could not be interpreted \(schema SCH_1200000_12006, which the file does not describe; \d known layouts were tried\)/);
  // a face whose surface points to something that is not a surface stops its body
  const wrong = example(NODES.replace('14 9 2 13 ?0 0 11 3 4 +', '14 9 2 13 ?0 0 11 3 5 +'));
  assert.throws(() => parseParasolid({ name: 'w.x_t', text: wrong }), /w\.x_t: no geometry could be read from the Parasolid file\. w was not read: FACE\.surface points to a CIRCLE node\. Nothing of this body is shown/);
  assert.equal(formatFor('part.x_t').status, 'native');
  assert.equal(formatFor('part.x_t').via, 'importers/parasolid.js');
  assert.throws(() => parseParasolid({ name: 'a.x_t', text: EXAMPLE }, { angularStep: 0 }), /angularStep: must be within 0…90 degrees/);
});
