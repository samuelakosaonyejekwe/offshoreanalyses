import test from 'node:test';
import assert from 'node:assert/strict';
import { parseStl, parseObj, parsePly, parseOff } from '../../src/geometry/importers/surfaces.js';
import { parseGltf, parse3mf } from '../../src/geometry/importers/gltf.js';
import { makeZip } from '../../src/geometry/zip.js';
import { enc, close, CUBE_V, CUBE_Q, CUBE_T, surfaceVolume, surfaceArea, binaryStl, concat, assertContract } from './fixtures/helpers.js';

const asciiStl = (name, tris, verts) => `solid ${name}\n` + tris.map((t) => ` facet normal 0 0 0\n  outer loop\n${t.map((k) => `   vertex ${verts[k].join(' ')}`).join('\n')}\n  endloop\n endfacet`).join('\n') + `\nendsolid ${name}\n`;

test('STL binary: unit cube welds to 8 vertices, 12 triangles, volume 1', () => {
  const r = parseStl({ name: 'cube.stl', bytes: binaryStl(CUBE_T, CUBE_V) });
  assertContract(r, ['surface']);
  const m = r.models[0];
  assert.equal(m.positions.length, 24);
  assert.equal(m.indices.length, 36);
  close(surfaceVolume(m), 1, 1e-6);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [1, 1, 1] });
  assert.equal(r.meta.encoding, 'binary');
  assert.equal(r.meta.units.basis, 'assumed');
  assert.deepEqual(r.meta.elementTypes, { tri3: 12 });
});

test('STL ascii: CRLF, two solids become two groups; units option scales', () => {
  const half = CUBE_T.slice(0, 6), rest = CUBE_T.slice(6);
  const text = (asciiStl('lower', half, CUBE_V) + asciiStl('upper', rest, CUBE_V)).replace(/\n/g, '\r\n');
  const r = parseStl({ name: 'two.stl', text }, { units: 'mm' });
  assertContract(r, ['surface']);
  const m = r.models[0];
  assert.equal(m.groups.length, 2);
  assert.deepEqual(m.groups.map((g) => g.name), ['lower', 'upper']);
  assert.equal(m.groups[0].count + m.groups[1].count, 36);
  close(surfaceVolume(m), 1e-9, 1e-15);
  assert.equal(r.meta.units.unit, 'mm');
  assert.equal(r.meta.encoding, 'ascii');
});

test('STL binary whose header starts with "solid" is still read as binary', () => {
  const r = parseStl({ name: 'x.stl', bytes: binaryStl(CUBE_T, CUBE_V, 'solid exported by a CAD system') });
  assert.equal(r.meta.encoding, 'binary');
  assert.equal(r.meta.counts.triangles, 12);
});

test('STL: implausibly large unit-less model gets a millimetre hint; degenerate facets are dropped', () => {
  const big = CUBE_V.map((v) => v.map((c) => c * 20000));
  const r = parseStl({ name: 'big.stl', bytes: binaryStl([...CUBE_T, [0, 0, 1]], big) });
  assert.equal(r.meta.counts.triangles, 12);
  assert.equal(r.meta.skipped.degenerateFacets, 1);
  assert.ok(r.warnings.some((w) => /millimetres/.test(w)));
});

test('STL: malformed input errors say what was expected', () => {
  assert.throws(() => parseStl({ name: 'a.stl', bytes: enc('hello') }), /too short/);
  assert.throws(() => parseStl({ name: 'a.stl', text: 'solid x\n facet normal 0 0 0\n outer loop\n vertex 1 2\n' }), /expected "vertex x y z"/);
  assert.throws(() => parseStl({ name: 'a.stl', text: 'solid empty\nendsolid empty\n' + ' '.repeat(80) }), /no "vertex x y z" records/);
});

test('OBJ: quads, negative indices, groups, polyline and y-up option', () => {
  const text = [
    '# units: mm',
    ...CUBE_V.map((v) => `v ${v.map((c) => c * 1000).join(' ')}`),
    'vn 0 0 1',
    'g bottom', `f ${CUBE_Q[0].map((i) => `${i + 1}//1`).join(' ')}`,
    'g rest', ...CUBE_Q.slice(1).map((q) => `f ${q.map((i) => i - 8).join(' ')}`),
    'v 0 0 2000', 'v 0 0 3000', 'l -2 -1',
    'curv 0 1 1 2',
  ].join('\n');
  const r = parseObj({ name: 'cube.obj', text });
  assertContract(r, ['surface', 'polyline']);
  const s = r.models.find((m) => m.kind === 'surface');
  assert.equal(s.positions.length, 24, 'curve vertices are not surface vertices');
  close(surfaceVolume(s), 1, 1e-6);
  assert.deepEqual(s.groups.map((g) => [g.name, g.count]), [['bottom', 6], ['rest', 30]]);
  const l = r.models.find((m) => m.kind === 'polyline');
  assert.deepEqual([...l.points], [0, 0, 2, 0, 0, 3]);
  assert.equal(r.meta.units.unit, 'mm');
  assert.equal(r.meta.units.basis, 'file');
  assert.deepEqual(r.meta.elementTypes, { quad4: 6, polyline: 1 });
  assert.equal(r.meta.skipped.curv, 1);
  const y = parseObj({ name: 'cube.obj', text }, { upAxis: 'y', units: 'm' });
  assert.deepEqual(y.meta.bbox.min.slice(0, 3), [0, -3000, 0]);
});

test('OBJ: concave polygon is triangulated without changing its area', () => {
  // L-shaped hexagon of area 3
  const text = 'v 0 0 0\nv 2 0 0\nv 2 1 0\nv 1 1 0\nv 1 2 0\nv 0 2 0\nf 1 2 3 4 5 6\n';
  const r = parseObj({ name: 'l.obj', text });
  close(surfaceArea(r.models[0]), 3, 1e-9);
  assert.equal(r.models[0].indices.length, 12);
  assert.throws(() => parseObj({ name: 'e.obj', text: 'f 1 2 3\n' }), /no vertex records/);
});

function plyCube(format) {
  const header = `ply\nformat ${format} 1.0\ncomment units: m\nelement vertex 8\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face 6\nproperty list uchar int vertex_indices\nelement edge 1\nproperty int vertex1\nproperty int vertex2\nend_header\n`;
  if (format === 'ascii') return enc(header + CUBE_V.map((v, i) => `${v.join(' ')} ${i} 0 255`).join('\n') + '\n' + CUBE_Q.map((q) => `4 ${q.join(' ')}`).join('\n') + '\n0 1\n');
  const le = format === 'binary_little_endian';
  const body = new Uint8Array(8 * 15 + 6 * 17 + 8);
  const dv = new DataView(body.buffer);
  let o = 0;
  CUBE_V.forEach((v, i) => { v.forEach((c) => { dv.setFloat32(o, c, le); o += 4; }); dv.setUint8(o++, i); dv.setUint8(o++, 0); dv.setUint8(o++, 255); });
  CUBE_Q.forEach((q) => { dv.setUint8(o++, 4); q.forEach((k) => { dv.setInt32(o, k, le); o += 4; }); });
  dv.setInt32(o, 0, le); dv.setInt32(o + 4, 1, le);
  return concat(enc(header), body);
}

for (const fmt of ['ascii', 'binary_little_endian', 'binary_big_endian']) {
  test(`PLY ${fmt}: cube with colours, extra element skipped`, () => {
    const r = parsePly({ name: 'cube.ply', bytes: plyCube(fmt) });
    assertContract(r, ['surface']);
    const m = r.models[0];
    close(surfaceVolume(m), 1, 1e-6);
    assert.equal(m.attrs.rgb.length, 24);
    assert.deepEqual([...m.attrs.rgb.slice(3, 6)], [1, 0, 255]);
    assert.equal(r.meta.skipped.edge, 1);
    assert.equal(r.meta.elementTypes.quad4, 6);
    assert.equal(r.meta.encoding, fmt);
  });
}

test('PLY without faces is a point cloud; malformed headers are rejected', () => {
  const text = 'ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nproperty float intensity\nend_header\n0 0 0 0.5\n1 0 0 0.25\n0 2 0 1\n';
  const r = parsePly({ name: 'pts.ply', text });
  assertContract(r, ['points']);
  assert.equal(r.models[0].positions.length, 9);
  assert.deepEqual([...r.models[0].attrs.intensity], [0.5, 0.25, 1]);
  assert.throws(() => parsePly({ name: 'x.ply', text: 'plx\n' }), /first line to be "ply"/);
  assert.throws(() => parsePly({ name: 'x.ply', text: 'ply\nformat ascii 1.0\nelement vertex 1\n' }), /end_header/);
  assert.throws(() => parsePly({ name: 'x.ply', text: 'ply\nformat ascii 1.0\nelement vertex 2\nproperty float x\nproperty float y\nproperty float z\nend_header\n0 0 0\n' }), /ends inside element "vertex"/);
});

test('OFF: counts on the header line, comments, quads and a pentagon', () => {
  const text = `OFF 8 6 12 # cube\n${CUBE_V.map((v) => v.join(' ')).join('\n')}\n# faces\n${CUBE_Q.map((q) => `4 ${q.join(' ')} 255 0 0`).join('\n')}\n`;
  const r = parseOff({ name: 'cube.off', text });
  assertContract(r, ['surface']);
  close(surfaceVolume(r.models[0]), 1, 1e-6);
  assert.equal(r.meta.elementTypes.quad4, 6);
  const c = parseOff({ name: 'p.off', text: 'COFF\n5 1 0\n0 0 0 1 0 0 1\n2 0 0 1 0 0 1\n2 1 0 1 0 0 1\n1 2 0 1 0 0 1\n0 1 0 1 0 0 1\n5 0 1 2 3 4\n' });
  close(surfaceArea(c.models[0]), 3, 1e-9);
  assert.throws(() => parseOff({ name: 'x.off', text: 'cube\n' }), /expected the keyword OFF/);
  assert.throws(() => parseOff({ name: 'x.off', text: '4OFF\n1 0 0\n' }), /not supported/);
});

function gltfCube() {
  // y-up source: the importer must rotate (x, y, z) → (x, −z, y)
  const pos = new Float32Array(CUBE_V.flat());
  const idx = new Uint16Array(CUBE_T.flat());
  const bin = concat(new Uint8Array(pos.buffer), new Uint8Array(idx.buffer));
  const json = {
    asset: { version: '2.0', generator: 'test' }, scene: 0, scenes: [{ nodes: [0] }],
    nodes: [{ name: 'root', translation: [10, 0, 0], children: [1] }, { name: 'part', mesh: 0, scale: [2, 2, 2] }],
    meshes: [{ name: 'cube', primitives: [{ attributes: { POSITION: 0 }, indices: 1 }] }],
    accessors: [{ bufferView: 0, componentType: 5126, count: 8, type: 'VEC3' }, { bufferView: 1, componentType: 5123, count: 36, type: 'SCALAR' }],
    bufferViews: [{ buffer: 0, byteOffset: 0, byteLength: 96 }, { buffer: 0, byteOffset: 96, byteLength: 72 }],
    buffers: [{ byteLength: bin.length }],
  };
  return { json, bin };
}

test('glTF with base64 buffer: node transforms applied, y-up rotated to z-up', () => {
  const { json, bin } = gltfCube();
  json.buffers[0].uri = 'data:application/octet-stream;base64,' + Buffer.from(bin).toString('base64');
  const r = parseGltf({ name: 'cube.gltf', text: JSON.stringify(json) });
  assertContract(r, ['surface']);
  close(surfaceVolume(r.models[0]), 8, 1e-5);
  assert.deepEqual(r.meta.bbox, { min: [10, -2, 0], max: [12, 0, 2] });
  assert.equal(r.meta.units.unit, 'm');
  assert.equal(r.meta.upAxis.source, 'y');
});

test('GLB container and .gltf with side-car .bin', () => {
  const { json, bin } = gltfCube();
  const j = enc(JSON.stringify(json));
  const jp = concat(j, enc(' '.repeat((4 - (j.length % 4)) % 4)));
  const glb = new Uint8Array(12 + 8 + jp.length + 8 + bin.length);
  const dv = new DataView(glb.buffer);
  dv.setUint32(0, 0x46546c67, true); dv.setUint32(4, 2, true); dv.setUint32(8, glb.length, true);
  dv.setUint32(12, jp.length, true); dv.setUint32(16, 0x4e4f534a, true); glb.set(jp, 20);
  dv.setUint32(20 + jp.length, bin.length, true); dv.setUint32(24 + jp.length, 0x004e4942, true); glb.set(bin, 28 + jp.length);
  const r = parseGltf({ name: 'cube.glb', bytes: glb });
  assert.equal(r.meta.format, 'glb');
  close(surfaceVolume(r.models[0]), 8, 1e-5);
  json.buffers[0].uri = 'cube.bin';
  const s = parseGltf({ name: 'cube.gltf', text: JSON.stringify(json), files: new Map([['cube.bin', bin]]) });
  assert.equal(s.meta.counts.triangles, 12);
  assert.throws(() => parseGltf({ name: 'cube.gltf', text: JSON.stringify(json) }), /"cube\.bin" was not supplied/);
  assert.throws(() => parseGltf({ name: 'x.gltf', text: '{"meshes":[]}' }), /missing the "asset" member/);
  json.meshes[0].primitives[0].extensions = { KHR_draco_mesh_compression: {} };
  assert.throws(() => parseGltf({ name: 'cube.gltf', text: JSON.stringify(json), files: new Map([['cube.bin', bin]]) }), /Draco\/meshopt compressed/);
});

const model3mf = `<?xml version="1.0" encoding="UTF-8"?>
<model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">
 <resources>
  <object id="1" type="model" name="cube"><mesh>
   <vertices>${CUBE_V.map((v) => `<vertex x="${v[0] * 1000}" y="${v[1] * 1000}" z="${v[2] * 1000}"/>`).join('')}</vertices>
   <triangles>${CUBE_T.map((t) => `<triangle v1="${t[0]}" v2="${t[1]}" v3="${t[2]}"/>`).join('')}</triangles>
  </mesh></object>
  <object id="2" type="model" name="pair"><components><component objectid="1"/><component objectid="1" transform="1 0 0 0 1 0 0 0 1 0 0 5000"/></components></object>
 </resources>
 <build><item objectid="2" transform="1 0 0 0 1 0 0 0 1 2000 0 0"/></build>
</model>`;

test('3MF package: millimetres to metres, components and build transform', async () => {
  const zip = makeZip([
    { name: '[Content_Types].xml', data: enc('<Types/>') },
    { name: '_rels/.rels', data: enc('<Relationships><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>') },
    { name: '3D/3dmodel.model', data: enc(model3mf) },
  ]);
  const r = await parse3mf({ name: 'pair.3mf', bytes: zip });
  assertContract(r, ['surface']);
  assert.equal(r.meta.counts.triangles, 24);
  close(surfaceVolume(r.models[0]), 2, 1e-5);
  assert.deepEqual(r.meta.bbox, { min: [2, 0, 0], max: [3, 1, 6] });
  assert.equal(r.meta.units.unit, 'mm');
  assert.equal(r.meta.units.basis, 'file');
  const bare = await parse3mf({ name: 'pair.model', text: model3mf });
  assert.equal(bare.meta.counts.triangles, 24);
  await assert.rejects(parse3mf({ name: 'x.3mf', text: '<svg/>' }), /expected root element <model>/);
});
