import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseFluent } from '../../src/geometry/importers/fluent.js';
import * as hdf5 from '../../src/geometry/adapters/hdf5.js';
import { enc, close, concat, assertContract, assertTruncationSafe, surfaceArea, surfaceVolume, CUBE_V, CUBE_Q } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));
const hx = (v) => v.toString(16);
const groupArea = (m, g) => surfaceArea({ positions: m.positions, indices: m.indices.subarray(g.start, g.start + g.count) });
const edgeLength = (m, group) => m.members.filter((e) => !group || e.group === group).reduce((s, e) => s + Math.hypot(m.nodes[3 * e.a] - m.nodes[3 * e.b], m.nodes[3 * e.a + 1] - m.nodes[3 * e.b + 1], m.nodes[3 * e.a + 2] - m.nodes[3 * e.b + 2]), 0);

// ---- hand-built meshes ----------------------------------------------------------------------------------
// A block of nx × ny × nz unit hexahedra. Fluent lists the nodes of a face so that its right-hand normal points
// into the first cell (c0); here every face is written with its normal along +x, +y or +z, so the cell behind
// it is c1 and the far faces of the block have c0 = 0.
function block(nx, ny, nz) {
  const N = (i, j, k) => 1 + i + (nx + 1) * (j + (ny + 1) * k);
  const C = (i, j, k) => (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz ? 0 : 1 + i + nx * (j + ny * k));
  const nodes = [];
  for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) nodes.push([i, j, k]);
  const zones = [
    { id: 2, bc: 2, type: 'interior', name: 'interior-fluid', faces: [] },
    { id: 3, bc: 10, type: 'velocity-inlet', name: 'inlet', faces: [] },
    { id: 4, bc: 5, type: 'pressure-outlet', name: 'outlet', faces: [] },
    { id: 5, bc: 3, type: 'wall', name: 'walls', faces: [] },
  ];
  const put = (ids, c0, c1, zone) => zones[c0 && c1 ? 0 : zone].faces.push({ ids, c0, c1 });
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) put([N(i, j, k), N(i, j + 1, k), N(i, j + 1, k + 1), N(i, j, k + 1)], C(i, j, k), C(i - 1, j, k), i === 0 ? 1 : 2);
  for (let k = 0; k < nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) put([N(i, j, k), N(i, j, k + 1), N(i + 1, j, k + 1), N(i + 1, j, k)], C(i, j, k), C(i, j - 1, k), 3);
  for (let k = 0; k <= nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) put([N(i, j, k), N(i + 1, j, k), N(i + 1, j + 1, k), N(i, j + 1, k)], C(i, j, k), C(i, j, k - 1), 3);
  return { dim: 3, nodes, zones: zones.filter((z) => z.faces.length), cells: nx * ny * nz, cellType: 4 };
}

// nx × ny unit squares in the plane: the faces are the edges
function grid2(nx, ny) {
  const N = (i, j) => 1 + i + (nx + 1) * j;
  const C = (i, j) => (i < 0 || j < 0 || i >= nx || j >= ny ? 0 : 1 + i + nx * j);
  const nodes = [];
  for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) nodes.push([i, j]);
  const zones = [{ id: 2, bc: 2, type: 'interior', name: 'int_fluid', faces: [] }, { id: 3, bc: 3, type: 'wall', name: 'wall', faces: [] }, { id: 4, bc: 7, type: 'symmetry', name: 'sym', faces: [] }];
  const put = (ids, c0, c1, zone) => zones[c0 && c1 ? 0 : zone].faces.push({ ids, c0, c1 });
  for (let j = 0; j < ny; j++) for (let i = 0; i <= nx; i++) put([N(i, j), N(i, j + 1)], C(i - 1, j), C(i, j), 1);
  for (let j = 0; j <= ny; j++) for (let i = 0; i < nx; i++) put([N(i, j), N(i + 1, j)], C(i, j), C(i, j - 1), j === 0 ? 2 : 1);
  return { dim: 2, nodes, zones: zones.filter((z) => z.faces.length), cells: nx * ny, cellType: 3 };
}

const totalFaces = (m) => m.zones.reduce((s, z) => s + z.faces.length, 0);
const zoneText = (m, extra = '') => `(45 (1 fluid fluid)())\n${m.zones.map((z, i) => (i % 2 ? `(39 (${z.id} ${z.type} ${z.name} 1)(\n (roughness . 0.5) (note . "a ) in a string"))\n)` : `(45 (${z.id} ${z.type} ${z.name})())`)).join('\n')}\n${extra}`;

/** ASCII file. `mixed`: face rows start with their node count; `split`: nodes in two zones; `upper`: upper-case hexadecimal. */
function ascii(m, { mixed = false, split = false, upper = false, extra = '' } = {}) {
  const H = upper ? (v) => hx(v).toUpperCase() : hx;
  const n = m.nodes.length, nf = totalFaces(m);
  let s = `(0 "hand-built test mesh")\n(0 grid written by a test\n   nodes:       (10 (id start end type) (x y z ...))\n   faces:       (13 (id start end type elemType)\n                (v-0 v-1 .. v-n right-cell left-cell ...))\n)\n(2 ${m.dim})\n(10 (0 1 ${H(n)} 0 ${m.dim}))\n(13 (0 1 ${H(nf)} 0))\n(12 (0 1 ${H(m.cells)} 0))\n`;
  const cut = split ? Math.floor(n / 2) : n;
  const rows = (a, b) => m.nodes.slice(a, b).map((p) => p.map((v) => v.toExponential(10)).join(' ')).join('\n');
  s += `(10 (6 1 ${H(cut)} 1 ${m.dim})(\n${rows(0, cut)}\n))\n`;
  if (split) s += `(10 (7 ${H(cut + 1)} ${H(n)} 2 ${m.dim})\n(\n${rows(cut, n)}))\n`;
  s += `(12 (1 1 ${H(m.cells)} 1 ${m.cellType}))\n`;
  let first = 1;
  for (const z of m.zones) {
    const k = z.faces[0].ids.length;
    s += `(13 (${H(z.id)} ${H(first)} ${H(first + z.faces.length - 1)} ${H(z.bc)} ${mixed ? 0 : k})(\n${z.faces.map((f) => (mixed ? [f.ids.length] : []).concat(f.ids, [f.c0, f.c1]).map(H).join(' ')).join('\n')}\n))\n`;
    first += z.faces.length;
  }
  return s + zoneText(m, extra);
}

/** Binary sections: `real` 4 or 8 bytes per coordinate, `int` 4 or 8 bytes per integer. */
function binary(m, { real: rb = 8, int: ib = 4, mixed = false, extra = '' } = {}) {
  const n = m.nodes.length, nf = totalFaces(m), sid = rb === 8 ? 3000 : 2000;
  const parts = [enc(`(0 "hand-built test mesh")\n(4 (60 0 0 1 2 4 4 4 8 ${ib} 8 4))\n(2 ${m.dim})\n(10 (0 1 ${hx(n)} 0 ${m.dim}))\n(13 (0 1 ${hx(nf)} 0))\n(12 (0 1 ${hx(m.cells)} 0))\n(${sid + 10} (6 1 ${hx(n)} 1 ${m.dim})\n(`)];
  const xyz = new DataView(new ArrayBuffer(n * m.dim * rb));
  m.nodes.forEach((p, i) => p.forEach((v, a) => (rb === 8 ? xyz.setFloat64((i * m.dim + a) * 8, v, true) : xyz.setFloat32((i * m.dim + a) * 4, v, true))));
  parts.push(new Uint8Array(xyz.buffer), enc(`\n)\nEnd of Binary Section ${sid + 10})\n(${sid + 12} (1 1 ${hx(m.cells)} 1 ${m.cellType})\nEnd of Binary Section   ${sid + 12})\n`));
  let first = 1;
  for (const z of m.zones) {
    const k = z.faces[0].ids.length;
    const ints = z.faces.flatMap((f) => (mixed ? [f.ids.length] : []).concat(f.ids, [f.c0, f.c1]));
    const dv = new DataView(new ArrayBuffer(ints.length * ib));
    ints.forEach((v, i) => (ib === 8 ? dv.setBigUint64(8 * i, BigInt(v), true) : dv.setUint32(4 * i, v, true)));
    parts.push(enc(`(${sid + 13} (${hx(z.id)} ${hx(first)} ${hx(first + z.faces.length - 1)} ${hx(z.bc)} ${mixed ? 0 : k})\n(`), new Uint8Array(dv.buffer), enc(`)\nEnd of Binary Section   ${sid + 13})\n`));
    first += z.faces.length;
  }
  parts.push(enc(zoneText(m, extra)));
  return concat(...parts);
}

test('Fluent ASCII, 3-D: a 3 × 2 × 2 block of hexahedra → boundary surface grouped by face zone', () => {
  const m = block(3, 2, 2);
  const r = parseFluent({ name: 'duct.msh', text: ascii(m, { split: true }) });
  assertContract(r, ['surface']);
  const s = r.models[0];
  assert.equal(s.name, 'duct (boundary)');
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.deepEqual([r.meta.format, r.meta.dimension, r.meta.encoding, r.meta.precision, r.meta.points, r.meta.faces, r.meta.cells], ['fluent', 3, 'ascii', null, 36, 52, 12]);
  // 32 boundary quadrilaterals, 20 faces between two cells
  assert.deepEqual(r.meta.elementTypes, { quad4: 32 });
  assert.deepEqual(r.meta.skipped, { invalidFaces: 0, internalFaces: 20, twoSidedFaces: 0, parentFaces: 0, facesWithoutNodes: 0, cells: 12 });
  assert.deepEqual(r.meta.patches, [{ name: 'inlet', type: 'velocity-inlet', nFaces: 4, id: 3 }, { name: 'outlet', type: 'pressure-outlet', nFaces: 4, id: 4 }, { name: 'walls', type: 'wall', nFaces: 24, id: 5 }]);
  assert.deepEqual(r.meta.zones[0], { id: 2, name: 'interior-fluid', type: 'interior', nFaces: 20, boundary: 0, internal: 20 });
  assert.deepEqual(r.meta.cellZones, [{ id: 1, count: 12, shape: 'hexahedron', name: 'fluid', type: 'fluid' }]);
  assert.deepEqual(s.groups.map((g) => [g.name, g.type, g.count / 3]), [['inlet', 'velocity-inlet', 8], ['outlet', 'pressure-outlet', 8], ['walls', 'wall', 48]]);
  close(groupArea(s, s.groups[0]), 4, 1e-9);
  close(groupArea(s, s.groups[1]), 4, 1e-9);
  close(groupArea(s, s.groups[2]), 24, 1e-9);
  // faces with c0 = 0 and faces with c1 = 0 both come out with their normal pointing out of the mesh
  close(surfaceVolume(s), 12, 1e-9);
  // only the 34 boundary vertices are kept: the 2 interior nodes are gone
  assert.equal(s.positions.length / 3, 34);
  assert.equal(r.meta.counts.triangles, 64);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [3, 2, 2] });
  assert.equal(r.warnings.length, 1);
  assert.match(r.warnings[0], /cells only as the faces that bound them.*12 cells in 1 zone\(s\) and 20 interior faces/);
  assert.deepEqual(r.meta.comments.slice(0, 1), ['hand-built test mesh']);
  assert.deepEqual(r.meta.sections, { '0 comment': 2, '2 dimension': 1, '10 nodes': 3, '13 faces': 5, '12 cells': 2, '45 zone': 3, '39 zone': 2 });
});

test('Fluent ASCII variants give the same surface: mixed rows, upper-case hexadecimal, CR LF, units option', () => {
  const m = block(11, 2, 3);
  const ref = parseFluent({ name: 'a.msh', text: ascii(m) });
  for (const [label, text] of [['mixed', ascii(m, { mixed: true })], ['upper', ascii(m, { upper: true, split: true })], ['crlf', ascii(m).replace(/\n/g, '\r\n')], ['one line', ascii(m, { mixed: true }).replace(/\n/g, ' ')]]) {
    const r = parseFluent({ name: 'a.msh', bytes: enc(text) });
    assert.deepEqual(r.models[0].positions, ref.models[0].positions, label);
    assert.deepEqual(r.models[0].indices, ref.models[0].indices, label);
    assert.deepEqual(r.meta.patches, ref.meta.patches, label);
  }
  close(surfaceArea(ref.models[0]), 2 * (11 * 2 + 2 * 3 + 11 * 3), 1e-9);
  const mm = parseFluent({ name: 'a.msh', text: ascii(m) }, { units: 'mm' });
  assert.deepEqual(mm.meta.units, { unit: 'mm', scale: 0.001, basis: 'option' });
  close(surfaceArea(mm.models[0]), 122e-6, 1e-10);
});

test('Fluent binary sections: single and double precision, 32- and 64-bit integers, fixed and mixed face rows', () => {
  const m = block(3, 2, 2);
  const ref = parseFluent({ name: 'duct.msh', text: ascii(m) });
  for (const o of [{ real: 4, int: 4 }, { real: 8, int: 4 }, { real: 8, int: 8 }, { real: 4, int: 4, mixed: true }, { real: 8, int: 8, mixed: true }]) {
    const r = parseFluent({ name: 'duct.msh', bytes: binary(m, o) });
    assertContract(r, ['surface']);
    assert.deepEqual(r.models[0].positions, ref.models[0].positions, JSON.stringify(o));
    assert.deepEqual(r.models[0].indices, ref.models[0].indices, JSON.stringify(o));
    assert.deepEqual(r.models[0].groups, ref.models[0].groups);
    assert.deepEqual([r.meta.encoding, r.meta.precision, r.meta.integerBytes], ['binary', o.real === 8 ? 'double' : 'single', o.int]);
    assert.deepEqual(r.meta.skipped, ref.meta.skipped);
    assert.deepEqual(r.meta.machine, [60, 0, 0, 1, 2, 4, 4, 4, 8, o.int, 8, 4]);
  }
});

test('Fluent case file: variable, tree, partition and binary sections of other kinds are skipped and counted', () => {
  const m = block(2, 1, 1);
  const flags = concat(enc('(3041 (1 c)\n('), new Uint8Array(48).fill(41), enc(')\nEnd of Binary Section   3041)\n'));
  const extra = '(37 (\n(case-config ((rp-seg? . #t) (separator . #\\() (text . "unbalanced ( in a string \\" and a quote")))\n(list (1 2 (3 4)))))\n(38 ((x . 1)))\n(59 (1 2 3 4)(\n 2 5 6\n))\n(58 (1 2 1 1)(\n 2 1 2\n))\n(61 (1 2)(\n 3 4\n))\n(18 (1 2 3 4)(\n 1 2\n))\n(40 (1 1 2 1)(0 0))\n(64 ())\n';
  const text = ascii(m, { extra });
  const r = parseFluent({ name: 'run.cas', bytes: concat(enc(text), flags) });
  assertContract(r, ['surface']);
  close(surfaceArea(r.models[0]), 10, 1e-9);
  close(surfaceVolume(r.models[0]), 2, 1e-9);
  assert.deepEqual(r.meta.skipped, { invalidFaces: 0, internalFaces: 1, twoSidedFaces: 0, parentFaces: 0, facesWithoutNodes: 0, cells: 2, 'section 37 variables': 1, 'section 38 variables': 1, 'section 59 face tree': 1, 'section 58 cell tree': 1, 'section 61 interface face parents': 1, 'section 18 periodic shadow faces': 1, 'section 40 partition': 1, 'section 64': 1, 'section 3041 node flags': 1 });
  assert.match(r.warnings[0], /^Sections not read: 64, 37 variables, 38 variables, 59 face tree, 58 cell tree, 61 interface face parents, 18 periodic shadow faces, 40 partition, 3041 node flags\.$/);
  // the mesh is in ASCII, one other section is binary
  assert.equal(r.meta.encoding, 'mixed');
});

test('Fluent polygonal faces (face type 5): a cube cell bounded by six octagons', () => {
  // corners plus one node on the middle of every edge; each face lists corner, mid-edge, corner, …
  const nodes = CUBE_V.map((p) => p.slice());
  const mid = new Map();
  const M = (a, b) => { const key = a < b ? `${a}-${b}` : `${b}-${a}`; if (!mid.has(key)) { mid.set(key, nodes.length); nodes.push(CUBE_V[a].map((v, i) => (v + CUBE_V[b][i]) / 2)); } return mid.get(key); };
  // CUBE_Q is wound outward; Fluent wants the normal to point into c0, so the order is reversed
  const faces = CUBE_Q.map((q) => { const w = q.slice().reverse(); return w.flatMap((a, i) => [a + 1, M(a, w[(i + 1) % 4]) + 1]); });
  const text = `(1 "polyhedral test")\n(2 3)\n(10 (1 1 ${hx(nodes.length)} 1 3)(\n${nodes.map((p) => p.join(' ')).join('\n')}\n))\n(12 (2 1 1 1 7))\n(13 (3 1 6 3 5)(\n${faces.map((f) => `8 ${f.map(hx).join(' ')} 1 0`).join('\n')}\n))\n(39 (3 wall skin)())\n(39 (2 fluid cell)())\n`;
  const r = parseFluent({ name: 'poly.msh', text });
  assertContract(r, ['surface']);
  assert.deepEqual(r.meta.elementTypes, { polygon: 6 });
  assert.equal(r.meta.writer, 'polyhedral test');
  assert.equal(r.models[0].indices.length / 3, 36);
  close(surfaceArea(r.models[0]), 6, 1e-9);
  close(surfaceVolume(r.models[0]), 1, 1e-9);
  assert.deepEqual(r.meta.cellZones, [{ id: 2, count: 1, shape: 'polyhedron', name: 'cell', type: 'fluid' }]);
});

test('Fluent 2-D mesh: the boundary edges come back as a frame-kind model with the zones as member groups', () => {
  const m = grid2(4, 3);
  for (const bytes of [enc(ascii(m)), enc(ascii(m, { mixed: true, split: true })), binary(m, { real: 8 }), binary(m, { real: 4, mixed: true })]) {
    const r = parseFluent({ name: 'plate.msh', bytes });
    assertContract(r, ['frame']);
    const f = r.models[0];
    assert.equal(f.name, 'plate (boundary edges)');
    assert.deepEqual([r.meta.dimension, r.meta.points, r.meta.faces, r.meta.cells], [2, 20, 31, 12]);
    assert.deepEqual(r.meta.elementTypes, { line2: 14 });
    assert.deepEqual(r.meta.patches, [{ name: 'wall', type: 'wall', nFaces: 10, id: 3 }, { name: 'sym', type: 'symmetry', nFaces: 4, id: 4 }]);
    assert.equal(r.meta.skipped.internalFaces, 17);
    close(edgeLength(f), 14, 1e-12);
    close(edgeLength(f, 'sym'), 4, 1e-12);
    assert.equal(f.nodes.length / 3, 14);
    for (let i = 2; i < f.nodes.length; i += 3) assert.equal(f.nodes[i], 0);
    assert.match(r.warnings.at(-1), /only the boundary edges are returned \(the mesh has 12 cells/);
  }
});

test('Fluent zone kinds: two-sided walls, parent faces, faces without nodes, zones without a name', () => {
  // two unit cubes side by side; the face between them is a wall (a baffle), not an interior face
  const m = block(2, 1, 1);
  m.zones[0].bc = 3; m.zones[0].type = 'wall'; m.zones[0].name = 'baffle';
  // a parent-face zone (type 31 = 1f) repeating a boundary face, and a zone of virtual faces without nodes
  const f = totalFaces(m);
  const extra = `(13 (9 ${hx(f + 1)} ${hx(f + 1)} 1f 4)(\n1 4 a 7 1 0\n))\n(13 (a ${hx(f + 2)} ${hx(f + 3)} 3f4 0)(\n0 1 0\n0 2 0\n))\n`;
  const text = ascii(m, { extra }).replace('(39 (5 wall walls 1)', '(39 (77 wall other 1)');
  const r = parseFluent({ name: 'baffle.msh', text });
  assertContract(r, ['surface']);
  close(surfaceArea(r.models[0]), 10, 1e-9);
  assert.deepEqual(r.meta.skipped, { invalidFaces: 0, internalFaces: 0, twoSidedFaces: 1, parentFaces: 1, facesWithoutNodes: 2, cells: 2 });
  assert.deepEqual(r.models[0].groups.map((g) => [g.name, g.type]), [['inlet', 'velocity-inlet'], ['outlet', 'pressure-outlet'], ['zone 5', 'wall']]);
  assert.match(r.warnings.join('\n'), /1 face\(s\) lie between two cells in zones that are not of type interior.*baffle \(wall, 1\)\. Import with twoSided = true/);
  assert.match(r.warnings.join('\n'), /2 face\(s\) without nodes/);
  assert.match(r.warnings.join('\n'), /1 parent face\(s\) of a hanging-node mesh were left out/);
  const both = parseFluent({ name: 'baffle.msh', text }, { twoSided: true });
  close(surfaceArea(both.models[0]), 11, 1e-9);
  assert.equal(both.models[0].groups[0].name, 'baffle');
  assert.equal(both.meta.skipped.twoSidedFaces, 0);
  assert.match(both.warnings.join('\n'), /kept as two-sided surfaces: baffle/);
});

test('Fluent surface mesh: faces that name no cell are returned whole', () => {
  const text = `(0 "GAMBIT to Fluent File")\n(2 3)\n(10 (1 1 8 1 3)(\n${CUBE_V.map((p) => p.join(' ')).join('\n')}))\n(13(5 1 6 3 4)(\n${CUBE_Q.map((q) => q.slice().reverse().map((i) => i + 1).join(' ') + ' 0 0').join('\n')}\n))\n(12 (0 1 0 0))\n(45 (5 wall skin)())\n`;
  const r = parseFluent({ name: 'skin.msh', text });
  assertContract(r, ['surface']);
  assert.equal(r.models[0].name, 'skin (surface mesh)');
  assert.equal(r.meta.cells, null);
  close(surfaceArea(r.models[0]), 6, 1e-9);
  close(surfaceVolume(r.models[0]), 1, 1e-9);
  assert.match(r.warnings.at(-1), /surface mesh/);
});

test('Fluent files with defects that are tolerated and reported', () => {
  const m = block(2, 1, 1);
  const good = ascii(m);
  // a zone that holds fewer faces than its header says (the outlet loses its only face → one group less)
  const lines = good.split('\n');
  const at = lines.findIndex((l) => l.startsWith('(13 (4 '));
  const short = lines.slice(0, at + 1).concat(lines.slice(at + 2)).join('\n');
  const r1 = parseFluent({ name: 'short.msh', text: short });
  assert.match(r1.warnings.join('\n'), /face section of zone 4 declares 1 faces but holds 0/);
  assert.match(r1.warnings.join('\n'), /declares b faces|declares 11 faces but its face sections hold 10/);
  close(surfaceArea(r1.models[0]), 9, 1e-9);
  // a face that names a node outside the node sections
  const r2 = parseFluent({ name: 'bad.msh', text: good.replace(/\(13 \(3 [^\n]*\n\S+ /, (s) => s.replace(/\S+ $/, 'ff ')) });
  assert.equal(r2.meta.skipped.invalidFaces, 1);
  assert.match(r2.warnings.join('\n'), /1 boundary face\(s\) that name nodes outside the node sections/);
  close(surfaceArea(r2.models[0]), 9, 1e-9);
  // cell list left open and an empty list appended, as an OpenFOAM-written file has them
  const foam = good.replace('(12 (1 1 2 1 4))', '(12 (1 1 2 1 0)(\n 4 4\n)())');
  const r3 = parseFluent({ name: 'foam.msh', text: foam });
  close(surfaceArea(r3.models[0]), 10, 1e-9);
  assert.equal(r3.meta.cellZones[0].shape, 'mixed');
  assert.equal(r3.warnings.length, 1);
  // a list of variables that never closes: reading goes on at the next mesh section
  const r5 = parseFluent({ name: 'open.cas', text: `(0 "case")\n(2 3)\n(38 ((a . 1) (b (c d)\n${good.slice(good.indexOf('(10 (0'))}` });
  assert.match(r5.warnings[0], /^The parentheses of section 38 at byte 17 do not balance; reading went on at the next mesh section \(byte \d+\)\.$/);
  close(surfaceArea(r5.models[0]), 10, 1e-9);
  assert.equal(r5.meta.patches.length, 3);
  // text between sections
  const r4 = parseFluent({ name: 'stray.msh', text: good.replace('(2 3)\n', '(2 3)\nsome text here\n') });
  assert.match(r4.warnings.join('\n'), /1 run\(s\) of text outside any section/);
  close(surfaceArea(r4.models[0]), 10, 1e-9);
});

test('Fluent: unsupported and unreadable input is refused with a message that says what the file is', () => {
  const m = block(2, 1, 1);
  const good = ascii(m);
  const bad = (input, re, opts) => assert.throws(() => parseFluent(input, opts), (e) => e instanceof Error && !(e instanceof TypeError) && !(e instanceof RangeError) && re.test(e.message), `expected ${re}`);
  bad({ name: 'a.msh.h5', bytes: Uint8Array.of(0x89, 0x48, 0x44, 0x46, 13, 10, 26, 10, 0, 0) }, /a\.msh\.h5: this is a Fluent file in the HDF5-based common fluids format.*cff-files\? no/);
  bad({ name: 'a.cas.gz', bytes: Uint8Array.of(0x1f, 0x8b, 8, 0, 0, 0) }, /still gzip-compressed/);
  bad({ name: 'a.msh', text: '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n' }, /it is a Gmsh mesh, not a Fluent mesh/);
  bad({ name: 'a.msh', text: '  \n' }, /the file is empty/);
  bad({ name: 'a.msh', text: 'hello (2 3)' }, /not a Fluent mesh or case file \(expected the file to start with a parenthesised section/);
  bad({ name: 'a.msh', text: '(define x 1)' }, /the first list does not start with a section number/);
  bad({ name: 'a.msh', text: '(0 "only a comment")\n(2 3)\n' }, /no node section \(10\) was found \(sections present: 0 comment, 2 dimension\)/);
  bad({ name: 'a.dat', text: '(0 "fluent data")\n(2 3)\n(300 (1 1 1 0 0 1 2)(\n1.0\n2.0))\n' }, /looks like a Fluent data file \(\.dat\)/);
  // meshio writes cells by their nodes and no faces: nothing Fluent itself could read either
  bad({ name: 'a.msh', text: '(1 "meshio 5.3.5")\n(2 3)\n(10 (1 1 4 1 3)(\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n))\n(12 (1 1 1 1 2))\n' }, /holds 4 nodes but no face section \(13\)/);
  bad({ name: 'a.msh', text: '(2 4)\n' }, /dimension section \(2 …\) at byte 0: expected 2 or 3/);
  bad({ name: 'a.msh', text: good.replace('(13 (2 ', '(13 (2 zz ') }, /section 13 at byte \d+: the header list holds something other than hexadecimal integers/);
  bad({ name: 'a.msh', text: good.replace(/\(13 \(3 ([^)]*) 4\)/, '(13 (3 $1 9)') }, /face type 9 is not one of 0 \(mixed\), 2 \(line\)/);
  bad({ name: 'a.msh', text: good.replace('0.0000000000e+0 0.0000000000e+0 0.0000000000e+0', '0.0 0.0 abc') }, /node section of zone 6: "abc[^"]*" at byte \d+ is not a number/);
  bad({ name: 'a.msh', text: good.replace(/\(10 \(6 1 c /, '(10 (6 1 d ') }, /declares 13 nodes but holds 36 coordinates \(expected 26 or 39\)/);
  bad({ name: 'a.msh', text: good.replace(/\(10 \(6 1 c /, '(10 (6 1 ffffffffff ') }, /numbers its nodes up to 1099511627775, more than a file of \d+ bytes can hold/);
  bad({ name: 'a.msh', text: good.replace(/\(13 \(5 (\S+) \S+ /, '(13 (5 $1 fffffff ') }, /face section of zone 5 declares \d+ faces, more than the \d+ bytes left in the file can hold/);
  bad({ name: 'a.msh', text: good.replace(/(\(13 \(5 [^\n]*\n)(\S+)/, '$1$2g') }, /face section of zone 5: unexpected character "g" at byte \d+/);
  bad({ name: 'a.msh', text: good.replace(/(\(13 \(5 [^\n]*\n[^\n]*\n)/, '$11 2\n') }, /face section of zone 5 declares a faces but|more data follows them|closes in the middle of face/);
  bad({ name: 'a.msh', text: good.slice(0, good.indexOf('(13 (3')) }, /lie between two cells|no boundary faces could be built/);
  // binary: body shorter than the header says, and a body that is neither 32- nor 64-bit
  const bin = binary(m, { real: 8, int: 4 });
  const cutAt = (bytes, needle, plus) => { const s = new TextDecoder('latin1').decode(bytes); return bytes.slice(0, s.indexOf(needle) + plus); };
  bad({ name: 'b.msh', bytes: cutAt(bin, '(3010 (6', 40) }, /binary node section of zone 6 declares 12 nodes \(288 bytes\) but the file ends after \d+ bytes of it: the file is truncated/);
  bad({ name: 'b.msh', bytes: cutAt(bin, '(3013 (5', 60) }, /binary face section of zone 5 declares 8 faces but the file ends after \d+ bytes of it: the file is truncated/);
  const s = new TextDecoder('latin1').decode(bin);
  const i = s.indexOf('(3013 (5');
  const odd = concat(bin.slice(0, i + 40), Uint8Array.of(1, 2, 3), bin.slice(i + 40));
  bad({ name: 'b.msh', bytes: odd }, /binary face section of zone 5 at byte \d+: the end marker is not where 8 faces of 32-bit or of 64-bit integers end/);
  const j = s.indexOf('(3010 (6');
  bad({ name: 'b.msh', bytes: concat(bin.slice(0, j + 30), Uint8Array.of(1, 2, 3), bin.slice(j + 30)) }, /binary node section of zone 6 at byte \d+: the end marker is not where 12 nodes of 3 double-precision coordinates end/);
  bad({ name: 'b.cas', bytes: concat(enc(good), enc('(3041 (1 c)\n('), new Uint8Array(40)) }, /the file ends inside binary section 3041 that starts at byte \d+ \(no "End of Binary Section" marker follows/);
  bad({ name: 'b.cas', bytes: enc(`${good}(37 (\n(a "open string`) }, /the file ends inside section 37 that starts at byte \d+ \(unbalanced parentheses: the file is truncated\)/);
});

test('Fluent: truncated and corrupted files never fail with an internal error', async () => {
  const m3 = block(3, 2, 2), m2 = grid2(4, 3);
  // a result with indices in range, or an Error worded for the user that names the file
  const attempt = (b, label) => {
    try {
      const r = parseFluent({ name: 'damaged.msh', bytes: b });
      assert.ok(Array.isArray(r.models) && r.models.length > 0 && r.meta.counts, `${label}: result shape`);
      for (const mdl of r.models) if (mdl.kind === 'surface') { const nv = mdl.positions.length / 3; let hi = 0; for (let i = 0; i < mdl.indices.length; i++) if (mdl.indices[i] > hi) hi = mdl.indices[i]; assert.ok(hi < nv, `${label}: index in range`); }
    } catch (e) {
      assert.ok(e instanceof Error && !(e instanceof TypeError) && !(e instanceof RangeError) && !(e instanceof ReferenceError) && e.code !== 'ERR_ASSERTION' && /^damaged\.msh: .{15,}/.test(e.message), `${label}: ${e && e.constructor.name}: ${e && e.message}\n${e && e.stack}`);
    }
  };
  const samples = [
    ['ascii 3-D', enc(ascii(m3, { split: true, extra: '(37 ((a . "b")))\n(59 (1 2 3 4)(\n 2 5 6\n))\n' }))],
    ['ascii mixed 2-D', enc(ascii(m2, { mixed: true }))],
    ['binary double', binary(m3, { real: 8, int: 4 })],
    ['binary single mixed', binary(m3, { real: 4, int: 4, mixed: true })],
    ['binary 64-bit', binary(m2, { real: 8, int: 8, mixed: true })],
    ['real binary msh', real('fluent-hex-binary.msh')],
    ['real binary cas', real('fluent-cube-5tets-binary.cas')],
    ['real ascii cas', real('fluent-wedge-cylinder.cas')],
  ];
  for (const [label, bytes] of samples) {
    await assertTruncationSafe(bytes, (b) => parseFluent({ name: 'cut.msh', bytes: b }), label);
    // every prefix of the two smallest files, every third of the other small ones
    if (bytes.length < 4000) for (let c = 0; c < bytes.length; c += bytes.length < 1000 ? 1 : 3) attempt(bytes.subarray(0, c), `${label} cut at ${c}`);
    let seed = 12345 + bytes.length;
    const rnd = () => { seed = (Math.imul(seed, 1103515245) + 12345) >>> 0; return seed / 4294967296; };
    for (let round = 0; round < 90; round++) {
      let b = bytes.slice();
      const kind = round % 3;
      const at = Math.floor(rnd() * b.length);
      if (kind === 0) for (let k = 0; k < 1 + round % 4; k++) b[Math.floor(rnd() * b.length)] = [40, 41, 34, 48, 102, 32, 10, 255, 0, 69][Math.floor(rnd() * 10)];
      else if (kind === 1) b = concat(b.slice(0, at), b.slice(at + 1 + Math.floor(rnd() * 9)));
      else b = concat(b.slice(0, at), Uint8Array.from({ length: 1 + Math.floor(rnd() * 6) }, () => Math.floor(rnd() * 256)), b.slice(at));
      attempt(b, `${label} round ${round}`);
    }
  }
  // random bytes behind a plausible start
  let seed = 99;
  const rnd = () => { seed = (Math.imul(seed, 1103515245) + 12345) >>> 0; return seed / 4294967296; };
  for (let round = 0; round < 200; round++) {
    const junk = Uint8Array.from({ length: 20 + Math.floor(rnd() * 400) }, () => (round % 2 ? Math.floor(rnd() * 256) : '()0123af "\n'.charCodeAt(Math.floor(rnd() * 11))));
    assert.throws(() => parseFluent({ name: 'junk.msh', bytes: concat(enc(['(0 "x")\n(2 3)\n', '(10 (1 1 8 1 3)(', '(3013 (2 1 8 2 0)\n(', '(13 (2 1 8 2 0)('][round % 4]), junk) }), (e) => e instanceof Error && !(e instanceof TypeError) && !(e instanceof RangeError) && /^junk\.msh: .{15,}/.test(e.message));
  }
});

// ---- Fluent's HDF5 files (.msh.h5, .cas.h5), read by the HDF5 adapter with the same boundary builder --------
// The tree below follows files written by Fluent 2021 R2 – 2023 R1: a case file keeps all faces in one section,
// lists the boundary zones first and stores c1 for the range of interior faces only.

/** Write a mesh of `block` / `grid2` as an HDF5 file in Fluent's layout. */
async function cff(m, { names = true } = {}) {
  const { h5, FS } = await hdf5.load();
  const path = `/fluent-fixture-${m.nodes.length}-${names}.h5`;
  const f = new h5.File(path, 'w');
  const u64 = (a) => BigUint64Array.from(a, (v) => BigInt(v));
  const range = (o, lo, hi) => { o.create_attribute('minId', lo); o.create_attribute('maxId', hi); return o; };
  try {
    const mesh = f.create_group('meshes').create_group('1');
    const zones = m.zones.filter((z) => z.bc !== 2).concat(m.zones.filter((z) => z.bc === 2));
    const faces = zones.flatMap((z) => z.faces);
    const nInt = faces.filter((x) => x.c0 && x.c1).length, nBnd = faces.length - nInt;
    for (const [k, v] of [['dimension', m.dim], ['nodeCount', m.nodes.length], ['faceCount', faces.length], ['cellCount', m.cells]]) mesh.create_attribute(k, v);
    range(mesh.create_group('nodes').create_group('coords').create_dataset({ name: '1', data: Float64Array.from(m.nodes.flat()), shape: [m.nodes.length, m.dim] }), 1, m.nodes.length);
    const fg = mesh.create_group('faces');
    const sec = range(fg.create_group('nodes').create_group('1'), 1, faces.length);
    sec.create_dataset({ name: 'nnodes', data: Int16Array.from(faces, (x) => x.ids.length) });
    sec.create_dataset({ name: 'nodes', data: Uint32Array.from(faces.flatMap((x) => x.ids)) });
    // the boundary faces of the far side have their cell in c1: stored here as c0, with the order reversed
    range(fg.create_group('c0').create_dataset({ name: '1', data: Uint32Array.from(faces, (x) => x.c0 || x.c1) }), 1, faces.length);
    range(fg.create_group('c1').create_dataset({ name: '1', data: Uint32Array.from(faces.slice(nBnd), (x) => x.c1) }), nBnd + 1, faces.length);
    const zt = fg.create_group('zoneTopology');
    let first = 1;
    const lo = zones.map((z) => { const v = first; first += z.faces.length; return v; });
    zt.create_dataset({ name: 'id', data: Int32Array.from(zones, (z) => z.id) });
    zt.create_dataset({ name: 'zoneType', data: Int32Array.from(zones, (z) => z.bc) });
    zt.create_dataset({ name: 'minId', data: u64(lo) });
    zt.create_dataset({ name: 'maxId', data: u64(zones.map((z, i) => lo[i] + z.faces.length - 1)) });
    if (names) zt.create_dataset({ name: 'name', data: [zones.map((z) => z.name).join(';')] });
    const ct = mesh.create_group('cells').create_group('zoneTopology');
    ct.create_dataset({ name: 'id', data: Int32Array.of(1) });
    ct.create_dataset({ name: 'minId', data: u64([1]) });
    ct.create_dataset({ name: 'maxId', data: u64([m.cells]) });
    ct.create_dataset({ name: 'cellType', data: Int32Array.of(m.cellType) });
    ct.create_dataset({ name: 'zoneType', data: Int32Array.of(17) });
    ct.create_dataset({ name: 'name', data: ['steel'] });
    f.create_group('settings').create_dataset({ name: 'Origin', data: ['ANSYS_FLUENT 22.2 Build 10212'] });
  } finally { f.close(); }
  const bytes = FS.readFile(path).slice();
  FS.unlink(path);
  return bytes;
}

test('Fluent HDF5 (.msh.h5 / .cas.h5): the HDF5 adapter returns the same boundary as the legacy file', { skip: !hdf5.available() && 'the vendored HDF5 library is not present' }, async () => {
  const m = block(3, 2, 2);
  // the far faces carry their cell in c0 in this file, so their node order is the reverse of the ASCII one
  const flipped = { ...m, zones: m.zones.map((z) => ({ ...z, faces: z.faces.map((x) => (x.c0 ? x : { ids: x.ids.slice().reverse(), c0: x.c1, c1: 0 })) })) };
  const ref = parseFluent({ name: 'duct.msh', text: ascii(m) });
  const r = await hdf5.parse({ name: 'duct.cas.h5', bytes: await cff(flipped) });
  assertContract(r, ['surface']);
  assert.equal(r.models[0].name, 'duct (boundary)');
  assert.deepEqual([r.meta.format, r.meta.flavour, r.meta.encoding, r.meta.writer, r.meta.dimension], ['fluent', 'fluent', 'hdf5', 'ANSYS_FLUENT 22.2 Build 10212', 3]);
  assert.deepEqual([r.meta.points, r.meta.faces, r.meta.cells], [36, 52, 12]);
  assert.deepEqual(r.meta.skipped, ref.meta.skipped);
  assert.deepEqual(r.meta.patches, ref.meta.patches);
  assert.deepEqual(r.meta.cellZones, [{ id: 1, count: 12, shape: 'hexahedron', name: 'steel', type: 'solid' }]);
  assert.deepEqual(r.models[0].groups, ref.models[0].groups);
  close(surfaceArea(r.models[0]), 32, 1e-9);
  close(surfaceVolume(r.models[0]), 12, 1e-9);
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  // without the name table the zones are numbered, the types still come from zoneType
  const anon = await hdf5.parse({ name: 'duct.msh.h5', bytes: await cff(flipped, { names: false }) }, { units: 'mm' });
  assert.deepEqual(anon.models[0].groups.map((g) => [g.name, g.type]), [['zone 3', 'velocity-inlet'], ['zone 4', 'pressure-outlet'], ['zone 5', 'wall']]);
  close(surfaceArea(anon.models[0]), 32e-6, 1e-10);
  // 2-D
  const e = await hdf5.parse({ name: 'plate.msh.h5', bytes: await cff(grid2(4, 3)) });
  assertContract(e, ['frame']);
  assert.equal(e.models[0].name, 'plate (boundary edges)');
  assert.deepEqual(e.meta.patches, [{ name: 'wall', type: 'wall', nFaces: 10, id: 3 }, { name: 'sym', type: 'symmetry', nFaces: 4, id: 4 }]);
  close(edgeLength(e.models[0]), 14, 1e-12);
  assert.equal(e.meta.skipped.internalFaces, 17);
});

// ---- files written by other programs (tests/geometry/fixtures/real, sources in its README) ----------------
// Node, face and zone counts and the areas below were confirmed with an independent reading of each file.

test('Fluent real files: PETSc cube of five tetrahedra, ASCII and binary (single precision) case files', () => {
  const a = parseFluent({ name: 'cube_5tets_ascii.cas', bytes: real('fluent-cube-5tets-ascii.cas') });
  const b = parseFluent({ name: 'cube_5tets.cas', bytes: real('fluent-cube-5tets-binary.cas') });
  for (const r of [a, b]) {
    assertContract(r, ['surface']);
    assert.deepEqual([r.meta.points, r.meta.faces, r.meta.cells, r.meta.writer], [8, 16, 5, 'cube_5tets test mesh']);
    assert.deepEqual(r.meta.elementTypes, { tri3: 12 });
    assert.equal(r.meta.skipped.internalFaces, 4);
    assert.deepEqual(r.meta.patches, [{ name: 'L1', type: 'wall', nFaces: 12, id: 12 }]);
    assert.deepEqual(r.meta.cellZones, [{ id: 10, count: 5, shape: 'mixed', name: 'Volume_1', type: 'fluid' }]);
    // a cube of side 9.999983
    close(surfaceArea(r.models[0]), 6 * 9.999983 ** 2, 1e-3);
    close(surfaceVolume(r.models[0]), 9.999983 ** 3, 1e-2);
  }
  assert.deepEqual([a.meta.encoding, a.meta.precision], ['ascii', null]);
  assert.deepEqual([b.meta.encoding, b.meta.precision, b.meta.integerBytes], ['binary', 'single', 4]);
  assert.deepEqual(a.models[0].positions, b.models[0].positions);
  assert.deepEqual(a.models[0].indices, b.models[0].indices);
});

test('Fluent real files: hexahedral mesh in binary double precision, wedge mesh in ASCII from Fluent 24.2', () => {
  const h = parseFluent({ name: 'hex.msh', bytes: real('fluent-hex-binary.msh') });
  assertContract(h, ['surface']);
  assert.deepEqual([h.meta.points, h.meta.faces, h.meta.cells, h.meta.encoding, h.meta.precision, h.meta.integerBytes], [343, 756, 216, 'binary', 'double', 4]);
  assert.deepEqual(h.meta.patches, [{ name: 'WALL_SOLID', type: 'wall', nFaces: 216, id: 11 }]);
  assert.deepEqual(h.meta.cellZones, [{ id: 9, count: 216, shape: 'hexahedron', name: 'SOLID', type: 'solid' }]);
  assert.deepEqual([h.meta.skipped.internalFaces, h.meta.elementTypes.quad4], [540, 216]);
  assert.deepEqual(h.meta.comments.slice(0, 2), ['Created by: Fluent Interface 23.1', 'Units: Meters']);
  // a block of 32 mm × 32 mm × 32 mm
  close(surfaceArea(h.models[0]), 6 * 0.032 ** 2, 1e-9);
  close(surfaceVolume(h.models[0]), 0.032 ** 3, 1e-11);
  assert.equal(h.models[0].positions.length / 3, 218);

  const w = parseFluent({ name: 'wedge_cylinder.cas', bytes: real('fluent-wedge-cylinder.cas') });
  assertContract(w, ['surface']);
  assert.deepEqual([w.meta.points, w.meta.faces, w.meta.cells, w.meta.encoding], [104, 264, 84, 'ascii']);
  assert.deepEqual(w.meta.patches.map((p) => [p.name, p.type, p.nFaces]), [['inlet', 'wall', 12], ['outlet', 'wall', 12], ['wall', 'wall', 42]]);
  assert.deepEqual(w.meta.elementTypes, { quad4: 54, tri3: 12 });
  assert.deepEqual([w.meta.skipped.internalFaces, w.meta.skipped['section 41 node flags']], [198, 1]);
  const g = w.models[0].groups;
  close(groupArea(w.models[0], g[0]), 0.000419043716, 1e-10);
  close(groupArea(w.models[0], g[1]), 0.000419043716, 1e-10);
  close(groupArea(w.models[0], g[2]), 0.003244595989, 1e-10);
  close(surfaceVolume(w.models[0]), 1.78428813e-5, 1e-11);
  assert.deepEqual(w.meta.comments[0], 'fluent24.2.0  build-id: 10192');
});

test('Fluent real files: 2-D meshes — CUBIT triangles in ASCII, a binary case file with four boundary zones', () => {
  const c = parseFluent({ name: 'square.cas', bytes: real('fluent-cubit-square.cas') });
  assertContract(c, ['frame']);
  assert.deepEqual([c.meta.dimension, c.meta.points, c.meta.faces, c.meta.cells], [2, 92, 241, 150]);
  assert.deepEqual(c.meta.patches.map((p) => [p.name, p.nFaces]), [['wall_0', 32]]);
  assert.equal(c.meta.skipped.internalFaces, 209);
  assert.equal(c.meta.cellZones[0].shape, 'triangle');
  // the unit square
  close(edgeLength(c.models[0]), 4, 1e-9);
  assert.deepEqual(c.meta.bbox, { min: [-0.5, -0.5, 0], max: [0.5, 0.5, 0] });

  const s = parseFluent({ name: 'pro_species.cas', bytes: real('fluent-species-2d.cas') });
  assertContract(s, ['frame']);
  assert.equal(s.models[0].name, 'pro_species (boundary edges)');
  assert.deepEqual([s.meta.dimension, s.meta.points, s.meta.faces, s.meta.cells, s.meta.precision], [2, 861, 1660, 800, 'double']);
  assert.deepEqual(s.meta.patches.map((p) => [p.name, p.type, p.nFaces]), [['PO', 'wall', 20], ['VI-OXID', 'wall', 10], ['VI-FUEL', 'wall', 10], ['SYMMETRY', 'symmetry', 80]]);
  assert.equal(s.meta.skipped.internalFaces, 1540);
  // a 40 mm × 10 mm rectangle: outlet 10 mm, two inlets of 5 mm, two symmetry sides of 20 mm
  close(edgeLength(s.models[0], 'PO'), 0.01, 1e-9);
  close(edgeLength(s.models[0], 'VI-OXID'), 0.005, 1e-9);
  close(edgeLength(s.models[0], 'VI-FUEL'), 0.005, 1e-9);
  close(edgeLength(s.models[0], 'SYMMETRY'), 0.04, 1e-9);
});
