import test from 'node:test';
import assert from 'node:assert/strict';
import { deflateSync } from 'node:zlib';
import { parseGmsh } from '../../src/geometry/importers/gmsh.js';
import { parseNastran } from '../../src/geometry/importers/nastran.js';
import { parseAbaqus, abaqusType } from '../../src/geometry/importers/abaqus.js';
import { parseAnsysCdb } from '../../src/geometry/importers/ansys.js';
import { parseUnv } from '../../src/geometry/importers/unv.js';
import { parseVtk, parseVtkXml } from '../../src/geometry/importers/vtk.js';
import { parseSu2, parseOpenFoam } from '../../src/geometry/importers/cfd.js';
import { enc, close, closeArr, concat, assertContract, REF, CANON_EDGES, withMids, checkElements, assertTruncationSafe } from './fixtures/helpers.js';

const byKind = (r, k) => r.models.find((m) => m.kind === k);
/** Node table builder: returns the 1-based ids of the points added. */
function nodeTable() {
  const pts = [];
  return { pts, add: (list) => list.map((p) => pts.push(p)) };
}

// ---- Gmsh ------------------------------------------------------------------------------------------------

// native Gmsh mid-side edge order (Gmsh reference manual, "Node ordering")
const GMSH_TET10 = [[0, 1], [1, 2], [2, 0], [3, 0], [2, 3], [1, 3]];
const GMSH_HEX20 = [[0, 1], [0, 3], [0, 4], [1, 2], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 7], [5, 6], [6, 7]];

test('Gmsh MSH 2.2: mixed elements, node-order permutation, physical groups', () => {
  const T = nodeTable();
  const hex = T.add(withMids(REF.hex, GMSH_HEX20));
  const tet = T.add(withMids(REF.tet, GMSH_TET10, [2, 0, 0]));
  const wedge = T.add(withMids(REF.wedge, [], [4, 0, 0]));
  const pyr = T.add(withMids(REF.pyr, [], [6, 0, 0]));
  const beam = T.add([[0, 0, 5], [0, 0, 9]]);
  const els = [
    `1 17 2 7 1 ${hex.join(' ')}`, `2 11 2 7 1 ${tet.join(' ')}`, `3 6 2 8 2 ${wedge.join(' ')}`, `4 7 2 8 2 ${pyr.join(' ')}`,
    `5 3 2 3 9 ${hex.slice(0, 4).join(' ')}`, `6 2 2 4 9 ${tet.slice(0, 3).join(' ')}`, `7 1 2 5 1 ${beam.join(' ')}`, `8 15 2 0 1 ${beam[0]}`,
  ];
  const text = `$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n3 7 "steel"\n3 8 "grout"\n2 3 "inlet"\n1 5 "brace"\n$EndPhysicalNames\n$Nodes\n${T.pts.length}\n${T.pts.map((p, i) => `${i + 1} ${p.join(' ')}`).join('\n')}\n$EndNodes\n$Elements\n${els.length}\n${els.join('\n')}\n$EndElements\n$NodeData\n1\n"t"\n$EndNodeData\n`;
  const r = parseGmsh({ name: 'mix.msh', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { hex20: 1, tet10: 1, wedge6: 1, pyr5: 1, quad4: 1, tri3: 1, line2: 1 });
  const vol = byKind(r, 'volume');
  checkElements(vol);
  assert.deepEqual(vol.cells.map((c) => c.type).sort(), ['hex20', 'pyr5', 'tet10', 'wedge6']);
  assert.deepEqual(Object.values(vol.tagNames).sort(), ['grout', 'steel']);
  assert.deepEqual([...vol.sets.cells.steel], [0, 1]);
  const surf = byKind(r, 'surface');
  assert.deepEqual(surf.groups.map((g) => [g.name, g.count]), [['inlet', 6], ['4', 3]]);
  const frame = byKind(r, 'frame');
  assert.equal(frame.members[0].group, 'brace');
  close(frame.nodes[5] - frame.nodes[2], 4, 1e-12);
  assert.equal(r.meta.version, 2.2);
  assert.ok(r.warnings.some((w) => /\$NodeData/.test(w)));
});

test('Gmsh MSH 4.1: entity blocks, physical tags through $Entities', () => {
  const text = `$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n2\n2 1 "wall"\n3 2 "fluid"\n$EndPhysicalNames\n$Entities\n0 0 1 1\n5 0 0 0 1 1 0 1 1 0\n9 0 0 0 1 1 1 1 2 0\n$EndEntities\n$Nodes\n1 4 1 4\n3 9 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n$Elements\n2 2 1 2\n2 5 2 1\n1 1 3 2\n3 9 4 1\n2 1 2 3 4\n$EndElements\n`;
  const r = parseGmsh({ name: 't.msh', text }, { units: 'mm' });
  assertContract(r, ['surface', 'volume']);
  checkElements(byKind(r, 'volume'));
  assert.deepEqual(byKind(r, 'volume').tagNames, { 3000002: 'fluid' });
  assert.equal(byKind(r, 'surface').groups[0].name, 'wall');
  close(r.meta.bbox.max[0], 0.001, 1e-9);
  assert.deepEqual(r.meta.elementSets.sort(), ['fluid', 'wall']);
});

test('Gmsh: binary, unsupported version and non-mesh input are rejected with the remedy', () => {
  assert.throws(() => parseGmsh({ name: 'b.msh', text: '$MeshFormat\n4.1 1 8\n' }), /binary Gmsh files are not read\. Re-save as ASCII/);
  assert.throws(() => parseGmsh({ name: 'b.msh', text: '$MeshFormat\n4 0 8\n$EndMeshFormat\n' }), /version 4 is not read/);
  assert.throws(() => parseGmsh({ name: 'b.msh', text: 'hello\n' }), /expected the first line to be \$MeshFormat/);
  assert.throws(() => parseGmsh({ name: 'b.msh', text: '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n1 0 0 0\n' }), /file ends inside Nodes/);
});

// ---- Nastran ---------------------------------------------------------------------------------------------

const sf = (...f) => f.map((v, i) => (i === 0 ? String(v).padEnd(8) : String(v).padStart(8))).join('');
// Nastran CHEXA: G9–G12 bottom edges, G13–G16 vertical edges, G17–G20 top edges
const NAS_HEX20 = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 5], [2, 6], [3, 7], [4, 5], [5, 6], [6, 7], [7, 4]];

test('Nastran BDF: small, large and free field, continuations, properties and materials', () => {
  const T = nodeTable();
  const hex = T.add(withMids(REF.hex, NAS_HEX20).map((p) => p.map((c) => c * 1000)));
  const tet = T.add(withMids(REF.tet, CANON_EDGES.tet10, [2, 0, 0]).map((p) => p.map((c) => c * 1000)));
  const n0 = T.pts.length;
  const lines = ['SOL 101', 'CEND', 'SUBCASE 1', '  LOAD = 1', 'BEGIN BULK', '$ comment'];
  T.pts.forEach((p, i) => lines.push(sf('GRID', i + 1, '', p[0].toFixed(1), p[1].toFixed(1), p[2].toFixed(1))));
  // large field, free field, and Nastran's exponent-without-E real
  lines.push(`GRID*   ${String(n0 + 1).padStart(16)}${''.padStart(16)}${'0.0'.padStart(16)}${'0.0'.padStart(16)}*`, `*       ${'5.+3'.padStart(16)}`);
  lines.push(`GRID,${n0 + 2},,0.0,0.0,9000.`, `GRID,${n0 + 3},,4000.,0.0,9000.`, `GRID,${n0 + 4},,4000.,3000.,9000.`, `GRID,${n0 + 5},,0.,3000.,9000.`);
  lines.push(sf('CBAR', 1, 10, n0 + 1, n0 + 2, '0.0', '1.0', '0.0'), sf('CROD', 2, 11, n0 + 2, n0 + 3));
  lines.push(sf('CQUAD4', 3, 20, n0 + 2, n0 + 3, n0 + 4, n0 + 5), `CTRIA3,4,20,${n0 + 2},${n0 + 3},${n0 + 4}`);
  lines.push(sf('CTETRA', 5, 30, ...tet.slice(0, 6)) + '+T1', sf('+T1', ...tet.slice(6)));
  lines.push(sf('CHEXA', 6, 30, ...hex.slice(0, 6)), sf('', ...hex.slice(6, 14)), sf('', ...hex.slice(14)));
  lines.push(sf('PBARL', 10, 1, '', 'TUBE'), sf('', '305.', '286.'), sf('PTUBE', 11, 1, '406.4', '12.7'), sf('PSHELL', 20, 1, '25.0'), sf('PSOLID', 30, 1));
  lines.push(sf('MAT1', 1, '2.1+5', '', '0.3', '7.85-9'), sf('SPC1', 1, 123456, 1), sf('FORCE', 1, 2, '', '1.0', '0.0', '0.0', '-1.0'), 'ENDDATA', 'GRID,999,,1.+9,0.,0.');
  const r = parseNastran({ name: 'model.bdf', bytes: enc(lines.join('\n')) }, { units: 'mm' });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { line2: 2, quad4: 1, tri3: 1, tet10: 1, hex20: 1 });
  const vol = byKind(r, 'volume');
  checkElements(vol);
  assert.equal(vol.tagNames[30], 'PSOLID 30');
  const frame = byKind(r, 'frame');
  closeArr([...frame.nodes], [0, 0, 5, 0, 0, 9, 4, 0, 9], 1e-12);
  assert.equal(frame.members[0].section.type, 'tube');
  closeArr([frame.members[0].section.D, frame.members[0].section.t], [0.61, 0.019], 1e-12);
  assert.deepEqual(frame.members[0].orient, [0, 1, 0]);
  close(frame.members[1].section.D, 0.4064, 1e-12);
  assert.equal(frame.members[0].material, 'MAT1 1');
  const surf = byKind(r, 'surface');
  close(surf.groups[0].thickness, 0.025, 1e-12);
  assert.equal(surf.elements.length, 2);
  close(r.meta.materials['MAT1 1'].E, 2.1e5, 1e-6);
  close(r.meta.materials['MAT1 1'].rho, 7.85e-9, 1e-18);
  assert.deepEqual(r.meta.skipped.cards, { SPC1: 1, FORCE: 1 });
  assert.equal(r.meta.bbox.max[0], 4, 'the GRID after ENDDATA is not read');
});

test('Nastran: cylindrical coordinate system, missing bulk data', () => {
  const text = ['CORD2C,5,0,10.,0.,0.,10.,0.,1.,+', '+,11.,0.,0.', 'GRID,1,5,2.,90.,3.', 'GRID,2,,0.,0.,0.', 'CROD,1,1,1,2'].join('\n');
  const r = parseNastran({ name: 'c.bdf', text });
  closeArr([...r.models[0].nodes].slice(0, 3), [10, 2, 3], 1e-9);
  assert.throws(() => parseNastran({ name: 'e.bdf', text: '$ nothing here\n' }), /no bulk-data cards were found/);
});

// ---- Abaqus ----------------------------------------------------------------------------------------------

test('Abaqus INP: sets, sections, materials, wrapped C3D20, 3-node beam ordering', () => {
  const T = nodeTable();
  const hex = T.add(withMids(REF.hex, CANON_EDGES.hex20));
  const tet = T.add(withMids(REF.tet, CANON_EDGES.tet10, [2, 0, 0]));
  const b = T.add([[0, 0, 5], [0, 0, 7], [0, 0, 9], [4, 0, 9], [4, 3, 9], [0, 3, 9]]);
  const text = [
    '*Heading', ' test deck', '** comment', '*Node, nset=ALL',
    ...T.pts.map((p, i) => `${i + 1}, ${p.join(', ')}`),
    '*Element, type=C3D20R, elset=BLOCK', `1, ${hex.slice(0, 15).join(', ')},`, ` ${hex.slice(15).join(', ')}`,
    '*Element, type=C3D10, elset=BLOCK', `2, ${tet.join(', ')}`,
    '*Element, type=B32, elset=LEG', `3, ${b[0]}, ${b[1]}, ${b[2]}`,
    '*Element, type=B31, elset=BRACE', `4, ${b[2]}, ${b[3]}`,
    '*Element, type=S4R, elset=DECK', `5, ${b[2]}, ${b[3]}, ${b[4]}, ${b[5]}`,
    '*Element, type=SPRINGA, elset=SPR', `6, ${b[0]}, ${b[1]}`,
    '*Nset, nset=BASE, generate', '1, 4, 1', '*Elset, elset=STEEL', 'LEG, BRACE,', '5',
    '*Beam Section, elset=LEG, material=S355, section=PIPE', '0.305, 0.019', '0., 1., 0.',
    '*Beam Section, elset=BRACE, material=S355, section=BOX', '0.2, 0.3, 0.01, 0.012, 0.01, 0.012',
    '*Shell Section, elset=DECK, material=S355', '0.025, 5', '*Solid Section, elset=BLOCK, material=GROUT', ',',
    '*Material, name=S355', '*Elastic', '2.1e11, 0.3', '*Density', '7850.,', '*Step', '*Static', '*Boundary', 'BASE, 1, 6', '*End Step',
  ].join('\r\n');
  const r = parseAbaqus({ name: 'deck.inp', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { hex20: 1, tet10: 1, line3: 1, line2: 1, quad4: 1 });
  checkElements(byKind(r, 'volume'));
  const f = byKind(r, 'frame');
  // the B32 element is split at its mid node, which Abaqus lists second
  assert.equal(f.members.length, 3);
  closeArr([f.nodes[3 * f.members[0].a + 2], f.nodes[3 * f.members[0].b + 2], f.nodes[3 * f.members[1].b + 2]], [5, 7, 9]);
  assert.deepEqual(f.members[0].section, { type: 'tube', D: 0.61, t: 0.019 });
  assert.equal(f.members[0].material, 'S355');
  assert.deepEqual(f.members[0].orient, [0, 1, 0]);
  assert.equal(f.members[2].section.shape, 'box');
  // thin-walled box 0.2 × 0.3, t = 10 / 12 mm: A = bh − (b − 2tw)(h − 2tf)
  close(f.members[2].section.A, 0.2 * 0.3 - 0.18 * 0.276, 1e-12);
  assert.deepEqual([...f.sets.members.STEEL], [0, 2]);
  const s = byKind(r, 'surface');
  close(s.groups[0].thickness, 0.025, 1e-12);
  assert.deepEqual([...byKind(r, 'volume').sets.nodes.BASE].length, 4);
  assert.deepEqual(r.meta.materials.S355, { units: 'as written in the deck (not converted)', E: 2.1e11, nu: 0.3, rho: 7850 });
  assert.equal(r.meta.skipped.keywords['element type SPRINGA'], 1);
  assert.equal(r.meta.heading, 'test deck');
  assert.deepEqual(abaqusType('CPE4R'), ['quad4', 4]);
  assert.equal(abaqusType('MASS'), null);
});

test('Abaqus: parts instanced twice with translation and rotation', () => {
  const text = ['*Part, name=Pile', '*Node', '1, 0, 0, 0', '2, 0, 0, -10', '*Element, type=B31, elset=ALL', '1, 1, 2', '*Beam Section, elset=ALL, material=M, section=PIPE', '1.0, 0.05', '*End Part',
    '*Assembly, name=A', '*Instance, name=P1, part=Pile', '*End Instance', '*Instance, name=P2, part=Pile', '20., 0., 0.', '20., 0., 0., 20., 1., 0., 90.', '*End Instance',
    '*Nset, nset=HEADS, instance=P2', '1,', '*End Assembly'].join('\n');
  const r = parseAbaqus({ name: 'piles.inp', text });
  const f = r.models[0];
  assert.equal(f.members.length, 2);
  // P2: translated to x = 20, then rotated 90° about the y-parallel axis through (20, 0, 0): −z → −x
  closeArr([...f.nodes], [0, 0, 0, 0, 0, -10, 20, 0, 0, 10, 0, 0], 1e-9);
  assert.deepEqual(f.members.map((m) => m.group), ['P1.PIPE ALL', 'P2.PIPE ALL']);
  assert.deepEqual([...f.sets.nodes.HEADS], [2]);
  assert.deepEqual(r.meta.instances, ['P1', 'P2']);
  assert.throws(() => parseAbaqus({ name: 'x.inp', text: 'just text\n' }), /no Abaqus keyword lines/);
  assert.throws(() => parseAbaqus({ name: 'x.inp', text: '*Heading\n x\n*Step\n*End Step\n' }), /no nodes were found \(expected a \*NODE block\)/);
});

// ---- ANSYS CDB -------------------------------------------------------------------------------------------

test('ANSYS CDB: NBLOCK/EBLOCK, degenerate solids, sections, CMBLOCK ranges', () => {
  const T = nodeTable();
  const hex = T.add(REF.hex);
  const tet = T.add(REF.tet.map((p) => [p[0] + 2, p[1], p[2]]));
  const b = T.add([[0, 0, 5], [0, 0, 9], [4, 0, 9], [4, 3, 9]]);
  const e = (v) => v.toExponential(13).replace(/e([+-])(\d+)$/, (m, s, d) => `E${s}${d.padStart(3, '0')}`).padStart(21);
  const i9 = (...v) => v.map((x) => String(x).padStart(9)).join('');
  const nb = T.pts.map((p, i) => i9(i + 1, 0, 0) + p.filter((c, k) => k < 2 || c !== 0).map(e).join(''));
  const el = (mat, type, sec, id, nodes) => i9(mat, type, 1, sec, 0, 0, 0, 0, nodes.length, 0, id, ...nodes);
  const text = ['/PREP7', 'ET,1,185', 'ET,2,SHELL181', 'ET,3,188', 'ET,4,CONTA174',
    `NBLOCK,6,SOLID,${T.pts.length},${T.pts.length}`, '(3i9,6e21.13e3)', ...nb, 'N,R5.3,LOC,-1,',
    'EBLOCK,19,SOLID,5,5', '(19i9)',
    el(1, 1, 0, 1, hex), el(1, 1, 0, 2, [tet[0], tet[1], tet[2], tet[2], tet[3], tet[3], tet[3], tet[3]]),
    el(2, 2, 7, 3, [b[1], b[2], b[3], b[3]]), el(2, 3, 5, 4, [b[0], b[1], b[3]]), el(1, 4, 0, 5, hex.slice(0, 4)), '-1',
    'CMBLOCK,BASE,NODE,2', '(8i10)', `${'1'.padStart(10)}${'-4'.padStart(10)}`, 'CMBLOCK,SOLIDS,ELEM,2', '(8i10)', `${'1'.padStart(10)}${'2'.padStart(10)}`,
    'SECTYPE,5,PIPE,,riser', 'SECDATA,0.61,0.019', 'SECTYPE,7,SHELL', 'SECDATA,0.02,2', 'MP,EX,2,2.1E11', 'MPDATA,DENS,1,,7850.', 'D,1,UX,0', 'FINISH'].join('\n');
  const r = parseAnsysCdb({ name: 'model.cdb', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { hex8: 1, tet4: 1, tri3: 1, line2: 1 });
  checkElements(byKind(r, 'volume'));
  const f = byKind(r, 'frame');
  assert.deepEqual(f.members[0].section, { type: 'tube', D: 0.61, t: 0.019 });
  assert.equal(f.members[0].group, 'riser');
  assert.ok(f.members[0].orientNode > 0, 'BEAM188 third node is the orientation node, not a mid node');
  closeArr([...f.nodes], [0, 0, 5, 0, 0, 9]);
  close(byKind(r, 'surface').groups[0].thickness, 0.02, 1e-12);
  assert.deepEqual([...byKind(r, 'volume').sets.nodes.BASE], [0, 1, 2, 3]);
  assert.deepEqual([...byKind(r, 'volume').sets.cells.SOLIDS].length, 2);
  assert.equal(r.meta.skipped.records['element routine 174'], 1);
  assert.equal(r.meta.materials['material 2'].E, 2.1e11);
  assert.throws(() => parseAnsysCdb({ name: 'x.cdb', text: '/PREP7\nET,1,185\n' }), /no NBLOCK or EBLOCK was found/);
  assert.throws(() => parseAnsysCdb({ name: 'x.cdb', text: 'NBLOCK,6,SOLID\nno format\n' }), /not followed by a Fortran format line/);
});

// ---- UNV -------------------------------------------------------------------------------------------------

test('UNV: datasets 164, 2411, 2412, 2467 with perimeter-ordered parabolic elements', () => {
  const T = nodeTable();
  // I-DEAS orders parabolic elements around the perimeter: corner, mid, corner, mid …
  const tetC = REF.tet.map((p) => [p[0] * 1000, p[1] * 1000, p[2] * 1000]);
  const mid = (a, b) => a.map((c, k) => (c + b[k]) / 2);
  const tet = T.add([tetC[0], mid(tetC[0], tetC[1]), tetC[1], mid(tetC[1], tetC[2]), tetC[2], mid(tetC[2], tetC[0]), mid(tetC[0], tetC[3]), mid(tetC[1], tetC[3]), mid(tetC[2], tetC[3]), tetC[3]]);
  const q = REF.quad.map((p) => [p[0] * 1000 + 3000, p[1] * 1000, 0]);
  const quad = T.add([q[0], mid(q[0], q[1]), q[1], mid(q[1], q[2]), q[2], mid(q[2], q[3]), q[3], mid(q[3], q[0])]);
  const beam = T.add([[0, 0, 5000], [0, 0, 9000]]);
  const I = (...v) => v.map((x) => String(x).padStart(10)).join('');
  const D = (v) => v.toExponential(16).replace('e', 'D').padStart(25);
  const text = ['    -1', '   164', `         5mm (milli-newton)             2`, '  1.0000000000000000D+03  1.0000000000000000D+03  1.0000000000000000D+00', '  2.7314999999999998D+02', '    -1',
    '    -1', '  2411', ...T.pts.flatMap((p, i) => [I(i + 1, 1, 1, 11), p.map(D).join('')]), '    -1',
    '    -1', '  2412', I(1, 118, 2, 1, 7, 10), I(...tet.slice(0, 8)), I(...tet.slice(8)), I(2, 95, 3, 1, 7, 8), I(...quad), I(3, 21, 4, 1, 7, 2), I(0, 0, 0), I(...beam), I(4, 136, 1, 1, 7, 1), I(beam[0]), '    -1',
    '    -1', '  2467', I(1, 0, 0, 0, 0, 0, 0, 3), 'FIXED', I(7, tet[0], 0, 0, 7, tet[2], 0, 0), I(8, 1, 0, 0), '    -1',
    '    -1', '    55', '         1         2         3', '    -1', ''].join('\n');
  const r = parseUnv({ name: 'part.unv', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { tet10: 1, quad8: 1, line2: 1 });
  assert.equal(r.meta.units.unit, 'mm');
  assert.equal(r.meta.units.basis, 'file');
  checkElements(byKind(r, 'volume'));
  checkElements(byKind(r, 'surface'));
  close(r.meta.bbox.max[2], 9, 1e-9);
  assert.deepEqual([...byKind(r, 'volume').sets.nodes.FIXED].length, 2);
  assert.deepEqual([...byKind(r, 'volume').sets.cells.FIXED], [0]);
  assert.equal(r.meta.skipped['dataset 55'], 1);
  assert.equal(r.meta.skipped['element descriptor 136'], 1);
  assert.throws(() => parseUnv({ name: 'x.unv', text: 'nothing\n' }), /no universal-file datasets were found/);
});

// ---- VTK -------------------------------------------------------------------------------------------------

function vtkMesh() {
  const T = nodeTable();
  const hex = T.add(REF.hex).map((i) => i - 1);
  // VTK wedge base normal points away from the top face: canonical order with 1↔2 and 4↔5 swapped
  const w = REF.wedge.map((p) => [p[0] + 2, p[1], p[2]]);
  const wedge = T.add([w[0], w[2], w[1], w[3], w[5], w[4]]).map((i) => i - 1);
  const tri = T.add([[5, 0, 0], [6, 0, 0], [5, 1, 0], [6, 1, 0]]).map((i) => i - 1);
  const cells = [[12, hex], [13, wedge], [5, tri.slice(0, 3)], [8, tri], [3, [tri[0], tri[1]]]];
  return { pts: T.pts, cells, scalar: T.pts.map((p, i) => i * 0.5) };
}

test('VTK legacy ASCII: unstructured grid, wedge/pixel reordering, nodal scalars', () => {
  const { pts, cells, scalar } = vtkMesh();
  const size = cells.reduce((s, c) => s + 1 + c[1].length, 0);
  const text = `# vtk DataFile Version 3.0\ntest mesh\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ${pts.length} double\n${pts.map((p) => p.join(' ')).join('\n')}\nCELLS ${cells.length} ${size}\n${cells.map((c) => `${c[1].length} ${c[1].join(' ')}`).join('\n')}\nCELL_TYPES ${cells.length}\n${cells.map((c) => c[0]).join('\n')}\nPOINT_DATA ${pts.length}\nSCALARS stress double 1\nLOOKUP_TABLE default\n${scalar.join(' ')}\nVECTORS disp float\n${pts.map(() => '0 0 1').join('\n')}\nCELL_DATA ${cells.length}\nSCALARS mat int 1\nLOOKUP_TABLE default\n${cells.map(() => 1).join(' ')}\n`;
  const r = parseVtk({ name: 'mesh.vtk', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'volume']);
  assert.deepEqual(r.meta.elementTypes, { hex8: 1, wedge6: 1, tri3: 1, quad4: 1, line2: 1 });
  const vol = byKind(r, 'volume');
  checkElements(vol);
  // node data follow the nodes through compaction
  for (let i = 0; i < vol.nodeIds.length; i++) close(vol.nodeData.stress[i], scalar[vol.nodeIds[i]], 1e-6);
  assert.equal(vol.nodeDataComponents.disp, 3);
  const surf = byKind(r, 'surface');
  assert.equal(surf.indices.length, 9);
  close(surf.nodeData.stress[3], scalar[surf.nodeIds[3]], 1e-6);
  assert.deepEqual(r.meta.pointData, ['stress', 'disp']);
  assert.ok(r.warnings.some((w) => /cell-data array/.test(w)));
  assert.equal(r.meta.title, 'test mesh');
});

test('VTK legacy BINARY polydata (big-endian) and 5.1 OFFSETS/CONNECTIVITY layout', () => {
  const pts = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]];
  const be = (arr, set, sz) => { const b = new Uint8Array(arr.length * sz); const dv = new DataView(b.buffer); arr.forEach((v, i) => dv[set](i * sz, v, false)); return b; };
  const bin = concat(enc(`# vtk DataFile Version 4.2\nbinary\nBINARY\nDATASET POLYDATA\nPOINTS 5 float\n`), be(pts.flat(), 'setFloat32', 4), enc(`\nPOLYGONS 2 9\n`), be([4, 0, 1, 2, 3, 3, 1, 4, 2], 'setInt32', 4), enc(`\nPOINT_DATA 5\nSCALARS h float\nLOOKUP_TABLE default\n`), be([1, 2, 3, 4, 5], 'setFloat32', 4), enc('\n'));
  const r = parseVtk({ name: 'poly.vtk', bytes: bin });
  assertContract(r, ['surface']);
  assert.deepEqual(r.meta.elementTypes, { quad4: 1, tri3: 1 });
  assert.deepEqual([...r.models[0].nodeData.h], [1, 2, 3, 4, 5]);
  assert.equal(r.meta.encoding, 'binary');
  const v51 = `# vtk DataFile Version 5.1\nv\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 4 float\n0 0 0 1 0 0 0 1 0 0 0 1\nCELLS 2 4\nOFFSETS vtktypeint64\n0 4\nCONNECTIVITY vtktypeint64\n0 1 2 3\nCELL_TYPES 1\n10\n`;
  const t = parseVtk({ name: 't.vtk', text: v51 });
  assert.deepEqual(t.meta.elementTypes, { tet4: 1 });
  checkElements(t.models[0]);
});

test('VTK legacy: structured points image → grid; structured grid → quads; bad headers', () => {
  const img = `# vtk DataFile Version 3.0\ndem\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS 3 2 1\nORIGIN 100 200 0\nSPACING 10 20 1\nPOINT_DATA 6\nSCALARS elevation float 1\nLOOKUP_TABLE default\n-1 -2 -3\n-4 -5 -6\n`;
  const g = parseVtk({ name: 'dem.vtk', text: img });
  assertContract(g, ['grid']);
  assert.deepEqual([g.models[0].nx, g.models[0].ny, g.models[0].x0, g.models[0].y0, g.models[0].dx, g.models[0].dy], [3, 2, 100, 200, 10, 20]);
  assert.deepEqual([...g.models[0].z], [-1, -2, -3, -4, -5, -6]);
  const sg = `# vtk DataFile Version 3.0\nsg\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS 3 2 1\nPOINTS 6 float\n0 0 0 1 0 0 2 0 0 0 1 0 1 1 0 2 1 0\n`;
  const s = parseVtk({ name: 'sg.vtk', text: sg });
  assert.deepEqual(s.meta.elementTypes, { quad4: 2 });
  assert.throws(() => parseVtk({ name: 'x.vtk', text: 'hello\n' }), /expected the first line "# vtk DataFile/);
  assert.throws(() => parseVtk({ name: 'x.vtk', text: '# vtk DataFile Version 3.0\nt\nUNICODE\n' }), /expected ASCII or BINARY/);
  assert.throws(() => parseVtk({ name: 'x.vtk', text: '# vtk DataFile Version 3.0\nt\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 3 float\n0 0 0 1\n' }), /file ends while reading 9 values/);
});

const le = (arr, T) => new Uint8Array(new T(arr).buffer);
const b64 = (u8) => Buffer.from(u8).toString('base64');
const withHeader = (u8) => concat(le([u8.length], Uint32Array), u8);

function vtuParts() {
  const { pts, cells, scalar } = vtkMesh();
  const conn = cells.flatMap((c) => c[1]);
  let o = 0;
  const offs = cells.map((c) => (o += c[1].length));
  return { pts, cells, scalar, conn, offs, types: cells.map((c) => c[0]) };
}

test('VTU ascii and inline base64', async () => {
  const m = vtuParts();
  const mk = (fmt, f64, i32, u8) => `<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n<UnstructuredGrid><Piece NumberOfPoints="${m.pts.length}" NumberOfCells="${m.cells.length}">\n<PointData Scalars="stress"><DataArray type="Float64" Name="stress" format="${fmt}">${f64(m.scalar)}</DataArray></PointData>\n<CellData><DataArray type="Int32" Name="mat" format="${fmt}">${i32(m.types)}</DataArray></CellData>\n<Points><DataArray type="Float64" NumberOfComponents="3" format="${fmt}">${f64(m.pts.flat())}</DataArray></Points>\n<Cells><DataArray type="Int32" Name="connectivity" format="${fmt}">${i32(m.conn)}</DataArray><DataArray type="Int32" Name="offsets" format="${fmt}">${i32(m.offs)}</DataArray><DataArray type="UInt8" Name="types" format="${fmt}">${u8(m.types)}</DataArray></Cells>\n</Piece></UnstructuredGrid></VTKFile>`;
  const variants = {
    ascii: mk('ascii', (a) => a.join(' '), (a) => a.join(' '), (a) => a.join(' ')),
    base64: mk('binary', (a) => b64(withHeader(le(a, Float64Array))), (a) => b64(withHeader(le(a, Int32Array))), (a) => b64(withHeader(le(a, Uint8Array)))),
    // header encoded on its own, as some writers do
    'base64 split header': mk('binary', (a) => b64(le([a.length * 8], Uint32Array)) + b64(le(a, Float64Array)), (a) => b64(le([a.length * 4], Uint32Array)) + b64(le(a, Int32Array)), (a) => b64(le([a.length], Uint32Array)) + b64(le(a, Uint8Array))),
  };
  for (const [label, text] of Object.entries(variants)) {
    const r = await parseVtkXml({ name: 'mesh.vtu', bytes: enc(text) });
    assertContract(r, ['frame', 'surface', 'volume']);
    assert.deepEqual(r.meta.elementTypes, { hex8: 1, wedge6: 1, tri3: 1, quad4: 1, line2: 1 }, label);
    const vol = byKind(r, 'volume');
    checkElements(vol);
    for (let i = 0; i < vol.nodeIds.length; i++) close(vol.nodeData.stress[i], m.scalar[vol.nodeIds[i]], 1e-6, label);
    assert.equal(r.meta.skipped.cellDataArrays, 1);
  }
});

test('VTU appended raw with zlib compression (UInt64 headers) and VTP polydata', async () => {
  const m = vtuParts();
  const blocks = [];
  let off = 0;
  const add = (u8) => {
    const c = new Uint8Array(deflateSync(Buffer.from(u8)));
    const head = new Uint8Array(new BigUint64Array([1n, BigInt(u8.length), BigInt(u8.length), BigInt(c.length)]).buffer);
    const start = off;
    blocks.push(head, c);
    off += head.length + c.length;
    return start;
  };
  const oP = add(le(m.pts.flat(), Float32Array)), oC = add(new Uint8Array(new BigInt64Array(m.conn.map(BigInt)).buffer)), oO = add(new Uint8Array(new BigInt64Array(m.offs.map(BigInt)).buffer)), oT = add(le(m.types, Uint8Array));
  const xml = `<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64" compressor="vtkZLibDataCompressor">\n<UnstructuredGrid><Piece NumberOfPoints="${m.pts.length}" NumberOfCells="${m.cells.length}">\n<Points><DataArray type="Float32" NumberOfComponents="3" format="appended" offset="${oP}"/></Points>\n<Cells><DataArray type="Int64" Name="connectivity" format="appended" offset="${oC}"/><DataArray type="Int64" Name="offsets" format="appended" offset="${oO}"/><DataArray type="UInt8" Name="types" format="appended" offset="${oT}"/></Cells>\n</Piece></UnstructuredGrid>\n<AppendedData encoding="raw">\n_`;
  const r = await parseVtkXml({ name: 'mesh.vtu', bytes: concat(enc(xml), ...blocks, enc('\n</AppendedData>\n</VTKFile>\n')) });
  assert.deepEqual(r.meta.elementTypes, { hex8: 1, wedge6: 1, tri3: 1, quad4: 1, line2: 1 });
  checkElements(byKind(r, 'volume'));
  assert.equal(r.meta.appended, 'raw');
  const vtp = `<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian"><PolyData><Piece NumberOfPoints="5" NumberOfPolys="2" NumberOfLines="1">\n<Points><DataArray type="Float32" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 1 1 0 0 1 0 2 0 0</DataArray></Points>\n<Polys><DataArray type="Int32" Name="connectivity" format="ascii">0 1 2 3 1 4 2</DataArray><DataArray type="Int32" Name="offsets" format="ascii">4 7</DataArray></Polys>\n<Lines><DataArray type="Int32" Name="connectivity" format="ascii">0 1 4</DataArray><DataArray type="Int32" Name="offsets" format="ascii">3</DataArray></Lines>\n</Piece></PolyData></VTKFile>`;
  const p = await parseVtkXml({ name: 'poly.vtp', text: vtp });
  assert.deepEqual(p.meta.elementTypes, { quad4: 1, tri3: 1, line2: 2 });
  await assert.rejects(parseVtkXml({ name: 'x.vtu', text: '<VTKFile type="ImageData"></VTKFile>' }), /type "ImageData" is not read/);
  await assert.rejects(parseVtkXml({ name: 'x.vtu', text: '<svg/>' }), /expected root element <VTKFile>/);
});

// ---- SU2 and OpenFOAM ------------------------------------------------------------------------------------

test('SU2: 3-D mesh with boundary markers as surface groups; 2-D mesh with marker edges', () => {
  const text = `% tet\nNDIME= 3\nNELEM= 1\n10 0 1 2 3 0\nNPOIN= 4\n0 0 0 0\n1 0 0 1\n0 1 0 2\n0 0 1 3\nNMARK= 2\nMARKER_TAG= wall\nMARKER_ELEMS= 3\n5 0 2 1\n5 0 1 3\n5 0 3 2\nMARKER_TAG= farfield\nMARKER_ELEMS= 1\n5 1 2 3\n`;
  const r = parseSu2({ name: 'tet.su2', bytes: enc(text) });
  assertContract(r, ['surface', 'volume']);
  checkElements(byKind(r, 'volume'));
  const s = byKind(r, 'surface');
  assert.deepEqual(s.groups.map((g) => [g.name, g.count]), [['wall', 9], ['farfield', 3]]);
  assert.deepEqual(r.meta.markers, [{ name: 'wall', elements: 3 }, { name: 'farfield', elements: 1 }]);
  const two = `NDIME= 2\nNELEM= 2\n5 0 1 2\n5 0 2 3\nNPOIN= 4\n0 0\n1 0\n1 1\n0 1\nNMARK= 1\nMARKER_TAG= airfoil\nMARKER_ELEMS= 2\n3 0 1\n3 1 2\n`;
  const t = parseSu2({ name: 'sq.su2', text: two });
  assert.deepEqual(t.meta.elementTypes, { tri3: 2, line2: 2 });
  assert.equal(byKind(t, 'frame').name, 'sq (boundary edges)');
  assert.equal(byKind(t, 'frame').members[0].group, 'airfoil');
  assert.throws(() => parseSu2({ name: 'x.su2', text: 'hello\n' }), /expected NDIME=, NELEM= and NPOIN=/);
});

test('OpenFOAM polyMesh: boundary patches as groups, cell count reported', async () => {
  const head = (cls, obj) => `/*--------------------------------*- C++ -*----------------------------------*\\\n| OpenFOAM |\n\\*---------------------------------------------------------------------------*/\nFoamFile\n{\n    version 2.0;\n    format ascii;\n    class ${cls};\n    object ${obj};\n}\n// * * * //\n`;
  const files = new Map([
    ['constant/polyMesh/points', enc(head('vectorField', 'points') + `8\n(\n${REF.hex.map((p) => `(${p.join(' ')})`).join('\n')}\n)\n`)],
    ['constant/polyMesh/faces', enc(head('faceList', 'faces') + '6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(1 2 6 5)\n4(2 3 7 6)\n4(3 0 4 7)\n)\n')],
    ['constant/polyMesh/owner', enc(head('labelList', 'owner') + '6\n(\n0\n0\n0\n0\n0\n0\n)\n')],
    ['constant/polyMesh/boundary', enc(head('polyBoundaryMesh', 'boundary') + '2\n(\n    seabed\n    {\n        type wall;\n        nFaces 1;\n        startFace 0;\n    }\n    sides\n    {\n        type patch;\n        nFaces 5;\n        startFace 1;\n    }\n)\n')],
  ]);
  const r = await parseOpenFoam({ name: 'case', files });
  assertContract(r, ['surface']);
  assert.deepEqual(r.models[0].groups.map((g) => [g.name, g.count, g.type]), [['seabed', 6, 'wall'], ['sides', 30, 'patch']]);
  assert.equal(r.meta.cells, 1);
  assert.equal(r.meta.elementTypes.quad4, 6);
  assert.ok(r.warnings.some((w) => /not reconstructed/.test(w)));
  files.delete('constant/polyMesh/faces');
  await assert.rejects(parseOpenFoam({ name: 'case', files }), /missing polyMesh file\(s\): faces/);
  files.set('constant/polyMesh/faces', enc(head('faceList', 'faces').replace('ascii', 'binary') + '6\n('));
  await assert.rejects(parseOpenFoam({ name: 'case', files }), /binary OpenFOAM files are not read/);
});

test('truncated mesh files fail with a worded error, never an internal one', async () => {
  const pts = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]];
  const be = (arr, set, sz) => { const b = new Uint8Array(arr.length * sz); const dv = new DataView(b.buffer); arr.forEach((v, i) => dv[set](i * sz, v, false)); return b; };
  const bin = concat(enc('# vtk DataFile Version 4.2\nbinary\nBINARY\nDATASET POLYDATA\nPOINTS 5 float\n'), be(pts.flat(), 'setFloat32', 4), enc('\nPOLYGONS 2 9\n'), be([4, 0, 1, 2, 3, 3, 1, 4, 2], 'setInt32', 4), enc('\nPOINT_DATA 5\nSCALARS h float\nLOOKUP_TABLE default\n'), be([1, 2, 3, 4, 5], 'setFloat32', 4), enc('\n'));
  await assertTruncationSafe(bin, (b) => parseVtk({ name: 'p.vtk', bytes: b }), 'VTK binary');
  const m = vtuParts();
  const blocks = [];
  let off = 0;
  const add = (u8) => { const c = new Uint8Array(deflateSync(Buffer.from(u8))); const head = new Uint8Array(new Uint32Array([1, u8.length, u8.length, c.length]).buffer); const start = off; blocks.push(head, c); off += head.length + c.length; return start; };
  const oP = add(le(m.pts.flat(), Float32Array)), oC = add(le(m.conn, Int32Array)), oO = add(le(m.offs, Int32Array)), oT = add(le(m.types, Uint8Array));
  const xml = `<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor"><UnstructuredGrid><Piece NumberOfPoints="${m.pts.length}" NumberOfCells="${m.cells.length}"><Points><DataArray type="Float32" NumberOfComponents="3" format="appended" offset="${oP}"/></Points><Cells><DataArray type="Int32" Name="connectivity" format="appended" offset="${oC}"/><DataArray type="Int32" Name="offsets" format="appended" offset="${oO}"/><DataArray type="UInt8" Name="types" format="appended" offset="${oT}"/></Cells></Piece></UnstructuredGrid><AppendedData encoding="raw">_`;
  const vtu = concat(enc(xml), ...blocks, enc('</AppendedData></VTKFile>'));
  const full = await parseVtkXml({ name: 'm.vtu', bytes: vtu });
  assert.deepEqual(full.meta.elementTypes, { hex8: 1, wedge6: 1, tri3: 1, quad4: 1, line2: 1 });
  await assertTruncationSafe(vtu, (b) => parseVtkXml({ name: 'm.vtu', bytes: b }), 'VTU appended raw + zlib');
  const b64 = (u8) => Buffer.from(u8).toString('base64');
  const inline = `<VTKFile type="UnstructuredGrid" byte_order="LittleEndian"><UnstructuredGrid><Piece><Points><DataArray type="Float64" NumberOfComponents="3" format="binary">${b64(withHeader(le(m.pts.flat(), Float64Array)))}</DataArray></Points><Cells><DataArray type="Int32" Name="connectivity" format="binary">${b64(withHeader(le(m.conn, Int32Array)))}</DataArray><DataArray type="Int32" Name="offsets" format="binary">${b64(withHeader(le(m.offs, Int32Array)))}</DataArray><DataArray type="UInt8" Name="types" format="binary">${b64(withHeader(le(m.types, Uint8Array)))}</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>`;
  await assertTruncationSafe(enc(inline), (b) => parseVtkXml({ name: 'm.vtu', bytes: b }), 'VTU inline base64');
});

test('truncated solver decks fail with a worded error, never an internal one', async () => {
  const msh2 = '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n3 7 "steel"\n$EndPhysicalNames\n$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n$EndNodes\n$Elements\n2\n1 4 2 7 1 1 2 3 4\n2 2 2 3 9 1 2 3\n$EndElements\n';
  await assertTruncationSafe(enc(msh2), (b) => parseGmsh({ name: 'a.msh', bytes: b }), 'MSH 2.2');
  const msh4 = '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Entities\n0 0 1 1\n5 0 0 0 1 1 0 1 1 0\n9 0 0 0 1 1 1 1 2 0\n$EndEntities\n$Nodes\n1 4 1 4\n3 9 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n$Elements\n2 2 1 2\n2 5 2 1\n1 1 3 2\n3 9 4 1\n2 1 2 3 4\n$EndElements\n';
  await assertTruncationSafe(enc(msh4), (b) => parseGmsh({ name: 'a.msh', bytes: b }), 'MSH 4.1');
  const bdf = ['BEGIN BULK', 'CORD2C,5,0,10.,0.,0.,10.,0.,1.,+', '+,11.,0.,0.', 'GRID,1,5,2.,90.,3.', 'GRID,2,,0.,0.,0.', 'GRID,3,,1.,0.,0.', 'GRID,4,,0.,1.,0.', 'GRID*                  5                             0.0             0.0*', '*                    1.0', sf('CBAR', 1, 10, 1, 2, '0.0', '1.0', '0.0'), sf('CTETRA', 2, 30, 2, 3, 4, 5), sf('CQUAD4', 3, 20, 1, 2, 3, 4), sf('PBARL', 10, 1, '', 'TUBE'), sf('', '0.3', '0.28'), sf('PSHELL', 20, 1, '0.02'), sf('MAT1', 1, '2.1+11', '', '0.3', '7850.'), 'ENDDATA'].join('\n');
  await assertTruncationSafe(enc(bdf), (b) => parseNastran({ name: 'a.bdf', bytes: b }), 'BDF');
  const inp = ['*Heading', ' t', '*Part, name=P', '*Node', '1, 0, 0, 0', '2, 1, 0, 0', '3, 0, 1, 0', '4, 0, 0, 1', '*Element, type=C3D4, elset=E', '1, 1, 2, 3, 4', '*Element, type=B31, elset=B', '2, 1, 2', '*Nset, nset=N, generate', '1, 4, 1', '*Elset, elset=ALL', 'E, B', '*Beam Section, elset=B, material=M, section=PIPE', '0.3, 0.02', '0., 1., 0.', '*Solid Section, elset=E, material=M', ',', '*End Part', '*Assembly, name=A', '*Instance, name=I, part=P', '1., 0., 0.', '1., 0., 0., 1., 1., 0., 90.', '*End Instance', '*End Assembly', '*Material, name=M', '*Elastic', '2.1e11, 0.3'].join('\n');
  await assertTruncationSafe(enc(inp), (b) => parseAbaqus({ name: 'a.inp', bytes: b }), 'INP');
  const e21 = (v) => v.toExponential(13).replace(/e([+-])(\d+)$/, (m, s, d) => `E${s}${d.padStart(3, '0')}`).padStart(21);
  const i9 = (...v) => v.map((x) => String(x).padStart(9)).join('');
  const cdb = ['ET,1,185', 'ET,2,188', 'NBLOCK,6,SOLID,4,4', '(3i9,6e21.13e3)', ...REF.tet.map((p, i) => i9(i + 1, 0, 0) + p.map(e21).join('')), 'N,R5.3,LOC,-1,', 'EBLOCK,19,SOLID,2,2', '(19i9)', i9(1, 1, 1, 0, 0, 0, 0, 0, 8, 0, 1, 1, 2, 3, 3, 4, 4, 4, 4), i9(1, 2, 1, 5, 0, 0, 0, 0, 3, 0, 2, 1, 2, 3), '-1', 'CMBLOCK,BASE,NODE,2', '(8i10)', `${'1'.padStart(10)}${'-3'.padStart(10)}`, 'SECTYPE,5,PIPE,,r', 'SECDATA,0.6,0.02', 'MP,EX,1,2.1E11'].join('\n');
  await assertTruncationSafe(enc(cdb), (b) => parseAnsysCdb({ name: 'a.cdb', bytes: b }), 'CDB');
  const I = (...v) => v.map((x) => String(x).padStart(10)).join('');
  const D = (v) => v.toExponential(16).replace('e', 'D').padStart(25);
  const unv = ['    -1', '   164', '         1SI', '  1.0D+00  1.0D+00  1.0D+00', '  2.7D+02', '    -1', '    -1', '  2411', ...REF.tet.flatMap((p, i) => [I(i + 1, 1, 1, 11), p.map(D).join('')]), '    -1', '    -1', '  2412', I(1, 111, 2, 1, 7, 4), I(1, 2, 3, 4), I(2, 21, 4, 1, 7, 2), I(0, 0, 0), I(1, 2), '    -1', '    -1', '  2467', I(1, 0, 0, 0, 0, 0, 0, 2), 'G', I(7, 1, 0, 0, 8, 1, 0, 0), '    -1'].join('\n');
  await assertTruncationSafe(enc(unv), (b) => parseUnv({ name: 'a.unv', bytes: b }), 'UNV');
  const su2 = 'NDIME= 3\nNELEM= 1\n10 0 1 2 3 0\nNPOIN= 4\n0 0 0 0\n1 0 0 1\n0 1 0 2\n0 0 1 3\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 2\n5 0 2 1\n5 0 1 3\n';
  await assertTruncationSafe(enc(su2), (b) => parseSu2({ name: 'a.su2', bytes: b }), 'SU2');
  const head = 'FoamFile\n{\n format ascii;\n}\n';
  const pointsFile = enc(`${head}8\n(\n${REF.hex.map((p) => `(${p.join(' ')})`).join('\n')}\n)\n`);
  const facesFile = enc(`${head}6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(1 2 6 5)\n4(2 3 7 6)\n4(3 0 4 7)\n)\n`);
  const boundaryFile = enc(`${head}1\n(\n w\n {\n  type wall;\n  nFaces 6;\n  startFace 0;\n }\n)\n`);
  await assertTruncationSafe(facesFile, (b) => parseOpenFoam({ name: 'c', files: new Map([['points', pointsFile], ['faces', b], ['boundary', boundaryFile]]) }), 'polyMesh faces');
  await assertTruncationSafe(pointsFile, (b) => parseOpenFoam({ name: 'c', files: new Map([['points', b], ['faces', facesFile], ['boundary', boundaryFile]]) }), 'polyMesh points');
  await assertTruncationSafe(boundaryFile, (b) => parseOpenFoam({ name: 'c', files: new Map([['points', pointsFile], ['faces', facesFile], ['boundary', b]]) }), 'polyMesh boundary');
});
