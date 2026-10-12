import test from 'node:test';
import assert from 'node:assert/strict';
import { parseDxf, bulgeArc, ocsMatrix } from '../../src/geometry/importers/dxf.js';
import { parseSdnf } from '../../src/geometry/importers/sdnf.js';
import { parseStep, parseCis2, P21, Ref, Typed } from '../../src/geometry/importers/step.js';
import { parseIfc } from '../../src/geometry/importers/ifc.js';
import { parseIges, igesFields } from '../../src/geometry/importers/iges.js';
import { parseSvg, parseAirfoil, flattenPath } from '../../src/geometry/importers/profiles.js';
import { parseGdf, parseNemoh } from '../../src/geometry/importers/panels.js';
import { sampleCurve, tessellateSurface } from '../../src/geometry/importers/nurbs.js';
import { ringArea, M4 } from '../../src/geometry/importers/common.js';
import { resolveSection } from '../../src/geometry/importers/sections.js';
import { enc, close, closeArr, assertContract, surfaceArea, surfaceVolume, CUBE_V, CUBE_Q, CUBE_T, assertTruncationSafe } from './fixtures/helpers.js';

const byKind = (r, k) => r.models.find((m) => m.kind === k);

// ---- NURBS ----------------------------------------------------------------------------------------------

test('NURBS: rational quadratic quarter circle stays on the circle; bilinear patch area', () => {
  const w = Math.SQRT1_2;
  const pts = sampleCurve(2, [0, 0, 0, 1, 1, 1], [1, 0, 0, 1, 1, 0, 0, 1, 0], [1, w, 1], 16);
  for (let i = 0; i < pts.length; i += 3) close(Math.hypot(pts[i], pts[i + 1]), 1, 1e-12);
  closeArr([...pts.slice(0, 3), ...pts.slice(-3)], [1, 0, 0, 0, 1, 0], 1e-12);
  // without weights the same control polygon is a parabola: its mid-point is at (0.75, 0.75)
  const par = sampleCurve(2, [0, 0, 0, 1, 1, 1], [1, 0, 0, 1, 1, 0, 0, 1, 0], null, 2);
  closeArr([par[3], par[4]], [0.75, 0.75], 1e-12);
  const s = tessellateSurface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 0, 2, 0, 0, 0, 3, 0, 2, 3, 0], null, 1, 1, 4);
  close(surfaceArea({ positions: s.coords, indices: s.indices }), 6, 1e-12);
  assert.equal(sampleCurve(3, [0, 0, 1, 1], [0, 0, 0, 1, 0, 0], null), null, 'too few control points for the degree');
});

// ---- DXF ------------------------------------------------------------------------------------------------

const dxf = (...pairs) => { let s = ''; for (let i = 0; i < pairs.length; i += 2) s += `${String(pairs[i]).padStart(3)}\n${pairs[i + 1]}\n`; return s; };
const dxfFile = ({ header = [], tables = [], blocks = [], entities }) => dxf(0, 'SECTION', 2, 'HEADER', ...header, 0, 'ENDSEC', 0, 'SECTION', 2, 'TABLES', ...tables, 0, 'ENDSEC', 0, 'SECTION', 2, 'BLOCKS', ...blocks, 0, 'ENDSEC', 0, 'SECTION', 2, 'ENTITIES', ...entities, 0, 'ENDSEC', 0, 'EOF');
const LINE = (layer, a, b) => [0, 'LINE', 8, layer, 10, a[0], 20, a[1], 30, a[2], 11, b[0], 21, b[1], 31, b[2]];

test('DXF: 3-D LINE model in millimetres → welded frame with layers as groups', () => {
  const text = dxfFile({
    header: [9, '$ACADVER', 1, 'AC1027', 9, '$INSUNITS', 70, 4],
    tables: [0, 'TABLE', 2, 'LAYER', 0, 'LAYER', 2, 'LEGS', 70, 0, 0, 'LAYER', 2, 'BRACES', 70, 0, 0, 'ENDTAB'],
    entities: [...LINE('LEGS', [0, 0, 0], [0, 0, 10000]), ...LINE('LEGS', [8000, 0, 0], [8000, 0, 10000]), ...LINE('BRACES', [0, 0, 10000], [8000, 0, 10000]), ...LINE('BRACES', [0, 0, 0], [8000, 0, 10000]), 0, 'TEXT', 8, '0', 10, 0, 20, 0, 1, 'label', 0, '3DSOLID', 8, '0', 1, 'acis data'],
  }).replace(/\n/g, '\r\n');
  const r = parseDxf({ name: 'jacket.dxf', bytes: enc(text) });
  assertContract(r, ['frame']);
  const f = r.models[0];
  assert.equal(f.nodes.length / 3, 4);
  assert.equal(f.members.length, 4);
  assert.deepEqual(f.members.map((m) => m.group), ['LEGS', 'LEGS', 'BRACES', 'BRACES']);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [8, 0, 10] });
  assert.equal(r.meta.units.unit, 'mm');
  assert.equal(r.meta.units.basis, 'file');
  assert.deepEqual(r.meta.layers, ['LEGS', 'BRACES']);
  assert.deepEqual(r.meta.elementTypes, { LINE: 4 });
  assert.equal(r.meta.skipped.TEXT, 1);
  assert.ok(r.warnings.some((w) => /ACIS data.*Export them as STEP or IGES/.test(w)));
});

test('DXF: closed LWPOLYLINE with bulge and CIRCLE → section loops; bulge geometry', () => {
  const out = [0, 0];
  bulgeArc(0, 0, 0, 2, -1, out); // clockwise semicircle from (0,0) to (0,2): passes through (−1, 1)
  const mid = out.length / 4;
  closeArr([out[2 * Math.floor(mid)], out[2 * Math.floor(mid) + 1]], [-1, 1], 1e-9);
  const text = dxfFile({ entities: [
    0, 'LWPOLYLINE', 8, 'SHEATH', 90, 4, 70, 1, 10, 0, 20, 0, 10, 2, 20, 0, 42, 1, 10, 2, 20, 1, 10, 0, 20, 1,
    0, 'CIRCLE', 8, 'CORE', 10, 0.5, 20, 0.5, 30, 0, 40, 0.25,
    0, 'LWPOLYLINE', 8, 'OPEN', 90, 2, 70, 0, 10, 5, 20, 5, 10, 6, 20, 6,
  ] });
  const r = parseDxf({ name: 'section.dxf', text });
  assertContract(r, ['section', 'polyline']);
  const s = byKind(r, 'section');
  assert.equal(s.loops.length, 2);
  // 2 × 1 rectangle with a semicircular end of radius 0.5 (polygonal approximation)
  close(Math.abs(ringArea(s.loops[0])), 2 + (Math.PI * 0.25) / 2, 2e-3);
  close(Math.abs(ringArea(s.loops[1])), Math.PI * 0.0625, 3e-4);
  assert.deepEqual(s.layers, [{ name: 'SHEATH', loops: [0] }, { name: 'CORE', loops: [1] }]);
  assert.equal(byKind(r, 'polyline').layer, 'OPEN');
  assert.ok(r.warnings.some((w) => /declares no length unit/.test(w)));
});

test('DXF: separate LINE + ARC entities closing a flat outline are assembled into a section', () => {
  const ents = [...LINE('P', [0, 0, 0], [1, 0, 0]), 0, 'ARC', 8, 'P', 10, 1, 20, 0.5, 30, 0, 40, 0.5, 50, 270, 51, 90, ...LINE('P', [1, 1, 0], [0, 1, 0]), ...LINE('P', [0, 0, 0], [0, 1, 0])];
  const r = parseDxf({ name: 'd.dxf', text: dxfFile({ entities: ents }) });
  assertContract(r, ['section']);
  close(Math.abs(ringArea(r.models[0].loops[0])), 1 + Math.PI / 8, 1e-3);
  assert.ok(r.warnings.some((w) => /assembled from separate LINE/.test(w)));
  const f = parseDxf({ name: 'd.dxf', text: dxfFile({ entities: ents }) }, { dxfLines: 'frame' });
  assert.deepEqual(f.models.map((m) => m.kind).sort(), ['frame', 'polyline']);
  // an extra dangling line means the drawing is not purely outlines: lines stay members
  const g = parseDxf({ name: 'd.dxf', text: dxfFile({ entities: [...ents, ...LINE('P', [5, 5, 0], [6, 5, 0])] }) });
  assert.equal(byKind(g, 'frame').members.length, 4);
});

test('DXF: 3DFACE, polyface mesh, INSERT with scale and rotation, extrusion direction, SPLINE', () => {
  const face = (v, layer = 'SKIN') => [0, '3DFACE', 8, layer, ...v.flatMap((p, k) => [10 + k, p[0], 20 + k, p[1], 30 + k, p[2]])];
  const poly = [0, 'POLYLINE', 8, 'MESH', 66, 1, 70, 64, 71, 4, 72, 1,
    ...[[0, 0, 5], [1, 0, 5], [1, 1, 5], [0, 1, 5]].flatMap((p) => [0, 'VERTEX', 8, 'MESH', 10, p[0], 20, p[1], 30, p[2], 70, 192]),
    0, 'VERTEX', 8, 'MESH', 10, 0, 20, 0, 30, 0, 70, 128, 71, 1, 72, 2, 73, 3, 74, -4, 0, 'SEQEND'];
  const text = dxfFile({
    blocks: [0, 'BLOCK', 2, 'PAD', 10, 0, 20, 0, 30, 0, ...face([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], '0'), 0, 'ENDBLK'],
    entities: [
      ...face([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0]]), ...poly,
      0, 'INSERT', 8, 'PADS', 2, 'PAD', 10, 10, 20, 0, 30, 0, 41, 2, 42, 3, 43, 1, 50, 90,
      0, 'CIRCLE', 8, 'C', 10, 5, 20, 0, 30, 0, 40, 1, 210, 0, 220, 0, 230, -1,
      0, 'SPLINE', 8, 'S', 70, 4, 71, 2, 72, 6, 73, 3, 40, 0, 40, 0, 40, 0, 40, 1, 40, 1, 40, 1, 10, 1, 20, 0, 30, 9, 10, 1, 20, 1, 30, 9, 10, 0, 20, 1, 30, 9, 41, 1, 41, Math.SQRT1_2, 41, 1,
    ],
  });
  const r = parseDxf({ name: 'm.dxf', text }, { units: 'm' });
  assertContract(r, ['surface', 'section', 'polyline']);
  const s = byKind(r, 'surface');
  assert.deepEqual(s.groups.map((g) => [g.name, g.count / 3]), [['SKIN', 1], ['MESH', 2], ['PADS', 2]]);
  // the block's face is on layer 0, so it takes the layer of the INSERT
  // inserted unit square: scaled 2 × 3, rotated 90° about z, moved to x = 10 → spans x 7…10, y 0…2
  const g = s.groups[2];
  let lo = [9e9, 9e9], hi = [-9e9, -9e9];
  for (let i = g.start; i < g.start + g.count; i++) { const k = 3 * s.indices[i]; lo = [Math.min(lo[0], s.positions[k]), Math.min(lo[1], s.positions[k + 1])]; hi = [Math.max(hi[0], s.positions[k]), Math.max(hi[1], s.positions[k + 1])]; }
  closeArr([...lo, ...hi], [7, 0, 10, 2], 1e-6);
  close(surfaceArea(s), 0.5 + 1 + 6, 1e-6);
  // extrusion (0, 0, −1): the OCS x axis is world −x, so the circle centre (5, 0) lies at world x = −5
  const circ = byKind(r, 'section').loops[0];
  let cx = 0;
  for (let i = 0; i < circ.length; i += 2) cx += circ[i];
  close(cx / (circ.length / 2), -5, 1e-9);
  closeArr([...M4.apply(ocsMatrix(0, 0, -1), 1, 0, 0)], [-1, 0, 0], 1e-12);
  const sp = byKind(r, 'polyline');
  for (let i = 0; i < sp.points.length; i += 3) { close(Math.hypot(sp.points[i], sp.points[i + 1]), 1, 1e-9); close(sp.points[i + 2], 9, 1e-12); }
  assert.equal(r.meta.elementTypes.POLYFACE, 1);
  assert.equal(r.meta.blocks, 1);
});

test('DXF: binary and non-DXF input are rejected with the remedy', () => {
  assert.throws(() => parseDxf({ name: 'b.dxf', bytes: enc('AutoCAD Binary DXF\r\n\x1a\0binarybinary') }), /binary DXF is not read\. Save the drawing as ASCII DXF/);
  assert.throws(() => parseDxf({ name: 'x.dxf', text: 'this is\nnot a\ndrawing\nfile\n' }), /not an ASCII DXF drawing with entities/);
  assert.throws(() => parseDxf({ name: 't.dxf', text: dxfFile({ entities: [0, 'TEXT', 8, '0', 1, 'only text'] }) }), /contains no geometry this reader keeps \(found only: TEXT\)/);
});

// ---- SDNF -----------------------------------------------------------------------------------------------

test('SDNF 3.0: members with resolved sections, plates with thickness, other packets reported', () => {
  const text = `# SDNF export
Packet 00
"Engineering Co" "Client"
"Topsides module M1" "Project X"
Packet 10
"millimeters" 3
1 0 0 0 "BEAM" "B-101" 0
"W12X26" "A36" 0.0 0 0
0.0 0.0 1.0   0.0 0.0 5000.0   6000.0 0.0 5000.0   0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
0 0 0 0 0 0 0 0 0 0 1 1
2 0 0 0 "COLUMN" "C-1" 0
"PIPE 610x19" "S355" 90.0 0 0
1.0 0.0 0.0   0.0 0.0 0.0   0.0 0.0 5000.0   0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
0 0 0 0 0 0 0 0 0 0 0 0
3 0 0 0 "BRACE" "X-1" 0
"L4X4X1/2" "A36" 0.0 0 0
0.0 1.0 0.0   0.0 0.0 0.0   6000.0 0.0 5000.0   0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
0 0 0 0 0 0 0 0 0 0 0 0
Packet 20
"millimeters" "millimeters" 1
7 0 0 0 "PLATE" "PL-1" "S355" 25.0 4
0.0 0.0 5000.0
6000.0 0.0 5000.0
6000.0 3000.0 5000.0
0.0 3000.0 5000.0
Packet 30
1 "load"
`;
  const r = parseSdnf({ name: 'module.sdnf', bytes: enc(text) });
  assertContract(r, ['frame', 'surface']);
  const f = byKind(r, 'frame');
  assert.equal(f.nodes.length / 3, 3, 'shared ends are welded');
  assert.equal(f.members.length, 3);
  const [beam, col, brace] = f.members;
  assert.equal(beam.group, 'BEAM');
  assert.equal(beam.material, 'A36');
  assert.equal(beam.mark, 'B-101');
  assert.equal(beam.section.designation, 'W12X26');
  // W12X26: A = 26 lb/ft ÷ 3.403 = 7.64 in² (AISC table 7.65 in²)
  close(beam.section.A / 0.0254 ** 2, 7.65, 0.02);
  assert.deepEqual(beam.releases, { a: [0, 0, 0, 0, 0, 0], b: [0, 0, 0, 0, 1, 1] });
  assert.deepEqual(col.section, { type: 'tube', D: 0.61, t: 0.019, name: 'PIPE 610x19' });
  assert.equal(col.rotationDeg, 90);
  assert.equal(brace.sectionName, 'L4X4X1/2');
  closeArr([f.nodes[3 * beam.a], f.nodes[3 * beam.a + 2], f.nodes[3 * beam.b], f.nodes[3 * beam.b + 2]], [0, 5, 6, 5], 1e-12);
  const s = byKind(r, 'surface');
  close(surfaceArea(s), 18, 1e-9);
  close(s.groups[0].thickness, 0.025, 1e-12);
  assert.equal(s.groups[0].name, 'PL-1');
  assert.equal(r.meta.units.unit, 'mm');
  assert.deepEqual(r.meta.title.slice(0, 2), ['Engineering Co', 'Client']);
  assert.ok(r.meta.skipped['Packet 30 (member loads)']);
  assert.ok(r.warnings.some((w) => /L4X4X1\/2/.test(w)));
  assert.throws(() => parseSdnf({ name: 'x.sdnf', text: 'nothing\n' }), /expected lines "Packet 00"/);
  assert.equal(resolveSection('HE 300 B').section.designation, 'HE 300 B');
});

// ---- ISO 10303-21 ---------------------------------------------------------------------------------------

const p21 = (schema, data) => `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test'),'2;1');\nFILE_NAME('part.stp','2024-01-01T00:00:00',('a'),('o'),'pre','Test System 1.0','');\nFILE_SCHEMA(('${schema}'));\nENDSEC;\nDATA;\n${data}\nENDSEC;\nEND-ISO-10303-21;\n`;
const MM_CONTEXT = `#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n#2=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n#3=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n#4=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(())GLOBAL_UNIT_ASSIGNED_CONTEXT((#1,#2,#3))REPRESENTATION_CONTEXT('',''));`;

test('Part 21 tokenizer: strings, escapes, typed values, nested lists, complex instances, comments', () => {
  const text = p21('TEST', `/* a comment ; with a semicolon */\n#10 = THING('it''s; here', $, *, .T., 1.5E+2, -3, (#11, #12), (1, (2, 3)), LENGTH_MEASURE(25.4), 'caf\\X2\\00E9\\X0\\');\n#11=POINT('',(1.,2.,3.));\n#12 = ( A_PART ( 7 ) B_PART ( 'x' ) ) ;\n#13=EMPTY();`);
  const p = new P21(enc(text), 't.stp');
  assert.equal(p.count, 4);
  assert.deepEqual(p.header.schema, ['TEST']);
  assert.equal(p.header.originatingSystem, 'Test System 1.0');
  const e = p.get(10);
  assert.equal(e.type, 'THING');
  assert.equal(e.args[0], "it's; here");
  assert.deepEqual(e.args.slice(1, 6), [null, undefined, '.T.', 150, -3]);
  assert.ok(e.args[6][0] instanceof Ref && e.args[6][1].id === 12);
  assert.deepEqual(e.args[7], [1, [2, 3]]);
  assert.ok(e.args[8] instanceof Typed && e.args[8].type === 'LENGTH_MEASURE' && e.args[8].value === 25.4);
  assert.equal(e.args[9], 'café');
  assert.deepEqual(p.xyz(new Ref(11)), [1, 2, 3]);
  assert.deepEqual(p.get(12).parts, [{ type: 'A_PART', args: [7] }, { type: 'B_PART', args: ['x'] }]);
  assert.deepEqual(p.get(13).args, []);
  assert.equal(p.get(99), null);
  assert.deepEqual(p.idsOf('POINT'), [11]);
  assert.throws(() => new P21(enc('hello'), 'x.stp'), /expected the first line "ISO-10303-21;"/);
});

test('STEP AP242 tessellated face (millimetres) → surface in metres', () => {
  const data = `${MM_CONTEXT}\n#20=COORDINATES_LIST('',8,(${CUBE_V.map((v) => `(${v.map((c) => (c * 1000).toFixed(1)).join(',')})`).join(',')}));\n#21=TRIANGULATED_FACE('cube',#20,8,(),$,(),(${CUBE_T.map((t) => `(${t.map((k) => k + 1).join(',')})`).join(',')}));\n#22=TESSELLATED_SHELL('hull',(#21),$);\n#23=COMPLEX_TRIANGULATED_FACE('lid',#20,4,(),$,(5,6,7,8),((1,2,4,3)),());\n#24=PRODUCT('p1','Pontoon','',(#99));`;
  const r = parseStep({ name: 'pontoon.stp', bytes: enc(p21('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }', data)) });
  assertContract(r, ['surface']);
  const s = r.models[0];
  assert.deepEqual(s.groups.map((g) => [g.name, g.count / 3]), [['hull', 12], ['lid', 2]]);
  close(surfaceVolume({ positions: s.positions, indices: s.indices.subarray(0, 36) }), 1, 1e-6);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [1, 1, 1] });
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.equal(r.meta.needsAdapter, null);
  assert.deepEqual(r.meta.products, ['Pontoon']);
});

test('STEP faceted B-rep (poly loops) and wireframe curves', () => {
  let n = 100;
  const pts = CUBE_V.map((v) => `#${++n}=CARTESIAN_POINT('',(${v.join('.,')}.));`);
  const faces = CUBE_Q.map((q, i) => `#${200 + 3 * i}=POLY_LOOP('',(${q.map((k) => `#${101 + k}`).join(',')}));\n#${201 + 3 * i}=FACE_OUTER_BOUND('',#${200 + 3 * i},.T.);\n#${202 + 3 * i}=FACE_SURFACE('',(#${201 + 3 * i}),$,.T.);`);
  const f = parseStep({ name: 'box.step', text: p21('CONFIG_CONTROL_DESIGN', [...pts, ...faces, `#300=CLOSED_SHELL('',(${CUBE_Q.map((_, i) => `#${202 + 3 * i}`).join(',')}));`, `#301=FACETED_BREP('box',#300);`].join('\n')) }, { units: 'm' });
  assertContract(f, ['surface']);
  assert.equal(f.models[0].positions.length / 3, 8);
  close(surfaceVolume(f.models[0]), 1, 1e-6);
  assert.equal(f.meta.needsAdapter, null);

  const wire = `${MM_CONTEXT}
#10=CARTESIAN_POINT('',(0.,0.,0.));\n#11=CARTESIAN_POINT('',(1000.,0.,0.));\n#12=CARTESIAN_POINT('',(1000.,1000.,0.));
#13=POLYLINE('pl',(#10,#11,#12));
#14=DIRECTION('',(0.,0.,1.));\n#15=VECTOR('',#14,1000.);\n#16=LINE('',#10,#15);
#17=TRIMMED_CURVE('riser',#16,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(2.5)),.T.,.PARAMETER.);
#18=DIRECTION('',(1.,0.,0.));\n#19=AXIS2_PLACEMENT_3D('',#10,#14,#18);\n#20=CIRCLE('',#19,500.);
#21=TRIMMED_CURVE('bend',#20,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.5707963267948966)),.T.,.PARAMETER.);
#22=B_SPLINE_CURVE_WITH_KNOTS('spl',1,(#10,#11,#12),.UNSPECIFIED.,.F.,.F.,(2,1,2),(0.,1.,2.),.UNSPECIFIED.);
#23=GEOMETRIC_CURVE_SET('wires',(#13,#17,#21,#22,#12));`;
  const w = parseStep({ name: 'wire.stp', text: p21('AUTOMOTIVE_DESIGN', wire) });
  assertContract(w, ['polyline', 'points']);
  const [pl, riser, bend, spl] = w.models;
  closeArr([...pl.points], [0, 0, 0, 1, 0, 0, 1, 1, 0], 1e-12);
  closeArr([...riser.points], [0, 0, 0, 0, 0, 2.5], 1e-12);
  assert.equal(riser.name, 'riser');
  for (let i = 0; i < bend.points.length; i += 3) close(Math.hypot(bend.points[i], bend.points[i + 1]), 0.5, 1e-12);
  closeArr([...bend.points.slice(-3)], [0, 0.5, 0], 1e-12);
  closeArr([...spl.points.slice(-3)], [1, 1, 0], 1e-12);
  assert.equal(byKind(w, 'points').positions.length, 3);
});

test('STEP with exact B-rep only: error names the adapter; allowEmpty reports needsAdapter', () => {
  const data = `${MM_CONTEXT}\n#10=CARTESIAN_POINT('',(0.,0.,0.));\n#11=VERTEX_POINT('',#10);\n#12=EDGE_CURVE('',#11,#11,#13,.T.);\n#13=CIRCLE('',#14,5.);\n#14=AXIS2_PLACEMENT_3D('',#10,$,$);\n#15=ORIENTED_EDGE('',*,*,#12,.T.);\n#16=EDGE_LOOP('',(#15));\n#17=FACE_BOUND('',#16,.T.);\n#18=PLANE('',#14);\n#19=ADVANCED_FACE('',(#17),#18,.T.);\n#20=CLOSED_SHELL('',(#19));\n#21=MANIFOLD_SOLID_BREP('disc',#20);`;
  const text = p21('AUTOMOTIVE_DESIGN', data);
  assert.throws(() => parseStep({ name: 'solid.stp', text }), /exact B-rep geometry only \(1 MANIFOLD_SOLID_BREP, 1 ADVANCED_FACE, 1 CLOSED_SHELL\).*adapter "occt".*tessellated STEP \(AP242/);
  const r = parseStep({ name: 'solid.stp', text }, { allowEmpty: true });
  assert.equal(r.models.length, 0);
  assert.equal(r.meta.needsAdapter, 'occt');
  assert.equal(r.meta.skipped['B-rep faces (exact geometry)'], 1);
});

test('CIS/2: analysis elements on nodes and design located parts → frame', () => {
  const data = `#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#10=CARTESIAN_POINT('',(0.,0.,0.));\n#11=CARTESIAN_POINT('',(0.,0.,4000.));\n#12=CARTESIAN_POINT('',(6000.,0.,4000.));
#20=NODE('N1',#10,$,#90);\n#21=NODE('N2',#11,$,#90);\n#22=NODE('N3',#12,$,#90);
#30=SECTION_PROFILE(1,'HEB300',$,$,5,.T.);\n#31=SECTION_PROFILE(2,'W12X26',$,$,5,.T.);
#40=ELEMENT_CURVE_SIMPLE('COL1',$,#90,1,#30,(#99),8);\n#41=ELEMENT_CURVE_SIMPLE('BM1',$,#90,1,#31,(#99),8);
#50=ELEMENT_NODE_CONNECTIVITY(1,'start',#20,#40,$,$);\n#51=ELEMENT_NODE_CONNECTIVITY(2,'end',#21,#40,$,$);
#52=ELEMENT_NODE_CONNECTIVITY(2,'end',#22,#41,$,$);\n#53=ELEMENT_NODE_CONNECTIVITY(1,'start',#21,#41,$,$);
#60=POSITIVE_LENGTH_MEASURE_WITH_UNIT(POSITIVE_LENGTH_MEASURE(3000.),#1);
#61=PART_PRISMATIC_SIMPLE(1,'brace','',$,.ROLLED.,$,#62,#60,$,$);\n#62=SECTION_PROFILE(3,'PIPE 219.1x8',$,$,5,.T.);
#70=CARTESIAN_POINT('',(6000.,0.,4000.));\n#71=DIRECTION('',(1.,0.,0.));\n#72=DIRECTION('',(0.,0.,-1.));
#73=AXIS2_PLACEMENT_3D('',#70,#71,#72);\n#74=COORD_SYSTEM_CARTESIAN_3D('cs','',$,3,#73);\n#75=LOCATED_PART(1,'BR1',$,#74,#61,#98);`;
  const r = parseCis2({ name: 'frame.stp', text: p21('STRUCTURAL_FRAME_SCHEMA', data) });
  assertContract(r, ['frame']);
  const f = r.models[0];
  assert.equal(f.members.length, 3);
  assert.equal(f.nodes.length / 3, 4);
  const [col, beam, brace] = f.members;
  assert.equal(col.name, 'COL1');
  closeArr([f.nodes[3 * col.a + 2], f.nodes[3 * col.b + 2]], [0, 4], 1e-12);
  assert.equal(col.section.designation, 'HE 300 B');
  closeArr([f.nodes[3 * beam.a], f.nodes[3 * beam.b]], [0, 6], 1e-12);
  // the located part runs 3 m along its local x axis (the placement's ref direction), here straight down
  assert.equal(brace.a, beam.b, 'the part starts on the beam end node');
  closeArr([f.nodes[3 * brace.b], f.nodes[3 * brace.b + 2]], [6, 1], 1e-12);
  close(brace.section.D, 0.2191, 1e-12);
  assert.deepEqual(f.nodeNames.slice(0, 3), ['N1', 'N2', 'N3']);
  assert.equal(r.meta.units.unit, 'mm');
  assert.throws(() => parseCis2({ name: 'x.stp', text: p21('STRUCTURAL_FRAME_SCHEMA', '#1=CARTESIAN_POINT(\'\',(0.,0.,0.));') }), /no members were found/);
});

// ---- IFC ------------------------------------------------------------------------------------------------

const IFC = `#1=IFCPROJECT('2O2Fr$t4X7Zf8NOew3FLOH',$,'Jetty',$,$,$,$,(#20),#7);
#5=IFCSIUNIT(*,.LENGTHUNIT.,.MILLI.,.METRE.);\n#6=IFCSIUNIT(*,.PLANEANGLEUNIT.,$,.RADIAN.);\n#7=IFCUNITASSIGNMENT((#5,#6));
#10=IFCCARTESIANPOINT((0.,0.,0.));\n#11=IFCDIRECTION((0.,0.,1.));\n#12=IFCDIRECTION((1.,0.,0.));\n#13=IFCAXIS2PLACEMENT3D(#10,#11,#12);
#20=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-5,#13,$);
#30=IFCLOCALPLACEMENT($,#13);
#31=IFCCARTESIANPOINT((0.,0.,3000.));\n#32=IFCAXIS2PLACEMENT3D(#31,$,$);\n#33=IFCLOCALPLACEMENT(#30,#32);
#34=IFCBUILDINGSTOREY('3h9Lmn0000000000000001',$,'Deck level',$,$,#33,$,$,.ELEMENT.,3000.);
#40=IFCCARTESIANPOINT((1000.,2000.,0.));\n#41=IFCAXIS2PLACEMENT3D(#40,$,$);\n#42=IFCLOCALPLACEMENT(#33,#41);
#50=IFCCARTESIANPOINT((0.,0.));\n#51=IFCCARTESIANPOINT((6000.,0.));\n#52=IFCPOLYLINE((#50,#51));
#53=IFCSHAPEREPRESENTATION(#20,'Axis','Curve2D',(#52));
#54=IFCISHAPEPROFILEDEF(.AREA.,'IPE300',$,150.,300.,7.1,10.7,15.,$,$);
#55=IFCDIRECTION((0.,1.,0.));\n#56=IFCAXIS2PLACEMENT3D(#10,#12,#55);
#57=IFCEXTRUDEDAREASOLID(#54,#56,#11,6000.);
#58=IFCSHAPEREPRESENTATION(#20,'Body','SweptSolid',(#57));
#59=IFCPRODUCTDEFINITIONSHAPE($,$,(#53,#58));
#60=IFCBEAM('0BTBFw6f90Nfh9rP1dlXr2',$,'B-1','deck beam','Beam:IPE300',#42,#59,'T-100',.BEAM.);
#70=IFCCIRCLEHOLLOWPROFILEDEF(.AREA.,'CHS 610x19',$,305.,19.);
#71=IFCDIRECTION((0.,0.,-1.));\n#72=IFCEXTRUDEDAREASOLID(#70,#13,#71,12000.);
#73=IFCSHAPEREPRESENTATION(#20,'Body','SweptSolid',(#72));\n#74=IFCPRODUCTDEFINITIONSHAPE($,$,(#73));
#75=IFCPILE('1kTvXnbbzCWw8lcMd1dR4o',$,'P-1',$,$,#42,#74,$,.DRIVEN.,$);
#80=IFCMATERIAL('S355',$,$);\n#81=IFCRELASSOCIATESMATERIAL('2YBqaV_8L15eWJ9DA1sGmT',$,$,$,(#60,#75),#80);
#82=IFCRELCONTAINEDINSPATIALSTRUCTURE('3Sa3dTJGn0H8TQIGiuGQd5',$,$,$,(#60,#75,#90),#34);
#90=IFCPLATE('2gRXFgjRn2HPE$YoDLX3FV',$,'PL-1',$,$,#42,$,$,.NOTDEFINED.);
#95=IFCCARTESIANPOINT((0.,0.,0.));\n#96=IFCCARTESIANPOINT((100000.,0.,0.));\n#97=IFCCARTESIANPOINT((100000.,50000.,0.));\n#98=IFCPOLYLINE((#95,#96,#97));
#99=IFCSHAPEREPRESENTATION(#20,'Axis','Curve3D',(#98));\n#100=IFCPRODUCTDEFINITIONSHAPE($,$,(#99));
#101=IFCALIGNMENT('0Hk8sMX6b4T9qhv3Ljd6eP',$,'Approach trestle',$,$,#30,#100,$);`;

test('IFC: beam (Axis + profile), pile (extruded Body), placements, material, container, alignment', () => {
  const r = parseIfc({ name: 'jetty.ifc', bytes: enc(p21('IFC4', IFC)) });
  assertContract(r, ['frame', 'polyline']);
  const f = byKind(r, 'frame');
  assert.equal(f.members.length, 2);
  const [beam, pile] = f.members;
  // storey at z = 3000 mm, element placement (1000, 2000, 0) mm
  closeArr([...f.nodes].slice(3 * beam.a, 3 * beam.a + 3), [1, 2, 3], 1e-12);
  closeArr([...f.nodes].slice(3 * beam.b, 3 * beam.b + 3), [7, 2, 3], 1e-12);
  assert.equal(beam.group, 'IfcBeam');
  assert.equal(beam.name, 'B-1');
  assert.equal(beam.tag, 'T-100');
  assert.equal(beam.globalId, '0BTBFw6f90Nfh9rP1dlXr2');
  assert.equal(beam.material, 'S355');
  assert.equal(beam.container, 'Deck level');
  assert.equal(beam.section.shape, 'I');
  // IPE 300: A = 53.81 cm², Iy = 8356 cm⁴ (EN 10365 tables)
  close(beam.section.A * 1e4, 53.81, 0.05);
  close(beam.section.Iy * 1e8, 8356, 10);
  assert.equal(pile.group, 'IfcPile');
  assert.equal(pile.a, beam.a, 'pile head welds to the beam start');
  closeArr([...f.nodes].slice(3 * pile.b, 3 * pile.b + 3), [1, 2, -9], 1e-12);
  assert.deepEqual(pile.section, { type: 'tube', D: 0.61, t: 0.019, name: 'CHS 610x19' });
  const al = byKind(r, 'polyline');
  assert.equal(al.name, 'Approach trestle');
  closeArr([...al.attrs.kp], [0, 100, 150], 1e-9);
  assert.equal(r.meta.project, 'Jetty');
  assert.deepEqual(r.meta.classes, { IfcBeam: 1, IfcPile: 1 });
  assert.deepEqual(r.meta.otherProducts, { IFCPLATE: 1 });
  assert.equal(r.meta.needsAdapter, 'web-ifc');
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.ok(r.warnings.some((w) => /1 other product\(s\) \(IFCPLATE\)/.test(w)));
});

test('IFC: files without linear elements or with another schema explain themselves', () => {
  const only = `#5=IFCSIUNIT(*,.LENGTHUNIT.,$,.METRE.);\n#7=IFCUNITASSIGNMENT((#5));\n#90=IFCSLAB('2gRXFgjRn2HPE$YoDLX3FV',$,'S-1',$,$,$,$,$,.FLOOR.);`;
  assert.throws(() => parseIfc({ name: 's.ifc', text: p21('IFC2X3', only) }), /no beams, columns, members, piles.*1 other products.*adapter "web-ifc"/);
  assert.equal(parseIfc({ name: 's.ifc', text: p21('IFC2X3', only) }, { allowEmpty: true }).meta.needsAdapter, 'web-ifc');
  assert.throws(() => parseIfc({ name: 's.ifc', text: p21('AUTOMOTIVE_DESIGN', only) }), /not an IFC file \(FILE_SCHEMA is AUTOMOTIVE_DESIGN/);
});

// ---- IGES -----------------------------------------------------------------------------------------------

function iges(entities, unitFlag = 2) {
  const lines = [];
  const rec = (body, sec, n) => body.padEnd(72) + sec + String(n).padStart(7);
  lines.push(rec('unit test file', 'S', 1));
  const g = ['1H,', '1H;', '4Htest', '8Htest.igs', '4Hunit', '3H1.0', '32', '38', '6', '308', '15', '4Htest', '1.0', String(unitFlag), unitFlag === 2 ? '2HMM' : '1HM', '1', '0.01', '15H20240101.000000', '0.001', '1000.0', '1Ha', '1Ho', '11', '0'].join(',') + ';';
  for (let i = 0, n = 1; i < g.length; i += 72, n++) lines.push(rec(g.slice(i, i + 72), 'G', n));
  const dlines = [], plines = [];
  const f8 = (v) => String(v).padStart(8);
  entities.forEach((e, k) => {
    const seq = 2 * k + 1;
    const text = `${e.type},${e.params.join(',')};`;
    const chunks = [];
    // break only at commas so that numbers are never split across records
    let cur = '';
    for (const part of text.split(/(?<=,)/)) { if (cur.length + part.length > 64) { chunks.push(cur); cur = ''; } cur += part; }
    chunks.push(cur);
    const pstart = plines.length + 1;
    chunks.forEach((c) => plines.push(c.padEnd(64) + ' ' + String(seq).padStart(7) + 'P' + String(plines.length + 1).padStart(7)));
    dlines.push(f8(e.type) + f8(pstart) + f8(0) + f8(0) + f8(e.level || 0) + f8(0) + f8(e.xform ? 2 * e.xform - 1 : 0) + f8(0) + (e.status || '00000000') + 'D' + String(seq).padStart(7));
    dlines.push(f8(e.type) + f8(0) + f8(0) + f8(chunks.length) + f8(e.form || 0) + f8('') + f8('') + (e.label || '').padStart(8) + f8(0) + 'D' + String(seq + 1).padStart(7));
  });
  lines.push(...dlines, ...plines, rec(`S${String(1).padStart(7)}G${String(3).padStart(7)}D${String(dlines.length).padStart(7)}P${String(plines.length).padStart(7)}`, 'T', 1));
  return lines.join('\n') + '\n';
}

test('IGES fields: Hollerith strings may contain delimiters', () => {
  assert.deepEqual(igesFields('1H,,1H;,5Ha,b;c,12,3.5D2;rest').fields, [',', ';', 'a,b;c', '12', '3.5D2']);
});

test('IGES: lines, arc with transformation matrix, points, B-spline curve and surface, units', () => {
  const w = Math.SQRT1_2;
  const text = iges([
    { type: 110, params: [0, 0, 0, 0, 0, '1.0D4'], level: 3, label: 'LEG' },
    { type: 110, params: [0, 0, 10000, 5000, 0, 10000], level: 4 },
    { type: 124, params: [1, 0, 0, 1000, 0, 1, 0, 0, 0, 0, 1, 2000], status: '00010000' },
    { type: 100, params: [0, 0, 0, 500, 0, 0, 500], xform: 3, label: 'BEND' },
    { type: 116, params: [1000, 2000, 3000] },
    { type: 126, params: [2, 2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, w, 1, 1000, 0, 0, 1000, 1000, 0, 0, 1000, 0, 0, 1, 0, 0, 1], label: 'ARC' },
    { type: 128, params: [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2000, 0, 0, 0, 3000, 0, 2000, 3000, 0, 0, 1, 0, 1], label: 'PATCH' },
    { type: 144, params: [13, 0, 0, 0] },
    { type: 126, params: [1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1], status: '00010500' },
    { type: 212, params: [1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '4Hnote'], status: '00000100' },
  ]);
  const r = parseIges({ name: 'model.igs', bytes: enc(text) });
  assertContract(r, ['frame', 'surface', 'polyline', 'points']);
  assert.equal(r.meta.units.unit, 'mm');
  const f = byKind(r, 'frame');
  assert.equal(f.nodes.length / 3, 3);
  assert.deepEqual(f.members.map((m) => m.group), ['level 3', 'level 4']);
  closeArr([...f.nodes], [0, 0, 0, 0, 0, 10, 5, 0, 10], 1e-12);
  const [bend, arc] = r.models.filter((m) => m.kind === 'polyline');
  // quarter circle of radius 0.5 m, shifted by the matrix to (1, 0, 2)
  closeArr([...bend.points.slice(0, 3), ...bend.points.slice(-3)], [1.5, 0, 2, 1, 0.5, 2], 1e-9);
  for (let i = 0; i < arc.points.length; i += 3) close(Math.hypot(arc.points[i], arc.points[i + 1]), 1, 1e-9);
  assert.equal(arc.name, 'ARC');
  close(surfaceArea(byKind(r, 'surface')), 6, 1e-9);
  closeArr([...byKind(r, 'points').positions], [1, 2, 3], 1e-9);
  assert.equal(r.meta.needsAdapter, 'occt');
  assert.deepEqual(r.meta.kernelEntities, { '144 trimmed surface': 1 });
  assert.equal(r.meta.skipped['126 (dependent or non-geometry use)'], 1);
  assert.equal(r.meta.skipped['entity 212'], 1);
  assert.ok(r.warnings.some((x) => /1 entity needs a geometry kernel/.test(x)));
});

test('IGES: closed planar curves also give a section; kernel-only and malformed files', () => {
  const sq = iges([{ type: 106, form: 63, params: [1, 5, 0, 0, 0, 2, 0, 2, 1, 0, 1, 0, 0], label: 'OUTER' }], 6);
  const r = parseIges({ name: 'p.iges', text: sq });
  assert.deepEqual(r.models.map((m) => m.kind), ['section', 'polyline']);
  close(Math.abs(ringArea(r.models[0].loops[0])), 2, 1e-12);
  assert.equal(r.meta.units.unit, 'm');
  const k = iges([{ type: 144, params: [3, 1, 0, 5] }, { type: 120, params: [1, 3, 0, 6.28] }]);
  assert.throws(() => parseIges({ name: 'k.igs', text: k }), /only entities that need a geometry kernel \(144 trimmed surface, 120 surface of revolution\).*adapter "occt"/);
  assert.equal(parseIges({ name: 'k.igs', text: k }, { allowEmpty: true }).meta.needsAdapter, 'occt');
  assert.throws(() => parseIges({ name: 'x.igs', text: 'not an iges file\n' + ' '.repeat(80) + '\n' }), /expected section letter "S" in column 73/);
});

// ---- SVG and airfoils -----------------------------------------------------------------------------------

test('SVG: physical size from width + viewBox, y flipped, arcs and transforms, layers, open path', () => {
  const sub = flattenPath('M0 0 h10 a5 5 0 0 1 0 10 h-10 z m20 0 l5 5');
  assert.equal(sub.length, 2);
  assert.ok(sub[0].closed && !sub[1].closed);
  const svg = `<?xml version="1.0"?><svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" width="200mm" height="100mm" viewBox="0 0 100 50">
  <defs><circle id="unused" r="99"/></defs>
  <g inkscape:label="sheath" id="g1"><path d="M 0,0 H 10 A 5,5 0 0 1 10,10 H 0 Z"/></g>
  <g id="cores" transform="translate(50,0) scale(2)"><circle cx="5" cy="5" r="2"/><rect x="10" y="0" width="4" height="3"/></g>
  <polygon points="80,40 90,40 90,50 80,50"/><path d="M0 40 L10 45"/><text x="1" y="1">note</text></svg>`;
  const r = parseSvg({ name: 'umbilical.svg', text: svg });
  assertContract(r, ['section', 'polyline']);
  const s = byKind(r, 'section');
  assert.deepEqual(s.layers, [{ name: 'sheath', loops: [0] }, { name: 'cores', loops: [1, 2] }, { name: 'drawing', loops: [3] }]);
  // 1 user unit = 2 mm
  close(Math.abs(ringArea(s.loops[0])), (100 + (Math.PI * 25) / 2) * 4e-6, 2e-7);
  close(Math.abs(ringArea(s.loops[1])), Math.PI * 16 * 4e-6, 3e-7, 'circle of radius 2, scaled by 2 (72-gon)');
  close(Math.abs(ringArea(s.loops[2])), 4 * 3 * 4 * 4e-6, 1e-12);
  // the polygon at the bottom of the drawing ends up at the bottom (y ≈ 0) after the flip
  const yMin = Math.min(...s.loops[3].filter((_, i) => i % 2 === 1));
  close(yMin, 0, 1e-12);
  close(Math.max(...s.loops[3].filter((_, i) => i % 2 === 0)), 0.18, 1e-12);
  assert.equal(r.meta.units.basis, 'file');
  close(r.meta.units.scale, 0.002, 1e-15);
  assert.equal(r.meta.skipped['<text>'], 1);
  const px = parseSvg({ name: 'p.svg', text: '<svg><rect width="96" height="96"/></svg>' });
  close(Math.abs(ringArea(px.models[0].loops[0])), 0.0254 ** 2, 1e-12);
  assert.ok(px.warnings.some((w) => /one user unit was taken as one CSS pixel/.test(w)));
  close(Math.abs(ringArea(parseSvg({ name: 'p.svg', text: '<svg><rect width="96" height="96"/></svg>' }, { units: 'mm' }).models[0].loops[0])), 0.096 ** 2, 1e-12);
  assert.throws(() => parseSvg({ name: 'x.svg', text: '<kml/>' }), /expected root element <svg>/);
  assert.throws(() => parseSvg({ name: 'x.svg', text: '<svg><text>a</text></svg>' }), /no outlines were found/);
});

// NACA 4-digit symmetric thickness distribution (Abbott & von Doenhoff), closed trailing edge variant not used
const nacaT = (x, t = 0.12) => 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1015 * x ** 4);

test('Airfoils: Selig and Lednicer layouts give the same loop; NACA 0012 thickness; percent chord', () => {
  const xs = Array.from({ length: 41 }, (_, i) => (1 - Math.cos((Math.PI * i) / 40)) / 2);
  const selig = ['NACA 0012', ...xs.slice().reverse().map((x) => `${x.toFixed(6)} ${nacaT(x).toFixed(6)}`), ...xs.slice(1).map((x) => `${x.toFixed(6)} ${(-nacaT(x)).toFixed(6)}`)].join('\n');
  const led = ['NACA 0012 (Lednicer)', '41. 41.', '', ...xs.map((x) => `${x.toFixed(6)} ${nacaT(x).toFixed(6)}`), '', ...xs.map((x) => `${x.toFixed(6)} ${(-nacaT(x)).toFixed(6)}`)].join('\r\n');
  const a = parseAirfoil({ name: 'n0012.dat', text: selig });
  const b = parseAirfoil({ name: 'n0012-il.dat', text: led });
  assertContract(a, ['section']);
  assert.equal(a.meta.airfoil.layout, 'selig');
  assert.equal(b.meta.airfoil.layout, 'lednicer');
  assert.equal(a.models[0].name, 'NACA 0012');
  closeArr([...a.models[0].loops[0]], [...b.models[0].loops[0]], 1e-12);
  close(a.meta.airfoil.maxThickness, 0.12, 2e-4);
  close(a.meta.airfoil.maxThicknessAt, 0.30, 0.02);
  close(a.meta.airfoil.maxCamber, 0, 1e-9);
  assert.ok(ringArea(a.models[0].loops[0]) > 0, 'upper surface first is counter-clockwise');
  const scaled = parseAirfoil({ name: 'n0012.dat', text: selig.replace(/^(\S+) (\S+)$/gm, (m, x, y) => `${x * 100} ${y * 100}`) }, { chord: 3 });
  close(scaled.meta.airfoil.maxThickness, 0.12, 2e-4);
  close(Math.max(...scaled.models[0].loops[0]), 3, 1e-9);
  assert.ok(scaled.warnings.some((w) => /percent of chord/.test(w)));
  assert.throws(() => parseAirfoil({ name: 'x.dat', text: 'name only\n1 2\n' }), /at least five "x y" coordinate pairs/);
});

// ---- panel meshes ---------------------------------------------------------------------------------------

function outward(m, centre) {
  const p = m.positions, t = m.indices;
  for (let i = 0; i < t.length; i += 3) {
    const a = 3 * t[i], b = 3 * t[i + 1], c = 3 * t[i + 2];
    const n = [(p[b + 1] - p[a + 1]) * (p[c + 2] - p[a + 2]) - (p[b + 2] - p[a + 2]) * (p[c + 1] - p[a + 1]), (p[b + 2] - p[a + 2]) * (p[c] - p[a]) - (p[b] - p[a]) * (p[c + 2] - p[a + 2]), (p[b] - p[a]) * (p[c + 1] - p[a + 1]) - (p[b + 1] - p[a + 1]) * (p[c] - p[a])];
    const g = [(p[a] + p[b] + p[c]) / 3 - centre[0], (p[a + 1] + p[b + 1] + p[c + 1]) / 3 - centre[1], (p[a + 2] + p[b + 2] + p[c + 2]) / 3 - centre[2]];
    assert.ok(n[0] * g[0] + n[1] * g[1] + n[2] * g[2] > 0, `triangle ${i / 3} faces outward`);
  }
}

// quarter of the wetted surface of a 2 × 2 × 1 box (x, y ≥ 0), normals pointing out of the body
const QUARTER = [[[0, 0, -1], [0, 1, -1], [1, 1, -1], [1, 0, -1]], [[1, 0, -1], [1, 1, -1], [1, 1, 0], [1, 0, 0]], [[0, 1, -1], [0, 1, 0], [1, 1, 0], [1, 1, -1]]];

test('WAMIT GDF: two planes of symmetry expanded with outward normals; unit from GRAV', () => {
  const text = `Box quarter\n1.0  9.80665\n1  1\n3\n${QUARTER.map((q) => q.map((v) => v.join(' ')).join('\n')).join('\n')}\n`;
  const r = parseGdf({ name: 'box.gdf', bytes: enc(text) });
  assertContract(r, ['surface']);
  const m = r.models[0];
  assert.equal(m.indices.length / 3, 24);
  assert.equal(m.quads.length / 4, 12);
  close(surfaceArea(m), 12, 1e-9);
  assert.deepEqual(r.meta.bbox, { min: [-1, -1, -1], max: [1, 1, 0] });
  outward(m, [0, 0, -0.5]);
  assert.deepEqual(m.groups.map((g) => g.name), ['as stored', 'mirror x', 'mirror y', 'mirror xy']);
  assert.equal(m.positions.length / 3, 17, 'vertices on the symmetry planes are shared: 9 on the bottom, 8 on the waterline');
  assert.deepEqual(r.meta.symmetry, { x: true, y: true, expanded: true, note: 'x: mirror plane x = 0; y: mirror plane y = 0' });
  assert.equal(r.meta.units.unit, 'm');
  assert.equal(r.meta.units.basis, 'file');
  const half = parseGdf({ name: 'box.gdf', text }, { expandSymmetry: false });
  assert.equal(half.models[0].indices.length / 3, 6);
  const ft = parseGdf({ name: 'box.gdf', text: text.replace('9.80665', '32.174') });
  close(ft.meta.bbox.max[0], 0.3048, 1e-7);
  assert.throws(() => parseGdf({ name: 'ho.gdf', text: 'hdr\n1.0 9.80665\n0 0\n4 1\n' }), /higher-order GDF/);
  assert.throws(() => parseGdf({ name: 'x.gdf', text: 'hdr\nabc def\n' }), /file ends before the "ISX ISY" line|not a WAMIT low-order GDF/);
});

test('NEMOH mesh: node table, panels, xOz symmetry', () => {
  const nodes = [];
  const id = (v) => { let k = nodes.findIndex((n) => n.every((c, i) => c === v[i])); if (k < 0) { nodes.push(v); k = nodes.length - 1; } return k + 1; };
  // half body (y ≥ 0): the quarter and its mirror in x
  const panels = [...QUARTER, ...QUARTER.map((q) => q.map((v) => [-v[0], v[1], v[2]]).reverse())].map((q) => q.map(id));
  const text = `    2    1\n${nodes.map((n, i) => `${i + 1} ${n.map((c) => c.toFixed(4)).join(' ')}`).join('\n')}\n    0   0.00 0.00 0.00\n${panels.map((p) => p.join(' ')).join('\n')}\n   0   0   0   0\n`;
  const r = parseNemoh({ name: 'box.dat', text });
  assertContract(r, ['surface']);
  close(surfaceArea(r.models[0]), 12, 1e-9);
  outward(r.models[0], [0, 0, -0.5]);
  assert.deepEqual(r.meta.symmetry.y, true);
  assert.equal(r.meta.panelsStored, 6);
  assert.deepEqual(r.meta.elementTypes, { quad4: 12 });
  assert.throws(() => parseNemoh({ name: 'x.dat', text: '3 0\n' }), /expected the first line "2 0" or "2 1"/);
});

test('truncated CAD / structural text files fail with a worded error, never an internal one', async () => {
  const d = dxfFile({
    header: [9, '$INSUNITS', 70, 4],
    blocks: [0, 'BLOCK', 2, 'B', 10, 0, 20, 0, 30, 0, ...LINE('0', [0, 0, 0], [1, 0, 0]), 0, 'ENDBLK'],
    entities: [...LINE('A', [0, 0, 0], [0, 0, 5]), 0, 'LWPOLYLINE', 8, 'P', 90, 3, 70, 1, 10, 0, 20, 0, 42, 0.5, 10, 2, 20, 0, 10, 2, 20, 1, 0, 'INSERT', 8, 'I', 2, 'B', 10, 5, 20, 5, 30, 0,
      0, 'POLYLINE', 8, 'M', 66, 1, 70, 64, 0, 'VERTEX', 8, 'M', 10, 0, 20, 0, 30, 0, 70, 192, 0, 'VERTEX', 8, 'M', 10, 1, 20, 0, 30, 0, 70, 192, 0, 'VERTEX', 8, 'M', 10, 1, 20, 1, 30, 0, 70, 192, 0, 'VERTEX', 8, 'M', 70, 128, 71, 1, 72, 2, 73, 3, 0, 'SEQEND',
      0, 'SPLINE', 8, 'S', 70, 0, 71, 2, 40, 0, 40, 0, 40, 0, 40, 1, 40, 1, 40, 1, 10, 1, 20, 0, 30, 0, 10, 1, 20, 1, 30, 0, 10, 0, 20, 1, 30, 0],
  });
  await assertTruncationSafe(enc(d), (b) => parseDxf({ name: 'd.dxf', bytes: b }), 'DXF');
  const sdnf = 'Packet 00\n"a" "b"\nPacket 10\n"millimeters" 1\n1 0 0 0 "BEAM" "B1" 0\n"W12X26" "A36" 0.0 0 0\n0 0 1 0 0 0 6000 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0 0 0 0 0 0 0\nPacket 20\n"millimeters" "millimeters" 1\n7 0 0 0 "PLATE" "PL-1" "S355" 25.0 4\n0 0 0\n1000 0 0\n1000 500 0\n0 500 0\n';
  await assertTruncationSafe(enc(sdnf), (b) => parseSdnf({ name: 'm.sdnf', bytes: b }), 'SDNF');
  await assertTruncationSafe(enc(p21('IFC4', IFC)), (b) => parseIfc({ name: 'j.ifc', bytes: b }), 'IFC');
  const tess = `${MM_CONTEXT}\n#20=COORDINATES_LIST('',8,(${CUBE_V.map((v) => `(${v.join('.,')}.)`).join(',')}));\n#21=TRIANGULATED_FACE('cube',#20,8,(),$,(),(${CUBE_T.map((t) => `(${t.map((k) => k + 1).join(',')})`).join(',')}));\n#30=CARTESIAN_POINT('',(0.,0.,0.));\n#31=CARTESIAN_POINT('',(1.,0.,0.));\n#32=POLYLINE('',(#30,#31));\n#33=DIRECTION('',(0.,0.,1.));\n#34=AXIS2_PLACEMENT_3D('',#30,#33,$);\n#35=CIRCLE('',#34,5.);\n#36=TRIMMED_CURVE('',#35,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.);\n#37=GEOMETRIC_CURVE_SET('',(#32,#36));`;
  await assertTruncationSafe(enc(p21('AP242', tess)), (b) => parseStep({ name: 's.stp', bytes: b }), 'STEP');
  const cis = `#1=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n#10=CARTESIAN_POINT('',(0.,0.,0.));\n#11=CARTESIAN_POINT('',(0.,0.,4000.));\n#20=NODE('N1',#10,$,#90);\n#21=NODE('N2',#11,$,#90);\n#30=SECTION_PROFILE(1,'HEB300',$,$,5,.T.);\n#40=ELEMENT_CURVE_SIMPLE('C',$,#90,1,#30,(#99),8);\n#50=ELEMENT_NODE_CONNECTIVITY(1,'s',#20,#40,$,$);\n#51=ELEMENT_NODE_CONNECTIVITY(2,'e',#21,#40,$,$);`;
  await assertTruncationSafe(enc(p21('STRUCTURAL_FRAME_SCHEMA', cis)), (b) => parseCis2({ name: 'c.stp', bytes: b }), 'CIS/2');
  const ig = iges([{ type: 110, params: [0, 0, 0, 0, 0, 1] }, { type: 124, params: [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 2], status: '00010000' }, { type: 100, params: [0, 0, 0, 5, 0, 0, 5], xform: 2 },
    { type: 126, params: [2, 2, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0.7, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1] }, { type: 128, params: [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 2, 0, 0, 0, 3, 0, 2, 3, 0, 0, 1, 0, 1] }, { type: 102, params: [2, 1, 5] }]);
  await assertTruncationSafe(enc(ig), (b) => parseIges({ name: 'm.igs', bytes: b }), 'IGES');
  const svg = '<svg width="10mm" height="10mm" viewBox="0 0 10 10"><g id="a" transform="translate(1,1) rotate(30)"><path d="M0 0 h5 a2 2 0 0 1 0 4 q1 1 2 2 c1 1 2 2 3 3 z"/><circle cx="2" cy="2" r="1"/></g><polygon points="0,0 1,0 1,1"/></svg>';
  await assertTruncationSafe(enc(svg), (b) => parseSvg({ name: 'p.svg', bytes: b }), 'SVG');
  const gdf = `box\n1.0 9.80665\n1 1\n3\n${QUARTER.map((q) => q.map((v) => v.join(' ')).join('\n')).join('\n')}\n`;
  await assertTruncationSafe(enc(gdf), (b) => parseGdf({ name: 'b.gdf', bytes: b }), 'GDF');
  await assertTruncationSafe(enc('2 1\n1 0 0 -1\n2 1 0 -1\n3 1 1 -1\n4 0 1 -1\n0 0 0 0\n1 2 3 4\n0 0 0 0\n'), (b) => parseNemoh({ name: 'b.dat', bytes: b }), 'NEMOH');
  await assertTruncationSafe(enc('NACA\n1 0\n0.5 0.06\n0 0\n0.5 -0.06\n1 0\n'), (b) => parseAirfoil({ name: 'a.dat', bytes: b }), 'airfoil');
});
