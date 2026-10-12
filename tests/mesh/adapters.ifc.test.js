import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import * as ifc from '../../src/geometry/adapters/web-ifc.js';
import { importFile, GeometryFormatError } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { repair } from '../../src/mesh/repair.js';
import { close, rel } from './helpers.js';

// ---------------------------------------------------------------------------------------------------------
// IFC models written out by hand (ISO 10303-21), small enough for every expected number to be worked out on
// paper: extruded rectangles and a circle, a wall with an opening, and one model in feet.
// ---------------------------------------------------------------------------------------------------------

const head = (schema) => `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1');\nFILE_NAME('m.ifc','2024-01-01T00:00:00',(''),(''),'','','');\nFILE_SCHEMA(('${schema}'));\nENDSEC;\nDATA;\n`;
const SI = (prefix) => `#1=IFCSIUNIT(*,.LENGTHUNIT.,${prefix ? `.${prefix}.` : '$'},.METRE.);\n#3=IFCUNITASSIGNMENT((#1));\n`;
const FOOT = "#1=IFCSIUNIT(*,.LENGTHUNIT.,$,.METRE.);\n#2=IFCDIMENSIONALEXPONENTS(1,0,0,0,0,0,0);\n#12=IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(0.3048),#1);\n#13=IFCCONVERSIONBASEDUNIT(#2,.LENGTHUNIT.,'FOOT',#12);\n#3=IFCUNITASSIGNMENT((#13));\n";
// world placement #11, z axis #5, x axis #6, 2-D profile placement #32
const CTX = "#4=IFCCARTESIANPOINT((0.,0.,0.));\n#5=IFCDIRECTION((0.,0.,1.));\n#6=IFCDIRECTION((1.,0.,0.));\n#7=IFCAXIS2PLACEMENT3D(#4,#5,#6);\n#8=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-5,#7,$);\n#10=IFCPROJECT('0YvctVUKr0kugbFTf53O9L',$,'Project',$,$,$,$,(#8),#3);\n#11=IFCLOCALPLACEMENT($,#7);\n#30=IFCCARTESIANPOINT((0.,0.));\n#31=IFCDIRECTION((1.,0.));\n#32=IFCAXIS2PLACEMENT2D(#30,#31);\n";
const TAIL = 'ENDSEC;\nEND-ISO-10303-21;\n';

// millimetres: a 200 × 400 beam, 6000 long, standing at (1000, 2000, 3000); a Ø1000 column, 2000 high, at the origin
const FRAME = `${head('IFC4')}${SI('MILLI')}${CTX}`
  + "#20=IFCCARTESIANPOINT((1000.,2000.,3000.));\n#21=IFCAXIS2PLACEMENT3D(#20,#5,#6);\n#22=IFCLOCALPLACEMENT(#11,#21);\n#33=IFCRECTANGLEPROFILEDEF(.AREA.,'R200x400',#32,200.,400.);\n#34=IFCEXTRUDEDAREASOLID(#33,#7,#5,6000.);\n#35=IFCSHAPEREPRESENTATION(#8,'Body','SweptSolid',(#34));\n#36=IFCPRODUCTDEFINITIONSHAPE($,$,(#35));\n#40=IFCBEAM('3cUkl32yn9qRSPvBJVyWYp',$,'B-101',$,$,#22,#36,'tag-1',.BEAM.);\n"
  + "#50=IFCCIRCLEPROFILEDEF(.AREA.,'D1000',#32,500.);\n#51=IFCEXTRUDEDAREASOLID(#50,#7,#5,2000.);\n#52=IFCSHAPEREPRESENTATION(#8,'Body','SweptSolid',(#51));\n#53=IFCPRODUCTDEFINITIONSHAPE($,$,(#52));\n#54=IFCCOLUMN('1cUkl32yn9qRSPvBJVyWYq',$,'C-1',$,$,#11,#53,$,.COLUMN.);\n"
  + TAIL;

// metres: a wall 4 long, 0.2 thick, 3 high, with an opening 1 wide and 1 high cut through it from z = 1
const WALL = `${head('IFC2X3')}${SI('')}${CTX}`
  + "#33=IFCRECTANGLEPROFILEDEF(.AREA.,'W',#32,4.,0.2);\n#34=IFCEXTRUDEDAREASOLID(#33,#7,#5,3.);\n#35=IFCSHAPEREPRESENTATION(#8,'Body','SweptSolid',(#34));\n#36=IFCPRODUCTDEFINITIONSHAPE($,$,(#35));\n#40=IFCWALLSTANDARDCASE('2cUkl32yn9qRSPvBJVyWYp',$,'W-1',$,$,#11,#36,$);\n"
  + "#60=IFCCARTESIANPOINT((0.,0.,1.));\n#61=IFCAXIS2PLACEMENT3D(#60,#5,#6);\n#62=IFCLOCALPLACEMENT(#11,#61);\n#63=IFCRECTANGLEPROFILEDEF(.AREA.,'O',#32,1.,0.6);\n#64=IFCEXTRUDEDAREASOLID(#63,#7,#5,1.);\n#65=IFCSHAPEREPRESENTATION(#8,'Body','SweptSolid',(#64));\n#66=IFCPRODUCTDEFINITIONSHAPE($,$,(#65));\n#67=IFCOPENINGELEMENT('0cUkl32yn9qRSPvBJVyWYr',$,'O-1',$,$,#62,#66,$);\n#68=IFCRELVOIDSELEMENT('3dUkl32yn9qRSPvBJVyWYs',$,$,$,#40,#67);\n"
  + TAIL;

// feet (conversion-based unit): a slab profile of 2 × 4 ft extruded 10 ft
const FEET = `${head('IFC2X3')}${FOOT}${CTX}`
  + "#33=IFCRECTANGLEPROFILEDEF(.AREA.,'R',#32,2.,4.);\n#34=IFCEXTRUDEDAREASOLID(#33,#7,#5,10.);\n#35=IFCSHAPEREPRESENTATION(#8,'Body','SweptSolid',(#34));\n#36=IFCPRODUCTDEFINITIONSHAPE($,$,(#35));\n#40=IFCSLAB('3cUkl32yn9qRSPvBJVyWYp',$,'S-1',$,$,#11,#36,$,.FLOOR.);\n"
  + TAIL;

const enc = (s) => new TextEncoder().encode(s);
const byName = (r, name) => r.models.find((m) => m.name === name);
const box = (m) => inspect(m).stats.bbox;

test('vendored IFC engine matches the record in docs/vendor.md; only the files that are loaded are shipped', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  const dir = new URL('../../vendor/web-ifc/', import.meta.url);
  const sizes = {};
  for (const f of ['web-ifc-api.js', 'web-ifc.wasm', 'node/web-ifc-api-node.js', 'node/web-ifc-node.wasm', 'LICENSE.md']) {
    const buf = await readFile(new URL(f, dir));
    sizes[f] = buf.length;
    const sha = createHash('sha256').update(buf).digest('hex');
    assert.ok(doc.includes(sha), `${f}: SHA-256 ${sha} is recorded`);
    assert.ok(doc.includes(String(buf.length)), `${f}: size ${buf.length} is recorded`);
  }
  // browser runtime: ES-module API + single-threaded WebAssembly
  assert.equal(ifc.sizeBytes, sizes['web-ifc-api.js'] + sizes['web-ifc.wasm']);
  // no multi-threaded build, source maps, type definitions or bundles for other loaders
  assert.deepEqual((await readdir(dir)).sort(), ['LICENSE.md', 'node', 'web-ifc-api.js', 'web-ifc.wasm']);
  assert.deepEqual((await readdir(new URL('node/', dir))).sort(), ['package.json', 'web-ifc-api-node.js', 'web-ifc-node.wasm']);
  assert.deepEqual(ifc.library, { name: 'web-ifc', version: '0.0.78', licence: 'MPL-2.0' });
  assert.ok(doc.includes('web-ifc') && doc.includes('0.0.78') && doc.includes('MPL-2.0'));
  assert.equal(ifc.id, 'web-ifc');
  assert.deepEqual(ifc.formats, ['ifc']);
  assert.equal(ifc.available(), true);
  const first = await ifc.load();
  assert.equal(typeof first.api.OpenModel, 'function');
  assert.equal(await ifc.load(), first);
});

test('extruded beam in millimetres → closed surface in metres, z up: 0.2 × 0.4 × 6 m at (1, 2, 3)', async () => {
  const r = await ifc.parse({ name: 'models/frame.ifc', text: FRAME });
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.format, 'ifc'); assert.equal(r.meta.schema, 'IFC4');
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.equal(r.meta.engine, 'web-ifc 0.0.78');
  assert.deepEqual(r.meta.classes, { IfcColumn: 1, IfcBeam: 1 });
  const m = byName(r, 'IfcBeam');
  assert.equal(m.kind, 'surface');
  assert.ok(m.positions instanceof Float32Array && m.indices instanceof Uint32Array);
  assert.equal(m.ifcClass, 'IfcBeam'); assert.equal(m.source, 'frame');
  assert.equal(m.indices.length / 3, 12);
  assert.deepEqual(m.groups.map((g) => ({ name: g.name, start: g.start, count: g.count, expressId: g.expressId, globalId: g.globalId, tag: g.tag })),
    [{ name: 'B-101', start: 0, count: 36, expressId: 40, globalId: '3cUkl32yn9qRSPvBJVyWYp', tag: 'tag-1' }]);
  const i = inspect(m);
  assert.ok(i.ok && i.watertight && i.stats.outward);
  assert.equal(i.stats.components, 1);
  rel(i.stats.volume, 0.2 * 0.4 * 6, 1e-5);
  rel(i.stats.area, 2 * (0.2 * 0.4 + 0.2 * 6 + 0.4 * 6), 1e-5);
  // the extrusion runs along the placement z axis, which stays vertical
  const b = i.stats.bbox;
  [0.9, 1.8, 3].forEach((v, k) => close(b.min[k], v, 1e-5));
  [1.1, 2.2, 9].forEach((v, k) => close(b.max[k], v, 1e-5));
  // faces are tessellated separately; welding gives the eight corners back
  assert.equal(repair(m).model.positions.length / 3, 8);
  assert.deepEqual(r.meta.elements.find((e) => e.class === 'IfcBeam'), { class: 'IfcBeam', name: 'B-101', globalId: '3cUkl32yn9qRSPvBJVyWYp', expressId: 40, model: r.models.indexOf(m), group: 0, triangles: 12 });
  assert.equal(r.meta.counts.triangles, r.models.reduce((s, q) => s + q.indices.length / 3, 0));
});

test('circular column: the facet count follows circleSegments and the volume converges on π r² h', async () => {
  const exact = Math.PI * 0.5 ** 2 * 2;
  const vol = async (circleSegments) => {
    const r = await ifc.parse({ name: 'frame.ifc', text: FRAME }, { circleSegments, classes: ['ifccolumn'] });
    assert.deepEqual(Object.keys(r.meta.classes), ['IfcColumn'], 'the class filter is case-insensitive and drops the beam');
    assert.equal(r.meta.circleSegments, Math.min(128, Math.max(4, circleSegments)));
    const i = inspect(r.models[0]);
    assert.ok(i.watertight && i.stats.outward);
    close(i.stats.bbox.min[2], 0, 1e-6); close(i.stats.bbox.max[2], 2, 1e-6);
    return { v: i.stats.volume, tris: r.models[0].indices.length / 3 };
  };
  const coarse = await vol(12), fine = await vol(48), finest = await vol(500);
  // an inscribed polygon is always smaller than the circle
  assert.ok(coarse.v < fine.v && fine.v < finest.v && finest.v < exact);
  rel(coarse.v, exact, 0.07); rel(fine.v, exact, 0.005); rel(finest.v, exact, 0.001);
  assert.ok(coarse.tris < fine.tris && fine.tris < finest.tris);
});

test('wall with an opening: the void is subtracted (2.4 − 0.2 = 2.2 m³); opening volumes on request only', async () => {
  const r = await ifc.parse({ name: 'wall.ifc', text: WALL });
  assert.equal(r.meta.schema, 'IFC2X3');
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.deepEqual(r.meta.classes, { IfcWallStandardCase: 1 });
  assert.deepEqual(r.warnings, ['1 opening / space volume(s) were left out (pass openings = true to keep them).']);
  const i = inspect(r.models[0]);
  assert.ok(i.ok && i.watertight && i.stats.outward);
  rel(i.stats.volume, 4 * 0.2 * 3 - 1 * 0.2 * 1, 1e-6);
  // two wall faces with a hole, the outer rim and the four reveals of the opening
  rel(i.stats.area, 2 * (4 * 3 - 1) + 2 * 0.2 * (4 + 3) + 4 * 0.2 * 1, 1e-6);
  assert.equal(i.stats.genus, 1);
  const o = await ifc.parse({ name: 'wall.ifc', text: WALL }, { openings: true });
  assert.deepEqual(o.meta.classes, { IfcWallStandardCase: 1, IfcOpeningElement: 1 });
  assert.deepEqual(o.warnings, []);
  const op = inspect(byName(o, 'IfcOpeningElement'));
  rel(op.stats.volume, 1 * 0.6 * 1, 1e-6);
  close(op.stats.bbox.min[2], 1, 1e-6); close(op.stats.bbox.max[2], 2, 1e-6);
  assert.equal(byName(o, 'IfcOpeningElement').groups[0].name, 'O-1');
  // asking for the class by name works as well
  const c = await ifc.parse({ name: 'wall.ifc', text: WALL }, { classes: ['IfcOpeningElement'] });
  assert.deepEqual(c.meta.classes, { IfcOpeningElement: 1 });
});

test('length units: conversion-based feet are scaled by the file; opts.units overrides the declared unit', async () => {
  const r = await ifc.parse({ name: 'feet.ifc', bytes: enc(FEET) });
  assert.deepEqual(r.meta.units, { unit: 'ft', scale: 0.3048, basis: 'file' });
  rel(inspect(r.models[0]).stats.volume, 2 * 4 * 10 * 0.3048 ** 3, 1e-5);
  close(box(r.models[0]).max[2], 3.048, 1e-5);
  const m = await ifc.parse({ name: 'feet.ifc', bytes: enc(FEET) }, { units: 'm' });
  assert.deepEqual(m.meta.units, { unit: 'm', scale: 1, basis: 'option' });
  rel(inspect(m.models[0]).stats.volume, 80, 1e-6);
  const mm = await ifc.parse({ name: 'feet.ifc', bytes: enc(FEET) }, { units: 'mm' });
  rel(inspect(mm.models[0]).stats.volume, 80e-9, 1e-5);
  // no unit assignment at all: metres, and it is said
  const bare = await ifc.parse({ name: 'bare.ifc', text: FEET.replace(FOOT, '').replace(',(#8),#3);', ',(#8),$);') });
  assert.deepEqual(bare.meta.units, { unit: 'm', scale: 1, basis: 'assumed' });
  assert.ok(bare.warnings.some((w) => /No length unit was found in IfcUnitAssignment; metres assumed/.test(w)));
  rel(inspect(bare.models[0]).stats.volume, 80, 1e-6);
});

test('through the registry: fullGeometry routes an IFC file to the engine; without it the direct reader answers', async () => {
  const r = await importFile({ name: 'wall.ifc', bytes: enc(WALL) }, { fullGeometry: true });
  assert.equal(r.format, 'ifc-full');
  assert.equal(r.meta.via, 'adapter:web-ifc');
  assert.deepEqual(r.models.map((m) => [m.kind, m.name]), [['surface', 'IfcWallStandardCase']]);
  rel(inspect(r.models[0]).stats.volume, 2.2, 1e-6);
  await assert.rejects(importFile({ name: 'wall.ifc', bytes: enc(WALL) }), /no beams, columns, members, piles.*adapter "web-ifc"/);
  await assert.rejects(importFile({ name: 'wall.ifc', bytes: enc(WALL) }, { fullGeometry: true, adapters: { 'web-ifc': null } }), (e) => e instanceof GeometryFormatError && e.reason === 'adapter-missing');
});

test('content that is not an IFC STEP file, or has no surfaces, is reported precisely', async () => {
  assert.equal(ifc.sniff({ text: FRAME }), true);
  assert.equal(ifc.sniff({ text: FRAME.replace("FILE_SCHEMA(('IFC4'))", "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'))") }), false);
  await assert.rejects(() => ifc.parse({ name: 'notes.ifc', text: 'hello' }), /notes\.ifc: not an IFC STEP file \(expected "ISO-10303-21"/);
  await assert.rejects(() => ifc.parse({ name: 'pack.ifc', bytes: Uint8Array.from([0x50, 0x4b, 3, 4, 0, 0, 0, 0]) }), /ZIP archive: unpack the \.ifc inside an \.ifczip/);
  await assert.rejects(() => ifc.parse({ name: 'x.ifc', text: '<?xml version="1.0"?><ifcXML/>' }), /ifcXML is not read: export IFC-SPF/);
  const none = await ifc.parse({ name: 'frame.ifc', text: FRAME }, { classes: ['IfcPlate'] });
  assert.deepEqual(none.models, []);
  assert.deepEqual(none.warnings, ['No product of the requested classes (IfcPlate) has surface geometry.']);
  // products without any representation
  const empty = await ifc.parse({ name: 'empty.ifc', text: `${head('IFC4')}${SI('')}${CTX}#40=IFCBEAM('3cUkl32yn9qRSPvBJVyWYp',$,'B',$,$,#11,$,$,.BEAM.);\n${TAIL}` });
  assert.deepEqual(empty.models, []);
  assert.ok(empty.warnings.some((w) => /none of its products has surface geometry/.test(w)));
  // the engine keeps working after all of the above
  const again = await ifc.parse({ name: 'feet.ifc', text: FEET });
  assert.equal(again.models.length, 1);
});
