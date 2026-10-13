import test from 'node:test';
import assert from 'node:assert/strict';
import { exportStl, exportObj, exportPly, exportVtk, exportCsv, exportGeoJson, exportEsriAscii, exportJson, exportModels, EXPORTERS, gridToSurface, worldCoords } from '../../src/geometry/exporters.js';
import { importFile } from '../../src/geometry/registry.js';
import { parseStl, parseObj, parsePly } from '../../src/geometry/importers/surfaces.js';
import { parseVtk } from '../../src/geometry/importers/vtk.js';
import { parseTable } from '../../src/geometry/importers/tables.js';
import { parseGeoJson } from '../../src/geometry/importers/vector.js';
import { parseEsriAscii } from '../../src/geometry/importers/grids.js';
import { parseOaJson } from '../../src/geometry/importers/json.js';
import { surfaceModel, pointsModel } from '../../src/geometry/importers/common.js';
import { enc, close, closeArr, CUBE_V, CUBE_T, surfaceVolume, surfaceArea, checkElements, REF } from './fixtures/helpers.js';

const cube = (shift = [0, 0, 0]) => surfaceModel('cube', Float64Array.from(CUBE_V.flatMap((v) => v.map((c, k) => c + shift[k]))), Uint32Array.from(CUBE_T.flat()), [{ name: 'bottom', start: 0, count: 6 }, { name: 'rest', start: 6, count: 30 }]);
const GRID = { kind: 'grid', name: 'dem', nx: 3, ny: 2, x0: 500005, y0: 6000005, dx: 10, dy: 10, z: Float32Array.from([-1, -2, -3, -4, NaN, -6]), crs: 'EPSG:32631' };
const FRAME = { kind: 'frame', name: 'portal', nodes: Float64Array.from([0, 0, 0, 0, 0, 5, 4, 0, 5]), members: [{ a: 0, b: 1, section: { type: 'tube', D: 0.61, t: 0.019 }, material: 'S355', group: 'leg', id: 'M1' }, { a: 1, b: 2, section: { A: 0.01, Iy: 2e-4, Iz: 1e-5, J: 3e-6 }, group: 'beam', id: 'M2' }] };
const ROUTE = { kind: 'polyline', name: 'export cable', points: Float64Array.from([500000, 6000000, -30, 500100, 6000000, -31.5, 500200, 6000050, -33]), attrs: { kp: Float64Array.from([0, 100, 211.8]), burial: Float64Array.from([1, 1.5, 2]) } };

test('STL: binary export re-imports to the same mesh; UTM-sized models keep their precision through ORIGIN', () => {
  const m = cube();
  const b = exportStl(m);
  assert.equal(b.length, 84 + 50 * 12);
  const r = parseStl({ name: 'cube.stl', bytes: b });
  assert.equal(r.models[0].positions.length, 24);
  close(surfaceVolume(r.models[0]), 1, 1e-6);
  assert.equal(r.meta.units.basis, 'file', 'the header states UNITS=m');
  // outward normals are written
  const dv = new DataView(b.buffer);
  closeArr([dv.getFloat32(84, true), dv.getFloat32(88, true), dv.getFloat32(92, true)], [0, 0, -1], 1e-6);
  const far = cube([500000.25, 6000000.5, -40]);
  assert.ok(far.origin, 'stored relative to an origin');
  const rf = parseStl({ name: 'far.stl', bytes: exportStl(far) });
  assert.deepEqual(rf.models[0].origin, far.origin);
  close(rf.meta.bbox.min[0], 500000.25, 1e-6);
  close(rf.meta.bbox.max[1], 6000001.5, 1e-6);
  close(rf.meta.bbox.min[2], -40, 1e-6);
  close(surfaceVolume(rf.models[0]), 1, 1e-4);
  const g = parseStl({ name: 'g.stl', bytes: exportStl({ ...GRID, x0: 5, y0: 5 }) });
  assert.equal(g.meta.counts.triangles, 1, 'cells touching the no-data node are left out; one triangle of the other cell survives');
  assert.throws(() => exportStl(FRAME), /STL export does not handle a "frame" model/);
});

test('OBJ: groups, world coordinates, frames as line elements', () => {
  const far = cube([500000.25, 6000000.5, -40]);
  const text = exportObj(far);
  assert.match(text, /^# units: m\no cube\nv 500000\.25 6000000\.5 -40\n/);
  const r = parseObj({ name: 'far.obj', text });
  assert.deepEqual(r.models[0].groups.map((g) => [g.name, g.count]), [['bottom', 6], ['rest', 30]]);
  close(r.meta.bbox.max[0], 500001.25, 1e-6);
  close(surfaceVolume(r.models[0]), 1, 1e-4);
  const two = parseObj({ name: 'both.obj', text: exportObj([cube(), FRAME]) });
  assert.deepEqual(two.models.map((m) => m.kind), ['surface', 'polyline', 'polyline']);
  closeArr([...two.models[2].points], [0, 0, 5, 4, 0, 5]);
});

test('PLY: binary and ASCII with vertex colours', () => {
  const m = cube();
  m.attrs = { rgb: Uint8Array.from({ length: 24 }, (_, i) => i * 10) };
  for (const binary of [true, false]) {
    const out = exportPly(m, { binary });
    const r = parsePly(binary ? { name: 'c.ply', bytes: out } : { name: 'c.ply', text: out });
    close(surfaceVolume(r.models[0]), 1, 1e-6);
    assert.deepEqual([...r.models[0].attrs.rgb], [...m.attrs.rgb]);
    assert.equal(r.meta.units.basis, 'file');
  }
  const pts = pointsModel('cloud', Float64Array.from([500000.5, 6000000.5, 1, 500001.5, 6000002.5, 2]));
  const rp = parsePly({ name: 'p.ply', bytes: exportPly(pts) });
  assert.equal(rp.models[0].kind, 'points');
  closeArr([...worldCoords(rp.models[0])], [500000.5, 6000000.5, 1, 500001.5, 6000002.5, 2], 1e-6);
});

test('VTK: volume cells with nodal scalars and vectors survive a round trip', () => {
  const pts = [...REF.hex, ...REF.wedge.map((p) => [p[0] + 2, p[1], p[2]]), ...REF.tet.map((p) => [p[0] + 4, p[1], p[2]])];
  const vol = { kind: 'volume', name: 'block', positions: Float64Array.from(pts.flat()), cells: [{ type: 'hex8', conn: Uint32Array.from([0, 1, 2, 3, 4, 5, 6, 7]) }, { type: 'wedge6', conn: Uint32Array.from([8, 9, 10, 11, 12, 13]) }, { type: 'tet4', conn: Uint32Array.from([14, 15, 16, 17]) }] };
  checkElements(vol);
  const temp = Float64Array.from(pts, (p) => 10 + p[0]);
  const disp = Float64Array.from(pts.flatMap((p) => [0, 0, 0.001 * p[0]]));
  const text = exportVtk(vol, { scalars: { temperature: temp, 'disp vec': disp }, title: 'block\nresults' });
  assert.match(text, /^# vtk DataFile Version 3\.0\nblock results\nASCII\nDATASET UNSTRUCTURED_GRID\n/);
  const r = parseVtk({ name: 'block.vtk', text });
  const v = r.models[0];
  assert.deepEqual(r.meta.elementTypes, { hex8: 1, wedge6: 1, tet4: 1 });
  checkElements(v);
  for (let i = 0; i < v.nodeIds.length; i++) { close(v.nodeData.temperature[i], temp[v.nodeIds[i]], 1e-5); close(v.nodeData.disp_vec[3 * i + 2], disp[3 * v.nodeIds[i] + 2], 1e-9); }
  assert.throws(() => exportVtk(vol, { scalars: { bad: [1, 2, 3] } }), /nodal field "bad" has 3 values for 18 nodes/);
});

test('VTK: surface, frame, polyline and grid forms', () => {
  const s = parseVtk({ name: 's.vtk', text: exportVtk(cube()) });
  close(surfaceVolume(s.models[0]), 1, 1e-6);
  const f = parseVtk({ name: 'f.vtk', text: exportVtk(FRAME) });
  assert.equal(f.models[0].kind, 'frame');
  assert.equal(f.models[0].members.length, 2);
  const l = parseVtk({ name: 'l.vtk', text: exportVtk(ROUTE) });
  assert.equal(l.models[0].members.length, 2);
  close(l.meta.bbox.max[1], 6000050, 1e-9);
  const g = parseVtk({ name: 'g.vtk', text: exportVtk(GRID) });
  assert.equal(g.models[0].kind, 'grid');
  assert.deepEqual([g.models[0].nx, g.models[0].ny, g.models[0].x0, g.models[0].y0, g.models[0].dx], [3, 2, 500005, 6000005, 10]);
  assert.deepEqual([...g.models[0].z].map((v) => (Number.isNaN(v) ? null : v)), [-1, -2, -3, -4, null, -6]);
  assert.throws(() => exportVtk({ kind: 'section', loops: [] }), /VTK export does not handle a "section" model/);
});

test('CSV: frame, route, points and grid tables re-import to the same models', () => {
  const f = parseTable({ name: 'portal.csv', text: exportCsv(FRAME) });
  const m = f.models[0];
  assert.equal(m.kind, 'frame');
  closeArr([...m.nodes], [...FRAME.nodes]);
  assert.deepEqual(m.members[0], { a: 0, b: 1, section: { type: 'tube', D: 0.61, t: 0.019 }, material: 'S355', group: 'leg', id: 'M1' });
  assert.deepEqual(m.members[1].section, { A: 0.01, Iy: 2e-4, Iz: 1e-5, J: 3e-6 });
  const r = parseTable({ name: 'export_cable.csv', text: exportCsv(ROUTE) });
  assert.equal(r.models[0].kind, 'polyline');
  closeArr([...r.models[0].points], [...ROUTE.points]);
  closeArr([...r.models[0].attrs.kp], [...ROUTE.attrs.kp]);
  closeArr([...r.models[0].attrs.burial], [1, 1.5, 2]);
  const pts = pointsModel('wtg', Float64Array.from([0, 0, 0, 800, 0, 0, 0, 900, 1]), { names: ['A, 1', 'A"2"', 'B1'], rating: Float64Array.from([15, 15, NaN]) });
  const text = exportCsv(pts);
  assert.match(text, /"A, 1"/);
  const p = parseTable({ name: 'wtg.csv', text });
  assert.deepEqual(p.models[0].attrs.names, ['A, 1', 'A"2"', 'B1']);
  assert.ok(Number.isNaN(p.models[0].attrs.rating[2]));
  const g = parseTable({ name: 'dem.csv', text: exportCsv({ ...GRID, z: Float32Array.from([-1, -2, -3, -4, -5, -6]) }) });
  assert.equal(g.models[0].kind, 'grid');
  assert.deepEqual([g.models[0].nx, g.models[0].ny, g.models[0].dx, ...g.models[0].z], [3, 2, 10, -1, -2, -3, -4, -5, -6]);
  assert.match(exportCsv({ kind: 'section', name: 's', loops: [Float64Array.from([0, 0, 1, 0, 0, 1])] }), /^loop,x \(m\),y \(m\)\n1,0,0\n1,1,0\n1,0,1\n$/);
  assert.throws(() => exportCsv(cube()), /CSV export does not handle a "surface" model/);
});

test('GeoJSON: projected output with a crs member, and WGS 84 output through the tangent-plane origin', () => {
  const text = exportGeoJson([ROUTE, pointsModel('wtg', Float64Array.from([500010, 6000020, 0]), { names: ['A1'] })], { crs: 'EPSG:32631' });
  const doc = JSON.parse(text);
  assert.equal(doc.crs.properties.name, 'urn:ogc:def:crs:EPSG::32631');
  assert.deepEqual(doc.features[0].geometry.coordinates[1], [500100, 6000000, -31.5]);
  const r = parseGeoJson({ name: 'r.geojson', text });
  closeArr([...r.models.find((m) => m.kind === 'polyline').points], [...ROUTE.points]);
  assert.deepEqual(r.models.find((m) => m.kind === 'points').attrs.names, ['A1']);
  // local metres ↔ longitude/latitude about an origin
  const origin = { lon: 3, lat: 54 };
  const local = { kind: 'polyline', name: 'cable', points: Float64Array.from([0, 0, -30, 1500, 800, -32]) };
  const sec = { kind: 'section', name: 'lease', loops: [Float64Array.from([0, 0, 1000, 0, 1000, 500, 0, 500]), Float64Array.from([100, 100, 100, 200, 200, 200])] };
  const wgs = exportGeoJson([local, sec, FRAME], { origin });
  const d2 = JSON.parse(wgs);
  assert.equal(d2.crs, undefined);
  close(d2.features[0].geometry.coordinates[0][0], 3, 1e-12);
  assert.equal(d2.features[1].geometry.type, 'Polygon');
  assert.equal(d2.features[1].geometry.coordinates[0].length, 5, 'rings are closed');
  const back = parseGeoJson({ name: 'b.geojson', text: wgs }, { origin });
  closeArr([...back.models.find((m) => m.kind === 'polyline' && m.name === 'cable').points], [...local.points], 1e-6);
  closeArr([...back.models.find((m) => m.kind === 'section').loops[0]], [...sec.loops[0]], 1e-6);
  assert.equal(back.models.filter((m) => m.kind === 'polyline').length, 3, 'two frame members come back as lines');
  assert.throws(() => exportGeoJson(GRID), /GeoJSON export does not handle a "grid" model/);
});

test('ESRI ASCII grid: no-data, cell registration and unequal spacing round-trip', () => {
  const text = exportEsriAscii(GRID);
  assert.match(text, /^ncols 3\nnrows 2\nxllcorner 500000\nyllcorner 6000000\ncellsize 10\nNODATA_value -9999\n-4 -9999 -6\n-1 -2 -3\n$/);
  const r = parseEsriAscii({ name: 'dem.asc', text });
  const g = r.models[0];
  assert.deepEqual([g.nx, g.ny, g.x0, g.y0, g.dx, g.dy], [3, 2, 500005, 6000005, 10, 10]);
  assert.deepEqual([...g.z].map((v) => (Number.isNaN(v) ? null : v)), [-1, -2, -3, -4, null, -6]);
  const rect = parseEsriAscii({ name: 'r.asc', text: exportEsriAscii({ ...GRID, dy: 20 }) });
  assert.deepEqual([rect.models[0].dx, rect.models[0].dy, rect.models[0].y0], [10, 20, 6000005]);
  assert.throws(() => exportEsriAscii(FRAME), /ESRI ASCII grid export does not handle a "frame" model/);
});

test('application JSON: every model kind round-trips losslessly, typed arrays and NaN included', async () => {
  const far = cube([500000.25, 6000000.5, -40]);
  const vol = { kind: 'volume', name: 'v', positions: Float64Array.from(REF.tet.flat()), cells: [{ type: 'tet4', conn: Uint32Array.from([0, 1, 2, 3]), tags: Int32Array.from([-7]) }], sets: { nodes: { base: Uint32Array.from([0, 1, 2]) } } };
  const sec = { kind: 'section', name: 's', loops: [Float64Array.from([0, 0, 1, 0, 0, 1])], layers: [{ name: 'steel', loops: [0] }] };
  const pts = pointsModel('p', Float64Array.from([1, 2, 3]), { names: ['a'], rgb: Uint8Array.from([1, 2, 3]), flag: Float64Array.from([NaN]) });
  const models = [far, vol, FRAME, ROUTE, GRID, pts, sec];
  const text = exportJson(models, { crs: 'EPSG:32631', source: 'unit test' });
  const doc = JSON.parse(text);
  assert.equal(doc.format, 'offshore-analyses-geometry');
  assert.equal(doc.version, 1);
  const r = parseOaJson({ name: 'all.json', text });
  assert.deepEqual(r.models, models);
  assert.ok(r.models[4].z instanceof Float32Array && Number.isNaN(r.models[4].z[4]));
  assert.ok(r.models[1].cells[0].tags instanceof Int32Array);
  assert.equal(r.meta.crs, 'EPSG:32631');
  assert.deepEqual(r.meta.counts.byKind, { surface: 1, volume: 1, frame: 1, polyline: 1, grid: 1, points: 1, section: 1 });
  // through the registry, by content
  const viaRegistry = await importFile({ name: 'whatever.json', bytes: enc(text) });
  assert.equal(viaRegistry.format, 'oa-json');
  assert.deepEqual(viaRegistry.models[2], FRAME);
  assert.throws(() => parseOaJson({ name: 'x.json', text: '{"format":"offshore-analyses-geometry","version":99,"models":[{}]}' }), /version 99 is newer than this build reads/);
  assert.throws(() => parseOaJson({ name: 'x.json', text: '{"format":"offshore-analyses-geometry","version":1,"models":[{"kind":"blob"}]}' }), /unknown kind "blob"/);
});

test('exportModels: file names, kind checks, format list', () => {
  assert.deepEqual(EXPORTERS.map((e) => e.id), ['stl', 'obj', 'ply', 'vtk', 'csv', 'geojson', 'esri-ascii', 'hst', 'json']);
  const e = exportModels('stl', cube(), { fileName: 'my hull/v2.step' });
  assert.equal(e.fileName, 'my hull_v2.stl');
  assert.equal(e.mime, 'model/stl');
  assert.ok(e.data instanceof Uint8Array);
  assert.equal(exportModels('esri-ascii', GRID).fileName, 'dem.asc');
  assert.equal(typeof exportModels('json', [cube(), FRAME]).data, 'string');
  assert.throws(() => exportModels('stl', FRAME), /STL \(binary\) export does not handle a "frame" model \(supported kinds: surface, grid\)/);
  assert.throws(() => exportModels('stl', [cube(), cube()]), /holds one model per file; 2 were given/);
  assert.throws(() => exportModels('dwg', cube()), /unknown format "dwg" \(available: stl, obj, ply, vtk, csv, geojson, esri-ascii, hst, json\)/);
  close(surfaceArea(gridToSurface({ ...GRID, z: new Float32Array(6) })), 200, 1e-9);
});
