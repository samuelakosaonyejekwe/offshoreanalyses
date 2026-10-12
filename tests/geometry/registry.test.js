import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { FORMATS, importFile, formatFor, acceptString, adapterStatus, GeometryFormatError } from '../../src/geometry/registry.js';
import { detectFormat, detectFileSet, EXTENSION_MAP } from '../../src/geometry/detect.js';
import { makeZip } from '../../src/geometry/zip.js';
import { enc, close, binaryStl, CUBE_T, CUBE_V, CUBE_Q, concat } from './fixtures/helpers.js';

const root = fileURLToPath(new URL('../../', import.meta.url));
const KINDS = ['surface', 'volume', 'frame', 'polyline', 'grid', 'points', 'section'];
const CATEGORIES = ['cad', 'structural', 'mesh', 'centreline', 'geospatial', 'pointcloud', 'panel', 'section', 'airfoil', 'table'];

test('FORMATS: every entry is complete, consistent and truthful about how it is read', async () => {
  const ids = new Set();
  for (const f of FORMATS) {
    assert.ok(!ids.has(f.id), `duplicate id ${f.id}`);
    ids.add(f.id);
    assert.equal(typeof f.name, 'string');
    assert.ok(Array.isArray(f.extensions) && f.extensions.every((e) => /^[a-z0-9_]+$/.test(e)), `${f.id} extensions`);
    assert.ok(f.kinds.every((k) => KINDS.includes(k)), `${f.id} kinds`);
    assert.ok(CATEGORIES.includes(f.category), `${f.id} category ${f.category}`);
    assert.ok(Array.isArray(f.suites) && f.suites.every((s) => Number.isInteger(s) && s >= 1 && s <= 16), `${f.id} suites`);
    assert.ok(typeof f.notes === 'string' && f.notes.length > 5, `${f.id} notes`);
    if (f.status === 'native') {
      assert.match(f.via, /^importers\/[a-z0-9]+\.js$/, f.id);
      assert.ok(existsSync(`${root}src/geometry/${f.via}`), `${f.id}: ${f.via} exists`);
      const mod = await f.load();
      assert.equal(typeof mod.parse, 'function', `${f.id} load() gives parse`);
      assert.ok(f.kinds.length > 0, `${f.id} declares kinds`);
    } else if (f.status === 'adapter') {
      assert.match(f.via, /^adapter:(occt|web-ifc|hdf5|laz|sqlite)$/, f.id);
      assert.ok(typeof f.missing === 'string' && f.missing.length > 20, `${f.id} says what to do without the adapter`);
    } else if (f.status === 'convert') {
      assert.equal(f.via, null);
      assert.equal(f.load, null);
      assert.deepEqual(f.kinds, []);
      assert.ok(typeof f.convert === 'string' && /[Ee]xport|[Ss]ave|[Cc]onvert/.test(f.convert), `${f.id} says what to export`);
    } else assert.fail(`${f.id}: unknown status ${f.status}`);
  }
  // every id the detector can return has an entry
  for (const id of new Set(Object.values(EXTENSION_MAP))) assert.ok(ids.has(id), `detector id ${id} is registered`);
  for (const id of ['panel-csv', 'airfoil', 'xyzrgb', 'openfoam', 'oa-json', 'landxml', 'exodus-hdf5', 'dxf-binary', 'fluent', 'cis2']) assert.ok(ids.has(id), id);
});

test('FORMATS covers every format family named in the suite specifications', () => {
  const spec = [];
  for (const f of readdirSync(`${root}src/spec`)) if (/^suite-\d+\.json$/.test(f)) for (const g of JSON.parse(readFileSync(`${root}src/spec/${f}`, 'utf8')).geometry || []) spec.push(g.text);
  const text = spec.join(' ');
  const need = { STEP: 'step', IGES: 'iges', Parasolid: 'parasolid', ACIS: 'acis', DXF: 'dxf', DWG: 'dwg', IFC: 'ifc', SDNF: 'sdnf', 'CIS/2': 'cis2', Nastran: 'nastran', Abaqus: 'abaqus', ANSYS: 'ansys-cdb', Gmsh: 'msh', MED: 'med', UNV: 'unv', VTK: 'vtk', Exodus: 'exodus', CGNS: 'cgns', SU2: 'su2', OpenFOAM: 'openfoam', STL: 'stl', OBJ: 'obj', PLY: 'ply', '3MF': '3mf', GLB: 'glb', LAS: 'las', LAZ: 'laz', E57: 'e57', GeoTIFF: 'geotiff', 'ESRI ASCII': 'esri-ascii', LandXML: 'landxml', GML: 'gml', GeoJSON: 'geojson', GeoPackage: 'geopackage', Shapefile: 'shapefile', NetCDF: 'netcdf', HDF5: 'hdf5', CSV: 'table', JSON: 'table-json', SVG: 'svg', GDF: 'gdf', NEMOH: 'nemoh', Hydrostar: 'hydrostar', Fluent: 'fluent' };
  for (const [word, id] of Object.entries(need)) {
    assert.ok(text.includes(word), `spec mentions ${word}`);
    const f = FORMATS.find((x) => x.id === id);
    assert.ok(f, `${word} → ${id} registered`);
    assert.ok(f.suites.length > 0, `${id} lists the suites that ask for it`);
  }
});

test('detection: magic bytes, content sniffing and extensions', () => {
  const d = (name, content) => { const r = detectFormat(typeof content === 'string' ? { name, bytes: enc(content) } : { name, bytes: content }); return r ? r.id : null; };
  assert.equal(d('a.bin', binaryStl(CUBE_T, CUBE_V)), 'stl');
  assert.equal(d('a.txt', 'solid x\n facet normal 0 0 1\n endfacet\nendsolid x\n'), 'stl');
  assert.equal(d('a', 'ply\nformat ascii 1.0\n'), 'ply');
  assert.equal(d('a.dat', '# vtk DataFile Version 3.0\n'), 'vtk');
  assert.equal(d('a.dat', '$MeshFormat\n2.2 0 8\n'), 'msh');
  assert.equal(d('a.msh', '(0 "Fluent mesh")\n(2 3)\n'), 'fluent');
  assert.equal(d('m.stp', 'ISO-10303-21;\nHEADER;\nFILE_SCHEMA((\'IFC4\'));\nENDSEC;\n'), 'ifc');
  assert.equal(d('m.stp', 'ISO-10303-21;\nHEADER;\nFILE_SCHEMA((\'STRUCTURAL_FRAME_SCHEMA\'));\nENDSEC;\n'), 'cis2');
  assert.equal(d('m.ifc', 'ISO-10303-21;\nHEADER;\nFILE_SCHEMA((\'AUTOMOTIVE_DESIGN\'));\nENDSEC;\n'), 'step');
  assert.equal(d('s.xml', '<?xml version="1.0"?>\n<!-- c -->\n<LandXML version="1.2"></LandXML>'), 'landxml');
  assert.equal(d('s.xml', '<gml:FeatureCollection xmlns:gml="http://www.opengis.net/gml"/>'), 'gml');
  assert.equal(d('s.xml', '<Model><Nodes/></Model>'), 'xml-model');
  assert.equal(d('s.xml', '<VTKFile type="UnstructuredGrid"/>'), 'vtu');
  assert.equal(d('s.json', '{"type":"FeatureCollection","features":[]}'), 'geojson');
  assert.equal(d('s.json', '{"format":"offshore-analyses-geometry","version":1,"models":[]}'), 'oa-json');
  assert.equal(d('s.json', '{"asset":{"version":"2.0"},"meshes":[]}'), 'gltf');
  assert.equal(d('s.json', '{"nodes":[],"members":[]}'), 'table-json');
  assert.equal(d('g.txt', 'ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\n'), 'esri-ascii');
  assert.equal(d('pts.asc', '1.0 2.0 3.0\n4.0 5.0 6.0\n7 8 9\n'), 'table');
  assert.equal(d('d.inp', '*Heading\n*Node\n1, 0, 0, 0\n'), 'abaqus');
  assert.equal(d('d.dat', 'BEGIN BULK\nGRID,1,,0.,0.,0.\n'), 'nastran');
  assert.equal(d('d.dat', '    -1\n  2411\n         1         1         1        11\n'), 'unv');
  assert.equal(d('d.txt', 'Packet 00\n"a"\nPacket 10\n'), 'sdnf');
  assert.equal(d('d.txt', 'NDIME= 2\nNELEM= 1\n'), 'su2');
  assert.equal(d('d.txt', 'LINESTRING Z (0 0 0, 1 1 1)'), 'wkt');
  assert.equal(d('d.txt', '  0\nSECTION\n  2\nHEADER\n'), 'dxf');
  assert.equal(d('d.dxf', 'AutoCAD Binary DXF\r\n\x1a\0'), 'dxf-binary');
  assert.equal(d('foil.dat', `NACA 0012\n${Array.from({ length: 12 }, (_, i) => `${(1 - i / 11).toFixed(4)} ${(0.05 * Math.sin((Math.PI * i) / 11)).toFixed(4)}`).join('\n')}\n`), 'airfoil');
  assert.equal(d('mesh.dat', '2 1\n1 0.0 0.0 -1.0\n2 1.0 0.0 -1.0\n0 0 0 0\n1 2 2 1\n0 0 0 0\n'), 'nemoh');
  assert.equal(d('c.txt', '0 0 0 255 0 0\n1 0 0 0 255 0\n1 1 0 0 0 255\n'), 'xyzrgb');
  assert.equal(d('p.csv', '0 0 0 1 0 0 1 1 0 0 1 0\n'), 'panel-csv');
  assert.equal(d('x.nc', concat(Uint8Array.from([0x43, 0x44, 0x46, 2]), new Uint8Array(40))), 'netcdf');
  assert.equal(d('x.nc', concat(Uint8Array.from([0x89, 0x48, 0x44, 0x46, 13, 10, 26, 10]), new Uint8Array(40))), 'netcdf4');
  assert.equal(d('x.e', concat(Uint8Array.from([0x89, 0x48, 0x44, 0x46, 13, 10, 26, 10]), new Uint8Array(40))), 'exodus-hdf5');
  assert.equal(d('x.tif', concat(Uint8Array.from([0x4d, 0x4d, 0, 42]), new Uint8Array(40))), 'geotiff');
  const las = new Uint8Array(230); las.set(enc('LASF'));
  assert.equal(d('x.bin', las), 'las');
  las[104] = 0x83;
  assert.equal(d('x.las', las), 'laz');
  assert.equal(d('x.x_t', '**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**\n'), 'parasolid');
  assert.equal(d('x.sldprt', new Uint8Array(64)), 'native-cad');
  assert.equal(d('x.unknown', 'just some words\n'), null);
  assert.equal(detectFileSet(['case/constant/polyMesh/points', 'case/constant/polyMesh/faces', 'case/constant/polyMesh/owner']).id, 'openfoam');
  assert.equal(detectFileSet(['a/roads.shp', 'a/roads.dbf', 'a/roads.prj']).main, 'a/roads.shp');
});

test('formatFor and acceptString', () => {
  assert.equal(formatFor('hull.stl').id, 'stl');
  assert.equal(formatFor('hull.STL', binaryStl(CUBE_T, CUBE_V)).status, 'native');
  assert.equal(formatFor('part.x_t').status, 'convert');
  assert.equal(formatFor('scan.laz').via, 'adapter:laz');
  assert.equal(formatFor('model.stp.gz', Uint8Array.from([0x1f, 0x8b, 8, 0])).id, 'step');
  assert.equal(formatFor('set.zip', makeZip([{ name: 'r.shp', data: new Uint8Array(4) }, { name: 'r.dbf', data: new Uint8Array(4) }])).id, 'shapefile');
  assert.equal(formatFor('notes.docx'), null);
  const acc = acceptString().split(',');
  for (const e of ['.stl', '.step', '.ifc', '.tif', '.shp', '.dbf', '.zip', '.gz', '.x_t', '.dwg', '.e57', '.3mf']) assert.ok(acc.includes(e), e);
  assert.ok(acc.every((e) => /^\.[a-z0-9_]+$/.test(e)));
});

test('importFile: detection, result envelope, kinds within the registered kinds', async () => {
  const r = await importFile({ name: 'cube.stl', bytes: binaryStl(CUBE_T, CUBE_V) });
  assert.equal(r.format, 'stl');
  assert.equal(r.models[0].kind, 'surface');
  assert.equal(r.meta.detected.id, 'stl');
  assert.equal(r.meta.fileName, 'cube.stl');
  assert.equal(r.meta.counts.triangles, 12);
  const cases = [
    ['a.obj', 'v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n'],
    ['a.xyz', '0 0 1\n1 0 2\n0 1 3\n1 1 4\n'],
    ['route.csv', 'kp,x,y,z\n0,0,0,-10\n1,1000,0,-12\n'],
    ['a.geojson', '{"type":"Feature","properties":{},"geometry":{"type":"LineString","coordinates":[[3,54],[3.1,54]]}}'],
    ['a.wkt', 'POINT (1 2)'],
    ['a.asc', 'ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 1\n1 2\n3 4\n'],
    ['a.su2', 'NDIME= 2\nNELEM= 1\n5 0 1 2\nNPOIN= 3\n0 0\n1 0\n0 1\n'],
    ['a.svg', '<svg width="10mm" height="10mm" viewBox="0 0 10 10"><rect width="5" height="5"/></svg>'],
    ['a.gdf', 'h\n1 9.80665\n0 0\n1\n0 0 -1 1 0 -1 1 1 -1 0 1 -1\n'],
    ['a.json', '{"nodes":[{"id":1,"x":0,"y":0,"z":0},{"id":2,"x":0,"y":0,"z":5}],"members":[{"a":1,"b":2,"D":0.6,"t":0.02}]}'],
  ];
  for (const [name, text] of cases) {
    const x = await importFile({ name, bytes: enc(text) });
    const f = FORMATS.find((q) => q.id === x.format);
    assert.ok(f && f.status === 'native', `${name} → ${x.format}`);
    for (const m of x.models) assert.ok(f.kinds.includes(m.kind), `${name}: kind ${m.kind} is declared for ${f.id}`);
    assert.ok(x.meta.bbox && x.meta.counts && x.meta.elementTypes && x.meta.skipped, `${name} report`);
  }
  await assert.rejects(importFile({ name: 'x' }), /expected \{ name, bytes \}/);
  await assert.rejects(importFile({ name: 'x.stl', bytes: binaryStl(CUBE_T, CUBE_V) }, { format: 'nope' }), /unknown format id "nope"/);
});

test('importFile: opts.format overrides detection; options reach the importer', async () => {
  const text = '0 0 0\n10 0 0\n20 5 0\n';
  assert.equal((await importFile({ name: 'data.txt', text })).models[0].kind, 'points');
  const r = await importFile({ name: 'data.txt', text }, { format: 'table', schema: 'route', units: 'mm' });
  assert.equal(r.models[0].kind, 'polyline');
  close(r.meta.bbox.max[0], 0.02, 1e-12);
  assert.equal(r.meta.detected.reason, 'chosen by the caller');
  // the schema-specific table entries of the matrix
  const route = await importFile({ name: 'data.txt', text }, { format: 'route-table' });
  assert.equal(route.format, 'route-table');
  assert.equal(route.models[0].kind, 'polyline');
  const hull = await importFile({ name: 'h.csv', text: 'x,z,y\n0,0,0\n0,2,3\n10,0,0\n10,2,3\n' }, { format: 'hull-offsets', mapping: { station: 'x', waterline: 'z', halfBreadth: 'y' } });
  assert.equal(hull.models[0].kind, 'surface');
  const frame = await importFile({ name: 'f.csv', text: 'x1,y1,z1,x2,y2,z2\n0,0,0,0,0,5\n' }, { format: 'frame-table' });
  assert.equal(frame.models[0].kind, 'frame');
  const blade = await importFile({ name: 'b.csv', text: 'a,b\n1,3\n20,2\n' }, { format: 'blade-table', mapping: { r: 'a', chord: 'b' } });
  assert.deepEqual([...blade.models[0].attrs.chord], [3, 2]);
});

test('importFile: gzip and zip containers, multi-file sets', async () => {
  const stl = binaryStl(CUBE_T, CUBE_V);
  const gz = await importFile({ name: 'cube.stl.gz', bytes: new Uint8Array(gzipSync(Buffer.from(stl))) });
  assert.equal(gz.format, 'stl');
  assert.deepEqual(gz.meta.detected.containers, ['gzip']);
  assert.equal(gz.meta.fileName, 'cube.stl');
  const zip = await importFile({ name: 'bundle.zip', bytes: makeZip([{ name: '__MACOSX/._cube.stl', data: new Uint8Array(3) }, { name: 'models/cube.stl', data: stl }, { name: 'readme.md', data: enc('notes') }]) });
  assert.equal(zip.format, 'stl');
  assert.deepEqual(zip.meta.detected.containers, ['zip']);
  // glTF with a side-car buffer, zipped
  const pos = new Float32Array(CUBE_V.flat()), idx = new Uint16Array(CUBE_T.flat());
  const bin = concat(new Uint8Array(pos.buffer), new Uint8Array(idx.buffer));
  const gltf = { asset: { version: '2.0' }, scene: 0, scenes: [{ nodes: [0] }], nodes: [{ mesh: 0 }], meshes: [{ primitives: [{ attributes: { POSITION: 0 }, indices: 1 }] }], accessors: [{ bufferView: 0, componentType: 5126, count: 8, type: 'VEC3' }, { bufferView: 1, componentType: 5123, count: 36, type: 'SCALAR' }], bufferViews: [{ buffer: 0, byteLength: 96 }, { buffer: 0, byteOffset: 96, byteLength: 72 }], buffers: [{ uri: 'cube.bin', byteLength: 168 }] };
  const g = await importFile({ name: 'scene.zip', bytes: makeZip([{ name: 'scene/cube.gltf', data: enc(JSON.stringify(gltf)) }, { name: 'scene/cube.bin', data: bin }]) });
  assert.equal(g.format, 'gltf');
  assert.equal(g.meta.counts.triangles, 12);
  // OpenFOAM case as a file map
  const head = 'FoamFile\n{\n format ascii;\n}\n';
  const foam = await importFile({ name: 'case', files: new Map([
    ['constant/polyMesh/points', enc(`${head}8\n(\n${CUBE_V.map((p) => `(${p.join(' ')})`).join('\n')}\n)\n`)],
    ['constant/polyMesh/faces', enc(`${head}6\n(\n${CUBE_Q.map((q) => `4(${q.join(' ')})`).join('\n')}\n)\n`)],
    ['constant/polyMesh/owner', enc(`${head}6\n(\n0\n0\n0\n0\n0\n0\n)\n`)],
    ['constant/polyMesh/boundary', enc(`${head}1\n(\n walls\n {\n  type wall;\n  nFaces 6;\n  startFace 0;\n }\n)\n`)],
  ]) });
  assert.equal(foam.format, 'openfoam');
  assert.equal(foam.models[0].groups[0].name, 'walls');
  await assert.rejects(importFile({ name: 'empty.zip', bytes: makeZip([{ name: 'notes.md', data: enc('x') }]) }), (e) => e instanceof GeometryFormatError && e.reason === 'unrecognised' && /Files supplied: notes\.md/.test(e.message));
});

test('importFile: convert-only formats say exactly what to export', async () => {
  const cases = [
    ['part.x_t', '**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz**\n', /Parasolid files are not read\. Export STEP \(AP214 or AP242\)/],
    ['part.sat', '700 0 1 0\n', /ACIS SAT \/ SAB files are not read\. Export STEP/],
    ['plan.dwg', 'AC1032\0\0\0\0\0', /DWG drawings are not read\..*save the drawing as ASCII DXF/],
    ['plan.dxf', 'AutoCAD Binary DXF\r\n\x1a\0', /Binary DXF is not read\. Save the drawing as ASCII DXF/],
    ['scan.e57', 'ASTM-E57\0\0\0\0', /E57 scans are not read\. Export the cloud as LAS/],
    ['part.sldprt', '\0\0\0\0binary', /Native CAD files are not read\..*export STEP/],
    ['mesh.hst', 'HYDROSTAR\n', /Export the panel mesh as WAMIT GDF or NEMOH/],
    ['flow.cas', '(0 "fluent")\n(2 3)\n', /Convert the mesh to VTK, SU2 or Gmsh/],
  ];
  for (const [name, text, re] of cases) await assert.rejects(importFile({ name, bytes: enc(text) }), (e) => e instanceof GeometryFormatError && e.reason === 'convert' && re.test(e.message), name);
});

const P21_HEAD = (schema) => `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('a','',(''),(''),'','','');\nFILE_SCHEMA(('${schema}'));\nENDSEC;\nDATA;\n`;
const BREP = `${P21_HEAD('AUTOMOTIVE_DESIGN')}#10=CARTESIAN_POINT('',(0.,0.,0.));\n#14=AXIS2_PLACEMENT_3D('',#10,$,$);\n#16=EDGE_LOOP('',());\n#17=FACE_BOUND('',#16,.T.);\n#18=PLANE('',#14);\n#19=ADVANCED_FACE('',(#17),#18,.T.);\n#20=CLOSED_SHELL('',(#19));\n#21=MANIFOLD_SOLID_BREP('s',#20);\nENDSEC;\nEND-ISO-10303-21;\n`;

test('importFile: adapter-backed formats — missing adapter is explained, supplied adapter is used', async () => {
  const hdf = concat(Uint8Array.from([0x89, 0x48, 0x44, 0x46, 13, 10, 26, 10]), new Uint8Array(64));
  // `null` declares an adapter unavailable, whatever this build happens to ship
  const none = { adapters: { occt: null, 'web-ifc': null, hdf5: null, laz: null, sqlite: null } };
  const missing = (re) => (e) => e instanceof GeometryFormatError && e.reason === 'adapter-missing' && re.test(e.message);
  await assert.rejects(importFile({ name: 'bathy.nc', bytes: hdf }, none), missing(/NetCDF-4 needs the "hdf5" reader.*not available in this build.*nccopy -k classic/));
  await assert.rejects(importFile({ name: 'mesh.med', bytes: hdf }, none), missing(/MED \(SALOME\) needs the "hdf5" reader.*Export the mesh from SALOME as UNV/));
  await assert.rejects(importFile({ name: 'mesh.cgns', bytes: hdf }, none), missing(/Convert the mesh to VTK, SU2 or Gmsh/));
  const laz = new Uint8Array(230); laz.set(enc('LASF')); laz[104] = 0x80;
  await assert.rejects(importFile({ name: 'scan.laz', bytes: laz }, none), missing(/needs the "laz" reader.*laszip -i in\.laz -o out\.las/));
  await assert.rejects(importFile({ name: 'site.gpkg', bytes: concat(enc('SQLite format 3\0'), new Uint8Array(100)) }, none), missing(/needs the "sqlite" reader.*ogr2ogr -f GeoJSON/));
  await assert.rejects(importFile({ name: 'shape.brep', bytes: enc('DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\n') }, none), missing(/OpenCascade BREP needs the "occt" reader/));
  await assert.rejects(importFile({ name: 'solid.stp', bytes: enc(BREP) }, none), missing(/only exact geometry \(1 MANIFOLD_SOLID_BREP, 1 ADVANCED_FACE, 1 CLOSED_SHELL\).*needs the "occt" reader.*tessellated STEP/));
  const ifc = `${P21_HEAD('IFC4')}#5=IFCSIUNIT(*,.LENGTHUNIT.,$,.METRE.);\n#7=IFCUNITASSIGNMENT((#5));\n#90=IFCSLAB('2gRXFgjRn2HPE$YoDLX3FV',$,'S-1',$,$,$,$,$,.FLOOR.);\nENDSEC;\nEND-ISO-10303-21;\n`;
  await assert.rejects(importFile({ name: 'slab.ifc', bytes: enc(ifc) }, { ...none, fullGeometry: true }), missing(/full geometry.*needs the "web-ifc" reader/));
  // an adapter module that reports itself as not built in counts as missing too
  // a supplied adapter receives the input and its result is passed through
  let seen = null;
  const fake = { parse: async (input, o) => { seen = { name: input.name, n: input.bytes.length, units: o.units, format: o.format }; return { models: [{ kind: 'surface', name: 'solid', positions: new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), indices: new Uint32Array([0, 1, 2]) }], warnings: ['from adapter'], meta: { tolerance: 0.1 } }; } };
  const r = await importFile({ name: 'solid.stp', bytes: enc(BREP) }, { adapters: { occt: fake }, units: 'mm' });
  assert.equal(r.format, 'step');
  assert.deepEqual(seen, { name: 'solid.stp', n: BREP.length, units: 'mm', format: 'step' });
  assert.equal(r.meta.via, 'adapter:occt');
  assert.equal(r.meta.tolerance, 0.1);
  assert.deepEqual(r.warnings, ['from adapter']);
  assert.equal(r.meta.counts.triangles, 1);
  const h = await importFile({ name: 'mesh.med', bytes: hdf }, { adapters: { hdf5: fake } });
  assert.equal(h.format, 'med');
  assert.equal(seen.format, 'med');
  await assert.rejects(importFile({ name: 'mesh.med', bytes: hdf }, { adapters: { hdf5: { parse: async () => { throw new Error('corrupt dataset /ENS_MAA'); } } } }), /corrupt dataset/);
  const st = await adapterStatus();
  assert.deepEqual(Object.keys(st).sort(), ['hdf5', 'laz', 'occt', 'sqlite', 'web-ifc']);
  for (const k in st) assert.equal(typeof st[k].available, 'boolean');
});

test('importFile: STEP subsets are read natively; the kernel is used only when the file needs it', async () => {
  const tess = `${P21_HEAD('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF')}#20=COORDINATES_LIST('',3,((0.,0.,0.),(1.,0.,0.),(0.,1.,0.)));\n#21=TRIANGULATED_FACE('f',#20,3,(),$,(),((1,2,3)));\nENDSEC;\nEND-ISO-10303-21;\n`;
  let called = 0;
  const spy = { parse: async () => { called++; return { models: [{ kind: 'surface', name: 'exact', positions: new Float32Array(9), indices: new Uint32Array([0, 1, 2]) }], warnings: [], meta: {} }; } };
  const r = await importFile({ name: 'tess.stp', bytes: enc(tess) }, { units: 'm', adapters: { occt: spy } });
  assert.equal(called, 0);
  assert.equal(r.models[0].kind, 'surface');
  assert.equal(r.meta.via, undefined);
  // mixed file: B-rep plus a curve set
  const mixed = BREP.replace('ENDSEC;\nEND-ISO', "#110=CARTESIAN_POINT('',(0.,0.,0.));\n#111=CARTESIAN_POINT('',(5.,0.,0.));\n#112=POLYLINE('axis',(#110,#111));\n#113=GEOMETRIC_CURVE_SET('w',(#112));\nENDSEC;\nEND-ISO");
  // … without the kernel: the curve is returned with a warning about what is missing
  const m = await importFile({ name: 'mixed.stp', bytes: enc(mixed) }, { units: 'm', adapters: { occt: null } });
  assert.deepEqual(m.models.map((x) => x.kind), ['polyline']);
  assert.ok(m.warnings.some((w) => /"occt" reader is not available in this build, so only the part of the file that can be read without it is shown/.test(w)));
  // … with the kernel: its surfaces plus the curve from the built-in reader
  const k = await importFile({ name: 'mixed.stp', bytes: enc(mixed) }, { units: 'm', adapters: { occt: spy } });
  assert.equal(called, 1);
  assert.deepEqual(k.models.map((x) => [x.kind, x.name]), [['surface', 'exact'], ['polyline', 'axis']]);
  assert.equal(k.meta.via, 'adapter:occt');
  // … with a kernel that fails on this file: fall back to what was read, and say why
  const f = await importFile({ name: 'mixed.stp', bytes: enc(mixed) }, { units: 'm', adapters: { occt: { parse: async () => { throw new Error('kernel ran out of memory'); } } } });
  assert.deepEqual(f.models.map((x) => x.kind), ['polyline']);
  assert.ok(f.warnings.some((w) => /could not process the file \(kernel ran out of memory\)/.test(w)));
});
