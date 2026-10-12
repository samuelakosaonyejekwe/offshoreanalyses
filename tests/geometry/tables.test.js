import test from 'node:test';
import assert from 'node:assert/strict';
import { parseTable, readDelimited, splitHeader, gridFromPoints } from '../../src/geometry/importers/tables.js';
import { parseTableJson, parseXmlModel } from '../../src/geometry/importers/json.js';
import { enc, close, closeArr, surfaceArea, assertContract } from './fixtures/helpers.js';

test('header parsing: units in brackets and suffixes', () => {
  assert.deepEqual(splitHeader('Depth (m)'), { name: 'Depth', unit: 'm' });
  assert.deepEqual(splitHeader('x[mm]'), { name: 'x', unit: 'mm' });
  assert.deepEqual(splitHeader('KP_km'), { name: 'KP', unit: 'km' });
  assert.deepEqual(splitHeader('"Easting"'), { name: 'Easting', unit: null });
});

test('XYZ soundings on a regular lattice become a grid (row 0 = south, missing node = NaN)', () => {
  // 3 × 3 lattice, pitch 10 m, one node missing, rows deliberately shuffled
  const pts = [];
  for (let iy = 0; iy < 3; iy++) for (let ix = 0; ix < 3; ix++) if (!(ix === 1 && iy === 1)) pts.push([100 + 10 * ix, 200 + 10 * iy, -(ix + 10 * iy)]);
  pts.reverse();
  const r = parseTable({ name: 'bathy.xyz', bytes: enc(pts.map((p) => p.join(' ')).join('\n')) });
  assertContract(r, ['grid']);
  const g = r.models[0];
  assert.deepEqual([g.nx, g.ny, g.x0, g.y0, g.dx, g.dy], [3, 3, 100, 200, 10, 10]);
  assert.equal(g.z[0], 0);
  assert.equal(g.z[2], -2);
  assert.equal(g.z[3 * 2 + 1], -21);
  assert.ok(Number.isNaN(g.z[4]));
  assert.equal(r.meta.schema.type, 'grid');
  assert.deepEqual(r.meta.mapping.z, { column: 'c3', index: 2, unit: null });
});

test('scattered soundings stay points; depth column is negated; grid:false keeps points', () => {
  const text = 'Easting,Northing,Depth\n500000.5,6000000.25,35.2\n500013.1,6000004,36.0\n500021.7,6000019.5,34.1\n';
  const r = parseTable({ name: 'soundings.csv', text });
  assertContract(r, ['points']);
  const m = r.models[0];
  assert.deepEqual(m.origin.length, 3, 'UTM-sized coordinates are stored relative to an origin');
  close(m.positions[2] + m.origin[2], -35.2, 1e-5);
  close(m.positions[0] + m.origin[0], 500000.5, 1e-2);
  assert.equal(r.meta.vertical.source, 'depth-positive-down');
  assert.equal(r.meta.mapping.depth.column, 'Depth');
  const lattice = '0 0 1\n1 0 2\n0 1 3\n1 1 4\n';
  assert.equal(parseTable({ name: 'a.xyz', text: lattice }).models[0].kind, 'grid');
  assert.equal(parseTable({ name: 'a.xyz', text: lattice }, { grid: false }).models[0].kind, 'points');
  assert.equal(gridFromPoints('g', new Float64Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 1.3, 1, 0])), null);
});

test('route table: KP in km, easting/northing, depth → polyline with attrs; mapping reported', () => {
  const text = '# crs: EPSG:32631\nKP,Easting (m),Northing (m),Depth (m),OD (mm),Burial\n0,500000,6000000,100,323.9,yes\n0.5,500500,6000000,102,323.9,no\n1.0,501000,6000000,104,323.9,no\n';
  const r = parseTable({ name: 'export_cable.csv', text });
  assertContract(r, ['polyline']);
  const m = r.models[0];
  closeArr([...m.points], [500000, 6000000, -100, 500500, 6000000, -102, 501000, 6000000, -104]);
  closeArr([...m.attrs.kp], [0, 500, 1000]);
  closeArr([...m.attrs.OD], [323.9, 323.9, 323.9]);
  assert.deepEqual(m.attrs.Burial, ['yes', 'no', 'no']);
  assert.equal(m.attrUnits.OD, 'mm');
  assert.equal(r.meta.schema.type, 'route');
  assert.equal(r.meta.mapping.kp.scale, 1000);
  assert.equal(r.meta.crs.source, 'EPSG:32631');
  assert.equal(r.meta.units.basis, 'file');
});

test('route: several lines in one table, explicit mapping and schema, computed chainage', () => {
  const text = 'leg\tE\tN\tlvl\nML1\t0\t0\t-20\nML1\t30\t40\t-80\nML2\t0\t0\t-20\nML2\t-300\t0\t-80\nML2\t-600\t0\t-80\n';
  const r = parseTable({ name: 'moorings.tsv', text }, { schema: 'route', mapping: { x: 'E', y: 'N', z: 'lvl' } });
  assert.deepEqual(r.models.map((m) => [m.kind, m.name, m.points.length / 3]), [['polyline', 'ML1', 2], ['polyline', 'ML2', 3]]);
  close(r.models[0].attrs.kp[1], Math.hypot(30, 40, 60), 1e-9);
  assert.equal(r.meta.mapping.z.column, 'lvl');
  assert.equal(r.meta.delimiter, 'tab');
  assert.throws(() => parseTable({ name: 'm.tsv', text }, { mapping: { x: 'nope' } }), /column "nope" given for role "x" does not exist/);
});

test('profile-only table (chainage + depth) becomes a polyline in the x–z plane', () => {
  const r = parseTable({ name: 'profile.csv', text: 'chainage,depth\n0,10\n100,12\n250,11\n' });
  closeArr([...r.models[0].points], [0, 0, -10, 100, 0, -12, 250, 0, -11]);
  assert.equal(r.meta.schema.type, 'route-profile');
});

test('node + member tables (semicolon, decimal comma) → frame with tube sections in metres', () => {
  const text = '[nodes]\nid;x;y;z\nN1;0,0;0;0\nN2;0;0;10,5\nN3;5;0;10,5\n[members]\nid;start;end;D (mm);t (mm);material;group\nM1;N1;N2;610;19;S355;leg\nM2;N2;N3;406,4;12,7;S355;brace\nM3;N3;N9;100;5;S355;brace\n';
  const r = parseTable({ name: 'jacket.csv', text });
  assertContract(r, ['frame']);
  const f = r.models[0];
  closeArr([...f.nodes], [0, 0, 0, 0, 0, 10.5, 5, 0, 10.5]);
  assert.equal(f.members.length, 2);
  assert.deepEqual(f.members[0], { a: 0, b: 1, section: { type: 'tube', D: 0.61, t: 0.019 }, material: 'S355', group: 'leg', id: 'M1' });
  close(f.members[1].section.D, 0.4064, 1e-12);
  assert.deepEqual(f.nodeNames, ['N1', 'N2', 'N3']);
  assert.equal(r.meta.skipped.membersWithMissingNodes, 1);
  assert.ok(r.warnings.some((w) => /reference node ids that are not in the node table/.test(w)));
});

test('members with inline end coordinates and named sections; nodes are welded', () => {
  const text = 'x1,y1,z1,x2,y2,z2,section\n0,0,0,0,0,5000,PIPE 610x19\n0,0,5000,4000,0,5000,HEB300\n4000,0,5000,4000,0,0,UNKNOWN 5\n';
  const r = parseTable({ name: 'frame.csv', text }, { units: 'mm' });
  const f = r.models[0];
  assert.equal(f.nodes.length / 3, 4);
  assert.deepEqual(f.members[0].section, { type: 'tube', D: 0.61, t: 0.019, name: 'PIPE 610x19' });
  close(f.members[1].section.A, 0.01491, 2e-5, 'HE 300 B area 149.1 cm²');
  assert.equal(f.members[2].sectionName, 'UNKNOWN 5');
  assert.ok(r.warnings.some((w) => /not resolved/.test(w)));
  assert.throws(() => parseTable({ name: 'm.csv', text: 'id,n1,n2\n1,1,2\n' }), /no node table/);
});

test('hull offset matrix → mirrored lofted surface with outward normals', () => {
  // prismatic chine hull: keel (0,0), chine (5,2), sheer (5,4); 30 m long
  const r = parseTable({ name: 'barge.csv', text: 'station,0,2,4\n0,0,5,5\n10,0,5,5\n30,0,5,5\n' });
  assertContract(r, ['surface']);
  const m = r.models[0];
  close(surfaceArea(m), 2 * 30 * (Math.hypot(5, 2) + 2), 1e-4);
  assert.deepEqual(m.groups.map((g) => g.name), ['starboard', 'port']);
  assert.deepEqual(r.meta.bbox, { min: [0, -5, 0], max: [30, 5, 4] });
  // keel vertices are shared between the two halves: 3 stations × (3 + 2) vertices
  assert.equal(m.positions.length / 3, 15);
  const p = m.positions, t = m.indices;
  for (const g of m.groups) {
    let ny = 0;
    for (let i = g.start; i < g.start + g.count; i += 3) {
      const a = 3 * t[i], b = 3 * t[i + 1], c = 3 * t[i + 2];
      ny += (p[b + 2] - p[a + 2]) * (p[c] - p[a]) - (p[b] - p[a]) * (p[c + 2] - p[a + 2]);
    }
    assert.ok(g.name === 'starboard' ? ny > 0 : ny < 0, `${g.name} normals point outboard`);
  }
  assert.equal(m.stations.length, 3);
  const long = parseTable({ name: 'offsets.csv', text: 'station,waterline,half_breadth\n0,0,0\n0,2,5\n0,4,5\n30,0,0\n30,2,5\n30,4,5\n' });
  close(surfaceArea(long.models[0]), 2 * 30 * (Math.hypot(5, 2) + 2), 1e-4);
  assert.equal(long.meta.schema.type, 'hull');
});

test('blade table → span polyline with chord and twist; layout table → named points in local metres', () => {
  const b = parseTable({ name: 'blade.dat', text: 'r (m)  chord (mm)  twist (rad)  airfoil\n61.5 1400 0 NACA64\n2 3500 0.2 cylinder\n30 3000 0.1 DU30\n' });
  const m = b.models[0];
  assert.equal(m.kind, 'polyline');
  closeArr([...m.attrs.r], [2, 30, 61.5]);
  closeArr([...m.attrs.chord], [3.5, 3, 1.4]);
  close(m.attrs.twistDeg[0], (0.2 * 180) / Math.PI, 1e-12);
  assert.deepEqual(m.attrs.airfoil, ['cylinder', 'DU30', 'NACA64']);
  assert.equal(b.meta.schema.type, 'blade');
  const l = parseTable({ name: 'wtg.csv', text: 'name,lon,lat,rating\nA01,3.0,54.0,15\nA02,3.02,54.0,15\n' }, { origin: { lon: 3.0, lat: 54.0 } });
  assertContract(l, ['points']);
  assert.deepEqual(l.models[0].attrs.names, ['A01', 'A02']);
  // 0.02° of longitude at 54°N on WGS 84: N·cos φ·Δλ with N = a / sqrt(1 − e² sin² φ), e² = f(2 − f)
  const f = 1 / 298.257223563, phi = (54 * Math.PI) / 180;
  const east = (6378137 / Math.sqrt(1 - f * (2 - f) * Math.sin(phi) ** 2)) * Math.cos(phi) * ((0.02 * Math.PI) / 180);
  close(east, 1311.5, 0.05, 'hand value');
  close(l.models[0].positions[3], east, 1e-3);
  assert.equal(l.meta.crs.projection, 'local-tangent-plane');
  assert.equal(l.meta.schema.type, 'layout');
});

test('cross-section component table → section loops with layers; panel rows → surface', () => {
  const s = parseTable({ name: 'umbilical.csv', text: 'name,x,y,od,id,material\ntube 1,10,0,12.7,10.2,duplex\ncable,-10,0,8,,copper\n' }, { schema: 'section', units: 'mm' });
  assertContract(s, ['section']);
  assert.equal(s.models[0].loops.length, 3);
  assert.deepEqual(s.models[0].layers.map((l) => [l.name, l.loops.length]), [['tube 1', 2], ['cable', 1]]);
  close(s.models[0].loops[0][0], 0.01 + 0.00635, 1e-12);
  const p = parseTable({ name: 'panels.csv', text: '0 0 0 1 0 0 1 1 0 0 1 0\n1 0 0 2 0 0 2 1 0 2 1 0\n' });
  assertContract(p, ['surface']);
  assert.deepEqual(p.meta.elementTypes, { quad4: 1, tri3: 1 });
  assert.equal(p.models[0].positions.length / 3, 6);
  close(surfaceArea(p.models[0]), 1.5, 1e-12);
});

test('units row, quoted fields, ragged rows and unreadable tables', () => {
  const r = readDelimited({ name: 't.csv', text: '"Name, long","x","y","z"\n,(mm),(mm),(mm)\n"A, 1",1000,0,0\n"B",0,2000\n' });
  const t = r.tables[0];
  assert.equal(t.cols[0].name, 'Name, long');
  assert.equal(t.cols[1].unit, 'mm');
  assert.equal(t.cols[0].str[0], 'A, 1');
  assert.equal(t.ragged, 1);
  assert.ok(Number.isNaN(t.cols[3].num[1]));
  assert.throws(() => parseTable({ name: 'x.csv', text: 'alpha,beta\n1,2\n3,4\n' }), /could not find coordinate columns\. Columns present: alpha, beta/);
  assert.throws(() => parseTable({ name: 'x.csv', text: '# nothing\n\n' }), /no data rows found/);
});

test('JSON: nodes + members with nested positions and units; point rows; hull stations', () => {
  const doc = { units: 'mm', nodes: [{ id: 1, position: [0, 0, 0] }, { id: 2, position: { x: 0, y: 0, z: 5000 } }], members: [{ id: 'a', n1: 1, n2: 2, section: 'PIPE 610x19', material: 'S355' }] };
  const r = parseTableJson({ name: 'm.json', text: JSON.stringify(doc) });
  assertContract(r, ['frame']);
  closeArr([...r.models[0].nodes], [0, 0, 0, 0, 0, 5]);
  assert.deepEqual(r.models[0].members[0].section, { type: 'tube', D: 0.61, t: 0.019, name: 'PIPE 610x19' });
  assert.equal(r.meta.units.unit, 'mm');
  const pts = parseTableJson({ name: 'p.json', text: '{"points":[[0,0,1],[5,0,2],[9,3,4]]}' });
  assert.equal(pts.models[0].kind, 'points');
  const route = parseTableJson({ name: 'riser.json', text: '[{"x":0,"y":0,"z":-100},{"x":50,"y":0,"z":-60},{"x":80,"y":0,"z":0}]' });
  assert.equal(route.models[0].kind, 'polyline', 'file name hints a centreline');
  const hull = parseTableJson({ name: 'hull.json', text: JSON.stringify({ stations: [{ x: 0, offsets: [[0, 0], [5, 2], [5, 4]] }, { x: 30, offsets: [[0, 0], [5, 2], [5, 4]] }] }) });
  close(surfaceArea(hull.models[0]), 2 * 30 * (Math.hypot(5, 2) + 2), 1e-4);
  assert.throws(() => parseTableJson({ name: 'x.json', text: '{"a":1}' }), /no arrays of records were found/);
  assert.throws(() => parseTableJson({ name: 'x.json', text: '{"a":' }), /not valid JSON/);
});

test('XML node/member model → frame', () => {
  const xml = '<?xml version="1.0"?><Model units="m"><Nodes><Node id="1" x="0" y="0" z="0"/><Node id="2" x="0" y="0" z="5"/><Node id="3" x="4" y="0" z="5"/></Nodes><Members><Member id="1" start="1" end="2" section="HEB300"/><Member id="2" start="2" end="3"><section>IPE 300</section></Member></Members></Model>';
  const r = parseXmlModel({ name: 'm.xml', text: xml });
  assertContract(r, ['frame']);
  assert.equal(r.models[0].members.length, 2);
  assert.equal(r.models[0].members[1].section.designation, 'IPE 300');
  assert.equal(r.meta.rootElement, 'Model');
  assert.throws(() => parseXmlModel({ name: 'x.xml', text: '<a><b>text</b></a>' }), /no repeated record elements/);
});
