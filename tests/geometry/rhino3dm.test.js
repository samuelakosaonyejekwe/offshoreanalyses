import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { deflateSync } from 'node:zlib';
import { parseRhino3dm, readRhino3dm, meshVertices, isRhino3dm } from '../../src/geometry/importers/rhino3dm.js';
import { crc32 } from '../../src/geometry/zip.js';
import { close, concat, assertContract, assertTruncationSafe, surfaceArea, surfaceVolume } from './fixtures/helpers.js';

const real = (name) => new Uint8Array(readFileSync(new URL(`./fixtures/real/${name}`, import.meta.url)));

// ---------------------------------------------------------------------------------------------------------
// A 3DM writer for the tests: chunks, tables and the class data of the objects the reader knows, laid out as
// the openNURBS toolkit writes them. `v` is the archive version: 4-byte chunk lengths below 50, 8-byte from 50.
// ---------------------------------------------------------------------------------------------------------

const TC = {
  COMMENT: 0x00000001, EOF: 0x00007fff, END: 0xffffffff, ANON: 0x40008000, MODEL_ATTRS: 0x40008002,
  PROPERTIES: 0x10000014, SETTINGS: 0x10000015, LAYERS: 0x10000011, IDEFS: 0x10000021, OBJECTS: 0x10000013, LIGHTS: 0x10000012,
  ONVERSION: 0xa0000026, APPLICATION: 0x20008024, UNITS: 0x20008031, LAYER_REC: 0x20008050, IDEF_REC: 0x20008076, OBJECT_REC: 0x20008070,
  OBJ_TYPE: 0x82000071, OBJ_ATTRS: 0x02008072, OBJ_END: 0x8200007f,
  CLASS: 0x00027ffa, CLASS_UUID: 0x0002fffb, CLASS_DATA: 0x0002fffc, CLASS_USERDATA: 0x00027ffd, CLASS_UD_HEADER: 0x0002fff9, CLASS_END: 0x80027fff,
};
const ID = {
  mesh: '4ED7D4E4-E947-11D3-BFE5-0010830122F0', brep: '60B5DBC5-E660-11D3-BFE4-0010830122F0', extrusion: '36F53175-72B8-4D47-BF1F-B4E6FC24F4B9',
  nurbsCurve: '4ED7D4DD-E947-11D3-BFE5-0010830122F0', polyline: '4ED7D4E6-E947-11D3-BFE5-0010830122F0', line: '4ED7D4DB-E947-11D3-BFE5-0010830122F0',
  arc: 'CF33BE2A-09B4-11D4-BFFB-0010830122F0', polycurve: '4ED7D4E0-E947-11D3-BFE5-0010830122F0', nurbsSurface: '4ED7D4DE-E947-11D3-BFE5-0010830122F0',
  plane: '4ED7D4DF-E947-11D3-BFE5-0010830122F0', rev: 'A16220D3-163B-11D4-8000-0010830122F0', sum: 'C4CD5359-446D-4690-9FF5-29059732472B',
  point: 'C3101A1D-F157-11D3-BFE7-0010830122F0', cloud: '2488F347-F8FA-11D3-BFEC-0010830122F0', iref: 'F9CFB638-B9D4-4340-87E3-C56E7865D96A',
  idef: '26F8BFF6-2618-417F-A158-153D64A94989', layer: '95809813-E985-11D3-BFE5-0010830122F0', subd: 'F09BA4D9-455B-42C3-BA3B-E6CCACEF853B',
  doubleVertices: '17F24E75-21BE-4A7B-9F3D-7F85225247E3', tlBrep: 'F06FC243-A32A-4608-9DD8-A7D2C4CE2A36',
};
const u8 = (...v) => Uint8Array.from(v);
const i32 = (...v) => { const b = new Uint8Array(4 * v.length), d = new DataView(b.buffer); v.forEach((x, i) => d.setInt32(4 * i, x | 0, true)); return b; };
const u32 = (...v) => { const b = new Uint8Array(4 * v.length), d = new DataView(b.buffer); v.forEach((x, i) => d.setUint32(4 * i, x >>> 0, true)); return b; };
const f32 = (...v) => { const b = new Uint8Array(4 * v.length), d = new DataView(b.buffer); v.forEach((x, i) => d.setFloat32(4 * i, x, true)); return b; };
const f64 = (...v) => { const b = new Uint8Array(8 * v.length), d = new DataView(b.buffer); v.forEach((x, i) => d.setFloat64(8 * i, x, true)); return b; };
const zeros = (n) => new Uint8Array(n);
const wstr = (s) => { if (!s) return u32(0); const b = new Uint8Array(4 + 2 * (s.length + 1)), d = new DataView(b.buffer); d.setUint32(0, s.length + 1, true); for (let i = 0; i < s.length; i++) d.setUint16(4 + 2 * i, s.charCodeAt(i), true); return b; };
const uuid = (s) => { const h = s.replace(/-/g, ''), b = new Uint8Array(16), d = new DataView(b.buffer); d.setUint32(0, parseInt(h.slice(0, 8), 16), true); d.setUint16(4, parseInt(h.slice(8, 12), 16), true); d.setUint16(6, parseInt(h.slice(12, 16), 16), true); for (let i = 0; i < 8; i++) b[8 + i] = parseInt(h.slice(16 + 2 * i, 18 + 2 * i), 16); return b; };
const objectId = (n) => `00000000-0000-4000-8000-${n.toString(16).toUpperCase().padStart(12, '0')}`;

/** Writer for one archive version. */
function writer(v) {
  const L = v >= 50 ? 8 : 4;
  const val = (x) => { if (L === 4) return u32(x); const b = new Uint8Array(8); new DataView(b.buffer).setBigInt64(0, BigInt(x), true); return b; };
  const w = { v, L };
  /** Long chunk; chunks whose typecode has bit 0x8000 end with a CRC-32 of their body. */
  w.chunk = (tc, ...body) => { const data = concat(...body); const crc = tc & 0x8000 ? u32(crc32(data)) : zeros(0); return concat(u32(tc), val(data.length + crc.length), data, crc); };
  w.short = (tc, x = 0) => concat(u32(tc), val(x));
  w.anon = (...body) => w.chunk(TC.ANON, ...body);
  w.versioned = (major, minor, ...body) => w.chunk(TC.ANON, i32(major, minor), ...body);
  w.table = (tc, ...records) => w.chunk(tc, ...records, w.short(TC.END));
  w.object = (id, data, ...userData) => w.chunk(TC.CLASS, w.chunk(TC.CLASS_UUID, uuid(id)), w.chunk(TC.CLASS_DATA, data), ...userData, w.short(TC.CLASS_END));
  /** Compressed buffer: size, CRC-32, method, then the bytes or a zlib stream in a chunk of its own. */
  w.cbuf = (data, zip) => (data.length ? concat(u32(data.length, crc32(data)), u8(zip ? 1 : 0), zip ? w.anon(new Uint8Array(deflateSync(data))) : data) : u32(0));
  return w;
}

// Attributes: the item list of Rhino 5 and later, or the fixed record of earlier versions (the reader chooses
// by the archive version and the openNURBS version of the properties table).
const attributes = (w, modern, { id = 0, layer = 0, name = '', visible = true, mode = 0 } = {}) => (modern
  ? concat(u8(0x2d), uuid(objectId(id)), i32(layer), name ? concat(u8(1), wstr(name)) : zeros(0), u8(6), zeros(4), u8(11, visible ? 1 : 0), u8(12, mode), u8(0))
  : concat(u8(0x12), uuid(objectId(id)), i32(layer, -1), zeros(4 + 2 + 2 + 8 + 8), i32(1), u8(mode, 0, 0, 0), wstr(name), wstr(''), i32(0), u8(visible ? 1 : 0)));
const record = (w, modern, type, obj, at) => w.chunk(TC.OBJECT_REC, w.short(TC.OBJ_TYPE, type), obj, w.chunk(TC.OBJ_ATTRS, attributes(w, modern, at)), w.short(TC.OBJ_END));
const layerRecord = (w, index, name, { visible = true, id = 100 + index, parent = 0 } = {}) => w.chunk(TC.LAYER_REC, w.object(ID.layer, concat(
  u8(0x16), i32(visible ? 0 : 1, index, 0, -1, 0), zeros(4 + 2 + 2 + 8 + 8), wstr(name), u8(visible ? 1 : 0), i32(-1), zeros(4 + 8), u8(0), uuid(objectId(id)), parent ? uuid(objectId(parent)) : zeros(16), u8(1))));

/** The file: header, comment, properties, settings (unit system), layers, block definitions, objects, end mark. */
function archive(w, { onv, units = 2, custom = null, layers = [], idefs = [], objects = [], application = 'Test writer', end = true }) {
  const head = new TextEncoder().encode(`3D Geometry File Format ${String(w.v).padStart(8, ' ')}`);
  const parts = [
    head, w.chunk(TC.COMMENT, new TextEncoder().encode('written by the test suite\n')),
    w.table(TC.PROPERTIES, w.short(TC.ONVERSION, onv), w.chunk(TC.APPLICATION, u8(0x10), wstr(application), wstr(''), wstr(''))),
    w.table(TC.SETTINGS, w.chunk(TC.UNITS, i32(102, units), f64(0.001, Math.PI / 180, 0.01), i32(0, 3), f64(custom ? custom.scale : 1), wstr(custom ? custom.name : ''))),
    w.table(TC.LAYERS, ...layers), w.table(TC.LIGHTS), w.table(TC.IDEFS, ...idefs), w.table(TC.OBJECTS, ...objects),
  ];
  if (end) { const size = parts.reduce((s, p) => s + p.length, 0) + 4 + 2 * w.L; parts.push(w.chunk(TC.EOF, w.L === 8 ? concat(u32(size), u32(0)) : u32(size))); }
  return concat(...parts);
}

// ---- class data ----
const meshHead = (nv, nf) => concat(i32(nv, nf), zeros(64 + 16), f32(0, 0, 0, 0, 0, 0), zeros(24 + 16), i32(0), u8(0), u8(0, 0, 0, 0));
/** ON_Mesh, chunk version 3.x: faces as quadruples (a triangle repeats its third vertex), compressed vertex buffers. */
function meshData(w, verts, faces, { zip = false, minor = 2, doubles = false, indexSize = 4 } = {}) {
  const nv = verts.length / 3, nf = faces.length / 4;
  const F = indexSize === 4 ? i32(...faces) : indexSize === 2 ? (() => { const b = new Uint8Array(2 * faces.length), d = new DataView(b.buffer); faces.forEach((x, i) => d.setUint16(2 * i, x, true)); return b; })() : u8(...faces);
  const parts = [u8(0x30 | minor), meshHead(nv, nf), i32(indexSize), F, w.cbuf(f32(...verts), zip), u32(0), u32(0), u32(0), u32(0), i32(0)];
  if (minor >= 7) parts.push(zeros(16), u32(0), w.versioned(1, 1, zeros(16 + 4 + 128 + 4)), u8(0, 0, 0), u8(0), u8(doubles ? 1 : 0), doubles ? w.versioned(1, 0, u32(nv), w.cbuf(f64(...verts), zip)) : zeros(0), f64(0, 0, 0, 0, 0, 0));
  return concat(...parts);
}
/** ON_Mesh, chunk version 1.x (Rhino 2): plain arrays. */
const meshDataV1 = (verts, faces) => concat(u8(0x10), meshHead(verts.length / 3, faces.length / 4), i32(4), i32(...faces), i32(verts.length / 3), f32(...verts), i32(0), i32(0), i32(0), i32(0));
const planeData = (o, x, y) => { const z = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]; return f64(...o, ...x, ...y, ...z, ...z, -(z[0] * o[0] + z[1] * o[1] + z[2] * o[2])); };
const lineData = (a, b, t0 = 0, t1 = 1, dim = 3) => concat(u8(0x10), f64(...a, ...b, t0, t1), i32(dim));
const arcData = (o, x, y, r, a0, a1, t0 = a0, t1 = a1) => concat(u8(0x10), planeData(o, x, y), f64(r), zeros(72), f64(a0, a1, t0, t1), i32(3));
const polylineData = (pts) => concat(u8(0x10), i32(pts.length / 3), f64(...pts), i32(pts.length / 3), f64(...Array.from({ length: pts.length / 3 }, (_, i) => i)), i32(3));
/** ON_NurbsCurve: `knots` in openNURBS form (order + count − 2 values), control points homogeneous when rational. */
const nurbsCurveData = (dim, rat, order, knots, cvs) => concat(u8(0x10), i32(dim, rat ? 1 : 0, order, cvs.length / (dim + (rat ? 1 : 0)), 0, 0), zeros(48), i32(knots.length), f64(...knots), i32(cvs.length / (dim + (rat ? 1 : 0))), f64(...cvs));
const nurbsSurfaceData = (order0, order1, n0, n1, k0, k1, cvs, rat = false) => concat(u8(0x10), i32(3, rat ? 1 : 0, order0, order1, n0, n1, 0, 0), zeros(48), i32(k0.length), f64(...k0), i32(k1.length), f64(...k1), i32(n0 * n1), f64(...cvs));
const planeSurfaceData = (o, x, y, du, dv) => concat(u8(0x11), planeData(o, x, y), f64(...du, ...dv, ...du, ...dv));
const polycurveData = (w, segments, params) => concat(u8(0x10), i32(segments.length, 0, 0), zeros(48), i32(params.length), f64(...params), ...segments);
const xform = (m) => f64(...m); // row by row

const QUAD = [0, 0, 0, 2, 0, 0, 2, 3, 0, 0, 3, 0]; // 2 × 3 rectangle in the x–y plane
const QUAD_F = [0, 1, 2, 3];

async function parse(bytes, opts) { return parseRhino3dm({ name: 'model.3dm', bytes }, opts); }
const areaOf = (r) => r.models.filter((m) => m.kind === 'surface').reduce((s, m) => s + surfaceArea(m), 0);

// ---------------------------------------------------------------------------------------------------------
// Meshes, units, layers, attributes
// ---------------------------------------------------------------------------------------------------------

for (const [v, onv] of [[2, 0], [3, 200304080], [4, 200612060], [50, 201004095], [60, 2348834708], [70, 2382428588], [80, 2415953804]]) {
  test(`3DM archive ${v}: mesh with a quad and a triangle, layers, names, hidden objects, millimetres`, async () => {
    const w = writer(v), modern = v >= 5 && onv >= 200712190;
    const verts = [...QUAD, 1, 5, 1], faces = [...QUAD_F, 3, 2, 4, 4];
    const mesh = (zip) => w.object(ID.mesh, v === 2 ? meshDataV1(verts, faces) : meshData(w, verts, faces, { zip }));
    const bytes = archive(w, {
      onv, units: 2,
      layers: [layerRecord(w, 0, 'Hull'), layerRecord(w, 1, 'Plating', { parent: 100 }), layerRecord(w, 2, 'Construction', { visible: false })],
      objects: [
        record(w, modern, 32, mesh(false), { id: 1, layer: 0, name: 'bottom' }),
        record(w, modern, 32, mesh(true), { id: 2, layer: 1, name: 'side' }),
        record(w, modern, 32, mesh(true), { id: 3, layer: 0, name: 'spare', visible: false }),
        record(w, modern, 32, mesh(false), { id: 4, layer: 2, name: 'jig' }),
        record(w, modern, 1, w.object(ID.point, concat(u8(0x10), f64(1000, 2000, 3000))), { id: 5, layer: 0 }),
      ],
    });
    assert.ok(isRhino3dm(bytes));
    const r = await parse(bytes);
    assertContract(r, ['surface', 'points']);
    assert.equal(r.meta.format, 'rhino-3dm');
    assert.equal(r.meta.archiveVersion, v);
    assert.deepEqual(r.meta.units, { unit: 'mm', scale: 1e-3, basis: 'file' });
    assert.deepEqual(r.meta.layers, ['Hull', onv > 200505110 ? 'Hull::Plating' : 'Plating', 'Construction']);
    assert.equal(r.meta.application, 'Test writer');
    const surf = r.models.filter((m) => m.kind === 'surface');
    assert.deepEqual(surf.map((m) => m.name), ['Hull', onv > 200505110 ? 'Hull::Plating' : 'Plating']);
    // quad 2 × 3 mm and a triangle with base 2 and apex (1, 5, 1): 6 + √5 mm² each
    const one = (6 + Math.sqrt(5)) * 1e-6;
    for (const m of surf) { assert.equal(m.indices.length, 9); close(surfaceArea(m), one, 1e-12); }
    assert.deepEqual(r.meta.elementTypes, { quad4: 2, tri3: 2, point: 1 });
    assert.deepEqual(r.meta.skipped, { 'hidden objects': 1, 'objects on hidden layers': 1 });
    assert.deepEqual(r.meta.hidden, { objects: 1, onHiddenLayers: 1, kept: false });
    assert.ok(r.warnings.some((x) => /1 hidden object\(s\) were left out/.test(x)) && r.warnings.some((x) => /on hidden layers/.test(x)));
    const pts = r.models.find((m) => m.kind === 'points');
    assert.deepEqual(Array.from(pts.positions), [1, 2, 3]);
    assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [1, 2, 3] });
    assert.deepEqual(r.meta.objects, { ON_Mesh: 4, ON_Point: 1 });

    // hidden: true keeps everything; units can be overridden
    const all = await parse(bytes, { hidden: true, units: 'm' });
    assert.deepEqual(all.models.filter((m) => m.kind === 'surface').map((m) => m.name), ['Hull', onv > 200505110 ? 'Hull::Plating' : 'Plating', 'Construction']);
    assert.deepEqual(all.models[0].groups.map((g) => g.name), ['bottom', 'spare']);
    close(areaOf(all), 4 * (6 + Math.sqrt(5)), 1e-6);
    assert.deepEqual(all.meta.skipped, {});
    assert.equal(all.meta.units.basis, 'option');
    await assertTruncationSafe(bytes, (b) => parse(b), `3DM archive ${v}`);
  });
}

test('3DM: unit systems of the settings table', async () => {
  const w = writer(60), obj = [record(w, true, 32, w.object(ID.mesh, meshData(w, QUAD, QUAD_F)), { id: 1 })];
  const scaleOf = async (units, custom) => (await parse(archive(w, { onv: 2348834708, units, custom, objects: obj }))).meta.units;
  for (const [code, unit, scale] of [[1, 'um', 1e-6], [2, 'mm', 1e-3], [3, 'cm', 1e-2], [4, 'm', 1], [5, 'km', 1e3], [7, 'mil', 2.54e-5], [8, 'in', 0.0254], [9, 'ft', 0.3048], [10, 'mi', 1609.344], [14, 'dm', 0.1], [19, 'yd', 0.9144], [22, 'nmi', 1852]]) {
    assert.deepEqual(await scaleOf(code), { unit, scale, basis: 'file' }, `unit system ${code}`);
  }
  assert.deepEqual(await scaleOf(6), { unit: 'custom', scale: 2.54e-8, basis: 'file', name: 'microinches' });
  assert.deepEqual(await scaleOf(11, { scale: 0.5, name: 'half metres' }), { unit: 'custom', scale: 0.5, basis: 'file', name: 'half metres' });
  // no unit system: metres are assumed and the warning says so
  const none = await parse(archive(w, { onv: 2348834708, units: 0, objects: obj }));
  assert.equal(none.meta.units.basis, 'assumed');
  assert.ok(none.warnings.some((x) => /declares no length unit; m was assumed/.test(x)));
  close(areaOf(none), 6, 1e-12);
  close(areaOf(await parse(archive(w, { onv: 2348834708, units: 8, objects: obj }))), 6 * 0.0254 ** 2, 1e-10);
});

test('3DM: double-precision vertices (Rhino 6 and later in the mesh, Rhino 5 as user data), index sizes', async () => {
  // a vertex pair that single precision cannot tell apart at a UTM-sized offset
  const E = 500000, N = 6000000, verts = [E, N, 0, E + 0.125, N, 0, E + 0.125, N + 0.0625, 0, E, N + 0.0625, 0];
  const w6 = writer(70);
  const r6 = await parse(archive(w6, { onv: 2382428588, units: 4, objects: [record(w6, true, 32, w6.object(ID.mesh, meshData(w6, verts, QUAD_F, { zip: true, minor: 8, doubles: true })), { id: 1 })] }));
  const m6 = r6.models[0];
  assert.deepEqual(m6.origin, [500000, 6000000, 0]);
  assert.deepEqual(Array.from(m6.positions), [0, 0, 0, 0.125, 0, 0, 0.125, 0.0625, 0, 0, 0.0625, 0]);
  assert.equal(r6.meta.elementTypes['meshes with double-precision vertices'], 1);

  // Rhino 5: ON_V5_MeshDoubleVertices user data behind the class data
  const w5 = writer(50);
  const ud = w5.chunk(TC.CLASS_USERDATA, u8(0x22), w5.chunk(TC.CLASS_UD_HEADER, uuid(ID.doubleVertices), zeros(16 + 4 + 128 + 16 + 1 + 8)), w5.anon(w5.versioned(1, 0, i32(4, 4, 0, 0), i32(4), f64(...verts))));
  const r5 = await parse(archive(w5, { onv: 201004095, units: 4, objects: [record(w5, true, 32, w5.object(ID.mesh, meshData(w5, verts, QUAD_F, { minor: 5 }), ud), { id: 1 })] }));
  assert.deepEqual(Array.from(r5.models[0].positions), [0, 0, 0, 0.125, 0, 0, 0.125, 0.0625, 0, 0, 0.0625, 0]);
  // without them the single-precision vertices are all there is (0.5 m resolution at this northing)
  const rf = await parse(archive(w5, { onv: 201004095, units: 4, objects: [record(w5, true, 32, w5.object(ID.mesh, meshData(w5, verts, QUAD_F, { minor: 5 })), { id: 1 })] }));
  assert.notDeepEqual(Array.from(rf.models[0].positions), Array.from(r5.models[0].positions));

  // 1- and 2-byte face indices; a face that points past the vertex list is dropped and counted
  for (const indexSize of [1, 2]) {
    const r = await parse(archive(w5, { onv: 201004095, objects: [record(w5, true, 32, w5.object(ID.mesh, meshData(w5, QUAD, [...QUAD_F, 0, 1, 9, 9], { indexSize })), { id: 1 })] }));
    assert.equal(r.models[0].indices.length, 6);
    assert.deepEqual(r.meta.skipped, { 'mesh faces with invalid vertex indices': 1 });
  }
  // the decoded archive gives the same vertices
  const ar = readRhino3dm(archive(w6, { onv: 2382428588, objects: [record(w6, true, 32, w6.object(ID.mesh, meshData(w6, QUAD, QUAD_F, { zip: true })), { id: 1 })] }));
  assert.deepEqual(Array.from(meshVertices(ar.objects[0].geom)), QUAD);
});

// ---------------------------------------------------------------------------------------------------------
// B-reps and extrusions
// ---------------------------------------------------------------------------------------------------------

const emptyArray = (w) => w.anon(u8(0x10), i32(0));
/** ON_Brep 3.1 with only its faces and their meshes (what a reader needs when render meshes are stored). */
function meshedBrep(w, faceMeshes, analysis = false) {
  const nf = faceMeshes.length, faces = [];
  for (let i = 0; i < nf; i++) faces.push(i32(i), i32(0), i32(-1, 0, 0));
  const flags = (use) => w.anon(...faceMeshes.map((m) => (use && m ? concat(u8(1), w.object(ID.mesh, m)) : u8(0))));
  return concat(u8(0x31), ...Array.from({ length: 7 }, () => emptyArray(w)), w.anon(u8(0x10), i32(nf), ...faces), f64(0, 0, 0, 1, 1, 1), flags(!analysis), flags(analysis));
}
const BOX_FACES = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]];
const BOX_V = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

test('3DM: the render meshes of a B-rep make a closed solid; analysis meshes are used when there is no render mesh', async () => {
  for (const [v, onv, cls] of [[4, 200612060, ID.brep], [50, 201004095, ID.tlBrep], [80, 2415953804, ID.brep]]) {
    const w = writer(v), modern = v >= 5;
    // one mesh per face of the unit cube, each with its own four vertices (as Rhino stores them)
    const faceMeshes = BOX_FACES.map((q) => meshData(w, q.flatMap((i) => BOX_V[i]), QUAD_F, { zip: true }));
    for (const analysis of [false, true]) {
      const r = await parse(archive(w, { onv, units: 4, objects: [record(w, modern, 16, w.object(cls, meshedBrep(w, faceMeshes, analysis)), { id: 1, name: 'tank' })] }));
      assertContract(r, ['surface']);
      const m = r.models[0];
      assert.equal(m.positions.length / 3, 8, 'the vertices of the face meshes are welded');
      close(surfaceArea(m), 6, 1e-12);
      close(surfaceVolume(m), 1, 1e-12);
      assert.equal(r.meta.elementTypes['B-rep faces from stored meshes'], 6);
      assert.equal(r.meta.tessellation, undefined);
    }
  }
});

/** A planar face trimmed to the rectangle [0, a] × [0, b] of the plane (O, X, Y), with full topology and no mesh. */
function planarBrep(w, a, b, { rev = false } = {}) {
  const P = [[0, 0, 0], [a, 0, 0], [a, b, 0], [0, b, 0]];
  const c2 = [], c3 = [], V = [], E = [], T = [];
  for (let i = 0; i < 4; i++) {
    const p = P[i], q = P[(i + 1) % 4];
    c2.push(concat(i32(1), w.object(ID.line, lineData(p, q, 0, 1, 2))));
    c3.push(concat(i32(1), w.object(ID.line, lineData(p, q, 0, 1, 3))));
    V.push(concat(i32(i), f64(...p), i32(2, (i + 3) % 4, i), f64(0)));
    E.push(concat(i32(i, i, 0), f64(0, 1), i32(i, (i + 1) % 4), i32(1, i), f64(0), f64(0, 1)));
    T.push(concat(i32(i, i), f64(0, 1), i32(i, i, (i + 1) % 4, 0, 1, 0, 0), f64(0, 0), f64(0, 1), zeros(8 + 24), f64(0, 0)));
  }
  return concat(u8(0x31),
    w.anon(u8(0x10), i32(4), ...c2), w.anon(u8(0x10), i32(4), ...c3),
    w.anon(u8(0x10), i32(1), i32(1), w.object(ID.plane, planeSurfaceData([0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, a + 1], [-1, b + 1]))),
    w.anon(u8(0x10), i32(4), ...V), w.anon(u8(0x10), i32(4), ...E), w.anon(u8(0x10), i32(4), ...T),
    w.anon(u8(0x10), i32(1), i32(0), i32(4, 0, 1, 2, 3), i32(1, 0)),
    w.anon(u8(0x10), i32(1), i32(0), i32(1, 0), i32(0, rev ? 1 : 0, 0)),
    f64(0, 0, 0, a, b, 0), w.anon(u8(0)), w.anon(u8(0)), i32(0));
}

test('3DM: a B-rep face without a stored mesh is triangulated between its trimming curves', async () => {
  const w = writer(60);
  const bytes = archive(w, { onv: 2348834708, units: 4, objects: [record(w, true, 16, w.object(ID.brep, planarBrep(w, 4, 2.5)), { id: 1, name: 'plate' })] });
  const r = await parse(bytes);
  assertContract(r, ['surface']);
  // the plane surface is larger than the face: only the trimmed rectangle 4 × 2.5 comes out
  close(areaOf(r), 10, 1e-9);
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [4, 2.5, 0] });
  assert.deepEqual(r.meta.tessellation, { faces: 1, facesRead: 1, facesRelaxed: 0, chord: 0.01 * Math.hypot(4, 2.5), angularStepDeg: 15 });
  assert.ok(r.warnings.some((x) => /1 face\(s\) of 1 object\(s\) carry no render mesh .* triangulated here .* 1 read/.test(x)), r.warnings.join('\n'));
  // the face normal follows the surface (+z); a reversed face is wound the other way
  const nz = (m) => { const p = m.positions, t = m.indices; return (p[3 * t[1]] - p[3 * t[0]]) * (p[3 * t[2] + 1] - p[3 * t[0] + 1]) - (p[3 * t[1] + 1] - p[3 * t[0] + 1]) * (p[3 * t[2]] - p[3 * t[0]]); };
  assert.ok(nz(r.models[0]) > 0);
  const flipped = await parse(archive(w, { onv: 2348834708, units: 4, objects: [record(w, true, 16, w.object(ID.brep, planarBrep(w, 4, 2.5, { rev: true })), { id: 1 })] }));
  assert.ok(nz(flipped.models[0]) < 0);
  // a stored mesh is used when there is one, and can be ignored on request
  const meshed = planarBrep(w, 4, 2.5), cut = meshed.length - 4 - 2 * (4 + w.L + 1 + 4);
  const withMesh = concat(meshed.slice(0, cut), w.anon(u8(1), w.object(ID.mesh, meshData(w, [0, 0, 0, 4, 0, 0, 4, 1, 0, 0, 1, 0], QUAD_F))), w.anon(u8(0)), i32(0));
  const stored = archive(w, { onv: 2348834708, units: 4, objects: [record(w, true, 16, w.object(ID.brep, withMesh), { id: 1 })] });
  close(areaOf(await parse(stored)), 4, 1e-9);
  const again = await parse(stored, { renderMeshes: false });
  close(areaOf(again), 10, 1e-9);
  assert.ok(again.warnings.some((x) => /1 face\(s\) of 1 object\(s\) were triangulated here/.test(x)));
  await assertTruncationSafe(bytes, (b) => parse(b), '3DM B-rep without meshes');
});

test('3DM: surfaces as objects — NURBS patch, plane, surface of revolution, sum surface', async () => {
  const w = writer(70), onv = 2382428588;
  const one = async (id, data, opts) => parse(archive(w, { onv, units: 4, objects: [record(w, true, 8, w.object(id, data), { id: 1 })] }), opts);
  // bilinear patch on the corners of a 3 × 2 rectangle lifted to z = 1: control points with v running fastest
  const patch = await one(ID.nurbsSurface, nurbsSurfaceData(2, 2, 2, 2, [0, 1], [0, 1], [0, 0, 1, 0, 2, 1, 3, 0, 1, 3, 2, 1]));
  close(areaOf(patch), 6, 1e-9);
  assert.deepEqual(patch.meta.bbox, { min: [0, 0, 1], max: [3, 2, 1] });
  const plane = await one(ID.plane, planeSurfaceData([1, 1, 1], [0, 1, 0], [0, 0, 1], [0, 2], [0, 5]));
  close(areaOf(plane), 10, 1e-9);
  assert.deepEqual(plane.meta.bbox, { min: [1, 1, 1], max: [1, 3, 6] });
  // cylinder: the line (2, 0, 0) → (2, 0, 3) revolved about z; area 2πrh from the inscribed side
  const cyl = (transposed) => concat(u8(0x20), f64(0, 0, 0, 0, 0, 1), f64(0, 2 * Math.PI), f64(0, 2 * Math.PI), zeros(48), i32(transposed ? 1 : 0), u8(1), w.object(ID.line, lineData([2, 0, 0], [2, 0, 3])));
  for (const transposed of [false, true]) {
    const coarse = await one(ID.rev, cyl(transposed)), fine = await one(ID.rev, cyl(transposed), { deflection: 1e-4, angularStep: 2 });
    const exact = 2 * Math.PI * 2 * 3;
    assert.ok(areaOf(coarse) < exact && areaOf(coarse) > 0.97 * exact, `coarse cylinder area ${areaOf(coarse)}`);
    assert.ok(areaOf(fine) < exact && areaOf(fine) > 0.9995 * exact && areaOf(fine) > areaOf(coarse), `fine cylinder area ${areaOf(fine)}`);
    for (let i = 0; i < fine.models[0].positions.length; i += 3) close(Math.hypot(fine.models[0].positions[i], fine.models[0].positions[i + 1]), 2, 1e-6);
  }
  // sum surface: a quarter arc of radius 1 in the x–y plane swept along the line (0, 0, 0) → (0, 0, 4)
  const sum = await one(ID.sum, concat(u8(0x10), f64(0, 0, 0), zeros(48), w.object(ID.arc, arcData([0, 0, 0], [1, 0, 0], [0, 1, 0], 1, 0, Math.PI / 2)), w.object(ID.line, lineData([0, 0, 0], [0, 0, 4]))), { deflection: 1e-4, angularStep: 2 });
  assert.ok(Math.abs(areaOf(sum) / (2 * Math.PI) - 1) < 1e-3, `sum surface area ${areaOf(sum)}`);
  assert.deepEqual(sum.meta.elementTypes['surfaces triangulated: sum surface'], 1);
});

test('3DM: extrusion — cached mesh when stored, else walls and caps from the profile', async () => {
  const w = writer(60), onv = 2348834708;
  // closed square profile 2 × 2 in the profile plane, path of length 5 along z, up = y
  const square = w.object(ID.polyline, polylineData([0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0, 0]));
  const body = (cache) => w.versioned(1, 3, square, f64(0, 0, 0, 0, 0, 5), f64(0, 1), f64(0, 1, 0), u8(0, 0), f64(0, 0, 0, 0, 0, 0), f64(0, 5), u8(0), i32(1), u8(1, 1), cache);
  const noCache = w.versioned(1, 0, u8(0));
  const r = await parse(archive(w, { onv, units: 4, objects: [record(w, true, 0x40000000, w.object(ID.extrusion, body(noCache)), { id: 1 })] }));
  const m = r.models[0];
  close(surfaceArea(m), 4 * 2 * 5 + 2 * 4, 1e-9);
  close(Math.abs(surfaceVolume(m)), 20, 1e-9);
  assert.equal(r.meta.elementTypes['extrusions triangulated'], 1);
  // with a cached mesh, that mesh is what comes out
  const cache = w.versioned(1, 0, u8(1), w.versioned(1, 0, uuid('F5E3BAA9-A7A2-49FD-B8A1-66EB274A5F91'), w.object(ID.mesh, meshData(w, QUAD, QUAD_F, { zip: true }))), u8(0));
  const c = await parse(archive(w, { onv, units: 4, objects: [record(w, true, 0x40000000, w.object(ID.extrusion, body(cache)), { id: 1 })] }));
  close(areaOf(c), 6, 1e-12);
  assert.equal(c.meta.elementTypes['extrusions from stored meshes'], 1);
  // mitred ends are not built: counted, with the advice to store meshes
  const mitred = w.versioned(1, 3, square, f64(0, 0, 0, 0, 0, 5), f64(0, 1), f64(0, 1, 0), u8(1, 0), f64(0.6, 0, 0.8, 0, 0, 0), f64(0, 5), u8(0), i32(1), u8(1, 1), noCache);
  const both = await parse(archive(w, { onv, units: 4, objects: [record(w, true, 0x40000000, w.object(ID.extrusion, mitred), { id: 1 }), record(w, true, 0x40000000, w.object(ID.extrusion, body(noCache)), { id: 2 })] }));
  assert.deepEqual(both.meta.skipped, { 'ON_Extrusion without a stored mesh (mitred ends)': 1 });
  assert.ok(both.warnings.some((x) => /Skipped 1 ON_Extrusion without a stored mesh \(mitred ends\); re-save the model in Rhino with render meshes or export STEP/.test(x)));
});

// ---------------------------------------------------------------------------------------------------------
// Curves, points, blocks
// ---------------------------------------------------------------------------------------------------------

test('3DM: curves become polylines — line, arc, polyline, NURBS, polycurve', async () => {
  const w = writer(50), onv = 201004095;
  const s = Math.SQRT1_2;
  const objects = [
    record(w, true, 4, w.object(ID.line, lineData([0, 0, 0], [3, 4, 0])), { id: 1, name: 'keel' }),
    // half circle of radius 2 about (0, 0, 1), parameter 0 … 1
    record(w, true, 4, w.object(ID.arc, arcData([0, 0, 1], [1, 0, 0], [0, 1, 0], 2, 0, Math.PI, 0, 1)), { id: 2 }),
    record(w, true, 4, w.object(ID.polyline, polylineData([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0])), { id: 3, name: 'frame' }),
    // rational quadratic quarter circle of radius 1 (homogeneous control points) and a cubic Bézier
    record(w, true, 4, w.object(ID.nurbsCurve, nurbsCurveData(3, true, 3, [0, 0, 1, 1], [1, 0, 0, 1, s, s, 0, s, 0, 1, 0, 1])), { id: 4 }),
    record(w, true, 4, w.object(ID.nurbsCurve, nurbsCurveData(3, false, 4, [0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 2, 0, 3, 2, 0, 4, 0, 0])), { id: 5 }),
    // polycurve: a line, then a quarter arc that keeps its own parameters
    record(w, true, 4, w.object(ID.polycurve, polycurveData(w, [w.object(ID.line, lineData([0, 0, 0], [1, 0, 0], 5, 7)), w.object(ID.arc, arcData([1, 1, 0], [0, -1, 0], [1, 0, 0], 1, 0, Math.PI / 2))], [0, 1, 2])), { id: 6 }),
  ];
  const r = await parse(archive(w, { onv, units: 4, layers: [layerRecord(w, 0, 'Lines')], objects }));
  assertContract(r, ['polyline']);
  assert.equal(r.models.length, 6);
  const [line, arc, pl, quarter, bez, poly] = r.models;
  assert.equal(line.name, 'keel'); assert.equal(line.layer, 'Lines'); assert.equal(line.entity, 'ON_LineCurve');
  assert.deepEqual(Array.from(line.points), [0, 0, 0, 3, 4, 0]);
  assert.equal(arc.name, 'Lines curve 2');
  for (let i = 0; i < arc.points.length; i += 3) { close(Math.hypot(arc.points[i], arc.points[i + 1]), 2, 1e-12); close(arc.points[i + 2], 1, 1e-12); }
  close(arc.points[0], 2, 1e-12); close(arc.points[arc.points.length - 3], -2, 1e-12);
  assert.equal(pl.closed, true); assert.equal(pl.points.length, 12);
  for (let i = 0; i < quarter.points.length; i += 3) close(Math.hypot(quarter.points[i], quarter.points[i + 1]), 1, 1e-12, 'rational arc');
  // the Bézier at t = 1/2: (0 + 3·1 + 3·3 + 4) / 8, (3·2 + 3·2) / 8
  const mid = (bez.points.length / 3 - 1) / 2;
  close(bez.points[3 * mid], 2, 1e-12); close(bez.points[3 * mid + 1], 1.5, 1e-12);
  assert.deepEqual(Array.from(poly.points.slice(0, 6)), [0, 0, 0, 1, 0, 0]);
  close(poly.points[poly.points.length - 3], 2, 1e-12); close(poly.points[poly.points.length - 2], 1, 1e-12);
  assert.equal(r.meta.elementTypes.curve, 6);
  // curves can be left out
  await assert.rejects(parse(archive(w, { onv, units: 4, objects }), { curves: false }), /none of the 6 object\(s\).*6 × curves \(left out on request\)/);
});

test('3DM: points per layer, point clouds, block instances with transforms and mirroring', async () => {
  for (const [v, onv] of [[50, 201004095], [70, 2382428588]]) {
    const w = writer(v);
    const tetra = meshData(w, [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1], [0, 2, 1, 1, 0, 1, 3, 3, 1, 2, 3, 3, 2, 0, 3, 3]);
    // block definition: Rhino 5 record, or the component form of Rhino 6 and later
    const members = concat(i32(2), uuid(objectId(1)), uuid(objectId(2)));
    const idef = v >= 60
      ? w.versioned(1, 0, w.chunk(TC.MODEL_ATTRS, i32(1, 0), u8(0), u8(1), uuid(objectId(50)), u8(0), u8(0), u8(1), wstr('Bracket')), u32(1), w.versioned(1, 0, i32(4), f64(1), wstr('')), wstr(''), wstr(''), wstr(''), zeros(48), u8(1), members, u8(0))
      : concat(u8(0x16), uuid(objectId(50)), members, wstr('Bracket'), wstr(''), wstr(''), wstr(''), zeros(48), u32(0), wstr(''));
    const ref = (m) => concat(u8(0x10), uuid(objectId(50)), xform(m), zeros(48));
    const objects = [
      record(w, true, 32, w.object(ID.mesh, tetra), { id: 1, mode: 3, layer: 0 }),
      record(w, true, 1, w.object(ID.point, concat(u8(0x10), f64(0.5, 0.5, 0.5))), { id: 2, mode: 3, layer: 0 }),
      record(w, true, 0x1000, w.object(ID.iref, ref([1, 0, 0, 10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])), { id: 3, layer: 1, name: 'port' }),
      // mirrored in x and scaled by 2
      record(w, true, 0x1000, w.object(ID.iref, ref([-2, 0, 0, -10, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1])), { id: 4, layer: 1, name: 'starboard' }),
      record(w, true, 2, w.object(ID.cloud, concat(u8(0x10), i32(3), f64(0, 0, 0, 1, 1, 1, 2, 2, 2), zeros(128 + 48), i32(0))), { id: 5, layer: 1, name: 'scan' }),
      record(w, true, 0x1000, w.object(ID.iref, concat(u8(0x10), uuid(objectId(77)), xform([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]), zeros(48))), { id: 6, layer: 1 }),
    ];
    const r = await parse(archive(w, { onv, units: 4, layers: [layerRecord(w, 0, 'Parts'), layerRecord(w, 1, 'Assembly')], idefs: [w.chunk(TC.IDEF_REC, w.object(ID.idef, idef))], objects }));
    assertContract(r, ['surface', 'points']);
    const s = r.models.find((m) => m.kind === 'surface');
    // the members keep their own layer; each instance is a group named after the instance
    assert.equal(s.layer, 'Parts');
    assert.deepEqual(s.groups.map((g) => g.name), ['port', 'starboard']);
    const part = (g) => ({ positions: s.positions, indices: s.indices.subarray(g.start, g.start + g.count) });
    close(surfaceVolume(part(s.groups[0])), 1 / 6, 1e-9, 'instance volume');
    close(surfaceVolume(part(s.groups[1])), 8 / 6, 1e-9, 'mirrored instance keeps outward orientation');
    const xs = Array.from(s.positions).filter((_, i) => i % 3 === 0);
    assert.equal(Math.min(...xs), -12); assert.equal(Math.max(...xs), 11);
    const pts = r.models.filter((m) => m.kind === 'points');
    assert.deepEqual(pts.map((m) => m.name), ['Parts points', 'scan']);
    assert.deepEqual(Array.from(pts[0].positions), [10.5, 0.5, 0.5, -11, 1, 1]);
    assert.equal(pts[1].positions.length, 9);
    assert.equal(r.meta.blockDefinitions, 1);
    assert.deepEqual(r.meta.skipped, { 'ON_InstanceRef to a block definition that is not in the file': 1 });
    assert.equal(r.meta.elementTypes['block instance'], 2);
  }
});

// ---------------------------------------------------------------------------------------------------------
// What is not read, and damaged input
// ---------------------------------------------------------------------------------------------------------

test('3DM: classes that are not read are counted by name; a file with nothing readable says why', async () => {
  const w = writer(80), onv = 2415953804;
  const subd = record(w, true, 0x40000, w.object(ID.subd, zeros(40)), { id: 1 });
  const text = record(w, true, 0x200, w.object('57376349-62A9-4A16-B411-A46BCD544790', zeros(20)), { id: 2 });
  const unknown = record(w, true, 0, w.object('12345678-1234-1234-1234-123456789ABC', zeros(8)), { id: 3 });
  const broken = record(w, true, 32, w.object(ID.mesh, concat(u8(0x30), i32(4, 1), zeros(20))), { id: 4 });
  const mesh = record(w, true, 32, w.object(ID.mesh, meshData(w, QUAD, QUAD_F)), { id: 5 });
  const r = await parse(archive(w, { onv, objects: [subd, subd, text, unknown, broken, mesh] }));
  assert.equal(r.models.length, 1);
  assert.equal(r.meta.skipped.ON_SubD, 2);
  assert.equal(r.meta.skipped.ON_Text, 1);
  assert.equal(r.meta.skipped['class 12345678-1234-1234-1234-123456789ABC'], 1);
  assert.equal(Object.keys(r.meta.skipped).filter((k) => /^ON_Mesh that could not be decoded \(.*runs past the end of its chunk/.test(k)).length, 1);
  assert.ok(r.warnings.some((x) => /Skipped 2 SubD object\(s\).*ToNURBS or Mesh/.test(x)));
  assert.deepEqual(r.meta.objects, { ON_SubD: 2, ON_Text: 1, 'class 12345678-1234-1234-1234-123456789ABC': 1, ON_Mesh: 2 });
  await assert.rejects(parse(archive(w, { onv, objects: [subd, text] })), /model\.3dm: none of the 2 object\(s\) of the 3DM file could be turned into geometry \(1 × ON_SubD; 1 × ON_Text\)/);
  await assert.rejects(parse(archive(w, { onv, objects: [] })), /model\.3dm: the 3DM file holds no objects \(tables found: properties, settings, layer, light, instance definition, object\)/);
});

test('3DM: wrong signature, version 1, unknown versions, truncated and damaged files', async () => {
  const w = writer(60), onv = 2348834708;
  const good = archive(w, { onv, objects: [record(w, true, 32, w.object(ID.mesh, meshData(w, QUAD, QUAD_F, { zip: true })), { id: 1 })] });
  assert.equal(isRhino3dm(good), true);
  assert.equal(isRhino3dm(good.slice(0, 20)), false);
  await assert.rejects(parse(new Uint8Array(10)), /not a Rhino 3DM file: it holds 10 byte\(s\)/);
  await assert.rejects(parse(new TextEncoder().encode('ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((\'\'),\'2;1\');\nENDSEC;\nDATA;\nENDSEC;\nEND-ISO-10303-21;\n')), /the signature "3D Geometry File Format" is missing/);
  const head = (ver) => { const b = good.slice(); b.set(new TextEncoder().encode(String(ver).padStart(8, ' ')), 24); return b; };
  await assert.rejects(parse(head(1)), /version 1 archive \(Rhino 1\).*save it as a Rhino 2 or later file/);
  await assert.rejects(parse(head(37)), /archive version 37; archives 2, 3, 4, 50, 60, 70 and 80 are known/);
  await assert.rejects(parse(head('abc')), /archive version "abc", which is not a number/);
  // a 64-bit archive read with the version of a 32-bit one: the chunk lengths make no sense
  await assert.rejects(parse(head(4)), /model\.3dm: (the 3DM file is truncated or damaged|the 3DM file holds no objects|none of)/);
  // the file may sit behind a wrapper (OLE storage): the signature is searched for
  const wrapped = await parse(concat(new Uint8Array(512).fill(0xd0), good));
  close(areaOf(wrapped), 6e-6, 1e-12);
  // no end mark
  const open = await parse(archive(w, { onv, end: false, objects: [record(w, true, 32, w.object(ID.mesh, meshData(w, QUAD, QUAD_F)), { id: 1 })] }));
  assert.ok(open.warnings.some((x) => /no end-of-file mark/.test(x)));
  // a chunk that claims more bytes than the file holds
  const cut = good.slice(0, good.length - 60);
  await assert.rejects(parse(cut), /model\.3dm: the 3DM file is truncated or damaged in the section after the instance definition table: the chunk 0x10000013 at offset \d+ claims \d+ byte\(s\), but \d+ are left in the file/);
  // a corrupt zlib stream: the mesh is given up with a reason, not the program
  const zipAt = (() => { for (let i = good.length - 8; i > 0; i--) if (good[i] === 0x78 && good[i + 1] === 0x9c) return i; return -1; })();
  assert.ok(zipAt > 0);
  const bad = good.slice(); for (let i = 2; i < 10; i++) bad[zipAt + i] ^= 0xa5;
  await assert.rejects(parse(bad), /none of the 1 object\(s\).*ON_Mesh that could not be decoded \(the vertex buffer of a mesh could not be inflated/);
  // a buffer that claims an impossible size is refused before anything is allocated
  const huge = good.slice(); new DataView(huge.buffer).setUint32(zipAt - 4 - w.L - 1 - 4 - 4, 0x7fffff00, true);
  await assert.rejects(parse(huge), /none of the 1 object\(s\).*vertex buffer holds 2147483392 bytes where 48 are expected/);
  await assertTruncationSafe(good, (b) => parse(b), '3DM archive 60');
  // garbage behind a valid header, and random damage anywhere: a worded error or a result, never an internal error
  let seed = 7;
  const rnd = () => { seed = (Math.imul(seed, 1103515245) + 12345) >>> 0; return seed / 4294967296; };
  for (let k = 0; k < 200; k++) {
    const b = k < 20 ? concat(good.slice(0, 32), Uint8Array.from({ length: 400 }, () => Math.floor(rnd() * 256))) : good.slice();
    if (k >= 20) for (let n = 0; n < 1 + (k % 5); n++) b[32 + Math.floor(rnd() * (b.length - 32))] = Math.floor(rnd() * 256);
    try { const r = await parse(b); assert.ok(Array.isArray(r.models) && r.meta); } catch (e) {
      assert.ok(e instanceof Error && !(e instanceof TypeError) && !(e instanceof ReferenceError) && /^model\.3dm: /.test(e.message), `damage ${k}: ${e && e.constructor.name}: ${e && e.message}\n${e && e.stack}`);
    }
  }
});

// ---------------------------------------------------------------------------------------------------------
// Files written by Rhino and by the openNURBS bindings (tests/geometry/fixtures/real)
// ---------------------------------------------------------------------------------------------------------

test('3DM written by Rhino 5: mesh, B-rep with its render mesh, 220 polyline curves', async () => {
  const bytes = real('rhino5-mesh-brep-polylines.3dm');
  const ar = readRhino3dm(bytes, 'tutorial.3dm');
  assert.equal(ar.version, 50);
  assert.equal(ar.openNurbs, 201705225);
  assert.equal(ar.application, 'Rhinoceros 5.0');
  assert.equal(ar.unitSystem, 2);
  assert.deepEqual(ar.layers.map((l) => [l.name, l.visible]), [['Default', true]]);
  assert.equal(ar.objects.length, 222);
  // the values openNURBS (rhino3dm 8.35) reports for the same file
  const [mesh, brep] = ar.objects.map((o) => o.geom);
  assert.deepEqual([mesh.t, mesh.nv, mesh.nf], ['mesh', 121, 100]);
  assert.deepEqual([brep.t, brep.nf, brep.meshes.length, brep.meshes[0].nv, brep.meshes[0].nf], ['brep', 1, 1, 25, 16]);
  assert.deepEqual(brep.box, [[15, 0, 0], [25, 10.000000000000002, 0]]);
  assert.equal(meshVertices(brep.meshes[0]).length, 75);

  const r = await parseRhino3dm({ name: 'tutorial.3dm', bytes });
  assertContract(r, ['surface', 'polyline']);
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 1e-3, basis: 'file' });
  assert.deepEqual(r.meta.objects, { ON_Mesh: 1, ON_Brep: 1, ON_PolylineCurve: 220 });
  assert.deepEqual(r.meta.counts.byKind, { surface: 1, polyline: 220 });
  assert.deepEqual(r.meta.skipped, {});
  assert.deepEqual(r.warnings, []);
  const s = r.models[0];
  assert.equal(s.name, 'tutorial');
  assert.deepEqual(s.groups.map((g) => [g.name, g.count / 3]), [['Mesh', 200], ['Brep', 32]]);
  // mesh 10 × 10 mm (area 100 mm²) and the planar face 10 × 10 mm at x = 15 … 25
  close(surfaceArea(s), 200e-6, 1e-12);
  const g = s.groups[1], xs = [];
  for (let i = g.start; i < g.start + g.count; i++) xs.push(s.positions[3 * s.indices[i]]);
  close(Math.min(...xs), 0.015, 1e-9); close(Math.max(...xs), 0.025, 1e-9);
  assert.equal(r.meta.elementTypes.quad4, 116);
  assert.equal(r.meta.elementTypes['B-rep faces from stored meshes'], 1);
  // first and last curve: (40, 9, 0) → (40, 10, 0) and (30, 1, 0) → (31, 1, 0) mm
  const lines = r.models.filter((m) => m.kind === 'polyline');
  for (const [m, want] of [[lines[0], [0.04, 0.009, 0, 0.04, 0.01, 0]], [lines[219], [0.03, 0.001, 0, 0.031, 0.001, 0]]]) Array.from(m.points).forEach((x, i) => close(x, want[i], 1e-12));
  assert.equal(lines[0].name, 'Default curve 1');
  assert.deepEqual(r.meta.bbox.min, [0, 0, 0]);
  await assertTruncationSafe(bytes, (b) => parseRhino3dm({ name: 'tutorial.3dm', bytes: b }), '3DM Rhino 5');
});

test('3DM written through rhino3dm (archive 80): a sphere as a surface of revolution, triangulated here', async () => {
  const bytes = real('rhino3dm-revolved-sphere.3dm');
  const ar = readRhino3dm(bytes);
  assert.equal(ar.version, 80);
  const g = ar.objects[0].geom;
  assert.deepEqual([g.t, g.k, g.curve.k, g.curve.radius, g.transposed], ['surface', 'rev', 'arc', 10, false]);
  assert.deepEqual(g.ang, [0, 2 * Math.PI]);
  const R = 0.01, exactA = 4 * Math.PI * R * R, exactV = (4 / 3) * Math.PI * R ** 3;
  const check = async (opts, tolA, tolV) => {
    const r = await parseRhino3dm({ name: 'sphere.3dm', bytes }, opts);
    assertContract(r, ['surface']);
    assert.deepEqual(r.meta.units, { unit: 'mm', scale: 1e-3, basis: 'file' });
    const m = r.models[0];
    for (let i = 0; i < m.positions.length; i += 3) close(Math.hypot(m.positions[i], m.positions[i + 1], m.positions[i + 2]), R, 1e-8, 'vertex on the sphere');
    const a = surfaceArea(m), v = Math.abs(surfaceVolume(m));
    // inscribed facets: below the exact values, within the tolerance
    assert.ok(a < exactA && a > (1 - tolA) * exactA, `area ${a} of ${exactA}`);
    assert.ok(v < exactV && v > (1 - tolV) * exactV, `volume ${v} of ${exactV}`);
    assert.equal(r.meta.tessellation.facesRead, 1);
    assert.ok(r.warnings.some((x) => /carry no render mesh/.test(x)));
    return { a, v };
  };
  const coarse = await check({}, 0.02, 0.04);
  const fine = await check({ deflection: 5e-6, angularStep: 6 }, 0.004, 0.008);
  assert.ok(fine.a > coarse.a && fine.v > coarse.v);
  await assertTruncationSafe(bytes, (b) => parseRhino3dm({ name: 'sphere.3dm', bytes: b }), '3DM rhino3dm');
});
