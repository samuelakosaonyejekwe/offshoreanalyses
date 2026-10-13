// Rhino 3DM — the openNURBS archive written by Rhinoceros 2 and later, by the openNURBS toolkit and its bindings
// (rhino3dm), by MoI and by other programs that use the toolkit.
//
// The file is a tree of typed chunks:
//   header     32 bytes: "3D Geometry File Format " and the archive version as a right-justified number
//              (1, 2, 3, 4, 50, 60, 70, 80; 5 was used by early Rhino 5 work-in-progress builds)
//   chunk      32-bit typecode, then a value of 4 bytes (archives 1 – 4) or 8 bytes (archives 50 and later).
//              Typecodes with the top bit set are "short": the value is the datum. Otherwise it is the length of
//              the body that follows; typecodes with bit 0x8000 end their body with a CRC-32 of it.
//   tables     start section (comment), properties (application, openNURBS version), settings (unit system and
//              tolerances), bitmap, texture mapping, material, linetype, layer, group, font, dimension style,
//              light, hatch pattern, instance definition, object, history record and user tables, end mark.
//              Each table is a chunk of record chunks closed by the short chunk 0xFFFFFFFF.
//   object     record chunk = object type (short), class chunk { class UUID, class data, user data … }, attributes
//              (object id, layer index, name, visibility, mode …), end mark. The class data starts with a
//              version byte (major in the high nibble) and is laid out per class; nested objects (segments of a
//              polycurve, curves and surfaces of a B-rep, the meshes stored with its faces) are class chunks again.
//   buffers    the vertex arrays of meshes are written through ON_BinaryArchive::WriteCompressedBuffer: size,
//              CRC-32 of the data, method (0 stored, 1 zlib stream in a chunk of its own).
// The layouts follow the openNURBS toolkit of Robert McNeel & Associates (opennurbs_archive.cpp, opennurbs_3dm.h,
// opennurbs_3dm_attributes.cpp, opennurbs_layer.cpp, opennurbs_mesh.cpp, opennurbs_brep_io.cpp, opennurbs_beam.cpp,
// opennurbs_instance.cpp, opennurbs_zlib.cpp and the curve and surface classes), which is the published reference
// for the format. The reader below is written for this application from that description.
//
// Read:
//   ON_Mesh            triangles and quads (quads are split along the shorter diagonal); the double-precision
//                      vertices of Rhino 5 (user data) and of Rhino 6 and later when stored
//   ON_Brep            the render mesh stored with each face, else its analysis mesh, welded into one group per
//                      object. A face without a mesh ("Save Small", files of the toolkit and of other programs) is
//                      triangulated in the parameter plane of its surface between its 2-D trimming loops by the
//                      shared B-rep tessellator (brep.js); the 3-D boundary points are those of the edge curves,
//                      which neighbouring faces share. Surfaces: ON_NurbsSurface, ON_PlaneSurface, ON_RevSurface,
//                      ON_SumSurface. The chord error is opts.deflection (default 1 % of the diagonal of the
//                      object) and never more than 2 % of the smaller extent of the face; one facet spans at most
//                      opts.angularStep. The trimmed surfaces of Rhino 2 (B-rep chunk version 2) are read the same way
//   ON_Extrusion       its cached mesh, else walls and caps from the profile and the path (not with mitred ends)
//   ON_NurbsSurface, ON_PlaneSurface, ON_RevSurface, ON_SumSurface as objects: the whole patch
//   ON_LineCurve, ON_ArcCurve, ON_PolylineCurve, ON_NurbsCurve, ON_PolyCurve, the arrows of Rhino 2 → polylines
//   ON_Point, ON_PointCloud → point sets
//   ON_InstanceRef     block instances, nested, with their transforms
//   layers (full path "parent::child"; one surface model per layer, one group per object name), object names,
//   visibility of objects and layers (hidden ones are left out and counted unless opts.hidden), the unit system
// Not read (counted per class in the warnings and in meta.skipped): ON_SubD, annotations, text dots, hatches,
// lights, details, clipping planes, cages, linked (external) block definitions, version 1 archives (refused),
// materials, textures, groups, other user data, history; normals, texture coordinates and colours of meshes.

import { inputBytes, baseName, Grow, surfaceModel, pointsModel, finish, tally, cleanName, M4, V3, Welder, pushPolygon, earclip, groupsFromBuckets, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits } from '../units.js';
import { inflateZlib, inflateZlibSync, crc32 } from '../zip.js';
import { sampleCurve, arcPoints, findSpan, basisFuns } from './nurbs.js';
import { tessellateFace, BodyMesh, CROSSING } from './brep.js';

const SIGNATURE = '3D Geometry File Format ';

const TC = {
  COMMENT: 0x00000001, EOF: 0x00007fff, ENDOFTABLE: 0xffffffff, ANON: 0x40008000, MODEL_ATTRS: 0x40008002,
  PROPERTIES: 0x10000014, SETTINGS: 0x10000015, LAYERS: 0x10000011, IDEFS: 0x10000021, OBJECTS: 0x10000013,
  PROP_APPLICATION: 0x20008024, PROP_ONVERSION: 0xa0000026, SET_UNITS: 0x20008031,
  LAYER_REC: 0x20008050, IDEF_REC: 0x20008076, OBJECT_REC: 0x20008070,
  OBJ_TYPE: 0x82000071, OBJ_ATTRS: 0x02008072, OBJ_END: 0x8200007f,
  CLASS: 0x00027ffa, CLASS_UUID: 0x0002fffb, CLASS_DATA: 0x0002fffc, CLASS_USERDATA: 0x00027ffd, CLASS_UD_HEADER: 0x0002fff9, CLASS_END: 0x80027fff,
};
// short chunks whose value is unsigned (all others are signed)
const UNSIGNED_SHORT = new Set([TC.PROP_ONVERSION, TC.OBJ_TYPE]);
const TABLE_NAMES = { 0x10000010: 'material', 0x10000011: 'layer', 0x10000012: 'light', 0x10000013: 'object', 0x10000014: 'properties', 0x10000015: 'settings', 0x10000016: 'bitmap', 0x10000017: 'user', 0x10000018: 'group', 0x10000019: 'font', 0x10000020: 'dimension style', 0x10000021: 'instance definition', 0x10000022: 'hatch pattern', 0x10000023: 'linetype', 0x10000025: 'texture mapping', 0x10000026: 'history record' };

const NIL = '00000000-0000-0000-0000-000000000000';
const CLS = {
  '4ED7D4E4-E947-11D3-BFE5-0010830122F0': 'ON_Mesh', '60B5DBC5-E660-11D3-BFE4-0010830122F0': 'ON_Brep', '36F53175-72B8-4D47-BF1F-B4E6FC24F4B9': 'ON_Extrusion',
  '4ED7D4DD-E947-11D3-BFE5-0010830122F0': 'ON_NurbsCurve', '4ED7D4E6-E947-11D3-BFE5-0010830122F0': 'ON_PolylineCurve', '4ED7D4DB-E947-11D3-BFE5-0010830122F0': 'ON_LineCurve',
  'CF33BE2A-09B4-11D4-BFFB-0010830122F0': 'ON_ArcCurve', '4ED7D4E0-E947-11D3-BFE5-0010830122F0': 'ON_PolyCurve', '4ED7D4D8-E947-11D3-BFE5-0010830122F0': 'ON_CurveOnSurface',
  '4ED7D4DE-E947-11D3-BFE5-0010830122F0': 'ON_NurbsSurface', '4ED7D4DF-E947-11D3-BFE5-0010830122F0': 'ON_PlaneSurface', 'A16220D3-163B-11D4-8000-0010830122F0': 'ON_RevSurface',
  'C4CD5359-446D-4690-9FF5-29059732472B': 'ON_SumSurface', 'DBC5A584-CE3F-4170-98A8-497069CA5C36': 'ON_ClippingPlaneSurface',
  'C3101A1D-F157-11D3-BFE7-0010830122F0': 'ON_Point', '2488F347-F8FA-11D3-BFEC-0010830122F0': 'ON_PointCloud', '4ED7D4E5-E947-11D3-BFE5-0010830122F0': 'ON_PointGrid',
  'F9CFB638-B9D4-4340-87E3-C56E7865D96A': 'ON_InstanceRef', '26F8BFF6-2618-417F-A158-153D64A94989': 'ON_InstanceDefinition', '95809813-E985-11D3-BFE5-0010830122F0': 'ON_Layer',
  'F09BA4D9-455B-42C3-BA3B-E6CCACEF853B': 'ON_SubD', '0559733B-5332-49D1-A936-0532AC76ADE5': 'ON_Hatch', '85A08513-F383-11D3-BFE7-0010830122F0': 'ON_Light',
  'C8C66EFA-B3CB-4E00-9440-2AD66203379E': 'ON_DetailView', '74198302-CDF4-4F95-9609-6D684F22AB37': 'ON_TextDot', '8BD94E19-59E1-11D4-8018-0010830122F0': 'ON_TextDot (Rhino 2)',
  '06936AFB-3D3C-41AC-BF70-C9319FA480A1': 'ON_NurbsCage', 'D379E6D8-7C31-4407-A913-E3B7040D034A': 'ON_MorphControl', 'D66E5CCF-EA39-11D3-BFE5-0010830122F0': 'ON_Viewport',
  '57376349-62A9-4A16-B411-A46BCD544790': 'ON_Text', '945BF594-6FF9-4F5C-BFC0-B3AF528F29D2': 'ON_Leader', 'E550882B-F44D-4154-A1EF-6E50CBBBF543': 'ON_DimLinear',
  'D417786B-F6CD-4F12-9E1F-063F414DBEB6': 'ON_DimAngular', 'FC749C2F-4C00-41FD-9840-26D94F047AD3': 'ON_DimRadial', '03124828-4C9B-4D28-9A82-664DDDE7A14F': 'ON_DimOrdinate',
  'D46767BA-7E8F-4D9D-9A92-66050219A5B9': 'ON_Centermark', '46F75541-F46B-48BE-AA7E-B353BBE068A7': 'ON_TextObject (Rhino 5)', '14922B7A-5B65-4F11-8345-D415A9637129': 'ON_Leader (Rhino 5)',
  'BD57F33B-A1B2-46E9-9C6E-AF09D30FFDDE': 'ON_DimLinear (Rhino 5)', '841BC40B-A971-4A8E-94E5-BBA26D67348E': 'ON_DimAngular (Rhino 5)', 'B2B683FC-7964-4E96-B1F9-9B356A76B08B': 'ON_DimRadial (Rhino 5)',
  'C8288D69-5BD8-4F50-9BAF-525A0086B0C3': 'ON_DimOrdinate (Rhino 5)', '5DE6B20D-486B-11D4-8014-0010830122F0': 'ON_DimLinear (Rhino 2)', '5DE6B20E-486B-11D4-8014-0010830122F0': 'ON_DimRadial (Rhino 2)',
  '5DE6B20F-486B-11D4-8014-0010830122F0': 'ON_DimAngular (Rhino 2)', '5DE6B210-486B-11D4-8014-0010830122F0': 'ON_TextObject (Rhino 2)', '5DE6B211-486B-11D4-8014-0010830122F0': 'ON_Leader (Rhino 2)',
  '8BD94E1A-59E1-11D4-8018-0010830122F0': 'ON_AnnotationArrow (Rhino 2)',
  // the classes Rhino itself writes for curves, surfaces and B-reps (its internal "TL" classes and their
  // predecessors) share the layout of the openNURBS class
  '5EAF1119-0B51-11D4-BFFE-0010830122F0': 'ON_NurbsCurve', '76A709D5-1550-11D4-8000-0010830122F0': 'ON_NurbsCurve', '4760C817-0BE3-11D4-BFFE-0010830122F0': 'ON_NurbsSurface',
  'FA4FD4B5-1613-11D4-8000-0010830122F0': 'ON_NurbsSurface', 'EF638317-154B-11D4-8000-0010830122F0': 'ON_PolyCurve', '0705FDEF-3E2A-11D4-800E-0010830122F0': 'ON_Brep',
  '2D4CFEDB-3E2A-11D4-800E-0010830122F0': 'ON_Brep', 'F06FC243-A32A-4608-9DD8-A7D2C4CE2A36': 'ON_Brep', '0A8401B6-4D34-4B99-8615-1B4E723DC4E5': 'ON_RevSurface',
  '665F6331-2A66-4CCE-81D0-B5EEBD9B5417': 'ON_SumSurface',
};
const UD_DOUBLE_VERTICES = '17F24E75-21BE-4A7B-9F3D-7F85225247E3', UD_EXTRUSION_MESHES = 'A8130A3E-E4F3-4CB0-BB8A-F10A473912D0';
const MESH_RENDER = 'F5E3BAA9-A7A2-49FD-B8A1-66EB274A5F91', MESH_ANALYSIS = 'AC12F955-A29F-437B-A4C9-06C85B7FF57C';

// ON::LengthUnitSystem → unit name of units.js, or [metres per unit, name] where units.js has none
const UNITS = {
  1: 'um', 2: 'mm', 3: 'cm', 4: 'm', 5: 'km', 6: [2.54e-8, 'microinches'], 7: 'mil', 8: 'in', 9: 'ft', 10: 'mi', 12: [1e-10, 'angstroms'], 13: [1e-9, 'nanometres'], 14: 'dm',
  15: [10, 'dekametres'], 16: [100, 'hectometres'], 17: [1e6, 'megametres'], 18: [1e9, 'gigametres'], 19: 'yd', 20: [0.0254 / 72, 'printer points'], 21: [0.0254 / 6, 'printer picas'], 22: 'nmi',
  23: [1.4959787e11, 'astronomical units'], 24: [9.4607304725808e15, 'light years'], 25: [3.08567758e16, 'parsecs'],
};

/** A flaw in the data of one chunk: the object it belongs to is given up, the file is read on. */
class Flaw extends Error {}
const hex = (tc) => `0x${(tc >>> 0).toString(16).toUpperCase().padStart(8, '0')}`;
const MAX_DEPTH = 24;
// facets across the parameter box of one face (before trimming) at the default facet size: at most MAX_CELLS,
// fewer when a file has so many faces to triangulate that CELL_BUDGET would be passed, never fewer than MIN_CELLS
const MAX_CELLS = 1024, MIN_CELLS = 36, CELL_BUDGET = 100000;

/** True when the bytes start with the 3DM signature. */
export function isRhino3dm(bytes) {
  if (!bytes || bytes.length < 32) return false;
  for (let i = 0; i < 24; i++) if (bytes[i] !== SIGNATURE.charCodeAt(i)) return false;
  return true;
}

// ---------------------------------------------------------------------------------------------------------
// Byte reader with a chunk limit: nothing is read past the end of the chunk that is open
// ---------------------------------------------------------------------------------------------------------

class Reader {
  constructor(bytes) {
    this.b = bytes;
    this.dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    this.pos = 0; this.lim = bytes.length;
    this.L = 4;       // bytes of a chunk value
    this.v = 0;       // archive version
    this.onv = 0;     // openNURBS version that wrote the file
  }
  need(n, what) {
    if (!(n >= 0) || this.pos + n > this.lim) throw new Flaw(`${what || 'a value'} of ${n} byte(s) at offset ${this.pos} runs past the end of ${this.lim === this.b.length ? 'the file' : 'its chunk'} (offset ${this.lim})`);
  }
  skip(n, what) { this.need(n, what); this.pos += n; }
  u8() { this.need(1); return this.b[this.pos++]; }
  bool() { return this.u8() !== 0; }
  i32() { this.need(4); const v = this.dv.getInt32(this.pos, true); this.pos += 4; return v; }
  u32() { this.need(4); const v = this.dv.getUint32(this.pos, true); this.pos += 4; return v; }
  f64() { this.need(8); const v = this.dv.getFloat64(this.pos, true); this.pos += 8; return v; }
  /** n doubles as a Float64Array (copied: the file bytes need not be aligned). */
  f64s(n, what) {
    this.need(8 * n, what);
    const out = new Float64Array(n);
    for (let i = 0, p = this.pos; i < n; i++, p += 8) out[i] = this.dv.getFloat64(p, true);
    this.pos += 8 * n;
    return out;
  }
  pt() { return [this.f64(), this.f64(), this.f64()]; }
  interval() { return [this.f64(), this.f64()]; }
  /** Element count of an array whose elements take `size` bytes each; refused when it cannot fit. */
  count(size, what) {
    const n = this.i32();
    if (n < 0 || n * size > this.lim - this.pos) throw new Flaw(`${what} claims ${n} element(s), more than the ${this.lim - this.pos} byte(s) left in its chunk`);
    return n;
  }
  uuid() {
    this.need(16, 'an id');
    const d = this.dv, p = this.pos, h = (v, n) => v.toString(16).toUpperCase().padStart(n, '0');
    let s = `${h(d.getUint32(p, true), 8)}-${h(d.getUint16(p + 4, true), 4)}-${h(d.getUint16(p + 6, true), 4)}-${h(this.b[p + 8], 2)}${h(this.b[p + 9], 2)}-`;
    for (let i = 10; i < 16; i++) s += h(this.b[p + i], 2);
    this.pos += 16;
    return s;
  }
  /** UTF-16 string: element count including the terminating zero, then the elements. */
  wstr() {
    const n = this.u32();
    if (n > (this.lim - this.pos) / 2) throw new Flaw(`a string claims ${n} characters, more than its chunk holds`);
    let s = '';
    for (let i = 0, p = this.pos; i < n; i++, p += 2) { const c = this.dv.getUint16(p, true); if (c === 0) break; s += String.fromCharCode(c); }
    this.pos += 2 * n;
    return s;
  }
  /** Version byte of a class: major in the high nibble. */
  ver() { const v = this.u8(); return { major: v >> 4, minor: v & 15 }; }
  /** Open the chunk at the current position. A long chunk becomes the limit of all reads until it is closed. */
  open() {
    this.need(4 + this.L, 'a chunk header');
    const at = this.pos, tc = this.dv.getUint32(at, true);
    const short = tc >= 0x80000000;
    let val;
    if (this.L === 8) val = Number(this.dv.getBigInt64(at + 4, true));
    else val = short && !UNSIGNED_SHORT.has(tc) ? this.dv.getInt32(at + 4, true) : this.dv.getUint32(at + 4, true);
    this.pos = at + 4 + this.L;
    const c = { tc, val, at, end: this.pos, short, prev: this.lim };
    if (!short) {
      if (!(val >= 0) || this.pos + val > this.lim) { this.pos = at; throw new Flaw(`the chunk ${hex(tc)} at offset ${at} claims ${val} byte(s), but ${this.lim - at - 4 - this.L} are left ${this.lim === this.b.length ? 'in the file' : 'in the chunk around it'}`); }
      c.end = this.pos + val;
      this.lim = c.end;
    }
    return c;
  }
  /** Leave a chunk (and every chunk opened inside it that was not closed). */
  close(c) { this.pos = c.end; this.lim = c.prev; }
  /** Open a chunk that must have the given typecode. */
  expect(tc, what) {
    const c = this.open();
    if (c.tc !== tc) { this.close(c); throw new Flaw(`${what}: found the chunk ${hex(c.tc)} where ${hex(tc)} should be`); }
    return c;
  }
  /** Open an anonymous chunk that starts with two 32-bit version numbers. */
  versioned(what, tc = TC.ANON) { const c = this.expect(tc, what); c.major = this.i32(); c.minor = this.i32(); return c; }
  skipChunk(what) { this.close(this.expect(TC.ANON, what)); }
  /**
   * Compressed buffer (ON_BinaryArchive::WriteCompressedBuffer). Returns { size, crc, raw | zip } without
   * inflating, or null for an empty buffer. `expect` is the size the buffer must have (null: any).
   */
  cbuf(expect, what, keep = true) {
    const size = this.u32();
    if (!size) return null;
    if (expect !== null && size !== expect) throw new Flaw(`the ${what} buffer holds ${size} bytes where ${expect} are expected`);
    const crc = this.u32(), method = this.u8();
    if (method === 0) { this.need(size, `the ${what} buffer`); const raw = this.b.subarray(this.pos, this.pos + size); this.pos += size; return keep ? { size, crc, raw } : null; }
    if (method !== 1) throw new Flaw(`the ${what} buffer uses the unknown compression method ${method}`);
    const c = this.expect(TC.ANON, `the ${what} buffer`);
    const zip = this.b.subarray(this.pos, Math.max(this.pos, c.end - 4));
    this.close(c);
    // deflate cannot expand more than about 1032 : 1
    if (size > 1040 * zip.length + 64) throw new Flaw(`the ${what} buffer claims ${size} bytes from ${zip.length} compressed bytes, which deflate cannot produce`);
    return keep ? { size, crc, zip } : null;
  }
}

// ---------------------------------------------------------------------------------------------------------
// Class data
// ---------------------------------------------------------------------------------------------------------

function readPlane(r) {
  const o = r.pt(), x = r.pt(), y = r.pt(), z = r.pt();
  r.skip(32, 'a plane equation');
  return { o, x, y, z };
}

function readNurbsCurve(r) {
  const v = r.ver();
  if (v.major !== 1) throw new Flaw(`ON_NurbsCurve chunk version ${v.major}.${v.minor} is not known`);
  const dim = r.i32(), rat = r.i32() ? 1 : 0, order = r.i32(), n = r.i32();
  r.skip(8 + 48, 'a NURBS curve header');
  if (!(dim >= 1 && dim <= 16) || !(order >= 2 && order <= 64) || !(n >= order)) throw new Flaw(`ON_NurbsCurve with dimension ${dim}, order ${order} and ${n} control points`);
  const nk = r.count(8, 'the knot vector of a NURBS curve');
  if (nk !== order + n - 2) throw new Flaw(`ON_NurbsCurve of order ${order} with ${n} control points has ${nk} knots instead of ${order + n - 2}`);
  const knots = r.f64s(nk);
  const size = dim + rat, nc = r.count(8 * size, 'the control points of a NURBS curve');
  if (nc !== n) throw new Flaw(`ON_NurbsCurve lists ${nc} control points instead of ${n}`);
  const cv = r.f64s(n * size);
  return { t: 'curve', k: 'nurbs', dim, rat, order, n, knots, cv, dom: [knots[order - 2], knots[n - 1]] };
}

function readNurbsSurface(r) {
  const v = r.ver();
  if (v.major !== 1) throw new Flaw(`ON_NurbsSurface chunk version ${v.major}.${v.minor} is not known`);
  const dim = r.i32(), rat = r.i32() ? 1 : 0, o0 = r.i32(), o1 = r.i32(), n0 = r.i32(), n1 = r.i32();
  r.skip(8 + 48, 'a NURBS surface header');
  if (!(dim >= 1 && dim <= 16) || !(o0 >= 2 && o0 <= 64 && o1 >= 2 && o1 <= 64) || !(n0 >= o0 && n1 >= o1)) throw new Flaw(`ON_NurbsSurface with dimension ${dim}, orders ${o0} × ${o1} and ${n0} × ${n1} control points`);
  const k0 = r.count(8, 'a knot vector of a NURBS surface');
  if (k0 !== o0 + n0 - 2) throw new Flaw(`ON_NurbsSurface has ${k0} knots in u instead of ${o0 + n0 - 2}`);
  const U = r.f64s(k0);
  const k1 = r.count(8, 'a knot vector of a NURBS surface');
  if (k1 !== o1 + n1 - 2) throw new Flaw(`ON_NurbsSurface has ${k1} knots in v instead of ${o1 + n1 - 2}`);
  const V = r.f64s(k1);
  const size = dim + rat, nc = r.count(8 * size, 'the control points of a NURBS surface');
  if (nc !== n0 * n1) throw new Flaw(`ON_NurbsSurface lists ${nc} control points instead of ${n0 * n1}`);
  return { t: 'surface', k: 'nurbs', dim, rat, order: [o0, o1], n: [n0, n1], knots: [U, V], cv: r.f64s(nc * size) };
}

/** Mesh: faces and the (still compressed) vertex buffers. */
function readMesh(r) {
  const v = r.ver();
  if (v.major !== 1 && v.major !== 3) throw new Flaw(`ON_Mesh chunk version ${v.major}.${v.minor} is not known`);
  const nv = r.i32(), nf = r.i32();
  if (nv < 0 || nf < 0) throw new Flaw(`ON_Mesh with ${nv} vertices and ${nf} faces`);
  // texture and surface domains, scales, vertex / normal / texture boxes, closed flag
  r.skip(64 + 16 + 24 + 24 + 16 + 4, 'a mesh header');
  if (r.u8()) r.skipChunk('the meshing parameters of a mesh');
  for (let i = 0; i < 4; i++) if (r.u8()) r.skipChunk('the curvature statistics of a mesh');
  const isz = r.i32();
  if (isz !== 1 && isz !== 2 && isz !== 4) { if (nf) throw new Flaw(`ON_Mesh face indices of ${isz} byte(s)`); }
  const F = new Uint32Array(4 * nf);
  if (nf) {
    r.need(4 * nf * isz, 'the faces of a mesh');
    const d = r.dv;
    let p = r.pos;
    if (isz === 1) for (let i = 0; i < 4 * nf; i++) F[i] = r.b[p++];
    else if (isz === 2) for (let i = 0; i < 4 * nf; i++, p += 2) F[i] = d.getUint16(p, true);
    else for (let i = 0; i < 4 * nf; i++, p += 4) F[i] = d.getUint32(p, true);
    r.pos = p;
  }
  const m = { t: 'mesh', nv, nf, F, vf: null, vd: null, V: null, precision: 'single' };
  if (v.major === 1) {
    const n = r.count(12, 'the vertices of a mesh');
    m.vf = { size: 12 * n, crc: null, raw: r.b.subarray(r.pos, r.pos + 12 * n) };
    r.pos += 12 * n;
    m.nv = n;
    // normals, texture coordinates, curvatures, colours
    for (const s of [12, 8, 16, 4]) { const k = r.count(s, 'a vertex array of a mesh'); r.pos += s * k; }
  } else {
    if (nv) {
      m.vf = r.cbuf(12 * nv, 'vertex');
      r.cbuf(12 * nv, 'vertex normal', false); r.cbuf(8 * nv, 'texture coordinate', false); r.cbuf(16 * nv, 'vertex curvature', false); r.cbuf(4 * nv, 'vertex colour', false);
    }
  }
  if (nv && !m.vf) throw new Flaw(`ON_Mesh declares ${nv} vertices but stores none`);
  // what follows the vertex arrays only matters for the double-precision vertices of Rhino 6 and later
  try {
    if (v.minor >= 2) r.i32();
    if (v.major === 3 && v.minor >= 3) {
      r.skip(16);
      if (nv) r.cbuf(null, 'surface parameter', false);
      if (v.minor >= 4 && r.onv >= 200606010) {
        r.skipChunk('the texture mapping tag of a mesh');
        if (v.minor >= 5) {
          r.skip(3);
          if (v.minor >= 6) {
            if (r.bool()) r.skipChunk('the n-gons of a mesh');
            if (v.minor >= 7 && r.bool()) {
              const c = r.versioned('the double-precision vertices of a mesh');
              if (c.major === 1) { const n = r.u32(); if (n && n === nv) m.vd = r.cbuf(24 * n, 'double-precision vertex'); }
              r.close(c);
            }
          }
        }
      }
    }
  } catch (e) { if (!(e instanceof Flaw)) throw e; m.vd = null; }
  return m;
}

// ---------------------------------------------------------------------------------------------------------
// Archive
// ---------------------------------------------------------------------------------------------------------

/**
 * Decode a 3DM archive into its settings, layers, block definitions and objects. Mesh vertex buffers stay
 * compressed (see meshVertices); the curves and surfaces of a B-rep are decoded when a face needs them.
 * @param {Uint8Array} bytes @param {string} [fname]
 * @returns {{ version: number, openNurbs: number, comment: string, application: string, unitSystem: number|null,
 *   metresPerUnit: number|null, unitName: string, layers: object[], idefs: Map<string, object>, objects: object[],
 *   tables: string[], reader: object }}
 */
export function readRhino3dm(bytes, fname = '3DM') {
  if (!bytes || bytes.length < 32) throw new Error(`${fname}: not a Rhino 3DM file: it holds ${bytes ? bytes.length : 0} byte(s), fewer than the 32-byte header "3D Geometry File Format".`);
  let start = 0;
  if (!isRhino3dm(bytes)) {
    // an OLE wrapper may precede the archive (files saved through some shell extensions)
    start = -1;
    const lim = Math.min(bytes.length - 32, 1 << 20);
    for (let i = 1; i <= lim && start < 0; i++) { if (bytes[i] !== 0x33 || bytes[i + 1] !== 0x44 || bytes[i + 2] !== 0x20) continue; let ok = true; for (let k = 3; k < 24 && ok; k++) ok = bytes[i + k] === SIGNATURE.charCodeAt(k); if (ok) start = i; }
    if (start < 0) throw new Error(`${fname}: not a Rhino 3DM file: the signature "3D Geometry File Format" is missing.`);
  }
  let vs = '';
  for (let i = 24; i < 32; i++) vs += String.fromCharCode(bytes[start + i]);
  vs = vs.trim().replace(/X$/, '2');
  if (!/^\d+$/.test(vs) || !(Number(vs) > 0)) throw new Error(`${fname}: the 3DM header gives the archive version "${vs}", which is not a number.`);
  const version = Number(vs);
  if (version === 1) throw new Error(`${fname}: this is a version 1 archive (Rhino 1), which stores its objects in an older record format that is not read. Open it in Rhino and save it as a Rhino 2 or later file.`);
  if (!(version >= 2 && version <= 5) && !(version >= 50 && version % 10 === 0)) throw new Error(`${fname}: the 3DM header gives the archive version ${version}; archives 2, 3, 4, 50, 60, 70 and 80 are known.`);
  const r = new Reader(start ? bytes.subarray(start) : bytes);
  r.v = version; r.L = version >= 50 ? 8 : 4; r.pos = 32;
  const ar = { version, openNurbs: 0, comment: '', application: '', unitSystem: null, metresPerUnit: null, unitName: '', layers: [], idefs: new Map(), objects: [], tables: [], flaws: [], reader: r, newer: version > 80 };
  const damaged = (e, where) => new Error(`${fname}: the 3DM file is truncated or damaged in ${where}: ${e.message}.`);

  /** Records of a table: `fn` is called inside each record chunk of type `recTc`. */
  const table = (c, recTc, what, fn) => {
    while (r.pos + 4 + r.L <= c.end) {
      let k;
      try { k = r.open(); } catch (e) { if (e instanceof Flaw) throw damaged(e, `the ${what} table`); throw e; }
      if (k.tc === TC.ENDOFTABLE) { r.close(k); break; }
      try { if (k.tc === recTc || recTc === null) fn(k); } catch (e) { if (!(e instanceof Flaw)) throw e; ar.flaws.push(`${what} table: ${e.message}`); } finally { r.close(k); }
    }
  };

  let sawEnd = false;
  while (r.pos + 4 + r.L <= r.lim) {
    let c;
    try { c = r.open(); } catch (e) { if (e instanceof Flaw) throw damaged(e, ar.tables.length ? `the section after the ${ar.tables[ar.tables.length - 1]} table` : 'its first section'); throw e; }
    if (TABLE_NAMES[c.tc]) ar.tables.push(TABLE_NAMES[c.tc]);
    switch (c.tc) {
      case TC.COMMENT: { let s = ''; for (let p = r.pos; p < c.end && s.length < 2000; p++) { const ch = r.b[p]; if (ch === 0 || ch === 26) break; s += String.fromCharCode(ch); } ar.comment = s.trim(); break; }
      case TC.PROPERTIES:
        table(c, null, 'properties', (k) => {
          if (k.tc === TC.PROP_ONVERSION) { ar.openNurbs = k.val > 0 ? k.val : 0; r.onv = ar.openNurbs; } else if (k.tc === TC.PROP_APPLICATION) { r.ver(); ar.application = r.wstr(); }
        });
        break;
      case TC.SETTINGS:
        table(c, TC.SET_UNITS, 'settings', () => {
          const ver = r.i32();
          if (!(ver >= 100 && ver < 200)) return;
          ar.unitSystem = r.i32();
          r.skip(24);
          if (ver >= 101) { r.skip(8); if (ver >= 102) { ar.metresPerUnit = r.f64(); ar.unitName = r.wstr(); } }
        });
        break;
      case TC.LAYERS:
        table(c, TC.LAYER_REC, 'layer', () => { const l = readObject(r, ar, 0); if (l && l.t === 'layer') ar.layers.push(l); });
        break;
      case TC.IDEFS:
        table(c, TC.IDEF_REC, 'instance definition', () => { const d = readObject(r, ar, 0); if (d && d.t === 'idef') ar.idefs.set(d.id, d); });
        break;
      case TC.OBJECTS:
        table(c, TC.OBJECT_REC, 'object', (k) => {
          const t = r.expect(TC.OBJ_TYPE, 'an object record');
          const type = t.val;
          r.close(t);
          const geom = readObject(r, ar, 0);
          const attrs = { id: NIL, layer: 0, name: '', visible: true, mode: 0 };
          while (r.pos + 4 + r.L <= k.end) {
            const a = r.open();
            try { if (a.tc === TC.OBJ_ATTRS) readAttributes(r, attrs); } catch (e) { if (!(e instanceof Flaw)) throw e; } finally { r.close(a); }
            if (a.tc === TC.OBJ_END) break;
          }
          ar.objects.push({ type, geom: geom || { t: 'none', cls: 'empty object record' }, attrs });
        });
        break;
      case TC.EOF: sawEnd = true; break;
      default: break;
    }
    r.close(c);
    if (sawEnd) break;
  }
  ar.complete = sawEnd;
  return ar;
}

/** Object attributes: id, layer, name, visibility and mode (0 normal, 1 hidden, 2 locked, 3 member of a block definition). */
function readAttributes(r, a) {
  if (r.v >= 5 && r.onv >= 200712190) {
    const v = r.ver();
    if (v.major !== 2) throw new Flaw(`object attributes of chunk version ${v.major}.${v.minor}`);
    a.id = r.uuid(); a.layer = r.i32();
    // items in ascending order of their id, closed by 0; nothing after the mode is needed
    let id = r.u8();
    const next = () => { id = r.u8(); return id; };
    if (id === 1) { a.name = r.wstr(); next(); }
    if (id === 2) { r.wstr(); next(); }
    if (id === 3) { r.skip(4); next(); }
    if (id === 4) { r.skip(4); next(); }
    if (id === 5) { r.skipChunk('rendering attributes'); next(); }
    if (id === 6) { r.skip(4); next(); }
    if (id === 7) { r.skip(4); next(); }
    if (id === 8) { r.skip(8); next(); }
    if (id === 9) { r.skip(1); next(); }
    if (id === 10) { r.skip(4); next(); }
    if (id === 11) { a.visible = r.bool(); next(); }
    if (id === 12) { a.mode = r.u8(); }
    return;
  }
  const v = r.ver();
  if (v.major !== 1) throw new Flaw(`object attributes of chunk version ${v.major}.${v.minor}`);
  a.id = r.uuid(); a.layer = r.i32();
  r.skip(4 + 4 + 2 + 2 + 8 + 8 + 4, 'object attributes');
  a.mode = r.u8();
  r.skip(3);
  a.name = r.wstr(); r.wstr();
  a.visible = a.mode !== 1;
  if (v.minor >= 1) { const n = r.count(4, 'the group list of an object'); r.pos += 4 * n; if (v.minor >= 2) a.visible = r.bool(); }
}

function readLayer(r) {
  const v = r.ver();
  if (v.major !== 1) throw new Flaw(`ON_Layer chunk version ${v.major}.${v.minor} is not known`);
  const mode = r.i32(), index = r.i32();
  r.skip(4 + 4 + 4 + 4 + 2 + 2 + 8 + 8, 'a layer record');
  const l = { t: 'layer', index, name: r.wstr(), visible: mode !== 1, id: NIL, parent: NIL };
  if (v.minor >= 1) {
    l.visible = r.bool();
    if (v.minor >= 2) { r.skip(4); if (v.minor >= 3) { r.skip(12); if (v.minor >= 4) { r.skip(1); if (v.minor >= 5) { l.id = r.uuid(); if (v.minor >= 6 && r.onv > 200505110) l.parent = r.uuid(); } } } }
  }
  return l;
}

function readInstanceDefinition(r) {
  const uuids = () => { const n = r.count(16, 'the object list of a block definition'); const out = []; for (let i = 0; i < n; i++) out.push(r.uuid()); return out; };
  const v6 = r.v > 60 || (r.v === 60 && r.onv > 2348834153 && r.pos + 4 <= r.lim && r.dv.getUint32(r.pos, true) === TC.ANON);
  if (!v6) {
    const v = r.ver();
    if (v.major !== 1) throw new Flaw(`ON_InstanceDefinition chunk version ${v.major}.${v.minor} is not known`);
    const d = { t: 'idef', id: r.uuid(), members: uuids(), name: r.wstr(), linked: false };
    try { r.wstr(); r.wstr(); r.wstr(); r.skip(48); const type = r.u32(); d.linked = (type === 2 || type === 3) && !d.members.length; } catch (e) { if (!(e instanceof Flaw)) throw e; }
    return d;
  }
  const c = r.versioned('a block definition');
  if (c.major !== 1) throw new Flaw(`ON_InstanceDefinition chunk version ${c.major}.${c.minor} is not known`);
  const d = { t: 'idef', id: NIL, members: [], name: '', linked: false };
  const m = r.versioned('the component attributes of a block definition', TC.MODEL_ATTRS);
  if (r.u8() === 1) r.skip(12);
  if (r.u8() === 1) d.id = r.uuid();
  if (r.u8() === 1) r.skip(4);
  if (r.u8() === 1) r.skip(4);
  if (r.u8() === 1) d.name = r.wstr();
  r.close(m);
  const type = r.u32();
  r.skipChunk('the unit system of a block definition');
  r.wstr(); r.wstr(); r.wstr(); r.skip(48);
  if (r.bool()) d.members = uuids();
  d.linked = type === 3 && !d.members.length;
  r.close(c);
  return d;
}

function readPolyCurve(r, ar, depth) {
  r.ver();
  const n = r.i32();
  r.skip(8 + 48, 'a polycurve header');
  const nt = r.count(8, 'the parameters of a polycurve'), T = r.f64s(nt);
  if (n < 0 || n > (r.lim - r.pos) / 8 || nt !== n + 1) throw new Flaw(`ON_PolyCurve with ${n} segments and ${nt} parameters`);
  const segs = [];
  for (let i = 0; i < n; i++) { const s = readObject(r, ar, depth + 1); if (!s || s.t !== 'curve') throw new Flaw(`segment ${i + 1} of an ON_PolyCurve is ${s ? s.cls || s.t : 'missing'}`); segs.push(s); }
  return { t: 'curve', k: 'poly', segs, T, dom: [T[0], T[n]] };
}

/**
 * Trimmed surface of Rhino 2 files (B-rep chunk version 2): counts, the 2-D and 3-D curves as polycurves, the
 * surfaces as NURBS surfaces, faces with their loops and trims, then the face meshes.
 */
function readLegacyBrep(r, ar, depth, minor) {
  const nf = r.i32(), ne = r.i32(), nl = r.i32(), nt = r.i32();
  if (!(nf >= 1 && ne >= 1 && nl >= 1 && nt >= 1) || nf > r.lim - r.pos || ne > r.lim - r.pos || nt > r.lim - r.pos) throw new Flaw(`legacy ON_Brep with ${nf} faces, ${ne} edges, ${nl} loops and ${nt} trims`);
  r.skip(4);
  const b = { t: 'brep', nf, meshes: new Array(nf).fill(null), box: [r.pt(), r.pt()], solid: 0, topo: null };
  const t = { C2: [], C3: [], S: [], V: [], E: [], T: [], L: [], F: [] };
  for (let i = 0; i < nt; i++) t.C2.push(readPolyCurve(r, ar, depth + 1));
  for (let i = 0; i < ne; i++) { const c = readPolyCurve(r, ar, depth + 1); t.C3.push(c); t.E.push({ c3: i, rev: false, dom: c.dom, v: [-1, -1], pts: null }); }
  for (let i = 0; i < nf; i++) t.S.push(readNurbsSurface(r));
  for (let fi = 0; fi < nf; fi++) {
    r.skip(8);
    const f = { loops: [], s: fi, rev: r.i32() !== 0 };
    r.skip(4 + 48);
    const nb = r.count(44, 'the loops of a legacy B-rep face');
    for (let k = 0; k < nb; k++) {
      r.skip(4);
      const type = r.i32(), l = { trims: [], type: type === 0 ? 1 : type === 1 ? 2 : 3 };
      r.skip(32);
      const n = r.count(41, 'the trims of a legacy B-rep loop');
      for (let j = 0; j < n; j++) { r.skip(4 + 4 + 1); const e = r.i32(), rev = r.i32() !== 0; r.skip(4 + 4 + 16); const c2 = t.T.length; l.trims.push(c2); t.T.push({ c2, dom: c2 < nt ? t.C2[c2].dom : [0, 0], e, v: [-1, -1], rev, flip: false }); }
      f.loops.push(t.L.length); t.L.push(l);
    }
    t.F.push(f);
  }
  b.topo = t;
  try {
    for (let pass = 0; pass < (minor >= 1 ? 2 : 1); pass++) for (let i = 0; i < nf; i++) if (r.u8()) { const m = readObject(r, ar, depth + 1); if (m && m.t === 'mesh' && !b.meshes[i]) b.meshes[i] = m; }
  } catch (e) { if (!(e instanceof Flaw)) throw e; }
  return b;
}

/** Array of polymorphic curves or surfaces of a B-rep. */
function readObjectArray(r, ar, depth, what) {
  const c = r.expect(TC.ANON, what);
  const v = r.ver();
  if (v.major !== 1) throw new Flaw(`${what} of chunk version ${v.major}.${v.minor}`);
  const n = r.count(4, what), out = new Array(n).fill(null);
  for (let i = 0; i < n; i++) if (r.i32() === 1) out[i] = readObject(r, ar, depth + 1);
  r.close(c);
  return out;
}

/** B-rep: the chunks of its geometry and topology are located, its face meshes are read. */
function readBrep(r, ar, depth) {
  const v = r.ver();
  if (v.major === 2) return readLegacyBrep(r, ar, depth, v.minor);
  if (v.major !== 3) throw new Flaw(`ON_Brep chunk version ${v.major}.${v.minor} is not known`);
  const span = (what) => { const c = r.expect(TC.ANON, what); const s = { at: c.at, end: c.end }; r.close(c); return s; };
  const b = { t: 'brep', c2: span('the 2-D curves of a B-rep'), c3: span('the 3-D curves of a B-rep'), s: span('the surfaces of a B-rep'), v: span('the vertices of a B-rep'), e: span('the edges of a B-rep'), tr: span('the trims of a B-rep'), l: span('the loops of a B-rep'), nf: 0, meshes: [], topo: null, solid: 0 };
  const f = r.expect(TC.ANON, 'the faces of a B-rep');
  b.f = { at: f.at, end: f.end };
  r.ver();
  b.nf = r.count(16, 'the faces of a B-rep');
  r.close(f);
  b.box = [r.pt(), r.pt()];
  b.meshes = new Array(b.nf).fill(null);
  if (v.minor >= 1) {
    // render meshes, then analysis meshes: one flag byte per face, a mesh object where it is set
    for (let pass = 0; pass < 2; pass++) {
      const c = r.expect(TC.ANON, 'the face meshes of a B-rep');
      try {
        for (let i = 0; i < b.nf; i++) if (r.u8()) { const m = readObject(r, ar, depth + 1); if (m && m.t === 'mesh' && !b.meshes[i]) b.meshes[i] = m; }
      } catch (e) { if (!(e instanceof Flaw)) throw e; } finally { r.close(c); }
    }
    if (v.minor >= 2 && r.pos + 4 <= r.lim) b.solid = r.i32();
  }
  return b;
}

function readExtrusion(r, ar, depth) {
  const c = r.versioned('an extrusion');
  if (c.major !== 1) throw new Flaw(`ON_Extrusion chunk version ${c.major}.${c.minor} is not known`);
  const x = { t: 'extrusion', profile: readObject(r, ar, depth + 1), path: [r.pt(), r.pt()], pt: r.interval(), up: r.pt(), haveN: [r.bool(), r.bool()], N: [r.pt(), r.pt()], profiles: 1, cap: [false, false], capsKnown: false, mesh: null };
  r.skip(16 + 1, 'an extrusion');
  if (c.minor >= 1) {
    x.profiles = r.i32();
    if (c.minor >= 2) {
      x.cap = [r.bool(), r.bool()]; x.capsKnown = true;
      if (c.minor >= 3) {
        try {
          const mc = r.versioned('the mesh cache of an extrusion');
          let any = null, analysis = null;
          while (mc.major === 1 && r.u8() === 1) {
            const it = r.versioned('a cached mesh of an extrusion');
            const id = r.uuid(), m = readObject(r, ar, depth + 1);
            r.close(it);
            if (!m || m.t !== 'mesh') continue;
            if (id === MESH_RENDER) x.mesh = m; else if (id === MESH_ANALYSIS) analysis = m; else any = any || m;
          }
          r.close(mc);
          x.mesh = x.mesh || analysis || any;
        } catch (e) { if (!(e instanceof Flaw)) throw e; }
      }
    }
  }
  r.close(c);
  return x;
}

function readClass(r, ar, id, depth) {
  const name = CLS[id];
  switch (name) {
    case 'ON_Mesh': return readMesh(r);
    case 'ON_Brep': return readBrep(r, ar, depth);
    case 'ON_Extrusion': return readExtrusion(r, ar, depth);
    case 'ON_NurbsCurve': return readNurbsCurve(r);
    case 'ON_NurbsSurface': return readNurbsSurface(r);
    case 'ON_LineCurve': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_LineCurve chunk version ${v.major}.${v.minor} is not known`);
      return { t: 'curve', k: 'line', a: r.pt(), b: r.pt(), dom: r.interval() };
    }
    case 'ON_ArcCurve': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_ArcCurve chunk version ${v.major}.${v.minor} is not known`);
      const plane = readPlane(r), radius = r.f64();
      r.skip(72, 'an arc');
      return { t: 'curve', k: 'arc', plane, radius, ang: r.interval(), dom: r.interval() };
    }
    case 'ON_PolylineCurve': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_PolylineCurve chunk version ${v.major}.${v.minor} is not known`);
      const n = r.count(24, 'the points of a polyline'), P = r.f64s(3 * n);
      const nt = r.count(8, 'the parameters of a polyline'), T = r.f64s(nt);
      if (n < 2 || nt !== n) throw new Flaw(`ON_PolylineCurve with ${n} points and ${nt} parameters`);
      return { t: 'curve', k: 'polyline', P, T, dom: [T[0], T[n - 1]] };
    }
    case 'ON_PolyCurve': return readPolyCurve(r, ar, depth);
    case 'ON_PlaneSurface': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_PlaneSurface chunk version ${v.major}.${v.minor} is not known`);
      const s = { t: 'surface', k: 'plane', plane: readPlane(r), dom: [r.interval(), r.interval()], ext: null };
      s.ext = v.minor >= 1 ? [r.interval(), r.interval()] : s.dom;
      return s;
    }
    case 'ON_RevSurface': {
      const v = r.ver();
      if (v.major !== 1 && v.major !== 2) throw new Flaw(`ON_RevSurface chunk version ${v.major}.${v.minor} is not known`);
      const s = { t: 'surface', k: 'rev', axis: [r.pt(), r.pt()], ang: r.interval(), par: null, transposed: false, curve: null };
      s.par = v.major === 2 ? r.interval() : [Math.min(s.ang[0], s.ang[1]), Math.max(s.ang[0], s.ang[1])];
      r.skip(48);
      s.transposed = r.i32() !== 0;
      if (r.u8()) s.curve = readObject(r, ar, depth + 1);
      return s;
    }
    case 'ON_SumSurface': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_SumSurface chunk version ${v.major}.${v.minor} is not known`);
      const base = r.pt();
      r.skip(48);
      return { t: 'surface', k: 'sum', base, curves: [readObject(r, ar, depth + 1), readObject(r, ar, depth + 1)] };
    }
    case 'ON_Point': { const v = r.ver(); if (v.major !== 1) throw new Flaw(`ON_Point chunk version ${v.major}.${v.minor} is not known`); return { t: 'point', p: r.pt() }; }
    case 'ON_PointCloud': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_PointCloud chunk version ${v.major}.${v.minor} is not known`);
      const n = r.count(24, 'the points of a point cloud');
      return { t: 'cloud', P: r.f64s(3 * n) };
    }
    case 'ON_InstanceRef': {
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`ON_InstanceRef chunk version ${v.major}.${v.minor} is not known`);
      const def = r.uuid(), x = r.f64s(16), M = new Float64Array(16);
      // stored row by row
      for (let row = 0; row < 4; row++) for (let col = 0; col < 4; col++) M[col * 4 + row] = x[row * 4 + col];
      return { t: 'iref', def, M };
    }
    case 'ON_AnnotationArrow (Rhino 2)': {
      // an arrow of Rhino 2: a line from tail to head (openNURBS reads it as a line curve too)
      const v = r.ver();
      if (v.major !== 1) throw new Flaw(`annotation arrow of chunk version ${v.major}.${v.minor}`);
      return { t: 'curve', k: 'line', a: r.pt(), b: r.pt(), dom: [0, 1] };
    }
    case 'ON_InstanceDefinition': return readInstanceDefinition(r);
    case 'ON_Layer': return readLayer(r);
    default: return { t: 'other', cls: name || `class ${id}` };
  }
}

/**
 * One polymorphic object: class chunk { class id, class data, user data …, end mark }. Class data that cannot be
 * decoded gives { t: 'bad', cls, why } and the archive is read on behind the chunk.
 */
function readObject(r, ar, depth) {
  if (depth > MAX_DEPTH) throw new Flaw(`objects are nested more than ${MAX_DEPTH} deep`);
  const c = r.expect(TC.CLASS, 'an object');
  try {
    const u = r.expect(TC.CLASS_UUID, 'an object');
    const id = r.uuid();
    r.close(u);
    if (id === NIL) return null;
    const d = r.expect(TC.CLASS_DATA, 'an object');
    let g;
    try { g = readClass(r, ar, id, depth); } catch (e) {
      if (!(e instanceof Flaw)) { if (e instanceof RangeError || e instanceof TypeError) g = { t: 'bad', cls: CLS[id] || `class ${id}`, why: `its data could not be decoded (${e.message})` }; else throw e; } else g = { t: 'bad', cls: CLS[id] || `class ${id}`, why: e.message };
    } finally { r.close(d); }
    // user data: double-precision vertices of Rhino 5 meshes, cached meshes of Rhino 5 extrusions
    const wantUd = (g.t === 'mesh' && r.v <= 50) || (g.t === 'extrusion' && !g.mesh);
    while (wantUd && r.pos + 4 + r.L <= c.end) {
      const ud = r.open();
      try { if (ud.tc === TC.CLASS_USERDATA) readUserData(r, ar, g, depth); } catch (e) { if (!(e instanceof Flaw)) throw e; } finally { r.close(ud); }
      if (ud.tc === TC.CLASS_END) break;
    }
    return g;
  } finally { r.close(c); }
}

function readUserData(r, ar, g, depth) {
  const v = r.ver();
  if (v.major !== 2) return;
  const h = r.expect(TC.CLASS_UD_HEADER, 'user data');
  const cls = r.uuid();
  r.close(h);
  if (cls !== UD_DOUBLE_VERTICES && cls !== UD_EXTRUSION_MESHES) return;
  const a = r.expect(TC.ANON, 'user data');
  if (cls === UD_DOUBLE_VERTICES && g.t === 'mesh') {
    const c = r.versioned('the double-precision vertices of a mesh');
    const nfl = r.i32(), nd = r.i32();
    r.skip(8);
    const n = r.count(24, 'the double-precision vertices of a mesh');
    if (n === g.nv && nd === n && nfl === n && n) g.vd = { size: 24 * n, crc: null, raw: r.b.subarray(r.pos, r.pos + 24 * n) };
    r.close(c);
  } else if (cls === UD_EXTRUSION_MESHES && g.t === 'extrusion') {
    for (let i = 0; i < 2 && !g.mesh; i++) { const m = readObject(r, ar, depth + 1); if (m && m.t === 'mesh') g.mesh = m; }
  }
  r.close(a);
}

// ---------------------------------------------------------------------------------------------------------
// Mesh vertices
// ---------------------------------------------------------------------------------------------------------

const inflated = (buf, out, what, notes) => {
  if (out.length !== buf.size) throw new Flaw(`the ${what} buffer inflates to ${out.length} bytes instead of ${buf.size}`);
  if (buf.crc !== null && crc32(out) !== buf.crc) notes.crc++;
  return out;
};
function floatsOf(bytes, Type) {
  const n = bytes.length / Type.BYTES_PER_ELEMENT;
  if (bytes.byteOffset % Type.BYTES_PER_ELEMENT === 0) return new Type(bytes.buffer, bytes.byteOffset, n);
  return new Type(bytes.slice().buffer, 0, n);
}

/**
 * Vertices of a mesh decoded by readRhino3dm, as a Float64Array (double-precision vertices) or Float32Array.
 * The buffers are inflated on first use.
 */
export function meshVertices(m, notes = { crc: 0 }) {
  if (m.V) return m.V;
  const get = (buf, what) => (buf.raw ? buf.raw : buf.out ? buf.out : inflated(buf, inflateZlibSync(buf.zip, buf.size), what, notes));
  let V = null;
  if (m.vd) {
    try { const d = floatsOf(get(m.vd, 'double-precision vertex'), Float64Array); if (d.length === 3 * m.nv) { V = d; m.precision = 'double'; } } catch (e) { if (!(e instanceof Flaw) && !/[Dd]eflate|zlib|Huffman/.test(e.message)) throw e; }
  }
  if (!V) {
    if (!m.vf) V = new Float32Array(0);
    else {
      let bytes;
      try { bytes = get(m.vf, 'vertex'); } catch (e) { if (e instanceof Flaw) throw e; throw new Flaw(`the vertex buffer of a mesh could not be inflated (${e.message})`); }
      V = floatsOf(bytes, Float32Array);
    }
  }
  m.V = V;
  return V;
}

/** Inflate the larger vertex buffers of the meshes with the platform's decompressor. */
async function inflateMeshes(meshes, notes) {
  const jobs = [];
  for (const m of meshes) {
    if (m.V) continue;
    const buf = m.vd || m.vf;
    if (!buf || !buf.zip || buf.out || buf.zip.length < 16384) continue;
    jobs.push(inflateZlib(buf.zip.slice(), buf.size).then((out) => { if (out.length === buf.size) { if (buf.crc !== null && crc32(out) !== buf.crc) notes.crc++; buf.out = out; } }, () => {}));
    if (jobs.length >= 64) await Promise.all(jobs.splice(0));
  }
  await Promise.all(jobs);
}

// ---------------------------------------------------------------------------------------------------------
// Curves
// ---------------------------------------------------------------------------------------------------------

/** Knots and control points of a NURBS curve as nurbs.js takes them (full knot vector, weights apart). */
function nurbsArrays(c) {
  if (c.arrays) return c.arrays;
  const { dim, rat, order, n, knots, cv } = c, size = dim + rat;
  const U = new Float64Array(knots.length + 2);
  U.set(knots, 1); U[0] = knots[0]; U[U.length - 1] = knots[knots.length - 1];
  const P = new Float64Array(3 * n), W = rat ? new Float64Array(n) : null;
  for (let i = 0; i < n; i++) {
    const w = rat ? cv[i * size + dim] : 1, iw = w !== 0 ? 1 / w : 1;
    for (let k = 0; k < 3 && k < dim; k++) P[3 * i + k] = cv[i * size + k] * iw;
    if (W) W[i] = w;
  }
  return (c.arrays = { p: order - 1, U, P, W, n });
}

/**
 * Points of a curve between the parameters t0 and t1 (in that order along the result).
 * @param {(radius: number) => number} arcStep angular step [rad] for an arc of the given radius
 * @returns {number[][]|null}
 */
function curvePoints(c, t0, t1, arcStep, depth = 0) {
  if (!c || c.t !== 'curve' || depth > MAX_DEPTH) return null;
  const [d0, d1] = c.dom;
  if (t0 === undefined) { t0 = d0; t1 = d1; }
  const rev = t0 > t1, lo = Math.max(Math.min(t0, t1), Math.min(d0, d1)), hi = Math.min(Math.max(t0, t1), Math.max(d0, d1));
  if (!(hi >= lo)) return null;
  let out = null;
  const span = d1 - d0;
  switch (c.k) {
    case 'line': {
      const at = (t) => { const f = span !== 0 ? (t - d0) / span : 0; return [c.a[0] + f * (c.b[0] - c.a[0]), c.a[1] + f * (c.b[1] - c.a[1]), c.a[2] + f * (c.b[2] - c.a[2])]; };
      out = [at(lo), at(hi)];
      break;
    }
    case 'arc': {
      if (!(span !== 0)) return null;
      const [a0, a1] = c.ang, ang = (t) => a0 + ((t - d0) / span) * (a1 - a0);
      const s0 = ang(lo), sweep = ang(hi) - s0;
      const flat = arcPoints(c.plane.o, c.plane.x, c.plane.y, c.radius, s0, sweep, (arcStep(Math.abs(c.radius)) * 180) / Math.PI);
      out = [];
      for (let i = 0; i < flat.length; i += 3) out.push([flat[i], flat[i + 1], flat[i + 2]]);
      break;
    }
    case 'polyline': {
      const { P, T } = c, n = T.length;
      const at = (t) => {
        let i = 0;
        while (i + 2 < n && T[i + 1] <= t) i++;
        const f = T[i + 1] !== T[i] ? (t - T[i]) / (T[i + 1] - T[i]) : 0;
        return [P[3 * i] + f * (P[3 * i + 3] - P[3 * i]), P[3 * i + 1] + f * (P[3 * i + 4] - P[3 * i + 1]), P[3 * i + 2] + f * (P[3 * i + 5] - P[3 * i + 2])];
      };
      out = [at(lo)];
      for (let i = 0; i < n; i++) if (T[i] > lo && T[i] < hi) out.push([P[3 * i], P[3 * i + 1], P[3 * i + 2]]);
      out.push(at(hi));
      break;
    }
    case 'nurbs': {
      const a = nurbsArrays(c);
      const flat = sampleCurve(a.p, a.U, a.P, a.W, a.p === 1 ? 1 : a.p === 2 ? 6 : 8, lo, hi);
      if (!flat) return null;
      out = [];
      for (let i = 0; i < flat.length; i += 3) out.push([flat[i], flat[i + 1], flat[i + 2]]);
      break;
    }
    case 'poly': {
      out = [];
      for (let i = 0; i < c.segs.length; i++) {
        const a = c.T[i], b = c.T[i + 1], s = c.segs[i];
        if (!(b > a) || hi <= a || lo >= b) continue;
        // the segment keeps its own parameters; the polycurve maps [a, b] onto them linearly
        const map = (t) => s.dom[0] + ((t - a) / (b - a)) * (s.dom[1] - s.dom[0]);
        const part = curvePoints(s, map(Math.max(lo, a)), map(Math.min(hi, b)), arcStep, depth + 1);
        if (!part) return null;
        for (let k = out.length ? 1 : 0; k < part.length; k++) out.push(part[k]);
      }
      if (out.length < 2) return null;
      break;
    }
    default: return null;
  }
  for (const p of out) if (!(Number.isFinite(p[0]) && Number.isFinite(p[1]) && Number.isFinite(p[2]))) return null;
  return rev ? out.reverse() : out;
}

/** Parameters between lo and hi at which a curve may have a corner or changes its piece (ends included). */
function curveBreaks(c, lo, hi, depth = 0) {
  const out = [lo];
  const add = (t) => { if (t > out[out.length - 1] && t < hi) out.push(t); };
  if (c.k === 'polyline') for (const t of c.T) add(t);
  else if (c.k === 'nurbs') {
    // knots of full multiplicity (possible corners), and every fourth span so that no piece is long in parameter
    const K = c.knots, p = c.order - 1;
    for (let i = p, spans = 0; i < c.n - 1; i++) {
      if (!(K[i] > K[i - 1])) continue;
      let mult = 1;
      while (i + mult < K.length && K[i + mult] === K[i]) mult++;
      if (mult >= p || ++spans % 4 === 0) add(K[i]);
    }
  }
  else if (c.k === 'arc') for (let i = 1; i < 4; i++) add(c.dom[0] + ((c.dom[1] - c.dom[0]) * i) / 4);
  else if (c.k === 'poly' && depth < MAX_DEPTH) {
    for (let i = 0; i < c.segs.length; i++) {
      const a = c.T[i], b = c.T[i + 1], s = c.segs[i];
      if (!(b > a) || hi <= a || lo >= b) continue;
      const to = (t) => s.dom[0] + ((t - a) / (b - a)) * (s.dom[1] - s.dom[0]), back = (ts) => a + ((ts - s.dom[0]) / (s.dom[1] - s.dom[0] || 1)) * (b - a);
      const s0 = to(Math.max(lo, a)), s1 = to(Math.min(hi, b));
      for (const ts of curveBreaks(s, Math.min(s0, s1), Math.max(s0, s1), depth + 1)) add(back(ts));
      add(b);
    }
  }
  out.push(hi);
  return out;
}

/**
 * Points of a curve from t0 to t1 for the boundary of a face: between two points the curve stays within `tol` of
 * the chord and turns by no more than `maxTurn` [rad].
 */
function edgePolyline(c, t0, t1, tol, maxTurn) {
  if (!c || c.t !== 'curve') return null;
  const rev = t0 > t1, lo = Math.max(Math.min(t0, t1), Math.min(c.dom[0], c.dom[1])), hi = Math.min(Math.max(t0, t1), Math.max(c.dom[0], c.dom[1]));
  if (!(hi > lo)) return null;
  const out = [curveAt(c, lo)], cosTurn = Math.cos(maxTurn);
  // a piece is flat when the curve at 1/4, 1/2 and 3/4 of it lies within `tol` of the chord, in order along it,
  // and the chords to and from the middle turn by no more than `maxTurn`
  const walk = (a, pa, b, pb, depth) => {
    const m = (a + b) / 2, pm = curveAt(c, m);
    const ex = pb[0] - pa[0], ey = pb[1] - pa[1], ez = pb[2] - pa[2], l2 = ex * ex + ey * ey + ez * ez;
    const probe = (q, f0, f1) => { const dx = q[0] - pa[0], dy = q[1] - pa[1], dz = q[2] - pa[2], f = l2 > 0 ? (dx * ex + dy * ey + dz * ez) / l2 : -1; return f > f0 && f < f1 && Math.hypot(dx - f * ex, dy - f * ey, dz - f * ez) <= tol; };
    let flat = depth >= 12 || !(pm[0] === pm[0]);
    if (!flat && probe(pm, 0.2, 0.8)) {
      const dx = pm[0] - pa[0], dy = pm[1] - pa[1], dz = pm[2] - pa[2], vx = pb[0] - pm[0], vy = pb[1] - pm[1], vz = pb[2] - pm[2], l = Math.hypot(dx, dy, dz) * Math.hypot(vx, vy, vz);
      flat = (l === 0 || (dx * vx + dy * vy + dz * vz) / l >= cosTurn) && probe(curveAt(c, (3 * a + b) / 4), 0, 0.5) && probe(curveAt(c, (a + 3 * b) / 4), 0.5, 1);
    }
    if (flat) { out.push(pb); return; }
    walk(a, pa, m, pm, depth + 1);
    walk(m, pm, b, pb, depth + 1);
  };
  const br = curveBreaks(c, lo, hi);
  let pa = out[0];
  for (let i = 1; i < br.length; i++) {
    const pb = curveAt(c, br[i]);
    // a straight piece needs no point in between
    if (c.k === 'line' || (c.k === 'nurbs' && c.order === 2) || c.k === 'polyline') out.push(pb); else walk(br[i - 1], pa, br[i], pb, 0);
    pa = pb;
  }
  for (const p of out) if (!(Number.isFinite(p[0]) && Number.isFinite(p[1]) && Number.isFinite(p[2]))) return null;
  return rev ? out.reverse() : out;
}

/** Point of a curve at the parameter t. */
function curveAt(c, t, depth = 0) {
  const [d0, d1] = c.dom, span = d1 - d0;
  switch (c.k) {
    case 'line': { const f = span !== 0 ? (t - d0) / span : 0; return [c.a[0] + f * (c.b[0] - c.a[0]), c.a[1] + f * (c.b[1] - c.a[1]), c.a[2] + f * (c.b[2] - c.a[2])]; }
    case 'arc': {
      const ang = c.ang[0] + (span !== 0 ? (t - d0) / span : 0) * (c.ang[1] - c.ang[0]), x = c.radius * Math.cos(ang), y = c.radius * Math.sin(ang), pl = c.plane;
      return [pl.o[0] + x * pl.x[0] + y * pl.y[0], pl.o[1] + x * pl.x[1] + y * pl.y[1], pl.o[2] + x * pl.x[2] + y * pl.y[2]];
    }
    case 'polyline': {
      const { P, T } = c, n = T.length;
      let lo = 0, hi = n - 2;
      while (lo < hi) { const m = (lo + hi + 1) >> 1; if (T[m] <= t) lo = m; else hi = m - 1; }
      const f = T[lo + 1] !== T[lo] ? (t - T[lo]) / (T[lo + 1] - T[lo]) : 0;
      return [P[3 * lo] + f * (P[3 * lo + 3] - P[3 * lo]), P[3 * lo + 1] + f * (P[3 * lo + 4] - P[3 * lo + 1]), P[3 * lo + 2] + f * (P[3 * lo + 5] - P[3 * lo + 2])];
    }
    case 'nurbs': {
      const a = nurbsArrays(c), p = a.p, U = a.U, n = a.n;
      if (!a.N) { a.N = new Float64Array(p + 1); a.L = new Float64Array(p + 1); a.R = new Float64Array(p + 1); }
      const u = Math.max(U[p], Math.min(t, U[n])), s = findSpan(n - 1, p, u, U);
      basisFuns(s, u, p, U, a.N, a.L, a.R);
      let x = 0, y = 0, z = 0, w = 0;
      for (let j = 0; j <= p; j++) { const i = s - p + j, k = a.W ? a.N[j] * a.W[i] : a.N[j]; x += k * a.P[3 * i]; y += k * a.P[3 * i + 1]; z += k * a.P[3 * i + 2]; w += k; }
      if (!a.W || w === 0) w = 1;
      return [x / w, y / w, z / w];
    }
    case 'poly': {
      const T = c.T, n = c.segs.length;
      if (!n || depth > MAX_DEPTH) return [NaN, NaN, NaN];
      let i = 0;
      while (i + 1 < n && T[i + 1] <= t) i++;
      const s = c.segs[i], a = T[i], b = T[i + 1];
      return curveAt(s, s.dom[0] + (b !== a ? (t - a) / (b - a) : 0) * (s.dom[1] - s.dom[0]), depth + 1);
    }
    default: return [NaN, NaN, NaN];
  }
}

// ---------------------------------------------------------------------------------------------------------
// Surfaces
// ---------------------------------------------------------------------------------------------------------

/**
 * Evaluator of a surface: { at(u, v), dom: [[u0, u1], [v0, v1]], spans: [nu, nv], how } in the parameters of
 * the file, so that the 2-D trimming curves of a B-rep apply as they are; or { why } when it cannot be evaluated.
 */
function surfaceEval(s) {
  if (!s) return { why: 'surface missing from the B-rep' };
  if (s.t !== 'surface') return { why: `surface of class ${s.cls || s.t}${s.why ? ` (${s.why})` : ''}` };
  if (s.ev) return s.ev;
  const spansOf = (K, a, b) => { let n = 0; for (let i = a; i < b; i++) if (K[i + 1] > K[i]) n++; return Math.max(1, n); };
  const curveSpans = (c) => (c.k === 'nurbs' ? spansOf(c.knots, c.order - 2, c.n - 1) : c.k === 'polyline' ? c.T.length - 1 : c.k === 'poly' ? c.segs.reduce((n, x) => n + curveSpans(x), 0) : c.k === 'arc' ? 4 : 1);
  switch (s.k) {
    case 'plane': {
      const pl = s.plane, [dx, dy] = s.dom, [ex, ey] = s.ext;
      const map = (t, d, e) => (d[1] !== d[0] ? e[0] + ((t - d[0]) / (d[1] - d[0])) * (e[1] - e[0]) : e[0]);
      return (s.ev = { how: 'plane', flat: true, dom: [dx, dy], spans: [1, 1], at: (u, v) => { const x = map(u, dx, ex), y = map(v, dy, ey); return [pl.o[0] + x * pl.x[0] + y * pl.y[0], pl.o[1] + x * pl.x[1] + y * pl.y[1], pl.o[2] + x * pl.x[2] + y * pl.y[2]]; } });
    }
    case 'nurbs': {
      const { dim, rat, order, n, knots, cv } = s, size = dim + rat, nu = n[0], nv = n[1], pu = order[0] - 1, pv = order[1] - 1;
      if (dim < 3) return (s.ev = { why: `ON_NurbsSurface of dimension ${dim}` });
      const full = (K) => { const U = new Float64Array(K.length + 2); U.set(K, 1); U[0] = K[0]; U[U.length - 1] = K[K.length - 1]; return U; };
      const U = full(knots[0]), V = full(knots[1]);
      if (!(U[nu] > U[pu]) || !(V[nv] > V[pv])) return (s.ev = { why: 'ON_NurbsSurface with an empty parameter range' });
      const Nu = new Float64Array(pu + 1), Nv = new Float64Array(pv + 1), L = new Float64Array(Math.max(pu, pv) + 1), R = new Float64Array(Math.max(pu, pv) + 1);
      const at = (u, v) => {
        u = Math.max(U[pu], Math.min(u, U[nu])); v = Math.max(V[pv], Math.min(v, V[nv]));
        const su = findSpan(nu - 1, pu, u, U);
        basisFuns(su, u, pu, U, Nu, L, R);
        const sv = findSpan(nv - 1, pv, v, V);
        basisFuns(sv, v, pv, V, Nv, L, R);
        let x = 0, y = 0, z = 0, w = 0;
        for (let a = 0; a <= pu; a++) for (let b = 0; b <= pv; b++) {
          // control points are stored with v running fastest, in homogeneous form when rational
          const i = ((su - pu + a) * nv + sv - pv + b) * size, k = Nu[a] * Nv[b];
          x += k * cv[i]; y += k * cv[i + 1]; z += k * cv[i + 2]; w += rat ? k * cv[i + dim] : k;
        }
        if (w === 0) w = 1;
        return [x / w, y / w, z / w];
      };
      return (s.ev = { how: 'NURBS surface', at, dom: [[U[pu], U[nu]], [V[pv], V[nv]]], spans: [spansOf(U, pu, nu), spansOf(V, pv, nv)] });
    }
    case 'rev': {
      const c = s.curve;
      if (!c || c.t !== 'curve') return (s.ev = { why: 'surface of revolution without a profile curve' });
      const A = s.axis[0], D = V3.norm(V3.sub(s.axis[1], s.axis[0])), [a0, a1] = s.ang, [t0, t1] = s.par;
      if (!(V3.len(D) > 0) || !(a1 !== a0) || !(t1 !== t0)) return (s.ev = { why: 'surface of revolution without an axis or an angle' });
      const rot = (ta, tc) => {
        const p = curveAt(c, tc), ang = a0 + ((ta - t0) / (t1 - t0)) * (a1 - a0), q = V3.sub(p, A), h = V3.dot(q, D);
        const X = [q[0] - h * D[0], q[1] - h * D[1], q[2] - h * D[2]], Y = V3.cross(D, X), cs = Math.cos(ang), sn = Math.sin(ang);
        return [A[0] + h * D[0] + cs * X[0] + sn * Y[0], A[1] + h * D[1] + cs * X[1] + sn * Y[1], A[2] + h * D[2] + cs * X[2] + sn * Y[2]];
      };
      const na = Math.max(1, Math.ceil(Math.abs(a1 - a0) / (Math.PI / 2) - 1e-9)), nc = curveSpans(c);
      return (s.ev = s.transposed ? { how: 'surface of revolution', at: (u, v) => rot(v, u), dom: [c.dom, [t0, t1]], spans: [nc, na] } : { how: 'surface of revolution', at: rot, dom: [[t0, t1], c.dom], spans: [na, nc] });
    }
    case 'sum': {
      const [a, b] = s.curves, o = s.base;
      if (!a || !b || a.t !== 'curve' || b.t !== 'curve') return (s.ev = { why: 'sum surface without its two curves' });
      return (s.ev = { how: 'sum surface', at: (u, v) => { const p = curveAt(a, u), q = curveAt(b, v); return [p[0] + q[0] + o[0], p[1] + q[1] + o[1], p[2] + q[2] + o[2]]; }, dom: [a.dom, b.dom], spans: [curveSpans(a), curveSpans(b)] });
    }
    default: return (s.ev = { why: `surface of kind ${s.k}` });
  }
}

/**
 * Facets a surface needs within the parameter box [u0, u1] × [v0, v1]: from the sag and the turn of sampled
 * iso-lines, as brep.js does for its own spline surfaces. Returns the surface object tessellateFace takes; its
 * parameters are those of the file and `toParam` is filled in per face (the trims give the parameters).
 */
function facetting(ev, box, defl, angStep, cells) {
  const [u0, u1, v0, v1] = box;
  if (ev.flat) return { type: 'spline', inward: false, perU: false, perV: false, poleRows: null, fromParam: ev.at, toParam: null, nu: 1, nv: 1, relax: 1, relaxed: false, scales: () => ({ su: 1 / (u1 - u0 || 1), sv: 1 / (v1 - v0 || 1), stepU: Infinity, stepV: Infinity }) };
  const frac = (d, a, b) => Math.max(1e-3, Math.min(1, Math.abs(b - a) / Math.abs(d[1] - d[0])));
  const gu = Math.max(6, Math.min(40, Math.ceil(4 * ev.spans[0] * frac(ev.dom[0], u0, u1)))), gv = Math.max(6, Math.min(40, Math.ceil(4 * ev.spans[1] * frac(ev.dom[1], v0, v1))));
  const G = [];
  for (let j = 0; j <= gv; j++) for (let i = 0; i <= gu; i++) { const p = ev.at(u0 + ((u1 - u0) * i) / gu, v0 + ((v1 - v0) * j) / gv); if (!(p[0] === p[0] && p[1] === p[1] && p[2] === p[2])) return null; G.push(p); }
  const g = (i, j) => G[i + (gu + 1) * j];
  const need = (along) => {
    const n = along ? gu : gv, m = along ? gv : gu;
    let sag = 0, turn = 0;
    for (let j = 0; j <= m; j++) for (let i = 1; i < n; i++) {
      const a = along ? g(i - 1, j) : g(j, i - 1), b = along ? g(i, j) : g(j, i), c = along ? g(i + 1, j) : g(j, i + 1);
      sag = Math.max(sag, Math.hypot(b[0] - (a[0] + c[0]) / 2, b[1] - (a[1] + c[1]) / 2, b[2] - (a[2] + c[2]) / 2));
      const e1 = V3.sub(b, a), e2 = V3.sub(c, b), l = V3.len(e1) * V3.len(e2);
      if (l > 0 && V3.len(e1) > 1e-3 * defl && V3.len(e2) > 1e-3 * defl) turn = Math.max(turn, Math.acos(Math.max(-1, Math.min(1, V3.dot(e1, e2) / l))));
    }
    return Math.max(1, Math.min(512, Math.ceil(n * Math.max(Math.sqrt(sag / (4 * defl)), turn / angStep) - 1e-9)));
  };
  let nu = need(true), nv = need(false), relax = 1;
  // a small feature of high curvature sets the facet size of the whole face: the number of facets is capped
  if (nu * nv > cells) { relax = Math.sqrt(cells / (nu * nv)); nu = Math.max(1, Math.round(nu * relax)); nv = Math.max(1, Math.round(nv * relax)); }
  const stepU = Math.abs(u1 - u0) / nu, stepV = Math.abs(v1 - v0) / nv;
  return { type: 'spline', inward: false, perU: false, perV: false, poleRows: null, fromParam: ev.at, toParam: null, nu, nv, relax, relaxed: relax < 1, scales: () => ({ su: 1 / stepU, sv: 1 / stepV, stepU, stepV }) };
}

/**
 * Extent of a surface within a parameter box: the mean length of its iso-lines in u and in v. The smaller of the
 * two bounds the chord error a face may have: facets that are coarse against the width of a narrow face (a fillet)
 * would stand at an angle to its surface.
 */
function patchSize(ev, box) {
  const [u0, u1, v0, v1] = box, n = 6;
  const G = [];
  for (let j = 0; j <= n; j++) for (let i = 0; i <= n; i++) G.push(ev.at(u0 + ((u1 - u0) * i) / n, v0 + ((v1 - v0) * j) / n));
  const d = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  let lu = 0, lv = 0;
  for (let j = 0; j <= n; j++) for (let i = 0; i < n; i++) { lu += d(G[i + (n + 1) * j], G[i + 1 + (n + 1) * j]); lv += d(G[j + (n + 1) * i], G[j + (n + 1) * (i + 1)]); }
  lu /= n + 1; lv /= n + 1;
  return lu === lu && lv === lv ? [lu, lv] : [0, 0];
}
// chord error of a face as a fraction of its smaller extent; of a boundary chord as a fraction of the facet size
const FACE_SAG = 0.02, EDGE_SAG = 0.1;

/** Triangulate a face given by loops of { uv, p } in the parameter plane of its surface. */
function liftFace(surf, rings) {
  const uv = new Map();
  for (const ring of rings) for (const q of ring) uv.set(q.p, q.uv);
  surf.toParam = (p) => { const q = uv.get(p); return q ? [q[0], q[1], 0] : [NaN, NaN, 0]; };
  try { return tessellateFace(surf, rings.map((ring) => ring.map((q) => q.p)), false); } catch (e) { return `triangulation failed (${e.message})`; }
}

/** The whole patch of a surface as one face. */
function patchMesh(ev, defl, angStep, diag, name, cells) {
  const [[u0, u1], [v0, v1]] = ev.dom, size = Math.min(...patchSize(ev, [u0, u1, v0, v1]));
  if (size > 0) defl = Math.min(defl, FACE_SAG * size);
  const surf = facetting(ev, [u0, u1, v0, v1], defl, angStep, cells);
  if (!surf) return 'the surface cannot be evaluated';
  const ring = [], add = (u, v) => ring.push({ uv: [u, v], p: ev.at(u, v) });
  const { nu, nv } = surf;
  for (let i = 0; i < nu; i++) add(u0 + ((u1 - u0) * i) / nu, v0);
  for (let j = 0; j < nv; j++) add(u1, v0 + ((v1 - v0) * j) / nv);
  for (let i = nu; i > 0; i--) add(u0 + ((u1 - u0) * i) / nu, v1);
  for (let j = nv; j > 0; j--) add(u0, v0 + ((v1 - v0) * j) / nv);
  const res = liftFace(surf, [ring]);
  if (typeof res === 'string') return res;
  const bm = new BodyMesh(diag);
  if (!bm.addFace(res.pts, res.indices, false, name)) return 'the surface is degenerate';
  return { mesh: bm.finish(null), relaxed: surf.relaxed };
}

// ---------------------------------------------------------------------------------------------------------
// B-rep topology and tessellation
// ---------------------------------------------------------------------------------------------------------

/** Decode the curves, surfaces and topology of a B-rep located by readBrep. */
function brepTopology(b, ar) {
  if (b.topo) return b.topo;
  const r = ar.reader, save = [r.pos, r.lim];
  const within = (s, fn) => { r.pos = s.at; r.lim = s.end; return fn(); };
  const records = (s, what, size, one) => within(s, () => {
    const c = r.expect(TC.ANON, what);
    const v = r.ver();
    if (v.major !== 1) throw new Flaw(`${what} of chunk version ${v.major}.${v.minor}`);
    const n = r.count(size, what), out = [];
    for (let i = 0; i < n; i++) out.push(one());
    r.close(c);
    return out;
  });
  const ints = (what) => { const n = r.count(4, what), out = []; for (let i = 0; i < n; i++) out.push(r.i32()); return out; };
  const newer = r.v >= 3 && r.onv >= 200206180;
  try {
    const t = {};
    t.C2 = within(b.c2, () => readObjectArray(r, ar, 1, 'the 2-D curves of a B-rep'));
    t.C3 = within(b.c3, () => readObjectArray(r, ar, 1, 'the 3-D curves of a B-rep'));
    t.S = within(b.s, () => readObjectArray(r, ar, 1, 'the surfaces of a B-rep'));
    t.V = records(b.v, 'the vertices of a B-rep', 44, () => { r.skip(4); const p = r.pt(); ints('the edges of a B-rep vertex'); r.skip(8); return p; });
    t.E = records(b.e, 'the edges of a B-rep', 52, () => {
      r.skip(4);
      const e = { c3: r.i32(), rev: r.i32() !== 0, dom: r.interval(), v: [r.i32(), r.i32()], pts: null };
      ints('the trims of a B-rep edge'); r.skip(8);
      if (newer) r.skip(16);
      return e;
    });
    t.T = records(b.tr, 'the trims of a B-rep', 68, () => {
      r.skip(4);
      // 2-D curve and the part of it the trim uses, edge, vertices, direction of the edge against the trim
      const tr = { c2: r.i32(), dom: r.interval(), e: r.i32(), v: [r.i32(), r.i32()], rev: r.i32() !== 0, flip: false };
      r.skip(4 + 4 + 4 + 16);
      if (newer) { r.skip(16); tr.flip = r.u8() === 1; r.skip(7 + 24); } else r.skip(48);
      r.skip(16);
      return tr;
    });
    t.L = records(b.l, 'the loops of a B-rep', 16, () => { r.skip(4); const l = { trims: ints('the trims of a B-rep loop'), type: r.i32() }; r.skip(4); return l; });
    t.F = records(b.f, 'the faces of a B-rep', 16, () => { r.skip(4); const f = { loops: ints('the loops of a B-rep face'), s: r.i32(), rev: r.i32() !== 0 }; r.skip(4); return f; });
    b.topo = t;
  } catch (e) {
    if (!(e instanceof Flaw)) throw e;
    b.topo = { why: e.message };
  } finally { r.pos = save[0]; r.lim = save[1]; }
  return b.topo;
}

/**
 * Triangulate the faces of a B-rep that have no stored mesh. Each face is triangulated in the parameter plane
 * of its surface between its trimming loops; the 3-D points of the boundary are those of the edge curves, which
 * neighbouring faces share, so that the faces meet without gaps.
 * @returns {{ mesh: BodyMesh, read: number, skips: object, hist: object, crossing: number, relaxed: number }}
 */
function tessellateBrep(b, ar, faces, defl, angStep, diag, cells) {
  const t = brepTopology(b, ar);
  const out = { mesh: new BodyMesh(diag), read: 0, skips: {}, hist: {}, crossing: 0, relaxed: 0 };
  if (t.why) { tally(out.skips, `B-rep topology could not be decoded (${t.why})`, faces.length); return out; }
  const fine = () => (2 * Math.PI) / 180;
  // pass 1: the trims of every face in the parameter plane, the chord error each face may have
  const jobs = [], edgeTol = new Map();
  for (const fi of faces) {
    const f = t.F[fi];
    if (!f) { tally(out.skips, 'face missing from the B-rep'); continue; }
    const ev = surfaceEval(f.s >= 0 && f.s < t.S.length ? t.S[f.s] : null);
    if (ev.why) { tally(out.skips, ev.why); continue; }
    const loops = [];
    let why = '';
    const lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
    for (const li of f.loops) {
      const l = t.L[li];
      if (!l || l.type >= 3) continue;
      const parts = [];
      for (const ti of l.trims) {
        const tr = t.T[ti], c2 = tr && tr.c2 >= 0 && tr.c2 < t.C2.length ? t.C2[tr.c2] : null;
        const uv = c2 ? curvePoints(c2, tr.dom[0], tr.dom[1], fine) : null;
        if (!uv || uv.length < 2) { why = 'a trim has no 2-D curve that can be evaluated'; break; }
        if (tr.flip) uv.reverse();
        for (const q of uv) for (let k = 0; k < 2; k++) { if (q[k] < lo[k]) lo[k] = q[k]; if (q[k] > hi[k]) hi[k] = q[k]; }
        parts.push({ tr, uv });
      }
      if (why) break;
      if (parts.length) loops.push(parts);
    }
    if (why) { tally(out.skips, `${ev.how}: ${why}`); continue; }
    if (!loops.length || !(hi[0] > lo[0]) || !(hi[1] > lo[1])) { tally(out.skips, `${ev.how}: the face has no boundary`); continue; }
    const box = [lo[0], hi[0], lo[1], hi[1]], [lu, lv] = patchSize(ev, box), size = Math.min(lu, lv), sag = size > 0 ? Math.min(defl, FACE_SAG * size) : defl;
    const surf = facetting(ev, box, sag, angStep, cells);
    if (!surf) { tally(out.skips, `${ev.how}: the surface cannot be evaluated`); continue; }
    // An edge is sampled once, to the finer of the faces it bounds. Its chords must stay close to the curve
    // against the size of the facets beside it, or the triangles along the boundary would stand at an angle to
    // the surface.
    // Where the facets of a face had to be made larger than asked for, its boundary is sampled as coarsely: fine
    // chords along the boundary of coarse facets would give fans of slender triangles across the curvature.
    const coarse = sag / (4 * surf.relax * surf.relax), tol = ev.flat ? sag / 4 : Math.min(coarse, EDGE_SAG * Math.min(lu / surf.nu, lv / surf.nv) || coarse), turn = Math.min(angStep / surf.relax, 1);
    for (const parts of loops) for (const { tr } of parts) if (tr.e >= 0 && tr.e < t.E.length) { const have = edgeTol.get(tr.e); edgeTol.set(tr.e, have ? { tol: Math.min(have.tol, tol), turn: Math.min(have.turn, turn) } : { tol, turn }); }
    jobs.push({ fi, f, ev, loops, box, surf });
  }
  for (const e of t.E) e.pts = null;
  const edgePts = (ei) => {
    const e = t.E[ei];
    if (e.pts) return e.pts;
    const { tol, turn } = edgeTol.get(ei) ?? { tol: defl / 4, turn: angStep }, va = t.V[e.v[0]], vb = t.V[e.v[1]];
    let pts = e.c3 >= 0 && e.c3 < t.C3.length ? edgePolyline(t.C3[e.c3], e.dom[0], e.dom[1], tol, turn) : null;
    if (pts && e.rev) pts.reverse();
    if (!pts || pts.length < 2) { if (!va || !vb) return (e.pts = []); pts = [va, vb]; }
    // the ends are the vertices of the B-rep, so that neighbouring edges meet exactly
    const near = (p, q) => q && Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) <= 1e-3 * diag;
    if (near(pts[0], va)) pts[0] = va;
    if (near(pts[pts.length - 1], vb)) pts[pts.length - 1] = vb;
    return (e.pts = pts);
  };
  for (const { fi, f, ev, loops, box, surf } of jobs) {
    let why = '';
    const sc = surf.scales(), eu = 1e-9 * (box[1] - box[0]), ev2 = 1e-9 * (box[3] - box[2]);
    // pass 2: 3-D points of the edges, each with the parameters of the nearest point of its trim
    const rings = [];
    for (const parts of loops) {
      const ring = [];
      for (const { tr, uv } of parts) {
        if (tr.e < 0 || tr.e >= t.E.length) {
          // singular trim: a side of the parameter box that is a single point in space
          const a = uv[0], z = uv[uv.length - 1], p = ev.at(a[0], a[1]);
          const k = Math.max(1, Math.ceil(Math.max(Number.isFinite(sc.stepU) ? Math.abs(z[0] - a[0]) / sc.stepU : 0, Number.isFinite(sc.stepV) ? Math.abs(z[1] - a[1]) / sc.stepV : 0) - 1e-9));
          for (let i = 0; i < k; i++) ring.push({ uv: [a[0] + ((z[0] - a[0]) * i) / k, a[1] + ((z[1] - a[1]) * i) / k], p: [p[0], p[1], p[2]] });
          continue;
        }
        const e = edgePts(tr.e), n = e.length;
        if (n < 2) { why = 'an edge has no curve'; break; }
        // the trim as a dense polyline on the surface
        const sub = Math.max(1, Math.ceil((4 * n) / (uv.length - 1))), duv = [], dp = [];
        for (let i = 0; i + 1 < uv.length; i++) for (let k = 0; k < sub; k++) { const q = [uv[i][0] + ((uv[i + 1][0] - uv[i][0]) * k) / sub, uv[i][1] + ((uv[i + 1][1] - uv[i][1]) * k) / sub]; duv.push(q); dp.push(ev.at(q[0], q[1])); }
        duv.push(uv[uv.length - 1]); dp.push(ev.at(uv[uv.length - 1][0], uv[uv.length - 1][1]));
        const m = duv.length - 1, win = Math.max(8, m >> 1);
        let s = 0;
        for (let j = 0; j + 1 < n; j++) {
          const q = tr.rev ? e[n - 1 - j] : e[j];
          let best = duv[0];
          if (j > 0) {
            let bd = Infinity;
            for (let k = s; k < m && k <= s + win; k++) {
              const a = dp[k], c = dp[k + 1], ex = c[0] - a[0], ey = c[1] - a[1], ez = c[2] - a[2], l2 = ex * ex + ey * ey + ez * ez;
              const fr = l2 > 0 ? Math.max(0, Math.min(1, ((q[0] - a[0]) * ex + (q[1] - a[1]) * ey + (q[2] - a[2]) * ez) / l2)) : 0;
              const d = Math.hypot(q[0] - a[0] - fr * ex, q[1] - a[1] - fr * ey, q[2] - a[2] - fr * ez);
              if (d < bd) { bd = d; s = k; best = [duv[k][0] + fr * (duv[k + 1][0] - duv[k][0]), duv[k][1] + fr * (duv[k + 1][1] - duv[k][1])]; }
            }
          }
          ring.push({ uv: best, p: [q[0], q[1], q[2]] });
        }
      }
      if (why) break;
      // drop points that repeat their predecessor in the parameter plane (to rounding)
      const clean = [], same = (a, b) => Math.abs(a.uv[0] - b.uv[0]) <= eu && Math.abs(a.uv[1] - b.uv[1]) <= ev2;
      for (const q of ring) if (!clean.length || !same(clean[clean.length - 1], q)) clean.push(q);
      while (clean.length > 1 && same(clean[0], clean[clean.length - 1])) clean.pop();
      if (clean.length >= 3) rings.push(clean);
    }
    if (why) { tally(out.skips, `${ev.how}: ${why}`); continue; }
    if (!rings.length) { tally(out.skips, `${ev.how}: its boundary is degenerate`); continue; }
    const res = liftFace(surf, rings);
    if (typeof res === 'string') { if (res === CROSSING) out.crossing++; tally(out.skips, `${ev.how}: ${res}`); continue; }
    if (!out.mesh.addFace(res.pts, res.indices, f.rev, `face ${fi + 1}`)) { tally(out.skips, `${ev.how}: its boundary is degenerate`); continue; }
    out.read++;
    if (surf.relaxed) out.relaxed++;
    tally(out.hist, ev.how);
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------------------------------------

const det3 = (M) => M[0] * (M[5] * M[10] - M[9] * M[6]) - M[4] * (M[1] * M[10] - M[9] * M[2]) + M[8] * (M[1] * M[6] - M[5] * M[2]);

/**
 * Rhino 3DM (openNURBS archive), versions 2 to 8.
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ units?: string, scale?: number, hidden?: boolean, deflection?: number, angularStep?: number, curves?: boolean, renderMeshes?: boolean }} [opts]
 *   `units` / `scale` override the unit system of the file; `hidden: true` keeps hidden objects and the objects of
 *   hidden layers (skipped and counted by default); `deflection`: largest distance [m] between the triangles and
 *   a surface that has to be triangulated here (default 1 % of the diagonal of the object; a face is never given
 *   more than 2 % of its smaller extent), `angularStep`: largest
 *   angle one facet spans [deg], default 15 — given either, a face gets as many facets as that asks for; by default
 *   the number per face is bounded; `curves: false` leaves curve objects out; `renderMeshes: false` ignores the
 *   meshes stored with B-rep faces and extrusions and triangulates their surfaces here instead.
 */
export async function parseRhino3dm(input, opts = {}) {
  const fname = input.name || '3DM';
  if (opts.deflection !== undefined && !(opts.deflection > 0)) throw new RangeError('deflection: must be a positive length in metres.');
  if (opts.angularStep !== undefined && !(opts.angularStep > 0 && opts.angularStep <= 90)) throw new RangeError('angularStep: must be within 0…90 degrees.');
  const ar = readRhino3dm(inputBytes(input), fname);
  const warnings = [], skipped = {}, hist = {}, classes = {};
  const notes = { crc: 0 };

  // ---- units ----
  const u = UNITS[ar.unitSystem];
  let units;
  if (opts.units || (Number.isFinite(opts.scale) && opts.scale > 0)) units = resolveUnits(opts, null, 'm');
  else if (typeof u === 'string') units = resolveUnits(opts, u, 'm');
  else if (Array.isArray(u)) units = { unit: 'custom', scale: u[0], basis: 'file', name: u[1] };
  else if (ar.unitSystem === 11 && ar.metresPerUnit > 0 && Number.isFinite(ar.metresPerUnit)) units = { unit: 'custom', scale: ar.metresPerUnit, basis: 'file', name: cleanName(ar.unitName, 'custom unit') };
  else units = resolveUnits(opts, null, 'm');
  const S = units.scale;

  // ---- layers ----
  const byId = new Map(), byIndex = new Map();
  ar.layers.forEach((l, i) => { if (l.id !== NIL) byId.set(l.id, l); byIndex.set(l.index, l); if (!byIndex.has(i)) byIndex.set(i, l); });
  const pathOf = (l, depth = 0) => {
    if (l.path) return l.path;
    const p = l.parent !== NIL && depth < 64 ? byId.get(l.parent) : null;
    return (l.path = cleanName(p && p !== l ? `${pathOf(p, depth + 1)}::${l.name}` : l.name, `layer ${l.index}`));
  };
  const layerOf = (i) => byIndex.get(i) || null;
  const layerName = (i) => { const l = layerOf(i); return l ? pathOf(l) : `layer ${i}`; };

  // ---- which objects are drawn ----
  const byObjectId = new Map();
  for (const o of ar.objects) { tally(classes, classOf(o.geom)); if (o.attrs.id !== NIL) byObjectId.set(o.attrs.id, o); }
  const members = new Set();
  for (const d of ar.idefs.values()) for (const id of d.members) members.add(id);
  let hidden = 0, hiddenLayer = 0;
  const top = [];
  for (const o of ar.objects) {
    if (o.attrs.mode === 3 || members.has(o.attrs.id)) continue;
    if (!opts.hidden) {
      if (!o.attrs.visible || o.attrs.mode === 1) { hidden++; continue; }
      const l = layerOf(o.attrs.layer);
      if (l && !l.visible) { hiddenLayer++; continue; }
    }
    top.push(o);
  }

  // ---- meshes to inflate: those of the objects drawn, through block instances ----
  const wanted = new Set();
  const stored = opts.renderMeshes !== false;
  const meshesOf = (g) => (g.t === 'mesh' ? [g] : !stored ? [] : g.t === 'brep' ? g.meshes.filter(Boolean) : g.t === 'extrusion' && g.mesh ? [g.mesh] : []);
  const seenDefs = new Set();
  let toTessellate = 0;
  const collect = (g) => {
    for (const m of meshesOf(g)) wanted.add(m);
    if (g.t === 'brep') { for (const m of g.meshes) if (!stored || !m || !m.nf) toTessellate++; } else if (g.t === 'surface') toTessellate++;
    if (g.t === 'iref' && !seenDefs.has(g.def)) { seenDefs.add(g.def); const d = ar.idefs.get(g.def); if (d) for (const id of d.members) { const o = byObjectId.get(id); if (o) collect(o.geom); } }
  };
  for (const o of top) collect(o.geom);
  await inflateMeshes(wanted, notes);
  // a facet size given by the caller is kept (up to 65536 facets per face); the default one is bounded per file
  const cells = opts.deflection !== undefined || opts.angularStep !== undefined ? 65536 : Math.max(MIN_CELLS, Math.min(MAX_CELLS, Math.floor(CELL_BUDGET / Math.max(1, toTessellate))));

  // ---- collectors, one surface per layer ----
  const surfaces = new Map(), pointSets = new Map(), polylines = [], clouds = [];
  const surfaceOf = (layer) => { let s = surfaces.get(layer); if (!s) { s = { coords: new Grow(Float64Array, 3 * 4096), buckets: new Map() }; surfaces.set(layer, s); } return s; };
  const bucketOf = (s, name) => { let g = s.buckets.get(name); if (!g) { g = new Grow(Uint32Array, 3072); s.buckets.set(name, g); } return g; };
  /** Append triangles (indices into `V`, three coordinates per vertex) under a transform. */
  const addTriangles = (layer, group, V, T, M) => {
    if (!T.length) return;
    const s = surfaceOf(layer), g = bucketOf(s, group), base = s.coords.n / 3, mirror = det3(M) < 0;
    s.coords.ensure(V.length);
    for (let i = 0; i + 2 < V.length; i += 3) { M4.apply(M, V[i], V[i + 1], V[i + 2], s.coords.a, s.coords.n); s.coords.n += 3; }
    g.ensure(T.length);
    for (let i = 0; i + 2 < T.length; i += 3) { g.a[g.n++] = base + T[i]; g.a[g.n++] = base + T[i + (mirror ? 2 : 1)]; g.a[g.n++] = base + T[i + (mirror ? 1 : 2)]; }
  };
  /** Triangles of a stored mesh: quads split along the shorter diagonal, faces with invalid indices dropped. */
  const meshTriangles = (m, weld) => {
    const V = meshVertices(m, notes), nv = V.length / 3, F = m.F, T = [];
    if (m.precision === 'double') tally(hist, 'meshes with double-precision vertices');
    let ids = null;
    if (weld) { ids = new Uint32Array(nv); for (let i = 0; i < nv; i++) ids[i] = weld.add(V[3 * i], V[3 * i + 1], V[3 * i + 2]); }
    const q = [0, 0, 0, 0];
    let nt = 0, nq = 0, bad = 0;
    for (let f = 0; f + 3 < F.length; f += 4) {
      const a = F[f], b = F[f + 1], c = F[f + 2], d = F[f + 3];
      if (a >= nv || b >= nv || c >= nv || d >= nv) { bad++; continue; }
      if (c === d) { nt++; T.push(a, b, c); } else { nq++; q[0] = a; q[1] = b; q[2] = c; q[3] = d; pushPolygon(T, q, V); }
    }
    if (nt) tally(hist, 'tri3', nt);
    if (nq) tally(hist, 'quad4', nq);
    if (bad) tally(skipped, 'mesh faces with invalid vertex indices', bad);
    if (ids) for (let i = 0; i < T.length; i++) T[i] = ids[T[i]];
    return { V, T };
  };
  const angStep = ((opts.angularStep ?? 15) * Math.PI) / 180;
  const curveStep = () => (5 * Math.PI) / 180;
  const tess = { faces: 0, facesRead: 0, objects: 0, chord: 0, relaxed: 0 };
  const faceSkips = {};
  const flatOf = (pts) => { const out = new Float64Array(3 * pts.length); pts.forEach((p, i) => { out[3 * i] = p[0]; out[3 * i + 1] = p[1]; out[3 * i + 2] = p[2]; }); return out; };
  const boxDiag = (flat) => { let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity]; for (let i = 0; i + 2 < flat.length; i += 3) for (let k = 0; k < 3; k++) { const v = flat[i + k]; if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; } return lo[0] <= hi[0] ? Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) : 0; };
  const deflOf = (diag) => (opts.deflection !== undefined ? opts.deflection / S : 0.01 * diag);
  let nCurve = 0;

  const emit = (g, M, layer, label, depth) => {
    switch (g.t) {
      case 'mesh': {
        if (!g.nf) { tally(skipped, 'ON_Mesh without faces'); return; }
        const { V, T } = meshTriangles(g, null);
        addTriangles(layer, label || 'Mesh', V, T, M);
        return;
      }
      case 'brep': {
        if (!g.nf) { tally(skipped, 'ON_Brep without faces'); return; }
        const group = label || 'Brep';
        const missing = [];
        const weld = new Welder(0, 1024), T = [];
        for (let i = 0; i < g.nf; i++) {
          const m = stored ? g.meshes[i] : null;
          if (!m || !m.nf) { missing.push(i); continue; }
          const r = meshTriangles(m, weld);
          for (let k = 0; k < r.T.length; k++) T.push(r.T[k]);
        }
        if (T.length) { addTriangles(layer, group, weld.finish(), T, M); tally(hist, 'B-rep faces from stored meshes', g.nf - missing.length); }
        if (!missing.length) return;
        // faces without a stored mesh: triangulate them from the surface and the edge curves
        const b = g.box, diag = b && b[0].every(Number.isFinite) && b[1].every(Number.isFinite) && b[1][0] >= b[0][0] ? Math.hypot(b[1][0] - b[0][0], b[1][1] - b[0][1], b[1][2] - b[0][2]) : 0;
        let d = diag;
        if (!(d > 0)) { const t = brepTopology(g, ar); if (!t.why) d = boxDiag(flatOf(t.V)); }
        if (!(d > 0)) { tally(skipped, 'ON_Brep faces without a stored mesh (B-rep of no extent)', missing.length); return; }
        let defl = deflOf(d), pass = tessellateBrep(g, ar, missing, defl, angStep, d, cells);
        // boundary loops that touch at this facet size come apart at a finer one: one more try
        if (pass.crossing) { const finer = tessellateBrep(g, ar, missing, defl / 4, angStep, d, cells); if (finer.read > pass.read) { pass = finer; defl /= 4; } }
        tess.faces += missing.length; tess.facesRead += pass.read; tess.objects++; tess.relaxed += pass.relaxed;
        tess.chord = Math.max(tess.chord, defl);
        for (const k in pass.skips) tally(faceSkips, k, pass.skips[k]);
        for (const k in pass.hist) tally(hist, `B-rep faces triangulated: ${k}`, pass.hist[k]);
        if (!pass.mesh.empty) { const fm = pass.mesh.finish(null); addTriangles(layer, group, fm.coords, fm.indices, M); tally(hist, 'tri3', fm.indices.length / 3); }
        return;
      }
      case 'surface': {
        tess.faces++; tess.objects++;
        const ev = surfaceEval(g);
        if (ev.why) { tally(faceSkips, ev.why); return; }
        // extent from a coarse sample of the patch
        const [[u0, u1], [v0, v1]] = ev.dom, sample = [];
        for (let j = 0; j <= 8; j++) for (let i = 0; i <= 8; i++) { const p = ev.at(u0 + ((u1 - u0) * i) / 8, v0 + ((v1 - v0) * j) / 8); sample.push(p[0], p[1], p[2]); }
        const diag = boxDiag(sample);
        if (!(diag > 0)) { tally(faceSkips, `${ev.how} of no extent`); return; }
        const defl = deflOf(diag), res = patchMesh(ev, defl, angStep, diag, 'surface', cells);
        if (typeof res === 'string') { tally(faceSkips, `${ev.how}: ${res}`); return; }
        tess.facesRead++; tess.chord = Math.max(tess.chord, defl);
        if (res.relaxed) tess.relaxed++;
        tally(hist, `surfaces triangulated: ${ev.how}`); tally(hist, 'tri3', res.mesh.indices.length / 3);
        addTriangles(layer, label || 'Surface', res.mesh.coords, res.mesh.indices, M);
        return;
      }
      case 'extrusion': {
        const group = label || 'Extrusion';
        if (stored && g.mesh && g.mesh.nf) { const { V, T } = meshTriangles(g.mesh, null); addTriangles(layer, group, V, T, M); tally(hist, 'extrusions from stored meshes'); return; }
        const res = extrusionMesh(g, curveStep);
        if (typeof res === 'string') { tally(skipped, `ON_Extrusion without a stored mesh (${res})`); return; }
        addTriangles(layer, group, res.V, res.T, M);
        tally(hist, 'extrusions triangulated'); tally(hist, 'tri3', res.T.length / 3);
        return;
      }
      case 'curve': {
        if (opts.curves === false) { tally(skipped, 'curves (left out on request)'); return; }
        const pts = curvePoints(g, undefined, undefined, curveStep);
        if (!pts || pts.length < 2) { tally(skipped, `${classOf(g)} that could not be evaluated`); return; }
        const p = new Float64Array(3 * pts.length);
        pts.forEach((q, i) => M4.apply(M, q[0], q[1], q[2], p, 3 * i));
        nCurve++;
        const a = pts[0], z = pts[pts.length - 1];
        const m = { kind: 'polyline', name: label || `${layer} curve ${nCurve}`, points: p, layer, entity: classOf(g) };
        if (pts.length > 2 && Math.hypot(a[0] - z[0], a[1] - z[1], a[2] - z[2]) <= 1e-9 * (Math.abs(a[0]) + Math.abs(a[1]) + Math.abs(a[2]) + 1)) m.closed = true;
        polylines.push(m);
        tally(hist, 'curve');
        return;
      }
      case 'point': {
        let s = pointSets.get(layer);
        if (!s) { s = new Grow(Float64Array, 96); pointSets.set(layer, s); }
        s.ensure(3); M4.apply(M, g.p[0], g.p[1], g.p[2], s.a, s.n); s.n += 3;
        tally(hist, 'point');
        return;
      }
      case 'cloud': {
        if (!g.P.length) { tally(skipped, 'ON_PointCloud without points'); return; }
        const p = new Float64Array(g.P.length);
        for (let i = 0; i + 2 < p.length; i += 3) M4.apply(M, g.P[i], g.P[i + 1], g.P[i + 2], p, i);
        clouds.push({ name: label || `${layer} point cloud ${clouds.length + 1}`, p, layer });
        tally(hist, 'point', p.length / 3);
        return;
      }
      case 'iref': {
        const d = ar.idefs.get(g.def);
        if (!d) { tally(skipped, 'ON_InstanceRef to a block definition that is not in the file'); return; }
        if (d.linked) { tally(skipped, 'ON_InstanceRef to a linked block (its geometry is in another file)'); return; }
        if (depth > 16) { tally(skipped, 'ON_InstanceRef nested more than 16 deep'); return; }
        const Mi = M4.mul(M, g.M), name = label || cleanName(d.name, 'Block');
        tally(hist, 'block instance');
        for (const id of d.members) {
          const o = byObjectId.get(id);
          if (!o) { tally(skipped, 'block member that is not in the object table'); continue; }
          if (!opts.hidden && (!o.attrs.visible || o.attrs.mode === 1)) continue;
          emit(o.geom, Mi, layerName(o.attrs.layer), o.attrs.name ? cleanName(o.attrs.name) : name, depth + 1);
        }
        return;
      }
      case 'bad': tally(skipped, `${g.cls} that could not be decoded (${g.why})`); return;
      default: tally(skipped, g.cls || g.t);
    }
  };

  const I = M4.identity();
  for (const o of top) {
    try { emit(o.geom, I, layerName(o.attrs.layer), o.attrs.name ? cleanName(o.attrs.name) : '', 0); } catch (e) {
      if (e instanceof Flaw) tally(skipped, `${classOf(o.geom)} that could not be decoded (${e.message})`); else throw e;
    }
  }

  // ---- models ----
  const name = baseName(fname);
  const scaleIn = (c) => { if (S !== 1) for (let i = 0; i < c.length; i++) c[i] *= S; return c; };
  const models = [];
  const single = surfaces.size === 1;
  for (const [layer, s] of surfaces) {
    if (!s.coords.n) continue;
    const { indices, groups } = groupsFromBuckets(s.buckets);
    if (!indices.length) continue;
    models.push(surfaceModel(single ? name : layer, scaleIn(s.coords.finish()), indices, groups.length > 1 ? groups : null, { layer }));
  }
  for (const l of polylines) { scaleIn(l.points); models.push(l); }
  for (const [layer, s] of pointSets) models.push(pointsModel(pointSets.size === 1 && !clouds.length ? `${name} points` : `${layer} points`, scaleIn(s.finish()), null, { layer }));
  for (const c of clouds) models.push(pointsModel(c.name, scaleIn(c.p), null, { layer: c.layer }));

  // ---- report ----
  if (hidden) { skipped['hidden objects'] = hidden; warnings.push(`${hidden} hidden object(s) were left out; import with hidden = true to keep them.`); }
  if (hiddenLayer) { skipped['objects on hidden layers'] = hiddenLayer; warnings.push(`${hiddenLayer} object(s) on hidden layers were left out; import with hidden = true to keep them.`); }
  const nFaceSkips = Object.values(faceSkips).reduce((a, b) => a + b, 0);
  if (tess.faces) {
    const list = Object.entries(faceSkips).sort((a, b) => b[1] - a[1]).map(([k, n]) => `${n} × ${k}`);
    warnings.push(`${tess.faces} face(s) of ${tess.objects} object(s) ${stored ? 'carry no render mesh (the file was saved with "Save Small", or written by a program that stores none) and ' : ''}were triangulated here from their surfaces at a chord error of ${(tess.chord * S).toPrecision(3)} m: ${tess.facesRead} read${tess.relaxed ? ` (${tess.relaxed} of them with larger facets than that, to keep the triangle count in bounds: their shape is coarser and their area can be off by a few per cent — import with deflection or angularStep set to have them triangulated as asked)` : ''}${nFaceSkips ? `, ${nFaceSkips} not (${list.slice(0, 6).join('; ')}${list.length > 6 ? '; …' : ''})` : ''}.${nFaceSkips ? ' For the faces left out, re-save the model in Rhino with render meshes (shade the view, then save without "Save Small") or export STEP.' : ''}`);
    if (nFaceSkips) skipped['B-rep faces without a stored mesh that could not be triangulated'] = nFaceSkips;
  }
  for (const k in skipped) {
    if (/^hidden objects$|^objects on hidden layers$|^B-rep faces without/.test(k)) continue;
    if (k === 'ON_SubD') warnings.push(`Skipped ${skipped[k]} SubD object(s): subdivision surfaces are not evaluated here. Convert them in Rhino (ToNURBS or Mesh) and save again.`);
    else if (/^ON_Extrusion without a stored mesh/.test(k)) warnings.push(`Skipped ${skipped[k]} ${k}; re-save the model in Rhino with render meshes or export STEP.`);
    else warnings.push(`Skipped ${skipped[k]} ${k}.`);
  }
  if (notes.crc) warnings.push(`${notes.crc} mesh vertex buffer(s) fail their CRC-32 check; the file may be damaged.`);
  if (ar.flaws.length) warnings.push(`${ar.flaws.length} record(s) could not be decoded and were passed over (first: ${ar.flaws[0]}).`);
  if (!ar.complete) warnings.push('The file has no end-of-file mark: it may be incomplete.');
  if (ar.newer) warnings.push(`Archive version ${ar.version} is newer than the versions this reader was checked with (2 – 80).`);
  if (!models.length) {
    const n = ar.objects.length;
    if (!n) throw new Error(`${fname}: the 3DM file holds no objects (tables found: ${ar.tables.join(', ') || 'none'}).`);
    const why = [...Object.entries(skipped).filter(([k]) => !/^B-rep faces without/.test(k)), ...Object.entries(faceSkips).map(([k, c]) => [`face(s) not triangulated — ${k}`, c])].map(([k, c]) => `${c} × ${k}`).slice(0, 8).join('; ');
    throw new Error(`${fname}: none of the ${n} object(s) of the 3DM file could be turned into geometry (${why || 'nothing readable'}).${hidden + hiddenLayer ? ' Import with hidden = true to keep hidden objects.' : ''}${tess.faces ? ' Re-save the model in Rhino with render meshes (not "Save Small") or export STEP.' : ''}`);
  }
  const meta = {
    format: 'rhino-3dm', archiveVersion: ar.version, openNurbsVersion: ar.openNurbs, application: ar.application || null, comment: ar.comment || null,
    units, upAxis: { source: 'z', basis: 'file' }, crs: null, elementTypes: hist, skipped, objects: classes,
    layers: ar.layers.map((l) => pathOf(l)), blockDefinitions: ar.idefs.size, hidden: { objects: hidden, onHiddenLayers: hiddenLayer, kept: !!opts.hidden },
  };
  if (tess.faces) meta.tessellation = { faces: tess.faces, facesRead: tess.facesRead, facesRelaxed: tess.relaxed, chord: tess.chord * S, angularStepDeg: opts.angularStep ?? 15 };
  const r = finish(models, warnings, meta);
  if (units.basis === 'assumed') unitAdvice(units, r.meta.bbox, suggestUnits, warnings);
  return r;
}

function classOf(g) {
  switch (g.t) {
    case 'mesh': return 'ON_Mesh';
    case 'brep': return 'ON_Brep';
    case 'extrusion': return 'ON_Extrusion';
    case 'curve': return { nurbs: 'ON_NurbsCurve', line: 'ON_LineCurve', arc: 'ON_ArcCurve', polyline: 'ON_PolylineCurve', poly: 'ON_PolyCurve' }[g.k];
    case 'surface': return { nurbs: 'ON_NurbsSurface', plane: 'ON_PlaneSurface', rev: 'ON_RevSurface', sum: 'ON_SumSurface' }[g.k];
    case 'point': return 'ON_Point';
    case 'cloud': return 'ON_PointCloud';
    case 'iref': return 'ON_InstanceRef';
    default: return g.cls || g.t;
  }
}

/**
 * Walls and caps of an extrusion: the profile (in the x–y plane) is carried along the path, with y along the
 * `up` vector and x along up × path. Mitred ends are not handled.
 * @returns {{ V: Float64Array, T: number[] }|string}
 */
function extrusionMesh(x, arcStep) {
  if (!x.profile || x.profile.t !== 'curve') return 'no profile curve';
  if ((x.haveN[0] && Math.abs(x.N[0][2]) < 1 - 1e-9) || (x.haveN[1] && Math.abs(x.N[1][2]) < 1 - 1e-9)) return 'mitred ends';
  const A = x.path[0], B = x.path[1], Tn = V3.norm(V3.sub(B, A)), up = V3.norm(x.up);
  if (!(V3.len(Tn) > 0) || !(V3.len(up) > 0)) return 'no path';
  const X = V3.norm(V3.cross(up, Tn));
  const at = (t) => [A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]), A[2] + t * (B[2] - A[2])];
  const P0 = at(x.pt[0]), P1 = at(x.pt[1]);
  // profiles: the curve itself, or the segments of a polycurve when several profiles are stored
  const profs = x.profiles > 1 && x.profile.k === 'poly' && x.profile.segs.length === x.profiles ? x.profile.segs : [x.profile];
  const rings = [];
  for (const c of profs) { const pts = curvePoints(c, undefined, undefined, arcStep); if (!pts || pts.length < 2) return 'its profile could not be evaluated'; rings.push(pts); }
  const V = [], T = [];
  const closedRings = [];
  for (const pts of rings) {
    const a = pts[0], z = pts[pts.length - 1], closed = pts.length > 3 && Math.hypot(a[0] - z[0], a[1] - z[1]) <= 1e-9 * (1 + Math.abs(a[0]) + Math.abs(a[1]));
    const ring = closed ? pts.slice(0, -1) : pts, n = ring.length, base = V.length / 3;
    for (const O of [P0, P1]) for (const p of ring) V.push(O[0] + p[0] * X[0] + p[1] * up[0], O[1] + p[0] * X[1] + p[1] * up[1], O[2] + p[0] * X[2] + p[1] * up[2]);
    for (let i = 0; i + 1 < n || (closed && i < n); i++) { const j = (i + 1) % n; T.push(base + i, base + j, base + n + j, base + i, base + n + j, base + n + i); }
    if (closed) closedRings.push({ ring, base, n });
  }
  const capped = x.capsKnown ? x.cap : [closedRings.length > 0, closedRings.length > 0];
  if ((capped[0] || capped[1]) && closedRings.length && closedRings.length === rings.length) {
    // outer profile first, the others are holes
    const flat = [], holes = [], bottom = [], top = [];
    closedRings.forEach((cr, k) => { if (k) holes.push(flat.length / 2); cr.ring.forEach((p, i) => { flat.push(p[0], p[1]); bottom.push(cr.base + i); top.push(cr.base + cr.n + i); }); });
    const tri = earclip(flat, holes, 2);
    for (let i = 0; i + 2 < tri.length; i += 3) {
      if (capped[0]) T.push(bottom[tri[i]], bottom[tri[i + 2]], bottom[tri[i + 1]]);
      if (capped[1]) T.push(top[tri[i]], top[tri[i + 1]], top[tri[i + 2]]);
    }
  }
  if (!T.length) return 'its profile is degenerate';
  return { V: Float64Array.from(V), T };
}
