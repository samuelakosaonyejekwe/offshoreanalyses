// ANSYS Fluent mesh and case files (.msh, .cas): the legacy parenthesised format, ASCII and binary.
//
// Format references: the "Mesh File Format" appendix of the Fluent User's Guide (section ids, zone and
// element type numbers), and the section layouts as described by the meshio `ansys` reader (MIT) and the
// OpenFOAM manual pages of fluentMeshToFoam / fluent3DMeshToFoam. No code was taken from any of them.
//
// A file is a sequence of sections "(id …)". Integers in the section headers and in the face lists are
// hexadecimal; zone ids in the zone sections are decimal.
//
//   (0 "text")                                   comment (also written without quotes, over several lines)
//   (1 "TGrid 3D 3.6.8")                         header: the writer
//   (2 3)                                        dimension
//   (4 (60 0 0 1 2 4 4 4 8 4 8 4))               machine configuration
//   (10 (zone first last type ND) (x y [z] …))   nodes; zone 0 only declares the total
//   (12 (zone first last type elem-type) (…))    cells; the body exists for mixed zones only
//   (13 (zone first last bc-type face-type) (n0 n1 … c0 c1 …))
//                                                faces: node numbers, then the two adjacent cells (0 = none).
//                                                face-type 2 line, 3 triangle, 4 quadrilateral; 0 mixed and
//                                                5 polygonal rows start with their node count
//   (39 (id type name [domain]) (…))  (45 (id type name) ())      zone type and name
//
// Binary sections carry the id + 2000 (single precision reals) or + 3000 (double precision). Their body
// follows the header list as "(" raw little-endian data ")" and the section ends with the text
// "End of Binary Section   3013)". Integers are 32-bit, or 64-bit in files of builds with 64-bit indices:
// the width is the one for which the end marker is found where the body should end.
//
// What is read: the nodes and the faces. A cell is not stored in the file — it is the set of faces that
// name it — so the interior is not reconstructed. The result is the boundary of the mesh: every face with
// exactly one adjacent cell (and every face with none, which is how a surface mesh is written), grouped by
// face zone with the zone name and type. Faces are wound so that their normal points out of the cell.
//   3-D mesh → one `surface` model "(boundary)", triangles; polygonal faces are ear-clipped
//   2-D mesh → one `frame` model "(boundary edges)" whose member groups are the zones (they are
//              boundaries, not beams)
// Faces between two cells are counted (`meta.skipped.internalFaces`); those that lie in a zone that is not
// of type interior — two-sided walls, fans, porous jumps, radiators — are counted separately
// (`twoSidedFaces`) and returned as well with `opts.twoSided: true`. Zones of type parent-face (the unrefined
// faces of a hanging-node mesh, whose children are in the file too) are left out.
//
// Not read, counted in `meta.sections` / `meta.skipped`: cells (12), edges (11), periodic shadow faces (18),
// cell and face trees (58, 59), interface face parents (61), partition and node flags (40, 41), the
// variable sections of case files (37, 38, 64 …) and any other section. Data files (.dat) are refused.
// Fluent's HDF5 files (.msh.h5, .cas.h5) hold the same tables and are read by adapters/hdf5.js, which
// hands them to `fluentBoundary` below.
// Coordinates are taken as metres: Fluent stores SI whatever unit is shown on screen. A mesh that was never
// scaled in Fluent is still in the unit of the mesh generator — `opts.units` says so.

import { inputBytes, decodeText, baseName, Grow, finish, surfaceModel, pushPolygon, tally, cleanName } from './common.js';
import { resolveUnits } from '../units.js';

/** Boundary-condition numbers of face-section headers (the zone sections give the same as text). */
const BC = { 2: 'interior', 3: 'wall', 4: 'pressure-inlet', 5: 'pressure-outlet', 7: 'symmetry', 8: 'periodic-shadow', 9: 'pressure-far-field', 10: 'velocity-inlet', 12: 'periodic', 14: 'fan', 20: 'mass-flow-inlet', 24: 'interface', 31: 'parent-face', 36: 'outflow', 37: 'axis' };
const CELL_ZONE = { 0: 'dead', 1: 'fluid', 17: 'solid', 32: 'inactive' };
const CELL_SHAPE = { 0: 'mixed', 1: 'triangle', 2: 'tetrahedron', 3: 'quadrilateral', 4: 'hexahedron', 5: 'pyramid', 6: 'wedge', 7: 'polyhedron' };
const SECTION = { 0: 'comment', 1: 'header', 2: 'dimension', 4: 'machine configuration', 10: 'nodes', 11: 'edges', 12: 'cells', 13: 'faces', 18: 'periodic shadow faces', 37: 'variables', 38: 'variables', 39: 'zone', 40: 'partition', 41: 'node flags', 45: 'zone', 58: 'cell tree', 59: 'face tree', 61: 'interface face parents' };
const READ = new Set([0, 1, 2, 4, 10, 13, 39, 45]);
const MARK = [69, 110, 100, 32, 111, 102, 32, 66, 105, 110, 97, 114, 121, 32, 83, 101, 99, 116, 105, 111, 110]; // "End of Binary Section"
const MAX_POLY = 4096;
// the node numbers and the two cells of the face being read (doubles: 64-bit files)
const ids = new Float64Array(MAX_POLY + 2);
const P10 = [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];
// byte → hexadecimal digit value, −1 otherwise
const HEX = new Int8Array(256).fill(-1);
for (let i = 0; i < 10; i++) HEX[48 + i] = i;
for (let i = 0; i < 6; i++) { HEX[65 + i] = 10 + i; HEX[97 + i] = 10 + i; }

const sectionName = (id) => { const base = id >= 2000 && id < 4000 ? id % 1000 : id; return SECTION[base] ? `${id} ${SECTION[base]}` : String(id); };

/**
 * Fluent mesh or case file → the boundary faces, one group per face zone.
 * @param {{name?: string, bytes?: Uint8Array, text?: string}} input
 * @param {object} [opts] `units` / `scale` (the file is in metres), `twoSided: true` to also return the faces
 *   of non-interior zones that lie between two cells
 */
export function parseFluent(input, opts = {}) {
  const b = inputBytes(input);
  const n = b.length;
  const fname = input.name || 'Fluent mesh';
  const dv = new DataView(b.buffer, b.byteOffset, n);
  let p = n >= 3 && b[0] === 0xef && b[1] === 0xbb && b[2] === 0xbf ? 3 : 0;
  const fail = (msg) => { throw new Error(`${fname}: ${msg}`); };

  if (b[p] === 0x89 && b[p + 1] === 0x48 && b[p + 2] === 0x44 && b[p + 3] === 0x46) fail('this is a Fluent file in the HDF5-based common fluids format (.msh.h5 / .cas.h5). It is read by the HDF5 reader, not by this one: import it as HDF5, or write the legacy format from Fluent (enter "/file/cff-files? no" in the console and write the mesh or case again).');
  if (b[p] === 0x1f && b[p + 1] === 0x8b) fail('the file is still gzip-compressed. Decompress it (gunzip) and import the .msh / .cas inside.');
  while (p < n && b[p] <= 32) p++;
  if (p >= n) fail('the file is empty (expected Fluent mesh sections such as (2 3), (10 …) and (13 …)).');
  if (b[p] === 36) fail('this .msh file starts with a "$" keyword: it is a Gmsh mesh, not a Fluent mesh. Import it as Gmsh MSH.');
  if (b[p] !== 40) fail('not a Fluent mesh or case file (expected the file to start with a parenthesised section such as (0 "comment"), (1 "header") or (2 3)).');

  // ---- low-level scanning -------------------------------------------------------------------------------
  const skipWs = () => { while (p < n && b[p] <= 32) p++; };
  /** Skip the list that opens at p, honouring nested lists and "strings"; p ends after its ")". */
  const skipList = (what, start) => {
    let depth = 0;
    while (p < n) {
      const c = b[p++];
      if (c === 40) depth++;
      else if (c === 41) { if (--depth === 0) return; } else if (c === 34) {
        while (p < n && b[p] !== 34) p += b[p] === 92 ? 2 : 1;
        p++;
      } else if (c === 35 && b[p] === 92) p += 2; // a Scheme character such as #\( in the variables of a case file
    }
    // The list never closes. When a mesh section starts on a later line the fault is in this list (a quote or
    // parenthesis the scanner misread, in text that is not needed) and reading goes on there.
    for (let q = start + 1; q + 4 < n; q++) {
      if (b[q] !== 40 || (b[q - 1] !== 10 && b[q - 1] !== 13)) continue;
      let r = q + 1, v = 0;
      while (r < n && r - q < 6 && b[r] >= 48 && b[r] <= 57) v = v * 10 + b[r++] - 48;
      const base = v >= 2000 && v < 4000 ? v % 1000 : v;
      if ((base === 10 || base === 12 || base === 13 || v === 39 || v === 45) && (b[r] === 32 || b[r] === 40)) { notes.push(`The parentheses of ${what} at byte ${start} do not balance; reading went on at the next mesh section (byte ${q}).`); p = q; return; }
    }
    p = n;
    fail(`the file ends inside ${what} that starts at byte ${start} (unbalanced parentheses: the file is truncated).`);
  };
  /** Hexadecimal (or decimal) integer token at p; −1 when there is none. */
  const int = (radix) => {
    let v = 0, any = false;
    while (p < n) {
      const d = HEX[b[p]];
      if (d < 0 || d >= radix) break;
      v = v * radix + d; any = true; p++;
    }
    return any ? v : -1;
  };
  /** The integers of a header list "(a b c …)" at p, hexadecimal. */
  const headerInts = (id) => {
    skipWs();
    if (b[p] !== 40) fail(`section ${id} at byte ${p}: expected a header list "(zone first last …)".`);
    p++;
    const out = [];
    for (;;) {
      skipWs();
      if (p >= n) fail(`the file ends inside the header of section ${id}.`);
      if (b[p] === 41) { p++; return out; }
      const v = int(16);
      if (v < 0 || out.length > 16 || (p < n && b[p] > 32 && b[p] !== 41)) fail(`section ${id} at byte ${p}: the header list holds something other than hexadecimal integers.`);
      out.push(v);
    }
  };
  /** Decimal real at p (Clinger's fast path; anything longer goes through Number()). */
  const real = () => {
    const s = p;
    let c = b[p], neg = false;
    if (c === 45) { neg = true; c = b[++p]; } else if (c === 43) c = b[++p];
    let m = 0, nd = 0, e = 0, any = false, exact = true;
    while (c >= 48 && c <= 57) { any = true; if (nd < 15) { m = m * 10 + (c - 48); if (m) nd++; } else { e++; if (c !== 48) exact = false; } c = b[++p]; }
    if (c === 46) {
      c = b[++p];
      while (c >= 48 && c <= 57) { any = true; if (nd < 15) { m = m * 10 + (c - 48); e--; if (m) nd++; } else if (c !== 48) exact = false; c = b[++p]; }
    }
    if (!any) return NaN;
    if (c === 101 || c === 69 || c === 100 || c === 68) {
      c = b[++p];
      let en = false, x = 0, ed = 0;
      if (c === 45) { en = true; c = b[++p]; } else if (c === 43) c = b[++p];
      while (c >= 48 && c <= 57) { if (ed < 6) { x = x * 10 + (c - 48); ed++; } c = b[++p]; }
      if (!ed) return NaN;
      e += en ? -x : x;
    }
    if (p > n) { p = n; return NaN; }
    let v;
    if (exact && e >= -22 && e <= 22) v = e < 0 ? m / P10[-e] : m * P10[e];
    else v = Number(String.fromCharCode.apply(null, b.subarray(s + (neg || b[s] === 43 ? 1 : 0), p)).replace(/[dD]/, 'e'));
    return neg ? -v : v;
  };
  /** Position after the end marker of a binary section whose body ends at q, or −1 when it is not there. */
  const markerAt = (q) => {
    // at most "\n)\n" before the text: bytes of the data that happen to be blanks must not be skipped
    let w = 0;
    while (w < 2 && q < n && (b[q] === 10 || b[q] === 13 || b[q] === 32)) { q++; w++; }
    if (b[q] === 41) q++;
    for (w = 0; w < 2 && q < n && (b[q] === 10 || b[q] === 13 || b[q] === 32); w++) q++;
    if (q + MARK.length > n) return -1;
    for (let k = 0; k < MARK.length; k++) if (b[q + k] !== MARK[k]) return -1;
    q += MARK.length;
    while (q < n && b[q] !== 41 && b[q] !== 10) q++;
    return q < n && b[q] === 41 ? q + 1 : -1;
  };
  /** Skip a binary section of unknown layout by searching its end marker from p. */
  const skipBinary = (id, start) => {
    for (let q = p; q + MARK.length <= n; q++) {
      if (b[q] !== 69 || b[q + 1] !== 110 || b[q + 2] !== 100) continue;
      const e = markerAt(q);
      if (e > 0) { p = e; return; }
    }
    p = n;
    fail(`the file ends inside binary section ${id} that starts at byte ${start} (no "End of Binary Section" marker follows: the file is truncated).`);
  };
  /** After the body of an ASCII section: empty lists some writers append, then the ")" of the section. */
  const closeSection = () => {
    for (;;) {
      skipWs();
      if (b[p] !== 40) break;
      let q = p + 1;
      while (q < n && b[q] <= 32) q++;
      if (b[q] !== 41) return;
      p = q + 1;
    }
    if (b[p] === 41) p++;
  };
  const i64 = (q) => { const hi = dv.getUint32(q + 4, true); return hi > 0x1fffff ? -1 : hi * 4294967296 + dv.getUint32(q, true); };

  // ---- what is collected --------------------------------------------------------------------------------
  let dim = 0, coordDim = 0, declNodes = 0, declFaces = 0, declCells = 0, nNodes = 0, nodesRead = 0;
  let coords = new Float64Array(0);
  const sections = {}, skippedSections = {}, comments = [], notes = [], zoneNames = new Map(), cellZones = [], faceZones = [];
  let writer = null, machine = null, sawBinary = false, sawAscii = false, intBytes = 0, single = false, double = false, cellsMax = 0, stray = 0;
  const wantTwoSided = !!opts.twoSided;
  const growNodes = (last) => {
    if (last <= nNodes) return;
    const c = new Float64Array(last * 3).fill(NaN);
    c.set(coords);
    coords = c; nNodes = last;
  };
  // ---- sections -----------------------------------------------------------------------------------------
  const nodeSection = (id, start) => {
    const h = headerInts(id);
    if (h.length < 3) fail(`node section at byte ${start}: the header holds ${h.length} value(s); expected (zone first last type [dimension]).`);
    const [zone, first, last] = h;
    const bin = id !== 10;
    skipWs();
    if (zone === 0 || (!bin && b[p] !== 40)) { if (zone === 0) declNodes = Math.max(declNodes, last); if (bin) skipBinary(id, start); else closeSection(); return; }
    const cnt = last - first + 1;
    if (!(first >= 1) || !(cnt >= 1)) fail(`node section of zone ${zone} at byte ${start}: the range ${first} – ${last} is not a valid node range.`);
    // a node takes four bytes at least ("0 0\n", or two single-precision reals)
    if (last > n / 4 + 1) fail(`node section of zone ${zone} at byte ${start} numbers its nodes up to ${last}, more than a file of ${n} bytes can hold: the header is damaged or the file is truncated.`);
    if (bin) {
      const nd = h[4] === 2 || h[4] === 3 ? h[4] : dim || coordDim || 3;
      const sz = id >= 3000 ? 8 : 4;
      if (b[p] !== 40) fail(`binary node section of zone ${zone} at byte ${start}: no data follows the header.`);
      p++;
      const end = p + cnt * nd * sz;
      if (end > n) fail(`binary node section of zone ${zone} declares ${cnt} nodes (${cnt * nd * sz} bytes) but the file ends after ${n - p} bytes of it: the file is truncated.`);
      const after = markerAt(end);
      if (after < 0) fail(`binary node section of zone ${zone} at byte ${start}: the end marker is not where ${cnt} nodes of ${nd} ${sz === 8 ? 'double' : 'single'}-precision coordinates end. The section uses a layout that is not read.`);
      growNodes(last);
      let o = (first - 1) * 3, q = p;
      for (let i = 0; i < cnt; i++, o += 3) {
        if (sz === 8) { coords[o] = dv.getFloat64(q, true); coords[o + 1] = dv.getFloat64(q + 8, true); coords[o + 2] = nd === 3 ? dv.getFloat64(q + 16, true) : 0; } else { coords[o] = dv.getFloat32(q, true); coords[o + 1] = dv.getFloat32(q + 4, true); coords[o + 2] = nd === 3 ? dv.getFloat32(q + 8, true) : 0; }
        q += nd * sz;
      }
      if (sz === 8) double = true; else single = true;
      coordDim = Math.max(coordDim, nd);
      nodesRead += cnt;
      p = after;
      return;
    }
    p++;
    sawAscii = true;
    // two characters at least per coordinate, two coordinates at least per node
    if (cnt > (n - p) / 4 + 1) fail(`node section of zone ${zone} declares ${cnt} nodes, more than the ${n - p} bytes left in the file can hold: the file is truncated or the header is damaged.`);
    const tmp = new Float64Array(cnt * 3);
    let k = 0;
    for (;;) {
      while (p < n && b[p] <= 32) p++;
      if (p >= n) fail(`the file ends inside the node section of zone ${zone} after ${k} of ${cnt * (dim || 2)} coordinates: the file is truncated.`);
      if (b[p] === 41) { p++; break; }
      const at = p, v = real();
      if (v !== v || (p < n && b[p] > 32 && b[p] !== 41)) fail(`node section of zone ${zone}: "${cleanName(String.fromCharCode.apply(null, b.subarray(at, Math.min(n, at + 24))), '?')}" at byte ${at} is not a number.`);
      if (k >= tmp.length) fail(`node section of zone ${zone} declares ${cnt} nodes but holds more than ${tmp.length} coordinates.`);
      tmp[k++] = v;
    }
    const nd = k === cnt * 3 ? 3 : k === cnt * 2 ? 2 : 0;
    if (!nd) fail(`node section of zone ${zone} declares ${cnt} nodes but holds ${k} coordinates (expected ${cnt * 2} or ${cnt * 3}).`);
    growNodes(last);
    for (let i = 0, o = (first - 1) * 3, s = 0; i < cnt; i++, o += 3) { coords[o] = tmp[s++]; coords[o + 1] = tmp[s++]; coords[o + 2] = nd === 3 ? tmp[s++] : 0; }
    coordDim = Math.max(coordDim, nd);
    nodesRead += cnt;
    closeSection();
  };

  const cellSection = (id, start) => {
    const h = headerInts(id);
    if (h.length >= 3) {
      const [zone, first, last, type, shape] = h;
      if (zone === 0) declCells = Math.max(declCells, last);
      else { cellZones.push({ id: zone, count: Math.max(0, last - first + 1), typeCode: type ?? 1, shapeCode: shape ?? 0 }); cellsMax = Math.max(cellsMax, last); }
    }
    if (id !== 12) { skipBinary(id, start); return; }
    skipWs();
    if (b[p] === 40) {
      // the cell types of a mixed zone; a "(" inside means the writer left the list open
      let q = p + 1;
      while (q < n && (b[q] <= 32 || HEX[b[q]] >= 0)) q++;
      if (q >= n) { p = n; fail(`the file ends inside the cell section that starts at byte ${start}: the file is truncated.`); }
      if (b[q] !== 41) { p = q; if (b[q] !== 40) fail(`cell section at byte ${start}: unexpected character "${cleanName(String.fromCharCode(b[q]), '?')}" at byte ${q} (expected cell type numbers).`); notes.push(`The cell section at byte ${start} is not closed; reading went on at the next section.`); return; }
      p = q + 1;
    }
    closeSection();
  };

  /** Walk the binary rows of a face section; returns the offset after the last row, or −1 when they do not fit. */
  const binaryFaces = (z, cnt, fixed, w, store) => {
    let q = p;
    for (let f = 0; f < cnt; f++) {
      let k = fixed;
      if (!k) {
        if (q + w > n) return -1;
        k = w === 4 ? dv.getUint32(q, true) : i64(q);
        q += w;
        if (k < 0 || k > MAX_POLY) return -1;
      }
      const e = q + (k + 2) * w;
      if (e > n) return -1;
      if (store) {
        if (w === 4) { for (let i = 0; i < k; i++) ids[i] = dv.getUint32(q + 4 * i, true); addFace(z, ids, k, dv.getUint32(e - 8, true), dv.getUint32(e - 4, true), wantTwoSided); } else { for (let i = 0; i < k; i++) ids[i] = i64(q + 8 * i); addFace(z, ids, k, i64(e - 16), i64(e - 8), wantTwoSided); }
      }
      q = e;
    }
    return q;
  };

  const faceSection = (id, start) => {
    const h = headerInts(id);
    if (h.length < 3) fail(`face section at byte ${start}: the header holds ${h.length} value(s); expected (zone first last bc-type face-type).`);
    const [zone, first, last, bc] = h;
    const bin = id !== 13;
    skipWs();
    if (zone === 0 || (!bin && b[p] !== 40)) { if (zone === 0) declFaces = Math.max(declFaces, last); if (bin) skipBinary(id, start); else closeSection(); return; }
    const cnt = last - first + 1;
    const ft = h.length > 4 ? h[4] : 0;
    if (!(first >= 1) || !(cnt >= 1)) fail(`face section of zone ${zone} at byte ${start}: the range ${first} – ${last} is not a valid face range.`);
    if (ft !== 0 && ft !== 5 && !(ft >= 2 && ft <= 4)) fail(`face section of zone ${zone} at byte ${start}: face type ${ft} is not one of 0 (mixed), 2 (line), 3 (triangle), 4 (quadrilateral), 5 (polygon).`);
    const fixed = ft >= 2 && ft <= 4 ? ft : 0;
    const z = faceZone(zone, bc);
    faceZones.push(z);
    if (bin) {
      if (b[p] !== 40) fail(`binary face section of zone ${zone} at byte ${start}: no data follows the header.`);
      p++;
      let w = 0, after = -1;
      // rows with a node count are tried as 64-bit first: 32-bit rows never pass as 64-bit ones (a node number
      // would sit in the high word of the count), whereas 64-bit rows read as 32-bit can by accident
      for (const t of fixed ? [4, 8] : [8, 4]) {
        const end = fixed ? p + cnt * (fixed + 2) * t : binaryFaces(z, cnt, 0, t, false);
        if (end < 0 || end > n) continue;
        after = markerAt(end);
        if (after > 0) { w = t; break; }
      }
      if (!w) {
        if (p + cnt * ((fixed || 1) + 2) * 4 > n) fail(`binary face section of zone ${zone} declares ${cnt} faces but the file ends after ${n - p} bytes of it: the file is truncated.`);
        fail(`binary face section of zone ${zone} at byte ${start}: the end marker is not where ${cnt} faces of 32-bit or of 64-bit integers end. The section is damaged or uses a layout that is not read.`);
      }
      binaryFaces(z, cnt, fixed, w, true);
      intBytes = Math.max(intBytes, w);
      p = after;
      return;
    }
    p++;
    sawAscii = true;
    // "1 2 0 1" — a face takes eight characters at least
    if (cnt > (n - p) / 8 + 1) fail(`face section of zone ${zone} declares ${cnt} faces, more than the ${n - p} bytes left in the file can hold: the file is truncated or the header is damaged.`);
    const bad = (q) => { p = q; fail(`face section of zone ${zone}: unexpected character "${cleanName(String.fromCharCode(b[q]), '?')}" at byte ${q} (expected hexadecimal node and cell numbers).`); };
    const i0 = fixed ? 0 : -1;
    let q = p, f = 0, short = false;
    rows: for (; f < cnt; f++) {
      let k = fixed;
      for (let i = i0; ; i++) {
        while (q < n && b[q] <= 32) q++;
        if (q >= n) { p = n; fail(`the file ends inside the face section of zone ${zone} after ${f} of ${cnt} faces: the file is truncated.`); }
        let v = 0, d = HEX[b[q]];
        if (d < 0) {
          if (b[q] !== 41) bad(q);
          if (i !== i0) { p = q; fail(`face section of zone ${zone}: the list closes in the middle of face ${f + 1} at byte ${q}.`); }
          short = true;
          break rows;
        }
        do { v = v * 16 + d; d = HEX[b[++q]]; } while (d >= 0);
        if (q < n && b[q] > 32 && b[q] !== 41) bad(q);
        if (i < 0) { k = v; if (k > MAX_POLY) { p = q; fail(`face section of zone ${zone}: face ${f + 1} declares ${k} nodes (a mixed or polygonal face row starts with its node count, at most ${MAX_POLY}).`); } } else if (i <= k) ids[i] = v; else { addFace(z, ids, k, ids[k], v, wantTwoSided); break; }
      }
    }
    p = q;
    skipWs();
    if (short) notes.push(`The face section of zone ${zone} declares ${cnt} faces but holds ${f}.`);
    else if (p < n && b[p] !== 41) fail(`face section of zone ${zone} declares ${cnt} faces but more data follows them at byte ${p}${fixed ? '' : ' (or a row does not start with its node count)'}.`);
    if (b[p] === 41) p++;
    closeSection();
  };

  const zoneSection = (id, start) => {
    skipWs();
    if (b[p] === 40) {
      let q = p + 1;
      while (q < n && b[q] !== 41 && b[q] !== 40 && q - p < 600) q++;
      if (b[q] === 41) {
        const t = String.fromCharCode.apply(null, b.subarray(p + 1, q)).trim().split(/\s+/);
        const zid = Number(t[0]);
        if (t.length >= 3 && Number.isInteger(zid) && zid > 0) zoneNames.set(zid, { type: t[1], name: cleanName(t[2], `zone ${zid}`) });
      }
    }
    p = start;
    skipList(`zone section ${id}`, start);
  };

  for (;;) {
    skipWs();
    if (p >= n) break;
    if (b[p] !== 40) {
      // text between sections: not part of the format; skipped up to the next "(" that starts a line
      stray++;
      while (p < n && !(b[p] === 40 && (b[p - 1] === 10 || b[p - 1] === 13))) p++;
      continue;
    }
    const start = p++;
    skipWs();
    const id = int(10);
    if (id < 0 || (p < n && b[p] > 32 && b[p] !== 40 && b[p] !== 41 && b[p] !== 34)) {
      if (!Object.keys(sections).length) fail('not a Fluent mesh or case file (the first list does not start with a section number such as (0 "comment") or (2 3)).');
      p = start; skipList('an unnumbered list', start); tally(skippedSections, 'unnumbered list');
      continue;
    }
    tally(sections, sectionName(id));
    const binary = id >= 2000 && id < 4000;
    const base = binary ? id % 1000 : id;
    if (binary) sawBinary = true;
    if (base === 10) nodeSection(id, start);
    else if (base === 13) faceSection(id, start);
    else if (base === 12) cellSection(id, start);
    else if (binary) { skipBinary(id, start); tally(skippedSections, sectionName(id)); } else if (id === 39 || id === 45) zoneSection(id, start);
    else if (id === 2) { skipWs(); const d = int(10); if (d === 2 || d === 3) dim = d; else fail(`dimension section (2 …) at byte ${start}: expected 2 or 3.`); p = start; skipList('the dimension section', start); } else {
      const s = p;
      p = start; skipList(`section ${id}`, start);
      if (id === 0 || id === 1) {
        const t = cleanName(decodeText(b.subarray(s, Math.min(p - 1, s + 200))).replace(/^\s*["']|["']\s*$/g, ''), '');
        if (id === 1) writer = t; else if (t && comments.length < 12) comments.push(t.length > 120 ? `${t.slice(0, 120)}…` : t);
      } else if (id === 4) machine = String.fromCharCode.apply(null, b.subarray(s, Math.min(p - 1, s + 200))).replace(/[()]/g, ' ').trim().split(/\s+/).map(Number).filter(Number.isFinite);
      if (!READ.has(id)) tally(skippedSections, sectionName(id));
    }
  }

  if (!nodesRead) {
    const seen = Object.keys(sections);
    fail(seen.some((k) => /^(300|2300|3300)$/.test(k)) && !faceZones.length ? 'no mesh was found: this looks like a Fluent data file (.dat), which holds the solution and not the mesh. Import the case (.cas) or mesh (.msh) file instead.' : `no node section (10) was found${seen.length ? ` (sections present: ${seen.slice(0, 12).join(', ')})` : ''}: not a Fluent mesh or case file, or one without a mesh.`);
  }
  if (!faceZones.length) fail(`the file holds ${nodesRead} nodes but no face section (13), so there is no mesh to show.`);
  if (stray) notes.push(`${stray} run(s) of text outside any section were ignored.`);
  return fluentBoundary({ coords, nNodes, nodesRead, faceZones, zoneNames, cellZones, dim: dim || coordDim || 3, declNodes, declFaces, nCells: declCells || cellsMax || null, notes, skippedSections }, input, opts, {
    encoding: sawBinary ? (sawAscii ? 'mixed' : 'binary') : 'ascii', precision: sawBinary ? (double ? 'double' : single ? 'single' : null) : null, integerBytes: intBytes || null, writer, comments, machine, sections,
  });
}

/** A face zone of the collector: counters and the kept faces as [node count, node numbers …] runs. */
export function faceZone(id, bc) {
  return { id, bc, parent: bc === 31, faces: 0, boundary: 0, internal: 0, free: 0, empty: 0, kept: new Grow(Uint32Array, 1024) };
}

/**
 * Hand one face to its zone; it is kept when it bounds the mesh.
 * @param {ArrayLike<number>} ids the node numbers (1-based) of the face, `k` of them
 * @param {number} c0 @param {number} c1 the adjacent cells, 0 for none
 */
export function addFace(z, ids, k, c0, c1, twoSided) {
  z.faces++;
  if (z.parent) return;
  // faces without nodes: the virtual faces of a non-conformal periodic or mesh interface
  if (!k) { z.empty++; return; }
  let flip = false;
  if (c0 && c1) { z.internal++; if (!(twoSided && z.bc !== 2)) return; } else if (c0 || c1) { z.boundary++; flip = !c0; } else z.free++;
  const g = z.kept;
  g.ensure(k + 1);
  const a = g.a;
  let o = g.n;
  a[o++] = k;
  // the right-hand normal of the stored order points from c1 into c0; reversed so that it leaves the mesh
  if (flip) for (let i = 0; i < k; i++) a[o++] = ids[i] > 4294967295 ? 0 : ids[i]; else for (let i = k - 1; i >= 0; i--) a[o++] = ids[i] > 4294967295 ? 0 : ids[i];
  g.n = o;
}

/**
 * Models and report from the collected nodes and face zones. Also used by the HDF5 adapter, which reads the
 * same tables from Fluent's common fluids format.
 * @param {object} m `{ coords, nNodes, nodesRead, faceZones, zoneNames: Map(id → { name, type? }),
 *   cellZones: [{ id, count, typeCode, shapeCode }], dim, declNodes, declFaces, nCells, notes, skippedSections }`;
 *   `coords` holds x, y, z per node number − 1, NaN where the file gave none
 * @param {object} extra format-specific report fields
 */
export function fluentBoundary(m, input, opts = {}, extra = {}) {
  const { coords, nNodes, nodesRead, faceZones, zoneNames, cellZones, dim, declNodes, declFaces, nCells, skippedSections } = m;
  const fail = (msg) => { throw new Error(`${input.name || 'Fluent mesh'}: ${msg}`); };
  const wantTwoSided = !!opts.twoSided;
  const units = resolveUnits(opts, 'm', 'm');
  const scale = units.scale;
  const warnings = m.notes || [], hist = {}, patches = [], zones = [];
  const map = new Int32Array(nNodes).fill(-1);
  const used = new Grow(Float64Array, 3072);
  const vertex = (k) => { let v = map[k]; if (v < 0) { v = used.n / 3; map[k] = v; used.push3(coords[3 * k] * scale, coords[3 * k + 1] * scale, coords[3 * k + 2] * scale); } return v; };
  const tris = new Grow(Uint32Array, 3072), groups = [], members = [];
  const tri = [0, 0, 0], quad = [0, 0, 0, 0];
  let invalid = 0, internal = 0, twoSided = 0, parentFaces = 0, total = 0, free = 0, bnd = 0, empty = 0;
  const twoSidedZones = [];
  for (const z of faceZones) {
    const zn = zoneNames.get(z.id);
    const name = zn ? zn.name : `zone ${z.id}`;
    const type = (zn && zn.type) || BC[z.bc] || `type ${z.bc}`;
    total += z.faces; free += z.free; bnd += z.boundary; empty += z.empty;
    if (z.parent) parentFaces += z.faces;
    else if (z.bc === 2) internal += z.internal;
    else if (z.internal) { twoSided += z.internal; twoSidedZones.push(`${name} (${type}, ${z.internal})`); }
    const a = z.kept.a, end = z.kept.n, t0 = tris.n, m0 = members.length;
    let bad = 0;
    for (let o = 0; o < end;) {
      const k = a[o++];
      let ok = true;
      for (let i = 0; i < k; i++) { const v = a[o + i]; if (v < 1 || v > nNodes || coords[3 * v - 3] !== coords[3 * v - 3]) { ok = false; break; } }
      if (!ok || k < 2 || (k === 2 && a[o] === a[o + 1])) { bad++; o += k; continue; }
      if (k === 2) { members.push({ a: vertex(a[o] - 1), b: vertex(a[o + 1] - 1), group: name }); tally(hist, 'line2'); } else if (k === 3) { tri[0] = a[o] - 1; tri[1] = a[o + 1] - 1; tri[2] = a[o + 2] - 1; pushPolygon(tris, tri, coords); tally(hist, 'tri3'); } else if (k === 4) { for (let i = 0; i < 4; i++) quad[i] = a[o + i] - 1; pushPolygon(tris, quad, coords); tally(hist, 'quad4'); } else {
        const poly = new Array(k);
        for (let i = 0; i < k; i++) poly[i] = a[o + i] - 1;
        pushPolygon(tris, poly, coords); tally(hist, 'polygon');
      }
      o += k;
    }
    invalid += bad;
    const rec = { id: z.id, name, type, nFaces: z.faces, boundary: z.boundary + z.free, internal: z.internal };
    zones.push(rec);
    if (tris.n > t0) { groups.push({ name, start: t0, count: tris.n - t0, type }); patches.push({ name, type, nFaces: z.boundary + z.free + (wantTwoSided && z.bc !== 2 ? z.internal : 0) - bad, id: z.id }); } else if (members.length > m0) patches.push({ name, type, nFaces: members.length - m0, id: z.id });
    z.kept = null;
  }
  const cells = cellZones.map((c) => { const zn = zoneNames.get(c.id); return { id: c.id, count: c.count, shape: CELL_SHAPE[c.shapeCode] || `type ${c.shapeCode}`, name: zn ? zn.name : `zone ${c.id}`, type: (zn && zn.type) || CELL_ZONE[c.typeCode] || `type ${c.typeCode}` }; });
  const indices = tris.finish();
  for (let i = 0; i < indices.length; i++) indices[i] = vertex(indices[i]);
  if (!indices.length && !members.length) fail(`no boundary faces could be built: of the ${total} face(s) in ${faceZones.length} zone(s), ${internal + twoSided} lie between two cells${parentFaces ? `, ${parentFaces} are parent faces` : ''}${invalid ? ` and ${invalid} name nodes that are not in the file` : ''}.`);

  const base = baseName(input.name);
  const models = [];
  if (indices.length) models.push(surfaceModel(`${base} (${free && !bnd ? 'surface mesh' : 'boundary'})`, used.view(), indices, groups));
  if (members.length) models.push({ kind: 'frame', name: `${base} (boundary edges)`, nodes: used.finish(), members });
  if (indices.length && members.length) warnings.push(`The file mixes ${members.length} two-node faces with faces of three or more nodes; the former are returned as a separate model of edges.`);
  if (declNodes && nodesRead < declNodes) warnings.push(`The file declares ${declNodes} nodes but its node sections hold ${nodesRead}.`);
  if (declFaces && total !== declFaces) warnings.push(`The file declares ${declFaces} faces but its face sections hold ${total}.`);
  if (invalid) warnings.push(`${invalid} boundary face(s) that name nodes outside the node sections, or the same node twice, were skipped.`);
  if (twoSided) warnings.push(wantTwoSided ? `${twoSided} face(s) between two cells in zones that are not interior were kept as two-sided surfaces: ${twoSidedZones.join(', ')}.` : `${twoSided} face(s) lie between two cells in zones that are not of type interior (two-sided walls, fans, porous jumps) and were left out: ${twoSidedZones.join(', ')}. Import with twoSided = true to keep them.`);
  if (empty) warnings.push(`${empty} face(s) without nodes (virtual faces of periodic or non-conformal interfaces) were left out.`);
  if (parentFaces) warnings.push(`${parentFaces} parent face(s) of a hanging-node mesh were left out (their child faces are in the file).`);
  const sk = Object.keys(skippedSections || {});
  if (sk.length) warnings.push(`Sections not read: ${sk.map((k) => `${k}${skippedSections[k] > 1 ? ` ×${skippedSections[k]}` : ''}`).join(', ')}.`);
  if (free && !bnd) warnings.push('No face names a cell: the file is a surface mesh, returned whole.');
  else warnings.push(`Fluent stores cells only as the faces that bound them, and they are not reconstructed: only the boundary ${dim === 2 && !indices.length ? 'edges' : 'faces'} are returned${nCells ? ` (the mesh has ${nCells} cells in ${cellZones.length || 1} zone(s) and ${internal} interior faces)` : ''}.`);
  const skipped = { invalidFaces: invalid, internalFaces: internal, twoSidedFaces: wantTwoSided ? 0 : twoSided, parentFaces, facesWithoutNodes: empty, cells: nCells ?? 0 };
  for (const k of sk) skipped[`section ${k}`] = skippedSections[k];
  return finish(models, warnings, {
    format: 'fluent', units, crs: null, elementTypes: hist, skipped,
    points: nodesRead, faces: total, cells: nCells, patches, zones, cellZones: cells,
    dimension: dim, ...extra,
  });
}
