// VTK files: the legacy format (.vtk, ASCII and big-endian binary; versions 2–5.1) and the XML formats for
// unstructured grids and poly data (.vtu / .vtp; ASCII, inline base64 and appended base64/raw arrays, with or
// without zlib compression).
//
// Legacy datasets read: UNSTRUCTURED_GRID, POLYDATA, STRUCTURED_GRID, RECTILINEAR_GRID and STRUCTURED_POINTS
// (a single-layer STRUCTURED_POINTS image with a point scalar becomes a `grid`, the scalar being z).
// Point-data arrays are kept on the models as `nodeData[name]` (Float32Array, `nodeDataComponents[name]`
// values per node). Cell data, polyhedral cells (type 42) and the XML image/structured/parallel variants
// (.vti, .vts, .vtr, .pvtu) are not read.

import { inputBytes, baseName, FeMesh, LineCursor, decodeText, base64Decode, triangulate3D, finish } from './common.js';
import { feResult } from './fe.js';
import { resolveUnits } from '../units.js';
import { parseXml, kids, kid, attr, desc, numberList } from '../xml.js';
import { inflateZlib } from '../zip.js';

// VTK cell type → [canonical type, node count, canonical-from-native order]
const VTK_CELLS = {
  3: ['line2', 2], 21: ['line3', 3], 5: ['tri3', 3], 22: ['tri6', 6], 9: ['quad4', 4], 8: ['quad4', 4, [0, 1, 3, 2]], 23: ['quad8', 8], 28: ['quad9', 9],
  10: ['tet4', 4], 24: ['tet10', 10], 12: ['hex8', 8], 11: ['hex8', 8, [0, 1, 3, 2, 4, 5, 7, 6]], 25: ['hex20', 20], 29: ['hex27', 27],
  // VTK numbers the wedge base so that its normal points away from the top face; the canonical wedge has it pointing in
  13: ['wedge6', 6, [0, 2, 1, 3, 5, 4]], 26: ['wedge15', 15, [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13]], 14: ['pyr5', 5], 27: ['pyr13', 13],
};

function addCell(fe, type, ids, id, xyz) {
  const def = VTK_CELLS[type];
  if (def) {
    if (ids.length < def[1]) { fe.skip(`vtk cell type ${type} with too few points`); return; }
    fe.addElement(def[0], def[2] ? def[2].map((k) => ids[k]) : ids, id, 0);
  } else if (type === 4) { for (let i = 0; i + 1 < ids.length; i++) fe.addElement('line2', [ids[i], ids[i + 1]], 0, 0); } else if (type === 6) { for (let i = 0; i + 2 < ids.length; i++) fe.addElement('tri3', i & 1 ? [ids[i + 1], ids[i], ids[i + 2]] : [ids[i], ids[i + 1], ids[i + 2]], 0, 0); } else if (type === 7) {
    if (ids.length === 3) fe.addElement('tri3', ids, id, 0);
    else if (ids.length === 4) fe.addElement('quad4', ids, id, 0);
    else {
      const flat = new Float64Array(ids.length * 3);
      ids.forEach((k, i) => { flat[3 * i] = xyz[3 * k]; flat[3 * i + 1] = xyz[3 * k + 1]; flat[3 * i + 2] = xyz[3 * k + 2]; });
      const t = triangulate3D(flat);
      for (let i = 0; i + 2 < t.length; i += 3) fe.addElement('tri3', [ids[t[i]], ids[t[i + 1]], ids[t[i + 2]]], 0, 0);
    }
  } else if (type === 1 || type === 2) fe.skip('vertex cell', 1);
  else fe.skip(`vtk cell type ${type}`);
}

function attachNodeData(r, pointData, comps) {
  const names = Object.keys(pointData);
  if (!names.length) return;
  for (const m of r.models) {
    if (!m.nodeIds) continue;
    m.nodeData = {}; m.nodeDataComponents = {};
    for (const nm of names) {
      const nc = comps[nm], src = pointData[nm];
      const out = new Float32Array(m.nodeIds.length * nc);
      for (let i = 0; i < m.nodeIds.length; i++) for (let c = 0; c < nc; c++) out[i * nc + c] = src[m.nodeIds[i] * nc + c];
      m.nodeData[nm] = out; m.nodeDataComponents[nm] = nc;
    }
  }
  r.meta.pointData = names;
}

// ---------------------------------------------------------------------------------------------------------
// Legacy
// ---------------------------------------------------------------------------------------------------------

const LEGACY_TYPES = {
  bit: null, unsigned_char: [1, 'getUint8'], char: [1, 'getInt8'], unsigned_short: [2, 'getUint16'], short: [2, 'getInt16'],
  unsigned_int: [4, 'getUint32'], int: [4, 'getInt32'], unsigned_long: [8, 'getBigUint64'], long: [8, 'getBigInt64'],
  float: [4, 'getFloat32'], double: [8, 'getFloat64'], vtktypeint64: [8, 'getBigInt64'], vtktypeint32: [4, 'getInt32'], vtkidtype: [4, 'getInt32'],
  vtktypeuint8: [1, 'getUint8'], vtktypeuint64: [8, 'getBigUint64'],
};

/** Legacy VTK. */
export function parseVtk(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  const cur = new LineCursor(bytes);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const l0 = cur.next();
  if (l0 === null || !/^# vtk DataFile/i.test(l0)) throw new Error(`${input.name || 'VTK'}: not a legacy VTK file (expected the first line "# vtk DataFile Version x.y").`);
  const title = (cur.next() || '').trim();
  const enc = (cur.next() || '').trim().toUpperCase();
  if (enc !== 'ASCII' && enc !== 'BINARY') throw new Error(`${input.name || 'VTK'}: expected ASCII or BINARY on line 3, found "${enc.slice(0, 40)}".`);
  const binary = enc === 'BINARY';
  let queue = [], qi = 0;
  const token = () => {
    while (qi >= queue.length) { const l = cur.nextNonBlank(); if (l === null) return null; queue = l.split(/\s+/); qi = 0; }
    return queue[qi++];
  };
  const read = (n, type) => {
    const out = new Float64Array(n);
    const t = LEGACY_TYPES[String(type).toLowerCase()];
    if (t === undefined || t === null) throw new Error(`${input.name || 'VTK'}: data type "${type}" is not supported.`);
    if (!binary) { for (let i = 0; i < n; i++) { const s = token(); if (s === null) throw new Error(`${input.name || 'VTK'}: file ends while reading ${n} values (got ${i}).`); out[i] = Number(s); } return out; }
    let o = cur.offset;
    if (o + n * t[0] > bytes.length) throw new Error(`${input.name || 'VTK'}: binary block of ${n} × ${type} runs past the end of the file.`);
    const big = t[1].includes('Big');
    for (let i = 0; i < n; i++, o += t[0]) out[i] = big ? Number(dv[t[1]](o, false)) : dv[t[1]](o, false);
    cur.offset = o;
    queue = []; qi = 0;
    return out;
  };
  const keyLine = () => { queue = []; qi = 0; return cur.nextNonBlank(); };
  /** Cell list in either the classic "n id id …" layout or the 5.1 OFFSETS / CONNECTIVITY layout. */
  const cellList = (n, size) => {
    const save = { off: cur.offset, ln: cur.lineNo };
    const nl = cur.nextNonBlank();
    if (nl !== null && /^OFFSETS\b/i.test(nl)) {
      const offs = read(n, nl.split(/\s+/)[1]);
      const cl = keyLine();
      const conn = read(size, cl.split(/\s+/)[1]);
      const cells = [];
      for (let i = 0; i + 1 < n; i++) cells.push(Array.from(conn.subarray(offs[i], offs[i + 1])));
      return cells;
    }
    cur.offset = save.off; cur.lineNo = save.ln;
    const raw = read(size, 'int');
    const cells = [];
    for (let p = 0; p < size;) { const k = raw[p++]; cells.push(Array.from(raw.subarray(p, p + k))); p += k; }
    return cells;
  };

  const fe = new FeMesh(name);
  const warnings = [];
  const pointData = {}, comps = {};
  let dataset = null, nPoints = 0, mode = null, dims = null, origin = [0, 0, 0], spacing = [1, 1, 1], axes = {}, cells = null, types = null, skippedCellData = 0;
  const poly = {};
  const keep = (nm, nc, vals) => { if (mode === 'POINT' && vals.length === nPoints * nc) { pointData[nm] = vals; comps[nm] = nc; } else if (mode === 'CELL') skippedCellData++; };

  for (let line = keyLine(); line !== null; line = keyLine()) {
    const p = line.split(/\s+/);
    const kw = p[0].toUpperCase();
    if (kw === 'DATASET') dataset = (p[1] || '').toUpperCase();
    else if (kw === 'POINTS') { nPoints = Number(p[1]); const c = read(nPoints * 3, p[2]); for (let i = 0; i < nPoints; i++) fe.addNode(i, c[3 * i], c[3 * i + 1], c[3 * i + 2]); } else if (kw === 'DIMENSIONS') { dims = [Number(p[1]), Number(p[2]), Number(p[3] || 1)]; if (!nPoints) nPoints = dims[0] * dims[1] * dims[2]; }
    else if (kw === 'ORIGIN') origin = [Number(p[1]), Number(p[2]), Number(p[3] || 0)];
    else if (kw === 'SPACING' || kw === 'ASPECT_RATIO') spacing = [Number(p[1]), Number(p[2]), Number(p[3] || 1)];
    else if (/^[XYZ]_COORDINATES$/.test(kw)) axes[kw[0]] = read(Number(p[1]), p[2]);
    else if (kw === 'CELLS') cells = cellList(Number(p[1]), Number(p[2]));
    else if (kw === 'CELL_TYPES') types = read(Number(p[1]), 'int');
    else if (['VERTICES', 'LINES', 'POLYGONS', 'TRIANGLE_STRIPS'].includes(kw)) poly[kw] = cellList(Number(p[1]), Number(p[2]));
    else if (kw === 'POINT_DATA') mode = 'POINT';
    else if (kw === 'CELL_DATA') mode = 'CELL';
    else if (kw === 'SCALARS') {
      const nc = Number(p[3]) || 1;
      const n = mode === 'CELL' ? null : nPoints;
      const save = { off: cur.offset, ln: cur.lineNo };
      const lt = cur.nextNonBlank();
      if (lt === null || !/^LOOKUP_TABLE\b/i.test(lt)) { cur.offset = save.off; cur.lineNo = save.ln; }
      const count = (n ?? cellCount()) * nc;
      keep(p[1], nc, read(count, p[2]));
    } else if (kw === 'COLOR_SCALARS') { const nc = Number(p[2]); const vals = read((mode === 'CELL' ? cellCount() : nPoints) * nc, binary ? 'unsigned_char' : 'float'); keep(p[1], nc, vals); } else if (kw === 'VECTORS' || kw === 'NORMALS') keep(p[1], 3, read((mode === 'CELL' ? cellCount() : nPoints) * 3, p[2]));
    else if (kw === 'TENSORS') read((mode === 'CELL' ? cellCount() : nPoints) * 9, p[2]);
    else if (kw === 'TEXTURE_COORDINATES') read((mode === 'CELL' ? cellCount() : nPoints) * Number(p[2]), p[3]);
    else if (kw === 'LOOKUP_TABLE') read(Number(p[2]) * 4, binary ? 'unsigned_char' : 'float');
    else if (kw === 'FIELD') {
      const k = Number(p[2]);
      for (let i = 0; i < k; i++) {
        let h = keyLine();
        while (h !== null && /^METADATA\b/i.test(h)) { for (let x = cur.next(); x !== null && x.trim(); x = cur.next()); h = keyLine(); }
        if (h === null) break;
        const q = h.split(/\s+/);
        if (String(q[3]).toLowerCase() === 'string') { for (let j = 0; j < Number(q[1]) * Number(q[2]); j++) cur.next(); continue; }
        const vals = read(Number(q[1]) * Number(q[2]), q[3]);
        if (mode) keep(q[0], Number(q[1]), vals);
      }
    } else if (kw === 'METADATA') { for (let x = cur.next(); x !== null && x.trim(); x = cur.next()); } else if (kw === 'CELL_TYPES' || kw === 'OFFSETS' || kw === 'CONNECTIVITY') { /* handled with CELLS */ } else fe.skip(`section ${kw}`);
  }
  function cellCount() {
    if (cells) return cells.length;
    let n = 0;
    for (const k in poly) n += poly[k].length;
    if (!n && dims) n = Math.max(1, dims[0] - 1) * Math.max(1, dims[1] - 1) * Math.max(1, dims[2] - 1);
    return n;
  }
  if (!dataset) throw new Error(`${input.name || 'VTK'}: no DATASET line found (expected DATASET UNSTRUCTURED_GRID, POLYDATA, STRUCTURED_GRID, RECTILINEAR_GRID or STRUCTURED_POINTS).`);

  const structured = dataset === 'STRUCTURED_POINTS' || dataset === 'STRUCTURED_GRID' || dataset === 'RECTILINEAR_GRID';
  if (structured) {
    if (!dims) throw new Error(`${input.name || 'VTK'}: ${dataset} needs a DIMENSIONS line.`);
    const [nx, ny, nz] = dims;
    if (dataset === 'STRUCTURED_POINTS') {
      const zname = nz === 1 ? Object.keys(pointData).find((k) => comps[k] === 1) : null;
      if (zname && spacing[0] > 0 && spacing[1] > 0) {
        // a single-layer image with a scalar is a raster surface: the scalar is the elevation
        const units = resolveUnits(opts, null, 'm');
        const z = new Float32Array(nx * ny);
        const v = pointData[zname];
        for (let i = 0; i < z.length; i++) z[i] = v[i] * units.scale;
        const g = { kind: 'grid', name, nx, ny, x0: origin[0] * units.scale, y0: origin[1] * units.scale, dx: spacing[0] * units.scale, dy: spacing[1] * units.scale, z };
        return finish([g], [`STRUCTURED_POINTS image read as a grid with scalar "${zname}" as elevation.`], { format: 'vtk', encoding: enc.toLowerCase(), dataset, title, units, crs: null, elementTypes: { gridCell: nx * ny }, skipped: {}, pointData: Object.keys(pointData) });
      }
      nPoints = nx * ny * nz;
      for (let k = 0, i = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let ii = 0; ii < nx; ii++, i++) fe.addNode(i, origin[0] + ii * spacing[0], origin[1] + j * spacing[1], origin[2] + k * spacing[2]);
    } else if (dataset === 'RECTILINEAR_GRID') {
      const X = axes.X || [0], Y = axes.Y || [0], Z = axes.Z || [0];
      for (let k = 0, i = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let ii = 0; ii < nx; ii++, i++) fe.addNode(i, X[ii] || 0, Y[j] || 0, Z[k] || 0);
    }
    const at = (i, j, k) => i + nx * (j + ny * k);
    let id = 0;
    if (nz > 1 && ny > 1 && nx > 1) { for (let k = 0; k + 1 < nz; k++) for (let j = 0; j + 1 < ny; j++) for (let i = 0; i + 1 < nx; i++) fe.addElement('hex8', [at(i, j, k), at(i + 1, j, k), at(i + 1, j + 1, k), at(i, j + 1, k), at(i, j, k + 1), at(i + 1, j, k + 1), at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1)], ++id, 0); } else {
      // one layer thin in one direction: quads in the remaining two
      const A = nx > 1 ? nx : ny, B = nx > 1 ? (ny > 1 ? ny : nz) : nz;
      const idx = (a, b) => (nx > 1 ? (ny > 1 ? at(a, b, 0) : at(a, 0, b)) : at(0, a, b));
      for (let b = 0; b + 1 < B; b++) for (let a = 0; a + 1 < A; a++) fe.addElement('quad4', [idx(a, b), idx(a + 1, b), idx(a + 1, b + 1), idx(a, b + 1)], ++id, 0);
    }
  } else if (dataset === 'UNSTRUCTURED_GRID') {
    if (cells && !types) throw new Error(`${input.name || 'VTK'}: CELLS without CELL_TYPES.`);
    if (cells) cells.forEach((c, i) => addCell(fe, types[i], c, i + 1, fe.coords.a));
  } else if (dataset === 'POLYDATA') {
    let id = 0;
    for (const c of poly.POLYGONS || []) addCell(fe, 7, c, ++id, fe.coords.a);
    for (const c of poly.TRIANGLE_STRIPS || []) addCell(fe, 6, c, ++id, fe.coords.a);
    for (const c of poly.LINES || []) addCell(fe, 4, c, ++id, fe.coords.a);
    if (poly.VERTICES && poly.VERTICES.length) fe.skip('vertex cell', poly.VERTICES.length);
  } else throw new Error(`${input.name || 'VTK'}: DATASET ${dataset} is not supported (expected UNSTRUCTURED_GRID, POLYDATA, STRUCTURED_GRID, RECTILINEAR_GRID or STRUCTURED_POINTS).`);
  if (skippedCellData) warnings.push(`${skippedCellData} cell-data array(s) were not kept (only point data is carried on the models).`);
  if (fe.skipped['vertex cell'] && (fe.blocks.size === 0)) delete fe.skipped['vertex cell'];
  const r = feResult(fe, 'vtk', opts, { name: input.name || 'VTK', expected: 'expected a POINTS section', warnings, skipped: skippedCellData ? { cellDataArrays: skippedCellData } : {}, extra: { encoding: enc.toLowerCase(), dataset, title } });
  attachNodeData(r, pointData, comps);
  return r;
}

// ---------------------------------------------------------------------------------------------------------
// XML (.vtu / .vtp)
// ---------------------------------------------------------------------------------------------------------

const XML_TYPES = {
  Int8: [1, 'getInt8'], UInt8: [1, 'getUint8'], Int16: [2, 'getInt16'], UInt16: [2, 'getUint16'], Int32: [4, 'getInt32'], UInt32: [4, 'getUint32'],
  Int64: [8, 'getBigInt64'], UInt64: [8, 'getBigUint64'], Float32: [4, 'getFloat32'], Float64: [8, 'getFloat64'],
};

function findSeq(bytes, text, from = 0) {
  const n = text.length;
  outer: for (let i = from; i + n <= bytes.length; i++) {
    for (let k = 0; k < n; k++) if (bytes[i + k] !== text.charCodeAt(k)) continue outer;
    return i;
  }
  return -1;
}

/** VTK XML unstructured grid (.vtu) or poly data (.vtp). */
export async function parseVtkXml(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  // the appended block may be raw binary, so the XML is parsed only up to it
  const ap = findSeq(bytes, '<AppendedData');
  let xmlText, appended = null;
  if (ap >= 0) {
    const tagEnd = bytes.indexOf(62, ap);
    const tag = decodeText(bytes.subarray(ap, tagEnd + 1));
    const us = bytes.indexOf(95, tagEnd);
    const m = /encoding\s*=\s*["'](\w+)["']/.exec(tag);
    const encoding = m ? m[1].toLowerCase() : 'base64';
    const close = findSeq(bytes, '</AppendedData', us);
    appended = { encoding, data: bytes.subarray(us + 1, encoding === 'raw' || close < 0 ? bytes.length : close) };
    if (encoding !== 'raw') appended.text = decodeText(appended.data);
    xmlText = decodeText(bytes.subarray(0, ap)) + '</VTKFile>';
  } else xmlText = decodeText(bytes);
  const root = parseXml(xmlText);
  if (root.name !== 'VTKFile') throw new Error(`${input.name || 'VTU'}: not a VTK XML file (expected root element <VTKFile>, found <${root.tag}>).`);
  const type = attr(root, 'type', '');
  if (type !== 'UnstructuredGrid' && type !== 'PolyData') throw new Error(`${input.name || 'VTU'}: VTK XML type "${type}" is not read (supported: UnstructuredGrid .vtu and PolyData .vtp). In ParaView use File → Save Data and choose .vtu, or apply "Merge Blocks" / "Append Datasets" first for parallel or multi-block data.`);
  const little = attr(root, 'byte_order', 'LittleEndian') !== 'BigEndian';
  const H = attr(root, 'header_type', 'UInt32') === 'UInt64' ? 8 : 4;
  const compressor = attr(root, 'compressor', '');
  if (compressor && !/ZLib/i.test(compressor)) throw new Error(`${input.name || 'VTU'}: compressor "${compressor}" is not supported (only vtkZLibDataCompressor). Re-save with ZLib compression or none.`);

  const headNum = (dv, o) => (H === 8 ? Number(dv.getBigUint64(o, little)) : dv.getUint32(o, little));
  const b64len = (n) => Math.ceil(n / 3) * 4;
  /** Raw bytes of one binary array starting at character/byte `start` of `src`. */
  const binaryPayload = async (src, start, isText) => {
    const view = (b) => new DataView(b.buffer, b.byteOffset, b.byteLength);
    if (!compressor) {
      if (!isText) { const n = headNum(view(src), start); return src.subarray(start + H, start + H + n); }
      const hd = base64Decode(src, start, start + b64len(H));
      const n = headNum(view(hd), 0);
      // some writers encode the length header on its own (padded), most encode header + data as one stream
      if (src[start + b64len(H) - 1] === '=') return base64Decode(src, start + b64len(H), start + b64len(H) + b64len(n));
      return base64Decode(src, start, start + b64len(H + n)).subarray(H, H + n);
    }
    let head, dataStart;
    if (isText) {
      const h3 = base64Decode(src, start, start + b64len(3 * H));
      const nb = headNum(view(h3), 0);
      const hl = b64len((3 + nb) * H);
      head = base64Decode(src, start, start + hl);
      dataStart = start + hl;
    } else { const nb = headNum(view(src), start); head = src.subarray(start, start + (3 + nb) * H); dataStart = start + (3 + nb) * H; }
    const hv = view(head);
    const nb = headNum(hv, 0), bs = headNum(hv, H), last = headNum(hv, 2 * H);
    let total = 0;
    const sizes = [];
    for (let i = 0; i < nb; i++) { sizes.push(headNum(hv, (3 + i) * H)); total += sizes[i]; }
    const comp = isText ? base64Decode(src, dataStart, dataStart + b64len(total)) : src.subarray(dataStart, dataStart + total);
    const out = new Uint8Array(nb ? (nb - 1) * bs + (last || bs) : 0);
    let ci = 0, oi = 0;
    for (let i = 0; i < nb; i++) { const blk = await inflateZlib(comp.subarray(ci, ci + sizes[i]).slice()); out.set(blk, oi); oi += blk.length; ci += sizes[i]; }
    return out.subarray(0, oi);
  };
  const readArray = async (el) => {
    if (!el) return null;
    const fmt = attr(el, 'format', 'ascii');
    if (fmt === 'ascii') return numberList(el.text);
    const t = XML_TYPES[attr(el, 'type', 'Float32')];
    if (!t) throw new Error(`${input.name || 'VTU'}: DataArray type "${attr(el, 'type')}" is not supported.`);
    let raw;
    if (fmt === 'binary') { const s = el.text.replace(/\s+/g, ''); raw = await binaryPayload(s, 0, true); } else if (fmt === 'appended') {
      if (!appended) throw new Error(`${input.name || 'VTU'}: a DataArray refers to appended data but the file has no <AppendedData> block.`);
      const off = Number(attr(el, 'offset', 0));
      raw = appended.encoding === 'raw' ? await binaryPayload(appended.data, off, false) : await binaryPayload(appended.text, off, true);
    } else throw new Error(`${input.name || 'VTU'}: DataArray format "${fmt}" is not supported.`);
    const n = Math.floor(raw.length / t[0]);
    const out = new Float64Array(n);
    const dv = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);
    const big = t[1].includes('Big');
    for (let i = 0, o = 0; i < n; i++, o += t[0]) out[i] = big ? Number(dv[t[1]](o, little)) : dv[t[1]](o, little);
    return out;
  };
  const named = (parent, nm) => (parent ? kids(parent, 'DataArray').find((d) => attr(d, 'Name') === nm) : null);

  const fe = new FeMesh(name);
  const warnings = [];
  const pieces = desc(root, 'Piece');
  if (!pieces.length) throw new Error(`${input.name || 'VTU'}: no <Piece> element found in the ${type} file.`);
  const pointData = {}, comps = {};
  let base = 0, eid = 0, cellArrays = 0;
  for (const pc of pieces) {
    const pts = kid(pc, 'Points');
    const c = await readArray(pts && kid(pts, 'DataArray'));
    if (!c) continue;
    const n = c.length / 3;
    for (let i = 0; i < n; i++) fe.addNode(base + i, c[3 * i], c[3 * i + 1], c[3 * i + 2]);
    const xyz = () => fe.coords.a;
    const walk = async (parent, typeFor) => {
      const conn = await readArray(named(parent, 'connectivity')), offs = await readArray(named(parent, 'offsets'));
      if (!conn || !offs) return;
      const types = typeof typeFor === 'number' ? null : await readArray(named(parent, 'types'));
      let s = 0;
      for (let i = 0; i < offs.length; i++) {
        const ids = new Array(offs[i] - s);
        for (let k = s; k < offs[i]; k++) ids[k - s] = conn[k] + base;
        s = offs[i];
        addCell(fe, types ? types[i] : typeFor, ids, ++eid, xyz());
      }
    };
    if (type === 'UnstructuredGrid') await walk(kid(pc, 'Cells'), null);
    else { await walk(kid(pc, 'Polys'), 7); await walk(kid(pc, 'Strips'), 6); await walk(kid(pc, 'Lines'), 4); await walk(kid(pc, 'Verts'), 1); }
    const pd = kid(pc, 'PointData');
    for (const d of pd ? kids(pd, 'DataArray') : []) {
      const nm = attr(d, 'Name', `array ${Object.keys(pointData).length + 1}`), nc = Number(attr(d, 'NumberOfComponents', 1)) || 1;
      const vals = await readArray(d);
      if (vals.length !== n * nc) continue;
      if (!pointData[nm]) { pointData[nm] = []; comps[nm] = nc; }
      pointData[nm].push({ base, vals });
    }
    const cd = kid(pc, 'CellData');
    if (cd) cellArrays += kids(cd, 'DataArray').length;
    base += n;
  }
  const merged = {};
  for (const nm in pointData) { const a = new Float64Array(base * comps[nm]); for (const part of pointData[nm]) a.set(part.vals, part.base * comps[nm]); merged[nm] = a; }
  if (cellArrays) warnings.push(`${cellArrays} cell-data array(s) were not kept (only point data is carried on the models).`);
  if (fe.skipped['vertex cell'] && fe.blocks.size === 0) delete fe.skipped['vertex cell'];
  const r = feResult(fe, 'vtu', opts, { name: input.name || 'VTU', expected: 'expected <Points> with a DataArray', warnings, skipped: cellArrays ? { cellDataArrays: cellArrays } : {}, extra: { dataset: type, pieces: pieces.length, compressor: compressor || null, appended: appended ? appended.encoding : null } });
  attachNodeData(r, merged, comps);
  return r;
}
