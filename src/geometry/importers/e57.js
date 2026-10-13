// ASTM E57 (ASTM E2807) — the vendor-neutral exchange file of terrestrial and mobile laser scanners (".e57").
//
// Layout of the file:
//   header     48 bytes: "ASTM-E57", version major / minor (uint32), file length, physical offset and logical
//              length of the XML section, page size (uint64; always 1024). All integers are little-endian.
//   pages      the whole file, header included, is a sequence of 1024-byte pages: 1020 bytes of data followed by
//              the CRC-32C (Castagnoli) of those bytes, stored most significant byte first. "Logical" positions
//              count data bytes only; the offsets in the header and in the XML are physical.
//   XML        <e57Root> with data3D (a vector of scans) and images2D. A scan has guid, name, pose (rotation
//              quaternion w x y z and translation, scan → file coordinates), intensityLimits, colorLimits and
//              points: a CompressedVector with fileOffset, recordCount and a prototype that lists the fields of a
//              record — Float (precision single | double), Integer (minimum, maximum) or ScaledInteger (minimum,
//              maximum, scale, offset; value = raw · scale + offset).
//   binary     one section per CompressedVector: a 32-byte header (id 1, logical length, physical offsets of the
//              first data packet and of the index), then packets of at most 64 kB, each a multiple of 4 bytes:
//                data packet   type 1, flags, length − 1 (uint16), bytestream count (uint16), one uint16 length
//                              per bytestream, then the bytestreams — one per prototype field, in prototype order
//                index packet  type 0 (record number → packet offset; not needed for sequential reading)
//                empty packet  type 2 (padding)
//              A bytestream is the field's values end to end: floats as IEEE 754, integers as (value − minimum)
//              in ceil(log2(maximum − minimum + 1)) bits, least significant bit first, with no padding between
//              values (a field with minimum = maximum takes no bits). A value may straddle two packets, and the
//              packets need not hold the same number of values of every field, so each field is decoded into
//              its own queue and points are formed from what all queues hold.
// The description follows ASTM E2807 as documented by libE57Format (https://github.com/asmaloney/libE57Format,
// Boost Software License 1.0: its format notes and the packet and codec definitions in its source); dropping
// the points with a non-zero invalid state is what pye57 (https://github.com/davidcaron/pye57, MIT) does, and
// pye57 served as the independent check of the results. The reader below is written for this application from
// that description.
//
// Read: cartesianX/Y/Z, or sphericalRange / Azimuth / Elevation when a scan has no cartesian fields (x = r cos e
// cos a, y = r cos e sin a, z = r sin e); cartesianInvalidState / sphericalInvalidState (points with a non-zero
// state are dropped and counted, as are points with a non-finite coordinate); intensity (kept as written, see
// meta.scans[].intensityLimits for its range); colorRed/Green/Blue scaled to 0–255 by colorLimits (or by the
// range of the field). The pose of every scan is applied and all scans are merged into one `points` model;
// meta.scans lists them. Lengths are metres by definition. A projected system in coordinateMetadata (WKT)
// is reported as the CRS. Clouds larger than `opts.maxPoints` (default 5 000 000) are thinned by keeping every
// n-th valid point. The packets are decoded one at a time; nothing of the size of the file is copied.
//
// Not read: images2D (counted), every other point field — rowIndex, columnIndex, returnIndex, returnCount,
// timeStamp, isIntensityInvalid, isColorInvalid, isTimeStampInvalid, normals and other extension fields (their
// bytestreams are skipped; names in meta.ignoredFields), pointGroupingSchemes, the bounds and the index packets.
// Refused, per scan, with the reason in the warning (and as the error when no scan can be read): a coordinate,
// colour or intensity field of a type other than Float / Integer / ScaledInteger, integer fields wider than 64
// bits or with maximum < minimum, and a non-empty codecs list (extension codecs). A page of the XML section
// that fails its checksum is an error; data pages that fail are counted in a warning and read as they are.

import { inputBytes, decodeText, baseName, cleanName, M4 } from './common.js';
import { parseXml, kids, kid, attr, kidText } from '../xml.js';
import { parsePrj } from '../units.js';
import { cloudResult } from './pointcloud.js';

const DEFAULT_MAX_POINTS = 5000000;
const PAGE = 1024, DATA = 1020;
const USED = ['cartesianX', 'cartesianY', 'cartesianZ', 'cartesianInvalidState', 'sphericalRange', 'sphericalAzimuth', 'sphericalElevation', 'sphericalInvalidState', 'intensity', 'colorRed', 'colorGreen', 'colorBlue'];

/** A problem with one scan: reported as a warning when other scans can be read, as the error otherwise. */
class ScanError extends Error {}
const bad = (msg) => { throw new ScanError(msg); };

// CRC-32C (polynomial 0x1EDC6F41, reflected), four bytes per step
const CRC = (() => {
  const t = new Int32Array(1024);
  for (let i = 0; i < 256; i++) { let c = i; for (let k = 0; k < 8; k++) c = c & 1 ? (c >>> 1) ^ 0x82f63b78 : c >>> 1; t[i] = c; }
  for (let i = 0; i < 256; i++) for (let k = 1; k < 4; k++) { const p = t[(k - 1) * 256 + i]; t[k * 256 + i] = (p >>> 8) ^ t[p & 255]; }
  return t;
})();
function crc32c(b, off, len) {
  let c = -1, i = off;
  for (const e = off + (len & ~3); i < e; i += 4) {
    c ^= b[i] | (b[i + 1] << 8) | (b[i + 2] << 16) | (b[i + 3] << 24);
    c = CRC[768 + (c & 255)] ^ CRC[512 + ((c >>> 8) & 255)] ^ CRC[256 + ((c >>> 16) & 255)] ^ CRC[c >>> 24];
  }
  for (const e = off + len; i < e; i++) c = (c >>> 8) ^ CRC[(c ^ b[i]) & 255];
  return ~c >>> 0;
}

/** Unsigned 64-bit little-endian integer as a Number; null when it exceeds 2^53 (no real offset does). */
function u64(b, o) {
  const hi = (b[o + 4] | (b[o + 5] << 8) | (b[o + 6] << 16) | (b[o + 7] << 24)) >>> 0;
  return hi >= 0x200000 ? null : hi * 4294967296 + ((b[o] | (b[o + 1] << 8) | (b[o + 2] << 16) | (b[o + 3] << 24)) >>> 0);
}

/** Page-aware access: copies logical byte ranges out of the paged file and checks each page once. */
class Pages {
  constructor(bytes) {
    this.b = bytes;
    this.seen = new Uint8Array(Math.ceil(bytes.length / PAGE));
    this.checked = 0; this.failed = 0;
  }
  /** Logical position of a physical offset. */
  logical(phys, what) {
    const r = phys % PAGE;
    if (!(phys >= 0 && phys < this.b.length) || r >= DATA) bad(`${what} is at offset ${phys}, which is ${r >= DATA ? 'inside a page checksum' : `outside the file (${this.b.length} bytes)`} — the file is truncated or corrupt.`);
    return (phys - r) / PAGE * DATA + r;
  }
  /** Copy `len` logical bytes from logical position `pos` into `dst`. */
  read(pos, len, dst, what) {
    const b = this.b;
    let page = Math.floor(pos / DATA), r = pos - page * DATA, o = 0;
    while (o < len) {
      const p0 = page * PAGE;
      if (p0 + PAGE > b.length) bad(`${what} runs past the end of the file (${b.length} bytes) — the file is truncated or corrupt.`);
      if (!this.seen[page]) {
        this.seen[page] = 1; this.checked++;
        if (crc32c(b, p0, DATA) !== ((b[p0 + DATA] << 24) | (b[p0 + DATA + 1] << 16) | (b[p0 + DATA + 2] << 8) | b[p0 + DATA + 3]) >>> 0) this.failed++;
      }
      const n = Math.min(DATA - r, len - o);
      dst.set(b.subarray(p0 + r, p0 + r + n), o);
      o += n; r = 0; page++;
    }
    return dst;
  }
}

/** Value of a numeric child element (an empty element is 0; a ScaledInteger element holds the raw integer). */
function num(node, name, dflt = 0) {
  const c = node ? kid(node, name) : null, t = c ? c.text.trim() : '';
  let v = t === '' ? (c ? 0 : dflt) : Number(t);
  if (c && attr(c, 'type') === 'ScaledInteger') v = v * Number(attr(c, 'scale', '1')) + Number(attr(c, 'offset', '0'));
  return Number.isFinite(v) ? v : dflt;
}

function bigAttr(node, name, field) {
  const s = String(attr(node, name, '0')).trim();
  if (!/^[+-]?\d{1,30}(\.0*)?$/.test(s)) bad(`field ${field}: ${name}="${s.slice(0, 40)}" is not an integer.`);
  return BigInt(s.replace(/\..*$/, ''));
}

/** Terminal fields of a prototype in bytestream order (nested structures are numbered depth-first). */
function protoFields(node, out, depth) {
  if (depth > 8) bad('the prototype is nested more than 8 levels deep.');
  for (const c of node.children) {
    const type = attr(c, 'type', '');
    if (type === 'Structure') { protoFields(c, out, depth + 1); continue; }
    const f = { name: cleanName(c.tag, '?').slice(0, 60), type: String(type).slice(0, 30), kind: '', bits: 0, min: 0, minBig: 0n, big: false, scale: 1, offset: 0, lo: null, hi: null, error: null, q: null, qn: 0, got: 0, carry: new Uint8Array(16), carryN: 0, bit: 0, constant: false };
    try {
      if (type === 'Float') {
        const p = attr(c, 'precision', 'double');
        if (p !== 'single' && p !== 'double') bad(`field ${f.name}: Float precision "${p}" is not known (expected single or double).`);
        f.kind = p === 'single' ? 'f32' : 'f64'; f.bits = p === 'single' ? 32 : 64;
        const lo = Number(attr(c, 'minimum', 'x')), hi = Number(attr(c, 'maximum', 'x')); // optional bounds
        if (hi > lo && hi - lo < 1e30) { f.lo = lo; f.hi = hi; }
      } else if (type === 'Integer' || type === 'ScaledInteger') {
        const lo = bigAttr(c, 'minimum', f.name), hi = bigAttr(c, 'maximum', f.name);
        if (hi < lo) bad(`field ${f.name}: maximum (${hi}) is below minimum (${lo}).`);
        f.bits = hi === lo ? 0 : (hi - lo).toString(2).length;
        if (f.bits > 64) bad(`field ${f.name}: the range ${lo} … ${hi} needs ${f.bits} bits per value; E57 integers have at most 64.`);
        f.kind = 'int'; f.minBig = lo; f.min = Number(lo);
        // values that a double cannot hold exactly are assembled with BigInt
        f.big = f.bits > 53 || !Number.isSafeInteger(Number(lo)) || !Number.isSafeInteger(Number(hi));
        if (type === 'ScaledInteger') {
          f.scale = Number(attr(c, 'scale', '1')); f.offset = Number(attr(c, 'offset', '0'));
          if (!Number.isFinite(f.scale) || !Number.isFinite(f.offset)) bad(`field ${f.name}: scale="${String(attr(c, 'scale', '1')).slice(0, 30)}" and offset="${String(attr(c, 'offset', '0')).slice(0, 30)}" are not both numbers.`);
        }
        f.lo = Number(lo) * f.scale + f.offset; f.hi = Number(hi) * f.scale + f.offset;
        f.constant = f.bits === 0;
      } else bad(`field ${f.name}: element type "${type}" is not one this reader decodes (Float, Integer, ScaledInteger).`);
    } catch (e) { if (!(e instanceof ScanError)) throw e; f.error = e.message; f.bits = 0; }
    out.push(f);
  }
  return out;
}

/** Append the values of one bytestream piece to the field's queue (at most `want` values). */
function feed(f, src, off, len, want, work, wdv) {
  const total = f.carryN + len;
  work.set(f.carry.subarray(0, f.carryN));
  work.set(src.subarray(off, off + len), f.carryN);
  work.fill(0, total, total + 8);
  const bits = f.bits;
  let n = Math.floor((total * 8 - f.bit) / bits);
  if (n > want) n = want;
  if (f.qn + n > f.q.length) { let cap = f.q.length * 2; while (cap < f.qn + n) cap *= 2; const q = new Float64Array(cap); q.set(f.q.subarray(0, f.qn)); f.q = q; }
  const q = f.q;
  let k = f.qn, bp = f.bit;
  if (f.kind === 'f32') { for (let i = 0; i < n; i++) q[k++] = wdv.getFloat32(4 * i, true); bp = 32 * n; } else if (f.kind === 'f64') { for (let i = 0; i < n; i++) q[k++] = wdv.getFloat64(8 * i, true); bp = 64 * n; } else {
    const w = work, min = f.min, sc = f.scale, of = f.offset;
    if (bits <= 24) {
      const mask = (1 << bits) - 1;
      for (let j = 0; j < n; j++, bp += bits) { const i = bp >>> 3; q[k++] = ((((w[i] | (w[i + 1] << 8) | (w[i + 2] << 16) | (w[i + 3] << 24)) >>> (bp & 7)) & mask) + min) * sc + of; }
    } else {
      // 32 bits starting at any bit position
      const get = (p) => { const i = p >>> 3, s = p & 7, lo = (w[i] | (w[i + 1] << 8) | (w[i + 2] << 16) | (w[i + 3] << 24)) >>> s; return s ? (lo | (w[i + 4] << (32 - s))) >>> 0 : lo; };
      if (bits <= 32) {
        const mask = bits === 32 ? 0xffffffff : 2 ** bits - 1;
        for (let j = 0; j < n; j++, bp += bits) q[k++] = (((get(bp) & mask) >>> 0) + min) * sc + of;
      } else {
        const hmask = bits === 64 ? 0xffffffff : 2 ** (bits - 32) - 1;
        if (!f.big) for (let j = 0; j < n; j++, bp += bits) q[k++] = (((get(bp + 32) & hmask) >>> 0) * 4294967296 + get(bp) + min) * sc + of;
        else { const mb = f.minBig; for (let j = 0; j < n; j++, bp += bits) q[k++] = Number(((BigInt((get(bp + 32) & hmask) >>> 0) << 32n) | BigInt(get(bp))) + mb) * sc + of; }
      }
    }
  }
  f.qn = k; f.got += n;
  if (n === want) { f.carryN = 0; f.bit = 0; return; } // complete: whatever follows is padding
  const used = bp >>> 3;
  f.bit = bp & 7; f.carryN = total - used;
  f.carry.set(work.subarray(used, total));
}

/** ASTM E57 point cloud. */
export function parseE57(input, opts = {}) {
  const bytes = inputBytes(input);
  const fname = input.name || 'E57';
  if (bytes.length < 48 || decodeText(bytes.subarray(0, 8)) !== 'ASTM-E57') throw new Error(`${fname}: not an E57 file (expected the signature "ASTM-E57" and a 48-byte header).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const vMajor = dv.getUint32(8, true), vMinor = dv.getUint32(12, true);
  const fileLen = u64(bytes, 16), xmlOff = u64(bytes, 24), xmlLen = u64(bytes, 32), pageSize = u64(bytes, 40);
  if (vMajor !== 1) throw new Error(`${fname}: E57 format version ${vMajor}.${vMinor} is not known (expected 1.0).`);
  if (pageSize !== PAGE) throw new Error(`${fname}: the header gives a page size of ${pageSize ?? 'more than 2^53'} bytes; E57 files use 1024 — the header is corrupt.`);
  const warnings = [];
  const short = fileLen === null || fileLen > bytes.length;
  if (xmlOff === null || xmlLen === null) throw new Error(`${fname}: the header places the XML section beyond 2^53 bytes — the header is corrupt.`);
  if (!xmlOff || !xmlLen) throw new Error(`${fname}: the header gives no XML section (offset ${xmlOff}, length ${xmlLen}) — the file was not closed by the program that wrote it.`);
  if (xmlOff + xmlLen > bytes.length || xmlOff < 48) throw new Error(`${fname}: the XML section (offset ${xmlOff}, ${xmlLen} bytes) lies outside the file (${bytes.length} bytes${short ? `; the header declares ${fileLen ?? 'more than 2^53'}` : ''}) — the file is truncated or corrupt.`);
  if (short) warnings.push(`The header declares a file of ${fileLen ?? 'more than 2^53'} bytes but ${bytes.length} were supplied; scans that lie beyond the end are reported below.`);

  const pages = new Pages(bytes);
  let root;
  try {
    const xml = pages.read(pages.logical(xmlOff, 'the XML section'), xmlLen, new Uint8Array(xmlLen), 'the XML section');
    if (pages.failed) bad(`${pages.failed} of the ${pages.checked} pages of the XML section fail their CRC-32C checksum — the file is corrupt.`);
    try { root = parseXml(decodeText(xml)); } catch { bad('the XML section holds no XML element — the file is corrupt.'); }
    if (root.name !== 'e57Root') bad(`the XML section starts with <${String(root.tag).slice(0, 40)}> where <e57Root> was expected.`);
  } catch (e) { if (e instanceof ScanError) throw new Error(`${fname}: ${e.message}`); throw e; }
  const xmlPages = pages.checked;

  const d3 = kid(root, 'data3D'), i2 = kid(root, 'images2D');
  const scanNodes = d3 ? kids(d3, 'vectorChild') : [];
  const nImages = i2 ? kids(i2, 'vectorChild').length : 0;
  if (!scanNodes.length) throw new Error(`${fname}: the file holds no scans (data3D is ${d3 ? 'empty' : 'missing'}${nImages ? `; it holds ${nImages} image(s), which are not read` : ''}).`);

  // ---- pass 1: prototypes and section headers of all scans, so that the thinning stride is known up front
  const head = new Uint8Array(32);
  const scans = [], problems = [], ignored = {};
  let declared = 0, anyI = false, anyRgb = false;
  scanNodes.forEach((node, si) => {
    const name = cleanName(kidText(node, 'name', ''), `scan ${si + 1}`).slice(0, 80);
    const pts = kid(node, 'points');
    const s = { node, label: `Scan ${si + 1}${kidText(node, 'name', '') ? ` ("${name}")` : ''}`, ok: false, meta: { name, guid: kidText(node, 'guid', '') || null, records: 0, points: 0, invalid: 0, coordinates: null, fields: [], pose: null, description: kidText(node, 'description', '') || null, sensor: [kidText(node, 'sensorVendor', ''), kidText(node, 'sensorModel', '')].filter(Boolean).join(' ') || null } };
    scans.push(s);
    try {
      if (!pts || attr(pts, 'type') !== 'CompressedVector') bad('it has no points element of type CompressedVector.');
      const proto = kid(pts, 'prototype');
      if (!proto || !proto.children.length) bad('its points have no prototype (the list of fields of a record).');
      const count = Number(attr(pts, 'recordCount', '0')), off = Number(attr(pts, 'fileOffset', ''));
      if (!(count >= 0) || !Number.isSafeInteger(count)) bad(`recordCount="${String(attr(pts, 'recordCount')).slice(0, 40)}" is not a count.`);
      const fields = protoFields(proto, [], 0);
      s.fields = fields;
      s.meta.fields = fields.map((f) => f.name);
      const by = {};
      for (const f of fields) if (USED.includes(f.name) && !by[f.name]) by[f.name] = f; else ignored[f.name] = (ignored[f.name] || 0) + 1;
      s.by = by;
      s.cart = !!(by.cartesianX && by.cartesianY && by.cartesianZ);
      const sph = !!(by.sphericalRange && by.sphericalAzimuth && by.sphericalElevation);
      if (!s.cart && !sph) bad(`its records have no complete coordinate triple (fields: ${s.meta.fields.join(', ').slice(0, 200)}).`);
      s.meta.coordinates = s.cart ? 'cartesian' : 'spherical';
      s.rgb = !!(by.colorRed && by.colorGreen && by.colorBlue);
      s.used = (s.cart ? [by.cartesianX, by.cartesianY, by.cartesianZ, by.cartesianInvalidState] : [by.sphericalRange, by.sphericalAzimuth, by.sphericalElevation, by.sphericalInvalidState]).concat(by.intensity, s.rgb ? [by.colorRed, by.colorGreen, by.colorBlue] : []).filter(Boolean);
      for (const f of s.used) if (f.error) bad(f.error);
      const codecs = kid(pts, 'codecs');
      if (codecs && codecs.children.length) bad(`its points name ${codecs.children.length} codec(s) in the codecs list; only the standard encoding (bit-packed integers, IEEE floats) is read, not E57 extension codecs.`);
      s.meta.records = count;
      if (!count) { s.ok = true; s.count = 0; return; }
      if (!Number.isSafeInteger(off)) bad(`fileOffset="${String(attr(pts, 'fileOffset')).slice(0, 40)}" is not an offset.`);
      const start = pages.logical(off, 'its binary section');
      pages.read(start, 32, head, 'the header of its binary section');
      const secLen = u64(head, 8), dataOff = u64(head, 16);
      if (head[0] !== 1 || secLen === null || dataOff === null) bad(`its binary section at offset ${off} does not start with a CompressedVector header (id ${head[0]}) — the file is corrupt.`);
      s.end = start + Math.min(secLen, bytes.length);
      s.pos = pages.logical(dataOff, 'its first data packet');
      if (s.pos < start + 32 || s.pos > s.end) bad(`its first data packet (offset ${dataOff}) lies outside its binary section — the file is corrupt.`);
      // a record cannot take less than the bits of its fields: header counts larger than the data are not trusted
      const recBits = fields.reduce((a, f) => a + f.bits, 0);
      if (!recBits && fields.every((f) => !f.error)) bad(`all fields of its ${count} records are constants (minimum = maximum); there is nothing to read.`);
      const room = recBits ? Math.floor((s.end - s.pos) * 8 / recBits) : count;
      s.count = Math.min(count, room);
      if (s.count < count) warnings.push(`${s.label} declares ${count} records but its binary section has room for ${room}; read ${s.count}.`);
      s.ok = true;
      declared += s.count;
      if (by.intensity) anyI = true;
      if (s.rgb) anyRgb = true;
    } catch (e) { if (!(e instanceof ScanError)) throw e; s.problem = true; problems.push(`${s.label} was not read: ${e.message}`); }
  });

  const maxPts = opts.maxPoints ?? DEFAULT_MAX_POINTS;
  const stride = maxPts > 0 && declared > maxPts ? Math.ceil(declared / maxPts) : 1;
  const cap = Math.ceil(declared / stride);
  const coords = new Float64Array(cap * 3);
  const intensity = anyI ? new Float32Array(cap) : null, rgb = anyRgb ? new Uint8ClampedArray(cap * 3) : null;
  let n = 0, valid = 0, invalid = 0, noI = 0, noRgb = 0;

  // ---- pass 2: packets
  const pkt = new Uint8Array(65536), work = new Uint8Array(65536 + 32), wdv = new DataView(work.buffer);
  for (const s of scans) {
    if (!s.ok) continue;
    const sm = s.meta, by = s.by, count = s.count, nStreams = s.fields.length;
    // pose: scan → file coordinates
    const pose = kid(s.node, 'pose');
    let M = null;
    if (pose) {
      const rot = kid(pose, 'rotation'), tr = kid(pose, 'translation');
      let q = [num(rot, 'x'), num(rot, 'y'), num(rot, 'z'), num(rot, 'w', rot ? 0 : 1)];
      const t = [num(tr, 'x'), num(tr, 'y'), num(tr, 'z')];
      const l = Math.hypot(q[0], q[1], q[2], q[3]);
      if (!(l > 0)) { warnings.push(`${s.label}: the rotation of its pose is a zero quaternion; no rotation was applied.`); q = [0, 0, 0, 1]; } else q = q.map((v) => v / l);
      sm.pose = { rotation: [q[3], q[0], q[1], q[2]], translation: t };
      M = M4.fromTRS(t, q);
      if (M4.isIdentity(M)) M = null;
    }
    const lim = (parent, a, b, f) => {
      const L = kid(s.node, parent);
      if (L && kid(L, a) && kid(L, b)) return [num(L, a), num(L, b)];
      if (f.kind === 'int') return [f.lo, f.hi];
      if (f.lo === null) return null;
      // bounds of a Float field are often those of the data: for colour only the nominal range is taken from them
      return parent === 'colorLimits' ? [0, f.hi <= 1 ? 1 : f.hi <= 255 ? 255 : 65535] : [f.lo, f.hi];
    };
    if (by.intensity) { const l = lim('intensityLimits', 'intensityMinimum', 'intensityMaximum', by.intensity); if (l) sm.intensityLimits = l; }
    // colour → 0–255 by colorLimits, else by the range of the field, else 0–1 floats
    const cs = [0, 0, 0], c0 = [0, 0, 0];
    if (s.rgb) ['Red', 'Green', 'Blue'].forEach((c, k) => { const l = lim('colorLimits', `color${c}Minimum`, `color${c}Maximum`, by[`color${c}`]) || [0, 1]; c0[k] = l[0]; cs[k] = l[1] > l[0] ? 255 / (l[1] - l[0]) : 1; });
    const F = s.cart ? [by.cartesianX, by.cartesianY, by.cartesianZ] : [by.sphericalRange, by.sphericalAzimuth, by.sphericalElevation];
    const fInv = s.cart ? by.cartesianInvalidState : by.sphericalInvalidState, fI = by.intensity, fR = by.colorRed, fG = by.colorGreen, fB = by.colorBlue;
    for (const f of s.used) { f.q = new Float64Array(4096); if (f.constant) f.q.fill(f.lo); }
    const live = s.used.filter((f) => !f.constant);
    const n0 = n;
    let done = 0;

    const emit = (m) => {
      for (const f of s.used) if (f.constant && f.q.length < m) f.q = new Float64Array(m).fill(f.lo);
      const A = F[0].q, B = F[1].q, C = F[2].q, V = fInv ? fInv.q : null, I = fI ? fI.q : null, R = s.rgb ? fR.q : null, G = s.rgb ? fG.q : null, Bl = s.rgb ? fB.q : null, cart = s.cart;
      let k = n, nv = valid, dropped = 0;
      for (let i = 0; i < m; i++) {
        if (V !== null && V[i] !== 0) { dropped++; continue; }
        let x, y, z;
        if (cart) { x = A[i]; y = B[i]; z = C[i]; } else { const r = A[i], az = B[i], el = C[i], h = r * Math.cos(el); x = h * Math.cos(az); y = h * Math.sin(az); z = r * Math.sin(el); }
        if (!(x - x === 0 && y - y === 0 && z - z === 0)) { dropped++; continue; }
        if (nv++ % stride) continue;
        const o = 3 * k;
        if (M !== null) { coords[o] = M[0] * x + M[4] * y + M[8] * z + M[12]; coords[o + 1] = M[1] * x + M[5] * y + M[9] * z + M[13]; coords[o + 2] = M[2] * x + M[6] * y + M[10] * z + M[14]; } else { coords[o] = x; coords[o + 1] = y; coords[o + 2] = z; }
        if (I !== null) intensity[k] = I[i];
        if (R !== null) { rgb[o] = (R[i] - c0[0]) * cs[0]; rgb[o + 1] = (G[i] - c0[1]) * cs[1]; rgb[o + 2] = (Bl[i] - c0[2]) * cs[2]; } // the clamped array rounds
        k++;
      }
      sm.invalid += dropped; sm.points += nv - valid;
      n = k; valid = nv;
      for (const f of live) { f.q.copyWithin(0, m, f.qn); f.qn -= m; }
      done += m;
    };

    try {
      if (!live.length) emit(count);
      let pos = s.pos;
      while (done < count) {
        const left = `after ${done} of its ${count} records`;
        if (pos + 4 > s.end) bad(`its binary section ends ${left} — the file is truncated or corrupt.`);
        pages.read(pos, 4, pkt, 'a packet header');
        const type = pkt[0], plen = (pkt[2] | (pkt[3] << 8)) + 1;
        if (type > 2) bad(`a packet of unknown type ${type} follows ${left} — the file is corrupt.`);
        if (plen % 4 || pos + plen > s.end) bad(`a packet of ${plen} bytes ${plen % 4 ? '(not a multiple of 4)' : 'that runs past the end of the binary section'} follows ${left} — the file is corrupt.`);
        if (type === 1) {
          pages.read(pos, plen, pkt, 'a data packet');
          const ns = pkt[4] | (pkt[5] << 8);
          if (ns !== nStreams || 6 + 2 * ns > plen) bad(`a data packet holds ${ns} bytestreams where the prototype has ${nStreams} fields (${left}) — the file is corrupt.`);
          let o = 6 + 2 * ns, m = count - done;
          for (let k = 0; k < ns; k++) {
            const len = pkt[6 + 2 * k] | (pkt[7 + 2 * k] << 8), f = s.fields[k];
            if (o + len > plen) bad(`the bytestreams of a data packet are longer than the packet (${left}) — the file is corrupt.`);
            if (f.q !== null && !f.constant && f.got < count) feed(f, pkt, o, len, count - f.got, work, wdv);
            o += len;
          }
          for (const f of live) if (f.qn < m) m = f.qn;
          if (m > 0) emit(m);
        }
        pos += plen;
      }
    } catch (e) {
      if (!(e instanceof ScanError)) throw e;
      s.problem = true;
      problems.push(`${s.label}: ${e.message}${done ? ` The ${done} records before that were kept.` : ''}`);
    }
    for (const f of s.used) f.q = null;
    invalid += sm.invalid;
    // scans without a colour or intensity of their own, merged with scans that have one
    if (rgb && !s.rgb && n > n0) { rgb.fill(160, 3 * n0, 3 * n); noRgb++; }
    if (intensity && !fI && n > n0) noI++;
  }

  if (!n) {
    if (problems.length) throw new Error(`${fname}: ${problems[0]}${problems.length > 1 ? ` (and ${problems.length - 1} more scan(s) with problems)` : ''}`);
    if (invalid) throw new Error(`${fname}: all ${invalid} points of the ${scans.length} scan(s) are flagged invalid or have non-finite coordinates.`);
  }
  for (const p of problems) warnings.push(p);
  const dataPages = pages.checked - xmlPages;
  if (pages.failed) warnings.push(`${pages.failed} of the ${dataPages} data pages read fail their CRC-32C checksum; the points stored in them may be wrong.`);
  if (invalid) warnings.push(`${invalid} point(s) flagged invalid by the scanner (or with non-finite coordinates) were dropped.`);
  if (nImages) warnings.push(`${nImages} image(s) (images2D) were not read.`);
  if (noRgb) warnings.push(`${noRgb} of the ${scans.length} scans carry no colour; their points are grey in attrs.rgb.`);
  if (noI) warnings.push(`${noI} of the ${scans.length} scans carry no intensity; their points have intensity 0.`);
  if (stride > 1) warnings.push(`The cloud has ${valid} points; every ${stride}${stride === 2 ? 'nd' : stride === 3 ? 'rd' : 'th'} point was kept (${n}) to stay under maxPoints.`);

  const attrs = {};
  if (intensity) attrs.intensity = n < cap ? intensity.slice(0, n) : intensity;
  if (rgb) attrs.rgb = n < cap ? new Uint8Array(rgb.subarray(0, 3 * n)) : new Uint8Array(rgb.buffer);
  // E57 lengths are metres whatever unit the coordinate system names
  const wkt = kidText(root, 'coordinateMetadata', '');
  let crs = null;
  if (wkt) { const info = parsePrj(wkt); if (info.kind === 'projected' || info.kind === 'local') crs = { ...info, unit: 'm', unitScale: 1 }; }
  const failedScans = scans.filter((s) => s.problem).length;
  return cloudResult(baseName(input.name), n < cap ? coords.subarray(0, 3 * n) : coords, attrs, opts, {
    format: 'e57', fileName: fname, crs, unitHint: 'm', warnings, expected: `the ${scans.length} scan(s) hold no point records`,
    skipped: { thinnedPoints: valid - n, invalidPoints: invalid, images2D: nImages, scans: failedScans },
    extra: { version: `${vMajor}.${vMinor}`, pointCount: valid, thinning: stride, scans: scans.map((s) => s.meta), attributes: Object.keys(attrs), ignoredFields: ignored, guid: kidText(root, 'guid', '') || null, library: kidText(root, 'e57LibraryVersion', '') || null, coordinateMetadata: wkt || null, checksums: { pages: pages.checked, failed: pages.failed } },
  });
}
