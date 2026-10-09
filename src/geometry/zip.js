// ZIP container reader (central directory; stored and deflate entries) plus the inflate helpers the importers
// share (3MF, KMZ, zipped shapefiles, zipped multi-file uploads, VTK XML zlib blocks, TIFF Deflate).
// Decompression uses the platform DecompressionStream when present (browsers, workers, Node ≥ 18) and falls
// back to a compact pure-JavaScript inflate otherwise.
//
// Not supported: encrypted entries, ZIP64 archives (> 4 GiB or > 65 535 entries), multi-disk archives,
// compression methods other than 0 (stored) and 8 (deflate).

const SIG_EOCD = 0x06054b50, SIG_CEN = 0x02014b50, SIG_LOC = 0x04034b50;

/** True when the buffer starts with a ZIP local-file or empty-archive signature. */
export function isZip(bytes) {
  return bytes && bytes.length >= 4 && bytes[0] === 0x50 && bytes[1] === 0x4b && ((bytes[2] === 3 && bytes[3] === 4) || (bytes[2] === 5 && bytes[3] === 6));
}

/**
 * List the entries of a ZIP archive from its central directory.
 * @param {Uint8Array} bytes
 * @returns {{ name: string, method: number, compressedSize: number, size: number, offset: number, crc: number, isDirectory: boolean, encrypted: boolean }[]}
 * @throws {Error} when the end-of-central-directory record is missing
 */
export function zipEntries(bytes) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let eocd = -1;
  const lowest = Math.max(0, bytes.length - 22 - 65535);
  for (let i = bytes.length - 22; i >= lowest; i--) if (dv.getUint32(i, true) === SIG_EOCD) { eocd = i; break; }
  if (eocd < 0) throw new Error('Not a ZIP archive: expected an end-of-central-directory record (PK\\x05\\x06) near the end of the file.');
  const total = dv.getUint16(eocd + 10, true);
  const cenSize = dv.getUint32(eocd + 12, true);
  const cenOff = dv.getUint32(eocd + 16, true);
  if (total === 0xffff || cenOff === 0xffffffff) throw new Error('ZIP64 archives are not supported: re-zip the files as a standard ZIP below 4 GiB, or upload the files unzipped.');
  // Tolerate archives with bytes prepended (self-extractors): shift when the directory is not where it says.
  let shift = 0;
  if (cenOff + cenSize !== eocd && eocd - cenSize >= 0 && dv.getUint32(eocd - cenSize, true) === SIG_CEN) shift = eocd - cenSize - cenOff;
  const out = [];
  const dec = new TextDecoder('utf-8');
  let p = cenOff + shift;
  for (let i = 0; i < total; i++) {
    if (p + 46 > bytes.length || dv.getUint32(p, true) !== SIG_CEN) throw new Error(`Corrupt ZIP central directory at entry ${i + 1}: expected signature PK\\x01\\x02.`);
    const flags = dv.getUint16(p + 8, true);
    const method = dv.getUint16(p + 10, true);
    const crc = dv.getUint32(p + 16, true);
    const compressedSize = dv.getUint32(p + 20, true);
    const size = dv.getUint32(p + 24, true);
    const nl = dv.getUint16(p + 28, true), xl = dv.getUint16(p + 30, true), cl = dv.getUint16(p + 32, true);
    const offset = dv.getUint32(p + 42, true) + shift;
    const name = dec.decode(bytes.subarray(p + 46, p + 46 + nl));
    out.push({ name, method, compressedSize, size, offset, crc, isDirectory: name.endsWith('/'), encrypted: (flags & 1) === 1 });
    p += 46 + nl + xl + cl;
  }
  return out;
}

/** Decompress one entry returned by zipEntries(). */
export async function readZipEntry(bytes, entry) {
  if (entry.encrypted) throw new Error(`ZIP entry "${entry.name}" is encrypted; extract it and upload the plain file.`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const p = entry.offset;
  if (p + 30 > bytes.length || dv.getUint32(p, true) !== SIG_LOC) throw new Error(`Corrupt ZIP: no local header (PK\\x03\\x04) for "${entry.name}".`);
  const start = p + 30 + dv.getUint16(p + 26, true) + dv.getUint16(p + 28, true);
  const data = bytes.subarray(start, start + entry.compressedSize);
  if (entry.method === 0) return data.slice();
  if (entry.method === 8) return inflateRaw(data, entry.size);
  throw new Error(`ZIP entry "${entry.name}" uses compression method ${entry.method}; only stored (0) and deflate (8) are supported.`);
}

/**
 * Read all (or selected) files of a ZIP archive.
 * @param {Uint8Array} bytes
 * @param {object} [opts]
 * @param {(name: string) => boolean} [opts.filter] keep only entries for which this returns true
 * @returns {Promise<Map<string, Uint8Array>>} file name (as stored, forward slashes) → content
 */
export async function readZip(bytes, opts = {}) {
  const out = new Map();
  for (const e of zipEntries(bytes)) {
    if (e.isDirectory) continue;
    if (opts.filter && !opts.filter(e.name)) continue;
    out.set(e.name, await readZipEntry(bytes, e));
  }
  return out;
}

async function streamInflate(bytes, format) {
  const ds = new DecompressionStream(format);
  const writer = ds.writable.getWriter();
  writer.write(bytes).catch(() => {});
  writer.close().catch(() => {});
  const reader = ds.readable.getReader();
  const chunks = [];
  let total = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    total += value.length;
  }
  if (chunks.length === 1) return chunks[0];
  const out = new Uint8Array(total);
  let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

/** Inflate a raw deflate stream (RFC 1951). `sizeHint` is the expected output length, if known. */
export async function inflateRaw(bytes, sizeHint) {
  if (typeof DecompressionStream === 'function') {
    try { return await streamInflate(bytes, 'deflate-raw'); } catch { /* fall through for a precise error */ }
  }
  return inflateSync(bytes, sizeHint);
}

/** Inflate a zlib stream (RFC 1950: 2-byte header, deflate data, Adler-32). */
export async function inflateZlib(bytes, sizeHint) {
  if (typeof DecompressionStream === 'function') {
    try { return await streamInflate(bytes, 'deflate'); } catch { /* fall through */ }
  }
  return inflateSync(zlibBody(bytes), sizeHint);
}

/** Decompress a gzip member (RFC 1952). */
export async function gunzip(bytes, sizeHint) {
  if (typeof DecompressionStream === 'function') {
    try { return await streamInflate(bytes, 'gzip'); } catch { /* fall through */ }
  }
  if (bytes[0] !== 0x1f || bytes[1] !== 0x8b) throw new Error('Not gzip data: expected magic bytes 1F 8B.');
  const flg = bytes[3];
  let p = 10;
  if (flg & 4) p += 2 + (bytes[p] | (bytes[p + 1] << 8));
  if (flg & 8) { while (bytes[p++] !== 0 && p < bytes.length); }
  if (flg & 16) { while (bytes[p++] !== 0 && p < bytes.length); }
  if (flg & 2) p += 2;
  return inflateSync(bytes.subarray(p, bytes.length - 8), sizeHint);
}

function zlibBody(bytes) {
  if (bytes.length < 6 || (bytes[0] & 0x0f) !== 8 || ((bytes[0] << 8) | bytes[1]) % 31 !== 0) throw new Error('Not zlib data: expected a deflate stream with a valid 2-byte zlib header.');
  if (bytes[1] & 0x20) throw new Error('zlib preset dictionaries are not supported.');
  return bytes.subarray(2, bytes.length - 4);
}

/** Synchronous zlib inflate (pure JavaScript). Use where an async call is not possible. */
export function inflateZlibSync(bytes, sizeHint) {
  return inflateSync(zlibBody(bytes), sizeHint);
}

const LEN_BASE = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEN_EXTRA = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DIST_BASE = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DIST_EXTRA = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];
const CL_ORDER = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

function huffman(lengths, off, n) {
  const count = new Uint16Array(16), symbol = new Uint16Array(n), offs = new Uint16Array(16);
  for (let i = 0; i < n; i++) count[lengths[off + i]]++;
  count[0] = 0;
  for (let i = 1; i < 16; i++) offs[i] = offs[i - 1] + count[i - 1];
  for (let i = 0; i < n; i++) { const l = lengths[off + i]; if (l) symbol[offs[l]++] = i; }
  return { count, symbol };
}

let FIXED = null;
function fixedTables() {
  if (FIXED) return FIXED;
  const l = new Uint8Array(288 + 30);
  for (let i = 0; i < 144; i++) l[i] = 8;
  for (let i = 144; i < 256; i++) l[i] = 9;
  for (let i = 256; i < 280; i++) l[i] = 7;
  for (let i = 280; i < 288; i++) l[i] = 8;
  for (let i = 0; i < 30; i++) l[288 + i] = 5;
  FIXED = { lit: huffman(l, 0, 288), dist: huffman(l, 288, 30) };
  return FIXED;
}

/**
 * Synchronous raw-deflate decoder (RFC 1951), canonical-Huffman bit-serial decoding. Slower than the platform
 * stream but dependency-free and usable anywhere.
 * @param {Uint8Array} src raw deflate data
 * @param {number} [sizeHint] expected output size
 * @returns {Uint8Array}
 * @throws {Error} on a malformed stream
 */
export function inflateSync(src, sizeHint) {
  let out = new Uint8Array(sizeHint > 0 ? sizeHint : Math.max(1024, src.length * 4));
  let op = 0, pos = 0, bitBuf = 0, bitCnt = 0;
  const need = (extra) => {
    if (op + extra <= out.length) return;
    let cap = out.length * 2;
    while (cap < op + extra) cap *= 2;
    const b = new Uint8Array(cap);
    b.set(out.subarray(0, op));
    out = b;
  };
  const bits = (n) => {
    while (bitCnt < n) {
      if (pos >= src.length) throw new Error('Deflate stream ends unexpectedly (truncated data).');
      bitBuf |= src[pos++] << bitCnt;
      bitCnt += 8;
    }
    const v = bitBuf & ((1 << n) - 1);
    bitBuf >>>= n;
    bitCnt -= n;
    return v;
  };
  const decode = (h) => {
    let code = 0, first = 0, index = 0;
    for (let len = 1; len <= 15; len++) {
      code |= bits(1);
      const c = h.count[len];
      if (code - c < first) return h.symbol[index + (code - first)];
      index += c;
      first += c;
      first <<= 1;
      code <<= 1;
    }
    throw new Error('Invalid Huffman code in deflate stream.');
  };
  let last = 0;
  while (!last) {
    last = bits(1);
    const type = bits(2);
    if (type === 0) {
      bitBuf = 0; bitCnt = 0;
      if (pos + 4 > src.length) throw new Error('Deflate stream ends unexpectedly (stored block header).');
      const len = src[pos] | (src[pos + 1] << 8);
      const nlen = src[pos + 2] | (src[pos + 3] << 8);
      pos += 4;
      if ((len ^ 0xffff) !== nlen) throw new Error('Corrupt stored block in deflate stream.');
      if (pos + len > src.length) throw new Error('Deflate stream ends unexpectedly (stored block).');
      need(len);
      out.set(src.subarray(pos, pos + len), op);
      op += len; pos += len;
      continue;
    }
    if (type === 3) throw new Error('Invalid block type in deflate stream.');
    let lit, dist;
    if (type === 1) ({ lit, dist } = fixedTables());
    else {
      const nlen = bits(5) + 257, ndist = bits(5) + 1, ncode = bits(4) + 4;
      const lengths = new Uint8Array(320);
      for (let i = 0; i < ncode; i++) lengths[CL_ORDER[i]] = bits(3);
      const cl = huffman(lengths, 0, 19);
      const ll = new Uint8Array(nlen + ndist);
      let i = 0;
      while (i < nlen + ndist) {
        const sym = decode(cl);
        if (sym < 16) ll[i++] = sym;
        else {
          let rep, val = 0;
          if (sym === 16) { if (!i) throw new Error('Corrupt code lengths in deflate stream.'); val = ll[i - 1]; rep = 3 + bits(2); } else if (sym === 17) rep = 3 + bits(3); else rep = 11 + bits(7);
          if (i + rep > nlen + ndist) throw new Error('Corrupt code lengths in deflate stream.');
          while (rep--) ll[i++] = val;
        }
      }
      lit = huffman(ll, 0, nlen);
      dist = huffman(ll, nlen, ndist);
    }
    for (;;) {
      const sym = decode(lit);
      if (sym < 256) { if (op === out.length) need(1); out[op++] = sym; continue; }
      if (sym === 256) break;
      const li = sym - 257;
      if (li >= 29) throw new Error('Invalid length symbol in deflate stream.');
      const len = LEN_BASE[li] + bits(LEN_EXTRA[li]);
      const di = decode(dist);
      if (di >= 30) throw new Error('Invalid distance symbol in deflate stream.');
      const d = DIST_BASE[di] + bits(DIST_EXTRA[di]);
      if (d > op) throw new Error('Deflate back-reference before start of output.');
      need(len);
      for (let k = 0; k < len; k++, op++) out[op] = out[op - d];
    }
  }
  return op === out.length ? out : out.slice(0, op);
}

let CRC_TABLE = null;
/** CRC-32 (IEEE 802.3) of a byte array. */
export function crc32(bytes) {
  if (!CRC_TABLE) {
    CRC_TABLE = new Uint32Array(256);
    for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; CRC_TABLE[n] = c >>> 0; }
  }
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

/**
 * Write a ZIP archive with stored (uncompressed) entries — enough for packaging exports and fixtures.
 * @param {{ name: string, data: Uint8Array }[]} files
 * @returns {Uint8Array}
 */
export function makeZip(files) {
  const enc = new TextEncoder();
  const recs = files.map((f) => ({ name: enc.encode(f.name), data: f.data, crc: crc32(f.data) }));
  let size = 22;
  for (const r of recs) size += 30 + r.name.length + r.data.length + 46 + r.name.length;
  const out = new Uint8Array(size);
  const dv = new DataView(out.buffer);
  let p = 0;
  for (const r of recs) {
    r.offset = p;
    dv.setUint32(p, SIG_LOC, true); dv.setUint16(p + 4, 20, true); dv.setUint16(p + 6, 0x0800, true);
    dv.setUint32(p + 14, r.crc, true); dv.setUint32(p + 18, r.data.length, true); dv.setUint32(p + 22, r.data.length, true);
    dv.setUint16(p + 26, r.name.length, true);
    out.set(r.name, p + 30);
    out.set(r.data, p + 30 + r.name.length);
    p += 30 + r.name.length + r.data.length;
  }
  const cen = p;
  for (const r of recs) {
    dv.setUint32(p, SIG_CEN, true); dv.setUint16(p + 4, 20, true); dv.setUint16(p + 6, 20, true); dv.setUint16(p + 8, 0x0800, true);
    dv.setUint32(p + 16, r.crc, true); dv.setUint32(p + 20, r.data.length, true); dv.setUint32(p + 24, r.data.length, true);
    dv.setUint16(p + 28, r.name.length, true);
    dv.setUint32(p + 42, r.offset, true);
    out.set(r.name, p + 46);
    p += 46 + r.name.length;
  }
  dv.setUint32(p, SIG_EOCD, true);
  dv.setUint16(p + 8, recs.length, true); dv.setUint16(p + 10, recs.length, true);
  dv.setUint32(p + 12, p - cen, true); dv.setUint32(p + 16, cen, true);
  return out;
}
