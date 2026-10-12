// Compressed LiDAR / multibeam point clouds (LAZ, i.e. LAS with LASzip-compressed point records) through the
// vendored laz-perf decoder (see docs/vendor.md).
//
// A LAZ file is a LAS file whose public header and variable-length records are stored as usual; only the point
// records are compressed, in chunks, and a "laszip encoded" VLR (record 22204) describes how. This module
//   1. reads the header and that VLR and refuses, with the reason, what the decoder cannot open;
//   2. lets the decoder walk the records in file order (the arithmetic coder cannot skip) and keeps every
//      n-th one so that at most `opts.maxPoints` (default 2 000 000) are kept;
//   3. hands the kept records, as an uncompressed LAS image with the original header and VLRs, to the LAS
//      importer (importers/pointcloud.js). Scale and offset, intensity, classification, colour, the coordinate
//      reference system, unit handling and the lon/lat projection are therefore exactly those of a LAS file.
//
// Memory: the decoder needs the whole compressed file inside its WebAssembly heap, so a file occupies about
// twice its size while it is read, plus the kept records. Files of 2 GB and more cannot be opened that way.
//
// The runtime (JavaScript glue + WebAssembly) is 0.3 MB and is loaded on first use only, relative to this
// module.

import { inputBytes, decodeText } from '../importers/common.js';
import { parseLas } from '../importers/pointcloud.js';

export const id = 'laz';
export const title = 'Compressed LiDAR point clouds (LAZ, LASzip)';
export const formats = ['laz'];
export const extensions = ['laz'];
export const library = { name: 'laz-perf', version: '0.0.7', licence: 'Apache-2.0' };
/** Download size of the runtime (JavaScript glue + WebAssembly) [bytes]. */
export const sizeBytes = 89055 + 214351;
/** What to convert a file to so that it can be imported when this adapter cannot run. */
export const convertTo = 'uncompressed LAS (las2las -i in.laz -o out.las, or pdal translate in.laz out.las)';

const DEFAULT_MAX_POINTS = 2000000;
const BASE = new URL('../../../vendor/laz-perf/', import.meta.url);
const isNode = typeof process !== 'undefined' && !!(process.versions && process.versions.node) && typeof window === 'undefined';

/**
 * True when the decoder can run here: WebAssembly is supported and, under Node (where it can be checked
 * without a request), the vendored files are on disk.
 */
export function available() {
  if (typeof WebAssembly !== 'object' || typeof WebAssembly.instantiate !== 'function') return false;
  if (!isNode || typeof process.getBuiltinModule !== 'function') return true;
  try { const fs = process.getBuiltinModule('node:fs'); return ['laz-perf.js', 'laz-perf.wasm'].every((f) => fs.existsSync(new URL(f, BASE))); } catch { return true; }
}

let loading = null;

async function loadFactory() {
  const src = new URL('laz-perf.js', BASE);
  if (isNode) {
    const { createRequire } = await import('node:module');
    const { fileURLToPath } = await import('node:url');
    return createRequire(import.meta.url)(fileURLToPath(src));
  }
  const g = globalThis;
  if (typeof g.createLazPerf === 'function') return g.createLazPerf;
  if (typeof g.document !== 'undefined') {
    await new Promise((resolve, reject) => {
      const s = g.document.createElement('script');
      s.src = src.href;
      s.onload = resolve;
      s.onerror = () => reject(new Error('the decoder script could not be loaded'));
      g.document.head.append(s);
    });
  } else if (typeof g.importScripts === 'function') {
    try { g.importScripts(src.href); } catch { /* module workers have no importScripts: fall through */ }
  }
  if (typeof g.createLazPerf !== 'function') {
    // module worker: evaluate the classic script text in function scope
    const text = await (await fetch(src.href)).text();
    g.createLazPerf = new Function(`${text}\nreturn createLazPerf;`)();
  }
  return g.createLazPerf;
}

/**
 * Load and initialise the decoder (once; later calls return the same instance).
 * @returns {Promise<object>} the laz-perf module (LASZip, ChunkDecoder, _malloc, _free, HEAPU8)
 */
export function load() {
  if (!loading) {
    loading = (async () => {
      const factory = await loadFactory();
      const overrides = { locateFile: (f) => new URL(f, BASE).href };
      if (isNode) {
        const { readFile } = await import('node:fs/promises');
        overrides.wasmBinary = await readFile(new URL('laz-perf.wasm', BASE));
      }
      return factory(overrides);
    })().catch((err) => { loading = null; throw new Error(`The LAZ decoder could not be started: ${err && err.message ? err.message : err}. LAZ import needs the files in vendor/laz-perf/ and WebAssembly support; otherwise convert the file to ${convertTo}.`); });
  }
  return loading;
}

const cstr = (bytes, a, b) => decodeText(bytes.subarray(a, b)).replace(/\u0000[\s\S]*$/, '').trim();

/** True when the content is a LAS file whose point records are LASzip-compressed. */
export function sniff(input) {
  const b = inputBytes(input);
  return b.length >= 227 && b[0] === 0x4c && b[1] === 0x41 && b[2] === 0x53 && b[3] === 0x46 && (b[104] & 0xc0) !== 0;
}

const COMPRESSORS = { 0: 'none', 1: 'pointwise', 2: 'pointwise, chunked', 3: 'layered, chunked' };
const ITEM_TYPES = { 0: 'BYTE', 1: 'SHORT', 2: 'INT', 3: 'LONG', 4: 'FLOAT', 5: 'DOUBLE', 6: 'POINT10', 7: 'GPSTIME11', 8: 'RGB12', 9: 'WAVEPACKET13', 10: 'POINT14', 11: 'RGB14', 12: 'RGBNIR14', 13: 'WAVEPACKET14', 14: 'BYTE14' };
/** Shortest record of each LAS point data format [bytes]. */
const MIN_RECORD = [20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67];

/**
 * Read the public header of a LAS / LAZ file and its "laszip encoded" variable-length record.
 * Nothing is decompressed, so this works without the runtime.
 * @param {Uint8Array} bytes the file
 * @param {string} [fileName] used in error messages
 * @returns {{ version: string, versionMinor: number, headerSize: number, dataOffset: number, vlrCount: number,
 *   pointFormat: number, compressed: boolean, recordLength: number, pointCount: number, scale: number[],
 *   offset: number[], min: number[], max: number[], systemId: string, software: string,
 *   vlrs: { userId: string, recordId: number, offset: number, length: number, description: string }[],
 *   laszip: null | { compressor: number, compressorName: string, coder: number, version: string,
 *     options: number, chunkSize: number, variableChunks: boolean, items: { type: number, name: string,
 *     size: number, version: number }[], recordLength: number } }}
 *   `vlrs[i].offset` is where the record's payload starts; `chunkSize` is in points
 */
export function readLazHeader(bytes, fileName = 'LAZ') {
  if (bytes.length < 4 || cstr(bytes, 0, 4) !== 'LASF') throw new Error(`${fileName}: not a LAS / LAZ file (the file does not start with the signature "LASF").`);
  if (bytes.length < 227) throw new Error(`${fileName}: the file is truncated (${bytes.length} bytes; the LAS public header alone is 227 bytes).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const vMajor = bytes[24], vMinor = bytes[25];
  const headerSize = dv.getUint16(94, true), dataOffset = dv.getUint32(96, true), vlrCount = dv.getUint32(100, true);
  if (headerSize < 227 || headerSize > bytes.length) throw new Error(`${fileName}: the header size in the file (${headerSize}) is not possible for a ${bytes.length}-byte LAS file — the file is truncated or corrupt.`);
  if (dataOffset < headerSize || dataOffset > bytes.length) throw new Error(`${fileName}: the point data offset in the header (${dataOffset}) lies outside the file (${bytes.length} bytes) — the file is truncated or corrupt.`);
  const fmtByte = bytes[104];
  let pointCount = dv.getUint32(107, true);
  if (vMajor === 1 && vMinor >= 4 && headerSize >= 375) { const n64 = Number(dv.getBigUint64(247, true)); if (n64 > 0) pointCount = n64; }
  const vlrs = [];
  let laszip = null;
  let p = headerSize;
  for (let i = 0; i < vlrCount && p + 54 <= dataOffset; i++) {
    const userId = cstr(bytes, p + 2, p + 18), recordId = dv.getUint16(p + 18, true), length = dv.getUint16(p + 20, true);
    const body = p + 54;
    if (body + length > bytes.length) break;
    vlrs.push({ userId, recordId, offset: body, length, description: cstr(bytes, p + 22, p + 54) });
    if (recordId === 22204 && /^laszip encoded/i.test(userId) && length >= 34) {
      const nItems = dv.getUint16(body + 32, true);
      const items = [];
      let recordLength = 0;
      for (let k = 0; k < nItems && body + 34 + 6 * k + 6 <= body + length; k++) {
        const type = dv.getUint16(body + 34 + 6 * k, true), size = dv.getUint16(body + 36 + 6 * k, true);
        items.push({ type, name: ITEM_TYPES[type] || `type ${type}`, size, version: dv.getUint16(body + 38 + 6 * k, true) });
        recordLength += size;
      }
      const compressor = dv.getUint16(body, true), chunkSize = dv.getUint32(body + 12, true);
      laszip = {
        compressor, compressorName: COMPRESSORS[compressor] || `unknown (${compressor})`, coder: dv.getUint16(body + 2, true),
        version: `${bytes[body + 4]}.${bytes[body + 5]}.${dv.getUint16(body + 6, true)}`, options: dv.getUint32(body + 8, true),
        chunkSize, variableChunks: chunkSize === 0xffffffff, items, recordLength,
      };
    }
    p = body + length;
  }
  return {
    version: `${vMajor}.${vMinor}`, versionMinor: vMinor, headerSize, dataOffset, vlrCount,
    pointFormat: fmtByte & 0x3f, compressed: (fmtByte & 0xc0) !== 0, recordLength: dv.getUint16(105, true), pointCount,
    scale: [dv.getFloat64(131, true), dv.getFloat64(139, true), dv.getFloat64(147, true)],
    offset: [dv.getFloat64(155, true), dv.getFloat64(163, true), dv.getFloat64(171, true)],
    min: [dv.getFloat64(187, true), dv.getFloat64(203, true), dv.getFloat64(219, true)],
    max: [dv.getFloat64(179, true), dv.getFloat64(195, true), dv.getFloat64(211, true)],
    systemId: cstr(bytes, 26, 58), software: cstr(bytes, 58, 90), vlrs, laszip,
  };
}

const ordinal = (n) => `${n}${n % 100 >= 11 && n % 100 <= 13 ? 'th' : n % 10 === 1 ? 'st' : n % 10 === 2 ? 'nd' : n % 10 === 3 ? 'rd' : 'th'}`;

/**
 * Import a LAZ file.
 * @param {{ name?: string, bytes?: Uint8Array }} input file content
 * @param {object} [opts] options of the LAS importer (`units`, `vertical`, `projection`, `origin`) and
 * @param {number} [opts.maxPoints] keep at most this many points by taking every n-th record
 *   (default 2 000 000; 0 keeps all)
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} one `points` model in metres with
 *   `attrs.classification`, `attrs.intensity` (when any is non-zero) and `attrs.rgb` (formats with colour);
 *   meta as for LAS plus `compression: { compressor, version, chunkSize, items }` and `decoder`
 */
export async function parse(input, opts = {}) {
  const fname = (input && input.name) || 'LAZ';
  const bytes = inputBytes(input || {});
  const h = readLazHeader(bytes, fname);
  // an uncompressed LAS file with a .laz name needs no decoder
  if (!h.compressed && !h.laszip) return parseLas({ name: fname, bytes }, { ...opts, maxPoints: opts.maxPoints ?? DEFAULT_MAX_POINTS });
  if (!h.laszip) throw new Error(`${fname}: the header marks the point records as compressed but the file has no "laszip encoded" record (VLR 22204) describing the compression, so they cannot be decoded. The file was probably truncated or rewritten by a tool that dropped the record; export it again or convert the original to ${convertTo}.`);
  const z = h.laszip;
  if (z.compressor === 0) throw new Error(`${fname}: the LASzip record declares no compression although the header says the points are compressed — the file is inconsistent. Convert it to ${convertTo}.`);
  if (z.compressor === 1) throw new Error(`${fname}: the file uses the original pointwise LASzip stream without chunks (written by LAStools before 2011), which this decoder does not read. Rewrite it with a current tool ("laszip -i in.laz -o out.laz") or convert it to ${convertTo}.`);
  if (z.compressor > 3) throw new Error(`${fname}: unknown LASzip compressor ${z.compressor} (expected 2 or 3). Convert the file to ${convertTo}.`);
  if (z.coder !== 0) throw new Error(`${fname}: unknown LASzip entropy coder ${z.coder} (only the arithmetic coder, 0, exists in the published format).`);
  if (h.pointFormat > 10) throw new Error(`${fname}: LAS point data format ${h.pointFormat} is not known (expected 0–10).`);
  if ([4, 5, 9, 10].includes(h.pointFormat) || z.items.some((it) => it.type === 9 || it.type === 13)) throw new Error(`${fname}: point data format ${h.pointFormat} carries full-waveform packets, which this decoder does not read. Convert the file to ${convertTo} (the LAS importer reads the point attributes of formats 4, 5, 9 and 10).`);
  if (h.pointFormat >= 6 && z.compressor !== 3) throw new Error(`${fname}: point data format ${h.pointFormat} must use the layered LASzip compressor (3), but the file declares ${z.compressor}. Convert it to ${convertTo}.`);
  if (h.recordLength < MIN_RECORD[h.pointFormat]) throw new Error(`${fname}: point record length ${h.recordLength} is shorter than format ${h.pointFormat} requires (${MIN_RECORD[h.pointFormat]} bytes).`);
  if (z.recordLength !== h.recordLength) throw new Error(`${fname}: the LASzip record describes points of ${z.recordLength} bytes but the header says ${h.recordLength} — the file is inconsistent and cannot be decoded safely.`);
  if (bytes.length >= 0x7fffffff) throw new Error(`${fname}: the file is ${(bytes.length / 2 ** 30).toFixed(1)} GB; the decoder needs the whole file in its 32-bit memory and cannot open files of 2 GB or more. Split it into tiles (lassplit / pdal tile) first.`);
  if (bytes.length < h.dataOffset + 8) throw new Error(`${fname}: the file ends before the compressed point data (truncated).`);

  const count = h.pointCount;
  const maxPts = opts.maxPoints ?? DEFAULT_MAX_POINTS;
  const stride = maxPts > 0 && count > maxPts ? Math.ceil(count / maxPts) : 1;
  const n = Math.ceil(count / stride);
  const recLen = h.recordLength;
  const las = new Uint8Array(h.dataOffset + n * recLen);
  las.set(bytes.subarray(0, h.dataOffset));
  const out = new DataView(las.buffer);
  las[104] = h.pointFormat;
  out.setUint32(107, h.pointFormat < 6 && n <= 0xffffffff ? n : 0, true);
  if (h.versionMinor >= 4 && h.headerSize >= 375) out.setBigUint64(247, BigInt(n), true);

  let decoded = 0;
  if (count > 0) {
    const mod = await load();
    let filePtr = 0, pointPtr = 0, zip = null, failed = false;
    try {
      filePtr = mod._malloc(bytes.length);
      pointPtr = mod._malloc(recLen);
      if (!filePtr || !pointPtr) throw new Error('out of memory');
      mod.HEAPU8.set(bytes, filePtr);
      zip = new mod.LASZip();
      zip.open(filePtr, bytes.length);
      const len = zip.getPointLength();
      if (len !== recLen) throw new Error(`the decoder reports ${len}-byte records, the header ${recLen}`);
      for (let i = 0, k = 0; i < count; i++) {
        zip.getPoint(pointPtr);
        // the heap may have been reallocated by the call: take the view afterwards
        if (i % stride === 0) { las.set(mod.HEAPU8.subarray(pointPtr, pointPtr + recLen), h.dataOffset + k * recLen); k++; }
        decoded++;
      }
    } catch (err) {
      // the decoder is built without C++ exception support: a failure inside it leaves the instance unusable,
      // so it is dropped and a fresh one is created for the next file
      failed = true;
      loading = null;
      const msg = err instanceof Error ? err.message : typeof err === 'string' ? err : '';
      const why = !msg || /exception catching is disabled|unreachable|out of bounds|aborted/i.test(msg) ? 'the compressed stream is damaged or ends early' : msg;
      throw new Error(`${fname}: the point records could not be decompressed after ${decoded} of ${count} points (${why}). The file is probably truncated or corrupt; try "lasinfo -i in.laz" to check it, or convert a good copy to ${convertTo}.`);
    } finally {
      if (!failed) {
        try { if (zip) zip.delete(); } catch { /* the instance is gone either way */ }
        if (pointPtr) mod._free(pointPtr);
        if (filePtr) mod._free(filePtr);
      }
    }
  }

  // The compressed stream has no checksum: damaged bytes decode into wrong coordinates without any error.
  // The bounding box of the header is the only independent record to compare with.
  let outside = 0;
  const boxOk = [0, 1, 2].every((k) => Number.isFinite(h.min[k]) && Number.isFinite(h.max[k]) && h.max[k] >= h.min[k]) && (h.max[0] > h.min[0] || h.max[1] > h.min[1]);
  if (boxOk) {
    const tol = [0, 1, 2].map((k) => 2 * Math.abs(h.scale[k]) + 1e-9 * Math.max(Math.abs(h.min[k]), Math.abs(h.max[k])));
    for (let k = 0; k < n; k++) {
      const o = h.dataOffset + k * recLen;
      for (let a = 0; a < 3; a++) { const q = out.getInt32(o + 4 * a, true) * h.scale[a] + h.offset[a]; if (q < h.min[a] - tol[a] || q > h.max[a] + tol[a]) { outside++; break; } }
    }
  }

  const r = parseLas({ name: fname, bytes: las }, { ...opts, maxPoints: 0 });
  r.meta.format = 'laz';
  r.meta.pointCount = count;
  r.meta.thinning = stride;
  r.meta.skipped = { ...(r.meta.skipped || {}), thinnedPoints: count - n };
  r.meta.compression = { compressor: z.compressorName, version: z.version, chunkSize: z.variableChunks ? 'variable' : z.chunkSize, items: z.items.map((it) => `${it.name} v${it.version}${it.type === 0 || it.type === 14 ? ` × ${it.size}` : ''}`) };
  r.meta.decoder = `${library.name} ${library.version}`;
  if (outside) r.warnings.push(`${outside} of the ${n} points read lie outside the bounding box stated in the file header. Either the header is wrong or the compressed data is damaged (LAZ has no checksum); compare with "lasinfo -i in.laz".`);
  if (stride > 1) r.warnings.unshift(`The cloud has ${count} points; every ${ordinal(stride)} point was kept (${n}) to stay under maxPoints.`);
  return r;
}
