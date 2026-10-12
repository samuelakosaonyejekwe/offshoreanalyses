// GeoTIFF elevation / bathymetry rasters (classic TIFF and BigTIFF, either byte order).
//
// Read: stripped and tiled layouts; compression none (1), LZW (5), Deflate (8 and the legacy 32946) and
// PackBits (32773); horizontal (2) and floating-point (3) predictors; 8/16/32-bit integer and 32/64-bit float
// samples; chunky or planar samples (the first band is used); GDAL_NODATA; georeferencing from
// ModelTiepoint + ModelPixelScale or an unrotated ModelTransformation; the CRS from the GeoKey directory.
// Large rasters are block-averaged while decoding so that the full-resolution image is never held in memory:
// the result has at most `opts.maxCells` cells (default 4 000 000).
//
// Not read: JPEG, LERC, ZSTD, WebP, CCITT compressions; rotated or sheared transforms; palette and 1-bit
// images; internal overviews and masks (the first image directory is used).

import { inputBytes, baseName } from './common.js';
import { geoKeysToCrs } from '../units.js';
import { inflateZlib } from '../zip.js';
import { finishGrid, decimationFactor, DEFAULT_MAX_CELLS } from './grids.js';

const TYPE_SIZE = [0, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8, 0, 0, 0, 8, 8, 8];

/** TIFF-flavoured LZW (MSB-first codes, 9–12 bits, early change). */
export function lzwDecode(src, expected) {
  let out = new Uint8Array(Math.max(expected, 16));
  const prefix = new Uint16Array(4097), suffix = new Uint8Array(4097), len = new Uint16Array(4097);
  for (let i = 0; i < 256; i++) len[i] = 1;
  let op = 0, pos = 0, bitBuf = 0, bitCnt = 0, width = 9, next = 258, old = -1;
  const emit = (code) => {
    const n = len[code];
    if (op + n > out.length) { const b = new Uint8Array(Math.max(out.length * 2, op + n)); b.set(out.subarray(0, op)); out = b; }
    let p = op + n - 1, c = code;
    while (c >= 258) { out[p--] = suffix[c]; c = prefix[c]; }
    out[p] = c;
    op += n;
  };
  for (;;) {
    while (bitCnt < width) { if (pos >= src.length) return out.subarray(0, op); bitBuf = ((bitBuf << 8) | src[pos++]) >>> 0; bitCnt += 8; }
    const code = (bitBuf >>> (bitCnt - width)) & ((1 << width) - 1);
    bitCnt -= width;
    bitBuf &= (1 << bitCnt) - 1;
    if (code === 257) break;
    if (code === 256) { width = 9; next = 258; old = -1; continue; }
    if (old < 0) { if (code > 255) throw new Error('Corrupt LZW stream: the first code after a clear is not a literal.'); emit(code); old = code; continue; }
    if (code > next || next > 4096) throw new Error('Corrupt LZW stream: code outside the string table.');
    const start = op;
    if (code < next) {
      emit(code);
      if (next < 4096) { prefix[next] = old; suffix[next] = out[start]; len[next] = len[old] + 1; }
    } else {
      // the code being defined right now: previous string + its own first byte
      const first = out[op - len[old]];
      prefix[next] = old; suffix[next] = first; len[next] = len[old] + 1;
      emit(next);
    }
    next++;
    if (next === (1 << width) - 1 && width < 12) width++;
    old = code;
  }
  return out.subarray(0, op);
}

/** Apple PackBits run-length decoding. */
export function packBitsDecode(src, expected) {
  const out = new Uint8Array(expected);
  let op = 0, p = 0;
  while (p < src.length && op < expected) {
    const n = src[p++];
    if (n < 128) { const c = Math.min(n + 1, expected - op, src.length - p); out.set(src.subarray(p, p + c), op); op += c; p += n + 1; } else if (n > 128) { const c = Math.min(257 - n, expected - op); out.fill(src[p++], op, op + c); op += c; }
  }
  return out;
}

function readIfd(bytes, dv, le, big, off) {
  const tags = new Map();
  const u64 = (o) => Number(dv.getBigUint64(o, le));
  const n = big ? u64(off) : dv.getUint16(off, le);
  const esz = big ? 20 : 12;
  let p = off + (big ? 8 : 2);
  for (let i = 0; i < n; i++, p += esz) {
    if (p + esz > bytes.length) break;
    const tag = dv.getUint16(p, le), type = dv.getUint16(p + 2, le);
    const count = big ? u64(p + 4) : dv.getUint32(p + 4, le);
    const size = (TYPE_SIZE[type] || 0) * count;
    const inline = big ? 8 : 4;
    const vo = size <= inline ? p + (big ? 12 : 8) : (big ? u64(p + 12) : dv.getUint32(p + 8, le));
    tags.set(tag, { type, count, off: vo });
  }
  return tags;
}

function tagValues(bytes, dv, le, t) {
  if (!t) return null;
  const { type, count, off } = t;
  if (off + (TYPE_SIZE[type] || 0) * count > bytes.length) return null;
  if (type === 2) { let s = ''; for (let i = 0; i < count; i++) { const c = bytes[off + i]; if (c === 0) break; s += String.fromCharCode(c); } return s; }
  const out = new Array(count);
  for (let i = 0; i < count; i++) {
    const o = off + i * TYPE_SIZE[type];
    switch (type) {
      case 1: case 7: out[i] = bytes[o]; break;
      case 6: out[i] = dv.getInt8(o); break;
      case 3: out[i] = dv.getUint16(o, le); break;
      case 8: out[i] = dv.getInt16(o, le); break;
      case 4: out[i] = dv.getUint32(o, le); break;
      case 9: out[i] = dv.getInt32(o, le); break;
      case 5: out[i] = dv.getUint32(o, le) / dv.getUint32(o + 4, le); break;
      case 10: out[i] = dv.getInt32(o, le) / dv.getInt32(o + 4, le); break;
      case 11: out[i] = dv.getFloat32(o, le); break;
      case 12: out[i] = dv.getFloat64(o, le); break;
      case 16: case 18: out[i] = Number(dv.getBigUint64(o, le)); break;
      case 17: out[i] = Number(dv.getBigInt64(o, le)); break;
      default: out[i] = 0;
    }
  }
  return out;
}

/** GeoTIFF raster → grid. */
export async function parseGeoTiff(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  const fname = input.name || 'TIFF';
  if (bytes.length < 16 || !((bytes[0] === 0x49 && bytes[1] === 0x49) || (bytes[0] === 0x4d && bytes[1] === 0x4d))) throw new Error(`${fname}: not a TIFF file (expected the byte-order mark II or MM at the start).`);
  const le = bytes[0] === 0x49;
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const magic = dv.getUint16(2, le);
  if (magic !== 42 && magic !== 43) throw new Error(`${fname}: not a TIFF file (expected version 42 or 43 after the byte-order mark, found ${magic}).`);
  const big = magic === 43;
  const ifdOff = big ? Number(dv.getBigUint64(8, le)) : dv.getUint32(4, le);
  if (!(ifdOff > 0 && ifdOff < bytes.length)) throw new Error(`${fname}: the first image directory lies outside the file (truncated TIFF).`);
  const tags = readIfd(bytes, dv, le, big, ifdOff);
  const val = (id) => tagValues(bytes, dv, le, tags.get(id));
  const one = (id, d) => { const v = val(id); return v && v.length ? v[0] : d; };
  const width = one(256, 0), height = one(257, 0);
  if (!(width > 0 && height > 0)) throw new Error(`${fname}: the TIFF has no ImageWidth / ImageLength tags.`);
  const bits = one(258, 1), comp = one(259, 1), spp = one(277, 1), planar = one(284, 1), predictor = one(317, 1), fmt = one(339, 1), photometric = one(262, 1);
  if (![8, 16, 32, 64].includes(bits) || (fmt === 3 && bits < 32) || photometric === 3) throw new Error(`${fname}: ${bits}-bit ${photometric === 3 ? 'palette' : fmt === 3 ? 'float' : 'integer'} samples are not an elevation raster this reader handles (expected 8/16/32-bit integers or 32/64-bit floats).`);
  if (![1, 5, 8, 32946, 32773].includes(comp)) {
    const nm = { 7: 'JPEG', 6: 'old JPEG', 34887: 'LERC', 50000: 'ZSTD', 50001: 'WebP', 34925: 'LZMA', 2: 'CCITT', 3: 'CCITT', 4: 'CCITT', 34712: 'JPEG 2000' }[comp] || `code ${comp}`;
    throw new Error(`${fname}: TIFF compression ${nm} is not decoded here. Re-write the file with LZW or Deflate, for example: gdal_translate -co COMPRESS=DEFLATE in.tif out.tif`);
  }
  const tiled = tags.has(322);
  const bw = tiled ? one(322) : width, bh = tiled ? one(323) : Math.min(one(278, height), height);
  const offsets = val(tiled ? 324 : 273), counts = val(tiled ? 325 : 279);
  if (!offsets || !counts || !offsets.length) throw new Error(`${fname}: the TIFF has no ${tiled ? 'tile' : 'strip'} offset table.`);
  const bps = bits / 8;
  const chunky = planar !== 2;
  const pix = chunky ? spp * bps : bps;
  const across = Math.ceil(width / bw), down = Math.ceil(height / bh);
  const nodataTag = val(42113);
  const nodata = nodataTag !== null && String(nodataTag).trim() !== '' ? Number(String(nodataTag).trim()) : null;

  // georeferencing
  const scale = val(33550), tie = val(33922), xf = val(34264);
  const warnings = [];
  let sx = 1, sy = 1, X = 0, Y = height, I = 0, J = 0, georef = 'none';
  if (scale && tie && tie.length >= 6) { sx = scale[0]; sy = scale[1]; I = tie[0]; J = tie[1]; X = tie[3]; Y = tie[4]; georef = 'tiepoint + pixel scale'; } else if (xf && xf.length >= 16) {
    if (Math.abs(xf[1]) > 1e-12 * Math.abs(xf[0]) || Math.abs(xf[4]) > 1e-12 * Math.abs(xf[5])) throw new Error(`${fname}: the raster is rotated or sheared (ModelTransformation with off-diagonal terms), which a north-up grid cannot represent. Warp it first, for example: gdalwarp -of GTiff in.tif out.tif`);
    sx = xf[0]; sy = -xf[5]; X = xf[3]; Y = xf[7]; georef = 'model transformation';
  } else warnings.push('The TIFF carries no georeferencing tags; pixel coordinates are used (1 unit per pixel, origin at the lower-left corner).');
  if (!(sx > 0) || sy === 0) throw new Error(`${fname}: invalid pixel scale (${sx}, ${sy}).`);

  // GeoKeys
  const dir = val(34735), dbl = val(34736) || [], asc = val(34737) || '';
  const keys = {};
  if (dir && dir.length >= 4) for (let k = 0; k < dir[3]; k++) {
    const o = 4 + 4 * k, id = dir[o], loc = dir[o + 1], cnt = dir[o + 2], vo = dir[o + 3];
    if (loc === 0) keys[id] = vo; else if (loc === 34736) keys[id] = dbl[vo]; else if (loc === 34737) keys[id] = String(asc).substr(vo, cnt).replace(/\|$/, '');
  }
  const crs = geoKeysToCrs(keys);
  const point = crs.rasterType === 'point';

  // decode, block-averaging into the output lattice as blocks arrive
  const f = decimationFactor(width, height, opts.maxCells ?? DEFAULT_MAX_CELLS);
  const mx = Math.ceil(width / f), my = Math.ceil(height / f);
  const sum = new Float64Array(mx * my), cnt = new Uint32Array(mx * my);
  const getter = fmt === 3 ? (bits === 32 ? 'getFloat32' : 'getFloat64') : fmt === 2 ? ({ 8: 'getInt8', 16: 'getInt16', 32: 'getInt32' })[bits] : ({ 8: 'getUint8', 16: 'getUint16', 32: 'getUint32' })[bits];
  if (!getter) throw new Error(`${fname}: ${bits}-bit ${fmt === 2 ? 'signed' : 'unsigned'} integer samples are not supported.`);
  const nBlocks = across * down;
  for (let b = 0; b < nBlocks; b++) {
    if (b >= offsets.length) break;
    const bx = (b % across) * bw, by = Math.floor(b / across) * bh;
    const rows = tiled ? bh : Math.min(bh, height - by);
    const expected = bw * rows * pix;
    if (!counts[b]) continue; // sparse file: block not written
    if (offsets[b] + counts[b] > bytes.length) throw new Error(`${fname}: ${tiled ? 'tile' : 'strip'} ${b + 1} of ${nBlocks} lies outside the file (truncated TIFF).`);
    const raw = bytes.subarray(offsets[b], offsets[b] + counts[b]);
    let data;
    if (comp === 1) data = raw;
    else if (comp === 5) data = lzwDecode(raw, expected);
    else if (comp === 32773) data = packBitsDecode(raw, expected);
    else data = await inflateZlib(raw.slice(), expected);
    if (data.length < expected) { const pad = new Uint8Array(expected); pad.set(data); data = pad; warnings.push(`${tiled ? 'Tile' : 'Strip'} ${b + 1} decoded short; missing pixels are zero-filled.`); }
    if (data === raw && predictor !== 1) data = raw.slice();
    const rowBytes = bw * pix, ns = chunky ? spp : 1;
    if (predictor === 2) {
      const d = new DataView(data.buffer, data.byteOffset, data.byteLength);
      for (let r = 0; r < rows; r++) {
        const o = r * rowBytes;
        if (bps === 1) for (let i = ns; i < rowBytes; i++) data[o + i] = (data[o + i] + data[o + i - ns]) & 0xff;
        else if (bps === 2) for (let i = ns; i < bw * ns; i++) d.setUint16(o + 2 * i, (d.getUint16(o + 2 * i, le) + d.getUint16(o + 2 * (i - ns), le)) & 0xffff, le);
        else if (bps === 4) for (let i = ns; i < bw * ns; i++) d.setUint32(o + 4 * i, (d.getUint32(o + 4 * i, le) + d.getUint32(o + 4 * (i - ns), le)) >>> 0, le);
      }
    } else if (predictor === 3) {
      // floating-point predictor: byte-wise differences over a row whose bytes are grouped by significance
      const tmp = new Uint8Array(rowBytes), nv = bw * ns;
      for (let r = 0; r < rows; r++) {
        const o = r * rowBytes;
        for (let i = ns; i < rowBytes; i++) data[o + i] = (data[o + i] + data[o + i - ns]) & 0xff;
        for (let v = 0; v < nv; v++) for (let k = 0; k < bps; k++) tmp[v * bps + (le ? bps - 1 - k : k)] = data[o + k * nv + v];
        data.set(tmp, o);
      }
    } else if (predictor !== 1) throw new Error(`${fname}: TIFF predictor ${predictor} is not supported.`);
    const d = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const cols = Math.min(bw, width - bx);
    for (let r = 0; r < rows && by + r < height; r++) {
      const orow = Math.floor((by + r) / f) * mx;
      let o = r * rowBytes;
      for (let c = 0; c < cols; c++, o += pix) {
        const v = bps === 1 ? (fmt === 2 ? d.getInt8(o) : data[o]) : d[getter](o, le);
        if (v !== v || v === nodata) continue;
        const k = orow + Math.floor((bx + c) / f);
        sum[k] += v; cnt[k]++;
      }
    }
  }
  // image rows run north → south (for a positive y scale); the grid stores south → north
  const z = new Float32Array(mx * my);
  const flip = sy > 0;
  for (let r = 0; r < my; r++) for (let c = 0; c < mx; c++) { const k = r * mx + c; z[(flip ? my - 1 - r : r) * mx + c] = cnt[k] ? sum[k] / cnt[k] : NaN; }
  const edgeX = X - I * sx - (point ? sx / 2 : 0), edgeY = Y + J * sy + (point ? sy / 2 : 0);
  const dx = f * sx, dyAbs = f * Math.abs(sy);
  const x0 = edgeX + dx / 2;
  const y0 = flip ? edgeY - (my - 0.5) * dyAbs : edgeY + dyAbs / 2;
  if (spp > 1) warnings.push(`The TIFF has ${spp} samples per pixel; the first band is used as elevation.`);
  return finishGrid({ name, nx: mx, ny: my, x0, y0, dx, dy: dyAbs, z }, { ...opts, maxCells: 0 }, {
    format: 'geotiff', crs: crs.kind === 'unknown' ? null : crs, zUnit: crs.verticalUnit !== 'm' ? crs.verticalUnit : null, warnings, decimated: f, sourceSize: [width, height],
    extra: {
      tiff: { bigTiff: big, byteOrder: le ? 'little' : 'big', layout: tiled ? `tiles ${bw} × ${bh}` : `strips of ${bh} rows`, compression: { 1: 'none', 5: 'LZW', 8: 'Deflate', 32946: 'Deflate', 32773: 'PackBits' }[comp], predictor, bitsPerSample: bits, sampleFormat: { 1: 'uint', 2: 'int', 3: 'float' }[fmt] || 'uint', samplesPerPixel: spp },
      georeferencing: georef, registration: point ? 'pixel is point' : 'pixel is area (cell centres)', nodata, geoKeys: keys, verticalEpsg: crs.verticalEpsg,
    },
  });
}
