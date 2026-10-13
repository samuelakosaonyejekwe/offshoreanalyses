// Point clouds: ASPRS LAS 1.0–1.4 (uncompressed, point formats 0–10), Leica PTS and PTX, and header-less
// XYZ + colour / intensity / normal text files.
//
// LAS: X/Y/Z are scaled and offset per the header; intensity, classification and RGB are kept as attributes.
// The CRS is taken from the GeoTIFF-key or WKT variable-length records when present. Clouds larger than
// `opts.maxPoints` (default 5 000 000) are thinned by keeping every n-th record. LAZ (compressed LAS) is
// not decoded here.
//
// PTX: each scan's 4 × 4 registration matrix is applied; invalid returns (0 0 0) are dropped.

import { inputBytes, lineCursor, baseName, Grow, finish, pointsModel, decodeText, M4, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits, geoKeysToCrs, parsePrj, crsLabel, geoProjector } from '../units.js';

const DEFAULT_MAX_POINTS = 5000000;

/** Common tail of the point-cloud readers (also used by e57.js): units, CRS, report. Scales `coords` in place. */
export function cloudResult(name, coords, attrs, opts, o) {
  const warnings = o.warnings || [];
  const n = coords.length / 3;
  if (!n) throw new Error(`${o.fileName || name}: no points were found (${o.expected}).`);
  let units, crsMeta = null;
  const crs = o.crs && o.crs.kind !== 'unknown' ? o.crs : null;
  if (crs && crs.kind === 'geographic') {
    let lo = Infinity, hi = -Infinity, la = Infinity, lb = -Infinity;
    for (let i = 0; i < coords.length; i += 3) { const x = coords[i], y = coords[i + 1]; if (x < lo) lo = x; if (x > hi) hi = x; if (y < la) la = y; if (y > lb) lb = y; }
    const P = geoProjector(opts, (lo + hi) / 2, (la + lb) / 2);
    for (let i = 0; i < coords.length; i += 3) { const p = P.forward(coords[i], coords[i + 1]); coords[i] = p[0]; coords[i + 1] = p[1]; }
    units = { unit: 'deg', scale: null, basis: 'file' };
    crsMeta = { source: crsLabel(crs), kind: 'geographic', epsg: crs.epsg ?? null, name: crs.name || null, basis: 'file', ...P.crs };
  } else {
    units = crs && crs.unitScale && crs.unitScale !== 1 && !opts.units ? { unit: crs.unit, scale: crs.unitScale, basis: 'file' } : resolveUnits(opts, o.unitHint || (crs ? 'm' : null), 'm');
    const s = units.scale, zs = s * (opts.vertical === 'depth' ? -1 : 1);
    if (s !== 1 || zs !== 1) for (let i = 0; i < coords.length; i += 3) { coords[i] *= s; coords[i + 1] *= s; coords[i + 2] *= zs; }
    if (crs) crsMeta = { source: crsLabel(crs), kind: crs.kind, epsg: crs.epsg ?? null, name: crs.name || null, basis: 'file', units: 'm' };
  }
  const r = finish([pointsModel(name, coords, attrs)], warnings, { format: o.format, units, crs: crsMeta, elementTypes: { point1: n }, skipped: o.skipped || {}, ...(o.extra || {}) });
  if (!crs) unitAdvice(units, r.meta.bbox, (b) => { const sg = suggestUnits(b); return sg && sg.unit === 'km' ? sg : null; }, r.warnings);
  return r;
}

/** ASPRS LAS (uncompressed). */
export function parseLas(input, opts = {}) {
  const bytes = inputBytes(input);
  const fname = input.name || 'LAS';
  if (bytes.length < 227 || decodeText(bytes.subarray(0, 4)) !== 'LASF') throw new Error(`${fname}: not a LAS file (expected the signature "LASF" and a public header of at least 227 bytes).`);
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const vMajor = bytes[24], vMinor = bytes[25];
  const headerSize = dv.getUint16(94, true), dataOff = dv.getUint32(96, true), nVlr = dv.getUint32(100, true);
  const fmtByte = bytes[104], fmt = fmtByte & 0x3f, recLen = dv.getUint16(105, true);
  if (fmtByte & 0xc0) throw new Error(`${fname}: the point records are LASzip-compressed (LAZ). Decompress first, for example "laszip -i in.laz -o out.las" or "pdal translate in.laz out.las".`);
  if (fmt > 10) throw new Error(`${fname}: LAS point data format ${fmt} is not known (expected 0–10).`);
  let count = dv.getUint32(107, true);
  if (vMajor === 1 && vMinor >= 4 && headerSize >= 375) { const n64 = Number(dv.getBigUint64(247, true)); if (n64 > 0) count = n64; }
  const sc = [dv.getFloat64(131, true), dv.getFloat64(139, true), dv.getFloat64(147, true)];
  const of = [dv.getFloat64(155, true), dv.getFloat64(163, true), dv.getFloat64(171, true)];
  const minLen = [20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67][fmt];
  if (recLen < minLen) throw new Error(`${fname}: point record length ${recLen} is shorter than format ${fmt} requires (${minLen} bytes).`);
  if (dataOff < headerSize || dataOff > bytes.length) throw new Error(`${fname}: the point data offset in the header (${dataOff}) lies outside the file (${bytes.length} bytes) — the file is truncated or corrupt.`);
  const avail = Math.floor((bytes.length - dataOff) / recLen);
  const warnings = [];
  if (avail < count) { warnings.push(`The header declares ${count} points but the file holds ${avail}; read ${avail}.`); count = avail; }

  // variable-length records: coordinate system
  let crs = null;
  {
    let p = headerSize;
    const keys = {};
    let dir = null, dbl = [], asc = '', wkt = null;
    for (let i = 0; i < nVlr && p + 54 <= bytes.length; i++) {
      const user = decodeText(bytes.subarray(p + 2, p + 18)).replace(/\u0000.*$/, ''), id = dv.getUint16(p + 18, true), len = dv.getUint16(p + 20, true);
      const body = p + 54;
      if (user === 'LASF_Projection') {
        if (id === 34735) { dir = []; for (let k = 0; k + 1 < len; k += 2) dir.push(dv.getUint16(body + k, true)); } else if (id === 34736) { for (let k = 0; k + 7 < len; k += 8) dbl.push(dv.getFloat64(body + k, true)); } else if (id === 34737) asc = decodeText(bytes.subarray(body, body + len));
        else if (id === 2112) wkt = decodeText(bytes.subarray(body, body + len)).replace(/\u0000+$/, '');
      }
      p = body + len;
    }
    if (wkt) { const info = parsePrj(wkt); if (info.kind !== 'unknown') crs = info; }
    if (!crs && dir && dir.length >= 4) {
      for (let k = 0; k < dir[3]; k++) { const o = 4 + 4 * k, id = dir[o], loc = dir[o + 1], cnt = dir[o + 2], vo = dir[o + 3]; if (loc === 0) keys[id] = vo; else if (loc === 34736) keys[id] = dbl[vo]; else if (loc === 34737) keys[id] = asc.substr(vo, cnt).replace(/\|$/, ''); }
      const info = geoKeysToCrs(keys);
      if (info.kind !== 'unknown') crs = info;
    }
  }
  const maxPts = opts.maxPoints ?? DEFAULT_MAX_POINTS;
  const stride = maxPts > 0 && count > maxPts ? Math.ceil(count / maxPts) : 1;
  const n = Math.ceil(count / stride);
  const coords = new Float64Array(n * 3);
  const intensity = new Uint16Array(n), cls = new Uint8Array(n);
  const rgbOff = { 2: 20, 3: 28, 5: 28, 7: 30, 8: 30, 10: 30 }[fmt];
  const rgb16 = rgbOff !== undefined ? new Uint16Array(n * 3) : null;
  const clsOff = fmt >= 6 ? 16 : 15;
  let anyI = false, maxC = 0, withheld = 0;
  for (let i = 0, k = 0; i < count; i += stride, k++) {
    const o = dataOff + i * recLen;
    coords[3 * k] = dv.getInt32(o, true) * sc[0] + of[0];
    coords[3 * k + 1] = dv.getInt32(o + 4, true) * sc[1] + of[1];
    coords[3 * k + 2] = dv.getInt32(o + 8, true) * sc[2] + of[2];
    const it = dv.getUint16(o + 12, true);
    intensity[k] = it;
    if (it) anyI = true;
    // formats 0–5 pack flags into the top three bits of the classification byte
    cls[k] = fmt >= 6 ? bytes[o + clsOff] : bytes[o + clsOff] & 0x1f;
    if ((fmt >= 6 ? bytes[o + 15] & 0x04 : bytes[o + clsOff] & 0x80)) withheld++;
    if (rgb16) { const r = dv.getUint16(o + rgbOff, true), g = dv.getUint16(o + rgbOff + 2, true), b = dv.getUint16(o + rgbOff + 4, true); rgb16[3 * k] = r; rgb16[3 * k + 1] = g; rgb16[3 * k + 2] = b; if (r > maxC) maxC = r; if (g > maxC) maxC = g; if (b > maxC) maxC = b; }
  }
  const attrs = { classification: cls };
  if (anyI) attrs.intensity = intensity;
  if (rgb16) { const rgb = new Uint8Array(n * 3); const sh = maxC > 255 ? 8 : 0; for (let i = 0; i < rgb.length; i++) rgb[i] = rgb16[i] >> sh; attrs.rgb = rgb; }
  if (stride > 1) warnings.push(`The cloud has ${count} points; every ${stride}${stride === 2 ? 'nd' : stride === 3 ? 'rd' : 'th'} point was kept (${n}) to stay under maxPoints.`);
  if (withheld) warnings.push(`${withheld} of the points read carry the "withheld" flag (kept; filter on classification if needed).`);
  const classes = {};
  for (let i = 0; i < n; i++) classes[cls[i]] = (classes[cls[i]] || 0) + 1;
  return cloudResult(baseName(input.name), coords, attrs, opts, {
    format: 'las', fileName: fname, crs, warnings, expected: 'the header declares zero point records', skipped: { thinnedPoints: count - n },
    extra: { version: `${vMajor}.${vMinor}`, pointFormat: fmt, pointCount: count, thinning: stride, scale: sc, offset: of, classes, systemId: decodeText(bytes.subarray(26, 58)).replace(/\u0000.*$/, '').trim(), verticalCrs: crs && crs.vertical ? crs.vertical : null },
  });
}

function textCloud(input, opts, layout) {
  const cur = lineCursor(input);
  const fname = input.name || layout.format;
  const xyz = new Grow(Float64Array, 3 * 4096), inten = new Grow(Float32Array, 4096), rgb = new Grow(Uint8Array, 3 * 4096), nrm = new Grow(Float32Array, 3 * 4096);
  let bad = 0, dropped = 0, scans = 0, declared = null, cols = 0;
  let M = null, hasI = false, hasRgb = false, hasN = false, rgbScale = 1;
  const push = (p) => {
    const n = p.length;
    if (!cols) {
      cols = n;
      hasI = n === 4 || n === 7 || n === 10;
      hasRgb = n === 6 || n === 7 || (n >= 9 && layout.format !== 'xyzn');
      hasN = n === 9 || n === 10 || (n === 6 && layout.normals);
      if (hasN && n === 6) hasRgb = false;
    }
    if (n < 3 || !(p[0] === p[0] && p[1] === p[1] && p[2] === p[2])) { bad++; return; }
    if (layout.format === 'ptx' && p[0] === 0 && p[1] === 0 && p[2] === 0) { dropped++; return; }
    if (M) { const q = M4.apply(M, p[0], p[1], p[2]); xyz.push3(q[0], q[1], q[2]); } else xyz.push3(p[0], p[1], p[2]);
    let k = 3;
    if (hasI) inten.push(p[k++] ?? 0);
    if (hasRgb) { rgb.push((p[k] ?? 0) * rgbScale); rgb.push((p[k + 1] ?? 0) * rgbScale); rgb.push((p[k + 2] ?? 0) * rgbScale); k += 3; }
    if (hasN) { nrm.push(p[k] ?? 0); nrm.push(p[k + 1] ?? 0); nrm.push(p[k + 2] ?? 0); }
  };
  const nums = (l) => l.split(/[\s,;]+/).filter(Boolean).map(Number);
  if (layout.format === 'ptx') {
    for (;;) {
      const a = cur.nextNonBlank();
      if (a === null) break;
      const b = cur.nextNonBlank();
      const nc = Number(a), nr = Number(b);
      if (!(nc >= 0 && nr >= 0) || b === null) throw new Error(`${fname}: not a PTX file (expected the number of columns and rows on the first two lines of each scan).`);
      for (let i = 0; i < 4; i++) cur.nextNonBlank(); // scanner position and axes; the 4 × 4 matrix below carries the same registration
      const m = [];
      for (let i = 0; i < 4; i++) { const l = cur.nextNonBlank(); if (l === null) throw new Error(`${fname}: PTX scan ${scans + 1} ends inside its header.`); m.push(...nums(l)); }
      // rows of the file are the columns of the transform (points are row vectors in PTX)
      M = Float64Array.from(m);
      scans++;
      cols = 0;
      for (let i = 0; i < nc * nr; i++) { const l = cur.nextNonBlank(); if (l === null) break; push(nums(l)); }
    }
  } else {
    for (let l = cur.nextNonBlank(); l !== null; l = cur.nextNonBlank()) {
      if (l[0] === '#' || l.startsWith('//')) continue;
      const p = nums(l);
      if (declared === null && p.length === 1 && Number.isInteger(p[0]) && !xyz.n) { declared = p[0]; continue; }
      if (!cols && p.length >= 6 && layout.format !== 'pts') {
        // colours given as 0–1 floats rather than 0–255 integers
        const c0 = p.length === 7 || p.length === 10 ? 4 : 3;
        if (!layout.normals && p.length !== 9 && p[c0] <= 1 && p[c0 + 1] <= 1 && p[c0 + 2] <= 1 && !(Number.isInteger(p[c0]) && Number.isInteger(p[c0 + 1]) && Number.isInteger(p[c0 + 2]))) rgbScale = 255;
      }
      push(p);
    }
  }
  const n = xyz.n / 3;
  const attrs = {};
  if (hasI && inten.n === n) attrs.intensity = inten.finish();
  if (hasRgb && rgb.n === 3 * n) attrs.rgb = rgb.finish();
  if (hasN && nrm.n === 3 * n) attrs.normals = nrm.finish();
  const warnings = [];
  if (bad) warnings.push(`${bad} line(s) without three numeric coordinates were skipped.`);
  if (declared !== null && declared !== n + bad) warnings.push(`The first line declares ${declared} points; ${n} were read.`);
  return cloudResult(baseName(input.name), xyz.finish(), attrs, opts, {
    format: layout.format === 'xyzn' ? 'xyzrgb' : layout.format, fileName: fname, warnings, expected: 'expected lines of "x y z [intensity] [r g b]"', skipped: { invalidLines: bad, emptyReturns: dropped },
    extra: { columns: cols, scans: scans || undefined, attributes: Object.keys(attrs) },
  });
}

/** Leica PTS: optional point count, then "x y z intensity r g b" rows. */
export function parsePts(input, opts = {}) { return textCloud(input, opts, { format: 'pts' }); }
/** Leica PTX: one or more registered scans. */
export function parsePtx(input, opts = {}) { return textCloud(input, opts, { format: 'ptx' }); }
/** Header-less x y z plus colour / intensity / normals (4, 6, 7, 9 or 10 columns). */
export function parseXyzRgb(input, opts = {}) { return textCloud(input, opts, { format: /\.xyzn$/i.test(input.name || '') ? 'xyzn' : 'xyzrgb', normals: /\.xyzn$/i.test(input.name || '') || opts.columns === 'normals' }); }
