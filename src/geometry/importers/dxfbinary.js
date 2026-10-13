// AutoCAD DXF, binary form (AutoCAD Release 10 onwards).
//
// Layout (Autodesk DXF reference, "Binary DXF files"): the 22-byte sentinel "AutoCAD Binary DXF\r\n\x1a\0",
// then (group code, value) pairs to the end of the file. The group code is a 16-bit little-endian integer in
// files written by R13 and later, and a single byte in R12 and earlier, where the byte 255 announces a 16-bit
// code (used for extended data, codes 1000 and above). The value type follows from the group code exactly as
// in the ASCII form: null-terminated strings, 16-, 32- and 64-bit little-endian integers, IEEE doubles,
// one-byte booleans (290–299) and binary chunks with a leading length byte (310–319 and 1004).
//
// The pairs are decoded into the (code, trimmed string) stream of the ASCII form and handed to the entity
// handling of dxf.js, so both encodings of a drawing give the same models — including the solids stored as
// ACIS data (SAT text in groups 1 / 3, or SAB chunks in the ACDSDATA section, see dxf.js). Numbers are written with the
// shortest decimal text that reads back to the same double.

import { dxfModels } from './dxf.js';
import { inputBytes } from './common.js';

const SENTINEL = 'AutoCAD Binary DXF\r\n\x1a\0';
const STR = 0, F64 = 1, I16 = 2, I32 = 3, I64 = 4, BOOL = 5, CHUNK = 6;

/** Value type of a group code (DXF reference, "Group code value types"); −1 for codes outside the tables. */
export function groupType(code) {
  if (code < 0) return code >= -5 ? STR : -1; // application-specific markers, never written to a file
  if (code <= 9) return STR;
  if (code <= 59) return F64;
  if (code <= 79) return I16;
  if (code >= 90 && code <= 99) return I32;
  // 101 is the embedded-object marker, met in the ACDSDATA section that holds the binary ACIS data of solids
  if (code === 100 || code === 101 || code === 102 || code === 105) return STR;
  if (code >= 110 && code <= 149) return F64;
  if (code >= 160 && code <= 169) return I64;
  if (code >= 170 && code <= 179) return I16;
  if (code >= 210 && code <= 239) return F64;
  if (code >= 270 && code <= 289) return I16;
  if (code >= 290 && code <= 299) return BOOL;
  if (code >= 300 && code <= 309) return STR;
  if (code >= 310 && code <= 319) return CHUNK;
  if (code >= 320 && code <= 369) return STR;
  if (code >= 370 && code <= 389) return I16;
  if (code >= 390 && code <= 399) return STR;
  if (code >= 400 && code <= 409) return I16;
  if (code >= 410 && code <= 419) return STR;
  if (code >= 420 && code <= 429) return I32;
  if (code >= 430 && code <= 439) return STR;
  if (code >= 440 && code <= 459) return I32;
  if (code >= 460 && code <= 469) return F64;
  if (code >= 470 && code <= 481) return STR;
  if (code === 999) return STR;
  if (code === 1004) return CHUNK;
  if (code >= 1000 && code <= 1009) return STR;
  if (code >= 1010 && code <= 1059) return F64;
  if (code >= 1060 && code <= 1070) return I16;
  if (code === 1071) return I32;
  return -1;
}

class Truncated extends Error {}

/**
 * Decode the pair stream. `wide` selects 16-bit group codes (R13+) or single-byte codes with the 255 escape.
 * `limit` stops after that many pairs (used to probe the code width).
 * @returns {{ pairs: Array, eof: boolean, count: number }}
 */
function decode(bytes, wide, limit = Infinity) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const end = bytes.length;
  const utf8 = new TextDecoder('utf-8');
  const pairs = [];
  let o = SENTINEL.length, count = 0, eof = false;
  const need = (n, code) => { if (o + n > end) throw new Truncated(`the file ends inside the value of group ${code} at byte ${o} (pair ${count + 1})`); };
  while (o < end && count < limit) {
    const at = o;
    let code;
    if (wide) { if (o + 2 > end) throw new Truncated(`the file ends inside a group code at byte ${o} (pair ${count + 1})`); code = dv.getInt16(o, true); o += 2; } else {
      code = bytes[o++];
      if (code === 255) { if (o + 2 > end) throw new Truncated(`the file ends inside a group code at byte ${at} (pair ${count + 1})`); code = dv.getInt16(o, true); o += 2; }
    }
    let v;
    switch (groupType(code)) {
      case STR: {
        let e = o;
        while (e < end && bytes[e] !== 0) e++;
        if (e >= end) throw new Truncated(`the file ends inside the string of group ${code} at byte ${o} (pair ${count + 1})`);
        // groups 1 and 3 keep their blanks: embedded ACIS text is split across them at any character
        v = utf8.decode(bytes.subarray(o, e));
        if (code !== 1 && code !== 3) v = v.trim();
        o = e + 1;
        break;
      }
      case F64: need(8, code); v = String(dv.getFloat64(o, true)); o += 8; break;
      case I16: need(2, code); v = String(dv.getInt16(o, true)); o += 2; break;
      case I32: need(4, code); v = String(dv.getInt32(o, true)); o += 4; break;
      case I64: need(8, code); v = String(dv.getBigInt64(o, true)); o += 8; break;
      case BOOL: need(1, code); v = String(bytes[o]); o += 1; break;
      case CHUNK: {
        need(1, code);
        const n = bytes[o];
        need(1 + n, code);
        let hex = '';
        for (let k = 0; k < n; k++) hex += bytes[o + 1 + k].toString(16).padStart(2, '0').toUpperCase();
        v = hex;
        o += 1 + n;
        break;
      }
      default: throw new RangeError(`group code ${code} at byte ${at} (pair ${count + 1}) is not in the DXF group-code tables, so the length of its value is unknown`);
    }
    pairs.push(code, v);
    count++;
    if (code === 0 && v === 'EOF') { eof = true; break; }
  }
  return { pairs, eof, count };
}

/** True when the first pairs read as "0 SECTION / 2 <name>" with this code width. */
function plausible(bytes, wide) {
  try {
    const { pairs } = decode(bytes, wide, 2);
    return pairs.length === 4 && pairs[0] === 0 && pairs[1] === 'SECTION' && pairs[2] === 2 && /^[A-Z_]+$/.test(pairs[3]);
  } catch { return false; }
}

/**
 * Decode a binary DXF file into its group pairs.
 * @returns {{ pairs: Array, codeBytes: 1|2 }} flat [code, value, …] with values as the ASCII form would hold them
 */
export function readBinaryDxfPairs(input) {
  const fname = input.name || 'DXF';
  const bytes = inputBytes(input);
  let ok = bytes.length >= SENTINEL.length;
  for (let i = 0; ok && i < SENTINEL.length; i++) if (bytes[i] !== SENTINEL.charCodeAt(i)) ok = false;
  if (!ok) throw new Error(`${fname}: not a binary DXF file (the 22-byte sentinel "AutoCAD Binary DXF" + CR LF SUB NUL is missing or damaged).`);
  if (bytes.length === SENTINEL.length) throw new Error(`${fname}: the binary DXF file is truncated: it ends right after the sentinel, before the first group.`);
  const wide = plausible(bytes, true) ? true : plausible(bytes, false) ? false : null;
  if (wide === null) {
    if (bytes.length < SENTINEL.length + 24) throw new Error(`${fname}: the binary DXF file is truncated: it ends inside the first groups, before any section.`);
    throw new Error(`${fname}: the data after the binary DXF sentinel do not start with the groups "0 SECTION / 2 <name>" in either the 16-bit (R13 and later) or the single-byte (R12 and earlier) group-code layout.`);
  }
  let res;
  try { res = decode(bytes, wide); } catch (e) {
    if (e instanceof Truncated) throw new Error(`${fname}: the binary DXF file is truncated: ${e.message}.`);
    throw new Error(`${fname}: binary DXF could not be decoded: ${e.message}.`);
  }
  if (!res.eof) throw new Error(`${fname}: the binary DXF file is truncated: it ends after ${res.count} groups without the closing "0 EOF" group.`);
  return { pairs: res.pairs, codeBytes: wide ? 2 : 1 };
}

/** AutoCAD DXF (binary). */
export function parseDxfBinary(input, opts = {}) {
  const fname = input.name || 'DXF';
  const { pairs, codeBytes } = readBinaryDxfPairs(input);
  const sec = { HEADER: [], TABLES: [], BLOCKS: [], ENTITIES: [], ACDSDATA: [] };
  let target = null;
  for (let i = 0; i < pairs.length; i += 2) {
    const code = pairs[i], v = pairs[i + 1];
    if (code === 0 && v === 'SECTION') { i += 2; target = i < pairs.length ? sec[String(pairs[i + 1]).trim().toUpperCase()] || null : null; continue; }
    if (code === 0 && v === 'ENDSEC') { target = null; continue; }
    if (code === 0 && v === 'EOF') break;
    if (target) target.push(code, v);
  }
  if (!sec.ENTITIES.length && !sec.BLOCKS.length) throw new Error(`${fname}: the binary DXF drawing has no ENTITIES or BLOCKS section with content.`);
  const r = dxfModels(sec, input, opts);
  r.meta.format = 'dxf-binary';
  r.meta.encoding = `binary, ${codeBytes === 2 ? '16-bit group codes (R13 and later)' : 'single-byte group codes (R12 and earlier)'}`;
  return r;
}
