// ACIS SAB — the binary save format of the ACIS / Autodesk ShapeManager kernel (".sab", and the solid data of
// 3DSOLID, REGION, BODY and SURFACE entities in DXF R2013 and later).
//
// The file holds the same records as the SAT text form, written as a stream of tagged values:
//   signature  15 bytes: "ACIS BinaryFile", "ASM BinaryFile4" or "ASM BinaryFile8" (the last with 64-bit integers)
//   header     four integers without tags (version, record count, body count, flags), then tagged values:
//              product, kernel release and date as strings, millimetres per unit, resabs, resnor
//   records    one tag byte per value:
//                0x02 char   0x03 16-bit integer   0x04 integer (32 or 64 bit, as the signature says)
//                0x05 float  0x06 double           0x17 64-bit integer
//                0x07 / 0x08 / 0x09 string with an 8-, 16- or 32-bit length;  0x12 string with a 32-bit length
//                0x0A true, 0x0B false — the logical value of a keyword pair of the text form: reversed / forward,
//                     double / single, in / out, reverse_v / forward_v, and "F" (a finite interval end, followed
//                     by its value) / "I"
//                0x0C pointer to a record (index, −1 = null)
//                0x0E part of an entity name, 0x0D its last part ("0E plane 0D surface" = plane-surface)
//                0x0F opens and 0x10 closes a sub-type block ("{" and "}"), 0x11 ends the record ("#")
//                0x13 position, 0x14 direction (three doubles), 0x16 parameter pair (two doubles)
//                0x15 value of an enumeration (integer; the text form writes a keyword)
//              The data ends with the record name End-of-ACIS-data / End-of-ASM-data, or where a history
//              section begins.
// The tag values and the header layout follow the description of the format in the ezdxf project
// (https://github.com/mozman/ezdxf, src/ezdxf/acis/sab.py, MIT licence) and, for the 64-bit variant and the
// string and vector tags ezdxf does not meet, the notes in InventorLoader (https://github.com/jmplonka/
// InventorLoader). The decoder below is written for this application from that description.
//
// The records are converted to the token lists of acis.js (readSat) — logical values and enumerations stay what
// they are, since their keyword depends on the field — and the geometry is then read by the same code as SAT.

import { inputBytes } from './common.js';
import { acisModels, END_MARK } from './acis.js';

const SIGNATURES = ['ACIS BinaryFile', 'ASM BinaryFile4', 'ASM BinaryFile8'];
const IDENT = /^[A-Za-z_][A-Za-z0-9_]*$/;

/** True when the bytes start with one of the SAB signatures. */
export function isSab(bytes) {
  if (!bytes || bytes.length < 15) return false;
  let s = '';
  for (let i = 0; i < 15; i++) s += String.fromCharCode(bytes[i]);
  return SIGNATURES.includes(s);
}

/**
 * Decode a SAB file into the header and records that acis.js reads.
 * @param {Uint8Array} bytes @param {string} [fname]
 * @returns {{ header: object, records: Array<{ type: string, kind: string, tok: object[], idx: number }>, wide: boolean }}
 */
export function readSab(bytes, fname = 'SAB') {
  if (!isSab(bytes)) throw new Error(`${fname}: not an ACIS binary (SAB) file: the signature "ACIS BinaryFile" / "ASM BinaryFile4" is missing.`);
  const wide = bytes[14] === 0x38;
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const end = bytes.length;
  const utf8 = new TextDecoder('utf-8');
  let o = 15;
  const truncated = (what) => new Error(`${fname}: the SAB file is truncated: it ends inside ${what} at byte ${o}.`);
  const need = (n, what) => { if (o + n > end) throw truncated(what); };
  const int = (what) => {
    if (wide) { need(8, what); const v = Number(dv.getBigInt64(o, true)); o += 8; return v; }
    need(4, what); const v = dv.getInt32(o, true); o += 4; return v;
  };
  const f64 = (what) => { need(8, what); const v = dv.getFloat64(o, true); o += 8; return v; };
  const str = (lenBytes, what) => {
    need(lenBytes, what);
    const n = lenBytes === 1 ? bytes[o] : lenBytes === 2 ? dv.getUint16(o, true) : dv.getUint32(o, true);
    o += lenBytes;
    need(n, what);
    const s = utf8.decode(bytes.subarray(o, o + n));
    o += n;
    return s;
  };

  const header = { version: int('the header'), records: int('the header'), bodies: int('the header'), flag: int('the header'), product: null, acis: null, date: null, mmPerUnit: null, resabs: null, resnor: null };
  if (!(header.version >= 100 && header.version <= 99999)) throw new Error(`${fname}: the SAB header gives the version ${header.version}, which is not an ACIS save version.`);
  const tagged = (want, what) => {
    need(1, what);
    const tag = bytes[o];
    if (want === 's' && (tag === 0x07 || tag === 0x08 || tag === 0x09 || tag === 0x12)) { o++; return str(tag === 0x07 ? 1 : tag === 0x08 ? 2 : 4, what).trim(); }
    if (want === 'd' && tag === 0x06) { o++; return f64(what); }
    throw new Error(`${fname}: the SAB header is not laid out as expected (tag 0x${tag.toString(16)} at byte ${o} where ${want === 's' ? 'a string' : 'a double'} should be).`);
  };
  header.product = tagged('s', 'the header'); header.acis = tagged('s', 'the header'); header.date = tagged('s', 'the header');
  header.mmPerUnit = tagged('d', 'the header'); header.resabs = tagged('d', 'the header'); header.resnor = tagged('d', 'the header');

  const records = [];
  let stop = false;
  while (o < end && !stop) {
    // optional stream identifier some writers put after the header: a logical value, then a string when true
    const parts = [];
    const tok = [];
    let type = null, depth = 0, closed = false;
    const what = () => `record ${records.length}${type ? ` (${type})` : ''}`;
    while (o < end) {
      const at = o, tag = bytes[o++];
      if (type === null) {
        if (tag === 0x0e) { parts.push(str(1, what())); continue; }
        if (tag === 0x0d) { parts.push(str(1, what())); type = parts.join('-'); if (END_MARK.test(type)) { stop = true; closed = true; break; } continue; }
        // values between the header and the first record (stream identifier) are passed over
        if (!records.length && !parts.length && (tag === 0x0a || tag === 0x0b)) continue;
        if (!records.length && !parts.length && (tag === 0x07 || tag === 0x08 || tag === 0x09 || tag === 0x12)) { str(tag === 0x07 ? 1 : tag === 0x08 ? 2 : 4, 'the header'); continue; }
        if (!records.length && !parts.length && tag === 0x04) { int('the header'); continue; }
        throw new Error(`${fname}: record ${records.length} of the SAB file does not start with an entity name (tag 0x${tag.toString(16)} at byte ${at}).`);
      }
      if (tag === 0x11) { closed = true; break; }
      switch (tag) {
        case 0x02: need(1, what()); tok.push({ t: 'n', v: bytes[o] }); o += 1; break;
        case 0x03: need(2, what()); tok.push({ t: 'n', v: dv.getInt16(o, true) }); o += 2; break;
        case 0x04: tok.push({ t: 'n', v: int(what()) }); break;
        case 0x05: need(4, what()); tok.push({ t: 'n', v: dv.getFloat32(o, true) }); o += 4; break;
        case 0x06: tok.push({ t: 'n', v: f64(what()) }); break;
        case 0x17: need(8, what()); tok.push({ t: 'n', v: Number(dv.getBigInt64(o, true)) }); o += 8; break;
        case 0x07: case 0x08: case 0x09: case 0x12: {
          const s = str(tag === 0x07 ? 1 : tag === 0x08 ? 2 : 4, what());
          if (tag === 0x12 && type === 'transform') {
            // the matrix is stored as the text of the SAT form
            for (const w of s.trim().split(/\s+/)) if (w) tok.push(/^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/.test(w) ? { t: 'n', v: Number(w) } : { t: 'w', v: w });
          } else tok.push(IDENT.test(s) && s.length <= 40 ? { t: 'w', v: s } : { t: 's', v: s });
          break;
        }
        case 0x0a: tok.push({ t: 'b', v: true }); break;
        case 0x0b: tok.push({ t: 'b', v: false }); break;
        case 0x0c: tok.push({ t: 'p', v: int(what()) }); break;
        case 0x0d: case 0x0e: tok.push({ t: 'w', v: str(1, what()) }); break;
        case 0x0f: depth++; tok.push({ t: '{' }); break;
        case 0x10: depth--; tok.push({ t: '}' }); break;
        case 0x13: case 0x14: tok.push({ t: 'n', v: f64(what()) }, { t: 'n', v: f64(what()) }, { t: 'n', v: f64(what()) }); break;
        case 0x16: tok.push({ t: 'n', v: f64(what()) }, { t: 'n', v: f64(what()) }); break;
        case 0x15: tok.push({ t: 'e', v: int(what()) }); break;
        default: throw new Error(`${fname}: unknown value tag 0x${tag.toString(16)} at byte ${at} in ${what()} of the SAB file; the length of its value is not known, so the rest cannot be read.`);
      }
    }
    if (type === null) { if (parts.length) throw truncated(`the name of record ${records.length}`); break; }
    if (stop) break;
    if (!closed) throw truncated(what());
    records.push({ type, kind: type.slice(type.lastIndexOf('-') + 1), tok, idx: records.length });
  }
  if (!records.length) throw new Error(`${fname}: the SAB file holds a header but no entity records.`);
  return { header, records, wide };
}

/**
 * ACIS SAB (binary).
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ units?: string, scale?: number, deflection?: number, angularStep?: number }} [opts] as for SAT (acis.js)
 */
export function parseSab(input, opts = {}) {
  const fname = input.name || 'SAB';
  const sat = readSab(inputBytes(input), fname);
  const r = acisModels(sat, input, opts, 'acis-sab');
  r.meta.encoding = `binary, ${sat.wide ? '64' : '32'}-bit integers`;
  return r;
}
