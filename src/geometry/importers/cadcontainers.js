// Native CAD part files that are containers around a kernel stream this application reads: SolidWorks parts
// (.sldprt) and Solid Edge parts (.par, .psm, .pwd) around Parasolid, Autodesk Inventor parts (.ipt) around ACIS.
// The container is opened here and the stream is handed to the Parasolid reader (parasolid.js) or the SAB reader
// (sab.js), so what is read of the geometry — and what is not (rolling-ball blends, offsets …) — is exactly what
// those readers read of an .x_b or .sab file. Nothing else of the file is interpreted: no feature tree,
// no sketches, no configurations other than the stored bodies, no colours or materials.
//
// SolidWorks, files up to 2014: a compound (OLE2) file. Each configuration that has been rebuilt keeps its bodies
//   in the stream "Contents/Config-<n>-Partition": one or more blocks of
//       [16-byte marker 23 1D D5 71 DA 81 48 A2 A8 58 98 B2 1B 89 EF 99] { [uncompressed size][compressed size][zlib data] }
//   (large bodies in several pieces, the first of 1 MiB) whose joined data are neutral-binary Parasolid transmit files without the "**" header — first the partition (the
//   bodies), then its deltas (rollback states, not read). "…-GhostPartition" holds bodies that are not part of the
//   model and is not read.
// SolidWorks 2015 and later: not a compound file but a run of records, each
//       [14 00 06 00 08 00][4 bytes][CRC-32][compressed size][uncompressed size][name length][name][raw deflate data]
//   (32-bit little-endian numbers; every byte of the name has its two halves swapped). The names are the stream
//   paths of the older files and the data of "Contents/Config-<n>-Partition" are the blocks described above.
//   Records whose data do not inflate to the announced size are filler and are passed over.
//   A part made of imported geometry has an empty partition; its bodies are kept, one per import feature and not
//   knitted, in "Config-<n>-FeatureBodies/LocalBodies". That stream is recognised and the file refused with this
//   explanation: which of those bodies make up the model is not known to this reader.
// Solid Edge: a compound file with the stream "PARASOLID/STREAM<n>.P_B" — a neutral-binary partition file with its
//   "**" header, uncompressed; "….D_B" holds the deltas.
//
// Autodesk Inventor part (.ipt): a compound file whose storage "RSeStorage" holds one stream "B<id>" per segment of
//   the document and a stream "M<id>" that names it. The segment named "…BRepSegment" keeps the bodies of the
//   model as a standard binary ACIS stream ("ASM BinaryFile…" to its end marker) inside other records; it is found
//   by its signature and handed to the SAB reader (sab.js / acis.js). A "B" stream is a 16-byte identifier, the
//   bytes 04 and the compression (01: zlib, in files of Inventor 2013 – 2025; 02: Zstandard, in files of Inventor
//   2026, which is not decoded here), then the compressed segment. Sketches, work geometry and feature bodies of
//   the other segments
//   are not read.
//
// None of these layouts is published by its vendor. They were worked out from files and are checked against the
// Parasolid or STEP export of the same parts (docs/geometry.md, "Checks on real files"). Assemblies (.sldasm,
// .asm) hold references to their part files, not bodies, and are refused with that explanation.

import { inputBytes, finish } from './common.js';
import { isCfb, readCfb } from '../cfb.js';
import { inflateRaw, inflateZlib } from '../zip.js';
import { parseParasolid } from './parasolid.js';
import { parseSab } from './sab.js';

const MARKER = [0x23, 0x1d, 0xd5, 0x71, 0xda, 0x81, 0x48, 0xa2, 0xa8, 0x58, 0x98, 0xb2, 0x1b, 0x89, 0xef, 0x99];
const RECORD = [0x14, 0x00, 0x06, 0x00, 0x08, 0x00];
const u32 = (b, o) => (b[o] | (b[o + 1] << 8) | (b[o + 2] << 16)) + b[o + 3] * 16777216;
const matchAt = (b, o, sig) => { for (let i = 0; i < sig.length; i++) if (b[o + i] !== sig[i]) return false; return true; };
const indexOf = (b, sig, from) => { for (let o = from, e = b.length - sig.length; o <= e; o++) if (b[o] === sig[0] && matchAt(b, o, sig)) return o; return -1; };
const isPs = (b) => b.length > 8 && b[0] === 0x50 && b[1] === 0x53 && b[2] === 0 && b[3] === 0;
/** "partition", "deltas" or "" (a part file) from the modeller-version string that opens a binary transmit file. */
const psGuise = (b) => { const n = Math.min(b.length, 6 + ((b[4] << 8) | b[5])); let s = ''; for (let i = 6; i < n; i++) s += String.fromCharCode(b[i]); const m = /\((\w+)\)/.exec(s); return m ? m[1] : ''; };

/** The Parasolid transmit files in a SolidWorks "…Partition" stream, inflated: `[{ guise, bytes }]`. */
async function partitionBlocks(stream, what) {
  const out = [];
  if (isPs(stream)) return [{ guise: psGuise(stream), bytes: stream }];
  let p = indexOf(stream, MARKER, 0);
  while (p >= 0) {
    // after the marker: pieces of [uncompressed size][compressed size][zlib data] (the first of up to 1 MiB, then
    // of 4 KiB), closed by a piece of size zero
    const pieces = [];
    let q = p + 16, total = 0;
    while (q + 8 <= stream.length) {
      const size = u32(stream, q), packed = u32(stream, q + 4);
      if (!packed || !size || q + 8 + packed > stream.length || stream[q + 8] !== 0x78) break;
      let data;
      try { data = await inflateZlib(stream.subarray(q + 8, q + 8 + packed), size); } catch { throw new Error(`${what}: a compressed block of the body data could not be inflated — the file is damaged.`); }
      if (data.length !== size) throw new Error(`${what}: a compressed block of the body data inflates to ${data.length} bytes where ${size} are announced — the file is damaged.`);
      pieces.push(data);
      total += size;
      q += 8 + packed;
    }
    if (pieces.length) {
      let data = pieces[0];
      if (pieces.length > 1) { data = new Uint8Array(total); let o = 0; for (const c of pieces) { data.set(c, o); o += c.length; } }
      if (isPs(data)) out.push({ guise: psGuise(data), bytes: data });
    }
    p = indexOf(stream, MARKER, Math.max(q, p + 16));
  }
  return out;
}

/** Streams of a SolidWorks 2015+ file whose name passes `want`: `Map<name, Uint8Array>` (last record of a name wins). */
async function sw2015Streams(bytes, want) {
  const found = new Map();
  let records = 0;
  for (let p = indexOf(bytes, RECORD, 0); p >= 0; p = indexOf(bytes, RECORD, p + 6)) {
    if (p + 26 > bytes.length) break;
    const packed = u32(bytes, p + 14), size = u32(bytes, p + 18), nlen = u32(bytes, p + 22);
    if (nlen === 0 || nlen > 512 || p + 26 + nlen + packed > bytes.length) continue;
    let name = '';
    for (let i = 0; i < nlen; i++) { const c = bytes[p + 26 + i]; name += String.fromCharCode(((c << 4) & 0xf0) | (c >> 4)); }
    if (!/^[\x20-\x7e]+$/.test(name)) continue;
    records++;
    if (!want(name)) continue;
    let data = null;
    try { data = await inflateRaw(bytes.subarray(p + 26 + nlen, p + 26 + nlen + packed), size); } catch { data = null; }
    // the same name occurs in filler records, whose data are not a deflate stream of the announced size
    if (data && data.length === size) found.set(name, data);
  }
  return { found, records };
}

function tag(result, fname, extra) {
  for (const k in extra) result.meta[k] = extra[k];
  result.meta.fileName = fname;
  return result;
}

const PARTITION = /^Contents\/Config-(\d+)-Partition$/;
const LOCAL = /(^|\/)Config-\d+-FeatureBodies\/LocalBodies$/;

/**
 * SolidWorks part (.sldprt).
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ configuration?: number, deflection?: number, angularStep?: number }} [opts] `configuration`: number n of
 *   the "Config-n" whose bodies are wanted (default: the lowest that has body data); the rest as for Parasolid
 */
export async function parseSolidWorks(input, opts = {}) {
  const bytes = inputBytes(input);
  const fname = input.name || 'SLDPRT';
  const asm = /\.sldasm$/i.test(fname), drw = /\.slddrw$/i.test(fname);
  let streams, container, names, local = false;
  if (isCfb(bytes)) {
    const cfb = readCfb(bytes, fname);
    names = cfb.entries.filter((e) => !e.storage).map((e) => e.path);
    if (!names.some((n) => /^Contents\//.test(n) || n === 'Header2' || /^_MO_VERSION_/.test(n))) throw new Error(`${fname}: this compound file has none of the streams of a SolidWorks document (Contents/…, Header2); it was written by another program.`);
    streams = new Map(cfb.find(PARTITION).map((e) => [e.path, cfb.read(e.path)]));
    local = names.some((n) => LOCAL.test(n));
    container = 'SolidWorks compound file (2014 or earlier)';
  } else {
    if (indexOf(bytes, RECORD, 0) < 0) throw new Error(`${fname}: not a SolidWorks file — neither a compound (OLE2) file nor the record stream of SolidWorks 2015 and later was found. If this is another program's part file, export STEP from that program.`);
    const r = await sw2015Streams(bytes, (n) => PARTITION.test(n) || LOCAL.test(n));
    for (const n of [...r.found.keys()]) if (LOCAL.test(n)) { local = true; r.found.delete(n); }
    if (!r.records) throw new Error(`${fname}: not a SolidWorks file — no named records were found. If this is another program's part file, export STEP from that program.`);
    streams = r.found;
    container = 'SolidWorks record file (2015 or later)';
  }
  const configs = [...streams.keys()].map((n) => Number(PARTITION.exec(n)[1])).sort((a, b) => a - b);
  if (!configs.length) {
    if (asm) throw new Error(`${fname}: a SolidWorks assembly stores references to its part files and no bodies of its own. Import the part files (.sldprt), or export the assembly as STEP or Parasolid from SolidWorks.`);
    if (drw) throw new Error(`${fname}: a SolidWorks drawing holds views, not bodies. Export the drawing as DXF, or import the part it shows.`);
    throw new Error(`${fname}: the file holds no body data (no "Contents/Config-n-Partition" stream) — it may have been saved without its bodies. Open it in SolidWorks, rebuild and save it, or export STEP or Parasolid.`);
  }
  const cfg = opts.configuration ?? configs[0];
  if (!configs.includes(cfg)) throw new RangeError(`configuration: the file has body data for configuration${configs.length > 1 ? 's' : ''} ${configs.join(', ')}, not for ${cfg}.`);
  const blocks = await partitionBlocks(streams.get(`Contents/Config-${cfg}-Partition`), fname);
  const parts = blocks.filter((b) => b.guise !== 'deltas');
  if (!parts.length) throw new Error(`${fname}: the body data of configuration ${cfg} hold no Parasolid partition (${blocks.length} block(s) found) — the layout of this SolidWorks version is not known to this reader. Export STEP or Parasolid from SolidWorks.`);
  let r;
  try { r = parseParasolid({ name: fname, bytes: parts[0].bytes }, opts); } catch (e) {
    // parts made of imported geometry keep one Parasolid body per import feature in "…FeatureBodies/LocalBodies"
    // and an empty partition; the bodies are there, but not which of them make up the model, nor knitted together
    if (/no geometry could be read from the Parasolid file\.$/.test(e.message) && local) throw new Error(`${fname}: this part is made of imported geometry: SolidWorks keeps its bodies with the import features (stream "Config-${cfg}-FeatureBodies/LocalBodies"), not as a finished model, and which of them make up the part cannot be told from the file by this reader. Export STEP or Parasolid from SolidWorks.`);
    if (/no geometry could be read from the Parasolid file\.$/.test(e.message)) throw new Error(`${fname}: the body data of configuration ${cfg} hold no bodies — the part was saved without solid or surface bodies (all of them deleted or suppressed, or only display data kept). Open it in SolidWorks, rebuild and save it, or export STEP or Parasolid.`);
    throw e;
  }
  // the partition holds every body of the configuration; whether one is hidden is kept in the feature data, not read
  if (r.meta.bodies > 1) r.warnings.push(`The part holds ${r.meta.bodies} bodies. Bodies hidden in SolidWorks cannot be told from shown ones in the file and are included; each body is a model of its own and can be removed.`);
  if (configs.length > 1) r.warnings.push(`The file stores bodies for ${configs.length} configurations (${configs.join(', ')}); configuration ${cfg} was read (opts.configuration chooses another).`);
  if (parts.length > 1) r.warnings.push(`The body data hold ${parts.length} Parasolid partitions; the first was read.`);
  return tag(r, fname, { format: 'solidworks', container, kernelStream: `Contents/Config-${cfg}-Partition`, configuration: cfg, configurations: configs });
}

const SE_STREAM = /^PARASOLID\/STREAM(\d+)\.P_B$/i;

/**
 * Solid Edge part, sheet-metal part or weldment (.par, .psm, .pwd).
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ deflection?: number, angularStep?: number }} [opts] as for Parasolid
 */
export function parseSolidEdge(input, opts = {}) {
  const bytes = inputBytes(input);
  const fname = input.name || 'PAR';
  if (!isCfb(bytes)) throw new Error(`${fname}: not a Solid Edge file (a Solid Edge document is a compound (OLE2) file). If this is another program's part file, export STEP from that program.`);
  const cfb = readCfb(bytes, fname);
  const found = cfb.find(SE_STREAM).sort((a, b) => Number(SE_STREAM.exec(a.path)[1]) - Number(SE_STREAM.exec(b.path)[1]));
  if (!found.length) {
    const names = cfb.entries.map((e) => e.path);
    if (!names.some((n) => /^PSM|^JSitesList$|^BuildVersions$/.test(n))) throw new Error(`${fname}: this compound file has none of the streams of a Solid Edge document; it was written by another program.`);
    if (/\.asm$/i.test(fname)) throw new Error(`${fname}: a Solid Edge assembly stores references to its part files and no bodies of its own. Import the part files (.par, .psm), or export the assembly as STEP or Parasolid from Solid Edge.`);
    throw new Error(`${fname}: the file holds no body data (no "PARASOLID/STREAMn.P_B" stream) — a draft, an assembly or an empty part. Export STEP or Parasolid from Solid Edge.`);
  }
  const results = [];
  const failed = [];
  for (const e of found) {
    try { results.push([e.path, parseParasolid({ name: fname, bytes: cfb.read(e.path) }, opts)]); } catch (err) { failed.push(`${e.path}: ${err.message.replace(`${fname}: `, '')}`); }
  }
  if (!results.length) throw new Error(`${fname}: ${failed.join(' ')}`);
  const first = results[0][1];
  if (results.length === 1 && !failed.length) return tag(first, fname, { format: 'solid-edge', container: 'Solid Edge compound file', kernelStream: results[0][0] });
  // several partitions: the models side by side, the counts of the first kept and the rest listed
  const models = [], warnings = [];
  for (const [path, r] of results) { models.push(...r.models); for (const w of r.warnings) warnings.push(results.length > 1 ? `${path}: ${w}` : w); }
  for (const f of failed) warnings.push(`Not read — ${f}`);
  const sum = (k) => results.reduce((a, x) => a + (x[1].meta[k] || 0), 0);
  const elementTypes = {}, skipped = {};
  for (const [, r] of results) { for (const k in r.meta.elementTypes) elementTypes[k] = (elementTypes[k] || 0) + r.meta.elementTypes[k]; for (const k in r.meta.skipped) skipped[k] = (skipped[k] || 0) + r.meta.skipped[k]; }
  if (failed.length) skipped['partition that could not be interpreted'] = failed.length;
  return tag(finish(models, warnings, { ...first.meta, counts: undefined, bbox: undefined, elementTypes, skipped, bodies: sum('bodies'), instances: sum('instances'), solidBodies: sum('solidBodies'), closedBodies: sum('closedBodies'), orientedBodies: sum('orientedBodies'), faces: sum('faces'), faceTotal: sum('faceTotal'), facesSkipped: sum('facesSkipped'), volumes: results.flatMap((x) => x[1].meta.volumes || []) }), fname, { format: 'solid-edge', container: 'Solid Edge compound file', kernelStream: results.map((x) => x[0]) });
}

const ASM = [0x41, 0x53, 0x4d, 0x20, 0x42, 0x69, 0x6e, 0x61, 0x72, 0x79, 0x46, 0x69, 0x6c, 0x65]; // "ASM BinaryFile"
const ACIS = [0x41, 0x43, 0x49, 0x53, 0x20, 0x42, 0x69, 0x6e, 0x61, 0x72, 0x79, 0x46, 0x69, 0x6c, 0x65]; // "ACIS BinaryFile"

/**
 * Autodesk Inventor part (.ipt).
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ deflection?: number, angularStep?: number }} [opts] as for ACIS SAB
 */
export async function parseInventor(input, opts = {}) {
  const bytes = inputBytes(input);
  const fname = input.name || 'IPT';
  if (!isCfb(bytes)) throw new Error(`${fname}: not an Inventor file (an Inventor document is a compound (OLE2) file). If this is another program's part file, export STEP from that program.`);
  const cfb = readCfb(bytes, fname);
  const segs = cfb.find(/^RSeStorage\/M[^/]+$/);
  if (!segs.length) throw new Error(`${fname}: this compound file has none of the segment streams of an Inventor document (RSeStorage/M…); it was written by another program or by a release before Inventor 11.`);
  // the name of a segment is UTF-16 text in its "M" stream
  const named = [];
  for (const e of segs) {
    const d = cfb.read(e.path);
    let t = '';
    for (let i = 0; i + 1 < d.length;) { if (d[i + 1] === 0 && d[i] >= 0x20 && d[i] < 0x7f) { t += String.fromCharCode(d[i]); i += 2; } else { t += ' '; i++; } }
    const m = /[A-Za-z]{2}[A-Za-z]*Segment\b/.exec(t);
    named.push({ id: e.name.slice(1), name: m ? m[0] : '' });
  }
  const brep = named.filter((x) => /BRepSegment$/i.test(x.name));
  if (!brep.length) {
    if (/\.(iam|ipn|idw)$/i.test(fname) || named.some((x) => /^Am/.test(x.name))) throw new Error(`${fname}: an Inventor assembly, presentation or drawing stores references to its part files and no bodies of its own. Import the part files (.ipt), or use Save Copy As → STEP or SAT in Inventor.`);
    throw new Error(`${fname}: the file has no B-rep segment (segments found: ${named.map((x) => x.name || '?').join(', ')}). Use Save Copy As → SAT or STEP in Inventor.`);
  }
  const results = [];
  let zstd = 0, empty = 0;
  for (const seg of brep) {
    const d = cfb.read(`RSeStorage/B${seg.id}`);
    if (!d || d.length < 20) continue;
    const method = d[17];
    if (d[16] === 4 && method === 2) { zstd++; continue; }
    let data;
    if (d[16] === 4 && method === 1) { try { data = await inflateZlib(d.subarray(18)); } catch { throw new Error(`${fname}: the B-rep segment could not be inflated — the file is damaged.`); } } else data = d;
    for (const sig of [ASM, ACIS]) for (let p = indexOf(data, sig, 0); p >= 0; p = indexOf(data, sig, p + sig.length)) {
      try { results.push(parseSab({ name: fname, bytes: data.subarray(p) }, opts)); } catch (e) {
        // a part that holds sketches or parameters only has an ACIS stream without bodies
        if (/no geometry could be read from the SAB file\.$/.test(e.message)) empty++; else throw e;
      }
    }
  }
  if (!results.length) {
    if (zstd) throw new Error(`${fname}: the segments of this Inventor file are compressed with Zstandard, as Inventor 2026 writes them, which this reader does not decode. Use Save Copy As → SAT (ACIS text) or STEP in Inventor; both are read.`);
    if (empty) throw new Error(`${fname}: the part has no solid or surface bodies — its B-rep segment is empty (a part of sketches, work geometry or parameters only).`);
    throw new Error(`${fname}: the B-rep segment holds no ACIS body data — the part has no solid or surface bodies, or keeps them in a form this reader does not know. Use Save Copy As → SAT or STEP in Inventor.`);
  }
  const extra = { format: 'inventor', container: 'Inventor compound file', kernelStream: brep.map((x) => `RSeStorage/B${x.id} (${x.name})`).join(', '), segments: named.map((x) => x.name).filter(Boolean) };
  if (results.length === 1) return tag(results[0], fname, extra);
  const models = [], warnings = [];
  for (const r of results) { models.push(...r.models); for (const w of r.warnings) if (!warnings.includes(w)) warnings.push(w); }
  const first = results[0];
  const sum = (k) => results.reduce((a, r) => a + (r.meta[k] || 0), 0);
  const elementTypes = {}, skipped = {};
  for (const r of results) { for (const k in r.meta.elementTypes) elementTypes[k] = (elementTypes[k] || 0) + r.meta.elementTypes[k]; for (const k in r.meta.skipped) skipped[k] = (skipped[k] || 0) + r.meta.skipped[k]; }
  return tag(finish(models, warnings, { ...first.meta, counts: undefined, bbox: undefined, elementTypes, skipped, bodies: sum('bodies'), solidBodies: sum('solidBodies'), closedBodies: sum('closedBodies'), orientedBodies: sum('orientedBodies'), faces: sum('faces'), faceTotal: sum('faceTotal'), volumes: results.flatMap((r) => r.meta.volumes || []) }), fname, extra);
}
