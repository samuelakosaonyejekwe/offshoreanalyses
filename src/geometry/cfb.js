// Compound File Binary reader ("OLE2 structured storage", Microsoft [MS-CFB]): the container of SolidWorks
// files up to 2014, Solid Edge and Inventor files. Read-only, synchronous, no dependencies.
//
// Layout: a 512-byte header, then sectors of 512 bytes (version 3) or 4096 bytes (version 4). The FAT chains
// sectors into streams; its own sectors are listed in the DIFAT (109 entries in the header, more in chained
// DIFAT sectors). The directory is a stream of 128-byte entries forming a tree (child / left / right
// siblings). Streams shorter than the cut-off (4096 bytes) live in 64-byte mini sectors inside the "mini
// stream" owned by the root entry, chained by the mini FAT.
//
// Every chain is followed with a bound on its length, so a damaged or hostile file ends in an Error with the
// reason, never in a loop.

const SIGNATURE = [0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1];
const END_OF_CHAIN = 0xfffffffe, FREE = 0xffffffff, MAX_REGULAR = 0xfffffffa;
const NO_STREAM = 0xffffffff;

/** Whether the bytes start with the compound-file signature. */
export function isCfb(bytes) {
  if (!bytes || bytes.length < 512) return false;
  for (let i = 0; i < 8; i++) if (bytes[i] !== SIGNATURE[i]) return false;
  return true;
}

/**
 * Open a compound file.
 * @param {Uint8Array} bytes
 * @param {string} [fname] name used in error messages
 * @returns {{ entries: { path: string, name: string, size: number, storage: boolean }[],
 *   read: (path: string) => Uint8Array|null, find: (test: RegExp) => object[] }} `path` joins the names of the
 *   storages with "/" (the root entry left out); `read` is case-sensitive on the full path
 */
export function readCfb(bytes, fname = 'file') {
  const fail = (why) => new Error(`${fname}: ${why}`);
  if (!isCfb(bytes)) throw fail('not a compound (OLE2) file: the signature D0 CF 11 E0 A1 B1 1A E1 is missing or the file is shorter than its 512-byte header.');
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const u16 = (o) => dv.getUint16(o, true), u32 = (o) => dv.getUint32(o, true);
  const major = u16(0x1a), shift = u16(0x1e), miniShift = u16(0x20);
  if (u16(0x1c) !== 0xfffe) throw fail('the compound-file header has an unknown byte-order mark.');
  if (!((major === 3 && shift === 9) || (major === 4 && shift === 12)) || miniShift !== 6) throw fail(`the compound-file header declares version ${major} with sector size 2^${shift} and mini-sector size 2^${miniShift}; only version 3 (512 bytes) and version 4 (4096 bytes) with 64-byte mini sectors exist.`);
  const size = 1 << shift, per = size >> 2;
  const total = Math.ceil(bytes.length / size);
  const offsetOf = (sec) => {
    const o = (sec + 1) * size;
    if (sec > MAX_REGULAR || o >= bytes.length) throw fail(`a sector number (${sec}) points beyond the end of the file — the file is truncated or damaged.`);
    return o;
  };
  const nFat = u32(0x2c), dirStart = u32(0x30), cutoff = u32(0x38), miniFatStart = u32(0x3c), nMiniFat = u32(0x40);
  let difatSec = u32(0x44);
  const nDifat = u32(0x48);
  if (nFat > total + 1 || nDifat > total + 1) throw fail('the compound-file header declares more allocation-table sectors than the file has sectors — the file is truncated or damaged.');

  // ---- FAT ----
  const fatSectors = [];
  for (let i = 0; i < 109 && fatSectors.length < nFat; i++) { const s = u32(0x4c + 4 * i); if (s <= MAX_REGULAR) fatSectors.push(s); }
  for (let k = 0; fatSectors.length < nFat && difatSec <= MAX_REGULAR; k++) {
    if (k > nDifat + 1) throw fail('the chain of DIFAT sectors does not end — the file is damaged.');
    const o = offsetOf(difatSec);
    for (let i = 0; i < per - 1 && fatSectors.length < nFat; i++) { if (o + 4 * i + 4 > bytes.length) break; const s = u32(o + 4 * i); if (s <= MAX_REGULAR) fatSectors.push(s); }
    if (o + size > bytes.length) break;
    difatSec = u32(o + size - 4);
  }
  const fat = new Uint32Array(fatSectors.length * per).fill(FREE);
  fatSectors.forEach((s, k) => { const o = offsetOf(s); for (let i = 0; i < per && o + 4 * i + 4 <= bytes.length; i++) fat[k * per + i] = u32(o + 4 * i); });

  /** Bytes of the sector chain starting at `start`; `len` < 0 reads to the end of the chain. */
  const chain = (start, len, what) => {
    const secs = [];
    for (let s = start; s !== END_OF_CHAIN && s !== FREE; s = fat[s]) {
      if (s > MAX_REGULAR || s >= fat.length) throw fail(`the sector chain of ${what} leaves the allocation table — the file is truncated or damaged.`);
      secs.push(s);
      if (secs.length > total + 1) throw fail(`the sector chain of ${what} runs in a loop — the file is damaged.`);
      if (len >= 0 && secs.length * size >= len) break;
    }
    const n = len >= 0 ? len : secs.length * size;
    if (secs.length * size < n) throw fail(`${what} is declared with ${n} bytes but its sector chain holds ${secs.length * size} — the file is truncated or damaged.`);
    const out = new Uint8Array(n);
    for (let i = 0, p = 0; p < n; i++, p += size) {
      const o = offsetOf(secs[i]);
      const k = Math.min(size, n - p);
      if (o + k > bytes.length) throw fail(`${what} runs beyond the end of the file — the file is truncated.`);
      out.set(bytes.subarray(o, o + k), p);
    }
    return out;
  };

  // ---- directory ----
  const dir = chain(dirStart, -1, 'the directory');
  const ddv = new DataView(dir.buffer);
  const count = Math.floor(dir.length / 128);
  const raw = [];
  for (let i = 0; i < count; i++) {
    const o = 128 * i;
    const type = dir[o + 0x42];
    let nlen = ddv.getUint16(o + 0x40, true);
    if (nlen > 64) nlen = 64;
    let name = '';
    for (let k = 0; k + 2 < nlen; k += 2) name += String.fromCharCode(ddv.getUint16(o + k, true));
    const lo = ddv.getUint32(o + 0x78, true), hi = major === 4 ? ddv.getUint32(o + 0x7c, true) : 0;
    raw.push({ name, type, left: ddv.getUint32(o + 0x44, true), right: ddv.getUint32(o + 0x48, true), child: ddv.getUint32(o + 0x4c, true), start: ddv.getUint32(o + 0x74, true), size: hi * 4294967296 + lo });
  }
  if (!raw.length || raw[0].type !== 5) throw fail('the directory of the compound file has no root entry — the file is damaged.');
  const entries = [];
  const byPath = new Map();
  {
    const seen = new Uint8Array(count);
    // siblings form a binary tree; walk it with an explicit stack
    const stack = raw[0].child !== NO_STREAM ? [[raw[0].child, '']] : [];
    while (stack.length) {
      const [id, prefix] = stack.pop();
      if (id === NO_STREAM || id >= count || seen[id]) continue;
      seen[id] = 1;
      const e = raw[id];
      if (e.left !== NO_STREAM) stack.push([e.left, prefix]);
      if (e.right !== NO_STREAM) stack.push([e.right, prefix]);
      if (e.type !== 1 && e.type !== 2) continue;
      const path = prefix + e.name;
      const entry = { path, name: e.name, size: e.type === 2 ? e.size : 0, storage: e.type === 1 };
      entries.push(entry);
      if (e.type === 2) byPath.set(path, e);
      else if (e.child !== NO_STREAM) stack.push([e.child, `${path}/`]);
    }
  }
  entries.sort((a, b) => (a.path < b.path ? -1 : a.path > b.path ? 1 : 0));

  // ---- mini stream ----
  let mini = null, miniFat = null;
  const loadMini = () => {
    if (mini) return;
    mini = raw[0].start <= MAX_REGULAR ? chain(raw[0].start, -1, 'the mini stream') : new Uint8Array(0);
    const mf = miniFatStart <= MAX_REGULAR && nMiniFat ? chain(miniFatStart, -1, 'the mini FAT') : new Uint8Array(0);
    miniFat = new Uint32Array(mf.length >> 2);
    const mdv = new DataView(mf.buffer);
    for (let i = 0; i < miniFat.length; i++) miniFat[i] = mdv.getUint32(4 * i, true);
  };
  const readEntry = (e, path) => {
    if (e.size > bytes.length) throw fail(`the stream "${path}" is declared with ${e.size} bytes, more than the file holds — the file is truncated or damaged.`);
    if (e.size === 0) return new Uint8Array(0);
    if (e.size >= cutoff) return chain(e.start, e.size, `the stream "${path}"`);
    loadMini();
    const out = new Uint8Array(e.size);
    let s = e.start, p = 0;
    for (let k = 0; p < e.size; k++) {
      if (s >= miniFat.length || k > miniFat.length) throw fail(`the mini-sector chain of the stream "${path}" is broken — the file is truncated or damaged.`);
      const o = s * 64, n = Math.min(64, e.size - p);
      if (o + n > mini.length) throw fail(`the stream "${path}" runs beyond the mini stream — the file is truncated or damaged.`);
      out.set(mini.subarray(o, o + n), p);
      p += n;
      s = miniFat[s];
    }
    return out;
  };
  return {
    entries,
    read(path) { const e = byPath.get(path); return e ? readEntry(e, path) : null; },
    find(test) { return entries.filter((e) => !e.storage && test.test(e.path)); },
  };
}

/**
 * Which program's document a compound file is, judged by its stream names: 'solidworks', 'solid-edge',
 * 'inventor', 'revit', or '' when none of these (or when the directory cannot be read).
 */
export function cfbProgram(bytes) {
  let names;
  try { names = readCfb(bytes).entries.map((e) => e.path); } catch { return ''; }
  if (names.some((n) => /^Contents\/|^Header2$|^_MO_VERSION_/.test(n))) return 'solidworks';
  if (names.some((n) => /^PARASOLID\/|^PSMcluster|^JSitesList$/.test(n))) return 'solid-edge';
  if (names.some((n) => /^RSeStorage\/|^RSeDb$|^UFRxDoc$/.test(n))) return 'inventor';
  if (names.some((n) => /^BasicFileInfo$|^Formats\/|^Global\/|^Partitions\//.test(n))) return 'revit';
  return '';
}
