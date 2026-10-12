// Registry of the solver hand-offs offered to the user interface, and a ZIP packer so that a whole case folder
// downloads as one file.
//
// Every writer returns `{ files: [{ path, text | bytes }], notes: string[] }`; a path ending in '/' is an empty
// folder. Readers take the text of a result file (or a list of files) and return plain objects; see the JSDoc
// of each module and docs/interop.md for what was run against the real solvers.

import { writeNemohCase, readNemohResults, writeWamitGdf, writeCapytaineScript, readCapytaineCsv } from './nemoh.js';
import { writeCalculixCase, readCalculixFrd, readCalculixDat } from './calculix.js';
import { writeOpenfoamWaveTank, writeOpenfoamCurrentCase, readOpenfoamForces } from './openfoam.js';
import { hydroDbFromWamit } from '../engine/hydro/hydrodb.js';

export * from './nemoh.js';
export * from './calculix.js';
export * from './openfoam.js';

/**
 * Solver hand-offs. `write(model, options)` takes the geometry model as its first argument and the remaining
 * options of the underlying writer as the second; `read(input, options)` parses what the solver wrote
 * (`reads` names the files and the shape of `input`). `executed` tells whether a case written by this entry
 * was run through the real solver during development (see docs/interop.md).
 */
export const SOLVER_EXPORTS = Object.freeze([
  {
    id: 'nemoh', solver: 'NEMOH v3.0', purpose: 'Frequency-domain radiation/diffraction: added mass, wave damping and excitation of a floating hull',
    kinds: ['surface'], extension: 'zip', executed: false,
    write: (model, o = {}) => writeNemohCase({ ...o, surface: model, name: o.name || model.name }),
    read: (files, o = {}) => readNemohResults({ ...o, files }),
    reads: ['input: [{ path, text }]', 'results/RadiationCoefficients.tec (or CM.dat + CA.dat)', 'results/ExcitationForce.tec', 'Nemoh.cal (optional: rho, g, depth)'],
    notes: 'Closed hull, z up, waterline at z = 0. Options: rho, g, depth, omegas (rad/s), headingsDeg, cog. Follows the NEMOH v3 manual; not executed against NEMOH.',
  },
  {
    id: 'capytaine', solver: 'Capytaine 2.x (Python)', purpose: 'Frequency-domain radiation/diffraction through a ready-to-run Python script',
    kinds: ['surface'], extension: 'zip', executed: false,
    write: (model, o = {}) => writeCapytaineScript({ ...o, surface: model, name: o.name || model.name }),
    read: (text, o = {}) => readCapytaineCsv(text, o),
    reads: ['input: text of <name>_hydro.csv'],
    notes: 'Same options as NEMOH; frequencies and headings need not be uniformly spaced. The script also requests the infinite-frequency added mass. Not executed against Capytaine.',
  },
  {
    id: 'wamit-gdf', solver: 'WAMIT v7 (low-order GDF)', purpose: 'Panel geometry of the wetted hull for WAMIT and other codes that read GDF',
    kinds: ['surface'], extension: 'gdf', executed: false,
    write: (model, o = {}) => writeWamitGdf({ ...o, surface: model, name: o.name || model.name }),
    read: (input, o = {}) => hydroDbFromWamit({ ...o, text1: input.text1, text3: input.text3 }),
    reads: ['input: { text1, text3 }', '<run>.1 (added mass and damping) as text1', '<run>.3 (exciting forces) as text3; pass the same ULEN as option L'],
    notes: 'Options: ulen, g. Only the geometry is written; the WAMIT control files (.pot, .frc, .cfg) are the user\'s. Not executed against WAMIT.',
  },
  {
    id: 'calculix', solver: 'CalculiX 2.21', purpose: 'Linear and non-linear finite elements: static, large displacement, plasticity, buckling, natural frequencies, implicit dynamics',
    kinds: ['frame', 'surface', 'volume'], extension: 'inp', executed: false,
    write: (model, o = {}) => writeCalculixCase(model, o),
    read: (text, o = {}) => (/^\s{4}[12]C/m.test(text) ? readCalculixFrd(text, o) : readCalculixDat(text, o)),
    reads: ['input: text of <job>.frd (ASCII nodal results) or of <job>.dat (printed tables)'],
    notes: 'Options: material { E, nu, rho, fy, hardening }, supports, loads, analysis, steps, elementOrder, thickness. Follows the CalculiX manual; not executed against ccx.',
  },
  {
    id: 'openfoam-wave-tank', solver: 'OpenFOAM v2406 (interFoam)', purpose: 'Wave loads on a fixed structure in a numerical wave tank (volume of fluid)',
    kinds: ['surface'], extension: 'zip', executed: true,
    write: (model, o = {}) => writeOpenfoamWaveTank({ ...o, surface: model }),
    read: (text, o = {}) => readOpenfoamForces(text, o),
    reads: ['input: text of postProcessing/forces/0/force.dat', 'option momentText: text of postProcessing/forces/0/moment.dat'],
    notes: 'Options: domain { length, width, height, waterDepth }, wave { H, T, theory }, cellSize, refinement, endTime, turbulence. A current is not superposed on the waves.',
  },
  {
    id: 'openfoam-current', solver: 'OpenFOAM v2406 (simpleFoam)', purpose: 'Steady current loads on a fully submerged fixed body (RANS, k-omega SST)',
    kinds: ['surface'], extension: 'zip', executed: true,
    write: (model, o = {}) => writeOpenfoamCurrentCase({ ...o, surface: model }),
    read: (text, o = {}) => readOpenfoamForces(text, o),
    reads: ['input: text of postProcessing/forces/0/force.dat', 'option momentText: text of postProcessing/forces/0/moment.dat'],
    notes: 'Options: domain { length, width, height }, current (m/s), cellSize, refinement, iterations, turbulence. No free surface.',
  },
]);

/**
 * Writes a case for one of the registered solvers.
 * @param {string} id entry of SOLVER_EXPORTS
 * @param {object} model geometry model (CONTRACT §4)
 * @param {object} [options] options of the writer
 * @returns {{files:{path:string, text?:string, bytes?:Uint8Array}[], notes:string[], solver:string}}
 */
export function exportCase(id, model, options = {}) {
  const ex = SOLVER_EXPORTS.find((e) => e.id === id);
  if (!ex) throw new RangeError(`solver export: unknown id "${id}" (available: ${SOLVER_EXPORTS.map((e) => e.id).join(', ')}).`);
  if (!model || !ex.kinds.includes(model.kind)) throw new RangeError(`${ex.solver} export does not handle a "${model && model.kind}" model (supported kinds: ${ex.kinds.join(', ')}).`);
  return { ...ex.write(model, options), solver: ex.solver };
}

let CRC_TABLE = null;
/**
 * CRC-32 (IEEE 802.3, reflected polynomial 0xEDB88320) of a byte array, as stored in ZIP headers.
 * @param {Uint8Array} bytes
 * @returns {number} unsigned 32-bit checksum
 */
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
 * Packs the files of a case into a ZIP archive (method 0, stored; UTF-8 names; local headers, central directory
 * and end-of-central-directory record per the PKWARE APPNOTE 6.3). Paths ending in '/' become folder entries;
 * `Allrun`, `Allclean` and `.sh` / `.py` files are marked executable (Unix mode 0755) so that the scripts run
 * after extraction with unzip.
 * Validity: fewer than 65 535 entries and less than 4 GiB in total (no ZIP64).
 * @param {{path:string, text?:string, bytes?:Uint8Array}[]} files case files
 * @param {{root?:string, date?:Date}} [opts] root: folder name prefixed to every path; date: time stamp of the
 *   entries (default 1980-01-01, which keeps the archive reproducible)
 * @returns {Uint8Array} the archive
 */
export function packCase(files, { root = '', date } = {}) {
  if (!Array.isArray(files) || files.length === 0) throw new RangeError('packCase: files must be a non-empty list');
  if (files.length >= 0xffff) throw new RangeError('packCase: too many files for a ZIP archive without ZIP64');
  const enc = new TextEncoder();
  const prefix = root ? `${String(root).replace(/[\\/]+$/, '')}/` : '';
  const seen = new Set();
  const recs = files.map((f) => {
    const path = prefix + String(f.path || '').replace(/\\/g, '/').replace(/^\/+/, '');
    if (!f.path || path.split('/').includes('..')) throw new RangeError(`packCase: invalid path "${f.path}"`);
    if (seen.has(path)) throw new RangeError(`packCase: duplicate path "${path}"`);
    seen.add(path);
    const isDir = path.endsWith('/');
    const data = isDir ? new Uint8Array(0) : f.bytes ? f.bytes : enc.encode(f.text ?? '');
    const base = path.split('/').pop();
    const exec = /^(Allrun|Allclean)$/.test(base) || /\.(sh|py)$/.test(base);
    const mode = isDir ? 0o40755 : exec ? 0o100755 : 0o100644;
    return { name: enc.encode(path), data, crc: crc32(data), attr: ((mode << 16) | (isDir ? 0x10 : 0)) >>> 0 };
  });
  let size = 22;
  for (const r of recs) size += 30 + 46 + 2 * r.name.length + r.data.length;
  if (size > 0xffffffff) throw new RangeError('packCase: archive would exceed 4 GiB');
  const d = date || new Date(1980, 0, 1);
  const dosTime = (d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1);
  const dosDate = ((Math.max(1980, d.getFullYear()) - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate();
  const out = new Uint8Array(size);
  const dv = new DataView(out.buffer);
  let p = 0;
  for (const r of recs) {
    r.offset = p;
    dv.setUint32(p, 0x04034b50, true); dv.setUint16(p + 4, 20, true); dv.setUint16(p + 6, 0x0800, true); dv.setUint16(p + 8, 0, true);
    dv.setUint16(p + 10, dosTime, true); dv.setUint16(p + 12, dosDate, true);
    dv.setUint32(p + 14, r.crc, true); dv.setUint32(p + 18, r.data.length, true); dv.setUint32(p + 22, r.data.length, true);
    dv.setUint16(p + 26, r.name.length, true); dv.setUint16(p + 28, 0, true);
    out.set(r.name, p + 30); out.set(r.data, p + 30 + r.name.length);
    p += 30 + r.name.length + r.data.length;
  }
  const cen = p;
  for (const r of recs) {
    dv.setUint32(p, 0x02014b50, true);
    dv.setUint16(p + 4, (3 << 8) | 20, true); // made by: Unix, so that the external attributes carry the file mode
    dv.setUint16(p + 6, 20, true); dv.setUint16(p + 8, 0x0800, true); dv.setUint16(p + 10, 0, true);
    dv.setUint16(p + 12, dosTime, true); dv.setUint16(p + 14, dosDate, true);
    dv.setUint32(p + 16, r.crc, true); dv.setUint32(p + 20, r.data.length, true); dv.setUint32(p + 24, r.data.length, true);
    dv.setUint16(p + 28, r.name.length, true);
    dv.setUint32(p + 38, r.attr, true); dv.setUint32(p + 42, r.offset, true);
    out.set(r.name, p + 46);
    p += 46 + r.name.length;
  }
  dv.setUint32(p, 0x06054b50, true);
  dv.setUint16(p + 8, recs.length, true); dv.setUint16(p + 10, recs.length, true);
  dv.setUint32(p + 12, p - cen, true); dv.setUint32(p + 16, cen, true);
  return out;
}
