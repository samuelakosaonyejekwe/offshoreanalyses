// CAD boundary-representation import — STEP (ISO 10303-21, AP203 / AP214 / AP242), IGES 5.3 and OpenCascade BREP —
// through the vendored occt-import-js build of OpenCASCADE Technology (see docs/vendor.md).
//
// The kernel reads the exact geometry and tessellates every face; one `surface` model is returned per solid or
// shell, with one group per B-rep face. Lengths come out in metres: STEP and IGES declare their unit and the
// kernel converts; BREP files carry no unit, so millimetres (the OpenCascade convention) are assumed unless
// `opts.units` says otherwise.
//
// Tessellation is face by face, so vertices on shared edges are duplicated; `topology.inspect` and
// `repair` weld them. The chordal tolerance governs how well curved faces are followed — see `opts` of `parse`.
//
// The runtime is 7.7 MB of WebAssembly and is loaded on first use only, relative to this module.

import { inputBytes, baseName, cleanName, surfaceModel, finish, extOf } from '../importers/common.js';
import { resolveUnits, yUpToZUp } from '../units.js';

export const id = 'occt';
export const title = 'CAD B-rep (STEP, IGES, BREP)';
export const formats = ['step', 'iges', 'brep'];
export const extensions = ['step', 'stp', 'p21', 'iges', 'igs', 'brep', 'brp'];
export const library = { name: 'occt-import-js', version: '0.0.23', licence: 'LGPL-2.1' };
/** Download size of the runtime (JavaScript glue + WebAssembly) [bytes]. */
export const sizeBytes = 96871 + 7604031;

const BASE = new URL('../../../vendor/occt/', import.meta.url);
const isNode = typeof process !== 'undefined' && !!(process.versions && process.versions.node) && typeof window === 'undefined';

/** True when the runtime is part of this build. */
export function available() { return true; }

let loading = null;

async function loadFactory() {
  const src = new URL('occt-import-js.js', BASE);
  if (isNode) {
    const { createRequire } = await import('node:module');
    const { fileURLToPath } = await import('node:url');
    return createRequire(import.meta.url)(fileURLToPath(src));
  }
  const g = globalThis;
  if (typeof g.occtimportjs === 'function') return g.occtimportjs;
  if (typeof g.document !== 'undefined') {
    await new Promise((resolve, reject) => {
      const s = g.document.createElement('script');
      s.src = src.href;
      s.onload = resolve;
      s.onerror = () => reject(new Error('The CAD kernel script could not be loaded.'));
      g.document.head.append(s);
    });
  } else if (typeof g.importScripts === 'function') {
    try { g.importScripts(src.href); } catch { /* module workers have no importScripts: fall through */ }
  }
  if (typeof g.occtimportjs !== 'function') {
    // module worker: evaluate the classic script text in function scope
    const text = await (await fetch(src.href)).text();
    g.occtimportjs = new Function(`${text}\nreturn occtimportjs;`)();
  }
  return g.occtimportjs;
}

/**
 * Load and initialise the kernel (once; later calls return the same instance).
 * @returns {Promise<object>} the occt-import-js module (ReadStepFile, ReadIgesFile, ReadBrepFile)
 */
export function load() {
  if (!loading) {
    loading = (async () => {
      const factory = await loadFactory();
      const overrides = { locateFile: (f) => new URL(f, BASE).href };
      if (isNode) {
        const { readFile } = await import('node:fs/promises');
        overrides.wasmBinary = await readFile(new URL('occt-import-js.wasm', BASE));
      }
      return factory(overrides);
    })().catch((err) => { loading = null; throw new Error(`The CAD kernel could not be started: ${err && err.message ? err.message : err}. STEP, IGES and BREP import needs the files in vendor/occt/ and WebAssembly support.`); });
  }
  return loading;
}

/** Format of a CAD file from its content, falling back to its extension: 'step' | 'iges' | 'brep' | null. */
export function sniff(input) {
  const bytes = inputBytes(input);
  const head = String.fromCharCode.apply(null, bytes.subarray(0, Math.min(bytes.length, 400)));
  if (/ISO-10303-21/.test(head)) return 'step';
  if (/DBRep_DrawableShape|CASCADE Topology/.test(head)) return 'brep';
  const line = head.split(/\r?\n/)[0] || '';
  if (line.length >= 73 && /[SC]/.test(line[72])) return 'iges';
  const ext = extOf(input.name || '');
  if (['step', 'stp', 'p21'].includes(ext)) return 'step';
  if (['iges', 'igs'].includes(ext)) return 'iges';
  if (['brep', 'brp'].includes(ext)) return 'brep';
  return null;
}

const DEFLECTION = { coarse: [0.005, 0.8], normal: [0.001, 0.5], fine: [0.0002, 0.2] };

/**
 * Import a STEP, IGES or BREP file.
 * @param {{ name?: string, bytes?: Uint8Array, text?: string }} input file content
 * @param {object} [opts]
 * @param {'step'|'iges'|'brep'} [opts.format] default: detected from the content, then the extension
 * @param {'coarse'|'normal'|'fine'} [opts.quality] tessellation: chordal tolerance 0.5 % / 0.1 % / 0.02 % of the
 *   mean bounding-box size and angular tolerance 0.8 / 0.5 / 0.2 rad (default 'normal')
 * @param {number} [opts.linearDeflection] absolute chordal tolerance [m], overrides `quality`
 * @param {number} [opts.angularDeflection] angular tolerance [rad], overrides `quality`
 * @param {string} [opts.units] length unit of a BREP file (default 'mm'); ignored for STEP and IGES, which
 *   declare their unit
 * @param {'z'|'y'} [opts.upAxis] 'y' converts a y-up model to z up (default 'z')
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} surface models in metres;
 *   meta: { format, units, kernel, tessellation, assembly, counts, bbox }
 */
export async function parse(input, opts = {}) {
  const bytes = inputBytes(input);
  const format = opts.format || sniff(input);
  if (!format || !formats.includes(format)) throw new Error(`${input.name || 'input'}: not recognised as a STEP, IGES or BREP file.`);
  const occt = await load();
  const [ratio, ang] = DEFLECTION[opts.quality || 'normal'] || DEFLECTION.normal;
  const units = format === 'brep' ? resolveUnits(opts, null, 'mm') : { unit: 'm', scale: 1, basis: 'file' };
  const params = { linearUnit: 'meter', angularDeflection: opts.angularDeflection > 0 ? opts.angularDeflection : ang };
  if (opts.linearDeflection > 0) { params.linearDeflectionType = 'absolute_value'; params.linearDeflection = opts.linearDeflection / units.scale; }
  else { params.linearDeflectionType = 'bounding_box_ratio'; params.linearDeflection = ratio; }
  const read = format === 'step' ? occt.ReadStepFile : format === 'iges' ? occt.ReadIgesFile : occt.ReadBrepFile;
  let res;
  try { res = read.call(occt, bytes, params); } catch (err) {
    throw new Error(`${input.name || 'input'}: the CAD kernel failed while reading the file (${err && err.message ? err.message : 'out of memory or corrupt data'}).`);
  }
  if (!res || !res.success) throw new Error(`${input.name || 'input'}: the ${format.toUpperCase()} file could not be read — it is damaged, empty or uses entities the kernel does not support.`);
  const warnings = [];
  const fallback = baseName(input.name || format);
  const models = [];
  let empty = 0;
  (res.meshes || []).forEach((mesh, mi) => {
    const pos = mesh.attributes && mesh.attributes.position && mesh.attributes.position.array;
    const idx = mesh.index && mesh.index.array;
    if (!pos || !idx || !idx.length) { empty++; return; }
    const coords = Float64Array.from(pos);
    if (units.scale !== 1) for (let i = 0; i < coords.length; i++) coords[i] *= units.scale;
    if (opts.upAxis === 'y') yUpToZUp(coords);
    const groups = (mesh.brep_faces || []).map((f, fi) => ({ name: `face-${fi + 1}`, start: 3 * f.first, count: 3 * (f.last - f.first + 1), ...(f.color ? { color: f.color } : {}) }));
    const name = cleanName(mesh.name, res.meshes.length > 1 ? `${fallback}-${mi + 1}` : fallback);
    models.push(surfaceModel(name, coords, Uint32Array.from(idx), groups, mesh.color ? { color: mesh.color } : undefined));
  });
  if (empty) warnings.push(`${empty} shape(s) produced no triangles (wire or point geometry, or faces the kernel could not mesh) and were left out.`);
  if (!models.length) warnings.push('The file was read but holds no surfaces or solids.');
  if (format === 'brep' && units.basis === 'assumed') warnings.push('BREP files carry no length unit: millimetres were assumed. Set the unit if the model comes out 1000 times too small.');
  const tree = (node) => ({ name: node.name || '', meshes: node.meshes || [], children: (node.children || []).map(tree) });
  return finish(models, warnings, {
    format, units, kernel: `${library.name} ${library.version}`,
    tessellation: { linearDeflectionType: params.linearDeflectionType, linearDeflection: params.linearDeflection, angularDeflection: params.angularDeflection },
    assembly: res.root ? tree(res.root) : null,
  });
}
