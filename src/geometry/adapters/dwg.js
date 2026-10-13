// AutoCAD DWG — read through LibreDWG, which is NOT part of this application.
//
// The only open reader for DWG is LibreDWG (GNU GPL, version 3). To keep that licence out of this code base the
// library is neither vendored nor imported by any module of the application. Instead:
//   * dwg.worker.js, a module Web Worker that contains nothing but a message loop, imports the library at run
//     time from the address in LIBRARY_URL — the library's own distribution on a public CDN, at one exact
//     version — when, and only when, the caller confirms that the user has agreed (opts.consent === true);
//   * the DWG bytes go into the worker, and a plain description of the drawing comes out (header variables,
//     layer names, blocks, entities with their coordinates). Nothing else crosses the boundary;
//   * that description is turned into DXF group codes here and read by the application's own DXF reader
//     (importers/dxf.js), so a DWG drawing gives exactly the models its DXF export would give.
// Consequences: a network connection is needed the first time (about 11 MB; afterwards the browser cache may
// serve it), the worker must be allowed to load from the CDN, and it runs in a browser only. This describes the
// technical separation; it is not legal advice (docs/vendor.md).

import { inputBytes } from '../importers/common.js';
import { dxfModels } from '../importers/dxf.js';

/**
 * The owner's switch. Set to false and DWG reading is off for every user of this copy of the app: nothing is
 * downloaded, the worker is never started, and a DWG file is answered with "save it as DXF". Nothing else needs
 * to change.
 */
export const ENABLED = true;
const OFF = 'DWG reading is switched off in this copy of the app. Save the drawing as DXF in the originating program (in AutoCAD: SAVEAS → DXF) and import that.';

export const id = 'dwg';
export const title = 'AutoCAD DWG (through LibreDWG, downloaded on request)';
export const formats = ['dwg'];
export const extensions = ['dwg'];
/** Where the worker fetches the library from: its npm distribution on jsDelivr, at one exact version. */
export const LIBRARY_URL = 'https://cdn.jsdelivr.net/npm/@mlightcad/libredwg-web@0.7.17/dist/libredwg-web.js';
export const library = { name: '@mlightcad/libredwg-web (LibreDWG compiled to WebAssembly)', version: '0.7.17', licence: 'GPL-3.0', bundled: false, url: LIBRARY_URL };
/** What the browser downloads on first use (wrapper, glue code and WebAssembly) [bytes]. */
export const sizeBytes = 279544 + 110192 + 10816837;
export const requirements = { consent: true, network: 'first use', environment: 'browser (module worker)' };
/** Text to show the user before opts.consent is set. */
export const CONSENT = 'Reading a DWG drawing downloads the LibreDWG library (GNU GPL version 3, about 11 MB) from the public CDN cdn.jsdelivr.net into your browser and runs it there in a separate worker. The library is not part of this application. Confirm to continue (pass consent: true), or save the drawing as DXF instead.';

const FETCHED_KEY = 'offshore-analyses.dwg.library';
const store = () => { try { return typeof localStorage !== 'undefined' ? localStorage : null; } catch { return null; } };
/** True when this browser fetched this version of the library before (its cache may still hold it). */
export function fetchedBefore() { const s = store(); try { return !!s && s.getItem(FETCHED_KEY) === library.version; } catch { return false; } }

/**
 * Whether a DWG drawing can be read here and now: in a browser with module workers, and either on line or with
 * the library fetched before (the browser cache may then serve it; if it does not, parse() says so). The
 * user's consent is a separate matter: parse() refuses without opts.consent === true.
 */
export function available() {
  if (!ENABLED) return false;
  if (typeof Worker !== 'function' || typeof window === 'undefined') return false;
  const online = typeof navigator === 'undefined' || navigator.onLine !== false;
  return online || fetchedBefore();
}

/** Why available() is false, for the message of the registry. */
export function unavailableReason() {
  if (!ENABLED) return 'is switched off in this copy of the app';
  if (typeof Worker !== 'function' || typeof window === 'undefined') return 'runs in a web browser only (it loads the DWG library into a Web Worker at run time)';
  return 'needs a network connection the first time: the DWG library is downloaded on request and this device is off line';
}

const S = (v) => String(v);
const xyz = (out, base, p) => { if (p) out.push(10 + base, S(p.x ?? 0), 20 + base, S(p.y ?? 0), 30 + base, S(p.z ?? 0)); };
const extrusion = (out, n) => { if (n && (n.x || n.y || n.z !== 1) && Number.isFinite(n.x + n.y + n.z)) out.push(210, S(n.x), 220, S(n.y), 230, S(n.z)); };
const DEG = 180 / Math.PI;

/** Group codes of one entity of the drawing description (flat [code, value, …] starting at group 0), or null. */
function entityGroups(e) {
  const o = [0, e.type, 8, e.layer || '0'];
  if (e.isInPaperSpace) o.push(67, '1');
  switch (e.type) {
    case 'LINE': xyz(o, 0, e.startPoint); xyz(o, 1, e.endPoint); return o;
    case 'POINT': xyz(o, 0, e.position); return o;
    case 'CIRCLE': xyz(o, 0, e.center); o.push(40, S(e.radius)); extrusion(o, e.extrusionDirection); return o;
    // angles arrive in radians; DXF holds degrees
    case 'ARC': xyz(o, 0, e.center); o.push(40, S(e.radius), 50, S(e.startAngle * DEG), 51, S(e.endAngle * DEG)); extrusion(o, e.extrusionDirection); return o;
    case 'ELLIPSE': xyz(o, 0, e.center); xyz(o, 1, e.majorAxisEndPoint); o.push(40, S(e.axisRatio), 41, S(e.startAngle), 42, S(e.endAngle)); extrusion(o, e.extrusionDirection); return o;
    case '3DFACE': xyz(o, 0, e.corner1); xyz(o, 1, e.corner2); xyz(o, 2, e.corner3); xyz(o, 3, e.corner4 || e.corner3); return o;
    case 'SOLID': case 'TRACE': xyz(o, 0, e.corner1); xyz(o, 1, e.corner2); xyz(o, 2, e.corner3); xyz(o, 3, e.corner4 || e.corner3); extrusion(o, e.extrusionDirection); return o;
    case 'LWPOLYLINE':
      o.push(70, S(e.flag || 0), 38, S(e.elevation || 0));
      for (const v of e.vertices || []) { o.push(10, S(v.x), 20, S(v.y)); if (v.bulge) o.push(42, S(v.bulge)); }
      extrusion(o, e.extrusionDirection);
      return o;
    case 'POLYLINE2D': case 'POLYLINE3D': {
      o[1] = 'POLYLINE';
      o.push(70, S((e.flag || 0) | (e.type === 'POLYLINE3D' ? 8 : 0)), 30, S(e.elevation || 0));
      if (e.meshMVertexCount) o.push(71, S(e.meshMVertexCount), 72, S(e.meshNVertexCount || 0));
      extrusion(o, e.extrusionDirection);
      for (const v of e.vertices || []) {
        o.push(0, 'VERTEX', 8, e.layer || '0', 10, S(v.x ?? 0), 20, S(v.y ?? 0), 30, S(v.z ?? 0), 70, S(v.flag || 0));
        if (v.bulge) o.push(42, S(v.bulge));
        ['polyfaceIndex0', 'polyfaceIndex1', 'polyfaceIndex2', 'polyfaceIndex3'].forEach((k, i) => { if (v[k]) o.push(71 + i, S(v[k])); });
      }
      o.push(0, 'SEQEND');
      return o;
    }
    case 'SPLINE':
      o.push(70, S(e.flag || 0), 71, S(e.degree || 3));
      for (const k of e.knots || []) o.push(40, S(k));
      if (e.weights && e.weights.length === (e.controlPoints || []).length) for (const w of e.weights) o.push(41, S(w));
      for (const p of e.controlPoints || []) xyz(o, 0, p);
      for (const p of e.fitPoints || []) xyz(o, 1, p);
      return o;
    case 'INSERT':
      o.push(2, e.name || '');
      xyz(o, 0, e.insertionPoint);
      o.push(41, S(e.xScale ?? 1), 42, S(e.yScale ?? 1), 43, S(e.zScale ?? 1), 50, S((e.rotation || 0) * DEG), 70, S(e.columnCount || 1), 71, S(e.rowCount || 1), 44, S(e.columnSpacing || 0), 45, S(e.rowSpacing || 0));
      extrusion(o, e.extrusionDirection);
      return o;
    default:
      // solids and surfaces: the library passes on the ACIS text when the drawing stores it as text
      if (typeof e.data === 'string' && e.data.length) { for (const line of e.data.split(/\r?\n/)) if (line.length) o.push(1, line); }
      return o;
  }
}

/**
 * Turn the drawing description the worker returns into the group-code sections of the DXF reader.
 * @param {{ header: object, layers: string[], blocks: object[], entities: object[] }} drawing
 * @returns {{ HEADER: Array, TABLES: Array, BLOCKS: Array, ENTITIES: Array }}
 */
export function drawingToDxfSections(drawing) {
  const sec = { HEADER: [], TABLES: [], BLOCKS: [], ENTITIES: [] };
  const h = drawing.header || {};
  if (Number.isFinite(h.INSUNITS)) sec.HEADER.push(9, '$INSUNITS', 70, S(h.INSUNITS));
  if (typeof h.ACADVER === 'string') sec.HEADER.push(9, '$ACADVER', 1, h.ACADVER);
  for (const name of drawing.layers || []) sec.TABLES.push(0, 'LAYER', 2, S(name));
  for (const b of drawing.blocks || []) {
    // model and paper space are block records too: their entities are the drawing itself
    if (/^\*(model|paper)_space/i.test(b.name || '')) continue;
    sec.BLOCKS.push(0, 'BLOCK', 2, S(b.name));
    const bp = b.basePoint || {};
    sec.BLOCKS.push(10, S(bp.x ?? 0), 20, S(bp.y ?? 0), 30, S(bp.z ?? 0));
    for (const e of b.entities || []) sec.BLOCKS.push(...entityGroups(e));
    sec.BLOCKS.push(0, 'ENDBLK');
  }
  for (const e of drawing.entities || []) sec.ENTITIES.push(...entityGroups(e));
  return sec;
}

let serial = 0;

/**
 * Ask the worker to read a drawing.
 * @param {Uint8Array} bytes @param {{ worker?: () => object, timeoutMs?: number }} opts `worker`: factory of the
 *   worker (tests supply a stand-in)
 * @returns {Promise<object>} the drawing description
 */
function readInWorker(bytes, opts) {
  const make = opts.worker || (() => new Worker(new URL('./dwg.worker.js', import.meta.url), { type: 'module' }));
  return new Promise((resolve, reject) => {
    let w;
    try { w = make(); } catch (e) { reject(new Error(`The DWG worker could not be started (${e.message}).`)); return; }
    const idNo = ++serial;
    const done = (fn, v) => { clearTimeout(timer); try { w.terminate(); } catch { /* already gone */ } fn(v); };
    const timer = setTimeout(() => done(reject, new Error('The DWG library did not answer in time. Check the network connection (it is downloaded from cdn.jsdelivr.net on first use) and try again, or save the drawing as DXF.')), opts.timeoutMs > 0 ? opts.timeoutMs : 180000);
    w.onmessage = (ev) => {
      const m = ev.data;
      if (!m || m.id !== idNo) return;
      if (m.ok) done(resolve, m.drawing);
      else done(reject, new Error(/import|fetch|load|network|module/i.test(m.error) ? `The DWG library could not be loaded from ${LIBRARY_URL} (${m.error}). A network connection is needed the first time; alternatively save the drawing as DXF.` : `The DWG library could not read the drawing (${m.error}). Save it as DXF in the originating program.`));
    };
    w.onerror = (ev) => done(reject, new Error(`The DWG worker failed (${(ev && ev.message) || 'the library could not be loaded'}). A network connection to cdn.jsdelivr.net is needed the first time; alternatively save the drawing as DXF.`));
    const buf = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    w.postMessage({ id: idNo, cmd: 'read', consent: true, url: opts.libraryUrl || LIBRARY_URL, bytes: buf }, [buf]);
  });
}

/**
 * Read a DWG drawing.
 * @param {{ name?: string, bytes: Uint8Array }} input
 * @param {{ consent?: boolean, worker?: Function, timeoutMs?: number }} [opts] `consent` must be true: the user has
 *   agreed to download and run the GPL-3.0-licensed library; all options of the DXF reader apply as well
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>}
 */
export async function parse(input, opts = {}) {
  const fname = input.name || 'DWG';
  if (!ENABLED) throw new Error(`${fname}: ${OFF}`);
  if (opts.consent !== true) throw new Error(`${fname}: ${CONSENT}`);
  if (!opts.worker && !available()) throw new Error(`${fname}: DWG reading ${unavailableReason()}. Save the drawing as DXF in the originating program instead.`);
  const bytes = inputBytes(input);
  if (bytes.length < 6 || bytes[0] !== 0x41 || bytes[1] !== 0x43) throw new Error(`${fname}: not a DWG drawing (the "AC10xx" version string at the start is missing).`);
  const drawing = await readInWorker(bytes, opts);
  const s = store();
  try { if (s) s.setItem(FETCHED_KEY, library.version); } catch { /* storage disabled */ }
  if (!drawing || !Array.isArray(drawing.entities)) throw new Error(`${fname}: the DWG library returned no drawing.`);
  const sec = drawingToDxfSections(drawing);
  if (!sec.ENTITIES.length && !sec.BLOCKS.length) throw new Error(`${fname}: the drawing holds no entities.`);
  const r = dxfModels(sec, input, opts);
  r.meta.format = 'dwg';
  r.meta.dwgVersion = String.fromCharCode(...bytes.subarray(0, 6));
  r.meta.library = { name: library.name, version: library.version, licence: library.licence, loadedFrom: opts.libraryUrl || LIBRARY_URL };
  r.warnings.push(`The drawing was read with ${library.name} ${library.version} (${library.licence}), downloaded into the browser at your request; it is not part of this application.`);
  return r;
}
