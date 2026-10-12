// IFC (ISO 16739: IFC2x3, IFC4, IFC4x3) full product geometry through the vendored web-ifc engine
// (see docs/vendor.md).
//
// The engine evaluates every representation it supports (swept and extruded solids, B-reps, boolean clipping,
// openings, mapped items, tessellations) into triangles. One `surface` model is returned per IFC class
// (IfcBeam, IfcWall, IfcSlab …) with one group per element, so a whole discipline can be shown, hidden or
// measured at once and a single element can still be picked. Groups carry the element name, GlobalId and
// express id; `meta.elements` lists them all.
//
// Units and axes: the engine places everything with a 4 × 4 matrix that already contains the project length
// unit (IfcUnitAssignment, SI prefixes and conversion-based units) and a rotation to its internal y-up frame.
// The matrices are applied here and the result is rotated back to z-up, which is the IFC convention, so
// coordinates come out as in the file, in metres. `opts.units` overrides the unit declared in the file.
//
// Linear elements with a profile and a material are better imported without this adapter (importers/ifc.js
// gives a frame model with section properties); this one is for plates, slabs, walls, equipment and anything
// else that has to be seen or meshed as a surface.
//
// The runtime (JavaScript API with the IFC schema tables, plus WebAssembly) is about 7.5 MB and is loaded on
// first use only, relative to this module.

import { inputBytes, baseName, cleanName, surfaceModel, finish, Grow } from '../importers/common.js';
import { resolveUnits, normaliseUnit, UNIT_SCALE } from '../units.js';

export const id = 'web-ifc';
export const title = 'IFC full geometry (IFC2x3 / IFC4 / IFC4x3)';
export const formats = ['ifc'];
export const extensions = ['ifc'];
export const library = { name: 'web-ifc', version: '0.0.78', licence: 'MPL-2.0' };
/** Download size of the browser runtime (JavaScript API + WebAssembly) [bytes]. */
export const sizeBytes = 5909050 + 1595268;
/** What to convert a file to so that it can be imported when this adapter cannot run. */
export const convertTo = 'STEP, glTF or OBJ exported from the BIM tool';

const BASE = new URL('../../../vendor/web-ifc/', import.meta.url);
const NODE_BASE = new URL('node/', BASE);
const isNode = typeof process !== 'undefined' && !!(process.versions && process.versions.node) && typeof window === 'undefined';

function filesPresent(urls) {
  if (!isNode || typeof process.getBuiltinModule !== 'function') return true;
  try { const fs = process.getBuiltinModule('node:fs'); return urls.every((u) => fs.existsSync(u)); } catch { return true; }
}

/**
 * True when the engine can run here: WebAssembly is supported and, under Node (where it can be checked
 * without a request), the vendored files are on disk.
 */
export function available() {
  if (typeof WebAssembly !== 'object' || typeof WebAssembly.instantiate !== 'function') return false;
  return filesPresent(isNode ? [new URL('web-ifc-api-node.js', NODE_BASE), new URL('web-ifc-node.wasm', NODE_BASE)] : []);
}

let loading = null;

/**
 * Load and initialise the engine (once; later calls return the same instance).
 * @returns {Promise<{ api: object, lib: object }>} the IfcAPI instance and the web-ifc module (type codes)
 */
export function load() {
  if (!loading) {
    loading = (async () => {
      let lib, locate;
      if (isNode) {
        const { createRequire } = await import('node:module');
        const { fileURLToPath } = await import('node:url');
        lib = createRequire(import.meta.url)(fileURLToPath(new URL('web-ifc-api-node.js', NODE_BASE)));
        locate = (file) => fileURLToPath(new URL(file, NODE_BASE));
      } else {
        lib = await import(new URL('web-ifc-api.js', BASE).href);
        locate = (file) => new URL(file, BASE).href;
      }
      const api = new lib.IfcAPI();
      // single-threaded build: the threaded one needs cross-origin isolation and a second WebAssembly file
      await api.Init(locate, true);
      if (!api.wasmModule) throw new Error('the WebAssembly module did not initialise');
      return { api, lib };
    })().catch((err) => { loading = null; throw new Error(`The IFC geometry engine could not be started: ${err && err.message ? err.message : err}. Full IFC geometry needs the files in vendor/web-ifc/ and WebAssembly support; without it, import the file without fullGeometry (beams, columns, members and piles are read directly) or convert it to ${convertTo}.`); });
  }
  return loading;
}

/** True when the content is an IFC STEP physical file (ISO-10303-21 with an IFC schema). */
export function sniff(input) {
  const bytes = inputBytes(input);
  const head = String.fromCharCode.apply(null, bytes.subarray(0, Math.min(bytes.length, 2000)));
  return /ISO-10303-21/.test(head) && /FILE_SCHEMA\s*\(\s*\(\s*'IFC/i.test(head);
}

const SI_PREFIX = { EXA: 1e18, PETA: 1e15, TERA: 1e12, GIGA: 1e9, MEGA: 1e6, KILO: 1e3, HECTO: 100, DECA: 10, DECI: 0.1, CENTI: 0.01, MILLI: 1e-3, MICRO: 1e-6, NANO: 1e-9 };
const val = (v) => (v && typeof v === 'object' ? (v.value ?? v._representationValue ?? v._internalValue) : v);

/** Length unit declared by IfcUnitAssignment: `{ name, scale }` (metres per unit) or null. */
function declaredLengthUnit(api, lib, model) {
  let ids;
  try { ids = api.GetLineIDsWithType(model, lib.IFCUNITASSIGNMENT); } catch { return null; }
  const si = (u) => (u && val(u.UnitType) === 'LENGTHUNIT' && /METRE|METER/.test(String(val(u.Name))) ? (SI_PREFIX[String(val(u.Prefix) || '')] || 1) : null);
  for (let i = 0; i < ids.size(); i++) {
    let line;
    try { line = api.GetLine(model, ids.get(i), true); } catch { continue; }
    for (const u of (line && line.Units) || []) {
      if (!u || val(u.UnitType) !== 'LENGTHUNIT') continue;
      const s = si(u);
      if (s !== null && !u.ConversionFactor) return { name: { 1: 'm', 0.001: 'mm', 0.01: 'cm', 1000: 'km' }[s] || 'custom', scale: s };
      const cf = u.ConversionFactor;
      if (cf) {
        const f = Number(val(cf.ValueComponent)), b = si(cf.UnitComponent);
        const nm = String(val(u.Name) || '');
        const known = normaliseUnit(nm.toLowerCase());
        if (f > 0 && b !== null) return { name: known || nm.toLowerCase() || 'custom', scale: f * b };
        if (known) return { name: known, scale: UNIT_SCALE[known] };
      }
    }
  }
  return null;
}

/**
 * Import the tessellated geometry of every product of an IFC file.
 * @param {{ name?: string, bytes?: Uint8Array, text?: string }} input file content (STEP physical file)
 * @param {object} [opts]
 * @param {string} [opts.units] length unit, overriding the one in IfcUnitAssignment
 * @param {number} [opts.circleSegments] segments per full circle for curved profiles (default 12, 4–128)
 * @param {string[]} [opts.classes] keep only these IFC classes (case-insensitive, e.g. ['IfcPlate', 'IfcSlab'])
 * @param {boolean} [opts.openings] also return IfcOpeningElement and IfcSpace volumes (default false: they
 *   are voids and rooms, not structure)
 * @returns {Promise<{ models: object[], warnings: string[], meta: object }>} one surface per IFC class, metres,
 *   z up; meta: { format, schema, units, engine, classes: { IfcBeam: n, … }, elements: [{ class, name,
 *   globalId, expressId, model, group, triangles }], counts, bbox }
 */
export async function parse(input, opts = {}) {
  const fname = (input && input.name) || 'IFC';
  const bytes = inputBytes(input || {});
  if (!sniff({ bytes })) {
    const zip = bytes.length > 3 && bytes[0] === 0x50 && bytes[1] === 0x4b;
    const xml = bytes.length > 1 && bytes[0] === 0x3c;
    throw new Error(`${fname}: not an IFC STEP file (expected "ISO-10303-21" and FILE_SCHEMA(('IFC…')) in the header).${zip ? ' This is a ZIP archive: unpack the .ifc inside an .ifczip first.' : xml ? ' ifcXML is not read: export IFC-SPF (.ifc) from the BIM tool.' : ''}`);
  }
  const { api, lib } = await load();
  const segs = Number.isFinite(opts.circleSegments) ? Math.min(128, Math.max(4, Math.round(opts.circleSegments))) : 12;
  let model = -1;
  try { model = api.OpenModel(bytes, { COORDINATE_TO_ORIGIN: false, CIRCLE_SEGMENTS: segs }); } catch (err) {
    throw new Error(`${fname}: the IFC file could not be opened (${err && err.message ? err.message : 'corrupt data or out of memory'}). Check that the file is complete — it must end with END-ISO-10303-21;.`);
  }
  if (!(model >= 0)) throw new Error(`${fname}: the IFC file could not be opened — it is damaged or its schema is not one of IFC2x3, IFC4, IFC4x3.`);
  try {
    const warnings = [];
    let schema = null;
    try { schema = api.GetModelSchema(model) || null; } catch { /* reported as unknown */ }
    const lu = declaredLengthUnit(api, lib, model);
    let units;
    const override = !!(opts.units || (Number.isFinite(opts.scale) && opts.scale > 0));
    if (override) units = resolveUnits(opts, null, 'm');
    else if (lu) units = { unit: normaliseUnit(lu.name) || lu.name, scale: lu.scale, basis: 'file' };
    else { units = { unit: 'm', scale: 1, basis: 'assumed' }; warnings.push('No length unit was found in IfcUnitAssignment; metres assumed.'); }
    const want = Array.isArray(opts.classes) && opts.classes.length ? new Set(opts.classes.map((c) => String(c).toUpperCase())) : null;
    const VOIDS = new Set(['IFCOPENINGELEMENT', 'IFCOPENINGSTANDARDCASE', 'IFCSPACE', 'IFCVIRTUALELEMENT']);

    const byClass = new Map();
    const elements = [];
    let empty = 0, hidden = 0, failed = 0;
    const p = [0, 0, 0];
    api.StreamAllMeshes(model, (mesh) => {
      const eid = mesh.expressID;
      let cls = 'IfcProduct';
      try { cls = api.GetNameFromTypeCode(api.GetLineType(model, eid)) || cls; } catch { /* keep the generic name */ }
      const key = cls.toUpperCase();
      if ((want && !want.has(key)) || (!want && !opts.openings && VOIDS.has(key))) { hidden++; return; }
      let b = byClass.get(cls);
      if (!b) { b = { coords: new Grow(Float64Array, 3 * 1024), idx: new Grow(Uint32Array, 3 * 1024), groups: [] }; byClass.set(cls, b); }
      const start = b.idx.n;
      let colour = null;
      const n = mesh.geometries.size();
      for (let gi = 0; gi < n; gi++) {
        const pg = mesh.geometries.get(gi);
        let g;
        try { g = api.GetGeometry(model, pg.geometryExpressID); } catch { failed++; continue; }
        // views into the engine's memory: consumed before the next call into it
        const v = api.GetVertexArray(g.GetVertexData(), g.GetVertexDataSize());
        const ix = api.GetIndexArray(g.GetIndexData(), g.GetIndexDataSize());
        const m = pg.flatTransformation;
        const nv = v.length / 6;
        if (nv && ix.length) {
          const base = b.coords.n / 3;
          b.coords.ensure(3 * nv);
          const out = b.coords.a;
          let o = b.coords.n;
          for (let k = 0; k < nv; k++) {
            const x = v[6 * k], y = v[6 * k + 1], z = v[6 * k + 2];
            p[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
            p[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
            p[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
            // engine frame is y-up; IFC (and this application) is z-up
            out[o++] = p[0]; out[o++] = -p[2]; out[o++] = p[1];
          }
          b.coords.n = o;
          b.idx.ensure(ix.length);
          const io = b.idx.a;
          let q = b.idx.n;
          for (let k = 0; k < ix.length; k++) io[q++] = base + ix[k];
          b.idx.n = q;
          if (!colour && pg.color) colour = [pg.color.x, pg.color.y, pg.color.z, pg.color.w];
        }
        if (g && typeof g.delete === 'function') g.delete();
      }
      const count = b.idx.n - start;
      if (!count) { empty++; return; }
      let name = '', globalId = null, tag = null;
      try {
        const line = api.GetLine(model, eid);
        name = String(val(line.Name) ?? '');
        globalId = line.GlobalId ? String(val(line.GlobalId)) : null;
        tag = line.Tag ? String(val(line.Tag)) : null;
        if (typeof api.DecodeText === 'function' && /\\X2?\\|\\S\\/.test(name)) name = api.DecodeText(name);
      } catch { /* an element without readable attributes still has geometry */ }
      const grp = { name: cleanName(name, `#${eid}`), start, count, expressId: eid };
      if (globalId) grp.globalId = globalId;
      if (tag) grp.tag = tag;
      if (colour) grp.color = colour;
      b.groups.push(grp);
    }, !override);

    const fallback = baseName(fname);
    const models = [];
    const classes = {};
    for (const [cls, b] of byClass) {
      if (!b.idx.n) continue;
      const coords = b.coords.view();
      if (override && units.scale !== 1) for (let i = 0; i < coords.length; i++) coords[i] *= units.scale;
      const mi = models.length;
      models.push(surfaceModel(cls, coords, b.idx.finish(), b.groups, { ifcClass: cls, source: fallback }));
      classes[cls] = b.groups.length;
      b.groups.forEach((g, gi) => elements.push({ class: cls, name: g.name, globalId: g.globalId || null, expressId: g.expressId, model: mi, group: gi, triangles: g.count / 3 }));
    }
    if (failed) warnings.push(`${failed} geometry item(s) could not be evaluated by the engine and were left out.`);
    if (empty) warnings.push(`${empty} product(s) produced no triangles (curve or point representations, or geometry the engine does not support) and were left out.`);
    if (hidden && !want) warnings.push(`${hidden} opening / space volume(s) were left out (pass openings = true to keep them).`);
    if (!models.length) warnings.push(want ? `No product of the requested classes (${opts.classes.join(', ')}) has surface geometry.` : 'The file was read but none of its products has surface geometry. Linear elements defined by an axis only are read by importing without fullGeometry.');
    return finish(models, warnings, { format: 'ifc', schema, units, engine: `${library.name} ${library.version}`, circleSegments: segs, classes, elements });
  } finally {
    try { api.CloseModel(model); } catch { /* the model is gone either way */ }
  }
}
