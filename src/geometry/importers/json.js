// Structured-text geometry: generic JSON and XML documents holding node / member / route / offset tables, and
// the application's own JSON exchange format (written by exporters.js).
//
// Generic JSON and XML have no geometry schema of their own. The reader looks for arrays of records (or
// repeated XML elements), flattens each into a table and hands the tables to the same role mapping used for
// CSV files (see tables.js): what was recognised is returned in `meta.schema` and `meta.mapping`, and
// `opts.mapping` / `opts.schema` override it.

import { inputText, baseName, finish, base64Decode, cleanName } from './common.js';
import { makeTable, modelsFromTables } from './tables.js';
import { parseXml } from '../xml.js';

const isObj = (v) => v !== null && typeof v === 'object' && !Array.isArray(v) && !ArrayBuffer.isView(v);
const isNumArr = (v) => Array.isArray(v) && v.length > 0 && v.every((x) => typeof x === 'number');
const POS_KEYS = /^(position|pos|coords?|coordinates?|xyz|location|loc|point|p|centre|center|origin)$/i;
const START_KEYS = /^(start|a|from|end1|p1|starpoint|startpoint|i)$/i;
const END_KEYS = /^(end|b|to|end2|p2|endpoint|j)$/i;

/** Flatten one record: nested objects are merged, [x, y, z] arrays under position-like keys become columns. */
function flatten(rec, out = {}, prefix = '') {
  for (const k in rec) {
    const v = rec[k];
    if (isNumArr(v) && v.length >= 2 && v.length <= 3) {
      const ax = ['x', 'y', 'z'];
      if (POS_KEYS.test(k)) v.forEach((c, i) => { out[ax[i]] = c; });
      else if (START_KEYS.test(k)) v.forEach((c, i) => { out[`${ax[i]}1`] = c; });
      else if (END_KEYS.test(k)) v.forEach((c, i) => { out[`${ax[i]}2`] = c; });
      else v.forEach((c, i) => { out[`${prefix}${k}_${ax[i]}`] = c; });
    } else if (isObj(v)) {
      if (POS_KEYS.test(k)) flatten(v, out, '');
      else if (START_KEYS.test(k) && ('x' in v || 'X' in v)) { for (const a in v) out[`${a.toLowerCase()}1`] = v[a]; } else if (END_KEYS.test(k) && ('x' in v || 'X' in v)) { for (const a in v) out[`${a.toLowerCase()}2`] = v[a]; } else flatten(v, out, `${prefix}${k}_`);
    } else if (v === null || typeof v !== 'object') out[prefix + k] = v;
  }
  return out;
}

function tableFromRecords(name, recs) {
  const flat = recs.map((r) => flatten(r));
  const headers = [];
  const seen = new Set();
  for (const f of flat) for (const k in f) if (!seen.has(k)) { seen.add(k); headers.push(k); }
  if (!headers.length) return null;
  return makeTable(name, headers, headers.map((h) => flat.map((f) => f[h])));
}

function tableFromRows(name, rows) {
  const w = Math.max(...rows.map((r) => r.length));
  const t = makeTable(name, Array.from({ length: w }, (_, i) => `c${i + 1}`), Array.from({ length: w }, (_, i) => rows.map((r) => r[i])));
  t.headerless = true;
  return t;
}

/** Collect every array of records in a JSON value as a table. */
export function tablesFromValue(root, name) {
  const tables = [];
  const visit = (v, nm, depth) => {
    if (depth > 6 || v === null || v === undefined) return;
    if (Array.isArray(v)) {
      if (!v.length) return;
      if (v.every(isObj)) {
        // hull stations: [{ x, offsets: [[y, z], …] }]
        const offKey = ['offsets', 'points', 'halfBreadths', 'section'].find((k) => v.every((r) => Array.isArray(r[k])));
        const xKey = ['x', 'station', 'position', 'X'].find((k) => v.every((r) => typeof r[k] === 'number'));
        if (offKey && xKey && v.every((r) => r[offKey].every((p) => isNumArr(p) && p.length === 2 || (isObj(p) && 'y' in p && 'z' in p)))) {
          const st = [], hb = [], wl = [];
          for (const r of v) for (const p of r[offKey]) { st.push(r[xKey]); hb.push(Array.isArray(p) ? p[0] : p.y); wl.push(Array.isArray(p) ? p[1] : p.z); }
          tables.push(makeTable(nm, ['station', 'half_breadth', 'waterline'], [st, hb, wl]));
          return;
        }
        const t = tableFromRecords(nm, v);
        if (t) tables.push(t);
        // records may themselves hold arrays of records (a route with its points, a line with its nodes)
        if (v.length <= 200) v.forEach((r, i) => { for (const k in r) if (Array.isArray(r[k]) && r[k].length && (r[k].every(isObj) || r[k].every(isNumArr))) visit(r[k], cleanName(r.name ?? r.id ?? `${nm} ${i + 1}`), depth + 1); });
        return;
      }
      if (v.every(isNumArr)) { tables.push(tableFromRows(nm, v)); return; }
      return;
    }
    if (!isObj(v)) return;
    const keys = Object.keys(v);
    // columnar: { x: […], y: […], z: […] }
    const arrs = keys.filter((k) => Array.isArray(v[k]) && v[k].length && v[k].every((x) => x === null || typeof x !== 'object'));
    if (arrs.length >= 2 && arrs.every((k) => v[k].length === v[arrs[0]].length) && arrs.some((k) => /^(x|lon|easting|kp|r|station|longitude)$/i.test(k))) {
      tables.push(makeTable(nm, arrs, arrs.map((k) => v[k])));
      return;
    }
    // id → coordinates map: { "N1": [0, 0, 0], … }
    if (keys.length >= 2 && keys.every((k) => isNumArr(v[k]) && v[k].length >= 2 && v[k].length <= 3)) {
      tables.push(makeTable(nm, ['id', 'x', 'y', 'z'], [keys, keys.map((k) => v[k][0]), keys.map((k) => v[k][1]), keys.map((k) => v[k][2] ?? 0)]));
      return;
    }
    const arrayKeys = keys.filter((k) => Array.isArray(v[k]) && v[k].length && typeof v[k][0] === 'object');
    for (const k of keys) {
      const child = v[k];
      if (Array.isArray(child)) visit(child, arrayKeys.length === 1 && depth > 0 && !/s$/i.test(k) ? nm : k, depth + 1);
      else if (isObj(child)) visit(child, k, depth + 1);
    }
  };
  visit(root, name, 0);
  return tables;
}

function docHints(root) {
  const o = isObj(root) ? root : {};
  const find = (...names) => { for (const n of names) for (const k in o) if (k.toLowerCase() === n && (typeof o[k] === 'string' || typeof o[k] === 'number')) return o[k]; return null; };
  let unit = find('units', 'unit', 'lengthunit', 'lengthunits', 'length_unit');
  if (!unit && isObj(o.units)) unit = o.units.length || o.units.Length || null;
  let crs = find('crs', 'epsg', 'srs', 'srid');
  if (typeof crs === 'number') crs = `EPSG:${crs}`;
  if (crs && /^\d{4,5}$/.test(String(crs))) crs = `EPSG:${crs}`;
  return { unitHint: unit ? String(unit) : null, crsHint: crs ? String(crs) : null };
}

/** Generic JSON holding tables of geometry records. */
export function parseTableJson(input, opts = {}) {
  let root;
  try { root = JSON.parse(inputText(input)); } catch (e) { throw new Error(`${input.name || 'JSON'}: not valid JSON (${e.message}).`); }
  if (isObj(root) && root.format === OA_FORMAT) return parseOaJson(input, opts, root);
  const name = baseName(input.name);
  const tables = tablesFromValue(root, name);
  if (!tables.length) throw new Error(`${input.name || 'JSON'}: no arrays of records were found. Expected, for example, { "nodes": [{ "id", "x", "y", "z" }], "members": [{ "a", "b", "D", "t" }] }, { "points": [[x, y, z], …] } or an array of { "kp", "easting", "northing", "depth" } objects.`);
  return modelsFromTables(tables, opts, { name, format: 'table-json', ...docHints(root) });
}

function xmlToValue(el) {
  const keys = Object.keys(el.attrs);
  if (!el.children.length && !keys.length) {
    const t = el.text.trim();
    if (t === '') return null;
    if (/^[-+.\d][-+.\deE\s,;]*$/.test(t)) {
      const p = t.split(/[\s,;]+/).filter(Boolean).map(Number);
      if (p.every((v) => v === v)) return p.length === 1 ? p[0] : p;
    }
    return t;
  }
  const o = {};
  for (const k of keys) { if (/^xmlns|:schemaLocation$/.test(k)) continue; const v = el.attrs[k]; const n = Number(v); o[k.includes(':') ? k.slice(k.indexOf(':') + 1) : k] = v !== '' && n === n ? n : v; }
  const groups = new Map();
  for (const c of el.children) { let g = groups.get(c.name); if (!g) groups.set(c.name, (g = [])); g.push(c); }
  for (const [nm, list] of groups) o[nm] = list.length > 1 ? list.map(xmlToValue) : xmlToValue(list[0]);
  const t = el.text.trim();
  if (t && !el.children.length) {
    const p = t.split(/[\s,;]+/).map(Number);
    o[p.length >= 2 && p.length <= 3 && p.every((v) => v === v) ? 'position' : 'text'] = p.length >= 2 && p.length <= 3 && p.every((v) => v === v) ? p : t;
  }
  return o;
}

/** Generic XML holding repeated record elements (nodes, members, route points, stations …). */
export function parseXmlModel(input, opts = {}) {
  const root = parseXml(inputText(input));
  const name = baseName(input.name);
  const value = xmlToValue(root);
  const tables = tablesFromValue(isObj(value) ? value : { [root.name]: value }, name);
  if (!tables.length) throw new Error(`${input.name || 'XML'}: root element <${root.tag}> is not a format this application knows, and no repeated record elements (for example <Node id x y z/> … <Member start end section/>) were found in it.`);
  return modelsFromTables(tables, opts, { name, format: 'xml-model', ...docHints(value), extra: { rootElement: root.tag } });
}

// ---------------------------------------------------------------------------------------------------------
// Application JSON
// ---------------------------------------------------------------------------------------------------------

export const OA_FORMAT = 'offshore-analyses-geometry';
export const OA_VERSION = 1;
const TYPED = { Float32Array, Float64Array, Uint32Array, Int32Array, Uint16Array, Int16Array, Uint8Array, Int8Array };

function revive(v) {
  if (Array.isArray(v)) return v.map(revive);
  if (!isObj(v)) return v;
  if (typeof v.$typed === 'string') {
    const T = TYPED[v.$typed];
    if (!T) throw new Error(`application JSON: unknown typed-array type "${v.$typed}".`);
    if (typeof v.base64 === 'string') {
      const b = base64Decode(v.base64);
      if (b.length % T.BYTES_PER_ELEMENT) throw new Error(`application JSON: ${v.$typed} payload length ${b.length} is not a multiple of ${T.BYTES_PER_ELEMENT} bytes.`);
      // copy: the decoded bytes may sit at an unaligned offset
      return new T(b.slice().buffer);
    }
    // JSON has no NaN: null marks no-data in value lists
    return T.from(v.values || [], (x) => (x === null ? NaN : x));
  }
  const o = {};
  for (const k in v) o[k] = revive(v[k]);
  return o;
}

/** The application's own exchange format: `{ format, version, models: [...] }` with tagged typed arrays. */
export function parseOaJson(input, opts = {}, parsed) {
  let root = parsed;
  if (!root) { try { root = JSON.parse(inputText(input)); } catch (e) { throw new Error(`${input.name || 'JSON'}: not valid JSON (${e.message}).`); } }
  if (!isObj(root) || root.format !== OA_FORMAT) throw new Error(`${input.name || 'JSON'}: not an application geometry file (expected "format": "${OA_FORMAT}").`);
  if (!(root.version >= 1) || root.version > OA_VERSION) throw new Error(`${input.name || 'JSON'}: application geometry version ${root.version} is newer than this build reads (${OA_VERSION}).`);
  if (!Array.isArray(root.models) || !root.models.length) throw new Error(`${input.name || 'JSON'}: the file holds no models (expected a non-empty "models" array).`);
  const kinds = ['surface', 'volume', 'frame', 'polyline', 'grid', 'points', 'section'];
  const models = root.models.map((m, i) => {
    const r = revive(m);
    if (!kinds.includes(r.kind)) throw new Error(`${input.name || 'JSON'}: model ${i + 1} has unknown kind "${r.kind}".`);
    return r;
  });
  const hist = {};
  for (const m of models) hist[m.kind] = (hist[m.kind] || 0) + 1;
  return finish(models, [], { format: 'oa-json', units: { unit: 'm', scale: 1, basis: 'file' }, crs: root.crs ?? null, elementTypes: hist, skipped: {}, version: root.version, source: root.source ?? null });
}
