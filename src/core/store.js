// Project state and persistence. Everything lives in the browser (IndexedDB) so the app works offline
// and nothing depends on a server: projects, imported geometry, cached live data.

const DB_NAME = 'offshore-analyses';
const STORES = ['kv', 'cache', 'geometry'];
let dbp;

function db() {
  if (dbp) return dbp;
  dbp = new Promise((resolve, reject) => {
    const req = indexedDB.open(DB_NAME, 1);
    req.onupgradeneeded = () => { for (const s of STORES) if (!req.result.objectStoreNames.contains(s)) req.result.createObjectStore(s); };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
  return dbp;
}
function tx(store, mode, fn) {
  return db().then((d) => new Promise((resolve, reject) => {
    const t = d.transaction(store, mode);
    const r = fn(t.objectStore(store));
    t.oncomplete = () => resolve(r?.result);
    t.onerror = () => reject(t.error);
    t.onabort = () => reject(t.error);
  }));
}
// Private-mode browsers may refuse IndexedDB; fall back to memory so the app still runs for the session.
const mem = new Map();
const guard = (p, fallback) => p.catch(() => fallback());
export const idb = {
  get: (s, k) => guard(tx(s, 'readonly', (o) => o.get(k)), () => mem.get(`${s}:${k}`)),
  set: (s, k, v) => guard(tx(s, 'readwrite', (o) => o.put(v, k)), () => { mem.set(`${s}:${k}`, v); }),
  del: (s, k) => guard(tx(s, 'readwrite', (o) => o.delete(k)), () => { mem.delete(`${s}:${k}`); }),
  keys: (s) => guard(tx(s, 'readonly', (o) => o.getAllKeys()), () => [...mem.keys()].filter((k) => k.startsWith(`${s}:`)).map((k) => k.slice(s.length + 1))),
};

const listeners = new Map();
export function on(evt, fn) { if (!listeners.has(evt)) listeners.set(evt, new Set()); listeners.get(evt).add(fn); return () => listeners.get(evt).delete(fn); }
export function emit(evt, data) { for (const fn of listeners.get(evt) || []) { try { fn(data); } catch (e) { console.error(e); } } }

const blank = (name = 'Untitled project') => ({
  id: `p${Date.now().toString(36)}`, name, created: new Date().toISOString(), currency: 'USD',
  site: { name: '', lat: null, lon: null, pulled: null, data: {} },
  bus: {},        // key -> { value, from, at }
  inputs: {},     // moduleId -> { key: value }
  runs: {},       // moduleId -> { at, status, kpis, quantities, maxUC, icbc }
  geometry: [],   // [{ id, name, kind, format, stats }]
  materials: [],  // user-defined material records
  factors: {},    // 'moduleId.key' -> { value, edition, by, at }: code factors confirmed against the governing standard
  evidence: {},   // moduleId -> { calibration?, validation? }: what measurements have demonstrated, and over what range
  notes: '',
});

let project = blank();
let saveTimer;

export const getProject = () => project;

export async function boot() {
  const current = await idb.get('kv', 'currentProject');
  const saved = current && await idb.get('kv', `project:${current}`);
  if (saved) project = { ...blank(), ...saved };
  else await persist();
  emit('project', project);
  return project;
}
async function persist() {
  await idb.set('kv', `project:${project.id}`, project);
  await idb.set('kv', 'currentProject', project.id);
}
export function save() { clearTimeout(saveTimer); saveTimer = setTimeout(() => persist().catch(console.error), 300); }

export async function listProjects() {
  const keys = (await idb.keys('kv')).filter((k) => String(k).startsWith('project:'));
  const out = [];
  for (const k of keys) { const p = await idb.get('kv', k); if (p) out.push({ id: p.id, name: p.name, created: p.created, runs: Object.keys(p.runs || {}).length }); }
  return out.sort((a, b) => (a.created < b.created ? 1 : -1));
}
export async function newProject(name) { await persist(); project = blank(name); await persist(); emit('project', project); return project; }
export async function openProject(id) { const p = await idb.get('kv', `project:${id}`); if (!p) throw new Error('Project not found'); project = { ...blank(), ...p }; await idb.set('kv', 'currentProject', id); emit('project', project); return project; }
export async function deleteProject(id) { await idb.del('kv', `project:${id}`); if (project.id === id) { project = blank(); await persist(); emit('project', project); } }
export function importProject(obj) { project = { ...blank(), ...obj, id: `p${Date.now().toString(36)}` }; save(); emit('project', project); return project; }

/** Publish values on the project bus. Values are SI (docs/CONTRACT.md §3). */
export function publish(values, from) {
  const at = new Date().toISOString();
  let changed = false;
  for (const [k, v] of Object.entries(values || {})) {
    if (typeof v !== 'number' || !Number.isFinite(v)) continue;
    project.bus[k] = { value: v, from, at }; changed = true;
  }
  if (changed) { save(); emit('bus', project.bus); }
}
export const busValues = () => Object.fromEntries(Object.entries(project.bus).map(([k, v]) => [k, v.value]));

/** Confirmed code factors of one module as { key: value }. */
export const confirmedFactors = (moduleId) => Object.fromEntries(Object.entries(project.factors || {}).filter(([k]) => k.startsWith(`${moduleId}.`)).map(([k, v]) => [k.slice(moduleId.length + 1), v.value]));
export function confirmFactor(moduleId, key, rec) { project.factors ||= {}; if (rec) project.factors[`${moduleId}.${key}`] = { ...rec, at: new Date().toISOString() }; else delete project.factors[`${moduleId}.${key}`]; save(); emit('factors', project.factors); }
export function recordEvidence(moduleId, kind, rec) { project.evidence ||= {}; (project.evidence[moduleId] ||= {})[kind] = { ...rec, at: new Date().toISOString() }; save(); emit('evidence', project.evidence); }
export const evidenceOf = (moduleId) => (project.evidence || {})[moduleId] || {};

export function setInputs(moduleId, values) { project.inputs[moduleId] = values; save(); }
export function recordRun(moduleId, summary) { project.runs[moduleId] = { at: new Date().toISOString(), ...summary }; save(); emit('run', { moduleId, summary }); }

export async function addGeometry(entry, models) {
  const id = `g${Date.now().toString(36)}${Math.floor(performance.now() % 1000)}`;
  await idb.set('geometry', id, models);
  project.geometry.push({ id, ...entry });
  save(); emit('geometry', project.geometry);
  return id;
}
export const loadGeometry = (id) => idb.get('geometry', id);
export async function removeGeometry(id) { await idb.del('geometry', id); project.geometry = project.geometry.filter((g) => g.id !== id); save(); emit('geometry', project.geometry); }

export const prefs = {
  get: (k, d) => { try { const v = localStorage.getItem(`oa:${k}`); return v == null ? d : JSON.parse(v); } catch { return d; } },
  set: (k, v) => { try { localStorage.setItem(`oa:${k}`, JSON.stringify(v)); } catch { /* storage unavailable */ } },
};
