// Suite and module registry. Suite metadata comes from the specification data; module code is loaded on
// demand so the first paint never waits for the engines.

const pad = (n) => String(n).padStart(2, '0');
let indexP;
const specs = new Map();
const suites = new Map();

export const SHORT = ['', 'Fixed platforms & jackets', 'Topsides & decks', 'Wind turbine foundations', 'Subsea structures & mudmats', 'Wharves, jetties & quays',
  'Dolphins', 'Other marine structures', 'Pipelines & flowlines', 'Risers', 'Moorings & anchors', 'Floating structures', 'Umbilicals',
  'Power & telecom cables', 'Flexible pipes', 'Installation & marine ops', 'Renewable-energy systems', 'Economics'];

const fetchJson = (rel) => fetch(new URL(rel, import.meta.url)).then((r) => { if (!r.ok) throw new Error(`Could not load ${rel}`); return r.json(); });

/** [{ n, title, subtitle, scope, analyses[] }] for all 17 suites. */
export function suiteIndex() { return (indexP ||= fetchJson('../spec/index.json')); }
export function suiteSpec(n) { if (!specs.has(n)) specs.set(n, fetchJson(`../spec/suite-${pad(n)}.json`)); return specs.get(n); }
export const commonSpec = (() => { let p; return () => (p ||= fetchJson('../spec/common.json')); })();

let manifestP;
const manifest = () => (manifestP ||= fetchJson('../suites/manifest.json').catch(() => ({ suites: [] })));

/** Catalogue entries of one suite's modules (titles, coverage, wiring) — no calculation code is loaded. */
export async function suiteModules(n) {
  if (!suites.has(n)) suites.set(n, manifest().then((m) => m.suites.find((s) => s.n === n)?.modules || []));
  return suites.get(n);
}

const described = new Map();
/**
 * Full declaration of one module (inputs, mesh ladders, verification cases) as plain data. The module's code
 * stays in the analysis worker, so the page never loads the engines itself.
 */
export async function findModule(id) {
  const n = Number(id.slice(1, 3));
  if (!(await suiteModules(n)).some((m) => m.id === id)) return null;
  if (!described.has(id)) {
    const { task } = await import('./runner.js');
    described.set(id, task('describe', id, {}).catch((err) => { described.delete(id); throw err; }));
  }
  return described.get(id);
}

/** For each analysis in the specification, the modules that cover it. */
export async function coverage(n) {
  const [idx, mods] = await Promise.all([suiteIndex(), suiteModules(n)]);
  const suite = idx.find((s) => s.n === n);
  return suite.analyses.map((name, i) => ({ index: i + 1, name, modules: mods.filter((m) => (m.covers || []).includes(i + 1)) }));
}

/** Cross-suite wiring: who publishes and who consumes each bus key, over every loaded suite. */
export async function wiring() {
  const out = new Map();
  const idx = await suiteIndex();
  await Promise.all(idx.map(async (s) => {
    for (const m of await suiteModules(s.n)) {
      for (const inp of m.inputs || []) if (inp.bind) { const e = out.get(inp.bind) || { key: inp.bind, publishers: [], consumers: [] }; e.consumers.push({ suite: s.n, id: m.id, title: m.title }); out.set(inp.bind, e); }
      for (const key of m.publishes || []) { const e = out.get(key) || { key, publishers: [], consumers: [] }; e.publishers.push({ suite: s.n, id: m.id, title: m.title }); out.set(key, e); }
    }
  }));
  return [...out.values()].sort((a, b) => a.key.localeCompare(b.key));
}
