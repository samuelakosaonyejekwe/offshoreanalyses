// Analysis worker: runs module computations and generic studies off the main thread so the UI stays responsive.
import * as studies from './studies.js';
import { material, setCustomMaterials } from '../library/materials.js';

const suites = new Map();
async function findModule(id) {
  const n = Number(id.slice(1, 3));
  if (!suites.has(n)) suites.set(n, import(`../suites/s${String(n).padStart(2, '0')}/index.js`).then((m) => m.default));
  const suite = await suites.get(n);
  const mod = suite.modules.find((m) => m.id === id);
  if (!mod) throw new Error(`Module ${id} is not available.`);
  return mod;
}

// the module declaration without its functions: what the page needs to build the form and the studies
const describe = (mod) => JSON.parse(JSON.stringify(mod, (k, v) => (typeof v === 'function' ? undefined : v === Infinity ? 1e308 : v === -Infinity ? -1e308 : v)));

const aborts = new Map();
self.onmessage = async ({ data }) => {
  const { id, type, moduleId, values, options, ctxData } = data;
  if (type === 'abort') { aborts.get(data.target)?.abort(); return; }
  const ac = new AbortController();
  aborts.set(id, ac);
  const progress = (fraction, message) => self.postMessage({ id, progress: fraction, message });
  try {
    setCustomMaterials(ctxData?.materials);
    const ctx = { site: ctxData?.site || {}, bus: ctxData?.bus || {}, material, progress, signal: ac.signal };
    const mod = await findModule(moduleId);
    let result;
    if (type === 'describe') result = describe(mod);
    else if (type === 'run') result = await mod.compute(values, ctx);
    else if (type === 'verify') result = await studies.verify(mod, ctx);
    else if (typeof studies[type] === 'function') result = await studies[type](mod, values, ctx, { ...options, onProgress: progress });
    else throw new Error(`Unknown task ${type}`);
    self.postMessage({ id, result });
  } catch (err) {
    self.postMessage({ id, error: err?.message || String(err) });
  } finally {
    aborts.delete(id);
  }
};
