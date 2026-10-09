// Main-thread side of the analysis worker, with an in-page fallback where module workers are unavailable.
import { getProject, busValues } from './store.js';

let worker, seq = 0, broken = false;
const pending = new Map();

function ensureWorker() {
  if (worker || broken) return worker;
  try {
    worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = ({ data }) => {
      const p = pending.get(data.id);
      if (!p) return;
      if (data.progress != null) { p.onProgress?.(data.progress, data.message); return; }
      pending.delete(data.id);
      if (data.error) p.reject(new Error(data.error)); else p.resolve(data.result);
    };
    worker.onerror = (e) => {
      // a worker that fails to load (old browser, blocked module workers) must not strand the queued jobs
      broken = true; worker = null;
      for (const [, p] of pending) p.retry();
      pending.clear();
      console.warn('Analysis worker unavailable, computing in the page instead.', e.message);
    };
  } catch { broken = true; worker = null; }
  return worker;
}

async function inPage(type, moduleId, values, options, ctxData, onProgress) {
  const [studies, lib] = await Promise.all([import('./studies.js'), import('../library/materials.js')]);
  const n = moduleId.slice(1, 3);
  const suite = (await import(`../suites/s${n}/index.js`)).default;
  const mod = suite.modules.find((m) => m.id === moduleId);
  if (!mod) throw new Error(`Module ${moduleId} is not available.`);
  lib.setCustomMaterials(ctxData.materials);
  const ctx = { site: ctxData.site, bus: ctxData.bus, material: lib.material, progress: onProgress || (() => {}), signal: new AbortController().signal };
  if (type === 'run') return mod.compute(values, ctx);
  if (type === 'verify') return studies.verify(mod, ctx);
  return studies[type](mod, values, ctx, { ...options, onProgress });
}

/** task: 'run' | 'verify' | 'sensitivity' | 'sweep' | 'monteCarlo' | 'convergence' | 'calibrate' | 'validation' | 'advise' */
export function task(type, moduleId, values, options = {}, onProgress) {
  const p = getProject();
  const ctxData = { site: { ...p.site.data, lat: p.site.lat, lon: p.site.lon, name: p.site.name }, bus: busValues(), materials: p.materials };
  const w = ensureWorker();
  if (!w) return inPage(type, moduleId, values, options, ctxData, onProgress);
  const id = ++seq;
  let cancel;
  const promise = new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject, onProgress, retry: () => inPage(type, moduleId, values, options, ctxData, onProgress).then(resolve, reject) });
    cancel = () => { w.postMessage({ type: 'abort', target: id }); pending.delete(id); reject(new Error('Cancelled')); };
    w.postMessage({ id, type, moduleId, values, options, ctxData });
  });
  promise.cancel = cancel;
  return promise;
}
