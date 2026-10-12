// Writes src/suites/manifest.json: the catalogue of every analysis module (titles, coverage, wiring, data inputs)
// so list pages never have to load calculation code. Run before publishing: `npm run manifest`.
import { writeFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';

const suites = [];
for (let n = 1; n <= 17; n++) {
  const file = new URL(`../src/suites/s${String(n).padStart(2, '0')}/index.js`, import.meta.url);
  if (!existsSync(file)) { suites.push({ n, modules: [] }); continue; }
  const { default: suite } = await import(file);
  suites.push({ n, modules: suite.modules.map((m) => ({
    id: m.id, suite: n, title: m.title, summary: m.summary, covers: m.covers, fidelity: m.fidelity,
    standards: m.standards, models: m.models, publishes: m.publishes || [], searchTerms: m.searchTerms,
    hasMesh: !!(m.mesh || m.meshTime), verifyCases: (m.verify || []).filter((v) => !v.pending).length,
    // only what list pages need: wiring (bind), the portal's table/series/geometry list, and one search string
    inputText: m.inputs.map((i) => i.label).join(' · '),
    inputs: m.inputs.filter((i) => i.bind || ['table', 'series', 'geometry'].includes(i.type)).map((i) => ({ label: i.label, type: i.type || 'number', ...(i.bind ? { bind: i.bind } : {}), ...(i.kinds ? { kinds: i.kinds } : {}),
      ...(i.columns ? { columns: i.columns.map((c) => ({ label: c.label })) } : {}), ...(i.xLabel ? { xLabel: i.xLabel, yLabel: i.yLabel } : {}) })),
  })) });
}
await writeFile(new URL('../src/suites/manifest.json', import.meta.url), JSON.stringify({ generated: new Date().toISOString(), suites }));
console.log(`manifest: ${suites.reduce((s, x) => s + x.modules.length, 0)} modules in ${suites.filter((x) => x.modules.length).length} suites`);
