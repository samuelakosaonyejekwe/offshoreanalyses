// The packed suite files in dist/ are what users actually run. They must be built from the current sources
// and behave identically to them.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { sourceHash } from '../tools/bundle.mjs';
import { defaults } from '../src/core/studies.js';
import { material } from '../src/library/materials.js';

const ctx = { site: {}, bus: {}, material, progress() {}, signal: new AbortController().signal };
const have = existsSync(new URL('../dist/bundles.json', import.meta.url));

test('dist/ is built from the current sources (run `npm run build` if this fails)', { skip: !have && 'dist/ has not been built' }, async () => {
  const meta = JSON.parse(await readFile(new URL('../dist/bundles.json', import.meta.url), 'utf8'));
  assert.equal(meta.source, await sourceHash(), 'sources changed since dist/ was built');
});

for (let n = 1; n <= 17; n++) {
  const nn = String(n).padStart(2, '0');
  test(`suite ${n}: packed file matches its sources`, { skip: !have && 'dist/ has not been built' }, async () => {
    const [{ default: src }, { default: packed }] = await Promise.all([import(`../src/suites/s${nn}/index.js`), import(`../dist/s${nn}.js`)]);
    assert.deepEqual(packed.modules.map((m) => m.id), src.modules.map((m) => m.id));
    // spot-check: the first and last module give identical headline numbers from both builds
    for (const i of [0, src.modules.length - 1]) {
      const values = defaults(src.modules[i]).values;
      const [a, b] = await Promise.all([src.modules[i].compute(structuredClone(values), ctx), packed.modules[i].compute(structuredClone(values), ctx)]);
      assert.deepEqual(b.kpis.map((k) => [k.key, k.value]), a.kpis.map((k) => [k.key, k.value]), src.modules[i].id);
    }
  });
}
