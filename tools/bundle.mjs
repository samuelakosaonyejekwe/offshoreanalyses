// Packs each suite (its modules and every engine they use) into one file under dist/, plus the analysis
// worker, so opening a suite costs one request instead of a few hundred. Run before publishing: `npm run build`.
// The bundler is fetched on demand with npx; nothing is added to the repository's dependencies.
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { readdir, readFile, writeFile, mkdir, rm, stat } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { join } from 'node:path';

const ESBUILD = process.env.ESBUILD || null;
const run = (args) => (ESBUILD ? execFileSync(ESBUILD, args, { stdio: ['ignore', 'pipe', 'pipe'] }) : execFileSync('npx', ['--yes', 'esbuild@0.24.0', ...args], { stdio: ['ignore', 'pipe', 'pipe'] }));

async function walk(dir, out = []) {
  for (const e of await readdir(dir, { withFileTypes: true })) {
    const p = join(dir, e.name);
    if (e.isDirectory()) await walk(p, out); else if (e.name.endsWith('.js')) out.push(p);
  }
  return out;
}
/** Hash of every source file a bundle can contain: the stamp that tells a stale bundle from a fresh one. */
export async function sourceHash() {
  const files = (await Promise.all(['src/engine', 'src/suites', 'src/core', 'src/library', 'src/mesh'].map((d) => walk(d)))).flat().sort();
  const h = createHash('sha256');
  for (const f of files) { h.update(f); h.update(await readFile(f)); }
  return h.digest('hex').slice(0, 16);
}

if (process.argv[1] && import.meta.url.endsWith(process.argv[1].split('/').pop())) {
  await rm('dist', { recursive: true, force: true });
  await mkdir('dist', { recursive: true });
  const entries = [];
  for (let n = 1; n <= 17; n++) { const f = `src/suites/s${String(n).padStart(2, '0')}/index.js`; if (existsSync(f)) entries.push(f); }
  const common = ['--bundle', '--format=esm', '--minify', '--target=es2020', '--legal-comments=none', '--log-level=warning'];
  for (const f of entries) run([f, ...common, `--outfile=dist/${f.split('/')[2]}.js`]);
  run(['src/core/worker.js', ...common, '--outfile=dist/worker.js']);
  const sizes = {};
  for (const f of await readdir('dist')) sizes[f] = (await stat(join('dist', f))).size;
  await writeFile('dist/bundles.json', JSON.stringify({ source: await sourceHash(), built: new Date().toISOString(), files: sizes }));
  console.log(`dist: ${Object.keys(sizes).length} bundles, ${(Object.values(sizes).reduce((s, v) => s + v, 0) / 1e6).toFixed(1)} MB`);
}
