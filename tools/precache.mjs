// Writes precache.json: every file the installed app needs offline, with content hashes so an update
// re-downloads only what changed. Run before publishing: `npm run precache`.
import { readdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { join, relative, sep } from 'node:path';

const root = process.cwd();
const include = ['index.html', 'manifest.webmanifest', 'mirrors.json', 'css', 'assets', 'src', 'vendor'];
// Node-only vendor builds serve the test suite, not the browser
const skip = /(^|\/)(\.|LICENSE|README)|\.map$|\.d\.ts$|^vendor\/[^/]+\/node\/|package\.json$/;

async function walk(p, out) {
  for (const e of await readdir(p, { withFileTypes: true }).catch(() => [])) {
    const full = join(p, e.name);
    if (e.isDirectory()) await walk(full, out); else out.push(full);
  }
  return out;
}
const files = [];
for (const item of include) {
  const full = join(root, item);
  const list = item.includes('.') && !item.endsWith('/') && !['css', 'assets', 'src', 'vendor'].includes(item) ? [full] : await walk(full, []);
  for (const f of list) {
    const path = relative(root, f).split(sep).join('/');
    if (skip.test(path)) continue;
    const buf = await readFile(f).catch(() => null);
    if (!buf) continue;
    files.push({ path, hash: createHash('sha256').update(buf).digest('hex').slice(0, 16), size: buf.length });
  }
}
files.sort((a, b) => a.path.localeCompare(b.path));
const version = createHash('sha256').update(files.map((f) => f.path + f.hash).join('|')).digest('hex').slice(0, 12);
await writeFile(join(root, 'precache.json'), JSON.stringify({ version, generated: new Date().toISOString(), files }));
console.log(`precache.json: ${files.length} files, ${(files.reduce((s, f) => s + f.size, 0) / 1e6).toFixed(1)} MB, version ${version}`);
