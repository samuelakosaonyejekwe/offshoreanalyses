// Publishes a second copy of the app on a host that does not depend on GitHub, and records it in mirrors.json so
// installed copies fall back to it when the primary host is unreachable.
//
//   node tools/mirror.mjs cloudflare [project-name]   Cloudflare Pages   (opens a browser sign-in the first time)
//   node tools/mirror.mjs netlify                     Netlify            (opens a browser sign-in the first time)
//   node tools/mirror.mjs gitlab <git-remote-url>     GitLab Pages       (pushes the repository; .gitlab-ci.yml publishes it)
//   node tools/mirror.mjs stage [folder]              only write the publishable folder, for any other static host
//   node tools/mirror.mjs add <https://address/>      record a mirror that is already online
//
// The app is a folder of static files, so any static host works: upload the staged folder and run `add`.
import { cp, mkdir, rm, readFile, writeFile, readdir, stat } from 'node:fs/promises';
import { spawnSync } from 'node:child_process';
import { join, resolve } from 'node:path';
import { tmpdir } from 'node:os';

const root = process.cwd();
// what a browser needs; tests, tools and Node-only vendor builds stay behind
const PUBLISH = ['index.html', 'sw.js', 'manifest.webmanifest', 'mirrors.json', 'precache.json', 'css', 'assets', 'src', 'dist', 'vendor', 'profiles', 'docs', 'README.md'];
const SKIP = /(^|[\\/])(node|node_modules)([\\/]|$)|\.map$|\.d\.ts$/;

export async function stage(dir) {
  const out = resolve(dir || join(tmpdir(), 'offshoreanalyses-site'));
  if (out === root || root.startsWith(out)) throw new Error('refusing to stage over the repository itself');
  await rm(out, { recursive: true, force: true });
  await mkdir(out, { recursive: true });
  for (const item of PUBLISH) await cp(join(root, item), join(out, item), { recursive: true, filter: (src) => !SKIP.test(src.slice(root.length)) }).catch((e) => { if (e.code !== 'ENOENT') throw e; });
  // hosts that honour a _headers file (Cloudflare Pages, Netlify): let the service worker update promptly and
  // let other copies of the app use this one as a fallback source
  await writeFile(join(out, '_headers'), '/*\n  Access-Control-Allow-Origin: *\n/sw.js\n  Cache-Control: no-cache\n/precache.json\n  Cache-Control: no-cache\n/mirrors.json\n  Cache-Control: no-cache\n');
  let n = 0, bytes = 0, largest = ['', 0];
  const walk = async (p) => { for (const e of await readdir(p, { withFileTypes: true })) { const f = join(p, e.name); if (e.isDirectory()) await walk(f); else { const s = (await stat(f)).size; n++; bytes += s; if (s > largest[1]) largest = [f.slice(out.length + 1), s]; } } };
  await walk(out);
  console.log(`Staged ${n} files, ${(bytes / 1e6).toFixed(1)} MB in ${out} (largest: ${largest[0]}, ${(largest[1] / 1e6).toFixed(1)} MB)`);
  return out;
}

export async function add(url, host) {
  let u;
  try { u = new URL(url); } catch { throw new Error(`not an address: ${url}`); }
  if (u.protocol !== 'https:') throw new Error('a mirror must be served over https (service workers and installation require it)');
  const base = u.href.replace(/index\.html$/, '').replace(/\/?$/, '/');
  const file = join(root, 'mirrors.json');
  const m = JSON.parse(await readFile(file, 'utf8'));
  if (!m.mirrors.some((x) => x.url.replace(/index\.html$/, '') === base)) m.mirrors.push({ url: base, host: host || u.hostname, note: 'independent of GitHub' });
  if (!m.assetBases.includes(base)) m.assetBases.push(base);
  await writeFile(file, `${JSON.stringify(m, null, 2)}\n`);
  console.log(`Recorded ${base} in mirrors.json. Run "npm run precache", commit and push so installed copies learn about it.`);
}

const run = (cmd, args, opts = {}) => { const r = spawnSync(cmd, args, { stdio: 'inherit', shell: process.platform === 'win32', ...opts }); if (r.status !== 0) throw new Error(`${cmd} ${args.join(' ')} did not finish (exit ${r.status ?? r.error?.message})`); };

if (import.meta.url === `file://${process.argv[1]}` || process.argv[1]?.endsWith('mirror.mjs')) {
  const [what, arg] = process.argv.slice(2);
  try {
    if (what === 'stage') await stage(arg);
    else if (what === 'add') { if (!arg) throw new Error('give the address of the mirror'); await add(arg); }
    else if (what === 'cloudflare') {
      const name = arg || 'offshoreanalyses';
      const dir = await stage();
      run('npx', ['--yes', 'wrangler@latest', 'pages', 'deploy', dir, '--project-name', name, '--branch', 'main', '--commit-dirty=true']);
      console.log(`\nIf the address printed above is https://${name}.pages.dev, record it with:\n  node tools/mirror.mjs add https://${name}.pages.dev/`);
    } else if (what === 'netlify') {
      const dir = await stage();
      run('npx', ['--yes', 'netlify-cli@latest', 'deploy', '--prod', '--dir', dir]);
      console.log('\nRecord the "Website URL" printed above with:\n  node tools/mirror.mjs add <that address>');
    } else if (what === 'gitlab') {
      if (!arg) throw new Error('give the GitLab repository address, e.g. git@gitlab.com:<user>/offshoreanalyses.git');
      const has = spawnSync('git', ['remote', 'get-url', 'gitlab'], { encoding: 'utf8' });
      if (has.status !== 0) run('git', ['remote', 'add', 'gitlab', arg]); else if (has.stdout.trim() !== arg) run('git', ['remote', 'set-url', 'gitlab', arg]);
      run('git', ['push', 'gitlab', 'HEAD:main']);
      const mt = arg.match(/gitlab\.com[:/]([^/]+)\/(.+?)(\.git)?$/);
      console.log(`\nGitLab builds the site from .gitlab-ci.yml in a minute or two.${mt ? ` Then record it with:\n  node tools/mirror.mjs add https://${mt[1]}.gitlab.io/${mt[2]}/` : ''}`);
    } else {
      console.log((await readFile(new URL(import.meta.url), 'utf8')).split('\n').slice(0, 10).map((l) => l.replace(/^\/\/ ?/, '')).join('\n'));
    }
  } catch (e) { console.error(`mirror: ${e.message}`); process.exit(1); }
}
