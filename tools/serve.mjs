// Minimal static file server for local use: `npm run serve` then open http://localhost:8080
import { createServer } from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { extname, join, normalize, resolve, sep } from 'node:path';

const root = resolve(process.argv[2] || '.');
const port = Number(process.env.PORT || 8080);
const types = { '.html': 'text/html; charset=utf-8', '.js': 'text/javascript; charset=utf-8', '.mjs': 'text/javascript; charset=utf-8', '.css': 'text/css; charset=utf-8', '.json': 'application/json', '.webmanifest': 'application/manifest+json', '.svg': 'image/svg+xml', '.png': 'image/png', '.wasm': 'application/wasm', '.txt': 'text/plain; charset=utf-8' };

createServer(async (req, res) => {
  try {
    let p = normalize(join(root, decodeURIComponent(new URL(req.url, 'http://x').pathname)));
    if (p !== root && !p.startsWith(root + sep)) { res.writeHead(403).end(); return; }
    if ((await stat(p)).isDirectory()) p = join(p, 'index.html');
    const body = await readFile(p);
    res.writeHead(200, { 'content-type': types[extname(p)] || 'application/octet-stream', 'cache-control': 'no-cache' }).end(body);
  } catch { res.writeHead(404, { 'content-type': 'text/plain' }).end('Not found'); }
}).listen(port, () => console.log(`Offshore Analyses on http://localhost:${port}`));
