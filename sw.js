// Service worker: makes the app installable and fully usable offline, keeps it current, and falls back to
// mirror hosts when the primary host cannot be reached.
const CACHE = 'offshore-analyses-v1';
const CORE = ['./', 'index.html', 'css/app.css', 'manifest.webmanifest', 'mirrors.json', 'precache.json', 'assets/icon.svg', 'assets/icon-192.png', 'assets/icon-512.png', 'assets/world.json',
  'src/app.js', 'src/ui/shell.js', 'src/ui/dom.js', 'src/ui/plot.js', 'src/ui/form.js', 'src/ui/results.js', 'src/core/store.js', 'src/core/registry.js', 'src/core/runner.js', 'src/core/studies.js', 'src/core/worker.js',
  'src/library/materials.js', 'src/data/connectors.js', 'src/ui/views/home.js', 'src/ui/views/suite.js', 'src/ui/views/module.js', 'src/ui/views/pages.js', 'src/spec/index.json'];

const scope = self.registration.scope;
const rel = (url) => (url.startsWith(scope) ? url.slice(scope.length) : null);

self.addEventListener('install', (e) => {
  e.waitUntil((async () => {
    const c = await caches.open(CACHE);
    // one missing file must not block installation; whatever arrives is cached
    await Promise.all(CORE.map((u) => c.add(new Request(u, { cache: 'reload' })).catch(() => {})));
    await self.skipWaiting();
  })());
});

self.addEventListener('activate', (e) => {
  e.waitUntil((async () => {
    for (const k of await caches.keys()) if (k !== CACHE && k.startsWith('offshore-analyses-')) await caches.delete(k);
    await self.clients.claim();
    refreshChanged().catch(() => {});
  })());
});

let mirrorBases;
async function mirrors() {
  if (mirrorBases) return mirrorBases;
  try { const r = await caches.match('mirrors.json'); mirrorBases = ((await r.json()).assetBases || []).filter((b) => !scope.startsWith(b)); } catch { mirrorBases = []; }
  return mirrorBases;
}

async function fromNetwork(request, path) {
  try {
    const r = await fetch(request, { cache: 'no-cache' });
    if (r.ok) return r;
    if (r.status < 500) return r;
    throw new Error(String(r.status));
  } catch (err) {
    // primary host unreachable: the same file is published on the mirror hosts
    for (const base of await mirrors()) {
      try {
        const r = await fetch(base + path, { mode: 'cors', cache: 'no-cache' });
        if (r.ok) {
          // some mirrors serve scripts as text/plain; restore the type so modules still load
          const type = path.endsWith('.js') || path.endsWith('.mjs') ? 'text/javascript' : path.endsWith('.json') || path.endsWith('.webmanifest') ? 'application/json' : path.endsWith('.css') ? 'text/css' : path.endsWith('.html') || path === '' ? 'text/html' : path.endsWith('.svg') ? 'image/svg+xml' : path.endsWith('.wasm') ? 'application/wasm' : r.headers.get('content-type') || 'application/octet-stream';
          return new Response(await r.arrayBuffer(), { status: 200, headers: { 'content-type': type } });
        }
      } catch { /* try the next mirror */ }
    }
    throw err;
  }
}

self.addEventListener('fetch', (e) => {
  const { request } = e;
  if (request.method !== 'GET') return;
  const path = rel(request.url.split('#')[0].split('?')[0]);
  if (path == null) return; // third-party data requests are cached by the page itself, with their own lifetimes
  e.respondWith((async () => {
    const c = await caches.open(CACHE);
    const key = request.mode === 'navigate' ? 'index.html' : request;
    const hit = await c.match(key, { ignoreSearch: true });
    const update = fromNetwork(request.mode === 'navigate' ? new Request('index.html') : request, request.mode === 'navigate' ? 'index.html' : path)
      .then((r) => { if (r.ok) c.put(key, r.clone()).catch(() => {}); return r; });
    if (hit) { e.waitUntil(update.catch(() => {})); return hit; }
    try { return await update; } catch {
      if (request.mode === 'navigate') return (await c.match('index.html')) || new Response('Offline — open the app once while connected to store it on this device.', { status: 503, headers: { 'content-type': 'text/plain' } });
      return new Response('', { status: 504, statusText: 'Offline and not stored on this device' });
    }
  })());
});

async function cacheAll(port) {
  const c = await caches.open(CACHE);
  const list = await fetch('precache.json', { cache: 'no-cache' }).then((r) => r.json()).catch(async () => (await c.match('precache.json'))?.json());
  const files = list?.files || CORE;
  let done = 0, failed = 0;
  const queue = files.slice();
  const worker = async () => {
    while (queue.length) {
      const f = queue.shift();
      try { const hit = await c.match(f.path || f); if (!hit) { const r = await fromNetwork(new Request(f.path || f), f.path || f); if (r.ok) await c.put(f.path || f, r); else failed++; } } catch { failed++; }
      done++;
      if (done % 5 === 0 || done === files.length) port?.postMessage({ progress: done / files.length, done, total: files.length, failed });
    }
  };
  await Promise.all(Array.from({ length: 6 }, worker));
  await c.put('precache-state', new Response(JSON.stringify({ version: list?.version, hashes: Object.fromEntries(files.map((f) => [f.path, f.hash])) })));
}

// Compare the published file list with what was stored and re-fetch only what changed.
async function refreshChanged() {
  const c = await caches.open(CACHE);
  const list = await fetch('precache.json', { cache: 'no-cache' }).then((r) => (r.ok ? r.json() : null));
  if (!list) return;
  const state = await c.match('precache-state').then((r) => r?.json()).catch(() => null);
  if (state?.version === list.version) return;
  const old = state?.hashes || {};
  for (const f of list.files) {
    if (old[f.path] === f.hash) continue;
    if (!(await c.match(f.path))) continue; // never stored: leave it for on-demand loading
    try { const r = await fetch(f.path, { cache: 'reload' }); if (r.ok) await c.put(f.path, r); } catch { /* next time */ }
  }
  await c.put('precache-state', new Response(JSON.stringify({ version: list.version, hashes: Object.fromEntries(list.files.map((f) => [f.path, f.hash])) })));
}

self.addEventListener('message', (e) => {
  if (e.data?.type === 'cache-all') e.waitUntil(cacheAll(e.ports[0]));
});

// Background refresh while the app is closed (browsers that support periodic background sync).
self.addEventListener('periodicsync', (e) => {
  if (e.tag === 'refresh-data') e.waitUntil(Promise.all([refreshChanged().catch(() => {}), refreshSiteForecast().catch(() => {})]));
});

function idbOpen() { return new Promise((res, rej) => { const r = indexedDB.open('offshore-analyses', 1); r.onsuccess = () => res(r.result); r.onerror = () => rej(r.error); r.onupgradeneeded = () => { r.transaction.abort(); rej(new Error('store not initialised')); }; }); }
const idbGet = (db, store, key) => new Promise((res, rej) => { const q = db.transaction(store).objectStore(store).get(key); q.onsuccess = () => res(q.result); q.onerror = () => rej(q.error); });
async function refreshSiteForecast() {
  const db = await idbOpen();
  const id = await idbGet(db, 'kv', 'currentProject');
  const p = id && await idbGet(db, 'kv', `project:${id}`);
  if (!p?.site || p.site.lat == null) return;
  const r4 = (v) => Math.round(v * 1e4) / 1e4;
  const vars = 'wave_height,wave_peak_period,wave_direction,wind_wave_height,swell_wave_height,swell_wave_peak_period,ocean_current_velocity,ocean_current_direction,sea_level_height_msl,sea_surface_temperature';
  const r = await fetch(`https://marine-api.open-meteo.com/v1/marine?latitude=${p.site.lat}&longitude=${p.site.lon}&current=${vars}&hourly=${vars}&forecast_days=8&timezone=GMT&cell_selection=sea`);
  if (!r.ok) return;
  const data = await r.json();
  await new Promise((res, rej) => { const t = db.transaction('cache', 'readwrite'); t.objectStore('cache').put({ at: Date.now(), data }, `marine:now:${r4(p.site.lat)},${r4(p.site.lon)}`); t.oncomplete = res; t.onerror = () => rej(t.error); });
}
