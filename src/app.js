// Entry point: boot the project store, register routes (each view is loaded on first visit), start the shell.
import { boot } from './core/store.js';
import { route, startShell } from './ui/shell.js';
import { toast } from './ui/dom.js';

const view = (file, fn = 'default') => () => import(`./ui/views/${file}.js`).then((m) => ({ default: m[fn] }));
route('/', view('home'));
route('/site', view('site'));
route('/portal', view('portal'));
route('/geometry', view('geometry'));
route('/geometry/:id', view('geometry'));
route('/materials', view('pages', 'materials'));
route('/suite/:n', view('suite'));
route('/m/:id', () => import('./ui/views/module.js').then((m) => ({ default: (root, p) => m.moduleView(root, p.id) })));
route('/decision', view('decision'));
route('/verification', view('pages', 'verification'));
route('/wiring', view('pages', 'wiring'));
route('/data', view('pages', 'data'));
route('/help', view('pages', 'help'));
route('/search/:q', view('pages', 'search'));

async function registerServiceWorker() {
  if (!('serviceWorker' in navigator) || location.protocol === 'file:') return;
  try {
    const reg = await navigator.serviceWorker.register('sw.js');
    reg.addEventListener('updatefound', () => {
      const sw = reg.installing;
      sw?.addEventListener('statechange', () => { if (sw.state === 'installed' && navigator.serviceWorker.controller) toast('A new version is ready — it loads the next time you open the app.', 6000); });
    });
    // ask the browser to wake the worker periodically so cached data stays fresh while the app is closed (where supported)
    if ('periodicSync' in reg) {
      const perm = await navigator.permissions.query({ name: 'periodic-background-sync' }).catch(() => null);
      if (perm?.state === 'granted') await reg.periodicSync.register('refresh-data', { minInterval: 6 * 3600e3 }).catch(() => {});
    }
    setInterval(() => reg.update().catch(() => {}), 30 * 60e3);
  } catch (err) { console.warn('Offline support could not be enabled:', err.message); }
}

boot().catch((e) => console.warn('Project store unavailable, using memory only:', e)).finally(() => {
  startShell();
  registerServiceWorker();
  if (navigator.storage?.persist) navigator.storage.persist().catch(() => {});
});
