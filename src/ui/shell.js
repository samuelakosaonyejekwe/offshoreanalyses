// Application shell: header with back / forward arrows, breadcrumbs, navigation drawer, hash router.
import { h, clear, icon, toast } from './dom.js';
import { suiteIndex, SHORT } from '../core/registry.js';
import { getProject, on, prefs } from '../core/store.js';

const routes = [];
export function route(pattern, loader) { routes.push({ re: new RegExp(`^${pattern.replace(/:(\w+)/g, '(?<$1>[^/]+)')}$`), loader }); }

// The browser does not expose whether "forward" is possible, so the shell keeps its own position in the
// session history: every entry carries an index in history.state.
let pos = 0, maxPos = 0;
const els = {};

export function pager(prev, next) {
  return h('nav.pager', { 'aria-label': 'Previous and next page' },
    prev ? h('a', { href: prev.href, rel: 'prev' }, h('span.arrow', icon('back')), h('div', h('small.muted', 'Previous'), h('b', prev.label))) : h('span'),
    next ? h('a.next', { href: next.href, rel: 'next' }, h('div', h('small.muted', 'Next'), h('b', next.label)), h('span.arrow', icon('fwd'))) : null);
}

const NAV = [
  ['Workspace', [['#/', 'home', 'Home'], ['#/site', 'globe', 'Project & site data'], ['#/portal', 'upload', 'Input data portal'], ['#/geometry', 'cube', 'Geometry & mesh'], ['#/materials', 'layers', 'Materials library']]],
  ['suites'],
  ['Assurance & decisions', [['#/decision', 'bulb', 'Decision support'], ['#/verification', 'check', 'Verification & validation'], ['#/factors', 'book', 'Code-factor register'], ['#/wiring', 'link', 'Suite wiring & coverage'], ['#/data', 'refresh', 'Live data sources'], ['#/help', 'help', 'Guide & install']]],
];

async function buildNav() {
  const idx = await suiteIndex();
  clear(els.nav);
  const search = h('input.plain', { type: 'search', placeholder: 'Find an analysis…', 'aria-label': 'Find an analysis', onkeydown: (e) => { if (e.key === 'Enter' && search.value.trim()) { location.hash = `#/search/${encodeURIComponent(search.value.trim())}`; } } });
  els.nav.append(h('div.navsearch', search));
  for (const [title, items] of NAV) {
    if (title === 'suites') {
      els.nav.append(h('h4', 'The 17 suites'));
      for (const s of idx) els.nav.append(h('a', { href: `#/suite/${s.n}`, title: s.title }, h('span.n', s.n), SHORT[s.n]));
      continue;
    }
    els.nav.append(h('h4', title));
    for (const [href, ic, label] of items) els.nav.append(h('a', { href }, icon(ic), label));
  }
  markCurrent();
}
function markCurrent() {
  const hash = location.hash || '#/';
  const suite = /^#\/m\/s(\d\d)/.exec(hash)?.[1];
  for (const a of els.nav.querySelectorAll('a')) {
    const href = a.getAttribute('href');
    const cur = href === hash || (suite && href === `#/suite/${Number(suite)}`) || (href !== '#/' && hash.startsWith(`${href}/`));
    if (cur) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current');
  }
}

function setCrumbs(crumbs = []) {
  clear(els.crumbs);
  els.crumbs.append(h('a', { href: '#/' }, getProject().name));
  for (const c of crumbs) els.crumbs.append(h('span.sep', '›'), h('a', { href: c.href }, c.label));
}

let navToken = 0;
async function render() {
  const token = ++navToken;
  const path = (location.hash || '#/').slice(1) || '/';
  if (history.state?.i == null) history.replaceState({ i: pos }, '');
  else pos = history.state.i;
  maxPos = Math.max(maxPos, pos);
  els.back.disabled = pos <= 0;
  els.fwd.disabled = pos >= maxPos;
  closeNav();
  const view = h('div');
  let meta = {};
  try {
    let matched = false;
    for (const r of routes) {
      const m = r.re.exec(path);
      if (!m) continue;
      matched = true;
      const mod = await r.loader();
      if (token !== navToken) return;
      meta = (await mod.default(view, Object.fromEntries(Object.entries(m.groups || {}).map(([k, v]) => [k, decodeURIComponent(v)])))) || {};
      break;
    }
    if (!matched) view.append(h('div.empty', h('h2', 'Page not found'), h('a.btn', { href: '#/' }, 'Go to the home page')));
  } catch (err) {
    console.error(err);
    view.append(h('div.banner.fail', icon('warn'), h('div', h('b', 'This page could not be shown. '), err.message)));
  }
  if (token !== navToken) return;
  clear(els.main).append(view);
  setCrumbs(meta.crumbs);
  document.title = `${meta.title ? `${meta.title} — ` : ''}Offshore Analyses`;
  markCurrent();
  window.scrollTo(0, 0);
  els.main.focus({ preventScroll: true });
}

const closeNav = () => { els.nav.classList.remove('open'); els.scrim.classList.remove('open'); };

export function applyTheme(t) {
  if (t === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', t);
  window.dispatchEvent(new Event('themechange'));
}

let installEvt = null;
export function installApp() {
  if (installEvt) { installEvt.prompt(); installEvt.userChoice.finally(() => { installEvt = null; }); return; }
  location.hash = '#/help';
  toast('Follow the install steps for your device on this page.');
}

export function startShell() {
  els.main = document.getElementById('main');
  els.nav = document.getElementById('nav');
  els.crumbs = document.getElementById('crumbs');
  els.back = document.getElementById('nav-back');
  els.fwd = document.getElementById('nav-fwd');
  els.scrim = document.getElementById('scrim');
  els.net = document.getElementById('net');

  els.back.addEventListener('click', () => history.back());
  els.fwd.addEventListener('click', () => history.forward());
  document.getElementById('menu-btn').addEventListener('click', () => { els.nav.classList.toggle('open'); els.scrim.classList.toggle('open'); });
  els.scrim.addEventListener('click', closeNav);
  document.getElementById('install-btn').addEventListener('click', installApp);
  const themes = ['auto', 'light', 'dark'];
  document.getElementById('theme-btn').addEventListener('click', () => { const t = themes[(themes.indexOf(prefs.get('theme', 'auto')) + 1) % 3]; prefs.set('theme', t); applyTheme(t); toast(`Theme: ${t === 'auto' ? 'follow device' : t}`); });
  applyTheme(prefs.get('theme', 'auto'));

  // a link click creates a new history entry one past the current position; anything beyond it is discarded
  window.addEventListener('hashchange', () => {
    if (history.state?.i == null) { pos += 1; maxPos = pos; history.replaceState({ i: pos }, ''); }
    render();
  });
  window.addEventListener('keydown', (e) => {
    if (e.altKey && e.key === 'ArrowLeft') history.back();
    if (e.altKey && e.key === 'ArrowRight') history.forward();
  });
  const net = () => { els.net.classList.toggle('off', !navigator.onLine); els.net.lastChild.textContent = navigator.onLine ? 'Online — live data' : 'Offline — cached data'; };
  window.addEventListener('online', net); window.addEventListener('offline', net); net();
  window.addEventListener('beforeinstallprompt', (e) => { e.preventDefault(); installEvt = e; });
  window.addEventListener('appinstalled', () => toast('Offshore Analyses is installed on this device.'));
  on('project', () => { setCrumbs(); });
  buildNav();
  render();
}
