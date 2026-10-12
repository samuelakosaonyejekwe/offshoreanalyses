// Smaller pages: materials library, verification matrix, suite wiring, live-data status, guide, search.
import { h, clear, fmt, icon, table, toast, tabs } from '../dom.js';
import { suiteIndex, suiteModules, coverage, wiring as wiringMap, commonSpec, SHORT } from '../../core/registry.js';
import { getProject, save, idb } from '../../core/store.js';
import { MATERIALS, CATEGORIES, setCustomMaterials } from '../../library/materials.js';
import { SOURCES, health, lastFetch } from '../../data/connectors.js';
import { task } from '../../core/runner.js';
import { pager, installApp } from '../shell.js';

// ---------------------------------------------------------------- materials
const PROPS = [['fy', 'Yield strength', 'MPa', 1e-6], ['fu', 'Tensile strength', 'MPa', 1e-6], ['fck', 'Characteristic compressive strength', 'MPa', 1e-6], ['E', 'Elastic modulus', 'GPa', 1e-9], ['nu', 'Poisson ratio', '-', 1],
  ['rho', 'Density', 'kg/m³', 1], ['alpha', 'Thermal expansion', '1/K', 1], ['k', 'Thermal conductivity', 'W/m·K', 1], ['cp', 'Specific heat', 'J/kg·K', 1], ['su', 'Undrained shear strength', 'kPa', 1e-3], ['suGrad', 'Strength gradient', 'kPa/m', 1e-3],
  ['phiDeg', 'Friction angle', '°', 1], ['gammaSub', 'Submerged unit weight', 'kN/m³', 1e-3], ['d50', 'Median grain size', 'mm', 1e3], ['mu', 'Dynamic viscosity', 'mPa·s', 1e3], ['rho20', 'Resistivity at 20 °C', 'nΩ·m', 1e9], ['epsR', 'Relative permittivity', '-', 1], ['rhoT', 'Thermal resistivity', 'K·m/W', 1]];

export async function materials(root) {
  const p = getProject();
  setCustomMaterials(p.materials);
  const body = h('div');
  const filter = h('input.plain', { type: 'search', placeholder: 'Filter materials…', 'aria-label': 'Filter materials', style: { maxWidth: '280px' } });
  let cat = 'steel';
  const draw = () => {
    clear(body);
    const list = [...MATERIALS, ...p.materials].filter((m) => m.cat === cat && (!filter.value || m.name.toLowerCase().includes(filter.value.toLowerCase())));
    const cols = PROPS.filter(([k]) => list.some((m) => m[k] != null));
    body.append(table([{ key: 'name', label: 'Material', left: true }, ...cols.map(([k, l, u]) => ({ key: k, label: l, unit: u })), { key: 'src', label: 'Source', left: true }],
      list.map((m) => ({ name: p.materials.includes(m) ? h('span', m.name, ' ', h('span.chip.info', 'yours')) : m.name, ...Object.fromEntries(cols.map(([k, , , f]) => [k, m[k] != null ? m[k] * f : null])), src: m.src || '' })), { max: 200 }));
  };
  filter.addEventListener('input', draw);
  const add = () => {
    const f = Object.fromEntries(PROPS.map(([k, l, u]) => [k, h('input.plain', { type: 'number', step: 'any', 'aria-label': l, placeholder: u })]));
    const nm = h('input.plain', { type: 'text', placeholder: 'Name, e.g. “Mill cert heat 4471 — X65”' });
    const dlg = h('dialog', h('h3', `Add your own ${cat} record`), h('p.muted', 'Leave blank what does not apply. Your records are stored with the project and appear in every material picker.'), h('div.field', h('label', 'Name'), nm),
      h('div.fields', { style: { marginTop: '10px' } }, PROPS.map(([k, l, u]) => h('div.field', h('label', `${l} [${u}]`), f[k]))),
      h('div.row', { style: { marginTop: '12px', justifyContent: 'flex-end' } }, h('button.btn', { type: 'button', onclick: () => dlg.close() }, 'Cancel'), h('button.btn.primary', { type: 'button', onclick: () => {
        if (!nm.value.trim()) { toast('Give the material a name.'); return; }
        const rec = { id: `user-${Date.now().toString(36)}`, cat, name: nm.value.trim(), src: 'project-specific' };
        for (const [k, , , fac] of PROPS) if (f[k].value !== '') rec[k] = Number(f[k].value) / fac;
        p.materials.push(rec); save(); setCustomMaterials(p.materials); dlg.close(); draw(); toast('Material added to this project.');
      } }, 'Add')));
    document.body.append(dlg); dlg.addEventListener('close', () => dlg.remove()); dlg.showModal();
  };
  root.append(h('h1', 'Materials library'), h('p.lede', `${MATERIALS.length} built-in records across steels, concretes, mooring components, polymers, conductors, soils, rock and fluids. Modules pick from here; add certified project values and they are offered everywhere.`),
    h('div.row', filter, h('button.btn', { type: 'button', onclick: add }, '+ Add your own')),
    tabs(CATEGORIES.map((c) => ({ label: c[0].toUpperCase() + c.slice(1), render: () => { cat = c; draw(); return body; } }))),
    h('p.muted', 'Strengths are specified minimums from the cited standard; soils, ropes and generic items are typical values and are labelled as such. Certificates and site investigation data take precedence.'),
    pager({ href: '#/geometry', label: 'Geometry & mesh' }, { href: '#/suite/1', label: `1. ${SHORT[1]}` }));
  return { title: 'Materials library', crumbs: [{ href: '#/materials', label: 'Materials library' }] };
}

// ---------------------------------------------------------------- verification
export async function verification(root) {
  const idx = await suiteIndex();
  const common = await commonSpec();
  const out = h('div');
  const summary = h('div.kpis');
  const prog = h('div');
  root.append(h('h1', 'Verification & validation'), h('p.lede', 'What has actually been demonstrated for each capability in this build — and what has not.'),
    h('div.card', h('h3', 'How to read the status'), h('ul.prose',
      h('li', h('b', 'Implemented'), ' — the model is coded with documented inputs, outputs and assumptions.'),
      h('li', h('b', 'Numerically verified'), ' — every built-in benchmark case of the module reproduces an independent closed-form or published value within its tolerance. This shows the equations are solved correctly.'),
      h('li', h('b', 'Calibrated / validated'), ' — only you can establish this, for your asset, with your measurements: each module has a Calibrate & validate tab that keeps calibration and validation data apart and reports error, bias and peak behaviour. Nothing here is labelled experimentally or field validated by default.'))),
    h('div.row', { style: { margin: '12px 0' } }, h('button.btn.primary', { type: 'button', onclick: runAll }, icon('play'), 'Run every verification case now'), h('small.muted', 'Runs on this device; takes a minute or two.')), prog, summary, out,
    h('details.fold', h('summary', 'Common calibration, verification and validation requirements'), h('div.prose', (common.cvv || []).map((t) => h('p', t)))),
    pager({ href: '#/decision', label: 'Decision support' }, { href: '#/wiring', label: 'Suite wiring & coverage' }));

  async function runAll() {
    clear(out); clear(summary);
    const bar = h('i'), msg = h('small.muted');
    clear(prog).append(h('div.progress', bar), msg);
    const all = [];
    for (const s of idx) for (const m of await suiteModules(s.n)) all.push([s, m]);
    let done = 0, pass = 0, fail = 0, pending = 0, cases = 0;
    const bySuite = new Map();
    for (const [s, m] of all) {
      msg.textContent = `Suite ${s.n}: ${m.title}`;
      let res;
      try { res = await task('verify', m.id, {}); } catch (e) { res = [{ name: 'Module failed to run', status: 'fail', error: e.message, rows: [] }]; }
      const st = !res.length || res.every((c) => c.status === 'pending') ? 'pending' : res.some((c) => c.status === 'fail') ? 'fail' : 'pass';
      cases += res.filter((c) => c.status !== 'pending').length;
      if (st === 'pass') pass++; else if (st === 'fail') fail++; else pending++;
      if (!bySuite.has(s.n)) bySuite.set(s.n, []);
      bySuite.get(s.n).push({ m, res, st });
      bar.style.width = `${Math.round(100 * ++done / all.length)}%`;
    }
    clear(prog);
    summary.append(h('div.kpi.ok', h('span', 'Modules numerically verified'), h('b', pass), h('small.muted', `of ${all.length}`)), h('div.kpi', { class: fail ? 'fail' : '' }, h('span', 'Modules with a failing case'), h('b', fail)),
      h('div.kpi', h('span', 'Implemented, no independent benchmark yet'), h('b', pending)), h('div.kpi', h('span', 'Benchmark cases executed'), h('b', cases)));
    await idb.set('kv', 'verification', { at: Date.now(), pass, fail, pending, total: all.length, cases });
    for (const s of idx) {
      const rows = bySuite.get(s.n) || [];
      out.append(h('details.fold', { open: rows.some((r) => r.st === 'fail') }, h('summary', `${s.n}. ${SHORT[s.n]} — ${rows.filter((r) => r.st === 'pass').length}/${rows.length} verified`),
        h('div', table([{ key: 'title', label: 'Module', left: true }, { key: 'status', label: 'Status', left: true }, { key: 'cases', label: 'Benchmark cases', left: true }],
          rows.map((r) => ({ title: h('a', { href: `#/m/${r.m.id}` }, r.m.title), status: h('span', { class: `chip ${r.st === 'pass' ? 'ok' : r.st === 'fail' ? 'fail' : 'warn'}` }, r.st === 'pass' ? '✓ numerically verified' : r.st === 'fail' ? '✗ a case fails' : 'implemented'),
            cases: h('span', r.res.map((c) => h('div', h('small', `${c.status === 'pass' ? '✓' : c.status === 'pending' ? '·' : '✗'} ${c.name}${c.error ? ` — ${c.error}` : ''}${c.rows?.some((x) => !x.pass) ? ` — got ${c.rows.filter((x) => !x.pass).map((x) => `${fmt(x.got)} vs ${fmt(x.expected)}`).join(', ')}` : ''}`)))) }))))));
    }
  }
  const prev = await idb.get('kv', 'verification');
  if (prev) summary.append(h('div.kpi.ok', h('span', `Last full run ${new Date(prev.at).toLocaleDateString()}`), h('b', `${prev.pass}/${prev.total}`), h('small.muted', 'modules numerically verified')), h('div.kpi', { class: prev.fail ? 'fail' : '' }, h('span', 'Failing'), h('b', prev.fail)), h('div.kpi', h('span', 'Benchmark cases'), h('b', prev.cases)));
  return { title: 'Verification & validation', crumbs: [{ href: '#/verification', label: 'Verification & validation' }] };
}

// ---------------------------------------------------------------- wiring & coverage
export async function wiring(root) {
  const idx = await suiteIndex();
  const [map, common] = await Promise.all([wiringMap(), commonSpec()]);
  const bus = getProject().bus;
  const cov = await Promise.all(idx.map((s) => coverage(s.n)));
  const total = cov.reduce((s, c) => s + c.length, 0), done = cov.reduce((s, c) => s + c.filter((x) => x.modules.length).length, 0);
  const mods = await Promise.all(idx.map((s) => suiteModules(s.n)));
  // suite-to-suite links: a publisher in one suite feeding a consumer in another
  const links = new Map();
  for (const e of map) for (const a of e.publishers) for (const b of e.consumers) if (a.suite !== b.suite) { const k = `${a.suite}→${b.suite}`; if (!links.has(k)) links.set(k, new Set()); links.get(k).add(e.key); }
  const matrix = idx.map((a) => idx.map((b) => links.get(`${a.n}→${b.n}`)?.size || 0));

  root.append(h('h1', 'Suite wiring & coverage'), h('p.lede', 'The 17 suites are connected through named project values: what one analysis computes becomes the default input of the analyses that depend on it — loads, motions, reactions, stiffnesses, quantities, energy and cost.'),
    h('div.kpis', h('div.kpi', h('span', 'Analysis modules'), h('b', mods.reduce((s, m) => s + m.length, 0))), h('div.kpi', h('span', 'Specified analyses covered'), h('b', `${done}/${total}`)), h('div.kpi', h('span', 'Shared project values'), h('b', map.length)), h('div.kpi', h('span', 'Suite-to-suite links'), h('b', links.size))),
    h('h2', { style: { marginTop: '20px' } }, 'Who feeds whom'),
    (await import('../plot.js')).renderPlot({ type: 'heatmap', title: 'Values passed from suite (row) to suite (column)', xs: idx.map((s) => s.n), ys: idx.map((s) => s.n), z: matrix, x: { label: 'Receiving suite' }, y: { label: 'Providing suite' }, zLabel: 'Shared values' }),
    h('h2', { style: { marginTop: '20px' } }, 'Project values'),
    table([{ key: 'key', label: 'Value', left: true }, { key: 'now', label: 'Current (SI)' }, { key: 'pub', label: 'Computed by', left: true }, { key: 'con', label: 'Used by', left: true }],
      map.map((e) => ({ key: h('code', e.key), now: bus[e.key]?.value ?? null, pub: h('span', e.key.startsWith('site.') && !e.publishers.length ? h('a', { href: '#/site' }, 'site data') : e.publishers.map((x, i) => [i ? ', ' : '', h('a', { href: `#/m/${x.id}` }, `${x.suite}: ${x.title}`)])), con: h('span', e.consumers.slice(0, 6).map((x, i) => [i ? ', ' : '', h('a', { href: `#/m/${x.id}` }, `${x.suite}: ${x.title}`)]), e.consumers.length > 6 ? ` +${e.consumers.length - 6} more` : '') })), { max: 600 }),
    h('h2', { style: { marginTop: '20px' } }, 'Coverage of the specified analyses'),
    idx.map((s, i) => h('details.fold', h('summary', `${s.n}. ${SHORT[s.n]} — ${cov[i].filter((c) => c.modules.length).length}/${cov[i].length}`), h('div', h('ul.ticks', cov[i].map((c) => h('li', { class: c.modules.length ? 'on' : '' }, `${c.name} `, c.modules.map((m) => h('a', { href: `#/m/${m.id}`, style: { marginRight: '8px' } }, `→ ${m.title}`)))))))),
    h('details.fold', h('summary', 'Integration requirements from the specification'), h('div.prose', (common.integration || []).slice(0, 5).map((t) => h('p', t)))),
    pager({ href: '#/verification', label: 'Verification & validation' }, { href: '#/data', label: 'Live data sources' }));
  return { title: 'Suite wiring & coverage', crumbs: [{ href: '#/wiring', label: 'Suite wiring & coverage' }] };
}

// ---------------------------------------------------------------- live data
export async function data(root) {
  const body = h('div');
  const draw = async (probe) => {
    clear(body).append(h('small.muted', probe ? 'Contacting every source…' : 'Reading cache…'));
    const hl = probe ? await health() : {};
    const rows = [];
    for (const s of SOURCES) {
      const at = await lastFetch(s.id);
      const st = hl[s.id];
      rows.push({ name: h('span', h('a', { href: s.home, target: '_blank', rel: 'noopener' }, s.name), h('br'), h('small.muted', s.provider)), gives: s.gives,
        status: st ? h('span', { class: `chip ${st.ok ? 'ok' : 'fail'}` }, st.ok ? `✓ reachable · ${st.ms} ms` : `✗ ${st.error || 'no answer'}`) : h('span.chip', probe ? 'on demand' : '—'),
        last: at ? new Date(at).toLocaleString() : 'never', suites: s.suites.length === 17 ? 'all' : s.suites.join(', ') || 'app', lic: s.licence });
    }
    clear(body).append(table([{ key: 'name', label: 'Source', left: true }, { key: 'gives', label: 'What it provides', left: true }, { key: 'status', label: 'Right now', left: true }, { key: 'last', label: 'Last fetched here', left: true }, { key: 'suites', label: 'Suites', left: true }, { key: 'lic', label: 'Licence', left: true }], rows));
  };
  root.append(h('h1', 'Live data sources'), h('p.lede', 'Every source is contacted directly by the browser on the device you are holding. Nothing passes through a server of ours, so data stays current wherever the app is opened and whoever opens it. Offline, the last answers are served from the device cache and marked as cached.'),
    h('div.row', { style: { marginBottom: '12px' } }, h('button.btn.primary', { type: 'button', onclick: () => draw(true) }, icon('refresh'), 'Check all sources now'), h('span', { class: `chip ${navigator.onLine ? 'ok' : 'warn'}` }, navigator.onLine ? 'online' : 'offline')), body,
    h('div.card', { style: { marginTop: '14px' } }, h('h3', 'How freshness works'), h('ul.prose',
      h('li', 'Each answer is stored with its fetch time and a lifetime that matches how fast the source changes — an hour for forecasts and grid intensity, a day for prices, months for hindcasts and bathymetry.'),
      h('li', 'When the lifetime has passed and the device is online, the next use fetches again. “Refresh” on the site page forces it.'),
      h('li', 'Installed as an app, supporting browsers also wake it in the background to keep the cache warm.'),
      h('li', 'Engineering judgement still applies: open global models are for screening and early design. Measured site data entered on the site page always overrides them.'))),
    pager({ href: '#/wiring', label: 'Suite wiring & coverage' }, { href: '#/help', label: 'Guide & install' }));
  draw(false);
  return { title: 'Live data sources', crumbs: [{ href: '#/data', label: 'Live data sources' }] };
}

// ---------------------------------------------------------------- guide & install
export async function help(root) {
  const mirrors = await fetch(new URL('../../../mirrors.json', import.meta.url)).then((r) => r.json()).catch(() => ({ mirrors: [] }));
  const offline = h('div');
  const cacheAll = async () => {
    const reg = await navigator.serviceWorker?.ready.catch(() => null);
    if (!reg?.active) { toast('Offline storage is not available in this browser or on this address.'); return; }
    const bar = h('i'), msg = h('small.muted', 'Starting…');
    clear(offline).append(h('div.progress', bar), msg);
    const ch = new MessageChannel();
    ch.port1.onmessage = ({ data: d }) => {
      if (d.progress != null) { bar.style.width = `${Math.round(d.progress * 100)}%`; msg.textContent = `${d.done} of ${d.total} files stored`; }
      if (d.done === d.total) { msg.textContent = `All ${d.total} files are stored on this device${d.failed ? ` (${d.failed} could not be fetched)` : ''}. The app now runs fully offline.`; }
    };
    reg.active.postMessage({ type: 'cache-all' }, [ch.port2]);
  };
  const est = navigator.storage?.estimate ? await navigator.storage.estimate().catch(() => null) : null;
  root.append(h('h1', 'Guide & install'),
    h('div.card', h('h3', 'Install on this device'), h('p', 'Installed, the app opens from your home screen or desktop like any other, in its own window, and works without a connection.'),
      h('div.row', h('button.btn.primary', { type: 'button', onclick: installApp }, icon('install'), 'Install now'), h('button.btn', { type: 'button', onclick: cacheAll }, icon('download'), 'Store everything for offline use')), offline,
      h('div.grid.c3', { style: { marginTop: '12px' } },
        h('div', h('h4', 'Windows, macOS, Linux, ChromeOS'), h('p.muted', 'Chrome or Edge: press “Install now”, or use the install icon at the right of the address bar. Safari on macOS: File → Add to Dock.')),
        h('div', h('h4', 'Android'), h('p.muted', 'Chrome, Edge, Samsung Internet or Firefox: press “Install now”, or open the browser menu and choose “Install app” / “Add to Home screen”.')),
        h('div', h('h4', 'iPhone and iPad'), h('p.muted', 'Safari: tap the Share button, then “Add to Home Screen”. (Apple does not let a web page trigger this, so the button cannot do it for you.)'))),
      est ? h('small.muted', `Stored on this device: ${fmt((est.usage || 0) / 1e6, 3)} MB of ${fmt((est.quota || 0) / 1e9, 3)} GB available.`) : null),
    h('div.card', h('h3', 'Working offline and in flight mode'), h('ul.prose',
      h('li', 'After the first visit the application shell is on the device. Press “Store everything for offline use” once to keep every suite, the geometry readers and the specification data too.'),
      h('li', 'All calculations run on the device. Projects, imported geometry and results are saved in the browser’s own database — nothing is uploaded.'),
      h('li', 'Live data shows the last fetched values with their date while offline, and refreshes by itself when the connection returns.'))),
    h('div.card', h('h3', 'Always reachable'), h('p', 'The same application is published at more than one independent address. If one host is down, open another — they are identical, and an installed copy keeps working regardless.'),
      mirrors.mirrors?.length ? h('ul.prose', mirrors.mirrors.map((m) => h('li', h('a', { href: m.url, target: '_blank', rel: 'noopener' }, m.url), ` — ${m.host}${m.note ? ` (${m.note})` : ''}`))) : h('p.muted', 'No mirror list is published with this copy.')),
    h('div.card', h('h3', 'A typical session'), h('ol.prose',
      h('li', h('b', 'Site. '), 'Pick the location; design waves, wind, current, depth and seismicity are derived and shared with every suite. Override anything with measured data.'),
      h('li', h('b', 'Data and geometry. '), 'Drop CAD, mesh, route, survey or table files into the portal. They are recognised, checked, repaired and meshed, and become selectable inside the analyses that can use them.'),
      h('li', h('b', 'Analyse. '), 'Open a suite, open a module, edit the inputs. Results update as you type for fast analyses. Linked fields show where a value came from.'),
      h('li', h('b', 'Explore. '), 'Rank the inputs by influence, sweep a parameter, propagate uncertainty, and prove discretisation convergence with an automatic mesh selection.'),
      h('li', h('b', 'Calibrate and validate. '), 'Paste measurements to tune uncertain parameters inside bounds, then test against independent data.'),
      h('li', h('b', 'Decide. '), 'The decision page ranks what to fix, shows where reserve can be given back, and puts cost and carbon next to the engineering result.'))),
    h('div.card', h('h3', 'Keyboard'), h('p.muted', 'Alt + ← / → back and forward · ↑ / ↓ in a number field nudges it by 1 % (10 % with Shift) · in a 3-D view: arrows orbit, + / − zoom, F fits.')),
    h('div.card', h('h3', 'Engineering responsibility'), h('p', 'This software performs calculations; it does not replace a qualified engineer or an independent design verification. Every module states its method, its standard and its limits on its Method tab, and its verification status on the Verification page. Use results for design decisions only inside those limits and after review.')),
    pager({ href: '#/data', label: 'Live data sources' }, { href: '#/', label: 'Home' }));
  return { title: 'Guide & install', crumbs: [{ href: '#/help', label: 'Guide & install' }] };
}

// ---------------------------------------------------------------- search
export async function search(root, { q }) {
  const idx = await suiteIndex();
  const words = q.toLowerCase().split(/\s+/).filter(Boolean);
  const hits = [];
  for (const s of idx) {
    for (const m of await suiteModules(s.n)) {
      const hay = `${m.title} ${m.summary} ${(m.models || []).join(' ')} ${(m.standards || []).join(' ')} ${(m.inputs || []).map((i) => i.label).join(' ')}`.toLowerCase();
      const score = words.reduce((sc, w) => sc + (m.title.toLowerCase().includes(w) ? 3 : 0) + (hay.includes(w) ? 1 : 0), 0);
      if (words.every((w) => hay.includes(w))) hits.push({ m, s, score });
    }
    s.analyses.forEach((a, i) => { if (words.every((w) => a.toLowerCase().includes(w))) hits.push({ analysis: a, s, score: 2, i }); });
  }
  hits.sort((a, b) => b.score - a.score);
  root.append(h('h1', `Search: ${q}`), hits.length ? h('div.grid.c2', hits.slice(0, 60).map((x) => (x.m
    ? h('a.card.tile', { href: `#/m/${x.m.id}` }, h('small.muted', `${x.s.n}. ${SHORT[x.s.n]}`), h('h3', x.m.title), h('p', x.m.summary))
    : h('a.card.tile', { href: `#/suite/${x.s.n}` }, h('small.muted', `${x.s.n}. ${SHORT[x.s.n]} — specified analysis`), h('h3', x.analysis))))) : h('div.empty', 'Nothing matches. Try a shorter term such as “fatigue”, “catenary”, “burst” or “LCOE”.'));
  return { title: `Search: ${q}`, crumbs: [{ href: `#/search/${encodeURIComponent(q)}`, label: 'Search' }] };
}
