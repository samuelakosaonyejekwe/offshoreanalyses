// Suite page: its analysis modules, how they cover the specified analyses, and the full specification
// (models, inputs/outputs, initial/boundary conditions, CVV plan, geometry formats, live data).
import { h, icon, tabs, fmt } from '../dom.js';
import { suiteIndex, suiteSpec, suiteModules, coverage, SHORT } from '../../core/registry.js';
import { getProject } from '../../core/store.js';
import { SOURCES, literature } from '../../data/connectors.js';
import { pager } from '../shell.js';

const para = (items) => h('div.prose', items.map((t) => h('p', t)));
const catList = (items) => h('div.prose', items.map((x) => h('p', h('b.g', `${x.cat || x.group}${x.cat || x.group ? '. ' : ''}`), x.text)));

export default async function suite(root, { n: nStr }) {
  const n = Number(nStr);
  const idx = await suiteIndex();
  const meta = idx.find((s) => s.n === n);
  if (!meta) { root.append(h('div.empty', 'No such suite.')); return { title: 'Suite' }; }
  const [spec, mods, cov] = await Promise.all([suiteSpec(n), suiteModules(n), coverage(n)]);
  const runs = getProject().runs;
  const covered = cov.filter((c) => c.modules.length).length;

  const modulesTab = () => {
    if (!mods.length) return h('div.empty', 'The analysis modules of this suite are not installed in this build.');
    return h('div.grid.c2', mods.map((m) => {
      const r = runs[m.id];
      return h('a.card.tile', { href: `#/m/${m.id}` },
        h('div.row.spread', h('h3', m.title), r ? h('span', { class: `chip ${r.status === 'fail' ? 'fail' : 'ok'}` }, r.status === 'fail' ? `✗ UC ${fmt(r.maxUC, 3)}` : '✓ run') : h('span.chip', 'not run')),
        h('p', m.summary),
        h('div.row', { style: { marginTop: '8px' } }, h('span.chip.info', m.fidelity), (m.standards || []).slice(0, 2).map((s) => h('span.chip', s))));
    }));
  };
  const coverageTab = () => h('div.card', h('p.muted', `${covered} of ${cov.length} analyses in the specification have at least one working module. Each line links to the modules that perform it.`),
    h('ul.ticks', cov.map((c) => h('li', { class: c.modules.length ? 'on' : '' }, h('b', `${c.name}. `), c.modules.length ? c.modules.map((m, i) => [i ? ' · ' : '', h('a', { href: `#/m/${m.id}` }, m.title)]) : h('span.muted', 'not yet covered by a module')))));
  const liveTab = () => {
    const box = h('div');
    const srcs = SOURCES.filter((s) => s.suites.includes(n));
    box.append(h('div.card', h('h3', 'Live sources wired to this suite'), h('p.muted', 'Fetched directly by your browser and cached for offline use. Site-dependent values arrive through the project site page.'),
      h('ul.clean', srcs.map((s) => h('li', { style: { padding: '8px 0', borderBottom: '1px solid var(--line)' } }, h('a', { href: s.home, target: '_blank', rel: 'noopener' }, s.name), h('small.muted', ` — ${s.provider}`), h('br'), h('small', s.gives)))),
      h('div.row', { style: { marginTop: '10px' } }, h('a.btn', { href: '#/site' }, icon('globe'), 'Pull site data'), h('a.btn', { href: '#/data' }, 'Source status'))));
    const lit = h('div.card', h('h3', 'Newest publications for this suite'), h('p.muted', 'Loading…'));
    box.append(lit);
    literature(`${meta.title} offshore engineering analysis`.replace(/,| and /g, ' '), { n: 8 }).then((r) => {
      lit.replaceChildren(h('h3', 'Newest publications for this suite'), h('small.muted', `Live from OpenAlex · ${new Date(r.at).toLocaleString()}${r.stale ? ' · cached (offline)' : ''}`),
        h('ul.clean', r.data.map((w) => h('li', { style: { padding: '8px 0', borderBottom: '1px solid var(--line)' } }, h('a', { href: w.oa || w.doi || '#', target: '_blank', rel: 'noopener' }, w.title || 'Untitled'), h('br'), h('small.muted', `${w.authors} · ${w.venue} · ${w.date || w.year}`)))));
    }).catch((e) => lit.replaceChildren(h('h3', 'Newest publications for this suite'), h('p.muted', `Not available right now (${e.message}).`)));
    return box;
  };

  root.append(
    h('div.row', h('span.chip.info', `Suite ${n} of 17`), h('span.chip', `${mods.length} modules`), h('span.chip', `${covered}/${cov.length} analyses covered`)),
    h('h1', { style: { marginTop: '8px' } }, meta.title), h('p.lede', meta.subtitle || ''), meta.scope ? h('p.muted', meta.scope) : null,
    tabs([
      { label: 'Analysis modules', render: modulesTab },
      { label: 'Coverage', render: coverageTab },
      { label: 'Equations & models', render: () => h('div.card', catList(spec.models || [])) },
      { label: 'Inputs & outputs', render: () => h('div.grid.c2', h('div.card', h('h3', 'Input data accepted'), catList(spec.inputs || [])), h('div.card', h('h3', 'Output data generated'), h('ul.ticks', (spec.outputs || []).map((o) => h('li.on', o))))) },
      { label: 'Initial & boundary conditions', render: () => h('div.grid.c2', h('div.card', h('h3', 'Initial conditions'), para(spec.ic || [])), h('div.card', h('h3', n === 17 ? 'Boundary conditions and financial constraints' : 'Boundary conditions'), para(spec.bc || []))) },
      { label: 'Calibration, verification, validation', render: () => h('div', h('div.banner', icon('info'), h('div', 'This is the plan each capability is held to. What has actually been demonstrated in this build is on the ', h('a', { href: '#/verification' }, 'Verification & validation'), ' page — nothing is called validated without evidence.')), h('div.grid.c3', ['calibration', 'verification', 'validation'].map((k) => h('div.card', h('h3', k[0].toUpperCase() + k.slice(1)), para(spec.cvv?.[k] || []))))) },
      ...(n <= 16 ? [{ label: 'Geometry formats', render: () => h('div.card', catList(spec.geometry || []), h('div.row', h('a.btn', { href: '#/geometry' }, icon('cube'), 'Import geometry'))) }] : []),
      { label: 'Live data', render: liveTab },
    ]),
    pager(n > 1 ? { href: `#/suite/${n - 1}`, label: `${n - 1}. ${SHORT[n - 1]}` } : { href: '#/materials', label: 'Materials library' },
      mods.length ? { href: `#/m/${mods[0].id}`, label: mods[0].title } : (n < 17 ? { href: `#/suite/${n + 1}`, label: `${n + 1}. ${SHORT[n + 1]}` } : { href: '#/decision', label: 'Decision support' })));
  return { title: meta.title, crumbs: [{ href: `#/suite/${n}`, label: `${n}. ${SHORT[n]}` }] };
}
