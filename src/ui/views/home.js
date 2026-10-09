// Home: what the app is, where the project stands, and the 17 suites.
import { h, fmt, icon } from '../dom.js';
import { suiteIndex, suiteModules, SHORT } from '../../core/registry.js';
import { getProject } from '../../core/store.js';
import { pager, installApp } from '../shell.js';

export default async function home(root) {
  const p = getProject();
  const idx = await suiteIndex();
  const runs = Object.values(p.runs);
  const failing = runs.filter((r) => r.status === 'fail').length;
  const site = p.site.pulled ? `${p.site.name} · ${p.bus['site.depth'] ? `${fmt(p.bus['site.depth'].value, 3)} m water depth` : 'on land'}` : 'No site chosen yet';

  root.append(h('section.hero',
    h('h1', 'Analyse any offshore, marine or subsea structure — from one place'),
    h('p', 'Seventeen connected suites cover fixed and floating structures, foundations, marine terminals, pipelines, risers, moorings, cables, installation, offshore renewables and project economics. Pick a site anywhere on Earth, feed in your geometry and data, run the analyses, and get checks, charts and recommendations — online or offline.'),
    h('div.row', h('a.btn', { href: '#/site' }, icon('globe'), p.site.pulled ? 'Review site data' : 'Start with your site'), h('a.btn.ghost', { href: '#/portal' }, icon('upload'), 'Feed in data & geometry'), h('button.btn.ghost', { type: 'button', onclick: installApp }, icon('install'), 'Install the app')),
    h('svg.wave', { viewBox: '0 0 420 160', 'aria-hidden': 'true' }, h('path', { d: 'M0 90c35-40 70-40 105 0s70 40 105 0 70-40 105 0 70 40 105 0v70H0z', fill: '#fff' }))));

  root.append(h('div.grid.c4', { style: { margin: '16px 0' } },
    h('a.card.tile', { href: '#/site' }, h('div.stat', h('span', 'Project'), h('b', { style: { fontSize: '1.05rem' } }, p.name), h('span', site))),
    h('a.card.tile', { href: '#/decision' }, h('div.stat', h('span', 'Analyses run'), h('b', runs.length), h('span', failing ? `${failing} with a criterion not met` : runs.length ? 'all criteria met' : 'none yet'))),
    h('a.card.tile', { href: '#/geometry' }, h('div.stat', h('span', 'Imported geometry'), h('b', p.geometry.length), h('span', 'files in this project'))),
    h('a.card.tile', { href: '#/wiring' }, h('div.stat', h('span', 'Shared project values'), h('b', Object.keys(p.bus).length), h('span', 'passed between suites')))));

  root.append(h('h2', { style: { marginTop: '22px' } }, 'How it works'),
    h('div.grid.c4',
      [['1', 'Locate', 'Choose the site on the map or by name. Waves, wind, currents, tide, depth, seismicity and solar resource are pulled live and reduced to design values.', '#/site'],
        ['2', 'Feed in', 'Upload CAD, mesh, route, bathymetry or tabular data in the portal. Geometry is checked, repaired and meshed for you.', '#/portal'],
        ['3', 'Analyse', 'Open a suite and run its modules. Results pass to the suites that depend on them automatically.', '#/suite/1'],
        ['4', 'Decide', 'See every check ranked, the smallest design change that fixes a failure, and the cost and carbon consequence.', '#/decision']]
        .map(([n, t, d, href]) => h('a.card.tile', { href }, h('span.n', n), h('h3', t), h('p', d)))));

  const grid = h('div.grid.c3');
  root.append(h('h2', { style: { marginTop: '26px' } }, 'The 17 suites'), grid);
  for (const s of idx) {
    const count = h('small.muted', ' ');
    grid.append(h('a.card.tile', { href: `#/suite/${s.n}` }, h('span.n', s.n), h('h3', SHORT[s.n]), h('p', s.subtitle || s.title), h('p', { style: { marginTop: '8px' } }, count)));
    suiteModules(s.n).then((mods) => {
      const done = mods.filter((m) => p.runs[m.id]).length;
      count.textContent = mods.length ? `${mods.length} analysis modules · ${s.analyses.length} specified analyses${done ? ` · ${done} run` : ''}` : 'modules loading…';
    });
  }
  root.append(pager(null, { href: '#/site', label: 'Project & site data' }));
  return { title: 'Home' };
}
