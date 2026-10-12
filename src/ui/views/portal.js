// Input data portal: one place to feed a case in. Whatever is dropped here is recognised and routed —
// geometry to the geometry bench, saved projects and run records back into the project, tables of site and
// design values onto the project bus where every suite picks them up.
import { h, clear, icon, toast, table, download, fmt } from '../dom.js';
import { getProject, importProject, setInputs, publish, save, busValues } from '../../core/store.js';
import { suiteIndex, suiteSpec, suiteModules, findModule, SHORT } from '../../core/registry.js';
import { SITE_KEYS } from '../../data/site.js';
import { importFiles, dropZone } from './geometry.js';
import { pager } from '../shell.js';

const TEMPLATES = {
  'site-values.csv': `key,value,unit,note\n${SITE_KEYS.map(([k, l, u]) => `${k},,${u},${l}`).join('\n')}\n`,
  'route-centreline.csv': 'KP,easting,northing,depth,OD,WT,burial\n0.000,455000.0,6250000.0,92.1,0.3239,0.0159,0.0\n0.500,455480.2,6250139.6,93.4,0.3239,0.0159,0.0\n1.000,455960.3,6250279.3,95.0,0.3239,0.0159,1.0\n',
  'frame-nodes-members.csv': '[NODES]\nid,x,y,z\n1,0,0,-90\n2,20,0,-90\n3,2,0,10\n4,18,0,10\n\n[MEMBERS]\nid,a,b,D,t,material,group\n1,1,3,1.2,0.040,s355,leg\n2,2,4,1.2,0.040,s355,leg\n3,3,4,0.6,0.020,s355,brace\n',
  'soil-layers.csv': 'top,thickness,type,gammaSub,su,suGrad,phi\n0,8,clay,7.0,20,1.5,\n8,15,sand,9.5,,,35\n23,40,clay,9.0,120,2.0,\n',
  'wave-scatter.csv': 'Hs,Tp,probability\n0.5,6,0.21\n1.5,7,0.34\n2.5,8,0.24\n3.5,9,0.13\n4.5,10,0.06\n6.0,11,0.02\n',
  'time-series.csv': 'time_s,value\n0,0.00\n0.5,0.12\n1.0,0.31\n',
};

function parseKeyValues(text) {
  const out = {};
  for (const line of text.split(/\r?\n/)) {
    const m = /^\s*([A-Za-z][\w.]*)\s*[,=:\t;]\s*(-?[\d.]+(?:[eE][-+]?\d+)?)/.exec(line);
    if (m) out[m[1]] = Number(m[2]);
  }
  return out;
}

export default async function portal(root) {
  const log = h('div');
  const say = (ok, text, link) => log.prepend(h('div.banner', { class: ok ? '' : 'fail' }, icon(ok ? 'check' : 'warn'), h('div', text, link ? [' ', link] : null)));

  async function route(files) {
    const geo = [];
    for (const f of [...files]) {
      const name = f.name.toLowerCase();
      if (name.endsWith('.json')) {
        let obj = null;
        try { obj = JSON.parse(await f.text()); } catch { /* not JSON after all: let the geometry reader try */ }
        if (obj && obj.application === 'Offshore Analyses' && obj.module?.id && obj.inputs) {
          const mod = await findModule(obj.module.id);
          if (!mod) { say(false, `${f.name}: run record for “${obj.module.title}”, which is not installed in this build.`); continue; }
          setInputs(mod.id, obj.inputs);
          say(true, `${f.name}: inputs of “${mod.title}” restored.`, h('a', { href: `#/m/${mod.id}` }, 'Open the analysis →'));
          continue;
        }
        if (obj && obj.bus && obj.inputs && obj.site) { importProject(obj); say(true, `${f.name}: project “${obj.name}” imported and opened.`, h('a', { href: '#/site' }, 'Review →')); continue; }
        if (obj && !Array.isArray(obj) && Object.keys(obj).every((k) => /^s\d\d-/.test(k))) {
          let n = 0;
          for (const [id, vals] of Object.entries(obj)) if (await findModule(id)) { setInputs(id, { ...(getProject().inputs[id] || {}), ...vals }); n++; }
          say(n > 0, `${f.name}: inputs set for ${n} analysis module${n === 1 ? '' : 's'}.`);
          continue;
        }
        if (obj && !Array.isArray(obj) && Object.values(obj).every((v) => typeof v === 'number') && Object.keys(obj).some((k) => k.includes('.'))) {
          publish(obj, 'user'); say(true, `${f.name}: ${Object.keys(obj).length} project values set.`, h('a', { href: '#/wiring' }, 'See where they are used →')); continue;
        }
      }
      if (/\.(csv|txt)$/.test(name)) {
        const text = await f.text();
        const kv = parseKeyValues(text);
        const keys = Object.keys(kv).filter((k) => k.includes('.'));
        if (keys.length && keys.length >= text.trim().split(/\r?\n/).length * 0.5 - 1) {
          publish(Object.fromEntries(keys.map((k) => [k, kv[k]])), 'user');
          say(true, `${f.name}: ${keys.length} project values set (${keys.slice(0, 4).join(', ')}${keys.length > 4 ? '…' : ''}).`, h('a', { href: '#/site' }, 'Review →'));
          continue;
        }
      }
      geo.push(f);
    }
    if (geo.length) {
      const res = await importFiles(geo);
      for (const r of res) say(r.ok, r.ok ? `${r.file}: geometry imported (${r.models} model${r.models > 1 ? 's' : ''}).` : `${r.file}: ${r.error}`, r.ok ? h('a', { href: `#/geometry/${r.id}` }, 'Open in the geometry bench →') : null);
    }
  }

  const paste = h('textarea.plain', { rows: 6, placeholder: 'site.depth = 92\nsite.Hs100 = 13.4\nsite.Tp100 = 15.8\nsite.su0 = 20000', 'aria-label': 'Project values as key = value lines' });
  const applyPaste = () => { const kv = parseKeyValues(paste.value); const n = Object.keys(kv).length; if (!n) { toast('No “key = value” lines found.'); return; } publish(kv, 'user'); toast(`${n} project value${n > 1 ? 's' : ''} set.`); paste.value = ''; drawBus(); };

  const busBox = h('div');
  const drawBus = () => {
    const bus = getProject().bus; const keys = Object.keys(bus).sort();
    clear(busBox).append(keys.length ? table([{ key: 'k', label: 'Project value', left: true }, { key: 'v', label: 'Value (SI)' }, { key: 'from', label: 'Origin', left: true }, { key: 'at', label: 'Set', left: true }], keys.map((k) => ({ k: h('code', k), v: bus[k].value, from: bus[k].from === 'user' ? 'entered by you' : bus[k].from === 'site' ? 'live site data' : bus[k].from, at: new Date(bus[k].at).toLocaleDateString() })), { max: 300 }) : h('div.empty', 'No project values yet.'));
  };

  // which modules take which kind of input, per suite
  const idx = await suiteIndex();
  const bySuite = h('div');
  for (const s of idx) {
    const d = h('details.fold', h('summary', `${s.n}. ${SHORT[s.n]}`), h('div', h('small.muted', 'Open to load…')));
    d.addEventListener('toggle', async () => {
      if (!d.open || d.dataset.done) return;
      d.dataset.done = '1';
      const [spec, mods] = await Promise.all([suiteSpec(s.n), suiteModules(s.n)]);
      const tableInputs = mods.flatMap((m) => (m.inputs || []).filter((i) => ['table', 'series', 'geometry'].includes(i.type)).map((i) => ({ what: i.label, type: i.type === 'geometry' ? `imported geometry (${(i.kinds || []).join(' / ')})` : i.type === 'table' ? `table: ${i.columns.map((c) => c.label).join(', ')}` : `two-column series: ${i.xLabel || 'x'}, ${i.yLabel || 'y'}`, mod: h('a', { href: `#/m/${m.id}` }, m.title) })));
      clear(d.lastChild).append(h('h4', 'Data this suite is specified to accept'), h('div.prose', (spec.inputs || []).map((x) => h('p', h('b.g', `${x.cat}. `), x.text))),
        h('h4', { style: { marginTop: '12px' } }, 'Tables, series and geometry its modules take'), tableInputs.length ? table([{ key: 'what', label: 'Input', left: true }, { key: 'type', label: 'Form', left: true }, { key: 'mod', label: 'Module', left: true }], tableInputs, { max: 200 }) : h('p.muted', 'Modules of this suite take scalar inputs only, or are not installed.'));
    });
    bySuite.append(d);
  }

  root.append(h('h1', 'Input data portal'), h('p.lede', 'Feed the case in here. Drop any mix of files — geometry, routes, bathymetry, saved projects, run records, tables of design values — and each is recognised and sent where it belongs.'),
    dropZone(route, '', 'Drop any input files here'), log,
    h('div.grid.c2', { style: { marginTop: '16px' } },
      h('div.card', h('h3', 'What is recognised'), h('ul.ticks',
        h('li.on', h('b', 'Geometry and meshes'), ' — CAD, structural models, finite-element and CFD meshes, panel meshes, cross-sections, airfoils. ', h('a', { href: '#/geometry' }, 'Full format list')),
        h('li.on', h('b', 'Routes and layouts'), ' — pipeline, cable, riser and mooring centrelines; turbine and anchor coordinates (CSV, GeoJSON, Shapefile, DXF, KML…).'),
        h('li.on', h('b', 'Seabed and terrain'), ' — GeoTIFF, ASCII grids, NetCDF, LandXML, XYZ soundings, point clouds.'),
        h('li.on', h('b', 'Project values'), ' — a CSV or text file of “key, value” lines, or a JSON object, using the names listed below (SI units).'),
        h('li.on', h('b', 'Run records'), ' — the JSON a module exports: restores that module’s inputs exactly.'),
        h('li.on', h('b', 'Whole projects'), ' — a project export from this app on any device.'),
        h('li.on', h('b', 'Module inputs in bulk'), ' — a JSON object { "module-id": { "input": value } }.'))),
      h('div.card', h('h3', 'Type or paste project values'), h('p.muted', 'One “key = value” per line, SI units. These become the default in every analysis field linked to them.'), paste,
        h('div.row', { style: { marginTop: '8px' } }, h('button.btn.primary', { type: 'button', onclick: applyPaste }, 'Set values'), h('a.btn', { href: '#/site' }, icon('globe'), 'Or pull them live for a site')))),
    h('h2', { style: { marginTop: '22px' } }, 'Templates'), h('p.muted', 'Download, fill in a spreadsheet, drop back here or paste into the matching table of a module.'),
    h('div.row', Object.keys(TEMPLATES).map((name) => h('button.btn', { type: 'button', onclick: () => download(name, TEMPLATES[name], 'text/csv') }, icon('download'), name))),
    h('h2', { style: { marginTop: '22px' } }, 'Project values now'), busBox,
    h('div.row', { style: { marginTop: '8px' } }, h('button.btn', { type: 'button', onclick: () => download('project-values.json', JSON.stringify(busValues(), null, 1), 'application/json') }, icon('download'), 'Export values'),
      h('button.btn', { type: 'button', onclick: () => download('module-inputs.json', JSON.stringify(getProject().inputs, null, 1), 'application/json') }, icon('download'), 'Export all module inputs'),
      h('button.btn', { type: 'button', onclick: () => { if (confirm('Clear every project value? Module inputs you typed are kept.')) { getProject().bus = {}; save(); drawBus(); } } }, icon('trash'), 'Clear values')),
    h('h2', { style: { marginTop: '22px' } }, 'What each suite accepts'), bySuite,
    pager({ href: '#/site', label: 'Project & site data' }, { href: '#/geometry', label: 'Geometry & mesh' }));
  drawBus();
  void fmt;
  return { title: 'Input data portal', crumbs: [{ href: '#/portal', label: 'Input data portal' }] };
}
