// Project & site: choose any location on Earth, pull live metocean / seabed / seismic / resource data,
// review the derived design values and override any of them with measured data.
import { h, clear, fmt, icon, toast, table, download } from '../dom.js';
import { renderPlot } from '../plot.js';
import { getProject, save, listProjects, newProject, openProject, deleteProject, importProject, emit, publish } from '../../core/store.js';
import { pullSite, SITE_KEYS } from '../../data/site.js';
import { searchPlace, bathyGrid } from '../../data/connectors.js';
import { fieldChecks } from './fieldcheck.js';
import { pager } from '../shell.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

/** Equirectangular world map drawn from bundled coastline data: works offline, click or tap to pick a point. */
async function worldMap(onPick, initial) {
  const W = 720, H = 360;
  let view = { cx: initial?.lon ?? 0, cy: initial?.lat ?? 10, z: initial ? 3 : 1 };
  const svg = h('svg', { viewBox: `0 0 ${W} ${H}`, role: 'application', 'aria-label': 'World map: click or tap to choose the site', style: { width: '100%', borderRadius: '10px', background: 'var(--accent-soft)', cursor: 'crosshair', touchAction: 'none' } });
  const land = h('path', { fill: 'var(--surface-3)', stroke: 'var(--ink-3)', 'stroke-width': 0.4, 'vector-effect': 'non-scaling-stroke' });
  const grat = h('path', { fill: 'none', stroke: 'var(--line)', 'stroke-width': 0.5, 'vector-effect': 'non-scaling-stroke' });
  const pin = h('g', { visibility: 'hidden' }, h('circle', { r: 7, fill: 'var(--critical)', stroke: '#fff', 'stroke-width': 2 }), h('circle', { r: 14, fill: 'none', stroke: 'var(--critical)', 'stroke-width': 1.5, opacity: 0.6 }));
  const g = h('g', grat, land);
  svg.append(g, pin);
  const world = await fetch(new URL('../../../assets/world.json', import.meta.url)).then((r) => r.json()).catch(() => ({ rings: [] }));
  let d = '';
  for (const ring of world.rings) { d += `M${ring[0]},${-ring[1]}`; for (let i = 2; i < ring.length; i += 2) d += `L${ring[i]},${-ring[i + 1]}`; d += 'Z'; }
  land.setAttribute('d', d);
  let gd = '';
  for (let lon = -180; lon <= 180; lon += 30) gd += `M${lon},-90L${lon},90`;
  for (let lat = -60; lat <= 60; lat += 30) gd += `M-180,${lat}L180,${lat}`;
  grat.setAttribute('d', gd);
  let picked = initial ? { lat: initial.lat, lon: initial.lon } : null;
  const apply = () => {
    const s = (W / 360) * view.z;
    view.cx = Math.max(-180 + 180 / view.z, Math.min(180 - 180 / view.z, view.cx));
    view.cy = Math.max(-90 + 90 / view.z, Math.min(90 - 90 / view.z, view.cy));
    g.setAttribute('transform', `translate(${W / 2 - view.cx * s} ${H / 2 + view.cy * s}) scale(${s})`);
    if (picked) { pin.setAttribute('visibility', 'visible'); pin.setAttribute('transform', `translate(${W / 2 + (picked.lon - view.cx) * s} ${H / 2 - (picked.lat - view.cy) * s})`); }
  };
  const toGeo = (e) => { const r = svg.getBoundingClientRect(); const s = (W / 360) * view.z; return { lon: view.cx + ((e.clientX - r.left) / r.width * W - W / 2) / s, lat: view.cy - ((e.clientY - r.top) / r.height * H - H / 2) / s }; };
  let drag = null;
  svg.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY, cx: view.cx, cy: view.cy, moved: false }; svg.setPointerCapture(e.pointerId); });
  svg.addEventListener('pointermove', (e) => {
    if (!drag) return;
    const r = svg.getBoundingClientRect(); const s = (r.width / 360) * view.z;
    if (Math.hypot(e.clientX - drag.x, e.clientY - drag.y) > 5) drag.moved = true;
    if (drag.moved) { view.cx = drag.cx - (e.clientX - drag.x) / s; view.cy = drag.cy + (e.clientY - drag.y) / s; apply(); }
  });
  svg.addEventListener('pointerup', (e) => { if (drag && !drag.moved) { const p = toGeo(e); set(p.lat, p.lon, true); } drag = null; });
  svg.addEventListener('wheel', (e) => { e.preventDefault(); view.z = Math.max(1, Math.min(60, view.z * (e.deltaY < 0 ? 1.3 : 1 / 1.3))); apply(); }, { passive: false });
  function set(lat, lon, notify) { picked = { lat: Math.max(-89.9, Math.min(89.9, lat)), lon: ((lon + 540) % 360) - 180 }; apply(); if (notify) onPick(picked.lat, picked.lon); }
  apply();
  const zoom = (f) => { view.z = Math.max(1, Math.min(60, view.z * f)); apply(); };
  return { el: h('div', svg, h('div.row', { style: { marginTop: '6px' } }, h('button.btn.small', { type: 'button', onclick: () => zoom(1.6), 'aria-label': 'Zoom in' }, '+ Zoom in'), h('button.btn.small', { type: 'button', onclick: () => zoom(1 / 1.6), 'aria-label': 'Zoom out' }, '− Zoom out'), h('small.muted', 'Click or tap to place the site · drag to move · wheel to zoom'))), set: (lat, lon) => { view.cx = lon; view.cy = lat; view.z = Math.max(view.z, 4); set(lat, lon, false); } };
}

export default async function site(root) {
  const p = getProject();
  const lat = h('input.plain', { type: 'number', step: 'any', min: -90, max: 90, value: p.site.lat ?? '', placeholder: 'Latitude', 'aria-label': 'Latitude in decimal degrees', style: { width: '130px' } });
  const lon = h('input.plain', { type: 'number', step: 'any', min: -180, max: 180, value: p.site.lon ?? '', placeholder: 'Longitude', 'aria-label': 'Longitude in decimal degrees', style: { width: '130px' } });
  const place = h('input.plain', { type: 'search', placeholder: 'Search a place, port or field name…', 'aria-label': 'Search a place', style: { maxWidth: '320px' } });
  const hits = h('div.row');
  const years = h('select.plain', { style: { width: '120px' }, 'aria-label': 'Length of hindcast' }, [10, 20, 30, 40].map((y) => h('option', { value: y, selected: y === 20 }, `${y} years`)));
  const prog = h('div');
  const summary = h('div');
  const map = await worldMap((la, lo) => { lat.value = la.toFixed(4); lon.value = lo.toFixed(4); }, p.site.lat != null ? p.site : null);

  const find = async () => {
    if (!place.value.trim()) return;
    clear(hits).append(h('small.muted', 'Searching…'));
    try {
      const r = await searchPlace(place.value.trim());
      clear(hits);
      const list = r.data.results || [];
      if (!list.length) hits.append(h('small.muted', 'No place found — enter the coordinates or click the map.'));
      for (const x of list) hits.append(h('button.btn.small', { type: 'button', onclick: () => { lat.value = x.latitude.toFixed(4); lon.value = x.longitude.toFixed(4); map.set(x.latitude, x.longitude); getProject().site.name = `${x.name}${x.country ? `, ${x.country}` : ''}`; } }, `${x.name}${x.admin1 ? `, ${x.admin1}` : ''}${x.country_code ? ` (${x.country_code})` : ''}`));
    } catch (e) { clear(hits).append(h('small.muted', `Search unavailable: ${e.message}`)); }
  };
  place.addEventListener('keydown', (e) => { if (e.key === 'Enter') find(); });

  const pull = async (force) => {
    const la = Number(lat.value), lo = Number(lon.value);
    if (!Number.isFinite(la) || !Number.isFinite(lo) || Math.abs(la) > 90 || Math.abs(lo) > 180 || lat.value === '' || lon.value === '') { toast('Enter a latitude between −90 and 90 and a longitude between −180 and 180.'); return; }
    map.set(la, lo);
    const bar = h('i'), msg = h('small.muted', 'Starting…');
    clear(prog).append(h('div.progress', bar), msg);
    try {
      const r = await pullSite(la, lo, { years: Number(years.value), force, onProgress: (f, m) => { bar.style.width = `${Math.round(f * 100)}%`; msg.textContent = m; } });
      clear(prog);
      if (r.data.issues.length) prog.append(h('div.banner.warn', icon('warn'), h('div', h('b', 'Some sources did not answer: '), r.data.issues.join(' · '))));
      toast('Site data pulled and shared with every suite.');
      draw();
    } catch (e) { clear(prog).append(h('div.banner.fail', icon('warn'), e.message)); }
  };

  function draw() {
    clear(summary);
    const s = getProject().site, d = s.data || {}, bus = getProject().bus;
    if (!s.pulled && !Object.keys(bus).some((k) => k.startsWith('site.'))) { summary.append(h('div.empty', 'No site data yet. Choose a location and press “Pull live data”.')); return; }
    if (s.pulled) summary.append(h('p.muted', `${s.name} — ${Number(s.lat).toFixed(4)}°, ${Number(s.lon).toFixed(4)}° · pulled ${new Date(s.pulled).toLocaleString()}`));

    // design values table: editable, with provenance
    const rows = SITE_KEYS.map(([key, label, unit]) => {
      const b = bus[key];
      const inp = h('input', { type: 'text', inputmode: 'decimal', value: b ? Number(b.value.toPrecision(5)) : '', 'aria-label': label, style: { width: '110px', textAlign: 'right', border: '1px solid var(--line)', borderRadius: '6px', padding: '3px 6px', background: 'var(--surface)', color: 'var(--ink)' } });
      inp.addEventListener('change', () => { const v = Number(inp.value); if (Number.isFinite(v) && inp.value !== '') { publish({ [key]: v }, 'user'); toast(`${label} set by you — it now overrides the pulled estimate.`); draw(); } });
      return { label, v: inp, unit, from: b ? (b.from === 'user' ? 'entered by you' : b.from === 'site' ? 'derived from live data' : b.from) : '—' };
    });
    summary.append(h('div.card', h('h3', 'Design values shared with every suite'), h('p.muted', 'Derived from the live sources below. Type over any value to replace it with measured or project-specific data — yours always wins.'),
      table([{ key: 'label', label: 'Quantity', left: true }, { key: 'v', label: 'Value' }, { key: 'unit', label: 'Unit', left: true }, { key: 'from', label: 'Origin', left: true }], rows)));

    const charts = h('div.grid.c2', { style: { marginTop: '14px' } });
    if (d.waves?.extremes?.maxima) {
      const ex = d.waves.extremes, T = [1.5, 2, 5, 10, 20, 50, 100, 200];
      charts.append(renderPlot({ type: 'line', title: 'Significant wave height — return values', subtitle: ex.method, x: { label: 'Return period', unit: 'years', log: true }, y: { label: 'Hs', unit: 'm', zero: false },
        series: [{ name: 'Gumbel fit', x: T, y: T.map((t) => ex.fit.loc - ex.fit.scale * Math.log(-Math.log(1 - 1 / t))) },
          { name: 'Annual maxima', x: ex.maxima.slice().sort((a, b) => a - b).map((_, i, a) => 1 / (1 - (i + 1 - 0.44) / (a.length + 0.12))), y: ex.maxima.slice().sort((a, b) => a - b), markers: true, line: false }].map((q, i) => (i === 1 ? { ...q, dash: false } : q)), }));
      charts.append(renderPlot({ type: 'bar', title: 'Wave climate by month', categories: MONTHS, series: [{ name: 'Mean Hs', y: d.waves.monthlyMean }, { name: '90th percentile Hs', y: d.waves.monthlyP90 }], y: { label: 'Hs', unit: 'm' } }));
      if (d.waves.scatter?.prob) {
        const sc = d.waves.scatter;
        const mid = (e) => Array.from(e).slice(0, -1).map((v, i) => (v + e[i + 1]) / 2);
        charts.append(renderPlot({ type: 'heatmap', title: 'Wave scatter diagram (Hs – Tp)', xs: mid(sc.tpEdges), ys: mid(sc.hsEdges), z: sc.prob.map((r) => Array.from(r).map((v) => v * 100)), x: { label: 'Tp', unit: 's' }, y: { label: 'Hs', unit: 'm' }, zLabel: 'Occurrence', zUnit: '%' }));
      }
      charts.append(renderPlot({ type: 'line', title: 'Workability — time Hs is exceeded', x: { label: 'Hs threshold', unit: 'm' }, y: { label: 'Time exceeded', unit: '%' }, series: [{ name: 'Exceedance', x: d.waves.exceed.map((e) => e.hs), y: d.waves.exceed.map((e) => e.frac * 100), markers: true }] }));
    }
    if (d.wind?.rose) {
      charts.append(renderPlot({ type: 'polar', title: 'Wind rose at 100 m (direction the wind comes from)', rUnit: '%', series: [{ name: 'Frequency', thetaDeg: [...d.wind.rose.map((r) => r.dir), 360], r: [...d.wind.rose.map((r) => r.freq * 100), d.wind.rose[0].freq * 100], fill: true, closed: true }] }));
      charts.append(renderPlot({ type: 'bar', title: 'Mean wind speed at 100 m by month', categories: MONTHS, series: [{ name: 'Mean wind', y: d.wind.monthlyMean }], y: { label: 'Wind speed', unit: 'm/s' } }));
    }
    if (d.forecast?.time) {
      const t = d.forecast.time.map((_, i) => i / 24);
      charts.append(renderPlot({ type: 'line', title: 'Sea-state forecast — next 8 days', x: { label: 'Days from now' }, y: { label: 'Wave height', unit: 'm' }, series: [{ name: 'Total Hs', x: t, y: d.forecast.wave_height }, { name: 'Swell', x: t, y: d.forecast.swell_wave_height || [] }, { name: 'Wind sea', x: t, y: d.forecast.wind_wave_height || [] }] }));
    }
    summary.append(charts);
    if (d.seismic) summary.append(h('div.card', { style: { marginTop: '14px' } }, h('h3', 'Seismicity'), h('p', `${d.seismic.n} earthquakes of magnitude 4.5 or more within 300 km since 1970${d.seismic.maxMag ? `; largest M ${fmt(d.seismic.maxMag, 2)}` : ''}. Screening PGA (475-year) ≈ ${fmt(d.seismic.pga475, 2)} g.`), h('p.muted', d.seismic.note)));
    if (d.seabed) summary.append(h('div.card', { style: { marginTop: '14px' } }, h('h3', 'Seabed sediment'), d.seabed.cls
      ? h('p', h('b', `${d.seabed.label}. `), d.seabed.basis === 'map' ? `Mapped at the site as “${d.seabed.described}” on the EMODnet Geology substrate map at 1:${Number(d.seabed.scale).toLocaleString('en')}${d.seabed.holder ? ` (${d.seabed.holder})` : ''}. `
        : `${d.seabed.agree} of the ${d.seabed.n} nearest described bottom samples (${fmt(d.seabed.nearestKm, 2)}–${fmt(d.seabed.farthestKm, 2)} km away, surveyed ${d.seabed.years[0]}–${d.seabed.years[1]}) are of this class; nearest descriptions: ${d.seabed.described}. `,
        d.seabed.d50 ? `Representative median grain size ${fmt(d.seabed.d50 * 1000, 2)} mm, used as the default for scour and sediment-transport analyses. ` : 'No grain size is set. ', d.seabed.means)
      : h('p', 'No open seabed-sediment map or described sample covers this location (open coverage is European seas and United States waters). Enter the soil values from a site investigation or a regional geological map in the table above.'),
      h('p.why', 'This describes the top few centimetres of the seabed. It indicates what kind of ground to expect and sets a grain size; it does not give strength, layering or depth to rock, and no strength value is derived from it. Shear strength, friction angle and unit weight must come from a geotechnical investigation (boreholes, cone penetration tests) before any foundation is sized.')));
    if (d.kept?.length) summary.append(h('div.banner', icon('info'), `Kept as you entered them, not replaced by this pull: ${d.kept.map((k) => (SITE_KEYS.find((r) => r[0] === k) || [k, k])[1]).join(', ')}.`));
    if (d.provenance) summary.append(h('details.fold', h('summary', 'Where each number comes from'), h('div', table([{ key: 'what', label: 'Quantity', left: true }, { key: 'source', label: 'Source', left: true }, { key: 'method', label: 'Reduction', left: true }, { key: 'when', label: 'Fetched', left: true }],
      Object.entries(d.provenance).map(([k, v]) => ({ what: k, source: v.source, method: v.method || '—', when: v.at ? new Date(v.at).toLocaleString() : '—' }))))));

    const bat = h('div.card', { style: { marginTop: '14px' } }, h('h3', 'Seabed around the site'),
      h('div.row', h('button.btn', { type: 'button', onclick: async (e) => {
        e.target.disabled = true;
        try {
          const g = (await bathyGrid(Number(s.lat), Number(s.lon), { halfKm: 10, n: 25 })).data;
          const xs = Array.from({ length: g.nx }, (_, i) => (g.x0 + i * g.dx) / 1000), ys = Array.from({ length: g.ny }, (_, j) => (g.y0 + j * g.dy) / 1000);
          const z = ys.map((_, j) => xs.map((__, i) => g.z[j * g.nx + i]));
          bat.append(renderPlot({ type: 'contour', title: 'Bathymetry ±10 km (elevation, negative below sea level)', xs, ys, z, x: { label: 'East', unit: 'km' }, y: { label: 'North', unit: 'km' }, zLabel: 'Elevation', zUnit: 'm', diverging: false }));
          const { addGeometry } = await import('../../core/store.js');
          await addGeometry({ name: `Bathymetry at ${s.name}`, kind: 'grid', kinds: ['grid'], format: 'NOAA NCEI DEM', stats: `${g.nx}×${g.ny} cells` }, [{ ...g, z: Float32Array.from(g.z) }]);
          toast('Bathymetry grid added to the project geometry.');
        } catch (err) { toast(`Bathymetry not available: ${err.message}`); e.target.disabled = false; }
      } }, icon('download'), 'Fetch bathymetry grid')), h('p.muted', 'Samples the NOAA NCEI global mosaic on a 25 × 25 grid and stores it as project geometry for route, foundation and installation analyses. A site survey supersedes it.'));
    if (s.lat != null) summary.append(bat, fieldChecks({ lat: Number(s.lat), lon: Number(s.lon) }));
  }

  // ----- project management -----
  const name = h('input.plain', { type: 'text', value: p.name, 'aria-label': 'Project name', style: { maxWidth: '320px' } });
  name.addEventListener('change', () => { getProject().name = name.value || 'Untitled project'; save(); emit('project', getProject()); });
  const projList = h('div');
  const drawProjects = async () => {
    const list = await listProjects();
    clear(projList).append(table([{ key: 'name', label: 'Project', left: true }, { key: 'created', label: 'Created', left: true }, { key: 'runs', label: 'Analyses run' }, { key: 'act', label: '', left: true }],
      list.map((x) => ({ ...x, created: new Date(x.created).toLocaleDateString(), act: h('span.row', x.id === getProject().id ? h('span.chip.ok', 'open') : h('button.btn.small', { type: 'button', onclick: async () => { await openProject(x.id); location.reload(); } }, 'Open'),
        h('button.btn.small', { type: 'button', 'aria-label': `Delete ${x.name}`, onclick: async () => { if (confirm(`Delete project “${x.name}” and its results from this device? This cannot be undone.`)) { await deleteProject(x.id); drawProjects(); } } }, icon('trash'))) }))));
  };
  const file = h('input', { type: 'file', accept: '.json', style: { display: 'none' }, onchange: async () => { try { importProject(JSON.parse(await file.files[0].text())); location.reload(); } catch { toast('That file is not a project export.'); } } });

  root.append(h('h1', 'Project & site data'), h('p.lede', 'Everything the suites need to know about where the asset stands. Pulled live from public sources by your own browser, stored on your device, and shared with all 17 suites.'),
    h('div.card', h('h3', '1 · Choose the location'),
      h('div.row', place, h('button.btn', { type: 'button', onclick: find }, icon('search'), 'Search'), hits),
      h('div', { style: { margin: '10px 0' } }, map.el),
      h('div.row', h('label', 'Latitude'), lat, h('label', 'Longitude'), lon, h('label', 'Hindcast'), years,
        h('button.btn.primary', { type: 'button', onclick: () => pull(false) }, icon('globe'), 'Pull live data'), h('button.btn', { type: 'button', title: 'Ignore the cache and fetch everything again', onclick: () => pull(true) }, icon('refresh'), 'Refresh')),
      prog),
    h('h2', { style: { marginTop: '20px' } }, '2 · Site conditions'), summary,
    h('h2', { style: { marginTop: '24px' } }, 'Projects on this device'),
    h('div.card', h('div.row', h('label', 'Current project name'), name,
      h('button.btn', { type: 'button', onclick: async () => { const n = prompt('Name of the new project'); if (n) { await newProject(n); location.reload(); } } }, '+ New project'),
      h('button.btn', { type: 'button', onclick: () => download(`${getProject().name.replace(/[^\w-]+/g, '_')}.project.json`, JSON.stringify(getProject()), 'application/json') }, icon('download'), 'Export'),
      file, h('button.btn', { type: 'button', onclick: () => file.click() }, icon('upload'), 'Import')), h('div', { style: { marginTop: '10px' } }, projList)),
    pager({ href: '#/', label: 'Home' }, { href: '#/portal', label: 'Input data portal' }));
  draw(); drawProjects();
  return { title: 'Project & site data', crumbs: [{ href: '#/site', label: 'Project & site data' }] };
}
