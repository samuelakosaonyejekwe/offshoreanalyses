// Geometry & mesh workbench: import any supported file, see it, check it, repair it, refine or coarsen it,
// and prove that a measured quantity no longer depends on the mesh. This is not a CAD tool — geometry comes
// from the originating design system and is prepared here for analysis.
import { h, clear, fmt, icon, toast, table, download, tabs } from '../dom.js';
import { renderPlot } from '../plot.js';
import { getProject, addGeometry, loadGeometry, removeGeometry } from '../../core/store.js';
import { pager } from '../shell.js';

const KIND_LABEL = { surface: 'Surface mesh', volume: 'Volume mesh', frame: 'Frame (nodes and members)', polyline: 'Route / centreline', grid: 'Bathymetry / terrain grid', points: 'Points', section: 'Cross-section' };

function statsOf(m) {
  if (m.kind === 'surface') return `${fmt(m.indices.length / 3)} triangles`;
  if (m.kind === 'volume') return `${fmt((m.cells || []).reduce((s, c) => s + c.conn.length / ({ tet4: 4, tet10: 10, hex8: 8, hex20: 20, wedge6: 6, pyr5: 5 }[c.type] || 4), 0))} cells`;
  if (m.kind === 'frame') return `${m.members.length} members, ${m.nodes.length / 3} nodes`;
  if (m.kind === 'polyline') return `${m.points.length / 3} points`;
  if (m.kind === 'grid') return `${m.nx}×${m.ny} cells`;
  if (m.kind === 'points') return `${fmt(m.positions.length / 3)} points`;
  if (m.kind === 'section') return `${m.loops.length} loop${m.loops.length > 1 ? 's' : ''}`;
  return '';
}

/** Import files through the format registry and store every resulting model in the project. */
export async function importFiles(fileList, onStatus) {
  const { importFile } = await import('../../geometry/registry.js');
  const files = [...fileList];
  const results = [];
  // multi-file formats (shapefile sets, glTF with buffers, polyMesh folders) need their companions: pass the whole selection
  const all = new Map();
  for (const f of files) all.set(f.name, new Uint8Array(await f.arrayBuffer()));
  const companions = /\.(dbf|prj|shx|cpg|bin|mtl)$/i;
  // DWG is read by a GPL-licensed library that is not part of the app: it is downloaded only if the user agrees, once per batch
  let dwg = null;
  if (files.some((f) => /\.dwg$/i.test(f.name))) {
    const a = await import('../../geometry/adapters/dwg.js');
    dwg = { on: a.ENABLED, consent: a.ENABLED && window.confirm(`Reading a DWG drawing needs the LibreDWG library, which is not part of this app.\n\nIf you continue, your browser downloads it (${a.library.licence}, about ${Math.round(a.sizeBytes / 1e6)} MB${a.fetchedBefore() ? ', probably already in this browser\'s cache' : ''}) from cdn.jsdelivr.net and runs it in a separate worker on this device. Your drawing is not uploaded anywhere.\n\nCancel to skip the DWG files — you can save them as DXF and import those instead.`) };
  }
  for (const f of files) {
    if (companions.test(f.name) && files.length > 1) continue;
    onStatus?.(`Reading ${f.name}…`);
    try {
      const bytes = all.get(f.name);
      const isDwg = /\.dwg$/i.test(f.name);
      if (isDwg && dwg.on && !dwg.consent) throw new Error('Skipped: the DWG library was not downloaded. Save the drawing as DXF in the originating program and import that, or import again and agree to the download.');
      const res = await importFile({ name: f.name, bytes, files: all }, isDwg ? { consent: dwg.consent } : {});
      if (!res.models?.length) throw new Error('The file was read but contains no geometry this application can use.');
      const kinds = [...new Set(res.models.map((m) => m.kind))];
      const id = await addGeometry({ name: f.name, kind: kinds[0], kinds, format: res.format?.name || res.format || '', stats: res.models.map(statsOf).filter(Boolean).slice(0, 2).join('; '), warnings: (res.warnings || []).slice(0, 12), imported: new Date().toISOString() }, res.models);
      results.push({ file: f.name, ok: true, id, models: res.models.length, warnings: res.warnings || [] });
    } catch (err) {
      results.push({ file: f.name, ok: false, error: err.message });
    }
  }
  return results;
}

export function dropZone(onFiles, accept, text) {
  const input = h('input', { type: 'file', multiple: true, accept, style: { display: 'none' }, onchange: () => { const picked = [...input.files]; input.value = ''; if (picked.length) onFiles(picked); } });   // copy first: clearing the input empties the live FileList
  const zone = h('div.drop', { tabindex: 0, role: 'button', 'aria-label': 'Choose files to import', onclick: () => input.click(), onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); input.click(); } } },
    icon('upload'), h('p', h('b', text || 'Drop files here, or click to choose')), h('small.muted', 'Several files at once are fine. Files are read on this device and are not uploaded anywhere.'), input);
  zone.querySelector('svg').style.cssText = 'width:32px;height:32px';
  zone.addEventListener('dragover', (e) => { e.preventDefault(); zone.classList.add('over'); });
  zone.addEventListener('dragleave', () => zone.classList.remove('over'));
  zone.addEventListener('drop', (e) => { e.preventDefault(); zone.classList.remove('over'); if (e.dataTransfer.files.length) onFiles([...e.dataTransfer.files]); });
  return zone;
}

async function formatMatrix() {
  const { FORMATS } = await import('../../geometry/registry.js');
  const cats = [...new Set(FORMATS.map((f) => f.category))];
  return h('div', h('p.muted', 'Read natively means the reader is built in and works offline. Adapter means a larger reader that is fetched the first time it is needed. Convert first means the format needs a licensed geometry kernel: export the stated neutral format from the originating system.'),
    cats.map((c) => h('details.fold', { open: c === 'cad' }, h('summary', `${c[0].toUpperCase()}${c.slice(1)} — ${FORMATS.filter((f) => f.category === c).length} formats`), h('div',
      table([{ key: 'name', label: 'Format', left: true }, { key: 'ext', label: 'Extensions', left: true }, { key: 'status', label: 'Support', left: true }, { key: 'kinds', label: 'Becomes', left: true }, { key: 'suites', label: 'Suites', left: true }, { key: 'notes', label: 'Notes', left: true }],
        FORMATS.filter((f) => f.category === c).map((f) => ({ name: f.name, ext: (f.extensions || []).join(' '), status: h('span', { class: `chip ${f.status === 'native' ? 'ok' : f.status === 'adapter' ? 'info' : 'warn'}` }, f.status === 'native' ? '✓ read natively' : f.status === 'adapter' ? 'adapter' : 'convert first'), kinds: (f.kinds || []).join(', '), suites: (f.suites || []).length >= 16 ? 'all' : (f.suites || []).join(', '), notes: f.notes || '' })))))));
}

async function workbench(entry, box) {
  clear(box).append(h('small.muted', 'Loading model…'));
  const models = await loadGeometry(entry.id);
  if (!models?.length) { clear(box).append(h('div.banner.fail', 'This model is no longer stored on the device.')); return; }
  const [topology, measure, quality, repairM, refineM, conv, viewerM] = await Promise.all([
    import('../../mesh/topology.js'), import('../../mesh/measure.js'), import('../../mesh/quality.js'), import('../../mesh/repair.js'), import('../../mesh/refine.js'), import('../../mesh/convergence.js'), import('../viewer3d.js')]);
  let idx = 0;
  let model = models[0];
  const host = h('div.viewer');
  const info = h('div');
  const pickInfo = h('small.muted', 'Click the model to read coordinates.');
  clear(box);
  if (models.length > 1) box.append(h('div.row', h('label', 'Part'), h('select.plain', { style: { maxWidth: '360px' }, onchange: (e) => { idx = Number(e.target.value); model = models[idx]; show(); } }, models.map((m, i) => h('option', { value: i }, `${m.name || `part ${i + 1}`} — ${KIND_LABEL[m.kind] || m.kind}`)))));
  box.append(host, h('div.row', { style: { margin: '6px 0' } },
    h('button.btn.small', { type: 'button', onclick: () => download(`${entry.name}.png`, dataUrlBlob(viewer.snapshot())) }, icon('download'), 'Image'), pickInfo), info);
  const viewer = viewerM.createViewer(host, { edges: true });
  viewer.onPick((p) => { pickInfo.textContent = p ? `x ${fmt(p.point[0])}  y ${fmt(p.point[1])}  z ${fmt(p.point[2])} m${p.value != null ? ` · value ${fmt(p.value)}` : ''}` : 'Click the model to read coordinates.'; });
  let modelId;
  const save = async (label) => { models[idx] = model; const id = await addGeometry({ ...entry, id: undefined, name: `${entry.name} (${label})`, kind: model.kind, kinds: [...new Set(models.map((m) => m.kind))], stats: statsOf(model) }, models); toast(`Saved as a new project model: ${entry.name} (${label}).`); return id; };

  function show() {
    viewer.clear(); modelId = viewer.addModel(model); viewer.fit();
    clear(info);
    const rep = topology.inspect(model);
    const mes = safe(() => measure.measure(model)) || {};
    const kp = [];
    if (mes.bbox || rep.stats?.bbox) { const b = mes.bbox || rep.stats.bbox; const d = b.size || [b.max[0] - b.min[0], b.max[1] - b.min[1], b.max[2] - b.min[2]]; kp.push(['Overall size', `${fmt(d[0], 4)} × ${fmt(d[1], 4)} × ${fmt(d[2], 4)}`, 'm']); }
    if (mes.area != null) kp.push(['Surface area', fmt(mes.area), 'm²']);
    if (mes.volume != null && rep.watertight) kp.push(['Enclosed volume', fmt(Math.abs(mes.volume)), 'm³']);
    if (mes.length != null) kp.push(['Length', fmt(mes.length), 'm']);
    kp.push(['Content', statsOf(model), '']);
    info.append(h('div.kpis', kp.map(([l, v, u]) => h('div.kpi', h('span', l), h('b', { style: { fontSize: '1.05rem' } }, v), u ? h('em', u) : null))));

    const checkTab = () => {
      const el = h('div');
      el.append(h('div.banner', { class: rep.ok ? '' : 'warn' }, icon(rep.ok ? 'check' : 'warn'), h('div', h('b', rep.ok ? 'No problems found. ' : `${rep.issues.length} kind${rep.issues.length > 1 ? 's' : ''} of problem found. `),
        model.kind === 'surface' ? (rep.watertight ? 'The surface is closed (watertight) and consistently oriented, so it can be used for volume, hydrostatic and solid analyses.' : 'The surface is not closed. Area-based uses work as it is; volume, hydrostatics and solid meshing need a watertight surface — try Repair.') : '')));
      if (rep.issues.length) el.append(table([{ key: 'message', label: 'Finding', left: true }, { key: 'count', label: 'Count' }, { key: 'show', label: '', left: true }], rep.issues.map((i) => ({ ...i, show: (i.edges || i.elements || i.vertices) ? h('button.btn.small', { type: 'button', onclick: () => viewer.highlight(modelId, i) }, 'Show') : '' }))));
      el.append(h('div.row', { style: { marginTop: '10px' } }, h('button.btn.primary', { type: 'button', onclick: () => {
        const r = repairM.repair(model, {});
        model = r.model; show();
        toast(r.actions.length ? `Repair applied: ${r.actions.map((a) => a.message).join('; ')}`.slice(0, 220) : 'Nothing needed repairing.');
      } }, 'Repair automatically'), h('button.btn', { type: 'button', onclick: () => save('repaired') }, 'Save as new model')));
      if (entry.warnings?.length) el.append(h('details.fold', h('summary', `Notes from the import (${entry.warnings.length})`), h('div', h('ul.prose', entry.warnings.map((w) => h('li', w))))));
      return el;
    };

    const qualityTab = () => {
      if (!['surface', 'volume'].includes(model.kind)) return h('div.empty', 'Element quality applies to surface and volume meshes.');
      const q = quality.assess(model, {});
      const el = h('div');
      el.append(h('div.banner', { class: q.grade === 'A' || q.grade === 'B' ? '' : 'warn' }, icon('info'), h('div', h('b', `Mesh quality grade ${q.grade}. `), (q.diagnosis || []).join(' '))));
      const hist = q.histograms?.aspectRatio || q.histogram;
      if (q.perElement?.aspectRatio) el.append(h('div.grid.c2',
        renderPlot({ type: 'hist', title: 'Element aspect ratio', values: Array.from(q.perElement.aspectRatio).filter(Number.isFinite).slice(0, 200000), x: { label: 'Aspect ratio (1 is ideal)' }, bins: 24 }),
        q.perElement.minAngle ? renderPlot({ type: 'hist', title: 'Smallest corner angle', values: Array.from(q.perElement.minAngle).filter(Number.isFinite).slice(0, 200000), x: { label: 'Angle', unit: '°' }, bins: 24 }) : h('div')));
      void hist;
      if (q.summary) el.append(h('details.fold', h('summary', 'Quality statistics'), h('div', h('pre', { style: { fontSize: '.78rem', overflowX: 'auto' } }, JSON.stringify(q.summary, (k, v) => (typeof v === 'number' ? Number(v.toPrecision(4)) : v), 1)))));
      return el;
    };

    const meshTab = () => {
      const el = h('div');
      if (model.kind === 'surface' || model.kind === 'volume') {
        const size = h('input.plain', { type: 'number', step: 'any', min: 0, placeholder: 'target edge length [m]', style: { width: '190px' } });
        const count = h('input.plain', { type: 'number', min: 4, placeholder: 'target triangle count', style: { width: '190px' } });
        el.append(h('div.card', h('h3', 'Refine or coarsen'), h('p.muted', 'Refinement splits elements until edges are no longer than the target. Coarsening merges edges with the least geometric error until the target count is reached.'),
          h('div.row', size, h('button.btn', { type: 'button', onclick: () => { const s = Number(size.value); if (!(s > 0)) { toast('Enter a target edge length.'); return; } try { model = refineM.refineTo(model, s); show(); toast(`Refined: ${statsOf(model)}.`); } catch (e) { toast(e.message); } } }, 'Refine'),
            h('button.btn', { type: 'button', onclick: () => { try { model = refineM.refine(model); show(); toast(`Refined once: ${statsOf(model)}.`); } catch (e) { toast(e.message); } } }, 'Refine ×1 (every element split)')),
          model.kind === 'surface' ? h('div.row', { style: { marginTop: '8px' } }, count, h('button.btn', { type: 'button', onclick: () => { const n = Number(count.value); if (!(n >= 4)) { toast('Enter a target count.'); return; } try { model = refineM.coarsenTo(model, n); show(); toast(`Coarsened: ${statsOf(model)}.`); } catch (e) { toast(e.message); } } }, 'Coarsen'),
            h('button.btn', { type: 'button', onclick: () => { try { model = refineM.smooth(model, { iterations: 5, taubin: true }); show(); toast('Smoothed (volume-preserving).'); } catch (e) { toast(e.message); } } }, 'Smooth')) : null,
          h('div.row', { style: { marginTop: '8px' } }, h('button.btn', { type: 'button', onclick: () => save('remeshed') }, 'Save as new model'))));
      }
      if (model.kind === 'surface') {
        const res = h('div');
        const tol = h('input.plain', { type: 'number', value: 1, min: 0.01, step: 0.1, style: { width: '80px' }, 'aria-label': 'Tolerance in percent' });
        const run = async () => {
          clear(res).append(h('small.muted', 'Running the mesh ladder…'));
          try {
            const base = model.indices.length / 3;
            const levels = [Math.max(24, Math.round(base / 16)), Math.max(48, Math.round(base / 4)), base, base * 4].filter((v, i, a) => i === 0 || v > a[i - 1] * 1.5);
            const built = new Map([[base, model]]);
            const get = (n) => { if (!built.has(n)) built.set(n, n < base ? refineM.coarsenTo(model, n) : refineM.refine(model)); return built.get(n); };
            const rows = [];
            for (const n of levels) {
              const m = get(n);
              const t0 = performance.now();
              const ms = measure.measure(m);
              rows.push({ n: m.indices.length / 3, area: ms.area, volume: Math.abs(ms.volume ?? NaN), ms: performance.now() - t0, model: m });
            }
            const closed = rep.watertight;
            const key = closed ? 'volume' : 'area';
            const f = rows.map((r) => r[key]);
            const hRep = rows.map((r) => 1 / Math.sqrt(r.n));
            const last = rows.length - 1;
            const tg = safe(() => conv.threeGrid({ f1: f[last], f2: f[last - 1], f3: f[last - 2], r21: hRep[last - 1] / hRep[last], r32: hRep[last - 2] / hRep[last - 1] }, {})) || {};
            const ext = tg.extrapolated ?? tg.fExt ?? f[last];
            rows.forEach((r) => { r.err = Math.abs((r[key] - ext) / ext); });
            const pick = rows.find((r) => r.err <= Number(tol.value) / 100) || rows[last];
            clear(res).append(h('div.banner', icon('check'), h('div', h('b', `Use about ${fmt(pick.n)} triangles. `), `At that density the ${closed ? 'enclosed volume' : 'surface area'} is within ${fmt(100 * pick.err, 2)} % of the extrapolated value ${fmt(ext)} ${closed ? 'm³' : 'm²'}`, tg.p != null ? `; observed order ${fmt(tg.p, 3)}` : '', tg.gciFine != null ? `, fine-grid GCI ${fmt(100 * tg.gciFine, 2)} %` : '', '. Finer meshes cost run time without changing the answer.')),
              h('div.grid.c2',
                renderPlot({ type: 'line', title: `${closed ? 'Enclosed volume' : 'Surface area'} vs mesh density`, x: { label: 'Triangles', log: true }, y: { label: closed ? 'Volume' : 'Area', unit: closed ? 'm³' : 'm²', zero: false }, series: [{ name: 'Measured', x: rows.map((r) => r.n), y: f, markers: true }, { name: 'Extrapolated', x: [rows[0].n, rows[last].n], y: [ext, ext], dash: true }] }),
                renderPlot({ type: 'line', title: 'Discretisation error', x: { label: 'Triangles', log: true }, y: { label: 'Relative error', log: true }, series: [{ name: 'Error', x: rows.filter((r) => r.err > 0).map((r) => r.n), y: rows.filter((r) => r.err > 0).map((r) => r.err), markers: true }] })),
              table([{ key: 'n', label: 'Triangles' }, { key: 'area', label: 'Area', unit: 'm²', sig: 7 }, { key: 'volume', label: 'Volume', unit: 'm³', sig: 7 }, { key: 'errp', label: 'Error [%]' }, { key: 'pick', label: '', left: true }], rows.map((r) => ({ ...r, errp: 100 * r.err, pick: r === pick ? '◀ selected' : '' }))),
              h('div.row', { style: { marginTop: '8px' } }, h('button.btn', { type: 'button', onclick: () => { model = pick.model; show(); toast(`Mesh set to ${fmt(pick.n)} triangles.`); } }, 'Use the selected mesh')));
          } catch (e) { clear(res).append(h('div.banner.fail', e.message)); }
        };
        el.append(h('div.card', h('h3', 'Mesh sensitivity and selection'), h('p.muted', 'The model is measured on a ladder of coarser and finer meshes. The coarsest mesh whose result is within your tolerance of the extrapolated value is selected — the same test every numerical analysis module applies to its own discretisation in its Convergence tab.'),
          h('div.row', h('label', 'Tolerance'), tol, h('span', '%'), h('button.btn.primary', { type: 'button', onclick: run }, 'Run mesh study')), res));
      }
      if (model.kind === 'frame') {
        const len = h('input.plain', { type: 'number', step: 'any', min: 0, placeholder: 'max element length [m]', style: { width: '200px' } });
        el.append(h('div.card', h('h3', 'Member subdivision'), h('div.row', len, h('button.btn', { type: 'button', onclick: () => { const L = Number(len.value); if (!(L > 0)) return; model = refineM.subdivideMembers(model, { maxLength: L }); show(); toast(`Subdivided: ${statsOf(model)}.`); } }, 'Subdivide'), h('button.btn', { type: 'button', onclick: () => save('subdivided') }, 'Save as new model'))));
      }
      if (model.kind === 'grid') {
        const gs = safe(() => measure.gridStats(model));
        if (gs) el.append(h('div.card', h('h3', 'Seabed statistics'), h('pre', { style: { fontSize: '.8rem', overflowX: 'auto' } }, JSON.stringify(gs, (k, v) => (typeof v === 'number' ? Number(v.toPrecision(4)) : (ArrayBuffer.isView(v) ? undefined : v)), 1))));
      }
      if (!el.childNodes.length) el.append(h('div.empty', 'No meshing operation applies to this kind of model.'));
      return el;
    };

    const exportTab = () => {
      const el = h('div.card', h('h3', 'Export'), h('p.muted', 'Take the checked, repaired or remeshed model back out.'));
      import('../../geometry/exporters.js').then((ex) => {
        const row = h('div.row');
        for (const e of ex.EXPORTERS.filter((x) => x.kinds.includes(model.kind))) {
          row.append(h('button.btn', { type: 'button', onclick: () => { try { const out = ex.exportModels(e.id, model, { fileName: model.name || entry.name }); download(out.fileName, out.data instanceof Uint8Array ? new Blob([out.data], { type: out.mime }) : out.data, out.mime); } catch (err) { toast(err.message); } } }, icon('download'), e.name));
        }
        el.append(row.childNodes.length ? row : h('p.muted', 'No exporter applies to this kind of model.'));
      }).catch(() => el.append(h('p.muted', 'Exporters are not available in this build.')));
      return el;
    };

    // Hand-off to high-fidelity solvers, and the built-in panel-method added mass
    const solverTab = () => {
      const el = h('div');
      const num = (label, value, unit, w = 90) => { const i = h('input.plain', { type: 'number', step: 'any', value, 'aria-label': label, style: { width: `${w}px` } }); return { i, row: h('span.row', h('label', label), i, unit ? h('small.muted', unit) : null), get v() { return Number(i.value); } }; };
      const bb = topology.bounds ? safe(() => topology.bounds(model.positions || model.nodes)) : null;
      const size = bb ? [bb.max[0] - bb.min[0], bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]] : [10, 10, 10];
      if (model.kind === 'surface') {
        const res = h('div');
        const rho = num('Water density', 1025, 'kg/m³');
        const lim = h('select.plain', { style: { maxWidth: '300px' }, 'aria-label': 'Frequency limit' }, h('option', { value: 'infinite' }, 'High-frequency limit (A∞, for time-domain models)'), h('option', { value: 'zero' }, 'Low-frequency limit (rigid free surface)'), h('option', { value: 'unbounded' }, 'Deeply submerged (no free surface)'));
        el.append(h('div.card', h('h3', 'Added mass by the boundary-element method'), h('p.muted', 'Solves the potential flow around this hull with panels on its wetted surface (waterline at z = 0) and returns the 6 × 6 added-mass matrix. This is a true panel-method result for the stated frequency limit. It is not a frequency-dependent radiation and diffraction analysis: it gives no wave damping and no wave excitation. For those, export a case below and bring the results back.'),
          h('div.row', rho.row, lim, h('button.btn.primary', { type: 'button', onclick: async () => {
            clear(res).append(h('small.muted', 'Solving…'));
            try {
              const { addedMassBem } = await import('../../engine/hydro/bem.js');
              let hull = model;
              if (hull.indices.length / 3 > 2400) hull = refineM.coarsenTo(model, 2400);   // keeps the dense solve to a few seconds on a phone
              await new Promise((r) => setTimeout(r, 30));
              const t0 = performance.now();
              const r = addedMassBem({ surface: hull, rho: rho.v, limit: lim.value });
              const lab = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'];
              clear(res).append(h('p', `${r.panels} wetted panels, displaced volume ${fmt(r.volume)} m³, solved in ${fmt((performance.now() - t0) / 1000, 2)} s. Asymmetry removed from the matrix: ${fmt(100 * r.symmetryError, 2)} % (a measure of discretisation error — refine the mesh if it is more than a few per cent).`),
                h('p.why', 'Rows and columns are surge, sway, heave [kg], roll, pitch, yaw [kg·m²]; off-diagonal blocks couple translation and rotation [kg·m]. Compare heave added mass with the displaced mass ' + fmt(r.volume * rho.v) + ' kg.'),
                table([{ key: 'n', label: '', left: true }, ...lab.map((l, j) => ({ key: `c${j}`, label: l, sig: 4 }))], r.A.map((row, i) => ({ n: lab[i], ...Object.fromEntries(row.map((v, j) => [`c${j}`, v])) }))),
                h('div.row', { style: { marginTop: '8px' } }, h('button.btn.small', { type: 'button', onclick: () => download('added-mass.json', JSON.stringify({ limit: r.limit, rho: rho.v, panels: r.panels, A: r.A }, null, 1), 'application/json') }, icon('download'), 'Matrix (JSON)')));
            } catch (e) { clear(res).append(h('div.banner.fail', icon('warn'), e.message)); }
          } }, 'Compute added mass')), res));
      }
      const card = h('div.card', h('h3', 'Hand off to a high-fidelity solver'), h('p.muted', 'Writes a complete, ready-to-run case for an established open-source solver from this model, as one ZIP. Run it on a workstation or cluster, then bring the results back into the analyses. This is how questions beyond the built-in models — wave diffraction and radiation, nonlinear collapse, viscous flow — are answered without leaving the workflow.'));
      el.append(card);
      import('../../interop/index.js').then((io) => {
        const list = io.SOLVER_EXPORTS.filter((x) => x.kinds.includes(model.kind));
        if (!list.length) { card.append(h('p.muted', 'No solver hand-off applies to this kind of model.')); return; }
        const depth = num('Water depth', Math.max(20, Math.round(-((bb?.min?.[2]) ?? -20) * 3)), 'm'), H = num('Wave height', 2, 'm'), T = num('Wave period', 8, 's'), U = num('Current', 1, 'm/s');
        const E = num('E', 210, 'GPa', 80), fy = num('Yield', 355, 'MPa', 80), thk = num('Shell thickness', 20, 'mm', 80);
        card.append(h('div.row', { style: { margin: '8px 0' } }, depth.row, H.row, T.row, U.row), model.kind !== 'grid' ? h('div.row', { style: { margin: '8px 0' } }, E.row, fy.row, model.kind === 'surface' ? thk.row : null) : null);
        const L = Math.max(...size);
        const options = (id) => {
          const name = (model.name || entry.name).replace(/\.[^.]+$/, '').replace(/[^\w-]+/g, '_') || 'model';
          const omegas = Array.from({ length: 30 }, (_, i) => 0.2 + i * 0.06);
          if (id === 'nemoh' || id === 'capytaine') return { surface: model, rho: 1025, g: 9.80665, depth: depth.v, omegas, headingsDeg: [0, 45, 90, 135, 180], cog: [0, 0, 0], name };
          if (id === 'wamit-gdf') return { surface: model, ulen: 1, g: 9.80665, name };
          if (id === 'calculix') {
            // default restraint: every node on the lowest 1 % of the model height is fixed; default load: self-weight
            const P = model.positions || model.nodes; const zmin = bb.min[2], tol = 0.01 * (size[2] || 1) + 1e-9; const fixed = [];
            for (let i = 0; i < P.length / 3; i++) if (P[3 * i + 2] <= zmin + tol) fixed.push(i);
            return { material: { E: E.v * 1e9, nu: 0.3, rho: 7850, fy: fy.v * 1e6, hardening: { Et: E.v * 1e9 / 100, epsMax: 0.15 } }, supports: [{ nodes: fixed, dofs: 'fixed' }], loads: [{ type: 'gravity' }], analysis: 'static-nlgeom', steps: 10, thickness: thk.v / 1000, name };
          }
          if (id === 'openfoam-wave-tank') return { surface: model, domain: { length: 12 * L, width: 6 * L, height: depth.v + 3 * Math.max(H.v, size[2]), waterDepth: depth.v, xMin: bb.min[0] - 5 * L, yMin: (bb.min[1] + bb.max[1]) / 2 - 3 * L }, wave: { H: H.v, T: T.v, theory: 'stokes2' }, cellSize: L / 6, refinement: { structure: 2, box: 1, freeSurface: 1 }, endTime: 6 * T.v, turbulence: 'laminar', name };
          if (id === 'openfoam-current') return { surface: model, domain: { length: 14 * L, width: 6 * L, height: 6 * L, xMin: bb.min[0] - 4 * L, yMin: (bb.min[1] + bb.max[1]) / 2 - 3 * L, zMin: (bb.min[2] + bb.max[2]) / 2 - 3 * L }, current: U.v, cellSize: L / 6, refinement: { structure: 2, box: 1 }, name };
          return { name };
        };
        card.append(table([{ key: 'solver', label: 'Solver', left: true }, { key: 'purpose', label: 'Answers', left: true }, { key: 'state', label: 'Case files', left: true }, { key: 'act', label: '', left: true }], list.map((x) => ({ solver: x.solver, purpose: h('span', x.purpose, x.notes ? h('small.why', Array.isArray(x.notes) ? x.notes.join(' ') : x.notes) : null),
          state: h('span', { class: `chip ${x.executed ? 'ok' : 'warn'}`, title: x.executed ? 'A case written by this exporter was run end to end with the solver during development.' : 'Written to the solver’s documented input format; not executed during development because the solver was not installed. Check the first case you run.' }, x.executed ? '✓ run with the solver' : 'format-checked only'),
          act: h('button.btn.small', { type: 'button', onclick: () => { try { const out = io.exportCase(x.id, model, options(x.id)); download(`${x.id}-case.zip`, new Blob([io.packCase(out.files, { root: `${x.id}-case` })], { type: 'application/zip' })); toast(`${out.files.length} files written.${out.notes?.length ? ` ${out.notes[0]}` : ''}`.slice(0, 200)); } catch (e) { toast(e.message); } } }, icon('download'), 'Case (ZIP)') }))),
        h('p.why', 'Defaults are filled in from this model’s size (domain extents, cell size, restraint at its base, self-weight). They produce a case that runs; refine mesh, loads and run time in the solver’s own files for a production analysis. Results come back through the floating-structure RAO module (hydrodynamic coefficients) and the readers documented in docs/interop.md.'));
      }).catch((e) => card.append(h('p.muted', `Solver hand-off is not available in this build (${e.message}).`)));
      return el;
    };

    info.append(tabs([{ label: 'Check & repair', render: checkTab }, { label: 'Mesh quality', render: qualityTab }, { label: 'Refine, coarsen & mesh study', render: meshTab }, { label: 'High-fidelity solvers', render: solverTab }, { label: 'Export', render: exportTab }]));
  }
  show();
}

const safe = (fn) => { try { return fn(); } catch { return null; } };
function dataUrlBlob(url) { const [head, b64] = url.split(','); const bin = atob(b64); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return new Blob([u], { type: /:(.*?);/.exec(head)[1] }); }

export default async function geometry(root, params = {}) {
  const list = h('div');
  const status = h('div');
  const bench = h('div');
  let accept = '';
  try { accept = (await import('../../geometry/registry.js')).acceptString(); } catch { /* registry not present in this build */ }

  const onFiles = async (files) => {
    clear(status).append(h('small.muted', 'Importing…'));
    const res = await importFiles(files, (m) => { status.firstChild.textContent = m; });
    clear(status);
    for (const r of res) status.append(h('div.banner', { class: r.ok ? '' : 'fail' }, icon(r.ok ? 'check' : 'warn'), h('div', h('b', r.file), r.ok ? ` — imported (${r.models} model${r.models > 1 ? 's' : ''}).${r.warnings.length ? ` ${r.warnings.length} note${r.warnings.length > 1 ? 's' : ''}: ${r.warnings.slice(0, 2).join(' ')}` : ''}` : ` — ${r.error}`)));
    draw();
    const first = res.find((r) => r.ok);
    if (first) open(first.id);
  };
  const open = (id) => { const e = getProject().geometry.find((g) => g.id === id); if (e) { bench.scrollIntoView({ behavior: 'smooth', block: 'start' }); clear(bench).append(h('h2', e.name)); const b = h('div'); bench.append(b); workbench(e, b).catch((err) => clear(b).append(h('div.banner.fail', err.message))); } };
  function draw() {
    const g = getProject().geometry;
    clear(list);
    if (!g.length) { list.append(h('div.empty', 'No geometry in this project yet.')); return; }
    list.append(table([{ key: 'name', label: 'Model', left: true }, { key: 'kind', label: 'Kind', left: true }, { key: 'format', label: 'Format', left: true }, { key: 'stats', label: 'Content', left: true }, { key: 'act', label: '', left: true }],
      g.map((e) => ({ name: h('a', { href: `#/geometry/${e.id}` }, e.name), kind: (e.kinds || [e.kind]).map((k) => KIND_LABEL[k] || k).join(', '), format: e.format || '', stats: e.stats || '',
        act: h('span.row', h('button.btn.small', { type: 'button', onclick: () => open(e.id) }, 'Open'), h('button.btn.small', { type: 'button', 'aria-label': `Remove ${e.name}`, onclick: async () => { if (confirm(`Remove “${e.name}” from this project?`)) { await removeGeometry(e.id); draw(); clear(bench); } } }, icon('trash'))) }))));
  }

  root.append(h('h1', 'Geometry & mesh'), h('p.lede', 'Bring in geometry from your CAD, structural, survey or meshing system. It is recognised by content, converted to metres, checked, and made selectable inside every analysis that can use it. This is a preparation and meshing bench, not a modelling tool.'),
    dropZone(onFiles, accept, 'Drop CAD, mesh, structural, route, bathymetry or point-cloud files here'), status,
    h('h2', { style: { marginTop: '20px' } }, 'Models in this project'), list, bench,
    h('h2', { style: { marginTop: '24px' } }, 'Supported formats'), await formatMatrix().catch(() => h('div.empty', 'The format registry is not available in this build.')),
    pager({ href: '#/portal', label: 'Input data portal' }, { href: '#/materials', label: 'Materials library' }));
  draw();
  if (params.id) open(params.id);
  return { title: 'Geometry & mesh', crumbs: [{ href: '#/geometry', label: 'Geometry & mesh' }] };
}
