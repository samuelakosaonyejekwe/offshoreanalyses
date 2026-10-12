// Renders a module result (docs/CONTRACT.md §2.2): headline numbers, acceptance checks, charts, 3-D views,
// tables, remarks and recommendations, plus the exact set of conditions the run used.
import { h, fmt, table, icon, download, toCsv } from './dom.js';
import { renderPlot } from './plot.js';
import { kpiOk, summarise, conditionReport, ROLE_LABEL } from '../core/studies.js';

export function kpiCards(result) {
  return h('div.kpis', (result.kpis || []).map((k) => {
    const ok = kpiOk(k);
    return h('div.kpi', { class: ok == null ? '' : ok ? 'ok' : 'fail' },
      h('span', k.label), h('b', fmt(k.value)), k.unit && k.unit !== '-' ? h('em', k.unit) : null,
      k.limit != null ? h('small', { class: ok ? 'muted' : '' }, `${ok ? '✓ within' : '✗ outside'} limit ${k.sense === 'min' ? '≥' : '≤'} ${fmt(k.limit)}`) : null,
      k.explain ? h('small.why', k.explain) : null);
  }));
}

export function checksTable(result) {
  const checks = result.checks || [];
  if (!checks.length) return null;
  const sorted = checks.slice().sort((a, b) => (b.uc ?? 0) - (a.uc ?? 0));
  return h('div.tablewrap', { style: { maxHeight: '380px', overflowY: 'auto' } }, h('table.t',
    h('thead', h('tr', h('th.l', 'Check'), h('th', 'Demand'), h('th', 'Capacity'), h('th', 'Utilisation'), h('th.l', 'Reference'))),
    h('tbody', sorted.map((c) => {
      const cls = c.uc > 1 ? 'fail' : c.uc > 0.9 ? 'warn' : '';
      return h('tr',
        h('td.l', c.name, c.explain ? h('small.why', c.explain) : null), h('td', fmt(c.demand), c.unit ? h('small.muted', ` ${c.unit}`) : null), h('td', fmt(c.capacity), c.unit ? h('small.muted', ` ${c.unit}`) : null),
        h('td', h('span', { class: `ucbar ${cls}` }, h('i', { style: { width: `${Math.min(100, Math.max(2, (c.uc || 0) * 100))}%` } })), h('b', fmt(c.uc, 3)), ' ', h('span', { class: `chip ${c.uc > 1 ? 'fail' : 'ok'}` }, c.uc > 1 ? '✗ fails' : '✓ passes')),
        h('td.l', h('small', c.ref || '')));
    }))));
}

function modelView(p) {
  const box = h('div.viewer');
  const wrap = h('div.figure', h('h3', p.title || '3-D view'), box, h('small.muted', 'Drag to orbit · two fingers or right-drag to pan · wheel or pinch to zoom · F fits the view.'));
  import('./viewer3d.js').then(({ createViewer }) => {
    const viewer = createViewer(box, {});
    viewer.setScene({ models: [{ kind: 'model3d', ...p }] });
    viewer.fit();
  }).catch(() => { box.append(h('div.empty', '3-D view is not available on this device.')); });
  return wrap;
}

export function adviceList(items) {
  if (!items?.length) return null;
  const rank = { fail: 0, warn: 1, info: 2, ok: 3 };
  return h('div', items.slice().sort((a, b) => (rank[a.severity] ?? 2) - (rank[b.severity] ?? 2)).map((a) =>
    h('div.advice', { class: a.severity || 'info' }, h('span.ico', { 'aria-hidden': 'true' }, a.severity === 'fail' ? '!' : a.severity === 'warn' ? '!' : a.severity === 'ok' ? '✓' : 'i'),
      h('div', a.kind ? h('span', { class: `chip ${a.kind === 'economic' ? 'warn' : 'info'}`, style: { marginRight: '6px' } }, a.kind === 'economic' ? 'Economic practice' : 'Engineering practice') : null, h('b', a.text), a.action ? h('p', a.action) : null,
        a.sustainability ? h('p.sus', h('span', { 'aria-hidden': 'true' }, '♻ '), h('b', 'Sustainability: '), a.sustainability) : null,
        a.basis ? h('p', h('small', h('b', 'Basis: '), a.basis)) : null))));
}

export function conditionsTable(mod, values) {
  const rows = conditionReport(mod, values);
  const order = ['ic', 'bc', 'geometry', 'material', 'soil', 'environment', 'load', 'operation', 'criteria', 'control', 'cost'];
  rows.sort((a, b) => order.indexOf(a.role) - order.indexOf(b.role));
  return table([{ key: 'roleL', label: 'Type', left: true }, { key: 'label', label: 'Condition / input', left: true }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit', left: true }],
    rows.map((r) => ({ ...r, roleL: ROLE_LABEL[r.role] || r.role })));
}

/** Full result panel. */
export function renderResult(mod, values, result, { elapsed } = {}) {
  const s = summarise(result);
  const el = h('div');
  el.append(h('div.banner', { class: s.status === 'fail' ? 'fail' : '' }, icon(s.status === 'fail' ? 'warn' : 'check'),
    h('div', h('b', s.status === 'fail' ? `${s.fails} criterion${s.fails > 1 ? 'a are' : ' is'} not met` : (result.checks?.length || (result.kpis || []).some((k) => k.limit != null) ? 'All criteria are met' : 'Calculation complete')),
      s.maxUC ? ` — highest utilisation ${fmt(s.maxUC, 3)}${s.governing ? ` (${s.governing})` : ''}` : '', elapsed != null ? h('small.muted', ` · ${elapsed < 1000 ? `${Math.round(elapsed)} ms` : `${(elapsed / 1000).toFixed(1)} s`}`) : null)));
  if (result.kpis?.length) el.append(kpiCards(result));
  const adv = adviceList(result.advice);
  if (adv) el.append(h('h3', { style: { marginTop: '18px' } }, 'What this means and what to do'), adv);
  const ct = checksTable(result);
  if (ct) el.append(h('h3', { style: { marginTop: '18px' } }, 'Acceptance checks'), ct);
  if (result.plots?.length) {
    el.append(h('h3', { style: { marginTop: '18px' } }, 'Charts'));
    el.append(h('div.grid.c2', result.plots.map((p) => { const f = p.type === 'model3d' ? modelView(p) : renderPlot(p); if (p.explain) f.append(h('p.why', p.explain)); return f; })));
  }
  for (const t of result.tables || []) {
    el.append(h('div.row.spread', { style: { marginTop: '18px' } }, h('h3', { style: { margin: 0 } }, t.title || 'Table'),
      h('button.btn.small', { type: 'button', onclick: () => download(`${(t.title || 'table').replace(/[^\w-]+/g, '_')}.csv`, toCsv(t.columns, t.rows), 'text/csv') }, icon('download'), 'CSV')));
    if (t.explain) el.append(h('p.why', t.explain));
    el.append(table(t.columns, t.rows));
  }
  if (result.notes?.length) el.append(h('h3', { style: { marginTop: '18px' } }, 'Remarks'), h('ul.prose', result.notes.map((n) => h('li', n))));
  el.append(h('details.fold', h('summary', 'Initial conditions, boundary conditions and inputs used for this run'), h('div', conditionsTable(mod, values))));
  if (result.publish && Object.keys(result.publish).length) {
    el.append(h('details.fold', h('summary', `Values passed to other suites (${Object.keys(result.publish).length})`), h('div',
      table([{ key: 'k', label: 'Project value', left: true }, { key: 'v', label: 'Value (SI)' }], Object.entries(result.publish).map(([k, v]) => ({ k, v }))))));
  }
  return el;
}

/** Self-contained JSON record of a run: enough to reproduce it and to audit it. */
export function runRecord(mod, values, result) {
  const strip = (v) => (v && typeof v === 'object' && v.kind && (v.positions || v.nodes || v.points || v.z) ? { kind: v.kind, name: v.name, note: 'imported geometry omitted from the record' } : v);
  return {
    application: 'Offshore Analyses', module: { id: mod.id, title: mod.title, suite: mod.suite, fidelity: mod.fidelity, standards: mod.standards, models: mod.models, limits: mod.limits },
    generated: new Date().toISOString(), inputs: Object.fromEntries(Object.entries(values).map(([k, v]) => [k, strip(v)])), conditions: conditionReport(mod, values),
    explanations: 'Every kpi, check, plot and table carries an `explain` field; advice items carry `basis` and `sustainability`.',
    result: { kpis: result.kpis, checks: result.checks, tables: result.tables, notes: result.notes, advice: result.advice, publish: result.publish, quantities: result.quantities },
  };
}
