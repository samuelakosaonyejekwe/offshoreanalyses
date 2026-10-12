// Code-factor register: every standards-derived factor used by any module, in one place, with a sign-off.
// Built-in values were written from the standards as recalled; a factor counts as confirmed only when someone
// has checked it against the governing edition and said so here.
import { h, clear, fmt, icon, toast, table, download } from '../dom.js';
import { suiteIndex, suiteModules, SHORT } from '../../core/registry.js';
import { getProject, confirmFactor, save, emit, prefs } from '../../core/store.js';
import { pager } from '../shell.js';

export default async function factors(root, params = {}) {
  const idx = await suiteIndex();
  const mods = [];
  for (const s of idx) for (const m of await suiteModules(s.n)) if (m.criteria?.length) mods.push(m);
  const total = mods.reduce((n, m) => n + m.criteria.length, 0);
  const p = getProject();
  p.factors ||= {};
  const summary = h('div.kpis');
  const body = h('div');
  const who = h('input.plain', { type: 'text', value: prefs.get('checker', ''), placeholder: 'Your name or initials', 'aria-label': 'Checker name', style: { maxWidth: '220px' } });
  who.addEventListener('change', () => prefs.set('checker', who.value));
  const only = h('input', { type: 'checkbox', checked: false });
  const find = h('input.plain', { type: 'search', placeholder: 'Filter by factor, module or standard…', 'aria-label': 'Filter', style: { maxWidth: '300px' } });
  let suiteSel = params.suite ? Number(params.suite) : 0;
  const suitePick = h('select.plain', { style: { maxWidth: '260px' }, 'aria-label': 'Suite', onchange: (e) => { suiteSel = Number(e.target.value); draw(); } },
    h('option', { value: 0 }, 'All suites'), idx.map((s) => h('option', { value: s.n, selected: s.n === suiteSel }, `${s.n}. ${SHORT[s.n]}`)));

  const drawSummary = () => {
    const done = Object.keys(p.factors).length;
    clear(summary).append(
      h('div.kpi', h('span', 'Standards-derived factors'), h('b', total), h('small.muted', `in ${mods.length} modules`)),
      h('div.kpi', { class: done === total ? 'ok' : '' }, h('span', 'Confirmed against the governing edition'), h('b', done), h('small.muted', `${fmt(100 * done / total, 3)} %`)),
      h('div.kpi', { class: total - done ? 'fail' : 'ok' }, h('span', 'Still built-in defaults'), h('b', total - done), h('small.muted', 'check before design use')));
  };

  function draw() {
    clear(body);
    const q = find.value.trim().toLowerCase();
    const shown = mods.filter((m) => (!suiteSel || m.suite === suiteSel) && (!params.m || m.id === params.m));
    let rendered = 0;
    for (const m of shown) {
      const rows = m.criteria.filter((c) => {
        const conf = p.factors[`${m.id}.${c.key}`];
        if (only.checked && conf) return false;
        return !q || `${c.label} ${m.title} ${(m.standards || []).join(' ')}`.toLowerCase().includes(q);
      });
      if (!rows.length) continue;
      if (++rendered > 40 && !q && !params.m) { body.append(h('p.muted', `Showing the first 40 modules. Choose a suite or type a filter to see the rest of the ${shown.length}.`)); break; }
      const edition = h('input.plain', { type: 'text', placeholder: 'Edition checked, e.g. “ISO 19902:2020”', 'aria-label': `Edition for ${m.title}`, style: { maxWidth: '280px' }, value: (m.standards || [])[0] || '' });
      const card = h('div.card', { style: { marginTop: '12px' } },
        h('div.row.spread', h('h3', { style: { margin: 0 } }, h('a', { href: `#/m/${m.id}` }, m.title), h('small.muted', `  suite ${m.suite}`)),
          h('div.row', edition, h('button.btn.small', { type: 'button', title: 'Record every factor shown for this module, at its current value, as checked against the edition entered', onclick: () => {
            if (!who.value.trim()) { toast('Enter your name or initials first — a sign-off needs a person.'); who.focus(); return; }
            if (!edition.value.trim()) { toast('Enter the edition you checked against.'); edition.focus(); return; }
            for (const c of rows) { const inp = card.querySelector(`[data-k="${c.key}"]`); const v = Number(inp.value); if (Number.isFinite(v)) confirmFactor(m.id, c.key, { value: v, edition: edition.value.trim(), by: who.value.trim() }); }
            drawSummary(); draw(); toast(`${rows.length} factor${rows.length > 1 ? 's' : ''} of “${m.title}” recorded as confirmed.`);
          } }, icon('check'), 'Confirm all shown'))),
        h('small.muted', (m.standards || []).join(' · ')),
        table([{ key: 'label', label: 'Factor', left: true }, { key: 'def', label: 'Built-in' }, { key: 'val', label: 'Value to use' }, { key: 'unit', label: 'Unit', left: true }, { key: 'st', label: 'Status', left: true }, { key: 'act', label: '', left: true }],
          rows.map((c) => {
            const conf = p.factors[`${m.id}.${c.key}`];
            const inp = h('input', { type: 'text', inputmode: 'decimal', 'data-k': c.key, value: conf ? conf.value : c.default, 'aria-label': c.label, style: { width: '110px', textAlign: 'right', border: '1px solid var(--line)', borderRadius: '6px', padding: '3px 6px', background: 'var(--surface)', color: 'var(--ink)' } });
            return { label: h('span', c.label, c.help ? h('small.why', c.help) : null), def: c.default, val: inp, unit: c.unit,
              st: conf ? h('span.chip.ok', { title: `${conf.edition} — ${conf.by}, ${new Date(conf.at).toLocaleDateString()}` }, `✓ ${conf.edition} · ${conf.by}`) : h('span.chip.warn', 'built-in default, unconfirmed'),
              act: conf ? h('button.btn.small', { type: 'button', onclick: () => { confirmFactor(m.id, c.key, null); drawSummary(); draw(); } }, 'Withdraw')
                : h('button.btn.small', { type: 'button', onclick: () => {
                  if (!who.value.trim()) { toast('Enter your name or initials first — a sign-off needs a person.'); who.focus(); return; }
                  if (!edition.value.trim()) { toast('Enter the edition you checked against.'); edition.focus(); return; }
                  const v = Number(inp.value); if (!Number.isFinite(v)) { toast('Enter a number.'); return; }
                  confirmFactor(m.id, c.key, { value: v, edition: edition.value.trim(), by: who.value.trim() }); drawSummary(); draw();
                } }, 'Confirm') };
          })));
      body.append(card);
    }
    if (!rendered) body.append(h('div.empty', only.checked ? 'Nothing left to confirm in this selection.' : 'No factor matches the filter.'));
  }
  only.addEventListener('change', draw);
  find.addEventListener('input', draw);

  const file = h('input', { type: 'file', accept: '.json', style: { display: 'none' }, onchange: async () => {
    try {
      const obj = JSON.parse(await file.files[0].text());
      if (obj.kind !== 'offshore-analyses-code-factor-profile' || typeof obj.factors !== 'object') throw new Error('not a profile');
      let n = 0;
      for (const [k, v] of Object.entries(obj.factors)) if (Number.isFinite(v.value) && v.edition && v.by) { p.factors[k] = v; n++; }
      save(); emit('factors', p.factors); drawSummary(); draw(); toast(`${n} confirmed factors loaded from “${obj.name || 'profile'}”.`);
    } catch { toast('That file is not a code-factor profile.'); }
    file.value = '';
  } });

  root.append(h('h1', 'Code-factor register'),
    h('p.lede', 'Partial factors, utilisation limits, model coefficients and other numbers taken from standards are editable inputs in every module. Their built-in values were written from the standards as recalled, not checked against a licensed copy. This register is where each one gets checked once, against the edition that governs your project, and signed off.'),
    h('div.banner.warn', icon('warn'), h('div', h('b', 'Until a factor is confirmed here, treat it as unverified. '), 'Every analysis result states how many of the factors it used are unconfirmed. A confirmed value replaces the built-in default in that module from then on, and a confirmed set can be exported as a profile so one check serves a whole team.')),
    summary,
    h('div.row', { style: { margin: '14px 0 4px' } }, h('label', 'Checked by'), who, suitePick, find, h('label.check', only, h('span', 'Unconfirmed only')),
      h('button.btn', { type: 'button', onclick: () => download('code-factor-profile.json', JSON.stringify({ kind: 'offshore-analyses-code-factor-profile', name: p.name, exported: new Date().toISOString(), factors: p.factors }, null, 1), 'application/json') }, icon('download'), 'Export profile'),
      file, h('button.btn', { type: 'button', onclick: () => file.click() }, icon('upload'), 'Load profile')),
    params.m ? h('p', h('a', { href: '#/factors' }, '← All modules')) : null,
    body,
    pager({ href: '#/verification', label: 'Verification & validation' }, { href: '#/wiring', label: 'Suite wiring & coverage' }));
  drawSummary(); draw();
  return { title: 'Code-factor register', crumbs: [{ href: '#/factors', label: 'Code-factor register' }] };
}
