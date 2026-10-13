// Input form generated from a module's declared inputs: grouped fields with units, role tags (IC/BC…),
// live validation, project-bus links, editable tables and series, material pickers and geometry pickers.
import { h, clear, icon, fmt, toast } from './dom.js';
import { validate, hidden, fromSI, ROLE_LABEL } from '../core/studies.js';
import { materialsIn, material } from '../library/materials.js';
import { getProject, loadGeometry } from '../core/store.js';

/**
 * Inputs of a module that carry a value of the project's screening soil profile: bound to a project value that
 * was published by the profile and still equal to it. Found from the values themselves, so the flag survives a
 * reload and a saved input set.
 * @returns {{ key: string, label: string, bind: string }[]}
 */
export function screeningInputs(mod, values) {
  const bus = getProject().bus || {};
  return (mod.inputs || []).filter((inp) => {
    const b = inp.bind && bus[inp.bind];
    if (!b || b.from !== 'screening') return false;
    const v = fromSI(b.value, inp.unit);
    return Number.isFinite(values[inp.key]) && Math.abs(values[inp.key] - v) <= 1e-9 * Math.max(1, Math.abs(v));
  }).map((inp) => ({ key: inp.key, label: inp.label, bind: inp.bind }));
}

const ROLE_SHORT = { ic: 'IC', bc: 'BC', criteria: 'limit', control: 'solver' };

function parseNum(s) {
  if (typeof s === 'number') return s;
  const t = String(s).trim().replace(/\s/g, '').replace(',', '.');
  if (t === '') return NaN;
  return Number(t);
}

function numberField(inp, state, onChange) {
  const input = h('input', { type: 'text', inputmode: 'decimal', id: `f-${inp.key}`, value: state.values[inp.key] ?? '', autocomplete: 'off', 'aria-describedby': `h-${inp.key}` });
  input.addEventListener('input', () => { state.values[inp.key] = parseNum(input.value); delete state.linked[inp.key]; onChange(inp.key); });
  input.addEventListener('keydown', (e) => {
    // arrow keys nudge by 1 % (10 % with Shift): quick what-if without retyping
    if (e.key !== 'ArrowUp' && e.key !== 'ArrowDown') return;
    const v = parseNum(input.value); if (!Number.isFinite(v)) return;
    e.preventDefault();
    const step = (Math.abs(v) || 1) * (e.shiftKey ? 0.1 : 0.01) * (e.key === 'ArrowUp' ? 1 : -1);
    const nv = Number((v + step).toPrecision(6));
    input.value = nv; state.values[inp.key] = inp.type === 'int' ? Math.round(nv) : nv; delete state.linked[inp.key]; onChange(inp.key);
  });
  return h('div.inp', input, inp.unit && inp.unit !== '-' ? h('span.unit', inp.unit) : null);
}

function selectField(inp, state, onChange) {
  let options = inp.options || [];
  if (inp.source?.startsWith('materials:')) options = materialsIn(inp.source.split(':')[1]).map((m) => ({ value: m.id, label: m.name }));
  const sel = h('select', { id: `f-${inp.key}` }, options.map((o) => h('option', { value: o.value, selected: o.value === state.values[inp.key] }, o.label)));
  sel.addEventListener('change', () => {
    const raw = sel.value; const opt = options.find((o) => String(o.value) === raw);
    state.values[inp.key] = opt ? opt.value : raw; onChange(inp.key, true);
  });
  return h('div.inp', sel);
}

function tableField(inp, state, onChange) {
  const wrap = h('div');
  const draw = () => {
    clear(wrap);
    const rows = state.values[inp.key] || [];
    const cell = (row, col) => {
      if (col.type === 'select') {
        const s = h('select', { 'aria-label': col.label }, (col.options || []).map((o) => h('option', { value: o.value, selected: o.value === row[col.key] }, o.label)));
        s.addEventListener('change', () => { row[col.key] = s.value; onChange(inp.key); });
        return s;
      }
      const i = h('input', { type: 'text', inputmode: col.type === 'text' ? 'text' : 'decimal', value: row[col.key] ?? '', 'aria-label': col.label, style: col.type === 'text' ? { textAlign: 'left' } : null });
      i.addEventListener('input', () => { row[col.key] = col.type === 'text' ? i.value : parseNum(i.value); onChange(inp.key); });
      return i;
    };
    wrap.append(h('div.tablewrap', h('table.t',
      h('thead', h('tr', inp.columns.map((c) => h('th', c.label, c.unit ? h('small', ` [${c.unit}]`) : null)), h('th', ''))),
      h('tbody', rows.map((row, ri) => h('tr', inp.columns.map((c) => h('td', cell(row, c))),
        h('td', h('button.btn.small', { type: 'button', title: 'Remove this row', 'aria-label': 'Remove row', onclick: () => { rows.splice(ri, 1); onChange(inp.key); draw(); } }, '×'))))))),
    h('div.row', { style: { marginTop: '6px' } },
      h('button.btn.small', { type: 'button', onclick: () => { rows.push({ ...(rows[rows.length - 1] || Object.fromEntries(inp.columns.map((c) => [c.key, c.type === 'text' ? '' : 0]))) }); state.values[inp.key] = rows; onChange(inp.key); draw(); } }, '+ Add row'),
      h('button.btn.small', { type: 'button', title: 'Paste rows copied from a spreadsheet (tab- or comma-separated, same column order)', onclick: () => pasteRows(inp, state, () => { onChange(inp.key); draw(); }) }, 'Paste from spreadsheet')));
  };
  draw();
  return wrap;
}

function pasteRows(inp, state, done) {
  const ta = h('textarea.plain', { rows: 8, placeholder: inp.columns.map((c) => c.label).join('\t') });
  const dlg = h('dialog', h('h3', `Paste rows — ${inp.label}`), h('p.muted', `Columns in order: ${inp.columns.map((c) => c.label + (c.unit ? ` [${c.unit}]` : '')).join(', ')}. A header line is skipped automatically.`), ta,
    h('div.row', { style: { marginTop: '10px', justifyContent: 'flex-end' } },
      h('button.btn', { type: 'button', onclick: () => dlg.close() }, 'Cancel'),
      h('button.btn.primary', { type: 'button', onclick: () => {
        const rows = ta.value.trim().split(/\r?\n/).map((l) => l.split(/\t|,|;/).map((s) => s.trim())).filter((r) => r.length >= 1 && r.some(Boolean));
        const numericCols = inp.columns.filter((c) => c.type !== 'text' && c.type !== 'select').length;
        const body = rows.filter((r) => !numericCols || r.some((c, i) => inp.columns[i] && inp.columns[i].type !== 'text' && Number.isFinite(parseNum(c))));
        if (!body.length) { toast('No numeric rows found in the pasted text.'); return; }
        state.values[inp.key] = body.map((r) => Object.fromEntries(inp.columns.map((c, i) => [c.key, c.type === 'text' || c.type === 'select' ? (r[i] ?? '') : parseNum(r[i] ?? '')])));
        dlg.close(); done(); toast(`${body.length} rows imported.`);
      } }, 'Import')));
  document.body.append(dlg); dlg.addEventListener('close', () => dlg.remove()); dlg.showModal(); ta.focus();
}

function seriesField(inp, state, onChange) {
  const asText = (s) => (s?.x || []).map((x, i) => `${x}\t${s.y[i]}`).join('\n');
  const ta = h('textarea.plain', { id: `f-${inp.key}`, rows: 6, 'aria-describedby': `h-${inp.key}` }, asText(state.values[inp.key]));
  const count = h('small.muted');
  const parse = () => {
    const x = [], y = [];
    for (const line of ta.value.split(/\r?\n/)) { const p = line.trim().split(/[\t,; ]+/).map(parseNum); if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1])) { x.push(p[0]); y.push(p[1]); } }
    state.values[inp.key] = { x, y }; count.textContent = `${x.length} points`; onChange(inp.key);
  };
  ta.addEventListener('input', parse);
  count.textContent = `${(state.values[inp.key]?.x || []).length} points`;
  const file = h('input', { type: 'file', accept: '.csv,.txt,.tsv,.dat', style: { display: 'none' }, onchange: async () => { if (file.files[0]) { ta.value = await file.files[0].text(); parse(); } } });
  return h('div', h('small.muted', `Two columns: ${inp.xLabel || 'x'}${inp.xUnit ? ` [${inp.xUnit}]` : ''}, ${inp.yLabel || 'y'}${inp.yUnit ? ` [${inp.yUnit}]` : ''} — type, paste from a spreadsheet, or load a file.`), ta,
    h('div.row', { style: { marginTop: '4px' } }, count, file, h('button.btn.small', { type: 'button', onclick: () => file.click() }, icon('upload'), 'Load file')));
}

function geometryField(inp, state, onChange) {
  const geo = getProject().geometry.filter((g) => !inp.kinds || (g.kinds || [g.kind]).some((k) => inp.kinds.includes(k)));
  const sel = h('select', { id: `f-${inp.key}` },
    h('option', { value: '' }, 'Parametric model (from the fields in this form)'),
    geo.map((g) => h('option', { value: g.id, selected: state.geoIds?.[inp.key] === g.id }, `${g.name} — ${g.kind}${g.stats ? `, ${g.stats}` : ''}`)));
  sel.addEventListener('change', async () => {
    state.geoIds ||= {};
    if (!sel.value) { state.values[inp.key] = null; delete state.geoIds[inp.key]; onChange(inp.key, true); return; }
    const models = await loadGeometry(sel.value);
    const model = (models || []).find((m) => !inp.kinds || inp.kinds.includes(m.kind));
    if (!model) { toast('That file holds no model of the kind this analysis needs.'); sel.value = ''; return; }
    state.values[inp.key] = model; state.geoIds[inp.key] = sel.value; onChange(inp.key, true);
  });
  return h('div', h('div.inp', sel),
    h('small.hint', geo.length ? `Needs: ${inp.kinds?.join(' / ') || 'any'} model.` : h('span', `No imported ${inp.kinds?.join(' / ') || ''} model in this project yet — `, h('a', { href: '#/geometry' }, 'import one in Geometry & mesh'), '.')));
}

/**
 * Build the form. Returns { el, state, refresh } where state.values is the live input object.
 * onChange(key, structural) fires after every edit; `structural` means visibility of other fields may change.
 */
export function buildForm(mod, initial, { onChange, onLinkClick } = {}) {
  const state = { values: initial.values, linked: { ...(initial.linked || {}) }, errors: {}, geoIds: {} };
  const el = h('form.form', { novalidate: true, onsubmit: (e) => e.preventDefault() });
  const fieldEls = new Map();

  const mark = () => {
    state.errors = validate(mod, state.values);
    const screening = new Set(screeningInputs(mod, state.values).map((s) => s.key));
    for (const [key, f] of fieldEls) {
      const err = state.errors[key];
      f.wrap.classList.toggle('bad', !!err);
      f.wrap.classList.toggle('linked', !!state.linked[key]);
      f.wrap.classList.toggle('screening', screening.has(key));
      const inp = f.inp;
      f.hint.textContent = err || (screening.has(key) ? `Screening value — no site investigation yet (project value ${inp.bind}). Replace it with site data.`
        : state.linked[key] ? `Linked to project value ${state.linked[key]} — edit to override.` : (inp.help || rangeHint(inp)));
    }
  };
  const changed = (key, structural) => { if (structural) draw(); else mark(); onChange?.(key, state); };

  function draw() {
    clear(el); fieldEls.clear();
    const groups = new Map();
    for (const inp of mod.inputs || []) {
      if (hidden(inp, state.values)) continue;
      const g = inp.group || 'Inputs';
      if (!groups.has(g)) groups.set(g, []);
      groups.get(g).push(inp);
    }
    for (const [name, inputs] of groups) {
      const fields = h('div.fields');
      for (const inp of inputs) {
        const type = inp.type || 'number';
        let control;
        if (type === 'number' || type === 'int') control = numberField(inp, state, changed);
        else if (type === 'select') control = selectField(inp, state, changed);
        else if (type === 'bool') {
          const cb = h('input', { type: 'checkbox', id: `f-${inp.key}`, checked: !!state.values[inp.key] });
          cb.addEventListener('change', () => { state.values[inp.key] = cb.checked; changed(inp.key, true); });
          control = h('div.check', cb, h('span', inp.checkLabel || 'Yes'));
        } else if (type === 'text') {
          const t = h('input', { type: 'text', id: `f-${inp.key}`, value: state.values[inp.key] ?? '' });
          t.addEventListener('input', () => { state.values[inp.key] = t.value; changed(inp.key); });
          control = h('div.inp', t);
        } else if (type === 'table') control = tableField(inp, state, changed);
        else if (type === 'series') control = seriesField(inp, state, changed);
        else if (type === 'geometry') control = geometryField(inp, state, changed);
        const hint = h('span.hint', { id: `h-${inp.key}` });
        const role = inp.role && inp.role !== 'operation' ? h('span', { class: `role ${inp.role}`, title: ROLE_LABEL[inp.role] }, ROLE_SHORT[inp.role] || inp.role) : null;
        const matInfo = inp.source?.startsWith('materials:') ? h('span.hint', materialLine(state.values[inp.key])) : null;
        const wrap = h('div.field', { class: ['table', 'series', 'geometry'].includes(type) ? 'wide' : '' },
          h('label', { for: `f-${inp.key}` }, inp.label, role, inp.design ? h('span.role', { title: 'Design variable: the advisor may vary this to satisfy failed checks' }, 'design') : null),
          control, type === 'geometry' ? null : hint, matInfo);
        if (type !== 'geometry') fieldEls.set(inp.key, { wrap, hint, inp });
        fields.append(wrap);
      }
      el.append(h('fieldset', h('legend', name), fields));
    }
    mark();
  }
  draw();
  return { el, state, refresh: draw, isValid: () => Object.keys(validate(mod, state.values)).length === 0 };
}

function rangeHint(inp) {
  if (inp.type && inp.type !== 'number' && inp.type !== 'int') return '';
  if (inp.min != null && inp.max != null) return `Range ${fmt(inp.min)} – ${fmt(inp.max)}${inp.unit && inp.unit !== '-' ? ' ' + inp.unit : ''}`;
  if (inp.min != null) return `≥ ${fmt(inp.min)}`;
  if (inp.max != null) return `≤ ${fmt(inp.max)}`;
  return '';
}
function materialLine(id) {
  const m = material(id);
  if (!m) return '';
  const bits = [];
  if (m.fy) bits.push(`fy ${fmt(m.fy / 1e6)} MPa`);
  if (m.fu) bits.push(`fu ${fmt(m.fu / 1e6)} MPa`);
  if (m.fck) bits.push(`fck ${fmt(m.fck / 1e6)} MPa`);
  if (m.E) bits.push(`E ${fmt(m.E / 1e9)} GPa`);
  if (m.rho) bits.push(`ρ ${fmt(m.rho)} kg/m³`);
  if (m.su) bits.push(`su ${fmt(m.su / 1e3)} kPa`);
  if (m.phiDeg) bits.push(`φ′ ${m.phiDeg}°`);
  return `${bits.join(' · ')}${m.src ? ` — ${m.src}` : ''}`;
}
