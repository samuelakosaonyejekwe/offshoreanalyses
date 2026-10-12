// Analysis module page: inputs → results, plus the studies every module gets for free through the contract.
import { h, clear, fmt, icon, toast, download, tabs, table, debounce } from '../dom.js';
import { buildForm } from '../form.js';
import { renderResult, runRecord, adviceList } from '../results.js';
import { renderPlot } from '../plot.js';
import { task } from '../../core/runner.js';
import { defaults, summarise, numericInputs, factorReport, ROLE_LABEL } from '../../core/studies.js';
import { findModule, suiteModules, suiteIndex, SHORT } from '../../core/registry.js';
import { getProject, busValues, setInputs, recordRun, publish, prefs, confirmedFactors, recordEvidence, evidenceOf } from '../../core/store.js';
import { literature } from '../../data/connectors.js';
import { pager } from '../shell.js';

const FIDELITY = {
  analytical: ['Closed-form', 'Closed-form or code-equation calculation.'], numerical: ['Numerical', 'Discretised numerical solution; check convergence in the Convergence tab.'],
  'reduced-order': ['Reduced-order', 'Simplified physical model that captures the governing behaviour; confirm critical cases with a higher-fidelity model.'],
  screening: ['Screening', 'Conservative screening estimate for early design and ranking.'],
};

const progressBar = () => { const bar = h('i'); const msg = h('small.muted'); const el = h('div', h('div.progress', bar), msg); return { el, set: (f, m) => { bar.style.width = `${Math.round(100 * Math.min(1, f || 0))}%`; if (m) msg.textContent = m; } }; };
const kpiPick = (kpis, id) => h('select.plain', { id, style: { maxWidth: '320px' } }, kpis.map((k) => h('option', { value: k.key }, `${k.label}${k.unit && k.unit !== '-' ? ` [${k.unit}]` : ''}`)));

export async function moduleView(root, id) {
  const mod = await findModule(id);
  if (!mod) { root.append(h('div.empty', h('h2', 'Module not found'), h('p', `No analysis module with id ${id} is installed.`), h('a.btn', { href: '#/' }, 'Home'))); return { title: 'Not found' }; }
  const project = getProject();
  const init = defaults(mod, busValues(), project.inputs[mod.id], confirmedFactors(mod.id));
  let lastResult = null, lastValues = null;

  const out = h('div', { 'aria-live': 'polite' }, h('div.empty', 'Results appear here. Press Run — or just edit a value.'));
  const status = h('small.muted');
  const live = h('input', { type: 'checkbox', checked: prefs.get('live', true) });
  live.addEventListener('change', () => prefs.set('live', live.checked));
  let running = false, queued = false, lastMs = 0;

  async function run() {
    if (running) { queued = true; return; }
    if (!form.isValid()) { status.textContent = 'Fix the highlighted inputs first.'; return; }
    running = true; status.textContent = 'Computing…';
    const values = structuredCloneSafe(form.state.values);
    const t0 = performance.now();
    try {
      const result = await task('run', mod.id, values, {}, (f, m) => { status.textContent = `Computing… ${Math.round(100 * f)} %${m ? ` — ${m}` : ''}`; });
      lastMs = performance.now() - t0; lastResult = result; lastValues = values;
      clear(out).append(assurance(mod, values), renderResult(mod, values, result, { elapsed: lastMs }), practiceFeed(mod));
      status.textContent = '';
      const s = summarise(result);
      setInputs(mod.id, stripGeometry(values));
      recordRun(mod.id, { status: s.status, maxUC: s.maxUC, governing: s.governing, kpis: (result.kpis || []).map((k) => ({ key: k.key, label: k.label, value: k.value, unit: k.unit })), quantities: result.quantities || null, advice: result.advice || [] });
      if (result.publish) publish(result.publish, mod.id);
    } catch (err) {
      clear(out).append(h('div.banner.fail', icon('warn'), h('div', h('b', 'This case could not be computed. '), err.message)));
      status.textContent = '';
    } finally {
      running = false;
      if (queued) { queued = false; run(); }
    }
  }
  const autoRun = debounce(() => { if (live.checked && lastMs < 1500) run(); }, 350);
  const form = buildForm(mod, init, { onChange: () => autoRun() });

  const analyse = () => h('div.split',
    h('div.sticky', form.el,
      h('div.runbar',
        h('button.btn.primary', { type: 'button', onclick: run }, icon('play'), 'Run analysis'),
        h('button.btn', { type: 'button', title: 'Restore every input to the module default and current project values', onclick: () => { const d = defaults(mod, busValues(), {}, confirmedFactors(mod.id)); Object.assign(form.state.values, d.values); form.state.linked = d.linked; form.refresh(); run(); } }, icon('refresh'), 'Reset'),
        h('label.check', { title: 'Recompute automatically while you edit (fast analyses only)' }, live, h('span', 'Live')),
        status)),
    h('div', out,
      h('div.row', { style: { marginTop: '14px' } },
        h('button.btn.small', { type: 'button', onclick: () => (lastResult ? download(`${mod.id}.json`, JSON.stringify(runRecord(mod, lastValues, lastResult), null, 1), 'application/json') : toast('Run the analysis first.')) }, icon('download'), 'Run record (JSON)'),
        h('button.btn.small', { type: 'button', onclick: () => window.print() }, 'Print / save as PDF'))));

  // ---------- Explore: sensitivity + sweeps ----------
  const explore = () => {
    const box = h('div');
    const kpis = lastResult?.kpis || [];
    if (!kpis.length) return h('div.empty', 'Run the analysis once, then explore how each input moves the result.');
    const nums = numericInputs(mod, form.state.values);
    const kSel = kpiPick([...kpis, { key: 'maxUC', label: 'Highest utilisation' }], 'ex-kpi');
    const pct = h('input.plain', { type: 'number', value: 10, min: 1, max: 50, style: { width: '80px' }, 'aria-label': 'Variation in percent' });
    const res1 = h('div');
    const pb = progressBar();
    const doSens = async () => {
      clear(res1).append(pb.el);
      try {
        const r = await task('sensitivity', mod.id, structuredCloneSafe(form.state.values), { pct: Number(pct.value) || 10 }, pb.set);
        const k = kSel.value, base = r.base[k];
        const items = r.rows.filter((x) => x.out.low && x.out.high).map((x) => ({ name: x.name, low: x.out.low[k], high: x.out.high[k] })).filter((x) => Number.isFinite(x.low) && Number.isFinite(x.high) && (x.low !== base || x.high !== base));
        clear(res1);
        if (!items.length) { res1.append(h('div.empty', 'No input changes this result within the chosen variation.')); return; }
        const lab = kSel.selectedOptions[0].textContent;
        res1.append(renderPlot({ type: 'tornado', title: `What drives ${lab} (inputs ±${r.pct} %)`, base, items, x: { label: lab } }));
        const ranked = items.map((x) => ({ ...x, swing: Math.abs(x.high - x.low) })).sort((a, b) => b.swing - a.swing);
        res1.append(h('p.muted', { style: { marginTop: '8px' } }, `Most influential: ${ranked.slice(0, 3).map((x) => x.name).join(', ')}. The longer the bar, the more that input deserves measurement or margin.`));
      } catch (e) { clear(res1).append(h('div.banner.fail', e.message)); }
    };
    const xSel = h('select.plain', { style: { maxWidth: '260px' }, 'aria-label': 'Input to sweep' }, nums.map((i) => h('option', { value: i.key }, i.label)));
    const from = h('input.plain', { type: 'number', style: { width: '110px' }, 'aria-label': 'From' }), to = h('input.plain', { type: 'number', style: { width: '110px' }, 'aria-label': 'To' });
    const setRange = () => { const inp = nums.find((i) => i.key === xSel.value); const v = form.state.values[inp.key]; from.value = Number(Math.max(inp.min ?? -Infinity, v * 0.5).toPrecision(4)); to.value = Number(Math.min(inp.max ?? Infinity, v * 1.5 || 1).toPrecision(4)); };
    xSel.addEventListener('change', setRange); if (nums.length) setRange();
    const res2 = h('div');
    const doSweep = async () => {
      const pb2 = progressBar(); clear(res2).append(pb2.el);
      try {
        const inp = nums.find((i) => i.key === xSel.value);
        const r = await task('sweep', mod.id, structuredCloneSafe(form.state.values), { key: inp.key, from: Number(from.value), to: Number(to.value), n: 25 }, pb2.set);
        const k = kSel.value, kp = [...kpis, { key: 'maxUC', label: 'Highest utilisation', limit: 1 }].find((q) => q.key === k);
        const ys = r.grid[0].map((g) => (g ? g[k] : NaN));
        const series = [{ name: kp.label, x: r.xs, y: ys }];
        if (kp.limit != null) series.push({ name: 'Limit', x: [r.xs[0], r.xs[r.xs.length - 1]], y: [kp.limit, kp.limit], dash: true });
        clear(res2).append(renderPlot({ type: 'line', title: `${kp.label} vs ${inp.label}`, x: { label: inp.label, unit: inp.unit }, y: { label: kp.label, unit: kp.unit }, series }));
      } catch (e) { clear(res2).append(h('div.banner.fail', e.message)); }
    };
    box.append(h('div.card', h('h3', 'Sensitivity'), h('p.muted', 'Each numeric input is moved up and down in turn while the others stay fixed.'),
      h('div.row', h('label', { for: 'ex-kpi' }, 'Result'), kSel, h('label', 'Variation ±'), pct, h('span', '%'), h('button.btn.primary', { type: 'button', onclick: doSens }, 'Rank the inputs')), res1),
    h('div.card', h('h3', 'Parameter sweep'), h('div.row', h('span', 'Vary'), xSel, h('span', 'from'), from, h('span', 'to'), to, h('button.btn.primary', { type: 'button', onclick: doSweep }, 'Plot')), res2));
    return box;
  };

  // ---------- Uncertainty ----------
  const uncertainty = () => {
    const nums = numericInputs(mod, form.state.values).filter((i) => i.role !== 'control' && i.role !== 'criteria');
    if (!lastResult) return h('div.empty', 'Run the analysis once, then propagate input uncertainty through it.');
    const rows = nums.map((i) => {
      const on = h('input', { type: 'checkbox', checked: !!i.uncertain, 'aria-label': `Treat ${i.label} as uncertain` });
      const d = h('select', { 'aria-label': 'Distribution' }, ['normal', 'lognormal', 'uniform', 'triangular'].map((x) => h('option', { value: x, selected: (i.uncertain?.dist || 'normal') === x }, x)));
      const cov = h('input', { type: 'text', inputmode: 'decimal', value: (i.uncertain?.cov ?? 0.1) * 100, 'aria-label': 'Coefficient of variation in percent', style: { width: '70px' } });
      return { i, on, d, cov };
    });
    const n = h('input.plain', { type: 'number', value: 400, min: 50, max: 5000, style: { width: '100px' }, 'aria-label': 'Number of samples' });
    const res = h('div');
    const go = async () => {
      const vars = {};
      for (const r of rows) if (r.on.checked) {
        const base = form.state.values[r.i.key], c = (Number(r.cov.value) || 10) / 100;
        vars[r.i.key] = r.d.value === 'uniform' || r.d.value === 'triangular' ? { dist: r.d.value, min: base * (1 - 1.732 * c), max: base * (1 + 1.732 * c), mode: base } : { dist: r.d.value, cov: c };
      }
      const pb = progressBar(); clear(res).append(pb.el);
      try {
        const r = await task('monteCarlo', mod.id, structuredCloneSafe(form.state.values), { n: Number(n.value) || 400, seed: 1, vars }, pb.set);
        clear(res);
        const kpis = [...(lastResult.kpis || []), { key: 'maxUC', label: 'Highest utilisation', unit: '-', limit: 1 }].filter((k) => r.stats[k.key]);
        res.append(h('div.banner', { class: r.pFail > 0.05 ? 'fail' : '' }, icon(r.pFail > 0.05 ? 'warn' : 'info'), h('div', h('b', `${fmt(100 * r.pFail, 3)} % of the ${r.n} sampled cases miss at least one criterion`), ` (± ${fmt(100 * r.pFailSE, 2)} % sampling error). `,
          'This is the spread caused by the input uncertainty you selected — not a calibrated failure probability.', r.failedRuns ? ` ${r.failedRuns} samples could not be computed.` : '')));
        res.append(table([{ key: 'label', label: 'Result', left: true }, { key: 'mean', label: 'Mean' }, { key: 'std', label: 'Std dev' }, { key: 'p10', label: '10th percentile' }, { key: 'p50', label: 'Median' }, { key: 'p90', label: '90th percentile' }, { key: 'unit', label: 'Unit', left: true }],
          kpis.map((k) => ({ label: k.label, unit: k.unit || '', ...r.stats[k.key] }))));
        res.append(h('p.muted', 'Percentiles are non-exceedance values: 90 % of sampled cases fall at or below the 90th percentile.'));
        const g = h('div.grid.c2', { style: { marginTop: '12px' } });
        for (const k of kpis.slice(0, 4)) {
          g.append(renderPlot({ type: 'hist', title: k.label, values: Array.from(r.samples.outputs[k.key]), x: { label: k.label, unit: k.unit }, markers: [{ value: r.stats[k.key].p50, label: 'median' }, ...(k.limit != null ? [{ value: k.limit, label: 'limit' }] : [])] }));
        }
        const first = kpis[0];
        if (first && r.influence[first.key]) g.append(renderPlot({ type: 'bar', title: `Which uncertainty drives ${first.label}`, categories: r.influence[first.key].map((x) => nums.find((i) => i.key === x.key)?.label || x.key), series: [{ name: 'Rank correlation', y: r.influence[first.key].map((x) => x.rho) }], y: { label: 'Spearman rank correlation', min: -1, max: 1 } }));
        res.append(g);
      } catch (e) { clear(res).append(h('div.banner.fail', e.message)); }
    };
    return h('div.card', h('h3', 'Uncertainty propagation (Latin hypercube sampling)'),
      h('p.muted', 'Tick the inputs that are uncertain and give each a coefficient of variation. Sampling is seeded, so a rerun reproduces the same numbers.'),
      h('div.tablewrap', { style: { maxHeight: '300px', overflowY: 'auto' } }, h('table.t', h('thead', h('tr', h('th.l', 'Uncertain'), h('th.l', 'Input'), h('th', 'Value'), h('th.l', 'Distribution'), h('th', 'CoV [%]'))),
        h('tbody', rows.map((r) => h('tr', h('td.l', r.on), h('td.l', r.i.label), h('td', fmt(form.state.values[r.i.key]), ' ', h('small.muted', r.i.unit || '')), h('td.l', r.d), h('td', r.cov)))))),
      h('div.row', { style: { marginTop: '10px' } }, h('label', 'Samples'), n, h('button.btn.primary', { type: 'button', onclick: go }, 'Run uncertainty study')), res);
  };

  // ---------- Convergence ----------
  const convergenceTab = () => {
    const specs = [['mesh', mod.mesh], ['meshTime', mod.meshTime]].filter(([, s]) => s);
    if (!specs.length) return h('div.card', h('h3', 'No discretisation in this module'), h('p.muted', 'This analysis evaluates closed-form or code equations, so there is no mesh, segment count or time step to converge. Its numerical accuracy is that of floating-point arithmetic; see the Method tab for the verification cases.'));
    const box = h('div');
    for (const [which, spec] of specs) {
      const lv = h('input.plain', { type: 'text', value: spec.levels.join(', '), 'aria-label': 'Refinement levels', style: { maxWidth: '260px' } });
      const tol = h('input.plain', { type: 'number', value: 2, min: 0.01, max: 20, step: 0.1, style: { width: '80px' }, 'aria-label': 'Tolerance in percent' });
      const res = h('div');
      const go = async () => {
        const pb = progressBar(); clear(res).append(pb.el);
        try {
          const levels = lv.value.split(/[,\s]+/).map(Number).filter((x) => Number.isFinite(x) && x > 0);
          const r = await task('convergence', mod.id, structuredCloneSafe(form.state.values), { which, levels, tol: (Number(tol.value) || 2) / 100 }, pb.set);
          clear(res);
          res.append(h('div.banner', { class: r.met ? '' : 'warn' }, icon(r.met ? 'check' : 'warn'), h('div',
            r.met ? h('b', `Use ${r.controlLabel} = ${r.selected}. `) : h('b', r.trustworthy ? 'The tolerance is not reached on this ladder — add finer levels. ' : `The sequence is ${r.behaviour}, so the extrapolation is not reliable — add levels or check the model. `),
            `Estimated discretisation error at that level ${fmt(100 * r.selectedErr, 2)} % against the Richardson-extrapolated value ${fmt(r.fExt)} ${r.unit || ''}. Observed order of accuracy ${Number.isFinite(r.observedOrder) ? fmt(r.observedOrder, 3) : '—'}; fine-grid GCI ${fmt(100 * r.gciFine, 2)} %; asymptotic-range ratio ${fmt(r.asymptotic, 3)} (1 is ideal).`)));
          res.append(h('div.grid.c2',
            renderPlot({ type: 'line', title: `${r.qoiLabel} vs ${r.controlLabel}`, x: { label: r.controlLabel, log: true }, y: { label: r.qoiLabel, unit: r.unit, zero: false }, series: [{ name: 'Computed', x: r.runs.map((q) => q.level), y: r.runs.map((q) => q.value), markers: true }, { name: 'Richardson extrapolation', x: [r.runs[0].level, r.runs[r.runs.length - 1].level], y: [r.fExt, r.fExt], dash: true }] }),
            renderPlot({ type: 'line', title: 'Estimated discretisation error', x: { label: r.controlLabel, log: true }, y: { label: 'Relative error', log: true }, series: [{ name: 'Error estimate', x: r.runs.filter((q) => q.errEst > 0).map((q) => q.level), y: r.runs.filter((q) => q.errEst > 0).map((q) => q.errEst), markers: true }, { name: 'Tolerance', x: [r.runs[0].level, r.runs[r.runs.length - 1].level], y: [r.tol, r.tol], dash: true }] })));
          res.append(table([{ key: 'level', label: r.controlLabel }, { key: 'value', label: r.qoiLabel, unit: r.unit, sig: 7 }, { key: 'err', label: 'Error estimate [%]' }, { key: 'ms', label: 'Run time [ms]' }, { key: 'pick', label: '', left: true }],
            r.runs.map((q) => ({ level: q.level, value: q.value, err: 100 * q.errEst, ms: q.ms, pick: q.level === r.selected ? '◀ selected' : '' }))));
          res.append(h('div.row', { style: { marginTop: '8px' } }, h('button.btn', { type: 'button', onclick: () => { form.state.values[r.key] = r.selected; form.refresh(); toast(`${r.controlLabel} set to ${r.selected}.`); run(); } }, `Apply ${r.controlLabel} = ${r.selected}`)));
        } catch (e) { clear(res).append(h('div.banner.fail', e.message)); }
      };
      box.append(h('div.card', h('h3', which === 'mesh' ? `Spatial / discretisation convergence — ${spec.label || spec.key}` : `Time-step convergence — ${spec.label || spec.key}`),
        h('p.muted', 'The analysis is repeated on a ladder of refinement levels. Observed order, Richardson extrapolation and the Grid Convergence Index (Roache; Celik et al. 2008) follow, and the coarsest level that meets your tolerance is selected so accuracy is bought at the lowest cost.'),
        h('div.row', h('label', 'Levels'), lv, h('label', 'Tolerance'), tol, h('span', '%'), h('button.btn.primary', { type: 'button', onclick: go }, 'Run convergence study')), res));
    }
    return box;
  };

  // ---------- Calibration & validation ----------
  const cvvTab = () => {
    if (!lastResult) return h('div.empty', 'Run the analysis once, then calibrate it against your measurements.');
    const nums = numericInputs(mod, form.state.values);
    const kpis = lastResult.kpis || [];
    const vary = h('select.plain', { style: { maxWidth: '260px' }, 'aria-label': 'Input that changes between measurements' }, h('option', { value: '' }, '(none — repeated measurements of one case)'), nums.map((i) => h('option', { value: i.key }, i.label)));
    const kSel = kpiPick(kpis, 'cv-kpi');
    const calTa = h('textarea.plain', { rows: 6, placeholder: 'input value <tab> measured result   (or one measured result per line)', 'aria-label': 'Calibration data' });
    const valTa = h('textarea.plain', { rows: 6, placeholder: 'independent data, same columns — never reuse the calibration rows', 'aria-label': 'Validation data' });
    const pRows = nums.filter((i) => i.role !== 'control').map((i) => {
      const v = form.state.values[i.key];
      return { i, on: h('input', { type: 'checkbox', 'aria-label': `Calibrate ${i.label}` }), lo: h('input', { type: 'text', value: Number((i.min != null ? Math.max(i.min, v * 0.5) : v * 0.5).toPrecision(4)), style: { width: '90px' } }), hi: h('input', { type: 'text', value: Number((i.max != null ? Math.min(i.max, v * 1.5) : v * 1.5).toPrecision(4)), style: { width: '90px' } }) };
    });
    const parse = (txt) => txt.trim().split(/\r?\n/).map((l) => l.trim().split(/[\t,; ]+/).map(Number)).filter((r) => r.every(Number.isFinite) && r.length).map((r) => (vary.value && r.length >= 2 ? { set: { [vary.value]: r[0] }, kpi: kSel.value, value: r[1] } : { set: {}, kpi: kSel.value, value: r[r.length - 1] }));
    const res = h('div'), res2 = h('div');
    const parity = (title, pred, obs) => { const lo = Math.min(...pred, ...obs), hi = Math.max(...pred, ...obs); return renderPlot({ type: 'scatter', title, x: { label: 'Measured', zero: false }, y: { label: 'Predicted', zero: false }, series: [{ name: 'Cases', x: obs, y: pred }, { name: 'Perfect agreement', x: [lo, hi], y: [lo, hi], line: true, dash: true, markers: false }] }); };
    const doCal = async () => {
      const pb = progressBar(); clear(res).append(pb.el);
      try {
        const params = pRows.filter((r) => r.on.checked).map((r) => ({ key: r.i.key, min: Number(r.lo.value), max: Number(r.hi.value) }));
        const r = await task('calibrate', mod.id, structuredCloneSafe(form.state.values), { params, data: parse(calTa.value) }, pb.set);
        clear(res);
        recordEvidence(mod.id, 'calibration', { n: r.observed.length, kpi: kSel.value, vary: vary.value || null, rmse: r.metricsAfter.rmse, bias: r.metricsAfter.bias, params: r.params.map((q) => ({ key: q.key, fitted: q.fitted })), dataHash: hashRows(parse(calTa.value)), warnings: r.warnings.length });
        for (const w of r.warnings) res.append(h('div.banner.warn', icon('warn'), w));
        res.append(table([{ key: 'label', label: 'Parameter', left: true }, { key: 'initial', label: 'Before' }, { key: 'fitted', label: 'Calibrated' }, { key: 'min', label: 'Lower bound' }, { key: 'max', label: 'Upper bound' }], r.params));
        res.append(h('p', `Root-mean-square error ${fmt(r.metricsBefore.rmse)} → `, h('b', fmt(r.metricsAfter.rmse)), `; bias ${fmt(r.metricsBefore.bias)} → `, h('b', fmt(r.metricsAfter.bias)), ` over ${r.observed.length} measurements.`));
        res.append(parity('Calibration fit', r.after, r.observed));
        res.append(h('div.row', { style: { marginTop: '8px' } }, h('button.btn', { type: 'button', onclick: () => { for (const p of r.params) form.state.values[p.key] = Number(p.fitted.toPrecision(6)); form.refresh(); run(); toast('Calibrated values applied to the inputs.'); } }, 'Apply calibrated values')));
      } catch (e) { clear(res).append(h('div.banner.fail', e.message)); }
    };
    const doVal = async () => {
      clear(res2).append(h('small.muted', 'Comparing…'));
      try {
        const vData = parse(valTa.value);
        // validation only means something on data the model has not been tuned to
        const cal = evidenceOf(mod.id).calibration;
        const calRows = new Set(parse(calTa.value).map((d) => JSON.stringify([d.set, d.value])));
        const shared = vData.filter((d) => calRows.has(JSON.stringify([d.set, d.value]))).length;
        if (shared > 0 || (cal && cal.dataHash === hashRows(vData))) { clear(res2).append(h('div.banner.fail', icon('warn'), h('div', h('b', `${shared || vData.length} of these rows are the calibration data. `), 'Agreement with data the parameters were fitted to proves nothing about prediction. Use measurements that were held back.'))); return; }
        const r = await task('validation', mod.id, structuredCloneSafe(form.state.values), { data: vData, acceptance: 0.1 });
        clear(res2);
        const xs = vary.value ? vData.map((d) => d.set[vary.value]) : [];
        recordEvidence(mod.id, 'validation', { n: vData.length, kpi: kSel.value, kpiLabel: kSel.selectedOptions[0].textContent, vary: vary.value || null, varyLabel: vary.value ? vary.selectedOptions[0].textContent : null,
          range: xs.length ? [Math.min(...xs), Math.max(...xs)] : null, measuredRange: [Math.min(...r.observed), Math.max(...r.observed)], rmse: r.metrics.rmse, bias: r.metrics.bias, r2: r.metrics.r2, within: r.within, peakBias: r.peakBias, acceptance: 0.1, pass: r.pass, dataHash: hashRows(vData) });
        res2.append(h('div.banner', { class: r.pass ? '' : 'warn' }, icon(r.pass ? 'check' : 'warn'), h('div', h('b', r.pass ? 'Predictions agree with the independent data within ±10 %. ' : 'Agreement is outside the ±10 % acceptance band. '),
          `${fmt(100 * r.within, 3)} % of cases within the band; RMSE ${fmt(r.metrics.rmse)}, mean bias ${fmt(r.metrics.bias)}, R² ${fmt(r.metrics.r2, 3)}. Mean error on the largest fifth of measured values ${fmt(100 * r.peakBias, 3)} % (negative = under-prediction of peaks, which matters more than the average). This outcome is now recorded in the project as evidence for this module, with the range it covers.`)));
        res2.append(parity('Validation against independent data', r.predicted, r.observed));
        res2.append(h('p.muted', 'This statement holds for these measurements and this range of conditions only.'));
      } catch (e) { clear(res2).append(h('div.banner.fail', e.message)); }
    };
    return h('div',
      h('div.card', h('h3', 'Calibration against measurements'), h('p.muted', 'Calibration adjusts uncertain parameters inside physically justified bounds. It must not be used to compensate for wrong geometry, boundary conditions or physics — and the data used here must be kept apart from validation data.'),
        h('div.row', h('label', 'Measured quantity'), kSel, h('label', 'Varies with'), vary),
        h('div.grid.c2', { style: { marginTop: '10px' } },
          h('div', h('h4', 'Measured data'), calTa),
          h('div', h('h4', 'Parameters to adjust'), h('div.tablewrap', { style: { maxHeight: '190px', overflowY: 'auto' } }, h('table.t', h('thead', h('tr', h('th.l', ''), h('th.l', 'Parameter'), h('th', 'Lower'), h('th', 'Upper'))), h('tbody', pRows.map((r) => h('tr', h('td.l', r.on), h('td.l', r.i.label, ' ', h('small.muted', r.i.unit || '')), h('td', r.lo), h('td', r.hi)))))))),
        h('div.row', { style: { marginTop: '10px' } }, h('button.btn.primary', { type: 'button', onclick: doCal }, 'Calibrate')), res),
      h('div.card', h('h3', 'Validation against independent data'), h('p.muted', 'Paste measurements that were not used for calibration. The comparison reports error, bias and how the largest values are predicted.'), valTa,
        h('div.row', { style: { marginTop: '10px' } }, h('button.btn.primary', { type: 'button', onclick: doVal }, 'Validate')), res2));
  };

  // ---------- Advisor ----------
  const advisorTab = () => {
    const box = h('div');
    if (!lastResult) return h('div.empty', 'Run the analysis once to get recommendations.');
    const s = summarise(lastResult);
    box.append(h('div.card', h('h3', 'Recommendations for this result'), adviceList(lastResult.advice) || h('p.muted', s.status === 'ok' ? 'Every criterion is met; no corrective action is needed for this case.' : 'See the design options below.')));
    const res = h('div.card', h('h3', 'Design options'), h('p.muted', 'Searching the design variables…'));
    box.append(res);
    task('advise', mod.id, structuredCloneSafe(form.state.values), {}).then((r) => {
      clear(res).append(h('h3', r.status === 'fail' ? 'Smallest single changes that satisfy every criterion' : 'Margin available on each design variable'));
      if (!r.designVars) { res.append(h('p.muted', 'This module declares no design variables.')); return; }
      if (!r.options.length) { res.append(h('p.muted', r.status === 'fail' ? 'No single design variable can bring this case within its criteria inside the allowed input ranges — a combination of changes or a different concept is needed.' : 'Every criterion keeps passing across the whole allowed range of each design variable.')); return; }
      res.append(h('p.muted', r.status === 'fail' ? 'Each row changes one variable only and leaves the rest as entered; ranked by the size of the change.' : 'Each row shows how far one variable can move before the first criterion is reached — an indication of where the design carries reserve.'));
      res.append(table([{ key: 'label', label: 'Design variable', left: true }, { key: 'from', label: 'Now' }, { key: 'to', label: r.status === 'fail' ? 'Required' : 'Limit' }, { key: 'unit', label: 'Unit', left: true }, { key: 'pct', label: 'Change [%]' }, { key: 'apply', label: '', left: true }],
        r.options.map((o) => ({ ...o, pct: 100 * o.change, apply: h('button.btn.small', { type: 'button', onclick: () => { form.state.values[o.key] = Number(o.to.toPrecision(5)); form.refresh(); run(); toast(`${o.label} set to ${fmt(o.to)} ${o.unit || ''}.`); } }, 'Apply') }))));
    }).catch((e) => { clear(res).append(h('div.banner.fail', e.message)); });
    const lit = h('div.card', h('h3', 'Latest published work on this topic'), h('p.muted', 'Loading from OpenAlex…'));
    box.append(lit);
    const query = mod.searchTerms || `${mod.title} offshore ${SHORT[mod.suite]}`.replace(/[&,]/g, ' ');
    literature(query, { n: 6 }).then((r) => {
      clear(lit).append(h('h3', 'Latest published work on this topic'), h('small.muted', `Live from OpenAlex · fetched ${new Date(r.at).toLocaleString()}${r.stale ? ' (shown from cache — offline)' : ''}`));
      if (!r.data.length) { lit.append(h('p.muted', 'No recent publications found for this query.')); return; }
      lit.append(h('ul.clean', r.data.map((w) => h('li', { style: { padding: '8px 0', borderBottom: '1px solid var(--line)' } }, h('a', { href: w.oa || w.doi || '#', target: '_blank', rel: 'noopener' }, w.title || 'Untitled'), h('br'), h('small.muted', `${w.authors || ''} · ${w.venue || ''} · ${w.date || w.year}`)))));
    }).catch((e) => { clear(lit).append(h('h3', 'Latest published work on this topic'), h('p.muted', `Not available right now (${e.message}).`)); });
    return box;
  };

  // ---------- Method ----------
  const methodTab = async () => {
    const idx = await suiteIndex();
    const suite = idx.find((x) => x.n === mod.suite);
    const box = h('div');
    box.append(h('div.card', h('h3', 'What this module computes'), h('p', mod.summary),
      h('h4', 'Fidelity'), h('p', h('span.chip.info', FIDELITY[mod.fidelity]?.[0] || mod.fidelity), ' ', FIDELITY[mod.fidelity]?.[1] || ''),
      h('h4', 'Assumptions and limits'), h('p', mod.limits || '—'),
      h('h4', 'Equations and models used'), h('div.row', (mod.models || []).map((m) => h('span.chip', m))),
      h('h4', { style: { marginTop: '12px' } }, 'Standards and references'), h('ul.prose', (mod.standards || []).map((x) => h('li', x))),
      h('h4', 'Specified analyses covered'), h('ul.ticks', (mod.covers || []).map((c) => h('li.on', suite?.analyses[c - 1] || `#${c}`)))));
    const roles = {};
    for (const i of mod.inputs || []) (roles[i.role || 'operation'] ||= []).push(i.label);
    box.append(h('div.card', h('h3', 'Conditions this module needs'), ['ic', 'bc', 'geometry', 'material', 'soil', 'environment', 'load', 'operation', 'criteria', 'control', 'cost'].filter((r) => roles[r]).map((r) => h('p', h('b', `${ROLE_LABEL[r]}: `), roles[r].join(', '), '.'))));
    const ver = h('div.card', h('h3', 'Verification cases'), h('p.muted', 'Running the built-in benchmark cases…'));
    box.append(ver);
    task('verify', mod.id, {}).then((cases) => {
      clear(ver).append(h('h3', 'Verification cases'), h('p.muted', 'Each case compares this module with an independent closed-form or published value. Passing them shows the equations are implemented and solved correctly; it does not by itself show that the model predicts a real structure — use the Calibrate & validate tab with your own measurements for that.'));
      for (const c of cases) {
        ver.append(h('div', { style: { margin: '10px 0' } }, h('b', c.name), ' ', h('span', { class: `chip ${c.status === 'pass' ? 'ok' : c.status === 'pending' ? 'warn' : 'fail'}` }, c.status === 'pass' ? '✓ passes' : c.status === 'pending' ? 'no independent benchmark yet' : '✗ fails'), c.source ? h('div', h('small.muted', c.source)) : null, c.error ? h('div', h('small', c.error)) : null,
          c.rows.length ? table([{ key: 'label', label: 'Quantity', left: true }, { key: 'expected', label: 'Reference', sig: 6 }, { key: 'got', label: 'Computed', sig: 6 }, { key: 'errPct', label: 'Difference' }, { key: 'tolPct', label: 'Tolerance' }], c.rows.map((r) => ({ ...r, errPct: r.absolute ? fmt(r.error) : `${fmt(100 * r.error, 2)} %`, tolPct: r.absolute ? fmt(r.tol) : `${fmt(100 * r.tol, 2)} %` }))) : null));
      }
    }).catch((e) => { clear(ver).append(h('h3', 'Verification cases'), h('div.banner.fail', e.message)); });
    return box;
  };

  const idx = await suiteIndex();
  const mods = await suiteModules(mod.suite);
  const pos = mods.findIndex((m) => m.id === mod.id);
  root.append(
    h('div.row', h('span.chip.info', `Suite ${mod.suite}`), h('span.chip', { title: FIDELITY[mod.fidelity]?.[1] }, FIDELITY[mod.fidelity]?.[0] || mod.fidelity), (mod.standards || []).slice(0, 2).map((s) => h('span.chip', s))),
    h('h1', { style: { marginTop: '8px' } }, mod.title), h('p.lede', mod.summary),
    tabs([
      { label: 'Analyse', render: analyse }, { label: 'Explore', render: explore }, { label: 'Uncertainty', render: uncertainty }, { label: 'Convergence', render: convergenceTab },
      { label: 'Calibrate & validate', render: cvvTab }, { label: 'Advisor', render: advisorTab }, { label: 'Method', render: methodTab },
    ]),
    pager(pos > 0 ? { href: `#/m/${mods[pos - 1].id}`, label: mods[pos - 1].title } : { href: `#/suite/${mod.suite}`, label: `Suite ${mod.suite} overview` },
      pos < mods.length - 1 ? { href: `#/m/${mods[pos + 1].id}`, label: mods[pos + 1].title } : (mod.suite < 17 ? { href: `#/suite/${mod.suite + 1}`, label: `Suite ${mod.suite + 1} — ${SHORT[mod.suite + 1]}` } : { href: '#/decision', label: 'Decision support' })));
  run();
  return { title: mod.title, crumbs: [{ href: `#/suite/${mod.suite}`, label: `${mod.suite}. ${SHORT[mod.suite]}` }, { href: `#/m/${mod.id}`, label: mod.title }], suite: idx.find((s) => s.n === mod.suite) };
}

// What stands behind this result: the status of the standards-derived factors it used and the evidence
// (verification, calibration, validation) held in the project for this module.
function assurance(mod, values) {
  const fr = factorReport(mod, values, getProject().factors || {});
  const ev = evidenceOf(mod.id);
  const bits = [];
  if (fr.total) bits.push(h('span', { class: `chip ${fr.unconfirmed ? 'warn' : 'ok'}` }, fr.unconfirmed ? `${fr.unconfirmed} of ${fr.total} code factors unconfirmed` : `all ${fr.total} code factors confirmed`), ' ', h('a', { href: `#/factors/m/${mod.id}` }, fr.unconfirmed ? 'Check them against your governing edition →' : 'Register'));
  const v = ev.validation;
  bits.push(' ', v ? h('span', { class: `chip ${v.pass ? 'ok' : 'warn'}`, title: `${v.n} independent measurements, ${new Date(v.at).toLocaleDateString()}` }, v.pass ? `validated for ${v.kpiLabel || v.kpi}${v.range ? ` over ${v.varyLabel} ${fmt(v.range[0])}–${fmt(v.range[1])}` : ''}` : 'validation attempted — outside acceptance')
    : h('span.chip', { title: 'No independent measurements have been compared with this module in this project. Use the Calibrate & validate tab.' }, 'not validated against measurements'));
  if (ev.calibration) bits.push(' ', h('span.chip.info', `calibrated on ${ev.calibration.n} measurements`));
  return h('div', { style: { margin: '0 0 8px', fontSize: '.82rem', lineHeight: 2 } }, bits);
}
function hashRows(rows) { let hsh = 2166136261; for (const ch of JSON.stringify(rows.map((d) => [d.set, d.value]))) { hsh ^= ch.charCodeAt(0); hsh = Math.imul(hsh, 16777619); } return (hsh >>> 0).toString(16); }

// Live check that the recommendations stay in step with current practice: the newest publications on the
// module's topic, fetched by the browser (cached for a day, shown from cache when offline).
function practiceFeed(mod) {
  const box = h('details.fold', h('summary', 'Current practice watch — newest publications on this topic'), h('div', h('small.muted', 'Open to load…')));
  let loaded = false;
  box.addEventListener('toggle', () => {
    if (!box.open || loaded) return;
    loaded = true;
    const body = box.lastChild;
    const query = mod.searchTerms || `${mod.title} offshore ${SHORT[mod.suite]}`.replace(/[&,]/g, ' ');
    literature(query, { n: 5 }).then((r) => {
      clear(body).append(h('small.muted', `OpenAlex · fetched ${new Date(r.at).toLocaleString()}${r.stale ? ' · cached (offline)' : ''} · compare with the “Basis” cited under each recommendation`),
        h('ul.clean', r.data.map((w) => h('li', { style: { padding: '6px 0', borderBottom: '1px solid var(--line)' } }, h('a', { href: w.oa || w.doi || '#', target: '_blank', rel: 'noopener' }, w.title || 'Untitled'), h('br'), h('small.muted', `${w.authors || ''} · ${w.venue || ''} · ${w.date || w.year}`)))));
    }).catch((e) => { clear(body).append(h('small.muted', `Not available right now (${e.message}).`)); loaded = false; });
  });
  return box;
}

// Typed arrays survive structured cloning; the copy keeps the form state from being mutated by a module.
function structuredCloneSafe(v) { try { return structuredClone(v); } catch { return JSON.parse(JSON.stringify(v)); } }
function stripGeometry(values) { const o = {}; for (const [k, v] of Object.entries(values)) if (!(v && typeof v === 'object' && v.kind)) o[k] = v; return o; }
