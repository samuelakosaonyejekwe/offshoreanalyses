// Contract and verification gate for every analysis module of every suite.
// Run all: node --test tests/suites.test.js      One suite: SUITE=8 node --test tests/suites.test.js
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { existsSync } from 'node:fs';
import { defaults, summarise, kpiOk, verify, validate, convergence } from '../src/core/studies.js';
import { material } from '../src/library/materials.js';

const pad = (n) => String(n).padStart(2, '0');
const only = process.env.SUITE ? process.env.SUITE.split(',').map(Number) : null;
const ctx = { site: {}, bus: {}, material, progress() {}, signal: new AbortController().signal };
const ROLES = ['geometry', 'material', 'load', 'environment', 'soil', 'ic', 'bc', 'operation', 'control', 'criteria', 'cost'];
const TYPES = ['number', 'int', 'select', 'bool', 'text', 'table', 'series', 'geometry'];
const PLOTS = ['line', 'scatter', 'bar', 'area', 'hist', 'polar', 'heatmap', 'contour', 'tornado', 'model3d'];
const seen = new Set();

for (let n = 1; n <= 17; n++) {
  if (only && !only.includes(n)) continue;
  const file = new URL(`../src/suites/s${pad(n)}/index.js`, import.meta.url);
  if (!existsSync(file)) { test(`suite ${n} is installed`, { todo: 'suite not built yet' }, () => {}); continue; }
  const suite = (await import(file)).default;
  const spec = JSON.parse(await readFile(new URL(`../src/spec/suite-${pad(n)}.json`, import.meta.url), 'utf8'));

  test(`suite ${n}: covers every specified analysis`, () => {
    assert.equal(suite.n, n);
    const covered = new Set(suite.modules.flatMap((m) => m.covers || []));
    const missing = spec.analyses.map((a, i) => [i + 1, a]).filter(([i]) => !covered.has(i));
    assert.deepEqual(missing, [], `analyses without a module: ${missing.map(([i, a]) => `${i}. ${a}`).join('; ')}`);
  });

  for (const mod of suite.modules) {
    test(`${mod.id}: contract`, () => {
      assert.match(mod.id, new RegExp(`^s${pad(n)}-[a-z0-9-]+$`));
      assert.ok(!seen.has(mod.id), 'duplicate id'); seen.add(mod.id);
      assert.equal(mod.suite, n);
      for (const k of ['title', 'summary', 'limits']) assert.ok(typeof mod[k] === 'string' && mod[k].length > 3, `${k} missing`);
      assert.ok(['analytical', 'numerical', 'reduced-order', 'screening'].includes(mod.fidelity), 'fidelity');
      assert.ok(Array.isArray(mod.standards) && mod.standards.length, 'standards');
      assert.ok(Array.isArray(mod.models) && mod.models.length, 'models');
      assert.ok(Array.isArray(mod.covers) && mod.covers.length && mod.covers.every((c) => Number.isInteger(c) && c >= 1 && c <= spec.analyses.length), 'covers out of range');
      assert.ok(Array.isArray(mod.inputs) && mod.inputs.length, 'inputs');
      const keys = new Set();
      for (const inp of mod.inputs) {
        assert.ok(inp.key && !keys.has(inp.key), `duplicate or missing input key ${inp.key}`); keys.add(inp.key);
        assert.ok(inp.label, `${inp.key}: label`);
        assert.ok(TYPES.includes(inp.type || 'number'), `${inp.key}: type ${inp.type}`);
        assert.ok(ROLES.includes(inp.role), `${inp.key}: role '${inp.role}'`);
        assert.ok(inp.group, `${inp.key}: group`);
        if ((inp.type || 'number') === 'number' || inp.type === 'int') assert.ok(Number.isFinite(inp.default), `${inp.key}: numeric default`);
        if (inp.type === 'select' && !inp.source) assert.ok(inp.options?.some((o) => o.value === inp.default), `${inp.key}: default not among options`);
        if (inp.type === 'select' && inp.source) assert.ok(material(inp.default), `${inp.key}: default material '${inp.default}' not in library`);
        if (inp.type === 'table') assert.ok(Array.isArray(inp.columns) && Array.isArray(inp.default), `${inp.key}: table columns/default`);
      }
      assert.ok(mod.inputs.some((i) => i.role === 'bc' || i.role === 'ic') || n === 17, 'no input is tagged as an initial or boundary condition');
      assert.deepEqual(validate(mod, defaults(mod).values), {}, 'defaults violate their own min/max');
      for (const m of [mod.mesh, mod.meshTime].filter(Boolean)) { assert.ok(keys.has(m.key), `mesh key ${m.key}`); assert.ok(m.levels?.length >= 3, 'mesh levels'); assert.ok(m.qoi, 'mesh qoi'); }
      assert.ok(Array.isArray(mod.verify) && mod.verify.length, 'verify cases');
    });

    test(`${mod.id}: runs with defaults`, async () => {
      const { values } = defaults(mod);
      // Processor time, not wall-clock time: the limit must mean "this module is too heavy", never
      // "the machine happened to be busy" (a wall-clock limit failed s01-pushover at 8.1 s under load).
      const c0 = process.cpuUsage();
      const res = await mod.compute(structuredClone(values), ctx);
      const c = process.cpuUsage(c0);
      const ms = Math.round((c.user + c.system) / 1000);
      assert.ok(ms < 8000, `default run used ${ms} ms of processor time`);
      assert.ok(res && Array.isArray(res.kpis) && res.kpis.length, 'no kpis');
      const kk = new Set();
      for (const k of res.kpis) { assert.ok(k.key && !kk.has(k.key), `kpi key ${k.key}`); kk.add(k.key); assert.ok(k.label, 'kpi label'); assert.ok(Number.isFinite(k.value), `kpi ${k.key} = ${k.value}`); if (k.limit != null) assert.ok(['min', 'max'].includes(k.sense), `kpi ${k.key}: sense`); assert.ok(typeof k.explain === 'string' && k.explain.length >= 25, `kpi ${k.key}: explain missing`); }
      for (const c of res.checks || []) { assert.ok(c.name, 'check name'); assert.ok(Number.isFinite(c.uc) && c.uc >= 0, `check '${c.name}' uc = ${c.uc}`); assert.ok(typeof c.explain === 'string' && c.explain.length >= 25, `check '${c.name}': explain missing`); }
      for (const p of res.plots || []) {
        assert.ok(PLOTS.includes(p.type), `plot type ${p.type}`); assert.ok(p.title, 'plot title'); assert.ok(typeof p.explain === 'string' && p.explain.length >= 25, `plot '${p.title}': explain missing`);
        for (const s of p.series || []) { const len = (s.x || s.thetaDeg || s.y || []).length; assert.ok(len <= 4000, `series '${s.name}' has ${len} points`); }
      }
      for (const t of res.tables || []) { assert.ok(Array.isArray(t.columns) && Array.isArray(t.rows), 'table shape'); assert.ok(t.rows.length <= 2000, 'table too long'); assert.ok(typeof t.explain === 'string' && t.explain.length >= 25, `table '${t.title}': explain missing`); }
      assert.ok(Array.isArray(res.advice) && res.advice.length, 'advice must not be empty');
      for (const a of res.advice) { assert.ok(['info', 'warn', 'fail', 'ok'].includes(a.severity), 'advice severity'); assert.ok(a.text && a.action && a.basis, `advice needs text, action and basis: ${a.text}`); assert.ok(['engineering', 'economic'].includes(a.kind), `advice kind: ${a.text}`); }
      for (const kind of ['engineering', 'economic']) assert.ok(res.advice.some((a) => a.kind === kind), `no ${kind}-practice recommendation`);
      assert.ok(res.advice.some((a) => typeof a.sustainability === 'string' && a.sustainability.length >= 25), 'no advice item carries a sustainability consequence');
      for (const [k, v] of Object.entries(res.publish || {})) { assert.ok(Number.isFinite(v), `publish ${k}`); assert.ok((mod.publishes || []).includes(k), `publishes[] lacks ${k}`); }
      for (const m of [mod.mesh, mod.meshTime].filter(Boolean)) assert.ok(kk.has(m.qoi), `mesh qoi '${m.qoi}' is not a kpi`);
      assert.doesNotThrow(() => structuredClone(res), 'result must be structured-cloneable (no functions)');
      // The built-in example is the first thing a user sees: it must satisfy its own acceptance criteria.
      const sum = summarise(res);
      const bad = [...res.kpis.filter((k) => kpiOk(k) === false).map((k) => `kpi ${k.key} = ${k.value} (limit ${k.limit}, ${k.sense})`), ...(res.checks || []).filter((c) => c.uc > 1).map((c) => `check '${c.name}' uc = ${c.uc.toFixed(2)}`)];
      assert.equal(sum.status, 'ok', `the default case fails its own criteria (${bad.join('; ')}): supply a default design that passes, or report a screening result as a KPI without limit plus advice of severity 'warn'`);
      // same inputs → same outputs
      const again = await mod.compute(structuredClone(values), ctx);
      for (const k of res.kpis) assert.equal(again.kpis.find((q) => q.key === k.key).value, k.value, `kpi ${k.key} not reproducible`);
    });

    test(`${mod.id}: verification cases`, async () => {
      const cases = await verify(mod, ctx);
      const bad = cases.filter((c) => c.status === 'fail');
      assert.deepEqual(bad.map((c) => `${c.name}: ${c.error || c.rows.filter((r) => !r.pass).map((r) => `${r.kpi} got ${r.got}, expected ${r.expected}`).join(', ')}`), []);
    });

    if (mod.mesh) test(`${mod.id}: discretisation study runs`, async () => {
      const r = await convergence(mod, defaults(mod).values, ctx, { levels: mod.mesh.levels.slice(0, 3) });
      assert.ok(Number.isFinite(r.fExt), 'extrapolated value');
    });
  }
}
