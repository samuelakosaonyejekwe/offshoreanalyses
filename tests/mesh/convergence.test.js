import test from 'node:test';
import assert from 'node:assert/strict';
import { classify, observedOrder, richardson, gci, asymptoticCheck, threeGrid, fitOrder, selectLevel, runStudy, errorNorms, orderFromErrors } from '../../src/mesh/convergence.js';
import { close, rel } from './helpers.js';

// Celik et al. (2008), J. Fluids Eng. 130, 078001, Table 1, column "φ = dimensionless reattachment length"
const CELIK = { f1: 6.063, f2: 5.972, f3: 5.863, r21: 1.5, r32: 1.333 };

test('Celik et al. (2008) Table 1: p = 1.53, φ_ext = 6.1685, e_a = 1.5 %, e_ext = 1.7 %, GCI_fine = 2.2 %', () => {
  const o = observedOrder(CELIK);
  assert.ok(o.converged);
  close(o.p, 1.53, 0.01);
  const g = gci({ f1: CELIK.f1, f2: CELIK.f2, r21: CELIK.r21, p: o.p });
  close(g.extrapolated, 6.1685, 5e-4);
  close(100 * g.ea, 1.5, 0.05);
  close(100 * g.eext, 1.7, 0.05);
  close(100 * g.gci, 2.2, 0.05);
  close(richardson({ f1: CELIK.f1, f2: CELIK.f2, r21: CELIK.r21, p: o.p }), g.extrapolated, 1e-12);
  const t = threeGrid(CELIK);
  assert.equal(t.type, 'monotone');
  close(t.p, o.p, 1e-12); close(t.gciFine, g.gci, 1e-12); close(t.extrapolated, g.extrapolated, 1e-12);
  assert.equal(t.Fs, 1.25);
  assert.ok(t.inAsymptoticRange); close(t.asymptoticRatio, 1, 0.05);
});

test('constant refinement ratio reduces to p = ln(ε32 / ε21) / ln r', () => {
  const o = observedOrder({ f1: 1.01, f2: 1.04, f3: 1.16, r21: 2, r32: 2 });
  close(o.p, 2, 1e-12); assert.equal(o.iterations, 0); assert.equal(o.s, 1);
});

test('manufactured data f = f0 + C h^p: order, extrapolated value and asymptotic ratio recovered', () => {
  for (const p of [1, 2, 3.5]) {
    for (const [h1, h2, h3] of [[1, 2, 4], [1, 1.5, 2], [0.1, 0.13, 0.26]]) {
      const f = (h) => 7.25 - 0.003 * h ** p;
      const t = threeGrid({ f1: f(h1), f2: f(h2), f3: f(h3), h: [h1, h2, h3] });
      close(t.p, p, 1e-7, `p for h = ${h1}, ${h2}, ${h3}`);
      close(t.extrapolated, 7.25, 1e-8);
      // for exact power-law data the ratio is f1 / f2 (the two GCIs are normalised by different solutions)
      close(t.asymptoticRatio, f(h1) / f(h2), 1e-6);
      assert.ok(t.inAsymptoticRange);
      // GCI = 1.25 × exact relative error of the fine grid
      rel(t.gciFine, 1.25 * 0.003 * h1 ** p / f(h1), 1e-5);
    }
  }
});

test('classification: monotone, oscillatory, divergent, converged', () => {
  assert.equal(classify(1.01, 1.04, 1.16).type, 'monotone');
  close(classify(1.01, 1.04, 1.16).R, 0.25, 1e-12);
  assert.equal(classify(1.01, 0.98, 1.06).type, 'oscillatory');
  assert.equal(classify(1.3, 1.1, 1.0).type, 'divergent');
  assert.equal(classify(0.7, 1.1, 1.0).type, 'oscillatory-divergent');
  assert.equal(classify(2, 2, 2.1).type, 'converged');
  assert.equal(classify(2, 2.1, 2.1).type, 'indeterminate');
  // unequal ratios: equal differences still converge when the finer step is the larger one (first order here)
  assert.equal(classify(1, 1.5, 2).type, 'divergent');
  assert.equal(classify(1, 1.5, 2, { r21: 1.5, r32: 4 / 3 }).type, 'monotone');
  assert.throws(() => classify(1, NaN, 2), RangeError);
});

test('oscillatory convergence still yields an order; divergence yields none unless a formal order is given', () => {
  // f = 1 + 0.1 (−1/2)^k on grids halved each time
  const t = threeGrid({ f1: 1 + 0.1 / 16, f2: 1 - 0.1 / 8, f3: 1 + 0.1 / 4, r21: 2, r32: 2 });
  assert.equal(t.type, 'oscillatory');
  close(t.p, 1, 1e-9);
  assert.match(t.notes.join(' '), /oscillate/);
  const d = threeGrid({ f1: 1.3, f2: 1.1, f3: 1.0, r21: 2, r32: 2 });
  assert.equal(d.type, 'divergent'); assert.ok(Number.isNaN(d.gciFine)); assert.ok(Number.isNaN(d.extrapolated));
  assert.match(d.notes.join(' '), /not converging/);
  const withFormal = threeGrid({ f1: 1.3, f2: 1.1, f3: 1.0, r21: 2, r32: 2 }, { formalOrder: 2 });
  assert.equal(withFormal.Fs, 3); assert.equal(withFormal.pUsed, 2);
  rel(withFormal.gciFine, 3 * (0.2 / 1.3) / 3, 1e-12);
  // super-convergence is capped at twice the formal order
  const fast = threeGrid({ f1: 1 + 1e-6, f2: 1 + 1e-6 * 2 ** 6, f3: 1 + 1e-6 * 4 ** 6, r21: 2, r32: 2 }, { formalOrder: 2 });
  close(fast.p, 6, 1e-6); assert.equal(fast.pUsed, 4);
  const same = threeGrid({ f1: 5, f2: 5, f3: 5.1, r21: 2, r32: 2 });
  assert.equal(same.type, 'converged'); assert.equal(same.gciFine, 0); assert.equal(same.extrapolated, 5);
});

test('GCI safety factors and argument checks', () => {
  const a = gci({ f1: 10, f2: 10.4, r21: 2, p: 2 });
  close(a.gci, 1.25 * 0.04 / 3, 1e-15);
  close(gci({ f1: 10, f2: 10.4, r21: 2, p: 2, Fs: 3 }).gci, 3 * 0.04 / 3, 1e-15);
  close(a.extrapolated, (4 * 10 - 10.4) / 3, 1e-12);
  const c = asymptoticCheck({ f1: 10 + 0.1, f2: 10 + 0.4, f3: 10 + 1.6, r21: 2, r32: 2, p: 2 });
  close(c.ratio, 1 / (10.4 / 10.1) * 1, 0.05); assert.ok(c.gci32 > c.gci21);
  assert.throws(() => gci({ f1: 1, f2: 2, r21: 1, p: 2 }), RangeError);
  assert.throws(() => gci({ f1: 1, f2: 2, r21: 2, p: 0 }), RangeError);
  assert.throws(() => observedOrder({ f1: 1, f2: 2, f3: 3, r21: 0.5, r32: 2 }), RangeError);
  assert.ok(Number.isNaN(observedOrder({ f1: 1, f2: 1, f3: 3, r21: 2, r32: 2 }).p));
});

test('least-squares fit over more than three grids recovers f0, C and p; noisy data give a residual', () => {
  const h = [0.5, 1, 2, 4, 8];
  const f = h.map((x) => 3 + 0.01 * x ** 1.7);
  const fit = fitOrder(h, f);
  close(fit.p, 1.7, 1e-6); close(fit.f0, 3, 1e-8); close(fit.C, 0.01, 1e-8);
  assert.ok(fit.residualStd < 1e-9); assert.equal(fit.atBound, false);
  const noisy = fitOrder(h, f.map((v, i) => v + (i % 2 ? 1 : -1) * 2e-3));
  assert.ok(noisy.residualStd > 1e-3 && Math.abs(noisy.p - 1.7) < 0.3);
  // weighted towards the fine grids
  close(fitOrder(h, f, { weights: h.map((x) => 1 / x) }).p, 1.7, 1e-6);
  assert.equal(fitOrder([1, 2, 4, 8], [1, 1.5, 1.2, 1.9]).atBound, true);
  assert.throws(() => fitOrder([1, 2], [1, 2]), RangeError);
});

test('coarsest admissible level and its cost saving', () => {
  const f = (h) => 100 * (1 + 0.02 * h * h);
  const levels = [8, 4, 2, 1].map((h) => ({ h, value: f(h) }));
  // errors with Fs = 1.25: 160 %, 40 %, 10 %, 2.5 %
  const s = selectLevel({ levels, extrapolated: 100, tolerance: 0.12, costExponent: 3 });
  assert.equal(s.index, 2); assert.ok(s.met);
  close(s.errors[2], 0.1, 1e-12);
  close(s.relativeCost, 1 / 8, 1e-12); close(s.saving, 7 / 8, 1e-12);
  const none = selectLevel({ levels, extrapolated: 100, tolerance: 0.001 });
  assert.equal(none.index, -1); assert.equal(none.met, false); assert.equal(none.level.h, 1);
  // measured cost overrides the estimate
  const timed = selectLevel({ levels: levels.map((l) => ({ ...l, cost: 10 / l.h })), extrapolated: 100, tolerance: 0.5 });
  assert.equal(timed.index, 1); close(timed.relativeCost, 0.25, 1e-12);
});

test('runStudy: trapezoidal rule on ∫₀^π sin x dx — second order, extrapolates to 2, picks the coarsest adequate level', async () => {
  const trap = (n) => { let s = 0; const h = Math.PI / n; for (let i = 0; i <= n; i++) s += (i === 0 || i === n ? 0.5 : 1) * Math.sin(i * h); return s * h; };
  const seen = [];
  const r = await runStudy({ levels: [10, 20, 40, 80], evaluate: (n) => ({ value: trap(n), cost: n }), tolerance: 0.005, onProgress: (f, msg) => seen.push([f, msg]) });
  assert.equal(r.type, 'monotone');
  close(r.p, 2, 0.01);
  close(r.extrapolated, 2, 1e-6);
  assert.ok(r.converged && r.inAsymptoticRange);
  assert.deepEqual(r.levels.map((l) => l.level), [80, 40, 20, 10]);
  rel(r.levels[0].h, 1 / 80, 1e-12);
  // exact errors π²/(12 n²)·… relative: 0.82 % at n = 10, 0.21 % at n = 20; with Fs = 1.25 → 1.03 %, 0.26 %
  assert.equal(r.selected.level, 20);
  rel(r.selected.error, 1.25 * (2 - trap(20)) / 2, 1e-3);
  close(r.selected.relativeCost, 0.25, 1e-12);
  close(r.fit.p, 2, 0.01); close(r.fit.f0, 2, 1e-5);
  assert.equal(seen.length, 5); assert.equal(seen[4][0], 1);
  assert.deepEqual(r.notes, []);
});

test('runStudy: element size as the control, async evaluation, failure modes', async () => {
  // Simpson's rule is fourth order in the step
  const simpson = (h) => { const n = Math.round(1 / h); let s = 0; for (let i = 0; i <= n; i++) s += (i === 0 || i === n ? 1 : i % 2 ? 4 : 2) * Math.exp(i * h); return s * h / 3; };
  const r = await runStudy({ levels: [0.25, 0.125, 0.0625], evaluate: async (h) => simpson(h), isSize: true, tolerance: 1e-4 });
  close(r.p, 4, 0.05);
  close(r.extrapolated, Math.E - 1, 1e-8);
  assert.equal(r.selected.level, 0.25);
  assert.equal(r.fit, null);
  // total cell counts in three dimensions: h ∝ N^(−1/3)
  const r3 = await runStudy({ levels: [1000, 8000, 64000], evaluate: (N) => 5 + 2 / N ** (2 / 3), dimension: 3 });
  close(r3.p, 2, 1e-6); close(r3.extrapolated, 5, 1e-8);
  const stuck = await runStudy({ levels: [10, 20, 40], evaluate: (n) => 1 + 1 / n, tolerance: 1e-6 });
  assert.equal(stuck.selected, null); assert.equal(stuck.converged, false);
  assert.match(stuck.notes.join(' '), /Refine further/);
  await assert.rejects(() => runStudy({ levels: [1, 2], evaluate: () => 1 }), RangeError);
  await assert.rejects(() => runStudy({ levels: [1, 2, 3], evaluate: () => NaN }), /finite value/);
  await assert.rejects(() => runStudy({ levels: [1, 2, 2], evaluate: (n) => n }), RangeError);
  const ac = new AbortController(); ac.abort();
  await assert.rejects(() => runStudy({ levels: [1, 2, 3], evaluate: (n) => n, signal: ac.signal }), /cancelled/);
});

test('manufactured-solution error norms and the order of a central difference', () => {
  const e = errorNorms([1, 2, 3, 5], [1, 2, 4, 3]);
  close(e.l1, 0.75, 1e-12); close(e.l2, Math.sqrt(5 / 4), 1e-12); close(e.linf, 2, 1e-12); close(e.relL2, Math.sqrt(5 / 30), 1e-12);
  close(errorNorms([1, 3], [0, 0], [3, 1]).l1, 1.5, 1e-12); // weighted
  assert.throws(() => errorNorms([1], [1, 2]), RangeError);
  // d/dx sin x by central differences on n intervals of [0, 1]
  const hs = [], errs = [];
  for (const n of [16, 32, 64, 128]) {
    const h = 1 / n;
    const num = [], ex = [];
    for (let i = 1; i < n; i++) { num.push((Math.sin((i + 1) * h) - Math.sin((i - 1) * h)) / (2 * h)); ex.push(Math.cos(i * h)); }
    hs.push(h); errs.push(errorNorms(num, ex).l2);
  }
  const o = orderFromErrors(hs, errs);
  close(o.p, 2, 5e-3);
  assert.equal(o.pairwise.length, 3);
  o.pairwise.forEach((p) => close(p, 2, 1e-2));
  assert.ok(Math.abs(o.pairwise[0] - 2) < Math.abs(o.pairwise[2] - 2), 'order approaches 2 as the grid is refined');
  rel(o.coefficient * hs[3] ** o.p, errs[3], 1e-2);
  assert.throws(() => orderFromErrors([1, 1], [1, 2]), RangeError);
});
