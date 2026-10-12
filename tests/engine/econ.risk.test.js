import test from 'node:test';
import assert from 'node:assert/strict';
import * as R from '../../src/engine/econ/risk.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('quantile (type 7), normal pdf / cdf / inverse and Spearman correlation', () => {
  // PERCENTILE.INC of 1…5: q = 0.25 → 2; q = 0.9 → 4.6
  close(R.quantile([1, 2, 3, 4, 5], 0.25), 2, 1e-12); close(R.quantile([1, 2, 3, 4, 5], 0.9), 4.6, 1e-12);
  close(R.quantile([10, 20], 0.5), 15, 1e-12); close(R.quantile([7], 0.3), 7, 1e-12);
  close(R.normPdf(0), 0.3989422804014327, 1e-15); close(R.normPdf(1), 0.24197072451914337, 1e-15);
  close(R.normCdf(0), 0.5, 1e-15); close(R.normCdf(1.96), 0.9750021048517795, 1e-13); close(R.normCdf(-3), 0.0013498980316301, 1e-14);
  close(R.normCdf(-8), 6.220960574271785e-16, 1e-26); close(R.normCdf(1) + R.normCdf(-1), 1, 1e-15);
  close(R.normInv(0.975), 1.959963984540054, 1e-12); close(R.normInv(0.5), 0, 1e-14); close(R.normInv(1e-6), -4.753424308822899, 1e-10);
  for (const p of [0.001, 0.2, 0.7, 0.9999]) close(R.normCdf(R.normInv(p)), p, 1e-13);
  assert.throws(() => R.normInv(1), RangeError);
  // monotone transformation gives ρ = 1; reversed −1; textbook ranks with one swap: ρ = 1 − 6Σd²/(n(n² − 1)) = 1 − 12/120 = 0.9
  close(R.spearman([1, 2, 3, 4, 5], [1, 8, 27, 64, 125]), 1, 1e-12); close(R.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1, 1e-12);
  close(R.spearman([1, 2, 3, 4, 5], [1, 2, 3, 5, 4]), 0.9, 1e-12);
});

test('summary statistics state both percentile conventions: exceedance P90 = cumulative p10', () => {
  const x = Array.from({ length: 101 }, (_, i) => i);
  const s = R.summarise(x, { kind: 'resource' });
  close(s.mean, 50, 1e-12); close(s.std, Math.sqrt(101 * 102 / 12), 1e-9); close(s.standardError, s.std / Math.sqrt(101), 1e-12);
  close(s.cumulative.p10, 10, 1e-12); close(s.cumulative.p90, 90, 1e-12); close(s.exceedance.P90, 10, 1e-12); close(s.exceedance.P10, 90, 1e-12);
  close(s.headline.P90, 10, 1e-12); // resource convention: P90 is the low case
  close(R.summarise(x, { kind: 'cost' }).headline.P90, 90, 1e-12); // cost convention: P90 is the high case
  assert.equal(s.min, 0); assert.equal(s.max, 100); assert.equal(s.n, 101);
});

test('deterministic sensitivity: tornado, spider, two-way table and scenarios on NPV = price × volume − cost', () => {
  const model = (x) => x.price * x.volume - x.cost;
  const base = { price: 10, volume: 100, cost: 600 };
  const t = R.oneAtATime({ model, base, ranges: { price: { low: 8, high: 13 }, volume: { low: 90, high: 105 }, cost: { low: 500, high: 800 } } });
  close(t.base, 400, 1e-12);
  assert.deepEqual(t.items.map((i) => i.name), ['price', 'cost', 'volume']); // swings 500, 300, 150
  close(t.items[0].low, 200, 1e-12); close(t.items[0].high, 700, 1e-12); close(t.items[0].swing, 500, 1e-12);
  close(t.items[1].low, 500, 1e-12); close(t.items[1].high, 200, 1e-12);
  const sp = R.spider({ model, base, names: ['price', 'cost'], factors: [0.9, 1, 1.1] });
  assert.deepEqual(sp.series[0].y, [300, 400, 500]); assert.deepEqual(sp.series[1].y, [460, 400, 340]);
  const d = R.dataTable2({ model, base, x: { name: 'price', values: [8, 10] }, y: { name: 'volume', values: [50, 100, 150] } });
  assert.deepEqual(d.z, [[-200, -100], [200, 400], [600, 900]]);
  const sc = R.scenarios({ model, base, scenarios: { low: { price: 8 }, high: { price: 12, volume: 110 } } });
  assert.deepEqual(sc.map((s) => [s.name, s.result]), [['base', 400], ['low', 200], ['high', 720]]);
});

test('sampling: Latin hypercube reproduces the marginals; Iman–Conover imposes a rank correlation without changing them', () => {
  const s = R.sampleInputs({ vars: { a: { dist: 'normal', mean: 10, std: 2 }, b: { dist: 'uniform', min: 0, max: 1 } }, n: 2000, seed: 3 });
  assert.deepEqual(s.names, ['a', 'b']); assert.equal(s.samples.a.length, 2000);
  const st = R.summarise(s.samples.a);
  close(st.mean, 10, 0.01); rel(st.std, 2, 0.01); // stratified sampling: far tighter than plain Monte Carlo
  close(R.summarise(s.samples.b).mean, 0.5, 1e-3);
  assert.ok(Math.abs(R.spearman(s.samples.a, s.samples.b)) < 0.08);
  const c = R.sampleInputs({ vars: { a: { dist: 'normal', mean: 10, std: 2 }, b: { dist: 'uniform', min: 0, max: 1 } }, n: 2000, seed: 3, rankCorr: [{ a: 'a', b: 'b', rho: 0.7 }] });
  close(R.spearman(c.samples.a, c.samples.b), 0.7, 0.03);
  // the marginals are exactly the same numbers, only paired differently
  assert.deepEqual(Array.from(c.samples.a).sort((x, y) => x - y), Array.from(s.samples.a).sort((x, y) => x - y));
  assert.deepEqual(Array.from(c.samples.b).sort((x, y) => x - y), Array.from(s.samples.b).sort((x, y) => x - y));
  assert.deepEqual(R.sampleInputs({ vars: { a: { dist: 'normal', mean: 0, std: 1 } }, n: 50, seed: 9 }).samples, R.sampleInputs({ vars: { a: { dist: 'normal', mean: 0, std: 1 } }, n: 50, seed: 9 }).samples);
});

test('propagation: sum of independent normals has mean Σμ and variance Σσ²; correlation adds 2ρσ1σ2', () => {
  const vars = { a: { dist: 'normal', mean: 100, std: 10 }, b: { dist: 'normal', mean: 50, std: 20 } };
  const m = R.simulate({ model: (x) => x.a + x.b + x.k, vars, fixed: { k: 5 }, n: 4000, seed: 2 });
  close(m.stats.value.mean, 155, 0.05); rel(m.stats.value.std, Math.sqrt(500), 0.03);
  const c = R.simulate({ model: (x) => ({ sum: x.a + x.b, diff: x.a - x.b }), vars, n: 4000, seed: 2, rankCorr: [{ a: 'a', b: 'b', rho: 0.8 }] });
  // rank correlation 0.8 between normals ↔ Pearson 2 sin(π × 0.8/6) = 0.8135: Var = 500 + 2 × 0.8135 × 200 = 825.4
  rel(c.stats.sum.std, Math.sqrt(500 + 2 * 0.8135 * 200), 0.03); rel(c.stats.diff.std, Math.sqrt(500 - 2 * 0.8135 * 200), 0.05);
  // probability of a negative outcome and VaR against the normal closed form: X ~ N(20, 10)
  const x = R.simulate({ model: (v) => v.x, vars: { x: { dist: 'normal', mean: 20, std: 10 } }, n: 20000, seed: 5 }).outputs.value;
  const pb = R.probabilityBelow(x);
  close(pb.probability, R.normCdf(-2), 0.002); close(pb.standardError, Math.sqrt(pb.probability * (1 - pb.probability) / 20000), 1e-12);
  const nv = R.normalVaR({ mean: 20, std: 10, alpha: 0.99 });
  close(nv.z, 2.3263479, 1e-6); close(nv.var, -20 + 23.263479, 1e-5); close(nv.cvar, -20 + 10 * R.normPdf(nv.z) / 0.01, 1e-9);
  const v = R.valueAtRisk(x, 0.99);
  close(v.var, nv.var, 0.3); close(v.cvar, nv.cvar, 0.4); assert.ok(v.cvar >= v.var);
  close(R.valueAtRisk(x, 0.99, { relativeTo: 'mean' }).var, 23.26, 0.3);
});

test('expected monetary value, value of perfect and of imperfect information (drilling example)', () => {
  const e = R.expectedMonetaryValue([{ p: 0.4, value: 70 }, { p: 0.6, value: -30 }]);
  close(e.emv, 10, 1e-12); close(e.std, Math.sqrt(0.4 * 3600 + 0.6 * 1600), 1e-12);
  assert.throws(() => R.expectedMonetaryValue([{ p: 0.5, value: 1 }, { p: 0.3, value: 2 }]), RangeError);
  // actions [drill, walk]; states oil (0.4): [70, 0], dry (0.6): [−30, 0]
  const states = [{ p: 0.4, payoffs: [70, 0] }, { p: 0.6, payoffs: [-30, 0] }];
  const p = R.evpi(states);
  close(p.evWithoutInformation, 10, 1e-12); close(p.evWithPerfectInformation, 28, 1e-12); close(p.evpi, 18, 1e-12); assert.equal(p.bestAction, 0);
  // seismic: P(positive | oil) = 0.8, P(positive | dry) = 0.3
  //   P(+) = 0.32 + 0.18 = 0.5; P(oil | +) = 0.64 → drill: 0.64 × 70 − 0.36 × 30 = 34
  //   P(−) = 0.5; P(oil | −) = 0.16 → drill: 11.2 − 25.2 = −14 → walk (0). EV with information = 17; EVII = 7
  const i = R.evii({ states, likelihood: [[0.8, 0.2], [0.3, 0.7]] });
  close(i.signals[0].probability, 0.5, 1e-12); close(i.signals[0].posterior[0], 0.64, 1e-12); close(i.signals[0].value, 34, 1e-12);
  assert.equal(i.signals[1].bestAction, 1); close(i.evWithInformation, 17, 1e-12); close(i.evii, 7, 1e-12); close(i.evpi, 18, 1e-12);
  // an uninformative survey is worth nothing; a perfect one is worth EVPI
  close(R.evii({ states, likelihood: [[0.5, 0.5], [0.5, 0.5]] }).evii, 0, 1e-12);
  close(R.evii({ states, likelihood: [[1, 0], [0, 1]] }).evii, 18, 1e-12);
});

test('convergence diagnostics and importance measures for a linear model y = 3a + b', () => {
  const vars = { a: { dist: 'normal', mean: 0, std: 1 }, b: { dist: 'normal', mean: 0, std: 1 }, c: { dist: 'normal', mean: 0, std: 1 } };
  const sim = R.simulate({ model: (x) => 3 * x.a + x.b, vars, n: 4000, seed: 7 });
  const cv = R.convergence(sim.outputs.value, 10);
  const last = cv.n.length - 1;
  assert.ok(cv.n.length >= 10 && cv.n.length <= 11); assert.equal(cv.n[0], 2); assert.equal(cv.n[last], 4000);
  close(cv.runningMean[last], sim.stats.value.mean, 1e-12); close(cv.standardError[last], sim.stats.value.standardError, 1e-9);
  assert.ok(cv.standardError[last] < cv.standardError[1]);
  close(cv.relativeError, cv.standardError[last] / Math.abs(cv.runningMean[last]), 1e-12);
  // variance shares 9/10 and 1/10; the unused input gets none
  const vd = R.varianceDecomposition(sim.inputs, sim.outputs.value);
  close(vd.varianceShare.a, 0.9, 0.02); close(vd.varianceShare.b, 0.1, 0.02); close(vd.varianceShare.c, 0, 0.005); close(vd.r2, 1, 1e-9);
  close(vd.src.a, 3 / Math.sqrt(10), 0.02);
  const wm = R.whatMatters(sim.inputs, sim.outputs.value);
  assert.deepEqual(wm.map((w) => w.name), ['a', 'b', 'c']); close(wm.reduce((s, w) => s + w.share, 0), 1, 1e-12);
  assert.ok(wm[0].rankCorrelation > 0.9 && Math.abs(wm[2].rankCorrelation) < 0.06);
  // Sobol indices of the additive model: first = total = variance shares
  const sb = R.sobol({ model: (x) => 3 * x.a + x.b, vars, n: 1024, seed: 3 });
  close(sb.first.a, 0.9, 0.05); close(sb.total.b, 0.1, 0.05); close(sb.total.c, 0, 0.02); assert.equal(sb.ranking[0], 'a');
});

test('risk register: expected value Σ p × impact; a certain risk gives its fixed impact at every percentile', () => {
  const r = R.riskRegister({ risks: [{ name: 'weather', probability: 0.3, impact: 100 }, { name: 'soil', probability: 0.1, impact: { dist: 'triangular', min: 200, mode: 300, max: 700 } }], n: 20000, seed: 4 });
  close(r.expectedValue, 30 + 0.1 * 400, 1e-9);
  rel(r.mean, 70, 0.05);
  close(r.byRisk[0].expected + r.byRisk[1].expected, 70, 1e-9); close(r.byRisk.reduce((s, b) => s + b.share, 0), 1, 1e-12);
  // P(no risk occurs) = 0.7 × 0.9 = 0.63 > 0.5 → median 0; P80: weather alone occurs (0.27) → 100
  close(r.p50, 0, 1e-9); close(r.p80, 100, 1e-9); assert.ok(r.p90 >= 100);
  const c = R.riskRegister({ risks: [{ name: 'sure', probability: 1, impact: 55 }], n: 200 });
  close(c.p50, 55, 1e-12); close(c.p90, 55, 1e-12); close(c.std, 0, 1e-12);
});

test('cash-flow simulation: a price multiplier scales revenue; degenerate inputs reproduce the deterministic case', () => {
  const inputs = { timeline: { constructionYears: 1, operatingYears: 5 }, discountRate: 0.08, capex: [{ name: 'plant', amount: 1000 }],
    opex: { fixed: 50 }, production: 100, price: { type: 'fixed', price: 4 } };
  // deterministic: −1000 + (400 − 50)(P/A, 8 %, 5)/… with first operating year at t = 1 → NPV = −1000 + 350 × 3.99271 = 397.45
  const fixed = R.simulateCashflow({ inputs, vars: { priceFactor: { dist: 'uniform', min: 1, max: 1 + 1e-12 } }, n: 50, seed: 1 });
  close(fixed.npv.mean, -1000 + 350 * (1 - 1.08 ** -5) / 0.08, 1e-6); close(fixed.npv.std, 0, 1e-6);
  close(fixed.capex.mean, 1000, 1e-9); close(fixed.probNpvNegative, 0, 1e-12);
  // LCOE-type levelised cost = (1000 + 50 (P/A)) / (100 (P/A)) = 3.0046
  close(fixed.lcoe.mean, (1000 + 50 * 3.99271) / (100 * 3.99271), 1e-4);
  // price factor N(1, 0.2): NPV is linear in it with slope 400 (P/A) = 1597.08 → std 319.4
  const sim = R.simulateCashflow({ inputs, vars: { priceFactor: { dist: 'normal', mean: 1, std: 0.2 }, capexFactor: { dist: 'uniform', min: 0.9, max: 1.1 } }, n: 2000, seed: 2 });
  rel(sim.npv.std, Math.hypot(1597.08 * 0.2, 1000 * 0.2 / Math.sqrt(12)), 0.03);
  close(sim.npv.mean, 397.45, 5);
  assert.equal(sim.ranking[0].name, 'priceFactor');
  // P[NPV < 0]: NPV ≈ N(397.45, 324.6) → Φ(−1.2245) = 0.110
  close(sim.probNpvNegative, 0.110, 0.02);
  assert.ok(sim.var.cvar >= sim.var.var); assert.equal(sim.n, 2000);
  assert.ok(typeof sim.conventions === 'object');
});
