import test from 'node:test';
import assert from 'node:assert/strict';
import * as D from '../../src/engine/econ/decision.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const crit = [{ name: 'capacity', weight: 3, sense: 'max' }, { name: 'cost', weight: 1, sense: 'min' }];
const matrix = [[10, 100], [20, 400], [15, 200]]; // A, B, C

test('normalisation schemes by hand', () => {
  // min–max: capacity (x − 10)/10 → 0, 1, 0.5; cost (400 − x)/300 → 1, 0, 0.6667
  const mm = D.normaliseMatrix(matrix, crit);
  close(mm[0][0], 0, 1e-12); close(mm[1][0], 1, 1e-12); close(mm[2][0], 0.5, 1e-12);
  close(mm[0][1], 1, 1e-12); close(mm[1][1], 0, 1e-12); close(mm[2][1], 2 / 3, 1e-12);
  // 'max': x/max = 0.5, 1, 0.75; min/x = 1, 0.25, 0.5
  const mx = D.normaliseMatrix(matrix, crit, 'max');
  close(mx[0][0], 0.5, 1e-12); close(mx[2][0], 0.75, 1e-12); close(mx[1][1], 0.25, 1e-12);
  // 'sum': x/Σx = 10/45; (1/x)/Σ(1/x) with Σ = 0.0175 → 0.5714, 0.1429, 0.2857
  const sm = D.normaliseMatrix(matrix, crit, 'sum');
  close(sm[0][0], 10 / 45, 1e-12); close(sm[0][1], 0.01 / 0.0175, 1e-12);
  // 'vector': x/√(100 + 400 + 225) = 10/26.926
  close(D.normaliseMatrix(matrix, crit, 'vector')[0][0], 10 / Math.sqrt(725), 1e-12);
  assert.throws(() => D.normaliseMatrix(matrix, crit, 'rank'), RangeError);
});

test('weighted sum: weights 0.75 / 0.25 → scores A 0.25, B 0.75, C 0.5417; ranking B, C, A', () => {
  const w = D.weightedSum({ alternatives: ['A', 'B', 'C'], criteria: crit, matrix });
  close(w.weights[0], 0.75, 1e-12);
  close(w.scores[0], 0.25, 1e-12); close(w.scores[1], 0.75, 1e-12); close(w.scores[2], 0.375 + 0.25 * 2 / 3, 1e-12);
  assert.equal(w.best, 'B'); assert.deepEqual(w.ranking.map((r) => r.name), ['B', 'C', 'A']); assert.equal(w.ranking[0].rank, 1);
});

test('TOPSIS: closeness of the ideal alternative is 1, of the anti-ideal 0; hand value for the middle one', () => {
  // dominating alternative D = (20, 100), dominated E = (10, 400), middle C = (15, 200); equal weights
  const t = D.topsis({ alternatives: ['D', 'E', 'C'], criteria: [{ name: 'cap', weight: 1, sense: 'max' }, { name: 'cost', weight: 1, sense: 'min' }], matrix: [[20, 100], [10, 400], [15, 200]] });
  close(t.closeness[0], 1, 1e-12); close(t.closeness[1], 0, 1e-12); assert.equal(t.best, 'D');
  // vector norms √725 and √210000; weighted (0.5): C = (7.5/26.926, 100/458.26) = (0.27854, 0.21822)
  //   ideal (0.37139, 0.10911), anti-ideal (0.18570, 0.43644): d+ = 0.14325, d− = 0.23719 → C = 0.6235
  const n1 = Math.sqrt(725); const n2 = Math.sqrt(210000);
  const dp = Math.hypot(0.5 * (20 - 15) / n1, 0.5 * (200 - 100) / n2); const dm = Math.hypot(0.5 * (15 - 10) / n1, 0.5 * (400 - 200) / n2);
  close(t.closeness[2], dm / (dp + dm), 1e-12); close(t.closeness[2], 0.6235, 5e-4);
});

test('AHP: a perfectly consistent matrix returns its generating weights with CR = 0; Saaty example is consistent enough', () => {
  const w = [0.6, 0.3, 0.1];
  const A = w.map((a) => w.map((b) => a / b));
  const r = D.ahp(A);
  for (let i = 0; i < 3; i++) close(r.weights[i], w[i], 1e-9);
  close(r.lambdaMax, 3, 1e-9); close(r.ci, 0, 1e-9); close(r.cr, 0, 1e-8); assert.equal(r.consistent, true);
  // [[1, 3, 5], [1/3, 1, 3], [1/5, 1/3, 1]]: principal eigenvalue 3.0385, weights 0.637, 0.258, 0.105, CR = 0.033
  const s = D.ahp([[1, 3, 5], [1 / 3, 1, 3], [1 / 5, 1 / 3, 1]]);
  close(s.lambdaMax, 3.0385, 2e-4); close(s.weights[0], 0.637, 1e-3); close(s.weights[1], 0.2583, 1e-3); close(s.weights[2], 0.1047, 1e-3);
  close(s.cr, (3.0385 - 3) / 2 / 0.58, 5e-4); assert.equal(s.consistent, true);
  // a strongly intransitive matrix is flagged
  assert.equal(D.ahp([[1, 9, 1 / 9], [1 / 9, 1, 9], [9, 1 / 9, 1]]).consistent, false);
  assert.equal(D.SAATY_RANDOM_INDEX[3], 0.58); assert.equal(D.SAATY_RANDOM_INDEX[5], 1.12);
  assert.throws(() => D.ahp([[1, 2], [3, 1]]), RangeError);
});

test('ELECTRE I and Pareto filter', () => {
  // D dominates E on both criteria: concordance(D, E) = 1, discordance 0 → D outranks E; kernel = {D}
  const e = D.electre1({ alternatives: ['D', 'E'], criteria: [{ name: 'cap', weight: 1, sense: 'max' }, { name: 'cost', weight: 1, sense: 'min' }], matrix: [[20, 100], [10, 400]] });
  close(e.concordance[0][1], 1, 1e-12); close(e.concordance[1][0], 0, 1e-12); close(e.discordance[0][1], 0, 1e-12);
  assert.equal(e.outranks[0][1], true); assert.equal(e.outranks[1][0], false); assert.deepEqual(e.kernel, ['D']);
  // points (capacity max, cost min): (20, 400) and (10, 100) and (15, 200) are efficient; (12, 300) is dominated by (15, 200)
  const p = D.paretoFilter([[20, 400], [10, 100], [15, 200], [12, 300]], ['max', 'min']);
  assert.deepEqual(p.indices, [0, 1, 2]); assert.deepEqual(p.dominated, [3]);
});

test('weight stability: the capacity weight at which the winner changes from B to A', () => {
  // min–max scores: A = (0, 1), B = (1, 0), C = (0.5, 0.667). With capacity weight w: A = 1 − w, B = w,
  // C = 0.5w + 0.667(1 − w). B stays best while w ≥ 1 − w and w ≥ C → w ≥ 0.5 and w ≥ 0.5714: lower limit 4/7
  const s = D.weightStability({ alternatives: ['A', 'B', 'C'], criteria: crit, matrix });
  assert.equal(s.best, 'B');
  const cap = s.intervals.find((i) => i.criterion === 'capacity');
  close(cap.weight, 0.75, 1e-12); close(cap.lower, 4 / 7, 2e-3); close(cap.upper, 1, 1e-9); assert.equal(cap.robust, false);
  const cost = s.intervals.find((i) => i.criterion === 'cost');
  close(cost.upper, 3 / 7, 2e-3); close(cost.lower, 0, 1e-9);
});

test('decision tree rollback: expected value, policy, exponential utility and the value of flexibility', () => {
  // drill (cost 30): 40 % success worth 100, else 0 → EV = 10; walk away = 0
  const tree = { type: 'decision', name: 'drill?', options: [
    { name: 'drill', payoff: -30, node: { type: 'chance', name: 'result', branches: [{ name: 'success', p: 0.4, node: { type: 'terminal', value: 100 } }, { name: 'dry', p: 0.6, node: { type: 'terminal', value: 0 } }] } },
    { name: 'walk', node: { type: 'terminal', value: 0 } }] };
  const r = D.rollbackTree(tree);
  close(r.expectedValue, 10, 1e-12); close(r.certaintyEquivalent, 10, 1e-12); close(r.riskPremium, 0, 1e-12);
  assert.equal(r.policy[0].choice, 'drill');
  // risk tolerance R = 50: E[u] = 0.4 (1 − e^(−70/50)) + 0.6 (1 − e^(30/50)) = 0.30136 − 0.49327 < 0 → walking away (CE 0) is preferred
  const a = D.rollbackTree(tree, { riskTolerance: 50 });
  assert.equal(a.policy[0].choice, 'walk'); close(a.certaintyEquivalent, 0, 1e-12);
  // a very large tolerance approaches risk neutrality
  const n = D.rollbackTree(tree, { riskTolerance: 1e7 });
  assert.equal(n.policy[0].choice, 'drill'); close(n.certaintyEquivalent, 10, 1e-2);
  // forced drilling under R = 50: CE = −50 ln(1 − E[u])
  const forced = { type: 'chance', name: 'result', branches: [{ name: 'success', p: 0.4, payoff: -30, node: { type: 'terminal', value: 100 } }, { name: 'dry', p: 0.6, payoff: -30, node: { type: 'terminal', value: 0 } }] };
  const eu = 0.4 * (1 - Math.exp(-1.4)) + 0.6 * (1 - Math.exp(0.6));
  const f = D.rollbackTree(forced, { riskTolerance: 50 });
  close(f.certaintyEquivalent, -50 * Math.log(1 - eu), 1e-9); close(f.expectedValue, 10, 1e-12); close(f.riskPremium, 10 - f.certaintyEquivalent, 1e-9);
  // flexibility: appraise first (learn the outcome), then build only on success, against committing to build now
  const flexible = { type: 'chance', name: 'appraisal', branches: [
    { name: 'good', p: 0.4, node: { type: 'decision', name: 'build?', options: [{ name: 'build', payoff: -30, node: { type: 'terminal', value: 100 } }, { name: 'stop', node: { type: 'terminal', value: 0 } }] } },
    { name: 'bad', p: 0.6, node: { type: 'decision', name: 'build?', options: [{ name: 'build', payoff: -30, node: { type: 'terminal', value: 0 } }, { name: 'stop', node: { type: 'terminal', value: 0 } }] } }] };
  const v = D.valueOfFlexibility({ flexible, rigid: forced });
  close(v.flexible, 28, 1e-12); close(v.rigid, 10, 1e-12); close(v.valueOfFlexibility, 18, 1e-12);
  assert.throws(() => D.rollbackTree({ type: 'chance', name: 'x', branches: [{ name: 'a', p: 0.5, node: { type: 'terminal', value: 1 } }] }), RangeError);
});

test('comparison table and recommendation composer', () => {
  const t = D.comparisonTable({ alternatives: [{ name: 'base', values: { npv: 100, co2: 50 } }, { name: 'alt', values: { npv: 120, co2: 70 } }],
    measures: [{ key: 'npv', label: 'NPV', unit: 'M', sense: 'max' }, { key: 'co2', label: 'Carbon', unit: 'kt', sense: 'min' }] });
  assert.equal(t.rows.length, 2); assert.equal(t.wins.base, 1); assert.equal(t.wins.alt, 1);
  assert.ok(Array.isArray(t.columns) && typeof t.title === 'string');
  const rec = D.composeRecommendations({
    checks: [{ name: 'Member UC', uc: 1.25, ref: 'ISO 19902' }, { name: 'Joint UC', uc: 0.3 }, { name: 'Pile UC', uc: 0.97 }],
    economics: { npv: -5e6, irr: 0.06, hurdleRate: 0.08, currency: 'USD' } });
  assert.ok(rec.length >= 3);
  assert.equal(rec[0].severity, 'fail'); assert.match(rec[0].text, /Member UC/); assert.match(rec[0].text, /1\.25/);
  for (let i = 1; i < rec.length; i++) assert.ok(rec[i].priority <= rec[i - 1].priority);
  assert.deepEqual(rec.map((r) => r.rank), rec.map((_, i) => i + 1));
  assert.ok(rec.some((r) => r.severity === 'warn' && /Pile UC/.test(r.text)));
  assert.ok(rec.some((r) => r.severity === 'fail' && /Net present value|Internal rate of return/.test(r.text)));
  assert.deepEqual(D.composeRecommendations({ checks: [{ name: 'x', uc: 1.2 }] }), D.composeRecommendations({ checks: [{ name: 'x', uc: 1.2 }] }));
});
