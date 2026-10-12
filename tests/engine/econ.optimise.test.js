import test from 'node:test';
import assert from 'node:assert/strict';
import * as O from '../../src/engine/econ/optimise.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);

test('constrained minimum on an active inequality: x² + y² with x + y ≥ 1 → (½, ½), f = ½, multiplier 1; feasible', () => {
  for (const method of ['de', 'nm']) {
    const r = O.constrainedMinimise({ objective: (x) => x[0] ** 2 + x[1] ** 2, bounds: [[-3, 3], [-3, 3]], ineq: [(x) => 1 - x[0] - x[1]], method });
    close(r.x[0], 0.5, 1e-4); close(r.x[1], 0.5, 1e-4); close(r.f, 0.5, 1e-5);
    assert.equal(r.feasible, true); assert.equal(r.converged, true);
    assert.ok(r.maxViolation <= 1e-6);
    close(r.multipliers.ineq[0], 1, 1e-3); // KKT: 2x = λ
  }
});

test('active constraint with a large multiplier: feasibility is reported consistently with the stopping rule', () => {
  // conductor sizing shape: cost(A) = s (2 + 0.003 A + 20/A) is decreasing below A = 81.6, so a thermal minimum
  // A ≥ 400 is active. With s = 1000 the multiplier is 1000 × (0.003 − 20/400²) × 400 = 1150: the multiplier
  // method needs more than 8 outer iterations to bring the residual under 1e-6.
  const aMin = 400;
  const prob = { objective: (x) => 1000 * (2 + 0.003 * x[0] + 20 / x[0]), bounds: [[20, 2000]], ineq: [(x) => (aMin - x[0]) / aMin], maxGen: 60 };
  for (const scheme of ['auglag', 'penalty']) {
    const r = O.constrainedMinimise({ ...prob, scheme, method: 'nm' });
    assert.equal(r.feasible, true, `${scheme}: violation ${r.maxViolation}`);
    assert.equal(r.converged, true);
    assert.ok(r.maxViolation <= 1e-6);
    close(r.x[0], 400, 1e-3);
    close(r.f, 1000 * (2 + 1.2 + 0.05), 1e-3);
  }
  close(O.constrainedMinimise({ ...prob, method: 'nm' }).multipliers.ineq[0], 1150, 5);
  // an iteration budget that is too small is reported, not hidden: not converged and not feasible
  const short = O.constrainedMinimise({ ...prob, method: 'nm', outerIter: 3 });
  assert.equal(short.converged, false); assert.equal(short.feasible, short.maxViolation <= 1e-6);
  assert.equal(short.history.length, 3);
  // a looser feasibility tolerance is honoured by both the loop and the flag
  const loose = O.constrainedMinimise({ ...prob, method: 'nm', feasTol: 1e-3 });
  assert.equal(loose.feasible, true); assert.ok(loose.history.length <= O.constrainedMinimise({ ...prob, method: 'nm' }).history.length);
  assert.throws(() => O.constrainedMinimise({ ...prob, feasTol: 0 }), RangeError);
  assert.throws(() => O.constrainedMinimise({ ...prob, outerIter: 0 }), RangeError);
});

test('equality constraint and integer design: min (x − 2)² + (y − 1)² on x = 2y − ½ within the unit disc; integer rounding', () => {
  // Lagrange on the line alone: y = 1.2, x = 1.9 lies outside x² + y² ≤ 1, so both constraints are active:
  // (2y − ½)² + y² = 1 → 5y² − 2y − ¾ = 0 → y = (2 + √19)/10 = 0.63589, x = 0.77178
  const r = O.constrainedMinimise({ objective: (x) => (x[0] - 2) ** 2 + (x[1] - 1) ** 2, bounds: [[-3, 3], [-3, 3]],
    ineq: [(x) => x[0] ** 2 + x[1] ** 2 - 1], eq: [(x) => x[0] - 2 * x[1] + 0.5], method: 'nm' });
  const y = (2 + Math.sqrt(19)) / 10;
  close(r.x[1], y, 1e-4); close(r.x[0], 2 * y - 0.5, 1e-4); assert.equal(r.feasible, true);
  assert.equal(r.eqValues.length, 1); assert.equal(r.ineqValues.length, 1);
  // unconstrained with an integer variable: (n − 2.6)² is smallest at n = 3
  const i = O.constrainedMinimise({ objective: (x) => (x[0] - 2.6) ** 2 + (x[1] - 0.3) ** 2, bounds: [[0, 10], [0, 1]], integer: [true, false], method: 'de', maxGen: 60 });
  assert.equal(i.x[0], 3); close(i.x[1], 0.3, 1e-3); assert.equal(i.feasible, true);
});

test('knee point, Pareto front of two quadratics', () => {
  // extremes (0,1) and (1,0); (0.2, 0.2) lies farthest below the chord
  const k = O.kneePoint([{ f: [0, 1] }, { f: [0.2, 0.2] }, { f: [0.6, 0.1] }, { f: [1, 0] }]);
  assert.equal(k.index, 1);
  // f1 = x², f2 = (x − 2)² on [−1, 3]: Pareto set x ∈ [0, 2], knee at x = 1 (f = (1, 1))
  const m = O.multiObjective({ objectives: (x) => [x[0] ** 2, (x[0] - 2) ** 2], bounds: [[-1, 3]], pop: 40, generations: 40 });
  assert.ok(m.front.length > 5);
  for (const p of m.front) assert.ok(p.x[0] > -0.05 && p.x[0] < 2.05, `non-dominated point at x = ${p.x[0]}`);
  close(m.knee.x[0], 1, 0.15);
  assert.equal(m.extremes.length, 2);
});

test('knapsack and capital rationing: the profitability-index ranking fails on a lumpy project', () => {
  // classic: values 60, 100, 120; weights 10, 20, 30; capacity 50 → items 2 and 3, value 220
  const k = O.knapsack({ values: [60, 100, 120], weights: [10, 20, 30], capacity: 50 });
  assert.deepEqual(k.selected.slice().sort(), [1, 2]); close(k.value, 220, 1e-12); close(k.weight, 50, 1e-12);
  close(O.knapsack({ values: [5, 4], weights: [2.5, 2.4], capacity: 4.9, resolution: 0.1 }).value, 9, 1e-12);
  // budget 100: A (cost 60, NPV 30, PI 0.5), B (50, 22, 0.44), C (50, 21, 0.42): ranking funds A only (30);
  // the optimum is B + C (43)
  const c = O.capitalRationing({ projects: [{ name: 'A', cost: 60, npv: 30 }, { name: 'B', cost: 50, npv: 22 }, { name: 'C', cost: 50, npv: 21 }], budget: 100 });
  assert.equal(c.ranking[0].name, 'A'); close(c.ranking[0].pi, 0.5, 1e-12);
  assert.deepEqual(c.greedy.names, ['A']); close(c.greedy.npv, 30, 1e-12);
  assert.deepEqual(c.optimal.names.slice().sort(), ['B', 'C']); close(c.optimal.npv, 43, 1e-12);
  assert.equal(c.rankingIsOptimal, false);
});

test('portfolio MILP with budget rows, exclusion and dependency', () => {
  // values 10, 8, 7, 4; cost 5, 4, 3, 2; budget 9 → {0, 1} = 18 (cost 9) beats {0, 2} = 17 and {1, 2, 3} = 19? cost 9 → 19
  const base = { values: [10, 8, 7, 4], budgets: [{ use: [5, 4, 3, 2], limit: 9 }] };
  const p = O.portfolioMilp(base);
  close(p.value, 19, 1e-9); assert.deepEqual(p.selected.slice().sort(), [1, 2, 3]); close(p.usage[0], 9, 1e-9);
  // projects 1 and 2 mutually exclusive: best is {0, 1} = 18
  const ex = O.portfolioMilp({ ...base, exclusive: [[1, 2]] });
  close(ex.value, 18, 1e-9); assert.deepEqual(ex.selected.slice().sort(), [0, 1]);
  // project 3 requires project 0: {1, 2, 3} is no longer allowed → {0, 1} = 18
  const rq = O.portfolioMilp({ ...base, requires: [[3, 0]] });
  close(rq.value, 18, 1e-9);
});

test('decision criteria on a payoff table: minimax regret, maximin, maximax, Laplace, Hurwicz, expected value', () => {
  // payoffs (rows: build small / medium / large; columns: low / high demand)
  const r = O.minMaxRegret({ alternatives: ['small', 'medium', 'large'], scenarios: ['low', 'high'], payoffs: [[50, 60], [20, 120], [-40, 200]], probabilities: [0.6, 0.4] });
  // regret: best per scenario 50 and 200 → small [0, 140], medium [30, 80], large [90, 0]
  assert.deepEqual(r.regret, [[0, 140], [30, 80], [90, 0]]); assert.deepEqual(r.maxRegret, [140, 80, 90]);
  assert.equal(r.minimaxRegret, 'medium'); assert.equal(r.maximin, 'small'); assert.equal(r.maximax, 'large');
  // Laplace means 55, 70, 80 → large; Hurwicz α = 0.5: 55, 70, 80 → large; EV: 54, 60, 56 → medium
  assert.equal(r.laplace, 'large'); assert.equal(r.hurwicz, 'large'); assert.equal(r.expectedValue, 'medium');
  // as costs the roles reverse
  const c = O.minMaxRegret({ alternatives: ['a', 'b'], scenarios: ['s1', 's2'], payoffs: [[10, 30], [18, 20]], sense: 'min' });
  assert.deepEqual(c.regret, [[0, 10], [8, 0]]); assert.equal(c.minimaxRegret, 'b');
});

test('chance constraint: smallest capacity x with P[demand ≤ x] ≥ 0.9 for normal demand is the 90 % quantile', () => {
  const r = O.chanceConstrained({ objective: (x) => x[0], constraint: (x, xi) => xi.d - x[0], vars: { d: { dist: 'normal', mean: 100, std: 10 } },
    bounds: [[50, 200]], p: 0.9, n: 2000, seed: 4, maxGen: 60 });
  close(r.x[0], 100 + 1.2816 * 10, 0.8); assert.ok(r.probability >= 0.9 - 1e-9 && r.probability < 0.92); assert.equal(r.feasible, true);
});

test('two-stage stochastic LP (newsvendor): order at the critical fractile; EVPI = RP − WS', () => {
  // buy x at 1 per unit now; demand 10, 20 or 30 with probability 1/3 each; unmet demand is bought at 3 per unit:
  // T x + W y ≥ d  ⇔  −x − y ≤ −d. Optimal x = quantile at 1 − c/q = 2/3 → x = 20; cost = 20 + (1/3) × 3 × 10 = 30.
  const sc = [10, 20, 30].map((d) => ({ p: 1 / 3, q: [3], T: [[-1]], W: [[-1]], h: [-d] }));
  const r = O.twoStageStochastic({ c: [1], scenarios: sc });
  close(r.x[0], 20, 1e-6); close(r.objective, 30, 1e-6); close(r.firstStageCost, 20, 1e-6); close(r.expectedRecourse, 10, 1e-6);
  // wait and see: buy exactly the demand → 20 on average; EVPI = 10
  close(r.waitAndSee, 20, 1e-6); close(r.evpi, 10, 1e-6);
});

test('dynamic programming: staged path, equipment replacement and maintenance policy', () => {
  // two stages; from "start": invest (−10 → "built") or wait (0 → "start"); "built" earns 8 per stage
  const d = O.dpStaged({ stages: 3, initial: 'start', actions: (stage, state) => (state === 'start'
    ? [{ name: 'invest', reward: -10, next: 'built' }, { name: 'wait', reward: 0, next: 'start' }] : [{ name: 'operate', reward: 8, next: 'built' }]) });
  close(d.value, 6, 1e-12); assert.deepEqual(d.path.map((p) => p.action), ['invest', 'operate', 'operate']);
  // with a heavy discount the investment no longer pays: −10 + 0.4 × 8 + 0.16 × 8 = −5.52 < 0
  close(O.dpStaged({ stages: 3, initial: 'start', discount: 0.4, actions: (stage, state) => (state === 'start'
    ? [{ name: 'invest', reward: -10, next: 'built' }, { name: 'wait', reward: 0, next: 'start' }] : [{ name: 'operate', reward: 8, next: 'built' }]) }).value, 0, 1e-12);
  // replacement: O&M grows by 30 per year of age, new unit 100, no salvage, 6 years, no discounting.
  // Keeping a unit n years costs 30 (0 + 1 + … + n−1) O&M; cycles of 3 years: 100 + 0 + 30 + 60 = 190 per cycle
  const e = O.equipmentReplacementDp({ horizon: 6, maxAge: 6, price: 100, om: (age) => 30 * age, salvage: () => 0, initialAge: 0 });
  // never replace: 30 (0 + … + 5) = 450; replace once after 3 years: 90 + 100 + 90 = 280; twice (after 2, 4): 30 + 100 + 30 + 100 + 30 = 290
  close(e.pvCost, 280, 1e-9); assert.deepEqual(e.replacements, [4]); assert.equal(e.plan.length, 6);
  // maintenance MDP: two states (good, failed); failure probability 0.5 per period when good; failed costs 100 per
  // period; repair costs 20 and returns the asset to good before the period's deterioration
  const m = O.maintenanceScheduleDp({ horizon: 3, P: [[0.5, 0.5], [0, 1]], stateCost: [0, 100],
    actions: [{ name: 'nothing', cost: [0, 0], next: [0, 1] }, { name: 'repair', cost: [20, 20], next: [0, 0] }] });
  assert.equal(m.policy[0][1], 'repair'); assert.equal(m.policy[0][0], 'nothing'); assert.equal(m.firstAction, 'nothing');
  assert.ok(m.value[0][1] > m.value[0][0] && m.expectedCost === m.value[0][0]);
});

test('surrogate optimisation, production mix and transportation plan', () => {
  const s = O.surrogateOptimise({ f: (x) => (x[0] - 1.3) ** 2 + 2, bounds: [[-2, 4]], nInit: 5, nIter: 12, seed: 2 });
  close(s.x[0], 1.3, 0.1); close(s.f, 2, 0.02); assert.equal(s.evaluations, 17);
  for (let i = 1; i < s.history.length; i++) assert.ok(s.history[i] <= s.history[i - 1] + 1e-12);
  // product mix (Hillier & Lieberman, Wyndor Glass): max 3x1 + 5x2, x1 ≤ 4, 2x2 ≤ 12, 3x1 + 2x2 ≤ 18 → (2, 6), Z = 36
  const p = O.productionPlan({ margin: [3, 5], use: [[1, 0], [0, 2], [3, 2]], capacity: [4, 12, 18] });
  close(p.x[0], 2, 1e-8); close(p.x[1], 6, 1e-8); close(p.profit, 36, 1e-8); close(p.slack[0], 2, 1e-8); close(p.slack[2], 0, 1e-8);
  close(O.productionPlan({ margin: [3, 5], use: [[1, 0], [0, 2], [3, 2]], capacity: [4, 12, 18], maxDemand: [10, 3] }).profit, 3 * 4 + 5 * 3, 1e-8);
  // transportation: supplies 20, 30; demands 10, 25, 15; costs [[2, 3, 1], [5, 4, 8]].
  // Source 1 serves sink 3 (15) and 5 of sink 1; source 2 serves the rest: 15 + 10 + 25 + 100 = 150
  const t = O.transportPlan({ supply: [20, 30], demand: [10, 25, 15], cost: [[2, 3, 1], [5, 4, 8]] });
  close(t.cost, 150, 1e-7);
  close(t.flow[0][2], 15, 1e-7); close(t.flow[1][1], 25, 1e-7);
  for (let j = 0; j < 3; j++) close(t.flow[0][j] + t.flow[1][j], [10, 25, 15][j], 1e-7);
});
