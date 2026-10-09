import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  goldenSection,
  goalSeek,
  nelderMead,
  levenbergMarquardt,
  bfgs,
  differentialEvolution,
  particleSwarm,
  geneticAlgorithm,
  simulatedAnnealing,
  nsga2,
  nonDominatedSort,
  linprog,
  milp,
} from '../../src/engine/math/optim.js';
import { makeRng } from '../../src/engine/math/rng.js';

const rosen = (x) => 100 * (x[1] - x[0] * x[0]) ** 2 + (1 - x[0]) ** 2;
const himmel = (x) => (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2;
// the four global minima of Himmelblau's function (f = 0)
const HIMMEL_MINIMA = [
  [3, 2],
  [-2.805118, 3.131312],
  [-3.77931, -3.283186],
  [3.584428, -1.848126],
];
const distToHimmelMin = (x) => Math.min(...HIMMEL_MINIMA.map((m) => Math.hypot(x[0] - m[0], x[1] - m[1])));
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);
const box = [
  [-5, 5],
  [-5, 5],
];

test('goldenSection: minimum of (x-2)^2 + 1 on [0, 5] is x = 2, f = 1', () => {
  const r = goldenSection((x) => (x - 2) ** 2 + 1, 0, 5, { tol: 1e-10 });
  close(r.x, 2, 1e-7);
  close(r.f, 1, 1e-12);
});

test('goalSeek: x^3 = 10 on [0, 5] gives 10^(1/3)', () => {
  close(goalSeek((x) => x ** 3, 10, 0, 5), Math.cbrt(10), 1e-10);
  assert.throws(() => goalSeek((x) => x * x, -1, 0, 5), RangeError);
});

test('nelderMead: Rosenbrock minimum (1, 1) from the classical start (-1.2, 1)', () => {
  const r = nelderMead(rosen, [-1.2, 1]);
  assert.ok(r.converged);
  close(r.x[0], 1, 1e-5);
  close(r.x[1], 1, 1e-5);
  assert.ok(r.f < 1e-10);
});

test('nelderMead: Himmelblau from (0, 0) reaches a global minimum (f = 0)', () => {
  const r = nelderMead(himmel, [0, 0]);
  assert.ok(r.f < 1e-10);
  assert.ok(distToHimmelMin(r.x) < 1e-4);
});

test('nelderMead: bounds are respected (Rosenbrock with x0 <= 0.5 has its constrained minimum on the bound)', () => {
  let violated = false;
  const f = (x) => {
    if (x[0] > 0.5 + 1e-15 || x[0] < -2 || x[1] < -2 || x[1] > 2) violated = true;
    return rosen(x);
  };
  const r = nelderMead(f, [-1, 1], { bounds: [[-2, 0.5], [-2, 2]] });
  assert.equal(violated, false);
  // on x0 = 0.5 the minimum over x1 is x1 = 0.25, f = 0.25
  close(r.x[0], 0.5, 1e-6);
  close(r.x[1], 0.25, 1e-4);
  close(r.f, 0.25, 1e-8);
});

test('bfgs: Rosenbrock with numerical and analytical gradient', () => {
  const r1 = bfgs(rosen, [-1.2, 1]);
  close(r1.x[0], 1, 1e-5);
  close(r1.x[1], 1, 1e-5);
  const grad = (x) => [-400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]), 200 * (x[1] - x[0] * x[0])];
  const r2 = bfgs(rosen, [-1.2, 1], { grad, tol: 1e-10 });
  assert.ok(r2.converged);
  close(r2.x[0], 1, 1e-8);
  close(r2.x[1], 1, 1e-8);
  assert.ok(r2.iterations < 100);
});

test('bfgs: Himmelblau from (1, 1) converges to the minimum (3, 2)', () => {
  const r = bfgs(himmel, [1, 1]);
  assert.ok(r.f < 1e-10);
  assert.ok(distToHimmelMin(r.x) < 1e-5);
});

test('bfgs: bounded Rosenbrock, constrained minimum (0.5, 0.25) with f = 0.25', () => {
  let violated = false;
  const f = (x) => {
    if (x[0] > 0.5 || x[0] < -2 || x[1] < -2 || x[1] > 2) violated = true;
    return rosen(x);
  };
  const r = bfgs(f, [-1, 1], { bounds: [[-2, 0.5], [-2, 2]] });
  assert.equal(violated, false);
  close(r.x[0], 0.5, 1e-9);
  close(r.x[1], 0.25, 1e-5);
  close(r.f, 0.25, 1e-9);
});

test('levenbergMarquardt: Rosenbrock as residuals [10(y - x^2), 1 - x] and Himmelblau as residuals', () => {
  const r = levenbergMarquardt((x) => [10 * (x[1] - x[0] * x[0]), 1 - x[0]], [-1.2, 1]);
  assert.ok(r.converged);
  close(r.x[0], 1, 1e-8);
  close(r.x[1], 1, 1e-8);
  assert.ok(r.cost < 1e-16);
  const h = levenbergMarquardt((x) => [x[0] * x[0] + x[1] - 11, x[0] + x[1] * x[1] - 7], [1, 1]);
  assert.ok(distToHimmelMin(h.x) < 1e-6);
});

test('levenbergMarquardt: recovers y = a exp(-b t) + c from noisy data, covariance consistent with the errors', () => {
  const rng = makeRng(42);
  const truth = [5, 1.3, 0.7];
  const sigma = 0.05;
  const t = [];
  const y = [];
  for (let i = 0; i < 200; i++) {
    t.push((i * 4) / 199);
    y.push(truth[0] * Math.exp(-truth[1] * t[i]) + truth[2] + sigma * rng.normal());
  }
  const res = (p) => t.map((ti, i) => p[0] * Math.exp(-p[1] * ti) + p[2] - y[i]);
  const r = levenbergMarquardt(res, [1, 0.5, 0]);
  assert.ok(r.converged);
  assert.equal(r.residuals.length, 200);
  for (let j = 0; j < 3; j++) {
    const se = Math.sqrt(r.covariance[j][j]);
    assert.ok(se > 0 && se < 0.05, `se[${j}] = ${se}`);
    assert.ok(Math.abs(r.x[j] - truth[j]) < 4 * se, `parameter ${j}: ${r.x[j]} vs ${truth[j]} (se ${se})`);
    close(r.x[j], truth[j], 0.05);
  }
  // residual variance estimates the noise variance
  close(Math.sqrt(r.sigma2), sigma, 0.008);
  // covariance is symmetric
  close(r.covariance[0][1], r.covariance[1][0], 1e-12);
  // an analytical Jacobian gives the same answer
  const jac = (p) => t.map((ti) => [Math.exp(-p[1] * ti), -p[0] * ti * Math.exp(-p[1] * ti), 1]);
  const r2 = levenbergMarquardt(res, [1, 0.5, 0], { jac });
  for (let j = 0; j < 3; j++) close(r2.x[j], r.x[j], 1e-6);
  close(r2.covariance[1][1], r.covariance[1][1], 1e-6);
});

test('levenbergMarquardt: bounds are respected', () => {
  // unconstrained optimum is a = 2; with a <= 1.5 the solution sits on the bound
  const r = levenbergMarquardt((p) => [p[0] - 2, 0.5 * (p[0] - 2)], [0], { bounds: [[-1, 1.5]] });
  close(r.x[0], 1.5, 1e-12);
});

test('differentialEvolution: Rosenbrock and Himmelblau global minima, deterministic for a seed', () => {
  const r = differentialEvolution(rosen, box, { seed: 3 });
  close(r.x[0], 1, 1e-4);
  close(r.x[1], 1, 1e-4);
  assert.ok(r.f < 1e-8);
  const r2 = differentialEvolution(rosen, box, { seed: 3 });
  assert.deepEqual(r2.x, r.x);
  assert.equal(r.history.length, r.generations);
  for (let i = 1; i < r.history.length; i++) assert.ok(r.history[i] <= r.history[i - 1]);
  const h = differentialEvolution(himmel, box, { seed: 5 });
  assert.ok(h.f < 1e-8);
  assert.ok(distToHimmelMin(h.x) < 1e-4);
});

test('particleSwarm: Rosenbrock and Himmelblau global minima, deterministic for a seed', () => {
  const r = particleSwarm(rosen, box, { seed: 7, maxIter: 400 });
  close(r.x[0], 1, 1e-3);
  close(r.x[1], 1, 1e-3);
  assert.ok(r.f < 1e-6);
  assert.deepEqual(particleSwarm(rosen, box, { seed: 7, maxIter: 400 }).x, r.x);
  const h = particleSwarm(himmel, box, { seed: 2 });
  assert.ok(h.f < 1e-6);
  assert.ok(distToHimmelMin(h.x) < 1e-3);
});

test('geneticAlgorithm: Rosenbrock and Himmelblau global minima, deterministic for a seed', () => {
  const r = geneticAlgorithm(rosen, box, { seed: 11, pop: 80, generations: 300 });
  assert.ok(r.f < 1e-3, `f = ${r.f}`);
  close(r.x[0], 1, 0.05);
  close(r.x[1], 1, 0.1);
  assert.deepEqual(geneticAlgorithm(rosen, box, { seed: 11, pop: 80, generations: 300 }).x, r.x);
  const h = geneticAlgorithm(himmel, box, { seed: 4 });
  assert.ok(h.f < 1e-3, `f = ${h.f}`);
  assert.ok(distToHimmelMin(h.x) < 0.02);
});

test('geneticAlgorithm: integer variables stay integral ((x-3.2)^2 + (y+1.6)^2 -> (3, -2))', () => {
  const r = geneticAlgorithm((x) => (x[0] - 3.2) ** 2 + (x[1] + 1.6) ** 2, [[-10, 10], [-10, 10]], { seed: 1, integer: [true, true], generations: 60 });
  assert.deepEqual(r.x, [3, -2]);
});

test('simulatedAnnealing: Rosenbrock and Himmelblau global minima, deterministic for a seed', () => {
  const r = simulatedAnnealing(rosen, [-1.2, 1], box, { seed: 9 });
  assert.ok(r.f < 1e-4, `f = ${r.f}`);
  close(r.x[0], 1, 0.02);
  close(r.x[1], 1, 0.04);
  assert.deepEqual(simulatedAnnealing(rosen, [-1.2, 1], box, { seed: 9 }).x, r.x);
  const h = simulatedAnnealing(himmel, [0, 0], box, { seed: 6 });
  assert.ok(h.f < 1e-4, `f = ${h.f}`);
  assert.ok(distToHimmelMin(h.x) < 0.01);
});

test('nonDominatedSort: hand-checked ranks', () => {
  // (1,5) (2,3) (4,1) are mutually non-dominated; (3,4) is dominated by (2,3); (5,5) by all of the first front
  const { rank, fronts } = nonDominatedSort([[1, 5], [2, 3], [4, 1], [3, 4], [5, 5]]);
  assert.deepEqual(Array.from(rank), [0, 0, 0, 1, 2]);
  assert.deepEqual(fronts.map((f) => f.slice().sort()), [[0, 1, 2], [3], [4]]);
});

test('nsga2: ZDT1 (Zitzler, Deb & Thiele 2000) front f2 = 1 - sqrt(f1)', () => {
  const nVar = 10;
  const zdt1 = (x) => {
    let s = 0;
    for (let i = 1; i < nVar; i++) s += x[i];
    const g = 1 + (9 * s) / (nVar - 1);
    return [x[0], g * (1 - Math.sqrt(x[0] / g))];
  };
  const bounds = Array.from({ length: nVar }, () => [0, 1]);
  const r = nsga2(zdt1, bounds, { pop: 100, generations: 200, seed: 1 });
  assert.ok(r.front.length >= 80, `front size ${r.front.length}`);
  // distance to the analytical front; the extreme point f1 -> 0 is excluded from the max because there
  // f2 = g - sqrt(f1 g) is non-dominated for any g (boundary solutions are kept by the crowding operator)
  let worst = 0;
  let avg = 0;
  for (const p of r.front) {
    const d = Math.abs(p.f[1] - (1 - Math.sqrt(p.f[0])));
    avg += d / r.front.length;
    if (p.f[0] > 0.01) worst = Math.max(worst, d);
  }
  assert.ok(avg < 0.005, `mean distance to the true front ${avg}`);
  assert.ok(worst < 0.02, `max distance to the true front ${worst}`);
  // spread: the front covers f1 from ~0 to ~1 without large gaps
  assert.ok(r.front[0].f[0] < 0.02);
  assert.ok(r.front[r.front.length - 1].f[0] > 0.95);
  let gap = 0;
  for (let i = 1; i < r.front.length; i++) gap = Math.max(gap, r.front[i].f[0] - r.front[i - 1].f[0]);
  assert.ok(gap < 0.08, `largest gap ${gap}`);
  // members are mutually non-dominated
  for (const p of r.front) for (const q of r.front) assert.ok(!(q.f[0] < p.f[0] && q.f[1] < p.f[1]));
  const again = nsga2(zdt1, bounds, { pop: 100, generations: 200, seed: 1 });
  assert.deepEqual(again.front[5], r.front[5]);
});

test('linprog: Wyndor Glass (Hillier & Lieberman, Intro. to OR, §3.1): max 3x + 5y -> (2, 6), Z = 36', () => {
  const r = linprog({ c: [-3, -5], Aub: [[1, 0], [0, 2], [3, 2]], bub: [4, 12, 18] });
  assert.equal(r.status, 'optimal');
  close(r.x[0], 2, 1e-9);
  close(r.x[1], 6, 1e-9);
  close(r.f, -36, 1e-9);
});

test('linprog: >= and = constraints (two-phase): min 2x + 3y, x + y >= 4, x + 3y >= 6, x - y = 1', () => {
  // x = y + 1 -> 2y + 1 >= 4 -> y >= 1.5 ; 4y + 1 >= 6 -> y >= 1.25 ; optimum y = 1.5, x = 2.5, f = 9.5
  const r = linprog({ c: [2, 3], Aub: [[-1, -1], [-1, -3]], bub: [-4, -6], Aeq: [[1, -1]], beq: [1] });
  assert.equal(r.status, 'optimal');
  close(r.x[0], 2.5, 1e-9);
  close(r.x[1], 1.5, 1e-9);
  close(r.f, 9.5, 1e-9);
});

test('linprog: bounds, free and upper-bounded variables', () => {
  // min x + y with -3 <= x <= 5, y free, y >= x - 10, y >= -x - 2 -> x = 4 (kink), y = -6, f = -2
  const r = linprog({ c: [1, 1], Aub: [[1, -1], [-1, -1]], bub: [10, 2], bounds: [[-3, 5], [-Infinity, Infinity]] });
  assert.equal(r.status, 'optimal');
  close(r.f, -2, 1e-9);
  // max x with x <= 7.5 only via the bound
  const r2 = linprog({ c: [-1], bounds: [[null, 7.5]] });
  close(r2.x[0], 7.5, 1e-12);
});

test('linprog: infeasible and unbounded problems are reported', () => {
  assert.equal(linprog({ c: [1, 1], Aub: [[1, 1], [-1, -1]], bub: [1, -3] }).status, 'infeasible');
  assert.equal(linprog({ c: [-1, 0], Aub: [[-1, 1]], bub: [1] }).status, 'unbounded');
});

test('linprog: classical degenerate cycling LP terminates with Bland rule at x = (1, 0, 1, 0), f = -5/4', () => {
  // min -3/4 x1 + 20 x2 - 1/2 x3 + 6 x4; Dantzig's rule cycles on this problem.
  // Hand check: x = (1, 0, 1, 0): row1 = 0.25 - 1 <= 0, row2 = 0.5 - 0.5 <= 0, x3 <= 1, f = -0.75 - 0.5 = -1.25
  const r = linprog({
    c: [-0.75, 20, -0.5, 6],
    Aub: [
      [0.25, -8, -1, 9],
      [0.5, -12, -0.5, 3],
      [0, 0, 1, 0],
    ],
    bub: [0, 0, 1],
  });
  assert.equal(r.status, 'optimal');
  close(r.f, -1.25, 1e-9);
});

test('linprog: transportation problem (2 sources, 3 sinks) matches the hand optimum', () => {
  // supplies 20, 30; demands 10, 25, 15; costs [[2,3,1],[5,4,8]]
  // optimum: x13 = 15, x11 = 5, x21 = 5, x22 = 25 -> 15 + 10 + 25 + 100 = 150
  const r = linprog({
    c: [2, 3, 1, 5, 4, 8],
    Aeq: [
      [1, 1, 1, 0, 0, 0],
      [0, 0, 0, 1, 1, 1],
      [1, 0, 0, 1, 0, 0],
      [0, 1, 0, 0, 1, 0],
      [0, 0, 1, 0, 0, 1],
    ],
    beq: [20, 30, 10, 25, 15],
  });
  assert.equal(r.status, 'optimal');
  close(r.f, 150, 1e-8);
});

test('milp: 0/1 knapsack (values 60, 100, 120; weights 10, 20, 30; capacity 50) -> 220', () => {
  const r = milp({ c: [-60, -100, -120], Aub: [[10, 20, 30]], bub: [50], bounds: [[0, 1], [0, 1], [0, 1]], integer: [true, true, true] });
  assert.equal(r.status, 'optimal');
  assert.deepEqual(r.x, [0, 1, 1]);
  close(r.f, -220, 1e-9);
  assert.ok(r.nodes >= 1);
});

test('milp: 12-item knapsack equals exhaustive enumeration', () => {
  const v = [23, 31, 29, 44, 53, 38, 63, 85, 89, 82, 12, 47];
  const w = [92, 57, 49, 68, 60, 43, 67, 84, 87, 72, 15, 39];
  const cap = 300;
  let best = 0;
  for (let mask = 0; mask < 1 << 12; mask++) {
    let sv = 0;
    let sw = 0;
    for (let i = 0; i < 12; i++) {
      if (mask & (1 << i)) {
        sv += v[i];
        sw += w[i];
      }
    }
    if (sw <= cap && sv > best) best = sv;
  }
  const r = milp({ c: v.map((q) => -q), Aub: [w], bub: [cap], bounds: v.map(() => [0, 1]), integer: v.map(() => true) });
  assert.equal(r.status, 'optimal');
  close(-r.f, best, 1e-9);
  let sw = 0;
  r.x.forEach((xi, i) => {
    assert.ok(xi === 0 || xi === 1);
    sw += xi * w[i];
  });
  assert.ok(sw <= cap);
});

test('milp: mixed problem, max x + y with 2x + 2y <= 9, x integer, y continuous <= 2.3 -> x = 2, y = 2.3', () => {
  const r = milp({ c: [-1, -1], Aub: [[2, 2]], bub: [9], bounds: [[0, 10], [0, 2.3]], integer: [true, false] });
  assert.equal(r.status, 'optimal');
  close(r.f, -4.5, 1e-9);
  assert.ok(Number.isInteger(r.x[0]));
  // integer infeasible: 0.2 <= x <= 0.8
  assert.equal(milp({ c: [1], Aub: [[1], [-1]], bub: [0.8, -0.2], integer: [true] }).status, 'infeasible');
});
