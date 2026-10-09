// Techno-economic optimisation wrappers over the shared optimisation library: constrained single-objective design
// optimisation, multi-objective trade-offs with knee point, portfolio selection (knapsack by DP and MILP), capital
// rationing, robust (min–max regret) choice, chance-constrained optimisation, two-stage stochastic programming,
// dynamic programming for replacement and maintenance scheduling, surrogate-assisted optimisation and linear
// production / logistics planning. Every routine returns its convergence history where one exists.

import { nelderMead, differentialEvolution, particleSwarm, simulatedAnnealing, geneticAlgorithm, nsga2, linprog, milp }
  from '../math/optim.js';
import { bayesOptimise } from '../reliability/surrogate.js';
import { requireNumber, requireSeries } from './tvm.js';
import { sampleInputs } from './risk.js';

function checkBounds(bounds) {
  if (!Array.isArray(bounds) || bounds.length === 0) throw new RangeError('bounds must be a non-empty array of [lo, hi]');
  for (const b of bounds) if (!b || !(b[1] > b[0]) || !Number.isFinite(b[0]) || !Number.isFinite(b[1])) throw new RangeError('each bound must be finite with lo < hi');
}

/**
 * Constrained single-objective design optimisation: minimise f(x) subject to g_i(x) <= 0, h_j(x) = 0 and bounds,
 * by an augmented-Lagrangian outer loop (or a fixed quadratic penalty) around a derivative-free solver. The first
 * outer iteration uses the chosen global method (differential evolution, particle swarm, simulated annealing,
 * genetic algorithm) or Nelder–Mead; later iterations are warm-started with Nelder–Mead. Mixed-integer designs
 * are handled by `integer` flags (variables rounded inside the objective, i.e. a derivative-free MINLP heuristic).
 * Ref: Hestenes (1969) and Powell (1969) (method of multipliers); Nocedal & Wright, Numerical Optimization, 2nd
 * ed., ch. 17; Storn & Price (1997) for DE.
 * @param {{objective: (x: number[]) => number, bounds: [number, number][], ineq?: ((x: number[]) => number)[],
 *   eq?: ((x: number[]) => number)[], method?: 'de' | 'pso' | 'sa' | 'ga' | 'nm', scheme?: 'auglag' | 'penalty',
 *   integer?: boolean[], x0?: number[], seed?: number, outerIter?: number, penalty?: number, tol?: number,
 *   maxGen?: number}} o constraints in the objective's natural units (scale them to order one for best behaviour)
 * @returns {{x: number[], f: number, feasible: boolean, maxViolation: number, ineqValues: number[],
 *   eqValues: number[], multipliers: {ineq: number[], eq: number[]}, history: {iter: number, f: number,
 *   maxViolation: number, penalty: number}[], innerHistory: number[]}}
 */
export function constrainedMinimise({ objective, bounds, ineq = [], eq = [], method = 'de', scheme = 'auglag', integer, x0, seed = 1,
  outerIter = 8, penalty = 10, tol = 1e-6, maxGen = 150 }) {
  if (typeof objective !== 'function') throw new RangeError('objective must be a function');
  checkBounds(bounds);
  const n = bounds.length;
  const snap = (x) => (integer ? x.map((v, j) => (integer[j] ? Math.round(v) : v)) : x);
  let lam = new Array(ineq.length).fill(0); let nu = new Array(eq.length).fill(0); let mu = penalty;
  const merit = (xr) => {
    const x = snap(xr);
    let L = objective(x);
    for (let i = 0; i < ineq.length; i++) {
      const g = ineq[i](x);
      if (scheme === 'penalty') L += 0.5 * mu * Math.max(0, g) ** 2;
      else { const s = Math.max(0, g + lam[i] / mu); L += 0.5 * mu * (s * s - (lam[i] / mu) ** 2); }
    }
    for (let j = 0; j < eq.length; j++) { const h = eq[j](x); L += (scheme === 'penalty' ? 0 : nu[j] * h) + 0.5 * mu * h * h; }
    return L;
  };
  const viol = (x) => { let v = 0; for (const g of ineq) v = Math.max(v, g(x)); for (const h of eq) v = Math.max(v, Math.abs(h(x))); return Math.max(0, v); };
  let x = x0 ? Array.from(x0) : bounds.map((b) => 0.5 * (b[0] + b[1]));
  const history = []; let innerHistory = []; let prevViol = Infinity;
  for (let k = 0; k < outerIter; k++) {
    let r;
    if (k === 0 && method !== 'nm') {
      if (method === 'de') r = differentialEvolution(merit, bounds, { seed, maxGen, tol: 1e-12 });
      else if (method === 'pso') r = particleSwarm(merit, bounds, { seed, maxIter: maxGen });
      else if (method === 'sa') r = simulatedAnnealing(merit, x, bounds, { seed, maxIter: maxGen * 20 * n });
      else if (method === 'ga') r = geneticAlgorithm(merit, bounds, { seed, generations: maxGen });
      else throw new RangeError("method must be 'de', 'pso', 'sa', 'ga' or 'nm'");
      innerHistory = Array.from(r.history || []);
      r = nelderMead(merit, r.x, { bounds, tol: 1e-12, maxIter: 400 * n });
    } else r = nelderMead(merit, x, { bounds, tol: 1e-12, maxIter: 400 * n });
    x = Array.from(r.x);
    const xs = snap(x); const v = viol(xs);
    history.push({ iter: k + 1, f: objective(xs), maxViolation: v, penalty: mu });
    if (scheme === 'auglag') {
      for (let i = 0; i < ineq.length; i++) lam[i] = Math.max(0, lam[i] + mu * ineq[i](xs));
      for (let j = 0; j < eq.length; j++) nu[j] += mu * eq[j](xs);
    }
    if (v <= tol && k > 0 && Math.abs(history[k].f - history[k - 1].f) <= tol * (1 + Math.abs(history[k].f))) break;
    if (ineq.length + eq.length === 0) break;
    if (v > 0.25 * prevViol || scheme === 'penalty') mu *= scheme === 'penalty' ? 10 : 4;
    prevViol = v;
  }
  const xs = snap(x); const v = viol(xs);
  return { x: xs, f: objective(xs), feasible: v <= Math.max(tol, 1e-6), maxViolation: v, ineqValues: ineq.map((g) => g(xs)),
    eqValues: eq.map((h) => h(xs)), multipliers: { ineq: lam, eq: nu }, history, innerHistory };
}

/**
 * Knee point of a Pareto front (all objectives minimised): objectives are scaled to 0–1 over the front; for two
 * objectives the knee is the point farthest below the line joining the two extreme solutions, otherwise the point
 * closest to the ideal (utopia) point.
 * Ref: Das, I. (1999), On characterizing the "knee" of the Pareto curve based on normal-boundary intersection,
 * Structural Optimization 18(2–3); Branke et al. (2004), Finding knees in multi-objective optimization, PPSN VIII.
 * @param {{x?: number[], f: number[]}[]} front non-dominated points
 * @returns {{index: number, point: object, normalised: number[][]}}
 */
export function kneePoint(front) {
  if (!Array.isArray(front) || front.length === 0) throw new RangeError('front must be a non-empty array');
  const m = front[0].f.length;
  const lo = new Array(m).fill(Infinity); const hi = new Array(m).fill(-Infinity);
  for (const p of front) for (let j = 0; j < m; j++) { lo[j] = Math.min(lo[j], p.f[j]); hi[j] = Math.max(hi[j], p.f[j]); }
  const N = front.map((p) => p.f.map((v, j) => (hi[j] > lo[j] ? (v - lo[j]) / (hi[j] - lo[j]) : 0)));
  let best = 0; let bestD = -Infinity;
  for (let i = 0; i < N.length; i++) {
    // two objectives: extremes are (0, 1) and (1, 0), the line is f1 + f2 = 1; distance below it ∝ 1 − f1 − f2
    const d = m === 2 ? 1 - N[i][0] - N[i][1] : -Math.hypot(...N[i]);
    if (d > bestD) { bestD = d; best = i; }
  }
  return { index: best, point: front[best], normalised: N };
}

/**
 * Multi-objective optimisation (NSGA-II) for trade-offs such as cost – risk – carbon, returning the Pareto set,
 * the extreme solution for each objective and the knee point.
 * Ref: Deb, Pratap, Agarwal & Meyarivan (2002), A fast and elitist multiobjective genetic algorithm: NSGA-II,
 * IEEE Trans. Evolutionary Computation 6(2).
 * @param {{objectives: (x: number[]) => number[], bounds: [number, number][], names?: string[], pop?: number,
 *   generations?: number, seed?: number}} o all objectives are minimised (negate to maximise)
 * @returns {{front: {x: number[], f: number[]}[], knee: {x: number[], f: number[]}, kneeIndex: number,
 *   extremes: {objective: string, x: number[], f: number[]}[], evaluations: number}}
 */
export function multiObjective({ objectives, bounds, names, pop = 80, generations = 120, seed = 1 }) {
  if (typeof objectives !== 'function') throw new RangeError('objectives must be a function returning an array');
  checkBounds(bounds);
  const r = nsga2(objectives, bounds, { pop, generations, seed });
  const front = r.front.map((p) => ({ x: Array.from(p.x), f: Array.from(p.f) }));
  const k = kneePoint(front);
  const m = front[0].f.length;
  const extremes = Array.from({ length: m }, (_, j) => { const p = front.reduce((b, q) => (q.f[j] < b.f[j] ? q : b)); return { objective: names?.[j] ?? `f${j + 1}`, x: p.x, f: p.f }; });
  return { front, knee: front[k.index], kneeIndex: k.index, extremes, evaluations: r.evaluations ?? pop * (generations + 1) };
}

/**
 * 0/1 knapsack (project selection under one budget) by dynamic programming on integer weights; non-integer costs
 * are scaled with `resolution` (cost units per DP cell) and rounded up, which keeps the solution feasible.
 * Ref: Bellman, R. (1957), Dynamic Programming; Martello & Toth (1990), Knapsack Problems, ch. 2.
 * @param {{values: number[], weights: number[], capacity: number, resolution?: number}} o values (e.g. NPV),
 *   weights (e.g. capex) and capacity (budget) in consistent units
 * @returns {{selected: number[], value: number, weight: number}} indices of the chosen items
 */
export function knapsack({ values, weights, capacity, resolution = 1 }) {
  requireSeries('values', values); requireSeries('weights', weights); requireNumber('capacity', capacity, { min: 0 });
  requireNumber('resolution', resolution, { gt: 0 });
  if (values.length !== weights.length) throw new RangeError('values and weights must have the same length');
  const n = values.length; const C = Math.floor(capacity / resolution + 1e-9);
  if (C > 2e6) throw new RangeError('capacity / resolution is too large for the DP table; increase resolution');
  const w = Array.from(weights, (v) => { if (v < 0) throw new RangeError('weights must be >= 0'); return Math.ceil(v / resolution - 1e-9); });
  const best = new Float64Array(C + 1); const take = Array.from({ length: n }, () => new Uint8Array(C + 1));
  for (let i = 0; i < n; i++) {
    if (values[i] <= 0) continue;
    for (let c = C; c >= w[i]; c--) { const v = best[c - w[i]] + values[i]; if (v > best[c]) { best[c] = v; take[i][c] = 1; } }
  }
  const selected = []; let c = C;
  for (let i = n - 1; i >= 0; i--) if (take[i][c]) { selected.push(i); c -= w[i]; }
  selected.reverse();
  return { selected, value: selected.reduce((s, i) => s + values[i], 0), weight: selected.reduce((s, i) => s + weights[i], 0) };
}

/**
 * Portfolio selection as a mixed-integer linear programme: maximise Σ value_i x_i subject to one or more budget
 * rows (e.g. capex per year), optional mutual-exclusion groups and dependencies (j requires i), x binary.
 * Ref: Weingartner, H.M. (1963), Mathematical Programming and the Analysis of Capital Budgeting Problems;
 * branch and bound of the shared MILP solver.
 * @param {{values: number[], budgets: {use: number[], limit: number}[], exclusive?: number[][],
 *   requires?: [number, number][]}} o budgets[k].use[i] = consumption of resource k by project i; exclusive =
 *   groups of which at most one may be chosen; requires = [j, i] pairs meaning j needs i
 * @returns {{selected: number[], value: number, usage: number[], status: string, nodes: number}}
 */
export function portfolioMilp({ values, budgets, exclusive = [], requires = [] }) {
  requireSeries('values', values);
  const n = values.length; const Aub = []; const bub = [];
  for (const b of budgets || []) { requireSeries('budget use', b.use); Aub.push(Array.from(b.use)); bub.push(requireNumber('budget limit', b.limit)); }
  for (const g of exclusive) { const row = new Array(n).fill(0); for (const i of g) row[i] = 1; Aub.push(row); bub.push(1); }
  for (const [j, i] of requires) { const row = new Array(n).fill(0); row[j] = 1; row[i] = -1; Aub.push(row); bub.push(0); }
  const r = milp({ c: Array.from(values, (v) => -v), Aub, bub, bounds: Array.from({ length: n }, () => [0, 1]), integer: new Array(n).fill(true) });
  if (!r.x) return { selected: [], value: 0, usage: (budgets || []).map(() => 0), status: r.status, nodes: r.nodes };
  const selected = []; r.x.forEach((v, i) => { if (v > 0.5) selected.push(i); });
  return { selected, value: selected.reduce((s, i) => s + values[i], 0),
    usage: (budgets || []).map((b) => selected.reduce((s, i) => s + b.use[i], 0)), status: r.status, nodes: r.nodes };
}

/**
 * Capital rationing: rank projects by profitability index (NPV per unit of capital) and fund in that order within
 * the budget, then compare with the exact integer optimum — the ranking rule can fail when projects are lumpy.
 * Ref: Brealey, Myers & Allen, Principles of Corporate Finance, ch. 5 (capital rationing and the PI rule);
 * Lorie & Savage (1955), J. Business 28(4).
 * @param {{projects: {name: string, cost: number, npv: number}[], budget: number}} o currency
 * @returns {{ranking: {name: string, pi: number, cost: number, npv: number, funded: boolean}[], greedy: {names: string[],
 *   npv: number, cost: number}, optimal: {names: string[], npv: number, cost: number}, rankingIsOptimal: boolean}}
 *   pi = NPV / cost
 */
export function capitalRationing({ projects, budget }) {
  if (!Array.isArray(projects) || projects.length === 0) throw new RangeError('projects must be a non-empty array');
  requireNumber('budget', budget, { min: 0 });
  const ranking = projects.map((p) => ({ name: p.name, pi: requireNumber('npv', p.npv) / requireNumber('cost', p.cost, { gt: 0 }), cost: p.cost, npv: p.npv, funded: false }))
    .sort((a, b) => b.pi - a.pi);
  let left = budget; const g = { names: [], npv: 0, cost: 0 };
  for (const r of ranking) if (r.npv > 0 && r.cost <= left + 1e-9) { r.funded = true; left -= r.cost; g.names.push(r.name); g.npv += r.npv; g.cost += r.cost; }
  const m = portfolioMilp({ values: projects.map((p) => p.npv), budgets: [{ use: projects.map((p) => p.cost), limit: budget }] });
  const optimal = { names: m.selected.map((i) => projects[i].name), npv: m.value, cost: m.usage[0] };
  return { ranking, greedy: g, optimal, rankingIsOptimal: Math.abs(optimal.npv - g.npv) <= 1e-9 * Math.max(1, Math.abs(optimal.npv)) };
}

/**
 * Decision under deep uncertainty from a scenario table of payoffs: regret table and the min–max regret choice,
 * with maximin (Wald), maximax, Laplace (equal likelihood) and Hurwicz results for comparison, and expected value
 * when probabilities are given.
 * Ref: Savage, L.J. (1951), The theory of statistical decision, J. American Statistical Association 46(253);
 * Kouvelis & Yu (1997), Robust Discrete Optimization and Its Applications, ch. 1.
 * @param {{alternatives: string[], scenarios: string[], payoffs: number[][], sense?: 'max' | 'min',
 *   probabilities?: number[], hurwiczAlpha?: number}} o payoffs[alternative][scenario]; sense 'max' for value
 *   (NPV), 'min' for cost; hurwiczAlpha = optimism coefficient 0–1 (default 0.5)
 * @returns {{regret: number[][], maxRegret: number[], minimaxRegret: string, maximin: string, maximax: string,
 *   laplace: string, hurwicz: string, expectedValue: string | null}}
 */
export function minMaxRegret({ alternatives, scenarios, payoffs, sense = 'max', probabilities, hurwiczAlpha = 0.5 }) {
  const na = alternatives.length; const ns = scenarios.length;
  if (!Array.isArray(payoffs) || payoffs.length !== na) throw new RangeError('payoffs needs one row per alternative');
  const sg = sense === 'max' ? 1 : sense === 'min' ? -1 : NaN;
  if (Number.isNaN(sg)) throw new RangeError("sense must be 'max' or 'min'");
  const U = payoffs.map((r) => { if (r.length !== ns) throw new RangeError('every payoff row needs one value per scenario'); return r.map((v) => sg * requireNumber('payoff', v)); });
  const bestS = Array.from({ length: ns }, (_, s) => Math.max(...U.map((r) => r[s])));
  const regret = U.map((r) => r.map((v, s) => bestS[s] - v));
  const maxRegret = regret.map((r) => Math.max(...r));
  const arg = (a, f) => alternatives[a.indexOf(f(...a))];
  const worst = U.map((r) => Math.min(...r)); const bestA = U.map((r) => Math.max(...r));
  const avg = U.map((r) => r.reduce((s, v) => s + v, 0) / ns);
  requireNumber('hurwiczAlpha', hurwiczAlpha, { min: 0, max: 1 });
  const hz = U.map((_, a) => hurwiczAlpha * bestA[a] + (1 - hurwiczAlpha) * worst[a]);
  let ev = null;
  if (probabilities) {
    const sp = probabilities.reduce((s, v) => s + v, 0);
    if (probabilities.length !== ns || Math.abs(sp - 1) > 1e-9) throw new RangeError('probabilities must match the scenarios and sum to 1');
    ev = arg(U.map((r) => r.reduce((s, v, k) => s + v * probabilities[k], 0)), Math.max);
  }
  return { regret, maxRegret, minimaxRegret: arg(maxRegret, Math.min), maximin: arg(worst, Math.max), maximax: arg(bestA, Math.max),
    laplace: arg(avg, Math.max), hurwicz: arg(hz, Math.max), expectedValue: ev };
}

/**
 * Chance-constrained optimisation by sample-average approximation: minimise f(x) subject to
 * P[g(x, ξ) <= 0] >= p, with the probability estimated on a fixed (common-random-number) Latin-hypercube sample of
 * the uncertain inputs ξ and enforced by a penalty inside differential evolution.
 * Ref: Charnes, A. & Cooper, W.W. (1959), Chance-constrained programming, Management Science 6(1); Luedtke &
 * Ahmed (2008), A sample approximation approach for optimization with probabilistic constraints, SIAM J. Optim. 19(2).
 * @param {{objective: (x: number[]) => number, constraint: (x: number[], xi: object) => number,
 *   vars: Record<string, object>, bounds: [number, number][], p: number, n?: number, seed?: number,
 *   maxGen?: number}} o constraint satisfied when <= 0; vars = distributions of ξ; p = required probability
 * @returns {{x: number[], f: number, probability: number, standardError: number, feasible: boolean,
 *   history: number[]}} probability = estimated satisfaction probability at the solution
 */
export function chanceConstrained({ objective, constraint, vars, bounds, p, n = 400, seed = 1, maxGen = 120 }) {
  checkBounds(bounds); requireNumber('p', p, { gt: 0, max: 1 });
  const s = sampleInputs({ vars, n, seed });
  const xi = Array.from({ length: n }, (_, k) => { const o = {}; for (const nm of s.names) o[nm] = s.samples[nm][k]; return o; });
  const prob = (x) => { let c = 0; for (let k = 0; k < n; k++) if (constraint(x, xi[k]) <= 0) c++; return c / n; };
  // scale the penalty to the objective's magnitude at the centre of the box
  const f0 = Math.abs(objective(bounds.map((b) => 0.5 * (b[0] + b[1])))) + 1;
  const merit = (x) => { const pr = prob(x); return objective(x) + (pr < p ? f0 * 1e3 * (1 + (p - pr)) : 0); };
  const r = differentialEvolution(merit, bounds, { seed, maxGen });
  const pr = prob(r.x);
  return { x: Array.from(r.x), f: objective(r.x), probability: pr, standardError: Math.sqrt(pr * (1 - pr) / n), feasible: pr >= p,
    history: Array.from(r.history || []) };
}

/**
 * Two-stage stochastic linear programme with recourse solved in extensive form: minimise cᵀx + Σ_s p_s q_sᵀ y_s
 * subject to A x <= b and T_s x + W_s y_s <= h_s (or = h_s with sense 'eq') for every scenario, x, y >= 0.
 * Also returns the wait-and-see value and EVPI = RP − WS.
 * Ref: Dantzig, G.B. (1955), Linear programming under uncertainty, Management Science 1(3–4); Birge & Louveaux,
 * Introduction to Stochastic Programming, 2nd ed.
 * @param {{c: number[], A?: number[][], b?: number[], scenarios: {p: number, q: number[], T: number[][],
 *   W: number[][], h: number[], sense?: 'le' | 'eq'}[], xBounds?: [number, number][]}} o first-stage costs and
 *   constraints; per scenario: probability, recourse costs, technology and recourse matrices, right-hand side
 * @returns {{x: number[], y: number[][], objective: number, firstStageCost: number, expectedRecourse: number,
 *   status: string, waitAndSee: number | null, evpi: number | null}}
 */
export function twoStageStochastic({ c, A = [], b = [], scenarios, xBounds }) {
  requireSeries('c', c);
  if (!Array.isArray(scenarios) || scenarios.length === 0) throw new RangeError('scenarios must be a non-empty array');
  const sp = scenarios.reduce((s, sc) => s + requireNumber('scenario p', sc.p, { min: 0, max: 1 }), 0);
  if (Math.abs(sp - 1) > 1e-9) throw new RangeError('scenario probabilities must sum to 1');
  const n1 = c.length;
  const build = (list, weights) => {
    const offs = []; let nv = n1;
    for (const sc of list) { offs.push(nv); nv += sc.q.length; }
    const cc = new Array(nv).fill(0); for (let j = 0; j < n1; j++) cc[j] = c[j];
    const Aub = []; const bub = []; const Aeq = []; const beq = [];
    A.forEach((row, i) => { const r = new Array(nv).fill(0); for (let j = 0; j < n1; j++) r[j] = row[j]; Aub.push(r); bub.push(b[i]); });
    list.forEach((sc, s) => {
      for (let j = 0; j < sc.q.length; j++) cc[offs[s] + j] = weights[s] * sc.q[j];
      sc.h.forEach((hv, i) => {
        const r = new Array(nv).fill(0);
        for (let j = 0; j < n1; j++) r[j] = sc.T[i][j];
        for (let j = 0; j < sc.q.length; j++) r[offs[s] + j] = sc.W[i][j];
        if (sc.sense === 'eq') { Aeq.push(r); beq.push(hv); } else { Aub.push(r); bub.push(hv); }
      });
    });
    const bounds = Array.from({ length: nv }, (_, j) => (j < n1 && xBounds ? xBounds[j] : [0, Infinity]));
    const prob = { c: cc, bounds };
    if (Aub.length) { prob.Aub = Aub; prob.bub = bub; }
    if (Aeq.length) { prob.Aeq = Aeq; prob.beq = beq; }
    return { r: linprog(prob), offs };
  };
  const { r, offs } = build(scenarios, scenarios.map((s) => s.p));
  if (r.status !== 'optimal') return { x: [], y: [], objective: r.f, firstStageCost: NaN, expectedRecourse: NaN, status: r.status, waitAndSee: null, evpi: null };
  const x = Array.from(r.x).slice(0, n1);
  const y = scenarios.map((sc, s) => Array.from(r.x).slice(offs[s], offs[s] + sc.q.length));
  const first = x.reduce((s, v, j) => s + v * c[j], 0);
  let ws = 0; let wsOk = true;
  for (const sc of scenarios) { const one = build([sc], [1]).r; if (one.status !== 'optimal') { wsOk = false; break; } ws += sc.p * one.f; }
  return { x, y, objective: r.f, firstStageCost: first, expectedRecourse: r.f - first, status: r.status,
    waitAndSee: wsOk ? ws : null, evpi: wsOk ? r.f - ws : null };
}

/**
 * Deterministic finite-horizon dynamic programme by backward induction (staged development, phasing, expansion
 * timing): maximises the discounted sum of rewards. States are identified by string keys.
 * Ref: Bellman, R. (1957), Dynamic Programming (principle of optimality).
 * @param {{stages: number, initial: string, actions: (stage: number, state: string) => {name: string, reward: number,
 *   next: string}[], terminal?: (state: string) => number, discount?: number}} o discount = factor per stage
 *   (default 1); actions returns the feasible decisions at a stage/state (empty = dead end)
 * @returns {{value: number, path: {stage: number, state: string, action: string, reward: number}[],
 *   statesEvaluated: number}}
 */
export function dpStaged({ stages, initial, actions, terminal = () => 0, discount = 1 }) {
  requireNumber('stages', stages, { min: 1, integer: true }); requireNumber('discount', discount, { gt: 0 });
  const memo = new Map();
  const solve = (t, s) => {
    if (t === stages) return { v: terminal(s), a: null };
    const key = `${t}|${s}`;
    if (memo.has(key)) return memo.get(key);
    let best = { v: -Infinity, a: null };
    for (const a of actions(t, s) || []) {
      const v = requireNumber('reward', a.reward) + discount * solve(t + 1, a.next).v;
      if (v > best.v) best = { v, a };
    }
    memo.set(key, best);
    return best;
  };
  const root = solve(0, initial);
  const path = []; let s = initial;
  for (let t = 0; t < stages; t++) { const b = solve(t, s); if (!b.a) break; path.push({ stage: t, state: s, action: b.a.name, reward: b.a.reward }); s = b.a.next; }
  return { value: root.v, path, statesEvaluated: memo.size };
}

/**
 * Equipment replacement by dynamic programming: each year keep (pay O&M for the current age) or replace (buy new
 * less trade-in, pay O&M of a new unit); minimises the present value of cost over the horizon, with salvage of
 * the final unit credited at the end.
 * Ref: Bellman (1955), Equipment replacement policy, J. SIAM 3(3); Hillier & Lieberman, Introduction to
 * Operations Research (dynamic programming).
 * @param {{horizon: number, maxAge: number, price: number, om: (age: number) => number, salvage: (age: number) => number,
 *   rate?: number, initialAge?: number}} o horizon (years), maximum allowed age, price of a new unit (currency),
 *   O&M cost of a unit of given age during the year, trade-in value at a given age, discount rate
 * @returns {{pvCost: number, plan: {year: number, age: number, action: 'keep' | 'replace'}[], replacements: number[]}}
 *   replacements = years (1-based) at whose start a replacement is made
 */
export function equipmentReplacementDp({ horizon, maxAge, price, om, salvage, rate = 0, initialAge = 0 }) {
  requireNumber('horizon', horizon, { min: 1, integer: true }); requireNumber('maxAge', maxAge, { min: 1, integer: true });
  requireNumber('price', price, { min: 0 }); requireNumber('rate', rate, { gt: -1 });
  requireNumber('initialAge', initialAge, { min: 0, max: maxAge, integer: true });
  const d = 1 / (1 + rate);
  const r = dpStaged({ stages: horizon, initial: String(initialAge), discount: d,
    actions: (t, s) => {
      const age = Number(s); const acts = [];
      if (age < maxAge) acts.push({ name: 'keep', reward: -om(age) * d, next: String(age + 1) });
      if (age > 0) acts.push({ name: 'replace', reward: -(price - salvage(age)) - om(0) * d, next: '1' });
      return acts;
    },
    terminal: (s) => salvage(Number(s)) });
  const plan = r.path.map((p) => ({ year: p.stage + 1, age: Number(p.state), action: p.action }));
  return { pvCost: -r.value, plan, replacements: plan.filter((p) => p.action === 'replace').map((p) => p.year) };
}

/**
 * Maintenance / inspection scheduling over a horizon as a finite-horizon Markov decision process: the asset
 * deteriorates through condition states with transition matrix P; each period one action is chosen (e.g. do
 * nothing, repair, replace), which costs `cost[state]`, moves the asset to `next[state]` and is followed by one
 * period of deterioration and the state-dependent operating/risk cost. Backward induction gives the cost-optimal
 * policy for every period and state.
 * Ref: Howard, R.A. (1960), Dynamic Programming and Markov Processes; Puterman, Markov Decision Processes, ch. 4;
 * application to inspection and maintenance planning as in Frangopol, Kallen & van Noortwijk (2004), Progress in
 * Structural Engineering and Materials 6(4).
 * @param {{horizon: number, P: number[][], stateCost: number[], actions: {name: string, cost: number[],
 *   next: number[]}[], rate?: number, terminalCost?: number[], initialState?: number}} o P rows sum to 1;
 *   stateCost = expected cost per period of being in each state (downtime, failure risk)
 * @returns {{expectedCost: number, policy: string[][], firstAction: string, value: number[][]}} policy[t][state];
 *   value[t][state] = expected present value of cost from period t
 */
export function maintenanceScheduleDp({ horizon, P, stateCost, actions, rate = 0, terminalCost, initialState = 0 }) {
  requireNumber('horizon', horizon, { min: 1, integer: true });
  const n = P.length; const d = 1 / (1 + requireNumber('rate', rate, { gt: -1 }));
  for (const row of P) { const s = row.reduce((a, v) => a + requireNumber('P', v, { min: 0, max: 1 }), 0); if (row.length !== n || Math.abs(s - 1) > 1e-9) throw new RangeError('P must be square with rows summing to 1'); }
  requireSeries('stateCost', stateCost);
  if (!Array.isArray(actions) || actions.length === 0) throw new RangeError('actions must be a non-empty array');
  const value = Array.from({ length: horizon + 1 }, () => new Array(n).fill(0));
  const policy = Array.from({ length: horizon }, () => new Array(n).fill(''));
  if (terminalCost) for (let s = 0; s < n; s++) value[horizon][s] = terminalCost[s];
  for (let t = horizon - 1; t >= 0; t--) {
    for (let s = 0; s < n; s++) {
      let best = Infinity; let name = '';
      for (const a of actions) {
        const s1 = a.next[s];
        if (s1 === undefined || s1 === null || !(a.cost[s] < Infinity)) continue;
        let ev = 0; for (let j = 0; j < n; j++) ev += P[s1][j] * (stateCost[j] + value[t + 1][j]);
        const v = a.cost[s] + d * ev;
        if (v < best) { best = v; name = a.name; }
      }
      value[t][s] = best; policy[t][s] = name;
    }
  }
  return { expectedCost: value[0][initialState], policy, firstAction: policy[0][initialState], value };
}

/**
 * Surrogate-assisted (Bayesian) optimisation for expensive engineering models: a Gaussian-process surrogate with
 * expected-improvement infill (shared library routine), reporting the best-so-far history.
 * Ref: Jones, Schonlau & Welch (1998), Efficient global optimization of expensive black-box functions, J. Global
 * Optimization 13(4).
 * @param {{f: (x: number[]) => number, bounds: [number, number][], nInit?: number, nIter?: number, seed?: number}} o
 *   objective to minimise; nInit initial design points, nIter adaptive evaluations
 * @returns {{x: number[], f: number, history: number[], evaluations: number}}
 */
export function surrogateOptimise({ f, bounds, nInit, nIter = 25, seed = 1 }) {
  checkBounds(bounds);
  const opts = { nIter, seed }; if (nInit !== undefined) opts.nInit = nInit;
  const r = bayesOptimise(f, bounds, opts);
  return { x: Array.from(r.x), f: r.f, history: Array.from(r.history || []), evaluations: (r.history || []).length };
}

/**
 * Linear production planner: maximise Σ margin_j x_j subject to resource rows Σ use_ij x_j <= capacity_i and
 * optional demand caps (product-mix problem).
 * Ref: Dantzig (1963), Linear Programming and Extensions.
 * @param {{margin: number[], use: number[][], capacity: number[], maxDemand?: number[]}} o margin per unit of each
 *   product (currency), use[resource][product], capacity per resource
 * @returns {{x: number[], profit: number, slack: number[], status: string}} production quantities and unused capacity
 */
export function productionPlan({ margin, use, capacity, maxDemand }) {
  requireSeries('margin', margin); requireSeries('capacity', capacity);
  if (!Array.isArray(use) || use.length !== capacity.length) throw new RangeError('use needs one row per resource');
  const n = margin.length;
  const r = linprog({ c: Array.from(margin, (v) => -v), Aub: use.map((row) => Array.from(row)), bub: Array.from(capacity),
    bounds: Array.from({ length: n }, (_, j) => [0, maxDemand ? maxDemand[j] : Infinity]) });
  if (r.status !== 'optimal') return { x: [], profit: NaN, slack: [], status: r.status };
  const x = Array.from(r.x);
  return { x, profit: -r.f, slack: use.map((row, i) => capacity[i] - row.reduce((s, v, j) => s + v * x[j], 0)), status: r.status };
}

/**
 * Transportation / logistics plan: minimise Σ cost_ij x_ij shipping from sources (supply limits) to sinks (demand
 * to be met).
 * Ref: Hitchcock, F.L. (1941), The distribution of a product from several sources to numerous localities,
 * J. Mathematics and Physics 20; Dantzig (1963).
 * @param {{supply: number[], demand: number[], cost: number[][]}} o cost[source][sink] per unit
 * @returns {{flow: number[][], cost: number, status: string}}
 */
export function transportPlan({ supply, demand, cost }) {
  requireSeries('supply', supply); requireSeries('demand', demand);
  const m = supply.length; const n = demand.length;
  if (!Array.isArray(cost) || cost.length !== m) throw new RangeError('cost needs one row per source');
  const c = []; for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) c.push(requireNumber('cost', cost[i][j]));
  const Aub = []; const bub = []; const Aeq = []; const beq = [];
  for (let i = 0; i < m; i++) { const r = new Array(m * n).fill(0); for (let j = 0; j < n; j++) r[i * n + j] = 1; Aub.push(r); bub.push(supply[i]); }
  for (let j = 0; j < n; j++) { const r = new Array(m * n).fill(0); for (let i = 0; i < m; i++) r[i * n + j] = 1; Aeq.push(r); beq.push(demand[j]); }
  const r = linprog({ c, Aub, bub, Aeq, beq });
  if (r.status !== 'optimal') return { flow: [], cost: NaN, status: r.status };
  return { flow: Array.from({ length: m }, (_, i) => Array.from(r.x).slice(i * n, (i + 1) * n)), cost: r.f, status: r.status };
}
