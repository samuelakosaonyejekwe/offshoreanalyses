// Multi-criteria decision analysis (weighted sum, TOPSIS, AHP, ELECTRE I, Pareto filtering, weight stability),
// decision trees with expected-value / expected-utility rollback, alternative comparison tables and a
// deterministic, template-based recommendation composer.

import { requireNumber } from './tvm.js';

function checkMatrix(matrix, nCrit) {
  if (!Array.isArray(matrix) || matrix.length === 0) throw new RangeError('matrix must have one row per alternative');
  for (const row of matrix) {
    if (!row || row.length !== nCrit) throw new RangeError('every matrix row needs one value per criterion');
    for (const v of row) requireNumber('matrix value', v);
  }
}
function normWeights(criteria) {
  let s = 0;
  for (const c of criteria) s += requireNumber(`weight of '${c.name}'`, c.weight, { min: 0 });
  if (!(s > 0)) throw new RangeError('weights must sum to a positive number');
  return criteria.map((c) => c.weight / s);
}
const isMax = (c) => {
  if (c.sense !== 'max' && c.sense !== 'min') throw new RangeError(`sense of '${c.name}' must be 'max' or 'min'`);
  return c.sense === 'max';
};

/**
 * Normalise a decision matrix column by column to a 0–1 "higher is better" scale.
 * 'minmax': (x − min)/(max − min) or (max − x)/(max − min); 'max': x/max or min/x; 'sum': x/Σx or (1/x)/Σ(1/x);
 * 'vector': x/√Σx² (cost criteria: 1 − that).
 * Ref: Hwang, C.L. & Yoon, K. (1981), Multiple Attribute Decision Making (normalisation methods).
 * @param {number[][]} matrix matrix[alternative][criterion] in the criteria's own units
 * @param {{name: string, sense: 'max' | 'min'}[]} criteria
 * @param {'minmax' | 'max' | 'sum' | 'vector'} [method] default 'minmax'
 * @returns {number[][]} dimensionless scores
 */
export function normaliseMatrix(matrix, criteria, method = 'minmax') {
  checkMatrix(matrix, criteria.length);
  const out = matrix.map((r) => r.slice());
  criteria.forEach((c, j) => {
    const col = matrix.map((r) => r[j]); const mx = Math.max(...col); const mn = Math.min(...col); const up = isMax(c);
    for (let i = 0; i < matrix.length; i++) {
      const x = col[i]; let v;
      if (method === 'minmax') v = mx === mn ? 1 : up ? (x - mn) / (mx - mn) : (mx - x) / (mx - mn);
      else if (method === 'max') {
        if (!(mn > 0)) throw new RangeError("'max' normalisation needs positive values");
        v = up ? x / mx : mn / x;
      } else if (method === 'sum') {
        if (!(mn > 0)) throw new RangeError("'sum' normalisation needs positive values");
        v = up ? x / col.reduce((s, y) => s + y, 0) : (1 / x) / col.reduce((s, y) => s + 1 / y, 0);
      } else if (method === 'vector') {
        const nrm = Math.sqrt(col.reduce((s, y) => s + y * y, 0)) || 1;
        v = up ? x / nrm : 1 - x / nrm;
      } else throw new RangeError("method must be 'minmax', 'max', 'sum' or 'vector'");
      out[i][j] = v;
    }
  });
  return out;
}

/**
 * Weighted-sum (simple additive weighting) scores and ranking.
 * Ref: Fishburn, P.C. (1967), Additive utilities with incomplete product sets; Hwang & Yoon (1981).
 * @param {{alternatives: string[], criteria: {name: string, weight: number, sense: 'max' | 'min'}[],
 *   matrix: number[][], normalisation?: string}} o weights on any positive scale (normalised internally)
 * @returns {{scores: number[], ranking: {name: string, score: number, rank: number}[], weights: number[],
 *   normalised: number[][], best: string}}
 */
export function weightedSum({ alternatives, criteria, matrix, normalisation = 'minmax' }) {
  const w = normWeights(criteria);
  const N = normaliseMatrix(matrix, criteria, normalisation);
  if (alternatives.length !== matrix.length) throw new RangeError('alternatives and matrix rows must match');
  const scores = N.map((r) => r.reduce((s, v, j) => s + v * w[j], 0));
  const ranking = alternatives.map((name, i) => ({ name, score: scores[i], rank: 0 })).sort((a, b) => b.score - a.score);
  ranking.forEach((r, k) => { r.rank = k + 1; });
  return { scores, ranking, weights: w, normalised: N, best: ranking[0].name };
}

/**
 * TOPSIS: vector normalisation, weighting, distances to the ideal and anti-ideal solutions and relative closeness
 * C_i = d⁻ / (d⁺ + d⁻); the highest closeness is preferred.
 * Ref: Hwang, C.L. & Yoon, K. (1981), Multiple Attribute Decision Making: Methods and Applications.
 * @param {{alternatives: string[], criteria: {name: string, weight: number, sense: 'max' | 'min'}[],
 *   matrix: number[][]}} o
 * @returns {{closeness: number[], dPlus: number[], dMinus: number[], ideal: number[], antiIdeal: number[],
 *   ranking: {name: string, closeness: number, rank: number}[], best: string}}
 */
export function topsis({ alternatives, criteria, matrix }) {
  checkMatrix(matrix, criteria.length);
  if (alternatives.length !== matrix.length) throw new RangeError('alternatives and matrix rows must match');
  const w = normWeights(criteria); const m = matrix.length; const n = criteria.length;
  const V = matrix.map((r) => r.slice());
  const ideal = new Array(n); const antiIdeal = new Array(n);
  for (let j = 0; j < n; j++) {
    const nrm = Math.sqrt(matrix.reduce((s, r) => s + r[j] * r[j], 0)) || 1;
    for (let i = 0; i < m; i++) V[i][j] = matrix[i][j] / nrm * w[j];
    const col = V.map((r) => r[j]); const up = isMax(criteria[j]);
    ideal[j] = up ? Math.max(...col) : Math.min(...col); antiIdeal[j] = up ? Math.min(...col) : Math.max(...col);
  }
  const dPlus = V.map((r) => Math.sqrt(r.reduce((s, v, j) => s + (v - ideal[j]) ** 2, 0)));
  const dMinus = V.map((r) => Math.sqrt(r.reduce((s, v, j) => s + (v - antiIdeal[j]) ** 2, 0)));
  const closeness = dPlus.map((d, i) => (d + dMinus[i] > 0 ? dMinus[i] / (d + dMinus[i]) : 0.5));
  const ranking = alternatives.map((name, i) => ({ name, closeness: closeness[i], rank: 0 })).sort((a, b) => b.closeness - a.closeness);
  ranking.forEach((r, k) => { r.rank = k + 1; });
  return { closeness, dPlus, dMinus, ideal, antiIdeal, ranking, best: ranking[0].name };
}

/**
 * Saaty's random consistency index RI by matrix order (average consistency index of random reciprocal matrices).
 * Ref: Saaty, T.L. (1980), The Analytic Hierarchy Process (n = 1…10); extended values n = 11…15 from
 * Saaty's later tabulation.
 */
export const SAATY_RANDOM_INDEX = Object.freeze({ 1: 0, 2: 0, 3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24, 7: 1.32, 8: 1.41,
  9: 1.45, 10: 1.49, 11: 1.51, 12: 1.48, 13: 1.56, 14: 1.57, 15: 1.59 });

/**
 * Analytic hierarchy process priorities from a reciprocal pairwise-comparison matrix: principal right eigenvector
 * (power iteration), λmax, consistency index CI = (λmax − n)/(n − 1) and consistency ratio CR = CI/RI.
 * CR < 0.10 is the customary acceptance limit.
 * Ref: Saaty, T.L. (1977), A scaling method for priorities in hierarchical structures, J. Mathematical
 * Psychology 15(3); Saaty (1980), The Analytic Hierarchy Process.
 * @param {number[][]} pairwise n × n positive matrix with a_ji = 1/a_ij and a_ii = 1 (n <= 15)
 * @returns {{weights: number[], lambdaMax: number, ci: number, cr: number, consistent: boolean, n: number}}
 */
export function ahp(pairwise) {
  const n = Array.isArray(pairwise) ? pairwise.length : 0;
  if (n < 1 || n > 15) throw new RangeError('pairwise must be an n × n matrix with 1 <= n <= 15');
  for (let i = 0; i < n; i++) {
    if (!pairwise[i] || pairwise[i].length !== n) throw new RangeError('pairwise must be square');
    for (let j = 0; j < n; j++) {
      requireNumber('pairwise value', pairwise[i][j], { gt: 0 });
      if (Math.abs(pairwise[i][j] * pairwise[j][i] - 1) > 1e-6) throw new RangeError(`pairwise is not reciprocal at (${i + 1}, ${j + 1})`);
    }
  }
  let w = new Array(n).fill(1 / n); let lambda = n;
  for (let it = 0; it < 1000; it++) {
    const y = pairwise.map((row) => row.reduce((s, a, j) => s + a * w[j], 0));
    const sum = y.reduce((s, v) => s + v, 0);
    const next = y.map((v) => v / sum);
    let diff = 0; for (let i = 0; i < n; i++) diff = Math.max(diff, Math.abs(next[i] - w[i]));
    w = next; lambda = sum;                                   // with Σw = 1, Σ(Aw) converges to λmax
    if (diff < 1e-14) break;
  }
  const ci = n > 1 ? (lambda - n) / (n - 1) : 0;
  const ri = SAATY_RANDOM_INDEX[n];
  const cr = ri > 0 ? ci / ri : 0;
  return { weights: w, lambdaMax: lambda, ci: Math.abs(ci) < 1e-12 ? 0 : ci, cr: Math.abs(cr) < 1e-12 ? 0 : cr, consistent: cr < 0.1, n };
}

/**
 * ELECTRE I outranking: a outranks b when the concordance (share of weight on criteria where a is at least as
 * good) reaches the concordance threshold and the discordance (largest normalised amount by which b beats a) stays
 * below the discordance threshold. The kernel is the set of alternatives not outranked by any other.
 * Ref: Roy, B. (1968), Classement et choix en présence de points de vue multiples (la méthode ELECTRE), RIRO 2(8).
 * @param {{alternatives: string[], criteria: {name: string, weight: number, sense: 'max' | 'min'}[],
 *   matrix: number[][], concordanceThreshold?: number, discordanceThreshold?: number}} o thresholds in 0–1
 * @returns {{concordance: number[][], discordance: number[][], outranks: boolean[][], kernel: string[]}}
 */
export function electre1({ alternatives, criteria, matrix, concordanceThreshold = 0.65, discordanceThreshold = 0.35 }) {
  checkMatrix(matrix, criteria.length);
  const w = normWeights(criteria); const m = matrix.length; const n = criteria.length;
  const sgn = criteria.map((c) => (isMax(c) ? 1 : -1));
  const range = criteria.map((_, j) => { const col = matrix.map((r) => r[j]); return Math.max(...col) - Math.min(...col) || 1; });
  const C = Array.from({ length: m }, () => new Array(m).fill(0)); const D = Array.from({ length: m }, () => new Array(m).fill(0));
  const O = Array.from({ length: m }, () => new Array(m).fill(false));
  for (let a = 0; a < m; a++) for (let b = 0; b < m; b++) {
    if (a === b) continue;
    let c = 0; let d = 0;
    for (let j = 0; j < n; j++) {
      const diff = sgn[j] * (matrix[a][j] - matrix[b][j]);
      if (diff >= 0) c += w[j]; else d = Math.max(d, -diff / range[j]);
    }
    C[a][b] = c; D[a][b] = d; O[a][b] = c >= concordanceThreshold && d <= discordanceThreshold;
  }
  const kernel = alternatives.filter((_, b) => !O.some((row, a) => a !== b && row[b] && !O[b][a]));
  return { concordance: C, discordance: D, outranks: O, kernel };
}

/**
 * Pareto (non-dominated) filter: a point is dominated when another is at least as good in every objective and
 * strictly better in one.
 * Ref: Deb, K. (2001), Multi-Objective Optimization using Evolutionary Algorithms.
 * @param {number[][]} points points[i][objective]
 * @param {('max' | 'min')[]} senses sense of each objective
 * @returns {{indices: number[], dominated: number[]}} indices of non-dominated and dominated points
 */
export function paretoFilter(points, senses) {
  if (!Array.isArray(points) || points.length === 0) throw new RangeError('points must be a non-empty array');
  const s = senses.map((x) => { if (x !== 'max' && x !== 'min') throw new RangeError("senses must be 'max' or 'min'"); return x === 'max' ? 1 : -1; });
  const indices = []; const dominated = [];
  for (let i = 0; i < points.length; i++) {
    let dom = false;
    for (let k = 0; k < points.length && !dom; k++) {
      if (k === i) continue;
      let ge = true; let gt = false;
      for (let j = 0; j < s.length; j++) {
        const d = s[j] * (points[k][j] - points[i][j]);
        if (d < 0) { ge = false; break; }
        if (d > 0) gt = true;
      }
      dom = ge && gt;
    }
    (dom ? dominated : indices).push(i);
  }
  return { indices, dominated };
}

/**
 * Weight-stability intervals for the weighted-sum ranking: for each criterion, the range of its weight (the
 * other weights rescaled proportionally) over which the currently best alternative stays best.
 * Ref: Triantaphyllou, E. & Sánchez, A. (1997), A sensitivity analysis approach for some deterministic
 * multi-criteria decision-making methods, Decision Sciences 28(1).
 * @param {{alternatives: string[], criteria: object[], matrix: number[][], normalisation?: string}} o as weightedSum
 * @returns {{best: string, intervals: {criterion: string, weight: number, lower: number, upper: number,
 *   robust: boolean}[]}} weights normalised to sum 1; robust = the full 0–1 range keeps the winner
 */
export function weightStability({ alternatives, criteria, matrix, normalisation = 'minmax' }) {
  const ws = weightedSum({ alternatives, criteria, matrix, normalisation });
  const { weights: w, normalised: N, scores } = ws;
  const best = scores.indexOf(Math.max(...scores));
  const intervals = criteria.map((c, j) => {
    const rest = 1 - w[j];
    const S = N.map((r) => (rest > 0 ? r.reduce((s, v, k) => s + (k === j ? 0 : v * w[k]), 0) / rest : 0));
    let lo = 0; let hi = 1;
    for (let b = 0; b < N.length; b++) {
      if (b === best) continue;
      const A = N[best][j] - N[b][j]; const B = S[best] - S[b];      // lead(t) = B + (A − B) t must stay >= 0
      const slope = A - B;
      if (Math.abs(slope) < 1e-15) { if (B < 0) { lo = 1; hi = 0; } continue; }
      const t0 = -B / slope;
      if (slope > 0) lo = Math.max(lo, t0); else hi = Math.min(hi, t0);
    }
    lo = Math.min(Math.max(lo, 0), 1); hi = Math.max(Math.min(hi, 1), 0);
    return { criterion: c.name, weight: w[j], lower: Math.min(lo, w[j]), upper: Math.max(hi, w[j]), robust: lo <= 1e-12 && hi >= 1 - 1e-12 };
  });
  return { best: alternatives[best], intervals };
}

// ---------------------------------------------------------------------------------------------------------------
// Decision trees
// ---------------------------------------------------------------------------------------------------------------

/**
 * Roll back a decision tree. Nodes: { type: 'decision', name, options: [{ name, payoff?, node }] },
 * { type: 'chance', name, branches: [{ name, p, payoff?, node }] }, { type: 'terminal', value }.
 * Branch payoffs (negative for costs) accumulate along the path and are added to the terminal value.
 * Without `riskTolerance` decisions maximise expected value. With it, they maximise expected exponential utility
 * u(x) = 1 − exp(−x/R) of the total path payoff and the certainty equivalent CE = −R ln(1 − E[u]) is reported.
 * Ref: Raiffa, H. (1968), Decision Analysis; Clemen & Reilly, Making Hard Decisions, ch. 4 and 14 (exponential
 * utility and risk tolerance).
 * @param {object} tree root node
 * @param {{riskTolerance?: number}} [o] R in currency (> 0)
 * @returns {{expectedValue: number, certaintyEquivalent: number, riskPremium: number,
 *   policy: {node: string, choice: string, expectedValue: number}[]}} values in currency; policy lists the chosen
 *   option at every decision node reached under the optimal strategy
 */
export function rollbackTree(tree, { riskTolerance } = {}) {
  if (riskTolerance !== undefined) requireNumber('riskTolerance', riskTolerance, { gt: 0 });
  const R = riskTolerance;
  const u = (x) => (R ? 1 - Math.exp(-x / R) : x);
  const visit = (node, acc, path) => {
    if (!node || typeof node !== 'object') throw new RangeError(`tree node missing at ${path}`);
    if (node.type === 'terminal') {
      const x = acc + requireNumber(`terminal value at ${path}`, node.value);
      return { ev: x, eu: u(x), policy: [] };
    }
    if (node.type === 'chance') {
      if (!Array.isArray(node.branches) || node.branches.length === 0) throw new RangeError(`chance node ${path} needs branches`);
      let sp = 0; let ev = 0; let eu = 0; const policy = [];
      for (const b of node.branches) {
        requireNumber(`probability at ${path}`, b.p, { min: 0, max: 1 });
        const r = visit(b.node, acc + (b.payoff ?? 0), `${path}/${b.name ?? '?'}`);
        sp += b.p; ev += b.p * r.ev; eu += b.p * r.eu;
        if (b.p > 0) policy.push(...r.policy);
      }
      if (Math.abs(sp - 1) > 1e-9) throw new RangeError(`probabilities at chance node ${path} must sum to 1`);
      return { ev, eu, policy };
    }
    if (node.type === 'decision') {
      if (!Array.isArray(node.options) || node.options.length === 0) throw new RangeError(`decision node ${path} needs options`);
      let best = null; let bestName = '';
      for (const opt of node.options) {
        const r = visit(opt.node, acc + (opt.payoff ?? 0), `${path}/${opt.name ?? '?'}`);
        if (!best || r.eu > best.eu) { best = r; bestName = opt.name ?? '?'; }
      }
      return { ev: best.ev, eu: best.eu, policy: [{ node: node.name ?? path, choice: bestName, expectedValue: best.ev }, ...best.policy] };
    }
    throw new RangeError(`node type '${node.type}' at ${path} is not decision, chance or terminal`);
  };
  const r = visit(tree, 0, 'root');
  const ce = R ? -R * Math.log(1 - r.eu) : r.ev;
  return { expectedValue: r.ev, certaintyEquivalent: ce, riskPremium: r.ev - ce, policy: r.policy };
}

/**
 * Value of flexibility: rollback value of a tree that keeps downstream decisions, less the value of the rigid
 * tree in which those decisions are committed up front.
 * Ref: Trigeorgis (1996), Real Options, ch. 4 (expanded NPV); Clemen & Reilly ch. 12 (value of control).
 * @param {{flexible: object, rigid: object, riskTolerance?: number}} o two trees in rollbackTree format
 * @returns {{flexible: number, rigid: number, valueOfFlexibility: number}} certainty-equivalent values (currency)
 */
export function valueOfFlexibility({ flexible, rigid, riskTolerance }) {
  const a = rollbackTree(flexible, { riskTolerance }).certaintyEquivalent;
  const b = rollbackTree(rigid, { riskTolerance }).certaintyEquivalent;
  return { flexible: a, rigid: b, valueOfFlexibility: a - b };
}

// ---------------------------------------------------------------------------------------------------------------
// Comparison table and recommendations (deterministic text templates)
// ---------------------------------------------------------------------------------------------------------------

function fmt(x, digits = 2) {
  if (x === null || x === undefined || !Number.isFinite(x)) return 'n/a';
  const a = Math.abs(x);
  let s;
  if (a >= 1e9) s = `${(x / 1e9).toFixed(2)} bn`;
  else if (a >= 1e6) s = `${(x / 1e6).toFixed(2)} M`;
  else if (a >= 1e4) s = `${(x / 1e3).toFixed(1)} k`;
  else s = x.toFixed(digits);
  return s;
}
const pct = (x, d = 1) => (Number.isFinite(x) ? `${(100 * x).toFixed(d)} %` : 'n/a');

/**
 * Alternative-comparison table combining technical KPIs, cost, schedule, risk and carbon. Each row is one
 * measure with its sense; the best alternative per row is flagged and the difference to the first (base)
 * alternative is given. The result follows the app's table shape (columns/rows).
 * Ref: comparison of alternatives, Blank & Tarquin, Engineering Economy, ch. 10 (multiple attributes).
 * @param {{alternatives: {name: string, values: Record<string, number>}[], measures: {key: string, label: string,
 *   unit?: string, sense?: 'max' | 'min', group?: string}[]}} o values keyed by measure key (bus values, cost,
 *   schedule, risk and carbon figures supplied by the caller)
 * @returns {{title: string, columns: object[], rows: object[], wins: Record<string, number>}} wins = number of
 *   measures on which each alternative is best
 */
export function comparisonTable({ alternatives, measures }) {
  if (!Array.isArray(alternatives) || alternatives.length === 0) throw new RangeError('alternatives must be a non-empty array');
  if (!Array.isArray(measures) || measures.length === 0) throw new RangeError('measures must be a non-empty array');
  const columns = [{ key: 'measure', label: 'Measure' }, { key: 'unit', label: 'Unit' },
    ...alternatives.map((a, i) => ({ key: `a${i}`, label: a.name })), { key: 'best', label: 'Best' }];
  const wins = {}; for (const a of alternatives) wins[a.name] = 0;
  const rows = measures.map((m) => {
    const row = { measure: m.label, unit: m.unit ?? '-', group: m.group ?? '' };
    let best = -1;
    alternatives.forEach((a, i) => {
      const v = a.values?.[m.key];
      row[`a${i}`] = typeof v === 'number' && Number.isFinite(v) ? v : null;
      if (row[`a${i}`] === null || !m.sense) return;
      if (best < 0 || (m.sense === 'max' ? v > row[`a${best}`] : v < row[`a${best}`])) best = i;
    });
    row.best = best >= 0 ? alternatives[best].name : '';
    if (best >= 0) wins[alternatives[best].name]++;
    row.deltaToBase = alternatives.map((_, i) => (row[`a${i}`] !== null && row.a0 !== null ? row[`a${i}`] - row.a0 : null));
    return row;
  });
  return { title: 'Comparison of alternatives', columns, rows, wins };
}

/**
 * Rule-based recommendation composer. From engineering check results (utilisation uc against a limit), KPI deltas
 * between a base and a proposed design, economic metrics with their acceptance criteria and an optional ranked
 * list of alternatives, it produces ranked plain-language recommendations with the governing numbers quoted.
 * The text is generated from fixed templates and is fully deterministic: same input, same output.
 * Priority: failed checks (by overshoot) > failed financial criteria > marginal items > improvements > information.
 * Ref: none (rule base documented here); thresholds: uc > limit = fail, uc > 0.9 × limit = marginal,
 * uc < 0.5 × limit = optimisation potential.
 * @param {{checks?: {name: string, uc: number, limit?: number, ref?: string}[], kpiDeltas?: {name: string,
 *   base: number, value: number, unit?: string, sense?: 'max' | 'min'}[], economics?: {npv?: number, irr?: number,
 *   hurdleRate?: number, paybackDiscounted?: number | null, maxPayback?: number, lcoe?: number, lcoeTarget?: number,
 *   probNpvNegative?: number, maxProbLoss?: number, minDscr?: number, covenantDscr?: number, currency?: string},
 *   alternatives?: {name: string, npv?: number, capex?: number, carbon?: number, score?: number}[]}} o
 * @returns {{rank: number, severity: 'fail' | 'warn' | 'ok' | 'info', priority: number, text: string,
 *   action: string}[]} sorted by priority (highest first)
 */
export function composeRecommendations({ checks = [], kpiDeltas = [], economics, alternatives = [] } = {}) {
  const out = [];
  const add = (severity, priority, text, action) => out.push({ rank: 0, severity, priority, text, action });
  for (const c of checks) {
    const limit = c.limit ?? 1;
    requireNumber(`uc of '${c.name}'`, c.uc, { min: 0 }); requireNumber(`limit of '${c.name}'`, limit, { gt: 0 });
    const ratio = c.uc / limit; const ref = c.ref ? ` (${c.ref})` : '';
    if (ratio > 1) add('fail', 1000 + 100 * (ratio - 1), `${c.name} is not satisfied: utilisation ${c.uc.toFixed(2)} against a limit of ${limit.toFixed(2)}, ${pct(ratio - 1, 0)} over${ref}.`,
      `Increase capacity or reduce demand by at least ${pct(1 - 1 / ratio, 0)} before the design is priced.`);
    else if (ratio > 0.9) add('warn', 500 + 100 * (ratio - 0.9), `${c.name} is marginal: utilisation ${c.uc.toFixed(2)} leaves ${pct(1 - ratio, 0)} reserve${ref}.`,
      'Confirm the governing load case and input uncertainty; hold contingency for strengthening.');
    else if (ratio < 0.5) add('info', 100 + 100 * (0.5 - ratio), `${c.name} is lightly utilised (${c.uc.toFixed(2)} of ${limit.toFixed(2)}).`,
      `Up to ${pct(1 - ratio, 0)} of the capacity is unused; test a lighter or cheaper design in the optimisation study.`);
  }
  if (economics) {
    const e = economics; const cur = e.currency ? ` ${e.currency}` : '';
    if (typeof e.npv === 'number') {
      if (e.npv < 0) add('fail', 900, `Net present value is negative (${fmt(e.npv)}${cur}).`,
        'The project destroys value at the stated discount rate; review capex, production and price assumptions or reject.');
      else add('ok', 50, `Net present value is positive (${fmt(e.npv)}${cur}).`, 'Economically acceptable on NPV at the stated discount rate.');
    }
    if (typeof e.irr === 'number' && typeof e.hurdleRate === 'number') {
      const d = e.irr - e.hurdleRate;
      if (d < 0) add('fail', 880 + 100 * -d, `Internal rate of return ${pct(e.irr)} is below the hurdle rate ${pct(e.hurdleRate)} by ${(100 * -d).toFixed(1)} points.`,
        'Improve returns or do not sanction; check the break-even price and capex for the size of the gap.');
      else if (d < 0.02) add('warn', 450, `Internal rate of return ${pct(e.irr)} clears the hurdle rate ${pct(e.hurdleRate)} by only ${(100 * d).toFixed(1)} points.`,
        'Returns are sensitive to small changes; run the sensitivity and Monte Carlo studies before sanction.');
      else add('ok', 40, `Internal rate of return ${pct(e.irr)} exceeds the hurdle rate ${pct(e.hurdleRate)} by ${(100 * d).toFixed(1)} points.`, 'No action.');
    }
    if (typeof e.maxPayback === 'number') {
      if (e.paybackDiscounted === null || e.paybackDiscounted === undefined) add('fail', 860, 'The investment is not recovered within the evaluation period on a discounted basis.',
        `A discounted payback within ${e.maxPayback.toFixed(1)} years is required; the case does not meet it.`);
      else if (e.paybackDiscounted > e.maxPayback) add('fail', 850, `Discounted payback of ${e.paybackDiscounted.toFixed(1)} years exceeds the maximum of ${e.maxPayback.toFixed(1)} years.`,
        'Shorten construction, phase capex or improve early cash flow.');
      else add('ok', 30, `Discounted payback of ${e.paybackDiscounted.toFixed(1)} years is within the maximum of ${e.maxPayback.toFixed(1)} years.`, 'No action.');
    }
    if (typeof e.lcoe === 'number' && typeof e.lcoeTarget === 'number') {
      const r = e.lcoe / e.lcoeTarget - 1;
      if (r > 0) add('warn', 600 + 100 * r, `Levelised cost ${e.lcoe.toFixed(1)} is ${pct(r, 0)} above the target of ${e.lcoeTarget.toFixed(1)} per unit.`,
        'Target the largest LCOE components first (see the component breakdown) and the capacity-factor sensitivity.');
      else add('ok', 20, `Levelised cost ${e.lcoe.toFixed(1)} is ${pct(-r, 0)} below the target of ${e.lcoeTarget.toFixed(1)} per unit.`, 'No action.');
    }
    if (typeof e.probNpvNegative === 'number') {
      const lim = e.maxProbLoss ?? 0.2;
      if (e.probNpvNegative > lim) add('warn', 700 + 100 * (e.probNpvNegative - lim), `Probability of a negative NPV is ${pct(e.probNpvNegative, 0)}, above the tolerance of ${pct(lim, 0)}.`,
        'Reduce exposure on the top-ranked uncertainty drivers (contract price, capex fix, schedule float).');
      else add('ok', 25, `Probability of a negative NPV is ${pct(e.probNpvNegative, 0)}, within the tolerance of ${pct(lim, 0)}.`, 'No action.');
    }
    if (typeof e.minDscr === 'number' && typeof e.covenantDscr === 'number') {
      if (e.minDscr < e.covenantDscr) add('fail', 870, `Minimum debt-service cover ${e.minDscr.toFixed(2)} breaches the covenant of ${e.covenantDscr.toFixed(2)}.`,
        'Reduce gearing, lengthen the tenor or sculpt repayments to cash flow.');
      else add('ok', 20, `Minimum debt-service cover ${e.minDscr.toFixed(2)} meets the covenant of ${e.covenantDscr.toFixed(2)}.`, 'No action.');
    }
  }
  for (const k of kpiDeltas) {
    requireNumber(`base of '${k.name}'`, k.base); requireNumber(`value of '${k.name}'`, k.value);
    const d = k.value - k.base; const rel = k.base !== 0 ? d / Math.abs(k.base) : NaN;
    if (d === 0) continue;
    const unit = k.unit ? ` ${k.unit}` : '';
    const better = k.sense ? (k.sense === 'max' ? d > 0 : d < 0) : null;
    const change = `${k.name} changes from ${fmt(k.base)} to ${fmt(k.value)}${unit} (${d > 0 ? '+' : '−'}${Number.isFinite(rel) ? pct(Math.abs(rel), 1) : fmt(Math.abs(d))})`;
    if (better === true) add('ok', 200 + Math.min(100, 100 * Math.abs(rel || 0)), `${change}, an improvement.`, 'Retain the change.');
    else if (better === false) add('warn', 300 + Math.min(100, 100 * Math.abs(rel || 0)), `${change}, a deterioration.`, 'Justify the change against its benefits elsewhere or revert.');
    else add('info', 60, `${change}.`, 'For information.');
  }
  if (alternatives.length > 1) {
    const key = alternatives.every((a) => typeof a.score === 'number') ? 'score' : alternatives.every((a) => typeof a.npv === 'number') ? 'npv' : null;
    if (key) {
      const sorted = alternatives.slice().sort((a, b) => b[key] - a[key] || (a.name < b.name ? -1 : 1));
      const lead = sorted[0][key] - sorted[1][key];
      let text = `${sorted[0].name} ranks first on ${key === 'npv' ? 'NPV' : 'overall score'} (${fmt(sorted[0][key])}), ahead of ${sorted[1].name} by ${fmt(lead)}.`;
      if (typeof sorted[0].carbon === 'number' && typeof sorted[1].carbon === 'number') {
        const dc = sorted[0].carbon - sorted[1].carbon;
        text += dc <= 0 ? ` It also has the lower life-cycle emissions (${fmt(-dc)} t CO2e less).`
          : ` It carries ${fmt(dc)} t CO2e more life-cycle emissions; the implied abatement cost of choosing ${sorted[1].name} is ${key === 'npv' ? `${fmt(lead / dc)} per t CO2e` : 'to be evaluated'}.`;
      }
      add('info', 400, text, `Select ${sorted[0].name} unless non-monetised criteria reverse the ranking.`);
    }
  }
  out.sort((a, b) => b.priority - a.priority || (a.text < b.text ? -1 : a.text > b.text ? 1 : 0));
  out.forEach((r, k) => { r.rank = k + 1; });
  return out;
}
