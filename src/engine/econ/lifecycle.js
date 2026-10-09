// Life-cycle economics: life-cycle cost and total cost of ownership, economic service life, replacement analysis,
// maintenance and inspection economics on simple lifetime models (Weibull, exponential, lognormal, renewal and
// NHPP processes, Markov / semi-Markov availability, reliability block diagrams, fault and event trees),
// repair-versus-replace, life extension, decommissioning liability (estimate, ARO, funding, timing) and residual value.
// Time in years unless stated; money in project currency.

import { gamma } from '../math/special.js';
import { solveDense } from '../math/linalg.js';
import { requireNumber, requireSeries, factorAP, factorAF } from './tvm.js';
import { npv } from './cashflow.js';
import { normCdf } from './risk.js';

// ---------------------------------------------------------------------------------------------------------------
// Life-cycle cost
// ---------------------------------------------------------------------------------------------------------------

/**
 * Life-cycle cost by phase: undiscounted and discounted totals of each phase's annual cost series (index 0 = time
 * zero) and the equivalent annual cost over the horizon.
 * Ref: ISO 15663:2021 (petroleum and natural gas industries — life cycle costing); IEC 60300-3-3 (life cycle
 * costing).
 * @param {{phases: {name: string, costs: number[]}[], rate: number}} o costs per year (currency, positive = cost;
 *   negative entries are credits such as residual value), discount rate (fraction)
 * @returns {{total: number, discounted: number, eac: number, years: number, byPhase: {name: string,
 *   undiscounted: number, discounted: number, share: number}[]}} currency; share = fraction of discounted LCC
 */
export function lifeCycleCost({ phases, rate }) {
  if (!Array.isArray(phases) || phases.length === 0) throw new RangeError('phases must be a non-empty array');
  requireNumber('rate', rate, { gt: -1 });
  let total = 0; let discounted = 0; let years = 0;
  const byPhase = phases.map((p) => {
    requireSeries(`costs of '${p.name}'`, p.costs);
    const u = p.costs.reduce((s, v) => s + v, 0); const d = npv(p.costs, rate);
    total += u; discounted += d; years = Math.max(years, p.costs.length - 1);
    return { name: p.name, undiscounted: u, discounted: d, share: 0 };
  });
  for (const b of byPhase) b.share = discounted !== 0 ? b.discounted / discounted : 0;
  return { total, discounted, eac: years > 0 ? discounted * factorAP(rate, years) : discounted, years, byPhase };
}

/**
 * Total cost of ownership of an asset: acquisition plus the present value of operating, maintenance and downtime
 * costs (each escalating at its own rate), disposal cost and residual-value credit at the end of life.
 * Ref: ISO 15663:2021; Ellram, L.M. (1995), Total cost of ownership, Int. J. Physical Distribution & Logistics
 * Management 25(8).
 * @param {{acquisition: number, operating?: number, maintenance?: number, downtime?: number, disposal?: number,
 *   residual?: number, life: number, rate: number, escalation?: number}} o annual costs in currency per year
 *   (year-1 values), life in years, discount rate and common real escalation (fractions)
 * @returns {{tco: number, eac: number, components: Record<string, number>}} present values (currency)
 */
export function totalCostOfOwnership({ acquisition, operating = 0, maintenance = 0, downtime = 0, disposal = 0, residual = 0,
  life, rate, escalation = 0 }) {
  requireNumber('acquisition', acquisition, { min: 0 }); requireNumber('life', life, { min: 1, integer: true });
  requireNumber('rate', rate, { gt: -1 }); requireNumber('escalation', escalation, { gt: -1 });
  let f = 0;
  for (let k = 1; k <= life; k++) f += Math.pow(1 + escalation, k - 1) * Math.pow(1 + rate, -k);
  const end = Math.pow(1 + rate, -life);
  const components = { acquisition, operating: operating * f, maintenance: maintenance * f, downtime: downtime * f,
    disposal: disposal * end, residual: -residual * end };
  const tco = Object.values(components).reduce((s, v) => s + v, 0);
  return { tco, eac: tco * factorAP(rate, life), components };
}

/**
 * Equivalent annual cost comparison of alternatives with different lives:
 * EAC = (P − S)(A/P, i, n) + S i + annual cost; the lowest EAC is preferred under like-for-like replacement.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 6.
 * @param {{name: string, capex: number, annualCost?: number, salvage?: number, life: number}[]} alternatives
 * @param {number} rate discount rate (fraction)
 * @returns {{results: {name: string, eac: number, capitalRecovery: number}[], best: string}} currency per year
 */
export function eacComparison(alternatives, rate) {
  if (!Array.isArray(alternatives) || alternatives.length === 0) throw new RangeError('alternatives must be a non-empty array');
  const results = alternatives.map((a) => {
    const cr = (a.capex - (a.salvage ?? 0)) * factorAP(rate, a.life) + (a.salvage ?? 0) * rate;
    return { name: a.name, eac: cr + (a.annualCost ?? 0), capitalRecovery: cr };
  });
  return { results, best: results.reduce((b, r) => (r.eac < b.eac ? r : b)).name };
}

/**
 * Economic service life: the retention period n that minimises the equivalent annual cost
 * EAC(n) = [P − S_n (P/F, i, n) + Σ_{k<=n} OM_k (P/F, i, k)] (A/P, i, n) with rising O&M and falling salvage.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 11 (economic service life).
 * @param {{capex: number, om: number[], salvage: number[], rate: number}} o first cost (currency), O&M cost in
 *   year k (index k − 1), salvage value at the end of year k (index k − 1), discount rate (fraction)
 * @returns {{economicLife: number, minEac: number, eac: number[], capitalRecovery: number[], omEquivalent: number[]}}
 *   years and currency per year; arrays indexed by n − 1
 */
export function economicServiceLife({ capex, om, salvage, rate }) {
  requireNumber('capex', capex, { gt: 0 }); requireSeries('om', om); requireSeries('salvage', salvage);
  requireNumber('rate', rate, { gt: -1 });
  const N = Math.min(om.length, salvage.length);
  const eac = []; const capitalRecovery = []; const omEquivalent = [];
  let pvOm = 0; let best = 0;
  for (let n = 1; n <= N; n++) {
    pvOm += om[n - 1] * Math.pow(1 + rate, -n);
    const ap = factorAP(rate, n);
    const cr = (capex - salvage[n - 1] * Math.pow(1 + rate, -n)) * ap;
    capitalRecovery.push(cr); omEquivalent.push(pvOm * ap); eac.push(cr + pvOm * ap);
    if (eac[n - 1] < eac[best]) best = n - 1;
  }
  return { economicLife: best + 1, minEac: eac[best], eac, capitalRecovery, omEquivalent };
}

/**
 * Replacement analysis (defender versus challenger): the challenger's minimum EAC is compared with the defender's
 * marginal cost of one more year, MC_k = MV_(k−1)(1 + i) − MV_k + OM_k. Replace as soon as the defender's
 * marginal cost exceeds the challenger's minimum EAC. Sunk costs of the defender are excluded: its capital is
 * its current market value.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 11 (replacement study, marginal-cost approach).
 * @param {{defender: {marketValue: number, om: number[], salvage: number[]}, challenger: {capex: number, om: number[],
 *   salvage: number[]}, rate: number}} o defender market value now and at the end of each further year; O&M per year
 * @returns {{replaceNow: boolean, keepYears: number, defenderMarginalCost: number[], challengerMinEac: number,
 *   challengerLife: number, defenderMinEac: number}} currency per year; keepYears = further years to keep the defender
 */
export function replacementAnalysis({ defender, challenger, rate }) {
  if (!defender || !challenger) throw new RangeError('defender and challenger are required');
  requireNumber('defender.marketValue', defender.marketValue, { min: 0 });
  const ch = economicServiceLife({ capex: challenger.capex, om: challenger.om, salvage: challenger.salvage, rate });
  requireSeries('defender.om', defender.om); requireSeries('defender.salvage', defender.salvage);
  const N = Math.min(defender.om.length, defender.salvage.length);
  const mc = []; let prev = defender.marketValue; let keep = N;
  for (let k = 1; k <= N; k++) {
    mc.push(prev * (1 + rate) - defender.salvage[k - 1] + defender.om[k - 1]);
    prev = defender.salvage[k - 1];
  }
  for (let k = 0; k < N; k++) if (mc[k] > ch.minEac) { keep = k; break; }
  const def = defender.marketValue > 0
    ? economicServiceLife({ capex: defender.marketValue, om: defender.om, salvage: defender.salvage, rate }).minEac : Math.min(...mc);
  return { replaceNow: keep === 0, keepYears: keep, defenderMarginalCost: mc, challengerMinEac: ch.minEac,
    challengerLife: ch.economicLife, defenderMinEac: def };
}

// ---------------------------------------------------------------------------------------------------------------
// Lifetime models
// ---------------------------------------------------------------------------------------------------------------

/**
 * Two-parameter Weibull lifetime: reliability R(t) = exp[−(t/η)^β], hazard h(t) = (β/η)(t/η)^(β−1), unreliability
 * and mean life η Γ(1 + 1/β).
 * Ref: Weibull, W. (1951), J. Applied Mechanics 18; Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 2.
 * @param {{t: number, shape: number, scale: number}} o time, shape β (-), scale η (same time unit)
 * @returns {{reliability: number, unreliability: number, hazard: number, mttf: number}} hazard per time unit
 */
export function weibullLife({ t, shape, scale }) {
  requireNumber('t', t, { min: 0 }); requireNumber('shape', shape, { gt: 0 }); requireNumber('scale', scale, { gt: 0 });
  const R = Math.exp(-Math.pow(t / scale, shape));
  const hazard = t === 0 ? (shape < 1 ? Infinity : shape === 1 ? 1 / scale : 0) : shape / scale * Math.pow(t / scale, shape - 1);
  return { reliability: R, unreliability: 1 - R, hazard, mttf: scale * gamma(1 + 1 / shape) };
}

/**
 * Exponential lifetime (constant failure rate λ): R(t) = e^(−λt), MTTF = 1/λ.
 * Ref: Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 2.
 * @param {{t: number, rate: number}} o time and failure rate (1/time)
 * @returns {{reliability: number, unreliability: number, hazard: number, mttf: number}}
 */
export function exponentialLife({ t, rate }) {
  requireNumber('t', t, { min: 0 }); requireNumber('rate', rate, { gt: 0 });
  const R = Math.exp(-rate * t);
  return { reliability: R, unreliability: 1 - R, hazard: rate, mttf: 1 / rate };
}

/**
 * Lognormal lifetime: ln T ~ N(μ, σ²); R(t) = 1 − Φ[(ln t − μ)/σ], MTTF = exp(μ + σ²/2).
 * Ref: Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 2.
 * @param {{t: number, mu: number, sigma: number}} o time (> 0), log-mean and log-standard deviation
 * @returns {{reliability: number, unreliability: number, hazard: number, mttf: number}}
 */
export function lognormalLife({ t, mu, sigma }) {
  requireNumber('t', t, { gt: 0 }); requireNumber('mu', mu); requireNumber('sigma', sigma, { gt: 0 });
  const z = (Math.log(t) - mu) / sigma; const F = normCdf(z);
  const pdf = Math.exp(-0.5 * z * z) / (t * sigma * Math.sqrt(2 * Math.PI));
  return { reliability: 1 - F, unreliability: F, hazard: 1 - F > 0 ? pdf / (1 - F) : Infinity, mttf: Math.exp(mu + 0.5 * sigma * sigma) };
}

// ∫0^T R(t) dt for a Weibull life by composite Simpson
function weibullMeanUpTo(T, shape, scale, n = 400) {
  const h = T / n; let s = 1 + Math.exp(-Math.pow(T / scale, shape));
  for (let k = 1; k < n; k++) s += (k % 2 ? 4 : 2) * Math.exp(-Math.pow(k * h / scale, shape));
  return s * h / 3;
}

/**
 * Age-replacement policy with Weibull life: replace preventively at age T (cost cp) or on failure (cost cf > cp).
 * Long-run cost rate C(T) = [cp R(T) + cf F(T)] / ∫0^T R(t) dt; the optimum exists for shape > 1.
 * Ref: Barlow, R.E. & Hunter, L. (1960), Optimum preventive maintenance policies, Operations Research 8(1);
 * Barlow & Proschan (1965), Mathematical Theory of Reliability, ch. 4.
 * @param {{shape: number, scale: number, costPreventive: number, costFailure: number, maxAge?: number}} o Weibull
 *   β and η (η in the time unit of the result), costs per event (currency)
 * @returns {{optimalAge: number | null, costRate: number, runToFailureCostRate: number, saving: number,
 *   curve: {age: number[], costRate: number[]}}} cost per unit time; optimalAge null = run to failure is optimal
 */
export function ageReplacement({ shape, scale, costPreventive, costFailure, maxAge }) {
  requireNumber('shape', shape, { gt: 0 }); requireNumber('scale', scale, { gt: 0 });
  requireNumber('costPreventive', costPreventive, { gt: 0 }); requireNumber('costFailure', costFailure, { gt: 0 });
  const mttf = scale * gamma(1 + 1 / shape);
  const rtf = costFailure / mttf;
  const C = (T) => { const R = Math.exp(-Math.pow(T / scale, shape)); return (costPreventive * R + costFailure * (1 - R)) / weibullMeanUpTo(T, shape, scale); };
  const Tmax = maxAge ?? 3 * scale;
  const age = []; const costRate = []; let kBest = 0;
  const n = 120;
  for (let k = 1; k <= n; k++) { const T = Tmax * k / n; age.push(T); costRate.push(C(T)); if (costRate[k - 1] < costRate[kBest]) kBest = k - 1; }
  if (shape <= 1 || costFailure <= costPreventive || kBest === n - 1 || costRate[kBest] >= rtf) {
    return { optimalAge: null, costRate: rtf, runToFailureCostRate: rtf, saving: 0, curve: { age, costRate } };
  }
  let a = age[Math.max(0, kBest - 1)]; let b = age[Math.min(n - 1, kBest + 1)];
  const gr = 0.6180339887498949;
  let x1 = b - gr * (b - a); let x2 = a + gr * (b - a); let f1 = C(x1); let f2 = C(x2);
  for (let it = 0; it < 80; it++) {
    if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = C(x1); } else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = C(x2); }
  }
  const T = 0.5 * (a + b); const c = C(T);
  return { optimalAge: T, costRate: c, runToFailureCostRate: rtf, saving: rtf - c, curve: { age, costRate } };
}

/**
 * Renewal function M(t) (expected number of failures with perfect renewal at each failure) for a Weibull life,
 * from the discretised renewal equation M(t) = F(t) + ∫0^t M(t − x) dF(x) (midpoint rule).
 * Ref: Xie, M. (1989), On the solution of renewal-type integral equations, Communications in Statistics B18(1);
 * Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 7.
 * @param {{t: number, shape: number, scale: number, steps?: number}} o time horizon, Weibull β and η, grid steps
 * @returns {{expectedFailures: number, grid: {t: number[], M: number[]}}}
 */
export function renewalFunction({ t, shape, scale, steps = 400 }) {
  requireNumber('t', t, { gt: 0 }); requireNumber('shape', shape, { gt: 0 }); requireNumber('scale', scale, { gt: 0 });
  requireNumber('steps', steps, { min: 10, integer: true });
  const h = t / steps;
  const cdf = (x) => 1 - Math.exp(-Math.pow(x / scale, shape));
  // Riemann–Stieltjes scheme: M_i = F_i + Σ_j F(t_i − t_(j−1/2)) (M_j − M_(j−1)), solved for M_i
  const Fmid = new Float64Array(steps + 1);
  for (let k = 1; k <= steps; k++) Fmid[k] = cdf((k - 0.5) * h);
  const M = new Float64Array(steps + 1);
  for (let i = 1; i <= steps; i++) {
    let s = cdf(i * h) - Fmid[1] * M[i - 1];
    for (let j = 1; j < i; j++) s += Fmid[i - j + 1] * (M[j] - M[j - 1]);
    M[i] = s / (1 - Fmid[1]);
  }
  return { expectedFailures: M[steps], grid: { t: Array.from(M, (_, k) => k * h), M: Array.from(M) } };
}

/**
 * Block-replacement policy: replace at fixed intervals T regardless of age, with renewal at failures in between.
 * Cost rate C(T) = [cp + cf M(T)] / T with M the renewal function.
 * Ref: Barlow & Proschan (1965), Mathematical Theory of Reliability, ch. 4.
 * @param {{interval: number, shape: number, scale: number, costPreventive: number, costFailure: number}} o
 * @returns {{costRate: number, expectedFailures: number}} cost per unit time
 */
export function blockReplacement({ interval, shape, scale, costPreventive, costFailure }) {
  requireNumber('costPreventive', costPreventive, { min: 0 }); requireNumber('costFailure', costFailure, { min: 0 });
  const m = renewalFunction({ t: interval, shape, scale }).expectedFailures;
  return { costRate: (costPreventive + costFailure * m) / interval, expectedFailures: m };
}

/**
 * Non-homogeneous Poisson process with power-law intensity (Crow–AMSAA, minimal repair): expected failures
 * N(t) = λ t^β, intensity u(t) = λ β t^(β−1); with minimal repair cost cm and replacement cost cr the optimal
 * replacement age is T* = [cr / (cm λ (β − 1))]^(1/β) for β > 1.
 * Ref: Crow, L.H. (1974), Reliability analysis for complex repairable systems; Barlow & Hunter (1960) (minimal
 * repair policy).
 * @param {{lambda: number, beta: number, t: number, minimalRepairCost?: number, replacementCost?: number}} o scale
 *   λ (failures per time^β), shape β, time; costs in currency
 * @returns {{expectedFailures: number, intensity: number, expectedRepairCost: number, optimalReplacementAge: number | null,
 *   minCostRate: number | null}}
 */
export function nhppPowerLaw({ lambda, beta, t, minimalRepairCost = 0, replacementCost }) {
  requireNumber('lambda', lambda, { gt: 0 }); requireNumber('beta', beta, { gt: 0 }); requireNumber('t', t, { gt: 0 });
  const N = lambda * Math.pow(t, beta);
  let T = null; let c = null;
  if (replacementCost !== undefined && minimalRepairCost > 0 && beta > 1) {
    T = Math.pow(replacementCost / (minimalRepairCost * lambda * (beta - 1)), 1 / beta);
    c = (replacementCost + minimalRepairCost * lambda * Math.pow(T, beta)) / T;
  }
  return { expectedFailures: N, intensity: lambda * beta * Math.pow(t, beta - 1), expectedRepairCost: minimalRepairCost * N,
    optimalReplacementAge: T, minCostRate: c };
}

/**
 * Steady-state probabilities of a continuous-time Markov chain from its transition-rate matrix (off-diagonal
 * rates q_ij from state i to j; the diagonal is ignored and rebuilt), solving π Q = 0, Σπ = 1; availability is the
 * probability mass of the listed up states.
 * Ref: Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 8 (Markov processes).
 * @param {{rates: number[][], upStates?: number[]}} o rates per unit time; upStates = indices counted as available
 * @returns {{probabilities: number[], availability: number | null}}
 */
export function markovSteadyState({ rates, upStates }) {
  const n = Array.isArray(rates) ? rates.length : 0;
  if (n < 2) throw new RangeError('rates must be an n × n matrix with n >= 2');
  const A = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    if (!rates[i] || rates[i].length !== n) throw new RangeError('rates must be square');
    let out = 0;
    for (let j = 0; j < n; j++) if (j !== i) { requireNumber('rate', rates[i][j], { min: 0 }); out += rates[i][j]; A[j][i] = rates[i][j]; }
    A[i][i] = -out;
  }
  for (let j = 0; j < n; j++) A[n - 1][j] = 1;                 // replace one balance equation by normalisation
  const b = new Array(n).fill(0); b[n - 1] = 1;
  const p = Array.from(solveDense(A, b));
  return { probabilities: p, availability: upStates ? upStates.reduce((s, k) => s + p[k], 0) : null };
}

/**
 * Discrete-time Markov chain evolution p_(k+1) = p_k P (e.g. condition-state deterioration per inspection
 * interval).
 * Ref: Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 8.
 * @param {{P: number[][], p0: number[], steps: number}} o row-stochastic transition matrix, initial distribution
 * @returns {{history: number[][], final: number[]}} history[k] = distribution after k steps
 */
export function markovChainEvolve({ P, p0, steps }) {
  const n = P.length; requireNumber('steps', steps, { min: 0, integer: true });
  for (const row of P) { const s = row.reduce((a, v) => a + requireNumber('P', v, { min: 0, max: 1 }), 0); if (row.length !== n || Math.abs(s - 1) > 1e-9) throw new RangeError('P must be square with rows summing to 1'); }
  let p = Array.from(p0); const history = [p];
  for (let k = 0; k < steps; k++) {
    const q = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) q[j] += p[i] * P[i][j];
    p = q; history.push(p);
  }
  return { history, final: p };
}

/**
 * Semi-Markov steady state: with embedded-chain stationary distribution ν (ν = ν P) and mean sojourn times τ, the
 * long-run fraction of time in state i is ν_i τ_i / Σ ν_j τ_j — sojourn distributions need not be exponential.
 * Ref: Ross, S.M., Introduction to Probability Models, ch. 7 (semi-Markov processes).
 * @param {{P: number[][], meanSojourn: number[], upStates?: number[]}} o embedded transition matrix (rows sum to 1)
 *   and mean time spent in each state per visit
 * @returns {{probabilities: number[], embedded: number[], availability: number | null}}
 */
export function semiMarkovSteadyState({ P, meanSojourn, upStates }) {
  const n = P.length;
  const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => P[j][i] - (i === j ? 1 : 0)));
  for (let j = 0; j < n; j++) A[n - 1][j] = 1;
  const b = new Array(n).fill(0); b[n - 1] = 1;
  const nu = Array.from(solveDense(A, b));
  let s = 0; for (let i = 0; i < n; i++) s += nu[i] * requireNumber('meanSojourn', meanSojourn[i], { min: 0 });
  const p = nu.map((v, i) => v * meanSojourn[i] / s);
  return { probabilities: p, embedded: nu, availability: upStates ? upStates.reduce((a, k) => a + p[k], 0) : null };
}

/**
 * Reliability block diagram with independent blocks. A block is a number (its reliability or availability) or
 * { type: 'series' | 'parallel' | 'koon', k?, blocks: [...] } (k-out-of-n good, blocks need not be identical).
 * Ref: Rausand & Høyland, System Reliability Theory, 2nd ed., ch. 3–4.
 * @param {number | object} block nested block description
 * @returns {number} system reliability (0–1)
 */
export function rbdReliability(block) {
  if (typeof block === 'number') return requireNumber('block reliability', block, { min: 0, max: 1 });
  if (!block || !Array.isArray(block.blocks) || block.blocks.length === 0) throw new RangeError('block needs a non-empty blocks array');
  const r = block.blocks.map(rbdReliability);
  if (block.type === 'series') return r.reduce((p, v) => p * v, 1);
  if (block.type === 'parallel') return 1 - r.reduce((p, v) => p * (1 - v), 1);
  if (block.type === 'koon') {
    const k = requireNumber('k', block.k, { min: 0, max: r.length, integer: true });
    const dp = new Array(r.length + 1).fill(0); dp[0] = 1;       // dp[j] = P(exactly j of the blocks so far work)
    for (const v of r) for (let j = r.length; j >= 0; j--) dp[j] = dp[j] * (1 - v) + (j > 0 ? dp[j - 1] * v : 0);
    let s = 0; for (let j = k; j <= r.length; j++) s += dp[j];
    return s;
  }
  throw new RangeError("block.type must be 'series', 'parallel' or 'koon'");
}

/**
 * Fault-tree top-event probability for independent basic events. A node is a number (basic-event probability) or
 * { gate: 'AND' | 'OR' | 'VOTE', k?, inputs: [...] }. Repeated basic events are treated as independent, so trees
 * with shared events should be restructured (or bounded with minimal cut sets) by the caller.
 * Ref: Vesely et al. (1981), Fault Tree Handbook, NUREG-0492; IEC 61025:2006.
 * @param {number | object} node nested gate description
 * @returns {number} probability of the top event (0–1)
 */
export function faultTreeProbability(node) {
  if (typeof node === 'number') return requireNumber('basic event probability', node, { min: 0, max: 1 });
  if (!node || !Array.isArray(node.inputs) || node.inputs.length === 0) throw new RangeError('gate needs a non-empty inputs array');
  const p = node.inputs.map(faultTreeProbability);
  if (node.gate === 'AND') return p.reduce((a, v) => a * v, 1);
  if (node.gate === 'OR') return 1 - p.reduce((a, v) => a * (1 - v), 1);
  if (node.gate === 'VOTE') return rbdReliability({ type: 'koon', k: node.k, blocks: p });
  throw new RangeError("gate must be 'AND', 'OR' or 'VOTE'");
}

/**
 * Event tree: an initiating event of given frequency followed by barriers that each fail with a probability;
 * returns every sequence with its frequency and, when `cost(failed)` is supplied, the consequence cost and the
 * expected annual loss.
 * Ref: IEC 62502:2010 (event tree analysis); ISO 31010:2019 Annex B.
 * @param {{frequency: number, barriers: {name: string, pFail: number}[], cost?: (failed: boolean[]) => number}} o
 *   frequency per year; cost returns the consequence (currency) of a sequence given which barriers failed
 * @returns {{sequences: {failed: boolean[], label: string, frequency: number, cost: number | null}[],
 *   expectedAnnualLoss: number | null}}
 */
export function eventTree({ frequency, barriers, cost }) {
  requireNumber('frequency', frequency, { min: 0 });
  if (!Array.isArray(barriers) || barriers.length === 0 || barriers.length > 16) throw new RangeError('barriers must have 1 to 16 entries');
  const n = barriers.length; const sequences = []; let eal = cost ? 0 : null;
  for (let mask = 0; mask < 1 << n; mask++) {
    let f = frequency; const failed = [];
    for (let j = 0; j < n; j++) {
      const fl = !!(mask & (1 << j)); failed.push(fl);
      const p = requireNumber(`pFail of '${barriers[j].name}'`, barriers[j].pFail, { min: 0, max: 1 });
      f *= fl ? p : 1 - p;
    }
    const c = cost ? cost(failed) : null;
    if (cost) eal += f * c;
    sequences.push({ failed, label: barriers.map((b, j) => `${b.name}: ${failed[j] ? 'fails' : 'works'}`).join(', '), frequency: f, cost: c });
  }
  return { sequences, expectedAnnualLoss: eal };
}

/**
 * Risk-based inspection economics: for each candidate inspection interval, the present value of inspection cost
 * plus expected failure cost over the horizon. Model (deliberately simple and explicit): the annual failure
 * probability is the conditional Weibull probability at the component's effective age; an inspection with
 * effectiveness e (probability that degradation is found and repaired) reduces the effective age by the factor
 * (1 − e); a failure costs `failureCost`. PoF(t) may instead be supplied as a function of effective age.
 * Ref: API RP 580 / 581 (risk-based inspection) and DNV-RP-G101 for the framework; age-reduction (virtual age)
 * model after Kijima, M. (1989), J. Applied Probability 26(1).
 * @param {{intervals: number[], horizon: number, inspectionCost: number, failureCost: number, effectiveness: number,
 *   rate: number, weibull?: {shape: number, scale: number}, pof?: (age: number) => number, initialAge?: number}} o
 *   intervals and horizon in whole years; pof = cumulative failure probability at an effective age
 * @returns {{results: {interval: number, pvInspection: number, pvRisk: number, pvTotal: number,
 *   expectedFailures: number}[], best: {interval: number, pvTotal: number}, noInspection: {pvRisk: number}}}
 */
export function rbiEconomics({ intervals, horizon, inspectionCost, failureCost, effectiveness, rate, weibull, pof, initialAge = 0 }) {
  requireSeries('intervals', intervals); requireNumber('horizon', horizon, { min: 1, integer: true });
  requireNumber('inspectionCost', inspectionCost, { min: 0 }); requireNumber('failureCost', failureCost, { min: 0 });
  requireNumber('effectiveness', effectiveness, { min: 0, max: 1 }); requireNumber('rate', rate, { gt: -1 });
  const F = pof ?? ((a) => 1 - Math.exp(-Math.pow(a / weibull.scale, weibull.shape)));
  if (!pof && !(weibull && weibull.shape > 0 && weibull.scale > 0)) throw new RangeError('give weibull {shape, scale} or a pof function');
  const run = (interval) => {
    let age = initialAge; let pvI = 0; let pvR = 0; let nf = 0;
    for (let y = 1; y <= horizon; y++) {
      const Fa = F(age); const Fb = F(age + 1);
      const p = Fa < 1 ? (Fb - Fa) / (1 - Fa) : 1;                // conditional probability of failure in the year
      const d = Math.pow(1 + rate, -y);
      pvR += p * failureCost * d; nf += p;
      age += 1;
      if (interval > 0 && y % interval === 0 && y < horizon) { pvI += inspectionCost * d; age *= 1 - effectiveness; }
    }
    return { interval, pvInspection: pvI, pvRisk: pvR, pvTotal: pvI + pvR, expectedFailures: nf };
  };
  const results = Array.from(intervals, (iv) => run(requireNumber('interval', iv, { min: 1, integer: true })));
  const best = results.reduce((b, r) => (r.pvTotal < b.pvTotal ? r : b));
  return { results, best: { interval: best.interval, pvTotal: best.pvTotal }, noInspection: { pvRisk: run(0).pvRisk } };
}

/**
 * Repair versus replace: equivalent annual cost of repairing (repair cost buys `repairLife` more years at the old
 * O&M) against replacing (new asset over its life at the new O&M), and the break-even repair cost.
 * Ref: Blank & Tarquin, Engineering Economy, ch. 11; EAC comparison under repeatability.
 * @param {{repairCost: number, repairLife: number, omOld: number, replaceCost: number, newLife: number, omNew: number,
 *   salvageOld?: number, salvageNew?: number, rate: number}} o costs in currency, lives in years, O&M per year;
 *   salvageOld = value of the old asset given up if repaired rather than sold (opportunity cost)
 * @returns {{eacRepair: number, eacReplace: number, decision: 'repair' | 'replace', breakEvenRepairCost: number}}
 */
export function repairVsReplace({ repairCost, repairLife, omOld, replaceCost, newLife, omNew, salvageOld = 0, salvageNew = 0, rate }) {
  requireNumber('repairCost', repairCost, { min: 0 }); requireNumber('replaceCost', replaceCost, { min: 0 });
  const eacRepair = (repairCost + salvageOld) * factorAP(rate, repairLife) + omOld;
  const eacReplace = (replaceCost - salvageNew) * factorAP(rate, newLife) + salvageNew * rate + omNew;
  return { eacRepair, eacReplace, decision: eacRepair <= eacReplace ? 'repair' : 'replace',
    breakEvenRepairCost: (eacReplace - omOld) / factorAP(rate, repairLife) - salvageOld };
}

/**
 * Life-extension business case including ageing risk: NPV of extending for n years = −upgrade capex + Σ (net
 * revenue − expected failure cost) DF − deferred decommissioning PV, compared with decommissioning now. The annual
 * failure probability is the conditional Weibull probability at the asset's age (an engineering life assessment
 * should supply the Weibull parameters or the remaining-life estimate).
 * Ref: incremental cash flows, Brealey, Myers & Allen ch. 6; ageing-risk treatment as in DNV-RP-G101 / NORSOK
 * N-006 life-extension assessments.
 * @param {{currentAge: number, extensionYears: number, upgradeCapex: number, netRevenue: number | number[],
 *   failureCost?: number, weibull?: {shape: number, scale: number}, decomCost: number, decomEscalation?: number,
 *   rate: number}} o years, currency; netRevenue per extension year; real discount rate
 * @returns {{npvExtend: number, npvDecommissionNow: number, npvIncremental: number, expectedFailureCost: number[],
 *   annualFailureProbability: number[], decision: 'extend' | 'decommission', optimalExtension: number}}
 *   optimalExtension = number of extension years (0…n) with the highest NPV
 */
export function lifeExtensionCase({ currentAge, extensionYears, upgradeCapex, netRevenue, failureCost = 0, weibull, decomCost,
  decomEscalation = 0, rate }) {
  requireNumber('currentAge', currentAge, { min: 0 }); requireNumber('extensionYears', extensionYears, { min: 1, integer: true });
  requireNumber('upgradeCapex', upgradeCapex, { min: 0 }); requireNumber('decomCost', decomCost, { min: 0 });
  requireNumber('rate', rate, { gt: -1 });
  const probs = []; const efc = [];
  for (let y = 1; y <= extensionYears; y++) {
    let p = 0;
    if (weibull) {
      const R0 = Math.exp(-Math.pow((currentAge + y - 1) / weibull.scale, weibull.shape));
      const R1 = Math.exp(-Math.pow((currentAge + y) / weibull.scale, weibull.shape));
      p = R0 > 0 ? 1 - R1 / R0 : 1;
    }
    probs.push(p); efc.push(p * failureCost);
  }
  const rev = (k) => (typeof netRevenue === 'number' ? netRevenue : netRevenue[Math.min(k, netRevenue.length - 1)]);
  const valueOf = (n) => {
    if (n === 0) return -decomCost;
    let v = -upgradeCapex;
    for (let y = 1; y <= n; y++) v += (rev(y - 1) - efc[y - 1]) * Math.pow(1 + rate, -y);
    return v - decomCost * Math.pow(1 + decomEscalation, n) * Math.pow(1 + rate, -n);
  };
  let best = 0; let bestV = valueOf(0);
  for (let n = 1; n <= extensionYears; n++) { const v = valueOf(n); if (v > bestV) { bestV = v; best = n; } }
  const ext = valueOf(extensionYears); const now = valueOf(0);
  return { npvExtend: ext, npvDecommissionNow: now, npvIncremental: ext - now, expectedFailureCost: efc,
    annualFailureProbability: probs, decision: ext >= now ? 'extend' : 'decommission', optimalExtension: best };
}

// ---------------------------------------------------------------------------------------------------------------
// Decommissioning liability and residual value
// ---------------------------------------------------------------------------------------------------------------

/**
 * Decommissioning cost estimate build-up: removal, well plugging and abandonment, subsea recovery, onshore
 * disposal, site restoration and post-decommissioning monitoring, less the recycling credit (recoverable steel
 * tonnage × scrap price), plus contingency on the gross cost.
 * Ref: cost categories of the OEUK (Oil & Gas UK) Decommissioning Cost Estimation guidelines work breakdown;
 * all amounts are caller inputs.
 * @param {{removal?: number, wells?: {count: number, costPerWell: number}, subseaRecovery?: number, disposal?: number,
 *   restoration?: number, monitoring?: {years: number, costPerYear: number}, scrap?: {mass_t: number,
 *   recoverableFraction: number, pricePerT: number}, contingency?: number}} o currency; contingency fraction of gross
 * @returns {{gross: number, contingency: number, recyclingCredit: number, net: number, breakdown: Record<string, number>}}
 */
export function decomEstimate({ removal = 0, wells, subseaRecovery = 0, disposal = 0, restoration = 0, monitoring, scrap, contingency = 0 }) {
  const nn = (name, v) => requireNumber(name, v ?? 0, { min: 0 });
  const b = { removal: nn('removal', removal), wellPA: wells ? nn('wells.count', wells.count) * nn('wells.costPerWell', wells.costPerWell) : 0,
    subseaRecovery: nn('subseaRecovery', subseaRecovery), disposal: nn('disposal', disposal), restoration: nn('restoration', restoration),
    monitoring: monitoring ? nn('monitoring.years', monitoring.years) * nn('monitoring.costPerYear', monitoring.costPerYear) : 0 };
  const gross = Object.values(b).reduce((s, v) => s + v, 0);
  const cont = gross * nn('contingency', contingency);
  const credit = scrap ? nn('scrap.mass_t', scrap.mass_t) * requireNumber('scrap.recoverableFraction', scrap.recoverableFraction, { min: 0, max: 1 }) * nn('scrap.pricePerT', scrap.pricePerT) : 0;
  return { gross, contingency: cont, recyclingCredit: credit, net: gross + cont - credit, breakdown: b };
}

/**
 * Asset retirement obligation: today's cost estimate inflated to the settlement date and discounted back at the
 * (credit-adjusted risk-free or other prescribed) rate, with the year-by-year accretion of the liability.
 * Ref: IAS 37 §45–47 (present value of provisions); US GAAP ASC 410-20 (asset retirement obligations).
 * @param {{costToday: number, inflation: number, years: number, discountRate: number}} o currency, fractions per
 *   year, whole years to settlement
 * @returns {{futureCost: number, presentValue: number, accretion: {year: number, opening: number, expense: number,
 *   closing: number}[]}} currency
 */
export function assetRetirementObligation({ costToday, inflation, years, discountRate }) {
  requireNumber('costToday', costToday, { min: 0 }); requireNumber('inflation', inflation, { gt: -1 });
  requireNumber('years', years, { min: 0, integer: true }); requireNumber('discountRate', discountRate, { gt: -1 });
  const futureCost = costToday * Math.pow(1 + inflation, years);
  const pv = futureCost * Math.pow(1 + discountRate, -years);
  const accretion = []; let bal = pv;
  for (let y = 1; y <= years; y++) { const e = bal * discountRate; accretion.push({ year: y, opening: bal, expense: e, closing: bal + e }); bal += e; }
  return { futureCost, presentValue: pv, accretion };
}

/**
 * Decommissioning funding (escrow / security) schedule to reach a target at the end of n years: level annual
 * contributions (sinking fund at the fund's return) or contributions in proportion to production (unit of
 * production, solved so that the compounded fund equals the target).
 * Ref: sinking-fund factor, Blank & Tarquin ch. 2; unit-of-production abandonment funding as in petroleum
 * fiscal agreements.
 * @param {{target: number, years: number, fundRate: number, method?: 'level' | 'unitOfProduction',
 *   production?: number[]}} o target in the money of the settlement date; fund return (fraction per year)
 * @returns {{contributions: number[], balance: number[], totalContributed: number, interestEarned: number}}
 */
export function decomFunding({ target, years, fundRate, method = 'level', production }) {
  requireNumber('target', target, { min: 0 }); requireNumber('years', years, { min: 1, integer: true });
  requireNumber('fundRate', fundRate, { gt: -1 });
  let contributions;
  if (method === 'level') contributions = new Array(years).fill(target * factorAF(fundRate, years));
  else if (method === 'unitOfProduction') {
    requireSeries('production', production);
    if (production.length < years) throw new RangeError('production must cover every funding year');
    let s = 0; for (let k = 0; k < years; k++) s += production[k] * Math.pow(1 + fundRate, years - 1 - k);
    if (!(s > 0)) throw new RangeError('production must be positive in at least one year');
    contributions = Array.from({ length: years }, (_, k) => target * production[k] / s);
  } else throw new RangeError("method must be 'level' or 'unitOfProduction'");
  const balance = []; let bal = 0; let tot = 0;
  for (let k = 0; k < years; k++) { bal = bal * (1 + fundRate) + contributions[k]; balance.push(bal); tot += contributions[k]; }
  return { contributions, balance, totalContributed: tot, interestEarned: bal - tot };
}

/**
 * Timing optionality of decommissioning: present value of executing in each candidate year, with cost escalation
 * and the annual cost of holding the asset idle (monitoring, insurance, integrity) until then.
 * Ref: deferral valuation by discounted cash flow; Brealey, Myers & Allen ch. 6 (timing decisions).
 * @param {{costToday: number, escalation: number, rate: number, holdingCostPerYear?: number, candidateYears: number[]}} o
 * @returns {{options: {year: number, pvCost: number}[], best: {year: number, pvCost: number}}} PV of total cost
 */
export function decomTiming({ costToday, escalation, rate, holdingCostPerYear = 0, candidateYears }) {
  requireNumber('costToday', costToday, { min: 0 }); requireNumber('escalation', escalation, { gt: -1 });
  requireNumber('rate', rate, { gt: -1 }); requireSeries('candidateYears', candidateYears);
  const options = Array.from(candidateYears, (y) => {
    requireNumber('candidate year', y, { min: 0, integer: true });
    let pv = costToday * Math.pow((1 + escalation) / (1 + rate), y);
    for (let k = 1; k <= y; k++) pv += holdingCostPerYear * Math.pow(1 + rate, -k);
    return { year: y, pvCost: pv };
  });
  return { options, best: options.reduce((b, o) => (o.pvCost < b.pvCost ? o : b)) };
}

/**
 * Residual value of an asset at a given age: 'linear' (straight line to the end-of-life value), 'declining'
 * (constant annual decline rate) or 'scrap' (recoverable mass × scrap price less removal cost).
 * Ref: residual-value conventions, IAS 16 §6 and §53; Blank & Tarquin ch. 16.
 * @param {{method: 'linear' | 'declining' | 'scrap', cost?: number, life?: number, age?: number, endValue?: number,
 *   declineRate?: number, mass_t?: number, scrapPricePerT?: number, recoverableFraction?: number,
 *   removalCost?: number}} o currency, years, fractions, tonnes
 * @returns {number} residual value (currency; may be negative for 'scrap' when removal exceeds the credit)
 */
export function residualValue({ method, cost, life, age, endValue = 0, declineRate, mass_t, scrapPricePerT, recoverableFraction = 1, removalCost = 0 }) {
  if (method === 'linear') {
    requireNumber('cost', cost, { min: 0 }); requireNumber('life', life, { gt: 0 }); requireNumber('age', age, { min: 0 });
    return age >= life ? endValue : cost - (cost - endValue) * age / life;
  }
  if (method === 'declining') {
    requireNumber('cost', cost, { min: 0 }); requireNumber('declineRate', declineRate, { min: 0, max: 1 }); requireNumber('age', age, { min: 0 });
    return cost * Math.pow(1 - declineRate, age);
  }
  if (method === 'scrap') {
    return requireNumber('mass_t', mass_t, { min: 0 }) * requireNumber('recoverableFraction', recoverableFraction, { min: 0, max: 1 })
      * requireNumber('scrapPricePerT', scrapPricePerT, { min: 0 }) - requireNumber('removalCost', removalCost, { min: 0 });
  }
  throw new RangeError("method must be 'linear', 'declining' or 'scrap'");
}
