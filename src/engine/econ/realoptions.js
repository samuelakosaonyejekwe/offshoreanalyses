// Real options: Black–Scholes–Merton with dividend yield and Greeks, CRR binomial lattice (European/American),
// real-option templates on a project-value lattice (defer, expand, contract, abandon, switch, sequential compound),
// Margrabe exchange option, project volatility from simulation, expanded NPV and Longstaff–Schwartz Monte Carlo.
// Rates are continuously compounded per year, volatility per √year, time in years.

import { makeRng } from '../math/rng.js';
import { solveDense } from '../math/linalg.js';
import { requireNumber } from './tvm.js';
import { normCdf, normPdf } from './risk.js';
import { rollbackTree } from './decision.js';

const common = (S, K, r, sigma, T, q) => {
  requireNumber('S', S, { gt: 0 }); requireNumber('K', K, { gt: 0 }); requireNumber('r', r);
  requireNumber('sigma', sigma, { gt: 0 }); requireNumber('T', T, { gt: 0 }); requireNumber('q', q);
};

/**
 * Black–Scholes–Merton European call and put on an asset paying a continuous yield q (for real options, q is the
 * value leakage or cost of delay), with Greeks.
 * Assumptions: lognormal value, constant volatility and rates, European exercise, frictionless replication —
 * use only where these are defensible for the underlying project value.
 * Ref: Black, F. & Scholes, M. (1973), J. Political Economy 81(3); Merton, R.C. (1973), Bell J. Economics 4(1);
 * Hull, Options, Futures and Other Derivatives, ch. 15, 17 and 19.
 * @param {{S: number, K: number, r: number, sigma: number, T: number, q?: number}} o underlying value and strike
 *   (currency), risk-free rate and yield (1/yr, continuous), volatility (1/√yr), maturity (yr)
 * @returns {{call: number, put: number, d1: number, d2: number, greeks: {call: object, put: object}}} Greeks per
 *   unit: delta (−), gamma (1/currency), vega (currency per unit volatility), theta (currency/yr), rho (currency
 *   per unit rate)
 */
export function blackScholes({ S, K, r, sigma, T, q = 0 }) {
  common(S, K, r, sigma, T, q);
  const sq = sigma * Math.sqrt(T);
  const d1 = (Math.log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / sq; const d2 = d1 - sq;
  const eq = Math.exp(-q * T); const er = Math.exp(-r * T);
  const call = S * eq * normCdf(d1) - K * er * normCdf(d2);
  const put = K * er * normCdf(-d2) - S * eq * normCdf(-d1);
  const gamma = eq * normPdf(d1) / (S * sq); const vega = S * eq * normPdf(d1) * Math.sqrt(T);
  const decay = -S * eq * normPdf(d1) * sigma / (2 * Math.sqrt(T));
  return { call, put, d1, d2, greeks: {
    call: { delta: eq * normCdf(d1), gamma, vega, theta: decay - r * K * er * normCdf(d2) + q * S * eq * normCdf(d1), rho: K * T * er * normCdf(d2) },
    put: { delta: -eq * normCdf(-d1), gamma, vega, theta: decay + r * K * er * normCdf(-d2) - q * S * eq * normCdf(-d1), rho: -K * T * er * normCdf(-d2) } } };
}

/**
 * Cox–Ross–Rubinstein binomial lattice for a European or American call/put: u = e^(σ√Δt), d = 1/u,
 * p = (e^((r−q)Δt) − d)/(u − d).
 * Ref: Cox, J.C., Ross, S.A. & Rubinstein, M. (1979), Option pricing: a simplified approach, J. Financial
 * Economics 7(3).
 * @param {{S: number, K: number, r: number, sigma: number, T: number, q?: number, steps?: number,
 *   type?: 'call' | 'put', american?: boolean}} o as for blackScholes; steps = number of time steps (default 200)
 * @returns {{value: number, delta: number, earlyExercisePremium: number | null, steps: number}} currency
 */
export function binomialCRR({ S, K, r, sigma, T, q = 0, steps = 200, type = 'call', american = false }) {
  common(S, K, r, sigma, T, q);
  if (type !== 'call' && type !== 'put') throw new RangeError("type must be 'call' or 'put'");
  const sign = type === 'call' ? 1 : -1;
  const res = valueLattice({ V0: S, sigma, r, T, q, steps,
    terminal: (V) => Math.max(sign * (V - K), 0),
    exercise: american ? (V, cont) => Math.max(cont, sign * (V - K)) : undefined });
  let prem = null;
  if (american) {
    prem = res.value - valueLattice({ V0: S, sigma, r, T, q, steps, terminal: (V) => Math.max(sign * (V - K), 0) }).value;
  }
  return { value: res.value, delta: res.delta, earlyExercisePremium: prem, steps };
}

/**
 * General recombining CRR lattice valuation of a claim on a value process V: terminal payoff `terminal(V)` and an
 * optional decision rule `exercise(V, continuation, time, step)` applied at every earlier node (return the node
 * value). This is the engine behind the real-option templates.
 * Ref: Cox, Ross & Rubinstein (1979); Trigeorgis, L. (1996), Real Options, ch. 3 and 5; Copeland & Antikarov
 * (2001), Real Options: A Practitioner's Guide, ch. 5.
 * @param {{V0: number, sigma: number, r: number, T: number, q?: number, steps?: number,
 *   terminal: (V: number) => number, exercise?: (V: number, cont: number, t: number, step: number) => number}} o
 *   V0 = present value of the underlying project (currency), q = value leakage (1/yr)
 * @returns {{value: number, delta: number, u: number, d: number, p: number}} value in currency; p = risk-neutral
 *   up probability
 */
export function valueLattice({ V0, sigma, r, T, q = 0, steps = 200, terminal, exercise }) {
  requireNumber('V0', V0, { gt: 0 }); requireNumber('sigma', sigma, { gt: 0 }); requireNumber('T', T, { gt: 0 });
  requireNumber('r', r); requireNumber('q', q); requireNumber('steps', steps, { min: 1, max: 20000, integer: true });
  if (typeof terminal !== 'function') throw new RangeError('terminal must be a function');
  const dt = T / steps; const u = Math.exp(sigma * Math.sqrt(dt)); const d = 1 / u;
  const p = (Math.exp((r - q) * dt) - d) / (u - d);
  if (!(p > 0 && p < 1)) throw new RangeError('risk-neutral probability outside (0, 1): increase steps or check r, q, sigma');
  const disc = Math.exp(-r * dt);
  const val = new Float64Array(steps + 1);
  for (let j = 0; j <= steps; j++) val[j] = terminal(V0 * Math.pow(u, 2 * j - steps));
  let delta = NaN;
  for (let i = steps - 1; i >= 0; i--) {
    for (let j = 0; j <= i; j++) {
      const cont = disc * (p * val[j + 1] + (1 - p) * val[j]);
      val[j] = exercise ? exercise(V0 * Math.pow(u, 2 * j - i), cont, i * dt, i) : cont;
    }
    if (i === 1) delta = (val[1] - val[0]) / (V0 * (u - d));
  }
  return { value: val[0], delta, u, d, p };
}

/**
 * Project value with embedded real options, valued on one lattice so interacting options are not double counted.
 * Options: { type: 'expand', factor, cost } → factor·V − cost; { type: 'contract', factor, saving } → factor·V +
 * saving; { type: 'abandon', salvage } → salvage; { type: 'switch', value, cost?, growth? } → alternative-use value
 * growing at `growth` less the switching cost. Each may be limited to a window [from, to] in years; `american`
 * false restricts exercise to maturity.
 * Option value = value with options − value of the same lattice without options.
 * Ref: Trigeorgis (1996), Real Options, ch. 5 and 7 (option interactions); Copeland & Antikarov (2001), ch. 5–6.
 * @param {{V0: number, sigma: number, r: number, T: number, q?: number, steps?: number, american?: boolean,
 *   options: object[]}} o V0 = PV of project operating cash flows (currency); option amounts in currency
 * @returns {{valueWithOptions: number, valueWithoutOptions: number, optionValue: number}} currency
 */
export function realOptionLattice({ V0, sigma, r, T, q = 0, steps = 200, american = true, options }) {
  if (!Array.isArray(options) || options.length === 0) throw new RangeError('options must be a non-empty array');
  const pay = options.map((o) => {
    const from = o.from ?? 0; const to = o.to ?? T;
    let f;
    if (o.type === 'expand') { requireNumber('factor', o.factor, { gt: 1 }); requireNumber('cost', o.cost, { min: 0 }); f = (V) => o.factor * V - o.cost; }
    else if (o.type === 'contract') { requireNumber('factor', o.factor, { gt: 0, lt: 1 }); requireNumber('saving', o.saving, { min: 0 }); f = (V) => o.factor * V + o.saving; }
    else if (o.type === 'abandon') { requireNumber('salvage', o.salvage, { min: 0 }); f = () => o.salvage; }
    else if (o.type === 'switch') { requireNumber('value', o.value, { min: 0 }); f = (V, t) => o.value * Math.exp((o.growth ?? 0) * t) - (o.cost ?? 0); }
    else throw new RangeError(`option type '${o.type}' is not supported`);
    return (V, t) => (t >= from - 1e-12 && t <= to + 1e-12 ? f(V, t) : -Infinity);
  });
  const best = (V, t, hold) => { let m = hold; for (const f of pay) { const v = f(V, t); if (v > m) m = v; } return m; };
  const withOpt = valueLattice({ V0, sigma, r, T, q, steps, terminal: (V) => best(V, T, V),
    exercise: american ? (V, cont, t) => best(V, t, cont) : undefined }).value;
  const without = valueLattice({ V0, sigma, r, T, q, steps, terminal: (V) => V }).value;
  return { valueWithOptions: withOpt, valueWithoutOptions: without, optionValue: withOpt - without };
}

/**
 * Option to defer (wait): an American call on the project value V with the investment cost as strike. `q` is the
 * cost of delay (cash flows forgone while waiting); with q = 0 early exercise is never optimal.
 * Ref: McDonald, R. & Siegel, D. (1986), The value of waiting to invest, Quarterly J. Economics 101(4);
 * Dixit & Pindyck (1994), Investment under Uncertainty, ch. 5.
 * @param {{V: number, I: number, r: number, sigma: number, T: number, q?: number, steps?: number}} o V = PV of
 *   operating cash flows, I = investment (currency), T = time the opportunity stays open (yr)
 * @returns {{optionValue: number, staticNpv: number, valueOfWaiting: number, investNow: boolean}} currency;
 *   valueOfWaiting = option value − max(static NPV, 0)
 */
export function optionToDefer({ V, I, r, sigma, T, q = 0, steps = 200 }) {
  const v = binomialCRR({ S: V, K: I, r, sigma, T, q, steps, type: 'call', american: true }).value;
  const st = V - I;
  return { optionValue: v, staticNpv: st, valueOfWaiting: v - Math.max(st, 0), investNow: st > 0 && v - st < 1e-9 * Math.max(1, V) };
}

/**
 * Option to expand the project scale by `factor` at a cost.
 * Ref: Trigeorgis (1996), Real Options, ch. 5.
 * @param {{V: number, factor: number, cost: number, r: number, sigma: number, T: number, q?: number, steps?: number,
 *   american?: boolean}} o V = PV of the base project (currency), factor > 1, cost (currency)
 * @returns {{optionValue: number, valueWithOptions: number, valueWithoutOptions: number}} currency
 */
export function optionToExpand({ V, factor, cost, r, sigma, T, q = 0, steps = 200, american = true }) {
  return realOptionLattice({ V0: V, sigma, r, T, q, steps, american, options: [{ type: 'expand', factor, cost }] });
}

/**
 * Option to contract the project scale to `factor` (< 1) and save `saving`.
 * Ref: Trigeorgis (1996), Real Options, ch. 5.
 * @param {{V: number, factor: number, saving: number, r: number, sigma: number, T: number, q?: number,
 *   steps?: number, american?: boolean}} o currency
 * @returns {{optionValue: number, valueWithOptions: number, valueWithoutOptions: number}} currency
 */
export function optionToContract({ V, factor, saving, r, sigma, T, q = 0, steps = 200, american = true }) {
  return realOptionLattice({ V0: V, sigma, r, T, q, steps, american, options: [{ type: 'contract', factor, saving }] });
}

/**
 * Option to abandon the project for its salvage value (an American put on the project value).
 * Ref: Myers, S.C. & Majd, S. (1990), Abandonment value and project life, Advances in Futures and Options
 * Research 4; Trigeorgis (1996), ch. 5.
 * @param {{V: number, salvage: number, r: number, sigma: number, T: number, q?: number, steps?: number,
 *   american?: boolean}} o currency
 * @returns {{optionValue: number, valueWithOptions: number, valueWithoutOptions: number}} currency
 */
export function optionToAbandon({ V, salvage, r, sigma, T, q = 0, steps = 200, american = true }) {
  return realOptionLattice({ V0: V, sigma, r, T, q, steps, american, options: [{ type: 'abandon', salvage }] });
}

/**
 * Option to switch the asset to an alternative use worth `value` (growing at `growth`) at a switching cost.
 * Ref: Trigeorgis (1996), Real Options, ch. 5 (switching options).
 * @param {{V: number, value: number, cost?: number, growth?: number, r: number, sigma: number, T: number, q?: number,
 *   steps?: number}} o currency; growth 1/yr (continuous)
 * @returns {{optionValue: number, valueWithOptions: number, valueWithoutOptions: number}} currency
 */
export function optionToSwitch({ V, value, cost = 0, growth = 0, r, sigma, T, q = 0, steps = 200 }) {
  return realOptionLattice({ V0: V, sigma, r, T, q, steps, options: [{ type: 'switch', value, cost, growth }] });
}

/**
 * Staged (sequential compound) investment: at each phase gate the next investment is made only if the
 * continuation value exceeds it; the last gate buys the operating project worth V. Valued on a lattice of the
 * project value with the gates mapped to the nearest time steps.
 * Ref: Geske, R. (1979), The valuation of compound options, J. Financial Economics 7(1); Copeland & Antikarov
 * (2001), ch. 6 (sequential compound options on a lattice).
 * @param {{V: number, phases: {time: number, cost: number}[], r: number, sigma: number, q?: number, steps?: number}} o
 *   V = PV today of the completed project's cash flows (currency); phases in increasing time (yr) with their
 *   investment (currency); the first phase may be at time 0
 * @returns {{value: number, staticNpv: number, flexibilityValue: number, proceed: boolean}} value = expanded NPV
 *   of the staged programme today (after paying a time-0 phase if it is worth proceeding)
 */
export function sequentialCompoundOption({ V, phases, r, sigma, q = 0, steps = 240 }) {
  if (!Array.isArray(phases) || phases.length === 0) throw new RangeError('phases must be a non-empty array');
  const ph = phases.slice().sort((a, b) => a.time - b.time);
  for (const p of ph) { requireNumber('phase time', p.time, { min: 0 }); requireNumber('phase cost', p.cost, { min: 0 }); }
  const T = ph[ph.length - 1].time;
  const staticNpv = V * Math.exp(-q * T) - ph.reduce((s, p) => s + p.cost * Math.exp(-r * p.time), 0);
  if (T === 0) {
    const tot = ph.reduce((s, p) => s + p.cost, 0);
    return { value: Math.max(V - tot, 0), staticNpv: V - tot, flexibilityValue: Math.max(V - tot, 0) - Math.max(V - tot, 0), proceed: V > tot };
  }
  const dt = T / steps;
  const gate = new Map();
  for (let k = 0; k < ph.length - 1; k++) { const s = Math.round(ph[k].time / dt); gate.set(s, (gate.get(s) ?? 0) + ph[k].cost); }
  const last = ph[ph.length - 1].cost;
  const res = valueLattice({ V0: V, sigma, r, T, q, steps, terminal: (v) => Math.max(v - last - (gate.get(steps) ?? 0), 0),
    exercise: (v, cont, t, step) => (gate.has(step) ? Math.max(cont - gate.get(step), 0) : cont) });
  return { value: res.value, staticNpv, flexibilityValue: res.value - Math.max(staticNpv, 0), proceed: res.value > 0 };
}

/**
 * Margrabe option to exchange asset 2 for asset 1 at T (payoff max(S1 − S2, 0)); the strike is itself uncertain,
 * e.g. switching between two outputs or fuels.
 * Ref: Margrabe, W. (1978), The value of an option to exchange one asset for another, J. Finance 33(1).
 * @param {{S1: number, S2: number, sigma1: number, sigma2: number, rho: number, T: number, q1?: number, q2?: number}} o
 *   values (currency), volatilities (1/√yr), correlation (−1…1), maturity (yr), yields (1/yr)
 * @returns {{value: number, sigma: number, d1: number, d2: number}} currency
 */
export function margrabe({ S1, S2, sigma1, sigma2, rho, T, q1 = 0, q2 = 0 }) {
  requireNumber('S1', S1, { gt: 0 }); requireNumber('S2', S2, { gt: 0 }); requireNumber('sigma1', sigma1, { min: 0 });
  requireNumber('sigma2', sigma2, { min: 0 }); requireNumber('rho', rho, { min: -1, max: 1 }); requireNumber('T', T, { gt: 0 });
  const sigma = Math.sqrt(sigma1 * sigma1 + sigma2 * sigma2 - 2 * rho * sigma1 * sigma2);
  if (!(sigma > 0)) return { value: Math.max(S1 * Math.exp(-q1 * T) - S2 * Math.exp(-q2 * T), 0), sigma: 0, d1: NaN, d2: NaN };
  const sq = sigma * Math.sqrt(T);
  const d1 = (Math.log(S1 / S2) + (q2 - q1 + 0.5 * sigma * sigma) * T) / sq; const d2 = d1 - sq;
  return { value: S1 * Math.exp(-q1 * T) * normCdf(d1) - S2 * Math.exp(-q2 * T) * normCdf(d2), sigma, d1, d2 };
}

/**
 * Project-value volatility by the consolidated (market-asset-disclaimer) approach: simulate the first-year
 * outcome of the cash-flow model, form the logarithmic return z = ln[(PV1 + CF1)/PV0] and take its standard
 * deviation as the annual volatility of project value.
 * Ref: Copeland, T. & Antikarov, V. (2001), Real Options: A Practitioner's Guide, ch. 9.
 * @param {{pv0: number, pv1: ArrayLike<number>, cf1?: ArrayLike<number>}} o pv0 = expected PV today (currency);
 *   pv1 = simulated PV at year 1 of the remaining cash flows; cf1 = simulated year-1 cash flow (currency)
 * @returns {{sigma: number, meanReturn: number, n: number}} volatility 1/√yr, mean log return
 */
export function projectVolatility({ pv0, pv1, cf1 }) {
  requireNumber('pv0', pv0, { gt: 0 });
  const n = pv1.length;
  if (n < 2) throw new RangeError('pv1 needs at least 2 samples');
  const z = new Float64Array(n); let m = 0;
  for (let k = 0; k < n; k++) {
    const v = pv1[k] + (cf1 ? cf1[k] : 0);
    if (!(v > 0)) throw new RangeError('PV1 + CF1 must be positive in every sample for a log return');
    z[k] = Math.log(v / pv0); m += z[k];
  }
  m /= n;
  let s = 0; for (let k = 0; k < n; k++) s += (z[k] - m) ** 2;
  return { sigma: Math.sqrt(s / (n - 1)), meanReturn: m, n };
}

/**
 * Expanded (strategic) NPV = static NPV + value of managerial flexibility.
 * Ref: Trigeorgis (1996), Real Options, ch. 4.
 * @param {{staticNpv: number, optionValue: number}} o currency
 * @returns {{expandedNpv: number, staticNpv: number, optionValue: number, decisionChanged: boolean}}
 */
export function expandedNpv({ staticNpv, optionValue }) {
  requireNumber('staticNpv', staticNpv); requireNumber('optionValue', optionValue, { min: 0 });
  const e = staticNpv + optionValue;
  return { expandedNpv: e, staticNpv, optionValue, decisionChanged: staticNpv < 0 && e >= 0 };
}

/**
 * Decision-tree rollback with chance nodes (shared with decision.js; re-exported for real-options work where a
 * tree with explicit probabilities is preferred over a lattice).
 * Ref: Magee, J.F. (1964), Decision trees for decision making, Harvard Business Review 42(4).
 * @param {object} tree see decision.rollbackTree
 * @param {object} [o] see decision.rollbackTree
 * @returns {object} rollback result
 */
export function decisionTreeRollback(tree, o) { return rollbackTree(tree, o); }

/**
 * Longstaff–Schwartz least-squares Monte Carlo for an American call or put under geometric Brownian motion,
 * with antithetic paths and a quadratic polynomial basis on in-the-money paths. Seeded.
 * Ref: Longstaff, F.A. & Schwartz, E.S. (2001), Valuing American options by simulation: a simple least-squares
 * approach, Review of Financial Studies 14(1).
 * @param {{S: number, K: number, r: number, sigma: number, T: number, q?: number, steps?: number, paths?: number,
 *   seed?: number, type?: 'call' | 'put'}} o as for blackScholes; paths = number of simulated paths (even)
 * @returns {{value: number, standardError: number, paths: number, steps: number}} currency
 */
export function longstaffSchwartz({ S, K, r, sigma, T, q = 0, steps = 50, paths = 20000, seed = 1, type = 'put' }) {
  common(S, K, r, sigma, T, q);
  requireNumber('steps', steps, { min: 1, integer: true }); requireNumber('paths', paths, { min: 100, integer: true });
  const sign = type === 'call' ? 1 : -1;
  const n = paths + (paths % 2); const dt = T / steps;
  const drift = (r - q - 0.5 * sigma * sigma) * dt; const vol = sigma * Math.sqrt(dt); const disc = Math.exp(-r * dt);
  const rng = makeRng(seed);
  const path = Array.from({ length: steps + 1 }, () => new Float64Array(n));
  path[0].fill(S);
  for (let i = 1; i <= steps; i++) {
    for (let k = 0; k < n; k += 2) {
      const z = rng.normal();
      path[i][k] = path[i - 1][k] * Math.exp(drift + vol * z);
      path[i][k + 1] = path[i - 1][k + 1] * Math.exp(drift - vol * z);
    }
  }
  const cash = new Float64Array(n);
  for (let k = 0; k < n; k++) cash[k] = Math.max(sign * (path[steps][k] - K), 0);
  for (let i = steps - 1; i >= 1; i--) {
    for (let k = 0; k < n; k++) cash[k] *= disc;
    // regress continuation value on 1, x, x² (x = S/K) for in-the-money paths
    const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; const b = [0, 0, 0]; let cnt = 0;
    for (let k = 0; k < n; k++) {
      const ex = sign * (path[i][k] - K);
      if (ex <= 0) continue;
      const x = path[i][k] / K; const basis = [1, x, x * x];
      for (let a = 0; a < 3; a++) { b[a] += basis[a] * cash[k]; for (let c = 0; c < 3; c++) A[a][c] += basis[a] * basis[c]; }
      cnt++;
    }
    if (cnt < 10) continue;
    let beta;
    try { beta = solveDense(A, b); } catch { continue; }
    for (let k = 0; k < n; k++) {
      const ex = sign * (path[i][k] - K);
      if (ex <= 0) continue;
      const x = path[i][k] / K;
      if (ex > beta[0] + beta[1] * x + beta[2] * x * x) cash[k] = ex;
    }
  }
  let m = 0; let s2 = 0;
  for (let k = 0; k < n; k += 2) { const v = 0.5 * (cash[k] + cash[k + 1]) * disc; m += v; s2 += v * v; }
  const np = n / 2; m /= np;
  const value = Math.max(m, Math.max(sign * (S - K), 0));
  return { value, standardError: Math.sqrt(Math.max(0, s2 / np - m * m) / (np - 1)), paths: n, steps };
}
