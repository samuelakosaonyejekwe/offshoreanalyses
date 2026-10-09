// Markov models for reliability and availability: continuous-time chains (steady state, transient, absorption),
// discrete-time chains, a simple semi-Markov model and standard redundancy models.
// Rates and times must use consistent units (e.g. 1/h and h).

import { solveDense } from '../math/linalg.js';

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function checkSquare(M, name) {
  need(Array.isArray(M) && M.length > 0, `${name} must be a non-empty square matrix`);
  for (const row of M) need(row.length === M.length, `${name} must be square`);
}

/**
 * Build the generator (transition-rate) matrix Q of a continuous-time Markov chain: q_ij = rate i → j for
 * i ≠ j, q_ii = −Σ_j q_ij.
 * @param {object} o
 * @param {string[]} o.states state names
 * @param {{from:string, to:string, rate:number}[]} o.transitions rates [1/time]; duplicates add
 * @returns {{Q:number[][], states:string[], index:Object<string, number>}}
 */
export function buildGenerator({ states, transitions }) {
  need(Array.isArray(states) && states.length > 0, 'states must be a non-empty array');
  const index = {};
  states.forEach((s, i) => {
    need(index[s] === undefined, `states: duplicate state "${s}"`);
    index[s] = i;
  });
  const n = states.length;
  const Q = Array.from({ length: n }, () => new Array(n).fill(0));
  for (const tr of transitions || []) {
    const i = index[tr.from];
    const j = index[tr.to];
    need(i !== undefined && j !== undefined, `transitions: unknown state in "${tr.from}" -> "${tr.to}"`);
    need(i !== j, 'transitions: self-transitions are not allowed');
    need(tr.rate >= 0, 'transitions: rate must be non-negative');
    Q[i][j] += tr.rate;
    Q[i][i] -= tr.rate;
  }
  return { Q, states: states.slice(), index };
}

/**
 * Steady-state probabilities of an irreducible continuous-time Markov chain: π·Q = 0, Σπ = 1.
 * Reference: Rausand & Høyland, System Reliability Theory (2nd ed.); IEC 61165.
 * @param {number[][]} Q generator matrix [1/time]
 * @returns {number[]} stationary probabilities [-]
 */
export function ctmcSteadyState(Q) {
  checkSquare(Q, 'Q');
  const n = Q.length;
  const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => Q[j][i]));
  const b = new Array(n).fill(0);
  A[n - 1] = new Array(n).fill(1);
  b[n - 1] = 1;
  const pi = Array.from(solveDense(A, b));
  for (let i = 0; i < n; i++) if (pi[i] < 0 && pi[i] > -1e-12) pi[i] = 0;
  return pi;
}

/**
 * Transient state probabilities p(t) = p(0)·exp(Q·t) by uniformisation (Jensen's method):
 * p(t) = Σ_k e^(−Λt)(Λt)^k/k! · p(0)·P^k with P = I + Q/Λ, Λ >= max|q_ii|. Long horizons are split into
 * segments so that the Poisson weights never underflow.
 * Reference: Jensen (1953); Grassmann (1977); Stewart, Introduction to the Numerical Solution of Markov Chains
 * (1994).
 * @param {number[][]} Q generator matrix [1/time]
 * @param {number[]} p0 initial distribution [-]
 * @param {number} t time [time]
 * @param {object} [o]
 * @param {number} [o.tol=1e-12] truncation error of the Poisson series [-]
 * @returns {number[]} state probabilities at t [-]
 */
export function ctmcTransient(Q, p0, t, { tol = 1e-12 } = {}) {
  checkSquare(Q, 'Q');
  const n = Q.length;
  need(Array.isArray(p0) || ArrayBuffer.isView(p0), 'p0 must be an array');
  need(p0.length === n, 'p0 must have one entry per state');
  need(t >= 0, 't must be non-negative');
  let lam = 0;
  for (let i = 0; i < n; i++) lam = Math.max(lam, -Q[i][i]);
  let p = Float64Array.from(p0);
  if (lam === 0 || t === 0) return Array.from(p);
  lam *= 1.02;
  const P = Array.from({ length: n }, (_, i) => Float64Array.from({ length: n }, (_, j) => (i === j ? 1 : 0) + Q[i][j] / lam));
  const segments = Math.max(1, Math.ceil((lam * t) / 100));
  const x = (lam * t) / segments;
  const kMax = Math.ceil(x + 12 * Math.sqrt(x) + 30);
  let v = new Float64Array(n);
  let w = new Float64Array(n);
  for (let s = 0; s < segments; s++) {
    const acc = new Float64Array(n);
    v.set(p);
    let weight = Math.exp(-x);
    let cum = 0;
    for (let k = 0; k <= kMax; k++) {
      for (let j = 0; j < n; j++) acc[j] += weight * v[j];
      cum += weight;
      if (1 - cum < tol && k > x) break;
      w.fill(0);
      for (let i = 0; i < n; i++) {
        const vi = v[i];
        if (vi === 0) continue;
        const row = P[i];
        for (let j = 0; j < n; j++) w[j] += vi * row[j];
      }
      const tmp = v; v = w; w = tmp;
      weight *= x / (k + 1);
    }
    // Renormalise by the captured Poisson mass so truncation does not leak probability across segments.
    for (let j = 0; j < n; j++) acc[j] /= cum;
    p = acc;
  }
  return Array.from(p);
}

/**
 * Absorption analysis of a continuous-time chain: mean time to absorption from each transient state, expected
 * time spent in and expected number of visits to each transient state, and absorption probabilities.
 *   N = (−Q_TT)⁻¹ (expected sojourn times),  τ = N·1,  visits_ij = N_ij·|q_jj|,  B = N·Q_TA.
 * Reference: Rausand & Høyland (mean time to system failure); Kemeny & Snell, Finite Markov Chains.
 * @param {number[][]} Q generator [1/time]
 * @param {number[]} absorbing indices of the absorbing (failed) states
 * @param {number[]} [p0] initial distribution; when given, `mean` is the mean time to absorption from it
 * @returns {{transient:number[], meanTime:number[], sojourn:number[][], visits:number[][], absorption:number[][],
 *   mean:(number|null)}} meanTime[i] from transient state transient[i] [time]; sojourn[i][j] time in transient j
 *   starting from transient i [time]; absorption[i][a] probability of ending in absorbing[a]
 */
export function ctmcAbsorption(Q, absorbing, p0) {
  checkSquare(Q, 'Q');
  const n = Q.length;
  need(Array.isArray(absorbing) && absorbing.length > 0, 'absorbing must be a non-empty array of state indices');
  const isAbs = new Array(n).fill(false);
  for (const a of absorbing) {
    need(Number.isInteger(a) && a >= 0 && a < n, 'absorbing: state index out of range');
    isAbs[a] = true;
  }
  const tr = [];
  for (let i = 0; i < n; i++) if (!isAbs[i]) tr.push(i);
  need(tr.length > 0, 'at least one transient state is required');
  const m = tr.length;
  const A = tr.map((i) => tr.map((j) => -Q[i][j]));
  const I = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : 0)));
  // Columns of the inverse: solve A·X = I (multi right-hand side).
  const N = solveDense(A, I).map((r) => Array.from(r));
  const meanTime = N.map((row) => row.reduce((s, v) => s + v, 0));
  const visits = N.map((row) => row.map((v, j) => v * -Q[tr[j]][tr[j]]));
  const absorption = N.map((row) => absorbing.map((a) => row.reduce((s, v, j) => s + v * Q[tr[j]][a], 0)));
  let mean = null;
  if (p0) {
    need(p0.length === n, 'p0 must have one entry per state');
    mean = tr.reduce((s, i, k) => s + p0[i] * meanTime[k], 0);
  }
  return { transient: tr, meanTime, sojourn: N, visits, absorption, mean };
}

/**
 * Availability measures of a repairable system modelled as a continuous-time chain.
 *   A = Σ_{i up} π_i,  failure frequency ω = Σ_{i up} Σ_{j down} π_i·q_ij,  MUT = A/ω,  MDT = (1 − A)/ω.
 * Reference: Rausand & Høyland; IEC 61165.
 * @param {number[][]} Q generator [1/time]
 * @param {number[]} upStates indices of the functioning states
 * @param {object} [o]
 * @param {number[]} [o.p0] initial distribution — with `t`, the point availability A(t) is also returned
 * @param {number} [o.t] time [time]
 * @returns {{availability:number, unavailability:number, failureFrequency:number, meanUpTime:number,
 *   meanDownTime:number, steadyState:number[], pointAvailability:(number|null)}}
 */
export function ctmcAvailability(Q, upStates, { p0, t } = {}) {
  checkSquare(Q, 'Q');
  const n = Q.length;
  need(Array.isArray(upStates) && upStates.length > 0, 'upStates must be a non-empty array');
  const up = new Array(n).fill(false);
  for (const i of upStates) {
    need(Number.isInteger(i) && i >= 0 && i < n, 'upStates: state index out of range');
    up[i] = true;
  }
  const pi = ctmcSteadyState(Q);
  let A = 0;
  let w = 0;
  for (let i = 0; i < n; i++) {
    if (!up[i]) continue;
    A += pi[i];
    for (let j = 0; j < n; j++) if (!up[j]) w += pi[i] * Q[i][j];
  }
  let pointAvailability = null;
  if (p0 && t !== undefined) {
    const p = ctmcTransient(Q, p0, t);
    pointAvailability = p.reduce((s, v, i) => s + (up[i] ? v : 0), 0);
  }
  return {
    availability: A, unavailability: 1 - A, failureFrequency: w, meanUpTime: w > 0 ? A / w : Infinity,
    meanDownTime: w > 0 ? (1 - A) / w : 0, steadyState: pi, pointAvailability,
  };
}

function checkStochastic(P) {
  checkSquare(P, 'P');
  P.forEach((row, i) => {
    let s = 0;
    for (const v of row) { need(v >= -1e-15, 'P entries must be non-negative'); s += v; }
    need(Math.abs(s - 1) < 1e-9, `P row ${i} must sum to 1`);
  });
}

/**
 * State distribution of a discrete-time Markov chain after n steps: p_n = p_0·Pⁿ.
 * @param {number[][]} P one-step transition matrix (rows sum to 1)
 * @param {number[]} p0 initial distribution
 * @param {number} steps number of steps, integer >= 0
 * @returns {number[]} distribution after `steps` steps
 */
export function dtmcStep(P, p0, steps) {
  checkStochastic(P);
  const n = P.length;
  need(p0.length === n, 'p0 must have one entry per state');
  need(Number.isInteger(steps) && steps >= 0, 'steps must be a non-negative integer');
  let p = Array.from(p0);
  for (let s = 0; s < steps; s++) {
    const q = new Array(n).fill(0);
    for (let i = 0; i < n; i++) if (p[i] !== 0) for (let j = 0; j < n; j++) q[j] += p[i] * P[i][j];
    p = q;
  }
  return p;
}

/**
 * Stationary distribution of an irreducible discrete-time chain: π = π·P, Σπ = 1.
 * @param {number[][]} P transition matrix
 * @returns {number[]} stationary distribution
 */
export function dtmcStationary(P) {
  checkStochastic(P);
  const n = P.length;
  const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => P[j][i] - (i === j ? 1 : 0)));
  const b = new Array(n).fill(0);
  A[n - 1] = new Array(n).fill(1);
  b[n - 1] = 1;
  return Array.from(solveDense(A, b));
}

/**
 * Absorbing discrete-time chain: fundamental matrix N = (I − Q)⁻¹ (expected visits), expected number of steps
 * to absorption and absorption probabilities B = N·R.
 * Reference: Kemeny & Snell, Finite Markov Chains (1960).
 * @param {number[][]} P transition matrix
 * @param {number[]} absorbing indices of absorbing states
 * @returns {{transient:number[], fundamental:number[][], expectedSteps:number[], absorption:number[][]}}
 */
export function dtmcAbsorbing(P, absorbing) {
  checkStochastic(P);
  const n = P.length;
  need(Array.isArray(absorbing) && absorbing.length > 0, 'absorbing must be a non-empty array');
  const isAbs = new Array(n).fill(false);
  for (const a of absorbing) {
    need(Number.isInteger(a) && a >= 0 && a < n, 'absorbing: state index out of range');
    isAbs[a] = true;
  }
  const tr = [];
  for (let i = 0; i < n; i++) if (!isAbs[i]) tr.push(i);
  need(tr.length > 0, 'at least one transient state is required');
  const m = tr.length;
  const A = tr.map((i, r) => tr.map((j, c) => (r === c ? 1 : 0) - P[i][j]));
  const I = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : 0)));
  const N = solveDense(A, I).map((r) => Array.from(r));
  return {
    transient: tr, fundamental: N, expectedSteps: N.map((row) => row.reduce((s, v) => s + v, 0)),
    absorption: N.map((row) => absorbing.map((a) => row.reduce((s, v, j) => s + v * P[tr[j]][a], 0))),
  };
}

/**
 * Semi-Markov process, steady state only: states change according to an embedded discrete chain P and the
 * sojourn time in state i has an arbitrary distribution with mean h_i. Long-run fraction of time in state i:
 *   φ_i = π_i·h_i / Σ_j π_j·h_j   with π the stationary distribution of P.
 * Only the mean sojourn times matter for the steady state; transient semi-Markov behaviour is not modelled.
 * Reference: Ross, Stochastic Processes (2nd ed.); Limnios & Oprişan, Semi-Markov Processes and Reliability.
 * @param {object} o
 * @param {number[][]} o.P embedded transition matrix (zero diagonal for a pure jump chain)
 * @param {number[]} o.meanSojourn mean sojourn time in each state [time]
 * @param {number[]} [o.upStates] indices of functioning states; when given, the availability is returned
 * @returns {{timeFraction:number[], embedded:number[], meanCycleTime:number, availability:(number|null)}}
 */
export function semiMarkovSteadyState({ P, meanSojourn, upStates }) {
  const pi = dtmcStationary(P);
  need(Array.isArray(meanSojourn) && meanSojourn.length === P.length, 'meanSojourn must have one entry per state');
  let den = 0;
  for (let i = 0; i < pi.length; i++) {
    need(meanSojourn[i] >= 0, 'meanSojourn must be non-negative');
    den += pi[i] * meanSojourn[i];
  }
  need(den > 0, 'meanSojourn must not be all zero');
  const frac = pi.map((v, i) => (v * meanSojourn[i]) / den);
  return { timeFraction: frac, embedded: pi, meanCycleTime: den, availability: upStates ? upStates.reduce((s, i) => s + frac[i], 0) : null };
}

/**
 * Redundant pair (1-out-of-2, active) of identical repairable units with β-factor common-cause failures.
 * States: 0 = both working, 1 = one failed, 2 = both failed (system down).
 *   0 → 1: 2(1 − β)λ,  0 → 2: βλ,  1 → 2: λ,  1 → 0: μ,  2 → 1: μ·crews (crews = 1 or 2).
 * MTTF is the mean time to first system failure from state 0 (state 2 absorbing); for β = 0 and one crew it
 * equals (3λ + μ)/(2λ²).
 * Reference: Rausand & Høyland (parallel structure of two components); β-factor model of Fleming (1975);
 * IEC 61508-6 Annex D (β-factor).
 * @param {object} o
 * @param {number} o.lambda failure rate of one unit [1/time]
 * @param {number} o.mu repair rate of one unit [1/time]
 * @param {number} [o.beta=0] common-cause fraction [-]
 * @param {number} [o.repairCrews=1] 1 or 2
 * @returns {{Q:number[][], states:string[], availability:number, unavailability:number, mttf:number,
 *   failureFrequency:number, steadyState:number[]}}
 */
export function model1oo2({ lambda, mu, beta = 0, repairCrews = 1 }) {
  need(lambda > 0, 'lambda must be positive');
  need(mu >= 0, 'mu must be non-negative');
  need(beta >= 0 && beta <= 1, 'beta must be in [0, 1]');
  need(repairCrews === 1 || repairCrews === 2, 'repairCrews must be 1 or 2');
  const { Q, states } = buildGenerator({
    states: ['2 up', '1 up', '0 up'],
    transitions: [
      { from: '2 up', to: '1 up', rate: 2 * (1 - beta) * lambda },
      { from: '2 up', to: '0 up', rate: beta * lambda },
      { from: '1 up', to: '0 up', rate: lambda },
      { from: '1 up', to: '2 up', rate: mu },
      { from: '0 up', to: '1 up', rate: mu * repairCrews },
    ],
  });
  const mttf = ctmcAbsorption(Q, [2]).meanTime[0];
  if (mu === 0) return { Q, states, availability: 0, unavailability: 1, mttf, failureFrequency: 0, steadyState: [0, 0, 1] };
  const av = ctmcAvailability(Q, [0, 1]);
  return { Q, states, availability: av.availability, unavailability: av.unavailability, mttf, failureFrequency: av.failureFrequency, steadyState: av.steadyState };
}

/**
 * Two-state repairable unit: point availability A(t) = μ/(λ + μ) + λ/(λ + μ)·exp(−(λ + μ)t) (starting up),
 * steady state μ/(λ + μ).
 * Reference: Rausand & Høyland
 * @param {object} o
 * @param {number} o.lambda failure rate [1/time]
 * @param {number} o.mu repair rate [1/time]
 * @param {number} [o.t] time [time]
 * @returns {{availability:number, pointAvailability:(number|null), Q:number[][]}}
 */
export function twoStateUnit({ lambda, mu, t }) {
  need(lambda >= 0 && mu >= 0 && lambda + mu > 0, 'lambda and mu must be non-negative and not both zero');
  const s = lambda + mu;
  return {
    availability: mu / s,
    pointAvailability: t === undefined ? null : mu / s + (lambda / s) * Math.exp(-s * t),
    Q: [[-lambda, lambda], [mu, -mu]],
  };
}

/**
 * Degradation chain: states 0 (new) → 1 → … → n (failed) with stage rates, optional repair from the failed state
 * to new and optional maintenance transitions from intermediate states back to new.
 * @param {object} o
 * @param {number[]} o.rates degradation rate out of each state 0..n−1 [1/time]
 * @param {number} [o.repairRate=0] rate failed → new [1/time]
 * @param {number[]} [o.maintenanceRates] rate from each state 0..n−1 back to new [1/time] (entry 0 ignored)
 * @returns {{Q:number[][], states:string[], meanTimeToFailure:number, steadyState:(number[]|null),
 *   availability:(number|null)}} meanTimeToFailure from new with the failed state absorbing [time]
 */
export function degradationChain({ rates, repairRate = 0, maintenanceRates }) {
  need(Array.isArray(rates) && rates.length > 0, 'rates must be a non-empty array');
  const n = rates.length;
  const states = Array.from({ length: n + 1 }, (_, i) => (i === n ? 'failed' : `stage ${i}`));
  const transitions = [];
  rates.forEach((r, i) => {
    need(r > 0, 'rates must be positive');
    transitions.push({ from: states[i], to: states[i + 1], rate: r });
  });
  if (maintenanceRates) {
    need(maintenanceRates.length === n, 'maintenanceRates must have one entry per degradation state');
    for (let i = 1; i < n; i++) if (maintenanceRates[i] > 0) transitions.push({ from: states[i], to: states[0], rate: maintenanceRates[i] });
  }
  need(repairRate >= 0, 'repairRate must be non-negative');
  const absorbingQ = buildGenerator({ states, transitions }).Q;
  const mttfVal = ctmcAbsorption(absorbingQ, [n]).meanTime[0];
  if (repairRate === 0) return { Q: absorbingQ, states, meanTimeToFailure: mttfVal, steadyState: null, availability: null };
  transitions.push({ from: states[n], to: states[0], rate: repairRate });
  const { Q } = buildGenerator({ states, transitions });
  const pi = ctmcSteadyState(Q);
  return { Q, states, meanTimeToFailure: mttfVal, steadyState: pi, availability: 1 - pi[n] };
}
