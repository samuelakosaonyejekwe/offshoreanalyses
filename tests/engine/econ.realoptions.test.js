import test from 'node:test';
import assert from 'node:assert/strict';
import {
  blackScholes, binomialCRR, valueLattice, realOptionLattice, optionToDefer, optionToExpand, optionToContract,
  optionToAbandon, optionToSwitch, sequentialCompoundOption, margrabe, projectVolatility, expandedNpv,
  decisionTreeRollback, longstaffSchwartz,
} from '../../src/engine/econ/realoptions.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} (±${tol})`);
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('Black–Scholes: S = K = 100, r = 5 %, σ = 20 %, T = 1 → call 10.4506, put 5.5735; Greeks and parity with yield', () => {
  // d1 = (0 + (0.05 + 0.02))/0.2 = 0.35, d2 = 0.15; N(0.35) = 0.63683, N(0.15) = 0.55962
  // call = 100 × 0.63683 − 100 e^(−0.05) × 0.55962 = 63.683 − 53.232 = 10.451
  const b = blackScholes({ S: 100, K: 100, r: 0.05, sigma: 0.2, T: 1 });
  close(b.d1, 0.35, 1e-12); close(b.d2, 0.15, 1e-12);
  close(b.call, 10.4506, 1e-4); close(b.put, 5.5735, 1e-4);
  close(b.greeks.call.delta, 0.63683, 1e-5); close(b.greeks.put.delta, 0.63683 - 1, 1e-5);
  // gamma = φ(d1)/(S σ √T) = 0.375240/20 = 0.018762; vega = S φ(d1) √T = 37.524
  close(b.greeks.call.gamma, 0.018762, 1e-6); close(b.greeks.call.vega, 37.524, 1e-3);
  // put–call parity with a yield: C − P = S e^(−qT) − K e^(−rT)
  const y = blackScholes({ S: 120, K: 100, r: 0.04, sigma: 0.3, T: 2.5, q: 0.03 });
  close(y.call - y.put, 120 * Math.exp(-0.075) - 100 * Math.exp(-0.1), 1e-10);
  // finite-difference check of delta, vega and rho with yield
  const f = (o) => blackScholes({ S: 120, K: 100, r: 0.04, sigma: 0.3, T: 2.5, q: 0.03, ...o }).call;
  close(y.greeks.call.delta, (f({ S: 120.01 }) - f({ S: 119.99 })) / 0.02, 1e-6);
  close(y.greeks.call.vega, (f({ sigma: 0.3001 }) - f({ sigma: 0.2999 })) / 0.0002, 1e-4);
  close(y.greeks.call.rho, (f({ r: 0.0401 }) - f({ r: 0.0399 })) / 0.0002, 1e-4);
  close(y.greeks.call.theta, -(f({ T: 2.5001 }) - f({ T: 2.4999 })) / 0.0002, 1e-4);
  assert.throws(() => blackScholes({ S: -1, K: 100, r: 0.05, sigma: 0.2, T: 1 }), RangeError);
});

test('CRR lattice: converges to Black–Scholes; American put of Hull (S = K = 50, r = 10 %, σ = 40 %, T = 5 months)', () => {
  const o = { S: 100, K: 95, r: 0.05, sigma: 0.25, T: 1.5, q: 0.02 };
  const bs = blackScholes(o);
  close(binomialCRR({ ...o, steps: 800, type: 'call' }).value, bs.call, 0.01);
  close(binomialCRR({ ...o, steps: 800, type: 'put' }).value, bs.put, 0.01);
  // Hull, Options, Futures and Other Derivatives, American put on a five-step tree: 4.49; fine trees approach 4.28
  const h = { S: 50, K: 50, r: 0.1, sigma: 0.4, T: 5 / 12, type: 'put', american: true };
  close(binomialCRR({ ...h, steps: 5 }).value, 4.49, 0.005);
  close(binomialCRR({ ...h, steps: 500 }).value, 4.283, 0.005);
  const a = binomialCRR({ ...h, steps: 500 });
  const e = binomialCRR({ ...h, american: false, steps: 500 });
  close(a.earlyExercisePremium, a.value - e.value, 1e-12);
  assert.ok(a.earlyExercisePremium > 0.1);
  close(e.value, blackScholes(h).put, 0.01);
  // independent five-step backward induction written out
  const dt = (5 / 12) / 5; const u = Math.exp(0.4 * Math.sqrt(dt)); const d = 1 / u;
  const p = (Math.exp(0.1 * dt) - d) / (u - d); const disc = Math.exp(-0.1 * dt);
  let v = Array.from({ length: 6 }, (_, j) => Math.max(50 - 50 * u ** (2 * j - 5), 0));
  for (let i = 4; i >= 0; i--) v = Array.from({ length: i + 1 }, (_, j) => Math.max(disc * (p * v[j + 1] + (1 - p) * v[j]), 50 - 50 * u ** (2 * j - i)));
  close(binomialCRR({ ...h, steps: 5 }).value, v[0], 1e-12);
  // an American call without yield is never exercised early
  close(binomialCRR({ S: 50, K: 50, r: 0.1, sigma: 0.4, T: 5 / 12, steps: 200, american: true }).earlyExercisePremium, 0, 1e-12);
  assert.throws(() => binomialCRR({ ...h, type: 'straddle' }), RangeError);
});

test('owner options on a leaking project value: the value without options is V for any yield (cash received while holding)', () => {
  // V = 100, q = 6 %/yr, T = 5: a lattice carrying the terminal value only would give 100 e^(−0.3) = 74.08
  for (const q of [0, 0.03, 0.06]) {
    const r = realOptionLattice({ V0: 100, sigma: 0.3, r: 0.05, T: 5, q, steps: 100, options: [{ type: 'abandon', salvage: 1e-9 }] });
    close(r.valueWithoutOptions, 100, 1e-9, `q = ${q}`);
    close(r.optionValue, 0, 1e-6);
  }
  // valueLattice itself: payout false is the claim on V_T (prepaid forward), payout true the asset
  const hold = { V0: 100, sigma: 0.3, r: 0.05, T: 5, q: 0.06, steps: 50, terminal: (V) => V };
  close(valueLattice(hold).value, 100 * Math.exp(-0.3), 1e-9);
  close(valueLattice({ ...hold, payout: true }).value, 100, 1e-9);
});

test('abandonment option: European limit = Black–Scholes–Merton put with yield; q = 0 American = CRR American put', () => {
  const o = { V: 100, salvage: 90, r: 0.05, sigma: 0.3, T: 3, q: 0.04 };
  // European: V(1 − e^(−qT)) + e^(−rT) E[max(V_T, X)] = V + put(S = V, K = X, q)
  const eu = optionToAbandon({ ...o, steps: 1000, american: false });
  const bs = blackScholes({ S: 100, K: 90, r: 0.05, sigma: 0.3, T: 3, q: 0.04 });
  close(eu.optionValue, bs.put, 0.02);
  close(eu.valueWithOptions, 100 + bs.put, 0.02);
  // no yield: the owner's American abandonment option is the American put
  const am0 = optionToAbandon({ ...o, q: 0, steps: 400 });
  close(am0.optionValue, binomialCRR({ S: 100, K: 90, r: 0.05, sigma: 0.3, T: 3, steps: 400, type: 'put', american: true }).value, 1e-9);
  // with yield the American option lies between the European value and the value of the lattice that omits the
  // cash flow received while holding (which overstates early exercise)
  const am = optionToAbandon({ ...o, steps: 400 });
  const noCash = valueLattice({ V0: 100, sigma: 0.3, r: 0.05, T: 3, q: 0.04, steps: 400, terminal: (V) => Math.max(V, 90),
    exercise: (V, cont) => Math.max(cont, 90) }).value - 100 * Math.exp(-0.12);
  assert.ok(am.optionValue >= eu.optionValue - 0.02 && am.optionValue < noCash - 0.3, `American ${am.optionValue}, European ${eu.optionValue}, without holding cash ${noCash}`);
  // independent backward induction with the holding cash flow written out (50 steps)
  const n = 50; const dt = 3 / n; const u = Math.exp(0.3 * Math.sqrt(dt)); const d = 1 / u;
  const p = (Math.exp(0.01 * dt) - d) / (u - d); const disc = Math.exp(-0.05 * dt); const pay = 1 - Math.exp(-0.04 * dt);
  let v = Array.from({ length: n + 1 }, (_, j) => Math.max(100 * u ** (2 * j - n), 90));
  for (let i = n - 1; i >= 0; i--) {
    v = Array.from({ length: i + 1 }, (_, j) => { const V = 100 * u ** (2 * j - i); return Math.max(V * pay + disc * (p * v[j + 1] + (1 - p) * v[j]), 90); });
  }
  close(optionToAbandon({ ...o, steps: n }).valueWithOptions, v[0], 1e-9);
});

test('expand and contract options at maturity reduce to scaled Black–Scholes calls and puts', () => {
  // expand by 1.5 for 40 at T: max(V, 1.5V − 40) = V + 0.5 max(V − 80, 0) → 0.5 × call(K = 80)
  const base = { V: 100, r: 0.04, sigma: 0.35, T: 2, q: 0.05, steps: 1000, american: false };
  const ex = optionToExpand({ ...base, factor: 1.5, cost: 40 });
  close(ex.optionValue, 0.5 * blackScholes({ S: 100, K: 80, r: 0.04, sigma: 0.35, T: 2, q: 0.05 }).call, 0.02);
  // contract to 0.7 saving 25: max(V, 0.7V + 25) = V + 0.3 max(83.33 − V, 0) → 0.3 × put(K = 83.33)
  const co = optionToContract({ ...base, factor: 0.7, saving: 25 });
  close(co.optionValue, 0.3 * blackScholes({ S: 100, K: 25 / 0.3, r: 0.04, sigma: 0.35, T: 2, q: 0.05 }).put, 0.02);
  // American exercise is worth at least as much; several options on one lattice are worth less than their sum
  const exA = optionToExpand({ ...base, american: true, factor: 1.5, cost: 40, steps: 300 });
  assert.ok(exA.optionValue >= ex.optionValue - 0.02);
  const both = realOptionLattice({ V0: 100, sigma: 0.35, r: 0.04, T: 2, q: 0.05, steps: 300,
    options: [{ type: 'expand', factor: 1.5, cost: 40 }, { type: 'abandon', salvage: 70 }] });
  const ab = optionToAbandon({ V: 100, salvage: 70, r: 0.04, sigma: 0.35, T: 2, q: 0.05, steps: 300 });
  assert.ok(both.optionValue <= exA.optionValue + ab.optionValue + 1e-9 && both.optionValue >= Math.max(exA.optionValue, ab.optionValue) - 1e-9);
  // switching to an alternative use worth more than the project today is exercised at once
  const sw = optionToSwitch({ V: 100, value: 130, cost: 10, r: 0.04, sigma: 0.2, T: 2, q: 0.05, steps: 100 });
  assert.ok(sw.valueWithOptions >= 120 - 1e-9);
  assert.throws(() => realOptionLattice({ V0: 100, sigma: 0.3, r: 0.05, T: 1, options: [] }), RangeError);
  assert.throws(() => realOptionLattice({ V0: 100, sigma: 0.3, r: 0.05, T: 1, options: [{ type: 'grow' }] }), RangeError);
  assert.throws(() => optionToExpand({ ...base, factor: 0.9, cost: 1 }), RangeError);
});

test('option to defer: without cost of delay it equals the European call and waiting has value; with delay cost early investment can be optimal', () => {
  const d = optionToDefer({ V: 100, I: 100, r: 0.05, sigma: 0.2, T: 1, steps: 600 });
  close(d.optionValue, 10.4506, 0.01); close(d.staticNpv, 0, 1e-12); close(d.valueOfWaiting, d.optionValue, 1e-12);
  assert.equal(d.investNow, false);
  const deep = optionToDefer({ V: 300, I: 100, r: 0.02, sigma: 0.1, T: 2, q: 0.15, steps: 300 });
  close(deep.optionValue, 200, 1e-6); assert.equal(deep.investNow, true);
});

test('sequential compound option: a single gate at T is the Black–Scholes call; an earlier gate lowers the value', () => {
  const one = sequentialCompoundOption({ V: 100, phases: [{ time: 2, cost: 90 }], r: 0.05, sigma: 0.3, steps: 800 });
  close(one.value, blackScholes({ S: 100, K: 90, r: 0.05, sigma: 0.3, T: 2 }).call, 0.03);
  close(one.staticNpv, 100 - 90 * Math.exp(-0.1), 1e-10);
  const two = sequentialCompoundOption({ V: 100, phases: [{ time: 1, cost: 10 }, { time: 2, cost: 90 }], r: 0.05, sigma: 0.3, steps: 400 });
  assert.ok(two.value < one.value && two.value > one.value - 10 * Math.exp(-0.05) - 1e-9);
  // all phases today: plain NPV decision
  const now = sequentialCompoundOption({ V: 100, phases: [{ time: 0, cost: 60 }], r: 0.05, sigma: 0.3 });
  close(now.value, 40, 1e-12); assert.equal(now.proceed, true);
  assert.throws(() => sequentialCompoundOption({ V: 100, phases: [], r: 0.05, sigma: 0.3 }), RangeError);
});

test('Margrabe exchange option: equal assets give S[2N(σ√T/2) − 1]; zero volatility of the ratio gives the intrinsic value', () => {
  // σ1 = 0.3, σ2 = 0.4, ρ = 0: σ = 0.5; T = 1 → d1 = 0.25, value = 100 (N(0.25) − N(−0.25)) = 100 × 0.19741 = 19.741
  const m = margrabe({ S1: 100, S2: 100, sigma1: 0.3, sigma2: 0.4, rho: 0, T: 1 });
  close(m.sigma, 0.5, 1e-12); close(m.d1, 0.25, 1e-12); close(m.value, 19.741, 1e-3);
  // against a riskless second asset it is a Black–Scholes call with strike S2 e^(rT) and r = 0
  close(margrabe({ S1: 100, S2: 90, sigma1: 0.25, sigma2: 0, rho: 0, T: 2 }).value, blackScholes({ S: 100, K: 90, r: 1e-12, sigma: 0.25, T: 2 }).call, 1e-6);
  const z = margrabe({ S1: 110, S2: 100, sigma1: 0.3, sigma2: 0.3, rho: 1, T: 1 });
  close(z.value, 10, 1e-12); assert.equal(z.sigma, 0);
});

test('project volatility, expanded NPV and tree rollback', () => {
  // log returns ln(PV1/PV0) = ±0.2, ±0.1: mean 0, sample variance (2 × 0.04 + 2 × 0.01)/3 = 0.033333
  const pv1 = [Math.exp(0.2), Math.exp(-0.2), Math.exp(0.1), Math.exp(-0.1)].map((v) => 50 * v);
  const v = projectVolatility({ pv0: 50, pv1 });
  close(v.sigma, Math.sqrt(0.1 / 3), 1e-12); close(v.meanReturn, 0, 1e-12); assert.equal(v.n, 4);
  close(projectVolatility({ pv0: 50, pv1: pv1.map((x) => x - 5), cf1: [5, 5, 5, 5] }).sigma, v.sigma, 1e-12);
  assert.throws(() => projectVolatility({ pv0: 50, pv1: [10, -60] }), RangeError);
  const e = expandedNpv({ staticNpv: -5, optionValue: 12 });
  close(e.expandedNpv, 7, 1e-12); assert.equal(e.decisionChanged, true);
  assert.equal(expandedNpv({ staticNpv: 5, optionValue: 1 }).decisionChanged, false);
  assert.equal(typeof decisionTreeRollback, 'function');
});

test('Longstaff–Schwartz Monte Carlo agrees with the lattice for an American put within sampling error', () => {
  const o = { S: 36, K: 40, r: 0.06, sigma: 0.2, T: 1 };
  // Longstaff & Schwartz (2001) Table 1: finite-difference value 4.478 for this case
  const lat = binomialCRR({ ...o, steps: 1000, type: 'put', american: true }).value;
  close(lat, 4.478, 0.01);
  const mc = longstaffSchwartz({ ...o, steps: 50, paths: 20000, seed: 3 });
  assert.ok(Math.abs(mc.value - lat) < 4 * mc.standardError + 0.03, `LSM ${mc.value} ± ${mc.standardError} vs ${lat}`);
  assert.equal(longstaffSchwartz({ ...o, steps: 20, paths: 2000, seed: 3 }).value, longstaffSchwartz({ ...o, steps: 20, paths: 2000, seed: 3 }).value);
});
