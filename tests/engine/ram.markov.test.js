import test from 'node:test';
import assert from 'node:assert/strict';
import {
  buildGenerator, ctmcSteadyState, ctmcTransient, ctmcAbsorption, ctmcAvailability, dtmcStep, dtmcStationary,
  dtmcAbsorbing, semiMarkovSteadyState, model1oo2, twoStateUnit, degradationChain,
} from '../../src/engine/ram/markov.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + 1e-15, `${msg || ''} got ${a}, expected ${b}`);

test('two-state repairable unit: steady availability μ/(λ + μ) and transient A(t) = μ/(λ+μ) + λ/(λ+μ)·e^(−(λ+μ)t)', () => {
  const lambda = 1e-3;
  const mu = 5e-2;
  const { Q } = buildGenerator({ states: ['up', 'down'], transitions: [{ from: 'up', to: 'down', rate: lambda }, { from: 'down', to: 'up', rate: mu }] });
  assert.deepEqual(Q, [[-lambda, lambda], [mu, -mu]]);
  const pi = ctmcSteadyState(Q);
  close(pi[0], mu / (lambda + mu), 1e-12);
  close(pi[0], 0.980392, 1e-6);
  for (const t of [1, 10, 50, 200, 5000]) {
    const p = ctmcTransient(Q, [1, 0], t);
    close(p[0], mu / (lambda + mu) + (lambda / (lambda + mu)) * Math.exp(-(lambda + mu) * t), 1e-10, `t = ${t}`);
    close(p[0] + p[1], 1, 1e-12);
    close(twoStateUnit({ lambda, mu, t }).pointAvailability, p[0], 1e-10);
  }
  // Starting down: A(t) = μ/(λ+μ)·(1 − e^(−(λ+μ)t))
  close(ctmcTransient(Q, [0, 1], 20)[0], (mu / (lambda + mu)) * (1 - Math.exp(-(lambda + mu) * 20)), 1e-10);
  const av = ctmcAvailability(Q, [0], { p0: [1, 0], t: 10 });
  close(av.availability, mu / (lambda + mu), 1e-12);
  close(av.failureFrequency, (lambda * mu) / (lambda + mu), 1e-12);
  close(av.meanUpTime, 1 / lambda, 1e-10);
  close(av.meanDownTime, 1 / mu, 1e-10);
});

test('uniformisation on a stiff chain agrees with the analytic solution of a pure-death (Erlang) process', () => {
  // 0 → 1 → 2 with the same rate r: p0 = e^(−rt), p1 = rt·e^(−rt), p2 = 1 − e^(−rt)(1 + rt); r·t = 800 forces segmentation.
  const r = 4;
  const Q = [[-r, r, 0], [0, -r, r], [0, 0, 0]];
  for (const t of [0.1, 1, 3]) {
    const p = ctmcTransient(Q, [1, 0, 0], t);
    close(p[0], Math.exp(-r * t), 1e-9);
    close(p[1], r * t * Math.exp(-r * t), 1e-9);
    close(p[2], 1 - Math.exp(-r * t) * (1 + r * t), 1e-9);
  }
  const far = ctmcTransient(Q, [1, 0, 0], 200);
  close(far[2], 1, 1e-12);
});

test('mean time to absorption: redundant pair with repair, MTTF = (3λ + μ)/(2λ²)', () => {
  const lambda = 1e-3;
  const mu = 0.1;
  const m = model1oo2({ lambda, mu });
  close(m.mttf, (3 * lambda + mu) / (2 * lambda * lambda), 1e-10);
  // Steady state of the birth–death chain (one repair crew): π ∝ [1, 2λ/μ, 2λ²/μ²]
  const norm = 1 + (2 * lambda) / mu + (2 * lambda * lambda) / (mu * mu);
  close(m.unavailability, (2 * lambda * lambda) / (mu * mu) / norm, 1e-10);
  // Two crews: π2 ∝ λ²/μ² → unavailability (λ/(λ+μ))² (independent units)
  close(model1oo2({ lambda, mu, repairCrews: 2 }).unavailability, (lambda / (lambda + mu)) ** 2, 1e-10);
  // No repair: MTTF = 3/(2λ)
  close(model1oo2({ lambda, mu: 0 }).mttf, 1.5 / lambda, 1e-10);
  // Sojourn times and visits: from state 0, time in state 1 = MTTF − time in state 0.
  const ab = ctmcAbsorption(m.Q, [2], [1, 0, 0]);
  close(ab.mean, m.mttf, 1e-12);
  close(ab.sojourn[0][0] + ab.sojourn[0][1], m.mttf, 1e-12);
  close(ab.absorption[0][0], 1, 1e-12);
  // Expected visits to state 1 before failure: 1/P(fail from 1) = (λ + μ)/λ
  close(ab.visits[0][1], (lambda + mu) / lambda, 1e-10);
});

test('β-factor common cause: MTTF with β against the closed form of the three-state chain', () => {
  const lambda = 1e-3;
  const mu = 0.1;
  const beta = 0.1;
  const m = model1oo2({ lambda, mu, beta });
  // τ0 = [1 + a·(1/(λ+μ))]/(a + b − a·μ/(λ+μ)) with a = 2(1−β)λ, b = βλ
  const a = 2 * (1 - beta) * lambda;
  const b = beta * lambda;
  close(m.mttf, (1 + a / (lambda + mu)) / (a + b - (a * mu) / (lambda + mu)), 1e-10);
  assert.ok(m.mttf < model1oo2({ lambda, mu }).mttf / 5, 'common cause dominates a well-repaired pair');
  assert.ok(m.unavailability > model1oo2({ lambda, mu }).unavailability);
});

test('discrete-time chains: n-step distribution, stationary distribution, gambler\'s ruin absorption', () => {
  const P = [[0.9, 0.1], [0.5, 0.5]];
  const pi = dtmcStationary(P);
  close(pi[0], 5 / 6, 1e-12);
  // Two-state closed form: p_n(0) = π0 + (1 − π0)·(1 − a − b)^n with a = 0.1, b = 0.5
  close(dtmcStep(P, [1, 0], 5)[0], 5 / 6 + (1 / 6) * 0.4 ** 5, 1e-12);
  // Gambler's ruin on 0..4 with p = 0.5: from state i, P(reach 4) = i/4 and expected duration i·(4 − i).
  const G = [[1, 0, 0, 0, 0], [0.5, 0, 0.5, 0, 0], [0, 0.5, 0, 0.5, 0], [0, 0, 0.5, 0, 0.5], [0, 0, 0, 0, 1]];
  const r = dtmcAbsorbing(G, [0, 4]);
  assert.deepEqual(r.transient, [1, 2, 3]);
  [1, 2, 3].forEach((i, k) => {
    close(r.absorption[k][1], i / 4, 1e-12);
    close(r.expectedSteps[k], i * (4 - i), 1e-12);
  });
});

test('semi-Markov steady state: alternating up/down with arbitrary sojourn distributions → MUT/(MUT + MDT)', () => {
  const r = semiMarkovSteadyState({ P: [[0, 1], [1, 0]], meanSojourn: [900, 100], upStates: [0] });
  close(r.availability, 0.9, 1e-12);
  // Three states: operation → (0.7 minor repair | 0.3 major repair) → operation; means 1000, 10, 200 h.
  const s = semiMarkovSteadyState({ P: [[0, 0.7, 0.3], [1, 0, 0], [1, 0, 0]], meanSojourn: [1000, 10, 200], upStates: [0] });
  close(s.availability, 1000 / (1000 + 7 + 60), 1e-12);
  close(s.timeFraction.reduce((a, b) => a + b, 0), 1, 1e-12);
});

test('degradation chain: mean time to failure is the sum of stage means; maintenance extends it', () => {
  const d = degradationChain({ rates: [0.1, 0.2, 0.5] });
  close(d.meanTimeToFailure, 10 + 5 + 2, 1e-10);
  const withRepair = degradationChain({ rates: [0.1, 0.2, 0.5], repairRate: 1 });
  // Renewal argument: availability = MTTF/(MTTF + MTTR) = 17/18
  close(withRepair.availability, 17 / 18, 1e-10);
  const maintained = degradationChain({ rates: [0.1, 0.2, 0.5], maintenanceRates: [0, 0, 2] });
  assert.ok(maintained.meanTimeToFailure > 17);
  // From stage 2 the chance of failing before maintenance is 0.5/2.5 = 0.2 → 5 passes on average:
  // MTTF = 5·(10 + 5 + 1/2.5) = 77
  close(maintained.meanTimeToFailure, 77, 1e-10);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => buildGenerator({ states: ['a', 'a'], transitions: [] }), /duplicate/);
  assert.throws(() => buildGenerator({ states: ['a', 'b'], transitions: [{ from: 'a', to: 'c', rate: 1 }] }), /unknown state/);
  assert.throws(() => ctmcTransient([[-1, 1], [1, -1]], [1, 0], -1), /t must/);
  assert.throws(() => dtmcStep([[0.5, 0.4], [0, 1]], [1, 0], 1), /row 0/);
  assert.throws(() => model1oo2({ lambda: 1, mu: 1, beta: 2 }), /beta/);
});
