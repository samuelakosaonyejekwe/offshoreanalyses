import test from 'node:test';
import assert from 'node:assert/strict';
import {
  basicEventProbability, minimalCutSets, inclusionExclusion, evaluateFaultTree, evaluateEventTree, fmeca,
} from '../../src/engine/ram/faulttree.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + 1e-300, `${msg || ''} got ${a}, expected ${b}`);
const B = (id, p) => ({ type: 'basic', id, p });

test('basic-event models: fixed, mission time, repairable (steady and transient), periodically tested', () => {
  assert.equal(basicEventProbability({ p: 0.3 }), 0.3);
  close(basicEventProbability({ lambda: 1e-4, t: 1000 }), 1 - Math.exp(-0.1), 1e-14);
  close(basicEventProbability({ lambda: 1e-3, mu: 0.1 }), 1e-3 / 0.101, 1e-14);
  close(basicEventProbability({ lambda: 1e-3, mu: 0.1, t: 10 }), (1e-3 / 0.101) * (1 - Math.exp(-1.01)), 1e-14);
  // Mean unavailability over a test interval: 1 − (1 − e^−λτ)/(λτ) ≈ λτ/2 for small λτ
  close(basicEventProbability({ lambda: 1e-5, testInterval: 8760 }), 1 - (1 - Math.exp(-0.0876)) / 0.0876, 1e-12);
  close(basicEventProbability({ lambda: 1e-9, testInterval: 1 }), 0.5e-9, 1e-6);
});

test('repeated event: TOP = (A or B) and (A or C) = A or (B and C), hand calculation', () => {
  const pA = 0.01; const pB = 0.1; const pC = 0.2;
  const tree = { type: 'and', children: [{ type: 'or', children: [B('A', pA), B('B', pB)] }, { type: 'or', children: [B('A', pA), B('C', pC)] }] };
  const r = evaluateFaultTree(tree);
  // Exact: pA + pB·pC − pA·pB·pC = 0.01 + 0.02 − 0.0002 = 0.0298
  close(r.probability, 0.0298, 1e-12);
  // Treating the two OR gates as independent would give 0.109·0.208 = 0.02267 — wrong.
  assert.ok(Math.abs(r.probability - 0.109 * 0.208) > 5e-3);
  assert.deepEqual(minimalCutSets(tree), [['A'], ['B', 'C']]);
  close(r.rareEvent, 0.03, 1e-12);
  close(r.minCutUpperBound, 1 - 0.99 * 0.98, 1e-12);
  assert.ok(r.probability <= r.minCutUpperBound && r.minCutUpperBound <= r.rareEvent);
  close(inclusionExclusion(minimalCutSets(tree), { A: pA, B: pB, C: pC }), 0.0298, 1e-12);
  // Importance, by hand: I_B(A) = 1 − pB·pC = 0.98; I_B(B) = pC·(1 − pA) = 0.198; I_B(C) = pB·(1 − pA) = 0.099
  const imp = Object.fromEntries(r.importance.map((x) => [x.id, x]));
  close(imp.A.birnbaum, 0.98, 1e-12);
  close(imp.B.birnbaum, 0.198, 1e-12);
  close(imp.C.birnbaum, 0.099, 1e-12);
  close(imp.A.criticality, (0.98 * 0.01) / 0.0298, 1e-12);
  // Fussell–Vesely: P(cut sets containing the event)/P(top): A → 0.01/0.0298; B and C → 0.02/0.0298
  close(imp.A.fussellVesely, 0.01 / 0.0298, 1e-12);
  close(imp.B.fussellVesely, 0.02 / 0.0298, 1e-12);
  close(imp.A.raw, 1 / 0.0298, 1e-12);
  close(imp.A.rrw, 0.0298 / 0.02, 1e-12);
  assert.equal(r.cutSets[0].events.join(), 'B,C');
});

test('2-out-of-3 voting gate: 3p² − 2p³, cut sets are the three pairs', () => {
  const p = 0.1;
  const tree = { type: 'vote', k: 2, children: [B('X1', p), B('X2', p), B('X3', p)] };
  const r = evaluateFaultTree(tree);
  close(r.probability, 3 * p * p - 2 * p ** 3, 1e-12);
  assert.deepEqual(minimalCutSets(tree), [['X1', 'X2'], ['X1', 'X3'], ['X2', 'X3']]);
  // Non-identical components, by enumeration: p1p2 + p1p3 + p2p3 − 2p1p2p3
  const t2 = { type: 'vote', k: 2, children: [B('X1', 0.1), B('X2', 0.2), B('X3', 0.3)] };
  close(evaluateFaultTree(t2).probability, 0.02 + 0.03 + 0.06 - 2 * 0.006, 1e-12);
  // 3-out-of-5 identical: Σ_{k>=3} C(5,k) p^k (1−p)^(5−k)
  const t5 = { type: 'vote', k: 3, children: [1, 2, 3, 4, 5].map((i) => B(`Y${i}`, 0.2)) };
  close(evaluateFaultTree(t5).probability, 10 * 0.008 * 0.64 + 5 * 0.0016 * 0.8 + 0.00032, 1e-12);
  assert.equal(minimalCutSets(t5).length, 10);
});

test('exact BDD probability against brute-force enumeration on a tree with shared events and nested gates', () => {
  const p = { A: 0.05, B: 0.1, C: 0.15, D: 0.2, E: 0.25, F: 0.3 };
  const b = (id) => ({ type: 'basic', id });
  const tree = {
    type: 'or', children: [
      { type: 'and', children: [b('A'), { type: 'or', children: [b('B'), b('C')] }] },
      { type: 'vote', k: 2, children: [b('C'), b('D'), { type: 'and', children: [b('E'), b('A')] }] },
      { type: 'and', children: [b('F'), b('B'), b('D')] },
    ],
  };
  const events = Object.fromEntries(Object.entries(p).map(([k, v]) => [k, { p: v }]));
  const r = evaluateFaultTree(tree, { events });
  const ids = Object.keys(p);
  const evalNode = (n, s) => {
    if (n.type === 'basic') return s[n.id];
    const v = n.children.map((c) => evalNode(c, s));
    if (n.type === 'and') return v.every(Boolean);
    if (n.type === 'or') return v.some(Boolean);
    return v.filter(Boolean).length >= n.k;
  };
  let exact = 0;
  for (let mask = 0; mask < 64; mask++) {
    const s = {};
    let pr = 1;
    ids.forEach((id, i) => { s[id] = !!(mask & (1 << i)); pr *= s[id] ? p[id] : 1 - p[id]; });
    if (evalNode(tree, s)) exact += pr;
  }
  close(r.probability, exact, 1e-12);
  close(inclusionExclusion(minimalCutSets(tree), p), exact, 1e-12);
  // Every cut set is minimal: removing any event stops it from failing the top event.
  for (const cs of minimalCutSets(tree)) {
    const s = Object.fromEntries(ids.map((id) => [id, cs.includes(id)]));
    assert.equal(evalNode(tree, s), true);
    for (const id of cs) assert.equal(evalNode(tree, { ...s, [id]: false }), false);
  }
});

test('NOT gate (non-coherent tree): P(A and not B) = pA·(1 − pB); cut sets are refused', () => {
  const tree = { type: 'and', children: [B('A', 0.3), { type: 'not', child: B('B', 0.4) }] };
  const r = evaluateFaultTree(tree);
  close(r.probability, 0.18, 1e-12);
  assert.equal(r.coherent, false);
  assert.equal(r.cutSets, null);
  assert.throws(() => minimalCutSets(tree), /coherent/);
  // Exclusive or built from gates: pA(1 − pB) + pB(1 − pA)
  const xor = { type: 'or', children: [tree, { type: 'and', children: [B('B', 0.4), { type: 'not', child: B('A', 0.3) }] }] };
  close(evaluateFaultTree(xor).probability, 0.3 * 0.6 + 0.4 * 0.7, 1e-12);
});

test('event tree: outcome frequencies and consequence-weighted risk by hand', () => {
  // Leak 1e-2/yr; ignition 0.1; if ignited, deluge fails 0.2.
  const tree = {
    name: 'ignition', branches: [
      { label: 'ignited', probability: 0.1, next: { name: 'deluge', branches: [
        { label: 'works', probability: 0.8, outcome: 'controlled fire', consequence: 1e6 },
        { label: 'fails', outcome: 'escalated fire', consequence: 5e7 }] } },
      { label: 'not ignited', probability: 0.9, outcome: 'release', consequence: 1e4 },
    ],
  };
  const r = evaluateEventTree({ initiatingFrequency: 1e-2, tree });
  close(r.outcomes['controlled fire'].frequency, 8e-4, 1e-12);
  close(r.outcomes['escalated fire'].frequency, 2e-4, 1e-12);
  close(r.outcomes.release.frequency, 9e-3, 1e-12);
  close(r.totalRisk, 800 + 10000 + 90, 1e-12);
  close(r.sequences.reduce((s, q) => s + q.probability, 0), 1, 1e-12);
  assert.deepEqual(r.sequences[1].path, ['ignited', 'fails']);
  assert.throws(() => evaluateEventTree({ initiatingFrequency: 1, tree: { name: 'x', branches: [{ probability: 0.5, outcome: 'a' }, { probability: 0.6, outcome: 'b' }] } }), /sum to 1/);
});

test('FMECA: RPN and MIL-STD-1629A criticality numbers and ranking', () => {
  const r = fmeca([
    { item: 'pump', mode: 'seal leak', severity: 6, occurrence: 5, detection: 4, failureRate: 2e-5, modeRatio: 0.6, effectProbability: 0.5, operatingTime: 8760 },
    { item: 'pump', mode: 'bearing seizure', severity: 9, occurrence: 2, detection: 6, failureRate: 2e-5, modeRatio: 0.4, effectProbability: 1, operatingTime: 8760 },
    { item: 'valve', mode: 'fails closed', severity: 8, occurrence: 3, detection: 2 },
  ]);
  assert.equal(r.rows[0].rpn, 120);
  assert.equal(r.rows[1].rpn, 108);
  assert.equal(r.rows[2].rpn, 48);
  assert.equal(r.rows[0].rpnRank, 1);
  // C_m = β·α·λ·t: 0.5·0.6·2e-5·8760 = 0.05256; 1·0.4·2e-5·8760 = 0.07008
  close(r.rows[0].criticality, 0.05256, 1e-12);
  close(r.rows[1].criticality, 0.07008, 1e-12);
  assert.equal(r.rows[1].criticalityRank, 1);
  assert.equal(r.rows[2].criticality, null);
  close(r.itemCriticality.pump, 0.12264, 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => evaluateFaultTree({ type: 'or', children: [{ type: 'basic', id: 'A' }] }), /no probability model/);
  assert.throws(() => evaluateFaultTree({ type: 'nand', children: [B('A', 0.1)] }), /unknown gate/);
  assert.throws(() => basicEventProbability({ p: 1.5 }), /p must/);
});
