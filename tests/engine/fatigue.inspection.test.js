import test from 'node:test';
import assert from 'node:assert/strict';
import {
  pod, podPreset, crackPopulation, crackSizesAt, updateAfterInspection, detectionProbability,
  populationFailureProbability, planInspections, valueOfInformation,
} from '../../src/engine/fatigue/inspection.js';
import { parisConstantToSI, parisCrackSizeClosedForm } from '../../src/engine/fatigue/fracture.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('POD curve: POD(X0) = 0.5 and hand values of 1 − 1/(1 + (a/X0)^b)', () => {
  assert.equal(pod(0.45e-3, { x0: 0.45e-3, b: 0.9 }), 0.5);
  // a = 2 mm, X0 = 0.4 mm, b = 1.43: 5^1.43 = 9.989 → POD = 1 − 1/10.989 = 0.9090
  close(pod(2e-3, { x0: 0.4e-3, b: 1.43 }), 1 - 1 / (1 + 5 ** 1.43), 1e-14);
  close(pod(2e-3, { x0: 0.4e-3, b: 1.43 }), 0.9090, 1e-3);
  assert.equal(pod(0, { x0: 1e-3, b: 1 }), 0);
  assert.equal(pod(Infinity, { x0: 1e-3, b: 1 }), 1);
  assert.equal(podPreset('visual-easy').measure, 'length');
  assert.ok(pod(1e-3, podPreset('ec-good')) > pod(1e-3, podPreset('ec-underwater')));
});

const model = {
  a0: { dist: 'exponential', mean: 0.11e-3 },
  C: parisConstantToSI(5.21e-13, 3),
  m: 3, stressRange: 30e6, Y: 1.12, cyclesPerYear: 5e6, aCrit: 0.02, n: 40000, seed: 3,
};

test('prior failure probability of the crack population against the analytic result (only a0 random, exponential)', () => {
  const pop = crackPopulation(model);
  // Failure by time t ⇔ a0 >= a*(t), with a*(t) from the Paris closed form run backwards:
  // a*^(−1/2) = aCrit^(−1/2) + 0.5·k·t, k = C·(Y·S·√π)^3·ν. P_f = exp(−a*/mean).
  const k = model.C * (1.12 * 30e6 * Math.sqrt(Math.PI)) ** 3 * 5e6;
  const times = [10, 20, 30];
  const r = populationFailureProbability(pop, times);
  times.forEach((t, i) => {
    const aStar = (0.02 ** -0.5 + 0.5 * k * t) ** -2;
    const exact = Math.exp(-aStar / 0.11e-3);
    const se = Math.sqrt((exact * (1 - exact)) / model.n);
    assert.ok(Math.abs(r.pf[i] - exact) < 4 * se + 1e-12, `t=${t}: ${r.pf[i]} vs ${exact}`);
  });
  // Crack sizes follow the closed form.
  const a = crackSizesAt(pop, 5);
  close(a[0], parisCrackSizeClosedForm({ a0: pop.a0[0], cycles: 5 * 5e6, stressRange: 30e6, Y: 1.12, A: model.C, m: 3 }), 1e-10);
});

test('Bayesian updating: no detection lowers, a large measured crack raises the failure probability; Bayes identity holds', () => {
  const pop = crackPopulation({ ...model, C: { dist: 'lognormal', mean: model.C, cov: 0.5 }, stressFactor: { dist: 'lognormal', mean: 1, cov: 0.15 } });
  const podPar = { x0: 0.45e-3, b: 0.9 };
  const prior = populationFailureProbability(pop, [25], { survivedTo: 10 }).pf[0];
  const clean = updateAfterInspection(pop, { time: 10, pod: podPar, detected: false });
  const postClean = populationFailureProbability(clean, [25], { survivedTo: 10 }).pf[0];
  const found = updateAfterInspection(pop, { time: 10, pod: podPar, detected: true });
  const postFound = populationFailureProbability(found, [25], { survivedTo: 10 }).pf[0];
  assert.ok(postClean < prior, `${postClean} < ${prior}`);
  assert.ok(postFound > prior, `${postFound} > ${prior}`);
  // Law of total probability over the two inspection outcomes.
  const pDet = detectionProbability(pop, 10, podPar);
  close(pDet * postFound + (1 - pDet) * postClean, prior, 1e-9);
  // A measured crack of 5 mm concentrates the posterior around that size.
  const sized = updateAfterInspection(pop, { time: 10, pod: podPar, detected: true, measuredSize: 5e-3, sizingStd: 0.5e-3 });
  const a = crackSizesAt(pop, 10);
  let sw = 0;
  let sa = 0;
  for (let i = 0; i < pop.n; i++) { if (sized.weights[i] > 0) { sw += sized.weights[i]; sa += sized.weights[i] * a[i]; } }
  close(sa / sw, 5e-3, 0.1);
  assert.ok(populationFailureProbability(sized, [25], { survivedTo: 10 }).pf[0] > postFound);
  // The original population is not modified.
  assert.equal(pop.weights[0], 1);
});

test('risk-based inspection plan keeps the annual failure probability at or below the target', () => {
  const pop = crackPopulation({ ...model, C: { dist: 'lognormal', mean: model.C, cov: 0.5 }, stressFactor: { dist: 'lognormal', mean: 1, cov: 0.15 }, n: 60000 });
  const none = planInspections({ population: pop, pod: { x0: 0.45e-3, b: 0.9 }, targetAnnualPf: 0.5, serviceLife: 30 });
  assert.equal(none.inspectionTimes.length, 0);
  const target = none.maxAnnualPf / 3;
  const plan = planInspections({ population: pop, pod: { x0: 0.45e-3, b: 0.9 }, targetAnnualPf: target, serviceLife: 30 });
  assert.ok(plan.inspectionTimes.length >= 1);
  assert.equal(plan.table.length, 30);
  assert.equal(plan.table.filter((r) => r.inspection).length, plan.inspectionTimes.length);
  assert.ok(plan.maxAnnualPf < none.maxAnnualPf);
  if (plan.meetsTarget) assert.ok(plan.maxAnnualPf <= target);
  // A poorer inspection method needs at least as many inspections.
  const poor = planInspections({ population: pop, pod: { x0: 5e-3, b: 0.9 }, targetAnnualPf: target, serviceLife: 30 });
  assert.ok(poor.inspectionTimes.length >= plan.inspectionTimes.length);
});

test('value of information: hand calculation of the pre-posterior decision analysis', () => {
  // Prior: P(damage) = 0.1, P(fail | damage) = 0.5, failure cost 1000, repair 60 → do nothing costs 50, repair 60.
  // Inspection: POD 0.9, false alarm 0.05 → P(indication) = 0.09 + 0.045 = 0.135,
  // P(damage | indication) = 0.6667 → min(333.3, 60) = 60; P(damage | clear) = 0.01/0.865 = 0.011561 → 5.780.
  // Cost with inspection = 0.135·60 + 0.865·5.780 = 8.1 + 5.0 = 13.1; VoI = 50 − 13.1 = 36.9.
  const r = valueOfInformation({ pDamage: 0.1, pFailGivenDamage: 0.5, pod: 0.9, falseAlarm: 0.05, costFailure: 1000, costRepair: 60, costInspection: 10 });
  close(r.costWithout, 50, 1e-12);
  assert.equal(r.actionWithout, 'do-nothing');
  close(r.posteriorDamageGivenIndication, 0.09 / 0.135, 1e-12);
  close(r.costWith, 13.1, 1e-12);
  close(r.valueOfInformation, 36.9, 1e-12);
  close(r.netValue, 26.9, 1e-12);
  assert.equal(r.worthwhile, true);
  // A useless inspection (POD = false-alarm rate) has no value.
  const u = valueOfInformation({ pDamage: 0.1, pFailGivenDamage: 0.5, pod: 0.3, falseAlarm: 0.3, costFailure: 1000, costRepair: 60 });
  assert.ok(Math.abs(u.valueOfInformation) < 1e-9);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => pod(1, { x0: 0, b: 1 }), /x0/);
  assert.throws(() => crackPopulation({ ...model, m: -1 }), /m must/);
  assert.throws(() => valueOfInformation({ pDamage: 1.2, pFailGivenDamage: 0.5, pod: 0.9, costFailure: 1, costRepair: 1 }), /pDamage/);
});
