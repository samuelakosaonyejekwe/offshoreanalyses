import test from 'node:test';
import assert from 'node:assert/strict';
import * as sch from '../../src/engine/econ/schedule.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

// Six-activity textbook network: A(3) -> B(4), C(2); B -> D(5); C -> E(1); D, E -> F(2). Longest path A-B-D-F = 14.
const net = [
  { id: 'A', duration: 3 }, { id: 'B', duration: 4, predecessors: ['A'] }, { id: 'C', duration: 2, predecessors: ['A'] },
  { id: 'D', duration: 5, predecessors: ['B'] }, { id: 'E', duration: 1, predecessors: ['C'] }, { id: 'F', duration: 2, predecessors: ['D', 'E'] },
];

test('CPM on a small textbook network: duration 14, critical path A-B-D-F, float 6 on C and E', () => {
  const r = sch.cpm(net);
  assert.equal(r.duration, 14);
  assert.deepEqual(r.criticalPath, ['A', 'B', 'D', 'F']);
  const by = Object.fromEntries(r.activities.map((a) => [a.id, a]));
  assert.deepEqual([by.C.ES, by.C.EF, by.C.LS, by.C.LF], [3, 5, 9, 11]);
  assert.equal(by.C.totalFloat, 6); assert.equal(by.C.freeFloat, 0);
  assert.equal(by.E.totalFloat, 6); assert.equal(by.E.freeFloat, 6);
  assert.deepEqual([by.D.ES, by.D.LF], [7, 12]);
  for (const id of ['A', 'B', 'D', 'F']) assert.equal(by[id].critical, true);
});

test('CPM with lags and start-to-start links; cycles and unknown predecessors are rejected', () => {
  const r = sch.cpm([{ id: 'A', duration: 5 }, { id: 'B', duration: 3, predecessors: [{ id: 'A', type: 'SS', lag: 2 }] },
    { id: 'C', duration: 2, predecessors: [{ id: 'A', lag: 4 }] }]);
  const by = Object.fromEntries(r.activities.map((a) => [a.id, a]));
  assert.equal(by.B.ES, 2); assert.equal(by.C.ES, 9); assert.equal(r.duration, 11);
  assert.throws(() => sch.cpm([{ id: 'A', duration: 1, predecessors: ['B'] }, { id: 'B', duration: 1, predecessors: ['A'] }]), /cycle/);
  assert.throws(() => sch.cpm([{ id: 'A', duration: 1, predecessors: ['Z'] }]), /unknown predecessor/);
});

test('PERT three-point estimate (2, 4, 12): mean 5, standard deviation 5/3', () => {
  const p = sch.pertEstimate({ a: 2, m: 4, b: 12 });
  close(p.mean, 5); close(p.std, 10 / 6, 1e-14);
  const n = sch.pertNetwork([{ id: 'A', a: 2, m: 4, b: 12 }, { id: 'B', a: 3, m: 3, b: 3, predecessors: ['A'] }], 8);
  close(n.expectedDuration, 8); close(n.std, 10 / 6, 1e-14); close(n.probabilityByTarget, 0.5, 1e-12);
  close(sch.pertNetwork([{ id: 'A', a: 2, m: 4, b: 12 }], 5 + 10 / 6).probabilityByTarget, 0.8413447460685429, 1e-12);
});

test('schedule risk: reproducible for a seed, exact for fixed durations, P80 above P50', () => {
  const fixed = sch.scheduleRisk({ activities: net, n: 50 });
  for (const d of fixed.durations) assert.equal(d, 14);
  assert.equal(fixed.criticalityIndex.A, 1); assert.equal(fixed.criticalityIndex.C, 0);
  const risky = net.map((a) => ({ id: a.id, predecessors: a.predecessors, a: a.duration * 0.8, m: a.duration, b: a.duration * 1.6 }));
  const r1 = sch.scheduleRisk({ activities: risky, n: 1500, seed: 11 }); const r2 = sch.scheduleRisk({ activities: risky, n: 1500, seed: 11 });
  assert.deepEqual(Array.from(r1.durations), Array.from(r2.durations));
  assert.notEqual(sch.scheduleRisk({ activities: risky, n: 1500, seed: 12 }).mean, r1.mean);
  assert.ok(r1.p80 > r1.p50 && r1.p50 > r1.p10 && r1.p90 > r1.p80);
  // PERT mean of each activity is d (0.8 + 4 + 1.6)/6; the critical path A-B-D-F sums to 14 x 1.0667
  close(r1.mean, 14 * 6.4 / 6, 0.01);
  assert.equal(r1.criticalityIndex.A, 1); assert.ok(r1.criticalityIndex.C < 0.01);
  assert.match(r1.convention, /not exceeded/);
});

test('schedule risk: two identical parallel paths share criticality and the mean exceeds the deterministic duration (merge bias)', () => {
  const acts = [{ id: 'P', dist: { dist: 'normal', mean: 10, std: 2 } }, { id: 'Q', dist: { dist: 'normal', mean: 10, std: 2 } },
    { id: 'End', duration: 0, predecessors: ['P', 'Q'] }];
  const r = sch.scheduleRisk({ activities: acts, n: 4000, seed: 5 });
  close(r.criticalityIndex.P, 0.5, 0.05); close(r.criticalityIndex.P + r.criticalityIndex.Q, 1, 1e-9);
  // E[max of two iid normals] = mu + sigma / sqrt(pi)
  close(r.mean, 10 + 2 / Math.sqrt(Math.PI), 0.01);
  assert.equal(r.deterministic, 10);
});

test('weather adjustment: constant operability 0.5 doubles the duration; seasonal operability is integrated month by month', () => {
  close(sch.weatherAdjustedDuration({ netDays: 10, operabilityByMonth: new Array(12).fill(0.5) }).grossDays, 20, 1e-12);
  const month = 365.25 / 12; const op = new Array(12).fill(1); op[0] = 0.5;
  // start 1 January, 20 net days: January yields 0.5 x 30.44 = 15.22 net days, the remaining 4.78 take 4.78 days in February
  const r = sch.weatherAdjustedDuration({ netDays: 20, operabilityByMonth: op });
  close(r.grossDays, month + (20 - 0.5 * month), 1e-12);
  const hooked = sch.scheduleRisk({ activities: [{ id: 'lay', duration: 10, operabilityByMonth: new Array(12).fill(0.5) }], n: 5 });
  close(hooked.mean, 20, 1e-12);
});

test('cost loading: uniform spread over activity dates, S-curve ends at the total, fractions sum to 1', () => {
  const s = sch.cpm([{ id: 'A', duration: 10 }, { id: 'B', duration: 10, predecessors: ['A'] }]);
  const c = sch.costLoading({ schedule: s, costs: { A: 100, B: 300 }, periodDays: 5 });
  assert.deepEqual(c.cost, [50, 50, 150, 150]); assert.deepEqual(c.cumulative, [50, 100, 250, 400]);
  close(c.fractions.reduce((a, v) => a + v, 0), 1, 1e-14);
  const full = sch.cpm(net);
  const early = sch.costLoading({ schedule: full, costs: { A: 10, B: 10, C: 40, D: 10, E: 10, F: 10 }, periodDays: 7 });
  const late = sch.costLoading({ schedule: full, costs: { A: 10, B: 10, C: 40, D: 10, E: 10, F: 10 }, periodDays: 7, timing: 'late' });
  close(early.total, 90); close(early.cumulative[1], 90, 1e-12); close(late.cumulative[1], 90, 1e-12);
  assert.ok(early.cumulative[0] > late.cumulative[0]);      // early-start curve leads the late-start curve
});

test('delay cost and liquidated damages by hand', () => {
  const d = sch.delayCost({ delayDays: 30, standbyDayRates: [200000, 50000], revenuePerDay: 400000 });
  close(d.standby, 7.5e6); close(d.revenueImpact, 1.2e7); close(d.total, 1.95e7);
  const def = sch.delayCost({ delayDays: 365.25, mode: 'deferred', pvOfOperations: 1e9, discountRate: 0.08 });
  close(def.revenueImpact, 1e9 * (1 - 1 / 1.08), 1e-12);
  const ld = sch.liquidatedDamages({ delayDays: 60, ratePerDay: 100000, graceDays: 10, contractValue: 4e7, capFraction: 0.1 });
  close(ld.damages, 4e6); assert.equal(ld.capped, true); assert.equal(ld.chargeableDays, 50);
  close(sch.liquidatedDamages({ delayDays: 20, ratePerDay: 100000, graceDays: 10 }).damages, 1e6);
});
