import test from 'node:test';
import assert from 'node:assert/strict';
import * as Wn from '../../src/engine/metocean/windows.js';
import { makeRng } from '../../src/engine/math/rng.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

// 3-hourly hand-built record, 16 steps (48 h)
const hs = [1.0, 1.2, 2.6, 2.8, 1.4, 1.1, 0.9, 1.3, 1.6, 3.0, 3.2, 1.9, 1.0, 0.8, 0.7, 2.4];
const wind = [5, 6, 9, 12, 8, 7, 11, 6, 5, 14, 15, 9, 6, 5, 4, 8];
const tp = [7, 7, 9, 10, 8, 8, 7, 13, 8, 11, 11, 9, 8, 7, 7, 9];

test('workability mask with several simultaneous limits (Hs AND wind AND Tp range), NaN is not workable', () => {
  const mask = Wn.workabilityMask({ hs, wind, tp }, [{ key: 'hs', max: 2.0 }, { key: 'wind', max: 10 }, { key: 'tp', min: 6, max: 12 }]);
  //                         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  assert.deepEqual(Array.from(mask), [1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0]);
  const m2 = Wn.workabilityMask({ hs: [1, NaN, 1] }, [{ key: 'hs', max: 2 }]);
  assert.deepEqual(Array.from(m2), [1, 0, 1]);
  assert.throws(() => Wn.workabilityMask({ hs }, [{ key: 'nope', max: 2 }]), RangeError);
  assert.throws(() => Wn.workabilityMask({ hs, wind: [1, 2] }, [{ key: 'hs', max: 2 }, { key: 'wind', max: 2 }]), RangeError);
});

test('alpha factor reduces the forecast limit (OPWF = α · OPLIM)', () => {
  assert.ok(rel(Wn.forecastLimit(2.5, 0.76), 1.9) < 1e-15);
  const full = Wn.workabilityMask({ hs }, [{ key: 'hs', max: 2.0 }]);
  const red = Wn.workabilityMask({ hs }, [{ key: 'hs', max: 2.0, alpha: 0.6 }]); // forecast limit 1.2 m
  assert.equal(full.reduce((a, b) => a + b), 11);
  assert.deepEqual(Array.from(red), [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0]);
  assert.throws(() => Wn.forecastLimit(2, 1.2), RangeError);
});

test('operation reference period TR = TPOP + TC and weather-restricted classification', () => {
  const a = Wn.operationReferencePeriod({ plannedHours: 24, contingencyHours: 12 });
  assert.equal(a.TR, 36); assert.equal(a.weatherRestricted, true);
  assert.equal(Wn.operationReferencePeriod({ plannedHours: 60, contingencyHours: 48 }).weatherRestricted, false); // TR = 108 h
  assert.equal(Wn.operationReferencePeriod({ plannedHours: 80, contingencyHours: 10 }).weatherRestricted, false); // TPOP > 72 h
});

test('weather windows on the hand-built series: list, durations and minimum length filter', () => {
  const mask = Wn.workabilityMask({ hs }, [{ key: 'hs', max: 2.0 }]);
  //                                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  assert.deepEqual(Array.from(mask), [1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0]);
  const w = Wn.findWindows(mask, 3);
  assert.deepEqual(w, [
    { startIndex: 0, endIndex: 1, startHours: 0, durationHours: 6 },
    { startIndex: 4, endIndex: 8, startHours: 12, durationHours: 15 },
    { startIndex: 11, endIndex: 14, startHours: 33, durationHours: 12 },
  ]);
  assert.equal(Wn.findWindows(mask, 3, 12).length, 2);
  assert.equal(Wn.findWindows(mask, 3, 13).length, 1);
});

test('operability: workable and startable fractions and waiting time, hand calculation', () => {
  // required 9 h = 3 steps. Startable steps: 4, 5, 6 (window 4..8) and 11, 12 (window 11..14) -> 5/16
  const r = Wn.operability({ series: { hs }, limits: [{ key: 'hs', max: 2.0 }], dtHours: 3, requiredHours: 9 });
  assert.equal(r.requiredSteps, 3);
  assert.ok(rel(r.workableFraction, 11 / 16) < 1e-15);
  assert.ok(rel(r.startableFraction, 5 / 16) < 1e-15);
  assert.equal(r.windows.length, 2);
  // waiting steps: i=0..3 -> 4,3,2,1 ; 4,5,6 -> 0 ; 7..10 -> 4,3,2,1 ; 11,12 -> 0 ; 13..15 no later window (excluded)
  // mean = (10 + 10) / 13 steps = 4.615 h ; max = 12 h
  assert.ok(rel(r.meanWaitingHours, (20 / 13) * 3) < 1e-14);
  assert.equal(r.maxWaitingHours, 12);
  // with requiredHours omitted every workable step is startable
  const r1 = Wn.operability({ series: { hs }, limits: [{ key: 'hs', max: 2.0 }], dtHours: 3 });
  assert.equal(r1.startableFraction, r1.workableFraction);
  // nothing workable
  const r0 = Wn.operability({ series: { hs }, limits: [{ key: 'hs', max: 0.1 }], dtHours: 3 });
  assert.equal(r0.workableFraction, 0); assert.ok(Number.isNaN(r0.meanWaitingHours)); assert.equal(r0.windows.length, 0);
});

test('operability: monthly and seasonal percentages from time stamps', () => {
  // daily samples over 2023: January calm on odd days only, July always calm, everything else rough
  const n = 365; const t = new Array(n); const h = new Array(n);
  const t0 = Date.UTC(2023, 0, 1, 12);
  for (let i = 0; i < n; i++) {
    t[i] = t0 + i * 86400000;
    const dte = new Date(t[i]);
    const mo = dte.getUTCMonth();
    h[i] = mo === 6 ? 0.5 : mo === 0 && dte.getUTCDate() % 2 === 1 ? 0.5 : 3;
  }
  const r = Wn.operability({ series: { hs: h }, limits: [{ key: 'hs', max: 1.5 }], dtHours: 24, requiredHours: 48, timeMs: t });
  assert.equal(r.monthly.length, 12);
  assert.equal(r.monthly[0].n, 31);
  assert.ok(rel(r.monthly[0].workableFraction, 16 / 31) < 1e-15);
  assert.equal(r.monthly[0].startableFraction, 0); // never two calm days in a row
  assert.equal(r.monthly[6].workableFraction, 1);
  assert.ok(rel(r.monthly[6].startableFraction, 30 / 31) < 1e-15); // the last day of July has only one day left
  assert.equal(r.monthly[3].workableFraction, 0);
  const jja = r.seasonal.find((s) => s.season === 'JJA');
  assert.ok(rel(jja.workableFraction, 31 / 92) < 1e-15);
  assert.ok(rel(r.seasonal.find((s) => s.season === 'DJF').workableFraction, 16 / 90) < 1e-15);
  // waiting in June: from 1 June (30 days to 1 July) down to 1 day -> mean 15.5 days
  assert.ok(rel(r.monthly[5].meanWaitingHours, 15.5 * 24) < 1e-14);
});

test('persistence table of calms and storms on the hand-built series', () => {
  const rows = Wn.persistenceTable({ values: hs, dtHours: 3, thresholds: [2.0, 1.0], durations: [6, 12, 15, 18] });
  const a = rows[0]; // calms Hs <= 2.0: 6 h, 15 h, 12 h
  assert.equal(a.nEvents, 3);
  assert.ok(rel(a.occurrence, 11 / 16) < 1e-15);
  assert.ok(rel(a.meanDurationHours, 11) < 1e-15);
  assert.equal(a.maxDurationHours, 15);
  assert.deepEqual(a.eventExceedance, [1, 2 / 3, 1 / 3, 0]);
  assert.deepEqual(a.timeFraction, [33 / 48, 27 / 48, 15 / 48, 0]);
  const b = rows[1]; // Hs <= 1.0: steps 0, 6, 12-14 -> 3 h, 3 h, 9 h
  assert.equal(b.nEvents, 3); assert.equal(b.maxDurationHours, 9);
  const st = Wn.persistenceTable({ values: hs, dtHours: 3, thresholds: [2.0], durations: [3, 6], sense: 'above' })[0];
  // storms Hs > 2.0: steps 2-3 (6 h), 9-10 (6 h), 15 (3 h)
  assert.equal(st.nEvents, 3);
  assert.deepEqual(st.eventExceedance, [1, 2 / 3]);
  assert.ok(rel(st.occurrence + a.occurrence, 1) < 1e-15);
});

test('Markov chain: stationary distribution of a two-state chain is [b, a]/(a + b); periodic chain handled', () => {
  const a = 0.2; const b = 0.5;
  const pi = Wn.markovStationary([[1 - a, a], [b, 1 - b]]);
  assert.ok(Math.abs(pi[0] - b / (a + b)) < 1e-11 && Math.abs(pi[1] - a / (a + b)) < 1e-11);
  const per = Wn.markovStationary([[0, 1], [1, 0]]);
  assert.ok(Math.abs(per[0] - 0.5) < 1e-11);
});

test('Markov chain: transition counts on a hand-built sequence', () => {
  // states by Hs bins [0,1), [1,2), [2,∞): sequence 0 0 1 2 1 0 0 1
  const m = Wn.markovFit([0.5, 0.6, 1.5, 2.5, 1.2, 0.3, 0.2, 1.9], [0, 1, 2, 10]);
  assert.deepEqual(m.counts, [[2, 2, 0], [1, 0, 1], [0, 1, 0]]);
  assert.deepEqual(m.P, [[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 1, 0]]);
  // stationary: π = π P -> π = [0.4, 0.4, 0.2]
  for (const [i, v] of [0.4, 0.4, 0.2].entries()) assert.ok(Math.abs(m.stationary[i] - v) < 1e-10);
  assert.throws(() => Wn.markovFit([1, 2, 3], [0, 0]), RangeError);
});

test('Markov chain: simulate from a known matrix, refit, and recover matrix and stationary distribution', () => {
  const P = [[0.90, 0.09, 0.01], [0.15, 0.75, 0.10], [0.02, 0.28, 0.70]];
  const edges = [0, 1.5, 3, 6];
  const truth = { edges, P, stationary: Wn.markovStationary(P) };
  // synthetic record: 30 years of 3-hourly sea states
  const sim = Wn.markovSimulate(truth, { nSteps: 30 * 2922, seed: 11 });
  const fit = Wn.markovFit(sim.values, edges);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) assert.ok(Math.abs(fit.P[i][j] - P[i][j]) < 0.012, `P[${i}][${j}] = ${fit.P[i][j]}`);
  for (let i = 0; i < 3; i++) {
    assert.ok(Math.abs(fit.stationary[i] - truth.stationary[i]) < 0.02);
    assert.ok(Math.abs(fit.occupancy[i] - truth.stationary[i]) < 0.02);
    const s = fit.P[i].reduce((x, y) => x + y);
    assert.ok(Math.abs(s - 1) < 1e-12);
  }
  // analytic check of the stationary vector: π P = π
  for (let j = 0; j < 3; j++) {
    let v = 0;
    for (let i = 0; i < 3; i++) v += truth.stationary[i] * P[i][j];
    assert.ok(Math.abs(v - truth.stationary[j]) < 1e-11);
  }
  // mean sojourn time in state 0 is dt/(1 − P00) = 10 steps
  const calm = Wn.persistenceTable({ values: sim.values, dtHours: 3, thresholds: [1.5 - 1e-9], durations: [] })[0];
  assert.ok(rel(calm.meanDurationHours, 30) < 0.05);
  // deterministic for a seed, values inside their bins
  const again = Wn.markovSimulate(truth, { nSteps: 500, seed: 11 });
  assert.deepEqual(Array.from(again.states), Array.from(sim.states.slice(0, 500)));
  for (let i = 0; i < 500; i++) assert.ok(again.values[i] >= edges[again.states[i]] && again.values[i] < edges[again.states[i] + 1]);
  assert.equal(Wn.markovSimulate(truth, { nSteps: 3, seed: 1, startState: 2, sample: 'mid' }).values[0], 4.5);
});

test('campaign Monte Carlo: always-workable weather gives zero downtime; periodic weather gives the hand result', () => {
  const calm = { hs: new Array(240).fill(0.5) };
  const ops = [
    { name: 'transit', limits: [{ key: 'hs', max: 3 }], netHours: 12 },
    { name: 'lift', limits: [{ key: 'hs', max: 1.5, alpha: 0.8 }], netHours: 6, contingencyHours: 6 },
  ];
  const r0 = Wn.campaignMonteCarlo({ series: calm, dtHours: 3, operations: ops, nTrials: 50, seed: 3 });
  assert.equal(r0.netHours, 18);
  assert.equal(r0.completionHours.p10, 18); assert.equal(r0.completionHours.p90, 18); assert.equal(r0.downtimeHours.mean, 0);
  // periodic weather with a 24 h cycle at 3 h steps: 4 calm steps (12 h) then 4 rough steps (12 h).
  // one operation: net 6 h + 3 h contingency = 3 steps. Startable phases 0 and 1 (remaining >= 3);
  // waiting steps by phase 0..7: 0, 0, 6, 5, 4, 3, 2, 1  -> mean 21/8 steps = 7.875 h
  const cyc = { hs: Array.from({ length: 8 * 40 }, (_, i) => (i % 8 < 4 ? 0.5 : 3)) };
  const one = [{ limits: [{ key: 'hs', max: 1.5 }], netHours: 6, contingencyHours: 3 }];
  const r = Wn.campaignMonteCarlo({ series: cyc, dtHours: 3, operations: one, nTrials: 8000, seed: 5 });
  assert.ok(rel(r.downtimeHours.mean, 7.875) < 0.03);
  assert.ok(rel(r.operations[0].startableFraction, 0.25) < 1e-15);
  assert.equal(r.completionHours.min, 6); assert.equal(r.completionHours.max, 6 + 18);
  assert.ok(r.completionHours.p10 <= r.completionHours.p50 && r.completionHours.p50 <= r.completionHours.p90);
  // two such operations in sequence: the first ends 2 steps after its start phase (0 or 1), i.e. at phase 2 or 3,
  // so the second always waits 6 or 5 steps
  const two = Wn.campaignMonteCarlo({ series: cyc, dtHours: 3, operations: [one[0], one[0]], nTrials: 2000, seed: 5 });
  assert.ok(two.operations[1].meanWaitHours >= 15 && two.operations[1].meanWaitHours <= 18);
  assert.ok(two.completionHours.min >= 12 + 15);
  // same seed -> same samples; start month restriction and impossible operations
  const again = Wn.campaignMonteCarlo({ series: cyc, dtHours: 3, operations: one, nTrials: 8000, seed: 5 });
  assert.deepEqual(Array.from(again.samples.slice(0, 100)), Array.from(r.samples.slice(0, 100)));
  assert.throws(() => Wn.campaignMonteCarlo({ series: cyc, dtHours: 3, operations: [{ limits: [{ key: 'hs', max: 1.5 }], netHours: 24 }] }), Error);
});

test('campaign Monte Carlo: start months restrict the draw; wrap-around of the record is periodic', () => {
  // one year of daily data: workable only in the first 10 days of each month; operation needs 5 days
  const n = 365; const t = new Array(n); const h = new Array(n);
  for (let i = 0; i < n; i++) { t[i] = Date.UTC(2023, 0, 1 + i, 12); h[i] = new Date(t[i]).getUTCDate() <= 10 ? 0.5 : 3; }
  const ops = [{ limits: [{ key: 'hs', max: 1 }], netHours: 120 }];
  const dec = Wn.campaignMonteCarlo({ series: { hs: h }, dtHours: 24, operations: ops, nTrials: 400, seed: 2, timeMs: t, startMonths: [12] });
  // starting on 1–6 Dec: no wait; 7–31 Dec: wait until 1 Jan of the wrapped record (25 … 1 days)
  // mean wait = (0·6 + Σ_{k=1}^{25} k) / 31 = 325/31 days
  assert.ok(rel(dec.downtimeHours.mean, (325 / 31) * 24) < 0.08);
  assert.equal(dec.completionHours.max, 120 + 25 * 24);
  assert.throws(() => Wn.campaignMonteCarlo({ series: { hs: h }, dtHours: 24, operations: ops, startMonths: [12] }), RangeError);
});
