import test from 'node:test';
import assert from 'node:assert/strict';
import {
  extractReversals, rainflowThreePoint, rainflowFourPoint, closeResidue, stitchResidues, rainflowSegments,
  rangeHistogram, rangeMeanMatrix, fromToMatrix, damageEquivalentLoad, levelCrossings, cycleSummary,
} from '../../src/engine/fatigue/rainflow.js';
import { makeRng } from '../../src/engine/math/rng.js';

const ASTM = [-2, 1, -3, 5, -1, 3, -4, 4, -2];

function tally(cycles) {
  const m = new Map();
  for (const c of cycles) m.set(c.range, (m.get(c.range) || 0) + c.count);
  return Object.fromEntries([...m.entries()].sort((a, b) => a[0] - b[0]));
}

test('ASTM E1049-85 §5.4.4 worked example (Fig. 6): ranges 3, 4, 6, 8, 9 with counts 0.5, 1.5, 0.5, 1.0, 0.5', () => {
  const { cycles } = rainflowThreePoint(ASTM);
  assert.deepEqual(tally(cycles), { 3: 0.5, 4: 1.5, 6: 0.5, 8: 1.0, 9: 0.5 });
  // The single full cycle is E–F (−1 → 3), mean 1.
  const full = cycles.filter((c) => c.count === 1);
  assert.equal(full.length, 1);
  assert.equal(full[0].range, 4);
  assert.equal(full[0].mean, 1);
  assert.equal(cycleSummary(cycles).totalCycles, 4);
});

test('four-point counting with half-cycle residue equals the ASTM three-point count', () => {
  const r = rainflowFourPoint(ASTM, { residue: 'half' });
  assert.deepEqual(tally(r.cycles), { 3: 0.5, 4: 1.5, 6: 0.5, 8: 1.0, 9: 0.5 });
  assert.deepEqual(r.residue, [-2, 1, -3, 5, -4, 4, -2]);
});

test('repeat-residue method = counting the repeating history from its largest peak (hand count: 3, 4, 7, 9)', () => {
  // Hand count of the block rearranged to start at the maximum, 5, −1, 3, −4, 4, −2, 1, −3, 5:
  // −1→3 closes (4); −2→1 closes (3); 4→−3 closes (7); 5→−4→5 closes (9). Four full cycles.
  const r = rainflowFourPoint(ASTM, { residue: 'repeat' });
  assert.deepEqual(tally(r.cycles), { 3: 1, 4: 1, 7: 1, 9: 1 });
  const rotated = rainflowThreePoint([5, -1, 3, -4, 4, -2, 1, -3, 5]);
  const full = {};
  for (const c of rotated.cycles) full[c.range] = (full[c.range] || 0) + c.count;
  // The three-point rule sees the outermost loop as two half cycles (5 → −4 and −4 → 5), i.e. one cycle of 9.
  assert.deepEqual(full, { 3: 1, 4: 1, 7: 1, 9: 1 });
  assert.deepEqual(tally(closeResidue([-2, 1, -3, 5, -4, 4, -2], 'repeat')), { 3: 1, 7: 1, 9: 1 });
});

test('peak–valley extraction: plateaus and monotonic points removed, hysteresis gate removes small reversals', () => {
  assert.deepEqual(extractReversals([0, 1, 2, 2, 3, 2, 1, 1, 4]).values, [0, 3, 1, 4]);
  // The 0.2 ripple at the top is below a 0.5 gate.
  assert.deepEqual(extractReversals([0, 5, 4.8, 5.1, 0, 3], { gate: 0.5 }).values, [0, 5.1, 0, 3]);
  assert.deepEqual(extractReversals([0, 5, 4.8, 5.1, 0, 3]).values, [0, 5, 4.8, 5.1, 0, 3]);
  assert.deepEqual(extractReversals([0, 5, 4.8, 5.1, 0, 3]).indices, [0, 1, 2, 3, 4, 5]);
});

test('constant-amplitude signal: n cycles of range 2A', () => {
  const x = [];
  for (let i = 0; i <= 2000; i++) x.push(3 * Math.sin((2 * Math.PI * 10 * i) / 2000 + Math.PI / 2));
  const { cycles } = rainflowFourPoint(x, { residue: 'half' });
  const s = cycleSummary(cycles);
  assert.ok(Math.abs(s.totalCycles - 10) < 1e-9);
  assert.ok(Math.abs(s.maxRange - 6) < 1e-4);
});

test('segment stitching recovers exactly the cycles of the unsplit history', () => {
  const rng = makeRng(42);
  const x = [];
  let v = 0;
  for (let i = 0; i < 5000; i++) { v = 0.97 * v + rng.normal(); x.push(v + 5 * Math.sin(i / 400)); }
  const whole = rainflowFourPoint(x, { residue: 'repeat' }).cycles;
  const parts = rainflowSegments([x.slice(0, 1234), x.slice(1234, 3999), x.slice(3999)], { residue: 'repeat' }).cycles;
  const key = (cs) => cs.map((c) => c.range).sort((a, b) => a - b);
  const a = key(whole);
  const b = key(parts);
  assert.equal(a.length, b.length);
  for (let i = 0; i < a.length; i++) assert.ok(Math.abs(a[i] - b[i]) < 1e-12);
  // Without stitching the largest cycle is lost: per-segment residues counted separately give a smaller maximum.
  const st = stitchResidues([[0, 5, 1], [1, 4, -6, 2]], { residue: 'half' });
  assert.equal(cycleSummary(st.cycles).maxRange, 11);
});

test('damage-equivalent load: hand calculation, m = 3', () => {
  // (10·2³ + 1·10³)/100 = 10.8 → DEL = 10.8^(1/3) = 2.2104
  const del = damageEquivalentLoad({ ranges: [2, 10], counts: [10, 1], m: 3, nRef: 100 });
  assert.ok(Math.abs(del - 10.8 ** (1 / 3)) < 1e-12);
  // A constant-amplitude history with nRef equal to its own count returns its range.
  const cyc = [{ range: 7, count: 50 }];
  assert.ok(Math.abs(damageEquivalentLoad({ cycles: cyc, m: 4, nRef: 50 }) - 7) < 1e-12);
});

test('histogram, range–mean and from–to matrices conserve the cycle count', () => {
  const { cycles } = rainflowThreePoint(ASTM);
  const total = cycleSummary(cycles).totalCycles;
  const h = rangeHistogram(cycles, { bins: 9, maxRange: 9 });
  assert.equal(h.counts.reduce((a, b) => a + b, 0), total);
  assert.equal(h.counts[3], 0.5); // range 3 falls in [3, 4)
  assert.equal(h.counts[4], 1.5);
  assert.equal(h.counts[8], 1.5); // ranges 8 (1.0) and 9 (0.5, top edge) share the last bin
  const rm = rangeMeanMatrix(cycles, { rangeBins: 4, meanBins: 3 });
  assert.equal(rm.matrix.flat().reduce((a, b) => a + b, 0), total);
  const ft = fromToMatrix(cycles, { bins: 9, min: -4, max: 5 });
  assert.equal(ft.matrix.flat().reduce((a, b) => a + b, 0), total);
  // Cycle −1 → 3: from class floor((−1+4)/1) = 3 to class 7.
  assert.equal(ft.matrix[3][7], 1);
});

test('level crossings of the ASTM sequence (positive slope)', () => {
  const r = levelCrossings(ASTM, [0, 2, 4.5, 6]);
  assert.deepEqual(Array.from(r.upCrossings), [4, 3, 1, 0]);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => extractReversals([1, 2], { gate: -1 }), /gate/);
  assert.throws(() => damageEquivalentLoad({ ranges: [1], counts: [1], m: 0, nRef: 1 }), /m must/);
  assert.throws(() => closeResidue([1, 2, 1], 'bogus'), RangeError);
});
