import test from 'node:test';
import assert from 'node:assert/strict';
import {
  lifeReliability, minimalPaths, systemReliability, reliabilityCurve, mttf, availabilityOf, seriesRepairable,
  parallelRepairable,
} from '../../src/engine/ram/rbd.js';
import { gamma } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);
const comp = (R) => ({ type: 'component', life: { dist: 'constant', reliability: R } });
const expo = (rate) => ({ type: 'component', life: { dist: 'exponential', rate } });

test('component life distributions', () => {
  close(lifeReliability({ dist: 'exponential', rate: 0.01 }, 100), Math.exp(-1), 1e-14);
  close(lifeReliability({ dist: 'weibull', shape: 2, scale: 100 }, 50), Math.exp(-0.25), 1e-14);
  close(lifeReliability({ dist: 'lognormal', mu: Math.log(100), sigma: 0.5 }, 100), 0.5, 1e-12);
  assert.equal(lifeReliability({ R: (t) => 1 / (1 + t) }, 3), 0.25);
});

test('series, parallel and 2-out-of-3: R = 3R² − 2R³', () => {
  close(systemReliability({ type: 'series', children: [comp(0.9), comp(0.8), comp(0.95)] }, 1), 0.684, 1e-12);
  close(systemReliability({ type: 'parallel', children: [comp(0.9), comp(0.8)] }, 1), 0.98, 1e-12);
  const R = 0.9;
  close(systemReliability({ type: 'koon', k: 2, children: [comp(R), comp(R), comp(R)] }, 1), 3 * R * R - 2 * R ** 3, 1e-12);
  close(3 * R * R - 2 * R ** 3, 0.972, 1e-12);
  // Non-identical 2oo3: R1R2 + R1R3 + R2R3 − 2R1R2R3 = 0.72 + 0.63 + 0.56 − 1.008 = 0.902
  close(systemReliability({ type: 'koon', k: 2, children: [comp(0.9), comp(0.8), comp(0.7)] }, 1), 0.902, 1e-12);
  // Time-dependent 2oo3 with exponential units: 3e^(−2λt) − 2e^(−3λt)
  const sys = { type: 'koon', k: 2, children: [expo(1e-3), expo(1e-3), expo(1e-3)] };
  close(systemReliability(sys, 500), 3 * Math.exp(-1) - 2 * Math.exp(-1.5), 1e-12);
  // Nested: series of a parallel pair and a single block.
  close(systemReliability({ type: 'series', children: [{ type: 'parallel', children: [comp(0.9), comp(0.9)] }, comp(0.95)] }, 1), 0.99 * 0.95, 1e-12);
});

test('bridge network by path enumeration: R = 2R² + 2R³ − 5R⁴ + 2R⁵', () => {
  const R = 0.9;
  const bridge = {
    type: 'network', source: 's', sink: 't', edges: [
      { from: 's', to: 'a', block: comp(R) }, { from: 's', to: 'b', block: comp(R) },
      { from: 'a', to: 't', block: comp(R) }, { from: 'b', to: 't', block: comp(R) },
      { from: 'a', to: 'b', block: comp(R) },
    ],
  };
  assert.equal(minimalPaths(bridge).length, 4);
  close(systemReliability(bridge, 1), 2 * R ** 2 + 2 * R ** 3 - 5 * R ** 4 + 2 * R ** 5, 1e-12);
  close(2 * R ** 2 + 2 * R ** 3 - 5 * R ** 4 + 2 * R ** 5, 0.97848, 1e-12);
  // Unequal blocks against conditioning on the bridge element (factoring theorem).
  const r = [0.9, 0.8, 0.7, 0.6, 0.5];
  const net = { ...bridge, edges: bridge.edges.map((e, i) => ({ ...e, block: comp(r[i]) })) };
  const up = (1 - (1 - r[0]) * (1 - r[1])) * (1 - (1 - r[2]) * (1 - r[3]));
  const down = 1 - (1 - r[0] * r[2]) * (1 - r[1] * r[3]);
  close(systemReliability(net, 1), r[4] * up + (1 - r[4]) * down, 1e-12);
});

test('cold standby: identical exponential units R = e^(−λt)(1 + λt); imperfect switch; three units', () => {
  const lam = 1e-3;
  const t = 1500;
  const two = { type: 'standby', children: [expo(lam), expo(lam)] };
  close(systemReliability(two, t), Math.exp(-lam * t) * (1 + lam * t), 1e-5);
  const imperfect = { type: 'standby', switchReliability: 0.9, children: [expo(lam), expo(lam)] };
  close(systemReliability(imperfect, t), Math.exp(-lam * t) * (1 + 0.9 * lam * t), 1e-5);
  const three = { type: 'standby', children: [expo(lam), expo(lam), expo(lam)] };
  close(systemReliability(three, t), Math.exp(-lam * t) * (1 + lam * t + (lam * t) ** 2 / 2), 1e-5);
  // Different rates: R = e^(−λ1 t) + λ1/(λ1 − λ2)·(e^(−λ2 t) − e^(−λ1 t))
  const l2 = 3e-3;
  const mixed = { type: 'standby', children: [expo(lam), expo(l2)] };
  close(systemReliability(mixed, t), Math.exp(-lam * t) + (lam / (lam - l2)) * (Math.exp(-l2 * t) - Math.exp(-lam * t)), 1e-4);
  // Standby beats active redundancy.
  assert.ok(systemReliability(two, t) > systemReliability({ type: 'parallel', children: [expo(lam), expo(lam)] }, t));
});

test('MTTF by integration: exponential structures and Weibull against closed forms', () => {
  const lam = 2e-3;
  close(mttf(expo(lam), { scale: 100 }), 1 / lam, 1e-6);
  close(mttf({ type: 'parallel', children: [expo(lam), expo(lam)] }, { scale: 100 }), 1.5 / lam, 1e-6);
  close(mttf({ type: 'koon', k: 2, children: [expo(lam), expo(lam), expo(lam)] }, { scale: 100 }), 5 / (6 * lam), 1e-6);
  close(mttf({ type: 'series', children: [expo(lam), expo(3 * lam)] }, { scale: 100 }), 1 / (4 * lam), 1e-6);
  for (const shape of [0.8, 1.5, 3]) {
    close(mttf({ type: 'component', life: { dist: 'weibull', shape, scale: 1000 } }), 1000 * gamma(1 + 1 / shape), 1e-4, `shape ${shape}`);
  }
  close(mttf({ type: 'standby', children: [expo(lam), expo(lam)] }, { scale: 100 }), 2 / lam, 2e-3);
  const c = reliabilityCurve(expo(lam), [0, 100, 500]);
  close(c[2], Math.exp(-1), 1e-12);
  assert.equal(c[0], 1);
});

test('availability of repairable series and parallel systems', () => {
  // A = μ/(λ + μ): 0.1/(0.001 + 0.1) = 0.990099
  const a = { type: 'component', lambda: 1e-3, mu: 0.1 };
  const b = { type: 'component', mtbf: 2000, mttr: 20 };
  close(availabilityOf(a), 0.1 / 0.101, 1e-14);
  close(availabilityOf({ type: 'series', children: [a, b] }), (0.1 / 0.101) * (2000 / 2020), 1e-14);
  close(availabilityOf({ type: 'parallel', children: [a, a] }), 1 - (0.001 / 0.101) ** 2, 1e-14);
  const A = 0.1 / 0.101;
  close(availabilityOf({ type: 'koon', k: 2, children: [a, a, a] }), 3 * A * A - 2 * A ** 3, 1e-12);
  // Frequency–duration: series MUT = 1/Σλ; parallel MDT = 1/Σμ; A = MUT/(MUT + MDT) in both.
  const s = seriesRepairable([{ lambda: 1e-3, mu: 0.1 }, { lambda: 2e-3, mu: 0.05 }]);
  close(s.availability, (0.1 / 0.101) * (0.05 / 0.052), 1e-14);
  close(s.meanUpTime, 1 / 3e-3, 1e-12);
  close(s.meanUpTime / (s.meanUpTime + s.meanDownTime), s.availability, 1e-12);
  const p = parallelRepairable([{ lambda: 1e-3, mu: 0.1 }, { lambda: 2e-3, mu: 0.05 }]);
  close(p.meanDownTime, 1 / 0.15, 1e-12);
  close(p.meanUpTime / (p.meanUpTime + p.meanDownTime), p.availability, 1e-12);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => systemReliability({ type: 'koon', k: 4, children: [comp(0.9)] }, 1), /koon/);
  assert.throws(() => systemReliability({ type: 'ring', children: [] }, 1), /unknown block type/);
  assert.throws(() => lifeReliability({ dist: 'gamma' }, 1), /unknown distribution/);
  assert.throws(() => systemReliability(expo(1), -1), /t must/);
});
