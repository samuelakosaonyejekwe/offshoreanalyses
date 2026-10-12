import test from 'node:test';
import assert from 'node:assert/strict';
import * as p from '../../src/engine/econ/production.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);
const sum = (a) => a.reduce((s, v) => s + v, 0);

test('Arps (1945) closed forms: exponential, harmonic and hyperbolic rate and cumulative at t = 5', () => {
  const a = { qi: 1000, Di: 0.2, t: 5 };
  close(p.arpsRate({ ...a, b: 0 }), 1000 * Math.exp(-1), 1e-14);
  close(p.arpsCumulative({ ...a, b: 0 }), (1000 - 1000 * Math.exp(-1)) / 0.2, 1e-13);       // 3160.6
  close(p.arpsRate({ ...a, b: 1 }), 500, 1e-14);
  close(p.arpsCumulative({ ...a, b: 1 }), 5000 * Math.LN2, 1e-13);                             // 3465.7
  close(p.arpsRate({ ...a, b: 0.5 }), 1000 / 2.25, 1e-13);
  close(p.arpsCumulative({ ...a, b: 0.5 }), 10000 * (1 - 1 / 1.5), 1e-12);                    // 3333.3
});

test('Arps cumulative equals the numerical integral of the rate for several exponents', () => {
  for (const b of [0, 0.3, 0.7, 1, 1.4]) {
    const n = 20000; const T = 8; const h = T / n; let s = 0;
    for (let k = 0; k <= n; k++) s += (k === 0 || k === n ? 0.5 : 1) * p.arpsRate({ qi: 750, Di: 0.35, b, t: k * h });
    close(p.arpsCumulative({ qi: 750, Di: 0.35, b, t: T }), s * h, 1e-7, `b=${b}`);
  }
});

test('time to economic limit and EUR: exponential (qi - qlim)/Di = 4500 at ln(10)/0.2 years', () => {
  const e = p.arpsEur({ qi: 1000, Di: 0.2, b: 0, qLimit: 100 });
  close(e.eur, 4500, 1e-13); close(e.timeToLimit, Math.log(10) / 0.2, 1e-14);
  const h = p.arpsEur({ qi: 1000, Di: 0.2, b: 0.5, qLimit: 100 });
  close(h.timeToLimit, (Math.sqrt(10) - 1) / 0.1, 1e-13);
  close(p.arpsRate({ qi: 1000, Di: 0.2, b: 0.5, t: h.timeToLimit }), 100, 1e-12);
  close(h.eur, 10000 * (1 - Math.sqrt(0.1)), 1e-12);
  assert.throws(() => p.arpsTimeToRate({ qi: 100, Di: 0.1, qLimit: 200 }), RangeError);
});

test('nominal and effective decline conversion', () => {
  close(p.declineConversion({ Di: 0.2 }).De, 1 - Math.exp(-0.2), 1e-15);
  close(p.declineConversion({ De: 0.15 }).Di, -Math.log(0.85), 1e-15);
  close(p.declineConversion({ De: p.declineConversion({ Di: 0.4, b: 0.6 }).De, b: 0.6 }).Di, 0.4, 1e-13);
});

test('decline-curve fitting recovers qi, Di and b from synthetic rate-time data', () => {
  const t = Array.from({ length: 24 }, (_, k) => k * 0.25);
  const hyp = p.arpsFit({ t, q: t.map((x) => p.arpsRate({ qi: 1200, Di: 0.45, b: 0.6, t: x })) });
  close(hyp.qi, 1200, 1e-6); close(hyp.Di, 0.45, 1e-5); close(hyp.b, 0.6, 1e-5); assert.ok(hyp.r2 > 0.999999);
  const ex = p.arpsFit({ t, q: t.map((x) => 800 * Math.exp(-0.3 * x)), model: 'exponential' });
  close(ex.qi, 800, 1e-8); close(ex.Di, 0.3, 1e-8); assert.equal(ex.b, 0);
  const har = p.arpsFit({ t, q: t.map((x) => 500 / (1 + 0.5 * x)), model: 'harmonic' });
  close(har.Di, 0.5, 1e-7);
});

test('plateau + decline profile: annual volumes integrate each segment exactly', () => {
  const r = p.plateauDeclineProfile({ plateauRate: 100, rampYears: 2, plateauYears: 3, Di: 0.2, years: 10 });
  close(r.annual[0], 25); close(r.annual[1], 75); close(r.annual[2], 100); close(r.annual[4], 100);
  close(r.annual[5], 100 * (1 - Math.exp(-0.2)) / 0.2, 1e-13);
  close(r.total, 100 + 300 + 100 * (1 - Math.exp(-1)) / 0.2, 1e-13);
  const lim = p.plateauDeclineProfile({ plateauRate: 100, plateauYears: 1, Di: 0.5, years: 12, qLimit: 20, uptime: 0.9 });
  close(lim.total, 0.9 * (100 + (100 - 20) / 0.5), 1e-12);            // production stops at the economic limit
  assert.equal(lim.lastEconomicYear, Math.ceil(1 + Math.log(5) / 0.5));
});

test('multi-well build-up: unconstrained profile equals the sum of type curves; capacity defers but does not lose volume', () => {
  const well = { qi: 100, Di: 0.3, b: 0 };
  const free = p.multiWellProfile({ wellsOnline: [2, 2, 0, 0, 0, 0], well, stepsPerYear: 400 });
  const yr = (k) => 100 * (Math.exp(-0.3 * k) - Math.exp(-0.3 * (k + 1))) / 0.3;       // one well's volume in its k-th year
  close(free.constrained[0], 2 * yr(0), 2e-3); close(free.constrained[1], 2 * yr(1) + 2 * yr(0), 2e-3);
  close(sum(free.deferred), 0, 1e-9);
  const cap = p.multiWellProfile({ wellsOnline: [2, 2, ...new Array(48).fill(0)], well, facilityCapacity: 250, stepsPerYear: 48 });
  for (const v of cap.constrained) assert.ok(v <= 250 + 1e-9);
  close(cap.constrained[1], 250, 1e-9); assert.ok(cap.deferred[1] > 0 && cap.deferred[2] < 0);
  // reserves are conserved: 4 wells x qi/Di
  close(cap.cumulative[49], 4 * 100 / 0.3, 1e-3);
  assert.deepEqual(cap.wellCount.slice(0, 3), [2, 4, 4]);
});

test('associated streams, gas shrinkage, throughput take-or-pay and ramp-up', () => {
  const s = p.associatedStreams({ oil: [100, 100], waterCut: [0, 0.5], gor: 800 });
  assert.deepEqual(s.water, [0, 100]); assert.deepEqual(s.gas, [80000, 80000]); assert.deepEqual(s.liquids, [100, 200]);
  const g = p.gasSales({ gross: [1000], fuelFraction: 0.05, flareFraction: 0.01, shrinkage: 0.1 });
  close(g.sales[0], 1000 * 0.94 * 0.9, 1e-12); close(g.sales[0] + g.fuel[0] + g.flare[0] + g.shrink[0], 1000, 1e-12);
  const t = p.throughputRevenue({ volumes: [50, 120, 80], tariff: 2, escalation: 0.1, capacity: 100, takeOrPayFraction: 0.8, contractedVolume: 100 });
  assert.deepEqual(t.billedVolume, [80, 100, 80]); close(t.revenue[1], 220, 1e-12); close(t.shortfallPayments[0], 60, 1e-12);
  const r = p.rampUpFactors({ years: 4, rampYears: 2, start: 0.5 });
  [0.625, 0.875, 1, 1].forEach((v, k) => close(r[k], v, 1e-12));
});

test('renewable generation: availability, losses, degradation and P90 = P50 (1 - 1.2816 sigma)', () => {
  const r = p.renewableGeneration({ capacityMW: 500, capacityFactor: 0.5, availability: 0.96, lossesFrac: 0.03, degradation: 0.005, years: 3, sigma: 0.08,
    monthlyShape: [3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 3, 3] });
  const y1 = 500 * 0.5 * 8760 * 0.96 * 0.97;
  close(r.p50[0], y1, 1e-12); close(r.p50[2], y1 * 0.995 ** 2, 1e-12);
  close(r.cases.P90[0], y1 * (1 - 1.2815515655446004 * 0.08), 1e-12);
  assert.deepEqual(r.cases.P50.map((v) => Math.round(v)), r.p50.map((v) => Math.round(v)));
  assert.ok(r.cases.P99[0] < r.cases.P90[0] && r.cases.P90[0] < r.cases.P75[0]);
  close(sum(r.monthlyP50Year1), y1, 1e-12); close(r.netCapacityFactor, 0.5 * 0.96 * 0.97, 1e-12);
  assert.match(r.convention, /exceeded with/);
  close(p.renewableGeneration({ grossAep_MWh: 1000, curtailment: [0, 0.1], years: 2 }).p50[1], 900);
});
