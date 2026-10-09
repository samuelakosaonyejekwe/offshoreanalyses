import test from 'node:test';
import assert from 'node:assert/strict';
import {
  ageReplacementCostRate, optimalAgeReplacement, blockReplacementCostRate, optimalBlockReplacement, costRateCurves,
  kijimaExpectedFailures, cbmGammaProcess, optimiseCbmThreshold, delayTimeModel, optimalInspectionInterval,
  expectedDowntime, genericPowerFraction, simulateWindOM,
} from '../../src/engine/ram/maintenance.js';
import { weibullRenewals } from '../../src/engine/ram/lifedata.js';
import { gamma } from '../../src/engine/math/special.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

// Independent evaluation of the age-replacement cost rate with a fine trapezoidal integral of R(t).
function bruteCostRate(T, shape, scale, cp, cf) {
  const n = 20000;
  const h = T / n;
  let integral = 0;
  for (let i = 0; i <= n; i++) integral += (i === 0 || i === n ? 0.5 : 1) * Math.exp(-(((i * h) / scale) ** shape)) * h;
  const R = Math.exp(-((T / scale) ** shape));
  return (cp * R + cf * (1 - R)) / integral;
}

test('age replacement: cost rate against direct numerical integration; exponential life gives a constant-hazard result', () => {
  const r = ageReplacementCostRate({ interval: 600, shape: 2.5, scale: 1000, costPreventive: 1, costFailure: 10 });
  close(r.costRate, bruteCostRate(600, 2.5, 1000, 1, 10), 1e-7);
  close(r.preventiveShare + r.correctiveShare, r.costRate, 1e-14);
  // β = 1: ∫R = η(1 − e^(−T/η)) → C(T) = [cp·e^(−x) + cf·(1 − e^(−x))]/(η(1 − e^(−x)))
  const e = ageReplacementCostRate({ interval: 500, shape: 1, scale: 1000, costPreventive: 1, costFailure: 10 });
  const q = 1 - Math.exp(-0.5);
  close(e.costRate, (Math.exp(-0.5) + 10 * q) / (1000 * q), 1e-10);
});

test('optimal age-replacement interval against a numerically minimised cost rate and the optimality condition', () => {
  const shape = 2.5; const scale = 1000; const cp = 1; const cf = 10;
  const opt = optimalAgeReplacement({ shape, scale, costPreventive: cp, costFailure: cf });
  // Brute-force minimum on a 1-hour grid.
  let best = Infinity;
  let tBest = 0;
  for (let T = 100; T <= 1500; T += 1) {
    const c = bruteCostRate(T, shape, scale, cp, cf);
    if (c < best) { best = c; tBest = T; }
  }
  assert.ok(Math.abs(opt.interval - tBest) <= 1, `${opt.interval} vs ${tBest}`);
  close(opt.costRate, best, 1e-5);
  assert.ok(opt.costRate <= best);
  // First-order condition (Barlow & Proschan): h(T*)·∫₀^T* R dt − F(T*) = cp/(cf − cp)
  const T = opt.interval;
  const h = (shape / scale) * (T / scale) ** (shape - 1);
  const cyc = ageReplacementCostRate({ interval: T, shape, scale, costPreventive: cp, costFailure: cf }).meanCycle;
  close(h * cyc - (1 - Math.exp(-((T / scale) ** shape))), cp / (cf - cp), 1e-5);
  // ... and at the optimum the cost rate equals (cf − cp)·h(T*).
  close(opt.costRate, (cf - cp) * h, 1e-5);
  close(opt.runToFailureCostRate, cf / (scale * gamma(1 + 1 / shape)), 1e-12);
  assert.ok(opt.saving > 0.3);
  // No finite optimum for a constant or decreasing hazard, or when prevention costs as much as failure.
  assert.equal(optimalAgeReplacement({ shape: 1, scale, costPreventive: cp, costFailure: cf }).interval, Infinity);
  assert.equal(optimalAgeReplacement({ shape: 0.7, scale, costPreventive: cp, costFailure: cf }).interval, Infinity);
  assert.equal(optimalAgeReplacement({ shape: 3, scale, costPreventive: 5, costFailure: 5 }).interval, Infinity);
  const curves = costRateCurves({ intervals: [200, 400, T, 900, 1400], shape, scale, costPreventive: cp, costFailure: cf });
  assert.ok(curves.total[2] <= Math.min(curves.total[0], curves.total[1], curves.total[3], curves.total[4]));
  assert.ok(curves.preventive[0] > curves.preventive[4] && curves.corrective[0] < curves.corrective[4]);
});

test('block replacement: cost rate uses the renewal function; optimum satisfies T·m(T) − M(T) = cp/cf', () => {
  const shape = 3; const scale = 1000; const cp = 1; const cf = 8;
  const c = blockReplacementCostRate({ interval: 500, shape, scale, costPreventive: cp, costFailure: cf });
  close(c.expectedFailures, weibullRenewals({ shape, scale, t: 500, steps: 600 }), 1e-12);
  close(c.costRate, (cp + cf * c.expectedFailures) / 500, 1e-14);
  const opt = optimalBlockReplacement({ shape, scale, costPreventive: cp, costFailure: cf });
  assert.ok(Number.isFinite(opt.interval));
  // Optimality: d/dT[(cp + cf·M)/T] = 0 ⇔ T·M'(T) − M(T) = cp/cf, with M' by central difference.
  const d = 5;
  const M = (t) => weibullRenewals({ shape, scale, t, steps: 1500 });
  const lhs = (opt.interval * (M(opt.interval + d) - M(opt.interval - d))) / (2 * d) - M(opt.interval);
  close(lhs, cp / cf, 0.03);
  // Age replacement is never worse than block replacement.
  assert.ok(optimalAgeReplacement({ shape, scale, costPreventive: cp, costFailure: cf }).costRate <= opt.costRate * (1 + 1e-9));
  // Exponential life: M(T) = T/η, cost rate strictly decreasing → no finite optimum.
  assert.equal(optimalBlockReplacement({ shape: 1, scale, costPreventive: cp, costFailure: cf }).interval, Infinity);
});

test('Kijima type I imperfect repair: limits q = 1 (minimal repair, E[N] = (t/η)^β) and q = 0 (renewal process)', () => {
  const bad = kijimaExpectedFailures({ shape: 2, scale: 100, q: 1, horizon: 300, n: 4000, seed: 2 });
  assert.ok(Math.abs(bad.expectedFailures - 9) < 4 * bad.stdError + 1e-9, `${bad.expectedFailures}`);
  const good = kijimaExpectedFailures({ shape: 2, scale: 100, q: 0, horizon: 300, n: 4000, seed: 2 });
  const M = weibullRenewals({ shape: 2, scale: 100, t: 300, steps: 1500 });
  assert.ok(Math.abs(good.expectedFailures - M) < 4 * good.stdError + 1e-3, `${good.expectedFailures} vs ${M}`);
  const mid = kijimaExpectedFailures({ shape: 2, scale: 100, q: 0.5, horizon: 300, n: 4000, seed: 2 });
  assert.ok(mid.expectedFailures > good.expectedFailures && mid.expectedFailures < bad.expectedFailures);
});

test('condition-based maintenance on a gamma process: mean cycle against the renewal estimate; interior optimum of the threshold', () => {
  const base = { shapeRate: 2, scale: 0.5, failureLevel: 20, inspectionInterval: 1, costInspection: 1, costPreventive: 20, costFailure: 200, n: 6000, seed: 4 };
  // Threshold at the failure level: every cycle ends in failure; mean degradation rate 1/yr → first passage of
  // level 20 at about 20 + (mean overshoot ≈ (σ²/μ + μ)/2 = 0.75) years, found at the next inspection.
  const rtf = cbmGammaProcess({ ...base, preventiveThreshold: 20 });
  assert.equal(rtf.failureFraction, 1);
  close(rtf.meanCycleLength, 20.75, 0.03);
  const opt = optimiseCbmThreshold(base);
  assert.ok(opt.preventiveThreshold > 5 && opt.preventiveThreshold < 20, `threshold ${opt.preventiveThreshold}`);
  assert.ok(opt.costRate < rtf.costRate);
  assert.ok(opt.costRate < cbmGammaProcess({ ...base, preventiveThreshold: 2 }).costRate);
  assert.equal(opt.grid.length, 20);
  // A lower threshold reduces the share of cycles ending in failure.
  assert.ok(cbmGammaProcess({ ...base, preventiveThreshold: 10 }).failureFraction < cbmGammaProcess({ ...base, preventiveThreshold: 18 }).failureFraction);
  // Deterministic for a seed.
  assert.equal(cbmGammaProcess({ ...base, preventiveThreshold: 15 }).costRate, cbmGammaProcess({ ...base, preventiveThreshold: 15 }).costRate);
});

test('delay-time model: exponential delay time gives b(T) = 1 − (1 − e^(−αT))/(αT); optimum against a grid search', () => {
  const alpha = 0.05;
  const args = { defectRate: 0.2, delayCdf: (h) => 1 - Math.exp(-alpha * h), costInspection: 50, costFailure: 1000, costDefectRepair: 100 };
  const r = delayTimeModel({ ...args, interval: 30 });
  const b = 1 - (1 - Math.exp(-1.5)) / 1.5;
  close(r.failureFraction, b, 1e-8);
  close(r.costRate, (6 * (b * 1000 + (1 - b) * 100) + 50) / 30, 1e-8);
  close(r.failuresPerInterval, 6 * b, 1e-8);
  const opt = optimalInspectionInterval({ ...args, minInterval: 0.5, maxInterval: 200 });
  let best = Infinity;
  let tBest = 0;
  for (let T = 0.5; T <= 200; T += 0.01) {
    const bb = 1 - (1 - Math.exp(-alpha * T)) / (alpha * T);
    const c = (0.2 * T * (bb * 1000 + (1 - bb) * 100) + 50) / T;
    if (c < best) { best = c; tBest = T; }
  }
  close(opt.interval, tBest, 5e-3);
  close(opt.costRate, best, 1e-6);
  assert.equal(opt.atBound, false);
  // Downtime objective with inspection downtime.
  const d = optimalInspectionInterval({ ...args, downtimeInspection: 0.2, downtimeFailure: 3, minInterval: 0.5, maxInterval: 200, objective: 'downtime' });
  const dr = (T) => delayTimeModel({ ...args, downtimeInspection: 0.2, downtimeFailure: 3, interval: T }).downtimeRate;
  assert.ok(dr(d.interval) <= dr(d.interval * 1.1) && dr(d.interval) <= dr(d.interval * 0.9));
});

test('expected downtime and downtime cost: hand calculation', () => {
  // λ = 1e-3 /h, MTTR 20 h, delay 80 h → A = 1/(1 + 0.1) = 0.90909; failures/yr = 1e-3·8760·0.90909 = 7.9636
  const r = expectedDowntime({ failureRate: 1e-3, mttr: 20, logisticDelay: 80, downtimeCostPerHour: 500, repairCost: 10000 });
  close(r.availability, 1 / 1.1, 1e-12);
  close(r.failuresPerYear, 8.76 / 1.1, 1e-12);
  close(r.downtimeHoursPerYear, 876 / 1.1, 1e-12);
  close(r.totalCostPerYear, (876 / 1.1) * 500 + (8.76 / 1.1) * 10000, 1e-12);
  close(r.downtimeHoursPerYear, 8760 * (1 - r.availability), 1e-12);
});

test('offshore-wind O&M simulator: unrestricted access reproduces A = 1/(1 + λ·MTTR); weather and vessel limits reduce it', () => {
  // Constant calm weather at rated wind.
  const n = 8760;
  const calm = { hs: new Float64Array(n).fill(1), wind: new Float64Array(n).fill(12), dtHours: 1 };
  const components = [
    { name: 'minor', failureRate: 6, repairHours: 8, vessel: 'ctv', partsCost: 2000, technicians: 3 },
    { name: 'major', failureRate: 1, repairHours: 60, vessel: 'jackup', partsCost: 150000, leadTimeHours: 100, technicians: 6 },
  ];
  const vessels = { ctv: { hsLimit: 1.5, dayRate: 3000, count: 3 }, jackup: { hsLimit: 2.0, windLimit: 15, dayRate: 150000, charter: 'on-demand', mobilisationHours: 300, mobilisationCost: 500000 } };
  const base = { turbines: 40, ratedPower: 8e6, years: 5, weather: calm, components, vessels, energyPricePerMWh: 80, labourRatePerHour: 100, simulations: 6, seed: 11 };
  const r = simulateWindOM(base);
  // Downtime per failure: minor 8 h, major 300 h (mobilisation governs the lead time) + 60 h.
  // Failures per operating year: 7 → mean downtime per operating hour = (6·8 + 1·360)/8760 = 0.046575 → A = 0.95550
  const expected = 1 / (1 + (6 * 8 + 360) / 8760);
  close(r.timeAvailability, expected, 0.004);
  close(r.energyAvailability, r.timeAvailability, 1e-12); // constant wind: energy- and time-based availability coincide
  close(r.failuresPerYear, 40 * 7 * expected, 0.05);
  close(r.byComponent[1].meanDowntimeHours, 360, 0.02);
  close(r.capacityFactor, r.energyAvailability, 1e-9);
  // Costs: parts ≈ failures × part cost; lost revenue = lost energy × price; CTVs paid every day.
  close(r.costPerYear.parts, r.byComponent[0].failuresPerYear * 2000 + r.byComponent[1].failuresPerYear * 150000, 1e-9);
  close(r.costPerYear.lostRevenue, r.energyLostMWhPerYear * 80, 1e-9);
  close(r.energyLostMWhPerYear, 40 * 8 * 8760 * (1 - r.energyAvailability), 1e-6);
  close(r.costPerYear.labour, (r.byComponent[0].failuresPerYear * 3 * 8 + r.byComponent[1].failuresPerYear * 6 * 60) * 100, 1e-9);
  const jackupJobs = r.byComponent[1].failuresPerYear;
  close(r.costPerYear.vessels, 3 * 3000 * 365 + jackupJobs * 500000 + (r.vesselUtilisation.jackup * 8760 / 24) * 150000, 0.02);
  close(r.costPerYear.total, r.costPerYear.vessels + r.costPerYear.parts + r.costPerYear.labour + r.costPerYear.lostRevenue, 1e-12);
  // CTV utilisation = busy hours / (3 vessels × hours): 6·8 h per turbine-year × 40 / (3·8760)
  close(r.vesselUtilisation.ctv, (r.byComponent[0].failuresPerYear * 8) / (3 * 8760), 0.02);
  assert.equal(r.accessibility.ctv, 1);
  // Same seed → identical result.
  assert.equal(simulateWindOM(base).timeAvailability, r.timeAvailability);

  // Seasonal weather: half of each day above the CTV limit → access 50 %, availability drops.
  const rough = { hs: Float64Array.from({ length: n }, (_, i) => (i % 24 < 12 ? 1 : 2.5)), wind: calm.wind, dtHours: 1 };
  const w = simulateWindOM({ ...base, weather: rough });
  assert.equal(w.accessibility.ctv, 0.5);
  assert.ok(w.timeAvailability < r.timeAvailability - 0.002);
  // A single CTV with long repairs queues the work.
  const q = simulateWindOM({ ...base, components: [{ ...components[0], repairHours: 30 }, components[1]], vessels: { ...vessels, ctv: { hsLimit: 1.5, dayRate: 3000, count: 1 } } });
  const u = simulateWindOM({ ...base, components: [{ ...components[0], repairHours: 30 }, components[1]] });
  assert.ok(q.timeAvailability < u.timeAvailability);
  assert.ok(q.vesselUtilisation.ctv > 0.7 && q.vesselUtilisation.ctv <= 1);
});

test('energy-based availability differs from time-based when downtime correlates with wind; generic power curve', () => {
  assert.equal(genericPowerFraction(2), 0);
  assert.equal(genericPowerFraction(12), 1);
  assert.equal(genericPowerFraction(26), 0);
  close(genericPowerFraction(8), (512 - 27) / (1728 - 27), 1e-12);
  // Windy days are also the inaccessible ones, so short repairs wait through high-production hours: a failure on a
  // rough day loses on average 12 h at full power, one on a calm day about 2 h at 6 % power.
  const n = 8760;
  const hs = Float64Array.from({ length: n }, (_, i) => (i % 48 < 24 ? 0.8 : 2.5));
  const wind = Float64Array.from({ length: n }, (_, i) => (i % 48 < 24 ? 5 : 14));
  const r = simulateWindOM({
    turbines: 30, ratedPower: 5e6, years: 4, weather: { hs, wind, dtHours: 1 },
    components: [{ name: 'c', failureRate: 20, repairHours: 2, vessel: 'ctv' }], vessels: { ctv: { hsLimit: 1.5 } }, simulations: 4, seed: 3,
  });
  assert.ok(r.energyAvailability < r.timeAvailability - 0.003, `${r.energyAvailability} < ${r.timeAvailability}`);
  // Time-based: mean downtime ≈ 0.5·(12 + 2) + 0.5·2 = 8 h per failure → A ≈ 1/(1 + 20·8/8760) = 0.982
  close(r.timeAvailability, 1 / (1 + 160 / 8760), 0.004);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => ageReplacementCostRate({ interval: 0, shape: 2, scale: 1, costPreventive: 1, costFailure: 2 }), /interval/);
  assert.throws(() => cbmGammaProcess({ shapeRate: 1, scale: 1, failureLevel: 10, preventiveThreshold: 12, inspectionInterval: 1, costInspection: 1, costPreventive: 1, costFailure: 1 }), /preventiveThreshold/);
  assert.throws(() => simulateWindOM({ turbines: 1, ratedPower: 1, years: 1, weather: { hs: [1], wind: [1], dtHours: 1 }, components: [{ name: 'x', failureRate: 1, repairHours: 1, vessel: 'sov' }], vessels: {} }), /vessel/);
});
