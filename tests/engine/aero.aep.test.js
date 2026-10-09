import test from 'node:test';
import assert from 'node:assert/strict';
import { gamma, gammaP } from '../../src/engine/math/special.js';
import {
  curvePower, cubicPowerCurve, powerCurveFromRotor, parseTurbineCsv, airDensityCorrection, airDensity, weibull,
  weibullScaleFromMean, extrapolateWind, aepWeibull, aepBins, normaliseRose, farmAep, lossCascade,
  normalQuantile, aepUncertainty, mcp, timeSeriesProduction, capacityValue, turbulenceSummary, shearSummary,
  extremeWind, generationProfile,
} from '../../src/engine/aero/aep.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg ?? ''} got ${a}, expected ${b}`);

test('Weibull mean A Γ(1 + 1/k); Rayleigh mean A √π/2', () => {
  close(weibull(10, 2).mean, (10 * Math.sqrt(Math.PI)) / 2, 1e-12);
  close(weibull(9, 2.3).mean, 9 * gamma(1 + 1 / 2.3), 1e-12);
  close(weibullScaleFromMean(8.86226925, 2), 10, 1e-8);
  // cdf at u = A is 1 − 1/e for any k
  close(weibull(7, 1.7).cdf(7), 1 - Math.exp(-1), 1e-12);
  // E[U³] for Rayleigh: A³ Γ(2.5) = 1.3293404 A³
  close(weibull(10, 2).cubeMean, 1329.3404, 1e-7);
});

test('AEP of an ideal cubic power curve under a Rayleigh climate: closed form', () => {
  const c = 0.5 * 1.225 * ((Math.PI * 100 * 100) / 4) * 0.45; // P = c U³
  const A = 9;
  // uncapped: mean power = c A³ Γ(5/2)
  const r = aepWeibull({ curve: (U) => c * U ** 3, A, k: 2, uMin: 0, uMax: 80, n: 8000 });
  close(r.meanPower, c * A ** 3 * gamma(2.5), 1e-8);
  close(r.aepWh, 8760 * c * A ** 3 * gamma(2.5), 1e-8);
});

test('AEP of the capped cubic curve: integration and bin method against the incomplete-gamma closed form', () => {
  const curve = cubicPowerCurve({ ratedPower: 3e6, D: 100, Cp: 0.45, cutIn: 3, cutOut: 25, step: 0.05 });
  const c = 0.5 * 1.225 * ((Math.PI * 100 * 100) / 4) * 0.45;
  const A = 9.5;
  const ur = Math.cbrt(3e6 / c);
  close(curve.ratedSpeed, ur, 1e-12);
  // ∫_a^b c u³ f du = c A³ Γ(2.5) [P(2.5, (b/A)²) − P(2.5, (a/A)²)] for k = 2
  const below = c * A ** 3 * gamma(2.5) * (gammaP(2.5, (ur / A) ** 2) - gammaP(2.5, (3 / A) ** 2));
  const w = weibull(A, 2);
  const exact = below + 3e6 * (w.cdf(25) - w.cdf(ur));
  close(aepWeibull({ curve, A, k: 2 }).meanPower, exact, 5e-4);
  close(aepBins({ curve, A, k: 2 }).meanPower, exact, 5e-4);
  // a coarse 1 m/s table with the bin method stays within 1 %
  const coarse = cubicPowerCurve({ ratedPower: 3e6, D: 100, Cp: 0.45, cutIn: 3, cutOut: 25, step: 1 });
  close(aepBins({ curve: coarse, A, k: 2 }).meanPower, exact, 1e-2);
});

test('generic power curve: Betz thrust at the Betz power coefficient, shedding above rated, smoothing', () => {
  const betz = cubicPowerCurve({ ratedPower: 5e6, D: 126, Cp: 16 / 27, efficiency: 1 });
  close(betz.ct[2], 8 / 9, 1e-9);
  const c = cubicPowerCurve({ ratedPower: 5e6, D: 126, Cp: 0.44, efficiency: 0.94 });
  const iLast = c.u.length - 1;
  assert.equal(c.p[iLast], 5e6);
  assert.ok(c.ct[iLast] < 0.15 && c.ct[2] > 0.5 && c.ct[2] < 8 / 9);
  // momentum check below rated: Cp/η = 4a(1 − a)², Ct = 4a(1 − a)
  const a = (1 - Math.sqrt(1 - c.ct[2])) / 2;
  close(4 * a * (1 - a) ** 2, 0.44 / 0.94, 1e-9);
  const s = cubicPowerCurve({ ratedPower: 5e6, D: 126, Cp: 0.44, smoothing: 10 });
  assert.ok(curvePower(s, c.ratedSpeed) < 5e6 && curvePower(s, c.ratedSpeed) > 0.9 * 5e6);
  assert.throws(() => cubicPowerCurve({ ratedPower: 5e6, D: 126, Cp: 0.7 }), RangeError);
});

test('power curve from a rotor callback and from an NREL-layout CSV', () => {
  const D = 100;
  const A = (Math.PI * D * D) / 4;
  const pc = powerCurveFromRotor({
    rotor: (U) => ({ power: 0.5 * 1.225 * A * 0.4 * U ** 3, thrust: 0.5 * 1.225 * A * 0.7 * U * U }),
    speeds: [4, 6, 8, 20],
    D,
    ratedPower: 2e6,
  });
  close(pc.cp[1], 0.4, 1e-12);
  close(pc.ct[1], 0.7, 1e-12);
  assert.equal(pc.p[3], 2e6);
  const csv = 'Wind Speed [m/s],Power [kW],Cp [-],Thrust [kN],Ct [-]\n3,40,0.20,170,0.80\n4,180,0.38,220,0.78\n5,400,0.44,340,0.77\n';
  const t = parseTurbineCsv(csv);
  assert.deepEqual(t.u, [3, 4, 5]);
  assert.deepEqual(t.p, [40e3, 180e3, 400e3]);
  assert.deepEqual(t.ct, [0.8, 0.78, 0.77]);
  assert.equal(t.thrust[2], 340e3);
  // Ct derived from thrust when the column is missing
  const t2 = parseTurbineCsv('Wind Speed [m/s],Power [MW],Thrust [kN]\n5,0.4,340\n10,3.2,900\n', { D: 100, rho: 1.225 });
  close(t2.ct[0], 340e3 / (0.5 * 1.225 * A * 25), 1e-12);
  close(t2.p[1], 3.2e6, 1e-12);
  assert.throws(() => parseTurbineCsv('a,b\n1,2\n'), RangeError);
});

test('air density and IEC 61400-12-1 density correction', () => {
  close(airDensity({ pressure: 101325, temperatureC: 15 }), 1.225, 1e-3);
  assert.ok(airDensity({ pressure: 101325, temperatureC: 15, relativeHumidity: 1 }) < 1.225);
  const curve = cubicPowerCurve({ ratedPower: 3e6, D: 100, Cp: 0.45, step: 0.25 });
  const site = airDensityCorrection(curve, { rho: 1.15 });
  // in the cubic region the power at a given speed scales with ρ/ρ0; rated power is unchanged
  close(curvePower(site, 8) / curvePower(curve, 8), 1.15 / 1.225, 2e-3);
  assert.equal(Math.max(...site.p), 3e6);
  const stall = airDensityCorrection(curve, { rho: 1.15, regulation: 'stall' });
  close(curvePower(stall, 20), (3e6 * 1.15) / 1.225, 1e-12);
});

test('vertical extrapolation by hand', () => {
  // power law: 8 (100/10)^0.14 = 8 · 1.380384 = 11.04307
  close(extrapolateWind({ U: 8, z1: 10, z2: 100, alpha: 0.14 }), 11.04307, 1e-6);
  // log law: 8 ln(100/0.0002)/ln(10/0.0002) = 8 · 13.122363/10.819778 = 9.70250
  close(extrapolateWind({ U: 8, z1: 10, z2: 100, z0: 0.0002 }), 9.7025, 1e-5);
});

test('loss cascade is multiplicative and reports capacity factor and full-load hours', () => {
  const r = lossCascade({
    grossWh: 1000e9,
    losses: { wake: 0.1, electrical: 0.03, availabilityTurbine: 0.04, curtailment: 0.01 },
    ratedPower: 250e6,
  });
  const net = 1000e9 * 0.9 * 0.97 * 0.96 * 0.99;
  close(r.netWh, net, 1e-12);
  close(r.totalLoss, 1 - 0.9 * 0.97 * 0.96 * 0.99, 1e-12);
  close(r.fullLoadHours, net / 250e6, 1e-12);
  close(r.capacityFactor, net / 250e6 / 8760, 1e-12);
  close(r.steps.reduce((s, x) => s + x.lostWh, 0), 1000e9 - net, 1e-9);
  assert.throws(() => lossCascade({ grossWh: 1, losses: { wake: 1.2 } }), RangeError);
});

test('uncertainty: root-sum-square combination and P90 = P50 (1 − 1.2816 σ)', () => {
  close(normalQuantile(0.9), 1.2815516, 1e-6);
  close(normalQuantile(0.75), 0.6744898, 1e-6);
  close(normalQuantile(0.99), 2.3263479, 1e-6);
  const r = aepUncertainty({
    p50Wh: 800e9,
    windSpeed: { measurement: 0.02, longTerm: 0.02, vertical: 0.01 },
    sensitivity: 2,
    energy: { wakeModel: 0.03, powerCurve: 0.04 },
    interAnnual: 0.06,
    years: 10,
  });
  // (2·0.02)² + (2·0.02)² + (2·0.01)² + 0.03² + 0.04² + 0.06²/10 = 0.0016+0.0016+0.0004+0.0009+0.0016+0.00036
  const sigma = Math.sqrt(0.00646);
  close(r.sigma, sigma, 1e-12);
  close(r.p90Wh, 800e9 * (1 - 1.2815516 * sigma), 1e-7);
  close(r.p75Wh, 800e9 * (1 - 0.6744898 * sigma), 1e-7);
  close(r.p99Wh, 800e9 * (1 - 2.3263479 * sigma), 1e-7);
  assert.ok(r.p99Wh < r.p90Wh && r.p90Wh < r.p75Wh && r.p75Wh < r.p50Wh);
  close(r.exceedance.p90, r.p90Wh, 1e-12);
  // a one-year horizon carries the full inter-annual variability
  assert.ok(aepUncertainty({ p50Wh: 1, interAnnual: 0.06, years: 1 }).sigma > aepUncertainty({ p50Wh: 1, interAnnual: 0.06, years: 20 }).sigma);
});

test('measure–correlate–predict: exact linear relation recovered; variance-ratio preserves the variance', () => {
  const ref = Array.from({ length: 200 }, (_, i) => 4 + 6 * (0.5 + 0.5 * Math.sin(i * 0.37)) + 0.01 * (i % 7));
  const site = ref.map((x) => 1.2 + 0.85 * x);
  const lt = ref.map((x) => x * 1.05);
  const r = mcp({ site, reference: ref, longTermReference: lt });
  close(r.slope, 0.85, 1e-10);
  close(r.intercept, 1.2, 1e-9);
  close(r.r2, 1, 1e-12);
  const ltMean = lt.reduce((s, x) => s + x, 0) / lt.length;
  close(r.longTermMean, 1.2 + 0.85 * ltMean, 1e-10);
  // scatter added: OLS shrinks the slope by r, variance ratio does not
  const noisy = site.map((y, i) => y + 0.8 * Math.sin(i * 12.9898) * Math.cos(i * 4.1));
  const ols = mcp({ site: noisy, reference: ref, longTermReference: ref });
  const vr = mcp({ site: noisy, reference: ref, longTermReference: ref, method: 'varianceRatio' });
  close(ols.slope, vr.slope * ols.r, 1e-10);
  const sd = (a) => {
    const m = a.reduce((s, x) => s + x, 0) / a.length;
    return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / a.length);
  };
  close(sd(Array.from(vr.predicted)), sd(noisy), 1e-9);
  close(vr.longTermMean, noisy.reduce((s, x) => s + x, 0) / noisy.length, 1e-10);
});

test('farm AEP: a single turbine reproduces the bin sum; an aligned pair loses energy only in the aligned sectors', () => {
  const curve = cubicPowerCurve({ ratedPower: 3e6, D: 100, Cp: 0.45, step: 0.5 });
  const rose = [0, 90, 180, 270].map((dirDeg) => ({ dirDeg, frequency: 0.25, A: 9, k: 2 }));
  const one = farmAep({ turbines: [{ x: 0, y: 0, D: 100, hubHeight: 90, curve }], rose });
  close(one.arrayEfficiency, 1, 1e-12);
  // independent aggregation: Σ over 1 m/s bins of P(u) [F(u + ½) − F(u − ½)]
  const w = weibull(9, 2);
  let s = 0;
  for (let u = 0.5; u <= 30; u += 1) s += curvePower(curve, u) * (w.cdf(u + 0.5) - w.cdf(Math.max(0, u - 0.5)));
  close(one.aepWh, s * 8760, 1e-10);
  close(one.capacityFactor, s / 3e6, 1e-10);
  const two = farmAep({
    turbines: [{ x: 0, y: 0, D: 100, hubHeight: 90, curve }, { x: 0, y: 600, D: 100, hubHeight: 90, curve }],
    rose,
    wake: { model: 'jensen' },
  });
  assert.ok(two.arrayEfficiency < 1 && two.arrayEfficiency > 0.8);
  close(two.grossWh, 2 * one.grossWh, 1e-10);
  close(two.perSector.reduce((a, x) => a + x.aepWh, 0), two.aepWh, 1e-10);
  close(two.perTurbine.reduce((a, x) => a + x.aepWh, 0), two.aepWh, 1e-10);
  // north and south sectors are waked, east and west are not
  assert.ok(two.perSector[0].efficiency < 0.95 && two.perSector[2].efficiency < 0.95);
  close(two.perSector[1].efficiency, 1, 1e-12);
  close(two.perSector[0].efficiency, two.perSector[2].efficiency, 1e-10);
  // northerly wind: the southern turbine (index 0) is the waked one; symmetric rose → equal annual totals
  close(two.perTurbine[0].aepWh, two.perTurbine[1].aepWh, 1e-10);
  close(two.wakeLoss, 1 - two.arrayEfficiency, 1e-15);
  assert.equal(normaliseRose([{ dirDeg: 0, frequency: 3, A: 8, k: 2 }])[0].frequency, 1);
});

test('time-series production: totals, monthly and diurnal aggregation', () => {
  const curve = { u: [0, 3, 12, 25], p: [0, 0, 2e6, 2e6] };
  const t0 = Date.UTC(2023, 0, 31, 0);
  const n = 48;
  const times = Array.from({ length: n }, (_, i) => t0 + i * 3600e3);
  const speeds = Array.from({ length: n }, (_, i) => (i % 24 < 12 ? 7.5 : 12));
  const r = timeSeriesProduction({ speeds, curve, times });
  // 7.5 m/s → 1 MW, 12 m/s → 2 MW: 24 h at each → 72 MWh
  close(r.energyWh, 72e6, 1e-12);
  close(r.capacityFactor, 0.75, 1e-12);
  assert.equal(r.monthly.length, 2);
  close(r.monthly[0].energyWh + r.monthly[1].energyWh, r.energyWh, 1e-12);
  assert.equal(r.daily.length, 2);
  close(r.diurnal[3], 1e6, 1e-12);
  close(r.diurnal[20], 2e6, 1e-12);
  const cv = capacityValue({ generation: r.power, demand: speeds.map((_, i) => (i % 24 >= 18 ? 10 : 5)), ratedPower: 2e6, topFraction: 0.25 });
  close(cv.capacityValue, 1, 1e-12);
  const g = generationProfile({ power: [1e6, 3e6], prices: [50, 100], exportLimit: 2e6 });
  close(g.totalMWh, 3, 1e-12);
  close(g.curtailedMWh, 1, 1e-12);
  close(g.revenue, 250, 1e-12);
  close(g.capturePrice, 250 / 3, 1e-12);
});

test('turbulence and shear summaries by hand', () => {
  const t = turbulenceSummary({ mean: [15.2, 14.9, 15.1, 8], std: [1.52, 1.788, 1.3288, 1.0] });
  // TI = 0.10, 0.12, 0.088 in the 15 m/s bin: mean 0.102667, sample std 0.016166, representative 0.123359
  const b = t.bins.find((x) => x.u === 15);
  close(b.meanTI, 0.1026667, 1e-6);
  close(b.representativeTI, 0.1026667 + 1.28 * 0.0161658, 1e-5);
  close(t.ti15, b.representativeTI, 1e-15);
  // exact power-law pair: α = 0.2
  const u1 = [6, 8, 10];
  const u2 = u1.map((u) => u * 2 ** 0.2);
  const s = shearSummary({ u1, u2, z1: 50, z2: 100 });
  close(s.alpha, 0.2, 1e-12);
  close(s.alphaMeanOfSamples, 0.2, 1e-12);
  // exact log-law pair recovers z0
  const z0 = 0.03;
  const v1 = [7];
  const v2 = [(7 * Math.log(100 / z0)) / Math.log(50 / z0)];
  close(shearSummary({ u1: v1, u2: v2, z1: 50, z2: 100 }).z0, z0, 1e-9);
});

test('50-year extreme wind from Gumbel-distributed annual maxima', () => {
  const loc = 28;
  const scale = 2.5;
  const n = 60;
  // exact Gumbel quantiles at Gringorten plotting positions
  const maxima = Array.from({ length: n }, (_, i) => loc - scale * Math.log(-Math.log((i + 1 - 0.44) / (n + 0.12))));
  const v50 = loc - scale * Math.log(-Math.log(1 - 1 / 50)); // 28 + 2.5·3.90194 = 37.755
  close(v50, 37.75485, 1e-6);
  for (const method of ['mom', 'lmom', 'mle']) {
    const r = extremeWind({ annualMaxima: maxima, method });
    close(r.value, v50, 0.03, method);
  }
  const g = extremeWind({ annualMaxima: maxima, gustFactor: 1.4, method: 'lmom' });
  close(g.value, 1.4 * extremeWind({ annualMaxima: maxima, method: 'lmom' }).value, 1e-12);
  assert.throws(() => extremeWind({ annualMaxima: [20, 21] }), RangeError);
});
