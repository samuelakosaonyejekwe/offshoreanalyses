import test from 'node:test';
import assert from 'node:assert/strict';
import {
  expIntegralE1, vanWormer, cableLadder, ladderStepResponse, twoLoopResponse, soilStepResponse, transientModel,
  responseTime, temperatureHistory, exceedanceStatistics, cyclicRatingFactor, cyclicRatingFactorIEC, emergencyRating,
  adiabaticFinalTemperature, solveTrench2D,
} from '../../src/engine/electrical/thermaltransient.js';
import { cableGeometry, t4Buried, mutualHeatingRise } from '../../src/engine/electrical/ampacity.js';

const rel = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${a} vs ${b} (tol ${tol})`);

const single = {
  cores: 1,
  conductor: { material: 'Cu', area: 1000e-6, diameter: 0.0376 },
  conductorScreenThickness: 0.0012,
  insulation: { material: 'xlpe', thickness: 0.016 },
  insulationScreenThickness: 0.0012,
  sheath: { material: 'lead', thickness: 0.003 },
  serving: { thickness: 0.004, rhoT: 3.5 },
};
const g = cableGeometry(single);
const buried = { type: 'buried', rhoSoil: 1.0, diffusivity: 0.5e-6, depth: 1.0 };

test('exponential integral E1 against tabulated values (Abramowitz & Stegun Table 5.1)', () => {
  rel(expIntegralE1(0.1), 1.8229239585, 1e-9);
  rel(expIntegralE1(1), 0.2193839344, 1e-9);
  rel(expIntegralE1(2), 0.0489005107, 1e-9);
  rel(expIntegralE1(5), 0.0011482956, 1e-8);
  // small-argument limit E1(x) → −γ − ln x
  rel(expIntegralE1(1e-8), -0.5772156649 - Math.log(1e-8), 1e-8);
  assert.throws(() => expIntegralE1(0), RangeError);
});

test('Van Wormer coefficient by hand for Di/dc = 2: p = 1/(2 ln 2) − 1/3 = 0.38801; short duration p* = 1/ln 2 − 1 = 0.44270', () => {
  rel(vanWormer({ outerDiameter: 0.08, innerDiameter: 0.04 }), 0.388014, 1e-5);
  rel(vanWormer({ outerDiameter: 0.08, innerDiameter: 0.04, shortDuration: true }), 0.442695, 1e-5);
  // thin layer limit: half the capacitance at each side
  rel(vanWormer({ outerDiameter: 1.0001, innerDiameter: 1 }), 0.5, 1e-4);
});

test('single thermal capacitance: θ(t) = W·T·(1 − exp(−t/(T·Q)))', () => {
  const r = ladderStepResponse({ Q: [5000], T: [0.4], w: [10] });
  rel(r.timeConstants[0], 2000, 1e-12);
  for (const t of [100, 2000, 10000]) rel(r.rise(t), 10 * 0.4 * (1 - Math.exp(-t / 2000)), 1e-12);
  rel(r.steady[0], 4, 1e-12);
});

test('two-node ladder: eigen-solution equals the IEC 60853-2 two-loop closed form; steady state = W·(TA + TB)', () => {
  const p = { TA: 0.39, TB: 0.052, QA: 6600, QB: 7100 };
  const two = twoLoopResponse(p);
  const lad = ladderStepResponse({ Q: [p.QA, p.QB], T: [p.TA, p.TB], w: [1, 0] });
  for (const t of [60, 600, 3600, 36000]) rel(lad.rise(t), two.rise(t), 1e-10);
  rel(two.Ta + two.Tb, p.TA + p.TB, 1e-12);
  rel(lad.steady[0], p.TA + p.TB, 1e-12);
  rel(lad.steady[1], p.TB, 1e-12);
  rel(lad.timeConstants[0], 1 / two.a, 1e-10); rel(lad.timeConstants[1], 1 / two.b, 1e-10);
  // initial slope of the conductor response is W/QA
  rel(lad.rise(1e-3) / 1e-3, 1 / p.QA, 1e-4);
});

test('three-node ladder with losses at every node: steady state by hand and energy conservation in the transient', () => {
  const Q = [4000, 6000, 9000]; const T = [0.3, 0.1, 0.05]; const w = [20, 5, 8];
  const r = ladderStepResponse({ Q, T, w });
  // steady: node 3 = 33·0.05, node 2 = + 25·0.1, node 1 = + 20·0.3
  rel(r.steady[2], 1.65, 1e-12); rel(r.steady[1], 4.15, 1e-12); rel(r.steady[0], 10.15, 1e-12);
  // stored heat + heat conducted out = heat injected (trapezoidal quadrature of the outflow θ3/T3)
  const tEnd = 5000; const n = 20000; let out = 0;
  for (let i = 0; i <= n; i++) { const t = (tEnd * i) / n; out += (i === 0 || i === n ? 0.5 : 1) * r.rise(t, 2) / T[2]; }
  out *= tEnd / n;
  const stored = Q[0] * r.rise(tEnd, 0) + Q[1] * r.rise(tEnd, 1) + Q[2] * r.rise(tEnd, 2);
  rel(stored + out, 33 * tEnd, 1e-6);
});

test('cable ladder from geometry: capacitances add up to the cable heat capacity less the outer part of the covering', () => {
  const lad = cableLadder(g);
  assert.equal(lad.Q.length, 2);
  const cap = (name) => { const l = g.layers.find((x) => x.name === name); return l.heatCapacity * l.area; };
  const p = vanWormer({ outerDiameter: g.dUnderSheath, innerDiameter: g.dc });
  rel(lad.Q[0], cap('conductor') + p * cap('insulation'), 1e-12);
  const pc = vanWormer({ outerDiameter: g.De, innerDiameter: g.sheathOuterDiameter });
  rel(lad.Q[1], (1 - p) * cap('insulation') + cap('sheath') + pc * cap('serving'), 1e-12);
  rel(cap('conductor'), 3.45e6 * 1000e-6, 1e-12);
  assert.deepEqual(lad.T, [g.T1, g.T3]);
});

test('soil step response tends to ρ/(2π)·ln(4L/De), the steady T4 with ln(2u)', () => {
  const big = soilStepResponse({ t: 1e14, rhoSoil: 1.2, diffusivity: 0.5e-6, depth: 1.5, De: 0.12 });
  rel(big, (1.2 / (2 * Math.PI)) * Math.log((4 * 1.5) / 0.12), 1e-5);
  rel(big, t4Buried({ rhoT: 1.2, depth: 1.5, De: 0.12 }), 2e-4);
  // early time: no image effect, θ = ρ/(4π)·E1(De²/(16δt))
  const t = 3600;
  rel(soilStepResponse({ t, rhoSoil: 1.2, diffusivity: 0.5e-6, depth: 1.5, De: 0.12 }), (1.2 / (4 * Math.PI)) * expIntegralE1(0.12 ** 2 / (16 * 0.5e-6 * t)), 1e-9);
  // group: mutual term tends to ρ/(2π)·ln(d′/d)
  const grp = soilStepResponse({ t: 1e14, rhoSoil: 1.2, diffusivity: 0.5e-6, depth: 1.5, De: 0.12, others: [{ d: 0.3, dImage: Math.hypot(0.3, 3) }] });
  rel(grp - big, (1.2 / (2 * Math.PI)) * Math.log(Math.hypot(0.3, 3) / 0.3), 1e-5);
});

test('complete transient response is monotonic and tends to the steady state T1 + T3 + ρ/(2π)·ln(4L/De)', () => {
  const m = transientModel({ geometry: g, external: buried });
  const steady = g.T1 + g.T3 + (1 / (2 * Math.PI)) * Math.log(4 / g.De);
  rel(m.steady, steady, 1e-12);
  let prev = 0;
  for (const t of [60, 600, 3600, 86400, 30 * 86400, 3.15e7, 3.15e9]) { const s = m.step(t); assert.ok(s > prev); prev = s; }
  rel(m.step(1e13), steady, 1e-4);
  // within 0.3 % of the Kennelly-based steady resistance used by the rating engine
  rel(m.steady, g.T1 + g.T3 + t4Buried({ rhoT: 1, depth: 1, De: g.De }), 3e-3);
  // very early: adiabatic heating of the conductor node, θ ≈ t/QA
  rel(m.step(1) / 1, 1 / m.ladder.Q[0], 0.01);
  const t63 = responseTime({ model: m });
  rel(m.step(t63), (1 - Math.exp(-1)) * m.steady, 1e-9);
  assert.ok(t63 > 3600 && t63 < 30 * 86400);
});

test('sheath and armour losses enter at their nodes: steady state = T1 + (1 + λ1)·T3 + (1 + λ1)·T4 for an unarmoured cable', () => {
  const m = transientModel({ geometry: g, lambda1: 0.3, external: { type: 'resistance', T4: 0.5 } });
  rel(m.steady, g.T1 + 1.3 * g.T3 + 1.3 * 0.5, 1e-12);
  rel(m.dielectricSteady, 0.5 * g.T1 + g.T3 + 0.5, 1e-12);
});

test('superposition: constant current from cold follows I²R·S(t); from a steady preload at the same current nothing changes', () => {
  const m = transientModel({ geometry: g, external: buried });
  const R = 2.8e-5; const I = 900; const dt = 1800; const N = 48;
  const cold = temperatureHistory({ model: m, current: new Array(N).fill(I), dt, R, ambientC: 15 });
  for (const k of [0, 10, 47]) rel(cold.conductorC[k] - 15, I * I * R * m.step((k + 1) * dt), 1e-12);
  const hot = temperatureHistory({ model: m, current: new Array(N).fill(I), dt, R, ambientC: 15, initialCurrent: I });
  for (const k of [0, 47]) rel(hot.conductorC[k], 15 + I * I * R * m.steady, 1e-12);
  // a load step up then down: the response is the difference of two shifted step responses
  const cur = new Array(N).fill(0); for (let k = 10; k < 20; k++) cur[k] = I;
  const pulse = temperatureHistory({ model: m, current: cur, dt, R, ambientC: 15 });
  rel(pulse.conductorC[30] - 15, I * I * R * (m.step(21 * dt) - m.step(11 * dt)), 1e-10);
  assert.equal(pulse.tMax, 20 * dt);
});

test('dynamic rating of an intermittent load: peak temperature stays well below the steady value for the peak current; statistics', () => {
  const m = transientModel({ geometry: g, external: buried });
  const R = 2.8e-5; const N = 24 * 60;
  // deterministic synthetic wind-like current: slow swings between 15 % and 100 % of 1400 A
  const current = Array.from({ length: N }, (_, k) => 1400 * (0.575 + 0.425 * Math.sin(k / 37) * Math.cos(k / 211)));
  const h = temperatureHistory({ model: m, current, dt: 3600, R, ambientC: 12 });
  const steadyPeak = 12 + 1400 * 1400 * R * m.steady;
  assert.ok(h.maxC < steadyPeak - 5, `max ${h.maxC} vs steady ${steadyPeak}`);
  const st = exceedanceStatistics({ conductorC: h.conductorC, limitC: h.maxC - 1e-9, dt: 3600 });
  rel(st.probability, 1 / N, 1e-9); assert.equal(st.durationAbove, 3600);
  assert.equal(exceedanceStatistics({ conductorC: h.conductorC, limitC: 90 }).probability, h.conductorC.filter((x) => x > 90).length / N);
  // temperature-dependent resistance feedback lowers the estimate (R given at 90 °C, cable cooler than that)
  const fb = temperatureHistory({ model: m, current, dt: 3600, R, ambientC: 12, feedback: true });
  assert.ok(fb.maxC < h.maxC && fb.maxC > 12);
});

test('cyclic rating factor: flat load gives M = 1; the direct numerical value agrees with the IEC 60853 six-ordinate formula within 2 %', () => {
  const m = transientModel({ geometry: g, external: buried });
  rel(cyclicRatingFactor({ model: m, cycle: new Array(24).fill(1) }).M, 1, 1e-9);
  // 12 h at peak, 12 h at half current; the conductor is hottest at the end of the peak period
  const cycle = Array.from({ length: 24 }, (_, hgt) => (hgt >= 8 && hgt < 20 ? 1 : 0.5));
  const num = cyclicRatingFactor({ model: m, cycle });
  rel(num.lossLoadFactor, (12 + 12 * 0.25) / 24, 1e-12);
  const iec = cyclicRatingFactorIEC({ model: m, Y: [1, 1, 1, 1, 1, 1], mu: num.lossLoadFactor });
  assert.ok(num.M > 1.03 && num.M < 1.2);
  rel(num.M, iec, 0.02);
  // more repetitions do not change the periodic state
  rel(cyclicRatingFactor({ model: m, cycle, nCycles: 60 }).M, num.M, 2e-3);
});

test('emergency rating by hand: I2 = sqrt(I1² + (θmax − θ1)/(R·S(t))); preload at the continuous rating leaves no overload margin', () => {
  const m = transientModel({ geometry: g, external: buried });
  const R = 2.8e-5;
  const e = emergencyRating({ model: m, R, preloadCurrent: 700, duration: 6 * 3600, thetaMaxC: 105, ambientC: 15 });
  const th1 = 15 + 700 * 700 * R * m.steady;
  rel(e.preloadC, th1, 1e-12);
  rel(e.I, Math.sqrt(700 * 700 + (105 - th1) / (R * m.step(6 * 3600))), 1e-12);
  assert.ok(e.ratio > 1);
  const full = emergencyRating({ model: m, R, preloadCurrent: e.steadyRating, duration: 3600, thetaMaxC: 105, ambientC: 15 });
  rel(full.I, e.steadyRating, 1e-9);
  // consistency with the history calculation: applying I2 after the preload reaches θmax at the end
  const hist = temperatureHistory({ model: m, current: new Array(12).fill(e.I), dt: 1800, R, ambientC: 15, initialCurrent: 700 });
  rel(hist.conductorC[11], 105, 1e-9);
});

test('adiabatic check hook: 630 mm² Cu from 90 °C carrying 90.14 kA for 1 s ends at 250 °C (IEC 60949)', () => {
  const I = 226 * 630 * Math.sqrt(Math.log((250 + 234.5) / (90 + 234.5)));
  rel(adiabaticFinalTemperature({ I, duration: 1, area: 630e-6, material: 'Cu', thetaInitialC: 90 }), 250, 1e-9);
});

test('2-D trench solver vs the image solution: cable-surface temperature within 0.5 % of Kennelly, heat balance closed', () => {
  const W = 50; const L = 1; const De = 0.2;
  const s = solveTrench2D({ cables: [{ x: 0, depth: L, diameter: De, W }], rhoSoil: 1, ambientC: 8, h: De / 4, refinement: 1, width: 400, height: 200 });
  const kennelly = W * t4Buried({ rhoT: 1, depth: L, De });
  rel(s.cables[0].surfaceMeanC - 8, kennelly, 5e-3);
  rel(s.cables[0].T4, t4Buried({ rhoT: 1, depth: L, De }), 5e-3);
  assert.ok(s.heatBalanceError < 1e-8);
  // a highly conductive cable interior approximates the isothermal surface assumed by Kennelly; the stair-step
  // boundary makes this first-order accurate (about 2 % at this grid), so the tolerance is 3 %
  const iso = solveTrench2D({ cables: [{ x: 0, depth: L, diameter: De, W, rhoT: 0.001 }], rhoSoil: 1, h: De / 4, refinement: 1, width: 400, height: 200 });
  rel(iso.cables[0].surfaceMeanC, kennelly, 0.03);
  assert.ok(iso.cables[0].maxC - iso.cables[0].meanC < 0.02 * kennelly);
});

test('2-D trench solver grid convergence: error against ρW/(2π)·ln(d′/d) at a probe decreases at about second order', () => {
  const W = 50; const L = 1; const De = 0.2;
  const probe = { x: 0.5, z: 0.6 };
  const exact = (W / (2 * Math.PI)) * Math.log(Math.hypot(0.5, 1.6) / Math.hypot(0.5, 0.4));
  const err = [0, 1, 2].map((refinement) => {
    const s = solveTrench2D({ cables: [{ x: 0, depth: L, diameter: De, W }], rhoSoil: 1, h: De / 4, refinement, probes: [probe], width: 800, height: 400 });
    return Math.abs(s.probes[0] - exact);
  });
  assert.ok(err[1] < err[0] && err[2] < err[1], `errors ${err}`);
  const order = Math.log2(err[0] / err[2]) / 2;
  assert.ok(order > 1.5 && order < 3.2, `observed order ${order}`);
  assert.ok(err[2] < 2e-3 * exact);
});

test('2-D trench solver: two cables superpose as in the image method; low-resistivity backfill and convective seabed behave physically', () => {
  const cables = [{ x: -0.25, depth: 1.0, diameter: 0.15, W: 40 }, { x: 0.25, depth: 1.0, diameter: 0.15, W: 40 }];
  const base = solveTrench2D({ cables, rhoSoil: 1.2, h: 0.15 / 6, refinement: 1, width: 400, height: 200, probes: [{ x: 0, z: 0.5 }] });
  const img = mutualHeatingRise({ rhoT: 1.2, target: { x: 0, depth: 0.5 }, sources: cables.map((c) => ({ x: c.x, depth: c.depth, W: c.W })) });
  rel(base.probes[0], img, 5e-3);
  rel(base.cables[0].surfaceMeanC, base.cables[1].surfaceMeanC, 1e-6);
  const single = (1.2 / (2 * Math.PI)) * 40 * (Math.log(4 / 0.15) + Math.log(Math.hypot(0.5, 2) / 0.5));
  rel(base.cables[0].surfaceMeanC, single, 5e-3);
  const bf = solveTrench2D({ cables, rhoSoil: 1.2, regions: [{ x0: -0.6, x1: 0.6, z0: 0.6, z1: 1.3, rhoT: 0.6 }], h: 0.15 / 6, refinement: 1, width: 400, height: 200 });
  assert.ok(bf.cables[0].surfaceMeanC < base.cables[0].surfaceMeanC);
  assert.ok(bf.heatBalanceError < 1e-8);
  // a finite surface coefficient adds resistance; a very large one recovers the isothermal surface
  const conv = solveTrench2D({ cables, rhoSoil: 1.2, surface: { type: 'convective', h: 5 }, h: 0.15 / 6, refinement: 1, width: 400, height: 200 });
  const convBig = solveTrench2D({ cables, rhoSoil: 1.2, surface: { type: 'convective', h: 1e7 }, h: 0.15 / 6, refinement: 1, width: 400, height: 200 });
  assert.ok(conv.cables[0].surfaceMeanC > base.cables[0].surfaceMeanC);
  rel(convBig.cables[0].surfaceMeanC, base.cables[0].surfaceMeanC, 1e-5);
});
