import test from 'node:test';
import assert from 'node:assert/strict';
import {
  makeGrid, flowlineRiserProfile, simulateTransient, analyseOscillation, detectSlugArrivals, rampUp, shutdownRestart,
  severeSluggingCase, trackSlugs,
} from '../../src/engine/flow/transient.js';
import { marchPipeline } from '../../src/engine/flow/multiphase.js';
import { G } from '../../src/engine/flow/driftflux.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
const gasRT = 287 * 293; // air at 20 °C
const fluid = { rhoL: 1000, aL: 300, gasRT };

test('grid: uniform cells along the route with the inclination of each segment', () => {
  const g = makeGrid([{ x: 0, z: 0 }, { x: 30, z: 0 }, { x: 30, z: 10 }], 40);
  close(g.length, 40, 1e-14); close(g.ds, 1, 1e-14);
  assert.equal(g.sinTheta[10], 0); close(g.sinTheta[35], 1, 1e-14);
  close(g.z[39], 9.5, 1e-13); close(g.x[0], 0.5, 1e-13);
  const p = flowlineRiserProfile({ flowlineLength: 100, flowlineAngleDeg: -3, riserHeight: 20 });
  close(p[1].z, -100 * Math.sin(3 * Math.PI / 180), 1e-13);
  close(p[2].z - p[1].z, 20, 1e-13); close(p[2].x, p[1].x, 1e-13);
});

test('mass conservation to round-off with closed ends (sloshing in a V-shaped closed pipe)', () => {
  const r = simulateTransient({
    profile: [{ x: 0, z: 5 }, { x: 20, z: 0 }, { x: 40, z: 8 }], D: 0.1, nCells: 60, tEnd: 6, fluid,
    inlet: { type: 'closed' }, outlet: { type: 'closed' }, initial: { type: 'uniform', P: 2e5, alphaG: 0.5 },
  });
  assert.ok(r.steps > 1000);
  assert.ok(Math.abs(r.conservation.residualL) < 1e-12, `liquid residual ${r.conservation.residualL}`);
  assert.ok(Math.abs(r.conservation.residualG) < 1e-12, `gas residual ${r.conservation.residualG}`);
  assert.equal(r.conservation.netInflowL, 0); assert.equal(r.conservation.netInflowG, 0);
  assert.ok(r.conservation.clippedL < 1e-12 * r.conservation.massL0 && r.conservation.clippedG < 1e-12 * r.conservation.massG0);
  // the liquid has moved towards the low point and the gas towards the high ends
  const h = r.final.holdup;
  assert.ok(h[30] > 0.6 && h[59] < 0.2 && h[0] < 0.4, `holdup low point ${h[30]}, high end ${h[59]}`);
  for (const a of r.final.alphaG) assert.ok(a >= 0 && a <= 1);
});

test('mass conservation to round-off with through-flow, a gas-lift source and a choke outlet', () => {
  const r = simulateTransient({
    profile: [{ x: 0, z: 0 }, { x: 100, z: 5 }], D: 0.1, nCells: 40, tEnd: 20, fluid, roughness: 1e-5, slip: { C0: 1.2, vd: 0.2 },
    inlet: { type: 'massflow', massRateL: (t) => 4 + 2 * Math.sin(t), massRateG: { t: [0, 10, 20], v: [0.01, 0.03, 0.02] } },
    outlet: { type: 'choke', pressure: 2e5, opening: 0.3 }, sources: [{ s: 50, massRateG: 0.01 }],
    initial: { type: 'steady', massRateL: 4, massRateG: 0.01, pOut: 2.5e5 },
  });
  assert.ok(Math.abs(r.conservation.residualL) < 1e-12, `liquid residual ${r.conservation.residualL}`);
  assert.ok(Math.abs(r.conservation.residualG) < 1e-12, `gas residual ${r.conservation.residualG}`);
  assert.ok(r.conservation.netInflowG !== 0 && r.conservation.netInflowL !== 0);
  // choke: the pipe-side pressure sits above the downstream pressure by the orifice drop ρ v²/(2 z²)
  assert.ok(r.final.P[39] > 2e5);
});

test('hydrostatic equilibrium is preserved: settled liquid below gas in a closed vertical pipe stays at rest', () => {
  const r = simulateTransient({
    profile: [{ x: 0, z: 0 }, { x: 0, z: 20 }], D: 0.1, nCells: 40, tEnd: 2, fluid,
    inlet: { type: 'closed' }, outlet: { type: 'closed' }, initial: { type: 'settled', liquidFraction: 0.5, pRef: 2e5 },
  });
  assert.ok(r.steps > 2000);
  const p0 = r.profiles.P[0];
  let vmax = 0;
  let dp = 0;
  for (let j = 0; j < 40; j++) {
    if (r.final.holdup[j] > 0.5) vmax = Math.max(vmax, Math.abs(r.final.vL[j])); else vmax = Math.max(vmax, Math.abs(r.final.vG[j]));
    dp = Math.max(dp, Math.abs(r.final.P[j] - p0[j]) / p0[j]);
  }
  assert.ok(vmax < 1e-8, `max velocity ${vmax} m/s`);
  assert.ok(dp < 1e-10, `max relative pressure drift ${dp}`);
  // the initial state is the hand-calculated hydrostatic column: bottom pressure = top + gas head + liquid head
  const bottom = p0[0] + 0.5 * G * 0.5 * (1000 + (p0[0] - 101325) / 9e4);
  assert.ok(Math.abs(bottom - (2e5 + 1000 * G * 10)) / bottom < 2e-3, `bottom pressure ${bottom}`);
  close(r.history.liquidInventory[r.history.t.length - 1], Math.PI * 0.01 / 4 * 10, 1e-9);
});

test('gas-fraction front with no slip (C0 = 1, vd = 0) travels at the mixture velocity', () => {
  const L = 100; const v = 2; const tEnd = 20; const N = 200;
  const r = simulateTransient({
    profile: [{ x: 0, z: 0 }, { x: L, z: 0 }], D: 0.1, nCells: N, tEnd, fluid, slip: { C0: 1, vd: 0 }, friction: false,
    inlet: { type: 'periodic' }, outlet: { type: 'periodic' },
    initial: { type: 'profile', P: 5e5, alphaG: (s) => (s > 20 && s < 40 ? 0.6 : 0.2), vL: v },
  });
  // the block [20, 40] m must be at [60, 80] m: locate the mid-height crossings of the holdup profile
  const h = r.final.holdup;
  const s = r.grid.s;
  const cross = [];
  for (let j = 1; j < N; j++) {
    if ((h[j - 1] - 0.6) * (h[j] - 0.6) < 0) cross.push(s[j - 1] + (0.6 - h[j - 1]) / (h[j] - h[j - 1]) * (s[j] - s[j - 1]));
  }
  assert.equal(cross.length, 2);
  assert.ok(Math.abs(cross[0] - 60) < 0.5, `rear front at ${cross[0]} m`);
  assert.ok(Math.abs(cross[1] - 80) < 0.5, `leading front at ${cross[1]} m`);
  // pressure and velocity stay uniform across the contact; no over/undershoot of the volume fraction
  for (let j = 0; j < N; j++) {
    assert.ok(Math.abs(r.final.P[j] - 5e5) < 1e-3 && Math.abs(r.final.vL[j] - v) < 1e-9);
    assert.ok(h[j] > 0.4 - 1e-9 && h[j] < 0.8 + 1e-9);
  }
  assert.ok(Math.abs(r.conservation.residualL) < 1e-12 && Math.abs(r.conservation.residualG) < 1e-12);
});

test('observed convergence order on a smooth advected pulse: ≈ 2 (MUSCL + SSP-RK2), ≈ 1 first order', () => {
  const L = 10; const v = 1; const tEnd = 2;
  const a0 = (s) => 0.3 + 0.2 * Math.exp(-(((s - 4) / 0.8) ** 2));
  const error = (N, limiter) => {
    const r = simulateTransient({
      profile: [{ x: 0, z: 0 }, { x: L, z: 0 }], D: 0.1, nCells: N, tEnd, fluid, slip: { C0: 1, vd: 0 }, friction: false, limiter,
      inlet: { type: 'periodic' }, outlet: { type: 'periodic' }, initial: { type: 'profile', P: 2e5, alphaG: a0, vL: v },
    });
    let e = 0;
    for (let j = 0; j < N; j++) e += Math.abs(r.final.alphaG[j] - a0(((r.grid.s[j] - v * tEnd) % L + L) % L));
    return e / N;
  };
  const second = [40, 80, 160].map((N) => error(N, 'none'));
  const order = Math.log2(second[1] / second[2]);
  assert.ok(order > 1.8 && order < 2.3, `unlimited MUSCL order ${order}`);
  assert.ok(Math.log2(second[0] / second[1]) > 1.8);
  const limited = [40, 80].map((N) => error(N, 'vanleer'));
  assert.ok(Math.log2(limited[0] / limited[1]) > 1.5, 'van Leer limiter order');
  const first = [40, 80].map((N) => error(N, 'first'));
  const o1 = Math.log2(first[0] / first[1]);
  assert.ok(o1 > 0.5 && o1 < 1.2, `first-order scheme order ${o1}`);
  assert.ok(second[1] < first[1] / 5);
});

test('steady state of the transient solver agrees with the steady drift-flux marching solution', () => {
  const D = 0.1; const mL = 6; const mG = 0.03; const pOut = 3e5;
  const profile = [{ x: 0, z: 0 }, { x: 300, z: 12 }];
  const r = simulateTransient({
    profile, D, nCells: 60, tEnd: 400, fluid, slip: { C0: 1.2, vd: 0.3 }, roughness: 1e-5,
    inlet: { type: 'massflow', massRateL: mL, massRateG: mG }, outlet: { type: 'pressure', pressure: pOut },
    initial: { type: 'uniform', P: pOut, alphaG: 0.5 },
  });
  // same equation of state, friction law and slip closure in the independent marching solver
  const fluidAt = (P) => ({ rhoL: 1000 + (P - 101325) / 300 ** 2, rhoG: P / gasRT, muL: 1e-3, muG: 1.8e-5, sigma: 0.07 });
  const m = marchPipeline({ profile, D, roughness: 1e-5, fluidAt, temperature: 293, pOut, massRateL: mL, massRateG: mG, method: 'drift-flux', closure: 'slug', C0: 1.2, vd: 0.3, maxStep: 1 });
  const at = (s, key) => {
    let k = 1;
    while (m.stations[k].s < s) k++;
    const a = m.stations[k - 1]; const b = m.stations[k];
    return a[key] + (b[key] - a[key]) * (s - a.s) / (b.s - a.s);
  };
  // steady: outlet rates equal inlet rates, inlet pressure no longer changes
  const n = r.history.t.length;
  close(r.history.outletLiquidRate[n - 1], mL, 1e-6);
  close(r.history.outletGasRate[n - 1], mG, 1e-6);
  close(r.history.inletPressure[n - 1], r.history.inletPressure[n - 20], 1e-8);
  const dpRef = at(r.grid.s[0], 'P') - pOut;
  assert.ok(dpRef > 1e5);
  // (the first and last cells carry the first-order boundary treatment and are compared on pressure only)
  for (const j of [0, 59]) assert.ok(Math.abs(r.final.P[j] - at(r.grid.s[j], 'P')) < 0.01 * dpRef, `P at boundary cell ${j}`);
  for (const j of [5, 15, 30, 45, 55]) {
    assert.ok(Math.abs(r.final.P[j] - at(r.grid.s[j], 'P')) < 0.01 * dpRef, `P at cell ${j}: ${r.final.P[j]} vs ${at(r.grid.s[j], 'P')}`);
    assert.ok(Math.abs(r.final.holdup[j] - at(r.grid.s[j], 'holdup')) < 2e-3, `holdup at cell ${j}`);
    close(r.final.vL[j], at(r.grid.s[j], 'vL'), 5e-3);
    close(r.final.vG[j], at(r.grid.s[j], 'vG'), 5e-3);
  }
  // the 'steady' initial condition is that same state: starting from it nothing moves
  const s2 = simulateTransient({
    profile, D, nCells: 60, tEnd: 5, fluid, slip: { C0: 1.2, vd: 0.3 }, roughness: 1e-5,
    inlet: { type: 'massflow', massRateL: mL, massRateG: mG }, outlet: { type: 'pressure', pressure: pOut },
    initial: { type: 'steady', massRateL: mL, massRateG: mG, pOut },
  });
  assert.ok(Math.abs(s2.profiles.P[0][0] - r.final.P[0]) < 0.01 * dpRef);
});

const D = 0.05;
const A = Math.PI * D * D / 4;
// Laboratory-scale flowline–riser loop (cf. Schmidt, Brill & Beggs 1980; Taitel 1986): 60 m flowline at −2°,
// 10 m vertical riser, 50 mm bore, air–water, atmospheric separator, vsl = 0.1 m/s.
const loop = {
  flowlineLength: 60, flowlineAngleDeg: -2, riserHeight: 10, D, massRateL: 1000 * 0.1 * A, pSep: 1e5,
  fluid: { rhoL: 1000, aL: 100, gasRT }, slip: { C0: 1.2, vd: 'buoyant' }, nCells: 42, nHistory: 400,
};

test('severe slugging: classic unstable case (low gas rate, Bøe criterion violated) gives a cyclic build-up / blow-out', () => {
  const r = severeSluggingCase({ ...loop, massRateG: 1.2 * 0.4 * A, tEnd: 420, discardFraction: 0.25 });
  assert.equal(r.boe.severeSlugging, true, `ΠSS = ${r.boe.piSS}`);
  assert.equal(r.taitel.stable, false);
  assert.equal(r.slugging, true);
  // riser-base pressure swings over a large part of the full-riser liquid head (98 kPa)
  assert.ok(r.headAmplitude > 0.4, `head amplitude ${r.headAmplitude}`);
  assert.ok(r.pressure.max > 1.6e5 && r.pressure.max < 1e5 + 1000 * G * 10 * 1.1, `peak pressure ${r.pressure.max}`);
  // liquid production is intermittent: long periods without production, then a surge many times the mean rate
  assert.ok(r.liquid.max > 5 * loop.massRateL, `peak liquid rate ${r.liquid.max}`);
  assert.ok(r.liquid.min < 0.1 * loop.massRateL);
  assert.ok(r.pressure.cycles >= 2, `${r.pressure.cycles} cycles`);
  assert.ok(r.period > 30 && r.period < 200, `period ${r.period} s`);
  const surges = detectSlugArrivals({ t: r.history.t, rate: r.history.outletLiquidRate, thresholdFactor: 3, discardFraction: 0.25 });
  assert.ok(surges.count >= 2 && surges.surgeFraction > 0.5);
  assert.ok(Math.abs(r.conservation.residualL) < 1e-11 && Math.abs(r.conservation.residualG) < 1e-11);
  // mean production equals the inlet rate over whole cycles (within the cycle-phase error of a finite record)
  close(r.liquid.mean, loop.massRateL, 0.35);
});

test('severe slugging: stable flow when the Bøe criterion says stable (high gas rate)', () => {
  const r = severeSluggingCase({ ...loop, massRateG: 1.2 * 3 * A, tEnd: 90 });
  assert.equal(r.boe.severeSlugging, false, `ΠSS = ${r.boe.piSS}`);
  assert.ok(r.boe.piSS > 2);
  assert.equal(r.slugging, false);
  assert.ok(r.headAmplitude < 0.01, `head amplitude ${r.headAmplitude}`);
  assert.ok(Number.isNaN(r.period));
  close(r.liquid.max, loop.massRateL, 0.01);
  close(r.liquid.min, loop.massRateL, 0.01);
});

test('severe slugging mitigation: riser-base gas lift and topside choking suppress the cycle', () => {
  const lift = severeSluggingCase({ ...loop, massRateG: 1.2 * 0.4 * A, tEnd: 260, gasLiftRate: 1.2 * 3 * A });
  assert.equal(lift.slugging, false);
  assert.ok(lift.headAmplitude < 0.02, `gas-lift head amplitude ${lift.headAmplitude}`);
  assert.ok(lift.liquid.max < 1.2 * loop.massRateL);
  const choke = severeSluggingCase({ ...loop, massRateG: 1.2 * 0.4 * A, tEnd: 260, chokeOpening: 0.03 });
  assert.ok(choke.headAmplitude < 0.2, `choked head amplitude ${choke.headAmplitude}`);
  assert.ok(choke.liquid.max < 2 * loop.massRateL, `choked peak liquid rate ${choke.liquid.max}`);
  // the price of choking is back-pressure
  assert.ok(choke.pressure.mean > 1.5e5);
});

test('ramp-up and shutdown/restart scenarios: liquid swept out on ramp-up, zero flow when shut in, recovery after restart', () => {
  const base = { profile: [{ x: 0, z: 0 }, { x: 150, z: 6 }], D: 0.1, nCells: 30, fluid, slip: { C0: 1.2, vd: 0.3 }, roughness: 1e-5, pOut: 3e5 };
  const up = rampUp({ ...base, tEnd: 300, from: { massRateL: 3, massRateG: 0.02 }, to: { massRateL: 3, massRateG: 0.08 }, rampStart: 20, rampDuration: 20 });
  const n = up.history.t.length;
  close(up.history.outletLiquidRate[n - 1], 3, 1e-3);
  close(up.history.outletGasRate[n - 1], 0.08, 1e-3);
  // more gas → lower holdup: liquid is swept out, so the outlet liquid rate overshoots during the transient
  assert.ok(up.liquidSwept > 0);
  assert.ok(up.surge.peakRate > 3.05);
  let cum = 0;
  for (let i = 0; i < n; i++) cum += (up.history.outletLiquidRate[i] - 3) * (300 / n);
  close(cum / 1000, up.liquidSwept, 0.02);
  const sd = shutdownRestart({ ...base, tEnd: 300, massRateL: 3, massRateG: 0.05, tShutdown: 20, tRestart: 150, valveTime: 10 });
  const t = sd.history.t;
  const iShut = t.findIndex((x) => x > 140);
  assert.ok(Math.abs(sd.history.outletLiquidRate[iShut]) < 0.02 * 3, `liquid still flowing when shut in: ${sd.history.outletLiquidRate[iShut]}`);
  assert.ok(sd.settledInventory > 0);
  close(sd.history.outletLiquidRate[t.length - 1], 3, 0.02);
  assert.ok(Math.abs(sd.conservation.residualL) < 1e-12);
});

test('boundary specification is checked: over-constrained and inconsistent combinations throw', () => {
  const base = { profile: [{ x: 0, z: 0 }, { x: 10, z: 0 }], D: 0.1, nCells: 10, tEnd: 1, fluid, initial: { type: 'uniform', P: 1e5, alphaG: 0.5 } };
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'massflow', massRateL: 1, massRateG: 0.01, pressure: 2e5 }, outlet: { type: 'pressure', pressure: 1e5 } }), /Over-constrained inlet/);
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'closed' }, outlet: { type: 'pressure', pressure: 1e5, massRateL: 1 } }), /Over-constrained outlet/);
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'periodic' }, outlet: { type: 'closed' } }), /periodic/);
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'closed' }, outlet: { type: 'closed' }, slip: { C0: 1.4, alphaBlend: 0.75 } }), /slip/);
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'closed' }, outlet: { type: 'closed' }, limiter: 'superbee' }), /limiter/);
  assert.throws(() => simulateTransient({ ...base, inlet: { type: 'closed' }, outlet: { type: 'closed' }, initial: { type: 'uniform', P: 1e5, alphaG: 1.5 } }), /alphaG/);
});

test('surge detection and oscillation statistics on synthetic signals', () => {
  const n = 2000;
  const t = Float64Array.from({ length: n }, (_, i) => (i + 0.5) * 0.1);
  // square wave: 2 s surges of rate 10 every 20 s on a base rate of 1
  const rate = Float64Array.from(t, (x) => (x % 20 < 2 ? 10 : 1));
  const d = detectSlugArrivals({ t, rate, thresholdFactor: 2 });
  assert.equal(d.count, 10);
  close(d.meanPeriod, 20, 1e-9);
  close(d.maxAmount, 20, 1e-9);
  close(d.meanRate, 1.9, 1e-9);
  close(d.surgeFraction, 200 / 380, 1e-9);
  assert.ok(d.periodCov < 1e-9);
  const y = Float64Array.from(t, (x) => 5 + 2 * Math.sin(2 * Math.PI * x / 25));
  const o = analyseOscillation({ t, y, discardFraction: 0.25 });
  close(o.period, 25, 5e-3); close(o.mean, 5, 2e-2); close(o.relativeAmplitude, 0.8, 2e-2);
  assert.ok(o.cycles >= 4 && o.cycles <= 6);
});

test('Lagrangian slug tracking: prescribed moving fronts, growth/decay, arrival and merging are exact', () => {
  // vF = (HLS vm − HLF vLF)/(HLS − HLF) = (0.9 × 3 − 0.3 × 0.5)/0.6 = 4.25 m/s; tail 3.8 m/s → growth 0.45 m/s
  const base = { length: 500, D: 0.1, vm: 3, vt: 3.8, slugHoldup: 0.9, filmHoldup: 0.3, vLF: 0.5 };
  const one = trackSlugs({ ...base, slugs: [{ tail: 0, front: 5 }], tEnd: 40, dt: 0.5 });
  close(one.frontVelocity, 4.25, 1e-13);
  close(one.slugs[0].front, 5 + 4.25 * 40, 1e-12);
  close(one.slugs[0].tail, 3.8 * 40, 1e-12);
  close(one.slugs[0].length, 5 + 0.45 * 40, 1e-12);
  // time-step refinement does not change the result (velocities are constant between events)
  close(trackSlugs({ ...base, slugs: [{ tail: 0, front: 5 }], tEnd: 40, dt: 0.05 }).slugs[0].length, 23, 1e-11);
  // arrival at the outlet: front reaches 500 m at t = 495/4.25, with length 5 + 0.45 t
  const arr = trackSlugs({ ...base, slugs: [{ tail: 0, front: 5 }], tEnd: 140, dt: 0.5 });
  close(arr.arrivals[0].t, 495 / 4.25, 1e-12);
  close(arr.arrivals[0].length, 5 + 0.45 * 495 / 4.25, 1e-12);
  assert.equal(arr.slugs.length, 0);
  // decay over a thin film: vF = (2.7 − 0.025)/0.85 = 3.147 < 3.8 → a 5 m slug dies after 5/0.653 = 7.66 s
  const thin = trackSlugs({ ...base, filmHoldup: 0.05, slugs: [{ tail: 0, front: 5 }], tEnd: 20, dt: 0.01 });
  assert.equal(thin.dissipated, 1);
  const tDeath = thin.t[thin.count.findIndex((c) => c === 0)];
  close(tDeath, 5 / (3.8 - 2.675 / 0.85), 2e-3);
  // merging: a short slug with a wake-accelerated tail is caught by the slug behind it
  const two = trackSlugs({ ...base, slugs: [{ tail: 20, front: 21 }, { tail: 0, front: 10 }], wakeA: 0, tEnd: 30, dt: 0.01 });
  // gap 10 m closes at vF − vt = 0.45 m/s → merge at t = 22.2 s into one slug of the combined extent
  assert.equal(two.merges.length, 1);
  close(two.merges[0].t, 10 / 0.45, 1e-3);
  assert.equal(two.slugs.length, 1);
  close(two.slugs[0].length, 21 + 4.25 * 30 - 3.8 * 30, 1e-9);
  // periodic slug generation at the inlet
  const gen = trackSlugs({ ...base, inletPeriod: 10, inletSlugLength: 3, tEnd: 100, dt: 0.1 });
  assert.ok(gen.count[gen.count.length - 1] >= 9);
});
