import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCatenary } from '../../src/engine/lines/catenary.js';
import { spreadLayout, brokenLineAnalysis } from '../../src/engine/lines/mooringsystem.js';
import {
  buildLumpedLine, relaxLine, simulateLine, detectSnapLoads, lineBreakTransient, airyWaveNumber,
} from '../../src/engine/lines/lumpedmass.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const rho = 1025, g = 9.80665;
const chain = { EA: 8e8, mass: 200, d: 0.18 };
const wChain = (chain.mass - rho * Math.PI * chain.d ** 2 / 4) * g;

test('static state (dynamic relaxation) converges to the elastic catenary at second order in segment count', () => {
  const cat = solveCatenary({ X: 800, Z: 300, L: 900, w: wChain, EA: chain.EA, seabed: false });
  const err = [8, 16, 32].map((n) => {
    const m = buildLumpedLine({ anchor: [0, 0, -300], fairlead: [800, 0, 0], sections: [{ L: 900, ...chain, nSeg: n }], seabed: false });
    const r = relaxLine(m, { tol: 1e-9 });
    assert.ok(r.converged);
    return { H: Math.abs(r.H / cat.H - 1), T: Math.abs(r.fairleadTension / cat.T - 1) };
  });
  assert.ok(err[2].H < 1e-3 && err[2].T < 1e-3);
  // observed order p = log2(e(h) / e(h/2)) ~ 2
  for (const k of ['H', 'T']) {
    const p1 = Math.log2(err[0][k] / err[1][k]), p2 = Math.log2(err[1][k] / err[2][k]);
    assert.ok(Math.abs(p1 - 2) < 0.15 && Math.abs(p2 - 2) < 0.1, `order ${k}: ${p1}, ${p2}`);
  }
});

test('static state with seabed contact matches the grounded catenary within discretisation error', () => {
  const cat = solveCatenary({ X: 830, Z: 100, L: 850, w: wChain, EA: chain.EA });
  const m = buildLumpedLine({ anchor: [0, 0, -100], fairlead: [830, 0, 0], sections: [{ L: 850, ...chain, nSeg: 40 }] });
  const r = relaxLine(m);
  assert.ok(r.converged);
  assert.ok(rel(r.fairleadTension, cat.T) < 2e-3 && rel(r.H, cat.H) < 2e-3);
  assert.ok(Math.abs(r.touchdownS - cat.groundedLength) < 850 / 40); // within one segment
  assert.ok(rel(r.anchorTension, cat.anchorT) < 5e-3);
});

test('uniform current on a taut neutrally buoyant line: midspan deflection = f L^2 / (8 T)', () => {
  const d = 0.1, L = 200, EA = 5e7, U = 1.2, Cdn = 1.2;
  const mass = rho * Math.PI * d * d / 4; // neutrally buoyant
  const m = buildLumpedLine({
    anchor: [0, 0, -50], fairlead: [L * 1.004, 0, -50], sections: [{ L, EA, mass, d, nSeg: 40, Cdn, Cdt: 0 }], seabed: false,
    current: { speed: U, direction: Math.PI / 2 },
  });
  const r = relaxLine(m, { tol: 1e-8 });
  assert.ok(r.converged);
  const f = 0.5 * rho * Cdn * d * U * U; // N/m, normal to the line
  const T = r.segmentTension[20];
  const yMid = m.r[3 * 20 + 1];
  assert.ok(yMid > 0.5);
  assert.ok(rel(yMid, f * (L * 1.004) ** 2 / (8 * T)) < 5e-3);
});

test('bending stiffness: tensioned beam under self weight, midspan sag vs closed form', () => {
  // y_mid = w L^2/(8T) - (w/(T k^2)) (1 - 1/cosh(k L/2)), k = sqrt(T/EI) (pinned ends; Timoshenko tie-rod)
  const L = 100, EA = 1e9, EI = 1e8, w = 500, mass = w / g;
  const build = (ei) => buildLumpedLine({
    anchor: [0, 0, 0], fairlead: [L * 1.001, 0, 0], sections: [{ L, EA, mass, d: 0.3, nSeg: 40, EI: ei }], seabed: false, rho: 0,
  });
  const mb = build(EI);
  const rb = relaxLine(mb, { tol: 1e-9 });
  assert.ok(rb.converged);
  const T = rb.segmentTension[20], k = Math.sqrt(T / EI), Ls = L * 1.001;
  const sagBeam = w * Ls * Ls / (8 * T) - (w / (T * k * k)) * (1 - 1 / Math.cosh(k * Ls / 2));
  assert.ok(rel(-mb.r[3 * 20 + 2], sagBeam) < 5e-3);
  // without EI the same model gives the cable sag w L^2 / (8 T)
  const mc = build(0);
  const rc = relaxLine(mc, { tol: 1e-9 });
  const sagCable = w * Ls * Ls / (8 * rc.segmentTension[20]);
  assert.ok(rel(-mc.r[3 * 20 + 2], sagCable) < 5e-3);
  assert.ok(-mb.r[3 * 20 + 2] < -mc.r[3 * 20 + 2] * 0.97);
});

test('taut line, slow axial fairlead oscillation: tension amplitude = EA a / L (quasi-static limit)', () => {
  const d = 0.1, L = 500, EA = 1e8, a = 0.5;
  const mass = rho * Math.PI * d * d / 4;
  const m = buildLumpedLine({ anchor: [0, 0, -100], fairlead: [L * 1.01, 0, -100], sections: [{ L, EA, mass, d, nSeg: 10, Cdt: 0, zeta: 0.05 }], seabed: false });
  const res = simulateLine(m, { motion: { type: 'harmonic', amplitude: [a, 0, 0], period: 60 }, duration: 180, dt: 0.25, statsFrom: 60, segments: [5] });
  assert.ok(rel(res.static.fairleadTension, EA * 0.01) < 1e-6);
  const amp = 0.5 * (res.stats.fairleadMax - res.stats.fairleadMin);
  assert.ok(rel(amp, EA * a / L) < 5e-3, `amplitude ${amp} vs ${EA * a / L}`);
  assert.ok(rel(res.stats.fairleadMean, EA * 0.01) < 5e-3);
  // quasi-static (catenary) tension follows the same amplitude -> DAF ~ 1
  assert.ok(Math.abs(res.stats.daf - 1) < 0.01);
  // the stored mid-line segment history carries the same amplitude
  let mx = -Infinity;
  for (let i = 240; i < res.t.length; i++) mx = Math.max(mx, res.segmentTension[0][i]);
  assert.ok(rel(mx - EA * 0.01, EA * a / L) < 5e-3);
});

const catenaryModel = (n, extra = {}) => buildLumpedLine({
  anchor: [0, 0, -100], fairlead: [830, 0, 0], sections: [{ L: 850, ...chain, nSeg: n, Cdn: 2.4, Cdt: 1.15, Can: 1, Cat: 0.5, zeta: 0.02 }], ...extra,
});

test('dynamic fairlead tension converges with time step (RK4) and with segment count; integrators agree', () => {
  const motion = { type: 'harmonic', amplitude: [2, 0, 1], period: 10 };
  const run = (n, o = {}) => simulateLine(catenaryModel(n), { motion, duration: 30, dt: 0.1, statsFrom: 10, quasiStatic: false, ...o }).stats.fairleadMax;
  // time-step refinement at fixed mesh (dtCritical ~ 0.01 s for 20 segments)
  const a = run(20, { dtInternal: 0.01 }), b = run(20, { dtInternal: 0.005 }), c = run(20, { dtInternal: 0.0025 });
  assert.ok(Math.abs(b - c) <= Math.abs(a - b) + 1e-9 * c);
  assert.ok(rel(a, c) < 1e-3);
  const si = run(20, { integrator: 'semi-implicit', dtInternal: 0.0025 });
  assert.ok(rel(si, c) < 5e-3);
  // mesh refinement: successive differences shrink
  const m10 = run(10), m20 = run(20), m40 = run(40);
  assert.ok(Math.abs(m40 - m20) < Math.abs(m20 - m10));
  assert.ok(rel(m20, m40) < 0.03);
  // dynamic tension exceeds the static value
  assert.ok(m40 > relaxLine(catenaryModel(40)).fairleadTension);
});

test('energy balance without dissipation: kinetic + strain + potential - fairlead work stays constant', () => {
  const m = buildLumpedLine({ anchor: [0, 0, -300], fairlead: [800, 0, 0], sections: [{ L: 900, ...chain, nSeg: 16, Cdn: 0, Cdt: 0 }], seabed: false, rho: 0 });
  const res = simulateLine(m, { motion: { type: 'harmonic', amplitude: [3, 0, 2], period: 12 }, duration: 36, dt: 0.1, dtInternal: 0.004, quasiStatic: false });
  const E = (i) => res.energy.kinetic[i] + res.energy.strain[i] + res.energy.potential[i] - res.energy.fairleadWork[i];
  let swing = 0, drift = 0;
  for (let i = 0; i < res.t.length; i++) {
    swing = Math.max(swing, Math.abs(res.energy.fairleadWork[i]));
    drift = Math.max(drift, Math.abs(E(i) - E(0)));
  }
  assert.ok(swing > 1e5);
  assert.ok(drift < 2e-3 * swing, `drift ${drift} vs work ${swing}`);
});

test('dynamic amplification and envelopes on a catenary line in waves and current', () => {
  const m = catenaryModel(20, {
    seabed: { mu: 0.5 }, current: { z: [-100, 0], speed: [0.2, 1.0], direction: 0 },
    wave: { H: 4, T: 9, direction: 0, depth: 100 },
  });
  const res = simulateLine(m, { motion: { type: 'harmonic', amplitude: [3, 0, 1.5], period: 9 }, duration: 45, dt: 0.1, statsFrom: 18 });
  assert.ok(res.stats.daf > 1); // drag and inertia amplify the wave-frequency tension
  assert.ok(res.stats.fairleadMax > res.static.fairleadTension && res.stats.fairleadMin < res.static.fairleadTension);
  assert.ok(res.stats.touchdownMax > res.stats.touchdownMin);
  assert.ok(res.envelopes.Tmax[19] >= res.envelopes.Tmin[19]);
  assert.ok(res.envelopes.zMax[10] > res.envelopes.zMin[10]);
  for (let i = 0; i < res.t.length; i++) assert.ok(Number.isFinite(res.fairleadTension[i]));
  // deep-water dispersion check: k = w^2/g
  assert.ok(rel(airyWaveNumber(2 * Math.PI / 5, 2000), (2 * Math.PI / 5) ** 2 / g) < 1e-10);
  // finite depth: w^2 = g k tanh(k h)
  const k = airyWaveNumber(0.5, 20);
  assert.ok(rel(g * k * Math.tanh(k * 20), 0.25) < 1e-12);
});

test('snap-load detection: slack then re-tension spike on a synthetic history and on a slack-taut simulation', () => {
  const t = [], T = [];
  for (let i = 0; i <= 400; i++) {
    const x = i * 0.05;
    t.push(x);
    // tension dips to zero around t = 5-6 s and rebounds to 2.6 x pretension
    T.push(x < 4 ? 100 + 20 * Math.sin(x) : x < 5 ? 100 * (5 - x) : x < 6 ? 0 : x < 6.5 ? 520 * (x - 6) : Math.max(100, 260 - 300 * (x - 6.5)));
  }
  const r = detectSnapLoads({ t, T, Tref: 100 });
  assert.equal(r.slack, true); assert.equal(r.snap, true); assert.equal(r.events.length, 1);
  assert.ok(Math.abs(r.events[0].peak - 260) < 1e-9 && Math.abs(r.events[0].ratio - 2.6) < 1e-9);
  assert.ok(Math.abs(r.events[0].tSlack - 4.9) < 0.051 && r.events[0].slackDuration > 1);
  assert.equal(detectSnapLoads({ t, T: T.map(() => 100), Tref: 100 }).slack, false);
  // taut light line driven hard along its axis goes slack each cycle and snaps back
  const d = 0.05, L = 100, EA = 2e7, mass = 8;
  const m = buildLumpedLine({ anchor: [0, 0, -80], fairlead: [L * 1.002, 0, -80], sections: [{ L, EA, mass, d, nSeg: 10, zeta: 0.05 }], seabed: false });
  const res = simulateLine(m, { motion: { type: 'harmonic', amplitude: [1.0, 0, 0], period: 4, ramp: 4 }, duration: 16, dt: 0.01, quasiStatic: false });
  const sn = detectSnapLoads({ t: res.t, T: res.fairleadTension, Tref: res.static.fairleadTension });
  assert.equal(sn.slack, true);
  assert.equal(sn.snap, true);
  assert.ok(sn.minTension <= sn.threshold);
});

test('line-break transient of a spread-moored floater: overshoot consistent with the energy method', () => {
  const sys = spreadLayout({ nLines: 8, anchorRadius: 880, fairleadRadius: 30, fairleadZ: -10, depth: 150, segments: [{ L: 900, w: 1500, EA: 7e8 }] });
  const Fext = [1.5e6, 0, 0, 0, 0, 0];
  const und = lineBreakTransient({ system: sys, brokenLine: 4, Fext, mass: [3e7, 3e7, 2e10], duration: 400, dt: 0.5 });
  assert.ok(und.staticExcursion > 0.5);
  const en = brokenLineAnalysis(sys, { Fext, lines: [4] }).damaged[0];
  // undamped time-domain peak excursion vs work-energy estimate (motion is along x by symmetry)
  assert.ok(rel(und.overshootFactor, en.overshootFactor) < 0.02, `${und.overshootFactor} vs ${en.overshootFactor}`);
  assert.ok(rel(und.peakTension, en.maxTransientTension) < 0.02);
  // with damping the overshoot is smaller and the floater settles at the damaged equilibrium
  const dmp = lineBreakTransient({ system: sys, brokenLine: 4, Fext, mass: [3e7, 3e7, 2e10], damping: [8e5, 8e5, 1e9], duration: 1500, dt: 0.5 });
  assert.ok(dmp.overshootFactor < und.overshootFactor && dmp.overshootFactor > 1);
  const n = dmp.t.length - 1;
  assert.ok(Math.abs(dmp.x[n] - dmp.damaged.pose[0]) < 0.02 * dmp.staticExcursion);
  // before the break nothing moves
  const late = lineBreakTransient({ system: sys, brokenLine: 4, Fext, mass: [3e7, 3e7, 2e10], t0: 50, duration: 60, dt: 0.5 });
  assert.ok(late.offset[90] < 1e-6 && late.offset[119] > 1e-3);
});

test('invalid input throws RangeError', () => {
  assert.throws(() => buildLumpedLine({ anchor: [0, 0, 0], fairlead: [1, 0, 0], sections: [] }), /sections/);
  assert.throws(() => buildLumpedLine({ anchor: [0, 0, 0], fairlead: [1, 0, 0], sections: [{ L: 1, EA: 1, mass: 1, d: 1, nSeg: 0 }] }), /nSeg/);
  const m = buildLumpedLine({ anchor: [0, 0, 0], fairlead: [10, 0, 5], sections: [{ L: 12, EA: 1e6, mass: 5, d: 0.05, nSeg: 4 }], seabed: false });
  assert.throws(() => simulateLine(m, { duration: -1, dt: 0.1 }), /duration/);
  assert.throws(() => simulateLine(m, { duration: 1, dt: 0.1, motion: { type: 'spin' } }), /motion/);
});
