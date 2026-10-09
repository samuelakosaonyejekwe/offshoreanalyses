import test from 'node:test';
import assert from 'node:assert/strict';
import * as td from '../../src/engine/hydro/timedomain.js';
import { polyRoots } from '../../src/engine/math/roots.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);

// Analytic radiation model: K(t) = k0 exp(-a t)  <=>  B(w) = k0 a/(a² + w²),  A(w) = Ainf - k0/(a² + w²)
const k0 = 4e5; const a = 0.8; const Ainf = 6e5;
const Bw = (w) => (k0 * a) / (a * a + w * w);
const Aw = (w) => Ainf - k0 / (a * a + w * w);

test('Retardation kernel by cosine transform reproduces K(t) = k0 exp(-a t) from B(w) = k0 a/(a² + w²)', () => {
  // (the analytic B does not vanish at w = 0, so the table starts very close to zero)
  const omega = []; for (let w = 1e-4; w <= 400; w *= 1.01) omega.push(w);
  const ker = td.retardationKernel({ omega, B: omega.map(Bw), dt: 0.1, tMax: 10 });
  for (const l of [0, 5, 10, 30, 60, 100]) near(ker.K[l][0], k0 * Math.exp(-a * ker.t[l]), 0.004 * k0, `t = ${ker.t[l]}`);
  // Ogilvie relation recovers the infinite-frequency added mass from A(w) and K(t)
  const long = td.retardationKernel({ omega, B: omega.map(Bw), dt: 0.05, tMax: 25 });
  const ws = [0.4, 0.8, 1.2, 2.0];
  const est = td.infiniteFrequencyAddedMass({ omega: ws, A: ws.map(Aw), kernel: long });
  close(est.Ainf, Ainf, 5e-3); assert.ok(est.spread < 0.01);
});

test('Cummins free decay reproduces the frequency-domain natural period (complex pole of the transfer function)', () => {
  const M = 1.4e6; const C = 2.0e6; const dt = 0.02;
  // frequency domain: [-w²(M + A(w)) + i w B(w) + C] = 0 with s = i w  <=>  (M + Ainf) s²(s + a) + k0 s + C (s + a) = 0
  const mt = M + Ainf;
  const roots = polyRoots([mt, mt * a, k0 + C, C * a]);
  let wd = 0; let sig = 0;
  roots.im.forEach((im, i) => { if (im > 1e-9) { wd = im; sig = roots.re[i]; } });
  assert.ok(wd > 0 && sig < 0);
  const nLag = Math.round(30 / dt) + 1;
  const kernel = { dt, n: 1, K: Array.from({ length: nLag }, (_, l) => Float64Array.of(k0 * Math.exp(-a * l * dt))) };
  const r = td.freeDecay({ mass: mt, stiffness: C, kernel, x0: 1, dt, duration: 60 });
  const da = td.decayAnalysis({ t: r.t, x: r.x });
  close(da.period, (2 * Math.PI) / wd, 2e-3, 'damped period');
  close(da.logDecrement, (-sig * 2 * Math.PI) / wd, 0.03, 'log decrement');
  // consistent with the undamped frequency-domain estimate w² (M + A(w)) = C within the damping shift
  let w = 1; for (let i = 0; i < 50; i++) w = Math.sqrt(C / (M + Aw(w)));
  close(da.period, (2 * Math.PI) / w, 0.02, 'frequency-domain natural period');
});

test('Cummins with a kernel computed from tabulated B(w) gives the same decay as the analytic kernel', () => {
  const omega = []; for (let w = 0.02; w <= 400; w *= 1.01) omega.push(w);
  const dt = 0.05;
  const ker = td.retardationKernel({ omega, B: omega.map(Bw), dt, tMax: 20 });
  const mt = 1.4e6 + Ainf;
  const r = td.freeDecay({ mass: mt, stiffness: 2e6, kernel: ker, x0: 1, dt, duration: 40 });
  const roots = polyRoots([mt, mt * a, k0 + 2e6, 2e6 * a]);
  const wd = Math.max(...roots.im);
  close(td.decayAnalysis({ t: r.t, x: r.x }).period, (2 * Math.PI) / wd, 4e-3);
});

test('Constant-coefficient fallback: damped period and damping ratio of a linear oscillator', () => {
  const m = 2e6; const c = 5e6; const zeta = 0.05; const b = 2 * zeta * Math.sqrt(m * c);
  const wn = Math.sqrt(c / m);
  const r = td.freeDecay({ mass: m, stiffness: c, b1: b, x0: 2, dt: 0.01, duration: 60 });
  const da = td.decayAnalysis({ t: r.t, x: r.x, inertia: m });
  close(da.period, (2 * Math.PI) / (wn * Math.sqrt(1 - zeta * zeta)), 1e-4);
  close(da.dampingRatio, zeta, 1e-3);
  close(da.b1, b, 5e-3); near(da.p2, 0, 1e-3);
  // exact solution at t = 10 s
  const wdd = wn * Math.sqrt(1 - zeta * zeta); const t = 10;
  const exact = 2 * Math.exp(-zeta * wn * t) * (Math.cos(wdd * t) + ((zeta * wn) / wdd) * Math.sin(wdd * t));
  near(r.x[1000], exact, 2e-4);
});

test('Decay analysis separates linear and quadratic damping (p–q fit)', () => {
  const m = 1e6; const c = 1e6; const b1 = 4e4; const b2 = 3e5;
  // initial equivalent damping ratio (b1 + 8/(3 pi) b2 w X0)/(2 sqrt(m c)) = 0.08: light damping as the method assumes
  const r = td.freeDecay({ mass: m, stiffness: c, b1, b2, x0: 0.5, dt: 0.01, duration: 200 });
  const da = td.decayAnalysis({ t: r.t, x: r.x, inertia: m });
  close(da.b2, b2, 0.05, 'quadratic'); close(da.b1, b1, 0.05, 'linear');
  assert.ok(da.r2 > 0.99);
});

test('Forced response: steady-state amplitude F/|C - w²M + i w B| and nonlinear force callback', () => {
  const m = 1e5; const c = 4e5; const b = 6e4; const F = 1e4; const w = 1.5;
  const r = td.simulateCummins({ M: [[m]], Blin: [[b]], C: [[c]], dt: 0.01, nSteps: 12000, force: (t) => [F * Math.cos(w * t)] });
  let amp = 0; for (let i = 10000; i <= 12000; i++) amp = Math.max(amp, Math.abs(r.x[0][i]));
  close(amp, F / Math.hypot(c - w * w * m, w * b), 2e-3);
  // two-DOF system with a nonlinear (cubic) spring through the callback conserves energy when undamped
  const E = (x, v) => 0.5 * (v[0] ** 2 + 2 * v[1] ** 2) + 0.5 * (x[0] ** 2 + (x[1] - x[0]) ** 2) + 0.25 * x[1] ** 4;
  const s = td.simulateCummins({ M: [[1, 0], [0, 2]], C: [[2, -1], [-1, 1]], dt: 0.002, nSteps: 20000, x0: [1, 0], force: (t, x) => [0, -(x[1] ** 3)] });
  const e0 = E([1, 0], [0, 0]); const e1 = E([s.x[0][20000], s.x[1][20000]], [s.v[0][20000], s.v[1][20000]]);
  close(e1, e0, 1e-4);
});

test('Wave components carry the spectrum variance; excitation superposition returns Re{a X exp(i(wt + phase))}', () => {
  const omega = Array.from({ length: 200 }, (_, i) => 0.2 + i * 0.01);
  const S = omega.map((w) => Math.exp(-((w - 0.8) ** 2) / 0.02));
  const comp = td.waveComponents({ omega, S, seed: 7 });
  let m0 = 0; for (let i = 1; i < omega.length; i++) m0 += 0.5 * (S[i] + S[i - 1]) * 0.01;
  let v = 0; for (const ai of comp.amp) v += 0.5 * ai * ai;
  close(v, m0, 1e-3);
  assert.deepEqual(Array.from(td.waveComponents({ omega, S, seed: 7 }).phase), Array.from(comp.phase));
  const f = td.waveExcitationForce({ components: { omega: [0.5], amp: [2], phase: [0.3] }, Fre: [[3]], Fim: [[4]] });
  close(f(1.7)[0], 2 * (3 * Math.cos(0.5 * 1.7 + 0.3) - 4 * Math.sin(0.5 * 1.7 + 0.3)), 1e-12);
});
