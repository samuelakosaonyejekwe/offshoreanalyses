import test from 'node:test';
import assert from 'node:assert/strict';
import { rk4, rk45, implicitEuler, bdf2, newmark, newmarkNonlinear } from '../../src/engine/math/ode.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

const sho = (t, y) => [y[1], -y[0]]; // unit harmonic oscillator, y = (cos t, −sin t) from (1, 0)

test('rk4: exponential decay y′ = −y, fourth-order convergence (error ratio 16 per halving)', () => {
  const errs = [10, 20, 40, 80].map((n) => {
    const r = rk4((t, y) => [-y[0]], 0, [1], 2, n);
    assert.equal(r.t.length, n + 1);
    assert.equal(r.t[n], 2);
    return Math.abs(r.y[n][0] - Math.exp(-2));
  });
  for (let k = 1; k < errs.length; k++) {
    const order = Math.log2(errs[k - 1] / errs[k]);
    assert.ok(order > 3.8 && order < 4.2, `observed order ${order}`);
  }
});

test('rk4: harmonic oscillator over 10 periods vs cos/sin; onStep callback and early stop', () => {
  const T = 20 * Math.PI;
  const r = rk4(sho, 0, [1, 0], T, 4000);
  for (const k of [1000, 2500, 4000]) {
    // RK4 phase error of the oscillator: h⁴ t / 120 = 3.2e-8 at t = 20π with h = π/200
    close(r.y[k][0], Math.cos(r.t[k]), 5e-8);
    close(r.y[k][1], -Math.sin(r.t[k]), 5e-8);
  }
  assert.ok(r.y[0] instanceof Float64Array);
  let calls = 0;
  const s = rk4(sho, 0, [1, 0], 1, 100, (t, y, k) => {
    calls++;
    return k < 30;
  });
  assert.equal(calls, 30);
  assert.equal(s.t.length, 31);
  close(s.t[30], 0.3, 1e-12);
  assert.throws(() => rk4(sho, 0, [1, 0], 1, 0), RangeError);
});

test('rk4: right-hand side returning a reused buffer is handled', () => {
  const buf = new Float64Array(2);
  const f = (t, y) => {
    buf[0] = y[1];
    buf[1] = -y[0];
    return buf;
  };
  const r = rk4(f, 0, [1, 0], 3, 300);
  close(r.y[300][0], Math.cos(3), 1e-9);
});

test('rk45: harmonic oscillator over 10 periods, global error consistent with the tolerance', () => {
  const T = 20 * Math.PI;
  const r = rk45(sho, 0, [1, 0], T, { rtol: 1e-9, atol: 1e-12 });
  const n = r.t.length - 1;
  assert.equal(r.t[n], T);
  assert.equal(r.steps, n);
  let emax = 0;
  for (let k = 0; k <= n; k++) {
    emax = Math.max(emax, Math.abs(r.y[k][0] - Math.cos(r.t[k])), Math.abs(r.y[k][1] + Math.sin(r.t[k])));
  }
  assert.ok(emax < 2e-7, `max error ${emax}`);
  // energy drift of the non-symplectic scheme stays at tolerance level
  const yn = r.y[n];
  close(yn[0] * yn[0] + yn[1] * yn[1], 1, 4e-7);
  assert.ok(r.steps < 3000, `efficient: ${r.steps} steps`);
  assert.ok(r.nfev <= 6 * (r.steps + r.rejected) + 3, 'FSAL: six evaluations per step');
});

test('rk45: error decreases in proportion to the tolerance (tolerance proportionality)', () => {
  const errs = [1e-4, 1e-6, 1e-8, 1e-10].map((rtol) => {
    const r = rk45(sho, 0, [1, 0], 10, { rtol, atol: rtol * 1e-3 });
    return Math.abs(r.y[r.y.length - 1][0] - Math.cos(10));
  });
  for (let k = 1; k < errs.length; k++) assert.ok(errs[k] < errs[k - 1] / 10, `errors ${errs}`);
  assert.ok(errs[3] < 1e-8);
});

test('rk45: dense output at prescribed times matches the analytic solution to the integration accuracy', () => {
  const times = [];
  for (let i = 0; i <= 400; i++) times.push(i * 0.05);
  const r = rk45(sho, 0, [1, 0], 20, { rtol: 1e-9, atol: 1e-12, dense: times });
  assert.equal(r.t.length, 401);
  assert.equal(r.steps, rk45(sho, 0, [1, 0], 20, { rtol: 1e-9, atol: 1e-12 }).steps, 'dense output does not alter the step sequence');
  let e = 0;
  for (let k = 0; k <= 400; k++) {
    assert.equal(r.t[k], times[k]);
    e = Math.max(e, Math.abs(r.y[k][0] - Math.cos(times[k])), Math.abs(r.y[k][1] + Math.sin(times[k])));
  }
  assert.ok(e < 1e-7, `dense error ${e}`);
  // with loose tolerance (large steps) the 4th-order interpolant must still follow the solution
  const c = rk45(sho, 0, [1, 0], 20, { rtol: 1e-5, atol: 1e-8, dense: times });
  let ec = 0;
  for (let k = 0; k <= 400; k++) ec = Math.max(ec, Math.abs(c.y[k][0] - Math.cos(times[k])));
  assert.ok(ec < 5e-4, `coarse dense error ${ec}`);
  assert.throws(() => rk45(sho, 0, [1, 0], 20, { dense: [1, 0.5] }), RangeError);
  assert.throws(() => rk45(sho, 0, [1, 0], 20, { dense: [21] }), RangeError);
});

test('rk45: stiff-ish relaxation y′ = −50 (y − cos t), analytic y = (2500 cos t + 50 sin t − 2500 e^(−50t))/2501', () => {
  const exact = (t) => (2500 * Math.cos(t) + 50 * Math.sin(t) - 2500 * Math.exp(-50 * t)) / 2501;
  const r = rk45((t, y) => [-50 * (y[0] - Math.cos(t))], 0, [0], 3, { rtol: 1e-8, atol: 1e-10 });
  let e = 0;
  for (let k = 0; k < r.t.length; k++) e = Math.max(e, Math.abs(r.y[k][0] - exact(r.t[k])));
  assert.ok(e < 5e-8, `max error ${e}`);
  // once the transient has died the step is limited by stability (h·50 <~ 3.3), not accuracy
  assert.ok(r.steps > 3 * 50 / 3.4, 'stability-limited step count');
  assert.ok(r.steps < 400, `steps ${r.steps}`);
});

test('rk45: backward integration, hMax, zero span, scalar-like systems and step-limit error', () => {
  const b = rk45(sho, 5, [Math.cos(5), -Math.sin(5)], 0, { rtol: 1e-9, atol: 1e-12 });
  close(b.y[b.y.length - 1][0], 1, 1e-7);
  close(b.y[b.y.length - 1][1], 0, 1e-7);
  assert.ok(b.t[1] < b.t[0]);
  const h = rk45((t, y) => [1], 0, [0], 10, { hMax: 0.5 });
  assert.ok(h.steps >= 20);
  for (let k = 1; k < h.t.length; k++) assert.ok(h.t[k] - h.t[k - 1] <= 0.5 + 1e-12);
  close(h.y[h.y.length - 1][0], 10, 1e-12);
  const z = rk45(sho, 1, [1, 0], 1);
  assert.equal(z.steps, 0);
  assert.deepEqual(Array.from(z.y[0]), [1, 0]);
  assert.throws(() => rk45(sho, 0, [1, 0], 1000, { rtol: 1e-10, atol: 1e-12, maxSteps: 50 }), /maxSteps/);
  assert.throws(() => rk45(sho, 0, [1, 0], 1, { rtol: 0, atol: 0 }), RangeError);
});

test('rk45: Kepler orbit with eccentricity 0.6 — energy and angular momentum conserved over one period', () => {
  const e = 0.6;
  const f = (t, y) => {
    const r3 = Math.hypot(y[0], y[1]) ** 3;
    return [y[2], y[3], -y[0] / r3, -y[1] / r3];
  };
  const y0 = [1 - e, 0, 0, Math.sqrt((1 + e) / (1 - e))];
  const r = rk45(f, 0, y0, 2 * Math.PI, { rtol: 1e-10, atol: 1e-12 });
  const y = r.y[r.y.length - 1];
  for (let i = 0; i < 4; i++) close(y[i], y0[i], 1e-7, 'periodic orbit closes');
  const energy = 0.5 * (y[2] * y[2] + y[3] * y[3]) - 1 / Math.hypot(y[0], y[1]);
  close(energy, -0.5, 1e-8);
  close(y[0] * y[3] - y[1] * y[2], Math.sqrt(1 - e * e), 1e-8);
  assert.ok(r.rejected >= 0 && r.steps < 1500);
});

test('implicitEuler (order 1) and bdf2 (order 2) on the stiff Prothero–Robinson problem y′ = −1000 (y − cos t) − sin t', () => {
  // exact solution y = cos t for y(0) = 1; explicit methods would need h < 0.002
  const f = (t, y) => [-1000 * (y[0] - Math.cos(t)) - Math.sin(t)];
  const ie = [20, 40, 80].map((n) => Math.abs(implicitEuler(f, 0, [1], 2, n).y[n][0] - Math.cos(2)));
  const bd = [20, 40, 80].map((n) => Math.abs(bdf2(f, 0, [1], 2, n).y[n][0] - Math.cos(2)));
  for (let k = 1; k < 3; k++) {
    const o1 = Math.log2(ie[k - 1] / ie[k]);
    const o2 = Math.log2(bd[k - 1] / bd[k]);
    assert.ok(o1 > 0.8 && o1 < 1.3, `implicit Euler order ${o1}`);
    assert.ok(o2 > 1.7 && o2 < 2.4, `BDF2 order ${o2}`);
  }
  assert.ok(bd[2] < ie[2] / 20);
  assert.ok(ie[0] < 1e-3, 'stable with h = 0.1 although h·λ = −100');
});

test('bdf2: linear stiff system with eigenvalues −1 and −1000, analytic solution, analytic vs finite-difference Jacobian', () => {
  // y1′ = −y1, y2′ = −1000 (y2 − y1): y1 = e^(−t), y2 = (1000/999) e^(−t) + C e^(−1000t)
  const f = (t, y) => [-y[0], -1000 * (y[1] - y[0])];
  const jac = () => [[-1, 0], [1000, -1000]];
  const n = 200;
  const a = bdf2(f, 0, [1, 0], 4, n, { jac });
  const b = bdf2(f, 0, [1, 0], 4, n);
  close(a.y[n][0], Math.exp(-4), 2e-5);
  close(a.y[n][1], (1000 / 999) * Math.exp(-4), 2e-5);
  close(b.y[n][0], a.y[n][0], 1e-9);
  close(b.y[n][1], a.y[n][1], 1e-9);
  const ie = implicitEuler(f, 0, [1, 0], 4, n, { jac });
  close(ie.y[n][0], Math.exp(-4), 2e-3);
  assert.equal(a.t.length, n + 1);
  assert.equal(a.t[n], 4);
});

test('bdf2 / implicitEuler: nonlinear logistic equation y′ = y(1 − y) (Newton inside), second-order accuracy', () => {
  const exact = (t) => 1 / (1 + 9 * Math.exp(-t));
  const f = (t, y) => [y[0] * (1 - y[0])];
  const e1 = Math.abs(bdf2(f, 0, [0.1], 5, 50).y[50][0] - exact(5));
  const e2 = Math.abs(bdf2(f, 0, [0.1], 5, 100).y[100][0] - exact(5));
  const order = Math.log2(e1 / e2);
  assert.ok(order > 1.7 && order < 2.4, `order ${order}`);
  assert.ok(e2 < 1e-4);
  close(implicitEuler(f, 0, [0.1], 5, 2000).y[2000][0], exact(5), 5e-4);
});

// ---------------------------------------------------------------------------------------------------------------

test('newmark (average acceleration), undamped SDOF: energy is conserved to round-off over 100 periods', () => {
  const m = 3;
  const k = 3 * (2 * Math.PI) ** 2; // T = 1 s
  const r = newmark({ M: m, K: k, u0: 0.2, v0: -0.5, dt: 0.02, nSteps: 5000 });
  const E0 = 0.5 * m * 0.25 + 0.5 * k * 0.04;
  let dev = 0;
  for (let i = 0; i <= 5000; i++) {
    const E = 0.5 * m * r.v[i][0] ** 2 + 0.5 * k * r.u[i][0] ** 2;
    dev = Math.max(dev, Math.abs(E - E0) / E0);
  }
  assert.ok(dev < 1e-11, `relative energy deviation ${dev}`);
  assert.equal(r.t.length, 5001);
  close(r.t[5000], 100, 1e-9);
  // initial acceleration from equilibrium: a0 = −k u0 / m
  close(r.a[0][0], (-k * 0.2) / m, 1e-12);
});

test('newmark (average acceleration): period elongation — discrete solution is cos(ω̄ t) with tan(ω̄Δt/2) = ωΔt/2, ΔT/T ≈ (ωΔt)²/12', () => {
  const w = 2 * Math.PI;
  for (const dt of [0.05, 0.1]) {
    const n = Math.round(20 / dt);
    const r = newmark({ M: 1, K: w * w, u0: 1, v0: 0, dt, nSteps: n });
    const wbar = (2 / dt) * Math.atan((w * dt) / 2);
    for (const i of [1, 7, n >> 1, n]) close(r.u[i][0], Math.cos(wbar * i * dt), 1e-10, `dt=${dt}, step ${i}`);
    // period measured from the last upward zero crossing pair
    const zc = [];
    for (let i = 1; i <= n; i++) {
      const a = r.u[i - 1][0];
      const b = r.u[i][0];
      if (a < 0 && b >= 0) zc.push((i - 1 + a / (a - b)) * dt);
    }
    const Tnum = (zc[zc.length - 1] - zc[0]) / (zc.length - 1);
    const elong = Tnum - 1;
    close(elong, (2 * Math.PI) / wbar - 1, 2e-4, 'measured vs exact discrete period');
    close(elong / ((w * dt) ** 2 / 12), 1, 0.08, 'leading-order elongation (ωΔt)²/12');
  }
});

test('newmark: damped SDOF step response vs closed form; second-order convergence in Δt', () => {
  const m = 2;
  const zeta = 0.05;
  const w = 4;
  const k = m * w * w;
  const c = 2 * zeta * w * m;
  const F0 = 10;
  const wd = w * Math.sqrt(1 - zeta * zeta);
  const exact = (t) => (F0 / k) * (1 - Math.exp(-zeta * w * t) * (Math.cos(wd * t) + (zeta / Math.sqrt(1 - zeta * zeta)) * Math.sin(wd * t)));
  const errs = [0.04, 0.02, 0.01].map((dt) => {
    const n = Math.round(3 / dt);
    const r = newmark({ M: [[m]], C: [[c]], K: [[k]], f: () => [F0], u0: [0], v0: [0], dt, nSteps: n });
    return Math.abs(r.u[n][0] - exact(3));
  });
  for (let i = 1; i < 3; i++) {
    const order = Math.log2(errs[i - 1] / errs[i]);
    assert.ok(order > 1.8 && order < 2.2, `order ${order}`);
  }
  assert.ok(errs[2] < 2e-4);
  // long-time limit is the static deflection F0/k
  const s = newmark({ M: m, C: 20 * c, K: k, f: () => F0, u0: 0, v0: 0, dt: 0.05, nSteps: 2000 });
  close(s.u[2000][0], F0 / k, 1e-9);
});

test('newmark: 2-DOF chain released in its first mode shape vibrates harmonically at ω₁² = (k/m)(3 − √5)/2; energy conserved', () => {
  const k = 500;
  const m = 2;
  const K = [[2 * k, -k], [-k, k]];
  const M = [[m, 0], [0, m]];
  const w1 = Math.sqrt(((k / m) * (3 - Math.sqrt(5))) / 2);
  const phi = [1, (1 + Math.sqrt(5)) / 2]; // (2 − λ m/k) from the first row
  const dt = 0.0005;
  const n = 4000;
  const r = newmark({ M, K, u0: phi, v0: [0, 0], dt, nSteps: n });
  for (const i of [500, 2000, 4000]) {
    const c = Math.cos(w1 * i * dt);
    // phase error ≈ ω₁ t (ω₁Δt)²/12 ≈ 4e-5 at t = 2 s
    close(r.u[i][0], phi[0] * c, 1e-4);
    close(r.u[i][1], phi[1] * c, 1e-4);
  }
  const energy = (i) => {
    const u = r.u[i];
    const v = r.v[i];
    return 0.5 * m * (v[0] * v[0] + v[1] * v[1]) + 0.5 * (2 * k * u[0] * u[0] - 2 * k * u[0] * u[1] + k * u[1] * u[1]);
  };
  close(energy(n) / energy(0), 1, 1e-11);
  // second mode released likewise
  const w2 = Math.sqrt(((k / m) * (3 + Math.sqrt(5))) / 2);
  const phi2 = [1, (1 - Math.sqrt(5)) / 2];
  const r2 = newmark({ M, K, u0: phi2, v0: [0, 0], dt, nSteps: n });
  // the scheme decouples modally, so the discrete response is exactly cos(ω̄₂ t) with tan(ω̄Δt/2) = ωΔt/2
  const w2bar = (2 / dt) * Math.atan((w2 * dt) / 2);
  close(r2.u[n][0], Math.cos(w2bar * n * dt), 1e-8);
  close(r2.u[n][1], phi2[1] * Math.cos(w2bar * n * dt), 1e-8);
  close(r2.u[n][0], Math.cos(w2 * n * dt), 2e-3);
});

test('newmark: resonant harmonic forcing of an undamped SDOF grows linearly, u = (F0/2mω) t sin ωt', () => {
  const w = 5;
  const F0 = 3;
  const dt = 0.001;
  const n = 6000;
  const r = newmark({ M: 1, K: w * w, f: (t) => F0 * Math.sin(w * t), u0: 0, v0: 0, dt, nSteps: n });
  // u = F0/(2ω²) (sin ωt − ωt cos ωt)
  const exact = (t) => (F0 / (2 * w * w)) * (Math.sin(w * t) - w * t * Math.cos(w * t));
  for (const i of [1500, 4000, 6000]) close(r.u[i][0], exact(i * dt), 2e-4);
});

test('newmark HHT-α: α = 0 reproduces Newmark; α < 0 damps high frequencies yet stays second-order accurate', () => {
  const w = 2 * Math.PI;
  const base = newmark({ M: 1, K: w * w, u0: 1, v0: 0, dt: 0.01, nSteps: 100 });
  const hht0 = newmark({ M: 1, K: w * w, u0: 1, v0: 0, dt: 0.01, nSteps: 100, alphaHHT: 0 });
  assert.deepEqual(hht0.u, base.u);
  // high-frequency response (ωΔt ≈ 31): average acceleration keeps the amplitude, HHT-α dissipates it
  const big = { M: 1, K: w * w, u0: 1, v0: 0, dt: 5, nSteps: 60 };
  const undamped = newmark(big);
  const damped = newmark({ ...big, alphaHHT: -0.3 });
  const amp = (r) => Math.max(...r.u.slice(40).map((u) => Math.abs(u[0])));
  assert.ok(amp(undamped) > 0.9, `average acceleration amplitude ${amp(undamped)}`);
  assert.ok(amp(damped) < 0.05, `HHT amplitude ${amp(damped)}`);
  // low-frequency accuracy: error ratio ≈ 4 when Δt is halved
  // sampled away from a displacement peak so that the (second-order) phase error is visible
  const errs = [0.02, 0.01].map((dt) => {
    const n = Math.round(2.3 / dt);
    const r = newmark({ M: 1, K: w * w, u0: 1, v0: 0, dt, nSteps: n, alphaHHT: -0.1 });
    return Math.abs(r.u[n][0] - Math.cos(w * 2.3));
  });
  const order = Math.log2(errs[0] / errs[1]);
  assert.ok(order > 1.7 && order < 2.3, `HHT order ${order}`);
  assert.throws(() => newmark({ M: 1, K: 1, u0: 1, v0: 0, dt: 0.1, nSteps: 2, alphaHHT: 0.1 }), RangeError);
  assert.throws(() => newmark({ M: 1, K: 1, u0: 1, v0: 0, dt: 0.1, nSteps: 2, alphaHHT: -0.5 }), RangeError);
});

test('newmark: linear-acceleration variant (β = 1/6) is more accurate than average acceleration; argument checks', () => {
  const w = 2 * Math.PI;
  const run = (beta) => {
    const r = newmark({ M: 1, K: w * w, u0: 1, v0: 0, dt: 0.02, nSteps: 500, beta, gamma: 0.5 });
    return Math.abs(r.u[500][0] - Math.cos(w * 10));
  };
  assert.ok(run(1 / 6) < run(0.25) / 1.5);
  assert.throws(() => newmark({ M: 0, K: 1, u0: 1, v0: 0, dt: 0.1, nSteps: 2 }), RangeError);
  assert.throws(() => newmark({ M: 1, K: 1, u0: 1, v0: 0, dt: 0, nSteps: 2 }), RangeError);
  assert.throws(() => newmark({ M: [[1, 0], [0, 1]], K: [[1]], u0: [0, 0], v0: [0, 0], dt: 0.1, nSteps: 2 }), RangeError);
  assert.throws(() => newmark({ M: 1, K: 1, f: () => [1, 2], u0: 1, v0: 0, dt: 0.1, nSteps: 2 }), RangeError);
});

test('newmarkNonlinear: with a linear internal force it reproduces the linear Newmark solution', () => {
  const K = [[30, -10], [-10, 10]];
  const M = [[2, 0], [0, 1]];
  const C = [[0.3, -0.1], [-0.1, 0.1]];
  const load = (t) => [0, Math.sin(3 * t)];
  const lin = newmark({ M, C, K, f: load, u0: [0.1, 0], v0: [0, 0.2], dt: 0.01, nSteps: 300 });
  const nl = newmarkNonlinear({
    M, C,
    fint: (u) => ({ f: [K[0][0] * u[0] + K[0][1] * u[1], K[1][0] * u[0] + K[1][1] * u[1]], Kt: K }),
    fext: load, u0: [0.1, 0], v0: [0, 0.2], dt: 0.01, nSteps: 300,
  });
  for (const i of [1, 100, 300]) {
    for (let d = 0; d < 2; d++) {
      close(nl.u[i][d], lin.u[i][d], 1e-10);
      close(nl.v[i][d], lin.v[i][d], 1e-9);
      close(nl.a[i][d], lin.a[i][d], 1e-7);
    }
  }
  // a linear problem converges in one Newton iteration per step
  assert.equal(nl.iterations, 300);
  assert.equal(nl.iterationsPerStep.length, 300);
});

test('newmarkNonlinear: hardening Duffing oscillator ü + u + u³ = 0 vs an RK45 reference; energy nearly conserved', () => {
  const dt = 0.002;
  const n = 5000;
  const r = newmarkNonlinear({
    M: 1,
    fint: (u) => ({ f: u[0] + u[0] ** 3, Kt: 1 + 3 * u[0] * u[0] }),
    u0: 1.5, v0: 0, dt, nSteps: n,
  });
  const ref = rk45((t, y) => [y[1], -y[0] - y[0] ** 3], 0, [1.5, 0], n * dt, { rtol: 1e-11, atol: 1e-13 });
  const yr = ref.y[ref.y.length - 1];
  close(r.u[n][0], yr[0], 2e-4);
  close(r.v[n][0], yr[1], 4e-4);
  const E = (u, v) => 0.5 * v * v + 0.5 * u * u + 0.25 * u ** 4;
  const E0 = E(1.5, 0);
  let dev = 0;
  for (let i = 0; i <= n; i++) dev = Math.max(dev, Math.abs(E(r.u[i][0], r.v[i][0]) - E0) / E0);
  assert.ok(dev < 1e-5, `energy deviation ${dev}`);
  assert.ok(r.iterations >= n);
});

test('newmarkNonlinear: velocity-dependent internal force (quadratic drag) with tangent Ct — terminal velocity sqrt(F/c)', () => {
  // m v̇ + c v|v| = F, no stiffness: v -> sqrt(F/c) = 2, exact v(t) = 2 tanh(t·sqrt(F c)/m)
  const c = 2.5;
  const F = 10;
  const r = newmarkNonlinear({
    M: 4,
    fint: (u, v) => ({ f: c * v[0] * Math.abs(v[0]), Kt: 0, Ct: 2 * c * Math.abs(v[0]) }),
    fext: () => F, u0: 0, v0: 0, dt: 0.01, nSteps: 1000,
  });
  close(r.v[1000][0], 2, 1e-6);
  close(r.v[100][0], 2 * Math.tanh((1 * Math.sqrt(F * c)) / 4), 1e-4);
  close(r.a[1000][0], 0, 1e-5);
});

test('newmarkNonlinear: non-convergence and invalid arguments are reported', () => {
  const duffing = { M: 1, fint: (u) => ({ f: u[0] + 50 * u[0] ** 3, Kt: 1 + 150 * u[0] * u[0] }), u0: 2, v0: 0, dt: 0.2, nSteps: 10 };
  assert.throws(() => newmarkNonlinear({ ...duffing, tol: 1e-15, maxIter: 1 }), /did not converge at step 1/);
  assert.doesNotThrow(() => newmarkNonlinear({ ...duffing, tol: 1e-9, maxIter: 30 }));
  assert.throws(() => newmarkNonlinear({ ...duffing, beta: 0 }), RangeError);
  assert.throws(() => newmarkNonlinear({ ...duffing, fint: null }), RangeError);
  assert.throws(() => newmarkNonlinear({ ...duffing, fint: () => ({ f: 1 }) }), RangeError);
});
