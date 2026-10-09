import test from 'node:test';
import assert from 'node:assert/strict';
import * as dp from '../../src/engine/hydro/dp.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);

const fleet = () => [
  { name: 'bow tunnel', x: 40, y: 0, type: 'tunnel', Tmax: 2.0e5 },
  { name: 'fwd azi', x: 30, y: 0, type: 'azimuth', Tmax: 3.0e5 },
  { name: 'aft azi P', x: -40, y: 8, type: 'azimuth', Tmax: 4.0e5 },
  { name: 'aft azi S', x: -40, y: -8, type: 'azimuth', Tmax: 4.0e5 },
];
const sumForces = (ths, res) => {
  let Fx = 0; let Fy = 0; let Mz = 0;
  res.thrusters.forEach((t, i) => { Fx += t.Tx; Fy += t.Ty; Mz += ths[i].x * t.Ty - ths[i].y * t.Tx; });
  return { Fx, Fy, Mz };
};

test('Thrust allocation satisfies the demand within limits (pseudo-inverse and bounded)', () => {
  const ths = fleet();
  const demand = { Fx: 1.5e5, Fy: 2.0e5, Mz: 1.0e6 };
  for (const method of ['pinv', 'bounded']) {
    const r = dp.allocateThrust({ thrusters: ths, demand, method });
    const s = sumForces(ths, r);
    close(s.Fx, demand.Fx, 1e-9); close(s.Fy, demand.Fy, 1e-9); close(s.Mz, demand.Mz, 1e-9);
    assert.ok(r.feasible && r.maxUtilisation <= 1);
    near(r.thrusters[0].Tx, 0, 1e-6);    // a tunnel thruster only pushes sideways
  }
});

test('Pseudo-inverse is the minimum weighted-norm solution: equal thrusters share a pure surge demand equally', () => {
  const ths = [{ x: -40, y: 8, Tmax: 4e5 }, { x: -40, y: -8, Tmax: 4e5 }];
  const r = dp.allocateThrust({ thrusters: ths, demand: { Fx: 3e5, Fy: 0, Mz: 0 }, method: 'pinv' });
  close(r.thrusters[0].Tx, 1.5e5, 1e-9); close(r.thrusters[1].Tx, 1.5e5, 1e-9);
  near(r.thrusters[0].Ty, 0, 1e-6);
  // any other exact solution has a larger sum of squared utilisations
  const base = r.thrusters.reduce((s, t) => s + t.utilisation ** 2, 0);
  assert.ok(base < (2e5 / 4e5) ** 2 + (1e5 / 4e5) ** 2);
});

test('Bounded allocation redistributes a saturated thruster and still meets the demand', () => {
  const ths = fleet();
  const demand = { Fx: 0, Fy: 9.0e5, Mz: 4.0e6 };
  const free = dp.allocateThrust({ thrusters: ths, demand, method: 'pinv' });
  assert.ok(free.maxUtilisation > 1, 'unconstrained solution should exceed a limit for this test');
  const r = dp.allocateThrust({ thrusters: ths, demand });
  const s = sumForces(ths, r);
  close(s.Fy, demand.Fy, 1e-8); close(s.Mz, demand.Mz, 1e-8); near(s.Fx, 0, 1e-3);
  assert.ok(r.feasible);
  for (const t of r.thrusters) assert.ok(t.utilisation <= 1 + 1e-9, `${t.name} ${t.utilisation}`);
  // pure sway beyond the capability. Hand optimum: bow units at full sway thrust give +17 MN·m of yaw moment; the
  // aft pair (each |T| = 0.4 MN, sway t, opposite surge components ±x at y = ±8 m) must cancel it:
  //   17e6 - 80 t + 16 x = 0 with x² + t² = (4e5)²  ->  6656 t² - 2.72e9 t + 2.4804e14 = 0  ->  t = 2.7129e5 N
  const tAft = (2.72e9 + Math.sqrt(2.72e9 ** 2 - 4 * 6656 * 2.4804e14)) / (2 * 6656);
  const over = dp.allocateThrust({ thrusters: ths, demand: { Fx: 0, Fy: 2e6, Mz: 0 } });
  assert.ok(!over.feasible);
  for (const t of over.thrusters) assert.ok(t.utilisation <= 1 + 1e-9);
  close(over.achieved.Fy, 5e5 + 2 * tAft, 2e-3); close(over.demandFraction, (5e5 + 2 * tAft) / 2e6, 2e-3);
  near(over.achieved.Fx, 0, 1); near(over.achieved.Mz, 0, 10);
});

test('Forbidden azimuth sectors are respected and failed thrusters carry no load', () => {
  const ths = fleet();
  ths[2].forbidden = [[-20, 20]];
  ths[3].forbidden = [[-20, 20]];
  const demand = { Fx: 3.0e5, Fy: 0, Mz: 0 };
  const r = dp.allocateThrust({ thrusters: ths, demand });
  const s = sumForces(ths, r);
  close(s.Fx, demand.Fx, 1e-8); near(s.Fy, 0, 1e-2); near(s.Mz, 0, 1);
  assert.ok(r.feasible);
  for (const i of [2, 3]) {
    const t = r.thrusters[i];
    if (t.T > 1) { const a = ((t.angleDeg + 180) % 360) - 180; assert.ok(Math.abs(a) >= 20 - 1e-6, `${t.name} at ${a} deg`); }
  }
  const f = dp.allocateThrust({ thrusters: fleet(), demand: { Fx: 1e5, Fy: 1e5, Mz: 0 }, disabled: [1] });
  assert.equal(f.thrusters[1].T, 0); assert.ok(f.feasible);
  close(f.achieved.Fy, 1e5, 1e-8);
});

test('DP capability: symmetric azimuth layout balances F = c V² up to V = sqrt(sum Tmax / c); failures reduce it', () => {
  const ths = [[30, 10], [30, -10], [-30, 10], [-30, -10]].map(([x, y]) => ({ x, y, Tmax: 2.5e5 }));
  const c = 800;
  const env = (V, dir) => ({ Fx: c * V * V * Math.cos((dir * Math.PI) / 180), Fy: c * V * V * Math.sin((dir * Math.PI) / 180), Mz: 0 });
  const cap = dp.dpCapability({ thrusters: ths, envForce: env, directionsDeg: [0, 45, 90, 180], vMax: 60 });
  for (const v of cap.intact) close(v, Math.sqrt(1e6 / c), 2e-3);
  cap.failure.forEach((v, i) => { assert.ok(v < cap.intact[i] && v > 0.5 * cap.intact[i]); assert.equal(cap.failureCase[i].length, 1); });
  // with 80 % utilisation the limit drops by sqrt(0.8)
  const c80 = dp.dpCapability({ thrusters: ths, envForce: env, directionsDeg: [90], failures: 'none', maxUtilisation: 0.8 });
  close(c80.intact[0], Math.sqrt(0.8e6 / c), 2e-3); assert.equal(c80.failure, null);
});

test('LQR by the sign-function Riccati solver: scalar closed form and the double integrator K = [1, sqrt 3]', () => {
  const a = 0.5; const b = 2; const q = 3; const r = 0.5;
  const s = dp.lqr({ A: [[a]], B: [[b]], Q: [[q]], R: [[r]] });
  close(s.P[0][0], (r * (a + Math.sqrt(a * a + (b * b * q) / r))) / (b * b), 1e-10);
  const d = dp.lqr({ A: [[0, 1], [0, 0]], B: [[0], [1]], Q: [[1, 0], [0, 1]], R: [[1]] });
  close(d.K[0][0], 1, 1e-9); close(d.K[0][1], Math.sqrt(3), 1e-9);
  assert.ok(d.residual < 1e-9);
  // DP model: Riccati residual and closed-loop stability of the decoupled surge axis (Kp, Kd > 0)
  const M = [[6e6, 0, 0], [0, 9e6, 0], [0, 0, 4e9]]; const D = [[5e4, 0, 0], [0, 2e5, 0], [0, 0, 8e7]];
  const Q = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j ? (i < 3 ? 1 : 100) : 0)));
  const R = [[1e-10, 0, 0], [0, 1e-10, 0], [0, 0, 1e-14]];
  const gains = dp.dpLqrGains({ M, D, Q, R });
  assert.ok(gains.residual < 1e-6);
  close(gains.Kp[0][0], Math.sqrt(1 / 1e-10), 1e-6);    // Kp = sqrt(q/r) for a damped double integrator
  assert.ok(gains.Kd[0][0] > 0 && gains.Kd[2][2] > 0);
});

test('PID pole placement (Fossen): Kp = m wn², Kd = 2 zeta wn m - d, Ki = wn Kp / 10', () => {
  const g = dp.dpPidGains({ M: [[6e6, 0, 0], [0, 9e6, 0], [0, 0, 4e9]], D: [[5e4, 0, 0], [0, 2e5, 0], [0, 0, 8e7]], omegaN: [0.05, 0.05, 0.06] });
  close(g.Kp[0], 6e6 * 0.0025, 1e-12); close(g.Kd[1], 2 * 0.7 * 0.05 * 9e6 - 2e5, 1e-12); close(g.Ki[2], 0.006 * 4e9 * 0.0036, 1e-12);
});

test('Closed-loop DP: holds station against a constant force with thrust equal to the disturbance; wave filtering', () => {
  const M = [[6e6, 0, 0], [0, 9e6, 0], [0, 0, 4e9]]; const D = [[5e4, 0, 0], [0, 2e5, 0], [0, 0, 8e7]];
  const ths = fleet();
  const gains = dp.dpPidGains({ M, D, omegaN: [0.06, 0.06, 0.08] });
  const run = (extra) => dp.simulateDp({ M, D, thrusters: ths, gains, envForce: () => ({ Fx: -1.0e5, Fy: 1.5e5, Mz: 0 }), duration: 1500, dt: 0.25, statsStart: 900, ...extra });
  const calm = run({});
  // integral action removes the steady offset
  near(calm.footprint.meanX, 0, 0.1); near(calm.footprint.meanY, 0, 0.1);
  assert.ok(calm.footprint.maxRadius < 0.5, `radius ${calm.footprint.maxRadius}`);
  assert.ok(calm.footprint.maxPsiDevDeg < 0.5);
  // mean thrust balances the disturbance: sum(u_i Tmax_i) >= |F|, and no thruster is overloaded
  const total = calm.utilisation.reduce((s, u, i) => s + u.mean * ths[i].Tmax, 0);
  assert.ok(total >= Math.hypot(1.0e5, 1.5e5) * 0.999 && total < 1.5 * Math.hypot(1.0e5, 1.5e5), `total ${total}`);
  for (const u of calm.utilisation) assert.ok(u.max <= 1 + 1e-9);
  assert.equal(calm.t.length, calm.x.length);
  // first-order wave motion (1 m std) and sensor noise in the measurements: the wave filter keeps the vessel
  // from chasing it — low-frequency excursions stay far below the wave-frequency motion
  const waves = run({ waveMotion: { sigma: [1.0, 1.0, 0.01], period: 9 }, noiseStd: [0.1, 0.1, 0.001] });
  assert.ok(waves.footprint.stdX < 0.3 && waves.footprint.stdY < 0.3, `std ${waves.footprint.stdX}, ${waves.footprint.stdY}`);
  // deterministic for a seed
  const again = run({ waveMotion: { sigma: [1.0, 1.0, 0.01], period: 9 }, noiseStd: [0.1, 0.1, 0.001] });
  assert.equal(again.footprint.maxRadius, waves.footprint.maxRadius);
});
