import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCatenary } from '../../src/engine/lines/catenary.js';
import * as ms from '../../src/engine/lines/mooringsystem.js';

const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
const chain = { L: 900, w: 1500, EA: 7e8 };
const spread = (n = 8, extra = {}) => ms.spreadLayout({
  nLines: n, anchorRadius: 880, fairleadRadius: 30, fairleadZ: -10, depth: 150, segments: [chain], MBL: 8e6, ...extra,
});

test('symmetric spread: zero net load at zero offset; stiffness matrix symmetric and equal to the line-sum formula', () => {
  const sys = spread(8);
  const r = ms.lineForces(sys);
  const Tsum = r.lines.reduce((a, l) => a + l.T, 0);
  for (const i of [0, 1, 5]) assert.ok(Math.abs(r.F[i]) < 1e-8 * Tsum * (i === 5 ? 30 : 1));
  assert.ok(Math.abs(r.F[3]) < 1e-8 * Tsum * 30 && Math.abs(r.F[4]) < 1e-8 * Tsum * 30);
  // vertical pull = sum of V
  assert.ok(rel(-r.F[2], r.lines.reduce((a, l) => a + l.V, 0)) < 1e-12);
  const K = ms.stiffnessMatrix(sys);
  const sF = K[0][0], sM = K[5][5];
  for (let i = 0; i < 6; i++) for (let j = 0; j < i; j++) {
    const sc = (i < 3 ? 1 : 30) * (j < 3 ? 1 : 30) * sF;
    assert.ok(Math.abs(K[i][j] - K[j][i]) < 2e-5 * sc, `K[${i}][${j}] vs K[${j}][${i}]`);
  }
  // hand formula: Kxx = sum(k cos^2 + (H/X) sin^2) = (n/2)(dH/dX + H/X) for n >= 3 identical lines
  const one = solveCatenary({ X: 850, Z: 140, ...chain });
  assert.ok(rel(r.lines[0].T, one.T) < 1e-9);
  assert.ok(rel(K[0][0], 4 * (one.stiffness[0][0] + one.H / 850)) < 1e-5);
  assert.ok(rel(K[1][1], K[0][0]) < 1e-6);
  // yaw: fairleads on radius r, lines radial: Kyaw = n (H/X) r (r + X)
  assert.ok(rel(sM, 8 * (one.H / 850) * 30 * 880) < 1e-4);
  assert.ok(rel(K[2][2], 8 * one.stiffness[1][1]) < 1e-5);
});

test('restoring curve is monotonic, odd for a symmetric layout, and consistent with the equilibrium solver', () => {
  const sys = spread(8);
  const offsets = [];
  for (let d = -20; d <= 20; d += 2) offsets.push(d);
  const c = ms.restoringCurve(sys, { direction: 0.3, offsets });
  for (let i = 1; i < offsets.length; i++) assert.ok(c.restoring[i] > c.restoring[i - 1]);
  for (let i = 0; i < offsets.length; i++) {
    assert.ok(Math.abs(c.restoring[i] + c.restoring[offsets.length - 1 - i]) < 1e-6 * Math.abs(c.restoring[0]));
  }
  // equilibrium under the restoring force found at 12 m offset returns the 12 m offset
  const k = offsets.indexOf(12);
  const F = c.restoring[k];
  const eq = ms.equilibrium(sys, { Fext: [F * Math.cos(0.3), F * Math.sin(0.3), 0, 0, 0, 0] });
  assert.ok(Math.abs(eq.offset - 12) < 1e-4 && Math.abs(eq.offsetDirection - 0.3) < 1e-5);
  // sum of forces = 0 on the active DOFs
  assert.ok(Math.abs(eq.F[0] + F * Math.cos(0.3)) < 1e-6 * F && Math.abs(eq.F[1] + F * Math.sin(0.3)) < 1e-6 * F);
  assert.ok(Math.abs(eq.F[5]) < 1e-6 * F * 30);
  assert.ok(eq.maxTension > ms.lineForces(sys).maxTension);
});

test('6-DOF equilibrium with hydrostatic stiffness balances all six components', () => {
  const sys = spread(6);
  const hydroK = Array.from({ length: 6 }, () => new Array(6).fill(0));
  hydroK[2][2] = 5e6; hydroK[3][3] = 8e8; hydroK[4][4] = 8e8;
  const Fext = [6e5, 2e5, 0, 0, 5e6, 1e6];
  const eq = ms.equilibrium(sys, { Fext, dofs: [0, 1, 2, 3, 4, 5], hydroK });
  const q = eq.pose;
  for (let i = 0; i < 6; i++) {
    let v = eq.F[i] + Fext[i];
    for (let j = 0; j < 6; j++) v -= hydroK[i][j] * q[j];
    assert.ok(Math.abs(v) < 1e-5 * 6e5 * (i < 3 ? 1 : 30));
  }
  assert.ok(q[2] < 0 && q[0] > 0); // lines pull the floater down
});

test('pretension set-up by line length and by anchor radius', () => {
  const sys = spread(4);
  const a = ms.setPretension(sys, { target: 9e5 });
  for (const t of a.tensions) assert.ok(rel(t, 9e5) < 1e-8);
  for (const l of ms.lineForces(a.system).lines) assert.ok(rel(l.T, 9e5) < 1e-7);
  const b = ms.setPretension(sys, { target: [5e5, 6e5, 7e5, 8e5], mode: 'anchorRadius' });
  ms.lineForces(b.system).lines.forEach((l, i) => assert.ok(rel(l.T, [5e5, 6e5, 7e5, 8e5][i]) < 1e-7));
  assert.equal(sys.lines[0].segments[0].L, 900); // input not modified
  assert.throws(() => ms.setPretension(sys, { target: -1 }), /target/);
});

test('one-line-broken: larger offset and tension; energy-method overshoot factor -> 2 for a linear system', () => {
  const sys = spread(8);
  const res = ms.brokenLineAnalysis(sys, { Fext: [1.5e6, 0, 0, 0, 0, 0], lines: [4, 3] });
  assert.equal(res.damaged.length, 2);
  const d = res.damaged[0]; // line 5 (azimuth 180 deg) is the most loaded for +x load
  assert.equal(res.intact.maxTensionLine, 4);
  assert.ok(d.equilibrium.pose[0] > res.intact.pose[0]);
  assert.ok(d.overshootFactor > 1 && d.overshootFactor <= 2.0001);
  assert.ok(d.transientPose[0] > d.equilibrium.pose[0]);
  assert.ok(d.maxTransientTension > d.equilibrium.maxTension);
  // work-energy check: net work of the unbalanced load from the intact equilibrium to the transient pose is zero
  const q0 = res.intact.pose, q1 = d.transientPose, n = 400;
  let W = 0;
  for (let k = 0; k < n; k++) {
    const s = (k + 0.5) / n, q = q0.map((v, j) => v + s * (q1[j] - v));
    const F = ms.lineForces(sys, q, { exclude: [4] }).F;
    W += ((F[0] + 1.5e6) * (q1[0] - q0[0]) + F[1] * (q1[1] - q0[1]) + F[5] * (q1[5] - q0[5])) / n;
  }
  const scale = 1.5e6 * Math.abs(d.equilibrium.pose[0] - q0[0]);
  assert.ok(Math.abs(W) < 2e-3 * scale);
  // nearly linear taut system (weightless springs, small load): overshoot factor ~ 2
  const taut = ms.spreadLayout({ nLines: 6, anchorRadius: 500, fairleadRadius: 10, depth: 500, segments: [{ L: 690, w: 0, EA: 2e8 }] });
  const lin = ms.brokenLineAnalysis(taut, { Fext: [2e4, 0, 0, 0, 0, 0], lines: [3], dofs: [0, 1] });
  assert.ok(Math.abs(lin.damaged[0].overshootFactor - 2) < 0.02);
});

test('API RP 2SK safety factors and DNV-OS-E301 partial-factor utilisation (hand calculations)', () => {
  // 4000 kN on 8000 kN MBL: SF = 2.0; intact dynamic needs 1.67 -> UC = 0.835
  const a = ms.checkApi2sk({ tension: 4e6, MBL: 8e6 });
  assert.ok(rel(a.safetyFactor, 2) < 1e-12 && rel(a.uc, 0.835) < 1e-12 && a.ok);
  assert.equal(ms.checkApi2sk({ tension: 4e6, MBL: 8e6, analysis: 'quasi-static' }).required, 2.0);
  assert.equal(ms.checkApi2sk({ tension: 4e6, MBL: 8e6, condition: 'damaged' }).required, 1.25);
  assert.equal(ms.checkApi2sk({ tension: 4e6, MBL: 8e6, condition: 'damaged', analysis: 'quasi-static' }).required, 1.43);
  assert.equal(ms.checkApi2sk({ tension: 5e6, MBL: 8e6 }).ok, false);
  // DNV ULS class 2: (2000*1.40 + 1500*2.10)/(0.95*8000) = 5950/7600 = 0.78289
  const d = ms.checkDnvE301({ Tmean: 2e6, Tdyn: 1.5e6, MBL: 8e6, consequenceClass: 2 });
  assert.ok(rel(d.uc, 5950 / 7600) < 1e-12);
  // ALS class 1: (2000*1.00 + 1500*1.10)/7600
  assert.ok(rel(ms.checkDnvE301({ Tmean: 2e6, Tdyn: 1.5e6, MBL: 8e6, limitState: 'ALS' }).uc, 3650 / 7600) < 1e-12);
  assert.deepEqual(ms.DNV_E301_FACTORS.ULS[1].dynamic, [1.10, 1.50]);
  assert.deepEqual(ms.DNV_E301_FACTORS.ALS[2].dynamic, [1.00, 1.25]);
  assert.throws(() => ms.checkDnvE301({ Tmean: 1, Tdyn: 1, MBL: 1, consequenceClass: 3 }), /consequenceClass/);
  assert.equal(ms.checkAnchorUplift({ upliftAngle: 0.02 }).ok, false);
  assert.equal(ms.checkGroundedLength({ groundedLength: 120, required: 100 }).margin, 20);
  assert.ok(rel(ms.checkOffset({ offset: 12, depth: 150, limitFraction: 0.1 }).uc, 0.8) < 1e-12);
});

test('API RP 2SK mean + LF + WF combination and quasi-dynamic tensions', () => {
  // sigma = 2 m, N = 100: sig = 4 m, max = 2 sqrt(2 ln 100) = 6.0697 m
  const st = ms.responseStatistics({ sigma: 2, N: 100 });
  assert.ok(rel(st.max, 6.0697) < 1e-5 && st.sig === 4);
  const c = ms.combineApi2sk({ mean: 10, lfSig: 4, lfMax: 6, wfSig: 2, wfMax: 3.5 });
  assert.equal(c.case1, 18); assert.equal(c.case2, 17.5); assert.equal(c.max, 18);
  const sys = spread(8);
  const eq = ms.equilibrium(sys, { Fext: [1e6, 0, 0, 0, 0, 0] });
  const lf = { sig: 4, max: 6 }, wf = { sig: 2, max: 3.5 };
  const qd = ms.quasiDynamicTension(sys, { meanPose: eq.pose, direction: 0, lf, wf });
  const l5 = qd.lines.find((l) => l.index === 4);
  const at = (d) => ms.lineForces(sys, [eq.pose[0] + d, eq.pose[1], 0, 0, 0, eq.pose[5]]).lines[4].T;
  assert.ok(rel(l5.case1, at(8)) < 1e-9 && rel(l5.case2, at(7.5)) < 1e-9 && rel(l5.Tmean, at(0)) < 1e-9);
  // DAF 1.5 amplifies only the wave-frequency increment
  const qd2 = ms.quasiDynamicTension(sys, { meanPose: eq.pose, direction: 0, lf, wf, daf: 1.5 });
  assert.ok(rel(qd2.lines.find((l) => l.index === 4).case1, at(6) + 1.5 * (at(8) - at(6))) < 1e-9);
});

test('TLP tendons: inextensible set-down d = L - sqrt(L^2 - x^2), F = T0 tan(a); stiffness and tension distribution', () => {
  const tendons = [[30, 30], [-30, 30], [-30, -30], [30, -30]].map(([x, y]) => ({ x, y, L: 400, EA: 1e16, T0: 5e6 }));
  const r = ms.tendonOffset({ tendons, offset: 40 });
  assert.ok(rel(r.setDown, 400 - Math.sqrt(400 * 400 - 40 * 40)) < 1e-6);
  // no waterplane stiffness: T cos(a) = T0 -> F = 4 T0 tan(a)
  assert.ok(rel(r.restoring, 4 * 5e6 * 40 / Math.sqrt(400 * 400 - 40 * 40)) < 1e-6);
  // waterplane stiffness raises the tension with set-down: sum(T cos a) = sum(T0) + rho g Aw d
  const el = tendons.map((t) => ({ ...t, EA: 2e9 }));
  const r2 = ms.tendonOffset({ tendons: el, offset: 40, waterplaneStiffness: 1e6 });
  assert.ok(rel(r2.tensions.reduce((a, T, i) => a + T * Math.cos(r2.angles[i]), 0), 2e7 + 1e6 * r2.setDown) < 1e-9);
  assert.ok(r2.setDown < r.setDown && r2.tensions[0] > 5e6);
  const k = ms.tendonStiffness(el);
  assert.equal(k.heave, 4 * 2e9 / 400);
  assert.equal(k.surge, 4 * 5e6 / 400);
  assert.equal(k.pitch, 4 * (2e9 / 400) * 900);
  assert.equal(k.yaw, 4 * (5e6 / 400) * 1800);
  // pitch moment My = 6e7 N m: dT = -+ My/(4*30) = 5e5 N; positive My lowers x > 0 corners -> tension drops there
  const tt = ms.tendonTensions({ tendons: el, My: 6e7 });
  assert.ok(rel(tt.maxTension, 5.5e6) < 1e-10 && rel(tt.minTension, 4.5e6) < 1e-10);
  assert.ok(rel(tt.tensions[0], 4.5e6) < 1e-10);
  assert.ok(rel(ms.tendonTensions({ tendons: el, Fz: 4e6 }).tensions[2], 6e6) < 1e-10);
  // natural period: T = 2 pi sqrt((M + A)/K)
  assert.ok(rel(ms.naturalPeriod({ mass: 4e7, addedMass: 1e7, stiffness: k.surge }), 2 * Math.PI * Math.sqrt(5e7 / 5e4)) < 1e-12);
});

test('layout generators, turret weathervaning and CALM hawser hand calculation', () => {
  const cl = ms.clusterLayout({ nGroups: 3, linesPerGroup: 3, groupSpread: 5 * Math.PI / 180, anchorRadius: 880, fairleadRadius: 30, depth: 150, segments: [chain] });
  assert.equal(cl.lines.length, 9);
  const F = ms.lineForces(cl).F;
  assert.ok(Math.hypot(F[0], F[1]) < 1e-6 * 1e6);
  assert.ok(rel(Math.atan2(cl.lines[2].anchor[1], cl.lines[2].anchor[0]), 5 * Math.PI / 180) < 1e-12);
  const tur = ms.turretLayout({ nLines: 6, anchorRadius: 900, fairleadRadius: 5, depth: 150, segments: [chain], turretCenter: [80, 0, -15] });
  assert.ok(tur.lines.every((l) => l.turret));
  // yawing the hull about the turret axis leaves the line tensions unchanged
  const yaw = 0.6;
  const p = [80 - 80 * Math.cos(yaw), -80 * Math.sin(yaw), 0, 0, 0, yaw];
  const a = ms.lineForces(tur), b = ms.lineForces(tur, p);
  a.lines.forEach((l, i) => assert.ok(rel(b.lines[i].T, l.T) < 1e-9));
  // CALM hawser: Fmean 500 kN, kh = 200 kN/m, kb = 300 kN/m -> ks = 120 kN/m; M = 1.2e8 kg, v = 0.1 m/s
  const h = ms.calmHawser({ Fmean: 5e5, hawserStiffness: 2e5, buoyStiffness: 3e5, mass: 1e8, addedMass: 2e7, surgeVelocity: 0.1, MBL: 4e6 });
  assert.ok(rel(h.seriesStiffness, 1.2e5) < 1e-12);
  assert.ok(rel(h.dT, 0.1 * Math.sqrt(1.2e8 * 1.2e5)) < 1e-12); // 379.5 kN
  assert.ok(rel(h.Tmax, 5e5 + 379473.3192) < 1e-9);
});

test('clearances: skew segments, polylines, adjacent mooring lines and a pipeline crossing', () => {
  // skew lines along x at z = 0 and along y at z = 3: distance 3
  const d = ms.segmentDistance([[-1, 0, 0], [1, 0, 0]], [[0, -1, 3], [0, 1, 3]]);
  assert.ok(rel(d.distance, 3) < 1e-14 && d.s === 0.5 && d.t === 0.5);
  // end-point governed: distance from (2,0,0) to segment end (1,0,0)... = sqrt(1 + 4)
  assert.ok(rel(ms.segmentDistance([[0, 0, 0], [1, 0, 0]], [[2, 2, 0], [3, 5, 0]]).distance, Math.sqrt(5)) < 1e-14);
  assert.ok(rel(ms.polylineMinDistance([[0, 0, 0], [10, 0, 0], [10, 10, 0]], [[5, 5, 4]]).distance, Math.sqrt(41)) < 1e-14);
  const sys = spread(8);
  const c = ms.clearances(sys, { assets: [{ name: 'pipeline', points: [[400, -500, -150], [400, 500, -150]], radius: 0.5 }] });
  // radial lines are closest at the fairleads: 2 r sin(pi/n)
  assert.ok(rel(c.minLineLine, 2 * 30 * Math.sin(Math.PI / 8)) < 1e-9);
  // line 1 lies along +x: it crosses the pipeline at x = 400, where it is grounded or just above the seabed
  const first = c.lineAsset.find((a) => a.line === 0);
  const sol = ms.lineForces(sys).lines[0];
  const xFromAnchor = 880 - 400;
  if (xFromAnchor < sol.solution.touchdownX) assert.ok(first.distance < 1e-6);
  assert.ok(first.clearance === first.distance - 0.5);
});
