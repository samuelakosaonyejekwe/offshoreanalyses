import test from 'node:test';
import assert from 'node:assert/strict';
import * as Lw from '../../src/engine/lines/lowering.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);
const g = 9.80665;

test('static hook load vs depth: hand calculation (object weight in water + submerged cable weight)', () => {
  const s = Lw.staticTension({ mass: 100e3, volume: 30, cableWeight: 250, depths: [0, 500, 2000] });
  const W = (100e3 - 1025 * 30) * g;
  assert.ok(rel(s.objectWeightInWater, W) < 1e-12);
  assert.ok(rel(s.rows[1].tension, W + 250 * 500) < 1e-12 && rel(s.rows[2].tension, W + 250 * 2000) < 1e-12);
});

test('axial resonance: exact transcendental equation vs spring–mass with one-third of the cable mass, and limiting cases', () => {
  const EA = 6e8, L = 1500, m = 30, M = 150e3, A33 = 50e3;
  const r = Lw.axialResonance({ mass: M, addedMass: A33, EA, length: L, cableMass: m });
  // closed-form lumped value written out
  const wl = Math.sqrt((EA / L) / (M + A33 + m * L / 3));
  assert.ok(rel(r.omegaLumped, wl) < 1e-12);
  // exact root satisfies (ωL/c) tan(ωL/c) = mL/M'
  const x = r.omegaExact * L / Math.sqrt(EA / m);
  assert.ok(Math.abs(x * Math.tan(x) - m * L / (M + A33)) < 1e-9);
  // the one-third rule is accurate to better than 1 % for cable mass ≤ 25 % of the end mass
  assert.ok(rel(r.periodLumped, r.periodExact) < 0.002, `${r.periodLumped} vs ${r.periodExact}`);
  // massless cable → plain spring–mass
  const r0 = Lw.axialResonance({ mass: M, addedMass: A33, EA, length: L, cableMass: 0 });
  assert.ok(rel(r0.periodExact, 2 * Math.PI * Math.sqrt((M + A33) * L / EA)) < 1e-12);
  // very light end mass → fixed–free rod, ωL/c → π/2
  const rr = Lw.axialResonance({ mass: 1e-3, EA, length: L, cableMass: m });
  assert.ok(rel(rr.omegaExact, (Math.PI / 2) * Math.sqrt(EA / m) / L) < 1e-4);
  // soft spring in series: 1/K = 1/kc + L/EA
  const rs = Lw.axialResonance({ mass: M, addedMass: A33, EA, length: L, cableMass: m, kc: 2e5 });
  assert.ok(rel(rs.stiffness, 1 / (1 / 2e5 + L / EA)) < 1e-12 && rs.periodLumped > r.periodLumped);
});

test('crane-tip heave response: undamped SDOF closed form, drag-limited resonance and slack-wire detection', () => {
  const base = { mass: 150e3, addedMass: 50e3, EA: 6e8, length: 1500, cableMass: 30, heaveAmp: 1.0, staticTension: 1.2e6 };
  const K = 6e8 / 1500, Mp = 200e3 + 15e3;
  const T = 9, om = 2 * Math.PI / T;
  const u = Lw.heaveResponse({ ...base, period: T });
  // undamped: x/η = 1/(1 − (ω/ωn)²), Fd = K η (ω/ωn)²/|1 − (ω/ωn)²|
  const r2 = Mp * om * om / K;
  assert.ok(rel(u.motionRatio, 1 / Math.abs(1 - r2)) < 1e-9);
  assert.ok(rel(u.dynamicForce, K * r2 / Math.abs(1 - r2)) < 1e-9);
  assert.ok(rel(u.daf, 1 + u.dynamicForce / 1.2e6) < 1e-12);
  // at resonance with quadratic drag the linearised amplitude solves K η = ω c x, c = (8/3π)·½ρ Cd A ω x
  const Tn = u.naturalPeriod, on = 2 * Math.PI / Tn;
  const d = Lw.heaveResponse({ ...base, period: Tn, Cd: 2.0, area: 40 });
  const xRes = Math.sqrt(K * 1.0 / (on * on * (8 / (3 * Math.PI)) * 0.5 * 1025 * 2.0 * 40));
  assert.ok(rel(d.motionAmp, xRes) < 1e-6, `${d.motionAmp} vs ${xRes}`);
  assert.ok(d.slack === (d.dynamicForce >= 1.2e6) && rel(d.minTension, 1.2e6 - d.dynamicForce) < 1e-12);
  // long period: object follows the crane tip, negligible dynamic force
  const q = Lw.heaveResponse({ ...base, period: 200 });
  assert.ok(Math.abs(q.motionRatio - 1) < 1e-3 && q.dynamicForce < 0.01 * 1.2e6 && !q.slack);
});

test('heave compensator transmissibility (base-excited spring–damper–mass closed form) and pendulum period', () => {
  const M = 200e3, k = 5e4, c = 3e4, T = 8, om = 2 * Math.PI / T;
  const h = Lw.heaveCompensator({ mass: M, k, c, period: T });
  assert.ok(rel(h.transmissibility, Math.sqrt((k * k + (c * om) ** 2) / ((k - M * om * om) ** 2 + (c * om) ** 2))) < 1e-12);
  assert.ok(h.transmissibility < 1 && rel(h.efficiency, 1 - h.transmissibility) < 1e-12);
  // at ω = sqrt(2)·ωn the transmissibility is exactly 1 for any damping
  const Tcross = 2 * Math.PI / Math.sqrt(2 * k / M);
  assert.ok(Math.abs(Lw.heaveCompensator({ mass: M, k, c, period: Tcross }).transmissibility - 1) < 1e-12);
  assert.ok(rel(Lw.heaveCompensator({ mass: M, k, c, period: T, activeResidual: 0.2 }).transmissibility, 0.2 * h.transmissibility) < 1e-12);
  assert.ok(rel(Lw.pendulumPeriod({ length: 40 }), 2 * Math.PI * Math.sqrt(40 / g)) < 1e-12);
  assert.ok(rel(Lw.pendulumPeriod({ length: 40, amplitudeDeg: 20 }) / Lw.pendulumPeriod({ length: 40 }), 1 + (20 * Math.PI / 180) ** 2 / 16) < 1e-12);
});

test('horizontal offset in uniform current vs closed form x = (q/w)L + ((FD0 − qW/w)/w) ln(1 + wL/W)', () => {
  const L = 1000, W = 8e5, w = 200, dc = 0.08, U = 0.8, CdA = 30;
  const o = Lw.currentOffset({ length: L, objectWeight: W, cableWeight: w, cableDiameter: dc, objectDragArea: CdA, current: U });
  const q = 0.5 * 1025 * 1.2 * dc * U * U, FD0 = 0.5 * 1025 * CdA * U * U;
  const exact = (q / w) * L + ((FD0 - q * W / w) / w) * Math.log(1 + w * L / W);
  assert.ok(rel(o.offset, exact) < 1e-5, `${o.offset} vs ${exact}`);
  assert.ok(rel(o.objectDrag, FD0) < 1e-12);
  assert.ok(rel(Math.tan(o.bottomAngleDeg * Math.PI / 180), FD0 / W) < 1e-9);
  // sheared current (zero below 500 m) gives a smaller offset and no object drag
  const sh = Lw.currentOffset({ length: L, objectWeight: W, cableWeight: w, cableDiameter: dc, objectDragArea: CdA, current: (d) => (d < 500 ? U : 0) });
  assert.ok(sh.offset < o.offset && sh.objectDrag === 0);
  const lv = Lw.landingVelocity({ loweringSpeed: 0.3, objectMotionAmp: 0.5, period: 10 });
  assert.ok(rel(lv.velocity, 0.3 + 2 * Math.PI / 10 * 0.5) < 1e-12);
});

test('DNV-RP-N103 Simplified Method: crane-tip motion, wave kinematics, force components and slack criterion (hand calculation)', () => {
  const ct = Lw.craneTipMotion({ heaveAmp: 1.0, rollAmpDeg: 3, pitchAmpDeg: 2, b: 25, l: 40, T3: 8, T4: 12, T5: 7 });
  const r = 25 * Math.sin(3 * Math.PI / 180), p = 40 * Math.sin(2 * Math.PI / 180);
  assert.ok(rel(ct.eta, Math.sqrt(1 + r * r + p * p)) < 1e-12);
  assert.ok(rel(ct.velocity, 2 * Math.PI * Math.sqrt((1 / 8) ** 2 + (r / 12) ** 2 + (p / 7) ** 2)) < 1e-12);
  assert.ok(rel(ct.acceleration, 4 * Math.PI ** 2 * Math.sqrt((1 / 64) ** 2 + (r / 144) ** 2 + (p / 49) ** 2)) < 1e-12);
  const wk = Lw.waveKinematicsN103({ Hs: 2.0, Tz: 6, d: 3 });
  assert.ok(rel(wk.zetaA, 1.8) < 1e-12);
  const dec = Math.exp(-4 * Math.PI ** 2 * 3 / (36 * g));
  assert.ok(rel(wk.velocity, 1.8 * (2 * Math.PI / 6) * dec) < 1e-12 && rel(wk.acceleration, 1.8 * (2 * Math.PI / 6) ** 2 * dec) < 1e-12);
  assert.equal(wk.tzInRange, 6 >= 8.9 * Math.sqrt(2 / g));
  const inp = { mass: 60e3, volume: 20, addedMass: 40e3, slamArea: 12, Cs: 5, dragArea: 25, CD: 2.5, waterplaneArea: 6, vc: 0.5,
    etaCt: ct.eta, vCt: ct.velocity, aCt: ct.acceleration, zetaA: wk.zetaA, vW: wk.velocity, aW: wk.acceleration };
  const s = Lw.splashZoneN103(inp);
  const vs = 0.5 + Math.hypot(ct.velocity, wk.velocity);
  const Fslam = 0.5 * 1025 * 5 * 12 * vs * vs, FD = 0.5 * 1025 * 2.5 * 25 * vs * vs;
  const Frho = 1025 * 6 * Math.hypot(1.8, ct.eta) * g;
  const FM = Math.hypot(100e3 * ct.acceleration, (1025 * 20 + 40e3) * wk.acceleration);
  const Fhyd = Math.hypot(FD + Fslam, FM - Frho), stat = (60e3 - 1025 * 20) * g;
  assert.ok(rel(s.Fslam, Fslam) < 1e-12 && rel(s.FD, FD) < 1e-12 && rel(s.Frho, Frho) < 1e-12 && rel(s.FM, FM) < 1e-12);
  assert.ok(rel(s.Fhyd, Fhyd) < 1e-12);
  assert.ok(rel(s.slackUc, Fhyd / (0.9 * stat)) < 1e-12 && s.slackOk === (Fhyd <= 0.9 * stat));
  assert.ok(rel(s.dafConv, (stat + Fhyd) / (60e3 * g)) < 1e-12);
  // calm sea, no vessel motion: only the lowering-speed drag and slam remain
  const calm = Lw.splashZoneN103({ ...inp, etaCt: 0, vCt: 0, aCt: 0, zetaA: 0, vW: 0, aW: 0 });
  assert.ok(rel(calm.Fhyd, 0.5 * 1025 * 0.25 * (5 * 12 + 2.5 * 25)) < 1e-12);
  // snap load: vff/vr = 0.45 → C = cos(π·0.25)
  const sn = Lw.snapLoadN103({ mass: 60e3, addedMass: 40e3, K: 2e6, vff: 0.9, vr: 2.0 });
  assert.ok(rel(sn.C, Math.cos(Math.PI * 0.25)) < 1e-12);
  assert.ok(rel(sn.Fsnap, (0.9 + sn.C * 2.0) * Math.sqrt(2e6 * 100e3)) < 1e-12);
  assert.equal(Lw.snapLoadN103({ mass: 1, K: 1, vff: 0.1, vr: 1 }).C, 1);
  assert.equal(Lw.snapLoadN103({ mass: 1, K: 1, vff: 0.9, vr: 1 }).C, 0);
  assert.throws(() => Lw.splashZoneN103({ ...inp, mass: -1 }), /argument mass/);
});
