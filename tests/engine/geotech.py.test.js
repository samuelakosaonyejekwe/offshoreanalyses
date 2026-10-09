import test from 'node:test';
import assert from 'node:assert/strict';
import {
  matlockPu, matlockTransitionDepth, pySoftClay, pyStiffClayNoFreeWater, sandCoefficients, sandInitialModulus,
  sandPu, sandTransitionDepth, pySand, pyWeakRock, pyJeanjean, pyUser, pyLinear, scalePy, samplePy,
  initialStiffness, pyCurveAt, zeroCurve,
} from '../../src/engine/geotech/py.js';
import { createProfile, applyScour } from '../../src/engine/geotech/soilprofile.js';

const close = (a, b, rel = 1e-9) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-12), `${a} vs ${b}`);

test('Matlock (1970) static curve: p/pu = 0.5 (y/yc)^(1/3), pu = 9 su D at depth, yc = 2.5 eps50 D', () => {
  const su = 40e3, D = 1.2, eps50 = 0.01;
  const c = pySoftClay({ su, sigmaV: 30 * 8000, z: 30, D, eps50 });
  close(c.pu, 9 * su * D);
  const yc = 2.5 * eps50 * D;
  close(c.meta.yc, yc);
  close(c.p(yc), 0.5 * c.pu);
  close(c.p(2 * yc), 0.5 * Math.cbrt(2) * c.pu);
  close(c.p(8 * yc), c.pu);
  close(c.p(20 * yc), c.pu);
  close(c.p(-yc), -0.5 * c.pu); // odd
  // API tabulated points p/pu = 0.23, 0.33, 0.50, 0.72, 1.00 at y/yc = 0.1, 0.3, 1, 3, 8 (two-decimal rounding)
  for (const [a, r] of [[0.1, 0.23], [0.3, 0.33], [1, 0.5], [3, 0.72], [8, 1.0]]) assert.ok(Math.abs(c.p(a * yc) / c.pu - r) < 0.006);
  // tangent equals the numerical derivative
  const h = 1e-7;
  close(c.k(2 * yc), (c.p(2 * yc + h) - c.p(2 * yc - h)) / (2 * h), 1e-5);
});

test('Matlock shallow resistance and transition depth XR = 6D/(gamma D/su + J)', () => {
  const su = 20e3, g = 6000, D = 1, J = 0.5;
  const XR = matlockTransitionDepth({ su, gammaSub: g, D, J });
  close(XR, 7.5); // 6 / (0.3 + 0.5)
  close(matlockPu({ su, sigmaV: 0, z: 0, D, J }).pu, 3 * su * D);
  // hand: at z = 3 m, Np = 3 + 18/20 + 0.5·3 = 5.4
  const u = matlockPu({ su, sigmaV: g * 3, z: 3, D, J });
  close(u.Np, 5.4);
  assert.equal(u.shallow, true);
  close(u.zOverZR, 3 / XR);
  close(matlockPu({ su, sigmaV: g * XR, z: XR, D, J }).Np, 9);
  assert.equal(matlockPu({ su, sigmaV: g * 12, z: 12, D, J }).shallow, false);
});

test('Matlock cyclic curve: 0.72 pu plateau below XR, residual 0.72 pu X/XR at 15 yc above XR', () => {
  const su = 20e3, g = 6000, D = 1, eps50 = 0.02;
  const yc = 2.5 * eps50 * D, peak = 0.5 * Math.cbrt(3);
  assert.ok(Math.abs(peak - 0.72) < 0.002);
  const deep = pySoftClay({ su, sigmaV: g * 10, z: 10, D, eps50, loading: 'cyclic' });
  close(deep.p(3 * yc), peak * deep.pu);
  close(deep.p(30 * yc), peak * deep.pu);
  const sh = pySoftClay({ su, sigmaV: g * 3, z: 3, D, eps50, loading: 'cyclic' });
  close(sh.p(3 * yc), peak * sh.pu);
  close(sh.p(15 * yc), peak * sh.pu * (3 / 7.5));
  close(sh.p(40 * yc), peak * sh.pu * (3 / 7.5));
  close(sh.p(9 * yc), peak * sh.pu * (1 - 0.5 * (1 - 3 / 7.5))); // halfway down the linear branch
  assert.ok(sh.k(9 * yc) < 0);
  // same as static below 3 yc
  const st = pySoftClay({ su, sigmaV: g * 3, z: 3, D, eps50 });
  close(sh.p(2 * yc), st.p(2 * yc));
});

test('API sand coefficients from the Reese wedge expressions: hand values at 30 deg', () => {
  const { C1, C2, C3 } = sandCoefficients(30);
  // beta = 60°, alpha = 15°, Ka = 1/3, K0 = 0.4
  close(C2, 3 - 1 / 3, 1e-12); // tan60/tan30 − Ka
  close(C3, (81 - 1) / 3 + 0.4 * Math.tan(Math.PI / 6) * 9, 1e-12);
  close(C1, 0.358630 + 1.392305 + 0.160770, 2e-6);
  // read-off values quoted from the API figure: C1 ≈ 1.9, C2 ≈ 2.7, C3 ≈ 28 at 30°; ≈ 4.6, 4.3, 100 at 40°
  const c40 = sandCoefficients(40);
  assert.ok(Math.abs(c40.C1 - 4.6) < 0.15 && Math.abs(c40.C2 - 4.3) < 0.15 && Math.abs(c40.C3 - 100) < 6, JSON.stringify(c40));
  assert.throws(() => sandCoefficients(50), RangeError);
});

test('API sand: pu continuous at the transition depth; tanh initial slope = k z; A factor', () => {
  const phiDeg = 35, D = 1.5, g = 10e3;
  const zt = sandTransitionDepth(phiDeg, D);
  const u = sandPu({ phiDeg, sigmaV: g * zt, z: zt, D });
  close(u.pus, u.pud, 1e-12);
  const a = sandPu({ phiDeg, sigmaV: g * zt * 0.999, z: zt * 0.999, D }), b = sandPu({ phiDeg, sigmaV: g * zt * 1.001, z: zt * 1.001, D });
  assert.equal(a.shallow, true); assert.equal(b.shallow, false);
  assert.ok(Math.abs(a.pu - b.pu) / u.pu < 0.01);
  const k = 20e6, z = 4;
  const c = pySand({ phiDeg, sigmaV: g * z, z, D, k });
  close(c.k(0), k * z);
  close(c.p(1e-9) / 1e-9, k * z, 1e-6);
  close(c.meta.A, 0.9); // 3 − 0.8·4/1.5 = 0.87 is below the floor of 0.9
  close(pySand({ phiDeg, sigmaV: g * 1, z: 1, D, k }).meta.A, 3 - 0.8 / 1.5);
  close(pySand({ phiDeg, sigmaV: g * 20, z: 20, D, k }).meta.A, 0.9);
  close(pySand({ phiDeg, sigmaV: g * 1, z: 1, D, k, loading: 'cyclic' }).meta.A, 0.9);
  // asymptote A pu, closed-form value at one point
  close(c.p(10), c.meta.A * c.meta.puApi, 1e-9);
  close(c.p(0.01), c.pu * Math.tanh((k * z * 0.01) / c.pu));
  // zero resistance at the surface
  assert.equal(pySand({ phiDeg, sigmaV: 0, z: 0, D, k }).p(0.1), 0);
});

test('k(phi) fit: closed form, and the medium / dense submerged values of Reese et al. (16.3, 34 MN/m3)', () => {
  assert.ok(Math.abs(sandInitialModulus(33) / 1e6 - 16.3) < 0.3);
  assert.ok(Math.abs(sandInitialModulus(38) / 1e6 - 34) < 0.3);
  close(sandInitialModulus(35), (0.008085 * 35 ** 2.45 - 26.09) * 1e6);
  for (let phi = 29; phi < 45; phi++) assert.ok(sandInitialModulus(phi + 1) > sandInitialModulus(phi));
  assert.throws(() => sandInitialModulus(25), RangeError);
});

test('stiff clay without free water (Reese & Welch): 0.5 pu at y50, pu at 16 y50, cyclic shift', () => {
  const o = { su: 120e3, sigmaV: 200e3, z: 20, D: 1, eps50: 0.005 };
  const c = pyStiffClayNoFreeWater(o);
  const y50 = 2.5 * 0.005 * 1;
  close(c.p(y50), 0.5 * c.pu);
  close(c.p(16 * y50), c.pu);
  // cyclic: y_c = y_s + y50 · 9.6 (p/pu)^4 log10 N; at p = 0.5 pu, y_s = y50 → y_c = y50 (1 + 0.6 log10 N)
  const cy = pyStiffClayNoFreeWater({ ...o, loading: 'cyclic', nCycles: 1000 });
  close(cy.p(y50 * (1 + 0.6 * 3)), 0.5 * c.pu);
});

test('weak rock (Reese 1997): continuity at yA and ultimate resistance', () => {
  const o = { qur: 2e6, Eir: 500e6, D: 1.5, krm: 0.0005 };
  const sh = pyWeakRock({ ...o, z: 1.5 });
  close(sh.pu, 2e6 * 1.5 * (1 + 1.4)); // xr = b
  close(sh.meta.Kir, (100 + 400 / 3) * 500e6);
  const yA = sh.meta.yA;
  close(sh.p(yA * (1 - 1e-9)), sh.p(yA * (1 + 1e-9)), 1e-6);
  close(pyWeakRock({ ...o, z: 10 }).pu, 5.2 * 2e6 * 1.5);
  close(sh.p(1), sh.pu);
});

test('Jeanjean (2009): Np from 8 at the mudline to 12 at depth; tanh form', () => {
  const D = 1, su0 = 5e3, su1 = 1.5e3;
  const top = pyJeanjean({ su: su0, su0, su1, Gmax: 550 * su0, z: 0, D });
  close(top.meta.Np, 8);
  close(top.meta.xi, 0.25 + 0.05 * (su0 / (su1 * D)));
  const deep = pyJeanjean({ su: su0 + 100 * su1, su0, su1, Gmax: 550 * (su0 + 100 * su1), z: 100, D });
  assert.ok(Math.abs(deep.meta.Np - 12) < 1e-6);
  close(deep.p(0.04), deep.pu * Math.tanh(5.5 * Math.sqrt(0.04)));
});

test('user curves, linear springs, multipliers, sampling and initial stiffness', () => {
  const c = pyUser({ y: [0, 0.01, 0.05], p: [0, 100, 200] });
  close(c.p(0.005), 50); close(c.p(0.03), 150); close(c.p(1), 200); close(c.p(-0.03), -150);
  close(c.k(0.03), 2500); assert.equal(c.k(1), 0);
  const s = scalePy(c, { pMult: 0.5, yMult: 2 });
  close(s.p(0.06), 0.5 * 150); close(s.k(0.06), (0.5 / 2) * 2500); close(s.pu, 100);
  close(pyLinear(3e6).p(0.02), 6e4);
  close(initialStiffness(c), 1e4, 1e-9);
  const smp = samplePy(c, { yMax: 0.05, n: 11 });
  assert.equal(smp.y.length, 11); close(smp.p[10], 200);
  assert.equal(zeroCurve().p(1), 0);
  assert.throws(() => pyUser({ y: [0, 1, 1], p: [0, 1, 2] }), RangeError);
});

test('curve selection from a layered profile, with scour', () => {
  const p = createProfile([
    { type: 'sand', top: 0, thickness: 10, gammaSub: 10e3, phiDeg: 35 },
    { type: 'clay', top: 10, thickness: 30, gammaSub: 8e3, suTop: 150e3, suGrad: 0, eps50: 0.005 },
  ]);
  const D = 2;
  const sand = pyCurveAt(p, 5, { D });
  assert.equal(sand.type, 'sand');
  close(sand.meta.kz, sandInitialModulus(35) * 5);
  const clay = pyCurveAt(p, 30, { D });
  close(clay.pu, 9 * 150e3 * D);
  // stiff clay under cyclic loading: explicit degradation factor on the Matlock curve
  const deg = pyCurveAt(p, 30, { D, loading: 'cyclic', stiffClayCyclicFactor: 0.5 });
  close(deg.pu, 0.5 * clay.pu);
  const ps = applyScour(p, { general: 0, D });
  assert.equal(pyCurveAt(ps, 2, { D }).pu, 0);
  const sc = pyCurveAt(ps, 5, { D });
  close(sc.meta.kz, sandInitialModulus(35) * 2); // depth measured from the scour base at 3 m
  assert.ok(sc.pu < sand.pu);
  close(pyCurveAt(p, 5, { D, pMult: 0.4 }).pu, 0.4 * sand.pu);
});
