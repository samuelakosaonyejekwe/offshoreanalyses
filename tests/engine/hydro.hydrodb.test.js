import test from 'node:test';
import assert from 'node:assert/strict';
import * as db from '../../src/engine/hydro/hydrodb.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);
const near = (a, b, abs, msg) => assert.ok(Math.abs(a - b) <= abs, `${msg || ''} got ${a}, expected ${b}`);
const z6 = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);

test('WAMIT .1: A = rho L^k Abar, B = rho L^k omega Bbar with k = 3, 4, 5; zero/infinite period records', () => {
  const text = [
    '  -1.0  3 3  0.90',
    '   0.0  3 3  0.40',
    '  10.0  1 1  0.20  0.05',
    '  10.0  3 3  0.50  0.10',
    '  10.0  1 5  0.03  0.01',
    '  10.0  5 5  0.07  0.02',
    '   5.0  3 3  0.45  0.30',
  ].join('\n');
  const r = db.parseWamit1(text, { rho: 1000, L: 2 });
  assert.equal(r.omega.length, 2);
  close(r.omega[0], (2 * Math.PI) / 10, 1e-12);       // sorted by frequency: T = 10 s first
  close(r.A[0][2][2], 1000 * 8 * 0.5, 1e-12);
  close(r.B[0][2][2], 1000 * 8 * r.omega[0] * 0.1, 1e-12);
  close(r.A[0][0][4], 1000 * 16 * 0.03, 1e-12);
  close(r.B[0][4][4], 1000 * 32 * r.omega[0] * 0.02, 1e-12);
  close(r.A[1][2][2], 1000 * 8 * 0.45, 1e-12);
  close(r.Ainf[2][2], 1000 * 8 * 0.4, 1e-12); close(r.A0[2][2], 1000 * 8 * 0.9, 1e-12);
});

test('WAMIT .3 and .4: X = rho g L^m Xbar (m = 2 forces, 3 moments); RAO rotations divided by L', () => {
  const t3 = ['10.0  0.0 3  1.0 0.0  0.6 -0.8', '10.0  0.0 5  0.5 90.0  0.0 0.5', '10.0 90.0 3  1.0 0.0  1.0 0.0', '10.0 90.0 5 0 0 0 0'].join('\n');
  const f = db.parseWamit3(t3, { rho: 1000, g: 10, L: 2 });
  assert.deepEqual(f.headingsDeg, [0, 90]);
  close(f.Fre[0][0][2], 1000 * 10 * 4 * 0.6, 1e-12); close(f.Fim[0][0][2], -1000 * 10 * 4 * 0.8, 1e-12);
  close(f.Fim[0][0][4], 1000 * 10 * 8 * 0.5, 1e-12);
  const r4 = db.parseWamit4('8.0 0.0 3 1 0 0.9 0.1\n8.0 0.0 5 1 0 0.4 0.0', { L: 2 });
  close(r4.re[0][0][2], 0.9, 1e-12); close(r4.re[0][0][4], 0.2, 1e-12);
});

function heaveDb() {
  // single-DOF (heave) data, axisymmetric: excitation independent of heading
  const rho = 1025; const g = 9.80665;
  const omega = [0.3, 0.5, 0.7, 0.9, 1.1];
  const headingsDeg = [0, 45, 90, 135, 180, 225, 270, 315];
  const X0 = omega.map((w) => 4e6 * Math.exp(-w));
  const A = omega.map((w) => { const m = z6(); m[2][2] = 5e6 + 1e6 * w; return m; });
  // axisymmetric heave identity (deep water): B33 = omega³ |X3|² / (2 rho g³)
  const B = omega.map((w, i) => { const m = z6(); m[2][2] = (w ** 3 * X0[i] ** 2) / (2 * rho * g ** 3); return m; });
  const Fre = headingsDeg.map(() => omega.map((w, i) => [0, 0, X0[i] * Math.cos(0.3), 0, 0, 0]));
  const Fim = headingsDeg.map(() => omega.map((w, i) => [0, 0, X0[i] * Math.sin(0.3), 0, 0, 0]));
  return db.createHydroDb({ omega, A, B, headingsDeg, Fre, Fim, rho, g });
}

test('RAO solve reproduces the single-degree-of-freedom closed form xi = F / (C - w²(M + A) + i w B)', () => {
  const d = heaveDb();
  const M = z6(); const C = z6();
  for (let i = 0; i < 6; i++) { M[i][i] = 1e7; C[i][i] = 1e7; }
  C[2][2] = 6e6;
  const r = db.raoSolve({ db: d, M, C, headingDeg: 20 });
  d.omega.forEach((w, i) => {
    const re = C[2][2] - w * w * (M[2][2] + d.A[i][2][2]); const im = w * d.B[i][2][2];
    const X = 4e6 * Math.exp(-w);
    close(r.amp[2][i], X / Math.hypot(re, im), 1e-10);
    near(r.phase[2][i], 0.3 - Math.atan2(im, re), 1e-10);
  });
  // interpolated frequency: coefficients are linear between tabulated points
  const h = db.hydroAt(d, 0.6, 20);
  close(h.A[2][2], 5e6 + 1e6 * 0.6, 1e-12);
  close(h.B[2][2], 0.5 * (d.B[1][2][2] + d.B[2][2][2]), 1e-12);
});

test('Heading interpolation: linear in real/imaginary parts and periodic over the full circle', () => {
  const omega = [0.5, 1.0];
  const headingsDeg = [0, 90, 180, 270];
  const mk = (v) => omega.map(() => [v, 0, 0, 0, 0, 0]);
  const d = db.createHydroDb({ omega, A: omega.map(z6), B: omega.map(z6), headingsDeg, Fre: [mk(1), mk(0), mk(-1), mk(0)], Fim: [mk(0), mk(2), mk(0), mk(-2)] });
  assert.ok(d.periodicHeading);
  close(db.hydroAt(d, 0.7, 45).F.re[0], 0.5, 1e-12); close(db.hydroAt(d, 0.7, 45).F.im[0], 1, 1e-12);
  close(db.hydroAt(d, 0.7, 315).F.re[0], 0.5, 1e-12); close(db.hydroAt(d, 0.7, -45).F.im[0], -1, 1e-12);
  close(db.hydroAt(d, 0.7, 405).F.re[0], 0.5, 1e-12);
});

test('Consistency checks: Haskind damping–excitation identity (axisymmetric heave), negative damping, asymmetry', () => {
  const d = heaveDb();
  const ok = db.checkHydroDb(d);
  assert.ok(ok.ok, JSON.stringify(ok.haskind));
  assert.ok(ok.haskind.find((h) => h.dof === 3).maxRelError < 1e-10);
  // halve the damping: violates the identity
  const bad = heaveDb(); bad.B.forEach((m) => { m[2][2] *= 0.5; });
  assert.ok(!db.checkHydroDb(bad).ok);
  const neg = heaveDb(); neg.B[1][2][2] = -100; neg.A[0][0][0] = 1e6; neg.A[0][0][2] = 3e5; neg.A[0][2][0] = 1e5;
  const c = db.checkHydroDb(neg);
  assert.equal(c.negativeDamping.length, 1); assert.equal(c.negativeDamping[0].dof, 3);
  close(c.asymmetryA, 2e5 / Math.sqrt(1e6 * 5.3e6), 1e-9);
});

test('JSON and CSV schemas round-trip to the same database', () => {
  const d = heaveDb();
  const j = db.parseHydroJson(JSON.stringify({ rho: d.rho, g: d.g, omega: d.omega, headingsDeg: d.headingsDeg, A: d.A, B: d.B, Fre: d.Fre, Fim: d.Fim }));
  close(j.B[3][2][2], d.B[3][2][2], 1e-15); close(j.Fim[5][2][2], d.Fim[5][2][2], 1e-15);
  const rows = ['kind,omega,headingDeg,i,j,re,im', 'AINF,,,3,3,4.5e6,'];
  d.omega.forEach((w, i) => {
    rows.push(`A,${w},,3,3,${d.A[i][2][2]},`); rows.push(`B,${w},,3,3,${d.B[i][2][2]},`);
    d.headingsDeg.forEach((h, ih) => rows.push(`F,${w},${h},3,,${d.Fre[ih][i][2]},${d.Fim[ih][i][2]}`));
  });
  const c = db.parseHydroCsv(rows.join('\n'), { rho: d.rho, g: d.g });
  close(c.A[2][2][2], d.A[2][2][2], 1e-15); close(c.Fre[7][4][2], d.Fre[7][4][2], 1e-15); close(c.Ainf[2][2], 4.5e6, 1e-15);
  assert.ok(db.checkHydroDb(c).ok);
});

test('hydroDbFromWamit assembles .1 and .3 into a database', () => {
  const t1 = '10.0 3 3 0.5 0.1\n5.0 3 3 0.45 0.3';
  const t3 = '10.0 0.0 3 1 0 1.0 0.0\n5.0 0.0 3 1 0 0.5 0.5';
  const d = db.hydroDbFromWamit({ text1: t1, text3: t3, rho: 1000, g: 10, L: 1 });
  assert.equal(d.omega.length, 2); close(d.Fre[0][1][2], 5000, 1e-12); close(d.A[1][2][2], 450, 1e-12);
  assert.throws(() => db.hydroDbFromWamit({ text1: t1, text3: '9.0 0.0 3 1 0 1 0\n5.0 0.0 3 1 0 1 0' }), RangeError);
});
