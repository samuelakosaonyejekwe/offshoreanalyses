import test from 'node:test';
import assert from 'node:assert/strict';
import { streamFunctionWave, streamFunctionKinematics, streamFunctionEta, streamFunctionProfile } from '../../src/engine/metocean/streamfunction.js';
import { airyKinematics, wavelength, wavenumber, stokes5Wave, stokes5Kinematics, G } from '../../src/engine/metocean/waves.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

test('stream function reduces to Airy theory for a small wave (length, surface, kinematics, pressure)', () => {
  const H = 0.002; const T = 8; const d = 20;
  const s = streamFunctionWave({ H, T, d });
  assert.ok(rel(s.L, wavelength(T, d)) < 1e-7);
  assert.ok(rel(s.T, T) < 1e-12);
  assert.ok(Math.abs(s.currentEulerian) < 1e-12);
  for (const [x, z, t] of [[0, -1, 0], [7, -5, 1.3], [-20, -19, 3], [33, -0.5, 6.1]]) {
    const a = airyKinematics({ H, T, d, x, z, t });
    const f = streamFunctionKinematics(s, { x, z, t });
    // differences are of relative order kH/2 = 1e-4
    for (const key of ['eta', 'u', 'w', 'ax', 'az', 'pDyn']) {
      assert.ok(Math.abs(f[key] - a[key]) < 5e-4 * Math.max(Math.abs(a[key]), 1e-3 * Math.abs(a.u) + (key === 'pDyn' ? 1 : 0)), `${key} at x=${x}`);
    }
  }
});

test('stream function satisfies both free-surface boundary conditions to better than 1e-6', () => {
  const s = streamFunctionWave({ H: 5, T: 10, d: 30, order: 16 });
  assert.equal(s.report.converged, true);
  assert.ok(s.report.residual < 1e-10);
  assert.ok(s.report.kinematicResidual < 1e-6);
  assert.ok(s.report.dynamicResidual < 1e-6);
  // independent check through the public kinematics, at points that are not collocation points
  const rho = 1025; const n = 37;
  for (let i = 0; i < n; i++) {
    const x = s.L * (i / (n - 1) - 0.5) * 0.999;
    const eta = streamFunctionEta(s, x, 0);
    const k = streamFunctionKinematics(s, { x, z: eta, rho });
    // dynamic condition: zero gauge pressure on the free surface
    assert.ok(Math.abs(k.p) / (rho * G * s.H) < 1e-6, `pressure at x=${x}`);
    // kinematic condition: the surface is a streamline in the wave frame, w = (u − c) dη/dx
    const h = 1e-4 * s.L;
    const ex = (streamFunctionEta(s, x + h, 0) - streamFunctionEta(s, x - h, 0)) / (2 * h);
    assert.ok(Math.abs(k.w - (k.u - s.c) * ex) / s.c < 1e-6, `kinematic at x=${x}`);
  }
  // bed is impermeable
  assert.ok(Math.abs(streamFunctionKinematics(s, { x: 11, z: -30 }).w) < 1e-12);
});

test('stream function: wave height, mean level, period and field equations', () => {
  const s = streamFunctionWave({ H: 6, T: 11, d: 18, order: 20 });
  assert.ok(rel(s.crest - s.trough, 6) < 1e-11);
  assert.ok(rel(s.T, 11) < 1e-11);
  assert.ok(s.crest > 3 && s.crest < 6);
  const p = streamFunctionProfile(s, 401);
  let mean = 0;
  for (let i = 0; i < 400; i++) mean += 0.5 * (p.eta[i] + p.eta[i + 1]);
  assert.ok(Math.abs(mean / 400) < 1e-6 * s.H);
  assert.ok(Math.abs(p.eta[200] - s.crest) < 1e-9 && Math.abs(p.eta[0] - s.trough) < 1e-9);
  // continuity and irrotationality, and local acceleration = ∂u/∂t
  const h = 1e-3;
  const f = (dx, dz, dt = 0) => streamFunctionKinematics(s, { x: 12 + dx, z: -6 + dz, t: 0.7 + dt });
  assert.ok(Math.abs((f(h, 0).u - f(-h, 0).u + f(0, h).w - f(0, -h).w) / (2 * h)) < 1e-7);
  assert.ok(Math.abs((f(0, h).u - f(0, -h).u - f(h, 0).w + f(-h, 0).w) / (2 * h)) < 1e-7);
  assert.ok(Math.abs((f(0, 0, h).u - f(0, 0, -h).u) / (2 * h) - f(0, 0).ax) < 1e-6);
  assert.ok(Math.abs((f(0, 0, h).w - f(0, 0, -h).w) / (2 * h) - f(0, 0).az) < 1e-6);
  // above the surface: dry
  assert.equal(streamFunctionKinematics(s, { x: s.L / 2, z: 0 }).wet, false);
  assert.throws(() => streamFunctionKinematics(s, { x: 0, z: -19 }), RangeError);
});

test('stream function agrees with the fifth-order Stokes theory (Fenton 1985) for a moderately steep wave', () => {
  const o = { H: 5, T: 10, d: 30 };
  const s = streamFunctionWave({ ...o, order: 16 });
  const s5 = stokes5Wave(o);
  assert.ok(rel(s.L, s5.L) < 1e-5);
  assert.ok(rel(s.crest, s5.crest) < 2e-4);
  assert.ok(rel(s.currentStokes, s5.currentStokes) < 3e-3);
  const a = streamFunctionKinematics(s, { x: 0, z: 0 });
  const b = stokes5Kinematics(s5, { x: 0, z: 0 });
  assert.ok(rel(a.u, b.u) < 2e-4 && rel(a.azTotal, b.azTotal) < 2e-4 && rel(a.pDyn, b.pDyn) < 2e-4);
});

test('stream function: mean-current closures (Eulerian current, zero mass transport)', () => {
  const e = streamFunctionWave({ H: 4, T: 9, d: 25, currentEulerian: 1.0 });
  assert.ok(rel(e.currentEulerian, 1.0) < 1e-11);
  assert.ok(e.currentStokes > 1.0);
  // time-mean velocity at a point below the trough equals the Eulerian current
  let m = 0; const n = 64;
  for (let i = 0; i < n; i++) m += streamFunctionKinematics(e, { x: 0, z: -15, t: (i * e.T) / n }).u;
  assert.ok(Math.abs(m / n - 1.0) < 1e-9);
  // a following current lengthens the wave; linear Doppler estimate within the nonlinear correction
  const e0 = streamFunctionWave({ H: 4, T: 9, d: 25 });
  assert.ok(e.L > e0.L);
  assert.ok(rel(e.L / e0.L, wavenumber({ T: 9, d: 25 }) / wavenumber({ T: 9, d: 25, U: 1 })) < 5e-3);
  const z = streamFunctionWave({ H: 4, T: 9, d: 25, currentStokes: 0 });
  assert.ok(Math.abs(z.currentStokes) < 1e-11);
  assert.ok(z.currentEulerian < 0);
  // depth-integrated mean flux equals currentStokes · d: check Q = (c − cS) d
  assert.ok(rel(z.Q, z.c * 25) < 1e-11);
});

test('stream function, deep water: c²k/g = 1 + ε² + ε⁴/2 + ε⁶/4 (Stokes series) at ε = ak = 0.2', () => {
  const L = 100; const k = (2 * Math.PI) / L; const eps = 0.2;
  const s = streamFunctionWave({ H: (2 * eps) / k, L, d: 3 * L, order: 16 });
  assert.ok(s.report.dynamicResidual < 1e-9);
  const c2 = (s.c * s.c * k) / G;
  assert.ok(Math.abs(c2 - (1 + eps ** 2 + eps ** 4 / 2 + eps ** 6 / 4)) < 2e-5);
});

test('stream function, near-limiting deep-water wave: phase speed peaks near ak = 0.436 at about 1.093 sqrt(g/k) '
  + '(Longuet-Higgins 1975, Proc. R. Soc. A 342; Cokelet 1977, Phil. Trans. A 286)', () => {
  const L = 100; const k = (2 * Math.PI) / L;
  const c = (ak, order) => streamFunctionWave({ H: (2 * ak) / k, L, d: 3 * L, order }).c / Math.sqrt(G / k);
  const cPeak = c(0.436, 24);
  assert.ok(Math.abs(cPeak - 1.0929) < 5e-4, `c/c0 = ${cPeak}`);
  assert.ok(cPeak > c(0.42, 24) && c(0.42, 24) > c(0.40, 24));
  // crest is markedly higher than half the wave height for a steep wave
  const s = streamFunctionWave({ H: (2 * 0.42) / k, L, d: 3 * L, order: 24 });
  assert.ok(s.crest / s.H > 0.64 && s.crest / s.H < 0.66);
});

test('stream function: shallow-water steep wave converges with height stepping; invalid input is rejected', () => {
  const s = streamFunctionWave({ H: 6, T: 12, d: 10, order: 20 });
  assert.ok(s.report.heightSteps > 1);
  assert.ok(s.report.residual < 1e-10);
  assert.ok(s.crest / s.H > 0.75); // strongly peaked crest in shallow water
  assert.throws(() => streamFunctionWave({ H: 12, T: 12, d: 10 }), Error); // beyond breaking
  assert.throws(() => streamFunctionWave({ H: -1, T: 12, d: 10 }), RangeError);
  assert.throws(() => streamFunctionWave({ H: 1, T: 12, d: 10, order: 31 }), RangeError);
  assert.throws(() => streamFunctionWave({ H: 1, d: 10 }), RangeError);
});
