import test from 'node:test';
import assert from 'node:assert/strict';
import { G, dispersion, jonswap, spectrumMoments, csolve, boundedLeastSquares, prepWave, airyPoint, airyEta } from '../../src/engine/hydro/_wavehelp.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg || ''} got ${a}, expected ${b}`);

test('dispersion: satisfies omega² = g k tanh(kd); deep and shallow limits', () => {
  for (const [T, d] of [[10, 30], [6, 5], [15, 1000], [20, 3]]) {
    const w = (2 * Math.PI) / T;
    const k = dispersion(w, d);
    close(G * k * Math.tanh(k * d), w * w, 1e-12);
  }
  close(dispersion(1.0), 1 / G, 1e-14, 'deep water');
  // shallow water: c = sqrt(g d)
  const w = (2 * Math.PI) / 200;
  close(w / dispersion(w, 2), Math.sqrt(G * 2), 2e-3, 'shallow celerity');
});

test('JONSWAP: m0 = Hs²/16 (exact for gamma = 1, within 1 % for gamma = 3.3 by the DNV normalisation)', () => {
  const n = 4000;
  const om = Array.from({ length: n }, (_, i) => 0.05 + (i * 5.95) / (n - 1));
  for (const [gamma, tol] of [[1, 2e-3], [3.3, 1e-2]]) {
    const S = om.map((w) => jonswap(w, { Hs: 4, Tp: 10, gamma }));
    close(spectrumMoments(om, S).m0, 1.0, tol, `gamma ${gamma}`);
  }
});

test('Airy kinematics: deep-water orbital velocity a·omega at the surface, w = d(eta)/dt at z = 0', () => {
  const w = prepWave({ H: 2, T: 8, depth: Infinity });
  const o = {};
  airyPoint(w, 0, 0, 0, 0, o);
  close(o.ux, 1 * w.omega, 1e-12);
  const dt = 1e-5;
  airyPoint(w, 3, 0, 0, 1.3, o);
  close(o.uz, (airyEta(w, 3, 0, 1.3 + dt) - airyEta(w, 3, 0, 1.3 - dt)) / (2 * dt), 1e-7);
});

test('csolve: complex 3×3 system residual', () => {
  const Ar = [[4, 1, 0], [1, 3, -1], [0.5, 0, 2]]; const Ai = [[0.5, 0, 1], [0, -1, 0], [0.2, 0.3, 0]];
  const br = [1, 2, 3]; const bi = [0, -1, 0.5];
  const x = csolve(Ar.map((r) => r.slice()), Ai.map((r) => r.slice()), br.slice(), bi.slice());
  for (let i = 0; i < 3; i++) {
    let rr = -br[i]; let ri = -bi[i];
    for (let j = 0; j < 3; j++) { rr += Ar[i][j] * x.re[j] - Ai[i][j] * x.im[j]; ri += Ar[i][j] * x.im[j] + Ai[i][j] * x.re[j]; }
    assert.ok(Math.hypot(rr, ri) < 1e-12);
  }
});

test('boundedLeastSquares: hand-solved 2-variable problem with an active bound', () => {
  // minimise (x1 + x2 - 3)² + (x1 - x2)² with 0 <= x <= 1: unconstrained optimum (1.5, 1.5) -> clipped to (1, 1)
  const r = boundedLeastSquares([[1, 1], [1, -1]], [3, 0], [0, 0], [1, 1]);
  close(r.x[0], 1, 1e-9); close(r.x[1], 1, 1e-9);
  // inactive bounds reproduce the unconstrained solution
  const u = boundedLeastSquares([[1, 1], [1, -1]], [3, 1], [-10, -10], [10, 10]);
  close(u.x[0], 2, 1e-8); close(u.x[1], 1, 1e-8);
});
