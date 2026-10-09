import test from 'node:test';
import assert from 'node:assert/strict';
import { brent, bisect, newton, findBracket, newtonSystem, fixedPoint, polyRoots } from '../../src/engine/math/roots.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

const DOTTIE = 0.7390851332151607; // fixed point of cos x

test('brent: cos x = x (Dottie number 0.7390851332151607)', () => {
  close(brent((x) => Math.cos(x) - x, 0, 1), DOTTIE, 1e-12);
  close(brent((x) => Math.cos(x) - x, 1, 0), DOTTIE, 1e-12, 'reversed bracket');
});

test("brent: Wallis' cubic x³ − 2x − 5 = 0, root 2.0945514815423265", () => {
  let evals = 0;
  const r = brent((x) => {
    evals++;
    return x * x * x - 2 * x - 5;
  }, 2, 3);
  close(r, 2.0945514815423265, 1e-12);
  assert.ok(evals < 20, `super-linear convergence expected, used ${evals} evaluations`);
});

test('brent: root at an end point, flat multiple root, tolerance honoured, no sign change throws', () => {
  assert.equal(brent((x) => x - 1, 1, 2), 1);
  close(brent((x) => (x - 0.3) ** 3, -1, 2, { tol: 1e-10 }), 0.3, 1e-6);
  close(brent(Math.sin, 3, 4, { tol: 1e-3 }), Math.PI, 2e-3);
  assert.throws(() => brent((x) => x * x + 1, -1, 1), RangeError);
  assert.throws(() => brent(Math.sin, 0.1, Infinity), RangeError);
});

test('bisect: √2 and a discontinuous sign change', () => {
  close(bisect((x) => x * x - 2, 0, 2), Math.SQRT2, 1e-11);
  close(bisect((x) => (x < 0.3 ? -1 : 1), 0, 1, { tol: 1e-9 }), 0.3, 1e-8);
  assert.throws(() => bisect((x) => x * x + 1, -1, 1), RangeError);
});

test('newton: √2 with analytic and numerical derivative; Kepler equation E − e sin E = M', () => {
  close(newton((x) => x * x - 2, 1, { df: (x) => 2 * x }), Math.SQRT2, 1e-14);
  close(newton((x) => x * x - 2, 1), Math.SQRT2, 1e-12);
  const e = 0.6;
  const M = 1.2;
  const E = newton((x) => x - e * Math.sin(x) - M, M);
  close(E - e * Math.sin(E), M, 1e-12);
  assert.throws(() => newton((x) => x * x + 1, 0.5, { maxIter: 20 }), Error);
});

test('findBracket + brent: linear wave dispersion ω² = g k tanh(k h) for T = 10 s, h = 30 m', () => {
  const g = 9.80665;
  const h = 30;
  const w = (2 * Math.PI) / 10;
  const f = (k) => g * k * Math.tanh(k * h) - w * w;
  const [a, b] = findBracket(f, w * w / g, { min: 0 });
  assert.ok(a >= 0 && a < b);
  assert.ok(f(a) * f(b) <= 0);
  const k = brent(f, a, b);
  close(f(k), 0, 1e-12);
  // shallower than deep water: k > ω²/g, and the wavelength is between the shallow- and deep-water limits
  assert.ok(k > (w * w) / g);
  const L = (2 * Math.PI) / k;
  assert.ok(L < (g * 100) / (2 * Math.PI) && L > 10 * Math.sqrt(g * h) * 0.7);
});

test('findBracket: expands outwards, respects limits, throws when no sign change exists', () => {
  const [a, b] = findBracket((x) => x * x - 50, 1, { min: 0 });
  assert.ok(a >= 0 && (a * a - 50) * (b * b - 50) <= 0);
  const [c, d] = findBracket((x) => x + 1000, 0);
  assert.ok(c <= -1000 && d >= -1000);
  assert.throws(() => findBracket((x) => x * x + 1, 0, { min: -10, max: 10 }), Error);
  assert.throws(() => findBracket(Math.sin, 0, { factor: 1 }), RangeError);
});

test('newtonSystem: circle x² + y² = 4 intersected with y = x -> (√2, √2), finite-difference Jacobian', () => {
  const r = newtonSystem((x) => [x[0] * x[0] + x[1] * x[1] - 4, x[0] - x[1]], [1, 0.5]);
  assert.ok(r.converged);
  close(r.x[0], Math.SQRT2, 1e-9);
  close(r.x[1], Math.SQRT2, 1e-9);
  assert.ok(r.residual < 1e-9);
  assert.ok(r.iterations <= 8);
});

test('newtonSystem: analytic Jacobian converges quadratically; Rosenbrock stationarity from a far start', () => {
  const F = (x) => [x[0] * x[0] + x[1] * x[1] - 4, x[0] - x[1]];
  const jac = (x) => [[2 * x[0], 2 * x[1]], [1, -1]];
  const r = newtonSystem(F, [1, 0.5], { jac, tol: 1e-14 });
  assert.ok(r.converged);
  close(r.x[0], Math.SQRT2, 1e-14);
  // gradient of the Rosenbrock function vanishes at (1, 1); needs the line search from (−1.2, 1)
  const G = (x) => [-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0] * x[0]), 200 * (x[1] - x[0] * x[0])];
  const rr = newtonSystem(G, [-1.2, 1], { maxIter: 200 });
  assert.ok(rr.converged, `residual ${rr.residual} after ${rr.iterations} iterations`);
  close(rr.x[0], 1, 1e-7);
  close(rr.x[1], 1, 1e-7);
});

test('newtonSystem: catenary parameters from span and sag (3 unknowns), residual check', () => {
  // unknowns: a (catenary parameter), x1, x2 (end abscissae); ends 100 m apart horizontally, 20 m vertically,
  // unstretched length 120 m
  const F = (p) => {
    const [a, x1, x2] = p;
    return [
      x2 - x1 - 100,
      a * (Math.cosh(x2 / a) - Math.cosh(x1 / a)) - 20,
      a * (Math.sinh(x2 / a) - Math.sinh(x1 / a)) - 120,
    ];
  };
  const r = newtonSystem(F, [50, -40, 60]);
  assert.ok(r.converged);
  const res = F(r.x);
  for (const v of res) close(v, 0, 1e-8);
  // independent closed form: sqrt(L² − v²) = 2a sinh(H/(2a))
  const a = r.x[0];
  close(2 * a * Math.sinh(100 / (2 * a)), Math.sqrt(120 * 120 - 20 * 20), 1e-7);
});

test('newtonSystem: reports non-convergence instead of throwing', () => {
  const r = newtonSystem((x) => [x[0] * x[0] + 1], [3], { maxIter: 15 });
  assert.equal(r.converged, false);
  assert.ok(r.residual >= 1);
});

test('fixedPoint: x = cos x (scalar), 2-D contraction, under-relaxation of a divergent map', () => {
  const r = fixedPoint(Math.cos, 1, { tol: 1e-13, maxIter: 500 });
  assert.ok(r.converged);
  close(r.x, DOTTIE, 1e-11);
  // x = 0.5 cos y, y = 0.5 sin x + 1
  const v = fixedPoint((p) => [0.5 * Math.cos(p[1]), 0.5 * Math.sin(p[0]) + 1], [0, 0], { tol: 1e-13 });
  assert.ok(v.converged);
  close(v.x[0], 0.5 * Math.cos(v.x[1]), 1e-11);
  close(v.x[1], 0.5 * Math.sin(v.x[0]) + 1, 1e-11);
  // g(x) = 3 − 2x has |g'| = 2: plain iteration diverges, relax = 0.4 converges to x = 1
  const bad = fixedPoint((x) => 3 - 2 * x, 0, { maxIter: 30 });
  assert.equal(bad.converged, false);
  const good = fixedPoint((x) => 3 - 2 * x, 0, { relax: 0.4 });
  assert.ok(good.converged);
  close(good.x, 1, 1e-9);
});

function sortedPairs(r) {
  return r.re.map((v, i) => [v, r.im[i]]);
}

test('polyRoots: (x−1)(x−2)(x−3)(x−4)(x−5)(x−6) expanded', () => {
  const r = polyRoots([1, -21, 175, -735, 1624, -1764, 720]);
  for (let k = 0; k < 6; k++) {
    close(r.re[k], k + 1, 1e-9);
    assert.equal(r.im[k], 0);
  }
});

test('polyRoots: complex roots — x² + 1, x⁴ + 4 (±1 ± i), x³ − 1 (cube roots of unity)', () => {
  const a = polyRoots([1, 0, 1]);
  assert.deepEqual(a.re.map((v) => Math.abs(v) < 1e-15 ? 0 : v), [0, 0]);
  close(a.im[0], -1, 1e-14);
  close(a.im[1], 1, 1e-14);
  const b = sortedPairs(polyRoots([1, 0, 0, 0, 4]));
  const expected = [[-1, -1], [-1, 1], [1, -1], [1, 1]];
  for (let k = 0; k < 4; k++) {
    close(b[k][0], expected[k][0], 1e-13);
    close(b[k][1], expected[k][1], 1e-13);
  }
  const c = sortedPairs(polyRoots([1, 0, 0, -1]));
  close(c[0][0], -0.5, 1e-14);
  close(c[0][1], -Math.sqrt(3) / 2, 1e-14);
  close(c[1][1], Math.sqrt(3) / 2, 1e-14);
  close(c[2][0], 1, 1e-14);
  assert.equal(c[2][1], 0);
});

test('polyRoots: Chebyshev T8 roots are cos((2k−1)π/16)', () => {
  // T8(x) = 128x⁸ − 256x⁶ + 160x⁴ − 32x² + 1
  const r = polyRoots([128, 0, -256, 0, 160, 0, -32, 0, 1]);
  const exact = [];
  for (let k = 1; k <= 8; k++) exact.push(Math.cos(((2 * k - 1) * Math.PI) / 16));
  exact.sort((p, q) => p - q);
  for (let k = 0; k < 8; k++) {
    close(r.re[k], exact[k], 1e-12);
    assert.equal(r.im[k], 0);
  }
});

test('polyRoots: double root, zero roots, leading zeros, widely scaled coefficients, linear', () => {
  // (x − 1)²(x + 2) = x³ − 3x + 2: a double root is resolved to ~sqrt(eps)
  const d = polyRoots([1, 0, -3, 2]);
  close(d.re[0], -2, 1e-10);
  close(d.re[1], 1, 1e-6);
  close(d.re[2], 1, 1e-6);
  assert.deepEqual(d.im, [0, 0, 0]);
  // x²(x − 3)
  const z = polyRoots([1, -3, 0, 0]);
  assert.deepEqual(z.re.slice(0, 2), [0, 0]);
  close(z.re[2], 3, 1e-14);
  // leading zero coefficients are ignored
  const l = polyRoots([0, 0, 2, -8]);
  assert.deepEqual(l.re, [4]);
  // roots 1e-6 and 1e6: x² − (1e6 + 1e-6)x + 1
  const w = polyRoots([1, -(1e6 + 1e-6), 1]);
  close(w.re[0] / 1e-6, 1, 1e-9);
  close(w.re[1] / 1e6, 1, 1e-12);
  assert.throws(() => polyRoots([0, 0, 5]), RangeError);
  assert.throws(() => polyRoots([]), RangeError);
});

test('polyRoots: characteristic polynomial of a damped oscillator s² + 2ζω s + ω² (ζ = 0.05, ω = 3)', () => {
  const zeta = 0.05;
  const w = 3;
  const r = polyRoots([1, 2 * zeta * w, w * w]);
  close(r.re[0], -zeta * w, 1e-14);
  close(r.re[1], -zeta * w, 1e-14);
  close(Math.abs(r.im[0]), w * Math.sqrt(1 - zeta * zeta), 1e-14);
  close(r.im[0] + r.im[1], 0, 1e-14);
});
