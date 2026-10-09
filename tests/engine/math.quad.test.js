import test from 'node:test';
import assert from 'node:assert/strict';
import { trapz, cumtrapz, simpson, gaussLegendre, integrate } from '../../src/engine/math/quad.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

function sample(f, a, b, n) {
  const x = new Float64Array(n + 1);
  const y = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    x[i] = a + ((b - a) * i) / n;
    y[i] = f(x[i]);
  }
  return { x, y, dx: (b - a) / n };
}

test('trapz: exact for linear data, uniform and non-uniform abscissae', () => {
  // ∫₀⁴ (2x + 1) dx = 20
  assert.equal(trapz([1, 3, 5, 7, 9], 1), 20);
  assert.equal(trapz([1, 3, 9], [0, 1, 4]), 20);
  assert.equal(trapz([5], 1), 0);
  assert.equal(trapz([1, 2, 3]), 4, 'default spacing 1');
  assert.throws(() => trapz([1, 2], [0, 1, 2]), RangeError);
});

test('trapz: second-order convergence on ∫₀^π sin x dx = 2 (error ratio 4 per halving)', () => {
  const errs = [16, 32, 64, 128].map((n) => {
    const s = sample(Math.sin, 0, Math.PI, n);
    return Math.abs(trapz(s.y, s.dx) - 2);
  });
  for (let k = 1; k < errs.length; k++) {
    const order = Math.log2(errs[k - 1] / errs[k]);
    close(order, 2, 0.02, 'observed order');
  }
  // leading error term −(b−a)h²/12·f″: for sin on [0,π] the error is −h²/6 (Euler–Maclaurin)
  const h = Math.PI / 128;
  close(errs[3], (h * h) / 6, 1e-7);
});

test('simpson: exact for cubics with even and odd numbers of intervals', () => {
  const f = (x) => x * x * x - 2 * x * x + 3 * x - 1;
  const exact = 5 ** 4 / 4 - (2 * 5 ** 3) / 3 + (3 * 25) / 2 - 5; // ∫₀⁵
  for (const n of [2, 3, 4, 5, 6, 7, 10, 11]) {
    const s = sample(f, 0, 5, n);
    close(simpson(s.y, s.dx), exact, 1e-12, `n=${n}`);
  }
  assert.equal(simpson([1, 3], 2), 4, 'two samples fall back to the trapezoid');
  assert.equal(simpson([7], 1), 0);
});

test('simpson: fourth-order convergence on ∫₀^π sin x dx = 2, even and odd interval counts', () => {
  for (const base of [16, 17]) {
    const errs = [base, 2 * base, 4 * base, 8 * base].map((n) => {
      const s = sample(Math.sin, 0, Math.PI, n);
      return Math.abs(simpson(s.y, s.dx) - 2);
    });
    for (let k = 1; k < errs.length; k++) {
      const order = Math.log2(errs[k - 1] / errs[k]);
      assert.ok(order > 3.7 && order < 4.5, `base ${base}: observed order ${order}`);
    }
  }
  // composite Simpson error ≈ (b−a) h⁴ f⁗/180: for sin on [0, π] -> h⁴/90
  const s = sample(Math.sin, 0, Math.PI, 64);
  const h = Math.PI / 64;
  close(Math.abs(simpson(s.y, s.dx) - 2), h ** 4 / 90, 1e-9);
});

test('simpson: non-uniform abscissae are exact for quadratics (even and odd interval counts)', () => {
  const f = (x) => 3 * x * x - x + 2;
  const F = (x) => x * x * x - (x * x) / 2 + 2 * x;
  for (const x of [[0, 0.5, 2, 2.2, 4], [0, 0.5, 2, 2.2, 4, 4.3], [1, 3, 4]]) {
    const y = x.map(f);
    close(simpson(y, x), F(x[x.length - 1]) - F(x[0]), 1e-12);
  }
});

test('cumtrapz: running integral, last value equals trapz, non-uniform grid', () => {
  const c = cumtrapz([1, 3, 5, 7], 1);
  assert.deepEqual(Array.from(c), [0, 2, 6, 12]);
  const x = [0, 1, 3, 6];
  const y = [2, 2, 4, 4];
  const cx = cumtrapz(y, x);
  assert.deepEqual(Array.from(cx), [0, 2, 8, 20]);
  assert.equal(cx[3], trapz(y, x));
  // ∫₀ˣ cos = sin x
  const s = sample(Math.cos, 0, 2, 2000);
  const ci = cumtrapz(s.y, s.dx);
  close(ci[1000], Math.sin(1), 1e-7);
});

test('gaussLegendre: tabulated nodes and weights for n = 1, 2, 3 (A&S Table 25.4)', () => {
  const g1 = gaussLegendre(1);
  assert.deepEqual(g1, { x: [0], w: [2] });
  const g2 = gaussLegendre(2);
  close(g2.x[0], -1 / Math.sqrt(3), 1e-15);
  close(g2.x[1], 1 / Math.sqrt(3), 1e-15);
  close(g2.w[0], 1, 1e-15);
  const g3 = gaussLegendre(3);
  close(g3.x[0], -Math.sqrt(0.6), 1e-15);
  assert.equal(g3.x[1], 0);
  close(g3.w[0], 5 / 9, 1e-15);
  close(g3.w[1], 8 / 9, 1e-15);
  // n = 5 outer node 0.9061798459386640, weight 0.2369268850561891
  const g5 = gaussLegendre(5);
  close(g5.x[4], 0.906179845938664, 1e-15);
  close(g5.w[4], 0.2369268850561891, 1e-15);
});

test('gaussLegendre: exact to degree 2n − 1 and not beyond, weights sum to 2, n up to 64', () => {
  for (const n of [4, 10, 20, 33, 64]) {
    const { x, w } = gaussLegendre(n);
    assert.equal(x.length, n);
    let sw = 0;
    for (const v of w) sw += v;
    close(sw, 2, 1e-13, `n=${n} weight sum`);
    for (let i = 1; i < n; i++) assert.ok(x[i] > x[i - 1]);
    // ∫ x^(2n−2) dx = 2/(2n−1) exactly integrated
    let s = 0;
    for (let i = 0; i < n; i++) s += w[i] * x[i] ** (2 * n - 2);
    close(s * (2 * n - 1), 2, 1e-12, `n=${n} degree 2n−2`);
    // odd monomial integrates to zero by symmetry
    let o = 0;
    for (let i = 0; i < n; i++) o += w[i] * x[i] ** 3;
    close(o, 0, 1e-15);
  }
  // degree 2n is not exact: for n = 2, ∫x⁴ = 2/5 but the rule gives 2/9
  const g = gaussLegendre(2);
  close(g.w[0] * g.x[0] ** 4 + g.w[1] * g.x[1] ** 4, 2 / 9, 1e-15);
  assert.throws(() => gaussLegendre(0), RangeError);
  assert.throws(() => gaussLegendre(2.5), RangeError);
});

test('gaussLegendre: returned arrays are independent copies', () => {
  const a = gaussLegendre(6);
  a.x[0] = 99;
  assert.notEqual(gaussLegendre(6).x[0], 99);
});

test('integrate: smooth integrands to tolerance — ∫₀^π sin = 2, ∫₀¹ 4/(1+x²) = π, ∫₀¹ eˣ = e − 1', () => {
  close(integrate(Math.sin, 0, Math.PI), 2, 1e-12);
  close(integrate((x) => 4 / (1 + x * x), 0, 1), Math.PI, 1e-12);
  close(integrate(Math.exp, 0, 1), Math.E - 1, 1e-12);
  close(integrate(Math.sin, Math.PI, 0), -2, 1e-12, 'reversed limits change the sign');
  assert.equal(integrate(Math.sin, 1, 1), 0);
});

test('integrate: one 15-point panel suffices when both rules are exact (degree 12); degree 22 integrated to round-off', () => {
  let calls = 0;
  const v = integrate((x) => {
    calls++;
    return x ** 12;
  }, -1, 1);
  close(v, 2 / 13, 1e-15);
  assert.equal(calls, 15, 'Gauss(7) is exact to degree 13, so the error estimate vanishes');
  close(integrate((x) => x ** 22, -1, 1), 2 / 23, 1e-15);
});

test('integrate: sharp Lorentzian peak 1/((x−0.3)² + 10⁻⁴) on [0, 1] vs arctan closed form', () => {
  const eps = 0.01;
  const exact = (Math.atan((1 - 0.3) / eps) - Math.atan((0 - 0.3) / eps)) / eps;
  const v = integrate((x) => 1 / ((x - 0.3) ** 2 + eps * eps), 0, 1);
  close(v / exact, 1, 1e-9);
});

test('integrate: oscillatory ∫₀¹⁰ cos(20x) dx = sin(200)/20 with heavy cancellation', () => {
  close(integrate((x) => Math.cos(20 * x), 0, 10, { tol: 1e-10 }), Math.sin(200) / 20, 1e-9);
});

test('integrate: infinite ranges — Gaussian, Γ(4) = ∫₀^∞ x³e⁻ˣ = 6, ∫_{−∞}^0 eˣ = 1', () => {
  close(integrate((x) => Math.exp(-x * x), -Infinity, Infinity), Math.sqrt(Math.PI), 1e-9);
  close(integrate((x) => x ** 3 * Math.exp(-x), 0, Infinity), 6, 1e-8);
  close(integrate(Math.exp, -Infinity, 0), 1, 1e-9);
  close(integrate((x) => 1 / (x * x), 1, Infinity), 1, 1e-9);
});

test('integrate: integrable end-point singularities — ∫₀¹ ln x = −1, ∫₀¹ x^(−1/2) = 2', () => {
  close(integrate(Math.log, 0, 1, { tol: 1e-10 }), -1, 1e-8);
  close(integrate((x) => 1 / Math.sqrt(x), 0, 1, { maxDepth: 80 }), 2, 1e-7);
});

test('integrate: nested (re-entrant) double integral ∫₀¹∫₀ˣ xy dy dx = 1/8', () => {
  const v = integrate((x) => integrate((y) => x * y, 0, x), 0, 1);
  close(v, 1 / 8, 1e-12);
  // area of the unit quarter disc, inner limit with a square-root end point
  const a = integrate((x) => integrate(() => 1, 0, Math.sqrt(1 - x * x)), 0, 1, { tol: 1e-10 });
  close(a, Math.PI / 4, 1e-8);
});

test('integrate: tighter tolerance gives a smaller error; invalid arguments throw', () => {
  const f = (x) => Math.sqrt(x) * Math.cos(3 * x);
  const ref = integrate(f, 0, 2, { tol: 1e-13, maxDepth: 60 });
  const loose = Math.abs(integrate(f, 0, 2, { tol: 1e-3 }) - ref);
  const tight = Math.abs(integrate(f, 0, 2, { tol: 1e-9 }) - ref);
  assert.ok(loose <= 1e-3 * Math.abs(ref) + 1e-12);
  assert.ok(tight <= 1e-9);
  assert.throws(() => integrate(5, 0, 1), RangeError);
  assert.throws(() => integrate(Math.sin, NaN, 1), RangeError);
  assert.throws(() => integrate(Math.sin, 0, 1, { tol: 0 }), RangeError);
});
