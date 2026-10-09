import test from 'node:test';
import assert from 'node:assert/strict';
import {
  cx, cadd, csub, cmul, cdiv, cabs, carg, cpolar, csqrt, cexp, clog, csinh, ccosh, ctanh, cinv, cconj, cparallel,
  csolve, abcdMul,
} from '../../src/engine/electrical/complex.js';

const near = (a, b, tol = 1e-12) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${a} vs ${b}`);
const cnear = (z, re, im, tol = 1e-12) => { near(z.re, re, tol); near(z.im, im, tol); };

test('arithmetic by hand: (1+2j)(3−j) = 5+5j, (1+2j)/(3−j) = 0.1+0.7j', () => {
  const a = cx(1, 2); const b = cx(3, -1);
  cnear(cmul(a, b), 5, 5);
  cnear(cdiv(a, b), 0.1, 0.7);
  cnear(cadd(a, b), 4, 1);
  cnear(csub(a, b), -2, 3);
  cnear(cmul(cinv(b), b), 1, 0);
  cnear(cconj(a), 1, -2);
});

test('polar form: |3+4j| = 5, arg(j) = π/2, polar round trip', () => {
  near(cabs(cx(3, 4)), 5);
  near(carg(cx(0, 1)), Math.PI / 2);
  cnear(cpolar(2, Math.PI / 3), 1, Math.sqrt(3));
});

test('principal square root: sqrt(−4) = 2j, sqrt(−4 − 0j·ε) branch, sqrt(3+4j) = 2+j', () => {
  cnear(csqrt(cx(-4, 0)), 0, 2);
  cnear(csqrt(cx(3, 4)), 2, 1);
  cnear(csqrt(cx(-3, -4)), 1, -2);
  const z = cx(-0.3, 0.7); const r = csqrt(z);
  cnear(cmul(r, r), z.re, z.im);
  assert.ok(r.re >= 0);
});

test('Euler identity e^(jπ) = −1 and log/exp inverse', () => {
  cnear(cexp(cx(0, Math.PI)), -1, 0, 1e-15);
  const z = cx(0.4, -1.1);
  cnear(clog(cexp(z)), z.re, z.im);
});

test('hyperbolic identities: cosh² − sinh² = 1, cosh(jx) = cos x, sinh(jx) = j sin x, tanh = sinh/cosh', () => {
  const z = cx(0.7, 1.9);
  const c = ccosh(z); const s = csinh(z);
  cnear(csub(cmul(c, c), cmul(s, s)), 1, 0);
  cnear(ccosh(cx(0, 0.3)), Math.cos(0.3), 0);
  cnear(csinh(cx(0, 0.3)), 0, Math.sin(0.3));
  const t = ctanh(z); const q = cdiv(s, c);
  cnear(t, q.re, q.im);
});

test('parallel impedances: 4 Ω ∥ 4j Ω = 2+2j Ω', () => {
  cnear(cparallel(cx(4, 0), cx(0, 4)), 2, 2);
});

test('complex linear solve, 2×2 by hand: [[1, j], [−j, 2]] x = [1, 0] → x = [2, j]', () => {
  const x = csolve([[cx(1, 0), cx(0, 1)], [cx(0, -1), cx(2, 0)]], [cx(1, 0), cx(0, 0)]);
  cnear(x[0], 2, 0); cnear(x[1], 0, 1);
});

test('complex linear solve needs pivoting and reproduces the right-hand side (4×4 residual)', () => {
  const A = [
    [cx(0, 0), cx(2, 1), cx(1, 0), cx(0, -1)],
    [cx(1, 1), cx(0, 0), cx(3, -2), cx(1, 0)],
    [cx(2, 0), cx(1, -1), cx(0, 0), cx(4, 1)],
    [cx(0, 3), cx(1, 0), cx(1, 1), cx(0, 0)],
  ];
  const b = [cx(1, 0), cx(0, 2), cx(-1, 1), cx(3, -3)];
  const x = csolve(A, b);
  for (let i = 0; i < 4; i++) {
    let s = cx(0, 0);
    for (let j = 0; j < 4; j++) s = cadd(s, cmul(A[i][j], x[j]));
    cnear(s, b[i].re, b[i].im, 1e-11);
  }
  assert.throws(() => csolve([[cx(1, 0), cx(2, 0)], [cx(2, 0), cx(4, 0)]], [cx(1, 0), cx(1, 0)]), RangeError);
  assert.throws(() => cdiv(cx(1, 0), cx(0, 0)), RangeError);
});

test('ABCD product: two series impedances add, determinant stays 1', () => {
  const s = (z) => ({ A: cx(1, 0), B: z, C: cx(0, 0), D: cx(1, 0) });
  const m = abcdMul(s(cx(1, 2)), s(cx(3, -1)));
  cnear(m.B, 4, 1);
  const y = { A: cx(1, 0), B: cx(0, 0), C: cx(0, 0.01), D: cx(1, 0) };
  const p = abcdMul(abcdMul(s(cx(1, 2)), y), s(cx(0.5, 0.2)));
  cnear(csub(cmul(p.A, p.D), cmul(p.B, p.C)), 1, 0);
});
