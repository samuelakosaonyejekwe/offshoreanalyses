import test from 'node:test';
import assert from 'node:assert/strict';
import {
  erf, erfc, erfinv, gamma, lgamma, gammaP, gammaQ, betaI, besselJ0, besselJ1, besselY0, besselY1, besselJn,
  besselI0, besselI1, besselK0, besselK1,
} from '../../src/engine/math/special.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

function rel(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg} got ${a}, expected ${b} (rel tol ${tol})`);
}

test('erf: tabulated values (A&S Table 7.1) — erf(0.5) = 0.5204998778, erf(1) = 0.8427007929, erf(2) = 0.9953222650', () => {
  close(erf(0.5), 0.5204998778, 1e-10);
  close(erf(0.5), 0.5204998778130465, 1e-15);
  close(erf(1), 0.8427007929497149, 1e-15);
  close(erf(2), 0.9953222650189527, 1e-15);
  close(erf(0.1), 0.1124629160182849, 1e-15);
  close(erf(3), 0.9999779095030014, 1e-15);
  assert.equal(erf(0), 0);
  assert.equal(erf(Infinity), 1);
  assert.equal(erf(-Infinity), -1);
  assert.equal(erf(-0.7), -erf(0.7), 'odd symmetry');
  rel(erf(1e-10), (2 / Math.sqrt(Math.PI)) * 1e-10, 1e-15, 'small-argument limit 2x/√π');
  assert.throws(() => erf(NaN), RangeError);
});

test('erfc: far tail keeps relative accuracy — erfc(3), erfc(5), erfc(10) = 2.088487583762545e-45', () => {
  rel(erfc(0.5), 0.4795001221869535, 1e-15);
  rel(erfc(3), 2.209049699858544e-5, 1e-14);
  rel(erfc(5), 1.537459794428035e-12, 1e-14);
  rel(erfc(10), 2.088487583762545e-45, 1e-13);
  assert.equal(erfc(0), 1);
  assert.equal(erfc(40), 0);
  close(erfc(-1), 2 - erfc(1), 1e-15, 'erfc(−x) = 2 − erfc(x)');
  assert.equal(erfc(-30), 2);
});

test('erf + erfc = 1 and continuity across the internal switch at |x| = 0.5', () => {
  for (const x of [-3, -0.5, -0.2, 0.3, 0.4999999, 0.5, 0.5000001, 1.7, 4]) close(erf(x) + erfc(x), 1, 2e-16 * 2);
  close(erf(0.5 + 1e-9) - erf(0.5 - 1e-9), (2e-9 * 2 * Math.exp(-0.25)) / Math.sqrt(Math.PI), 1e-15);
});

test('erfc agrees with the independent continued fraction Q(1/2, x²) over 0.5 <= x <= 25', () => {
  for (const x of [0.5, 0.8, 1.3, 2, 3.5, 6, 6.6, 9, 15, 25]) rel(erfc(x), gammaQ(0.5, x * x), 1e-12, `x=${x}`);
});

test('erfinv: erfinv(0.5) = 0.4769362762044699, round trips incl. the tails, end points', () => {
  close(erfinv(0.5), 0.4769362762044699, 1e-15);
  for (const p of [-0.999999, -0.9, -0.3, 1e-12, 0.1, 0.7, 0.95, 0.999, 1 - 1e-10]) {
    const y = erfinv(p);
    if (Math.abs(p) > 0.7) rel(erfc(Math.abs(y)), 1 - Math.abs(p), 1e-12, `p=${p}`);
    else close(erf(y), p, 1e-15, `p=${p}`);
  }
  assert.equal(erfinv(0), 0);
  assert.equal(erfinv(1), Infinity);
  assert.equal(erfinv(-1), -Infinity);
  assert.equal(erfinv(-0.4), -erfinv(0.4));
  // standard normal 97.5 % quantile: √2 erfinv(0.95) = 1.959963984540054
  close(Math.SQRT2 * erfinv(0.95), 1.959963984540054, 1e-14);
  assert.throws(() => erfinv(1.0001), RangeError);
});

test('gamma: Γ(1/2) = √π, factorials, Γ(1/3) = 2.678938534707747, Γ(−1/2) = −2√π, Γ(−3/2) = 4√π/3', () => {
  rel(gamma(0.5), Math.sqrt(Math.PI), 1e-15);
  assert.equal(gamma(1), 1);
  assert.equal(gamma(5), 24);
  assert.equal(gamma(11), 3628800);
  rel(gamma(1 / 3), 2.678938534707747, 1e-14);
  rel(gamma(1.5), Math.sqrt(Math.PI) / 2, 1e-14);
  rel(gamma(4.5), 11.631728396567448, 1e-14);
  rel(gamma(-0.5), -2 * Math.sqrt(Math.PI), 1e-14);
  rel(gamma(-1.5), (4 * Math.sqrt(Math.PI)) / 3, 1e-14);
  rel(gamma(171), 7.257415615307994e306, 1e-14);
  assert.equal(gamma(172), Infinity);
  assert.throws(() => gamma(0), RangeError);
  assert.throws(() => gamma(-3), RangeError);
});

test('gamma: recurrence Γ(x+1) = x Γ(x) and reflection Γ(x)Γ(1−x) = π/sin(πx) on a grid', () => {
  for (let x = 0.05; x < 60; x += 0.731) rel(gamma(x + 1), x * gamma(x), 4e-14, `x=${x}`);
  for (const x of [0.1, 0.25, 0.37, 0.8]) rel(gamma(x) * gamma(1 - x), Math.PI / Math.sin(Math.PI * x), 1e-13);
  // Legendre duplication: Γ(x)Γ(x+½) = 2^(1−2x) √π Γ(2x)
  for (const x of [0.3, 1.7, 6.2]) rel(gamma(x) * gamma(x + 0.5), 2 ** (1 - 2 * x) * Math.sqrt(Math.PI) * gamma(2 * x), 1e-13);
});

test('lgamma: ln Γ(100) = ln 99! = 359.1342053695754, ln Γ(1/2) = ln √π, Stirling limit, negative arguments', () => {
  rel(lgamma(100), 359.1342053695754, 1e-15);
  rel(lgamma(0.5), 0.5 * Math.log(Math.PI), 1e-14);
  assert.equal(lgamma(1), 0);
  assert.equal(lgamma(2), 0);
  close(lgamma(3), Math.log(2), 1e-15);
  // Stirling series at x = 1e6: (x−½)ln x − x + ½ln 2π + 1/(12x)
  const x = 1e6;
  rel(lgamma(x), (x - 0.5) * Math.log(x) - x + 0.5 * Math.log(2 * Math.PI) + 1 / (12 * x), 1e-15);
  rel(lgamma(-0.5), Math.log(2 * Math.sqrt(Math.PI)), 1e-14);
  for (const y of [0.2, 3.3, 17.5, 150.5]) rel(Math.exp(lgamma(y) - lgamma(y + 1)), 1 / y, 1e-13);
  assert.throws(() => lgamma(-2), RangeError);
});

test('gammaP/gammaQ: closed forms P(1,x) = 1 − e⁻ˣ, P(2,x) = 1 − (1+x)e⁻ˣ, P(1/2,x) = erf √x', () => {
  for (const x of [0.01, 0.5, 1, 1.9, 2.1, 5, 30]) {
    rel(gammaP(1, x), 1 - Math.exp(-x), 1e-13, `P(1, ${x})`);
    rel(gammaQ(1, x), Math.exp(-x), 1e-12, `Q(1, ${x})`);
    // absolute tolerance: the closed form itself cancels for small x
    close(gammaP(2, x), 1 - (1 + x) * Math.exp(-x), 1e-15, `P(2, ${x})`);
    rel(gammaP(0.5, x), erf(Math.sqrt(x)), 1e-13, `P(1/2, ${x})`);
    close(gammaP(3.7, x) + gammaQ(3.7, x), 1, 1e-14);
  }
  assert.equal(gammaP(2, 0), 0);
  assert.equal(gammaQ(2, 0), 1);
  assert.equal(gammaP(2, Infinity), 1);
  assert.throws(() => gammaP(0, 1), RangeError);
  assert.throws(() => gammaQ(1, -1), RangeError);
});

test('gammaP: chi-square 95 % critical values (A&S Table 26.8) — ν = 2: 5.99146, ν = 10: 18.3070, ν = 30: 43.7730', () => {
  close(gammaP(1, 5.991464547107979 / 2), 0.95, 1e-13);
  close(gammaP(5, 18.307038053275146 / 2), 0.95, 1e-12);
  close(gammaP(15, 43.77297182574219 / 2), 0.95, 1e-12);
  // median of the exponential distribution
  close(gammaP(1, Math.LN2), 0.5, 1e-15);
});

test('betaI: closed forms I_x(1,1) = x, I_x(a,1) = xᵃ, I_x(2,3) = 6x² − 8x³ + 3x⁴, arcsine law', () => {
  for (const x of [0.05, 0.3, 0.4, 0.5, 0.77, 0.99]) {
    close(betaI(1, 1, x), x, 1e-14);
    rel(betaI(3.5, 1, x), x ** 3.5, 1e-13);
    close(betaI(2, 3, x), 6 * x * x - 8 * x ** 3 + 3 * x ** 4, 1e-13);
    close(betaI(0.5, 0.5, x), (2 / Math.PI) * Math.asin(Math.sqrt(x)), 1e-13);
    close(betaI(2.2, 7.1, x), 1 - betaI(7.1, 2.2, 1 - x), 1e-13, 'symmetry');
  }
  close(betaI(2, 3, 0.4), 0.5248, 1e-13);
  close(betaI(6.3, 6.3, 0.5), 0.5, 1e-13);
  assert.equal(betaI(2, 3, 0), 0);
  assert.equal(betaI(2, 3, 1), 1);
  assert.throws(() => betaI(2, 3, 1.1), RangeError);
  assert.throws(() => betaI(-1, 3, 0.5), RangeError);
});

test("betaI: Student's t two-sided 5 % point for ν = 10, t = 2.228138852 (A&S Table 26.10)", () => {
  const nu = 10;
  const t = 2.228138852; // tabulated to 10 significant figures
  close(betaI(nu / 2, 0.5, nu / (nu + t * t)), 0.05, 1e-10);
  // binomial tail: P(X >= 3 | n = 10, p = 0.2) = I_0.2(3, 8) = 1 − Σ_{k<=2} C(10,k) 0.2^k 0.8^(10−k)
  const tail = 1 - (0.8 ** 10 + 10 * 0.2 * 0.8 ** 9 + 45 * 0.04 * 0.8 ** 8);
  close(betaI(3, 8, 0.2), tail, 1e-14);
});

test('Bessel J0, J1, Y0, Y1: tabulated values at x = 1, 10, 50 (A&S Tables 9.1, 9.2)', () => {
  close(besselJ0(1), 0.7651976865579666, 1e-15);
  close(besselJ1(1), 0.4400505857449335, 1e-15);
  close(besselY0(1), 0.08825696421567696, 1e-15);
  close(besselY1(1), -0.7812128213002887, 1e-15);
  close(besselJ0(10), -0.2459357644513483, 1e-15);
  close(besselJ1(10), 0.04347274616886144, 1e-15);
  close(besselY0(10), 0.05567116728359939, 1e-15);
  close(besselY1(10), 0.24901542420695388, 1e-15);
  close(besselJ0(50), 0.05581232766925181, 1e-15);
  close(besselJ1(50), -0.09751182812517514, 1e-15);
  close(besselY0(50), -0.09806499547007708, 1e-15);
  close(besselY1(50), -0.05679566856201477, 1e-15);
  assert.equal(besselJ0(0), 1);
  assert.equal(besselJ1(0), 0);
  assert.equal(besselJ0(-3.3), besselJ0(3.3));
  assert.equal(besselJ1(-3.3), -besselJ1(3.3));
  assert.throws(() => besselY0(0), RangeError);
  assert.throws(() => besselY1(-1), RangeError);
});

test('Bessel zeros (A&S Table 9.5): j0,1 = 2.404825557695773, j0,2 = 5.520078110286311, j1,1 = 3.831705970207512, y0,1 = 0.893576966279168', () => {
  close(besselJ0(2.404825557695773), 0, 1e-15);
  close(besselJ0(5.520078110286311), 0, 1e-15);
  close(besselJ0(8.653727912911013), 0, 1e-15);
  close(besselJ1(3.8317059702075125), 0, 1e-15);
  close(besselJ1(7.015586669815619), 0, 1e-15);
  close(besselY0(0.8935769662791675), 0, 1e-15);
  close(besselY1(2.197141326031017), 0, 1e-15);
  // large-order zero in the asymptotic regime: j0,10 = 30.63460646843198
  close(besselJ0(30.634606468431976), 0, 1e-15);
});

test('Bessel Wronskian J1 Y0 − J0 Y1 = 2/(πx) across the series / recurrence / asymptotic regimes', () => {
  for (const x of [1e-3, 0.01, 0.5, 1.5, 1.999999, 2, 2.000001, 3.7, 8, 12.3, 20, 24.999999, 25, 25.000001, 40, 100, 1000, 12345.6]) {
    const w = besselJ1(x) * besselY0(x) - besselJ0(x) * besselY1(x);
    rel(w, 2 / (Math.PI * x), 4e-15 * Math.max(1, x), `x=${x}`);
  }
});

test('Bessel functions are continuous across the internal switches at x = 2 and x = 25', () => {
  for (const x0 of [2, 25]) {
    const h = 1e-7;
    for (const [f, df] of [
      [besselJ0, (x) => -besselJ1(x)],
      [besselY0, (x) => -besselY1(x)],
      [besselJ1, (x) => besselJ0(x) - besselJ1(x) / x],
      [besselY1, (x) => besselY0(x) - besselY1(x) / x],
    ]) {
      const jump = f(x0 + h) - f(x0 - h);
      close(jump, 2 * h * df(x0), 1e-13, `x0=${x0}`);
    }
  }
});

test('Bessel small-argument limits: J0 ≈ 1 − x²/4, Y0 ≈ (2/π)(ln(x/2) + γ), Y1 ≈ −2/(πx)', () => {
  const x = 1e-6;
  const g = 0.5772156649015329;
  close(besselJ0(x), 1 - (x * x) / 4, 1e-16);
  rel(besselJ1(x), x / 2, 1e-12);
  rel(besselY0(x), (2 / Math.PI) * (Math.log(x / 2) + g), 1e-12);
  rel(besselY1(x), -2 / (Math.PI * x), 1e-11);
});

test('besselJn: tabulated J2(1), J5(10), J10(1), three-term recurrence, parity, Σ identity', () => {
  close(besselJn(2, 1), 0.1149034849319005, 1e-15);
  close(besselJn(5, 10), -0.2340615281867936, 1e-15);
  rel(besselJn(10, 1), 2.630615123687453e-10, 1e-13);
  close(besselJn(3, 2.5), 0.2166003910391135, 1e-15);
  assert.equal(besselJn(0, 1.3), besselJ0(1.3));
  assert.equal(besselJn(1, 1.3), besselJ1(1.3));
  for (const x of [0.3, 2, 7.7, 30, 120]) {
    for (const n of [1, 2, 5, 9, 40]) {
      const lhs = besselJn(n - 1, x) + besselJn(n + 1, x);
      close(lhs, ((2 * n) / x) * besselJn(n, x), 1e-13, `n=${n}, x=${x}`);
    }
    // 1 = J0 + 2 Σ J_2k (A&S 9.1.46)
    let s = besselJ0(x);
    for (let k = 1; k < 120; k++) s += 2 * besselJn(2 * k, x);
    close(s, 1, 1e-13, `sum identity at x=${x}`);
  }
  assert.equal(besselJn(3, -2), -besselJn(3, 2));
  assert.equal(besselJn(-3, 2), -besselJn(3, 2));
  assert.equal(besselJn(-4, 2), besselJn(4, 2));
  assert.equal(besselJn(4, 0), 0);
  // leading term for small x: (x/2)^n / n!
  rel(besselJn(6, 1e-3), (5e-4) ** 6 / 720, 1e-6);
  assert.throws(() => besselJn(1.5, 1), RangeError);
});

test('modified Bessel I0, I1, K0, K1: tabulated values at x = 1, 2, 10 (A&S Table 9.8)', () => {
  rel(besselI0(1), 1.2660658777520084, 1e-15);
  rel(besselI1(1), 0.5651591039924851, 1e-15);
  rel(besselK0(1), 0.42102443824070834, 1e-15);
  rel(besselK1(1), 0.6019072301972346, 1e-15);
  rel(besselK0(2), 0.11389387274953344, 1e-15);
  rel(besselK1(2), 0.13986588181652243, 1e-15);
  rel(besselI0(10), 2815.716628466254, 1e-15);
  rel(besselI1(10), 2670.988303701255, 1e-15);
  rel(besselK0(10), 1.778006231616765e-5, 1e-14);
  rel(besselK1(10), 1.8648773453825585e-5, 1e-14);
  assert.equal(besselI0(0), 1);
  assert.equal(besselI1(0), 0);
  assert.equal(besselI0(-2), besselI0(2));
  assert.equal(besselI1(-2), -besselI1(2));
  assert.throws(() => besselK0(0), RangeError);
  assert.throws(() => besselK1(-2), RangeError);
});

test('modified Bessel Wronskian I0 K1 + I1 K0 = 1/x across all internal regimes (1e-3 .. 600)', () => {
  for (const x of [1e-3, 0.02, 0.5, 0.999999, 1, 1.000001, 1.7, 3, 9, 20, 29.999, 30, 30.001, 75, 300, 600]) {
    rel(besselI0(x) * besselK1(x) + besselI1(x) * besselK0(x), 1 / x, 2e-14, `x=${x}`);
  }
});

test('modified Bessel asymptotics (A&S 9.7.1–9.7.2, three terms, next term ~1e-8 at x = 200); small-x limits; derivative identities', () => {
  const x = 200;
  rel(besselK0(x), Math.sqrt(Math.PI / (2 * x)) * Math.exp(-x) * (1 - 1 / (8 * x) + 9 / (128 * x * x)), 5e-8);
  rel(besselI0(x), (Math.exp(x) / Math.sqrt(2 * Math.PI * x)) * (1 + 1 / (8 * x) + 9 / (128 * x * x)), 5e-8);
  rel(besselK1(x), Math.sqrt(Math.PI / (2 * x)) * Math.exp(-x) * (1 + 3 / (8 * x) - 15 / (128 * x * x)), 5e-8);
  const s = 1e-5;
  rel(besselK0(s), -(Math.log(s / 2) + 0.5772156649015329), 1e-9);
  rel(besselK1(s), 1 / s, 1e-9);
  // derivative identities I0' = I1 and K0' = −K1 by central differences
  const h = 1e-5;
  close((besselI0(2.3 + h) - besselI0(2.3 - h)) / (2 * h), besselI1(2.3), 1e-9);
  close((besselK0(2.3 + h) - besselK0(2.3 - h)) / (2 * h), -besselK1(2.3), 1e-9);
});
