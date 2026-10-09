import test from 'node:test';
import assert from 'node:assert/strict';
import { eigSym, eigGenSym, subspaceIteration, eigGeneral } from '../../src/engine/math/eigen.js';
import { SkylineSym } from '../../src/engine/math/sparse.js';
import { det, zeros } from '../../src/engine/math/linalg.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

/** Linear congruential generator (Numerical Recipes constants), uniform in [0, 1). */
function lcg(seed) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function tridiag(n, lo, di) {
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) {
    A[i][i] = di;
    if (i > 0) {
      A[i][i - 1] = lo;
      A[i - 1][i] = lo;
    }
  }
  return A;
}

function mulMV(A, x) {
  return A.map((row) => row.reduce((s, v, j) => s + v * x[j], 0));
}

function dotv(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

test('eigSym: 2×2 hand example [[2,1],[1,2]] -> λ = 1, 3 with vectors (1,−1)/√2 and (1,1)/√2', () => {
  const { values, vectors } = eigSym([[2, 1], [1, 2]]);
  close(values[0], 1, 1e-15);
  close(values[1], 3, 1e-15);
  close(Math.abs(vectors[0][0]), Math.SQRT1_2, 1e-15);
  close(vectors[0][0] + vectors[0][1], 0, 1e-15);
  close(vectors[1][0], Math.SQRT1_2, 1e-15);
  close(vectors[1][1], Math.SQRT1_2, 1e-15);
});

test('eigSym: tridiagonal [−1, 2, −1] has λ_k = 2 − 2 cos(kπ/(n+1)) and sine eigenvectors', () => {
  for (const n of [1, 2, 7, 30, 100]) {
    const { values, vectors } = eigSym(tridiag(n, -1, 2));
    assert.equal(values.length, n);
    for (let k = 1; k <= n; k++) close(values[k - 1], 2 - 2 * Math.cos((k * Math.PI) / (n + 1)), 1e-13, `n=${n}, k=${k}`);
    if (n >= 7) {
      // first mode shape sqrt(2/(n+1)) sin(iπ/(n+1)), positive by the sign convention
      const v = vectors[0];
      for (let i = 0; i < n; i += 3) close(v[i], Math.sqrt(2 / (n + 1)) * Math.sin(((i + 1) * Math.PI) / (n + 1)), 1e-11);
    }
  }
});

test('eigSym: random symmetric matrix (seeded LCG) — residual, orthonormality, trace and Frobenius invariants', () => {
  const rnd = lcg(123);
  const n = 60;
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    A[i][j] = rnd() - 0.5;
    A[j][i] = A[i][j];
  }
  const copy = A.map((r) => r.slice());
  const { values, vectors } = eigSym(A);
  assert.deepEqual(A, copy, 'input untouched');
  let tr = 0;
  let fro = 0;
  for (let i = 0; i < n; i++) {
    tr += A[i][i];
    for (let j = 0; j < n; j++) fro += A[i][j] * A[i][j];
  }
  close(values.reduce((s, v) => s + v, 0), tr, 1e-12);
  close(values.reduce((s, v) => s + v * v, 0), fro, 1e-11);
  for (let k = 0; k < n; k++) {
    if (k > 0) assert.ok(values[k] >= values[k - 1]);
    const Av = mulMV(A, vectors[k]);
    for (let i = 0; i < n; i++) close(Av[i], values[k] * vectors[k][i], 1e-12, `residual of pair ${k}`);
    for (let l = 0; l <= k; l++) close(dotv(vectors[k], vectors[l]), k === l ? 1 : 0, 1e-12);
  }
});

test('eigSym: repeated eigenvalues, diagonal and rank-one matrices; Hilbert(5) extreme eigenvalues', () => {
  const d = eigSym([[3, 0, 0], [0, 1, 0], [0, 0, 3]]);
  assert.deepEqual(d.values, [1, 3, 3]);
  close(Math.abs(d.vectors[0][1]), 1, 1e-15);
  close(dotv(d.vectors[1], d.vectors[2]), 0, 1e-15);
  // all-ones 4×4: eigenvalues 0, 0, 0, 4
  const ones = eigSym([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]);
  for (let k = 0; k < 3; k++) close(ones.values[k], 0, 1e-15);
  close(ones.values[3], 4, 1e-14);
  for (const v of ones.vectors[3]) close(v, 0.5, 1e-14);
  // Hilbert(5): λ_min = 3.2879e-6, λ_max = 1.5670507 (Gregory & Karney, test matrix 3.8)
  const H = zeros(5, 5);
  for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) H[i][j] = 1 / (i + j + 1);
  const h = eigSym(H);
  close(h.values[4], 1.5670506910982311, 1e-13);
  close(h.values[0] / 3.2879287721718e-6, 1, 1e-8);
  assert.throws(() => eigSym([[1, 2], [3, 4]]), RangeError);
  assert.throws(() => eigSym([[1, 2, 3], [2, 1, 3]]), RangeError);
});

test('eigGenSym: 2-DOF shear building K = k[[2,−1],[−1,1]], M = m I — ω² = (k/m)(3 ∓ √5)/2', () => {
  const k = 1.2e6;
  const m = 800;
  const K = [[2 * k, -k], [-k, k]];
  const M = [[m, 0], [0, m]];
  const { values, vectors } = eigGenSym(K, M);
  close(values[0] / (((k / m) * (3 - Math.sqrt(5))) / 2), 1, 1e-13);
  close(values[1] / (((k / m) * (3 + Math.sqrt(5))) / 2), 1, 1e-13);
  // first mode shape ratio φ2/φ1 = golden ratio
  close(vectors[0][1] / vectors[0][0], (1 + Math.sqrt(5)) / 2, 1e-12);
  // M-orthonormal: φᵀ M φ = 1 -> |φ|² = 1/m
  close(dotv(vectors[0], vectors[0]) * m, 1, 1e-12);
  close(dotv(vectors[0], vectors[1]), 0, 1e-15);
  // lumped mass given as a vector
  const v = eigGenSym(K, [m, m]);
  close(v.values[0], values[0], 1e-9 * values[0]);
});

test('eigGenSym: fixed–fixed rod with consistent mass — K = tridiag(−1,2,−1), M = tridiag(1,4,1)/6, λ_k = 6(1−c)/(2+c), c = cos(kπ/(n+1))', () => {
  const n = 40;
  const K = tridiag(n, -1, 2);
  const M = tridiag(n, 1 / 6, 4 / 6);
  const { values, vectors } = eigGenSym(K, M, { nModes: 6 });
  assert.equal(values.length, 6);
  assert.equal(vectors.length, 6);
  for (let k = 1; k <= 6; k++) {
    const c = Math.cos((k * Math.PI) / (n + 1));
    close(values[k - 1], (6 * (1 - c)) / (2 + c), 1e-13, `mode ${k}`);
  }
  // residual K φ = λ M φ and M-orthonormality
  for (let a = 0; a < 6; a++) {
    const Kv = mulMV(K, vectors[a]);
    const Mv = mulMV(M, vectors[a]);
    for (let i = 0; i < n; i++) close(Kv[i], values[a] * Mv[i], 1e-12);
    for (let b = 0; b <= a; b++) close(dotv(vectors[b], Mv), a === b ? 1 : 0, 1e-12);
  }
  // consistent mass overestimates the continuum eigenvalue (kπ h)² — upper-bound property of Rayleigh–Ritz
  const h = 1 / (n + 1);
  assert.ok(values[0] > (Math.PI * h) ** 2);
  close(values[0] / (Math.PI * h) ** 2, 1, 1e-3);
  assert.throws(() => eigGenSym(K, tridiag(n, 1, 1)), RangeError, 'indefinite M');
  assert.throws(() => eigGenSym(K, M, { nModes: 0 }), RangeError);
});

function springChain(n, k) {
  // fixed–fixed chain of n masses and n+1 springs
  const K = new SkylineSym(n, new Int32Array(n).fill(1));
  for (let i = 0; i < n; i++) {
    K.add(i, i, 2 * k);
    if (i > 0) K.add(i - 1, i, -k);
  }
  return K.factor();
}

test('subspaceIteration: 200-dof spring–mass chain vs closed form ω_j² = (2k/m)(1 − cos(jπ/(n+1))), 8 modes', () => {
  const n = 200;
  const k = 4.0e5;
  const m = 25;
  const K = springChain(n, k);
  const out = new Float64Array(n);
  const r = subspaceIteration({
    n, nModes: 8,
    solveK: (b) => K.solve(b),
    mulM: (x) => {
      for (let i = 0; i < n; i++) out[i] = m * x[i]; // deliberately reuses its buffer
      return out;
    },
    tol: 1e-10,
  });
  assert.ok(r.converged, `not converged after ${r.iterations} sweeps`);
  assert.ok(r.iterations <= 30, `sweeps ${r.iterations}`);
  assert.equal(r.values.length, 8);
  for (let j = 1; j <= 8; j++) {
    const exact = ((2 * k) / m) * (1 - Math.cos((j * Math.PI) / (n + 1)));
    close(r.values[j - 1] / exact, 1, 1e-8, `mode ${j}`);
  }
  // mode shapes: M-orthonormal sines, sqrt(2/((n+1) m)) sin(i j π/(n+1)) up to sign
  for (let a = 0; a < 8; a++) {
    assert.ok(r.vectors[a] instanceof Float64Array);
    for (let b = 0; b <= a; b++) close(m * dotv(r.vectors[a], r.vectors[b]), a === b ? 1 : 0, 1e-9);
  }
  for (const j of [1, 2, 5]) {
    const v = r.vectors[j - 1];
    const amp = Math.sqrt(2 / ((n + 1) * m));
    const sgn = Math.sign(v[0]);
    for (let i = 0; i < n; i += 17) close(sgn * v[i], amp * Math.sin(((i + 1) * j * Math.PI) / (n + 1)), 1e-6 * amp * 100, `shape ${j}`);
  }
  // fundamental frequency in Hz for reference: f1 = ω1/2π
  close(Math.sqrt(r.values[0]) / (2 * Math.PI), (Math.sqrt((2 * k) / m) * Math.sqrt(1 - Math.cos(Math.PI / 201))) / (2 * Math.PI), 1e-8);
});

test('subspaceIteration: non-uniform consistent mass — agrees with the dense eigGenSym solution (60 dof, 5 modes)', () => {
  const n = 60;
  const rnd = lcg(77);
  // tapered rod: element stiffness and mass vary along the length
  const ke = Array.from({ length: n + 1 }, (_, e) => 1 + 2 * (e / n) + 0.2 * rnd());
  const me = Array.from({ length: n + 1 }, (_, e) => 2 - e / n + 0.2 * rnd());
  const Kd = zeros(n, n);
  const Md = zeros(n, n);
  for (let e = 0; e <= n; e++) {
    const a = e - 1;
    const b = e;
    if (a >= 0) {
      Kd[a][a] += ke[e];
      Md[a][a] += me[e] / 3;
    }
    if (b < n) {
      Kd[b][b] += ke[e];
      Md[b][b] += me[e] / 3;
    }
    if (a >= 0 && b < n) {
      Kd[a][b] -= ke[e];
      Kd[b][a] -= ke[e];
      Md[a][b] += me[e] / 6;
      Md[b][a] += me[e] / 6;
    }
  }
  const K = new SkylineSym(n, new Int32Array(n).fill(1));
  for (let i = 0; i < n; i++) {
    K.add(i, i, Kd[i][i]);
    if (i > 0) K.add(i - 1, i, Kd[i - 1][i]);
  }
  K.factor();
  const ref = eigGenSym(Kd, Md, { nModes: 5 });
  const r = subspaceIteration({ n, nModes: 5, solveK: (b) => K.solve(b), mulM: (x) => mulMV(Md, x), tol: 1e-11, seed: 5 });
  assert.ok(r.converged);
  for (let j = 0; j < 5; j++) {
    close(r.values[j] / ref.values[j], 1, 1e-9, `mode ${j + 1}`);
    const s = Math.sign(dotv(r.vectors[j], mulMV(Md, ref.vectors[j])));
    for (let i = 0; i < n; i += 7) close(s * r.vectors[j][i], ref.vectors[j][i], 1e-6);
  }
});

test('subspaceIteration: deterministic for a seed; all modes of a small system; argument checks', () => {
  const n = 12;
  const K = springChain(n, 1);
  const args = { n, solveK: (b) => K.solve(b), mulM: (x) => x };
  const a = subspaceIteration({ ...args, nModes: 3, seed: 9 });
  const b = subspaceIteration({ ...args, nModes: 3, seed: 9 });
  assert.deepEqual(a.values, b.values);
  assert.deepEqual(a.vectors, b.vectors);
  // subspace spans the whole space: exact after the first Rayleigh–Ritz projection
  const full = subspaceIteration({ ...args, nModes: n });
  assert.ok(full.converged);
  assert.ok(full.iterations <= 3);
  // the first projection uses the raw (ill-conditioned) basis K⁻¹X, hence the slightly relaxed tolerance
  for (let j = 1; j <= n; j++) close(full.values[j - 1], 2 - 2 * Math.cos((j * Math.PI) / (n + 1)), 1e-9);
  assert.throws(() => subspaceIteration({ ...args, nModes: 13 }), RangeError);
  assert.throws(() => subspaceIteration({ ...args, nModes: 2, solveK: null }), RangeError);
  assert.throws(() => subspaceIteration({ ...args, nModes: 2, mulM: () => [1, 2] }), RangeError);
  const capped = subspaceIteration({ n: 200, nModes: 6, solveK: (v) => springChain(200, 1).solve(v), mulM: (x) => x, maxIter: 1 });
  assert.equal(capped.converged, false);
  assert.equal(capped.iterations, 1);
});

test('subspaceIteration: massless degrees of freedom (singular M) are tolerated as long as rank(M) >= nModes', () => {
  // 30-dof chain where every third dof carries no mass
  const n = 30;
  const K = springChain(n, 1);
  const mass = Float64Array.from({ length: n }, (_, i) => (i % 3 === 1 ? 0 : 1));
  const r = subspaceIteration({ n, nModes: 3, solveK: (b) => K.solve(b), mulM: (x) => x.map((v, i) => v * mass[i]), tol: 1e-10 });
  assert.ok(r.converged);
  // verify K φ = λ M φ directly
  for (let j = 0; j < 3; j++) {
    const v = r.vectors[j];
    for (let i = 0; i < n; i++) {
      const Kv = 2 * v[i] - (i > 0 ? v[i - 1] : 0) - (i < n - 1 ? v[i + 1] : 0);
      close(Kv, r.values[j] * mass[i] * v[i], 1e-6);
    }
  }
});

test('eigGeneral: rotation generator [[0,1],[−1,0]] -> ±i; upper-triangular matrix returns its diagonal', () => {
  const r = eigGeneral([[0, 1], [-1, 0]]);
  assert.deepEqual(r.re, [0, 0]);
  assert.deepEqual(r.im, [-1, 1]);
  const t = eigGeneral([[3, 5, -2], [0, -1, 7], [0, 0, 2]]);
  assert.deepEqual(t.re, [-1, 2, 3]);
  assert.deepEqual(t.im, [0, 0, 0]);
  assert.deepEqual(eigGeneral([[7]]), { re: [7], im: [0] });
});

test('eigGeneral: companion matrix of (x−1)(x−2)(x−3)(x−4); matrix [[1,2,3],[4,5,6],[7,8,10]]', () => {
  const c = eigGeneral([[0, 0, 0, -24], [1, 0, 0, 50], [0, 1, 0, -35], [0, 0, 1, 10]]);
  for (let k = 0; k < 4; k++) {
    close(c.re[k], k + 1, 1e-11);
    assert.equal(c.im[k], 0);
  }
  // characteristic polynomial λ³ − 16λ² − 12λ + 3 = 0
  const m = eigGeneral([[1, 2, 3], [4, 5, 6], [7, 8, 10]]);
  for (let k = 0; k < 3; k++) {
    const l = m.re[k];
    close(l * l * l - 16 * l * l - 12 * l + 3, 0, 1e-10);
  }
  close(m.re[2], 16.707493316124740, 1e-12);
});

test('eigGeneral: state matrix of a damped oscillator — eigenvalues −ζω ± iω√(1−ζ²); overdamped pair is real', () => {
  const w = 7;
  const zeta = 0.03;
  const r = eigGeneral([[0, 1], [-w * w, -2 * zeta * w]]);
  close(r.re[0], -zeta * w, 1e-13);
  close(r.re[1], -zeta * w, 1e-13);
  close(r.im[1], w * Math.sqrt(1 - zeta * zeta), 1e-13);
  close(r.im[0], -r.im[1], 0);
  const o = eigGeneral([[0, 1], [-w * w, -2 * 2 * w]]);
  assert.deepEqual(o.im, [0, 0]);
  close(o.re[0], -w * (2 + Math.sqrt(3)), 1e-12);
  close(o.re[1], -w * (2 - Math.sqrt(3)), 1e-12);
  // 2-DOF non-proportionally damped system: 4 eigenvalues in two conjugate pairs with negative real parts
  const K = [[2, -1], [-1, 1]];
  const C = [[0.1, 0], [0, 0]];
  const A = [
    [0, 0, 1, 0], [0, 0, 0, 1],
    [-K[0][0], -K[0][1], -C[0][0], -C[0][1]], [-K[1][0], -K[1][1], -C[1][0], -C[1][1]],
  ];
  const s = eigGeneral(A);
  let sumRe = 0;
  for (let k = 0; k < 4; k++) {
    assert.ok(s.re[k] < 0, 'stable');
    sumRe += s.re[k];
  }
  close(sumRe, -0.1, 1e-13, 'trace');
  close(s.im.reduce((p, v) => p + v, 0), 0, 1e-13);
});

test('eigGeneral: random unsymmetric matrices (seeded LCG) — Σλ = trace, Πλ = det, conjugate pairs', () => {
  const rnd = lcg(2025);
  for (const n of [5, 12, 30]) {
    const A = zeros(n, n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) A[i][j] = rnd() - 0.5;
    const copy = A.map((r) => r.slice());
    const { re, im } = eigGeneral(A);
    assert.deepEqual(A, copy);
    let tr = 0;
    for (let i = 0; i < n; i++) tr += A[i][i];
    close(re.reduce((s, v) => s + v, 0), tr, 1e-11, `n=${n} trace`);
    close(im.reduce((s, v) => s + v, 0), 0, 1e-11);
    // product of eigenvalues (complex arithmetic) equals the determinant
    let pr = 1;
    let pi = 0;
    for (let k = 0; k < n; k++) {
      const a = pr * re[k] - pi * im[k];
      pi = pr * im[k] + pi * re[k];
      pr = a;
    }
    const d = det(A);
    close(pr / d, 1, 1e-8, `n=${n} determinant`);
    close(pi / d, 0, 1e-8);
    // Σλ² = trace(A²)
    let tr2 = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) tr2 += A[i][j] * A[j][i];
    close(re.reduce((s, v, k) => s + v * v - im[k] * im[k], 0), tr2, 1e-10);
  }
});

test('eigGeneral: symmetric input reproduces eigSym; badly scaled matrix is handled by balancing', () => {
  const A = tridiag(15, -1, 2);
  const g = eigGeneral(A);
  const s = eigSym(A);
  for (let k = 0; k < 15; k++) {
    close(g.re[k], s.values[k], 1e-12);
    assert.equal(g.im[k], 0);
  }
  // similarity transform D⁻¹ B D with D = diag(1, 1e6, 1e-6) keeps the eigenvalues 1, 2, 3 of the companion-like B
  const B = [[6, -11, 6], [1, 0, 0], [0, 1, 0]];
  const dsc = [1, 1e6, 1e-6];
  const S = B.map((row, i) => row.map((v, j) => (v * dsc[j]) / dsc[i]));
  const r = eigGeneral(S);
  for (let k = 0; k < 3; k++) close(r.re[k], k + 1, 1e-9);
  assert.throws(() => eigGeneral([[1, 2], [3]]), RangeError);
  assert.throws(() => eigGeneral([[1, NaN], [3, 4]]), RangeError);
});
