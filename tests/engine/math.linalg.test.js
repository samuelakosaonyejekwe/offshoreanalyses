import test from 'node:test';
import assert from 'node:assert/strict';
import {
  zeros, identity, matMul, matVec, transpose, dot, norm, axpy, solveDense, luFactor, luSolve, cholesky, cholSolve,
  inverse, det, leastSquares, svd, pinv, solveTridiagonal, solveBandedSym, factorBandedSym, solveBandedSymFactored,
} from '../../src/engine/math/linalg.js';

const EPS = 2.220446049250313e-16;

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

function hilbert(n) {
  const H = zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) H[i][j] = 1 / (i + j + 1);
  return H;
}

function randomSpd(n, rnd) {
  // A = Bᵀ B + n I
  const B = zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) B[i][j] = rnd() - 0.5;
  const A = matMul(transpose(B), B);
  for (let i = 0; i < n; i++) A[i][i] += 0.1 * n;
  return A;
}

function maxAbsDiff(A, B) {
  let e = 0;
  for (let i = 0; i < A.length; i++) for (let j = 0; j < A[i].length; j++) e = Math.max(e, Math.abs(A[i][j] - B[i][j]));
  return e;
}

test('zeros, identity, transpose, matMul, matVec: hand-checked 2×3 example', () => {
  assert.deepEqual(zeros(3), [0, 0, 0]);
  assert.deepEqual(zeros(2, 3), [[0, 0, 0], [0, 0, 0]]);
  assert.deepEqual(identity(2), [[1, 0], [0, 1]]);
  const A = [[1, 2, 3], [4, 5, 6]];
  assert.deepEqual(transpose(A), [[1, 4], [2, 5], [3, 6]]);
  // A Aᵀ = [[14, 32], [32, 77]]
  assert.deepEqual(matMul(A, transpose(A)), [[14, 32], [32, 77]]);
  assert.deepEqual(matVec(A, [1, 0, -1]), [-2, -2]);
  assert.deepEqual(matVec(A, new Float64Array([1, 1, 1])), [6, 15]);
  assert.throws(() => matMul(A, A), RangeError);
  assert.throws(() => matVec(A, [1, 2]), RangeError);
});

test('dot, norm, axpy; norm is overflow- and underflow-safe', () => {
  assert.equal(dot([1, 2, 3], new Float64Array([4, 5, 6])), 32);
  assert.equal(norm([3, 4]), 5);
  close(norm([3e200, 4e200]) / 5e200, 1, 1e-15);
  close(norm([3e-200, 4e-200]) / 5e-200, 1, 1e-15);
  assert.equal(norm([0, 0]), 0);
  assert.deepEqual(axpy(2, [1, 2], [10, 20]), [12, 24]);
  const x = [1, 2];
  axpy(3, x, x);
  assert.deepEqual(x, [1, 2], 'inputs are not mutated');
});

test('solveDense: Hilbert matrices, error bounded by cond(H)·eps (cond from Todd 1954: 4.8e5, 1.5e10, 1.6e13)', () => {
  const conds = { 5: 4.766e5, 8: 1.526e10, 10: 1.603e13 };
  for (const n of [5, 8, 10]) {
    const H = hilbert(n);
    const ones = new Array(n).fill(1);
    const b = matVec(H, ones);
    const x = solveDense(H, b);
    let err = 0;
    for (let i = 0; i < n; i++) err = Math.max(err, Math.abs(x[i] - 1));
    assert.ok(err <= conds[n] * EPS * 4, `n=${n}: error ${err} exceeds cond·eps bound ${conds[n] * EPS * 4}`);
    // backward stability: the residual is at round-off level regardless of conditioning
    const r = matVec(H, x);
    let res = 0;
    for (let i = 0; i < n; i++) res = Math.max(res, Math.abs(r[i] - b[i]));
    assert.ok(res < 1e-14, `n=${n}: residual ${res}`);
  }
});

test('solveDense: needs pivoting; multi-RHS; inputs untouched; singular throws', () => {
  // zero leading pivot: [[0,1],[1,0]] x = [2,3] -> x = [3,2]
  assert.deepEqual(solveDense([[0, 1], [1, 0]], [2, 3]), [3, 2]);
  const A = [[4, -2, 1], [3, 6, -4], [2, 1, 8]];
  const copy = A.map((r) => r.slice());
  // right-hand sides built from known solutions [1,2,3] and [-1,0,1]
  const B = [[3, -3], [3, -7], [28, 6]];
  const X = solveDense(A, B);
  assert.deepEqual(A, copy);
  close(X[0][0], 1, 1e-14);
  close(X[1][0], 2, 1e-14);
  close(X[2][0], 3, 1e-14);
  close(X[0][1], -1, 1e-14);
  close(X[1][1], 0, 1e-14);
  close(X[2][1], 1, 1e-14);
  assert.throws(() => solveDense([[1, 2], [2, 4]], [1, 1]), RangeError);
  assert.throws(() => solveDense([[1, 2, 3], [4, 5, 6]], [1, 1]), RangeError);
});

test('luFactor/luSolve reuse, inverse and det (det of Hilbert(3) = 1/2160)', () => {
  const rnd = lcg(11);
  const n = 25;
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) A[i][j] = rnd() - 0.5;
  const f = luFactor(A);
  for (let k = 0; k < 3; k++) {
    const xt = Array.from({ length: n }, () => rnd());
    const x = luSolve(f, matVec(A, xt));
    for (let i = 0; i < n; i++) close(x[i], xt[i], 1e-11);
  }
  const Ainv = inverse(A);
  assert.ok(maxAbsDiff(matMul(A, Ainv), identity(n)) < 1e-11);
  close(det(hilbert(3)), 1 / 2160, 1e-17);
  close(det([[1, 2], [3, 4]]), -2, 1e-15);
  // odd permutation of the identity
  assert.equal(det([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), -1);
  assert.equal(det([[1, 2], [2, 4]]), 0);
  // det(A) det(A⁻¹) = 1
  close(det(A) * det(Ainv), 1, 1e-10);
  assert.throws(() => luFactor([[1, 2], [2, 4]]), RangeError);
});

test('cholesky: random SPD systems (seeded LCG) — L Lᵀ = A and residual at round-off', () => {
  const rnd = lcg(2024);
  for (const n of [5, 30, 80]) {
    const A = randomSpd(n, rnd);
    const L = cholesky(A);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) assert.equal(L[i][j], 0);
    assert.ok(maxAbsDiff(matMul(L, transpose(L)), A) < 1e-12 * n);
    const b = Array.from({ length: n }, () => rnd() - 0.5);
    const x = cholSolve(L, b);
    const r = matVec(A, x);
    let res = 0;
    for (let i = 0; i < n; i++) res = Math.max(res, Math.abs(r[i] - b[i]));
    assert.ok(res < 1e-12, `n=${n}: residual ${res}`);
    // agrees with LU
    const xl = solveDense(A, b);
    for (let i = 0; i < n; i++) close(x[i], xl[i], 1e-11);
  }
});

test('cholesky: hand example [[4,2],[2,3]] -> L = [[2,0],[1,√2]]; indefinite matrix throws', () => {
  const L = cholesky([[4, 2], [2, 3]]);
  close(L[0][0], 2, 0);
  close(L[1][0], 1, 0);
  close(L[1][1], Math.SQRT2, 1e-15);
  assert.throws(() => cholesky([[1, 2], [2, 1]]), RangeError);
  assert.throws(() => cholesky([[0, 0], [0, 1]]), RangeError);
});

test('leastSquares: straight-line fit through (0,1), (1,2), (2,4) — slope 3/2, intercept 5/6, residual 1/√6', () => {
  const { x, residual, rank } = leastSquares([[1, 0], [1, 1], [1, 2]], [1, 2, 4]);
  close(x[0], 5 / 6, 1e-14);
  close(x[1], 1.5, 1e-14);
  close(residual, 1 / Math.sqrt(6), 1e-14);
  assert.equal(rank, 2);
});

test('leastSquares: consistent overdetermined system is solved exactly; normal equations satisfied for noisy data', () => {
  const rnd = lcg(5);
  const m = 40;
  const n = 6;
  const A = zeros(m, n);
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) A[i][j] = rnd() - 0.5;
  const xt = [1, -2, 3, -4, 5, -6];
  const exact = leastSquares(A, matVec(A, xt));
  for (let j = 0; j < n; j++) close(exact.x[j], xt[j], 1e-12);
  assert.ok(exact.residual < 1e-12);
  const b = matVec(A, xt).map((v) => v + 0.1 * (rnd() - 0.5));
  const ls = leastSquares(A, b);
  // Aᵀ (A x − b) = 0
  const r = matVec(A, ls.x).map((v, i) => v - b[i]);
  const g = matVec(transpose(A), r);
  for (let j = 0; j < n; j++) close(g[j], 0, 1e-13);
  close(ls.residual, norm(r), 1e-13);
  assert.throws(() => leastSquares([[1, 2, 3]], [1]), RangeError);
});

test('leastSquares: rank-deficient matrix returns the minimum-norm solution', () => {
  // columns identical: x1 + x2 = t; minimum norm splits evenly
  const { x, rank } = leastSquares([[1, 1], [1, 1], [1, 1]], [2, 2, 2]);
  close(x[0], 1, 1e-12);
  close(x[1], 1, 1e-12);
  assert.equal(rank, 1);
});

test('svd: singular values of [[1,2],[3,4],[5,6]] = sqrt(eig(AᵀA)) = 9.525518091565, 0.514300580659', () => {
  const { U, S, V } = svd([[1, 2], [3, 4], [5, 6]]);
  // eig of [[35,44],[44,56]]: (91 ± sqrt(21² + 4·44²))/2
  const disc = Math.sqrt(21 * 21 + 4 * 44 * 44);
  close(S[0], Math.sqrt((91 + disc) / 2), 1e-13);
  close(S[1], Math.sqrt((91 - disc) / 2), 1e-13);
  assert.equal(U.length, 3);
  assert.equal(U[0].length, 2);
  assert.equal(V.length, 2);
});

test('svd: reconstruction and orthogonality for tall, wide and rank-deficient matrices', () => {
  const rnd = lcg(99);
  for (const [m, n] of [[7, 4], [3, 6], [5, 5]]) {
    const A = zeros(m, n);
    for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) A[i][j] = rnd() - 0.5;
    const { U, S, V } = svd(A);
    const k = Math.min(m, n);
    assert.equal(S.length, k);
    for (let i = 1; i < k; i++) assert.ok(S[i] <= S[i - 1], 'descending');
    const US = U.map((row) => row.map((v, j) => v * S[j]));
    assert.ok(maxAbsDiff(matMul(US, transpose(V)), A) < 1e-13);
    assert.ok(maxAbsDiff(matMul(transpose(U), U), identity(k)) < 1e-13);
    assert.ok(maxAbsDiff(matMul(transpose(V), V), identity(k)) < 1e-13);
  }
  // rank 1: outer product [1,2,3]ᵀ[1,1] has the single singular value sqrt(14)·sqrt(2)
  const { S } = svd([[1, 1], [2, 2], [3, 3]]);
  close(S[0], Math.sqrt(28), 1e-13);
  close(S[1], 0, 1e-14);
});

test('svd: 2-norm condition number of Hilbert(6) = 1.4951e7 (Todd 1954)', () => {
  const { S } = svd(hilbert(6));
  close(S[0] / S[5] / 1.4951058642e7, 1, 1e-6);
});

test('pinv: Moore–Penrose conditions; equals inverse for a regular matrix; rank-deficient case', () => {
  const A = [[1, 2], [3, 4], [5, 6], [7, 9]];
  const P = pinv(A);
  assert.ok(maxAbsDiff(matMul(matMul(A, P), A), A) < 1e-12);
  assert.ok(maxAbsDiff(matMul(matMul(P, A), P), P) < 1e-12);
  const AP = matMul(A, P);
  assert.ok(maxAbsDiff(AP, transpose(AP)) < 1e-12);
  const B = [[2, 1], [1, 3]];
  assert.ok(maxAbsDiff(pinv(B), inverse(B)) < 1e-14);
  // pinv of [[1,1],[1,1]] is [[1,1],[1,1]]/4
  const Q = pinv([[1, 1], [1, 1]]);
  for (const row of Q) for (const v of row) close(v, 0.25, 1e-14);
});

test('solveTridiagonal: −u″ = 1 discretised, tridiag(−1,2,−1) x = 1 has x_i = i(n+1−i)/2', () => {
  const n = 50;
  const lo = new Array(n).fill(-1);
  const di = new Array(n).fill(2);
  const up = new Array(n).fill(-1);
  const x = solveTridiagonal(lo, di, up, new Array(n).fill(1));
  assert.ok(x instanceof Float64Array);
  for (let i = 1; i <= n; i++) close(x[i - 1], (i * (n + 1 - i)) / 2, 1e-10);
  // off-diagonals of length n−1 are accepted too
  const x2 = solveTridiagonal(new Float64Array(n - 1).fill(-1), di, new Float64Array(n - 1).fill(-1), new Float64Array(n).fill(1));
  for (let i = 0; i < n; i++) close(x2[i], x[i], 1e-12);
});

test('solveTridiagonal: unsymmetric system against dense LU; zero pivot throws', () => {
  const rnd = lcg(3);
  const n = 12;
  const lo = Array.from({ length: n }, () => rnd() - 0.5);
  const up = Array.from({ length: n }, () => rnd() - 0.5);
  const di = Array.from({ length: n }, () => 2 + rnd());
  const b = Array.from({ length: n }, () => rnd());
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) {
    A[i][i] = di[i];
    if (i > 0) A[i][i - 1] = lo[i];
    if (i < n - 1) A[i][i + 1] = up[i];
  }
  const x = solveTridiagonal(lo, di, up, b);
  const xd = solveDense(A, b);
  for (let i = 0; i < n; i++) close(x[i], xd[i], 1e-13);
  assert.throws(() => solveTridiagonal([0, 1], [0, 1], [1, 0], [1, 1]), RangeError);
});

test('banded SPD Cholesky: agrees with dense solve on a random banded SPD matrix', () => {
  const rnd = lcg(77);
  const n = 40;
  const bw = 5;
  const bands = [];
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) {
    const row = [2 * bw + 1 + rnd()];
    for (let k = 1; k <= bw && i + k < n; k++) row.push(rnd() - 0.5);
    bands.push(row);
    for (let k = 0; k < row.length; k++) {
      A[i][i + k] = row[k];
      A[i + k][i] = row[k];
    }
  }
  const b = Array.from({ length: n }, () => rnd() - 0.5);
  const x = solveBandedSym(bands, b);
  assert.ok(x instanceof Float64Array);
  const xd = solveDense(A, b);
  for (let i = 0; i < n; i++) close(x[i], xd[i], 1e-13);
  // factor once, solve twice
  const f = factorBandedSym(bands);
  const b2 = Array.from({ length: n }, () => rnd());
  const x2 = solveBandedSymFactored(f, b2);
  const xd2 = solveDense(A, b2);
  for (let i = 0; i < n; i++) close(x2[i], xd2[i], 1e-13);
});

test('banded SPD Cholesky: beam-like pentadiagonal system of 20 000 unknowns, residual at round-off', () => {
  // fourth-difference operator [1,−4,6,−4,1] plus a spring foundation: SPD
  const n = 20000;
  const bands = new Array(n);
  for (let i = 0; i < n; i++) bands[i] = [6.5, -4, 1];
  const b = new Float64Array(n).fill(1);
  const x = solveBandedSym(bands, b);
  let res = 0;
  for (let i = 0; i < n; i++) {
    let s = 6.5 * x[i];
    if (i >= 1) s -= 4 * x[i - 1];
    if (i >= 2) s += x[i - 2];
    if (i + 1 < n) s -= 4 * x[i + 1];
    if (i + 2 < n) s += x[i + 2];
    res = Math.max(res, Math.abs(s - 1));
  }
  assert.ok(res < 1e-11, `residual ${res}`);
  // far from the ends the solution tends to the uniform foundation value 1/0.5
  close(x[n >> 1], 2, 1e-9);
});

test('banded SPD Cholesky: non-positive-definite matrix throws', () => {
  assert.throws(() => factorBandedSym([[1, 2], [1]]), RangeError);
  assert.throws(() => factorBandedSym([]), RangeError);
});
