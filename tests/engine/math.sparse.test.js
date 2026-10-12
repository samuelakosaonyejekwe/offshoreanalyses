import test from 'node:test';
import assert from 'node:assert/strict';
import { CooBuilder, csrMatVec, pcg, SkylineSym, reverseCuthillMcKee, skylineHeights } from '../../src/engine/math/sparse.js';
import { solveDense, zeros } from '../../src/engine/math/linalg.js';

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

/** Five-point Laplacian on an nx × ny grid with Dirichlet boundaries (SPD), plus its coupling edges. */
function laplace2d(nx, ny) {
  const n = nx * ny;
  const coo = new CooBuilder(n);
  const edges = [];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const k = j * nx + i;
      coo.add(k, k, 4);
      if (i > 0) {
        coo.add(k, k - 1, -1);
        coo.add(k - 1, k, -1);
        edges.push([k - 1, k]);
      }
      if (j > 0) {
        coo.add(k, k - nx, -1);
        coo.add(k - nx, k, -1);
        edges.push([k - nx, k]);
      }
    }
  }
  return { A: coo.toCsr(), edges, n };
}

function residualInf(A, x, b) {
  const r = csrMatVec(A, x);
  let e = 0;
  for (let i = 0; i < b.length; i++) e = Math.max(e, Math.abs(r[i] - b[i]));
  return e;
}

test('CooBuilder: duplicates are summed, columns sorted, typed CSR arrays, growth beyond initial capacity', () => {
  const coo = new CooBuilder(3);
  coo.add(2, 2, 5);
  coo.add(0, 1, 1);
  coo.add(0, 0, 2);
  coo.add(0, 1, 3);
  coo.add(2, 0, -1);
  coo.add(2, 2, 0.5);
  const A = coo.toCsr();
  assert.equal(A.n, 3);
  assert.ok(A.rowPtr instanceof Int32Array && A.colIdx instanceof Int32Array && A.val instanceof Float64Array);
  assert.deepEqual(Array.from(A.rowPtr), [0, 2, 2, 4]);
  assert.deepEqual(Array.from(A.colIdx), [0, 1, 0, 2]);
  assert.deepEqual(Array.from(A.val), [2, 4, -1, 5.5]);
  assert.throws(() => coo.add(3, 0, 1), RangeError);
  assert.throws(() => coo.add(0, -1, 1), RangeError);
  // 5000 entries force several reallocations
  const big = new CooBuilder(50);
  for (let k = 0; k < 5000; k++) big.add(k % 50, (k * 7) % 50, 1);
  const B = big.toCsr();
  let sum = 0;
  for (const v of B.val) sum += v;
  assert.equal(sum, 5000);
  for (let i = 0; i < 50; i++) for (let k = B.rowPtr[i] + 1; k < B.rowPtr[i + 1]; k++) assert.ok(B.colIdx[k] > B.colIdx[k - 1]);
});

test('csrMatVec: agrees with the dense product for a random sparse matrix; out buffer is used', () => {
  const rnd = lcg(17);
  const n = 40;
  const D = zeros(n, n);
  const coo = new CooBuilder(n);
  for (let k = 0; k < 300; k++) {
    const i = Math.floor(rnd() * n);
    const j = Math.floor(rnd() * n);
    const v = rnd() - 0.5;
    D[i][j] += v;
    coo.add(i, j, v);
  }
  const A = coo.toCsr();
  const x = Array.from({ length: n }, () => rnd());
  const out = new Float64Array(n);
  const y = csrMatVec(A, x, out);
  assert.equal(y, out);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += D[i][j] * x[j];
    close(y[i], s, 1e-13);
  }
  assert.throws(() => csrMatVec(A, [1, 2]), RangeError);
});

test('pcg: 2-D Poisson problem (40 × 40) — all preconditioners reach the tolerance; IC(0) and SSOR beat Jacobi', () => {
  const { A, n } = laplace2d(40, 40);
  const b = new Float64Array(n).fill(1);
  const its = {};
  for (const precond of ['none', 'jacobi', 'ssor', 'ic0']) {
    const r = pcg(A, b, { precond, tol: 1e-10 });
    assert.ok(r.converged, `${precond} converged`);
    assert.ok(r.residual <= 1e-10);
    assert.ok(residualInf(A, r.x, b) < 1e-8, `${precond} true residual`);
    its[precond] = r.iterations;
  }
  assert.equal(its.none, its.jacobi, 'constant diagonal: Jacobi equals plain CG');
  assert.ok(its.ssor < 0.7 * its.jacobi, `ssor ${its.ssor} vs jacobi ${its.jacobi}`);
  assert.ok(its.ic0 < 0.6 * its.jacobi, `ic0 ${its.ic0} vs jacobi ${its.jacobi}`);
  // CG on an n×n grid converges in O(n_grid) iterations
  assert.ok(its.jacobi < 150);
});

test('pcg: discrete Poisson solution with a manufactured right-hand side (eigenvector sin·sin)', () => {
  // u_ij = sin(π i/(N+1)) sin(π j/(N+1)) is an eigenvector of the 5-point Laplacian with λ = 4 − 4 cos(π/(N+1))
  const N = 30;
  const { A, n } = laplace2d(N, N);
  const u = new Float64Array(n);
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) u[j * N + i] = Math.sin((Math.PI * (i + 1)) / (N + 1)) * Math.sin((Math.PI * (j + 1)) / (N + 1));
  const lam = 4 - 4 * Math.cos(Math.PI / (N + 1));
  const b = u.map((v) => lam * v);
  const r = pcg(A, b, { precond: 'ic0', tol: 1e-12 });
  assert.ok(r.converged);
  for (let k = 0; k < n; k += 37) close(r.x[k], u[k], 1e-9);
});

test('pcg: random sparse SPD system with strongly varying diagonal (seeded LCG) — Jacobi scaling pays off', () => {
  const rnd = lcg(321);
  const n = 400;
  const coo = new CooBuilder(n);
  const diag = new Float64Array(n);
  for (let k = 0; k < 1600; k++) {
    const i = Math.floor(rnd() * n);
    const j = Math.floor(rnd() * n);
    if (i === j) continue;
    const scale = 10 ** (3 * ((i + j) / (2 * n)));
    const v = -rnd() * scale;
    coo.add(i, j, v);
    coo.add(j, i, v);
    diag[i] -= v;
    diag[j] -= v;
  }
  for (let i = 0; i < n; i++) coo.add(i, i, diag[i] * 1.05 + 10 ** ((3 * i) / n));
  const A = coo.toCsr();
  const xt = Float64Array.from({ length: n }, () => rnd() - 0.5);
  const b = csrMatVec(A, xt);
  const jac = pcg(A, b, { precond: 'jacobi', tol: 1e-12 });
  const none = pcg(A, b, { precond: 'none', tol: 1e-12 });
  const ic = pcg(A, b, { precond: 'ic0', tol: 1e-12 });
  const ss = pcg(A, b, { precond: 'ssor', tol: 1e-12, omega: 1.2 });
  for (const r of [jac, none, ic, ss]) {
    assert.ok(r.converged);
    for (let i = 0; i < n; i += 13) close(r.x[i], xt[i], 1e-8);
  }
  assert.ok(jac.iterations < none.iterations, `jacobi ${jac.iterations} vs none ${none.iterations}`);
  assert.ok(ic.iterations <= jac.iterations);
});

test('pcg: warm start, zero right-hand side, iteration cap, custom preconditioner, indefinite matrix', () => {
  const { A, n } = laplace2d(20, 20);
  const b = new Float64Array(n).fill(1);
  const ref = pcg(A, b, { tol: 1e-12 });
  const warm = pcg(A, b, { x0: ref.x, tol: 1e-10 });
  assert.equal(warm.iterations, 0);
  assert.ok(warm.converged);
  const zero = pcg(A, new Float64Array(n));
  assert.deepEqual(Array.from(zero.x), new Array(n).fill(0));
  assert.ok(zero.converged);
  const capped = pcg(A, b, { maxIter: 3 });
  assert.equal(capped.converged, false);
  assert.equal(capped.iterations, 3);
  const custom = pcg(A, b, { precond: (r, z) => { for (let i = 0; i < n; i++) z[i] = r[i] / 4; } });
  assert.equal(custom.iterations, pcg(A, b, { precond: 'jacobi' }).iterations);
  // negative diagonal: not SPD
  const coo = new CooBuilder(2);
  coo.add(0, 0, 1);
  coo.add(1, 1, -1);
  assert.throws(() => pcg(coo.toCsr(), [1, 1]), RangeError);
  assert.throws(() => pcg(A, b, { precond: 'ilu' }), RangeError);
});

test('SkylineSym: random SPD profile matrix (seeded LCG) agrees with dense LU; mulVec; b is not modified', () => {
  const rnd = lcg(55);
  const n = 60;
  const heights = new Int32Array(n);
  for (let j = 0; j < n; j++) heights[j] = Math.min(j, Math.floor(rnd() * 9));
  const K = new SkylineSym(n, heights);
  const D = zeros(n, n);
  for (let j = 0; j < n; j++) {
    for (let i = j - heights[j]; i < j; i++) {
      const v = rnd() - 0.5;
      K.add(i, j, v);
      D[i][j] = v;
      D[j][i] = v;
    }
    K.add(j, j, 12);
    D[j][j] = 12;
  }
  const b = Float64Array.from({ length: n }, () => rnd() - 0.5);
  const b0 = Float64Array.from(b);
  const y = K.mulVec(b);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += D[i][j] * b[j];
    close(y[i], s, 1e-13);
  }
  close(K.get(3, 10), D[3][10], 0);
  close(K.get(10, 3), D[3][10], 0);
  const x = K.solve(b);
  assert.ok(x instanceof Float64Array);
  assert.deepEqual(b, b0);
  const xd = solveDense(D, Array.from(b));
  for (let i = 0; i < n; i++) close(x[i], xd[i], 1e-13);
  assert.equal(K.negativePivots, 0);
  // second right-hand side reuses the factor
  const x2 = K.solve(y);
  for (let i = 0; i < n; i++) close(x2[i], b[i], 1e-12);
});

test('SkylineSym: full-matrix assembly loops add each coefficient once (i > j ignored); profile violations throw', () => {
  const K = new SkylineSym(3, new Int32Array([0, 1, 1]));
  const ke = [[2, -1], [-1, 2]];
  for (const [a, b] of [[0, 1], [1, 2]]) {
    const dofs = [a, b];
    for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) K.add(dofs[p], dofs[q], ke[p][q]);
  }
  assert.equal(K.get(0, 0), 2);
  assert.equal(K.get(1, 1), 4);
  assert.equal(K.get(0, 1), -1);
  assert.equal(K.get(1, 0), -1);
  assert.equal(K.get(0, 2), 0);
  assert.throws(() => K.add(0, 2, 1), RangeError);
  assert.throws(() => K.add(0, 3, 1), RangeError);
  assert.throws(() => new SkylineSym(3, new Int32Array([0, 1])), RangeError);
  K.factor();
  assert.throws(() => K.add(0, 0, 1), RangeError);
});

test('SkylineSym: axial bar under end load — fixed–free chain of springs, u_i = P·i/k exactly', () => {
  const n = 50;
  const k = 2.1e8; // EA/L per element [N/m]
  const P = 1.5e5;
  const K = new SkylineSym(n + 1, new Int32Array(n + 1).fill(1));
  for (let e = 0; e < n; e++) {
    K.add(e, e, k);
    K.add(e, e + 1, -k);
    K.add(e + 1, e + 1, k);
  }
  const f = new Float64Array(n + 1);
  f[n] = P;
  K.applyBc(0);
  const u = K.solve(f);
  for (let i = 0; i <= n; i++) close(u[i], (P * i) / k, 1e-15 * n);
  assert.equal(u[0], 0);
});

test('SkylineSym: prescribed non-zero displacement via applyBc(dof, value, rhs) — linear field between supports', () => {
  const n = 20;
  const K = new SkylineSym(n + 1, new Int32Array(n + 1).fill(1));
  for (let e = 0; e < n; e++) {
    K.add(e, e, 3);
    K.add(e, e + 1, -3);
    K.add(e + 1, e + 1, 3);
  }
  const rhs = new Float64Array(n + 1);
  K.applyBc(0, 0, rhs);
  K.applyBc(n, 0.04, rhs);
  const u = K.solve(rhs);
  for (let i = 0; i <= n; i++) close(u[i], (0.04 * i) / n, 1e-15);
});

test('SkylineSym: unconstrained structure throws a zero-pivot error naming the dof', () => {
  const n = 6;
  const K = new SkylineSym(n, new Int32Array(n).fill(1));
  for (let e = 0; e < n - 1; e++) {
    K.add(e, e, 1);
    K.add(e, e + 1, -1);
    K.add(e + 1, e + 1, 1);
  }
  assert.throws(() => K.factor(), (err) => {
    assert.match(err.message, /zero pivot at dof 5/);
    assert.equal(err.dof, 5);
    return true;
  });
  // an isolated massless dof in the middle
  const Z = new SkylineSym(3, new Int32Array(3));
  Z.add(0, 0, 1);
  Z.add(2, 2, 1);
  assert.throws(() => Z.factor(), /dof 1/);
});

test('SkylineSym: indefinite shift K − σI — number of negative pivots equals the Sturm count of eigenvalues below σ', () => {
  // tridiag(−1, 2, −1), n = 20: λ_k = 2 − 2 cos(kπ/21)
  const n = 20;
  for (const sigma of [0.01, 0.5, 1.7, 3.3, 3.99]) {
    const K = new SkylineSym(n, new Int32Array(n).fill(1));
    for (let i = 0; i < n; i++) {
      K.add(i, i, 2 - sigma);
      if (i > 0) K.add(i - 1, i, -1);
    }
    K.factor();
    let below = 0;
    for (let k = 1; k <= n; k++) if (2 - 2 * Math.cos((k * Math.PI) / (n + 1)) < sigma) below++;
    assert.equal(K.negativePivots, below, `sigma = ${sigma}`);
  }
});

test('reverseCuthillMcKee: valid permutation; restores a narrow band for a randomly renumbered 2-D grid', () => {
  const nx = 12;
  const ny = 40;
  const { edges, n } = laplace2d(nx, ny);
  // random renumbering (Fisher–Yates with the LCG)
  const rnd = lcg(8);
  const shuffle = Int32Array.from({ length: n }, (_, i) => i);
  for (let i = n - 1; i > 0; i--) {
    const j = Math.floor(rnd() * (i + 1));
    const t = shuffle[i];
    shuffle[i] = shuffle[j];
    shuffle[j] = t;
  }
  const scrambled = edges.map(([a, b]) => [shuffle[a], shuffle[b]]);
  const bandwidth = (h) => Math.max(...h);
  const before = skylineHeights(n, scrambled);
  assert.ok(bandwidth(before) > n / 2, 'scrambled numbering has a huge bandwidth');
  const perm = reverseCuthillMcKee(n, scrambled);
  assert.ok(perm instanceof Int32Array);
  assert.deepEqual(Array.from(perm).sort((a, b) => a - b), Array.from({ length: n }, (_, i) => i));
  const after = skylineHeights(n, scrambled, perm);
  // the natural optimum is nx = 12; RCM levels on a grid are diagonals, at most ~nx + 1 wide
  assert.ok(bandwidth(after) <= 2 * nx, `bandwidth after RCM ${bandwidth(after)}`);
  const profile = (h) => h.reduce((s, v) => s + v, 0);
  assert.ok(profile(after) < profile(before) / 10);
});

test('reverseCuthillMcKee + SkylineSym: reordered solve equals the PCG solution of the original system', () => {
  const { A, edges, n } = laplace2d(25, 18);
  const perm = reverseCuthillMcKee(n, edges);
  const inv = new Int32Array(n);
  for (let k = 0; k < n; k++) inv[perm[k]] = k;
  const K = new SkylineSym(n, skylineHeights(n, edges, perm));
  for (let i = 0; i < n; i++) for (let k = A.rowPtr[i]; k < A.rowPtr[i + 1]; k++) K.add(inv[i], inv[A.colIdx[k]], A.val[k]);
  const rnd = lcg(4);
  const b = Float64Array.from({ length: n }, () => rnd());
  const bp = new Float64Array(n);
  for (let k = 0; k < n; k++) bp[k] = b[perm[k]];
  const xp = K.solve(bp);
  const x = new Float64Array(n);
  for (let k = 0; k < n; k++) x[perm[k]] = xp[k];
  assert.ok(residualInf(A, x, b) < 1e-12);
  const it = pcg(A, b, { precond: 'ic0', tol: 1e-13 });
  for (let i = 0; i < n; i += 11) close(x[i], it.x[i], 1e-10);
});

test('reverseCuthillMcKee: disconnected graph, isolated nodes, duplicate edges and self-loops; invalid edges throw', () => {
  const perm = reverseCuthillMcKee(7, [[0, 1], [1, 0], [1, 2], [2, 2], [4, 5]]);
  assert.deepEqual(Array.from(perm).sort((a, b) => a - b), [0, 1, 2, 3, 4, 5, 6]);
  const inv = new Int32Array(7);
  perm.forEach((o, k) => { inv[o] = k; });
  // members of a component are numbered contiguously
  assert.equal(Math.abs(inv[4] - inv[5]), 1);
  assert.ok(Math.max(inv[0], inv[1], inv[2]) - Math.min(inv[0], inv[1], inv[2]) === 2);
  assert.deepEqual(Array.from(reverseCuthillMcKee(3, [])).sort(), [0, 1, 2]);
  assert.throws(() => reverseCuthillMcKee(3, [[0, 3]]), RangeError);
  assert.deepEqual(Array.from(skylineHeights(4, [[0, 2], [3, 1]])), [0, 0, 2, 2]);
});

test('SkylineSym performance: 30 000 dof, half-bandwidth 300 (9 million stored coefficients) factorised and solved', () => {
  const n = 30000;
  const bw = 300;
  const heights = new Int32Array(n);
  for (let j = 0; j < n; j++) heights[j] = Math.min(j, bw);
  const K = new SkylineSym(n, heights);
  const rnd = lcg(1);
  for (let j = 0; j < n; j++) {
    K.add(j, j, bw + 2);
    for (let i = Math.max(0, j - bw); i < j; i++) K.add(i, j, -rnd());
  }
  const xt = Float64Array.from({ length: n }, () => rnd());
  const b = K.mulVec(xt);
  const c0 = process.cpuUsage();
  K.factor();
  const cu = process.cpuUsage(c0);
  const tf = (cu.user + cu.system) / 1000;   // processor time: independent of how busy the machine is
  const x = K.solve(b);
  let e = 0;
  for (let i = 0; i < n; i++) e = Math.max(e, Math.abs(x[i] - xt[i]));
  assert.ok(e < 1e-11, `solution error ${e}`);
  // "a few seconds" on a laptop; generous bound so that a loaded CI machine does not fail spuriously
  assert.ok(tf < 30000, `factorisation took ${tf.toFixed(0)} ms`);
});
