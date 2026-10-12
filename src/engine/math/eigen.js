/**
 * Eigenvalue problems: dense symmetric, dense generalised symmetric-definite, large generalised problems by
 * subspace iteration, and eigenvalues of general real matrices.
 *
 * Unit-agnostic. For structural dynamics K φ = λ M φ with K in N/m and M in kg gives λ = ω² in (rad/s)².
 *
 * References: Wilkinson & Reinsch, "Handbook for Automatic Computation II: Linear Algebra" (1971): tred2/tql2
 * (Householder tridiagonalisation + implicit QL), balanc/elmhes/hqr (balancing, elimination to Hessenberg form,
 * shifted QR); Golub & Van Loan, "Matrix Computations", 4th ed., §8.3, §7.5, §8.7.2 (Cholesky reduction of the
 * generalised problem); Bathe, "Finite Element Procedures" (1996) §11.6 (subspace iteration with Rayleigh–Ritz).
 */

import { makeRng } from './rng.js';

const EPS = 2.220446049250313e-16;

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

/** Validate a square symmetric matrix and return its symmetrised flat row-major copy. */
function symFlat(A, name) {
  if (!Array.isArray(A) || A.length < 1) throw new RangeError(`${name} must be a non-empty square matrix`);
  const n = A.length;
  const a = new Float64Array(n * n);
  let amax = 0;
  for (let i = 0; i < n; i++) {
    if (!isVec(A[i]) || A[i].length !== n) throw new RangeError(`${name} must be square (row ${i} has wrong length)`);
    for (let j = 0; j < n; j++) {
      const v = A[i][j];
      if (!Number.isFinite(v)) throw new RangeError(`${name} contains a non-finite entry at (${i}, ${j})`);
      a[i * n + j] = v;
      if (Math.abs(v) > amax) amax = Math.abs(v);
    }
  }
  const tol = 1e-8 * amax;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const u = a[i * n + j];
      const l = a[j * n + i];
      if (Math.abs(u - l) > tol) throw new RangeError(`${name} must be symmetric (entries (${i}, ${j}) and (${j}, ${i}) differ)`);
      const m = 0.5 * (u + l);
      a[i * n + j] = m;
      a[j * n + i] = m;
    }
  }
  return a;
}

/**
 * Symmetric eigen-decomposition of a flat row-major matrix (overwritten with the eigenvectors as columns).
 * Householder tridiagonalisation (tred2) followed by the implicit QL algorithm (tql2).
 * Returns the eigenvalues in ascending order; column k of V is the unit eigenvector of d[k].
 * Internally the transformation matrix is held transposed (eigenvectors as rows): the Householder accumulation
 * and the QL plane rotations then run along contiguous memory instead of striding by n, which is what governs the
 * run time once the matrix no longer fits the first-level cache.
 */
function symEigFlat(V, n) {
  const d = new Float64Array(n);
  const e = new Float64Array(n);
  // --- tred2 ---
  for (let j = 0; j < n; j++) d[j] = V[j * n + n - 1];
  for (let i = n - 1; i > 0; i--) {
    let scale = 0;
    let h = 0;
    for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
    if (scale === 0) {
      e[i] = d[i - 1];
      for (let j = 0; j < i; j++) {
        d[j] = V[j * n + i - 1];
        V[j * n + i] = 0;
        V[i * n + j] = 0;
      }
    } else {
      for (let k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      let f = d[i - 1];
      let g = Math.sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (let j = 0; j < i; j++) e[j] = 0;
      for (let j = 0; j < i; j++) {
        f = d[j];
        V[i * n + j] = f;
        g = e[j] + V[j * n + j] * f;
        for (let k = j + 1; k <= i - 1; k++) {
          g += V[j * n + k] * d[k];
          e[k] += V[j * n + k] * f;
        }
        e[j] = g;
      }
      f = 0;
      for (let j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      const hh = f / (h + h);
      for (let j = 0; j < i; j++) e[j] -= hh * d[j];
      for (let j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (let k = j; k <= i - 1; k++) V[j * n + k] -= f * e[k] + g * d[k];
        d[j] = V[j * n + i - 1];
        V[j * n + i] = 0;
      }
    }
    d[i] = h;
  }
  for (let i = 0; i < n - 1; i++) {
    V[i * n + n - 1] = V[i * n + i];
    V[i * n + i] = 1;
    const h = d[i + 1];
    if (h !== 0) {
      for (let k = 0; k <= i; k++) d[k] = V[(i + 1) * n + k] / h;
      for (let j = 0; j <= i; j++) {
        let g = 0;
        for (let k = 0; k <= i; k++) g += V[(i + 1) * n + k] * V[j * n + k];
        for (let k = 0; k <= i; k++) V[j * n + k] -= g * d[k];
      }
    }
    for (let k = 0; k <= i; k++) V[(i + 1) * n + k] = 0;
  }
  for (let j = 0; j < n; j++) {
    d[j] = V[j * n + n - 1];
    V[j * n + n - 1] = 0;
  }
  V[(n - 1) * n + n - 1] = 1;
  e[0] = 0;
  // --- tql2 ---
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  let f = 0;
  let tst1 = 0;
  for (let l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    let m = l;
    while (m < n - 1) {
      if (Math.abs(e[m]) <= EPS * tst1) break;
      m++;
    }
    if (m > l) {
      let iter = 0;
      do {
        if (++iter > 200) throw new Error('eigSym: QL iteration did not converge');
        let g = d[l];
        let p = (d[l + 1] - g) / (2 * e[l]);
        let r = Math.hypot(p, 1);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const dl1 = d[l + 1];
        let h = g - d[l];
        for (let i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        p = d[m];
        let c = 1;
        let c2 = c;
        let c3 = c;
        const el1 = e[l + 1];
        let s = 0;
        let s2 = 0;
        for (let i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = Math.hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          const r0 = i * n;
          const r1 = r0 + n;
          for (let k = 0; k < n; k++) {
            h = V[r1 + k];
            V[r1 + k] = s * V[r0 + k] + c * h;
            V[r0 + k] = c * V[r0 + k] - s * h;
          }
        }
        p = (-s * s2 * c3 * el1 * e[l]) / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (Math.abs(e[l]) > EPS * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }
  // --- sort ascending (selection sort with column swaps) ---
  for (let i = 0; i < n - 1; i++) {
    let k = i;
    let p = d[i];
    for (let j = i + 1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k !== i) {
      d[k] = d[i];
      d[i] = p;
      for (let r = 0; r < n; r++) {
        const t = V[i * n + r];
        V[i * n + r] = V[k * n + r];
        V[k * n + r] = t;
      }
    }
  }
  // back to the documented layout: eigenvectors as columns
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const t = V[i * n + j];
      V[i * n + j] = V[j * n + i];
      V[j * n + i] = t;
    }
  }
  return d;
}

/** Copy column k of a flat matrix (rows × cols) into a plain array, with the largest component made positive. */
function columnOf(V, rows, cols, k) {
  const v = new Array(rows);
  let big = 0;
  for (let i = 0; i < rows; i++) {
    const x = V[i * cols + k];
    v[i] = x;
    if (Math.abs(x) > Math.abs(big)) big = x;
  }
  if (big < 0) for (let i = 0; i < rows; i++) v[i] = -v[i];
  return v;
}

/**
 * Eigenvalues and eigenvectors of a dense real symmetric matrix: Householder reduction to tridiagonal form and
 * implicit QL iterations (Wilkinson & Reinsch tred2/tql2). Cost ≈ 9n³ flops; eigenvalues accurate to ~eps·‖A‖.
 * Validity: A symmetric (asymmetry above 1e-8·max|A| throws RangeError); intended for n up to ~1500.
 * @param {number[][]} A n×n symmetric matrix
 * @returns {{ values: number[], vectors: number[][] }} values ascending; vectors[k] is the unit-norm eigenvector
 *   of values[k] (orthonormal set; sign fixed so that the largest component is positive)
 */
export function eigSym(A) {
  const a = symFlat(A, 'eigSym: A');
  const n = A.length;
  const d = symEigFlat(a, n);
  const vectors = new Array(n);
  for (let k = 0; k < n; k++) vectors[k] = columnOf(a, n, n, k);
  return { values: Array.from(d), vectors };
}

/**
 * Generalised symmetric-definite eigenproblem K φ = λ M φ for dense matrices, by Cholesky reduction
 * M = L Lᵀ, C = L⁻¹ K L⁻ᵀ, C y = λ y, φ = L⁻ᵀ y (GVL §8.7.2), then eigSym.
 * Validity: K symmetric, M symmetric positive definite (throws RangeError otherwise; massless degrees of freedom
 * must be condensed out by the caller). M may also be given as a vector holding a diagonal (lumped) mass matrix.
 * @param {number[][]} K n×n symmetric matrix (stiffness)
 * @param {number[][] | ArrayLike<number>} M n×n SPD matrix (mass), or its diagonal as a vector of length n
 * @param {{ nModes?: number }} [opts] nModes: number of lowest modes returned (default all n)
 * @returns {{ values: number[], vectors: number[][] }} values ascending (λ = ω² for K/M in SI); vectors[k] is the
 *   mode shape of values[k], M-orthonormal (φᵢᵀ M φⱼ = δᵢⱼ)
 */
export function eigGenSym(K, M, { nModes } = {}) {
  const k = symFlat(K, 'eigGenSym: K');
  const n = K.length;
  const p = nModes === undefined ? n : nModes;
  if (!Number.isInteger(p) || p < 1 || p > n) throw new RangeError('eigGenSym: nModes must be an integer in 1..n');
  // Cholesky factor of M (flat, lower)
  const L = new Float64Array(n * n);
  const diagonalM = isVec(M) && M.length === n && typeof M[0] === 'number';
  if (diagonalM) {
    for (let i = 0; i < n; i++) {
      if (!(M[i] > 0) || !Number.isFinite(M[i])) throw new RangeError(`eigGenSym: M must be positive definite (diagonal entry ${i})`);
      L[i * n + i] = Math.sqrt(M[i]);
    }
  } else {
    const m = symFlat(M, 'eigGenSym: M');
    if (M.length !== n) throw new RangeError('eigGenSym: K and M must have the same order');
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let s = m[i * n + j];
        for (let q = 0; q < j; q++) s -= L[i * n + q] * L[j * n + q];
        if (i === j) {
          if (!(s > 0) || !Number.isFinite(s)) throw new RangeError(`eigGenSym: M must be positive definite (pivot ${s} at row ${i})`);
          L[i * n + i] = Math.sqrt(s);
        } else {
          L[i * n + j] = s / L[j * n + j];
        }
      }
    }
  }
  // B ← L⁻¹ B, row-oriented forward substitution on all columns at once
  const forward = (B) => {
    for (let i = 0; i < n; i++) {
      const ri = i * n;
      for (let q = 0; q < i; q++) {
        const l = L[ri + q];
        if (l === 0) continue;
        const rq = q * n;
        for (let c = 0; c < n; c++) B[ri + c] -= l * B[rq + c];
      }
      const inv = 1 / L[ri + i];
      for (let c = 0; c < n; c++) B[ri + c] *= inv;
    }
  };
  forward(k); // L⁻¹ K
  // transpose in place, then L⁻¹ (L⁻¹ K)ᵀ = L⁻¹ K L⁻ᵀ (symmetric)
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const t = k[i * n + j];
      k[i * n + j] = k[j * n + i];
      k[j * n + i] = t;
    }
  }
  forward(k);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const mid = 0.5 * (k[i * n + j] + k[j * n + i]);
      k[i * n + j] = mid;
      k[j * n + i] = mid;
    }
  }
  const d = symEigFlat(k, n);
  // φ = L⁻ᵀ y for the requested columns (back substitution)
  const vectors = new Array(p);
  const y = new Float64Array(n);
  for (let c = 0; c < p; c++) {
    for (let i = 0; i < n; i++) y[i] = k[i * n + c];
    let big = 0;
    for (let i = n - 1; i >= 0; i--) {
      let s = y[i];
      for (let q = i + 1; q < n; q++) s -= L[q * n + i] * y[q];
      y[i] = s / L[i * n + i];
      if (Math.abs(y[i]) > Math.abs(big)) big = y[i];
    }
    const v = new Array(n);
    const sg = big < 0 ? -1 : 1;
    for (let i = 0; i < n; i++) v[i] = sg * y[i];
    vectors[c] = v;
  }
  return { values: Array.from(d.subarray(0, p)), vectors };
}

/**
 * Reduced generalised problem Kr z = λ Mr z of order q with possibly rank-deficient Mr: Mr = W diag(μ) Wᵀ, keep
 * μ > 1e-12 μ_max, T = W μ^(−½), eigen-decompose Tᵀ Kr T. Returns { r, values, Q } with Q (q×r, flat) such that
 * Qᵀ Mr Q = I and Qᵀ Kr Q = diag(values).
 */
function reducedEig(Kr, Mr, q) {
  const W = Float64Array.from(Mr);
  const mu = symEigFlat(W, q);
  const mmax = mu[q - 1];
  if (!(mmax > 0)) return { r: 0, values: new Float64Array(0), Q: new Float64Array(0) };
  let start = 0;
  while (start < q && !(mu[start] > 1e-12 * mmax)) start++;
  const r = q - start;
  const T = new Float64Array(q * r);
  for (let c = 0; c < r; c++) {
    const s = 1 / Math.sqrt(mu[start + c]);
    for (let i = 0; i < q; i++) T[i * r + c] = W[i * q + start + c] * s;
  }
  // G = Tᵀ Kr T
  const KT = new Float64Array(q * r);
  for (let i = 0; i < q; i++) {
    for (let c = 0; c < r; c++) {
      let s = 0;
      for (let j = 0; j < q; j++) s += Kr[i * q + j] * T[j * r + c];
      KT[i * r + c] = s;
    }
  }
  const G = new Float64Array(r * r);
  for (let a = 0; a < r; a++) {
    for (let b = a; b < r; b++) {
      let s = 0;
      for (let i = 0; i < q; i++) s += T[i * r + a] * KT[i * r + b];
      G[a * r + b] = s;
      G[b * r + a] = s;
    }
  }
  const values = symEigFlat(G, r);
  const Q = new Float64Array(q * r);
  for (let i = 0; i < q; i++) {
    for (let c = 0; c < r; c++) {
      let s = 0;
      for (let a = 0; a < r; a++) s += T[i * r + a] * G[a * r + c];
      Q[i * r + c] = s;
    }
  }
  return { r, values, Q };
}

/**
 * Lowest eigenpairs of the large generalised problem K φ = λ M φ by subspace (simultaneous inverse) iteration
 * with a Rayleigh–Ritz projection in every sweep (Bathe 1996 §11.6). The caller supplies the action of K⁻¹
 * (e.g. SkylineSym.solve) and of M. The subspace starts with q = min(n, max(2·nModes, nModes + 8)) vectors;
 * convergence is declared when every requested eigenvalue changes by less than tol (relative) between two sweeps.
 *
 * Clustered spectra. The error of mode k shrinks by (λ_k/λ_(q+1))² per sweep, so a group of nearly equal
 * eigenvalues that reaches beyond the subspace (repeated members, symmetric structures) stalls the iteration:
 * the Ritz values inside the group cannot be told apart until the subspace holds the whole group. After every
 * sweep the rate is estimated from the Ritz values as θ_nModes/θ_q; when it exceeds `growRatio` (default 0.7,
 * i.e. fewer than a factor two gained per sweep) the subspace is doubled with fresh random guard vectors, up to
 * `maxSubspace` (default min(n, max(8·nModes, nModes + 64))). A group wider than that limit is not resolved:
 * `converged` is then false and the values are accurate to about the width of the group only — raise
 * `maxSubspace` or use a dense solver.
 * Validity: K symmetric positive definite (shift K by σM first if rigid-body modes are present, then subtract σ),
 * M symmetric positive semi-definite with rank >= nModes. Does not throw on non-convergence: check `converged`.
 * @param {{ n: number, nModes: number, solveK: (b: Float64Array) => ArrayLike<number>,
 *   mulM: (x: Float64Array) => ArrayLike<number>, tol?: number, maxIter?: number, seed?: number,
 *   maxSubspace?: number, growRatio?: number }} p
 *   n: problem order; nModes: number of lowest modes wanted; solveK(b): returns K⁻¹ b; mulM(x): returns M x
 *   (both may reuse their output buffer; the arguments must not be modified); tol: relative eigenvalue tolerance
 *   (default 1e-8); maxIter: sweep cap (default 60); seed: seed of the random starting vectors (default 1);
 *   maxSubspace: largest number of iteration vectors (default min(n, max(8·nModes, nModes + 64))); growRatio:
 *   Ritz-value ratio θ_nModes/θ_q above which the subspace is enlarged, in (0, 1] (default 0.7; 1 never enlarges)
 * @returns {{ values: number[], vectors: Float64Array[], iterations: number, converged: boolean,
 *   subspaceSize: number }} values: the nModes lowest λ ascending (ω² for K/M in SI); vectors[k]: mode shape of
 *   values[k], M-orthonormal; iterations: sweeps performed; subspaceSize: iteration vectors in the last sweep
 */
export function subspaceIteration({ n, nModes, solveK, mulM, tol = 1e-8, maxIter = 60, seed = 1, maxSubspace, growRatio = 0.7 } = {}) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('subspaceIteration: n must be an integer >= 1');
  if (!Number.isInteger(nModes) || nModes < 1 || nModes > n) throw new RangeError('subspaceIteration: nModes must be an integer in 1..n');
  if (typeof solveK !== 'function') throw new RangeError('subspaceIteration: solveK must be a function');
  if (typeof mulM !== 'function') throw new RangeError('subspaceIteration: mulM must be a function');
  if (!(tol > 0)) throw new RangeError('subspaceIteration: tol must be > 0');
  if (!Number.isInteger(maxIter) || maxIter < 1) throw new RangeError('subspaceIteration: maxIter must be an integer >= 1');
  const p = nModes;
  let q = Math.min(n, Math.max(2 * p, p + 8));
  if (!(growRatio > 0 && growRatio <= 1)) throw new RangeError('subspaceIteration: growRatio must lie in (0, 1]');
  if (maxSubspace !== undefined && (!Number.isInteger(maxSubspace) || maxSubspace < p)) {
    throw new RangeError('subspaceIteration: maxSubspace must be an integer >= nModes');
  }
  const qMax = Math.min(n, maxSubspace !== undefined ? Math.max(maxSubspace, q) : Math.max(8 * p, p + 64));
  let lastGrow = 0;
  const rng = makeRng(seed);
  const apply = (fn, x, name) => {
    const y = fn(x);
    if (!y || y.length !== n) throw new RangeError(`subspaceIteration: ${name} must return a vector of length n`);
    return Float64Array.from(y);
  };
  // starting vectors: one uniform vector (excites the smooth modes) and random vectors
  let X = new Array(q);
  for (let j = 0; j < q; j++) {
    const x = new Float64Array(n);
    if (j === 0) x.fill(1);
    else for (let i = 0; i < n; i++) x[i] = 2 * rng.next() - 1;
    X[j] = x;
  }
  let Y = X.map((x) => apply(mulM, x, 'mulM'));
  let prev = null;
  let values = null;
  let converged = false;
  let iterations = 0;
  while (iterations < maxIter) {
    iterations++;
    const Xb = new Array(q);
    const Z = new Array(q);
    for (let j = 0; j < q; j++) {
      Xb[j] = apply(solveK, Y[j], 'solveK');
      Z[j] = apply(mulM, Xb[j], 'mulM');
    }
    const Kr = new Float64Array(q * q);
    const Mr = new Float64Array(q * q);
    for (let i = 0; i < q; i++) {
      const xi = Xb[i];
      for (let j = i; j < q; j++) {
        const yj = Y[j];
        const yi = Y[i];
        const xj = Xb[j];
        const zj = Z[j];
        let sk = 0;
        let sk2 = 0;
        let sm = 0;
        for (let r = 0; r < n; r++) {
          sk += xi[r] * yj[r];
          sk2 += xj[r] * yi[r];
          sm += xi[r] * zj[r];
        }
        if (!Number.isFinite(sk) || !Number.isFinite(sm)) throw new Error('subspaceIteration: non-finite projection (check solveK and mulM)');
        Kr[i * q + j] = 0.5 * (sk + sk2);
        Kr[j * q + i] = 0.5 * (sk + sk2);
        Mr[i * q + j] = sm;
        Mr[j * q + i] = sm;
      }
    }
    const red = reducedEig(Kr, Mr, q);
    const r = red.r;
    if (r < p) {
      throw new Error(`subspaceIteration: the iteration subspace has M-rank ${r} < nModes = ${p} (M must have rank >= nModes)`);
    }
    // Ritz vectors X = Xb Q and their M-images Y = Z Q
    const Xn = new Array(r);
    const Yn = new Array(r);
    for (let c = 0; c < r; c++) {
      const x = new Float64Array(n);
      const y = new Float64Array(n);
      for (let i = 0; i < q; i++) {
        const w = red.Q[i * r + c];
        if (w === 0) continue;
        const xb = Xb[i];
        const z = Z[i];
        for (let k = 0; k < n; k++) {
          x[k] += w * xb[k];
          y[k] += w * z[k];
        }
      }
      Xn[c] = x;
      Yn[c] = y;
    }
    X = Xn;
    Y = Yn;
    q = r;
    values = red.values;
    if (prev && prev.length >= p) {
      let ok = true;
      for (let k = 0; k < p; k++) {
        if (!(Math.abs(values[k] - prev[k]) <= tol * Math.abs(values[k]))) {
          ok = false;
          break;
        }
      }
      if (ok) {
        converged = true;
        break;
      }
    }
    prev = values;
    // Slow rate (a cluster reaching beyond the subspace): add random guard vectors. Not judged in the first sweep
    // or right after an enlargement, while the top Ritz values are still far from the eigenvalues they bound.
    if (q < qMax && iterations >= 2 && iterations - lastGrow >= 2 && values[q - 1] > 0
      && Math.abs(values[p - 1]) > growRatio * Math.abs(values[q - 1])) {
      const add = Math.min(q, qMax - q);
      for (let j = 0; j < add; j++) {
        const x = new Float64Array(n);
        for (let i = 0; i < n; i++) x[i] = 2 * rng.next() - 1;
        X.push(x);
        Y.push(apply(mulM, x, 'mulM'));
      }
      q += add;
      lastGrow = iterations;
    }
  }
  const vectors = new Array(p);
  for (let k = 0; k < p; k++) {
    const v = X[k];
    let big = 0;
    for (let i = 0; i < n; i++) if (Math.abs(v[i]) > Math.abs(big)) big = v[i];
    if (big < 0) for (let i = 0; i < n; i++) v[i] = -v[i];
    vectors[k] = v;
  }
  return { values: Array.from(values.subarray(0, p)), vectors, iterations, converged, subspaceSize: X.length };
}

/**
 * Eigenvalues of a general (unsymmetric) real matrix: balancing, Gaussian elimination with pivoting to upper
 * Hessenberg form and the Francis double-shift QR algorithm (Wilkinson & Reinsch balanc/elmhes/hqr).
 * Complex eigenvalues come in conjugate pairs. Eigenvectors are not computed.
 * Validity: square real matrix; throws Error if the QR iteration fails to converge (rare; 60 iterations per
 * eigenvalue with exceptional shifts).
 * @param {number[][]} A n×n real matrix
 * @returns {{ re: number[], im: number[] }} real and imaginary parts, sorted by real part then imaginary part
 */
export function eigGeneral(A) {
  if (!Array.isArray(A) || A.length < 1) throw new RangeError('eigGeneral: A must be a non-empty square matrix');
  const n = A.length;
  const a = new Array(n);
  for (let i = 0; i < n; i++) {
    if (!isVec(A[i]) || A[i].length !== n) throw new RangeError(`eigGeneral: A must be square (row ${i} has wrong length)`);
    a[i] = Float64Array.from(A[i]);
    for (let j = 0; j < n; j++) if (!Number.isFinite(a[i][j])) throw new RangeError(`eigGeneral: A contains a non-finite entry at (${i}, ${j})`);
  }
  // --- balance ---
  const RADIX = 2;
  const sqrdx = RADIX * RADIX;
  let done = false;
  while (!done) {
    done = true;
    for (let i = 0; i < n; i++) {
      let r = 0;
      let c = 0;
      for (let j = 0; j < n; j++) {
        if (j !== i) {
          c += Math.abs(a[j][i]);
          r += Math.abs(a[i][j]);
        }
      }
      if (c !== 0 && r !== 0) {
        let g = r / RADIX;
        let f = 1;
        const s = c + r;
        while (c < g) {
          f *= RADIX;
          c *= sqrdx;
        }
        g = r * RADIX;
        while (c > g) {
          f /= RADIX;
          c /= sqrdx;
        }
        if ((c + r) / f < 0.95 * s) {
          done = false;
          g = 1 / f;
          for (let j = 0; j < n; j++) a[i][j] *= g;
          for (let j = 0; j < n; j++) a[j][i] *= f;
        }
      }
    }
  }
  // --- reduction to Hessenberg form by stabilised elementary similarity transformations ---
  for (let m = 1; m < n - 1; m++) {
    let x = 0;
    let i = m;
    for (let j = m; j < n; j++) {
      if (Math.abs(a[j][m - 1]) > Math.abs(x)) {
        x = a[j][m - 1];
        i = j;
      }
    }
    if (i !== m) {
      for (let j = m - 1; j < n; j++) {
        const t = a[i][j];
        a[i][j] = a[m][j];
        a[m][j] = t;
      }
      for (let j = 0; j < n; j++) {
        const t = a[j][i];
        a[j][i] = a[j][m];
        a[j][m] = t;
      }
    }
    if (x !== 0) {
      for (let r = m + 1; r < n; r++) {
        let y = a[r][m - 1];
        if (y !== 0) {
          y /= x;
          a[r][m - 1] = y;
          for (let j = m; j < n; j++) a[r][j] -= y * a[m][j];
          for (let j = 0; j < n; j++) a[j][m] += y * a[j][r];
        }
      }
    }
  }
  for (let i = 2; i < n; i++) for (let j = 0; j < i - 1; j++) a[i][j] = 0;
  // --- shifted QR on the Hessenberg matrix (eigenvalues only) ---
  const wr = new Float64Array(n);
  const wi = new Float64Array(n);
  let anorm = 0;
  for (let i = 0; i < n; i++) for (let j = Math.max(i - 1, 0); j < n; j++) anorm += Math.abs(a[i][j]);
  let nn = n - 1;
  let t = 0;
  let p = 0;
  let q = 0;
  let r = 0;
  let s = 0;
  let w = 0;
  let x = 0;
  let y = 0;
  let z = 0;
  while (nn >= 0) {
    let its = 0;
    let l;
    do {
      for (l = nn; l >= 1; l--) {
        s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
        if (s === 0) s = anorm;
        if (Math.abs(a[l][l - 1]) <= EPS * s) {
          a[l][l - 1] = 0;
          break;
        }
      }
      x = a[nn][nn];
      if (l === nn) {
        // one real root
        wr[nn] = x + t;
        wi[nn--] = 0;
      } else {
        y = a[nn - 1][nn - 1];
        w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l === nn - 1) {
          // a 2×2 block: two roots
          p = 0.5 * (y - x);
          q = p * p + w;
          z = Math.sqrt(Math.abs(q));
          x += t;
          if (q >= 0) {
            z = p + (p >= 0 ? Math.abs(z) : -Math.abs(z));
            wr[nn - 1] = x + z;
            wr[nn] = wr[nn - 1];
            if (z !== 0) wr[nn] = x - w / z;
            wi[nn - 1] = 0;
            wi[nn] = 0;
          } else {
            wr[nn - 1] = x + p;
            wr[nn] = x + p;
            wi[nn - 1] = z;
            wi[nn] = -z;
          }
          nn -= 2;
        } else {
          if (its >= 60) throw new Error('eigGeneral: QR iteration did not converge');
          if (its > 0 && its % 10 === 0) {
            // exceptional shift
            t += x;
            for (let i = 0; i <= nn; i++) a[i][i] -= x;
            s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
            x = 0.75 * s;
            y = x;
            w = -0.4375 * s * s;
          }
          ++its;
          let m;
          for (m = nn - 2; m >= l; m--) {
            z = a[m][m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
            q = a[m + 1][m + 1] - z - r - s;
            r = a[m + 2][m + 1];
            s = Math.abs(p) + Math.abs(q) + Math.abs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m === l) break;
            const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
            const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
            if (u <= EPS * v) break;
          }
          for (let i = m + 2; i <= nn; i++) {
            a[i][i - 2] = 0;
            if (i !== m + 2) a[i][i - 3] = 0;
          }
          for (let k = m; k <= nn - 1; k++) {
            if (k !== m) {
              p = a[k][k - 1];
              q = a[k + 1][k - 1];
              r = 0;
              if (k !== nn - 1) r = a[k + 2][k - 1];
              x = Math.abs(p) + Math.abs(q) + Math.abs(r);
              if (x !== 0) {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            const nrm = Math.sqrt(p * p + q * q + r * r);
            s = p >= 0 ? nrm : -nrm;
            if (s !== 0) {
              if (k === m) {
                if (l !== m) a[k][k - 1] = -a[k][k - 1];
              } else {
                a[k][k - 1] = -s * x;
              }
              p += s;
              x = p / s;
              y = q / s;
              z = r / s;
              q /= p;
              r /= p;
              for (let j = k; j <= nn; j++) {
                p = a[k][j] + q * a[k + 1][j];
                if (k !== nn - 1) {
                  p += r * a[k + 2][j];
                  a[k + 2][j] -= p * z;
                }
                a[k + 1][j] -= p * y;
                a[k][j] -= p * x;
              }
              const mmin = nn < k + 3 ? nn : k + 3;
              for (let i = l; i <= mmin; i++) {
                p = x * a[i][k] + y * a[i][k + 1];
                if (k !== nn - 1) {
                  p += z * a[i][k + 2];
                  a[i][k + 2] -= p * r;
                }
                a[i][k + 1] -= p * q;
                a[i][k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  const idx = Array.from({ length: n }, (_, i) => i).sort((i, j) => wr[i] - wr[j] || wi[i] - wi[j]);
  return { re: idx.map((i) => wr[i]), im: idx.map((i) => wi[i]) };
}
