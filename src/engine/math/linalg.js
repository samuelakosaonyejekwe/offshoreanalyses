/**
 * Dense and banded linear algebra.
 *
 * Dense matrices are arrays of rows (`number[][]`; rows may also be Float64Array). Vectors may be plain arrays or
 * Float64Array. All quantities are dimensionless as far as this file is concerned: units follow the caller's
 * system (A in [a], b in [b] gives x in [b]/[a]). Inputs are never mutated.
 *
 * References: Golub & Van Loan, "Matrix Computations", 4th ed. (GVL): LU §3.4, Cholesky §4.2, banded Cholesky
 * §4.3.5, Householder QR §5.2, one-sided Jacobi SVD §8.6.3; Thomas algorithm: GVL §4.3.6.
 */

const EPS = 2.220446049250313e-16;

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

function isMat(A) {
  return Array.isArray(A) && A.length > 0 && isVec(A[0]);
}

function checkSquare(A, name) {
  if (!isMat(A)) throw new RangeError(`${name} must be a non-empty matrix (array of rows)`);
  const n = A.length;
  for (let i = 0; i < n; i++) {
    if (!isVec(A[i]) || A[i].length !== n) throw new RangeError(`${name} must be square (row ${i} has wrong length)`);
  }
  return n;
}

function checkRect(A, name) {
  if (!isMat(A)) throw new RangeError(`${name} must be a non-empty matrix (array of rows)`);
  const m = A.length;
  const n = A[0].length;
  for (let i = 0; i < m; i++) {
    if (!isVec(A[i]) || A[i].length !== n) throw new RangeError(`${name} must be rectangular (row ${i} has wrong length)`);
  }
  return [m, n];
}

/**
 * Zero vector or zero matrix.
 * Validity: n >= 0, m >= 0 integers.
 * @param {number} n number of rows (or vector length) [-]
 * @param {number} [m] number of columns; omitted -> vector [-]
 * @returns {number[] | number[][]} zero vector (length n) or n×m zero matrix
 */
export function zeros(n, m) {
  if (!(n >= 0) || !Number.isInteger(n)) throw new RangeError('zeros: n must be a non-negative integer');
  if (m === undefined) return new Array(n).fill(0);
  if (!(m >= 0) || !Number.isInteger(m)) throw new RangeError('zeros: m must be a non-negative integer');
  const A = new Array(n);
  for (let i = 0; i < n; i++) A[i] = new Array(m).fill(0);
  return A;
}

/**
 * Identity matrix.
 * @param {number} n order [-]
 * @returns {number[][]} n×n identity
 */
export function identity(n) {
  const A = zeros(n, n);
  for (let i = 0; i < n; i++) A[i][i] = 1;
  return A;
}

/**
 * Matrix product C = A B.
 * Validity: A is m×p, B is p×n.
 * @param {number[][]} A m×p matrix
 * @param {number[][]} B p×n matrix
 * @returns {number[][]} m×n product
 */
export function matMul(A, B) {
  const [m, p] = checkRect(A, 'matMul: A');
  const [pb, n] = checkRect(B, 'matMul: B');
  if (p !== pb) throw new RangeError('matMul: inner dimensions of A and B differ');
  const C = new Array(m);
  for (let i = 0; i < m; i++) {
    const Ai = A[i];
    const Ci = new Array(n).fill(0);
    for (let k = 0; k < p; k++) {
      const a = Ai[k];
      if (a === 0) continue;
      const Bk = B[k];
      for (let j = 0; j < n; j++) Ci[j] += a * Bk[j];
    }
    C[i] = Ci;
  }
  return C;
}

/**
 * Matrix–vector product y = A x.
 * @param {number[][]} A m×n matrix
 * @param {ArrayLike<number>} x vector, length n
 * @returns {number[]} vector, length m
 */
export function matVec(A, x) {
  const [m, n] = checkRect(A, 'matVec: A');
  if (!isVec(x) || x.length !== n) throw new RangeError('matVec: x length must equal the number of columns of A');
  const y = new Array(m);
  for (let i = 0; i < m; i++) {
    const Ai = A[i];
    let s = 0;
    for (let j = 0; j < n; j++) s += Ai[j] * x[j];
    y[i] = s;
  }
  return y;
}

/**
 * Transpose.
 * @param {number[][]} A m×n matrix
 * @returns {number[][]} n×m matrix
 */
export function transpose(A) {
  const [m, n] = checkRect(A, 'transpose: A');
  const T = new Array(n);
  for (let j = 0; j < n; j++) {
    const Tj = new Array(m);
    for (let i = 0; i < m; i++) Tj[i] = A[i][j];
    T[j] = Tj;
  }
  return T;
}

/**
 * Euclidean inner product.
 * @param {ArrayLike<number>} a vector
 * @param {ArrayLike<number>} b vector of the same length
 * @returns {number} Σ a_i b_i
 */
export function dot(a, b) {
  if (a.length !== b.length) throw new RangeError('dot: a and b must have the same length');
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

/**
 * Euclidean (2-) norm, computed with scaling so that it neither overflows nor underflows.
 * @param {ArrayLike<number>} a vector
 * @returns {number} sqrt(Σ a_i²)
 */
export function norm(a) {
  let scale = 0;
  for (let i = 0; i < a.length; i++) {
    const v = Math.abs(a[i]);
    if (v > scale) scale = v;
  }
  if (scale === 0 || !Number.isFinite(scale)) return scale;
  let s = 0;
  for (let i = 0; i < a.length; i++) {
    const v = a[i] / scale;
    s += v * v;
  }
  return scale * Math.sqrt(s);
}

/**
 * alpha·x + y (new vector; inputs untouched).
 * @param {number} alpha scalar
 * @param {ArrayLike<number>} x vector
 * @param {ArrayLike<number>} y vector of the same length
 * @returns {number[]} alpha·x + y
 */
export function axpy(alpha, x, y) {
  if (x.length !== y.length) throw new RangeError('axpy: x and y must have the same length');
  const r = new Array(x.length);
  for (let i = 0; i < x.length; i++) r[i] = alpha * x[i] + y[i];
  return r;
}

/** LU with partial pivoting; never throws on singularity (flag returned). */
function luRaw(A, name) {
  const n = checkSquare(A, name);
  const lu = new Array(n);
  let amax = 0;
  for (let i = 0; i < n; i++) {
    const r = Float64Array.from(A[i]);
    for (let j = 0; j < n; j++) {
      const v = Math.abs(r[j]);
      if (!(v <= Infinity)) throw new RangeError(`${name} contains NaN at (${i}, ${j})`);
      if (v > amax) amax = v;
    }
    lu[i] = r;
  }
  const piv = new Int32Array(n);
  let sign = 1;
  let singular = !(amax > 0) || !Number.isFinite(amax);
  const tiny = EPS * amax;
  for (let k = 0; k < n; k++) {
    let p = k;
    let pmax = Math.abs(lu[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(lu[i][k]);
      if (v > pmax) {
        pmax = v;
        p = i;
      }
    }
    piv[k] = p;
    if (p !== k) {
      const t = lu[p];
      lu[p] = lu[k];
      lu[k] = t;
      sign = -sign;
    }
    if (!(pmax > tiny)) {
      singular = true;
      if (pmax === 0) continue;
    }
    const rk = lu[k];
    const inv = 1 / rk[k];
    for (let i = k + 1; i < n; i++) {
      const ri = lu[i];
      const f = ri[k] * inv;
      ri[k] = f;
      if (f !== 0) for (let j = k + 1; j < n; j++) ri[j] -= f * rk[j];
    }
  }
  return { lu, piv, n, sign, singular };
}

/**
 * LU factorisation with partial (row) pivoting, P A = L U (GVL Alg. 3.4.1).
 * Validity: A square and non-singular to working precision (pivot > eps·max|A|); otherwise throws RangeError.
 * @param {number[][]} A n×n matrix
 * @returns {{ lu: Float64Array[], piv: Int32Array, n: number, sign: number }} lu: packed factors (unit-lower L below
 *   the diagonal, U on and above), rows already permuted; piv[k]: row swapped with row k at step k; sign: ±1 parity of
 *   the permutation
 */
export function luFactor(A) {
  const f = luRaw(A, 'luFactor: A');
  if (f.singular) throw new RangeError('luFactor: A is singular to working precision');
  return f;
}

function luSolveInPlace(f, x) {
  const { lu, piv, n } = f;
  for (let k = 0; k < n; k++) {
    const p = piv[k];
    if (p !== k) {
      const t = x[k];
      x[k] = x[p];
      x[p] = t;
    }
  }
  for (let i = 1; i < n; i++) {
    const ri = lu[i];
    let s = x[i];
    for (let j = 0; j < i; j++) s -= ri[j] * x[j];
    x[i] = s;
  }
  for (let i = n - 1; i >= 0; i--) {
    const ri = lu[i];
    let s = x[i];
    for (let j = i + 1; j < n; j++) s -= ri[j] * x[j];
    x[i] = s / ri[i];
  }
  return x;
}

/**
 * Solve A x = b from an LU factorisation.
 * @param {{ lu: Float64Array[], piv: Int32Array, n: number }} f result of luFactor
 * @param {ArrayLike<number>} b right-hand side, length n
 * @returns {number[]} solution x, length n
 */
export function luSolve(f, b) {
  if (!f || !f.lu || !f.piv) throw new RangeError('luSolve: f must be a factorisation returned by luFactor');
  if (!isVec(b) || b.length !== f.n) throw new RangeError('luSolve: b must have length n');
  return luSolveInPlace(f, Array.from(b));
}

/**
 * Solve A x = b by LU with partial pivoting (GVL §3.4).
 * Validity: A square, non-singular (throws RangeError otherwise). O(n³): intended for n up to ~1500.
 * @param {number[][]} A n×n matrix
 * @param {ArrayLike<number> | number[][]} b right-hand side vector (length n) or n×k matrix of k right-hand sides
 * @returns {number[] | number[][]} x (length n) or n×k solution matrix
 */
export function solveDense(A, b) {
  const f = luRaw(A, 'solveDense: A');
  if (f.singular) throw new RangeError('solveDense: A is singular to working precision');
  const n = f.n;
  if (isMat(b)) {
    if (b.length !== n) throw new RangeError('solveDense: b must have n rows');
    const k = b[0].length;
    const X = zeros(n, k);
    const col = new Array(n);
    for (let c = 0; c < k; c++) {
      for (let i = 0; i < n; i++) col[i] = b[i][c];
      luSolveInPlace(f, col);
      for (let i = 0; i < n; i++) X[i][c] = col[i];
    }
    return X;
  }
  if (!isVec(b) || b.length !== n) throw new RangeError('solveDense: b must have length n');
  return luSolveInPlace(f, Array.from(b));
}

/**
 * Cholesky factorisation A = L Lᵀ (GVL Alg. 4.2.1). Only the lower triangle of A is read.
 * Validity: A symmetric positive definite; throws RangeError naming the failing row otherwise.
 * @param {number[][]} A n×n SPD matrix
 * @returns {number[][]} lower-triangular L (upper part zero)
 */
export function cholesky(A) {
  const n = checkSquare(A, 'cholesky: A');
  const L = zeros(n, n);
  for (let i = 0; i < n; i++) {
    const Li = L[i];
    const Ai = A[i];
    for (let j = 0; j <= i; j++) {
      const Lj = L[j];
      let s = Ai[j];
      for (let k = 0; k < j; k++) s -= Li[k] * Lj[k];
      if (i === j) {
        if (!(s > 0) || !Number.isFinite(s)) {
          throw new RangeError(`cholesky: A is not positive definite (pivot ${s} at row ${i})`);
        }
        Li[i] = Math.sqrt(s);
      } else {
        Li[j] = s / Lj[j];
      }
    }
  }
  return L;
}

/**
 * Solve L Lᵀ x = b given the Cholesky factor.
 * @param {number[][]} L lower-triangular factor from cholesky()
 * @param {ArrayLike<number>} b right-hand side, length n
 * @returns {number[]} solution x
 */
export function cholSolve(L, b) {
  const n = checkSquare(L, 'cholSolve: L');
  if (!isVec(b) || b.length !== n) throw new RangeError('cholSolve: b must have length n');
  const x = Array.from(b);
  for (let i = 0; i < n; i++) {
    const Li = L[i];
    let s = x[i];
    for (let k = 0; k < i; k++) s -= Li[k] * x[k];
    x[i] = s / Li[i];
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = x[i];
    for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k];
    x[i] = s / L[i][i];
  }
  return x;
}

/**
 * Matrix inverse by LU with partial pivoting.
 * Validity: A square and non-singular (throws RangeError otherwise). Prefer solveDense when only A⁻¹b is needed.
 * @param {number[][]} A n×n matrix
 * @returns {number[][]} A⁻¹
 */
export function inverse(A) {
  const f = luRaw(A, 'inverse: A');
  if (f.singular) throw new RangeError('inverse: A is singular to working precision');
  const n = f.n;
  const X = zeros(n, n);
  const col = new Array(n);
  for (let c = 0; c < n; c++) {
    col.fill(0);
    col[c] = 1;
    luSolveInPlace(f, col);
    for (let i = 0; i < n; i++) X[i][c] = col[i];
  }
  return X;
}

/**
 * Determinant by LU with partial pivoting. Returns exactly 0 when an exactly zero pivot column is met.
 * Validity: A square; beware overflow/underflow for large n (product of pivots).
 * @param {number[][]} A n×n matrix
 * @returns {number} det(A)
 */
export function det(A) {
  const f = luRaw(A, 'det: A');
  let d = f.sign;
  for (let i = 0; i < f.n; i++) d *= f.lu[i][i];
  return d + 0; // maps −0 to +0
}

/**
 * Linear least squares min ‖A x − b‖₂ by Householder QR (GVL Alg. 5.2.1 and §5.3.2).
 * If A is numerically rank deficient the minimum-norm solution from the pseudo-inverse is returned instead.
 * Validity: A is m×n with m >= n.
 * @param {number[][]} A m×n matrix
 * @param {ArrayLike<number>} b right-hand side, length m
 * @returns {{ x: number[], residual: number, rank: number }} x: solution (length n); residual: ‖A x − b‖₂ (units of
 *   b); rank: numerical rank of A
 */
export function leastSquares(A, b) {
  const [m, n] = checkRect(A, 'leastSquares: A');
  if (m < n) throw new RangeError('leastSquares: A must have at least as many rows as columns');
  if (!isVec(b) || b.length !== m) throw new RangeError('leastSquares: b must have length m');
  const R = new Array(m);
  for (let i = 0; i < m; i++) R[i] = Float64Array.from(A[i]);
  const y = Float64Array.from(b);
  let rmax = 0;
  let rmin = Infinity;
  for (let k = 0; k < n; k++) {
    let nrm = 0;
    let scale = 0;
    for (let i = k; i < m; i++) scale = Math.max(scale, Math.abs(R[i][k]));
    if (scale > 0) {
      for (let i = k; i < m; i++) {
        const v = R[i][k] / scale;
        nrm += v * v;
      }
      nrm = scale * Math.sqrt(nrm);
    }
    if (nrm > 0) {
      const alpha = R[k][k] > 0 ? -nrm : nrm;
      // v = x − alpha e1, stored in place; H = I − v vᵀ / (vᵀv / 2)
      R[k][k] -= alpha;
      let vtv = 0;
      for (let i = k; i < m; i++) vtv += R[i][k] * R[i][k];
      const beta = 2 / vtv;
      for (let j = k + 1; j < n; j++) {
        let s = 0;
        for (let i = k; i < m; i++) s += R[i][k] * R[i][j];
        s *= beta;
        for (let i = k; i < m; i++) R[i][j] -= s * R[i][k];
      }
      let s = 0;
      for (let i = k; i < m; i++) s += R[i][k] * y[i];
      s *= beta;
      for (let i = k; i < m; i++) y[i] -= s * R[i][k];
      R[k][k] = alpha;
    }
    const d = Math.abs(R[k][k]);
    if (d > rmax) rmax = d;
    if (d < rmin) rmin = d;
  }
  if (!(rmin > m * EPS * rmax)) {
    const P = pinv(A);
    const x = matVec(P, b);
    const Ax = matVec(A, x);
    let s = 0;
    for (let i = 0; i < m; i++) s += (Ax[i] - b[i]) * (Ax[i] - b[i]);
    const S = svd(A).S;
    const tol = Math.max(m, n) * EPS * S[0];
    let rank = 0;
    for (let i = 0; i < S.length; i++) if (S[i] > tol) rank++;
    return { x, residual: Math.sqrt(s), rank };
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = y[i];
    for (let j = i + 1; j < n; j++) s -= R[i][j] * x[j];
    x[i] = s / R[i][i];
  }
  let res = 0;
  for (let i = n; i < m; i++) res += y[i] * y[i];
  return { x, residual: Math.sqrt(res), rank: n };
}

/**
 * Thin singular value decomposition A = U diag(S) Vᵀ by one-sided (Hestenes) Jacobi rotations (GVL §8.6.3).
 * Singular values are computed to high relative accuracy. Columns of U belonging to exactly zero singular values
 * are returned as zero vectors.
 * Validity: small/medium dense matrices (cost O(m n²) per sweep, typically 6–10 sweeps).
 * @param {number[][]} A m×n matrix
 * @returns {{ U: number[][], S: number[], V: number[][] }} with k = min(m, n): U is m×k (orthonormal columns),
 *   S has k singular values in descending order, V is n×k (orthonormal columns)
 */
export function svd(A) {
  const [m0, n0] = checkRect(A, 'svd: A');
  const flip = m0 < n0;
  const m = flip ? n0 : m0;
  const n = flip ? m0 : n0;
  // column storage of the (possibly transposed) matrix
  const u = new Array(n);
  const v = new Array(n);
  for (let j = 0; j < n; j++) {
    const c = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      const a = flip ? A[j][i] : A[i][j];
      if (!Number.isFinite(a)) throw new RangeError('svd: A contains a non-finite entry');
      c[i] = a;
    }
    u[j] = c;
    const w = new Float64Array(n);
    w[j] = 1;
    v[j] = w;
  }
  const maxSweeps = 60;
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let rotated = false;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const up = u[p];
        const uq = u[q];
        let alpha = 0;
        let beta = 0;
        let gamma = 0;
        for (let i = 0; i < m; i++) {
          alpha += up[i] * up[i];
          beta += uq[i] * uq[i];
          gamma += up[i] * uq[i];
        }
        if (gamma === 0 || Math.abs(gamma) <= EPS * Math.sqrt(alpha * beta)) continue;
        rotated = true;
        const zeta = (beta - alpha) / (2 * gamma);
        const t = (zeta >= 0 ? 1 : -1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
        const c = 1 / Math.sqrt(1 + t * t);
        const s = c * t;
        for (let i = 0; i < m; i++) {
          const a = up[i];
          const b = uq[i];
          up[i] = c * a - s * b;
          uq[i] = s * a + c * b;
        }
        const vp = v[p];
        const vq = v[q];
        for (let i = 0; i < n; i++) {
          const a = vp[i];
          const b = vq[i];
          vp[i] = c * a - s * b;
          vq[i] = s * a + c * b;
        }
      }
    }
    if (!rotated) break;
  }
  const sv = new Array(n);
  const order = new Array(n);
  for (let j = 0; j < n; j++) {
    sv[j] = norm(u[j]);
    order[j] = j;
  }
  order.sort((a, b) => sv[b] - sv[a]);
  const S = new Array(n);
  const Um = zeros(m, n);
  const Vm = zeros(n, n);
  for (let k = 0; k < n; k++) {
    const j = order[k];
    const s = sv[j];
    S[k] = s;
    const inv = s > 0 ? 1 / s : 0;
    for (let i = 0; i < m; i++) Um[i][k] = u[j][i] * inv;
    for (let i = 0; i < n; i++) Vm[i][k] = v[j][i];
  }
  return flip ? { U: Vm, S, V: Um } : { U: Um, S, V: Vm };
}

/**
 * Moore–Penrose pseudo-inverse from the SVD: A⁺ = V diag(1/S) Uᵀ, dropping singular values <= tol.
 * @param {number[][]} A m×n matrix
 * @param {number} [tol] absolute singular-value cut-off; default max(m, n)·eps·S_max
 * @returns {number[][]} n×m pseudo-inverse
 */
export function pinv(A, tol) {
  const [m, n] = checkRect(A, 'pinv: A');
  const { U, S, V } = svd(A);
  const k = S.length;
  const cut = tol === undefined ? Math.max(m, n) * EPS * (S[0] || 0) : tol;
  if (!(cut >= 0)) throw new RangeError('pinv: tol must be >= 0');
  const P = zeros(n, m);
  for (let r = 0; r < k; r++) {
    if (!(S[r] > cut)) continue;
    const inv = 1 / S[r];
    for (let i = 0; i < n; i++) {
      const f = V[i][r] * inv;
      if (f === 0) continue;
      const Pi = P[i];
      for (let j = 0; j < m; j++) Pi[j] += f * U[j][r];
    }
  }
  return P;
}

/**
 * Tridiagonal solve by the Thomas algorithm (no pivoting; GVL §4.3.6).
 * Row i reads lower[i]·x[i−1] + diag[i]·x[i] + upper[i]·x[i+1] = rhs[i]; lower[0] and upper[n−1] are ignored.
 * `lower` and `upper` may also be given with length n−1 (then lower[i−1] multiplies x[i−1] in row i).
 * Validity: matrix diagonally dominant or SPD (otherwise a zero pivot throws RangeError).
 * @param {ArrayLike<number>} lower sub-diagonal, length n (or n−1)
 * @param {ArrayLike<number>} diag diagonal, length n
 * @param {ArrayLike<number>} upper super-diagonal, length n (or n−1)
 * @param {ArrayLike<number>} rhs right-hand side, length n
 * @returns {Float64Array} solution, length n
 */
export function solveTridiagonal(lower, diag, upper, rhs) {
  const n = diag.length;
  if (!(n >= 1)) throw new RangeError('solveTridiagonal: diag must be non-empty');
  if (rhs.length !== n) throw new RangeError('solveTridiagonal: rhs must have length n');
  if (lower.length !== n && lower.length !== n - 1) throw new RangeError('solveTridiagonal: lower must have length n or n-1');
  if (upper.length !== n && upper.length !== n - 1) throw new RangeError('solveTridiagonal: upper must have length n or n-1');
  const off = lower.length === n ? 0 : -1;
  const c = new Float64Array(n);
  const x = new Float64Array(n);
  let piv = diag[0];
  if (piv === 0 || !Number.isFinite(piv)) throw new RangeError('solveTridiagonal: zero pivot at row 0');
  c[0] = n > 1 ? upper[0] / piv : 0;
  x[0] = rhs[0] / piv;
  for (let i = 1; i < n; i++) {
    const l = lower[i + off];
    piv = diag[i] - l * c[i - 1];
    if (piv === 0 || !Number.isFinite(piv)) throw new RangeError(`solveTridiagonal: zero pivot at row ${i}`);
    if (i < n - 1) c[i] = upper[i] / piv;
    x[i] = (rhs[i] - l * x[i - 1]) / piv;
  }
  for (let i = n - 2; i >= 0; i--) x[i] -= c[i] * x[i + 1];
  return x;
}

/**
 * Banded Cholesky factorisation A = Uᵀ U of a symmetric positive-definite band matrix (GVL Alg. 4.3.5).
 * Storage: bands[i][k] = A[i][i+k] for k = 0..bandwidth (upper band, row-wise). Rows may be shorter near the end of
 * the matrix or wherever trailing entries are zero. Cost O(n·bw²).
 * Validity: A SPD; throws RangeError with the failing row otherwise.
 * @param {ArrayLike<number>[]} bands n rows of upper-band coefficients
 * @returns {{ n: number, bw: number, u: Float64Array }} factor: u[i·(bw+1)+k] = U[i][i+k]
 */
export function factorBandedSym(bands) {
  if (!Array.isArray(bands) || bands.length === 0) throw new RangeError('factorBandedSym: bands must be a non-empty array of rows');
  const n = bands.length;
  let bw = 0;
  for (let i = 0; i < n; i++) {
    if (!isVec(bands[i]) || bands[i].length < 1) throw new RangeError(`factorBandedSym: bands[${i}] must hold at least the diagonal`);
    bw = Math.max(bw, Math.min(bands[i].length - 1, n - 1 - i));
  }
  const w = bw + 1;
  const u = new Float64Array(n * w);
  for (let i = 0; i < n; i++) {
    const row = bands[i];
    const kmax = Math.min(row.length - 1, n - 1 - i);
    for (let k = 0; k <= kmax; k++) u[i * w + k] = row[k];
  }
  for (let i = 0; i < n; i++) {
    const ri = i * w;
    const d = u[ri];
    if (!(d > 0) || !Number.isFinite(d)) {
      throw new RangeError(`factorBandedSym: matrix is not positive definite (pivot ${d} at row ${i})`);
    }
    const s = Math.sqrt(d);
    const inv = 1 / s;
    u[ri] = s;
    const kmax = Math.min(bw, n - 1 - i);
    for (let k = 1; k <= kmax; k++) u[ri + k] *= inv;
    for (let k = 1; k <= kmax; k++) {
      const f = u[ri + k];
      if (f === 0) continue;
      const rj = (i + k) * w - k;
      for (let l = k; l <= kmax; l++) u[rj + l] -= f * u[ri + l];
    }
  }
  return { n, bw, u };
}

/**
 * Solve A x = b from a banded Cholesky factor.
 * @param {{ n: number, bw: number, u: Float64Array }} factor result of factorBandedSym
 * @param {ArrayLike<number>} b right-hand side, length n
 * @returns {Float64Array} solution x
 */
export function solveBandedSymFactored(factor, b) {
  if (!factor || !factor.u) throw new RangeError('solveBandedSymFactored: factor must come from factorBandedSym');
  const { n, bw, u } = factor;
  if (!isVec(b) || b.length !== n) throw new RangeError('solveBandedSymFactored: b must have length n');
  const w = bw + 1;
  const x = Float64Array.from(b);
  for (let i = 0; i < n; i++) {
    const ri = i * w;
    const xi = x[i] / u[ri];
    x[i] = xi;
    if (xi === 0) continue;
    const kmax = Math.min(bw, n - 1 - i);
    for (let k = 1; k <= kmax; k++) x[i + k] -= u[ri + k] * xi;
  }
  for (let i = n - 1; i >= 0; i--) {
    const ri = i * w;
    let s = x[i];
    const kmax = Math.min(bw, n - 1 - i);
    for (let k = 1; k <= kmax; k++) s -= u[ri + k] * x[i + k];
    x[i] = s / u[ri];
  }
  return x;
}

/**
 * Solve a symmetric positive-definite banded system (factor + solve).
 * @param {ArrayLike<number>[]} bands bands[i][k] = A[i][i+k], k = 0..bandwidth
 * @param {ArrayLike<number>} b right-hand side, length n
 * @returns {Float64Array} solution x
 */
export function solveBandedSym(bands, b) {
  return solveBandedSymFactored(factorBandedSym(bands), b);
}
