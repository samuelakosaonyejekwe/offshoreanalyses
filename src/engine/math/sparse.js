/**
 * Sparse linear algebra for finite-element sized systems: COO assembly to CSR, preconditioned conjugate
 * gradients, a symmetric skyline (profile) LDLᵀ direct solver and reverse Cuthill–McKee reordering.
 *
 * Everything is stored in typed arrays. Unit-agnostic (K in N/m and b in N give x in m).
 *
 * References: Saad, "Iterative Methods for Sparse Linear Systems", 2nd ed. (2003): CSR §3.4, PCG Alg. 9.1,
 * SSOR §10.2, IC(0) §10.3.2; Manteuffel, Math. Comp. 34 (1980) 473–497 (shifted incomplete Cholesky);
 * Bathe, "Finite Element Procedures" (1996) §8.2.3 and App. (active-column skyline LDLᵀ, COLSOL);
 * Cuthill & McKee, Proc. 24th ACM Nat. Conf. (1969) 157–172; George & Liu, "Computer Solution of Large Sparse
 * Positive Definite Systems" (1981) §4.3 (reverse ordering, pseudo-peripheral start node).
 */

function isVec(v) {
  return Array.isArray(v) || ArrayBuffer.isView(v);
}

function checkOrder(n, name) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError(`${name}: n must be an integer >= 1`);
}

/**
 * Coordinate-format builder for a square sparse matrix. Entries may be added in any order; duplicates are summed
 * when converting to CSR (the usual finite-element assembly).
 */
export class CooBuilder {
  /**
   * @param {number} n matrix order (number of rows = number of columns)
   */
  constructor(n) {
    checkOrder(n, 'CooBuilder');
    this.n = n;
    this.count = 0;
    this._i = new Int32Array(1024);
    this._j = new Int32Array(1024);
    this._v = new Float64Array(1024);
  }

  /**
   * Add v to entry (i, j). Validity: 0 <= i, j < n (throws RangeError otherwise).
   * @param {number} i row index
   * @param {number} j column index
   * @param {number} v value to accumulate
   */
  add(i, j, v) {
    const n = this.n;
    if (!(i >= 0 && i < n && j >= 0 && j < n)) throw new RangeError(`CooBuilder.add: index (${i}, ${j}) outside 0..${n - 1}`);
    let c = this.count;
    if (c === this._v.length) {
      const cap = c * 2;
      const ni = new Int32Array(cap);
      const nj = new Int32Array(cap);
      const nv = new Float64Array(cap);
      ni.set(this._i);
      nj.set(this._j);
      nv.set(this._v);
      this._i = ni;
      this._j = nj;
      this._v = nv;
    }
    this._i[c] = i;
    this._j[c] = j;
    this._v[c] = v;
    this.count = c + 1;
  }

  /**
   * Convert to compressed sparse row storage with duplicates summed and column indices ascending in every row
   * (two stable counting sorts, O(nnz + n)). Explicit zeros are kept so that the sparsity pattern is preserved.
   * @returns {{ n: number, rowPtr: Int32Array, colIdx: Int32Array, val: Float64Array }} CSR matrix
   */
  toCsr() {
    const n = this.n;
    const nnz = this.count;
    const I = this._i;
    const J = this._j;
    const V = this._v;
    // pass 1: stable sort by column
    const cnt = new Int32Array(n + 1);
    for (let k = 0; k < nnz; k++) cnt[J[k] + 1]++;
    for (let c = 0; c < n; c++) cnt[c + 1] += cnt[c];
    const ord1 = new Int32Array(nnz);
    for (let k = 0; k < nnz; k++) ord1[cnt[J[k]]++] = k;
    // pass 2: stable sort by row
    const rowStart = new Int32Array(n + 1);
    for (let k = 0; k < nnz; k++) rowStart[I[k] + 1]++;
    for (let r = 0; r < n; r++) rowStart[r + 1] += rowStart[r];
    const next = rowStart.slice(0, n);
    const ord2 = new Int32Array(nnz);
    for (let q = 0; q < nnz; q++) {
      const k = ord1[q];
      ord2[next[I[k]]++] = k;
    }
    // merge duplicates
    const rowPtr = new Int32Array(n + 1);
    const colTmp = new Int32Array(nnz);
    const valTmp = new Float64Array(nnz);
    let w = 0;
    for (let r = 0; r < n; r++) {
      rowPtr[r] = w;
      let lastCol = -1;
      for (let q = rowStart[r]; q < rowStart[r + 1]; q++) {
        const k = ord2[q];
        const c = J[k];
        if (c === lastCol) {
          valTmp[w - 1] += V[k];
        } else {
          colTmp[w] = c;
          valTmp[w] = V[k];
          lastCol = c;
          w++;
        }
      }
    }
    rowPtr[n] = w;
    return { n, rowPtr, colIdx: colTmp.slice(0, w), val: valTmp.slice(0, w) };
  }
}

function checkCsr(A, name) {
  if (!A || !Number.isInteger(A.n) || !A.rowPtr || !A.colIdx || !A.val) {
    throw new RangeError(`${name}: A must be a CSR matrix { n, rowPtr, colIdx, val }`);
  }
}

/**
 * CSR matrix–vector product y = A x.
 * @param {{ n: number, rowPtr: Int32Array, colIdx: Int32Array, val: Float64Array }} A CSR matrix
 * @param {ArrayLike<number>} x vector, length n
 * @param {Float64Array} [out] optional result buffer (length n, must not alias x)
 * @returns {Float64Array} y = A x
 */
export function csrMatVec(A, x, out) {
  checkCsr(A, 'csrMatVec');
  const n = A.n;
  if (!isVec(x) || x.length !== n) throw new RangeError('csrMatVec: x must have length n');
  const y = out || new Float64Array(n);
  const { rowPtr, colIdx, val } = A;
  for (let i = 0; i < n; i++) {
    let s = 0;
    const e = rowPtr[i + 1];
    for (let k = rowPtr[i]; k < e; k++) s += val[k] * x[colIdx[k]];
    y[i] = s;
  }
  return y;
}

function csrDiagonal(A, name) {
  const { n, rowPtr, colIdx, val } = A;
  const d = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) if (colIdx[k] === i) s += val[k];
    if (!(s > 0)) throw new RangeError(`${name}: diagonal entry of row ${i} is not positive (matrix is not SPD)`);
    d[i] = s;
  }
  return d;
}

/** Returns A itself when every row has strictly ascending columns, otherwise a sorted, duplicate-free copy. */
function canonicalCsr(A) {
  const { n, rowPtr, colIdx, val } = A;
  let ok = true;
  for (let i = 0; i < n && ok; i++) {
    for (let k = rowPtr[i] + 1; k < rowPtr[i + 1]; k++) {
      if (colIdx[k] <= colIdx[k - 1]) {
        ok = false;
        break;
      }
    }
  }
  if (ok) return A;
  const coo = new CooBuilder(n);
  for (let i = 0; i < n; i++) for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) coo.add(i, colIdx[k], val[k]);
  return coo.toCsr();
}

/**
 * Zero-fill incomplete Cholesky A ≈ L Lᵀ on the pattern of the lower triangle (Saad §10.3.2), with the diagonal
 * scaled by (1 + shift). Rows of L are stored with ascending columns, diagonal last. Returns null on breakdown.
 */
function ic0Factor(A0, shift) {
  const A = canonicalCsr(A0);
  const { n, rowPtr, colIdx, val } = A;
  const lptr = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    let c = 0;
    let hasDiag = false;
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
      if (colIdx[k] < i) c++;
      else if (colIdx[k] === i) hasDiag = true;
    }
    if (!hasDiag) return null;
    lptr[i + 1] = lptr[i] + c + 1;
  }
  const nl = lptr[n];
  const lcol = new Int32Array(nl);
  const lval = new Float64Array(nl);
  for (let i = 0; i < n; i++) {
    let w = lptr[i];
    for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
      const c = colIdx[k];
      if (c > i) break;
      lcol[w] = c;
      lval[w] = c === i ? val[k] * (1 + shift) : val[k];
      w++;
    }
  }
  for (let i = 0; i < n; i++) {
    const rs = lptr[i];
    const dpos = lptr[i + 1] - 1;
    for (let kk = rs; kk <= dpos; kk++) {
      const j = lcol[kk];
      let s = lval[kk];
      let p = rs;
      let q = lptr[j];
      const qEnd = j === i ? kk : lptr[j + 1] - 1;
      while (p < kk && q < qEnd) {
        const cp = lcol[p];
        const cq = lcol[q];
        if (cp === cq) {
          s -= lval[p] * lval[q];
          p++;
          q++;
        } else if (cp < cq) {
          p++;
        } else {
          q++;
        }
      }
      if (j < i) {
        lval[kk] = s / lval[lptr[j + 1] - 1];
      } else {
        if (!(s > 0) || !Number.isFinite(s)) return null;
        lval[kk] = Math.sqrt(s);
      }
    }
  }
  return { lptr, lcol, lval };
}

function makePreconditioner(A, precond, omega) {
  const { n, rowPtr, colIdx, val } = A;
  if (typeof precond === 'function') return precond;
  if (precond === 'none') {
    return (r, z) => {
      z.set(r);
    };
  }
  const d = csrDiagonal(A, 'pcg');
  if (precond === 'jacobi') {
    const inv = new Float64Array(n);
    for (let i = 0; i < n; i++) inv[i] = 1 / d[i];
    return (r, z) => {
      for (let i = 0; i < n; i++) z[i] = r[i] * inv[i];
    };
  }
  if (precond === 'ssor') {
    if (!(omega > 0 && omega < 2)) throw new RangeError('pcg: omega must be in (0, 2)');
    const sc = (2 - omega) / omega;
    return (r, z) => {
      // (D/ω + L) y = r
      for (let i = 0; i < n; i++) {
        let s = r[i];
        for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
          const c = colIdx[k];
          if (c < i) s -= val[k] * z[c];
        }
        z[i] = (s * omega) / d[i];
      }
      // y ← ((2 − ω)/ω) D y, then (D/ω + U) z = y
      for (let i = n - 1; i >= 0; i--) {
        let s = z[i] * d[i] * sc;
        for (let k = rowPtr[i]; k < rowPtr[i + 1]; k++) {
          const c = colIdx[k];
          if (c > i) s -= val[k] * z[c];
        }
        z[i] = (s * omega) / d[i];
      }
    };
  }
  if (precond === 'ic0') {
    // Manteuffel shift: retry with A + α·diag(A) until the incomplete factorisation exists
    let fac = null;
    let shift = 0;
    for (let attempt = 0; attempt < 8 && !fac; attempt++) {
      fac = ic0Factor(A, shift);
      shift = shift === 0 ? 1e-3 : shift * 10;
    }
    if (!fac) throw new RangeError('pcg: incomplete Cholesky factorisation failed (matrix is not SPD)');
    const { lptr, lcol, lval } = fac;
    const y = new Float64Array(n);
    return (r, z) => {
      for (let i = 0; i < n; i++) {
        let s = r[i];
        const dpos = lptr[i + 1] - 1;
        for (let k = lptr[i]; k < dpos; k++) s -= lval[k] * y[lcol[k]];
        y[i] = s / lval[dpos];
      }
      for (let i = n - 1; i >= 0; i--) {
        const dpos = lptr[i + 1] - 1;
        const zi = y[i] / lval[dpos];
        z[i] = zi;
        for (let k = lptr[i]; k < dpos; k++) y[lcol[k]] -= lval[k] * zi;
      }
    };
  }
  throw new RangeError("pcg: precond must be 'jacobi', 'ssor', 'ic0', 'none' or a function (r, z)");
}

/**
 * Preconditioned conjugate gradients for a symmetric positive-definite CSR system A x = b (Saad Alg. 9.1).
 * Preconditioners: 'jacobi' (diagonal scaling), 'ssor' (symmetric successive over-relaxation, ω = 1 by default),
 * 'ic0' (zero-fill incomplete Cholesky with Manteuffel diagonal shift on breakdown), 'none', or a user function
 * (r, z) writing z = M⁻¹ r.
 * Validity: A symmetric positive definite with the full matrix stored (both triangles). Does not throw on
 * non-convergence: check `converged`.
 * @param {{ n: number, rowPtr: Int32Array, colIdx: Int32Array, val: Float64Array }} A CSR matrix
 * @param {ArrayLike<number>} b right-hand side, length n
 * @param {{ x0?: ArrayLike<number>, tol?: number, maxIter?: number,
 *   precond?: 'jacobi' | 'ssor' | 'ic0' | 'none' | ((r: Float64Array, z: Float64Array) => void), omega?: number }} [opts]
 *   x0: starting guess (default 0); tol: relative residual tolerance ‖b − A x‖₂ <= tol·‖b‖₂ (default 1e-10);
 *   maxIter: iteration cap (default max(200, 4n)); precond (default 'jacobi'); omega: SSOR relaxation in (0, 2)
 * @returns {{ x: Float64Array, iterations: number, residual: number, converged: boolean }} x: solution;
 *   residual: final relative residual ‖b − A x‖₂ / ‖b‖₂ (recurrence value)
 */
export function pcg(A, b, { x0, tol = 1e-10, maxIter, precond = 'jacobi', omega = 1 } = {}) {
  checkCsr(A, 'pcg');
  const n = A.n;
  if (!isVec(b) || b.length !== n) throw new RangeError('pcg: b must have length n');
  if (x0 !== undefined && x0 !== null && (!isVec(x0) || x0.length !== n)) throw new RangeError('pcg: x0 must have length n');
  if (!(tol > 0)) throw new RangeError('pcg: tol must be > 0');
  const itMax = maxIter !== undefined ? maxIter : Math.max(200, 4 * n);
  const x = x0 ? Float64Array.from(x0) : new Float64Array(n);
  const r = new Float64Array(n);
  let bnorm = 0;
  for (let i = 0; i < n; i++) bnorm += b[i] * b[i];
  bnorm = Math.sqrt(bnorm);
  if (!Number.isFinite(bnorm)) throw new RangeError('pcg: b contains a non-finite entry');
  if (bnorm === 0) return { x: new Float64Array(n), iterations: 0, residual: 0, converged: true };
  const apply = makePreconditioner(A, precond, omega);
  const z = new Float64Array(n);
  const p = new Float64Array(n);
  const Ap = new Float64Array(n);
  if (x0) {
    csrMatVec(A, x, Ap);
    for (let i = 0; i < n; i++) r[i] = b[i] - Ap[i];
  } else {
    for (let i = 0; i < n; i++) r[i] = b[i];
  }
  let rnorm = 0;
  for (let i = 0; i < n; i++) rnorm += r[i] * r[i];
  rnorm = Math.sqrt(rnorm);
  let iterations = 0;
  let converged = rnorm <= tol * bnorm;
  if (!converged) {
    apply(r, z);
    p.set(z);
    let rz = 0;
    for (let i = 0; i < n; i++) rz += r[i] * z[i];
    while (iterations < itMax) {
      csrMatVec(A, p, Ap);
      let pAp = 0;
      for (let i = 0; i < n; i++) pAp += p[i] * Ap[i];
      if (!(pAp > 0) || !Number.isFinite(pAp)) break; // not positive definite (or breakdown)
      const alpha = rz / pAp;
      let rr = 0;
      for (let i = 0; i < n; i++) {
        x[i] += alpha * p[i];
        const ri = r[i] - alpha * Ap[i];
        r[i] = ri;
        rr += ri * ri;
      }
      iterations++;
      rnorm = Math.sqrt(rr);
      if (rnorm <= tol * bnorm) {
        converged = true;
        break;
      }
      apply(r, z);
      let rzNew = 0;
      for (let i = 0; i < n; i++) rzNew += r[i] * z[i];
      const beta = rzNew / rz;
      rz = rzNew;
      for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
    }
  }
  return { x, iterations, residual: rnorm / bnorm, converged };
}

/**
 * Symmetric skyline (variable-band, active-column) matrix with in-place LDLᵀ factorisation for finite-element
 * systems (Bathe 1996 §8.2.3). Only the upper triangle inside the column profile is stored, column by column, in
 * one Float64Array.
 *
 * Profile: colHeights[j] = j − (row index of the first non-zero in column j), i.e. the number of stored entries
 * above the diagonal of column j (0 for a diagonal-only column); values larger than j are clamped.
 *
 * Assembly convention: `add(i, j, v)` accumulates into the upper triangle only. Calls with i > j are ignored, so
 * looping over a full symmetric element matrix assembles each coefficient exactly once. Do not assemble only the
 * lower triangle.
 *
 * Typical use: K = new SkylineSym(n, heights); K.add(...); K.applyBc(dof); K.factor(); x = K.solve(b).
 * Cost: factor O(n·h²/2) for mean column height h, solve O(2·n·h).
 */
export class SkylineSym {
  /**
   * @param {number} n matrix order
   * @param {ArrayLike<number>} colHeights column heights above the diagonal, length n (see class description)
   */
  constructor(n, colHeights) {
    checkOrder(n, 'SkylineSym');
    if (!isVec(colHeights) || colHeights.length !== n) throw new RangeError('SkylineSym: colHeights must have length n');
    this.n = n;
    /** first stored row of each column */
    this.first = new Int32Array(n);
    /** column j occupies val[ptr[j] .. ptr[j+1]−1], rows first[j]..j (diagonal last) */
    this.ptr = new Int32Array(n + 1);
    let total = 0;
    for (let j = 0; j < n; j++) {
      let h = colHeights[j];
      if (!(h >= 0) || !Number.isInteger(h)) throw new RangeError(`SkylineSym: colHeights[${j}] must be a non-negative integer`);
      if (h > j) h = j;
      this.first[j] = j - h;
      this.ptr[j] = total;
      total += h + 1;
      if (total > 2147483647) throw new RangeError('SkylineSym: profile too large');
    }
    this.ptr[n] = total;
    this.val = new Float64Array(total);
    this.factored = false;
    /** number of negative pivots found by factor() (0 for a positive-definite matrix; Sturm count otherwise) */
    this.negativePivots = 0;
  }

  /**
   * Accumulate v into K[i][j] (upper triangle). Calls with i > j are ignored (see class description).
   * Validity: indices in range and (i, j) inside the profile (throws RangeError otherwise); not after factor().
   * @param {number} i row index
   * @param {number} j column index
   * @param {number} v value to add
   */
  add(i, j, v) {
    if (i > j) {
      if (!(i < this.n && j >= 0)) throw new RangeError(`SkylineSym.add: index (${i}, ${j}) outside 0..${this.n - 1}`);
      return;
    }
    if (!(i >= 0 && j < this.n)) throw new RangeError(`SkylineSym.add: index (${i}, ${j}) outside 0..${this.n - 1}`);
    if (this.factored) throw new RangeError('SkylineSym.add: matrix is already factorised');
    const f = this.first[j];
    if (i < f) throw new RangeError(`SkylineSym.add: entry (${i}, ${j}) lies outside the column profile (first row ${f})`);
    this.val[this.ptr[j] + i - f] += v;
  }

  /**
   * Read K[i][j] (symmetric access; zero outside the profile). After factor() this returns factor entries.
   * @param {number} i row index
   * @param {number} j column index
   * @returns {number} stored coefficient
   */
  get(i, j) {
    if (i > j) {
      const t = i;
      i = j;
      j = t;
    }
    if (!(i >= 0 && j < this.n)) throw new RangeError(`SkylineSym.get: index (${i}, ${j}) outside 0..${this.n - 1}`);
    const f = this.first[j];
    return i < f ? 0 : this.val[this.ptr[j] + i - f];
  }

  /**
   * Impose a prescribed value on a degree of freedom: row and column `dof` are zeroed and the diagonal set to 1,
   * so that the solution returns x[dof] = rhs[dof]. When `rhs` is given it is updated in place for the eliminated
   * column (rhs[i] −= K[i][dof]·value) and rhs[dof] is set to `value`.
   * Validity: before factor() only.
   * @param {number} dof degree-of-freedom index
   * @param {number} [value] prescribed value (default 0)
   * @param {Float64Array | number[]} [rhs] right-hand side to modify in place
   */
  applyBc(dof, value = 0, rhs) {
    const n = this.n;
    if (!Number.isInteger(dof) || dof < 0 || dof >= n) throw new RangeError(`SkylineSym.applyBc: dof ${dof} outside 0..${n - 1}`);
    if (this.factored) throw new RangeError('SkylineSym.applyBc: matrix is already factorised');
    if (rhs !== undefined && rhs !== null && (!isVec(rhs) || rhs.length !== n)) throw new RangeError('SkylineSym.applyBc: rhs must have length n');
    const { first, ptr, val } = this;
    const useRhs = !!rhs && value !== 0;
    const f = first[dof];
    const base = ptr[dof] - f;
    for (let i = f; i < dof; i++) {
      if (useRhs) rhs[i] -= val[base + i] * value;
      val[base + i] = 0;
    }
    for (let j = dof + 1; j < n; j++) {
      if (first[j] <= dof) {
        const q = ptr[j] + dof - first[j];
        if (useRhs) rhs[j] -= val[q] * value;
        val[q] = 0;
      }
    }
    val[ptr[dof + 1] - 1] = 1;
    if (rhs) rhs[dof] = value;
  }

  /**
   * y = K x with the assembled (unfactorised) matrix.
   * @param {ArrayLike<number>} x vector, length n
   * @param {Float64Array} [out] optional result buffer (must not alias x)
   * @returns {Float64Array} K x
   */
  mulVec(x, out) {
    const n = this.n;
    if (this.factored) throw new RangeError('SkylineSym.mulVec: matrix is already factorised');
    if (!isVec(x) || x.length !== n) throw new RangeError('SkylineSym.mulVec: x must have length n');
    const y = out || new Float64Array(n);
    y.fill(0);
    const { first, ptr, val } = this;
    for (let j = 0; j < n; j++) {
      const f = first[j];
      const base = ptr[j] - f;
      const xj = x[j];
      let s = 0;
      for (let i = f; i < j; i++) {
        const a = val[base + i];
        s += a * x[i];
        y[i] += a * xj;
      }
      y[j] += s + val[base + j] * xj;
    }
    return y;
  }

  /**
   * In-place LDLᵀ factorisation (active-column reduction, no pivoting). Negative pivots are allowed (indefinite
   * matrices factor as long as no pivot vanishes) and counted in `negativePivots`.
   * Validity: throws Error naming the degree of freedom when a pivot is zero or negligible relative to the
   * original diagonal (|d| <= 1e-12·|K_jj|): unconstrained rigid-body mode or mechanism. The error carries `.dof`.
   * @returns {SkylineSym} this
   */
  factor() {
    if (this.factored) return this;
    const { n, first, ptr, val: a } = this;
    let neg = 0;
    for (let j = 0; j < n; j++) {
      const mj = first[j];
      const pj = ptr[j] - mj;
      const kjj = a[pj + j];
      for (let i = mj + 1; i < j; i++) {
        const mi = first[i];
        const lo = mi > mj ? mi : mj;
        const pi = ptr[i] - mi;
        let s = 0;
        let r = lo;
        const r4 = i - 3;
        for (; r < r4; r += 4) {
          s += a[pi + r] * a[pj + r] + a[pi + r + 1] * a[pj + r + 1] + a[pi + r + 2] * a[pj + r + 2] + a[pi + r + 3] * a[pj + r + 3];
        }
        for (; r < i; r++) s += a[pi + r] * a[pj + r];
        a[pj + i] -= s;
      }
      let d = kjj;
      for (let i = mj; i < j; i++) {
        const g = a[pj + i];
        const l = g / a[ptr[i + 1] - 1];
        a[pj + i] = l;
        d -= l * g;
      }
      if (!Number.isFinite(d) || d === 0 || Math.abs(d) <= 1e-12 * Math.abs(kjj)) {
        const err = new Error(`SkylineSym.factor: zero pivot at dof ${j} (matrix singular: unconstrained or mechanism)`);
        err.dof = j;
        throw err;
      }
      if (d < 0) neg++;
      a[pj + j] = d;
    }
    this.negativePivots = neg;
    this.factored = true;
    return this;
  }

  /**
   * Solve K x = b with the factorised matrix (factor() is called automatically if needed). b is not modified.
   * @param {ArrayLike<number>} b right-hand side, length n
   * @returns {Float64Array} solution x
   */
  solve(b) {
    const n = this.n;
    if (!isVec(b) || b.length !== n) throw new RangeError('SkylineSym.solve: b must have length n');
    if (!this.factored) this.factor();
    const { first, ptr, val: a } = this;
    const x = Float64Array.from(b);
    for (let j = 1; j < n; j++) {
      const mj = first[j];
      const pj = ptr[j] - mj;
      let s = 0;
      for (let r = mj; r < j; r++) s += a[pj + r] * x[r];
      x[j] -= s;
    }
    for (let j = 0; j < n; j++) x[j] /= a[ptr[j + 1] - 1];
    for (let j = n - 1; j > 0; j--) {
      const xj = x[j];
      if (xj === 0) continue;
      const mj = first[j];
      const pj = ptr[j] - mj;
      for (let r = mj; r < j; r++) x[r] -= a[pj + r] * xj;
    }
    return x;
  }
}

/** Adjacency structure (CSR, neighbours unique and ascending) of an undirected graph given as an edge list. */
function buildAdjacency(n, edges, name) {
  const m = edges.length;
  const deg = new Int32Array(n + 1);
  for (let e = 0; e < m; e++) {
    const i = edges[e][0];
    const j = edges[e][1];
    if (!(i >= 0 && i < n && j >= 0 && j < n) || !Number.isInteger(i) || !Number.isInteger(j)) {
      throw new RangeError(`${name}: edge ${e} = [${i}, ${j}] has a node outside 0..${n - 1}`);
    }
    if (i === j) continue;
    deg[i + 1]++;
    deg[j + 1]++;
  }
  for (let i = 0; i < n; i++) deg[i + 1] += deg[i];
  const fill = deg.slice(0, n);
  const adjRaw = new Int32Array(deg[n]);
  for (let e = 0; e < m; e++) {
    const i = edges[e][0];
    const j = edges[e][1];
    if (i === j) continue;
    adjRaw[fill[i]++] = j;
    adjRaw[fill[j]++] = i;
  }
  const xadj = new Int32Array(n + 1);
  const adj = new Int32Array(deg[n]);
  let w = 0;
  for (let i = 0; i < n; i++) {
    xadj[i] = w;
    const seg = adjRaw.subarray(deg[i], deg[i + 1]);
    seg.sort();
    let last = -1;
    for (let k = 0; k < seg.length; k++) {
      if (seg[k] !== last) {
        last = seg[k];
        adj[w++] = last;
      }
    }
  }
  xadj[n] = w;
  return { xadj, adj };
}

/**
 * Reverse Cuthill–McKee ordering of an undirected graph, to reduce the bandwidth/profile of a symmetric sparse
 * matrix before a banded or skyline factorisation. Each connected component is started from a pseudo-peripheral
 * node (George & Liu 1981 §4.3) and traversed breadth-first with neighbours taken in order of increasing degree;
 * the final order is reversed.
 * Validity: node indices 0..n−1; self-loops and duplicate edges are ignored; disconnected graphs are handled.
 * @param {number} n number of nodes (matrix order)
 * @param {ArrayLike<ArrayLike<number>>} edges list of [i, j] couplings (non-zero off-diagonal positions)
 * @returns {Int32Array} permutation new -> old: perm[k] is the original index of the node numbered k
 */
export function reverseCuthillMcKee(n, edges) {
  checkOrder(n, 'reverseCuthillMcKee');
  if (!edges || typeof edges.length !== 'number') throw new RangeError('reverseCuthillMcKee: edges must be an array of [i, j] pairs');
  const { xadj, adj } = buildAdjacency(n, edges, 'reverseCuthillMcKee');
  const degree = new Int32Array(n);
  for (let i = 0; i < n; i++) degree[i] = xadj[i + 1] - xadj[i];
  const order = new Int32Array(n);
  const visited = new Uint8Array(n);
  const level = new Int32Array(n);
  const mark = new Int32Array(n).fill(-1);
  const queue = new Int32Array(n);
  let stamp = 0;

  // breadth-first search from s inside the unvisited part; returns [eccentricity, min-degree node of last level]
  const bfsLevels = (s) => {
    stamp++;
    let head = 0;
    let tail = 0;
    queue[tail++] = s;
    mark[s] = stamp;
    level[s] = 0;
    let ecc = 0;
    while (head < tail) {
      const u = queue[head++];
      ecc = level[u];
      for (let k = xadj[u]; k < xadj[u + 1]; k++) {
        const v = adj[k];
        if (mark[v] !== stamp && !visited[v]) {
          mark[v] = stamp;
          level[v] = ecc + 1;
          queue[tail++] = v;
        }
      }
    }
    let best = queue[tail - 1];
    for (let q = tail - 1; q >= 0 && level[queue[q]] === ecc; q--) {
      if (degree[queue[q]] < degree[best]) best = queue[q];
    }
    return [ecc, best];
  };

  let count = 0;
  // nodes by increasing degree give the candidate start of every component
  const byDegree = Array.from({ length: n }, (_, i) => i).sort((p, q) => degree[p] - degree[q] || p - q);
  for (let c = 0; c < n; c++) {
    let start = byDegree[c];
    if (visited[start]) continue;
    let [ecc, cand] = bfsLevels(start);
    for (let it = 0; it < 8; it++) {
      const [ecc2, cand2] = bfsLevels(cand);
      if (ecc2 > ecc) {
        start = cand;
        ecc = ecc2;
        cand = cand2;
      } else {
        break;
      }
    }
    // Cuthill–McKee traversal of this component
    let head = count;
    order[count++] = start;
    visited[start] = 1;
    while (head < count) {
      const u = order[head++];
      const from = count;
      for (let k = xadj[u]; k < xadj[u + 1]; k++) {
        const v = adj[k];
        if (!visited[v]) {
          visited[v] = 1;
          order[count++] = v;
        }
      }
      // insertion sort of the newly added neighbours by degree
      for (let p = from + 1; p < count; p++) {
        const v = order[p];
        const dv = degree[v];
        let q = p - 1;
        while (q >= from && degree[order[q]] > dv) {
          order[q + 1] = order[q];
          q--;
        }
        order[q + 1] = v;
      }
    }
  }
  order.reverse();
  return order;
}

/**
 * Column heights of the skyline profile for a given coupling graph, optionally after renumbering.
 * heights[j] = j − min(i) over all couplings (i, j) with i < j, in the (new) numbering: the argument expected by
 * SkylineSym.
 * @param {number} n matrix order
 * @param {ArrayLike<ArrayLike<number>>} edges list of [i, j] couplings in the original numbering
 * @param {ArrayLike<number>} [perm] permutation new -> old (e.g. from reverseCuthillMcKee); identity if omitted
 * @returns {Int32Array} column heights, length n
 */
export function skylineHeights(n, edges, perm) {
  checkOrder(n, 'skylineHeights');
  let inv = null;
  if (perm !== undefined && perm !== null) {
    if (!isVec(perm) || perm.length !== n) throw new RangeError('skylineHeights: perm must have length n');
    inv = new Int32Array(n);
    for (let k = 0; k < n; k++) inv[perm[k]] = k;
  }
  const h = new Int32Array(n);
  for (let e = 0; e < edges.length; e++) {
    let i = edges[e][0];
    let j = edges[e][1];
    if (!(i >= 0 && i < n && j >= 0 && j < n)) throw new RangeError(`skylineHeights: edge ${e} has a node outside 0..${n - 1}`);
    if (inv) {
      i = inv[i];
      j = inv[j];
    }
    if (i > j) {
      const t = i;
      i = j;
      j = t;
    }
    if (j - i > h[j]) h[j] = j - i;
  }
  return h;
}
