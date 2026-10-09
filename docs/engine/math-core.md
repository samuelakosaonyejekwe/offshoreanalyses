# Numerical core — `src/engine/math`

One line per export: signature, shapes/units, algorithm. Dense matrices are arrays of rows; vectors may be plain
arrays or `Float64Array`. All routines are unit-agnostic unless a unit is stated; inputs are never mutated unless
stated. Invalid input throws `RangeError`; iteration failures throw `Error` or return `converged: false` as noted.
Tests: `tests/engine/math.<file>.test.js`.

## linalg.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `zeros(n, m?) -> number[] \| number[][]` | vector of length n, or n×m matrix | fill |
| `identity(n) -> number[][]` | n×n | — |
| `matMul(A, B) -> number[][]` | (m×p)(p×n) -> m×n | row-oriented triple loop |
| `matVec(A, x) -> number[]` | (m×n)(n) -> m | — |
| `transpose(A) -> number[][]` | m×n -> n×m | — |
| `dot(a, b) -> number` | equal lengths | — |
| `norm(a) -> number` | Euclidean norm | scaled sum of squares (no overflow/underflow) |
| `axpy(alpha, x, y) -> number[]` | alpha·x + y (new vector) | — |
| `solveDense(A, b) -> number[] \| number[][]` | A n×n; b length n, or n×k (multi-RHS -> n×k) | LU with partial pivoting; throws if singular |
| `luFactor(A) -> { lu, piv, n, sign }` | lu: packed rows (`Float64Array[]`), piv: row swaps | LU with partial pivoting; throws if singular |
| `luSolve(f, b) -> number[]` | f from `luFactor`, b length n | forward/back substitution |
| `cholesky(A) -> number[][]` | lower L with A = L Lᵀ (reads lower triangle) | Cholesky–Banachiewicz; throws if not SPD |
| `cholSolve(L, b) -> number[]` | L from `cholesky` | two triangular solves |
| `inverse(A) -> number[][]` | n×n | LU, n unit right-hand sides |
| `det(A) -> number` | n×n | product of LU pivots with permutation sign (0 if singular) |
| `leastSquares(A, b) -> { x, residual, rank }` | A m×n (m ≥ n); residual = ‖Ax − b‖₂ | Householder QR; minimum-norm (pseudo-inverse) fallback when rank deficient |
| `svd(A) -> { U, S, V }` | thin: U m×k, S (k, descending), V n×k, k = min(m, n) | one-sided (Hestenes) Jacobi |
| `pinv(A, tol?) -> number[][]` | n×m; tol = absolute singular-value cut-off (default max(m,n)·eps·S₀) | SVD |
| `solveTridiagonal(lower, diag, upper, rhs) -> Float64Array` | row i: lower[i]·x[i−1] + diag[i]·x[i] + upper[i]·x[i+1]; off-diagonals of length n or n−1 | Thomas algorithm (no pivoting) |
| `factorBandedSym(bands) -> { n, bw, u }` | bands[i][k] = A[i][i+k], k = 0..bw (ragged rows allowed); u flat `Float64Array` | banded Cholesky A = UᵀU, O(n·bw²); throws if not SPD |
| `solveBandedSymFactored(factor, b) -> Float64Array` | b length n | banded triangular solves |
| `solveBandedSym(bands, b) -> Float64Array` | as above | factor + solve |

## sparse.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `class CooBuilder(n)` | `.add(i, j, v)`, `.toCsr() -> Csr`, `.count` | growable typed arrays; two stable counting sorts, duplicates summed, columns ascending |
| `Csr` (plain object) | `{ n, rowPtr: Int32Array(n+1), colIdx: Int32Array, val: Float64Array }` | — |
| `csrMatVec(A, x, out?) -> Float64Array` | y = A x; `out` must not alias x | row loop |
| `pcg(A, b, { x0?, tol = 1e-10, maxIter?, precond = 'jacobi', omega = 1 }) -> { x, iterations, residual, converged }` | A SPD CSR with both triangles stored; residual = ‖b − Ax‖₂/‖b‖₂; maxIter default max(200, 4n); precond `'jacobi' \| 'ssor' \| 'ic0' \| 'none' \| (r, z) => void` | preconditioned CG; SSOR(ω); IC(0) with Manteuffel shift on breakdown; no throw on non-convergence |
| `class SkylineSym(n, colHeights)` | colHeights[j] = j − (first non-zero row of column j) = stored entries above the diagonal; upper profile in one `Float64Array` | active-column skyline storage |
| `SkylineSym#add(i, j, v)` | accumulates into the upper triangle; **calls with i > j are ignored** (full element-matrix loops assemble each term once); outside the profile throws | — |
| `SkylineSym#get(i, j) -> number` | symmetric read, 0 outside the profile | — |
| `SkylineSym#applyBc(dof, value = 0, rhs?)` | zeroes row/column, unit diagonal; optional in-place rhs update for a prescribed value; before `factor()` only | — |
| `SkylineSym#mulVec(x, out?) -> Float64Array` | K x with the unfactorised matrix | — |
| `SkylineSym#factor() -> this` | in place; sets `.negativePivots` (Sturm count); throws `Error` with `.dof` on a zero pivot | LDLᵀ column reduction, O(n·h²/2) |
| `SkylineSym#solve(b) -> Float64Array` | b untouched; factors on first use | forward / diagonal / back substitution |
| `reverseCuthillMcKee(n, edges) -> Int32Array` | edges: `[i, j][]`; result perm[new] = old | RCM with George–Liu pseudo-peripheral start per component |
| `skylineHeights(n, edges, perm?) -> Int32Array` | column heights for `SkylineSym`, optionally in the renumbered system (perm new -> old) | max over couplings |

## eigen.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `eigSym(A) -> { values, vectors }` | A symmetric n×n; values ascending (`number[]`); vectors[k] unit eigenvector (`number[]`), largest component positive | Householder tridiagonalisation + implicit QL (tred2/tql2) |
| `eigGenSym(K, M, { nModes? }) -> { values, vectors }` | K φ = λ M φ; M SPD matrix or diagonal vector; vectors M-orthonormal; λ = ω² [rad²/s²] for K [N/m], M [kg] | Cholesky reduction to standard form + `eigSym` |
| `subspaceIteration({ n, nModes, solveK, mulM, tol = 1e-8, maxIter = 60, seed = 1 }) -> { values, vectors, iterations, converged }` | solveK(b) = K⁻¹b, mulM(x) = Mx (may reuse their buffers); vectors `Float64Array[]`, M-orthonormal; lowest nModes | subspace iteration with q = min(n, max(2p, p+8)) vectors, Rayleigh–Ritz each sweep, convergence on relative eigenvalue change; rank-deficient M tolerated |
| `eigGeneral(A) -> { re, im }` | real n×n; eigenvalues sorted by real then imaginary part | balancing + elimination to Hessenberg + Francis double-shift QR |

## roots.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `brent(f, a, b, { tol = 1e-12, maxIter = 200 }) -> number` | tol absolute on x; throws `RangeError` if no sign change | Brent (IQI / secant / bisection) |
| `bisect(f, a, b, { tol = 1e-12, maxIter = 200 }) -> number` | as `brent` | bisection |
| `newton(f, x0, { df?, tol = 1e-12, maxIter = 50 }) -> number` | converged on \|dx\| ≤ tol(1+\|x\|); throws `Error` on failure | Newton–Raphson, central-difference derivative if `df` omitted |
| `findBracket(f, x0, { factor = 1.6, maxIter = 60, min?, max?, step? }) -> [a, b]` | a < b with f(a)·f(b) ≤ 0 inside [min, max]; throws `Error` if none | geometric expansion of the end with the smaller \|f\| |
| `newtonSystem(F, x0, { jac?, tol = 1e-10, maxIter = 50, lineSearch = true, fdStep? }) -> { x, iterations, residual, converged }` | F: n -> n; residual = ‖F‖₂; converged on ‖F‖ ≤ tol·max(1, ‖F₀‖) or full step < tol(1+‖x‖) | damped Newton, Armijo backtracking, forward-difference Jacobian, pseudo-inverse step if singular |
| `fixedPoint(g, x0, { tol = 1e-10, maxIter = 200, relax = 1 }) -> { x, iterations, converged }` | scalar or vector x (same type returned) | relaxed Picard iteration |
| `polyRoots(coeffs) -> { re, im }` | real coefficients, highest degree first; sorted by re then im; \|im\| < 6e-8\|z\| snapped to 0 | Aberth–Ehrlich simultaneous iteration on the scaled monic polynomial |

## ode.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `rk4(f, t0, y0, t1, nSteps, onStep?) -> { t, y }` | f(t, y) -> dy; t `Float64Array(nSteps+1)` [s]; y[k] `Float64Array` state at t[k]; `onStep(t, y, k)` may return false to stop | classical RK4, fixed step |
| `rk45(f, t0, y0, t1, { rtol = 1e-6, atol = 1e-9, hMax?, hInit?, maxSteps = 1e6, dense? }) -> { t, y, steps, rejected, nfev }` | atol scalar or per component; output at accepted steps, or exactly at `dense` times | Dormand–Prince 5(4), FSAL, PI step control, 4th-order dense output; throws on maxSteps / step underflow |
| `implicitEuler(f, t0, y0, t1, nSteps, { jac?, tol = 1e-10 }) -> { t, y }` | jac(t, y) -> rows of ∂f/∂y (finite differences if omitted) | backward Euler, simplified Newton with LU |
| `bdf2(f, t0, y0, t1, nSteps, { jac?, tol = 1e-10 }) -> { t, y }` | as `implicitEuler` | BDF2 (first step backward Euler), simplified Newton with LU |
| `newmark({ M, C?, K, f?, u0, v0, dt, nSteps, beta = 1/4, gamma = 1/2, alphaHHT = 0 }) -> { t, u, v, a }` | M [kg], C [N·s/m], K [N/m] n×n rows (or scalars for SDOF); f(t) [N]; u [m], v [m/s], a [m/s²] as `Float64Array[]` of length nSteps+1; t [s] | Newmark-β / HHT-α in acceleration form, one LU; with alphaHHT ≠ 0 the defaults become γ = ½ − α, β = (1 − α)²/4 |
| `newmarkNonlinear({ M, C?, fint, fext?, u0, v0, dt, nSteps, beta = 1/4, gamma = 1/2, tol = 1e-8, maxIter = 25 }) -> { t, u, v, a, iterations, iterationsPerStep }` | fint(u, v) -> { f [N], Kt [N/m], Ct? [N·s/m] }; iterations = total Newton iterations (number); iterationsPerStep `Int32Array` | implicit Newmark with full Newton–Raphson per step; throws `Error` naming the step if not converged |

## quad.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `trapz(y, x \| dx = 1) -> number` | samples y with abscissae x or uniform spacing dx | composite trapezoid |
| `cumtrapz(y, x \| dx = 1) -> Float64Array` | running integral, out[0] = 0 | cumulative trapezoid |
| `simpson(y, dx = 1) -> number` | uniform spacing (an abscissa array is also accepted) | composite Simpson 1/3; last three intervals by the 3/8 rule when the interval count is odd; non-uniform three-point rule for arrays |
| `gaussLegendre(n) -> { x, w }` | nodes ascending on [−1, 1], Σw = 2 | Newton iteration on Legendre recurrence (cached) |
| `integrate(f, a, b, { tol = 1e-9, maxDepth = 40, absTol = 0 }) -> number` | tol relative; a, b may be ±Infinity; re-entrant (nested integrals) | globally adaptive Gauss–Kronrod G7–K15 (QUADPACK QAG strategy), variable transformation for infinite ranges |

## interp.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `interp1(xs, ys, x, { extrapolate = 'clamp' }) -> number` | xs strictly ascending; `'clamp' \| 'linear' \| 'error'` | binary search + linear interpolation |
| `makeLinear(xs, ys, opts) -> (x) => number` | as `interp1`, table validated and copied once | — |
| `makePchip(xs, ys, { extrapolate = 'clamp' }) -> fn` with `fn.d1(x)` | shape-preserving C¹; `'clamp' \| 'linear' \| 'extend' \| 'error'` | Fritsch–Carlson / Fritsch–Butland slopes |
| `makeSpline(xs, ys, { bc = 'natural', d0?, dn?, extrapolate = 'linear' }) -> fn` with `fn.d1(x)`, `fn.d2(x)` | C² cubic; bc `'natural' \| 'not-a-knot' \| 'clamped'` (end slopes d0, dn); extrapolate `'linear' \| 'clamp' \| 'extend' \| 'error'` | tridiagonal system for knot second derivatives (Thomas) |
| `interp2(xs, ys, z, x, y) -> number` | z[iy][ix]; clamped to the grid | bilinear |
| `linspace(a, b, n) -> Float64Array` | n points a..b inclusive | — |
| `logspace(a, b, n) -> Float64Array` | n points 10^a..10^b | — |
| `resample(xs, ys, n) -> { x, y }` | n uniform abscissae over [xs[0], xs[last]], `Float64Array` | linear interpolation, single sweep |
| `decimate(xs, ys, maxPoints) -> { x, y }` | ≤ maxPoints samples (≥ 4), `Float64Array`; end points kept | per-bucket min/max in original order |

## fft.js

| export | shapes / units | algorithm |
| --- | --- | --- |
| `fft(re, im?) -> { re, im }` | X[k] = Σ x[n] e^{−2πikn/N}, any N, `Float64Array` | iterative radix-2; Bluestein chirp-z otherwise |
| `ifft(re, im?) -> { re, im }` | inverse with 1/N | conjugate + forward transform |
| `windowFunction(name, n) -> Float64Array` | `'hann' \| 'hamming' \| 'blackman' \| 'rect'`, DFT-periodic form | — |
| `psdWelch(x, fs, { nfft = 1024, overlap = 0.5, window = 'hann', detrend = 'mean' }) -> { f, psd, nSegments }` | fs [Hz]; f [Hz] 0..fs/2 (nfft/2+1 bins); psd one-sided [units²/Hz], Σpsd·Δf = variance | Welch averaged periodograms, per-segment mean removal; single shorter segment if the record is shorter than nfft |
| `spectralMoments(f, S, orders = [0, 1, 2, 4]) -> { m0, m1, m2, m4 }` | f [Hz], S [units²/Hz]; m_n = ∫ fⁿ S df [units²·Hzⁿ]; order −k keyed `m_k` | trapezoidal rule |
| `synthesizeFromSpectrum(f, S, { duration, dt, seed = 1, randomAmplitude = false }) -> { t, x }` | t = 0, dt, … ≤ duration [s]; x [units]; variance → m0 | random-phase harmonic superposition, a_i = sqrt(2 S_i Δf_i), phases from `rng.js`; phasor recurrence with periodic re-synchronisation |

## special.js

| export | validity | algorithm |
| --- | --- | --- |
| `erf(x)` | all real x | positive-term series for \|x\| < 0.5, else 1 − erfc |
| `erfc(x)` | all real x, relative accuracy kept in the tail | trapezoidal rule with pole correction (Matta & Reichel), h = 1/2 |
| `erfinv(p)` | −1 ≤ p ≤ 1 (±1 -> ±∞) | Winitzki start + Halley iterations (on erfc in the tails) |
| `gamma(x)` | x not a non-positive integer; exact factorials for integer x ≤ 171 | Lanczos (g = 7, 9 terms) + reflection |
| `lgamma(x)` | ln\|Γ(x)\|, x not a non-positive integer | Lanczos + reflection |
| `gammaP(a, x)`, `gammaQ(a, x)` | regularised lower / upper incomplete gamma; a > 0, x ≥ 0 | series for x < a + 1, Lentz continued fraction otherwise |
| `betaI(a, b, x)` | regularised incomplete beta; a, b > 0, 0 ≤ x ≤ 1 | Lentz continued fraction with symmetry switch |
| `besselJ0(x)`, `besselJ1(x)` | all real x | series (\|x\| < 2), Miller backward recurrence (2 ≤ \|x\| < 25), Hankel asymptotic expansion (\|x\| ≥ 25) |
| `besselY0(x)`, `besselY1(x)` | x > 0 | series, Neumann series on the Miller sequence, Hankel asymptotic expansion (same ranges) |
| `besselJn(n, x)` | integer n (negative allowed), all real x | forward recurrence for x > n, Miller backward recurrence otherwise |
| `besselI0(x)`, `besselI1(x)` | all real x (overflow beyond \|x\| ≈ 713) | ascending series (\|x\| ≤ 30), asymptotic expansion above |
| `besselK0(x)`, `besselK1(x)` | x > 0 | series (x ≤ 1), trapezoidal rule on ∫ exp(−x cosh t) cosh(νt) dt otherwise |
