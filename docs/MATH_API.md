# `src/engine/math` and `src/engine/reliability` — fixed public API

These signatures are a contract: implementers must export exactly these names with exactly these shapes;
consumers may rely on them before the files exist. Extra exports are welcome. Dense matrices are arrays of row
arrays (`number[][]`) unless a name says `Flat`. Vectors are `number[]` or `Float64Array` (accept both, return
`Float64Array` for n > 16 where practical, otherwise plain arrays are fine).

## math/linalg.js
```
zeros(n, m?) -> number[][] | number[]            // m omitted -> vector
identity(n) -> number[][]
matMul(A, B) -> number[][]        matVec(A, x) -> number[]        transpose(A) -> number[][]
dot(a, b) -> number               norm(a) -> number               axpy(alpha, x, y) -> number[]  // alpha*x + y
solveDense(A, b) -> number[]                     // LU with partial pivoting; b may be number[] or number[][] (multi-RHS -> number[][])
luFactor(A) -> { lu, piv, n }     luSolve(f, b) -> number[]
cholesky(A) -> number[][]                        // lower L with A = L Lᵀ; throws if not SPD
cholSolve(L, b) -> number[]
inverse(A) -> number[][]          det(A) -> number
leastSquares(A, b) -> { x, residual }            // QR (Householder), A is m×n, m >= n
svd(A) -> { U, S, V }                            // thin SVD, small matrices (one-sided Jacobi)
pinv(A, tol?) -> number[][]
solveTridiagonal(lower, diag, upper, rhs) -> Float64Array        // Thomas
solveBandedSym(bands, b) / factorBandedSym(bands) -> factor; solveBandedSymFactored(factor, b)
    // symmetric positive-definite banded Cholesky. bands[i][k] = A[i][i+k], k = 0..bandwidth
```

## math/sparse.js
```
class CooBuilder(n) { add(i, j, v); toCsr() -> Csr }      // duplicates summed
Csr = { n, rowPtr: Int32Array, colIdx: Int32Array, val: Float64Array }
csrMatVec(A, x, out?) -> Float64Array
pcg(A, b, { x0?, tol = 1e-10, maxIter?, precond = 'jacobi' | 'ssor' | 'ic0' }) -> { x, iterations, residual, converged }
class SkylineSym(n, colHeights: Int32Array) { add(i, j, v); factor(); solve(b) -> Float64Array; applyBc(dof) }
    // symmetric skyline LDLᵀ direct solver for FEM (upper profile storage); factor() throws on zero pivot with dof index
reverseCuthillMcKee(n, edges: [i, j][]) -> Int32Array     // permutation new->old
```

## math/eigen.js
```
eigSym(A) -> { values: number[] (ascending), vectors: number[][] (columns as vectors[k]) }   // dense symmetric, Jacobi or QL
eigGenSym(K, M, { nModes? }) -> { values, vectors }   // dense K φ = λ M φ, M SPD (or diagonal incl. zeros handled by static condensation is NOT required)
subspaceIteration({ n, nModes, solveK(b) -> x, mulM(x) -> y, tol = 1e-8, maxIter = 60, seed? }) -> { values, vectors: Float64Array[] , iterations, converged }
    // lowest modes of K φ = λ M φ for large systems, given a factorised-K solve and an M product; vectors M-orthonormal
eigGeneral(A) -> { re: number[], im: number[] }       // unsymmetric real matrix, eigenvalues only (Hessenberg + shifted QR)
```

## math/roots.js
```
brent(f, a, b, { tol = 1e-12, maxIter = 200 }) -> number           // throws if no sign change
bisect(f, a, b, opts) -> number
newton(f, x0, { df?, tol, maxIter }) -> number                      // numerical derivative if df omitted
findBracket(f, x0, { factor = 1.6, maxIter = 60, min?, max? }) -> [a, b]
newtonSystem(F, x0, { jac?, tol = 1e-10, maxIter = 50, lineSearch = true, fdStep? }) -> { x, iterations, residual, converged }
    // F: number[] -> number[]; damped Newton with finite-difference Jacobian when jac omitted
fixedPoint(g, x0, { tol, maxIter, relax = 1 }) -> { x, iterations, converged }
polyRoots(coeffs) -> { re, im }                                     // highest degree first
```

## math/ode.js
```
rk4(f, t0, y0, t1, nSteps, onStep?) -> { t: Float64Array, y: Float64Array[] }      // f(t, y) -> dy (number[]); y[k] is the state at t[k]
rk45(f, t0, y0, t1, { rtol = 1e-6, atol = 1e-9, hMax?, hInit?, maxSteps = 1e6, dense?: number[] (output times) }) -> { t, y, steps, rejected }
implicitEuler / bdf2 (f, t0, y0, t1, nSteps, { jac?, tol }) -> { t, y }             // stiff systems, Newton inside
newmark({ M, C, K, f(t) -> number[], u0, v0, dt, nSteps, beta = 0.25, gamma = 0.5, alphaHHT = 0 }) -> { t, u: Float64Array[], v, a }
    // linear MDOF, dense matrices; alphaHHT in [-1/3, 0] gives Hilber–Hughes–Taylor
newmarkNonlinear({ M, C, fint(u, v) -> { f, Kt }, fext(t), u0, v0, dt, nSteps, beta, gamma, tol, maxIter }) -> { t, u, v, a, iterations }
```

## math/quad.js
```
trapz(y, x | dx) -> number        cumtrapz(y, x | dx) -> Float64Array        simpson(y, dx) -> number
gaussLegendre(n) -> { x: number[], w: number[] }   // on [-1, 1], n <= 64
integrate(f, a, b, { tol = 1e-9, maxDepth = 40 }) -> number      // adaptive Gauss–Kronrod or Simpson
```

## math/interp.js
```
interp1(xs, ys, x, { extrapolate = 'clamp' | 'linear' | 'error' }) -> number    // xs ascending
makeLinear(xs, ys, opts) -> (x) => number      makePchip(xs, ys) -> (x) => number      makeSpline(xs, ys, { bc = 'natural' }) -> fn with fn.d1(x), fn.d2(x)
interp2(xs, ys, z /* z[iy][ix] */, x, y) -> number        // bilinear, clamped
linspace(a, b, n) -> Float64Array        logspace(a, b, n) -> Float64Array        // logspace(a,b): 10^a..10^b
resample(xs, ys, n) -> { x, y }          decimate(xs, ys, maxPoints) -> { x, y }  // min/max preserving
```

## math/fft.js
```
fft(re, im?) -> { re, im }        ifft(re, im) -> { re, im }        // any length (radix-2 fast path, Bluestein otherwise)
psdWelch(x, fs, { nfft = 1024, overlap = 0.5, window = 'hann' }) -> { f, psd }     // one-sided, units²/Hz
spectralMoments(f, S, orders = [0, 1, 2, 4]) -> { m0, m1, m2, m4 }                  // f in Hz
synthesizeFromSpectrum(f, S, { duration, dt, seed }) -> { t, x }                    // random-phase sum, deterministic for a seed
```

## math/special.js
```
erf(x), erfc(x), erfinv(p), gamma(x), lgamma(x), gammaP(a, x) /* regularised lower */, gammaQ(a, x), betaI(a, b, x) /* regularised */,
besselJ0, besselJ1, besselY0, besselY1, besselJn(n, x), besselI0, besselI1, besselK0, besselK1
```

## math/rng.js
```
makeRng(seed: number) -> { next() -> [0,1), normal() -> N(0,1), int(n) -> 0..n-1, shuffle(arr) }
```

## math/stats.js
```
mean, variance, std, min, max, median, percentile(arr, p /* 0..100 */), skewness, kurtosis, histogram(arr, bins) -> { edges, counts },
covariance(a, b), correlation(a, b), rankCorrelation(a, b)
dist(spec) -> { pdf(x), cdf(x), icdf(p), mean, std, sample(rng) }
   // spec: { dist: 'normal', mean, std } | { dist: 'lognormal', mean, std /* of X, not ln X */ } | { dist: 'uniform', min, max } |
   //       { dist: 'triangular', min, mode, max } | { dist: 'pert', min, mode, max, lambda = 4 } | { dist: 'beta', a, b, min = 0, max = 1 } |
   //       { dist: 'weibull', scale, shape, loc = 0 } | { dist: 'gumbel', loc, scale } | { dist: 'gev', loc, scale, shape } |
   //       { dist: 'gpd', loc, scale, shape } | { dist: 'exponential', rate } | { dist: 'gamma', shape, scale } | { dist: 'rayleigh', sigma } |
   //       { dist: 'deterministic', value }
   // convenience: any spec may give { mean, cov } instead of { mean, std }
linearRegression(x, y) -> { slope, intercept, r2, seSlope, seIntercept }
polyfit(x, y, degree) -> number[] /* highest degree first */       polyval(c, x)
multipleRegression(X /* rows */, y) -> { beta, r2, residuals, se }  // intercept added automatically as beta[0]
metrics(pred, obs) -> { rmse, mae, bias, nrmse, r2, maxAbs, n }     // r2 = 1 - SSE/SST; bias = mean(pred - obs)
ksTest(sample, cdf) -> { D, p }
```

## math/optim.js
```
goldenSection(f, a, b, { tol }) -> { x, f }
nelderMead(f, x0, { bounds?: [lo, hi][], tol = 1e-8, maxIter = 2000, step? }) -> { x, f, iterations, converged }
levenbergMarquardt(residuals /* x -> number[] */, x0, { bounds?, tol = 1e-10, maxIter = 200, jac? }) -> { x, cost, iterations, converged, covariance /* (JᵀJ)⁻¹ σ² */, residuals }
bfgs(f, x0, { grad?, bounds?, tol, maxIter }) -> { x, f, iterations, converged }
differentialEvolution(f, bounds, { pop = 15*n, F = 0.7, CR = 0.9, maxGen = 200, seed, tol }) -> { x, f, generations, history }
particleSwarm(f, bounds, { particles, maxIter, seed, w, c1, c2 }) -> { x, f, history }
geneticAlgorithm(f, bounds, { pop, generations, seed, mutation, crossover, integer?: boolean[] }) -> { x, f, history }
simulatedAnnealing(f, x0, bounds, { T0, cooling, maxIter, seed }) -> { x, f, history }
nsga2(fs /* x -> number[] objectives (minimise) */, bounds, { pop, generations, seed }) -> { front: [{ x, f }] }
linprog({ c, Aub?, bub?, Aeq?, beq?, bounds? }) -> { x, f, status: 'optimal' | 'infeasible' | 'unbounded' }     // minimise cᵀx, two-phase simplex
milp({ c, Aub, bub, Aeq, beq, bounds, integer: boolean[] }) -> { x, f, status, nodes }                           // branch and bound on linprog
goalSeek(f, target, a, b) -> number
```

## math/timeseries.js
```
acf(x, maxLag), pacf(x, maxLag), difference(x, d), movingAverage(x, w), expSmoothing(x, { alpha, beta?, gamma?, season? }) -> { fitted, forecast(h) }
arima(x, { p, d, q }) -> { phi, theta, sigma2, aic, forecast(h) -> { mean, lo95, hi95 } }      // conditional sum of squares
kalmanFilter({ F, H, Q, R, x0, P0 }, observations) -> { x: number[][], P, innovations }
rainfallLikeMarkov?  // not required
```

## reliability/sampling.js
```
monteCarlo(g /* x: object -> number */, vars /* { name: distSpec } */, { n = 10000, seed = 1, corr? }) -> { samples: Float64Array (g values), inputs?: , mean, std, pf (P[g<0]), ciPf }
latinHypercube(vars, n, { seed, corr? /* name pairs -> rho, Iman–Conover */ }) -> { [name]: Float64Array }
propagate(model /* x: object -> { [out]: number } */, vars, { n, seed, method = 'lhs' | 'mc' }) -> { outputs: { [out]: Float64Array }, inputs: { [name]: Float64Array }, stats: { [out]: { mean, std, p10, p50, p90, min, max } } }
sobolIndices(model /* x: object -> number */, vars, { n = 1024, seed }) -> { first: { [name]: number }, total: { [name]: number } }     // Saltelli
tornado(model, vars /* { name: { low, base, high } } */) -> [{ name, low, high, swing }] sorted by swing
```

## reliability/form.js
```
form(g /* x: object -> number, failure when g < 0 */, vars, { corr?, tol = 1e-6, maxIter = 100 }) -> { beta, pf, alpha: { [name]: number }, designPoint: { [name]: number }, iterations, converged }
sorm(g, vars, opts) -> { beta, pfBreitung, betaGeneralised, curvatures }
importanceSampling(g, vars, { designPoint, n, seed }) -> { pf, cov }
```

## reliability/extreme.js
```
fitGumbel(x, { method = 'mle' | 'mom' | 'lmom' }) -> { loc, scale, se?, logLik }
fitWeibull(x, { threeParam = false }) -> { scale, shape, loc }
fitGev(x) -> { loc, scale, shape }            // L-moments (Hosking)
fitGpd(excess) -> { scale, shape }            // L-moments / MLE
fitLognormal(x) -> { mu, sigma }             fitNormal(x)
annualMaxima(times /* ms epoch or ISO strings */, values) -> { years, maxima }
peaksOverThreshold(times, values, { threshold, minSeparationHours = 48 }) -> { peaks, times, ratePerYear }
returnValue(fit /* { dist, ...params } */, T /* years */, { ratePerYear = 1 }) -> number
bootstrapReturnValue(x, fitFn, T, { n = 500, seed }) -> { lo, median, hi }
iformContour({ hs: Float64Array, tp: Float64Array }, T /* years */, { seaStateHours = 1, nPoints = 72 }) -> { hs: number[], tp: number[], fit }
    // 3-p Weibull marginal for Hs, conditional lognormal Tp | Hs with fitted mean/std functions
scatterDiagram(hs, tp, { hsBin = 0.5, tpBin = 1 }) -> { hsEdges, tpEdges, counts: number[][], prob: number[][] }
```

## reliability/surrogate.js
```
gaussianProcess(X /* rows */, y, { kernel = 'rbf' | 'matern52', lengthScale?, noise = 1e-8, optimise = true }) -> { predict(x) -> { mean, std }, logLik, params }
metropolis(logPost, x0, { n, stepSizes, seed, burn }) -> { chain: number[][], acceptance }
polynomialChaos?  // optional
bayesOptimise(f, bounds, { nInit, nIter, seed }) -> { x, f, history }      // GP + expected improvement
```
