# Statistics, optimisation, time series and reliability engines

One line per export: signature → return shape — algorithm / reference. Signatures of the names listed in
`docs/MATH_API.md` are exactly as specified there; fields and exports marked *(extra)* are additions.
Vectors may be plain or typed arrays. Random variables are `{ name: distSpec }`; models take `{ name: value }`.

## `src/engine/math/rng.js`

- `makeRng(seed = 1)` → `{ next() ∈ [0,1), normal(), int(n) ∈ 0..n-1, shuffle(arr) }` + *(extra)* `uniform(a, b)`, `uint32()`, `exponential(rate)`, `fork()`, `seed` — sfc32 generator seeded through splitmix32, 53-bit uniforms; normal deviates by the Marsaglia polar method; Fisher–Yates shuffle.
- `splitmix32(seed)` → `() => uint32` *(extra)* — seed expander.

## `src/engine/math/stats.js`

- `mean(arr)`, `variance(arr, population = false)`, `std(arr, population = false)`, `min(arr)`, `max(arr)` → number — definitions; variance/std use the n−1 divisor unless `population`.
- `median(arr)`, `percentile(arr, p /* 0..100 */)` → number — linear interpolation of order statistics (Hyndman & Fan type 7). `percentileSorted(sorted, p)` *(extra)* skips the sort.
- `skewness(arr)` → g1 = m3/m2^1.5; `kurtosis(arr, excess = true)` → **excess** kurtosis m4/m2² − 3 (0 for a normal sample; pass `false` for the Pearson value).
- `histogram(arr, bins /* count or edges */)` → `{ edges, counts }` + *(extra)* `centres`, `density` — equal-width bins over [min, max], last bin closed on the right.
- `covariance(a, b)`, `correlation(a, b)` (Pearson), `rankCorrelation(a, b)` (Spearman, average ranks for ties) → number; `ranks(arr)` → Float64Array *(extra)*.
- `normPdf(z)`, `normCdf(z)`, `normInv(p)` *(extra)* — Φ by erfc; Φ⁻¹ by Acklam's rational approximation + Halley refinement (tails included).
- `dist(spec)` → `{ pdf(x), cdf(x), icdf(p), mean, std, sample(rng) }` + *(extra)* `variance`, `support: [lo, hi]`, `params` (resolved native parameters), `name`, `spec`, `fromZ(z)` (x = F⁻¹(Φ(z))), `toZ(x)`.
  - `normal { mean, std }`; `lognormal { mean, std }` — moments of X itself (σ_ln² = ln(1+V²), μ_ln = ln m − σ_ln²/2), or `{ mu, sigma }` of ln X;
    `uniform { min, max }`; `triangular { min, mode, max }`; `pert { min, mode, max, lambda = 4 }` — beta-PERT, Beta(1+λ(mode−min)/(max−min), 1+λ(max−mode)/(max−min)) on [min, max];
    `beta { a, b, min = 0, max = 1 }`; `weibull { scale, shape, loc = 0 }`; `gumbel { loc, scale }` (maxima); `gev { loc, scale, shape }`; `gpd { loc, scale, shape }`;
    `exponential { rate }`; `gamma { shape, scale }`; `rayleigh { sigma }`; `deterministic { value }`.
  - **GEV/GPD shape sign:** ξ > 0 heavy (Fréchet/Pareto) tail, ξ < 0 bounded upper tail, ξ = 0 Gumbel/exponential (Coles 2001); Hosking's k = −ξ. `fitGev`/`fitGpd` return the same ξ.
  - `{ mean, cov }` may replace `{ mean, std }`; `uniform`, `triangular` (symmetric), `beta`, `weibull`, `gumbel`, `gamma` also accept `{ mean, std | cov }` (moment matching), `exponential`/`rayleigh`/`deterministic` accept `{ mean }`.
  - icdf: closed form, except beta/pert/gamma (safeguarded Newton on the regularised incomplete beta/gamma function). Sampling: inverse transform; normal/lognormal polar; gamma and beta by Marsaglia–Tsang.
  - Invalid parameters or p ∉ [0, 1] throw `RangeError` naming the argument. `mean`/`std` are `Infinity` where the moment does not exist.
- `DISTRIBUTIONS` (name list), `EULER_GAMMA` *(extra)*.
- `linearRegression(x, y)` → `{ slope, intercept, r2, seSlope, seIntercept }` + *(extra)* `n`, `sigma` — ordinary least squares (Draper & Smith ch. 1).
- `polyfit(x, y, degree)` → `number[]` highest degree first — Householder QR on centred/scaled abscissae; `polyval(c, x)` → number (or array for array x) — Horner.
- `multipleRegression(X /* rows */, y)` → `{ beta /* [0] = intercept */, r2, residuals, se }` + *(extra)* `sigma`, `adjR2` — QR least squares, se = √(diag((AᵀA)⁻¹)·SSE/(n−p)).
- `metrics(pred, obs)` → `{ rmse, mae, bias, nrmse, r2, maxAbs, n }` — bias = mean(pred − obs); nrmse = rmse/(max(obs) − min(obs)); r2 = 1 − SSE/SST.
- `ksTest(sample, cdf)` → `{ D, p }` + *(extra)* `n` — one-sample Kolmogorov–Smirnov, asymptotic p-value with Stephens' correction.

## `src/engine/math/optim.js` (all minimise)

- `goldenSection(f, a, b, { tol })` → `{ x, f }` + *(extra)* `iterations` — golden-section search (unimodal f).
- `nelderMead(f, x0, { bounds?, tol = 1e-8, maxIter = 2000, step? })` → `{ x, f, iterations, converged }` + *(extra)* `evaluations` — downhill simplex (1, 2, ½, ½); trial points projected onto the bounds.
- `levenbergMarquardt(residuals, x0, { bounds?, tol = 1e-10, maxIter = 200, jac? })` → `{ x, cost /* ½Σr² */, iterations, converged, covariance /* (JᵀJ)⁻¹σ², σ² = Σr²/(m−n) */, residuals }` + *(extra)* `sigma2`, `se` — Marquardt-scaled damping, forward-difference Jacobian unless `jac`, projection onto bounds.
- `bfgs(f, x0, { grad?, bounds?, tol = 1e-8, maxIter = 500 })` → `{ x, f, iterations, converged }` + *(extra)* `gradNorm`, `evaluations` — BFGS inverse-Hessian update, Armijo backtracking, projected active-set treatment of bounds, central-difference gradient unless `grad`.
- `differentialEvolution(f, bounds, { pop = 15n, F = 0.7, CR = 0.9, maxGen = 200, seed = 1, tol = 1e-12 })` → `{ x, f, generations, history /* best f per generation */ }` — DE/rand/1/bin (Storn & Price 1997), reflection at bounds.
- `particleSwarm(f, bounds, { particles = max(30, 10n), maxIter = 200, seed = 1, w = 0.7298, c1 = c2 = 1.49618 })` → `{ x, f, history }` — global-best PSO with constriction coefficients (Clerc & Kennedy 2002), absorbing walls.
- `geneticAlgorithm(f, bounds, { pop = max(40, 10n), generations = 200, seed = 1, mutation, crossover = 0.9, integer? })` → `{ x, f, history }` — real-coded GA: tournament selection, SBX crossover, polynomial mutation, 2-elitism; `integer[j]` keeps variable j integral.
- `simulatedAnnealing(f, x0, bounds, { T0, cooling, maxIter = 20000, seed = 1 })` → `{ x, f, history }` — Metropolis acceptance, geometric cooling per iteration (default reaches 1e-10·T0), adaptive Gaussian step (Corana et al. 1987); T0 defaults to the mean |Δf| of 20 trial moves.
- `nsga2(fs, bounds, { pop = 100, generations = 200, seed = 1 })` → `{ front: [{ x, f }] }` sorted by f[0] — fast non-dominated sort + crowding distance, crowded tournament, SBX (η = 15), polynomial mutation (η = 20) (Deb et al. 2002). `nonDominatedSort(F)`, `crowdingDistance(F, front, out)` *(extra)*.
- `linprog({ c, Aub?, bub?, Aeq?, beq?, bounds? })` → `{ x, f, status: 'optimal' | 'infeasible' | 'unbounded' }` — dense two-phase primal simplex with Bland's rule; default bounds x ≥ 0; `null`/±Infinity bounds allowed (free variables are split). `x` is `null` when not optimal.
- `milp({ c, Aub, bub, Aeq, beq, bounds, integer, maxNodes = 1e5 })` → `{ x, f, status, nodes }` — depth-first branch and bound on `linprog`, most-fractional branching; *(extra)* status `'node-limit'`.
- `goalSeek(f, target, a, b)` → number — Brent's method on f − target (throws `RangeError` without a sign change).
- All stochastic methods are deterministic for a given `seed`; `history` arrays hold the best objective so far.

## `src/engine/math/timeseries.js`

- `acf(x, maxLag)` → `number[maxLag+1]` (r₀ = 1) — biased sample autocorrelation (Box–Jenkins §2.1).
- `pacf(x, maxLag)` → `number[maxLag+1]` (first entry 1) — Durbin–Levinson recursion; `yuleWalker(x, p)` → φ *(extra)*.
- `difference(x, d = 1)` → `number[n−d]`; `movingAverage(x, w)` → `number[n]` — centred window, truncated symmetrically at the ends, 2×w average for even w.
- `expSmoothing(x, { alpha, beta?, gamma?, season?, trend?, seasonal = 'additive' | 'multiplicative' })` → `{ fitted, forecast(h) → number[] }` + *(extra)* `level`, `trend`, `seasonals`, `alpha`, `beta`, `gamma`, `sse`, `rmse` — simple / Holt / Holt–Winters; omitted smoothing constants are fitted by minimising the one-step SSE.
- `arima(x, { p, d, q, includeMean? })` → `{ phi, theta, sigma2, aic, forecast(h) → { mean, lo95, hi95 } }` + *(extra)* `mean`, `residuals`, `se`, `converged`, `stationary`, `forecast().se` — conditional sum of squares minimised by Levenberg–Marquardt; model (w−μ) = Σφᵢ(w₋ᵢ−μ) + e + Σθⱼe₋ⱼ on w = ∇ᵈx; intervals from ψ-weights.
- `kalmanFilter({ F, H, Q, R, x0, P0 }, observations)` → `{ x: number[][], P: number[][][], innovations }` + *(extra)* `xPred`, `S`, `logLik` — predict/update with Joseph-form covariance; `null`/NaN observations are treated as missing.

## `src/engine/reliability/sampling.js`

- `monteCarlo(g, vars, { n = 10000, seed = 1, corr?, keepInputs? })` → `{ samples: Float64Array, inputs?, mean, std, pf, ciPf: [lo, hi] }` + *(extra)* `covPf`, `nFail`, `n`, `beta` — crude Monte Carlo; pf = P[g < 0] with the 95 % Wilson interval; Nataf transformation when `corr` is given.
- `latinHypercube(vars, n, { seed = 1, corr? })` → `{ [name]: Float64Array }` — one sample per equiprobable stratum (McKay et al. 1979); `corr` imposes rank correlation by Iman–Conover (1982).
- `propagate(model, vars, { n = 1000, seed = 1, method = 'lhs' | 'mc', corr? })` → `{ outputs, inputs, stats: { [out]: { mean, std, p10, p50, p90, min, max } } }` + *(extra)* `p05`, `p95`.
- `sobolIndices(model, vars, { n = 1024, seed = 1 })` → `{ first, total }` + *(extra)* `variance`, `mean`, `evaluations` — Saltelli scheme, n(d+2) runs on a randomised Halton base sample; first order by Saltelli et al. (2010), total effect by Jansen (1999); independent inputs.
- `tornado(model, vars /* { name: { low, base, high } } */)` → `[{ name, low, high, swing }]` sorted by swing + *(extra)* `base`, `xLow`, `xHigh`.
- *(extra)* `resolveVars(vars)`, `corrMatrix(names, corr)`, `natafRho(di, dj, rho)`, `makeTransform(names, dists, R)` → `{ toX(u), toU(x), L, R0 }`, `scrambledHalton(n, dim, rng)`.
- `corr` forms: `{ 'a,b': rho }` (separators `,` `|` `:` `~`), `{ a: { b: rho } }`, `[[a, b, rho]]` / `[{ a, b, rho }]`, or a full matrix ordered like `Object.keys(vars)`.

## `src/engine/reliability/form.js`

- `form(g, vars, { corr?, tol = 1e-6, maxIter = 100, fdStep?, u0? })` → `{ beta, pf, alpha: { [name] }, designPoint: { [name] }, iterations, converged }` + *(extra)* `u`, `importance` (α²), `gradient`, `evaluations` — HL-RF iteration with merit-function line search (improved HL-RF) in standard normal space; isoprobabilistic transform of arbitrary marginals; Nataf model (Liu & Der Kiureghian 1986) for `corr`. Sign: u* = β·α, so α < 0 for resistance-type and α > 0 for load-type variables.
- `sorm(g, vars, opts)` → `{ beta, pfBreitung, betaGeneralised, curvatures }` + *(extra)* `pfForm`, `pfHohenbichler`, `valid`, `form` — finite-difference Hessian at the design point, principal curvatures, Breitung (1984).
- `importanceSampling(g, vars, { designPoint, n = 2000, seed = 1, corr? })` → `{ pf, cov }` + *(extra)* `beta`, `nFail`, `n` — unit-variance normal sampling density centred on the design point (found by FORM when omitted).

## `src/engine/reliability/extreme.js`

Every fit result also carries `dist` (so it can be passed to `returnValue` / `dist`), `logLik` and `n`.

- `fitGumbel(x, { method = 'mle' | 'mom' | 'lmom' })` → `{ loc, scale, se? /* { loc, scale }, MLE only */, logLik }`.
- `fitWeibull(x, { threeParam = false })` → `{ scale, shape, loc }` — MLE (Newton on the shape equation); three-parameter fit by profile likelihood over the location.
- `fitGev(x)` → `{ loc, scale, shape /* ξ */ }` — L-moments (Hosking, Wallis & Wood 1985), shape equation solved exactly.
- `fitGpd(excess, { method = 'lmom' | 'mle' })` → `{ scale, shape /* ξ */ }` (+ `loc: 0`; set `loc` to the threshold for return values).
- `fitLognormal(x)` → `{ mu, sigma }` of ln X (+ `mean`, `std` of X); `fitNormal(x)` → `{ mean, std }` (+ aliases `mu`, `sigma`).
- `lMoments(x)` → `{ l1, l2, l3, l4, t3, t4 }` *(extra)* — unbiased probability-weighted moments (Hosking 1990).
- `annualMaxima(times, values)` → `{ years, maxima }` + *(extra)* `counts`, `timesOfMax` — UTC calendar years; times in ms epoch, ISO strings or Dates.
- `peaksOverThreshold(times, values, { threshold, minSeparationHours = 48 })` → `{ peaks, times, ratePerYear }` + *(extra)* `excesses`, `threshold`, `years` — runs declustering, one peak per cluster.
- `returnValue(fit /* { dist, ...params } */, T, { ratePerYear = 1 })` → number — F⁻¹(1 − 1/(λT)).
- `bootstrapReturnValue(x, fitFn, T, { n = 500, seed = 1, ratePerYear = 1, level = 0.95 })` → `{ lo, median, hi }` + *(extra)* `estimate`, `nValid`, `values` — percentile bootstrap.
- `iformContour({ hs, tp }, T, { seaStateHours = 1, nPoints = 72, hsBin?, minCount? })` → `{ hs: number[], tp: number[], fit }` + *(extra)* `hsMax`, `tpAtHsMax`; `fit = { weibull: { scale, shape, loc }, mu: { a0, a1, a2 }, sigma: { b0, b1, b2 }, beta, nSeaStates, bins }` — 3-parameter Weibull Hs marginal, lognormal Tp | Hs with μ(h) = a0 + a1·h^a2 and σ(h) = b0 + b1·exp(b2·h) of ln Tp; contour = image of the β-circle (Winterstein et al. 1993). Points start at θ = 0 (largest Hs) and are not repeated at the end.
- `scatterDiagram(hs, tp, { hsBin = 0.5, tpBin = 1 })` → `{ hsEdges, tpEdges, counts /* [iHs][iTp] */, prob }` + *(extra)* `n`.

## `src/engine/reliability/surrogate.js`

- `gaussianProcess(X, y, { kernel = 'rbf' | 'matern52', lengthScale?, noise = 1e-8, optimise = true })` → `{ predict(x) → { mean, std }, logLik, params }` + *(extra)* `n`, `dim` — constant-mean GP on scaled inputs, Cholesky solves, signal variance profiled, length scales (one per input for d ≤ 8) by maximising the concentrated marginal likelihood with `nelderMead`; `predict` accepts one point or an array of points; `params = { kernel, lengthScale[], signalVariance, noise, mean }`.
- `metropolis(logPost, x0, { n = 5000, stepSizes, seed = 1, burn = n/5, thin = 1 })` → `{ chain: number[][], acceptance }` + *(extra)* `logPost`, `mean`, `std` — random-walk Metropolis, Gaussian proposals.
- `bayesOptimise(f, bounds, { nInit = max(5, 2d+1), nIter = 25, seed = 1, kernel = 'matern52' })` → `{ x, f, history /* best so far per evaluation */ }` + *(extra)* `X`, `y` — GP surrogate + expected improvement (Jones et al. 1998). `expectedImprovement(mean, std, fBest)` *(extra)*.
- Not implemented: `polynomialChaos` (optional in the API).
