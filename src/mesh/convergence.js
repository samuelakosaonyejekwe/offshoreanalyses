/**
 * Discretisation-error estimation and mesh-convergence studies.
 *
 * Implements the procedure of Celik, Ghia, Roache, Freitas, Coleman & Raad, "Procedure for estimation and
 * reporting of uncertainty due to discretization in CFD applications", J. Fluids Eng. 130 (2008) 078001 —
 * observed order from three grids with non-constant refinement ratio, Richardson extrapolation and the Grid
 * Convergence Index of Roache ("Perspective: a method for uniform reporting of grid refinement studies",
 * J. Fluids Eng. 116 (1994)) — together with the least-squares variant of Eça & Hoekstra ("A procedure for the
 * estimation of the numerical uncertainty of CFD calculations based on grid refinement studies", J. Comput.
 * Phys. 262 (2014)) for more than three grids, and error norms for manufactured solutions.
 *
 * Index convention (Celik et al.): 1 = finest grid, 2 = medium, 3 = coarsest; h is the representative cell
 * size, r21 = h2 / h1 > 1, r32 = h3 / h2 > 1. All quantities are dimensionless unless they carry the units of
 * the solution value.
 */
import { goldenSection } from '../engine/math/optim.js';

const finite = (v, name) => { if (!Number.isFinite(v)) throw new RangeError(`${name}: must be a finite number`); return v; };
const ratio = (v, name) => { if (!(v > 1) || !Number.isFinite(v)) throw new RangeError(`${name}: refinement ratio must be > 1`); return v; };

/**
 * Convergence behaviour of three solutions from R = (f1 − f2) / (f2 − f3)
 * (Stern et al., "Comprehensive approach to verification and validation of CFD simulations", J. Fluids Eng.
 * 123 (2001)): 0 < R < R_lim monotone convergence, −1 < R < 0 oscillatory convergence, R ≥ R_lim monotone
 * divergence, R ≤ −1 oscillatory divergence. For a constant refinement ratio R_lim = 1; for unequal ratios the
 * limit of R as the order tends to zero is ln r21 / ln r32, which is used instead.
 * @param {number} f1 finest-grid value
 * @param {number} f2 medium-grid value
 * @param {number} f3 coarsest-grid value
 * @param {{ r21?: number, r32?: number, tol?: number }} [opts] refinement ratios (default equal); relative
 *   difference below which two solutions count as equal (default 1e-12)
 * @returns {{ type: 'monotone'|'oscillatory'|'divergent'|'oscillatory-divergent'|'converged'|'indeterminate', R: number }}
 *   'converged': the two finest solutions agree; 'indeterminate': the two coarsest agree but the finest differs.
 */
export function classify(f1, f2, f3, { r21 = 2, r32 = 2, tol = 1e-12 } = {}) {
  finite(f1, 'f1'); finite(f2, 'f2'); finite(f3, 'f3'); ratio(r21, 'r21'); ratio(r32, 'r32');
  const scale = Math.max(Math.abs(f1), Math.abs(f2), Math.abs(f3), 1e-300);
  const e21 = f2 - f1, e32 = f3 - f2;
  if (Math.abs(e21) <= tol * scale) return { type: 'converged', R: 0 };
  if (Math.abs(e32) <= tol * scale) return { type: 'indeterminate', R: Infinity };
  const R = e21 / e32;
  const lim = Math.log(r21) / Math.log(r32);
  return { type: R > 0 ? (R < lim * (1 - 1e-12) ? 'monotone' : 'divergent') : (R > -1 ? 'oscillatory' : 'oscillatory-divergent'), R };
}

/**
 * Observed (apparent) order of accuracy from three grids, Celik et al. (2008) Eq. (3):
 *   p = |ln|ε32 / ε21| + q(p)| / ln r21,  q(p) = ln[(r21^p − s) / (r32^p − s)],  s = sign(ε32 / ε21),
 * with ε32 = f3 − f2 and ε21 = f2 − f1, solved by fixed-point iteration (q = 0 for a constant ratio).
 * @param {{ f1: number, f2: number, f3: number, r21: number, r32: number, tol?: number, maxIter?: number }} o
 * @returns {{ p: number, s: number, iterations: number, converged: boolean }} p is NaN when either difference
 *   vanishes (order undefined)
 */
export function observedOrder({ f1, f2, f3, r21, r32, tol = 1e-10, maxIter = 200 }) {
  finite(f1, 'f1'); finite(f2, 'f2'); finite(f3, 'f3'); ratio(r21, 'r21'); ratio(r32, 'r32');
  const e32 = f3 - f2, e21 = f2 - f1;
  if (e21 === 0 || e32 === 0) return { p: NaN, s: 0, iterations: 0, converged: false };
  const s = Math.sign(e32 / e21);
  const base = Math.log(Math.abs(e32 / e21));
  const l21 = Math.log(r21);
  if (Math.abs(r21 - r32) < 1e-12 * r21) return { p: Math.abs(base) / l21, s, iterations: 0, converged: true };
  const rhs = (p) => {
    const num = r21 ** p - s, den = r32 ** p - s;
    return num / den > 0 ? Math.abs(base + Math.log(num / den)) / l21 : NaN;
  };
  // the constant-ratio estimate starts the iteration; it is zero when the two differences happen to be equal,
  // and p = 0 is a singular point of q(p), so start from first order then
  let p = Math.abs(base) / l21;
  if (!(p > 0.05)) p = 1;
  for (let it = 1; it <= maxIter; it++) {
    const next = rhs(p);
    if (!(next > 0) || !Number.isFinite(next)) break;
    if (Math.abs(next - p) <= tol * Math.max(1, next)) return { p: next, s, iterations: it, converged: true };
    // light under-relaxation keeps the iteration stable for oscillatory data
    p = s < 0 ? 0.5 * (p + next) : next;
  }
  // fixed-point iteration failed: bracket a root of p − rhs(p) on a logarithmic scan and bisect
  const g = (x) => x - rhs(x);
  let lo = 0.01, glo = g(lo);
  for (let i = 1; i <= 400; i++) {
    const hi = 0.01 * 2000 ** (i / 400), ghi = g(hi);
    if (Number.isFinite(glo) && Number.isFinite(ghi) && glo * ghi <= 0) {
      let a = lo, b = hi, ga = glo;
      for (let k = 0; k < 200 && b - a > tol * Math.max(1, b); k++) {
        const m = (a + b) / 2, gm = g(m);
        if (!Number.isFinite(gm)) break;
        if (ga * gm <= 0) b = m; else { a = m; ga = gm; }
      }
      return { p: (a + b) / 2, s, iterations: maxIter, converged: true };
    }
    lo = hi; glo = ghi;
  }
  return { p: NaN, s, iterations: maxIter, converged: false };
}

/**
 * Richardson extrapolation to zero cell size, Celik et al. (2008) Eq. (4): f_ext = (r^p f1 − f2) / (r^p − 1).
 * @param {{ f1: number, f2: number, r21: number, p: number }} o finer and coarser values, their size ratio and
 *   the order of accuracy (observed or formal)
 * @returns {number} extrapolated value (units of f)
 */
export function richardson({ f1, f2, r21, p }) {
  finite(f1, 'f1'); finite(f2, 'f2'); ratio(r21, 'r21');
  if (!(p > 0)) throw new RangeError('p: order must be > 0');
  const rp = r21 ** p;
  return (rp * f1 - f2) / (rp - 1);
}

/**
 * Grid Convergence Index of the finer of two grids, Roache (1994) / Celik et al. (2008) Eq. (7):
 *   GCI = Fs · e_a / (r^p − 1),  e_a = |(f1 − f2) / f1|.
 * Fs = 1.25 when the order is observed from three or more grids, 3.0 for a two-grid study with an assumed order.
 * @param {{ f1: number, f2: number, r21: number, p: number, Fs?: number }} o
 * @returns {{ gci: number, ea: number, eext: number, extrapolated: number }} GCI and approximate relative error
 *   (fractions, not per cent), extrapolated relative error |(f_ext − f1) / f_ext| and the extrapolated value.
 *   When f1 = 0 the errors are returned as absolute differences.
 */
export function gci({ f1, f2, r21, p, Fs = 1.25 }) {
  finite(f1, 'f1'); finite(f2, 'f2'); ratio(r21, 'r21');
  if (!(p > 0)) throw new RangeError('p: order must be > 0');
  if (!(Fs >= 1)) throw new RangeError('Fs: safety factor must be >= 1');
  const rp = r21 ** p;
  const ext = (rp * f1 - f2) / (rp - 1);
  const ea = f1 !== 0 ? Math.abs((f1 - f2) / f1) : Math.abs(f1 - f2);
  const eext = ext !== 0 ? Math.abs((ext - f1) / ext) : Math.abs(ext - f1);
  return { gci: Fs * ea / (rp - 1), ea, eext, extrapolated: ext };
}

/**
 * Asymptotic-range check (Roache 1998, "Verification and Validation in Computational Science and Engineering",
 * §5.10): for solutions in the asymptotic range GCI32 = r21^p · GCI21, so the ratio below is close to one.
 * @param {{ f1, f2, f3, r21, r32, p }} o
 * @returns {{ ratio: number, inRange: boolean, gci21: number, gci32: number }} inRange when |ratio − 1| ≤ `tol`
 *   (default 0.05)
 */
export function asymptoticCheck({ f1, f2, f3, r21, r32, p, tol = 0.05 }) {
  const g21 = gci({ f1, f2, r21, p }).gci;
  const g32 = gci({ f1: f2, f2: f3, r21: r32, p }).gci;
  const rr = g32 / (r21 ** p * g21);
  return { ratio: rr, inRange: Math.abs(rr - 1) <= tol, gci21: g21, gci32: g32 };
}

/**
 * Complete three-grid analysis (Celik et al. 2008, steps 3–5).
 * @param {{ f1, f2, f3, r21, r32 }} o solutions on the fine, medium and coarse grids and the size ratios; or give
 *   `h: [h1, h2, h3]` instead of the ratios
 * @param {{ formalOrder?: number, Fs?: number }} [opts] formal order of the scheme, used (with Fs = 3) when the
 *   observed order cannot be established, and to cap an implausibly high observed order at twice its value
 * @returns {{ type: string, R: number, p: number, pUsed: number, extrapolated: number, gciFine: number,
 *   gciMedium: number, ea21: number, eext21: number, asymptoticRatio: number, inAsymptoticRange: boolean,
 *   Fs: number, notes: string[] }} errors and GCI as fractions of the fine-grid value
 */
export function threeGrid({ f1, f2, f3, r21, r32, h }, { formalOrder, Fs } = {}) {
  if (h) { r21 = h[1] / h[0]; r32 = h[2] / h[1]; }
  ratio(r21, 'r21'); ratio(r32, 'r32');
  const cls = classify(f1, f2, f3, { r21, r32 });
  const notes = [];
  const out = { type: cls.type, R: cls.R, p: NaN, pUsed: NaN, extrapolated: f1, gciFine: 0, gciMedium: 0, ea21: 0, eext21: 0, asymptoticRatio: NaN, inAsymptoticRange: false, Fs: Fs ?? 1.25, notes };
  if (cls.type === 'converged') {
    notes.push('The two finest solutions are identical to working precision: the result is mesh independent.');
    out.inAsymptoticRange = true;
    return out;
  }
  const ord = observedOrder({ f1, f2, f3, r21, r32 });
  out.p = ord.p;
  let p = ord.p, fs = Fs ?? 1.25;
  const usable = Number.isFinite(p) && p > 0 && ord.converged && (cls.type === 'monotone' || cls.type === 'oscillatory');
  if (!usable) {
    if (!(formalOrder > 0)) {
      notes.push(cls.type.includes('divergent') ? 'The differences between successive solutions grow as the mesh is refined: the solutions are not converging and no error estimate can be given.' : 'The observed order could not be established from these three solutions.');
      out.ea21 = gci({ f1, f2, r21, p: 1 }).ea;
      out.gciFine = NaN; out.gciMedium = NaN; out.extrapolated = NaN; out.eext21 = NaN;
      return out;
    }
    p = formalOrder; fs = Fs ?? 3;
    notes.push(`The observed order is not usable; the formal order ${formalOrder} was assumed with safety factor ${fs}.`);
  } else if (formalOrder > 0 && p > 2 * formalOrder) {
    notes.push(`The observed order ${p.toFixed(2)} is more than twice the formal order ${formalOrder}; it was limited to ${2 * formalOrder} (super-convergence on three grids is usually fortuitous).`);
    p = 2 * formalOrder;
  } else if (p < 0.5) {
    notes.push(`The observed order ${p.toFixed(2)} is very low: the solutions are probably not yet in the asymptotic range.`);
  }
  if (cls.type === 'oscillatory') notes.push('The solutions oscillate about the limit; the error band is indicative only.');
  out.pUsed = p; out.Fs = fs;
  const g = gci({ f1, f2, r21, p, Fs: fs });
  out.gciFine = g.gci; out.ea21 = g.ea; out.eext21 = g.eext; out.extrapolated = g.extrapolated;
  out.gciMedium = gci({ f1: f2, f2: f3, r21: r32, p, Fs: fs }).gci;
  const a = asymptoticCheck({ f1, f2, f3, r21, r32, p });
  out.asymptoticRatio = a.ratio; out.inAsymptoticRange = usable && a.inRange;
  return out;
}

/**
 * Least-squares fit of f(h) = f0 + C·h^p to four or more grids (Eça & Hoekstra 2014, Eq. 13): for each trial
 * order the two linear coefficients follow from a 2 × 2 normal system; the order minimising the residual is
 * found by a scan refined with a golden-section search.
 * @param {ArrayLike<number>} h representative cell sizes (any order, at least three distinct values, four for a
 *   meaningful residual)
 * @param {ArrayLike<number>} f solution values
 * @param {{ pMin?: number, pMax?: number, weights?: ArrayLike<number> }} [opts] search range for the order
 *   (default 0.1 … 8); weights (default 1; Eça & Hoekstra suggest 1/h to favour the fine grids)
 * @returns {{ f0: number, C: number, p: number, residualStd: number, atBound: boolean }} extrapolated value,
 *   coefficient, order, standard deviation of the fit (units of f; n − 3 degrees of freedom) and whether the
 *   order ran into a limit of the search range (the data then do not show clean power-law convergence)
 */
export function fitOrder(h, f, { pMin = 0.1, pMax = 8, weights } = {}) {
  const n = h.length;
  if (n < 3 || f.length !== n) throw new RangeError('h/f: need at least three grids and equal lengths');
  for (let i = 0; i < n; i++) { if (!(h[i] > 0)) throw new RangeError(`h[${i}]: must be > 0`); finite(f[i], `f[${i}]`); }
  const hs = Math.max(...h);
  const solve = (p) => {
    let sw = 0, sx = 0, sxx = 0, sy = 0, sxy = 0;
    for (let i = 0; i < n; i++) { const w = weights ? weights[i] : 1, x = (h[i] / hs) ** p; sw += w; sx += w * x; sxx += w * x * x; sy += w * f[i]; sxy += w * x * f[i]; }
    const det = sw * sxx - sx * sx;
    if (!(Math.abs(det) > 1e-300)) return { f0: sy / sw, C: 0, sse: Infinity };
    const f0 = (sy * sxx - sx * sxy) / det, c = (sw * sxy - sx * sy) / det;
    let sse = 0;
    for (let i = 0; i < n; i++) { const w = weights ? weights[i] : 1; sse += w * (f[i] - f0 - c * (h[i] / hs) ** p) ** 2; }
    return { f0, C: c / hs ** p, sse };
  };
  let best = pMin, bestS = Infinity;
  const N = 80;
  for (let i = 0; i <= N; i++) { const p = pMin + (pMax - pMin) * i / N; const s = solve(p).sse; if (s < bestS) { bestS = s; best = p; } }
  const step = (pMax - pMin) / N;
  const lo = Math.max(pMin, best - step), hi = Math.min(pMax, best + step);
  const g = goldenSection((p) => solve(p).sse, lo, hi, { tol: 1e-10 });
  const p = g.x;
  const r = solve(p);
  return { f0: r.f0, C: r.C, p, residualStd: n > 3 ? Math.sqrt(r.sse / (n - 3)) : 0, atBound: p <= pMin + 1e-6 || p >= pMax - 1e-6 };
}

/**
 * Choose the coarsest level whose estimated discretisation error meets a tolerance.
 * The error of level i is estimated as Fs · |f_i − f_ext| relative to |f_ext|.
 * @param {{ levels: Array<{ h: number, value: number, cost?: number }>, extrapolated: number, tolerance: number,
 *   Fs?: number, costExponent?: number }} o levels in any order; tolerance as a fraction (0.01 = 1 %); safety
 *   factor on the error estimate (default 1.25); when a level has no measured cost its cost is taken as
 *   h^−costExponent (default 1: cost proportional to the element count along one direction — use the number of
 *   space dimensions, plus one for explicit time stepping, when h is a cell size)
 * @returns {{ index: number, level: object|null, errors: number[], relativeCost: number, saving: number, met: boolean }}
 *   index into `levels` (−1 when no level meets the tolerance; the finest is then reported with met = false),
 *   estimated relative errors per level, cost of the chosen level relative to the finest and the saving 1 − that
 */
export function selectLevel({ levels, extrapolated, tolerance, Fs = 1.25, costExponent = 1 }) {
  if (!Array.isArray(levels) || !levels.length) throw new RangeError('levels: at least one level is required');
  if (!(tolerance > 0)) throw new RangeError('tolerance: must be > 0');
  const scale = Math.abs(extrapolated) > 0 ? Math.abs(extrapolated) : 1;
  const errors = levels.map((l) => Fs * Math.abs(l.value - extrapolated) / scale);
  const cost = (l) => (l.cost > 0 ? l.cost : l.h ** -costExponent);
  let finest = 0;
  levels.forEach((l, i) => { if (l.h < levels[finest].h) finest = i; });
  let pick = -1;
  levels.forEach((l, i) => { if (errors[i] <= tolerance && (pick < 0 || l.h > levels[pick].h)) pick = i; });
  const met = pick >= 0;
  const idx = met ? pick : finest;
  const rel = cost(levels[idx]) / cost(levels[finest]);
  return { index: met ? idx : -1, level: levels[idx], errors, relativeCost: rel, saving: 1 - rel, met };
}

/**
 * Run a refinement study: evaluate the quantity of interest on every level, estimate the order, extrapolated
 * value and GCI from the three finest levels (plus a least-squares fit when there are more), and pick the
 * coarsest level that meets the tolerance.
 *
 * @param {object} o
 * @param {number[]} o.levels discretisation control per level, at least three (any order)
 * @param {(level: number, index: number) => number | { value: number, cost?: number } | Promise<…>} o.evaluate
 *   returns the quantity of interest (and optionally its cost); wall-clock time is used as the cost otherwise
 * @param {number} [o.tolerance] acceptable relative discretisation error (fraction), default 0.01
 * @param {(fraction: number, message: string) => void} [o.onProgress]
 * @param {boolean} [o.isSize] true when a level is an element size or time step; false (default) when it is a
 *   count, in which case h ∝ level^(−1 / dimension)
 * @param {number} [o.dimension] number of directions the count is spread over (default 1: the count is per
 *   direction, or the model is one-dimensional)
 * @param {number} [o.formalOrder] formal order of the method (see `threeGrid`)
 * @param {AbortSignal} [o.signal] aborts between evaluations
 * @returns {Promise<{ levels: Array<{ level, h, value, cost, error }>, type, p, pUsed, extrapolated, gciFine,
 *   asymptoticRatio, inAsymptoticRange, fit, selected: { level, index, error, relativeCost, saving } | null,
 *   converged: boolean, tolerance: number, notes: string[] }>} levels sorted fine → coarse with the estimated
 *   relative error of each; `converged` is true when the fine-grid GCI is within the tolerance
 */
export async function runStudy({ levels, evaluate, tolerance = 0.01, onProgress, isSize = false, dimension = 1, formalOrder, signal }) {
  if (!Array.isArray(levels) || levels.length < 3) throw new RangeError('levels: at least three refinement levels are required');
  if (typeof evaluate !== 'function') throw new RangeError('evaluate: a function is required');
  if (!(tolerance > 0)) throw new RangeError('tolerance: must be > 0');
  if (!(dimension >= 1)) throw new RangeError('dimension: must be >= 1');
  const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
  const rows = [];
  for (let i = 0; i < levels.length; i++) {
    if (signal && signal.aborted) throw new Error('The convergence study was cancelled.');
    const level = levels[i];
    if (!(level > 0)) throw new RangeError(`levels[${i}]: must be > 0`);
    if (onProgress) onProgress(i / levels.length, `Level ${i + 1} of ${levels.length} (${level})`);
    const t0 = now();
    const r = await evaluate(level, i);
    const value = typeof r === 'number' ? r : r && r.value;
    if (!Number.isFinite(value)) throw new Error(`The evaluation at level ${level} did not return a finite value.`);
    const measured = now() - t0;
    rows.push({ level, h: isSize ? level : level ** (-1 / dimension), value, cost: r && r.cost > 0 ? r.cost : measured > 0 ? measured : NaN });
  }
  if (onProgress) onProgress(1, 'Analysing convergence');
  rows.sort((a, b) => a.h - b.h);
  for (let i = 1; i < rows.length; i++) if (!(rows[i].h > rows[i - 1].h)) throw new RangeError('levels: every level must be different');
  const [a, b, c] = rows;
  const tg = threeGrid({ f1: a.value, f2: b.value, f3: c.value, h: [a.h, b.h, c.h] }, { formalOrder });
  const notes = tg.notes.slice();
  let fit = null;
  if (rows.length > 3) {
    try {
      fit = fitOrder(rows.map((r) => r.h), rows.map((r) => r.value));
      if (fit.atBound) notes.push('The least-squares fit over all levels ran into the limit of its order range: the coarser levels are outside the asymptotic range.');
    } catch { fit = null; }
  }
  const ext = tg.type === 'converged' ? a.value : tg.extrapolated;
  let selected = null;
  if (Number.isFinite(ext)) {
    const timed = rows.every((r) => r.cost > 0);
    const sel = selectLevel({ levels: rows.map((r) => ({ h: r.h, value: r.value, cost: timed ? r.cost : undefined })), extrapolated: ext, tolerance, Fs: tg.Fs, costExponent: isSize ? 1 : dimension });
    rows.forEach((r, i) => { r.error = sel.errors[i]; });
    if (sel.met) selected = { level: rows[sel.index].level, index: sel.index, error: sel.errors[sel.index], relativeCost: sel.relativeCost, saving: sel.saving };
    else notes.push(`No level meets the ${(100 * tolerance).toPrecision(2)} % tolerance; the finest level has an estimated error of ${(100 * sel.errors[0]).toPrecision(2)} %. Refine further.`);
  } else rows.forEach((r) => { r.error = NaN; });
  const converged = tg.type === 'converged' || (Number.isFinite(tg.gciFine) && tg.gciFine <= tolerance);
  if (converged && !tg.inAsymptoticRange && tg.type !== 'converged') notes.push('The tolerance is met, but the three finest solutions are not yet in the asymptotic range; treat the error estimate with caution.');
  return {
    levels: rows, type: tg.type, R: tg.R, p: tg.p, pUsed: tg.pUsed, extrapolated: ext, gciFine: tg.gciFine, gciMedium: tg.gciMedium, ea21: tg.ea21,
    asymptoticRatio: tg.asymptoticRatio, inAsymptoticRange: tg.inAsymptoticRange, Fs: tg.Fs, fit, selected, converged, tolerance, notes,
  };
}

/**
 * Discrete error norms of a numerical solution against an exact (manufactured) one.
 * @param {ArrayLike<number>} numerical computed values
 * @param {ArrayLike<number>} exact exact values at the same points
 * @param {ArrayLike<number>} [weights] quadrature weights (cell volumes, lengths); default equal weights
 * @returns {{ l1: number, l2: number, linf: number, relL2: number, n: number }} L1 = Σw|e| / Σw,
 *   L2 = sqrt(Σw e² / Σw), L∞ = max|e| (units of the solution); relL2 = L2 / L2-norm of the exact solution
 */
export function errorNorms(numerical, exact, weights) {
  const n = numerical.length;
  if (!n || exact.length !== n || (weights && weights.length !== n)) throw new RangeError('numerical/exact/weights: arrays must be non-empty and of equal length');
  let sw = 0, s1 = 0, s2 = 0, se = 0, linf = 0;
  for (let i = 0; i < n; i++) {
    const w = weights ? weights[i] : 1, e = Math.abs(numerical[i] - exact[i]);
    if (!Number.isFinite(e)) throw new RangeError(`numerical[${i}]: not finite`);
    sw += w; s1 += w * e; s2 += w * e * e; se += w * exact[i] * exact[i];
    if (e > linf) linf = e;
  }
  const l2 = Math.sqrt(s2 / sw);
  return { l1: s1 / sw, l2, linf, relL2: se > 0 ? Math.sqrt(s2 / se) : NaN, n };
}

/**
 * Observed order from known errors (manufactured or exact solutions): slope of log(error) against log(h).
 * @param {ArrayLike<number>} h cell sizes (at least two)
 * @param {ArrayLike<number>} errors error norm on each grid (> 0)
 * @returns {{ p: number, pairwise: number[], coefficient: number }} least-squares order over all grids, the order
 *   between each successive pair (sorted fine → coarse) and the constant C in error ≈ C·h^p
 */
export function orderFromErrors(h, errors) {
  const n = h.length;
  if (n < 2 || errors.length !== n) throw new RangeError('h/errors: need at least two grids and equal lengths');
  const idx = Array.from({ length: n }, (_, i) => i).sort((i, j) => h[i] - h[j]);
  const x = idx.map((i) => { if (!(h[i] > 0)) throw new RangeError('h: must be > 0'); return Math.log(h[i]); });
  const y = idx.map((i) => { if (!(errors[i] > 0)) throw new RangeError('errors: must be > 0'); return Math.log(errors[i]); });
  const mx = x.reduce((s, v) => s + v, 0) / n, my = y.reduce((s, v) => s + v, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sxx += (x[i] - mx) ** 2; sxy += (x[i] - mx) * (y[i] - my); }
  if (!(sxx > 0)) throw new RangeError('h: the grids must differ in size');
  const p = sxy / sxx;
  const pairwise = [];
  for (let i = 0; i + 1 < n; i++) pairwise.push((y[i + 1] - y[i]) / (x[i + 1] - x[i]));
  return { p, pairwise, coefficient: Math.exp(my - p * mx) };
}
